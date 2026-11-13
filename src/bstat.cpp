// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of crimlab's count and table statistics: count
// regression (cntreg), rates by group (grprr) and standardised rates
// (dastd, sirex), and contingency tables (ctabx). Each checks its arguments, builds the kernel's
// arrays, calls the kernels and builds the R result. Distribution
// functions come from R's maths library. Nothing calls back into the R
// evaluator.
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

#include "crimlab.h"
#include "rchk.h"
#include "spin.h"

// Older R headers declare these under their short names only when
// R_NO_REMAP_RMATH is set; the library has them as Rf_*, and newer headers
// declare the same names identically (rchk.h does this for Rf_qnorm5).
extern "C" double Rf_pnorm5(double, double, double, int, int);
extern "C" double Rf_qgamma(double, double, double, int, int);
extern "C" double Rf_qbeta(double, double, double, int, int);
extern "C" double Rf_qchisq(double, double, int, int);
extern "C" double Rf_qpois(double, double, int, int);
extern "C" double Rf_pchisq(double, double, int, int);
extern "C" double Rf_pnchisq(double, double, double, int, int);

using clab::fail;

namespace {

const double qnan = std::numeric_limits<double>::quiet_NaN();

int len(SEXP x, const char *nm) { return clab::ilen(static_cast<double>(XLENGTH(x)), nm); }

// Columns of a result, all protected until the frame is built.
struct cols {
  std::vector<std::string> nm;
  SEXP v = R_NilValue;
  int k = 0;
  explicit cols(int ncol) {
    v = PROTECT(Rf_allocVector(VECSXP, ncol));
    nm.reserve(static_cast<std::size_t>(ncol));
  }
  double *dbl(const char *name, int n) {
    SEXP c = Rf_allocVector(REALSXP, n);
    SET_VECTOR_ELT(v, k++, c);
    nm.emplace_back(name);
    return REAL(c);
  }
  int *lgc(const char *name, int n) {
    SEXP c = Rf_allocVector(LGLSXP, n);
    SET_VECTOR_ELT(v, k++, c);
    nm.emplace_back(name);
    return LOGICAL(c);
  }
  SEXP str(const char *name, int n) {
    SEXP c = Rf_allocVector(STRSXP, n);
    SET_VECTOR_ELT(v, k++, c);
    nm.emplace_back(name);
    return c;
  }
  void put(const char *name, SEXP c) {
    SET_VECTOR_ELT(v, k++, c);
    nm.emplace_back(name);
  }
  // the data frame; unprotects the list
  SEXP frame() {
    std::vector<const char *> p;
    for (const auto &s : nm) p.push_back(s.c_str());
    // the columns added, however many were reserved
    SEXP w = PROTECT(Rf_allocVector(VECSXP, k));
    for (int j = 0; j < k; ++j) SET_VECTOR_ELT(w, j, VECTOR_ELT(v, j));
    SEXP f = clab::frame(w, p.data());
    UNPROTECT(2);
    return f;
  }
};

// A named list of the given elements (protected by the caller).
SEXP nlist(std::initializer_list<std::pair<const char *, SEXP>> e) {
  SEXP v = PROTECT(Rf_allocVector(VECSXP, static_cast<R_xlen_t>(e.size())));
  SEXP n = PROTECT(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(e.size())));
  R_xlen_t i = 0;
  for (const auto &x : e) {
    SET_VECTOR_ELT(v, i, x.second);
    SET_STRING_ELT(n, i, Rf_mkChar(x.first));
    ++i;
  }
  Rf_setAttrib(v, R_NamesSymbol, n);
  UNPROTECT(2);
  return v;
}

double zq(SEXP level) {
  const double lv = clab::one(level, "level", 0.0, 1.0);
  return Rf_qnorm5(0.5 + 0.5 * lv, 0.0, 1.0, 1, 0);
}

// ------------------------------------------------------------ cntreg

// One column of the design: its values and its term label.
struct term {
  std::vector<double> x;
  std::string nm;
};

// The adjustment columns: a data frame or list of columns, a numeric
// matrix (its columns), or one vector. Labels (factor, text, logical)
// become indicators of every level but the first; numbers enter as they
// are.
void adjcols(SEXP adj, int n, std::vector<term> &out) {
  if (Rf_isNull(adj)) return;
  std::vector<std::pair<std::string, SEXP>> cs;
  if (TYPEOF(adj) == VECSXP) {
    SEXP nms = Rf_getAttrib(adj, R_NamesSymbol);
    for (R_xlen_t j = 0; j < XLENGTH(adj); ++j) {
      std::string nm = !Rf_isNull(nms) && STRING_ELT(nms, j) != NA_STRING &&
                               CHAR(STRING_ELT(nms, j))[0] != '\0'
                           ? Rf_translateCharUTF8(STRING_ELT(nms, j))
                           : "adj" + std::to_string(j + 1);
      cs.emplace_back(nm, VECTOR_ELT(adj, j));
    }
  } else if (clab::isnum(adj) && Rf_isMatrix(adj)) {
    const int nr = Rf_nrows(adj), nc = Rf_ncols(adj);
    if (nr != n) fail("`adj` must have one row per count, %d.", n);
    SEXP dn = Rf_getAttrib(adj, R_DimNamesSymbol);
    SEXP cn = Rf_isNull(dn) ? R_NilValue : VECTOR_ELT(dn, 1);
    SEXP v = PROTECT(clab::dbl(adj));
    for (int j = 0; j < nc; ++j) {
      term t;
      t.nm = !Rf_isNull(cn) && STRING_ELT(cn, j) != NA_STRING
                 ? Rf_translateCharUTF8(STRING_ELT(cn, j))
                 : "adj" + std::to_string(j + 1);
      t.x.assign(REAL(v) + static_cast<std::size_t>(j) * static_cast<std::size_t>(n),
                 REAL(v) + static_cast<std::size_t>(j + 1) * static_cast<std::size_t>(n));
      for (double e : t.x) {
        if (!std::isfinite(e)) fail("`adj` must hold finite numbers with no missing values.");
      }
      out.push_back(t);
    }
    UNPROTECT(1);
    return;
  } else {
    cs.emplace_back("adj", adj);
  }
  for (auto &c : cs) {
    SEXP x = c.second;
    const std::string an = "adj " + c.first;
    clab::unlab(x, an.c_str());
    if (!Rf_isVectorAtomic(x) || XLENGTH(x) != n) {
      fail("`adj` column %s must be a vector with one value per count, %d.", c.first.c_str(), n);
    }
    if (clab::isnum(x) && !Rf_isFactor(x)) {
      SEXP v = PROTECT(clab::num(x, an.c_str(), n));
      term t;
      t.nm = c.first;
      t.x.assign(REAL(v), REAL(v) + n);
      out.push_back(t);
      UNPROTECT(1);
      continue;
    }
    csp::lab l = csp::codes(x, an.c_str(), n);
    for (int k = 1; k < l.g; ++k) {
      term t;
      t.nm = c.first + "=" + Rf_translateCharUTF8(STRING_ELT(l.lev, k));
      t.x.resize(static_cast<std::size_t>(n));
      for (int i = 0; i < n; ++i) t.x[static_cast<std::size_t>(i)] = INTEGER(l.code)[i] == k ? 1.0 : 0.0;
      out.push_back(t);
    }
    UNPROTECT(2);
  }
}

// A kernel's verdict: CLAB_OK, or stop with msg.
void kok(int rc, const char *msg) {
  if (rc != CLAB_OK) fail("%s", msg);
}

// the Poisson fit's failure: its arguments were checked, so the IRLS
// left the finite range
const char *const cntmsg =
    "the fit failed: the deviance or a coefficient left the finite range "
    "(a group or cell with no counts can do this).";

// ------------------------------------------------------------ rates

double na(double v) { return std::isnan(v) ? NA_REAL : v; }

// Clopper-Pearson limits of x of n at tail a / 2 each side
double cplo(double x, double n, double a) { return x > 0 ? Rf_qbeta(a / 2, x, n - x + 1, 1, 0) : 0.0; }
double cphi(double x, double n, double a) { return x < n ? Rf_qbeta(1 - a / 2, x + 1, n - x, 1, 0) : 1.0; }
// the exact Poisson limits of a count x
double gmlo(double x, double a) { return x > 0 ? Rf_qgamma(a / 2, x, 1.0, 1, 0) : 0.0; }
double gmhi(double x, double a) { return Rf_qgamma(1 - a / 2, x + 1, 1.0, 1, 0); }

// Per-group sums of x and n over rows labelled by grp (NULL: each row its
// own group, labelled by the names of x or 1, 2, ...). lev gets the
// labels; protects two (lev and one scratch) or one.
struct grpd {
  std::vector<double> x, n;
  std::vector<int> code;  // the group of each row
  SEXP lev;
  int g;
};

grpd groups(SEXP x, SEXP grp, int m, const double *xv, const double *nv, int *np) {
  grpd r;
  if (Rf_isNull(grp)) {
    SEXP nm = Rf_getAttrib(x, R_NamesSymbol);
    r.lev = PROTECT(Rf_allocVector(STRSXP, m));
    ++*np;
    for (int i = 0; i < m; ++i) {
      if (!Rf_isNull(nm) && STRING_ELT(nm, i) != NA_STRING) {
        SET_STRING_ELT(r.lev, i, STRING_ELT(nm, i));
      } else {
        SET_STRING_ELT(r.lev, i, Rf_mkChar(std::to_string(i + 1).c_str()));
      }
    }
    r.g = m;
    r.code.resize(static_cast<std::size_t>(m));
    for (int i = 0; i < m; ++i) r.code[static_cast<std::size_t>(i)] = i;
  } else {
    csp::lab l = csp::codes(grp, "grp", m);
    *np += 2;
    r.lev = l.lev;
    r.g = l.g;
    r.code.assign(INTEGER(l.code), INTEGER(l.code) + m);
  }
  r.x.assign(static_cast<std::size_t>(r.g), 0.0);
  r.n.assign(static_cast<std::size_t>(r.g), 0.0);
  for (int i = 0; i < m; ++i) {
    const std::size_t c = static_cast<std::size_t>(r.code[static_cast<std::size_t>(i)]);
    r.x[c] += xv[i];
    r.n[c] += nv ? nv[i] : 1.0;
  }
  return r;
}

// Weights or rates given per stratum: named ones matched to the strata by
// label, unnamed ones in the strata's order. Returns ns values.
std::vector<double> perstr(SEXP v, const char *nm, SEXP lev, int ns) {
  SEXP d = PROTECT(clab::num(v, nm, -1, 0.0));
  SEXP vn = Rf_getAttrib(v, R_NamesSymbol);
  std::vector<double> out(static_cast<std::size_t>(ns));
  if (Rf_isNull(vn)) {
    if (XLENGTH(d) != ns) fail("`%s` must have one value per stratum, %d.", nm, ns);
    for (int s = 0; s < ns; ++s) out[static_cast<std::size_t>(s)] = REAL(d)[s];
  } else {
    for (int s = 0; s < ns; ++s) {
      int hit = -1;
      for (R_xlen_t j = 0; j < XLENGTH(d); ++j) {
        if (clab::seql(STRING_ELT(vn, j), STRING_ELT(lev, s))) {
          if (hit >= 0) fail("`%s` names stratum %s twice.", nm, Rf_translateChar(STRING_ELT(lev, s)));
          hit = static_cast<int>(j);
        }
      }
      if (hit < 0) fail("`%s` has no value for stratum %s.", nm, Rf_translateChar(STRING_ELT(lev, s)));
      out[static_cast<std::size_t>(s)] = REAL(d)[hit];
    }
  }
  UNPROTECT(1);
  return out;
}

const char *const rtype[] = {"rate", "prop"};

// ------------------------------------------------------------ tables

// The noncentrality at which a chi-square statistic st on df degrees of
// freedom sits at lower-tail probability q: pnchisq(st, df, ncp) falls as
// ncp grows; 0 when even ncp = 0 leaves it below q.
double ncp(double st, double df, double q) {
  if (Rf_pchisq(st, df, 1, 0) <= q) return 0.0;
  double a = 0, b = std::fmax(1.0, st);
  while (Rf_pnchisq(st, df, b, 1, 0) > q) b *= 2;
  for (int it = 0; it < 200 && b - a > 1e-12 * b; ++it) {
    const double m = (a + b) / 2;
    if (Rf_pnchisq(st, df, m, 1, 0) > q) a = m;
    else b = m;
  }
  return (a + b) / 2;
}

const char *const aggs[] = {"n", "sum", "max", "min", "mean", "dist"};

SEXP labmat(const double *v, int r, int c, SEXP rl, SEXP cl) {
  SEXP m = PROTECT(Rf_allocMatrix(REALSXP, r, c));
  std::copy(v, v + static_cast<std::ptrdiff_t>(r) * c, REAL(m));
  SEXP dn = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dn, 0, rl);
  SET_VECTOR_ELT(dn, 1, cl);
  Rf_setAttrib(m, R_DimNamesSymbol, dn);
  UNPROTECT(2);
  return m;
}

// labels 1, 2, ..., n
SEXP seqlab(int n) {
  SEXP s = PROTECT(Rf_allocVector(STRSXP, n));
  for (int i = 0; i < n; ++i) SET_STRING_ELT(s, i, Rf_mkChar(std::to_string(i + 1).c_str()));
  UNPROTECT(1);
  return s;
}

}  // namespace

extern "C" {

//' @srrstats {G2.6} cntreg() reads its count, group, exposure and cluster
//'   arguments as plain vectors of one length, named ones matched by name;
//'   adjustment columns come as a data frame, a list, a matrix or one
//'   vector.
//' @srrstats {G2.4c} Groups, clusters and labelled adjustment columns are
//'   read as labels (as.character(), a factor by its levels).
//' @noRd
SEXP C_clab_r_cntreg(SEXP y, SEXP grp, SEXP adj, SEXP expo, SEXP clus, SEXP ref,
                     SEXP gee, SEXP cadj, SEXP level) {
  clab::who = "cntreg";
  SEXP uv[4] = {y, grp, expo, clus};
  static const char *const un[] = {"y", "grp", "expo", "clus"};
  int np = csp::byname(uv, un, 4);
  SEXP yv = PROTECT(clab::num(uv[0], "y", -1, 0.0));
  ++np;
  const int n = len(yv, "y");
  const double *yy = REAL(yv);
  double tot = 0.0;
  for (int i = 0; i < n; ++i) tot += yy[i];
  if (!(tot > 0)) fail("`y` has no events: every count is 0.");
  if (Rf_isNull(uv[1])) fail("`grp` must give the group of each count.");
  csp::lab g = csp::codes(uv[1], "grp", n);
  np += 2;
  if (g.g < 2) fail("`grp` has one group: a rate ratio needs two.");
  const int r0 = Rf_isNull(ref) ? 0 : clab::ref(ref, g.lev, g.g) - 1;
  const double z = zq(level);
  const bool dogee = csp::flag(gee, "gee"), adjc = csp::flag(cadj, "cadj");
  std::vector<double> off;
  if (!Rf_isNull(uv[2])) {
    SEXP e = PROTECT(clab::num(uv[2], "expo", n, 0.0, true));
    ++np;
    off.resize(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) off[static_cast<std::size_t>(i)] = std::log(REAL(e)[i]);
  }
  int ncl = 0;
  const int *cl = nullptr;
  if (!Rf_isNull(uv[3])) {
    csp::lab c = csp::codes(uv[3], "clus", n);
    np += 2;
    if (c.g < 2) fail("`clus` has one cluster: a cluster-robust variance needs two or more.");
    ncl = c.g;
    cl = INTEGER(c.code);
  }
  if (dogee && !cl) fail("`gee = TRUE` needs `clus`.");
  // the design: intercept, the groups but the reference, then adj
  std::vector<term> ts;
  {
    term t;
    t.nm = "intercept";
    t.x.assign(static_cast<std::size_t>(n), 1.0);
    ts.push_back(t);
  }
  std::vector<int> gt;  // the group of each group term
  for (int k = 0; k < g.g; ++k) {
    if (k == r0) continue;
    term t;
    t.nm = Rf_translateCharUTF8(STRING_ELT(g.lev, k));
    t.x.resize(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) t.x[static_cast<std::size_t>(i)] = INTEGER(g.code)[i] == k ? 1.0 : 0.0;
    ts.push_back(t);
    gt.push_back(k);
  }
  adjcols(adj, n, ts);
  const int p = static_cast<int>(ts.size());
  const std::size_t sn = static_cast<std::size_t>(n), sp = static_cast<std::size_t>(p);
  std::vector<double> x(sn * sp);
  for (int j = 0; j < p; ++j) {
    std::copy(ts[static_cast<std::size_t>(j)].x.begin(), ts[static_cast<std::size_t>(j)].x.end(),
              x.begin() + static_cast<std::ptrdiff_t>(static_cast<std::size_t>(j) * sn));
  }
  const std::size_t sc = static_cast<std::size_t>(ncl);
  std::vector<double> work(sn * sp + 6 * sn + sp + sp * sp + 2 * sc * sp + 2 * sn);
  std::vector<int> iwork(sp + sc + 1 + sn);
  // fits: 0 Poisson, 1 negative binomial, 2 GEE
  std::vector<double> beta(3 * sp, qnan), vb(2 * sp * sp, qnan), vc(3 * sp * sp, qnan),
      st(3 * 8, qnan);
  std::vector<int> alias(2 * sp, 1);
  const double *of = off.empty() ? nullptr : off.data();
  bool ok[3] = {false, false, false};
  for (int f = 0; f < 2; ++f) {
    const std::size_t fo = static_cast<std::size_t>(f);
    const int rc = clab_cntreg(yy, x.data(), n, p, of, cl, ncl, f, 25, 1e-8, &beta[fo * sp],
                               &vb[fo * sp * sp], &vc[fo * sp * sp], &st[fo * 8],
                               &alias[fo * sp], work.data(), iwork.data());
    if (f == 0) kok(rc, cntmsg);
    if (rc == CLAB_OK) {
      ok[f] = true;
    } else {
      Rf_warningcall(R_NilValue,
                     "cntreg: the negative binomial fit failed (theta has no finite "
                     "estimate); its rows are NA.");
    }
  }
  if (dogee) {
    std::copy(beta.begin(), beta.begin() + static_cast<std::ptrdiff_t>(sp), beta.begin() + static_cast<std::ptrdiff_t>(2 * sp));
    const int rc = clab_geeex(yy, x.data(), n, p, of, cl, ncl, alias.data(), 100, 1e-10,
                              &beta[2 * sp], &vc[2 * sp * sp], &st[16], work.data(),
                              iwork.data());
    if (rc == CLAB_OK) {
      ok[2] = true;
    } else {
      Rf_warningcall(R_NilValue,
                     "cntreg: the GEE fit failed (the exchangeable working correlation "
                     "is not positive definite, or the scoring diverged); its rows are NA.");
    }
  }
  const double fc = adjc && ncl > 1 ? static_cast<double>(ncl) / (ncl - 1) : 1.0;
  // the rows of coef: (fit, variance) pairs
  struct row {
    const char *lab;
    int f;     // the fit
    int v;     // 0 model-based (vb), 1 sandwich (vc)
  };
  std::vector<row> rows = {{"pois", 0, 0}, {"nb", 1, 0}};
  if (cl) {
    rows.push_back({"poiscl", 0, 1});
    rows.push_back({"nbcl", 1, 1});
  }
  if (dogee) rows.push_back({"gee", 2, 1});
  const int nr = static_cast<int>(rows.size());
  const int ng = static_cast<int>(gt.size());
  // separation: a group term whose group, or the reference, has no counts
  std::vector<double> gsum(static_cast<std::size_t>(g.g), 0.0);
  for (int i = 0; i < n; ++i) gsum[static_cast<std::size_t>(INTEGER(g.code)[i])] += yy[i];
  auto est = [&](const row &w, int j, double *b, double *s) {
    const std::size_t fo = static_cast<std::size_t>(w.f), jj = static_cast<std::size_t>(j);
    *b = ok[w.f] ? beta[fo * sp + jj] : qnan;
    const double v = !ok[w.f] ? qnan
                     : w.v == 0 ? vb[fo * sp * sp + jj * sp + jj]
                                : vc[fo * sp * sp + jj * sp + jj] * (w.f == 2 ? 1.0 : fc);
    *s = std::sqrt(v);
  };
  cols cf(6);
  SEXP cfit = cf.str("fit", nr * p), cterm = cf.str("term", nr * p);
  double *ce = cf.dbl("est", nr * p), *cs = cf.dbl("se", nr * p), *cz = cf.dbl("z", nr * p),
         *cp = cf.dbl("p", nr * p);
  cols ir(7);
  SEXP ifit = ir.str("fit", nr * ng), iterm = ir.str("term", nr * ng);
  double *ii = ir.dbl("irr", nr * ng), *il = ir.dbl("lo", nr * ng), *ih = ir.dbl("hi", nr * ng),
         *ip = ir.dbl("p", nr * ng);
  int *is = ir.lgc("sep", nr * ng);
  for (int a = 0; a < nr; ++a) {
    for (int j = 0; j < p; ++j) {
      const int q = a * p + j;
      double b, s;
      est(rows[static_cast<std::size_t>(a)], j, &b, &s);
      SET_STRING_ELT(cfit, q, Rf_mkChar(rows[static_cast<std::size_t>(a)].lab));
      SET_STRING_ELT(cterm, q, Rf_mkCharCE(ts[static_cast<std::size_t>(j)].nm.c_str(), CE_UTF8));
      ce[q] = std::isnan(b) ? NA_REAL : b;
      cs[q] = std::isnan(s) ? NA_REAL : s;
      cz[q] = std::isnan(b / s) ? NA_REAL : b / s;
      cp[q] = std::isnan(b / s) ? NA_REAL : 2 * Rf_pnorm5(-std::fabs(b / s), 0.0, 1.0, 1, 0);
    }
    for (int k = 0; k < ng; ++k) {
      const int q = a * ng + k, j = k + 1;
      double b, s;
      est(rows[static_cast<std::size_t>(a)], j, &b, &s);
      SET_STRING_ELT(ifit, q, Rf_mkChar(rows[static_cast<std::size_t>(a)].lab));
      SET_STRING_ELT(iterm, q, Rf_mkCharCE(ts[static_cast<std::size_t>(j)].nm.c_str(), CE_UTF8));
      const bool fin = !std::isnan(b) && !std::isnan(s);
      ii[q] = std::isnan(b) ? NA_REAL : std::exp(b);
      il[q] = fin ? std::exp(b - z * s) : NA_REAL;
      ih[q] = fin ? std::exp(b + z * s) : NA_REAL;
      ip[q] = fin ? 2 * Rf_pnorm5(-std::fabs(b / s), 0.0, 1.0, 1, 0) : NA_REAL;
      is[q] = gsum[static_cast<std::size_t>(gt[static_cast<std::size_t>(k)])] == 0 ||
                      gsum[static_cast<std::size_t>(r0)] == 0
                  ? 1
                  : 0;
    }
  }
  SEXP coef = PROTECT(cf.frame());
  SEXP irr = PROTECT(ir.frame());
  np += 2;
  // the fits
  const int nf = dogee ? 3 : 2;
  static const char *const fl[] = {"pois", "nb", "gee"};
  cols fs(9);
  SEXP ffit = fs.str("fit", nf);
  double *fll = fs.dbl("loglik", nf), *fa = fs.dbl("aic", nf), *ft = fs.dbl("theta", nf),
         *fst = fs.dbl("setheta", nf), *fal = fs.dbl("alpha", nf), *fit = fs.dbl("iter", nf);
  int *fcv = fs.lgc("conv", nf);
  double *frk = fs.dbl("rank", nf);
  for (int f = 0; f < nf; ++f) {
    const double *s = &st[static_cast<std::size_t>(f) * 8];
    SET_STRING_ELT(ffit, f, Rf_mkChar(fl[f]));
    auto na = [&](double v) { return ok[f] && !std::isnan(v) ? v : NA_REAL; };
    if (f < 2) {
      fll[f] = na(s[0]);
      fa[f] = na(s[1]);
      ft[f] = na(s[2]);
      fst[f] = na(s[3]);
      fal[f] = NA_REAL;
      fit[f] = na(s[4]);
      fcv[f] = ok[f] ? (s[5] != 0 ? 1 : 0) : NA_LOGICAL;
      frk[f] = na(s[6]);
    } else {
      fll[f] = fa[f] = ft[f] = fst[f] = NA_REAL;
      fal[f] = na(s[0]);
      fit[f] = na(s[2]);
      fcv[f] = ok[f] ? (s[3] != 0 ? 1 : 0) : NA_LOGICAL;
      frk[f] = na(st[6]);
    }
  }
  SEXP fits = PROTECT(fs.frame());
  SEXP nn = PROTECT(Rf_ScalarInteger(n)), nc = PROTECT(Rf_ScalarInteger(ncl)),
       rf = PROTECT(Rf_ScalarString(STRING_ELT(g.lev, r0))), lv = PROTECT(Rf_ScalarReal(REAL_RO(clab::dbl(level))[0]));
  np += 5;
  SEXP out = PROTECT(nlist({{"irr", irr}, {"coef", coef}, {"fit", fits}, {"ref", rf},
                            {"n", nn}, {"ncl", nc}, {"level", lv}}));
  ++np;
  Rf_classgets(out, PROTECT(Rf_mkString("cntreg")));
  UNPROTECT(np + 1);
  return out;
}

SEXP C_clab_r_grprr(SEXP x, SEXP n, SEXP grp, SEXP ref, SEXP type, SEXP per,
                    SEXP level) {
  clab::who = "grprr";
  SEXP uv[3] = {x, n, grp};
  static const char *const un[] = {"x", "n", "grp"};
  int np = csp::byname(uv, un, 3);
  SEXP xv = PROTECT(clab::num(uv[0], "x", -1, 0.0));
  ++np;
  const int m = len(xv, "x");
  const int prop = csp::pick(type, "type", rtype, 2);
  const double* nv = nullptr;
  if (!Rf_isNull(uv[1])) {
    SEXP nn = PROTECT(clab::num(uv[1], "n", m, 0.0, true));
    ++np;
    nv = REAL(nn);
  }
  const double pr = clab::one(per, "per", 0.0, R_PosInf);
  const double z = zq(level), a = 1 - REAL_RO(PROTECT(clab::dbl(level)))[0];
  ++np;
  grpd g = groups(uv[0], uv[2], m, REAL(xv), nv, &np);
  for (int k = 0; k < g.g; ++k) {
    if (prop && g.x[static_cast<std::size_t>(k)] > g.n[static_cast<std::size_t>(k)]) {
      fail("group %s has more events than trials: `x` must not exceed `n` for proportions.",
           Rf_translateChar(STRING_ELT(g.lev, k)));
    }
  }
  if (!prop && !nv) fail("`n` must give the exposure of each count for `type = \"rate\"`.");
  const int r0 = Rf_isNull(ref) ? 0 : clab::ref(ref, g.lev, g.g) - 1;
  const int k = g.g;
  const std::size_t sk = static_cast<std::size_t>(k);
  std::vector<double> o(10 * sk);
  double *rt = &o[0], *lo = &o[sk], *hi = &o[2 * sk], *rp = &o[3 * sk], *ra = &o[4 * sk],
         *rl = &o[5 * sk], *rh = &o[6 * sk], *df = &o[7 * sk], *dl = &o[8 * sk], *dh = &o[9 * sk];
  kok(clab_grprr(g.x.data(), g.n.data(), k, r0, prop, z, rt, lo, hi, rp, ra, rl, rh, df, dl, dh),
      "the kernel refused the input.");
  cols rs(9);
  rs.put("group", g.lev);
  double *cx = rs.dbl("x", k), *cn = rs.dbl("n", k), *cr = rs.dbl("rate", k),
         *cl = rs.dbl("lo", k), *ch = rs.dbl("hi", k), *el = rs.dbl("elo", k),
         *eh = rs.dbl("ehi", k), *cp = rs.dbl("rep", k);
  for (int i = 0; i < k; ++i) {
    const double xi = g.x[static_cast<std::size_t>(i)], ni = g.n[static_cast<std::size_t>(i)];
    cx[i] = xi;
    cn[i] = ni;
    cr[i] = pr * rt[i];
    cl[i] = na(pr * lo[i]);
    ch[i] = na(pr * hi[i]);
    el[i] = pr * (prop ? cplo(xi, ni, a) : gmlo(xi, a) / ni);
    eh[i] = pr * (prop ? cphi(xi, ni, a) : gmhi(xi, a) / ni);
    cp[i] = na(rp[i]);
  }
  SEXP rates = PROTECT(rs.frame());
  ++np;
  // the comparisons, every group but the reference
  SEXP cg = PROTECT(Rf_allocVector(STRSXP, k - 1));
  ++np;
  cols cm(9);
  cm.put("group", cg);
  double *mr = cm.dbl("ratio", k - 1), *ml = cm.dbl("lo", k - 1), *mh = cm.dbl("hi", k - 1),
         *me = cm.dbl("elo", k - 1), *mf = cm.dbl("ehi", k - 1), *md = cm.dbl("diff", k - 1),
         *mdl = cm.dbl("dlo", k - 1), *mdh = cm.dbl("dhi", k - 1);
  const double x0 = g.x[static_cast<std::size_t>(r0)], n0 = g.n[static_cast<std::size_t>(r0)];
  for (int i = 0, q = 0; i < k; ++i) {
    if (i == r0) continue;
    const double x1 = g.x[static_cast<std::size_t>(i)], n1 = g.n[static_cast<std::size_t>(i)];
    SET_STRING_ELT(cg, q, STRING_ELT(g.lev, i));
    mr[q] = na(ra[i]);
    ml[q] = na(rl[i]);
    mh[q] = na(rh[i]);
    if (!prop && x1 + x0 > 0) {
      // conditional on the total: x1 is binomial with p = n1 R / (n1 R + n0)
      const double pl = cplo(x1, x1 + x0, a), pu = cphi(x1, x1 + x0, a);
      me[q] = pl * n0 / ((1 - pl) * n1);
      mf[q] = pu < 1 ? pu * n0 / ((1 - pu) * n1) : R_PosInf;
    } else {
      me[q] = mf[q] = NA_REAL;
    }
    md[q] = pr * df[i];
    mdl[q] = na(pr * dl[i]);
    mdh[q] = na(pr * dh[i]);
    ++q;
  }
  SEXP cmp = PROTECT(cm.frame());
  SEXP rf = PROTECT(Rf_ScalarString(STRING_ELT(g.lev, r0)));
  SEXP ty = PROTECT(Rf_mkString(rtype[prop]));
  SEXP pe = PROTECT(Rf_ScalarReal(pr)), lv = PROTECT(Rf_ScalarReal(1 - a));
  np += 5;
  SEXP out = PROTECT(nlist({{"rates", rates}, {"cmp", cmp}, {"ref", rf}, {"type", ty},
                            {"per", pe}, {"level", lv}}));
  Rf_classgets(out, PROTECT(Rf_mkString("grprr")));
  UNPROTECT(np + 2);
  return out;
}

//' @srrstats {G2.3a} grprr() takes `type` as one of "rate" or "prop" and
//'   refuses anything else, naming the choices.
//' @noRd
SEXP C_clab_r_dastd(SEXP x, SEXP n, SEXP strat, SEXP grp, SEXP std, SEXP ref,
                    SEXP per, SEXP level) {
  clab::who = "dastd";
  SEXP uv[4] = {x, n, strat, grp};
  static const char *const un[] = {"x", "n", "strat", "grp"};
  int np = csp::byname(uv, un, 4);
  SEXP xv = PROTECT(clab::num(uv[0], "x", -1, 0.0));
  const int m = len(xv, "x");
  SEXP nv = PROTECT(clab::num(uv[1], "n", m, 0.0, true));
  np += 2;
  csp::lab s = csp::codes(uv[2], "strat", m);
  np += 2;
  csp::lab g{R_NilValue, R_NilValue, 1};
  std::vector<int> gc(static_cast<std::size_t>(m), 0);
  if (!Rf_isNull(uv[3])) {
    g = csp::codes(uv[3], "grp", m);
    np += 2;
    gc.assign(INTEGER(g.code), INTEGER(g.code) + m);
  } else {
    g.lev = PROTECT(Rf_mkString("all"));
    ++np;
  }
  std::vector<double> w = perstr(std, "std", s.lev, s.g);
  const int r0 = Rf_isNull(ref) ? 0 : clab::ref(ref, g.lev, g.g) - 1;
  const double pr = clab::one(per, "per", 0.0, R_PosInf);
  const double z = zq(level), a = 1 - REAL_RO(PROTECT(clab::dbl(level)))[0];
  ++np;
  const int k = g.g;
  const std::size_t sk = static_cast<std::size_t>(k);
  std::vector<double> o(10 * sk);
  std::vector<int> seen(sk * static_cast<std::size_t>(s.g));
  double *rt = &o[0], *vr = &o[sk], *wm = &o[2 * sk], *ev = &o[3 * sk], *ra = &o[4 * sk],
         *rl = &o[5 * sk], *rh = &o[6 * sk], *df = &o[7 * sk], *dl = &o[8 * sk], *dh = &o[9 * sk];
  const int rc = clab_dastd(REAL(xv), REAL(nv), gc.data(), INTEGER(s.code), m, k, s.g, w.data(),
                            r0, z, rt, vr, wm, ev, ra, rl, rh, df, dl, dh, seen.data());
  if (rc == CLAB_EARG) fail("each group must have one row per stratum.");
  if (rc != CLAB_OK) {
    fail("every group needs a row for each stratum of positive weight, and `std` a "
         "positive sum.");
  }
  cols rs(6);
  rs.put("group", g.lev);
  double *cx = rs.dbl("x", k), *cr = rs.dbl("rate", k), *cl = rs.dbl("lo", k),
         *ch = rs.dbl("hi", k), *cv = rs.dbl("var", k);
  for (int i = 0; i < k; ++i) {
    const double y = rt[i], v = vr[i], wx = wm[i];
    cx[i] = ev[i];
    cr[i] = pr * y;
    // Fay and Feuer's gamma interval
    cl[i] = y > 0 ? pr * v / (2 * y) * Rf_qchisq(a / 2, 2 * y * y / v, 1, 0) : 0.0;
    ch[i] = pr * (v + wx * wx) / (2 * (y + wx)) *
            Rf_qchisq(1 - a / 2, 2 * (y + wx) * (y + wx) / (v + wx * wx), 1, 0);
    cv[i] = pr * pr * v;
  }
  SEXP rates = PROTECT(rs.frame());
  SEXP cg = PROTECT(Rf_allocVector(STRSXP, k - 1));
  np += 2;
  cols cm(7);
  cm.put("group", cg);
  double *mr = cm.dbl("ratio", k - 1), *ml = cm.dbl("lo", k - 1), *mh = cm.dbl("hi", k - 1),
         *md = cm.dbl("diff", k - 1), *mdl = cm.dbl("dlo", k - 1), *mdh = cm.dbl("dhi", k - 1);
  for (int i = 0, q = 0; i < k; ++i) {
    if (i == r0) continue;
    SET_STRING_ELT(cg, q, STRING_ELT(g.lev, i));
    mr[q] = na(ra[i]);
    ml[q] = na(rl[i]);
    mh[q] = na(rh[i]);
    md[q] = pr * df[i];
    mdl[q] = pr * dl[i];
    mdh[q] = pr * dh[i];
    ++q;
  }
  SEXP cmp = PROTECT(cm.frame());
  SEXP rf = PROTECT(Rf_ScalarString(STRING_ELT(g.lev, r0)));
  SEXP pe = PROTECT(Rf_ScalarReal(pr)), lv = PROTECT(Rf_ScalarReal(1 - a));
  np += 4;
  SEXP out = PROTECT(nlist({{"rates", rates}, {"cmp", cmp}, {"ref", rf}, {"per", pe},
                            {"level", lv}}));
  Rf_classgets(out, PROTECT(Rf_mkString("dastd")));
  UNPROTECT(np + 2);
  return out;
}

SEXP C_clab_r_sirex(SEXP x, SEXP n, SEXP strat, SEXP grp, SEXP rates, SEXP level) {
  clab::who = "sirex";
  SEXP uv[4] = {x, n, strat, grp};
  static const char *const un[] = {"x", "n", "strat", "grp"};
  int np = csp::byname(uv, un, 4);
  SEXP xv = PROTECT(clab::num(uv[0], "x", -1, 0.0));
  const int m = len(xv, "x");
  SEXP nv = PROTECT(clab::num(uv[1], "n", m, 0.0, true));
  np += 2;
  csp::lab s = csp::codes(uv[2], "strat", m);
  if (Rf_isNull(uv[3])) fail("`grp` must give the group (area) of each row.");
  csp::lab g = csp::codes(uv[3], "grp", m);
  np += 4;
  std::vector<double> sr;
  if (!Rf_isNull(rates)) sr = perstr(rates, "rates", s.lev, s.g);
  const double a = 1 - REAL_RO(PROTECT(clab::dbl(level)))[0];
  ++np;
  (void)zq(level);
  const int k = g.g;
  const std::size_t sk = static_cast<std::size_t>(k);
  std::vector<double> o(5 * sk), wk(2 * static_cast<std::size_t>(s.g));
  double prior[2];
  const int rc = clab_sirex(REAL(xv), REAL(nv), INTEGER(g.code), INTEGER(s.code), m, k, s.g,
                            sr.empty() ? nullptr : sr.data(), &o[0], &o[sk], &o[2 * sk],
                            &o[3 * sk], &o[4 * sk], prior, wk.data());
  if (rc != CLAB_OK) {
    fail("every group needs a positive expected count: a stratum rate of 0 or a group "
         "with no exposure leaves none.");
  }
  cols cs(11);
  cs.put("group", g.lev);
  double *co = cs.dbl("obs", k), *ce = cs.dbl("expd", k), *ci = cs.dbl("sir", k),
         *cl = cs.dbl("lo", k), *ch = cs.dbl("hi", k), *cb = cs.dbl("eb", k),
         *cw = cs.dbl("shr", k), *fl = cs.dbl("flo", k), *fh = cs.dbl("fhi", k);
  int *cf = cs.lgc("out", k);
  for (int i = 0; i < k; ++i) {
    const double ob = o[static_cast<std::size_t>(i)], e = o[sk + static_cast<std::size_t>(i)];
    co[i] = ob;
    ce[i] = e;
    ci[i] = o[2 * sk + static_cast<std::size_t>(i)];
    cl[i] = gmlo(ob, a) / e;
    ch[i] = gmhi(ob, a) / e;
    cb[i] = o[3 * sk + static_cast<std::size_t>(i)];
    cw[i] = o[4 * sk + static_cast<std::size_t>(i)];
    // funnel limits for a ratio of 1 at this expected count
    fl[i] = Rf_qpois(a / 2, e, 1, 0) / e;
    fh[i] = Rf_qpois(1 - a / 2, e, 1, 0) / e;
    cf[i] = ci[i] < fl[i] || ci[i] > fh[i] ? 1 : 0;
  }
  SEXP tab = PROTECT(cs.frame());
  const char *const pn[] = {"mean", "var"};
  SEXP pv = PROTECT(clab::nvec(2, pn, prior));
  SEXP lv = PROTECT(Rf_ScalarReal(1 - a));
  np += 3;
  SEXP out = PROTECT(nlist({{"sir", tab}, {"prior", pv}, {"level", lv}}));
  Rf_classgets(out, PROTECT(Rf_mkString("sirex")));
  UNPROTECT(np + 2);
  return out;
}

//' @srrstats {G2.14} ctabx() states its rule for missing values: none
//'   are accepted; a table with an empty row or column has it dropped
//'   (counted in `drop`) before the statistics.
//' @noRd
SEXP C_clab_r_ctabx(SEXP row, SEXP col, SEXP w, SEXP agg, SEXP correct, SEXP exact,
                    SEXP level) {
  clab::who = "ctabx";
  int np = 0;
  const bool yates = csp::flag(correct, "correct"), ex = csp::flag(exact, "exact");
  const double z = zq(level);
  (void)z;
  const double conf = REAL_RO(PROTECT(clab::dbl(level)))[0];
  ++np;
  int r = 0, c = 0, ag = 0;
  std::vector<double> t;
  SEXP rl, cl;
  if (Rf_isNull(col)) {
    // a table given whole
    if (!clab::isnum(row) || !Rf_isMatrix(row)) {
      fail("`row` must be the row labels (with `col`) or a numeric matrix or 2-d table.");
    }
    clab::unlab(row, "row");
    if (!Rf_isNull(w)) fail("`w` goes with `row` and `col` labels, not with a table.");
    r = Rf_nrows(row);
    c = Rf_ncols(row);
    SEXP v = PROTECT(clab::dbl(row));
    ++np;
    t.assign(REAL(v), REAL(v) + static_cast<std::ptrdiff_t>(r) * c);
    for (double e : t) {
      if (!(e >= 0) || !std::isfinite(e)) fail("`row` must hold counts: finite numbers from 0.");
    }
    SEXP dn = Rf_getAttrib(row, R_DimNamesSymbol);
    rl = !Rf_isNull(dn) && !Rf_isNull(VECTOR_ELT(dn, 0)) ? VECTOR_ELT(dn, 0) : seqlab(r);
    PROTECT(rl);
    cl = !Rf_isNull(dn) && !Rf_isNull(VECTOR_ELT(dn, 1)) ? VECTOR_ELT(dn, 1) : seqlab(c);
    PROTECT(cl);
    np += 2;
    ag = Rf_isNull(agg) ? 0 : csp::pick(agg, "agg", aggs, 6);
    if (ag > 1) fail("`agg` applies to `row` and `col` labels, not to a table.");
  } else {
    SEXP uv[3] = {row, col, w};
    static const char *const un[] = {"row", "col", "w"};
    np += csp::byname(uv, un, 3);
    const int n = clab::ilen(static_cast<double>(Rf_xlength(uv[0])), "row");
    if (n < 1) fail("`row` must have at least one label.");
    csp::lab a = csp::codes(uv[0], "row", n);
    csp::lab b = csp::codes(uv[1], "col", n);
    np += 4;
    r = a.g;
    c = b.g;
    rl = a.lev;
    cl = b.lev;
    ag = Rf_isNull(agg) ? (Rf_isNull(uv[2]) ? 0 : 1) : csp::pick(agg, "agg", aggs, 6);
    const std::size_t rc = static_cast<std::size_t>(r) * static_cast<std::size_t>(c);
    t.assign(rc, 0.0);
    std::vector<double> cnt(rc, 0.0);
    const int *ia = INTEGER(a.code), *ib = INTEGER(b.code);
    if (ag == 0) {
      if (!Rf_isNull(uv[2])) fail("`agg = \"n\"` counts rows: leave `w` out.");
      for (int i = 0; i < n; ++i) t[static_cast<std::size_t>(ib[i]) * static_cast<std::size_t>(r) + static_cast<std::size_t>(ia[i])] += 1;
    } else if (Rf_isNull(uv[2])) {
      fail("`agg = \"%s\"` needs `w`, the values to aggregate.", aggs[ag]);
    } else if (ag == 5) {
      // distinct values of w (an id) in each cell
      csp::lab d = csp::codes(uv[2], "w", n);
      np += 2;
      std::vector<std::vector<int>> seen(rc);
      for (int i = 0; i < n; ++i) {
        seen[static_cast<std::size_t>(ib[i]) * static_cast<std::size_t>(r) + static_cast<std::size_t>(ia[i])].push_back(INTEGER(d.code)[i]);
      }
      for (std::size_t q = 0; q < rc; ++q) {
        std::sort(seen[q].begin(), seen[q].end());
        t[q] = static_cast<double>(std::unique(seen[q].begin(), seen[q].end()) - seen[q].begin());
      }
    } else {
      SEXP wv = PROTECT(clab::num(uv[2], "w", n, ag == 1 ? 0.0 : clab::ninf));
      ++np;
      for (int i = 0; i < n; ++i) {
        const std::size_t q = static_cast<std::size_t>(ib[i]) * static_cast<std::size_t>(r) + static_cast<std::size_t>(ia[i]);
        const double v = REAL(wv)[i];
        if (ag == 1 || ag == 4) t[q] += v;
        else if (cnt[q] == 0) t[q] = v;
        else if (ag == 2) t[q] = std::fmax(t[q], v);
        else t[q] = std::fmin(t[q], v);
        cnt[q] += 1;
      }
      if (ag == 4) {
        for (std::size_t q = 0; q < rc; ++q) t[q] = cnt[q] > 0 ? t[q] / cnt[q] : 0.0;
      }
    }
  }
  // the counts with the empty rows and columns dropped
  std::vector<int> kr, kc;
  for (int i = 0; i < r; ++i) {
    double s = 0;
    for (int j = 0; j < c; ++j) s += t[static_cast<std::size_t>(j) * static_cast<std::size_t>(r) + static_cast<std::size_t>(i)];
    if (s > 0) kr.push_back(i);
  }
  for (int j = 0; j < c; ++j) {
    double s = 0;
    for (int i = 0; i < r; ++i) s += t[static_cast<std::size_t>(j) * static_cast<std::size_t>(r) + static_cast<std::size_t>(i)];
    if (s > 0) kc.push_back(j);
  }
  const int r2 = static_cast<int>(kr.size()), c2 = static_cast<int>(kc.size());
  const bool counts = ag == 0 || ag == 1 || ag == 5;
  std::vector<double> u(static_cast<std::size_t>(r2) * static_cast<std::size_t>(c2)), e(u.size()),
      rs(u.size()), st(9, qnan);
  std::vector<int> dg(static_cast<std::size_t>(r2), -1);
  for (int j = 0; j < c2; ++j) {
    for (int i = 0; i < r2; ++i) {
      u[static_cast<std::size_t>(j) * static_cast<std::size_t>(r2) + static_cast<std::size_t>(i)] =
          t[static_cast<std::size_t>(kc[static_cast<std::size_t>(j)]) * static_cast<std::size_t>(r) +
            static_cast<std::size_t>(kr[static_cast<std::size_t>(i)])];
    }
  }
  for (int i = 0; i < r2; ++i) {
    for (int j = 0; j < c2; ++j) {
      if (clab::seql(STRING_ELT(rl, kr[static_cast<std::size_t>(i)]), STRING_ELT(cl, kc[static_cast<std::size_t>(j)]))) {
        dg[static_cast<std::size_t>(i)] = j;
      }
    }
  }
  bool any = false;
  for (int i = 0; i < r2; ++i) any = any || dg[static_cast<std::size_t>(i)] >= 0;
  double fs[5] = {qnan, qnan, qnan, qnan, qnan};
  const bool ok = counts && r2 >= 1 && c2 >= 1;
  if (ok) {
    kok(clab_ctabx(u.data(), r2, c2, yates ? 1 : 0, any ? dg.data() : nullptr, e.data(),
                   rs.data(), st.data()),
        "the kernel refused the table.");
    if (ex && r2 >= 2 && c2 >= 2) {
      std::vector<double> rem(static_cast<std::size_t>(r2 + c2));
      const int rc = clab_fshx(u.data(), r2, c2, 1e7, conf, fs, rem.data());
      if (rc != CLAB_OK) {
        for (double v : u) {
          if (v != std::floor(v)) fail("the exact test needs whole counts.");
        }
        fail("the exact test would visit more than 1e7 tables: leave `exact` FALSE.");
      }
    }
  } else if (ex) {
    fail("the exact test needs counts: `agg` \"n\", \"sum\" or \"dist\".");
  }
  const bool two = ok && r2 >= 2 && c2 >= 2;
  SEXP rl2 = PROTECT(Rf_allocVector(STRSXP, r2)), cl2 = PROTECT(Rf_allocVector(STRSXP, c2));
  np += 2;
  for (int i = 0; i < r2; ++i) SET_STRING_ELT(rl2, i, STRING_ELT(rl, kr[static_cast<std::size_t>(i)]));
  for (int j = 0; j < c2; ++j) SET_STRING_ELT(cl2, j, STRING_ELT(cl, kc[static_cast<std::size_t>(j)]));
  SEXP tab = PROTECT(labmat(t.data(), r, c, rl, cl));
  ++np;
  for (double &v : e) v = ok ? v : NA_REAL;
  for (double &v : rs) v = ok ? v : NA_REAL;
  SEXP em = PROTECT(labmat(e.data(), r2, c2, rl2, cl2));
  SEXP rm = PROTECT(labmat(rs.data(), r2, c2, rl2, cl2));
  np += 2;
  cols sc(17);
  auto put = [&](const char *nm, double v) { sc.dbl(nm, 1)[0] = std::isnan(v) ? NA_REAL : v; };
  const double df = st[2], ch = st[1], nn = st[0];
  const int kk = r2 < c2 ? r2 : c2;
  put("n", st[0]);
  put("chisq", two ? ch : qnan);
  put("df", two ? df : qnan);
  put("p", two ? Rf_pchisq(ch, df, 0, 0) : qnan);
  put("g2", two ? st[3] : qnan);
  put("g2p", two ? Rf_pchisq(st[3], df, 0, 0) : qnan);
  put("v", two ? st[4] : qnan);
  // DescTools' noncentral chi-square interval: sqrt(ncp / (n (k - 1)))
  const double chu = two ? (yates && r2 == 2 && c2 == 2 ? st[4] * st[4] * nn * (kk - 1) : ch) : qnan;
  put("vlo", two ? std::fmin(1.0, std::sqrt(ncp(chu, df, 1 - (1 - conf) / 2) / (nn * (kk - 1)))) : qnan);
  put("vhi", two ? std::fmin(1.0, std::sqrt(ncp(chu, df, (1 - conf) / 2) / (nn * (kk - 1)))) : qnan);
  put("vbc", two ? st[5] : qnan);
  put("cc", two ? st[6] : qnan);
  put("phi", st[7]);
  put("diag", st[8]);
  put("fp", fs[0]);
  put("or", fs[1]);
  put("orlo", fs[2]);
  put("orhi", fs[3]);
  SEXP stat = PROTECT(sc.frame());
  SEXP dr = PROTECT(Rf_ScalarInteger((r - r2) + (c - c2)));
  SEXP am = PROTECT(Rf_mkString(aggs[ag]));
  SEXP cr = PROTECT(Rf_ScalarLogical(yates ? 1 : 0));
  SEXP lv = PROTECT(Rf_ScalarReal(conf));
  np += 5;
  SEXP out = PROTECT(nlist({{"tab", tab}, {"expd", em}, {"res", rm}, {"stat", stat},
                            {"drop", dr}, {"agg", am}, {"correct", cr}, {"level", lv}}));
  Rf_classgets(out, PROTECT(Rf_mkString("ctabx")));
  UNPROTECT(np + 2);
  return out;
}

}  // extern "C"
