// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of crimlab's count and table statistics: count
// regression (cntreg). Each checks its arguments, builds the kernel's
// arrays, calls the kernels and builds the R result. Distribution
// functions come from R's maths library. Nothing calls back into the R
// evaluator.
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
    SEXP f = clab::frame(v, p.data());
    UNPROTECT(1);
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

const char *cntmsg(int rc) {
  return rc == CLAB_EDOM ? "the fit failed: the deviance or a coefficient left the finite "
                           "range (a group or cell with no counts can do this)."
                         : "the kernel refused the input.";
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
    if (rc == CLAB_OK) {
      ok[f] = true;
    } else if (f == 0) {
      fail("%s", cntmsg(rc));
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

}  // extern "C"
