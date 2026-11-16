// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the environmental-crime family: the LMDI
// decomposition of a count over the factors of a Kaya identity (lmdix),
// the two-way fixed-effects regression of a count on continuous drivers
// with an interaction (fepnl), Pesaran's CD test (pcdtst), and least
// squares with its diagnostics (olsdx). Each checks
// its arguments, builds the kernel's arrays, calls the kernels (env.cpp,
// and did.cpp for the regression) and builds the R result. Nothing calls
// back into the R evaluator.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

#include "crimlab.h"
#include "bmat.h"
#include "pevr.h"
#include "rchk.h"
#include "spin.h"

extern "C" double Rf_pt(double, double, int, int);
extern "C" double Rf_qt(double, double, int, int);

using clab::fail;
using pevr::cols;
using bmat::mat;
using bmat::readm;
using pevr::nlist;
using pevr::z;

namespace {

// The periods at: NULL for def (0 the first and last, 1 all), else
// numbers that are periods of the panel, increasing; their codes.
std::vector<int> periods(SEXP at, const pevr::pan &pn, int def, int need) {
  std::vector<int> c;
  if (Rf_isNull(at)) {
    if (def == 0) {
      c.push_back(0);
      c.push_back(pn.nt - 1);
    } else {
      for (int s = 0; s < pn.nt; ++s) c.push_back(s);
    }
    return c;
  }
  SEXP a = PROTECT(clab::num(at, "at"));
  for (R_xlen_t q = 0; q < XLENGTH(a); ++q) {
    int s = -1;
    for (int r = 0; r < pn.nt; ++r) {
      if (pn.tv[z(r)] == REAL(a)[q]) s = r;
    }
    if (s < 0) fail("`at` must hold periods of `time`.");
    if (!c.empty() && s <= c.back()) fail("`at` must be increasing.");
    c.push_back(s);
  }
  UNPROTECT(1);
  if (static_cast<int>(c.size()) < need) fail("`at` must hold at least %d periods.", need);
  return c;
}

// The factors of row values y (1) and q (m): y / q1, q1 / q2, ..., qm.
void factors(double y, const double *q, int m, std::size_t step, double *f) {
  f[0] = y / q[0];
  for (int j = 1; j < m; ++j) f[z(j)] = q[z(j - 1) * step] / q[z(j) * step];
  f[z(m)] = q[z(m - 1) * step];
}

// The factor names: given (k of them) or from the quantities' names.
std::vector<std::string> fnames(SEXP fnm, const mat &q) {
  const int k = q.p + 1;
  std::vector<std::string> nm;
  if (!Rf_isNull(fnm)) {
    if (TYPEOF(fnm) != STRSXP || XLENGTH(fnm) != k) {
      fail("`fnm` must be %d names, one per factor.", k);
    }
    for (int j = 0; j < k; ++j) {
      if (STRING_ELT(fnm, j) == NA_STRING) fail("`fnm` must not hold NA.");
      nm.emplace_back(Rf_translateCharUTF8(STRING_ELT(fnm, j)));
    }
    return nm;
  }
  nm.push_back("y/" + q.nm[0]);
  for (int j = 1; j < q.p; ++j) nm.push_back(q.nm[z(j - 1)] + "/" + q.nm[z(j)]);
  nm.push_back(q.nm[z(q.p - 1)]);
  return nm;
}

// One decomposition: base (y0, F0) against target (y1, F1) per row.
struct dec {
  std::vector<int> unit, s0, s1;  // unit codes; base and target period codes
  std::vector<double> y0, y1, dy, f0, f1, c;
};

bool good(double v) { return std::isfinite(v) && v > 0; }

// The contributions of the rows whose values are all positive; NA in
// the others (a zero left by `zero = 0`).
void decompose(dec &d, int k) {
  const int r = static_cast<int>(d.unit.size());
  d.c.assign(z(r) * z(k), NA_REAL);
  std::vector<int> ok;
  for (int i = 0; i < r; ++i) {
    bool g = good(d.y0[z(i)]) && good(d.y1[z(i)]);
    for (int j = 0; j < k; ++j) {
      g = g && good(d.f0[z(i) * z(k) + z(j)]) && good(d.f1[z(i) * z(k) + z(j)]);
    }
    if (g) ok.push_back(i);
  }
  const int m = static_cast<int>(ok.size());
  if (m == 0) return;
  // f0 and f1 are held row by row; the kernel wants columns
  std::vector<double> y0(z(m)), y1(z(m)), a(z(m) * z(k)), b(z(m) * z(k)), c(z(m) * z(k));
  for (int h = 0; h < m; ++h) {
    const std::size_t i = z(ok[z(h)]);
    y0[z(h)] = d.y0[i];
    y1[z(h)] = d.y1[i];
    for (int j = 0; j < k; ++j) {
      a[z(j) * z(m) + z(h)] = d.f0[i * z(k) + z(j)];
      b[z(j) * z(m) + z(h)] = d.f1[i * z(k) + z(j)];
    }
  }
  clab_lmdi(y0.data(), y1.data(), a.data(), b.data(), m, k, c.data());
  for (int h = 0; h < m; ++h) {
    for (int j = 0; j < k; ++j) d.c[z(j) * z(r) + z(ok[z(h)])] = c[z(j) * z(m) + z(h)];
  }
}

// The rows of the decomposition as a data frame.
SEXP decframe(const dec &d, const pevr::pan &pn, const std::vector<std::string> &fn, bool sp) {
  const int r = static_cast<int>(d.unit.size()), k = static_cast<int>(fn.size());
  cols c(5 + k + (sp ? 0 : 1));
  SEXP cu = c.str("unit", r);
  double *a = c.dbl(sp ? "time" : "from", r), *b = sp ? nullptr : c.dbl("to", r);
  double *v0 = c.dbl(sp ? "ref" : "y0", r), *v1 = c.dbl(sp ? "y" : "y1", r);
  double *dy = c.dbl("dy", r);
  std::vector<double *> cf;
  for (int j = 0; j < k; ++j) cf.push_back(c.dbl(fn[z(j)].c_str(), r));
  for (int i = 0; i < r; ++i) {
    SET_STRING_ELT(cu, i, STRING_ELT(pn.ulev, d.unit[z(i)]));
    a[i] = pn.tv[z(d.s0[z(i)])];
    if (b) b[i] = pn.tv[z(d.s1[z(i)])];
    v0[i] = d.y0[z(i)];
    v1[i] = d.y1[z(i)];
    dy[i] = d.dy[z(i)];
    for (int j = 0; j < k; ++j) cf[z(j)][i] = d.c[z(j) * z(r) + z(i)];
  }
  return c.frame();
}

// The unit-period values of y and the factors, as given (zeros kept).
SEXP facframe(const pevr::pan &pn, const mat &q, const std::vector<std::string> &fn) {
  const int n = pn.n, k = q.p + 1;
  cols c(3 + k);
  SEXP cu = c.str("unit", n);
  double *t = c.dbl("time", n), *y = c.dbl("y", n);
  std::vector<double *> cf;
  for (int j = 0; j < k; ++j) cf.push_back(c.dbl(fn[z(j)].c_str(), n));
  std::vector<double> f(z(k));
  for (int i = 0; i < n; ++i) {
    SET_STRING_ELT(cu, i, STRING_ELT(pn.ulev, pn.u[z(i)]));
    t[i] = pn.tv[z(pn.t[z(i)])];
    y[i] = pn.y[i];
    factors(pn.y[i], q.v.data() + i, q.p, z(n), f.data());
    for (int j = 0; j < k; ++j) cf[z(j)][i] = pevr::na(f[z(j)]);
  }
  return c.frame();
}

// The reference region of period s: from the means (gm 0) or geometric
// means (gm 1) of y and of each quantity over the units present, with
// zeros replaced by zr. Out: ry, rf (k).
void refer(const pevr::pan &pn, const mat &q, const std::vector<int> &rows, int gm, double zr,
           double *ry, double *rf) {
  const int m = q.p, n = pn.n;
  std::vector<double> a(z(m) + 1, 0.0);
  for (int i : rows) {
    for (int j = 0; j <= m; ++j) {
      double v = j == 0 ? pn.y[i] : q.v[z(j - 1) * z(n) + z(i)];
      if (gm) v = std::log(v > 0 ? v : zr);  // -Inf, so a 0 mean, when zr is 0
      a[z(j)] += v / static_cast<double>(rows.size());
    }
  }
  for (int j = 0; j <= m; ++j) {
    if (gm) a[z(j)] = std::exp(a[z(j)]);
    if (!gm && a[z(j)] <= 0) a[z(j)] = zr;
  }
  *ry = a[0];
  factors(a[0], a.data() + 1, m, 1, rf);
}

}  // namespace

extern "C" {

//' @srrstats {G2.6} lmdix(), fepnl() and pcdtst() read the outcome, units
//'   and periods as plain vectors of one length; quantities and regressors
//'   as a numeric matrix with one row per row.
//' @noRd
SEXP C_clab_r_lmdix(SEXP y, SEXP q, SEXP unit, SEXP time, SEXP type, SEXP at, SEXP ref,
                    SEXP zero, SEXP fnm) {
  clab::who = "lmdix";
  static const char *const tyo[] = {"time", "space"};
  static const char *const rfo[] = {"mean", "gmean"};
  const int sp = csp::pick(type, "type", tyo, 2);
  const int gm = csp::pick(ref, "ref", rfo, 2);
  const double zr = clab::one(zero, "zero", -1.0, R_PosInf);
  if (zr < 0) fail("`zero` must be 0 (leave zeros out) or a positive number.");
  int np = 0;
  pevr::pan pn = pevr::panel(y, unit, time, 0.0, &np);
  const int n = pn.n;
  const mat qm = readm(q, "q", "q", n, 0.0);
  const int m = qm.p, k = m + 1;
  const std::vector<std::string> fn = fnames(fnm, qm);
  // the factors of every row with the zeros replaced
  std::vector<double> yz(z(n)), fz(z(n) * z(k)), qz(qm.v);
  for (double &v : qz) v = v > 0 ? v : zr;  // 0 stays 0 when zr is 0
  for (int i = 0; i < n; ++i) {
    yz[z(i)] = pn.y[i] > 0 ? pn.y[i] : zr;
    factors(yz[z(i)], qz.data() + i, m, z(n), fz.data() + z(i) * z(k));
  }
  const double cells = static_cast<double>(pn.nu) * pn.nt;
  if (cells > 67108864.0) fail("the units by periods (%.0f cells) are too many.", cells);
  std::vector<int> row(z(pn.nu) * z(pn.nt), -1);
  for (int i = 0; i < n; ++i) row[z(pn.u[z(i)]) * z(pn.nt) + z(pn.t[z(i)])] = i;
  const std::vector<int> ps = periods(at, pn, sp, sp ? 1 : 2);
  dec d;
  SEXP rf = R_NilValue;
  auto add = [&](int a, int s0, int s1, double y0, double y1, double dy, const double *f0,
                 const double *f1) {
    d.unit.push_back(a);
    d.s0.push_back(s0);
    d.s1.push_back(s1);
    d.y0.push_back(y0);
    d.y1.push_back(y1);
    d.dy.push_back(dy);
    d.f0.insert(d.f0.end(), f0, f0 + k);
    d.f1.insert(d.f1.end(), f1, f1 + k);
  };
  if (!sp) {
    for (std::size_t h = 1; h < ps.size(); ++h) {
      for (int a = 0; a < pn.nu; ++a) {
        const int i0 = row[z(a) * z(pn.nt) + z(ps[h - 1])], i1 = row[z(a) * z(pn.nt) + z(ps[h])];
        if (i0 < 0 || i1 < 0) continue;
        add(a, ps[h - 1], ps[h], pn.y[i0], pn.y[i1], pn.y[i1] - pn.y[i0],
            fz.data() + z(i0) * z(k), fz.data() + z(i1) * z(k));
      }
    }
  } else {
    cols rc(2 + k);
    const int nr = static_cast<int>(ps.size());
    double *rt = rc.dbl("time", nr), *ry = rc.dbl("y", nr);
    std::vector<double *> rcf;
    for (int j = 0; j < k; ++j) rcf.push_back(rc.dbl(fn[z(j)].c_str(), nr));
    std::vector<double> f(z(k));
    for (int h = 0; h < nr; ++h) {
      std::vector<int> rows;
      for (int a = 0; a < pn.nu; ++a) {
        const int i = row[z(a) * z(pn.nt) + z(ps[z(h)])];
        if (i >= 0) rows.push_back(i);
      }
      double yr = 0;
      refer(pn, qm, rows, gm, zr, &yr, f.data());
      rt[h] = pn.tv[z(ps[z(h)])];
      ry[h] = yr;
      for (int j = 0; j < k; ++j) rcf[z(j)][h] = f[z(j)];
      for (int i : rows) {
        add(pn.u[z(i)], ps[z(h)], ps[z(h)], yr, pn.y[i], pn.y[i] - yr, f.data(),
            fz.data() + z(i) * z(k));
      }
    }
    rf = PROTECT(rc.frame());
    ++np;
  }
  // the kernel needs positive values: y as the factors saw it
  std::vector<double> y0(d.y0), y1(d.y1);
  for (double &v : d.y0) v = v > 0 ? v : zr;
  for (double &v : d.y1) v = v > 0 ? v : zr;  // as the factors were made
  decompose(d, k);
  d.y0.swap(y0);
  d.y1.swap(y1);
  SEXP de = PROTECT(decframe(d, pn, fn, sp));
  SEXP fa = PROTECT(facframe(pn, qm, fn));
  SEXP fnv = PROTECT(Rf_allocVector(STRSXP, k));
  for (int j = 0; j < k; ++j) SET_STRING_ELT(fnv, j, Rf_mkCharCE(fn[z(j)].c_str(), CE_UTF8));
  SEXP ty = PROTECT(Rf_mkString(tyo[sp])), rs = PROTECT(Rf_mkString(rfo[gm]));
  np += 5;
  SEXP out = PROTECT(nlist({{"dec", de}, {"fac", fa}, {"ref", rf}, {"factors", fnv},
                            {"type", ty}, {"refkind", rs}}));
  Rf_classgets(out, PROTECT(Rf_mkString("lmdix")));
  UNPROTECT(np + 2);
  return out;
}

}  // extern "C"

namespace {

// The cluster of each row: by unit (clus NULL) or the labels given.
int clust(SEXP clus, const pevr::pan &pn, std::vector<int> &c) {
  if (Rf_isNull(clus)) {
    c = pn.u;
    return pn.nu;
  }
  csp::lab l = csp::codes(clus, "clus", pn.n);
  c.assign(INTEGER(l.code), INTEGER(l.code) + pn.n);
  const int g = l.g;
  UNPROTECT(2);
  if (g < 2) fail("`clus` has one cluster: a cluster-robust variance needs two or more.");
  return g;
}

// The column of x named by element j of int, or a refusal.
int which(const mat &x, SEXP in, int j) {
  const char *s = Rf_translateCharUTF8(STRING_ELT(in, j));
  for (int h = 0; h < x.p; ++h) {
    if (x.nm[z(h)] == s) return h;
  }
  fail("`int` names %s, which is not a column of `x`.", s);
}

double sdev(const double *v, int n) {
  double m = 0, s = 0;
  for (int i = 0; i < n; ++i) m += v[i] / n;
  for (int i = 0; i < n; ++i) s += (v[i] - m) * (v[i] - m);
  return std::sqrt(s / (n - 1));
}

double mean(const double *v, int n) {
  double m = 0;
  for (int i = 0; i < n; ++i) m += v[i] / n;
  return m;
}

// The simple slopes of column a at values of the moderator b, through
// the interaction column c: est b_a + b_c m, se from the variance.
SEXP slopes(SEXP at, const std::vector<double> &mv, int a, int c, const std::vector<double> &b,
            const std::vector<double> &vc, int p, double df, double lv) {
  std::vector<double> m;
  const int n = static_cast<int>(mv.size());
  if (Rf_isNull(at)) {
    const double mu = mean(mv.data(), n), s = sdev(mv.data(), n);
    m = {mu - s, mu, mu + s};
  } else {
    SEXP v = PROTECT(clab::num(at, "at"));
    m.assign(REAL(v), REAL(v) + XLENGTH(v));
    UNPROTECT(1);
  }
  const int r = static_cast<int>(m.size());
  const double tq = Rf_qt(1 - (1 - lv) / 2, df, 1, 0);
  cols o(7);
  double *oa = o.dbl("at", r), *oe = o.dbl("slope", r), *os = o.dbl("se", r),
         *ol = o.dbl("lo", r), *oh = o.dbl("hi", r), *ot = o.dbl("t", r), *op = o.dbl("p", r);
  const std::size_t sp = z(p);
  for (int h = 0; h < r; ++h) {
    const double x = m[z(h)];
    const double e = b[z(a)] + b[z(c)] * x;
    const double v = vc[z(a) * sp + z(a)] + x * x * vc[z(c) * sp + z(c)] +
                     2 * x * vc[z(a) * sp + z(c)];
    oa[h] = x;
    oe[h] = e;
    os[h] = std::sqrt(v);
    ol[h] = e - tq * os[h];
    oh[h] = e + tq * os[h];
    ot[h] = e / os[h];
    op[h] = 2 * Rf_pt(-std::fabs(ot[h]), df, 1, 0);
  }
  return o.frame();
}

}  // namespace

extern "C" {

SEXP C_clab_r_fepnl(SEXP y, SEXP x, SEXP unit, SEXP time, SEXP in, SEXP clus, SEXP at,
                    SEXP ssc, SEXP level) {
  clab::who = "fepnl";
  static const char *const sso[] = {"cr1", "none"};
  const int raw = csp::pick(ssc, "ssc", sso, 2);
  const double lv = clab::one(level, "level", 0.0, 1.0);
  int np = 0;
  pevr::pan pn = pevr::panel(y, unit, time, clab::ninf, &np);
  const int n = pn.n;
  mat xm = readm(x, "x", "x", n, clab::ninf);
  int ia = -1, ib = -1;
  if (!Rf_isNull(in)) {
    if (TYPEOF(in) != STRSXP || XLENGTH(in) != 2 || STRING_ELT(in, 0) == NA_STRING ||
        STRING_ELT(in, 1) == NA_STRING) {
      fail("`int` must be two names of columns of `x`.");
    }
    ia = which(xm, in, 0);
    ib = which(xm, in, 1);
    if (ia == ib) fail("`int` must name two different columns of `x`.");
    for (int i = 0; i < n; ++i) {
      xm.v.push_back(xm.v[z(ia) * z(n) + z(i)] * xm.v[z(ib) * z(n) + z(i)]);
    }
    xm.nm.push_back(xm.nm[z(ia)] + ":" + xm.nm[z(ib)]);
    ++xm.p;
  } else if (!Rf_isNull(at)) {
    fail("`at` gives moderator values for the simple slopes: give `int` too.");
  }
  const int p = xm.p;
  std::vector<int> cl;
  const int ncl = clust(clus, pn, cl);
  if (ncl < 10) {
    Rf_warningcall(R_NilValue,
                   "%s: fewer than ten clusters; cluster-robust errors are unreliable with so "
                   "few (wcbx() is made for this case).",
                   clab::who);
  }
  const std::size_t sp = z(p), sn = z(n);
  std::vector<double> b(sp), vc(sp * sp), work(sn * (sp + 2) + sp * sp + 2 * z(pn.nu + pn.nt) +
                                                (sn + 2 * sp) * sp + sp);
  std::vector<int> iwork(sn + z(std::max(pn.nu, std::max(pn.nt, ncl))));
  double st[5];
  if (n <= p + pn.nu + pn.nt ||
      clab_didfe(pn.y, xm.v.data(), n, p, pn.u.data(), pn.nu, pn.t.data(), pn.nt, 1, nullptr,
                 cl.data(), ncl, nullptr, 0, 1e-10, 10000, b.data(), vc.data(), st, work.data(),
                 iwork.data()) != CLAB_OK) {
    fail("a column of `x` is absorbed by the unit and period effects, the columns are "
         "collinear, or there are too few rows for the effects.");
  }
  if (raw) {
    for (double &v : vc) v /= st[3];
  }
  // the within outcome and the residuals stay in work (did.cpp)
  const double *yw = work.data(), *res = work.data() + sn * (sp + 1);
  std::vector<double> vif(sp), wv(sn * (sp + 1) + 2 * sp * (sp + 1));
  if (clab_vif(xm.v.data(), n, p, vif.data(), wv.data()) != CLAB_OK) vif.assign(sp, NA_REAL);
  double jb[4];
  const bool jok = clab_jb(res, n, jb) == CLAB_OK;
  const double df = st[1] - 1, tq = Rf_qt(1 - (1 - lv) / 2, df, 1, 0);
  const double sy = sdev(pn.y, n);
  cols co(9);
  SEXP ct = co.str("term", p);
  double *ce = co.dbl("est", p), *cs = co.dbl("se", p), *cl0 = co.dbl("lo", p),
         *ch = co.dbl("hi", p), *cz = co.dbl("t", p), *cp = co.dbl("p", p),
         *cb = co.dbl("beta", p), *cv = co.dbl("vif", p);
  for (int j = 0; j < p; ++j) {
    SET_STRING_ELT(ct, j, Rf_mkCharCE(xm.nm[z(j)].c_str(), CE_UTF8));
    ce[j] = b[z(j)];
    cs[j] = std::sqrt(vc[z(j) * sp + z(j)]);
    cl0[j] = ce[j] - tq * cs[j];
    ch[j] = ce[j] + tq * cs[j];
    cz[j] = ce[j] / cs[j];
    cp[j] = 2 * Rf_pt(-std::fabs(cz[j]), df, 1, 0);
    cb[j] = ce[j] * sdev(xm.v.data() + z(j) * sn, n) / sy;
    cv[j] = pevr::na(vif[z(j)]);
  }
  SEXP coef = PROTECT(co.frame());
  ++np;
  SEXP sl = R_NilValue;
  if (ia >= 0) {
    std::vector<double> mv(xm.v.begin() + static_cast<std::ptrdiff_t>(z(ib) * sn),
                           xm.v.begin() + static_cast<std::ptrdiff_t>(z(ib + 1) * sn));
    sl = PROTECT(slopes(at, mv, ia, p - 1, b, vc, p, df, lv));
    ++np;
  }
  // R-squared of the levels and of the within fit
  double tss = 0, tw = 0;
  const double my = mean(pn.y, n);
  for (int i = 0; i < n; ++i) {
    tss += (pn.y[i] - my) * (pn.y[i] - my);
    tw += yw[i] * yw[i];
  }
  cols fs(8);
  double *fn = fs.dbl("n", 1), *fu = fs.dbl("units", 1), *ft = fs.dbl("periods", 1),
         *fg = fs.dbl("ncl", 1), *fk = fs.dbl("k", 1), *fr = fs.dbl("r2", 1),
         *fw = fs.dbl("r2w", 1), *frs = fs.dbl("rss", 1);
  fn[0] = n;
  fu[0] = pn.nu;
  ft[0] = pn.nt;
  fg[0] = st[1];
  fk[0] = st[2];
  fr[0] = pevr::na(1 - st[4] / tss);
  fw[0] = pevr::na(1 - st[4] / tw);
  frs[0] = st[4];
  SEXP info = PROTECT(fs.frame());
  SEXP vm = PROTECT(Rf_allocMatrix(REALSXP, p, p));
  SEXP rv = PROTECT(Rf_allocVector(REALSXP, n));
  np += 3;
  for (std::size_t q = 0; q < sp * sp; ++q) REAL(vm)[q] = vc[q];
  for (int i = 0; i < n; ++i) REAL(rv)[i] = res[i];
  static const char *const jn[] = {"skew", "kurt", "jb", "p"};
  const double na4[4] = {NA_REAL, NA_REAL, NA_REAL, NA_REAL};
  SEXP jv = PROTECT(clab::nvec(4, jn, jok ? jb : na4));
  SEXP sv = PROTECT(Rf_mkString(sso[raw])), lvs = PROTECT(Rf_ScalarReal(lv));
  np += 3;
  SEXP out = PROTECT(nlist({{"coef", coef}, {"slopes", sl}, {"fit", info}, {"jb", jv},
                            {"vcov", vm}, {"resid", rv}, {"ssc", sv}, {"level", lvs}}));
  Rf_classgets(out, PROTECT(Rf_mkString("fepnl")));
  UNPROTECT(np + 2);
  return out;
}

SEXP C_clab_r_pcdtst(SEXP x, SEXP unit, SEXP time) {
  clab::who = "pcdtst";
  int np = 0;
  pevr::pan pn = pevr::panel(x, unit, time, clab::ninf, &np);
  const double cells = static_cast<double>(pn.nu) * pn.nt;
  if (cells > 67108864.0) fail("the units by periods (%.0f cells) are too many.", cells);
  std::vector<double> e(z(pn.nu) * z(pn.nt), std::numeric_limits<double>::quiet_NaN());
  for (int i = 0; i < pn.n; ++i) e[z(pn.u[z(i)]) * z(pn.nt) + z(pn.t[z(i)])] = pn.y[i];
  double o[7];
  if (clab_pcd(e.data(), pn.nu, pn.nt, o) != CLAB_OK) {
    fail("no pair of units shares two periods with variation in both: CD is undefined.");
  }
  o[5] = pn.nu;
  o[6] = pn.nt;
  static const char *const nm[] = {"cd", "p", "rho", "absrho", "pairs", "units", "periods"};
  UNPROTECT(np);
  return clab::nvec(7, nm, o);
}

SEXP C_clab_r_olsdx(SEXP y, SEXP x, SEXP at, SEXP icp, SEXP pred, SEXP level) {
  clab::who = "olsdx";
  static const char *const po[] = {"mean", "new"};
  const int pn = csp::pick(pred, "pred", po, 2);
  const double lv = clab::one(level, "level", 0.0, 1.0);
  clab::lgl(icp, "int", 1);
  const int ic = LOGICAL(icp)[0] ? 1 : 0;
  int np = 0;
  SEXP yv = PROTECT(clab::num(y, "y"));
  ++np;
  const int n = static_cast<int>(XLENGTH(yv));
  mat xm = readm(x, "x", "x", n, clab::ninf);
  const int p = xm.p, k = p + ic;
  if (n < k + 2) {
    fail("`y` must have at least %d values: the %d coefficients and two more.", k + 2, k);
  }
  bool cn = true;
  for (int i = 1; i < n; ++i) cn = cn && REAL(yv)[i] == REAL(yv)[0];
  if (cn) fail("`y` is constant: there is nothing to explain.");
  int m = 0;
  std::vector<double> av;
  if (!Rf_isNull(at)) {
    if (!clab::isnum(at) || !Rf_isMatrix(at) || Rf_ncols(at) != p || Rf_nrows(at) < 1) {
      fail("`at` must be a numeric matrix with the %d columns of `x`, one row per point.", p);
    }
    m = Rf_nrows(at);
    mat am = readm(at, "at", "x", m, clab::ninf);
    SEXP dn = Rf_getAttrib(at, R_DimNamesSymbol);
    if (!Rf_isNull(dn) && !Rf_isNull(VECTOR_ELT(dn, 1)) && am.nm != xm.nm) {
      std::string s;
      for (const std::string &v : xm.nm) s += (s.empty() ? "" : ", ") + v;
      fail("`at` must name its columns as `x` does, in its order: %s.", s.c_str());
    }
    av = am.v;
  }
  const std::size_t sn = z(n), sp = z(p), sk = z(k), sm = z(m);
  std::vector<double> co(sk * 9), vc(sk * sk), cr(sp * sp), cs(sn * 7), pr(sm * 4);
  double st[10];
  if (clab_olsdx(REAL(yv), xm.v.data(), n, p, ic, m ? av.data() : nullptr, m, pn, lv, co.data(),
                 vc.data(), cr.data(), cs.data(), pr.data(), st) != CLAB_OK) {
    fail("the columns of `x` are collinear%s: least squares has no unique solution.",
         ic ? " (with the constant)" : "");
  }
  if (!(st[1] < 1.0)) {
    fail("the fit is exact (an R-squared of 1 to rounding): `y` is a linear function of `x`.");
  }
  static const char *const cn9[] = {"est", "se", "t", "p", "lo", "hi", "beta", "tol", "vif"};
  cols cc(10);
  SEXP ct = cc.str("term", k);
  if (ic) SET_STRING_ELT(ct, 0, Rf_mkChar("(Intercept)"));
  for (int j = 0; j < p; ++j) SET_STRING_ELT(ct, j + ic, Rf_mkCharCE(xm.nm[z(j)].c_str(), CE_UTF8));
  for (std::size_t c = 0; c < 9; ++c) {
    double *o = cc.dbl(cn9[c], k);
    for (std::size_t j = 0; j < sk; ++j) o[j] = c < 6 ? co[c * sk + j] : pevr::na(co[c * sk + j]);
  }
  SEXP coef = PROTECT(cc.frame());
  ++np;
  static const char *const fn[] = {"n", "k", "r", "r2", "adjr2", "sigma", "f",
                                   "df1", "df2", "p", "dw", "rss"};
  const double fv[12] = {static_cast<double>(n), static_cast<double>(k), st[0], st[1], st[2],
                         st[3], st[4], st[5], st[6], st[7], st[8], st[9]};
  cols fc(12);
  for (int q = 0; q < 12; ++q) *fc.dbl(fn[q], 1) = fv[q];
  SEXP fit = PROTECT(fc.frame());
  ++np;
  static const char *const cn7[] = {"fit", "resid", "zres", "sres", "tres", "hat", "cook"};
  cols ca(7);
  for (std::size_t c = 0; c < 7; ++c) {
    double *o = ca.dbl(cn7[c], n);
    for (std::size_t i = 0; i < sn; ++i) o[i] = cs[c * sn + i];
  }
  SEXP cases = PROTECT(ca.frame());
  ++np;
  SEXP cm = PROTECT(Rf_allocMatrix(REALSXP, p, p));
  SEXP vm = PROTECT(Rf_allocMatrix(REALSXP, k, k));
  np += 2;
  for (std::size_t q = 0; q < sp * sp; ++q) REAL(cm)[q] = cr[q];
  for (std::size_t q = 0; q < sk * sk; ++q) REAL(vm)[q] = vc[q];
  SEXP pf = R_NilValue;
  if (m) {
    static const char *const cn4[] = {"fit", "se", "lo", "hi"};
    cols pc(4);
    for (std::size_t c = 0; c < 4; ++c) {
      double *o = pc.dbl(cn4[c], m);
      for (std::size_t a = 0; a < sm; ++a) o[a] = pr[c * sm + a];
    }
    pf = PROTECT(pc.frame());
    ++np;
  }
  SEXP iv = PROTECT(Rf_ScalarLogical(ic)), ps = PROTECT(Rf_mkString(po[pn])),
       lvs = PROTECT(Rf_ScalarReal(lv));
  np += 3;
  SEXP out = PROTECT(nlist({{"coef", coef}, {"fit", fit}, {"cases", cases}, {"cor", cm},
                            {"vcov", vm}, {"pred", pf}, {"int", iv}, {"ptype", ps},
                            {"level", lvs}}));
  Rf_classgets(out, PROTECT(Rf_mkString("olsdx")));
  UNPROTECT(np + 2);
  return out;
}

}  // extern "C"
