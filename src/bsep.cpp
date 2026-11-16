// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the spatial epidemiology methods
// (src/sep.cpp). Each reads its points through csp::xy() and its weights
// through csp::wts() (spin.h), checks the other arguments, calls the
// kernels and builds the R result. Nothing here calls back into the R
// evaluator.
#include <cmath>
#include <cstdio>
#include <cstring>

#include "bkh.h"
#include "crimlab.h"
#include "kutil.h"
#include "rchk.h"
#include "spin.h"

using namespace bkh;
using clab::num;

namespace {

// A design column that does not vary cannot be told from the intercept:
// the fit would be singular, or not, by the rounding of its sums.
void varies(SEXP m, int n, int p, const char *nm) {
  for (int j = 0; j < p && n > 1; ++j) {
    if (kutil::cnst(REAL(m) + static_cast<R_xlen_t>(j) * n, n)) {
      fail("column %d of `%s` is constant: it cannot be told from the intercept.", j + 1, nm);
    }
  }
}

// the region of each subject: whole numbers from 1 to nr
int *regions(SEXP reg, int n, int nr) {
  SEXP r = PROTECT(num(reg, "reg", n, 1.0));
  int *o = ial(n);
  for (int i = 0; i < n; ++i) {
    const double v = REAL(r)[i];
    if (v != std::floor(v) || v > nr) fail("`reg` must be area numbers from 1 to %d, the rows of `w`.", nr);
    o[i] = static_cast<int>(v) - 1;
  }
  UNPROTECT(1);
  return o;
}

// the times and events of a survival model, n of each
void survdat(SEXP time, SEXP event, int *n, SEXP *tv, SEXP *dv) {
  *tv = PROTECT(num(time, "time", -1, 0.0, true));
  *n = clab::ilen(static_cast<double>(XLENGTH(*tv)), "time");
  *dv = PROTECT(csp::bin01(event, "event", *n));
  SEXP d = PROTECT(real(*n));
  for (int i = 0; i < *n; ++i) REAL(d)[i] = INTEGER(*dv)[i];
  *dv = d;
}

// receptors and sources: planar, or longitude and latitude (raw)
void rsrc(SEXP x, SEXP y, SEXP src, SEXP lonlat, csp::pts *r, csp::pts *s, bool *ll) {
  *ll = csp::flag(lonlat, "lonlat");
  *r = csp::xy(x, y, R_NilValue, *ll ? csp::raw : csp::planar);
  *s = csp::xy(src, R_NilValue, R_NilValue, *ll ? csp::raw : csp::planar);
}

const double *swts(SEXP w, int m) {
  if (Rf_isNull(w)) return nullptr;
  return REAL(PROTECT(num(w, "w", m, 0.0)));
}

// the number of areas a weights argument w is for: its matrix's rows, its
// list's n, or a listw's neighbour lists
int nareas(SEXP w) {
  if (Rf_isMatrix(w)) return Rf_nrows(w);
  if (TYPEOF(w) == VECSXP) {
    SEXP nb = elt(w, "neighbours");
    if (Rf_inherits(w, "listw") && TYPEOF(nb) == VECSXP) return static_cast<int>(XLENGTH(nb));
    SEXP nn = elt(w, "n");
    const double v = clab::val1(nn);
    if (v >= 1 && v <= 1e9 && v == std::floor(v)) return static_cast<int>(v);
  }
  fail("`w` must be the list spwt() returns, an spdep listw or a matrix.");
}

// the frailty (kind 0) and cure (kind 1) models
SEXP surv(SEXP time, SEXP event, SEXP x, SEXP reg, SEXP w, SEXP tau, SEXP tol, int kind) {
  int n = 0, p = 0;
  SEXP tv, dv;
  survdat(time, event, &n, &tv, &dv);
  SEXP xm = zmat(x, "x", n, false, &p);
  const int nr = nareas(w);
  SEXP c = csp::wts(w, nr);
  const int *rg = regions(reg, n, nr);
  SEXP tu = PROTECT(num(tau, "tau", 1, 0.0));
  const double tl = clab::one(tol, "tol", 0.0, 1.0);
  const int k = (kind == 0 ? 1 : 2) + p + nr;
  double *th = dal(k), *se = dal(k), f = 0;
  int it = 0;
  SEXP cure = PROTECT(real(kind == 1 ? n : 0));
  const int *rp = INTEGER(VECTOR_ELT(c, 0)), *ci = INTEGER(VECTOR_ELT(c, 1));
  const double *wv = REAL(VECTOR_ELT(c, 2));
  const int rc = kind == 0 ? clab_sfrail(REAL(tv), REAL(dv), n, REAL(xm), p, rg, rp, ci, wv, nr, REAL(tu)[0], tl, th,
                                         se, &f, &it)
                           : clab_sfcure(REAL(tv), REAL(dv), n, REAL(xm), p, rg, rp, ci, wv, nr, REAL(tu)[0], tl, th,
                                         se, &f, REAL(cure), &it);
  kc(rc, "the Hessian is singular: the model is not identified by these data.");
  const int h = kind == 0 ? 1 : 2;
  SEXP cf = PROTECT(real(p)), u = PROTECT(real(nr)), sv = PROTECT(real(k));
  for (int j = 0; j < p; ++j) REAL(cf)[j] = th[h + j];
  for (int j = 0; j < nr; ++j) REAL(u)[j] = th[h + p + j];
  for (int j = 0; j < k; ++j) REAL(sv)[j] = se[j];
  SEXP r = PROTECT(Rf_allocVector(VECSXP, kind == 0 ? 7 : 9));
  SET_VECTOR_ELT(r, 0, Rf_ScalarReal(f));
  SET_VECTOR_ELT(r, 1, Rf_ScalarReal(th[0]));
  if (kind == 0) {
    SET_VECTOR_ELT(r, 2, cf);
    SET_VECTOR_ELT(r, 3, u);
    SET_VECTOR_ELT(r, 4, sv);
    SET_VECTOR_ELT(r, 5, Rf_ScalarInteger(it));
    SET_VECTOR_ELT(r, 6, Rf_ScalarReal(REAL(tu)[0]));
    static const char *const nm[] = {"pll", "shape", "coef", "frail", "se", "iter", "tau"};
    lst(r, nm, "sfrail");
  } else {
    SET_VECTOR_ELT(r, 2, Rf_ScalarReal(th[1]));
    SET_VECTOR_ELT(r, 3, cf);
    SET_VECTOR_ELT(r, 4, u);
    SET_VECTOR_ELT(r, 5, sv);
    SET_VECTOR_ELT(r, 6, cure);
    SET_VECTOR_ELT(r, 7, Rf_ScalarInteger(it));
    SET_VECTOR_ELT(r, 8, Rf_ScalarReal(REAL(tu)[0]));
    static const char *const nm[] = {"pll", "shape", "lscale", "coef", "effect", "se", "cure", "iter", "tau"};
    lst(r, nm, "sfcure");
  }
  UNPROTECT(11);
  return r;
}

}  // namespace

extern "C" {

// ----------------------------------------------------- spatial epidemiology

SEXP C_clab_r_bayout(SEXP x, SEXP freq, SEXP b, SEXP w, SEXP acty, SEXP alpha, SEXP at) {
  clab::who = "bayout";
  clab::unlab(x, "x");
  if (!clab::isnum(x) || clab::xlen(x, "x") < 1) fail("`x` must be a numeric vector of counts.");
  clab::vec(x, "x");
  const int n = clab::xlen(x, "x");
  SEXP xd = PROTECT(clab::dbl(x));
  SEXP xv = PROTECT(clab::dcopy(xd));
  for (int i = 0; i < n; ++i) {
    const double e = REAL(xv)[i];
    if (ISNAN(e)) {
      REAL(xv)[i] = NAN;
    } else if (!(e >= 0) || !std::isfinite(e)) {
      fail("`x` must be counts from 0, or NA.");
    }
  }
  const int fq = clab::whole(freq, "freq", 1), bb = clab::whole(b, "b", 0), ww = clab::whole(w, "w", 0);
  const bool ay = csp::flag(acty, "acty");
  const double al = clab::one(alpha, "alpha", 0.0, 1.0);
  const double need = static_cast<double>(bb) * fq + ww + 1;
  int m = 0;
  int *t = nullptr;
  if (Rf_isNull(at)) {
    if (n < need) fail("`x` has %d values; the reference window needs at least %.0f.", n, need);
    m = n - static_cast<int>(need) + 1;
    t = ial(m);
    for (int q = 0; q < m; ++q) t[q] = static_cast<int>(need) + q;
  } else {
    SEXP a = PROTECT(num(at, "at", -1, 1.0));
    m = clab::xlen(a, "at");
    t = ial(m);
    for (int q = 0; q < m; ++q) {
      const double e = REAL(a)[q];
      if (e != std::floor(e) || e > n) fail("`at` must be times from 1 to %d.", n);
      if (e < need) fail("time %.0f is too early for the reference window, which needs %.0f.", e, need);
      t[q] = static_cast<int>(e);
    }
    UNPROTECT(1);
  }
  SEXP c = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(c, 0, ints(m));
  SET_VECTOR_ELT(c, 1, real(m));
  SET_VECTOR_ELT(c, 2, Rf_allocVector(LGLSXP, m));
  int *al2 = ial(m);
  const int rc = clab_bayout(REAL(xv), n, fq, bb, ww, ay, al, t, m, REAL(VECTOR_ELT(c, 1)), al2);
  if (rc == CLAB_EARG) fail("the counts of a reference window are too large for the negative binomial bound.");
  kc(rc, "a reference window holds no observed count.");
  for (int q = 0; q < m; ++q) {
    INTEGER(VECTOR_ELT(c, 0))[q] = t[q];
    LOGICAL(VECTOR_ELT(c, 2))[q] = al2[q] < 0 ? NA_LOGICAL : al2[q];
  }
  static const char *const nm[] = {"t", "ub", "alarm"};
  SEXP f = PROTECT(clab::frame(c, nm));
  UNPROTECT(4);
  return f;
}



SEXP C_clab_r_bufexp(SEXP x, SEXP y, SEXP src, SEXP r, SEXP w, SEXP lonlat) {
  clab::who = "bufexp";
  csp::pts a, s;
  bool ll = false;
  rsrc(x, y, src, lonlat, &a, &s, &ll);
  SEXP rv = PROTECT(num(r, "r", 1, 0.0));
  const double *wv = swts(w, s.n);
  SEXP c = PROTECT(Rf_allocVector(VECSXP, 3));
  for (int j = 0; j < 3; ++j) SET_VECTOR_ELT(c, j, real(a.n));
  kc(clab_bufexp(a.x, a.y, a.n, s.x, s.y, s.n, REAL(rv)[0], wv, ll, REAL(VECTOR_ELT(c, 0)),
                 REAL(VECTOR_ELT(c, 1)), REAL(VECTOR_ELT(c, 2))),
     "the kernel refused the input.");
  static const char *const nm[] = {"count", "wsum", "near"};
  SEXP f = PROTECT(clab::frame(c, nm));
  UNPROTECT(5 + (wv != nullptr));
  return f;
}

SEXP C_clab_r_kerexp(SEXP x, SEXP y, SEXP src, SEXP h, SEXP w, SEXP kern) {
  clab::who = "kerexp";
  csp::pts a = csp::xy(x, y, R_NilValue, csp::planar);
  csp::pts s = csp::xy(src, R_NilValue, R_NilValue, csp::planar);
  const double hh = clab::one(h, "h", 0.0, HUGE_VAL);
  const double *wv = swts(w, s.n);
  static const char *const ko[] = {"gaussian", "quartic"};
  const int k = csp::pick(kern, "kern", ko, 2);
  SEXP o = PROTECT(real(a.n));
  kc(clab_kerexp(a.x, a.y, a.n, s.x, s.y, s.n, hh, wv, k, REAL(o)), "the kernel refused the input.");
  UNPROTECT(3 + (wv != nullptr));
  return o;
}


SEXP C_clab_r_lerprc(SEXP w, SEXP rho, SEXP tau) {
  clab::who = "lerprc";
  const int n = nareas(w);
  SEXP c = csp::wts(w, n);
  const double r = clab::p01(rho, "rho"), t = clab::one(tau, "tau", 0.0, HUGE_VAL);
  SEXP q = PROTECT(Rf_allocMatrix(REALSXP, n, n));
  kc(clab_lerprc(INTEGER(VECTOR_ELT(c, 0)), INTEGER(VECTOR_ELT(c, 1)), REAL(VECTOR_ELT(c, 2)), n, r, t, REAL(q)),
     "the kernel refused the input.");
  UNPROTECT(2);
  return q;
}

SEXP C_clab_r_bym2sc(SEXP w, SEXP tau, SEXP phi, SEXP tol) {
  clab::who = "bym2sc";
  const int n = nareas(w);
  if (n < 2) fail("`w` must link at least 2 areas.");
  SEXP c = csp::wts(w, n);
  const double t = clab::one(tau, "tau", 0.0, HUGE_VAL), ph = clab::p01(phi, "phi");
  const double tl = clab::one(tol, "tol", 0.0, 1.0);
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(r, 1, Rf_allocMatrix(REALSXP, n, n));
  SET_VECTOR_ELT(r, 2, Rf_allocMatrix(REALSXP, n, n));
  double s = 0;
  kc(clab_bym2(INTEGER(VECTOR_ELT(c, 0)), INTEGER(VECTOR_ELT(c, 1)), REAL(VECTOR_ELT(c, 2)), n, t, ph, tl, &s,
               REAL(VECTOR_ELT(r, 1)), REAL(VECTOR_ELT(r, 2))),
     "every area needs a neighbour.");
  SET_VECTOR_ELT(r, 0, Rf_ScalarReal(s));
  static const char *const nm[] = {"scale", "ginv", "cov"};
  lst(r, nm, nullptr);
  UNPROTECT(2);
  return r;
}

SEXP C_clab_r_ecolrg(SEXP y, SEXP e, SEXP x, SEXP family, SEXP z, SEXP tol, SEXP maxit) {
  clab::who = "ecolrg";
  SEXP yv = PROTECT(num(y, "y", -1, 0.0));
  const int n = clab::xlen(yv, "y");
  SEXP ev = PROTECT(num(e, "e", n, 0.0, true));
  int p = 0, q = 0;
  SEXP xm = zmat(x, "x", n, false, &p);
  varies(xm, n, p, "x");
  static const char *const fo[] = {"poisson", "negbin", "zip"};
  const int fam = csp::pick(family, "family", fo, 3);
  SEXP zm = zmat(z, "z", n, false, &q);
  varies(zm, n, q, "z");
  const double tl = clab::one(tol, "tol", 0.0, 1.0);
  const int mi = clab::whole(maxit, "maxit", 1);
  if (n <= p + 1) fail("`y` must have more values than the model has coefficients.");
  SEXP r = PROTECT(Rf_allocVector(VECSXP, fam == 0 ? 3 : 4));
  SEXP cf = PROTECT(real(p + 1)), fit = PROTECT(real(n));
  double *zc = dal(q + 1), th = 0, ll = 0;
  kc(clab_ecolrg(REAL(yv), REAL(ev), n, REAL(xm), p, fam, REAL(zm), q, tl, mi,
                 REAL(cf), zc, &th, REAL(fit), &ll),
     "the design is singular, or the fit diverged.");
  SET_VECTOR_ELT(r, 0, cf);
  if (fam == 0) {
    SET_VECTOR_ELT(r, 1, fit);
    SET_VECTOR_ELT(r, 2, Rf_ScalarReal(ll));
    static const char *const nm[] = {"coef", "fitted", "loglik"};
    lst(r, nm, "ecolrg");
  } else if (fam == 1) {
    SET_VECTOR_ELT(r, 1, Rf_ScalarReal(th));
    SET_VECTOR_ELT(r, 2, fit);
    SET_VECTOR_ELT(r, 3, Rf_ScalarReal(ll));
    static const char *const nm[] = {"coef", "theta", "fitted", "loglik"};
    lst(r, nm, "ecolrg");
  } else {
    SEXP zv = PROTECT(real(q + 1));
    for (int j = 0; j <= q; ++j) REAL(zv)[j] = zc[j];
    SET_VECTOR_ELT(r, 1, zv);
    SET_VECTOR_ELT(r, 2, fit);
    SET_VECTOR_ELT(r, 3, Rf_ScalarReal(ll));
    static const char *const nm[] = {"coef", "zcoef", "fitted", "loglik"};
    lst(r, nm, "ecolrg");
    UNPROTECT(1);
  }
  UNPROTECT(7);
  return r;
}

SEXP C_clab_r_spbias(SEXP d, SEXP b, SEXP ub) {
  clab::who = "spbias";
  clab::unlab(d, "d");
  if (!clab::isnum(d) || !Rf_isMatrix(d)) fail("`d` must be a numeric matrix of hop distances.");
  const int nr = Rf_nrows(d), nc = Rf_ncols(d);
  SEXP dv = PROTECT(clab::dbl(d));
  const R_xlen_t m = XLENGTH(dv);
  int *di = ial(m);
  for (R_xlen_t q = 0; q < m; ++q) {
    const double e = REAL(dv)[q];
    if (!std::isfinite(e) || e != std::floor(e) || e < -1 || e > 2147483647.0) {
      fail("`d` must hold whole numbers from 0, or -1 for an unreachable pair.");
    }
    di[q] = static_cast<int>(e);
  }
  SEXP bv = PROTECT(num(b, "b"));
  const int nb = clab::xlen(bv, "b");
  SEXP uv = PROTECT(num(ub, "ub", 1));
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(r, 0, Rf_allocMatrix(REALSXP, nr, nc));
  kc(clab_spbias(di, nr, nc, REAL(bv), nb, REAL(uv)[0], REAL(VECTOR_ELT(r, 0))), "the kernel refused the input.");
  SET_VECTOR_ELT(r, 1, Rf_ScalarReal(REAL(uv)[0]));
  static const char *const nm[] = {"bias", "ub"};
  lst(r, nm, nullptr);
  UNPROTECT(4);
  return r;
}

SEXP C_clab_r_spsir(SEXP i0, SEXP w, SEXP pop, SEXP beta, SEXP gamma, SEXP kappa, SEXP tmax, SEXP dt) {
  clab::who = "spsir";
  SEXP iv = PROTECT(num(i0, "i0", -1, 0.0));
  const int n = clab::xlen(iv, "i0");
  SEXP pv;
  if (Rf_isNull(pop)) {
    pv = PROTECT(real(n));
    for (int i = 0; i < n; ++i) REAL(pv)[i] = 1000;
  } else {
    pv = PROTECT(num(pop, "pop", n, 0.0, true));
  }
  for (int i = 0; i < n; ++i) {
    if (REAL(iv)[i] > REAL(pv)[i]) fail("`i0` must not exceed `pop`.");
  }
  SEXP c = Rf_isNull(w) ? PROTECT(R_NilValue) : csp::wts(w, n);
  SEXP bv = PROTECT(num(beta, "beta", 1, 0.0)), gv = PROTECT(num(gamma, "gamma", 1, 0.0));
  const double k = clab::p01(kappa, "kappa");
  const double tm = clab::one(tmax, "tmax", 0.0, HUGE_VAL), h = clab::one(dt, "dt", 0.0, HUGE_VAL);
  const double st = std::nearbyint(tm / h);
  if (st > 1e6) fail("`tmax` / `dt` must be at most 1e6 steps.");
  const int steps = static_cast<int>(st);
  SEXP S = PROTECT(Rf_allocMatrix(REALSXP, steps + 1, n)), I = PROTECT(Rf_allocMatrix(REALSXP, steps + 1, n)),
       R = PROTECT(Rf_allocMatrix(REALSXP, steps + 1, n));
  const bool hw = !Rf_isNull(w);
  kc(clab_spsir(REAL(iv), n, REAL(pv), hw ? INTEGER(VECTOR_ELT(c, 0)) : nullptr,
                hw ? INTEGER(VECTOR_ELT(c, 1)) : nullptr, hw ? REAL(VECTOR_ELT(c, 2)) : nullptr, REAL(bv)[0],
                REAL(gv)[0], k, h, steps, REAL(S), REAL(I), REAL(R)),
     "the kernel refused the input.");
  SEXP tv = PROTECT(real(steps + 1));
  for (int q = 0; q <= steps; ++q) REAL(tv)[q] = q * h;
  double rs = 0, ns = 0;
  for (int i = 0; i < n; ++i) {
    rs += REAL(R)[steps + static_cast<R_xlen_t>(steps + 1) * i];
    ns += REAL(pv)[i];
  }
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 6));
  SET_VECTOR_ELT(r, 0, Rf_ScalarReal(rs / ns));
  SET_VECTOR_ELT(r, 1, Rf_ScalarReal(REAL(bv)[0] / REAL(gv)[0]));
  SET_VECTOR_ELT(r, 2, tv);
  SET_VECTOR_ELT(r, 3, S);
  SET_VECTOR_ELT(r, 4, I);
  SET_VECTOR_ELT(r, 5, R);
  static const char *const nm[] = {"attack", "r0", "t", "S", "I", "R"};
  lst(r, nm, "spsir");
  UNPROTECT(10);
  return r;
}


SEXP C_clab_r_sfrail(SEXP time, SEXP event, SEXP x, SEXP reg, SEXP w, SEXP tau, SEXP tol) {
  clab::who = "sfrail";
  return surv(time, event, x, reg, w, tau, tol, 0);
}

SEXP C_clab_r_sfcure(SEXP time, SEXP event, SEXP x, SEXP reg, SEXP w, SEXP tau, SEXP tol) {
  clab::who = "sfcure";
  return surv(time, event, x, reg, w, tau, tol, 1);
}

SEXP C_clab_r_evfilt(SEXP v, SEXP w, SEXP x, SEXP method, SEXP cand, SEXP tol, SEXP maxv) {
  clab::who = "evfilt";
  SEXP yv = PROTECT(num(v, "v"));
  const int n = clab::xlen(yv, "v");
  if (n < 3) fail("`v` must have at least 3 values.");
  if (kutil::cnst(REAL(yv), n)) fail("`v` is constant: there is no autocorrelation to filter.");
  SEXP c = csp::wts(w, n);
  int p = 0;
  SEXP xm = zmat(x, "x", n, false, &p);
  varies(xm, n, p, "x");
  static const char *const mo[] = {"moran", "aic", "press", "r2"};
  static const char *const co[] = {"positive", "negative", "all"};
  const int cr = csp::pick(method, "method", mo, 4), ca = csp::pick(cand, "cand", co, 3);
  const double tl = clab::one(tol, "tol", 0.0, HUGE_VAL);
  const int mv = Rf_isNull(maxv) ? -1 : clab::whole(maxv, "maxv", 0);
  int *sel = ial(n), ns = 0;
  double *cf = dal(p + n), st[4];
  SEXP res = PROTECT(real(n)), mi = PROTECT(real(n - 1));
  kc(clab_evfilt(REAL(yv), n, REAL(xm), p, INTEGER(VECTOR_ELT(c, 0)), INTEGER(VECTOR_ELT(c, 1)),
                 REAL(VECTOR_ELT(c, 2)), cr, ca, tl, mv, sel, &ns, cf, st, REAL(res), REAL(mi)),
     "the weights have no links, or the design is singular.");
  SEXP sv = PROTECT(ints(ns)), cv = PROTECT(real(p + 1 + ns));
  for (int q = 0; q < ns; ++q) INTEGER(sv)[q] = sel[q];
  for (int q = 0; q < p + 1 + ns; ++q) REAL(cv)[q] = cf[q];
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 8));
  SET_VECTOR_ELT(r, 0, sv);
  SET_VECTOR_ELT(r, 1, cv);
  for (int j = 0; j < 4; ++j) SET_VECTOR_ELT(r, 2 + j, Rf_ScalarReal(st[j]));
  SET_VECTOR_ELT(r, 6, res);
  SET_VECTOR_ELT(r, 7, mi);
  static const char *const nm[] = {"sel", "coef", "imoran", "aic", "r2", "fr2", "resid", "mi"};
  lst(r, nm, "evfilt");
  UNPROTECT(8);
  return r;
}

}  // extern "C"
