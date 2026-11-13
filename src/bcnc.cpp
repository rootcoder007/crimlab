// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the concentration family (cconc, hilltl): each
// checks its arguments, calls the kernels and builds the R result.
// Nothing calls back into the R evaluator.
#include <algorithm>
#include <cmath>
#include <vector>

#include "bcol.h"
#include "rchk.h"

using clab::fail;

namespace {

// The Gini of each group of units: units sorted by group, the kernel run
// on each group's values and weights.
SEXP bygini(const double *x, const double *w, int n, csp::lab g) {
  bj::cols f(4, g.g);
  f.put("group", g.lev);
  double *un = f.dbl("units"), *tot = f.dbl("total"), *gi = f.dbl("gini");
  SEXP sc = PROTECT(Rf_allocVector(REALSXP, 4 * static_cast<R_xlen_t>(n) + 2));
  SEXP ic = PROTECT(Rf_allocVector(INTSXP, 2 * static_cast<R_xlen_t>(n) + g.g + 1));
  double *xs = REAL(sc), *wsub = xs + n, *lp = wsub + n, *lv = lp + n + 1;
  int *start = INTEGER(ic), *ord = start + g.g + 1, *o = ord + n;
  const int *code = INTEGER(g.code);
  // a counting sort of the units by group
  for (int k = 0; k <= g.g; ++k) start[k] = 0;
  for (int i = 0; i < n; ++i) ++start[code[i] + 1];
  for (int k = 0; k < g.g; ++k) start[k + 1] += start[k];
  for (int i = 0; i < n; ++i) ord[start[code[i]]++] = i;
  for (int k = g.g; k > 0; --k) start[k] = start[k - 1];
  start[0] = 0;
  double st[16];
  for (int k = 0; k < g.g; ++k) {
    const int a = start[k], m = start[k + 1] - a;
    for (int j = 0; j < m; ++j) {
      xs[j] = x[ord[a + j]];
      wsub[j] = w ? w[ord[a + j]] : 1.0;
    }
    if (clab_cconc(xs, wsub, m, nullptr, 0, st, lp, lv, nullptr, nullptr, o) == CLAB_OK) {
      un[k] = st[0];
      tot[k] = st[1];
      gi[k] = st[10];
    } else {
      // no events in the group: no concentration to measure
      double tw = 0;
      for (int j = 0; j < m; ++j) tw += wsub[j];
      un[k] = tw;
      tot[k] = 0;
      gi[k] = NA_REAL;
    }
  }
  SEXP out = f.done();
  UNPROTECT(3);
  return out;
}

}  // namespace

extern "C" {

SEXP C_clab_r_cconc(SEXP x, SEXP w, SEXP by, SEXP fracs, SEXP pop, SEXP nsim, SEXP seed) {
  clab::who = "cconc";
  int np = 0;
  const bool labs = !clab::isnum(x) || Rf_isFactor(x);
  SEXP xv, wv = R_NilValue, lev = R_NilValue, ev = R_NilValue;
  int n, ne = 0;
  if (labs) {
    // events: the place of each, counted per place
    ne = bj::len(x, "x");
    if (ne < 2) fail("`x` must hold at least two events when it gives their places.");
    if (!Rf_isNull(w)) fail("`w` goes with counts: with places of events, each event counts once.");
    if (!Rf_isNull(pop)) fail("`pop` goes with counts per unit, not with places of events.");
    csp::lab p = bj::labs(x, "x", ne);
    np += 2;
    ev = p.code;
    lev = p.lev;
    n = p.g;
    xv = PROTECT(Rf_allocVector(REALSXP, n));
    ++np;
    std::fill(REAL(xv), REAL(xv) + n, 0.0);
    for (int i = 0; i < ne; ++i) REAL(xv)[INTEGER(ev)[i]] += 1;
  } else {
    xv = PROTECT(bj::nums(x, "x", -1, 0.0));
    ++np;
    n = bj::len(xv, "x");
    if (!Rf_isNull(w)) {
      wv = PROTECT(bj::nums(w, "w", n, 0.0, true));
      ++np;
    }
  }
  SEXP fr = PROTECT(bj::nums(fracs, "fracs", -1, 0.0));
  ++np;
  const int nf = bj::len(fr, "fracs");
  for (int f = 0; f < nf; ++f) {
    if (REAL(fr)[f] > 1) fail("`fracs` must be shares from 0 to 1.");
  }
  const int ns = clab::whole(nsim, "nsim", 0);
  const double *wp = Rf_isNull(wv) ? nullptr : REAL(wv);
  SEXP st = PROTECT(Rf_allocVector(REALSXP, 16));
  ++np;
  bj::cols lc(2, n + 1);
  double *lp = lc.dbl("pop"), *lvv = lc.dbl("share");
  bj::cols tc(3, nf);
  double *fv = tc.dbl("frac"), *tn = tc.dbl("units"), *ts = tc.dbl("share");
  std::copy(REAL(fr), REAL(fr) + nf, fv);
  SEXP ow = PROTECT(Rf_allocVector(INTSXP, n));
  ++np;
  np += 2;
  if (clab_cconc(REAL(xv), wp, n, REAL(fr), nf, REAL(st), lp, lvv, ts, tn, INTEGER(ow)) !=
      CLAB_OK) {
    fail("`x` has no events: every count is 0.");
  }
  static const char *const sn[] = {"units", "total", "mean",  "var",  "disp", "ivar",
                                   "isd",   "zero",  "znull", "zgap", "gini", "gpos",
                                   "gnull", "vpop",  "cvb",   "nullb"};
  SEXP stat = PROTECT(clab::nvec(16, sn, REAL(st)));
  ++np;
  for (int k = 0; k < 16; ++k) REAL(stat)[k] = bj::na(REAL(stat)[k]);
  // the Poisson null at the observed mean, over the same number of units
  SEXP nul = R_NilValue;
  if (ns > 0) {
    const double un = REAL(st)[0];
    if (std::fabs(un - std::round(un)) > 1e-9 || un > 1e8) {
      fail("`nsim` needs whole-number weights (and at most 1e8 units).");
    }
    const unsigned int sd = csp::seed(seed);
    const int nu = static_cast<int>(std::round(un));
    SEXP gs = PROTECT(Rf_allocVector(REALSXP, ns));
    SEXP wk = PROTECT(Rf_allocVector(INTSXP, nu));
    np += 2;
    bj::kc(clab_cnull(REAL(st)[2], nu, ns, sd, REAL(gs), INTEGER(wk)),
           "the mean count is above 1e8: no null to simulate.");
    double *g = REAL(gs);
    const double gobs = REAL(st)[10];
    int ge = 0;
    double sum = 0;
    for (int s = 0; s < ns; ++s) {
      ge += g[s] >= gobs - 1e-12;
      sum += g[s];
    }
    std::sort(g, g + ns);
    // the 0.95 quantile, type 7
    const double h = 0.95 * (ns - 1);
    const int lo = static_cast<int>(std::floor(h));
    const double q = g[lo] + (h - lo) * (g[std::min(lo + 1, ns - 1)] - g[lo]);
    const double nv[5] = {static_cast<double>(ns), sum / ns, q, (1.0 + ge) / (ns + 1.0),
                          gobs - sum / ns};
    static const char *const nn[] = {"nsim", "mean", "q95", "p", "excess"};
    nul = PROTECT(clab::nvec(5, nn, nv));
    ++np;
  }
  SEXP byf = R_NilValue;
  if (!Rf_isNull(by)) {
    if (labs) {
      // groups of events: the places counted within each group
      csp::lab g = bj::labs(by, "by", ne);
      np += 2;
      std::vector<long long> key(static_cast<std::size_t>(ne));
      for (int i = 0; i < ne; ++i) {
        key[static_cast<std::size_t>(i)] =
            static_cast<long long>(INTEGER(g.code)[i]) * n + INTEGER(ev)[i];
      }
      std::sort(key.begin(), key.end());
      int m = 0;
      for (int i = 0; i < ne; ++i) m += i == 0 || key[i] != key[i - 1];
      SEXP cx = PROTECT(Rf_allocVector(REALSXP, m));
      SEXP cg = PROTECT(Rf_allocVector(INTSXP, m));
      np += 2;
      int j = -1;
      for (int i = 0; i < ne; ++i) {
        if (i == 0 || key[i] != key[i - 1]) {
          ++j;
          REAL(cx)[j] = 0;
          INTEGER(cg)[j] = static_cast<int>(key[i] / n);
        }
        REAL(cx)[j] += 1;
      }
      csp::lab gg = g;
      gg.code = cg;
      byf = PROTECT(bygini(REAL(cx), nullptr, m, gg));
      ++np;
    } else {
      csp::lab g = bj::labs(by, "by", n);
      np += 2;
      byf = PROTECT(bygini(REAL(xv), wp, n, g));
      ++np;
    }
  }
  SEXP scl = R_NilValue;
  if (!Rf_isNull(pop)) {
    SEXP pv = PROTECT(bj::nums(pop, "pop", n, 0.0, true));
    ++np;
    double sc[5];
    const int rc = clab_urbsc(REAL(xv), REAL(pv), n, sc);
    if (rc != CLAB_OK) {
      for (int k = 0; k < 4; ++k) sc[k] = NA_REAL;
    }
    static const char *const cn[] = {"beta", "pref", "se", "r2", "used"};
    scl = PROTECT(clab::nvec(5, cn, sc));
    ++np;
  }
  SEXP grw = R_NilValue;
  if (labs) {
    bj::cols gc(5, ne);
    int *ix = gc.itg("i");
    double *di = gc.dbl("dist"), *ex = gc.dbl("expd"), *glo = gc.dbl("lo"), *ghi = gc.dbl("hi");
    SEXP sw = PROTECT(Rf_allocVector(INTSXP, n));
    np += 2;
    double mh, sl;
    clab_cgrow(INTEGER(ev), ne, n, di, &mh, ex, glo, ghi, &sl, INTEGER(sw));
    for (int i = 0; i < ne; ++i) {
      ix[i] = i + 1;
      ex[i] = bj::na(ex[i]);
      glo[i] = bj::na(glo[i]);
      ghi[i] = bj::na(ghi[i]);
    }
    grw = PROTECT(Rf_allocVector(VECSXP, 3));
    ++np;
    SET_VECTOR_ELT(grw, 0, Rf_ScalarReal(mh));
    SET_VECTOR_ELT(grw, 1, Rf_ScalarReal(bj::na(sl)));
    SET_VECTOR_ELT(grw, 2, gc.done());
    static const char *const gn[] = {"mhat", "slope", "curve"};
    grw = bj::named(grw, gn);
    PROTECT(grw);
    ++np;
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 8));
  ++np;
  SET_VECTOR_ELT(out, 0, stat);
  SET_VECTOR_ELT(out, 1, tc.done());
  SET_VECTOR_ELT(out, 2, lc.done());
  SET_VECTOR_ELT(out, 3, nul);
  SET_VECTOR_ELT(out, 4, byf);
  SET_VECTOR_ELT(out, 5, scl);
  SET_VECTOR_ELT(out, 6, grw);
  SET_VECTOR_ELT(out, 7, labs ? lev : R_NilValue);
  static const char *const nm[] = {"stat", "top",   "lorenz", "null",
                                   "by",   "scale", "grow",   "places"};
  UNPROTECT(np);
  return bj::named(out, nm);
}

SEXP C_clab_r_hilltl(SEXP x, SEXP xmin, SEXP discrete, SEXP approx, SEXP mintl) {
  clab::who = "hilltl";
  int np = 0;
  SEXP xv = PROTECT(bj::nums(x, "x", -1, 0.0));
  ++np;
  const int n0 = bj::len(xv, "x");
  // the positive values; zeros have no place on a power law
  int n = 0;
  for (int i = 0; i < n0; ++i) n += REAL(xv)[i] > 0;
  if (n == 0) fail("`x` has no values above 0 to fit a tail to.");
  SEXP pv = PROTECT(Rf_allocVector(REALSXP, n));
  ++np;
  for (int i = 0, j = 0; i < n0; ++i) {
    if (REAL(xv)[i] > 0) REAL(pv)[j++] = REAL(xv)[i];
  }
  const bool dis = csp::flag(discrete, "discrete"), apx = csp::flag(approx, "approx");
  const int mt = std::max(2, clab::whole(mintl, "mintl", 1));
  double xm;
  if (Rf_isNull(xmin)) {
    // the 50th largest value, or the smallest when there are fewer
    std::vector<double> s(REAL(pv), REAL(pv) + n);
    std::sort(s.begin(), s.end(), [](double a, double b) { return a > b; });
    xm = n >= 50 ? s[49] : s[n - 1];
  } else {
    xm = clab::one(xmin, "xmin", 0.0, R_PosInf);
  }
  SEXP wk = PROTECT(Rf_allocVector(REALSXP, n));
  ++np;
  double o[5];
  const int rc = clab_hilltl(REAL(pv), n, xm, dis, apx, o, REAL(wk));
  const int ntl = static_cast<int>(o[3]);
  const char *meth;
  if (ntl < mt) {
    meth = "too few tail observations";
  } else if (rc != CLAB_OK) {
    meth = "not estimable: every tail value is the same";
  } else if (dis && !apx) {
    meth = o[4] > 0 ? "discrete maximum likelihood (alpha at an end of the search)"
                    : "discrete maximum likelihood";
  } else {
    meth = dis ? "continuity-corrected Hill estimator (approximate)" : "continuous Hill estimator";
  }
  const bool fit = ntl >= mt && rc == CLAB_OK;
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 7));
  ++np;
  SET_VECTOR_ELT(out, 0, Rf_ScalarReal(fit ? o[0] : NA_REAL));
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(fit ? bj::na(o[1]) : NA_REAL));
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(xm));
  SET_VECTOR_ELT(out, 3, Rf_ScalarInteger(ntl));
  SET_VECTOR_ELT(out, 4, Rf_ScalarReal(fit ? o[2] : NA_REAL));
  SET_VECTOR_ELT(out, 5, Rf_ScalarLogical(fit && ntl >= 50 && o[4] == 0));
  SET_VECTOR_ELT(out, 6, Rf_mkString(meth));
  static const char *const nm[] = {"alpha", "se", "xmin", "ntl", "ks", "reliable", "method"};
  UNPROTECT(np);
  return bj::named(out, nm);
}

SEXP C_clab_r_hzta(SEXP s, SEXP q) {
  clab::who = "hzta";
  SEXP sv = PROTECT(bj::nums(s, "s", -1, 1.0, true));
  const double qv = clab::one(q, "q", 0.0, R_PosInf);
  const R_xlen_t n = XLENGTH(sv);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  for (R_xlen_t i = 0; i < n; ++i) clab_hzta(REAL(sv)[i], qv, REAL(out) + i);
  UNPROTECT(2);
  return out;
}

}  // extern "C"
