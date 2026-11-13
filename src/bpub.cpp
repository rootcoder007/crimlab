// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the published-table family (bandx, pbbnd,
// yoypx): each checks its arguments, calls the kernels and builds the R
// result. The band labels and the text cells are read in R; exact
// binomial tails and beta quantiles come from R's maths library. Nothing
// calls back into the R evaluator.
#include <algorithm>
#include <cmath>
#include <vector>

#include "bcol.h"
#include "rchk.h"

extern "C" double Rf_qbeta(double, double, double, int, int);
extern "C" double Rf_dbinom(double, double, double, int);
extern "C" double Rf_pbinom(double, double, double, int, int);

using clab::fail;

namespace {

// One number, or NaN for NULL.
double opt(SEXP x, const char *nm, double lo) {
  return Rf_isNull(x) ? NAN : bj::inr(x, nm, lo, R_PosInf);
}

// The two-sided exact binomial p-value of x successes in n trials at
// chance p, as binom.test(): the tail beyond x and the far tail's points
// no more likely than x (relative 1e-7). The far tail's edge is found by
// bisection: the probabilities fall monotonically away from the mode.
double binp(double x, double n, double p) {
  const double d = Rf_dbinom(x, n, p, 0), m = n * p, re = 1 + 1e-7;
  if (x == m) return 1.0;
  if (x < m) {
    // the first k from ceil(m) with dbinom(k) <= d re
    double a = std::ceil(m), b = n + 1;
    while (a < b) {
      const double c = std::floor((a + b) / 2);
      if (Rf_dbinom(c, n, p, 0) <= d * re)
        b = c;
      else
        a = c + 1;
    }
    return std::fmin(1.0, Rf_pbinom(x, n, p, 1, 0) + (a <= n ? Rf_pbinom(a - 1, n, p, 0, 0) : 0));
  }
  // the last k up to floor(m) with dbinom(k) <= d re
  double a = -1, b = std::floor(m);
  while (a < b) {
    const double c = std::ceil((a + b) / 2);
    if (Rf_dbinom(c, n, p, 0) <= d * re)
      a = c;
    else
      b = c - 1;
  }
  return std::fmin(1.0, (a >= 0 ? Rf_pbinom(a, n, p, 1, 0) : 0) + Rf_pbinom(x - 1, n, p, 0, 0));
}

}  // namespace

extern "C" {

SEXP C_clab_r_bandx(SEXP lo, SEXP hi, SEXP olo, SEXP oup, SEXP rule, SEXP cap, SEXP fac, SEXP flr) {
  clab::who = "bandx";
  // lo, hi, olo, oup come from the R reader of the labels
  const int n = bj::len(lo, "x");
  static const char *const rules[] = {"midpoint", "lower", "upper", "geometric"};
  const int r = csp::pick(rule, "rule", rules, 4);
  const double cp = Rf_isNull(cap) ? NAN : bj::inr(cap, "cap", -INFINITY, R_PosInf);
  const double fc = clab::one(fac, "fac", 0.0, R_PosInf);
  const double fl = Rf_isNull(flr) ? NAN : bj::inr(flr, "floor", -INFINITY, R_PosInf);
  bj::cols f(3, n);
  double *v = f.dbl("value"), *bl = f.dbl("lo"), *bh = f.dbl("hi");
  if (clab_bandx(REAL(lo), REAL(hi), LOGICAL(olo), LOGICAL(oup), n, r, cp, fc, fl, v, bl, bh) !=
      CLAB_OK) {
    fail("the open band's assumed end (`cap` or `floor`) lies on the wrong side of its other end.");
  }
  for (int i = 0; i < n; ++i) v[i] = bj::na(v[i]);
  SEXP out = f.done();
  UNPROTECT(1);
  return out;
}

SEXP C_clab_r_pbbnd(SEXP val, SEXP st, SEXP lim, SEXP rnd, SEXP kind, SEXP sup, SEXP pval, SEXP pst,
                    SEXP plim, SEXP pop, SEXP ppop, SEXP per) {
  clab::who = "pbbnd";
  int np = 0;
  const int n = bj::len(val, "x");
  const double rn = opt(rnd, "rnd", 0.0), sp = opt(sup, "suplim", 0.0);
  if (rn == 0) fail("`rnd` must be above 0.");
  static const char *const kinds[] = {"nearest", "random"};
  const int kd = csp::pick(kind, "kind", kinds, 2);
  for (int i = 0; i < n; ++i) {
    if (REAL(val)[i] < 0) fail("`x` must hold counts from 0.");
  }
  bj::cols f(4, n);
  ++np;
  f.put("value", val);
  double *lo = f.dbl("lower"), *hi = f.dbl("upper");
  SEXP so = PROTECT(Rf_allocVector(INTSXP, n));
  ++np;
  if (clab_pbbnd(REAL(val), INTEGER(st), REAL(lim), n, rn, kd, sp, lo, hi, INTEGER(so)) !=
      CLAB_OK) {
    fail("a suppressed cell has no limit: give `suplim` or write the cell as \"<k\".");
  }
  static const char *const sn[] = {"exact", "rounded", "suppressed", "missing"};
  SEXP ss = PROTECT(Rf_allocVector(STRSXP, n));
  ++np;
  for (int i = 0; i < n; ++i) {
    SET_STRING_ELT(ss, i, Rf_mkChar(sn[INTEGER(so)[i]]));
    lo[i] = bj::na(lo[i]);
    hi[i] = bj::na(hi[i]);
  }
  f.put("status", ss);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  ++np;
  SET_VECTOR_ELT(out, 0, f.done());
  if (!Rf_isNull(pval)) {
    if (XLENGTH(pval) != n) fail("`prev` must have one cell per cell of `x`, %d.", n);
    for (int i = 0; i < n; ++i) {
      if (REAL(pval)[i] < 0) fail("`prev` must hold counts from 0.");
    }
    SEXP w = PROTECT(Rf_allocVector(REALSXP, 10 * static_cast<R_xlen_t>(n)));
    SEXP so2 = PROTECT(Rf_allocVector(INTSXP, n));
    np += 2;
    double *plo = REAL(w), *phi = plo + n, *o = phi + n;
    if (clab_pbbnd(REAL(pval), INTEGER(pst), REAL(plim), n, rn, kd, sp, plo, phi, INTEGER(so2)) !=
        CLAB_OK) {
      fail("a suppressed cell of `prev` has no limit: give `suplim` or write it as \"<k\".");
    }
    SEXP pv = R_NilValue, qv = R_NilValue;
    if (!Rf_isNull(pop)) {
      pv = PROTECT(bj::nums(pop, "pop", -1, 0.0, true));
      ++np;
      if (XLENGTH(pv) != 1 && XLENGTH(pv) != n)
        fail("`pop` must be one number or one per cell, %d.", n);
    }
    if (!Rf_isNull(ppop)) {
      qv = PROTECT(bj::nums(ppop, "prepop", -1, 0.0, true));
      ++np;
      if (XLENGTH(qv) != 1 && XLENGTH(qv) != n) {
        fail("`prepop` must be one number or one per cell, %d.", n);
      }
    }
    // populations recycled to one per cell
    SEXP pr = PROTECT(Rf_allocVector(REALSXP, 2 * static_cast<R_xlen_t>(n)));
    ++np;
    for (int i = 0; i < n; ++i) {
      if (!Rf_isNull(pv)) REAL(pr)[i] = REAL(pv)[XLENGTH(pv) == 1 ? 0 : i];
      if (!Rf_isNull(qv)) REAL(pr)[n + i] = REAL(qv)[XLENGTH(qv) == 1 ? 0 : i];
    }
    const double pe = clab::one(per, "per", 0.0, R_PosInf);
    clab_pbchg(lo, hi, plo, phi, n, Rf_isNull(pv) ? nullptr : REAL(pr),
               Rf_isNull(qv) ? nullptr : REAL(pr) + n, pe, o);
    static const char *const cn[] = {"chglo",  "chgup",  "pctlo",     "pctup",
                                     "ratelo", "rateup", "preratelo", "prerateup"};
    const int nc = Rf_isNull(pv) ? 4 : Rf_isNull(qv) ? 6 : 8;
    bj::cols g(8, n);
    ++np;
    for (int c = 0; c < nc; ++c) {
      double *d = g.dbl(cn[c]);
      for (int i = 0; i < n; ++i) d[i] = bj::na(o[c * n + i]);
    }
    SET_VECTOR_ELT(out, 1, g.done());
  }
  static const char *const nm[] = {"bnd", "chg"};
  UNPROTECT(np);
  return bj::named(out, nm);
}

SEXP C_clab_r_yoypx(SEXP x, SEXP prd, SEXP by, SEXP lag, SEXP units, SEXP minbase, SEXP level,
                    SEXP rnd, SEXP kind, SEXP adj, SEXP alpha) {
  clab::who = "yoypx";
  int np = 0;
  static const char *const un[] = {"count", "continuous", "percent"};
  const int u = csp::pick(units, "units", un, 3);
  SEXP xv = PROTECT(bj::nums(x, "x", -1, u == 0 ? 0.0 : clab::ninf));
  ++np;
  const int n = bj::len(xv, "x");
  SEXP pv = PROTECT(bj::nums(prd, "prd", n));
  ++np;
  const double lg = clab::one(lag, "lag", 0.0, R_PosInf);
  const double mb = bj::inr(minbase, "minbase", 0.0, R_PosInf);
  const double lv = clab::one(level, "level", 0.0, 1.0);
  const double rn = opt(rnd, "rnd", 0.0);
  if (rn == 0) fail("`rnd` must be above 0.");
  if (u != 0 && !std::isnan(rn)) fail("`rnd` goes with counts (`units = \"count\"`).");
  static const char *const kinds[] = {"nearest", "random"};
  const int kd = csp::pick(kind, "kind", kinds, 2);
  static const char *const am[] = {"none", "bonferroni", "holm", "hochberg", "BH", "BY"};
  const int ad = csp::pick(adj, "adj", am, 6);
  const double al = clab::one(alpha, "alpha", 0.0, 1.0);
  // the groups, then the cells (group, period): rows of a cell summed
  // (counts) or averaged
  SEXP gc = R_NilValue, glev = R_NilValue;
  if (!Rf_isNull(by)) {
    csp::lab g = bj::labs(by, "by", n);
    np += 2;
    gc = g.code;
    glev = g.lev;
  }
  std::vector<int> o(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) o[static_cast<std::size_t>(i)] = i;
  const int *gi = Rf_isNull(gc) ? nullptr : INTEGER(gc);
  const double *pp = REAL(pv);
  std::stable_sort(o.begin(), o.end(), [&](int a, int b) {
    const int ga = gi ? gi[a] : 0, gb = gi ? gi[b] : 0;
    return ga < gb || (ga == gb && pp[a] < pp[b]);
  });
  int m = 0;
  for (int i = 0; i < n; ++i) {
    const int a = o[static_cast<std::size_t>(i)];
    const int b = i ? o[static_cast<std::size_t>(i - 1)] : -1;
    m += i == 0 || (gi ? gi[a] != gi[b] : false) || pp[a] != pp[b];
  }
  SEXP cv = PROTECT(Rf_allocVector(REALSXP, 2 * static_cast<R_xlen_t>(m)));
  SEXP cg = PROTECT(Rf_allocVector(INTSXP, 2 * static_cast<R_xlen_t>(m)));
  np += 2;
  double *cx = REAL(cv), *cp = cx + m;
  int *cgr = INTEGER(cg), *cn = cgr + m;
  for (int i = 0, j = -1; i < n; ++i) {
    const int a = o[static_cast<std::size_t>(i)];
    const int b = i ? o[static_cast<std::size_t>(i - 1)] : -1;
    if (i == 0 || (gi ? gi[a] != gi[b] : false) || pp[a] != pp[b]) {
      ++j;
      cx[j] = 0;
      cp[j] = pp[a];
      cgr[j] = gi ? gi[a] : 0;
      cn[j] = 0;
    }
    cx[j] += REAL(xv)[a];
    ++cn[j];
  }
  o.clear();
  o.shrink_to_fit();
  if (u != 0) {
    for (int j = 0; j < m; ++j) cx[j] /= cn[j];
  }
  bj::cols f(22, m);
  ++np;
  if (!Rf_isNull(glev)) {
    SEXP gl = PROTECT(Rf_allocVector(STRSXP, m));
    ++np;
    for (int j = 0; j < m; ++j) SET_STRING_ELT(gl, j, STRING_ELT(glev, cgr[j]));
    f.put("group", gl);
  }
  double *fp = f.dbl("prd"), *val = f.dbl("value"), *prev = f.dbl("prev"), *chg = f.dbl("change");
  double *pct = f.dbl(u == 2 ? "ppchg" : "pctchg");
  std::copy(cp, cp + m, fp);
  std::copy(cx, cx + m, val);
  clab_yoypx(cx, cp, cgr, m, lg, prev, chg, pct);
  for (int j = 0; j < m; ++j) {
    // a change in percentage points has no percent change
    if (u == 2) pct[j] = chg[j];
  }
  SEXP fl = PROTECT(Rf_allocVector(STRSXP, m));
  ++np;
  for (int j = 0; j < m; ++j) {
    const char *s = nullptr;
    if (std::isnan(prev[j]))
      s = "no earlier period";
    else if (u != 2 && prev[j] == 0)
      s = "earlier period is 0";
    else if (u != 2 && prev[j] < mb)
      s = "small base";
    SET_STRING_ELT(fl, j, s ? Rf_mkChar(s) : NA_STRING);
  }
  for (int j = 0; j < m; ++j) {
    prev[j] = bj::na(prev[j]);
    chg[j] = bj::na(chg[j]);
    pct[j] = bj::na(pct[j]);
  }
  f.put("flag", fl);
  if (u == 0) {
    // the percent change's exact interval, from the conditional binomial
    // of now among now and before, and the exact p-value of no change
    double *clo = f.dbl("pctlo"), *chi = f.dbl("pctup"), *pv2 = f.dbl("p"), *pa = f.dbl("padj");
    int *sg = f.lgl("signif");
    const double a = 1 - lv;
    SEXP po = PROTECT(Rf_allocVector(REALSXP, 2 * static_cast<R_xlen_t>(m)));
    SEXP oi = PROTECT(Rf_allocVector(INTSXP, m));
    np += 2;
    double *pk = REAL(po), *qk = pk + m;
    int k = 0;
    for (int j = 0; j < m; ++j) {
      const double v = val[j], b = prev[j];
      clo[j] = chi[j] = pv2[j] = pa[j] = NA_REAL;
      sg[j] = NA_LOGICAL;
      if (ISNAN(b) || v + b == 0) continue;
      const double l = v == 0 ? 0 : Rf_qbeta(a / 2, v, b + 1, 1, 0);
      const double h = b == 0 ? 1 : Rf_qbeta(1 - a / 2, v + 1, b, 1, 0);
      clo[j] = 100 * (l / (1 - l) - 1);
      chi[j] = h >= 1 ? R_PosInf : 100 * (h / (1 - h) - 1);
      pv2[j] = binp(std::round(v), std::round(v + b), 0.5);
      pk[k++] = pv2[j];
    }
    clab_pbadj(pk, k, ad, qk, INTEGER(oi));
    for (int j = 0, i = 0; j < m; ++j) {
      if (ISNAN(pv2[j])) continue;
      pa[j] = qk[i++];
      sg[j] = pa[j] < al;
    }
    if (!std::isnan(rn)) {
      // what the rounding of both values allows
      SEXP bw = PROTECT(Rf_allocVector(REALSXP, 12 * static_cast<R_xlen_t>(m)));
      SEXP bs = PROTECT(Rf_allocVector(INTSXP, 3 * static_cast<R_xlen_t>(m)));
      np += 2;
      double *nlo = REAL(bw), *nhi = nlo + m, *plo = nhi + m, *phi = plo + m, *eo = phi + m;
      int *zst = INTEGER(bs), *s1 = zst + m, *s2 = s1 + m;
      std::fill(zst, zst + m, 0);
      std::vector<double> lim(static_cast<std::size_t>(m), NAN), pz(static_cast<std::size_t>(m));
      for (int j = 0; j < m; ++j)
        pz[static_cast<std::size_t>(j)] = std::isnan(prev[j]) ? NAN : prev[j];
      clab_pbbnd(val, zst, lim.data(), m, rn, kd, NAN, nlo, nhi, s1);
      clab_pbbnd(pz.data(), zst, lim.data(), m, rn, kd, NAN, plo, phi, s2);
      clab_pbchg(nlo, nhi, plo, phi, m, nullptr, nullptr, 1, eo);
      double *el = f.dbl("envlo"), *eh = f.dbl("envup"), *ql = f.dbl("envplo"),
             *qh = f.dbl("envpup");
      // the percent change over both the sampling interval and the envelope
      double *cl = f.dbl("cpctlo"), *ch = f.dbl("cpctup");
      for (int j = 0; j < m; ++j) {
        el[j] = bj::na(eo[j]);
        eh[j] = bj::na(eo[m + j]);
        ql[j] = bj::na(eo[2 * m + j]);
        qh[j] = bj::na(eo[3 * m + j]);
        cl[j] = ISNAN(clo[j]) || ISNAN(ql[j]) ? NA_REAL : std::fmin(clo[j], ql[j]);
        ch[j] = ISNAN(chi[j]) || ISNAN(qh[j]) ? NA_REAL : std::fmax(chi[j], qh[j]);
      }
    }
  }
  SEXP out = PROTECT(f.done());
  ++np;
  UNPROTECT(np);
  return out;
}

}  // extern "C"
