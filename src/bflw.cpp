// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the custody and court-flow family (stkflw,
// lostay, bklog, dktflw, lftab, trnmat): each checks its arguments, calls
// the kernels and builds the R result. Quantiles of the t, beta and
// chi-square distributions come from R's maths library. Nothing calls
// back into the R evaluator.
#include <algorithm>
#include <cmath>
#include <vector>

#include "bcol.h"
#include "rchk.h"

extern "C" double Rf_qt(double, double, int, int);
extern "C" double Rf_qbeta(double, double, double, int, int);
extern "C" double Rf_qchisq(double, double, int, int);
extern "C" double Rf_pchisq(double, double, int, int);

using clab::fail;

namespace {

// Period labels: as.character(prd), or 1, 2, ... (protected by the
// caller).
SEXP plab(SEXP prd, int n) {
  if (Rf_isNull(prd)) {
    SEXP s = PROTECT(Rf_allocVector(STRSXP, n));
    for (int i = 0; i < n; ++i) SET_STRING_ELT(s, i, Rf_mkChar(std::to_string(i + 1).c_str()));
    UNPROTECT(1);
    return s;
  }
  clab::unlab(prd, "prd");
  clab::vec(prd, "prd");
  if (XLENGTH(prd) != n) fail("`prd` must have one label per period, %d.", n);
  SEXP s = clab::labs(prd);
  if (Rf_isNull(s)) fail("`prd` must be labels: numbers, text or a factor.");
  return s;
}

// A count per period, or NULL: numbers from 0, n of them.
SEXP per0(SEXP x, const char *nm, int n) {
  return Rf_isNull(x) ? R_NilValue : bj::nums(x, nm, n, 0.0);
}

// The exact interval of a ratio a / b of counts, from the conditional
// binomial of a among a + b, at tail p each side; NaN for no cases.
void rci(double a, double b, double p, double *lo, double *hi) {
  if (std::isnan(a) || std::isnan(b) || a + b == 0) {
    *lo = *hi = NA_REAL;
    return;
  }
  const double l = a == 0 ? 0 : Rf_qbeta(p, a, b + 1, 1, 0);
  const double h = b == 0 ? 1 : Rf_qbeta(1 - p, a + 1, b, 1, 0);
  *lo = l / (1 - l);
  *hi = h >= 1 ? R_PosInf : h / (1 - h);
}

}  // namespace

extern "C" {

SEXP C_clab_r_stkflw(SEXP people, SEXP days, SEXP adp, SEXP alos, SEXP t, SEXP prd, SEXP expo,
                     SEXP per, SEXP base, SEXP rel, SEXP early) {
  clab::who = "stkflw";
  int np = 0;
  const int given = !Rf_isNull(people) + !Rf_isNull(days) + !Rf_isNull(adp) + !Rf_isNull(alos);
  if (given != 2 || (!Rf_isNull(days) && !Rf_isNull(adp))) {
    fail("give two of `people`, `days`, `adp` and `alos` (not `days` with `adp`).");
  }
  // headcounts per period (a list of vectors) give the average daily
  // population as their means
  SEXP av = R_NilValue;
  if (!Rf_isNull(adp) && TYPEOF(adp) == VECSXP) {
    const int k = bj::len(adp, "adp");
    av = PROTECT(Rf_allocVector(REALSXP, k));
    ++np;
    for (int j = 0; j < k; ++j) {
      SEXP c = PROTECT(bj::nums(VECTOR_ELT(adp, j), "adp", -1, 0.0));
      REAL(av)[j] = clab::nmean(REAL(c), XLENGTH(c));
      UNPROTECT(1);
    }
  } else if (!Rf_isNull(adp)) {
    av = PROTECT(bj::nums(adp, "adp", -1, 0.0));
    ++np;
  }
  SEXP first = !Rf_isNull(people) ? people : !Rf_isNull(days) ? days : av;
  const int n = bj::len(first, "people");
  SEXP pv = PROTECT(Rf_isNull(people) ? R_NilValue : bj::nums(people, "people", n, 0.0, true));
  SEXP dv = PROTECT(per0(days, "days", n));
  SEXP lv = PROTECT(Rf_isNull(alos) ? R_NilValue : bj::nums(alos, "alos", n, 0.0, true));
  np += 3;
  if (!Rf_isNull(av) && XLENGTH(av) != n) fail("`adp` must have one value per period, %d.", n);
  SEXP t0 = PROTECT(bj::nums(t, "t", -1, 0.0, true));
  ++np;
  if (XLENGTH(t0) != 1 && XLENGTH(t0) != n)
    fail("`t` must be one number or one per period, %d.", n);
  SEXP tv = PROTECT(Rf_allocVector(REALSXP, n));
  ++np;
  for (int i = 0; i < n; ++i) REAL(tv)[i] = REAL(t0)[XLENGTH(t0) == 1 ? 0 : i];
  SEXP xv = R_NilValue, rv = R_NilValue;
  if (!Rf_isNull(expo)) {
    SEXP e0 = PROTECT(bj::nums(expo, "expo", -1, 0.0, true));
    ++np;
    if (XLENGTH(e0) != 1 && XLENGTH(e0) != n)
      fail("`expo` must be one number or one per period, %d.", n);
    xv = PROTECT(Rf_allocVector(REALSXP, n));
    ++np;
    for (int i = 0; i < n; ++i) REAL(xv)[i] = REAL(e0)[XLENGTH(e0) == 1 ? 0 : i];
  }
  const double pe = clab::one(per, "per", 0.0, R_PosInf);
  static const char *const bases[] = {"first", "previous"};
  const int b = csp::pick(base, "base", bases, 2);
  if (Rf_isNull(rel) != Rf_isNull(early)) fail("`rel` and `early` go together.");
  double ea = NAN;
  if (!Rf_isNull(rel)) {
    rv = PROTECT(bj::nums(rel, "rel", n, 0.0));
    ++np;
    ea = bj::inr(early, "early", 0.0, R_PosInf);
  }
  SEXP pl = PROTECT(plab(prd, n));
  ++np;
  SEXP ow = PROTECT(Rf_allocVector(REALSXP, 16 * static_cast<R_xlen_t>(n) + 4));
  ++np;
  const double *o = REAL(ow);
  clab_stkflw(Rf_isNull(pv) ? nullptr : REAL(pv), Rf_isNull(dv) ? nullptr : REAL(dv),
              Rf_isNull(av) ? nullptr : REAL(av), Rf_isNull(lv) ? nullptr : REAL(lv), REAL(tv), n,
              b, Rf_isNull(xv) ? nullptr : REAL(xv), pe, Rf_isNull(rv) ? nullptr : REAL(rv), ea,
              REAL(ow));
  bj::cols f(20, n);
  ++np;
  f.put("period", pl);
  f.put("t", tv);
  static const char *const cn[] = {"people",  "days",    "adp",    "alos",   "peochg",
                                   "aloschg", "dayschg", "adpchg", "pshare", "flowrate",
                                   "stkrate", "frchg",   "srchg",  "cut",    "cutpct"};
  for (int c = 0; c < 15; ++c) {
    if (c >= 9 && c <= 12 && Rf_isNull(xv)) continue;
    if (c >= 13 && Rf_isNull(rv)) continue;
    double *d = f.dbl(cn[c]);
    for (int i = 0; i < n; ++i) d[i] = bj::na(o[c * n + i]);
  }
  static const char *const mn[] = {"people", "days", "adp", "alos"};
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  ++np;
  SET_VECTOR_ELT(out, 0, f.done());
  SET_VECTOR_ELT(out, 1, clab::nvec(4, mn, o + 15 * n));
  SET_VECTOR_ELT(out, 2, Rf_mkString(bases[b]));
  SET_VECTOR_ELT(out, 3, Rf_ScalarReal(Rf_isNull(xv) ? NA_REAL : pe));
  static const char *const nm[] = {"flow", "mean", "base", "per"};
  UNPROTECT(np);
  return bj::named(out, nm);
}

SEXP C_clab_r_lostay(SEXP days, SEXP w, SEXP breaks, SEXP thr, SEXP flag, SEXP level) {
  clab::who = "lostay";
  int np = 0;
  SEXP xv = PROTECT(bj::nums(days, "days", -1, 0.0));
  ++np;
  const int n = bj::len(xv, "days");
  SEXP wv = R_NilValue;
  if (!Rf_isNull(w)) {
    wv = PROTECT(bj::nums(w, "w", n, 0.0));
    ++np;
    double tot = 0;
    for (int i = 0; i < n; ++i) {
      if (REAL(wv)[i] != std::floor(REAL(wv)[i])) fail("`w` must be whole-number counts.");
      tot += REAL(wv)[i];
    }
    if (!(tot > 0)) fail("`w` must have a positive total.");
  }
  SEXP bv = R_NilValue;
  int nb = 0;
  if (!Rf_isNull(breaks)) {
    clab::unlab(breaks, "breaks");
    if (!clab::isnum(breaks) || XLENGTH(breaks) < 2)
      fail("`breaks` must be at least two band edges.");
    bv = PROTECT(clab::dbl(breaks));
    ++np;
    nb = bj::len(bv, "breaks");
    for (int k = 0; k < nb; ++k) {
      if (std::isnan(REAL(bv)[k]) || (k > 0 && !(REAL(bv)[k] > REAL(bv)[k - 1]))) {
        fail("`breaks` must be increasing band edges with no missing values.");
      }
    }
  }
  SEXP tv = R_NilValue;
  int nt = 0;
  if (!Rf_isNull(thr)) {
    tv = PROTECT(bj::nums(thr, "thr", -1, 0.0));
    ++np;
    nt = bj::len(tv, "thr");
  }
  SEXP fv = R_NilValue;
  if (!Rf_isNull(flag)) {
    fv = PROTECT(Rf_allocVector(LGLSXP, n));
    ++np;
    clab::lgl(flag, "flag", n);
    std::copy(LOGICAL(flag), LOGICAL(flag) + n, LOGICAL(fv));
  }
  const double lv = clab::one(level, "level", 0.0, 1.0);
  double st[9], fit[6];
  SEXP bw = PROTECT(Rf_allocVector(REALSXP, 2 * static_cast<R_xlen_t>(nb) + nt + 2));
  SEXP ow = PROTECT(Rf_allocVector(INTSXP, n));
  np += 2;
  double *bt = REAL(bw), *tt = bt + 2 * nb;
  clab_lostay(REAL(xv), Rf_isNull(wv) ? nullptr : REAL(wv), n, Rf_isNull(bv) ? nullptr : REAL(bv),
              nb, Rf_isNull(tv) ? nullptr : REAL(tv), nt, Rf_isNull(fv) ? nullptr : LOGICAL(fv), st,
              bt, tt, fit, INTEGER(ow));
  const double nn = st[0];
  const double q = nn > 1 ? Rf_qt(0.5 + lv / 2, nn - 1, 1, 0) : NA_REAL;
  const double sv[11] = {nn,
                         st[1],
                         st[2],
                         bj::na(st[3]),
                         st[4],
                         st[5],
                         st[6],
                         st[7],
                         bj::na(st[8]),
                         nn > 1 ? st[2] - q * st[8] : NA_REAL,
                         nn > 1 ? st[2] + q * st[8] : NA_REAL};
  static const char *const sn[] = {"n",  "days", "mean", "sd", "median", "q1",
                                   "q3", "max",  "se",   "lo", "hi"};
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
  ++np;
  SET_VECTOR_ELT(out, 0, clab::nvec(11, sn, sv));
  if (nb >= 2) {
    bj::cols g(5, nb - 1);
    ++np;
    double *gl = g.dbl("lo"), *gh = g.dbl("hi"), *gn = g.dbl("n"), *gs = g.dbl("share");
    double *gf = Rf_isNull(fv) ? nullptr : g.dbl("flagged");
    for (int k = 0; k + 1 < nb; ++k) {
      gl[k] = REAL(bv)[k];
      gh[k] = REAL(bv)[k + 1];
      gn[k] = bt[k];
      gs[k] = bt[k] / nn;
      if (gf) gf[k] = bt[k] > 0 ? bt[nb - 1 + k] / bt[k] : NA_REAL;
    }
    SET_VECTOR_ELT(out, 1, g.done());
  }
  if (nt > 0) {
    bj::cols h(3, nt);
    ++np;
    h.put("thr", tv);
    double *at = h.dbl("atmost"), *ab = h.dbl("above");
    for (int i = 0; i < nt; ++i) {
      at[i] = tt[i];
      ab[i] = 1 - tt[i];
    }
    SET_VECTOR_ELT(out, 2, h.done());
  }
  bj::cols m(3, 3);
  ++np;
  SEXP ml = PROTECT(Rf_allocVector(STRSXP, 3));
  ++np;
  static const char *const mo[] = {"lognormal", "exponential", "pareto"};
  for (int k = 0; k < 3; ++k) SET_STRING_ELT(ml, k, Rf_mkChar(mo[k]));
  m.put("model", ml);
  double *ma = m.dbl("aic");
  int *mb = m.lgl("best");
  int best = -1;
  for (int k = 0; k < 3; ++k) {
    ma[k] = bj::na(fit[k]);
    if (!std::isnan(fit[k]) && (best < 0 || fit[k] < fit[best])) best = k;
  }
  for (int k = 0; k < 3; ++k) mb[k] = k == best;
  SET_VECTOR_ELT(out, 3, m.done());
  static const char *const pn[] = {"meanlog", "sdlog", "alpha"};
  const double pv[3] = {bj::na(fit[3]), bj::na(fit[4]), bj::na(fit[5])};
  SET_VECTOR_ELT(out, 4, clab::nvec(3, pn, pv));
  static const char *const nm[] = {"stay", "bands", "thr", "fit", "par"};
  UNPROTECT(np);
  return bj::named(out, nm);
}

SEXP C_clab_r_bklog(SEXP arrive, SEXP dispose, SEXP horizon, SEXP target) {
  clab::who = "bklog";
  int np = 0;
  SEXP av = PROTECT(bj::nums(arrive, "arrive", -1, 0.0));
  ++np;
  const int n = bj::len(av, "arrive");
  clab::unlab(dispose, "dispose");
  clab::vec(dispose, "dispose");
  if (!clab::isnum(dispose) || XLENGTH(dispose) != n) {
    fail("`dispose` must be numbers (NA for a pending case), one per case, %d.", n);
  }
  SEXP dv = PROTECT(clab::dcopy(PROTECT(clab::dbl(dispose))));
  np += 2;
  double h = 0;
  for (int i = 0; i < n; ++i) {
    const double d = REAL(dv)[i];
    if (std::isinf(d)) fail("`dispose` must be finite or NA.");
    if (!std::isnan(d) && d < REAL(av)[i])
      fail("a disposition precedes its arrival (case %d).", i + 1);
    h = std::fmax(h, std::fmax(REAL(av)[i], std::isnan(d) ? 0 : d));
  }
  if (!Rf_isNull(horizon)) h = clab::one(horizon, "horizon", 0.0, R_PosInf);
  double o[7];
  if (clab_bklog(REAL(av), REAL(dv), n, h, o) != CLAB_OK) {
    fail("every arrival must lie inside the horizon, after 0.");
  }
  for (int k = 0; k < 7; ++k) o[k] = bj::na(o[k]);
  static const char *const on[] = {"cases", "pending", "rate", "wait", "backlog", "occup", "resid"};
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  ++np;
  SET_VECTOR_ELT(out, 0, clab::nvec(7, on, o));
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(h));
  // the mean time a backlog target needs at this filing rate: W = L / lambda
  if (!Rf_isNull(target)) {
    const double tg = bj::inr(target, "target", 0.0, R_PosInf);
    SET_VECTOR_ELT(out, 2, Rf_ScalarReal(tg / o[2]));
  }
  static const char *const nm[] = {"est", "horizon", "need"};
  UNPROTECT(np);
  return bj::named(out, nm);
}

SEXP C_clab_r_dktflw(SEXP recv, SEXP disp, SEXP pend, SEXP prd, SEXP win, SEXP level, SEXP span,
                     SEXP mtime, SEXP shares, SEXP floors, SEXP mix, SEXP scores) {
  clab::who = "dktflw";
  int np = 0;
  if (Rf_isNull(recv) && Rf_isNull(mix)) fail("give `recv` and `disp`, or `mix`.");
  const int w = clab::whole(win, "win", 1);
  const double lv = clab::one(level, "level", 0.0, 1.0);
  const double a = 1 - lv;
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  ++np;
  if (!Rf_isNull(recv)) {
    SEXP rv = PROTECT(bj::nums(recv, "recv", -1, 0.0));
    ++np;
    const int n = bj::len(rv, "recv");
    if (Rf_isNull(disp)) fail("`disp` must come with `recv`.");
    SEXP dv = PROTECT(bj::nums(disp, "disp", n, 0.0));
    SEXP pv = PROTECT(per0(pend, "pend", n));
    SEXP pl = PROTECT(plab(prd, n));
    np += 3;
    SEXP ow = PROTECT(Rf_allocVector(REALSXP, 8 * static_cast<R_xlen_t>(n)));
    ++np;
    const double *o = REAL(ow);
    clab_dktflw(REAL(rv), REAL(dv), Rf_isNull(pv) ? nullptr : REAL(pv), n, w, REAL(ow));
    bj::cols f(20, n);
    ++np;
    f.put("period", pl);
    f.put("recv", rv);
    f.put("disp", dv);
    double *cl = f.dbl("clear"), *cl0 = f.dbl("lo"), *ch0 = f.dbl("hi"), *nt = f.dbl("net");
    double *rc = f.dbl("rclear"), *rl = f.dbl("rlo"), *rh = f.dbl("rhi");
    for (int i = 0; i < n; ++i) {
      cl[i] = bj::na(o[i]);
      rci(REAL(dv)[i], REAL(rv)[i], a / 2, cl0 + i, ch0 + i);
      nt[i] = o[n + i];
      rc[i] = bj::na(o[4 * n + i]);
      rci(o[3 * n + i], o[2 * n + i], a / 2, rl + i, rh + i);
    }
    if (!Rf_isNull(pv)) {
      f.put("pend", pv);
      double *pc = f.dbl("pchg"), *tc = f.dbl("toclear");
      for (int i = 0; i < n; ++i) {
        pc[i] = bj::na(o[5 * n + i]);
        tc[i] = bj::na(o[6 * n + i]);
      }
    } else {
      double *cn = f.dbl("cumnet");
      std::copy(o + 7 * n, o + 8 * n, cn);
    }
    // Little's law on published figures: the filing rate over span days,
    // the mean time the pending implies, and the published mean's backlog
    if (!Rf_isNull(span)) {
      if (Rf_isNull(pv)) fail("`span` needs `pend`.");
      SEXP sp = PROTECT(bj::nums(span, "span", -1, 0.0, true));
      ++np;
      if (XLENGTH(sp) != 1 && XLENGTH(sp) != n)
        fail("`span` must be one number or one per period, %d.", n);
      double *fr = f.dbl("frate"), *im = f.dbl("impl");
      for (int i = 0; i < n; ++i) {
        fr[i] = REAL(rv)[i] / REAL(sp)[XLENGTH(sp) == 1 ? 0 : i];
        im[i] = fr[i] > 0 ? REAL(pv)[i] / fr[i] : NA_REAL;
      }
      if (!Rf_isNull(mtime)) {
        SEXP mt = PROTECT(bj::nums(mtime, "mtime", n, 0.0, true));
        ++np;
        double *ra = f.dbl("ratio"), *lb = f.dbl("ltlbk");
        for (int i = 0; i < n; ++i) {
          ra[i] = ISNAN(im[i]) ? NA_REAL : im[i] / REAL(mt)[i];
          lb[i] = fr[i] * REAL(mt)[i];
        }
      }
    } else if (!Rf_isNull(mtime)) {
      fail("`mtime` needs `span` and `pend`.");
    }
    // the floors of published duration bands bound the mean from below
    if (!Rf_isNull(shares)) {
      SEXP dm = Rf_getAttrib(shares, R_DimSymbol);
      if (!clab::isnum(shares) || TYPEOF(dm) != INTSXP || XLENGTH(dm) != 2 || INTEGER(dm)[0] != n) {
        fail("`shares` must be a numeric matrix with one row per period, %d.", n);
      }
      const int k = INTEGER(dm)[1];
      SEXP fl = PROTECT(bj::nums(floors, "floors", k, 0.0));
      SEXP sv = PROTECT(clab::dbl(shares));
      np += 2;
      for (int j = 1; j < k; ++j) {
        if (REAL(fl)[j] < REAL(fl)[j - 1]) fail("`floors` must be increasing.");
      }
      double *bf = f.dbl("floorbd");
      for (int i = 0; i < n; ++i) {
        double s = 0, sw = 0;
        for (int j = 0; j < k; ++j) {
          const double v = REAL(sv)[i + static_cast<R_xlen_t>(n) * j];
          if (!(v >= 0) || !std::isfinite(v)) fail("`shares` must be numbers from 0.");
          s += v * REAL(fl)[j];
          sw += v;
        }
        if (std::fabs(sw - 1) > 0.01) fail("each row of `shares` must sum to 1 (within rounding).");
        bf[i] = s / sw;
      }
    }
    SET_VECTOR_ELT(out, 0, f.done());
  }
  if (!Rf_isNull(mix)) {
    SEXP dm = Rf_getAttrib(mix, R_DimSymbol);
    if (!clab::isnum(mix) || TYPEOF(dm) != INTSXP || XLENGTH(dm) != 2) {
      fail("`mix` must be a numeric matrix of periods by outcomes.");
    }
    const int nt = INTEGER(dm)[0], k = INTEGER(dm)[1];
    if (nt < 1 || k < 1) fail("`mix` must have a period and an outcome.");
    SEXP xv = PROTECT(clab::dbl(mix));
    ++np;
    SEXP sc = R_NilValue;
    if (Rf_isNull(scores)) {
      sc = PROTECT(Rf_allocVector(REALSXP, nt));
      for (int t = 0; t < nt; ++t) REAL(sc)[t] = t + 1;
    } else {
      sc = PROTECT(bj::nums(scores, "scores", nt));
    }
    ++np;
    const long long c = static_cast<long long>(nt) * k;
    SEXP mw = PROTECT(Rf_allocVector(REALSXP, 3 * c + 3 + 2 * static_cast<R_xlen_t>(k)));
    ++np;
    double *sh = REAL(mw), *lo = sh + c, *hi = lo + c, *st = hi + c, *tr = st + 3;
    if (clab_dkmix(REAL(xv), nt, k, REAL(sc), Rf_qchisq(1 - a / k, 1.0, 1, 0), sh, lo, hi, st,
                   tr) != CLAB_OK) {
      fail("`mix` must hold counts from 0, every period with a case.");
    }
    SEXP dn = Rf_getAttrib(mix, R_DimNamesSymbol);
    SEXP rn = Rf_isNull(dn) ? R_NilValue : VECTOR_ELT(dn, 0);
    SEXP cn = Rf_isNull(dn) ? R_NilValue : VECTOR_ELT(dn, 1);
    bj::cols f(6, static_cast<int>(c));
    ++np;
    SEXP fp = PROTECT(Rf_allocVector(STRSXP, c)), fo = PROTECT(Rf_allocVector(STRSXP, c));
    np += 2;
    for (int j = 0; j < k; ++j) {
      for (int t = 0; t < nt; ++t) {
        const long long q = t + static_cast<long long>(nt) * j;
        SET_STRING_ELT(
            fp, q, Rf_isNull(rn) ? Rf_mkChar(std::to_string(t + 1).c_str()) : STRING_ELT(rn, t));
        SET_STRING_ELT(fo, q,
                       Rf_isNull(cn) ? Rf_mkChar(("outcome" + std::to_string(j + 1)).c_str())
                                     : STRING_ELT(cn, j));
      }
    }
    f.put("period", fp);
    f.put("outcome", fo);
    double *fc = f.dbl("count"), *fs = f.dbl("share"), *fl = f.dbl("lo"), *fh = f.dbl("hi");
    std::copy(REAL(xv), REAL(xv) + c, fc);
    std::copy(sh, sh + c, fs);
    std::copy(lo, lo + c, fl);
    std::copy(hi, hi + c, fh);
    // the trend tests, Holm-adjusted over the outcomes that have one
    bj::cols g(6, k);
    ++np;
    SEXP go = PROTECT(Rf_allocVector(STRSXP, k));
    ++np;
    for (int j = 0; j < k; ++j)
      SET_STRING_ELT(go, j, STRING_ELT(fo, static_cast<long long>(nt) * j));
    g.put("outcome", go);
    double *gs = g.dbl("stat"), *gp = g.dbl("p"), *gh = g.dbl("holm"), *gc = g.dbl("change");
    SEXP pw = PROTECT(Rf_allocVector(REALSXP, 2 * static_cast<R_xlen_t>(k)));
    SEXP iw = PROTECT(Rf_allocVector(INTSXP, k));
    np += 2;
    int m = 0;
    for (int j = 0; j < k; ++j) {
      gs[j] = bj::na(tr[j]);
      gp[j] = ISNAN(gs[j]) ? NA_REAL : Rf_pchisq(gs[j], 1.0, 0, 0);
      gc[j] = tr[k + j];
      if (!ISNAN(gp[j])) REAL(pw)[m++] = gp[j];
    }
    clab_pbadj(REAL(pw), m, 2, REAL(pw) + k, INTEGER(iw));
    for (int j = 0, i = 0; j < k; ++j) gh[j] = ISNAN(gp[j]) ? NA_REAL : REAL(pw)[k + i++];
    const double hv[4] = {st[0], st[1], nt > 1 ? Rf_pchisq(st[0], st[1], 0, 0) : NA_REAL,
                          bj::na(st[2])};
    static const char *const hn[] = {"chisq", "df", "p", "cramerv"};
    SEXP mx = PROTECT(Rf_allocVector(VECSXP, 3));
    ++np;
    SET_VECTOR_ELT(mx, 0, f.done());
    SET_VECTOR_ELT(mx, 1, clab::nvec(4, hn, hv));
    SET_VECTOR_ELT(mx, 2, g.done());
    static const char *const xn[] = {"shares", "homog", "trend"};
    SET_VECTOR_ELT(out, 1, bj::named(mx, xn));
  }
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(lv));
  static const char *const nm[] = {"flow", "mix", "level"};
  UNPROTECT(np);
  return bj::named(out, nm);
}

SEXP C_clab_r_lftab(SEXP lo, SEXP hi, SEXP events, SEXP cens, SEXP level, SEXP ci) {
  clab::who = "lftab";
  int np = 0;
  const int k = bj::len(lo, "breaks");
  SEXP dv = PROTECT(bj::nums(events, "events", k, 0.0));
  SEXP c0 = PROTECT(bj::nums(cens, "cens", -1, 0.0));
  np += 2;
  if (XLENGTH(c0) != 1 && XLENGTH(c0) != k)
    fail("`cens` must be one number or one per band, %d.", k);
  SEXP cv = PROTECT(Rf_allocVector(REALSXP, k));
  ++np;
  for (int i = 0; i < k; ++i) REAL(cv)[i] = REAL(c0)[XLENGTH(c0) == 1 ? 0 : i];
  const double lv = clab::one(level, "level", 0.0, 1.0);
  static const char *const cis[] = {"plain", "log", "log-log"};
  const int c = csp::pick(ci, "ci", cis, 3);
  SEXP ow = PROTECT(Rf_allocVector(REALSXP, 8 * static_cast<R_xlen_t>(k) + 1));
  ++np;
  const int rc = clab_lftab(REAL(lo), REAL(hi), REAL(dv), REAL(cv), k,
                            Rf_qnorm5(0.5 + lv / 2, 0.0, 1.0, 1, 0), c, REAL(ow));
  if (rc != CLAB_OK) {
    fail("the bands must be contiguous, only the last open, with cases in them.");
  }
  const double *o = REAL(ow);
  bj::cols f(13, k);
  ++np;
  f.put("lower", lo);
  f.put("upper", hi);
  static const char *const cn[] = {"entering", "atrisk", "q", "surv", "se", "lo", "hi", "hazard"};
  double *e = f.dbl(cn[0]), *ar = f.dbl(cn[1]);
  f.put("events", dv);
  f.put("cens", cv);
  for (int i = 0; i < k; ++i) {
    e[i] = o[i];
    ar[i] = o[k + i];
  }
  for (int j = 2; j < 8; ++j) {
    double *d = f.dbl(cn[j]);
    for (int i = 0; i < k; ++i) d[i] = bj::na(o[j * k + i]);
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  ++np;
  SET_VECTOR_ELT(out, 0, f.done());
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(o[8 * k]));
  SET_VECTOR_ELT(out, 2, Rf_mkString(cis[c]));
  static const char *const nm[] = {"tab", "median", "ci"};
  UNPROTECT(np);
  return bj::named(out, nm);
}

SEXP C_clab_r_trnmat(SEXP state, SEXP id, SEXP time) {
  clab::who = "trnmat";
  int np = 0;
  const int n = bj::len(state, "state");
  csp::lab s = bj::labs(state, "state", n);
  csp::lab p = bj::labs(id, "id", n);
  np += 4;
  SEXP tv = R_NilValue;
  if (!Rf_isNull(time)) {
    tv = PROTECT(bj::nums(time, "time", n));
    ++np;
  }
  // each person's rows in time order (the row order breaks ties)
  std::vector<int> o(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) o[static_cast<std::size_t>(i)] = i;
  const int *pc = INTEGER(p.code);
  const double *tt = Rf_isNull(tv) ? nullptr : REAL(tv);
  std::stable_sort(o.begin(), o.end(), [&](int a, int b) {
    return pc[a] < pc[b] || (pc[a] == pc[b] && tt && tt[a] < tt[b]);
  });
  SEXP tw = PROTECT(Rf_allocVector(INTSXP, 2 * static_cast<R_xlen_t>(n)));
  ++np;
  int *fr = INTEGER(tw), *to = fr + n, m = 0;
  for (int i = 1; i < n; ++i) {
    const int a = o[static_cast<std::size_t>(i - 1)], b = o[static_cast<std::size_t>(i)];
    if (pc[a] != pc[b]) continue;
    fr[m] = INTEGER(s.code)[a];
    to[m] = INTEGER(s.code)[b];
    ++m;
  }
  o.clear();
  o.shrink_to_fit();
  const int k = s.g;
  if (static_cast<long long>(k) * k > 25000000LL) fail("`state` has more than 5000 states.");
  SEXP cnt = PROTECT(Rf_allocMatrix(REALSXP, k, k)), pm = PROTECT(Rf_allocMatrix(REALSXP, k, k));
  SEXP pv = PROTECT(Rf_allocVector(REALSXP, k));
  SEXP wk = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(k) * (k + 1)));
  np += 4;
  double st[3];
  clab_trnmat(fr, to, m, k, REAL(cnt), REAL(pm), REAL(pv), st, REAL(wk));
  SEXP dn = PROTECT(Rf_allocVector(VECSXP, 2));
  ++np;
  SET_VECTOR_ELT(dn, 0, s.lev);
  SET_VECTOR_ELT(dn, 1, s.lev);
  Rf_setAttrib(cnt, R_DimNamesSymbol, dn);
  Rf_setAttrib(pm, R_DimNamesSymbol, dn);
  for (int a = 0; a < k; ++a) REAL(pv)[a] = bj::na(REAL(pv)[a]);
  Rf_setAttrib(pv, R_NamesSymbol, s.lev);
  for (int j = 0; j < 3; ++j) st[j] = bj::na(st[j]);
  static const char *const sn[] = {"moves", "stay", "theil"};
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  ++np;
  SET_VECTOR_ELT(out, 0, cnt);
  SET_VECTOR_ELT(out, 1, pm);
  SET_VECTOR_ELT(out, 2, pv);
  SET_VECTOR_ELT(out, 3, clab::nvec(3, sn, st));
  static const char *const nm[] = {"count", "prob", "pi", "summ"};
  UNPROTECT(np);
  return bj::named(out, nm);
}

}  // extern "C"
