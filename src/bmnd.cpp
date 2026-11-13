// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry point of the Mandela-Rules classification (mndla): checks
// its arguments, calls the kernels and builds the R result. Beta
// quantiles come from R's maths library. Nothing calls back into the R
// evaluator.
#include <cctype>
#include <cmath>
#include <cstring>

#include "bcol.h"
#include "rchk.h"

extern "C" double Rf_qbeta(double, double, double, int, int);

using clab::fail;

namespace {

// "yes" 1, "no" 0 (any case, spaces around ignored), anything else NA.
int word(const char *s) {
  while (*s == ' ') ++s;
  char b[4] = {0, 0, 0, 0};
  int k = 0;
  for (; *s && *s != ' ' && k < 4; ++s, ++k) b[k] = static_cast<char>(std::tolower(*s));
  while (*s == ' ') ++s;
  if (*s || k == 4) return NA_INTEGER;
  if (std::strcmp(b, "yes") == 0) return 1;
  return std::strcmp(b, "no") == 0 ? 0 : NA_INTEGER;
}

// TRUE/FALSE, 0/1 or "yes"/"no" (any case) per row, as 0 or 1 into out.
void yes(SEXP x, const char *nm, int n, int *out) {
  clab::unlab(x, nm);
  if (XLENGTH(x) != n) fail("`%s` must have one value per stay, %d.", nm, n);
  for (int i = 0; i < n; ++i) {
    int v = NA_INTEGER;
    if (TYPEOF(x) == LGLSXP) {
      v = LOGICAL(x)[i];
    } else if (TYPEOF(x) == INTSXP && !Rf_isFactor(x)) {
      v = INTEGER(x)[i] == NA_INTEGER ? NA_INTEGER : INTEGER(x)[i] != 0;
    } else if (TYPEOF(x) == REALSXP) {
      v = std::isnan(REAL(x)[i]) ? NA_INTEGER : REAL(x)[i] != 0;
    } else if (TYPEOF(x) == STRSXP && STRING_ELT(x, i) != NA_STRING) {
      v = word(CHAR(STRING_ELT(x, i)));
    }
    if (v == NA_INTEGER) {
      fail("`%s` must be TRUE or FALSE, 0 or 1, or \"yes\" or \"no\", with no missing values.", nm);
    }
    out[i] = v;
  }
}

}  // namespace

extern "C" {

SEXP C_clab_r_mndla(SEXP days, SEXP id, SEXP by, SEXP w, SEXP thr, SEXP hours, SEXP miss,
                    SEXP alert, SEXP amin, SEXP elig, SEXP excl, SEXP level) {
  clab::who = "mndla";
  int np = 0;
  SEXP dv = PROTECT(bj::nums(days, "days", -1, 0.0));
  ++np;
  const int n = bj::len(dv, "days");
  const double th = bj::inr(thr, "thr", 0.0, R_PosInf);
  const double lv = clab::one(level, "level", 0.0, 1.0);
  SEXP hv = R_NilValue, mv = R_NilValue;
  if (Rf_isNull(hours) != Rf_isNull(miss)) fail("`hours` and `miss` go together.");
  if (!Rf_isNull(hours)) {
    hv = PROTECT(bj::nums(hours, "hours", n, 0.0));
    mv = PROTECT(bj::nums(miss, "miss", n, 0.0));
    np += 2;
    for (int i = 0; i < n; ++i) {
      if (REAL(mv)[i] > 100) fail("`miss` must be percentages from 0 to 100.");
    }
  }
  SEXP iw = PROTECT(Rf_allocVector(INTSXP, 4 * static_cast<R_xlen_t>(n)));
  ++np;
  int *nal = INTEGER(iw), *el = nal + n, *ex = el + n, *cls = ex + n;
  const int am = clab::whole(amin, "amin", 1);
  bool hasal = false;
  if (!Rf_isNull(alert)) {
    if (!Rf_isNull(hours)) fail("`alert` goes with the duration rule, not with `hours`.");
    hasal = true;
    for (int i = 0; i < n; ++i) nal[i] = 0;
    // one alert, or several: a data frame, list or matrix of columns
    SEXP dm = Rf_getAttrib(alert, R_DimSymbol);
    const bool m = TYPEOF(dm) == INTSXP && XLENGTH(dm) == 2;
    const int k = TYPEOF(alert) == VECSXP ? bj::len(alert, "alert") : m ? INTEGER(dm)[1] : 1;
    for (int j = 0; j < k; ++j) {
      SEXP c = TYPEOF(alert) == VECSXP ? VECTOR_ELT(alert, j) : alert;
      if (m) {
        if (INTEGER(dm)[0] != n) fail("`alert` must have one row per stay, %d.", n);
        SEXP s = PROTECT(Rf_allocVector(TYPEOF(alert), n));
        for (int i = 0; i < n; ++i) {
          const R_xlen_t q = i + static_cast<R_xlen_t>(n) * j;
          if (TYPEOF(alert) == LGLSXP)
            LOGICAL(s)[i] = LOGICAL(alert)[q];
          else if (TYPEOF(alert) == REALSXP)
            REAL(s)[i] = REAL(alert)[q];
          else if (TYPEOF(alert) == STRSXP)
            SET_STRING_ELT(s, i, STRING_ELT(alert, q));
          else
            INTEGER(s)[i] = INTEGER(alert)[q];
        }
        yes(s, "alert", n, cls);
        UNPROTECT(1);
      } else {
        yes(c, "alert", n, cls);
      }
      for (int i = 0; i < n; ++i) nal[i] += cls[i];
    }
  }
  if (!Rf_isNull(elig)) yes(elig, "elig", n, el);
  if (!Rf_isNull(excl)) yes(excl, "excl", n, ex);
  clab_mndla(REAL(dv), Rf_isNull(hv) ? nullptr : REAL(hv), Rf_isNull(mv) ? nullptr : REAL(mv),
             hasal ? nal : nullptr, am, Rf_isNull(elig) ? nullptr : el,
             Rf_isNull(excl) ? nullptr : ex, n, th, cls);
  // the groups, then all
  SEXP gcode = PROTECT(Rf_allocVector(INTSXP, n));
  ++np;
  int ng = 1;
  SEXP glev = R_NilValue;
  if (!Rf_isNull(by)) {
    csp::lab g = bj::labs(by, "by", n);
    np += 2;
    ng = g.g;
    glev = g.lev;
    std::memcpy(INTEGER(gcode), INTEGER(g.code), sizeof(int) * static_cast<std::size_t>(n));
  } else {
    for (int i = 0; i < n; ++i) INTEGER(gcode)[i] = 0;
  }
  const int *idc = nullptr;
  if (!Rf_isNull(id)) {
    csp::lab p = bj::labs(id, "id", n);
    np += 2;
    idc = INTEGER(p.code);
  }
  SEXP wv = R_NilValue;
  if (!Rf_isNull(w)) {
    wv = PROTECT(bj::nums(w, "w", n, 0.0));
    ++np;
  }
  SEXP ow = PROTECT(Rf_allocVector(REALSXP, 7 * (static_cast<R_xlen_t>(ng) + 1)));
  SEXP oi = PROTECT(Rf_allocVector(INTSXP, 2 * static_cast<R_xlen_t>(n)));
  np += 2;
  clab_mndtab(cls, idc, INTEGER(gcode), REAL(dv), Rf_isNull(wv) ? nullptr : REAL(wv),
              Rf_isNull(elig) ? nullptr : el, n, ng, th, REAL(ow), INTEGER(oi));
  // one row per group and denominator: the groups (none without by),
  // then all
  const int nd = idc ? 3 : 1, g0 = Rf_isNull(by) ? 1 : 0;
  const int nr = (ng + 1 - g0) * nd;
  bj::cols f(9, nr);
  ++np;
  SEXP gs = PROTECT(Rf_allocVector(STRSXP, nr)), ds = PROTECT(Rf_allocVector(STRSXP, nr));
  np += 2;
  f.put("group", gs);
  f.put("den", ds);
  double *fn = f.dbl("n"), *fk = f.dbl("k"), *fs = f.dbl("share"), *fl = f.dbl("lo"),
         *fh = f.dbl("hi");
  double *fsol = Rf_isNull(hv) ? nullptr : f.dbl("sol"),
         *ftor = Rf_isNull(hv) ? nullptr : f.dbl("tor");
  static const char *const dn[] = {"stays", "people", "cumul"};
  const double a = 1 - lv;
  for (int g = g0, r = 0; g <= ng; ++g) {
    const double *o = REAL(ow) + 7 * (g == ng ? ng : g);
    for (int d = 0; d < nd; ++d, ++r) {
      SET_STRING_ELT(gs, r, g == ng ? Rf_mkChar("all") : STRING_ELT(glev, g));
      SET_STRING_ELT(ds, r, Rf_mkChar(dn[d]));
      const double nn = d == 0 ? o[0] : o[4], kk = d == 0 ? o[1] : d == 1 ? o[5] : o[6];
      fn[r] = nn;
      fk[r] = kk;
      fs[r] = nn > 0 ? kk / nn : NA_REAL;
      // the exact (Clopper-Pearson) interval, for whole counts
      const bool whole = nn > 0 && nn == std::floor(nn) && kk == std::floor(kk);
      fl[r] = !whole ? NA_REAL : kk > 0 ? Rf_qbeta(a / 2, kk, nn - kk + 1, 1, 0) : 0.0;
      fh[r] = !whole ? NA_REAL : kk < nn ? Rf_qbeta(1 - a / 2, kk + 1, nn - kk, 1, 0) : 1.0;
      if (fsol) {
        fsol[r] = d == 0 ? o[2] : NA_REAL;
        ftor[r] = d == 0 ? o[3] : NA_REAL;
      }
    }
  }
  SEXP cl = PROTECT(Rf_allocVector(INTSXP, n));
  ++np;
  std::memcpy(INTEGER(cl), cls, sizeof(int) * static_cast<std::size_t>(n));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  ++np;
  SET_VECTOR_ELT(out, 0, f.done());
  SET_VECTOR_ELT(out, 1, cl);
  if (hasal) {
    SEXP na = PROTECT(Rf_allocVector(INTSXP, n));
    ++np;
    std::memcpy(INTEGER(na), nal, sizeof(int) * static_cast<std::size_t>(n));
    SET_VECTOR_ELT(out, 2, na);
  }
  SET_VECTOR_ELT(out, 3, Rf_mkString(Rf_isNull(hv) ? "duration" : "rule44"));
  static const char *const nm[] = {"tab", "cls", "alerts", "rule"};
  UNPROTECT(np);
  return bj::named(out, nm);
}

}  // extern "C"
