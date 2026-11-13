// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the life-course and pipeline family (agemix,
// cumrsk, jpipe): each checks its arguments, calls the kernels and builds
// the R result. Normal quantiles come from R's maths library. Nothing
// calls back into the R evaluator.
#include <algorithm>
#include <cmath>
#include <vector>

#include "bcol.h"
#include "rchk.h"

using clab::fail;

namespace {

// A numeric matrix or data frame of numeric columns as a column-major
// double matrix, with its dimensions and column names (protected; the
// caller unprotects one).
SEXP mat(SEXP x, const char *nm, int *nr, int *nc, SEXP *cn) {
  *cn = R_NilValue;
  if (TYPEOF(x) == VECSXP) {
    *nc = bj::len(x, nm);
    if (*nc < 1) fail("`%s` must have at least one column.", nm);
    *nr = bj::len(VECTOR_ELT(x, 0), nm);
    SEXP m = PROTECT(Rf_allocMatrix(REALSXP, *nr, *nc));
    for (int j = 0; j < *nc; ++j) {
      SEXP c = PROTECT(bj::nums(VECTOR_ELT(x, j), nm, *nr));
      std::copy(REAL(c), REAL(c) + *nr, REAL(m) + static_cast<R_xlen_t>(*nr) * j);
      UNPROTECT(1);
    }
    *cn = Rf_getAttrib(x, R_NamesSymbol);
    UNPROTECT(1);
    return PROTECT(m);
  }
  SEXP dm = Rf_getAttrib(x, R_DimSymbol);
  if (!clab::isnum(x) || TYPEOF(dm) != INTSXP || XLENGTH(dm) != 2) {
    fail("`%s` must be a numeric matrix or a data frame of numbers.", nm);
  }
  *nr = INTEGER(dm)[0];
  *nc = INTEGER(dm)[1];
  SEXP dn = Rf_getAttrib(x, R_DimNamesSymbol);
  if (!Rf_isNull(dn)) *cn = VECTOR_ELT(dn, 1);
  SEXP v = PROTECT(clab::dcopy(PROTECT(clab::dbl(x))));
  for (R_xlen_t k = 0; k < XLENGTH(v); ++k) {
    if (!std::isfinite(REAL(v)[k])) fail("`%s` must hold finite numbers.", nm);
  }
  Rf_setAttrib(v, R_DimSymbol, dm);
  UNPROTECT(2);
  return PROTECT(v);
}

}  // namespace

extern "C" {

SEXP C_clab_r_agemix(SEXP curves, SEXP shares, SEXP ages, SEXP panel, SEXP nperm, SEXP seed) {
  clab::who = "agemix";
  int np = 0;
  if (Rf_isNull(curves) && Rf_isNull(panel)) fail("give `curves` and `shares`, or `panel`.");
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  ++np;
  if (!Rf_isNull(curves)) {
    int na, nt;
    SEXP cn;
    SEXP cv = mat(curves, "curves", &na, &nt, &cn);
    ++np;
    if (na < 1) fail("`curves` must have at least one age.");
    if (nt > 30) fail("`curves` may have at most 30 types.");
    SEXP w = PROTECT(bj::nums(shares, "shares", nt, 0.0));
    ++np;
    SEXP av = R_NilValue;
    if (!Rf_isNull(ages)) {
      av = PROTECT(bj::nums(ages, "ages", na));
    } else {
      av = PROTECT(Rf_allocVector(REALSXP, na));
      for (int a = 0; a < na; ++a) REAL(av)[a] = a + 1;
    }
    ++np;
    bj::cols f(nt + 2, na);
    ++np;
    f.put("age", av);
    double *ag = f.dbl("agg");
    SEXP cov = PROTECT(Rf_allocMatrix(REALSXP, na, na));
    ++np;
    if (clab_agemix(REAL(cv), na, nt, REAL(w), ag, REAL(cov)) != CLAB_OK) {
      fail("`shares` must sum to 1.");
    }
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, nt));
    ++np;
    for (int k = 0; k < nt; ++k) {
      if (!Rf_isNull(cn) && STRING_ELT(cn, k) != NA_STRING && CHAR(STRING_ELT(cn, k))[0]) {
        SET_STRING_ELT(nms, k, STRING_ELT(cn, k));
      } else {
        SET_STRING_ELT(nms, k, Rf_mkChar(("type" + std::to_string(k + 1)).c_str()));
      }
    }
    // the type columns keep the names that are free
    for (int k = 0; k < nt; ++k) {
      SEXP c = PROTECT(Rf_allocVector(REALSXP, na));
      std::copy(REAL(cv) + static_cast<R_xlen_t>(na) * k,
                REAL(cv) + static_cast<R_xlen_t>(na) * (k + 1), REAL(c));
      f.put(CHAR(STRING_ELT(nms, k)), c);
      UNPROTECT(1);
    }
    SEXP dn = PROTECT(Rf_allocVector(VECSXP, 2));
    ++np;
    SEXP al = PROTECT(Rf_coerceVector(av, STRSXP));
    ++np;
    SET_VECTOR_ELT(dn, 0, al);
    SET_VECTOR_ELT(dn, 1, al);
    Rf_setAttrib(cov, R_DimNamesSymbol, dn);
    SET_VECTOR_ELT(out, 0, f.done());
    SET_VECTOR_ELT(out, 1, cov);
  }
  if (!Rf_isNull(panel)) {
    int n, k;
    SEXP cn;
    SEXP pv = mat(panel, "panel", &n, &k, &cn);
    ++np;
    if (k != 2) fail("`panel` must have two columns, the outcomes at two ages.");
    if (n < 3) fail("`panel` must have at least three individuals.");
    const int nps = clab::whole(nperm, "nperm", 0);
    const unsigned int sd = csp::seed(seed);
    SEXP wk = PROTECT(Rf_allocVector(REALSXP, n));
    ++np;
    double o[2];
    clab_agepan(REAL(pv), REAL(pv) + n, n, nps, sd, o, REAL(wk));
    static const char *const pn[] = {"cov", "p", "n", "nperm"};
    const double pr[4] = {o[0], nps > 0 ? o[1] : NA_REAL, static_cast<double>(n),
                          static_cast<double>(nps)};
    SET_VECTOR_ELT(out, 2, clab::nvec(4, pn, pr));
  }
  static const char *const nm[] = {"curve", "cov", "panel"};
  UNPROTECT(np);
  return bj::named(out, nm);
}

SEXP C_clab_r_cumrsk(SEXP age, SEXP ev, SEXP pop, SEXP atrisk, SEXP dth, SEXP level) {
  clab::who = "cumrsk";
  int np = 0;
  SEXP av = PROTECT(bj::nums(age, "age"));
  ++np;
  const int n = bj::len(av, "age");
  for (int i = 1; i < n; ++i) {
    if (!(REAL(av)[i] > REAL(av)[i - 1])) fail("`age` must be increasing.");
  }
  SEXP ev2 = PROTECT(bj::nums(ev, "ev", n, 0.0));
  SEXP pv = PROTECT(bj::nums(pop, "pop", n, 0.0, true));
  np += 2;
  SEXP rv = R_NilValue, dv = R_NilValue;
  if (!Rf_isNull(atrisk)) {
    rv = PROTECT(bj::nums(atrisk, "atrisk", n, 0.0, true));
    ++np;
  }
  if (!Rf_isNull(dth)) {
    dv = PROTECT(bj::nums(dth, "dth", n, 0.0));
    ++np;
    for (int i = 0; i < n; ++i) {
      if (REAL(dv)[i] >= 1) fail("`dth` must be death rates from 0 to below 1.");
    }
  }
  const double lv = clab::one(level, "level", 0.0, 1.0);
  const double z = Rf_qnorm5(0.5 + lv / 2, 0.0, 1.0, 1, 0);
  bj::cols f(6, n);
  ++np;
  f.put("age", av);
  double *h = f.dbl("hazard"), *nr = f.dbl("atrisk"), *cr = f.dbl("cumrisk"), *lo = f.dbl("lo"),
         *hi = f.dbl("hi");
  SEXP vw = PROTECT(Rf_allocVector(REALSXP, n));
  ++np;
  if (clab_cumrsk(REAL(av), REAL(ev2), REAL(pv), Rf_isNull(rv) ? nullptr : REAL(rv),
                  Rf_isNull(dv) ? nullptr : REAL(dv), n, h, nr, cr, REAL(vw)) != CLAB_OK) {
    fail("more first events (and deaths) than people at risk at some age.");
  }
  for (int i = 0; i < n; ++i) {
    // the interval of 1 - S on the log scale of S
    const double s = 1 - cr[i], sv = std::sqrt(REAL(vw)[i]);
    lo[i] = std::isnan(sv) ? NA_REAL : 1 - s * std::exp(z * sv);
    hi[i] = std::isnan(sv) ? NA_REAL : 1 - s * std::exp(-z * sv);
  }
  SEXP out = PROTECT(f.done());
  ++np;
  UNPROTECT(np);
  return out;
}

SEXP C_clab_r_jpipe(SEXP num, SEXP den, SEXP win, SEXP lag, SEXP grp) {
  clab::who = "jpipe";
  int np = 0;
  SEXP nv = PROTECT(bj::nums(num, "num", -1, 0.0));
  ++np;
  const int n = bj::len(nv, "num");
  SEXP dv = PROTECT(bj::nums(den, "den", n, 0.0));
  ++np;
  const int w = clab::whole(win, "win", 1), lg = clab::whole(lag, "lag", 0);
  // the rows of each group in their order, the groups in label order
  SEXP code = R_NilValue, lev = R_NilValue;
  if (!Rf_isNull(grp)) {
    csp::lab g = bj::labs(grp, "grp", n);
    np += 2;
    code = g.code;
    lev = g.lev;
  }
  std::vector<int> o(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) o[static_cast<std::size_t>(i)] = i;
  if (!Rf_isNull(code)) {
    const int *c = INTEGER(code);
    std::stable_sort(o.begin(), o.end(), [&](int a, int b) { return c[a] < c[b]; });
  }
  SEXP ws = PROTECT(Rf_allocVector(REALSXP, 2 * static_cast<R_xlen_t>(n)));
  ++np;
  double *sn = REAL(ws), *sd = sn + n;
  for (int i = 0; i < n; ++i) {
    sn[i] = REAL(nv)[o[static_cast<std::size_t>(i)]];
    sd[i] = REAL(dv)[o[static_cast<std::size_t>(i)]];
  }
  bj::cols f(5, n);
  ++np;
  SEXP gl = PROTECT(Rf_allocVector(STRSXP, n));
  ++np;
  f.put("group", gl);
  int *ix = f.itg("index");
  double *on = f.dbl("num"), *od = f.dbl("den"), *rt = f.dbl("rate");
  for (int a = 0; a < n;) {
    int b = a;
    const int ca = Rf_isNull(code) ? 0 : INTEGER(code)[o[static_cast<std::size_t>(a)]];
    while (b < n && (Rf_isNull(code) ? 0 : INTEGER(code)[o[static_cast<std::size_t>(b)]]) == ca)
      ++b;
    clab_jpipe(sn + a, sd + a, b - a, w, lg, on + a, od + a, rt + a);
    for (int i = a; i < b; ++i) {
      SET_STRING_ELT(gl, i, Rf_isNull(lev) ? Rf_mkChar("all") : STRING_ELT(lev, ca));
      ix[i] = i - a + 1;
      on[i] = bj::na(on[i]);
      od[i] = bj::na(od[i]);
      rt[i] = bj::na(rt[i]);
    }
    a = b;
  }
  SEXP out = PROTECT(f.done());
  ++np;
  UNPROTECT(np);
  return out;
}

}  // extern "C"
