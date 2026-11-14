// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry point of the split-population survival model (spcure,
// src/spc.cpp): it checks the arguments, calls the kernel and builds the
// R result. Nothing calls back into the R evaluator.
#include <cmath>
#include <string>
#include <vector>

#include "bcol.h"
#include "rchk.h"
#include "spin.h"

extern "C" double Rf_pnorm5(double, double, double, int, int);

using clab::fail;

namespace {

// The covariates: NULL, a numeric matrix, a data frame or list of numeric
// columns, or one numeric vector; n rows, finite; their names.
int xcov(SEXP x, int n, std::vector<double> &v, std::vector<std::string> &nm) {
  v.clear();
  nm.clear();
  if (Rf_isNull(x)) return 0;
  if (Rf_isMatrix(x)) {
    if (!clab::isnum(x) || Rf_nrows(x) != n) fail("`x` must be a numeric matrix with one row per person, %d.", n);
    SEXP dn = Rf_getAttrib(x, R_DimNamesSymbol);
    SEXP cn = Rf_isNull(dn) ? R_NilValue : VECTOR_ELT(dn, 1);
    SEXP xd = PROTECT(clab::dbl(x));
    v.assign(REAL(xd), REAL(xd) + XLENGTH(xd));
    UNPROTECT(1);
    for (double e : v) {
      if (!std::isfinite(e)) fail("`x` must hold finite numbers with no missing values.");
    }
    for (int j = 0; j < Rf_ncols(x); ++j) {
      const bool has = !Rf_isNull(cn) && STRING_ELT(cn, j) != NA_STRING && CHAR(STRING_ELT(cn, j))[0] != '\0';
      nm.push_back(has ? Rf_translateCharUTF8(STRING_ELT(cn, j)) : "x" + std::to_string(j + 1));
    }
    return Rf_ncols(x);
  }
  std::vector<SEXP> cs;
  SEXP names = R_NilValue;
  if (TYPEOF(x) == VECSXP) {
    for (R_xlen_t j = 0; j < XLENGTH(x); ++j) cs.push_back(VECTOR_ELT(x, j));
    names = Rf_getAttrib(x, R_NamesSymbol);
  } else {
    cs.push_back(x);
  }
  for (std::size_t j = 0; j < cs.size(); ++j) {
    SEXP c = PROTECT(clab::num(cs[j], "x", n));
    v.insert(v.end(), REAL(c), REAL(c) + n);
    UNPROTECT(1);
    const R_xlen_t jj = static_cast<R_xlen_t>(j);
    const bool has = !Rf_isNull(names) && STRING_ELT(names, jj) != NA_STRING && CHAR(STRING_ELT(names, jj))[0] != '\0';
    nm.push_back(has ? Rf_translateCharUTF8(STRING_ELT(names, jj)) : "x" + std::to_string(j + 1));
  }
  return static_cast<int>(cs.size());
}

}  // namespace

extern "C" {

SEXP C_clab_r_spcure(SEXP time, SEXP event, SEXP x, SEXP dist, SEXP level, SEXP maxit) {
  clab::who = "spcure";
  SEXP tv = PROTECT(bj::nums(time, "time", -1, 0.0, true));
  const int n = bj::len(tv, "time");
  SEXP ev = PROTECT(csp::bin01(event, "event", n));
  int ne = 0;
  for (int i = 0; i < n; ++i) ne += INTEGER(ev)[i];
  if (ne < 1 || ne >= n) fail("`event` must hold failures and censored times both.");
  std::vector<double> xv;
  std::vector<std::string> xn;
  const int p = xcov(x, n, xv, xn);
  static const char *const dn[] = {"weibull", "exp", "lognormal"};
  const int d = csp::pick(dist, "dist", dn, 3);
  const double lv = clab::one(level, "level", 0.0, 1.0);
  const int mi = clab::whole(maxit, "maxit", 1);
  std::vector<double> evd(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) evd[static_cast<std::size_t>(i)] = INTEGER(ev)[i];
  const int q = p + 1, np = 2 * q + (d != 1);
  std::vector<double> par(static_cast<std::size_t>(np)), cov(static_cast<std::size_t>(np) * np), st(5);
  const int rc = clab_spcure(REAL(tv), evd.data(), p ? xv.data() : nullptr, n, p, d, mi, par.data(), cov.data(),
                             st.data());
  clab::ichk(rc);
  bj::kc(rc, "the likelihood could not be maximised from its start.");
  const double z = Rf_qnorm5(0.5 + 0.5 * lv, 0.0, 1.0, 1, 0);
  bj::cols f(8, np);
  SEXP part = PROTECT(Rf_allocVector(STRSXP, np)), term = PROTECT(Rf_allocVector(STRSXP, np));
  for (int k = 0; k < np; ++k) {
    const int j = k < 2 * q ? k % q : 0;
    SET_STRING_ELT(part, k, Rf_mkChar(k < q ? "fail" : k < 2 * q ? "time" : "shape"));
    const std::string tn = k >= 2 * q ? (d == 2 ? "log(sigma)" : "log(shape)") : j == 0 ? "(Intercept)" : xn[static_cast<std::size_t>(j - 1)];
    SET_STRING_ELT(term, k, Rf_mkCharCE(tn.c_str(), CE_UTF8));
  }
  f.put("part", part);
  f.put("term", term);
  double *es = f.dbl("est"), *se = f.dbl("se"), *lo = f.dbl("lo"), *hi = f.dbl("hi"), *zz = f.dbl("z"), *pv = f.dbl("p");
  for (int k = 0; k < np; ++k) {
    const double s = std::sqrt(cov[static_cast<std::size_t>(k) * np + k]);
    es[k] = par[static_cast<std::size_t>(k)];
    se[k] = bj::na(s);
    lo[k] = bj::na(es[k] - z * s);
    hi[k] = bj::na(es[k] + z * s);
    zz[k] = bj::na(es[k] / s);
    pv[k] = std::isnan(s) ? NA_REAL : 2 * Rf_pnorm5(-std::fabs(es[k] / s), 0.0, 1.0, 1, 0);
  }
  static const char *const sn[] = {"loglik", "aic", "iter", "conv", "pfail"};
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out, 0, f.done());
  SET_VECTOR_ELT(out, 1, clab::nvec(5, sn, st.data()));
  SET_VECTOR_ELT(out, 2, Rf_mkString(dn[d]));
  SET_VECTOR_ELT(out, 3, Rf_ScalarInteger(n));
  static const char *const on[] = {"coef", "fit", "dist", "n"};
  bj::named(out, on);
  Rf_classgets(out, PROTECT(Rf_mkString("spcure")));
  UNPROTECT(7);
  return out;
}

}  // extern "C"
