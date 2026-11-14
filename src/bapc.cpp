// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry point of the age-period-cohort model (apcie, src/apc.cpp):
// it checks the arguments, calls the kernel and builds the R result.
// Nothing calls back into the R evaluator.
#include <cmath>
#include <string>
#include <vector>

#include "bcol.h"
#include "rchk.h"

using clab::fail;

namespace {

// The labels of dimension d of the matrix m, or 1, 2, ...
std::string dlab(SEXP m, int d, int k) {
  SEXP dn = Rf_getAttrib(m, R_DimNamesSymbol);
  if (!Rf_isNull(dn) && !Rf_isNull(VECTOR_ELT(dn, d))) {
    SEXP s = VECTOR_ELT(dn, d);
    if (STRING_ELT(s, k) != NA_STRING) return Rf_translateCharUTF8(STRING_ELT(s, k));
  }
  return std::to_string(k + 1);
}

// The cells of a numeric matrix as doubles: finite, from 0 (above 0 when
// open), with num()'s words.
SEXP cells(SEXP m, const char *nm, bool open) {
  if (!clab::isnum(m)) fail("`%s` must be finite numbers with no missing values.", nm);
  SEXP v = PROTECT(clab::dbl(m));
  for (R_xlen_t i = 0; i < XLENGTH(v); ++i) {
    if (!std::isfinite(REAL(v)[i])) fail("`%s` must be finite numbers with no missing values.", nm);
  }
  for (R_xlen_t i = 0; i < XLENGTH(v); ++i) {
    const double e = REAL(v)[i];
    if (open ? e <= 0 : e < 0) fail("`%s` must be %s0.", nm, open ? "greater than " : "at least ");
  }
  UNPROTECT(1);
  return v;
}

}  // namespace

extern "C" {

SEXP C_clab_r_apcie(SEXP y, SEXP expo, SEXP level, SEXP maxit) {
  clab::who = "apcie";
  if (!Rf_isMatrix(y) || !clab::isnum(y)) fail("`y` must be a numeric matrix of counts, ages by periods.");
  const int na = Rf_nrows(y), np = Rf_ncols(y);
  if (na < 3 || np < 3) fail("`y` must have three ages and three periods at least.");
  if (static_cast<double>(na) * np > 1e6) fail("`y` has more than 1e6 cells.");
  SEXP yv = PROTECT(cells(y, "y", false));
  double tot = 0;
  for (R_xlen_t i = 0; i < XLENGTH(yv); ++i) tot += REAL(yv)[i];
  if (!(tot > 0)) fail("`y` has no events: every count is 0.");
  if (!Rf_isMatrix(expo) || Rf_nrows(expo) != na || Rf_ncols(expo) != np) {
    fail("`expo` must be a matrix of the dimensions of `y`, %d by %d.", na, np);
  }
  SEXP ev = PROTECT(cells(expo, "expo", true));
  const double lv = clab::one(level, "level", 0.0, 1.0);
  const int mi = clab::whole(maxit, "maxit", 1);
  const int nc = na + np - 1, ne = 1 + na + np + nc;
  std::vector<double> est(static_cast<std::size_t>(ne)), se(static_cast<std::size_t>(ne)), st(5);
  const int rc = clab_apcie(REAL(yv), REAL(ev), na, np, mi, est.data(), se.data(), st.data());
  clab::ichk(rc);
  bj::kc(rc, "the Poisson model could not be fitted to these counts.");
  const double z = Rf_qnorm5(0.5 + 0.5 * lv, 0.0, 1.0, 1, 0);
  bj::cols f(6, ne);
  SEXP tm = PROTECT(Rf_allocVector(STRSXP, ne)), lb = PROTECT(Rf_allocVector(STRSXP, ne));
  SET_STRING_ELT(tm, 0, Rf_mkChar("(Intercept)"));
  SET_STRING_ELT(lb, 0, Rf_mkChar(""));
  for (int k = 1; k < ne; ++k) {
    const bool ag = k <= na, pe = !ag && k <= na + np;
    SET_STRING_ELT(tm, k, Rf_mkChar(ag ? "age" : pe ? "period" : "cohort"));
    const std::string s = ag ? dlab(y, 0, k - 1) : pe ? dlab(y, 1, k - 1 - na) : std::to_string(k - na - np);
    SET_STRING_ELT(lb, k, Rf_mkCharCE(s.c_str(), CE_UTF8));
  }
  f.put("term", tm);
  f.put("level", lb);
  double *es = f.dbl("est"), *ss = f.dbl("se"), *lo = f.dbl("lo"), *hi = f.dbl("hi");
  for (int k = 0; k < ne; ++k) {
    es[k] = est[static_cast<std::size_t>(k)];
    ss[k] = se[static_cast<std::size_t>(k)];
    lo[k] = es[k] - z * ss[k];
    hi[k] = es[k] + z * ss[k];
  }
  static const char *const sn[] = {"deviance", "df", "loglik", "iter", "conv"};
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, f.done());
  SET_VECTOR_ELT(out, 1, clab::nvec(5, sn, st.data()));
  const double dm[3] = {static_cast<double>(na), static_cast<double>(np), static_cast<double>(nc)};
  static const char *const dn[] = {"ages", "periods", "cohorts"};
  SET_VECTOR_ELT(out, 2, clab::nvec(3, dn, dm));
  static const char *const on[] = {"coef", "fit", "dims"};
  bj::named(out, on);
  Rf_classgets(out, PROTECT(Rf_mkString("apcie")));
  UNPROTECT(7);
  return out;
}

}  // extern "C"
