// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of data description (src/dsc.cpp). The R wrappers
// (R/dsc.R) check the arguments and take the data frames apart. Nothing
// calls back into the R evaluator.
#include <cmath>
#include <vector>

#include "bcol.h"
#include "rchk.h"

namespace {

int nr(SEXP x) { return Rf_isMatrix(x) ? Rf_nrows(x) : static_cast<int>(XLENGTH(x)); }
int nc(SEXP x) { return Rf_isMatrix(x) ? Rf_ncols(x) : 1; }

SEXP list2(SEXP a, SEXP b) {
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(r, 0, a);
  SET_VECTOR_ELT(r, 1, b);
  UNPROTECT(1);
  return r;
}

}  // namespace

extern "C" {

SEXP C_clab_r_hbin(SEXP x, SEXP bins) {
  const int b = Rf_asInteger(bins);
  SEXP o = PROTECT(Rf_allocVector(INTSXP, b));
  clab_hbin(REAL(x), static_cast<int>(clab::xlen(x, "x")), b, INTEGER(o));
  UNPROTECT(1);
  return o;
}

SEXP C_clab_r_freq(SEXP code, SEXP nlev, SEXP nna, SEXP nall, SEXP sort, SEXP maxlvl, SEXP fold) {
  const int n = static_cast<int>(clab::xlen(code, "code")), nl = Rf_asInteger(nlev), na = Rf_asInteger(nna);
  const double all = Rf_asReal(nall), valid = all - na;
  std::vector<int> lev(static_cast<size_t>(nl) + 2), cnt(static_cast<size_t>(nl) + 2);
  int k = clab_freq(INTEGER(code), n, nl, Rf_asLogical(sort), Rf_asInteger(maxlvl), Rf_asLogical(fold), lev.data(),
                    cnt.data());
  if (na > 0) {
    lev[static_cast<size_t>(k)] = NA_INTEGER;
    cnt[static_cast<size_t>(k++)] = na;
  }
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 4));
  SEXP l = Rf_allocVector(INTSXP, k);
  SET_VECTOR_ELT(r, 0, l);
  SEXP c = Rf_allocVector(INTSXP, k);
  SET_VECTOR_ELT(r, 1, c);
  SEXP p = Rf_allocVector(REALSXP, k);
  SET_VECTOR_ELT(r, 2, p);
  SEXP v = Rf_allocVector(REALSXP, k);
  SET_VECTOR_ELT(r, 3, v);
  for (int i = 0; i < k; ++i) {
    INTEGER(l)[i] = lev[static_cast<size_t>(i)];
    INTEGER(c)[i] = cnt[static_cast<size_t>(i)];
    REAL(p)[i] = 100 * static_cast<double>(cnt[static_cast<size_t>(i)]) / all;
    REAL(v)[i] = lev[static_cast<size_t>(i)] == NA_INTEGER || !(valid > 0)
                     ? NA_REAL
                     : 100 * static_cast<double>(cnt[static_cast<size_t>(i)]) / valid;
  }
  UNPROTECT(1);
  return r;
}

SEXP C_clab_r_pcor(SEXP x, SEXP pear, SEXP minprs) {
  const int n = Rf_nrows(x), p = Rf_ncols(x), k = p * (p - 1) / 2;
  SEXP c = PROTECT(Rf_allocVector(REALSXP, k));
  SEXP m = PROTECT(Rf_allocVector(INTSXP, k));
  const int mp = Rf_asInteger(minprs);
  clab_pcor(REAL(x), n, p, Rf_asLogical(pear), mp, REAL(c), INTEGER(m));
  // too few pairs is NA; no spread NaN, as the reference gives them
  for (int i = 0; i < k; ++i) {
    if (INTEGER(m)[i] < mp) REAL(c)[i] = NA_REAL;
  }
  SEXP r = list2(c, m);
  UNPROTECT(2);
  return r;
}

SEXP C_clab_r_rowany(SEXP na) {
  const int n = nr(na), p = nc(na);
  SEXP o = PROTECT(Rf_allocVector(INTSXP, n));
  for (int i = 0; i < n; ++i) {
    int c = 0;
    for (int j = 0; j < p; ++j) c += LOGICAL(na)[static_cast<R_xlen_t>(j) * n + i] == 0;
    INTEGER(o)[i] = c;
  }
  UNPROTECT(1);
  return o;
}

SEXP C_clab_r_nrun(SEXP na, SEXP minrun) {
  const int n = static_cast<int>(clab::xlen(na, "na"));
  std::vector<int> st(static_cast<size_t>(n) + 1), en(static_cast<size_t>(n) + 1);
  const int k = clab_nrun(LOGICAL(na), n, Rf_asInteger(minrun), st.data(), en.data());
  SEXP a = PROTECT(Rf_allocVector(INTSXP, k));
  SEXP b = PROTECT(Rf_allocVector(INTSXP, k));
  for (int i = 0; i < k; ++i) {
    INTEGER(a)[i] = st[static_cast<size_t>(i)];
    INTEGER(b)[i] = en[static_cast<size_t>(i)];
  }
  SEXP r = list2(a, b);
  UNPROTECT(2);
  return r;
}

SEXP C_clab_r_msum(SEXP na) {
  SEXP o = PROTECT(Rf_allocVector(REALSXP, 8));
  clab_msum(LOGICAL(na), nr(na), nc(na), REAL(o));
  UNPROTECT(1);
  return o;
}

SEXP C_clab_r_pstat(SEXP x) {
  SEXP o = PROTECT(Rf_allocVector(REALSXP, 7));
  clab_pstat(REAL(x), static_cast<int>(clab::xlen(x, "x")), REAL(o));
  UNPROTECT(1);
  return o;
}

SEXP C_clab_r_q7(SEXP x, SEXP pr) {
  SEXP o = PROTECT(Rf_allocVector(REALSXP, clab::xlen(pr, "pr")));
  clab_q7(REAL(x), static_cast<int>(clab::xlen(x, "x")), REAL(pr), static_cast<int>(clab::xlen(pr, "pr")), REAL(o));
  UNPROTECT(1);
  return o;
}

// NULL for collinear columns.
SEXP C_clab_r_maha(SEXP x, SEXP robust, SEXP p) {
  clab::who = "mahaoutl";
  const int m = Rf_nrows(x), pp = Rf_asInteger(p);
  std::vector<double> d2(static_cast<size_t>(m));
  if (clab_maha(REAL(x), m, pp, Rf_asLogical(robust), d2.data()) != CLAB_OK) return R_NilValue;
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 3));
  for (int k = 0; k < 3; ++k) SET_VECTOR_ELT(r, k, Rf_allocVector(REALSXP, m));
  for (int i = 0; i < m; ++i) {
    const double d = d2[static_cast<size_t>(i)];
    REAL(VECTOR_ELT(r, 0))[i] = std::sqrt(d);
    REAL(VECTOR_ELT(r, 1))[i] = clab_pchisq(d, pp, 0);
    REAL(VECTOR_ELT(r, 2))[i] = clab_lpchisq(d, pp);
  }
  UNPROTECT(1);
  return r;
}

// NULL for collinear columns; an error when EM cannot invert a block.
SEXP C_clab_r_mcar(SEXP x, SEXP maxit, SEXP tol) {
  clab::who = "mcartst";
  const int n = Rf_nrows(x), p = Rf_ncols(x);
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(r, 0, Rf_allocVector(REALSXP, 8));
  SET_VECTOR_ELT(r, 1, Rf_allocVector(REALSXP, p));
  SET_VECTOR_ELT(r, 2, Rf_allocMatrix(REALSXP, p, p));
  double *o = REAL(VECTOR_ELT(r, 0));
  const int rc = clab_mcar(REAL(x), n, p, Rf_asInteger(maxit), Rf_asReal(tol), o, REAL(VECTOR_ELT(r, 1)),
                           REAL(VECTOR_ELT(r, 2)));
  clab::ichk(rc);
  if (rc != CLAB_OK) {
    if (o[0] == 1) {
      UNPROTECT(1);
      return R_NilValue;
    }
    clab::fail("the EM estimates did not converge to a usable covariance: the columns may be collinear, or a "
               "pattern may have too few rows.");
  }
  UNPROTECT(1);
  return r;
}

SEXP C_clab_r_mgrid(SEXP na, SEXP bins) {
  const int n = nr(na), p = nc(na), b = Rf_asInteger(bins);
  std::vector<int> st(static_cast<size_t>(b) + 1), sh((static_cast<size_t>(b) + 1) * static_cast<size_t>(p));
  const int k = clab_mgrid(LOGICAL(na), n, p, b, st.data(), sh.data());
  SEXP s = PROTECT(Rf_allocVector(INTSXP, k));
  SEXP m = PROTECT(Rf_allocMatrix(INTSXP, k, p));
  for (int i = 0; i < k; ++i) INTEGER(s)[i] = st[static_cast<size_t>(i)];
  for (R_xlen_t i = 0; i < static_cast<R_xlen_t>(k) * p; ++i) INTEGER(m)[i] = sh[static_cast<size_t>(i)];
  SEXP r = list2(s, m);
  UNPROTECT(2);
  return r;
}

SEXP C_clab_r_dups(SEXP code, SEXP nlev) {
  const int n = static_cast<int>(clab::xlen(code, "code"));
  std::vector<int> rw(static_cast<size_t>(n)), ct(static_cast<size_t>(n));
  const int k = clab_dups(INTEGER(code), n, Rf_asInteger(nlev), rw.data(), ct.data());
  SEXP a = PROTECT(Rf_allocVector(INTSXP, k));
  SEXP b = PROTECT(Rf_allocVector(INTSXP, k));
  for (int i = 0; i < k; ++i) {
    INTEGER(a)[i] = rw[static_cast<size_t>(i)];
    INTEGER(b)[i] = ct[static_cast<size_t>(i)];
  }
  SEXP r = list2(a, b);
  UNPROTECT(2);
  return r;
}

SEXP C_clab_r_cochq(SEXP y, SEXP s) {
  SEXP o = PROTECT(Rf_allocVector(REALSXP, 6));
  clab_cochq(REAL(y), REAL(s), static_cast<int>(clab::xlen(y, "y")), REAL(o));
  UNPROTECT(1);
  return o;
}

}  // extern "C"
