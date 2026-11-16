// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of data description (src/dsc.cpp). The R wrappers
// (R/dsc.R) check the arguments and take the data frames apart. Nothing
// calls back into the R evaluator.
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

SEXP C_clab_r_msum(SEXP na) {
  SEXP o = PROTECT(Rf_allocVector(REALSXP, 8));
  clab_msum(LOGICAL(na), nr(na), nc(na), REAL(o));
  UNPROTECT(1);
  return o;
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
