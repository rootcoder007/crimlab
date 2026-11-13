// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Helpers of the R entry points of the profiling, toxicology, taphonomy
// and deterrence families (brtm.cpp, btox.cpp, btph.cpp, bdtr.cpp): the
// reading of matrices and the building of matrices. Header only.
#ifndef CLAB_BFX_H
#define CLAB_BFX_H

#include <cmath>

#include "bcol.h"

namespace bfx {

// A numeric matrix of finite numbers from lo with nc columns (any when nc
// < 0) and at least one row; a vector is one column. Its number of rows in
// nr. Not protected.
inline SEXP mat(SEXP x, const char *nm, int nc, int *nr, double lo = clab::ninf) {
  clab::unlab(x, nm);
  SEXP d = Rf_getAttrib(x, R_DimSymbol);
  const bool m = !Rf_isNull(d) && XLENGTH(d) == 2;
  const R_xlen_t r = m ? INTEGER(d)[0] : XLENGTH(x), c = m ? INTEGER(d)[1] : 1;
  if (!clab::isnum(x) || r < 1 || (nc >= 0 && c != nc)) {
    if (nc >= 0) clab::fail("`%s` must be a numeric matrix of %d columns.", nm, nc);
    clab::fail("`%s` must be a numeric matrix.", nm);
  }
  SEXP v = clab::dbl(x);
  for (R_xlen_t i = 0; i < XLENGTH(v); ++i) {
    if (!std::isfinite(REAL(v)[i]) || REAL(v)[i] < lo) {
      clab::fail("`%s` must hold finite numbers from %g.", nm, lo);
    }
  }
  *nr = clab::ilen(static_cast<double>(r), nm);
  return v;
}

// Numbers or NA (no infinite values), n of them or one recycled to n
// (n < 0: any length, at least one). A fresh vector, not protected.
inline SEXP nna(SEXP x, const char *nm, R_xlen_t n) {
  clab::unlab(x, nm);
  const R_xlen_t m = XLENGTH(x);
  if (!(clab::isnum(x) || TYPEOF(x) == LGLSXP) || m < 1 || (n >= 0 && m != n && m != 1)) {
    clab::fail("`%s` must be numbers or NA, one per case or one for all.", nm);
  }
  clab::vec(x, nm);
  const R_xlen_t k = n < 0 ? m : n;
  SEXP v = PROTECT(Rf_allocVector(REALSXP, k));
  SEXP d = PROTECT(Rf_coerceVector(x, REALSXP));
  for (R_xlen_t i = 0; i < k; ++i) {
    const double e = REAL(d)[m == 1 ? 0 : i];
    if (std::isinf(e)) clab::fail("`%s` must be numbers or NA, one per case or one for all.", nm);
    REAL(v)[i] = e;
  }
  UNPROTECT(2);
  return v;
}

// An r by c double matrix (protected by the caller) holding v.
inline SEXP dmat(int r, int c, const double *v) {
  SEXP m = Rf_allocMatrix(REALSXP, r, c);
  for (R_xlen_t i = 0; i < static_cast<R_xlen_t>(r) * c; ++i) REAL(m)[i] = v[i];
  return m;
}

// dimnames list(rn, cn) on m (either may be R_NilValue).
inline void dnm(SEXP m, SEXP rn, SEXP cn) {
  SEXP dn = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dn, 0, rn);
  SET_VECTOR_ELT(dn, 1, cn);
  Rf_setAttrib(m, R_DimNamesSymbol, dn);
  UNPROTECT(1);
}

// A character vector of the n C strings s.
inline SEXP strs(const char *const *s, int n) {
  SEXP v = PROTECT(Rf_allocVector(STRSXP, n));
  for (int i = 0; i < n; ++i) SET_STRING_ELT(v, i, Rf_mkCharCE(s[i], CE_UTF8));
  UNPROTECT(1);
  return v;
}

}  // namespace bfx

#endif
