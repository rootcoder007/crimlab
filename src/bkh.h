// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Helpers of the R entry points of the geostatistics and spatial
// epidemiology families (bkrg.cpp, bsep.cpp): kernel codes, list
// elements, classed and named results, and the reading of a design
// matrix. Scratch comes from R_alloc(), freed by R on an error. Header
// only.
#ifndef CLAB_BKH_H
#define CLAB_BKH_H

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#include <cmath>
#include <cstring>

#include "crimlab.h"
#include "rchk.h"

namespace bkh {

using clab::fail;

inline SEXP real(R_xlen_t n) { return Rf_allocVector(REALSXP, n); }
inline SEXP ints(R_xlen_t n) { return Rf_allocVector(INTSXP, n); }
inline double *dal(R_xlen_t n) { return reinterpret_cast<double *>(R_alloc(n > 0 ? n : 1, sizeof(double))); }
inline int *ial(R_xlen_t n) { return reinterpret_cast<int *>(R_alloc(n > 0 ? n : 1, sizeof(int))); }

// a kernel's code: interrupted, refused with msg, or OK
inline void kc(int rc, const char *msg) {
  clab::ichk(rc);
  if (rc != CLAB_OK) fail("%s", msg);
}

inline SEXP elt(SEXP x, const char *nm) {
  SEXP n = Rf_getAttrib(x, R_NamesSymbol);
  for (R_xlen_t j = 0; j < Rf_xlength(n); ++j) {
    if (std::strcmp(CHAR(STRING_ELT(n, j)), nm) == 0) return VECTOR_ELT(x, j);
  }
  return R_NilValue;
}

inline SEXP klass(SEXP x, const char *cl) {
  Rf_classgets(x, PROTECT(Rf_mkString(cl)));
  UNPROTECT(1);
  return x;
}

// a named list of the k elements of v (protected by the caller)
inline SEXP lst(SEXP v, const char *const *nm, const char *cl) {
  clab::nlist(v, nm);
  return cl ? klass(v, cl) : v;
}

// The design argument z: NULL (an intercept alone when one), a numeric
// vector of n values (one column) or a numeric matrix of n rows; finite.
// Returned protected; *p its columns.
inline SEXP zmat(SEXP z, const char *nm, int n, bool one, int *p) {
  if (Rf_isNull(z)) {
    *p = one ? 1 : 0;
    SEXP o = PROTECT(real(one ? n : 0));
    for (int i = 0; one && i < n; ++i) REAL(o)[i] = 1;
    return o;
  }
  clab::unlab(z, nm);
  if (!clab::isnum(z)) fail("`%s` must be a numeric matrix or vector.", nm);
  if (Rf_isMatrix(z) ? Rf_nrows(z) != n : clab::xlen(z, nm) != n) fail("`%s` must have %d rows.", nm, n);
  *p = Rf_isMatrix(z) ? Rf_ncols(z) : 1;
  SEXP o = PROTECT(clab::dbl(z));
  for (R_xlen_t q = 0; q < XLENGTH(o); ++q) {
    if (!std::isfinite(REAL(o)[q])) fail("`%s` must hold finite numbers.", nm);
  }
  return o;
}





}  // namespace bkh

#endif
