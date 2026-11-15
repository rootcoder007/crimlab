// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the exposure-and-crime family (src/apx.cpp and
// clab_ppmlfe, clab_fesep of src/cnt.cpp). The R wrappers (R/apx.R) check
// the arguments and pass codes from 0, numeric matrices (with no columns
// for none) and NULL for what is left out; these call the kernels and
// build the results. Scratch comes from R_alloc, so an error leaks
// nothing. Nothing calls back into the R evaluator.
#include <cmath>

#include "bcol.h"
#include "rchk.h"

namespace {

int n32(SEXP x) { return clab::ilen(static_cast<double>(XLENGTH(x)), nullptr); }

const char *const absorbed =
    "the fit failed: a regressor or an instrument does not vary within the fixed effects (they "
    "absorb it), or the columns are collinear";

void kc(int rc) {
  clab::ichk(rc);
  if (rc == CLAB_EDOM) clab::fail("%s.", absorbed);
  if (rc != CLAB_OK) clab::fail("%s", clab_strerror(rc));
}

const int *icodes(SEXP x) { return Rf_isNull(x) ? nullptr : INTEGER(x); }
const double *dnul(SEXP x) { return Rf_isNull(x) ? nullptr : REAL(x); }

SEXP named(SEXP v, const char *const *nm, int k) {
  SEXP s = PROTECT(Rf_allocVector(STRSXP, k));
  for (int j = 0; j < k; ++j) SET_STRING_ELT(s, j, Rf_mkChar(nm[j]));
  Rf_setAttrib(v, R_NamesSymbol, s);
  UNPROTECT(1);
  return v;
}

}  // namespace

extern "C" {

SEXP C_clab_r_fesep(SEXP y, SEXP f1, SEXP n1, SEXP f2, SEXP n2) {
  clab::who = "ppmlfe";
  const int n = n32(y), a = Rf_asInteger(n1), b = Rf_isNull(f2) ? 0 : Rf_asInteger(n2);
  SEXP k = PROTECT(Rf_allocVector(LGLSXP, n));
  double *w = reinterpret_cast<double *>(R_alloc(static_cast<size_t>(a + b + 1), sizeof(double)));
  if (clab_fesep(REAL(y), n, INTEGER(f1), a, icodes(f2), b, LOGICAL(k), w) < 0) {
    clab::fail("%s", clab_strerror(CLAB_EARG));
  }
  UNPROTECT(1);
  return k;
}

SEXP C_clab_r_ppmlfe(SEXP y, SEXP x, SEXP off, SEXP f1, SEXP n1, SEXP f2, SEXP n2, SEXP cl,
                     SEXP ncl, SEXP cl2, SEXP ncl2, SEXP tol, SEXP maxit) {
  clab::who = "ppmlfe";
  const int n = n32(y), p = Rf_ncols(x), a = Rf_asInteger(n1);
  const int b = Rf_isNull(f2) ? 0 : Rf_asInteger(n2);
  const int g1 = Rf_isNull(cl) ? 0 : Rf_asInteger(ncl), g2 = Rf_isNull(cl2) ? 0 : Rf_asInteger(ncl2);
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(r, 0, Rf_allocVector(REALSXP, p));
  SET_VECTOR_ELT(r, 1, Rf_allocMatrix(REALSXP, p, p));
  SET_VECTOR_ELT(r, 2, Rf_allocVector(REALSXP, 6));
  const size_t sn = static_cast<size_t>(n), sp = static_cast<size_t>(p);
  double *w = reinterpret_cast<double *>(
      R_alloc(sn * (2 * sp + 8) + 3 * sp * sp + sp + 4 * static_cast<size_t>(a + b), sizeof(double)));
  int *iw = reinterpret_cast<int *>(
      R_alloc(sn + static_cast<size_t>(g1) * static_cast<size_t>(g2) + 1, sizeof(int)));
  kc(clab_ppmlfe(REAL(y), REAL(x), n, p, dnul(off), nullptr, INTEGER(f1), a, icodes(f2), b,
                 icodes(cl), g1, icodes(cl2), g2, Rf_asInteger(maxit), 1e-10, Rf_asReal(tol),
                 REAL(VECTOR_ELT(r, 0)), REAL(VECTOR_ELT(r, 1)), REAL(VECTOR_ELT(r, 2)), w, iw));
  static const char *const nm[] = {"beta", "vc", "st"};
  named(r, nm, 3);
  UNPROTECT(1);
  return r;
}

}  // extern "C"
