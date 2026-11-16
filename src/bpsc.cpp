// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry point of pscred() (src/psc.cpp). The R wrapper (R/psc.R)
// checks the arguments and passes doubles, NULL for a sentence left out;
// this calls the kernel and returns its arrays. Nothing calls back into
// the R evaluator.
#include "crimlab.h"
#include "rchk.h"

extern "C" {

SEXP C_clab_r_pscred(SEXP days, SEXP ratio, SEXP sent) {
  clab::who = "pscred";
  const R_xlen_t n = XLENGTH(days);
  const bool hs = !Rf_isNull(sent);
  SEXP cr = PROTECT(Rf_allocVector(REALSXP, n));
  SEXP rm = PROTECT(Rf_allocVector(REALSXP, hs ? n : 0));
  const int rc = clab_pscred(REAL(days), n, REAL(ratio), XLENGTH(ratio), hs ? REAL(sent) : nullptr,
                             hs ? XLENGTH(sent) : 0, REAL(cr), hs ? REAL(rm) : nullptr);
  clab::ichk(rc);
  if (rc != CLAB_OK) clab::fail("%s", clab_strerror(rc));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, cr);
  SET_VECTOR_ELT(out, 1, rm);
  SEXP nm = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(nm, 0, Rf_mkChar("credit"));
  SET_STRING_ELT(nm, 1, Rf_mkChar("remain"));
  Rf_setAttrib(out, R_NamesSymbol, nm);
  UNPROTECT(4);
  return out;
}

}  // extern "C"
