// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the byte checks of vrfpol()'s exposure file
// (src/vrfp.cpp). Raw vectors in, a flag or a code out.
#include <cstdint>

#include "bcol.h"
#include "rchk.h"

extern "C" {

SEXP C_clab_r_bz2end(SEXP t) {
  return Rf_ScalarLogical(clab_bz2end(RAW(t), XLENGTH(t)) == 1);
}

SEXP C_clab_r_gzend(SEXP d, SEXP t) {
  return Rf_ScalarLogical(clab_gzend(RAW(d), XLENGTH(d), RAW(t), XLENGTH(t)) == 1);
}

// The encoding (0 none found, 1 to 4 UTF-32LE, UTF-32BE, UTF-16LE,
// UTF-16BE, -1 NUL bytes that are not UTF-16) and the length of the
// byte-order mark.
SEXP C_clab_r_txtenc(SEXP r) {
  int bom = 0;
  const int e = clab_txtenc(RAW(r), XLENGTH(r), &bom);
  SEXP o = PROTECT(Rf_allocVector(INTSXP, 2));
  INTEGER(o)[0] = e < -1 ? -1 : e;
  INTEGER(o)[1] = bom;
  UNPROTECT(1);
  return o;
}

}  // extern "C"
