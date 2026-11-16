// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry point of toascii() (src/txt.cpp): text in, its ASCII
// folding out, NA as "". Nothing calls back into the R evaluator.
#include <vector>

#include "bcol.h"
#include "rchk.h"

extern "C" {

SEXP C_clab_r_ascii(SEXP x) {
  const R_xlen_t n = clab::xlen(x, "x");
  SEXP out = PROTECT(Rf_allocVector(STRSXP, n));
  std::vector<char> b;
  for (R_xlen_t i = 0; i < n; ++i) {
    SEXP e = STRING_ELT(x, i);
    if (e == NA_STRING) {
      SET_STRING_ELT(out, i, Rf_mkChar(""));
      continue;
    }
    const int len = LENGTH(e);
    b.resize(static_cast<std::size_t>(len) * 2 + 1);
    const int m = clab_ascii(CHAR(e), len, b.data());
    SET_STRING_ELT(out, i, Rf_mkCharLen(b.data(), m));
  }
  UNPROTECT(1);
  return out;
}

}  // extern "C"
