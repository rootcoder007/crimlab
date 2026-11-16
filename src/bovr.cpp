// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the oversight report functions (src/ovr.cpp).
// The R wrappers (R/ovr.R) check the arguments and pass text as UTF-8,
// with NA for what is missing; these call the kernels one text at a time
// and return their results. Nothing calls back into the R evaluator.
#include <cstdint>
#include <cstring>

#include "crimlab.h"
#include "rchk.h"

namespace {

// The bytes of element i of the character vector x, as UTF-8.
const char *u8(SEXP x, R_xlen_t i, int64_t *n) {
  const char *s = Rf_translateCharUTF8(STRING_ELT(x, i));
  *n = static_cast<int64_t>(std::strlen(s));
  return s;
}

void kc(int rc) {
  clab::ichk(rc);
  if (rc != CLAB_OK) clab::fail("%s", clab_strerror(rc));
}

}  // namespace

extern "C" {

// ovrtxt(): each text as the rules read it; html -1 to detect markup in
// each, 0 or 1 to say.
SEXP C_clab_r_ovrtxt(SEXP x, SEXP html) {
  clab::who = "ovrtxt";
  const R_xlen_t n = XLENGTH(x);
  const int h = Rf_asLogical(html) == NA_LOGICAL ? -1 : Rf_asLogical(html);
  SEXP out = PROTECT(Rf_allocVector(STRSXP, n));
  for (R_xlen_t i = 0; i < n; ++i) {
    if (STRING_ELT(x, i) == NA_STRING) {
      SET_STRING_ELT(out, i, NA_STRING);
      continue;
    }
    int64_t m = 0, len = 0;
    const char *s = u8(x, i, &m);
    char *b = R_alloc(static_cast<std::size_t>(m) + 1, 1);
    kc(clab_ovrtxt(s, m, h, b, &len));
    SET_STRING_ELT(out, i, Rf_mkCharLenCE(b, static_cast<int>(len), CE_UTF8));
  }
  UNPROTECT(1);
  return out;
}

}  // extern "C"
