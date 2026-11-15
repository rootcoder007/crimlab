// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry point of the writers' numbers (src/ynum.h): mode 0 every
// digit, 1 fixed with commas, 2 a signed percentage, 3 signed percentage
// points; NA stays NA. Nothing calls back into the R evaluator.
#include "bcol.h"
#include "rchk.h"
#include "ynum.h"

extern "C" {

SEXP C_clab_r_fmt(SEXP x, SEXP mode, SEXP digits) {
  const R_xlen_t n = XLENGTH(x);
  const int m = Rf_asInteger(mode), d = Rf_asInteger(digits);
  SEXP out = PROTECT(Rf_allocVector(STRSXP, n));
  for (R_xlen_t i = 0; i < n; ++i) {
    const double v = REAL(x)[i];
    if (ISNA(v) || (m > 0 && std::isnan(v))) {
      SET_STRING_ELT(out, i, NA_STRING);
      continue;
    }
    const std::string s = m == 0 ? ynum::all(v) : m == 1 ? ynum::fix(v, d) : ynum::pct(v, d, m == 2 ? "%" : "pp");
    SET_STRING_ELT(out, i, Rf_mkCharLen(s.data(), static_cast<int>(s.size())));
  }
  UNPROTECT(1);
  return out;
}

}  // extern "C"
