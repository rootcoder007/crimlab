// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the one-pass summaries (src/skch.cpp). The R
// wrappers (R/skch.R) check the arguments. Nothing calls back into the R
// evaluator.
#include <cstring>
#include <vector>

#include "bcol.h"
#include "rchk.h"

extern "C" {

// the CRC-32 of the bytes of the raw vector x from position from (zero
// based), as a double
SEXP C_clab_r_crc32(SEXP x, SEXP from) {
  const R_xlen_t n = clab::xlen(x, "x"), f = static_cast<R_xlen_t>(Rf_asReal(from));
  const R_xlen_t k = f < 0 || f > n ? n : f;
  return Rf_ScalarReal(clab_crc32(RAW(x) + k, static_cast<std::size_t>(n - k), 0));
}

SEXP C_clab_r_momupd(SEXP st, SEXP x) {
  clab::who = "smryupd";
  double b[5];
  clab_momacc(REAL(x), static_cast<double>(clab::xlen(x, "x")), b);
  SEXP o = PROTECT(Rf_allocVector(REALSXP, 5));
  clab_mommrg(REAL(st), b, REAL(o));
  UNPROTECT(1);
  return o;
}

SEXP C_clab_r_mommrg(SEXP a, SEXP b) {
  SEXP o = PROTECT(Rf_allocVector(REALSXP, 5));
  clab_mommrg(REAL(a), REAL(b), REAL(o));
  UNPROTECT(1);
  return o;
}

SEXP C_clab_r_momstat(SEXP s) {
  SEXP o = PROTECT(Rf_allocVector(REALSXP, 6));
  clab_momstat(REAL(s), REAL(o));
  UNPROTECT(1);
  return o;
}

SEXP C_clab_r_rsv(SEXP n, SEXP k, SEXP seed) {
  clab::who = "rsvidx";
  const double nd = Rf_asReal(n), kd = Rf_asReal(k);
  SEXP o = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(kd < nd ? kd : nd)));
  const int rc = clab_rsv(nd, kd, Rf_asReal(seed), REAL(o));
  clab::ichk(rc);
  if (rc != CLAB_OK) clab::fail("%s", clab_strerror(rc));
  UNPROTECT(1);
  return o;
}

SEXP C_clab_r_hlladd(SEXP x, SEXP p, SEXP reg) {
  clab::who = "dstnskch";
  const int pp = Rf_asInteger(p);
  const R_xlen_t m = static_cast<R_xlen_t>(1) << pp;
  SEXP o = PROTECT(Rf_allocVector(INTSXP, m));
  for (R_xlen_t i = 0; i < m; ++i) INTEGER(o)[i] = Rf_isNull(reg) ? 0 : INTEGER(reg)[i];
  for (R_xlen_t i = 0; i < clab::xlen(x, "x"); ++i) {
    SEXP e = STRING_ELT(x, i);
    if (e == NA_STRING) continue;
    const char *s = CHAR(e);
    clab_hlladd(s, std::strlen(s), pp, INTEGER(o));
  }
  UNPROTECT(1);
  return o;
}

SEXP C_clab_r_hllcnt(SEXP reg) {
  return Rf_ScalarReal(clab_hllcnt(INTEGER(reg), static_cast<std::size_t>(clab::xlen(reg, "reg"))));
}

}  // extern "C"
