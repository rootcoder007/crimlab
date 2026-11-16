// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the footprint family (src/fp.cpp). The R wrappers
// (R/fp.R) check the arguments and read the counters; these look up the
// tables and compute the energy. Nothing calls back into the R evaluator.
#include <vector>

#include "bcol.h"
#include "rchk.h"

extern "C" {

SEXP C_clab_r_citab(void) {
  const int n = clab_cin();
  SEXP v = PROTECT(Rf_allocVector(VECSXP, 7));
  SEXP s[5];
  for (int j : {0, 1, 2, 5, 6}) SET_VECTOR_ELT(v, j, Rf_allocVector(STRSXP, n));
  SET_VECTOR_ELT(v, 3, Rf_allocVector(REALSXP, n));
  SET_VECTOR_ELT(v, 4, Rf_allocVector(INTSXP, n));
  s[0] = VECTOR_ELT(v, 0);
  s[1] = VECTOR_ELT(v, 1);
  s[2] = VECTOR_ELT(v, 2);
  s[3] = VECTOR_ELT(v, 5);
  s[4] = VECTOR_ELT(v, 6);
  for (int i = 0; i < n; ++i) {
    const char *loc, *iso3, *name, *src;
    double g;
    int yr;
    clab_cirow(i, &loc, &iso3, &name, &g, &yr, &src);
    SET_STRING_ELT(s[0], i, Rf_mkCharCE(loc, CE_UTF8));
    SET_STRING_ELT(s[1], i, Rf_mkCharCE(iso3, CE_UTF8));
    SET_STRING_ELT(s[2], i, Rf_mkCharCE(name, CE_UTF8));
    SET_STRING_ELT(s[3], i, Rf_mkChar("lifecycle"));
    SET_STRING_ELT(s[4], i, Rf_mkCharCE(src, CE_UTF8));
    REAL(VECTOR_ELT(v, 3))[i] = g;
    INTEGER(VECTOR_ELT(v, 4))[i] = yr;
  }
  static const char *const nm[] = {"location", "iso3", "name", "g_per_kwh", "year", "basis", "source"};
  SEXP f = clab::frame(v, nm);
  UNPROTECT(1);
  return f;
}

SEXP C_clab_r_proxy(SEXP terr) { return Rf_mkString(clab_terrproxy(CHAR(STRING_ELT(terr, 0)))); }

SEXP C_clab_r_fpe(SEXP in, SEXP r0, SEXP r1, SEXP rmx, SEXP dw, SEXP de, SEXP dk) {
  clab::who = "co2fp";
  const int nr = static_cast<int>(clab::xlen(r0, "r0")), nd = static_cast<int>(clab::xlen(dw, "dw"));
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(r, 0, Rf_allocVector(REALSXP, 10));
  SET_VECTOR_ELT(r, 1, Rf_allocVector(REALSXP, nd));
  const int rc = clab_fpe(REAL(in), nr, REAL(r0), REAL(r1), REAL(rmx), nd, REAL(dw), REAL(de), INTEGER(dk),
                          REAL(VECTOR_ELT(r, 0)), REAL(VECTOR_ELT(r, 1)));
  if (rc != CLAB_OK) clab::fail("%s", clab_strerror(rc));
  UNPROTECT(1);
  return r;
}

SEXP C_clab_r_fpeqv(SEXP g) {
  SEXP o = PROTECT(Rf_allocVector(REALSXP, 2));
  clab_fpeqv(Rf_asReal(g), REAL(o));
  UNPROTECT(1);
  return o;
}

}  // extern "C"
