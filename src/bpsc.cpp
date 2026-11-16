// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of pscred() and synfix() (src/psc.cpp,
// src/sfx.cpp). The R wrappers (R/psc.R, R/sfx.R) check the arguments
// and pass doubles, NULL for what is left out; these call the kernels and
// return their arrays. Nothing calls back into the R evaluator.
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

// The presets of synfix(): a list of kind, column, type, lo, hi and lev,
// one value per column of every preset.
SEXP C_clab_r_sfxtab(void) {
  clab::who = "synfix";
  int m = 0;
  const char *k, *nm, *lv;
  int t, lo, hi;
  while (clab_sfxcol(m, &k, &nm, &t, &lo, &hi, &lv) == CLAB_OK) ++m;
  SEXP v[6];
  for (int j = 0; j < 6; ++j) v[j] = PROTECT(Rf_allocVector(j == 2 || j == 3 || j == 4 ? INTSXP : STRSXP, m));
  for (int i = 0; i < m; ++i) {
    clab_sfxcol(i, &k, &nm, &t, &lo, &hi, &lv);
    SET_STRING_ELT(v[0], i, Rf_mkCharCE(k, CE_UTF8));
    SET_STRING_ELT(v[1], i, Rf_mkCharCE(nm, CE_UTF8));
    INTEGER(v[2])[i] = t;
    INTEGER(v[3])[i] = lo;
    INTEGER(v[4])[i] = hi;
    SET_STRING_ELT(v[5], i, Rf_mkCharCE(lv, CE_UTF8));
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 6));
  SEXP s = PROTECT(Rf_allocVector(STRSXP, 6));
  static const char *const nms[] = {"kind", "column", "type", "lo", "hi", "lev"};
  for (int j = 0; j < 6; ++j) {
    SET_VECTOR_ELT(out, j, v[j]);
    SET_STRING_ELT(s, j, Rf_mkChar(nms[j]));
  }
  Rf_setAttrib(out, R_NamesSymbol, s);
  UNPROTECT(8);
  return out;
}

// n draws of each of the length(dn) columns: an n x k integer matrix of
// indices from 0.
SEXP C_clab_r_synfix(SEXP dn, SEXP n, SEXP seed) {
  clab::who = "synfix";
  const int k = static_cast<int>(XLENGTH(dn));
  const R_xlen_t nn = static_cast<R_xlen_t>(Rf_asReal(n));
  SEXP out = PROTECT(clab::xalloc(INTSXP, nn * k, "n"));
  const int rc = clab_synfix(REAL(dn), k, nn, static_cast<uint32_t>(Rf_asReal(seed)), INTEGER(out));
  clab::ichk(rc);
  if (rc != CLAB_OK) clab::fail("%s", clab_strerror(rc));
  SEXP d = PROTECT(Rf_allocVector(INTSXP, 2));
  INTEGER(d)[0] = static_cast<int>(nn);
  INTEGER(d)[1] = k;
  Rf_setAttrib(out, R_DimSymbol, d);
  UNPROTECT(2);
  return out;
}

}  // extern "C"
