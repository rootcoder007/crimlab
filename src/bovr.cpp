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

extern "C" {

// ovrdte(): the first date in each text as YYYY-MM-DD, NA for none.
SEXP C_clab_r_ovrdte(SEXP x) {
  clab::who = "ovrdte";
  const R_xlen_t n = XLENGTH(x);
  SEXP out = PROTECT(Rf_allocVector(STRSXP, n));
  char b[16];
  for (R_xlen_t i = 0; i < n; ++i) {
    SET_STRING_ELT(out, i, NA_STRING);
    if (STRING_ELT(x, i) == NA_STRING) continue;
    int64_t m = 0;
    const char *s = u8(x, i, &m);
    kc(clab_ovrdte(s, m, b));
    if (b[0]) SET_STRING_ELT(out, i, Rf_mkChar(b));
  }
  UNPROTECT(1);
  return out;
}

}  // extern "C"

extern "C" {

// ovrsch(): the fields of a preset, a list of six character columns.
SEXP C_clab_r_ovrsch(SEXP preset) {
  clab::who = "ovrsch";
  const char *p = Rf_translateCharUTF8(STRING_ELT(preset, 0));
  const char *v[6];
  int m = 0;
  while (clab_ovrsch(p, m, &v[0], &v[1], &v[2], &v[3], &v[4], &v[5]) == CLAB_OK) ++m;
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 6));
  SEXP nm = PROTECT(Rf_allocVector(STRSXP, 6));
  static const char *const nms[] = {"name", "rule", "cue", "alt", "sect", "desc"};
  for (int k = 0; k < 6; ++k) {
    SET_VECTOR_ELT(out, k, Rf_allocVector(STRSXP, m));
    SET_STRING_ELT(nm, k, Rf_mkChar(nms[k]));
  }
  for (int i = 0; i < m; ++i) {
    clab_ovrsch(p, i, &v[0], &v[1], &v[2], &v[3], &v[4], &v[5]);
    for (int k = 0; k < 6; ++k) SET_STRING_ELT(VECTOR_ELT(out, k), i, Rf_mkCharCE(v[k], CE_UTF8));
  }
  Rf_setAttrib(out, R_NamesSymbol, nm);
  UNPROTECT(2);
  return out;
}

}  // extern "C"

extern "C" {

// ovrpar(): the fields of each text, a list of character columns (NA
// where a field is not found or the text is NA); rule the codes of the
// rules, cue, alt and sect their strings.
SEXP C_clab_r_ovrpar(SEXP x, SEXP html, SEXP rule, SEXP cue, SEXP alt, SEXP sect) {
  clab::who = "ovrpar";
  const R_xlen_t n = XLENGTH(x), k = XLENGTH(rule);
  const int h = Rf_asLogical(html) == NA_LOGICAL ? -1 : Rf_asLogical(html);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, k));
  for (R_xlen_t j = 0; j < k; ++j) SET_VECTOR_ELT(out, j, Rf_allocVector(STRSXP, n));
  for (R_xlen_t i = 0; i < n; ++i) {
    for (R_xlen_t j = 0; j < k; ++j) SET_STRING_ELT(VECTOR_ELT(out, j), i, NA_STRING);
    if (STRING_ELT(x, i) == NA_STRING) continue;
    int64_t m = 0, tl = 0, len = 0;
    const char *s = u8(x, i, &m);
    char *t = R_alloc(static_cast<std::size_t>(m) + 1, 1);
    kc(clab_ovrtxt(s, m, h, t, &tl));
    char *b = R_alloc(2 * static_cast<std::size_t>(tl) + 64, 1);
    for (R_xlen_t j = 0; j < k; ++j) {
      kc(clab_ovrfld(t, tl, INTEGER(rule)[j], Rf_translateCharUTF8(STRING_ELT(cue, j)),
                     Rf_translateCharUTF8(STRING_ELT(alt, j)), Rf_translateCharUTF8(STRING_ELT(sect, j)), b, &len));
      if (len >= 0) SET_STRING_ELT(VECTOR_ELT(out, j), i, Rf_mkCharLenCE(b, static_cast<int>(len), CE_UTF8));
    }
  }
  UNPROTECT(1);
  return out;
}

}  // extern "C"
