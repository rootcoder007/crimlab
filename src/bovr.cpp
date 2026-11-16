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

extern "C" {

// ovrsoc(): the count, the rule and the reason of each text; tag, term,
// wtag and head as '|'-joined strings.
SEXP C_clab_r_ovrsoc(SEXP x, SEXP html, SEXP tag, SEXP term, SEXP wtag, SEXP head) {
  clab::who = "ovrsoc";
  const R_xlen_t n = XLENGTH(x);
  const int h = Rf_asLogical(html) == NA_LOGICAL ? -1 : Rf_asLogical(html);
  SEXP cnt = PROTECT(Rf_allocVector(INTSXP, n));
  SEXP rul = PROTECT(Rf_allocVector(INTSXP, n));
  SEXP why = PROTECT(Rf_allocVector(STRSXP, n));
  const char *a[4];
  SEXP in[4] = {tag, term, wtag, head};
  for (int k = 0; k < 4; ++k) a[k] = Rf_translateCharUTF8(STRING_ELT(in[k], 0));
  char w[1024];
  for (R_xlen_t i = 0; i < n; ++i) {
    INTEGER(cnt)[i] = NA_INTEGER;
    INTEGER(rul)[i] = NA_INTEGER;
    SET_STRING_ELT(why, i, NA_STRING);
    if (STRING_ELT(x, i) == NA_STRING) continue;
    int64_t m = 0, tl = 0;
    const char *s = u8(x, i, &m);
    char *t = R_alloc(static_cast<std::size_t>(m) + 1, 1);
    kc(clab_ovrtxt(s, m, h, t, &tl));
    int c = 0, r = 0;
    kc(clab_ovrsoc(t, tl, a[0], a[1], a[2], a[3], &c, &r, w, sizeof w));
    INTEGER(cnt)[i] = c < 0 ? NA_INTEGER : c;
    INTEGER(rul)[i] = r;
    SET_STRING_ELT(why, i, Rf_mkCharCE(w, CE_UTF8));
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP nm = PROTECT(Rf_allocVector(STRSXP, 3));
  static const char *const nms[] = {"count", "rule", "reason"};
  SEXP v[3] = {cnt, rul, why};
  for (int k = 0; k < 3; ++k) {
    SET_VECTOR_ELT(out, k, v[k]);
    SET_STRING_ELT(nm, k, Rf_mkChar(nms[k]));
  }
  Rf_setAttrib(out, R_NamesSymbol, nm);
  UNPROTECT(5);
  return out;
}

}  // extern "C"

extern "C" {

// The keyword groups of a preset: a list of group and term.
SEXP C_clab_r_ovrvoc(SEXP preset) {
  clab::who = "ovrkwd";
  const char *p = Rf_translateCharUTF8(STRING_ELT(preset, 0)), *g, *t;
  int m = 0;
  while (clab_ovrvoc(p, m, &g, &t) == CLAB_OK) ++m;
  SEXP gs = PROTECT(Rf_allocVector(STRSXP, m));
  SEXP ts = PROTECT(Rf_allocVector(STRSXP, m));
  for (int i = 0; i < m; ++i) {
    clab_ovrvoc(p, i, &g, &t);
    SET_STRING_ELT(gs, i, Rf_mkCharCE(g, CE_UTF8));
    SET_STRING_ELT(ts, i, Rf_mkCharCE(t, CE_UTF8));
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP nm = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_VECTOR_ELT(out, 0, gs);
  SET_VECTOR_ELT(out, 1, ts);
  SET_STRING_ELT(nm, 0, Rf_mkChar("group"));
  SET_STRING_ELT(nm, 1, Rf_mkChar("term"));
  Rf_setAttrib(out, R_NamesSymbol, nm);
  UNPROTECT(4);
  return out;
}

// ovrkwd(): the times each term occurs in each text, an integer matrix
// (texts by terms), NA for a text that is NA.
SEXP C_clab_r_ovrkwd(SEXP x, SEXP html, SEXP terms) {
  clab::who = "ovrkwd";
  const R_xlen_t n = XLENGTH(x), k = XLENGTH(terms);
  const int h = Rf_asLogical(html) == NA_LOGICAL ? -1 : Rf_asLogical(html);
  SEXP out = PROTECT(Rf_allocMatrix(INTSXP, static_cast<int>(n), static_cast<int>(k)));
  for (R_xlen_t i = 0; i < n; ++i) {
    const bool na = STRING_ELT(x, i) == NA_STRING;
    int64_t m = 0, tl = 0;
    char *t = nullptr;
    if (!na) {
      const char *s = u8(x, i, &m);
      t = R_alloc(static_cast<std::size_t>(m) + 1, 1);
      kc(clab_ovrtxt(s, m, h, t, &tl));
    }
    for (R_xlen_t j = 0; j < k; ++j) {
      int64_t c = 0;
      if (!na) kc(clab_ovrkwd(t, tl, Rf_translateCharUTF8(STRING_ELT(terms, j)), &c));
      INTEGER(out)[i + j * n] = na ? NA_INTEGER : static_cast<int>(c);
    }
  }
  UNPROTECT(1);
  return out;
}

}  // extern "C"

extern "C" {

// ovrchk(): for each value of x (text, NA where missing), 1 when it is of
// the kind (0 a date, 1 a count, 2 a flag, 3 of a shape), else 0; NA for
// a missing one.
SEXP C_clab_r_ovrchk(SEXP x, SEXP kind, SEXP shape) {
  clab::who = "ovrchk";
  const R_xlen_t n = XLENGTH(x);
  const int k = Rf_asInteger(kind);
  const char *sh = Rf_translateCharUTF8(STRING_ELT(shape, 0));
  SEXP out = PROTECT(Rf_allocVector(LGLSXP, n));
  for (R_xlen_t i = 0; i < n; ++i) {
    if (STRING_ELT(x, i) == NA_STRING) {
      LOGICAL(out)[i] = NA_LOGICAL;
      continue;
    }
    int64_t m = 0;
    int ok = 0;
    const char *s = u8(x, i, &m);
    kc(clab_ovrchk(s, m, k, sh, &ok));
    LOGICAL(out)[i] = ok;
  }
  UNPROTECT(1);
  return out;
}

}  // extern "C"
