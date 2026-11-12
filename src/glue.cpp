// SPDX-License-Identifier: AGPL-3.0-or-later
//
// The R binding of libcrimlab: each entry point unpacks R vectors, calls one
// kernel and packs the result. The R wrappers have checked every argument, so
// a non-zero code here means the two layers disagree; it is reported in the
// kernel's own words. No C++ object with a destructor is alive when
// Rf_error() jumps.
#include <R.h>
#include <Rinternals.h>

#include "crimlab.h"

namespace {

void check(int code) {
  if (code != CLAB_OK) Rf_error("crimlab: %s", clab_strerror(code));
}

SEXP named_list(int n, const char *const *names) {
  SEXP out = PROTECT(Rf_allocVector(VECSXP, n));
  SEXP nm = PROTECT(Rf_allocVector(STRSXP, n));
  for (int i = 0; i < n; ++i) SET_STRING_ELT(nm, i, Rf_mkChar(names[i]));
  Rf_setAttrib(out, R_NamesSymbol, nm);
  UNPROTECT(2);
  return out;
}

}  // namespace

extern "C" {

SEXP C_clab_vod(SEXP day, SEXP night, SEXP ref, SEXP shift) {
  const int n = static_cast<int>(XLENGTH(day));
  const bool has = !Rf_isNull(shift);
  static const char *const names[] = {"share_day", "share_night", "ratio",
                                      "sign", "adj"};
  SEXP out = PROTECT(named_list(5, names));
  for (int i = 0; i < 3; ++i) SET_VECTOR_ELT(out, i, Rf_allocVector(REALSXP, n));
  SET_VECTOR_ELT(out, 3, Rf_allocVector(INTSXP, n));
  SET_VECTOR_ELT(out, 4, Rf_allocVector(REALSXP, has ? n : 0));
  check(clab_vod(REAL(day), REAL(night), n, Rf_asInteger(ref),
                 has ? REAL(shift) : nullptr, REAL(VECTOR_ELT(out, 0)),
                 REAL(VECTOR_ELT(out, 1)), REAL(VECTOR_ELT(out, 2)),
                 INTEGER(VECTOR_ELT(out, 3)), REAL(VECTOR_ELT(out, 4))));
  UNPROTECT(1);
  return out;
}

SEXP C_clab_orrisk(SEXP orr, SEXP q, SEXP s) {
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 2));
  check(clab_orrisk(Rf_asReal(orr), Rf_asReal(q), Rf_asReal(s), REAL(out),
                    REAL(out) + 1));
  UNPROTECT(1);
  return out;
}

SEXP C_clab_outtst(SEXP q, SEXP w, SEXP t) {
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 2));
  check(clab_outtst(REAL(q), REAL(w), static_cast<int>(XLENGTH(q)),
                    Rf_asReal(t), REAL(out), REAL(out) + 1));
  UNPROTECT(1);
  return out;
}

SEXP C_clab_dspbm(SEXP pop, SEXP contact, SEXP force, SEXP ref) {
  const int n = static_cast<int>(XLENGTH(pop));
  static const char *const names[] = {"cd", "fd", "rd"};
  SEXP out = PROTECT(named_list(3, names));
  for (int i = 0; i < 3; ++i) SET_VECTOR_ELT(out, i, Rf_allocVector(REALSXP, n));
  check(clab_dspbm(REAL(pop), REAL(contact), REAL(force), n, Rf_asInteger(ref),
                   REAL(VECTOR_ELT(out, 0)), REAL(VECTOR_ELT(out, 1)),
                   REAL(VECTOR_ELT(out, 2))));
  UNPROTECT(1);
  return out;
}

SEXP C_clab_dyadrt(SEXP count, SEXP off, SEXP vic, SEXP pop) {
  const int m = static_cast<int>(XLENGTH(count));
  static const char *const names[] = {"rate_off", "rate_pair", "null_off",
                                      "ratio"};
  SEXP out = PROTECT(named_list(4, names));
  for (int i = 0; i < 4; ++i) SET_VECTOR_ELT(out, i, Rf_allocVector(REALSXP, m));
  check(clab_dyadrt(REAL(count), INTEGER(off), INTEGER(vic), m, REAL(pop),
                    static_cast<int>(XLENGTH(pop)), REAL(VECTOR_ELT(out, 0)),
                    REAL(VECTOR_ELT(out, 1)), REAL(VECTOR_ELT(out, 2)),
                    REAL(VECTOR_ELT(out, 3))));
  UNPROTECT(1);
  return out;
}

SEXP C_clab_gtmp(SEXP p, SEXP q, SEXP hom, SEXP theta) {
  const int n = static_cast<int>(XLENGTH(p));
  static const char *const names[] = {"locus", "log10mp"};
  SEXP out = PROTECT(named_list(2, names));
  SET_VECTOR_ELT(out, 0, Rf_allocVector(REALSXP, n));
  SET_VECTOR_ELT(out, 1, Rf_allocVector(REALSXP, 1));
  check(clab_gtmp(REAL(p), REAL(q), LOGICAL(hom), n, Rf_asReal(theta),
                  REAL(VECTOR_ELT(out, 0)), REAL(VECTOR_ELT(out, 1))));
  UNPROTECT(1);
  return out;
}

SEXP C_clab_cllr(SEXP same, SEXP diff) {
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 1));
  check(clab_cllr(REAL(same), static_cast<int>(XLENGTH(same)), REAL(diff),
                  static_cast<int>(XLENGTH(diff)), REAL(out)));
  UNPROTECT(1);
  return out;
}

SEXP C_clab_lrcal(SEXP lr, SEXP same) {
  const int n = static_cast<int>(XLENGTH(lr));
  static const char *const names[] = {"post", "min"};
  SEXP out = PROTECT(named_list(2, names));
  SET_VECTOR_ELT(out, 0, Rf_allocVector(REALSXP, n));
  SET_VECTOR_ELT(out, 1, Rf_allocVector(REALSXP, 1));
  SEXP work = PROTECT(Rf_allocVector(REALSXP, 2 * static_cast<R_xlen_t>(n)));
  check(clab_lrcal(REAL(lr), LOGICAL(same), n, REAL(VECTOR_ELT(out, 0)),
                   REAL(work), REAL(VECTOR_ELT(out, 1))));
  UNPROTECT(2);
  return out;
}

SEXP C_clab_mixp(SEXP p, SEXP need, SEXP x) {
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 1));
  check(clab_mixp(REAL(p), LOGICAL(need), static_cast<int>(XLENGTH(p)),
                  Rf_asInteger(x), REAL(out)));
  UNPROTECT(1);
  return out;
}

SEXP C_clab_kinidx(SEXP g1, SEXP g2, SEXP f, SEXP k) {
  const int n = static_cast<int>(XLENGTH(g1) / 2);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  check(clab_kinidx(INTEGER(g1), INTEGER(g2), REAL(f), n, REAL(k), REAL(out)));
  UNPROTECT(1);
  return out;
}

SEXP C_clab_patidx(SEXP c, SEXP m, SEXP af, SEXP f) {
  const int n = static_cast<int>(XLENGTH(c) / 2);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  check(clab_patidx(INTEGER(c), INTEGER(m), INTEGER(af), REAL(f), n,
                    REAL(out)));
  UNPROTECT(1);
  return out;
}

SEXP C_clab_clrsmp(SEXP clra, SEXP reca, SEXP clrb, SEXP recb) {
  const int n = static_cast<int>(XLENGTH(clra));
  static const char *const names[] = {"ra", "rb", "wa", "wb", "reff",
                                      "meff", "tot", "cnt"};
  SEXP out = PROTECT(named_list(8, names));
  for (int i = 0; i < 6; ++i) SET_VECTOR_ELT(out, i, Rf_allocVector(REALSXP, n));
  SET_VECTOR_ELT(out, 6, Rf_allocVector(REALSXP, 5));
  SET_VECTOR_ELT(out, 7, Rf_allocVector(INTSXP, 3));
  check(clab_clrsmp(REAL(clra), REAL(reca), REAL(clrb), REAL(recb), n,
                    REAL(VECTOR_ELT(out, 0)), REAL(VECTOR_ELT(out, 1)),
                    REAL(VECTOR_ELT(out, 2)), REAL(VECTOR_ELT(out, 3)),
                    REAL(VECTOR_ELT(out, 4)), REAL(VECTOR_ELT(out, 5)),
                    REAL(VECTOR_ELT(out, 6)), INTEGER(VECTOR_ELT(out, 7))));
  UNPROTECT(1);
  return out;
}

SEXP C_clab_expbnd(SEXP y, SEXP m, SEXP ref, SEXP gamma) {
  const int n = static_cast<int>(XLENGTH(y));
  static const char *const names[] = {"rate", "lo", "hi", "ratio", "rlo",
                                      "rhi", "sign"};
  SEXP out = PROTECT(named_list(7, names));
  for (int i = 0; i < 6; ++i) SET_VECTOR_ELT(out, i, Rf_allocVector(REALSXP, n));
  SET_VECTOR_ELT(out, 6, Rf_allocVector(INTSXP, n));
  check(clab_expbnd(REAL(y), REAL(m), n, Rf_asInteger(ref), Rf_asReal(gamma),
                    REAL(VECTOR_ELT(out, 0)), REAL(VECTOR_ELT(out, 1)),
                    REAL(VECTOR_ELT(out, 2)), REAL(VECTOR_ELT(out, 3)),
                    REAL(VECTOR_ELT(out, 4)), REAL(VECTOR_ELT(out, 5)),
                    INTEGER(VECTOR_ELT(out, 6))));
  UNPROTECT(1);
  return out;
}

SEXP C_clab_misrd(SEXP r1, SEXP r0, SEXP p, SEXP se, SEXP sp) {
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 5));
  check(clab_misrd(Rf_asReal(r1), Rf_asReal(r0), Rf_asReal(p), Rf_asReal(se),
                   Rf_asReal(sp), REAL(out)));
  UNPROTECT(1);
  return out;
}

SEXP C_clab_rtmhs(SEXP x1, SEXP x2, SEXP w, SEXP c) {
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 7));
  check(clab_rtmhs(REAL(x1), REAL(x2), REAL(w), static_cast<int>(XLENGTH(x1)),
                   Rf_asReal(c), REAL(out)));
  UNPROTECT(1);
  return out;
}

SEXP C_clab_incap(SEXP lam, SEXP q, SEXP s, SEXP share) {
  const int n = static_cast<int>(XLENGTH(lam));
  const bool has = !Rf_isNull(share);
  static const char *const names[] = {"rate", "prev", "marg", "agg"};
  SEXP out = PROTECT(named_list(4, names));
  for (int i = 0; i < 3; ++i) SET_VECTOR_ELT(out, i, Rf_allocVector(REALSXP, n));
  SET_VECTOR_ELT(out, 3, Rf_allocVector(REALSXP, has ? 3 : 0));
  check(clab_incap(REAL(lam), REAL(q), REAL(s), n, has ? REAL(share) : nullptr,
                   REAL(VECTOR_ELT(out, 0)), REAL(VECTOR_ELT(out, 1)),
                   REAL(VECTOR_ELT(out, 2)), REAL(VECTOR_ELT(out, 3))));
  UNPROTECT(1);
  return out;
}

SEXP C_clab_incapr(SEXP lam, SEXP t0, SEXP s, SEXP r) {
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 5));
  check(clab_incapr(REAL(lam), static_cast<int>(XLENGTH(lam)), Rf_asInteger(t0),
                    Rf_asInteger(s), Rf_asReal(r), REAL(out)));
  UNPROTECT(1);
  return out;
}

SEXP C_clab_sellab(SEXP y, SEXP rel, SEXP rule, SEXP w) {
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 6));
  check(clab_sellab(REAL(y), LOGICAL(rel), LOGICAL(rule), REAL(w),
                    static_cast<int>(XLENGTH(y)), REAL(out)));
  UNPROTECT(1);
  return out;
}

SEXP C_clab_losbnd(SEXP t, SEXP a) {
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 4));
  check(clab_losbnd(REAL(t), static_cast<int>(XLENGTH(t)), REAL(a),
                    static_cast<int>(XLENGTH(a)), REAL(out)));
  UNPROTECT(1);
  return out;
}

SEXP C_clab_sntbnd(SEXP y, SEXP z, SEXP w) {
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 11));
  check(clab_sntbnd(REAL(y), INTEGER(z), REAL(w), static_cast<int>(XLENGTH(y)),
                    REAL(out)));
  UNPROTECT(1);
  return out;
}

SEXP C_clab_hzdsel(SEXP s, SEXP h, SEXP l, SEXP ah, SEXP al) {
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 5));
  check(clab_hzdsel(Rf_asReal(s), Rf_asReal(h), Rf_asReal(l), REAL(ah),
                    REAL(al), REAL(out)));
  UNPROTECT(1);
  return out;
}

// The self-test (selftest.cpp): the number of checks run, the number that
// failed and their report, one line each.
int clab_selftest(char *sink, int cap, int *nfail);

SEXP C_clab_selftest(void) {
  static char sink[1 << 16];
  int nf = 0;
  const int n = clab_selftest(sink, static_cast<int>(sizeof sink), &nf);
  static const char *const names[] = {"n", "nfail", "report"};
  SEXP out = PROTECT(named_list(3, names));
  SET_VECTOR_ELT(out, 0, Rf_ScalarInteger(n));
  SET_VECTOR_ELT(out, 1, Rf_ScalarInteger(nf));
  SET_VECTOR_ELT(out, 2, Rf_mkString(sink));
  UNPROTECT(1);
  return out;
}

static const R_CallMethodDef calls[] = {
    {"C_clab_vod", (DL_FUNC)&C_clab_vod, 4},
    {"C_clab_orrisk", (DL_FUNC)&C_clab_orrisk, 3},
    {"C_clab_outtst", (DL_FUNC)&C_clab_outtst, 3},
    {"C_clab_dspbm", (DL_FUNC)&C_clab_dspbm, 4},
    {"C_clab_dyadrt", (DL_FUNC)&C_clab_dyadrt, 4},
    {"C_clab_gtmp", (DL_FUNC)&C_clab_gtmp, 4},
    {"C_clab_cllr", (DL_FUNC)&C_clab_cllr, 2},
    {"C_clab_lrcal", (DL_FUNC)&C_clab_lrcal, 2},
    {"C_clab_mixp", (DL_FUNC)&C_clab_mixp, 3},
    {"C_clab_kinidx", (DL_FUNC)&C_clab_kinidx, 4},
    {"C_clab_patidx", (DL_FUNC)&C_clab_patidx, 4},
    {"C_clab_clrsmp", (DL_FUNC)&C_clab_clrsmp, 4},
    {"C_clab_expbnd", (DL_FUNC)&C_clab_expbnd, 4},
    {"C_clab_misrd", (DL_FUNC)&C_clab_misrd, 5},
    {"C_clab_rtmhs", (DL_FUNC)&C_clab_rtmhs, 4},
    {"C_clab_incap", (DL_FUNC)&C_clab_incap, 4},
    {"C_clab_incapr", (DL_FUNC)&C_clab_incapr, 4},
    {"C_clab_sellab", (DL_FUNC)&C_clab_sellab, 4},
    {"C_clab_losbnd", (DL_FUNC)&C_clab_losbnd, 2},
    {"C_clab_sntbnd", (DL_FUNC)&C_clab_sntbnd, 3},
    {"C_clab_hzdsel", (DL_FUNC)&C_clab_hzdsel, 5},
    {"C_clab_selftest", (DL_FUNC)&C_clab_selftest, 0},
    {nullptr, nullptr, 0}};

void R_init_crimlab(DllInfo *dll) {
  R_registerRoutines(dll, nullptr, calls, nullptr, nullptr);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}

}  // extern "C"
