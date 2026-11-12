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
#include "rchk.h"

namespace {

// An entry point as R's registration table holds it. The cast goes through
// void (*)(void), which compilers accept as matching every function type,
// so -Wcast-function-type (in -Wextra) stays quiet on (DL_FUNC) &f.
template <class F>
DL_FUNC fn(F f) {
  return reinterpret_cast<DL_FUNC>(reinterpret_cast<void (*)(void)>(f));
}

// The length of x for a kernel, which counts in int: a vector longer than
// INT_MAX is refused with the kernels' out-of-range error, not truncated.
int len(SEXP x) { return clab::ilen(static_cast<double>(XLENGTH(x)), nullptr); }

void check(int code) {
  if (code != CLAB_OK) Rf_error("%s: %s", clab::who, clab_strerror(code));
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
  clab::who = "clab_vod";
  const int n = len(day);
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
  clab::who = "clab_orrisk";
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 2));
  check(clab_orrisk(Rf_asReal(orr), Rf_asReal(q), Rf_asReal(s), REAL(out),
                    REAL(out) + 1));
  UNPROTECT(1);
  return out;
}

SEXP C_clab_outtst(SEXP q, SEXP w, SEXP t) {
  clab::who = "clab_outtst";
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 2));
  check(clab_outtst(REAL(q), REAL(w), len(q),
                    Rf_asReal(t), REAL(out), REAL(out) + 1));
  UNPROTECT(1);
  return out;
}

SEXP C_clab_dspbm(SEXP pop, SEXP contact, SEXP force, SEXP ref) {
  clab::who = "clab_dspbm";
  const int n = len(pop);
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
  clab::who = "clab_dyadrt";
  const int m = len(count);
  static const char *const names[] = {"rate_off", "rate_pair", "null_off",
                                      "ratio"};
  SEXP out = PROTECT(named_list(4, names));
  for (int i = 0; i < 4; ++i) SET_VECTOR_ELT(out, i, Rf_allocVector(REALSXP, m));
  check(clab_dyadrt(REAL(count), INTEGER(off), INTEGER(vic), m, REAL(pop),
                    len(pop), REAL(VECTOR_ELT(out, 0)),
                    REAL(VECTOR_ELT(out, 1)), REAL(VECTOR_ELT(out, 2)),
                    REAL(VECTOR_ELT(out, 3))));
  UNPROTECT(1);
  return out;
}

SEXP C_clab_gtmp(SEXP p, SEXP q, SEXP hom, SEXP theta) {
  clab::who = "clab_gtmp";
  const int n = len(p);
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
  clab::who = "clab_cllr";
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 1));
  check(clab_cllr(REAL(same), len(same), REAL(diff),
                  len(diff), REAL(out)));
  UNPROTECT(1);
  return out;
}

SEXP C_clab_lrcal(SEXP lr, SEXP same) {
  clab::who = "clab_lrcal";
  const int n = len(lr);
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
  clab::who = "clab_mixp";
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 1));
  check(clab_mixp(REAL(p), LOGICAL(need), len(p),
                  Rf_asInteger(x), REAL(out)));
  UNPROTECT(1);
  return out;
}

SEXP C_clab_kinidx(SEXP g1, SEXP g2, SEXP f, SEXP k) {
  clab::who = "clab_kinidx";
  const int n = clab::ilen(static_cast<double>(XLENGTH(g1) / 2), nullptr);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  check(clab_kinidx(INTEGER(g1), INTEGER(g2), REAL(f), n, REAL(k), REAL(out)));
  UNPROTECT(1);
  return out;
}

SEXP C_clab_patidx(SEXP c, SEXP m, SEXP af, SEXP f) {
  clab::who = "clab_patidx";
  const int n = clab::ilen(static_cast<double>(XLENGTH(c) / 2), nullptr);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  check(clab_patidx(INTEGER(c), INTEGER(m), INTEGER(af), REAL(f), n,
                    REAL(out)));
  UNPROTECT(1);
  return out;
}

SEXP C_clab_clrsmp(SEXP clra, SEXP reca, SEXP clrb, SEXP recb) {
  clab::who = "clab_clrsmp";
  const int n = len(clra);
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
  clab::who = "clab_expbnd";
  const int n = len(y);
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
  clab::who = "clab_misrd";
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 5));
  check(clab_misrd(Rf_asReal(r1), Rf_asReal(r0), Rf_asReal(p), Rf_asReal(se),
                   Rf_asReal(sp), REAL(out)));
  UNPROTECT(1);
  return out;
}

SEXP C_clab_rtmhs(SEXP x1, SEXP x2, SEXP w, SEXP c) {
  clab::who = "clab_rtmhs";
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 7));
  check(clab_rtmhs(REAL(x1), REAL(x2), REAL(w), len(x1),
                   Rf_asReal(c), REAL(out)));
  UNPROTECT(1);
  return out;
}

SEXP C_clab_incap(SEXP lam, SEXP q, SEXP s, SEXP share) {
  clab::who = "clab_incap";
  const int n = len(lam);
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
  clab::who = "clab_incapr";
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 5));
  check(clab_incapr(REAL(lam), len(lam), Rf_asInteger(t0),
                    Rf_asInteger(s), Rf_asReal(r), REAL(out)));
  UNPROTECT(1);
  return out;
}

SEXP C_clab_sellab(SEXP y, SEXP rel, SEXP rule, SEXP w) {
  clab::who = "clab_sellab";
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 6));
  check(clab_sellab(REAL(y), LOGICAL(rel), LOGICAL(rule), REAL(w),
                    len(y), REAL(out)));
  UNPROTECT(1);
  return out;
}

SEXP C_clab_losbnd(SEXP t, SEXP a) {
  clab::who = "clab_losbnd";
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 4));
  check(clab_losbnd(REAL(t), len(t), REAL(a),
                    len(a), REAL(out)));
  UNPROTECT(1);
  return out;
}

SEXP C_clab_sntbnd(SEXP y, SEXP z, SEXP w) {
  clab::who = "clab_sntbnd";
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 11));
  check(clab_sntbnd(REAL(y), INTEGER(z), REAL(w), len(y),
                    REAL(out)));
  UNPROTECT(1);
  return out;
}

SEXP C_clab_hzdsel(SEXP s, SEXP h, SEXP l, SEXP ah, SEXP al) {
  clab::who = "clab_hzdsel";
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
// The entry points of the exported functions (bind.cpp).
SEXP C_clab_r_vodtst(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_orrisk(SEXP, SEXP, SEXP);
SEXP C_clab_r_outtst(SEXP, SEXP, SEXP);
SEXP C_clab_r_dspbm(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_dyadrt(SEXP, SEXP);
SEXP C_clab_r_gtmp(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_cllr(SEXP, SEXP);
SEXP C_clab_r_lrpost(SEXP, SEXP);
SEXP C_clab_r_mixlr(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_dbsrch(SEXP, SEXP, SEXP);
SEXP C_clab_r_cllrmin(SEXP, SEXP);
SEXP C_clab_r_lrcal(SEXP, SEXP);
SEXP C_clab_r_tipplt(SEXP, SEXP);
SEXP C_clab_r_eceplt(SEXP, SEXP, SEXP);
SEXP C_clab_r_kinidx(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_patidx(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_clrsmp(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_expbnd(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_misrd(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_rtmhs(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_incap(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_incapr(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_sellab(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_losbnd(SEXP, SEXP);
SEXP C_clab_r_sntbnd(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_hzdsel(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_proof(SEXP);
SEXP C_clab_tlen(SEXP, SEXP);

static const R_CallMethodDef calls[] = {
    {"C_clab_vod", fn(&C_clab_vod), 4},
    {"C_clab_orrisk", fn(&C_clab_orrisk), 3},
    {"C_clab_outtst", fn(&C_clab_outtst), 3},
    {"C_clab_dspbm", fn(&C_clab_dspbm), 4},
    {"C_clab_dyadrt", fn(&C_clab_dyadrt), 4},
    {"C_clab_gtmp", fn(&C_clab_gtmp), 4},
    {"C_clab_cllr", fn(&C_clab_cllr), 2},
    {"C_clab_lrcal", fn(&C_clab_lrcal), 2},
    {"C_clab_mixp", fn(&C_clab_mixp), 3},
    {"C_clab_kinidx", fn(&C_clab_kinidx), 4},
    {"C_clab_patidx", fn(&C_clab_patidx), 4},
    {"C_clab_clrsmp", fn(&C_clab_clrsmp), 4},
    {"C_clab_expbnd", fn(&C_clab_expbnd), 4},
    {"C_clab_misrd", fn(&C_clab_misrd), 5},
    {"C_clab_rtmhs", fn(&C_clab_rtmhs), 4},
    {"C_clab_incap", fn(&C_clab_incap), 4},
    {"C_clab_incapr", fn(&C_clab_incapr), 4},
    {"C_clab_sellab", fn(&C_clab_sellab), 4},
    {"C_clab_losbnd", fn(&C_clab_losbnd), 2},
    {"C_clab_sntbnd", fn(&C_clab_sntbnd), 3},
    {"C_clab_hzdsel", fn(&C_clab_hzdsel), 5},
    {"C_clab_selftest", fn(&C_clab_selftest), 0},
    {"C_clab_r_vodtst", fn(&C_clab_r_vodtst), 4},
    {"C_clab_r_orrisk", fn(&C_clab_r_orrisk), 3},
    {"C_clab_r_outtst", fn(&C_clab_r_outtst), 3},
    {"C_clab_r_dspbm", fn(&C_clab_r_dspbm), 4},
    {"C_clab_r_dyadrt", fn(&C_clab_r_dyadrt), 2},
    {"C_clab_r_gtmp", fn(&C_clab_r_gtmp), 4},
    {"C_clab_r_cllr", fn(&C_clab_r_cllr), 2},
    {"C_clab_r_lrpost", fn(&C_clab_r_lrpost), 2},
    {"C_clab_r_mixlr", fn(&C_clab_r_mixlr), 5},
    {"C_clab_r_dbsrch", fn(&C_clab_r_dbsrch), 3},
    {"C_clab_r_cllrmin", fn(&C_clab_r_cllrmin), 2},
    {"C_clab_r_lrcal", fn(&C_clab_r_lrcal), 2},
    {"C_clab_r_tipplt", fn(&C_clab_r_tipplt), 2},
    {"C_clab_r_eceplt", fn(&C_clab_r_eceplt), 3},
    {"C_clab_r_kinidx", fn(&C_clab_r_kinidx), 4},
    {"C_clab_r_patidx", fn(&C_clab_r_patidx), 4},
    {"C_clab_r_clrsmp", fn(&C_clab_r_clrsmp), 4},
    {"C_clab_r_expbnd", fn(&C_clab_r_expbnd), 4},
    {"C_clab_r_misrd", fn(&C_clab_r_misrd), 5},
    {"C_clab_r_rtmhs", fn(&C_clab_r_rtmhs), 4},
    {"C_clab_r_incap", fn(&C_clab_r_incap), 4},
    {"C_clab_r_incapr", fn(&C_clab_r_incapr), 4},
    {"C_clab_r_sellab", fn(&C_clab_r_sellab), 4},
    {"C_clab_r_losbnd", fn(&C_clab_r_losbnd), 2},
    {"C_clab_r_sntbnd", fn(&C_clab_r_sntbnd), 4},
    {"C_clab_r_hzdsel", fn(&C_clab_r_hzdsel), 5},
    {"C_clab_r_proof", fn(&C_clab_r_proof), 1},
    {"C_clab_tlen", fn(&C_clab_tlen), 2},
    {nullptr, nullptr, 0}};

void R_init_crimlab(DllInfo *dll) {
  R_registerRoutines(dll, nullptr, calls, nullptr, nullptr);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}

}  // extern "C"
