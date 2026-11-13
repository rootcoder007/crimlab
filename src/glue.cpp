// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
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
  static const char *const names[] = {"shday", "shnight", "ratio",
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
  static const char *const names[] = {"rtoff", "rtpair", "nulloff",
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
SEXP C_clab_r_ntdeg(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_ntbtw(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_ntcls(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_ntevc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_ntcore(SEXP, SEXP, SEXP);
SEXP C_clab_r_ntcmp(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_ntlvn(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_ntmod(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_ntbrk(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_ntchg(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_coofg(SEXP, SEXP);
// The entry points of the spatial functions (bsp.cpp).
SEXP C_clab_r_hvdist(SEXP, SEXP, SEXP);
SEXP C_clab_r_eqproj(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_eqinv(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_cgstat(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_cgmedn(SEXP, SEXP, SEXP);
SEXP C_clab_r_webpt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_sdell(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_cvhull(SEXP, SEXP, SEXP);
SEXP C_clab_r_nnidx(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_ripk(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_ripenv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_fzmode(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_kmpart(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                     SEXP);
SEXP C_clab_r_nnhclu(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                     SEXP, SEXP, SEXP);
SEXP C_clab_r_stacel(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                     SEXP, SEXP);
SEXP C_clab_r_mantst(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_stmavg(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_jtdecay(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_jtsurf(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                     SEXP);
SEXP C_clab_r_jtcal(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_jtacc(SEXP, SEXP, SEXP);
SEXP C_clab_r_bjtc(SEXP, SEXP, SEXP);
// The entry points of the dossier functions (bdos.cpp).
SEXP C_clab_r_fair(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_fairdf(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_tmpaud(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_lipaud(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_nsbnd(SEXP, SEXP, SEXP);
SEXP C_clab_r_impfpr(SEXP, SEXP, SEXP);
SEXP C_clab_r_rnkres(SEXP, SEXP);
SEXP C_clab_r_calaud(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_scrdsp(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_nordet(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_fairsim(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_oaxbd(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_tlvc(SEXP, SEXP);
SEXP C_clab_r_tlvlr(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_aorwt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                    SEXP);
SEXP C_clab_r_bdig(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_bnft(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_skp(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_rptc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_rswp(SEXP, SEXP, SEXP);
SEXP C_clab_r_rshl(SEXP, SEXP);
SEXP C_clab_r_zgll(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_zgem(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_mvtll(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_mvtem(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_tfrp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_lwss(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_chmb(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_gstr(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_sptst(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
// The entry points of the weights and autocorrelation functions (bsp2.cpp).
SEXP C_clab_r_spwt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_gearc(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_gorg(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_lmorn(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_locgi(SEXP, SEXP, SEXP);
SEXP C_clab_r_dirmn(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_relchk(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_catchk(SEXP, SEXP, SEXP);
// The count and table statistics (bstat.cpp).
SEXP C_clab_r_cntreg(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_grprr(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_dastd(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_sirex(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_ctabx(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
// The learning and causal families (rlearn.cpp).
SEXP C_clab_r_enpath(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                     SEXP);
SEXP C_clab_r_lgtreg(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_rfor(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                   SEXP);
SEXP C_clab_r_dmlx(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                   SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                   SEXP, SEXP, SEXP);
SEXP C_clab_r_hdbang(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_crwalk(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_rnnhcl(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                     SEXP);
// The policy-evaluation families (bdid.cpp).
SEXP C_clab_r_didfe(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                    SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_evtst(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_stgdid(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                     SEXP);
SEXP C_clab_r_synthc(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_sdidx(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_wcbx(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                   SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_wystep(SEXP, SEXP);
SEXP C_clab_r_knxnr(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_kdsurf(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_dbscn(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_hsper(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_ebshr(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_kscan(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_gmoran(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_bvmor(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_sarml(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_spill(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_crossk(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_qdrtst(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_morsta(SEXP, SEXP);
SEXP C_clab_r_dissix(SEXP, SEXP, SEXP);
SEXP C_clab_r_pinpl(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_arwt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_bufas(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
// Time, contagion and operations (btm.cpp).
SEXP C_clab_r_hwkfit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_hwkres(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_hwksim(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_trndrr(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_stpcp(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_cusrv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_xcser(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_kmrec(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_hypq(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_sqrtlw(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_sdbsim(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                     SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_sdbthr(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_etasf(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_mvhwk(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_sthwk(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
// The design-based causal estimators (bcba.cpp).
SEXP C_clab_r_cbatt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_smdbal(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_evlx(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
// Survey estimates (bsvy.cpp).
SEXP C_clab_r_svymx(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
// The custody, court and published-table families (bsev.cpp and siblings).
SEXP C_clab_r_crsev(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_zcomp(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_cconc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_hilltl(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_hzta(SEXP, SEXP);
SEXP C_clab_r_crcap(SEXP, SEXP, SEXP);
SEXP C_clab_r_dkfig(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_bandx(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_pbbnd(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                    SEXP, SEXP);
SEXP C_clab_r_yoypx(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                    SEXP);
SEXP C_clab_r_jdgiv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_jmono(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_mmshr(SEXP, SEXP, SEXP);
SEXP C_clab_r_agemix(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_cumrsk(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_jpipe(SEXP, SEXP, SEXP, SEXP, SEXP);

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
    {"C_clab_r_ntdeg", fn(&C_clab_r_ntdeg), 5},
    {"C_clab_r_ntbtw", fn(&C_clab_r_ntbtw), 5},
    {"C_clab_r_ntcls", fn(&C_clab_r_ntcls), 4},
    {"C_clab_r_ntevc", fn(&C_clab_r_ntevc), 6},
    {"C_clab_r_ntcore", fn(&C_clab_r_ntcore), 3},
    {"C_clab_r_ntcmp", fn(&C_clab_r_ntcmp), 4},
    {"C_clab_r_ntlvn", fn(&C_clab_r_ntlvn), 6},
    {"C_clab_r_ntmod", fn(&C_clab_r_ntmod), 6},
    {"C_clab_r_ntbrk", fn(&C_clab_r_ntbrk), 5},
    {"C_clab_r_ntchg", fn(&C_clab_r_ntchg), 6},
    {"C_clab_r_coofg", fn(&C_clab_r_coofg), 2},
    {"C_clab_r_hvdist", fn(&C_clab_r_hvdist), 3},
    {"C_clab_r_eqproj", fn(&C_clab_r_eqproj), 5},
    {"C_clab_r_eqinv", fn(&C_clab_r_eqinv), 5},
    {"C_clab_r_cgstat", fn(&C_clab_r_cgstat), 5},
    {"C_clab_r_cgmedn", fn(&C_clab_r_cgmedn), 3},
    {"C_clab_r_webpt", fn(&C_clab_r_webpt), 7},
    {"C_clab_r_sdell", fn(&C_clab_r_sdell), 6},
    {"C_clab_r_cvhull", fn(&C_clab_r_cvhull), 3},
    {"C_clab_r_nnidx", fn(&C_clab_r_nnidx), 8},
    {"C_clab_r_ripk", fn(&C_clab_r_ripk), 6},
    {"C_clab_r_ripenv", fn(&C_clab_r_ripenv), 7},
    {"C_clab_r_fzmode", fn(&C_clab_r_fzmode), 6},
    {"C_clab_r_kmpart", fn(&C_clab_r_kmpart), 10},
    {"C_clab_r_nnhclu", fn(&C_clab_r_nnhclu), 12},
    {"C_clab_r_stacel", fn(&C_clab_r_stacel), 11},
    {"C_clab_r_mantst", fn(&C_clab_r_mantst), 9},
    {"C_clab_r_stmavg", fn(&C_clab_r_stmavg), 5},
    {"C_clab_r_jtdecay", fn(&C_clab_r_jtdecay), 4},
    {"C_clab_r_jtsurf", fn(&C_clab_r_jtsurf), 10},
    {"C_clab_r_jtcal", fn(&C_clab_r_jtcal), 6},
    {"C_clab_r_jtacc", fn(&C_clab_r_jtacc), 3},
    {"C_clab_r_bjtc", fn(&C_clab_r_bjtc), 3},
    {"C_clab_r_fair", fn(&C_clab_r_fair), 6},
    {"C_clab_r_fairdf", fn(&C_clab_r_fairdf), 7},
    {"C_clab_r_tmpaud", fn(&C_clab_r_tmpaud), 5},
    {"C_clab_r_lipaud", fn(&C_clab_r_lipaud), 4},
    {"C_clab_r_nsbnd", fn(&C_clab_r_nsbnd), 3},
    {"C_clab_r_impfpr", fn(&C_clab_r_impfpr), 3},
    {"C_clab_r_rnkres", fn(&C_clab_r_rnkres), 2},
    {"C_clab_r_calaud", fn(&C_clab_r_calaud), 7},
    {"C_clab_r_scrdsp", fn(&C_clab_r_scrdsp), 4},
    {"C_clab_r_nordet", fn(&C_clab_r_nordet), 6},
    {"C_clab_r_fairsim", fn(&C_clab_r_fairsim), 7},
    {"C_clab_r_oaxbd", fn(&C_clab_r_oaxbd), 5},
    {"C_clab_r_tlvc", fn(&C_clab_r_tlvc), 2},
    {"C_clab_r_tlvlr", fn(&C_clab_r_tlvlr), 5},
    {"C_clab_r_aorwt", fn(&C_clab_r_aorwt), 11},
    {"C_clab_r_bdig", fn(&C_clab_r_bdig), 4},
    {"C_clab_r_bnft", fn(&C_clab_r_bnft), 4},
    {"C_clab_r_skp", fn(&C_clab_r_skp), 5},
    {"C_clab_r_rptc", fn(&C_clab_r_rptc), 6},
    {"C_clab_r_rswp", fn(&C_clab_r_rswp), 3},
    {"C_clab_r_rshl", fn(&C_clab_r_rshl), 2},
    {"C_clab_r_zgll", fn(&C_clab_r_zgll), 6},
    {"C_clab_r_zgem", fn(&C_clab_r_zgem), 8},
    {"C_clab_r_mvtll", fn(&C_clab_r_mvtll), 5},
    {"C_clab_r_mvtem", fn(&C_clab_r_mvtem), 7},
    {"C_clab_r_tfrp", fn(&C_clab_r_tfrp), 8},
    {"C_clab_r_lwss", fn(&C_clab_r_lwss), 5},
    {"C_clab_r_chmb", fn(&C_clab_r_chmb), 5},
    {"C_clab_r_gstr", fn(&C_clab_r_gstr), 5},
    {"C_clab_r_sptst", fn(&C_clab_r_sptst), 6},
    {"C_clab_r_spwt", fn(&C_clab_r_spwt), 10},
    {"C_clab_r_gearc", fn(&C_clab_r_gearc), 4},
    {"C_clab_r_gorg", fn(&C_clab_r_gorg), 4},
    {"C_clab_r_lmorn", fn(&C_clab_r_lmorn), 4},
    {"C_clab_r_locgi", fn(&C_clab_r_locgi), 3},
    {"C_clab_r_dirmn", fn(&C_clab_r_dirmn), 6},
    {"C_clab_r_relchk", fn(&C_clab_r_relchk), 8},
    {"C_clab_r_catchk", fn(&C_clab_r_catchk), 3},
    {"C_clab_r_cntreg", fn(&C_clab_r_cntreg), 9},
    {"C_clab_r_grprr", fn(&C_clab_r_grprr), 7},
    {"C_clab_r_dastd", fn(&C_clab_r_dastd), 8},
    {"C_clab_r_sirex", fn(&C_clab_r_sirex), 6},
    {"C_clab_r_ctabx", fn(&C_clab_r_ctabx), 7},
    {"C_clab_r_enpath", fn(&C_clab_r_enpath), 11},
    {"C_clab_r_lgtreg", fn(&C_clab_r_lgtreg), 4},
    {"C_clab_r_rfor", fn(&C_clab_r_rfor), 11},
    {"C_clab_r_dmlx", fn(&C_clab_r_dmlx), 23},
    {"C_clab_r_hdbang", fn(&C_clab_r_hdbang), 9},
    {"C_clab_r_crwalk", fn(&C_clab_r_crwalk), 6},
    {"C_clab_r_rnnhcl", fn(&C_clab_r_rnnhcl), 10},
    {"C_clab_r_didfe", fn(&C_clab_r_didfe), 16},
    {"C_clab_r_evtst", fn(&C_clab_r_evtst), 18},
    {"C_clab_r_stgdid", fn(&C_clab_r_stgdid), 13},
    {"C_clab_r_synthc", fn(&C_clab_r_synthc), 5},
    {"C_clab_r_sdidx", fn(&C_clab_r_sdidx), 7},
    {"C_clab_r_wcbx", fn(&C_clab_r_wcbx), 18},
    {"C_clab_r_wystep", fn(&C_clab_r_wystep), 2},
    {"C_clab_r_knxnr", fn(&C_clab_r_knxnr), 8},
    {"C_clab_r_kdsurf", fn(&C_clab_r_kdsurf), 10},
    {"C_clab_r_dbscn", fn(&C_clab_r_dbscn), 6},
    {"C_clab_r_hsper", fn(&C_clab_r_hsper), 7},
    {"C_clab_r_ebshr", fn(&C_clab_r_ebshr), 4},
    {"C_clab_r_kscan", fn(&C_clab_r_kscan), 10},
    {"C_clab_r_gmoran", fn(&C_clab_r_gmoran), 4},
    {"C_clab_r_bvmor", fn(&C_clab_r_bvmor), 5},
    {"C_clab_r_sarml", fn(&C_clab_r_sarml), 7},
    {"C_clab_r_spill", fn(&C_clab_r_spill), 6},
    {"C_clab_r_crossk", fn(&C_clab_r_crossk), 10},
    {"C_clab_r_qdrtst", fn(&C_clab_r_qdrtst), 5},
    {"C_clab_r_morsta", fn(&C_clab_r_morsta), 2},
    {"C_clab_r_dissix", fn(&C_clab_r_dissix), 3},
    {"C_clab_r_pinpl", fn(&C_clab_r_pinpl), 4},
    {"C_clab_r_arwt", fn(&C_clab_r_arwt), 6},
    {"C_clab_r_bufas", fn(&C_clab_r_bufas), 7},
    {"C_clab_r_hwkfit", fn(&C_clab_r_hwkfit), 9},
    {"C_clab_r_hwkres", fn(&C_clab_r_hwkres), 6},
    {"C_clab_r_hwksim", fn(&C_clab_r_hwksim), 9},
    {"C_clab_r_trndrr", fn(&C_clab_r_trndrr), 6},
    {"C_clab_r_stpcp", fn(&C_clab_r_stpcp), 4},
    {"C_clab_r_cusrv", fn(&C_clab_r_cusrv), 7},
    {"C_clab_r_xcser", fn(&C_clab_r_xcser), 5},
    {"C_clab_r_kmrec", fn(&C_clab_r_kmrec), 6},
    {"C_clab_r_hypq", fn(&C_clab_r_hypq), 5},
    {"C_clab_r_sqrtlw", fn(&C_clab_r_sqrtlw), 5},
    {"C_clab_r_sdbsim", fn(&C_clab_r_sdbsim), 17},
    {"C_clab_r_sdbthr", fn(&C_clab_r_sdbthr), 8},
    {"C_clab_r_etasf", fn(&C_clab_r_etasf), 6},
    {"C_clab_r_mvhwk", fn(&C_clab_r_mvhwk), 6},
    {"C_clab_r_sthwk", fn(&C_clab_r_sthwk), 9},
    {"C_clab_r_cbatt", fn(&C_clab_r_cbatt), 12},
    {"C_clab_r_smdbal", fn(&C_clab_r_smdbal), 5},
    {"C_clab_r_evlx", fn(&C_clab_r_evlx), 6},
    {"C_clab_r_svymx", fn(&C_clab_r_svymx), 10},
    {"C_clab_r_crsev", fn(&C_clab_r_crsev), 9},
    {"C_clab_r_zcomp", fn(&C_clab_r_zcomp), 4},
    {"C_clab_r_cconc", fn(&C_clab_r_cconc), 7},
    {"C_clab_r_hilltl", fn(&C_clab_r_hilltl), 5},
    {"C_clab_r_hzta", fn(&C_clab_r_hzta), 2},
    {"C_clab_r_crcap", fn(&C_clab_r_crcap), 3},
    {"C_clab_r_dkfig", fn(&C_clab_r_dkfig), 10},
    {"C_clab_r_bandx", fn(&C_clab_r_bandx), 8},
    {"C_clab_r_pbbnd", fn(&C_clab_r_pbbnd), 12},
    {"C_clab_r_yoypx", fn(&C_clab_r_yoypx), 11},
    {"C_clab_r_jdgiv", fn(&C_clab_r_jdgiv), 6},
    {"C_clab_r_jmono", fn(&C_clab_r_jmono), 10},
    {"C_clab_r_mmshr", fn(&C_clab_r_mmshr), 3},
    {"C_clab_r_agemix", fn(&C_clab_r_agemix), 6},
    {"C_clab_r_cumrsk", fn(&C_clab_r_cumrsk), 6},
    {"C_clab_r_jpipe", fn(&C_clab_r_jpipe), 5},
    {nullptr, nullptr, 0}};

void R_init_crimlab(DllInfo *dll) {
  R_registerRoutines(dll, nullptr, calls, nullptr, nullptr);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}

}  // extern "C"
