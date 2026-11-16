// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R binding of libcrimlab: each entry point unpacks R vectors, calls one
// kernel and packs the result. The R wrappers have checked every argument, so
// a non-zero code here means the two layers disagree; it is reported in the
// kernel's own words. No C++ object with a destructor is alive when
// Rf_error() jumps.
//
// Interrupts. The kernels poll rintr() (clab_onintr), which runs R's check
// inside R_ToplevelExec(): an interrupt there ends that check only, never
// jumping over a kernel's C++ frames, and the kernel returns CLAB_EINTR.
// The binding's error paths (check(), clab::fail(), clab::kchk()) then
// throw clab::Intr instead of raising an error, which unwinds the
// binding's frames (destructors run) up to the entry point. There Safe
// catches it, or finds the flag set after a binding that did not look at
// the code, drops the half-made result and raises R's interrupt, with no
// C++ object alive.
#include <R.h>
#include <Rinternals.h>

#include <cstdio>
#include <exception>
#include <new>
#include <stdexcept>

#include "crimlab.h"
#include "rchk.h"

namespace {

void rchk1(void *) { R_CheckUserInterrupt(); }
int rintr() { return R_ToplevelExec(rchk1, nullptr) ? 0 : 1; }
int yes() { return 1; }

// What a C++ exception that reached an entry point said: copied here, a
// plain buffer, so that the error is raised after the try has unwound.
char xmsg[256];

// An entry point F wrapped as above. No C++ exception crosses into R: a
// failed allocation (std::bad_alloc, or std::length_error for a vector
// past its largest size) and any other std::exception are caught here and
// raised, once the binding's frames are gone, as an ordinary error named
// after the function (clab::who).
template <auto F>
struct Safe;
template <class... A, SEXP (*F)(A...)>
struct Safe<F> {
  static SEXP call(A... a) {
    clab_intrd(1);
    bool stop = false;
    const char *bad = nullptr;
    SEXP r = R_NilValue;
    try {
      r = F(a...);
    } catch (const clab::Intr &) {
      stop = true;
    } catch (const std::bad_alloc &) {
      bad = "not enough memory for this size.";
    } catch (const std::length_error &) {
      bad = "not enough memory for this size.";
    } catch (const std::exception &e) {
      std::snprintf(xmsg, sizeof xmsg, "%s", e.what());
      bad = xmsg;
    } catch (...) {
      bad = "an unexpected internal error.";
    }
    const bool hit = clab_intrd(1) != 0 || stop;
    if (hit) Rf_onintr();  // returns only while R holds interrupts back
    if (bad != nullptr) Rf_errorcall(R_NilValue, "%s: %s", clab::who, bad);
    return hit ? R_NilValue : r;
  }
};

// An entry point as R's registration table holds it. The cast goes through
// void (*)(void), which compilers accept as matching every function type,
// so -Wcast-function-type (in -Wextra) stays quiet on (DL_FUNC) &f.
template <auto F>
DL_FUNC fn() {
  return reinterpret_cast<DL_FUNC>(reinterpret_cast<void (*)(void)>(&Safe<F>::call));
}

// The length of x for a kernel, which counts in int: a vector longer than
// INT_MAX is refused with the kernels' out-of-range error, not truncated.
int len(SEXP x) { return clab::ilen(static_cast<double>(XLENGTH(x)), nullptr); }

void check(int code) {
  if (code == CLAB_EINTR || clab_intrd(0) != 0) throw clab::Intr{};
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
SEXP C_clab_r_pedlr(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_gmixlr(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_cmpsc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
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
// The entry points of the specialist methods of R/spec.R (bdos.cpp).
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
                   SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_suplrn(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                     SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_hdbang(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_crwalk(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_rnnhcl(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                     SEXP);
// The environmental-crime family (benv.cpp).
SEXP C_clab_r_lmdix(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_fepnl(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_gmmpnl(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                     SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_wsctst(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_ccepnl(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                     SEXP, SEXP);
SEXP C_clab_r_pcdtst(SEXP, SEXP, SEXP);
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
SEXP C_clab_r_evlx(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
// The summary statistics (bprim.cpp).
SEXP C_clab_r_prim(SEXP, SEXP, SEXP, SEXP);
// Text to plain ASCII (btxt.cpp).
SEXP C_clab_r_ascii(SEXP);
// The numbers of the change-table writers (bfmt.cpp).
SEXP C_clab_r_fmt(SEXP, SEXP, SEXP);
// Bounds under monotone response and selection (bmtr.cpp).
SEXP C_clab_r_mtrbnd(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
// The Indigenous health, published-table and census functions (bigh.cpp).
SEXP C_clab_r_igh(SEXP, SEXP);
// The power and energy readers (bpwr.cpp).
SEXP C_clab_r_pwrs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_pwrd(SEXP, SEXP);
SEXP C_clab_r_pwrp(SEXP, SEXP, SEXP);
// The accuracy of a hot-spot map (bhsa.cpp).
SEXP C_clab_r_hsacc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
// Scaling, ideal points and counts over time (bsvt.cpp).
SEXP C_clab_r_svk(SEXP, SEXP);
// The age-period-cohort model (bapc.cpp).
SEXP C_clab_r_apcie(SEXP, SEXP, SEXP, SEXP);
// The split-population survival model (bspc.cpp).
SEXP C_clab_r_spcure(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
// Crime linkage (blnk.cpp).
SEXP C_clab_r_lnkpr(SEXP, SEXP, SEXP, SEXP, SEXP);
// Prospective hot-spotting (bprh.cpp).
SEXP C_clab_r_prohs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
// The weighted displacement quotient (bwdq.cpp).
SEXP C_clab_r_wdqx(SEXP, SEXP, SEXP);
// Spatial epidemiology (bsep.cpp).
SEXP C_clab_r_bayout(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_bufexp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_kerexp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_lerprc(SEXP, SEXP, SEXP);
SEXP C_clab_r_bym2sc(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_ecolrg(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_spbias(SEXP, SEXP, SEXP);
SEXP C_clab_r_spsir(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_sfrail(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_sfcure(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_evfilt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
// Geostatistics (bkrg.cpp).
SEXP C_clab_r_vgemp(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_vgfit(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_okrge(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_bkrig(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_krigcv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_trsurf(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_sglmx(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_smixm(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_lfspec(SEXP);
SEXP C_clab_r_spsamp(SEXP, SEXP, SEXP, SEXP);
// Repeat victimisation (brpv.cpp).
SEXP C_clab_r_rptvic(SEXP, SEXP, SEXP, SEXP);
// Eyewitness lineups (blnr.cpp).
SEXP C_clab_r_lnroc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_lnsize(SEXP, SEXP, SEXP, SEXP);
// Internal benchmarking of officers (bibn.cpp).
SEXP C_clab_r_ibench(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
// Discrimination in searches (bthr.cpp).
SEXP C_clab_r_thrtst(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
// Search theory for search and rescue (bsar.cpp).
SEXP C_clab_r_sarlr(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_sarsw(SEXP, SEXP);
SEXP C_clab_r_sarpod(SEXP, SEXP, SEXP);
SEXP C_clab_r_sarpsw(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_sarscs(SEXP, SEXP);
SEXP C_clab_r_sarmap(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_sarcrc(SEXP, SEXP, SEXP);
SEXP C_clab_r_saropt(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_sarsq(SEXP, SEXP, SEXP);
SEXP C_clab_r_sarfig(void);
// Air pollution (bapol.cpp).
SEXP C_clab_r_pgsig(SEXP, SEXP, SEXP);
SEXP C_clab_r_briggs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_gplume(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_gspuff(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_advdif(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_laggrid(SEXP, SEXP, SEXP);
SEXP C_clab_r_lagr(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_lagstep(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_crf(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_levin(SEXP, SEXP);
SEXP C_clab_r_hif(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_concidx(SEXP, SEXP);
// Exposure and crime (bapx.cpp).
SEXP C_clab_r_fesep(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_ppmlfe(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                     SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_ivfe(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                   SEXP, SEXP, SEXP);
SEXP C_clab_r_ppmlcf(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                     SEXP, SEXP, SEXP);
SEXP C_clab_r_demean(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_expgmm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                     SEXP, SEXP, SEXP);
SEXP C_clab_r_ivoid(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                    SEXP, SEXP, SEXP);
SEXP C_clab_r_rdck(SEXP);
SEXP C_clab_r_ikbw(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_mccr(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_windbin(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_dwnexp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                     SEXP, SEXP);
SEXP C_clab_r_tailwd(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_nbfit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_idwint(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_gtwreg(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_asrule(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_dshift(SEXP, SEXP, SEXP);
SEXP C_clab_r_updwn(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_vcx(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_runlen(SEXP, SEXP, SEXP, SEXP);
// One-pass summaries (bskch.cpp).
SEXP C_clab_r_crc32(SEXP, SEXP);
SEXP C_clab_r_momupd(SEXP, SEXP);
SEXP C_clab_r_mommrg(SEXP, SEXP);
SEXP C_clab_r_momstat(SEXP);
SEXP C_clab_r_rsv(SEXP, SEXP, SEXP);
SEXP C_clab_r_hlladd(SEXP, SEXP, SEXP);
SEXP C_clab_r_hllcnt(SEXP);
// The footprint of a computation (bfp.cpp).
SEXP C_clab_r_citab(void);
SEXP C_clab_r_cimatch(SEXP);
SEXP C_clab_r_tzterr(SEXP);
SEXP C_clab_r_proxy(SEXP);
SEXP C_clab_r_fpe(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_fpeqv(SEXP);
SEXP C_clab_r_fploc(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_cizone(SEXP);
SEXP C_clab_r_rapldom(SEXP);
SEXP C_clab_r_raplrd(SEXP);
SEXP C_clab_r_cpuacc(SEXP);
SEXP C_clab_r_memgb(void);
SEXP C_clab_r_fpnote(SEXP, SEXP);
// The exposure-response estimate (bplr.cpp).
SEXP C_clab_r_plr(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_bootci(SEXP, SEXP);
// Data description (bdsc.cpp).
SEXP C_clab_r_hbin(SEXP, SEXP);
SEXP C_clab_r_freq(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_pcor(SEXP, SEXP, SEXP);
SEXP C_clab_r_rowany(SEXP);
SEXP C_clab_r_nrun(SEXP, SEXP);
SEXP C_clab_r_msum(SEXP);
SEXP C_clab_r_pstat(SEXP);
SEXP C_clab_r_q7(SEXP, SEXP);
SEXP C_clab_r_maha(SEXP, SEXP, SEXP);
SEXP C_clab_r_mcar(SEXP, SEXP, SEXP);
SEXP C_clab_r_mgrid(SEXP, SEXP);
SEXP C_clab_r_dups(SEXP, SEXP);
SEXP C_clab_r_cochq(SEXP, SEXP);
// Distributions fitted by maximum likelihood (bdfit.cpp).
SEXP C_clab_r_dfit(SEXP, SEXP);
SEXP C_clab_r_dfgof(SEXP, SEXP, SEXP);
SEXP C_clab_r_dfboot(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_dfhill(SEXP, SEXP, SEXP);
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
SEXP C_clab_r_mndla(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                    SEXP, SEXP);
SEXP C_clab_r_stkflw(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                     SEXP);
SEXP C_clab_r_lostay(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_bklog(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_dktflw(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                     SEXP, SEXP);
SEXP C_clab_r_lftab(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_trnmat(SEXP, SEXP, SEXP);
// Areas (bxw.cpp).
SEXP C_clab_r_bndxw(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_synexp(SEXP, SEXP, SEXP, SEXP);
// Deterrence theory (bdtr.cpp).
SEXP C_clab_r_detrsp(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_dtrid(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_inspg(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_insphs(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_colarr(SEXP, SEXP, SEXP);
// Taphonomy (btph.cpp).
SEXP C_clab_r_clrilr(SEXP, SEXP, SEXP);
SEXP C_clab_r_preslr(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_tpdmc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_tpbhm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_dirsim(SEXP, SEXP, SEXP, SEXP);
// Forensic toxicology (btox.cpp).
SEXP C_clab_r_txcal(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_lodimp(SEXP, SEXP, SEXP);
SEXP C_clab_r_pmrrat(SEXP, SEXP, SEXP);
SEXP C_clab_r_etohlr(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_mtxrel(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
// Geographic profiling and risk terrain (brtm.cpp).
SEXP C_clab_r_crchyp(SEXP, SEXP, SEXP);
SEXP C_clab_r_hitscr(SEXP, SEXP);
SEXP C_clab_r_rtmfit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_crfpar(SEXP, SEXP);
SEXP C_clab_r_pbarea(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_rocpts(SEXP, SEXP);
// Data description between arguments and results (bdscd.cpp).
SEXP C_clab_r_spark(SEXP, SEXP, SEXP);
SEXP C_clab_r_cnames(SEXP, SEXP, SEXP);
SEXP C_clab_r_mgridl(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_mpat(SEXP, SEXP, SEXP);
SEXP C_clab_r_prof(SEXP, SEXP, SEXP);
SEXP C_clab_r_mhtab(SEXP, SEXP, SEXP);
SEXP C_clab_r_mruns(SEXP, SEXP, SEXP);
SEXP C_clab_r_topcor(SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_dtsel(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_clab_r_mcar2(SEXP, SEXP, SEXP);
// The exposure file of vrfpol() (bvrfp.cpp).
SEXP C_clab_r_bz2end(SEXP);
SEXP C_clab_r_gzend(SEXP, SEXP);
SEXP C_clab_r_txtenc(SEXP);
SEXP C_clab_r_rdfit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

// For the tests: on TRUE every poll of a kernel reads as an interrupt
// (the stop paths run on small inputs, in this process); on FALSE R's
// check is back, at its interval.
SEXP C_clab_r_intrt(SEXP on) {
  if (Rf_asLogical(on) == TRUE) {
    clab_onintr(&yes, 0.0);
  } else {
    clab_onintr(&rintr, -1.0);
  }
  return R_NilValue;
}

static const R_CallMethodDef calls[] = {
    {"C_clab_vod", fn<&C_clab_vod>(), 4},
    {"C_clab_orrisk", fn<&C_clab_orrisk>(), 3},
    {"C_clab_outtst", fn<&C_clab_outtst>(), 3},
    {"C_clab_dspbm", fn<&C_clab_dspbm>(), 4},
    {"C_clab_dyadrt", fn<&C_clab_dyadrt>(), 4},
    {"C_clab_gtmp", fn<&C_clab_gtmp>(), 4},
    {"C_clab_cllr", fn<&C_clab_cllr>(), 2},
    {"C_clab_lrcal", fn<&C_clab_lrcal>(), 2},
    {"C_clab_mixp", fn<&C_clab_mixp>(), 3},
    {"C_clab_kinidx", fn<&C_clab_kinidx>(), 4},
    {"C_clab_patidx", fn<&C_clab_patidx>(), 4},
    {"C_clab_clrsmp", fn<&C_clab_clrsmp>(), 4},
    {"C_clab_expbnd", fn<&C_clab_expbnd>(), 4},
    {"C_clab_misrd", fn<&C_clab_misrd>(), 5},
    {"C_clab_rtmhs", fn<&C_clab_rtmhs>(), 4},
    {"C_clab_incap", fn<&C_clab_incap>(), 4},
    {"C_clab_incapr", fn<&C_clab_incapr>(), 4},
    {"C_clab_sellab", fn<&C_clab_sellab>(), 4},
    {"C_clab_losbnd", fn<&C_clab_losbnd>(), 2},
    {"C_clab_sntbnd", fn<&C_clab_sntbnd>(), 3},
    {"C_clab_hzdsel", fn<&C_clab_hzdsel>(), 5},
    {"C_clab_selftest", fn<&C_clab_selftest>(), 0},
    {"C_clab_r_vodtst", fn<&C_clab_r_vodtst>(), 4},
    {"C_clab_r_orrisk", fn<&C_clab_r_orrisk>(), 3},
    {"C_clab_r_outtst", fn<&C_clab_r_outtst>(), 3},
    {"C_clab_r_dspbm", fn<&C_clab_r_dspbm>(), 4},
    {"C_clab_r_dyadrt", fn<&C_clab_r_dyadrt>(), 2},
    {"C_clab_r_gtmp", fn<&C_clab_r_gtmp>(), 4},
    {"C_clab_r_cllr", fn<&C_clab_r_cllr>(), 2},
    {"C_clab_r_lrpost", fn<&C_clab_r_lrpost>(), 2},
    {"C_clab_r_mixlr", fn<&C_clab_r_mixlr>(), 5},
    {"C_clab_r_dbsrch", fn<&C_clab_r_dbsrch>(), 3},
    {"C_clab_r_cllrmin", fn<&C_clab_r_cllrmin>(), 2},
    {"C_clab_r_lrcal", fn<&C_clab_r_lrcal>(), 2},
    {"C_clab_r_tipplt", fn<&C_clab_r_tipplt>(), 2},
    {"C_clab_r_eceplt", fn<&C_clab_r_eceplt>(), 3},
    {"C_clab_r_kinidx", fn<&C_clab_r_kinidx>(), 4},
    {"C_clab_r_patidx", fn<&C_clab_r_patidx>(), 4},
    {"C_clab_r_pedlr", fn<&C_clab_r_pedlr>(), 7},
    {"C_clab_r_gmixlr", fn<&C_clab_r_gmixlr>(), 13},
    {"C_clab_r_cmpsc", fn<&C_clab_r_cmpsc>(), 6},
    {"C_clab_r_clrsmp", fn<&C_clab_r_clrsmp>(), 4},
    {"C_clab_r_expbnd", fn<&C_clab_r_expbnd>(), 4},
    {"C_clab_r_misrd", fn<&C_clab_r_misrd>(), 5},
    {"C_clab_r_rtmhs", fn<&C_clab_r_rtmhs>(), 4},
    {"C_clab_r_incap", fn<&C_clab_r_incap>(), 4},
    {"C_clab_r_incapr", fn<&C_clab_r_incapr>(), 4},
    {"C_clab_r_sellab", fn<&C_clab_r_sellab>(), 4},
    {"C_clab_r_losbnd", fn<&C_clab_r_losbnd>(), 2},
    {"C_clab_r_sntbnd", fn<&C_clab_r_sntbnd>(), 4},
    {"C_clab_r_hzdsel", fn<&C_clab_r_hzdsel>(), 5},
    {"C_clab_r_proof", fn<&C_clab_r_proof>(), 1},
    {"C_clab_tlen", fn<&C_clab_tlen>(), 2},
    {"C_clab_r_ntdeg", fn<&C_clab_r_ntdeg>(), 5},
    {"C_clab_r_ntbtw", fn<&C_clab_r_ntbtw>(), 5},
    {"C_clab_r_ntcls", fn<&C_clab_r_ntcls>(), 4},
    {"C_clab_r_ntevc", fn<&C_clab_r_ntevc>(), 6},
    {"C_clab_r_ntcore", fn<&C_clab_r_ntcore>(), 3},
    {"C_clab_r_ntcmp", fn<&C_clab_r_ntcmp>(), 4},
    {"C_clab_r_ntlvn", fn<&C_clab_r_ntlvn>(), 6},
    {"C_clab_r_ntmod", fn<&C_clab_r_ntmod>(), 6},
    {"C_clab_r_ntbrk", fn<&C_clab_r_ntbrk>(), 5},
    {"C_clab_r_ntchg", fn<&C_clab_r_ntchg>(), 6},
    {"C_clab_r_coofg", fn<&C_clab_r_coofg>(), 2},
    {"C_clab_r_hvdist", fn<&C_clab_r_hvdist>(), 3},
    {"C_clab_r_eqproj", fn<&C_clab_r_eqproj>(), 5},
    {"C_clab_r_eqinv", fn<&C_clab_r_eqinv>(), 5},
    {"C_clab_r_cgstat", fn<&C_clab_r_cgstat>(), 5},
    {"C_clab_r_cgmedn", fn<&C_clab_r_cgmedn>(), 3},
    {"C_clab_r_webpt", fn<&C_clab_r_webpt>(), 7},
    {"C_clab_r_sdell", fn<&C_clab_r_sdell>(), 6},
    {"C_clab_r_cvhull", fn<&C_clab_r_cvhull>(), 3},
    {"C_clab_r_nnidx", fn<&C_clab_r_nnidx>(), 8},
    {"C_clab_r_ripk", fn<&C_clab_r_ripk>(), 6},
    {"C_clab_r_ripenv", fn<&C_clab_r_ripenv>(), 7},
    {"C_clab_r_fzmode", fn<&C_clab_r_fzmode>(), 6},
    {"C_clab_r_kmpart", fn<&C_clab_r_kmpart>(), 10},
    {"C_clab_r_nnhclu", fn<&C_clab_r_nnhclu>(), 12},
    {"C_clab_r_stacel", fn<&C_clab_r_stacel>(), 11},
    {"C_clab_r_mantst", fn<&C_clab_r_mantst>(), 9},
    {"C_clab_r_stmavg", fn<&C_clab_r_stmavg>(), 5},
    {"C_clab_r_jtdecay", fn<&C_clab_r_jtdecay>(), 4},
    {"C_clab_r_jtsurf", fn<&C_clab_r_jtsurf>(), 10},
    {"C_clab_r_jtcal", fn<&C_clab_r_jtcal>(), 6},
    {"C_clab_r_jtacc", fn<&C_clab_r_jtacc>(), 3},
    {"C_clab_r_bjtc", fn<&C_clab_r_bjtc>(), 3},
    {"C_clab_r_fair", fn<&C_clab_r_fair>(), 6},
    {"C_clab_r_fairdf", fn<&C_clab_r_fairdf>(), 7},
    {"C_clab_r_tmpaud", fn<&C_clab_r_tmpaud>(), 5},
    {"C_clab_r_lipaud", fn<&C_clab_r_lipaud>(), 4},
    {"C_clab_r_nsbnd", fn<&C_clab_r_nsbnd>(), 3},
    {"C_clab_r_impfpr", fn<&C_clab_r_impfpr>(), 3},
    {"C_clab_r_rnkres", fn<&C_clab_r_rnkres>(), 2},
    {"C_clab_r_calaud", fn<&C_clab_r_calaud>(), 7},
    {"C_clab_r_scrdsp", fn<&C_clab_r_scrdsp>(), 4},
    {"C_clab_r_nordet", fn<&C_clab_r_nordet>(), 6},
    {"C_clab_r_fairsim", fn<&C_clab_r_fairsim>(), 7},
    {"C_clab_r_oaxbd", fn<&C_clab_r_oaxbd>(), 5},
    {"C_clab_r_tlvc", fn<&C_clab_r_tlvc>(), 2},
    {"C_clab_r_tlvlr", fn<&C_clab_r_tlvlr>(), 5},
    {"C_clab_r_aorwt", fn<&C_clab_r_aorwt>(), 11},
    {"C_clab_r_bdig", fn<&C_clab_r_bdig>(), 4},
    {"C_clab_r_bnft", fn<&C_clab_r_bnft>(), 4},
    {"C_clab_r_skp", fn<&C_clab_r_skp>(), 5},
    {"C_clab_r_rptc", fn<&C_clab_r_rptc>(), 6},
    {"C_clab_r_rswp", fn<&C_clab_r_rswp>(), 3},
    {"C_clab_r_rshl", fn<&C_clab_r_rshl>(), 2},
    {"C_clab_r_zgll", fn<&C_clab_r_zgll>(), 6},
    {"C_clab_r_zgem", fn<&C_clab_r_zgem>(), 8},
    {"C_clab_r_mvtll", fn<&C_clab_r_mvtll>(), 5},
    {"C_clab_r_mvtem", fn<&C_clab_r_mvtem>(), 7},
    {"C_clab_r_tfrp", fn<&C_clab_r_tfrp>(), 8},
    {"C_clab_r_lwss", fn<&C_clab_r_lwss>(), 5},
    {"C_clab_r_chmb", fn<&C_clab_r_chmb>(), 5},
    {"C_clab_r_gstr", fn<&C_clab_r_gstr>(), 5},
    {"C_clab_r_sptst", fn<&C_clab_r_sptst>(), 6},
    {"C_clab_r_spwt", fn<&C_clab_r_spwt>(), 10},
    {"C_clab_r_gearc", fn<&C_clab_r_gearc>(), 4},
    {"C_clab_r_gorg", fn<&C_clab_r_gorg>(), 4},
    {"C_clab_r_lmorn", fn<&C_clab_r_lmorn>(), 4},
    {"C_clab_r_locgi", fn<&C_clab_r_locgi>(), 3},
    {"C_clab_r_dirmn", fn<&C_clab_r_dirmn>(), 6},
    {"C_clab_r_relchk", fn<&C_clab_r_relchk>(), 8},
    {"C_clab_r_catchk", fn<&C_clab_r_catchk>(), 3},
    {"C_clab_r_cntreg", fn<&C_clab_r_cntreg>(), 9},
    {"C_clab_r_grprr", fn<&C_clab_r_grprr>(), 7},
    {"C_clab_r_dastd", fn<&C_clab_r_dastd>(), 8},
    {"C_clab_r_sirex", fn<&C_clab_r_sirex>(), 6},
    {"C_clab_r_ctabx", fn<&C_clab_r_ctabx>(), 7},
    {"C_clab_r_enpath", fn<&C_clab_r_enpath>(), 11},
    {"C_clab_r_lgtreg", fn<&C_clab_r_lgtreg>(), 4},
    {"C_clab_r_rfor", fn<&C_clab_r_rfor>(), 11},
    {"C_clab_r_dmlx", fn<&C_clab_r_dmlx>(), 26},
    {"C_clab_r_suplrn", fn<&C_clab_r_suplrn>(), 14},
    {"C_clab_r_hdbang", fn<&C_clab_r_hdbang>(), 9},
    {"C_clab_r_crwalk", fn<&C_clab_r_crwalk>(), 6},
    {"C_clab_r_rnnhcl", fn<&C_clab_r_rnnhcl>(), 10},
    {"C_clab_r_lmdix", fn<&C_clab_r_lmdix>(), 9},
    {"C_clab_r_fepnl", fn<&C_clab_r_fepnl>(), 9},
    {"C_clab_r_gmmpnl", fn<&C_clab_r_gmmpnl>(), 22},
    {"C_clab_r_wsctst", fn<&C_clab_r_wsctst>(), 6},
    {"C_clab_r_ccepnl", fn<&C_clab_r_ccepnl>(), 14},
    {"C_clab_r_pcdtst", fn<&C_clab_r_pcdtst>(), 3},
    {"C_clab_r_didfe", fn<&C_clab_r_didfe>(), 16},
    {"C_clab_r_evtst", fn<&C_clab_r_evtst>(), 18},
    {"C_clab_r_stgdid", fn<&C_clab_r_stgdid>(), 13},
    {"C_clab_r_synthc", fn<&C_clab_r_synthc>(), 5},
    {"C_clab_r_sdidx", fn<&C_clab_r_sdidx>(), 7},
    {"C_clab_r_wcbx", fn<&C_clab_r_wcbx>(), 18},
    {"C_clab_r_wystep", fn<&C_clab_r_wystep>(), 2},
    {"C_clab_r_knxnr", fn<&C_clab_r_knxnr>(), 8},
    {"C_clab_r_kdsurf", fn<&C_clab_r_kdsurf>(), 10},
    {"C_clab_r_dbscn", fn<&C_clab_r_dbscn>(), 6},
    {"C_clab_r_hsper", fn<&C_clab_r_hsper>(), 7},
    {"C_clab_r_ebshr", fn<&C_clab_r_ebshr>(), 4},
    {"C_clab_r_kscan", fn<&C_clab_r_kscan>(), 10},
    {"C_clab_r_gmoran", fn<&C_clab_r_gmoran>(), 4},
    {"C_clab_r_bvmor", fn<&C_clab_r_bvmor>(), 5},
    {"C_clab_r_sarml", fn<&C_clab_r_sarml>(), 7},
    {"C_clab_r_spill", fn<&C_clab_r_spill>(), 6},
    {"C_clab_r_crossk", fn<&C_clab_r_crossk>(), 10},
    {"C_clab_r_qdrtst", fn<&C_clab_r_qdrtst>(), 5},
    {"C_clab_r_morsta", fn<&C_clab_r_morsta>(), 2},
    {"C_clab_r_dissix", fn<&C_clab_r_dissix>(), 3},
    {"C_clab_r_pinpl", fn<&C_clab_r_pinpl>(), 4},
    {"C_clab_r_arwt", fn<&C_clab_r_arwt>(), 6},
    {"C_clab_r_bufas", fn<&C_clab_r_bufas>(), 7},
    {"C_clab_r_hwkfit", fn<&C_clab_r_hwkfit>(), 9},
    {"C_clab_r_hwkres", fn<&C_clab_r_hwkres>(), 6},
    {"C_clab_r_hwksim", fn<&C_clab_r_hwksim>(), 9},
    {"C_clab_r_trndrr", fn<&C_clab_r_trndrr>(), 6},
    {"C_clab_r_stpcp", fn<&C_clab_r_stpcp>(), 4},
    {"C_clab_r_cusrv", fn<&C_clab_r_cusrv>(), 7},
    {"C_clab_r_xcser", fn<&C_clab_r_xcser>(), 5},
    {"C_clab_r_kmrec", fn<&C_clab_r_kmrec>(), 6},
    {"C_clab_r_hypq", fn<&C_clab_r_hypq>(), 5},
    {"C_clab_r_sqrtlw", fn<&C_clab_r_sqrtlw>(), 5},
    {"C_clab_r_sdbsim", fn<&C_clab_r_sdbsim>(), 17},
    {"C_clab_r_sdbthr", fn<&C_clab_r_sdbthr>(), 8},
    {"C_clab_r_etasf", fn<&C_clab_r_etasf>(), 6},
    {"C_clab_r_mvhwk", fn<&C_clab_r_mvhwk>(), 6},
    {"C_clab_r_sthwk", fn<&C_clab_r_sthwk>(), 9},
    {"C_clab_r_cbatt", fn<&C_clab_r_cbatt>(), 12},
    {"C_clab_r_smdbal", fn<&C_clab_r_smdbal>(), 5},
    {"C_clab_r_evlx", fn<&C_clab_r_evlx>(), 7},
    {"C_clab_r_mtrbnd", fn<&C_clab_r_mtrbnd>(), 6},
    {"C_clab_r_igh", fn<&C_clab_r_igh>(), 2},
    {"C_clab_r_fmt", fn<&C_clab_r_fmt>(), 3},
    {"C_clab_r_ascii", fn<&C_clab_r_ascii>(), 1},
    {"C_clab_r_prim", fn<&C_clab_r_prim>(), 4},
    {"C_clab_r_hsacc", fn<&C_clab_r_hsacc>(), 6},
    {"C_clab_r_svk", fn<&C_clab_r_svk>(), 2},
    {"C_clab_r_pwrs", fn<&C_clab_r_pwrs>(), 6},
    {"C_clab_r_pwrd", fn<&C_clab_r_pwrd>(), 2},
    {"C_clab_r_pwrp", fn<&C_clab_r_pwrp>(), 3},
    {"C_clab_r_apcie", fn<&C_clab_r_apcie>(), 4},
    {"C_clab_r_spcure", fn<&C_clab_r_spcure>(), 6},
    {"C_clab_r_lnkpr", fn<&C_clab_r_lnkpr>(), 5},
    {"C_clab_r_prohs", fn<&C_clab_r_prohs>(), 11},
    {"C_clab_r_wdqx", fn<&C_clab_r_wdqx>(), 3},
    {"C_clab_r_bayout", fn<&C_clab_r_bayout>(), 7},
    {"C_clab_r_bufexp", fn<&C_clab_r_bufexp>(), 6},
    {"C_clab_r_kerexp", fn<&C_clab_r_kerexp>(), 6},
    {"C_clab_r_lerprc", fn<&C_clab_r_lerprc>(), 3},
    {"C_clab_r_bym2sc", fn<&C_clab_r_bym2sc>(), 4},
    {"C_clab_r_ecolrg", fn<&C_clab_r_ecolrg>(), 7},
    {"C_clab_r_spbias", fn<&C_clab_r_spbias>(), 3},
    {"C_clab_r_spsir", fn<&C_clab_r_spsir>(), 8},
    {"C_clab_r_sfrail", fn<&C_clab_r_sfrail>(), 7},
    {"C_clab_r_sfcure", fn<&C_clab_r_sfcure>(), 7},
    {"C_clab_r_evfilt", fn<&C_clab_r_evfilt>(), 7},
    {"C_clab_r_vgemp", fn<&C_clab_r_vgemp>(), 5},
    {"C_clab_r_vgfit", fn<&C_clab_r_vgfit>(), 4},
    {"C_clab_r_okrge", fn<&C_clab_r_okrge>(), 5},
    {"C_clab_r_bkrig", fn<&C_clab_r_bkrig>(), 8},
    {"C_clab_r_krigcv", fn<&C_clab_r_krigcv>(), 6},
    {"C_clab_r_trsurf", fn<&C_clab_r_trsurf>(), 4},
    {"C_clab_r_sglmx", fn<&C_clab_r_sglmx>(), 4},
    {"C_clab_r_smixm", fn<&C_clab_r_smixm>(), 4},
    {"C_clab_r_lfspec", fn<&C_clab_r_lfspec>(), 1},
    {"C_clab_r_spsamp", fn<&C_clab_r_spsamp>(), 4},
    {"C_clab_r_rptvic", fn<&C_clab_r_rptvic>(), 4},
    {"C_clab_r_lnroc", fn<&C_clab_r_lnroc>(), 9},
    {"C_clab_r_lnsize", fn<&C_clab_r_lnsize>(), 4},
    {"C_clab_r_ibench", fn<&C_clab_r_ibench>(), 7},
    {"C_clab_r_thrtst", fn<&C_clab_r_thrtst>(), 14},
    {"C_clab_r_sarlr", fn<&C_clab_r_sarlr>(), 5},
    {"C_clab_r_sarsw", fn<&C_clab_r_sarsw>(), 2},
    {"C_clab_r_sarpod", fn<&C_clab_r_sarpod>(), 3},
    {"C_clab_r_sarpsw", fn<&C_clab_r_sarpsw>(), 7},
    {"C_clab_r_sarscs", fn<&C_clab_r_sarscs>(), 2},
    {"C_clab_r_sarmap", fn<&C_clab_r_sarmap>(), 8},
    {"C_clab_r_sarcrc", fn<&C_clab_r_sarcrc>(), 3},
    {"C_clab_r_saropt", fn<&C_clab_r_saropt>(), 5},
    {"C_clab_r_sarsq", fn<&C_clab_r_sarsq>(), 3},
    {"C_clab_r_sarfig", fn<&C_clab_r_sarfig>(), 0},
    {"C_clab_r_pgsig", fn<&C_clab_r_pgsig>(), 3},
    {"C_clab_r_briggs", fn<&C_clab_r_briggs>(), 9},
    {"C_clab_r_gplume", fn<&C_clab_r_gplume>(), 8},
    {"C_clab_r_gspuff", fn<&C_clab_r_gspuff>(), 9},
    {"C_clab_r_advdif", fn<&C_clab_r_advdif>(), 10},
    {"C_clab_r_laggrid", fn<&C_clab_r_laggrid>(), 3},
    {"C_clab_r_lagr", fn<&C_clab_r_lagr>(), 11},
    {"C_clab_r_lagstep", fn<&C_clab_r_lagstep>(), 9},
    {"C_clab_r_crf", fn<&C_clab_r_crf>(), 4},
    {"C_clab_r_levin", fn<&C_clab_r_levin>(), 2},
    {"C_clab_r_hif", fn<&C_clab_r_hif>(), 4},
    {"C_clab_r_concidx", fn<&C_clab_r_concidx>(), 2},
    {"C_clab_r_crc32", fn<&C_clab_r_crc32>(), 2},
    {"C_clab_r_fesep", fn<&C_clab_r_fesep>(), 5},
    {"C_clab_r_ppmlfe", fn<&C_clab_r_ppmlfe>(), 16},
    {"C_clab_r_ivfe", fn<&C_clab_r_ivfe>(), 16},
    {"C_clab_r_ppmlcf", fn<&C_clab_r_ppmlcf>(), 15},
    {"C_clab_r_demean", fn<&C_clab_r_demean>(), 8},
    {"C_clab_r_expgmm", fn<&C_clab_r_expgmm>(), 14},
    {"C_clab_r_ivoid", fn<&C_clab_r_ivoid>(), 15},
    {"C_clab_r_rdck", fn<&C_clab_r_rdck>(), 1},
    {"C_clab_r_ikbw", fn<&C_clab_r_ikbw>(), 4},
    {"C_clab_r_mccr", fn<&C_clab_r_mccr>(), 4},
    {"C_clab_r_windbin", fn<&C_clab_r_windbin>(), 7},
    {"C_clab_r_dwnexp", fn<&C_clab_r_dwnexp>(), 14},
    {"C_clab_r_tailwd", fn<&C_clab_r_tailwd>(), 6},
    {"C_clab_r_nbfit", fn<&C_clab_r_nbfit>(), 7},
    {"C_clab_r_idwint", fn<&C_clab_r_idwint>(), 11},
    {"C_clab_r_gtwreg", fn<&C_clab_r_gtwreg>(), 8},
    {"C_clab_r_asrule", fn<&C_clab_r_asrule>(), 7},
    {"C_clab_r_dshift", fn<&C_clab_r_dshift>(), 3},
    {"C_clab_r_updwn", fn<&C_clab_r_updwn>(), 10},
    {"C_clab_r_vcx", fn<&C_clab_r_vcx>(), 10},
    {"C_clab_r_runlen", fn<&C_clab_r_runlen>(), 4},
    {"C_clab_r_momupd", fn<&C_clab_r_momupd>(), 2},
    {"C_clab_r_mommrg", fn<&C_clab_r_mommrg>(), 2},
    {"C_clab_r_momstat", fn<&C_clab_r_momstat>(), 1},
    {"C_clab_r_rsv", fn<&C_clab_r_rsv>(), 3},
    {"C_clab_r_hlladd", fn<&C_clab_r_hlladd>(), 3},
    {"C_clab_r_hllcnt", fn<&C_clab_r_hllcnt>(), 1},
    {"C_clab_r_citab", fn<&C_clab_r_citab>(), 0},
    {"C_clab_r_cimatch", fn<&C_clab_r_cimatch>(), 1},
    {"C_clab_r_tzterr", fn<&C_clab_r_tzterr>(), 1},
    {"C_clab_r_proxy", fn<&C_clab_r_proxy>(), 1},
    {"C_clab_r_fpe", fn<&C_clab_r_fpe>(), 7},
    {"C_clab_r_fpeqv", fn<&C_clab_r_fpeqv>(), 1},
    {"C_clab_r_fploc", fn<&C_clab_r_fploc>(), 4},
    {"C_clab_r_cizone", fn<&C_clab_r_cizone>(), 1},
    {"C_clab_r_rapldom", fn<&C_clab_r_rapldom>(), 1},
    {"C_clab_r_raplrd", fn<&C_clab_r_raplrd>(), 1},
    {"C_clab_r_cpuacc", fn<&C_clab_r_cpuacc>(), 1},
    {"C_clab_r_memgb", fn<&C_clab_r_memgb>(), 0},
    {"C_clab_r_fpnote", fn<&C_clab_r_fpnote>(), 2},
    {"C_clab_r_plr", fn<&C_clab_r_plr>(), 6},
    {"C_clab_r_bootci", fn<&C_clab_r_bootci>(), 2},
    {"C_clab_r_hbin", fn<&C_clab_r_hbin>(), 2},
    {"C_clab_r_freq", fn<&C_clab_r_freq>(), 7},
    {"C_clab_r_pcor", fn<&C_clab_r_pcor>(), 3},
    {"C_clab_r_rowany", fn<&C_clab_r_rowany>(), 1},
    {"C_clab_r_nrun", fn<&C_clab_r_nrun>(), 2},
    {"C_clab_r_msum", fn<&C_clab_r_msum>(), 1},
    {"C_clab_r_pstat", fn<&C_clab_r_pstat>(), 1},
    {"C_clab_r_q7", fn<&C_clab_r_q7>(), 2},
    {"C_clab_r_maha", fn<&C_clab_r_maha>(), 3},
    {"C_clab_r_mcar", fn<&C_clab_r_mcar>(), 3},
    {"C_clab_r_mgrid", fn<&C_clab_r_mgrid>(), 2},
    {"C_clab_r_dups", fn<&C_clab_r_dups>(), 2},
    {"C_clab_r_cochq", fn<&C_clab_r_cochq>(), 2},
    {"C_clab_r_dfit", fn<&C_clab_r_dfit>(), 2},
    {"C_clab_r_dfgof", fn<&C_clab_r_dfgof>(), 3},
    {"C_clab_r_dfboot", fn<&C_clab_r_dfboot>(), 4},
    {"C_clab_r_dfhill", fn<&C_clab_r_dfhill>(), 3},
    {"C_clab_r_svymx", fn<&C_clab_r_svymx>(), 10},
    {"C_clab_r_crsev", fn<&C_clab_r_crsev>(), 9},
    {"C_clab_r_zcomp", fn<&C_clab_r_zcomp>(), 4},
    {"C_clab_r_cconc", fn<&C_clab_r_cconc>(), 7},
    {"C_clab_r_hilltl", fn<&C_clab_r_hilltl>(), 5},
    {"C_clab_r_hzta", fn<&C_clab_r_hzta>(), 2},
    {"C_clab_r_crcap", fn<&C_clab_r_crcap>(), 3},
    {"C_clab_r_dkfig", fn<&C_clab_r_dkfig>(), 10},
    {"C_clab_r_bandx", fn<&C_clab_r_bandx>(), 8},
    {"C_clab_r_pbbnd", fn<&C_clab_r_pbbnd>(), 12},
    {"C_clab_r_yoypx", fn<&C_clab_r_yoypx>(), 11},
    {"C_clab_r_jdgiv", fn<&C_clab_r_jdgiv>(), 6},
    {"C_clab_r_jmono", fn<&C_clab_r_jmono>(), 10},
    {"C_clab_r_mmshr", fn<&C_clab_r_mmshr>(), 3},
    {"C_clab_r_agemix", fn<&C_clab_r_agemix>(), 6},
    {"C_clab_r_cumrsk", fn<&C_clab_r_cumrsk>(), 6},
    {"C_clab_r_jpipe", fn<&C_clab_r_jpipe>(), 5},
    {"C_clab_r_mndla", fn<&C_clab_r_mndla>(), 12},
    {"C_clab_r_stkflw", fn<&C_clab_r_stkflw>(), 11},
    {"C_clab_r_lostay", fn<&C_clab_r_lostay>(), 6},
    {"C_clab_r_bklog", fn<&C_clab_r_bklog>(), 4},
    {"C_clab_r_dktflw", fn<&C_clab_r_dktflw>(), 12},
    {"C_clab_r_lftab", fn<&C_clab_r_lftab>(), 6},
    {"C_clab_r_trnmat", fn<&C_clab_r_trnmat>(), 3},
    {"C_clab_r_intrt", fn<&C_clab_r_intrt>(), 1},
    {"C_clab_r_bndxw", fn<&C_clab_r_bndxw>(), 5},
    {"C_clab_r_synexp", fn<&C_clab_r_synexp>(), 4},
    {"C_clab_r_detrsp", fn<&C_clab_r_detrsp>(), 4},
    {"C_clab_r_dtrid", fn<&C_clab_r_dtrid>(), 6},
    {"C_clab_r_inspg", fn<&C_clab_r_inspg>(), 4},
    {"C_clab_r_insphs", fn<&C_clab_r_insphs>(), 5},
    {"C_clab_r_colarr", fn<&C_clab_r_colarr>(), 3},
    {"C_clab_r_clrilr", fn<&C_clab_r_clrilr>(), 3},
    {"C_clab_r_preslr", fn<&C_clab_r_preslr>(), 5},
    {"C_clab_r_tpdmc", fn<&C_clab_r_tpdmc>(), 7},
    {"C_clab_r_tpbhm", fn<&C_clab_r_tpbhm>(), 6},
    {"C_clab_r_dirsim", fn<&C_clab_r_dirsim>(), 4},
    {"C_clab_r_txcal", fn<&C_clab_r_txcal>(), 5},
    {"C_clab_r_lodimp", fn<&C_clab_r_lodimp>(), 3},
    {"C_clab_r_pmrrat", fn<&C_clab_r_pmrrat>(), 3},
    {"C_clab_r_etohlr", fn<&C_clab_r_etohlr>(), 5},
    {"C_clab_r_mtxrel", fn<&C_clab_r_mtxrel>(), 8},
    {"C_clab_r_crchyp", fn<&C_clab_r_crchyp>(), 3},
    {"C_clab_r_hitscr", fn<&C_clab_r_hitscr>(), 2},
    {"C_clab_r_rtmfit", fn<&C_clab_r_rtmfit>(), 6},
    {"C_clab_r_crfpar", fn<&C_clab_r_crfpar>(), 2},
    {"C_clab_r_pbarea", fn<&C_clab_r_pbarea>(), 6},
    {"C_clab_r_rocpts", fn<&C_clab_r_rocpts>(), 2},
    {"C_clab_r_spark", fn<&C_clab_r_spark>(), 3},
    {"C_clab_r_cnames", fn<&C_clab_r_cnames>(), 3},
    {"C_clab_r_mgridl", fn<&C_clab_r_mgridl>(), 5},
    {"C_clab_r_mpat", fn<&C_clab_r_mpat>(), 3},
    {"C_clab_r_prof", fn<&C_clab_r_prof>(), 3},
    {"C_clab_r_mhtab", fn<&C_clab_r_mhtab>(), 3},
    {"C_clab_r_mruns", fn<&C_clab_r_mruns>(), 3},
    {"C_clab_r_topcor", fn<&C_clab_r_topcor>(), 4},
    {"C_clab_r_dtsel", fn<&C_clab_r_dtsel>(), 7},
    {"C_clab_r_mcar2", fn<&C_clab_r_mcar2>(), 3},
    {"C_clab_r_bz2end", fn<&C_clab_r_bz2end>(), 1},
    {"C_clab_r_gzend", fn<&C_clab_r_gzend>(), 2},
    {"C_clab_r_txtenc", fn<&C_clab_r_txtenc>(), 1},
    {"C_clab_r_rdfit", fn<&C_clab_r_rdfit>(), 8},
    {nullptr, nullptr, 0}};

void R_init_crimlab(DllInfo *dll) {
  clab_onintr(&rintr, -1.0);
  R_registerRoutines(dll, nullptr, calls, nullptr, nullptr);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}

}  // extern "C"
