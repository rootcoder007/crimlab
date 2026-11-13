// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the dark-figure family (crcap, dkfig): each checks
// its arguments, calls the kernels and builds the R result. Nothing calls
// back into the R evaluator.
#include <cmath>
#include <cstring>

#include "bcol.h"
#include "rchk.h"

using clab::fail;

namespace {

// x in the order of the names nm (k of them) when every name is one of
// them once, otherwise as given (need: names required).
void order(SEXP x, SEXP xv, const char *const *nm, int k, bool need, double *out) {
  SEXP xn = Rf_getAttrib(x, R_NamesSymbol);
  if (Rf_isNull(xn)) {
    if (need) fail("`x` must name its seven cells 100, 010, 001, 110, 101, 011 and 111.");
    for (int j = 0; j < k; ++j) out[j] = REAL(xv)[j];
    return;
  }
  for (int j = 0; j < k; ++j) {
    int hit = -1;
    for (int i = 0; i < k; ++i) {
      if (STRING_ELT(xn, i) != NA_STRING && std::strcmp(bj::str(xn, i), nm[j]) == 0) {
        if (hit >= 0) fail("`x` names %s twice.", nm[j]);
        hit = i;
      }
    }
    if (hit < 0) {
      if (need) fail("`x` must name its seven cells 100, 010, 001, 110, 101, 011 and 111.");
      fail("`x` must be named n1, n2 and m, or not named at all.");
    }
    out[j] = REAL(xv)[hit];
  }
}

}  // namespace

extern "C" {

SEXP C_clab_r_crcap(SEXP x, SEXP kappa, SEXP it) {
  clab::who = "crcap";
  int np = 0;
  SEXP xv = PROTECT(bj::nums(x, "x", -1, 0.0));
  ++np;
  const int k = bj::len(xv, "x");
  if (k != 3 && k != 7) fail("`x` must be c(n1, n2, m) for two lists or the seven cells of three.");
  const double kp = bj::inr(kappa, "kappa", 1.0, R_PosInf);
  double c[7];
  static const char *const two[] = {"n1", "n2", "m"};
  static const char *const three[] = {"100", "010", "001", "110", "101", "011", "111"};
  order(x, xv, k == 3 ? two : three, k, k == 7, c);
  int ni = 0;
  SEXP iv = R_NilValue;
  if (!Rf_isNull(it)) {
    if (k == 3) fail("`int` goes with three lists.");
    iv = PROTECT(bj::nums(it, "int"));
    ++np;
    ni = bj::len(iv, "int");
  }
  double st[6], pw[18];
  SEXP sv = PROTECT(Rf_allocVector(REALSXP, 2 * static_cast<R_xlen_t>(ni)));
  ++np;
  const int rc = clab_crcap(c, k, kp, Rf_isNull(iv) ? nullptr : REAL(iv), ni, st, pw, REAL(sv));
  if (rc != CLAB_OK) {
    fail(k == 3 ? "`x` must have an overlap m above 0 and no larger than either list."
                : "the three cells of two lists only (110, 101, 011) must be above 0.");
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  ++np;
  if (k == 3) {
    static const char *const en[] = {"petersen", "lower", "upper", "chapman", "se", "floor"};
    SET_VECTOR_ELT(out, 0, Rf_mkString("two"));
    SET_VECTOR_ELT(out, 1, clab::nvec(6, en, st));
  } else {
    static const char *const en[] = {"observed", "missing", "total"};
    SET_VECTOR_ELT(out, 0, Rf_mkString("three"));
    SET_VECTOR_ELT(out, 1, clab::nvec(3, en, st));
    bj::cols p(7, 3);
    SEXP ls = PROTECT(Rf_allocVector(STRSXP, 3));
    np += 2;
    static const char *const pl[] = {"1-2", "1-3", "2-3"};
    for (int j = 0; j < 3; ++j) SET_STRING_ELT(ls, j, Rf_mkChar(pl[j]));
    p.put("lists", ls);
    static const char *const pn[] = {"n1", "n2", "m", "floor", "petersen", "chapman"};
    for (int q = 0; q < 6; ++q) {
      double *d = p.dbl(pn[q]);
      for (int j = 0; j < 3; ++j) d[j] = pw[3 * q + j];
    }
    SET_VECTOR_ELT(out, 2, p.done());
    if (ni > 0) {
      bj::cols s(3, ni);
      ++np;
      s.put("int", iv);
      double *ms = s.dbl("missing"), *tt = s.dbl("total");
      for (int j = 0; j < ni; ++j) {
        ms[j] = REAL(sv)[j];
        tt[j] = REAL(sv)[ni + j];
      }
      SET_VECTOR_ELT(out, 3, s.done());
    }
  }
  static const char *const nm[] = {"type", "est", "pairs", "sens"};
  UNPROTECT(np);
  return bj::named(out, nm);
}

SEXP C_clab_r_dkfig(SEXP v, SEXP r, SEXP amax, SEXP bmax, SEXP thr, SEXP va, SEXP ra, SEXP rb,
                    SEXP lam, SEXP k) {
  clab::who = "dkfig";
  int np = 0;
  if (Rf_isNull(v) && Rf_isNull(va) && Rf_isNull(k)) {
    fail("give `v` and `r` (survey and recorded rates), `va`, `ra` and `rb`, or `k`.");
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  ++np;
  if (!Rf_isNull(v)) {
    SEXP vv = PROTECT(bj::nums(v, "v", -1, 0.0));
    ++np;
    const int n = bj::len(vv, "v");
    if (Rf_isNull(r)) fail("`r` must give the recorded rate with `v`.");
    SEXP r0 = PROTECT(bj::nums(r, "r", -1, 0.0));
    ++np;
    if (XLENGTH(r0) != 1 && XLENGTH(r0) != n)
      fail("`r` must have length 1 or the length of `v`, %d.", n);
    for (int i = 0; i < n; ++i) {
      if (REAL(vv)[i] > 1) fail("`v` must be rates from 0 to 1.");
      if (REAL(r0)[XLENGTH(r0) == 1 ? 0 : i] > 1) fail("`r` must be rates from 0 to 1.");
    }
    const double a = bj::inr(amax, "amax", 0.0, 1.0), b = bj::inr(bmax, "bmax", 0.0, 1.0);
    if (a + b >= 1) fail("`amax + bmax` must be below 1.");
    double th = NAN;
    if (!Rf_isNull(thr)) {
      th = clab::one(thr, "thr", 0.0, R_PosInf);
      for (int i = 0; i < n; ++i) {
        if (!(REAL(vv)[i] > 0 && REAL(vv)[i] < th)) {
          fail("`thr` must lie above every `v`, and each `v` above 0.");
        }
      }
    }
    bj::cols f(9, n);
    ++np;
    f.put("v", vv);
    double *rr = f.dbl("r"), *lo = f.dbl("lo"), *hi = f.dbl("hi"), *dlo = f.dbl("dlo"),
           *dhi = f.dbl("dhi"), *ru = f.dbl("ratio");
    for (int i = 0; i < n; ++i) rr[i] = REAL(r0)[XLENGTH(r0) == 1 ? 0 : i];
    SEXP bk = PROTECT(Rf_allocVector(REALSXP, n));
    ++np;
    clab_dkfig(REAL(vv), rr, n, a, b, th, lo, hi, REAL(bk));
    for (int i = 0; i < n; ++i) {
      dlo[i] = lo[i] - rr[i];
      dhi[i] = hi[i] - rr[i];
      ru[i] = rr[i] > 0 ? hi[i] / rr[i] : NA_REAL;
    }
    if (!Rf_isNull(thr)) f.put("brk", bk);
    SET_VECTOR_ELT(out, 0, f.done());
  }
  if (!Rf_isNull(va)) {
    if (Rf_isNull(ra) || Rf_isNull(rb)) fail("`ra` and `rb` must come with `va`.");
    const double a = clab::one(va, "va", 0.0, R_PosInf);
    const double b = clab::one(ra, "ra", 0.0, R_PosInf);
    const double c = bj::inr(rb, "rb", 0.0, R_PosInf);
    const double l = bj::inr(lam, "lam", 0.0, 1.0);
    if (l == 0) fail("`lam` must be above 0.");
    if (b > a) fail("`ra` must not exceed `va`: a record is an offence.");
    double o[7];
    clab_dkcmp(a, b, c, l, o);
    static const char *const cn[] = {"rate", "rateup", "totlo", "tothi", "dklo", "dkhi", "width"};
    SET_VECTOR_ELT(out, 1, clab::nvec(7, cn, o));
  }
  if (!Rf_isNull(k)) {
    SEXP kv = PROTECT(bj::nums(k, "k", -1, 1.0));
    ++np;
    double o[6];
    if (clab_dkhier(REAL(kv), bj::len(kv, "k"), o) != CLAB_OK) {
      fail("`k` must be whole numbers of offences from 1.");
    }
    static const char *const hn[] = {"incid", "offen", "ratio", "extra", "lower", "upper"};
    SET_VECTOR_ELT(out, 2, clab::nvec(6, hn, o));
  }
  static const char *const nm[] = {"bnd", "cmp", "hier"};
  UNPROTECT(np);
  return bj::named(out, nm);
}

}  // extern "C"
