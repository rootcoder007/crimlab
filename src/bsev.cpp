// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the severity and composite indices (crsev,
// zcomp): each checks its arguments, calls the kernel and builds the R
// result. Nothing calls back into the R evaluator.
#include <algorithm>
#include <vector>

#include "bcol.h"
#include "rchk.h"

using clab::fail;

namespace {

// The cells of the rows: the pairs (unit, time) of codes that occur, in
// order of unit then time, and each row's cell.
int cells(const int *u, const int *t, int nt, int n, int *key, int *cu, int *ct) {
  std::vector<long long> k(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i)
    k[static_cast<std::size_t>(i)] = static_cast<long long>(u[i]) * nt + t[i];
  std::vector<long long> s(k);
  std::sort(s.begin(), s.end());
  s.erase(std::unique(s.begin(), s.end()), s.end());
  for (int i = 0; i < n; ++i) {
    key[i] = static_cast<int>(std::lower_bound(s.begin(), s.end(), k[static_cast<std::size_t>(i)]) -
                              s.begin());
  }
  for (std::size_t j = 0; j < s.size(); ++j) {
    cu[j] = static_cast<int>(s[j] / nt);
    ct[j] = static_cast<int>(s[j] % nt);
  }
  return static_cast<int>(s.size());
}

// The labels lev[code[j]] of nk cells.
SEXP pick(SEXP lev, const int *code, int nk) {
  SEXP s = PROTECT(Rf_allocVector(STRSXP, nk));
  for (int j = 0; j < nk; ++j) SET_STRING_ELT(s, j, STRING_ELT(lev, code[j]));
  UNPROTECT(1);
  return s;
}

}  // namespace

extern "C" {

SEXP C_clab_r_crsev(SEXP x, SEXP off, SEXP w, SEXP unit, SEXP time, SEXP pop, SEXP per, SEXP base,
                    SEXP sent) {
  clab::who = "crsev";
  int np = 0;
  SEXP xv = PROTECT(bj::nums(x, "x", -1, 0.0));
  ++np;
  const int n = bj::len(xv, "x");
  csp::lab o = bj::labs(off, "off", n);
  np += 2;
  SEXP wv = PROTECT(Rf_allocVector(REALSXP, o.g));
  ++np;
  bj::perlev(w, "w", o.lev, o.g, REAL(wv));
  if (!Rf_isNull(sent)) {
    SEXP sv = PROTECT(Rf_allocVector(REALSXP, o.g));
    ++np;
    bj::perlev(sent, "sent", o.lev, o.g, REAL(sv));
    for (int j = 0; j < o.g; ++j) REAL(wv)[j] *= REAL(sv)[j];
  }
  SEXP zero = PROTECT(Rf_allocVector(INTSXP, n));
  ++np;
  for (int i = 0; i < n; ++i) INTEGER(zero)[i] = 0;
  const int *uc = INTEGER(zero), *tc = INTEGER(zero);
  int nt = 1;
  SEXP ulev = R_NilValue, tlev = R_NilValue;
  if (!Rf_isNull(unit)) {
    csp::lab u = bj::labs(unit, "unit", n);
    np += 2;
    uc = INTEGER(u.code);
    ulev = u.lev;
  }
  if (!Rf_isNull(time)) {
    csp::lab t = bj::labs(time, "time", n);
    np += 2;
    tc = INTEGER(t.code);
    tlev = t.lev;
    nt = t.g;
  }
  SEXP kv = PROTECT(Rf_allocVector(INTSXP, 3 * static_cast<R_xlen_t>(n)));
  ++np;
  int *key = INTEGER(kv), *cu = key + n, *ct = cu + n;
  const int nk = cells(uc, tc, nt, n, key, cu, ct);
  // one population per cell, read from its rows
  SEXP pv = R_NilValue;
  if (!Rf_isNull(pop)) {
    SEXP pr = PROTECT(bj::nums(pop, "pop", n, 0.0, true));
    pv = PROTECT(Rf_allocVector(REALSXP, nk));
    np += 2;
    for (int j = 0; j < nk; ++j) REAL(pv)[j] = -1.0;
    for (int i = 0; i < n; ++i) {
      double &c = REAL(pv)[key[i]];
      if (c >= 0 && c != REAL(pr)[i])
        fail("`pop` must be the same on every row of a unit and period.");
      c = REAL(pr)[i];
    }
  }
  const double pe = clab::one(per, "per", 0.0, R_PosInf);
  const int b = Rf_isNull(base) ? -1 : clab::whole(base, "base", 1) - 1;
  if (b >= nk) fail("`base` must be a row of the result, from 1 to %d.", nk);
  bj::cols f(9, nk);
  if (!Rf_isNull(ulev)) f.put("unit", pick(ulev, cu, nk));
  if (!Rf_isNull(tlev)) f.put("time", pick(tlev, ct, nk));
  double *cnt = f.dbl("count"), *sev = f.dbl("sev");
  SEXP ws = PROTECT(Rf_allocVector(REALSXP, 3 * static_cast<R_xlen_t>(nk)));
  ++np;
  double *rate = REAL(ws), *srate = rate + nk, *idx = srate + nk;
  bj::kc(clab_crsev(REAL(xv), INTEGER(o.code), key, n, o.g, nk, REAL(wv),
                    Rf_isNull(pv) ? nullptr : REAL(pv), pe, b, cnt, sev, rate, srate, idx),
         "the kernel refused the input.");
  if (!Rf_isNull(pv)) {
    f.put("pop", pv);
    double *r = f.dbl("rate"), *s = f.dbl("srate");
    std::copy(rate, rate + nk, r);
    std::copy(srate, srate + nk, s);
  }
  if (b >= 0) {
    double *ix = f.dbl("idx");
    for (int j = 0; j < nk; ++j) ix[j] = bj::na(idx[j]);
  }
  bj::cols g(2, o.g);
  g.put("off", o.lev);
  g.put("w", wv);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  ++np;
  SET_VECTOR_ELT(out, 0, f.done());
  SET_VECTOR_ELT(out, 1, g.done());
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(pe));
  SET_VECTOR_ELT(out, 3, Rf_ScalarInteger(b < 0 ? NA_INTEGER : b + 1));
  static const char *const nm[] = {"tab", "w", "per", "base"};
  UNPROTECT(np + 2);
  return bj::named(out, nm);
}

SEXP C_clab_r_zcomp(SEXP x, SEXP unit, SEXP cat, SEXP w) {
  clab::who = "zcomp";
  int np = 0;
  SEXP xv = PROTECT(bj::nums(x, "x", -1, 0.0));
  ++np;
  const int n = bj::len(xv, "x");
  csp::lab u = bj::labs(unit, "unit", n);
  csp::lab c = bj::labs(cat, "cat", n);
  np += 4;
  if (u.g < 2) fail("`unit` must hold at least two units.");
  const long long cells = static_cast<long long>(u.g) * c.g;
  if (cells > 100000000LL) fail("`unit` and `cat` make more than 1e8 cells.");
  // counts per unit and category, 0 where no row gives one
  SEXP xm = PROTECT(Rf_allocMatrix(REALSXP, u.g, c.g));
  SEXP zm = PROTECT(Rf_allocMatrix(REALSXP, u.g, c.g));
  np += 2;
  std::fill(REAL(xm), REAL(xm) + cells, 0.0);
  for (int i = 0; i < n; ++i) {
    REAL(xm)[INTEGER(u.code)[i] + static_cast<long long>(u.g) * INTEGER(c.code)[i]] += REAL(xv)[i];
  }
  SEXP wv = PROTECT(Rf_allocVector(REALSXP, c.g));
  ++np;
  if (Rf_isNull(w)) {
    for (int j = 0; j < c.g; ++j) REAL(wv)[j] = 1.0;
  } else {
    bj::perlev(w, "w", c.lev, c.g, REAL(wv));
  }
  bj::cols f(3, u.g);
  f.put("unit", u.lev);
  double *comp = f.dbl("comp");
  int *rk = f.itg("rank");
  bj::cols g(4, c.g);
  g.put("cat", c.lev);
  double *mu = g.dbl("mean"), *sd = g.dbl("sd");
  g.put("w", wv);
  clab_zcomp(REAL(xm), u.g, c.g, REAL(wv), REAL(zm), comp, mu, sd);
  // rank 1 for the highest composite, ties sharing the lowest rank
  std::vector<int> ord(static_cast<std::size_t>(u.g));
  for (int j = 0; j < u.g; ++j) ord[static_cast<std::size_t>(j)] = j;
  std::stable_sort(ord.begin(), ord.end(), [&](int a, int b) { return comp[a] > comp[b]; });
  for (int j = 0; j < u.g; ++j) {
    const int k = ord[static_cast<std::size_t>(j)];
    rk[k] = j > 0 && comp[ord[static_cast<std::size_t>(j - 1)]] == comp[k]
                ? rk[ord[static_cast<std::size_t>(j - 1)]]
                : j + 1;
  }
  SEXP dn = PROTECT(Rf_allocVector(VECSXP, 2));
  ++np;
  SET_VECTOR_ELT(dn, 0, u.lev);
  SET_VECTOR_ELT(dn, 1, c.lev);
  Rf_setAttrib(zm, R_DimNamesSymbol, dn);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  ++np;
  SET_VECTOR_ELT(out, 0, f.done());
  SET_VECTOR_ELT(out, 1, g.done());
  SET_VECTOR_ELT(out, 2, zm);
  static const char *const nm[] = {"tab", "cat", "z"};
  UNPROTECT(np + 2);
  return bj::named(out, nm);
}

}  // extern "C"
