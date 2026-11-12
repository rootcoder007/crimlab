// SPDX-License-Identifier: AGPL-3.0-or-later
//
// The R entry points of crimlab, one per exported function: each checks its
// arguments in the order and words of the R wrapper it replaces, calls the
// kernel and builds the R result (names, data frames, lists). The R
// function is its documentation and one .Call. Inputs are read in place;
// each result is allocated once, at its final size; nothing calls back
// into the R evaluator.
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstring>

#include "crimlab.h"
#include "rchk.h"

using clab::dcopy;
using clab::fail;
using clab::ilen;
using clab::kchk;
using clab::Key;
using clab::num;

namespace {

SEXP real(R_xlen_t n) { return Rf_allocVector(REALSXP, n); }

// rep(1, n)
SEXP ones(R_xlen_t n) {
  SEXP w = real(n);
  double *p = REAL(w);
  for (R_xlen_t i = 0; i < n; ++i) p[i] = 1.0;
  return w;
}

// TRUE when some p[i] > 0 (with keep[i] when keep is given). For weights
// already checked to be finite and at least 0 this is sum(p[keep]) > 0,
// exactly: a sum of such numbers is 0 only when every one is.
bool pos(const double *p, R_xlen_t n, const int *keep = nullptr) {
  for (R_xlen_t i = 0; i < n; ++i) {
    if ((keep == nullptr || keep[i]) && p[i] > 0) return true;
  }
  return false;
}

// The keys of a character vector and the positions in key order, in R's
// transient memory: string sets are matched by sorting, never pairwise.
struct Keys {
  Key *k;
  R_xlen_t *o;
  R_xlen_t n;
};

Keys keys(SEXP s) {
  const R_xlen_t n = XLENGTH(s);
  Keys r{reinterpret_cast<Key *>(R_alloc(static_cast<std::size_t>(n), sizeof(Key))),
         reinterpret_cast<R_xlen_t *>(R_alloc(static_cast<std::size_t>(n), sizeof(R_xlen_t))), n};
  for (R_xlen_t i = 0; i < n; ++i) {
    r.k[i] = clab::key(STRING_ELT(s, i));
    r.o[i] = i;
  }
  const Key *k = r.k;
  std::stable_sort(r.o, r.o + n, [k](R_xlen_t a, R_xlen_t b) {
    return clab::kcmp(k[a], k[b]) < 0;
  });
  return r;
}

// x %in% the set
bool has(const Keys &s, const Key &x) {
  const Key *k = s.k;
  const R_xlen_t *p = std::lower_bound(
      s.o, s.o + s.n, x,
      [k](R_xlen_t a, const Key &v) { return clab::kcmp(k[a], v) < 0; });
  return p != s.o + s.n && clab::kcmp(k[*p], x) == 0;
}

// the number of distinct values
R_xlen_t ndist(const Keys &s) {
  R_xlen_t d = 0;
  for (R_xlen_t i = 0; i < s.n; ++i) {
    if (i == 0 || clab::kcmp(s.k[s.o[i - 1]], s.k[s.o[i]]) != 0) ++d;
  }
  return d;
}

const char *const vodnm[] = {"group", "share_day", "share_night", "ratio",
                             "sign", "adj"};

}  // namespace

extern "C" {

// Disparity.

SEXP C_clab_r_vodtst(SEXP day, SEXP night, SEXP ref, SEXP shift) {
  SEXP nms = PROTECT(clab::grp(day));
  SEXP d = PROTECT(num(day, "day", -1, 0.0));
  const int n = ilen(XLENGTH(d), "day");
  SEXP ni = PROTECT(num(night, "night", n, 0.0, true));
  const bool has = !Rf_isNull(shift);
  SEXP sh = PROTECT(has ? num(shift, "shift", n, 0.0, true) : R_NilValue);
  const int i = clab::ref(ref, nms, n);
  if (REAL_RO(d)[i - 1] <= 0) {
    fail("the reference group has no daylight stops.");
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, has ? 6 : 5));
  SET_VECTOR_ELT(out, 0, nms);
  for (int k = 1; k < 4; ++k) SET_VECTOR_ELT(out, k, real(n));
  SET_VECTOR_ELT(out, 4, Rf_allocVector(INTSXP, n));
  if (has) SET_VECTOR_ELT(out, 5, real(n));
  kchk(clab_vod(REAL_RO(d), REAL_RO(ni), n, i - 1, has ? REAL_RO(sh) : nullptr,
                REAL(VECTOR_ELT(out, 1)), REAL(VECTOR_ELT(out, 2)),
                REAL(VECTOR_ELT(out, 3)), INTEGER(VECTOR_ELT(out, 4)),
                has ? REAL(VECTOR_ELT(out, 5)) : nullptr),
       nullptr);
  // out[["adj"]] <- ... on a data frame leaves row.names before class
  clab::frame(out, vodnm, has);
  UNPROTECT(5);
  return out;
}

SEXP C_clab_r_orrisk(SEXP orr, SEXP rate, SEXP share) {
  SEXP o = PROTECT(num(orr, "or", 1, 0.0, true));
  const double v = REAL_RO(o)[0];
  if (!std::isfinite(1 / v)) {
    fail("`or` is too small: its reciprocal is not finite.");
  }
  const double q = clab::one(rate, "rate", 0, 1);
  const double s = clab::one(share, "share", 0, 1);
  double r[4];
  kchk(clab_orrisk(v, q, s, r, r + 1), nullptr);
  r[2] = r[0] / r[1];
  r[3] = v;
  static const char *const nm[] = {"exposed", "unexposed", "rr", "or"};
  UNPROTECT(1);
  return clab::nvec(4, nm, r);
}

SEXP C_clab_r_outtst(SEXP q, SEXP t, SEXP w) {
  SEXP qv = PROTECT(num(q, "q", -1, 0.0));
  const int n = ilen(XLENGTH(qv), "q");
  const double *p = REAL_RO(qv);
  for (int i = 0; i < n; ++i) {
    if (p[i] > 1) fail("`q` must be at most 1.");
  }
  SEXP tv = PROTECT(num(t, "t", 1));
  SEXP wv = PROTECT(Rf_isNull(w) ? ones(n) : num(w, "w", n, 0.0));
  if (!pos(REAL_RO(wv), n)) fail("`w` must have a positive total.");
  double r[2];
  kchk(clab_outtst(p, REAL_RO(wv), n, REAL_RO(tv)[0], r, r + 1), nullptr);
  static const char *const nm[] = {"stop", "hit"};
  UNPROTECT(3);
  return clab::nvec(2, nm, r);
}

SEXP C_clab_r_dspbm(SEXP pop, SEXP contact, SEXP force, SEXP ref) {
  SEXP nms = PROTECT(clab::grp(pop));
  SEXP p = PROTECT(num(pop, "pop", -1, 0.0, true));
  const int n = ilen(XLENGTH(p), "pop");
  SEXP c = PROTECT(num(contact, "contact", n, 0.0, true));
  SEXP f = PROTECT(num(force, "force", n, 0.0, true));
  const int i = clab::ref(ref, nms, n);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out, 0, nms);
  for (int k = 1; k < 4; ++k) SET_VECTOR_ELT(out, k, real(n));
  kchk(clab_dspbm(REAL_RO(p), REAL_RO(c), REAL_RO(f), n, i - 1,
                  REAL(VECTOR_ELT(out, 1)), REAL(VECTOR_ELT(out, 2)),
                  REAL(VECTOR_ELT(out, 3))),
       nullptr);
  static const char *const nm[] = {"group", "cd", "fd", "rd"};
  clab::frame(out, nm);
  UNPROTECT(5);
  return out;
}

SEXP C_clab_r_dyadrt(SEXP count, SEXP pop) {
  SEXP nms = PROTECT(clab::grp(pop));
  SEXP p = PROTECT(num(pop, "pop", -1, 0.0, true));
  const int g = ilen(XLENGTH(p), "pop");
  // !is.matrix(count) || nrow(count) != g || ncol(count) != g
  bool ok = Rf_isMatrix(count);
  if (ok) {
    const int *dm = INTEGER_RO(Rf_getAttrib(count, R_DimSymbol));
    ok = dm[0] == g && dm[1] == g;
  }
  if (!ok) fail("`count` must be a %d by %d matrix.", g, g);
  // .num(as.vector(count), ...): as.vector() drops the class of a numeric
  // matrix (a Date or difftime one is read as its numbers) but turns a
  // factor into its labels, which are not numbers.
  const int ty = TYPEOF(count);
  const bool data = (ty == INTSXP || ty == REALSXP) && !Rf_isFactor(count);
  if (!data) fail("`count` must be finite numbers with no missing values.");
  SEXP cv = PROTECT(Rf_isObject(count) ? (ty == REALSXP ? dcopy(count) : clab::dbl(count)) : count);
  SEXP x = PROTECT(num(cv, "count", -1, 0.0));
  const R_xlen_t m = XLENGTH(x);
  if (!pos(REAL_RO(x), m)) fail("`count` holds no offences.");
  // off <- rep(seq_len(g), times = g); vic <- rep(seq_len(g), each = g)
  int *off = reinterpret_cast<int *>(R_alloc(static_cast<std::size_t>(m), sizeof(int)));
  int *vic = reinterpret_cast<int *>(R_alloc(static_cast<std::size_t>(m), sizeof(int)));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 7));
  SEXP so = Rf_allocVector(STRSXP, m);
  SET_VECTOR_ELT(out, 0, so);
  SEXP sv = Rf_allocVector(STRSXP, m);
  SET_VECTOR_ELT(out, 1, sv);
  for (R_xlen_t k = 0; k < m; ++k) {
    off[k] = static_cast<int>(k % g);
    vic[k] = static_cast<int>(k / g);
    SET_STRING_ELT(so, k, STRING_ELT(nms, off[k]));
    SET_STRING_ELT(sv, k, STRING_ELT(nms, vic[k]));
  }
  SET_VECTOR_ELT(out, 2, dcopy(x));
  for (int k = 3; k < 7; ++k) SET_VECTOR_ELT(out, k, real(m));
  kchk(clab_dyadrt(REAL_RO(x), off, vic, ilen(m, "count"), REAL_RO(p), g,
                   REAL(VECTOR_ELT(out, 3)), REAL(VECTOR_ELT(out, 4)),
                   REAL(VECTOR_ELT(out, 5)), REAL(VECTOR_ELT(out, 6))),
       nullptr);
  static const char *const nm[] = {"off",       "vic",      "count", "rate_off",
                                   "rate_pair", "null_off", "ratio"};
  clab::frame(out, nm);
  UNPROTECT(5);
  return out;
}

}  // extern "C"
