// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
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
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <initializer_list>
#include <utility>

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

// A string set: the keys of a character vector in an open-addressing hash
// table, in R's transient memory, so matching is one pass (as R's match()
// hashes), never pairwise. first[i] marks the first appearance of a value;
// nd counts the distinct values.
struct Keys {
  Key *k;
  R_xlen_t *slot;
  char *first;
  R_xlen_t cap, n, nd;
};

std::uint64_t hash(const Key &x) {
  // FNV-1a over the UTF-8 bytes; NA and "bytes" strings hash apart
  std::uint64_t h = x.s == nullptr ? 1u : (x.b ? 2u : 3u);
  for (const unsigned char *p = reinterpret_cast<const unsigned char *>(x.s); p && *p; ++p) {
    h = (h ^ *p) * 1099511628211u;
  }
  return h;
}

// the slot that holds x, or the empty slot where x would go
R_xlen_t probe(const Keys &s, const Key &x) {
  R_xlen_t i = static_cast<R_xlen_t>(hash(x) & static_cast<std::uint64_t>(s.cap - 1));
  while (s.slot[i] >= 0 && clab::kcmp(s.k[s.slot[i]], x) != 0) i = (i + 1) & (s.cap - 1);
  return i;
}

Keys keys(SEXP v) {
  const R_xlen_t n = XLENGTH(v);
  Keys s{nullptr, nullptr, nullptr, 2, n, 0};
  while (s.cap < 2 * n) s.cap *= 2;
  s.k = reinterpret_cast<Key *>(R_alloc(static_cast<std::size_t>(n), sizeof(Key)));
  s.first = R_alloc(static_cast<std::size_t>(n), 1);
  s.slot = reinterpret_cast<R_xlen_t *>(R_alloc(static_cast<std::size_t>(s.cap), sizeof(R_xlen_t)));
  for (R_xlen_t i = 0; i < s.cap; ++i) s.slot[i] = -1;
  for (R_xlen_t i = 0; i < n; ++i) {
    s.k[i] = clab::key(STRING_ELT(v, i));
    const R_xlen_t j = probe(s, s.k[i]);
    s.first[i] = s.slot[j] < 0;
    if (s.first[i]) {
      s.slot[j] = i;
      ++s.nd;
    }
  }
  return s;
}

// x %in% the set
bool has(const Keys &s, const Key &x) { return s.n > 0 && s.slot[probe(s, x)] >= 0; }

// Group arguments.
//
// One rule for the arguments of a call that hold one value per group (or
// per unit, per locus, per arm), applied by every entry point that takes
// two or more of them. Names are the identity of a group:
//   - all named, the same set of names: every later argument is put in
//     the order of the first, by name;
//   - all named, different sets: an error that lists the names missing
//     and extra; nothing is dropped, padded or matched by position;
//   - some named, some not: positions, as before; the groups are named by
//     the first named argument;
//   - none named: positions; the groups are "1", "2", ...;
//   - a name given twice in one argument: an error.
// Only arguments of the length of the first take part (one value to be
// recycled, or a wrong length that the checks after refuse, keep their
// place). A matrix is matched by its rows, and only with other matrices.

const char *const forms =
    "a numeric vector, a 1-d table, a list of single numbers or a one-row "
    "data frame";

// The labels of the rows of x: the row names of a matrix, otherwise
// names(x), which for a 1-d table or array are its dimnames.
SEXP rnames(SEXP x) {
  if (Rf_isMatrix(x)) {
    SEXP dn = Rf_getAttrib(x, R_DimNamesSymbol);
    return Rf_isNull(dn) ? R_NilValue : VECTOR_ELT(dn, 0);
  }
  return Rf_getAttrib(x, R_NamesSymbol);
}

R_xlen_t rows(SEXP x) { return Rf_isMatrix(x) ? Rf_nrows(x) : Rf_xlength(x); }

// The types take() reorders; any other is left in place, to be refused.
bool movable(SEXP x) {
  const int t = TYPEOF(x);
  return t == REALSXP || t == INTSXP || t == LGLSXP || t == STRSXP || t == VECSXP;
}

// x with its rows (the elements of a vector) in the order o, n rows; its
// attributes but names and dimnames are kept (a factor keeps its levels).
SEXP take(SEXP x, const R_xlen_t *o, R_xlen_t n) {
  const R_xlen_t len = XLENGTH(x), nc = len / n;
  const int t = TYPEOF(x);
  SEXP v = PROTECT(Rf_allocVector(t, len));
  for (R_xlen_t c = 0; c < nc; ++c) {
    for (R_xlen_t i = 0; i < n; ++i) {
      const R_xlen_t a = c * n + i, b = c * n + o[i];
      if (t == REALSXP) REAL(v)[a] = REAL_RO(x)[b];
      else if (t == INTSXP) INTEGER(v)[a] = INTEGER_RO(x)[b];
      else if (t == LGLSXP) LOGICAL(v)[a] = LOGICAL_RO(x)[b];
      else if (t == STRSXP) SET_STRING_ELT(v, a, STRING_ELT(x, b));
      else SET_VECTOR_ELT(v, a, VECTOR_ELT(x, b));
    }
  }
  Rf_copyMostAttrib(x, v);
  if (Rf_isMatrix(x)) Rf_setAttrib(v, R_DimSymbol, Rf_getAttrib(x, R_DimSymbol));
  UNPROTECT(1);
  return v;
}

// The names nm (keys k) that are not in the set other, joined by ", ".
const char *absent(SEXP nm, const Keys &k, const Keys &other) {
  std::size_t len = 1;
  for (R_xlen_t i = 0; i < k.n; ++i) {
    if (!has(other, k.k[i])) len += std::strlen(Rf_translateChar(STRING_ELT(nm, i))) + 2;
  }
  char *s = R_alloc(len, 1), *e = s;
  for (R_xlen_t i = 0; i < k.n; ++i) {
    if (has(other, k.k[i])) continue;
    if (e != s) {
      std::memcpy(e, ", ", 2);
      e += 2;
    }
    const char *c = Rf_translateChar(STRING_ELT(nm, i));
    const std::size_t m = std::strlen(c);
    std::memcpy(e, c, m);
    e += m;
  }
  *e = '\0';
  return s;
}

// The keys of the names nm of argument an, refused when a name repeats
// (strict) or returned with nd < n (not strict).
Keys uniq(SEXP nm, const char *an, bool strict) {
  Keys k = keys(nm);
  for (R_xlen_t i = 0; strict && i < k.n; ++i) {
    if (!k.first[i]) {
      fail("%s names a group twice: %s.", an, Rf_translateChar(STRING_ELT(nm, i)));
    }
  }
  return k;
}

// o such that b[o[i]] is a[i], for two sets of unique names (keys ka, kb)
// of one length; refused, listing the names, when the sets differ. Returns
// whether o is the identity.
bool order(SEXP a, const Keys &ka, SEXP b, const Keys &kb, const char *an,
           const char *bn, R_xlen_t *o) {
  bool id = true;
  for (R_xlen_t i = 0; i < ka.n; ++i) {
    const R_xlen_t j = kb.slot[probe(kb, ka.k[i])];
    if (j < 0) {
      fail("%s must name the same groups as %s; missing: %s; extra: %s.", bn,
           an, absent(a, ka, kb), absent(b, kb, ka));
    }
    o[i] = j;
    id = id && j == i;
  }
  return id;
}

// list(NULL, x...): the arguments of one call, in a list that align()
// fills in place (slot 0 for the group names).
SEXP args(std::initializer_list<SEXP> xs) {
  SEXP a = Rf_allocVector(VECSXP, static_cast<R_xlen_t>(xs.size()) + 1);
  R_xlen_t j = 0;
  for (SEXP x : xs) SET_VECTOR_ELT(a, ++j, x);
  return a;
}

// The rule above on the arguments in a (from args()), named nm: later
// arguments are reordered in place and slot 0 receives the group names.
// strict: a repeated name is an error; otherwise (per-unit arguments,
// whose names may be labels that repeat) it leaves positions in use.
SEXP align(SEXP a, const char *const *nm, bool strict) {
  PROTECT(a);
  const int k = static_cast<int>(XLENGTH(a)) - 1;
  SEXP x0 = VECTOR_ELT(a, 1);
  const R_xlen_t n = rows(x0);
  const bool mat = Rf_isMatrix(x0);
  Keys ks[5];
  char an[5][64];
  bool use[5];
  int nu = 0, nn = 0, first = -1;
  bool dup = false;
  for (int j = 0; j < k; ++j) {
    SEXP x = VECTOR_ELT(a, j + 1);
    use[j] = !Rf_isNull(x) && movable(x) && Rf_isMatrix(x) == mat && rows(x) == n;
    std::snprintf(an[j], sizeof an[j], "`%s`", nm[j]);
    if (!use[j]) continue;
    ++nu;
    SEXP r = rnames(x);
    if (Rf_isNull(r)) continue;
    ++nn;
    if (first < 0) first = j;
    ks[j] = uniq(r, an[j], strict);
    dup = dup || ks[j].nd != ks[j].n;
  }
  if (nn == nu && nu > 1 && !dup) {
    R_xlen_t *o = reinterpret_cast<R_xlen_t *>(R_alloc(static_cast<std::size_t>(n), sizeof(R_xlen_t)));
    for (int j = 1; j < k; ++j) {
      if (!use[j]) continue;
      SEXP x = VECTOR_ELT(a, j + 1);
      if (!order(rnames(x0), ks[0], rnames(x), ks[j], an[0], an[j], o)) {
        SET_VECTOR_ELT(a, j + 1, take(x, o, n));
      }
    }
  }
  SET_VECTOR_ELT(a, 0, first < 0 ? clab::seqs(n) : rnames(VECTOR_ELT(a, first + 1)));
  UNPROTECT(1);
  return a;
}

// A per-group argument in any of its forms, as a vector with names: a
// list of single numbers (or NA) and a one-row data frame of numbers give
// their values with their names, a one-row or one-column numeric matrix
// its values named by its column or row names (a 1-d table already has
// names). Any other input is returned as it is; one that is not numeric is
// refused with the forms listed.
SEXP gvec(SEXP x, const char *nm) {
  if (Rf_isNull(x)) return x;
  SEXP v = x;
  if (TYPEOF(x) == VECSXP) {
    // each element (each column of a data frame) one number
    const R_xlen_t n = XLENGTH(x);
    v = PROTECT(Rf_allocVector(REALSXP, n));
    for (R_xlen_t i = 0; i < n; ++i) {
      SEXP e = VECTOR_ELT(x, i);
      clab::unlab(e, nm);
      const bool na = TYPEOF(e) == LGLSXP && XLENGTH(e) == 1 && LOGICAL_RO(e)[0] == NA_LOGICAL;
      if (!na && (!clab::isnum(e) || XLENGTH(e) != 1)) {
        fail("`%s` must be finite numbers, one per group: %s.", nm, forms);
      }
      REAL(v)[i] = clab::val1(e);
    }
    Rf_setAttrib(v, R_NamesSymbol, Rf_getAttrib(x, R_NamesSymbol));
    UNPROTECT(1);
  } else if (Rf_isMatrix(x) && clab::isnum(x) && (Rf_nrows(x) == 1 || Rf_ncols(x) == 1)) {
    SEXP dn = Rf_getAttrib(x, R_DimNamesSymbol);
    clab::unlab(x, nm);
    v = PROTECT(TYPEOF(x) == REALSXP ? dcopy(x) : clab::dbl(x));
    if (!Rf_isNull(dn)) Rf_setAttrib(v, R_NamesSymbol, VECTOR_ELT(dn, Rf_nrows(x) == 1 && Rf_ncols(x) != 1 ? 1 : 0));
    UNPROTECT(1);
  } else if (!clab::isnum(x)) {
    fail("`%s` must be finite numbers, one per group: %s.", nm, forms);
  }
  return v;
}

// gvec() on each argument in a (from args()), then align() strictly.
SEXP group(SEXP a, const char *const *nm) {
  PROTECT(a);
  for (R_xlen_t j = 1; j < XLENGTH(a); ++j) {
    SET_VECTOR_ELT(a, j, gvec(VECTOR_ELT(a, j), nm[j - 1]));
  }
  UNPROTECT(1);
  return align(a, nm, true);
}

// A logical column, TRUE at the reference group i (one-based) of n.
SEXP refcol(int n, int i) {
  SEXP v = Rf_allocVector(LGLSXP, n);
  for (int k = 0; k < n; ++k) LOGICAL(v)[k] = k == i - 1;
  return v;
}

// Doubles (no NaN) as unsigned integers in the same order: the IEEE 754
// bit pattern (R requires IEEE 754), sign bit flipped for positives and
// all bits for negatives; -0 is read as 0 so that the two tie, as they
// compare equal.
constexpr std::uint64_t top = std::uint64_t{1} << 63;

std::uint64_t ukey(double x) {
  const double v = x == 0 ? 0.0 : x;
  std::uint64_t u;
  std::memcpy(&u, &v, sizeof u);
  return (u & top) ? ~u : (u | top);
}

// LSD radix sort of n keys, 16 bits a pass, carrying the int tags o along
// when o is given (stable); a pass on which every key has the same digit
// is skipped. The result is left in k (and o).
void radix(std::uint64_t *k, int *o, int n) {
  std::uint64_t *k2 = reinterpret_cast<std::uint64_t *>(R_alloc(static_cast<std::size_t>(n), sizeof(std::uint64_t)));  // NOLINT
  int *o2 = o ? reinterpret_cast<int *>(R_alloc(static_cast<std::size_t>(n), sizeof(int))) : nullptr;
  int *cnt = reinterpret_cast<int *>(R_alloc(65536, sizeof(int)));
  std::uint64_t *ka = k;
  int *oa = o;
  for (int sh = 0; sh < 64; sh += 16) {
    std::memset(cnt, 0, 65536 * sizeof(int));
    for (int i = 0; i < n; ++i) ++cnt[(ka[i] >> sh) & 0xffff];
    if (n == 0 || cnt[(ka[0] >> sh) & 0xffff] == n) continue;
    for (int d = 0, t = 0; d < 65536; ++d) {
      const int c = cnt[d];
      cnt[d] = t;
      t += c;
    }
    for (int i = 0; i < n; ++i) {
      const int p = cnt[(ka[i] >> sh) & 0xffff]++;
      k2[p] = ka[i];
      if (o) o2[p] = oa[i];
    }
    std::swap(ka, k2);
    std::swap(oa, o2);
  }
  if (ka != k) {
    std::memcpy(k, ka, static_cast<std::size_t>(n) * sizeof(std::uint64_t));
    if (o) std::memcpy(o, oa, static_cast<std::size_t>(n) * sizeof(int));
  }
}

// order(x) for n doubles without NaN: stable, as R's radix order is
void rorder(const double *x, int n, int *o) {
  std::uint64_t *k = reinterpret_cast<std::uint64_t *>(R_alloc(static_cast<std::size_t>(n), sizeof(std::uint64_t)));  // NOLINT
  for (int i = 0; i < n; ++i) {
    k[i] = ukey(x[i]);
    o[i] = i;
  }
  radix(k, o, n);
}

// sort(x) in place for n doubles without NaN
void rsort(double *x, int n) {
  std::uint64_t *k = reinterpret_cast<std::uint64_t *>(R_alloc(static_cast<std::size_t>(n), sizeof(std::uint64_t)));  // NOLINT
  for (int i = 0; i < n; ++i) k[i] = ukey(x[i]);
  radix(k, nullptr, n);
  for (int i = 0; i < n; ++i) {
    const std::uint64_t u = (k[i] & top) ? (k[i] & ~top) : ~k[i];
    std::memcpy(x + i, &u, sizeof u);
  }
}

const char *const vodnm[] = {"group", "ref", "shday", "shnight",
                             "ratio", "sign", "adj"};

// vodtst(day = <table>) with night left out: one table of groups by
// {day, night}, a matrix (a 2-d table) or a data frame whose two columns
// are named day and night, in either order; never taken by position. The
// rows are named by the row names (character row names of a data frame).
// Returns list(day, night).
SEXP split2(SEXP x) {
  const bool df = Rf_inherits(x, "data.frame");
  const bool mt = !df && Rf_isMatrix(x);
  SEXP cn = R_NilValue, rn = R_NilValue;
  if (df) {
    cn = Rf_getAttrib(x, R_NamesSymbol);
    rn = Rf_getAttrib(x, R_RowNamesSymbol);
  } else if (mt) {
    SEXP dn = Rf_getAttrib(x, R_DimNamesSymbol);
    if (!Rf_isNull(dn)) {
      rn = VECTOR_ELT(dn, 0);
      cn = VECTOR_ELT(dn, 1);
    }
  }
  int col[2] = {-1, -1};
  if (Rf_xlength(cn) == 2) {
    static const char *const want[] = {"day", "night"};
    for (int c = 0; c < 2; ++c) {
      for (int w = 0; w < 2; ++w) {
        if (clab::kcmp(clab::key(STRING_ELT(cn, c)), Key{want[w], 0}) == 0) col[w] = c;
      }
    }
  }
  if (col[0] < 0 || col[1] < 0) {
    fail("`night` is missing: give `night`, or `day` as one table (a matrix "
         "or a data frame) with two columns named day and night.");
  }
  clab::unlab(x, "day");
  if (mt && !clab::isnum(x)) fail("`day` must be finite numbers, one per group: %s.", forms);
  const R_xlen_t n = mt ? Rf_nrows(x) : Rf_xlength(VECTOR_ELT(x, 0));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  for (int w = 0; w < 2; ++w) {
    SEXP v;
    if (df) {
      v = Rf_duplicate(VECTOR_ELT(x, col[w]));
    } else {
      SEXP m = PROTECT(clab::dbl(x));
      v = real(n);
      std::memcpy(REAL(v), REAL_RO(m) + col[w] * n, static_cast<std::size_t>(n) * sizeof(double));
      UNPROTECT(1);
    }
    SET_VECTOR_ELT(out, w, v);
    if (TYPEOF(rn) == STRSXP) Rf_setAttrib(v, R_NamesSymbol, rn);
  }
  UNPROTECT(1);
  return out;
}

}  // namespace

extern "C" {

// Disparity.

SEXP C_clab_r_vodtst(SEXP day, SEXP night, SEXP ref, SEXP shift) {
  clab::who = "vodtst";
  SEXP dn = PROTECT(Rf_isNull(night) ? split2(day) : R_NilValue);
  if (!Rf_isNull(dn)) {
    day = VECTOR_ELT(dn, 0);
    night = VECTOR_ELT(dn, 1);
  }
  static const char *const an[] = {"day", "night", "shift"};
  SEXP a = PROTECT(group(args({day, night, shift}), an));
  SEXP nms = VECTOR_ELT(a, 0);
  SEXP d = PROTECT(num(VECTOR_ELT(a, 1), "day", -1, 0.0));
  const int n = ilen(static_cast<double>(XLENGTH(d)), "day");
  SEXP ni = PROTECT(num(VECTOR_ELT(a, 2), "night", n, 0.0, true));
  const bool has = !Rf_isNull(shift);
  SEXP sh = PROTECT(has ? num(VECTOR_ELT(a, 3), "shift", n, 0.0, true) : R_NilValue);
  const int i = clab::ref(ref, nms, n);
  if (REAL_RO(d)[i - 1] <= 0) {
    fail("the reference group has no daylight stops.");
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, has ? 7 : 6));
  SET_VECTOR_ELT(out, 0, nms);
  SET_VECTOR_ELT(out, 1, refcol(n, i));
  for (int k = 2; k < 5; ++k) SET_VECTOR_ELT(out, k, real(n));
  SET_VECTOR_ELT(out, 5, Rf_allocVector(INTSXP, n));
  if (has) SET_VECTOR_ELT(out, 6, real(n));
  kchk(clab_vod(REAL_RO(d), REAL_RO(ni), n, i - 1, has ? REAL_RO(sh) : nullptr,
                REAL(VECTOR_ELT(out, 2)), REAL(VECTOR_ELT(out, 3)),
                REAL(VECTOR_ELT(out, 4)), INTEGER(VECTOR_ELT(out, 5)),
                has ? REAL(VECTOR_ELT(out, 6)) : nullptr),
       nullptr);
  // out[["adj"]] <- ... on a data frame leaves row.names before class
  clab::frame(out, vodnm, has);
  UNPROTECT(6);
  return out;
}

SEXP C_clab_r_orrisk(SEXP orr, SEXP rate, SEXP share) {
  clab::who = "orrisk";
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
  clab::who = "outtst";
  static const char *const an[] = {"q", "w"};
  SEXP a = PROTECT(align(args({q, w}), an, false));
  q = VECTOR_ELT(a, 1);
  w = VECTOR_ELT(a, 2);
  SEXP qv = PROTECT(num(q, "q", -1, 0.0));
  const int n = ilen(static_cast<double>(XLENGTH(qv)), "q");
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
  UNPROTECT(4);
  return clab::nvec(2, nm, r);
}

SEXP C_clab_r_dspbm(SEXP pop, SEXP contact, SEXP force, SEXP ref) {
  clab::who = "dspbm";
  static const char *const an[] = {"pop", "contact", "force"};
  SEXP a = PROTECT(group(args({pop, contact, force}), an));
  SEXP nms = VECTOR_ELT(a, 0);
  SEXP p = PROTECT(num(VECTOR_ELT(a, 1), "pop", -1, 0.0, true));
  const int n = ilen(static_cast<double>(XLENGTH(p)), "pop");
  SEXP c = PROTECT(num(VECTOR_ELT(a, 2), "contact", n, 0.0, true));
  SEXP f = PROTECT(num(VECTOR_ELT(a, 3), "force", n, 0.0, true));
  const int i = clab::ref(ref, nms, n);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
  SET_VECTOR_ELT(out, 0, nms);
  SET_VECTOR_ELT(out, 1, refcol(n, i));
  for (int k = 2; k < 5; ++k) SET_VECTOR_ELT(out, k, real(n));
  kchk(clab_dspbm(REAL_RO(p), REAL_RO(c), REAL_RO(f), n, i - 1,
                  REAL(VECTOR_ELT(out, 2)), REAL(VECTOR_ELT(out, 3)),
                  REAL(VECTOR_ELT(out, 4))),
       nullptr);
  static const char *const nm[] = {"group", "ref", "cd", "fd", "rd"};
  clab::frame(out, nm);
  UNPROTECT(5);
  return out;
}

SEXP C_clab_r_dyadrt(SEXP count, SEXP pop) {
  clab::who = "dyadrt";
  SEXP pv = PROTECT(gvec(pop, "pop"));
  SEXP p = PROTECT(num(pv, "pop", -1, 0.0, true));
  const int g = ilen(static_cast<double>(XLENGTH(p)), "pop");
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
  clab::unlab(count, "count");
  // the copy has no dim, as as.vector() leaves none, so num() reads it as
  // the vector it is
  SEXP cv = PROTECT(ty == REALSXP ? dcopy(count) : clab::dbl(count));
  SEXP x = PROTECT(num(cv, "count", -1, 0.0));
  const R_xlen_t m = XLENGTH(x);
  if (!pos(REAL_RO(x), m)) fail("`count` holds no offences.");
  // The groups by name: when count names both its rows and its columns,
  // the columns are put in the order of the rows, and both in the order of
  // pop when pop is named too. The groups are named by pop, else by the
  // rows, else by the columns.
  SEXP dn = Rf_getAttrib(count, R_DimNamesSymbol);
  SEXP rn = Rf_isNull(dn) ? R_NilValue : VECTOR_ELT(dn, 0);
  SEXP cn = Rf_isNull(dn) ? R_NilValue : VECTOR_ELT(dn, 1);
  SEXP pn = Rf_getAttrib(pv, R_NamesSymbol);
  const char *rd = "the row names of `count`", *cd = "the column names of `count`";
  const Keys kp = Rf_isNull(pn) ? Keys{} : uniq(pn, "`pop`", true);
  const Keys kr = Rf_isNull(rn) ? Keys{} : uniq(rn, rd, true);
  const Keys kc = Rf_isNull(cn) ? Keys{} : uniq(cn, cd, true);
  R_xlen_t *orow = reinterpret_cast<R_xlen_t *>(R_alloc(static_cast<std::size_t>(g), sizeof(R_xlen_t)));
  R_xlen_t *ocol = reinterpret_cast<R_xlen_t *>(R_alloc(static_cast<std::size_t>(g), sizeof(R_xlen_t)));
  for (int i = 0; i < g; ++i) orow[i] = ocol[i] = i;
  if (!Rf_isNull(rn) && !Rf_isNull(cn)) {
    if (!Rf_isNull(pn)) order(pn, kp, rn, kr, "`pop`", rd, orow);
    if (Rf_isNull(pn)) order(rn, kr, cn, kc, rd, cd, ocol);
    else order(pn, kp, cn, kc, "`pop`", cd, ocol);
    double *y = reinterpret_cast<double *>(R_alloc(static_cast<std::size_t>(m), sizeof(double)));
    for (R_xlen_t k = 0; k < m; ++k) y[k] = REAL_RO(x)[orow[k % g] + g * ocol[k / g]];
    std::memcpy(REAL(x), y, static_cast<std::size_t>(m) * sizeof(double));
  }
  SEXP nms = PROTECT(!Rf_isNull(pn) ? pn : !Rf_isNull(rn) ? rn : !Rf_isNull(cn) ? cn : clab::seqs(g));
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
  kchk(clab_dyadrt(REAL_RO(x), off, vic, ilen(static_cast<double>(m), "count"), REAL_RO(p), g,
                   REAL(VECTOR_ELT(out, 3)), REAL(VECTOR_ELT(out, 4)),
                   REAL(VECTOR_ELT(out, 5)), REAL(VECTOR_ELT(out, 6))),
       nullptr);
  static const char *const nm[] = {"off",       "vic",      "count", "rtoff",
                                   "rtpair", "nulloff", "ratio"};
  clab::frame(out, nm);
  UNPROTECT(6);
  return out;
}

}  // extern "C"

// Forensics.

namespace {

// all(is.na(x)) of a q that is not numeric (is.na() on each type as R's
// do_isna has it, the warning included); a list q never comes here, gvec()
// has made it numbers
bool allna(SEXP x) {
  const R_xlen_t n = Rf_xlength(x);
  switch (TYPEOF(x)) {
    case LGLSXP: {
      const int *p = LOGICAL_RO(x);
      for (R_xlen_t i = 0; i < n; ++i) if (p[i] != NA_LOGICAL) return false;
      return true;
    }
    case STRSXP:
      for (R_xlen_t i = 0; i < n; ++i) if (STRING_ELT(x, i) != NA_STRING) return false;
      return true;
    case INTSXP: {  // a factor, or a Date stored as integers
      const int *p = INTEGER_RO(x);
      for (R_xlen_t i = 0; i < n; ++i) if (p[i] != NA_INTEGER) return false;
      return true;
    }
    case REALSXP: {  // a Date, POSIXct or difftime
      const double *p = REAL_RO(x);
      for (R_xlen_t i = 0; i < n; ++i) if (!std::isnan(p[i])) return false;
      return true;
    }
    case CPLXSXP: {
      const Rcomplex *p = COMPLEX_RO(x);
      for (R_xlen_t i = 0; i < n; ++i) {
        if (!std::isnan(p[i].r) && !std::isnan(p[i].i)) return false;
      }
      return true;
    }
    case RAWSXP:
      return n == 0;
    default:
      Rf_warningcall(PROTECT(clab::lang("is.na(q)")),
                     "is.na() applied to non-(list or vector) of type '%s'",
                     Rf_type2char(static_cast<SEXPTYPE>(TYPEOF(x))));
      UNPROTECT(1);
      return n == 0;
  }
}

}  // namespace

extern "C" {

SEXP C_clab_r_gtmp(SEXP p, SEXP q, SEXP hom, SEXP theta) {
  clab::who = "gtmp";
  SEXP p0 = PROTECT(gvec(p, "p"));
  // q may be all NA of any type; only a list (or data frame) is converted
  SEXP q0 = PROTECT(TYPEOF(q) == VECSXP ? gvec(q, "q") : q);
  static const char *const an[] = {"p", "q", "hom"};
  SEXP a = PROTECT(align(args({p0, q0, hom}), an, true));
  SEXP nms = VECTOR_ELT(a, 0);
  q = VECTOR_ELT(a, 2);
  hom = VECTOR_ELT(a, 3);
  SEXP pv = PROTECT(num(VECTOR_ELT(a, 1), "p", -1, 0.0, true));
  const int n = ilen(static_cast<double>(XLENGTH(pv)), "p");
  const double *pp = REAL_RO(pv);
  for (int i = 0; i < n; ++i) {
    if (pp[i] > 1) fail("`p` must be at most 1.");
  }
  // !is.logical(hom) || anyNA(hom) || !length(hom) %in% c(1L, n);
  // rep_len(hom, n)
  const R_xlen_t nh = Rf_xlength(hom);
  if (TYPEOF(hom) != LGLSXP || clab::anyna(hom) || (nh != 1 && nh != n)) {
    fail("`hom` must be TRUE or FALSE, once or once per locus.");
  }
  clab::vec(hom, "hom");
  SEXP h = PROTECT(Rf_allocVector(LGLSXP, n));
  int *hh = LOGICAL(h);
  const int *h0 = LOGICAL_RO(hom);
  for (int i = 0; i < n; ++i) hh[i] = h0[nh == 1 ? 0 : i];
  // if (is.null(q)) q <- rep(NA_real_, n)
  // if (!is.numeric(q) && !all(is.na(q)) || length(q) != n) stop(...)
  // q <- as.double(q); q[hom] <- 0
  SEXP qv = PROTECT(real(n));
  double *qq = REAL(qv);
  bool nq = false;
  if (!Rf_isNull(q)) {
    clab::unlab(q, "q");
    nq = clab::isnum(q);
    if ((!nq && !allna(q)) || Rf_xlength(q) != n) {
      fail("`q` must be numbers, one per locus.");
    }
    clab::vec(q, "q");
  }
  if (nq && TYPEOF(q) == REALSXP) {
    std::memcpy(qq, REAL_RO(q), static_cast<std::size_t>(n) * sizeof(double));
  } else if (nq) {
    const int *qi = INTEGER_RO(q);
    for (int i = 0; i < n; ++i) qq[i] = qi[i] == NA_INTEGER ? NA_REAL : qi[i];
  } else {
    for (int i = 0; i < n; ++i) qq[i] = NA_REAL;
  }
  // anyNA(q) || any(q[het] <= 0) || any(p[het] + q[het] > 1 + 1e-12)
  bool bad = false;
  for (int i = 0; i < n; ++i) {
    if (hh[i]) qq[i] = 0;
    bad = bad || std::isnan(qq[i]) ||
          (!hh[i] && (qq[i] <= 0 || pp[i] + qq[i] > 1 + 1e-12));
  }
  if (bad) {
    fail("`q` must be above 0 with `p + q` at most 1 at every heterozygous "
         "locus.");
  }
  // !is.numeric(theta) || length(theta) != 1L || is.na(theta) ||
  //   theta < 0 || theta >= 1
  const double th = clab::val1(theta);
  if (std::isnan(th) || th < 0 || th >= 1) {
    fail("`theta` must be one number from 0 up to, not including, 1.");
  }
  SEXP loc = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(loc, 0, nms);
  SET_VECTOR_ELT(loc, 1, h);
  SET_VECTOR_ELT(loc, 2, real(n));
  double l10 = 0;
  kchk(clab_gtmp(pp, qq, hh, n, th, REAL(VECTOR_ELT(loc, 2)), &l10), nullptr);
  static const char *const lnm[] = {"locus", "hom", "mp"};
  clab::frame(loc, lnm);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, loc);
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(l10));
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(-l10));
  static const char *const nm[] = {"locus", "log10mp", "log10lr"};
  UNPROTECT(8);
  return clab::nlist(out, nm);
}

SEXP C_clab_r_cllr(SEXP same, SEXP diff) {
  clab::who = "cllr";
  SEXP s = PROTECT(num(same, "same", -1, 0.0, true));
  SEXP d = PROTECT(num(diff, "diff", -1, 0.0, true));
  double r = 0;
  kchk(clab_cllr(REAL_RO(s), ilen(static_cast<double>(XLENGTH(s)), "same"), REAL_RO(d),
                 ilen(static_cast<double>(XLENGTH(d)), "diff"), &r),
       nullptr);
  UNPROTECT(2);
  return Rf_ScalarReal(r);
}

SEXP C_clab_r_lrpost(SEXP lr, SEXP prior) {
  clab::who = "lrpost";
  static const char *const an[] = {"lr", "prior"};
  SEXP a = PROTECT(align(args({lr, prior}), an, false));
  lr = VECTOR_ELT(a, 1);
  prior = VECTOR_ELT(a, 2);
  SEXP l = PROTECT(num(lr, "lr", -1, 0.0, true));
  SEXP p = PROTECT(num(prior, "prior", -1, 0.0, true));
  const R_xlen_t nl = XLENGTH(l), np = XLENGTH(p);
  const double *lv = REAL_RO(l), *pv = REAL_RO(p);
  // any(prior >= 1) || !length(prior) %in% c(1L, length(lr))
  bool bad = np != 1 && np != nl;
  for (R_xlen_t i = 0; i < np && !bad; ++i) bad = pv[i] >= 1;
  if (bad) fail("`prior` must be below 1, once or once per ratio.");
  // 1 / (1 + (1 - prior) / (prior * lr)): no overflow for a huge ratio
  SEXP out = PROTECT(real(nl));
  double *o = REAL(out);
  for (R_xlen_t i = 0; i < nl; ++i) {
    const double q = pv[np == 1 ? 0 : i];
    o[i] = 1 / (1 + (1 - q) / (q * lv[i]));
  }
  UNPROTECT(4);
  return out;
}

}  // extern "C"

namespace {

// .mix(): one side of a mixture ratio, the chance that x unknowns carry
// only alleles of `p` with every allele of `need`. src is the old call.
double mix(SEXP pv, const Keys &nk, SEXP need, SEXP x, const char *a,
           const char *b, const char *src) {
  // !is.character(need) || anyNA(need) || !all(need %in% nms) ||
  //   length(unique(need)) > 20L
  bool bad = TYPEOF(need) != STRSXP || clab::anyna(need);
  Keys nd{};
  if (!bad) {
    nd = keys(need);
    for (R_xlen_t i = 0; i < nd.n && !bad; ++i) bad = !has(nk, nd.k[i]);
    bad = bad || nd.nd > 20;
  }
  if (bad) fail("`%s` must be names of `p`, at most 20.", a);
  clab::vec(need, a);
  // !is.numeric(x) || length(x) != 1L || is.na(x) || x < 0 ||
  //   x != round(x) || x > 50
  const double v = clab::val1(x);
  if (std::isnan(v) || v < 0 || v != std::floor(v) || v > 50) {
    fail("`%s` must be a whole number from 0 to 50.", b);
  }
  // nms %in% need
  const R_xlen_t n = XLENGTH(pv);
  int *in = reinterpret_cast<int *>(R_alloc(static_cast<std::size_t>(n), sizeof(int)));
  for (R_xlen_t i = 0; i < n; ++i) in[i] = has(nd, nk.k[i]);
  double r = 0;
  kchk(clab_mixp(REAL_RO(pv), in, static_cast<int>(n), static_cast<int>(v), &r),
       src);
  return r;
}

}  // namespace

extern "C" {

SEXP C_clab_r_mixlr(SEXP p, SEXP np, SEXP xp, SEXP nd, SEXP xd) {
  clab::who = "mixlr";
  SEXP nms = PROTECT(Rf_getAttrib(p, R_NamesSymbol));
  SEXP pv = PROTECT(num(p, "p", -1, 0.0, true));
  const int n = ilen(static_cast<double>(XLENGTH(pv)), "p");
  // is.null(nms) || anyDuplicated(nms) || sum(p) > 1 + 1e-12
  bool bad = Rf_isNull(nms);
  Keys nk{};
  if (!bad) {
    nk = keys(nms);
    bad = nk.nd != nk.n;
  }
  if (bad || clab::nsum(REAL_RO(pv), n) > 1 + 1e-12) {
    fail("`p` must be uniquely named with a sum of at most 1.");
  }
  double r[3];
  r[0] = mix(pv, nk, np, xp, "np", "xp",
             ".mix(p, nms, np, xp, \"np\", \"xp\")");
  r[1] = mix(pv, nk, nd, xd, "nd", "xd",
             ".mix(p, nms, nd, xd, \"nd\", \"xd\")");
  r[2] = r[0] / r[1];
  static const char *const nm[] = {"num", "den", "lr"};
  UNPROTECT(2);
  return clab::nvec(3, nm, r);
}

SEXP C_clab_r_dbsrch(SEXP mp, SEXP n, SEXP pop) {
  clab::who = "dbsrch";
  SEXP m = PROTECT(num(mp, "mp", 1, 0.0, true));
  const double mv = REAL_RO(m)[0];
  if (mv > 1) fail("`mp` must be at most 1.");
  SEXP nn = PROTECT(num(n, "n", 1, 1.0));
  const double nv = REAL_RO(nn)[0];
  double r[3];
  // c(np = min(1, n * mp), any = -expm1(n * log1p(-mp)))
  const double e = nv * mv;
  r[0] = e < 1 ? e : 1;
  r[1] = -std::expm1(nv * std::log1p(-mv));
  const bool has = !Rf_isNull(pop);
  if (has) {
    SEXP pp = PROTECT(num(pop, "pop", 1, nv));
    // 1 / (1 + (pop - n) * mp)
    r[2] = 1 / (1 + clab::rd((REAL_RO(pp)[0] - nv) * mv));
    UNPROTECT(1);
  }
  static const char *const nm[] = {"np", "any", "post"};
  UNPROTECT(2);
  return clab::nvec(has ? 3 : 2, nm, r);
}

}  // extern "C"

namespace {

// .cal(): the ratios checked, and calibrated in the order given. Leaves
// four objects PROTECTed (s, d, cs, cd).
struct Cal {
  SEXP s, d, cs, cd;
  double min;
};

Cal cal(SEXP same, SEXP diff) {
  Cal c{};
  c.s = PROTECT(num(same, "same", -1, 0.0, true));
  c.d = PROTECT(num(diff, "diff", -1, 0.0, true));
  const int ns = ilen(static_cast<double>(XLENGTH(c.s)), "same");
  const int nd = ilen(static_cast<double>(XLENGTH(c.d)), "diff");
  const int n = ilen(static_cast<double>(ns) + nd, "same");
  const double *sp = REAL_RO(c.s), *dp = REAL_RO(c.d);
  // lr <- c(same, diff); o <- order(lr): a stable order, as R's radix sort
  double *lr = reinterpret_cast<double *>(R_alloc(static_cast<std::size_t>(n), sizeof(double)));
  std::memcpy(lr, sp, static_cast<std::size_t>(ns) * sizeof(double));
  std::memcpy(lr + ns, dp, static_cast<std::size_t>(nd) * sizeof(double));
  int *o = reinterpret_cast<int *>(R_alloc(static_cast<std::size_t>(n), sizeof(int)));
  rorder(lr, n, o);
  int *is = reinterpret_cast<int *>(R_alloc(static_cast<std::size_t>(n), sizeof(int)));
  for (int k = 0; k < n; ++k) is[k] = o[k] < ns;
  // lr[o], in place through the stable order: copy the sorted values
  double *ls = reinterpret_cast<double *>(R_alloc(static_cast<std::size_t>(n), sizeof(double)));
  for (int k = 0; k < n; ++k) ls[k] = lr[o[k]];
  lr = ls;
  double *post = reinterpret_cast<double *>(R_alloc(static_cast<std::size_t>(n), sizeof(double)));
  double *work = reinterpret_cast<double *>(R_alloc(2 * static_cast<std::size_t>(n), sizeof(double)));
  kchk(clab_lrcal(lr, is, n, post, work, &c.min), ".cal(same, diff)");
  // post / (1 - post) * length(diff) / length(same), back in input order
  c.cs = PROTECT(real(ns));
  c.cd = PROTECT(real(nd));
  double *cs = REAL(c.cs), *cd = REAL(c.cd);
  for (int k = 0; k < n; ++k) {
    const double x = post[k] / (1 - post[k]) * nd / ns;
    if (o[k] < ns) cs[o[k]] = x; else cd[o[k] - ns] = x;
  }
  return c;
}

// one curve of tipplt(): sorted log10 ratios and the share at or above
SEXP curve(SEXP x) {
  const R_xlen_t n = XLENGTH(x);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP l = real(n);
  SET_VECTOR_ELT(out, 0, l);
  double *lp = REAL(l);
  const double *xp = REAL_RO(x);
  // sort(log10(x)); x passed ilen() in the caller
  for (R_xlen_t i = 0; i < n; ++i) lp[i] = std::log10(xp[i]);
  rsort(lp, static_cast<int>(n));
  SEXP s = real(n);
  SET_VECTOR_ELT(out, 1, s);
  double *sp = REAL(s);
  // rev(seq_along(x)) / length(x)
  for (R_xlen_t i = 0; i < n; ++i) {
    sp[i] = static_cast<double>(n - i) / static_cast<double>(n);
  }
  static const char *const nm[] = {"log10lr", "share"};
  clab::frame(out, nm);
  UNPROTECT(1);
  return out;
}

// log(1 + e^z), without overflow for large z
double softplus(double z) {
  return z > 0 ? z + std::log1p(std::exp(-z)) : std::log1p(std::exp(z));
}

// One empirical cross-entropy:
// (pr * mean(log(1 + 1 / (s * odds))) + (1 - pr) * mean(log1p(d * odds)))
//   / log(2), lo the prior log10 odds, odds = 10^lo.
// Where s * odds (or d * odds) is a positive finite number the terms are
// computed as before. Where it is not (a calibrated ratio of Inf or 0,
// odds that overflow to Inf or underflow to 0, or a product out of range)
// the term is taken in logs, log x = log(s) + lo * log(10), so that no
// 0 * Inf is formed: log(1 + 1 / x) = softplus(-log x) and
// log(1 + x) = softplus(log x). A same-source ratio of Inf and a
// different-source ratio of 0 then add 0, their limit, at any odds.
double ece(const double *s, R_xlen_t ns, const double *d, R_xlen_t nd,
           double lo, double odds, double pr, double *t) {
  const double lg = lo * std::log(10.0);
  const bool fin = std::isfinite(odds) && odds > 0;
  for (R_xlen_t k = 0; k < ns; ++k) {
    const double x = fin && std::isfinite(s[k]) && s[k] > 0 ? s[k] * odds : 0;
    if (std::isfinite(x) && x > 0) {
      // log(1 + 1 / x) without forming 1 / x where it overflows
      t[k] = x < 1 ? std::log1p(x) - std::log(x) : std::log1p(1 / x);
    } else {
      t[k] = softplus(-(std::log(s[k]) + lg));
    }
  }
  const double a = clab::nmean(t, ns);
  for (R_xlen_t k = 0; k < nd; ++k) {
    const double x = fin && std::isfinite(d[k]) && d[k] > 0 ? d[k] * odds : 0;
    t[k] = std::isfinite(x) && x > 0 ? std::log1p(x) : softplus(std::log(d[k]) + lg);
  }
  const double b = clab::nmean(t, nd);
  return (clab::rd(pr * a) + clab::rd((1 - pr) * b)) / std::log(2.0);
}

}  // namespace

extern "C" {

SEXP C_clab_r_cllrmin(SEXP same, SEXP diff) {
  clab::who = "cllrmin";
  Cal c = cal(same, diff);
  double r[3];
  // cllr(r[["same"]], r[["diff"]]): its .Call ran in that call's frame
  kchk(clab_cllr(REAL_RO(c.s), static_cast<int>(XLENGTH(c.s)), REAL_RO(c.d),
                 static_cast<int>(XLENGTH(c.d)), r),
       "cllr(r[[\"same\"]], r[[\"diff\"]])");
  r[1] = c.min;
  r[2] = r[0] - c.min;
  static const char *const nm[] = {"cllr", "min", "cal"};
  UNPROTECT(4);
  return clab::nvec(3, nm, r);
}

SEXP C_clab_r_lrcal(SEXP same, SEXP diff) {
  clab::who = "lrcal";
  Cal c = cal(same, diff);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, c.cs);
  SET_VECTOR_ELT(out, 1, c.cd);
  static const char *const nm[] = {"same", "diff"};
  UNPROTECT(5);
  return clab::nlist(out, nm);
}

SEXP C_clab_r_tipplt(SEXP same, SEXP diff) {
  clab::who = "tipplt";
  SEXP s = PROTECT(num(same, "same", -1, 0.0, true));
  SEXP d = PROTECT(num(diff, "diff", -1, 0.0, true));
  ilen(static_cast<double>(XLENGTH(s)), "same");
  ilen(static_cast<double>(XLENGTH(d)), "diff");
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, curve(s));
  SET_VECTOR_ELT(out, 1, curve(d));
  static const char *const nm[] = {"same", "diff"};
  UNPROTECT(3);
  return clab::nlist(out, nm);
}

SEXP C_clab_r_eceplt(SEXP same, SEXP diff, SEXP lo) {
  clab::who = "eceplt";
  Cal c = cal(same, diff);
  SEXP l = PROTECT(num(lo, "lo"));
  const R_xlen_t nl = XLENGTH(l), ns = XLENGTH(c.s), nd = XLENGTH(c.d);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out, 0, dcopy(l));
  for (int k = 1; k < 4; ++k) SET_VECTOR_ELT(out, k, real(nl));
  double *e = REAL(VECTOR_ELT(out, 1)), *m = REAL(VECTOR_ELT(out, 2)),
         *r = REAL(VECTOR_ELT(out, 3));
  double *t = reinterpret_cast<double *>(
      R_alloc(static_cast<std::size_t>(ns > nd ? ns : nd), sizeof(double)));
  const double *lv = REAL_RO(l), one = 1.0;
  for (R_xlen_t i = 0; i < nl; ++i) {
    const double odds = clab::rpow(10.0, lv[i]);
    // the prior: odds / (1 + odds), 1 where the odds overflow
    const double pr = std::isinf(odds) ? 1.0 : odds / (1 + odds);
    e[i] = ece(REAL_RO(c.s), ns, REAL_RO(c.d), nd, lv[i], odds, pr, t);
    m[i] = ece(REAL_RO(c.cs), ns, REAL_RO(c.cd), nd, lv[i], odds, pr, t);
    r[i] = ece(&one, 1, &one, 1, lv[i], odds, pr, t);
  }
  static const char *const nm[] = {"lo", "ece", "min", "ref"};
  clab::frame(out, nm);
  UNPROTECT(6);
  return out;
}

}  // extern "C"

namespace {

// matrix(x, 1L) for a vector of two alleles: a factor gives its labels,
// another class its data (matrix() calls as.vector()); a type that cannot
// hold alleles is returned unchanged, to be refused as before.
SEXP mat12(SEXP x) {
  SEXP v = x;
  if (Rf_isFactor(x)) v = clab::labs(x);
  PROTECT(v);
  const int t = TYPEOF(v);
  if (t != INTSXP && t != REALSXP && t != STRSXP) {
    UNPROTECT(1);
    return x;
  }
  SEXP m = PROTECT(Rf_allocMatrix(t, 1, 2));
  for (R_xlen_t i = 0; i < 2; ++i) {
    if (t == STRSXP) SET_STRING_ELT(m, i, STRING_ELT(v, i));
    else if (t == INTSXP) INTEGER(m)[i] = INTEGER_RO(v)[i];
    else REAL(m)[i] = REAL_RO(v)[i];
  }
  UNPROTECT(2);
  return m;
}

int nrow(SEXP x) { return INTEGER_RO(Rf_getAttrib(x, R_DimSymbol))[0]; }

// .gt(): a genotype table, one row per locus, two alleles; n rows if n >= 0
SEXP gt(SEXP x, const char *nm, int n) {
  // if (is.null(dim(x)) && length(x) == 2L) x <- matrix(x, 1L)
  clab::unlab(x, nm);
  if (Rf_isNull(Rf_getAttrib(x, R_DimSymbol)) && Rf_xlength(x) == 2) x = mat12(x);
  PROTECT(x);
  // !is.matrix(x) || ncol(x) != 2L || anyNA(x) ||
  //   !(is.numeric(x) || is.character(x))
  const bool ok = Rf_isMatrix(x) &&
                  INTEGER_RO(Rf_getAttrib(x, R_DimSymbol))[1] == 2 &&
                  !clab::anyna(x) && (clab::isnum(x) || TYPEOF(x) == STRSXP);
  if (!ok) {
    fail("`%s` must be a two-column matrix of alleles, one row per locus, "
         "with no missing values.", nm);
  }
  if (n >= 0 && nrow(x) != n) {
    fail("`%s` must have %d rows, not %d.", nm, n, nrow(x));
  }
  UNPROTECT(1);
  return x;
}

// .gf(): the frequencies of the alleles of g, same shape, row by row; a
// pair of frequencies is already a one-row table (loci() makes it one)
const double *gf(SEXP freq, SEXP g, const char *nm) {
  PROTECT(freq);
  // !is.matrix(freq) || !is.numeric(freq) || !identical(dim(freq), dim(g))
  //   || anyNA(freq) || any(freq <= 0) || any(freq > 1)
  const int nr = nrow(g);
  bool bad = !Rf_isMatrix(freq) || !clab::isnum(freq);
  if (!bad) {
    const int *df = INTEGER_RO(Rf_getAttrib(freq, R_DimSymbol));
    bad = df[0] != nr || df[1] != 2;
  }
  SEXP fv = PROTECT(bad ? R_NilValue : clab::dbl(freq));
  const double *f = bad ? nullptr : REAL_RO(fv);
  for (R_xlen_t i = 0; !bad && i < 2 * static_cast<R_xlen_t>(nr); ++i) {
    bad = std::isnan(f[i]) || f[i] <= 0 || f[i] > 1;
  }
  if (bad) {
    fail("`freq` must hold the frequencies of the alleles of `%s`, same "
         "shape, above 0 and at most 1.", nm);
  }
  // any(g[, 1L] == g[, 2L] & freq[, 1L] != freq[, 2L])
  SEXP gv = PROTECT(TYPEOF(g) == STRSXP ? g : clab::dbl(g));
  for (int i = 0; i < nr; ++i) {
    const bool same = TYPEOF(gv) == STRSXP
                          ? clab::seql(STRING_ELT(gv, i), STRING_ELT(gv, i + nr))
                          : REAL_RO(gv)[i] == REAL_RO(gv)[i + nr];
    if (same && f[i] != f[i + nr]) {
      fail("`freq` gives two frequencies for one allele of a homozygote.");
    }
  }
  // as.double(t(freq))
  double *out = reinterpret_cast<double *>(
      R_alloc(2 * static_cast<std::size_t>(nr), sizeof(double)));
  for (int i = 0; i < nr; ++i) {
    out[2 * i] = f[i];
    out[2 * i + 1] = f[i + nr];
  }
  UNPROTECT(3);
  return out;
}

//' @srrstats {G2.3} The character arguments rel (kinidx()), ref (vodtst(),
//'   dspbm(), expbnd()), treat (sntbnd()) and fn (proof()) are checked
//'   against the values allowed, and refused otherwise.
//' @srrstats {G2.3a} kin() matches rel against its four names, as
//'   match.arg() would, and lists them when it refuses.
//' @srrstats {G2.3b} rel, ref and treat are matched exactly, case
//'   included, as their help pages say.
//' @noRd
// .kin(): the identity-by-descent probabilities of a relationship
void kin(SEXP rel, double *k) {
  static const char *const nm[] = {"parent", "sibling", "half", "cousin"};
  static const double ks[4][3] = {
      {0, 1, 0}, {0.25, 0.5, 0.25}, {0.5, 0.5, 0}, {0.75, 0.25, 0}};
  // is.character(rel) && length(rel) == 1L && rel %in% names(ks)
  if (TYPEOF(rel) == STRSXP && XLENGTH(rel) == 1) {
    const Key r = clab::key(STRING_ELT(rel, 0));
    for (int i = 0; i < 4; ++i) {
      if (clab::kcmp(r, Key{nm[i], 0}) == 0) {
        std::memcpy(k, ks[i], sizeof ks[i]);
        return;
      }
    }
  }
  // is.numeric(rel) && length(rel) == 3L && !anyNA(rel) && all(rel >= 0) &&
  //   abs(sum(rel) - 1) < 1e-9
  if (clab::isnum(rel) && XLENGTH(rel) == 3) {
    SEXP v = PROTECT(clab::dbl(rel));
    const double *p = REAL_RO(v);
    bool ok = true;
    for (int i = 0; i < 3; ++i) ok = ok && !std::isnan(p[i]) && p[i] >= 0;
    if (ok && std::fabs(clab::nsum(p, 3) - 1) < 1e-9) {
      std::memcpy(k, p, 3 * sizeof(double));
      UNPROTECT(1);
      return;
    }
    UNPROTECT(1);
  }
  fail("`rel` must be one of parent, sibling, half, cousin, or three "
       "probabilities that sum to 1.");
}

// .code(): the alleles of m tables of nr rows as integer codes, row by row,
// one code per distinct label of as.character(); the kernels only compare
// codes, so a code is the label's rank in code-point order (R numbered
// labels by first appearance), found by one sort, not pairwise.
void codes(SEXP *g, int m, int nr, int **out, const char *nm) {
  // the hash table below holds up to 8 * m * nr slots, all of them ints
  ilen(8.0 * m * nr, nm);
  const int w = 2 * nr, all = m * w;
  for (int j = 0; j < m; ++j) {
    out[j] = reinterpret_cast<int *>(R_alloc(static_cast<std::size_t>(w), sizeof(int)));
  }
  // u <- unique(c(as.character(g1), ...)); match(as.character(g), u), by
  // one open-addressing hash table: a label's code is one plus the flat
  // position of its first appearance (R numbered the labels 1, 2, ...;
  // the kernels only compare codes). Allele p of the flat order (table j,
  // column c, row i) gets its code in out[j][2 * i + c].
  int cap = 2, bits = 1;
  while (cap < 2 * all) {
    cap *= 2;
    ++bits;
  }
  int *slot = reinterpret_cast<int *>(R_alloc(static_cast<std::size_t>(cap), sizeof(int)));
  for (int i = 0; i < cap; ++i) slot[i] = -1;
  // Whole numbers below 1e15 in size have one text each, and two of them
  // have the same text exactly when they are equal: such tables are coded
  // by value, without making a string.
  bool whole = true;
  for (int j = 0; j < m && whole; ++j) whole = TYPEOF(g[j]) != STRSXP;
  auto val = [g, w](int p) {
    SEXP x = g[p / w];
    return TYPEOF(x) == INTSXP ? INTEGER_RO(x)[p % w] : REAL_RO(x)[p % w];
  };
  for (int p = 0; p < all && whole; ++p) {
    const double v = val(p);
    whole = v == std::floor(v) && std::fabs(v) < 1e15;
  }
  Key *k = nullptr;
  if (!whole) {
    k = reinterpret_cast<Key *>(R_alloc(static_cast<std::size_t>(all), sizeof(Key)));
    for (int j = 0; j < m; ++j) {
      SEXP s = PROTECT(clab::labs(g[j]));
      for (int i = 0; i < w; ++i) k[j * w + i] = clab::key(STRING_ELT(s, i));
      UNPROTECT(1);
    }
  }
  for (int p = 0; p < all; ++p) {
    // Fibonacci hashing: the top bits of the product mix every key bit
    const std::uint64_t h = (whole ? ukey(val(p)) : hash(k[p])) * 0x9e3779b97f4a7c15u;
    int i = static_cast<int>(h >> (64 - bits));
    while (slot[i] >= 0 && (whole ? val(slot[i]) != val(p) : clab::kcmp(k[slot[i]], k[p]) != 0)) {
      i = (i + 1) & (cap - 1);
    }
    if (slot[i] < 0) slot[i] = p;
    const int j = p / w, c = (p % w) / nr, r = p % nr;
    out[j][2 * r + c] = slot[i] + 1;
  }
}

// .idx(): the per-locus indices with names(g[, 1L]) as loci, their log10
// total and the posterior at even prior odds
SEXP idx(SEXP g, SEXP ix) {
  const int nr = nrow(g);
  // g[, 1L] keeps the row names; a one-row table drops to one value, which
  // keeps a name only when exactly one of its row and column names exists
  SEXP dn = Rf_getAttrib(g, R_DimNamesSymbol);
  SEXP rn = Rf_isNull(dn) ? R_NilValue : VECTOR_ELT(dn, 0);
  SEXP cn = Rf_isNull(dn) ? R_NilValue : VECTOR_ELT(dn, 1);
  SEXP nm;
  if (nr == 1 && Rf_isNull(rn) != Rf_isNull(cn)) {
    nm = Rf_ScalarString(STRING_ELT(Rf_isNull(rn) ? cn : rn, 0));
  } else if (nr > 1 && !Rf_isNull(rn)) {
    nm = rn;
  } else {
    nm = clab::seqs(nr);
  }
  PROTECT(nm);
  // l10 <- sum(log10(idx)); post = 1 / (1 + 10^(-l10))
  double *lg = reinterpret_cast<double *>(R_alloc(static_cast<std::size_t>(nr), sizeof(double)));
  for (int i = 0; i < nr; ++i) lg[i] = std::log10(REAL_RO(ix)[i]);
  const double l10 = clab::nsum(lg, nr);
  SEXP loc = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(loc, 0, nm);
  SET_VECTOR_ELT(loc, 1, ix);
  static const char *const lnm[] = {"locus", "idx"};
  clab::frame(loc, lnm);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, loc);
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(l10));
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(1 / (1 + clab::rpow(10.0, -l10))));
  static const char *const onm[] = {"locus", "log10idx", "post"};
  UNPROTECT(3);
  return clab::nlist(out, onm);
}

// The genotype tables and freq, matched by locus (row names) under the
// rule of align(): list(loci, tables..., freq), each table from gt() and
// freq made a table as gf() makes it. src receives the table whose row
// names name the loci (the first one named, else the first).
SEXP loci(std::initializer_list<SEXP> gs, SEXP freq, const char *const *nm,
          SEXP *src) {
  if (Rf_isNull(Rf_getAttrib(freq, R_DimSymbol)) && Rf_xlength(freq) == 2) {
    freq = mat12(freq);
  }
  PROTECT(freq);
  SEXP a = PROTECT(args(gs));
  const R_xlen_t k = XLENGTH(a);
  SEXP b = PROTECT(Rf_allocVector(VECSXP, k + 1));
  for (R_xlen_t j = 1; j < k; ++j) SET_VECTOR_ELT(b, j, VECTOR_ELT(a, j));
  SET_VECTOR_ELT(b, k, freq);
  align(b, nm, true);
  *src = VECTOR_ELT(b, 1);
  for (R_xlen_t j = k; j >= 1; --j) {
    if (!Rf_isNull(rnames(VECTOR_ELT(b, j)))) *src = VECTOR_ELT(b, j);
  }
  UNPROTECT(3);
  return b;
}

SEXP kinidx(SEXP g1, SEXP g2, SEXP freq, SEXP rel) {
  SEXP a = PROTECT(gt(g1, "g1", -1));
  const int nr = nrow(a);
  SEXP b = PROTECT(gt(g2, "g2", nr));
  static const char *const an[] = {"g1", "g2", "freq"};
  SEXP src;
  SEXP al = PROTECT(loci({a, b}, freq, an, &src));
  b = VECTOR_ELT(al, 2);
  const double *f = gf(VECTOR_ELT(al, 3), a, "g1");
  double k[3];
  kin(rel, k);
  SEXP g[2] = {a, b};
  int *c[2];
  codes(g, 2, nr, c, "g1");
  SEXP ix = PROTECT(real(nr));
  kchk(clab_kinidx(c[0], c[1], f, nr, k, REAL(ix)),
       ".idx(g1, .Call(C_clab_kinidx, .code(g1, u), .code(g2, u), f, k))");
  SEXP out = idx(src, ix);
  UNPROTECT(4);
  return out;
}

}  // namespace

extern "C" {

SEXP C_clab_r_kinidx(SEXP g1, SEXP g2, SEXP freq, SEXP rel) {
  clab::who = "kinidx";
  return kinidx(g1, g2, freq, rel);
}

SEXP C_clab_r_patidx(SEXP child, SEXP father, SEXP freq, SEXP mother) {
  clab::who = "patidx";
  // if (is.null(mother)) return(kinidx(child, father, freq, "parent"))
  if (Rf_isNull(mother)) {
    SEXP rel = PROTECT(Rf_mkString("parent"));
    SEXP out = kinidx(child, father, freq, rel);
    UNPROTECT(1);
    return out;
  }
  SEXP c = PROTECT(gt(child, "child", -1));
  const int nr = nrow(c);
  SEXP fa = PROTECT(gt(father, "father", nr));
  SEXP mo = PROTECT(gt(mother, "mother", nr));
  static const char *const an[] = {"child", "father", "mother", "freq"};
  SEXP src;
  SEXP al = PROTECT(loci({c, fa, mo}, freq, an, &src));
  fa = VECTOR_ELT(al, 2);
  mo = VECTOR_ELT(al, 3);
  const double *f = gf(VECTOR_ELT(al, 4), c, "child");
  SEXP g[3] = {c, fa, mo};
  int *k[3];
  codes(g, 3, nr, k, "child");
  SEXP ix = PROTECT(real(nr));
  kchk(clab_patidx(k[0], k[2], k[1], f, nr, REAL(ix)),
       ".idx(child, .Call(C_clab_patidx, .code(child, u), .code(mother, u), "
       ".code(father, u), f))");
  SEXP out = idx(src, ix);
  UNPROTECT(5);
  return out;
}

}  // extern "C"

// Policing.

extern "C" {

SEXP C_clab_r_clrsmp(SEXP clra, SEXP reca, SEXP clrb, SEXP recb) {
  clab::who = "clrsmp";
  static const char *const an[] = {"clra", "reca", "clrb", "recb"};
  SEXP al = PROTECT(group(args({clra, reca, clrb, recb}), an));
  clra = VECTOR_ELT(al, 1);
  reca = VECTOR_ELT(al, 2);
  clrb = VECTOR_ELT(al, 3);
  recb = VECTOR_ELT(al, 4);
  SEXP nms = PROTECT(VECTOR_ELT(al, 0));
  SEXP a = PROTECT(num(clra, "clra", -1, 0.0));
  const int n = ilen(static_cast<double>(XLENGTH(a)), "clra");
  SEXP ra = PROTECT(num(reca, "reca", n, 0.0, true));
  SEXP b = PROTECT(num(clrb, "clrb", n, 0.0));
  SEXP rb = PROTECT(num(recb, "recb", n, 0.0, true));
  SEXP parts = PROTECT(Rf_allocVector(VECSXP, 7));
  SET_VECTOR_ELT(parts, 0, nms);
  for (int k = 1; k < 7; ++k) SET_VECTOR_ELT(parts, k, real(n));
  double tot[8];
  int cnt[3];
  kchk(clab_clrsmp(REAL_RO(a), REAL_RO(ra), REAL_RO(b), REAL_RO(rb), n,
                   REAL(VECTOR_ELT(parts, 1)), REAL(VECTOR_ELT(parts, 2)),
                   REAL(VECTOR_ELT(parts, 3)), REAL(VECTOR_ELT(parts, 4)),
                   REAL(VECTOR_ELT(parts, 5)), REAL(VECTOR_ELT(parts, 6)), tot,
                   cnt),
       nullptr);
  static const char *const pnm[] = {"part",   "ratea",   "rateb", "sharea",
                                    "shareb", "rateeff", "mixeff"};
  clab::frame(parts, pnm);
  // c(..., up = cnt[1L], ...) makes the counts doubles; they are never NA
  for (int k = 0; k < 3; ++k) tot[5 + k] = cnt[k];
  static const char *const tnm[] = {"ratea",   "rateb", "change", "rateeff",
                                    "mixeff", "up",    "down",   "simpson"};
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, parts);
  SET_VECTOR_ELT(out, 1, clab::nvec(8, tnm, tot));
  static const char *const nm[] = {"parts", "total"};
  UNPROTECT(8);
  return clab::nlist(out, nm);
}

SEXP C_clab_r_expbnd(SEXP y, SEXP m, SEXP gamma, SEXP ref) {
  clab::who = "expbnd";
  static const char *const an[] = {"y", "m"};
  SEXP a = PROTECT(group(args({y, m}), an));
  y = VECTOR_ELT(a, 1);
  m = VECTOR_ELT(a, 2);
  SEXP nms = PROTECT(VECTOR_ELT(a, 0));
  SEXP yv = PROTECT(num(y, "y", -1, 0.0, true));
  const int n = ilen(static_cast<double>(XLENGTH(yv)), "y");
  SEXP mv = PROTECT(num(m, "m", n, 0.0, true));
  SEXP g = PROTECT(num(gamma, "gamma", 1, 1.0));
  const int i = clab::ref(ref, nms, n);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 9));
  SET_VECTOR_ELT(out, 0, nms);
  SET_VECTOR_ELT(out, 1, refcol(n, i));
  for (int k = 2; k < 8; ++k) SET_VECTOR_ELT(out, k, real(n));
  SET_VECTOR_ELT(out, 8, Rf_allocVector(INTSXP, n));
  kchk(clab_expbnd(REAL_RO(yv), REAL_RO(mv), n, i - 1, REAL_RO(g)[0],
                   REAL(VECTOR_ELT(out, 2)), REAL(VECTOR_ELT(out, 3)),
                   REAL(VECTOR_ELT(out, 4)), REAL(VECTOR_ELT(out, 5)),
                   REAL(VECTOR_ELT(out, 6)), REAL(VECTOR_ELT(out, 7)),
                   INTEGER(VECTOR_ELT(out, 8))),
       nullptr);
  static const char *const nm[] = {"group", "ref", "rate", "lo",  "hi",
                                   "ratio", "rlo", "rhi",  "sign"};
  clab::frame(out, nm);
  UNPROTECT(6);
  return out;
}

SEXP C_clab_r_misrd(SEXP r1, SEXP r0, SEXP prev, SEXP se, SEXP sp) {
  clab::who = "misrd";
  const double a = clab::p01(r1, "r1");
  const double b = clab::p01(r0, "r0");
  const double p = clab::one(prev, "prev", 0, 1);
  const double s = clab::p01(se, "se");
  const double c = clab::p01(sp, "sp");
  // d1 <- prev * se + (1 - prev) * (1 - sp)
  const double d1 = clab::rd(p * s) + clab::rd((1 - p) * (1 - c));
  if (d1 <= 0 || d1 >= 1) {
    fail("`se` and `sp` leave a recorded group empty.");
  }
  double r[5];
  kchk(clab_misrd(a, b, p, s, c, r), nullptr);
  static const char *const nm[] = {"true", "measured", "factor", "d1", "d0"};
  return clab::nvec(5, nm, r);
}

SEXP C_clab_r_rtmhs(SEXP x1, SEXP x2, SEXP cut, SEXP w) {
  clab::who = "rtmhs";
  static const char *const an[] = {"x1", "x2", "w"};
  SEXP al = PROTECT(align(args({x1, x2, w}), an, false));
  x1 = VECTOR_ELT(al, 1);
  x2 = VECTOR_ELT(al, 2);
  w = VECTOR_ELT(al, 3);
  SEXP a = PROTECT(num(x1, "x1"));
  const int n = ilen(static_cast<double>(XLENGTH(a)), "x1");
  SEXP b = PROTECT(num(x2, "x2", n));
  SEXP c = PROTECT(num(cut, "cut", 1));
  SEXP wv = PROTECT(Rf_isNull(w) ? ones(n) : num(w, "w", n, 0.0));
  // sum(w[x1 > cut]) <= 0
  const double cv = REAL_RO(c)[0], *ap = REAL_RO(a), *wp = REAL_RO(wv);
  bool any = false;
  for (int i = 0; i < n && !any; ++i) any = ap[i] > cv && wp[i] > 0;
  if (!any) fail("no place with weight exceeds `cut`.");
  double r[7];
  kchk(clab_rtmhs(ap, REAL_RO(b), wp, n, cv, r), nullptr);
  static const char *const nm[] = {"nsel",   "change", "mirror", "sym",
                                   "excess", "low",    "lowsym"};
  UNPROTECT(5);
  return clab::nvec(7, nm, r);
}

// Corrections.

SEXP C_clab_r_incap(SEXP lam, SEXP q, SEXP s, SEXP share) {
  clab::who = "incap";
  static const char *const an[] = {"lam", "q", "s", "share"};
  SEXP a = PROTECT(group(args({lam, q, s, share}), an));
  lam = VECTOR_ELT(a, 1);
  q = VECTOR_ELT(a, 2);
  s = VECTOR_ELT(a, 3);
  share = VECTOR_ELT(a, 4);
  SEXP nms = PROTECT(VECTOR_ELT(a, 0));
  SEXP l = PROTECT(num(lam, "lam", -1, 0.0, true));
  const int n = ilen(static_cast<double>(XLENGTH(l)), "lam");
  SEXP qv = PROTECT(num(q, "q", -1, 0.0, true));
  const R_xlen_t nq = XLENGTH(qv);
  for (R_xlen_t i = 0; i < nq; ++i) {
    if (REAL_RO(qv)[i] > 1) fail("`q` must be at most 1.");
  }
  SEXP sv = PROTECT(num(s, "s", -1, 0.0));
  const R_xlen_t ns = XLENGTH(sv);
  if (nq != 1 && nq != n) fail("`q` must have length 1 or %d.", n);
  if (ns != 1 && ns != n) fail("`s` must have length 1 or %d.", n);
  const bool has = !Rf_isNull(share);
  SEXP sh = PROTECT(has ? num(share, "share", n, 0.0) : R_NilValue);
  // abs(sum(share) - 1) > 1e-9, the sum compensated
  if (has && std::fabs(clab::nsum(REAL_RO(sh), n) - 1) > 1e-9) {
    fail("`share` must sum to 1.");
  }
  // data.frame(group, lam, q = rep_len(q, n), s = rep_len(s, n), rate, prev,
  // marg), and with share the row rbind() adds: group "all", lam, rate and
  // prev of the whole, NA for the rest
  const R_xlen_t nr = has ? n + 1 : n;
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 7));
  // Without share the labels are the column as they stand; copying them
  // would turn the labels R defers for an unnamed lam into strings one by
  // one, which costs more than the whole computation.
  SEXP g = has ? Rf_allocVector(STRSXP, nr) : nms;
  SET_VECTOR_ELT(out, 0, g);
  if (has) {
    for (int i = 0; i < n; ++i) SET_STRING_ELT(g, i, STRING_ELT(nms, i));
    SET_STRING_ELT(g, n, Rf_mkChar("all"));
  }
  for (int k = 1; k < 7; ++k) SET_VECTOR_ELT(out, k, real(nr));
  double *c[7];
  for (int k = 1; k < 7; ++k) c[k] = REAL(VECTOR_ELT(out, k));
  for (int i = 0; i < n; ++i) {
    c[1][i] = REAL_RO(l)[i];
    c[2][i] = REAL_RO(qv)[nq == 1 ? 0 : i];
    c[3][i] = REAL_RO(sv)[ns == 1 ? 0 : i];
  }
  double agg[3];
  kchk(clab_incap(c[1], c[2], c[3], n, has ? REAL_RO(sh) : nullptr, c[4], c[5],
                  c[6], agg),
       nullptr);
  if (has) {
    c[1][n] = agg[0];
    c[2][n] = NA_REAL;
    c[3][n] = NA_REAL;
    c[4][n] = agg[1];
    c[5][n] = agg[2];
    c[6][n] = NA_REAL;
  }
  static const char *const nm[] = {"group", "lam",  "q",   "s",
                                   "rate",  "prev", "marg"};
  clab::frame(out, nm, has);
  UNPROTECT(7);
  return out;
}

SEXP C_clab_r_incapr(SEXP lam, SEXP t0, SEXP s, SEXP r) {
  clab::who = "incapr";
  SEXP l = PROTECT(num(lam, "lam", -1, 0.0));
  const R_xlen_t n = XLENGTH(l);
  const double *lp = REAL_RO(l);
  // any(diff(lam) > 0)
  for (R_xlen_t i = 1; i < n; ++i) {
    if (lp[i] - lp[i - 1] > 0) fail("`lam` must never increase.");
  }
  const int sv = clab::whole(s, "s", 1);
  const int tv = clab::whole(t0, "t0", 0);
  // t0 + s > length(lam) - 1, summed in double: t0 and s are each at most
  // INT_MAX, so the sum is exact and cannot overflow
  if (static_cast<double>(tv) + sv > static_cast<double>(n) - 1) {
    fail("`lam` must have length at least t0 + s + 1.");
  }
  const double rv = clab::p01(r, "r");
  double o[5];
  kchk(clab_incapr(lp, ilen(static_cast<double>(n), "lam"), tv, sv, rv, o), nullptr);
  static const char *const nm[] = {"prev", "net", "upper", "lower", "later"};
  UNPROTECT(1);
  return clab::nvec(5, nm, o);
}

SEXP C_clab_r_sellab(SEXP y, SEXP rel, SEXP rule, SEXP w) {
  clab::who = "sellab";
  static const char *const an[] = {"y", "rel", "rule", "w"};
  SEXP a = PROTECT(align(args({y, rel, rule, w}), an, false));
  y = VECTOR_ELT(a, 1);
  rel = VECTOR_ELT(a, 2);
  rule = VECTOR_ELT(a, 3);
  w = VECTOR_ELT(a, 4);
  const R_xlen_t n = Rf_xlength(rel);
  clab::lgl(rel, "rel", n);
  clab::lgl(rule, "rule", n);
  clab::unlab(y, "y");
  if (!clab::isnum(y) || XLENGTH(y) != n) {
    fail("`y` must be numbers, one per defendant.");
  }
  clab::vec(y, "y");
  SEXP yv = PROTECT(TYPEOF(y) == REALSXP ? dcopy(y) : clab::dbl(y));
  double *yp = REAL(yv);
  const int *rp = LOGICAL_RO(rel), *up = LOGICAL_RO(rule);
  // anyNA(y[rel]) || !all(y[rel] %in% c(0, 1)); y[!rel] <- 0
  for (R_xlen_t i = 0; i < n; ++i) {
    if (rp[i] && !(yp[i] == 0 || yp[i] == 1)) {
      fail("`y` must be 0 or 1 for every released defendant.");
    }
  }
  for (R_xlen_t i = 0; i < n; ++i) {
    if (!rp[i]) yp[i] = 0;
  }
  SEXP wv = PROTECT(Rf_isNull(w) ? ones(n) : num(w, "w", n, 0.0));
  if (!pos(REAL_RO(wv), n, rp)) fail("the judge released nobody.");
  if (!pos(REAL_RO(wv), n, up)) fail("the rule releases nobody.");
  double r[6];
  kchk(clab_sellab(yp, rp, up, REAL_RO(wv), ilen(static_cast<double>(n), "rel"), r), nullptr);
  static const char *const nm[] = {"obs",   "unobs", "lower",
                                   "upper", "width", "naive"};
  UNPROTECT(3);
  return clab::nvec(6, nm, r);
}

SEXP C_clab_r_losbnd(SEXP done, SEXP open) {
  clab::who = "losbnd";
  SEXP d = PROTECT(num(done, "done", -1, 0.0));
  SEXP o = PROTECT(num(open, "open", -1, 0.0));
  double r[4];
  kchk(clab_losbnd(REAL_RO(d), ilen(static_cast<double>(XLENGTH(d)), "done"), REAL_RO(o),
                   ilen(static_cast<double>(XLENGTH(o)), "open"), r),
       nullptr);
  static const char *const nm[] = {"mean", "age", "lower", "bias"};
  UNPROTECT(2);
  return clab::nvec(4, nm, r);
}

//' @srrstats {G2.4c} z and treat are converted as as.character() does.
//' @srrstats {G2.4e} A factor z is read by its labels; a factor given to a
//'   numeric argument fails isnum() and is refused.
//' @srrstats {G2.5} z is the one argument that may be a factor; its help
//'   page says level order and ordering are ignored.
//' @noRd
SEXP C_clab_r_sntbnd(SEXP y, SEXP z, SEXP treat, SEXP w) {
  clab::who = "sntbnd";
  static const char *const an[] = {"y", "z", "w"};
  SEXP a = PROTECT(align(args({y, z, w}), an, false));
  y = VECTOR_ELT(a, 1);
  z = VECTOR_ELT(a, 2);
  w = VECTOR_ELT(a, 3);
  SEXP yv = PROTECT(num(y, "y", -1, 0.0));
  const int n = ilen(static_cast<double>(XLENGTH(yv)), "y");
  const double *yp = REAL_RO(yv);
  for (int i = 0; i < n; ++i) {
    if (!(yp[i] == 0 || yp[i] == 1)) fail("`y` must be 0 or 1.");
  }
  if (Rf_xlength(z) != n || clab::anyna(z)) {
    fail("`z` must have length %d and no missing values.", n);
  }
  clab::vec(z, "z");
  clab::unlab(z, "z");
  if (!Rf_isNull(treat)) clab::unlab(treat, "treat");
  // z <- as.character(z); lv <- sort(unique(z)): z must take two values,
  // ordered by code point (not by the locale's collation)
  SEXP zc = PROTECT(Rf_isVectorAtomic(z) || TYPEOF(z) == VECSXP ? clab::labs(z) : R_NilValue);  // NOLINT
  if (Rf_isNull(zc)) {
    fail("`z` must be labels: a character, numeric or logical vector or a "
         "factor.");
  }
  SEXP lv[2] = {STRING_ELT(zc, 0), nullptr};
  for (int i = 1; i < n; ++i) {
    SEXP e = STRING_ELT(zc, i);
    if (clab::seql(e, lv[0])) continue;
    if (lv[1] == nullptr) {
      lv[1] = e;
    } else if (!clab::seql(e, lv[1])) {
      lv[1] = R_NilValue;
      break;
    }
  }
  if (lv[1] == nullptr || lv[1] == R_NilValue) {
    fail("`z` must take exactly two values.");
  }
  {
    const void *v = vmaxget();
    if (clab::kcmp(clab::key(lv[1]), clab::key(lv[0])) < 0) std::swap(lv[0], lv[1]);
    vmaxset(v);
  }
  // treat: one of lv, never chosen by sorting; by default TRUE for a
  // logical z and 1 for a numeric z of 0 and 1, otherwise required
  SEXP tr = lv[1];
  if (Rf_isNull(treat)) {
    const bool l = TYPEOF(z) == LGLSXP;
    const bool b = clab::isnum(z) && clab::kcmp(clab::key(lv[0]), Key{"0", 0}) == 0 &&
                   clab::kcmp(clab::key(lv[1]), Key{"1", 0}) == 0;
    if (!l && !b) {
      fail("`treat` must name the treatment, %s or %s: it is chosen by "
           "default only for a logical `z` (TRUE) or a `z` of 0 and 1 (1).",
           Rf_translateChar(lv[0]), Rf_translateChar(lv[1]));
    }
  }
  SEXP tc = PROTECT(!Rf_isNull(treat) && Rf_isVectorAtomic(treat) && XLENGTH(treat) == 1 ? clab::labs(treat) : R_NilValue);  // NOLINT
  if (!Rf_isNull(treat)) {
    // !is.atomic(treat) || length(treat) != 1L || !as.character(treat) %in% lv
    tr = Rf_isNull(tc) ? NA_STRING : STRING_ELT(tc, 0);
    if (!clab::seql(tr, lv[0]) && !clab::seql(tr, lv[1])) {
      fail("`treat` must be one of the values of `z`.");
    }
  }
  SEXP wv = PROTECT(Rf_isNull(w) ? ones(n) : num(w, "w", n, 0.0));
  if (!pos(REAL_RO(wv), n)) fail("`w` must have a positive total.");
  // as.integer(z == treat)
  int *zi = reinterpret_cast<int *>(R_alloc(static_cast<std::size_t>(n), sizeof(int)));
  for (int i = 0; i < n; ++i) zi[i] = clab::seql(STRING_ELT(zc, i), tr);
  double r[11];
  kchk(clab_sntbnd(yp, zi, REAL_RO(wv), n, r), nullptr);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
  SEXP term = Rf_allocVector(STRSXP, 3);
  SET_VECTOR_ELT(out, 0, term);
  SET_STRING_ELT(term, 0, clab::seql(lv[0], tr) ? lv[1] : lv[0]);
  SET_STRING_ELT(term, 1, tr);
  SET_STRING_ELT(term, 2, Rf_mkChar("effect"));
  const double col[4][3] = {{r[2], r[3], NA_REAL},
                            {r[4], r[6], r[8]},
                            {r[5], r[7], r[9]},
                            {r[0] / r[2], r[1] / r[3], r[10]}};
  for (int k = 0; k < 4; ++k) {
    SEXP v = Rf_allocVector(REALSXP, 3);
    SET_VECTOR_ELT(out, k + 1, v);
    for (int j = 0; j < 3; ++j) REAL(v)[j] = col[k][j];
  }
  static const char *const nm[] = {"term", "share", "lower", "upper",
                                   "observed"};
  clab::frame(out, nm);
  UNPROTECT(6);
  return out;
}

SEXP C_clab_r_hzdsel(SEXP s, SEXP h, SEXP l, SEXP surh, SEXP surl) {
  clab::who = "hzdsel";
  static const char *const an[] = {"surh", "surl"};
  SEXP al = PROTECT(group(args({surh, surl}), an));
  surh = VECTOR_ELT(al, 1);
  surl = VECTOR_ELT(al, 2);
  const double sv = clab::one(s, "s", 0, 1);
  const double hv = clab::p01(h, "h");
  const double lv = clab::p01(l, "l");
  // both are numbers in [0, 1], never NaN: this is the negation of l < h
  if (lv >= hv) fail("`l` must be below `h`.");
  SEXP a = PROTECT(num(surh, "surh", 2, 0.0, true));
  SEXP b = PROTECT(num(surl, "surl", 2, 0.0, true));
  for (int i = 0; i < 2; ++i) {
    if (REAL_RO(a)[i] > 1) fail("`surh` must be at most 1.");
  }
  for (int i = 0; i < 2; ++i) {
    if (REAL_RO(b)[i] > 1) fail("`surl` must be at most 1.");
  }
  // arms named control and treated (both arguments named) are read by
  // name, in either order; otherwise the first is the control arm
  SEXP g = VECTOR_ELT(al, 0);
  const bool sw = !Rf_isNull(rnames(surh)) && !Rf_isNull(rnames(surl)) &&
                  clab::kcmp(clab::key(STRING_ELT(g, 0)), Key{"treated", 0}) == 0 &&
                  clab::kcmp(clab::key(STRING_ELT(g, 1)), Key{"control", 0}) == 0;
  const double ah[2] = {REAL_RO(a)[sw], REAL_RO(a)[!sw]};
  const double bl[2] = {REAL_RO(b)[sw], REAL_RO(b)[!sw]};
  double r[5];
  kchk(clab_hzdsel(sv, hv, lv, ah, bl, r), nullptr);
  static const char *const nm[] = {"share0", "share1", "haz0", "haz1", "ratio"};
  UNPROTECT(3);
  return clab::nvec(5, nm, r);
}

// A hook for the tests, on lengths no test can allocate: with a name, the
// length check of the binding; with FALSE, that of the kernel entry
// points; with NULL, the length as an error message would print it.
SEXP C_clab_tlen(SEXP n, SEXP nm) {
  clab::who = "clab_tlen";
  const double v = REAL_RO(n)[0];
  if (Rf_isNull(nm)) {
    char b[48];
    clab::lfmt(v, b, sizeof b);
    return Rf_mkString(b);
  }
  const char *s = TYPEOF(nm) == STRSXP ? CHAR(STRING_ELT(nm, 0)) : nullptr;
  return Rf_ScalarInteger(ilen(v, s));
}

// The machine-checked theorems behind each function, by module of
// inst/lean, in the order proof() lists them.
namespace {
struct Thm {
  const char *fn, *mod, *thm;
};
const Thm thms[] = {
    {"vodtst", "Veil", "odds_ratio_identifies_bias"},
    {"vodtst", "Veil", "share_sign"},
    {"vodtst", "Veil", "exposure_shift_confounds"},
    {"orrisk", "RelRisk", "or_eq_rr_mul"},
    {"orrisk", "RelRisk", "rr_between"},
    {"orrisk", "RelRisk", "or_overstates"},
    {"outtst", "OutcomeTest", "infra_marginality"},
    {"outtst", "OutcomeTest", "equal_thresholds_unequal_hits"},
    {"dspbm", "Benchmark", "benchmark_product"},
    {"dspbm", "Benchmark", "benchmark_not_additive"},
    {"dyadrt", "Dyad", "null_slope_positive"},
    {"dyadrt", "Dyad", "pair_exposure_rate_constant"},
    {"gtmp", "Forensics", "bn_hom_zero"},
    {"gtmp", "Forensics", "bn_het_zero"},
    {"gtmp", "Forensics", "bn_hom_ge_product"},
    {"gtmp", "Forensics", "bn_het_not_ge_product"},
    {"gtmp", "Forensics", "bn_hom_mem"},
    {"gtmp", "Forensics", "bn_het_mem"},
    {"gtmp", "Forensics", "bn_hom_one"},
    {"gtmp", "Forensics", "bn_hom_tendsto_one"},
    {"lrpost", "Forensics", "post_eq"},
    {"lrpost", "Forensics", "post_lt_half"},
    {"cllr", "Forensics", "cllr_neutral"},
    {"cllr", "Forensics", "cllr_pos"},
    {"cllr", "Forensics", "cllr_perfect_tendsto"},
    {"dbsrch", "Forensics", "np_rule"},
    {"dbsrch", "Forensics", "any_match_eq"},
    {"dbsrch", "Forensics", "any_match_mem"},
    {"dbsrch", "Forensics", "any_match_kernel"},
    {"dbsrch", "Forensics", "balding_donnelly"},
    {"dbsrch", "Forensics", "post_ge_no_search"},
    {"mixlr", "Forensics", "mix_two_alleles"},
    {"mixlr", "Forensics", "mix_one_unknown_one_allele"},
    {"mixlr", "Forensics", "mix_one_unknown_two"},
    {"mixlr", "Forensics", "mix_one_unknown_three_zero"},
    {"mixlr", "Forensics", "mix_two_unknowns_four"},
    {"kinidx", "Forensics", "po_identical_het"},
    {"kinidx", "Forensics", "kin_mix"},
    {"kinidx", "Forensics", "r1_hom_hom"},
    {"kinidx", "Forensics", "r2_hom_hom"},
    {"kinidx", "Forensics", "r1_hom_het"},
    {"kinidx", "Forensics", "r1_het_hom"},
    {"kinidx", "Forensics", "r1_het_share_one"},
    {"kinidx", "Forensics", "r1_het_identical"},
    {"kinidx", "Forensics", "r2_het_identical"},
    {"kinidx", "Forensics", "r2_not_identical"},
    {"patidx", "Forensics", "post_even_prior"},
    {"patidx", "Forensics", "pat_hom_father_hom"},
    {"patidx", "Forensics", "pat_hom_father_het"},
    {"patidx", "Forensics", "pat_het_trio"},
    {"lrcal", "Calib", "push_pos"},
    {"lrcal", "Calib", "push_total"},
    {"lrcal", "Calib", "pav_increasing"},
    {"lrcal", "Calib", "above_iff"},
    {"lrcal", "Calib", "pool_weighted_mean"},
    {"lrcal", "Calib", "pool_between"},
    {"cllrmin", "Calib", "gibbs"},
    {"cllrmin", "Calib", "block_mean_optimal"},
    {"clrsmp", "Clearance", "agg_eq_weighted"},
    {"clrsmp", "Clearance", "kitagawa"},
    {"clrsmp", "Clearance", "change_ge_mix"},
    {"clrsmp", "Clearance", "fixed_mix_no_reversal"},
    {"clrsmp", "Clearance", "simpson_witness"},
    {"expbnd", "Exposure", "rate_bounds"},
    {"expbnd", "Exposure", "rate_lower_attained"},
    {"expbnd", "Exposure", "rate_upper_attained"},
    {"expbnd", "Exposure", "disparity_bounds"},
    {"expbnd", "Exposure", "disparity_sign_identified"},
    {"misrd", "Misclass", "binary_attenuation"},
    {"misrd", "Misclass", "binary_sign"},
    {"misrd", "Misclass", "factor_abs_le_one"},
    {"misrd", "Misclass", "measured_abs_le"},
    {"misrd", "Misclass", "binary_flip"},
    {"rtmhs", "RegMean", "selected_change_nonpos"},
    {"rtmhs", "RegMean", "low_selected_change_nonneg"},
    {"incap", "Incap", "steady_state_rate"},
    {"incap", "Incap", "prevented_share_eq"},
    {"incap", "Incap", "prevented_share_lt_one"},
    {"incap", "Incap", "rate_antitone_in_S"},
    {"incap", "Incap", "rate_antitone_in_q"},
    {"incap", "Incap", "marginal_prevention_eq"},
    {"incap", "Incap", "high_rate_more_prevented"},
    {"incap", "Incap", "kernel_marginal_eq"},
    {"incap", "Incap", "kernel_marginal_is"},
    {"incap", "Incap", "marginal_strict_anti"},
    {"incap", "Incap", "marginal_anti"},
    {"incapr", "Incap", "prevented_le_const"},
    {"incapr", "Incap", "prevented_ge_const"},
    {"incapr", "Incap", "later_sentence_prevents_less"},
    {"incapr", "Incap", "prevented_net_le"},
    {"sellab", "Select", "observed_rate_is_conditional"},
    {"sellab", "Select", "nested_rate_identified"},
    {"sellab", "Select", "unobserved_bounds"},
    {"sellab", "Select", "unobserved_width"},
    {"sellab", "Select", "unobserved_ends_attained"},
    {"sellab", "Select", "naive_between"},
    {"losbnd", "Censor", "true_mean_ge"},
    {"losbnd", "Censor", "lower_bound_sub"},
    {"losbnd", "Censor", "bias_lower"},
    {"losbnd", "Censor", "disposed_understates"},
    {"losbnd", "Censor", "no_upper_bound"},
    {"sntbnd", "Sentence", "outcome_bounds"},
    {"sntbnd", "Sentence", "lower_attained"},
    {"sntbnd", "Sentence", "upper_attained"},
    {"sntbnd", "Sentence", "ate_width_one"},
    {"sntbnd", "Sentence", "ate_contains_zero"},
    {"hzdsel", "Hazard", "survivor_hazard_mono"},
    {"hzdsel", "Hazard", "hr2_gt_one_of_depletion"},
    {"ntdeg", "Network", "handshake"},
    {"ntmod", "Network", "modularity_le_one"},
    {"ntmod", "Network", "modularity_ge_neg_half"},
    {"ntcore", "Network", "core_dense"},
    {"ntcore", "Network", "core_nested"},
    {"ntchg", "Network", "cheeger_lower"},
    {"relchk", "RelChk", "or_swap_ref"},
    {"relchk", "RelChk", "or_swap_event"},
    {"relchk", "RelChk", "or_swap_both"},
    {"relchk", "RelChk", "or_transpose"},
    {"relchk", "RelChk", "rr_swap_ref"},
    {"relchk", "RelChk", "rd_swap_ref"},
    {"relchk", "RelChk", "rd_swap_event"},
    {"cntreg", "Rates", "rr_swap"},
    {"cntreg", "Rates", "ci_swap"},
    {"grprr", "Rates", "rd_swap"},
    {"dastd", "Rates", "dsr_le_max"},
    {"dastd", "Rates", "dsr_ge_min"},
    {"sirex", "Rates", "eb_convex"},
    {"cgstat", "Spatial", "mean_centre_minimises"},
    {"sdell", "Spatial", "parallel_axis"},
    {"lmorn", "Spatial", "lisa_sum"},
    {"gearc", "Spatial", "geary_mean_one"},
    {"gorg", "Spatial", "getis_mean"},
    {"dirmn", "Spatial", "resultant_le_one"},
};
constexpr int nthm = sizeof thms / sizeof thms[0];
}  // namespace

SEXP C_clab_r_proof(SEXP fn) {
  clab::who = "proof";
  // if (!is.character(fn) || anyNA(fn)) stop("`fn` must be function names.")
  const bool all = Rf_isNull(fn);
  if (!all && (TYPEOF(fn) != STRSXP || clab::anyna(fn))) {
    fail("`fn` must be function names.");
  }
  Keys fk{};
  if (!all) {
    fk = keys(fn);
    // miss <- setdiff(fn, t[["fn"]]): the first appearance of each unknown
    // name, in the order given
    char *miss = R_alloc(static_cast<std::size_t>(fk.n), 1);
    std::size_t len = 0;
    for (R_xlen_t i = 0; i < fk.n; ++i) {
      bool known = !fk.first[i];
      for (int t = 0; t < nthm && !known; ++t) {
        known = clab::kcmp(fk.k[i], Key{thms[t].fn, 0}) == 0;
      }
      miss[i] = !known;
      if (!known) len += std::strlen(Rf_translateChar(STRING_ELT(fn, i))) + 2;
    }
    if (len > 0) {
      // toString(miss): joined by ", "
      char *msg = R_alloc(len + 1, 1), *e = msg;
      for (R_xlen_t i = 0; i < fk.n; ++i) {
        if (!miss[i]) continue;
        if (e != msg) {
          std::memcpy(e, ", ", 2);
          e += 2;
        }
        const char *c = Rf_translateChar(STRING_ELT(fn, i));
        const std::size_t m = std::strlen(c);
        std::memcpy(e, c, m);
        e += m;
      }
      *e = '\0';
      fail("no proof recorded for: %s.", msg);
    }
  }
  // t[t[["fn"]] %in% fn, , drop = FALSE], row names reset
  int k = 0;
  bool keep[nthm];
  for (int t = 0; t < nthm; ++t) {
    keep[t] = all || has(fk, Key{thms[t].fn, 0});
    k += keep[t];
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  for (int c = 0; c < 3; ++c) SET_VECTOR_ELT(out, c, Rf_allocVector(STRSXP, k));
  char b[128];
  for (int t = 0, i = 0; t < nthm; ++t) {
    if (!keep[t]) continue;
    SET_STRING_ELT(VECTOR_ELT(out, 0), i, Rf_mkChar(thms[t].fn));
    std::snprintf(b, sizeof b, "Crimlab.%s.%s", thms[t].mod, thms[t].thm);
    SET_STRING_ELT(VECTOR_ELT(out, 1), i, Rf_mkChar(b));
    std::snprintf(b, sizeof b, "lean/%s.lean", thms[t].mod);
    SET_STRING_ELT(VECTOR_ELT(out, 2), i, Rf_mkChar(b));
    ++i;
  }
  static const char *const nm[] = {"fn", "theorem", "file"};
  // rbind() and `rownames<-` leave row.names before class
  clab::frame(out, nm, true);
  UNPROTECT(1);
  return out;
}

}  // extern "C"

// Networks.

namespace {

// node labels: as.character() of an atomic vector (a factor by its labels),
// no NA; Date, POSIXt and difftime refused as elsewhere
SEXP nlab(SEXP x, const char *nm) {
  SEXP s = R_NilValue;
  if (Rf_isVectorAtomic(x)) {
    clab::vec(x, nm);
    s = clab::labs(x);
  }
  if (Rf_isNull(s) || clab::anyna(s)) {
    fail("`%s` must be labels (a vector) with no missing values.", nm);
  }
  return s;
}

// codes 0, 1, ... of the strings of v by first appearance; nd the number
// of distinct strings, first[c] the position of code c's first appearance
int *fcode(SEXP v, int *nd, R_xlen_t **first) {
  const Keys s = keys(v);
  const R_xlen_t n = s.n;
  int *c = reinterpret_cast<int *>(R_alloc(static_cast<std::size_t>(n), sizeof(int)));
  R_xlen_t *f = reinterpret_cast<R_xlen_t *>(R_alloc(static_cast<std::size_t>(s.nd), sizeof(R_xlen_t)));
  int k = 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    if (s.first[i]) {
      f[k] = i;
      c[i] = k++;
    } else {
      c[i] = c[s.slot[probe(s, s.k[i])]];
    }
  }
  *nd = k;
  if (first) *first = f;
  return c;
}

// a graph read from its labels: n nodes (labels lab, left PROTECTed), m
// edges a[k] - b[k] by node code
struct Net {
  SEXP lab;
  int n, m;
  int *a, *b;
};

//' @srrstats {G2.4c} Node, offender, incident and group labels are
//'   converted as as.character() does.
//' @srrstats {G2.4e} A factor given as labels is read by its labels, never
//'   its codes; levels that do not occur are not nodes.
//' @noRd
Net net(SEXP from, SEXP to, SEXP nodes) {
  Net g{R_NilValue, 0, 0, nullptr, nullptr};
  SEXP f = PROTECT(nlab(from, "from"));
  SEXP t = PROTECT(nlab(to, "to"));
  const R_xlen_t m = XLENGTH(f);
  if (XLENGTH(t) != m) {
    char x[32], y[32];
    clab::lfmt(static_cast<double>(m), x, sizeof x);
    clab::lfmt(static_cast<double>(XLENGTH(t)), y, sizeof y);
    fail("`to` must have the length of `from`, %s, not %s.", x, y);
  }
  ilen(2.0 * m, "from");
  g.m = static_cast<int>(m);
  g.a = reinterpret_cast<int *>(R_alloc(static_cast<std::size_t>(m), sizeof(int)));
  g.b = reinterpret_cast<int *>(R_alloc(static_cast<std::size_t>(m), sizeof(int)));
  SEXP lab;
  if (Rf_isNull(nodes)) {
    // first appearance, reading each edge's from, then its to
    SEXP v = PROTECT(Rf_allocVector(STRSXP, 2 * m));
    for (R_xlen_t k = 0; k < m; ++k) {
      SET_STRING_ELT(v, 2 * k, STRING_ELT(f, k));
      SET_STRING_ELT(v, 2 * k + 1, STRING_ELT(t, k));
    }
    R_xlen_t *fst;
    const int *c = fcode(v, &g.n, &fst);
    for (R_xlen_t k = 0; k < m; ++k) {
      g.a[k] = c[2 * k];
      g.b[k] = c[2 * k + 1];
    }
    lab = Rf_allocVector(STRSXP, g.n);
    for (int i = 0; i < g.n; ++i) SET_STRING_ELT(lab, i, STRING_ELT(v, fst[i]));
    UNPROTECT(3);
  } else {
    SEXP nd = PROTECT(nlab(nodes, "nodes"));
    g.n = ilen(static_cast<double>(XLENGTH(nd)), "nodes");
    const Keys s = keys(nd);
    for (R_xlen_t i = 0; i < s.n; ++i) {
      if (!s.first[i]) {
        fail("`nodes` must not repeat a label: %s.",
             Rf_translateChar(STRING_ELT(nd, i)));
      }
    }
    for (int e = 0; e < 2; ++e) {
      SEXP x = e ? t : f;
      int *o = e ? g.b : g.a;
      for (R_xlen_t k = 0; k < m; ++k) {
        const R_xlen_t j = s.n > 0 ? s.slot[probe(s, clab::key(STRING_ELT(x, k)))] : -1;
        if (j < 0) {
          fail("`%s` names a node not in `nodes`: %s.", e ? "to" : "from",
               Rf_translateChar(STRING_ELT(x, k)));
        }
        o[k] = static_cast<int>(j);
      }
    }
    // a fresh vector: the labels without the attributes `nodes` had
    lab = Rf_allocVector(STRSXP, g.n);
    for (int i = 0; i < g.n; ++i) SET_STRING_ELT(lab, i, STRING_ELT(nd, i));
    UNPROTECT(3);
  }
  if (g.n == 0) fail("the graph has no nodes: give edges or `nodes`.");
  g.lab = PROTECT(lab);
  return g;
}

// TRUE or FALSE
bool flag(SEXP x, const char *nm) {
  if (TYPEOF(x) != LGLSXP || XLENGTH(x) != 1 || LOGICAL_RO(x)[0] == NA_LOGICAL) {
    fail("`%s` must be TRUE or FALSE.", nm);
  }
  return LOGICAL_RO(x)[0] != 0;
}

// the edge weights: NULL, or one positive number per edge
const double *wts(SEXP w, const Net &g) {
  if (Rf_isNull(w)) return nullptr;
  // the copy lives in R's transient memory to the end of the call; it is
  // made before num(), whose result is read at once
  double *p = reinterpret_cast<double *>(R_alloc(static_cast<std::size_t>(g.m), sizeof(double)));
  SEXP v = num(w, "w", g.m, 0.0, true);
  for (int k = 0; k < g.m; ++k) p[k] = REAL_RO(v)[k];
  return p;
}

// a graph in compressed sparse rows (crimlab.h); rev reverses each edge
struct Csr {
  int *off, *adj;
  double *wt;
};
Csr csr(const Net &g, const double *w, bool dir, bool rev = false) {
  const std::size_t cap = static_cast<std::size_t>(g.m) * (dir ? 1 : 2);
  Csr c;
  c.off = reinterpret_cast<int *>(R_alloc(static_cast<std::size_t>(g.n) + 1, sizeof(int)));
  c.adj = reinterpret_cast<int *>(R_alloc(cap, sizeof(int)));
  c.wt = reinterpret_cast<double *>(R_alloc(cap, sizeof(double)));
  int *iw = reinterpret_cast<int *>(R_alloc(static_cast<std::size_t>(g.n), sizeof(int)));
  int nl, nd;
  kchk(clab_ntcsr(rev ? g.b : g.a, rev ? g.a : g.b, w, g.m, g.n, dir, c.off,
                  c.adj, c.wt, iw, &nl, &nd),
       nullptr);
  return c;
}

int *ibuf(double n) {
  return reinterpret_cast<int *>(R_alloc(static_cast<std::size_t>(n), sizeof(int)));
}
double *dbuf(double n) {
  return reinterpret_cast<double *>(R_alloc(static_cast<std::size_t>(n), sizeof(double)));
}

// one value per node of x, coded 0, 1, ... by first appearance in node
// order: from names(x) matched to the labels, or, unnamed, x in node order.
// The values as text go to *val (left PROTECTed).
int *pernode(SEXP x, const char *nm, const Net &g, SEXP *val) {
  SEXP s = PROTECT(nlab(x, nm));
  SEXP nms = Rf_getAttrib(x, R_NamesSymbol);
  SEXP v = PROTECT(Rf_allocVector(STRSXP, g.n));
  if (Rf_isNull(nms)) {
    if (XLENGTH(s) != g.n) {
      fail("`%s` must have one value per node (%d) or names.", nm, g.n);
    }
    for (int i = 0; i < g.n; ++i) SET_STRING_ELT(v, i, STRING_ELT(s, i));
  } else {
    const Keys k = keys(nms);
    for (R_xlen_t i = 0; i < k.n; ++i) {
      if (!k.first[i]) {
        fail("`%s` must not repeat a name: %s.", nm,
             Rf_translateChar(STRING_ELT(nms, i)));
      }
    }
    for (int i = 0; i < g.n; ++i) {
      const R_xlen_t j = k.n > 0 ? k.slot[probe(k, clab::key(STRING_ELT(g.lab, i)))] : -1;
      if (j < 0) {
        fail("`%s` has no value for node %s.", nm,
             Rf_translateChar(STRING_ELT(g.lab, i)));
      }
      SET_STRING_ELT(v, i, STRING_ELT(s, j));
    }
  }
  int nd;
  int *c = fcode(v, &nd, nullptr);
  UNPROTECT(2);
  *val = PROTECT(v);
  return c;
}

// a data frame: node, then the columns in cols (k of them)
SEXP nframe(const Net &g, SEXP *cols, int k, const char *const *nm) {
  SEXP out = PROTECT(Rf_allocVector(VECSXP, k + 1));
  SET_VECTOR_ELT(out, 0, g.lab);
  for (int j = 0; j < k; ++j) SET_VECTOR_ELT(out, j + 1, cols[j]);
  clab::frame(out, nm);
  UNPROTECT(1);
  return out;
}

SEXP ints(int n, const int *x, int add = 0) {
  SEXP v = Rf_allocVector(INTSXP, n);
  for (int i = 0; i < n; ++i) INTEGER(v)[i] = x[i] + add;
  return v;
}

}  // namespace

extern "C" {

SEXP C_clab_r_ntdeg(SEXP from, SEXP to, SEXP w, SEXP nodes, SEXP dir) {
  clab::who = "ntdeg";
  const Net g = net(from, to, nodes);
  const double *wp = wts(w, g);
  const bool d = flag(dir, "dir");
  const Csr o = csr(g, wp, d);
  SEXP c[6];
  for (int j = 0; j < 6; ++j) c[j] = R_NilValue;
  SEXP od = PROTECT(Rf_allocVector(INTSXP, g.n));
  SEXP os = PROTECT(real(g.n));
  kchk(clab_ntdeg(g.n, o.off, o.wt, INTEGER(od), REAL(os)), nullptr);
  if (!d) {
    c[0] = od;
    c[1] = os;
    static const char *const nm[] = {"node", "deg", "wdeg"};
    SEXP r = nframe(g, c, 2, nm);
    UNPROTECT(3);
    return r;
  }
  const Csr ic = csr(g, wp, true, true);
  SEXP id = PROTECT(Rf_allocVector(INTSXP, g.n));
  SEXP is = PROTECT(real(g.n));
  kchk(clab_ntdeg(g.n, ic.off, ic.wt, INTEGER(id), REAL(is)), nullptr);
  SEXP td = PROTECT(Rf_allocVector(INTSXP, g.n));
  SEXP ts = PROTECT(real(g.n));
  for (int i = 0; i < g.n; ++i) {
    INTEGER(td)[i] = INTEGER(id)[i] + INTEGER(od)[i];
    REAL(ts)[i] = REAL(is)[i] + REAL(os)[i];
  }
  SEXP cc[6] = {id, od, td, is, os, ts};
  static const char *const nm[] = {"node", "in", "out", "deg", "win", "wout", "wdeg"};
  SEXP r = nframe(g, cc, 6, nm);
  UNPROTECT(7);
  return r;
}

SEXP C_clab_r_ntbtw(SEXP from, SEXP to, SEXP nodes, SEXP dir, SEXP norm) {
  clab::who = "ntbtw";
  const Net g = net(from, to, nodes);
  const bool d = flag(dir, "dir"), nz = flag(norm, "norm");
  const Csr o = csr(g, nullptr, d);
  SEXP b = PROTECT(real(g.n));
  kchk(clab_ntbtw(g.n, o.off, o.adj, d, nz, REAL(b), dbuf(2.0 * g.n),
                  ibuf(2.0 * g.n)),
       nullptr);
  static const char *const nm[] = {"node", "btw"};
  SEXP r = nframe(g, &b, 1, nm);
  UNPROTECT(2);
  return r;
}

SEXP C_clab_r_ntcls(SEXP from, SEXP to, SEXP nodes, SEXP dir) {
  clab::who = "ntcls";
  const Net g = net(from, to, nodes);
  const bool d = flag(dir, "dir");
  const Csr o = csr(g, nullptr, d);
  SEXP c[4];
  c[0] = PROTECT(Rf_allocVector(INTSXP, g.n));
  for (int j = 1; j < 4; ++j) c[j] = PROTECT(real(g.n));
  kchk(clab_ntcls(g.n, o.off, o.adj, REAL(c[1]), REAL(c[2]), REAL(c[3]),
                  INTEGER(c[0]), ibuf(2.0 * g.n)),
       nullptr);
  static const char *const nm[] = {"node", "reach", "cls", "ncls", "hrm"};
  SEXP r = nframe(g, c, 4, nm);
  UNPROTECT(5);
  return r;
}

SEXP C_clab_r_ntevc(SEXP from, SEXP to, SEXP w, SEXP nodes, SEXP tol,
                    SEXP maxit) {
  clab::who = "ntevc";
  const Net g = net(from, to, nodes);
  const double *wp = wts(w, g);
  const double tl = clab::one(tol, "tol", 0, 1);
  const int mi = clab::whole(maxit, "maxit", 1);
  const Csr o = csr(g, wp, false);
  SEXP x = PROTECT(real(g.n));
  double val;
  int it, cv;
  kchk(clab_ntevc(g.n, o.off, o.adj, o.wt, tl, mi, REAL(x), &val, &it, &cv,
                  dbuf(4.0 * g.n), ibuf(2.0 * g.n)),
       nullptr);
  static const char *const fn[] = {"node", "evc"};
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out, 0, nframe(g, &x, 1, fn));
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(val));
  SET_VECTOR_ELT(out, 2, Rf_ScalarInteger(it));
  SET_VECTOR_ELT(out, 3, Rf_ScalarLogical(cv));
  static const char *const nm[] = {"cent", "value", "iter", "conv"};
  clab::nlist(out, nm);
  UNPROTECT(3);
  return out;
}

SEXP C_clab_r_ntcore(SEXP from, SEXP to, SEXP nodes) {
  clab::who = "ntcore";
  const Net g = net(from, to, nodes);
  const Csr o = csr(g, nullptr, false);
  SEXP c = PROTECT(Rf_allocVector(INTSXP, g.n));
  kchk(clab_ntcore(g.n, o.off, o.adj, INTEGER(c), ibuf(3.0 * g.n)), nullptr);
  static const char *const nm[] = {"node", "core"};
  SEXP r = nframe(g, &c, 1, nm);
  UNPROTECT(2);
  return r;
}

SEXP C_clab_r_ntcmp(SEXP from, SEXP to, SEXP nodes, SEXP dir) {
  clab::who = "ntcmp";
  const Net g = net(from, to, nodes);
  const bool d = flag(dir, "dir");
  const Csr o = csr(g, nullptr, d);
  const Csr ic = d ? csr(g, nullptr, true, true) : Csr{nullptr, nullptr, nullptr};
  int *wc = ibuf(g.n), *sc = ibuf(g.n), *iw = ibuf(6.0 * g.n), nc;
  kchk(clab_ntwcc(g.n, o.off, o.adj, ic.off, ic.adj, wc, &nc, iw), nullptr);
  kchk(clab_ntscc(g.n, o.off, o.adj, sc, &nc, iw), nullptr);
  SEXP c[2];
  c[0] = PROTECT(ints(g.n, wc, 1));
  c[1] = PROTECT(ints(g.n, sc, 1));
  static const char *const nm[] = {"node", "weak", "strong"};
  SEXP r = nframe(g, c, 2, nm);
  UNPROTECT(3);
  return r;
}

SEXP C_clab_r_ntlvn(SEXP from, SEXP to, SEXP w, SEXP nodes, SEXP res,
                    SEXP seed) {
  clab::who = "ntlvn";
  const Net g = net(from, to, nodes);
  const double *wp = wts(w, g);
  const double rs = REAL_RO(num(res, "res", 1, 0.0, true))[0];
  const bool use = !Rf_isNull(seed);
  const unsigned sd = use ? static_cast<unsigned>(clab::whole(seed, "seed", 0)) : 0u;
  const Csr o = csr(g, wp, false);
  double len, il, q;
  kchk(clab_ntlvnw(g.n, o.off[g.n], &len, &il), nullptr);
  int *m = ibuf(g.n), lev;
  kchk(clab_ntlvn(g.n, o.off, o.adj, o.wt, rs, sd, use, m, &q, &lev, dbuf(len),
                  ibuf(il)),
       nullptr);
  SEXP c = PROTECT(ints(g.n, m, 1));
  static const char *const fn[] = {"node", "comm"};
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, nframe(g, &c, 1, fn));
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(q));
  SET_VECTOR_ELT(out, 2, Rf_ScalarInteger(lev));
  static const char *const nm[] = {"memb", "q", "levels"};
  clab::nlist(out, nm);
  UNPROTECT(3);
  return out;
}

SEXP C_clab_r_ntmod(SEXP from, SEXP to, SEXP memb, SEXP w, SEXP nodes,
                    SEXP res) {
  clab::who = "ntmod";
  const Net g = net(from, to, nodes);
  SEXP val;
  const int *m = pernode(memb, "memb", g, &val);
  const double *wp = wts(w, g);
  const double rs = REAL_RO(num(res, "res", 1, 0.0))[0];
  const Csr o = csr(g, wp, false);
  double q;
  kchk(clab_ntmod(g.n, o.off, o.adj, o.wt, m, rs, &q, dbuf(2.0 * g.n)), nullptr);
  UNPROTECT(2);
  return Rf_ScalarReal(q);
}

SEXP C_clab_r_ntbrk(SEXP from, SEXP to, SEXP grp, SEXP nodes, SEXP dir) {
  clab::who = "ntbrk";
  const Net g = net(from, to, nodes);
  SEXP val;
  const int *gr = pernode(grp, "grp", g, &val);
  const bool d = flag(dir, "dir");
  const Csr o = csr(g, nullptr, d);
  const Csr ic = d ? csr(g, nullptr, true, true) : o;
  double *cnt = dbuf(5.0 * g.n);
  kchk(clab_ntbrk(g.n, o.off, o.adj, ic.off, ic.adj, gr, cnt, ibuf(g.n)),
       nullptr);
  SEXP c[7];
  c[0] = val;
  for (int j = 1; j < 7; ++j) c[j] = PROTECT(real(g.n));
  for (int i = 0; i < g.n; ++i) {
    double t = 0.0;
    for (int r = 0; r < 5; ++r) {
      REAL(c[r + 1])[i] = cnt[r * g.n + i];
      t += cnt[r * g.n + i];
    }
    REAL(c[6])[i] = t;
  }
  static const char *const nm[] = {"node",  "grp",  "coord", "itin",
                                   "repr",  "gate", "liais", "total"};
  SEXP r = nframe(g, c, 7, nm);
  UNPROTECT(8);
  return r;
}

SEXP C_clab_r_ntchg(SEXP from, SEXP to, SEXP w, SEXP nodes, SEXP tol,
                    SEXP maxit) {
  clab::who = "ntchg";
  const Net g = net(from, to, nodes);
  const double *wp = wts(w, g);
  const double tl = clab::one(tol, "tol", 0, 1);
  const int mi = clab::whole(maxit, "maxit", 1);
  if (g.n < 2) fail("the graph needs at least two nodes.");
  const Csr o = csr(g, wp, false);
  int *iw = ibuf(g.n), *set = ibuf(g.n), nc, it, cv;
  kchk(clab_ntwcc(g.n, o.off, o.adj, nullptr, nullptr, set, &nc, iw), nullptr);
  if (nc > 1) {
    fail("the graph must be connected, not in %d components (see ntcmp()).",
         nc);
  }
  double b[4];
  kchk(clab_ntchg(g.n, o.off, o.adj, o.wt, tl, mi, b, set, &it, &cv,
                  dbuf(4.0 * g.n), iw),
       nullptr);
  int k = 0;
  for (int i = 0; i < g.n; ++i) k += set[i];
  SEXP s = PROTECT(Rf_allocVector(STRSXP, k));
  for (int i = 0, j = 0; i < g.n; ++i) {
    if (set[i]) SET_STRING_ELT(s, j++, STRING_ELT(g.lab, i));
  }
  static const char *const bn[] = {"lambda2", "lower", "upper", "sweep"};
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out, 0, clab::nvec(4, bn, b));
  SET_VECTOR_ELT(out, 1, s);
  SET_VECTOR_ELT(out, 2, Rf_ScalarInteger(it));
  SET_VECTOR_ELT(out, 3, Rf_ScalarLogical(cv));
  static const char *const nm[] = {"bound", "set", "iter", "conv"};
  clab::nlist(out, nm);
  UNPROTECT(3);
  return out;
}

SEXP C_clab_r_coofg(SEXP inc, SEXP off) {
  clab::who = "coofg";
  SEXP si = PROTECT(nlab(inc, "inc"));
  SEXP so = PROTECT(nlab(off, "off"));
  const R_xlen_t k = XLENGTH(si);
  if (XLENGTH(so) != k) {
    char x[32], y[32];
    clab::lfmt(static_cast<double>(k), x, sizeof x);
    clab::lfmt(static_cast<double>(XLENGTH(so)), y, sizeof y);
    fail("`off` must have the length of `inc`, %s, not %s.", x, y);
  }
  if (k == 0) fail("`inc` must name at least one incident.");
  ilen(2.0 * k, "inc");
  int ni, no;
  R_xlen_t *fo;
  const int *ci = fcode(si, &ni, nullptr);
  const int *co = fcode(so, &no, &fo);
  double np, il;
  int *iw = ibuf(static_cast<double>(ni) + 1 + k + no);
  kchk(clab_coofw(ci, co, static_cast<int>(k), ni, no, &np, &il, iw), nullptr);
  const int cap = static_cast<int>(np);
  int *pa = ibuf(cap), *pb = ibuf(cap), n;
  double *pw = dbuf(cap);
  kchk(clab_coofg(ci, co, static_cast<int>(k), ni, no, cap, pa, pb, pw, &n,
                  ibuf(il)),
       nullptr);
  SEXP e = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP f = Rf_allocVector(STRSXP, n);
  SET_VECTOR_ELT(e, 0, f);
  SEXP t = Rf_allocVector(STRSXP, n);
  SET_VECTOR_ELT(e, 1, t);
  SEXP w = Rf_allocVector(REALSXP, n);
  SET_VECTOR_ELT(e, 2, w);
  for (int p = 0; p < n; ++p) {
    SET_STRING_ELT(f, p, STRING_ELT(so, fo[pa[p]]));
    SET_STRING_ELT(t, p, STRING_ELT(so, fo[pb[p]]));
    REAL(w)[p] = pw[p];
  }
  static const char *const en[] = {"from", "to", "w"};
  clab::frame(e, en);
  SEXP nd = PROTECT(Rf_allocVector(STRSXP, no));
  for (int j = 0; j < no; ++j) SET_STRING_ELT(nd, j, STRING_ELT(so, fo[j]));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, e);
  SET_VECTOR_ELT(out, 1, nd);
  static const char *const nm[] = {"edges", "nodes"};
  clab::nlist(out, nm);
  UNPROTECT(5);
  return out;
}

}  // extern "C"

// Checks.

namespace {

// One side of relchk()'s table: lab when it is character, with no label
// missing or given twice (either is an error naming the side); "1", "2",
// ... when there is none (the integer row names of a data frame count as
// none).
SEXP tlab(SEXP lab, R_xlen_t n, const char *side) {
  if (TYPEOF(lab) != STRSXP || XLENGTH(lab) != n) return clab::seqs(n);
  const Keys k = keys(lab);
  for (R_xlen_t i = 0; i < n; ++i) {
    if (STRING_ELT(lab, i) == NA_STRING) fail("`x` has a missing %s label.", side);
    if (!k.first[i]) {
      fail("`x` names a %s twice: %s.", side, Rf_translateChar(STRING_ELT(lab, i)));
    }
  }
  return lab;
}

const char *const tforms =
    "a matrix, a 2-d table, a data frame of numeric columns or a list of "
    "numeric vectors of one length";

// element i of an integer or double vector, as a double (NA as NA)
double celd(SEXP s, R_xlen_t i) {
  if (TYPEOF(s) == REALSXP) return REAL_RO(s)[i];
  const int v = INTEGER_RO(s)[i];
  return v == NA_INTEGER ? NA_REAL : static_cast<double>(v);
}

//' @srrstats {G2.7} relchk() takes its table as a matrix, a 2-d table, a
//'   data frame or a list of columns; catchk() takes a vector, a data
//'   frame or a list of vectors.
//' @srrstats {G2.8} Each form becomes one double matrix of counts with its
//'   row and column labels before the kernel sees it.
//' @srrstats {G2.9} No information is lost in the conversion: the labels
//'   are kept, integer row names of a data frame become "1", "2", ... as
//'   R prints them, and nothing is added.
//' @srrstats {G2.10} Columns are read by position from the list that holds
//'   them, the same way for a data frame and a plain list.
//' @srrstats {G2.11} A column with a class is read as its numbers when
//'   is.numeric() is TRUE for it and refused otherwise (a factor, a Date).
//' @srrstats {G2.12} A list column is refused by relchk(); catchk()
//'   reports it as type "list" and does not audit it.
//' @srrstats {G2.13} A missing count is an error, never dropped.
//' @noRd
// relchk()'s table: list(counts, row labels, column labels), the counts as
// doubles column-major, r rows by c columns, every check of relchk()'s
// help page made: at least 2 x 2, no NA, finite, non-negative, labels
// neither missing nor repeated, no row or column without a count.
SEXP ctab(SEXP x, int *r, int *c) {
  clab::unlab(x, "x");
  R_xlen_t nr = 0, nc = 0;
  SEXP rn = R_NilValue, cn = R_NilValue;
  const bool lst = TYPEOF(x) == VECSXP;
  if (lst) {
    nc = XLENGTH(x);
    for (R_xlen_t j = 0; j < nc; ++j) {
      SEXP e = VECTOR_ELT(x, j);
      clab::unlab(e, "x");
      if (!clab::isnum(e) || Rf_isMatrix(e) || (j > 0 && XLENGTH(e) != nr)) {
        fail("`x` must be counts: %s.", tforms);
      }
      if (j == 0) nr = XLENGTH(e);
    }
    cn = Rf_getAttrib(x, R_NamesSymbol);
    if (Rf_inherits(x, "data.frame")) {
      rn = Rf_getAttrib(x, R_RowNamesSymbol);
    } else if (nc > 0) {
      rn = Rf_getAttrib(VECTOR_ELT(x, 0), R_NamesSymbol);
    }
  } else if (Rf_isMatrix(x) && clab::isnum(x)) {
    nr = Rf_nrows(x);
    nc = Rf_ncols(x);
    SEXP dn = Rf_getAttrib(x, R_DimNamesSymbol);
    if (!Rf_isNull(dn)) {
      rn = VECTOR_ELT(dn, 0);
      cn = VECTOR_ELT(dn, 1);
    }
  } else {
    fail("`x` must be counts: %s.", tforms);
  }
  PROTECT(rn);
  if (nr < 2 || nc < 2) fail("`x` must have at least two rows and two columns.");
  // 3 r c (r + c - 2) readings, each an int index
  if (3.0 * static_cast<double>(nr) * static_cast<double>(nc) *
          (static_cast<double>(nr) + static_cast<double>(nc) - 2) > INT_MAX) {
    fail("`x` is too large: its readings would outnumber the largest integer.");
  }
  const int ir = static_cast<int>(nr), ic = static_cast<int>(nc);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP v = real(nr * nc);
  SET_VECTOR_ELT(out, 0, v);
  double *p = REAL(v);
  for (R_xlen_t j = 0; j < nc; ++j) {
    for (R_xlen_t i = 0; i < nr; ++i) {
      p[i + nr * j] = lst ? celd(VECTOR_ELT(x, j), i) : celd(x, i + nr * j);
    }
  }
  for (R_xlen_t k = 0; k < nr * nc; ++k) {
    if (std::isnan(p[k])) fail("`x` has a missing count: the table must be complete.");
    if (!std::isfinite(p[k])) fail("`x` must be finite counts.");
    if (p[k] < 0) fail("`x` has a negative count.");
  }
  SET_VECTOR_ELT(out, 1, tlab(rn, nr, "row"));
  SET_VECTOR_ELT(out, 2, tlab(cn, nc, "column"));
  for (R_xlen_t i = 0; i < nr; ++i) {
    double s = 0;
    for (R_xlen_t j = 0; j < nc; ++j) s += p[i + nr * j];
    if (s == 0) {
      fail("`x` has no count in row %s: a group with no units has no contrast.",
           Rf_translateChar(STRING_ELT(VECTOR_ELT(out, 1), i)));
    }
  }
  for (R_xlen_t j = 0; j < nc; ++j) {
    double s = 0;
    for (R_xlen_t i = 0; i < nr; ++i) s += p[i + nr * j];
    if (s == 0) {
      fail("`x` has no count in column %s: an outcome no unit has gives no "
           "contrast.", Rf_translateChar(STRING_ELT(VECTOR_ELT(out, 2), j)));
    }
  }
  *r = ir;
  *c = ic;
  UNPROTECT(2);
  return out;
}

// The fewest decimals, 0 to 15, at which v prints as itself; -1 if none.
int ndec(double v) {
  char b[400];
  for (int d = 0; d <= 15; ++d) {
    std::snprintf(b, sizeof b, "%.*f", d, v);
    if (std::strtod(b, nullptr) == v) return d;
  }
  return -1;
}

const char *const kinds[] = {"or", "rr", "rd"};

}  // namespace

extern "C" {

SEXP C_clab_r_relchk(SEXP x, SEXP value, SEXP type, SEXP group, SEXP ref,
                     SEXP event, SEXP digits, SEXP tol) {
  clab::who = "relchk";
  int r, c;
  SEXP tb = PROTECT(ctab(x, &r, &c));
  SEXP vv = PROTECT(num(value, "value", 1));
  const double val = REAL_RO(vv)[0];
  clab::unlab(type, "type");
  int q = -1;
  if (TYPEOF(type) == STRSXP && XLENGTH(type) == 1) {
    for (int i = 0; i < 3; ++i) {
      if (clab::kcmp(clab::key(STRING_ELT(type, 0)), Key{kinds[i], 0}) == 0) q = i;
    }
  }
  if (q < 0) fail("`type` must be \"or\", \"rr\" or \"rd\".");
  if (!Rf_isNull(digits) && !Rf_isNull(tol)) fail("give `digits` or `tol`, not both.");
  int dg = NA_INTEGER;
  double tl = 0;
  if (!Rf_isNull(tol)) {
    SEXP t = PROTECT(num(tol, "tol", 1, 0.0));
    tl = REAL_RO(t)[0];
    UNPROTECT(1);
  } else {
    const int d = Rf_isNull(digits) ? ndec(val) : clab::whole(digits, "digits", 0);
    if (d > 15) fail("`digits` must be at most 15.");
    if (d >= 0) {
      dg = d;
      tl = 0.5 * std::pow(10.0, -d);
    }
  }
  SEXP rn = VECTOR_ELT(tb, 1), cn = VECTOR_ELT(tb, 2);
  const int gi = Rf_isNull(group) ? -1 : clab::ref(group, rn, r, "group", "row") - 1;
  const int hi = Rf_isNull(ref) ? -1 : clab::ref(ref, rn, r, "ref", "row") - 1;
  const int ei = Rf_isNull(event) ? -1 : clab::ref(event, cn, c, "event", "column") - 1;
  if (gi >= 0 && gi == hi) fail("`group` and `ref` must name two different rows.");
  const int m = 3 * r * c * (r + c - 2);
  int *grp = ibuf(m), *g = ibuf(m), *h = ibuf(m), *e = ibuf(m), *kd = ibuf(m), *mt = ibuf(m);
  double *v = dbuf(m);
  int nm = 0;
  kchk(clab_relchk(REAL_RO(VECTOR_ELT(tb, 0)), r, c, val, tl, grp, g, h, e, kd, v,
                   mt, &nm),
       nullptr);
  SEXP rd = PROTECT(Rf_allocVector(VECSXP, 8));
  for (int j = 0; j < 5; ++j) SET_VECTOR_ELT(rd, j, Rf_allocVector(STRSXP, m));
  SET_VECTOR_ELT(rd, 5, real(m));
  SET_VECTOR_ELT(rd, 6, Rf_allocVector(LGLSXP, m));
  SET_VECTOR_ELT(rd, 7, Rf_allocVector(LGLSXP, m));
  SEXP side = PROTECT(Rf_allocVector(STRSXP, 5));
  SET_STRING_ELT(side, 0, Rf_mkChar("rows"));
  SET_STRING_ELT(side, 1, Rf_mkChar("cols"));
  for (int j = 0; j < 3; ++j) SET_STRING_ELT(side, 2 + j, Rf_mkChar(kinds[j]));
  int hold = 0;
  for (int k = 0; k < m; ++k) {
    SEXP gl = grp[k] == 0 ? rn : cn, el = grp[k] == 0 ? cn : rn;
    SET_STRING_ELT(VECTOR_ELT(rd, 0), k, STRING_ELT(side, grp[k]));
    SET_STRING_ELT(VECTOR_ELT(rd, 1), k, STRING_ELT(side, 2 + kd[k]));
    SET_STRING_ELT(VECTOR_ELT(rd, 2), k, STRING_ELT(gl, g[k]));
    SET_STRING_ELT(VECTOR_ELT(rd, 3), k, STRING_ELT(gl, h[k]));
    SET_STRING_ELT(VECTOR_ELT(rd, 4), k, STRING_ELT(el, e[k]));
    REAL(VECTOR_ELT(rd, 5))[k] = v[k];
    LOGICAL(VECTOR_ELT(rd, 6))[k] = mt[k];
    const bool cl = grp[k] == 0 && kd[k] == q && (gi < 0 || g[k] == gi) &&
                    (hi < 0 || h[k] == hi) && (ei < 0 || e[k] == ei);
    LOGICAL(VECTOR_ELT(rd, 7))[k] = cl ? 1 : 0;
    if (cl && mt[k] != 0) ++hold;
  }
  static const char *const rdnm[] = {"grp", "type", "group", "ref", "event",
                                     "value", "match", "claim"};
  clab::frame(rd, rdnm);
  SEXP cv = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(cv, 0, gi < 0 ? NA_STRING : STRING_ELT(rn, gi));
  SET_STRING_ELT(cv, 1, hi < 0 ? NA_STRING : STRING_ELT(rn, hi));
  SET_STRING_ELT(cv, 2, ei < 0 ? NA_STRING : STRING_ELT(cn, ei));
  static const char *const cvnm[] = {"group", "ref", "event"};
  SEXP cvn = Rf_allocVector(STRSXP, 3);
  Rf_setAttrib(cv, R_NamesSymbol, cvn);
  for (int j = 0; j < 3; ++j) SET_STRING_ELT(cvn, j, Rf_mkChar(cvnm[j]));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 8));
  SET_VECTOR_ELT(out, 0, Rf_mkString(hold > 0 ? "holds" : (nm > 0 ? "other" : "none")));
  SET_VECTOR_ELT(out, 1, Rf_ScalarLogical(hold > 0 ? 1 : 0));
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(val));
  SET_VECTOR_ELT(out, 3, Rf_ScalarReal(tl));
  SET_VECTOR_ELT(out, 4, Rf_ScalarInteger(dg));
  SET_VECTOR_ELT(out, 5, Rf_mkString(kinds[q]));
  SET_VECTOR_ELT(out, 6, cv);
  SET_VECTOR_ELT(out, 7, rd);
  static const char *const nmo[] = {"verdict", "ok", "value", "tol", "digits",
                                    "type", "claim", "reads"};
  clab::nlist(out, nmo);
  Rf_classgets(out, PROTECT(Rf_mkString("relchk")));
  UNPROTECT(7);
  return out;
}

}  // extern "C"

namespace {

// The text of a label with white space (ASCII) trimmed at both ends and
// each run of it inside made one space, and ASCII letters in lower case:
// two labels with one fold differ only by case or white space.
SEXP lfold(SEXP s) {
  const char *p = Rf_translateCharUTF8(s);
  char *b = R_alloc(std::strlen(p) + 1, 1), *o = b;
  bool sp = false;
  for (; *p; ++p) {
    const unsigned char ch = static_cast<unsigned char>(*p);
    if (ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' || ch == '\f' || ch == '\v') {
      sp = o != b;
      continue;
    }
    if (sp) *o++ = ' ';
    sp = false;
    *o++ = ch >= 'A' && ch <= 'Z' ? static_cast<char>(ch - 'A' + 'a') : static_cast<char>(ch);
  }
  *o = '\0';
  return Rf_mkCharCE(b, CE_UTF8);
}

// The number of missing values of an atomic vector (NaN counts).
int nacnt(SEXP x) {
  const R_xlen_t n = XLENGTH(x);
  int k = 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    bool na = false;
    switch (TYPEOF(x)) {
      case REALSXP: na = std::isnan(REAL_RO(x)[i]); break;
      case INTSXP: na = INTEGER_RO(x)[i] == NA_INTEGER; break;
      case LGLSXP: na = LOGICAL_RO(x)[i] == NA_LOGICAL; break;
      case STRSXP: na = STRING_ELT(x, i) == NA_STRING; break;
      default: break;
    }
    if (na) ++k;
  }
  return k;
}

// What catchk() finds in one variable.
struct Var {
  const char *type;
  bool cat, lab, bin, codes;
  int n, k, dref;      // values, levels, the default reference (-1: none)
  const char *why;     // why that level, or why none
  SEXP lev, labl;      // levels and their value labels (NA where none)
  int *cnt, *twin, s[4], nolab, dlab;
};

// The levels of x in the order catchk() documents, their counts and
// twins from the kernel. kp, a list of 2 the caller keeps, receives the
// levels and their labels, so that they stay protected.
Var vaud(SEXP x, int nmin, SEXP kp) {
  Var v{};
  v.n = ilen(static_cast<double>(Rf_xlength(x)), "x");
  v.dref = -1;
  v.lab = !Rf_isNull(Rf_getAttrib(x, Rf_install("labels")));
  const int t = TYPEOF(x);
  v.cat = true;
  bool fac = false;
  if (!Rf_isVectorAtomic(x) || !Rf_isNull(Rf_getAttrib(x, R_DimSymbol))) {
    v.type = Rf_isVectorAtomic(x) ? "matrix" : Rf_type2char(static_cast<SEXPTYPE>(t));
    v.cat = false;
  } else if (v.lab) {
    v.type = "labelled";
  } else if (Rf_isFactor(x)) {
    v.type = Rf_inherits(x, "ordered") ? "ordered" : "factor";
    fac = true;
  } else if (Rf_inherits(x, "Date") || Rf_inherits(x, "POSIXt") ||
             Rf_inherits(x, "difftime")) {
    v.type = CHAR(STRING_ELT(Rf_getAttrib(x, R_ClassSymbol), 0));
    v.cat = false;
  } else if (t == LGLSXP || t == INTSXP || t == REALSXP || t == STRSXP) {
    static const char *const tn[] = {"logical", "integer", "double", "character"};
    v.type = tn[t == LGLSXP ? 0 : t == INTSXP ? 1 : t == REALSXP ? 2 : 3];
  } else {
    v.type = Rf_type2char(static_cast<SEXPTYPE>(t));
    v.cat = false;
  }
  if (v.cat && t == REALSXP && !v.lab) {
    // only whole numbers are categories
    for (int i = 0; i < v.n && v.cat; ++i) {
      const double d = REAL_RO(x)[i];
      v.cat = std::isnan(d) || (std::isfinite(d) && d == std::floor(d));
    }
  }
  if (!v.cat) {
    v.why = "not categorical";
    v.s[0] = Rf_isVectorAtomic(x) ? nacnt(x) : 0;
    SET_VECTOR_ELT(kp, 0, Rf_allocVector(STRSXP, 0));
    SET_VECTOR_ELT(kp, 1, Rf_allocVector(STRSXP, 0));
    v.lev = VECTOR_ELT(kp, 0);
    v.labl = VECTOR_ELT(kp, 1);
    return v;
  }
  int *code = ibuf(v.n);
  const bool num = (t == INTSXP || t == REALSXP) && !fac && !v.lab;
  bool b01 = num;
  int seen = 0;
  for (int i = 0; i < v.n && b01; ++i) {
    const double d = celd(x, i);
    if (std::isnan(d)) continue;
    ++seen;
    b01 = d == 0 || d == 1;
  }
  v.bin = (t == LGLSXP && !v.lab) || (b01 && seen > 0);
  SEXP lb = R_NilValue;
  if (fac) {
    SET_VECTOR_ELT(kp, 0, Rf_getAttrib(x, R_LevelsSymbol));
    v.k = static_cast<int>(Rf_xlength(VECTOR_ELT(kp, 0)));
    for (int i = 0; i < v.n; ++i) {
      const int c = INTEGER_RO(x)[i];
      code[i] = c == NA_INTEGER || c < 1 || c > v.k ? -1 : c - 1;
    }
  } else if (v.bin) {
    SEXP lv = Rf_allocVector(STRSXP, 2);
    SET_VECTOR_ELT(kp, 0, lv);
    SET_STRING_ELT(lv, 0, Rf_mkChar(t == LGLSXP ? "FALSE" : "0"));
    SET_STRING_ELT(lv, 1, Rf_mkChar(t == LGLSXP ? "TRUE" : "1"));
    v.k = 2;
    for (int i = 0; i < v.n; ++i) {
      double d;
      if (t == LGLSXP) {
        const int b = LOGICAL_RO(x)[i];
        d = b == NA_LOGICAL ? NA_REAL : static_cast<double>(b);
      } else {
        d = celd(x, i);
      }
      code[i] = std::isnan(d) ? -1 : static_cast<int>(d);
    }
  } else {
    // the distinct values as text, by first appearance; for a labelled
    // vector, then the labelled values that do not occur
    SEXP s = PROTECT(Rf_coerceVector(x, STRSXP));
    SEXP lv = Rf_getAttrib(x, Rf_install("labels"));
    SEXP ls = PROTECT(v.lab ? Rf_coerceVector(lv, STRSXP) : Rf_allocVector(STRSXP, 0));
    const R_xlen_t nl = XLENGTH(ls);
    SEXP all = PROTECT(Rf_allocVector(STRSXP, v.n + nl));
    // NaN is missing, as is.na() says, though as.character() gives "NaN"
    for (int i = 0; i < v.n; ++i) {
      const bool nan = t == REALSXP && std::isnan(REAL_RO(x)[i]);
      SET_STRING_ELT(all, i, nan ? NA_STRING : STRING_ELT(s, i));
    }
    for (R_xlen_t j = 0; j < nl; ++j) SET_STRING_ELT(all, v.n + j, STRING_ELT(ls, j));
    const Keys ks = keys(all);
    int *ac = ibuf(static_cast<double>(ks.n));
    int k = 0;
    for (R_xlen_t i = 0; i < ks.n; ++i) {
      if (ks.k[i].s == nullptr) {
        ac[i] = -1;
      } else if (ks.first[i]) {
        ac[i] = k++;
      } else {
        ac[i] = ac[ks.slot[probe(ks, ks.k[i])]];
      }
    }
    v.k = k;
    SEXP lev = Rf_allocVector(STRSXP, k);
    SET_VECTOR_ELT(kp, 0, lev);
    for (R_xlen_t i = 0; i < ks.n; ++i) {
      if (ac[i] >= 0 && ks.first[i]) SET_STRING_ELT(lev, ac[i], STRING_ELT(all, i));
    }
    for (int i = 0; i < v.n; ++i) code[i] = ac[i];
    SEXP ln = Rf_getAttrib(lv, R_NamesSymbol);
    if (v.lab && !Rf_isNull(ln)) {
      // each level's label: the first label of its value
      lb = Rf_allocVector(STRSXP, k);
      SET_VECTOR_ELT(kp, 1, lb);
      for (int j = 0; j < k; ++j) SET_STRING_ELT(lb, j, NA_STRING);
      for (R_xlen_t j = 0; j < nl; ++j) {
        const int a = ac[v.n + j];
        if (a >= 0 && STRING_ELT(lb, a) == NA_STRING) SET_STRING_ELT(lb, a, STRING_ELT(ln, j));
      }
      // labels whose text another value has too
      const Keys kl = keys(ln);
      int *nk = ibuf(static_cast<double>(kl.n));
      for (R_xlen_t j = 0; j < kl.n; ++j) nk[j] = 0;
      for (R_xlen_t j = 0; j < kl.n; ++j) ++nk[kl.slot[probe(kl, kl.k[j])]];
      for (R_xlen_t j = 0; j < kl.n; ++j) {
        if (kl.k[j].s != nullptr && nk[kl.slot[probe(kl, kl.k[j])]] > 1) ++v.dlab;
      }
    }
    UNPROTECT(3);
  }
  if (Rf_isNull(lb)) {
    lb = Rf_allocVector(STRSXP, v.k);
    SET_VECTOR_ELT(kp, 1, lb);
    for (int j = 0; j < v.k; ++j) SET_STRING_ELT(lb, j, NA_STRING);
  }
  v.lev = VECTOR_ELT(kp, 0);
  v.labl = lb;
  // fold keys, coded by first appearance
  SEXP fl = PROTECT(Rf_allocVector(STRSXP, v.k));
  for (int j = 0; j < v.k; ++j) SET_STRING_ELT(fl, j, lfold(STRING_ELT(v.lev, j)));
  int nd;
  const int *fk = fcode(fl, &nd, nullptr);
  UNPROTECT(1);
  v.cnt = ibuf(v.k);
  v.twin = ibuf(v.k);
  kchk(clab_catchk(code, v.n, v.k, fk, nmin, v.cnt, v.twin, v.s, ibuf(v.k)), nullptr);
  for (int j = 0; j < v.k; ++j) {
    if (v.lab && v.cnt[j] > 0 && STRING_ELT(v.labl, j) == NA_STRING) ++v.nolab;
  }
  v.codes = num && !v.bin && v.k >= 1 && v.k <= 10;
  if (t == LGLSXP && !v.lab) {
    v.dref = 0;
    v.why = "logical: TRUE is the treatment";
  } else if (v.bin) {
    v.dref = 0;
    v.why = "0 and 1: 1 is the treatment";
  } else if (fac) {
    v.dref = v.k > 0 ? 0 : -1;
    v.why = v.k > 0 ? "factor: the first level" : "factor with no levels";
  } else if (v.lab) {
    v.why = "labelled: decode it first";
  } else {
    v.why = "no order: state it";
  }
  return v;
}

// catchk()'s verdict on a stated reference w (text) for v.
const char *refst(const Var &v, SEXP w, int nmin) {
  for (int j = 0; j < v.k; ++j) {
    if (!clab::seql(STRING_ELT(v.lev, j), w)) continue;
    if (v.cnt[j] == 0) return "empty";
    if (v.cnt[j] < nmin) return "low";
    return j == v.dref ? "default" : "other";
  }
  SEXP fw = PROTECT(lfold(w));
  for (int j = 0; j < v.k; ++j) {
    SEXP fj = PROTECT(lfold(STRING_ELT(v.lev, j)));
    const bool eq = clab::seql(fj, fw);
    UNPROTECT(1);
    if (eq) {
      UNPROTECT(1);
      return "case";
    }
  }
  UNPROTECT(1);
  return "absent";
}

// one element of ref as text: an atomic value of length 1, a factor by
// its label; NA_STRING for anything else
SEXP ref1(SEXP e) {
  if (!Rf_isVectorAtomic(e) || XLENGTH(e) != 1) return NA_STRING;
  SEXP s = clab::labs(e);
  return Rf_isNull(s) ? NA_STRING : STRING_ELT(s, 0);
}

}  // namespace

extern "C" {

//' @srrstats {EA5.3} catchk() summarises each column: its storage type,
//'   its levels and their counts, missing values, and the problems found.
//' @srrstats {G2.4c} Levels are compared as as.character() gives them; a
//'   factor by its labels.
//' @noRd
SEXP C_clab_r_catchk(SEXP x, SEXP ref, SEXP nmin) {
  clab::who = "catchk";
  const int mn = clab::whole(nmin, "nmin", 1);
  SEXP vars, nms;
  if (TYPEOF(x) == VECSXP) {
    vars = PROTECT(x);
    SEXP n0 = Rf_getAttrib(x, R_NamesSymbol);
    nms = PROTECT(Rf_allocVector(STRSXP, XLENGTH(x)));
    char b[32];
    for (R_xlen_t i = 0; i < XLENGTH(x); ++i) {
      SEXP s = Rf_isNull(n0) ? NA_STRING : STRING_ELT(n0, i);
      if (s == NA_STRING || CHAR(s)[0] == '\0') {
        clab::lfmt(static_cast<double>(i + 1), b, sizeof b);
        s = Rf_mkChar(b);
      }
      SET_STRING_ELT(nms, i, s);
    }
  } else if (Rf_isVectorAtomic(x)) {
    clab::vec(x, "x");
    vars = PROTECT(Rf_allocVector(VECSXP, 1));
    SET_VECTOR_ELT(vars, 0, x);
    nms = PROTECT(Rf_mkString("x"));
  } else {
    fail("`x` must be a vector, a data frame or a list of vectors.");
  }
  const int nv = ilen(static_cast<double>(XLENGTH(vars)), "x");
  if (nv == 0) fail("`x` has no variables.");
  const Keys kn = keys(nms);
  for (int i = 0; i < nv; ++i) {
    if (!kn.first[i]) fail("`x` names a variable twice: %s.", Rf_translateChar(STRING_ELT(nms, i)));
  }
  // the stated reference of each variable (NA: none)
  SEXP want = PROTECT(Rf_allocVector(STRSXP, nv));
  for (int i = 0; i < nv; ++i) SET_STRING_ELT(want, i, NA_STRING);
  if (!Rf_isNull(ref)) {
    clab::unlab(ref, "ref");
    const bool l = TYPEOF(ref) == VECSXP;
    if (!l && !Rf_isVectorAtomic(ref)) fail("`ref` must be levels: a vector or a list.");
    SEXP rv = PROTECT(l ? ref : clab::labs(ref));
    if (Rf_isNull(rv)) fail("`ref` must be levels: a vector or a list.");
    const R_xlen_t nr = XLENGTH(ref);
    SEXP rn = Rf_getAttrib(ref, R_NamesSymbol);
    if (nr == 0) fail("`ref` must name at least one level.");
    for (R_xlen_t j = 0; j < nr; ++j) {
      SEXP w = l ? ref1(VECTOR_ELT(rv, j)) : STRING_ELT(rv, j);
      if (w == NA_STRING) fail("`ref` must be levels with no missing values, one per variable.");
      if (Rf_isNull(rn)) {
        if (nr != 1) fail("`ref` must be one level, or levels named by variable.");
        for (int i = 0; i < nv; ++i) SET_STRING_ELT(want, i, w);
        continue;
      }
      SEXP nj = STRING_ELT(rn, j);
      const R_xlen_t i = kn.slot[probe(kn, clab::key(nj))];
      if (i < 0) fail("`ref` names no variable: %s.", Rf_translateChar(nj));
      if (STRING_ELT(want, i) != NA_STRING) fail("`ref` names a variable twice: %s.", Rf_translateChar(nj));
      SET_STRING_ELT(want, i, w);
    }
    UNPROTECT(1);
  }
  // one row per variable; the levels of all, one row each
  static const char *const vnm[] = {"var", "type", "n", "nna", "nlev", "unused", "low",
                                    "twin", "bin", "codes", "nolab", "dlab", "ref",
                                    "why", "want", "state", "ok"};
  SEXP vd = PROTECT(Rf_allocVector(VECSXP, 17));
  static const int vty[] = {STRSXP, STRSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, LGLSXP,
                            LGLSXP, INTSXP, INTSXP, STRSXP, STRSXP, STRSXP, STRSXP, LGLSXP};
  for (int j = 0; j < 17; ++j) SET_VECTOR_ELT(vd, j, Rf_allocVector(static_cast<SEXPTYPE>(vty[j]), nv));
  SEXP keep = PROTECT(Rf_allocVector(VECSXP, nv));
  int **cnts = reinterpret_cast<int **>(R_alloc(static_cast<std::size_t>(nv), sizeof(int *)));
  int **tws = reinterpret_cast<int **>(R_alloc(static_cast<std::size_t>(nv), sizeof(int *)));
  int *drefs = ibuf(nv);
  double tot = 0;
  bool all = true;
  for (int i = 0; i < nv; ++i) {
    SEXP kp = Rf_allocVector(VECSXP, 2);
    SET_VECTOR_ELT(keep, i, kp);
    const Var v = vaud(VECTOR_ELT(vars, i), mn, kp);
    cnts[i] = v.cnt;
    tws[i] = v.twin;
    drefs[i] = v.dref;
    tot += v.k;
    SEXP w = STRING_ELT(want, i);
    const char *st = w == NA_STRING ? nullptr : (v.cat ? refst(v, w, mn) : "absent");
    const bool good = st == nullptr || std::strcmp(st, "default") == 0 || std::strcmp(st, "other") == 0;
    const bool ok = v.cat ? v.s[0] == 0 && v.s[1] == 0 && v.s[2] == 0 && v.s[3] == 0 &&
                                !v.lab && !v.codes && v.k - v.s[1] >= 2 && good
                          : v.s[0] == 0 && good;
    all = all && ok;
    const int na = NA_INTEGER;
    SET_STRING_ELT(VECTOR_ELT(vd, 0), i, STRING_ELT(nms, i));
    SET_STRING_ELT(VECTOR_ELT(vd, 1), i, Rf_mkChar(v.type));
    INTEGER(VECTOR_ELT(vd, 2))[i] = v.n;
    INTEGER(VECTOR_ELT(vd, 3))[i] = v.s[0];
    INTEGER(VECTOR_ELT(vd, 4))[i] = v.cat ? v.k : na;
    INTEGER(VECTOR_ELT(vd, 5))[i] = v.cat ? v.s[1] : na;
    INTEGER(VECTOR_ELT(vd, 6))[i] = v.cat ? v.s[2] : na;
    INTEGER(VECTOR_ELT(vd, 7))[i] = v.cat ? v.s[3] : na;
    LOGICAL(VECTOR_ELT(vd, 8))[i] = v.cat && v.bin ? 1 : 0;
    LOGICAL(VECTOR_ELT(vd, 9))[i] = v.cat && v.codes ? 1 : 0;
    INTEGER(VECTOR_ELT(vd, 10))[i] = v.nolab;
    INTEGER(VECTOR_ELT(vd, 11))[i] = v.dlab;
    SET_STRING_ELT(VECTOR_ELT(vd, 12), i, v.dref < 0 ? NA_STRING : STRING_ELT(v.lev, v.dref));
    SET_STRING_ELT(VECTOR_ELT(vd, 13), i, Rf_mkChar(v.why));
    SET_STRING_ELT(VECTOR_ELT(vd, 14), i, w);
    SET_STRING_ELT(VECTOR_ELT(vd, 15), i, st == nullptr ? NA_STRING : Rf_mkChar(st));
    LOGICAL(VECTOR_ELT(vd, 16))[i] = ok ? 1 : 0;
  }
  clab::frame(vd, vnm);
  const int nl = ilen(tot, "x");
  static const char *const lnm[] = {"var", "level", "label", "n", "low", "twin", "ref"};
  static const int lty[] = {STRSXP, STRSXP, STRSXP, INTSXP, LGLSXP, STRSXP, LGLSXP};
  SEXP ld = PROTECT(Rf_allocVector(VECSXP, 7));
  for (int j = 0; j < 7; ++j) SET_VECTOR_ELT(ld, j, Rf_allocVector(static_cast<SEXPTYPE>(lty[j]), nl));
  for (int i = 0, r = 0; i < nv; ++i) {
    SEXP lev = VECTOR_ELT(VECTOR_ELT(keep, i), 0), lab = VECTOR_ELT(VECTOR_ELT(keep, i), 1);
    for (int j = 0; j < XLENGTH(lev); ++j, ++r) {
      SET_STRING_ELT(VECTOR_ELT(ld, 0), r, STRING_ELT(nms, i));
      SET_STRING_ELT(VECTOR_ELT(ld, 1), r, STRING_ELT(lev, j));
      SET_STRING_ELT(VECTOR_ELT(ld, 2), r, STRING_ELT(lab, j));
      INTEGER(VECTOR_ELT(ld, 3))[r] = cnts[i][j];
      LOGICAL(VECTOR_ELT(ld, 4))[r] = cnts[i][j] > 0 && cnts[i][j] < mn ? 1 : 0;
      SET_STRING_ELT(VECTOR_ELT(ld, 5), r, tws[i][j] < 0 ? NA_STRING : STRING_ELT(lev, tws[i][j]));
      LOGICAL(VECTOR_ELT(ld, 6))[r] = j == drefs[i] ? 1 : 0;
    }
  }
  clab::frame(ld, lnm);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out, 0, Rf_ScalarLogical(all ? 1 : 0));
  SET_VECTOR_ELT(out, 1, Rf_ScalarInteger(mn));
  SET_VECTOR_ELT(out, 2, vd);
  SET_VECTOR_ELT(out, 3, ld);
  static const char *const onm[] = {"ok", "nmin", "vars", "levels"};
  clab::nlist(out, onm);
  Rf_classgets(out, PROTECT(Rf_mkString("catchk")));
  UNPROTECT(8);
  return out;
}

}  // extern "C"
