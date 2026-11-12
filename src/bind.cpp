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
#include <cstdint>
#include <cstring>
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

const char *const vodnm[] = {"group", "share_day", "share_night", "ratio",
                             "sign", "adj"};

}  // namespace

extern "C" {

// Disparity.

SEXP C_clab_r_vodtst(SEXP day, SEXP night, SEXP ref, SEXP shift) {
  SEXP nms = PROTECT(clab::grp(day));
  SEXP d = PROTECT(num(day, "day", -1, 0.0));
  const int n = ilen(static_cast<double>(XLENGTH(d)), "day");
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
  UNPROTECT(3);
  return clab::nvec(2, nm, r);
}

SEXP C_clab_r_dspbm(SEXP pop, SEXP contact, SEXP force, SEXP ref) {
  SEXP nms = PROTECT(clab::grp(pop));
  SEXP p = PROTECT(num(pop, "pop", -1, 0.0, true));
  const int n = ilen(static_cast<double>(XLENGTH(p)), "pop");
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
  kchk(clab_dyadrt(REAL_RO(x), off, vic, ilen(static_cast<double>(m), "count"), REAL_RO(p), g,
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

// Forensics.

namespace {

// all(is.na(x)) of a q that is not numeric (is.na() on each type as R's
// do_isna has it, the warning included)
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
    case VECSXP:
      for (R_xlen_t i = 0; i < n; ++i) {
        SEXP e = VECTOR_ELT(x, i);
        if (!(Rf_isVectorAtomic(e) && XLENGTH(e) == 1 && clab::anyna(e))) return false;
      }
      return true;
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
  SEXP nms = PROTECT(clab::grp(p));
  SEXP pv = PROTECT(num(p, "p", -1, 0.0, true));
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
    nq = clab::isnum(q);
    if ((!nq && !allna(q)) || Rf_xlength(q) != n) {
      fail("`q` must be numbers, one per locus.");
    }
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
  UNPROTECT(6);
  return clab::nlist(out, nm);
}

SEXP C_clab_r_cllr(SEXP same, SEXP diff) {
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
  UNPROTECT(3);
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

// One empirical cross-entropy:
// (pr * mean(.l1pinv(s * odds)) + (1 - pr) * mean(log1p(d * odds))) / log(2)
double ece(const double *s, R_xlen_t ns, const double *d, R_xlen_t nd,
           double odds, double pr, double *t) {
  // .l1pinv(x): log(1 + 1 / x) without forming 1 / x where it overflows:
  //   small <- x < 1; out <- log1p(1 / x)
  //   out[small] <- log1p(x[small]) - log(x[small])
  // A NaN x (0 * Inf) makes small NA, and R refuses that assignment when
  // more than one value is assigned.
  R_xlen_t nsm = 0, nna = 0;
  for (R_xlen_t k = 0; k < ns; ++k) {
    const double x = s[k] * odds;
    nna += std::isnan(x);
    nsm += x < 1;
    t[k] = x < 1 ? std::log1p(x) - std::log(x) : std::log1p(1 / x);
  }
  if (nna > 0 && nsm + nna > 1) {
    Rf_errorcall(PROTECT(clab::lang("out[small] <- log1p(x[small]) - log(x[small])")),  // NOLINT
                 "NAs are not allowed in subscripted assignments");
  }
  const double a = clab::nmean(t, ns);
  for (R_xlen_t k = 0; k < nd; ++k) t[k] = std::log1p(d[k] * odds);
  const double b = clab::nmean(t, nd);
  return (clab::rd(pr * a) + clab::rd((1 - pr) * b)) / std::log(2.0);
}

}  // namespace

extern "C" {

SEXP C_clab_r_cllrmin(SEXP same, SEXP diff) {
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
  Cal c = cal(same, diff);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, c.cs);
  SET_VECTOR_ELT(out, 1, c.cd);
  static const char *const nm[] = {"same", "diff"};
  UNPROTECT(5);
  return clab::nlist(out, nm);
}

SEXP C_clab_r_tipplt(SEXP same, SEXP diff) {
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
    const double pr = odds / (1 + odds);
    e[i] = ece(REAL_RO(c.s), ns, REAL_RO(c.d), nd, odds, pr, t);
    m[i] = ece(REAL_RO(c.cs), ns, REAL_RO(c.cd), nd, odds, pr, t);
    r[i] = ece(&one, 1, &one, 1, odds, pr, t);
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

// .gf(): the frequencies of the alleles of g, same shape, row by row
const double *gf(SEXP freq, SEXP g, const char *nm) {
  if (Rf_isNull(Rf_getAttrib(freq, R_DimSymbol)) && Rf_xlength(freq) == 2) {
    freq = mat12(freq);
  }
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

SEXP kinidx(SEXP g1, SEXP g2, SEXP freq, SEXP rel) {
  SEXP a = PROTECT(gt(g1, "g1", -1));
  const int nr = nrow(a);
  SEXP b = PROTECT(gt(g2, "g2", nr));
  const double *f = gf(freq, a, "g1");
  double k[3];
  kin(rel, k);
  SEXP g[2] = {a, b};
  int *c[2];
  codes(g, 2, nr, c, "g1");
  SEXP ix = PROTECT(real(nr));
  kchk(clab_kinidx(c[0], c[1], f, nr, k, REAL(ix)),
       ".idx(g1, .Call(C_clab_kinidx, .code(g1, u), .code(g2, u), f, k))");
  SEXP out = idx(a, ix);
  UNPROTECT(3);
  return out;
}

}  // namespace

extern "C" {

SEXP C_clab_r_kinidx(SEXP g1, SEXP g2, SEXP freq, SEXP rel) {
  return kinidx(g1, g2, freq, rel);
}

SEXP C_clab_r_patidx(SEXP child, SEXP father, SEXP freq, SEXP mother) {
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
  const double *f = gf(freq, c, "child");
  SEXP g[3] = {c, fa, mo};
  int *k[3];
  codes(g, 3, nr, k, "child");
  SEXP ix = PROTECT(real(nr));
  kchk(clab_patidx(k[0], k[2], k[1], f, nr, REAL(ix)),
       ".idx(child, .Call(C_clab_patidx, .code(child, u), .code(mother, u), "
       ".code(father, u), f))");
  SEXP out = idx(c, ix);
  UNPROTECT(4);
  return out;
}

}  // extern "C"

// Policing.

extern "C" {

SEXP C_clab_r_clrsmp(SEXP clra, SEXP reca, SEXP clrb, SEXP recb) {
  SEXP nms = PROTECT(clab::grp(clra));
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
  UNPROTECT(7);
  return clab::nlist(out, nm);
}

SEXP C_clab_r_expbnd(SEXP y, SEXP m, SEXP gamma, SEXP ref) {
  SEXP nms = PROTECT(clab::grp(y));
  SEXP yv = PROTECT(num(y, "y", -1, 0.0, true));
  const int n = ilen(static_cast<double>(XLENGTH(yv)), "y");
  SEXP mv = PROTECT(num(m, "m", n, 0.0, true));
  SEXP g = PROTECT(num(gamma, "gamma", 1, 1.0));
  const int i = clab::ref(ref, nms, n);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 8));
  SET_VECTOR_ELT(out, 0, nms);
  for (int k = 1; k < 7; ++k) SET_VECTOR_ELT(out, k, real(n));
  SET_VECTOR_ELT(out, 7, Rf_allocVector(INTSXP, n));
  kchk(clab_expbnd(REAL_RO(yv), REAL_RO(mv), n, i - 1, REAL_RO(g)[0],
                   REAL(VECTOR_ELT(out, 1)), REAL(VECTOR_ELT(out, 2)),
                   REAL(VECTOR_ELT(out, 3)), REAL(VECTOR_ELT(out, 4)),
                   REAL(VECTOR_ELT(out, 5)), REAL(VECTOR_ELT(out, 6)),
                   INTEGER(VECTOR_ELT(out, 7))),
       nullptr);
  static const char *const nm[] = {"group", "rate", "lo",  "hi",
                                   "ratio", "rlo",  "rhi", "sign"};
  clab::frame(out, nm);
  UNPROTECT(5);
  return out;
}

SEXP C_clab_r_misrd(SEXP r1, SEXP r0, SEXP prev, SEXP se, SEXP sp) {
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
  UNPROTECT(4);
  return clab::nvec(7, nm, r);
}

// Corrections.

SEXP C_clab_r_incap(SEXP lam, SEXP q, SEXP s, SEXP share) {
  SEXP nms = PROTECT(clab::grp(lam));
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
  SEXP g = Rf_allocVector(STRSXP, nr);
  SET_VECTOR_ELT(out, 0, g);
  for (int i = 0; i < n; ++i) SET_STRING_ELT(g, i, STRING_ELT(nms, i));
  if (has) SET_STRING_ELT(g, n, Rf_mkChar("all"));
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
  UNPROTECT(6);
  return out;
}

SEXP C_clab_r_incapr(SEXP lam, SEXP t0, SEXP s, SEXP r) {
  SEXP l = PROTECT(num(lam, "lam", -1, 0.0));
  const R_xlen_t n = XLENGTH(l);
  const double *lp = REAL_RO(l);
  // any(diff(lam) > 0)
  for (R_xlen_t i = 1; i < n; ++i) {
    if (lp[i] - lp[i - 1] > 0) fail("`lam` must never increase.");
  }
  const int sv = clab::whole(s, "s", 1);
  const int tv = clab::whole(t0, "t0", 0);
  if (tv > INT_MAX - sv) {
    // t0 + s overflows R's integers: NA with a warning, and if () stops
    // on the NA. Both conditions are rebuilt as R raised them.
    Rf_warningcall(PROTECT(clab::lang("t0 + s")),
                   "NAs produced by integer overflow");
    SEXP cl = PROTECT(clab::lang(
        "if (t0 + s > length(lam) - 1) {\n"
        "  stop(\"`lam` must have length at least t0 + s + 1.\", call. = FALSE)\n"
        "}"));
    Rf_errorcall(cl, "missing value where TRUE/FALSE needed");
  }
  // t0 + s > length(lam) - 1
  if (static_cast<double>(tv + sv) > static_cast<double>(n) - 1) {
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
  const R_xlen_t n = Rf_xlength(rel);
  clab::lgl(rel, "rel", n);
  clab::lgl(rule, "rule", n);
  if (!clab::isnum(y) || XLENGTH(y) != n) {
    fail("`y` must be numbers, one per defendant.");
  }
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
  UNPROTECT(2);
  return clab::nvec(6, nm, r);
}

SEXP C_clab_r_losbnd(SEXP done, SEXP open) {
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

SEXP C_clab_r_sntbnd(SEXP y, SEXP z, SEXP treat, SEXP w) {
  SEXP yv = PROTECT(num(y, "y", -1, 0.0));
  const int n = ilen(static_cast<double>(XLENGTH(yv)), "y");
  const double *yp = REAL_RO(yv);
  for (int i = 0; i < n; ++i) {
    if (!(yp[i] == 0 || yp[i] == 1)) fail("`y` must be 0 or 1.");
  }
  if (Rf_xlength(z) != n || clab::anyna(z)) {
    fail("`z` must have length %d and no missing values.", n);
  }
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
  // treat: one of lv, lv[2L] by default
  SEXP tr = lv[1];
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
  UNPROTECT(5);
  return out;
}

SEXP C_clab_r_hzdsel(SEXP s, SEXP h, SEXP l, SEXP surh, SEXP surl) {
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
  double r[5];
  kchk(clab_hzdsel(sv, hv, lv, REAL_RO(a), REAL_RO(b), r), nullptr);
  static const char *const nm[] = {"share0", "share1", "haz0", "haz1", "ratio"};
  UNPROTECT(2);
  return clab::nvec(5, nm, r);
}

// A hook for the tests, on lengths no test can allocate: with a name, the
// length check of the binding; with FALSE, that of the kernel entry
// points; with NULL, the length as an error message would print it.
SEXP C_clab_tlen(SEXP n, SEXP nm) {
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
    {"lrpost", "Forensics", "post_eq"},
    {"lrpost", "Forensics", "post_lt_half"},
    {"cllr", "Forensics", "cllr_neutral"},
    {"dbsrch", "Forensics", "np_rule"},
    {"mixlr", "Forensics", "mix_two_alleles"},
    {"kinidx", "Forensics", "po_identical_het"},
    {"patidx", "Forensics", "post_even_prior"},
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
    {"rtmhs", "RegMean", "selected_change_nonpos"},
    {"rtmhs", "RegMean", "low_selected_change_nonneg"},
    {"incap", "Incap", "steady_state_rate"},
    {"incap", "Incap", "prevented_share_eq"},
    {"incap", "Incap", "prevented_share_lt_one"},
    {"incap", "Incap", "rate_antitone_in_S"},
    {"incap", "Incap", "rate_antitone_in_q"},
    {"incap", "Incap", "marginal_prevention_eq"},
    {"incap", "Incap", "high_rate_more_prevented"},
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
};
constexpr int nthm = sizeof thms / sizeof thms[0];
}  // namespace

SEXP C_clab_r_proof(SEXP fn) {
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
