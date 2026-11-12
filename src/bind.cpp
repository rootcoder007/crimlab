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
  const int n = ilen(XLENGTH(pv), "p");
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
  kchk(clab_cllr(REAL_RO(s), ilen(XLENGTH(s), "same"), REAL_RO(d),
                 ilen(XLENGTH(d), "diff"), &r),
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
  Keys nd{nullptr, nullptr, 0};
  if (!bad) {
    nd = keys(need);
    for (R_xlen_t i = 0; i < nd.n && !bad; ++i) bad = !has(nk, nd.k[i]);
    bad = bad || ndist(nd) > 20;
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
  const int n = ilen(XLENGTH(pv), "p");
  // is.null(nms) || anyDuplicated(nms) || sum(p) > 1 + 1e-12
  bool bad = Rf_isNull(nms);
  Keys nk{nullptr, nullptr, 0};
  if (!bad) {
    nk = keys(nms);
    bad = ndist(nk) != nk.n;
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
  const int ns = ilen(XLENGTH(c.s), "same"), nd = ilen(XLENGTH(c.d), "diff");
  if (ns > INT_MAX - nd) fail("`same` and `diff` are too long for crimlab: more than %d values.", INT_MAX);  // NOLINT
  const int n = ns + nd;
  const double *sp = REAL_RO(c.s), *dp = REAL_RO(c.d);
  // lr <- c(same, diff); o <- order(lr): a stable order, as R's radix sort
  int *o = reinterpret_cast<int *>(R_alloc(static_cast<std::size_t>(n), sizeof(int)));
  for (int j = 0; j < n; ++j) o[j] = j;
  auto v = [sp, dp, ns](int j) { return j < ns ? sp[j] : dp[j - ns]; };
  std::stable_sort(o, o + n, [&v](int a, int b) { return v(a) < v(b); });
  double *lr = reinterpret_cast<double *>(R_alloc(static_cast<std::size_t>(n), sizeof(double)));
  int *is = reinterpret_cast<int *>(R_alloc(static_cast<std::size_t>(n), sizeof(int)));
  for (int k = 0; k < n; ++k) {
    lr[k] = v(o[k]);
    is[k] = o[k] < ns;
  }
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
  for (R_xlen_t i = 0; i < n; ++i) lp[i] = std::log10(xp[i]);
  std::sort(lp, lp + n);
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
void codes(SEXP *g, int m, int nr, int **out) {
  if (nr > INT_MAX / (2 * m)) fail("`%s` is too long for crimlab: more than %d values.", "g1", INT_MAX);  // NOLINT
  const R_xlen_t w = 2 * static_cast<R_xlen_t>(nr), all = m * w;
  Key *k = reinterpret_cast<Key *>(R_alloc(static_cast<std::size_t>(all), sizeof(Key)));
  for (int j = 0; j < m; ++j) {
    SEXP s = PROTECT(clab::labs(g[j]));
    for (R_xlen_t i = 0; i < w; ++i) k[j * w + i] = clab::key(STRING_ELT(s, i));
    UNPROTECT(1);
  }
  R_xlen_t *o = reinterpret_cast<R_xlen_t *>(R_alloc(static_cast<std::size_t>(all), sizeof(R_xlen_t)));
  for (R_xlen_t i = 0; i < all; ++i) o[i] = i;
  std::sort(o, o + all, [k](R_xlen_t a, R_xlen_t b) { return clab::kcmp(k[a], k[b]) < 0; });
  int *cls = reinterpret_cast<int *>(R_alloc(static_cast<std::size_t>(all), sizeof(int)));
  int id = 0;
  for (R_xlen_t i = 0; i < all; ++i) {
    if (i > 0 && clab::kcmp(k[o[i - 1]], k[o[i]]) != 0) ++id;
    cls[o[i]] = id + 1;
  }
  for (int j = 0; j < m; ++j) {
    out[j] = reinterpret_cast<int *>(R_alloc(static_cast<std::size_t>(w), sizeof(int)));
    for (int i = 0; i < nr; ++i) {
      out[j][2 * i] = cls[j * w + i];
      out[j][2 * i + 1] = cls[j * w + nr + i];
    }
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
  SEXP nm = PROTECT(Rf_allocVector(STRSXP, nr));
  if (nr == 1 && Rf_isNull(rn) != Rf_isNull(cn)) {
    SET_STRING_ELT(nm, 0, STRING_ELT(Rf_isNull(rn) ? cn : rn, 0));
  } else if (nr > 1 && !Rf_isNull(rn)) {
    for (int i = 0; i < nr; ++i) SET_STRING_ELT(nm, i, STRING_ELT(rn, i));
  } else {
    char b[16];
    for (int i = 0; i < nr; ++i) {
      std::snprintf(b, sizeof b, "%d", i + 1);
      SET_STRING_ELT(nm, i, Rf_mkChar(b));
    }
  }
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
  codes(g, 2, nr, c);
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
  codes(g, 3, nr, k);
  SEXP ix = PROTECT(real(nr));
  kchk(clab_patidx(k[0], k[2], k[1], f, nr, REAL(ix)),
       ".idx(child, .Call(C_clab_patidx, .code(child, u), .code(mother, u), "
       ".code(father, u), f))");
  SEXP out = idx(c, ix);
  UNPROTECT(4);
  return out;
}

}  // extern "C"
