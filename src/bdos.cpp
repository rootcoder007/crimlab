// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of crimlab's dossier functions (fairness, trace
// evidence, aoristic weights, Benford screens, Scott-Knott grouping,
// constant reporters, sensitivity bounds, trajectory models, glass
// transfer, lowess, the Chumbley score, grid Gi* and the point pattern
// test): each checks its arguments, calls one kernel and builds the R
// result. Nothing calls back into the R evaluator.
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstring>
#include <limits>
#include <vector>

#define R_NO_REMAP_RMATH
#include <Rmath.h>

#include "crimlab.h"
#include "rchk.h"
#include "spin.h"

using clab::fail;
using clab::num;

namespace {

SEXP real(R_xlen_t n) { return Rf_allocVector(REALSXP, n); }
SEXP ints(R_xlen_t n) { return Rf_allocVector(INTSXP, n); }
int len(SEXP x, const char *nm) { return clab::ilen(static_cast<double>(XLENGTH(x)), nm); }

void kc(int rc, const char *msg) {
  if (rc != CLAB_OK) fail("%s", msg);
}

// a data frame of k double columns of length n
SEXP dframe(int k, int n, const char *const *nm, double **col) {
  SEXP v = PROTECT(Rf_allocVector(VECSXP, k));
  for (int j = 0; j < k; ++j) {
    SET_VECTOR_ELT(v, j, real(n));
    col[j] = REAL(VECTOR_ELT(v, j));
  }
  SEXP f = clab::frame(v, nm);
  UNPROTECT(1);
  return f;
}

// prepend a character column (labels) to the data frame f, named first
SEXP withlab(SEXP f, SEXP lab, const char *first) {
  const R_xlen_t k = XLENGTH(f);
  SEXP v = PROTECT(Rf_allocVector(VECSXP, k + 1));
  SET_VECTOR_ELT(v, 0, lab);
  SEXP on = Rf_getAttrib(f, R_NamesSymbol);
  const char *nm[32];
  nm[0] = first;
  for (R_xlen_t j = 0; j < k; ++j) {
    SET_VECTOR_ELT(v, j + 1, VECTOR_ELT(f, j));
    nm[j + 1] = CHAR(STRING_ELT(on, j));
  }
  SEXP out = clab::frame(v, nm);
  UNPROTECT(1);
  return out;
}

// Labels to codes from 0: a factor by its levels that occur, in level
// order; anything else by its sorted distinct values (as.character()).
// lev receives the labels of the codes. No missing labels.
struct lab {
  SEXP code, lev;
  int g;
};

lab codes(SEXP x, const char *nm, R_xlen_t n) {
  if (Rf_xlength(x) != n) fail("`%s` must have one label per value, %ld.", nm, static_cast<long>(n));
  SEXP s = PROTECT(clab::labs(x));
  if (Rf_isNull(s) || clab::anyna(s)) {
    fail("`%s` must be labels (numbers, text or a factor) with no missing "
         "values.", nm);
  }
  SEXP o = PROTECT(ints(n));
  int *ord = INTEGER(o);
  for (R_xlen_t i = 0; i < n; ++i) ord[i] = static_cast<int>(i);
  std::vector<clab::Key> k(static_cast<std::size_t>(n));
  for (R_xlen_t i = 0; i < n; ++i) k[i] = clab::key(STRING_ELT(s, i));
  if (Rf_isFactor(x)) {
    const int *f = INTEGER_RO(x);
    std::stable_sort(ord, ord + n, [f](int a, int b) { return f[a] < f[b]; });
  } else {
    std::stable_sort(ord, ord + n, [&k](int a, int b) { return clab::kcmp(k[a], k[b]) < 0; });
  }
  lab r{PROTECT(ints(n)), R_NilValue, 0};
  int g = -1;
  for (R_xlen_t q = 0; q < n; ++q) {
    const int i = ord[q];
    if (q == 0 || clab::kcmp(k[ord[q - 1]], k[i]) != 0) ++g;
    INTEGER(r.code)[i] = g;
  }
  r.g = g + 1;
  r.lev = PROTECT(Rf_allocVector(STRSXP, r.g));
  for (R_xlen_t i = 0; i < n; ++i) SET_STRING_ELT(r.lev, INTEGER(r.code)[i], STRING_ELT(s, i));
  UNPROTECT(4);
  PROTECT(r.code);
  PROTECT(r.lev);
  return r;
}

// 0 or 1 outcomes: logical or numbers
SEXP bin01(SEXP y, const char *nm, R_xlen_t n) {
  SEXP v = PROTECT(ints(n));
  if (Rf_xlength(y) != n || (TYPEOF(y) != LGLSXP && !clab::isnum(y))) {
    fail("`%s` must be 0 or 1 (or TRUE or FALSE), one per case.", nm);
  }
  for (R_xlen_t i = 0; i < n; ++i) {
    const double d = TYPEOF(y) == LGLSXP
                         ? (LOGICAL_RO(y)[i] == NA_LOGICAL ? NA_REAL : LOGICAL_RO(y)[i])
                         : (TYPEOF(y) == INTSXP ? (INTEGER_RO(y)[i] == NA_INTEGER ? NA_REAL : INTEGER_RO(y)[i])
                                                : REAL_RO(y)[i]);
    if (!(d == 0 || d == 1)) fail("`%s` must be 0 or 1 (or TRUE or FALSE), one per case.", nm);
    INTEGER(v)[i] = static_cast<int>(d);
  }
  UNPROTECT(1);
  return v;
}

// a numeric array of nd dimensions (a matrix for nd = 2): its dimensions
// into dim and its values as doubles, finite (or NaN when na)
SEXP narr(SEXP x, const char *nm, int nd, int *dim, bool na = false) {
  SEXP d = Rf_getAttrib(x, R_DimSymbol);
  if (!clab::isnum(x) || TYPEOF(d) != INTSXP || XLENGTH(d) != nd) {
    fail(nd == 2 ? "`%s` must be a numeric matrix." : "`%s` must be a numeric array of 3 dimensions.", nm);
  }
  for (int j = 0; j < nd; ++j) dim[j] = INTEGER_RO(d)[j];
  SEXP v = PROTECT(clab::dcopy(PROTECT(clab::dbl(x))));
  for (R_xlen_t i = 0; i < XLENGTH(v); ++i) {
    const double e = REAL(v)[i];
    if (na ? std::isinf(e) : !std::isfinite(e)) {
      fail(na ? "`%s` must hold finite numbers or NA." : "`%s` must hold finite numbers.", nm);
    }
  }
  UNPROTECT(2);
  return v;
}

SEXP nmat(SEXP x, const char *nm, int *nr, int *nc, bool na = false) {
  int d[2];
  SEXP v = narr(x, nm, 2, d, na);
  *nr = d[0];
  *nc = d[1];
  return v;
}

// a seed for the kernels' generators (spin.h)
unsigned int sd(SEXP s) { return csp::seed(s); }

const double qnan = std::numeric_limits<double>::quiet_NaN();

}  // namespace

extern "C" {

// ---------------------------------------------------------------- fairness

SEXP C_clab_r_fair(SEXP score, SEXP y, SEXP grp, SEXP cut, SEXP base,
                   SEXP conf) {
  clab::who = "fair";
  SEXP uv[3] = {score, y, grp};
  static const char *const un[] = {"score", "y", "grp"};
  const int nb = csp::byname(uv, un, 3);
  score = uv[0];
  y = uv[1];
  grp = uv[2];
  SEXP s = PROTECT(num(score, "score"));
  const int n = len(s, "score");
  SEXP yv = PROTECT(bin01(y, "y", n));
  lab g = codes(grp, "grp", n);
  const double ct = clab::val1(cut);
  if (std::isnan(ct)) fail("`cut` must be one number.");
  const int b = clab::ref(base, g.lev, g.g) - 1;
  const double z = Rf_qnorm5(0.5 + 0.5 * clab::one(conf, "conf", 0.0, 1.0), 0.0, 1.0, 1, 0);
  const int G = g.g;
  SEXP w = PROTECT(real(5 * static_cast<R_xlen_t>(G) + 20 * G + 16 * G));
  double *cnt = REAL(w), *met = cnt + 5 * G, *rat = met + 10 * G, *lo = rat + 10 * G,
         *hi = lo + 8 * G;
  kc(clab_fair(REAL_RO(s), INTEGER(yv), INTEGER(g.code), n, G, ct, b, z, cnt, met,
               rat, lo, hi),
     "the kernel refused the input.");
  static const char *const cn[] = {"tp", "fp", "tn", "fn", "n"};
  static const char *const mn[] = {"acc", "tpr", "tnr", "fpr", "fnr", "ppv",
                                   "npv", "ppos", "npos", "mcc"};
  const double *src[5] = {cnt, met, rat, lo, hi};
  const int wd[5] = {5, 10, 10, 8, 8};
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
  for (int t = 0; t < 5; ++t) {
    double *c[10];
    SEXP f = PROTECT(dframe(wd[t], G, t == 0 ? cn : mn, c));
    for (int k = 0; k < wd[t]; ++k) {
      for (int j = 0; j < G; ++j) c[k][j] = src[t][k * G + j];
    }
    SET_VECTOR_ELT(out, t, withlab(f, g.lev, "group"));
    UNPROTECT(1);
  }
  static const char *const nm[] = {"cnt", "met", "rat", "lo", "hi"};
  UNPROTECT(6 + nb);
  return clab::nlist(out, nm);
}

// ----------------------------------------------------------- trace evidence

SEXP C_clab_r_tlvc(SEXP x, SEXP item) {
  clab::who = "tlvc";
  int n, p;
  SEXP xv = PROTECT(nmat(x, "x", &n, &p));
  lab g = codes(item, "item", n);
  if (g.g < 2 || n <= g.g) {
    fail("at least 2 items are needed, and more observations than items.");
  }
  SEXP mu = PROTECT(real(p));
  SEXP u = PROTECT(Rf_allocMatrix(REALSXP, p, p));
  SEXP c = PROTECT(Rf_allocMatrix(REALSXP, p, p));
  SEXP wk = PROTECT(real(static_cast<R_xlen_t>(g.g) * p + 2 * g.g));
  kc(clab_tlvc(REAL(xv), n, p, INTEGER(g.code), g.g, REAL(mu), REAL(u), REAL(c),
               REAL(wk)),
     "the kernel refused the input.");
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, mu);
  SET_VECTOR_ELT(out, 1, u);
  SET_VECTOR_ELT(out, 2, c);
  static const char *const nm[] = {"mu", "u", "c"};
  UNPROTECT(8);
  return clab::nlist(out, nm);
}

SEXP C_clab_r_tlvlr(SEXP xa, SEXP xb, SEXP na, SEXP nb, SEXP bg) {
  clab::who = "tlvlr";
  if (TYPEOF(bg) != VECSXP || XLENGTH(bg) != 3) {
    fail("`bg` must be the list tlvc() returns.");
  }
  SEXP mu = PROTECT(num(VECTOR_ELT(bg, 0), "bg$mu"));
  const int p = len(mu, "bg$mu");
  int r1, c1, r2, c2;
  SEXP u = PROTECT(nmat(VECTOR_ELT(bg, 1), "bg$u", &r1, &c1));
  SEXP c = PROTECT(nmat(VECTOR_ELT(bg, 2), "bg$c", &r2, &c2));
  if (r1 != p || c1 != p || r2 != p || c2 != p) {
    fail("`bg$u` and `bg$c` must be square with the length of `bg$mu`.");
  }
  SEXP a = PROTECT(num(xa, "xa", p));
  SEXP b = PROTECT(num(xb, "xb", p));
  const double ma = clab::one(na, "na", 0.0, std::numeric_limits<double>::infinity());
  const double mb = clab::one(nb, "nb", 0.0, std::numeric_limits<double>::infinity());
  SEXP wk = PROTECT(real(static_cast<R_xlen_t>(p) * p + 2 * p));
  double o[2];
  kc(clab_tlvlr(REAL(u), REAL(c), REAL_RO(mu), p, REAL_RO(a), ma, REAL_RO(b), mb,
                o, REAL(wk)),
     "a covariance matrix of the model is not positive definite.");
  static const char *const nm[] = {"ln", "log10"};
  UNPROTECT(6);
  return clab::nvec(2, nm, o);
}

}  // extern "C"

namespace {

// instants: numbers or date-times (POSIXct), as seconds; NA allowed, not
// infinite
SEXP inst(SEXP x, const char *nm) {
  const bool tm = Rf_inherits(x, "POSIXct");
  if (!(tm || clab::isnum(x))) fail("`%s` must be numbers of seconds or date-times (POSIXct).", nm);
  SEXP v = PROTECT(clab::dcopy(PROTECT(clab::dbl(x))));
  for (R_xlen_t i = 0; i < XLENGTH(v); ++i) {
    if (std::isinf(REAL(v)[i])) fail("`%s` must not be infinite.", nm);
  }
  UNPROTECT(2);
  return v;
}

}  // namespace

extern "C" {

// ------------------------------------------------------------ aoristic

SEXP C_clab_r_aorwt(SEXP from, SEXP to, SEXP bw, SEXP nb, SEXP origin,
                    SEXP off, SEXP tt, SEXP toff, SEXP lo, SEXP hi,
                    SEXP each) {
  clab::who = "aorwt";
  SEXP f = PROTECT(inst(from, "from"));
  const int n = len(f, "from");
  SEXP t = PROTECT(inst(to, "to"));
  if (XLENGTH(t) != n) fail("`to` must have the length of `from`.");
  const double w = clab::one(bw, "bw", 0.0, std::numeric_limits<double>::infinity());
  const int k = clab::whole(nb, "nb", 1);
  SEXP org = PROTECT(inst(origin, "origin"));
  if (XLENGTH(org) != 1 || std::isnan(REAL(org)[0])) fail("`origin` must be one instant.");
  SEXP o0 = PROTECT(num(off, "off", 1));
  if (Rf_isNull(tt) != Rf_isNull(toff)) fail("give both `tt` and `toff`, or neither.");
  SEXP ttv = PROTECT(Rf_isNull(tt) ? real(0) : inst(tt, "tt"));
  SEXP tov = PROTECT(Rf_isNull(toff) ? real(0) : num(toff, "toff", XLENGTH(ttv)));
  const int nt = len(ttv, "tt");
  const double l = clab::val1(lo), h = clab::val1(hi);
  if (!(l >= 0) || !(h > l)) fail("`lo` and `hi` must be numbers with 0 <= lo < hi.");
  const bool ea = csp::flag(each, "each");
  SEXP wm = PROTECT(ea ? Rf_allocMatrix(REALSXP, n, k) : R_NilValue);
  SEXP tot = PROTECT(real(k));
  int cnt[5];
  kc(clab_aorwt(REAL(f), REAL(t), n, REAL(ttv), REAL(tov), nt, REAL_RO(o0)[0],
                REAL(org)[0], w, k, l, h, ea ? REAL(wm) : nullptr, REAL(tot), cnt),
     "`tt` must be strictly increasing, with finite offsets.");
  double c[5];
  for (int j = 0; j < 5; ++j) c[j] = cnt[j];
  static const char *const cn[] = {"left", "noend", "reversed", "short", "long"};
  SEXP cv = PROTECT(clab::nvec(5, cn, c));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, tot);
  SET_VECTOR_ELT(out, 1, wm);
  SET_VECTOR_ELT(out, 2, cv);
  static const char *const nm[] = {"tot", "w", "cnt"};
  UNPROTECT(10);
  return clab::nlist(out, nm);
}

// ------------------------------------------------------------- Benford

SEXP C_clab_r_bdig(SEXP x, SEXP nd, SEXP dp, SEXP sgn) {
  clab::who = "bdig";
  if (!clab::isnum(x)) fail("`x` must be numbers.");
  SEXP v = PROTECT(clab::dbl(x));
  const int n = len(v, "x");
  const int d = clab::whole(nd, "nd", 1);
  if (d > 9) fail("`nd` must be from 1 to 9.");
  const double dv = clab::val1(dp);
  if (!(dv == std::floor(dv) && dv >= -1 && dv <= 22)) fail("`dp` must be a whole number from -1 to 22.");
  static const char *const so[] = {"pos", "neg", "both"};
  const int s = csp::pick(sgn, "sgn", so, 3);
  SEXP dig = PROTECT(ints(n));
  int cnt[4];
  kc(clab_bdig(REAL_RO(v), n, d, static_cast<int>(dv), s == 0 ? 1 : (s == 1 ? -1 : 0),
               INTEGER(dig), cnt),
     "the kernel refused the input.");
  double c[4];
  for (int j = 0; j < 4; ++j) c[j] = cnt[j];
  static const char *const cn[] = {"zero", "nonfinite", "sign", "range"};
  SEXP cv = PROTECT(clab::nvec(4, cn, c));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, dig);
  SET_VECTOR_ELT(out, 1, cv);
  static const char *const nm[] = {"dig", "cnt"};
  UNPROTECT(4);
  return clab::nlist(out, nm);
}

SEXP C_clab_r_bnft(SEXP dig, SEXP nd, SEXP nsim, SEXP seed) {
  clab::who = "bnft";
  SEXP v = PROTECT(num(dig, "dig", -1, 0.0));
  const int n = len(v, "dig");
  const int d = clab::whole(nd, "nd", 1);
  if (d > 6) fail("`nd` must be from 1 to 6.");
  const int lo = static_cast<int>(std::pow(10.0, d - 1));
  SEXP iv = PROTECT(ints(n));
  for (int i = 0; i < n; ++i) {
    const double e = REAL(v)[i];
    if (e != std::floor(e) || (e != 0 && (e < lo || e >= 10.0 * lo))) {
      fail("`dig` must hold %d-digit numbers (or 0, skipped), as bdig() gives.", d);
    }
    INTEGER(iv)[i] = static_cast<int>(e);
  }
  const int ns = clab::whole(nsim, "nsim", 0);
  const unsigned int s = ns > 0 ? sd(seed) : 0u;
  const int k = 9 * lo;
  SEXP obs = PROTECT(real(k));
  SEXP wk = PROTECT(real(k));
  double o[6];
  kc(clab_bnft(INTEGER(iv), n, d, ns, s, REAL(obs), REAL(wk), o), "the kernel refused the input.");
  static const char *const tn[] = {"digit", "count", "exp"};
  double *c[3];
  SEXP tab = PROTECT(dframe(3, k, tn, c));
  for (int j = 0; j < k; ++j) {
    c[0][j] = lo + j;
    c[1][j] = REAL(obs)[j];
    c[2][j] = o[0] * std::log10(1.0 + 1.0 / (lo + j));
  }
  static const char *const sn[] = {"m", "chisq", "df", "p", "pmc", "mad"};
  SEXP st = PROTECT(clab::nvec(6, sn, o));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, tab);
  SET_VECTOR_ELT(out, 1, st);
  static const char *const nm[] = {"obs", "stat"};
  UNPROTECT(7);
  return clab::nlist(out, nm);
}

// ---------------------------------------------------- Scott-Knott, panels

SEXP C_clab_r_skp(SEXP x, SEXP s2, SEXP alpha, SEXP nsim, SEXP seed) {
  clab::who = "skp";
  SEXP v = PROTECT(num(x, "x"));
  const int n = len(v, "x");
  const double s = clab::one(s2, "s2", 0.0, std::numeric_limits<double>::infinity());
  const double a = clab::one(alpha, "alpha", 0.0, 1.0);
  const int ns = clab::whole(nsim, "nsim", 1);
  const unsigned int sv = sd(seed);
  SEXP grp = PROTECT(ints(n));
  SEXP iw = PROTECT(ints(3 * static_cast<R_xlen_t>(n)));
  SEXP wk = PROTECT(real(2 * static_cast<R_xlen_t>(n)));
  int ng = 0;
  kc(clab_skp(REAL_RO(v), n, s, a, ns, sv, INTEGER(grp), &ng, INTEGER(iw), REAL(wk)),
     "the kernel refused the input.");
  UNPROTECT(4);
  return grp;
}

SEXP C_clab_r_rptc(SEXP unit, SEXP year, SEXP rep, SEXP y0, SEXP y1,
                   SEXP minrep) {
  clab::who = "rptc";
  SEXP uv[3] = {unit, year, rep};
  static const char *const un[] = {"unit", "year", "rep"};
  const int nb = csp::byname(uv, un, 3);
  unit = uv[0];
  year = uv[1];
  rep = uv[2];
  SEXP yv = PROTECT(num(year, "year"));
  const int n = len(yv, "year");
  lab u = codes(unit, "unit", n);
  SEXP rv = PROTECT(num(rep, "rep", n));
  const int a0 = clab::whole(y0, "y0", 0);
  const int a1 = clab::whole(y1, "y1", a0);
  const int mr = clab::whole(minrep, "minrep", 0);
  if (static_cast<double>(a1 - a0 + 1) * u.g > 1e9) fail("too many units times years.");
  SEXP yi = PROTECT(ints(n));
  SEXP ri = PROTECT(ints(n));
  for (int i = 0; i < n; ++i) {
    const double e = REAL(yv)[i], r = REAL(rv)[i];
    if (e != std::floor(e) || std::fabs(e) > 1e9 || r != std::floor(r) || std::fabs(r) > 1e9) {
      fail("`year` and `rep` must be whole numbers.");
    }
    INTEGER(yi)[i] = static_cast<int>(e);
    INTEGER(ri)[i] = static_cast<int>(r);
  }
  SEXP fl = PROTECT(ints(u.g));
  SEXP wk = PROTECT(ints(static_cast<R_xlen_t>(u.g) * (a1 - a0 + 1)));
  int nc = 0;
  kc(clab_rptc(INTEGER(u.code), INTEGER(yi), INTEGER(ri), n, u.g, a0, a1, mr,
               INTEGER(fl), &nc, INTEGER(wk)),
     "the kernel refused the input.");
  SEXP lg = PROTECT(Rf_allocVector(LGLSXP, u.g));
  for (int j = 0; j < u.g; ++j) LOGICAL(lg)[j] = INTEGER(fl)[j];
  SEXP v = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(v, 0, u.lev);
  SET_VECTOR_ELT(v, 1, lg);
  static const char *const nm[] = {"unit", "const"};
  SEXP out = clab::frame(v, nm);
  UNPROTECT(10 + nb);
  return out;
}

// ------------------------------------------------------ sensitivity bounds

SEXP C_clab_r_rswp(SEXP d, SEXP gam, SEXP alpha) {
  clab::who = "rswp";
  SEXP v = PROTECT(num(d, "d"));
  const int n = len(v, "d");
  SEXP g = PROTECT(num(gam, "gam", -1, 1.0));
  const int ng = len(g, "gam");
  const double z = Rf_qnorm5(clab::one(alpha, "alpha", 0.0, 0.5), 0.0, 1.0, 0, 0);
  double gc = 0, ts = 0;
  static const char *const nm[] = {"gam", "plo", "pup"};
  double *c[3];
  SEXP f = PROTECT(dframe(3, ng, nm, c));
  std::memcpy(c[0], REAL_RO(g), ng * sizeof(double));
  SEXP wk = PROTECT(real(n));
  SEXP iw = PROTECT(ints(n));
  kc(clab_rswp(REAL_RO(v), n, REAL_RO(g), ng, z, c[1], c[2], &gc, &ts, REAL(wk),
               INTEGER(iw)),
     "the kernel refused the input.");
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, f);
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(gc));
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(ts));
  static const char *const on[] = {"p", "gcrit", "t"};
  UNPROTECT(6);
  return clab::nlist(out, on);
}

SEXP C_clab_r_rshl(SEXP d, SEXP gam) {
  clab::who = "rshl";
  SEXP v = PROTECT(num(d, "d"));
  const int n = len(v, "d");
  SEXP g = PROTECT(num(gam, "gam", -1, 1.0));
  const int ng = len(g, "gam");
  static const char *const nm[] = {"gam", "lo", "hi"};
  double *c[3];
  SEXP f = PROTECT(dframe(3, ng, nm, c));
  std::memcpy(c[0], REAL_RO(g), ng * sizeof(double));
  SEXP wk = PROTECT(real(n));
  kc(clab_rshl(REAL_RO(v), n, REAL_RO(g), ng, c[1], c[2], REAL(wk)),
     "`d` must hold at most 2^26 values.");
  UNPROTECT(4);
  return f;
}

}  // extern "C"

extern "C" {

// ------------------------------------------------- lowess, tool marks

SEXP C_clab_r_lwss(SEXP x, SEXP y, SEXP f, SEXP iter, SEXP delta) {
  clab::who = "lwss";
  SEXP uv[2] = {x, y};
  static const char *const un[] = {"x", "y"};
  const int nb = csp::byname(uv, un, 2);
  x = uv[0];
  y = uv[1];
  SEXP xv = PROTECT(num(x, "x"));
  const int n = len(xv, "x");
  SEXP yv = PROTECT(num(y, "y", n));
  const double sp = clab::one(f, "f", 0.0, std::numeric_limits<double>::infinity());
  const int it = clab::whole(iter, "iter", 0);
  SEXP o = PROTECT(ints(n));
  int *ord = INTEGER(o);
  const double *px = REAL_RO(xv);
  for (int i = 0; i < n; ++i) ord[i] = i;
  std::stable_sort(ord, ord + n, [px](int a, int b) { return px[a] < px[b]; });
  SEXP sx = PROTECT(real(n));
  SEXP sy = PROTECT(real(n));
  for (int i = 0; i < n; ++i) {
    REAL(sx)[i] = px[ord[i]];
    REAL(sy)[i] = REAL_RO(yv)[ord[i]];
  }
  const double dl = Rf_isNull(delta) ? 0.01 * (REAL(sx)[n - 1] - REAL(sx)[0])
                                     : clab::val1(delta);
  if (!(dl >= 0) || !std::isfinite(dl)) fail("`delta` must be one number from 0.");
  SEXP fit = PROTECT(real(n));
  SEXP wk = PROTECT(real(3 * static_cast<R_xlen_t>(n)));
  kc(clab_lwss(REAL(sx), REAL(sy), n, sp, it, dl, REAL(fit), REAL(wk)),
     "the kernel refused the input.");
  SEXP out = PROTECT(real(n));
  for (int i = 0; i < n; ++i) REAL(out)[ord[i]] = REAL(fit)[i];
  UNPROTECT(8 + nb);
  return out;
}

SEXP C_clab_r_chmb(SEXP a, SEXP b, SEXP wopt, SEXP wval, SEXP coarse) {
  clab::who = "chmb";
  SEXP av = PROTECT(num(a, "a"));
  SEXP bv = PROTECT(num(b, "b"));
  const int na = len(av, "a"), nb = len(bv, "b");
  const int wo = clab::whole(wopt, "wopt", 2);
  const int wv = clab::whole(wval, "wval", 2);
  if (wv > wo) fail("`wval` must be at most `wopt`.");
  const double co = clab::one(coarse, "coarse", 0.0, std::numeric_limits<double>::infinity());
  SEXP wk = PROTECT(real(3 * static_cast<R_xlen_t>(na) + 3 * static_cast<R_xlen_t>(nb) +
                         7 * static_cast<R_xlen_t>(std::max(na, nb))));
  int ns = 0, nd = 0;
  double u = 0, p = 0;
  const int rc = clab_chmb(REAL_RO(av), na, REAL_RO(bv), nb, wo, wv, co, &ns, &nd, &u,
                           &p, REAL(wk));
  kc(rc, rc == CLAB_EARG ? "a trimmed profile is shorter than `wopt`."
                         : "every window is flat: the profiles carry no marks.");
  double o[4] = {static_cast<double>(ns), static_cast<double>(nd), u, p};
  static const char *const nm[] = {"nsame", "ndiff", "u", "p"};
  UNPROTECT(3);
  return clab::nvec(4, nm, o);
}

// ------------------------------------------ grid Gi*, point pattern test

SEXP C_clab_r_gstr(SEXP cnt, SEXP ex, SEXP rad, SEXP nsim, SEXP seed) {
  clab::who = "gstr";
  int nr, nc;
  SEXP cv = PROTECT(nmat(cnt, "cnt", &nr, &nc));
  SEXP ev;
  if (Rf_isNull(ex)) {
    ev = PROTECT(real(XLENGTH(cv)));
    for (R_xlen_t i = 0; i < XLENGTH(cv); ++i) REAL(ev)[i] = 1.0;
  } else {
    int er, ec;
    ev = PROTECT(nmat(ex, "ex", &er, &ec));
    if (er != nr || ec != nc) fail("`ex` must have the dimensions of `cnt`.");
  }
  if (nr > 46340 || nc > 46340) fail("`cnt` must have at most 46340 rows and columns.");
  const double rd = clab::val1(rad);
  if (!(rd >= 0 && rd <= 1000)) fail("`rad` must be one number from 0 to 1000.");
  const int ns = clab::whole(nsim, "nsim", 0);
  const unsigned int s = ns > 0 ? sd(seed) : 0u;
  const R_xlen_t m = static_cast<R_xlen_t>(nr) * nc;
  const R_xlen_t q = 2 * static_cast<R_xlen_t>(std::floor(rd)) + 1;
  SEXP o = PROTECT(real(7 * m));
  SEXP wk = PROTECT(real(m));
  SEXP iw = PROTECT(ints(2 * q * q + 3 * m));
  kc(clab_gstr(REAL(cv), REAL(ev), nr, nc, rd, ns, s, REAL(o), REAL(wk), INTEGER(iw)),
     "`cnt` must be whole numbers from 0, `ex` numbers from 0, with at least 2 "
     "cells of positive exposure.");
  static const char *const nm[] = {"row", "col", "gz", "pz", "qz", "pp", "ze", "pe", "qe"};
  double *c[9];
  SEXP f = PROTECT(dframe(9, static_cast<int>(m), nm, c));
  for (R_xlen_t i = 0; i < m; ++i) {
    c[0][i] = static_cast<double>(i % nr + 1);
    c[1][i] = static_cast<double>(i / nr + 1);
    for (int j = 0; j < 7; ++j) c[2 + j][i] = REAL(o)[j * m + i];
  }
  UNPROTECT(6);
  return f;
}

SEXP C_clab_r_sptst(SEXP a, SEXP b, SEXP nboot, SEXP seed, SEXP alpha,
                    SEXP scale) {
  clab::who = "sptst";
  a = PROTECT(csp::gvec(a, "a"));
  b = PROTECT(csp::gvec(b, "b"));
  SEXP uv[2] = {a, b};
  static const char *const un[] = {"a", "b"};
  const int nb = csp::byname(uv, un, 2);
  a = uv[0];
  b = uv[1];
  if (!clab::isnum(a) || !clab::isnum(b)) fail("`a` and `b` must be counts.");
  SEXP av = PROTECT(clab::dcopy(PROTECT(clab::dbl(a))));
  SEXP bv = PROTECT(clab::dcopy(PROTECT(clab::dbl(b))));
  const int n = len(av, "a");
  if (n < 1 || XLENGTH(bv) != n) fail("`a` and `b` must have the same length, from 1.");
  const int nbt = clab::whole(nboot, "nboot", 1);
  const unsigned int s = sd(seed);
  const double al = clab::one(alpha, "alpha", 0.0, 1.0);
  const double sc = clab::one(scale, "scale", 0.0, std::numeric_limits<double>::infinity());
  SEXP o = PROTECT(real(3 * static_cast<R_xlen_t>(n)));
  SEXP wk = PROTECT(real(6 * static_cast<R_xlen_t>(n)));
  SEXP iw = PROTECT(ints(4 * static_cast<R_xlen_t>(n)));
  double g[8];
  kc(clab_sptst(REAL(av), REAL(bv), n, nbt, s, al, sc, REAL(o), g, REAL(wk), INTEGER(iw)),
     "`a` and `b` must be whole numbers from 0 (NA leaves a unit out).");
  static const char *const cn3[] = {"diff", "p", "flag"};
  double *c[3];
  SEXP f = PROTECT(dframe(3, n, cn3, c));
  for (int j = 0; j < 3; ++j) std::memcpy(c[j], REAL(o) + static_cast<R_xlen_t>(j) * n, n * sizeof(double));
  static const char *const gn[] = {"left", "empty", "elig", "flagged", "s", "z", "p", "nboot"};
  SEXP gv = PROTECT(clab::nvec(8, gn, g));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, f);
  SET_VECTOR_ELT(out, 1, gv);
  static const char *const nm[] = {"units", "glob"};
  UNPROTECT(10 + nb + 2);
  return clab::nlist(out, nm);
}

// ------------------------------------------------------- glass transfer

SEXP C_clab_r_tfrp(SEXP par, SEXP deffect, SEXP tdist, SEXP ymax, SEXP nsim,
                   SEXP seed, SEXP hl, SEXP nq) {
  clab::who = "tfrp";
  SEXP pv = PROTECT(num(par, "par", 17));
  const int de = csp::flag(deffect, "deffect") ? 1 : 0;
  static const char *const to[] = {"nb", "nbtrunc", "unif"};
  const int td = csp::pick(tdist, "tdist", to, 3);
  const int ym = clab::whole(ymax, "ymax", 0);
  const int ns = clab::whole(nsim, "nsim", 1);
  const unsigned int s = sd(seed);
  const double h = clab::val1(hl);
  const int q = clab::whole(nq, "nq", 1);
  int wl = 0;
  kc(clab_tfrw(REAL_RO(pv), de, td, ym, h, q, &wl),
     "`par` is outside the model's domain, or `hl` is outside 1e-4 to 1, or "
     "`nq` above 100000.");
  SEXP wk = PROTECT(real(wl));
  static const char *const nm[] = {"y", "p", "se"};
  double *c[3];
  SEXP f = PROTECT(dframe(3, ym + 1, nm, c));
  double cut[2];
  kc(clab_tfrp(REAL_RO(pv), de, td, ym, ns, s, h, q, c[1], c[2], cut, REAL(wk)),
     "the kernel refused the input.");
  for (int j = 0; j <= ym; ++j) c[0][j] = j;
  static const char *const cn[] = {"time", "grid"};
  SEXP cv = PROTECT(clab::nvec(2, cn, cut));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, f);
  SET_VECTOR_ELT(out, 1, cv);
  static const char *const on[] = {"prob", "cut"};
  UNPROTECT(5);
  return clab::nlist(out, on);
}

}  // extern "C"

namespace {

// group shares: ng positive numbers summing to 1
SEXP shares(SEXP prob, int *ng) {
  SEXP v = PROTECT(clab::dcopy(PROTECT(num(prob, "prob", -1, 0.0, true))));
  *ng = len(v, "prob");
  if (std::fabs(clab::nsum(REAL_RO(v), *ng) - 1) > 1e-9) fail("`prob` must sum to 1.");
  UNPROTECT(2);
  return v;
}

// the list of a fitted model: par (or beta and sig), prob, post, ll, iter
SEXP fitlist(int k, SEXP *el, const char *const *nm) {
  SEXP out = PROTECT(Rf_allocVector(VECSXP, k));
  for (int j = 0; j < k; ++j) SET_VECTOR_ELT(out, j, el[j]);
  UNPROTECT(1);
  return clab::nlist(out, nm);
}

// the ZIP trajectory model's inputs
struct zin {
  SEXP y, x, z, off, par, prob;
  int n, nt, p, q, ng;
};

zin zread(SEXP y, SEXP x, SEXP z, SEXP par, SEXP prob, SEXP off) {
  zin a{};
  a.y = PROTECT(nmat(y, "y", &a.n, &a.nt, true));
  int r, c;
  a.x = PROTECT(nmat(x, "x", &r, &a.p));
  if (r != a.nt) fail("`x` must have one row per period (column of `y`).");
  a.q = 0;
  a.z = PROTECT(Rf_isNull(z) ? real(0) : nmat(z, "z", &r, &a.q));
  if (!Rf_isNull(z) && r != a.nt) fail("`z` must have one row per period (column of `y`).");
  a.off = PROTECT(Rf_isNull(off) ? R_NilValue : nmat(off, "off", &r, &c));
  if (!Rf_isNull(off) && (r != a.n || c != a.nt)) fail("`off` must have the dimensions of `y`.");
  a.prob = PROTECT(shares(prob, &a.ng));
  a.par = PROTECT(nmat(par, "par", &r, &c));
  if (r != a.p + a.q || c != a.ng) {
    fail("`par` must have one row per column of `x` and `z` and one column per group.");
  }
  return a;
}

}  // namespace

extern "C" {

// ----------------------------------------------------- trajectory models

SEXP C_clab_r_zgll(SEXP y, SEXP x, SEXP z, SEXP par, SEXP prob, SEXP off) {
  clab::who = "zgll";
  zin a = zread(y, x, z, par, prob, off);
  SEXP post = PROTECT(Rf_allocMatrix(REALSXP, a.n, a.ng));
  SEXP wk = PROTECT(real(2 * static_cast<R_xlen_t>(a.nt)));
  double ll = 0;
  kc(clab_zgll(REAL(a.y), a.n, a.nt, Rf_isNull(a.off) ? nullptr : REAL(a.off), REAL(a.x),
               a.p, REAL(a.z), a.q, a.ng, REAL(a.par), REAL(a.prob), REAL(post), &ll,
               REAL(wk)),
     "`y` must hold whole counts from 0 (NA or a negative value is missing).");
  SEXP el[2] = {PROTECT(Rf_ScalarReal(ll)), post};
  static const char *const nm[] = {"ll", "post"};
  SEXP out = fitlist(2, el, nm);
  UNPROTECT(9);
  return out;
}

SEXP C_clab_r_zgem(SEXP y, SEXP x, SEXP z, SEXP par, SEXP prob, SEXP off,
                   SEXP maxit, SEXP tol) {
  clab::who = "zgem";
  zin a = zread(y, x, z, par, prob, off);
  const int mi = clab::whole(maxit, "maxit", 1);
  const double tl = clab::one(tol, "tol", 0.0, std::numeric_limits<double>::infinity());
  const int pq = a.p + a.q;
  SEXP post = PROTECT(Rf_allocMatrix(REALSXP, a.n, a.ng));
  SEXP wk = PROTECT(real(2 * static_cast<R_xlen_t>(pq) * pq + 3 * pq + 7 * static_cast<R_xlen_t>(a.nt)));
  SEXP pm = PROTECT(Rf_allocMatrix(REALSXP, pq, a.ng));
  std::memcpy(REAL(pm), REAL(a.par), static_cast<std::size_t>(pq) * a.ng * sizeof(double));
  double ll = 0;
  int it = 0;
  kc(clab_zgem(REAL(a.y), a.n, a.nt, Rf_isNull(a.off) ? nullptr : REAL(a.off), REAL(a.x),
               a.p, REAL(a.z), a.q, a.ng, REAL(pm), REAL(a.prob), REAL(post), &ll, mi, tl,
               &it, REAL(wk)),
     "`y` must hold whole counts from 0 (NA or a negative value is missing).");
  SEXP el[5] = {pm, a.prob, post, PROTECT(Rf_ScalarReal(ll)), PROTECT(Rf_ScalarInteger(it))};
  static const char *const nm[] = {"par", "prob", "post", "ll", "iter"};
  SEXP out = fitlist(5, el, nm);
  UNPROTECT(11);
  return out;
}

SEXP C_clab_r_mvtll(SEXP y, SEXP x, SEXP beta, SEXP sig, SEXP prob) {
  clab::who = "mvtll";
  int dy[3], dx[2], db[3], ds[3], ng;
  SEXP yv = PROTECT(narr(y, "y", 3, dy, true));
  SEXP xv = PROTECT(narr(x, "x", 2, dx));
  SEXP pv = PROTECT(shares(prob, &ng));
  SEXP bv = PROTECT(narr(beta, "beta", 3, db));
  SEXP sv = PROTECT(narr(sig, "sig", 3, ds));
  const int n = dy[0], nt = dy[1], r = dy[2], p = dx[1];
  if (dx[0] != nt || db[0] != p || db[1] != r || db[2] != ng || ds[0] != r || ds[1] != r ||
      ds[2] != ng) {
    fail("the dimensions must agree: `y` n x nt x r, `x` nt x p, `beta` p x r x "
         "groups, `sig` r x r x groups, `prob` one per group.");
  }
  SEXP post = PROTECT(Rf_allocMatrix(REALSXP, n, ng));
  SEXP wk = PROTECT(real(static_cast<R_xlen_t>(r) * r + r));
  double ll = 0;
  kc(clab_mvtll(REAL(yv), n, nt, r, REAL(xv), p, ng, REAL(bv), REAL(sv), REAL(pv),
                REAL(post), &ll, REAL(wk)),
     "a covariance in `sig` is not positive definite.");
  SEXP el[2] = {PROTECT(Rf_ScalarReal(ll)), post};
  static const char *const nm[] = {"ll", "post"};
  SEXP out = fitlist(2, el, nm);
  UNPROTECT(8);
  return out;
}

SEXP C_clab_r_mvtem(SEXP y, SEXP x, SEXP beta, SEXP sig, SEXP prob, SEXP maxit,
                    SEXP tol) {
  clab::who = "mvtem";
  int dy[3], dx[2], db[3], ds[3], ng;
  SEXP yv = PROTECT(narr(y, "y", 3, dy, true));
  SEXP xv = PROTECT(narr(x, "x", 2, dx));
  SEXP pv = PROTECT(shares(prob, &ng));
  SEXP bv = PROTECT(narr(beta, "beta", 3, db));
  SEXP sv = PROTECT(narr(sig, "sig", 3, ds));
  const int n = dy[0], nt = dy[1], r = dy[2], p = dx[1];
  if (dx[0] != nt || db[0] != p || db[1] != r || db[2] != ng || ds[0] != r || ds[1] != r ||
      ds[2] != ng) {
    fail("the dimensions must agree: `y` n x nt x r, `x` nt x p, `beta` p x r x "
         "groups, `sig` r x r x groups, `prob` one per group.");
  }
  const int mi = clab::whole(maxit, "maxit", 1);
  const double tl = clab::one(tol, "tol", 0.0, std::numeric_limits<double>::infinity());
  Rf_setAttrib(bv, R_DimSymbol, Rf_getAttrib(beta, R_DimSymbol));
  Rf_setAttrib(sv, R_DimSymbol, Rf_getAttrib(sig, R_DimSymbol));
  SEXP post = PROTECT(Rf_allocMatrix(REALSXP, n, ng));
  SEXP wk = PROTECT(real(static_cast<R_xlen_t>(p) * p + static_cast<R_xlen_t>(p) * r +
                         static_cast<R_xlen_t>(r) * r + r));
  double ll = 0;
  int it = 0;
  kc(clab_mvtem(REAL(yv), n, nt, r, REAL(xv), p, ng, REAL(bv), REAL(sv), REAL(pv),
                REAL(post), &ll, mi, tl, &it, REAL(wk)),
     "a weighted design is singular, a residual variance is 0, or a covariance "
     "is not positive definite.");
  SEXP el[6] = {bv, sv, pv, post, PROTECT(Rf_ScalarReal(ll)), PROTECT(Rf_ScalarInteger(it))};
  static const char *const nm[] = {"beta", "sig", "prob", "post", "ll", "iter"};
  SEXP out = fitlist(6, el, nm);
  UNPROTECT(9);
  return out;
}

}  // extern "C"
