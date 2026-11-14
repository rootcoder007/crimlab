// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the specialist methods, R/spec.R (fairness and the
// audits of risk tools and predictive policing, trace evidence, aoristic weights, Benford screens, Scott-Knott grouping,
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

// Older R headers declare these under their short names only when
// R_NO_REMAP_RMATH is set; the library has them as Rf_*, and newer headers
// declare the same names identically (rchk.h does this for Rf_qnorm5).
extern "C" double Rf_pt(double, double, int, int);
extern "C" double Rf_qt(double, double, int, int);
extern "C" double Rf_pf(double, double, double, int, int);

using clab::fail;
using clab::num;
using csp::bin01;
using csp::codes;
using csp::lab;

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

// -------------------------------------------- helpers of the fairness audits

// A data frame column by name, R_NilValue when absent.
SEXP col(SEXP df, const char *nm) {
  SEXP n = Rf_getAttrib(df, R_NamesSymbol);
  for (R_xlen_t j = 0; j < Rf_xlength(n); ++j) {
    if (std::strcmp(CHAR(STRING_ELT(n, j)), nm) == 0) return VECTOR_ELT(df, j);
  }
  return R_NilValue;
}

// Group labels for codes(): a vector (no matrix) of labels, not haven's
// labelled codes.
lab glab(SEXP x, const char *nm, R_xlen_t n) {
  clab::unlab(x, nm);
  clab::vec(x, nm);
  return codes(x, nm, n);
}

// 0 or 1 decisions or outcomes, as bin01(), from a vector that is not a
// labelled one.
SEXP b01(SEXP x, const char *nm, R_xlen_t n) {
  clab::unlab(x, nm);
  clab::vec(x, nm);
  return bin01(x, nm, n);
}

// The normal quantile of a two-sided interval at level conf.
double zconf(SEXP conf) {
  return Rf_qnorm5(0.5 + 0.5 * clab::one(conf, "conf", 0.0, 1.0), 0.0, 1.0, 1, 0);
}

// Columns of a result data frame of n rows, built in order. The list is
// protected by the constructor (one PROTECT the caller counts); done()
// returns the frame unprotected.
struct cols {
  SEXP v;
  const char *nm[40];
  int k = 0, n;
  cols(int ncol, int nr) : n(nr) { v = PROTECT(Rf_allocVector(VECSXP, ncol)); }
  double *dbl(const char *name) {
    SEXP c = Rf_allocVector(REALSXP, n);
    SET_VECTOR_ELT(v, k, c);
    nm[k++] = name;
    return REAL(c);
  }
  int *in(const char *name, int type = INTSXP) {
    SEXP c = Rf_allocVector(type, n);
    SET_VECTOR_ELT(v, k, c);
    nm[k++] = name;
    return type == LGLSXP ? LOGICAL(c) : INTEGER(c);
  }
  void put(const char *name, SEXP c) {
    SET_VECTOR_ELT(v, k, c);
    nm[k++] = name;
  }
  // the frame, of the k columns added
  SEXP done() {
    SEXP w = PROTECT(Rf_allocVector(VECSXP, k));
    for (int j = 0; j < k; ++j) SET_VECTOR_ELT(w, j, VECTOR_ELT(v, j));
    SEXP f = clab::frame(w, nm);
    UNPROTECT(1);
    return f;
  }
};

// A named double vector of k values (NaN kept).
SEXP named(int k, const char *const *nm, const double *x) { return clab::nvec(k, nm, x); }

// Planar points: a two-column numeric matrix or a data frame or list with
// columns x and y (or two unnamed columns), finite; into px, py (protected
// by the caller via the returned list). Returns their number; a NULL or
// zero-row input gives 0 when empty is allowed.
// a noise box: one number a for [0, a], or two, lo <= hi, all from 0
void box(SEXP x, const char *nm, double *lo, double *hi) {
  SEXP v = PROTECT(num(x, nm, -1, 0.0));
  const R_xlen_t k = XLENGTH(v);
  if (k > 2) fail("`%s` must be one number (the box [0, %s]) or two (its ends).", nm, nm);
  *lo = k == 1 ? 0.0 : REAL(v)[0];
  *hi = REAL(v)[k - 1];
  if (*lo > *hi) fail("`%s` must give its lower end first.", nm);
  UNPROTECT(1);
}

// names(x), or "1", "2", ... (protected by the caller)
SEXP nmsof(SEXP x, R_xlen_t n) {
  SEXP nm = Rf_getAttrib(x, R_NamesSymbol);
  return Rf_isNull(nm) ? clab::seqs(n) : nm;
}


int pts2(SEXP x, const char *nm, SEXP *px, SEXP *py, bool empty) {
  if (Rf_isNull(x) && empty) {
    *px = *py = R_NilValue;
    return 0;
  }
  SEXP cx = R_NilValue, cy = R_NilValue;
  SEXP d = Rf_getAttrib(x, R_DimSymbol);
  if (clab::isnum(x) && TYPEOF(d) == INTSXP && XLENGTH(d) == 2 && INTEGER_RO(d)[1] == 2) {
    const int nr = INTEGER_RO(d)[0];
    SEXP v = PROTECT(clab::dbl(x));
    cx = PROTECT(Rf_allocVector(REALSXP, nr));
    cy = PROTECT(Rf_allocVector(REALSXP, nr));
    for (int i = 0; i < nr; ++i) {
      REAL(cx)[i] = REAL_RO(v)[i];
      REAL(cy)[i] = REAL_RO(v)[i + nr];
    }
    UNPROTECT(3);
  } else if (TYPEOF(x) == VECSXP) {
    // columns x and y by name, or the two columns of a two-column table
    cx = col(x, "x");
    cy = col(x, "y");
    if ((Rf_isNull(cx) || Rf_isNull(cy)) && XLENGTH(x) == 2) {
      cx = VECTOR_ELT(x, 0);
      cy = VECTOR_ELT(x, 1);
    }
    if (Rf_isNull(cx) || Rf_isNull(cy)) fail("`%s` must have columns x and y.", nm);
  } else {
    fail("`%s` must be a two-column numeric matrix or a data frame with "
         "columns x and y.", nm);
  }
  if (!clab::isnum(cx) || !clab::isnum(cy) || XLENGTH(cx) != XLENGTH(cy)) {
    fail("`%s` must hold numeric coordinates, x and y of one length.", nm);
  }
  const int n = len(cx, nm);
  if (n == 0) {
    if (!empty) fail("`%s` must hold at least one point.", nm);
    *px = *py = R_NilValue;
    return 0;
  }
  *px = PROTECT(clab::dcopy(PROTECT(clab::dbl(cx))));
  *py = PROTECT(clab::dcopy(PROTECT(clab::dbl(cy))));
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(REAL(*px)[i]) || !std::isfinite(REAL(*py)[i])) {
      fail("`%s` must hold finite coordinates.", nm);
    }
  }
  UNPROTECT(4);
  return n;
}

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

SEXP C_clab_r_fairdf(SEXP pred, SEXP grp, SEXP y, SEXP base, SEXP conf,
                     SEXP nboot, SEXP seed) {
  clab::who = "fairdf";
  SEXP uv[3] = {pred, grp, y};
  static const char *const un[] = {"pred", "grp", "y"};
  const int nb = csp::byname(uv, un, 3);
  pred = uv[0];
  grp = uv[1];
  y = uv[2];
  const int n = len(pred, "pred");
  if (n < 1) fail("`pred` must hold at least one case.");
  SEXP pv = PROTECT(b01(pred, "pred", n));
  lab g = glab(grp, "grp", n);
  if (g.g < 2) fail("`grp` must hold at least two groups.");
  const bool hy = !Rf_isNull(y);
  SEXP yv = PROTECT(hy ? b01(y, "y", n) : R_NilValue);
  const int b = Rf_isNull(base) ? -1 : clab::ref(base, g.lev, g.g, "base") - 1;
  const double z = zconf(conf);
  const int B = clab::whole(nboot, "nboot", 0);
  const unsigned int s = B > 0 ? sd(seed) : 0U;
  const int G = g.g;
  SEXP w = PROTECT(real(40 * static_cast<R_xlen_t>(G) + 11 + 2 * static_cast<R_xlen_t>(B)));
  double *cnt = REAL(w), *met = cnt + 6 * G, *sum = met + 23 * G, *wk = sum + 11;
  kc(clab_fairdf(INTEGER(pv), hy ? INTEGER(yv) : nullptr, INTEGER(g.code), n, G, b, z, B, s,
                 cnt, met, sum, wk),
     "the kernel refused the input.");
  static const char *const mn[] = {"n", "npos", "ppos", "lo", "hi", "dp", "dplo", "dphi",
                                   "di", "dilo", "dihi", "tpr", "fpr", "dtpr", "dtprlo",
                                   "dtprhi", "dfpr", "dfprlo", "dfprhi", "aod", "aodlo",
                                   "aodhi", "eo"};
  const int k = hy ? 23 : 11;
  double *c[23] = {};
  SEXP f = PROTECT(dframe(k, G, mn, c));
  for (int t = 0; t < k; ++t) {
    for (int j = 0; j < G; ++j) c[t][j] = met[t * G + j];
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, withlab(f, g.lev, "group"));
  static const char *const sn[] = {"di", "dp", "gini", "bas", "ginilo", "ginihi",
                                   "baslo", "bashi", "aod", "eo"};
  SET_VECTOR_ELT(out, 1, named(hy ? 10 : 8, sn, sum));
  SET_VECTOR_ELT(out, 2, Rf_ScalarString(STRING_ELT(g.lev, static_cast<R_xlen_t>(sum[10]))));
  static const char *const nm[] = {"tab", "sum", "base"};
  UNPROTECT(7 + nb);
  return clab::nlist(out, nm);
}

SEXP C_clab_r_tmpaud(SEXP pred, SEXP grp, SEXP period, SEXP city, SEXP base) {
  clab::who = "tmpaud";
  SEXP uv[4] = {pred, grp, period, city};
  static const char *const un[] = {"pred", "grp", "period", "city"};
  const int nb = csp::byname(uv, un, 4);
  pred = uv[0];
  grp = uv[1];
  period = uv[2];
  city = uv[3];
  const int n = len(pred, "pred");
  if (n < 1) fail("`pred` must hold at least one case.");
  SEXP pv = PROTECT(b01(pred, "pred", n));
  lab g = glab(grp, "grp", n);
  if (g.g < 2) fail("`grp` must hold at least two groups.");
  lab per = glab(period, "period", n);
  SEXP cc, cl;
  int ncity = 1;
  if (Rf_isNull(city)) {
    cc = PROTECT(ints(n));
    for (int i = 0; i < n; ++i) INTEGER(cc)[i] = 0;
    cl = PROTECT(Rf_mkString("all"));
  } else {
    lab ct = glab(city, "city", n);
    cc = ct.code;
    cl = ct.lev;
    ncity = ct.g;
  }
  int np = 7;  // pv, g, per, the city codes and labels
  // the cells: the (city, period) pairs that occur, city then period order
  std::vector<long long> key(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    key[static_cast<std::size_t>(i)] =
        static_cast<long long>(INTEGER(cc)[i]) * per.g + INTEGER(per.code)[i];
  }
  std::vector<long long> u(key);
  std::sort(u.begin(), u.end());
  u.erase(std::unique(u.begin(), u.end()), u.end());
  const int nc = static_cast<int>(u.size());
  SEXP cell = PROTECT(ints(n)), cy = PROTECT(ints(nc));
  np += 2;
  for (int i = 0; i < n; ++i) {
    INTEGER(cell)[i] = static_cast<int>(
        std::lower_bound(u.begin(), u.end(), key[static_cast<std::size_t>(i)]) - u.begin());
  }
  for (int c = 0; c < nc; ++c) INTEGER(cy)[c] = static_cast<int>(u[static_cast<std::size_t>(c)] / per.g);
  const int b = Rf_isNull(base) ? -1 : clab::ref(base, g.lev, g.g, "base") - 1;
  SEXP w = PROTECT(real(6 * static_cast<R_xlen_t>(nc) + 9 * static_cast<R_xlen_t>(ncity) + 3 +
                        2 * static_cast<R_xlen_t>(nc) * g.g + 2 * g.g));
  ++np;
  double *co = REAL(w), *yo = co + 6 * static_cast<R_xlen_t>(nc), *sum = yo + 9 * ncity,
         *wk = sum + 3;
  kc(clab_tmpaud(INTEGER(pv), INTEGER(g.code), INTEGER(cell), n, g.g, nc, INTEGER(cy), ncity, b,
                 co, yo, sum, wk),
     "the kernel refused the input.");
  cols a(8, nc);
  SEXP cn = PROTECT(Rf_allocVector(STRSXP, nc)), pn = PROTECT(Rf_allocVector(STRSXP, nc));
  np += 3;
  for (int c = 0; c < nc; ++c) {
    const long long kk = u[static_cast<std::size_t>(c)];
    SET_STRING_ELT(cn, c, STRING_ELT(cl, static_cast<R_xlen_t>(kk / per.g)));
    SET_STRING_ELT(pn, c, STRING_ELT(per.lev, static_cast<R_xlen_t>(kk % per.g)));
  }
  a.put("city", cn);
  a.put("period", pn);
  int *cnn = a.in("n"), *ok = a.in("ok", LGLSXP);
  double *v[4] = {a.dbl("di"), a.dbl("dp"), a.dbl("gini"), a.dbl("bas")};
  for (int c = 0; c < nc; ++c) {
    cnn[c] = static_cast<int>(co[6 * c]);
    ok[c] = co[6 * c + 1] > 0.0;
    for (int t = 0; t < 4; ++t) v[t][c] = co[6 * c + 2 + t];
  }
  SEXP fc = PROTECT(a.done());
  cols e(10, ncity);
  np += 2;
  e.put("city", cl);
  int *ncl = e.in("ncell");
  static const char *const yn[] = {"di", "dp", "gini", "bas", "dimin", "dimax", "dirng"};
  double *yv[7];
  for (int t = 0; t < 7; ++t) yv[t] = e.dbl(yn[t]);
  int *gt = e.in("ngt1");
  for (int y = 0; y < ncity; ++y) {
    ncl[y] = static_cast<int>(yo[9 * y]);
    for (int t = 0; t < 7; ++t) yv[t][y] = yo[9 * y + 1 + t];
    gt[y] = static_cast<int>(yo[9 * y + 8]);
  }
  SEXP fy = PROTECT(e.done());
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  np += 2;
  SET_VECTOR_ELT(out, 0, fc);
  SET_VECTOR_ELT(out, 1, fy);
  static const char *const sn[] = {"rng", "spread"};
  SET_VECTOR_ELT(out, 2, named(2, sn, sum));
  SET_VECTOR_ELT(out, 3, Rf_ScalarString(STRING_ELT(g.lev, static_cast<R_xlen_t>(sum[2]))));
  static const char *const nm[] = {"cells", "cities", "sum", "base"};
  UNPROTECT(np + nb);
  return clab::nlist(out, nm);
}

SEXP C_clab_r_lipaud(SEXP h, SEXP d, SEXP pairs, SEXP lip) {
  clab::who = "lipaud";
  SEXP hv = PROTECT(num(h, "h"));
  const int n = len(hv, "h");
  const double L = REAL(PROTECT(num(lip, "lip", 1, 0.0)))[0];
  int np = 2, m;
  SEXP pi, pj, dv;
  if (Rf_isNull(pairs)) {
    // every pair i < j of a dist object or a square matrix, in the dist
    // object's order: j = 1 with i = 2, ..., n, then j = 2, ...
    if (Rf_isNull(d)) {
      fail("give `d`, the distances: a dist object or a square matrix, or "
           "`pairs` with one distance each.");
    }
    clab::unlab(d, "d");
    const bool dist = Rf_inherits(d, "dist");
    SEXP dm = Rf_getAttrib(d, R_DimSymbol);
    if (!clab::isnum(d) ||
        (!dist && !(TYPEOF(dm) == INTSXP && XLENGTH(dm) == 2 &&
                    INTEGER_RO(dm)[0] == n && INTEGER_RO(dm)[1] == n)) ||
        (dist && Rf_asInteger(Rf_getAttrib(d, Rf_install("Size"))) != n)) {
      fail("`d` must be a dist object or a square matrix of distances "
           "between the %d cases.", n);
    }
    if (n < 2) fail("`h` must hold at least two cases.");
    if (0.5 * n * (n - 1.0) > INT_MAX) fail("`h` holds too many cases for every pair.");
    m = n * (n - 1) / 2;
    SEXP dd = PROTECT(clab::dbl(d));
    pi = PROTECT(ints(m));
    pj = PROTECT(ints(m));
    dv = PROTECT(real(m));
    np += 4;
    int k = 0;
    for (int j = 0; j < n; ++j) {
      for (int i = j + 1; i < n; ++i, ++k) {
        INTEGER(pi)[k] = j;
        INTEGER(pj)[k] = i;
        double e;
        if (dist) {
          e = REAL_RO(dd)[k];
        } else {
          e = REAL_RO(dd)[i + static_cast<R_xlen_t>(n) * j];
          if (e != REAL_RO(dd)[j + static_cast<R_xlen_t>(n) * i]) fail("`d` must be symmetric.");
        }
        if (!(e >= 0.0) || !std::isfinite(e)) fail("`d` must hold finite distances, at least 0.");
        REAL(dv)[k] = e;
      }
    }
  } else {
    int r, c;
    SEXP pm = PROTECT(nmat(pairs, "pairs", &r, &c));
    ++np;
    if (c != 2 || r < 1) fail("`pairs` must be a matrix of two columns, one row per pair.");
    m = r;
    pi = PROTECT(ints(m));
    pj = PROTECT(ints(m));
    np += 2;
    for (int k = 0; k < 2 * m; ++k) {
      const double e = REAL(pm)[k];
      if (e != std::floor(e) || e < 1 || e > n) {
        fail("`pairs` must hold case numbers from 1 to %d.", n);
      }
      INTEGER(k < m ? pi : pj)[k % m] = static_cast<int>(e) - 1;
    }
    dv = Rf_isNull(d) ? R_NilValue : num(d, "d", m, 0.0);
    PROTECT(dv);
    ++np;
  }
  cols f(6, m);
  ++np;
  int *ci = f.in("i"), *cj = f.in("j");
  double *cd = f.dbl("d"), *gap = f.dbl("gap"), *rat = f.dbl("ratio");
  int *vi = f.in("viol", LGLSXP);
  double sum[6];
  kc(clab_lipaud(REAL(hv), n, INTEGER(pi), INTEGER(pj), Rf_isNull(dv) ? nullptr : REAL(dv), m, L,
                 gap, rat, vi, sum),
     "the kernel refused the input.");
  for (int k = 0; k < m; ++k) {
    ci[k] = INTEGER(pi)[k] + 1;
    cj[k] = INTEGER(pj)[k] + 1;
    cd[k] = Rf_isNull(dv) ? 1.0 : REAL(dv)[k];
  }
  SEXP fr = PROTECT(f.done());
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 7));
  np += 2;
  const int kr = static_cast<int>(sum[5]);
  SET_VECTOR_ELT(out, 0, Rf_ScalarReal(sum[0]));
  SET_VECTOR_ELT(out, 1, Rf_ScalarInteger(static_cast<int>(sum[1])));
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(sum[2]));
  SET_VECTOR_ELT(out, 3, Rf_ScalarLogical(sum[1] == 0.0));
  SET_VECTOR_ELT(out, 4, Rf_ScalarReal(sum[3]));
  SEXP wp = Rf_allocVector(INTSXP, 2);
  SET_VECTOR_ELT(out, 5, wp);
  INTEGER(wp)[0] = INTEGER(pi)[kr] + 1;
  INTEGER(wp)[1] = INTEGER(pj)[kr] + 1;
  SET_VECTOR_ELT(out, 6, fr);
  static const char *const nm[] = {"lreq", "nviol", "rate", "ok", "maxgap", "worst", "pairs"};
  UNPROTECT(np);
  return clab::nlist(out, nm);
}


SEXP C_clab_r_nsbnd(SEXP pobs, SEXP a, SEXP b) {
  clab::who = "nsbnd";
  SEXP pv = PROTECT(num(pobs, "pobs", -1, 0.0));
  const int n = len(pv, "pobs");
  for (int i = 0; i < n; ++i) {
    if (REAL(pv)[i] > 1.0) fail("`pobs` must be rates from 0 to 1.");
  }
  double alo, ahi, blo, bhi;
  box(a, "a", &alo, &ahi);
  box(b, "b", &blo, &bhi);
  if (!(ahi + bhi < 1.0)) fail("the upper ends of `a` and `b` must add to less than 1.");
  if (0.5 * n * (n - 1.0) > INT_MAX) fail("`pobs` holds too many rates for every pair.");
  const int m = n * (n - 1) / 2;
  SEXP lab = PROTECT(nmsof(pobs, n));
  cols f(5, n);
  f.put("rate", lab);
  double *po = f.dbl("pobs"), *lo = f.dbl("lo"), *hi = f.dbl("hi"), *wd = f.dbl("width");
  SEXP ord = PROTECT(ints(m)), bdv = PROTECT(real(m));
  kc(clab_nsbnd(REAL(pv), n, alo, ahi, blo, bhi, lo, hi, INTEGER(ord), REAL(bdv)),
     "the kernel refused the input.");
  for (int i = 0; i < n; ++i) {
    po[i] = REAL(pv)[i];
    wd[i] = hi[i] - lo[i];
  }
  SEXP fb = PROTECT(f.done());
  cols e(4, m);
  SEXP li = PROTECT(Rf_allocVector(STRSXP, m)), lj = PROTECT(Rf_allocVector(STRSXP, m)),
       os = PROTECT(Rf_allocVector(STRSXP, m));
  static const char *const ow[] = {"higher", "undecided", "lower"};
  int k = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j, ++k) {
      SET_STRING_ELT(li, k, STRING_ELT(lab, i));
      SET_STRING_ELT(lj, k, STRING_ELT(lab, j));
      const int o = INTEGER(ord)[k];
      SET_STRING_ELT(os, k, o == 2 ? NA_STRING : Rf_mkChar(ow[o + 1]));
    }
  }
  e.put("i", li);
  e.put("j", lj);
  e.put("ord", os);
  e.put("bd", bdv);
  SEXP fp = PROTECT(e.done());
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, fb);
  SET_VECTOR_ELT(out, 1, fp);
  static const char *const nm[] = {"bnd", "pairs"};
  UNPROTECT(12);
  return clab::nlist(out, nm);
}

SEXP C_clab_r_impfpr(SEXP p, SEXP ppv, SEXP fnr) {
  clab::who = "impfpr";
  SEXP uv[3] = {p, ppv, fnr};
  static const char *const un[] = {"p", "ppv", "fnr"};
  const int nb = csp::byname(uv, un, 3);
  SEXP v[3];
  R_xlen_t n = 1;
  for (int t = 0; t < 3; ++t) {
    v[t] = PROTECT(num(uv[t], un[t]));
    n = std::max(n, XLENGTH(v[t]));
  }
  for (int t = 0; t < 3; ++t) {
    if (XLENGTH(v[t]) != 1 && XLENGTH(v[t]) != n) {
      fail("`p`, `ppv` and `fnr` must have one length (or length 1).");
    }
  }
  const int m = clab::ilen(static_cast<double>(n), "p");
  SEXP x[3];
  for (int t = 0; t < 3; ++t) {
    x[t] = PROTECT(real(m));
    for (int i = 0; i < m; ++i) REAL(x[t])[i] = REAL(v[t])[XLENGTH(v[t]) == 1 ? 0 : i];
  }
  for (int i = 0; i < m; ++i) {
    if (!(REAL(x[0])[i] > 0.0 && REAL(x[0])[i] < 1.0)) fail("`p` must be base rates in (0, 1).");
    if (!(REAL(x[1])[i] > 0.0 && REAL(x[1])[i] <= 1.0)) fail("`ppv` must be in (0, 1].");
    if (!(REAL(x[2])[i] >= 0.0 && REAL(x[2])[i] < 1.0)) fail("`fnr` must be in [0, 1).");
  }
  SEXP out = PROTECT(real(m));
  kc(clab_impfpr(REAL(x[0]), REAL(x[1]), REAL(x[2]), m, REAL(out)), "the kernel refused the input.");
  for (int t = 0; t < 3; ++t) {
    SEXP nm = Rf_getAttrib(uv[t], R_NamesSymbol);
    if (XLENGTH(v[t]) == m && !Rf_isNull(nm)) {
      Rf_setAttrib(out, R_NamesSymbol, nm);
      break;
    }
  }
  UNPROTECT(7 + nb);
  return out;
}

SEXP C_clab_r_rnkres(SEXP est, SEXP hw) {
  clab::who = "rnkres";
  SEXP ev = PROTECT(num(est, "est"));
  const int n = len(ev, "est");
  if (n < 2) fail("`est` must hold at least two estimates.");
  if (static_cast<double>(n) * n > INT_MAX) fail("`est` holds too many estimates.");
  SEXP hv = PROTECT(num(hw, "hw", -1, 0.0));
  if (XLENGTH(hv) != 1 && XLENGTH(hv) != n) fail("`hw` must be one half-width or one per estimate.");
  SEXP hh = PROTECT(real(n));
  for (int i = 0; i < n; ++i) REAL(hh)[i] = REAL(hv)[XLENGTH(hv) == 1 ? 0 : i];
  SEXP id = PROTECT(Rf_allocMatrix(LGLSXP, n, n));
  double o[4];
  kc(clab_rnkres(REAL(ev), REAL(hh), n, LOGICAL(id), o), "the kernel refused the input.");
  for (R_xlen_t k = 0; k < static_cast<R_xlen_t>(n) * n; ++k) {
    if (LOGICAL(id)[k] < 0) LOGICAL(id)[k] = NA_LOGICAL;
  }
  SEXP lab = PROTECT(nmsof(est, n));
  SEXP dn = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dn, 0, lab);
  SET_VECTOR_ELT(dn, 1, lab);
  Rf_setAttrib(id, R_DimNamesSymbol, dn);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
  SET_VECTOR_ELT(out, 0, id);
  SET_VECTOR_ELT(out, 1, Rf_ScalarInteger(static_cast<int>(o[0])));
  SET_VECTOR_ELT(out, 2, Rf_ScalarInteger(static_cast<int>(o[1])));
  SET_VECTOR_ELT(out, 3, Rf_ScalarReal(o[2]));
  SET_VECTOR_ELT(out, 4, Rf_ScalarReal(o[3]));
  static const char *const nm[] = {"id", "npair", "nunid", "share", "res"};
  UNPROTECT(7);
  return clab::nlist(out, nm);
}

SEXP C_clab_r_calaud(SEXP area, SEXP risk, SEXP outcome, SEXP group,
                     SEXP population, SEXP nperm, SEXP seed) {
  clab::who = "calaud";
  if (Rf_inherits(area, "data.frame")) {
    if (!Rf_isNull(risk) || !Rf_isNull(outcome) || !Rf_isNull(group) || !Rf_isNull(population)) {
      fail("give `area` as a data frame of the columns, or every argument "
           "as a vector, not both.");
    }
    SEXP df = area;
    area = col(df, "area");
    risk = col(df, "risk");
    outcome = col(df, "outcome");
    group = col(df, "group");
    population = col(df, "population");
    if (Rf_isNull(area) || Rf_isNull(risk) || Rf_isNull(outcome) || Rf_isNull(group)) {
      fail("`area` as a data frame needs columns area, risk, outcome and group.");
    }
  }
  if (Rf_isNull(risk) || Rf_isNull(outcome) || Rf_isNull(group)) {
    fail("`risk`, `outcome` and `group` are needed.");
  }
  SEXP uv[4] = {area, risk, outcome, group};
  static const char *const un[] = {"area", "risk", "outcome", "group"};
  const int nb = csp::byname(uv, un, 4);
  const int n = len(uv[0], "area");
  if (n < 1) fail("`area` must hold at least one record.");
  lab a = glab(uv[0], "area", n);
  SEXP rv = PROTECT(num(uv[1], "risk", n));
  SEXP ov = PROTECT(num(uv[2], "outcome", n));
  lab g = glab(uv[3], "group", n);
  if (a.g < 3) fail("at least 3 areas are needed.");
  SEXP pp = R_NilValue;
  if (!Rf_isNull(population)) {
    SEXP pv = PROTECT(num(population, "population", -1, 0.0, true));
    SEXP pn = Rf_getAttrib(population, R_NamesSymbol);
    pp = PROTECT(real(n));
    if (!Rf_isNull(pn) && XLENGTH(pv) == a.g) {
      // one per area, matched by the area labels
      std::vector<int> at(static_cast<std::size_t>(a.g), -1);
      for (int k = 0; k < a.g; ++k) {
        for (int t = 0; t < a.g; ++t) {
          if (clab::seql(STRING_ELT(pn, k), STRING_ELT(a.lev, t))) {
            if (at[static_cast<std::size_t>(t)] >= 0) fail("`population` names an area twice.");
            at[static_cast<std::size_t>(t)] = k;
          }
        }
      }
      for (int t = 0; t < a.g; ++t) {
        if (at[static_cast<std::size_t>(t)] < 0) {
          fail("`population` has no value for area %s.", Rf_translateChar(STRING_ELT(a.lev, t)));
        }
      }
      for (int i = 0; i < n; ++i) {
        REAL(pp)[i] = REAL(pv)[at[static_cast<std::size_t>(INTEGER(a.code)[i])]];
      }
    } else if (XLENGTH(pv) == n) {
      for (int i = 0; i < n; ++i) REAL(pp)[i] = REAL(pv)[i];
    } else {
      fail("`population` must be one number per record, or one per area "
           "named by area.");
    }
    UNPROTECT(2);
    PROTECT(pp);
  } else {
    PROTECT(pp);
  }
  const int B = clab::whole(nperm, "nperm", 0);
  const unsigned int s = B > 0 ? sd(seed) : 0U;
  const int na = a.g, G = g.g;
  SEXP iw = PROTECT(ints(n + 2 * static_cast<R_xlen_t>(na)));
  SEXP w = PROTECT(real(7 * static_cast<R_xlen_t>(na) + 3 * static_cast<R_xlen_t>(G) + 7 +
                        na + 3 * static_cast<R_xlen_t>(G)));
  double *ar = REAL(w), *gr = ar + 7 * static_cast<R_xlen_t>(na), *st = gr + 3 * G, *wk = st + 7;
  const int rc = clab_calaud(INTEGER(a.code), REAL(rv), REAL(ov), INTEGER(g.code),
                             Rf_isNull(pp) ? nullptr : REAL(pp), n, na, G, B, s, ar, gr, st,
                             INTEGER(iw), wk);
  if (rc == CLAB_EDOM) fail("`population` must take one value within each area.");
  kc(rc, "the kernel refused the input.");
  cols fa(8, na);
  fa.put("area", a.lev);
  int *an = fa.in("n");
  double *ac[2] = {fa.dbl("risk"), fa.dbl("rate")};
  SEXP ag = PROTECT(Rf_allocVector(STRSXP, na));
  fa.put("group", ag);
  double *ar2[3] = {fa.dbl("riskrk"), fa.dbl("raterk"), fa.dbl("gap")};
  for (int k = 0; k < na; ++k) {
    an[k] = static_cast<int>(ar[k]);
    ac[0][k] = ar[na + k];
    ac[1][k] = ar[2 * na + k];
    SET_STRING_ELT(ag, k, STRING_ELT(g.lev, static_cast<R_xlen_t>(ar[3 * na + k])));
    for (int t = 0; t < 3; ++t) ar2[t][k] = ar[(4 + t) * na + k];
  }
  SEXP f1 = PROTECT(fa.done());
  cols fg(4, G);
  fg.put("group", g.lev);
  int *gn = fg.in("n");
  double *gm = fg.dbl("gap"), *gp = fg.dbl("p");
  for (int j = 0; j < G; ++j) {
    gn[j] = static_cast<int>(gr[j]);
    gm[j] = gr[G + j];
    gp[j] = gr[2 * G + j];
  }
  SEXP f2 = PROTECT(fg.done());
  // the large-sample p-value: t = rho sqrt((n - 2) / (1 - rho^2)) on n - 2
  // degrees of freedom, two-sided, as cor.test() gives it without exact
  const double rho = st[1], df = na - 2.0;
  const double tt = rho * std::sqrt(df / ((1.0 - rho) * (1.0 + rho)));
  const double sp[5] = {rho, st[2], st[3], tt,
                        std::isnan(tt) ? qnan : std::fmin(1.0, 2.0 * Rf_pt(-std::fabs(tt), df, 1, 0))};
  static const char *const spn[] = {"rho", "s", "pex", "t", "pt"};
  static const char *const bn[] = {"ss", "p"};
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 6));
  SET_VECTOR_ELT(out, 0, f1);
  SET_VECTOR_ELT(out, 1, f2);
  SET_VECTOR_ELT(out, 2, named(5, spn, sp));
  SET_VECTOR_ELT(out, 3, Rf_ScalarLogical(st[4] > 0.0));
  SET_VECTOR_ELT(out, 4, named(2, bn, st + 5));
  SET_VECTOR_ELT(out, 5, Rf_ScalarInteger(B));
  static const char *const nm[] = {"areas", "groups", "spear", "ties", "between", "nperm"};
  UNPROTECT(15 + nb);
  return clab::nlist(out, nm);
}

SEXP C_clab_r_scrdsp(SEXP score, SEXP group, SEXP ref, SEXP conf) {
  clab::who = "scrdsp";
  SEXP uv[2] = {score, group};
  static const char *const un[] = {"score", "group"};
  const int nb = csp::byname(uv, un, 2);
  SEXP xv = PROTECT(num(uv[0], "score"));
  const int n = len(xv, "score");
  lab g = glab(uv[1], "group", n);
  if (g.g < 2) fail("`group` must hold at least two groups.");
  const int r = Rf_isNull(ref) ? -1 : clab::ref(ref, g.lev, g.g) - 1;
  const double cf = clab::one(conf, "conf", 0.0, 1.0);
  const int G = g.g;
  SEXP iw = PROTECT(ints(n));
  SEXP w = PROTECT(real(6 * static_cast<R_xlen_t>(G) + 9));
  double *gs = REAL(w), *st = gs + 6 * G;
  kc(clab_scrdsp(REAL(xv), INTEGER(g.code), n, G, r, gs, st, INTEGER(iw)),
     "the kernel refused the input.");
  const double q = st[2] > 0.0 ? Rf_qt(0.5 + 0.5 * cf, st[2], 1, 0) : qnan;
  cols f(9, G);
  f.put("group", g.lev);
  int *gn = f.in("n");
  static const char *const cn[] = {"mean", "median", "sd", "gap", "se"};
  double *c[5];
  for (int t = 0; t < 5; ++t) c[t] = f.dbl(cn[t]);
  double *lo = f.dbl("lo"), *hi = f.dbl("hi");
  for (int j = 0; j < G; ++j) {
    gn[j] = static_cast<int>(gs[j]);
    for (int t = 0; t < 5; ++t) c[t][j] = gs[(t + 1) * G + j];
    lo[j] = c[3][j] - q * c[4][j];
    hi[j] = c[3][j] + q * c[4][j];
  }
  SEXP fr = PROTECT(f.done());
  // upper tails of F; an infinite F (no variation within groups) gives 0
  auto pf = [](double F, double d1, double d2) {
    return std::isnan(F) ? qnan : Rf_pf(F, d1, d2, 0, 0);
  };
  const double an[4] = {st[0], st[1], st[2], pf(st[0], st[1], st[2])};
  const double wl[4] = {st[3], st[4], st[5], pf(st[3], st[4], st[5])};
  static const char *const fn[] = {"f", "df1", "df2", "p"};
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
  SET_VECTOR_ELT(out, 0, fr);
  SET_VECTOR_ELT(out, 1, named(4, fn, an));
  SET_VECTOR_ELT(out, 2, named(4, fn, wl));
  SET_VECTOR_ELT(out, 3, Rf_ScalarReal(st[6]));
  SET_VECTOR_ELT(out, 4, Rf_ScalarString(STRING_ELT(g.lev, static_cast<R_xlen_t>(st[7]))));
  static const char *const nm[] = {"tab", "anova", "welch", "spread", "ref"};
  UNPROTECT(8 + nb);
  return clab::nlist(out, nm);
}

SEXP C_clab_r_nordet(SEXP x, SEXP off, SEXP radius, SEXP p, SEXP draw, SEXP seed) {
  clab::who = "nordet";
  SEXP cx, cy, ox, oy;
  const int n = pts2(x, "x", &cx, &cy, false);
  PROTECT(cx);
  PROTECT(cy);
  const int m = pts2(off, "off", &ox, &oy, true);
  PROTECT(ox);
  PROTECT(oy);
  const double r = REAL(PROTECT(num(radius, "radius", 1, 0.0, true)))[0];
  const double pd = clab::val1(p);
  if (!(pd > 0.0 && pd <= 1.0)) fail("`p` must be one number in (0, 1].");
  const bool dr = csp::flag(draw, "draw");
  const unsigned int s = dr ? sd(seed) : 0U;
  cols f(3, n);
  int *k = f.in("k");
  double *pr = f.dbl("prob");
  int *det = dr ? f.in("det") : nullptr;
  kc(clab_nordet(REAL(cx), REAL(cy), n, m > 0 ? REAL(ox) : nullptr, m > 0 ? REAL(oy) : nullptr, m,
                 r, pd, dr ? 1 : 0, s, k, pr, det),
     "the kernel refused the input.");
  SEXP out = f.done();
  UNPROTECT(6);
  return out;
}

SEXP C_clab_r_fairsim(SEXP n, SEXP groups, SEXP props, SEXP nareas, SEXP base,
                      SEXP bias, SEXP seed) {
  clab::who = "fairsim";
  const int N = clab::whole(n, "n", 1);
  clab::unlab(groups, "groups");
  SEXP gl = PROTECT(clab::labs(groups));
  if (Rf_isNull(gl) || TYPEOF(groups) == VECSXP || clab::anyna(gl) || XLENGTH(gl) < 2) {
    fail("`groups` must be two or more labels with no missing values.");
  }
  const int G = len(gl, "groups");
  for (int i = 0; i < G; ++i) {
    for (int j = 0; j < i; ++j) {
      if (clab::seql(STRING_ELT(gl, i), STRING_ELT(gl, j))) fail("`groups` must not repeat a label.");
    }
  }
  SEXP pv = R_NilValue;
  if (!Rf_isNull(props)) {
    pv = num(props, "props", G, 0.0);
  }
  PROTECT(pv);
  if (!Rf_isNull(pv)) {
    double t = 0.0;
    for (int j = 0; j < G; ++j) t += REAL(pv)[j];
    if (!(t > 0.0)) fail("`props` must not all be 0.");
  }
  const int na = clab::whole(nareas, "nareas", 1);
  if (na < G) fail("`nareas` must be at least the number of groups, %d.", G);
  const double b0 = clab::p01(base, "base");
  const double bs = clab::val1(bias);
  if (!(bs >= -1.0 && bs <= 1.0)) fail("`bias` must be one number from -1 to 1.");
  const unsigned int s = sd(seed);
  SEXP ai = PROTECT(ints(N)), gi = PROTECT(ints(N));
  cols f(5, N);
  int *ar = f.in("area");
  SEXP gs = PROTECT(Rf_allocVector(STRSXP, N));
  f.put("group", gs);
  int *tr = f.in("truth"), *dt = f.in("det");
  double *rk = f.dbl("risk");
  kc(clab_fairsim(N, G, Rf_isNull(pv) ? nullptr : REAL(pv), na, b0, bs, s, INTEGER(ai), INTEGER(gi),
                  tr, dt, rk),
     "the kernel refused the input.");
  for (int i = 0; i < N; ++i) {
    ar[i] = INTEGER(ai)[i] + 1;
    SET_STRING_ELT(gs, i, STRING_ELT(gl, INTEGER(gi)[i]));
  }
  SEXP out = f.done();
  UNPROTECT(6);
  return out;
}

SEXP C_clab_r_oaxbd(SEXP y, SEXP x, SEXP grp, SEXP ref, SEXP conf) {
  clab::who = "oaxbd";
  SEXP uv[2] = {y, grp};
  static const char *const un[] = {"y", "grp"};
  const int nb = csp::byname(uv, un, 2);
  SEXP yv = PROTECT(num(uv[0], "y"));
  const int n = len(yv, "y");
  lab g = glab(uv[1], "grp", n);
  if (g.g != 2) fail("`grp` must hold exactly two groups.");
  const int r = clab::ref(ref, g.lev, 2) - 1;
  const double z = zconf(conf);
  // the covariates: none, a vector, a numeric matrix, or a data frame or
  // list of numeric columns
  int p = 0;
  SEXP xv = R_NilValue, tn = R_NilValue;
  static const char *const fmsg =
      "`x` must be numeric covariates, one row per case: code a factor as "
      "indicator columns (how the unexplained part splits across them then "
      "depends on the category left out).";
  if (!Rf_isNull(x)) {
    clab::unlab(x, "x");
    SEXP dm = Rf_getAttrib(x, R_DimSymbol);
    if (TYPEOF(x) == VECSXP) {
      p = len(x, "x");
      xv = PROTECT(real(static_cast<R_xlen_t>(n) * p));
      for (int j = 0; j < p; ++j) {
        SEXP c = VECTOR_ELT(x, j);
        if (!clab::isnum(c) || Rf_isFactor(c) || XLENGTH(c) != n) fail("%s", fmsg);
        clab::unlab(c, "x");
        SEXP cd = PROTECT(clab::dbl(c));
        for (int i = 0; i < n; ++i) REAL(xv)[i + static_cast<R_xlen_t>(n) * j] = REAL(cd)[i];
        UNPROTECT(1);
      }
      tn = Rf_getAttrib(x, R_NamesSymbol);
    } else if (clab::isnum(x) && !Rf_isFactor(x)) {
      if (TYPEOF(dm) == INTSXP && XLENGTH(dm) == 2) {
        if (INTEGER_RO(dm)[0] != n) fail("%s", fmsg);
        p = INTEGER_RO(dm)[1];
        SEXP dn = Rf_getAttrib(x, R_DimNamesSymbol);
        if (!Rf_isNull(dn)) tn = VECTOR_ELT(dn, 1);
      } else {
        if (XLENGTH(x) != n) fail("%s", fmsg);
        p = 1;
      }
      xv = PROTECT(clab::dcopy(PROTECT(clab::dbl(x))));
      UNPROTECT(2);
      PROTECT(xv);
    } else {
      fail("%s", fmsg);
    }
    for (R_xlen_t k = 0; k < static_cast<R_xlen_t>(n) * p; ++k) {
      if (!std::isfinite(REAL(xv)[k])) fail("`x` must hold finite numbers with no missing values.");
    }
  } else {
    PROTECT(xv);
  }
  const int q = p + 1;
  int nn[2] = {0, 0};
  for (int i = 0; i < n; ++i) ++nn[INTEGER(g.code)[i]];
  if (nn[0] <= q || nn[1] <= q) {
    fail("each group needs more cases than coefficients, %d.", q);
  }
  // the term labels: the constant, then the covariates' names or x1, x2, ...
  SEXP lab = PROTECT(Rf_allocVector(STRSXP, q));
  SET_STRING_ELT(lab, 0, Rf_mkChar("const"));
  for (int j = 0; j < p; ++j) {
    char b[24];
    if (!Rf_isNull(tn) && STRING_ELT(tn, j) != NA_STRING && CHAR(STRING_ELT(tn, j))[0] != '\0') {
      SET_STRING_ELT(lab, j + 1, STRING_ELT(tn, j));
    } else {
      std::snprintf(b, sizeof b, p == 1 && Rf_isNull(tn) ? "x" : "x%d", j + 1);
      SET_STRING_ELT(lab, j + 1, Rf_mkChar(b));
    }
  }
  SEXP w = PROTECT(real(2 * static_cast<R_xlen_t>(q) * 3 + 14 + 5 * static_cast<R_xlen_t>(q) +
                        static_cast<R_xlen_t>(q) * (n + 5 * q + 1) + n));
  double *xm = REAL(w), *bb = xm + 2 * q, *bse = bb + 2 * q, *dec = bse + 2 * q, *var = dec + 14,
         *wk = var + 5 * q;
  const int rc = clab_oaxbd(REAL(yv), Rf_isNull(xv) ? nullptr : REAL(xv), INTEGER(g.code), n, p, r,
                            xm, bb, bse, dec, var, wk);
  if (rc == CLAB_EDOM) fail("a covariate is constant, or a combination of the others, within a group.");
  kc(rc, "the kernel refused the input.");
  SEXP gl = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(gl, 0, STRING_ELT(g.lev, 1 - r));
  SET_STRING_ELT(gl, 1, STRING_ELT(g.lev, r));
  static const char *const gn[] = {"a", "b"};
  Rf_setAttrib(gl, R_NamesSymbol, PROTECT(Rf_allocVector(STRSXP, 2)));
  for (int k = 0; k < 2; ++k) SET_STRING_ELT(Rf_getAttrib(gl, R_NamesSymbol), k, Rf_mkChar(gn[k]));
  cols fm(3, q);
  fm.put("term", lab);
  double *ma = fm.dbl("a"), *mb = fm.dbl("b");
  cols fc(5, q);
  fc.put("term", lab);
  double *ca = fc.dbl("a"), *sa = fc.dbl("ase"), *cb = fc.dbl("b"), *sb = fc.dbl("bse");
  cols fv(6, q);
  fv.put("term", lab);
  static const char *const vn[] = {"exp", "expse", "unexp", "unexpse", "shift"};
  double *vc[5];
  for (int t = 0; t < 5; ++t) vc[t] = fv.dbl(vn[t]);
  for (int j = 0; j < q; ++j) {
    ma[j] = xm[j];
    mb[j] = xm[q + j];
    ca[j] = bb[j];
    sa[j] = bse[j];
    cb[j] = bb[q + j];
    sb[j] = bse[q + j];
    for (int t = 0; t < 5; ++t) vc[t][j] = var[t * q + j];
  }
  // the parts: the gap, twofold at B's coefficients, twofold at A's, the
  // threefold (explained at B's, unexplained at B's means, interaction)
  static const char *const pn[] = {"gap", "exp", "unexp", "expa", "unexpa", "endow", "coef", "inter"};
  const int ix[8] = {2, 4, 6, 8, 10, 4, 10, 12};
  cols fd(5, 8);
  SEXP ps = PROTECT(Rf_allocVector(STRSXP, 8));
  for (int k = 0; k < 8; ++k) SET_STRING_ELT(ps, k, Rf_mkChar(pn[k]));
  fd.put("part", ps);
  double *de = fd.dbl("est"), *ds = fd.dbl("se"), *dl = fd.dbl("lo"), *dh = fd.dbl("hi");
  for (int k = 0; k < 8; ++k) {
    de[k] = dec[ix[k]];
    ds[k] = dec[ix[k] + 1];
    dl[k] = de[k] - z * ds[k];
    dh[k] = de[k] + z * ds[k];
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
  SET_VECTOR_ELT(out, 0, gl);
  SET_VECTOR_ELT(out, 1, fm.done());
  SET_VECTOR_ELT(out, 2, fc.done());
  SET_VECTOR_ELT(out, 3, fd.done());
  SET_VECTOR_ELT(out, 4, fv.done());
  static const char *const nm[] = {"groups", "means", "coef", "dec", "by"};
  UNPROTECT(14 + nb);
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
  int r = 0, c = 0;
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
