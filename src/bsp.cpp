// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of crimlab's spatial functions, one per exported
// function: each reads its points through csp::xy() (spin.h), checks the
// other arguments, calls one kernel and builds the R result. Results that
// describe points keep the input order. Nothing calls back into the R
// evaluator.
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstring>
#include <limits>

#define R_NO_REMAP_RMATH
#include <Rmath.h>

#include "crimlab.h"
#include "rchk.h"
#include "spin.h"

using clab::fail;
using clab::num;

namespace {

constexpr double pi = 3.14159265358979323846;
const double qnan = std::numeric_limits<double>::quiet_NaN();

SEXP real(R_xlen_t n) { return Rf_allocVector(REALSXP, n); }
SEXP ints(R_xlen_t n) { return Rf_allocVector(INTSXP, n); }

// a kernel's code: CLAB_OK passes; otherwise the message msg, or the
// kernel's own words when msg is null
void kc(int rc, const char *msg = nullptr) {
  if (rc == CLAB_OK) return;
  if (msg == nullptr) msg = clab_strerror(rc);
  fail("%s", msg);
}

// weights: NULL for none, otherwise n finite numbers >= 0 with a positive
// sum
const double *wts(SEXP w, int n, SEXP *keep) {
  if (Rf_isNull(w)) return nullptr;
  *keep = PROTECT(num(w, "w", n, 0.0));
  const double *p = REAL_RO(*keep);
  bool any = false;
  for (int i = 0; i < n && !any; ++i) any = p[i] > 0;
  if (!any) fail("`w` must have a positive sum.");
  UNPROTECT(1);
  return p;
}

double wsum(const double *w, int n) {
  if (w == nullptr) return n;
  double s = 0;
  for (int i = 0; i < n; ++i) s += w[i];
  return s;
}

// the divisor offset k in {0, 1, 2}, below the total weight
int kdiv(SEXP k, double W) {
  const int v = clab::whole(k, "k", 0);
  if (v > 2) fail("`k` must be 0, 1 or 2.");
  if (!(W > v)) fail("the total weight must exceed `k`.");
  return v;
}

// a positive finite number
double posn(SEXP x, const char *nm) {
  return clab::one(x, nm, 0.0, std::numeric_limits<double>::infinity());
}

// one number from lo to hi
double in(SEXP x, const char *nm, double lo, double hi) {
  const double v = clab::val1(x);
  if (!(v >= lo && v <= hi)) fail("`%s` must be one number from %g to %g.", nm, lo, hi);
  return v;
}

// a named list of the elements of v
SEXP lst(SEXP v, const char *const *nm) { return clab::nlist(v, nm); }

// a data frame of k double columns of length n, named nm; col[j] gives
// each column's storage
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

}  // namespace

extern "C" {

// ------------------------------------------------------------ geography

SEXP C_clab_r_hvdist(SEXP from, SEXP to, SEXP r) {
  clab::who = "hvdist";
  const double rad = posn(r, "r");
  SEXP t = PROTECT(Rf_ScalarLogical(TRUE));
  csp::pts a = csp::xy(from, R_NilValue, t, csp::raw);
  csp::pts b = csp::xy(to, R_NilValue, t, csp::raw);
  if (a.n != b.n && a.n != 1 && b.n != 1) {
    fail("`from` and `to` must have the same number of points, or one.");
  }
  const int n = std::max(a.n, b.n);
  SEXP d = PROTECT(real(n));
  SEXP rep = PROTECT(real(4 * static_cast<R_xlen_t>(n)));
  double *q = REAL(rep);
  for (int i = 0; i < n; ++i) {
    q[i] = a.x[a.n == 1 ? 0 : i];
    q[n + i] = a.y[a.n == 1 ? 0 : i];
    q[2 * n + i] = b.x[b.n == 1 ? 0 : i];
    q[3 * n + i] = b.y[b.n == 1 ? 0 : i];
  }
  kc(clab_hvdist(q, q + n, q + 2 * n, q + 3 * n, n, rad, REAL(d)));
  UNPROTECT(5);
  return d;
}

SEXP C_clab_r_eqproj(SEXP x, SEXP y, SEXP lon0, SEXP lat0, SEXP r) {
  clab::who = "eqproj";
  const double rad = posn(r, "r");
  SEXP t = PROTECT(Rf_ScalarLogical(TRUE));
  csp::pts p = csp::xy(x, y, t, csp::raw);
  double c0, c1;
  if (Rf_isNull(lon0) != Rf_isNull(lat0)) {
    fail("give both `lon0` and `lat0`, or neither.");
  }
  if (Rf_isNull(lon0)) {
    c0 = p.lon0;
    c1 = p.lat0;
  } else {
    c0 = in(lon0, "lon0", -180.0, 180.0);
    c1 = in(lat0, "lat0", -90.0, 90.0);
  }
  static const char *const nm[] = {"x", "y"};
  double *col[2];
  SEXP f = PROTECT(dframe(2, p.n, nm, col));
  kc(clab_eqproj(p.x, p.y, p.n, c0, c1, rad, col[0], col[1]),
     "a point lies at the antipode of the centre.");
  SEXP at = PROTECT(Rf_allocVector(REALSXP, 2));
  REAL(at)[0] = c0;
  REAL(at)[1] = c1;
  Rf_setAttrib(f, Rf_install("centre"), at);
  UNPROTECT(4);
  return f;
}

SEXP C_clab_r_eqinv(SEXP x, SEXP y, SEXP lon0, SEXP lat0, SEXP r) {
  clab::who = "eqinv";
  const double rad = posn(r, "r");
  csp::pts p = csp::xy(x, y, R_NilValue, csp::planar);
  const double c0 = in(lon0, "lon0", -180.0, 180.0);
  const double c1 = in(lat0, "lat0", -90.0, 90.0);
  static const char *const nm[] = {"lon", "lat"};
  double *col[2];
  SEXP f = PROTECT(dframe(2, p.n, nm, col));
  kc(clab_eqinv(p.x, p.y, p.n, c0, c1, rad, col[0], col[1]),
     "a point lies farther than 2 r from the centre.");
  UNPROTECT(2);
  return f;
}

// ---------------------------------------------------------- centrography

//' @srrstats {G3.1} The covariances of cgstat() and sdell() are computed
//'   in C++ (compensated sums about the refined mean, cgeo.cpp), never by
//'   stats::cov(), with the divisor chosen by `k`.
//' @srrstats {G3.1a} The help pages of cgstat() and sdell() document `k`,
//'   and the examples of cgstat() use it.
//' @noRd
SEXP C_clab_r_cgstat(SEXP x, SEXP y, SEXP w, SEXP k, SEXP lonlat) {
  clab::who = "cgstat";
  SEXP uv[3] = {x, y, w};
  static const char *const un[] = {"x", "y", "w"};
  const int nb = csp::byname(uv, un, 3);
  x = uv[0];
  y = uv[1];
  w = uv[2];
  csp::pts p = csp::xy(x, y, lonlat);
  SEXP keep = R_NilValue;
  const double *wp = wts(w, p.n, &keep);
  PROTECT(keep);
  const int kk = kdiv(k, wsum(wp, p.n));
  double o[10];
  kc(clab_cgstat(p.x, p.y, wp, p.n, kk, o));
  if (p.ll) {
    // the means of projected coordinates depend on their origin: not given
    csp::back(p, o + 1, o + 2, 1);
    for (int j = 6; j < 10; ++j) o[j] = qnan;
  }
  static const char *const nm[] = {"w", "mx", "my", "sx", "sy", "sdd",
                                   "gx", "gy", "hx", "hy"};
  UNPROTECT(2 + nb);
  return clab::nvec(10, nm, o);
}

SEXP C_clab_r_cgmedn(SEXP x, SEXP y, SEXP lonlat) {
  clab::who = "cgmedn";
  csp::pts p = csp::xy(x, y, lonlat);
  SEXP work = PROTECT(real(p.n));
  double o[2];
  kc(clab_cgmedn(p.x, p.y, p.n, REAL(work), o));
  csp::back(p, o, o + 1, 1);
  static const char *const nm[] = {"x", "y"};
  UNPROTECT(2);
  return clab::nvec(2, nm, o);
}

SEXP C_clab_r_webpt(SEXP x, SEXP y, SEXP w, SEXP start, SEXP tol,
                    SEXP maxit, SEXP lonlat) {
  clab::who = "webpt";
  SEXP uv[3] = {x, y, w};
  static const char *const un[] = {"x", "y", "w"};
  const int nb = csp::byname(uv, un, 3);
  x = uv[0];
  y = uv[1];
  w = uv[2];
  csp::pts p = csp::xy(x, y, lonlat);
  SEXP keep = R_NilValue;
  const double *wp = wts(w, p.n, &keep);
  PROTECT(keep);
  double s[2];
  const bool has = !Rf_isNull(start);
  if (has) {
    SEXP sv = PROTECT(num(start, "start", 2));
    s[0] = REAL_RO(sv)[0];
    s[1] = REAL_RO(sv)[1];
    UNPROTECT(1);
    csp::fwd(p, s, s + 1, 1);
  }
  const double tl = posn(tol, "tol");
  const int mi = clab::whole(maxit, "maxit", 0);
  double o[6];
  kc(clab_webpt(p.x, p.y, wp, p.n, has ? s : nullptr, tl, mi, o));
  csp::back(p, o, o + 1, 1);
  static const char *const nm[] = {"x", "y", "obj", "gap", "iter", "code"};
  UNPROTECT(2 + nb);
  return clab::nvec(6, nm, o);
}

SEXP C_clab_r_sdell(SEXP x, SEXP y, SEXP w, SEXP k, SEXP c, SEXP lonlat) {
  clab::who = "sdell";
  SEXP uv[3] = {x, y, w};
  static const char *const un[] = {"x", "y", "w"};
  const int nb = csp::byname(uv, un, 3);
  x = uv[0];
  y = uv[1];
  w = uv[2];
  csp::pts p = csp::xy(x, y, lonlat);
  SEXP keep = R_NilValue;
  const double *wp = wts(w, p.n, &keep);
  PROTECT(keep);
  const int kk = kdiv(k, wsum(wp, p.n));
  const double cc = posn(c, "c");
  double o[6];
  kc(clab_sdell(p.x, p.y, wp, p.n, kk, cc, o));
  csp::back(p, o, o + 1, 1);
  static const char *const nm[] = {"mx", "my", "a", "b", "angle", "area"};
  UNPROTECT(2 + nb);
  return clab::nvec(6, nm, o);
}

SEXP C_clab_r_cvhull(SEXP x, SEXP y, SEXP lonlat) {
  clab::who = "cvhull";
  csp::pts p = csp::xy(x, y, lonlat);
  SEXP iw = PROTECT(ints(p.n));
  SEXP h = PROTECT(ints(2 * static_cast<R_xlen_t>(p.n)));
  int nh = 0;
  double o[2];
  kc(clab_cvhull(p.x, p.y, p.n, INTEGER(iw), INTEGER(h), &nh, o),
     "coordinates must be 0 or of magnitude from 1e-140 to 1e150.");
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP id = ints(nh);
  SET_VECTOR_ELT(out, 0, id);
  for (int i = 0; i < nh; ++i) INTEGER(id)[i] = INTEGER(h)[i] + 1;
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(o[0]));
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(o[1]));
  static const char *const nm[] = {"idx", "area", "perim"};
  UNPROTECT(4);
  return lst(out, nm);
}

}  // extern "C"

namespace {

// planar points inside the rectangle w, at least two of them
csp::pts inwin(SEXP x, SEXP y, SEXP w, double *wv) {
  csp::pts p = csp::xy(x, y, R_NilValue, csp::planar);
  csp::win(w, wv);
  if (p.n < 2) fail("at least 2 points are needed.");
  for (int i = 0; i < p.n; ++i) {
    if (!csp::inside(wv, p.x[i], p.y[i])) {
      fail("every point must lie inside `win`.");
    }
  }
  return p;
}

// radii: non-negative and strictly increasing
SEXP radii(SEXP r) {
  SEXP v = PROTECT(num(r, "r", -1, 0.0));
  const double *p = REAL_RO(v);
  for (R_xlen_t j = 1; j < XLENGTH(v); ++j) {
    if (!(p[j] > p[j - 1])) fail("`r` must be strictly increasing.");
  }
  UNPROTECT(1);
  return v;
}

const char *const kcor[] = {"none", "border", "modif", "iso"};

}  // namespace

extern "C" {

// ------------------------------------------------------- distance analysis

SEXP C_clab_r_nnidx(SEXP x, SEXP y, SEXP win, SEXP kmax, SEXP nsim,
                    SEXP seed, SEXP dist) {
  clab::who = "nnidx";
  double w[4];
  csp::pts p = inwin(x, y, win, w);
  const int n = p.n;
  const int km = clab::whole(kmax, "kmax", 1);
  if (km > n - 1) fail("`kmax` must be below the number of points.");
  const int ns = clab::whole(nsim, "nsim", 0);
  const bool dd = csp::flag(dist, "dist");
  const unsigned int sd = ns > 0 ? csp::seed(seed) : 0u;
  SEXP iw = PROTECT(ints(2 * static_cast<R_xlen_t>(n) + 2 + 2 * km));
  SEXP wk = PROTECT(real(3 * static_cast<R_xlen_t>(km) +
                         (ns > 0 ? 2 * static_cast<R_xlen_t>(n) : 0)));
  SEXP dv = PROTECT(dd ? real(static_cast<R_xlen_t>(n) * km) : R_NilValue);
  SEXP o = PROTECT(real(6 * static_cast<R_xlen_t>(km) + 4));
  kc(clab_nnidx(p.x, p.y, n, w, km, ns, sd, dd ? REAL(dv) : nullptr,
                INTEGER(iw), REAL(wk), REAL(o)));
  static const char *const on[] = {"k", "mean", "exp", "ratio", "pclu",
                                   "preg", "sd"};
  double *c[7];
  SEXP ord = PROTECT(dframe(7, km, on, c));
  for (int k = 0; k < km; ++k) {
    c[0][k] = k + 1;
    for (int j = 0; j < 6; ++j) c[j + 1][k] = REAL(o)[6 * k + j];
  }
  static const char *const cn[] = {"donnelly", "ratio", "z", "zc"};
  SEXP ce = PROTECT(clab::nvec(4, cn, REAL(o) + 6 * km));
  SEXP dm = R_NilValue;
  if (dd) {
    dm = PROTECT(Rf_allocMatrix(REALSXP, n, km));
    for (int i = 0; i < n; ++i) {
      for (int k = 0; k < km; ++k) {
        REAL(dm)[i + static_cast<R_xlen_t>(n) * k] =
            REAL(dv)[static_cast<R_xlen_t>(i) * km + k];
      }
    }
  } else {
    PROTECT(dm);
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, ord);
  SET_VECTOR_ELT(out, 1, ce);
  SET_VECTOR_ELT(out, 2, dm);
  static const char *const nm[] = {"order", "ce", "dist"};
  UNPROTECT(9);
  return lst(out, nm);
}

SEXP C_clab_r_ripk(SEXP x, SEXP y, SEXP win, SEXP r, SEXP corr) {
  clab::who = "ripk";
  double w[4];
  csp::pts p = inwin(x, y, win, w);
  SEXP rv = PROTECT(radii(r));
  const int m = clab::ilen(static_cast<double>(XLENGTH(rv)), "r");
  const int co = csp::pick(corr, "corr", kcor, 4);
  SEXP iw = PROTECT(ints(2 * static_cast<R_xlen_t>(p.n) + 2));
  SEXP wk = PROTECT(real(4 * static_cast<R_xlen_t>(m)));
  static const char *const nm[] = {"r", "k", "l"};
  double *c[3];
  SEXP f = PROTECT(dframe(3, m, nm, c));
  kc(clab_ripk(p.x, p.y, p.n, w, REAL_RO(rv), m, co, INTEGER(iw), REAL(wk),
               c[1]));
  for (int j = 0; j < m; ++j) {
    c[0][j] = REAL_RO(rv)[j];
    c[2][j] = std::sqrt(c[1][j] / pi);
  }
  UNPROTECT(5);
  return f;
}

SEXP C_clab_r_ripenv(SEXP x, SEXP y, SEXP win, SEXP r, SEXP corr, SEXP nsim,
                     SEXP seed) {
  clab::who = "ripenv";
  double w[4];
  csp::pts p = inwin(x, y, win, w);
  SEXP rv = PROTECT(radii(r));
  const int m = clab::ilen(static_cast<double>(XLENGTH(rv)), "r");
  const int co = csp::pick(corr, "corr", kcor, 4);
  const int ns = clab::whole(nsim, "nsim", 1);
  const unsigned int sd = csp::seed(seed);
  SEXP iw = PROTECT(ints(2 * static_cast<R_xlen_t>(p.n) + 2));
  SEXP wk = PROTECT(real(5 * static_cast<R_xlen_t>(m) + 2 * static_cast<R_xlen_t>(p.n)));
  SEXP o = PROTECT(real(4 * static_cast<R_xlen_t>(m) + 2));
  kc(clab_ripenv(p.x, p.y, p.n, w, REAL_RO(rv), m, co, ns, sd, INTEGER(iw),
                 REAL(wk), REAL(o)));
  static const char *const en[] = {"r", "obs", "lo", "hi", "mean"};
  double *c[5];
  SEXP env = PROTECT(dframe(5, m, en, c));
  for (int j = 0; j < m; ++j) {
    c[0][j] = REAL_RO(rv)[j];
    for (int b = 0; b < 4; ++b) c[b + 1][j] = REAL(o)[b * m + j];
  }
  static const char *const tn[] = {"t", "p"};
  SEXP te = PROTECT(clab::nvec(2, tn, REAL(o) + 4 * m));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, env);
  SET_VECTOR_ELT(out, 1, te);
  static const char *const nm[] = {"env", "test"};
  UNPROTECT(8);
  return lst(out, nm);
}

}  // extern "C"

namespace {

// group numbers from 0 (or -1) to R's from 1 (or NA)
SEXP grpno(const int *cl, R_xlen_t n) {
  SEXP v = ints(n);
  int *o = INTEGER(v);
  for (R_xlen_t i = 0; i < n; ++i) o[i] = cl[i] >= 0 ? cl[i] + 1 : NA_INTEGER;
  return v;
}

// back-projection of the centre columns a, b of a table of m rows (rows
// with NaN stay NaN)
void backcols(const csp::pts &p, double *a, double *b, int m) {
  csp::back(p, a, b, m);
}

// the cluster table of nnhclu and stacel: nc rows of w values from cs,
// columns named nm; the columns cx, cy and cx + 2, cy + 2 hold centres
SEXP ctab(const csp::pts &p, const double *cs, int nc, int w,
          const char *const *nm, int cx) {
  double *c[16];
  SEXP f = PROTECT(dframe(w, nc, nm, c));
  for (int i = 0; i < nc; ++i) {
    for (int j = 0; j < w; ++j) c[j][i] = cs[static_cast<R_xlen_t>(i) * w + j];
  }
  backcols(p, c[cx], c[cx + 1], nc);
  backcols(p, c[cx + 2], c[cx + 3], nc);
  UNPROTECT(1);
  return f;
}

}  // namespace

extern "C" {

// --------------------------------------------------------------- hot spots

SEXP C_clab_r_fzmode(SEXP x, SEXP y, SEXP r, SEXP lonlat) {
  clab::who = "fzmode";
  csp::pts p = csp::xy(x, y, lonlat);
  const double rad = in(r, "r", 0.0, std::numeric_limits<double>::max());
  SEXP iw = PROTECT(ints(2 * static_cast<R_xlen_t>(p.n) + 2));
  SEXP out = PROTECT(ints(p.n));
  kc(clab_fzmode(p.x, p.y, p.n, rad, INTEGER(iw), INTEGER(out)));
  UNPROTECT(3);
  return out;
}

SEXP C_clab_r_kmpart(SEXP x, SEXP y, SEXP k, SEXP init, SEXP cell, SEXP upd,
                     SEXP maxit, SEXP seed, SEXP lonlat) {
  clab::who = "kmpart";
  csp::pts p = csp::xy(x, y, lonlat);
  const int n = p.n;
  const int kk = clab::whole(k, "k", 1);
  if (kk > n) fail("`k` must be at most the number of points.");
  SEXP cen = PROTECT(real(2 * static_cast<R_xlen_t>(kk)));
  int in0 = 1;
  double par[2] = {0, 0};
  if (TYPEOF(init) == STRSXP) {
    static const char *const io[] = {"pp", "cells"};
    in0 = 1 + csp::pick(init, "init", io, 2);
  } else {
    // starting centres: a k x 2 matrix or data frame, read as points
    csp::pts c = csp::xy(init, R_NilValue, R_NilValue, p.ll ? csp::raw : csp::planar);
    if (c.n != kk) fail("`init` must give `k` starting centres.");
    std::memcpy(REAL(cen), c.x, kk * sizeof(double));
    std::memcpy(REAL(cen) + kk, c.y, kk * sizeof(double));
    UNPROTECT(1);
    csp::fwd(p, REAL(cen), REAL(cen) + kk, kk);
    in0 = 0;
  }
  if (in0 == 2) {
    SEXP cv = PROTECT(num(cell, "cell", 2, 0.0));
    par[0] = REAL_RO(cv)[0];
    par[1] = REAL_RO(cv)[1];
    UNPROTECT(1);
    if (!(par[0] > 0)) fail("`cell` must give a positive cell side.");
  }
  static const char *const uo[] = {"mean", "weber"};
  const int up = csp::pick(upd, "upd", uo, 2);
  const int mi = clab::whole(maxit, "maxit", 0);
  const unsigned int sd = in0 == 1 ? csp::seed(seed) : 0u;
  SEXP cl = PROTECT(ints(n));
  SEXP iw = PROTECT(ints(2 * static_cast<R_xlen_t>(n)));
  SEXP wk = PROTECT(real(2 * static_cast<R_xlen_t>(n)));
  SEXP o = PROTECT(real(4 + 3 * static_cast<R_xlen_t>(kk)));
  kc(clab_kmpart(p.x, p.y, n, kk, in0, par, sd, up, mi, REAL(cen),
                 INTEGER(cl), INTEGER(iw), REAL(wk), REAL(o)),
     in0 == 1 ? "fewer distinct points than `k`."
              : "fewer than `k` cells hold points far enough apart.");
  static const char *const cn[] = {"x", "y", "n", "ss", "sd"};
  double *c[5];
  SEXP ct = PROTECT(dframe(5, kk, cn, c));
  for (int j = 0; j < kk; ++j) {
    c[0][j] = REAL(cen)[j];
    c[1][j] = REAL(cen)[kk + j];
    for (int b = 0; b < 3; ++b) c[2 + b][j] = REAL(o)[4 + 3 * j + b];
  }
  backcols(p, c[0], c[1], kk);
  static const char *const fn[] = {"ss", "sd", "iter", "conv"};
  SEXP fit = PROTECT(clab::nvec(4, fn, REAL(o)));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, grpno(INTEGER(cl), n));
  SET_VECTOR_ELT(out, 1, ct);
  SET_VECTOR_ELT(out, 2, fit);
  static const char *const nm[] = {"cl", "cen", "fit"};
  UNPROTECT(9);
  return lst(out, nm);
}

SEXP C_clab_r_nnhclu(SEXP x, SEXP y, SEXP area, SEXP fixd, SEXP pv,
                     SEXP minpts, SEXP minhi, SEXP maxord, SEXP k, SEXP c,
                     SEXP lonlat) {
  clab::who = "nnhclu";
  csp::pts p = csp::xy(x, y, lonlat);
  const int n = p.n;
  if (n < 2) fail("at least 2 points are needed.");
  if (Rf_isNull(area) == Rf_isNull(fixd)) {
    fail("give one of `area` (for the random-expectation threshold) and "
         "`fixd` (a fixed threshold).");
  }
  const double ar = Rf_isNull(area) ? 0.0 : posn(area, "area");
  const double fd = Rf_isNull(fixd) ? 0.0 : posn(fixd, "fixd");
  const double z = Rf_qnorm5(clab::one(pv, "p", 0.0, 1.0), 0.0, 1.0, 1, 0);
  const int mp = clab::whole(minpts, "minpts", 2);
  const int mh = clab::whole(minhi, "minhi", 2);
  const int mo = clab::whole(maxord, "maxord", 1);
  if (mo > 20) fail("`maxord` must be at most 20.");
  const int kk = clab::whole(k, "k", 0);
  if (kk > 2) fail("`k` must be 0, 1 or 2.");
  const double cc = posn(c, "c");
  SEXP iw = PROTECT(ints(6 * static_cast<R_xlen_t>(n) + 2));
  SEXP wk = PROTECT(real(4 * static_cast<R_xlen_t>(n)));
  SEXP cl = PROTECT(Rf_allocMatrix(INTSXP, n, mo));
  SEXP cs = PROTECT(real(11 * static_cast<R_xlen_t>(n)));
  int nc = 0;
  kc(clab_nnhclu(p.x, p.y, n, ar, fd, z, mp, mh, mo, kk, cc, INTEGER(iw),
                 REAL(wk), INTEGER(cl), &nc, REAL(cs)));
  int *q = INTEGER(cl);
  for (R_xlen_t i = 0; i < static_cast<R_xlen_t>(n) * mo; ++i) {
    q[i] = q[i] >= 0 ? q[i] + 1 : NA_INTEGER;
  }
  static const char *const cn[] = {"order", "points", "units", "mx", "my",
                                   "wx", "wy", "a", "b", "angle", "thresh"};
  SEXP ct = PROTECT(ctab(p, REAL(cs), nc, 11, cn, 3));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, cl);
  SET_VECTOR_ELT(out, 1, ct);
  static const char *const nm[] = {"cl", "clusters"};
  UNPROTECT(7);
  return lst(out, nm);
}

SEXP C_clab_r_stacel(SEXP x, SEXP y, SEXP win, SEXP nx, SEXP ny, SEXP r,
                     SEXP top, SEXP minpts, SEXP k, SEXP c) {
  clab::who = "stacel";
  csp::pts p = csp::xy(x, y, R_NilValue, csp::planar);
  double w[4];
  csp::win(win, w);
  const int gx = clab::whole(nx, "nx", 2), gy = clab::whole(ny, "ny", 2);
  if (gx > 10000 || gy > 10000) fail("`nx` and `ny` must be at most 10000.");
  const double rad = posn(r, "r");
  const int tp = clab::whole(top, "top", 1);
  const int mp = clab::whole(minpts, "minpts", 2);
  const int kk = clab::whole(k, "k", 0);
  if (kk > 2) fail("`k` must be 0, 1 or 2.");
  const double cc = posn(c, "c");
  const R_xlen_t n = p.n, nn = static_cast<R_xlen_t>(gx) * gy;
  SEXP iw = PROTECT(ints(3 * n + 2 + nn + tp + std::max<R_xlen_t>(n, 2 * static_cast<R_xlen_t>(tp))));
  SEXP wk = PROTECT(real(nn + 2 * n));
  SEXP cl = PROTECT(ints(n));
  SEXP cs = PROTECT(real(8 * n));
  int nc = 0;
  kc(clab_stacel(p.x, p.y, p.n, w, gx, gy, rad, tp, mp, kk, cc, INTEGER(iw),
                 REAL(wk), INTEGER(cl), &nc, REAL(cs)));
  static const char *const cn[] = {"points", "mx", "my", "wx", "wy",
                                   "a", "b", "angle"};
  SEXP ct = PROTECT(ctab(p, REAL(cs), nc, 8, cn, 1));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, grpno(INTEGER(cl), n));
  SET_VECTOR_ELT(out, 1, ct);
  static const char *const nm[] = {"cl", "clusters"};
  UNPROTECT(7);
  return lst(out, nm);
}

}  // extern "C"

namespace {

// a matrix of finite numbers from 0, as a plain double vector
SEXP nummat(SEXP x, const char *nm) {
  if (!clab::isnum(x) || XLENGTH(x) < 1) {
    fail("`%s` must be finite numbers from 0 with no missing values.", nm);
  }
  SEXP v = PROTECT(clab::dcopy(PROTECT(clab::dbl(x))));
  const double *p = REAL_RO(v);
  for (R_xlen_t i = 0; i < XLENGTH(v); ++i) {
    if (!(p[i] >= 0) || !std::isfinite(p[i])) {
      fail("`%s` must be finite numbers from 0 with no missing values.", nm);
    }
  }
  UNPROTECT(2);
  return v;
}

// times: finite numbers, one per point
SEXP times(SEXP t, int n) { return num(t, "t", n); }

// a grid of nx by ny cells over the rectangle win: corner and cell size
struct grid {
  int nx, ny;
  double x0, y0, dx, dy;
};

grid mkgrid(SEXP win, int nx, int ny) {
  double w[4];
  csp::win(win, w);
  grid g{nx, ny, w[0], w[2], 0, 0};
  if (static_cast<double>(g.nx) * g.ny > 1e9) {
    fail("`nx` times `ny` must be at most 1e9 cells.");
  }
  g.dx = (w[1] - w[0]) / g.nx;
  g.dy = (w[3] - w[2]) / g.ny;
  return g;
}

const char *const dtype[] = {"linear", "negexp", "normal", "lognormal",
                             "trunc", "table"};
const int dlen[] = {2, 2, 3, 3, 3, 1};

// a decay function: its type, parameters and table (type "table")
struct decay {
  int type;
  SEXP par, tab;
};

decay mkdecay(SEXP type, SEXP par, SEXP tab) {
  decay d{csp::pick(type, "type", dtype, 6), R_NilValue, R_NilValue};
  d.par = PROTECT(num(par, "par", dlen[d.type]));
  if (d.type == 5) {
    d.tab = num(tab, "tab", -1, 0.0);
    if (XLENGTH(d.tab) < 2) fail("`tab` must hold at least 2 values.");
  }
  PROTECT(d.tab);
  return d;
}

}  // namespace

extern "C" {

// -------------------------------------------------------------- space-time

SEXP C_clab_r_mantst(SEXP x, SEXP y, SEXP t, SEXP cs, SEXP ct, SEXP nsim,
                     SEXP seed, SEXP lonlat) {
  clab::who = "mantst";
  SEXP uv[3] = {x, y, t};
  static const char *const un[] = {"x", "y", "t"};
  const int nb = csp::byname(uv, un, 3);
  x = uv[0];
  y = uv[1];
  t = uv[2];
  csp::pts p = csp::xy(x, y, lonlat);
  SEXP tv = PROTECT(times(t, p.n));
  if (p.n < 3 || p.n > 46340) fail("the number of events must be from 3 to 46340.");
  const double a = in(cs, "cs", 0.0, std::numeric_limits<double>::max());
  const double b = in(ct, "ct", 0.0, std::numeric_limits<double>::max());
  const int ns = clab::whole(nsim, "nsim", 0);
  const unsigned int sd = ns > 0 ? csp::seed(seed) : 0u;
  SEXP iw = PROTECT(ints(p.n));
  double o[5];
  kc(clab_mantst(p.x, p.y, REAL_RO(tv), p.n, a, b, ns, sd, INTEGER(iw), o),
     "the distances or the time intervals are all equal: r is undefined.");
  static const char *const nm[] = {"r", "pairs", "pup", "plo", "p2"};
  UNPROTECT(3 + nb);
  return clab::nvec(5, nm, o);
}

SEXP C_clab_r_stmavg(SEXP x, SEXP y, SEXP t, SEXP span, SEXP lonlat) {
  clab::who = "stmavg";
  SEXP uv[3] = {x, y, t};
  static const char *const un[] = {"x", "y", "t"};
  const int nb = csp::byname(uv, un, 3);
  x = uv[0];
  y = uv[1];
  t = uv[2];
  csp::pts p = csp::xy(x, y, lonlat);
  SEXP tv = PROTECT(times(t, p.n));
  const int sp = clab::whole(span, "span", 1);
  if (sp % 2 == 0) fail("`span` must be odd.");
  SEXP iw = PROTECT(ints(p.n));
  SEXP o = PROTECT(real(4 * static_cast<R_xlen_t>(p.n)));
  kc(clab_stmavg(p.x, p.y, REAL_RO(tv), p.n, sp, INTEGER(iw), REAL(o)));
  static const char *const nm[] = {"mx", "my", "n"};
  double *c[3];
  SEXP f = PROTECT(dframe(3, p.n, nm, c));
  for (int q = 0; q < p.n; ++q) {
    const double *r = REAL(o) + 4 * static_cast<R_xlen_t>(q);
    const int i = static_cast<int>(r[0]);
    c[0][i] = r[1];
    c[1][i] = r[2];
    c[2][i] = r[3];
  }
  csp::back(p, c[0], c[1], p.n);
  UNPROTECT(5 + nb);
  return f;
}

// --------------------------------------------------------- journey to crime

SEXP C_clab_r_jtdecay(SEXP d, SEXP type, SEXP par, SEXP tab) {
  clab::who = "jtdecay";
  SEXP dv = PROTECT(num(d, "d", -1, 0.0));
  const int n = clab::ilen(static_cast<double>(XLENGTH(dv)), "d");
  decay f = mkdecay(type, par, tab);
  SEXP out = PROTECT(real(n));
  kc(clab_jtdecay(REAL_RO(dv), n, f.type, REAL_RO(f.par),
                  f.type == 5 ? REAL_RO(f.tab) : nullptr,
                  f.type == 5 ? static_cast<int>(XLENGTH(f.tab)) : 0, REAL(out)),
     "`par` is outside the domain of this decay type.");
  UNPROTECT(4);
  return out;
}

SEXP C_clab_r_jtsurf(SEXP x, SEXP y, SEXP win, SEXP nx, SEXP ny, SEXP type,
                     SEXP par, SEXP tab, SEXP metric, SEXP mode) {
  clab::who = "jtsurf";
  csp::pts p = csp::xy(x, y, R_NilValue, csp::planar);
  grid g = mkgrid(win, clab::whole(nx, "nx", 1), clab::whole(ny, "ny", 1));
  decay f = mkdecay(type, par, tab);
  static const char *const mo[] = {"euclid", "manhattan"};
  const int me = csp::pick(metric, "metric", mo, 2);
  static const char *const so[] = {"sum", "lik"};
  const int md = csp::pick(mode, "mode", so, 2);
  SEXP raw = PROTECT(Rf_allocMatrix(REALSXP, g.nx, g.ny));
  SEXP prob = PROTECT(Rf_allocMatrix(REALSXP, g.nx, g.ny));
  double pk[3];
  kc(clab_jtsurf(p.x, p.y, p.n, g.nx, g.ny, g.x0, g.y0, g.dx, g.dy, me, md,
                 f.type, REAL_RO(f.par), f.type == 5 ? REAL_RO(f.tab) : nullptr,
                 f.type == 5 ? static_cast<int>(XLENGTH(f.tab)) : 0, REAL(raw),
                 REAL(prob), pk),
     "the surface is 0 everywhere: no cell is within reach of the decay "
     "function, or `par` is outside its domain.");
  pk[0] += 1;
  static const char *const pn[] = {"cell", "x", "y"};
  SEXP peak = PROTECT(clab::nvec(3, pn, pk));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, raw);
  SET_VECTOR_ELT(out, 1, prob);
  SET_VECTOR_ELT(out, 2, peak);
  static const char *const nm[] = {"raw", "prob", "peak"};
  UNPROTECT(7);
  return lst(out, nm);
}

SEXP C_clab_r_jtcal(SEXP d, SEXP t, SEXP h, SEXP w, SEXP kern, SEXP refl) {
  clab::who = "jtcal";
  SEXP uv[2] = {d, w};
  static const char *const un[] = {"d", "w"};
  const int nb = csp::byname(uv, un, 2);
  d = uv[0];
  w = uv[1];
  SEXP dv = PROTECT(num(d, "d", -1, 0.0));
  const int m = clab::ilen(static_cast<double>(XLENGTH(dv)), "d");
  const bool rf = csp::flag(refl, "refl");
  SEXP tv = PROTECT(num(t, "t", -1, rf ? 0.0 : clab::ninf));
  const int k = clab::ilen(static_cast<double>(XLENGTH(tv)), "t");
  const double bw = posn(h, "h");
  SEXP keep = R_NilValue;
  const double *wp = wts(w, m, &keep);
  PROTECT(keep);
  static const char *const ko[] = {"normal", "quartic", "triangular",
                                   "uniform", "laplace"};
  const int kk = csp::pick(kern, "kern", ko, 5);
  SEXP out = PROTECT(real(k));
  kc(clab_jtcal(REAL_RO(dv), wp, m, kk, bw, rf ? 1 : 0, REAL_RO(tv), k,
                REAL(out)));
  UNPROTECT(4 + nb);
  return out;
}

SEXP C_clab_r_jtacc(SEXP prob, SEXP win, SEXP base) {
  clab::who = "jtacc";
  SEXP d = Rf_getAttrib(prob, R_DimSymbol);
  if (TYPEOF(d) != INTSXP || XLENGTH(d) != 2) {
    fail("`prob` must be a matrix of cell values (nx by ny), as jtsurf() "
         "returns.");
  }
  SEXP pv = PROTECT(nummat(prob, "prob"));
  grid g = mkgrid(win, INTEGER_RO(d)[0], INTEGER_RO(d)[1]);
  SEXP bv = PROTECT(num(base, "base", 2));
  const double hx = REAL_RO(bv)[0], hy = REAL_RO(bv)[1];
  double w[4] = {g.x0, g.x0 + g.dx * g.nx, g.y0, g.y0 + g.dy * g.ny};
  if (!csp::inside(w, hx, hy)) fail("`base` must lie inside `win`.");
  double o[6];
  kc(clab_jtacc(REAL_RO(pv), g.nx, g.ny, g.x0, g.y0, g.dx, g.dy, hx, hy, o));
  o[4] += 1;
  o[5] += 1;
  static const char *const nm[] = {"p", "hitbest", "hitworst", "err",
                                   "peak", "cell"};
  UNPROTECT(2);
  return clab::nvec(6, nm, o);
}

SEXP C_clab_r_bjtc(SEXP pj, SEXP od, SEXP dest) {
  clab::who = "bjtc";
  pj = PROTECT(csp::gvec(pj, "pj"));
  SEXP pv = PROTECT(num(pj, "pj", -1, 0.0));
  const int no = clab::ilen(static_cast<double>(XLENGTH(pv)), "pj");
  SEXP d = Rf_getAttrib(od, R_DimSymbol);
  if (TYPEOF(d) != INTSXP || XLENGTH(d) != 2 || INTEGER_RO(d)[0] != no) {
    fail("`od` must be a matrix with one row per origin zone (the length "
         "of `pj`).");
  }
  const int nd = INTEGER_RO(d)[1];
  SEXP ov = PROTECT(nummat(od, "od"));
  // the kernel reads od row by row: od[o nd + j]
  SEXP tr = PROTECT(real(XLENGTH(ov)));
  for (int o = 0; o < no; ++o) {
    for (int j = 0; j < nd; ++j) {
      REAL(tr)[static_cast<R_xlen_t>(o) * nd + j] =
          REAL_RO(ov)[o + static_cast<R_xlen_t>(no) * j];
    }
  }
  SEXP dv = PROTECT(num(dest, "dest", -1, 1.0));
  const int m = clab::ilen(static_cast<double>(XLENGTH(dv)), "dest");
  SEXP di = PROTECT(ints(m));
  for (int k = 0; k < m; ++k) {
    const double v = REAL_RO(dv)[k];
    if (v != std::floor(v) || v > nd) {
      fail("`dest` must be whole numbers from 1 to %d, the columns of `od`.", nd);
    }
    INTEGER(di)[k] = static_cast<int>(v) - 1;
  }
  SEXP o = PROTECT(real(5 * static_cast<R_xlen_t>(no)));
  kc(clab_bjtc(REAL_RO(pv), no, REAL(tr), nd, INTEGER(di), m, REAL(o)),
     "a block sums to 0: no journey-to-crime value, no solved offence, or "
     "no solved offence in the series' zones.");
  static const char *const nm[] = {"pj", "po", "pod", "prod", "risk"};
  double *c[5];
  SEXP f = PROTECT(dframe(5, no, nm, c));
  for (int b = 0; b < 5; ++b) {
    std::memcpy(c[b], REAL(o) + static_cast<R_xlen_t>(b) * no, no * sizeof(double));
  }
  UNPROTECT(8);
  return f;
}

}  // extern "C"
