// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of head-bang smoothing, the correlated random walk and
// risk-adjusted nearest-neighbour clustering (kernels in src/sp2.cpp).
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

SEXP real(R_xlen_t n) { return Rf_allocVector(REALSXP, n); }
SEXP ints(R_xlen_t n) { return Rf_allocVector(INTSXP, n); }

void kc(int rc, const char *msg) {
  if (rc != CLAB_OK) fail("%s", msg);
}

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

SEXP C_clab_r_hdbang(SEXP x, SEXP y, SEXP v, SEXP w, SEXP nn, SEXP ang,
                     SEXP ntrip, SEXP iter, SEXP lonlat) {
  clab::who = "hdbang";
  SEXP uv[4] = {x, y, v, w};
  static const char *const un[] = {"x", "y", "v", "w"};
  const int nb = csp::byname(uv, un, 4);
  csp::pts p = csp::xy(uv[0], uv[1], lonlat);
  const int n = p.n;
  if (n < 3) fail("at least 3 points are needed.");
  SEXP vv = PROTECT(num(uv[2], "v", n));
  SEXP wv = PROTECT(Rf_isNull(uv[3]) ? R_NilValue : num(uv[3], "w", n, 0.0, true));
  const int k = clab::whole(nn, "nn", 2);
  if (k > n - 1) fail("`nn` must be below the number of points.");
  const double a = clab::val1(ang);
  if (!(a >= 90 && a <= 180)) fail("`ang` must be one number from 90 to 180.");
  const int nt = clab::whole(ntrip, "ntrip", 0);
  const int it = clab::whole(iter, "iter", 0);
  const R_xlen_t T = static_cast<R_xlen_t>(k) * (k - 1) / 2;
  SEXP iw = PROTECT(ints(2 * static_cast<R_xlen_t>(n) + 2 + static_cast<R_xlen_t>(n) * k + 3 * T));
  SEXP wk = PROTECT(real(static_cast<R_xlen_t>(n) * k + n + 5 * T));
  SEXP out = PROTECT(real(n));
  kc(clab_hdbang(p.x, p.y, REAL(vv), Rf_isNull(wv) ? nullptr : REAL(wv), n, k, a, nt, it,
                 INTEGER(iw), REAL(wk), REAL(out)),
     "the kernel refused the input.");
  UNPROTECT(6 + nb);
  return out;
}

SEXP C_clab_r_crwalk(SEXP x, SEXP y, SEXP h, SEXP nsim, SEXP seed,
                     SEXP lonlat) {
  clab::who = "crwalk";
  csp::pts p = csp::xy(x, y, lonlat);
  const int n = p.n;
  if (n < 3) fail("at least 3 locations are needed.");
  const int hh = clab::whole(h, "h", 1);
  const int ns = clab::whole(nsim, "nsim", 1);
  const unsigned int sd = csp::seed(seed);
  SEXP st = PROTECT(real(n - 1));
  SEXP tu = PROTECT(real(n - 2));
  double f[5];
  kc(clab_crwfit(p.x, p.y, n, REAL(st), REAL(tu), f),
     "two successive locations coincide: a repeated location has no heading; "
     "drop repeats first.");
  static const char *const mn[] = {"k", "obs", "exp", "x", "y", "sdx", "sdy"};
  double *c[7];
  SEXP m = PROTECT(dframe(7, hh, mn, c));
  kc(clab_crwmsd(p.x, p.y, n, f[0], f[1], f[2], hh, c[1], c[2]), "the kernel refused the input.");
  SEXP px = PROTECT(real(static_cast<R_xlen_t>(hh) * ns));
  SEXP py = PROTECT(real(static_cast<R_xlen_t>(hh) * ns));
  kc(clab_crwsim(REAL(st), n - 1, REAL(tu), n - 2, p.x[n - 1], p.y[n - 1], f[4], hh, ns, sd,
                 REAL(px), REAL(py)),
     "the kernel refused the input.");
  for (int k = 0; k < hh; ++k) {
    double sx = 0, sy = 0, qx = 0, qy = 0;
    const double *ax = REAL(px) + static_cast<R_xlen_t>(k) * ns;
    const double *ay = REAL(py) + static_cast<R_xlen_t>(k) * ns;
    for (int s = 0; s < ns; ++s) {
      sx += ax[s];
      sy += ay[s];
    }
    sx /= ns;
    sy /= ns;
    for (int s = 0; s < ns; ++s) {
      qx += (ax[s] - sx) * (ax[s] - sx);
      qy += (ay[s] - sy) * (ay[s] - sy);
    }
    c[0][k] = k + 1;
    c[3][k] = sx;
    c[4][k] = sy;
    c[5][k] = ns > 1 ? std::sqrt(qx / (ns - 1)) : std::numeric_limits<double>::quiet_NaN();
    c[6][k] = ns > 1 ? std::sqrt(qy / (ns - 1)) : std::numeric_limits<double>::quiet_NaN();
  }
  // the predicted centres back to degrees (their spreads stay in metres)
  csp::back(p, c[3], c[4], hh);
  const double fv[6] = {static_cast<double>(n - 1), f[0], f[1], f[2], f[3], f[4]};
  static const char *const fn[] = {"steps", "ml", "ml2", "c", "s", "head"};
  SEXP fit = PROTECT(clab::nvec(6, fn, fv));
  SEXP sim = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP mx = Rf_allocMatrix(REALSXP, ns, hh);
  SET_VECTOR_ELT(sim, 0, mx);
  std::memcpy(REAL(mx), REAL(px), sizeof(double) * static_cast<std::size_t>(hh) * static_cast<std::size_t>(ns));
  SEXP my = Rf_allocMatrix(REALSXP, ns, hh);
  SET_VECTOR_ELT(sim, 1, my);
  std::memcpy(REAL(my), REAL(py), sizeof(double) * static_cast<std::size_t>(hh) * static_cast<std::size_t>(ns));
  static const char *const sn[] = {"x", "y"};
  clab::nlist(sim, sn);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, fit);
  SET_VECTOR_ELT(out, 1, m);
  SET_VECTOR_ELT(out, 2, sim);
  static const char *const on[] = {"fit", "msd", "sim"};
  UNPROTECT(9);
  return clab::nlist(out, on);
}

SEXP C_clab_r_rnnhcl(SEXP x, SEXP y, SEXP base, SEXP bw, SEXP h, SEXP pv,
                     SEXP minpts, SEXP k, SEXP c, SEXP lonlat) {
  clab::who = "rnnhcl";
  csp::pts p = csp::xy(x, y, lonlat);
  const int n = p.n;
  if (n < 2) fail("at least 2 points are needed.");
  // the baseline in the same coordinates and projection as the points
  csp::pts b = csp::xy(base, R_NilValue, R_NilValue, p.ll ? csp::raw : csp::planar);
  SEXP bxy = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(bxy, 0, clab::dcopy(VECTOR_ELT(b.v, 0)));
  SET_VECTOR_ELT(bxy, 1, clab::dcopy(VECTOR_ELT(b.v, 1)));
  double *bx = REAL(VECTOR_ELT(bxy, 0)), *by = REAL(VECTOR_ELT(bxy, 1));
  csp::fwd(p, bx, by, b.n);
  SEXP bv = PROTECT(Rf_isNull(bw) ? R_NilValue : num(bw, "bw", b.n, 0.0));
  SEXP ones = PROTECT(real(b.n));
  for (int q = 0; q < b.n; ++q) REAL(ones)[q] = Rf_isNull(bv) ? 1.0 : REAL(bv)[q];
  const double hb = clab::one(h, "h", 0.0, std::numeric_limits<double>::infinity());
  const double z = Rf_qnorm5(clab::one(pv, "p", 0.0, 1.0), 0.0, 1.0, 1, 0);
  const int mp = clab::whole(minpts, "minpts", 2);
  const int kk = clab::whole(k, "k", 0);
  if (kk > 2) fail("`k` must be 0, 1 or 2.");
  const double cc = clab::one(c, "c", 0.0, std::numeric_limits<double>::infinity());
  SEXP iw = PROTECT(ints(4 * static_cast<R_xlen_t>(n)));
  SEXP wk = PROTECT(real(3 * static_cast<R_xlen_t>(n)));
  SEXP thr = PROTECT(real(n));
  SEXP cl = PROTECT(ints(n));
  SEXP cs = PROTECT(real(8 * static_cast<R_xlen_t>(n)));
  int nc = 0;
  kc(clab_rnnhcl(p.x, p.y, n, bx, by, REAL(ones), b.n, hb, z, mp, kk, cc, INTEGER(iw),
                 REAL(wk), REAL(thr), INTEGER(cl), &nc, REAL(cs)),
     "`bw` must have a positive sum, and `p` must leave a positive threshold.");
  for (int i = 0; i < n; ++i) INTEGER(cl)[i] = INTEGER(cl)[i] >= 0 ? INTEGER(cl)[i] + 1 : NA_INTEGER;
  static const char *const cn[] = {"points", "mx", "my", "wx", "wy", "a", "b", "angle"};
  double *col[8];
  SEXP ct = PROTECT(dframe(8, nc, cn, col));
  for (int i = 0; i < nc; ++i)
    for (int j = 0; j < 8; ++j) col[j][i] = REAL(cs)[8 * static_cast<R_xlen_t>(i) + j];
  csp::back(p, col[1], col[2], nc);
  csp::back(p, col[3], col[4], nc);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, cl);
  SET_VECTOR_ELT(out, 1, ct);
  SET_VECTOR_ELT(out, 2, thr);
  static const char *const on[] = {"cl", "clusters", "thr"};
  UNPROTECT(12);
  return clab::nlist(out, on);
}

}  // extern "C"
