// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the geostatistics (src/krg.cpp). Each reads its
// points through csp::xy() (spin.h), checks the other arguments, calls the
// kernels and builds the R result. Nothing here calls back into the R
// evaluator.
#include <cmath>
#include <cstdio>
#include <cstring>

#include "bkh.h"
#include "crimlab.h"
#include "kutil.h"
#include "rchk.h"
#include "spin.h"

using namespace bkh;
using clab::num;

namespace {

const char *const vmodels[] = {"spherical", "exponential", "gaussian"};

// the covariance parameters of the block methods
void covpar(SEXP nugget, SEXP sill, SEXP range, double *c) {
  SEXP a = PROTECT(num(nugget, "nugget", 1, 0.0));
  SEXP b = PROTECT(num(sill, "sill", 1));
  c[0] = REAL(a)[0];
  c[1] = REAL(b)[0];
  c[2] = clab::one(range, "range", 0.0, HUGE_VAL);
  if (c[1] < c[0]) fail("`sill` must be at least `nugget`.");
  UNPROTECT(2);
}

// the coordinates and EPSG code of a spobj() object
void sobj(SEXP s, const char *an, SEXP *xy, int *n, int *crs) {
  if (!Rf_inherits(s, "spobj")) fail("`%s` must be the object spobj() returns.", an);
  *xy = elt(s, "coords");
  if (!clab::isnum(*xy) || !Rf_isMatrix(*xy) || Rf_ncols(*xy) != 2) fail("`%s$coords` must be a matrix of two columns.", an);
  *n = Rf_nrows(*xy);
  *crs = static_cast<int>(clab::val1(elt(s, "crs")));
}

// the spatial linear models: sglmx (kind 0) and smixm (kind 1)
SEXP slm(SEXP v, SEXP z, SEXP x, SEXP y, int kind) {
  csp::pts p = csp::xy(x, y, R_NilValue, csp::planar);
  const int n = p.n;
  SEXP vv = PROTECT(num(v, "v", n));
  int k = 0;
  SEXP zz = zmat(z, "z", n, true, &k);
  if (n <= k) fail("`v` must have more values than `z` has columns.");
  double *b = dal(k), *s = dal(k), o[4];
  const int rc = kind == 0 ? clab_sglmx(REAL(zz), n, k, REAL(vv), p.x, p.y, b, s, o)
                           : clab_smixm(REAL(zz), n, k, REAL(vv), p.x, p.y, b, s, o);
  kc(rc, "the columns of `z` are collinear.");
  SEXP cn = PROTECT(Rf_allocVector(STRSXP, k));
  SEXP dn = Rf_isMatrix(z) ? Rf_getAttrib(z, R_DimNamesSymbol) : R_NilValue;
  SEXP zc = Rf_isNull(dn) ? R_NilValue : VECTOR_ELT(dn, 1);
  for (int j = 0; j < k; ++j) {
    char buf[24];
    std::snprintf(buf, sizeof buf, "z%d", j + 1);
    SET_STRING_ELT(cn, j, Rf_isNull(zc) ? Rf_mkChar(Rf_isNull(z) ? "1" : buf) : STRING_ELT(zc, j));
  }
  SEXP bv = PROTECT(real(k)), sv = PROTECT(real(k));
  for (int j = 0; j < k; ++j) {
    REAL(bv)[j] = b[j];
    REAL(sv)[j] = s[j];
  }
  Rf_namesgets(bv, cn);
  Rf_namesgets(sv, cn);
  SEXP r = PROTECT(Rf_allocVector(VECSXP, kind == 0 ? 6 : 7));
  SET_VECTOR_ELT(r, 0, bv);
  SET_VECTOR_ELT(r, 1, sv);
  if (kind == 0) {
    for (int j = 0; j < 3; ++j) SET_VECTOR_ELT(r, 2 + j, Rf_ScalarReal(o[j]));
    SET_VECTOR_ELT(r, 5, Rf_ScalarInteger(n));
    static const char *const nm[] = {"coef", "se", "sigma2", "phi", "tau2", "n"};
    lst(r, nm, "sglmx");
  } else {
    for (int j = 0; j < 3; ++j) SET_VECTOR_ELT(r, 2 + j, Rf_ScalarReal(o[j]));
    SET_VECTOR_ELT(r, 5, Rf_ScalarLogical(o[3] != 0));
    SET_VECTOR_ELT(r, 6, Rf_ScalarInteger(n));
    static const char *const nm[] = {"coef", "se", "sigma2", "tau2", "phi", "conv", "n"};
    lst(r, nm, "smixm");
  }
  UNPROTECT(7);
  return r;
}

}  // namespace

extern "C" {


// ------------------------------------------------------------ geostatistics

SEXP C_clab_r_vgemp(SEXP x, SEXP y, SEXP v, SEXP nbins, SEXP cutoff) {
  clab::who = "vgemp";
  csp::pts p = csp::xy(x, y, R_NilValue, csp::planar);
  const int n = p.n;
  if (n < 2) fail("`x` must hold at least 2 points.");
  SEXP vv = PROTECT(num(v, "v", n));
  const int nb = clab::whole(nbins, "nbins", 1);
  if (nb > kutil::maxcell) fail("`nbins` must be at most %.0f.", kutil::maxcell);
  double cut = NAN;
  if (!Rf_isNull(cutoff)) cut = clab::one(cutoff, "cutoff", 0.0, HUGE_VAL);
  double *d = dal(nb), *g = dal(nb), *np = dal(nb);
  int nk = 0;
  kc(clab_vgemp(p.x, p.y, n, REAL(vv), nb, cut, d, g, np, &nk), "the points all coincide.");
  SEXP c = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(c, 0, real(nk));
  SET_VECTOR_ELT(c, 1, real(nk));
  SET_VECTOR_ELT(c, 2, ints(nk));
  for (int k = 0; k < nk; ++k) {
    REAL(VECTOR_ELT(c, 0))[k] = d[k];
    REAL(VECTOR_ELT(c, 1))[k] = g[k];
    INTEGER(VECTOR_ELT(c, 2))[k] = static_cast<int>(np[k]);
  }
  static const char *const nm[] = {"dist", "gamma", "np"};
  SEXP f = PROTECT(clab::frame(c, nm));
  UNPROTECT(4);
  return f;
}

SEXP C_clab_r_vgfit(SEXP x, SEXP y, SEXP v, SEXP model) {
  clab::who = "vgfit";
  csp::pts p = csp::xy(x, y, R_NilValue, csp::planar);
  const int n = p.n;
  if (n < 3) fail("`x` must hold at least 3 points.");
  SEXP vv = PROTECT(num(v, "v", n));
  const int md = csp::pick(model, "model", vmodels, 3);
  double o[6];
  kc(clab_vgfit(p.x, p.y, n, REAL(vv), md, o), "`v` must vary, and the points must not all coincide.");
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 7));
  SET_VECTOR_ELT(r, 0, Rf_mkString(vmodels[md]));
  for (int k = 0; k < 4; ++k) SET_VECTOR_ELT(r, k + 1, Rf_ScalarReal(o[k]));
  SET_VECTOR_ELT(r, 5, Rf_ScalarLogical(o[4] != 0));
  SET_VECTOR_ELT(r, 6, Rf_ScalarInteger(static_cast<int>(o[5])));
  static const char *const nm[] = {"model", "nugget", "psill", "range", "loglik", "conv", "nfn"};
  lst(r, nm, "vgfit");
  UNPROTECT(3);
  return r;
}

SEXP C_clab_r_okrge(SEXP x, SEXP y, SEXP v, SEXP at, SEXP vgm) {
  clab::who = "okrge";
  csp::pts p = csp::xy(x, y, R_NilValue, csp::planar);
  const int n = p.n;
  SEXP vv = PROTECT(num(v, "v", n));
  csp::pts a = csp::xy(at, R_NilValue, R_NilValue, csp::planar);
  int md = 1;
  double par[3];
  if (Rf_isNull(vgm)) {
    if (n < 2) fail("`x` must hold at least 2 points to fit a variogram.");
    kc(clab_vgwls(p.x, p.y, n, REAL(vv), 1, par), "no pair of points is close enough for a variogram.");
  } else {
    if (TYPEOF(vgm) != VECSXP) fail("`vgm` must be a list with `model`, `nugget`, `psill` and `range`.");
    md = csp::pick(elt(vgm, "model"), "vgm$model", vmodels, 3);
    SEXP a1 = PROTECT(num(elt(vgm, "nugget"), "vgm$nugget", 1, 0.0));
    SEXP a2 = PROTECT(num(elt(vgm, "psill"), "vgm$psill", 1, 0.0));
    SEXP a3 = PROTECT(num(elt(vgm, "range"), "vgm$range", 1, 0.0, true));
    par[0] = REAL(a1)[0];
    par[1] = REAL(a2)[0];
    par[2] = REAL(a3)[0];
    UNPROTECT(3);
  }
  SEXP c = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(c, 0, real(a.n));
  SET_VECTOR_ELT(c, 1, real(a.n));
  kc(clab_okrge(p.x, p.y, n, REAL(vv), a.x, a.y, a.n, md, par[0], par[1], par[2], REAL(VECTOR_ELT(c, 0)),
                REAL(VECTOR_ELT(c, 1))),
     "the kriging system is singular: two points coincide without a nugget, or the variogram is flat.");
  static const char *const nm[] = {"pred", "var"};
  SEXP f = PROTECT(clab::frame(c, nm));
  UNPROTECT(5);
  return f;
}

SEXP C_clab_r_bkrig(SEXP x, SEXP y, SEXP v, SEXP blocks, SEXP nq, SEXP nugget, SEXP sill, SEXP range) {
  clab::who = "bkrig";
  csp::pts p = csp::xy(x, y, R_NilValue, csp::planar);
  const int n = p.n;
  SEXP vv = PROTECT(num(v, "v", n));
  const int q = clab::whole(nq, "nq", 1);
  double c[3];
  covpar(nugget, sill, range, c);
  if (TYPEOF(blocks) != VECSXP || XLENGTH(blocks) < 1) fail("`blocks` must be a list of blocks.");
  const int nb = clab::ilen(static_cast<double>(XLENGTH(blocks)), "blocks");
  // each block: a 2 x 2 box (rows lo, hi) or a matrix of points
  int *off = ial(nb + 1);
  off[0] = 0;
  const int k = static_cast<int>(std::nearbyint(std::sqrt(static_cast<double>(q))));
  double tot = 0;  // the points of all blocks, counted before any is made
  for (int b = 0; b < nb; ++b) {
    SEXP e = VECTOR_ELT(blocks, b);
    if (!clab::isnum(e) || !Rf_isMatrix(e) || Rf_ncols(e) != 2) {
      fail("`blocks[[%d]]` must be a numeric matrix of two columns.", b + 1);
    }
    const int r = Rf_nrows(e);
    tot += r == 2 ? static_cast<double>(k) * k : r;
    if (tot > kutil::maxcell) fail("the points of `blocks`, `nq` for each box, must be at most %.0f.", kutil::maxcell);
    off[b + 1] = off[b] + (r == 2 ? k * k : r);
  }
  double *bx = dal(off[nb]), *by = dal(off[nb]);
  for (int b = 0; b < nb; ++b) {
    SEXP e = PROTECT(clab::dbl(VECTOR_ELT(blocks, b)));
    const int r = Rf_nrows(VECTOR_ELT(blocks, b));
    const double *m = REAL(e);
    for (int i = 0; i < 2 * r; ++i) {
      if (!std::isfinite(m[i])) fail("`blocks[[%d]]` must hold finite numbers.", b + 1);
    }
    if (r == 2) {
      const double lo[2] = {m[0], m[2]}, hi[2] = {m[1], m[3]};
      clab_bkpts(lo, hi, q, bx + off[b], by + off[b]);
    } else {
      for (int i = 0; i < r; ++i) {
        bx[off[b] + i] = m[i];
        by[off[b] + i] = m[r + i];
      }
    }
    UNPROTECT(1);
  }
  SEXP o = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(o, 0, real(nb));
  SET_VECTOR_ELT(o, 1, real(nb));
  kc(clab_bkrig(p.x, p.y, n, REAL(vv), bx, by, off, nb, c[0], c[1], c[2], REAL(VECTOR_ELT(o, 0)),
                REAL(VECTOR_ELT(o, 1))),
     "the kriging system is singular: two points coincide without a nugget.");
  static const char *const nm[] = {"est", "se"};
  SEXP f = PROTECT(clab::frame(o, nm));
  UNPROTECT(4);
  return f;
}

SEXP C_clab_r_krigcv(SEXP x, SEXP y, SEXP v, SEXP nugget, SEXP sill, SEXP range) {
  clab::who = "krigcv";
  csp::pts p = csp::xy(x, y, R_NilValue, csp::planar);
  const int n = p.n;
  if (n < 3) fail("`x` must hold at least 3 points.");
  SEXP vv = PROTECT(num(v, "v", n));
  double c[3];
  covpar(nugget, sill, range, c);
  SEXP res = PROTECT(real(n));
  kc(clab_krigcv(p.x, p.y, n, REAL(vv), c[0], c[1], c[2], REAL(res)),
     "a kriging system is singular: two points coincide without a nugget.");
  double s2 = 0, s1 = 0;
  for (int i = 0; i < n; ++i) {
    s2 += REAL(res)[i] * REAL(res)[i];
    s1 += std::fabs(REAL(res)[i]);
  }
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(r, 0, Rf_ScalarReal(s2 / n));
  SET_VECTOR_ELT(r, 1, Rf_ScalarReal(std::sqrt(s2 / n)));
  SET_VECTOR_ELT(r, 2, Rf_ScalarReal(s1 / n));
  SET_VECTOR_ELT(r, 3, res);
  static const char *const nm[] = {"mspe", "rmspe", "mae", "resid"};
  lst(r, nm, "krigcv");
  UNPROTECT(4);
  return r;
}

SEXP C_clab_r_trsurf(SEXP x, SEXP y, SEXP v, SEXP order) {
  clab::who = "trsurf";
  csp::pts p = csp::xy(x, y, R_NilValue, csp::planar);
  const int n = p.n;
  SEXP vv = PROTECT(num(v, "v", n));
  const int o = clab::whole(order, "order", 0);
  if (o > 3) fail("`order` must be 0, 1, 2 or 3.");
  static const int np[] = {1, 3, 6, 10};
  const int k = np[o];
  if (n < k) fail("order %d needs at least %d points.", o, k);
  static const char *const cn[] = {"1", "x", "y", "x2", "y2", "xy", "x3", "y3", "x2y", "xy2"};
  double *b = dal(k), *s = dal(k), r2 = 0;
  kc(clab_trsurf(p.x, p.y, n, REAL(vv), o, b, s, &r2), "the monomials are collinear at these points.");
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 5));
  SET_VECTOR_ELT(r, 0, clab::nvec(k, cn, b));
  SET_VECTOR_ELT(r, 1, clab::nvec(k, cn, s));
  SET_VECTOR_ELT(r, 2, Rf_ScalarReal(r2));
  SET_VECTOR_ELT(r, 3, Rf_ScalarInteger(o));
  SET_VECTOR_ELT(r, 4, Rf_ScalarInteger(n));
  static const char *const nm[] = {"coef", "se", "r2", "order", "n"};
  lst(r, nm, "trsurf");
  UNPROTECT(3);
  return r;
}


SEXP C_clab_r_sglmx(SEXP v, SEXP z, SEXP x, SEXP y) {
  clab::who = "sglmx";
  return slm(v, z, x, y, 0);
}

SEXP C_clab_r_smixm(SEXP v, SEXP z, SEXP x, SEXP y) {
  clab::who = "smixm";
  return slm(v, z, x, y, 1);
}


SEXP C_clab_r_lfspec(SEXP s) {
  clab::who = "lfspec";
  SEXP xy;
  int n = 0, crs = 0;
  sobj(s, "x", &xy, &n, &crs);
  if (crs != 4326 && crs != 3857) fail("only EPSG:4326 and EPSG:3857 can be shown on a web map.");
  SEXP pts = PROTECT(Rf_allocMatrix(REALSXP, n, 2));
  double *px = REAL(pts), *py = REAL(pts) + n;
  SEXP d = PROTECT(clab::dbl(xy));
  if (crs == 3857) {
    kc(clab_wmerc(REAL(d), REAL(d) + n, n, 1, px, py), "the coordinates must be finite.");
  } else {
    std::memcpy(px, REAL(d), sizeof(double) * 2 * static_cast<std::size_t>(n));
  }
  static const char *const cn[] = {"lon", "lat"};
  SEXP dn = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP cv = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(cv, 0, Rf_mkChar(cn[0]));
  SET_STRING_ELT(cv, 1, Rf_mkChar(cn[1]));
  SET_VECTOR_ELT(dn, 1, cv);
  Rf_setAttrib(pts, R_DimNamesSymbol, dn);
  const double ctr[2] = {clab::nmean(px, n), clab::nmean(py, n)};
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 5));
  SET_VECTOR_ELT(r, 0, Rf_mkString("OpenStreetMap"));
  SET_VECTOR_ELT(r, 1, Rf_ScalarLogical(1));
  SET_VECTOR_ELT(r, 2, clab::nvec(2, cn, ctr));
  SET_VECTOR_ELT(r, 3, Rf_ScalarInteger(8));
  SET_VECTOR_ELT(r, 4, pts);
  static const char *const nm[] = {"provider", "interactive", "center", "zoom", "points"};
  lst(r, nm, nullptr);
  UNPROTECT(5);
  return r;
}

SEXP C_clab_r_spsamp(SEXP s, SEXP size, SEXP dens, SEXP sd) {
  clab::who = "spsamp";
  SEXP xy;
  int n = 0, crs = 0;
  sobj(s, "x", &xy, &n, &crs);
  const int k = clab::whole(size, "size", 0);
  if (k > n) fail("`size` must be at most the number of points, %d.", n);
  const bool bd = csp::flag(dens, "dens");
  const unsigned int seed = csp::seed(sd);
  SEXP d = PROTECT(clab::dbl(xy));
  SEXP pr = PROTECT(real(n));
  if (bd) {
    if (n < 2) fail("`x` must hold at least 2 points to sample by density.");
    kc(clab_spdens(REAL(d), REAL(d) + n, n, crs == 4326, REAL(pr)), "the points all coincide.");
  } else {
    for (int i = 0; i < n; ++i) REAL(pr)[i] = 1.0 / n;
  }
  SEXP idx = PROTECT(ints(k));
  kc(clab_spsamp(bd ? REAL(pr) : nullptr, n, k, seed, INTEGER(idx)), "the kernel refused the input.");
  Rf_setAttrib(idx, PROTECT(Rf_install("prob")), pr);
  UNPROTECT(4);
  return idx;
}

}  // extern "C"
