// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the spatial crime methods: near repeats (src/knx.cpp)
// and the rest that follow. Each reads its points through csp::xy()
// (spin.h), checks the other arguments, calls the kernels and builds the R
// result. Nothing here calls back into the R evaluator.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <vector>

#include "crimlab.h"
#include "rchk.h"
#include "spin.h"

// R's distribution functions under their library names (Rmath.h declares
// them so only in newer versions)
extern "C" double Rf_ppois(double, double, int, int);
extern "C" double Rf_pnorm5(double, double, double, int, int);
extern "C" double Rf_pchisq(double, double, int, int);
extern "C" double Rf_qchisq(double, double, int, int);

using clab::fail;
using clab::num;

namespace {

SEXP real(R_xlen_t n) { return Rf_allocVector(REALSXP, n); }

void kc(int rc, const char *msg) {
  if (rc != CLAB_OK) fail("%s", msg);
}

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

// band edges: positive, finite and increasing (returned protected)
SEXP edges(SEXP e, const char *nm) {
  SEXP v = PROTECT(num(e, nm, -1, 0.0, true));
  const R_xlen_t k = XLENGTH(v);
  for (R_xlen_t q = 1; q < k; ++q) {
    if (!(REAL(v)[q] > REAL(v)[q - 1])) fail("`%s` must increase.", nm);
  }
  return v;
}

const char *const xyn[] = {"x", "y"};

// a classed list
SEXP klass(SEXP x, const char *cl) {
  Rf_classgets(x, PROTECT(Rf_mkString(cl)));
  UNPROTECT(1);
  return x;
}

// the normal-reference bandwidth of one axis, MASS::bandwidth.nrd() / 4:
// 1.06 min(sd, IQR / 1.34) n^(-1/5), IQR from type-7 quantiles
double nrd(const double *v, int n) {
  std::vector<double> u(v, v + n);
  std::sort(u.begin(), u.end());
  auto q = [&](double pr) {
    const double h = (n - 1) * pr;
    const int lo = static_cast<int>(std::floor(h));
    return lo + 1 < n ? u[lo] + (h - lo) * (u[lo + 1] - u[lo]) : u[lo];
  };
  double m = 0, ss = 0;
  for (int i = 0; i < n; ++i) m += v[i];
  m /= n;
  for (int i = 0; i < n; ++i) ss += (v[i] - m) * (v[i] - m);
  const double sd = std::sqrt(ss / (n - 1));
  return 1.06 * std::min(sd, (q(0.75) - q(0.25)) / 1.34) * std::pow(n, -0.2);
}

}  // namespace

extern "C" {

// ---------------------------------------------------------- near repeats

SEXP C_clab_r_knxnr(SEXP x, SEXP y, SEXP t, SEXP ds, SEXP dt, SEXP nsim, SEXP seed, SEXP lonlat) {
  clab::who = "knxnr";
  csp::pts p = csp::xy(x, y, lonlat);
  const int n = p.n;
  if (n < 3) fail("at least 3 incidents are needed.");
  SEXP tv = PROTECT(num(t, "t", n));
  SEXP sv = edges(ds, "ds");
  SEXP bv = edges(dt, "dt");
  const int ks = static_cast<int>(XLENGTH(sv)), kt = static_cast<int>(XLENGTH(bv));
  const int ns = clab::whole(nsim, "nsim", 0);
  const unsigned int sd = ns > 0 ? csp::seed(seed) : 0u;
  const int m = ks * kt;
  SEXP o = PROTECT(real(3 * static_cast<R_xlen_t>(m) + 1));
  kc(clab_knxnr(p.x, p.y, REAL(tv), n, REAL(sv), ks, REAL(bv), kt, ns, sd, REAL(o)),
     "the kernel refused the input.");
  static const char *const tn[] = {"dlo", "dhi", "tlo", "thi", "obs", "exp", "ratio", "p", "ppois"};
  double *c[9];
  SEXP tab = PROTECT(dframe(9, m, tn, c));
  const double *ob = REAL(o), *ex = REAL(o) + m, *pp = REAL(o) + 2 * m;
  for (int b = 0; b < kt; ++b) {
    for (int s = 0; s < ks; ++s) {
      const int q = s + b * ks;
      c[0][q] = s ? REAL(sv)[s - 1] : 0.0;
      c[1][q] = REAL(sv)[s];
      c[2][q] = b ? REAL(bv)[b - 1] : 0.0;
      c[3][q] = REAL(bv)[b];
      c[4][q] = ob[q];
      c[5][q] = ex[q];
      c[6][q] = ex[q] > 0 ? ob[q] / ex[q] : std::numeric_limits<double>::quiet_NaN();
      c[7][q] = pp[q];
      // P(X >= obs) for X Poisson with the expected count as its mean
      c[8][q] = ob[q] > 0 ? Rf_ppois(ob[q] - 1, ex[q], 0, 0) : 1.0;
    }
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out, 0, tab);
  SET_VECTOR_ELT(out, 1, Rf_ScalarInteger(n));
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(REAL(o)[3 * m]));
  SET_VECTOR_ELT(out, 3, Rf_ScalarInteger(ns));
  static const char *const on[] = {"tab", "n", "pairs", "nsim"};
  clab::nlist(out, on);
  UNPROTECT(7);
  return klass(out, "knxnr");
}

// --------------------------------------------------------------- hot spots

SEXP C_clab_r_kdsurf(SEXP x, SEXP y, SEXP w, SEXP h, SEXP win, SEXP nx, SEXP ny, SEXP edge, SEXP at,
                     SEXP cover) {
  clab::who = "kdsurf";
  csp::pts p = csp::xy(x, y, R_NilValue, csp::planar);
  const int n = p.n;
  if (n < 2) fail("at least 2 points are needed.");
  SEXP wv = PROTECT(Rf_isNull(w) ? R_NilValue : num(w, "w", n, 0.0));
  double hh[2];
  if (Rf_isNull(h)) {
    hh[0] = nrd(p.x, n);
    hh[1] = nrd(p.y, n);
    if (!(hh[0] > 0 && hh[1] > 0)) {
      fail("the points do not spread along both axes: give `h`.");
    }
  } else {
    SEXP hv = PROTECT(num(h, "h", -1, 0.0, true));
    if (XLENGTH(hv) != 1 && XLENGTH(hv) != 2) fail("`h` must be one or two positive numbers.");
    hh[0] = REAL(hv)[0];
    hh[1] = REAL(hv)[XLENGTH(hv) - 1];
    UNPROTECT(1);
  }
  static const char *const eo[] = {"none", "uniform", "diggle"};
  const int ed = csp::pick(edge, "edge", eo, 3);
  double wn[4];
  const bool hw = !Rf_isNull(win);
  if (hw) csp::win(win, wn);
  if (ed > 0 && !hw) fail("an edge correction needs `win`.");
  const double *pw = Rf_isNull(wv) ? nullptr : REAL(wv);
  SEXP out;
  if (!Rf_isNull(at)) {
    const bool own = TYPEOF(at) == STRSXP;
    if (own && !(XLENGTH(at) == 1 && std::strcmp(CHAR(STRING_ELT(at, 0)), "points") == 0)) {
      fail("`at` must be NULL, \"points\" or coordinates.");
    }
    csp::pts q = own ? p : csp::xy(at, R_NilValue, R_NilValue, csp::planar);
    SEXP z = PROTECT(real(q.n));
    kc(clab_kdsurf(p.x, p.y, pw, n, hh[0], hh[1], q.x, q.n, q.y, q.n, 0, hw ? wn : nullptr, ed,
                   own ? 1 : 0, REAL(z)),
       "a point lies where its window holds none of its kernel.");
    out = PROTECT(Rf_allocVector(VECSXP, 5));
    SET_VECTOR_ELT(out, 0, clab::dcopy(VECTOR_ELT(q.v, 0)));
    SET_VECTOR_ELT(out, 1, clab::dcopy(VECTOR_ELT(q.v, 1)));
    SET_VECTOR_ELT(out, 2, z);
    SET_VECTOR_ELT(out, 3, clab::nvec(2, xyn, hh));
    SET_VECTOR_ELT(out, 4, R_NilValue);
    UNPROTECT(own ? 2 : 3);
  } else {
    const int ax = clab::whole(nx, "nx", 2), ay = clab::whole(ny, "ny", 2);
    if (static_cast<double>(ax) * ay > 1e8) fail("`nx` times `ny` must be at most 1e8.");
    SEXP cv = PROTECT(num(cover, "cover", -1, 0.0, true));
    for (R_xlen_t q = 0; q < XLENGTH(cv); ++q) {
      if (!(REAL(cv)[q] < 1)) fail("`cover` must hold shares between 0 and 1.");
    }
    double b[4];
    if (hw) {
      std::memcpy(b, wn, sizeof b);
    } else {
      b[0] = *std::min_element(p.x, p.x + n) - 3 * hh[0];
      b[1] = *std::max_element(p.x, p.x + n) + 3 * hh[0];
      b[2] = *std::min_element(p.y, p.y + n) - 3 * hh[1];
      b[3] = *std::max_element(p.y, p.y + n) + 3 * hh[1];
    }
    SEXP gx = PROTECT(real(ax)), gy = PROTECT(real(ay));
    for (int a = 0; a < ax; ++a) REAL(gx)[a] = b[0] + a * (b[1] - b[0]) / (ax - 1);
    for (int a = 0; a < ay; ++a) REAL(gy)[a] = b[2] + a * (b[3] - b[2]) / (ay - 1);
    SEXP z = PROTECT(Rf_allocMatrix(REALSXP, ax, ay));
    kc(clab_kdsurf(p.x, p.y, pw, n, hh[0], hh[1], REAL(gx), ax, REAL(gy), ay, 1, hw ? wn : nullptr,
                   ed, 0, REAL(z)),
       "a point lies where its window holds none of its kernel.");
    // the smallest sets of cells holding each share of the surface's mass
    const R_xlen_t m = static_cast<R_xlen_t>(ax) * ay;
    std::vector<double> zs(REAL(z), REAL(z) + m);
    std::sort(zs.begin(), zs.end(), [](double a, double c) { return a > c; });
    double tot = 0;
    for (R_xlen_t q = 0; q < m; ++q) tot += zs[q];
    const double ca = (b[1] - b[0]) / (ax - 1) * (b[3] - b[2]) / (ay - 1);
    const int nc = static_cast<int>(XLENGTH(cv));
    static const char *const rn[] = {"cover", "level", "area", "share"};
    double *rc[4];
    SEXP reg = PROTECT(dframe(4, nc, rn, rc));
    for (int k = 0; k < nc; ++k) {
      double run = 0;
      R_xlen_t q = 0;
      while (q < m && run < REAL(cv)[k] * tot) run += zs[q++];
      rc[0][k] = REAL(cv)[k];
      rc[1][k] = zs[q > 0 ? q - 1 : 0];
      rc[2][k] = q * ca;
      rc[3][k] = static_cast<double>(q) / m;
    }
    out = PROTECT(Rf_allocVector(VECSXP, 5));
    SET_VECTOR_ELT(out, 0, gx);
    SET_VECTOR_ELT(out, 1, gy);
    SET_VECTOR_ELT(out, 2, z);
    SET_VECTOR_ELT(out, 3, clab::nvec(2, xyn, hh));
    SET_VECTOR_ELT(out, 4, reg);
    UNPROTECT(6);
  }
  PROTECT(out);
  static const char *const on[] = {"x", "y", "z", "h", "regions"};
  clab::nlist(out, on);
  UNPROTECT(3);
  return klass(out, "kdsurf");
}

SEXP C_clab_r_dbscn(SEXP x, SEXP y, SEXP eps, SEXP minpts, SEXP border, SEXP lonlat) {
  clab::who = "dbscn";
  csp::pts p = csp::xy(x, y, lonlat);
  const int n = p.n;
  const double e = clab::one(eps, "eps", 0.0, std::numeric_limits<double>::infinity());
  const int mp = clab::whole(minpts, "minpts", 1);
  const bool bd = csp::flag(border, "border");
  SEXP cl = PROTECT(Rf_allocVector(INTSXP, n));
  SEXP co = PROTECT(Rf_allocVector(LGLSXP, n));
  int k = 0;
  kc(clab_dbscn(p.x, p.y, n, e, mp, bd ? 1 : 0, INTEGER(cl), LOGICAL(co), &k),
     "the kernel refused the input.");
  SEXP sz = PROTECT(Rf_allocVector(INTSXP, k));
  std::fill(INTEGER(sz), INTEGER(sz) + k, 0);
  for (int i = 0; i < n; ++i) {
    if (INTEGER(cl)[i] >= 0) {
      ++INTEGER(sz)[INTEGER(cl)[i]];
      ++INTEGER(cl)[i];
    } else {
      INTEGER(cl)[i] = NA_INTEGER;
    }
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, cl);
  SET_VECTOR_ELT(out, 1, co);
  SET_VECTOR_ELT(out, 2, sz);
  static const char *const on[] = {"cl", "core", "sizes"};
  clab::nlist(out, on);
  UNPROTECT(5);
  return klass(out, "dbscn");
}

SEXP C_clab_r_hsper(SEXP x, SEXP y, SEXP prd, SEXP cell, SEXP top, SEXP org, SEXP lonlat) {
  clab::who = "hsper";
  csp::pts p = csp::xy(x, y, lonlat);
  const int n = p.n;
  csp::lab g = csp::codes(prd, "prd", n);
  const double cs = clab::one(cell, "cell", 0.0, std::numeric_limits<double>::infinity());
  const double tp = clab::val1(top);
  if (!(tp > 0 && tp <= 1)) fail("`top` must be one share above 0, at most 1.");
  double o[2];
  if (Rf_isNull(org)) {
    o[0] = *std::min_element(p.x, p.x + n);
    o[1] = *std::min_element(p.y, p.y + n);
  } else {
    SEXP ov = PROTECT(num(org, "org", 2));
    o[0] = REAL(ov)[0];
    o[1] = REAL(ov)[1];
    UNPROTECT(1);
  }
  const int np = g.g;
  SEXP cid = PROTECT(Rf_allocVector(INTSXP, n));
  SEXP ix = PROTECT(real(n)), iy = PROTECT(real(n));
  SEXP hot = PROTECT(Rf_allocVector(INTSXP, static_cast<R_xlen_t>(n) * np));
  SEXP per = PROTECT(real(9 * static_cast<R_xlen_t>(np)));
  SEXP cv = PROTECT(real(3 * static_cast<R_xlen_t>(n)));
  int nc = 0;
  kc(clab_hsper(p.x, p.y, INTEGER(g.code), n, np, cs, o[0], o[1], tp, INTEGER(cid), REAL(ix),
                REAL(iy), &nc, INTEGER(hot), REAL(per), REAL(cv)),
     "the cells are too many for the grid: choose a larger `cell` or an `org` nearer the points.");
  static const char *const pn[] = {"period", "n",      "active",  "nhalf",  "shalf",
                                   "hot",    "tiedin", "tiedout", "hotshr", "jacc"};
  SEXP pc = PROTECT(Rf_allocVector(VECSXP, 10));
  SET_VECTOR_ELT(pc, 0, g.lev);
  for (int j = 0; j < 9; ++j) {
    SEXP c = real(np);
    SET_VECTOR_ELT(pc, j + 1, c);
    for (int q = 0; q < np; ++q) REAL(c)[q] = REAL(per)[9 * static_cast<R_xlen_t>(q) + j];
  }
  SEXP pf = PROTECT(clab::frame(pc, pn));
  static const char *const cn[] = {"x", "y", "n", "nhot", "run"};
  double *cc[5];
  SEXP cf = PROTECT(dframe(5, nc, cn, cc));
  for (int j = 0; j < nc; ++j) {
    cc[0][j] = o[0] + (REAL(ix)[j] + 0.5) * cs;
    cc[1][j] = o[1] + (REAL(iy)[j] + 0.5) * cs;
    for (int q = 0; q < 3; ++q) cc[2 + q][j] = REAL(cv)[3 * static_cast<R_xlen_t>(j) + q];
  }
  csp::back(p, cc[0], cc[1], nc);
  SEXP hm = PROTECT(Rf_allocMatrix(LGLSXP, nc, np));
  for (int q = 0; q < np; ++q) {
    for (int j = 0; j < nc; ++j) {
      LOGICAL(hm)
      [j + static_cast<R_xlen_t>(q) * nc] = INTEGER(hot)[j + static_cast<R_xlen_t>(q) * nc];
    }
  }
  for (int i = 0; i < n; ++i) ++INTEGER(cid)[i];
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
  SET_VECTOR_ELT(out, 0, pf);
  SET_VECTOR_ELT(out, 1, cf);
  SET_VECTOR_ELT(out, 2, hm);
  SET_VECTOR_ELT(out, 3, cid);
  SET_VECTOR_ELT(out, 4, Rf_ScalarInteger(std::max(1, static_cast<int>(std::ceil(tp * nc)))));
  static const char *const on[] = {"periods", "cells", "hot", "cell", "k"};
  clab::nlist(out, on);
  UNPROTECT(14);
  return klass(out, "hsper");
}

SEXP C_clab_r_ebshr(SEXP y, SEXP v, SEXP w, SEXP b) {
  clab::who = "ebshr";
  SEXP yv = PROTECT(num(y, "y"));
  const int n = clab::ilen(static_cast<double>(XLENGTH(yv)), "y");
  if (n < 2) fail("at least 2 places are needed.");
  SEXP wv = PROTECT(Rf_isNull(w) ? R_NilValue : num(w, "w", n, 0.0));
  const double *pw = Rf_isNull(wv) ? nullptr : REAL(wv);
  double nv;
  if (Rf_isNull(v)) {
    double sw = 0, sy = 0;
    for (int i = 0; i < n; ++i) {
      sw += pw ? pw[i] : 1.0;
      sy += (pw ? pw[i] : 1.0) * REAL(yv)[i];
    }
    if (!(sw > 0)) fail("`w` must have a positive sum.");
    nv = sy / sw;
    if (!(nv >= 0)) fail("the mean of `y` is negative: give the noise variance `v`.");
  } else {
    SEXP vv = PROTECT(num(v, "v", 1, 0.0));
    nv = REAL(vv)[0];
    UNPROTECT(1);
  }
  const double bb = Rf_isNull(b) ? std::numeric_limits<double>::quiet_NaN() : clab::p01(b, "b");
  double fit[6];
  SEXP sh = PROTECT(real(n)), fa = PROTECT(real(n));
  kc(clab_ebshr(REAL(yv), pw, n, nv, bb, fit, REAL(sh), REAL(fa)), "`w` must have a positive sum.");
  static const char *const fn[] = {"mean", "total", "noise", "signal", "bstar", "b"};
  static const char *const en[] = {"y", "shrunk", "fall"};
  SEXP ec = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(ec, 0, clab::dcopy(yv));
  SET_VECTOR_ELT(ec, 1, sh);
  SET_VECTOR_ELT(ec, 2, fa);
  SEXP est = PROTECT(clab::frame(ec, en));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, clab::nvec(6, fn, fit));
  SET_VECTOR_ELT(out, 1, est);
  static const char *const on[] = {"fit", "est"};
  clab::nlist(out, on);
  UNPROTECT(7);
  return klass(out, "ebshr");
}

// ------------------------------------------------------------------ scans

SEXP C_clab_r_kscan(SEXP x, SEXP y, SEXP t, SEXP r, SEXP tw, SEXP minc, SEXP top, SEXP nsim,
                    SEXP seed, SEXP lonlat) {
  clab::who = "kscan";
  csp::pts p = csp::xy(x, y, lonlat);
  const int n = p.n;
  if (n < 3) fail("at least 3 incidents are needed.");
  SEXP tv = PROTECT(num(t, "t", n));
  SEXP rv = edges(r, "r");
  SEXP wv = PROTECT(num(tw, "tw", -1, 0.0, true));
  const int mc = clab::whole(minc, "minc", 1);
  const int tp = clab::whole(top, "top", 1);
  const int ns = clab::whole(nsim, "nsim", 0);
  const unsigned int sd = ns > 0 ? csp::seed(seed) : 0u;
  const int nr = static_cast<int>(XLENGTH(rv)), nw = static_cast<int>(XLENGTH(wv));
  const int tq = tp < n ? tp : n;
  SEXP o = PROTECT(real(10 * static_cast<R_xlen_t>(tq)));
  int nc = 0;
  kc(clab_kscan(p.x, p.y, REAL(tv), n, REAL(rv), nr, REAL(wv), nw, mc, ns, sd, tq, REAL(o), &nc),
     "the kernel refused the input.");
  static const char *const cn[] = {"x",   "y",     "r",    "start", "end", "n",
                                   "exp", "ncirc", "nwin", "rr",    "llr", "p"};
  double *c[12];
  SEXP cl = PROTECT(dframe(12, nc, cn, c));
  SEXP mb = PROTECT(Rf_allocVector(INTSXP, n));
  for (int i = 0; i < n; ++i) INTEGER(mb)[i] = NA_INTEGER;
  for (int k = nc - 1; k >= 0; --k) {
    const double *v = REAL(o) + 10 * k;
    const int ci = static_cast<int>(v[0]);
    c[0][k] = p.x[ci];
    c[1][k] = p.y[ci];
    for (int j = 0; j < 7; ++j) c[2 + j][k] = v[1 + j];
    // the relative risk inside against outside the cylinder
    c[9][k] = (v[4] / v[5]) / ((n - v[4]) / (n - v[5]));
    c[10][k] = v[8];
    c[11][k] = v[9];
    for (int i = 0; i < n; ++i) {
      const double dx = p.x[i] - p.x[ci], dy = p.y[i] - p.y[ci];
      if (dx * dx + dy * dy <= v[1] * v[1] && REAL(tv)[i] >= v[2] && REAL(tv)[i] < v[3]) {
        INTEGER(mb)[i] = k + 1;
      }
    }
  }
  csp::back(p, c[0], c[1], nc);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out, 0, cl);
  SET_VECTOR_ELT(out, 1, mb);
  SET_VECTOR_ELT(out, 2, Rf_ScalarInteger(n));
  SET_VECTOR_ELT(out, 3, Rf_ScalarInteger(ns));
  static const char *const on[] = {"cl", "memb", "n", "nsim"};
  clab::nlist(out, on);
  UNPROTECT(8);
  return klass(out, "kscan");
}

}  // extern "C"
