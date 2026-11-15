// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/krg.cpp and src/gla.cpp: hand values, exact interpolation,
// round trips and the refusals.
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "gla.h"
#include "tchk.h"

namespace {

void t_gla(tk &T) {
  // [[2, 1], [1, 2]]: eigenvalues 3 and 1
  double a[4] = {2, 1, 1, 2}, val[2], vec[4];
  tk_true(T, "gla.jac", gla::jacobi(a, 2, val, vec));
  tk_eq(T, "gla.l1", val[0], 3, 1e-14);
  tk_eq(T, "gla.l2", val[1], 1, 1e-14);
  tk_eq(T, "gla.v", std::fabs(vec[0]), std::sqrt(0.5), 1e-14);
  double bad[4] = {1, NAN, NAN, 1};
  tk_true(T, "gla.jnan", !gla::jacobi(bad, 2, val, vec));
  double s[4] = {1, 2, 2, 4}, b[2] = {1, 1};
  int piv[2];
  tk_true(T, "gla.sing", !gla::lusolve(s, 2, b, 1, piv));
  double l[4] = {4, 2, 2, 3}, r[2] = {1, 2};
  tk_true(T, "gla.chol", gla::chol(l, 2));
  tk_eq(T, "gla.l11", l[0], 2, 1e-15);
  gla::fsolve(l, 2, r, 1);
  gla::bsolve(l, 2, r, 1);
  tk_eq(T, "gla.cs0", r[0], -0.125, 1e-15);
  tk_eq(T, "gla.cs1", r[1], 0.75, 1e-15);
  double np[4] = {1, 2, 2, 1};
  tk_true(T, "gla.npd", !gla::chol(np, 2));
  // the minimisers on a quadratic
  double x[2] = {1, 1}, f = 0;
  int nf = 0;
  const int rc = gla::nmmin(2, x, &f, [](const double *p) { return (p[0] - 2) * (p[0] - 2) + 3 * (p[1] + 1) * (p[1] + 1); },
                            1e-12, 2000, &nf);
  tk_ieq(T, "gla.nmrc", rc, 0);
  tk_eq(T, "gla.nm0", x[0], 2, 1e-5);
  tk_eq(T, "gla.nm1", x[1], -1, 1e-5);
  double y0[1] = {0};
  tk_ieq(T, "gla.nmmax", gla::nmmin(1, y0, &f, [](const double *p) { return -p[0]; }, 1e-12, 20, &nf), 1);
  double y1[1] = {0};
  tk_ieq(T, "gla.nmnan", gla::nmmin(1, y1, &f, [](const double *) { return NAN; }, 1e-8, 20, &nf), 2);
  double y2[1] = {0};
  gla::nmmin(1, y2, &f, [](const double *p) { return p[0] > 0.05 ? NAN : -p[0]; }, 1e-8, 200, &nf);
  tk_le(T, "gla.nmbig", y2[0], 0.05);
  tk_eq(T, "gla.fmin", gla::fmin([](double u) { return (u - 0.3) * (u - 0.3); }, 0, 1, 1e-10), 0.3, 1e-8);
}

}  // namespace

void t_krg(tk &T, int big) {
  (void)big;
  t_gla(T);
  // three points on a line, one pair per bin of (0, 1], (1, 2], (2, 3]
  const double x[3] = {0, 1, 3}, y[3] = {0, 0, 0}, v[3] = {0, 1, 3};
  double d[3], g[3], np[3];
  int nk = 0;
  tk_ieq(T, "vg.ok", clab_vgemp(x, y, 3, v, 3, 3, d, g, np, &nk), CLAB_OK);
  tk_ieq(T, "vg.nk", nk, 3);
  tk_eq(T, "vg.d0", d[0], 0.5, 1e-15);
  tk_eq(T, "vg.g0", g[0], 0.5, 1e-15);
  tk_eq(T, "vg.g1", g[1], 2, 1e-15);
  tk_eq(T, "vg.g2", g[2], 4.5, 1e-15);
  tk_ieq(T, "vg.def", clab_vgemp(x, y, 3, v, 2, NAN, d, g, np, &nk), CLAB_OK);
  tk_ieq(T, "vg.nk1", nk, 1);
  tk_ieq(T, "vg.n", clab_vgemp(x, y, 1, v, 3, 3, d, g, np, &nk), CLAB_EARG);
  const double vn[3] = {0, NAN, 1}, z0[3] = {0, 0, 0};
  tk_ieq(T, "vg.nan", clab_vgemp(x, y, 3, vn, 3, 3, d, g, np, &nk), CLAB_EDOM);
  tk_ieq(T, "vg.same", clab_vgemp(z0, z0, 3, v, 3, NAN, d, g, np, &nk), CLAB_EDOM);
  // kriging interpolates a data point exactly without a nugget
  const double px[4] = {0, 1, 0, 1}, py[4] = {0, 0, 1, 1}, pv[4] = {1, 2, 3, 5};
  double pr[2], vr[2];
  const double ax[2] = {1, 0.5}, ay[2] = {0, 0.5};
  for (int m = 0; m < 3; ++m) {
    tk_ieq(T, "ok.ok", clab_okrge(px, py, 4, pv, ax, ay, 2, m, 0, 1, 2, pr, vr), CLAB_OK);
    tk_eq(T, "ok.exact", pr[0], 2, 1e-12);
    tk_eq(T, "ok.var0", vr[0], 0, 1e-12);
    // the centre: equal weights by symmetry
    tk_eq(T, "ok.mid", pr[1], 2.75, 1e-12);
  }
  tk_ieq(T, "ok.model", clab_okrge(px, py, 4, pv, ax, ay, 2, 3, 0, 1, 2, pr, vr), CLAB_EARG);
  tk_ieq(T, "ok.rg", clab_okrge(px, py, 4, pv, ax, ay, 2, 1, 0, 1, 0, pr, vr), CLAB_EDOM);
  const double dx[2] = {0, 0};
  tk_ieq(T, "ok.sing", clab_okrge(dx, dx, 2, pv, ax, ay, 2, 1, 0, 1, 2, pr, vr), CLAB_EDOM);
  // block kriging of a block that is one data point: that value, no error
  const double lo[2] = {0, 0}, hi[2] = {1, 1};
  double bx[9], by[9];
  tk_ieq(T, "bk.pts", clab_bkpts(lo, hi, 9, bx, by), CLAB_OK);
  tk_eq(T, "bk.mid", bx[4], 0.5, 1e-15);
  tk_eq(T, "bk.last", by[8], 1, 1e-15);
  tk_ieq(T, "bk.nq", clab_bkpts(lo, hi, 0, bx, by), CLAB_EARG);
  double one1[1] = {0};
  tk_ieq(T, "bk.k1", clab_bkpts(lo, hi, 1, one1, one1 + 0), CLAB_OK);
  const double qx[1] = {1}, qy[1] = {1};
  const int off[3] = {0, 1, 9};
  double est[2], se[2];
  tk_ieq(T, "bk.ok", clab_bkrig(px, py, 4, pv, qx, qy, off, 1, 0, 1, 1, est, se), CLAB_OK);
  tk_eq(T, "bk.est", est[0], 5, 1e-12);
  tk_eq(T, "bk.se", se[0], 0, 1e-6);
  tk_ieq(T, "bk.box", clab_bkrig(px, py, 4, pv, bx, by, off, 1, 0, 1, 1, est, se), CLAB_OK);
  tk_ieq(T, "bk.empty", clab_bkrig(px, py, 4, pv, bx, by, off, 2, 0, 1, 1, est, se), CLAB_OK);
  const int bad[2] = {0, 0};
  tk_ieq(T, "bk.none", clab_bkrig(px, py, 4, pv, bx, by, bad, 1, 0, 1, 1, est, se), CLAB_EARG);
  tk_ieq(T, "bk.n", clab_bkrig(px, py, 0, pv, bx, by, off, 1, 0, 1, 1, est, se), CLAB_EARG);
  tk_ieq(T, "bk.rg", clab_bkrig(px, py, 4, pv, bx, by, off, 1, 0, 1, -1, est, se), CLAB_EDOM);
  tk_ieq(T, "bk.sing", clab_bkrig(dx, dx, 2, pv, bx, by, off, 1, 0, 1, 1, est, se), CLAB_EDOM);
  // leave-one-out: by symmetry of the square, the residual of a corner is
  // its value less the kriged value from the other three
  double res[4];
  tk_ieq(T, "cv.ok", clab_krigcv(px, py, 4, pv, 0, 1, 1, res), CLAB_OK);
  tk_true(T, "cv.fin", std::isfinite(res[0]) && std::isfinite(res[3]));
  tk_ieq(T, "cv.n", clab_krigcv(px, py, 1, pv, 0, 1, 1, res), CLAB_EARG);
  tk_ieq(T, "cv.rg", clab_krigcv(px, py, 4, pv, 0, 1, 0, res), CLAB_EDOM);
  const double cx[3] = {0, 0, 1};
  tk_ieq(T, "cv.sing", clab_krigcv(cx, cx, 3, pv, 0, 1, 1, res), CLAB_EDOM);
  // a plane is fitted exactly
  const double tv[4] = {1, 3, 0, 2};  // 1 + 2 x - y
  double cf[10], cs[10], r2 = 0;
  tk_ieq(T, "tr.ok", clab_trsurf(px, py, 4, tv, 1, cf, cs, &r2), CLAB_OK);
  tk_eq(T, "tr.b0", cf[0], 1, 1e-13);
  tk_eq(T, "tr.b1", cf[1], 2, 1e-13);
  tk_eq(T, "tr.b2", cf[2], -1, 1e-13);
  tk_eq(T, "tr.r2", r2, 1, 1e-13);
  const double c4[4] = {2, 2, 2, 2};
  tk_ieq(T, "tr.const", clab_trsurf(px, py, 4, c4, 0, cf, cs, &r2), CLAB_OK);
  tk_eq(T, "tr.r2c", r2, 1, 0);
  tk_ieq(T, "tr.ord", clab_trsurf(px, py, 4, tv, 4, cf, cs, &r2), CLAB_EARG);
  tk_ieq(T, "tr.np", clab_trsurf(px, py, 4, tv, 2, cf, cs, &r2), CLAB_EARG);
  const double lx[4] = {0, 1, 2, 3};
  tk_ieq(T, "tr.coll", clab_trsurf(lx, lx, 4, tv, 1, cf, cs, &r2), CLAB_EDOM);
  // spatial linear models: an intercept and a slope
  const int n = 8;
  double sx[n], sy[n], sv[n], z[2 * n];
  for (int i = 0; i < n; ++i) {
    sx[i] = i % 3 + 0.1 * i;
    sy[i] = i / 3 + 0.05 * i * i;
    z[i] = 1;
    z[n + i] = sx[i];
    sv[i] = 1 + 2 * sx[i] + 0.3 * std::sin(3.0 * i);
  }
  double b[2], sb[2], o[4];
  tk_ieq(T, "sg.ok", clab_sglmx(z, n, 2, sv, sx, sy, b, sb, o), CLAB_OK);
  tk_near(T, "sg.b1", b[1], 2, 0.5);
  tk_gt(T, "sg.phi", o[1], 0);
  tk_ieq(T, "sm.ok", clab_smixm(z, n, 2, sv, sx, sy, b, sb, o), CLAB_OK);
  tk_near(T, "sm.b1", b[1], 2, 0.5);
  double zc[2 * n];
  for (int i = 0; i < n; ++i) zc[i] = zc[n + i] = 1;
  tk_ieq(T, "sg.coll", clab_sglmx(zc, n, 2, sv, sx, sy, b, sb, o), CLAB_EDOM);
  tk_ieq(T, "sm.coll", clab_smixm(zc, n, 2, sv, sx, sy, b, sb, o), CLAB_EDOM);
  tk_ieq(T, "sg.n", clab_sglmx(z, 1, 1, sv, sx, sy, b, sb, o), CLAB_EARG);
  tk_ieq(T, "sm.n", clab_smixm(z, 1, 1, sv, sx, sy, b, sb, o), CLAB_EARG);
  // a perfect fit makes the REML residual variance 0: the guard
  double sp[n];
  for (int i = 0; i < n; ++i) sp[i] = 1 + 2 * sx[i];
  tk_ieq(T, "sm.exact", clab_smixm(z, n, 2, sp, sx, sy, b, sb, o), CLAB_OK);
  // variogram fits
  double vf[6], vw[3];
  tk_ieq(T, "vf.ok", clab_vgfit(sx, sy, n, sv, 1, vf), CLAB_OK);
  tk_true(T, "vf.fin", std::isfinite(vf[3]) && vf[2] > 0);
  tk_ieq(T, "vf.m0", clab_vgfit(sx, sy, n, sv, 0, vf), CLAB_OK);
  tk_ieq(T, "vf.n", clab_vgfit(sx, sy, 2, sv, 1, vf), CLAB_EARG);
  double xn[n];
  for (int i = 0; i < n; ++i) xn[i] = i == 0 ? NAN : sx[i];
  tk_ieq(T, "vf.pnan", clab_vgfit(xn, sy, n, sv, 1, vf), CLAB_EARG);
  double cv8[n];
  for (int i = 0; i < n; ++i) cv8[i] = 4;
  tk_ieq(T, "vf.const", clab_vgfit(sx, sy, n, cv8, 1, vf), CLAB_EDOM);
  double vnn[n];
  for (int i = 0; i < n; ++i) vnn[i] = i == 2 ? NAN : sv[i];
  tk_ieq(T, "vf.vnan", clab_vgfit(sx, sy, n, vnn, 1, vf), CLAB_EDOM);
  tk_ieq(T, "vw.ok", clab_vgwls(sx, sy, n, sv, 2, vw), CLAB_OK);
  tk_ieq(T, "vw.m", clab_vgwls(sx, sy, n, sv, 5, vw), CLAB_EARG);
  tk_ieq(T, "vw.same", clab_vgwls(z0, z0, 3, v, 1, vw), CLAB_EDOM);
  // Web Mercator: a round trip
  const double lon[2] = {-79.4, 120}, lat[2] = {43.7, -33};
  double mx[2], my[2], bx2[2], by2[2];
  tk_ieq(T, "wm.fwd", clab_wmerc(lon, lat, 2, 0, mx, my), CLAB_OK);
  tk_ieq(T, "wm.inv", clab_wmerc(mx, my, 2, 1, bx2, by2), CLAB_OK);
  tk_eq(T, "wm.lon", bx2[1], 120, 1e-12);
  tk_eq(T, "wm.lat", by2[0], 43.7, 1e-12);
  const double pole[1] = {90};
  tk_ieq(T, "wm.pole", clab_wmerc(lon, pole, 1, 0, mx, my), CLAB_EDOM);
  tk_ieq(T, "wm.nan", clab_wmerc(vn + 1, lat, 1, 1, mx, my), CLAB_EDOM);
  // density: the corners of a square are alike
  double pd[4];
  tk_ieq(T, "sd.ok", clab_spdens(px, py, 4, 0, pd), CLAB_OK);
  tk_eq(T, "sd.eq", pd[2], 0.25, 1e-15);
  tk_ieq(T, "sd.ll", clab_spdens(px, py, 4, 1, pd), CLAB_OK);
  tk_ieq(T, "sd.n", clab_spdens(px, py, 1, 0, pd), CLAB_EARG);
  tk_ieq(T, "sd.same", clab_spdens(z0, z0, 3, 0, pd), CLAB_EDOM);
  const double lo2[2] = {0, 1}, la2[2] = {0, 100};
  double pf[2];
  tk_ieq(T, "sd.lat", clab_spdens(lo2, la2, 2, 1, pf), CLAB_EDOM);
  // sampling: all of them, each once
  int idx[4];
  tk_ieq(T, "ss.ok", clab_spsamp(pd, 4, 4, 7, idx), CLAB_OK);
  tk_eq(T, "ss.sum", idx[0] + idx[1] + idx[2] + idx[3], 10, 0);
  tk_ieq(T, "ss.uni", clab_spsamp(nullptr, 4, 2, 7, idx), CLAB_OK);
  const double zp[3] = {0, 0, 0}, np3[3] = {1, -1, 1};
  tk_ieq(T, "ss.zero", clab_spsamp(zp, 3, 2, 7, idx), CLAB_OK);
  tk_ieq(T, "ss.neg", clab_spsamp(np3, 3, 2, 7, idx), CLAB_EDOM);
  tk_ieq(T, "ss.big", clab_spsamp(pd, 4, 5, 7, idx), CLAB_EARG);
}
