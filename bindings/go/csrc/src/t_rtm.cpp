// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for clab_crchyp, clab_hitscr and clab_rtmlay
// (src/rtm.cpp) against geometry done by hand. Part of the self-test
// (selftest.cpp).
#include <cmath>

#include "crimlab.h"
#include "tchk.h"

namespace {

void crchyp(tk &t) {
  // (0,0) and (4,0) are farthest apart: centre (2,0), radius 2
  const double x[3] = {0, 4, 2}, y[3] = {0, 0, 1}, hx[2] = {2, 5}, hy[2] = {-1, 5};
  double st[4];
  int pr[2], mar[2];
  tk_ieq(t, "crchyp.rc", clab_crchyp(x, y, 3, hx, hy, 2, st, pr, mar), CLAB_OK);
  tk_eq(t, "crchyp.cx", st[0], 2, 0);
  tk_eq(t, "crchyp.cy", st[1], 0, 0);
  tk_eq(t, "crchyp.r", st[2], 2, 0);
  tk_eq(t, "crchyp.share", st[3], 1, 0);
  tk_ieq(t, "crchyp.p0", pr[0], 0);
  tk_ieq(t, "crchyp.p1", pr[1], 1);
  tk_ieq(t, "crchyp.m0", mar[0], 1);
  tk_ieq(t, "crchyp.m1", mar[1], 0);
  // a site outside the circle; the first of two equal pairs is kept
  const double x2[4] = {0, 2, 0, 2}, y2[4] = {0, 0, 2, 2};
  clab_crchyp(x2, y2, 4, nullptr, nullptr, 0, st, pr, mar);
  tk_ieq(t, "crchyp.tie0", pr[0], 0);
  tk_ieq(t, "crchyp.tie1", pr[1], 3);
  tk_eq(t, "crchyp.r2", st[2], std::sqrt(2.0), 1e-15);
  const double x3[3] = {0, 4, 4}, y3[3] = {0, 0, 4};
  clab_crchyp(x3, y3, 3, nullptr, nullptr, 0, st, pr, mar);
  tk_eq(t, "crchyp.share3", st[3], 1, 0);
  const double x4[4] = {0, 10, 5, 5}, y4[4] = {0, 0, 5.1, 0};
  clab_crchyp(x4, y4, 4, nullptr, nullptr, 0, st, pr, mar);
  tk_eq(t, "crchyp.share4", st[3], 0.75, 0);
  const double bad[2] = {0, NAN};
  tk_ieq(t, "crchyp.n", clab_crchyp(x, y, 1, hx, hy, 0, st, pr, mar), CLAB_EARG);
  tk_ieq(t, "crchyp.nh", clab_crchyp(x, y, 2, hx, hy, -1, st, pr, mar), CLAB_EARG);
  tk_ieq(t, "crchyp.x", clab_crchyp(bad, y, 2, hx, hy, 0, st, pr, mar), CLAB_EDOM);
  tk_ieq(t, "crchyp.y", clab_crchyp(x, bad, 2, hx, hy, 0, st, pr, mar), CLAB_EDOM);
  tk_ieq(t, "crchyp.hx", clab_crchyp(x, y, 2, bad, hy, 2, st, pr, mar), CLAB_EDOM);
  tk_ieq(t, "crchyp.hy", clab_crchyp(x, y, 2, hx, bad, 2, st, pr, mar), CLAB_EDOM);
}

void hitscr(tk &t) {
  const double s[5] = {0.9, 0.5, 0.7, 0.5, 0.1};
  const int h[2] = {1, 0};
  int ab[2], ti[2];
  double hit[2], hitt[2];
  tk_ieq(t, "hitscr.rc", clab_hitscr(s, 5, h, 2, ab, ti, hit, hitt), CLAB_OK);
  tk_ieq(t, "hitscr.above", ab[0], 2);
  tk_ieq(t, "hitscr.ties", ti[0], 2);
  tk_eq(t, "hitscr.hit", hit[0], 40, 0);
  tk_eq(t, "hitscr.hitt", hitt[0], 80, 0);
  tk_eq(t, "hitscr.top", hit[1], 0, 0);
  tk_eq(t, "hitscr.topt", hitt[1], 20, 0);
  const int hb[1] = {5};
  const double sb[2] = {1, INFINITY};
  tk_ieq(t, "hitscr.n", clab_hitscr(s, 0, h, 1, ab, ti, hit, hitt), CLAB_EARG);
  tk_ieq(t, "hitscr.nh", clab_hitscr(s, 5, h, 0, ab, ti, hit, hitt), CLAB_EARG);
  tk_ieq(t, "hitscr.home", clab_hitscr(s, 5, hb, 1, ab, ti, hit, hitt), CLAB_EARG);
  tk_ieq(t, "hitscr.s", clab_hitscr(sb, 2, h, 1, ab, ti, hit, hitt), CLAB_EDOM);
}

void rtmlay(tk &t) {
  // cells (0,0), (3,0), (10,0); features (1,0) and (2,0)
  const double cx[3] = {0, 3, 10}, cy[3] = {0, 0, 0}, fx[2] = {1, 2}, fy[2] = {0, 0};
  double o[3];
  tk_ieq(t, "rtmlay.rc", clab_rtmlay(cx, cy, 3, fx, fy, 2, 1.5, 0, o), CLAB_OK);
  tk_eq(t, "rtmlay.p0", o[0], 1, 0);
  tk_eq(t, "rtmlay.p2", o[2], 0, 0);
  clab_rtmlay(cx, cy, 3, fx, fy, 2, 2, 1, o);
  tk_eq(t, "rtmlay.d0", o[0], 2, 0);
  tk_eq(t, "rtmlay.d1", o[1], 2, 0);
  tk_eq(t, "rtmlay.d2", o[2], 0, 0);
  // the kernel sums, standardised by hand
  clab_rtmlay(cx, cy, 3, fx, fy, 2, 1, 2, o);
  double k[3], m = 0, v = 0;
  for (int i = 0; i < 3; ++i) {
    const double a = cx[i] - 1, b = cx[i] - 2;
    k[i] = std::exp(-0.5 * a * a) + std::exp(-0.5 * b * b);
    m += k[i] / 3;
  }
  for (int i = 0; i < 3; ++i) v += (k[i] - m) * (k[i] - m) / 2;
  for (int i = 0; i < 3; ++i) tk_eq(t, "rtmlay.kde", o[i], (k[i] - m) / std::sqrt(v), 1e-14);
  // no variation: zeros
  const double c2[2] = {1, 1};
  clab_rtmlay(c2, c2, 2, fx, fy, 1, 1, 2, o);
  tk_eq(t, "rtmlay.flat", o[1], 0, 0);
  tk_ieq(t, "rtmlay.nc", clab_rtmlay(cx, cy, 0, fx, fy, 2, 1, 0, o), CLAB_EARG);
  tk_ieq(t, "rtmlay.nf", clab_rtmlay(cx, cy, 3, fx, fy, 0, 1, 0, o), CLAB_EARG);
  tk_ieq(t, "rtmlay.kind", clab_rtmlay(cx, cy, 3, fx, fy, 2, 1, 3, o), CLAB_EARG);
  tk_ieq(t, "rtmlay.kneg", clab_rtmlay(cx, cy, 3, fx, fy, 2, 1, -1, o), CLAB_EARG);
  tk_ieq(t, "rtmlay.kde1", clab_rtmlay(cx, cy, 1, fx, fy, 2, 1, 2, o), CLAB_EARG);
  tk_ieq(t, "rtmlay.rad", clab_rtmlay(cx, cy, 3, fx, fy, 2, 0, 0, o), CLAB_EDOM);
  const double bad[3] = {0, NAN, 1};
  tk_ieq(t, "rtmlay.cx", clab_rtmlay(bad, cy, 3, fx, fy, 2, 1, 0, o), CLAB_EDOM);
  tk_ieq(t, "rtmlay.cy", clab_rtmlay(cx, bad, 3, fx, fy, 2, 1, 0, o), CLAB_EDOM);
  tk_ieq(t, "rtmlay.fx", clab_rtmlay(cx, cy, 3, bad, fy, 2, 1, 0, o), CLAB_EDOM);
  tk_ieq(t, "rtmlay.fy", clab_rtmlay(cx, cy, 3, fx, bad, 2, 1, 0, o), CLAB_EDOM);
}

}  // namespace

void t_rtm(tk &t, int big) {
  (void)big;
  crchyp(t);
  hitscr(t);
  rtmlay(t);
}
