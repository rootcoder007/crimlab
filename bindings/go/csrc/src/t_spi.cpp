// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for idw.cpp, gtw.cpp and arl.cpp: inverse-distance
// weights of whole distances by hand, a target on a monitor, the k
// nearest and the radius, zonal means; the weighted regression reducing
// to least squares with a flat kernel and fitting a line exactly; rules
// of a four-row table by hand; and the failure paths. Part of the
// self-test (selftest.cpp).
#include <cmath>

#include "crimlab.h"
#include "tchk.h"

namespace {

void idw(tk &tt) {
  // monitors at 0, 1 and 3 on a line with values 6, 9, 12; the target at
  // 2 is 2, 1 and 1 away: weights 1/d^2 = 1/4, 1, 1, so (6/4 + 9 + 12) /
  // 2.25 = 10
  const double mx[3] = {0, 1, 3}, my[3] = {0, 0, 0}, mv[3] = {6, 9, 12};
  const double tx[3] = {2, 1, 50}, ty[3] = {0, 0, 0};
  double o[3], u[3];
  tk_ieq(tt, "spi.idw.rc", clab_idwint(3, mx, my, mv, 3, tx, ty, 2, 0, 10, 0, nullptr, 0, o, u,
                                       nullptr, nullptr),
         CLAB_OK);
  tk_near(tt, "spi.idw.v", o[0], 10, 1e-14);
  tk_near(tt, "spi.idw.used", u[0], 3, 0);
  tk_near(tt, "spi.idw.on", o[1], 9, 0);
  tk_ieq(tt, "spi.idw.far", std::isnan(o[2]) ? 1 : 0, 1);
  // the two nearest of the target at 2 (distances 1, 1; the farther 2
  // left out), plain mean with p = 0
  clab_idwint(3, mx, my, mv, 1, tx, ty, 0, 2, INFINITY, 0, nullptr, 0, o, u, nullptr, nullptr);
  tk_near(tt, "spi.idw.k", o[0], 10.5, 1e-14);
  // zones: the first two targets in zone 0, the third (NaN) alone in 1
  const int zn[3] = {0, 0, 1};
  double zo[2], zc[2];
  clab_idwint(3, mx, my, mv, 3, tx, ty, 2, 0, 10, 0, zn, 2, o, u, zo, zc);
  tk_near(tt, "spi.idw.z0", zo[0], 9.5, 1e-14);
  tk_near(tt, "spi.idw.zn", zc[0], 2, 0);
  tk_ieq(tt, "spi.idw.z1", std::isnan(zo[1]) ? 1 : 0, 1);
  tk_ieq(tt, "spi.idw.arg", clab_idwint(3, mx, my, mv, 1, tx, ty, -1, 0, 10, 0, nullptr, 0, o, u,
                                        nullptr, nullptr),
         CLAB_EARG);
  const double nx[3] = {0, NAN, 3};
  tk_ieq(tt, "spi.idw.dom", clab_idwint(3, nx, my, mv, 1, tx, ty, 2, 0, 10, 0, nullptr, 0, o, u,
                                        nullptr, nullptr),
         CLAB_EDOM);
}

void gtw(tk &tt) {
  // y = 1 + 2 x exactly: every local fit is the line, residuals 0
  const int n = 6;
  const double x[12] = {1, 1, 1, 1, 1, 1, 0, 1, 2, 3, 4, 5};
  const double y[6] = {1, 3, 5, 7, 9, 11}, u[6] = {0, 1, 2, 0, 1, 2}, v[6] = {0, 0, 0, 1, 1, 1};
  double b[12], se[12], f[6], h[6], l[6], st[7];
  tk_ieq(tt, "spi.gtw.rc",
         clab_gtwreg(n, 2, x, y, u, v, nullptr, 1.5, 0, 0, b, se, f, h, l, st), CLAB_OK);
  tk_near(tt, "spi.gtw.b0", b[3], 1, 1e-12);
  tk_near(tt, "spi.gtw.b1", b[9], 2, 1e-12);
  tk_near(tt, "spi.gtw.rss", st[0], 0, 1e-20);
  // a flat kernel (h huge): the hat diagonal of least squares on (1, x),
  // 1/6 + (x - 2.5)^2 / 17.5, summing to 2
  clab_gtwreg(n, 2, x, y, u, v, nullptr, 1e9, 0, 0, b, se, f, h, l, st);
  tk_near(tt, "spi.gtw.hat0", h[0], 1.0 / 6 + 6.25 / 17.5, 1e-9);
  tk_near(tt, "spi.gtw.trs", st[1], 2, 1e-8);
  tk_ieq(tt, "spi.gtw.k", clab_gtwreg(n, 6, x, y, u, v, nullptr, 1, 0, 0, b, se, f, h, l, st),
         CLAB_EARG);
  tk_ieq(tt, "spi.gtw.h", clab_gtwreg(n, 2, x, y, u, v, nullptr, 0, 0, 0, b, se, f, h, l, st),
         CLAB_EARG);
  // one regressor constant: X'WX singular
  const double xs[12] = {1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2};
  tk_ieq(tt, "spi.gtw.sing",
         clab_gtwreg(n, 2, xs, y, u, v, nullptr, 1, 0, 0, b, se, f, h, l, st), CLAB_EDOM);
}

void arl(tk &tt) {
  // x = (0, 0, 1, 1), y = (0, 1, 1, 1): x=0 => y=1 has support 1/4,
  // confidence 1/2, coverage 1/2, lift 1/2 / 3/4 = 2/3; x=1 => y=1
  // support 1/2, confidence 1, lift 4/3
  const int xc[4] = {0, 0, 1, 1}, yc[4] = {0, 1, 1, 1};
  int nr = 0, lh[12], rh[12];
  double o[60];
  tk_ieq(tt, "spi.arl.rc", clab_asrule(4, xc, 2, yc, 2, 0, 0, 0, &nr, lh, rh, o), CLAB_OK);
  tk_ieq(tt, "spi.arl.n", nr, 6);
  // the rules in order: x0=>y0, y0=>x0, x0=>y1, y1=>x0, x1=>y1, y1=>x1
  // (x1 and y0 never together); x0=>y1 third
  tk_ieq(tt, "spi.arl.lhs", lh[2], 0);
  tk_ieq(tt, "spi.arl.rhs", rh[2], 3);
  tk_near(tt, "spi.arl.sup", o[10], 0.25, 0);
  tk_near(tt, "spi.arl.conf", o[11], 0.5, 0);
  tk_near(tt, "spi.arl.lift", o[13], 2.0 / 3, 1e-15);
  // at least 30% support and the empty left-hand side: {} => each of
  // the four items with support 1/2, 1/2, 1/4 (out), 3/4; x0=>y1 out
  clab_asrule(4, xc, 2, yc, 2, 0.3, 0, 1, &nr, lh, rh, o);
  tk_ieq(tt, "spi.arl.n2", nr, 5);
  tk_ieq(tt, "spi.arl.empty", lh[0], -1);
  tk_ieq(tt, "spi.arl.code", clab_asrule(4, xc, 1, yc, 2, 0, 0, 0, &nr, lh, rh, o), CLAB_EARG);
}

void ifx(tk &tt) {
  // two controls (1 2 3; 3 4 5) and one unit treated at the third time
  // (2 3 10): mu 3, time effects -1 0 1, the unit's effect 0, its
  // counterfactual at the third time 4 and its effect 6
  const double y[9] = {1, 2, 3, 3, 4, 5, 2, 3, 10};
  const int a[3] = {3, 3, 2};
  double beta[1], y0[9], fac[3], cv[1], ef[6], bef[1], st[8];
  tk_ieq(tt, "spi.ifx.rc", clab_ifesc(y, y, 0, a, 3, 3, 0, 0, 100, 1e-8, nullptr, 0, beta, y0, fac, cv, ef, bef, st),
         CLAB_OK);
  tk_near(tt, "spi.ifx.y0", y0[8], 4.0, 1e-12);
  tk_near(tt, "spi.ifx.att", ef[0], 6.0, 1e-12);
  tk_near(tt, "spi.ifx.ne", st[7], 3.0, 0);
  // no control left: every unit treated but one
  const int b[3] = {3, 2, 2};
  tk_ieq(tt, "spi.ifx.ctl", clab_ifesc(y, y, 0, b, 3, 3, 0, 0, 100, 1e-8, nullptr, 0, beta, y0, fac, cv, ef, bef, st),
         CLAB_EARG);
}

void cem(tk &tt) {
  // one covariate, cells 0 0 0 1 1 2: cell 0 has treated 1, controls 2;
  // cell 1 treated 1, no control (dropped); cell 2 a control alone
  // (dropped). Kept: m_T 1, m_C 2; the controls' weight (1/2)(2/1) = 1
  const int code[6] = {0, 0, 0, 1, 1, 2}, tr[6] = {1, 0, 0, 1, 1, 0};
  int s[6];
  double w[6], st[6];
  tk_ieq(tt, "spi.cem.rc", clab_cemtch(code, 6, 1, tr, code, 1, s, w, st), CLAB_OK);
  tk_ieq(tt, "spi.cem.s", s[0] * 100 + s[3] * 10 + s[5], 100);
  tk_near(tt, "spi.cem.w", w[1], 1.0, 0);
  // L1 before: treated shares 1/3 2/3 0, controls 2/3 0 1/3: (1/3 + 2/3
  // + 1/3) / 2 = 2/3; after 0
  tk_near(tt, "spi.cem.l1", st[0], 2.0 / 3, 1e-15);
  tk_near(tt, "spi.cem.l1a", st[1], 0.0, 0);
  tk_near(tt, "spi.cem.nk", st[2], 1.0, 0);
}

}  // namespace

void t_spi(tk &t, int big) {
  (void)big;
  idw(t);
  gtw(t);
  arl(t);
  ifx(t);
  cem(t);
}
