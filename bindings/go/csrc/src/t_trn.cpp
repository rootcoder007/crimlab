// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for the series kernels (src/trn.cpp): Mann-Kendall's S,
// its variance and exact p-value by hand (with ties), the Theil-Sen slope,
// the step scan's exact permutation p-value and Pettitt's statistic on a
// series whose orderings can be counted, CUSUM by hand, and the
// cross-correlations of two lines. Part of the self-test (selftest.cpp).
#include <cmath>

#include "crimlab.h"
#include "tchk.h"

namespace {

void mk(tk &t) {
  static double w[400];
  static int iw[20];
  double st[10];
  const double x[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  // 1..5: S = 10, variance 5 4 15 / 18, one ordering in 120 each way
  const double y[5] = {1, 2, 3, 4, 5};
  tk_ieq(t, "mk.rc", clab_trndrr(x, y, 5, 0, 1.959963984540054, st, w, iw), CLAB_OK);
  tk_eq(t, "mk.s", st[0], 10, 0);
  tk_eq(t, "mk.var", st[1], 300.0 / 18, 1e-14);
  tk_eq(t, "mk.tau", st[2], 1, 0);
  tk_eq(t, "mk.p", st[3], 2.0 / 120, 1e-15);
  tk_eq(t, "mk.slope", st[4], 1, 0);
  tk_eq(t, "mk.int", st[5], 0, 0);
  tk_ieq(t, "mk.up", clab_trndrr(x, y, 5, 1, 1.96, st, w, iw), CLAB_OK);
  tk_eq(t, "mk.pup", st[3], 1.0 / 120, 1e-15);
  tk_ieq(t, "mk.down", clab_trndrr(x, y, 5, 2, 1.96, st, w, iw), CLAB_OK);
  tk_eq(t, "mk.pdown", st[3], 1, 0);
  // ties: 1, 1, 2 has S = 2; of its six orderings four reach |S| = 2
  const double yt[3] = {1, 1, 2};
  tk_ieq(t, "mk.tie.rc", clab_trndrr(x, yt, 3, 0, 1.96, st, w, iw), CLAB_OK);
  tk_eq(t, "mk.tie.s", st[0], 2, 0);
  tk_eq(t, "mk.tie.var", st[1], (66.0 - 18.0) / 18, 1e-14);
  tk_eq(t, "mk.tie.p", st[3], 4.0 / 6, 1e-15);
  tk_true(t, "mk.tie.clamp", st[8] == 1 && st[9] == 1);
  // ten values, one adjacent swap: S = 43, the normal approximation
  const double yn[10] = {1, 2, 3, 5, 4, 6, 7, 8, 9, 10};
  tk_ieq(t, "mk.n.rc", clab_trndrr(x, yn, 10, 0, 1.96, st, w, iw), CLAB_OK);
  tk_eq(t, "mk.n.s", st[0], 43, 0);
  tk_eq(t, "mk.n.p", st[3], std::erfc(42 / std::sqrt(125.0) / std::sqrt(2.0)), 1e-14);
  tk_ieq(t, "mk.n.up", clab_trndrr(x, yn, 10, 1, 1.96, st, w, iw), CLAB_OK);
  tk_eq(t, "mk.n.pu", st[3], 0.5 * std::erfc(42 / std::sqrt(125.0) / std::sqrt(2.0)), 1e-14);
  tk_ieq(t, "mk.n.dn", clab_trndrr(x, yn, 10, 2, 1.96, st, w, iw), CLAB_OK);
  tk_eq(t, "mk.n.pd", st[3], 0.5 * std::erfc(-42 / std::sqrt(125.0) / std::sqrt(2.0)), 1e-14);
  // all tied beyond 8: p is 1
  const double yc[10] = {3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
  tk_ieq(t, "mk.c.rc", clab_trndrr(x, yc, 10, 0, 1.96, st, w, iw), CLAB_OK);
  tk_eq(t, "mk.c.p", st[3], 1, 0);
  // an even number of slopes: the median of the middle two
  const double y4[4] = {0, 2, 3, 3};
  clab_trndrr(x, y4, 4, 0, 1.96, st, w, iw);
  tk_eq(t, "mk.med", st[4], 1, 1e-15);
  tk_ieq(t, "mk.arg", clab_trndrr(x, y, 2, 0, 1.96, st, w, iw), CLAB_EARG);
  const double xb[3] = {1, 1, 2};
  tk_ieq(t, "mk.x", clab_trndrr(xb, yt, 3, 0, 1.96, st, w, iw), CLAB_EDOM);
}

void step(tk &t) {
  static double w[20];
  static int iw[20];
  double st[12];
  // 0, 0, 10, 10: the split in the middle, T = 100; 8 of the 24 orderings
  // put the tens together at one end
  const double y[4] = {0, 0, 10, 10};
  tk_ieq(t, "stp.rc", clab_stpcp(y, 4, 1, 9999, 1, st, w, iw), CLAB_OK);
  tk_eq(t, "stp.k", st[0], 2, 0);
  tk_eq(t, "stp.before", st[1], 0, 0);
  tk_eq(t, "stp.after", st[2], 10, 0);
  tk_eq(t, "stp.t", st[3], 100, 1e-14);
  tk_eq(t, "stp.p", st[4], 8.0 / 24, 1e-15);
  tk_eq(t, "stp.n", st[5], 24, 0);
  tk_eq(t, "stp.exact", st[6], 1, 0);
  // Pettitt: U_2 = 4, the largest
  tk_eq(t, "pet.t", st[7], 2, 0);
  tk_eq(t, "pet.k", st[8], 4, 0);
  tk_eq(t, "pet.p", st[9], 2 * std::exp(-1.2), 1e-15);
  tk_eq(t, "pet.before", st[10], 0, 0);
  // random reorderings: seeded, (1 + count) / (1 + nperm)
  const double z[12] = {5, 6, 5, 7, 6, 5, 1, 2, 1, 2, 1, 0};
  tk_ieq(t, "stp.mc", clab_stpcp(z, 12, 2, 199, 5, st, w, iw), CLAB_OK);
  const double p1 = st[4];
  clab_stpcp(z, 12, 2, 199, 5, st, w, iw);
  tk_eq(t, "stp.seed", st[4], p1, 0);
  tk_lt(t, "stp.mc.p", p1, 0.05);
  tk_eq(t, "stp.mc.n", st[5], 199, 0);
  tk_eq(t, "stp.mc.exact", st[6], 0, 0);
  tk_ieq(t, "stp.arg", clab_stpcp(z, 3, 2, 99, 5, st, w, iw), CLAB_EARG);
  const double zn[4] = {1, 2, std::nan(""), 3};
  tk_ieq(t, "stp.nan", clab_stpcp(zn, 4, 1, 99, 5, st, w, iw), CLAB_EDOM);
}

void cus(tk &t) {
  const double y[3] = {5, 3, 9}, m[3] = {4, 4, 4};
  double z[3], s[3];
  int al[3];
  tk_ieq(t, "cus.rc", clab_cusrv(y, m, 3, 1, 1.4, 0, 0, z, s, al), CLAB_OK);
  tk_eq(t, "cus.z", z[2], 2.5, 1e-15);
  tk_eq(t, "cus.s", s[2], 1.5, 1e-15);
  tk_true(t, "cus.alarm", al[0] == 0 && al[1] == 0 && al[2] == 1);
  const double y2[3] = {9, 9, 4};
  tk_ieq(t, "cus.reset", clab_cusrv(y2, m, 3, 1, 1.4, 0, 1, z, s, al), CLAB_OK);
  tk_true(t, "cus.reset.s", al[0] == 1 && al[1] == 1 && s[1] == 1.5);
  clab_cusrv(y, m, 3, 1, 1.4, 1, 0, z, s, al);
  tk_eq(t, "cus.rossi", z[0], (5 - 12 + 2 * std::sqrt(20.0)) / 4, 1e-15);
  clab_cusrv(y, m, 3, 1, 1.4, 2, 0, z, s, al);
  tk_eq(t, "cus.ansc", z[2], 1.5 * (std::pow(81.0, 1.0 / 3) - std::pow(16.0, 1.0 / 3)) / std::pow(4.0, 1.0 / 6),
        1e-14);
  clab_cusrv(y, m, 3, 1, 1.4, 3, 0, z, s, al);
  tk_eq(t, "cus.none", z[1], 3, 0);
  tk_ieq(t, "cus.arg", clab_cusrv(y, m, 3, 1, 0, 0, 0, z, s, al), CLAB_EARG);
  const double mb[3] = {4, 0, 4};
  tk_ieq(t, "cus.dom", clab_cusrv(y, mb, 3, 1, 1, 0, 0, z, s, al), CLAB_EDOM);
}

void xc(tk &t) {
  const double x[4] = {1, 2, 3, 4}, y[4] = {2, 4, 6, 8};
  double r[3];
  tk_ieq(t, "xc.rc", clab_xcser(x, y, 4, 1, r), CLAB_OK);
  tk_eq(t, "xc.0", r[1], 1, 1e-15);
  tk_eq(t, "xc.p1", r[2], 0.25, 1e-15);
  tk_eq(t, "xc.m1", r[0], 0.25, 1e-15);
  const double c[4] = {1, 1, 1, 1};
  tk_ieq(t, "xc.const", clab_xcser(x, c, 4, 1, r), CLAB_EDOM);
  const double xn[4] = {1, 2, HUGE_VAL, 4};
  tk_ieq(t, "xc.inf", clab_xcser(xn, y, 4, 1, r), CLAB_EDOM);
  tk_ieq(t, "xc.lag", clab_xcser(x, y, 4, 4, r), CLAB_EARG);
}

}  // namespace

void t_trn(tk &t, int big) {
  (void)big;
  mk(t);
  step(t);
  cus(t);
  xc(t);
}
