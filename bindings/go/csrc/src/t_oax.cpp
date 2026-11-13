// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for clab_oaxbd (src/oax.cpp). Part of the self-test
// (selftest.cpp).

#include <cmath>

#include "crimlab.h"
#include "tchk.h"

namespace {

void run(tk &tt) {
  // 1. exact fits by hand: group A (code 0) y = 1 + 2 x at x = 0..3, group
  // B (code 1, the reference) y = 3 + x / 2 at x = 1..5. Means: x 1.5 and
  // 3, y 4 and 4.5. D = (0, -1.5), Delta = (-2, 1.5).
  {
    const double x[9] = {0, 1, 2, 3, 1, 2, 3, 4, 5};
    double y[9];
    int g[9];
    for (int i = 0; i < 9; ++i) {
      g[i] = i < 4 ? 0 : 1;
      y[i] = i < 4 ? 1 + 2 * x[i] : 3 + 0.5 * x[i];
    }
    double xm[4], b[4], bs[4], dec[14], var[10], wk[2 * (9 + 5 + 6) + 9];
    tk_true(tt, "oaxbd.rc", clab_oaxbd(y, x, g, 9, 1, 1, xm, b, bs, dec, var, wk) == CLAB_OK);
    const double bw[4] = {1, 2, 3, 0.5};
    for (int k = 0; k < 4; ++k) tk_near(tt, "oaxbd.b", b[k], bw[k], 1e-14);
    for (int k = 0; k < 4; ++k) tk_near(tt, "oaxbd.bse", bs[k], 0, 1e-14);
    tk_near(tt, "oaxbd.xa", xm[1], 1.5, 1e-15);
    tk_near(tt, "oaxbd.gap", dec[2], -0.5, 1e-14);
    // explained at B's -0.75, unexplained 0.25; at A's -3 and 2.5;
    // interaction -2.25
    const double want[5] = {-0.75, 0.25, -3, 2.5, -2.25};
    for (int k = 0; k < 5; ++k) tk_near(tt, "oaxbd.dec", dec[4 + 2 * k], want[k], 1e-13);
    tk_near(tt, "oaxbd.var.shift", var[4 * 2 + 1], 1.5, 1e-14);
    tk_near(tt, "oaxbd.var.exp", var[1], -1.5 * 0.5, 1e-14);
    // the other reference swaps the groups
    double d2[14];
    clab_oaxbd(y, x, g, 9, 1, 0, xm, b, bs, d2, var, wk);
    tk_near(tt, "oaxbd.swap.gap", d2[2], 0.5, 1e-14);
    tk_near(tt, "oaxbd.swap.exp", d2[4], 3, 1e-13);
  }
  // 2. random data: the three decompositions add up to the gap
  {
    static double x[600], y[300], wk[3 * (300 + 16) + 300];
    static int g[300];
    rng r = {2024};
    for (int i = 0; i < 300; ++i) {
      g[i] = runi(r) < 0.4 ? 1 : 0;
      x[i] = runi(r, 0, 4) + g[i];
      x[300 + i] = rpois(r, 2);
      y[i] = 1 + 0.5 * x[i] - 0.3 * x[300 + i] + g[i] * (0.4 + 0.2 * x[i]) + runi(r, -1, 1);
    }
    double xm[6], b[6], bs[6], dec[14], var[15];
    clab_oaxbd(y, x, g, 300, 2, 1, xm, b, bs, dec, var, wk);
    tk_near(tt, "oaxbd.two.b", dec[2], dec[4] + dec[6], 1e-13);
    tk_near(tt, "oaxbd.two.a", dec[2], dec[8] + dec[10], 1e-13);
    tk_near(tt, "oaxbd.three", dec[2], dec[4] + dec[10] + dec[12], 1e-13);
    tk_near(tt, "oaxbd.ref", dec[8] - dec[4], dec[12], 1e-13);
    double se = 0, su = 0;
    for (int j = 0; j < 3; ++j) {
      se += var[j];
      su += var[6 + j];
    }
    tk_near(tt, "oaxbd.by.exp", se, dec[4], 1e-13);
    tk_near(tt, "oaxbd.by.unexp", su, dec[6], 1e-13);
    tk_true(tt, "oaxbd.se", dec[3] > 0 && dec[5] > 0 && dec[13] > 0);
    // collinear covariates, and too few cases
    for (int i = 0; i < 300; ++i) x[300 + i] = 2 * x[i];
    tk_true(tt, "oaxbd.guard.col", clab_oaxbd(y, x, g, 300, 2, 1, xm, b, bs, dec, var, wk) == CLAB_EDOM);
    tk_true(tt, "oaxbd.guard.n", clab_oaxbd(y, x, g, 4, 2, 1, xm, b, bs, dec, var, wk) == CLAB_EDOM);
    tk_true(tt, "oaxbd.guard.ref", clab_oaxbd(y, x, g, 300, 2, 2, xm, b, bs, dec, var, wk) == CLAB_EARG);
  }
}

}  // namespace

void t_oax(tk &t, int) { run(t); }
