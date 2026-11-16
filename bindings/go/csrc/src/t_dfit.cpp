// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/dfit.cpp: the closed forms against their definitions, the
// minimisers against them where both apply (the gamma with shape 1 is the
// exponential), and the statistics of a perfect fit.
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "tchk.h"

void t_dfit(tk &T, int) {
  const double x[8] = {0.5, 1.2, 1.9, 2.4, 3.1, 4.0, 5.2, 7.7};
  double p[3], q[3];
  tk_ieq(T, "dfit.np", clab_dfnp(5), 3);
  clab_dfit(0, x, 8, p);
  tk_eq(T, "dfit.mean", p[0], 3.25, 1e-15);
  clab_dfit(2, x, 8, q);
  tk_eq(T, "dfit.rate", q[0], 1 / 3.25, 1e-15);
  // the gamma and the Weibull nest the exponential: their likelihood is at
  // least the exponential's
  tk_ieq(T, "dfit.gam", clab_dfit(3, x, 8, p), CLAB_OK);
  tk_ge(T, "dfit.gll", clab_dfll(3, x, 8, p), clab_dfll(2, x, 8, q) - 1e-9);
  tk_ieq(T, "dfit.wei", clab_dfit(4, x, 8, p), CLAB_OK);
  tk_ge(T, "dfit.wll", clab_dfll(4, x, 8, p), clab_dfll(2, x, 8, q) - 1e-9);
  // the gamma's ML mean is the sample mean
  clab_dfit(3, x, 8, p);
  tk_eq(T, "dfit.gmean", p[0] / p[1], 3.25, 1e-6);
  const double c[10] = {0, 0, 1, 1, 2, 3, 0, 4, 2, 1};
  tk_ieq(T, "dfit.nb", clab_dfit(9, c, 10, p), CLAB_OK);
  tk_eq(T, "dfit.nbmu", p[1], 1.4, 1e-15);
  tk_ieq(T, "dfit.zip", clab_dfit(10, c, 10, p), CLAB_OK);
  tk_true(T, "dfit.zippi", p[0] >= 0 && p[0] < 1);
  double g[2];
  clab_dfit(0, x, 8, p);
  clab_dfgof(0, x, 8, p, g);
  tk_true(T, "dfit.ks", g[0] > 0 && g[0] < 0.5 && std::isfinite(g[1]));
  clab_dfgof(8, c, 10, q, g);
  tk_true(T, "dfit.disc", std::isnan(g[1]));
  // the discrete KS over every whole number up to the largest value
  double ks = 0;
  for (int k = 0; k <= 4; ++k) {
    double m = 0;
    for (int i = 0; i < 10; ++i) m += c[i] <= k;
    ks = std::fmax(ks, std::fabs(m / 10 - clab_ppois(k, q[0], 1)));
  }
  tk_eq(T, "dfit.dks", g[0], ks, 0);
  // values far past 2^53 take no longer than small ones
  const double h[3] = {0, 3e50, 1e51}, hp[1] = {4e50};
  clab_dfgof(8, h, 3, hp, g);
  tk_true(T, "dfit.big", g[0] >= 0 && g[0] <= 1);
}
