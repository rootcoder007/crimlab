// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/xpl.cpp: designs and summaries worked out by hand on a few
// rows. Built with the self-test.
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "tchk.h"

namespace {

void imp(tk &t) {
  // x = (1, 2, 3) and (5, 5, 5); f = column 1, y = f: column 2 permuted
  // changes nothing, column 1 by its permutation
  const double x[6] = {1, 2, 3, 5, 5, 5}, y[3] = {1, 2, 3};
  const int n = 3, p = 2, k = 2;
  std::vector<double> des(n * p * k * p), fp(n * p * k);
  std::vector<int> pm(n * p * k);
  tk_ieq(t, "xp.perm.rc", clab_xpperm(x, n, p, k, 7, des.data(), pm.data()), CLAB_OK);
  int ok = 1;
  for (int b = 0; b < p * k; ++b) {
    int s = 0;
    for (int i = 0; i < n; ++i) s += 1 << pm[b * n + i];
    ok = ok && s == 7;
  }
  tk_ieq(t, "xp.perm.isperm", ok, 1);
  const int N = n * p * k;
  for (int r = 0; r < N; ++r) fp[r] = des[r];
  double rep[4], st[5];
  tk_ieq(t, "xp.imp.rc", clab_xpimp(y, x, fp.data(), n, p, k, 0, 0, rep, st), CLAB_OK);
  double want = 0;
  for (int i = 0; i < n; ++i) want += (y[i] - x[pm[i]]) * (y[i] - x[pm[i]]) / n;
  tk_near(t, "xp.imp.j0", rep[0], want, 1e-15);
  tk_near(t, "xp.imp.j1", rep[1], 0, 0);
  tk_near(t, "xp.imp.base", st[0], 0, 0);
  tk_ieq(t, "xp.imp.ratio0", clab_xpimp(y, x, fp.data(), n, p, k, 0, 1, rep, st), CLAB_EDOM);
  tk_ieq(t, "xp.imp.loss", clab_xpimp(y, x, fp.data(), n, p, k, 2, 0, rep, st), CLAB_EARG);
}

void pdp(tk &t) {
  const double v[5] = {4, 0, 1, 3, 2};
  double g[5];
  int m = 0;
  // seq(0, 4, length.out = 5) and the quantiles of 0..4 at 0, 1/4, ..., 1
  tk_ieq(t, "xp.grid.rc", clab_xpgrid(v, 5, 5, 0, g, &m), CLAB_OK);
  tk_ieq(t, "xp.grid.m", m, 5);
  tk_near(t, "xp.grid.e2", g[2], 2, 0);
  tk_ieq(t, "xp.grid.q", clab_xpgrid(v, 5, 3, 1, g, &m), CLAB_OK);
  tk_near(t, "xp.grid.q1", g[1], 2, 0);
  // quantile(c(0, 0, 0, 1), c(0, 0.5, 1)) is 0, 0, 1: one repeat dropped
  const double w[4] = {0, 0, 0, 1};
  tk_ieq(t, "xp.grid.u", clab_xpgrid(w, 4, 3, 1, g, &m), CLAB_OK);
  tk_ieq(t, "xp.grid.um", m, 2);
  const double s[4] = {1, 2, 10, 20}, gg[2] = {-1, 7};
  double des[16], pd[2];
  tk_ieq(t, "xp.cpd.rc", clab_xpcpd(s, 2, 2, 0, gg, 2, des), CLAB_OK);
  tk_near(t, "xp.cpd.a", des[2], 7, 0);
  tk_near(t, "xp.cpd.b", des[4 + 3], 20, 0);
  const double f[4] = {1, 3, 10, 20};
  tk_ieq(t, "xp.pd.rc", clab_xppd(f, 2, 2, pd), CLAB_OK);
  tk_near(t, "xp.pd.1", pd[1], 15, 0);
}

void ale(tk &t) {
  // x = 1..4, K = 2: edges 1 and the type-1 quantiles at 1/2 and 1, 2 and
  // 4; the values 1 and 2 in [1, 2], 3 and 4 in (2, 4]
  const double x[4] = {3, 1, 4, 2};
  double z[3], des[8];
  int nz = 0, bin[4];
  tk_ieq(t, "xp.aled.rc", clab_xpaled(x, 4, 1, 0, 2, z, &nz, bin, des), CLAB_OK);
  tk_ieq(t, "xp.aled.nz", nz, 3);
  tk_near(t, "xp.aled.z1", z[1], 2, 0);
  tk_ieq(t, "xp.aled.bin", bin[0] * 1000 + bin[1] * 100 + bin[2] * 10 + bin[3], 1010);
  // f = x^2: effects (4 - 1) and (16 - 4); ale 0, 3, 15 less the mean of
  // the midpoints 1.5, 1.5, 9, 9
  double f[8], eff[2], cnt[2], a[3];
  for (int i = 0; i < 8; ++i) f[i] = des[i] * des[i];
  tk_ieq(t, "xp.ale.rc", clab_xpale(f, 4, bin, 2, eff, cnt, a), CLAB_OK);
  tk_near(t, "xp.ale.e1", eff[1], 12, 0);
  tk_near(t, "xp.ale.a0", a[0], -5.25, 1e-15);
  tk_near(t, "xp.ale.a2", a[2], 9.75, 1e-15);
  const double c[3] = {2, 2, 2};
  tk_ieq(t, "xp.aled.const", clab_xpaled(c, 3, 1, 0, 2, z, &nz, bin, des), CLAB_EDOM);
}

void shap(tk &t) {
  // one reference row w = (0, 0), instance (1, 2), f = 3 x1 + x2: every
  // sample gives 3 for x1 and 2 for x2
  const double x[2] = {0, 0}, s[2] = {1, 2};
  const int m = 5;
  std::vector<double> des(2 * 2 * m * 2), f(2 * 2 * m);
  tk_ieq(t, "xp.shd.rc", clab_xpshd(x, 1, s, 1, 2, m, 3, des.data()), CLAB_OK);
  const int N = 2 * 2 * m;
  for (int r = 0; r < N; ++r) f[r] = 3 * des[r] + des[N + r];
  double phi[2], se[2];
  tk_ieq(t, "xp.sh.rc", clab_xpsh(f.data(), 1, 2, m, phi, se), CLAB_OK);
  tk_near(t, "xp.sh.p0", phi[0], 3, 0);
  tk_near(t, "xp.sh.p1", phi[1], 2, 0);
  tk_near(t, "xp.sh.se", se[0], 0, 0);
}

}  // namespace

void t_xpl(tk &t, int big) {
  (void)big;
  imp(t);
  pdp(t);
  ale(t);
  shap(t);
}
