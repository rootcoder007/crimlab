// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for the dynamic and heterogeneous panel kernels
// (dpn.cpp): GMM in a just-identified case, where it is instrumental
// variables and both steps agree, and in an overidentified case with one
// instrument per block of rows, where the weights are diagonal and every
// estimate is a ratio by hand; Wooldridge's statistic of three short
// series by hand; the common correlated effects of an exact relation; and
// the failure paths. Part of the self-test (selftest.cpp).
#include <cmath>

#include "crimlab.h"
#include "tchk.h"

namespace {

// clab_gmmpnl without differenced rows (no Arellano-Bond tests)
int g0(int n, int k, int m, int ng, const int *go, const double *y, const double *x,
       const double *z, const double *h, const int *eq, const int *dgo, const int *sub, int s,
       int rob, double *b, double *v, double *st) {
  double ar[1];
  return clab_gmmpnl(n, k, m, ng, go, y, x, z, h, eq, 0, dgo, nullptr, nullptr, nullptr, nullptr,
                     nullptr, sub, s, rob, 0, 1.0, b, v, st, ar);
}

void gmmjust(tk &tt) {
  // one row per group, H = 1: b = z'y / z'x = 5 / 3, e = (-2, 1, 2) / 3,
  // the scores z e = (-2, 0, 2) / 3, V = (8 / 9) / 3^2 at both steps and no
  // correction, since Z'e = 0; the conventional one-step variance is
  // sigma^2 (z'x (z'z)^-1 z'x)^-1 = (1 / 3) / (9 / 2), sigma^2 = e'e / 3
  const int go[4] = {0, 1, 2, 3}, eq[3] = {1, 1, 1}, dgo[4] = {0, 0, 0, 0};
  const double y[3] = {1, 2, 4}, x[3] = {1, 1, 2}, z[3] = {1, 0, 1}, h[3] = {1, 1, 1};
  double b[1], v[1], st[6];
  for (int s = 1; s <= 2; ++s) {
    tk_ieq(tt, "dpn.just.rc", g0(3, 1, 1, 3, go, y, x, z, h, eq, dgo, nullptr, s, 1, b, v, st),
           CLAB_OK);
    tk_near(tt, "dpn.just.b", b[0], 5.0 / 3, 1e-14);
    tk_near(tt, "dpn.just.v", v[0], 8.0 / 81, 1e-14);
    tk_near(tt, "dpn.just.j", st[1], 0, 1e-13);
    tk_near(tt, "dpn.just.s", st[0], 0, 1e-13);
    g0(3, 1, 1, 3, go, y, x, z, h, eq, dgo, nullptr, s, 0, b, v, st);
    tk_near(tt, "dpn.just.vn", v[0], s == 1 ? 2.0 / 27 : 8.0 / 81, 1e-14);
    tk_near(tt, "dpn.just.sig", st[3], 1.0 / 3, 1e-14);
  }
  tk_ieq(tt, "dpn.just.k", g0(3, 2, 1, 3, go, y, x, z, h, eq, dgo, nullptr, 1, 1, b, v, st),
         CLAB_EARG);
  tk_ieq(tt, "dpn.just.steps", g0(3, 1, 1, 3, go, y, x, z, h, eq, dgo, nullptr, 3, 1, b, v, st),
         CLAB_EARG);
  const int bad[4] = {0, 2, 1, 3};
  tk_ieq(tt, "dpn.just.go", g0(3, 1, 1, 3, bad, y, x, z, h, eq, dgo, nullptr, 1, 1, b, v, st),
         CLAB_EARG);
  const double z0[3] = {0, 0, 0}, xn[3] = {1, NAN, 2};
  tk_ieq(tt, "dpn.just.zero", g0(3, 1, 1, 3, go, y, x, z0, h, eq, dgo, nullptr, 1, 1, b, v, st),
         CLAB_EDOM);
  tk_ieq(tt, "dpn.just.nan", g0(3, 1, 1, 3, go, y, xn, z, h, eq, dgo, nullptr, 1, 1, b, v, st),
         CLAB_EDOM);
}

void gmmover(tk &tt) {
  // two instruments, each on two of the four one-row groups: Z'Z =
  // diag(2, 2), Z'x = (4, 3), Z'y = (9, 4)
  const int go[5] = {0, 1, 2, 3, 4}, eq[4] = {1, 1, 1, 1}, dgo[5] = {0, 0, 0, 0, 0};
  const double y[4] = {2, 3, 7, 1}, x[4] = {1, 2, 3, 1}, z[8] = {1, 0, 1, 0, 0, 1, 0, 1},
               h[4] = {1, 1, 1, 1};
  double b[1], v[1], st[6];
  tk_ieq(tt, "dpn.over.rc1", g0(4, 1, 2, 4, go, y, x, z, h, eq, dgo, nullptr, 1, 1, b, v, st),
         CLAB_OK);
  const double b1 = (4.0 * 9 / 2 + 3.0 * 4 / 2) / (16.0 / 2 + 9.0 / 2);
  tk_near(tt, "dpn.over.b1", b[0], b1, 1e-14);
  // the moment variance is diagonal: the sums of e1^2 over each block
  double s1 = 0, s2 = 0, ss = 0;
  for (int i = 0; i < 4; ++i) {
    const double e = y[i] - b1 * x[i];
    (i % 2 ? s2 : s1) += e * e;
    ss += e * e;
  }
  const double m1 = 16.0 / 2 + 9.0 / 2;
  const double v1 = (16.0 / 4 * s1 + 9.0 / 4 * s2) / (m1 * m1);
  tk_near(tt, "dpn.over.v1", v[0], v1, 1e-13);
  // Sargan: (Z'e1)' (Z'Z)^-1 (Z'e1) over sigma^2 = e1'e1 / 4
  const double z1 = 9 - 4 * b1, z2 = 4 - 3 * b1;
  tk_near(tt, "dpn.over.sargan", st[0], (z1 * z1 / 2 + z2 * z2 / 2) / (ss / 4), 1e-12);
  const double b2 = (4.0 * 9 / s1 + 3.0 * 4 / s2) / (16.0 / s1 + 9.0 / s2);
  const double j = (9 - 4 * b2) * (9 - 4 * b2) / s1 + (4 - 3 * b2) * (4 - 3 * b2) / s2;
  tk_near(tt, "dpn.over.hansen", st[1], j, 1e-12);
  tk_near(tt, "dpn.over.sigma", st[3], ss / 4, 1e-14);
  g0(4, 1, 2, 4, go, y, x, z, h, eq, dgo, nullptr, 1, 0, b, v, st);
  tk_near(tt, "dpn.over.v1n", v[0], ss / 4 / m1, 1e-14);
  g0(4, 1, 2, 4, go, y, x, z, h, eq, dgo, nullptr, 2, 0, b, v, st);
  tk_near(tt, "dpn.over.b2", b[0], b2, 1e-13);
  tk_near(tt, "dpn.over.v2", v[0], 1 / (16.0 / s1 + 9.0 / s2), 1e-13);
  tk_true(tt, "dpn.over.gen", st[4] == 0);
  // the corrected variance differs from the uncorrected one here
  double vc[1];
  g0(4, 1, 2, 4, go, y, x, z, h, eq, dgo, nullptr, 2, 1, b, vc, st);
  tk_true(tt, "dpn.over.wc", std::fabs(vc[0] - v[0]) > 1e-8);
  // leaving the second instrument out leaves one: both statistics without
  // it are 0
  const int sub[2] = {0, 1};
  g0(4, 1, 2, 4, go, y, x, z, h, eq, dgo, sub, 2, 1, b, vc, st);
  tk_near(tt, "dpn.over.hx", st[2], 0, 1e-12);
  tk_near(tt, "dpn.over.sx", st[5], 0, 1e-12);
}

void wsct(tk &tt) {
  // three units, four periods, no x: the differences e are (1, 2, -1),
  // (0, 1, 2), (2, -1, 0); sum e e_1 = 0 and sum e_1^2 = 11, so the
  // coefficient is 0; the cluster scores are 0, 2, -2, so V = 8 / 121
  // times 3 / 2 (n - 1) / (n - 1)
  const double y[12] = {0, 1, 3, 2, 1, 1, 2, 4, 0, 2, 1, 1};
  const int u[12] = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2},
            t[12] = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
  double o[6];
  tk_ieq(tt, "dpn.wsct.rc", clab_wsctst(12, 0, y, nullptr, u, t, 0, 0, nullptr, nullptr, o),
         CLAB_OK);
  tk_near(tt, "dpn.wsct.rho", o[0], 0, 1e-15);
  tk_near(tt, "dpn.wsct.se", o[1], std::sqrt(12.0) / 11, 1e-14);
  tk_near(tt, "dpn.wsct.f", o[2], 121.0 / 48, 1e-13);
  tk_near(tt, "dpn.wsct.g", o[3], 3, 0);
  tk_near(tt, "dpn.wsct.np", o[4], 6, 0);
  tk_near(tt, "dpn.wsct.nf", o[5], 9, 0);
  clab_wsctst(12, 0, y, nullptr, u, t, 0, 1, nullptr, nullptr, o);
  tk_near(tt, "dpn.wsct.raw", o[2], 121.0 / 32, 1e-13);
  tk_ieq(tt, "dpn.wsct.n", clab_wsctst(1, 0, y, nullptr, u, t, 0, 0, nullptr, nullptr, o),
         CLAB_EARG);
  // one unit: one cluster
  tk_ieq(tt, "dpn.wsct.one", clab_wsctst(4, 0, y, nullptr, u, t, 0, 0, nullptr, nullptr, o),
         CLAB_EDOM);
}

void cce(tk &tt) {
  // y = 2 x exactly: every unit estimate is 2 whatever the averages, and
  // the spreads are 0
  const int nu = 3, nt = 6;
  double x[18], y[18], w[3] = {1.0 / 3, 1.0 / 3, 1.0 / 3}, d[6] = {1, 1, 1, 1, 1, 1};
  for (int i = 0; i < 18; ++i) {
    x[i] = std::sin(1.0 + i * i);
    y[i] = 2 * x[i];
  }
  double bi[3], bmg[1], vmg[1], bp[1], vp[1];
  int bad = 0;
  const int av[2] = {1, 1}, nav[2] = {0, 0};
  tk_ieq(tt, "dpn.cce.rc", clab_ccepnl(nu, nt, 1, y, x, w, d, 1, av, bi, bmg, vmg, bp, vp, &bad),
         CLAB_OK);
  tk_near(tt, "dpn.cce.bi", bi[2], 2, 1e-12);
  tk_near(tt, "dpn.cce.bmg", bmg[0], 2, 1e-12);
  tk_near(tt, "dpn.cce.bp", bp[0], 2, 1e-12);
  tk_near(tt, "dpn.cce.vmg", vmg[0], 0, 1e-20);
  tk_near(tt, "dpn.cce.vp", vp[0], 0, 1e-20);
  // x the same in every unit is x's average: nothing is left of it
  double xs[18];
  for (int i = 0; i < 18; ++i) xs[i] = x[i % 6];
  tk_ieq(tt, "dpn.cce.flat", clab_ccepnl(nu, nt, 1, y, xs, w, d, 1, av, bi, bmg, vmg, bp, vp, &bad),
         CLAB_EDOM);
  tk_ieq(tt, "dpn.cce.badu", bad, 0);
  // unit 2 not seen in periods 0 and 5: still y = 2 x on its four periods;
  // period 0 then has two units, and without unit 1 there it has one
  double yu[18];
  for (int i = 0; i < 18; ++i) yu[i] = y[i];
  yu[12] = yu[17] = NAN;
  tk_ieq(tt, "dpn.cce.unb", clab_ccepnl(nu, nt, 1, yu, x, w, d, 1, av, bi, bmg, vmg, bp, vp, &bad),
         CLAB_OK);
  tk_near(tt, "dpn.cce.unb2", bi[2], 2, 1e-12);
  tk_near(tt, "dpn.cce.unbp", bp[0], 2, 1e-12);
  yu[6] = NAN;
  tk_ieq(tt, "dpn.cce.per", clab_ccepnl(nu, nt, 1, yu, x, w, d, 1, av, bi, bmg, vmg, bp, vp, &bad),
         CLAB_EDOM);
  tk_ieq(tt, "dpn.cce.badp", bad, -2);
  const double w2[3] = {0.5, 0.5, 0.5};
  tk_ieq(tt, "dpn.cce.w", clab_ccepnl(nu, nt, 1, y, x, w2, d, 1, av, bi, bmg, vmg, bp, vp, &bad),
         CLAB_EDOM);
  // no average and no factor: nothing to take out
  tk_ieq(tt, "dpn.cce.none", clab_ccepnl(nu, nt, 1, y, x, w, d, 0, nav, bi, bmg, vmg, bp, vp, &bad),
         CLAB_EARG);
  // the average of x alone still leaves y = 2 x exact
  const int ax[2] = {0, 1};
  clab_ccepnl(nu, nt, 1, y, x, w, d, 1, ax, bi, bmg, vmg, bp, vp, &bad);
  tk_near(tt, "dpn.cce.avx", bmg[0], 2, 1e-12);
  tk_ieq(tt, "dpn.cce.n", clab_ccepnl(1, nt, 1, y, x, w, d, 1, av, bi, bmg, vmg, bp, vp, &bad),
         CLAB_EARG);
}

}  // namespace

void t_dpn(tk &t, int big) {
  (void)big;
  gmmjust(t);
  gmmover(t);
  wsct(t);
  cce(t);
}
