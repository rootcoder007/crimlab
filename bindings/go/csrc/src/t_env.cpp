// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for the environmental-crime kernels (env.cpp): the
// LMDI contributions of two-factor identities by hand, their sum, the
// variance inflation factor of two columns from their correlation, the
// Jarque-Bera moments of small samples worked by hand, Pesaran's CD of
// three series by hand, and the failure paths. Part of the self-test
// (selftest.cpp).
#include <cmath>

#include "crimlab.h"
#include "tchk.h"

namespace {

void lmdi(tk &tt) {
  // 6 = 2 x 3 to 12 = 3 x 4: L = 6 / ln 2, the shares ln 1.5 and ln 4/3
  {
    const double y0[1] = {6}, y1[1] = {12}, f0[2] = {2, 3}, f1[2] = {3, 4};
    double c[2];
    tk_ieq(tt, "env.lmdi.rc", clab_lmdi(y0, y1, f0, f1, 1, 2, c), CLAB_OK);
    const double l = 6 / std::log(2.0);
    tk_near(tt, "env.lmdi.c1", c[0], l * std::log(1.5), 1e-14);
    tk_near(tt, "env.lmdi.c2", c[1], l * std::log(4.0 / 3), 1e-14);
    tk_near(tt, "env.lmdi.sum", c[0] + c[1], 6, 1e-14);
  }
  // no change in y: L is y and the shares cancel
  {
    const double y[1] = {5}, f0[2] = {5, 1}, f1[2] = {1, 5};
    double c[2];
    clab_lmdi(y, y, f0, f1, 1, 2, c);
    tk_near(tt, "env.lmdi.same1", c[0], -5 * std::log(5.0), 1e-14);
    tk_near(tt, "env.lmdi.same2", c[1], 5 * std::log(5.0), 1e-14);
  }
  // two rows, factors by column: 10 = 2 x 5 to 20 = 4 x 5 puts all of 10
  // on the first factor; 1 x 1 to 1 x 1 nothing
  {
    const double y0[2] = {10, 1}, y1[2] = {20, 1}, f0[4] = {2, 1, 5, 1},
                 f1[4] = {4, 1, 5, 1};
    double c[4];
    clab_lmdi(y0, y1, f0, f1, 2, 2, c);
    tk_near(tt, "env.lmdi.col1", c[0], 10, 1e-13);
    tk_near(tt, "env.lmdi.col2", c[2], 0, 1e-15);
    tk_near(tt, "env.lmdi.row2", c[1] + c[3], 0, 1e-15);
  }
  // close values: the log mean of 1 and 1 + 1e-12 is about 1 + 5e-13
  {
    const double y0[1] = {1}, y1[1] = {1 + 1e-12}, f0[1] = {1}, f1[1] = {1 + 1e-12};
    double c[1];
    clab_lmdi(y0, y1, f0, f1, 1, 1, c);
    tk_near(tt, "env.lmdi.close", c[0], 1e-12, 1e-9);
  }
  {
    const double y0[1] = {0}, y1[1] = {1}, f[1] = {1};
    double c[1];
    tk_ieq(tt, "env.lmdi.zero", clab_lmdi(y0, y1, f, f, 1, 1, c), CLAB_EDOM);
    tk_ieq(tt, "env.lmdi.n", clab_lmdi(y1, y1, f, f, 0, 1, c), CLAB_EARG);
  }
}

void vif(tk &tt) {
  // 1..5 and (1, 3, 2, 5, 4) correlate 0.8: both VIFs are 1 / 0.36
  double x[15] = {1, 2, 3, 4, 5, 1, 3, 2, 5, 4, 0, 0, 0, 0, 0}, o[3];
  double wk[5 * 4 + 2 * 3 * 4];
  tk_ieq(tt, "env.vif.rc", clab_vif(x, 5, 2, o, wk), CLAB_OK);
  tk_near(tt, "env.vif.a", o[0], 25.0 / 9, 1e-12);
  tk_near(tt, "env.vif.b", o[1], 25.0 / 9, 1e-12);
  clab_vif(x, 5, 1, o, wk);
  tk_near(tt, "env.vif.one", o[0], 1, 0);
  // a constant column has none; the sum of the two others is exact
  clab_vif(x, 5, 3, o, wk);
  tk_nan(tt, "env.vif.const", o[2]);
  for (int i = 0; i < 5; ++i) x[10 + i] = x[i] + x[5 + i];
  clab_vif(x, 5, 3, o, wk);
  tk_true(tt, "env.vif.exact", std::isinf(o[2]) || o[2] > 1e10);
  tk_ieq(tt, "env.vif.n", clab_vif(x, 3, 2, o, wk), CLAB_EARG);
}

void jb(tk &tt) {
  // (-1, 0, 0, 1): S 0, K 2, JB 4 / 6 x 1 / 4
  {
    const double r[4] = {-1, 0, 0, 1};
    double o[4];
    tk_ieq(tt, "env.jb.rc", clab_jb(r, 4, o), CLAB_OK);
    tk_near(tt, "env.jb.s", o[0], 0, 1e-15);
    tk_near(tt, "env.jb.k", o[1], 2, 1e-14);
    tk_near(tt, "env.jb.jb", o[2], 1.0 / 6, 1e-14);
    tk_near(tt, "env.jb.p", o[3], std::exp(-1.0 / 12), 1e-14);
  }
  // (0, 0, 0, 1): S 2 / sqrt(3), K 7 / 3, JB 26 / 27
  {
    const double r[4] = {0, 0, 0, 1};
    double o[4];
    clab_jb(r, 4, o);
    tk_near(tt, "env.jb.s2", o[0], 2 / std::sqrt(3.0), 1e-13);
    tk_near(tt, "env.jb.k2", o[1], 7.0 / 3, 1e-13);
    tk_near(tt, "env.jb.jb2", o[2], 26.0 / 27, 1e-13);
  }
  const double c[4] = {2, 2, 2, 2};
  double o[4];
  tk_ieq(tt, "env.jb.const", clab_jb(c, 4, o), CLAB_EDOM);
  tk_ieq(tt, "env.jb.n", clab_jb(c, 3, o), CLAB_EARG);
}

void pcd(tk &tt) {
  // (1, 2, 3), (1, 3, 2), (3, 2, 1): rho 0.5, -1, -0.5; CD -1
  double e[9] = {1, 2, 3, 1, 3, 2, 3, 2, 1}, o[5];
  tk_ieq(tt, "env.pcd.rc", clab_pcd(e, 3, 3, o), CLAB_OK);
  tk_near(tt, "env.pcd.cd", o[0], -1, 1e-14);
  tk_near(tt, "env.pcd.p", o[1], 0.31731050786291415, 1e-14);
  tk_near(tt, "env.pcd.rho", o[2], -1.0 / 3, 1e-14);
  tk_near(tt, "env.pcd.abs", o[3], 2.0 / 3, 1e-14);
  tk_near(tt, "env.pcd.m", o[4], 3, 0);
  // the third unit misses its last period: its pairs are over two
  // periods, both rho -1
  e[8] = NAN;
  clab_pcd(e, 3, 3, o);
  tk_near(tt, "env.pcd.unb", o[0],
          (std::sqrt(3.0) * 0.5 - 2 * std::sqrt(2.0)) / std::sqrt(3.0), 1e-14);
  // one period in common, or a constant series: no pair
  const double f[4] = {1, NAN, NAN, 2};
  tk_ieq(tt, "env.pcd.none", clab_pcd(f, 2, 2, o), CLAB_EDOM);
  const double g[4] = {1, 1, 1, 2};
  tk_ieq(tt, "env.pcd.flat", clab_pcd(g, 2, 2, o), CLAB_EDOM);
  tk_ieq(tt, "env.pcd.n", clab_pcd(g, 1, 2, o), CLAB_EARG);
}

}  // namespace

void t_env(tk &t, int) {
  lmdi(t);
  vif(t);
  jb(t);
  pcd(t);
}
