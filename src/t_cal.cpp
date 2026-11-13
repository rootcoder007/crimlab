// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for clab_calaud and clab_scrdsp (src/cal.cpp). Part of
// the self-test (selftest.cpp). References: by hand, and R 4.6.1
// (cor.test(), the C routine behind it for the tails of Spearman's S,
// rank()).

#include <cmath>

#include "crimlab.h"
#include "tchk.h"

namespace {

void run(tk &tt) {
  // 1. eight areas, one record each, groups X X Y X Y Y X Y
  {
    const int a[8] = {0, 1, 2, 3, 4, 5, 6, 7}, g[8] = {0, 0, 1, 0, 1, 1, 0, 1};
    const double r[8] = {90, 80, 70, 60, 30, 20, 10, 5}, o[8] = {10, 25, 20, 80, 70, 90, 85, 60};
    double ar[56], gr[6], st[7], wk[8 + 6];
    int iw[8 + 16];
    tk_true(tt, "calaud.rc",
            clab_calaud(a, r, o, g, nullptr, 8, 8, 2, 0, 0, ar, gr, st, iw, wk) == CLAB_OK);
    // rank(-o) - rank(-r)
    const double gap[8] = {7, 4, 4, -1, -1, -5, -5, -3};
    for (int k = 0; k < 8; ++k) tk_near(tt, "calaud.gap", ar[48 + k], gap[k], 0);
    // cor.test(r, o, method = "spearman")
    tk_near(tt, "calaud.rho", st[1], -0.69047619047619069, 1e-15);
    tk_near(tt, "calaud.s", st[2], 142, 1e-13);
    tk_near(tt, "calaud.pex", st[3], 0.069394841269841273, 1e-15);
    tk_true(tt, "calaud.ties", st[4] == 0);
    // mean gaps: X (7 + 4 - 1 - 5) / 4, Y (4 - 1 - 5 - 3) / 4
    tk_near(tt, "calaud.gx", gr[2], 1.25, 0);
    tk_near(tt, "calaud.gy", gr[3], -1.25, 0);
    tk_near(tt, "calaud.ss", st[5], 2 * 25.0 / 4, 1e-15);
    tk_true(tt, "calaud.noperm", std::isnan(st[6]) && std::isnan(gr[4]));
    // the permutation p-values: the exact ones are 40 / 70 for the
    // groups (|sum| >= 5 among the 70 ways to pick X's four areas)
    double s2[7];
    clab_calaud(a, r, o, g, nullptr, 8, 8, 2, 20000, 3, ar, gr, st, iw, wk);
    clab_calaud(a, r, o, g, nullptr, 8, 8, 2, 20000, 3, ar, gr, s2, iw, wk);
    tk_near(tt, "calaud.perm.same", s2[6], st[6], 0);
    tk_true(tt, "calaud.perm.range", gr[4] > 0 && gr[4] <= 1 && gr[4] == gr[5]);
  }
  // the exact count behind 40 / 70, by enumeration
  {
    const double gap[8] = {7, 4, 4, -1, -1, -5, -5, -3};
    int c = 0, tot = 0;
    for (int m = 0; m < 256; ++m) {
      int k = 0;
      double s = 0;
      for (int i = 0; i < 8; ++i) {
        if (m >> i & 1) {
          ++k;
          s += gap[i];
        }
      }
      if (k != 4) continue;
      ++tot;
      c += std::fabs(s) >= 5 ? 1 : 0;
    }
    tk_ieq(tt, "calaud.enum", c, 40);
    tk_ieq(tt, "calaud.enum.tot", tot, 70);
    const int a[8] = {0, 1, 2, 3, 4, 5, 6, 7}, g[8] = {0, 0, 1, 0, 1, 1, 0, 1};
    const double r[8] = {90, 80, 70, 60, 30, 20, 10, 5}, o[8] = {10, 25, 20, 80, 70, 90, 85, 60};
    double ar[56], gr[6], st[7], wk[14];
    int iw[24];
    clab_calaud(a, r, o, g, nullptr, 8, 8, 2, 20000, 3, ar, gr, st, iw, wk);
    const double pe = static_cast<double>(c) / tot, se = std::sqrt(pe * (1 - pe) / 20000);
    tk_lt(tt, "calaud.perm.exact", std::fabs(gr[4] - pe), 4 * se + 1e-4);
  }
  // 2. twelve areas (the Edgeworth series; the eight above were exact):
  // R's cor.test gives S = 74 and p = 0.008170639815899201
  {
    const double x[12] = {-0x1.40be8dcb2f8f7p-1, 0x1.7819fdbd4cf06p-3, -0x1.abd78373ce57ap-1,
                          0x1.986452991e538p+0,  0x1.516a7c3ee633p-2,  -0x1.a4146e9a60afdp-1,
                          0x1.f32099fd33d41p-2,  0x1.7a05b21cb7b8dp-1, 0x1.26ccd035ff319p-1,
                          -0x1.38b7bbbda8b07p-2, 0x1.83041735347ecp+0, 0x1.8f3310bbd0c95p-2};
    const double y[12] = {-0x1.3f68e64df7ce6p+0, -0x1.03f9a9551c1eap+1, 0x1.283edd171fa1cp-2,
                          0x1.8ce38dbbe823bp+0,  0x1.40d64e14ce06p-2,   0x1.f9508ac7c2e8p-4,
                          0x1.4f03b3df07689p+0,  0x1.550cc3ce0a0a8p+0,  0x1.7ea881f7d1dbep+0,
                          0x1.e8309ae4f55e9p-2,  0x1.961ac80a2f724p+0,  -0x1.997962eb057c9p+0};
    int a[12], g[12];
    for (int i = 0; i < 12; ++i) {
      a[i] = i;
      g[i] = 0;
    }
    double ar[84], gr[3], st[7], wk[15];
    int iw[36];
    clab_calaud(a, x, y, g, nullptr, 12, 12, 1, 0, 0, ar, gr, st, iw, wk);
    tk_near(tt, "calaud.s12", st[2], 74, 1e-12);
    tk_near(tt, "calaud.p12", st[3], 0.008170639815899201, 1e-14);
  }
  // 3. roll-up: records of areas 0 0 1 1 2 2 2, a population per area
  {
    const int a[7] = {0, 0, 1, 1, 2, 2, 2}, g[7] = {0, 1, 1, 1, 0, 0, 1};
    const double r[7] = {10, 20, 30, 40, 5, 6, 7}, o[7] = {1, 0, 1, 1, 0, 2, 1},
                 p[7] = {100, 100, 200, 200, 50, 50, 50};
    double ar[21], gr[6], st[7], wk[9];
    int iw[13];
    clab_calaud(a, r, o, g, p, 7, 3, 2, 0, 0, ar, gr, st, iw, wk);
    const double n[3] = {2, 2, 3}, m[3] = {15, 35, 6}, rt[3] = {0.01, 0.01, 0.06},
                 gg[3] = {0, 1, 0};
    for (int k = 0; k < 3; ++k) {
      tk_near(tt, "calaud.agg.n", ar[k], n[k], 0);
      tk_near(tt, "calaud.agg.risk", ar[3 + k], m[k], 0);
      tk_near(tt, "calaud.agg.rate", ar[6 + k], rt[k], 1e-16);
      tk_near(tt, "calaud.agg.grp", ar[9 + k], gg[k], 0);
    }
    // ties in the rates: average ranks, no exact p-value
    tk_near(tt, "calaud.agg.rk", ar[15], 2.5, 0);
    tk_true(tt, "calaud.agg.ties", st[4] == 1 && std::isnan(st[3]));
    const double pv[7] = {100, 90, 200, 200, 50, 50, 50};
    tk_true(tt, "calaud.guard.pop",
            clab_calaud(a, r, o, g, pv, 7, 3, 2, 0, 0, ar, gr, st, iw, wk) == CLAB_EDOM);
    tk_true(tt, "calaud.guard.na",
            clab_calaud(a, r, o, g, nullptr, 7, 2, 2, 0, 0, ar, gr, st, iw, wk) == CLAB_EARG);
  }
  // 4. scrdsp by hand: {1, 2, 3} and {4, 5, 6}, means 2 and 5
  {
    const double x[6] = {3, 1, 2, 6, 4, 5};
    const int g[6] = {0, 0, 0, 1, 1, 1};
    double gs[12], st[9];
    int iw[6];
    tk_true(tt, "scrdsp.rc", clab_scrdsp(x, g, 6, 2, -1, gs, st, iw) == CLAB_OK);
    tk_near(tt, "scrdsp.mean", gs[3], 5, 1e-15);
    tk_near(tt, "scrdsp.median", gs[4], 2, 0);
    tk_near(tt, "scrdsp.sd", gs[7], 1, 1e-15);
    tk_near(tt, "scrdsp.gap", gs[9], 3, 1e-15);
    tk_near(tt, "scrdsp.se", gs[11], std::sqrt(2.0 / 3), 1e-15);
    // SSB = 13.5 on 1, SSW = 4 on 4: F = 13.5; Welch's equals it here
    tk_near(tt, "scrdsp.f", st[0], 13.5, 1e-14);
    tk_near(tt, "scrdsp.df2", st[2], 4, 0);
    tk_near(tt, "scrdsp.wf", st[3], 13.5, 1e-14);
    tk_near(tt, "scrdsp.wdf", st[5], 4, 1e-14);
    tk_near(tt, "scrdsp.ref", st[7], 0, 0);
    // a group of one: no sd, no Welch
    const int g1[6] = {0, 0, 0, 0, 0, 1};
    clab_scrdsp(x, g1, 6, 2, 1, gs, st, iw);
    tk_true(tt, "scrdsp.one", std::isnan(gs[7]) && std::isnan(st[3]) && std::isnan(gs[11]));
    tk_near(tt, "scrdsp.med5", gs[4], 3, 0);
    // even counts: the median is the mean of the middle two
    const int g3[6] = {0, 0, 0, 0, 1, 1};
    clab_scrdsp(x, g3, 6, 2, -1, gs, st, iw);
    tk_near(tt, "scrdsp.med4", gs[4], 2.5, 0);
    const int g2[6] = {0, 0, 0, 0, 0, 0};
    tk_true(tt, "scrdsp.guard", clab_scrdsp(x, g2, 6, 2, -1, gs, st, iw) == CLAB_EDOM);
  }
}

}  // namespace

void t_calaud(tk &t, int) { run(t); }
