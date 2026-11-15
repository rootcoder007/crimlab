// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/dsc.cpp against their definitions: bins of an even grid,
// perfect and tied correlations, runs, tallies, quantiles, distances of an
// orthogonal design, a log tail beyond underflow, Cochran's Q by hand, and
// Little's test unchanged by a row with nothing observed.
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "tchk.h"

void t_dsc(tk &T, int) {
  const double x[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  int cnt[10];
  clab_hbin(x, 10, 10, cnt);
  int one = 1;
  for (int b = 1; b < 10; ++b) one = one && cnt[b] == 1;
  tk_true(T, "dsc.hbin", cnt[0] == 1 && one);
  tk_ieq(T, "dsc.hflat", clab_hbin(x, 1, 3, cnt), CLAB_EDOM);
  double xy[20], cor[1];
  int np[1];
  for (int i = 0; i < 10; ++i) {
    xy[i] = x[i];
    xy[10 + i] = std::exp(x[i]);
  }
  clab_pcor(xy, 10, 2, 0, 3, cor, np);
  tk_eq(T, "dsc.spear", cor[0], 1, 1e-15);
  clab_pcor(xy, 10, 2, 1, 11, cor, np);
  tk_true(T, "dsc.minprs", std::isnan(cor[0]) && np[0] == 10);
  const int na[8] = {0, 1, 1, 1, 0, 1, 0, 1};
  int st[8], en[8];
  tk_ieq(T, "dsc.runs", clab_nrun(na, 8, 2, st, en), 1);
  tk_true(T, "dsc.run1", st[0] == 2 && en[0] == 4);
  double ms[8];
  clab_msum(na, 4, 2, ms);
  tk_eq(T, "dsc.miss", ms[2], 5, 0);
  tk_eq(T, "dsc.comp", ms[4], 1, 0);
  const double pr[2] = {0.5, 0.25};
  double q[2];
  clab_q7(x, 10, pr, 2, q);
  tk_eq(T, "dsc.q50", q[0], 5.5, 0);
  tk_eq(T, "dsc.q25", q[1], 3.25, 0);
  double ps[7];
  clab_pstat(x, 10, ps);
  tk_eq(T, "dsc.sd", ps[1], std::sqrt(55.0 / 6), 1e-14);
  tk_eq(T, "dsc.skew", ps[2], 0, 1e-14);
  // the corners of a square from its centre: equal distances
  const double sq[10] = {-1, 1, -1, 1, 0, -1, -1, 1, 1, 0};
  double d2[5];
  tk_ieq(T, "dsc.maha", clab_maha(sq, 5, 2, 0, d2), CLAB_OK);
  tk_eq(T, "dsc.mahaeq", d2[0], d2[3], 1e-12);
  tk_eq(T, "dsc.mahac", d2[4], 0, 1e-15);
  tk_eq(T, "dsc.lp", clab_lpchisq(3, 2), std::log(clab_pchisq(3, 2, 0)), 1e-12);
  tk_eq(T, "dsc.lpfar", clab_lpchisq(4000, 2), -2000, 1e-9);
  const double y[2] = {1, 3}, s[2] = {1, 1};
  double co[6];
  clab_cochq(y, s, 2, co);
  tk_eq(T, "dsc.q", co[2], 2, 1e-15);
  tk_eq(T, "dsc.i2", co[5], 0.5, 1e-15);
  const int code[6] = {2, 1, 2, 3, 2, 1};
  int lev[4], ct[4];
  tk_ieq(T, "dsc.freq", clab_freq(code, 6, 3, 1, 2, 1, lev, ct), 3);
  tk_true(T, "dsc.freqo", lev[0] == 2 && ct[0] == 3 && lev[2] == 4 && ct[2] == 1);
  int rows[6], rc[6];
  tk_ieq(T, "dsc.dups", clab_dups(code, 6, 3, rows, rc), 5);
  tk_true(T, "dsc.dupo", rows[0] == 2 && rows[1] == 6 && rc[2] == 3);
  int gs[4], sh[8];
  tk_ieq(T, "dsc.grid", clab_mgrid(na, 8, 1, 4, gs, sh), 4);
  tk_true(T, "dsc.shade", sh[1] == 4 && sh[0] == 2);
  // a row with every value missing carries no information: the EM fixed
  // point and Little's statistic are those of the table without it, and
  // the row is one more pattern and one more row used
  const double m6[12] = {1, 2, 3, 4, 5, 6, 2, 1, 4, NAN, 6, NAN};
  const double m7[14] = {1, 2, 3, 4, 5, 6, NAN, 2, 1, 4, NAN, 6, NAN, NAN};
  double o6[8], o7[8], mu6[2], mu7[2], sg6[4], sg7[4];
  tk_ieq(T, "dsc.mcar6", clab_mcar(m6, 6, 2, 100000, 1e-13, o6, mu6, sg6), CLAB_OK);
  tk_ieq(T, "dsc.mcar7", clab_mcar(m7, 7, 2, 100000, 1e-13, o7, mu7, sg7), CLAB_OK);
  tk_eq(T, "dsc.mcarst", o7[0], o6[0], 1e-9);
  tk_eq(T, "dsc.mcardf", o7[1], o6[1], 0);
  tk_eq(T, "dsc.mcarmu", mu7[1], mu6[1], 1e-9);
  tk_eq(T, "dsc.mcarsg", sg7[3], sg6[3], 1e-9);
  tk_eq(T, "dsc.mcarnp", o7[4], o6[4] + 1, 0);
  tk_eq(T, "dsc.mcarnu", o7[5], o6[5] + 1, 0);
}
