// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/prim.cpp: small cases by hand.
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "tchk.h"

void t_prim(tk &T, int big) {
  (void)big;
  const std::vector<double> x = {2, 4, 4, 4, 5, 5, 7, 9};
  const int n = 8;
  tk_eq(T, "pr.mean", clab_mean(x.data(), n), 5, 1e-15);
  tk_eq(T, "pr.var", clab_var(x.data(), n, 1), 32.0 / 7, 1e-15);
  tk_eq(T, "pr.var0", clab_var(x.data(), n, 0), 4, 1e-15);
  tk_nan(T, "pr.var1", clab_var(x.data(), 1, 1));
  tk_nan(T, "pr.mean0", clab_mean(x.data(), 0));
  const double big2[2] = {1e308, 1e308};
  tk_eq(T, "pr.over", clab_mean(big2, 2), 1e308, 1e-15);
  const std::vector<double> y = {1, 2, 3, 4, 5, 6, 7, 8};
  const std::vector<double> z = {8, 7, 6, 5, 4, 3, 2, 1};
  tk_eq(T, "pr.cor", clab_cor(y.data(), z.data(), n), -1, 1e-15);
  tk_nan(T, "pr.cor0", clab_cor(x.data(), x.data(), 1));
  const double a[2] = {0, 0}, b[2] = {3, 4};
  tk_eq(T, "pr.dist", clab_edist(a, b, 2), 5, 1e-15);
  double d[2];
  clab_dnorm(a, 1, 0, 1, 0, d);
  tk_eq(T, "pr.dn", d[0], 1 / std::sqrt(2 * 3.14159265358979323846), 1e-15);
  clab_dnorm(b, 1, 0, 1, 1, d);
  tk_eq(T, "pr.ldn", d[0], -4.5 - 0.5 * std::log(2 * 3.14159265358979323846), 1e-15);
  double m[4];
  clab_moms(x.data(), n, m);
  tk_eq(T, "pr.m1", m[0], 5, 1e-15);
  tk_eq(T, "pr.m2", m[1], 32.0 / 7, 1e-14);
  tk_eq(T, "pr.m3", m[2], (42.0 / 8) / std::pow(4.0, 1.5), 1e-14);
  tk_eq(T, "pr.m4", m[3], (356.0 / 8) / 16 - 3, 1e-14);
  const double p[3] = {0.25, 0.5, 0.1};
  double q[3];
  clab_quant(x.data(), n, p, 3, q);
  tk_eq(T, "pr.q25", q[0], 4, 1e-15);
  tk_eq(T, "pr.q50", q[1], 4.5, 1e-15);
  tk_eq(T, "pr.q10", q[2], 3.4, 1e-14);
  tk_eq(T, "pr.med", clab_med(x.data(), n), 4.5, 1e-15);
  tk_eq(T, "pr.mad", clab_mad(x.data(), n, 1), 0.5, 1e-15);
  const std::vector<double> w = {1, 2, 3, 4, 5, 6, 7, 8, 9, 100};
  tk_eq(T, "pr.trim", clab_trim(w.data(), 10, 0.1), 5.5, 1e-15);
  tk_eq(T, "pr.wins", clab_wins(w.data(), 10, 0.1), 5.5, 1e-15);
  tk_eq(T, "pr.trim5", clab_trim(w.data(), 10, 0.5), 5.5, 1e-15);
  tk_nan(T, "pr.trimbad", clab_trim(w.data(), 10, 0.7));
  double o[2];
  const double wx[4] = {10, 20, 30, 40}, ww[4] = {1, 1, 1, 1};
  clab_wtd(wx, ww, 4, o);
  tk_eq(T, "pr.wm", o[0], 25, 1e-15);
  tk_eq(T, "pr.wv", o[1], 500.0 / 3, 1e-14);
  const double t1[5] = {1, 2, 2, 3, 3}, t2[5] = {1, 4, 9, 16, 25};
  double r[5];
  clab_midrk(t1, 5, r);
  tk_eq(T, "pr.rk", r[1], 2.5, 1e-15);
  tk_eq(T, "pr.sp", clab_spear(y.data(), y.data(), n), 1, 1e-15);
  tk_eq(T, "pr.sp2", clab_spear(t2, t2, 5), 1, 1e-15);
  const double X[6] = {1, 2, 3, 2, 4, 6};
  double C[4];
  clab_covm(X, 3, 2, C);
  tk_eq(T, "pr.cv", C[1], 2, 1e-15);
  tk_eq(T, "pr.cv22", C[3], 4, 1e-15);
  std::vector<double> bm(50);
  tk_ieq(T, "pr.boot", clab_bmean(x.data(), n, 50, 7, bm.data()), CLAB_OK);
  bool in = true;
  for (double v : bm) in = in && v >= 2 && v <= 9;
  tk_true(T, "pr.bootin", in);
  const double tr[2] = {1, 0}, ps[2] = {0.001, 0.5};
  clab_ipw(tr, ps, 2, 0.01, 0.99, o);
  tk_eq(T, "pr.ipw", o[0], 100, 1e-12);
  tk_eq(T, "pr.ipw0", o[1], 2, 1e-15);
  double u[3];
  clab_unif(0, 3, u);
  bool ok = true;
  for (double v : u) ok = ok && v >= 0 && v < 1;
  tk_true(T, "pr.unif", ok);
}
