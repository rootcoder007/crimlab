// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for clab_ctabx and clab_fshx (src/tab.cpp): a 2 x 2
// table by hand, Fisher's p-value by direct summation over the support,
// the odds ratio's defining equations, and the refusals. Part of the
// self-test (selftest.cpp).
#include <cmath>

#include "crimlab.h"
#include "tchk.h"

namespace {

double ldh(double x, double m, double n, double k) {
  return std::lgamma(m + 1) - std::lgamma(x + 1) - std::lgamma(m - x + 1) + std::lgamma(n + 1) -
         std::lgamma(k - x + 1) - std::lgamma(n - k + x + 1) - std::lgamma(m + n + 1) +
         std::lgamma(k + 1) + std::lgamma(m + n - k + 1);
}

void tab(tk &t) {
  // 1. [[12, 9], [5, 21]] by hand
  const double x[4] = {12, 5, 9, 21};
  double e[4], res[4], st[9];
  const int dg[2] = {0, 1};
  tk_ieq(t, "tab.rc", clab_ctabx(x, 2, 2, 0, dg, e, res, st), CLAB_OK);
  const double n = 47, e11 = 21.0 * 17 / n;
  tk_eq(t, "tab.e11", e[0], e11, 1e-15);
  double ch = 0;
  const double r[2] = {21, 26}, c[2] = {17, 30};
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      const double ex = r[i] * c[j] / n, o = x[j * 2 + i];
      ch += (o - ex) * (o - ex) / ex;
    }
  }
  tk_eq(t, "tab.chisq", st[1], ch, 1e-14);
  tk_eq(t, "tab.phi2", st[7] * st[7], ch / n, 1e-14);
  tk_eq(t, "tab.v", st[4], std::fabs(st[7]), 1e-14);
  tk_eq(t, "tab.diag", st[8], 33 / n, 1e-15);
  tk_eq(t, "tab.df", st[2], 1, 0);
  // Yates: |o - e| - 1/2 throughout
  tk_ieq(t, "tab.rcy", clab_ctabx(x, 2, 2, 1, nullptr, e, res, st), CLAB_OK);
  const double dy = std::fabs(12 - e11) - 0.5;
  double chy = 0;
  for (int q = 0; q < 4; ++q) chy += dy * dy / e[q];
  tk_eq(t, "tab.yates", st[1], chy, 1e-13);
  tk_nan(t, "tab.nodiag", st[8]);
  // 2. Fisher: the p-value by summing the support directly
  double fs[5], rem[6];
  tk_ieq(t, "fsh.rc", clab_fshx(x, 2, 2, 1e6, 0.95, fs, rem), CLAB_OK);
  const double m = 17, nn = 30, k = 21, lo = 0, hi = 17;
  const double l0 = ldh(12, m, nn, k);
  double p = 0;
  for (double v = lo; v <= hi; v += 1) {
    const double l = ldh(v, m, nn, k);
    if (l <= l0 + 1e-7) p += std::exp(l);
  }
  tk_eq(t, "fsh.p", fs[0], p, 1e-12);
  tk_eq(t, "fsh.cnt", fs[4], hi - lo + 1, 0);
  // the estimate solves E[X; or] = 12; the ends P(X >= 12; lo) = 0.025
  auto mom = [&](double th, int which) {
    double s = 0, a = 0;
    for (double v = lo; v <= hi; v += 1) {
      const double d = std::exp(ldh(v, m, nn, k) + v * std::log(th));
      s += d;
      a += which == 0 ? d * v : (which == 1 ? (v >= 12 ? d : 0) : (v <= 12 ? d : 0));
    }
    return a / s;
  };
  tk_eq(t, "fsh.mle", mom(fs[1], 0), 12, 1e-9);
  tk_eq(t, "fsh.lo", mom(fs[2], 1), 0.025, 1e-9);
  tk_eq(t, "fsh.hi", mom(fs[3], 2), 0.025, 1e-9);
  // 3. a 2 x 3 table: the p-value at most 1, at least the observed term
  const double y[6] = {3, 0, 1, 4, 2, 2};
  tk_ieq(t, "fsh.rc2", clab_fshx(y, 2, 3, 1e6, 0.95, fs, rem), CLAB_OK);
  tk_true(t, "fsh.range", fs[0] > 0 && fs[0] <= 1 && std::isnan(fs[1]));
  // 4. refusals
  const double z[4] = {0, 0, 1, 2};
  tk_ieq(t, "tab.emptyrow", clab_ctabx(z, 2, 2, 0, nullptr, e, res, st), CLAB_EDOM);
  const double h[4] = {1.5, 1, 1, 1};
  tk_ieq(t, "fsh.whole", clab_fshx(h, 2, 2, 1e6, 0.95, fs, rem), CLAB_EDOM);
  tk_ieq(t, "fsh.cap", clab_fshx(x, 2, 2, 3, 0.95, fs, rem), CLAB_EDOM);
  tk_ieq(t, "fsh.dim", clab_fshx(x, 1, 4, 1e6, 0.95, fs, rem), CLAB_EARG);
  const int bad[2] = {0, 5};
  tk_ieq(t, "tab.dg", clab_ctabx(x, 2, 2, 0, bad, e, res, st), CLAB_EARG);
}

}  // namespace

void t_tab(tk &t, int big) {
  (void)big;
  tab(t);
}
