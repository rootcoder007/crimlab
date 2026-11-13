// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for clab_hypq and clab_sqrtlw (src/hyq.cpp): whatever
// the dispatch orders, the number of busy units of the hypercube is the
// Erlang loss (M/M/N/N) or delay (M/M/N) system, so the full state has
// Erlang's B or C probability; one unit by hand; the dispatch shares add
// to the calls served; the square-root law's constants by hand. Part of
// the self-test (selftest.cpp).
#include <cmath>

#include "crimlab.h"
#include "tchk.h"

namespace {

double erlb(double a, int n) {
  double b = 1;
  for (int k = 1; k <= n; ++k) b = a * b / (k + a * b);
  return b;
}

void hyq(tk &t) {
  static double w[64 * 64 + 64];
  double p[64], wl[6], d[12], st[3];
  // three units, four atoms; the orders are the rows of 4 x 3 matrices
  // (column-major): all 0 1 2, or 0 1 2 / 2 1 0 / 2 0 1 / 1 2 0
  const double lam[4] = {0.4, 0.7, 0.3, 0.6};
  const int pref[12] = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2};
  const int alt[12] = {0, 2, 2, 1, 1, 1, 0, 2, 2, 0, 1, 0};
  const double a = 2.0 / 1.5;
  tk_ieq(t, "hyq.rc", clab_hypq(lam, alt, 4, 3, 1.5, 0, p, wl, d, st, w), CLAB_OK);
  tk_eq(t, "hyq.erlb", st[0], erlb(a, 3), 1e-13);
  double tot = 0, sw = 0, sd = 0;
  for (int b = 0; b < 8; ++b) tot += p[b];
  for (int u = 0; u < 3; ++u) sw += wl[u];
  for (int i = 0; i < 12; ++i) sd += d[i];
  tk_eq(t, "hyq.sum", tot, 1, 1e-14);
  tk_eq(t, "hyq.busy", sw, a * (1 - st[0]), 1e-13);
  tk_eq(t, "hyq.disp", sd, 1 - st[0], 1e-13);
  tk_eq(t, "hyq.rho", st[2], 2.0 / 4.5, 1e-15);
  // the delay system: Erlang C
  tk_ieq(t, "hyq.wait.rc", clab_hypq(lam, pref, 4, 3, 1.5, 1, p, wl, d, st, w), CLAB_OK);
  const double b3 = erlb(a, 3), c3 = 3 * b3 / (3 - a * (1 - b3));
  tk_eq(t, "hyq.erlc", st[0], c3, 1e-13);
  tk_eq(t, "hyq.lq", st[1], c3 * st[2] / (1 - st[2]), 1e-13);
  // one unit: busy with probability lam / (lam + mu)
  const double l1[1] = {2};
  const int p1[1] = {0};
  tk_ieq(t, "hyq.one.rc", clab_hypq(l1, p1, 1, 1, 3, 0, p, wl, d, st, w), CLAB_OK);
  tk_eq(t, "hyq.one", st[0], 0.4, 1e-15);
  tk_eq(t, "hyq.one.d", d[0], 0.6, 1e-15);
  // the refusals
  const int dup[12] = {0, 0, 0, 0, 1, 0, 1, 1, 2, 2, 2, 2};
  tk_ieq(t, "hyq.dup", clab_hypq(lam, dup, 4, 3, 1.5, 0, p, wl, d, st, w), CLAB_EDOM);
  const double ln[4] = {0.4, 0, 0.3, 0.6};
  tk_ieq(t, "hyq.lam", clab_hypq(ln, pref, 4, 3, 1.5, 0, p, wl, d, st, w), CLAB_EDOM);
  tk_ieq(t, "hyq.sat", clab_hypq(lam, pref, 4, 3, 0.5, 1, p, wl, d, st, w), CLAB_EDOM);
  tk_ieq(t, "hyq.arg", clab_hypq(lam, pref, 4, 11, 1.5, 0, p, wl, d, st, w), CLAB_EARG);
}

void sqr(tk &t) {
  double d, c;
  tk_ieq(t, "sqr.rc", clab_sqrtlw(100, 4, 0, 0, &d, &c), CLAB_OK);
  tk_eq(t, "sqr.e", d, 2.5, 1e-15);
  tk_ieq(t, "sqr.rc2", clab_sqrtlw(100, 8, 0.5, 1, &d, &c), CLAB_OK);
  tk_eq(t, "sqr.r", d, std::sqrt(2 * 3.141592653589793) / 4 * 5, 1e-14);
  tk_ieq(t, "sqr.dom", clab_sqrtlw(100, 8, 1, 1, &d, &c), CLAB_EDOM);
}

}  // namespace

void t_hyq(tk &t, int big) {
  (void)big;
  hyq(t);
  sqr(t);
}
