// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/kscn.cpp: the most likely cylinder against brute force over
// every cylinder, the ratio by hand, a planted cluster found and tested,
// non-overlapping secondary clusters, and determinism.
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "tchk.h"

namespace {

using vd = std::vector<double>;

double lr(double c, double e, double N) {
  if (!(c > e)) return 0;
  return c * std::log(c / e) + (N - c > 0 ? (N - c) * std::log((N - c) / (N - e)) : 0);
}

void famtest(tk &T, int big) {
  rng g = {4242};
  // brute force: every centre, radius, start and length
  {
    const int n = big ? 120 : 50;
    vd x(n), y(n), t(n);
    for (int i = 0; i < n; ++i) {
      x[i] = runi(g, 0, 10);
      y[i] = runi(g, 0, 10);
      t[i] = std::floor(runi(g, 0, 30));
    }
    const vd r = {1, 2.5}, tw = {3, 8};
    vd out(10);
    int ncl = 0;
    tk_ieq(T, "ks.ok",
           clab_kscan(x.data(), y.data(), t.data(), n, r.data(), 2, tw.data(), 2, 2, 0, 1, 1,
                      out.data(), &ncl),
           CLAB_OK);
    double bmax = 0;
    for (int c = 0; c < n; ++c)
      for (double rr : r)
        for (int a = 0; a < n; ++a)
          for (double len : tw) {
            double ns = 0, nt = 0, ci = 0;
            for (int j = 0; j < n; ++j) {
              const bool in =
                  (x[j] - x[c]) * (x[j] - x[c]) + (y[j] - y[c]) * (y[j] - y[c]) <= rr * rr;
              const bool tin = t[j] >= t[a] && t[j] < t[a] + len;
              ns += in;
              nt += tin;
              ci += in && tin;
            }
            if (ci >= 2) bmax = std::fmax(bmax, lr(ci, ns * nt / n, n));
          }
    tk_ieq(T, "ks.ncl", ncl, 1);
    tk_eq(T, "ks.brute", out[8], bmax, 1e-12);
    tk_eq(T, "ks.cell", out[8], lr(out[4], out[5], n), 1e-12);
    tk_eq(T, "ks.exp", out[5], out[6] * out[7] / n, 1e-15);
    tk_nan(T, "ks.nop", out[9]);
  }
  // the ratio by hand: 4 in a cylinder expecting 1, of 10
  {
    // ten incidents: four at one place in days 0 to 1, six spread out
    const vd x = {0, 0, 0, 0, 50, 60, 70, 80, 90, 100}, y = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    const vd t = {0, 0, 1, 1, 0, 5, 9, 13, 17, 21};
    const double r = 1, tw = 2;
    vd out(10);
    int ncl = 0;
    clab_kscan(x.data(), y.data(), t.data(), 10, &r, 1, &tw, 1, 2, 0, 1, 1, out.data(), &ncl);
    // circle of 4, window [0, 2) holds 5: e = 4 * 5 / 10 = 2
    tk_eq(T, "ks.hand", out[8], 4 * std::log(2.0) + 6 * std::log(6.0 / 8), 1e-14);
    tk_eq(T, "ks.hande", out[5], 2, 0);
    tk_eq(T, "ks.end", out[3], 2, 0);
  }
  // a planted cluster is found, tested and alone in its circle
  {
    const int n = 160;
    vd x(n), y(n), t(n);
    for (int i = 0; i < n; ++i) {
      x[i] = runi(g, 0, 20);
      y[i] = runi(g, 0, 20);
      t[i] = std::floor(runi(g, 0, 100));
    }
    for (int i = 0; i < 20; ++i) {
      x[i] = 5 + runi(g, -0.5, 0.5);
      y[i] = 15 + runi(g, -0.5, 0.5);
      t[i] = 40 + std::floor(runi(g, 0, 5));
    }
    const vd r = {1, 3}, tw = {5, 20};
    vd out(30), o2(30);
    int ncl = 0, n2 = 0;
    clab_kscan(x.data(), y.data(), t.data(), n, r.data(), 2, tw.data(), 2, 2, 99, 5, 3, out.data(),
               &ncl);
    clab_kscan(x.data(), y.data(), t.data(), n, r.data(), 2, tw.data(), 2, 2, 99, 5, 3, o2.data(),
               &n2);
    tk_true(
        T, "ks.found",
        ncl >= 1 && std::fabs(x[static_cast<int>(out[0])] - 5) < 3 && out[2] >= 30 && out[2] <= 45);
    tk_eq(T, "ks.p", out[9], 0.01, 0);
    tk_veq(T, "ks.seed", out.data(), o2.data(), 10 * ncl, 0);
    bool apart = true, desc = true;
    for (int a = 0; a < ncl; ++a)
      for (int b = a + 1; b < ncl; ++b) {
        const int i = static_cast<int>(out[10 * a]), j = static_cast<int>(out[10 * b]);
        const double d = std::hypot(x[i] - x[j], y[i] - y[j]);
        apart = apart && d > out[10 * a + 1] + out[10 * b + 1];
        desc = desc && out[10 * a + 8] >= out[10 * b + 8] && out[10 * a + 9] <= out[10 * b + 9];
      }
    tk_true(T, "ks.apart", apart && desc);
  }
  // refusals
  {
    const double x[3] = {0, 1, 2}, r2[2] = {2, 1}, one = 1;
    double o[10];
    int k;
    tk_ieq(T, "ks.n", clab_kscan(x, x, x, 2, &one, 1, &one, 1, 2, 0, 1, 1, o, &k), CLAB_EARG);
    tk_ieq(T, "ks.r", clab_kscan(x, x, x, 3, r2, 2, &one, 1, 2, 0, 1, 1, o, &k), CLAB_EDOM);
  }
}

}  // namespace

void t_kscn(tk &T, int big) { famtest(T, big); }
