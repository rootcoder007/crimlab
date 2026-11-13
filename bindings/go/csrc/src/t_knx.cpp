// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/knx.cpp: the Knox table against brute force over all pairs,
// its expectation against the mean over every permutation of the times,
// and the permutation p-values' range and determinism.
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "crimlab.h"
#include "tchk.h"

namespace {

using vd = std::vector<double>;

// brute force: every pair, bands (lo, hi]
void brute(const vd &x, const vd &y, const vd &t, const vd &ds, const vd &dt, vd &obs) {
  const int n = static_cast<int>(x.size()), ks = static_cast<int>(ds.size());
  const int kt = static_cast<int>(dt.size());
  obs.assign(ks * kt, 0.0);
  for (int i = 0; i < n; ++i)
    for (int j = i + 1; j < n; ++j) {
      const double d = std::sqrt((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]));
      const double e = std::fabs(t[i] - t[j]);
      int s = 0, b = 0;
      while (s < ks && d > ds[s]) ++s;
      while (b < kt && e > dt[b]) ++b;
      if (s < ks && b < kt) obs[s + b * ks] += 1;
    }
}

void famtest(tk &T, int big) {
  rng g = {2024};
  // counts against brute force, with ties in time and in place
  {
    const int n = big ? 400 : 120;
    vd x(n), y(n), t(n);
    for (int i = 0; i < n; ++i) {
      x[i] = std::floor(runi(g, 0, 50));
      y[i] = std::floor(runi(g, 0, 50));
      t[i] = std::floor(runi(g, 0, 60));
    }
    const vd ds = {2, 5, 9}, dt = {3, 7, 14, 30};
    vd out(3 * 12 + 1), obs;
    tk_ieq(
        T, "knx.ok",
        clab_knxnr(x.data(), y.data(), t.data(), n, ds.data(), 3, dt.data(), 4, 0, 1, out.data()),
        CLAB_OK);
    brute(x, y, t, ds, dt, obs);
    tk_veq(T, "knx.obs", out.data(), obs.data(), 12, 0);
    tk_nan(T, "knx.nop", out[24]);
    // the expectation: (pairs in the distance band) (pairs in the time
    // band) / all pairs, by brute force
    vd sc(3, 0), tc(4, 0);
    for (int i = 0; i < n; ++i)
      for (int j = i + 1; j < n; ++j) {
        const double d = std::sqrt((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]));
        const double e = std::fabs(t[i] - t[j]);
        int s = 0, b = 0;
        while (s < 3 && d > ds[s]) ++s;
        while (b < 4 && e > dt[b]) ++b;
        if (s < 3) sc[s] += 1;
        if (b < 4) tc[b] += 1;
      }
    const double np = n * (n - 1) / 2.0;
    vd ex(12);
    for (int b = 0; b < 4; ++b)
      for (int s = 0; s < 3; ++s) ex[s + b * 3] = sc[s] * tc[b] / np;
    tk_veq(T, "knx.exp", out.data() + 12, ex.data(), 12, 1e-15);
    tk_eq(T, "knx.pairs", out[36], sc[0] + sc[1] + sc[2], 0);
    // p-values: in [1 / (nsim + 1), 1], 1 where nothing is observed, the
    // same for the same seed
    vd o2(37), o3(37);
    clab_knxnr(x.data(), y.data(), t.data(), n, ds.data(), 3, dt.data(), 4, 99, 7, o2.data());
    clab_knxnr(x.data(), y.data(), t.data(), n, ds.data(), 3, dt.data(), 4, 99, 7, o3.data());
    tk_veq(T, "knx.seed", o2.data() + 24, o3.data() + 24, 12, 0);
    bool rng01 = true;
    for (int c = 0; c < 12; ++c) {
      rng01 = rng01 && o2[24 + c] >= 0.01 && o2[24 + c] <= 1;
      if (o2[c] == 0) rng01 = rng01 && o2[24 + c] == 1;
    }
    tk_true(T, "knx.prange", rng01);
  }
  // the expectation is exact: the mean count over all 720 orders of the
  // times of six incidents
  {
    const vd x = {0, 1, 0, 3, 7, 7.5}, y = {0, 0, 2, 1, 6, 6};
    vd t = {1, 2, 2, 5, 9, 20};
    const vd ds = {1.5, 4}, dt = {1, 4};
    vd out(13), mean(4, 0), o(13);
    clab_knxnr(x.data(), y.data(), t.data(), 6, ds.data(), 2, dt.data(), 2, 0, 1, out.data());
    std::vector<int> pm = {0, 1, 2, 3, 4, 5};
    int np = 0;
    do {
      vd tp(6);
      for (int i = 0; i < 6; ++i) tp[i] = t[pm[i]];
      clab_knxnr(x.data(), y.data(), tp.data(), 6, ds.data(), 2, dt.data(), 2, 0, 1, o.data());
      for (int c = 0; c < 4; ++c) mean[c] += o[c];
      ++np;
    } while (std::next_permutation(pm.begin(), pm.end()));
    for (int c = 0; c < 4; ++c) mean[c] /= np;
    tk_veq(T, "knx.exact", out.data() + 4, mean.data(), 4, 1e-14);
    // by hand: within 4 are (0,1) at 1, (0,2) at 2, (1,2) and (1,3) at
    // sqrt 5, (0,3) and (2,3) at sqrt 10, (4,5) at 0.5; the time gaps put
    // (0,1) in cell (1, 1), (0,2) and (1,2) in (2, 1), (0,3), (1,3) and
    // (2,3) in (2, 2), and (4,5) (11 days) in none
    tk_eq(T, "knx.hand", out[12], 7, 0);
    const vd hand = {1, 2, 0, 3};
    tk_veq(T, "knx.hand2", out.data(), hand.data(), 4, 0);
  }
  // refusals
  {
    const double x[3] = {0, 1, 2}, t[3] = {0, 1, 2}, bad[2] = {2, 1}, e[1] = {1};
    double o[4];
    tk_ieq(T, "knx.n", clab_knxnr(x, x, t, 2, e, 1, e, 1, 0, 1, o), CLAB_EARG);
    tk_ieq(T, "knx.edges", clab_knxnr(x, x, t, 3, bad, 2, e, 1, 0, 1, o), CLAB_EDOM);
    const double nan[3] = {0, std::numeric_limits<double>::quiet_NaN(), 2};
    tk_ieq(T, "knx.fin", clab_knxnr(x, x, nan, 3, e, 1, e, 1, 0, 1, o), CLAB_EDOM);
  }
}

}  // namespace

void t_knx(tk &T, int big) { famtest(T, big); }
