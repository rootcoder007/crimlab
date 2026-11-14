// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for clab_rswp and clab_rshl (src/sens.cpp). Part of the
// self-test (selftest.cpp).

#include <cmath>
#include <cstdint>

#include "crimlab.h"
#include "tchk.h"

namespace {

// Reference values, dev/ref/cpp/sens.R (rbounds 2.2, R 4.5.0), for the 60
// differences d_i = (((i * 37) %% 101) - 35) / 16, Gamma = 1, 1.5, ..., 3:
// ref_plo, ref_pup: psens(d, 0, Gamma = 3, GammaInc = 0.5)$bounds (rbounds
// rounds them to 4 digits); ref_pup_u, ref_plo_u: the same formulas in
// base R without rounding (pnorm(lower.tail = FALSE)); ref_hl1_*, ref_hl4_*:
// hlsens(d, 0, pr = 0.1 or 1e-4, Gamma = 3, GammaInc = 0.5)$bounds (signif
// 5); ref_ex_*: the bounds by sorting all 1,830 Walsh averages in base R.
const double ref_plo[5] = {0.00050000000000000001, 0, 0, 0, 0};
const double ref_pup[5] = {0.00050000000000000001, 0.0241, 0.13550000000000001, 0.33040000000000003, 0.53890000000000005};
const double ref_pup_u[5] = {0.00051288079178312893, 0.024072917724745201, 0.13545262279122386, 0.33039265117440919, 0.53893768680597998};
const double ref_plo_u[5] = {0.00051288079178312893, 1.1445926742216369e-06, 2.2580616428256174e-09, 4.2692755559615404e-12, 7.9319870552254062e-15};
const double ref_hl1_lo[5] = {0.90627000000000002, 0.50627, 0.30626999999999999, 0.10627, -0.093727000000000005};
const double ref_hl1_hi[5] = {0.90627000000000002, 1.3063, 1.5063, 1.7062999999999999, 1.9063000000000001};
const double ref_hl4_lo[5] = {0.90627000000000002, 0.59367000000000003, 0.31247000000000003, 0.12497, -0.031327000000000001};
const double ref_hl4_hi[5] = {0.90627000000000002, 1.2188000000000001, 1.5001, 1.6876, 1.8438000000000001};
const double ref_ex_lo[5] = {0.90625, 0.59375, 0.3125, 0.125, -0.03125};
const double ref_ex_hi[5] = {0.90625, 1.234375, 1.5, 1.6875, 1.84375};
const double G[5] = {1, 1.5, 2, 2.5, 3};

struct xrng {
  std::uint64_t s;
  double unif() {
    s ^= s >> 12;
    s ^= s << 25;
    s ^= s >> 27;
    return static_cast<double>((s * 2685821657736338717ULL) >> 11) * 0x1p-53;
  }
};

void dsample(double *d) {
  for (int i = 0; i < 60; ++i) d[i] = (((i * 37) % 101) - 35) / 16.0;
}

void run(tk &tt) {
  static double d[600], w[600];
  static int iw[600];
  double plo[5], pup[5], lo[5], hi[5], gc, ts;
  dsample(d);
  // 1. p-value bounds: parity with rbounds (its 4-digit rounding) and with
  // the unrounded formulas
  tk_true(tt, "sens.p.rc", clab_rswp(d, 60, G, 5, 1.6448536269514722, plo, pup, &gc, &ts, w, iw) == CLAB_OK);
  for (int k = 0; k < 5; ++k) {
    tk_near(tt, "sens.p.up4", pup[k], ref_pup[k], 5e-5);
    tk_near(tt, "sens.p.lo4", plo[k], ref_plo[k], 5e-5);
    tk_near(tt, "sens.p.up", pup[k] / ref_pup_u[k], 1, 1e-12);
    tk_near(tt, "sens.p.lo", plo[k] / ref_plo_u[k], 1, 1e-10);
  }
  tk_near(tt, "sens.p.g1", pup[0], plo[0], 0);
  // 2. the critical Gamma: the upper bound there equals 1 - Phi(z)
  for (double z = 0.5; z < 4; z += 0.7) {
    clab_rswp(d, 60, G, 0, z, plo, pup, &gc, &ts, w, iw);
    double p2[1], q2[1];
    clab_rswp(d, 60, &gc, gc >= 1 ? 1 : 0, z, q2, p2, &gc, &ts, w, iw);
    if (gc >= 1) tk_near(tt, "sens.gcrit", p2[0] / (0.5 * std::erfc(z / std::sqrt(2.0))), 1, 1e-12);
  }
  // 3. the rank sums against a direct O(n^2) midrank count
  {
    double t = 0;
    for (int i = 0; i < 60; ++i) {
      if (d[i] <= 0) continue;
      double less = 0, eq = 0;
      for (int j = 0; j < 60; ++j) {
        if (d[j] == 0) continue;
        less += std::fabs(d[j]) < std::fabs(d[i]);
        eq += std::fabs(d[j]) == std::fabs(d[i]);
      }
      t += less + (eq + 1) / 2;
    }
    clab_rswp(d, 60, G, 0, 1.0, plo, pup, &gc, &ts, w, iw);
    tk_near(tt, "sens.T", ts, t, 0);
  }
  // 4. Hodges-Lehmann bounds: exact equality with the sorted Walsh
  // averages; parity with rbounds at its fine grid where the bound is a
  // single Walsh average; at its default grid (pr = 0.1) its bounds lie within
  // one grid step of the exact ones; at Gamma = 1.5 the upper bound falls on
  // a flat of T(tau) = mu, where the estimate is the midpoint of the flat
  // (Hodges and Lehmann; Rosenbaum 1993): 1.234375, while the grid search
  // stops at its first point, 1.2188 (pr = 1e-4)
  tk_true(tt, "sens.hl.rc", clab_rshl(d, 60, G, 5, lo, hi, w) == CLAB_OK);
  for (int k = 0; k < 5; ++k) {
    tk_near(tt, "sens.hl.ex.lo", lo[k], ref_ex_lo[k], 0);
    tk_near(tt, "sens.hl.ex.hi", hi[k], ref_ex_hi[k], 0);
    tk_near(tt, "sens.hl.fine.lo", lo[k], ref_hl4_lo[k], 2e-4);
    if (k != 1) tk_near(tt, "sens.hl.fine.hi", hi[k], ref_hl4_hi[k], 2e-4);
  }
  tk_true(tt, "sens.hl.flat", std::fabs(hi[1] - ref_hl4_hi[1]) > 0.01);
  tk_true(tt, "sens.hl.grid", std::fabs(lo[1] - ref_hl1_lo[1]) > 0.08);
  // and every default-grid bound lies within one grid step (0.1) of the
  // exact one
  for (int k = 0; k < 5; ++k) {
    tk_true(tt, "sens.hl.grid.lo", std::fabs(lo[k] - ref_hl1_lo[k]) <= 0.1);
    tk_true(tt, "sens.hl.grid.hi", std::fabs(hi[k] - ref_hl1_hi[k]) <= 0.1);
  }
  // 5. brute force: random differences, all Walsh averages sorted
  // directly, Gamma including values that make a an integer
  {
    xrng r{21};
    int bad = 0;
    for (int rep = 0; rep < 30; ++rep) {
      const int n = 2 + static_cast<int>(r.unif() * 40);
      for (int i = 0; i < n; ++i) d[i] = std::floor(r.unif() * 40 - 12) / 4 + (rep % 2 ? r.unif() : 0);
      static double all[900];
      int m = 0;
      for (int i = 0; i < n; ++i)
        for (int j = i; j < n; ++j) all[m++] = 0.5 * d[i] + 0.5 * d[j];
      for (int i = 1; i < m; ++i)
        for (int j = i; j > 0 && all[j] < all[j - 1]; --j) {
          const double t = all[j]; all[j] = all[j - 1]; all[j - 1] = t;
        }
      const double gs[5] = {1, 2, 3, 1.25, 7.5};
      clab_rshl(d, n, gs, 5, lo, hi, w);
      for (int k = 0; k < 5; ++k) {
        // ceil and floor + 1 by integer arithmetic where Gamma is a
        // multiple of 1/4: a = 4 g M / (4 + 4 g)
        const std::int64_t num = static_cast<std::int64_t>(4 * gs[k]) * m, den = 4 + static_cast<std::int64_t>(4 * gs[k]);
        const std::int64_t cu = (num + den - 1) / den, fu = num / den + 1;
        const std::int64_t cl = (4LL * m + den - 1) / den, fl = 4LL * m / den + 1;
        bad += hi[k] != 0.5 * all[cu - 1] + 0.5 * all[fu - 1];
        bad += lo[k] != 0.5 * all[cl - 1] + 0.5 * all[fl - 1];
      }
    }
    tk_true(tt, "sens.hl.brute", bad == 0);
  }
  // 6. guards
  {
    const double g0[1] = {0.5};
    tk_true(tt, "sens.guard.g", clab_rshl(d, 5, g0, 1, lo, hi, w) == CLAB_EDOM);
    tk_true(tt, "sens.guard.z", clab_rswp(d, 5, G, 1, 0, plo, pup, &gc, &ts, w, iw) == CLAB_EDOM);
    const double z0[2] = {0, 0};
    tk_true(tt, "sens.zero", clab_rswp(z0, 2, G, 1, 1, plo, pup, &gc, &ts, w, iw) == CLAB_OK && std::isnan(pup[0]));
  }
}


}  // namespace

void t_sens(tk &t, int) { run(t); }

#ifdef CLAB_MAIN
#include <chrono>
#include <cstdio>

// timing: 50,000 pairs, Gamma 1 to 2 by 0.25 (rbounds psens took 0.024 s,
// hlsens 3.112 s, one run each, arm64 Linux)
void b_sens() {
  const int n = 50000;
  static double d[n], w[n];
  static int iw[n];
  xrng r{1};
  for (int i = 0; i < n; ++i) {
    const double a = 1.0 - r.unif(), b = r.unif(), c = 1.0 - r.unif(), e = r.unif();
    d[i] = 0.3 + std::sqrt(-2 * std::log(a)) * std::cos(6.283185307179586 * b) -
           std::sqrt(-2 * std::log(c)) * std::cos(6.283185307179586 * e);
  }
  const double g[5] = {1, 1.25, 1.5, 1.75, 2};
  double plo[5], pup[5], lo[5], hi[5], gc, ts, t[5], t2[5];
  for (int k = 0; k < 5; ++k) {
    auto a = std::chrono::steady_clock::now();
    clab_rswp(d, n, g, 5, 1.6448536269514722, plo, pup, &gc, &ts, w, iw);
    t[k] = std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count();
    a = std::chrono::steady_clock::now();
    clab_rshl(d, n, g, 5, lo, hi, w);
    t2[k] = std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count();
  }
  std::printf("BENCH clab_rswp n=5e4, 5 Gammas: median %.6f s (k=5)\n", tk_med5(t));
  std::printf("BENCH clab_rshl n=5e4, 5 Gammas: median %.6f s (k=5)\n", tk_med5(t2));
}
#endif
