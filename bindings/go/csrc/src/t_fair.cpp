// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for clab_fair (src/fair.cpp). Part of the self-test (selftest.cpp).

#include <cmath>
#include <cstdint>

#include "crimlab.h"
#include "tchk.h"

namespace {

// xorshift64*: a small fixed PRNG so that every platform draws the same data
struct xrng {
  std::uint64_t s;
  double unif() {
    s ^= s >> 12;
    s ^= s << 25;
    s ^= s >> 27;
    return static_cast<double>((s * 2685821657736338717ULL) >> 11) * 0x1p-53;
  }
};

const double Z = 1.9599639845400534;  // qnorm(0.975) in R 4.5.0

// fairness 1.2.3, R 4.5.0, dev/doss/fairness.R: with i = 0..999,
// y = ((i * 13) %% 7) < 3, g = i %% 3 (labels a, b, c), p = (((i * 37) %% 101)
// + 30 y) / 131, and for each metric the call
// fairness::<fn>(d, "y", "g", probs = "p", cutoff = 0.5)$Metric[1, ] with
// <fn> = acc_parity, equal_odds, spec_parity, fpr_parity, fnr_parity,
// pred_rate_parity, npv_parity, prop_parity, dem_parity, mcc_parity.
const double ref[10][3] = {
    {0.64071856287425155, 0.65765765765765771, 0.64564564564564564},
    {0.63888888888888884, 0.65492957746478875, 0.63380281690140849},
    {0.64210526315789473, 0.65968586387434558, 0.65445026178010468},
    {0.35789473684210527, 0.34031413612565448, 0.34554973821989526},
    {0.3611111111111111, 0.34507042253521125, 0.36619718309859156},
    {0.57499999999999996, 0.58860759493670889, 0.57692307692307687},
    {0.70114942528735624, 0.72000000000000008, 0.70621468926553677},
    {0.47904191616766467, 0.47447447447447449, 0.46846846846846846},
    {160, 158, 156},
    {0.27856125645327512, 0.31159703894873814, 0.28568397364023473}};

// stats::prop.test(x, 100, correct = FALSE)$conf.int, R 4.5.0, for
// x = 0, 1, 7, 50, 99, 100 (Wilson score interval).
const double wx[6] = {0, 1, 7, 50, 99, 100};
const double wlo[6] = {0, 0.0017674320641406492, 0.034319261067272677,
                       0.40383153036599567, 0.94551380382129468,
                       0.96300650179301428};
const double whi[6] = {0.036993498206985651, 0.054486196178705287,
                       0.13749514739073498, 0.59616846963400427,
                       0.99823256793585935, 1};

void detsample(double *p, int *y, int *g, int n) {
  for (int i = 0; i < n; ++i) {
    y[i] = ((i * 13) % 7) < 3 ? 1 : 0;
    g[i] = i % 3;
    p[i] = (((i * 37) % 101) + 30.0 * y[i]) / 131.0;
  }
}

void run(tk &tt) {
  // 1. hand-derived: group 0 has tp 2, fp 1, tn 3, fn 1; group 1 is
  // predicted all positive (tp 2, fp 1, tn 0, fn 0).
  {
    const double s[10] = {.9, .8, .7, .1, .2, .3, .4, .9, .9, .9};
    const int y[10] = {1, 1, 0, 0, 0, 0, 1, 1, 1, 0};
    const int g[10] = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1};
    double cnt[10], met[20], rat[20], lo[16], hi[16];
    tk_true(tt, "fair.hand.rc", clab_fair(s, y, g, 10, 2, 0.5, 0, Z, cnt, met, rat,
                                 lo, hi) == CLAB_OK);
    const double c0[5] = {2, 1, 3, 1, 7}, c1[5] = {2, 1, 0, 0, 3};
    for (int k = 0; k < 5; ++k) {
      tk_near(tt, "fair.hand.cnt0", cnt[2 * k], c0[k], 0);
      tk_near(tt, "fair.hand.cnt1", cnt[2 * k + 1], c1[k], 0);
    }
    const double nan = std::nan("");
    // acc tpr tnr fpr fnr ppv npv ppr npos mcc
    const double m0[10] = {5. / 7, 2. / 3, 3. / 4, 1. / 4, 1. / 3,
                           2. / 3, 3. / 4, 3. / 7, 3, 5. / 12};
    const double m1[10] = {2. / 3, 1, 0, 1, 0, 2. / 3, nan, 1, 3, nan};
    for (int k = 0; k < 10; ++k) {
      tk_near(tt, "fair.hand.met0", met[2 * k], m0[k], 1e-15);
      tk_near(tt, "fair.hand.met1", met[2 * k + 1], m1[k], 1e-15);
      tk_near(tt, "fair.hand.rat", rat[2 * k + 1], m1[k] / m0[k], 1e-15);
      tk_near(tt, "fair.hand.rat0", rat[2 * k], 1, 0);
    }
    // all-positive group: the npv interval is undefined, the ppr interval
    // reaches 1 exactly
    tk_true(tt, "fair.hand.npvci", std::isnan(lo[2 * 6 + 1]) && std::isnan(hi[2 * 6 + 1]));
    tk_near(tt, "fair.hand.pprhi", hi[2 * 7 + 1], 1, 0);
    tk_near(tt, "fair.hand.fnrlo", lo[2 * 4 + 1], 0, 0);
    // x = 0 closed form: hi = z^2 / (m + z^2), here m = 2 (fnr of group 1)
    tk_near(tt, "fair.hand.fnrhi", hi[2 * 4 + 1], Z * Z / (2 + Z * Z), 1e-15);
  }
  // 2. Wilson: parity with prop.test, and the defining identity
  // (x - m p)^2 = z^2 m p (1 - p) at each interior end.
  {
    for (int k = 0; k < 6; ++k) {
      // one group, m = 100, x positives: scores 1 for the first x cases
      double s[100];
      int y[100], g[100];
      for (int i = 0; i < 100; ++i) {
        s[i] = i < wx[k] ? 1.0 : 0.0;
        y[i] = 1;
        g[i] = 0;
      }
      double cnt[5], met[10], rat[10], lo[8], hi[8];
      clab_fair(s, y, g, 100, 1, 0.5, 0, Z, cnt, met, rat, lo, hi);
      // ppr = x / 100
      tk_near(tt, "fair.wilson.lo", lo[7], wlo[k], 1e-14);
      tk_near(tt, "fair.wilson.hi", hi[7], whi[k], 1e-14);
      const double x = wx[k], m = 100;
      for (int e = 0; e < 2; ++e) {
        const double p = e == 0 ? lo[7] : hi[7];
        if (p <= 0.0 || p >= 1.0) continue;
        const double l = (x - m * p) * (x - m * p), r = Z * Z * m * p * (1 - p);
        tk_near(tt, "fair.wilson.score", l, r, 1e-12 * r);
      }
    }
  }
  // 3. parity with fairness 1.2.3 on mixed predictions
  {
    static double p[1000];
    static int y[1000], g[1000];
    detsample(p, y, g, 1000);
    double cnt[15], met[30], rat[30];
    tk_true(tt, "fair.parity.rc", clab_fair(p, y, g, 1000, 3, 0.5, 0, Z, cnt, met, rat,
                                   nullptr, nullptr) == CLAB_OK);
    for (int k = 0; k < 10; ++k) {
      for (int j = 0; j < 3; ++j) tk_near(tt, "fair.parity", met[k * 3 + j], ref[k][j], 1e-14);
    }
    // 4. every score above the cutoff: every case is predicted positive, so
    // the sensitivity is 1 and the positives are the group sizes.
    clab_fair(p, y, g, 1000, 3, -1.0, 0, Z, cnt, met, rat, nullptr, nullptr);
    for (int j = 0; j < 3; ++j) {
      tk_near(tt, "fair.allpos.tpr", met[1 * 3 + j], 1, 0);
      tk_near(tt, "fair.allpos.npos", met[8 * 3 + j], cnt[4 * 3 + j], 0);
      tk_near(tt, "fair.allpos.size", cnt[4 * 3 + j], j == 0 ? 334 : 333, 0);
    }
  }
  // 5. seeded simulation: counts against a direct recount, 5 groups
  {
    static double s[5000];
    static int y[5000], g[5000];
    xrng r{0x9E3779B97F4A7C15ULL};
    for (int i = 0; i < 5000; ++i) {
      g[i] = static_cast<int>(r.unif() * 5);
      y[i] = r.unif() < 0.1 + 0.15 * g[i] ? 1 : 0;
      s[i] = r.unif() * 0.5 + 0.4 * y[i];
    }
    double cnt[25], met[50], rat[50], lo[40], hi[40];
    clab_fair(s, y, g, 5000, 5, 0.45, 2, Z, cnt, met, rat, lo, hi);
    for (int j = 0; j < 5; ++j) {
      double tp = 0, fp = 0, tn = 0, fn = 0;
      for (int i = 0; i < 5000; ++i) {
        if (g[i] != j) continue;
        tp += y[i] == 1 && s[i] > 0.45;
        fp += y[i] == 0 && s[i] > 0.45;
        tn += y[i] == 0 && !(s[i] > 0.45);
        fn += y[i] == 1 && !(s[i] > 0.45);
      }
      tk_near(tt, "fair.sim.tp", cnt[j], tp, 0);
      tk_near(tt, "fair.sim.fp", cnt[5 + j], fp, 0);
      tk_near(tt, "fair.sim.tn", cnt[10 + j], tn, 0);
      tk_near(tt, "fair.sim.fn", cnt[15 + j], fn, 0);
      const double mcc = (tp * tn - fp * fn) /
                         std::sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));
      tk_near(tt, "fair.sim.mcc", met[45 + j], mcc, 1e-14);
      tk_near(tt, "fair.sim.base", rat[5 + j], met[5 + j] / met[5 + 2], 0);
      for (int k = 0; k < 8; ++k) {
        tk_true(tt, "fair.sim.ci", lo[k * 5 + j] <= met[k * 5 + j] &&
                              met[k * 5 + j] <= hi[k * 5 + j]);
      }
    }
  }
  // 6. guards
  {
    const double s[2] = {0.1, 0.9};
    const int y[2] = {0, 1}, g[2] = {0, 1}, y2[2] = {0, 2}, g2[2] = {0, 2};
    double cnt[10], met[20], rat[20], lo[16], hi[16];
    tk_true(tt, "fair.guard.n", clab_fair(s, y, g, 0, 2, .5, 0, Z, cnt, met, rat, lo, hi) == CLAB_EARG);
    tk_true(tt, "fair.guard.base", clab_fair(s, y, g, 2, 2, .5, 2, Z, cnt, met, rat, lo, hi) == CLAB_EARG);
    tk_true(tt, "fair.guard.grp", clab_fair(s, y, g2, 2, 2, .5, 0, Z, cnt, met, rat, lo, hi) == CLAB_EARG);
    tk_true(tt, "fair.guard.y", clab_fair(s, y2, g, 2, 2, .5, 0, Z, cnt, met, rat, lo, hi) == CLAB_EDOM);
    tk_true(tt, "fair.guard.cut", clab_fair(s, y, g, 2, 2, std::nan(""), 0, Z, cnt, met, rat, lo, hi) == CLAB_EDOM);
    tk_true(tt, "fair.guard.z", clab_fair(s, y, g, 2, 2, .5, 0, 0, cnt, met, rat, lo, hi) == CLAB_EDOM);
    const double sn[2] = {std::nan(""), 0.9};
    tk_true(tt, "fair.guard.nan", clab_fair(sn, y, g, 2, 2, .5, 0, Z, cnt, met, rat, lo, hi) == CLAB_EDOM);
  }
}


}  // namespace

void t_fair(tk &t, int) { run(t); }

#ifdef CLAB_MAIN
#include <chrono>
#include <cstdio>

// timing: n = 1e5, 4 groups (fairness::equal_odds took 0.106 s, median
// of 5 on zeus)
void b_fair() {
  const int n = 100000;
  static double s[n];
  static int y[n], g[n];
  xrng r{12345};
  for (int i = 0; i < n; ++i) {
    y[i] = r.unif() < 0.4 ? 1 : 0;
    g[i] = static_cast<int>(r.unif() * 4);
    s[i] = r.unif();
  }
  double cnt[20], met[40], rat[40], lo[32], hi[32], t[5];
  for (int k = 0; k < 5; ++k) {
    const auto a = std::chrono::steady_clock::now();
    clab_fair(s, y, g, n, 4, 0.5, 0, Z, cnt, met, rat, lo, hi);
    t[k] = std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count();
  }
  std::printf("BENCH clab_fair n=1e5 g=4: median %.6f s (k=5)\n", tk_med5(t));
}
#endif
