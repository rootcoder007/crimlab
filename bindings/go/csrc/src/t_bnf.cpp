// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for clab_bdig and clab_bnft (src/bnf.cpp). Part of the
// self-test (selftest.cpp).

#include <cmath>
#include <cstdint>

#include "crimlab.h"
#include "tchk.h"

namespace {

// Reference values, dev/ref/cpp/bnf.R, R 4.5.0. With i = 0..999:
// BenfordTests 1.2.0 chisq.benftest(d1 + 0.5), d1 = 1 + ((i * 7 + (i^2 %%
// 13)) %% 9): statistic 400.72808136642323 (its p-value prints as 0).
const double ref_chi1 = 400.72808136642323;
// chisq.benftest(rep(1:9, nb) + 0.5), nb = round(10000 log10(1 + 1/d)):
// statistic, p-value
const double ref_chib[2] = {0.0013828772800778026, 0.99999999999999045};
// benford.analysis 0.1.5 benford(d2 + 0.5, number.of.digits = 2, sign =
// "both"), d2 = 10 + ((i * 37 + (i^2 %% 31)) %% 90): chisq statistic, df,
// p-value (pchisq upper tail), and MAD(); BenfordTests chisq.benftest(x2,
// digits = 2) gives the same statistic.
const double ref_chi2[3] = {501.68220140807284, 89, 8.2333645048997026e-59};
const double ref_mad2 = 0.0062786528295698191;

struct xrng {
  std::uint64_t s;
  double unif() {
    s ^= s >> 12;
    s ^= s << 25;
    s ^= s >> 27;
    return static_cast<double>((s * 2685821657736338717ULL) >> 11) * 0x1p-53;
  }
};

std::int64_t lead(std::int64_t k, int nd) {
  std::int64_t hi = 1;
  for (int j = 0; j < nd; ++j) hi *= 10;
  while (k >= hi) k /= 10;
  while (k < hi / 10) k *= 10;
  return k;
}

void run(tk &tt) {
  int cnt[4];
  // 1. the first two digits of 0.29 0.57 0.58 1.15 0.07 are 29 57 58 11 70
  // (the decimals as recorded, by both routes).
  {
    const double x[5] = {0.29, 0.57, 0.58, 1.15, 0.07};
    const int want[5] = {29, 57, 58, 11, 70};
    int d[5];
    for (int dp = -1; dp <= 2; dp += 3) {
      tk_true(tt, "bnf.dig.rc", clab_bdig(x, 5, 2, dp, 0, d, cnt) == CLAB_OK);
      for (int i = 0; i < 5; ++i) tk_true(tt, "bnf.dig.dec", d[i] == want[i]);
    }
  }
  // 2. exhaustive: every 2-decimal value 0.10 .. 999.99 as k / 100 (the
  // nearest double), and the doubles one ulp either side with dp = 2; the
  // first one to four digits of k by integer arithmetic.
  {
    int bad = 0, naive = 0;
    for (int k = 10; k <= 99999; ++k) {
      const double x = k / 100.0;
      const double xs[3] = {x, std::nextafter(x, 0.0), std::nextafter(x, 1e9)};
      for (int nd = 1; nd <= 4; ++nd) {
        int d[3], c[4];
        clab_bdig(xs, 3, nd, 2, 1, d, c);
        bad += d[0] != lead(k, nd) || d[1] != lead(k, nd) || d[2] != lead(k, nd);
        clab_bdig(xs, 1, nd, -1, 1, d, c);
        bad += d[0] != lead(k, nd);
      }
      // the truncation formula 10^(1 - floor(log10 x)) x, which is not exact
      naive += static_cast<std::int64_t>(std::trunc(std::pow(10.0, 1 - std::floor(std::log10(x))) * x)) != lead(k, 2);
    }
    tk_true(tt, "bnf.dig.exhaustive", bad == 0);
    tk_true(tt, "bnf.dig.naive_fails", naive > 0);
  }
  // 3. shortest-decimal route on random decimals k / 10^p (k with 1..15
  // digits, no shorter decimal rounds to the same double) and 17-digit
  // doubles
  {
    xrng r{31};
    int bad = 0, oor = 0;
    for (int t = 0; t < 100000; ++t) {
      const int s = 1 + static_cast<int>(r.unif() * 15);
      double hi = 1;
      for (int j = 0; j < s; ++j) hi *= 10;
      std::int64_t k = static_cast<std::int64_t>(hi / 10 + r.unif() * (hi - hi / 10));
      const int p = static_cast<int>(r.unif() * 21);
      double pw = 1;
      for (int j = 0; j < p; ++j) pw *= 10;
      const double x = static_cast<double>(k) / pw;
      int d, c[4];
      const int nd = 1 + static_cast<int>(r.unif() * 9);
      clab_bdig(&x, 1, nd, -1, 0, &d, c);
      // out of the exact range only when the scaling 10^q, q = e - nd + 1,
      // passes 10^-22 (e = s - 1 - p, the decimal exponent of x)
      const int q = (s - 1 - p) - nd + 1;
      if (c[3] == 1) {
        bad += d != 0 || q >= -22;
        ++oor;
      } else {
        bad += d != lead(k, nd);
      }
    }
    tk_true(tt, "bnf.short.random", bad == 0);
    tk_true(tt, "bnf.short.oor", oor > 0 && oor < 5000);
    const double x[3] = {0.1 + 0.2, 2.0 / 3.0, 1e21 / 3.0};
    const int want[3] = {300000000, 666666666, 333333333};
    int d[3];
    clab_bdig(x, 3, 9, -1, 0, d, cnt);
    for (int i = 0; i < 3; ++i) tk_true(tt, "bnf.short.17", d[i] == want[i]);
  }
  // 4. zeros, non-finite values and signs are set aside and counted
  {
    const double x[7] = {0.0, -0.0, std::nan(""), HUGE_VAL, -3.5, 4.5, 0.001};
    int d[7];
    clab_bdig(x, 7, 1, 2, 1, d, cnt);
    tk_true(tt, "bnf.cnt", cnt[0] == 3 && cnt[1] == 2 && cnt[2] == 1 && cnt[3] == 0);
    tk_true(tt, "bnf.cnt.dig", d[5] == 4 && d[4] == 0 && d[6] == 0);
    clab_bdig(x, 7, 1, -1, 0, d, cnt);
    tk_true(tt, "bnf.cnt.both", cnt[0] == 2 && cnt[2] == 0 && d[4] == 3 && d[6] == 1);
  }
  // 5. tests: parity, and zeros left out of the count
  {
    static int dg[30000];
    double obs[90], work[90], out[6];
    for (int i = 0; i < 1000; ++i) dg[i] = 1 + ((i * 7 + (i * i) % 13) % 9);
    tk_true(tt, "bnf.t.rc", clab_bnft(dg, 1000, 1, 0, 1, obs, work, out) == CLAB_OK);
    tk_near(tt, "bnf.t.chi1", out[1], ref_chi1, 1e-13);
    // df = 8 closed form e^-y (1 + y + y^2/2 + y^3/6), y = chi / 2: the
    // reference package prints 1 - pchisq = 0 here
    {
      const double y = ref_chi1 / 2;
      tk_near(tt, "bnf.t.p1", out[3] / (std::exp(-y) * (1 + y + y * y / 2 + y * y * y / 6)), 1, 1e-12);
    }
    // the near-Benford sample, then with 2000 zeros (digit 0 = excluded)
    int m = 0;
    for (int d = 1; d <= 9; ++d) {
      const int nb = static_cast<int>(std::floor(10000 * std::log10(1 + 1.0 / d) + 0.5));
      for (int j = 0; j < nb; ++j) dg[m++] = d;
    }
    clab_bnft(dg, m, 1, 0, 1, obs, work, out);
    // the statistic is a sum of squares of differences of order 1e-5:
    // relative agreement to 1e-9 is what binary64 allows either side
    tk_near(tt, "bnf.t.chib", out[1] / ref_chib[0], 1, 1e-9);
    tk_near(tt, "bnf.t.pb", out[3], ref_chib[1], 1e-15);
    for (int j = 0; j < 2000; ++j) dg[m + j] = 0;
    clab_bnft(dg, m + 2000, 1, 0, 1, obs, work, out);
    // zeros are set aside (digit 0): the statistic and p are those of the
    // sample without them
    tk_near(tt, "bnf.t.zeros.m", out[0], m, 0);
    tk_near(tt, "bnf.t.zeros.chi", out[1] / ref_chib[0], 1, 1e-9);
    tk_near(tt, "bnf.t.zeros.p", out[3], ref_chib[1], 1e-15);
    // two digits
    for (int i = 0; i < 1000; ++i) dg[i] = 10 + ((i * 37 + (i * i) % 31) % 90);
    clab_bnft(dg, 1000, 2, 0, 1, obs, work, out);
    tk_near(tt, "bnf.t.chi2", out[1], ref_chi2[0], 1e-13);
    tk_near(tt, "bnf.t.df2", out[2], ref_chi2[1], 0);
    tk_near(tt, "bnf.t.p2", out[3] / ref_chi2[2], 1, 1e-10);
    tk_near(tt, "bnf.t.mad2", out[5], ref_mad2, 1e-15);
  }
  // 6. Monte Carlo: never 0, calibrated under the null (400 Benford
  // samples of 50 drawn here by a linear search of the cumulative
  // probabilities, a sampler independent of the kernel's; 199 simulations
  // each; the rejection rate at 0.05 has sd 0.011, tolerance 4 sd), and
  // close to the asymptotic p for n = 2000 (sd sqrt(p (1 - p) / 999)).
  {
    xrng r{8};
    double cum[9], acc = 0;
    for (int d = 1; d <= 9; ++d) cum[d - 1] = acc += std::log10(1 + 1.0 / d);
    int rej = 0;
    static int dg[2000];
    double obs[9], work[9], out[6];
    for (int rep = 0; rep < 400; ++rep) {
      for (int i = 0; i < 50; ++i) {
        const double u = r.unif() * acc;
        int d = 0;
        while (d < 8 && u >= cum[d]) ++d;
        dg[i] = d + 1;
      }
      clab_bnft(dg, 50, 1, 199, static_cast<unsigned>(rep), obs, work, out);
      tk_true(tt, "bnf.mc.pos", out[4] >= 1.0 / 200);
      rej += out[4] <= 0.05;
    }
    tk_true(tt, "bnf.mc.size", rej >= 0.006 * 400 && rej <= 0.094 * 400);
    for (int i = 0; i < 2000; ++i) dg[i] = 1 + (i * 7919 % 2000 < 650 ? 0 : (i % 9));
    clab_bnft(dg, 2000, 1, 999, 5, obs, work, out);
    const double p = out[3], sd = std::sqrt(p * (1 - p) / 999);
    tk_true(tt, "bnf.mc.asym", std::fabs(out[4] - p) <= 4 * sd + 1.0 / 1000);
    // a gross departure: p is 1 / (B + 1), not 0
    for (int i = 0; i < 2000; ++i) dg[i] = 9;
    clab_bnft(dg, 2000, 1, 99, 5, obs, work, out);
    tk_near(tt, "bnf.mc.min", out[4], 1.0 / 100, 0);
  }
  // 7. guards and the empty sample
  {
    int d[2] = {0, 0};
    double obs[9], work[9], out[6];
    tk_true(tt, "bnf.empty", clab_bnft(d, 2, 1, 0, 1, obs, work, out) == CLAB_OK && out[0] == 0 &&
                        std::isnan(out[1]));
    d[0] = 12;
    tk_true(tt, "bnf.guard.dig", clab_bnft(d, 2, 1, 0, 1, obs, work, out) == CLAB_EDOM);
    tk_true(tt, "bnf.guard.nd", clab_bdig(nullptr, 0, 10, 0, 0, d, cnt) == CLAB_EARG);
    tk_true(tt, "bnf.guard.dp", clab_bdig(nullptr, 0, 1, 23, 0, d, cnt) == CLAB_EARG);
    const double big = 1e15;
    clab_bdig(&big, 1, 1, 2, 0, d, cnt);
    tk_true(tt, "bnf.range", cnt[3] == 1 && d[0] == 0);
  }
}


}  // namespace

void t_bnf(tk &t, int) { run(t); }

#ifdef CLAB_MAIN
#include <chrono>
#include <cstdio>

// timing: n = 1e6, first two digits
// (benford.analysis::benford took 1.622 s, BenfordTests::chisq.benftest
// 0.639 s, medians of 3 on zeus)
void b_bnf() {
  const int n = 1000000;
  static double x[n];
  static int d[n];
  double obs[90], work[90], out[6], t[5], t2[5];
  int cnt[4];
  xrng r{5};
  for (int i = 0; i < n; ++i) x[i] = std::floor(std::pow(10.0, 6 * r.unif()) * 100) / 100;
  for (int k = 0; k < 5; ++k) {
    auto a = std::chrono::steady_clock::now();
    clab_bdig(x, n, 2, -1, 0, d, cnt);
    clab_bnft(d, n, 2, 0, 1, obs, work, out);
    t[k] = std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count();
    a = std::chrono::steady_clock::now();
    clab_bdig(x, n, 2, 2, 0, d, cnt);
    clab_bnft(d, n, 2, 0, 1, obs, work, out);
    t2[k] = std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count();
  }
  std::printf("BENCH clab_bdig+bnft n=1e6 2 digits, shortest decimal: median %.6f s (k=5)\n", tk_med5(t));
  std::printf("BENCH clab_bdig+bnft n=1e6 2 digits, dp = 2: median %.6f s (k=5)\n", tk_med5(t2));
}
#endif
