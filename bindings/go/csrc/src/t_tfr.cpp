// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for clab_tfrp and clab_tfrw (src/tfr.cpp). Part of the
// self-test (selftest.cpp).

#include <cmath>

#include "crimlab.h"
#include "tchk.h"

namespace {

// tfer 1.3, R 4.5.0, dev/ref/cpp/tfr.R: set.seed(2026);
// tprob(transfer(N = 1e6))[1:81], the probabilities of 0..80 recovered
// fragments under the package defaults, by its own Monte Carlo.
const double tf_ref[81] = {0.097566, 0.061145999999999999, 0.051212000000000001, 0.045168, 0.040996999999999999, 0.038432000000000001, 0.035656, 0.033175999999999997, 0.03108, 0.028909000000000001, 0.027015000000000001, 0.025871999999999999, 0.024117, 0.022884000000000002, 0.021977, 0.020486000000000001, 0.019237000000000001, 0.018626, 0.017527999999999998, 0.016875999999999999, 0.016035000000000001, 0.015539000000000001, 0.014812000000000001, 0.014174000000000001, 0.013679999999999999, 0.013011999999999999, 0.012241999999999999, 0.011894, 0.011261999999999999, 0.010824, 0.010305999999999999, 0.0097529999999999995, 0.0095469999999999999, 0.0090340000000000004, 0.0085760000000000003, 0.0083260000000000001, 0.0079590000000000008, 0.0077099999999999998, 0.0071060000000000003, 0.0068009999999999998, 0.0063759999999999997, 0.0060850000000000001, 0.0057029999999999997, 0.0053969999999999999, 0.0051549999999999999, 0.0049649999999999998, 0.004607, 0.0045079999999999999, 0.0041200000000000004, 0.0038289999999999999, 0.00362, 0.0034510000000000001, 0.0032499999999999999, 0.0030590000000000001, 0.0028860000000000001, 0.0026809999999999998, 0.002457, 0.0023140000000000001, 0.0021840000000000002, 0.0020460000000000001, 0.0018979999999999999, 0.0018240000000000001, 0.001585, 0.001469, 0.0014790000000000001, 0.0013569999999999999, 0.0012570000000000001, 0.0011039999999999999, 0.0011119999999999999, 0.00097599999999999998, 0.00090499999999999999, 0.00088000000000000003, 0.00080699999999999999, 0.00071500000000000003, 0.00066600000000000003, 0.00056999999999999998, 0.00052800000000000004, 0.00049600000000000002, 0.00045100000000000001, 0.00043899999999999999, 0.00038699999999999997};

// the package defaults, in clab_tfer's parameter order
const double DEF[17] = {120, 60, 0.5, 0.05, 0.8, 0.9, 0.1, 0.15, 0.45, 0.7,
                        0.05, 0.1, 0.5, 0.7, 1, 2, 0.5};

double work[400000];

// the defaults with every rate interval a point and one hour
void fixed(double *p) {
  for (int k = 0; k < 17; ++k) p[k] = DEF[k];
  p[5] = p[4];
  p[7] = p[6];
  p[9] = p[8];
  p[11] = p[10];
  p[13] = p[12];
  p[14] = p[15] = 1;
}

// Poisson pmf by logs
double pois(int y, double v) { return std::exp(-v + y * std::log(v) - std::lgamma(y + 1.0)); }

// standard normal tail Phi(-z)
double ptail(double z) { return 0.5 * std::erfc(z / std::sqrt(2.0)); }

void run(tk &tt) {
  int len;
  double prob[401], se[401], cut[2];
  // 1. no random effects (sd 0, no distance effect, fixed rates, one
  // hour): the count is Poisson(lam s), s = R ((1 - Q)(1 - a) + Q (1 - b))
  {
    double p[17];
    for (int k = 0; k < 17; ++k) p[k] = DEF[k];
    p[1] = 0;
    p[5] = p[4];
    p[7] = p[6];
    p[9] = p[8];
    p[11] = p[10];
    p[13] = p[12];
    p[14] = p[15] = 1;
    tk_true(tt, "tfr.pois.work", clab_tfrw(p, 0, 0, 60, 0.02, 256, &len) == CLAB_OK && len <= 400000);
    tk_true(tt, "tfr.pois.rc", clab_tfrp(p, 0, 0, 60, 10, 1, 0.02, 256, prob, se, cut, work) == CLAB_OK);
    const double s = 0.5 * (0.95 * 0.2 + 0.05 * 0.9);
    for (int y = 0; y <= 60; ++y) {
      tk_near(tt, "tfr.pois", prob[y], pois(y, 120 * s), 1e-14);
      tk_true(tt, "tfr.pois.se", se[y] == 0);
    }
  }
  // 2. normal transfer amount (folded), fixed rates, no distance effect:
  // against a composite Simpson rule with 20,000 intervals on [0, 840]
  // in this file (a different rule, 100 times finer; 200,000 intervals in
  // the extended run give 3.3e-12)
  {
    double p[17];
    for (int k = 0; k < 17; ++k) p[k] = DEF[k];
    p[5] = p[4];
    p[7] = p[6];
    p[9] = p[8];
    p[11] = p[10];
    p[13] = p[12];
    p[14] = p[15] = 1;
    clab_tfrp(p, 0, 0, 80, 10, 1, 0.02, 256, prob, se, cut, work);
    const double s = 0.5 * (0.95 * 0.2 + 0.05 * 0.9);
    double worst = 0;
    for (int y = 0; y <= 80; y += 5) {
      const int n = 20000;
      const double h = 840.0 / n;
      double acc = 0;
      for (int i = 0; i <= n; ++i) {
        const double x = i * h, z1 = (x - 120) / 60, z2 = (x + 120) / 60;
        const double f = (x == 0 ? (y == 0 ? 1.0 : 0.0) : pois(y, x * s)) *
                         (std::exp(-0.5 * z1 * z1) + std::exp(-0.5 * z2 * z2)) / (60 * 2.5066282746310002);
        acc += f * (i == 0 || i == n ? 1 : (i % 2 ? 4 : 2));
      }
      acc *= h / 3;
      const double rel = std::fabs(prob[y] - acc) / acc;
      worst = rel > worst ? rel : worst;
    }
    tk_true(tt, "tfr.normal", worst < 1e-8);
  }
  // 3. with the distance effect (D ~ Gamma(0.5)), fixed rates: the mean in
  // closed form and the total (the comparison with a 2-D Simpson rule and
  // the self-convergence check are in quad(), full size in the ext mode)
  {
    double p[17];
    fixed(p);
    const double s = 0.5 * (0.95 * 0.2 + 0.05 * 0.9);
    double cf[2];
    // the mean in closed form: E|X| e (1 + 1/d)^{-d} s (counts to 400:
    // the law has a long right tail)
    static double pm[401], sm[401];
    clab_tfrp(p, 1, 0, 400, 10, 1, 0.02, 256, pm, sm, cf, work);
    double mean = 0, tot = 0;
    for (int y = 0; y <= 400; ++y) {
      mean += y * pm[y];
      tot += pm[y];
    }
    tk_near(tt, "tfr.total.fixed", tot, 1, 1e-12);
    const double ex = 60 * std::sqrt(2 / 3.141592653589793) * std::exp(-2.0) + 120 * (1 - 2 * ptail(2.0));
    const double em = ex * std::exp(1.0) * std::pow(3.0, -0.5) * s;
    tk_near(tt, "tfr.mean.fixed", mean / em, 1, 1e-9);
    tk_true(tt, "tfr.cut", cf[1] > 0.99 && cf[1] < 1 + 1e-9);
  }
  // 3b. a fixed transfer amount with the distance effect (the w rule
  // alone): mean |lam| e (1 + 1/d)^{-d} s in closed form; and intensities
  // above 600 (the log-space branch): Poisson(lam s) exactly
  {
    double p[17];
    fixed(p);
    p[1] = 0;
    const double s = 0.5 * (0.95 * 0.2 + 0.05 * 0.9);
    static double pm[401], sm[401];
    double cf[2];
    clab_tfrp(p, 1, 0, 400, 1, 1, 0.02, 64, pm, sm, cf, work);
    double mean = 0;
    for (int y = 0; y <= 400; ++y) mean += y * pm[y];
    tk_near(tt, "tfr.mean.sd0", mean / (120 * std::exp(1.0) * std::pow(3.0, -0.5) * s), 1, 1e-9);
    p[0] = 6000;
    static double pb[1001], sb[1001];
    clab_tfrp(p, 0, 0, 1000, 1, 1, 0.02, 64, pb, sb, cf, work);
    tk_near(tt, "tfr.big", pb[700] / pois(700, 6000 * s), 1, 1e-10);
  }
  // 4. the default model: mean against the closed form
  // E[Y] = E|X| e (1 + 1/d)^{-d} E[R] sum_h P(h) ((1-Q) E[A_h] + Q E[B_h]),
  // the folded-normal mean E|X| = sd sqrt(2/pi) e^{-mu^2/(2 sd^2)} +
  // mu (1 - 2 Phi(-mu/sd)); and parity with tfer's 10^6 draws within 4
  // standard errors of both Monte Carlo estimates
  {
    tk_true(tt, "tfr.def.work", clab_tfrw(DEF, 1, 0, 400, 0.02, 256, &len) == CLAB_OK && len <= 400000);
    tk_true(tt, "tfr.def.rc", clab_tfrp(DEF, 1, 0, 400, 4, 7, 0.05, 64, prob, se, cut, work) == CLAB_OK);
    double mean = 0, tot = 0, msd = 0;
    for (int y = 0; y <= 400; ++y) {
      mean += y * prob[y];
      tot += prob[y];
      msd += y * se[y];
    }
    const double ex = 60 * std::sqrt(2 / 3.141592653589793) * std::exp(-2.0) + 120 * (1 - 2 * ptail(2.0));
    const double em = ex * std::exp(1.0) * std::pow(3.0, -0.5);
    // hours: NB(1.5, 0.5), h = max(t, 1); E[A_h] = 0.15 0.425^{h-1}, E[B_h] = 0.875 0.925^{h-1}
    double es = 0, pk = std::pow(0.5, 1.5);
    for (int t = 0; t < 400; ++t) {
      const int h = t < 1 ? 1 : t;
      es += pk * (0.95 * 0.15 * std::pow(0.425, h - 1) + 0.05 * 0.875 * std::pow(0.925, h - 1));
      pk *= (t + 1.5) * 0.5 / (t + 1);
    }
    es *= 0.6;
    tk_true(tt, "tfr.def.mean", std::fabs(mean - em * es) <= 4 * msd);
    tk_true(tt, "tfr.def.total", std::fabs(tot - 1) < 1e-6);
    int bad = 0;
    for (int y = 0; y <= 80; ++y) {
      const double sd = std::sqrt(tf_ref[y] * (1 - tf_ref[y]) / 1e6 + se[y] * se[y]);
      bad += std::fabs(prob[y] - tf_ref[y]) > 4 * sd + 1e-6;
    }
    tk_true(tt, "tfr.parity", bad <= 1);
    // every count has its own positive probability (no count above the
    // largest simulated one gets 0)
    bool pos = true;
    for (int y = 0; y <= 400; ++y) pos = pos && prob[y] > 0;
    tk_true(tt, "tfr.positive", pos && prob[209] > 0 && prob[209] < 1e-6);
  }
  // 5. reproducible, and the other time distributions
  {
    static double p2[401], s2[401];
    clab_tfrp(DEF, 1, 0, 100, 3, 11, 0.05, 64, prob, se, cut, work);
    clab_tfrp(DEF, 1, 0, 100, 3, 11, 0.05, 64, p2, s2, cut, work);
    bool same = true;
    for (int y = 0; y <= 100; ++y) same = same && prob[y] == p2[y] && se[y] == s2[y];
    tk_true(tt, "tfr.repro", same);
    for (int td = 1; td <= 2; ++td) {
      double p[17];
      for (int k = 0; k < 17; ++k) p[k] = DEF[k];
      p[14] = 1;
      p[15] = 4;
      tk_true(tt, "tfr.tdist.rc", clab_tfrp(p, 1, td, 200, 2, 3, 0.05, 64, prob, se, cut, work) == CLAB_OK);
      double tot = 0;
      for (int y = 0; y <= 200; ++y) tot += prob[y];
      tk_true(tt, "tfr.tdist.total", std::fabs(tot - 1) < 1e-6 && cut[0] == 0);
    }
  }
  // 6. guards (tfer 1.3 checks no parameter: Q = -1 ends in "NAs produced")
  {
    double p[17];
    for (int k = 0; k < 17; ++k) p[k] = DEF[k];
    p[3] = -1;
    tk_true(tt, "tfr.guard.Q", clab_tfrp(p, 1, 0, 10, 10, 1, 0.02, 256, prob, se, cut, work) == CLAB_EDOM);
    p[3] = 0.05;
    p[4] = 0.95;
    tk_true(tt, "tfr.guard.lu", clab_tfrp(p, 1, 0, 10, 10, 1, 0.02, 256, prob, se, cut, work) == CLAB_EDOM);
    p[4] = 0.8;
    p[16] = 0;
    tk_true(tt, "tfr.guard.r", clab_tfrp(p, 1, 0, 10, 10, 1, 0.02, 256, prob, se, cut, work) == CLAB_EDOM);
    p[16] = 0.5;
    p[14] = p[15] = 1.5;
    tk_true(tt, "tfr.guard.t", clab_tfrp(p, 1, 0, 10, 10, 1, 0.02, 256, prob, se, cut, work) == CLAB_EDOM);
    tk_true(tt, "tfr.guard.nsim", clab_tfrp(DEF, 1, 0, 10, 0, 1, 0.02, 256, prob, se, cut, work) == CLAB_EARG);
    tk_true(tt, "tfr.guard.d", clab_tfrp(DEF, 2, 0, 10, 10, 1, 0.02, 256, prob, se, cut, work) == CLAB_EARG);
  }
}


// The quadrature against independent rules: full size in the "ext" mode
// (big = 1), reduced in the default self-test.
void quad(tk &tt, int big) {
  static double prob[81], se[81], cut[2];
  double p[17];
  fixed(p);
  // full grid in the ext mode, a coarse one by default
  clab_tfrp(p, 1, 0, 80, 10, 1, big ? 0.02 : 0.05, big ? 256 : 64, prob, se, cut, work);
  const double s = 0.5 * (0.95 * 0.2 + 0.05 * 0.9);
  // 2-D Simpson in x and u = sqrt(D) (density 2 e^{-u^2} / Gamma(1/2)),
  // 2,000 x 4,000 intervals (200 x 400 by default), for a few counts
  const int ys[4] = {0, 3, 20, 60};
  double worst = 0;
  for (int q = 0; q < 4; ++q) {
    const int y = ys[q], nx = big ? 2000 : 200, nu = big ? 4000 : 400;
    const double hx = 840.0 / nx, hu = 7.0 / nu;
    double acc = 0;
    for (int i = 0; i <= nx; ++i) {
      const double x = i * hx, z1 = (x - 120) / 60, z2 = (x + 120) / 60;
      const double fx = (std::exp(-0.5 * z1 * z1) + std::exp(-0.5 * z2 * z2)) / (60 * 2.5066282746310002);
      const double cx = i == 0 || i == nx ? 1 : (i % 2 ? 4 : 2);
      double inner = 0;
      for (int j = 0; j <= nu; ++j) {
        const double u = j * hu, dd = u * u;
        const double v = x * std::exp(1 - dd / 0.5) * s;
        const double f = (v == 0 ? (y == 0 ? 1.0 : 0.0) : pois(y, v)) * 2 * std::exp(-dd) / std::sqrt(3.141592653589793);
        inner += f * (j == 0 || j == nu ? 1 : (j % 2 ? 4 : 2));
      }
      acc += cx * fx * inner * hu / 3;
    }
    acc *= hx / 3;
    const double rel = std::fabs(prob[y] - acc) / acc;
    worst = rel > worst ? rel : worst;
  }
  tk_true(tt, "tfr.gamma", worst < (big ? 1e-6 : 1e-3));
  // self-convergence: halving the grid step and doubling the panels moves
  // no probability by more than 1e-9 relative, so the difference above is
  // the Simpson reference's own error
  static double pf[81], sf[81];
  double cf[2];
  clab_tfrp(p, 1, 0, 80, 10, 1, big ? 0.01 : 0.04, big ? 512 : 96, pf, sf, cf, work);
  double conv = 0;
  for (int y = 0; y <= 80; ++y)
    if (pf[y] > 1e-300) conv = std::fabs(prob[y] - pf[y]) / pf[y] > conv ? std::fabs(prob[y] - pf[y]) / pf[y] : conv;
  tk_true(tt, "tfr.conv", conv < (big ? 1e-9 : 1e-5));
  // the default model with 40 rate draws (4 by default): parity with tfer
  static double pd[401], sd2[401];
  clab_tfrp(DEF, 1, 0, 400, big ? 40 : 4, 7, big ? 0.02 : 0.05, big ? 256 : 64, pd, sd2, cf, work);
  int bad = 0;
  for (int y = 0; y <= 80; ++y) {
    const double sd = std::sqrt(tf_ref[y] * (1 - tf_ref[y]) / 1e6 + sd2[y] * sd2[y]);
    bad += std::fabs(pd[y] - tf_ref[y]) > 4 * sd + 1e-6;
  }
  tk_true(tt, "tfr.xparity", bad <= 1);
}

}  // namespace

void t_tfr(tk &t, int big) {
  run(t);
  quad(t, big);
}

#ifdef CLAB_MAIN
#include <chrono>
#include <cstdio>


// tfer's standard case: transfer(N = 1e5) 0.681 s (median of 5, arm64 Linux)
void b_tfr() {
  double prob[201], se[201], cut[2], t[5];
  for (int k = 0; k < 5; ++k) {
    const auto a = std::chrono::steady_clock::now();
    clab_tfrp(DEF, 1, 0, 200, 200, 1, 0.02, 256, prob, se, cut, work);
    t[k] = std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count();
  }
  double mx = 0;
  for (int y = 0; y <= 200; ++y) mx = se[y] > mx ? se[y] : mx;
  std::printf("BENCH clab_tfrp defaults, counts 0..200, 200 draws, grid 0.02, 256 panels: median %.3f s (k=5), max se %.2g\n",
              tk_med5(t), mx);
}
#endif
