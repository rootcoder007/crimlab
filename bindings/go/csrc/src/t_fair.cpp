// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for clab_fair and the other kernels of src/fair.cpp
// (clab_fairdf, clab_tmpaud, clab_lipaud, clab_nsbnd, clab_impfpr,
// clab_rnkres). Part of the self-test (selftest.cpp).

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


// Values by hand or from R 4.6.1 (prop.test(x, m, correct = FALSE) for the
// Wilson intervals; the formulas of Newcombe 1998 and Katz et al. 1978
// evaluated on them in R).
void run2(tk &tt) {
  // 1. fairdf by hand. Group 0: decisions 1 1 1 0, outcomes 1 1 0 0;
  // group 1: decisions 0 0 1 0 0, outcomes 1 0 1 0 0. Shares 3/4 and 1/5;
  // true positive rates 1 and 1/2, false positive rates 1/2 and 0.
  {
    const int pr[9] = {1, 1, 1, 0, 0, 0, 1, 0, 0};
    const int y[9] = {1, 1, 0, 0, 1, 0, 1, 0, 0};
    const int g[9] = {0, 0, 0, 0, 1, 1, 1, 1, 1};
    double cnt[12], met[46], sum[11], wk[6 + 400];
    tk_true(tt, "fairdf.hand.rc",
            clab_fairdf(pr, y, g, 9, 2, -1, Z, 0, 0, cnt, met, sum, wk) == CLAB_OK);
    // the base: the larger share, group 0
    tk_near(tt, "fairdf.hand.base", sum[10], 0, 0);
    const double c0[6] = {4, 3, 2, 2, 2, 1}, c1[6] = {5, 1, 2, 1, 3, 0};
    for (int k = 0; k < 6; ++k) {
      tk_near(tt, "fairdf.hand.cnt0", cnt[2 * k], c0[k], 0);
      tk_near(tt, "fairdf.hand.cnt1", cnt[2 * k + 1], c1[k], 0);
    }
    auto m = [&](int k, int j) { return met[2 * k + j]; };
    tk_near(tt, "fairdf.hand.ppos", m(2, 1), 0.2, 1e-16);
    tk_near(tt, "fairdf.hand.dp", m(5, 1), 0.2 - 0.75, 1e-16);
    tk_near(tt, "fairdf.hand.di", m(8, 1), 0.2 / 0.75, 1e-16);
    tk_near(tt, "fairdf.hand.wlo", m(3, 1), 0.036224108632430141, 1e-14);
    tk_near(tt, "fairdf.hand.whi", m(4, 0), 0.95441273919029945, 1e-14);
    tk_near(tt, "fairdf.hand.nlo", m(6, 1), -0.81192959079974802, 1e-14);
    tk_near(tt, "fairdf.hand.nhi", m(7, 1), 0.068137204970365128, 1e-14);
    tk_near(tt, "fairdf.hand.klo", m(9, 1), 0.04226297146417294, 1e-14);
    tk_near(tt, "fairdf.hand.khi", m(10, 1), 1.6825866390249733, 1e-14);
    tk_near(tt, "fairdf.hand.tpr", m(11, 1), 0.5, 0);
    tk_near(tt, "fairdf.hand.fpr", m(12, 0), 0.5, 0);
    tk_near(tt, "fairdf.hand.dtpr", m(13, 1), -0.5, 0);
    tk_near(tt, "fairdf.hand.dfpr", m(16, 1), -0.5, 0);
    tk_near(tt, "fairdf.hand.tlo", m(14, 1), -0.90546879426576932, 1e-14);
    tk_near(tt, "fairdf.hand.thi", m(15, 1), 0.27257278511318139, 1e-14);
    tk_near(tt, "fairdf.hand.flo", m(17, 1), -0.90546879426576921, 1e-14);
    tk_near(tt, "fairdf.hand.fhi", m(18, 1), 0.19259213090610783, 1e-14);
    tk_near(tt, "fairdf.hand.aod", m(19, 1), -0.5, 0);
    // MOVER: the sum's interval from the two Newcombe intervals, halved
    const double lo = -1 - std::sqrt(std::pow(-0.5 + 0.90546879426576932, 2) +
                                     std::pow(-0.5 + 0.90546879426576921, 2));
    tk_near(tt, "fairdf.hand.aodlo", m(20, 1), 0.5 * lo, 1e-14);
    // equal gaps in size: the true positive one
    tk_near(tt, "fairdf.hand.eo", m(22, 1), -0.5, 0);
    tk_true(tt, "fairdf.hand.basenan", std::isnan(m(6, 0)) && std::isnan(m(9, 0)));
    // Gini of (3/4, 1/5): |3/4 - 1/5| * 2 / (2 * 4 * mean) = 11 / 38
    tk_near(tt, "fairdf.hand.gini", sum[2], 11.0 / 38.0, 1e-15);
    tk_near(tt, "fairdf.hand.bas", sum[3], -0.55 * 11.0 / 38.0, 1e-15);
    tk_near(tt, "fairdf.hand.dmin", sum[0], 0.2 / 0.75, 1e-16);
    tk_true(tt, "fairdf.hand.noboot", std::isnan(sum[4]) && std::isnan(sum[7]));
    // without outcomes the error-rate columns are NaN
    clab_fairdf(pr, nullptr, g, 9, 2, 1, Z, 0, 0, cnt, met, sum, wk);
    tk_true(tt, "fairdf.noy", std::isnan(m(11, 0)) && std::isnan(sum[8]) && std::isnan(sum[9]));
    tk_near(tt, "fairdf.noy.base", m(8, 0), 0.75 / 0.2, 1e-15);
    tk_near(tt, "fairdf.noy.cnt", cnt[2 * 2], 0, 0);
    // the bootstrap: the same seed the same ends, ordered
    double s2[11];
    clab_fairdf(pr, y, g, 9, 2, -1, Z, 200, 7, cnt, met, sum, wk);
    clab_fairdf(pr, y, g, 9, 2, -1, Z, 200, 7, cnt, met, s2, wk);
    for (int k = 4; k < 8; ++k) tk_near(tt, "fairdf.boot.same", s2[k], sum[k], 0);
    tk_true(tt, "fairdf.boot.order", sum[4] <= sum[5] && sum[6] <= sum[7]);
    tk_true(tt, "fairdf.boot.range", sum[4] >= 0 && sum[5] <= 1);
    // guards
    const int bad[2] = {0, 2}, g2[2] = {0, 0};
    tk_true(tt, "fairdf.guard.g", clab_fairdf(pr, y, g, 9, 1, -1, Z, 0, 0, cnt, met, sum, wk) == CLAB_EARG);
    tk_true(tt, "fairdf.guard.pred", clab_fairdf(bad, nullptr, g, 2, 2, -1, Z, 0, 0, cnt, met, sum, wk) == CLAB_EDOM);
    tk_true(tt, "fairdf.guard.y", clab_fairdf(pr, bad, g, 2, 2, -1, Z, 0, 0, cnt, met, sum, wk) == CLAB_EDOM);
    tk_true(tt, "fairdf.guard.empty", clab_fairdf(pr, nullptr, g2, 2, 2, -1, Z, 0, 0, cnt, met, sum, wk) == CLAB_EDOM);
    tk_true(tt, "fairdf.guard.code", clab_fairdf(pr, nullptr, bad, 2, 2, -1, Z, 0, 0, cnt, met, sum, wk) == CLAB_EARG);
    tk_true(tt, "fairdf.guard.z", clab_fairdf(pr, nullptr, g, 9, 2, -1, 0, 0, 0, cnt, met, sum, wk) == CLAB_EDOM);
  }
  // 2. tmpaud by hand: one city, two cells. Cell 0: group 0 flagged, group
  // 1 not; cell 1: both flagged. The pooled shares are 1 and 1/2, so the
  // base is group 0.
  {
    const int pr[4] = {1, 0, 1, 1}, g[4] = {0, 1, 0, 1}, c[4] = {0, 0, 1, 1}, cy[2] = {0, 0};
    double co[12], yo[9], sum[3], wk[12];
    tk_true(tt, "tmpaud.rc", clab_tmpaud(pr, g, c, 4, 2, 2, cy, 1, -1, co, yo, sum, wk) == CLAB_OK);
    const double c0[6] = {2, 1, 0, -1, 0.5, -0.5}, c1[6] = {2, 1, 1, 0, 0, 0};
    for (int k = 0; k < 6; ++k) {
      tk_near(tt, "tmpaud.cell0", co[k], c0[k], 1e-15);
      tk_near(tt, "tmpaud.cell1", co[6 + k], c1[k], 1e-15);
    }
    const double y0[9] = {2, 0.5, -0.5, 0.25, -0.25, 0, 1, 1, 0};
    for (int k = 0; k < 9; ++k) tk_near(tt, "tmpaud.city", yo[k], y0[k], 1e-15);
    tk_near(tt, "tmpaud.rng", sum[0], 1, 0);
    tk_true(tt, "tmpaud.spread", std::isnan(sum[1]));
    tk_near(tt, "tmpaud.base", sum[2], 0, 0);
    // base group 1 is flagged nowhere in cell 0: no finite ratio there
    clab_tmpaud(pr, g, c, 4, 2, 2, cy, 1, 1, co, yo, sum, wk);
    tk_true(tt, "tmpaud.inf", std::isnan(co[2]) && co[8] == 1);
    // a cell with one group is skipped
    const int g1[4] = {0, 0, 0, 1};
    clab_tmpaud(pr, g1, c, 4, 2, 2, cy, 1, 0, co, yo, sum, wk);
    tk_true(tt, "tmpaud.skip", co[1] == 0 && std::isnan(co[2]) && co[7] == 1 && yo[0] == 1);
    const int cb[2] = {0, 3};
    tk_true(tt, "tmpaud.guard", clab_tmpaud(pr, g, c, 4, 2, 2, cb, 1, 0, co, yo, sum, wk) == CLAB_EARG);
  }
  // 3. lipaud by hand: scores 0, 1, 3; pairs (0, 1), (0, 2), (1, 2) at
  // distances 1, 0, 4: gaps 1, 3, 2, ratios 1, Inf, 1/2
  {
    const double h[3] = {0, 1, 3}, d[3] = {1, 0, 4};
    const int pi[3] = {0, 0, 1}, pj[3] = {1, 2, 2};
    double gap[3], rat[3], sum[6];
    int v[3];
    tk_true(tt, "lipaud.rc", clab_lipaud(h, 3, pi, pj, d, 3, 0.75, gap, rat, v, sum) == CLAB_OK);
    tk_near(tt, "lipaud.gap", gap[1], 3, 0);
    tk_true(tt, "lipaud.inf", std::isinf(rat[1]) && rat[2] == 0.5);
    tk_true(tt, "lipaud.viol", v[0] == 1 && v[1] == 1 && v[2] == 0);
    tk_true(tt, "lipaud.sum", std::isinf(sum[0]) && sum[1] == 2 && sum[3] == 3 && sum[4] == 1 && sum[5] == 1);
    // all distances 1, and a large enough constant: no violation
    clab_lipaud(h, 3, pi, pj, nullptr, 3, 3, gap, rat, v, sum);
    tk_true(tt, "lipaud.ok", sum[0] == 3 && sum[1] == 0 && sum[2] == 0);
    const int bad[1] = {3};
    tk_true(tt, "lipaud.guard", clab_lipaud(h, 3, bad, pj, d, 1, 1, gap, rat, v, sum) == CLAB_EARG);
    const double dn[1] = {-1};
    tk_true(tt, "lipaud.guard.d", clab_lipaud(h, 3, pi, pj, dn, 1, 1, gap, rat, v, sum) == CLAB_EDOM);
  }
  // 4. nsbnd: the corners by hand, and every admissible noise on a grid
  // stays between them
  {
    const double q[2] = {0.35, 0.55};
    double lo[2], hi[2], bd[1];
    int ord[1];
    tk_true(tt, "nsbnd.rc", clab_nsbnd(q, 2, 0, 0.1, 0, 0.2, lo, hi, ord, bd) == CLAB_OK);
    tk_near(tt, "nsbnd.lo", lo[0], 0.25 / 0.9, 1e-15);
    tk_near(tt, "nsbnd.hi", hi[1], 0.55 / 0.8, 1e-15);
    tk_true(tt, "nsbnd.ord", ord[0] == 1);
    tk_near(tt, "nsbnd.bd", bd[0], 1 - 0.35 / 0.5, 1e-15);
    // the last rate is above 1 - blo: no admissible noise records it
    for (int t = 0; t < 8; ++t) {
      const double r = t == 7 ? 0.98 : 0.02 + 0.15 * t;
      double l1, h1;
      clab_nsbnd(&r, 1, 0.01, 0.12, 0.03, 0.25, &l1, &h1, ord, bd);
      double mn = 2, mx = -1;
      for (int i = 0; i <= 111; ++i) {
        for (int j = 0; j <= 221; ++j) {
          // grid steps of 0.001; the points a = pobs and b = 1 - pobs added
          const double a = i == 111 ? r : 0.01 + 0.001 * i, b = j == 221 ? 1 - r : 0.03 + 0.001 * j;
          if (a > 0.12 + 1e-12 || b > 0.25 + 1e-12 || a < 0.01 || b < 0.03) continue;
          const double p = (r - a) / (1 - a - b);
          if (p < -1e-12 || p > 1 + 1e-12) continue;
          mn = std::fmin(mn, p);
          mx = std::fmax(mx, p);
        }
      }
      if (mx < 0) {
        tk_true(tt, "nsbnd.grid.none", std::isnan(l1) && std::isnan(h1));
      } else {
        tk_near(tt, "nsbnd.grid.lo", l1, mn, 1e-12);
        tk_near(tt, "nsbnd.grid.hi", h1, mx, 1e-12);
      }
    }
    // known noise: the true rate
    const double r = 0.4 * 0.8 + 0.6 * 0.1;
    clab_nsbnd(&r, 1, 0.1, 0.1, 0.2, 0.2, lo, hi, ord, bd);
    tk_near(tt, "nsbnd.known", lo[0], 0.4, 1e-15);
    tk_near(tt, "nsbnd.known.hi", hi[0], 0.4, 1e-15);
    tk_true(tt, "nsbnd.guard", clab_nsbnd(q, 2, 0, 0.6, 0, 0.4, lo, hi, ord, bd) == CLAB_EDOM);
    const double q2[2] = {0.3, 0.3};
    clab_nsbnd(q2, 2, 0, 0, 0, 0, lo, hi, ord, bd);
    tk_true(tt, "nsbnd.equal", ord[0] == 0 && std::isnan(bd[0]));
  }
  // 5. impfpr: Chouldechova's identity on random confusion tables
  {
    xrng r{77};
    for (int k = 0; k < 200; ++k) {
      const double tp = 1 + 99 * r.unif(), fp = 1 + 99 * r.unif(), fn = 1 + 99 * r.unif(),
                   tn = 1 + 99 * r.unif();
      const double p = (tp + fn) / (tp + fp + fn + tn), ppv = tp / (tp + fp), fnr = fn / (tp + fn);
      double f;
      clab_impfpr(&p, &ppv, &fnr, 1, &f);
      tk_near(tt, "impfpr.table", f, fp / (fp + tn), 1e-13);
    }
    const double p = 1, a = 0.5;
    double f;
    tk_true(tt, "impfpr.guard", clab_impfpr(&p, &a, &a, 1, &f) == CLAB_EDOM);
  }
  // 6. rnkres by hand: 0.2, 0.35, 0.8 with half-widths 0.1, 0.1, 0.05
  {
    const double e[3] = {0.2, 0.35, 0.8}, w[3] = {0.1, 0.1, 0.05};
    int id[9];
    double o[4];
    tk_true(tt, "rnkres.rc", clab_rnkres(e, w, 3, id, o) == CLAB_OK);
    const int want[9] = {-1, 0, 1, 0, -1, 1, 1, 1, -1};
    for (int k = 0; k < 9; ++k) tk_ieq(tt, "rnkres.id", id[k], want[k]);
    tk_near(tt, "rnkres.share", o[2], 1.0 / 3, 1e-16);
    tk_near(tt, "rnkres.res", o[3], 0.2, 0);
    const double wn[3] = {0.1, -1, 0};
    tk_true(tt, "rnkres.guard", clab_rnkres(e, wn, 3, id, o) == CLAB_EDOM);
  }
}

}  // namespace

void t_fair(tk &t, int) {
  run(t);
  run2(t);
}

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
