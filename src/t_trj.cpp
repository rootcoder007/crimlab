// SPDX-License-Identifier: AGPL-3.0-or-later
//
// Standalone tests for clab_zgll and clab_zgem (src/trj.cpp). Part of the
// self-test (selftest.cpp).

#include <cmath>
#include <cstdint>

#include "crimlab.h"
#include "tchk.h"

namespace {

// Reference values, dev/doss/crimCV.R (crimCV 1.0.0, R 4.5.0): Y, 100 units
// x 6 periods of simulated two-group zero-inflated counts (row-major; -1 is
// missing), Xb = bs(oyr, degree = 2, intercept = TRUE), Zb = bs(oyr, 1,
// intercept = TRUE); cbeta, cgamma, cprob, cllike, cgwt (first 5 rows):
// crimCV(Y, 2, dpolyp = 2, dpolyl = 1, model = "ZIP", init = 10) after
// set.seed(3); pbeta, pll: glm(y ~ X - 1, family = poisson) over the
// observed cells.
const int Y[600] = {0, 0, 4, 0, 0, 5, 0, 1, 1, 4, 3, 4, 0, 1, 2, -1, 0, 6, 1, 0, 1, 0, 0, 0, 0, 1, 2, 5, 2, 6, 0, 1, 0, 3, 5, 0, 1, 0, 0, 3, 4, 4, 1, 2, 0, 3, 6, 0, 1, 0, 1, 0, 5, 7, 2, 0, 5, 0, 4, 4, 0, 3, 3, 0, 8, 0, 0, 6, 0, 3, 4, 6, 0, 2, 1, 4, 5, 6, 0, 0, 1, 4, 2, 2, 2, 0, 1, 4, 0, 4, 0, 0, 0, 5, 3, 3, 0, 2, 0, 2, 4, 3, 2, 4, 3, 0, 0, 12, 3, 2, 3, 4, 0, 5, 1, 2, 2, 0, 4, 3, 2, 3, 2, 1, 7, 4, 0, 0, 2, 3, 2, 0, 1, 0, 0, 0, 0, 10, 2, 2, 0, 3, 0, 7, 4, 2, 2, 2, 5, 0, 0, 1, 1, 3, 6, 1, 1, 3, 2, 2, 0, 0, 1, 0, 0, 5, 0, 0, 1, 0, 4, 4, 5, 5, 0, 0, 3, 4, 0, 5, 0, 0, 6, 0, 3, 0, 1, 1, 2, 3, 3, 0, 5, 3, 5, 4, 2, 0, 3, 1, 0, 0, 6, 4, 1, 0, 1, 0, 2, 0, 1, 2, 4, 6, 0, 12, 0, 1, 2, 3, 0, 6, 1, 0, 2, 1, 0, 0, 0, 1, 0, 2, 4, 3, 0, 2, 1, 2, 6, 4, 5, 0, 5, 3, 0, 1, 2, 6, 2, 1, 0, 1, 2, 7, 3, 0, 1, 3, 0, 3, 2, 0, 0, 0, 4, 1, 2, 1, 0, 0, 5, 0, 1, 2, 1, 1, 2, 2, 0, 2, 0, 0, 6, 3, 0, 5, 4, 2, 3, 5, 0, 3, 2, 1, 2, 4, 1, 3, 2, 1, 0, 4, 6, 0, 2, 2, 4, 0, 0, 4, 1, 3, 3, 4, 2, 1, 2, 3, 2, 3, 1, 2, 1, 0, 6, 1, 2, 4, 1, 0, 4, 2, 2, 0, 0, 2, 0, 0, 0, 1, 2, 2, 2, 0, 0, 2, 1, 1, 3, 3, 4, 1, 1, 0, 4, 1, 1, 2, 2, 2, 5, 2, 0, 0, 1, 3, 1, 3, 3, 2, 3, 0, 3, 2, 0, 1, 0, 0, 3, 4, 4, 2, 2, 2, 5, 2, 2, 2, 3, 2, 7, 0, 3, 0, 2, 0, 3, 1, 0, 1, 0, 0, 4, 0, 0, 5, 0, 0, 0, 4, 5, 7, 3, 0, 2, 1, 3, 1, 2, 1, 0, 5, 3, 3, 0, 1, 4, 4, 2, 1, 2, 2, 4, 3, 0, 2, 1, 0, 5, 4, 3, 3, 1, 1, 6, 1, 1, 3, 1, 0, 0, 1, 4, 0, 1, 1, 1, 0, 4, 2, 0, 1, 2, 3, 0, 0, 0, 1, 1, 4, 3, 0, 2, 2, 7, 2, 3, 1, 2, 1, 3, 0, 5, 1, 3, 0, 1, 2, 3, 2, 1, 1, 3, 2, 0, 0, 0, 0, 1, 4, 3, 5, 3, 3, 4, 4, 2, 3, 0, 1, 3, 2, 3, 2, 2, 2, 4, 1, 0, 3, 1, 1, 8, 0, 0, 1, 5, 0, 2, 4, 2, 0, 1, 0, 5, 1, 1, 1, 2, 0, 4, 1, 2, 0, 1, 2, 7, 0, 1, 3, 1, 2, 3, 4, 0, 0, 6, 2, 4, 2, 1, 3, 2, 1, 7, 3, 2, 1, 0, 0, 1, 0, 4, 5, 1, 1, 1, 0, 1, 4, 0, 2, 0, 2, 3, 3, 2, 0, 7, 2, 2, 1, 3, 4, 2, 2, 3, 0, 1, 1};
const double Xb[18] = {1, 0.64000000000000012, 0.35999999999999999, 0.16000000000000003, 0.03999999999999998, 0, 0, 0.32000000000000006, 0.47999999999999998, 0.47999999999999998, 0.31999999999999995, 0, 0, 0.040000000000000008, 0.16000000000000003, 0.35999999999999999, 0.64000000000000012, 1};
const double Zb[12] = {1, 0.80000000000000004, 0.59999999999999998, 0.40000000000000002, 0.19999999999999996, 0, 0, 0.20000000000000001, 0.40000000000000002, 0.59999999999999998, 0.80000000000000004, 1};
const double cbeta[6] = {-0.80721687685316923, 1.1370642596975833, 1.5615165658744286, 1.0884266637803843, 0.67648259053554904, 0.88937775291607102};
const double cgamma[4] = {-21.062130434380762, -20.52827810438664, -1.7007834683805703, -0.86774487871818629};
const double cprob[2] = {0.13423431608948957, 0.86576568391051045};
const double cllike = -1120.0694216133402;
const double cgwt[10] = {0.030190941184994758, 0.90413223315639546, 0.42865048228242342, 0.00036932446118156698, 0.94747394716528555, 0.96980905881500534, 0.095867766843604499, 0.57134951771757647, 0.99963067553881846, 0.052526052834714364};
const double pbeta[3] = {0.78434452946132105, 0.48451707169193858, 0.74883243333068827};
const double pll = -1203.1411172776643;

struct xrng {
  std::uint64_t s;
  double unif() {
    s ^= s >> 12;
    s ^= s << 25;
    s ^= s >> 27;
    return static_cast<double>((s * 2685821657736338717ULL) >> 11) * 0x1p-53;
  }
};

double yc[600];  // column-major 100 x 6
void load() {
  for (int i = 0; i < 100; ++i)
    for (int t = 0; t < 6; ++t) yc[i + 100 * t] = Y[6 * i + t];
}

// the mixture log-likelihood written out directly (no log-sum-exp, no
// compensated sums): sum_i log sum_k prob_k prod_t f(y_it)
double brute(const double *y, int n, int nt, const double *x, int p, const double *z, int q,
             int ng, const double *par, const double *prob) {
  double ll = 0;
  for (int i = 0; i < n; ++i) {
    double li = 0;
    for (int k = 0; k < ng; ++k) {
      const double *b = par + k * (p + q), *g = b + p;
      double pr = prob[k];
      for (int t = 0; t < nt; ++t) {
        const double v = y[i + n * t];
        if (!(v >= 0)) continue;
        double eta = 0, a = 0;
        for (int j = 0; j < p; ++j) eta += x[t + nt * j] * b[j];
        for (int j = 0; j < q; ++j) a += z[t + nt * j] * g[j];
        const double lam = std::exp(eta), rho = q > 0 ? std::exp(a) / (1 + std::exp(a)) : 0;
        const double pois = std::exp(v * eta - lam - std::lgamma(v + 1));
        pr *= (v == 0 ? rho : 0) + (1 - rho) * pois;
      }
      li += pr;
    }
    ll += std::log(li);
  }
  return ll;
}

void run(tk &tt) {
  load();
  static double w[100 * 3], work[256];
  double par[10], prob[2], ll;
  int it;
  for (int k = 0; k < 2; ++k) {
    for (int j = 0; j < 3; ++j) par[5 * k + j] = cbeta[3 * k + j];
    for (int j = 0; j < 2; ++j) par[5 * k + 3 + j] = cgamma[2 * k + j];
    prob[k] = cprob[k];
  }
  // 1. parity with crimCV at its own estimates: log-likelihood and the
  // posterior memberships
  tk_true(tt, "trj.ll.rc", clab_zgll(yc, 100, 6, nullptr, Xb, 3, Zb, 2, 2, par, prob, w, &ll, work) == CLAB_OK);
  tk_near(tt, "trj.ll.crimcv", ll / cllike, 1, 1e-12);
  for (int i = 0; i < 5; ++i)
    for (int k = 0; k < 2; ++k) tk_near(tt, "trj.gwt", w[i + 100 * k], cgwt[i + 5 * k], 1e-9);
  // 2. the direct product formula on the same parameters
  tk_near(tt, "trj.ll.brute", ll / brute(yc, 100, 6, Xb, 3, Zb, 2, 2, par, prob), 1, 1e-12);
  // 3. the reference estimate is a stationary point (central differences
  // of the log-likelihood vanish) on a flat ridge (a zero-inflation
  // intercept of -21); EM started from it reaches a higher maximum (shares
  // 0.42 / 0.58), confirmed by the direct formula.
  {
    double gmax = 0;
    for (int j = 0; j < 10; ++j) {
      double pp[10], lp, lm;
      for (int l = 0; l < 10; ++l) pp[l] = par[l];
      pp[j] = par[j] + 1e-6;
      clab_zgll(yc, 100, 6, nullptr, Xb, 3, Zb, 2, 2, pp, prob, w, &lp, work);
      pp[j] = par[j] - 1e-6;
      clab_zgll(yc, 100, 6, nullptr, Xb, 3, Zb, 2, 2, pp, prob, w, &lm, work);
      gmax = std::fabs(lp - lm) / 2e-6 > gmax ? std::fabs(lp - lm) / 2e-6 : gmax;
    }
    tk_true(tt, "trj.crimcv.stationary", gmax < 1e-4);
    double p2[10], pr2[2];
    for (int j = 0; j < 10; ++j) p2[j] = par[j];
    pr2[0] = prob[0];
    pr2[1] = prob[1];
    tk_true(tt, "trj.em.rc", clab_zgem(yc, 100, 6, nullptr, Xb, 3, Zb, 2, 2, p2, pr2, w, &ll, 2000, 1e-14,
                               &it, work) == CLAB_OK);
    tk_true(tt, "trj.em.higher", ll > cllike + 30);
    tk_near(tt, "trj.em.brute", ll / brute(yc, 100, 6, Xb, 3, Zb, 2, 2, p2, pr2), 1, 1e-12);
    tk_true(tt, "trj.em.shares", std::fabs(pr2[0] - 0.42) < 0.01 || std::fabs(pr2[1] - 0.42) < 0.01);
  }
  // 4. one Poisson group (q = 0) is a Poisson regression: parity with
  // glm(poisson) in base R
  {
    double b[3] = {0, 0, 0}, pr1[1] = {1};
    clab_zgem(yc, 100, 6, nullptr, Xb, 3, nullptr, 0, 1, b, pr1, w, &ll, 50, 1e-14, &it, work);
    for (int j = 0; j < 3; ++j) tk_near(tt, "trj.glm.beta", b[j], pbeta[j], 1e-8);
    tk_near(tt, "trj.glm.ll", ll / pll, 1, 1e-12);
  }
  // 5. EM is monotone, and at its fixed point the observed log-likelihood
  // is stationary (central differences of clab_zgll, step 1e-5)
  {
    double p2[10] = {0.5, 0.5, 0.5, -1, -1, 1, 0.5, 0, -2, -2}, pr2[2] = {0.5, 0.5};
    double prev = -HUGE_VAL;
    bool mono = true;
    for (int r = 0; r < 300; ++r) {
      clab_zgem(yc, 100, 6, nullptr, Xb, 3, Zb, 2, 2, p2, pr2, w, &ll, 1, 1e-300, &it, work);
      mono = mono && ll >= prev - 1e-10 * std::fabs(ll);
      prev = ll;
    }
    tk_true(tt, "trj.monotone", mono);
    clab_zgem(yc, 100, 6, nullptr, Xb, 3, Zb, 2, 2, p2, pr2, w, &ll, 5000, 1e-15, &it, work);
    double gmax = 0;
    for (int j = 0; j < 10; ++j) {
      double pp[10], lp, lm;
      for (int l = 0; l < 10; ++l) pp[l] = p2[l];
      pp[j] = p2[j] + 1e-5;
      clab_zgll(yc, 100, 6, nullptr, Xb, 3, Zb, 2, 2, pp, pr2, w, &lp, work);
      pp[j] = p2[j] - 1e-5;
      clab_zgll(yc, 100, 6, nullptr, Xb, 3, Zb, 2, 2, pp, pr2, w, &lm, work);
      const double gr = (lp - lm) / 2e-5;
      gmax = std::fabs(gr) > gmax ? std::fabs(gr) : gmax;
    }
    tk_true(tt, "trj.stationary", gmax < 1e-3);
  }
  // 6. NaN is missing, exactly like a negative count, and 2 units fit.
  {
    double y2[600];
    for (int j = 0; j < 600; ++j) y2[j] = yc[j];
    y2[4 + 100 * 1] = std::nan("");
    double l1, l2;
    clab_zgll(y2, 100, 6, nullptr, Xb, 3, Zb, 2, 2, par, prob, w, &l1, work);
    y2[4 + 100 * 1] = -1;
    clab_zgll(y2, 100, 6, nullptr, Xb, 3, Zb, 2, 2, par, prob, w, &l2, work);
    tk_near(tt, "trj.nan", l1, l2, 0);
    double ys[12], p2[10], pr2[2] = {0.5, 0.5};
    for (int t = 0; t < 6; ++t) {
      ys[0 + 2 * t] = yc[0 + 100 * t] < 0 ? 0 : yc[0 + 100 * t];
      ys[1 + 2 * t] = yc[1 + 100 * t] < 0 ? 0 : yc[1 + 100 * t];
    }
    for (int j = 0; j < 10; ++j) p2[j] = 0.1 * j - 0.5;
    tk_true(tt, "trj.n2", clab_zgem(ys, 2, 6, nullptr, Xb, 3, Zb, 2, 2, p2, pr2, w, &ll, 100, 1e-10, &it,
                            work) == CLAB_OK && std::isfinite(ll));
  }
  // 7. offsets: log exposure 2 everywhere is the same as adding log 2 to
  // the intercept-spanning basis (the B-spline rows sum to one)
  {
    static double off[600];
    for (int j = 0; j < 600; ++j) off[j] = std::log(2.0);
    double p2[10], l1, l2;
    for (int j = 0; j < 10; ++j) p2[j] = par[j];
    clab_zgll(yc, 100, 6, off, Xb, 3, Zb, 2, 2, p2, prob, w, &l1, work);
    for (int k = 0; k < 2; ++k)
      for (int j = 0; j < 3; ++j) p2[5 * k + j] += std::log(2.0);
    clab_zgll(yc, 100, 6, nullptr, Xb, 3, Zb, 2, 2, p2, prob, w, &l2, work);
    tk_near(tt, "trj.offset", l1 / l2, 1, 1e-12);
  }
  // 8. guards
  {
    double y2[600];
    for (int j = 0; j < 600; ++j) y2[j] = yc[j];
    y2[7] = 1.5;
    tk_true(tt, "trj.guard.int", clab_zgll(y2, 100, 6, nullptr, Xb, 3, Zb, 2, 2, par, prob, w, &ll, work) == CLAB_EDOM);
    const double bad[2] = {0.5, 0.6};
    tk_true(tt, "trj.guard.prob", clab_zgll(yc, 100, 6, nullptr, Xb, 3, Zb, 2, 2, par, bad, w, &ll, work) == CLAB_EDOM);
    tk_true(tt, "trj.guard.p", clab_zgll(yc, 100, 6, nullptr, Xb, 0, Zb, 2, 2, par, prob, w, &ll, work) == CLAB_EARG);
  }
}


}  // namespace

void t_trj(tk &t, int) { run(t); }

#ifdef CLAB_MAIN
#include <chrono>
#include <cstdio>

// timing: 1,000 units x 10 periods, 2 groups, cubic bases for the counts
// and the zero inflation (crimCV ZIP, dpolyp = dpolyl = 3, init = 20, i.e.
// 40 starts, took 32.544 s, one run, zeus; its default ZIPt 22.262 s).
// Here 40 EM runs from random starts on a polynomial basis of the same
// size.
void b_trj() {
  const int n = 1000, nt = 10, ng = 2, p = 4, q = 4, m = p + q;
  static double y[n * nt], w[n * ng], work[256];
  double x[nt * p], z[nt * q];
  xrng r{1};
  for (int t = 0; t < nt; ++t) {
    const double o = t / 9.0;
    for (int j = 0; j < p; ++j) x[t + nt * j] = z[t + nt * j] = std::pow(o, j);
  }
  for (int i = 0; i < n; ++i) {
    const int g = r.unif() < 0.4;
    for (int t = 0; t < nt; ++t) {
      const double o = t / 9.0, lam = g ? std::exp(0.2 + 1.5 * o) : std::exp(1.2 - o);
      // Poisson by inversion
      double u = r.unif(), pk = std::exp(-lam), c = pk;
      int k = 0;
      while (u > c && k < 200) {
        ++k;
        pk *= lam / k;
        c += pk;
      }
      y[i + n * t] = r.unif() < 0.15 ? 0 : k;
    }
  }
  double t[5], best = -HUGE_VAL;
  for (int rep = 0; rep < 5; ++rep) {
    const auto a = std::chrono::steady_clock::now();
    for (int s = 0; s < 40; ++s) {
      double par[ng * m], prob[ng] = {0.5, 0.5}, ll;
      int it;
      for (int j = 0; j < ng * m; ++j) par[j] = (r.unif() - 0.5) * (j % m == 0 ? 2.0 : 1.0) - (j % m >= p ? 1.5 : 0.0);
      clab_zgem(y, n, nt, nullptr, x, p, z, q, ng, par, prob, w, &ll, 500, 1e-9, &it, work);
      best = ll > best ? ll : best;
    }
    t[rep] = std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count();
  }
  std::printf("BENCH clab_zgem n=1000x10 ng=2 p=q=4, 40 starts: median %.3f s (k=5), best ll %.4f\n",
              tk_med5(t), best);
}
#endif
