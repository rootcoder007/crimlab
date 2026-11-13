// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of the random forest (src/rf.cpp): error against a known
// regression function and a known probability, determinism under a seed,
// exact invariance to increasing transforms of the columns, the
// out-of-bag predictions, the edge cases of the splitting and every
// refusal. Built with the self-test.
#include <cmath>
#include <cstring>

#include "crimlab.h"
#include "learn.h"
#include "tchk.h"

namespace {

const int N = 600, NT = 300, P = 5;
double X[N * P], XT[NT * P], Y[N], YB[N], FT[NT], PT[NT], PR[NT], PR2[NT],
    OOB[N];
double W[P * 255 + N + 512 + NT + 2 * N + 3 * (2 * N + 1)];
int IW[(N + NT) * P / 4 + 1 + 3 * N + 1 + P + N + NT];

double rnorm(rng &g) {
  return std::sqrt(-2.0 * std::log(runi(g))) * std::cos(6.283185307179586 * runi(g));
}
double fx(const double *x, int ld, int i) {
  return 2.0 * std::sin(x[i]) + (x[i + ld] > 0 ? 1.0 : -1.0) + 0.5 * x[i + 2 * ld];
}

int rf(const double *x, const double *y, const double *xt, int ntree, int mtry,
       int minnode, double frac, int repl, unsigned int seed, double *pred,
       double *oob, int n = N, int p = P, int nt = NT) {
  return clab_rfor(x, n, p, y, xt, nt, ntree, mtry, minnode, frac, repl, seed,
                   pred, oob, W, IW);
}

}  // namespace

void t_rf(tk &t, int) {
  rng g = {23};
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < P; ++j) X[i + j * N] = rnorm(g);
    Y[i] = fx(X, N, i) + 0.5 * rnorm(g);
    YB[i] = runi(g) < 1.0 / (1.0 + std::exp(-fx(X, N, i))) ? 1.0 : 0.0;
  }
  for (int i = 0; i < NT; ++i) {
    for (int j = 0; j < P; ++j) XT[i + j * NT] = rnorm(g);
    FT[i] = fx(XT, NT, i);
    PT[i] = 1.0 / (1.0 + std::exp(-FT[i]));
  }
  // regression: mean squared error against the true function well below
  // the variance of the function (about 3.3) and below 0.6
  tk_ieq(t, "rf.rc", rf(X, Y, XT, 200, 2, 5, 1.0, 1, 1, PR, OOB), CLAB_OK);
  double mse = 0.0, vf = 0.0, mf = 0.0;
  for (int i = 0; i < NT; ++i) mf += FT[i] / NT;
  for (int i = 0; i < NT; ++i) {
    mse += (PR[i] - FT[i]) * (PR[i] - FT[i]) / NT;
    vf += (FT[i] - mf) * (FT[i] - mf) / NT;
  }
  tk_lt(t, "rf.mse", mse, 0.6);
  tk_gt(t, "rf.vf", vf, 2.0);
  // out of bag: every row left out by some of 200 trees; error near the
  // noise variance plus the test error
  double om = 0.0;
  for (int i = 0; i < N; ++i) om += (OOB[i] - Y[i]) * (OOB[i] - Y[i]) / N;
  tk_lt(t, "rf.oob", om, 0.25 + 0.8);
  // determinism, and a different seed gives a different forest
  rf(X, Y, XT, 200, 2, 5, 1.0, 1, 1, PR2, nullptr);
  tk_true(t, "rf.same", std::memcmp(PR, PR2, sizeof PR) == 0);
  rf(X, Y, XT, 200, 2, 5, 1.0, 1, 2, PR2, nullptr);
  tk_true(t, "rf.seed", std::memcmp(PR, PR2, sizeof PR) != 0);
  // increasing transforms of every column: the same trees exactly
  {
    static double xe[N * P], xte[NT * P];
    for (int k = 0; k < N * P; ++k) xe[k] = std::exp(X[k]) * 3.0 - 1.0;
    for (int k = 0; k < NT * P; ++k) xte[k] = std::exp(XT[k]) * 3.0 - 1.0;
    rf(xe, Y, xte, 200, 2, 5, 1.0, 1, 1, PR2, nullptr);
    tk_true(t, "rf.mono", std::memcmp(PR, PR2, sizeof PR) == 0);
  }
  // probability: in [0, 1], error against the true probability small
  rf(X, YB, XT, 200, 2, 10, 1.0, 1, 3, PR, nullptr);
  double pm = 0.0, lo = 1.0, hi = 0.0;
  for (int i = 0; i < NT; ++i) {
    pm += (PR[i] - PT[i]) * (PR[i] - PT[i]) / NT;
    lo = std::fmin(lo, PR[i]);
    hi = std::fmax(hi, PR[i]);
  }
  tk_lt(t, "rf.prob", pm, 0.05);
  tk_true(t, "rf.prob01", lo >= 0.0 && hi <= 1.0);
  // subsampling without replacement; a single tree leaves some rows in
  // every bag (no out-of-bag prediction: NaN)
  tk_ieq(t, "rf.sub", rf(X, Y, XT, 50, 5, 3, 0.5, 0, 4, PR, OOB), CLAB_OK);
  rf(X, Y, XT, 1, 2, 5, 1.0, 1, 5, PR, OOB);
  int nn = 0;
  for (int i = 0; i < N; ++i) nn += std::isnan(OOB[i]);
  tk_true(t, "rf.oobnan", nn > 100 && nn < N - 100);
  // few distinct values: splits at midpoints; a constant column and a
  // node too small to split
  {
    static double xs[N * P];
    for (int i = 0; i < N; ++i) {
      xs[i] = i % 3;
      for (int j = 1; j < P; ++j) xs[i + j * N] = 1.0;
    }
    double y3[N], xt1[P] = {0.4, 1, 1, 1, 1}, p1;
    for (int i = 0; i < N; ++i) y3[i] = (i % 3 == 0) ? 5.0 : 1.0;
    rf(xs, y3, xt1, 20, P, 1, 1.0, 1, 6, &p1, nullptr, N, P, 1);
    tk_eq(t, "rf.mid", p1, 5.0, 1e-15);
    xt1[0] = 0.6;
    rf(xs, y3, xt1, 20, P, 1, 1.0, 1, 6, &p1, nullptr, N, P, 1);
    tk_eq(t, "rf.mid2", p1, 1.0, 1e-15);
    for (int i = 0; i < N; ++i) xs[i] = 1.0;
    rf(xs, Y, xt1, 3, 1, 5, 1.0, 0, 7, &p1, nullptr, N, P, 1);
    double my = 0.0;
    for (int i = 0; i < N; ++i) my += Y[i] / N;
    tk_eq(t, "rf.const", p1, my, 1e-12);
    rf(X, Y, xt1, 1, 1, N, 1.0, 0, 7, &p1, nullptr, N, P, 1);
    tk_eq(t, "rf.minnode", p1, my, 1e-12);
  }
  // refusals
  tk_ieq(t, "rf.e.n", rf(X, Y, XT, 10, 2, 5, 1.0, 1, 1, PR, nullptr, 1), CLAB_EARG);
  tk_ieq(t, "rf.e.p", rf(X, Y, XT, 10, 2, 5, 1.0, 1, 1, PR, nullptr, N, 0), CLAB_EARG);
  tk_ieq(t, "rf.e.nt", rf(X, Y, XT, 10, 2, 5, 1.0, 1, 1, PR, nullptr, N, P, -1), CLAB_EARG);
  tk_ieq(t, "rf.e.tree", rf(X, Y, XT, 0, 2, 5, 1.0, 1, 1, PR, nullptr), CLAB_EARG);
  tk_ieq(t, "rf.e.mtry", rf(X, Y, XT, 10, 6, 5, 1.0, 1, 1, PR, nullptr), CLAB_EARG);
  tk_ieq(t, "rf.e.node", rf(X, Y, XT, 10, 2, 0, 1.0, 1, 1, PR, nullptr), CLAB_EARG);
  tk_ieq(t, "rf.e.frac", rf(X, Y, XT, 10, 2, 5, 1.5, 1, 1, PR, nullptr), CLAB_EARG);
  tk_ieq(t, "rf.e.m", rf(X, Y, XT, 10, 2, 5, 1e-6, 1, 1, PR, nullptr), CLAB_EARG);
  {
    static double xn[N * P], yn[N], xtn[NT * P];
    std::memcpy(xn, X, sizeof X);
    std::memcpy(yn, Y, sizeof Y);
    std::memcpy(xtn, XT, sizeof XT);
    xn[9] = NAN;
    yn[4] = NAN;
    xtn[2] = INFINITY;
    tk_ieq(t, "rf.e.x", rf(xn, Y, XT, 10, 2, 5, 1.0, 1, 1, PR, nullptr), CLAB_EDOM);
    tk_ieq(t, "rf.e.xt", rf(X, Y, xtn, 10, 2, 5, 1.0, 1, 1, PR, nullptr), CLAB_EDOM);
    tk_ieq(t, "rf.e.y", rf(X, yn, XT, 10, 2, 5, 1.0, 1, 1, PR, nullptr), CLAB_EDOM);
  }
  // the internal entry point refuses a training set of one row
  int tr[1] = {0}, te[1] = {0};
  const lrn::Forest f = {10, 2, 5, 1.0, 1, 1};
  tk_ieq(t, "rf.e.mt", lrn::forest(X, N, P, Y, tr, 1, X, N, te, 1, f, PR, nullptr, W, IW),
         CLAB_EARG);
}
