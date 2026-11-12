// SPDX-License-Identifier: AGPL-3.0-or-later
//
// Standalone tests for clab_mvtll and clab_mvtem (src/mvt.cpp). Built only
// with -DCLAB_MAIN.
#ifdef CLAB_MAIN

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>

#include "crimlab.h"
#include "tk.h"

namespace {

// Reference values, dev/doss/gbmt.R (gbmt 0.1.4, R 4.5.0): D, 60 units x 6
// times of (id, time, y1, y2), two groups, values multiples of 1/8;
// gbeta1/2 (rows intercept, t; columns y1, y2), gS1/2, gprob (mean
// posteriors), gpost (first 5 rows), gll: gbmt(c("y1", "y2"), "id",
// "time", ng = 2, d = 1, data = d, scaling = 0, pruning = FALSE).
const double D[1440] = {1, 1, 1, 3.375, 1, 2, 1.25, 2.625, 1, 3, 3.5, 2.125, 1, 4, 2.75, 2, 1, 5, 3.375, 2, 1, 6, 4.75, 1.5, 2, 1, 0.875, 2.75, 2, 2, 1.375, 2.5, 2, 3, 2.125, 1.5, 2, 4, 3.125, 2.125, 2, 5, 4, 2.375, 2, 6, 4.875, 2.125, 3, 1, 2, 2.625, 3, 2, 2.875, 3.25, 3, 3, 2.125, 2, 3, 4, 3.25, 2.625, 3, 5, 4.875, 2.5, 3, 6, 4.875, 2.125, 4, 1, 0.875, 2, 4, 2, 1, 2.5, 4, 3, 3.375, 2.125, 4, 4, 3, 3, 4, 5, 3.25, 2.25, 4, 6, 3.5, 1.625, 5, 1, 1.125, 2.75, 5, 2, 2.875, 2.625, 5, 3, 3.125, 2.125, 5, 4, 2.875, 1.875, 5, 5, 4, 2, 5, 6, 4, 1.875, 6, 1, 0.875, 3.125, 6, 2, 1.625, 2.75, 6, 3, 2, 2.25, 6, 4, 1.75, 2.125, 6, 5, 2.75, 1.875, 6, 6, 3.875, 1.75, 7, 1, 1.75, 2.875, 7, 2, 2.25, 2.5, 7, 3, 3.125, 2.5, 7, 4, 3.125, 2, 7, 5, 3.75, 1.25, 7, 6, 4.625, 1.75, 8, 1, 1.875, 2.5, 8, 2, 3.375, 2.375, 8, 3, 2.5, 2.25, 8, 4, 3.5, 1.75, 8, 5, 3.875, 2.375, 8, 6, 4.625, 2, 9, 1, 1.25, 2.625, 9, 2, 1.875, 2.625, 9, 3, 1.25, 2.875, 9, 4, 3.375, 2.25, 9, 5, 3.5, 2.375, 9, 6, 3.875, 1, 10, 1, 1, 3.25, 10, 2, 1.375, 3.25, 10, 3, 4, 2.5, 10, 4, 2.25, 2.5, 10, 5, 2.625, 2.125, 10, 6, 4, 1.75, 11, 1, 2.25, 3.125, 11, 2, 1.25, 2, 11, 3, 3.125, 2.5, 11, 4, 2.75, 2.75, 11, 5, 2.625, 1.375, 11, 6, 4.75, 1, 12, 1, 0.75, 2.75, 12, 2, 2.5, 2.75, 12, 3, 2.625, 2, 12, 4, 3.375, 2.375, 12, 5, 3, 2, 12, 6, 5.25, 2.125, 13, 1, 1.375, 2.625, 13, 2, 2, 2.5, 13, 3, 1.625, 2.25, 13, 4, 3.625, 2.5, 13, 5, 3.25, 2.375, 13, 6, 3.625, 2.25, 14, 1, 1.25, 3.25, 14, 2, 1.75, 2.625, 14, 3, 0.875, 3.25, 14, 4, 3, 2.875, 14, 5, 3.25, 1.75, 14, 6, 4, 1.75, 15, 1, 1.5, 2.625, 15, 2, 1.125, 2.625, 15, 3, 1.75, 1.75, 15, 4, 3, 2.375, 15, 5, 4.5, 1.625, 15, 6, 4.5, 1.5, 16, 1, 2.75, 2.5, 16, 2, 2.5, 3, 16, 3, 1.75, 2.25, 16, 4, 2.875, 2, 16, 5, 4, 1.5, 16, 6, 4.375, 2, 17, 1, 1.25, 3, 17, 2, 1.875, 3.125, 17, 3, 3, 2.375, 17, 4, 2.625, 2.5, 17, 5, 3.375, 1.875, 17, 6, 4.625, 1.5, 18, 1, 1.625, 3, 18, 2, 2.375, 3, 18, 3, 2.25, 2.75, 18, 4, 3.75, 2.5, 18, 5, 3.875, 1.75, 18, 6, 4.5, 1.625, 19, 1, 1.375, 3.375, 19, 2, 2.25, 2.875, 19, 3, 2.625, 2.375, 19, 4, 4.25, 2, 19, 5, 2.5, 2, 19, 6, 2.875, 1.625, 20, 1, 1, 1.375, 20, 2, 1.75, 3, 20, 3, 2.125, 2.125, 20, 4, 2.625, 3, 20, 5, 4.375, 1.75, 20, 6, 3.875, 2.25, 21, 1, 1.625, 2.75, 21, 2, 2, 3.25, 21, 3, 1.625, 2.125, 21, 4, 3.25, 2.125, 21, 5, 2.75, 2.125, 21, 6, 3.875, 1.75, 22, 1, 1.625, 2.625, 22, 2, 0.875, 2, 22, 3, 2.125, 2, 22, 4, 3.5, 2.25, 22, 5, 2.625, 2, 22, 6, 4.625, 1.5, 23, 1, 1.5, 3.625, 23, 2, 2.5, 2.375, 23, 3, 1.5, 2.125, 23, 4, 2.75, 1.75, 23, 5, 3, 2.625, 23, 6, 5.125, 1.75, 24, 1, 0.75, 2.875, 24, 2, 2.125, 2.5, 24, 3, 2.625, 2.125, 24, 4, 2.375, 2.625, 24, 5, 3.25, 1.625, 24, 6, 4.875, 1.75, 25, 1, 2.625, 2.5, 25, 2, 1.875, 2.875, 25, 3, 2.875, 2.375, 25, 4, 3, 1.875, 25, 5, 4.125, 2.375, 25, 6, 4.125, 1.875, 26, 1, 4.375, 1.125, 26, 2, 4.625, 1.25, 26, 3, 2.125, 1, 26, 4, 3.5, 1.25, 26, 5, 2.375, 1.75, 26, 6, 1.625, 2.125, 27, 1, 4.125, 0.625, 27, 2, 3.5, 1.875, 27, 3, 3.875, 1.125, 27, 4, 2, 1.5, 27, 5, 2.25, 2, 27, 6, 2.375, 1.875, 28, 1, 3.125, 0.625, 28, 2, 3.75, 1.125, 28, 3, 2.25, 1, 28, 4, 2.125, 1.25, 28, 5, 2.375, 1.125, 28, 6, 2.5, 1.25, 29, 1, 2.875, 1.125, 29, 2, 3.625, 1.375, 29, 3, 1.375, 0.375, 29, 4, 1.75, 1.75, 29, 5, 1.875, 1.125, 29, 6, 2, 1.125, 30, 1, 3.75, 1.625, 30, 2, 2.125, 1.25, 30, 3, 3.625, 1, 30, 4, 2.75, 1.625, 30, 5, 3.625, 0.75, 30, 6, 2.375, 1.75, 31, 1, 3.25, 1.875, 31, 2, 3.625, 1.625, 31, 3, 2.75, 1.75, 31, 4, 2.375, 1.375, 31, 5, 2.75, 1.375, 31, 6, 2.25, 2, 32, 1, 2.625, 1.25, 32, 2, 3.125, 1.625, 32, 3, 2.75, 1.75, 32, 4, 3.125, 1.25, 32, 5, 3.75, 2, 32, 6, 2.75, 2, 33, 1, 4.75, 1.625, 33, 2, 3.875, 0.875, 33, 3, 2.625, 1.25, 33, 4, 2.5, 1.5, 33, 5, 2, 1.75, 33, 6, 2.375, 0.5, 34, 1, 3.375, 1.25, 34, 2, 2.5, 1, 34, 3, 2.875, 1.5, 34, 4, 1.875, 1, 34, 5, 2.75, 1.25, 34, 6, 2, 2.375, 35, 1, 3.375, 0.25, 35, 2, 4, 0.875, 35, 3, 2.5, 1.625, 35, 4, 2.75, 2.25, 35, 5, 2.25, 1.625, 35, 6, 2.25, 0.875, 36, 1, 3.625, 1.875, 36, 2, 3.625, 0.75, 36, 3, 3.25, 1.25, 36, 4, 3.25, 0.875, 36, 5, 2.875, 1.75, 36, 6, 1.625, 1, 37, 1, 3.875, 1.25, 37, 2, 4.375, 1.375, 37, 3, 3.875, 2, 37, 4, 3.375, 1.625, 37, 5, 2.25, 0.875, 37, 6, 2.25, 2, 38, 1, 4.125, 1.25, 38, 2, 4.25, 1.125, 38, 3, 4.25, 1.25, 38, 4, 3.125, 1.125, 38, 5, 2.5, 1.5, 38, 6, 2.125, 1.625, 39, 1, 2.875, 1.375, 39, 2, 2.75, 1.375, 39, 3, 2.875, 0.75, 39, 4, 3.125, 0.875, 39, 5, 3.625, 1.75, 39, 6, 2.75, 1, 40, 1, 4, 0.75, 40, 2, 4.25, 1.125, 40, 3, 1.25, 1.25, 40, 4, 2.625, 1.875, 40, 5, 2.5, 1.625, 40, 6, 2.5, 1.875, 41, 1, 3.75, 1, 41, 2, 2.875, 1.25, 41, 3, 1.875, 1.125, 41, 4, 2.375, 1, 41, 5, 3, 1.625, 41, 6, 1.75, 1.375, 42, 1, 3.25, 0.75, 42, 2, 4.375, 1.25, 42, 3, 2.5, 1.875, 42, 4, 2.75, 0.5, 42, 5, 2.5, 1.625, 42, 6, 3, 1.625, 43, 1, 2.875, 1, 43, 2, 4.5, 1.625, 43, 3, 2.75, 1.625, 43, 4, 3.25, 0.875, 43, 5, 3.125, 1, 43, 6, 2.625, 1.75, 44, 1, 4.125, 0.875, 44, 2, 2.875, 1.5, 44, 3, 2.75, 1, 44, 4, 4.875, 1.25, 44, 5, 3.5, 1.25, 44, 6, 2.5, 1.625, 45, 1, 2.75, 1.25, 45, 2, 4.5, 1.5, 45, 3, 3.25, 1.25, 45, 4, 2.25, 1.625, 45, 5, 2, 1.875, 45, 6, 2.375, 2.375, 46, 1, 3.125, 1.25, 46, 2, 3.125, 1.5, 46, 3, 3.75, 1.25, 46, 4, 2.25, 1.25, 46, 5, 3.25, 2, 46, 6, 1.125, 1.625, 47, 1, 4.125, 2, 47, 2, 3.375, 1.25, 47, 3, 3.25, 1.125, 47, 4, 1.75, 1.625, 47, 5, 2.5, 1.375, 47, 6, 1.25, 1.625, 48, 1, 4, 1.25, 48, 2, 4.375, 1.625, 48, 3, 2.5, 1, 48, 4, 3.125, 1.625, 48, 5, 2.125, 0.375, 48, 6, 1.5, 1.375, 49, 1, 3.125, 1.25, 49, 2, 2.875, 1, 49, 3, 2.625, 0.875, 49, 4, 3.375, 1.125, 49, 5, 2.5, 1.625, 49, 6, 2.5, 1.75, 50, 1, 3.875, 1.375, 50, 2, 3.25, 1.75, 50, 3, 3.25, 1.375, 50, 4, 2.625, 1.75, 50, 5, 2.25, 1.5, 50, 6, 2.125, 1.625, 51, 1, 4.5, 1.5, 51, 2, 4.375, 0.375, 51, 3, 3.125, 1.25, 51, 4, 3, 1.875, 51, 5, 2.625, 1.25, 51, 6, 1.625, 2.375, 52, 1, 3.125, 0.75, 52, 2, 2.625, 1, 52, 3, 3.125, 1.875, 52, 4, 3.125, 0.875, 52, 5, 2, 1.5, 52, 6, 3.125, 1.75, 53, 1, 4.5, 2.125, 53, 2, 2.75, 1.5, 53, 3, 2.5, 1.875, 53, 4, 2.625, 1.25, 53, 5, 2.75, 1.875, 53, 6, 3.625, 1.875, 54, 1, 4, 1.125, 54, 2, 3.5, 0.625, 54, 3, 2.125, 2, 54, 4, 2.875, 1.875, 54, 5, 2.625, 1.875, 54, 6, 3, 1.75, 55, 1, 2.25, 2, 55, 2, 3, 1, 55, 3, 2.25, 1.125, 55, 4, 3, 1.25, 55, 5, 1.625, 1.625, 55, 6, 2.5, 2.5, 56, 1, 3.875, 0.625, 56, 2, 3.125, 1.125, 56, 3, 3.5, 1.625, 56, 4, 2.25, 1.25, 56, 5, 3, 1, 56, 6, 2.75, 1.875, 57, 1, 3.25, 1.375, 57, 2, 4.25, 1.625, 57, 3, 1.25, 1.25, 57, 4, 3.875, 0.75, 57, 5, 2.125, 1.625, 57, 6, 2.125, 1.5, 58, 1, 3.5, 2.25, 58, 2, 1.75, 1.875, 58, 3, 4.625, 1, 58, 4, 3, 1.625, 58, 5, 1.25, 1.5, 58, 6, 2, 1.625, 59, 1, 5.125, 1.75, 59, 2, 4.375, 2.125, 59, 3, 2.875, 1.5, 59, 4, 3.5, 1, 59, 5, 2.5, 1.375, 59, 6, 1.5, 1.875, 60, 1, 4.125, 0.75, 60, 2, 3.625, 1.375, 60, 3, 2.5, 1, 60, 4, 3.375, 1, 60, 5, 2.75, 1.375, 60, 6, 3, 1.75};
const double gbeta1[4] = {0.80100000486680367, 0.56042857021994985, 3.0319999971898604, -0.20842857095861347};
const double gS1[4] = {0.37825571664821439, 0.00024595137420400383, 0.00024595137420400383, 0.13980738186464992};
const double gbeta2[4] = {3.9123809532709211, -0.27632653092317683, 1.1271428564769812, 0.075918367410283516};
const double gS2[4] = {0.4433427274384667, 0.0033103132406257889, 0.0033103132406257889, 0.16560872202380714};
const double gprob[2] = {0.41666666666666669, 0.58333333333333337};
const double gpost[10] = {1, 1, 1, 1, 1, 0, 0, 0, 0, 0};
const double gll = -567.68434587517129;

struct rng {
  std::uint64_t s;
  double unif() {
    s ^= s >> 12;
    s ^= s << 25;
    s ^= s >> 27;
    return static_cast<double>((s * 2685821657736338717ULL) >> 11) * 0x1p-53;
  }
  double norm() {
    const double a = 1.0 - unif(), b = unif();
    return std::sqrt(-2.0 * std::log(a)) * std::cos(6.283185307179586 * b);
  }
};

double Y[60 * 6 * 2], X[12];
void load() {
  for (int row = 0; row < 360; ++row) {
    const int i = static_cast<int>(D[4 * row]) - 1, t = static_cast<int>(D[4 * row + 1]) - 1;
    Y[i + 60 * (t + 6 * 0)] = D[4 * row + 2];
    Y[i + 60 * (t + 6 * 1)] = D[4 * row + 3];
  }
  for (int t = 0; t < 6; ++t) {
    X[t] = 1;
    X[t + 6] = t + 1;
  }
}

// the bivariate normal mixture log-likelihood written out directly
double brute(const double *y, int n, int nt, const double *beta, const double *sig,
             const double *prob) {
  double ll = 0;
  for (int i = 0; i < n; ++i) {
    double li = 0;
    for (int k = 0; k < 2; ++k) {
      const double *b = beta + 4 * k, *s = sig + 4 * k;
      const double det = s[0] * s[3] - s[1] * s[2];
      double lp = std::log(prob[k]);
      for (int t = 0; t < nt; ++t) {
        const double e1 = y[i + n * t] - (b[0] + b[1] * (t + 1));
        const double e2 = y[i + n * (t + nt)] - (b[2] + b[3] * (t + 1));
        const double q = (s[3] * e1 * e1 - 2 * s[1] * e1 * e2 + s[0] * e2 * e2) / det;
        lp += -std::log(2 * 3.141592653589793) - 0.5 * std::log(det) - 0.5 * q;
      }
      li += std::exp(lp);
    }
    ll += std::log(li);
  }
  return ll;
}

void t_mvt() {
  load();
  static double w[2 * 400], work[64];
  double beta[8], sig[8], prob[2] = {gprob[0], gprob[1]}, ll;
  int it;
  for (int j = 0; j < 4; ++j) {
    beta[j] = gbeta1[j];
    beta[4 + j] = gbeta2[j];
    sig[j] = gS1[j];
    sig[4 + j] = gS2[j];
  }
  // 1. parity with gbmt at its estimates: log-likelihood and posteriors
  TK("mvt.ll.rc", clab_mvtll(Y, 60, 6, 2, X, 2, 2, beta, sig, prob, w, &ll, work) == CLAB_OK);
  TKN("mvt.ll.gbmt", ll / gll, 1, 1e-12);
  for (int i = 0; i < 5; ++i)
    for (int k = 0; k < 2; ++k) TKN("mvt.post", w[i + 60 * k], gpost[i + 5 * k], 1e-9);
  TKN("mvt.ll.brute", ll / brute(Y, 60, 6, beta, sig, prob), 1, 1e-12);
  // 2. EM from there: not lower (gbmt stops at tol 1e-4), and the closed
  // M-step reproduces its estimates closely
  {
    double b2[8], s2[8], p2[2] = {prob[0], prob[1]};
    for (int j = 0; j < 8; ++j) {
      b2[j] = beta[j];
      s2[j] = sig[j];
    }
    TK("mvt.em.rc", clab_mvtem(Y, 60, 6, 2, X, 2, 2, b2, s2, p2, w, &ll, 500, 1e-14, &it, work) == CLAB_OK);
    TK("mvt.em.notlower", ll >= gll - 1e-9 * std::fabs(gll));
    for (int j = 0; j < 8; ++j) TKN("mvt.em.beta", b2[j], beta[j], 1e-6);
    for (int j = 0; j < 8; ++j) TKN("mvt.em.sig", s2[j], sig[j], 1e-6);
    // stationary: central differences of the log-likelihood in beta
    double gmax = 0;
    for (int j = 0; j < 8; ++j) {
      double bb[8], lp, lm;
      for (int l = 0; l < 8; ++l) bb[l] = b2[l];
      bb[j] = b2[j] + 1e-6;
      clab_mvtll(Y, 60, 6, 2, X, 2, 2, bb, s2, p2, w, &lp, work);
      bb[j] = b2[j] - 1e-6;
      clab_mvtll(Y, 60, 6, 2, X, 2, 2, bb, s2, p2, w, &lm, work);
      gmax = std::fabs(lp - lm) / 2e-6 > gmax ? std::fabs(lp - lm) / 2e-6 : gmax;
    }
    TK("mvt.stationary", gmax < 1e-4);
  }
  // 3. EM is monotone from a poor start
  {
    double b2[8] = {2, 0, 2, 0, 2.5, 0, 1.5, 0}, s2[8] = {1, 0, 0, 1, 1, 0, 0, 1}, p2[2] = {0.5, 0.5};
    double prev = -HUGE_VAL;
    bool mono = true;
    for (int r = 0; r < 100; ++r) {
      clab_mvtem(Y, 60, 6, 2, X, 2, 2, b2, s2, p2, w, &ll, 1, 1e-300, &it, work);
      mono = mono && ll >= prev - 1e-10 * std::fabs(ll);
      prev = ll;
    }
    TK("mvt.monotone", mono);
    TKN("mvt.from.poor", ll / gll, 1, 1e-9);  // reaches the same maximum
  }
  // 4. scale equivariance (the disagreement): gbmt 0.1.4 clamps each unit's
  // log density to [-745, 709.8] before exp, so data in large units (x 1e5,
  // 20 periods x 4 indicators, dev/doss/gbmt.R) give posteriors 0.5 for
  // every unit and 50% agreement with the truth. Here multiplying the data
  // by 1e5 (and the parameters accordingly) leaves the posteriors unchanged
  // and shifts the log-likelihood by -(observed cells) log 1e5 exactly.
  {
    static double y5[720];
    for (int j = 0; j < 720; ++j) y5[j] = Y[j] * 1e5;
    double b5[8], s5[8], l5;
    static double w5[120];
    for (int j = 0; j < 8; ++j) {
      b5[j] = beta[j] * 1e5;
      s5[j] = sig[j] * 1e10;
    }
    clab_mvtll(Y, 60, 6, 2, X, 2, 2, beta, sig, prob, w, &ll, work);
    clab_mvtll(y5, 60, 6, 2, X, 2, 2, b5, s5, prob, w5, &l5, work);
    TKN("mvt.scale.ll", l5 - ll, -720 * std::log(1e5), 1e-9);
    double dmax = 0;
    for (int j = 0; j < 120; ++j) dmax = std::fabs(w5[j] - w[j]) > dmax ? std::fabs(w5[j] - w[j]) : dmax;
    TK("mvt.scale.post", dmax < 1e-9);
    // and a fit on large-unit data of the probe's shape: 40 units, 20
    // periods, 4 indicators, x 1e5, classifies every unit
    static double yb[40 * 20 * 4], xb[40], wb[80];
    rng r{6};
    for (int t = 0; t < 20; ++t) {
      xb[t] = 1;
      xb[t + 20] = t + 1;
    }
    for (int i = 0; i < 40; ++i)
      for (int t = 0; t < 20; ++t)
        for (int j = 0; j < 4; ++j) {
          const double mu = i < 20 ? (j + 1) + 0.1 * (t + 1) : (4 - j) - 0.1 * (t + 1);
          yb[i + 40 * (t + 20 * j)] = 1e5 * (mu + r.norm());
        }
    double bb[16], sb[32], pb[2] = {0.5, 0.5}, lb;
    for (int k = 0; k < 2; ++k) {
      for (int j = 0; j < 4; ++j) {
        bb[8 * k + 2 * j] = 1e5 * (k == 0 ? (j + 1) : (4 - j)) + 3e4;  // offset start
        bb[8 * k + 2 * j + 1] = 0;
      }
      for (int j = 0; j < 16; ++j) sb[16 * k + j] = (j % 5 == 0) ? 1e10 : 0;
    }
    TK("mvt.big.rc", clab_mvtem(yb, 40, 20, 4, xb, 2, 2, bb, sb, pb, wb, &lb, 200, 1e-12, &it, work) == CLAB_OK);
    int agree = 0;
    for (int i = 0; i < 40; ++i) agree += (wb[i] > 0.5) == (i < 20);
    TK("mvt.big.classify", agree == 40 || agree == 0);
  }
  // 5. a period with a missing indicator drops out
  {
    static double y2[720];
    for (int j = 0; j < 720; ++j) y2[j] = Y[j];
    y2[3 + 60 * (2 + 6 * 1)] = std::nan("");
    double l1;
    clab_mvtll(y2, 60, 6, 2, X, 2, 2, beta, sig, prob, w, &l1, work);
    TK("mvt.missing", std::isfinite(l1) && l1 > ll);
  }
  // 6. guards
  {
    double bad[8] = {1, 2, 2, 1, 1, 0, 0, 1};  // first covariance not positive definite
    TK("mvt.guard.pd", clab_mvtll(Y, 60, 6, 2, X, 2, 2, beta, bad, prob, w, &ll, work) == CLAB_EDOM);
    const double pr[2] = {0.3, 0.3};
    TK("mvt.guard.prob", clab_mvtll(Y, 60, 6, 2, X, 2, 2, beta, sig, pr, w, &ll, work) == CLAB_EDOM);
    TK("mvt.guard.r", clab_mvtll(Y, 60, 6, 0, X, 2, 2, beta, sig, prob, w, &ll, work) == CLAB_EARG);
  }
}

// the dossier's sizes: 200 units x 8 periods x 2 indicators, ng = 2, d = 1
// (gbmt 4.667 s, one run with its Ward start, zeus); and 1,000 units,
// ng = 3, d = 2 (gbmt did not finish in 600 s). Here one EM run from a
// rough start, tol 1e-10.
void b_mvt() {
  static double y[1000 * 8 * 2], x[24], w[3000], work[64];
  rng r{1};
  for (int t = 0; t < 8; ++t) {
    x[t] = 1;
    x[t + 8] = t + 1;
    x[t + 16] = (t + 1) * (t + 1);
  }
  for (int sz = 0; sz < 2; ++sz) {
    const int n = sz == 0 ? 200 : 1000, ng = sz == 0 ? 2 : 3, p = sz == 0 ? 2 : 3;
    for (int i = 0; i < n; ++i) {
      const int g = 1 + i % ng;
      for (int t = 0; t < 8; ++t) {
        y[i + n * t] = g * 0.5 * (t + 1) + r.norm();
        y[i + n * (t + 8)] = -g * 0.5 * (t + 1) + r.norm();
      }
    }
    double tm[5], ll = 0;
    int it = 0;
    for (int rep = 0; rep < 5; ++rep) {
      double beta[18] = {0}, sig[12], prob[3];
      for (int k = 0; k < ng; ++k) {
        prob[k] = 1.0 / ng;
        beta[p * 2 * k + 1] = 0.4 * (k + 1);
        beta[p * 2 * k + p + 1] = -0.4 * (k + 1);
        sig[4 * k] = sig[4 * k + 3] = 2;
        sig[4 * k + 1] = sig[4 * k + 2] = 0;
      }
      const auto a = std::chrono::steady_clock::now();
      clab_mvtem(y, n, 8, 2, x, p, ng, beta, sig, prob, w, &ll, 1000, 1e-10, &it, work);
      tm[rep] = std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count();
    }
    std::printf("BENCH clab_mvtem n=%d ng=%d d=%d: median %.6f s (k=5), %d iterations\n", n, ng, p - 1,
                tk::med5(tm), it);
  }
}

}  // namespace

TK_FAMILY("mvt", t_mvt, b_mvt);

#endif
