// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of the cross-fitted estimators (src/dml.cpp): each score's
// estimate and variance from given nuisance predictions against the
// formula written out here; the aggregation over repetitions; the
// clustered variance; the targeting step of TMLE against a bisection;
// known-truth simulations with learned nuisances; every refusal. Built
// with the self-test.
#include <algorithm>
#include <cmath>
#include <cstring>

#include "crimlab.h"
#include "tchk.h"

namespace {

const int N = 400, P = 4, R = 3;
double X[N * P], Y[N], D[N], Z[N], S[N], PR[5 * N * R], EST[2 + 2 * R],
    PSI[2 * N * R];
int F[N * R], CL[N];
double W[400000];
int IW[100000];

double rnorm(rng &g) {
  return std::sqrt(-2.0 * std::log(runi(g))) * std::cos(6.283185307179586 * runi(g));
}

int run(int model, int sc, int fit, int nrep = 1, const int *clu = nullptr,
        int nc = 0, double trim = 0.0, int kind = 0, int n = N, int nf = 5) {
  const clab_lrn l = {kind, 1.0, 20, 3, 7, 50, 0, 0};
  long niw;
  const long nw = clab_dmlxsz(model, n, P, nf, nrep, &l, &l, &niw);
  if (nw > 400000 || niw > 100000) return -1;
  return clab_dmlx(model, sc, X, n, P, Y, D, Z, S, F, nf, nrep, clu, nc, &l, &l,
                  trim, fit, PR, EST, PSI, W, IW);
}

// theta and its DoubleML variance from psi_a, psi_b
void dmlvar(const double *pa, const double *pb, int n, double *th, double *se) {
  double sa = 0.0, sb = 0.0;
  for (int i = 0; i < n; ++i) {
    sa += pa[i];
    sb += pb[i];
  }
  *th = -sb / sa;
  double s2 = 0.0;
  for (int i = 0; i < n; ++i) s2 += (pa[i] * *th + pb[i]) * (pa[i] * *th + pb[i]);
  *se = std::sqrt(s2 / n / ((sa / n) * (sa / n)) / n);
}

// the closed forms, nuisance predictions given
void t_forms(tk &t) {
  rng g = {5};
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < P; ++j) X[i + j * N] = rnorm(g);
    D[i] = runi(g) < 0.4 ? 1.0 : 0.0;
    Z[i] = runi(g) < 0.5 ? 1.0 : 0.0;
    S[i] = runi(g) < 0.7 ? 1.0 : 0.0;
    Y[i] = 1.0 + D[i] + rnorm(g);
    F[i] = i % 5;
  }
  for (int k = 0; k < 5 * N; ++k) PR[k] = 0.2 + 0.6 * runi(g);
  const double *p0 = PR, *p1 = PR + N, *p2 = PR + 2 * N, *p3 = PR + 3 * N,
               *p4 = PR + 4 * N;
  double pa[N], pb[N], th, se;
  // PLR, partialling out and IV-type (the treatment taken continuous)
  for (int sc = 0; sc < 2; ++sc) {
    tk_ieq(t, "dml.plr.rc", run(0, sc, 0), CLAB_OK);
    for (int i = 0; i < N; ++i) {
      const double v = D[i] - p1[i];
      pa[i] = sc ? -v * D[i] : -v * v;
      pb[i] = sc ? v * (Y[i] - p2[i]) : v * (Y[i] - p0[i]);
    }
    dmlvar(pa, pb, N, &th, &se);
    tk_eq(t, "dml.plr.th", EST[0], th, 1e-12);
    tk_eq(t, "dml.plr.se", EST[1], se, 1e-12);
    tk_veq(t, "dml.plr.psi", PSI + N, pb, N, 1e-14);
  }
  // PLIV
  for (int sc = 0; sc < 2; ++sc) {
    run(1, sc, 0);
    for (int i = 0; i < N; ++i) {
      const double v = Z[i] - p1[i];
      pa[i] = sc ? -D[i] * v : -(D[i] - p2[i]) * v;
      pb[i] = sc ? v * (Y[i] - p3[i]) : v * (Y[i] - p0[i]);
    }
    dmlvar(pa, pb, N, &th, &se);
    tk_eq(t, "dml.pliv.th", EST[0], th, 1e-12);
    tk_eq(t, "dml.pliv.se", EST[1], se, 1e-12);
  }
  // IRM: ATE, and ATTE with each test fold's treated share; trimming
  for (int sc = 0; sc < 2; ++sc) {
    run(2, sc, 0, 1, nullptr, 0, 0.3);
    double q[5] = {0, 0, 0, 0, 0}, c[5] = {0, 0, 0, 0, 0};
    for (int i = 0; i < N; ++i) {
      q[F[i]] += D[i];
      c[F[i]] += 1;
    }
    for (int i = 0; i < N; ++i) {
      const double m = std::min(std::max(p2[i], 0.3), 0.7), h = q[F[i]] / c[F[i]];
      pa[i] = sc ? -D[i] / h : -1.0;
      pb[i] = sc ? D[i] * (Y[i] - p0[i]) / h - m * (1 - D[i]) * (Y[i] - p0[i]) / (h * (1 - m))
                 : p1[i] - p0[i] + D[i] * (Y[i] - p1[i]) / m - (1 - D[i]) * (Y[i] - p0[i]) / (1 - m);
    }
    dmlvar(pa, pb, N, &th, &se);
    tk_eq(t, "dml.irm.th", EST[0], th, 1e-12);
    tk_eq(t, "dml.irm.se", EST[1], se, 1e-12);
  }
  // AIPW: the IRM ATE, its standard error from the standard deviation
  {
    run(2, 0, 0, 1, nullptr, 0, 0.3);
    std::memcpy(pb, PSI + N, sizeof pb);
    run(4, 0, 0, 1, nullptr, 0, 0.3);
    double m = 0.0, v = 0.0;
    for (int i = 0; i < N; ++i) m += pb[i] / N;
    for (int i = 0; i < N; ++i) v += (pb[i] - m) * (pb[i] - m);
    tk_eq(t, "dml.aipw.th", EST[0], m, 1e-12);
    tk_eq(t, "dml.aipw.se", EST[1], std::sqrt(v / (N - 1) / N), 1e-12);
  }
  // IIVM
  run(3, 0, 0, 1, nullptr, 0, 0.1);
  for (int i = 0; i < N; ++i) {
    const double m = std::min(std::max(p2[i], 0.1), 0.9);
    pb[i] = p1[i] - p0[i] + Z[i] * (Y[i] - p1[i]) / m - (1 - Z[i]) * (Y[i] - p0[i]) / (1 - m);
    pa[i] = -(p4[i] - p3[i] + Z[i] * (D[i] - p4[i]) / m - (1 - Z[i]) * (D[i] - p3[i]) / (1 - m));
  }
  dmlvar(pa, pb, N, &th, &se);
  tk_eq(t, "dml.iivm.th", EST[0], th, 1e-12);
  tk_eq(t, "dml.iivm.se", EST[1], se, 1e-12);
  // DiD on a panel: in-sample normalised weights
  run(6, 0, 0, 1, nullptr, 0, 0.01);
  {
    double md = 0.0, mw = 0.0;
    for (int i = 0; i < N; ++i) {
      md += D[i] / N;
      mw += (1 - D[i]) * p1[i] / (1 - p1[i]) / N;
    }
    for (int i = 0; i < N; ++i) {
      pa[i] = -D[i] / md;
      pb[i] = (D[i] / md - (1 - D[i]) * p1[i] / (1 - p1[i]) / mw) * (Y[i] - p0[i]);
    }
    dmlvar(pa, pb, N, &th, &se);
    tk_eq(t, "dml.did.th", EST[0], th, 1e-12);
    tk_eq(t, "dml.did.se", EST[1], se, 1e-12);
  }
  // DiD on repeated cross-sections (S is the period)
  run(7, 0, 0, 1, nullptr, 0, 0.01);
  {
    double md = 0, a11 = 0, a10 = 0, a01 = 0, a00 = 0;
    for (int i = 0; i < N; ++i) {
      const double r = p4[i] / (1 - p4[i]);
      md += D[i] / N;
      a11 += D[i] * S[i] / N;
      a10 += D[i] * (1 - S[i]) / N;
      a01 += (1 - D[i]) * S[i] * r / N;
      a00 += (1 - D[i]) * (1 - S[i]) * r / N;
    }
    for (int i = 0; i < N; ++i) {
      const double r = p4[i] / (1 - p4[i]);
      pa[i] = -D[i] / md;
      pb[i] = D[i] / md * (p3[i] - p2[i] - p1[i] + p0[i]) + D[i] * S[i] / a11 * (Y[i] - p3[i]) -
              D[i] * (1 - S[i]) / a10 * (Y[i] - p2[i]) -
              (1 - D[i]) * S[i] * r / a01 * (Y[i] - p1[i]) +
              (1 - D[i]) * (1 - S[i]) * r / a00 * (Y[i] - p0[i]);
    }
    dmlvar(pa, pb, N, &th, &se);
    tk_eq(t, "dml.didcs.th", EST[0], th, 1e-12);
    tk_eq(t, "dml.didcs.se", EST[1], se, 1e-12);
  }
  // sample selection, missing at random
  run(8, 0, 0, 1, nullptr, 0, 0.05);
  for (int i = 0; i < N; ++i) {
    const double m = std::min(std::max(p1[i], 0.05), 0.95);
    pa[i] = -1.0;
    pb[i] = D[i] * S[i] * (Y[i] - p3[i]) / (m * p0[i]) + p3[i] -
            (1 - D[i]) * S[i] * (Y[i] - p2[i]) / ((1 - m) * p0[i]) - p2[i];
  }
  dmlvar(pa, pb, N, &th, &se);
  tk_eq(t, "dml.ssm.th", EST[0], th, 1e-12);
  tk_eq(t, "dml.ssm.se", EST[1], se, 1e-12);
  // TMLE: the fluctuation found again by bisection on its score
  {
    tk_ieq(t, "dml.tmle.rc", run(5, 0, 0, 1, nullptr, 0, 0.25), CLAB_OK);
    double a = Y[0], b = Y[0];
    for (int i = 0; i < N; ++i) {
      a = std::min(a, Y[i]);
      b = std::max(b, Y[i]);
    }
    auto lq = [&](double q) {
      q = std::min(std::max(q, a - 0.01 * std::fabs(a)), b + 0.01 * std::fabs(b));
      const double u = std::min(std::max((q - a) / (b - a), 5e-4), 1 - 5e-4);
      return std::log(u / (1 - u));
    };
    double eps[2];
    for (int arm = 0; arm < 2; ++arm) {
      auto sc = [&](double e) {
        double u = 0.0;
        for (int i = 0; i < N; ++i) {
          if (D[i] != arm) continue;
          const double gg = std::max(arm ? p2[i] : 1 - p2[i], 0.25);
          u += ((Y[i] - a) / (b - a) - 1 / (1 + std::exp(-(lq(arm ? p1[i] : p0[i]) + e / gg)))) / gg;
        }
        return u;
      };
      double lo = -50, hi = 50;
      for (int it = 0; it < 200; ++it) {
        const double mid = 0.5 * (lo + hi);
        (sc(mid) > 0 ? lo : hi) = mid;
      }
      eps[arm] = 0.5 * (lo + hi);
    }
    double m1 = 0, m0 = 0, ic[N];
    for (int i = 0; i < N; ++i) {
      const double g1 = std::max(p2[i], 0.25), g0 = std::max(1 - p2[i], 0.25);
      const double q1 = a + (b - a) / (1 + std::exp(-(lq(p1[i]) + eps[1] / g1)));
      const double q0 = a + (b - a) / (1 + std::exp(-(lq(p0[i]) + eps[0] / g0)));
      m1 += q1 / N;
      m0 += q0 / N;
      ic[i] = (D[i] / g1 - (1 - D[i]) / g0) * (Y[i] - (D[i] ? q1 : q0)) + q1 - q0;
    }
    double mu = 0, v = 0;
    for (int i = 0; i < N; ++i) mu += ic[i] / N;
    for (int i = 0; i < N; ++i) v += (ic[i] - mu) * (ic[i] - mu);
    tk_eq(t, "dml.tmle.th", EST[0], m1 - m0, 1e-10);
    tk_eq(t, "dml.tmle.se", EST[1], std::sqrt(v / (N - 1) / N), 1e-9);
    tk_lt(t, "dml.tmle.ic", std::fabs(mu - (m1 - m0)), 1e-9);
  }
  // repetitions: the median, and the median of n se^2 + (theta_r -
  // theta)^2 over n
  {
    double ths[R], ses[R];
    for (int r = 0; r < R; ++r) {
      for (int i = 0; i < N; ++i) F[r * N + i] = (i + r) % 5;
      for (int k = 0; k < 2 * N; ++k) PR[r * 2 * N + k] = PR[k] + 0.05 * r * runi(g);
    }
    tk_ieq(t, "dml.rep.rc", run(0, 0, 0, R), CLAB_OK);
    for (int r = 0; r < R; ++r) {
      for (int i = 0; i < N; ++i) {
        const double v = D[i] - PR[r * 2 * N + N + i];
        pa[i] = -v * v;
        pb[i] = v * (Y[i] - PR[r * 2 * N + i]);
      }
      dmlvar(pa, pb, N, ths + r, ses + r);
    }
    tk_veq(t, "dml.rep.each", EST + 2, ths, R, 1e-12);
    double s[R];
    std::memcpy(s, ths, sizeof s);
    std::sort(s, s + R);
    const double md = s[1];
    for (int r = 0; r < R; ++r) s[r] = N * ses[r] * ses[r] + (ths[r] - md) * (ths[r] - md);
    std::sort(s, s + R);
    tk_eq(t, "dml.rep.th", EST[0], md, 1e-15);
    tk_eq(t, "dml.rep.se", EST[1], std::sqrt(s[1] / N), 1e-12);
    // an even number of repetitions takes the middle two
    tk_ieq(t, "dml.rep2.rc", run(0, 0, 0, 2), CLAB_OK);
    tk_eq(t, "dml.rep2.th", EST[0], 0.5 * (EST[2] + EST[3]), 1e-15);
  }
  // one-way clusters: 40 clusters of 10, each in one fold
  {
    for (int i = 0; i < N; ++i) {
      CL[i] = i / 10;
      F[i] = (i / 10) % 5;
    }
    tk_ieq(t, "dml.clu.rc", run(0, 0, 0, 1, CL, 40), CLAB_OK);
    double cs[40] = {0}, sa = 0, sb = 0;
    const double kk = 8.0;  // clusters per fold
    for (int i = 0; i < N; ++i) {
      const double v = D[i] - p1[i];
      pa[i] = -v * v;
      pb[i] = v * (Y[i] - p0[i]);
      sa += pa[i] / kk;
      sb += pb[i] / kk;
    }
    const double th0 = -sb / sa;
    for (int i = 0; i < N; ++i) cs[CL[i]] += pa[i] * th0 + pb[i];
    double gam = 0;
    for (int c = 0; c < 40; ++c) gam += cs[c] * cs[c] / kk;
    gam /= 5;
    const double j = sa / 5;
    tk_eq(t, "dml.clu.th", EST[0], th0, 1e-12);
    tk_eq(t, "dml.clu.se", EST[1], std::sqrt(gam / (j * j) / 40), 1e-12);
    for (int i = 0; i < N; ++i) F[i] = i % 5;
  }
}

// known truth: simulated data, nuisances learned
void t_sim(tk &t) {
  rng g = {17};
  // PLR: theta = 0.5, confounding through x1, x2
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < P; ++j) X[i + j * N] = rnorm(g);
    D[i] = 0.8 * X[i] - 0.5 * X[i + N] + rnorm(g);
    Y[i] = 0.5 * D[i] + X[i] + 0.7 * X[i + N] + rnorm(g);
    F[i] = i % 5;
  }
  for (int kind = 0; kind < 3; ++kind) {
    for (int sc = 0; sc < 2; ++sc) {
      tk_ieq(t, "sim.plr.rc", run(0, sc, 1, 1, nullptr, 0, 0.0, kind), CLAB_OK);
      tk_lt(t, "sim.plr", std::fabs(EST[0] - 0.5), 4 * EST[1]);
      tk_lt(t, "sim.plr.se", EST[1], 0.1);
    }
  }
  // PLIV: z shifts d; d endogenous through u
  for (int i = 0; i < N; ++i) {
    const double u = rnorm(g);
    Z[i] = X[i] + rnorm(g);
    D[i] = Z[i] + 0.5 * X[i + N] + u + rnorm(g);
    Y[i] = 0.5 * D[i] + X[i + N] + u + rnorm(g);
  }
  for (int sc = 0; sc < 2; ++sc) {
    run(1, sc, 1);
    tk_lt(t, "sim.pliv", std::fabs(EST[0] - 0.5), 4 * EST[1]);
  }
  // IRM, AIPW, TMLE: effect 1, propensity logistic in x1
  for (int i = 0; i < N; ++i) {
    D[i] = runi(g) < 1 / (1 + std::exp(-X[i])) ? 1.0 : 0.0;
    Y[i] = D[i] + X[i] + 0.5 * X[i + N] + rnorm(g);
  }
  for (int model : {2, 4, 5}) {
    tk_ieq(t, "sim.irm.rc", run(model, 0, 1, 1, nullptr, 0, model == 5 ? 0.02 : 0.01), CLAB_OK);
    tk_lt(t, "sim.irm", std::fabs(EST[0] - 1.0), 4 * EST[1]);
  }
  run(2, 1, 1, 1, nullptr, 0, 0.01, 1);
  tk_lt(t, "sim.atte", std::fabs(EST[0] - 1.0), 4 * EST[1]);
  // forests for the outcome and the propensity: probabilities in [0, 1]
  tk_ieq(t, "sim.irm.rf", run(2, 0, 1, 1, nullptr, 0, 0.01, 2), CLAB_OK);
  tk_lt(t, "sim.irm.rf2", std::fabs(EST[0] - 1.0), 4 * EST[1]);
  {
    double lo = 1, hi = 0;
    for (int i = 0; i < N; ++i) {
      lo = std::fmin(lo, PR[2 * N + i]);
      hi = std::fmax(hi, PR[2 * N + i]);
    }
    tk_true(t, "sim.irm.rf01", lo >= 0 && hi <= 1);
  }
  // DiD panel: y the change, effect 1 on the treated, trend in x1
  for (int i = 0; i < N; ++i) Y[i] = D[i] + 0.5 * X[i] + rnorm(g);
  run(6, 0, 1, 1, nullptr, 0, 0.01);
  tk_lt(t, "sim.did", std::fabs(EST[0] - 1.0), 4 * EST[1]);
  // DiD on cross-sections: S the period
  for (int i = 0; i < N; ++i) {
    S[i] = runi(g) < 0.5 ? 1.0 : 0.0;
    Y[i] = X[i] + D[i] * 0.5 + S[i] * (0.3 + 0.5 * X[i]) + D[i] * S[i] + rnorm(g);
  }
  run(7, 0, 1, 1, nullptr, 0, 0.01);
  tk_lt(t, "sim.didcs", std::fabs(EST[0] - 1.0), 4 * EST[1]);
  // selection: observed with probability rising in x2 and d
  for (int i = 0; i < N; ++i) {
    S[i] = runi(g) < 1 / (1 + std::exp(-(0.5 + X[i + N] + 0.5 * D[i]))) ? 1.0 : 0.0;
    Y[i] = D[i] + X[i] + X[i + N] + rnorm(g);
  }
  run(8, 0, 1, 1, nullptr, 0, 0.01);
  tk_lt(t, "sim.ssm", std::fabs(EST[0] - 1.0), 4 * EST[1]);
  // IIVM, one-sided noncompliance: nobody with z = 0 is treated, so r0
  // is learned from a constant and predicted as 0
  for (int i = 0; i < N; ++i) {
    Z[i] = runi(g) < 0.5 ? 1.0 : 0.0;
    D[i] = Z[i] == 1.0 && runi(g) < 1 / (1 + std::exp(-X[i])) ? 1.0 : 0.0;
    Y[i] = 2 * D[i] + X[i] + rnorm(g);
  }
  tk_ieq(t, "sim.iivm.rc", run(3, 0, 1, 1, nullptr, 0, 0.01), CLAB_OK);
  tk_lt(t, "sim.iivm", std::fabs(EST[0] - 2.0), 4 * EST[1]);
  tk_eq(t, "sim.iivm.r0", PR[3 * N + 5], 0.0, 0.0);
}

void t_bad(tk &t) {
  const clab_lrn l = {0, 1.0, 20, 3, 7, 0, 0, 0};
  auto call = [&](int model, int sc, int n, int p, int nf, int nrep, const int *clu,
                  int nc, double trim) {
    return clab_dmlx(model, sc, X, n, p, Y, D, Z, S, F, nf, nrep, clu, nc, &l, &l,
                    trim, 0, PR, EST, PSI, W, IW);
  };
  tk_ieq(t, "dml.e.model", call(9, 0, N, P, 5, 1, nullptr, 0, 0), CLAB_EARG);
  tk_ieq(t, "dml.e.model2", call(-1, 0, N, P, 5, 1, nullptr, 0, 0), CLAB_EARG);
  tk_ieq(t, "dml.e.sc", call(0, 2, N, P, 5, 1, nullptr, 0, 0), CLAB_EARG);
  tk_ieq(t, "dml.e.sc2", call(0, -1, N, P, 5, 1, nullptr, 0, 0), CLAB_EARG);
  tk_ieq(t, "dml.e.n", call(0, 0, 1, P, 5, 1, nullptr, 0, 0), CLAB_EARG);
  tk_ieq(t, "dml.e.p", call(0, 0, N, 0, 5, 1, nullptr, 0, 0), CLAB_EARG);
  tk_ieq(t, "dml.e.nf", call(0, 0, N, P, 1, 1, nullptr, 0, 0), CLAB_EARG);
  tk_ieq(t, "dml.e.nf2", call(0, 0, N, P, N + 1, 1, nullptr, 0, 0), CLAB_EARG);
  tk_ieq(t, "dml.e.rep", call(0, 0, N, P, 5, 0, nullptr, 0, 0), CLAB_EARG);
  tk_ieq(t, "dml.e.clu", call(2, 0, N, P, 5, 1, CL, 40, 0), CLAB_EARG);
  tk_ieq(t, "dml.e.nc", call(0, 0, N, P, 5, 1, CL, 1, 0), CLAB_EARG);
  tk_ieq(t, "dml.e.cid", call(0, 0, N, P, 5, 1, CL, 20, 0), CLAB_EARG);
  tk_ieq(t, "dml.e.fid", call(0, 0, N, P, 4, 1, nullptr, 0, 0), CLAB_EARG);
  tk_ieq(t, "dml.e.trim", call(2, 0, N, P, 5, 1, nullptr, 0, 0.5), CLAB_EDOM);
  tk_ieq(t, "dml.e.trim2", call(2, 0, N, P, 5, 1, nullptr, 0, -0.1), CLAB_EDOM);
  tk_ieq(t, "dml.k", clab_dmlxk(9, 0) + clab_dmlxk(-1, 0), 0);
  // no root: the treatment equals its prediction
  for (int i = 0; i < N; ++i) PR[N + i] = D[i];
  tk_ieq(t, "dml.e.root", call(0, 0, N, P, 5, 1, nullptr, 0, 0), CLAB_EDOM);
  for (int i = 0; i < N; ++i) CL[i] = i / 10;
  for (int i = 0; i < N; ++i) F[i] = (i / 10) % 5;
  tk_ieq(t, "dml.e.root2", call(0, 0, N, P, 5, 1, CL, 40, 0), CLAB_EDOM);
  for (int i = 0; i < N; ++i) F[i] = i % 5;
  // DiD with nobody treated, or a cell empty
  double dsave[N];
  std::memcpy(dsave, D, sizeof D);
  for (int i = 0; i < N; ++i) D[i] = 0.0;
  tk_ieq(t, "dml.e.did", call(6, 0, N, P, 5, 1, nullptr, 0, 0.01), CLAB_EDOM);
  tk_ieq(t, "dml.e.didcs", call(7, 0, N, P, 5, 1, nullptr, 0, 0.01), CLAB_EDOM);
  // learning: an empty training set (nobody untreated for g0)
  const clab_lrn lg = {0, 1.0, 20, 3, 7, 0, 0, 0};
  for (int i = 0; i < N; ++i) D[i] = 1.0;
  tk_ieq(t, "dml.e.empty",
         clab_dmlx(2, 0, X, N, P, Y, D, Z, S, F, 5, 1, nullptr, 0, &lg, &lg, 0.01, 1,
                  PR, EST, PSI, W, IW),
         CLAB_EDOM);
  std::memcpy(D, dsave, sizeof D);
  // a learner that fails: collinear columns for least squares, and the
  // elastic net with too few passes
  static double xc[N * 2];
  for (int i = 0; i < N; ++i) {
    xc[i] = X[i];
    xc[i + N] = 2 * X[i];
  }
  tk_ieq(t, "dml.e.lrn",
         clab_dmlx(0, 0, xc, N, 2, Y, D, Z, S, F, 5, 1, nullptr, 0, &lg, &lg, 0.0, 1, PR,
                  EST, PSI, W, IW),
         CLAB_EDOM);
  const clab_lrn le = {1, 1.0, 0, 3, 7, 0, 0, 0}, la = {1, 2.0, 5, 3, 7, 0, 0, 0},
                 lf = {2, 1.0, 5, 3, 7, 0, 0, 0};
  tk_ieq(t, "dml.e.net",
         clab_dmlx(0, 0, X, N, P, Y, D, Z, S, F, 5, 1, nullptr, 0, &lg, &le, 0.0, 1, PR,
                  EST, PSI, W, IW),
         CLAB_EARG);
  tk_ieq(t, "dml.e.tree",
         clab_dmlx(0, 0, X, N, P, Y, D, Z, S, F, 5, 1, nullptr, 0, &lf, &lg, 0.0, 1, PR,
                   EST, PSI, W, IW),
         CLAB_EARG);
  tk_ieq(t, "dml.e.alpha",
         clab_dmlx(0, 0, X, N, P, Y, D, Z, S, F, 5, 1, nullptr, 0, &la, &lg, 0.0, 1, PR,
                  EST, PSI, W, IW),
         CLAB_EDOM);
}

}  // namespace

void t_dml(tk &t, int) {
  t_forms(t);
  t_sim(t);
  t_bad(t);
}
