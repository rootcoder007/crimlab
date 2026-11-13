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
  const clab_lrn l = {kind, 1.0, 20, 3, 7, 50, 0, 0, 1e-7, 1, 0, nullptr, 0, 0, 0};
  long niw;
  const long nw = clab_dmlxsz(model, n, P, nf, nrep, &l, &l, &niw);
  if (nw > 400000 || niw > 100000) return -1;
  return clab_dmlx(model, sc, X, n, P, Y, D, Z, S, F, nf, nrep, clu, nc, 0, &l, &l,
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
  // the experimental scores: no propensity, normalised by the shares
  run(6, 1, 0);
  {
    double md = 0.0;
    for (int i = 0; i < N; ++i) md += D[i] / N;
    for (int i = 0; i < N; ++i) {
      const double q = D[i] / md;
      pa[i] = -1.0;
      pb[i] = (q - 1) * p0[i] + (1 - q) * p1[i] + (q - (1 - D[i]) / (1 - md)) * (Y[i] - p0[i]);
    }
    dmlvar(pa, pb, N, &th, &se);
    tk_eq(t, "dml.didx.th", EST[0], th, 1e-12);
    tk_eq(t, "dml.didx.se", EST[1], se, 1e-12);
    tk_ieq(t, "dml.didx.k", clab_dmlxk(6, 1) * 10 + clab_dmlxk(7, 1), 24);
  }
  run(7, 1, 0);
  {
    double c[4] = {0, 0, 0, 0};
    for (int i = 0; i < N; ++i) c[static_cast<int>(2 * D[i] + S[i])] += 1.0 / N;
    const double *g[4] = {p0, p1, p2, p3};
    const double sg[4] = {1, -1, -1, 1};
    for (int i = 0; i < N; ++i) {
      const int k = static_cast<int>(2 * D[i] + S[i]);
      pa[i] = -1.0;
      pb[i] = p3[i] - p2[i] - p1[i] + p0[i] + sg[k] * (Y[i] - g[k][i]) / c[k];
    }
    dmlvar(pa, pb, N, &th, &se);
    tk_eq(t, "dml.didcsx.th", EST[0], th, 1e-12);
    tk_eq(t, "dml.didcsx.se", EST[1], se, 1e-12);
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
// two-way clusters: 40 clusters of one kind (i / 10), 20 of the other
// (i % 20), folds of each kind; theta and the variance written out block
// by block as DoubleML's two-way form
void t_two(tk &t) {
  static int cl2[2 * N], f2[2 * N];
  for (int i = 0; i < N; ++i) {
    cl2[i] = i / 10;
    cl2[N + i] = i % 20;
    f2[i] = (i / 10) % 5;
    f2[N + i] = (i % 20) % 5;
  }
  const clab_lrn l = {0, 1.0, 20, 3, 7, 0, 0, 0, 1e-7, 1, 0, nullptr, 0, 0, 0};
  rng g = {41};
  for (int k = 0; k < 2 * N; ++k) PR[k] = runi(g);
  tk_ieq(t, "two.rc",
         clab_dmlx(0, 0, X, N, P, Y, D, Z, S, f2, 5, 1, cl2, 40, 20, &l, &l, 0.0, 0, PR, EST,
                   PSI, W, IW),
         CLAB_OK);
  const double ki = 8.0, kj = 4.0;  // clusters of each kind per fold
  double sa = 0.0, sb = 0.0, pa[N], pb[N];
  for (int i = 0; i < N; ++i) {
    const double v = D[i] - PR[N + i];
    pa[i] = -v * v;
    pb[i] = v * (Y[i] - PR[i]);
    sa += pa[i] / (ki * kj);
    sb += pb[i] / (ki * kj);
  }
  const double th = -sb / sa;
  double gam = 0.0;
  for (int k = 0; k < 5; ++k) {
    for (int l2 = 0; l2 < 5; ++l2) {
      const double c = std::min(ki, kj) / ((ki * kj) * (ki * kj));
      for (int c1 = 0; c1 < 40; ++c1) {
        if (c1 % 5 != k) continue;
        double s = 0.0;
        for (int i = 0; i < N; ++i) {
          if (i / 10 == c1 && (i % 20) % 5 == l2) s += pa[i] * th + pb[i];
        }
        gam += c * s * s;
      }
      for (int c2 = 0; c2 < 20; ++c2) {
        if (c2 % 5 != l2) continue;
        double s = 0.0;
        for (int i = 0; i < N; ++i) {
          if (i % 20 == c2 && (i / 10) % 5 == k) s += pa[i] * th + pb[i];
        }
        gam += c * s * s;
      }
    }
  }
  gam /= 25.0;
  const double j = sa / 25.0;
  tk_eq(t, "two.th", EST[0], th, 1e-12);
  tk_eq(t, "two.se", EST[1], std::sqrt(gam / (j * j) / 20.0), 1e-12);
  // learned: every test block's training rows share no cluster with it
  tk_ieq(t, "two.fit",
         clab_dmlx(0, 0, X, N, P, Y, D, Z, S, f2, 5, 1, cl2, 40, 20, &l, &l, 0.0, 1, PR, EST,
                   PSI, W, IW),
         CLAB_OK);
  tk_lt(t, "two.fin", std::fabs(EST[0]), 1e6);
  // refusals: repetitions, one kind only, an id out of range, no root
  tk_ieq(t, "two.e.rep",
         clab_dmlx(0, 0, X, N, P, Y, D, Z, S, f2, 5, 2, cl2, 40, 20, &l, &l, 0.0, 0, PR, EST,
                   PSI, W, IW),
         CLAB_EARG);
  tk_ieq(t, "two.e.nc2",
         clab_dmlx(0, 0, X, N, P, Y, D, Z, S, f2, 5, 1, nullptr, 0, 20, &l, &l, 0.0, 0, PR, EST,
                   PSI, W, IW),
         CLAB_EARG);
  tk_ieq(t, "two.e.id",
         clab_dmlx(0, 0, X, N, P, Y, D, Z, S, f2, 5, 1, cl2, 40, 10, &l, &l, 0.0, 0, PR, EST,
                   PSI, W, IW),
         CLAB_EARG);
  for (int i = 0; i < N; ++i) PR[N + i] = D[i];
  tk_ieq(t, "two.e.root",
         clab_dmlx(0, 0, X, N, P, Y, D, Z, S, f2, 5, 1, cl2, 40, 20, &l, &l, 0.0, 0, PR, EST,
                   PSI, W, IW),
         CLAB_EDOM);
}

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
  // experimental: a randomised treatment, learned nuisances
  {
    double ds[N];
    std::memcpy(ds, D, sizeof ds);
    for (int i = 0; i < N; ++i) {
      D[i] = runi(g) < 0.5 ? 1.0 : 0.0;
      Y[i] = X[i] + D[i] * 0.5 + S[i] * (0.3 + 0.5 * X[i]) + D[i] * S[i] + rnorm(g);
    }
    run(7, 1, 1);
    tk_lt(t, "sim.didcsx", std::fabs(EST[0] - 1.0), 4 * EST[1]);
    for (int i = 0; i < N; ++i) Y[i] = D[i] + 0.5 * X[i] + rnorm(g);
    run(6, 1, 1);
    tk_lt(t, "sim.didx", std::fabs(EST[0] - 1.0), 4 * EST[1]);
    std::memcpy(D, ds, sizeof ds);
  }
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
  const clab_lrn l = {0, 1.0, 20, 3, 7, 0, 0, 0, 1e-7, 1, 0, nullptr, 0, 0, 0};
  auto call = [&](int model, int sc, int n, int p, int nf, int nrep, const int *clu,
                  int nc, double trim) {
    return clab_dmlx(model, sc, X, n, p, Y, D, Z, S, F, nf, nrep, clu, nc, 0, &l, &l,
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
  tk_ieq(t, "dml.e.didx", call(6, 1, N, P, 5, 1, nullptr, 0, 0.01), CLAB_EDOM);
  tk_ieq(t, "dml.e.didcsx", call(7, 1, N, P, 5, 1, nullptr, 0, 0.01), CLAB_EDOM);
  // learning: an empty training set (nobody untreated for g0)
  const clab_lrn lg = {0, 1.0, 20, 3, 7, 0, 0, 0, 1e-7, 1, 0, nullptr, 0, 0, 0};
  for (int i = 0; i < N; ++i) D[i] = 1.0;
  tk_ieq(t, "dml.e.empty",
         clab_dmlx(2, 0, X, N, P, Y, D, Z, S, F, 5, 1, nullptr, 0, 0, &lg, &lg, 0.01, 1,
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
         clab_dmlx(0, 0, xc, N, 2, Y, D, Z, S, F, 5, 1, nullptr, 0, 0, &lg, &lg, 0.0, 1, PR,
                  EST, PSI, W, IW),
         CLAB_EDOM);
  const clab_lrn le = {1, 1.0, 0, 3, 7, 0, 0, 0, 1e-7, 1, 0, nullptr, 0, 0, 0}, la = {1, 2.0, 5, 3, 7, 0, 0, 0, 1e-7, 1, 0, nullptr, 0, 0, 0},
                 lf = {2, 1.0, 5, 3, 7, 0, 0, 0, 1e-7, 1, 0, nullptr, 0, 0, 0},
                 lt = {1, 1.0, 5, 3, 7, 0, 0, 0, 0.0, 1, 0, nullptr, 0, 0, 0};
  tk_ieq(t, "dml.e.net",
         clab_dmlx(0, 0, X, N, P, Y, D, Z, S, F, 5, 1, nullptr, 0, 0, &lg, &le, 0.0, 1, PR,
                  EST, PSI, W, IW),
         CLAB_EARG);
  tk_ieq(t, "dml.e.tree",
         clab_dmlx(0, 0, X, N, P, Y, D, Z, S, F, 5, 1, nullptr, 0, 0, &lf, &lg, 0.0, 1, PR,
                   EST, PSI, W, IW),
         CLAB_EARG);
  tk_ieq(t, "dml.e.tol",
         clab_dmlx(0, 0, X, N, P, Y, D, Z, S, F, 5, 1, nullptr, 0, 0, &lt, &lg, 0.0, 1, PR,
                   EST, PSI, W, IW),
         CLAB_EDOM);
  tk_ieq(t, "dml.e.alpha",
         clab_dmlx(0, 0, X, N, P, Y, D, Z, S, F, 5, 1, nullptr, 0, 0, &la, &lg, 0.0, 1, PR,
                  EST, PSI, W, IW),
         CLAB_EDOM);
}

// clab_kfold and clab_stack: their refusals, and a stack of two least
// squares learners, whose weights are the scaled NNLS fit.
void t_stk(tk &t) {
  int f[N], iw[N];
  tk_ieq(t, "stk.kfold.n", clab_kfold(0, 2, 1, f, iw), CLAB_EARG);
  tk_ieq(t, "stk.kfold.k", clab_kfold(3, 4, 1, f, iw), CLAB_EARG);
  tk_ieq(t, "stk.kfold", clab_kfold(N, 5, 1, F, iw), CLAB_OK);
  int cnt[5] = {0, 0, 0, 0, 0};
  for (int i = 0; i < N; ++i) ++cnt[F[i]];
  tk_true(t, "stk.kfold.size", cnt[0] == N / 5 && cnt[4] == N / 5);
  const clab_lrn g = {0, 1.0, 20, 3, 7, 0, 0, 0, 1e-7, 1, 0, nullptr, 0, 0, 0};
  const clab_lrn two[2] = {g, g};
  const clab_lrn s = {3, 1.0, 20, 3, 7, 0, 0, 0, 1e-7, 1, 0, two, 2, 5, 0};
  double oof[2 * N], rk[2], wt[2], pr[N];
  long niw;
  const long nw = clab_stksz(N, P, N, &s, &niw);
  tk_true(t, "stk.size", nw <= 400000 && niw <= 100000);
  tk_ieq(t, "stk.ok", clab_stack(X, N, P, Y, 0, &s, F, X, N, oof, rk, wt, pr, W, IW), CLAB_OK);
  tk_true(t, "stk.wt", wt[0] >= 0.0 && wt[1] >= 0.0 && std::fabs(wt[0] + wt[1] - 1.0) < 1e-15);
  auto bad = [&](const char *nm, const clab_lrn &l, int n, int fam, int nt, int code) {
    tk_ieq(t, nm, clab_stack(X, n, P, Y, fam, &l, F, X, nt, oof, rk, wt, pr, W, IW), code);
  };
  bad("stk.e.kind", g, N, 0, N, CLAB_EARG);
  bad("stk.e.n", s, 1, 0, N, CLAB_EARG);
  bad("stk.e.fam", s, N, 2, N, CLAB_EARG);
  bad("stk.e.nt", s, N, 0, -1, CLAB_EARG);
  bad("stk.e.big", s, N, 0, 2147483647, CLAB_EARG);
  clab_lrn e = s;
  e.sub = nullptr;
  bad("stk.e.sub", e, N, 0, N, CLAB_EARG);
  e = s;
  e.nsub = 0;
  bad("stk.e.nsub", e, N, 0, N, CLAB_EARG);
  e.nsub = 17;
  bad("stk.e.nsub2", e, N, 0, N, CLAB_EARG);
  e = s;
  e.nstk = 1;
  bad("stk.e.nstk", e, N, 0, N, CLAB_EARG);
  e = s;
  e.disc = 2;
  bad("stk.e.disc", e, N, 0, N, CLAB_EARG);
  e.disc = -1;
  bad("stk.e.disc2", e, N, 0, N, CLAB_EARG);
  const clab_lrn nest[1] = {s};
  e = s;
  e.sub = nest;
  e.nsub = 1;
  bad("stk.e.nest", e, N, 0, N, CLAB_EARG);
  clab_lrn ga = g;
  ga.alpha = 2.0;
  const clab_lrn sa[1] = {ga};
  e = s;
  e.sub = sa;
  e.nsub = 1;
  bad("stk.e.alpha", e, N, 0, N, CLAB_EDOM);
  tk_ieq(t, "stk.e.dml",
         clab_dmlx(0, 0, X, N, P, Y, D, Z, S, F, 5, 1, nullptr, 0, 0, &e, &g, 0.0, 1, PR,
                   EST, PSI, W, IW),
         CLAB_EDOM);
  F[0] = 5;
  bad("stk.e.fold", s, N, 0, N, CLAB_EARG);
  F[0] = 0;
  const double y0 = Y[0];
  Y[0] = 2.0;
  bad("stk.e.bin", s, N, 1, N, CLAB_EDOM);
  Y[0] = std::nan("");
  bad("stk.e.y", s, N, 0, N, CLAB_EDOM);
  Y[0] = y0;
  const double x0 = X[0];
  X[0] = std::nan("");
  bad("stk.e.x", s, N, 0, N, CLAB_EDOM);
  X[0] = x0;
  // three rows for five folds; a fold that holds every row
  bad("stk.e.few", s, 3, 0, 1, CLAB_EDOM);
  for (int i = 0; i < N; ++i) f[i] = F[i];
  for (int i = 0; i < N; ++i) F[i] = 0;
  bad("stk.e.one", s, N, 0, N, CLAB_EDOM);
  for (int i = 0; i < N; ++i) F[i] = f[i];
}

}  // namespace

void t_dml(tk &t, int) {
  t_forms(t);
  t_two(t);
  t_sim(t);
  t_bad(t);
  t_stk(t);
}
