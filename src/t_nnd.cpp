// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests of the nearest-neighbour kernel (src/nnd.cpp).
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "cgeo.h"
#include "crimlab.h"
#include "tchk.h"

namespace {

// a message of the extended run (standalone, "ext"); silent by default
template <class... A>
void info(int big, const char *f, A... a) {
  if (big) std::printf(f, a...);
}

// Reference values from dev/doss/dst.R; the first line names the versions.
// R 4.5.0 spatstat.explore 3.8.3 spatstat.geom 3.8.3 
// colMeans(nndist(XD, k = 1:4))
const double nnmean[] = {
    0.059550270948419783, 0.089740395404999909, 0.11583584160154756,
    0.13677241103812113};
// nndist(XD, k = 1:4)[1:10, ] by row
const double nnrows[] = {
    0, 0.19975729119355101, 0.29587596855415832, 0.30015668934128381, 0,
    0.02393436967260635, 0.09608287239235512, 0.10068832420290702, 0,
    0.057982195794100556, 0.11545819691101115, 0.14939130145198631, 0,
    0.047910366785379854, 0.078647500827646127, 0.089654823532436909, 0,
    0.025179024285643805, 0.13351015350209999, 0.13637604100512685,
    0.047443554594508511, 0.078810981741796474, 0.082776327956382498,
    0.14273916562575945, 0.091803061746022527, 0.10071664895875472,
    0.11681324043382421, 0.11951813853854618, 0.02635652377877808,
    0.044839003115889411, 0.11951813853854618, 0.12017581747983559,
    0.051370930553423878, 0.11491730292260414, 0.13622015618752892,
    0.13737633182980138, 0.12948353421818412, 0.13218896058786847,
    0.14323511792587523, 0.16789350449423091};
// clarkevans(X, correction = c("none", "Donnelly"))
const double ce[] = {
    1.1067868451491931, 1.0663286425305207};

const double PI = cgeo::pi;

// pattern P of dev/doss/dst.R
void patp(double *x, double *y, int n) {
  for (long long i = 0; i < n; ++i) {
    x[i] = static_cast<double>((i * i * 7919 + i * 31) % 1009) / 1009 * 2;
    y[i] = static_cast<double>((i * i * 104729 + i * 7) % 1013) / 1013;
  }
}

// brute-force k-th nearest-neighbour distances (sorted insertion)
void brute(const double *x, const double *y, int n, int K, double *d) {
  double b[64];
  for (int i = 0; i < n; ++i) {
    int got = 0;
    for (int j = 0; j < n; ++j) {
      if (j == i) continue;
      const double dx = x[j] - x[i], dy = y[j] - y[i];
      const double v = dx * dx + dy * dy;
      if (got < K) {
        int p = got++;
        while (p > 0 && b[p - 1] > v) {
          b[p] = b[p - 1];
          --p;
        }
        b[p] = v;
      } else if (v < b[K - 1]) {
        int p = K - 1;
        while (p > 0 && b[p - 1] > v) {
          b[p] = b[p - 1];
          --p;
        }
        b[p] = v;
      }
    }
    for (int k = 0; k < K; ++k) d[i * K + k] = std::sqrt(b[k]);
  }
}

void famtest(tk &T, int big) {
  const double win[4] = {0, 2, 0, 1};
  static int iw[4000];
  static double wk[4000], dist[64000], out[200];
  // parity with spatstat on P plus 5 coincident copies
  {
    double x[155], y[155];
    patp(x, y, 150);
    for (int i = 0; i < 5; ++i) {
      x[150 + i] = x[i];
      y[150 + i] = y[i];
    }
    tk_true(T, "nn.ok", clab_nnidx(x, y, 155, win, 4, 0, 1, dist, iw, wk, out) ==
                    CLAB_OK);
    for (int k = 0; k < 4; ++k) tk_eq(T, "nn.mean", out[6 * k], nnmean[k], 1e-13);
    for (int i = 0; i < 40; ++i) tk_eq(T, "nn.rows", dist[i], nnrows[i], 1e-14);
  }
  // Clark-Evans index, naive and with Donnelly's correction (spatstat
  // clarkevans on P)
  {
    double x[150], y[150];
    patp(x, y, 150);
    clab_nnidx(x, y, 150, win, 1, 0, 1, nullptr, iw, wk, out);
    tk_eq(T, "ce.naive", out[2], ce[0], 1e-13);
    tk_eq(T, "ce.donnelly", out[7], ce[1], 1e-13);
  }
  // Thompson's constants: k (2k)! / ((2^k k!)^2)
  {
    cgeo::rng r(2);
    double x[400], y[400];
    for (int i = 0; i < 400; ++i) {
      x[i] = 2 * r.unif();
      y[i] = r.unif();
    }
    clab_nnidx(x, y, 400, win, 5, 0, 1, nullptr, iw, wk, out);
    const double lam = 200;
    const double c[5] = {0.5, 0.75, 0.9375, 1.09375, 1.23046875};
    for (int k = 0; k < 5; ++k) {
      tk_eq(T, "nn.thompson", out[6 * k + 1], c[k] / std::sqrt(lam), 1e-15);
      tk_eq(T, "nn.sd", out[6 * k + 5],
          std::sqrt((k + 1) / (PI * lam) - c[k] * c[k] / lam), 1e-14);
    }
  }
  // a square lattice: every point's nearest neighbour is at the spacing
  {
    double x[100], y[100];
    for (int i = 0; i < 100; ++i) {
      x[i] = 0.05 + 0.1 * (i % 10);
      y[i] = 0.05 + 0.1 * (i / 10);
    }
    const double w1[4] = {0, 1, 0, 1};
    clab_nnidx(x, y, 100, w1, 1, 0, 1, nullptr, iw, wk, out);
    tk_eq(T, "nn.lattice", out[0], 0.1, 1e-14);
    tk_eq(T, "nn.latR", out[2], 0.1 / (0.5 * std::sqrt(0.01)), 1e-13);
  }
  // brute force on random, clustered and collinear sets
  {
    cgeo::rng r(9);
    static double x[1000], y[1000], bd[64000];
    for (int rep = 0; rep < 12; ++rep) {
      const int n = 20 + rep * 80, K = 1 + rep % 6;
      for (int i = 0; i < n; ++i) {
        if (rep % 3 == 0) {  // clusters, with exact duplicates
          const int c = r.below(5);
          x[i] = 0.2 + 0.3 * c + 0.01 * r.unif();
          y[i] = 0.5 + 0.01 * r.unif();
          if (i % 7 == 0 && i > 0) {
            x[i] = x[i - 1];
            y[i] = y[i - 1];
          }
        } else if (rep % 3 == 1) {  // collinear
          x[i] = 2 * r.unif();
          y[i] = 0.5 * x[i] / 2;
        } else {
          x[i] = 2 * r.unif();
          y[i] = r.unif();
        }
      }
      clab_nnidx(x, y, n, win, K, 0, 1, dist, iw, wk, out);
      brute(x, y, n, K, bd);
      bool eq = true;
      for (int q = 0; q < n * K; ++q) eq = eq && dist[q] == bd[q];
      tk_true(T, "nn.brute", eq);
    }
  }
  // invariance: translation, a quarter turn with the window, permutation
  {
    double x[150], y[150], u[150], v[150], o2[40];
    patp(x, y, 150);
    clab_nnidx(x, y, 150, win, 3, 0, 1, nullptr, iw, wk, out);
    for (int i = 0; i < 150; ++i) {
      u[149 - i] = -y[i] + 5e5;
      v[149 - i] = x[i] - 3e5;
    }
    const double w2[4] = {-1 + 5e5, 0 + 5e5, 0 - 3e5, 2 - 3e5};
    clab_nnidx(u, v, 150, w2, 3, 0, 1, nullptr, iw, wk, o2);
    for (int k = 0; k < 3; ++k) tk_eq(T, "nn.inv", o2[6 * k], out[6 * k], 1e-9);
    tk_eq(T, "nn.invce", o2[7], out[7], 1e-9);
  }
  // Monte Carlo test under CSR: p-values have the nominal size, and the
  // Donnelly-corrected ratio averages 1 in the rectangle
  {
    cgeo::rng r(31);
    double x[100], y[100];
    int rej = 0, rej2 = 0;
    double sr = 0;
    // default: 60 data sets and 19 simulations each; extended: 400 and 99
    const int reps = big ? 400 : 60, ns = big ? 99 : 19;
    for (int rep = 0; rep < reps; ++rep) {
      for (int i = 0; i < 100; ++i) {
        x[i] = 2 * r.unif();
        y[i] = r.unif();
      }
      clab_nnidx(x, y, 100, win, 2, ns, 1000 + rep, nullptr, iw, wk, out);
      if (out[3] <= 0.05) ++rej;
      if (out[6 + 4] <= 0.05) ++rej2;
      sr += out[13];
    }
    const double sd = std::sqrt(0.05 * 0.95 / reps);
    tk_true(T, "nn.size1", std::fabs(rej / double(reps) - 0.05) < 4 * sd);
    tk_true(T, "nn.size2", std::fabs(rej2 / double(reps) - 0.05) < 4 * sd);
    info(big, "  CSR: mean Donnelly ratio %.4f, size %.3f %.3f\n", sr / reps,
                rej / double(reps), rej2 / double(reps));
    tk_true(T, "nn.donnelly", std::fabs(sr / reps - 1) < (big ? 0.01 : 0.03));
    // a clustered pattern is detected
    for (int i = 0; i < 100; ++i) {
      const int c = i % 4;
      x[i] = 0.3 + 0.4 * c + 0.05 * r.unif();
      y[i] = 0.5 + 0.05 * r.unif();
    }
    clab_nnidx(x, y, 100, win, 1, 99, 5, nullptr, iw, wk, out);
    tk_true(T, "nn.cluster", out[3] == 0.01 && out[2] < 0.5);
  }
  // guards
  double x[3] = {0.5, 1, 1.5}, y[3] = {0.5, 0.5, 0.5};
  tk_true(T, "nn.n", clab_nnidx(x, y, 1, win, 1, 0, 1, nullptr, iw, wk, out) ==
                 CLAB_EARG);
  tk_true(T, "nn.k", clab_nnidx(x, y, 3, win, 3, 0, 1, nullptr, iw, wk, out) ==
                 CLAB_EARG);
  const double wbad[4] = {0, 1, 0, 1};
  tk_true(T, "nn.outside", clab_nnidx(x, y, 3, wbad, 1, 0, 1, nullptr, iw, wk, out) ==
                       CLAB_EDOM);
  const double whuge[4] = {0, 1e151, 0, 1};
  tk_true(T, "nn.huge", clab_nnidx(x, y, 3, whuge, 1, 0, 1, nullptr, iw, wk,
                                   out) == CLAB_EDOM);
  const double wflat[4] = {0, 2, 0.5, 0.5};
  tk_true(T, "nn.flat", clab_nnidx(x, y, 3, wflat, 1, 0, 1, nullptr, iw, wk, out) ==
                    CLAB_EDOM);
  tk_true(T, "nn.n2", clab_nnidx(x, y, 2, win, 1, 0, 1, nullptr, iw, wk, out) ==
                  CLAB_OK &&
                  out[0] == 0.5);
}

#ifdef CLAB_MAIN
double now() {
  return std::chrono::duration<double>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}

// a heap buffer for the timings (kept out of the program's static size)
template <class T>
struct hb {
  T *p;
  explicit hb(long n) : p(static_cast<T *>(std::malloc(sizeof(T) * n))) {
    if (p == nullptr) std::abort();
  }
  ~hb() { std::free(p); }
  hb(const hb &) = delete;
  hb &operator=(const hb &) = delete;
};

double med5(double *v) { return tk_q7(v, 5, 0.5); }

void fambench() {
  hb<double> x_(1000000), y_(1000000), wk_(2000100);
  double *x = x_.p, *y = y_.p, *wk = wk_.p;
  hb<int> iw_(2000100);
  int *iw = iw_.p;
  double out[64];
  cgeo::rng r(4);
  for (int i = 0; i < 1000000; ++i) {
    x[i] = 1e4 * r.unif();
    y[i] = 1e4 * r.unif();
  }
  const double win[4] = {0, 1e4, 0, 1e4};
  const int ns[3] = {1000, 100000, 1000000};
  for (int n : ns) {
    double t1[5], t5[5];
    for (int k = 0; k < 5; ++k) {
      double a = now();
      clab_nnidx(x, y, n, win, 1, 0, 1, nullptr, iw, wk, out);
      t1[k] = now() - a;
      a = now();
      clab_nnidx(x, y, n, win, 5, 0, 1, nullptr, iw, wk, out);
      t5[k] = now() - a;
    }
    std::printf("bench nnd n = %d: nnidx k = 1 %.5f s, k = 5 %.5f s "
                "(medians of 5)\n",
                n, med5(t1), med5(t5));
  }
  double a = now();
  clab_nnidx(x, y, 1000, win, 1, 999, 1, nullptr, iw, wk, out);
  std::printf("bench nnd n = 1000, 999 simulations: %.4f s\n", now() - a);
}

#endif

}  // namespace

void t_nnd(tk &T, int big) { famtest(T, big); }
#ifdef CLAB_MAIN
void b_nnd() { fambench(); }
#endif
