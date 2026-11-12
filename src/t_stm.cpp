// SPDX-License-Identifier: AGPL-3.0-or-later
//
// Standalone tests of the space-time kernels (src/stm.cpp).
#include <chrono>
#include <cmath>
#include <cstdio>

#include "cgeo.h"
#include "crimlab.h"
#include "tchk.h"

namespace {

// a message of the extended run (standalone, "ext"); silent by default
template <class... A>
void info(int big, const char *f, A... a) {
  if (big) std::printf(f, a...);
}

// Reference values from dev/doss/stm.R; the first line names the versions.
// R 4.5.0 vegan 2.6.10 
// mantel(dist(cbind(x, y)), dist(t), permutations = 0)$statistic
const double mr[] = {
    0.044297537295815606};
// the same with 1 / (d + 0.1) and 1 / (dt + 2)
const double mrr[] = {
    -0.0046212289154358962};

// Pearson correlation of the pair dissimilarities, written out directly
double rdirect(const double *x, const double *y, const double *t, int n,
               double cs, double ct) {
  double a[2000], b[2000];
  int m = 0;
  for (int i = 0; i < n; ++i)
    for (int j = i + 1; j < n; ++j) {
      const double d = std::hypot(x[i] - x[j], y[i] - y[j]);
      const double e = std::fabs(t[i] - t[j]);
      a[m] = cs > 0 ? 1 / (d + cs) : d;
      b[m++] = ct > 0 ? 1 / (e + ct) : e;
    }
  double ma = 0, mb = 0;
  for (int k = 0; k < m; ++k) {
    ma += a[k];
    mb += b[k];
  }
  ma /= m;
  mb /= m;
  double sab = 0, saa = 0, sbb = 0;
  for (int k = 0; k < m; ++k) {
    sab += (a[k] - ma) * (b[k] - mb);
    saa += (a[k] - ma) * (a[k] - ma);
    sbb += (b[k] - mb) * (b[k] - mb);
  }
  return sab / std::sqrt(saa * sbb);
}

void t_mantel(tk &T, int big) {
  static int iw[4000];
  double out[5];
  double x[60], y[60], t[60];
  for (long long i = 0; i < 60; ++i) {
    x[i] = static_cast<double>((i * i * 7919 + i * 31) % 1009) / 1009 * 2;
    y[i] = static_cast<double>((i * i * 104729 + i * 7) % 1013) / 1013;
    t[i] = static_cast<double>((i * 17) % 61) + x[i];
  }
  // vegan::mantel parity, raw and reciprocal dissimilarities
  tk_true(T, "man.ok", clab_mantst(x, y, t, 60, 0, 0, 0, 1, iw, out) == CLAB_OK);
  tk_eq(T, "man.vegan", out[0], mr[0], 1e-13);
  tk_eq(T, "man.pairs", out[1], 1770, 0);
  tk_true(T, "man.nop", std::isnan(out[2]));
  clab_mantst(x, y, t, 60, 0.1, 2, 0, 1, iw, out);
  tk_eq(T, "man.veganrec", out[0], mrr[0], 1e-13);
  // the direct formula on small random sets
  cgeo::rng r(12);
  for (int rep = 0; rep < 10; ++rep) {
    const int n = 3 + 5 * rep;
    double u[50], v[50], w[50];
    for (int i = 0; i < n; ++i) {
      u[i] = r.unif();
      v[i] = r.unif();
      w[i] = 100 * r.unif();
    }
    const double cs = rep % 2 ? 0.05 : 0, ct = rep % 3 ? 1.0 : 0;
    clab_mantst(u, v, w, n, cs, ct, 0, 1, iw, out);
    tk_eq(T, "man.direct", out[0], rdirect(u, v, w, n, cs, ct), 1e-12);
    // invariance: shift and scale times, shift coordinates, reverse order
    if (cs == 0 && ct == 0) {
      double u2[50], v2[50], w2[50], o2[5];
      for (int i = 0; i < n; ++i) {
        u2[n - 1 - i] = u[i] + 5e5;
        v2[n - 1 - i] = v[i] - 5e5;
        w2[n - 1 - i] = 3 * w[i] + 1e6;
      }
      clab_mantst(u2, v2, w2, n, 0, 0, 0, 1, iw, o2);
      tk_eq(T, "man.inv", o2[0], out[0], 1e-9);
    }
  }
  // the permutation test has its nominal size when times are independent
  // of places (300 data sets of 30 events, 99 permutations), and power
  // against space-time clustering
  {
    int rej = 0, rej2 = 0;
    const int reps = big ? 300 : 80;
    for (int rep = 0; rep < reps; ++rep) {
      double u[30], v[30], w[30];
      for (int i = 0; i < 30; ++i) {
        u[i] = r.unif();
        v[i] = r.unif();
        w[i] = 365 * r.unif();
      }
      clab_mantst(u, v, w, 30, 0, 0, 99, 100 + rep, iw, out);
      if (out[2] <= 0.05) ++rej;
      if (out[4] <= 0.05) ++rej2;
    }
    const double sd = std::sqrt(0.05 * 0.95 / reps);
    info(big, "  Mantel permutation test size: one-sided %.3f, two-sided %.3f\n",
                rej / double(reps), rej2 / double(reps));
    tk_true(T, "man.size", std::fabs(rej / double(reps) - 0.05) < 4 * sd);
    tk_true(T, "man.size2", std::fabs(rej2 / double(reps) - 0.05) < 4 * sd);
    double u[40], v[40], w[40];
    for (int i = 0; i < 40; ++i) {
      w[i] = 365 * r.unif();
      u[i] = w[i] / 365 + 0.05 * r.unif();  // place follows time
      v[i] = r.unif() * 0.1;
    }
    clab_mantst(u, v, w, 40, 0, 0, 199, 9, iw, out);
    tk_true(T, "man.power", out[2] == 1.0 / 200 && out[0] > 0.5);
  }
  // guards
  double c[3] = {1, 1, 1}, a[3] = {0, 1, 2};
  tk_true(T, "man.const", clab_mantst(a, a, c, 3, 0, 0, 0, 1, iw, out) == CLAB_EDOM);
  double nv[3] = {0, std::nan(""), 2};
  tk_true(T, "man.nan", clab_mantst(a, a, nv, 3, 0, 0, 0, 1, iw, out) == CLAB_EDOM);
  tk_true(T, "man.n", clab_mantst(a, a, a, 2, 0, 0, 0, 1, iw, out) == CLAB_EARG);
  tk_true(T, "man.c", clab_mantst(a, a, a, 3, -1, 0, 0, 1, iw, out) == CLAB_EDOM);
}

void t_stma(tk &T, int /*big*/) {
  // points (k, 2k) at times k, given in a scrambled order
  double x[9], y[9], t[9], out[36];
  int iw[9];
  const int perm[9] = {4, 0, 8, 2, 6, 1, 7, 3, 5};
  for (int i = 0; i < 9; ++i) {
    const int k = perm[i];
    x[i] = k;
    y[i] = 2 * k;
    t[i] = 10 + k;
  }
  tk_true(T, "sma.ok", clab_stmavg(x, y, t, 9, 5, iw, out) == CLAB_OK);
  // ends: windows of 3, 4; the rest: 5, centred on k
  const double ex[9] = {1, 1.5, 2, 3, 4, 5, 6, 6.5, 7};
  const int em[9] = {3, 4, 5, 5, 5, 5, 5, 4, 3};
  for (int q = 0; q < 9; ++q) {
    tk_true(T, "sma.idx", perm[static_cast<int>(out[4 * q])] == q);
    tk_eq(T, "sma.x", out[4 * q + 1], ex[q], 1e-15);
    tk_eq(T, "sma.y", out[4 * q + 2], 2 * ex[q], 1e-15);
    tk_true(T, "sma.m", out[4 * q + 3] == em[q]);
  }
  clab_stmavg(x, y, t, 9, 1, iw, out);
  bool same = true;
  for (int q = 0; q < 9; ++q) same = same && out[4 * q + 1] == q;
  tk_true(T, "sma.span1", same);
  // tied times keep input order
  double tt[3] = {1, 1, 0}, xx[3] = {5, 6, 7}, o3[12];
  int i3[3];
  clab_stmavg(xx, xx, tt, 3, 1, i3, o3);
  tk_true(T, "sma.ties", o3[0] == 2 && o3[4] == 0 && o3[8] == 1);
  tk_true(T, "sma.even", clab_stmavg(x, y, t, 9, 4, iw, out) == CLAB_EARG);
  double nanv = std::nan("");
  tk_true(T, "sma.nan", clab_stmavg(&nanv, y, t, 1, 1, iw, out) == CLAB_EDOM);
}

void famtest(tk &T, int big) {
  t_mantel(T, big);
  t_stma(T, big);
}

#ifdef CLAB_MAIN
double now() {
  return std::chrono::duration<double>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}

void fambench() {
  double x[2000], y[2000], t[2000];
  int iw[2000];
  double out[5];
  cgeo::rng r(2);
  for (int i = 0; i < 2000; ++i) {
    x[i] = r.unif();
    y[i] = r.unif();
    t[i] = r.unif();
  }
  const int ns[2] = {200, 1000};
  for (int n : ns) {
    const double a = now();
    clab_mantst(x, y, t, n, 0, 0, 99, 1, iw, out);
    std::printf("bench stm: mantst n = %d, 99 permutations: %.4f s\n", n,
                now() - a);
  }
}

#endif

}  // namespace

void t_stm(tk &T, int big) { famtest(T, big); }
#ifdef CLAB_MAIN
void b_stm() { fambench(); }
#endif
