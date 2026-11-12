// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests of the hot-spot kernels (src/hsp.cpp).
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "cgeo.h"
#include "crimlab.h"
#include "tchk.h"

namespace {

// Reference values from dev/doss/hsp.R (base R).
// R 4.5.0 
// kmeans(X, centers = X[1:3, ], algorithm = "Lloyd", iter.max = 1000)
const double km3cl[] = {
    0, 1, 2, 1, 1, 2, 1, 1, 2, 1, 0, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 1, 0, 2,
    2, 1, 0, 2, 1, 2, 2, 0, 1, 1, 1, 2, 2, 2, 0, 2, 2, 2, 1, 1, 1, 2, 2, 2,
    0, 2, 2, 0, 2, 1, 2, 1, 0, 0, 1, 1, 1, 2, 0, 0, 1, 0, 2, 2, 1, 0, 0, 1,
    0, 0, 1, 2, 2, 0, 1, 0, 1, 2, 1, 1, 1, 0, 0, 1, 2, 1, 2, 0, 2, 2, 0, 2,
    2, 0, 1, 1, 1, 0, 2, 2, 0, 2, 2, 0, 1, 1, 2, 0, 0, 0, 1, 0, 0, 1, 0, 2,
    0, 2, 1, 1, 0, 2, 2, 0, 2, 2, 2, 1, 1, 1, 2, 2, 2, 0, 2, 2, 0, 1, 1, 2,
    0, 2, 0, 1, 0, 0, 1, 0, 2, 0, 2, 1, 2, 0, 2, 2, 0, 2, 2, 0, 2, 1, 1, 0,
    2, 0, 1, 2, 2, 0, 2, 2, 2, 1, 0, 1, 2, 1, 1, 2, 1, 1, 1, 0, 2, 0, 1, 0,
    0, 1, 0, 0, 1, 0, 2, 0};
const double km3cen[] = {
    0.32369689741479235, 1.64286130705999, 0.92846572627177348,
    0.44210015560426325, 0.45027872945821973, 0.55131239097215645};
const double km3wss[] = {
    23.916916620598911, 6.4538614805891656, 8.7410247005878752,
    8.7220304394218715};
// kmeans(X, centers = X[1:6, ], algorithm = "Lloyd", iter.max = 1000)
const double km6cl[] = {
    0, 1, 2, 3, 4, 5, 3, 3, 2, 1, 5, 1, 2, 4, 3, 5, 3, 4, 5, 4, 1, 3, 5, 2,
    5, 1, 5, 5, 1, 5, 2, 0, 3, 1, 4, 5, 4, 2, 0, 3, 4, 2, 1, 1, 3, 5, 3, 2,
    0, 2, 2, 5, 2, 3, 4, 1, 5, 0, 3, 1, 1, 3, 0, 0, 1, 0, 2, 5, 1, 0, 5, 4,
    0, 0, 1, 2, 2, 0, 3, 0, 1, 3, 1, 1, 4, 0, 0, 1, 2, 4, 4, 5, 2, 2, 0, 2,
    2, 0, 4, 4, 3, 0, 2, 2, 5, 2, 2, 5, 4, 3, 4, 0, 0, 5, 4, 0, 0, 4, 5, 5,
    5, 3, 4, 3, 5, 2, 4, 0, 2, 3, 2, 3, 1, 4, 5, 2, 3, 5, 2, 2, 5, 4, 4, 3,
    0, 5, 5, 1, 5, 0, 1, 0, 5, 0, 2, 3, 4, 5, 2, 2, 5, 2, 5, 0, 3, 3, 4, 0,
    2, 0, 1, 5, 5, 5, 2, 2, 5, 1, 0, 1, 4, 4, 3, 2, 1, 1, 3, 0, 2, 5, 3, 0,
    5, 4, 0, 5, 4, 5, 5, 5};
const double km6cen[] = {
    0.3077012767445928, 1.859193187240759, 0.93531897136299613,
    1.4044632787669595, 1.39509927890366, 0.49305091386820926,
    0.23962023111317574, 0.5081349859237323, 0.45365511508287004,
    0.7155938319093168, 0.15008339857711814, 0.76135815790077821};
const double km6wss[] = {
    13.004052864494165, 1.9044089591344444, 2.0118032171982727,
    2.568693094545563, 1.6602181773413334, 1.8035677534063086,
    3.0553616628682421};

const double PI = cgeo::pi;
void patp(double *x, double *y, int n) {
  for (long long i = 0; i < n; ++i) {
    x[i] = static_cast<double>((i * i * 7919 + i * 31) % 1009) / 1009 * 2;
    y[i] = static_cast<double>((i * i * 104729 + i * 7) % 1013) / 1013;
  }
}
// blobs of given sizes about given centres, radius rad, then background
int blobs(cgeo::rng &r, const double *cx, const double *cy, const int *sz,
          int nb, double rad, int nbg, double *x, double *y) {
  int n = 0;
  for (int b = 0; b < nb; ++b)
    for (int t = 0; t < sz[b]; ++t) {
      const double a = 2 * PI * r.unif(), d = rad * std::sqrt(r.unif());
      x[n] = cx[b] + d * std::cos(a);
      y[n++] = cy[b] + d * std::sin(a);
    }
  for (int t = 0; t < nbg; ++t) {
    x[n] = r.unif();
    y[n++] = r.unif();
  }
  return n;
}

void t_fzmode(tk &T, int /*big*/) {
  cgeo::rng r(3);
  static double x[600], y[600];
  static int iw[2000], c[600];
  for (int rep = 0; rep < 6; ++rep) {
    const int n = 100 * (rep + 1);
    for (int i = 0; i < n; ++i) {
      x[i] = std::floor(r.unif() * 30) / 30;  // coincident points occur
      y[i] = std::floor(r.unif() * 30) / 30;
    }
    const double rad = rep == 0 ? 0.0 : 0.02 * rep;
    tk_true(T, "fz.ok", clab_fzmode(x, y, n, rad, iw, c) == CLAB_OK);
    bool eq = true;
    for (int i = 0; i < n; ++i) {
      int b = 0;
      for (int j = 0; j < n; ++j) {
        const double dx = x[j] - x[i], dy = y[j] - y[i];
        if (dx * dx + dy * dy <= rad * rad) ++b;
      }
      eq = eq && b == c[i];
    }
    tk_true(T, "fz.brute", eq);
  }
  double one = 0.5, mr = -1, nanv = std::nan("");
  tk_true(T, "fz.nan", clab_fzmode(&nanv, &one, 1, 1, iw, c) == CLAB_EDOM);
  tk_true(T, "fz.one", clab_fzmode(&one, &one, 1, 1, iw, c) == CLAB_OK && c[0] == 1);
  tk_true(T, "fz.neg", clab_fzmode(&one, &one, 1, mr, iw, c) == CLAB_EDOM);
}

void t_kmpart(tk &T, int /*big*/) {
  static double x[200], y[200], cen[40], work[400], out[100];
  static int cl[200], iw[400];
  patp(x, y, 200);
  // stats::kmeans(..., algorithm = "Lloyd") from the same starting centres
  const int ks[2] = {3, 6};
  const double *rcl[2] = {km3cl, km6cl}, *rcen[2] = {km3cen, km6cen},
               *rwss[2] = {km3wss, km6wss};
  for (int q = 0; q < 2; ++q) {
    const int k = ks[q];
    for (int j = 0; j < k; ++j) {
      cen[j] = x[j];
      cen[k + j] = y[j];
    }
    tk_true(T, "km.ok", clab_kmpart(x, y, 200, k, 0, nullptr, 0, 0, 1000, cen, cl, iw,
                            work, out) == CLAB_OK);
    tk_true(T, "km.conv", out[3] == 1);
    bool same = true;
    for (int i = 0; i < 200; ++i) same = same && cl[i] == rcl[q][i];
    tk_true(T, "km.cluster", same);
    for (int j = 0; j < 2 * k; ++j) tk_eq(T, "km.centre", cen[j], rcen[q][j], 1e-13);
    tk_eq(T, "km.wss", out[0], rwss[q][0], 1e-12);
    for (int j = 0; j < k; ++j) tk_eq(T, "km.wssj", out[5 + 3 * j], rwss[q][1 + j], 1e-12);
  }
  // one group: the mean, and the total sum of squares
  {
    clab_kmpart(x, y, 200, 1, 0, nullptr, 0, 0, 100, cen, cl, iw, work, out);
    double mx = 0, my = 0, ss = 0;
    for (int i = 0; i < 200; ++i) {
      mx += x[i];
      my += y[i];
    }
    mx /= 200;
    my /= 200;
    for (int i = 0; i < 200; ++i)
      ss += (x[i] - mx) * (x[i] - mx) + (y[i] - my) * (y[i] - my);
    tk_eq(T, "km.k1x", cen[0], mx, 1e-14);
    tk_eq(T, "km.k1ss", out[0], ss, 1e-12);
  }
  // Lloyd's steps never raise the sum of squares; Weber steps never raise
  // the sum of distances
  for (int upd = 0; upd < 2; ++upd) {
    for (int j = 0; j < 5; ++j) {
      cen[j] = x[10 * j];
      cen[5 + j] = y[10 * j];
    }
    double prev = 1e300;
    bool mono = true;
    for (int s = 0; s < 40; ++s) {
      clab_kmpart(x, y, 200, 5, 0, nullptr, 0, upd, 1, cen, cl, iw, work, out);
      const double v = out[upd];
      if (v > prev * (1 + 1e-14)) mono = false;
      prev = v;
    }
    tk_true(T, upd ? "km.monoweber" : "km.monolloyd", mono);
  }
  // k-means++ seeding: the D^2 law on three points (0, 0), (1, 0), (3, 0):
  // the unordered seed pairs {0,1}, {0,2}, {1,2} have probabilities 0.1,
  // 0.3 + 0.3 / 1.3 and 0.8 / 3 + 0.4 / 3.9 (summing to 1)
  {
    const double px[3] = {0, 1, 3}, py[3] = {0, 0, 0};
    int cnt[3] = {0, 0, 0};
    const int reps = 30000;
    for (int s = 0; s < reps; ++s) {
      double c2[4];
      int c3[3], i3[6];
      double w3[6], o3[10];
      clab_kmpart(px, py, 3, 2, 1, nullptr, s, 0, 0, c2, c3, i3, w3, o3);
      const double a = c2[0], b = c2[1];
      const double lo = a < b ? a : b, hi = a < b ? b : a;
      if (lo == 0 && hi == 1) ++cnt[0];
      if (lo == 0 && hi == 3) ++cnt[1];
      if (lo == 1 && hi == 3) ++cnt[2];
    }
    const double p[3] = {0.1, 0.3 + 0.3 / 1.3, 0.8 / 3 + 0.4 / 3.9};
    for (int t = 0; t < 3; ++t) {
      const double f = static_cast<double>(cnt[t]) / reps;
      tk_true(T, "km.pp", std::fabs(f - p[t]) < 4 * std::sqrt(p[t] * (1 - p[t]) / reps));
    }
    double c2[4];
    int c3[3], i3[6];
    double w3[6], o3[10];
    const double same[3] = {2, 2, 2};
    tk_true(T, "km.ppdup", clab_kmpart(same, same, 3, 2, 1, nullptr, 1, 0, 0, c2, c3,
                               i3, w3, o3) == CLAB_EDOM);
  }
  // separated cells: three blobs inside single cells of side 0.1
  {
    cgeo::rng r(9);
    const double bx[3] = {0.15, 0.55, 0.85}, by[3] = {0.15, 0.75, 0.25};
    const int sz[3] = {30, 20, 10};
    static double u[200], v[200];
    const int n = blobs(r, bx, by, sz, 3, 0.04, 15, u, v);
    // the bounding box starts at the smallest coordinates: put a point at
    // (0, 0) and (1, 1) so that cells align with tenths
    u[n] = 0;
    v[n] = 0;
    u[n + 1] = 1;
    v[n + 1] = 1;
    const double par[2] = {0.1, 0.3};
    tk_true(T, "km.sepok", clab_kmpart(u, v, n + 2, 3, 2, par, 0, 0, 0, cen, cl, iw,
                               work, out) == CLAB_OK);
    for (int j = 0; j < 3; ++j) {
      tk_eq(T, "km.sepx", cen[j], bx[j], 1e-12);
      tk_eq(T, "km.sepy", cen[3 + j], by[j], 1e-12);
    }
    const double par2[2] = {0.1, 5.0};
    tk_true(T, "km.sepfail", clab_kmpart(u, v, n + 2, 3, 2, par2, 0, 0, 0, cen, cl, iw,
                                 work, out) == CLAB_EDOM);
  }
  // translation and permutation of the input
  {
    static double u[200], v[200], c1[12], c2[12];
    static int l2[200];
    for (int i = 0; i < 200; ++i) {
      u[199 - i] = x[i] + 3e5;
      v[199 - i] = y[i] - 7e5;
    }
    for (int j = 0; j < 6; ++j) {
      c1[j] = x[j];
      c1[6 + j] = y[j];
      c2[j] = x[j] + 3e5;
      c2[6 + j] = y[j] - 7e5;
    }
    clab_kmpart(x, y, 200, 6, 0, nullptr, 0, 0, 1000, c1, cl, iw, work, out);
    double o2[40];
    clab_kmpart(u, v, 200, 6, 0, nullptr, 0, 0, 1000, c2, l2, iw, work, o2);
    bool same = true;
    for (int i = 0; i < 200; ++i) same = same && l2[199 - i] == cl[i];
    tk_true(T, "km.inv", same);
    tk_eq(T, "km.invwss", o2[0], out[0], 1e-8);
  }
  tk_true(T, "km.k", clab_kmpart(x, y, 200, 0, 0, nullptr, 0, 0, 10, cen, cl, iw, work,
                         out) == CLAB_EARG);
  tk_true(T, "km.kn", clab_kmpart(x, y, 3, 4, 1, nullptr, 0, 0, 10, cen, cl, iw, work,
                          out) == CLAB_EARG);
}

// the first order of NNH by brute force: degrees, the seed order, the
// takeover, the size filter and the numbering
void nnhbrute(const double *x, const double *y, int n, double thr, int minp,
              int *lab) {
  static int deg[2000], ord[2000], asg[2000], sz[2000], ok[2000], nid[2000];
  for (int i = 0; i < n; ++i) {
    deg[i] = 0;
    for (int j = 0; j < n; ++j)
      if (j != i && std::hypot(x[j] - x[i], y[j] - y[i]) < thr) ++deg[i];
    ord[i] = i;
    asg[i] = -1;
  }
  for (int a = 1; a < n; ++a)  // stable insertion sort by degree desc
    for (int b = a; b > 0 && deg[ord[b]] > deg[ord[b - 1]]; --b) {
      const int t = ord[b];
      ord[b] = ord[b - 1];
      ord[b - 1] = t;
    }
  int nf = 0;
  for (int t = 0; t < n && deg[ord[t]] > 0; ++t) {
    const int s = ord[t];
    if (asg[s] >= 0) continue;
    asg[s] = nf;
    for (int j = 0; j < n; ++j)
      if (asg[j] < 0 && std::hypot(x[j] - x[s], y[j] - y[s]) < thr) asg[j] = nf;
    ++nf;
  }
  for (int q = 0; q < nf; ++q) sz[q] = 0;
  for (int i = 0; i < n; ++i)
    if (asg[i] >= 0) ++sz[asg[i]];
  int kept = 0;
  for (int q = 0; q < nf; ++q) ok[q] = sz[q] >= minp;
  // number kept clusters by size desc, then seed order
  for (int pass = n; pass >= minp; --pass)
    for (int q = 0; q < nf; ++q)
      if (ok[q] && sz[q] == pass) nid[q] = kept++;
  for (int i = 0; i < n; ++i)
    lab[i] = (asg[i] >= 0 && ok[asg[i]]) ? nid[asg[i]] : -1;
}

void t_nnh(tk &T, int /*big*/) {
  static double x[2000], y[2000], work[8000], cs[22000];
  static int iw[12100], cl[6000], lab[2000];
  int ncl = 0;
  // brute force on random clustered sets, fixed threshold
  cgeo::rng r(41);
  for (int rep = 0; rep < 10; ++rep) {
    const double bx[4] = {0.2, 0.7, 0.4, 0.8}, by[4] = {0.3, 0.8, 0.6, 0.2};
    const int sz[4] = {15 + rep, 12, 9, 25 - rep};
    const int n = blobs(r, bx, by, sz, 4, 0.05, 60 + 10 * rep, x, y);
    const double thr = 0.02 + 0.005 * rep;
    tk_true(T, "nnh.ok", clab_nnhclu(x, y, n, 1, thr, 0, 5, 2, 1, 1, 1, iw, work, cl,
                             &ncl, cs) == CLAB_OK);
    nnhbrute(x, y, n, thr, 5, lab);
    bool eq = true;
    for (int i = 0; i < n; ++i) eq = eq && cl[i] == lab[i];
    tk_true(T, "nnh.brute", eq);
    // each row: the members' count and mean agree with the labels
    for (int q = 0; q < ncl; ++q) {
      int m = 0;
      double mx = 0;
      for (int i = 0; i < n; ++i)
        if (cl[i] == q) {
          ++m;
          mx += x[i];
        }
      tk_true(T, "nnh.count", cs[11 * q + 1] == m && cs[11 * q + 2] == m);
      tk_eq(T, "nnh.mean", cs[11 * q + 3], mx / m, 1e-13);
    }
  }
  // the random-expectation threshold, and a two-level hierarchy: 8 groups
  // of 10 points (radius 0.01) in two super-groups whose group centres are
  // 0.1 apart, the super-groups 0.6 apart. Order 1: n = 80 points, A = 1,
  // z = 0: threshold 0.5 sqrt(1 / 80) = 0.0559 links each group alone.
  // Order 2: 8 units, threshold 0.5 sqrt(1 / 8) = 0.177 links each
  // super-group. Order 3: 2 units, fewer than 4: stop.
  {
    double gx[8], gy[8];
    for (int g = 0; g < 8; ++g) {
      const int s = g / 4, w = g % 4;
      gx[g] = 0.15 + 0.6 * s + 0.1 * (w % 2);
      gy[g] = 0.3 + 0.1 * (w / 2);
    }
    int sz[8];
    for (int g = 0; g < 8; ++g) sz[g] = 10;
    const int n = blobs(r, gx, gy, sz, 8, 0.01, 0, x, y);
    tk_true(T, "nnh.h", clab_nnhclu(x, y, n, 1, 0, 0, 5, 2, 3, 1, 1, iw, work, cl,
                            &ncl, cs) == CLAB_OK);
    tk_true(T, "nnh.hn", ncl == 10);
    tk_eq(T, "nnh.thr1", cs[10], 0.5 * std::sqrt(1.0 / 80), 1e-15);
    tk_eq(T, "nnh.thr2", cs[11 * 8 + 10], 0.5 * std::sqrt(1.0 / 8), 1e-15);
    tk_true(T, "nnh.order", cs[11 * 7] == 1 && cs[11 * 8] == 2 && cs[11 * 9] == 2);
    tk_true(T, "nnh.units", cs[11 * 8 + 2] == 4 && cs[11 * 8 + 1] == 40);
    // every point of a group shares its first- and second-order clusters
    bool ok = true;
    for (int g = 0; g < 8; ++g)
      for (int t = 1; t < 10; ++t) {
        const int i = 10 * g + t;
        ok = ok && cl[i] == cl[10 * g] && cl[n + i] == cl[n + 10 * g];
      }
    // groups 0-3 share one second-order cluster, groups 4-7 the other
    for (int g = 0; g < 8; ++g)
      ok = ok && cl[n + 10 * g] == cl[n + 40 * (g / 4)];
    ok = ok && cl[n] != cl[n + 40] && cl[n] >= 8 && cl[n + 40] >= 8;
    tk_true(T, "nnh.tree", ok);
    tk_true(T, "nnh.none3", cl[2 * n] == -1);
    // z shifts the threshold by z sqrt((4 - pi) A / (4 pi)) / m
    clab_nnhclu(x, y, n, 1, 0, -1.645, 5, 2, 1, 1, 1, iw, work, cl, &ncl, cs);
    tk_eq(T, "nnh.z", cs[10],
        0.5 * std::sqrt(1.0 / 80) - 1.645 * std::sqrt((4 - PI) / (4 * PI)) / 80,
        1e-15);
  }
  // invariance: translation and permutation of the input
  {
    const double bx[3] = {0.2, 0.5, 0.8}, by[3] = {0.5, 0.2, 0.7};
    const int sz[3] = {12, 18, 9};
    const int n = blobs(r, bx, by, sz, 3, 0.04, 40, x, y);
    static double u[200], v[200];
    static int c2[200];
    for (int i = 0; i < n; ++i) {
      u[n - 1 - i] = x[i] + 1e5;
      v[n - 1 - i] = y[i] + 2e5;
    }
    clab_nnhclu(x, y, n, 1, 0.03, 0, 5, 2, 1, 1, 1, iw, work, cl, &ncl, cs);
    int nc2;
    static double cs2[2200];
    clab_nnhclu(u, v, n, 1, 0.03, 0, 5, 2, 1, 1, 1, iw, work, c2, &nc2, cs2);
    bool same = nc2 == ncl;
    for (int q = 0; q < ncl && same; ++q) same = cs[11 * q + 1] == cs2[11 * q + 1];
    tk_true(T, "nnh.inv", same);
  }
  double z2[2] = {0, 1};
  tk_true(T, "nnh.z", clab_nnhclu(z2, z2, 2, 1, 0, std::nan(""), 2, 2, 1, 1, 1,
                                  iw, work, cl, &ncl, cs) == CLAB_EDOM);
  // two members and divisor W - 2: the ellipse is not defined (NaN)
  tk_true(T, "nnh.ellnan", clab_nnhclu(z2, z2, 2, 1, 2, 0, 2, 2, 1, 2, 1, iw,
                                       work, cl, &ncl, cs) == CLAB_OK &&
                               std::isnan(cs[7]) && cs[1] == 2);
  tk_true(T, "nnh.minpts", clab_nnhclu(z2, z2, 2, 1, 0, 0, 1, 2, 1, 1, 1, iw, work, cl,
                               &ncl, cs) == CLAB_EARG);
  tk_true(T, "nnh.area", clab_nnhclu(z2, z2, 2, 0, 0, 0, 2, 2, 1, 1, 1, iw, work, cl,
                             &ncl, cs) == CLAB_EDOM);
  tk_true(T, "nnh.two", clab_nnhclu(z2, z2, 2, 1, 2, 0, 2, 2, 1, 0, 1, iw, work, cl,
                            &ncl, cs) == CLAB_OK && ncl == 1 && cl[0] == 0);
}

// the search-circle procedure by brute force
int stbrute(const double *x, const double *y, int n, const double *w, int nx,
            int ny, double r, int top, int minp, int *lab) {
  static int cnt[10000], ord[10000], par[10000], first[2000], sz[10000],
      nid[10000];
  const int nn = nx * ny;
  for (int q = 0; q < nn; ++q) {
    const double px = w[0] + (w[1] - w[0]) * (q % nx) / (nx - 1.0);
    const double py = w[2] + (w[3] - w[2]) * (q / nx) / (ny - 1.0);
    cnt[q] = 0;
    for (int i = 0; i < n; ++i) {
      const double dx = x[i] - px, dy = y[i] - py;
      if (dx * dx + dy * dy <= r * r) ++cnt[q];
    }
    ord[q] = q;
  }
  for (int a = 1; a < nn; ++a)
    for (int b = a; b > 0 && cnt[ord[b]] > cnt[ord[b - 1]]; --b) {
      const int t = ord[b];
      ord[b] = ord[b - 1];
      ord[b - 1] = t;
    }
  int nk = 0;
  while (nk < top && nk < nn && cnt[ord[nk]] >= 2) ++nk;
  for (int t = 0; t < nk; ++t) par[t] = t;
  auto root = [](int a) {
    while (par[a] != a) a = par[a];
    return a;
  };
  for (int i = 0; i < n; ++i) first[i] = -1;
  for (int t = 0; t < nk; ++t) {
    const int q = ord[t];
    const double px = w[0] + (w[1] - w[0]) * (q % nx) / (nx - 1.0);
    const double py = w[2] + (w[3] - w[2]) * (q / nx) / (ny - 1.0);
    for (int i = 0; i < n; ++i) {
      const double dx = x[i] - px, dy = y[i] - py;
      if (dx * dx + dy * dy > r * r) continue;
      if (first[i] < 0) {
        first[i] = t;
      } else {
        const int a = root(first[i]), b = root(t);
        if (a < b) par[b] = a;
        if (b < a) par[a] = b;
      }
    }
  }
  for (int t = 0; t < nk; ++t) sz[t] = 0;
  for (int i = 0; i < n; ++i)
    if (first[i] >= 0) ++sz[root(first[i])];
  int kept = 0;
  for (int t = 0; t < nk; ++t) nid[t] = -1;
  for (int s = n; s >= minp; --s)
    for (int t = 0; t < nk; ++t)
      if (root(t) == t && sz[t] == s) nid[t] = kept++;
  for (int i = 0; i < n; ++i) lab[i] = first[i] >= 0 ? nid[root(first[i])] : -1;
  return kept;
}

void t_stac(tk &T, int /*big*/) {
  static double x[2000], y[2000], work[20000], cs[16000];
  static int iw[30000], cl[2000], lab[2000];
  int ncl = 0;
  cgeo::rng r(55);
  const double win[4] = {0, 1, 0, 1};
  for (int rep = 0; rep < 8; ++rep) {
    const double bx[3] = {0.25, 0.7, 0.5}, by[3] = {0.3, 0.75, 0.55};
    const int sz[3] = {30, 25 + rep, 12};
    const int n = blobs(r, bx, by, sz, 3, 0.06, 80, x, y);
    const double rad = 0.04 + 0.01 * rep;
    const int top = 10 + 5 * rep;
    tk_true(T, "st.ok", clab_stacel(x, y, n, win, 20, 20, rad, top, 5, 1, 1, iw, work,
                            cl, &ncl, cs) == CLAB_OK);
    const int kb = stbrute(x, y, n, win, 20, 20, rad, top, 5, lab);
    bool eq = kb == ncl;
    for (int i = 0; i < n && eq; ++i) eq = cl[i] == lab[i];
    tk_true(T, "st.brute", eq);
    for (int q = 0; q < ncl; ++q) {
      int m = 0;
      for (int i = 0; i < n; ++i) m += cl[i] == q;
      tk_true(T, "st.count", cs[8 * q] == m);
    }
  }
  // two far blobs give two hot clusters holding them
  {
    const double bx[2] = {0.2, 0.8}, by[2] = {0.2, 0.8};
    const int sz[2] = {40, 40};
    const int n = blobs(r, bx, by, sz, 2, 0.03, 0, x, y);
    clab_stacel(x, y, n, win, 20, 20, 0.08, 25, 10, 1, 1, iw, work, cl, &ncl,
                cs);
    bool ok = ncl == 2;
    for (int i = 0; i < 80 && ok; ++i) ok = cl[i] == cl[(i / 40) * 40] && cl[i] >= 0;
    tk_true(T, "st.two", ok && cl[0] != cl[40]);
    tk_eq(T, "st.mean", cs[1] + cs[9], 1.0, 0.02);
  }
  tk_true(T, "st.grid", clab_stacel(x, y, 10, win, 1, 20, 0.1, 5, 2, 1, 1, iw, work, cl,
                            &ncl, cs) == CLAB_EARG);
  const double wrev[4] = {1, 0, 0, 1};
  tk_true(T, "st.win", clab_stacel(x, y, 10, wrev, 20, 20, 0.1, 5, 2, 1, 1, iw,
                                   work, cl, &ncl, cs) == CLAB_EDOM);
  tk_true(T, "st.r", clab_stacel(x, y, 10, win, 20, 20, 0, 5, 2, 1, 1, iw, work, cl,
                         &ncl, cs) == CLAB_EDOM);
}

void famtest(tk &T, int big) {
  t_fzmode(T, big);
  t_kmpart(T, big);
  t_nnh(T, big);
  t_stac(T, big);
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
  hb<double> x_(1000000), y_(1000000), work_(4000000), cs_(11000000);

  double *x = x_.p, *y = y_.p, *work = work_.p, *cs = cs_.p;
  hb<int> iw_(6000100), cl_(3000000), cnt_(1000000);

  int *iw = iw_.p, *cl = cl_.p, *cnt = cnt_.p;
  double cen[20], out[64];
  cgeo::rng r(6);
  // uniform points: with the radii scaled as 1 / sqrt(n) every search
  // meets about the same number of points, so the times show the
  // algorithms' growth; clustered data cost in proportion to the pairs
  // within the radius
  for (int i = 0; i < 1000000; ++i) {
    x[i] = r.unif();
    y[i] = r.unif();
  }
  const int ns[3] = {1000, 100000, 1000000};
  for (int n : ns) {
    double t[4][5];
    int ncl;
    const double win[4] = {0, 1, 0, 1};
    for (int s = 0; s < 5; ++s) {
      double a = now();
      clab_fzmode(x, y, n, 0.1 / std::sqrt(n / 1000.0), iw, cnt);
      t[0][s] = now() - a;
      a = now();
      clab_kmpart(x, y, n, 10, 1, nullptr, 1, 0, 100, cen, cl, iw, work, out);
      t[1][s] = now() - a;
      a = now();
      clab_nnhclu(x, y, n, 1, 0, 0, 10, 2, 3, 1, 1, iw, work, cl, &ncl, cs);
      t[2][s] = now() - a;
      a = now();
      clab_stacel(x, y, n, win, 20, 20, 0.05 * std::sqrt(1000.0 / n), 25, 10,
                  1, 1, iw, work, cl, &ncl, cs);
      t[3][s] = now() - a;
    }
    std::printf("bench hsp n = %d: fzmode %.5f kmpart(k = 10, <= 100 it) %.5f "
                "nnhclu %.5f stacel %.5f s (medians of 5)\n",
                n, med5(t[0]), med5(t[1]), med5(t[2]),
                med5(t[3]));
  }
}

#endif

}  // namespace

void t_hsp(tk &T, int big) { famtest(T, big); }
#ifdef CLAB_MAIN
void b_hsp() { fambench(); }
#endif
