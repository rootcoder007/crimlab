// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Weights on the simplex for synthetic controls (synthc()) and synthetic
// difference in differences (sdidx()). Each weight problem, min (1/m) ||A w
// + w0 - b||^2 + z2 ||w||^2 over w >= 0, sum w = 1 (w0 free when there is
// an intercept), is the nearest point to the origin of the convex hull of
// k points whose inner products are the Gram matrix below; Wolfe's (1976)
// minimum-norm-point algorithm finds it exactly in finitely many steps.
// Scratch is held in vectors.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "kutil.h"
#include "pev.h"

namespace pev {

namespace {

inline std::size_t z(int i) { return static_cast<std::size_t>(i); }

// Solves the bordered system [G_S 1; 1' 0] [a; mu] = [0; 1] for the corral
// s (ns points) by Gaussian elimination with partial pivoting: the affine
// minimiser a. false when the system is singular.
bool affine(const std::vector<double> &g, int k, const std::vector<int> &s, std::vector<double> &a) {
  const int n = static_cast<int>(s.size()) + 1;
  std::vector<double> m(z(n) * z(n + 1), 0.0);
  auto at = [&](int r, int c) -> double & { return m[z(r) * z(n + 1) + z(c)]; };
  for (int r = 0; r < n - 1; ++r) {
    for (int c = 0; c < n - 1; ++c) at(r, c) = g[z(s[z(r)]) * z(k) + z(s[z(c)])];
    at(r, n - 1) = 1.0;
    at(n - 1, r) = 1.0;
  }
  at(n - 1, n) = 1.0;
  double big = 0.0;
  for (double v : m) big = std::fmax(big, std::fabs(v));
  for (int c = 0; c < n; ++c) {
    int p = c;
    for (int r = c + 1; r < n; ++r) {
      if (std::fabs(at(r, c)) > std::fabs(at(p, c))) p = r;
    }
    if (!(std::fabs(at(p, c)) > 1e-14 * big)) return false;
    for (int q = 0; q <= n; ++q) std::swap(at(c, q), at(p, q));
    for (int r = 0; r < n; ++r) {
      if (r == c) continue;
      const double f = at(r, c) / at(c, c);
      for (int q = c; q <= n; ++q) at(r, q) -= f * at(c, q);
    }
  }
  a.resize(z(n - 1));
  for (int r = 0; r < n - 1; ++r) a[z(r)] = at(r, n) / at(r, r);
  return true;
}

}  // namespace

void gram(const double *a, const double *b, int m, int k, int icpt, double z2, double *g) {
  // the columns less b, centred over the m rows when there is an intercept
  std::vector<double> d(z(m) * z(k));
  for (int j = 0; j < k; ++j) {
    double mu = 0.0;
    for (int i = 0; i < m; ++i) d[z(j) * z(m) + z(i)] = a[z(j) * z(m) + z(i)] - b[i];
    if (icpt) {
      for (int i = 0; i < m; ++i) mu += d[z(j) * z(m) + z(i)] / m;
      for (int i = 0; i < m; ++i) d[z(j) * z(m) + z(i)] -= mu;
    }
  }
  for (int j = 0; j < k; ++j) {
    for (int l = j; l < k; ++l) {
      kutil::ksum s;
      for (int i = 0; i < m; ++i) s.add(d[z(j) * z(m) + z(i)] * d[z(l) * z(m) + z(i)]);
      g[z(j) * z(k) + z(l)] = g[z(l) * z(k) + z(j)] = s.val() / m + (j == l ? z2 : 0.0);
    }
  }
}

int mnp(const double *gp, int k, double *w, int maxit) {
  const std::vector<double> g(gp, gp + z(k) * z(k));
  double big = 0.0;
  int j0 = 0;
  for (int j = 0; j < k; ++j) {
    big = std::fmax(big, g[z(j) * z(k) + z(j)]);
    if (g[z(j) * z(k) + z(j)] < g[z(j0) * z(k) + z(j0)]) j0 = j;
  }
  const double eps = 1e-13 * (big > 0 ? big : 1.0);
  for (int j = 0; j < k; ++j) w[j] = 0.0;
  w[j0] = 1.0;
  std::vector<int> s = {j0};
  std::vector<double> gw(z(k)), a;
  for (int it = 0; it < maxit; ++it) {
    if (kutil::stop(static_cast<double>(k) * k)) return -1;
    // x'p_j = (G w)_j; x'x = w'G w
    for (int j = 0; j < k; ++j) {
      double v = 0.0;
      for (int i : s) v += w[i] * g[z(i) * z(k) + z(j)];
      gw[z(j)] = v;
    }
    double xx = 0.0;
    for (int i : s) xx += w[i] * gw[z(i)];
    int jn = 0;
    for (int j = 1; j < k; ++j) {
      if (gw[z(j)] < gw[z(jn)]) jn = j;
    }
    if (xx - gw[z(jn)] <= eps || std::find(s.begin(), s.end(), jn) != s.end()) return it;
    s.push_back(jn);
    // minor cycles: move to the affine minimiser of the corral, dropping
    // the points whose weight it would make negative
    for (;;) {
      if (!affine(g, k, s, a)) return -1;
      bool pos = true;
      for (double v : a) pos = pos && v > 1e-15;
      if (pos) {
        for (std::size_t q = 0; q < s.size(); ++q) w[s[q]] = a[q];
        break;
      }
      double th = 1.0;
      for (std::size_t q = 0; q < s.size(); ++q) {
        if (a[q] <= 1e-15) th = std::min(th, w[s[q]] / (w[s[q]] - a[q]));
      }
      std::vector<int> keep;
      for (std::size_t q = 0; q < s.size(); ++q) {
        const double v = (1.0 - th) * w[s[q]] + th * a[q];
        w[s[q]] = v > 1e-15 ? v : 0.0;
        if (w[s[q]] > 0.0) keep.push_back(s[q]);
      }
      s.swap(keep);
    }
  }
  return -1;
}

double sdtau(const double *yw, int nu, int nt, const std::vector<int> &ctl,
             const std::vector<int> &trd, int t0, bool fit, double *om, double *la) {
  const int n0 = static_cast<int>(ctl.size()), t1 = nt - t0;
  const int n1 = static_cast<int>(trd.size());
  // the treated mean path
  std::vector<double> y1(z(nt), 0.0);
  for (int s = 0; s < nt; ++s) {
    for (int i : trd) y1[z(s)] += yw[z(s) * z(nu) + z(i)] / n1;
  }
  if (fit) {
    // the noise level: the standard deviation of the comparison units'
    // first differences before adoption
    kutil::ksum s1, s2;
    double c = 0;
    for (int s = 0; s + 1 < t0; ++s) {
      for (int j : ctl) {
        const double d = yw[z(s + 1) * z(nu) + z(j)] - yw[z(s) * z(nu) + z(j)];
        s1.add(d);
        s2.add(d * d);
        c += 1;
      }
    }
    const double sg = c > 1 ? std::sqrt(std::max(0.0, (s2.val() - s1.val() * s1.val() / c) / (c - 1)))
                            : 0.0;
    const double zo = std::pow(static_cast<double>(n1) * t1, 0.25) * sg, zl = 1e-6 * sg;
    // unit weights: the comparison units' paths before adoption onto the
    // treated mean path, with an intercept and the ridge zo^2
    std::vector<double> a(z(t0) * z(n0)), b(y1.begin(), y1.begin() + t0), g(z(n0) * z(n0));
    for (int q = 0; q < n0; ++q) {
      for (int s = 0; s < t0; ++s) a[z(q) * z(t0) + z(s)] = yw[z(s) * z(nu) + z(ctl[z(q)])];
    }
    gram(a.data(), b.data(), t0, n0, 1, zo * zo, g.data());
    if (mnp(g.data(), n0, om, 10000) < 0) return kutil::nan();
    // time weights: the periods before adoption onto each comparison
    // unit's mean after it
    std::vector<double> at(z(n0) * z(t0)), bt(z(n0), 0.0), gt(z(t0) * z(t0));
    for (int q = 0; q < n0; ++q) {
      for (int s = 0; s < nt; ++s) {
        const double v = yw[z(s) * z(nu) + z(ctl[z(q)])];
        if (s < t0) at[z(s) * z(n0) + z(q)] = v;
        else bt[z(q)] += v / t1;
      }
    }
    gram(at.data(), bt.data(), n0, t0, 1, zl * zl, gt.data());
    if (mnp(gt.data(), t0, la, 10000) < 0) return kutil::nan();
  }
  // tau: the treated change from the weighted pre-period to the mean
  // post-period, less the weighted comparison units' change
  auto chg = [&](const double *path, int stride, int i) {
    double post = 0.0, pre = 0.0;
    for (int s = t0; s < nt; ++s) post += path[z(s) * z(stride) + z(i)] / t1;
    for (int s = 0; s < t0; ++s) pre += la[s] * path[z(s) * z(stride) + z(i)];
    return post - pre;
  };
  double tau = chg(y1.data(), 1, 0);
  for (int q = 0; q < n0; ++q) tau -= om[q] * chg(yw, nu, ctl[z(q)]);
  return tau;
}

}  // namespace pev
