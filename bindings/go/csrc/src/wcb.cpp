// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Inference for a difference in differences with few clusters: the
// restricted wild cluster bootstrap-t (Cameron, Gelbach and Miller 2008)
// of the first coefficient of a regression already within-transformed, and
// the Westfall-Young stepdown adjustment of several p-values from their
// joint null draws. Scratch is held in vectors.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <vector>

#include "kutil.h"
#include "pev.h"
#include "pevw.h"

namespace pev {

namespace {

inline std::size_t z(int i) { return static_cast<std::size_t>(i); }

}  // namespace

bool wsetup(wctx &c) {
  const int n = c.n, p = c.p;
  std::vector<double> xi(z(p) * z(p)), beta(z(p)), r(z(n)), work(z(p) * z(p) + z(p));
  if (ols(c.x.data(), c.y.data(), n, p, nullptr, beta.data(), xi.data(), r.data(), work.data()) !=
      0) {
    return false;
  }
  // h_i = x_i' (X'X)^-1 e_1: the first coefficient is sum h_i y_i
  c.h.assign(z(n), 0.0);
  for (int j = 0; j < p; ++j) {
    for (int i = 0; i < n; ++i) c.h[z(i)] += c.x[z(j) * z(n) + z(i)] * xi[z(j)];
  }
  c.xi = xi;
  // the projection on the other columns, for the restricted fits: a
  // subset of the columns of a nonsingular X'X is nonsingular too
  c.rxi.clear();
  if (p > 1) {
    std::vector<double> rb(z(p - 1)), rr(z(n));
    c.rxi.resize(z(p - 1) * z(p - 1));
    ols(c.x.data() + n, c.y.data(), n, p - 1, nullptr, rb.data(), c.rxi.data(), rr.data(),
        work.data());
  }
  return true;
}

double tstat(const wctx &c, const double *y, double b0, double *est) {
  const int n = c.n, p = c.p;
  std::vector<double> b(z(p), 0.0), xy(z(p), 0.0);
  for (int j = 0; j < p; ++j) {
    for (int i = 0; i < n; ++i) xy[z(j)] += c.x[z(j) * z(n) + z(i)] * y[i];
  }
  for (int j = 0; j < p; ++j) {
    for (int k = 0; k < p; ++k) b[z(j)] += c.xi[z(k) * z(p) + z(j)] * xy[z(k)];
  }
  std::vector<double> s(z(c.ng), 0.0);
  for (int i = 0; i < n; ++i) {
    double f = 0.0;
    for (int j = 0; j < p; ++j) f += c.x[z(j) * z(n) + z(i)] * b[z(j)];
    s[z(c.g[z(i)])] += c.h[z(i)] * (y[i] - f);
  }
  double v = 0.0;
  for (double q : s) v += q * q;
  if (est) *est = b[0];
  return (b[0] - b0) / std::sqrt(c.adj * v);
}

double wpval(const wctx &c, double b0, const std::vector<double> &v, int nd, int ties) {
  const int n = c.n, p = c.p;
  // the fit with b0 imposed: y - b0 x1 on the other columns
  std::vector<double> y0(z(n)), fr(z(n), 0.0), ur(z(n)), ys(z(n));
  for (int i = 0; i < n; ++i) y0[z(i)] = c.y[z(i)] - b0 * c.x[z(i)];
  if (p > 1) {
    std::vector<double> xy(z(p - 1), 0.0), b(z(p - 1), 0.0);
    for (int j = 1; j < p; ++j) {
      for (int i = 0; i < n; ++i) xy[z(j - 1)] += c.x[z(j) * z(n) + z(i)] * y0[z(i)];
    }
    for (int j = 0; j < p - 1; ++j) {
      for (int k = 0; k < p - 1; ++k) b[z(j)] += c.rxi[z(k) * z(p - 1) + z(j)] * xy[z(k)];
    }
    for (int i = 0; i < n; ++i) {
      for (int j = 1; j < p; ++j) fr[z(i)] += c.x[z(j) * z(n) + z(i)] * b[z(j - 1)];
    }
  }
  for (int i = 0; i < n; ++i) ur[z(i)] = y0[z(i)] - fr[z(i)];
  const double t0 = std::fabs(tstat(c, c.y.data(), b0, nullptr));
  int hit = 0;
  std::vector<double> wk(2 * z(c.nu + c.nt));
  for (int d = 0; d < nd; ++d) {
    if (kutil::stop(4.0 * n * p)) return kutil::nan();
    const double *w = v.data() + z(d) * z(c.ng);
    for (int i = 0; i < n; ++i) ys[z(i)] = ur[z(i)] * w[c.g[z(i)]];
    within(ys.data(), n, 1, c.u.data(), c.nu, c.t.data(), c.nt, 1, nullptr, 1e-12, 10000,
           wk.data());
    for (int i = 0; i < n; ++i) ys[z(i)] += fr[z(i)] + b0 * c.x[z(i)];
    const double ts = std::fabs(tstat(c, ys.data(), b0, nullptr));
    hit += ties ? ts >= t0 - 1e-12 : ts > t0 + 1e-12;
  }
  return static_cast<double>(hit) / nd;
}

double wend(const wctx &c, double bh, double s0, double alpha, int dir,
            const std::vector<double> &v, int nd, int ties) {
  double lo = bh, hi = bh + dir * s0;
  for (int k = 0; k < 12 && wpval(c, hi, v, nd, ties) > alpha; ++k) hi = bh + 2 * (hi - bh);
  if (wpval(c, hi, v, nd, ties) > alpha) return dir * HUGE_VAL;
  for (int k = 0; k < 40; ++k) {
    const double mid = 0.5 * (lo + hi);
    (wpval(c, mid, v, nd, ties) > alpha ? lo : hi) = mid;
  }
  return 0.5 * (lo + hi);
}

void wystep(const double *obs, const double *null, int nd, int k, double *p, double *pa) {
  // order by |obs| decreasing; successive maxima of the null draws over
  // the hypotheses from the j-th on
  std::vector<int> o(z(k));
  std::iota(o.begin(), o.end(), 0);
  std::stable_sort(o.begin(), o.end(),
                   [&](int a, int b) { return std::fabs(obs[a]) > std::fabs(obs[b]); });
  std::vector<double> cnt(z(k), 0.0), raw(z(k), 0.0);
  for (int d = 0; d < nd; ++d) {
    double q = 0.0;
    for (int j = k - 1; j >= 0; --j) {
      const int h = o[z(j)];
      const double a = std::fabs(null[z(h) * z(nd) + z(d)]);
      q = std::max(q, a);
      cnt[z(j)] += q >= std::fabs(obs[h]) - 1e-12 ? 1.0 : 0.0;
      raw[z(h)] += a >= std::fabs(obs[h]) - 1e-12 ? 1.0 : 0.0;
    }
  }
  double run = 0.0;
  for (int j = 0; j < k; ++j) {
    run = std::max(run, cnt[z(j)] / nd);
    pa[o[z(j)]] = run;
  }
  for (int h = 0; h < k; ++h) p[h] = raw[z(h)] / nd;
}

}  // namespace pev
