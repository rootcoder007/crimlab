// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Geographically and temporally weighted regression (crimlab.h,
// gtwreg()): at each observation the least-squares fit of all of them
// weighted by a Gaussian kernel of the space-time distance, d^2 = ds^2 +
// tau dt^2 and w = exp(-d^2 / h^2) (Kuo and Putra 2021, eqs 5-8, after
// Huang, Wu and Barry 2010); the fitted values, the hat matrix's diagonal
// and the trace of S'S, the leave-one-out residuals of the
// cross-validation score (eq 9), the corrected Akaike criterion and the
// local standard errors. Work: O(n^2 k) per fit, n the observations.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "crimlab.h"
#include "kutil.h"
#include "pev.h"

namespace {

using std::size_t;
constexpr double rdn = 3.14159265358979323846 / 180.0;

inline size_t z(int i) { return static_cast<size_t>(i); }

double sdist2(double x0, double y0, double x1, double y1, int ll) {
  if (!ll) return (x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0);
  const double p0 = y0 * rdn, p1 = y1 * rdn;
  const double a = std::sin((p1 - p0) / 2), c = std::sin((x1 - x0) * rdn / 2);
  const double h = a * a + std::cos(p0) * std::cos(p1) * c * c;
  const double d = 2.0 * 6371.0088 * std::asin(std::sqrt(std::min(1.0, h)));
  return d * d;
}

struct gin {
  int n, k, ll;
  const double *x, *y, *u, *v, *t;
  double h, tau;
};

// The weights of the fit at observation i.
void wts(const gin &g, int i, std::vector<double> &w) {
  for (int j = 0; j < g.n; ++j) {
    const double dt = g.t ? g.t[j] - g.t[i] : 0.0;
    const double d2 = sdist2(g.u[i], g.v[i], g.u[j], g.v[j], g.ll) + g.tau * dt * dt;
    w[z(j)] = std::exp(-d2 / (g.h * g.h));
  }
}

// One local fit: A = X'WX (Cholesky inverse into ai), b = A^-1 X'Wy and
// B = X'W^2X; false when A is singular.
bool fit(const gin &g, const std::vector<double> &w, std::vector<double> &ai,
         std::vector<double> &b, std::vector<double> &bb) {
  const int k = g.k, n = g.n;
  std::vector<double> a(z(k) * z(k), 0.0), xy(z(k), 0.0);
  std::fill(bb.begin(), bb.end(), 0.0);
  for (int j = 0; j < n; ++j) {
    const double wj = w[z(j)];
    if (wj == 0.0) continue;
    for (int c = 0; c < k; ++c) {
      const double xc = g.x[z(c) * z(n) + z(j)];
      xy[z(c)] += wj * xc * g.y[j];
      for (int r = 0; r < k; ++r) {
        const double p = xc * g.x[z(r) * z(n) + z(j)];
        a[z(c) * z(k) + z(r)] += wj * p;
        bb[z(c) * z(k) + z(r)] += wj * wj * p;
      }
    }
  }
  if (!pev::chol(a.data(), k)) return false;
  pev::cinv(a.data(), k, ai.data());
  for (int r = 0; r < k; ++r) {
    double s = 0.0;
    for (int c = 0; c < k; ++c) s += ai[z(c) * z(k) + z(r)] * xy[z(c)];
    b[z(r)] = s;
  }
  return true;
}

}  // namespace

extern "C" int clab_gtwreg(int n, int k, const double *x, const double *y, const double *u,
                           const double *v, const double *t, double h, double tau, int ll,
                           double *beta, double *se, double *fitv, double *hat, double *loo,
                           double *st) CLAB_TRY {
  if (n < 2 || k < 1 || k >= n || !(h > 0.0) || !std::isfinite(h) || !(tau >= 0.0) ||
      !std::isfinite(tau)) {
    return CLAB_EARG;
  }
  if (static_cast<double>(n) * k > kutil::maxcell) return CLAB_EARG;
  for (size_t i = 0; i < z(n) * z(k); ++i) {
    if (!std::isfinite(x[i])) return CLAB_EDOM;
  }
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(y[i]) || !std::isfinite(u[i]) || !std::isfinite(v[i]) ||
        (t && !std::isfinite(t[i]))) {
      return CLAB_EDOM;
    }
  }
  const gin g{n, k, ll, x, y, u, v, t, h, tau};
  std::vector<double> w(z(n)), ai(z(k) * z(k)), b(z(k)), bb(z(k) * z(k)), xi(z(k)), q(z(k));
  std::vector<double> vb(z(n) * z(k));  // the diagonal of A^-1 B A^-1, by observation
  double rss = 0.0, trs = 0.0, trss = 0.0;
  for (int i = 0; i < n; ++i) {
    if (kutil::stop(static_cast<double>(n) * k * (k + 1))) return CLAB_EINTR;
    wts(g, i, w);
    if (!fit(g, w, ai, b, bb)) return CLAB_EDOM;
    double f = 0.0;
    for (int c = 0; c < k; ++c) {
      xi[z(c)] = x[z(c) * z(n) + z(i)];
      beta[z(c) * z(n) + z(i)] = b[z(c)];
      f += xi[z(c)] * b[z(c)];
    }
    // q = A^-1 x_i: S_ij = (x_j . q) w_j
    for (int r = 0; r < k; ++r) {
      double s = 0.0;
      for (int c = 0; c < k; ++c) s += ai[z(c) * z(k) + z(r)] * xi[z(c)];
      q[z(r)] = s;
    }
    double sii = 0.0;
    for (int j = 0; j < n; ++j) {
      double s = 0.0;
      for (int c = 0; c < k; ++c) s += x[z(c) * z(n) + z(j)] * q[z(c)];
      s *= w[z(j)];
      trss += s * s;
      if (j == i) sii = s;
    }
    // the diagonal of A^-1 B A^-1
    for (int r = 0; r < k; ++r) {
      double s = 0.0;
      for (int c = 0; c < k; ++c) {
        for (int d = 0; d < k; ++d) {
          s += ai[z(c) * z(k) + z(r)] * bb[z(d) * z(k) + z(c)] * ai[z(r) * z(k) + z(d)];
        }
      }
      vb[z(r) * z(n) + z(i)] = s;
    }
    fitv[i] = f;
    hat[i] = sii;
    const double e = y[i] - f;
    loo[i] = sii < 1.0 ? e / (1.0 - sii) : NAN;
    rss += e * e;
    trs += sii;
  }
  const double edf = n - 2.0 * trs + trss;
  const double s2 = edf > 0.0 ? rss / edf : NAN;
  double cv = 0.0;
  for (int i = 0; i < n; ++i) cv += loo[i] * loo[i];
  for (size_t c = 0; c < z(n) * z(k); ++c) se[c] = std::sqrt(s2 * vb[c]);
  const double nn = n;
  st[0] = rss;
  st[1] = trs;
  st[2] = trss;
  st[3] = s2;
  st[4] = cv;
  // the corrected Akaike criterion of the GWR literature: n log(RSS / n) +
  // n log(2 pi) + n (n + tr S) / (n - 2 - tr S)
  const double l2p = std::log(2.0 * 3.14159265358979323846);
  st[5] = nn - 2.0 - trs > 0.0
              ? nn * std::log(rss / nn) + nn * l2p + nn * (nn + trs) / (nn - 2.0 - trs)
              : NAN;
  st[6] = edf;
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)
