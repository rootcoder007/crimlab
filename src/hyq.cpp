// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, police operations: Larson's hypercube queue (clab_hypq) and
// the Kolesar-Blum square-root law of response distance (clab_sqrtlw). No
// allocation: the caller passes the work arrays.
#include <cmath>
#include <cstddef>

#include "crimlab.h"
#include "kutil.h"

namespace {

constexpr double pi = 3.14159265358979323846;

// the first free unit in atom j's order (b is not the full state)
int first(unsigned b, const int *pref, int j, int m) {
  int r = 0;
  while ((b >> pref[j + r * m]) & 1u) ++r;
  return pref[j + r * m];
}

}  // namespace

extern "C" {

int clab_hypq(const double *lam, const int *pref, int m, int n, double mu, int wait,
              double *pst, double *wl, double *disp, double *st, double *work) {
  if (m < 1 || n < 1 || n > 10 || !kutil::pos(mu)) return CLAB_EARG;
  double big = 0;
  for (int j = 0; j < m; ++j) {
    if (!kutil::pos(lam[j])) return CLAB_EDOM;
    big += lam[j];
    // each atom's order: every unit once
    unsigned seen = 0;
    for (int r = 0; r < n; ++r) {
      const int u = pref[j + r * m];
      if (u < 0 || u >= n || ((seen >> u) & 1u)) return CLAB_EDOM;
      seen |= 1u << u;
    }
  }
  const double rho = big / (n * mu);
  if (wait && rho >= 1) return CLAB_EDOM;
  const int s = 1 << n;
  const unsigned full = static_cast<unsigned>(s - 1);
  const std::size_t ss = static_cast<std::size_t>(s);
  // the transposed generator Q' (column b holds the rates out of b), the
  // last row replaced by the normalisation
  double *a = work, *x = work + ss * ss;
  for (std::size_t i = 0; i < ss * ss; ++i) a[i] = 0;
  for (unsigned b = 0; b < full + 1; ++b) {
    double out = 0;
    if (b != full) {
      for (int j = 0; j < m; ++j) {
        const unsigned c = b | (1u << first(b, pref, j, m));
        a[c + b * ss] += lam[j];
        out += lam[j];
      }
    }
    // the infinite line lumps every state with a queue into the full one:
    // it empties a unit at rate n mu (1 - rho) in all (Larson 1974)
    const double sc = b == full && wait ? 1 - rho : 1;
    for (int u = 0; u < n; ++u) {
      if ((b >> u) & 1u) {
        a[(b & ~(1u << u)) + b * ss] += mu * sc;
        out += mu * sc;
      }
    }
    a[b + b * ss] -= out;
  }
  for (std::size_t j = 0; j < ss; ++j) {
    a[(ss - 1) + j * ss] = 1;
    x[j] = j == ss - 1 ? 1 : 0;
  }
  // Gaussian elimination with partial pivoting
  for (std::size_t k = 0; k < ss; ++k) {
    std::size_t pv = k;
    for (std::size_t i = k + 1; i < ss; ++i) {
      if (std::fabs(a[i + k * ss]) > std::fabs(a[pv + k * ss])) pv = i;
    }
    for (std::size_t j = k; j < ss; ++j) {
      const double t = a[k + j * ss];
      a[k + j * ss] = a[pv + j * ss];
      a[pv + j * ss] = t;
    }
    const double t = x[k];
    x[k] = x[pv];
    x[pv] = t;
    for (std::size_t i = k + 1; i < ss; ++i) {
      const double f = a[i + k * ss] / a[k + k * ss];
      for (std::size_t j = k; j < ss; ++j) a[i + j * ss] -= f * a[k + j * ss];
      x[i] -= f * x[k];
    }
  }
  for (std::size_t k = ss; k-- > 0;) {
    double v = x[k];
    for (std::size_t j = k + 1; j < ss; ++j) v -= a[k + j * ss] * pst[j];
    pst[k] = v / a[k + k * ss];
  }
  // workloads and dispatch shares (of all calls)
  for (int u = 0; u < n; ++u) wl[u] = 0;
  for (int i = 0; i < n * m; ++i) disp[i] = 0;
  for (unsigned b = 0; b < full + 1; ++b) {
    for (int u = 0; u < n; ++u) wl[u] += ((b >> u) & 1u) ? pst[b] : 0;
    if (b == full) continue;
    for (int j = 0; j < m; ++j) disp[first(b, pref, j, m) + j * n] += lam[j] * pst[b] / big;
  }
  st[0] = pst[full];
  st[1] = wait ? pst[full] * rho / (1 - rho) : 0;
  st[2] = rho;
  return CLAB_OK;
}

int clab_sqrtlw(double area, double units, double busy, int metric, double *d, double *c) {
  if (!kutil::pos(area) || !kutil::pos(units) || !(busy >= 0 && busy < 1) || metric < 0 ||
      metric > 1) {
    return CLAB_EDOM;
  }
  // the mean distance to the nearest of a Poisson field of available
  // units: 1 / 2 for straight lines, sqrt(2 pi) / 4 on a grid of streets
  *c = metric ? std::sqrt(2 * pi) / 4 : 0.5;
  *d = *c * std::sqrt(area / (units * (1 - busy)));
  return CLAB_OK;
}

}  // extern "C"
