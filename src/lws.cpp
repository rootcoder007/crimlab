// SPDX-License-Identifier: AGPL-3.0-or-later
//
// libcrimlab, core: Cleveland's robust locally weighted regression
// (lowess), written from Cleveland (1979) with the conventions of the
// classic implementation (tricube neighbourhood of f n points, linear local
// fits, delta skipping with linear interpolation, bisquare robustness
// weights from six times the median absolute residual).
#include "crimlab.h"
#include "kutil.h"

namespace {

// k-th smallest of a[0..n) (0-based), partially reordering a; iterative
// quickselect with the middle element as pivot
double kth(double *a, int n, int k) {
  int lo = 0, hi = n - 1;
  while (lo < hi) {
    const double piv = a[lo + (hi - lo) / 2];
    int i = lo, j = hi;
    while (i <= j) {
      while (a[i] < piv) ++i;
      while (a[j] > piv) --j;
      if (i <= j) {
        const double t = a[i];
        a[i] = a[j];
        a[j] = t;
        ++i;
        --j;
      }
    }
    if (k <= j) hi = j;
    else if (k >= i) lo = i;
    else return a[k];
  }
  return a[k];
}

// the local fit at xs from points nl..nr (0-based, inclusive; points past
// nr tied within the bandwidth are included); w: n scratch weights.
// Returns false when every weight is zero.
bool local(const double *x, const double *y, int n, double xs, int nl, int nr,
           const double *rw, bool userw, double *w, double *ys) {
  const double range = x[n - 1] - x[0];
  double h = xs - x[nl];
  if (x[nr] - xs > h) h = x[nr] - xs;
  const double h9 = 0.999 * h, h1 = 0.001 * h;
  double a = 0.0;
  int j = nl;
  for (; j < n; ++j) {
    w[j] = 0.0;
    const double r = std::fabs(x[j] - xs);
    if (r <= h9) {
      if (r <= h1) {
        w[j] = 1.0;
      } else {
        const double q = r / h, t = 1.0 - q * q * q;
        w[j] = t * t * t;
      }
      if (userw) w[j] *= rw[j];
      a += w[j];
    } else if (x[j] > xs) {
      break;
    }
  }
  const int nrt = j - 1;
  if (!(a > 0.0)) return false;
  for (j = nl; j <= nrt; ++j) w[j] /= a;
  if (h > 0.0) {
    a = 0.0;
    for (j = nl; j <= nrt; ++j) a += w[j] * x[j];
    double b = xs - a, c = 0.0;
    for (j = nl; j <= nrt; ++j) c += w[j] * (x[j] - a) * (x[j] - a);
    if (std::sqrt(c) > 0.001 * range) {
      b /= c;
      for (j = nl; j <= nrt; ++j) w[j] *= b * (x[j] - a) + 1.0;
    }
  }
  double s = 0.0;
  for (j = nl; j <= nrt; ++j) s += w[j] * y[j];
  *ys = s;
  return true;
}

}  // namespace

extern "C" {

int clab_lwss(const double *x, const double *y, int n, double f, int iter,
              double delta, double *ys, double *work) {
  if (n < 1 || iter < 0) return CLAB_EARG;
  if (!(f > 0.0) || !std::isfinite(f) || !(delta >= 0.0) || !std::isfinite(delta))
    return CLAB_EDOM;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(x[i]) || !std::isfinite(y[i])) return CLAB_EDOM;
    if (i > 0 && x[i] < x[i - 1]) return CLAB_EDOM;
  }
  if (n == 1) {
    ys[0] = y[0];
    return CLAB_OK;
  }
  double *rw = work, *res = work + n, *w = work + 2 * n;
  int ns = static_cast<int>(f * n + 1e-7);
  if (ns > n) ns = n;
  if (ns < 2) ns = 2;
  for (int it = 0; it <= iter; ++it) {
    int nl = 0, nr = ns - 1, last = -1, i = 0;
    for (;;) {
      // slide the neighbourhood while its right end is closer to x[i]
      while (nr < n - 1) {
        const double d1 = x[i] - x[nl], d2 = x[nr + 1] - x[i];
        if (d1 <= d2) break;
        ++nl;
        ++nr;
      }
      if (!local(x, y, n, x[i], nl, nr, rw, it > 0, w, &ys[i])) ys[i] = y[i];
      // linear interpolation over the points skipped since the last fit
      if (last < i - 1) {
        const double den = x[i] - x[last];
        for (int j = last + 1; j < i; ++j) {
          const double al = (x[j] - x[last]) / den;
          ys[j] = al * ys[i] + (1.0 - al) * ys[last];
        }
      }
      last = i;
      const double cut = x[last] + delta;
      for (i = last + 1; i < n; ++i) {
        if (x[i] > cut) break;
        if (x[i] == x[last]) {
          ys[i] = ys[last];
          last = i;
        }
      }
      i = last + 1 > i - 1 ? last + 1 : i - 1;
      if (last >= n - 1) break;
    }
    for (int k = 0; k < n; ++k) res[k] = y[k] - ys[k];
    if (it == iter) break;
    double sc = 0.0;
    for (int k = 0; k < n; ++k) sc += std::fabs(res[k]);
    sc /= n;
    for (int k = 0; k < n; ++k) rw[k] = std::fabs(res[k]);
    // six times the median absolute residual
    const int m1 = n / 2;
    double cmad;
    if (n % 2 == 0) {
      const double hi = kth(rw, n, m1);
      const double lo = kth(rw, n, n - m1 - 1);
      cmad = 3.0 * (hi + lo);
    } else {
      cmad = 6.0 * kth(rw, n, m1);
    }
    if (cmad < 1e-7 * sc) break;
    const double c9 = 0.999 * cmad, c1 = 0.001 * cmad;
    for (int k = 0; k < n; ++k) {
      const double r = std::fabs(res[k]);
      if (r <= c1) {
        rw[k] = 1.0;
      } else if (r > c9) {
        rw[k] = 0.0;
      } else {
        const double q = r / cmad, t = 1.0 - q * q;
        rw[k] = t * t;
      }
    }
  }
  return CLAB_OK;
}

}  // extern "C"
