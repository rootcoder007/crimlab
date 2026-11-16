// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Least squares with an intercept for the experiment and trial kernels
// (fls.cpp and the files after it): the projection on [1, x] by the
// normal equations of the weighted, centred columns (Cholesky, gla.h),
// and the coefficient of one more column with its variance by
// Frisch-Waugh-Lovell. Matrices are column-major. Weights are frequency
// weights from 0: a row of weight 2 counts as two rows.
#ifndef CLAB_XLS_H
#define CLAB_XLS_H

#include <cmath>
#include <cstddef>
#include <vector>

#include "gla.h"

namespace xls {

inline std::size_t z(int i) { return static_cast<std::size_t>(i); }

// The projection on [1, x[, cols]] (cols NULL: the first k columns of x,
// n rows), weights w (NULL: 1). ok is false when the weights sum to 0, or
// a column after the ones before it keeps less than 1e-10 of its centred
// sum of squares (a constant or collinear column).
struct Proj {
  int n = 0, k = 0;
  const double *w = nullptr;
  double sw = 0;
  std::vector<double> xc, mu, l;
  bool ok = false;

  bool make(const double *x, int nn, const int *cols, int kk, const double *ww) {
    n = nn;
    k = kk;
    w = ww;
    sw = 0;
    for (int i = 0; i < n; ++i) sw += w ? w[i] : 1.0;
    ok = sw > 0;
    if (!ok) return false;
    xc.assign(z(n) * z(k), 0.0);
    mu.assign(z(k), 0.0);
    l.assign(z(k) * z(k), 0.0);
    for (int j = 0; j < k; ++j) {
      const double *xj = x + z(cols ? cols[j] : j) * z(n);
      double s = 0;
      for (int i = 0; i < n; ++i) s += (w ? w[i] : 1.0) * xj[i];
      mu[z(j)] = s / sw;
      for (int i = 0; i < n; ++i) xc[z(j) * z(n) + z(i)] = xj[i] - mu[z(j)];
    }
    for (int a = 0; a < k; ++a) {
      for (int b = a; b < k; ++b) {
        double s = 0;
        const double *xa = &xc[z(a) * z(n)], *xb = &xc[z(b) * z(n)];
        for (int i = 0; i < n; ++i) s += (w ? w[i] : 1.0) * xa[i] * xb[i];
        l[z(a) + z(b) * z(k)] = s;
        l[z(b) + z(a) * z(k)] = s;
      }
    }
    if (k == 0) return true;
    std::vector<double> dg(z(k));
    for (int j = 0; j < k; ++j) dg[z(j)] = l[z(j) * (z(k) + 1)];
    if (!gla::chol(l.data(), k)) return ok = false;
    for (int j = 0; j < k; ++j) {
      const double lj = l[z(j) * (z(k) + 1)];
      if (!(dg[z(j)] > 0) || lj * lj < 1e-10 * dg[z(j)]) return ok = false;
    }
    return true;
  }

  // r = v less its fitted values; b (k + 1 values, the intercept first)
  // the coefficients, when given.
  void res(const double *v, double *r, double *b = nullptr) const {
    double vm = 0;
    for (int i = 0; i < n; ++i) vm += (w ? w[i] : 1.0) * v[i];
    vm /= sw;
    std::vector<double> c(z(k), 0.0);
    for (int j = 0; j < k; ++j) {
      const double *xj = &xc[z(j) * z(n)];
      double s = 0;
      for (int i = 0; i < n; ++i) s += (w ? w[i] : 1.0) * xj[i] * (v[i] - vm);
      c[z(j)] = s;
    }
    if (k > 0) {
      gla::fsolve(l.data(), k, c.data(), 1);
      gla::bsolve(l.data(), k, c.data(), 1);
    }
    for (int i = 0; i < n; ++i) {
      double f = vm;
      for (int j = 0; j < k; ++j) f += c[z(j)] * xc[z(j) * z(n) + z(i)];
      r[i] = v[i] - f;
    }
    if (b) {
      b[0] = vm;
      for (int j = 0; j < k; ++j) {
        b[0] -= c[z(j)] * mu[z(j)];
        b[j + 1] = c[z(j)];
      }
    }
  }
};

// The coefficient of d in the regression of y on [1, x, d] from rd and
// ry, the residuals of d and y on [1, x] (Frisch-Waugh-Lovell): b = sum w
// rd ry / sum w rd^2, the residuals e = ry - b rd, and the variance of b,
// HC1 (hc) n / (n - q) sum w rd^2 e^2 / (sum w rd^2)^2 or classical sum w
// e^2 / (n - q) / sum w rd^2, with n = sum w and q the coefficients with
// the intercept and d; df = n - q. False when d keeps less than 1e-10 of
// its centred sum of squares dss after x (b, v and df NaN).
inline bool coef(const double *rd, const double *ry, const double *w, int n, int q, bool hc, double dss,
                 double &b, double &v, double &df) {
  b = v = df = std::nan("");
  double sw = 0, srr = 0, sry = 0;
  for (int i = 0; i < n; ++i) {
    const double wi = w ? w[i] : 1.0;
    sw += wi;
    srr += wi * rd[i] * rd[i];
    sry += wi * rd[i] * ry[i];
  }
  if (!(srr >= 1e-10 * dss) || !(srr > 0)) return false;
  b = sry / srr;
  df = sw - q;
  double s = 0;
  for (int i = 0; i < n; ++i) {
    const double wi = w ? w[i] : 1.0, e = ry[i] - b * rd[i];
    s += hc ? wi * rd[i] * rd[i] * e * e : wi * e * e;
  }
  if (!(df > 0)) {
    v = std::nan("");
  } else {
    v = hc ? sw / df * s / (srr * srr) : s / df / srr;
  }
  return true;
}

// The weighted centred sum of squares of v.
inline double css(const double *v, const double *w, int n) {
  double sw = 0, m = 0, s = 0;
  for (int i = 0; i < n; ++i) {
    sw += w ? w[i] : 1.0;
    m += (w ? w[i] : 1.0) * v[i];
  }
  m /= sw;
  for (int i = 0; i < n; ++i) s += (w ? w[i] : 1.0) * (v[i] - m) * (v[i] - m);
  return s;
}

}  // namespace xls

#endif
