// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The drift battery (R/drift.R): has the distribution of a column moved
// between two releases of a table? The two-sample Kolmogorov-Smirnov
// statistic with its asymptotic p-value, the quantile bins of the
// population stability index, the index itself and the Jensen-Shannon
// divergence, and the chi-square sums of goodness of fit and of
// homogeneity. The p-values of the chi-square tests are left to the
// caller.
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

namespace {

const double qnan = std::numeric_limits<double>::quiet_NaN();
const double pi = 3.14159265358979323846;

// 10^k by repeated products, as R_pow_di()
double pow10i(int k) {
  double r = 1, b = 10;
  for (int n = k; n; n >>= 1, b *= b) {
    if (n & 1) r *= b;
  }
  return r;
}

// signif(x, 12), as R's fprec()
double sig12(double x) {
  if (x == 0 || !std::isfinite(x)) return x;
  const double sgn = x < 0 ? -1 : 1, a = std::fabs(x), l10 = std::log10(a);
  int e10 = static_cast<int>(11 - std::floor(l10));
  if (std::fabs(l10) >= 308 - 2) return x;
  double p10 = 1;
  if (e10 > 308) {
    p10 = pow10i(e10 - 308);
    e10 = 308;
  }
  if (e10 > 0) {
    const double p = pow10i(e10);
    return sgn * (std::nearbyint(a * p * p10) / p) / p10;
  }
  const double p = pow10i(-e10);
  return sgn * (std::nearbyint(a / p) * p);
}

double ksp(double d, double ne) {
  if (!(d > 0.0) || !(ne > 0.0)) return 1.0;
  const double t = std::sqrt(ne) * d;
  if (t < 1e-12) return 1.0;
  double p;
  if (t < 1.0) {
    const double z = -(pi * pi / 8.0) / (t * t), w = std::log(t);
    double s = 0.0;
    for (int k = 1; k < 200; k += 2) {
      const double term = std::exp(static_cast<double>(k) * static_cast<double>(k) * z - w);
      s += term;
      if (term < 1e-300) break;
    }
    p = 1.0 - std::sqrt(2.0 * pi) * s;
  } else {
    double s = 0.0;
    for (int k = 1; k <= 200; ++k) {
      const double term = std::exp(-2.0 * static_cast<double>(k) * static_cast<double>(k) * t * t);
      s += (k % 2 == 1 ? 1.0 : -1.0) * term;
      if (term < 1e-300) break;
    }
    p = 2.0 * s;
  }
  return std::min(1.0, std::max(0.0, p));
}

}  // namespace

extern "C" {

void clab_ks(const double *x, int nx, const double *y, int ny, double *out) CLAB_TRY {
  out[2] = nx > 0 && ny > 0 ? static_cast<double>(nx) * ny / (static_cast<double>(nx) + ny) : 0.0;
  bool nan = nx < 1 || ny < 1;
  for (int i = 0; i < nx && !nan; ++i) nan = std::isnan(x[i]);
  for (int i = 0; i < ny && !nan; ++i) nan = std::isnan(y[i]);
  if (nan) {
    out[0] = out[1] = qnan;
    return;
  }
  std::vector<double> a(x, x + nx), b(y, y + ny);
  std::sort(a.begin(), a.end());
  std::sort(b.begin(), b.end());
  int i = 0, j = 0;
  double d = 0.0;
  while (i < nx && j < ny) {
    const double v = std::min(a[static_cast<std::size_t>(i)], b[static_cast<std::size_t>(j)]);
    while (i < nx && a[static_cast<std::size_t>(i)] == v) ++i;
    while (j < ny && b[static_cast<std::size_t>(j)] == v) ++j;
    d = std::max(d, std::fabs(static_cast<double>(i) / nx - static_cast<double>(j) / ny));
  }
  out[0] = d;
  out[1] = ksp(d, out[2]);
} CLAB_CATCHV

int clab_psibin(const double *x, int nx, const double *y, int ny, int bins, double *px, double *py) CLAB_TRY {
  std::vector<double> a(static_cast<std::size_t>(nx)), b(static_cast<std::size_t>(ny));
  for (int i = 0; i < nx; ++i) a[static_cast<std::size_t>(i)] = sig12(x[i]);
  for (int i = 0; i < ny; ++i) b[static_cast<std::size_t>(i)] = sig12(y[i]);
  std::vector<double> p(static_cast<std::size_t>(bins) + 1), q(p.size());
  for (int k = 0; k <= bins; ++k) {
    // seq(0, 1, length.out = bins + 1)
    p[static_cast<std::size_t>(k)] = k == bins ? 1.0 : k * (1.0 / bins);
  }
  clab_quant(a.data(), nx, p.data(), bins + 1, q.data());
  std::vector<double> u;
  for (double v : q) {
    if (std::find(u.begin(), u.end(), v) == u.end()) u.push_back(v);
  }
  std::vector<double> e = {-HUGE_VAL};
  if (u.size() >= 3) {
    e.insert(e.end(), u.begin() + 1, u.end() - 1);
  } else {
    e.insert(e.end(), u.begin(), u.end());
  }
  e.push_back(HUGE_VAL);
  const int k = static_cast<int>(e.size()) - 1;
  for (int j = 0; j < k; ++j) px[j] = py[j] = 0;
  auto bin = [&e, k](double v) {
    for (int j = 0; j < k; ++j) {
      if (v > e[static_cast<std::size_t>(j)] && v <= e[static_cast<std::size_t>(j) + 1]) return j;
    }
    return k - 1;
  };
  for (double v : a) px[bin(v)] += 1;
  for (double v : b) py[bin(v)] += 1;
  for (int j = 0; j < k; ++j) {
    px[j] /= nx;
    py[j] /= ny;
  }
  return k;
} CLAB_CATCH(CLAB_ENOMEM)

void clab_psi(const double *p, const double *q, int k, double eps, double *out) CLAB_TRY {
  out[0] = out[1] = qnan;
  if (k < 1) return;
  double s = 0.0, js = 0.0;
  for (int i = 0; i < k; ++i) {
    const double a = p[i] > eps ? p[i] : eps, b = q[i] > eps ? q[i] : eps;
    s += (a - b) * std::log(a / b);
    const double m = 0.5 * (p[i] + q[i]);
    if (p[i] > 0.0) js += 0.5 * p[i] * std::log(p[i] / m);
    if (q[i] > 0.0) js += 0.5 * q[i] * std::log(q[i] / m);
  }
  out[0] = s;
  out[1] = js < 0.0 ? 0.0 : js;
} CLAB_CATCHV

void clab_chig(const double *o, const double *e, int k, double *out) CLAB_TRY {
  double so = 0, se = 0;
  for (int j = 0; j < k; ++j) {
    so += o[j];
    se += e[j];
  }
  double stat = 0, df = -1, bad = 0;
  for (int j = 0; j < k; ++j) {
    const double ex = so * e[j] / se;
    if (ex > 0) {
      df += 1;
      stat += (o[j] - ex) * (o[j] - ex) / ex;
    } else if (o[j] > 0) {
      bad = 1;
    }
  }
  out[0] = stat;
  out[1] = df;
  out[2] = bad;
} CLAB_CATCHV

void clab_chih(const double *x, const double *y, int k, double *out) CLAB_TRY {
  double nx = 0, ny = 0;
  for (int j = 0; j < k; ++j) {
    nx += x[j];
    ny += y[j];
  }
  double sx = 0, sy = 0, df = -1, mn = HUGE_VAL;
  for (int j = 0; j < k; ++j) {
    const double p = (x[j] + y[j]) / (nx + ny), ex = nx * p, ey = ny * p;
    if (p > 0) {
      df += 1;
      sx += (x[j] - ex) * (x[j] - ex) / ex;
      sy += (y[j] - ey) * (y[j] - ey) / ey;
      mn = std::min(mn, std::min(ex, ey));
    }
  }
  out[0] = sx + sy;
  out[1] = df;
  out[2] = mn;
  out[3] = df + 1;
} CLAB_CATCHV

}  // extern "C"
