// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Summary statistics in compiled form (R/prim.R): the mean, variance and
// correlation of base R, the normal density, moments, type-7 quantiles,
// the median and its absolute deviation, trimmed and winsorized means,
// weighted moments, ranks and Spearman's rho, the covariance matrix,
// bootstrap replicate means, trimmed inverse-probability weights and the
// splitmix64 uniforms. Missing values propagate: there is no na.rm. The
// rules for each one are those of crimlab.h.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

namespace {

const double qnan = std::numeric_limits<double>::quiet_NaN();

bool hasnan(const double *a, int n) {
  for (int i = 0; i < n; ++i) {
    if (std::isnan(a[i])) return true;
  }
  return false;
}

// A running mean: it cannot overflow on finite input (each term is at most
// the largest |a|).
double runmean(const double *a, int n) {
  double m = 0.0;
  for (int i = 0; i < n; ++i) {
    const double k = i + 1.0;
    m += a[i] / k - m / k;
  }
  return m;
}

// The type-7 quantile of the sorted v at p: the lower value alone when
// the weight on the upper one is 0, a weighted mean when their difference
// is not finite.
double q7s(const std::vector<double> &v, double p) {
  const int n = static_cast<int>(v.size());
  if (n == 0 || std::isnan(p)) return n == 0 ? qnan : p;
  if (n == 1) return v[0];
  const double h = (n - 1.0) * p, fl = std::floor(h);
  int lo = static_cast<int>(fl);
  if (lo < 0) lo = 0;
  if (lo >= n - 1) return v[static_cast<std::size_t>(n - 1)];
  const double g = h - fl, a = v[static_cast<std::size_t>(lo)], b = v[static_cast<std::size_t>(lo) + 1];
  if (g == 0.0) return a;
  const double d = b - a;
  return std::isfinite(d) ? a + g * d : (1.0 - g) * a + g * b;
}

std::vector<double> sorted(const double *a, int n) {
  std::vector<double> v(a, a + n);
  std::sort(v.begin(), v.end());
  return v;
}

// Midranks of a[0..n), none NaN: tied values share the mean of their ranks.
void midr(const double *a, int n, double *out) {
  std::vector<int> o(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) o[static_cast<std::size_t>(i)] = i;
  std::sort(o.begin(), o.end(), [a](int i, int j) { return a[i] < a[j]; });
  int i = 0;
  while (i < n) {
    int j = i;
    while (j + 1 < n && a[o[static_cast<std::size_t>(j) + 1]] == a[o[static_cast<std::size_t>(i)]]) ++j;
    const double r = (i + j + 2.0) / 2.0;
    for (int k = i; k <= j; ++k) out[o[static_cast<std::size_t>(k)]] = r;
    i = j + 1;
  }
}

// An unbiased draw from [0, s), s >= 1: Lemire's multiply-and-reject on the
// raw outputs of mt19937_64 (whose sequence the C++ standard fixes), the
// 64 by 64 bit product written out.
std::uint64_t bdraw(std::mt19937_64 &g, std::uint64_t s) {
  auto mul = [](std::uint64_t x, std::uint64_t y, std::uint64_t &hi, std::uint64_t &lo) {
    const std::uint64_t x0 = x & 0xFFFFFFFFULL, x1 = x >> 32, y0 = y & 0xFFFFFFFFULL, y1 = y >> 32;
    const std::uint64_t p00 = x0 * y0, p01 = x0 * y1, p10 = x1 * y0, p11 = x1 * y1;
    const std::uint64_t mid = (p00 >> 32) + (p01 & 0xFFFFFFFFULL) + (p10 & 0xFFFFFFFFULL);
    hi = p11 + (p01 >> 32) + (p10 >> 32) + (mid >> 32);
    lo = (mid << 32) | (p00 & 0xFFFFFFFFULL);
  };
  std::uint64_t hi = 0, lo = 0;
  mul(static_cast<std::uint64_t>(g()), s, hi, lo);
  if (lo < s) {
    const std::uint64_t t = (0 - s) % s;
    while (lo < t) mul(static_cast<std::uint64_t>(g()), s, hi, lo);
  }
  return hi;
}

}  // namespace

extern "C" {

double clab_mean(const double *x, int n) CLAB_TRY {
  if (n <= 0) return qnan;
  long double s = 0.0L;
  for (int i = 0; i < n; ++i) s += x[i];
  s /= static_cast<long double>(n);
  if (std::isfinite(static_cast<double>(s))) {
    long double t = 0.0L;
    for (int i = 0; i < n; ++i) t += (x[i] - s);
    s += t / static_cast<long double>(n);
    return static_cast<double>(s);
  }
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(x[i])) return static_cast<double>(s);
  }
  return runmean(x, n);
} CLAB_CATCH(kutil::nan())

double clab_var(const double *x, int n, int ddof) CLAB_TRY {
  if (static_cast<long long>(n) - ddof <= 0) return qnan;
  const double m = clab_mean(x, n);
  double sq = 0.0;
  for (int i = 0; i < n; ++i) {
    const double d = x[i] - m;
    sq += d * d;
  }
  return sq / (static_cast<double>(n) - static_cast<double>(ddof));
} CLAB_CATCH(kutil::nan())

double clab_cor(const double *x, const double *y, int n) CLAB_TRY {
  if (n < 2) return qnan;
  const double mx = clab_mean(x, n), my = clab_mean(y, n);
  if (!std::isfinite(mx) || !std::isfinite(my)) return qnan;
  long double sxy = 0.0L, sxx = 0.0L, syy = 0.0L;
  for (int i = 0; i < n; ++i) {
    const long double a = x[i] - mx, b = y[i] - my;
    sxy += a * b;
    sxx += a * a;
    syy += b * b;
  }
  if (sxx <= 0.0L || syy <= 0.0L) return qnan;
  long double r = sxy / std::sqrt(sxx * syy);
  if (r > 1.0L) r = 1.0L;
  if (r < -1.0L) r = -1.0L;
  return static_cast<double>(r);
} CLAB_CATCH(kutil::nan())

double clab_edist(const double *a, const double *b, int n) CLAB_TRY {
  double s = 0.0;
  for (int i = 0; i < n; ++i) {
    const double d = a[i] - b[i];
    s += d * d;
  }
  return std::sqrt(s);
} CLAB_CATCH(kutil::nan())

void clab_dnorm(const double *x, int n, double mu, double sd, int lg, double *out) CLAB_TRY {
  const double inv = 1.0 / sd;
  const double c = 0.39894228040143267794, base = -std::log(sd) - 0.91893853320467274178;
  for (int i = 0; i < n; ++i) {
    const double z = (x[i] - mu) * inv;
    out[i] = lg ? base - 0.5 * z * z : inv * c * std::exp(-0.5 * z * z);
  }
} CLAB_CATCHV

void clab_moms(const double *x, int n, double *out) CLAB_TRY {
  out[0] = out[1] = out[2] = out[3] = qnan;
  if (n < 1 || hasnan(x, n)) return;
  const double m = clab_mean(x, n);
  out[0] = m;
  if (n < 2) return;
  double s = 0.0;
  for (int i = 0; i < n; ++i) s = std::max(s, std::fabs(x[i] - m));
  if (s == 0.0) {
    out[1] = 0.0;
    return;
  }
  long double z2 = 0.0L, z3 = 0.0L, z4 = 0.0L;
  for (int i = 0; i < n; ++i) {
    const long double z = (x[i] - m) / s, zz = z * z;
    z2 += zz;
    z3 += zz * z;
    z4 += zz * zz;
  }
  const long double dn = n;
  out[1] = static_cast<double>(z2 / (dn - 1.0L)) * s * s;
  const long double m2 = z2 / dn;
  if (n > 2) out[2] = static_cast<double>((z3 / dn) / std::pow(m2, 1.5L));
  if (n > 3) out[3] = static_cast<double>((z4 / dn) / (m2 * m2) - 3.0L);
} CLAB_CATCHV

void clab_quant(const double *x, int n, const double *p, int np, double *out) CLAB_TRY {
  if (n < 1 || hasnan(x, n)) {
    for (int k = 0; k < np; ++k) out[k] = qnan;
    return;
  }
  const std::vector<double> v = sorted(x, n);
  for (int k = 0; k < np; ++k) out[k] = q7s(v, p[k]);
} CLAB_CATCHV

double clab_med(const double *x, int n) CLAB_TRY {
  double r = 0.0;
  const double h = 0.5;
  clab_quant(x, n, &h, 1, &r);
  return r;
} CLAB_CATCH(kutil::nan())

double clab_mad(const double *x, int n, double k) CLAB_TRY {
  if (n < 1 || hasnan(x, n)) return qnan;
  std::vector<double> v = sorted(x, n);
  const double med = q7s(v, 0.5);
  for (int i = 0; i < n; ++i) v[static_cast<std::size_t>(i)] = std::fabs(x[i] - med);
  std::sort(v.begin(), v.end());
  return k * q7s(v, 0.5);
} CLAB_CATCH(kutil::nan())

double clab_trim(const double *x, int n, double t) CLAB_TRY {
  if (n < 1 || hasnan(x, n) || !(t >= 0.0) || t > 0.5) return qnan;
  const std::vector<double> v = sorted(x, n);
  if (t >= 0.5) return q7s(v, 0.5);
  const int k = static_cast<int>(std::floor(n * t));
  double s = 0.0;
  for (int i = k; i < n - k; ++i) s += v[static_cast<std::size_t>(i)];
  return s / (n - 2.0 * k);
} CLAB_CATCH(kutil::nan())

double clab_wins(const double *x, int n, double t) CLAB_TRY {
  if (n < 1 || hasnan(x, n) || !(t >= 0.0) || t > 0.5) return qnan;
  const std::vector<double> v = sorted(x, n);
  const int k = static_cast<int>(std::floor(n * t));
  if (2 * k >= n) return q7s(v, 0.5);
  const double lo = v[static_cast<std::size_t>(k)], hi = v[static_cast<std::size_t>(n - k - 1)];
  double s = 0.0;
  for (double z : v) s += std::min(std::max(z, lo), hi);
  return s / n;
} CLAB_CATCH(kutil::nan())

void clab_wtd(const double *x, const double *w, int n, double *out) CLAB_TRY {
  out[0] = out[1] = qnan;
  if (n < 1 || hasnan(x, n) || hasnan(w, n)) return;
  double sw = 0.0, s = 0.0, sw2 = 0.0;
  for (int i = 0; i < n; ++i) {
    sw += w[i];
    s += w[i] * x[i];
    sw2 += w[i] * w[i];
  }
  if (sw == 0.0) return;
  const double m = s / sw;
  out[0] = m;
  if (n < 2 || sw <= 0.0 || sw * sw == sw2) return;
  double q = 0.0;
  for (int i = 0; i < n; ++i) {
    const double d = x[i] - m;
    q += w[i] * d * d;
  }
  out[1] = q / (sw - sw2 / sw);
} CLAB_CATCHV

double clab_spear(const double *x, const double *y, int n) CLAB_TRY {
  if (n < 2 || hasnan(x, n) || hasnan(y, n)) return qnan;
  std::vector<double> rx(static_cast<std::size_t>(n)), ry(static_cast<std::size_t>(n));
  midr(x, n, rx.data());
  midr(y, n, ry.data());
  double sx = 0.0, sy = 0.0, sxx = 0.0, syy = 0.0, sxy = 0.0;
  for (int i = 0; i < n; ++i) {
    const double a = rx[static_cast<std::size_t>(i)], b = ry[static_cast<std::size_t>(i)];
    sx += a;
    sy += b;
    sxx += a * a;
    syy += b * b;
    sxy += a * b;
  }
  const double num = n * sxy - sx * sy, den = (n * sxx - sx * sx) * (n * syy - sy * sy);
  return den <= 0.0 ? qnan : num / std::sqrt(den);
} CLAB_CATCH(kutil::nan())

void clab_midrk(const double *x, int n, double *out) CLAB_TRY { midr(x, n, out); } CLAB_CATCHV

void clab_covm(const double *x, int n, int p, double *out) CLAB_TRY {
  const std::size_t un = static_cast<std::size_t>(n), up = static_cast<std::size_t>(p);
  if (n < 2) {
    for (std::size_t k = 0; k < up * up; ++k) out[k] = qnan;
    return;
  }
  std::vector<double> mu(up, 0.0);
  for (std::size_t j = 0; j < up; ++j) {
    double s = 0.0;
    for (std::size_t i = 0; i < un; ++i) s += x[j * un + i];
    mu[j] = s / n;
  }
  for (std::size_t j = 0; j < up; ++j) {
    for (std::size_t k = j; k < up; ++k) {
      double s = 0.0;
      for (std::size_t i = 0; i < un; ++i) s += (x[j * un + i] - mu[j]) * (x[k * un + i] - mu[k]);
      out[k * up + j] = out[j * up + k] = s / (n - 1.0);
    }
  }
} CLAB_CATCHV

int clab_bmean(const double *x, int n, int b, double seed, double *out) CLAB_TRY {
  if (n < 1) {
    for (int k = 0; k < b; ++k) out[k] = qnan;
    return CLAB_OK;
  }
  std::mt19937_64 g(static_cast<std::uint64_t>(seed));
  const std::uint64_t un = static_cast<std::uint64_t>(n);
  for (int k = 0; k < b; ++k) {
    long double s = 0.0L;
    for (int i = 0; i < n; ++i) s += x[bdraw(g, un)];
    out[k] = static_cast<double>(s / static_cast<long double>(n));
    if (kutil::stop(n)) return CLAB_EINTR;
  }
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

void clab_ipw(const double *t, const double *e, int n, double lo, double hi, double *out) CLAB_TRY {
  for (int i = 0; i < n; ++i) {
    const double p = std::min(std::max(e[i], lo), hi);
    out[i] = t[i] == 1.0 ? 1.0 / p : 1.0 / (1.0 - p);
  }
} CLAB_CATCHV

void clab_unif(double seed, int n, double *out) CLAB_TRY {
  std::uint64_t s = static_cast<std::uint64_t>(seed);
  for (int i = 0; i < n; ++i) {
    s += 0x9E3779B97F4A7C15ULL;
    std::uint64_t z = s;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    z ^= z >> 31;
    out[i] = static_cast<double>(z >> 11) * 0x1.0p-53;
  }
} CLAB_CATCHV

}  // extern "C"

extern "C" int clab_shrs(const double *cnt, int k, double tot, double conf, double *share, double *lo, double *hi) CLAB_TRY {
  if (!(conf > 0 && conf < 1) || !(tot > 0)) return CLAB_EDOM;
  const double z = clab_qnorm(1 - (1 - conf) / 2), z2 = z * z;
  for (int i = 0; i < k; ++i) {
    const double c = cnt[i];
    if (c < 0 || c > tot) return CLAB_EDOM;
    const double p = c / tot, den = 1 + z2 / tot;
    const double cen = (p + z2 / (2 * tot)) / den;
    const double half = (z / den) * std::sqrt(p * (1 - p) / tot + z2 / (4 * (tot * tot)));
    share[i] = 100 * c / tot;
    lo[i] = 100 * (c == 0 ? 0 : std::max(0.0, cen - half));
    hi[i] = 100 * (c == tot ? 1 : std::min(1.0, cen + half));
  }
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

extern "C" void clab_ysum(const double *v, const double *chg, int n, double lag, double *out) CLAB_TRY {
  // out: first, last, totalpct, cagrpct, periods, up, down, flat over the
  // rows in period order; v NaN where missing
  int first = -1, last = -1, np = 0;
  for (int i = 0; i < n; ++i) {
    if (std::isnan(v[i])) continue;
    if (first < 0) first = i;
    last = i;
    ++np;
  }
  double up = 0, dn = 0, fl = 0;
  for (int i = 0; i < n; ++i) {
    up += chg[i] > 0;
    dn += chg[i] < 0;
    fl += chg[i] == 0;
  }
  out[5] = up;
  out[6] = dn;
  out[7] = fl;
  out[4] = np;
  if (first < 0) {
    out[0] = out[1] = out[2] = out[3] = qnan;
    return;
  }
  const double a = v[first], b = v[last], steps = (np - 1) / lag;
  out[0] = a;
  out[1] = b;
  out[2] = a > 0 ? 100 * (b / a - 1) : qnan;
  out[3] = a > 0 && steps > 0 ? 100 * (std::pow(b / a, 1 / steps) - 1) : qnan;
} CLAB_CATCHV
