// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, series: the Mann-Kendall test with the Theil-Sen slope
// (clab_trndrr), a single step change by a permutation scan and Pettitt's
// test (clab_stpcp), CUSUM surveillance of counts (clab_cusrv) and the
// cross-correlation of two series (clab_xcser). No allocation: the caller
// passes the work arrays.
#include <cmath>
#include <cstddef>

#include "crimlab.h"
#include "kutil.h"

namespace {

double sgn(double a) { return a > 0 ? 1.0 : (a < 0 ? -1.0 : 0.0); }

// Mann-Kendall's S of y in the order given
double mks(const double *y, int n) {
  double s = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) s += sgn(y[j] - y[i]);
  }
  return s;
}

double median(double *v, int n) {
  kutil::hsort(v, n, [](double a, double b) { return a < b; });
  return n % 2 ? v[n / 2] : (v[n / 2 - 1] + v[n / 2]) / 2;
}

double unorm(double z) { return 0.5 * std::erfc(z / std::sqrt(2.0)); }

// the next ordering of v by Heap's algorithm (iterative); c holds n
// counters from 0, i the position; false when every ordering was visited
bool heap(double *v, int n, int *c, int *i) {
  while (*i < n) {
    if (c[*i] < *i) {
      const int j = *i % 2 ? c[*i] : 0;
      const double s = v[j];
      v[j] = v[*i];
      v[*i] = s;
      ++c[*i];
      *i = 1;
      return true;
    }
    c[*i] = 0;
    ++*i;
  }
  return false;
}

// the step scan's largest between-segment sum of squares and where
double tmax(const double *y, int n, int m, int *at) {
  double tot = 0, best = -1, a = 0;
  for (int i = 0; i < n; ++i) tot += y[i];
  for (int k = 1; k < n; ++k) {
    a += y[k - 1];
    if (k < m || n - k < m) continue;
    const double d = a / k - (tot - a) / (n - k), t = static_cast<double>(k) * (n - k) / n * d * d;
    if (t > best) {
      best = t;
      *at = k;
    }
  }
  return best;
}

}  // namespace

extern "C" {

int clab_trndrr(const double *x, const double *y, int n, int alt, double z, double *st,
                double *work, int *iwork) {
  if (n < 3 || alt < 0 || alt > 2 || !(z > 0)) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(x[i]) || !std::isfinite(y[i]) || (i > 0 && !(x[i] > x[i - 1]))) {
      return CLAB_EDOM;
    }
  }
  const double s = mks(y, n), np = n * (n - 1) / 2.0;
  // the variance with ties: sum over groups of tied values
  double *v = work;
  for (int i = 0; i < n; ++i) v[i] = y[i];
  kutil::hsort(v, n, [](double a, double b) { return a < b; });
  double var = n * (n - 1.0) * (2 * n + 5) / 18;
  for (int i = 0; i < n;) {
    int j = i;
    while (j < n && v[j] == v[i]) ++j;
    const double g = j - i;
    var -= g * (g - 1) * (2 * g + 5) / 18;
    i = j;
  }
  double p;
  if (n <= 8) {
    // every ordering of the observed values, ties kept
    for (int i = 0; i < n; ++i) {
      v[i] = y[i];
      iwork[i] = 0;
    }
    double cnt = 0, tot = 0;
    int pos = 1;
    do {
      const double q = mks(v, n);
      cnt += alt == 0 ? std::fabs(q) >= std::fabs(s) : (alt == 1 ? q >= s : q <= s);
      ++tot;
    } while (heap(v, n, iwork, &pos));
    p = cnt / tot;
  } else {
    const double zz = var > 0 ? (s - sgn(s)) / std::sqrt(var) : 0;
    p = var <= 0 ? 1 : (alt == 0 ? 2 * unorm(std::fabs(zz)) : (alt == 1 ? unorm(zz) : unorm(-zz)));
  }
  // Theil-Sen: the median of the pairwise slopes, and Sen's interval
  // between the order statistics the variance of S places at level z
  double *sl = work + n;
  int m = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) sl[m++] = (y[j] - y[i]) / (x[j] - x[i]);
  }
  const double b = median(sl, m);
  const double ca = z * std::sqrt(var), lr = std::floor((m - ca) / 2), hr = std::ceil((m + ca) / 2) + 1;
  for (int i = 0; i < n; ++i) v[i] = y[i] - b * x[i];
  st[0] = s;
  st[1] = var;
  st[2] = s / np;
  st[3] = std::fmin(1.0, p);
  st[4] = b;
  st[5] = median(v, n);
  st[6] = sl[static_cast<int>(std::fmin(m, std::fmax(1, lr))) - 1];
  st[7] = sl[static_cast<int>(std::fmin(m, std::fmax(1, hr))) - 1];
  st[8] = lr <= 1;
  st[9] = hr >= m;
  return CLAB_OK;
}

int clab_stpcp(const double *y, int n, int minseg, int nperm, unsigned int seed, double *st,
               double *work, int *iwork) {
  if (minseg < 1 || n < 2 * minseg || nperm < 1) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(y[i])) return CLAB_EDOM;
  }
  int at = 0;
  const double tobs = tmax(y, n, minseg, &at);
  double a = 0;
  for (int i = 0; i < at; ++i) a += y[i];
  double tot = a;
  for (int i = at; i < n; ++i) tot += y[i];
  st[0] = at;
  st[1] = a / at;
  st[2] = (tot - a) / (n - at);
  st[3] = tobs;
  // the null of the largest statistic over reorderings: all of them when
  // n! is at most nperm, else nperm random ones
  double nf = 1;
  for (int i = 2; i <= n && nf <= nperm; ++i) nf *= i;
  double *v = work;
  for (int i = 0; i < n; ++i) v[i] = y[i];
  const double lim = tobs - 1e-12 * std::fabs(tobs);
  double cnt = 0, tt = 0;
  int k;
  if (nf <= nperm) {
    for (int i = 0; i < n; ++i) iwork[i] = 0;
    int pos = 1;
    do {
      cnt += tmax(v, n, minseg, &k) >= lim;
      ++tt;
    } while (heap(v, n, iwork, &pos));
    st[4] = cnt / tt;
    st[6] = 1;
  } else {
    kutil::rng r(seed);
    for (int b = 0; b < nperm; ++b) {
      for (int i = n - 1; i > 0; --i) {
        const int j = r.below(i + 1);
        const double s = v[i];
        v[i] = v[j];
        v[j] = s;
      }
      cnt += tmax(v, n, minseg, &k) >= lim;
    }
    tt = nperm;
    st[4] = (1 + cnt) / (1 + tt);
    st[6] = 0;
  }
  st[5] = tt;
  // Pettitt (1979): U_t = sum over i <= t < j of sign(y_j - y_i)
  double u = 0, best = -1;
  int pt = 1;
  for (int t = 1; t < n; ++t) {
    for (int j = t; j < n; ++j) u += sgn(y[j] - y[t - 1]);
    for (int i = 0; i < t - 1; ++i) u -= sgn(y[t - 1] - y[i]);
    if (std::fabs(u) > best) {
      best = std::fabs(u);
      pt = t;
    }
  }
  a = 0;
  for (int i = 0; i < pt; ++i) a += y[i];
  st[7] = pt;
  st[8] = best;
  st[9] = std::fmin(1.0, 2 * std::exp(-6 * best * best / (static_cast<double>(n) * n * n + static_cast<double>(n) * n)));
  st[10] = a / pt;
  st[11] = (tot - a) / (n - pt);
  return CLAB_OK;
}

int clab_cusrv(const double *y, const double *m, int n, double k, double h, int trans,
               int reset, double *z, double *s, int *alarm) {
  if (n < 1 || trans < 0 || trans > 3 || !std::isfinite(k) || !kutil::pos(h)) return CLAB_EARG;
  double c = 0;
  for (int i = 0; i < n; ++i) {
    const double v = y[i], e = m[i];
    if (!kutil::nonneg(v) || !kutil::pos(e)) return CLAB_EDOM;
    switch (trans) {
      case 0:
        z[i] = (v - e) / std::sqrt(e);
        break;
      case 1:  // Rossi et al. (1999)
        z[i] = (v - 3 * e + 2 * std::sqrt(v * e)) / (2 * std::sqrt(e));
        break;
      case 2:  // Anscombe
        z[i] = 1.5 * (std::cbrt(v * v) - std::cbrt(e * e)) / std::pow(e, 1.0 / 6);
        break;
      default:
        z[i] = v;
    }
    c = std::fmax(0.0, c + z[i] - k);
    alarm[i] = c >= h;
    s[i] = c;
    if (reset && alarm[i]) c = 0;
  }
  return CLAB_OK;
}

int clab_xcser(const double *x, const double *y, int n, int lag, double *r) {
  if (n < 2 || lag < 0 || lag >= n) return CLAB_EARG;
  double mx = 0, my = 0;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(x[i]) || !std::isfinite(y[i])) return CLAB_EDOM;
    mx += x[i] / n;
    my += y[i] / n;
  }
  double sx = 0, sy = 0;
  for (int i = 0; i < n; ++i) {
    sx += (x[i] - mx) * (x[i] - mx);
    sy += (y[i] - my) * (y[i] - my);
  }
  if (!(sx > 0 && sy > 0)) return CLAB_EDOM;
  const double d = std::sqrt(sx * sy);
  // r[lag + k] estimates cor(x[t + k], y[t]), as stats::ccf()
  for (int k = -lag; k <= lag; ++k) {
    double c = 0;
    for (int t = 0; t < n; ++t) {
      if (t + k >= 0 && t + k < n) c += (x[t + k] - mx) * (y[t] - my);
    }
    r[lag + k] = c / d;
  }
  return CLAB_OK;
}

}  // extern "C"
