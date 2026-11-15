// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, regression discontinuity: the kernel constant and the
// bandwidth of Imbens and Kalyanaraman (2012) for a local-linear fit at a
// cutoff (clab_rdck, clab_ikbw), and the density test of McCrary (2008)
// (clab_mccr). The local fits themselves are clab_ivfe with kernel weights
// (src/apx.cpp). Least squares is that of the policy-evaluation families
// (pev.h).
#include <cmath>
#include <cstddef>
#include <vector>

#include "crimlab.h"
#include "kutil.h"
#include "pev.h"

namespace {

using std::size_t;

// The least squares fit of y on the columns of x (n x p) into b; false
// when the columns are collinear or there are too few rows.
bool lsq(const std::vector<double> &x, const std::vector<double> &y, int p, double *b,
         double *mse) {
  const int n = static_cast<int>(y.size());
  if (n <= p) return false;
  std::vector<double> xi(static_cast<size_t>(p) * p), r(static_cast<size_t>(n));
  std::vector<double> wk(static_cast<size_t>(p) * p + p);
  if (pev::ols(x.data(), y.data(), n, p, nullptr, b, xi.data(), r.data(), wk.data()) != CLAB_OK) {
    return false;
  }
  if (mse) {
    double s = 0.0;
    for (int i = 0; i < n; ++i) s += r[static_cast<size_t>(i)] * r[static_cast<size_t>(i)];
    *mse = s / (n - p);
  }
  return true;
}

// Powers 0..p of (v - c) as the columns of a design, for the rows kept.
std::vector<double> poly(const std::vector<double> &v, double c, int p) {
  const size_t n = v.size();
  std::vector<double> x(n * static_cast<size_t>(p + 1));
  for (size_t i = 0; i < n; ++i) {
    double t = 1.0;
    for (int j = 0; j <= p; ++j) {
      x[static_cast<size_t>(j) * n + i] = t;
      t *= v[i] - c;
    }
  }
  return x;
}

}  // namespace

extern "C" {

double clab_rdck(int kern) {
  // the moments on [0, 1] of the kernel (any constant multiple: C_K does
  // not depend on it): nu_j of K, pi_j of K^2
  double nu[4], pi[3];
  for (int j = 0; j < 4; ++j) {
    const double a = j + 1.0;
    if (kern == 0) {
      nu[j] = 1 / (a * (a + 1));
      if (j < 3) pi[j] = 2 / (a * (a + 1) * (a + 2));
    } else if (kern == 1) {
      nu[j] = 1 / a;
      if (j < 3) pi[j] = 1 / a;
    } else if (kern == 2) {
      nu[j] = 1 / a - 1 / (a + 2);
      if (j < 3) pi[j] = 1 / a - 2 / (a + 2) + 1 / (a + 4);
    } else {
      return kutil::nan();
    }
  }
  // Imbens and Kalyanaraman (2012, Lemma 3.1 and p. 937-938): C1, C2 and
  // C_K = (C2 / (4 C1))^(1/5)
  const double d = nu[2] * nu[0] - nu[1] * nu[1];
  const double q = (nu[2] * nu[2] - nu[1] * nu[3]) / d;
  const double c1 = q * q / 4;
  const double c2 = (nu[2] * nu[2] * pi[0] - 2 * nu[1] * nu[2] * pi[1] + nu[1] * nu[1] * pi[2]) / (d * d);
  return std::pow(c2 / (4 * c1), 0.2);
}

int clab_ikbw(const double *y, const double *x, int n, double c, int kern, double *h, double *st) {
  if (n < 10 || kern < 0 || kern > 2 || !std::isfinite(c)) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(y[i]) || !std::isfinite(x[i])) return CLAB_EDOM;
  }
  const size_t sn = static_cast<size_t>(n);
  // step 1: the pilot bandwidth h1 = 1.84 S_X N^(-1/5), the density and
  // the variances on either side (their equations 11-13, p. 941-942)
  double mx = 0.0;
  for (size_t i = 0; i < sn; ++i) mx += x[i] / n;
  double sx = 0.0;
  for (size_t i = 0; i < sn; ++i) sx += (x[i] - mx) * (x[i] - mx);
  sx = std::sqrt(sx / (n - 1));
  const double h1 = 1.84 * sx * std::pow(n, -0.2);
  double nl = 0, nr = 0, yl = 0, yr = 0;
  for (size_t i = 0; i < sn; ++i) {
    if (x[i] >= c - h1 && x[i] < c) {
      ++nl;
      yl += y[i];
    } else if (x[i] >= c && x[i] <= c + h1) {
      ++nr;
      yr += y[i];
    }
  }
  if (nl < 2 || nr < 2) return CLAB_EDOM;
  yl /= nl;
  yr /= nr;
  double vl = 0, vr = 0;
  for (size_t i = 0; i < sn; ++i) {
    if (x[i] >= c - h1 && x[i] < c) {
      vl += (y[i] - yl) * (y[i] - yl);
    } else if (x[i] >= c && x[i] <= c + h1) {
      vr += (y[i] - yr) * (y[i] - yr);
    }
  }
  vl /= nl - 1;
  vr /= nr - 1;
  const double f = (nl + nr) / (2.0 * n * h1);
  // step 2: the third derivative from the global cubic with a jump
  // (equation 14), the pilot bandwidths h2 (equation 15) and the second
  // derivatives by local quadratics on either side
  std::vector<double> xv(x, x + n), yv(y, y + n), d5(sn * 5);
  double nlt = 0, nrt = 0;
  for (size_t i = 0; i < sn; ++i) {
    const double t = x[i] - c;
    d5[i] = 1;
    d5[sn + i] = x[i] >= c ? 1 : 0;
    d5[2 * sn + i] = t;
    d5[3 * sn + i] = t * t;
    d5[4 * sn + i] = t * t * t;
    (x[i] >= c ? nrt : nlt) += 1;
  }
  double g[5];
  if (!lsq(d5, yv, 5, g, nullptr)) return CLAB_EDOM;
  const double m3 = 6 * g[4];
  if (!(m3 != 0) || !(f > 0)) return CLAB_EDOM;
  const double h2r = 3.56 * std::pow(vr / (f * m3 * m3), 1.0 / 7) * std::pow(nrt, -1.0 / 7);
  const double h2l = 3.56 * std::pow(vl / (f * m3 * m3), 1.0 / 7) * std::pow(nlt, -1.0 / 7);
  std::vector<double> xr, yr2, xl, yl2;
  for (size_t i = 0; i < sn; ++i) {
    if (x[i] >= c && x[i] <= c + h2r) {
      xr.push_back(x[i]);
      yr2.push_back(y[i]);
    }
    if (x[i] >= c - h2l && x[i] < c) {
      xl.push_back(x[i]);
      yl2.push_back(y[i]);
    }
  }
  double lr[3], ll[3];
  if (!lsq(poly(xr, c, 2), yr2, 3, lr, nullptr) || !lsq(poly(xl, c, 2), yl2, 3, ll, nullptr)) {
    return CLAB_EDOM;
  }
  const double m2r = 2 * lr[2], m2l = 2 * ll[2];
  // step 3: the regularisation terms (equation 16) and the bandwidth
  // (equation 17)
  const double rr = 2160 * vr / (static_cast<double>(xr.size()) * std::pow(h2r, 4));
  const double rl = 2160 * vl / (static_cast<double>(xl.size()) * std::pow(h2l, 4));
  const double ck = clab_rdck(kern);
  const double dm = m2r - m2l;
  *h = ck * std::pow((vl + vr) / (f * (dm * dm + rr + rl)), 0.2) * std::pow(n, -0.2);
  if (!(*h > 0) || !std::isfinite(*h)) return CLAB_EDOM;
  st[0] = h1;
  st[1] = f;
  st[2] = vl;
  st[3] = vr;
  st[4] = m3;
  st[5] = h2l;
  st[6] = h2r;
  st[7] = m2l;
  st[8] = m2r;
  st[9] = ck;
  return CLAB_OK;
}

int clab_mccr(const double *r, int n, double c, double b, double h, double *st) {
  if (n < 10 || !std::isfinite(c) || !(b >= 0) || !(h >= 0)) return CLAB_EARG;
  double mr = 0.0, lo = r[0], hi = r[0];
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(r[i])) return CLAB_EDOM;
    mr += r[i] / n;
    lo = std::fmin(lo, r[i]);
    hi = std::fmax(hi, r[i]);
  }
  if (!(lo < c) || !(hi >= c)) return CLAB_EDOM;
  // the binsize: 2 s n^(-1/2) by default (McCrary 2008, p. 705, step 1)
  if (!(b > 0)) {
    double s = 0.0;
    for (int i = 0; i < n; ++i) s += (r[i] - mr) * (r[i] - mr);
    b = 2 * std::sqrt(s / (n - 1)) / std::sqrt(static_cast<double>(n));
  }
  // the first-step histogram: g(R) = floor((R - c) / b) b + b / 2 + c
  // (equation 2), bins of width b over the support of g, Y_j the share
  // of the rows in bin j over b
  const long jlo = static_cast<long>(std::floor((lo - c) / b));
  const long jhi = static_cast<long>(std::floor((hi - c) / b));
  const long nb = jhi - jlo + 1;
  if (nb < 4 || nb > 10000000L) return CLAB_EDOM;
  std::vector<double> xj(static_cast<size_t>(nb)), yj(static_cast<size_t>(nb), 0.0);
  for (long j = 0; j < nb; ++j) xj[static_cast<size_t>(j)] = (jlo + j) * b + b / 2 + c;
  for (int i = 0; i < n; ++i) {
    const long j = static_cast<long>(std::floor((r[i] - c) / b)) - jlo;
    yj[static_cast<size_t>(j)] += 1.0 / (n * b);
  }
  // the bandwidth, when not given: a quartic in the bin midpoints on each
  // side, k (s^2 (b - a) / sum f''(X_j)^2)^(1/5) with k = 3.348, averaged
  // over the sides (p. 705, step 2)
  if (!(h > 0)) {
    double hs = 0.0;
    for (int side = 0; side < 2; ++side) {
      std::vector<double> xs, ys;
      for (long j = 0; j < nb; ++j) {
        const double v = xj[static_cast<size_t>(j)];
        if ((side == 0 && v < c) || (side == 1 && v > c)) {
          xs.push_back(v);
          ys.push_back(yj[static_cast<size_t>(j)]);
        }
      }
      double a[5], s2 = 0.0;
      if (!lsq(poly(xs, c, 4), ys, 5, a, &s2)) return CLAB_EDOM;
      double f2 = 0.0;
      for (double v : xs) {
        const double t = v - c, d2 = 2 * a[2] + 6 * a[3] * t + 12 * a[4] * t * t;
        f2 += d2 * d2;
      }
      const double len = side == 0 ? c - xs.front() : xs.back() - c;
      if (!(f2 > 0)) return CLAB_EDOM;
      hs += 3.348 * std::pow(s2 * len / f2, 0.2) / 2;
    }
    h = hs;
  }
  // the second step: local linear fits of Y_j on X_j - c with the
  // triangle kernel on each side; theta = ln f+ - ln f- and its standard
  // error sqrt((1 / (n h)) (24 / 5) (1 / f+ + 1 / f-)) (equation 5)
  double fe[2];
  for (int side = 0; side < 2; ++side) {
    double s0 = 0, s1 = 0, s2 = 0, t0 = 0, t1 = 0;
    for (long j = 0; j < nb; ++j) {
      const double v = xj[static_cast<size_t>(j)];
      if ((side == 0 && !(v < c)) || (side == 1 && !(v > c))) continue;
      const double u = v - c, k = std::fmax(0.0, 1 - std::fabs(u / h));
      if (k <= 0) continue;
      s0 += k;
      s1 += k * u;
      s2 += k * u * u;
      t0 += k * yj[static_cast<size_t>(j)];
      t1 += k * u * yj[static_cast<size_t>(j)];
    }
    const double det = s0 * s2 - s1 * s1;
    if (!(det > 0)) return CLAB_EDOM;
    fe[side] = (s2 * t0 - s1 * t1) / det;
  }
  if (!(fe[0] > 0) || !(fe[1] > 0)) return CLAB_EDOM;
  st[0] = std::log(fe[1]) - std::log(fe[0]);
  st[1] = std::sqrt(24.0 / 5.0 / (n * h) * (1 / fe[1] + 1 / fe[0]));
  st[2] = fe[0];
  st[3] = fe[1];
  st[4] = b;
  st[5] = h;
  st[6] = static_cast<double>(nb);
  return CLAB_OK;
}

}  // extern "C"
