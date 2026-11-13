// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, taphonomy: log-ratio transforms (clab_clrilr), likelihood
// ratios of measured evidence (clab_preslr), the decay chain
// (clab_tpdmc, clab_tpdsim) and a conjugate Bayesian preservation model
// (clab_tpbhm).
#include <algorithm>
#include <cmath>

#include "crimlab.h"
#include "kutil.h"

namespace {

// In-place Cholesky factor (lower, column-major) of the p by p matrix a,
// then its inverse in a; false when a is not positive definite.
bool cinv(double *a, int p) {
  for (int j = 0; j < p; ++j) {
    double d = a[j * p + j];
    for (int k = 0; k < j; ++k) d -= a[k * p + j] * a[k * p + j];
    if (!(d > 1e-12 * std::fabs(a[j * p + j]))) return false;
    d = std::sqrt(d);
    a[j * p + j] = d;
    for (int i = j + 1; i < p; ++i) {
      double s = a[j * p + i];
      for (int k = 0; k < j; ++k) s -= a[k * p + i] * a[k * p + j];
      a[j * p + i] = s / d;
    }
  }
  // L^-1 in the lower triangle, then (L^-1)' L^-1
  for (int j = 0; j < p; ++j) {
    a[j * p + j] = 1 / a[j * p + j];
    for (int i = j + 1; i < p; ++i) {
      double s = 0;
      for (int k = j; k < i; ++k) s -= a[k * p + i] * a[j * p + k];
      a[j * p + i] = s / a[i * p + i];
    }
  }
  for (int j = 0; j < p; ++j) {
    for (int i = j; i < p; ++i) {
      double s = 0;
      for (int k = i; k < p; ++k) s += a[i * p + k] * a[j * p + k];
      a[j * p + i] = s;  // lower: row i >= j holds (i, j)
    }
  }
  for (int j = 0; j < p; ++j) {
    for (int i = 0; i < j; ++i) a[j * p + i] = a[i * p + j];
  }
  return true;
}

// x'v for column j of the design (1, x)
double xcol(const double *x, int n, int j, int i) { return j == 0 ? 1.0 : x[(j - 1) * n + i]; }

}  // namespace

int clab_clrilr(const double *x, int n, int d, double ps, int type, double *out) {
  if (n < 1 || d < 2 || type < 0 || type > 1) return CLAB_EARG;
  if (!kutil::nonneg(ps)) return CLAB_EDOM;
  for (int i = 0; i < n; ++i) {
    kutil::ksum s;
    for (int j = 0; j < d; ++j) {
      const double v = x[j * n + i];
      if (!kutil::nonneg(v)) return CLAB_EDOM;
      s.add(v + ps);
    }
    const double tot = s.val();
    double lm = 0;
    for (int j = 0; j < d; ++j) {
      const double v = (x[j * n + i] + ps) / tot;
      if (!(v > 0)) return CLAB_EDOM;
      const double l = std::log(v);
      if (type == 0) out[j * n + i] = l;
      lm += l;
      if (type == 1 && j > 0) {
        out[(j - 1) * n + i] = std::sqrt(static_cast<double>(j) / (j + 1)) * ((lm - l) / j - l);
      }
    }
    if (type == 0) {
      for (int j = 0; j < d; ++j) out[j * n + i] -= lm / d;
    }
  }
  return CLAB_OK;
}

int clab_preslr(const double *x, int n, const double *m1, const double *s1, int k1,
                const double *m2, const double *s2, int k2, double *out, int *band, int *side) {
  if (n < 1 || (k1 != 1 && k1 != n) || (k2 != 1 && k2 != n)) return CLAB_EARG;
  const double c = 0.91893853320467274178;  // log(2 pi) / 2
  kutil::ksum l1, l2;
  for (int i = 0; i < n; ++i) {
    const double a = m1[k1 == 1 ? 0 : i], b = s1[k1 == 1 ? 0 : i];
    const double e = m2[k2 == 1 ? 0 : i], f = s2[k2 == 1 ? 0 : i];
    if (!std::isfinite(x[i]) || !std::isfinite(a) || !std::isfinite(e) || !kutil::pos(b) ||
        !kutil::pos(f)) {
      return CLAB_EDOM;
    }
    l1.add(-c - std::log(b) - 0.5 * ((x[i] - a) / b) * ((x[i] - a) / b));
    l2.add(-c - std::log(f) - 0.5 * ((x[i] - e) / f) * ((x[i] - e) / f));
  }
  out[0] = l1.val();
  out[1] = l2.val();
  out[2] = out[0] - out[1];
  out[3] = std::exp(out[2]);
  out[4] = out[2] / std::log(10.0);
  const double a = std::fabs(out[4]);
  static const double cut[6] = {0, 1, 2, 3, 4, 6};
  int b = 6;
  while (b > 0 && a <= cut[b - 1]) --b;
  *band = b;
  *side = out[2] > 0 ? 1 : out[2] < 0 ? 2 : 0;
  return CLAB_OK;
}

int clab_tpdmc(double pres, double dr, double mr, int k, double *P, double *N, double *B,
               double *t) {
  if (k < 1) return CLAB_EARG;
  if (!(pres >= 0 && pres <= 1) || !(dr > 0 && dr <= 1) || !(mr >= 0 && mr <= 1)) {
    return CLAB_EDOM;
  }
  // the two never sum above max(dr, mr) <= 1
  const double prog = dr * (1 - pres), mum = mr * pres, tot = prog + mum;
  if (!(tot > 0)) return CLAB_EDOM;
  const int m = k + 2;
  std::fill(P, P + m * m, 0.0);
  for (int i = 0; i < k; ++i) {
    P[i * m + i] = 1 - tot;
    P[(i + 1) * m + i] += prog;  // the next stage, or skeletal after the last
    P[(k + 1) * m + i] += mum;
  }
  P[k * m + k] = 1;
  P[(k + 1) * m + k + 1] = 1;
  // I - Q is upper bidiagonal (tot on the diagonal, -prog above it): its
  // inverse by back substitution, column by column
  std::fill(N, N + k * k, 0.0);
  for (int j = 0; j < k; ++j) {
    N[j * k + j] = 1 / tot;
    for (int i = j - 1; i >= 0; --i) N[j * k + i] = prog * N[j * k + i + 1] / tot;
  }
  for (int i = 0; i < k; ++i) {
    kutil::ksum s, mm;
    for (int j = 0; j < k; ++j) {
      s.add(N[j * k + i]);
      mm.add(N[j * k + i] * mum);
    }
    t[i] = s.val();
    B[i] = N[(k - 1) * k + i] * prog;
    B[k + i] = mm.val();
  }
  return CLAB_OK;
}

int clab_tpdsim(const double *P, int m, int k, int start, int steps, unsigned long long seed,
                int *path, int *len) {
  if (m < 2 || k < 1 || k >= m || start < 0 || start >= k || steps < 0) return CLAB_EARG;
  kutil::rng r(seed);
  int s = start, n = 1;
  path[0] = s;
  for (int i = 0; i < steps && s < k; ++i) {
    const double u = r.unif();
    double c = 0;
    int nx = m - 1;
    for (int j = 0; j < m; ++j) {
      c += P[j * m + s];
      if (u < c) {
        nx = j;
        break;
      }
    }
    s = nx;
    path[n++] = s;
  }
  *len = n;
  return CLAB_OK;
}

int clab_tpbhm(const double *y, const double *x, int n, int q, const double *m0,
               const double *s0, const int *g, int ng, double *mu, double *sd, double *st,
               double *gm, double *lam, int *nj, double *work) {
  const int p = q + 1;
  if (n < 2 || q < 1 || (g && ng < 2)) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(y[i])) return CLAB_EDOM;
    for (int j = 1; j < p; ++j) {
      if (!std::isfinite(xcol(x, n, j, i))) return CLAB_EDOM;
    }
    if (g && (g[i] < 0 || g[i] >= ng)) return CLAB_EARG;
  }
  for (int j = 0; j < p; ++j) {
    if (!std::isfinite(m0[j]) || !kutil::pos(s0[j])) return CLAB_EDOM;
  }
  double *a = work, *c = work + p * p, *v = c + p * p;
  // least squares: (X'X)^-1 X'y
  for (int j = 0; j < p; ++j) {
    kutil::ksum b;
    for (int i = 0; i < n; ++i) b.add(xcol(x, n, j, i) * y[i]);
    v[j] = b.val();
    for (int l = 0; l < p; ++l) {
      kutil::ksum s;
      for (int i = 0; i < n; ++i) s.add(xcol(x, n, j, i) * xcol(x, n, l, i));
      a[l * p + j] = s.val();
    }
  }
  std::copy(a, a + p * p, c);
  if (n < p || !cinv(c, p)) return CLAB_EDOM;
  kutil::ksum rs, ym, yy;
  for (int i = 0; i < n; ++i) {
    double f = 0;
    for (int j = 0; j < p; ++j) {
      double bj = 0;
      for (int l = 0; l < p; ++l) bj += c[l * p + j] * v[l];
      f += xcol(x, n, j, i) * bj;
    }
    rs.add((y[i] - f) * (y[i] - f));
    ym.add(y[i]);
    yy.add(y[i] * y[i]);
  }
  // An exact fit leaves residuals of rounding only, and their squares do not sum to 0 on every
  // machine (a fused multiply-add leaves about 1e-30 here): residuals below 1e-12 of the size
  // of y are none.
  double s2 = rs.val() > 1e-24 * yy.val() ? rs.val() / std::max(1, n - p) : 0;
  if (!(s2 > 0)) {
    kutil::ksum vy;
    for (int i = 0; i < n; ++i) vy.add((y[i] - ym.val() / n) * (y[i] - ym.val() / n));
    s2 = vy.val() / (n - 1);
    if (!(s2 > 0)) return CLAB_EDOM;
  }
  // posterior: (X'X / s2 + D)^-1 (X'y / s2 + D m0), D = diag(1 / s0^2)
  for (int j = 0; j < p; ++j) {
    for (int l = 0; l < p; ++l) c[l * p + j] = a[l * p + j] / s2;
    c[j * p + j] += 1 / (s0[j] * s0[j]);
    v[j] = v[j] / s2 + m0[j] / (s0[j] * s0[j]);
  }
  cinv(c, p);  // positive definite: X'X / s2 is, and D adds to it
  for (int j = 0; j < p; ++j) {
    double s = 0;
    for (int l = 0; l < p; ++l) s += c[l * p + j] * v[l];
    mu[j] = s;
    sd[j] = std::sqrt(c[j * p + j]);
  }
  st[0] = s2;
  st[1] = kutil::nan();
  if (!g) return CLAB_OK;
  std::fill(gm, gm + ng, 0.0);
  std::fill(nj, nj + ng, 0);
  for (int i = 0; i < n; ++i) {
    double f = 0;
    for (int j = 0; j < p; ++j) f += xcol(x, n, j, i) * mu[j];
    gm[g[i]] += y[i] - f;
    ++nj[g[i]];
  }
  kutil::ksum sg, sw;
  for (int j = 0; j < ng; ++j) {
    if (!nj[j]) return CLAB_EARG;
    gm[j] /= nj[j];
    sg.add(gm[j]);
    sw.add(s2 / nj[j]);
  }
  kutil::ksum vg;
  for (int j = 0; j < ng; ++j) vg.add((gm[j] - sg.val() / ng) * (gm[j] - sg.val() / ng));
  const double tau2 = std::max(0.0, vg.val() / (ng - 1) - sw.val() / ng);
  for (int j = 0; j < ng; ++j) lam[j] = tau2 / (tau2 + s2 / nj[j]);
  st[1] = tau2;
  return CLAB_OK;
}
