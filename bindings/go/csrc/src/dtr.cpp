// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, deterrence theory: the offending response (clab_detrsp),
// identification of a sanction design (clab_dtrid), the inspection games
// (clab_inspg, clab_insphs) and collider bias of arrest data
// (clab_colarr).
#include <algorithm>
#include <cmath>

#include "crimlab.h"
#include "kutil.h"

namespace {

// The rank of the columns of the n by 4 matrix X (column-major) flagged in
// use, by modified Gram-Schmidt into q (n by 4); r (4 by 4) the triangle.
int mgs(const double *X, int n, const int *use, double *q, double *r, int *keep) {
  int rk = 0;
  for (int j = 0; j < 4; ++j) {
    keep[j] = 0;
    if (!use[j]) continue;
    double *v = q + j * n;
    double n0 = 0;
    for (int i = 0; i < n; ++i) {
      v[i] = X[j * n + i];
      n0 += v[i] * v[i];
    }
    for (int k = 0; k < j; ++k) {
      r[j * 4 + k] = 0;
      if (!keep[k]) continue;
      double d = 0;
      for (int i = 0; i < n; ++i) d += q[k * n + i] * v[i];
      r[j * 4 + k] = d;
      for (int i = 0; i < n; ++i) v[i] -= d * q[k * n + i];
    }
    double nv = 0;
    for (int i = 0; i < n; ++i) nv += v[i] * v[i];
    if (!(n0 > 0) || std::sqrt(nv) < 1e-7 * std::sqrt(n0)) continue;
    nv = std::sqrt(nv);
    r[j * 4 + j] = nv;
    for (int i = 0; i < n; ++i) v[i] /= nv;
    keep[j] = 1;
    ++rk;
  }
  return rk;
}

bool flat(const double *x, int n) {
  for (int i = 1; i < n; ++i) {
    if (x[i] != x[0]) return false;
  }
  return true;
}

}  // namespace

int clab_detrsp(const double *x, const double *b, const double *s, int n, const double *p,
                int np, double *xo, double *v) {
  if (n < 1 || np < 1) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(x[i]) || !std::isfinite(b[i]) || !std::isfinite(s[i])) return CLAB_EDOM;
    if (i > 0 && (x[i] <= x[i - 1] || s[i] <= s[i - 1])) return CLAB_EDOM;
  }
  for (int k = 0; k < np; ++k) {
    if (!kutil::nonneg(p[k])) return CLAB_EDOM;
    double m = -INFINITY;
    for (int i = 0; i < n; ++i) m = std::max(m, b[i] - p[k] * s[i]);
    for (int i = 0; i < n; ++i) {
      if (b[i] - p[k] * s[i] >= m - 1e-12) xo[k] = x[i];
    }
    v[k] = m;
  }
  return CLAB_OK;
}

int clab_dtrid(const double *p, const double *s, const double *c, int n, int lg, const double *y,
               int *rank, int *ident, int *cst, double *beta, double *work) {
  if (n < 1) return CLAB_EARG;
  double *X = work, *q = work + 4 * n, r[16];
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(p[i]) || !std::isfinite(s[i]) || !std::isfinite(c[i])) return CLAB_EDOM;
    if (lg && (!(p[i] > 0) || !(s[i] > 0))) return CLAB_EDOM;
    if (y && !std::isfinite(y[i])) return CLAB_EDOM;
    X[i] = 1;
    X[n + i] = lg ? std::log(p[i]) : p[i];
    X[2 * n + i] = lg ? std::log(s[i]) : s[i];
    X[3 * n + i] = c[i];
  }
  int use[4] = {1, 1, 1, 1}, keep[4];
  for (int j = 1; j < 4; ++j) cst[j - 1] = flat(X + j * n, n);
  for (int j = 1; j < 4; ++j) {
    use[j] = 0;
    ident[j - 1] = mgs(X, n, use, q, r, keep);  // the rank without column j
    use[j] = 1;
  }
  *rank = mgs(X, n, use, q, r, keep);
  for (int j = 0; j < 3; ++j) ident[j] = *rank == 4 || ident[j] < *rank;
  if (!y) return CLAB_OK;
  if (*rank < 4) return CLAB_EDOM;
  // Q'y by projecting out each column in turn, then R beta = Q'y
  double z[4];
  for (int i = 0; i < n; ++i) X[i] = y[i];
  for (int k = 0; k < 4; ++k) {
    double d = 0;
    for (int i = 0; i < n; ++i) d += q[k * n + i] * X[i];
    z[k] = d;
    for (int i = 0; i < n; ++i) X[i] -= d * q[k * n + i];
  }
  for (int j = 3; j >= 0; --j) {
    double v = z[j];
    for (int l = j + 1; l < 4; ++l) v -= r[l * 4 + j] * beta[l];
    beta[j] = v / r[j * 4 + j];
  }
  return CLAB_OK;
}

int clab_inspg(double a, const double *x0, int steps, double dt, double *traj, double *v,
               double *f0, double *dv) {
  if (steps < 0) return CLAB_EARG;
  if (!std::isfinite(a) || a == 0 || !kutil::pos(dt)) return CLAB_EDOM;
  if (!kutil::nonneg(x0[0]) || !kutil::nonneg(x0[1]) || !kutil::nonneg(x0[2]) ||
      std::fabs(x0[0] + x0[1] + x0[2] - 1) > 1e-10) {
    return CLAB_EDOM;
  }
  auto field = [a](const double *x, double *f) {
    f[0] = a * x[0] * (x[1] - x[2]);
    f[1] = a * x[1] * (x[2] - x[0]);
    f[2] = a * x[2] * (x[0] - x[1]);
  };
  const int m = steps + 1;
  double x[3] = {x0[0], x0[1], x0[2]};
  field(x, f0);
  *dv = x[1] * x[2] * f0[0] + x[0] * x[2] * f0[1] + x[0] * x[1] * f0[2];
  for (int k = 0; k < m; ++k) {
    if (k > 0) {
      double k1[3], k2[3], k3[3], k4[3], w[3];
      field(x, k1);
      for (int i = 0; i < 3; ++i) w[i] = x[i] + dt / 2 * k1[i];
      field(w, k2);
      for (int i = 0; i < 3; ++i) w[i] = x[i] + dt / 2 * k2[i];
      field(w, k3);
      for (int i = 0; i < 3; ++i) w[i] = x[i] + dt * k3[i];
      field(w, k4);
      for (int i = 0; i < 3; ++i) x[i] += dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }
    for (int i = 0; i < 3; ++i) traj[i * m + k] = x[i];
    v[k] = x[0] * x[1] * x[2];
  }
  return CLAB_OK;
}

int clab_insphs(const double *T, int nt, const double *g, int ng, int steps, double rate,
                const double *x0, double *out) {
  if (nt < 1 || ng < 1 || steps < 0) return CLAB_EARG;
  if (!kutil::pos(rate)) return CLAB_EDOM;
  const double s0 = x0[0] + x0[1] + x0[2];
  if (!kutil::nonneg(x0[0]) || !kutil::nonneg(x0[1]) || !kutil::nonneg(x0[2]) || !(s0 > 0)) {
    return CLAB_EDOM;
  }
  for (int i = 0; i < nt; ++i) {
    if (!std::isfinite(T[i])) return CLAB_EDOM;
  }
  for (int j = 0; j < ng; ++j) {
    if (!std::isfinite(g[j])) return CLAB_EDOM;
  }
  const long long cells = static_cast<long long>(nt) * ng;
  for (int j = 0; j < ng; ++j) {
    for (int i = 0; i < nt; ++i) {
      const double t = T[i], c = 1 - g[j];
      // payoffs of C, P, O against C, P, O (rows: the strategy played)
      const double A[9] = {1, 0, 1, t, 0, t - 1, c, c, c};
      double x[3] = {x0[0] / s0, x0[1] / s0, x0[2] / s0};
      for (int k = 0; k < steps; ++k) {
        double f[3], phi = 0, sum = 0;
        for (int r = 0; r < 3; ++r) {
          f[r] = A[3 * r] * x[0] + A[3 * r + 1] * x[1] + A[3 * r + 2] * x[2];
          phi += x[r] * f[r];
        }
        for (int r = 0; r < 3; ++r) {
          x[r] = std::max(0.0, x[r] + rate * x[r] * (f[r] - phi));
          sum += x[r];
        }
        // the step keeps the sum at 1 before the clamp, so sum >= 1
        for (int r = 0; r < 3; ++r) x[r] /= sum;
      }
      const long long at = static_cast<long long>(j) * nt + i;
      for (int r = 0; r < 3; ++r) out[r * cells + at] = x[r];
    }
  }
  return CLAB_OK;
}

int clab_colarr(const double *a, const double *e, const double *bg, int n, double *cells,
                double *orr) {
  if (n < 1) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (!(a[i] > 0 && a[i] < 1) || !(e[i] > 0 && e[i] < 1) || !(bg[i] >= 0 && bg[i] <= 1)) {
      return CLAB_EDOM;
    }
    const double c[4] = {a[i] * e[i], a[i] * (1 - e[i]), (1 - a[i]) * e[i],
                         (1 - a[i]) * (1 - e[i]) * bg[i]};
    const double s = c[0] + c[1] + c[2] + c[3];
    for (int k = 0; k < 4; ++k) cells[k * n + i] = c[k] / s;
    orr[i] = c[0] * c[3] / (c[1] * c[2]);
  }
  return CLAB_OK;
}
