// SPDX-License-Identifier: AGPL-3.0-or-later
//
// libcrimlab, corrections: group-based multivariate trajectories (several
// continuous indicators per unit and period, multivariate normal within
// each group), fitted by EM with closed-form M-steps.
#include <cstdint>

#include "crimlab.h"
#include "kutil.h"

namespace {

struct mdat {
  const double *y, *x;
  int n, nt, r, p;
  double yv(int i, int t, int j) const {
    return y[i + static_cast<std::int64_t>(n) * (t + static_cast<std::int64_t>(nt) * j)];
  }
  double xv(int t, int j) const { return x[t + static_cast<std::int64_t>(nt) * j]; }
  // a period is used only when all r indicators are present
  bool row(int i, int t) const {
    for (int j = 0; j < r; ++j)
      if (std::isnan(yv(i, t, j))) return false;
    return true;
  }
};

// Cholesky of the r x r matrix s into l (lower); returns log det or NaN
// when s is not positive definite
double chol(const double *s, int r, double *l) {
  double ld = 0.0;
  for (int j = 0; j < r; ++j) {
    double d = s[j + r * j];
    for (int k = 0; k < j; ++k) d -= l[j + r * k] * l[j + r * k];
    if (!(d > 0.0) || !std::isfinite(d)) return kutil::nan();
    d = std::sqrt(d);
    l[j + r * j] = d;
    ld += 2.0 * std::log(d);
    for (int i = j + 1; i < r; ++i) {
      double e = s[i + r * j];
      for (int k = 0; k < j; ++k) e -= l[i + r * k] * l[j + r * k];
      l[i + r * j] = e / d;
    }
  }
  return ld;
}

// E-step: w (n x ng), returns the log-likelihood, or NaN when a covariance
// is not positive definite. work: r^2 + r.
double estep(const mdat &D, const double *beta, const double *sig,
             const double *prob, int ng, double *w, double *work) {
  const int r = D.r, p = D.p;
  double *l = work, *v = work + r * r;
  const double l2pi = std::log(2.0 * 3.141592653589793);
  for (int k = 0; k < ng; ++k) {
    const double *b = beta + static_cast<std::int64_t>(k) * p * r;
    const double ld = chol(sig + static_cast<std::int64_t>(k) * r * r, r, l);
    if (std::isnan(ld)) return kutil::nan();
    const double lp = std::log(prob[k]);
    for (int i = 0; i < D.n; ++i) {
      kutil::ksum s;
      for (int t = 0; t < D.nt; ++t) {
        if (!D.row(i, t)) continue;
        // residual y - B' x_t, then L v = residual
        double q = 0.0;
        for (int j = 0; j < r; ++j) {
          double mu = 0.0;
          for (int a = 0; a < p; ++a) mu += D.xv(t, a) * b[a + p * j];
          double e = D.yv(i, t, j) - mu;
          for (int c = 0; c < j; ++c) e -= l[j + r * c] * v[c];
          v[j] = e / l[j + r * j];
          q += v[j] * v[j];
        }
        s.add(-0.5 * (r * l2pi + ld + q));
      }
      w[i + static_cast<std::int64_t>(D.n) * k] = lp + s.val();
    }
  }
  kutil::ksum ll;
  for (int i = 0; i < D.n; ++i) {
    double mx = -HUGE_VAL;
    for (int k = 0; k < ng; ++k) {
      const double v2 = w[i + static_cast<std::int64_t>(D.n) * k];
      if (v2 > mx) mx = v2;
    }
    double s = 0.0;
    for (int k = 0; k < ng; ++k) {
      double &v2 = w[i + static_cast<std::int64_t>(D.n) * k];
      v2 = std::exp(v2 - mx);
      s += v2;
    }
    for (int k = 0; k < ng; ++k) w[i + static_cast<std::int64_t>(D.n) * k] /= s;
    ll.add(mx + std::log(s));
  }
  return ll.val();
}

// M-step for group k: weighted least squares (the design is the same for
// every indicator, so the multivariate fit is one p x p solve with r right
// sides) and the weighted residual covariance. False when the weighted
// design is singular or a residual variance is zero to rounding. work:
// p^2 + p r + r^2 + r.
bool mstep(const mdat &D, const double *wk, double *b, double *s, double *work) {
  const int p = D.p, r = D.r, nt = D.nt;
  double *a = work, *bn = a + p * p, *sc = bn + p * r;
  for (int j = 0; j < p * p; ++j) a[j] = 0.0;
  for (int j = 0; j < p * r; ++j) bn[j] = 0.0;
  for (int t = 0; t < nt; ++t) {
    double wt = 0.0;
    for (int j = 0; j < r; ++j) sc[j] = 0.0;  // sum_i w y_itj
    for (int i = 0; i < D.n; ++i) {
      if (!(wk[i] > 0.0) || !D.row(i, t)) continue;
      wt += wk[i];
      for (int j = 0; j < r; ++j) sc[j] += wk[i] * D.yv(i, t, j);
    }
    for (int c = 0; c < p; ++c) {
      for (int e = 0; e < p; ++e) a[c + p * e] += wt * D.xv(t, c) * D.xv(t, e);
      for (int j = 0; j < r; ++j) bn[c + p * j] += D.xv(t, c) * sc[j];
    }
  }
  // Cholesky of a in place (lower), then r triangular solves
  for (int j = 0; j < p; ++j) {
    double d = a[j + p * j];
    for (int k = 0; k < j; ++k) d -= a[j + p * k] * a[j + p * k];
    if (!(d > 0.0) || !std::isfinite(d)) return false;
    d = std::sqrt(d);
    a[j + p * j] = d;
    for (int i = j + 1; i < p; ++i) {
      double e = a[i + p * j];
      for (int k = 0; k < j; ++k) e -= a[i + p * k] * a[j + p * k];
      a[i + p * j] = e / d;
    }
  }
  for (int j = 0; j < r; ++j) {
    double *c = bn + p * j;
    for (int i = 0; i < p; ++i) {
      double e = c[i];
      for (int k = 0; k < i; ++k) e -= a[i + p * k] * c[k];
      c[i] = e / a[i + p * i];
    }
    for (int i = p - 1; i >= 0; --i) {
      double e = c[i];
      for (int k = i + 1; k < p; ++k) e -= a[k + p * i] * c[k];
      c[i] = e / a[i + p * i];
    }
    for (int i = 0; i < p; ++i) b[i + p * j] = c[i];
  }
  // residual covariance, weighted by w_ik over the used periods
  double *m2 = sc + r * r;  // r: weighted mean squares of the indicators
  for (int j = 0; j < r * r; ++j) sc[j] = 0.0;
  for (int j = 0; j < r; ++j) m2[j] = 0.0;
  double tw = 0.0;
  for (int i = 0; i < D.n; ++i) {
    if (!(wk[i] > 0.0)) continue;
    for (int t = 0; t < nt; ++t) {
      if (!D.row(i, t)) continue;
      tw += wk[i];
      for (int j = 0; j < r; ++j) {
        double mj = 0.0;
        for (int c = 0; c < p; ++c) mj += D.xv(t, c) * b[c + p * j];
        const double ej = D.yv(i, t, j) - mj;
        m2[j] += wk[i] * D.yv(i, t, j) * D.yv(i, t, j);
        for (int l = 0; l <= j; ++l) {
          double ml = 0.0;
          for (int c = 0; c < p; ++c) ml += D.xv(t, c) * b[c + p * l];
          sc[j + r * l] += wk[i] * ej * (D.yv(i, t, l) - ml);
        }
      }
    }
  }
  for (int j = 0; j < r; ++j)
    for (int l = 0; l <= j; ++l) s[j + r * l] = s[l + r * j] = sc[j + r * l] / tw;
  // a residual variance below the rounding level of its indicator (sd under
  // 64 ulps of the root mean square) is zero in all but rounding: refused,
  // so the outcome does not depend on how the arithmetic was contracted
  const double e64 = 64.0 * std::numeric_limits<double>::epsilon();
  for (int j = 0; j < r; ++j)
    if (!(s[j + r * j] > e64 * e64 * (m2[j] / tw))) return false;
  return true;
}

int check(const mdat &D, int ng) {
  if (D.n < 1 || D.nt < 1 || D.r < 1 || D.p < 1 || ng < 1) return CLAB_EARG;
  for (int t = 0; t < D.nt; ++t)
    for (int j = 0; j < D.p; ++j)
      if (!std::isfinite(D.xv(t, j))) return CLAB_EDOM;
  for (std::int64_t k = 0; k < static_cast<std::int64_t>(D.n) * D.nt * D.r; ++k)
    if (std::isinf(D.y[k])) return CLAB_EDOM;
  return CLAB_OK;
}

}  // namespace

extern "C" {

int clab_mvtll(const double *y, int n, int nt, int r, const double *x, int p,
               int ng, const double *beta, const double *sig,
               const double *prob, double *w, double *ll, double *work) {
  const mdat D{y, x, n, nt, r, p};
  const int rc = check(D, ng);
  if (rc != CLAB_OK) return rc;
  double s = 0.0;
  for (int k = 0; k < ng; ++k) {
    if (!(prob[k] > 0.0)) return CLAB_EDOM;
    s += prob[k];
  }
  if (std::fabs(s - 1.0) > 1e-9) return CLAB_EDOM;
  const double v = estep(D, beta, sig, prob, ng, w, work);
  if (std::isnan(v)) return CLAB_EDOM;
  *ll = v;
  return CLAB_OK;
}

int clab_mvtem(const double *y, int n, int nt, int r, const double *x, int p,
               int ng, double *beta, double *sig, double *prob, double *w,
               double *ll, int maxit, double tol, int *iter, double *work) {
  if (maxit < 1 || !(tol > 0.0)) return CLAB_EARG;
  int rc = clab_mvtll(y, n, nt, r, x, p, ng, beta, sig, prob, w, ll, work);
  if (rc != CLAB_OK) return rc;
  const mdat D{y, x, n, nt, r, p};
  double cur = *ll;
  int it = 0;
  for (; it < maxit; ++it) {
    for (int k = 0; k < ng; ++k) {
      double s = 0.0;
      for (int i = 0; i < n; ++i) s += w[i + static_cast<std::int64_t>(n) * k];
      prob[k] = s / n > 1e-300 ? s / n : 1e-300;
      if (!(s > 0.0)) continue;  // an empty group keeps its trajectory
      if (!mstep(D, w + static_cast<std::int64_t>(n) * k, beta + static_cast<std::int64_t>(k) * p * r,
                 sig + static_cast<std::int64_t>(k) * r * r, work))
        return CLAB_EDOM;
    }
    const double nl = estep(D, beta, sig, prob, ng, w, work);
    if (std::isnan(nl)) return CLAB_EDOM;
    const double gain = nl - cur;
    cur = nl;
    if (gain <= tol * (1.0 + std::fabs(nl))) {
      ++it;
      break;
    }
  }
  *ll = cur;
  *iter = it;
  return CLAB_OK;
}

}  // extern "C"
