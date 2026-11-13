// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, forensics: the two-level normal model for multivariate trace
// evidence (Aitken and Lucy 2004): variance components from a background
// sample, and the likelihood ratio of two items in log space.
#include <cstdint>

#include "crimlab.h"
#include "kutil.h"

namespace {

// M = a U + b C (lower triangle of the p x p column-major inputs), factored
// in place in w (p * p) as M = L L'. With r in, writes *q = r' M^-1 r and
// *ld = log det M; v (p) is scratch. False when M is not positive definite
// to working precision (a pivot that is not positive and finite).
bool qform(const double *u, double a, const double *c, double b, int p,
           const double *r, double *w, double *v, double *q, double *ld) {
  for (int j = 0; j < p; ++j) {
    for (int i = j; i < p; ++i) w[i + p * j] = a * u[i + p * j] + b * c[i + p * j];
  }
  double l = 0.0;
  for (int j = 0; j < p; ++j) {
    double d = w[j + p * j];
    for (int k = 0; k < j; ++k) d -= w[j + p * k] * w[j + p * k];
    if (!(d > 0.0) || !std::isfinite(d)) return false;
    const double s = std::sqrt(d);
    w[j + p * j] = s;
    l += std::log(s);
    for (int i = j + 1; i < p; ++i) {
      double e = w[i + p * j];
      for (int k = 0; k < j; ++k) e -= w[i + p * k] * w[j + p * k];
      w[i + p * j] = e / s;
    }
  }
  // forward substitution L v = r; q = |v|^2
  double s2 = 0.0;
  for (int i = 0; i < p; ++i) {
    double e = r[i];
    for (int k = 0; k < i; ++k) e -= w[i + p * k] * v[k];
    v[i] = e / w[i + p * i];
    s2 += v[i] * v[i];
  }
  *q = s2;
  *ld = 2.0 * l;
  return std::isfinite(s2);
}

}  // namespace

extern "C" {

int clab_tlvc(const double *x, int n, int p, const int *item, int m,
              double *mu, double *u, double *c, double *work) {
  if (p < 1 || m < 2 || n <= m) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (item[i] < 0 || item[i] >= m) return CLAB_EARG;
  }
  for (int k = 0; k < p; ++k) {
    for (int i = 0; i < n; ++i) {
      if (!std::isfinite(x[i + static_cast<std::int64_t>(n) * k])) return CLAB_EDOM;
    }
  }
  double *mean = work;  // m x p item means
  double *cnt = work + static_cast<std::int64_t>(m) * p;  // m counts
  double *res = cnt + m;                                  // m residual sums
  for (int i = 0; i < m; ++i) cnt[i] = 0.0;
  for (int i = 0; i < n; ++i) cnt[item[i]] += 1.0;
  std::int64_t ssq = 0;
  for (int i = 0; i < m; ++i) {
    if (!(cnt[i] > 0.0)) return CLAB_EARG;  // an item without observations
    const std::int64_t t = static_cast<std::int64_t>(cnt[i]);
    ssq += t * t;
  }
  // item and overall means, each refined by one corrective pass (the mean
  // of the first pass's residuals is added back)
  for (int k = 0; k < p; ++k) {
    const double *xk = x + static_cast<std::int64_t>(n) * k;
    double *mk = mean + static_cast<std::int64_t>(m) * k;
    for (int i = 0; i < m; ++i) mk[i] = res[i] = 0.0;
    double t = 0.0;
    for (int i = 0; i < n; ++i) {
      mk[item[i]] += xk[i];
      t += xk[i];
    }
    for (int i = 0; i < m; ++i) mk[i] /= cnt[i];
    t /= n;
    double rt = 0.0;
    for (int i = 0; i < n; ++i) {
      rt += xk[i] - t;
      res[item[i]] += xk[i] - mk[item[i]];
    }
    mu[k] = t + rt / n;
    for (int i = 0; i < m; ++i) mk[i] += res[i] / cnt[i];
  }
  // pooled within cross-products / (n - m); between: sum n_i d_i d_i'
  for (int a = 0; a < p; ++a) {
    for (int b = 0; b <= a; ++b) {
      const double *xa = x + static_cast<std::int64_t>(n) * a;
      const double *xb = x + static_cast<std::int64_t>(n) * b;
      const double *ma = mean + static_cast<std::int64_t>(m) * a;
      const double *mb = mean + static_cast<std::int64_t>(m) * b;
      kutil::ksum w, s;
      for (int i = 0; i < n; ++i) w.add((xa[i] - ma[item[i]]) * (xb[i] - mb[item[i]]));
      for (int i = 0; i < m; ++i) s.add(cnt[i] * (ma[i] - mu[a]) * (mb[i] - mu[b]));
      const double uw = w.val() / (n - m);
      // E[S_B] = (m - 1) U + (n - sum n_i^2 / n) C (Searle et al. 1992)
      const std::int64_t nn = static_cast<std::int64_t>(n) * n;
      const double den = static_cast<double>(nn - ssq) / n;
      const double cb = (s.val() - (m - 1) * uw) / den;
      u[a + p * b] = u[b + p * a] = uw;
      c[a + p * b] = c[b + p * a] = cb;
    }
  }
  return CLAB_OK;
}

int clab_tlvlr(const double *u, const double *c, const double *mu, int p,
               const double *xa, double na, const double *xb, double nb,
               double *out, double *work) {
  if (p < 1) return CLAB_EARG;
  if (!kutil::pos(na) || !kutil::pos(nb)) return CLAB_EDOM;
  for (int k = 0; k < p; ++k) {
    if (!std::isfinite(mu[k]) || !std::isfinite(xa[k]) || !std::isfinite(xb[k]))
      return CLAB_EDOM;
    for (int i = k; i < p; ++i) {
      if (!std::isfinite(u[i + p * k]) || !std::isfinite(c[i + p * k])) return CLAB_EDOM;
    }
  }
  double *w = work;                                     // p x p factor
  double *v = work + static_cast<std::int64_t>(p) * p;  // p scratch
  double *res = v + p;                                  // p residual
  const double nt = na + nb;
  double q[4], ld[4];
  // numerator: difference of the item means, and the weighted mean
  for (int k = 0; k < p; ++k) res[k] = xa[k] - xb[k];
  if (!qform(u, 1.0 / na + 1.0 / nb, c, 0.0, p, res, w, v, &q[0], &ld[0]))
    return CLAB_EDOM;
  for (int k = 0; k < p; ++k) res[k] = (na * (xa[k] - mu[k]) + nb * (xb[k] - mu[k])) / nt;
  if (!qform(u, 1.0 / nt, c, 1.0, p, res, w, v, &q[1], &ld[1])) return CLAB_EDOM;
  // denominator: each item against the population
  for (int k = 0; k < p; ++k) res[k] = xa[k] - mu[k];
  if (!qform(u, 1.0 / na, c, 1.0, p, res, w, v, &q[2], &ld[2])) return CLAB_EDOM;
  for (int k = 0; k < p; ++k) res[k] = xb[k] - mu[k];
  if (!qform(u, 1.0 / nb, c, 1.0, p, res, w, v, &q[3], &ld[3])) return CLAB_EDOM;
  const double l = 0.5 * ((q[2] + q[3]) - (q[0] + q[1])) +
                   0.5 * ((ld[2] + ld[3]) - (ld[0] + ld[1]));
  out[0] = l;
  out[1] = l / std::log(10.0);
  return CLAB_OK;
}

}  // extern "C"
