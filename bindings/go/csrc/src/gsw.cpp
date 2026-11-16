// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Generalisability weights and diagnostics (crimlab.h, gensw()): the
// chance of being in the trial given the covariates (the sampling score,
// a logistic regression of s on x, clab_lgtreg), trimmed to [trim, 1 -
// trim]; the weights of the trial units that make them stand for a
// target, 1 / ps to generalise to the population the trial was drawn
// from (all rows) and (1 - ps) / ps to transport to the units outside it
// (Stuart et al. 2011; Westreich et al. 2017); and how far the trial is
// from the target: the standardised difference of the scores' logits
// (Stuart et al. 2011), Tipton's (2014) generalisability index over equal
// bins of the scores, Kish's effective size of the weights, and the
// standardised differences of each covariate before and after weighting.
#include <cmath>
#include <cstddef>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

extern "C" int clab_gensw(const double *x, int n, int p, const int *s, int target, double trim, int bins,
                          double *ps, double *w, double *coef, double *st, double *bal) CLAB_TRY {
  using std::size_t;
  if (n < 2 || p < 0 || target < 0 || target > 1 || !(trim >= 0 && trim < 0.5) || bins < 1) return CLAB_EARG;
  std::vector<double> sv(static_cast<size_t>(n));
  double n1 = 0;
  for (int i = 0; i < n; ++i) {
    if (s[i] != 0 && s[i] != 1) return CLAB_EARG;
    sv[static_cast<size_t>(i)] = s[i];
    n1 += s[i];
  }
  if (n1 < 1 || n1 > n - 1) return CLAB_EDOM;
  const size_t q = static_cast<size_t>(p) + 1;
  std::vector<double> se(q), work(2 * q * q + 4 * q);
  std::vector<int> iw(static_cast<size_t>(n));
  double dev;
  const int rc = clab_lgtreg(x, n, p, sv.data(), 1e-10, 50, coef, se.data(), &dev, work.data(), iw.data());
  if (rc != CLAB_OK) return rc;
  // the target: the rows outside the trial (transport) or all (generalise)
  auto tg = [&](int i) { return target == 1 ? s[i] == 0 : true; };
  double nt = 0;
  for (int i = 0; i < n; ++i) {
    double e = coef[0];
    for (int j = 0; j < p; ++j) e += coef[j + 1] * x[static_cast<size_t>(j) * static_cast<size_t>(n) + i];
    const double pr = std::fmin(1 - trim, std::fmax(trim, 1 / (1 + std::exp(-e))));
    ps[i] = pr;
    w[i] = s[i] ? (target == 1 ? (1 - pr) / pr : 1 / pr) : 0.0;
    if (tg(i)) nt += 1;
  }
  // the logits of the scores: mean in the trial, mean and sd in the target
  double l1 = 0, lt = 0, p1 = 0, pt = 0;
  for (int i = 0; i < n; ++i) {
    const double l = std::log(ps[i] / (1 - ps[i]));
    if (s[i]) {
      l1 += l;
      p1 += ps[i];
    }
    if (tg(i)) {
      lt += l;
      pt += ps[i];
    }
  }
  l1 /= n1;
  lt /= nt;
  double vt = 0;
  for (int i = 0; i < n; ++i) {
    if (!tg(i)) continue;
    const double l = std::log(ps[i] / (1 - ps[i])) - lt;
    vt += l * l;
  }
  const double sdt = nt > 1 ? std::sqrt(vt / (nt - 1)) : kutil::nan();
  // Tipton's index: sum over bins of sqrt(share in trial x share in
  // target), bins of equal width over the range of all scores
  double lo = ps[0], hi = ps[0];
  for (int i = 1; i < n; ++i) {
    lo = std::fmin(lo, ps[i]);
    hi = std::fmax(hi, ps[i]);
  }
  std::vector<double> f1(static_cast<size_t>(bins), 0.0), ft(static_cast<size_t>(bins), 0.0);
  for (int i = 0; i < n; ++i) {
    int b = hi > lo ? static_cast<int>(std::floor((ps[i] - lo) / (hi - lo) * bins)) : 0;
    if (b >= bins) b = bins - 1;
    if (s[i]) f1[static_cast<size_t>(b)] += 1 / n1;
    if (tg(i)) ft[static_cast<size_t>(b)] += 1 / nt;
  }
  double gi = 0, sw = 0, sw2 = 0;
  for (int b = 0; b < bins; ++b) gi += std::sqrt(f1[static_cast<size_t>(b)] * ft[static_cast<size_t>(b)]);
  for (int i = 0; i < n; ++i) {
    sw += w[i];
    sw2 += w[i] * w[i];
  }
  st[0] = n1;
  st[1] = nt;
  st[2] = (l1 - lt) / sdt;
  st[3] = p1 / n1 - pt / nt;
  st[4] = gi;
  st[5] = sw * sw / sw2;
  // each covariate: target mean, trial mean, weighted trial mean, and the
  // standardised differences of the two trial means from the target's
  for (int j = 0; j < p; ++j) {
    const double *xj = x + static_cast<size_t>(j) * static_cast<size_t>(n);
    double mt = 0, m1 = 0, mw = 0;
    for (int i = 0; i < n; ++i) {
      if (tg(i)) mt += xj[i];
      if (s[i]) {
        m1 += xj[i];
        mw += w[i] * xj[i];
      }
    }
    mt /= nt;
    m1 /= n1;
    mw /= sw;
    double v = 0;
    for (int i = 0; i < n; ++i) {
      if (tg(i)) v += (xj[i] - mt) * (xj[i] - mt);
    }
    const double sd = nt > 1 ? std::sqrt(v / (nt - 1)) : kutil::nan();
    bal[j] = mt;
    bal[p + j] = m1;
    bal[2 * p + j] = mw;
    bal[3 * p + j] = (m1 - mt) / sd;
    bal[4 * p + j] = (mw - mt) / sd;
  }
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)
