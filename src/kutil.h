// SPDX-License-Identifier: AGPL-3.0-or-later
//
// Small predicates and a compensated sum shared by the kernel families
// outside core.cpp. Header only; nothing here allocates.
#ifndef CLAB_KUTIL_H
#define CLAB_KUTIL_H

#include <cmath>
#include <limits>

namespace kutil {

inline bool pos(double x) noexcept { return std::isfinite(x) && x > 0.0; }
inline bool nonneg(double x) noexcept { return std::isfinite(x) && x >= 0.0; }
inline double nan() noexcept {
  return std::numeric_limits<double>::quiet_NaN();
}

// Neumaier's compensated sum: s + c carries the exact sum of the terms to
// within one rounding of the result for any order of magnitudes (Neumaier
// 1974, doi:10.1002/zamm.19740540106).
struct ksum {
  double s = 0.0, c = 0.0;
  void add(double x) noexcept {
    const double t = s + x;
    c += std::fabs(s) >= std::fabs(x) ? (s - t) + x : (x - t) + s;
    s = t;
  }
  double val() const noexcept { return s + c; }
};

// A tail of Binomial(n, p) starting at g: P(X >= g) when up, P(X <= g)
// otherwise, for g on the far side of the mean (the terms then decrease
// monotonically away from g). From the log pmf at g by the ratio
// recursion, summed until a term is below 1e-17 of the sum: O(sd) terms.
// Needs p > 0 (both callers have it); p >= 1 is a point mass at n.
inline double btail(double g, double n, double p, bool up) {
  if (p >= 1.0) return up ? 1.0 : (g >= n ? 1.0 : 0.0);
  const double lt = std::lgamma(n + 1.0) - std::lgamma(g + 1.0) - std::lgamma(n - g + 1.0) +
                    g * std::log(p) + (n - g) * std::log1p(-p);
  double t = std::exp(lt), s = 0.0;
  const double r = p / (1.0 - p);
  for (double k = g; up ? k <= n : k >= 0.0; k += up ? 1.0 : -1.0) {
    s += t;
    if (t < 1e-17 * s) break;
    t *= up ? (n - k) / (k + 1.0) * r : k / (n - k + 1.0) / r;
  }
  return s;
}

}  // namespace kutil

#endif
