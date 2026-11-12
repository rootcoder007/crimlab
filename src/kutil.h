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
// monotonically away from g); O(sd) terms. Needs p > 0 (both callers have
// it); p >= 1 is a point mass at n. Defined in kutil.cpp (one translation
// unit, so its coverage is counted once).
double btail(double g, double n, double p, bool up);

}  // namespace kutil

#endif
