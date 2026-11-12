// SPDX-License-Identifier: AGPL-3.0-or-later
//
// libcrimlab: the shared numerical helpers of kutil.h that are not inline.
#include "kutil.h"

namespace kutil {

// A tail of Binomial(n, p) starting at g: P(X >= g) when up, P(X <= g)
// otherwise, for g on the far side of the mean (the terms then decrease
// monotonically away from g). From the log pmf at g by the ratio
// recursion, summed until a term is below 1e-17 of the sum: O(sd) terms.
// Needs p > 0 (both callers have it); p >= 1 is a point mass at n.
double btail(double g, double n, double p, bool up) {
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
