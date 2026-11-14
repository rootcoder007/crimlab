// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, bounds on the mean outcome under each of ordered treatments
// (sentences of increasing severity) with no assumption, under monotone
// treatment response, under monotone treatment selection, and under both
// (Manski and Pepper 2000).
#include <algorithm>
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

int clab_mtrbnd(const double *y, const int *z, const double *w, int n, int k, double ylo, double yhi, int as,
                double *out) {
  if (n < 1 || k < 2 || as < 0 || as > 3) return CLAB_EARG;
  if (!std::isfinite(ylo) || !std::isfinite(yhi) || !(ylo < yhi)) return CLAB_EDOM;
  std::vector<double> P(static_cast<std::size_t>(k), 0.0), m(static_cast<std::size_t>(k), 0.0);
  double W = 0;
  for (int i = 0; i < n; ++i) {
    const double wi = w ? w[i] : 1.0;
    if (z[i] < 0 || z[i] >= k) return CLAB_EARG;
    if (!(y[i] >= ylo && y[i] <= yhi) || !kutil::nonneg(wi)) return CLAB_EDOM;
    P[static_cast<std::size_t>(z[i])] += wi;
    m[static_cast<std::size_t>(z[i])] += wi * y[i];
    W += wi;
  }
  if (!(W > 0)) return CLAB_EDOM;
  for (int t = 0; t < k; ++t) {
    const std::size_t u = static_cast<std::size_t>(t);
    if (!(P[u] > 0)) return CLAB_EDOM;
    m[u] /= P[u];
    P[u] /= W;
  }
  // out: k x 6, column-major: share, mean, lower and upper bound of the
  // mean outcome under the treatment, lower and upper bound of its
  // difference from the first treatment's
  std::vector<double> lo(static_cast<std::size_t>(k)), hi(static_cast<std::size_t>(k));
  for (int t = 0; t < k; ++t) {
    double below = 0, above = 0, mbelow = 0, mabove = 0;  // shares and sums of m P
    for (int u = 0; u < k; ++u) {
      const std::size_t q = static_cast<std::size_t>(u);
      if (u < t) {
        below += P[q];
        mbelow += m[q] * P[q];
      } else if (u > t) {
        above += P[q];
        mabove += m[q] * P[q];
      }
    }
    const std::size_t s = static_cast<std::size_t>(t);
    const double mt = m[s], pt = P[s];
    double l = 0, h = 0;
    switch (as) {
      case 0:  // none
        l = mt * pt + ylo * (1 - pt);
        h = mt * pt + yhi * (1 - pt);
        break;
      case 1:  // monotone response: y(t) >= y(u) for u < t
        l = mbelow + mt * pt + ylo * above;
        h = yhi * below + mt * pt + mabove;
        break;
      case 2:  // monotone selection: E[y(t) | z = u] rises with u
        l = ylo * below + mt * (pt + above);
        h = mt * (below + pt) + yhi * above;
        break;
      default:  // both
        l = mbelow + mt * (pt + above);
        h = mt * (below + pt) + mabove;
        break;
    }
    lo[s] = l;
    hi[s] = h;
  }
  for (int t = 0; t < k; ++t) {
    const std::size_t s = static_cast<std::size_t>(t);
    out[t] = P[s];
    out[k + t] = m[s];
    out[2 * k + t] = lo[s];
    out[3 * k + t] = hi[s];
    double dl = lo[s] - hi[0];
    if ((as == 1 || as == 3) && dl < 0) dl = 0;
    out[4 * k + t] = t == 0 ? kutil::nan() : dl;
    out[5 * k + t] = t == 0 ? kutil::nan() : hi[s] - lo[0];
  }
  return CLAB_OK;
}
