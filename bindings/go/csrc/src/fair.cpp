// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, fairness: group confusion-matrix metrics of a binary risk tool
// with Wilson score intervals.
#include <cstdint>

#include "crimlab.h"
#include "kutil.h"

namespace {

// Wilson score interval for x successes in m > 0 trials at normal quantile
// z > 0: the set of p with (x - m p)^2 <= z^2 m p (1 - p). With
// d = m + z^2, centre c = (x + z^2 / 2) / d and half-width
// h = (z / d) sqrt(x (m - x) / m + z^2 / 4), the ends are c -+ h. Since
// c^2 - h^2 = x^2 / (m d), the lower end is x^2 / (m d) / (c + h), with no
// cancellation, and by symmetry in x <-> m - x the upper end is 1 minus the
// same expression for m - x. x = 0 gives lo = 0 and x = m gives hi = 1
// exactly.
void wilson(double x, double m, double z, double *lo, double *hi) {
  if (!(m > 0.0)) {
    *lo = *hi = kutil::nan();
    return;
  }
  const double z2 = z * z;
  const double d = m + z2;
  const double y = m - x;
  const double h = (z / d) * std::sqrt(x * y / m + 0.25 * z2);
  const double cx = (x + 0.5 * z2) / d;
  const double cy = (y + 0.5 * z2) / d;
  *lo = x * x / (m * d) / (cx + h);
  *hi = 1.0 - y * y / (m * d) / (cy + h);
}

inline double ratio(double a, double b) {
  return b > 0.0 ? a / b : kutil::nan();
}

}  // namespace

extern "C" {

int clab_fair(const double *score, const int *y, const int *grp, int n, int g,
              double cut, int base, double z, double *cnt, double *met,
              double *rat, double *lo, double *hi) {
  if (n < 1 || g < 1 || base < 0 || base >= g) return CLAB_EARG;
  if (std::isnan(cut)) return CLAB_EDOM;
  if (lo != nullptr && hi != nullptr && !kutil::pos(z)) return CLAB_EDOM;
  for (int i = 0; i < n; ++i) {
    if (grp[i] < 0 || grp[i] >= g) return CLAB_EARG;
    if ((y[i] != 0 && y[i] != 1) || std::isnan(score[i])) return CLAB_EDOM;
  }
  for (int j = 0; j < 5 * g; ++j) cnt[j] = 0.0;
  // cnt[0..3][j]: tp, fp, tn, fn of group j; cnt[4][j]: its size. Counts
  // below 2^31 are exact in a double.
  for (int i = 0; i < n; ++i) {
    const int pr = score[i] > cut ? 1 : 0;
    const int k = y[i] == 1 ? (pr == 1 ? 0 : 3) : (pr == 1 ? 1 : 2);
    cnt[k * g + grp[i]] += 1.0;
    cnt[4 * g + grp[i]] += 1.0;
  }
  for (int j = 0; j < g; ++j) {
    const double tp = cnt[j], fp = cnt[g + j], tn = cnt[2 * g + j],
                 fn = cnt[3 * g + j], m = cnt[4 * g + j];
    // successes and trials of the eight proportions, in output order
    const double xs[8] = {tp + tn, tp, tn, fp, fn, tp, tn, tp + fp};
    const double ms[8] = {m, tp + fn, tn + fp, fp + tn, fn + tp,
                          tp + fp, tn + fn, m};
    for (int k = 0; k < 8; ++k) {
      met[k * g + j] = ratio(xs[k], ms[k]);
      if (lo != nullptr && hi != nullptr) {
        wilson(xs[k], ms[k], z, &lo[k * g + j], &hi[k * g + j]);
      }
    }
    met[8 * g + j] = tp + fp;
    // Matthews correlation: the numerator tp tn - fp fn in 64-bit integers
    // (each count < 2^31, each product < 2^62), so it is exact before its
    // one rounding; NaN when a margin is empty (the coefficient is 0 / 0).
    const std::int64_t a = static_cast<std::int64_t>(tp),
                       b = static_cast<std::int64_t>(tn),
                       c = static_cast<std::int64_t>(fp),
                       e = static_cast<std::int64_t>(fn);
    const double num = static_cast<double>(a * b - c * e);
    const double d1 = static_cast<double>((a + c) * (a + e));
    const double d2 = static_cast<double>((b + c) * (b + e));
    met[9 * g + j] = d1 > 0.0 && d2 > 0.0
                         ? num / (std::sqrt(d1) * std::sqrt(d2))
                         : kutil::nan();
  }
  // ratios to the base group: IEEE division (x / 0 is Inf, 0 / 0 is NaN)
  for (int k = 0; k < 10; ++k) {
    for (int j = 0; j < g; ++j) rat[k * g + j] = met[k * g + j] / met[k * g + base];
  }
  return CLAB_OK;
}

}  // extern "C"
