// SPDX-License-Identifier: AGPL-3.0-or-later
//
// Small predicates and a compensated sum shared by the kernel families
// outside core.cpp. Header only; nothing here allocates.
#ifndef CLAB_KUTIL_H
#define CLAB_KUTIL_H

#include <cmath>
#include <cstdint>
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

// In-place heapsort of a[0..n) under the strict order less(a, b): no
// recursion, O(n log n), the same order on every platform.
template <class T, class L>
void hsort(T *a, int n, L less) {
  auto sift = [&](int i, int m) {
    for (;;) {
      int c = 2 * i + 1;
      if (c >= m) return;
      if (c + 1 < m && less(a[c], a[c + 1])) ++c;
      if (!less(a[i], a[c])) return;
      const T t = a[i];
      a[i] = a[c];
      a[c] = t;
      i = c;
    }
  };
  for (int i = n / 2 - 1; i >= 0; --i) sift(i, n);
  for (int m = n - 1; m > 0; --m) {
    const T t = a[0];
    a[0] = a[m];
    a[m] = t;
    sift(0, m);
  }
}

// xorshift64* (Vigna 2016, doi:10.1145/2845077) seeded through splitmix64;
// 53-bit uniforms in [0, 1). Small, fast and identical everywhere.
struct rng {
  std::uint64_t s;
  explicit rng(std::uint64_t seed) {
    std::uint64_t z = seed + 0x9E3779B97F4A7C15ULL;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    s = (z ^ (z >> 31)) | 1ULL;  // never 0, the one fixed point
  }
  std::uint64_t next() {
    s ^= s >> 12;
    s ^= s << 25;
    s ^= s >> 27;
    return s * 0x2545F4914F6CDD1DULL;
  }
  double unif() { return static_cast<double>(next() >> 11) * 0x1.0p-53; }
  // an integer in [0, m), m >= 1, without modulo bias (Lemire's rejection
  // on the high bits would need 128-bit products; plain rejection instead)
  int below(int m) {
    const std::uint64_t mm = static_cast<std::uint64_t>(m);
    const std::uint64_t lim = (~0ULL) - (~0ULL) % mm;
    std::uint64_t v;
    do {
      v = next();
    } while (v >= lim);
    return static_cast<int>(v % mm);
  }
};

}  // namespace kutil

#endif
