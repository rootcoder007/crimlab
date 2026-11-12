// SPDX-License-Identifier: AGPL-3.0-or-later
//
// libcrimlab, forensics: Scott-Knott grouping of measurements (glass
// refractive indices and the like) with a known measurement variance and
// a Monte Carlo test at every split.
#include <cstdint>

#include "crimlab.h"
#include "kutil.h"

namespace {

// in-place heapsort of a[0..n) under the strict order less(a, b); no
// recursion, O(n log n)
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

// B0 of sorted v[0..k): the largest between-groups sum of squares over the
// k - 1 splits into a lower and an upper part, C^2 k / (n1 n2) with C the
// running sum of the centred values; *at = size of the lower part (first
// maximum)
double bzero(const double *v, int k, int *at) {
  double m = 0.0;
  for (int i = 0; i < k; ++i) m += v[i];
  m /= k;
  double r = 0.0;
  for (int i = 0; i < k; ++i) r += v[i] - m;
  m += r / k;
  double c = 0.0, best = -1.0;
  *at = 1;
  for (int i = 1; i < k; ++i) {
    c += v[i - 1] - m;
    const double b = c * c * k / (static_cast<double>(i) * (k - i));
    if (b > best) {
      best = b;
      *at = i;
    }
  }
  return best;
}

struct rng {
  std::uint64_t s;
  explicit rng(std::uint64_t seed) {
    std::uint64_t z = seed + 0x9E3779B97F4A7C15ULL;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    s = (z ^ (z >> 31)) | 1ULL;
  }
  double unif() {
    s ^= s >> 12;
    s ^= s << 25;
    s ^= s >> 27;
    return static_cast<double>((s * 2685821657736338717ULL) >> 11) * 0x1p-53;
  }
  double norm() {
    const double a = 1.0 - unif(), b = unif();
    return std::sqrt(-2.0 * std::log(a)) * std::cos(6.283185307179586 * b);
  }
};

}  // namespace

extern "C" {

int clab_skp(const double *x, int n, double s2, double alpha, int nsim,
             unsigned int seed, int *grp, int *ng, int *iwork, double *work) {
  if (n < 1 || nsim < 1) return CLAB_EARG;
  if (!kutil::pos(s2) || !(alpha > 0.0 && alpha < 1.0)) return CLAB_EDOM;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(x[i])) return CLAB_EDOM;
  }
  int *idx = iwork;          // n: the sorting permutation
  int *stk = iwork + n;      // 2n: segments [lo, hi) still to examine
  double *v = work;          // n: sorted values
  double *z = work + n;      // n: simulated samples
  for (int i = 0; i < n; ++i) idx[i] = i;
  hsort(idx, n, [x](int a, int b) { return x[a] < x[b] || (x[a] == x[b] && a < b); });
  for (int i = 0; i < n; ++i) v[i] = x[idx[i]];
  const double c0 = 3.141592653589793 / (2.0 * (3.141592653589793 - 2.0));
  const double tol = 1.0 - 64 * std::numeric_limits<double>::epsilon();
  rng r(seed);
  int top = 0, g = 0;
  stk[top++] = 0;
  stk[top++] = n;
  while (top > 0) {
    const int hi = stk[--top], lo = stk[--top];
    const int k = hi - lo;
    bool split = false;
    int at = 1;
    if (k > 1) {
      const double b0 = bzero(v + lo, k, &at);
      const double lam = c0 * b0 / s2;
      if (b0 > 0.0) {
        // p = (1 + #{lambda* >= lambda}) / (nsim + 1) under k independent
        // N(mu, s2) values (lambda does not depend on mu or the scale);
        // stop as soon as p must exceed alpha
        const double cut = alpha * (nsim + 1.0) - 1.0;
        int ge = 0;
        for (int b = 0; b < nsim && ge <= cut; ++b) {
          for (int i = 0; i < k; ++i) z[i] = r.norm();
          hsort(z, k, [](double a, double c) { return a < c; });
          int a2;
          if (c0 * bzero(z, k, &a2) >= lam * tol) ++ge;
        }
        split = (1.0 + ge) / (nsim + 1.0) <= alpha;
      }
    }
    if (split) {
      // upper part pushed first, so the lower part is finished first and
      // the groups are numbered in increasing order of their values
      stk[top++] = lo + at;
      stk[top++] = hi;
      stk[top++] = lo;
      stk[top++] = lo + at;
    } else {
      ++g;
      for (int i = lo; i < hi; ++i) grp[idx[i]] = g;
    }
  }
  *ng = g;
  return CLAB_OK;
}

}  // extern "C"
