// SPDX-License-Identifier: AGPL-3.0-or-later
//
// libcrimlab, causal: Rosenbaum sensitivity analysis for matched pairs.
// Bounds on the signed-rank p-value with the critical Gamma in closed
// form, and exact Hodges-Lehmann bounds from Walsh-average order
// statistics.
#include <cstdint>
#include <cstring>

#include "crimlab.h"
#include "kutil.h"

namespace {

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

// upper normal tail 1 - Phi(z) without cancellation
inline double qnormt(double z) { return 0.5 * std::erfc(z / std::sqrt(2.0)); }

// order-preserving map of doubles to unsigned integers and back
inline std::uint64_t key(double x) {
  std::uint64_t u;
  std::memcpy(&u, &x, sizeof u);
  return (u >> 63) ? ~u : u | 0x8000000000000000ULL;
}
inline double unkey(std::uint64_t k) {
  const std::uint64_t u = (k >> 63) ? k & 0x7FFFFFFFFFFFFFFFULL : ~k;
  double x;
  std::memcpy(&x, &u, sizeof x);
  return x;
}

// Walsh average of sorted v[i], v[j]: halves first, so no overflow; the
// rounding is monotone in each argument, which the counting relies on
inline double walsh(const double *v, int i, int j) { return 0.5 * v[i] + 0.5 * v[j]; }

// #{i <= j : walsh(i, j) <= t} for sorted v, by two pointers: O(n)
std::int64_t wcount(const double *v, int n, double t) {
  std::int64_t c = 0;
  int j = n - 1;
  for (int i = 0; i < n; ++i) {
    while (j >= i && walsh(v, i, j) > t) --j;
    if (j < i) break;
    c += j - i + 1;
  }
  return c;
}

// the k-th smallest Walsh average (1 <= k <= n (n + 1) / 2): the smallest
// double t with wcount(t) >= k, by bisection on the ordered keys of the
// doubles between v[0] and v[n - 1]; at most 64 steps of O(n). Every
// Walsh average is a double, so t is one of them.
double wkth(const double *v, int n, std::int64_t k) {
  std::uint64_t lo = key(v[0]), hi = key(v[n - 1]);
  while (lo < hi) {
    const std::uint64_t mid = lo + (hi - lo) / 2;
    if (wcount(v, n, unkey(mid)) >= k) hi = mid; else lo = mid + 1;
  }
  return unkey(lo);
}

// For a = g M / (1 + g) (up = true) or a = M / (1 + g) (up = false), the
// integers ceil(a) and floor(a) + 1, decided exactly: k >= g M / (1 + g)
// iff g (M - k) - k <= 0, and k >= M / (1 + g) iff g k + (k - M) >= 0;
// fma rounds each once, so the signs are exact (M < 2^53).
// cmp(k) > 0 when k > a, 0 when k = a, < 0 when k < a.
int cmpa(std::int64_t k, std::int64_t m, double g, bool up) {
  const double kd = static_cast<double>(k);
  const double s = up ? -std::fma(g, static_cast<double>(m - k), -kd)
                      : std::fma(g, kd, static_cast<double>(k - m));
  return s > 0.0 ? 1 : (s < 0.0 ? -1 : 0);
}

void ceilfl(std::int64_t m, double g, bool up, std::int64_t *ce, std::int64_t *f1) {
  const double a = up ? g * static_cast<double>(m) / (1.0 + g) : static_cast<double>(m) / (1.0 + g);
  std::int64_t k = static_cast<std::int64_t>(std::ceil(a));
  while (k > 1 && cmpa(k - 1, m, g, up) >= 0) --k;
  while (cmpa(k, m, g, up) < 0) ++k;
  *ce = k;                                  // smallest k >= a
  *f1 = cmpa(k, m, g, up) > 0 ? k : k + 1;  // smallest k > a
}

}  // namespace

extern "C" {

int clab_rswp(const double *d, int n, const double *gam, int ng, double z,
              double *plo, double *pup, double *gcrit, double *tsum,
              double *work, int *iwork) {
  if (n < 1 || ng < 0) return CLAB_EARG;
  if (!kutil::pos(z)) return CLAB_EDOM;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(d[i])) return CLAB_EDOM;
  }
  for (int k = 0; k < ng; ++k) {
    if (!(gam[k] >= 1.0) || !std::isfinite(gam[k])) return CLAB_EDOM;
  }
  // the non-zero differences, ordered by absolute value; midranks for ties
  int m = 0;
  for (int i = 0; i < n; ++i) {
    if (d[i] != 0.0) iwork[m++] = i;
  }
  hsort(iwork, m, [d](int a, int b) {
    const double x = std::fabs(d[a]), y = std::fabs(d[b]);
    return x < y || (x == y && a < b);
  });
  for (int i = 0; i < m;) {
    int j = i;
    while (j + 1 < m && std::fabs(d[iwork[j + 1]]) == std::fabs(d[iwork[i]])) ++j;
    const double r = 0.5 * (i + 1 + j + 1);
    for (int k = i; k <= j; ++k) work[k] = r;
    i = j + 1;
  }
  kutil::ksum s2, tp, tn;
  for (int k = 0; k < m; ++k) {
    s2.add(work[k] * work[k]);
    if (d[iwork[k]] > 0.0) tp.add(work[k]); else tn.add(work[k]);
  }
  const double t = tp.val(), tneg = tn.val(), s1 = t + tneg, v2 = s2.val();
  *tsum = t;
  for (int k = 0; k < ng; ++k) {
    const double pp = gam[k] / (1.0 + gam[k]), pm = 1.0 / (1.0 + gam[k]);
    // E and Var of the bounding sum: p s1 and p (1 - p) s2
    pup[k] = m > 0 ? qnormt((t - pp * s1) / std::sqrt(pp * pm * v2)) : kutil::nan();
    plo[k] = m > 0 ? qnormt((t - pm * s1) / std::sqrt(pp * pm * v2)) : kutil::nan();
  }
  // the Gamma at which the upper bound equals 1 - Phi(z): the smaller root
  // of (T - p S1)^2 = z^2 S2 p (1 - p), p = 2 T^2 / D with
  // D = 2 T S1 + z^2 S2 + z sqrt(S2 (z^2 S2 + 4 T (S1 - T))); then
  // Gamma = p / (1 - p) = 2 T^2 / (D - 2 T^2), and
  // D - 2 T^2 = 2 T (S1 - T) + z^2 S2 + z sqrt(...) has no cancellation.
  if (m > 0) {
    const double z2 = z * z;
    const double rt = z * std::sqrt(v2 * (z2 * v2 + 4.0 * t * tneg));
    *gcrit = 2.0 * t * t / (2.0 * t * tneg + z2 * v2 + rt);
  } else {
    *gcrit = kutil::nan();
  }
  return CLAB_OK;
}

int clab_rshl(const double *d, int n, const double *gam, int ng, double *lo,
              double *hi, double *work) {
  if (n < 1 || n > (1 << 26) || ng < 0) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(d[i])) return CLAB_EDOM;
    work[i] = d[i];
  }
  for (int k = 0; k < ng; ++k) {
    if (!(gam[k] >= 1.0) || !std::isfinite(gam[k])) return CLAB_EDOM;
  }
  hsort(work, n, [](double a, double b) { return a < b; });
  const std::int64_t m = static_cast<std::int64_t>(n) * (n + 1) / 2;
  for (int k = 0; k < ng; ++k) {
    // upper bound at mu = M / (1 + Gamma), a = M - mu; lower at
    // mu = Gamma M / (1 + Gamma), a = M - mu; each the midpoint of
    // W_(ceil a) and W_(floor a + 1)
    for (int side = 0; side < 2; ++side) {
      std::int64_t ce, f1;
      ceilfl(m, gam[k], side == 1, &ce, &f1);
      const double w1 = wkth(work, n, ce);
      const double w2 = f1 == ce ? w1 : wkth(work, n, f1);
      (side == 1 ? hi : lo)[k] = 0.5 * w1 + 0.5 * w2;
    }
  }
  return CLAB_OK;
}

}  // extern "C"
