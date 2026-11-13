// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, forensics: Benford digit screens. Leading digits extracted
// exactly (no truncation of a rounded product), zeros and non-finite values
// set aside and counted, and the goodness-of-fit tests with a closed-form
// chi-square tail and a +1-corrected Monte Carlo p-value.
#include <cstdint>

#include "crimlab.h"
#include "kutil.h"

namespace {

// 10^k for 0 <= k <= 22: every one is exact in binary64 (5^22 < 2^53)
double p10(int k) {
  double r = 1.0;
  for (int i = 0; i < k; ++i) r *= 10.0;
  return r;
}

std::int64_t i10(int k) {
  std::int64_t r = 1;
  for (int i = 0; i < k; ++i) r *= 10;
  return r;
}

int ndig(std::int64_t k) {
  int d = 0;
  while (k > 0) {
    k /= 10;
    ++d;
  }
  return d;
}

// the first nd digits of the positive integer k (padded with zeros)
std::int64_t lead(std::int64_t k, int nd) {
  const int d = ndig(k);
  return d >= nd ? k / i10(d - nd) : k * i10(nd - d);
}

// a * 10^p rounded once: a multiplication or a division by an exact power
double scl(double a, int p) { return p >= 0 ? a * p10(p) : a / p10(-p); }

// The first nd digits D of the shortest decimal S that rounds to a > 0,
// with no need to find S. D 10^q <= S < (D + 1) 10^q is decided by
// comparing a with the doubles nearest to the two boundaries: D and 10^|q|
// (|q| <= 22) are exact, so each boundary B is correctly rounded to r(B).
// If a > r(B) then a > B, and if a < r(B) then a < B (rounding is
// monotone and B lies within half an ulp of r(B)); S lies in the rounding
// interval of a, so it is on the same side of B as a. If r(B) = a, then B
// itself rounds to a and has at most nd digits, so S = B (two distinct
// decimals of at most 9 digits cannot round to one double), and B starts
// the digit block it bounds from below. Hence D is the one prefix with
// r(D 10^q) <= a < r((D + 1) 10^q). Returns D, or -1 when the scaling
// leaves the exact range (|a| below about 10^(nd - 23) or above 10^(nd + 21)).
std::int64_t prefix(double a, int nd) {
  const int e0 = static_cast<int>(std::floor(std::log10(a)));
  const double lo = p10(nd - 1), hi = p10(nd);
  for (int e = e0 - 1; e <= e0 + 1; ++e) {
    const int p = nd - 1 - e;  // D = a 10^p, roughly
    if (p > 22 || p < -22) continue;
    const double d0 = std::floor(scl(a, p));
    if (d0 > hi + 1.0) continue;
    for (int o = -1; o <= 1; ++o) {
      const double d = d0 + o;
      if (d < lo || d >= hi) continue;
      if (scl(d, -p) <= a && a < scl(d + 1.0, -p)) return static_cast<std::int64_t>(d);
    }
  }
  return -1;
}

// chi-square upper tail Q(df / 2, x / 2) in closed form: for even df,
// e^-y sum_{j < df/2} y^j / j!; for odd df, erfc(sqrt y) plus
// e^-y sum_{j <= (df - 3)/2} y^(j + 1/2) / Gamma(j + 3/2) (y = x / 2;
// Abramowitz and Stegun 26.4.4, 26.4.5). Terms are formed in logs and
// summed from the largest, so the tail keeps its relative accuracy far
// below 1e-300.
double chisqq(double x, int df) {
  if (!(x > 0.0)) return 1.0;
  const double y = 0.5 * x, ly = std::log(y);
  const bool odd = df % 2 == 1;
  const int nt = odd ? (df - 1) / 2 : df / 2;  // number of terms
  // log of term j: -y + (j + h) ly - lgam(j + h + 1), h = 0 or 1/2
  const double h = odd ? 0.5 : 0.0;
  double lg = odd ? std::log(0.886226925452758014) : 0.0;  // lnGamma(1 + h)
  double mx = -HUGE_VAL;
  for (int j = 0; j < nt; ++j) {
    const double l = -y + (j + h) * ly - lg;
    if (l > mx) mx = l;
    lg += std::log(j + h + 1.0);
  }
  kutil::ksum s;
  lg = odd ? std::log(0.886226925452758014) : 0.0;
  for (int j = 0; j < nt; ++j) {
    s.add(std::exp(-y + (j + h) * ly - lg - mx));
    lg += std::log(j + h + 1.0);
  }
  const double t = nt > 0 ? std::exp(mx) * s.val() : 0.0;
  return odd ? std::erfc(std::sqrt(y)) + t : t;
}

// Benford probability of the nd-digit prefix d: log10(1 + 1 / d)
inline double pb(std::int64_t d) {
  return std::log1p(1.0 / static_cast<double>(d)) / std::log(10.0);
}

double chistat(const double *obs, int k, std::int64_t lo, double m) {
  kutil::ksum c;
  for (int j = 0; j < k; ++j) {
    const double e = m * pb(lo + j);
    c.add((obs[j] - e) * (obs[j] - e) / e);
  }
  return c.val();
}

}  // namespace

extern "C" {

int clab_bdig(const double *x, int n, int nd, int dp, int sgn, int *dig,
              int *cnt) {
  if (n < 0 || nd < 1 || nd > 9 || dp < -1 || dp > 22 || sgn < -1 || sgn > 1)
    return CLAB_EARG;
  for (int k = 0; k < 4; ++k) cnt[k] = 0;
  const double sc = dp >= 0 ? p10(dp) : 0.0;
  for (int i = 0; i < n; ++i) {
    dig[i] = 0;
    const double v = x[i];
    if (!std::isfinite(v)) {
      ++cnt[1];
      continue;
    }
    if ((sgn == 1 && v < 0.0) || (sgn == -1 && v > 0.0)) {
      ++cnt[2];
      continue;
    }
    const double a = std::fabs(v);
    std::int64_t k;
    if (dp >= 0) {
      // recorded to dp decimals: k = round(a 10^dp). For a the double
      // nearest to k 10^-dp, a 10^dp = k (1 + d) with |d| <= 2^-52, so the
      // error is below 1/2 while k < 2^50; 2^49 leaves a margin.
      const double y = a * sc;
      if (y >= 0x1p49) {
        ++cnt[3];
        continue;
      }
      k = static_cast<std::int64_t>(std::floor(y + 0.5));
    } else {
      if (a == 0.0) {
        k = 0;
      } else {
        k = prefix(a, nd);
        if (k < 0) {
          ++cnt[3];
          continue;
        }
      }
    }
    if (k == 0) {
      ++cnt[0];
      continue;
    }
    dig[i] = static_cast<int>(lead(k, nd));
  }
  return CLAB_OK;
}

int clab_bnft(const int *dig, int n, int nd, int nsim, unsigned int seed,
              double *obs, double *work, double *out) {
  if (n < 0 || nd < 1 || nd > 6 || nsim < 0) return CLAB_EARG;
  const std::int64_t lo = i10(nd - 1);
  const int k = static_cast<int>(9 * lo);
  for (int j = 0; j < k; ++j) obs[j] = 0.0;
  double m = 0.0;
  for (int i = 0; i < n; ++i) {
    if (dig[i] == 0) continue;
    if (dig[i] < lo || dig[i] >= 10 * lo) return CLAB_EDOM;
    obs[dig[i] - lo] += 1.0;
    m += 1.0;
  }
  out[0] = m;
  out[2] = k - 1;
  if (m == 0.0) {
    out[1] = out[3] = out[4] = out[5] = kutil::nan();
    return CLAB_OK;
  }
  const double chi = chistat(obs, k, lo, m);
  out[1] = chi;
  out[3] = chisqq(chi, k - 1);
  kutil::ksum mad;
  for (int j = 0; j < k; ++j) mad.add(std::fabs(obs[j] / m - pb(lo + j)));
  out[5] = mad.val() / k;
  if (nsim == 0) {
    out[4] = kutil::nan();
    return CLAB_OK;
  }
  // Monte Carlo under Benford: a prefix is floor(10^(nd - 1 + u)), the
  // inverse of the prefix distribution function log10(d / 10^(nd - 1)).
  // p = (1 + #{chi* >= chi}) / (nsim + 1), never 0 (North et al. 2002);
  // ties are taken within 64 ulps.
  kutil::rng r(seed);
  const double hi = chi * (1.0 - 64 * std::numeric_limits<double>::epsilon());
  const std::int64_t mi = static_cast<std::int64_t>(m);
  int ge = 0;
  for (int b = 0; b < nsim; ++b) {
    for (int j = 0; j < k; ++j) work[j] = 0.0;
    for (std::int64_t i = 0; i < mi; ++i) {
      std::int64_t d = static_cast<std::int64_t>(std::floor(std::pow(10.0, nd - 1 + r.unif())));
      if (d < lo) d = lo;
      if (d >= 10 * lo) d = 10 * lo - 1;
      work[d - lo] += 1.0;
    }
    if (chistat(work, k, lo, m) >= hi) ++ge;
  }
  out[4] = (1.0 + ge) / (nsim + 1.0);
  return CLAB_OK;
}

}  // extern "C"
