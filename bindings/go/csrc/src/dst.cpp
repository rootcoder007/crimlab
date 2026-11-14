// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, distributions: the normal, chi-square (central and
// noncentral), Student t, F, gamma, beta, Poisson and binomial
// distribution functions, their upper tails and quantiles, for the front
// ends outside R that turn a level into a quantile or a statistic into a
// p-value, as the R binding does with R's own functions. A lower tail is
// lower != 0, an upper tail lower == 0; arguments outside the domain give
// NaN. The results follow R's definitions, including its conventions at
// the edges (df 0, a = 0, infinite degrees of freedom).
//
// The incomplete gamma ratio is a series (x below a + 1) or a continued
// fraction, the incomplete beta ratio a continued fraction, both after
// the prefactor x^a e^-x / Gamma(a + 1) or x^a (1 - x)^b / B(a, b) is
// computed by Loader's saddle-point form (C. Loader, Fast and accurate
// computation of binomial probabilities, 2000), which keeps its relative
// precision when a and b are large. Quantiles are found by Newton steps
// on the logarithm of the tail, kept inside a bracket.
#include <cfloat>
#include <cmath>
#include <limits>

#include "crimlab.h"

namespace {

const double qnan = std::numeric_limits<double>::quiet_NaN();
const double eps = DBL_EPSILON;
const double ln2pi = 1.837877066409345483560659472811;
const double lnsqrt2pi = 0.918938533204672741780329736406;
const double euler = 0.5772156649015328606065120900824024;
const double tiny = 1e-300;
const double pi = 3.14159265358979323846;

// zeta(k) - 1 for k = 2 .. 61: Euler-Maclaurin from N = 20.
struct Zeta {
  double z[60];
  Zeta() {
    static const double b[] = {1.0 / 6, -1.0 / 30, 1.0 / 42, -1.0 / 30, 5.0 / 66, -691.0 / 2730, 7.0 / 6};
    for (int i = 0; i < 60; ++i) {
      const double k = i + 2;
      double s = 0;
      for (int n = 19; n >= 2; --n) s += std::pow(n, -k);
      const double N = 20;
      double t = std::pow(N, 1 - k) / (k - 1) + 0.5 * std::pow(N, -k), f = 1, rf = k;
      for (int j = 1; j <= 7; ++j) {
        // B_2j / (2j)! times k (k + 1) ... (k + 2j - 2) N^(-k - 2j + 1)
        f /= (2.0 * j - 1) * (2.0 * j);
        t += b[j - 1] * f * rf * std::pow(N, -k - 2 * j + 1);
        rf *= (k + 2 * j - 1) * (k + 2 * j);
      }
      z[i] = s + t;
    }
  }
};

// log(Gamma(1 + a)), precise near a = 0.
double lgamma1p(double a) {
  if (std::fabs(a) >= 0.5) return std::lgamma(a + 1);
  static const Zeta zt;
  double s = 0, p = -a;
  for (int k = 2; k < 62; ++k) {
    p *= -a;
    const double t = p * zt.z[k - 2] / k;
    s += t;
    if (std::fabs(t) < eps * 1e-3 * std::fabs(s)) break;
  }
  return -euler * a + (a - std::log1p(a)) + s;
}

// lgamma(n + 1) - (n + 1/2) log(n) + n - log(sqrt(2 pi)).
double stirlerr(double n) {
  if (n <= 15) {
    const double l = n < 0.5 ? lgamma1p(n) : std::lgamma(n + 1);
    return l - (n + 0.5) * std::log(n) + n - lnsqrt2pi;
  }
  const double S0 = 1.0 / 12, S1 = 1.0 / 360, S2 = 1.0 / 1260, S3 = 1.0 / 1680, S4 = 1.0 / 1188;
  const double nn = n * n;
  if (n > 500) return (S0 - S1 / nn) / n;
  if (n > 80) return (S0 - (S1 - S2 / nn) / nn) / n;
  if (n > 35) return (S0 - (S1 - (S2 - S3 / nn) / nn) / nn) / n;
  return (S0 - (S1 - (S2 - (S3 - S4 / nn) / nn) / nn) / nn) / n;
}

// x log(x / np) + np - x, precise when x is near np.
double bd0(double x, double np) {
  if (std::fabs(x - np) < 0.1 * (x + np)) {
    double v = (x - np) / (x + np), s = (x - np) * v;
    if (std::fabs(s) < DBL_MIN) return s;
    double ej = 2 * x * v;
    v *= v;
    for (int j = 1; j < 1000; ++j) {
      ej *= v;
      const double s1 = s + ej / (2 * j + 1);
      if (s1 == s) return s1;
      s = s1;
    }
  }
  return x * std::log(x / np) + np - x;
}

// x^a e^-x / Gamma(a + 1) for real a >= 0.
double dpois(double a, double x) {
  if (x == 0) return a == 0 ? 1 : 0;
  if (!std::isfinite(x)) return 0;
  if (a <= x * DBL_MIN) return std::exp(-x);
  if (x < a * DBL_MIN) return std::exp(-x + a * std::log(x) - std::lgamma(a + 1));
  return std::exp(-stirlerr(a) - bd0(a, x)) / std::sqrt(2 * pi * a);
}

// choose(n, x) p^x q^(n - x) for real x and n, q = 1 - p.
double dbin(double x, double n, double p, double q) {
  if (p == 0) return x == 0 ? 1 : 0;
  if (q == 0) return x == n ? 1 : 0;
  if (x == 0) {
    if (n == 0) return 1;
    return std::exp(p < 0.1 ? -bd0(n, n * q) - n * p : n * std::log(q));
  }
  if (x == n) return std::exp(q < 0.1 ? -bd0(n, n * p) - n * q : n * std::log(p));
  if (x < 0 || x > n) return 0;
  const double lc = stirlerr(n) - stirlerr(x) - stirlerr(n - x) - bd0(x, n * p) - bd0(n - x, n * q);
  const double lf = ln2pi + std::log(x) + std::log1p(-x / n);
  return std::exp(lc - 0.5 * lf);
}

// log(Gamma(b) / Gamma(a + b)) for b >= 16, any a >= 0.
double lgdiff(double a, double b) {
  return stirlerr(b) - stirlerr(a + b) - (b - 0.5) * std::log1p(a / b) - a * std::log(a + b) + a;
}

// lgamma(a + b) - lgamma(b), precise for small a.
double lgd(double a, double b) {
  double s = 0;
  for (; b < 16; b += 1) s += std::log1p(a / b);
  return -s - lgdiff(a, b);
}

// x^a y^b / B(a, b), y = 1 - x given apart for its precision.
double bfac(double x, double y, double a, double b) {
  if (x == 0 || y == 0) return 0;
  if (a >= 1 && b >= 1) {
    // Loader's form of choose(n, a) x^a y^b, n = a + b, times a b / n
    const double n = a + b;
    const double lc = stirlerr(n) - stirlerr(a) - stirlerr(b) - bd0(a, n * x) - bd0(b, n * y);
    return std::exp(lc) * std::sqrt(a * b / (2 * pi * n));
  }
  const double lx = x < 0.5 ? std::log(x) : std::log1p(-y), ly = y < 0.5 ? std::log(y) : std::log1p(-x);
  const double s = std::fmin(a, b), l = std::fmax(a, b);
  const double lb = std::lgamma(s) - lgd(s, l);
  return std::exp(a * lx + b * ly - lb);
}

// The continued fraction of I_x(a, b) (modified Lentz).
double bcf(double x, double a, double b) {
  const double qab = a + b, qap = a + 1, qam = a - 1;
  double c = 1, d = 1 - qab * x / qap;
  if (std::fabs(d) < tiny) d = tiny;
  d = 1 / d;
  double h = d;
  for (long m = 1; m < 100000000; ++m) {
    const double m2 = 2.0 * m;
    double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
    d = 1 + aa * d;
    if (std::fabs(d) < tiny) d = tiny;
    c = 1 + aa / c;
    if (std::fabs(c) < tiny) c = tiny;
    d = 1 / d;
    h *= d * c;
    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
    d = 1 + aa * d;
    if (std::fabs(d) < tiny) d = tiny;
    c = 1 + aa / c;
    if (std::fabs(c) < tiny) c = tiny;
    d = 1 / d;
    const double del = d * c;
    h *= del;
    if (std::fabs(del - 1) < eps) break;
  }
  return h;
}

// The 20 nodes and weights of Gauss-Legendre on [-1, 1].
struct Gl {
  double x[20], w[20];
  Gl() {
    const int n = 20;
    for (int i = 0; i < n; ++i) {
      double z = std::cos(pi * (i + 0.75) / (n + 0.5)), dp = 0;
      for (int it = 0; it < 100; ++it) {
        double p0 = 1, p1 = z;
        for (int k = 2; k <= n; ++k) {
          const double p2 = ((2 * k - 1) * z * p1 - (k - 1) * p0) / k;
          p0 = p1;
          p1 = p2;
        }
        dp = n * (z * p1 - p0) / (z * z - 1);
        const double dz = p1 / dp;
        z -= dz;
        if (std::fabs(dz) < 1e-17) break;
      }
      x[i] = z;
      w[i] = 2 / ((1 - z * z) * dp * dp);
    }
  }
};

// I_x(a, b) for x near 1 and a large, b small, where the continued
// fraction cancels: x^a y^b / (B(a, b) a y) times the integral over r
// from 0 to a of (1 - r / a)^(a - 1) (1 + x r / (a y))^(b - 1), by
// Gauss-Legendre on panels.
double ibquad(double x, double y, double a, double b) {
  static const Gl g;
  // panels as wide as four of the integrand's scales of decay
  const double wd = 4 * std::fmax(1.0, 0.5 * (b + 1));
  double s = 0;
  for (double lo = 0; lo < a; lo += wd) {
    const double h = std::fmin(lo + wd, a), m = 0.5 * (lo + h), r = 0.5 * (h - lo);
    double t = 0;
    for (int i = 0; i < 20; ++i) {
      const double u = m + r * g.x[i];
      t += g.w[i] * std::exp((a - 1) * std::log1p(-u / a) + (b - 1) * std::log1p(x * u / (a * y)));
    }
    t *= r;
    s += t;
    if (t <= 1e-18 * s) break;
  }
  return bfac(x, y, a, b) * s / (a * y);
}

// I_x(a, b) where x < (a + 1) / (a + b + 2).
double ibet(double x, double y, double a, double b) {
  if (a > 20 && y < 0.1 && b < 50 && a > 10 * b) return ibquad(x, y, a, b);
  return bfac(x, y, a, b) * bcf(x, a, b) / a;
}

// I_x(a, b) by its power series for b <= 1 (positive terms):
// x^a / (a B(a, b)) (1 + a sum_j (1 - b) ... (j - b) / j! x^j / (a + j));
// NaN when it would need more than 1e6 terms.
double bpser(double x, double y, double a, double b) {
  if (b > 1 || x > 0.9999) return qnan;
  double s = 0, c = 1;
  for (int j = 1; j < 1000000; ++j) {
    c *= (j - b) / j * x;
    const double t = c / (a + j);
    s += t;
    if (t <= eps * 1e-2 * s) return bfac(x, y, a, b) / (a * std::pow(y, b)) * (1 + a * s);
  }
  return qnan;
}

// 1 - I_x(a, b) for a < 1, x <= 1/2 and b x <= 2, by the series of
// bpser() taken to its logarithm: -expm1(a log x + lgamma(a + b) -
// lgamma(b) - lgamma(1 + a) + log1p(a sum)); NaN outside.
double apser(double x, double a, double b) {
  if (a >= 1 || x > 0.5 || b * x > 2) return qnan;
  double s = 0, c = 1;
  for (int j = 1; j < 100000; ++j) {
    c *= (j - b) / j * x;
    const double t = c / (a + j);
    s += t;
    if (std::fabs(t) <= eps * 1e-2 * std::fabs(s)) break;
  }
  return -std::expm1(a * std::log(x) + lgd(a, b) - lgamma1p(a) + std::log1p(a * s));
}

// I_x(a, b) and its complement; y = 1 - x.
double pbet(double x, double y, double a, double b, int lower) {
  if (std::isnan(x) || std::isnan(a) || std::isnan(b) || a < 0 || b < 0) return qnan;
  if (x <= 0) return lower ? 0 : 1;
  if (y <= 0) return lower ? 1 : 0;
  double lo;
  if (a == 0 || b == 0 || std::isinf(a) || std::isinf(b)) {
    // the limits: point masses at 0, 1 or 1/2
    if (a == 0 && b == 0) lo = 0.5;
    else if (a == 0 || a / b == 0) lo = 1;
    else if (b == 0 || b / a == 0) lo = 0;
    else lo = x < 0.5 ? 0 : 1;
    return lower ? lo : 1 - lo;
  }
  // w = I_x(a, b) on the side where its continued fraction converges
  const bool sw = !(x < (a + 1) / (a + b + 2));
  const double w = sw ? ibet(y, x, b, a) : ibet(x, y, a, b);
  // the other tail, 1 - w, loses digits when w is near 1: by the series
  if (w > 0.5) {
    double v = sw ? bpser(x, y, a, b) : bpser(y, x, b, a);
    if (std::isnan(v)) v = sw ? apser(y, b, a) : apser(x, a, b);
    if (!std::isnan(v)) return (lower != 0) == sw ? v : 1 - v;
  }
  return (lower != 0) == sw ? 1 - w : w;
}

// The regularized incomplete gamma ratio P(a, x), or Q(a, x).
double pgam(double x, double a, int lower) {
  if (std::isnan(x) || std::isnan(a) || a < 0) return qnan;
  if (x <= 0) return lower ? 0 : 1;
  if (a == 0 || std::isinf(x)) return lower ? 1 : 0;
  if (std::isinf(a)) return lower ? 0 : 1;
  if (a < 1 && x < 1) {
    // x^a / Gamma(a + 1) (1 + a sum_{n >= 1} (-x)^n / (n! (a + n)))
    double s = 0, c = 1;
    for (int n = 1; n < 1000; ++n) {
      c *= -x / n;
      const double t = c / (a + n);
      s += t;
      if (std::fabs(t) <= eps * std::fabs(s)) break;
    }
    const double lf = a * std::log(x) - lgamma1p(a), f = std::exp(lf);
    return lower ? f * (1 + a * s) : -std::expm1(lf) - f * a * s;
  }
  if (a >= 1 && x < a + 1) {
    // x^a e^-x / Gamma(a + 1) sum_{n >= 0} x^n / ((a + 1) ... (a + n))
    double s = 1, cmp = 0, t = 1, ap = a;
    for (long n = 1; n < 100000000; ++n) {
      ap += 1;
      t *= x / ap;
      const double u = s + t;
      cmp += std::fabs(s) >= std::fabs(t) ? (s - u) + t : (t - u) + s;
      s = u;
      if (t < eps * 1e-2 * s) break;
    }
    const double p = dpois(a, x) * (s + cmp);
    return lower ? p : 1 - p;
  }
  // x^a e^-x / Gamma(a) times the continued fraction (modified Lentz)
  double b = x + 1 - a, c = 1 / tiny, d = 1 / b, h = d;
  for (long i = 1; i < 100000000; ++i) {
    const double an = -i * (i - a);
    b += 2;
    d = an * d + b;
    if (std::fabs(d) < tiny) d = tiny;
    c = b + an / c;
    if (std::fabs(c) < tiny) c = tiny;
    d = 1 / d;
    const double del = d * c;
    h *= del;
    if (std::fabs(del - 1) < eps) break;
  }
  const double q = a * dpois(a, x) * h;
  return lower ? 1 - q : q;
}

// The quantile of an increasing function: u with h(u) = 0, h increasing,
// h' its derivative; u0 a start, [umin, umax] the range; rel: stop at a
// step relative to |u|, else absolute.
template <class H, class D>
double solve(const H &h, const D &dh, double u0, double umin, double umax, bool rel) {
  double r = h(u0);
  if (r == 0) return u0;
  double lo = umin, hi = umax, step = 1;
  if (r < 0) {
    lo = u0;
    for (;;) {
      const double v = std::fmin(u0 + step, umax);
      if (h(v) >= 0 || v == umax) {
        hi = v;
        break;
      }
      lo = v;
      step *= 2;
    }
  } else {
    hi = u0;
    for (;;) {
      const double v = std::fmax(u0 - step, umin);
      if (h(v) <= 0 || v == umin) {
        lo = v;
        break;
      }
      hi = v;
      step *= 2;
    }
  }
  double u = u0;
  for (int it = 0; it < 400; ++it) {
    r = h(u);
    if (r == 0) return u;
    if (r < 0) lo = u;
    else hi = u;
    const double g = dh(u);
    double un = u - r / g;
    if (!(un > lo && un < hi) || !std::isfinite(un)) un = 0.5 * (lo + hi);
    const double tol = 2 * eps * (rel ? std::fabs(un) : 1.0);
    if (std::fabs(un - u) <= tol || hi - lo <= tol) return un;
    u = un;
  }
  return u;
}

}  // namespace

// Wichura's algorithm AS 241 (PPND16), Applied Statistics 37 (1988),
// 477-484, doi:10.2307/2347330: relative error about 1e-16, the algorithm
// and the operation order of R's qnorm() (whose tails beyond 1e-300 are
// refined further).
double clab_qnorm(double p) {
  if (std::isnan(p) || p < 0 || p > 1) return std::numeric_limits<double>::quiet_NaN();
  if (p == 0) return -HUGE_VAL;
  if (p == 1) return HUGE_VAL;
  const double q = p - 0.5;
  double r, v;
  if (std::fabs(q) <= 0.425) {
    r = 0.180625 - q * q;
    return q *
           (((((((r * 2509.0809287301226727 + 33430.575583588128105) * r + 67265.770927008700853) * r +
                45921.953931549871457) * r + 13731.693765509461125) * r + 1971.5909503065514427) * r +
             133.14166789178437745) * r + 3.387132872796366608) /
           (((((((r * 5226.495278852545561 + 28729.085735721942674) * r + 39307.89580009271061) * r +
                21213.794301586595867) * r + 5394.1960214247511077) * r + 687.1870074920579083) * r +
             42.313330701600911252) * r + 1.0);
  }
  r = q > 0 ? 0.5 - p + 0.5 : p;
  r = std::sqrt(-std::log(r));
  if (r <= 5.0) {
    r -= 1.6;
    v = (((((((r * 7.7454501427834140764e-4 + 0.0227238449892691845833) * r + 0.24178072517745061177) * r +
             1.27045825245236838258) * r + 3.64784832476320460504) * r + 5.7694972214606914055) * r +
          4.6303378461565452959) * r + 1.42343711074968357734) /
        (((((((r * 1.05075007164441684324e-9 + 5.475938084995344946e-4) * r + 0.0151986665636164571966) * r +
             0.14810397642748007459) * r + 0.68976733498510000455) * r + 1.6763848301838038494) * r +
          2.05319162663775882187) * r + 1.0);
  } else {
    r -= 5.0;
    v = (((((((r * 2.01033439929228813265e-7 + 2.71155556874348757815e-5) * r + 0.0012426609473880784386) * r +
             0.026532189526576123093) * r + 0.29656057182850489123) * r + 1.7848265399172913358) * r +
          5.4637849111641143699) * r + 6.6579046435011037772) /
        (((((((r * 2.04426310338993978564e-15 + 1.4215117583164458887e-7) * r + 1.8463183175100546818e-5) * r +
             7.868691311456132591e-4) * r + 0.0148753612908506148525) * r + 0.13692988092273580531) * r +
          0.59983220655588793769) * r + 1.0);
  }
  return q < 0 ? -v : v;
}

// From erfc, so the lower tail keeps its relative precision far out.
double clab_pnorm(double x) { return 0.5 * std::erfc(-x * 0.70710678118654752440); }

double clab_pgamma(double x, double a, int lower) { return pgam(x, a, lower); }

double clab_qgamma(double p, double a, int lower) {
  if (std::isnan(p) || std::isnan(a) || p < 0 || p > 1 || a < 0) return qnan;
  if (a == 0) return 0;
  if (p == (lower ? 0 : 1)) return 0;
  if (p == (lower ? 1 : 0)) return HUGE_VAL;
  // Wilson and Hilferty's start, or the small-x one
  const double z = lower ? clab_qnorm(p) : -clab_qnorm(p), c = 1 / (9 * a);
  double x0 = a * std::pow(1 - c + z * std::sqrt(c), 3);
  if (!(x0 > 0) || a < 1) x0 = std::exp((std::log(lower ? p : 1 - p) + lgamma1p(a)) / a);
  if (!(x0 > 0) || !std::isfinite(x0)) x0 = a;
  const double lp = std::log(p), sg = lower ? 1 : -1;
  auto h = [&](double u) { return sg * (std::log(pgam(std::exp(u), a, lower)) - lp); };
  auto dh = [&](double u) {
    const double x = std::exp(u);
    return a * dpois(a, x) / pgam(x, a, lower);
  };
  return std::exp(solve(h, dh, std::log(x0), -745, 710, false));
}

double clab_pbeta(double x, double a, double b, int lower) { return pbet(x, 1 - x, a, b, lower); }

double clab_qbeta(double p, double a, double b, int lower) {
  if (std::isnan(p) || std::isnan(a) || std::isnan(b) || p < 0 || p > 1 || a < 0 || b < 0) return qnan;
  if (p == (lower ? 0 : 1)) return 0;
  if (p == (lower ? 1 : 0)) return 1;
  if (a == 0 || b == 0 || std::isinf(a) || std::isinf(b)) {
    if (a == 0 && b == 0) return (lower ? p < 0.5 : p > 0.5) ? 0 : 1;
    if (a == 0 || a / b == 0) return 0;
    if (b == 0 || b / a == 0) return 1;
    return 0.5;
  }
  // on the logit scale: x = 1 / (1 + e^-u), 1 - x = 1 / (1 + e^u)
  const double lp = std::log(p), sg = lower ? 1 : -1;
  auto t = [&](double u) { return pbet(1 / (1 + std::exp(-u)), 1 / (1 + std::exp(u)), a, b, lower); };
  auto h = [&](double u) { return sg * (std::log(t(u)) - lp); };
  auto dh = [&](double u) { return bfac(1 / (1 + std::exp(-u)), 1 / (1 + std::exp(u)), a, b) / t(u); };
  const double u = solve(h, dh, std::log(a / b), -700, 700, false);
  return 1 / (1 + std::exp(-u));
}

double clab_pchisq(double x, double df, int lower) { return pgam(x / 2, df / 2, lower); }

double clab_qchisq(double p, double df, int lower) { return 2 * clab_qgamma(p, df / 2, lower); }

// The Poisson mixture of central chi-squares, summed out from its mode.
double clab_pnchisq(double x, double df, double ncp, int lower) {
  if (std::isnan(x) || std::isnan(df) || std::isnan(ncp) || df < 0 || ncp < 0 || std::isinf(ncp)) return qnan;
  if (ncp == 0) return clab_pchisq(x, df, lower);
  if (x <= 0) return lower ? 0 : 1;
  if (std::isinf(x)) return lower ? 1 : 0;
  const double lam = ncp / 2;
  const double j0 = std::floor(lam);
  double s = 0;
  for (int dir = 0; dir < 2; ++dir) {
    double prev = HUGE_VAL;
    for (double j = dir ? j0 - 1 : j0; j >= 0; j += dir ? -1 : 1) {
      const double w = dpois(j, lam), t = w * clab_pchisq(x, df + 2 * j, lower);
      s += t;
      if (w == 0 || (t <= 1e-17 * s && t <= prev && std::fabs(j - j0) > 2)) break;
      prev = t;
    }
  }
  return std::fmin(1.0, s);
}

double clab_pt(double x, double n, int lower) {
  if (std::isnan(x) || std::isnan(n) || n <= 0) return qnan;
  if (std::isinf(x)) return (x < 0) == (lower != 0) ? 0 : 1;
  if (std::isinf(n)) return clab_pnorm(lower ? x : -x);
  const double nx = 1 + (x / n) * x;
  double val;
  if (nx > 1e100) {
    const double lv = -0.5 * n * (2 * std::log(std::fabs(x)) - std::log(n)) -
                      (std::lgamma(0.5 * n) + std::lgamma(0.5) - std::lgamma(0.5 * n + 0.5)) - std::log(0.5 * n);
    val = std::exp(lv);
  } else {
    val = n > x * x ? clab_pbeta(x * x / (n + x * x), 0.5, n / 2, 0) : clab_pbeta(1 / nx, n / 2, 0.5, 1);
  }
  val /= 2;
  const bool up = x <= 0 ? !lower : lower;
  return up ? 0.5 - val + 0.5 : val;
}

double clab_qt(double p, double n, int lower) {
  if (std::isnan(p) || std::isnan(n) || p < 0 || p > 1 || n <= 0) return qnan;
  if (p == 0) return lower ? -HUGE_VAL : HUGE_VAL;
  if (p == 1) return lower ? HUGE_VAL : -HUGE_VAL;
  if (p == 0.5) return 0;
  if (n > 1e20) return lower ? clab_qnorm(p) : -clab_qnorm(p);
  // q: the upper tail of |t|'s side, below 1/2
  const bool neg = lower ? p < 0.5 : p > 0.5;
  const double q = p < 0.5 ? p : 1 - p;
  const double z = -clab_qnorm(q), t0 = z + (z * z * z + z) / (4 * n);
  const double lq = std::log(q);
  auto h = [&](double u) { return lq - std::log(clab_pt(std::exp(u), n, 0)); };
  auto dh = [&](double u) {
    const double t = std::exp(u);
    const double d = std::exp(std::lgamma((n + 1) / 2) - std::lgamma(n / 2) - 0.5 * std::log(n * pi) -
                              (n + 1) / 2 * std::log1p(t * t / n));
    return d * t / clab_pt(t, n, 0);
  };
  const double t = std::exp(solve(h, dh, std::log(t0 > 0 ? t0 : 1), -745, 710, false));
  return neg ? -t : t;
}

double clab_pf(double x, double m, double n, int lower) {
  if (std::isnan(x) || std::isnan(m) || std::isnan(n) || m <= 0 || n <= 0) return qnan;
  if (x <= 0) return lower ? 0 : 1;
  if (std::isinf(x)) return lower ? 1 : 0;
  if (std::isinf(m)) {
    if (std::isinf(n)) return (x < 1) == (lower != 0) ? 0 : x == 1 ? 0.5 : 1;
    return clab_pchisq(n / x, n, !lower);
  }
  if (std::isinf(n)) return clab_pchisq(x * m, m, lower);
  if (m * x > n) return clab_pbeta(n / (n + m * x), n / 2, m / 2, !lower);
  return clab_pbeta(m * x / (n + m * x), m / 2, n / 2, lower);
}

double clab_qf(double p, double m, double n, int lower) {
  if (std::isnan(p) || std::isnan(m) || std::isnan(n) || p < 0 || p > 1 || m <= 0 || n <= 0) return qnan;
  if (p == (lower ? 0 : 1)) return 0;
  if (p == (lower ? 1 : 0)) return HUGE_VAL;
  if (m <= n && n > 4e5) {
    if (std::isinf(m)) return 1;
    return clab_qchisq(p, m, lower) / m;
  }
  if (m > 4e5) return n / clab_qchisq(p, n, !lower);
  // R's (1 / qbeta(p, n / 2, m / 2, !lower) - 1) n / m, from the side
  // of the beta where 1 - y does not cancel
  const double z = clab_qbeta(p, m / 2, n / 2, lower);
  const double r = z < 0.5 ? z / (1 - z) * (n / m) : (1 / clab_qbeta(p, n / 2, m / 2, !lower) - 1) * (n / m);
  return std::isnan(r) ? qnan : r;
}

double clab_ppois(double x, double lam, int lower) {
  if (std::isnan(x) || std::isnan(lam) || lam < 0) return qnan;
  if (x < 0) return lower ? 0 : 1;
  if (lam == 0 || std::isinf(x)) return lower ? 1 : 0;
  x = std::floor(x + 1e-7);
  return pgam(lam, x + 1, !lower);
}

// The smallest whole y whose tail reaches p, with R's fuzz of 64 units
// of rounding for left continuity.
double clab_qpois(double p, double lam, int lower) {
  if (std::isnan(p) || std::isnan(lam) || p < 0 || p > 1 || lam < 0 || std::isinf(lam)) return qnan;
  if (lam == 0) return 0;
  if (p == (lower ? 0 : 1)) return 0;
  if (p == (lower ? 1 : 0)) return HUGE_VAL;
  const double s = std::sqrt(lam), z = lower ? clab_qnorm(p) : -clab_qnorm(p);
  double y = std::fmax(0.0, std::nearbyint(lam + s * (z + (z * z - 1) / (6 * s))));
  // in the lower tail F(y) >= p; in the upper P(X > y) <= p
  const double pf = lower ? p * (1 - 64 * eps) : p * (1 + 64 * eps);
  auto ok = [&](double v) { return lower ? clab_ppois(v, lam, 1) >= pf : clab_ppois(v, lam, 0) <= pf; };
  if (ok(y)) {
    while (y > 0 && ok(y - 1)) y -= 1;
  } else {
    do y += 1;
    while (!ok(y));
  }
  return y;
}

double clab_dbinom(double x, double n, double p) {
  if (std::isnan(x) || std::isnan(n) || std::isnan(p) || p < 0 || p > 1 || n < 0 || n != std::nearbyint(n)) {
    return qnan;
  }
  if (x < 0 || x != std::nearbyint(x) || std::isinf(x)) return 0;
  return dbin(x, n, p, 1 - p);
}

double clab_pbinom(double x, double n, double p, int lower) {
  if (std::isnan(x) || std::isnan(n) || std::isnan(p) || p < 0 || p > 1 || n < 0 || n != std::nearbyint(n)) {
    return qnan;
  }
  if (x < 0) return lower ? 0 : 1;
  x = std::floor(x + 1e-7);
  if (n <= x) return lower ? 1 : 0;
  return clab_pbeta(p, x + 1, n - x, !lower);
}
