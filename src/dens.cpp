// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Kernel density estimation as stats::density() computes it (R 4.4 and
// later, old.coords = FALSE), for the command line of plotrdg() and
// plotvln(): the bandwidth rules of ?bw.nrd (nrd0, nrd, ucv, bcv, and SJ by
// solving the equation or by direct plug-in, with R's counts of pairs by
// bin distance, its golden-section and parabolic minimiser of optimize()
// and the Brent root finder of uniroot()); linear binning of the (weighted)
// values onto n points; the kernel convolved by FFT over 2n points; linear
// interpolation onto the output points. The seven kernels of density() are
// offered. The FFT is a radix-2 transform of its own, so the values agree
// with R's to rounding, not bit for bit.

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

namespace {

const double kPi = 3.14159265358979323846;
const double kSqrtPi = 1.772453850905516027298;
const double kSqrt2Pi = 2.506628274631000502416;
const double k1Sqrt2Pi = 0.398942280401432677939946059934;

using Vec = std::vector<double>;
using Cvec = std::vector<std::complex<double>>;

// R's mean() and var(): a mean refined by a second pass, sums in long
// double.
double rmean(const Vec &x) {
  long double s = 0.0L;
  for (double v : x) s += v;
  const long double n = static_cast<long double>(x.size());
  long double m = s / n;
  if (std::isfinite(static_cast<double>(m))) {
    long double t = 0.0L;
    for (double v : x) t += (v - m);
    m = m + t / n;
  }
  return static_cast<double>(m);
}
double rvar(const Vec &x) {
  const double m = rmean(x);
  long double s = 0.0L;
  for (double v : x) s += static_cast<long double>(v - m) * (v - m);
  return static_cast<double>(s / static_cast<long double>(x.size() - 1));
}
// quantile(x, p, type = 7) of the sorted values s.
double q7(const Vec &s, double p) {
  const double idx = 1.0 + static_cast<double>(s.size() - 1) * p;
  const double lo = std::floor(idx), hi = std::ceil(idx);
  const double ql = s[static_cast<std::size_t>(lo) - 1];
  const double qh = s[static_cast<std::size_t>(hi) - 1];
  if (!(idx > lo) || qh == ql) return ql;
  const double h = idx - lo;
  return (1.0 - h) * ql + h * qh;
}
double iqr(const Vec &x) {
  Vec s = x;
  std::sort(s.begin(), s.end());
  return q7(s, 0.75) - q7(s, 0.25);
}
double npow(std::size_t n, double e) { return std::pow(static_cast<double>(n), e); }

// bw.nrd0() and bw.nrd().
double nrd0(const Vec &x) {
  const double hi = std::sqrt(rvar(x));
  double lo = std::min(hi, iqr(x) / 1.34);
  if (!(lo != 0.0)) {
    lo = hi;
    if (!(lo != 0.0)) lo = std::fabs(x[0]);
    if (!(lo != 0.0)) lo = 1.0;
  }
  return 0.9 * lo * npow(x.size(), -0.2);
}
double nrd(const Vec &x) {
  const double h = iqr(x) / 1.34;
  return 1.06 * std::min(std::sqrt(rvar(x)), h) * npow(x.size(), -1.0 / 5.0);
}

// bw_pair_cnts(x, nb, n > nb / 2): the bin width d and the counts of pairs
// by bin distance, pair by pair (bw_den) or from binned counts
// (bw_den_binned). False for constant data.
struct Pairs {
  double d = 0.0;
  Vec cnt;
};
bool binned(const Vec &x, std::size_t nb, double xmin, double xmax, Pairs &p) {
  // d <- diff(r) * 1.01 / nb; xx <- trunc(abs(x) / d) * sign(x);
  // xx <- xx - min(xx) + 1; tabulate(xx, nb)
  p.d = (xmax - xmin) * 1.01 / static_cast<double>(nb);
  if (!(p.d > 0)) return false;
  Vec xx(x.size(), 0.0);
  double mn = HUGE_VAL;
  for (std::size_t i = 0; i < x.size(); ++i) {
    const double s = x[i] > 0 ? 1.0 : (x[i] < 0 ? -1.0 : 0.0);
    xx[i] = std::trunc(std::fabs(x[i]) / p.d) * s;
    mn = std::min(mn, xx[i]);
  }
  Vec w(nb, 0.0);
  for (double v : xx) {
    const double k = v - mn + 1;
    if (k >= 1 && k <= static_cast<double>(nb)) w[static_cast<std::size_t>(k) - 1] += 1.0;
  }
  for (std::size_t i = 0; i < nb; ++i) {
    p.cnt[0] += w[i] * (w[i] - 1.0);
    for (std::size_t j = 0; j < i; ++j) p.cnt[i - j] += w[i] * w[j];
  }
  p.cnt[0] *= 0.5;
  return true;
}
bool pairs(const Vec &x, Pairs &p) {
  const std::size_t nb = 1000;
  const double xmin = *std::min_element(x.begin(), x.end());
  const double xmax = *std::max_element(x.begin(), x.end());
  p.cnt.assign(nb, 0.0);
  if (static_cast<double>(x.size()) > static_cast<double>(nb) / 2.0) return binned(x, nb, xmin, xmax, p);
  const double rang = (xmax - xmin) * 1.01;
  if (rang == 0) return false;
  p.d = rang / static_cast<double>(nb);
  for (std::size_t i = 1; i < x.size(); ++i) {
    const std::int64_t ii = static_cast<std::int64_t>(x[i] / p.d);
    for (std::size_t j = 0; j < i; ++j) {
      const std::int64_t jj = static_cast<std::int64_t>(x[j] / p.d);
      const std::int64_t k = ii > jj ? ii - jj : jj - ii;
      if (k < static_cast<std::int64_t>(nb)) p.cnt[static_cast<std::size_t>(k)] += 1.0;
    }
  }
  return true;
}

// The criteria of R's bandwidths.c (cut off at a squared distance of
// 1000): 0 ucv, 1 bcv, 2 phi4, 3 phi6.
double crit(int which, double n, const Pairs &p, double h) {
  double sum = 0.0;
  for (std::size_t i = 0; i < p.cnt.size(); ++i) {
    double delta = static_cast<double>(i) * p.d / h;
    delta *= delta;
    if (delta >= 1000.0) break;
    double term = 0.0;
    if (which == 0) {
      term = std::exp(-delta / 4.) - std::sqrt(8.0) * std::exp(-delta / 2.);
    } else if (which == 1) {
      term = std::exp(-delta / 4) * (delta * delta - 12 * delta + 12);
    } else if (which == 2) {
      term = std::exp(-delta / 2.) * (delta * delta - 6. * delta + 3.);
    } else {
      term = std::exp(-delta / 2) * (delta * delta * delta - 15 * delta * delta + 45 * delta - 15);
    }
    sum += term * p.cnt[i];
  }
  if (which == 0) return (0.5 + sum / n) / (n * h * kSqrtPi);
  if (which == 1) return (1 + sum / (32.0 * n)) / (2.0 * n * h * kSqrtPi);
  if (which == 2) return (2. * sum + n * 3.) / (n * (n - 1) * std::pow(h, 5.0) * kSqrt2Pi);
  return (2. * sum - 15. * n) / (n * (n - 1) * std::pow(h, 7.0) * kSqrt2Pi);
}

// Brent_fmin() of R's optimize() (a value that is not finite counts as
// DBL_MAX, as there).
struct Fmin {
  double a, b, d = 0., e = 0., v, w, x, fv, fw, fx;
};
// One step of the minimiser: the next point u, by golden section or by a
// parabola through x, w and v.
double fstep(Fmin &s, double xm, double tol1, double t2) {
  const double c = (3. - std::sqrt(5.)) * .5;
  double p = 0., q = 0., r = 0.;
  if (std::fabs(s.e) > tol1) {
    r = (s.x - s.w) * (s.fx - s.fv);
    q = (s.x - s.v) * (s.fx - s.fw);
    p = (s.x - s.v) * q - (s.x - s.w) * r;
    q = (q - r) * 2.;
    if (q > 0.) {
      p = -p;
    } else {
      q = -q;
    }
    r = s.e;
    s.e = s.d;
  }
  if (std::fabs(p) >= std::fabs(q * .5 * r) || p <= q * (s.a - s.x) || p >= q * (s.b - s.x)) {
    s.e = s.x < xm ? s.b - s.x : s.a - s.x;
    s.d = c * s.e;
  } else {
    s.d = p / q;
    const double u = s.x + s.d;
    if (u - s.a < t2 || s.b - u < t2) s.d = s.x >= xm ? -tol1 : tol1;
  }
  if (std::fabs(s.d) >= tol1) return s.x + s.d;
  return s.d > 0. ? s.x + tol1 : s.x - tol1;
}
template <class F>
double fmin(double ax, double bx, F f, double tol) {
  const double c = (3. - std::sqrt(5.)) * .5;
  const double eps = std::sqrt(DBL_EPSILON);
  auto ff = [&](double u) {
    const double r = f(u);
    return std::isfinite(r) ? r : DBL_MAX;
  };
  Fmin s;
  s.a = ax;
  s.b = bx;
  s.v = s.a + c * (s.b - s.a);
  s.w = s.x = s.v;
  s.fx = s.fv = s.fw = ff(s.x);
  const double tol3 = tol / 3.;
  for (;;) {
    const double xm = (s.a + s.b) * .5, tol1 = eps * std::fabs(s.x) + tol3, t2 = tol1 * 2.;
    if (std::fabs(s.x - xm) <= t2 - (s.b - s.a) * .5) break;
    const double u = fstep(s, xm, tol1, t2), fu = ff(u);
    if (fu <= s.fx) {
      if (u < s.x) {
        s.b = s.x;
      } else {
        s.a = s.x;
      }
      s.v = s.w;
      s.w = s.x;
      s.x = u;
      s.fv = s.fw;
      s.fw = s.fx;
      s.fx = fu;
      continue;
    }
    if (u < s.x) {
      s.a = u;
    } else {
      s.b = u;
    }
    if (fu <= s.fw || s.w == s.x) {
      s.v = s.w;
      s.fv = s.fw;
      s.w = u;
      s.fw = fu;
    } else if (fu <= s.fv || s.v == s.x || s.v == s.w) {
      s.v = u;
      s.fv = fu;
    }
  }
  return s.x;
}

// R_zeroin2() of R's uniroot() (maxiter 1000), given f at both ends.
template <class F>
double zeroin(double a, double b, double fa, double fb, F f, double tol) {
  if (fa == 0.0) return a;
  if (fb == 0.0) return b;
  double c = a, fc = fa;
  for (int it = 0; it <= 1000; ++it) {
    const double prev = b - a;
    if (std::fabs(fc) < std::fabs(fb)) {
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
    const double tola = 2 * DBL_EPSILON * std::fabs(b) + tol / 2;
    double step = (c - b) / 2;
    if (std::fabs(step) <= tola || fb == 0.0) return b;
    if (std::fabs(prev) >= tola && std::fabs(fa) > std::fabs(fb)) {
      double p = 0., q = 0.;
      const double cb = c - b;
      if (a == c) {
        const double t1 = fb / fa;
        p = cb * t1;
        q = 1.0 - t1;
      } else {
        const double qq = fa / fc, t1 = fb / fc, t2 = fb / fa;
        p = t2 * (cb * qq * (qq - t1) - (b - a) * (t1 - 1.0));
        q = (qq - 1.0) * (t1 - 1.0) * (t2 - 1.0);
      }
      if (p > 0.0) {
        q = -q;
      } else {
        p = -p;
      }
      if (p < (0.75 * cb * q - std::fabs(tola * q) / 2) && p < std::fabs(prev * q / 2)) step = p / q;
    }
    if (std::fabs(step) < tola) step = step > 0.0 ? tola : -tola;
    a = b;
    fa = fb;
    b += step;
    fb = f(b);
    if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
      c = a;
      fc = fa;
    }
  }
  return b;
}

// bw.ucv() (ucv) and bw.bcv(); 1 for constant data, 9 when the minimum
// is at one end of the range (R's warning).
int bwcv(const Vec &x, bool ucv, double *h) {
  const double n = static_cast<double>(x.size());
  const double hmax = 1.144 * std::sqrt(rvar(x)) * npow(x.size(), -1.0 / 5.0);
  Pairs p;
  if (!pairs(x, p)) return 1;
  const double lower = 0.1 * hmax, tol = 0.1 * lower;
  *h = fmin(lower, hmax, [&](double v) { return crit(ucv ? 0 : 1, n, p, v); }, tol);
  return *h < lower + tol || *h > hmax - tol ? 9 : 0;
}

// bw.SJ(): 1 constant data, 2 too sparse to find TD, 3 too sparse to find
// alph2, 4 no solution in the range of bandwidths.
int bwsj(const Vec &x, bool dpi, double *h) {
  const double n = static_cast<double>(x.size());
  Pairs p;
  if (!pairs(x, p)) return 1;
  const double scale = std::min(std::sqrt(rvar(x)), iqr(x) / 1.349);
  const double a = 1.24 * scale * npow(x.size(), -1.0 / 7.0);
  const double b = 1.23 * scale * npow(x.size(), -1.0 / 9.0);
  const double c1 = 1 / (2 * std::sqrt(kPi) * n);
  const double td = -crit(3, n, p, b);
  if (!std::isfinite(td) || td <= 0) return 2;
  if (dpi) {
    *h = std::pow(c1 / crit(2, n, p, std::pow(2.394 / (n * td), 1.0 / 7.0)), 1.0 / 5.0);
    return 0;
  }
  const double hmax = 1.144 * scale * npow(x.size(), -1.0 / 5.0);
  const double alph2 = 1.357 * std::pow(crit(2, n, p, a) / td, 1.0 / 7.0);
  if (!std::isfinite(alph2)) return 3;
  auto fsd = [&](double v) { return std::pow(c1 / crit(2, n, p, alph2 * std::pow(v, 5.0 / 7.0)), 1.0 / 5.0) - v; };
  double lower = 0.1 * hmax, upper = hmax;
  for (int itry = 1; fsd(lower) * fsd(upper) > 0; ++itry) {
    if (itry > 99) return 4;
    if (itry % 2) {
      upper = upper * 1.2;
    } else {
      lower = lower / 1.2;
    }
  }
  // tol = 0.1 * lower is a default argument, forced after the widening.
  *h = zeroin(lower, upper, fsd(lower), fsd(upper), fsd, 0.1 * lower);
  return 0;
}

// dnorm(x, 0, sd) as R's nmath computes it.
double dnorm0(double x, double sd) {
  x = std::fabs(x / sd);
  if (x >= 2 * std::sqrt(DBL_MAX)) return 0.0;
  if (x < 5) return k1Sqrt2Pi * std::exp(-0.5 * x * x) / sd;
  if (x > std::sqrt(-2 * 0.693147180559945309417 * (DBL_MIN_EXP + 1 - DBL_MANT_DIG))) return 0.;
  const double x1 = std::ldexp(std::nearbyint(std::ldexp(x, 16)), -16);
  const double x2 = x - x1;
  return k1Sqrt2Pi / sd * (std::exp(-0.5 * x1 * x1) * std::exp((-0.5 * x2 - x1) * x2));
}

// The kernel at distance k for bandwidth bw, as density()'s switch: 0
// gaussian, 1 epanechnikov, 2 rectangular, 3 triangular, 4 biweight, 5
// cosine, 6 optcosine.
double kval(int kern, double k, double bw) {
  const double ak = std::fabs(k);
  double a = 0.0, t = 0.0;
  switch (kern) {
    case 1:
      a = bw * std::sqrt(5.0);
      return ak < a ? 3. / 4 * (1 - (ak / a) * (ak / a)) / a : 0.0;
    case 2:
      a = bw * std::sqrt(3.0);
      return ak < a ? .5 / a : 0.0;
    case 3:
      a = bw * std::sqrt(6.0);
      return ak < a ? (1 - ak / a) / a : 0.0;
    case 4:
      a = bw * std::sqrt(7.0);
      t = 1 - (ak / a) * (ak / a);
      return ak < a ? 15. / 16 * t * t / a : 0.0;
    case 5:
      a = bw / std::sqrt(1. / 3 - 2 / (kPi * kPi));
      return ak < a ? (1 + std::cos(kPi * k / a)) / (2 * a) : 0.0;
    case 6:
      a = bw / std::sqrt(1 - 8 / (kPi * kPi));
      return ak < a ? kPi / 4 * std::cos(kPi * k / (2 * a)) / a : 0.0;
    default:
      return dnorm0(k, bw);
  }
}

// seq.int(from, to, length.out = n) as R's C code fills it.
Vec seqlen(double from, double to, std::size_t n) {
  Vec s(n, from);
  if (n > 1) s[n - 1] = to;
  if (n > 2) {
    const double by = (to - from) / static_cast<double>(n - 1);
    for (std::size_t i = 1; i + 1 < n; ++i) {
      s[i] = i < n / 2 ? from + static_cast<double>(i) * by : to - static_cast<double>(n - 1 - i) * by;
    }
  }
  return s;
}

// In-place radix-2 FFT of a power-of-2 length; inv: the unnormalised
// inverse, as fft(inverse = TRUE).
void fft(Cvec &a, bool inv) {
  const std::size_t n = a.size();
  for (std::size_t i = 1, j = 0; i < n; ++i) {
    std::size_t bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) std::swap(a[i], a[j]);
  }
  for (std::size_t len = 2; len <= n; len <<= 1) {
    const double ang = 2 * kPi / static_cast<double>(len) * (inv ? 1 : -1);
    for (std::size_t i = 0; i < n; i += len) {
      for (std::size_t k = 0; k < len / 2; ++k) {
        const double t = ang * static_cast<double>(k);
        const std::complex<double> w(std::cos(t), std::sin(t));
        const std::complex<double> u = a[i + k], v = a[i + k + len / 2] * w;
        a[i + k] = u + v;
        a[i + k + len / 2] = u - v;
      }
    }
  }
}

// BinDist(): linear binning of the values onto n points of [lo, up], in a
// vector of 2n.
Cvec bindist(const Vec &v, const Vec &wt, double lo, double up, std::size_t n) {
  Cvec y(2 * n);
  const double dx = (up - lo) / static_cast<double>(n - 1);
  const std::int64_t ixmax = static_cast<std::int64_t>(n) - 2;
  for (std::size_t i = 0; i < v.size(); ++i) {
    const double pos = (v[i] - lo) / dx;
    const std::int64_t ix = static_cast<std::int64_t>(std::floor(pos));
    const double fx = pos - static_cast<double>(ix);
    if (ix >= 0 && ix <= ixmax) {
      y[static_cast<std::size_t>(ix)] += wt[i] * (1 - fx);
      y[static_cast<std::size_t>(ix) + 1] += wt[i] * fx;
    } else if (ix == -1) {
      y[0] += wt[i] * fx;
    } else if (ix == ixmax + 1) {
      y[static_cast<std::size_t>(ix)] += wt[i] * (1 - fx);
    }
  }
  return y;
}

// approx(xs, ys, t): linear, found by bisection as R's approx1(); NaN
// outside.
double approx1(const Vec &xs, const Vec &ys, double t) {
  const std::size_t n = xs.size();
  if (t < xs[0] || t > xs[n - 1]) return std::nan("");
  std::size_t i = 0, j = n - 1;
  while (i < j - 1) {
    const std::size_t ij = (i + j) / 2;
    if (t < xs[ij]) {
      j = ij;
    } else {
      i = ij;
    }
  }
  if (t == xs[j]) return ys[j];
  if (t == xs[i]) return ys[i];
  return ys[i] + (ys[j] - ys[i]) * ((t - xs[i]) / (xs[j] - xs[i]));
}

// The bandwidth of a rule (1 to 6), or 2 to 5 for R's refusals, 9 for its
// warning.
int rulebw(const Vec &v, int rule, double *bw) {
  if (rule == 1) *bw = nrd0(v);
  if (rule == 2) *bw = nrd(v);
  if (rule == 3 || rule == 4) {
    const int rc = bwcv(v, rule == 3, bw);
    return rc == 1 ? 2 : rc;
  }
  if (rule == 5 || rule == 6) {
    const int rc = bwsj(v, rule == 6, bw);
    return rc == 0 ? 0 : rc + 1;
  }
  return 0;
}

}  // namespace

int clab_dens(const double *x, const double *w, int64_t nx, int rule, double bw, double adjust, int kern, int64_t nout,
              double cut, int hasfrom, double from, int hasto, double to, double *xo, double *yo, double *st,
              int *why) CLAB_TRY {
  *why = 0;
  if (nx < 1 || nout < 1 || nout > 67108864 || rule < 0 || rule > 6 || kern < 0 || kern > 6) return CLAB_EARG;
  Vec v, wt;
  double wsum = 0.0, wfin = 0.0;
  for (std::int64_t i = 0; i < nx; ++i) {
    const double wi = w ? w[i] : 0.0;
    if (w && (!std::isfinite(wi) || wi < 0)) return CLAB_EARG;
    wsum += wi;
    if (!std::isfinite(x[i])) continue;
    v.push_back(x[i]);
    wt.push_back(wi);
    wfin += wi;
  }
  if (v.empty()) return CLAB_EARG;
  // totMass of density(): the share of finite values, or of their weight
  double tot = static_cast<double>(v.size()) / static_cast<double>(nx);
  if (w) {
    tot = v.size() == static_cast<std::size_t>(nx) ? 1.0 : wfin / wsum;
  } else {
    wt.assign(v.size(), 1.0 / static_cast<double>(v.size()));
  }
  if (rule > 0) {
    if (v.size() < 2) {
      *why = 1;
      return CLAB_EDOM;
    }
    *why = rulebw(v, rule, &bw);
    if (*why != 0 && *why != 9) return CLAB_EDOM;
  }
  if (!std::isfinite(bw)) {
    *why = 6;
    return CLAB_EDOM;
  }
  bw = adjust * bw;
  if (!(bw > 0)) {
    *why = 7;
    return CLAB_EDOM;
  }
  if (!hasfrom) from = *std::min_element(v.begin(), v.end()) - cut * bw;
  if (!hasto) to = *std::max_element(v.begin(), v.end()) + cut * bw;
  if (!std::isfinite(from) || !std::isfinite(to)) {
    *why = 8;
    return CLAB_EDOM;
  }
  std::size_t n = static_cast<std::size_t>(std::max<std::int64_t>(nout, 512));
  if (n > 512) n = static_cast<std::size_t>(1) << static_cast<int>(std::ceil(std::log2(static_cast<double>(n))));
  const double lo = from - 4 * bw, up = to + 4 * bw;
  Cvec y = bindist(v, wt, lo, up, n), k(2 * n);
  for (auto &e : y) e *= tot;
  const Vec ko = seqlen(0, static_cast<double>(2 * n - 1) / static_cast<double>(n - 1) * (up - lo), 2 * n);
  for (std::size_t i = 0; i < 2 * n; ++i) k[i] = kval(kern, i <= n ? ko[i] : -ko[2 * n - i], bw);
  fft(y, false);
  fft(k, false);
  for (std::size_t i = 0; i < 2 * n; ++i) y[i] *= std::conj(k[i]);
  fft(y, true);
  Vec dens(n);
  for (std::size_t i = 0; i < n; ++i) dens[i] = std::max(0.0, y[i].real() / static_cast<double>(2 * n));
  const Vec xords = seqlen(lo, up, n), xs = seqlen(from, to, static_cast<std::size_t>(nout));
  for (std::size_t m = 0; m < xs.size(); ++m) {
    xo[m] = xs[m];
    yo[m] = approx1(xords, dens, xs[m]);
  }
  st[0] = bw;
  st[1] = static_cast<double>(n);
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)
