// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Candidate distributions for one variable fitted by maximum likelihood
// (dstfit()): the closed forms, and for the others R's own minimisers as
// optim() and optimize() run them -- the Nelder-Mead simplex of nmmin() and
// the variable-metric BFGS of vmmin() with optim()'s central differences
// (step 1e-3), and Brent's fmin() -- so the estimates follow the reference
// implementation's path. Log densities and distribution functions of the
// twelve families, the goodness-of-fit statistics and Hill's tail index.
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <functional>
#include <vector>

#include "crimlab.h"

namespace {

typedef long double ld;
using std::vector;
const double lnsqrt2pi = 0.918938533204672741780329736406;
const double pi = 3.14159265358979323846;

double rmean(const vector<double> &x) {
  const size_t n = x.size();
  if (n == 0) return NAN;
  ld s = 0;
  for (double v : x) s += v;
  s /= n;
  if (std::isfinite(static_cast<double>(s))) {
    ld t = 0;
    for (double v : x) t += v - s;
    s += t / n;
  }
  return static_cast<double>(s);
}

double rsum(const vector<double> &x) {
  ld s = 0;
  for (double v : x) s += v;
  return static_cast<double>(s);
}

double rmedian(vector<double> v) {
  std::sort(v.begin(), v.end());
  const size_t n = v.size(), h = (n + 1) / 2;
  if (n % 2) return v[h - 1];
  return rmean({v[h - 1], v[h]});
}

double lg(double p) { return std::log(p / (1 - p)); }
double il(double z) { return 1 / (1 + std::exp(-z)); }

double dpoisl(double x, double l) {
  if (l == 0) return x == 0 ? 0 : -INFINITY;
  if (x < 0 || x != std::floor(x)) return -INFINITY;
  return x * std::log(l) - l - std::lgamma(x + 1);
}

double dnbinoml(double x, double size, double mu) {
  if (x < 0 || x != std::floor(x)) return -INFINITY;
  if (x == 0) return size * (size < mu ? std::log(size / (size + mu)) : std::log1p(-mu / (size + mu)));
  return std::lgamma(x + size) - std::lgamma(size) - std::lgamma(x + 1) + size * std::log(size / (size + mu)) +
         x * std::log(mu / (size + mu));
}

double pnbin(double q, double size, double mu) {
  const double x = std::floor(q + 1e-7);
  if (x < 0) return 0;
  return clab_pbeta(size / (size + mu), size, x + 1, 1);
}

// The log density of family f at x.
double ld1(int f, double x, const double *p) {
  switch (f) {
    case 0: {
      const double z = (x - p[0]) / p[1];
      return -(lnsqrt2pi + 0.5 * z * z + std::log(p[1]));
    }
    case 1: {
      if (x <= 0) return -INFINITY;
      const double z = (std::log(x) - p[0]) / p[1];
      return -(lnsqrt2pi + 0.5 * z * z + std::log(x * p[1]));
    }
    case 2: {
      const double sc = 1 / p[0];
      return x < 0 ? -INFINITY : (-x / sc) - std::log(sc);
    }
    case 3: {
      if (x < 0) return -INFINITY;
      const double a = p[0], sc = 1 / p[1];
      if (x == 0) return a < 1 ? INFINITY : a > 1 ? -INFINITY : -std::log(sc);
      // dpois_raw(a - 1, x / scale) - log(scale), as the gamma density
      const double l = x / sc;
      if (a < 1) return a * std::log(l) - l - std::lgamma(a + 1) + std::log(a / x);
      const double k = a - 1;
      return k * std::log(l) - l - std::lgamma(k + 1) - std::log(sc);
    }
    case 4: {
      if (x < 0) return -INFINITY;
      const double t1 = std::pow(x / p[1], p[0] - 1), t2 = t1 * (x / p[1]);
      return -t2 + std::log(p[0] * t1 / p[1]);
    }
    case 5: {
      const double z = (x - p[0]) / p[1], n = p[2];
      return std::lgamma((n + 1) / 2) - std::lgamma(n / 2) - 0.5 * std::log(n * pi) -
             (n + 1) / 2 * std::log1p(z * z / n) - std::log(p[1]);
    }
    case 6:
      return x >= p[0] ? std::log(p[1]) + p[1] * std::log(p[0]) - (p[1] + 1) * std::log(x) : -INFINITY;
    case 7: {
      if (x < 0 || x > 1) return -INFINITY;
      const double a = p[0], b = p[1];
      return (a - 1) * std::log(x) + (b - 1) * std::log1p(-x) -
             (std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b));
    }
    case 8:
      return dpoisl(x, p[0]);
    case 9:
      return dnbinoml(x, p[0], p[1]);
    case 10:
      return x == 0 ? std::log(p[0] + (1 - p[0]) * std::exp(-p[1])) : std::log(1 - p[0]) + dpoisl(x, p[1]);
    default:
      return x == 0 ? std::log(p[0] + (1 - p[0]) * std::exp(dnbinoml(0, p[1], p[2])))
                    : std::log(1 - p[0]) + dnbinoml(x, p[1], p[2]);
  }
}

// The distribution function of family f at q.
double pf1(int f, double q, const double *p) {
  switch (f) {
    case 0:
      return clab_pnorm((q - p[0]) / p[1]);
    case 1:
      return q <= 0 ? 0 : clab_pnorm((std::log(q) - p[0]) / p[1]);
    case 2:
      return q <= 0 ? 0 : -std::expm1(-q / (1 / p[0]));
    case 3:
      return clab_pgamma(q / (1 / p[1]), p[0], 1);
    case 4:
      return q <= 0 ? 0 : -std::expm1(-std::pow(q / p[1], p[0]));
    case 5:
      return clab_pt((q - p[0]) / p[1], p[2], 1);
    case 6:
      return q < p[0] ? 0 : 1 - std::pow(p[0] / q, p[1]);
    case 7:
      return clab_pbeta(q, p[0], p[1], 1);
    case 8:
      return clab_ppois(q, p[0], 1);
    case 9:
      return pnbin(q, p[0], p[1]);
    case 10:
      return q < 0 ? 0 : p[0] + (1 - p[0]) * clab_ppois(q, p[1], 1);
    default:
      return q < 0 ? 0 : p[0] + (1 - p[0]) * pnbin(q, p[1], p[2]);
  }
}

typedef std::function<double(const double *)> Fn;

// optim(method = "Nelder-Mead") of R: nmmin() with alpha 1, beta 0.5,
// gamma 2, abstol -Inf. false when f is not finite at the start.
bool nmmin(int n, vector<double> &b, const Fn &fn, double intol, int maxit) {
  const double big = 1.0e+35, alpha = 1, bet = 0.5, gamm = 2;
  vector<vector<double>> P(static_cast<size_t>(n) + 1, vector<double>(static_cast<size_t>(n) + 2));
  double f = fn(b.data());
  if (!std::isfinite(f)) return false;
  int funcount = 1;
  const double convtol = intol * (std::fabs(f) + intol);
  const int n1 = n + 1, C = n + 2;
  P[static_cast<size_t>(n1 - 1)][0] = f;
  for (int i = 0; i < n; ++i) P[static_cast<size_t>(i)][0] = b[static_cast<size_t>(i)];
  int L = 1, H;
  double size = 0, step = 0;
  for (int i = 0; i < n; ++i) step = std::max(step, 0.1 * std::fabs(b[static_cast<size_t>(i)]));
  if (step == 0) step = 0.1;
  for (int j = 2; j <= n1; ++j) {
    for (int i = 0; i < n; ++i) P[static_cast<size_t>(i)][static_cast<size_t>(j - 1)] = b[static_cast<size_t>(i)];
    double trystep = step;
    while (P[static_cast<size_t>(j - 2)][static_cast<size_t>(j - 1)] == b[static_cast<size_t>(j - 2)]) {
      P[static_cast<size_t>(j - 2)][static_cast<size_t>(j - 1)] = b[static_cast<size_t>(j - 2)] + trystep;
      trystep *= 10;
    }
    size += trystep;
  }
  double oldsize = size;
  bool calcvert = true;
  auto at = [&](int i, int j) -> double & { return P[static_cast<size_t>(i)][static_cast<size_t>(j)]; };
  do {
    if (calcvert) {
      for (int j = 0; j < n1; ++j) {
        if (j + 1 != L) {
          for (int i = 0; i < n; ++i) b[static_cast<size_t>(i)] = at(i, j);
          f = fn(b.data());
          if (!std::isfinite(f)) f = big;
          ++funcount;
          at(n1 - 1, j) = f;
        }
      }
      calcvert = false;
    }
    double VL = at(n1 - 1, L - 1), VH = VL;
    H = L;
    for (int j = 1; j <= n1; ++j) {
      if (j != L) {
        f = at(n1 - 1, j - 1);
        if (f < VL) {
          L = j;
          VL = f;
        }
        if (f > VH) {
          H = j;
          VH = f;
        }
      }
    }
    if (VH <= VL + convtol || VL <= -INFINITY) break;
    for (int i = 0; i < n; ++i) {
      double t = -at(i, H - 1);
      for (int j = 0; j < n1; ++j) t += at(i, j);
      at(i, C - 1) = t / n;
    }
    for (int i = 0; i < n; ++i) b[static_cast<size_t>(i)] = (1.0 + alpha) * at(i, C - 1) - alpha * at(i, H - 1);
    f = fn(b.data());
    if (!std::isfinite(f)) f = big;
    ++funcount;
    const double VR = f;
    if (VR < VL) {
      at(n1 - 1, C - 1) = f;
      for (int i = 0; i < n; ++i) {
        f = gamm * b[static_cast<size_t>(i)] + (1 - gamm) * at(i, C - 1);
        at(i, C - 1) = b[static_cast<size_t>(i)];
        b[static_cast<size_t>(i)] = f;
      }
      f = fn(b.data());
      if (!std::isfinite(f)) f = big;
      ++funcount;
      if (f < VR) {
        for (int i = 0; i < n; ++i) at(i, H - 1) = b[static_cast<size_t>(i)];
        at(n1 - 1, H - 1) = f;
      } else {
        for (int i = 0; i < n; ++i) at(i, H - 1) = at(i, C - 1);
        at(n1 - 1, H - 1) = VR;
      }
    } else {
      if (VR < VH) {
        for (int i = 0; i < n; ++i) at(i, H - 1) = b[static_cast<size_t>(i)];
        at(n1 - 1, H - 1) = VR;
      }
      for (int i = 0; i < n; ++i) b[static_cast<size_t>(i)] = (1 - bet) * at(i, H - 1) + bet * at(i, C - 1);
      f = fn(b.data());
      if (!std::isfinite(f)) f = big;
      ++funcount;
      if (f < at(n1 - 1, H - 1)) {
        for (int i = 0; i < n; ++i) at(i, H - 1) = b[static_cast<size_t>(i)];
        at(n1 - 1, H - 1) = f;
      } else if (VR >= VH) {
        calcvert = true;
        size = 0;
        for (int j = 0; j < n1; ++j) {
          if (j + 1 != L) {
            for (int i = 0; i < n; ++i) {
              at(i, j) = bet * (at(i, j) - at(i, L - 1)) + at(i, L - 1);
              size += std::fabs(at(i, j) - at(i, L - 1));
            }
          }
        }
        if (size < oldsize) oldsize = size;
        else break;
      }
    }
  } while (funcount <= maxit);
  for (int i = 0; i < n; ++i) b[static_cast<size_t>(i)] = at(i, L - 1);
  return true;
}

// optim()'s central-difference gradient (ndeps 1e-3); false when a
// difference is not finite (optim() stops with an error there).
bool numgr(int n, vector<double> &p, const Fn &fn, vector<double> &g) {
  for (int i = 0; i < n; ++i) {
    const double e = 1e-3, p0 = p[static_cast<size_t>(i)];
    p[static_cast<size_t>(i)] = p0 + e;
    const double v1 = fn(p.data());
    p[static_cast<size_t>(i)] = p0 - e;
    const double v2 = fn(p.data());
    g[static_cast<size_t>(i)] = (v1 - v2) / (2 * e);
    p[static_cast<size_t>(i)] = p0;
    if (!std::isfinite(g[static_cast<size_t>(i)])) return false;
  }
  return true;
}

// optim(method = "BFGS") of R: vmmin() with abstol -Inf.
bool vmmin(int n, vector<double> &b, const Fn &fn, double reltol, int maxit) {
  const double stepredn = 0.2, acctol = 0.0001, reltest = 10.0;
  vector<double> g(static_cast<size_t>(n)), t(static_cast<size_t>(n)), X(static_cast<size_t>(n)),
      c(static_cast<size_t>(n));
  vector<vector<double>> B(static_cast<size_t>(n), vector<double>(static_cast<size_t>(n)));
  double f = fn(b.data());
  if (!std::isfinite(f)) return false;
  double Fmin = f;
  int funcount = 1, gradcount = 1, iter = 0, count = 0;
  if (!numgr(n, b, fn, g)) return false;
  ++iter;
  int ilast = gradcount;
  do {
    if (ilast == gradcount) {
      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) B[static_cast<size_t>(i)][static_cast<size_t>(j)] = 0.0;
        B[static_cast<size_t>(i)][static_cast<size_t>(i)] = 1.0;
      }
    }
    for (int i = 0; i < n; ++i) {
      X[static_cast<size_t>(i)] = b[static_cast<size_t>(i)];
      c[static_cast<size_t>(i)] = g[static_cast<size_t>(i)];
    }
    double gradproj = 0.0;
    for (int i = 0; i < n; ++i) {
      double s = 0.0;
      for (int j = 0; j <= i; ++j) s -= B[static_cast<size_t>(i)][static_cast<size_t>(j)] * g[static_cast<size_t>(j)];
      for (int j = i + 1; j < n; ++j) s -= B[static_cast<size_t>(j)][static_cast<size_t>(i)] * g[static_cast<size_t>(j)];
      t[static_cast<size_t>(i)] = s;
      gradproj += s * g[static_cast<size_t>(i)];
    }
    if (gradproj < 0.0) {
      double steplength = 1.0;
      bool accpoint = false;
      do {
        count = 0;
        for (int i = 0; i < n; ++i) {
          b[static_cast<size_t>(i)] = X[static_cast<size_t>(i)] + steplength * t[static_cast<size_t>(i)];
          if (reltest + X[static_cast<size_t>(i)] == reltest + b[static_cast<size_t>(i)]) ++count;
        }
        if (count < n) {
          f = fn(b.data());
          ++funcount;
          accpoint = std::isfinite(f) && (f <= Fmin + gradproj * steplength * acctol);
          if (!accpoint) steplength *= stepredn;
        }
      } while (!(count == n || accpoint));
      const bool enough = std::fabs(f - Fmin) > reltol * (std::fabs(Fmin) + reltol);
      if (!enough) {
        count = n;
        Fmin = f;
      }
      if (count < n) {
        Fmin = f;
        if (!numgr(n, b, fn, g)) return false;
        ++gradcount;
        ++iter;
        double D1 = 0.0;
        for (int i = 0; i < n; ++i) {
          t[static_cast<size_t>(i)] = steplength * t[static_cast<size_t>(i)];
          c[static_cast<size_t>(i)] = g[static_cast<size_t>(i)] - c[static_cast<size_t>(i)];
          D1 += t[static_cast<size_t>(i)] * c[static_cast<size_t>(i)];
        }
        if (D1 > 0) {
          double D2 = 0.0;
          for (int i = 0; i < n; ++i) {
            double s = 0.0;
            for (int j = 0; j <= i; ++j) s += B[static_cast<size_t>(i)][static_cast<size_t>(j)] * c[static_cast<size_t>(j)];
            for (int j = i + 1; j < n; ++j) s += B[static_cast<size_t>(j)][static_cast<size_t>(i)] * c[static_cast<size_t>(j)];
            X[static_cast<size_t>(i)] = s;
            D2 += s * c[static_cast<size_t>(i)];
          }
          D2 = 1.0 + D2 / D1;
          for (int i = 0; i < n; ++i) {
            for (int j = 0; j <= i; ++j) {
              B[static_cast<size_t>(i)][static_cast<size_t>(j)] +=
                  (D2 * t[static_cast<size_t>(i)] * t[static_cast<size_t>(j)] -
                   X[static_cast<size_t>(i)] * c[static_cast<size_t>(j)] -
                   t[static_cast<size_t>(i)] * X[static_cast<size_t>(j)]) /
                  D1;
            }
          }
        } else {
          ilast = gradcount;
        }
      } else if (ilast < gradcount) {
        count = 0;
        ilast = gradcount;
      }
    } else {
      count = 0;
      if (ilast == gradcount) count = n;
      else ilast = gradcount;
    }
    if (iter >= maxit) break;
    if (gradcount - ilast > 2 * n) ilast = gradcount;
  } while (count != n || ilast != gradcount);
  (void)funcount;
  return true;
}

// optimize() of R: Brent's fmin() on [ax, bx], tolerance eps^0.25; a value
// that is not finite counts as the largest double.
double fmin(double ax, double bx, const std::function<double(double)> &fn) {
  const double c = (3. - std::sqrt(5.)) * .5, tol = std::pow(DBL_EPSILON, 0.25);
  auto f = [&](double u) {
    const double v = fn(u);
    return std::isfinite(v) ? v : DBL_MAX;
  };
  double eps = std::sqrt(DBL_EPSILON), a = ax, b = bx, v = a + c * (b - a), w = v, x = v, d = 0., e = 0.;
  double fx = f(x), fv = fx, fw = fx, tol3 = tol / 3.;
  for (;;) {
    const double xm = (a + b) * .5, tol1 = eps * std::fabs(x) + tol3, t2 = tol1 * 2.;
    if (std::fabs(x - xm) <= t2 - (b - a) * .5) break;
    double p = 0., q = 0., r = 0., u;
    if (std::fabs(e) > tol1) {
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = (q - r) * 2.;
      if (q > 0.) p = -p;
      else q = -q;
      r = e;
      e = d;
    }
    if (std::fabs(p) >= std::fabs(q * .5 * r) || p <= q * (a - x) || p >= q * (b - x)) {
      e = x < xm ? b - x : a - x;
      d = c * e;
    } else {
      d = p / q;
      u = x + d;
      if (u - a < t2 || b - u < t2) {
        d = tol1;
        if (x >= xm) d = -d;
      }
    }
    if (std::fabs(d) >= tol1) u = x + d;
    else if (d > 0.) u = x + tol1;
    else u = x - tol1;
    const double fu = f(u);
    if (fu <= fx) {
      if (u < x) b = x;
      else a = x;
      v = w;
      w = x;
      x = u;
      fv = fw;
      fw = fx;
      fx = fu;
    } else {
      if (u < x) a = u;
      else b = u;
      if (fu <= fw || w == x) {
        v = w;
        fv = fw;
        w = u;
        fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v = u;
        fv = fu;
      }
    }
  }
  return x;
}

// .dst_ml(): Nelder-Mead (maxit 4000, reltol 1e-12), then BFGS (maxit 500,
// reltol 1e-14), on the transformed scale z.
bool ml(vector<double> &z, const Fn &nll) {
  const int n = static_cast<int>(z.size());
  if (!nmmin(n, z, nll, 1e-12, 4000)) return false;
  return vmmin(n, z, nll, 1e-14, 500);
}

double nllsum(int f, const vector<double> &x, const double *p) {
  ld s = 0;
  for (double v : x) s += ld1(f, v, p);
  return static_cast<double>(-s);
}

double p0(const vector<double> &x) {
  double z = 0;
  for (double v : x) z += v == 0;
  const double m = rmean(x);
  return std::min(std::max(z / static_cast<double>(x.size()) - std::exp(-m), 0.05), 0.9);
}

double posmean(const vector<double> &x) {
  vector<double> y;
  for (double v : x) {
    if (v > 0) y.push_back(v);
  }
  return std::max(rmean(y), 0.1);
}

}  // namespace

extern "C" {

int clab_dfnp(int f) {
  static const int np[12] = {2, 2, 1, 2, 2, 3, 2, 2, 1, 2, 2, 3};
  return f >= 0 && f < 12 ? np[f] : 0;
}

int clab_dfit(int f, const double *x0, int n, double *par) {
  if (f < 0 || f > 11 || n < 1) return CLAB_EARG;
  const vector<double> x(x0, x0 + n);
  const int k = clab_dfnp(f);
  vector<double> z;
  bool ok = true;
  switch (f) {
    case 0:
    case 1: {
      vector<double> y = x;
      if (f == 1) {
        for (double &v : y) v = std::log(v);
      }
      const double m = rmean(y);
      vector<double> d;
      for (double v : y) d.push_back((v - m) * (v - m));
      par[0] = m;
      par[1] = std::sqrt(rmean(d));
      break;
    }
    case 2:
      par[0] = 1 / rmean(x);
      break;
    case 3:
    case 7: {
      const double m = rmean(x);
      vector<double> d;
      for (double v : x) d.push_back((v - m) * (v - m));
      const double v = rmean(d);
      if (f == 3) {
        z = {std::log(m * m / v), std::log(m / v)};
      } else {
        const double kk = std::max(m * (1 - m) / v - 1, 0.1);
        z = {std::log(m * kk), std::log((1 - m) * kk)};
      }
      ok = ml(z, [&](const double *u) {
        const double p[2] = {std::exp(u[0]), std::exp(u[1])};
        return nllsum(f, x, p);
      });
      par[0] = std::exp(z[0]);
      par[1] = std::exp(z[1]);
      break;
    }
    case 4:
      z = {std::log(1.2), std::log(rmean(x))};
      ok = ml(z, [&](const double *u) {
        const double p[2] = {std::exp(u[0]), std::exp(u[1])};
        return nllsum(4, x, p);
      });
      par[0] = std::exp(z[0]);
      par[1] = std::exp(z[1]);
      break;
    case 5: {
      const double med = rmedian(x);
      vector<double> a;
      for (double v : x) a.push_back(std::fabs(v - med));
      const double madv = 1.4826 * rmedian(a), m = rmean(x);
      ld s = 0;
      for (double v : x) {
        const double e = v - m;
        s += e * e;
      }
      const double sd = std::sqrt(static_cast<double>(s / (n - 1)));
      z = {med, std::log(std::max(madv, sd / 10)), std::log(5.0)};
      ok = ml(z, [&](const double *u) {
        ld t = 0;
        const double sc = std::exp(u[1]), df = std::exp(u[2]), p[3] = {0, 1, df};
        for (double v : x) t += ld1(5, (v - u[0]) / sc, p) - u[1];
        return static_cast<double>(-t);
      });
      par[0] = z[0];
      par[1] = std::exp(z[1]);
      par[2] = std::exp(z[2]);
      break;
    }
    case 6: {
      const double mn = *std::min_element(x.begin(), x.end());
      vector<double> l;
      for (double v : x) l.push_back(std::log(v / mn));
      par[0] = mn;
      par[1] = n / rsum(l);
      break;
    }
    case 8:
      par[0] = rmean(x);
      break;
    case 9: {
      const double m = rmean(x);
      const double zz = fmin(-10, 15, [&](double u) {
        const double p[2] = {std::exp(u), m};
        return nllsum(9, x, p);
      });
      par[0] = std::exp(zz);
      par[1] = m;
      break;
    }
    case 10:
      z = {lg(p0(x)), std::log(posmean(x))};
      ok = ml(z, [&](const double *u) {
        const double p[2] = {il(u[0]), std::exp(u[1])};
        return nllsum(10, x, p);
      });
      par[0] = il(z[0]);
      par[1] = std::exp(z[1]);
      break;
    default:
      z = {lg(p0(x)), 0, std::log(posmean(x))};
      ok = ml(z, [&](const double *u) {
        const double p[3] = {il(u[0]), std::exp(u[1]), std::exp(u[2])};
        return nllsum(11, x, p);
      });
      par[0] = il(z[0]);
      par[1] = std::exp(z[1]);
      par[2] = std::exp(z[2]);
  }
  if (!ok) return CLAB_EDOM;
  for (int i = 0; i < k; ++i) {
    if (!std::isfinite(par[i])) return CLAB_EDOM;
  }
  if (!std::isfinite(-nllsum(f, x, par))) return CLAB_EDOM;
  return CLAB_OK;
}

double clab_dfll(int f, const double *x, int n, const double *par) {
  return -nllsum(f, vector<double>(x, x + n), par);
}

void clab_dfgof(int f, const double *xs, int n, const double *par, double *out) {
  if (f >= 8) {
    double mx = 0, ks = 0;
    for (int i = 0; i < n; ++i) mx = std::max(mx, xs[i]);
    for (double k = 0; k <= mx; ++k) {
      double c = 0;
      for (int i = 0; i < n; ++i) c += xs[i] <= k;
      ks = std::max(ks, std::fabs(c / n - pf1(f, k, par)));
    }
    out[0] = ks;
    out[1] = NAN;
    return;
  }
  vector<double> u(static_cast<size_t>(n));
  for (int i = 0; i < n; ++i) u[static_cast<size_t>(i)] = std::min(std::max(pf1(f, xs[i], par), 1e-300), 1 - 1e-16);
  double ks = -INFINITY;
  vector<double> t(static_cast<size_t>(n));
  for (int i = 0; i < n; ++i) {
    const double ui = u[static_cast<size_t>(i)];
    ks = std::max(ks, std::max((i + 1.0) / n - ui, ui - static_cast<double>(i) / n));
    t[static_cast<size_t>(i)] = (2.0 * (i + 1) - 1) * (std::log(ui) + std::log(1 - u[static_cast<size_t>(n - 1 - i)]));
  }
  out[0] = ks;
  out[1] = -n - rmean(t);
}

}  // extern "C"
