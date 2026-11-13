// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The numerical helpers of qn.h: BFGS with a backtracking line search, the
// inverse of a positive definite matrix by its Cholesky factor, the
// covariance from a numerical Hessian, and the Kolmogorov-Smirnov p-value
// as R computes it. No allocation.
#include "qn.h"

#include <cmath>
#include <cstddef>

#include "kutil.h"

namespace qn {

namespace {

constexpr double pi = 3.14159265358979323846;

std::size_t ix(int i, int j, int p) {
  return static_cast<std::size_t>(i) + static_cast<std::size_t>(j) * static_cast<std::size_t>(p);
}

double dot(const double *a, const double *b, int p) {
  double s = 0;
  for (int i = 0; i < p; ++i) s += a[i] * b[i];
  return s;
}

bool fin(const double *a, int p) {
  for (int i = 0; i < p; ++i) {
    if (!std::isfinite(a[i])) return false;
  }
  return true;
}

void eye(double *h, int p, double s) {
  for (int j = 0; j < p; ++j) {
    for (int i = 0; i < p; ++i) h[ix(i, j, p)] = i == j ? s : 0.0;
  }
}

// c = a b, m x m
void mmul(const double *a, const double *b, double *c, int m) {
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < m; ++j) {
      double s = 0;
      for (int k = 0; k < m; ++k) s += a[ix(i, k, m)] * b[ix(k, j, m)];
      c[ix(i, j, m)] = s;
    }
  }
}

void mcopy(const double *a, double *b, int m) {
  for (std::size_t i = 0; i < static_cast<std::size_t>(m) * static_cast<std::size_t>(m); ++i) b[i] = a[i];
}

// P(D < d) for n points (Marsaglia, Tsang and Wang 2003,
// doi:10.18637/jss.v008.i18), the algorithm of R's ks.test()
double kexact(double d, int n, double *w) {
  const int k = static_cast<int>(n * d) + 1, m = 2 * k - 1;
  const double h = k - n * d;
  double *a = w, *r = w + static_cast<std::size_t>(m) * static_cast<std::size_t>(m),
         *tmp = r + static_cast<std::size_t>(m) * static_cast<std::size_t>(m);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < m; ++j) a[ix(i, j, m)] = i - j + 1 < 0 ? 0.0 : 1.0;
  }
  for (int i = 0; i < m; ++i) {
    a[ix(i, 0, m)] -= std::pow(h, i + 1);
    a[ix(m - 1, i, m)] -= std::pow(h, m - i);
  }
  a[ix(m - 1, 0, m)] += 2 * h - 1 > 0 ? std::pow(2 * h - 1, m) : 0.0;
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < m; ++j) {
      for (int g = 1; g <= i - j + 1; ++g) a[ix(i, j, m)] /= g;
    }
  }
  // r = a^n by squaring; its entries stay below e^n (each row of a sums
  // to less than e), so no rescaling is needed below 100 points
  eye(r, m, 1.0);
  for (int e = n;;) {
    if (e & 1) {
      mmul(r, a, tmp, m);
      mcopy(tmp, r, m);
    }
    e >>= 1;
    if (!e) break;
    mmul(a, a, tmp, m);
    mcopy(tmp, a, m);
  }
  // times n! / n^n
  return r[ix(k - 1, k - 1, m)] * std::exp(std::lgamma(n + 1.0) - n * std::log(static_cast<double>(n)));
}

// The limiting distribution of sqrt(n) D (R's pkolmogorov_two_limit with
// tolerance 1e-6), x > 0
double klimit(double x) {
  const double tol = 1e-6;
  if (x < 1) {
    const int kmax = static_cast<int>(std::sqrt(2 - std::log(tol)));
    const double z = -(pi * pi / 8) / (x * x), w = std::log(x);
    double s = 0;
    for (int k = 1; k < kmax; k += 2) s += std::exp(k * k * z - w);
    return s * std::sqrt(2 * pi);
  }
  const double z = -2 * x * x;
  double s = -1, old = 0, now = 1;
  for (int k = 1; std::fabs(old - now) > tol; ++k) {
    old = now;
    now += 2 * s * std::exp(z * k * k);
    s *= -1;
  }
  return now;
}

}  // namespace

int min(int p, double *u, fn f, void *ctx, int maxit, double gtol, double *fv, int *iter,
        double *work) {
  double *h = work, *g = h + static_cast<std::size_t>(p) * static_cast<std::size_t>(p), *d = g + p,
         *un = d + p, *gn = un + p, *hy = gn + p;
  double fx = f(u, g, ctx);
  *iter = 0;
  *fv = fx;
  if (!std::isfinite(fx) || !fin(g, p)) return -1;
  eye(h, p, 1.0);
  bool first = true;
  int it = 0, conv = 0;
  for (; it < maxit; ++it) {
    double gmax = 0;
    for (int i = 0; i < p; ++i) gmax = std::fmax(gmax, std::fabs(g[i]));
    if (gmax <= gtol * (1 + std::fabs(fx))) {
      conv = 1;
      break;
    }
    for (int i = 0; i < p; ++i) {
      d[i] = 0;
      for (int j = 0; j < p; ++j) d[i] -= h[ix(i, j, p)] * g[j];
    }
    // H stays positive definite (updates need s'y > 0): a descent direction
    const double dg = dot(d, g, p);
    double a = 1, fnew = 0;
    bool found = false;
    for (int ls = 0; ls < 60 && !found; ++ls) {
      for (int i = 0; i < p; ++i) un[i] = u[i] + a * d[i];
      fnew = f(un, gn, ctx);
      found = std::isfinite(fnew) && fin(gn, p) && fnew <= fx + 1e-4 * a * dg;
      if (!found) a *= 0.5;
    }
    if (!found) break;
    // s = a d (in d), y = gn - g (in g)
    for (int i = 0; i < p; ++i) {
      d[i] *= a;
      g[i] = gn[i] - g[i];
    }
    const double sy = dot(d, g, p), yy = dot(g, g, p);
    if (sy > 1e-12 * std::sqrt(dot(d, d, p) * yy)) {
      if (first) eye(h, p, sy / yy);
      first = false;
      for (int i = 0; i < p; ++i) {
        hy[i] = 0;
        for (int j = 0; j < p; ++j) hy[i] += h[ix(i, j, p)] * g[j];
      }
      const double yhy = dot(g, hy, p);
      for (int j = 0; j < p; ++j) {
        for (int i = 0; i < p; ++i) {
          h[ix(i, j, p)] += (sy + yhy) * d[i] * d[j] / (sy * sy) - (hy[i] * d[j] + d[i] * hy[j]) / sy;
        }
      }
    }
    for (int i = 0; i < p; ++i) {
      u[i] = un[i];
      g[i] = gn[i];
    }
    fx = fnew;
  }
  *fv = fx;
  *iter = it;
  return conv;
}

int spdinv(double *a, int p) {
  // the Cholesky factor L (lower triangle)
  for (int j = 0; j < p; ++j) {
    double s = a[ix(j, j, p)];
    for (int k = 0; k < j; ++k) s -= a[ix(j, k, p)] * a[ix(j, k, p)];
    if (!(s > 0) || !std::isfinite(s)) return 1;
    a[ix(j, j, p)] = std::sqrt(s);
    for (int i = j + 1; i < p; ++i) {
      double r = a[ix(i, j, p)];
      for (int k = 0; k < j; ++k) r -= a[ix(i, k, p)] * a[ix(j, k, p)];
      a[ix(i, j, p)] = r / a[ix(j, j, p)];
    }
  }
  // L^-1, a column at a time from the left
  for (int j = 0; j < p; ++j) {
    a[ix(j, j, p)] = 1 / a[ix(j, j, p)];
    for (int i = j + 1; i < p; ++i) {
      double s = 0;
      for (int k = j; k < i; ++k) s += a[ix(i, k, p)] * a[ix(k, j, p)];
      a[ix(i, j, p)] = -s / a[ix(i, i, p)];
    }
  }
  // A^-1 = L^-T L^-1: the upper triangle first, then the diagonal
  for (int j = 0; j < p; ++j) {
    for (int i = 0; i < j; ++i) {
      double s = 0;
      for (int k = j; k < p; ++k) s += a[ix(k, i, p)] * a[ix(k, j, p)];
      a[ix(i, j, p)] = s;
    }
  }
  for (int j = 0; j < p; ++j) {
    double s = 0;
    for (int k = j; k < p; ++k) s += a[ix(k, j, p)] * a[ix(k, j, p)];
    a[ix(j, j, p)] = s;
  }
  for (int j = 0; j < p; ++j) {
    for (int i = j + 1; i < p; ++i) a[ix(i, j, p)] = a[ix(j, i, p)];
  }
  return 0;
}

int hcov(int p, double *th, fn gr, void *ctx, double *cov, double *work) {
  double *gp = work, *gm = work + p;
  bool ok = true;
  for (int k = 0; k < p; ++k) {
    const double t0 = th[k], h = 1e-5 * std::fmax(1.0, std::fabs(t0));
    th[k] = t0 + h;
    ok = std::isfinite(gr(th, gp, ctx)) && fin(gp, p) && ok;
    th[k] = t0 - h;
    ok = std::isfinite(gr(th, gm, ctx)) && fin(gm, p) && ok;
    th[k] = t0;
    for (int i = 0; i < p; ++i) cov[ix(i, k, p)] = (gp[i] - gm[i]) / (2 * h);
  }
  for (int j = 0; j < p; ++j) {
    for (int i = 0; i < j; ++i) {
      const double s = (cov[ix(i, j, p)] + cov[ix(j, i, p)]) / 2;
      cov[ix(i, j, p)] = s;
      cov[ix(j, i, p)] = s;
    }
  }
  if (ok && spdinv(cov, p) == 0) return 0;
  for (std::size_t i = 0; i < static_cast<std::size_t>(p) * static_cast<std::size_t>(p); ++i) cov[i] = kutil::nan();
  return 1;
}

double ksd(double *u, int n) {
  kutil::hsort(u, n, [](double a, double b) { return a < b; });
  double d = 0;
  for (int i = 0; i < n; ++i) {
    d = std::fmax(d, std::fmax(static_cast<double>(i + 1) / n - u[i], u[i] - static_cast<double>(i) / n));
  }
  return d;
}

double kspv(double d, int n, bool exact, double *work) {
  const double c = exact && n < 100 ? kexact(d, n, work) : klimit(std::sqrt(static_cast<double>(n)) * d);
  return std::fmin(1.0, std::fmax(0.0, 1 - c));
}

}  // namespace qn
