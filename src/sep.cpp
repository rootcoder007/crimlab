// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, spatial epidemiology and surveillance: the Bayes outbreak
// detector, exposure within a buffer and by kernel, ecological Poisson,
// negative binomial and zero-inflated Poisson regression with expected
// counts, the Leroux and BYM2 area structures, eigenvector spatial
// filtering, the metapopulation SIR model, Weibull frailty and
// promotion-time cure models with intrinsic CAR area effects, and the
// structural bias of graph attention.
#include <cstdint>
#include <algorithm>
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "gla.h"
#include "kutil.h"

namespace {

using vd = std::vector<double>;
using gla::at;

constexpr double pi = 3.14159265358979323846;

// The dense n x n matrix of a weights structure.
bool dense(const int *rp, const int *ci, const double *w, int n, vd &A) {
  A.assign(static_cast<std::size_t>(n) * n, 0.0);
  for (int i = 0; i < n; ++i)
    for (int q = rp[i]; q < rp[i + 1]; ++q) {
      const int j = ci[q];
      if (j < 0 || j >= n || !kutil::nonneg(w[q])) return false;
      A[at(i, j, n)] += w[q];
    }
  return true;
}

// The Laplacian diag(rowSums(A)) - A.
vd lap(const vd &A, int n) {
  vd Q(A.size());
  for (int i = 0; i < n; ++i) {
    double s = 0;
    for (int j = 0; j < n; ++j) s += A[at(i, j, n)];
    for (int j = 0; j < n; ++j) Q[at(i, j, n)] = (i == j ? s : 0.0) - A[at(i, j, n)];
  }
  return Q;
}

// digamma and trigamma for x > 0 (NaN otherwise): the recurrence up to
// 12, then the asymptotic series
double dig(double x) {
  if (!(x > 0)) return kutil::nan();
  double acc = 0;
  while (x < 12) {
    acc -= 1 / x;
    x += 1;
  }
  const double inv = 1 / x, i2 = inv * inv;
  return acc + std::log(x) - 0.5 * inv - i2 * (1.0 / 12 - i2 * (1.0 / 120 - i2 * (1.0 / 252 - i2 / 240)));
}
double trig(double x) {
  if (!(x > 0)) return kutil::nan();
  double acc = 0;
  while (x < 12) {
    acc += 1 / (x * x);
    x += 1;
  }
  const double inv = 1 / x, i2 = inv * inv;
  return acc + inv + i2 / 2 + inv * i2 * (1.0 / 6 - i2 * (1.0 / 30 - i2 * (1.0 / 42 - i2 / 30)));
}

// a p x p system a b = r by LU (a and r copied)
bool solve(const vd &a, const vd &r, int p, vd &b) {
  vd m(a);
  b = r;
  std::vector<int> piv(static_cast<std::size_t>(p));
  return gla::lusolve(m.data(), p, b.data(), 1, piv.data());
}

// The design: an intercept, then the p columns of x (n rows).
vd design(const double *x, int n, int p) {
  vd X(static_cast<std::size_t>(n) * (p + 1));
  for (int i = 0; i < n; ++i) X[i] = 1;
  if (p > 0) std::copy(x, x + static_cast<std::size_t>(n) * p, X.begin() + n);
  return X;
}

// weighted cross-products X'WX and X'Wv (w NULL: unit weights)
void xwx(const vd &X, int n, int k, const double *w, const double *v, vd &a, vd &r) {
  a.assign(static_cast<std::size_t>(k) * k, 0.0);
  r.assign(static_cast<std::size_t>(k), 0.0);
  for (int j = 0; j < k; ++j) {
    for (int l = 0; l < k; ++l) {
      double s = 0;
      for (int i = 0; i < n; ++i) s += X[at(i, j, n)] * (w ? w[i] : 1.0) * X[at(i, l, n)];
      a[at(j, l, k)] = s;
    }
    double s = 0;
    for (int i = 0; i < n; ++i) s += X[at(i, j, n)] * (w ? w[i] : 1.0) * v[i];
    r[j] = s;
  }
}

void lin(const vd &X, int n, int k, const double *b, const double *off, double *eta) {
  for (int i = 0; i < n; ++i) {
    double s = 0;
    for (int j = 0; j < k; ++j) s += X[at(i, j, n)] * b[j];
    eta[i] = s + (off ? off[i] : 0.0);
  }
}

// Poisson-type IRLS with working weights wf(mu): 100 steps at most, to a
// relative change of 1e-12.
template <class W>
bool irls(const double *y, const vd &X, int n, int k, const double *off, W wf, vd &beta) {
  beta.assign(static_cast<std::size_t>(k), 0.0);
  double sy = 0, se = 0;
  for (int i = 0; i < n; ++i) {
    sy += y[i];
    se += std::exp(off[i]);
  }
  beta[0] = std::log(std::max(sy, 0.5) / se);
  vd eta(static_cast<std::size_t>(n)), w(static_cast<std::size_t>(n)), z(static_cast<std::size_t>(n)), a, r, nb;
  for (int it = 0; it < 100; ++it) {
    if (kutil::stop(static_cast<double>(n) * k * k)) return false;
    lin(X, n, k, beta.data(), off, eta.data());
    for (int i = 0; i < n; ++i) {
      const double mu = std::exp(eta[i]);
      w[i] = wf(mu);
      z[i] = eta[i] - off[i] + (y[i] - mu) / mu;
    }
    xwx(X, n, k, w.data(), z.data(), a, r);
    if (!solve(a, r, k, nb)) return false;
    double dm = 0, bm = 0;
    for (int j = 0; j < k; ++j) {
      dm = std::max(dm, std::fabs(nb[j] - beta[j]));
      bm = std::max(bm, std::fabs(beta[j]));
    }
    beta = nb;
    if (dm < 1e-12 * (1 + bm)) break;
  }
  return true;
}

// Moran's I of e on the dense weights A (n x n) with total s0
double moran(const double *e, int n, const vd &A, double s0) {
  double m = 0;
  for (int i = 0; i < n; ++i) m += e[i];
  m /= n;
  double num = 0, den = 0;
  for (int i = 0; i < n; ++i) {
    const double di = e[i] - m;
    den += di * di;
    for (int j = 0; j < n; ++j) num += A[at(i, j, n)] * di * (e[j] - m);
  }
  return n / s0 * num / den;
}

// Newton-Raphson on the analytic gradient g of the objective f (to be
// maximised), the Hessian by central differences of g, steps halved until
// f does not fall; stops when a step is below tol in every coordinate or
// after 200 iterations. H: the Hessian of the last iteration (k x k).
template <class Fo, class Gr>
int newton(Fo f, Gr g, int k, double tol, vd &x, double *fx, vd &H, int *iters) {
  x.assign(static_cast<std::size_t>(k), 0.0);
  H.assign(static_cast<std::size_t>(k) * k, 0.0);
  double fv = f(x.data());
  vd gv(static_cast<std::size_t>(k)), gp(gv), gm(gv), xp(x), xm(x), st, xn(x);
  int it = 0;
  for (it = 1; it <= 200; ++it) {
    if (kutil::stop(static_cast<double>(k) * k * 100)) return CLAB_EINTR;
    g(x.data(), gv.data());
    for (int j = 0; j < k; ++j) {
      const double h = 1e-5 * std::max(1.0, std::fabs(x[j]));
      xp = x;
      xm = x;
      xp[j] += h;
      xm[j] -= h;
      g(xp.data(), gp.data());
      g(xm.data(), gm.data());
      for (int l = 0; l < k; ++l) H[at(j, l, k)] = (gp[l] - gm[l]) / (2 * h);
    }
    for (int j = 0; j < k; ++j)
      for (int l = j + 1; l < k; ++l) H[at(j, l, k)] = H[at(l, j, k)] = (H[at(j, l, k)] + H[at(l, j, k)]) / 2;
    vd nH(H);
    for (double &e : nH) e = -e;
    if (!solve(nH, gv, k, st)) return CLAB_EDOM;
    double s = 1, fn = 0;
    for (;;) {
      for (int j = 0; j < k; ++j) xn[j] = x[j] + s * st[j];
      fn = f(xn.data());
      if (!std::isnan(fn) && fn >= fv - 1e-12 * std::fabs(fv)) break;
      s /= 2;
      if (s < 1e-12) {
        xn = x;
        fn = fv;
        break;
      }
    }
    double mx = 0;
    for (int j = 0; j < k; ++j) mx = std::max(mx, std::fabs(s * st[j]));
    x = xn;
    fv = fn;
    if (mx < tol) break;
  }
  *fx = fv;
  *iters = std::min(it, 200);
  return CLAB_OK;
}

// the standard errors from the Hessian H (k x k): sqrt(diag(solve(-H))),
// NaN where not positive
bool ses(const vd &H, int k, double *se) {
  vd nH(H), I(static_cast<std::size_t>(k) * k, 0.0);
  for (double &e : nH) e = -e;
  for (int j = 0; j < k; ++j) I[at(j, j, k)] = 1;
  std::vector<int> piv(static_cast<std::size_t>(k));
  if (!gla::lusolve(nH.data(), k, I.data(), k, piv.data())) return false;
  for (int j = 0; j < k; ++j) {
    const double d = I[at(j, j, k)];
    se[j] = d > 0 ? std::sqrt(d) : kutil::nan();
  }
  return true;
}

// The survival data of the frailty and cure models.
struct Sv {
  const double *t, *d, *x;
  const int *reg;
  int n, p, R;
  double tau;
  vd lt, Q;
};

bool svok(Sv &s, const int *rp, const int *ci, const double *w) {
  if (s.n < 1 || s.p < 0 || s.R < 1 || !kutil::nonneg(s.tau)) return false;
  for (int i = 0; i < s.n; ++i) {
    if (!kutil::pos(s.t[i]) || !(s.d[i] == 0 || s.d[i] == 1) || s.reg[i] < 0 || s.reg[i] >= s.R) return false;
  }
  for (std::size_t q = 0; q < static_cast<std::size_t>(s.n) * s.p; ++q) {
    if (!std::isfinite(s.x[q])) return false;
  }
  vd A;
  if (!dense(rp, ci, w, s.R, A)) return false;
  s.Q = lap(A, s.R);
  s.lt.resize(static_cast<std::size_t>(s.n));
  for (int i = 0; i < s.n; ++i) s.lt[i] = std::log(s.t[i]);
  return true;
}

// 0.5 tau u'Qu and tau Q u
double pen(const Sv &s, const double *u, double *qu) {
  double q = 0;
  for (int r = 0; r < s.R; ++r) {
    double a = 0;
    for (int l = 0; l < s.R; ++l) a += s.Q[at(r, l, s.R)] * u[l];
    if (qu) qu[r] = s.tau * a;
    q += u[r] * a;
  }
  return 0.5 * s.tau * q;
}

double xb(const Sv &s, const double *b, int i) {
  double e = 0;
  for (int j = 0; j < s.p; ++j) e += s.x[at(i, j, s.n)] * b[j];
  return e;
}

}  // namespace

int clab_bayout(const double *x, int n, int freq, int b, int w, int acty, double alpha, const int *at, int m,
                double *ub, int *alarm) CLAB_TRY {
  if (n < 1 || freq < 0 || b < 0 || w < 0 || m < 0 || !(alpha > 0 && alpha < 1)) return CLAB_EARG;
  for (int q = 0; q < m; ++q) {
    const int t = at[q];
    if (t < 1 || t > n || static_cast<long long>(t) - static_cast<long long>(b) * freq - w < 1) return CLAB_EDOM;
    double s = 0, c = 0;
    auto add = [&](int k) {  // k from 1
      if (!std::isnan(x[k - 1])) {
        s += x[k - 1];
        c += 1;
      }
    };
    if (acty) {
      for (int k = t - w; k <= t - 1; ++k) add(k);
    }
    for (int i = 1; i <= b; ++i)
      for (int k = t - i * freq - w; k <= t - i * freq + w; ++k) add(k);
    if (c == 0) return CLAB_EDOM;
    if (kutil::stop(1e4)) return CLAB_EINTR;
    ub[q] = kutil::qnbin(1 - alpha, s + 0.5, c / (c + 1));
    if (std::isnan(ub[q])) return CLAB_EARG;  // counts too large
    alarm[q] = std::isnan(x[t - 1]) ? -1 : (x[t - 1] > ub[q] ? 1 : 0);
  }
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_bufexp(const double *rx, const double *ry, int n, const double *sx, const double *sy, int m, double r,
                const double *w, int ll, double *cnt, double *ws, double *near) CLAB_TRY {
  if (n < 1 || m < 1 || !kutil::nonneg(r)) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    double c = 0, s = 0, mn = HUGE_VAL;
    for (int j = 0; j < m; ++j) {
      double d;
      if (ll) {
        if (clab_hvdist(rx + i, ry + i, sx + j, sy + j, 1, 6371008.8, &d) != CLAB_OK) return CLAB_EDOM;
      } else {
        const double a = rx[i] - sx[j], b = ry[i] - sy[j];
        d = std::sqrt(a * a + b * b);
      }
      if (!std::isfinite(d)) return CLAB_EDOM;
      if (d <= r) {
        c += 1;
        s += w ? w[j] : 1.0;
      }
      mn = std::min(mn, d);
    }
    cnt[i] = c;
    ws[i] = s;
    near[i] = mn;
  }
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_kerexp(const double *rx, const double *ry, int n, const double *sx, const double *sy, int m, double h,
                const double *w, int kern, double *out) CLAB_TRY {
  if (n < 1 || m < 1 || kern < 0 || kern > 1) return CLAB_EARG;
  if (!kutil::pos(h)) return CLAB_EDOM;
  for (int i = 0; i < n; ++i) {
    double s = 0;
    for (int j = 0; j < m; ++j) {
      const double a = sx[j] - rx[i], b = sy[j] - ry[i], u = std::sqrt(a * a + b * b) / h;
      if (!std::isfinite(u)) return CLAB_EDOM;
      const double kv = kern == 0 ? std::exp(-u * u / 2) / (2 * pi)
                                  : (u < 1 ? 3 / pi * (1 - u * u) * (1 - u * u) : 0.0);
      s += (w ? w[j] : 1.0) * kv;
    }
    out[i] = s / (h * h);
  }
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_lerprc(const int *rp, const int *ci, const double *w, int n, double rho, double tau, double *q) CLAB_TRY {
  vd A;
  if (n < 1 || !dense(rp, ci, w, n, A)) return CLAB_EARG;
  if (!(rho >= 0 && rho <= 1) || !kutil::pos(tau)) return CLAB_EDOM;
  const vd L = lap(A, n);
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i) q[at(i, j, n)] = tau * (rho * L[at(i, j, n)] + (1 - rho) * (i == j ? 1.0 : 0.0));
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_bym2(const int *rp, const int *ci, const double *w, int n, double tau, double phi, double tol,
              double *scale, double *g, double *cov) CLAB_TRY {
  vd A;
  if (n < 2 || !dense(rp, ci, w, n, A)) return CLAB_EARG;
  if (!kutil::pos(tau) || !(phi >= 0 && phi <= 1) || !kutil::pos(tol)) return CLAB_EDOM;
  vd Q = lap(A, n), val(static_cast<std::size_t>(n)), vec(Q.size());
  if (!gla::jacobi(Q.data(), n, val.data(), vec.data())) return CLAB_EDOM;
  const double top = val[0];
  std::fill(g, g + static_cast<std::size_t>(n) * n, 0.0);
  for (int k = 0; k < n; ++k) {
    if (!(val[k] > tol * top)) continue;
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < n; ++i) g[at(i, j, n)] += vec[at(i, k, n)] * vec[at(j, k, n)] / val[k];
  }
  double sl = 0;
  for (int i = 0; i < n; ++i) {
    if (!(g[at(i, i, n)] > 0)) return CLAB_EDOM;
    sl += std::log(g[at(i, i, n)]);
  }
  const double s = std::exp(sl / n);
  *scale = s;
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i) cov[at(i, j, n)] = ((1 - phi) * (i == j ? 1.0 : 0.0) + phi * g[at(i, j, n)] / s) / tau;
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_ecolrg(const double *y, const double *e, int n, const double *x, int p, int fam, const double *z, int q,
                double tol, int maxit, double *coef, double *zcoef, double *theta, double *fit, double *ll) CLAB_TRY {
  if (n < 1 || p < 0 || q < 0 || fam < 0 || fam > 2 || maxit < 1) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (!kutil::nonneg(y[i]) || !kutil::pos(e[i])) return CLAB_EDOM;
  }
  const int k = p + 1;
  const vd X = design(x, n, p);
  vd off(static_cast<std::size_t>(n)), beta, eta(static_cast<std::size_t>(n)), mu(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) off[i] = std::log(e[i]);
  auto fitted = [&](const vd &b) {
    lin(X, n, k, b.data(), off.data(), eta.data());
    for (int i = 0; i < n; ++i) mu[i] = std::exp(eta[i]);
  };
  if (!irls(y, X, n, k, off.data(), [](double m) { return m; }, beta)) return kutil::stop(0) ? CLAB_EINTR : CLAB_EDOM;
  fitted(beta);
  double L = 0;
  *theta = kutil::nan();
  if (fam == 0) {
    for (int i = 0; i < n; ++i) L += y[i] * std::log(mu[i]) - mu[i] - std::lgamma(y[i] + 1);
  } else if (fam == 1) {
    double ss = 0;
    for (int i = 0; i < n; ++i) ss += (y[i] / mu[i] - 1) * (y[i] / mu[i] - 1);
    double th = n / ss;
    vd nb;
    for (int it = 0; it < maxit; ++it) {
      if (kutil::stop(static_cast<double>(n) * 100)) return CLAB_EINTR;
      for (int j = 0; j < 100; ++j) {
        // Counts with no overdispersion send theta off without bound. Past
        // a million the variance mu + mu^2 / theta is the Poisson's to six
        // digits, and the score and its derivative, which fall as 1 / th^2
        // and 1 / th^3, are soon no more than rounding of terms of order
        // one: a step taken from them is noise (and went to a wrong fit
        // with the extended precision of x87). The search ends at the cap.
        if (th >= 1e6) {
          th = 1e10;
          break;
        }
        double sc = 0, inf = 0;
        const double dth = dig(th), tth = trig(th);
        for (int i = 0; i < n; ++i) {
          sc += dig(y[i] + th) - dth + std::log(th) + 1 - std::log(th + mu[i]) - (y[i] + th) / (mu[i] + th);
          inf += -trig(y[i] + th) + tth - 1 / th + 2 / (mu[i] + th) - (y[i] + th) / ((mu[i] + th) * (mu[i] + th));
        }
        if (!(inf > 0)) {
          th = std::min(2 * th, 1e10);
          if (th >= 1e10) break;
          continue;
        }
        const double d = sc / inf;
        th = std::max(th + d, 1e-8);
        if (std::fabs(d) < tol * (1 + th)) break;
      }
      const double tt = th;
      if (!irls(y, X, n, k, off.data(), [tt](double m) { return m / (1 + m / tt); }, nb)) {
        return kutil::stop(0) ? CLAB_EINTR : CLAB_EDOM;
      }
      double dm = 0;
      for (int j = 0; j < k; ++j) dm = std::max(dm, std::fabs(nb[j] - beta[j]));
      beta = nb;
      fitted(beta);
      if (dm < 1e-10) break;
    }
    for (int i = 0; i < n; ++i) {
      L += kutil::nbll(y[i], mu[i], th);
    }
    *theta = th;
  } else {
    const int kz = q + 1;
    const vd Z = design(z, n, q);
    vd gam(static_cast<std::size_t>(kz), 0.0), pr(static_cast<std::size_t>(n)), zz(static_cast<std::size_t>(n)),
        wt(static_cast<std::size_t>(n)), a, r, st, ez(static_cast<std::size_t>(n)), res(static_cast<std::size_t>(n));
    double old = -HUGE_VAL;
    auto probs = [&](const vd &g) {
      lin(Z, n, kz, g.data(), nullptr, ez.data());
      for (int i = 0; i < n; ++i) pr[i] = 1 / (1 + std::exp(-ez[i]));
    };
    const std::int64_t lim = static_cast<std::int64_t>(maxit) * 20;
    for (std::int64_t it = 0; it < lim; ++it) {
      if (kutil::stop(static_cast<double>(n) * 100)) return CLAB_EINTR;
      fitted(beta);
      probs(gam);
      for (int i = 0; i < n; ++i) {
        zz[i] = y[i] == 0 ? pr[i] / (pr[i] + (1 - pr[i]) * std::exp(-mu[i])) : 0.0;
        wt[i] = 1 - zz[i];
      }
      vd bb(beta);
      for (int j = 0; j < 50; ++j) {
        fitted(bb);
        // X'(w mu)X st = X'(w (y - mu))
        for (int i = 0; i < n; ++i) {
          ez[i] = wt[i] * mu[i];
          res[i] = wt[i] * (y[i] - mu[i]);
        }
        xwx(X, n, k, ez.data(), y, a, r);
        for (int l = 0; l < k; ++l) {
          double s = 0;
          for (int i = 0; i < n; ++i) s += X[at(i, l, n)] * res[i];
          r[l] = s;
        }
        if (!solve(a, r, k, st)) return CLAB_EDOM;
        double mx = 0;
        for (int l = 0; l < k; ++l) {
          bb[l] += st[l];
          mx = std::max(mx, std::fabs(st[l]));
        }
        if (mx < 1e-13) break;
      }
      beta = bb;
      vd gm(gam);
      for (int j = 0; j < 50; ++j) {
        probs(gm);
        vd pw(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) pw[i] = pr[i] * (1 - pr[i]);
        xwx(Z, n, kz, pw.data(), zz.data(), a, r);
        for (int l = 0; l < kz; ++l) {
          double s = 0;
          for (int i = 0; i < n; ++i) s += Z[at(i, l, n)] * (zz[i] - pr[i]);
          r[l] = s;
        }
        if (!solve(a, r, kz, st)) return CLAB_EDOM;
        double mx = 0;
        for (int l = 0; l < kz; ++l) {
          gm[l] += st[l];
          mx = std::max(mx, std::fabs(st[l]));
        }
        if (mx < 1e-13) break;
      }
      gam = gm;
      fitted(beta);
      probs(gam);
      L = 0;
      for (int i = 0; i < n; ++i) {
        L += y[i] == 0 ? std::log(pr[i] + (1 - pr[i]) * std::exp(-mu[i]))
                       : std::log(1 - pr[i]) + y[i] * std::log(mu[i]) - mu[i] - std::lgamma(y[i] + 1);
      }
      if (std::fabs(L - old) < tol * (1 + std::fabs(L))) break;
      old = L;
    }
    std::copy(gam.begin(), gam.end(), zcoef);
  }
  std::copy(beta.begin(), beta.end(), coef);
  std::copy(mu.begin(), mu.end(), fit);
  *ll = L;
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_spbias(const int *d, int nr, int nc, const double *b, int nb, double ub, double *out) CLAB_TRY {
  if (nr < 1 || nc < 1 || nb < 1) return CLAB_EARG;
  for (std::size_t q = 0; q < static_cast<std::size_t>(nr) * nc; ++q) {
    if (d[q] < -1) return CLAB_EDOM;
    out[q] = d[q] == -1 ? ub : b[std::min(d[q], nb - 1)];
  }
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_spsir(const double *i0, int n, const double *pop, const int *rp, const int *ci, const double *w, double beta,
               double gamma, double kappa, double dt, int steps, double *S, double *I, double *R) CLAB_TRY {
  if (n < 1 || steps < 0) return CLAB_EARG;
  if (!kutil::nonneg(beta) || !kutil::nonneg(gamma) || !(kappa >= 0 && kappa <= 1) || !kutil::pos(dt)) {
    return CLAB_EDOM;
  }
  for (int i = 0; i < n; ++i) {
    if (!kutil::nonneg(i0[i]) || !kutil::pos(pop[i]) || i0[i] > pop[i]) return CLAB_EDOM;
  }
  vd A(static_cast<std::size_t>(n) * n, 0.0);
  if (rp != nullptr && !dense(rp, ci, w, n, A)) return CLAB_EARG;
  // the coupling (1 - kappa) I + kappa W~, W~ the row-standardised weights
  vd C(A.size());
  for (int i = 0; i < n; ++i) {
    double s = 0;
    for (int j = 0; j < n; ++j) s += A[at(i, j, n)];
    const double d = s > 0 ? s : 1;
    for (int j = 0; j < n; ++j) C[at(i, j, n)] = (i == j ? 1 - kappa : 0.0) + kappa * (A[at(i, j, n)] / d);
  }
  const std::size_t T = static_cast<std::size_t>(steps) + 1;
  vd s(i0, i0 + n), in(i0, i0 + n), r(static_cast<std::size_t>(n), 0.0);
  for (int i = 0; i < n; ++i) s[i] = pop[i] - i0[i];
  auto put = [&](int k) {
    for (int i = 0; i < n; ++i) {
      S[k + T * i] = s[i];
      I[k + T * i] = in[i];
      R[k + T * i] = r[i];
    }
  };
  put(0);
  vd fs[4], fi[4], fr[4], ts(static_cast<std::size_t>(n)), ti(ts), ip(ts), fo(ts);
  auto rhs = [&](const vd &a, const vd &b, int q) {
    for (int i = 0; i < n; ++i) ip[i] = b[i] / pop[i];
    for (int i = 0; i < n; ++i) {
      double f = 0;
      for (int j = 0; j < n; ++j) f += C[at(i, j, n)] * ip[j];
      fo[i] = f;
    }
    fs[q].resize(static_cast<std::size_t>(n));
    fi[q].resize(static_cast<std::size_t>(n));
    fr[q].resize(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
      fs[q][i] = -beta * a[i] * fo[i];
      fi[q][i] = beta * a[i] * fo[i] - gamma * b[i];
      fr[q][i] = gamma * b[i];
    }
  };
  for (int k = 1; k <= steps; ++k) {
    if (kutil::stop(static_cast<double>(n) * n * 4)) return CLAB_EINTR;
    rhs(s, in, 0);
    for (int i = 0; i < n; ++i) {
      ts[i] = s[i] + 0.5 * dt * fs[0][i];
      ti[i] = in[i] + 0.5 * dt * fi[0][i];
    }
    rhs(ts, ti, 1);
    for (int i = 0; i < n; ++i) {
      ts[i] = s[i] + 0.5 * dt * fs[1][i];
      ti[i] = in[i] + 0.5 * dt * fi[1][i];
    }
    rhs(ts, ti, 2);
    for (int i = 0; i < n; ++i) {
      ts[i] = s[i] + dt * fs[2][i];
      ti[i] = in[i] + dt * fi[2][i];
    }
    rhs(ts, ti, 3);
    for (int i = 0; i < n; ++i) {
      s[i] = s[i] + dt / 6 * (fs[0][i] + 2 * fs[1][i] + 2 * fs[2][i] + fs[3][i]);
      in[i] = in[i] + dt / 6 * (fi[0][i] + 2 * fi[1][i] + 2 * fi[2][i] + fi[3][i]);
      r[i] = r[i] + dt / 6 * (fr[0][i] + 2 * fr[1][i] + 2 * fr[2][i] + fr[3][i]);
    }
    put(k);
  }
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_sfrail(const double *t, const double *d, int n, const double *x, int p, const int *reg, const int *rp,
                const int *ci, const double *w, int nr, double tau, double tol, double *th, double *se, double *f,
                int *iters) CLAB_TRY {
  Sv s{t, d, x, reg, n, p, nr, tau, {}, {}};
  if (!svok(s, rp, ci, w)) return CLAB_EDOM;
  const int k = 1 + p + nr;
  vd eta(static_cast<std::size_t>(n)), hc(eta), qu(static_cast<std::size_t>(nr));
  auto parts = [&](const double *v) {
    const double rho = std::exp(v[0]);
    for (int i = 0; i < n; ++i) {
      eta[i] = xb(s, v + 1, i) + v[1 + p + reg[i]];
      hc[i] = std::exp(rho * s.lt[i] + eta[i]);
    }
    return rho;
  };
  auto obj = [&](const double *v) {
    parts(v);
    double a = 0;
    for (int i = 0; i < n; ++i) a += d[i] * (v[0] + (std::exp(v[0]) - 1) * s.lt[i] + eta[i]) - hc[i];
    return a - pen(s, v + 1 + p, nullptr);
  };
  auto grad = [&](const double *v, double *g) {
    const double rho = parts(v);
    double g0 = 0;
    for (int i = 0; i < n; ++i) g0 += d[i] * (1 + rho * s.lt[i]) - hc[i] * rho * s.lt[i];
    g[0] = g0;
    for (int j = 0; j < p; ++j) {
      double a = 0;
      for (int i = 0; i < n; ++i) a += x[at(i, j, n)] * (d[i] - hc[i]);
      g[1 + j] = a;
    }
    pen(s, v + 1 + p, qu.data());
    for (int r = 0; r < nr; ++r) g[1 + p + r] = -qu[r];
    for (int i = 0; i < n; ++i) g[1 + p + reg[i]] += d[i] - hc[i];
  };
  vd xo, H;
  const int rc = newton(obj, grad, k, tol, xo, f, H, iters);
  if (rc != CLAB_OK) return rc;
  std::copy(xo.begin(), xo.end(), th);
  th[0] = std::exp(xo[0]);
  return ses(H, k, se) ? CLAB_OK : CLAB_EDOM;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_sfcure(const double *t, const double *d, int n, const double *x, int p, const int *reg, const int *rp,
                const int *ci, const double *w, int nr, double tau, double tol, double *th, double *se, double *f,
                double *cure, int *iters) CLAB_TRY {
  Sv s{t, d, x, reg, n, p, nr, tau, {}, {}};
  if (!svok(s, rp, ci, w)) return CLAB_EDOM;
  const int k = 2 + p + nr;
  vd eta(static_cast<std::size_t>(n)), G(eta), qu(static_cast<std::size_t>(nr));
  auto parts = [&](const double *v) {
    const double rho = std::exp(v[0]);
    for (int i = 0; i < n; ++i) {
      eta[i] = xb(s, v + 2, i) + v[2 + p + reg[i]];
      G[i] = std::exp(v[1] + rho * s.lt[i]);
    }
    return rho;
  };
  auto obj = [&](const double *v) {
    const double rho = parts(v);
    double a = 0;
    for (int i = 0; i < n; ++i) {
      a += d[i] * (eta[i] + v[0] + v[1] + (rho - 1) * s.lt[i] - G[i]) - std::exp(eta[i]) * (1 - std::exp(-G[i]));
    }
    return a - pen(s, v + 2 + p, nullptr);
  };
  auto grad = [&](const double *v, double *g) {
    const double rho = parts(v);
    double g0 = 0, g1 = 0;
    for (int r = 0; r < nr; ++r) g[2 + p + r] = 0;
    for (int j = 0; j < p; ++j) g[2 + j] = 0;
    for (int i = 0; i < n; ++i) {
      const double tv = std::exp(eta[i]), dF = std::exp(-G[i]) * G[i], e = d[i] - tv * (1 - std::exp(-G[i]));
      g0 += d[i] * (1 + rho * s.lt[i] - G[i] * rho * s.lt[i]) - tv * dF * rho * s.lt[i];
      g1 += d[i] * (1 - G[i]) - tv * dF;
      for (int j = 0; j < p; ++j) g[2 + j] += x[at(i, j, n)] * e;
      g[2 + p + reg[i]] += e;
    }
    g[0] = g0;
    g[1] = g1;
    pen(s, v + 2 + p, qu.data());
    for (int r = 0; r < nr; ++r) g[2 + p + r] -= qu[r];
  };
  vd xo, H;
  const int rc = newton(obj, grad, k, tol, xo, f, H, iters);
  if (rc != CLAB_OK) return rc;
  std::copy(xo.begin(), xo.end(), th);
  th[0] = std::exp(xo[0]);
  for (int i = 0; i < n; ++i) cure[i] = std::exp(-std::exp(xb(s, xo.data() + 2, i) + xo[2 + p + reg[i]]));
  return ses(H, k, se) ? CLAB_OK : CLAB_EDOM;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_evfilt(const double *y, int n, const double *x, int p, const int *rp, const int *ci, const double *w,
                int crit, int cand, double tol, int maxv, int *sel, int *nsel, double *coef, double *st,
                double *res, double *mi) CLAB_TRY {
  if (n < 3 || p < 0 || crit < 0 || crit > 3 || cand < 0 || cand > 2) return CLAB_EARG;
  vd A;
  if (!dense(rp, ci, w, n, A)) return CLAB_EARG;
  double s0 = 0;
  for (double a : A) s0 += a;
  if (!(s0 > 0)) return CLAB_EDOM;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(y[i])) return CLAB_EDOM;
  }
  // Moran eigenvectors: the eigenvectors of the symmetrised weights in the
  // Helmert basis of the complement of the constant
  const int m = n - 1;
  vd Hq(static_cast<std::size_t>(n) * m, 0.0);
  for (int k = 0; k < m; ++k) {
    const double cc = 1 / std::sqrt(static_cast<double>(k + 1) * (k + 2));
    for (int i = 0; i <= k; ++i) Hq[at(i, k, n)] = cc;
    Hq[at(k + 1, k, n)] = -(k + 1) * cc;
  }
  vd SQ(static_cast<std::size_t>(n) * m, 0.0), M(static_cast<std::size_t>(m) * m, 0.0);
  for (int k = 0; k < m; ++k)
    for (int i = 0; i < n; ++i) {
      double a = 0;
      for (int j = 0; j < n; ++j) a += (A[at(i, j, n)] + A[at(j, i, n)]) / 2 * Hq[at(j, k, n)];
      SQ[at(i, k, n)] = a;
    }
  for (int k = 0; k < m; ++k)
    for (int l = 0; l < m; ++l) {
      double a = 0;
      for (int i = 0; i < n; ++i) a += Hq[at(i, k, n)] * SQ[at(i, l, n)];
      M[at(k, l, m)] = a;
    }
  vd val(static_cast<std::size_t>(m)), vec(M.size());
  if (!gla::jacobi(M.data(), m, val.data(), vec.data())) return CLAB_EDOM;
  vd E(static_cast<std::size_t>(n) * m, 0.0);
  for (int k = 0; k < m; ++k) {
    for (int i = 0; i < n; ++i) {
      double a = 0;
      for (int l = 0; l < m; ++l) a += Hq[at(i, l, n)] * vec[at(l, k, m)];
      E[at(i, k, n)] = a;
    }
    int f = 0;
    while (f < n && !(std::fabs(E[at(f, k, n)]) > 1e-12)) ++f;
    if (f < n && E[at(f, k, n)] < 0) {
      for (int i = 0; i < n; ++i) E[at(i, k, n)] = -E[at(i, k, n)];
    }
    mi[k] = n / s0 * val[k];
  }
  const double e0 = -1.0 / (n - 1);
  std::vector<int> pool;
  for (int k = 0; k < m; ++k) {
    if (cand == 2 || (cand == 0 && mi[k] > e0) || (cand == 1 && mi[k] < e0)) pool.push_back(k);
  }
  const int kx = p + 1;
  const vd Xb = design(x, n, p);
  double my = 0;
  for (int i = 0; i < n; ++i) my += y[i];
  my /= n;
  double tss = 0;
  for (int i = 0; i < n; ++i) tss += (y[i] - my) * (y[i] - my);
  struct Fit {
    vd beta, res;
    double aic, r2, I, press;
  };
  auto fit = [&](const std::vector<int> &sl, Fit &o) {
    const int k = kx + static_cast<int>(sl.size());
    vd X2(static_cast<std::size_t>(n) * k);
    std::copy(Xb.begin(), Xb.end(), X2.begin());
    for (std::size_t q = 0; q < sl.size(); ++q) {
      for (int i = 0; i < n; ++i) X2[at(i, kx + static_cast<int>(q), n)] = E[at(i, sl[q], n)];
    }
    vd a, r, ai(static_cast<std::size_t>(k) * k, 0.0);
    xwx(X2, n, k, nullptr, y, a, r);
    for (int j = 0; j < k; ++j) ai[at(j, j, k)] = 1;
    vd a2(a);
    std::vector<int> piv(static_cast<std::size_t>(k));
    if (!gla::lusolve(a2.data(), k, ai.data(), k, piv.data()) || !solve(a, r, k, o.beta)) return false;
    o.res.assign(static_cast<std::size_t>(n), 0.0);
    double rss = 0, press = 0;
    for (int i = 0; i < n; ++i) {
      double fi = 0, h = 0;
      for (int j = 0; j < k; ++j) {
        fi += X2[at(i, j, n)] * o.beta[j];
        double t2 = 0;
        for (int l = 0; l < k; ++l) t2 += X2[at(i, l, n)] * ai[at(l, j, k)];
        h += t2 * X2[at(i, j, n)];
      }
      o.res[i] = y[i] - fi;
      rss += o.res[i] * o.res[i];
      const double u = o.res[i] / (1 - h);
      press += u * u;
    }
    o.aic = n * std::log(rss / n) + 2 * (k + 1);
    o.r2 = 1 - rss / tss;
    o.I = moran(o.res.data(), n, A, s0);
    o.press = press;
    return true;
  };
  Fit base, cur, f;
  std::vector<int> cs;
  if (!fit(cs, base)) return CLAB_EDOM;
  cur = base;
  int limit = maxv < 0 ? static_cast<int>(pool.size()) : std::min(maxv, static_cast<int>(pool.size()));
  limit = std::min(limit, n - kx - 2);
  while (static_cast<int>(cs.size()) < limit) {
    if (kutil::stop(static_cast<double>(n) * n * m)) return CLAB_EINTR;
    if (crit == 0 && std::fabs(cur.I) < tol) break;
    int bk = -1;
    double bv = HUGE_VAL;
    Fit best;
    for (int k : pool) {
      if (std::find(cs.begin(), cs.end(), k) != cs.end()) continue;
      std::vector<int> t2(cs);
      t2.push_back(k);
      if (!fit(t2, f)) return CLAB_EDOM;
      const double key = crit == 0 ? std::fabs(f.I) : crit == 1 ? f.aic : crit == 2 ? f.press : -f.r2;
      if (key < bv) {  // ties: the first, the lowest index
        bv = key;
        bk = k;
        best = f;
      }
    }
    // no candidate with a finite criterion (values so large the sums
    // overflow): nothing to add
    if (bk < 0) break;
    if (crit == 1 && best.aic >= cur.aic) break;
    if (crit == 2 && best.press >= cur.press) break;
    if (crit == 3 && best.r2 - cur.r2 < tol) break;
    cs.push_back(bk);
    cur = best;
  }
  *nsel = static_cast<int>(cs.size());
  for (std::size_t q = 0; q < cs.size(); ++q) sel[q] = cs[q] + 1;
  std::copy(cur.beta.begin(), cur.beta.end(), coef);
  std::copy(cur.res.begin(), cur.res.end(), res);
  st[0] = cur.I;
  st[1] = cur.aic;
  st[2] = cur.r2;
  st[3] = cur.r2 - base.r2;
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)
