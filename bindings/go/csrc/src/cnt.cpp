// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, count regression: Poisson and negative binomial log-linear
// models by iteratively reweighted least squares (the steps of R's glm.fit
// and of MASS::glm.nb), with model-based and cluster-robust covariances
// (clab_cntreg), and the Poisson model by generalised estimating equations
// with an exchangeable working correlation (clab_geeex). No allocation:
// the caller passes the work arrays.
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <limits>

#include "crimlab.h"

namespace {

using std::size_t;

constexpr double qnan = std::numeric_limits<double>::quiet_NaN();
// R's glm.fit: a column is aliased when its norm left after the earlier
// columns is below this share of its own (min(1e-7, epsilon / 1000)).
constexpr double qrtol = 1e-11;

// a * b as a size, for array lengths
size_t sz(int a, int b) { return static_cast<size_t>(a) * static_cast<size_t>(b); }

size_t at(int i, int j, int n) {
  return static_cast<size_t>(i) + sz(j, n);
}

// digamma and trigamma for x > 0: the recurrence up to x >= 10, then the
// asymptotic series to the x^-12 term (error below 1e-14 there).
double digam(double x) {
  double s = 0.0;
  while (x < 10.0) {
    s -= 1.0 / x;
    x += 1.0;
  }
  const double r = 1.0 / (x * x);
  return s + std::log(x) - 0.5 / x -
         r * (1.0 / 12 - r * (1.0 / 120 - r * (1.0 / 252 - r * (1.0 / 240 - r * (1.0 / 132)))));
}

double trigam(double x) {
  double s = 0.0;
  while (x < 10.0) {
    s += 1.0 / (x * x);
    x += 1.0;
  }
  const double r = 1.0 / (x * x);
  return s + 1.0 / x + r / 2 +
         r / x * (1.0 / 6 - r * (1.0 / 30 - r * (1.0 / 42 - r * (1.0 / 30 - r * (5.0 / 66)))));
}

// The family: Poisson when th is 0, otherwise negative binomial with size
// th (variance mu + mu^2 / th). Log link throughout.
struct fam {
  double th;
  double var(double mu) const { return th > 0 ? mu + mu * mu / th : mu; }
  // the unit deviance of glm's dev.resids
  double dev(double y, double mu) const {
    if (th > 0) {
      return 2 * (y * std::log(std::fmax(1.0, y) / mu) -
                  (y + th) * std::log((y + th) / (mu + th)));
    }
    return 2 * ((y > 0 ? y * std::log(y / mu) : 0.0) - (y - mu));
  }
};

// mu from eta, as glm's poisson()$linkinv and negative.binomial()$linkinv:
// pmax(exp(eta), .Machine$double.eps)
double inv(double eta) { return std::fmax(std::exp(eta), DBL_EPSILON); }

// The arrays of one fit. a: n x p (the weighted design, then its QR); b,
// mu, eta, sw, z, mold: n; beta, bold: p; kept: p ints (the columns kept, in
// order, then -1); rk: their number.
struct fit {
  const double *y, *x, *off;
  int n, p;
  double *a, *b, *mu, *eta, *sw, *z, *mold, *beta, *bold;
  int *kept;
  int rk;
};

// One weighted least squares solve: Householder QR of a = sw * x with the
// limited pivoting of R's dqrdc2 (a column whose norm left is below qrtol
// of its own norm is aliased and skipped), b = sw * z rotated alongside.
// beta gets the solution on the kept columns, 0 on the aliased ones. The
// kept columns' R stays in the first rk columns of a.
void wls(fit &f) {
  const int n = f.n, p = f.p;
  for (int j = 0; j < p; ++j) {
    for (int i = 0; i < n; ++i) f.a[at(i, j, n)] = f.sw[i] * f.x[at(i, j, n)];
  }
  for (int i = 0; i < n; ++i) f.b[i] = f.sw[i] * f.z[i];
  int k = 0;
  for (int j = 0; j < p; ++j) {
    double *c = f.a + at(0, j, n);
    double o = 0.0, r = 0.0;
    for (int i = 0; i < n; ++i) o += c[i] * c[i];
    for (int i = k; i < n; ++i) r += c[i] * c[i];
    o = std::sqrt(o);
    r = std::sqrt(r);
    f.beta[j] = 0.0;
    if (k >= n || !(o > 0) || r < qrtol * o) continue;
    // compact: column j into slot k
    double *d = f.a + at(0, k, n);
    if (d != c) {
      for (int i = 0; i < n; ++i) d[i] = c[i];
    }
    f.kept[k] = j;
    // the reflection that zeroes d[k+1..n-1]
    const double al = d[k] > 0 ? -r : r;
    const double v0 = d[k] - al;
    // v = (v0, d[k+1..]); H = I - 2 v v' / (v'v); v'v = 2 r (r + |d[k]|)
    const double vv = v0 * v0 + (r * r - d[k] * d[k]);
    d[k] = al;
    if (vv > 0) {
      auto refl = [&](double *e) {
        double s = v0 * e[k];
        for (int i = k + 1; i < n; ++i) s += d[i] * e[i];
        s = 2 * s / vv;
        e[k] -= s * v0;
        for (int i = k + 1; i < n; ++i) e[i] -= s * d[i];
      };
      for (int jj = j + 1; jj < p; ++jj) refl(f.a + at(0, jj, n));
      refl(f.b);
    }
    for (int i = k + 1; i < n; ++i) d[i] = 0.0;
    ++k;
  }
  f.rk = k;
  for (int j = k; j < p; ++j) f.kept[j] = -1;
  // back substitution on R c = (Q'b)[0..k)
  for (int q = k - 1; q >= 0; --q) {
    double s = f.b[q];
    for (int r2 = q + 1; r2 < k; ++r2) s -= f.a[at(q, r2, n)] * f.b[r2];
    f.b[q] = s / f.a[at(q, q, n)];
  }
  for (int q = 0; q < k; ++q) f.beta[f.kept[q]] = f.b[q];
}

// eta = x beta + off and mu, with the deviance under family fm
double update(fit &f, const fam &fm) {
  double dv = 0.0;
  for (int i = 0; i < f.n; ++i) {
    double e = f.off ? f.off[i] : 0.0;
    for (int j = 0; j < f.p; ++j) e += f.x[at(i, j, f.n)] * f.beta[j];
    f.eta[i] = e;
    f.mu[i] = inv(e);
    dv += fm.dev(f.y[i], f.mu[i]);
  }
  return dv;
}

// glm.fit from the current eta and mu: IRLS until the relative change of
// the deviance is below eps, at most maxit steps, halving a step whose
// deviance is not finite. sw keeps the root weights of the last solve and
// a its R, as glm.fit's weights and qr. The number of steps, or -1
// without convergence; dv the deviance.
int irls(fit &f, const fam &fm, int maxit, double eps, double *dv) {
  double dold = 0.0;
  for (int i = 0; i < f.n; ++i) dold += fm.dev(f.y[i], f.mu[i]);
  for (int it = 1; it <= maxit; ++it) {
    for (int i = 0; i < f.n; ++i) {
      const double me = std::fmax(std::exp(f.eta[i]), DBL_EPSILON);
      const double o = f.off ? f.off[i] : 0.0;
      f.z[i] = (f.eta[i] - o) + (f.y[i] - f.mu[i]) / me;
      f.sw[i] = std::sqrt(me * me / fm.var(f.mu[i]));
    }
    wls(f);
    double d = update(f, fm);
    for (int h = 0; h < 60 && !std::isfinite(d); ++h) {
      for (int j = 0; j < f.p; ++j) f.beta[j] = (f.beta[j] + f.bold[j]) / 2;
      d = update(f, fm);
    }
    if (!std::isfinite(d)) return -2;
    *dv = d;
    if (std::fabs(d - dold) / (std::fabs(d) + 0.1) < eps) return it;
    dold = d;
    for (int j = 0; j < f.p; ++j) f.bold[j] = f.beta[j];
  }
  return -1;
}

// MASS's theta.ml: Newton steps on the profile score of the size th, from
// the moment start; the information at the last step into *info.
double thml(const double *y, const double *mu, int n, int limit, double *info,
            int *ok) {
  double s = 0.0;
  for (int i = 0; i < n; ++i) s += (y[i] / mu[i] - 1) * (y[i] / mu[i] - 1);
  double t0 = n / s, del = 1.0, inf = qnan;
  const double eps = std::pow(DBL_EPSILON, 0.25);
  int it = 0;
  // a theta that leaves the finite range ends the search: no estimate
  while (++it < limit && std::fabs(del) > eps && std::isfinite(t0)) {
    t0 = std::fabs(t0);
    double sc = 0.0, in = 0.0;
    const double dt = digam(t0), tt = trigam(t0);
    for (int i = 0; i < n; ++i) {
      const double a = t0 + y[i], m = mu[i] + t0;
      sc += digam(a) - dt + std::log(t0) + 1 - std::log(m) - a / m;
      in += -trigam(a) + tt - 1 / t0 + 2 / m - a / (m * m);
    }
    inf = in;
    del = sc / in;
    t0 += del;
  }
  *ok = (it < limit && t0 > 0 && std::isfinite(t0)) ? 1 : 0;
  *info = inf;
  return t0;
}

double llnb(const double *y, const double *mu, int n, double th) {
  double s = 0.0;
  for (int i = 0; i < n; ++i) {
    s += std::lgamma(th + y[i]) - std::lgamma(th) - std::lgamma(y[i] + 1) + th * std::log(th) +
         y[i] * std::log(mu[i] + (y[i] == 0 ? 1.0 : 0.0)) - (th + y[i]) * std::log(th + mu[i]);
  }
  return s;
}

double llpois(const double *y, const double *mu, int n) {
  double s = 0.0;
  for (int i = 0; i < n; ++i) {
    s += (y[i] > 0 ? y[i] * std::log(mu[i]) : 0.0) - mu[i] - std::lgamma(y[i] + 1);
  }
  return s;
}

// R^{-1} R^{-T} of the rk x rk upper triangle in a, spread over p x p cov
// (NaN in the rows and columns of aliased columns); ri: rk x rk scratch.
void rcov(const fit &f, double *ri, double *cov) {
  const int n = f.n, p = f.p, k = f.rk;
  for (int j = 0; j < k; ++j) {
    for (int i = 0; i < k; ++i) ri[at(i, j, k)] = 0.0;
    // column j of R^{-1}: solve R c = e_j
    for (int i = j; i >= 0; --i) {
      double s = i == j ? 1.0 : 0.0;
      for (int q = i + 1; q <= j; ++q) s -= f.a[at(i, q, n)] * ri[at(q, j, k)];
      ri[at(i, j, k)] = s / f.a[at(i, i, n)];
    }
  }
  for (size_t q = 0; q < sz(p, p); ++q) cov[q] = qnan;
  for (int u = 0; u < k; ++u) {
    for (int v = 0; v < k; ++v) {
      double s = 0.0;
      for (int q = (u > v ? u : v); q < k; ++q) s += ri[at(u, q, k)] * ri[at(v, q, k)];
      cov[at(f.kept[u], f.kept[v], p)] = s;
    }
  }
}

// bread U' U bread over clusters: u is ncl x p (the score of each cluster
// summed), cov the bread (p x p, NaN where aliased); vc gets the sandwich.
void sandw(const fit &f, const double *cov, const double *u, int ncl, double *t,
           double *vc) {
  const int p = f.p, k = f.rk;
  // t = bread_kept U' : k x ncl
  for (int g = 0; g < ncl; ++g) {
    for (int a2 = 0; a2 < k; ++a2) {
      double s = 0.0;
      for (int b2 = 0; b2 < k; ++b2) s += cov[at(f.kept[a2], f.kept[b2], p)] * u[at(g, f.kept[b2], ncl)];
      t[at(a2, g, k)] = s;
    }
  }
  for (size_t q = 0; q < sz(p, p); ++q) vc[q] = qnan;
  for (int a2 = 0; a2 < k; ++a2) {
    for (int b2 = 0; b2 < k; ++b2) {
      double s = 0.0;
      for (int g = 0; g < ncl; ++g) s += t[at(a2, g, k)] * t[at(b2, g, k)];
      vc[at(f.kept[a2], f.kept[b2], p)] = s;
    }
  }
}

bool chkin(const double *y, const double *x, int n, int p, const double *off,
           const int *cl, int ncl) {
  if (n < 1 || p < 1) return false;
  for (int i = 0; i < n; ++i) {
    if (!(y[i] >= 0) || !std::isfinite(y[i])) return false;
    if (off && !std::isfinite(off[i])) return false;
    if (cl && (cl[i] < 0 || cl[i] >= ncl)) return false;
  }
  for (size_t q = 0; q < sz(n, p); ++q) {
    if (!std::isfinite(x[q])) return false;
  }
  return true;
}

// The inverse of the k x k symmetric positive definite m, in place, by
// Cholesky; false when m is not positive definite.
bool spdinv(double *m, int k) {
  for (int j = 0; j < k; ++j) {
    double d = m[at(j, j, k)];
    for (int q = 0; q < j; ++q) d -= m[at(j, q, k)] * m[at(j, q, k)];
    if (!(d > 0)) return false;
    d = std::sqrt(d);
    m[at(j, j, k)] = d;
    for (int i = j + 1; i < k; ++i) {
      double s = m[at(i, j, k)];
      for (int q = 0; q < j; ++q) s -= m[at(i, q, k)] * m[at(j, q, k)];
      m[at(i, j, k)] = s / d;
    }
  }
  // L^{-1} in the lower triangle
  for (int j = 0; j < k; ++j) {
    m[at(j, j, k)] = 1 / m[at(j, j, k)];
    for (int i = j + 1; i < k; ++i) {
      double s = 0.0;
      for (int q = j; q < i; ++q) s -= m[at(i, q, k)] * m[at(q, j, k)];
      m[at(i, j, k)] = s / m[at(i, i, k)];
    }
  }
  // L^{-T} L^{-1}, symmetric, both triangles
  for (int j = 0; j < k; ++j) {
    for (int i = j; i < k; ++i) {
      double s = 0.0;
      for (int q = i; q < k; ++q) s += m[at(q, i, k)] * m[at(q, j, k)];
      m[at(j, i, k)] = s;
    }
  }
  for (int j = 0; j < k; ++j) {
    for (int i = j + 1; i < k; ++i) m[at(i, j, k)] = m[at(j, i, k)];
  }
  return true;
}

}  // namespace

extern "C" {

int clab_cntreg(const double *y, const double *x, int n, int p, const double *off,
                const int *cl, int ncl, int nb, int maxit, double eps, double *beta,
                double *vb, double *vc, double *st, int *alias, double *work,
                int *iwork) {
  if (maxit < 1 || !(eps > 0) || (cl && ncl < 1)) return CLAB_EARG;
  if (!chkin(y, x, n, p, off, cl, ncl)) return CLAB_EDOM;
  double ysum = 0.0;
  for (int i = 0; i < n; ++i) ysum += y[i];
  if (!(ysum > 0)) return CLAB_EDOM;
  const size_t np = sz(n, p);
  fit f{y, x, off, n, p, work, work + np, nullptr, nullptr, nullptr, nullptr, nullptr,
        beta, nullptr, iwork, 0};
  f.mu = f.b + n;
  f.eta = f.mu + n;
  f.sw = f.eta + n;
  f.z = f.sw + n;
  f.mold = f.z + n;
  f.bold = f.mold + n;
  double *ri = f.bold + p;                      // p x p
  double *u = ri + sz(p, p);  // ncl x p
  double *t = u + sz(cl ? ncl : 0, p);  // p x ncl
  for (int j = 0; j < p; ++j) f.bold[j] = 0.0;
  // glm.fit's start for the Poisson family: mu = y + 0.1
  for (int i = 0; i < n; ++i) {
    f.mu[i] = y[i] + 0.1;
    f.eta[i] = std::log(f.mu[i]);
  }
  double dv = 0.0;
  fam fm{0.0};
  int it = irls(f, fm, maxit, eps, &dv);
  if (it == -2) return CLAB_EDOM;
  int conv = it > 0 ? 1 : 0, iter = it > 0 ? it : maxit;
  double th = qnan, sth = qnan, ll;
  if (nb) {
    // glm.nb: theta by ML at the Poisson fit, then alternate an IRLS at
    // the current theta (from the current mu) with a new theta (at the mu
    // that IRLS started from) until the log-likelihood and theta settle;
    // the fit reported is the last IRLS, at the theta before the last
    double info;
    int ok;
    th = thml(y, f.mu, n, maxit, &info, &ok);
    if (!(th > 0) || !std::isfinite(th)) return CLAB_EDOM;
    const double d1 = std::sqrt(2 * std::fmax(1.0, static_cast<double>(n - f.rk)));
    double del = 1.0, lm = llnb(y, f.mu, n, th), lm0 = lm + 2 * d1;
    int outer = 0, tok = ok;
    while (++outer <= maxit && std::fabs(lm0 - lm) / d1 + std::fabs(del) > eps) {
      fm.th = th;
      for (int i = 0; i < n; ++i) {
        f.mold[i] = f.mu[i];
        f.eta[i] = std::log(f.mu[i]);
        f.mu[i] = inv(f.eta[i]);
      }
      it = irls(f, fm, maxit, eps, &dv);
      if (it == -2) return CLAB_EDOM;
      iter += it > 0 ? it : maxit;
      const double t0 = th;
      th = thml(y, f.mold, n, maxit, &info, &ok);
      if (!(th > 0) || !std::isfinite(th)) return CLAB_EDOM;
      del = t0 - th;
      tok = ok;
      lm0 = lm;
      lm = llnb(y, f.mu, n, th);
    }
    conv = (outer <= maxit && tok && it > 0) ? 1 : 0;
    sth = std::sqrt(1 / info);
    ll = lm;
  } else {
    ll = llpois(y, f.mu, n);
  }
  rcov(f, ri, vb);
  if (cl) {
    for (size_t q = 0; q < sz(ncl, p); ++q) u[q] = 0.0;
    for (int i = 0; i < n; ++i) {
      // the working residual (final mu) times the working weight (last
      // solve), as sandwich's estfun.glm
      const double me = std::fmax(std::exp(f.eta[i]), DBL_EPSILON);
      const double r = (y[i] - f.mu[i]) / me * f.sw[i] * f.sw[i];
      for (int j = 0; j < p; ++j) u[at(cl[i], j, ncl)] += r * x[at(i, j, n)];
    }
    sandw(f, vb, u, ncl, t, vc);
  }
  for (int j = 0; j < p; ++j) alias[j] = 1;
  for (int q = 0; q < f.rk; ++q) alias[f.kept[q]] = 0;
  for (int j = 0; j < p; ++j) {
    if (alias[j]) beta[j] = qnan;
  }
  st[0] = ll;
  st[1] = -2 * ll + 2 * f.rk + (nb ? 2 : 0);
  st[2] = th;
  st[3] = sth;
  st[4] = iter;
  st[5] = conv;
  st[6] = f.rk;
  st[7] = dv;
  return CLAB_OK;
}

int clab_geeex(const double *y, const double *x, int n, int p, const double *off,
               const int *cl, int ncl, const int *alias, int maxit, double eps,
               double *beta, double *vc, double *st, double *work, int *iwork) {
  if (maxit < 1 || !(eps > 0) || ncl < 1) return CLAB_EARG;
  if (!cl || !chkin(y, x, n, p, off, cl, ncl)) return CLAB_EDOM;
  int *kept = iwork, *beg = iwork + p, *ord = beg + ncl + 1;
  int k = 0;
  for (int j = 0; j < p; ++j) {
    if (!alias[j]) {
      if (!std::isfinite(beta[j])) return CLAB_EDOM;
      kept[k++] = j;
    }
  }
  if (k == 0) return CLAB_EARG;
  // the rows of each cluster together: a counting sort on cl
  for (int g = 0; g <= ncl; ++g) beg[g] = 0;
  for (int i = 0; i < n; ++i) ++beg[cl[i] + 1];
  for (int g = 0; g < ncl; ++g) beg[g + 1] += beg[g];
  // ord[beg[g]..beg[g + 1]) are the rows of cluster g; beg is advanced
  // while filling, then moved back one place
  for (int i = 0; i < n; ++i) ord[beg[cl[i]]++] = i;
  for (int g = ncl; g > 0; --g) beg[g] = beg[g - 1];
  beg[0] = 0;
  double *mu = work, *r = mu + n, *m = r + n, *u = m + sz(k, k),
         *sx = u + k, *uc = sx + k, *t = uc + sz(ncl, k);
  double al = 0.0, phi = 0.0;
  int it = 0, conv = 0;
  // one pass: mu, residuals, alpha at beta; then M and U (uc per cluster)
  auto pass = [&]() -> bool {
    double ss = 0.0, num = 0.0, cnt = 0.0;
    for (int i = 0; i < n; ++i) {
      double e = off ? off[i] : 0.0;
      for (int q = 0; q < k; ++q) e += x[at(i, kept[q], n)] * beta[kept[q]];
      mu[i] = inv(e);
      r[i] = (y[i] - mu[i]) / std::sqrt(mu[i]);
      ss += r[i] * r[i];
    }
    for (int g = 0; g < ncl; ++g) {
      double a = 0.0, b = 0.0;
      for (int q = beg[g]; q < beg[g + 1]; ++q) {
        a += r[ord[q]];
        b += r[ord[q]] * r[ord[q]];
      }
      const double ni = beg[g + 1] - beg[g];
      num += (a * a - b) / 2;
      cnt += ni * (ni - 1) / 2;
    }
    phi = ss / n;
    al = cnt > 0 && phi > 0 ? num / (cnt * phi) : 0.0;
    for (size_t q = 0; q < sz(k, k); ++q) m[q] = 0.0;
    for (int q = 0; q < k; ++q) u[q] = 0.0;
    for (int g = 0; g < ncl; ++g) {
      const int ni = beg[g + 1] - beg[g];
      if (!(al < 1) || !(1 + (ni - 1) * al > 0)) return false;
      const double c = 1 / (1 - al), d = al / (1 + (ni - 1) * al);
      double sr = 0.0;
      for (int a2 = 0; a2 < k; ++a2) sx[a2] = 0.0;
      double *ug = uc + g;
      for (int a2 = 0; a2 < k; ++a2) ug[at(0, a2, ncl)] = 0.0;
      for (int q = beg[g]; q < beg[g + 1]; ++q) {
        const int i = ord[q];
        const double sm = std::sqrt(mu[i]);
        sr += r[i];
        for (int a2 = 0; a2 < k; ++a2) {
          const double xa = sm * x[at(i, kept[a2], n)];
          sx[a2] += xa;
          ug[at(0, a2, ncl)] += xa * r[i];
          for (int b2 = 0; b2 <= a2; ++b2) m[at(a2, b2, k)] += c * xa * sm * x[at(i, kept[b2], n)];
        }
      }
      for (int a2 = 0; a2 < k; ++a2) {
        ug[at(0, a2, ncl)] = c * (ug[at(0, a2, ncl)] - d * sx[a2] * sr);
        u[a2] += ug[at(0, a2, ncl)];
        for (int b2 = 0; b2 <= a2; ++b2) m[at(a2, b2, k)] -= c * d * sx[a2] * sx[b2];
      }
    }
    for (int a2 = 0; a2 < k; ++a2) {
      for (int b2 = a2 + 1; b2 < k; ++b2) m[at(a2, b2, k)] = m[at(b2, a2, k)];
    }
    return spdinv(m, k);
  };
  while (it < maxit) {
    ++it;
    if (!pass()) return CLAB_EDOM;
    double big = 0.0;
    for (int a2 = 0; a2 < k; ++a2) {
      double s = 0.0;
      for (int b2 = 0; b2 < k; ++b2) s += m[at(a2, b2, k)] * u[b2];
      beta[kept[a2]] += s;
      big = std::fmax(big, std::fabs(s));
    }
    if (!std::isfinite(big)) return CLAB_EDOM;
    if (big <= eps) {
      conv = 1;
      break;
    }
  }
  // the sandwich at the final estimate
  if (!pass()) return CLAB_EDOM;
  for (int g = 0; g < ncl; ++g) {
    for (int a2 = 0; a2 < k; ++a2) {
      double s = 0.0;
      for (int b2 = 0; b2 < k; ++b2) s += m[at(a2, b2, k)] * uc[at(g, b2, ncl)];
      t[at(a2, g, k)] = s;
    }
  }
  for (size_t q = 0; q < sz(p, p); ++q) vc[q] = qnan;
  for (int a2 = 0; a2 < k; ++a2) {
    for (int b2 = 0; b2 < k; ++b2) {
      double s = 0.0;
      for (int g = 0; g < ncl; ++g) s += t[at(a2, g, k)] * t[at(b2, g, k)];
      vc[at(kept[a2], kept[b2], p)] = s;
    }
  }
  for (int j = 0; j < p; ++j) {
    if (alias[j]) beta[j] = qnan;
  }
  st[0] = al;
  st[1] = phi;
  st[2] = it;
  st[3] = conv;
  return CLAB_OK;
}

}  // extern "C"
