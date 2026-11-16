// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The generalized synthetic control (crimlab.h, ifesc()): Xu (2017) on
// the interactive fixed effects model of Bai (2009). The control units'
// outcomes, net of covariates, are unit and time effects plus r common
// factors with unit loadings; Bai's estimator alternates the factors
// (principal components of the doubly demeaned residuals) and the
// covariates' coefficients. Each treated unit's loadings and effect come
// from its periods before adoption, by least squares on the factors; its
// counterfactual is the model's value after. r is fixed or chosen by
// leaving out one pre-adoption period at a time (Xu 2017, sec. 3.2).
// Work: O(maxit (T^2 Nc + T^3)) per fit.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "crimlab.h"
#include "kutil.h"
#include "svk.h"

namespace {

using std::size_t;
inline size_t z(int i) { return static_cast<size_t>(i); }

// The panel: y (nt x n, a column per unit), x (nt x n x p), the units
// fitted (with repeats) as controls cu and treated tu, adoption a.
struct Pan {
  const double *y, *x;
  int p, n, nt;
  const int *a;
  double yv(int t, int u) const { return y[z(u) * z(nt) + z(t)]; }
  double xv(int t, int u, int k) const { return x[(z(k) * z(n) + z(u)) * z(nt) + z(t)]; }
};

struct Fit {
  std::vector<double> beta, xi, f;  // f: nt x r
  double mu = 0.0;
  int it = 0;
  bool conv = true;
};

// x(t, u) beta
double xb(const Pan &d, int t, int u, const std::vector<double> &b) {
  double s = 0.0;
  for (int k = 0; k < d.p; ++k) s += d.xv(t, u, k) * b[z(k)];
  return s;
}

// m (nt x nc) less its row and column means plus the grand mean.
void ddm(std::vector<double> &m, int nt, int nc) {
  std::vector<double> rm(z(nt), 0.0), cm(z(nc), 0.0);
  double g = 0.0;
  for (int j = 0; j < nc; ++j) {
    for (int t = 0; t < nt; ++t) {
      const double v = m[z(j) * z(nt) + z(t)];
      rm[z(t)] += v / nc;
      cm[z(j)] += v / nt;
      g += v;
    }
  }
  g /= static_cast<double>(nt) * nc;
  for (int j = 0; j < nc; ++j) {
    for (int t = 0; t < nt; ++t) m[z(j) * z(nt) + z(t)] += g - rm[z(t)] - cm[z(j)];
  }
}

// The common component F L' of the doubly demeaned w (nt x nc): F the
// first r eigenvectors of w w' times sqrt(nt), L = w' F / nt.
void comp(const std::vector<double> &w, int nt, int nc, int r, std::vector<double> &f, std::vector<double> &c) {
  c.assign(w.size(), 0.0);
  f.assign(z(nt) * z(r), 0.0);
  if (r == 0) return;
  svk::Mat a(nt, nt), e;
  for (int s = 0; s < nt; ++s) {
    for (int t = s; t < nt; ++t) {
      double v = 0.0;
      for (int j = 0; j < nc; ++j) v += w[z(j) * z(nt) + z(s)] * w[z(j) * z(nt) + z(t)];
      a(s, t) = v;
      a(t, s) = v;
    }
  }
  std::vector<double> ev;
  svk::eigsym(a, ev, e);
  const double sq = std::sqrt(static_cast<double>(nt));
  for (int k = 0; k < r; ++k) {
    for (int t = 0; t < nt; ++t) f[z(k) * z(nt) + z(t)] = e(t, k) * sq;
  }
  for (int j = 0; j < nc; ++j) {
    for (int k = 0; k < r; ++k) {
      double l = 0.0;
      for (int t = 0; t < nt; ++t) l += w[z(j) * z(nt) + z(t)] * f[z(k) * z(nt) + z(t)];
      l /= nt;
      for (int t = 0; t < nt; ++t) c[z(j) * z(nt) + z(t)] += f[z(k) * z(nt) + z(t)] * l;
    }
  }
}

// beta by least squares of dd(y) - c on dd(x) over the controls.
bool bols(const Pan &d, const std::vector<int> &cu, const std::vector<double> &c, std::vector<double> &b) {
  const int nc = static_cast<int>(cu.size()), nt = d.nt;
  const size_t m = z(nc) * z(nt);
  std::vector<double> yy(m);
  svk::Mat xx(static_cast<int>(m), d.p);
  for (int j = 0; j < nc; ++j) {
    for (int t = 0; t < nt; ++t) yy[z(j) * z(nt) + z(t)] = d.yv(t, cu[z(j)]);
  }
  ddm(yy, nt, nc);
  for (size_t i = 0; i < m; ++i) yy[i] -= c[i];
  std::vector<double> col(m);
  for (int k = 0; k < d.p; ++k) {
    for (int j = 0; j < nc; ++j) {
      for (int t = 0; t < nt; ++t) col[z(j) * z(nt) + z(t)] = d.xv(t, cu[z(j)], k);
    }
    ddm(col, nt, nc);
    std::copy(col.begin(), col.end(), xx.v.begin() + static_cast<std::ptrdiff_t>(z(k) * m));
  }
  return svk::lsq(xx, yy, b);
}

// Bai's estimator on the controls cu with r factors.
int bai(const Pan &d, const std::vector<int> &cu, int r, int maxit, double tol, Fit &o) {
  const int nc = static_cast<int>(cu.size()), nt = d.nt;
  std::vector<double> w(z(nc) * z(nt)), c(w.size(), 0.0), b(z(d.p), 0.0);
  o.it = 0;
  o.conv = true;
  if (d.p > 0 && !bols(d, cu, c, b)) return CLAB_EDOM;
  for (int it = 0;; ++it) {
    if (kutil::stop(static_cast<double>(nt) * nt * nc)) return CLAB_EINTR;
    for (int j = 0; j < nc; ++j) {
      for (int t = 0; t < nt; ++t) w[z(j) * z(nt) + z(t)] = d.yv(t, cu[z(j)]) - xb(d, t, cu[z(j)], b);
    }
    ddm(w, nt, nc);
    comp(w, nt, nc, r, o.f, c);
    o.it = it + 1;
    if (d.p == 0 || r == 0) break;
    std::vector<double> b1;
    if (!bols(d, cu, c, b1)) return CLAB_EDOM;
    double dm = 0.0;
    for (int k = 0; k < d.p; ++k) dm = std::max(dm, std::fabs(b1[z(k)] - b[z(k)]));
    b = b1;
    if (dm < tol) break;
    if (it + 1 >= maxit) {
      o.conv = false;
      break;
    }
  }
  o.beta = b;
  // mu and the time effects of the residuals net of covariates
  o.xi.assign(z(nt), 0.0);
  double g = 0.0;
  for (int j = 0; j < nc; ++j) {
    for (int t = 0; t < nt; ++t) {
      const double v = d.yv(t, cu[z(j)]) - xb(d, t, cu[z(j)], b);
      o.xi[z(t)] += v / nc;
      g += v;
    }
  }
  o.mu = g / (static_cast<double>(nc) * nt);
  for (double &v : o.xi) v -= o.mu;
  return CLAB_OK;
}

// A treated unit's intercept and loadings by least squares of its
// pre-adoption residuals on [1, F], period s left out (-1 none); false
// when the periods do not identify them.
bool load(const Pan &d, const Fit &f, int u, int r, int s, std::vector<double> &g) {
  const int a = d.a[u];
  const int m = a - (s >= 0 ? 1 : 0);
  if (m < r + 1) return false;
  svk::Mat xx(m, r + 1);
  std::vector<double> yy(z(m));
  int i = 0;
  for (int t = 0; t < a; ++t) {
    if (t == s) continue;
    yy[z(i)] = d.yv(t, u) - xb(d, t, u, f.beta) - f.mu - f.xi[z(t)];
    xx(i, 0) = 1.0;
    for (int k = 0; k < r; ++k) xx(i, k + 1) = f.f[z(k) * z(d.nt) + z(t)];
    ++i;
  }
  return svk::lsq(xx, yy, g);
}

double pred(const Pan &d, const Fit &f, int u, int r, int t, const std::vector<double> &g) {
  double v = xb(d, t, u, f.beta) + f.mu + f.xi[z(t)] + g[0];
  for (int k = 0; k < r; ++k) v += f.f[z(k) * z(d.nt) + z(t)] * g[z(k) + 1];
  return v;
}

// The mean squared error of predicting each treated unit's pre-adoption
// periods left out one at a time; NaN when some unit has too few.
double cvr(const Pan &d, const std::vector<int> &tu, int r, const Fit &f) {
  double s = 0.0, m = 0.0;
  std::vector<double> g;
  for (int u : tu) {
    for (int t = 0; t < d.a[u]; ++t) {
      if (!load(d, f, u, r, t, g)) return NAN;
      const double e = d.yv(t, u) - pred(d, f, u, r, t, g);
      s += e * e;
      m += 1.0;
    }
  }
  return s / m;
}

// The effects of a fit: the counterfactuals of the treated units into
// y0 (when not NULL), the mean effect over treated cells and over each
// event time e - emin (ne of them) into ef (1 + ne); false when a
// treated unit's loadings are not identified.
bool eff(const Pan &d, const std::vector<int> &tu, int r, const Fit &f, int emin, int ne, double *y0, double *ef) {
  std::vector<double> s(z(ne) + 1, 0.0), m(z(ne) + 1, 0.0), g;
  for (int u : tu) {
    if (!load(d, f, u, r, -1, g)) return false;
    for (int t = 0; t < d.nt; ++t) {
      const double v = pred(d, f, u, r, t, g);
      if (y0) y0[z(u) * z(d.nt) + z(t)] = v;
      const double e = d.yv(t, u) - v;
      const size_t k = z(t - d.a[u] - emin) + 1;
      s[k] += e;
      m[k] += 1.0;
      if (t >= d.a[u]) {
        s[0] += e;
        m[0] += 1.0;
      }
    }
  }
  for (size_t k = 0; k <= z(ne); ++k) ef[k] = m[k] > 0 ? s[k] / m[k] : NAN;
  return true;
}

}  // namespace

extern "C" int clab_ifesc(const double *y, const double *x, int p, const int *a, int n, int nt, int r, int rmax,
                          int maxit, double tol, const int *bm, int nb, double *beta, double *y0, double *fac,
                          double *cv, double *ef, double *bef, double *st) CLAB_TRY {
  if (n < 2 || nt < 2 || p < 0 || r < -1 || rmax < 0 || maxit < 1 || !(tol > 0.0) || nb < 0 || (nb > 0 && !bm)) {
    return CLAB_EARG;
  }
  if (static_cast<double>(n) * nt * (p + 1) > kutil::maxcell || static_cast<double>(nt) * nt > kutil::maxcell) {
    return CLAB_EARG;
  }
  for (size_t i = 0; i < z(n) * z(nt) * z(p + 1); ++i) {
    if (!std::isfinite(i < z(n) * z(nt) ? y[i] : x[i - z(n) * z(nt)])) return CLAB_EDOM;
  }
  const Pan d{y, x, p, n, nt, a};
  std::vector<int> cu, tu;
  int amin = nt, amax = 0;
  for (int u = 0; u < n; ++u) {
    if (a[u] < 0 || a[u] > nt) return CLAB_EARG;
    if (a[u] == nt) {
      cu.push_back(u);
    } else {
      tu.push_back(u);
      amin = std::min(amin, a[u]);
      amax = std::max(amax, a[u]);
    }
  }
  if (tu.empty() || cu.size() < 2) return CLAB_EARG;
  const int emin = -amax, ne = nt - amin + amax;
  // r: given, or the smallest mean squared error from 0 to rmax
  Fit f;
  int rr = r;
  std::fill(cv, cv + rmax + 1, NAN);
  if (r < 0) {
    double best = INFINITY;
    for (int k = 0; k <= rmax && k < static_cast<int>(cu.size()) && k < nt; ++k) {
      const int rc = bai(d, cu, k, maxit, tol, f);
      if (rc != CLAB_OK) return rc;
      cv[k] = cvr(d, tu, k, f);
      if (std::isnan(cv[k])) break;
      if (cv[k] < best) {
        best = cv[k];
        rr = k;
      }
    }
    if (rr < 0) return CLAB_EDOM;
  }
  if (rr >= static_cast<int>(cu.size()) || rr >= nt) return CLAB_EARG;
  const int rc = bai(d, cu, rr, maxit, tol, f);
  if (rc != CLAB_OK) return rc;
  std::fill(y0, y0 + z(n) * z(nt), NAN);
  if (!eff(d, tu, rr, f, emin, ne, y0, ef)) return CLAB_EDOM;
  for (int k = 0; k < p; ++k) beta[k] = f.beta[z(k)];
  for (size_t i = 0; i < z(nt) * z(rr); ++i) fac[i] = f.f[i];
  // the bootstrap: units drawn with the counts of bm (n x nb)
  int nf = 0;
  std::vector<double> e1(z(ne) + 1);
  for (int b = 0; b < nb; ++b) {
    std::vector<int> c1, t1;
    for (int u = 0; u < n; ++u) {
      for (int m = 0; m < bm[z(b) * z(n) + z(u)]; ++m) (a[u] == nt ? c1 : t1).push_back(u);
    }
    Fit g;
    bool ok = !t1.empty() && static_cast<int>(c1.size()) > rr && c1.size() >= 2;
    if (ok) {
      const int r1 = bai(d, c1, rr, maxit, tol, g);
      if (r1 == CLAB_EINTR) return r1;
      ok = r1 == CLAB_OK && eff(d, t1, rr, g, emin, ne, nullptr, e1.data());
    }
    for (size_t k = 0; k <= z(ne); ++k) bef[k * z(nb) + z(b)] = ok ? e1[k] : NAN;
    nf += ok ? 0 : 1;
  }
  st[0] = rr;
  st[1] = f.it;
  st[2] = f.conv ? 1 : 0;
  st[3] = static_cast<double>(tu.size());
  st[4] = static_cast<double>(cu.size());
  st[5] = nf;
  st[6] = emin;
  st[7] = ne;
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)
