// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, geostatistics: the empirical semivariogram, variogram models
// fitted by Gaussian maximum likelihood and by weighted least squares,
// ordinary point and block kriging with leave-one-out cross-validation,
// polynomial trend surfaces, the spatial linear model with exponential
// covariance by profile maximum likelihood and the spatial linear mixed
// model by REML, the Web Mercator projection and density-weighted
// sampling of points.
#include <algorithm>
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "gla.h"
#include "kutil.h"
#include "pev.h"

namespace {

using vd = std::vector<double>;
using gla::at;

constexpr double pi = 3.14159265358979323846;

double dxy(const double *x, const double *y, int i, int j) {
  const double a = x[i] - x[j], b = y[i] - y[j];
  return std::sqrt(a * a + b * b);
}

// the n x n distance matrix of the points
vd dmat(const double *x, const double *y, int n) {
  vd d(static_cast<std::size_t>(n) * n, 0.0);
  for (int j = 0; j < n; ++j)
    for (int i = j + 1; i < n; ++i) d[at(i, j, n)] = d[at(j, i, n)] = dxy(x, y, i, j);
  return d;
}

// the semivariance of a model (0 spherical, 1 exponential, 2 gaussian) at
// distance h; 0 at h = 0
double vgam(double h, int model, double nug, double ps, double rg) {
  if (h == 0) return 0;
  double s;
  if (model == 0) {
    const double u = h / rg;
    s = h >= rg ? 1 : 1.5 * u - 0.5 * u * u * u;
  } else if (model == 1) {
    s = 1 - std::exp(-h / rg);
  } else {
    const double u = h / rg;
    s = 1 - std::exp(-(u * u));
  }
  return nug + ps * s;
}

// the exponential covariance of the block methods: c1 exp(-h / a), plus
// the nugget c0 at h = 0
double ecov(double h, double c0, double c1, double a) { return c1 * std::exp(-h / a) + (h == 0 ? c0 : 0.0); }

bool okmodel(int model) { return model >= 0 && model <= 2; }

bool okpts(const double *x, const double *y, int n) {
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(x[i]) || !std::isfinite(y[i])) return false;
  }
  return true;
}

double mean(const double *v, int n) {
  kutil::ksum s;
  for (int i = 0; i < n; ++i) s.add(v[i]);
  return s.val() / n;
}

// The binned semivariogram (as clab_vgemp) into dist, gam, np; *nk bins.
void bins(const double *x, const double *y, int n, const double *v, int nb, double cut, double *dist,
          double *gam, double *np, int *nk) {
  vd br(static_cast<std::size_t>(nb) + 1);
  const double by = cut / nb;
  br[0] = 0;
  for (int k = 1; k < nb; ++k) br[k] = k * by;
  br[nb] = cut;
  std::vector<kutil::ksum> s(static_cast<std::size_t>(nb));
  vd cnt(static_cast<std::size_t>(nb), 0.0);
  for (int i = 0; i < n; ++i)
    for (int j = i + 1; j < n; ++j) {
      const double d = dxy(x, y, i, j);
      if (!(d > 0 && d <= cut)) continue;
      // the bin (br[k - 1], br[k]], the first closed below
      const int k = static_cast<int>(std::lower_bound(br.begin() + 1, br.end(), d) - br.begin()) - 1;
      const double e = v[i] - v[j];
      s[k].add(e * e / 2);
      cnt[k] += 1;
    }
  int m = 0;
  for (int k = 0; k < nb; ++k) {
    if (cnt[k] == 0) continue;
    dist[m] = (br[k + 1] + br[k]) / 2;
    gam[m] = s[k].val() / cnt[k];
    np[m] = cnt[k];
    ++m;
  }
  *nk = m;
}

double maxd(const double *x, const double *y, int n) {
  double m = 0;
  for (int i = 0; i < n; ++i)
    for (int j = i + 1; j < n; ++j) m = std::max(m, dxy(x, y, i, j));
  return m;
}

// Generalised least squares of v on the p columns of z with the
// covariance whose lower Cholesky factor is l: beta, (Z' S^-1 Z)^-1 into
// xi (may be NULL) and the whitened residual sum of squares. False when
// Z' S^-1 Z is singular.
bool gls(const double *l, int n, const double *z, int p, const double *v, double *beta, double *xi, double *rss,
         double *ld2) {
  vd zw(z, z + static_cast<std::size_t>(n) * p), yw(v, v + n);
  gla::fsolve(l, n, zw.data(), p);
  gla::fsolve(l, n, yw.data(), 1);
  vd a(static_cast<std::size_t>(p) * p), b(static_cast<std::size_t>(p));
  for (int j = 0; j < p; ++j) {
    for (int k = 0; k < p; ++k) {
      double s = 0;
      for (int i = 0; i < n; ++i) s += zw[at(i, j, n)] * zw[at(i, k, n)];
      a[at(j, k, p)] = s;
    }
    double s = 0;
    for (int i = 0; i < n; ++i) s += zw[at(i, j, n)] * yw[i];
    b[j] = s;
  }
  vd lc(a);
  if (!pev::chol(lc.data(), p)) return false;
  if (ld2 != nullptr) {
    double s = 0;
    for (int j = 0; j < p; ++j) s += std::log(lc[at(j, j, p)]);
    *ld2 = 2 * s;
  }
  vd sol(b);
  gla::fsolve(lc.data(), p, sol.data(), 1);
  gla::bsolve(lc.data(), p, sol.data(), 1);
  for (int j = 0; j < p; ++j) beta[j] = sol[j];
  if (xi != nullptr) {
    vd e(static_cast<std::size_t>(p) * p, 0.0);
    for (int j = 0; j < p; ++j) e[at(j, j, p)] = 1;
    gla::fsolve(lc.data(), p, e.data(), p);
    gla::bsolve(lc.data(), p, e.data(), p);
    std::copy(e.begin(), e.end(), xi);
  }
  double s = 0;
  for (int i = 0; i < n; ++i) {
    double r = yw[i];
    for (int j = 0; j < p; ++j) r -= zw[at(i, j, n)] * beta[j];
    s += r * r;
  }
  *rss = s;
  return true;
}

bool okz(const double *z, int n, int p, const double *v) {
  for (std::size_t q = 0; q < static_cast<std::size_t>(n) * p; ++q) {
    if (!std::isfinite(z[q])) return false;
  }
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(v[i])) return false;
  }
  return p >= 1 && n > p;
}

}  // namespace

int clab_vgemp(const double *x, const double *y, int n, const double *v, int nb, double cut, double *dist,
               double *gam, double *np, int *nk) {
  if (n < 2 || nb < 1 || !okpts(x, y, n)) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(v[i])) return CLAB_EDOM;
  }
  if (std::isnan(cut)) cut = maxd(x, y, n) / 3;
  if (!kutil::pos(cut)) return CLAB_EDOM;
  bins(x, y, n, v, nb, cut, dist, gam, np, nk);
  return CLAB_OK;
}

int clab_vgfit(const double *x, const double *y, int n, const double *v, int model, double *out) {
  if (n < 3 || !okmodel(model) || !okpts(x, y, n)) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(v[i])) return CLAB_EDOM;
  }
  const vd D = dmat(x, y, n);
  const double mu = mean(v, n);
  double ss = 0;
  for (int i = 0; i < n; ++i) ss += (v[i] - mu) * (v[i] - mu);
  const double v0 = ss / (n - 1);
  const double r0 = *std::max_element(D.begin(), D.end()) / 4;
  if (!kutil::pos(v0) || !kutil::pos(r0)) return CLAB_EDOM;
  vd C(D.size()), r(static_cast<std::size_t>(n));
  auto nll = [&](const double *p) {
    const double nug = std::exp(p[0]), ps = std::exp(p[1]), rg = std::exp(p[2]);
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < n; ++i) C[at(i, j, n)] = (nug + ps) - vgam(D[at(i, j, n)], model, nug, ps, rg);
    for (int i = 0; i < n; ++i) C[at(i, i, n)] += 1e-8 * v0;
    if (!gla::chol(C.data(), n)) return 1e10;
    for (int i = 0; i < n; ++i) r[i] = v[i] - mu;
    gla::fsolve(C.data(), n, r.data(), 1);
    double ld = 0, q = 0;
    for (int i = 0; i < n; ++i) {
      ld += std::log(C[at(i, i, n)]);
      q += r[i] * r[i];
    }
    return ld + 0.5 * q;
  };
  double p[3] = {std::log(v0 * 0.1), std::log(v0 * 0.9), std::log(r0)}, f = 0;
  int nf = 0;
  const int rc = gla::nmmin(3, p, &f, nll, 1.4901161193847656e-08, 800, &nf);
  if (kutil::stop(0)) return CLAB_EINTR;
  out[0] = std::exp(p[0]);
  out[1] = std::exp(p[1]);
  out[2] = std::exp(p[2]);
  out[3] = -f;
  out[4] = rc == 0 ? 1 : 0;
  out[5] = nf;
  return CLAB_OK;
}

int clab_vgwls(const double *x, const double *y, int n, const double *v, int model, double *out) {
  if (n < 2 || !okmodel(model) || !okpts(x, y, n)) return CLAB_EARG;
  double dist[15], gam[15], np[15];
  int nk = 0;
  const int rc = clab_vgemp(x, y, n, v, 15, kutil::nan(), dist, gam, np, &nk);
  if (rc != CLAB_OK) return rc;
  double v0 = 0, hm = 0;
  for (int k = 0; k < nk; ++k) {
    v0 = std::max(v0, gam[k]);
    hm = std::max(hm, dist[k]);
  }
  auto obj = [&](const double *p) {
    const double nug = std::exp(p[0]), ps = std::exp(p[1]), rg = std::exp(p[2]);
    double s = 0;
    for (int k = 0; k < nk; ++k) {
      const double w = np[k] / std::max(dist[k] * dist[k], 1e-12), e = gam[k] - vgam(dist[k], model, nug, ps, rg);
      s += w * e * e;
    }
    return s;
  };
  double p[3] = {std::log(v0 * 0.1 + 1e-8), std::log(v0 * 0.9 + 1e-8), std::log(hm / 3)}, f = 0;
  int nf = 0;
  gla::nmmin(3, p, &f, obj, 1.4901161193847656e-08, 500, &nf);
  if (kutil::stop(0)) return CLAB_EINTR;
  for (int k = 0; k < 3; ++k) out[k] = std::exp(p[k]);
  return CLAB_OK;
}

int clab_okrge(const double *x, const double *y, int n, const double *v, const double *ax, const double *ay,
               int m, int model, double nug, double ps, double rg, double *pred, double *var) {
  if (n < 1 || m < 1 || !okmodel(model) || !okpts(x, y, n) || !okpts(ax, ay, m)) return CLAB_EARG;
  if (!kutil::nonneg(nug) || !kutil::nonneg(ps) || !kutil::pos(rg)) return CLAB_EDOM;
  const int k = n + 1;
  vd A(static_cast<std::size_t>(k) * k, 0.0), B(static_cast<std::size_t>(k) * m);
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) A[at(i, j, k)] = vgam(dxy(x, y, i, j), model, nug, ps, rg);
    A[at(n, j, k)] = A[at(j, n, k)] = 1;
  }
  for (int q = 0; q < m; ++q) {
    for (int i = 0; i < n; ++i) {
      const double a = x[i] - ax[q], b = y[i] - ay[q];
      B[at(i, q, k)] = vgam(std::sqrt(a * a + b * b), model, nug, ps, rg);
    }
    B[at(n, q, k)] = 1;
  }
  vd W(B);
  std::vector<int> piv(static_cast<std::size_t>(k));
  if (!gla::lusolve(A.data(), k, W.data(), m, piv.data())) return CLAB_EDOM;
  for (int q = 0; q < m; ++q) {
    double s = 0, t = 0;
    for (int i = 0; i < n; ++i) s += W[at(i, q, k)] * v[i];
    for (int i = 0; i < k; ++i) t += W[at(i, q, k)] * B[at(i, q, k)];
    pred[q] = s;
    var[q] = std::max(t, 0.0);
  }
  return CLAB_OK;
}

int clab_bkpts(const double *lo, const double *hi, int nq, double *px, double *py) {
  if (nq < 1) return CLAB_EARG;
  const int k = static_cast<int>(std::nearbyint(std::sqrt(static_cast<double>(nq))));
  if (static_cast<double>(k) * k > kutil::maxcell) return CLAB_EARG;
  // seq(lo, hi, length.out = k) on each axis, the first axis fastest
  auto seq = [&](double a, double b, int i) {
    if (k == 1) return a;
    if (i == k - 1) return b;
    return a + i * ((b - a) / (k - 1));
  };
  for (int j = 0; j < k; ++j)
    for (int i = 0; i < k; ++i) {
      px[j * k + i] = seq(lo[0], hi[0], i);
      py[j * k + i] = seq(lo[1], hi[1], j);
    }
  return CLAB_OK;
}

int clab_bkrig(const double *x, const double *y, int n, const double *v, const double *bx, const double *by,
               const int *off, int nblk, double nug, double sill, double rg, double *est, double *se) {
  if (n < 1 || nblk < 1 || !okpts(x, y, n) || off[nblk] > kutil::maxcell) return CLAB_EARG;
  if (!kutil::nonneg(nug) || !std::isfinite(sill) || !kutil::pos(rg)) return CLAB_EDOM;
  const double c1 = sill - nug;
  const int k = n + 1;
  vd A(static_cast<std::size_t>(k) * k, 0.0), R(static_cast<std::size_t>(k) * nblk), Cb(static_cast<std::size_t>(nblk));
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) A[at(i, j, k)] = ecov(dxy(x, y, i, j), nug, c1, rg);
    A[at(n, j, k)] = A[at(j, n, k)] = 1;
  }
  for (int b = 0; b < nblk; ++b) {
    const int q0 = off[b], q1 = off[b + 1], nq = q1 - q0;
    if (nq < 1 || !okpts(bx + q0, by + q0, nq)) return CLAB_EARG;
    for (int i = 0; i < n; ++i) {
      if (kutil::stop(nq)) return CLAB_EINTR;
      double s = 0;
      for (int q = q0; q < q1; ++q) {
        const double a = x[i] - bx[q], c = y[i] - by[q];
        s += ecov(std::sqrt(a * a + c * c), nug, c1, rg);
      }
      R[at(i, b, k)] = s / nq;
    }
    R[at(n, b, k)] = 1;
    double s = 0;
    for (int q = q0; q < q1; ++q) {
      if (kutil::stop(nq)) return CLAB_EINTR;
      for (int r = q0; r < q1; ++r) s += ecov(dxy(bx, by, q, r), nug, c1, rg);
    }
    Cb[b] = s / (static_cast<double>(nq) * nq);
  }
  vd S(R);
  std::vector<int> piv(static_cast<std::size_t>(k));
  if (!gla::lusolve(A.data(), k, S.data(), nblk, piv.data())) return CLAB_EDOM;
  for (int b = 0; b < nblk; ++b) {
    double e = 0, l = 0;
    for (int i = 0; i < n; ++i) {
      e += S[at(i, b, k)] * v[i];
      l += S[at(i, b, k)] * R[at(i, b, k)];
    }
    est[b] = e;
    se[b] = std::sqrt(std::max(Cb[b] - l - S[at(n, b, k)], 0.0));
  }
  return CLAB_OK;
}

int clab_krigcv(const double *x, const double *y, int n, const double *v, double nug, double sill, double rg,
                double *res) {
  if (n < 2 || !okpts(x, y, n)) return CLAB_EARG;
  if (!kutil::nonneg(nug) || !std::isfinite(sill) || !kutil::pos(rg)) return CLAB_EDOM;
  const double c1 = sill - nug;
  const int k = n;  // n - 1 points and the constraint
  vd A(static_cast<std::size_t>(k) * k), b(static_cast<std::size_t>(k));
  std::vector<int> piv(static_cast<std::size_t>(k));
  for (int o = 0; o < n; ++o) {
    if (kutil::stop(static_cast<double>(n) * n * n)) return CLAB_EINTR;
    int jj = 0;
    for (int j = 0; j < n; ++j) {
      if (j == o) continue;
      int ii = 0;
      for (int i = 0; i < n; ++i) {
        if (i == o) continue;
        A[at(ii++, jj, k)] = ecov(dxy(x, y, i, j), nug, c1, rg);
      }
      A[at(k - 1, jj, k)] = A[at(jj, k - 1, k)] = 1;
      b[jj] = ecov(dxy(x, y, o, j), nug, c1, rg);
      ++jj;
    }
    A[at(k - 1, k - 1, k)] = 0;
    b[k - 1] = 1;
    if (!gla::lusolve(A.data(), k, b.data(), 1, piv.data())) return CLAB_EDOM;
    double z = 0;
    jj = 0;
    for (int j = 0; j < n; ++j) {
      if (j != o) z += b[jj++] * v[j];
    }
    res[o] = v[o] - z;
  }
  return CLAB_OK;
}

int clab_trsurf(const double *x, const double *y, int n, const double *v, int order, double *coef, double *se,
                double *r2) {
  if (order < 0 || order > 3 || !okpts(x, y, n)) return CLAB_EARG;
  static const int np[] = {1, 3, 6, 10};
  const int p = np[order];
  if (n < p) return CLAB_EARG;
  vd F(static_cast<std::size_t>(n) * p);
  for (int i = 0; i < n; ++i) {
    const double s = x[i], t = y[i];
    const double c[10] = {1, s, t, s * s, t * t, s * t, s * s * s, t * t * t, s * s * t, s * t * t};
    for (int j = 0; j < p; ++j) F[at(i, j, n)] = c[j];
  }
  vd xi(static_cast<std::size_t>(p) * p), r(static_cast<std::size_t>(n)),
      work(static_cast<std::size_t>(p) * p + p);
  if (pev::ols(F.data(), v, n, p, nullptr, coef, xi.data(), r.data(), work.data()) != CLAB_OK) return CLAB_EDOM;
  double rss = 0;
  for (int i = 0; i < n; ++i) rss += r[i] * r[i];
  const double s2 = rss / std::max(n - p, 1);
  for (int j = 0; j < p; ++j) se[j] = std::sqrt(std::max(s2 * xi[at(j, j, p)], 0.0));
  const double mu = mean(v, n);
  double tss = 0;
  for (int i = 0; i < n; ++i) tss += (v[i] - mu) * (v[i] - mu);
  *r2 = tss > 0 ? 1 - rss / tss : 1;
  return CLAB_OK;
}

int clab_sglmx(const double *z, int n, int p, const double *v, const double *x, const double *y, double *coef,
               double *se, double *out) {
  if (n < 2 || !okpts(x, y, n)) return CLAB_EARG;
  if (!okz(z, n, p, v)) return CLAB_EDOM;
  const vd D = dmat(x, y, n);
  const double hm = *std::max_element(D.begin(), D.end());
  vd L(D.size()), b(static_cast<std::size_t>(p));
  auto fill = [&](double phi) {
    for (std::size_t q = 0; q < D.size(); ++q) L[q] = std::exp(-D[q] / phi);
    for (int i = 0; i < n; ++i) L[at(i, i, n)] += 1e-8;
    return gla::chol(L.data(), n);
  };
  auto nll = [&](double lp) {
    if (!fill(std::exp(lp))) return 1e12;
    double rss = 0;
    if (!gls(L.data(), n, z, p, v, b.data(), nullptr, &rss, nullptr)) return 1e12;
    double ld = 0;
    for (int i = 0; i < n; ++i) ld += std::log(L[at(i, i, n)]);
    return 0.5 * (n * std::log(2 * pi * (rss / n)) + 2 * ld + n);
  };
  const double lo = std::log(std::max(hm / 100, 1e-3)), hi = std::log(std::max(hm * 3, 1.0));
  const double lp = gla::fmin(nll, lo, hi, std::pow(DBL_EPSILON, 0.25));
  if (kutil::stop(0)) return CLAB_EINTR;
  const double phi = std::exp(lp);
  vd xi(static_cast<std::size_t>(p) * p);
  double rss = 0;
  if (!fill(phi) || !gls(L.data(), n, z, p, v, coef, xi.data(), &rss, nullptr)) return CLAB_EDOM;
  const double s2 = rss / std::max(n - p, 1);
  for (int j = 0; j < p; ++j) se[j] = std::sqrt(s2 * xi[at(j, j, p)]);
  out[0] = s2;
  out[1] = phi;
  out[2] = 0;
  return CLAB_OK;
}

int clab_smixm(const double *z, int n, int p, const double *v, const double *x, const double *y, double *coef,
               double *se, double *out) {
  if (n < 2 || !okpts(x, y, n)) return CLAB_EARG;
  if (!okz(z, n, p, v)) return CLAB_EDOM;
  const vd D = dmat(x, y, n);
  const double hm = std::max(*std::max_element(D.begin(), D.end()), 1.0);
  vd L(D.size()), b(static_cast<std::size_t>(p)), yw(static_cast<std::size_t>(n));
  auto fill = [&](double phi, double nu) {
    for (std::size_t q = 0; q < D.size(); ++q) L[q] = std::exp(-D[q] / phi);
    for (int i = 0; i < n; ++i) L[at(i, i, n)] += nu;
    return gla::chol(L.data(), n);
  };
  auto nreml = [&](const double *th) {
    if (!fill(std::exp(th[0]), std::exp(th[1]))) return 1e12;
    double rss = 0, ldk = 0;
    if (!gls(L.data(), n, z, p, v, b.data(), nullptr, &rss, &ldk)) return 1e12;
    const double s2 = rss / std::max(n - p, 1);
    std::copy(v, v + n, yw.begin());
    gla::fsolve(L.data(), n, yw.data(), 1);
    double my = 0;
    for (int i = 0; i < n; ++i) my += yw[i] * yw[i];
    const double tol = std::sqrt(DBL_EPSILON) * std::max(1.0, my / n);
    if (!std::isfinite(s2) || s2 <= tol) return 1e12;
    double ld = 0;
    for (int i = 0; i < n; ++i) ld += std::log(L[at(i, i, n)]);
    return 0.5 * (2 * ld + ldk + (n - p) * std::log(2 * pi * s2) + (n - p));
  };
  double th[2] = {std::log(hm / 3), std::log(0.1)}, f = 0;
  int nf = 0;
  const int rc = gla::nmmin(2, th, &f, nreml, 1e-6, 400, &nf);
  if (kutil::stop(0)) return CLAB_EINTR;
  const double phi = std::exp(th[0]), nu = std::exp(th[1]);
  vd xi(static_cast<std::size_t>(p) * p);
  double rss = 0;
  if (!fill(phi, nu) || !gls(L.data(), n, z, p, v, coef, xi.data(), &rss, nullptr)) return CLAB_EDOM;
  const double s2 = rss / std::max(n - p, 1);
  for (int j = 0; j < p; ++j) se[j] = std::sqrt(s2 * xi[at(j, j, p)]);
  out[0] = s2;
  out[1] = nu * s2;
  out[2] = phi;
  out[3] = rc == 0 ? 1 : 0;
  return CLAB_OK;
}

int clab_wmerc(const double *x, const double *y, int n, int inv, double *ox, double *oy) {
  const double R = 6378137;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(x[i]) || !std::isfinite(y[i])) return CLAB_EDOM;
    if (inv) {
      ox[i] = (x[i] / R) * 180 / pi;
      oy[i] = (2 * std::atan(std::exp(y[i] / R)) - pi / 2) * 180 / pi;
    } else {
      if (!(std::fabs(y[i]) < 90)) return CLAB_EDOM;
      ox[i] = R * x[i] * pi / 180;
      oy[i] = R * std::log(std::tan(pi / 4 + (y[i] * pi / 180) / 2));
    }
  }
  return CLAB_OK;
}

int clab_spdens(const double *x, const double *y, int n, int ll, double *prob) {
  if (n < 2 || !okpts(x, y, n)) return CLAB_EARG;
  vd D(static_cast<std::size_t>(n) * n, 0.0);
  for (int j = 0; j < n; ++j)
    for (int i = j + 1; i < n; ++i) {
      double d;
      if (ll) {
        if (clab_hvdist(x + i, y + i, x + j, y + j, 1, 6371008.8, &d) != CLAB_OK) return CLAB_EDOM;
      } else {
        d = dxy(x, y, i, j);
      }
      D[at(i, j, n)] = D[at(j, i, n)] = d;
    }
  // the bandwidth: the median of the positive distances (each pair twice)
  vd pd;
  for (double d : D) {
    if (d > 0) pd.push_back(d);
  }
  if (pd.empty()) return CLAB_EDOM;
  std::sort(pd.begin(), pd.end());
  const std::size_t h = pd.size() / 2;
  const double bw = pd.size() % 2 ? pd[h] : (pd[h - 1] + pd[h]) / 2;
  kutil::ksum tot;
  for (int i = 0; i < n; ++i) {
    kutil::ksum s;
    for (int j = 0; j < n; ++j) {
      const double u = D[at(i, j, n)] / bw;
      s.add(std::exp(-(u * u)));
    }
    prob[i] = s.val();
    tot.add(prob[i]);
  }
  const double t = tot.val();
  for (int i = 0; i < n; ++i) prob[i] /= t;
  return CLAB_OK;
}

int clab_spsamp(const double *prob, int n, int size, unsigned int seed, int *idx) {
  if (n < 1 || size < 0 || size > n) return CLAB_EARG;
  kutil::rng g(seed);
  vd p(static_cast<std::size_t>(n), 1.0);
  if (prob != nullptr) {
    for (int i = 0; i < n; ++i) {
      if (!kutil::nonneg(prob[i])) return CLAB_EDOM;
      p[i] = prob[i];
    }
  }
  std::vector<int> left(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) left[i] = i;
  for (int k = 0; k < size; ++k) {
    // draw from the units left in proportion to their weight; a set of
    // weights all 0 is drawn from uniformly
    const int m = n - k;
    double tot = 0;
    for (int q = 0; q < m; ++q) tot += p[left[q]];
    int pick = m - 1;
    if (tot > 0) {
      const double u = g.unif() * tot;
      double c = 0;
      for (int q = 0; q < m; ++q) {
        c += p[left[q]];
        if (u < c) {
          pick = q;
          break;
        }
      }
    } else {
      pick = g.below(m);
    }
    idx[k] = left[pick] + 1;
    left[pick] = left[m - 1];
  }
  return CLAB_OK;
}
