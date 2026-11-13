// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, learning: elastic net paths by cyclical coordinate descent
// (Friedman, Hastie and Tibshirani 2010, doi:10.18637/jss.v033.i01), with
// covariance updates for the gaussian family and a proximal Newton outer
// loop with active-set cycling and a check of the optimality conditions
// for the binomial one; K-fold cross-validation over caller-supplied
// folds; generalised linear models by Newton steps. Rows are used through index lists, never copied.
#include "learn.h"

#include <cmath>

#include "crimlab.h"

namespace {

// The soft-thresholding operator: the proximal map of t |.| (Lean
// candidate: soft_threshold_prox).
inline double soft(double z, double t) {
  if (z > t) return z - t;
  if (z < -t) return z + t;
  return 0.0;
}

// The weighted means and population standard deviations of the columns
// of x over the rows idx[0..m), each row weighing 1 / m. The second pass
// sums centred squares, so a column of large mean loses no digits.
void stdz(const lrn::Dat &d, double *mu, double *sd) {
  const double w = 1.0 / d.m;
  for (int j = 0; j < d.p; ++j) {
    const double *xj = d.x + static_cast<long>(j) * d.n;
    double s = 0.0;
    for (int i = 0; i < d.m; ++i) s += xj[d.idx[i]];
    const double a = s * w;
    double v = 0.0;
    for (int i = 0; i < d.m; ++i) {
      const double e = xj[d.idx[i]] - a;
      v += e * e;
    }
    mu[j] = a;
    sd[j] = std::sqrt(v * w);
  }
}

// mean and population sd of y over the rows used
void ystat(const lrn::Dat &d, double *ym, double *ys) {
  double s = 0.0;
  for (int i = 0; i < d.m; ++i) s += d.y[d.idx[i]];
  const double a = s / d.m;
  double v = 0.0;
  for (int i = 0; i < d.m; ++i) {
    const double e = d.y[d.idx[i]] - a;
    v += e * e;
  }
  *ym = a;
  *ys = std::sqrt(v / d.m);
}

bool binary(const lrn::Dat &d) {
  for (int i = 0; i < d.m; ++i) {
    const double v = d.y[d.idx[i]];
    if (v != 0.0 && v != 1.0) return false;
  }
  return true;
}

// sum_i r[i] (x[idx[i], j] - mu) / sd over the rows used
inline double xdot(const lrn::Dat &d, int j, double mu, double sd,
                   const double *r) {
  const double *xj = d.x + static_cast<long>(j) * d.n;
  double s = 0.0;
  for (int i = 0; i < d.m; ++i) s += r[i] * (xj[d.idx[i]] - mu);
  return s / sd;
}

// sum_i a[i] b[i], in four running sums (independent, so they pipeline)
inline double dot(const double *a, const double *b, int m) {
  double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
  int i = 0;
  for (; i + 3 < m; i += 4) {
    s0 += a[i] * b[i];
    s1 += a[i + 1] * b[i + 1];
    s2 += a[i + 2] * b[i + 2];
    s3 += a[i + 3] * b[i + 3];
  }
  for (; i < m; ++i) s0 += a[i] * b[i];
  return (s0 + s1) + (s2 + s3);
}

// The standardised columns of the rows used, packed: xs[i + j m] =
// (x[idx[i], j] - mu[j]) / sd[j], zero for a constant column. Every pass
// of the descent then reads contiguous memory.
void pack(const lrn::Dat &d, const double *mu, const double *sd, double *xs) {
  for (int j = 0; j < d.p; ++j) {
    const double *xj = d.x + static_cast<long>(j) * d.n;
    double *o = xs + static_cast<long>(j) * d.m;
    const double f = sd[j] > 0.0 ? 1.0 / sd[j] : 0.0;
    for (int i = 0; i < d.m; ++i) o[i] = (xj[d.idx[i]] - mu[j]) * f;
  }
}

// glmnet's early stop of a path it chose itself: from the fifth lambda on,
// when the share of deviance explained (r, r0 at the lambda before) grew
// by less than 1e-5 of itself or passed 0.999
inline bool early(int k, double r, double r0) {
  return k >= 4 && (r - r0 < 1e-5 * r || r > 0.999);
}

// the gaussian path, covariance updates: g holds x~_j' (y~ - X~ b) / m,
// standardised y, and column j of C = X~' X~ / m is formed when b_j first
// moves. work: mu, sd, b, g, C (p p), t (n), xs (n p); iwork: have (p).
int gauss(const lrn::Dat &d, double alpha, const double *lam, int nl,
          double thresh, int maxit, int stop, int *nfit, double *a0,
          double *beta, double *work, int *iwork) {
  const int p = d.p, m = d.m;
  double *mu = work, *sd = mu + p, *b = sd + p, *g = b + p, *c = g + p;
  double *t = c + static_cast<long>(p) * p, *xs = t + d.n;
  int *have = iwork;
  double ym, ys;
  ystat(d, &ym, &ys);
  if (!(ys > 0.0)) return CLAB_EDOM;
  stdz(d, mu, sd);
  pack(d, mu, sd, xs);
  for (int i = 0; i < m; ++i) t[i] = (d.y[d.idx[i]] - ym) / (ys * m);
  for (int j = 0; j < p; ++j) {
    b[j] = 0.0;
    have[j] = 0;
    g[j] = dot(xs + static_cast<long>(j) * m, t, m);
  }
  // the columns of C for the group of four columns that holds j, in 4 x 4
  // tiles: one pass over eight columns gives sixteen inner products
  auto col = [&](int j) {
    const int j0 = j - j % 4, nj = p - j0 < 4 ? p - j0 : 4;
    for (int k0 = 0; k0 < p; k0 += 4) {
      const int nk = p - k0 < 4 ? p - k0 : 4;
      double acc[4][4] = {{0.0}};
      if (nj == 4 && nk == 4) {
        const double *a0 = xs + static_cast<long>(j0) * m, *a1 = a0 + m, *a2 = a1 + m,
                     *a3 = a2 + m;
        const double *b0 = xs + static_cast<long>(k0) * m, *b1 = b0 + m, *b2 = b1 + m,
                     *b3 = b2 + m;
        for (int i = 0; i < m; ++i) {
          const double x0 = a0[i], x1 = a1[i], x2 = a2[i], x3 = a3[i];
          const double y0 = b0[i], y1 = b1[i], y2 = b2[i], y3 = b3[i];
          acc[0][0] += x0 * y0; acc[0][1] += x0 * y1; acc[0][2] += x0 * y2; acc[0][3] += x0 * y3;
          acc[1][0] += x1 * y0; acc[1][1] += x1 * y1; acc[1][2] += x1 * y2; acc[1][3] += x1 * y3;
          acc[2][0] += x2 * y0; acc[2][1] += x2 * y1; acc[2][2] += x2 * y2; acc[2][3] += x2 * y3;
          acc[3][0] += x3 * y0; acc[3][1] += x3 * y1; acc[3][2] += x3 * y2; acc[3][3] += x3 * y3;
        }
      } else {
        for (int u = 0; u < nj; ++u) {
          for (int v = 0; v < nk; ++v) {
            acc[u][v] = dot(xs + static_cast<long>(j0 + u) * m,
                            xs + static_cast<long>(k0 + v) * m, m);
          }
        }
      }
      for (int u = 0; u < nj; ++u) {
        for (int v = 0; v < nk; ++v) c[static_cast<long>(j0 + u) * p + k0 + v] = acc[u][v] / m;
      }
    }
    for (int u = 0; u < nj; ++u) have[j0 + u] = 1;
  };
  // one coordinate: its change, squared; rsq, the share of the variance
  // of y explained, follows each change (y is standardised)
  double rsq = 0.0, rsq0 = 0.0;
  auto upd = [&](int j, double l1, double l2) {
    const double nb = soft(g[j] + b[j], l1) / (1.0 + l2);
    const double dl = nb - b[j];
    if (dl == 0.0) return 0.0;
    rsq += dl * (2.0 * g[j] - dl);
    if (!have[j]) col(j);
    const double *cj = c + static_cast<long>(j) * p;
    for (int k = 0; k < p; ++k) g[k] -= dl * cj[k];
    b[j] = nb;
    return dl * dl;
  };
  int it = 0;
  for (int k = 0; k < nl; ++k) {
    const double l1 = lam[k] / ys * alpha, l2 = lam[k] / ys * (1.0 - alpha);
    // a pass over every column, then passes over the nonzero ones until
    // they settle, then a pass over every column again, until one is quiet
    bool full = true;
    for (;;) {
      double dlx = 0.0;
      for (int j = 0; j < p; ++j) {
        if (full ? sd[j] > 0.0 : b[j] != 0.0) dlx = std::fmax(dlx, upd(j, l1, l2));
      }
      if (++it > maxit) return CLAB_EDOM;
      if (dlx >= thresh) {
        full = false;
      } else if (full) {
        break;
      } else {
        full = true;
      }
    }
    double s = ym, *bk = beta + static_cast<long>(k) * p;
    for (int j = 0; j < p; ++j) {
      bk[j] = b[j] != 0.0 ? b[j] * ys / sd[j] : 0.0;
      s -= mu[j] * bk[j];
    }
    a0[k] = s;
    *nfit = k + 1;
    if (stop && early(k, rsq, rsq0)) break;
    rsq0 = rsq;
  }
  return CLAB_OK;
}

// the binomial path: a proximal Newton outer loop (iteratively reweighted
// least squares) around coordinate descent over the active set, then a
// check of the optimality conditions over all columns, which adds the
// columns that violate them. work: mu, sd, b, g (unused), xv, bo (p
// each), eta, v, r, yy (n each), xs (n p); iwork: act (p).
int binom(const lrn::Dat &d, double alpha, const double *lam, int nl,
          double thresh, int maxit, int stop, int screen, int *nfit,
          double *a0, double *beta, double *work, int *iwork) {
  const int p = d.p, m = d.m;
  const double w = 1.0 / m;
  double *mu = work, *sd = mu + p, *b = sd + p, *xv = b + 2 * p;
  double *bo = xv + p, *eta = bo + p, *v = eta + d.n, *r = v + d.n;
  double *yy = r + d.n, *xs = yy + d.n;
  int *act = iwork;
  double *ga = b + p;  // |gradient| at the last check: the strong rule's input
  if (!binary(d)) return CLAB_EDOM;
  double ym, ys;
  ystat(d, &ym, &ys);
  if (!(ys > 0.0)) return CLAB_EDOM;
  stdz(d, mu, sd);
  pack(d, mu, sd, xs);
  for (int i = 0; i < m; ++i) yy[i] = d.y[d.idx[i]];
  double a = std::log(ym / (1.0 - ym));
  for (int i = 0; i < m; ++i) eta[i] = a;
  // the null deviance per row, over 2
  const double dev0 = -(ym * std::log(ym) + (1.0 - ym) * std::log1p(-ym));
  double r20 = 0.0;
  for (int i = 0; i < m; ++i) r[i] = w * (yy[i] - ym);
  double lp = 0.0;  // the lambda before, for the strong rule
  for (int j = 0; j < p; ++j) {
    b[j] = 0.0;
    act[j] = 0;
    ga[j] = std::fabs(dot(xs + static_cast<long>(j) * m, r, m));
    if (alpha > 0.0) lp = std::fmax(lp, ga[j] / alpha);
  }
  int it = 0;
  for (int k = 0; k < nl; ++k) {
    const double l1 = lam[k] * alpha, l2 = lam[k] * (1.0 - alpha);
    // sequential strong rule (Tibshirani et al. 2012): columns likely to
    // enter join the active set before the descent; the optimality check
    // below still covers every column
    if (screen) {
      const double cut = alpha * (2.0 * lam[k] - (k > 0 ? lam[k - 1] : std::fmax(lp, lam[k])));
      for (int j = 0; j < p; ++j) {
        if (!act[j] && sd[j] > 0.0 && ga[j] >= cut) act[j] = 1;
      }
    }
    for (;;) {
      // Newton steps on the active set
      for (;;) {
        double sv = 0.0;
        for (int i = 0; i < m; ++i) {
          const double pr = 1.0 / (1.0 + std::exp(-eta[i]));
          v[i] = w * pr * (1.0 - pr);
          r[i] = w * (yy[i] - pr);
          sv += v[i];
        }
        for (int j = 0; j < p; ++j) {
          bo[j] = b[j];
          if (!act[j]) continue;
          // sum v x^2, in four running sums
          const double *xj = xs + static_cast<long>(j) * m;
          double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
          int i = 0;
          for (; i + 3 < m; i += 4) {
            s0 += v[i] * xj[i] * xj[i];
            s1 += v[i + 1] * xj[i + 1] * xj[i + 1];
            s2 += v[i + 2] * xj[i + 2] * xj[i + 2];
            s3 += v[i + 3] * xj[i + 3] * xj[i + 3];
          }
          for (; i < m; ++i) s0 += v[i] * xj[i] * xj[i];
          xv[j] = (s0 + s1) + (s2 + s3);
        }
        const double ao = a;
        double dlx;
        // a pass over the active columns, then passes over the nonzero
        // ones until they settle, then the active ones again (as the
        // gaussian path), until a pass over all of them is quiet
        bool full = true;
        for (;;) {
          dlx = 0.0;
          for (int j = 0; j < p; ++j) {
            if (full ? !act[j] : b[j] == 0.0) continue;
            const double *xj = xs + static_cast<long>(j) * m;
            const double gj = dot(xj, r, m);
            const double nb = soft(gj + xv[j] * b[j], l1) / (xv[j] + l2);
            const double dl = nb - b[j];
            if (dl == 0.0) continue;
            for (int i = 0; i < m; ++i) r[i] -= dl * v[i] * xj[i];
            b[j] = nb;
            dlx = std::fmax(dlx, xv[j] * dl * dl);
          }
          double sr = 0.0;
          for (int i = 0; i < m; ++i) sr += r[i];
          const double da = sr / sv;
          a += da;
          for (int i = 0; i < m; ++i) r[i] -= da * v[i];
          dlx = std::fmax(dlx, sv * da * da);
          if (++it > maxit) return CLAB_EDOM;
          if (dlx >= thresh) {
            full = false;
          } else if (full) {
            break;
          } else {
            full = true;
          }
        }
        for (int i = 0; i < m; ++i) eta[i] = a;
        double dout = sv * (a - ao) * (a - ao);
        for (int j = 0; j < p; ++j) {
          if (b[j] != 0.0) {
            const double *xj = xs + static_cast<long>(j) * m;
            for (int i = 0; i < m; ++i) eta[i] += b[j] * xj[i];
          }
          if (act[j]) dout = std::fmax(dout, xv[j] * (b[j] - bo[j]) * (b[j] - bo[j]));
        }
        if (dout < thresh) break;
      }
      // the optimality conditions outside the active set
      for (int i = 0; i < m; ++i) {
        r[i] = w * (yy[i] - 1.0 / (1.0 + std::exp(-eta[i])));
      }
      bool viol = false;
      for (int j = 0; j < p; ++j) {
        if (!(sd[j] > 0.0)) continue;
        if (act[j]) continue;
        ga[j] = std::fabs(dot(xs + static_cast<long>(j) * m, r, m));
        if (ga[j] > l1) {
          act[j] = 1;
          viol = true;
        }
      }
      if (!viol) break;
    }
    double s = a, *bk = beta + static_cast<long>(k) * p;
    for (int j = 0; j < p; ++j) {
      bk[j] = b[j] != 0.0 ? b[j] / sd[j] : 0.0;
      s -= mu[j] * bk[j];
    }
    a0[k] = s;
    *nfit = k + 1;
    // the share of the null deviance explained
    double dev = 0.0;
    for (int i = 0; i < m; ++i) {
      const double pr = 1.0 / (1.0 + std::exp(-eta[i]));
      dev -= yy[i] == 1.0 ? std::log(pr) : std::log1p(-pr);
    }
    const double r2 = 1.0 - dev / m / dev0;
    if (stop && early(k, r2, r20)) break;
    r20 = r2;
  }
  return CLAB_OK;
}

// In-place Cholesky factor (lower) of the q x q matrix a; false when a
// pivot falls below 1e-10 of its diagonal (numerically singular).
bool chol(double *a, int q) {
  for (int j = 0; j < q; ++j) {
    double s = a[j + j * q];
    for (int k = 0; k < j; ++k) s -= a[j + k * q] * a[j + k * q];
    if (!(s > 1e-10 * a[j + j * q])) return false;
    const double l = std::sqrt(s);
    a[j + j * q] = l;
    for (int i = j + 1; i < q; ++i) {
      double t = a[i + j * q];
      for (int k = 0; k < j; ++k) t -= a[i + k * q] * a[j + k * q];
      a[i + j * q] = t / l;
    }
  }
  return true;
}

// Solve L L' x = b in place, L from chol().
void csolve(const double *l, int q, double *b) {
  for (int i = 0; i < q; ++i) {
    double s = b[i];
    for (int k = 0; k < i; ++k) s -= l[i + k * q] * b[k];
    b[i] = s / l[i + i * q];
  }
  for (int i = q - 1; i >= 0; --i) {
    double s = b[i];
    for (int k = i + 1; k < q; ++k) s -= l[k + i * q] * b[k];
    b[i] = s / l[i + i * q];
  }
}

}  // namespace

namespace lrn {

long pwork(int n, int p) {
  return static_cast<long>(p) * p + 6L * p + 4L * n + static_cast<long>(n) * p;
}

int lmax(const Dat &d, int fam, double alpha, double *lam, double *work) {
  const int p = d.p;
  double *mu = work, *sd = mu + p, *t = sd + p;
  double ym, ys;
  if (fam == 1 && !binary(d)) return CLAB_EDOM;
  ystat(d, &ym, &ys);
  if (!(ys > 0.0)) return CLAB_EDOM;
  stdz(d, mu, sd);
  // gaussian: in units of the standardised outcome, then back
  const double sc = fam == 0 ? ys : 1.0;
  for (int i = 0; i < d.m; ++i) t[i] = (d.y[d.idx[i]] - ym) / (sc * d.m);
  double g = 0.0;
  for (int j = 0; j < p; ++j) {
    if (sd[j] > 0.0) g = std::fmax(g, std::fabs(xdot(d, j, mu[j], sd[j], t)));
  }
  *lam = g * sc / std::fmax(alpha, 1e-3);
  return CLAB_OK;
}

int path(const Dat &d, int fam, double alpha, const double *lam, int nl,
         double thresh, int maxit, int stop, int screen, int *nfit,
         double *a0, double *beta, double *work, int *iwork) {
  return fam == 0 ? gauss(d, alpha, lam, nl, thresh, maxit, stop, nfit, a0, beta, work, iwork)
                  : binom(d, alpha, lam, nl, thresh, maxit, stop, screen, nfit, a0, beta, work,
                          iwork);
}

void lpred(const double *x, int n, int p, const int *idx, int m, double a0,
           const double *b, double *eta) {
  for (int i = 0; i < m; ++i) eta[i] = a0;
  for (int j = 0; j < p; ++j) {
    if (b[j] == 0.0) continue;
    const double *xj = x + static_cast<long>(j) * n;
    for (int i = 0; i < m; ++i) eta[i] += b[j] * xj[idx[i]];
  }
}

long gwork(int p) { return 2L * (p + 1) * (p + 1) + 4L * (p + 1); }

// Newton's method on the centred design [1, x - mu]: the intercept is then
// nearly orthogonal to the slopes, which keeps the normal equations well
// conditioned; the covariance is mapped back to the raw intercept at the
// end. work: a (q q), vinv (q q), mu (p, with room for q), u, del, zr (q).
int glm(const Dat &d, int fam, double tol, int maxit, double *coef,
        double *se, double *dev, double *work) {
  const int p = d.p, q = p + 1, m = d.m;
  double *a = work, *vi = a + static_cast<long>(q) * q, *mu = vi + static_cast<long>(q) * q;
  double *u = mu + q, *del = u + q, *zr = del + q;
  if (m <= q || (fam == 1 && !binary(d))) return CLAB_EDOM;
  double ym, ys;
  ystat(d, &ym, &ys);
  if (fam == 1 && !(ys > 0.0)) return CLAB_EDOM;
  for (int j = 0; j < p; ++j) {
    const double *xj = d.x + static_cast<long>(j) * d.n;
    double s = 0.0;
    for (int i = 0; i < m; ++i) s += xj[d.idx[i]];
    mu[j] = s / m;
  }
  // th: the centred coefficients, intercept first, kept in coef
  double *th = coef;
  th[0] = fam == 1 ? std::log(ym / (1.0 - ym)) : ym;
  for (int j = 1; j < q; ++j) th[j] = 0.0;
  double dv = 0.0;
  for (int it = 0;; ++it) {
    if (it >= maxit) return CLAB_EDOM;
    for (long k = 0; k < static_cast<long>(q) * q; ++k) a[k] = 0.0;
    for (int c = 0; c < q; ++c) u[c] = 0.0;
    dv = 0.0;
    // one pass over the rows: the gradient, the information, the deviance;
    // each row is gathered once into zr, then a rank-one update
    for (int i = 0; i < m; ++i) {
      zr[0] = 1.0;
      double eta = th[0];
      for (int c = 1; c < q; ++c) {
        zr[c] = d.x[d.idx[i] + static_cast<long>(c - 1) * d.n] - mu[c - 1];
        eta += th[c] * zr[c];
      }
      const double yi = d.y[d.idx[i]];
      double wi, ri;
      if (fam == 0) {
        wi = 1.0;
        ri = yi - eta;
        dv += ri * ri;
      } else {
        const double pr = 1.0 / (1.0 + std::exp(-eta));
        wi = pr * (1.0 - pr);
        ri = yi - pr;
        dv -= 2.0 * (yi == 1.0 ? std::log(pr) : std::log1p(-pr));
      }
      for (int c = 0; c < q; ++c) {
        u[c] += zr[c] * ri;
        const double wz = wi * zr[c];
        double *ac = a + static_cast<long>(c) * q;
        for (int e = c; e < q; ++e) ac[e] += wz * zr[e];
      }
    }
    for (int c = 0; c < q; ++c) {
      for (int e = c + 1; e < q; ++e) a[c + e * q] = a[e + c * q];
    }
    if (!chol(a, q)) return CLAB_EDOM;
    for (int c = 0; c < q; ++c) del[c] = u[c];
    csolve(a, q, del);
    double big = 0.0;
    for (int c = 0; c < q; ++c) {
      th[c] += del[c];
      big = std::fmax(big, std::fabs(del[c]) / (1.0 + std::fabs(th[c])));
    }
    if (fam == 0 || big < tol) {
      if (fam == 0) {
        // the residual sum of squares at the solution
        dv = 0.0;
        for (int i = 0; i < m; ++i) {
          double eta = th[0];
          for (int c = 1; c < q; ++c) {
            eta += th[c] * (d.x[d.idx[i] + static_cast<long>(c - 1) * d.n] - mu[c - 1]);
          }
          const double e = d.y[d.idx[i]] - eta;
          dv += e * e;
        }
      }
      break;
    }
  }
  // the inverse information, column by column
  for (int c = 0; c < q; ++c) {
    double *vc = vi + static_cast<long>(c) * q;
    for (int e = 0; e < q; ++e) vc[e] = e == c ? 1.0 : 0.0;
    csolve(a, q, vc);
  }
  const double phi = fam == 0 ? dv / (m - q) : 1.0;
  // the raw intercept: th0 - mu' th[1..p], its variance c' V c, c = (1, -mu)
  double v0 = vi[0], a0 = th[0];
  for (int j = 1; j < q; ++j) {
    a0 -= mu[j - 1] * th[j];
    v0 -= 2.0 * mu[j - 1] * vi[j];
    for (int e = 1; e < q; ++e) v0 += mu[j - 1] * mu[e - 1] * vi[j + e * q];
  }
  coef[0] = a0;
  se[0] = std::sqrt(phi * v0);
  for (int j = 1; j < q; ++j) se[j] = std::sqrt(phi * vi[j + j * q]);
  *dev = dv;
  return CLAB_OK;
}

}  // namespace lrn

namespace {

// the loss of one held-out prediction: squared error (gaussian) or the
// binomial deviance with the probability kept in [1e-5, 1 - 1e-5]
inline double loss(int fam, double y, double eta) {
  if (fam == 0) return (y - eta) * (y - eta);
  double pr = 1.0 / (1.0 + std::exp(-eta));
  pr = std::fmin(std::fmax(pr, 1e-5), 1.0 - 1e-5);
  return -2.0 * (y == 1.0 ? std::log(pr) : std::log1p(-pr));
}

}  // namespace

namespace lrn {

long cvwork(int n, int p, int nl, int nf) {
  return pwork(n, p) + n + static_cast<long>(nl) * (p + 1) +
         static_cast<long>(nf) * (nl + 1);
}

int cvnet(const Dat &d, int fam, double alpha, double *lam, int nl,
          double lmr, double thresh, int maxit, int screen, const int *fold, int nf,
          double *a0, double *beta, double *cv, int *best, int *nfit,
          double *work, int *iwork) {
  const int n = d.n, p = d.p;
  int *idx = iwork, *iw = idx + d.m;
  double *pw = work, *eta = pw + pwork(n, p), *fs = eta + n;
  double *fa = fs + static_cast<long>(nf) * nl, *fb = fa + nl;
  double *fn = fb + static_cast<long>(p) * nl;
  if (lmr > 0.0) {
    double lm;
    const int rc = lmax(d, fam, alpha, &lm, pw);
    if (rc != CLAB_OK) return rc;
    for (int k = 0; k < nl; ++k)
      lam[k] = nl == 1 ? lm : lm * std::exp(k * std::log(lmr) / (nl - 1));
  }
  for (int k = 0; k < nl; ++k) {
    if (!(lam[k] >= 0.0 && std::isfinite(lam[k])) || (k > 0 && lam[k] > lam[k - 1]))
      return CLAB_EDOM;
  }
  // a path chosen here stops early as glmnet's does; the folds then fit
  // the lambdas it reached
  int rc = path(d, fam, alpha, lam, nl, thresh, maxit, lmr > 0.0, screen, nfit, a0, beta, pw,
                iw);
  nl = *nfit;
  if (rc != CLAB_OK || fold == nullptr) return rc;
  // K-fold cross-validation over the caller's folds, on the same lambdas;
  // fold[i] belongs to row d.idx[i]
  for (int i = 0; i < d.m; ++i) {
    if (fold[i] < 0 || fold[i] >= nf) return CLAB_EARG;
  }
  for (int f = 0; f < nf; ++f) {
    int m = 0, t = d.m;
    for (int i = 0; i < d.m; ++i) {
      if (fold[i] != f) idx[m++] = d.idx[i]; else idx[--t] = d.idx[i];
    }
    if (t == d.m) return CLAB_EARG;  // an empty fold
    const Dat tr = {d.x, n, p, d.y, idx, m};
    int nk;
    rc = path(tr, fam, alpha, lam, nl, thresh, maxit, 0, screen, &nk, fa, fb, pw, iw);
    if (rc != CLAB_OK) return rc;
    for (int k = 0; k < nl; ++k) {
      lpred(d.x, n, p, idx + t, d.m - t, fa[k], fb + static_cast<long>(k) * p, eta);
      double s = 0.0;
      for (int i = 0; i < d.m - t; ++i) s += loss(fam, d.y[idx[t + i]], eta[i]);
      fs[static_cast<long>(f) * nl + k] = s;
    }
    fn[f] = d.m - t;
  }
  // the mean loss and its standard error over folds, as cv.glmnet: fold
  // means weighted by fold size; the first lambda (the largest) at the
  // minimum, and the largest within one standard error of it
  double *cm = cv, *cs = cv + nl;
  for (int k = 0; k < nl; ++k) {
    double s = 0.0;
    for (int f = 0; f < nf; ++f) s += fs[static_cast<long>(f) * nl + k];
    cm[k] = s / d.m;
    double v = 0.0;
    for (int f = 0; f < nf; ++f) {
      const double e = fs[static_cast<long>(f) * nl + k] / fn[f] - cm[k];
      v += fn[f] * e * e;
    }
    cs[k] = std::sqrt(v / d.m / (nf - 1));
  }
  int im = 0;
  for (int k = 1; k < nl; ++k) {
    if (cm[k] < cm[im]) im = k;
  }
  int i1 = 0;
  while (cm[i1] > cm[im] + cs[im]) ++i1;
  best[0] = im;
  best[1] = i1;
  return CLAB_OK;
}

}  // namespace lrn

extern "C" {

//' @srrstats {G3.0} Convergence of the coordinate descent is a threshold on
//'   the largest weighted squared change of a standardised coefficient in
//'   a pass, never an equality of doubles.
//' @noRd
int clab_enpath(const double *x, int n, int p, const double *y, int fam,
               double alpha, double *lam, int nl, double lmr, double thresh,
               int maxit, int scr, const int *fold, int nf, double *a0,
               double *beta, double *cv, int *best, int *nfit, double *work,
               int *iwork) {
  if (n < 2 || p < 1 || nl < 1 || maxit < 1 || (fam != 0 && fam != 1) ||
      (fold != nullptr && (nf < 2 || nf > n)))
    return CLAB_EARG;
  if (!(alpha >= 0.0 && alpha <= 1.0) || !(thresh > 0.0) || !(lmr < 1.0))
    return CLAB_EDOM;
  for (long k = 0; k < static_cast<long>(n) * p; ++k) {
    if (!std::isfinite(x[k])) return CLAB_EDOM;
  }
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(y[i])) return CLAB_EDOM;
    iwork[i] = i;
  }
  const lrn::Dat d = {x, n, p, y, iwork, n};
  return lrn::cvnet(d, fam, alpha, lam, nl, lmr, thresh, maxit, scr != 0, fold,
                    fold == nullptr ? 0 : nf, a0, beta, cv, best, nfit, work,
                    iwork + n);
}

int clab_lgtreg(const double *x, int n, int p, const double *y, double tol,
                int maxit, double *coef, double *se, double *dev,
                double *work, int *iwork) {
  if (n < 1 || p < 0 || maxit < 1) return CLAB_EARG;
  if (!(tol > 0.0)) return CLAB_EDOM;
  for (long k = 0; k < static_cast<long>(n) * p; ++k) {
    if (!std::isfinite(x[k])) return CLAB_EDOM;
  }
  for (int i = 0; i < n; ++i) iwork[i] = i;
  const lrn::Dat d = {x, n, p, y, iwork, n};
  return lrn::glm(d, 1, tol, maxit, coef, se, dev, work);
}

}  // extern "C"
