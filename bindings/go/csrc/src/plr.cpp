// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The exposure-response estimate of exprspplr(): the partially linear
// regression with cross-fitted ordinary least squares nuisances (the DML2
// estimator of Chernozhukov et al. 2018) on the rows and folds the caller
// drew, and the bootstrap summary of its replicates. Sums run in long
// double, as R's sum(), colMeans() and colSums() run them; the least squares
// fit is the Cholesky factor of the centred and scaled design, or, when that
// is ill conditioned, the Householder QR with limited pivoting of LINPACK's
// dqrdc2 (as R's lm.fit()).
#include <algorithm>
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

namespace {

typedef long double ld;

// The Cholesky factor (upper, a = r'r) of the q x q column-major a.
bool chol(const std::vector<double> &a, int q, std::vector<double> &r) {
  r.assign(static_cast<size_t>(q) * q, 0.0);
  for (int j = 0; j < q; ++j) {
    double s = a[static_cast<size_t>(j) * q + j];
    for (int k = 0; k < j; ++k) s -= r[static_cast<size_t>(j) * q + k] * r[static_cast<size_t>(j) * q + k];
    if (!(s > 0)) return false;
    const double rj = std::sqrt(s);
    r[static_cast<size_t>(j) * q + j] = rj;
    for (int i = j + 1; i < q; ++i) {
      double t = a[static_cast<size_t>(i) * q + j];
      for (int k = 0; k < j; ++k) t -= r[static_cast<size_t>(j) * q + k] * r[static_cast<size_t>(i) * q + k];
      r[static_cast<size_t>(i) * q + j] = t / rj;
    }
  }
  return true;
}

// The reciprocal condition number in the 1-norm of the SPD a with factor r:
// 1 / (|a|_1 |a^-1|_1), the inverse from the factor.
double rcond(const std::vector<double> &a, const std::vector<double> &r, int q) {
  double na = 0;
  for (int j = 0; j < q; ++j) {
    double s = 0;
    for (int i = 0; i < q; ++i) s += std::fabs(a[static_cast<size_t>(j) * q + i]);
    na = std::max(na, s);
  }
  double ni = 0;
  std::vector<double> e(static_cast<size_t>(q));
  for (int j = 0; j < q; ++j) {
    std::fill(e.begin(), e.end(), 0.0);
    e[static_cast<size_t>(j)] = 1;
    for (int i = 0; i < q; ++i) {  // r' w = e
      double t = e[static_cast<size_t>(i)];
      for (int k = 0; k < i; ++k) t -= r[static_cast<size_t>(i) * q + k] * e[static_cast<size_t>(k)];
      e[static_cast<size_t>(i)] = t / r[static_cast<size_t>(i) * q + i];
    }
    for (int i = q - 1; i >= 0; --i) {  // r v = w
      double t = e[static_cast<size_t>(i)];
      for (int k = i + 1; k < q; ++k) t -= r[static_cast<size_t>(k) * q + i] * e[static_cast<size_t>(k)];
      e[static_cast<size_t>(i)] = t / r[static_cast<size_t>(i) * q + i];
    }
    double s = 0;
    for (double v : e) s += std::fabs(v);
    ni = std::max(ni, s);
  }
  return na > 0 && ni > 0 ? 1 / (na * ni) : 0;
}

double nrm2(const double *x, int n) {
  double s = 0;
  for (int i = 0; i < n; ++i) s += x[i] * x[i];
  return std::sqrt(s);
}

// lm.fit() of the two responses yd (nt x 2) on xt (nt x q): LINPACK dqrdc2
// with tolerance 1e-7, aliased coefficients 0. g: q x 2.
void qrfit(std::vector<double> xt, int nt, int q, std::vector<double> yd, std::vector<double> &g) {
  std::vector<double> qra(static_cast<size_t>(q)), w1(static_cast<size_t>(q)), w2(static_cast<size_t>(q));
  std::vector<int> piv(static_cast<size_t>(q));
  auto X = [&](int i, int j) -> double & { return xt[static_cast<size_t>(j) * nt + i]; };
  for (int j = 0; j < q; ++j) {
    piv[static_cast<size_t>(j)] = j;
    qra[static_cast<size_t>(j)] = nrm2(&X(0, j), nt);
    w1[static_cast<size_t>(j)] = qra[static_cast<size_t>(j)];
    w2[static_cast<size_t>(j)] = qra[static_cast<size_t>(j)] == 0 ? 1 : qra[static_cast<size_t>(j)];
  }
  const double tol = 1e-7;
  const int lup = std::min(nt, q);
  int k = q + 1;
  for (int l = 0; l < lup; ++l) {
    while (l + 1 < k && qra[static_cast<size_t>(l)] < w2[static_cast<size_t>(l)] * tol) {
      for (int i = 0; i < nt; ++i) {
        const double t = X(i, l);
        for (int j = l + 1; j < q; ++j) X(i, j - 1) = X(i, j);
        X(i, q - 1) = t;
      }
      const int pi = piv[static_cast<size_t>(l)];
      const double t1 = qra[static_cast<size_t>(l)], t2 = w1[static_cast<size_t>(l)], t3 = w2[static_cast<size_t>(l)];
      for (int j = l + 1; j < q; ++j) {
        piv[static_cast<size_t>(j - 1)] = piv[static_cast<size_t>(j)];
        qra[static_cast<size_t>(j - 1)] = qra[static_cast<size_t>(j)];
        w1[static_cast<size_t>(j - 1)] = w1[static_cast<size_t>(j)];
        w2[static_cast<size_t>(j - 1)] = w2[static_cast<size_t>(j)];
      }
      piv[static_cast<size_t>(q - 1)] = pi;
      qra[static_cast<size_t>(q - 1)] = t1;
      w1[static_cast<size_t>(q - 1)] = t2;
      w2[static_cast<size_t>(q - 1)] = t3;
      --k;
    }
    if (l == nt - 1) break;
    double nx = nrm2(&X(l, l), nt - l);
    if (nx == 0) continue;
    if (X(l, l) != 0) nx = std::copysign(nx, X(l, l));
    for (int i = l; i < nt; ++i) X(i, l) *= 1 / nx;
    X(l, l) = 1 + X(l, l);
    for (int j = l + 1; j < q; ++j) {
      double t = 0;
      for (int i = l; i < nt; ++i) t += X(i, l) * X(i, j);
      t = -t / X(l, l);
      for (int i = l; i < nt; ++i) X(i, j) += t * X(i, l);
      if (qra[static_cast<size_t>(j)] == 0) continue;
      double tt = 1 - (std::fabs(X(l, j)) / qra[static_cast<size_t>(j)]) * (std::fabs(X(l, j)) / qra[static_cast<size_t>(j)]);
      tt = std::max(tt, 0.0);
      if (std::fabs(tt) < 1e-6) {
        qra[static_cast<size_t>(j)] = nrm2(&X(l + 1, j), nt - l - 1);
        w1[static_cast<size_t>(j)] = qra[static_cast<size_t>(j)];
      } else {
        qra[static_cast<size_t>(j)] *= std::sqrt(tt);
      }
    }
    qra[static_cast<size_t>(l)] = X(l, l);
    X(l, l) = -nx;
  }
  const int rank = std::min(k - 1, nt);
  g.assign(static_cast<size_t>(q) * 2, 0.0);
  for (int c = 0; c < 2; ++c) {
    double *y = &yd[static_cast<size_t>(c) * nt];
    // q'y
    for (int l = 0; l < std::min(rank, nt - 1); ++l) {
      if (qra[static_cast<size_t>(l)] == 0) continue;
      const double keep = X(l, l);
      X(l, l) = qra[static_cast<size_t>(l)];
      double t = 0;
      for (int i = l; i < nt; ++i) t += X(i, l) * y[i];
      t = -t / X(l, l);
      for (int i = l; i < nt; ++i) y[i] += t * X(i, l);
      X(l, l) = keep;
    }
    std::vector<double> b(static_cast<size_t>(rank));
    for (int i = rank - 1; i >= 0; --i) {
      double t = y[i];
      for (int j = i + 1; j < rank; ++j) t -= X(i, j) * b[static_cast<size_t>(j)];
      b[static_cast<size_t>(i)] = t / X(i, i);
    }
    for (int i = 0; i < rank; ++i) g[static_cast<size_t>(c) * q + piv[static_cast<size_t>(i)]] = b[static_cast<size_t>(i)];
  }
}

}  // namespace

extern "C" {

int clab_plr(const double *y, const double *d, const double *x, int n, int p, const int *rows, int m,
             const int *fold, int k, double *out) CLAB_TRY {
  out[0] = out[1] = NAN;
  if (n < 1 || p < 0 || m < 1 || k < 2) return CLAB_EARG;
  const int q = p + 1;
  const size_t M = static_cast<size_t>(m);
  // the design (an intercept, then the covariates) and its scaled form
  std::vector<double> X(M * q), Z(M * q), Y(M), D(M);
  for (int i = 0; i < m; ++i) {
    const int r = rows[i];
    if (r < 0 || r >= n) return CLAB_EARG;
    Y[static_cast<size_t>(i)] = y[r];
    D[static_cast<size_t>(i)] = d[r];
    X[static_cast<size_t>(i)] = Z[static_cast<size_t>(i)] = 1;
    for (int j = 0; j < p; ++j) X[(static_cast<size_t>(j) + 1) * M + i] = x[static_cast<size_t>(j) * n + r];
  }
  for (int j = 1; j < q; ++j) {
    const double *c = &X[static_cast<size_t>(j) * M];
    ld s = 0;
    for (int i = 0; i < m; ++i) s += c[i];
    const double ctr = static_cast<double>(s / m);
    ld ss = 0;
    for (int i = 0; i < m; ++i) {
      const double e = c[i] - ctr;
      ss += e * e;
    }
    double scl = std::sqrt(static_cast<double>(ss) / m);
    if (!std::isfinite(scl) || scl <= 0) scl = 1;
    for (int i = 0; i < m; ++i) Z[static_cast<size_t>(j) * M + i] = (c[i] - ctr) / scl;
  }
  std::vector<double> U(M), V(M);
  for (int f = 1; f <= k; ++f) {
    std::vector<int> tr, te;
    for (int i = 0; i < m; ++i) (fold[i] == f ? te : tr).push_back(i);
    if (te.empty()) continue;
    const int nt = static_cast<int>(tr.size());
    std::vector<double> A(static_cast<size_t>(q) * q, 0.0), c(static_cast<size_t>(q) * 2, 0.0), g, r;
    bool fin = true;
    for (int a = 0; a < q; ++a) {
      for (int b = 0; b < q; ++b) {
        double s = 0;
        for (int i : tr) s += Z[static_cast<size_t>(a) * M + i] * Z[static_cast<size_t>(b) * M + i];
        A[static_cast<size_t>(b) * q + a] = s;
        fin = fin && std::isfinite(s);
      }
      double s1 = 0, s2 = 0;
      for (int i : tr) {
        s1 += Z[static_cast<size_t>(a) * M + i] * Y[static_cast<size_t>(i)];
        s2 += Z[static_cast<size_t>(a) * M + i] * D[static_cast<size_t>(i)];
      }
      c[static_cast<size_t>(a)] = s1;
      c[static_cast<size_t>(q + a)] = s2;
    }
    const double *W = Z.data();
    if (fin && chol(A, q, r) && rcond(A, r, q) >= 1e-5) {
      g.assign(static_cast<size_t>(q) * 2, 0.0);
      for (int col = 0; col < 2; ++col) {
        double *b = &g[static_cast<size_t>(col) * q];
        for (int i = 0; i < q; ++i) {
          double t = c[static_cast<size_t>(col) * q + i];
          for (int kk = 0; kk < i; ++kk) t -= r[static_cast<size_t>(i) * q + kk] * b[kk];
          b[i] = t / r[static_cast<size_t>(i) * q + i];
        }
        for (int i = q - 1; i >= 0; --i) {
          double t = b[i];
          for (int kk = i + 1; kk < q; ++kk) t -= r[static_cast<size_t>(kk) * q + i] * b[kk];
          b[i] = t / r[static_cast<size_t>(i) * q + i];
        }
      }
    } else {
      std::vector<double> xt(static_cast<size_t>(nt) * q), yd(static_cast<size_t>(nt) * 2);
      for (int t = 0; t < nt; ++t) {
        const int i = tr[static_cast<size_t>(t)];
        for (int j = 0; j < q; ++j) xt[static_cast<size_t>(j) * nt + t] = X[static_cast<size_t>(j) * M + i];
        yd[static_cast<size_t>(t)] = Y[static_cast<size_t>(i)];
        yd[static_cast<size_t>(nt + t)] = D[static_cast<size_t>(i)];
      }
      qrfit(xt, nt, q, yd, g);
      W = X.data();
    }
    for (int i : te) {
      double f1 = 0, f2 = 0;
      for (int j = 0; j < q; ++j) {
        f1 += W[static_cast<size_t>(j) * M + i] * g[static_cast<size_t>(j)];
        f2 += W[static_cast<size_t>(j) * M + i] * g[static_cast<size_t>(q + j)];
      }
      U[static_cast<size_t>(i)] = Y[static_cast<size_t>(i)] - f1;
      V[static_cast<size_t>(i)] = D[static_cast<size_t>(i)] - f2;
    }
  }
  ld den = 0, sd = 0, num = 0;
  for (int i = 0; i < m; ++i) {
    den += V[static_cast<size_t>(i)] * V[static_cast<size_t>(i)];
    sd += D[static_cast<size_t>(i)];
  }
  ld mdl = sd / m;
  if (std::isfinite(static_cast<double>(mdl))) {
    ld t = 0;
    for (int i = 0; i < m; ++i) t += D[static_cast<size_t>(i)] - mdl;
    mdl += t / m;
  }
  const double md = static_cast<double>(mdl);
  ld ss = 0;
  for (int i = 0; i < m; ++i) {
    const double e = D[static_cast<size_t>(i)] - md;
    ss += e * e;
  }
  const double denom = static_cast<double>(den), scale = std::max(1.0, static_cast<double>(ss));
  if (!std::isfinite(denom) || denom <= 1e-10 * scale) return CLAB_EDOM;
  for (int i = 0; i < m; ++i) num += V[static_cast<size_t>(i)] * U[static_cast<size_t>(i)];
  const double th = static_cast<double>(num) / denom;
  ld s2 = 0;
  for (int i = 0; i < m; ++i) {
    const double v = V[static_cast<size_t>(i)], e = U[static_cast<size_t>(i)] - th * v;
    s2 += v * v * (e * e);
  }
  out[0] = th;
  out[1] = std::sqrt(static_cast<double>(s2)) / denom;
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

void clab_bootci(const double *v, int m, const double *pr, int np, double *out) CLAB_TRY {
  // sd(): the mean refined by a second pass, the squares summed in long double
  out[0] = NAN;
  if (m > 1) {
    ld s = 0;
    for (int i = 0; i < m; ++i) s += v[i];
    ld ml = s / m;
    if (std::isfinite(static_cast<double>(ml))) {
      ld t = 0;
      for (int i = 0; i < m; ++i) t += v[i] - ml;
      ml += t / m;
    }
    const double mu = static_cast<double>(ml);
    ld q = 0;
    for (int i = 0; i < m; ++i) {
      const double e = v[i] - mu;
      q += e * e;
    }
    out[0] = std::sqrt(static_cast<double>(q / (m - 1)));
  }
  // quantile(type = 7): (1 - h) x[lo] + h x[hi]
  std::vector<double> x(v, v + m);
  std::sort(x.begin(), x.end());
  for (int j = 0; j < np; ++j) {
    if (m < 1) {
      out[1 + j] = NAN;
      continue;
    }
    const double idx = 1 + (m - 1) * pr[j];
    const double lo = std::floor(idx), hi = std::ceil(idx);
    const double ql = x[static_cast<size_t>(lo) - 1], qh = x[static_cast<size_t>(hi) - 1];
    double r = ql;
    if (idx > lo && qh != ql) {
      const double h = idx - lo;
      r = (1 - h) * ql + h * qh;
    }
    out[1 + j] = r;
  }
} CLAB_CATCHV

}  // extern "C"
