// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, the Oaxaca-Blinder decomposition of a gap in mean outcomes
// between two groups (clab_oaxbd): least squares in each group by
// Householder QR, the twofold decompositions at either group's
// coefficients, the threefold one, the contribution of each covariate, and
// delta-method standard errors (Jann 2008).
#include <cmath>

#include "crimlab.h"
#include "kutil.h"

namespace {

// a' M b for q-vectors and the q x q matrix M (column-major)
double qf(const double *a, const double *m, const double *b, int q) {
  kutil::ksum s;
  for (int j = 0; j < q; ++j) {
    for (int i = 0; i < q; ++i) s.add(a[i] * m[i + q * j] * b[j]);
  }
  return s.val();
}

// Least squares of y on the m x q matrix a (column-major, overwritten by
// its Householder QR): the coefficients b, their covariance v = s^2 (X'X)^-1
// with s^2 = RSS / (m - q), and false when a column is (numerically) a
// combination of the earlier ones. yv: m doubles (overwritten), tau: q,
// ri: q x q.
bool ols(double *a, int m, int q, double *yv, double *tau, double *ri,
         double *b, double *v) {
  for (int k = 0; k < q; ++k) {
    double *c = a + static_cast<long long>(m) * k;
    kutil::ksum o, s;
    for (int i = 0; i < m; ++i) o.add(c[i] * c[i]);
    for (int i = k; i < m; ++i) s.add(c[i] * c[i]);
    const double nk = std::sqrt(s.val());
    // what is left of the column after the earlier ones: none, or rounding
    if (!(nk > 1e-9 * std::sqrt(o.val()))) return false;
    const double al = c[k], be = al > 0.0 ? -nk : nk;
    tau[k] = (be - al) / be;
    const double sc = 1.0 / (al - be);
    for (int i = k + 1; i < m; ++i) c[i] *= sc;
    c[k] = be;
    // H = I - tau v v', v = (1, c[k+1..m)), on the later columns and y
    auto refl = [&](double *z) {
      double d = z[k];
      for (int i = k + 1; i < m; ++i) d += c[i] * z[i];
      d *= tau[k];
      z[k] -= d;
      for (int i = k + 1; i < m; ++i) z[i] -= d * c[i];
    };
    for (int j = k + 1; j < q; ++j) refl(a + static_cast<long long>(m) * j);
    refl(yv);
  }
  // R b = (Q'y)[0..q), and R^-1 column by column
  for (int k = q - 1; k >= 0; --k) {
    double s = yv[k];
    for (int j = k + 1; j < q; ++j) s -= a[k + static_cast<long long>(m) * j] * b[j];
    b[k] = s / a[k + static_cast<long long>(m) * k];
  }
  for (int j = 0; j < q; ++j) {
    for (int k = q - 1; k >= 0; --k) {
      double s = k == j ? 1.0 : 0.0;
      for (int t = k + 1; t < q; ++t) s -= a[k + static_cast<long long>(m) * t] * ri[t + q * j];
      ri[k + q * j] = k > j ? 0.0 : s / a[k + static_cast<long long>(m) * k];
    }
  }
  kutil::ksum rss;
  for (int i = q; i < m; ++i) rss.add(yv[i] * yv[i]);
  const double s2 = rss.val() / (m - q);
  for (int j = 0; j < q; ++j) {
    for (int i = 0; i < q; ++i) {
      double s = 0.0;
      for (int t = 0; t < q; ++t) s += ri[i + q * t] * ri[j + q * t];
      v[i + q * j] = s2 * s;
    }
  }
  return true;
}

}  // namespace

extern "C" {

int clab_oaxbd(const double *y, const double *x, const int *grp, int n, int p,
               int ref, double *xm, double *b, double *bse, double *dec,
               double *var, double *work) {
  if (n < 1 || p < 0 || ref < 0 || ref > 1) return CLAB_EARG;
  const int q = p + 1;
  int nn[2] = {0, 0};
  for (int i = 0; i < n; ++i) {
    if (grp[i] != 0 && grp[i] != 1) return CLAB_EARG;
    if (!std::isfinite(y[i])) return CLAB_EDOM;
    for (int j = 0; j < p; ++j) {
      if (!std::isfinite(x[i + static_cast<long long>(n) * j])) return CLAB_EDOM;
    }
    ++nn[grp[i]];
  }
  // each group needs more cases than coefficients for its residual variance
  if (nn[0] <= q || nn[1] <= q) return CLAB_EDOM;
  // work: the design of one group (n q), its outcomes (n), tau (q), R^-1
  // (q q), then V(b) and V(xbar) of A and of B (4 q q)
  double *a = work, *yv = a + static_cast<long long>(n) * q, *tau = yv + n, *ri = tau + q;
  double *vb[2] = {ri + q * q, ri + 2 * q * q}, *vx[2] = {ri + 3 * q * q, ri + 4 * q * q};
  double ym[2], yvar[2];
  // group A (row 0 of the results) is the one that is not ref
  const int gid[2] = {1 - ref, ref};
  for (int s = 0; s < 2; ++s) {
    const int gs = gid[s], m = nn[gs];
    double *xs = xm + s * q, *bs = b + s * q;
    int r = 0;
    for (int i = 0; i < n; ++i) {
      if (grp[i] != gs) continue;
      a[r] = 1.0;
      for (int j = 0; j < p; ++j) a[r + static_cast<long long>(m) * (j + 1)] = x[i + static_cast<long long>(n) * j];
      yv[r] = y[i];
      ++r;
    }
    // means (one corrective pass) and V(xbar) = S / m from the centred
    // columns, before the QR overwrites them
    for (int j = 0; j < q; ++j) {
      const double *c = a + static_cast<long long>(m) * j;
      kutil::ksum t;
      for (int i = 0; i < m; ++i) t.add(c[i]);
      double mu = t.val() / m;
      kutil::ksum e;
      for (int i = 0; i < m; ++i) e.add(c[i] - mu);
      xs[j] = mu + e.val() / m;
    }
    {
      kutil::ksum t;
      for (int i = 0; i < m; ++i) t.add(yv[i]);
      double mu = t.val() / m;
      kutil::ksum e, v;
      for (int i = 0; i < m; ++i) e.add(yv[i] - mu);
      ym[s] = mu + e.val() / m;
      for (int i = 0; i < m; ++i) v.add((yv[i] - ym[s]) * (yv[i] - ym[s]));
      yvar[s] = v.val() / (m - 1.0);
    }
    for (int j = 0; j < q; ++j) {
      for (int k = 0; k < q; ++k) {
        kutil::ksum t;
        if (j > 0 && k > 0) {
          const double *cj = a + static_cast<long long>(m) * j, *ck = a + static_cast<long long>(m) * k;
          for (int i = 0; i < m; ++i) t.add((cj[i] - xs[j]) * (ck[i] - xs[k]));
        }
        vx[s][j + q * k] = t.val() / (m - 1.0) / m;
      }
    }
    if (!ols(a, m, q, yv, tau, ri, bs, vb[s])) return CLAB_EDOM;
    for (int j = 0; j < q; ++j) bse[s * q + j] = std::sqrt(vb[s][j + q * j]);
  }
  const double *xa = xm, *xb = xm + q, *ba = b, *bb = b + q;
  // D = xbar_A - xbar_B and Delta = b_A - b_B, in the work array's tail
  double *d = tau, *dl = yv;
  for (int j = 0; j < q; ++j) {
    d[j] = xa[j] - xb[j];
    dl[j] = ba[j] - bb[j];
  }
  auto dot = [q](const double *u, const double *w) {
    kutil::ksum s;
    for (int j = 0; j < q; ++j) s.add(u[j] * w[j]);
    return s.val();
  };
  // V(u'w) for independent estimates u and w with covariances U and W:
  // w'Uw + u'Wu (Jann 2008; the term tr(U W) is of a smaller order and
  // left out, as there)
  auto v2 = [q](const double *u, const double *uu, const double *w, const double *ww) {
    return qf(w, uu, w, q) + qf(u, ww, u, q);
  };
  // the sums of the two groups' covariances, in ri (no longer needed)
  double *sx = ri, *sb = a;
  for (int k = 0; k < q * q; ++k) {
    sx[k] = vx[0][k] + vx[1][k];
    sb[k] = vb[0][k] + vb[1][k];
  }
  dec[0] = ym[0];
  dec[1] = ym[1];
  dec[2] = ym[0] - ym[1];
  dec[3] = std::sqrt(yvar[0] / nn[gid[0]] + yvar[1] / nn[gid[1]]);
  // explained at B's coefficients, unexplained at A's means; the same with
  // the roles swapped; the interaction
  dec[4] = dot(d, bb);
  dec[5] = std::sqrt(v2(d, sx, bb, vb[1]));
  dec[6] = dot(xa, dl);
  dec[7] = std::sqrt(v2(xa, vx[0], dl, sb));
  dec[8] = dot(d, ba);
  dec[9] = std::sqrt(v2(d, sx, ba, vb[0]));
  dec[10] = dot(xb, dl);
  dec[11] = std::sqrt(v2(xb, vx[1], dl, sb));
  dec[12] = dot(d, dl);
  dec[13] = std::sqrt(v2(d, sx, dl, sb));
  // by covariate (the constant first): explained at B's coefficients,
  // unexplained at A's means, each with its standard error, and b_A - b_B
  for (int j = 0; j < q; ++j) {
    const long long jj = j + static_cast<long long>(q) * j;
    var[j] = d[j] * bb[j];
    var[q + j] = std::sqrt(d[j] * d[j] * vb[1][jj] + bb[j] * bb[j] * sx[jj]);
    var[2 * q + j] = xa[j] * dl[j];
    var[3 * q + j] = std::sqrt(xa[j] * xa[j] * sb[jj] + dl[j] * dl[j] * vx[0][jj]);
    var[4 * q + j] = dl[j];
  }
  return CLAB_OK;
}

}  // extern "C"
