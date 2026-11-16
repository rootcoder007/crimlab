// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Multiple imputation by predictive mean matching (crimlab.h, mipmm()):
// van Buuren (2018, sec. 3.4 and ch. 4). Each of m chains starts the
// missing values of every target column at draws from its observed ones
// and, maxit times, imputes each target in turn from the others and the
// complete predictors: the least squares fit on the observed rows (ridge
// added to the diagonal of the cross-products, in proportion), a draw of
// the residual scale (rss / chi-square) and of the coefficients about
// their estimate, and for each missing row a donor drawn among the k
// observed rows whose fitted values (with the estimate) are nearest its
// fitted value (with the draw), whose observed value it takes (type 1
// matching). The analyses of the m completed data sets are pooled by
// Rubin's rules with the degrees of freedom of Barnard and Rubin (1999).
// Chi-square and normal draws by inversion from the package's generator.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "crimlab.h"
#include "gla.h"
#include "kutil.h"

namespace {

using std::size_t;
inline size_t z(int i) { return static_cast<size_t>(i); }

double uo(kutil::rng &r) { return (static_cast<double>(r.next() >> 11) + 0.5) * 0x1.0p-53; }

// The least squares fit of v on the q columns of the n x q design a at
// the rows idx (m of them), the cross-products' diagonal raised by ridge
// times itself: b (q) and the Cholesky factor l (q x q). False when it is
// not positive definite.
bool lsq(const double *a, int n, int q, const int *idx, int m, const double *v, double ridge, double *b,
         double *l) {
  for (int u = 0; u < q; ++u) {
    double s = 0;
    for (int i = 0; i < m; ++i) s += a[z(u) * z(n) + z(idx[i])] * v[idx[i]];
    b[u] = s;
    for (int w = 0; w <= u; ++w) {
      double t = 0;
      for (int i = 0; i < m; ++i) t += a[z(u) * z(n) + z(idx[i])] * a[z(w) * z(n) + z(idx[i])];
      l[z(u) + z(w) * z(q)] = l[z(w) + z(u) * z(q)] = t;
    }
  }
  for (int u = 0; u < q; ++u) l[z(u) * (z(q) + 1)] *= 1 + ridge;
  if (!gla::chol(l, q)) return false;
  gla::fsolve(l, q, b, 1);
  gla::bsolve(l, q, b, 1);
  return true;
}

}  // namespace

extern "C" int clab_mipmm(const double *x0, int n, int q, const int *tg, int nt, int m, int k, int maxit,
                          double ridge, unsigned int seed, int fy, const int *fx, int nfx, double level,
                          double *imp, double *pool) CLAB_TRY {
  if (n < 2 || q < 1 || nt < 1 || m < 2 || k < 1 || maxit < 1 || !(ridge >= 0) || nfx < 0 ||
      !(level > 0 && level < 1)) {
    return CLAB_EARG;
  }
  // the missing and observed rows of each target
  std::vector<std::vector<int>> mis(z(nt)), obs(z(nt));
  for (int j = 0; j < nt; ++j) {
    if (tg[j] < 0 || tg[j] >= q) return CLAB_EARG;
    for (int i = 0; i < n; ++i) (std::isnan(x0[z(tg[j]) * z(n) + z(i)]) ? mis : obs)[z(j)].push_back(i);
  }
  for (int c = 0; c < q; ++c) {
    bool t = false;
    for (int j = 0; j < nt; ++j) t = t || tg[j] == c;
    for (int i = 0; !t && i < n; ++i) {
      if (!std::isfinite(x0[z(c) * z(n) + z(i)])) return CLAB_EARG;
    }
  }
  // each target's design: 1 and every other column
  const int qp = q;
  for (int j = 0; j < nt; ++j) {
    if (static_cast<int>(obs[z(j)].size()) <= qp) return CLAB_EDOM;
  }
  size_t off = 0;
  std::vector<size_t> at(z(nt));
  for (int j = 0; j < nt; ++j) {
    at[z(j)] = off;
    off += mis[z(j)].size() * z(m);
  }
  kutil::rng r(seed);
  std::vector<double> x(x0, x0 + z(n) * z(q)), a(z(n) * z(qp)), b(z(qp)), bs(z(qp)), l(z(qp) * z(qp)),
      fo(z(n)), fm(z(n));
  std::vector<int> ord;
  for (int c = 0; c < m; ++c) {
    // the start: draws from the observed values
    for (int j = 0; j < nt; ++j) {
      const std::vector<int> &o = obs[z(j)];
      double *xc = &x[z(tg[j]) * z(n)];
      for (int i : mis[z(j)]) xc[i] = xc[o[z(r.below(static_cast<int>(o.size())))]];
    }
    for (int it = 0; it < maxit; ++it) {
      for (int j = 0; j < nt; ++j) {
        if (kutil::stop(static_cast<double>(n) * qp * qp)) return CLAB_EINTR;
        const std::vector<int> &o = obs[z(j)], &ms = mis[z(j)];
        if (ms.empty()) continue;
        double *yc = &x[z(tg[j]) * z(n)];
        for (int i = 0; i < n; ++i) a[z(i)] = 1;
        for (int u = 0, w = 1; u < q; ++u) {
          if (u == tg[j]) continue;
          std::copy(&x[z(u) * z(n)], &x[z(u) * z(n)] + n, &a[z(w++) * z(n)]);
        }
        const int no = static_cast<int>(o.size());
        if (!lsq(a.data(), n, qp, o.data(), no, yc, ridge, b.data(), l.data())) return CLAB_EDOM;
        double rss = 0;
        for (int i = 0; i < n; ++i) {
          double f = 0;
          for (int u = 0; u < qp; ++u) f += a[z(u) * z(n) + z(i)] * b[z(u)];
          fo[z(i)] = f;
        }
        for (int i : o) rss += (yc[i] - fo[z(i)]) * (yc[i] - fo[z(i)]);
        const double sig = std::sqrt(rss / clab_qchisq(uo(r), no - qp, 1));
        // the coefficients drawn: b + sig L'^-1 e, L L' the cross-products
        for (int u = 0; u < qp; ++u) bs[z(u)] = clab_qnorm(uo(r));
        gla::bsolve(l.data(), qp, bs.data(), 1);
        for (int u = 0; u < qp; ++u) bs[z(u)] = b[z(u)] + sig * bs[z(u)];
        for (int i : ms) {
          double f = 0;
          for (int u = 0; u < qp; ++u) f += a[z(u) * z(n) + z(i)] * bs[z(u)];
          fm[z(i)] = f;
        }
        // the donors: the k observed rows nearest in fitted value, found
        // from the observed rows sorted by it (ties by row)
        ord = o;
        std::sort(ord.begin(), ord.end(), [&](int p1, int p2) {
          return fo[z(p1)] < fo[z(p2)] || (fo[z(p1)] == fo[z(p2)] && p1 < p2);
        });
        const int kk = std::min(k, no);
        for (int i : ms) {
          const double t = fm[z(i)];
          int hi = static_cast<int>(std::lower_bound(ord.begin(), ord.end(), t,
                                                     [&](int p1, double v) { return fo[z(p1)] < v; }) -
                                    ord.begin());
          int lo = hi - 1, got = 0, pick = r.below(kk), don = -1;
          while (got <= pick) {
            const bool left = hi >= no || (lo >= 0 && t - fo[z(ord[z(lo)])] <= fo[z(ord[z(hi)])] - t);
            don = left ? ord[z(lo--)] : ord[z(hi++)];
            ++got;
          }
          yc[i] = x0[z(tg[j]) * z(n) + z(don)];
        }
      }
    }
    for (int j = 0; j < nt; ++j) {
      const std::vector<int> &ms = mis[z(j)];
      for (size_t u = 0; u < ms.size(); ++u) {
        imp[at[z(j)] + z(c) * ms.size() + u] = x[z(tg[j]) * z(n) + z(ms[u])];
      }
    }
  }
  // Rubin's rules over the m completed data sets: the least squares
  // coefficients of column fy on [1, fx] with their classical variances,
  // or (fy < 0) the mean of each target with its variance s^2 / n
  const int np = fy >= 0 ? nfx + 1 : nt;
  const double dfc = fy >= 0 ? n - (nfx + 1) : n - 1;
  if (!(dfc > 0)) return CLAB_EDOM;
  std::vector<double> qm(z(np) * z(m)), um(z(np) * z(m)), d(z(n) * z(nfx + 1)), lv(z(nfx + 1) * z(nfx + 1)),
      bb(z(nfx + 1)), e(z(nfx + 1));
  std::vector<int> all(z(n));
  for (int i = 0; i < n; ++i) all[z(i)] = i;
  for (int c = 0; c < m; ++c) {
    // the completed data set c
    std::copy(x0, x0 + z(n) * z(q), x.begin());
    for (int j = 0; j < nt; ++j) {
      const std::vector<int> &ms = mis[z(j)];
      for (size_t u = 0; u < ms.size(); ++u) {
        x[z(tg[j]) * z(n) + z(ms[u])] = imp[at[z(j)] + z(c) * ms.size() + u];
      }
    }
    if (fy >= 0) {
      for (int i = 0; i < n; ++i) d[z(i)] = 1;
      for (int u = 0; u < nfx; ++u) std::copy(&x[z(fx[u]) * z(n)], &x[z(fx[u]) * z(n)] + n, &d[z(u + 1) * z(n)]);
      const double *yv = &x[z(fy) * z(n)];
      if (!lsq(d.data(), n, nfx + 1, all.data(), n, yv, 0, bb.data(), lv.data())) return CLAB_EDOM;
      double rss = 0;
      for (int i = 0; i < n; ++i) {
        double f = 0;
        for (int u = 0; u <= nfx; ++u) f += d[z(u) * z(n) + z(i)] * bb[z(u)];
        rss += (yv[i] - f) * (yv[i] - f);
      }
      const double s2 = rss / dfc;
      for (int u = 0; u <= nfx; ++u) {
        // the u-th diagonal element of (D'D)^-1: |L^-1 e_u|^2
        std::fill(e.begin(), e.end(), 0.0);
        e[z(u)] = 1;
        gla::fsolve(lv.data(), nfx + 1, e.data(), 1);
        double s = 0;
        for (double v : e) s += v * v;
        qm[z(u) * z(m) + z(c)] = bb[z(u)];
        um[z(u) * z(m) + z(c)] = s2 * s;
      }
    } else {
      for (int j = 0; j < nt; ++j) {
        const double *v = &x[z(tg[j]) * z(n)];
        double s = 0, ss = 0;
        for (int i = 0; i < n; ++i) s += v[i];
        s /= n;
        for (int i = 0; i < n; ++i) ss += (v[i] - s) * (v[i] - s);
        qm[z(j) * z(m) + z(c)] = s;
        um[z(j) * z(m) + z(c)] = ss / (n - 1) / n;
      }
    }
  }
  const double tq = level / 2 + 0.5;
  for (int u = 0; u < np; ++u) {
    double qb = 0, ub = 0, bv = 0;
    for (int c = 0; c < m; ++c) {
      qb += qm[z(u) * z(m) + z(c)];
      ub += um[z(u) * z(m) + z(c)];
    }
    qb /= m;
    ub /= m;
    for (int c = 0; c < m; ++c) bv += (qm[z(u) * z(m) + z(c)] - qb) * (qm[z(u) * z(m) + z(c)] - qb);
    bv /= m - 1;
    const double t = ub + (1 + 1.0 / m) * bv, lam = (1 + 1.0 / m) * bv / t, riv = (1 + 1.0 / m) * bv / ub;
    const double tmp = (1 - lam) * (1 + dfc) * dfc, df = (m - 1) * tmp / ((dfc + 3) * (m - 1) + lam * lam * tmp);
    const double se = std::sqrt(t), c = clab_qt(tq, df, 1);
    pool[u] = qb;
    pool[np + u] = se;
    pool[2 * np + u] = df;
    pool[3 * np + u] = qb - c * se;
    pool[4 * np + u] = qb + c * se;
    pool[5 * np + u] = 2 * clab_pt(-std::fabs(qb / se), df, 1);
    pool[6 * np + u] = riv;
    pool[7 * np + u] = (riv + 2 / (df + 3)) / (riv + 1);
  }
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)
