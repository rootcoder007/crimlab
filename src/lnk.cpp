// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, crime linkage: the distance, the time and the similarity of
// behaviour (Jaccard) of every pair of crimes, a logistic regression of
// whether two solved crimes belong to one series on them, and the area
// under the ROC curve of each (Bennell and Canter 2002).
#include <algorithm>
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

namespace {

// The Mann-Whitney area: the chance that a linked pair scores above an
// unlinked one, ties counting a half, over the pairs with lk >= 0.
double auc(const double *s, const int *lk, std::size_t m) {
  std::vector<std::size_t> o;
  for (std::size_t k = 0; k < m; ++k) {
    if (lk[k] >= 0) o.push_back(k);
  }
  std::sort(o.begin(), o.end(), [&](std::size_t a, std::size_t b) { return s[a] < s[b]; });
  double r1 = 0, n1 = 0, n0 = 0;
  for (std::size_t q = 0; q < o.size();) {
    std::size_t e = q;
    while (e < o.size() && s[o[e]] == s[o[q]]) ++e;
    const double rk = 0.5 * (static_cast<double>(q) + 1 + static_cast<double>(e));
    for (std::size_t u = q; u < e; ++u) {
      if (lk[o[u]]) {
        r1 += rk;
        n1 += 1;
      } else {
        n0 += 1;
      }
    }
    q = e;
  }
  return n1 > 0 && n0 > 0 ? (r1 - n1 * (n1 + 1) / 2) / (n1 * n0) : kutil::nan();
}

}  // namespace

int clab_lnkpr(const double *px, const double *py, const double *t, const int *mo, int q, const int *ser, int n,
               int *pi, int *pj, double *dist, double *dt, double *jac, int *lk, double *prob, double *coef,
               double *se, double *aucs, double *dev) CLAB_TRY {
  if (n < 2 || n > 20000 || q < 0) return CLAB_EARG;
  const std::size_t m = static_cast<std::size_t>(n) * (static_cast<std::size_t>(n) - 1) / 2;
  std::size_t k = 0;
  int n1 = 0, n0 = 0;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(px[i]) || !std::isfinite(py[i]) || !std::isfinite(t[i])) return CLAB_EDOM;
  }
  for (int i = 0; i < n; ++i) {
    if (kutil::stop(static_cast<double>(n) * (q + 4))) return CLAB_EINTR;
    for (int j = i + 1; j < n; ++j, ++k) {
      pi[k] = i;
      pj[k] = j;
      dist[k] = std::hypot(px[i] - px[j], py[i] - py[j]);
      dt[k] = std::fabs(t[i] - t[j]);
      if (q > 0) {
        int a = 0, u = 0;
        for (int b = 0; b < q; ++b) {
          const int x1 = mo[static_cast<std::size_t>(b) * n + i] != 0, x2 = mo[static_cast<std::size_t>(b) * n + j] != 0;
          a += x1 && x2;
          u += x1 || x2;
        }
        jac[k] = u > 0 ? static_cast<double>(a) / u : 0.0;
      } else {
        jac[k] = kutil::nan();
      }
      lk[k] = ser[i] < 0 || ser[j] < 0 ? -1 : ser[i] == ser[j];
      if (lk[k] == 1) ++n1;
      if (lk[k] == 0) ++n0;
    }
  }
  if (n1 < 1 || n0 < 1) return CLAB_EDOM;
  // the model on the pairs of solved crimes
  const int p = q > 0 ? 3 : 2, nt = n1 + n0;
  const std::size_t un = static_cast<std::size_t>(nt);
  std::vector<double> X(un * p), y(un);
  std::size_t r = 0;
  for (std::size_t s = 0; s < m; ++s) {
    if (lk[s] < 0) continue;
    X[r] = dist[s];
    X[un + r] = dt[s];
    if (p == 3) X[2 * un + r] = jac[s];
    y[r] = lk[s];
    ++r;
  }
  std::vector<double> w(2 * static_cast<std::size_t>(p + 1) * (p + 1) + 4 * static_cast<std::size_t>(p + 1));
  std::vector<int> iw(un);
  const int rc = clab_lgtreg(X.data(), nt, p, y.data(), 1e-10, 100, coef, se, dev, w.data(), iw.data());
  if (rc != CLAB_OK) return rc;
  std::vector<double> eta(m), nd(m), nt2(m);
  for (std::size_t s = 0; s < m; ++s) {
    double e = coef[0] + coef[1] * dist[s] + coef[2] * dt[s];
    if (p == 3) e += coef[3] * jac[s];
    eta[s] = e;
    prob[s] = 1.0 / (1.0 + std::exp(-e));
    nd[s] = -dist[s];
    nt2[s] = -dt[s];
  }
  // aucs: distance (smaller is linked), time (smaller), Jaccard (larger),
  // the model's linear predictor
  aucs[0] = auc(nd.data(), lk, m);
  aucs[1] = auc(nt2.data(), lk, m);
  aucs[2] = p == 3 ? auc(jac, lk, m) : kutil::nan();
  aucs[3] = auc(eta.data(), lk, m);
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)
