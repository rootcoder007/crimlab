// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, areas: crosswalks between boundary versions (clab_bndxw)
// and synthetic exposures (clab_synexp).
#include <cmath>

#include "crimlab.h"
#include "kutil.h"

int clab_bndxw(const double *x, const int *u, int n, const int *from, const int *to,
               const double *w, int m, int ns, int nt, double *out, double *wsum, int *best,
               double *bw) {
  if (n < 0 || m < 1 || ns < 1 || nt < 1) return CLAB_EARG;
  for (int k = 0; k < m; ++k) {
    if (from[k] < 0 || from[k] >= ns || to[k] < 0 || to[k] >= nt) return CLAB_EARG;
    if (!kutil::nonneg(w[k])) return CLAB_EDOM;
  }
  for (int i = 0; i < n; ++i) {
    if (u[i] < -1 || u[i] >= ns) return CLAB_EARG;
    if (!std::isfinite(x[i])) return CLAB_EDOM;
  }
  for (int s = 0; s < ns; ++s) {
    wsum[s] = 0;
    best[s] = -1;
    bw[s] = -1;
  }
  for (int t = 0; t < nt; ++t) out[t] = 0;
  // ponytail: the sources' totals in out's place would save ns doubles;
  // wsum doubles as the scratch until the links are read
  for (int i = 0; i < n; ++i) {
    if (u[i] >= 0) wsum[u[i]] += x[i];
  }
  for (int k = 0; k < m; ++k) out[to[k]] += wsum[from[k]] * w[k];
  for (int s = 0; s < ns; ++s) wsum[s] = 0;
  for (int k = 0; k < m; ++k) {
    const int s = from[k];
    wsum[s] += w[k];
    if (w[k] > bw[s]) {
      bw[s] = w[k];
      best[s] = to[k];
    }
  }
  for (int s = 0; s < ns; ++s) {
    if (best[s] < 0) bw[s] = kutil::nan();
  }
  return CLAB_OK;
}

int clab_synexp(const double *coef, const double *ax, int na, int p, const double *pop,
                double *rate, double *expo) {
  if (na < 1 || p < 0) return CLAB_EARG;
  for (int j = 0; j <= p; ++j) {
    if (!std::isfinite(coef[j])) return CLAB_EDOM;
  }
  for (int a = 0; a < na; ++a) {
    if (!kutil::nonneg(pop[a])) return CLAB_EDOM;
    double eta = coef[0];
    for (int j = 0; j < p; ++j) {
      const double v = ax[static_cast<long long>(j) * na + a];
      if (!std::isfinite(v)) return CLAB_EDOM;
      eta += coef[j + 1] * v;
    }
    rate[a] = 1 / (1 + std::exp(-eta));
    expo[a] = rate[a] * pop[a];
  }
  return CLAB_OK;
}
