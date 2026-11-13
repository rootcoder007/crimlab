// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, severity-weighted crime indices (clab_crsev) and composite
// indices of standardised counts (clab_zcomp).
#include <cmath>

#include "crimlab.h"
#include "kutil.h"

int clab_crsev(const double *x, const int *off, const int *key, int n, int no, int nk,
               const double *w, const double *pop, double per, int base, double *cnt, double *sev,
               double *rate, double *srate, double *idx) {
  if (n < 1 || no < 1 || nk < 1 || base < -1 || base >= nk) return CLAB_EARG;
  if (!kutil::pos(per)) return CLAB_EDOM;
  for (int o = 0; o < no; ++o) {
    if (!kutil::nonneg(w[o])) return CLAB_EDOM;
  }
  for (int k = 0; k < nk; ++k) {
    cnt[k] = 0.0;
    sev[k] = 0.0;
    if (pop && !kutil::pos(pop[k])) return CLAB_EDOM;
  }
  for (int i = 0; i < n; ++i) {
    if (off[i] < 0 || off[i] >= no || key[i] < 0 || key[i] >= nk) return CLAB_EARG;
    if (!kutil::nonneg(x[i])) return CLAB_EDOM;
    cnt[key[i]] += x[i];
    sev[key[i]] += x[i] * w[off[i]];
  }
  // the index is the severity rate (or the weighted sum without a
  // population) against the base cell's, times 100
  for (int k = 0; k < nk; ++k) {
    rate[k] = pop ? cnt[k] / pop[k] * per : kutil::nan();
    srate[k] = pop ? sev[k] / pop[k] * per : kutil::nan();
  }
  const double *s = pop ? srate : sev;
  for (int k = 0; k < nk; ++k) {
    idx[k] = base >= 0 && s[base] > 0 ? 100.0 * s[k] / s[base] : kutil::nan();
  }
  return CLAB_OK;
}

int clab_zcomp(const double *x, int nu, int nc, const double *w, double *z, double *comp,
               double *mu, double *sd) {
  if (nu < 2 || nc < 1) return CLAB_EARG;
  for (int c = 0; c < nc; ++c) {
    if (!std::isfinite(w[c])) return CLAB_EDOM;
    const double *xc = x + static_cast<long long>(nu) * c;
    kutil::ksum s;
    for (int u = 0; u < nu; ++u) {
      if (!std::isfinite(xc[u])) return CLAB_EDOM;
      s.add(xc[u]);
    }
    const double m = s.val() / nu;
    kutil::ksum q;
    for (int u = 0; u < nu; ++u) q.add((xc[u] - m) * (xc[u] - m));
    const double d = std::sqrt(q.val() / (nu - 1));
    mu[c] = m;
    // a category that does not vary contributes 0 to every unit
    sd[c] = d > 0 ? d : 1.0;
    for (int u = 0; u < nu; ++u) z[u + static_cast<long long>(nu) * c] = (xc[u] - m) / sd[c];
  }
  for (int u = 0; u < nu; ++u) {
    kutil::ksum s;
    for (int c = 0; c < nc; ++c) s.add(w[c] * z[u + static_cast<long long>(nu) * c]);
    comp[u] = s.val();
  }
  return CLAB_OK;
}
