// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, durations: the Kaplan-Meier estimate by group with
// Greenwood's variance and log-scale intervals, the median, and the
// log-rank test between groups (clab_kmrec). No allocation: the caller
// passes the work arrays.
#include <cmath>
#include <cstddef>

#include "crimlab.h"
#include "kutil.h"
#include "qn.h"

extern "C" {

int clab_kmrec(const double *time, const int *event, const int *grp, int n, int g, double z,
               int *m, int *gout, double *tt, double *nr, double *nd, double *nc, double *s,
               double *se, double *lo, double *hi, double *med, double *lr, double *work,
               int *iwork) {
  if (n < 1 || g < 1 || !(z > 0)) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (!kutil::nonneg(time[i]) || (event[i] != 0 && event[i] != 1) || grp[i] < 0 || grp[i] >= g) {
      return CLAB_EDOM;
    }
  }
  // the records by group, then time
  int *ix = iwork;
  for (int i = 0; i < n; ++i) ix[i] = i;
  kutil::hsort(ix, n, [&](int a, int b) {
    return grp[a] < grp[b] || (grp[a] == grp[b] && (time[a] < time[b] || (time[a] == time[b] && a < b)));
  });
  int k = 0;
  for (int a = 0; a < n;) {
    const int gg = grp[ix[a]];
    int b = a;
    while (b < n && grp[ix[b]] == gg) ++b;
    double risk = b - a, surv = 1, gw = 0;
    med[gg] = kutil::nan();
    bool half = false;
    for (int i = a; i < b;) {
      const double t = time[ix[i]];
      int d = 0, c = 0;
      while (i < b && time[ix[i]] == t) {
        d += event[ix[i]];
        c += 1 - event[ix[i]];
        ++i;
      }
      surv *= 1 - d / risk;
      gw += d < risk ? d / (risk * (risk - d)) : kutil::nan();
      gout[k] = gg;
      tt[k] = t;
      nr[k] = risk;
      nd[k] = d;
      nc[k] = c;
      s[k] = surv;
      se[k] = surv > 0 ? surv * std::sqrt(gw) : kutil::nan();
      lo[k] = surv > 0 ? surv * std::exp(-z * std::sqrt(gw)) : kutil::nan();
      hi[k] = surv > 0 ? std::fmin(1.0, surv * std::exp(z * std::sqrt(gw))) : kutil::nan();
      // the median: the first time the curve reaches 1/2, or, where it
      // sits at exactly 1/2, halfway to the time it next falls
      if (d > 0 && std::isnan(med[gg])) {
        if (half) {
          med[gg] = (med[gg + g] + t) / 2;
        } else if (std::fabs(surv - 0.5) < 1e-12) {
          half = true;
          med[gg + g] = t;
        } else if (surv < 0.5) {
          med[gg] = t;
        }
      }
      risk -= d + c;
      ++k;
    }
    if (half && std::isnan(med[gg])) med[gg] = med[gg + g];
    a = b;
  }
  *m = k;
  // the log-rank test over the distinct times of all groups
  lr[0] = lr[1] = kutil::nan();
  if (g < 2) return CLAB_OK;
  double *atr = work, *dg = atr + g, *oe = dg + g, *v = oe + g;
  const int q = g - 1;
  for (int j = 0; j < g; ++j) atr[j] = oe[j] = 0;
  for (int j = 0; j < q * q; ++j) v[j] = 0;
  for (int i = 0; i < n; ++i) atr[grp[i]] += 1;
  kutil::hsort(ix, n, [&](int a, int b) { return time[a] < time[b] || (time[a] == time[b] && a < b); });
  double tot = n;
  for (int i = 0; i < n;) {
    const double t = time[ix[i]];
    for (int j = 0; j < g; ++j) dg[j] = 0;
    int e = i;
    double d = 0;
    while (e < n && time[ix[e]] == t) {
      dg[grp[ix[e]]] += event[ix[e]];
      d += event[ix[e]];
      ++e;
    }
    if (d > 0) {
      const double f = tot > 1 ? d * (tot - d) / (tot - 1) : 0;
      for (int j = 0; j < g; ++j) oe[j] += dg[j] - d * atr[j] / tot;
      for (int j = 0; j < q; ++j) {
        for (int h = 0; h < q; ++h) {
          v[j + q * h] += f * atr[j] / tot * ((j == h) - atr[h] / tot);
        }
      }
    }
    for (; i < e; ++i) {
      atr[grp[ix[i]]] -= 1;
      tot -= 1;
    }
  }
  if (qn::spdinv(v, q) == 0) {
    double x2 = 0;
    for (int j = 0; j < q; ++j) {
      for (int h = 0; h < q; ++h) x2 += oe[j] * v[j + q * h] * oe[h];
    }
    lr[0] = x2;
    lr[1] = q;
  }
  return CLAB_OK;
}

}  // extern "C"
