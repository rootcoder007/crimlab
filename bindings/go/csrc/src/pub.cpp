// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, published tables: values for banded categories
// (clab_bandx), bounds implied by rounding and suppression (clab_pbbnd)
// and on a change between two published values (clab_pbchg), period over
// period change (clab_yoypx) and the adjustment of many p-values
// (clab_pbadj).
#include <cmath>

#include "crimlab.h"
#include "kutil.h"

int clab_bandx(const double *lo, const double *hi, const int *olo, const int *oup, int n, int rule,
               double cap, double fac, double flr, double *val, double *blo, double *bhi) {
  if (n < 0 || rule < 0 || rule > 3) return CLAB_EARG;
  if (!kutil::pos(fac)) return CLAB_EDOM;
  for (int i = 0; i < n; ++i) {
    double a = lo[i], b = hi[i];
    // an open top band ends at the cap, or at fac times its start; an open
    // bottom band starts at the floor, or at min(0, its end)
    if (oup[i]) b = std::isnan(cap) ? a * fac : cap;
    if (olo[i]) a = std::isnan(flr) ? std::fmin(0.0, b) : flr;
    if (oup[i] && b < a) return CLAB_EDOM;
    if (olo[i] && a > b) return CLAB_EDOM;
    blo[i] = a;
    bhi[i] = b;
    switch (rule) {
      case 0:
        val[i] = (a + b) / 2;
        break;
      case 1:
        val[i] = a;
        break;
      case 2:
        val[i] = b;
        break;
      default:
        // the geometric mean, or the midpoint where it is not defined
        val[i] = a > 0 && std::isfinite(std::sqrt(a * b)) ? std::sqrt(a * b) : (a + b) / 2;
    }
  }
  return CLAB_OK;
}

int clab_pbbnd(const double *val, const int *st, const double *lim, int n, double rnd, int kind,
               double sup, double *lo, double *hi, int *sto) {
  if (n < 0 || kind < 0 || kind > 1) return CLAB_EARG;
  if (!std::isnan(rnd) && !kutil::pos(rnd)) return CLAB_EDOM;
  // a cell rounded to the nearest multiple of rnd lies within rnd / 2 of
  // it; randomly rounded, within rnd - 1
  const double half = kind == 0 ? rnd / 2 : rnd - 1;
  for (int i = 0; i < n; ++i) {
    if (st[i] == 1) {
      const double l = std::isnan(lim[i]) ? sup : lim[i];
      if (!(l > 0) || !std::isfinite(l)) return CLAB_EDOM;
      lo[i] = 0;
      hi[i] = l - 1;
      sto[i] = 2;
    } else if (st[i] == 2 || std::isnan(val[i])) {
      lo[i] = hi[i] = kutil::nan();
      sto[i] = 3;
    } else {
      if (!kutil::nonneg(val[i])) return CLAB_EDOM;
      if (std::isnan(rnd)) {
        lo[i] = hi[i] = val[i];
        sto[i] = 0;
      } else {
        lo[i] = std::fmax(0.0, val[i] - half);
        hi[i] = val[i] + half;
        sto[i] = 1;
      }
    }
  }
  return CLAB_OK;
}

int clab_pbchg(const double *nlo, const double *nhi, const double *plo, const double *phi, int n,
               const double *pop, const double *ppop, double per, double *out) {
  if (n < 0) return CLAB_EARG;
  if (!kutil::pos(per)) return CLAB_EDOM;
  double *clo = out, *chi = out + n, *qlo = chi + n, *qhi = qlo + n, *rlo = qhi + n, *rhi = rlo + n,
         *slo = rhi + n, *shi = slo + n;
  for (int i = 0; i < n; ++i) {
    if ((pop && !kutil::pos(pop[i])) || (ppop && !kutil::pos(ppop[i]))) return CLAB_EDOM;
    // every true pair inside the two boxes: the change runs from the
    // lowest now less the highest before to the reverse
    clo[i] = nlo[i] - phi[i];
    chi[i] = nhi[i] - plo[i];
    qlo[i] = phi[i] > 0 ? 100 * (nlo[i] / phi[i] - 1) : kutil::nan();
    qhi[i] = plo[i] > 0 ? 100 * (nhi[i] / plo[i] - 1) : kutil::nan();
    rlo[i] = pop ? per * nlo[i] / pop[i] : kutil::nan();
    rhi[i] = pop ? per * nhi[i] / pop[i] : kutil::nan();
    slo[i] = ppop ? per * plo[i] / ppop[i] : kutil::nan();
    shi[i] = ppop ? per * phi[i] / ppop[i] : kutil::nan();
  }
  return CLAB_OK;
}

int clab_yoypx(const double *v, const double *prd, const int *grp, int n, double lag, double *prev,
               double *chg, double *pct) {
  if (n < 0) return CLAB_EARG;
  if (!kutil::pos(lag)) return CLAB_EDOM;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(prd[i])) return CLAB_EDOM;
    if (i > 0 && (grp[i] < grp[i - 1] || (grp[i] == grp[i - 1] && prd[i] <= prd[i - 1]))) {
      return CLAB_EARG;
    }
  }
  // rows sorted by group, then period: the earlier period by binary search
  // within the group's run
  int a = 0;
  for (int i = 0; i < n; ++i) {
    if (i > 0 && grp[i] != grp[i - 1]) a = i;
    const double want = prd[i] - lag;
    int l = a, h = i;
    while (l < h) {
      const int m = l + (h - l) / 2;
      if (prd[m] < want)
        l = m + 1;
      else
        h = m;
    }
    prev[i] =
        l < i && std::fabs(prd[l] - want) <= 1e-9 * (1 + std::fabs(want)) ? v[l] : kutil::nan();
    chg[i] = v[i] - prev[i];
    pct[i] = prev[i] != 0 ? 100 * (v[i] - prev[i]) / prev[i] : kutil::nan();
  }
  return CLAB_OK;
}

int clab_pbadj(const double *p, int n, int meth, double *q, int *o) {
  if (n < 0 || meth < 0 || meth > 5) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (!(p[i] >= 0 && p[i] <= 1)) return CLAB_EDOM;
    o[i] = i;
  }
  kutil::hsort(o, n, [&](int a, int b) { return p[a] == p[b] ? a < b : p[a] < p[b]; });
  double by = 0;
  for (int i = n; i >= 1; --i) by += 1.0 / i;
  // 0 none, 1 Bonferroni, 2 Holm, 3 Hochberg, 4 Benjamini-Hochberg, 5
  // Benjamini-Yekutieli, as p.adjust()
  if (meth <= 1) {
    for (int i = 0; i < n; ++i) q[i] = meth == 0 ? p[i] : std::fmin(1.0, n * p[i]);
  } else if (meth == 2) {
    double run = 0;
    for (int i = 0; i < n; ++i) {
      run = std::fmax(run, (n - i) * p[o[i]]);
      q[o[i]] = std::fmin(1.0, run);
    }
  } else {
    double run = 1;
    for (int i = n - 1; i >= 0; --i) {
      const double f = meth == 3   ? n - i
                       : meth == 4 ? static_cast<double>(n) / (i + 1)
                                   : by * n / (i + 1);
      run = std::fmin(run, f * p[o[i]]);
      q[o[i]] = run;
    }
  }
  return CLAB_OK;
}
