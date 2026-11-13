// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, segregation stays under the Mandela Rules: the class of
// each placement (clab_mndla) and the shares by group with three
// denominators, placements, people with any qualifying placement and
// people by their cumulative days (clab_mndtab).
#include <cmath>
#include <initializer_list>

#include "crimlab.h"
#include "kutil.h"

int clab_mndla(const double *days, const double *hours, const double *miss, const int *nal,
               int amin, const int *elig, const int *excl, int n, double thr, int *cls) {
  if (n < 0) return CLAB_EARG;
  if (!kutil::nonneg(thr)) return CLAB_EDOM;
  for (int i = 0; i < n; ++i) {
    if (!kutil::nonneg(days[i])) return CLAB_EDOM;
    if (hours) {
      if (!kutil::nonneg(hours[i]) || !(miss[i] >= 0 && miss[i] <= 100)) return CLAB_EDOM;
      // Rules 43 and 44: at most 2 hours out of cell, the 4 hours missed
      // every day, for more than thr days (prolonged) or up to it
      const bool meets = hours[i] <= 2 && miss[i] >= 100;
      cls[i] = meets ? (days[i] > thr ? 2 : 1) : 0;
    } else {
      cls[i] = days[i] > thr || (nal && nal[i] >= amin) ? 1 : 0;
    }
    if ((excl && excl[i]) || (elig && !elig[i])) cls[i] = 0;
  }
  return CLAB_OK;
}

int clab_mndtab(const int *cls, const int *idc, const int *gc, const double *days, const double *w,
                const int *elig, int n, int ng, double thr, double *out, int *o) {
  if (n < 1 || ng < 1) return CLAB_EARG;
  for (int k = 0; k < 7 * (ng + 1); ++k) out[k] = 0.0;
  for (int i = 0; i < n; ++i) {
    if (gc[i] < 0 || gc[i] >= ng || cls[i] < 0 || cls[i] > 2) return CLAB_EARG;
    const double wi = w ? w[i] : 1.0;
    if (!kutil::nonneg(wi) || !kutil::nonneg(days[i])) return CLAB_EDOM;
    for (int g : {gc[i], ng}) {
      double *r = out + 7 * g;
      r[0] += wi;
      r[1] += cls[i] > 0 ? wi : 0.0;
      r[2] += cls[i] == 1 ? wi : 0.0;
      r[3] += cls[i] == 2 ? wi : 0.0;
    }
  }
  if (!idc) return CLAB_OK;
  // people: the rows of each group by id (o), then of all by id (o + n)
  int *a = o, *b = o + n;
  for (int i = 0; i < n; ++i) a[i] = b[i] = i;
  kutil::hsort(a, n, [&](int p, int q) {
    return gc[p] < gc[q] || (gc[p] == gc[q] && (idc[p] < idc[q] || (idc[p] == idc[q] && p < q)));
  });
  kutil::hsort(b, n, [&](int p, int q) { return idc[p] < idc[q] || (idc[p] == idc[q] && p < q); });
  for (int pass = 0; pass < 2; ++pass) {
    const int *s = pass ? b : a;
    for (int i = 0; i < n;) {
      const int g = pass ? ng : gc[s[i]];
      int j = i;
      bool any = false, el = false;
      kutil::ksum cd;
      while (j < n && idc[s[j]] == idc[s[i]] && (pass || gc[s[j]] == gc[s[i]])) {
        any = any || cls[s[j]] > 0;
        el = el || !elig || elig[s[j]];
        cd.add(days[s[j]]);
        ++j;
      }
      double *r = out + 7 * g;
      r[4] += 1;
      r[5] += any;
      r[6] += cd.val() > thr && el;
      i = j;
    }
  }
  return CLAB_OK;
}
