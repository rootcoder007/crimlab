// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for clab_kmrec (src/kmr.cpp): the Kaplan-Meier curve,
// Greenwood's sum and the median by hand, and the log-rank statistic of
// two small groups by hand. Part of the self-test (selftest.cpp).
#include <cmath>

#include "crimlab.h"
#include "tchk.h"

namespace {

struct out {
  int m, go[10], iw[10];
  double tt[10], nr[10], nd[10], nc[10], s[10], se[10], lo[10], hi[10], med[6], lr[2], w[20];
};

int km(const double *t, const int *e, const int *g, int n, int ng, out &o) {
  return clab_kmrec(t, e, g, n, ng, 1.959963984540054, &o.m, o.go, o.tt, o.nr, o.nd, o.nc, o.s, o.se,
                    o.lo, o.hi, o.med, o.lr, o.w, o.iw);
}

void kmr(tk &t) {
  static out o;
  // 1, 2, 2+, 3, 4+: S = 0.8, 0.6, 0.3, 0.3
  const double tm[5] = {2, 1, 4, 3, 2};
  const int ev[5] = {1, 1, 0, 1, 0}, g0[5] = {0, 0, 0, 0, 0};
  tk_ieq(t, "km.rc", km(tm, ev, g0, 5, 1, o), CLAB_OK);
  tk_ieq(t, "km.m", o.m, 4);
  tk_eq(t, "km.s1", o.s[0], 0.8, 1e-15);
  tk_eq(t, "km.s2", o.s[1], 0.6, 1e-15);
  tk_eq(t, "km.s3", o.s[2], 0.3, 1e-15);
  tk_eq(t, "km.s4", o.s[3], 0.3, 1e-15);
  tk_eq(t, "km.nc", o.nc[1], 1, 0);
  const double gw = 1.0 / 20 + 1.0 / 12 + 1.0 / 2;
  tk_eq(t, "km.se", o.se[2], 0.3 * std::sqrt(gw), 1e-14);
  tk_eq(t, "km.lo", o.lo[2], 0.3 * std::exp(-1.959963984540054 * std::sqrt(gw)), 1e-14);
  tk_eq(t, "km.med", o.med[0], 3, 0);
  tk_nan(t, "km.lr", o.lr[0]);
  // two groups, one event each: O - E = 1/2, V = 1/4, chi-square 1
  const double t2[2] = {1, 2};
  const int e2[2] = {1, 1}, g2[2] = {0, 1};
  tk_ieq(t, "km.lr.rc", km(t2, e2, g2, 2, 2, o), CLAB_OK);
  tk_eq(t, "km.lr.x2", o.lr[0], 1, 1e-15);
  tk_eq(t, "km.lr.df", o.lr[1], 1, 0);
  tk_eq(t, "km.lr.s", o.s[0], 0, 0);
  tk_nan(t, "km.lr.se", o.se[0]);
  // the second group leaves before any event: the variance is 0
  const double t3[2] = {1, 0.5};
  const int e3[2] = {1, 0};
  tk_ieq(t, "km.sing.rc", km(t3, e3, g2, 2, 2, o), CLAB_OK);
  tk_nan(t, "km.sing", o.lr[0]);
  tk_ieq(t, "km.arg", km(t2, e2, g2, 0, 2, o), CLAB_EARG);
  const int eb[2] = {1, 2};
  tk_ieq(t, "km.dom", km(t2, eb, g2, 2, 2, o), CLAB_EDOM);
}

}  // namespace

void t_kmr(tk &t, int big) {
  (void)big;
  kmr(t);
}
