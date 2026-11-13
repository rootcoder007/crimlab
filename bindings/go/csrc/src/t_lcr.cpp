// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for the life-course and pipeline kernels (src/lcr.cpp)
// against values worked by hand. Part of the self-test (selftest.cpp).
#include <cmath>
#include <limits>

#include "crimlab.h"
#include "tchk.h"

namespace {

const double nan = std::numeric_limits<double>::quiet_NaN();

void mix(tk &t) {
  // two types, one offending only at the first age and one only at the
  // second, half each
  const double cv[4] = {1, 0, 0, 1}, w[2] = {0.5, 0.5};
  double agg[2], cov[4];
  tk_ieq(t, "agemix.rc", clab_agemix(cv, 2, 2, w, agg, cov), CLAB_OK);
  tk_eq(t, "agemix.agg", agg[1], 0.5, 0);
  tk_eq(t, "agemix.var", cov[0], 0.25, 0);
  tk_eq(t, "agemix.cov", cov[2], -0.25, 0);
  tk_eq(t, "agemix.sym", cov[1], cov[2], 0);
  // one type: no covariance
  const double w1[2] = {1, 0};
  clab_agemix(cv, 2, 2, w1, agg, cov);
  tk_eq(t, "agemix.one", cov[2], 0, 0);
  const double w9[2] = {0.5, 0.4}, wn[2] = {1.5, -0.5}, cn[4] = {nan, 0, 0, 1};
  tk_ieq(t, "agemix.n", clab_agemix(cv, 0, 2, w, agg, cov), CLAB_EARG);
  tk_ieq(t, "agemix.sum", clab_agemix(cv, 2, 2, w9, agg, cov), CLAB_EDOM);
  tk_ieq(t, "agemix.neg", clab_agemix(cv, 2, 2, wn, agg, cov), CLAB_EDOM);
  tk_ieq(t, "agemix.nan", clab_agemix(cn, 2, 2, w, agg, cov), CLAB_EDOM);
  // the panel covariance and its permutation test
  const double ya[4] = {1, 2, 3, 4}, yb[4] = {2, 4, 6, 8};
  double o[2], wk[4], o2[2];
  tk_ieq(t, "agepan.rc", clab_agepan(ya, yb, 4, 0, 1, o, wk), CLAB_OK);
  tk_eq(t, "agepan.cov", o[0], 2.5, 1e-15);
  tk_eq(t, "agepan.p0", o[1], 1, 0);
  clab_agepan(ya, yb, 4, 200, 5, o, wk);
  clab_agepan(ya, yb, 4, 200, 5, o2, wk);
  tk_true(t, "agepan.p", o[1] > 0 && o[1] < 1 && o[1] == o2[1]);
  const double yn[4] = {1, nan, 3, 4};
  tk_ieq(t, "agepan.n", clab_agepan(ya, yb, 2, 0, 1, o, wk), CLAB_EARG);
  tk_ieq(t, "agepan.seed", clab_agepan(ya, yb, 4, 0, -1, o, wk), CLAB_EDOM);
  tk_ieq(t, "agepan.nan", clab_agepan(yn, yb, 4, 0, 1, o, wk), CLAB_EDOM);
}

void risk(tk &t) {
  const double age[2] = {18, 19}, e[2] = {10, 9}, pop[2] = {100, 100};
  double h[2], nr[2], cr[2], v[2];
  tk_ieq(t, "cumrsk.rc", clab_cumrsk(age, e, pop, nullptr, nullptr, 2, h, nr, cr, v), CLAB_OK);
  tk_eq(t, "cumrsk.nr", nr[1], 90, 1e-15);
  tk_eq(t, "cumrsk.cr", cr[1], 0.19, 1e-15);
  tk_eq(t, "cumrsk.v", v[1], 0.1 / 90 + 0.1 / 81, 1e-15);
  const double ar[2] = {50, 40}, dt[2] = {0.01, 0.01};
  clab_cumrsk(age, e, pop, ar, nullptr, 2, h, nr, cr, v);
  tk_eq(t, "cumrsk.atrisk", h[1], 0.225, 1e-15);
  clab_cumrsk(age, e, pop, nullptr, dt, 2, h, nr, cr, v);
  tk_eq(t, "cumrsk.dth", cr[1], 0.1 + 0.89 * 0.1, 1e-15);
  tk_nan(t, "cumrsk.nov", v[0]);
  // every one at risk has an event: survival 0, an infinite variance
  const double ea[1] = {100};
  clab_cumrsk(age, ea, pop, nullptr, nullptr, 1, h, nr, cr, v);
  tk_true(t, "cumrsk.all", cr[0] == 1 && std::isinf(v[0]));
  const double e2[1] = {200}, a2[2] = {19, 18}, p0[1] = {0}, d1[2] = {1, 0}, d9[2] = {0.95, 0};
  tk_ieq(t, "cumrsk.n", clab_cumrsk(age, e, pop, nullptr, nullptr, 0, h, nr, cr, v), CLAB_EARG);
  tk_ieq(t, "cumrsk.h", clab_cumrsk(age, e2, pop, nullptr, nullptr, 1, h, nr, cr, v), CLAB_EDOM);
  tk_ieq(t, "cumrsk.age", clab_cumrsk(a2, e, pop, nullptr, nullptr, 2, h, nr, cr, v), CLAB_EDOM);
  tk_ieq(t, "cumrsk.pop", clab_cumrsk(age, e, p0, nullptr, nullptr, 1, h, nr, cr, v), CLAB_EDOM);
  tk_ieq(t, "cumrsk.ar", clab_cumrsk(age, e, pop, p0, nullptr, 1, h, nr, cr, v), CLAB_EDOM);
  tk_ieq(t, "cumrsk.dth1", clab_cumrsk(age, e, pop, nullptr, d1, 2, h, nr, cr, v), CLAB_EDOM);
  tk_ieq(t, "cumrsk.cohort", clab_cumrsk(age, e, pop, nullptr, d9, 2, h, nr, cr, v), CLAB_EDOM);
}

void pipe(tk &t) {
  const double num[4] = {5, 7, 6, 9}, den[4] = {20, 22, 19, 25};
  double on[4], od[4], r[4];
  tk_ieq(t, "jpipe.rc", clab_jpipe(num, den, 4, 2, 1, on, od, r), CLAB_OK);
  tk_nan(t, "jpipe.first", on[0]);
  tk_eq(t, "jpipe.num", on[3], 15, 0);
  tk_nan(t, "jpipe.den1", od[1]);
  tk_eq(t, "jpipe.den", od[3], 41, 0);
  tk_eq(t, "jpipe.rate", r[2], 13.0 / 42, 1e-15);
  tk_nan(t, "jpipe.nanrate", r[1]);
  const double d0[4] = {0, 0, 0, 0}, nn[1] = {-1};
  clab_jpipe(num, d0, 4, 1, 0, on, od, r);
  tk_nan(t, "jpipe.zero", r[0]);
  tk_ieq(t, "jpipe.win", clab_jpipe(num, den, 4, 0, 0, on, od, r), CLAB_EARG);
  tk_ieq(t, "jpipe.lag", clab_jpipe(num, den, 4, 1, -1, on, od, r), CLAB_EARG);
  tk_ieq(t, "jpipe.neg", clab_jpipe(nn, den, 1, 1, 0, on, od, r), CLAB_EDOM);
}

}  // namespace

void t_lcr(tk &t, int big) {
  (void)big;
  mix(t);
  risk(t);
  pipe(t);
}
