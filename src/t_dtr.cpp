// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for the deterrence kernels (src/dtr.cpp) against
// optima, ranks and steps worked by hand. Part of the self-test
// (selftest.cpp).
#include <cmath>

#include "crimlab.h"
#include "tchk.h"

namespace {

void detrsp(tk &t) {
  const double x[3] = {0, 1, 2}, b[3] = {0, 2, 3}, s[3] = {0, 1, 2}, p[4] = {0, 1, 2, 3};
  double xo[4], v[4];
  tk_ieq(t, "detrsp.rc", clab_detrsp(x, b, s, 3, p, 4, xo, v), CLAB_OK);
  tk_eq(t, "detrsp.x0", xo[0], 2, 0);
  tk_eq(t, "detrsp.v0", v[0], 3, 0);
  tk_eq(t, "detrsp.x1", xo[1], 2, 0);
  tk_eq(t, "detrsp.x2", xo[2], 1, 0);
  tk_eq(t, "detrsp.x3", xo[3], 0, 0);
  tk_eq(t, "detrsp.v3", v[3], 0, 0);
  const double xb[3] = {0, 0, 2}, sb[3] = {0, 2, 1}, bn[3] = {0, NAN, 1}, pn[1] = {-1};
  tk_ieq(t, "detrsp.n", clab_detrsp(x, b, s, 0, p, 1, xo, v), CLAB_EARG);
  tk_ieq(t, "detrsp.np", clab_detrsp(x, b, s, 3, p, 0, xo, v), CLAB_EARG);
  tk_ieq(t, "detrsp.xinc", clab_detrsp(xb, b, s, 3, p, 1, xo, v), CLAB_EDOM);
  tk_ieq(t, "detrsp.sinc", clab_detrsp(x, b, sb, 3, p, 1, xo, v), CLAB_EDOM);
  tk_ieq(t, "detrsp.nan", clab_detrsp(x, bn, s, 3, p, 1, xo, v), CLAB_EDOM);
  tk_ieq(t, "detrsp.p", clab_detrsp(x, b, s, 3, pn, 1, xo, v), CLAB_EDOM);
}

void dtrid(tk &t) {
  const double p[4] = {0, 1, 0, 0}, s[4] = {0, 0, 1, 0}, c[4] = {0, 0, 0, 1};
  const double y[4] = {1, 3, 4, 5};
  int rk = 0, id[3], cs[3];
  double b[4], w[32];
  tk_ieq(t, "dtrid.rc", clab_dtrid(p, s, c, 4, 0, y, &rk, id, cs, b, w), CLAB_OK);
  tk_ieq(t, "dtrid.rank", rk, 4);
  tk_ieq(t, "dtrid.id", id[0] && id[1] && id[2], 1);
  tk_eq(t, "dtrid.b0", b[0], 1, 1e-14);
  tk_eq(t, "dtrid.b1", b[1], 2, 1e-14);
  tk_eq(t, "dtrid.b3", b[3], 4, 1e-14);
  // certainty varied alone: rank 2, only certainty identified
  const double p2[3] = {0.1, 0.3, 0.5}, s2[3] = {2, 2, 2}, c2[3] = {30, 30, 30};
  clab_dtrid(p2, s2, c2, 3, 0, nullptr, &rk, id, cs, b, w);
  tk_ieq(t, "dtrid.rank2", rk, 2);
  tk_ieq(t, "dtrid.idp", id[0], 1);
  tk_ieq(t, "dtrid.ids", id[1], 0);
  tk_ieq(t, "dtrid.csts", cs[1], 1);
  tk_ieq(t, "dtrid.cstp", cs[0], 0);
  // the log design: log p and log s
  const double pl[4] = {1, std::exp(1.0), 1, 1}, sl[4] = {1, 1, std::exp(1.0), 1};
  tk_ieq(t, "dtrid.lg", clab_dtrid(pl, sl, c, 4, 1, y, &rk, id, cs, b, w), CLAB_OK);
  tk_eq(t, "dtrid.lgb", b[2], 3, 1e-13);
  // a zero column
  const double z[4] = {0, 0, 0, 0};
  clab_dtrid(z, s, c, 4, 0, nullptr, &rk, id, cs, b, w);
  tk_ieq(t, "dtrid.zero", rk, 3);
  const double bad[4] = {0, NAN, 0, 0};
  tk_ieq(t, "dtrid.n", clab_dtrid(p, s, c, 0, 0, y, &rk, id, cs, b, w), CLAB_EARG);
  tk_ieq(t, "dtrid.rk", clab_dtrid(p2, s2, c2, 3, 0, y, &rk, id, cs, b, w), CLAB_EDOM);
  tk_ieq(t, "dtrid.lgp", clab_dtrid(p, s, c, 4, 1, nullptr, &rk, id, cs, b, w), CLAB_EDOM);
  tk_ieq(t, "dtrid.lgs", clab_dtrid(pl, s, c, 4, 1, nullptr, &rk, id, cs, b, w), CLAB_EDOM);
  tk_ieq(t, "dtrid.nan", clab_dtrid(bad, s, c, 4, 0, nullptr, &rk, id, cs, b, w), CLAB_EDOM);
  tk_ieq(t, "dtrid.nans", clab_dtrid(p, bad, c, 4, 0, nullptr, &rk, id, cs, b, w), CLAB_EDOM);
  tk_ieq(t, "dtrid.nanc", clab_dtrid(p, s, bad, 4, 0, nullptr, &rk, id, cs, b, w), CLAB_EDOM);
  tk_ieq(t, "dtrid.nany", clab_dtrid(p, s, c, 4, 0, bad, &rk, id, cs, b, w), CLAB_EDOM);
}

void inspg(tk &t) {
  const double x0[3] = {0.5, 0.3, 0.2};
  double tr[603], v[201], f[3], dv = 1;
  tk_ieq(t, "inspg.rc", clab_inspg(1, x0, 200, 0.01, tr, v, f, &dv), CLAB_OK);
  tk_eq(t, "inspg.dv", dv, 0, 1e-17);
  tk_eq(t, "inspg.f0", f[0], 0.5 * 0.1, 1e-15);
  tk_eq(t, "inspg.start", tr[201], 0.3, 0);
  double drift = 0;
  for (int k = 0; k < 201; ++k) drift = std::fmax(drift, std::fabs(v[k] - 0.03));
  tk_lt(t, "inspg.drift", drift, 1e-10);
  tk_eq(t, "inspg.sum", tr[200] + tr[401] + tr[602], 1, 1e-12);
  const double xb[3] = {0.5, 0.5, 0.5}, xn[3] = {-0.5, 1, 0.5};
  tk_ieq(t, "inspg.steps", clab_inspg(1, x0, -1, 0.01, tr, v, f, &dv), CLAB_EARG);
  tk_ieq(t, "inspg.a", clab_inspg(0, x0, 1, 0.01, tr, v, f, &dv), CLAB_EDOM);
  tk_ieq(t, "inspg.ainf", clab_inspg(INFINITY, x0, 1, 0.01, tr, v, f, &dv), CLAB_EDOM);
  tk_ieq(t, "inspg.dt", clab_inspg(1, x0, 1, 0, tr, v, f, &dv), CLAB_EDOM);
  tk_ieq(t, "inspg.sum1", clab_inspg(1, xb, 1, 0.01, tr, v, f, &dv), CLAB_EDOM);
  tk_ieq(t, "inspg.neg", clab_inspg(1, xn, 1, 0.01, tr, v, f, &dv), CLAB_EDOM);
}

void insphs(tk &t) {
  // one step from the centroid, T 1.5, g 0.5, rate 1: (19, 19, 16) / 54
  const double T[1] = {1.5}, g[1] = {0.5}, x0[3] = {1, 1, 1};
  double o[3];
  tk_ieq(t, "insphs.rc", clab_insphs(T, 1, g, 1, 1, 1, x0, o), CLAB_OK);
  tk_eq(t, "insphs.c", o[0], 19.0 / 54, 1e-15);
  tk_eq(t, "insphs.p", o[1], 19.0 / 54, 1e-15);
  tk_eq(t, "insphs.o", o[2], 16.0 / 54, 1e-15);
  // a large step drives the inspectors below 0: set to 0
  clab_insphs(T, 1, g, 1, 1, 100, x0, o);
  tk_eq(t, "insphs.clamp", o[2], 0, 0);
  tk_eq(t, "insphs.half", o[0], 0.5, 1e-15);
  clab_insphs(T, 1, g, 1, 0, 1, x0, o);
  tk_eq(t, "insphs.none", o[1], 1.0 / 3, 1e-15);
  const double bad[1] = {NAN}, xz[3] = {0, 0, 0}, xn[3] = {-1, 1, 1};
  tk_ieq(t, "insphs.nt", clab_insphs(T, 0, g, 1, 1, 1, x0, o), CLAB_EARG);
  tk_ieq(t, "insphs.ng", clab_insphs(T, 1, g, 0, 1, 1, x0, o), CLAB_EARG);
  tk_ieq(t, "insphs.steps", clab_insphs(T, 1, g, 1, -1, 1, x0, o), CLAB_EARG);
  tk_ieq(t, "insphs.rate", clab_insphs(T, 1, g, 1, 1, 0, x0, o), CLAB_EDOM);
  tk_ieq(t, "insphs.x0", clab_insphs(T, 1, g, 1, 1, 1, xz, o), CLAB_EDOM);
  tk_ieq(t, "insphs.xneg", clab_insphs(T, 1, g, 1, 1, 1, xn, o), CLAB_EDOM);
  tk_ieq(t, "insphs.T", clab_insphs(bad, 1, g, 1, 1, 1, x0, o), CLAB_EDOM);
  tk_ieq(t, "insphs.g", clab_insphs(T, 1, bad, 1, 1, 1, x0, o), CLAB_EDOM);
}

void colarr(tk &t) {
  const double a[2] = {0.3, 0.5}, e[2] = {0.2, 0.5}, bg[2] = {0.1, 0};
  double c[8], o[2];
  tk_ieq(t, "colarr.rc", clab_colarr(a, e, bg, 2, c, o), CLAB_OK);
  tk_eq(t, "colarr.or", o[0], 0.1, 1e-14);
  tk_eq(t, "colarr.or0", o[1], 0, 0);
  // cells .06, .24, .14, .056 over .496
  tk_eq(t, "colarr.c11", c[0], 0.06 / 0.496, 1e-14);
  tk_eq(t, "colarr.c00", c[6], 0.056 / 0.496, 1e-14);
  tk_eq(t, "colarr.c01", c[5], 1.0 / 3, 1e-15);
  const double one[1] = {1}, two[1] = {2};
  tk_ieq(t, "colarr.n", clab_colarr(a, e, bg, 0, c, o), CLAB_EARG);
  tk_ieq(t, "colarr.a", clab_colarr(one, e, bg, 1, c, o), CLAB_EDOM);
  tk_ieq(t, "colarr.e", clab_colarr(a, one, bg, 1, c, o), CLAB_EDOM);
  tk_ieq(t, "colarr.bg", clab_colarr(a, e, two, 1, c, o), CLAB_EDOM);
}

}  // namespace

void t_dtr(tk &t, int big) {
  (void)big;
  detrsp(t);
  dtrid(t);
  inspg(t);
  insphs(t);
  colarr(t);
}
