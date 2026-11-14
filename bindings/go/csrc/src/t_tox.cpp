// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for the forensic toxicology kernels (src/tox.cpp)
// against sums done by hand. Part of the self-test (selftest.cpp).
#include <cmath>

#include "crimlab.h"
#include "tchk.h"

namespace {

void txcal(tk &t) {
  // y = 1 + 2 x exactly: zero residuals; unit weights
  const double x[4] = {1, 2, 3, 4}, y[4] = {3, 5, 7, 9}, w[4] = {1, 1, 1, 1};
  const double u[3] = {5, NAN, 1};
  double st[6], xh[3];
  int fl[3];
  tk_ieq(t, "txcal.rc", clab_txcal(x, y, w, 4, u, 3, 3.3, 10, st, xh, fl), CLAB_OK);
  tk_eq(t, "txcal.int", st[0], 1, 1e-14);
  tk_eq(t, "txcal.slope", st[1], 2, 1e-14);
  tk_eq(t, "txcal.r2", st[2], 1, 1e-14);
  tk_eq(t, "txcal.lod", st[4], 0, 1e-14);
  tk_eq(t, "txcal.xh", xh[0], 2, 1e-14);
  tk_ieq(t, "txcal.flagq", fl[0], 2);
  tk_ieq(t, "txcal.flagna", fl[1], -1);
  // Sxx 5, Sxy 2.5, Syy 2.75: slope 0.5, intercept 1.75 - 1.25, rss 1.5
  const double y2[4] = {1, 2, 1, 3};
  clab_txcal(x, y2, w, 4, u, 3, 3.3, 10, st, xh, fl);
  tk_eq(t, "txcal.slope2", st[1], 0.5, 1e-14);
  tk_eq(t, "txcal.int2", st[0], 0.5, 1e-14);
  tk_eq(t, "txcal.sigma2", st[3], std::sqrt(0.75), 1e-14);
  tk_eq(t, "txcal.r22", st[2], 1.25 / 2.75, 1e-14);
  tk_eq(t, "txcal.lod2", st[4], 3.3 * std::sqrt(0.75) / 0.5, 1e-14);
  tk_ieq(t, "txcal.below", fl[2], 0);
  const double u2[1] = {5.5};
  clab_txcal(x, y2, w, 4, u2, 1, 3.3, 10, st, xh, fl);
  tk_ieq(t, "txcal.loq", fl[0], 1);
  // flat response (slope 0) and flat concentrations refused
  const double yf[4] = {2, 2, 2, 2}, xf[4] = {1, 1, 1, 1}, wb[4] = {1, 0, 1, 1};
  tk_ieq(t, "txcal.flat", clab_txcal(x, yf, w, 4, u, 0, 3.3, 10, st, xh, fl), CLAB_EDOM);
  tk_ieq(t, "txcal.xflat", clab_txcal(xf, y, w, 4, u, 0, 3.3, 10, st, xh, fl), CLAB_EDOM);
  tk_ieq(t, "txcal.n", clab_txcal(x, y, w, 2, u, 0, 3.3, 10, st, xh, fl), CLAB_EARG);
  tk_ieq(t, "txcal.nu", clab_txcal(x, y, w, 4, u, -1, 3.3, 10, st, xh, fl), CLAB_EARG);
  tk_ieq(t, "txcal.k", clab_txcal(x, y, w, 4, u, 0, 0, 10, st, xh, fl), CLAB_EDOM);
  tk_ieq(t, "txcal.k2", clab_txcal(x, y, w, 4, u, 0, 3, 0, st, xh, fl), CLAB_EDOM);
  tk_ieq(t, "txcal.w", clab_txcal(x, y, wb, 4, u, 0, 3, 10, st, xh, fl), CLAB_EDOM);
  const double y3[4] = {3, 4, 6, 8};
  clab_txcal(x, y3, w, 4, u, 0, 3.3, 10, st, xh, fl);
  tk_eq(t, "txcal.r2b", st[2], 1 - 0.3 / 14.75, 1e-14);
}

void lodimp(tk &t) {
  const double x[5] = {0.4, 0.05, NAN, 1.2, 0.08}, l1[1] = {0.1};
  const double l5[5] = {0.5, 0.01, 0.1, 0.1, 0.1};
  double o[5], st[2];
  int c[5];
  tk_ieq(t, "lodimp.rc", clab_lodimp(x, 5, l1, 1, 0, o, c, st), CLAB_OK);
  tk_eq(t, "lodimp.half", o[1], 0.05, 0);
  tk_eq(t, "lodimp.keep", o[0], 0.4, 0);
  tk_nan(t, "lodimp.na", o[2]);
  tk_ieq(t, "lodimp.c2", c[2], 0);
  tk_eq(t, "lodimp.frac", st[0], 0.5, 0);
  tk_eq(t, "lodimp.miss", st[1], 1, 0);
  clab_lodimp(x, 5, l5, 5, 1, o, c, st);
  tk_eq(t, "lodimp.sqrt2", o[0], 0.5 / std::sqrt(2.0), 1e-16);
  tk_ieq(t, "lodimp.c1", c[1], 0);
  clab_lodimp(x, 5, l1, 1, 2, o, c, st);
  tk_eq(t, "lodimp.lod", o[4], 0.1, 0);
  const double xn[2] = {NAN, NAN};
  clab_lodimp(xn, 2, l1, 1, 2, o, c, st);
  tk_nan(t, "lodimp.nofrac", st[0]);
  const double lb[1] = {0};
  tk_ieq(t, "lodimp.n", clab_lodimp(x, 0, l1, 1, 0, o, c, st), CLAB_EARG);
  tk_ieq(t, "lodimp.nl", clab_lodimp(x, 5, l5, 2, 0, o, c, st), CLAB_EARG);
  tk_ieq(t, "lodimp.m", clab_lodimp(x, 5, l1, 1, 3, o, c, st), CLAB_EARG);
  tk_ieq(t, "lodimp.mneg", clab_lodimp(x, 5, l1, 1, -1, o, c, st), CLAB_EARG);
  tk_ieq(t, "lodimp.lod0", clab_lodimp(x, 5, lb, 1, 0, o, c, st), CLAB_EDOM);
}

void pmrrat(tk &t) {
  const double c[3] = {0.4, 1.5, 0.9}, p[3] = {0.5, 1, 0.2};
  double r[3];
  int k[3];
  tk_ieq(t, "pmrrat.rc", clab_pmrrat(c, p, 3, 1, 2, r, k), CLAB_OK);
  tk_eq(t, "pmrrat.r", r[2], 4.5, 1e-15);
  tk_ieq(t, "pmrrat.k0", k[0], 0);
  tk_ieq(t, "pmrrat.k1", k[1], 1);
  tk_ieq(t, "pmrrat.k2", k[2], 2);
  const double z[1] = {0};
  tk_ieq(t, "pmrrat.n", clab_pmrrat(c, p, 0, 1, 2, r, k), CLAB_EARG);
  tk_ieq(t, "pmrrat.c1", clab_pmrrat(c, p, 3, 0, 2, r, k), CLAB_EDOM);
  tk_ieq(t, "pmrrat.c2", clab_pmrrat(c, p, 3, 2, 1, r, k), CLAB_EDOM);
  tk_ieq(t, "pmrrat.c2inf", clab_pmrrat(c, p, 3, 1, INFINITY, r, k), CLAB_EDOM);
  tk_ieq(t, "pmrrat.c", clab_pmrrat(z, p, 1, 1, 2, r, k), CLAB_EDOM);
  tk_ieq(t, "pmrrat.p", clab_pmrrat(c, z, 1, 1, 2, r, k), CLAB_EDOM);
}

void etohlr(tk &t) {
  const double e[4] = {0.8, 0.5, 0.3, 0.2}, pr[4] = {0, 0.02, NAN, NAN};
  const double bu[4] = {0.01, 0.03, 0.05, NAN}, g[4] = {1.1, NAN, 0, NAN};
  double hi[4];
  int v[4];
  tk_ieq(t, "etohlr.rc", clab_etohlr(e, pr, bu, g, 4, 0, hi, v), CLAB_OK);
  tk_ieq(t, "etohlr.v0", v[0], 0);
  tk_ieq(t, "etohlr.v1", v[1], 1);
  tk_eq(t, "etohlr.hi1", hi[1], 0.03, 0);
  tk_eq(t, "etohlr.hi2", hi[2], 0.05, 0);
  tk_nan(t, "etohlr.hi3", hi[3]);
  tk_ieq(t, "etohlr.v3", v[3], 2);
  clab_etohlr(e, pr, bu, g, 4, 0.04, hi, v);
  tk_ieq(t, "etohlr.thr", v[1], 2);
  tk_eq(t, "etohlr.hi0", hi[0], 0.01, 0);
  const double eb[1] = {-1};
  tk_ieq(t, "etohlr.n", clab_etohlr(e, pr, bu, g, 0, 0, hi, v), CLAB_EARG);
  tk_ieq(t, "etohlr.thrinf", clab_etohlr(e, pr, bu, g, 1, NAN, hi, v), CLAB_EDOM);
  tk_ieq(t, "etohlr.eth", clab_etohlr(eb, pr, bu, g, 1, 0, hi, v), CLAB_EDOM);
}

void mtxrel(tk &t) {
  const double b[4] = {0.55, 0.95, 0.35, 0.8};
  const int d[4] = {1, 0, 1, 0};
  double r[4];
  int o[4];
  tk_ieq(t, "mtxrel.rc", clab_mtxrel(b, d, 4, 15, 2, 30, 0.8, 10, 0.6, r, o), CLAB_OK);
  tk_eq(t, "mtxrel.r0", r[0], 0.55 * 0.8 * 0.5, 1e-15);
  tk_eq(t, "mtxrel.r1", r[1], 0.95 * 0.8, 1e-15);
  tk_ieq(t, "mtxrel.o0", o[0], 1);
  tk_ieq(t, "mtxrel.o1", o[1], 3);
  tk_ieq(t, "mtxrel.o3", o[3], 2);
  // capped losses; ties keep their order
  const double bt[3] = {0.5, 0.5, 0.6};
  const int dt[3] = {0, 0, 1};
  clab_mtxrel(bt, dt, 3, 100, 100, 30, 0.8, 10, 0.6, r, o);
  tk_eq(t, "mtxrel.cap", r[2], 0.6 * 0.4 * 0.2, 1e-15);
  tk_ieq(t, "mtxrel.tie0", o[0], 0);
  tk_ieq(t, "mtxrel.tie1", o[1], 1);
  const double bn[1] = {-1};
  tk_ieq(t, "mtxrel.n", clab_mtxrel(b, d, 0, 0, 0, 30, 0.8, 10, 0.6, r, o), CLAB_EARG);
  tk_ieq(t, "mtxrel.days", clab_mtxrel(b, d, 4, -1, 0, 30, 0.8, 10, 0.6, r, o), CLAB_EDOM);
  tk_ieq(t, "mtxrel.stage", clab_mtxrel(b, d, 4, 0, -1, 30, 0.8, 10, 0.6, r, o), CLAB_EDOM);
  tk_ieq(t, "mtxrel.dmax", clab_mtxrel(b, d, 4, 0, 0, 0, 0.8, 10, 0.6, r, o), CLAB_EDOM);
  tk_ieq(t, "mtxrel.smax", clab_mtxrel(b, d, 4, 0, 0, 30, 0.8, 0, 0.6, r, o), CLAB_EDOM);
  tk_ieq(t, "mtxrel.dcap", clab_mtxrel(b, d, 4, 0, 0, 30, 1.5, 10, 0.6, r, o), CLAB_EDOM);
  tk_ieq(t, "mtxrel.dcap0", clab_mtxrel(b, d, 4, 0, 0, 30, -1, 10, 0.6, r, o), CLAB_EDOM);
  tk_ieq(t, "mtxrel.scap", clab_mtxrel(b, d, 4, 0, 0, 30, 0.8, 10, 2, r, o), CLAB_EDOM);
  tk_ieq(t, "mtxrel.scap0", clab_mtxrel(b, d, 4, 0, 0, 30, 0.8, 10, -1, r, o), CLAB_EDOM);
  tk_ieq(t, "mtxrel.base", clab_mtxrel(bn, d, 1, 0, 0, 30, 0.8, 10, 0.6, r, o), CLAB_EDOM);
}

}  // namespace

void t_tox(tk &t, int big) {
  (void)big;
  txcal(t);
  lodimp(t);
  pmrrat(t);
  etohlr(t);
  mtxrel(t);
}
