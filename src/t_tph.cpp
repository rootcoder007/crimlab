// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for the taphonomy kernels (src/tph.cpp) against
// algebra done by hand. Part of the self-test (selftest.cpp).
#include <cmath>

#include "crimlab.h"
#include "tchk.h"

namespace {

void clrilr(tk &t) {
  // one row (1, 2, 4): logs l1, l1 + log 2, l1 + 2 log 2
  const double x[3] = {1, 2, 4};
  double o[3];
  const double l2 = std::log(2.0);
  tk_ieq(t, "clrilr.rc", clab_clrilr(x, 1, 3, 0, 0, o), CLAB_OK);
  tk_eq(t, "clrilr.c0", o[0], -l2, 1e-15);
  tk_eq(t, "clrilr.c1", o[1], 0, 1e-15);
  tk_eq(t, "clrilr.c2", o[2], l2, 1e-15);
  clab_clrilr(x, 1, 3, 0, 1, o);
  tk_eq(t, "clrilr.i1", o[0], std::sqrt(0.5) * -l2, 1e-15);
  tk_eq(t, "clrilr.i2", o[1], std::sqrt(2.0 / 3) * (0.5 * l2 - 2 * l2), 1e-15);
  const double z[2] = {0, 1}, b[2] = {-1, 1};
  tk_ieq(t, "clrilr.zero", clab_clrilr(z, 1, 2, 0, 0, o), CLAB_EDOM);
  tk_ieq(t, "clrilr.zps", clab_clrilr(z, 1, 2, 1, 0, o), CLAB_OK);
  tk_eq(t, "clrilr.zv", o[1], 0.5 * std::log(2.0), 1e-15);
  tk_ieq(t, "clrilr.neg", clab_clrilr(b, 1, 2, 0, 0, o), CLAB_EDOM);
  tk_ieq(t, "clrilr.ps", clab_clrilr(x, 1, 3, -1, 0, o), CLAB_EDOM);
  tk_ieq(t, "clrilr.n", clab_clrilr(x, 0, 3, 0, 0, o), CLAB_EARG);
  tk_ieq(t, "clrilr.d", clab_clrilr(x, 1, 1, 0, 0, o), CLAB_EARG);
  tk_ieq(t, "clrilr.type", clab_clrilr(x, 1, 3, 0, 2, o), CLAB_EARG);
  tk_ieq(t, "clrilr.typen", clab_clrilr(x, 1, 3, 0, -1, o), CLAB_EARG);
}

void preslr(tk &t) {
  // x at the H1 mean: ll1 = -log(2 pi)/2 - log s; H2 one sd away
  const double x[1] = {2}, m1[1] = {2}, s1[1] = {1}, m2[1] = {1}, s2[1] = {1};
  double o[5];
  int b = 0, s = 0;
  tk_ieq(t, "preslr.rc", clab_preslr(x, 1, m1, s1, 1, m2, s2, 1, o, &b, &s), CLAB_OK);
  tk_eq(t, "preslr.llr", o[2], 0.5, 1e-15);
  tk_eq(t, "preslr.lr", o[3], std::exp(0.5), 1e-15);
  tk_ieq(t, "preslr.band", b, 1);
  tk_ieq(t, "preslr.side", s, 1);
  clab_preslr(x, 1, m2, s2, 1, m1, s1, 1, o, &b, &s);
  tk_ieq(t, "preslr.side2", s, 2);
  clab_preslr(x, 1, m1, s1, 1, m1, s1, 1, o, &b, &s);
  tk_ieq(t, "preslr.none", b, 0);
  tk_ieq(t, "preslr.side0", s, 0);
  // n values with per-value means: a ratio of 10^5 is very strong
  const double xs[2] = {0, 0}, ma[2] = {0, 0}, mb[2] = {3, 3.787157244}, one2[2] = {1, 1};
  clab_preslr(xs, 2, ma, one2, 2, mb, one2, 2, o, &b, &s);
  tk_eq(t, "preslr.l10", o[4], 0.5 * (9 + 3.787157244 * 3.787157244) / std::log(10.0), 1e-14);
  tk_ieq(t, "preslr.band5", b, 5);
  const double mc[2] = {10, 10};
  clab_preslr(xs, 2, ma, one2, 2, mc, one2, 2, o, &b, &s);
  tk_ieq(t, "preslr.band6", b, 6);
  const double bad[1] = {NAN}, z[1] = {0};
  tk_ieq(t, "preslr.n", clab_preslr(x, 0, m1, s1, 1, m2, s2, 1, o, &b, &s), CLAB_EARG);
  tk_ieq(t, "preslr.k1", clab_preslr(xs, 2, m1, s1, 3, m2, s2, 1, o, &b, &s), CLAB_EARG);
  tk_ieq(t, "preslr.k2", clab_preslr(xs, 2, m1, s1, 1, m2, s2, 3, o, &b, &s), CLAB_EARG);
  tk_ieq(t, "preslr.x", clab_preslr(bad, 1, m1, s1, 1, m2, s2, 1, o, &b, &s), CLAB_EDOM);
  tk_ieq(t, "preslr.m1", clab_preslr(x, 1, bad, s1, 1, m2, s2, 1, o, &b, &s), CLAB_EDOM);
  tk_ieq(t, "preslr.m2", clab_preslr(x, 1, m1, s1, 1, bad, s2, 1, o, &b, &s), CLAB_EDOM);
  tk_ieq(t, "preslr.s1", clab_preslr(x, 1, m1, z, 1, m2, s2, 1, o, &b, &s), CLAB_EDOM);
  tk_ieq(t, "preslr.s2", clab_preslr(x, 1, m1, s1, 1, m2, z, 1, o, &b, &s), CLAB_EDOM);
}

void tpdmc(tk &t) {
  // one stage, pres 0.5, dr = mr = 0.5: advance 0.25, mummify 0.25
  double P[9], N[1], B[2], tt[1];
  tk_ieq(t, "tpdmc.rc", clab_tpdmc(0.5, 0.5, 0.5, 1, P, N, B, tt), CLAB_OK);
  tk_eq(t, "tpdmc.stay", P[0], 0.5, 1e-15);
  tk_eq(t, "tpdmc.N", N[0], 2, 1e-15);
  tk_eq(t, "tpdmc.Bs", B[0], 0.5, 1e-15);
  tk_eq(t, "tpdmc.Bm", B[1], 0.5, 1e-15);
  tk_eq(t, "tpdmc.t", tt[0], 2, 1e-15);
  tk_eq(t, "tpdmc.abs", P[8], 1, 0);
  // advance 1/2, mummify 1/2: nothing stays
  double P2[16], N2[4], B2[4], t2[2];
  clab_tpdmc(0.5, 1, 1, 2, P2, N2, B2, t2);
  tk_eq(t, "tpdmc.norm", P2[0], 0, 1e-15);
  // from stage 1: mummify 1/2 at once, else stage 2 and mummify 1/2 again
  tk_eq(t, "tpdmc.B2m", B2[2], 0.75, 1e-15);
  tk_eq(t, "tpdmc.B2s", B2[0], 0.25, 1e-15);
  tk_eq(t, "tpdmc.t2", t2[0], 1.5, 1e-15);
  tk_ieq(t, "tpdmc.k", clab_tpdmc(0.5, 0.5, 0.5, 0, P, N, B, tt), CLAB_EARG);
  tk_ieq(t, "tpdmc.stuck", clab_tpdmc(1, 0.5, 0, 1, P, N, B, tt), CLAB_EDOM);
  tk_ieq(t, "tpdmc.pres", clab_tpdmc(2, 0.5, 0.5, 1, P, N, B, tt), CLAB_EDOM);
  tk_ieq(t, "tpdmc.dr", clab_tpdmc(0.5, 0, 0.5, 1, P, N, B, tt), CLAB_EDOM);
  tk_ieq(t, "tpdmc.mr", clab_tpdmc(0.5, 0.5, -1, 1, P, N, B, tt), CLAB_EDOM);
  // a path ends absorbed, its length at most steps + 1
  int path[200], len = 0;
  tk_ieq(t, "tpdsim.rc", clab_tpdsim(P2, 4, 2, 0, 199, 7, path, &len), CLAB_OK);
  tk_true(t, "tpdsim.end", path[len - 1] >= 2);
  tk_true(t, "tpdsim.len", len >= 2 && len <= 3);
  clab_tpdsim(P2, 4, 2, 0, 0, 7, path, &len);
  tk_ieq(t, "tpdsim.zero", len, 1);
  tk_ieq(t, "tpdsim.m", clab_tpdsim(P2, 1, 1, 0, 5, 1, path, &len), CLAB_EARG);
  tk_ieq(t, "tpdsim.k", clab_tpdsim(P2, 4, 0, 0, 5, 1, path, &len), CLAB_EARG);
  tk_ieq(t, "tpdsim.km", clab_tpdsim(P2, 4, 4, 0, 5, 1, path, &len), CLAB_EARG);
  tk_ieq(t, "tpdsim.start", clab_tpdsim(P2, 4, 2, 2, 5, 1, path, &len), CLAB_EARG);
  tk_ieq(t, "tpdsim.startn", clab_tpdsim(P2, 4, 2, -1, 5, 1, path, &len), CLAB_EARG);
  tk_ieq(t, "tpdsim.steps", clab_tpdsim(P2, 4, 2, 0, -1, 1, path, &len), CLAB_EARG);
}

void tpbhm(tk &t) {
  // y = 1 + 2 x exactly: rss 0, the noise variance falls back on var(y)
  const double x[4] = {0, 1, 2, 3}, y[4] = {1, 3, 5, 7}, m0[2] = {0, 0};
  const double big[2] = {1e8, 1e8};
  double mu[2], sd[2], st[2], w[10], gm[2], lam[2];
  int nj[2];
  tk_ieq(t, "tpbhm.rc", clab_tpbhm(y, x, 4, 1, m0, big, nullptr, 0, mu, sd, st, gm, lam, nj, w),
         CLAB_OK);
  tk_eq(t, "tpbhm.s2", st[0], 20.0 / 3, 1e-14);
  tk_eq(t, "tpbhm.b0", mu[0], 1, 1e-6);
  tk_eq(t, "tpbhm.b1", mu[1], 2, 1e-6);
  tk_nan(t, "tpbhm.tau", st[1]);
  // groups: residuals about the line
  const double y2[4] = {1.5, 2.5, 5.5, 6.5};
  const int g[4] = {0, 1, 0, 1};
  clab_tpbhm(y2, x, 4, 1, m0, big, g, 2, mu, sd, st, gm, lam, nj, w);
  // least squares: slope 1.8, intercept 1.3; residuals .2 -.6 .6 -.2
  tk_eq(t, "tpbhm.ls1", mu[1], 1.8, 1e-6);
  tk_eq(t, "tpbhm.gm0", gm[0], 0.4, 1e-6);
  tk_ieq(t, "tpbhm.nj", nj[1], 2);
  // var(gm) = 0.32, mean(s2 / nj) = 0.4 / 2 / 2... s2 = 0.8 / 2 = 0.4
  tk_eq(t, "tpbhm.tau2", st[1], 0.32 - 0.2, 1e-6);
  tk_eq(t, "tpbhm.lam", lam[0], 0.12 / 0.32, 1e-5);
  const double xd[8] = {0, 1, 2, 3, 0, 2, 4, 6}, yc[4] = {2, 2, 2, 2};
  const int gb[4] = {0, 2, 0, 0}, g1[4] = {0, 0, 0, 2};
  const double nan1[4] = {0, NAN, 1, 2}, s0b[2] = {0, 1}, m0b[2] = {NAN, 0};
  tk_ieq(t, "tpbhm.dep", clab_tpbhm(y, xd, 4, 2, m0, big, nullptr, 0, mu, sd, st, gm, lam, nj, w),
         CLAB_EDOM);
  tk_ieq(t, "tpbhm.flat", clab_tpbhm(yc, x, 4, 1, m0, big, nullptr, 0, mu, sd, st, gm, lam, nj, w),
         CLAB_EDOM);
  tk_ieq(t, "tpbhm.fewer", clab_tpbhm(y, x, 2, 2, m0, big, nullptr, 0, mu, sd, st, gm, lam, nj, w),
         CLAB_EDOM);
  tk_ieq(t, "tpbhm.n", clab_tpbhm(y, x, 1, 1, m0, big, nullptr, 0, mu, sd, st, gm, lam, nj, w),
         CLAB_EARG);
  tk_ieq(t, "tpbhm.q", clab_tpbhm(y, x, 4, 0, m0, big, nullptr, 0, mu, sd, st, gm, lam, nj, w),
         CLAB_EARG);
  tk_ieq(t, "tpbhm.ng", clab_tpbhm(y, x, 4, 1, m0, big, g, 1, mu, sd, st, gm, lam, nj, w),
         CLAB_EARG);
  tk_ieq(t, "tpbhm.g", clab_tpbhm(y, x, 4, 1, m0, big, gb, 2, mu, sd, st, gm, lam, nj, w),
         CLAB_EARG);
  tk_ieq(t, "tpbhm.unused", clab_tpbhm(y, x, 4, 1, m0, big, g1, 3, mu, sd, st, gm, lam, nj, w),
         CLAB_EARG);
  tk_ieq(t, "tpbhm.ynan", clab_tpbhm(nan1, x, 4, 1, m0, big, nullptr, 0, mu, sd, st, gm, lam, nj, w),
         CLAB_EDOM);
  tk_ieq(t, "tpbhm.xnan", clab_tpbhm(y, nan1, 4, 1, m0, big, nullptr, 0, mu, sd, st, gm, lam, nj, w),
         CLAB_EDOM);
  tk_ieq(t, "tpbhm.s0", clab_tpbhm(y, x, 4, 1, m0, s0b, nullptr, 0, mu, sd, st, gm, lam, nj, w),
         CLAB_EDOM);
  tk_ieq(t, "tpbhm.m0", clab_tpbhm(y, x, 4, 1, m0b, big, nullptr, 0, mu, sd, st, gm, lam, nj, w),
         CLAB_EDOM);
}

}  // namespace

void t_tph(tk &t, int big) {
  (void)big;
  clrilr(t);
  preslr(t);
  tpdmc(t);
  tpbhm(t);
}
