// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for the hot-spot kernels (src/sdb.cpp): without
// offenders the attractiveness decays by 1 - omega dt a step; offenders
// at a site of certain burglary all burgle; a lone offender who cannot
// burgle moves to a neighbour; the continuum model rests at its uniform
// steady state; suppression stops crime; and the stability coefficients
// by hand. Part of the self-test (selftest.cpp).
#include <cmath>

#include "crimlab.h"
#include "tchk.h"

namespace {

struct grid {
  double b[25], rho[25], burg[5], mb[5], tot[25], w[75];
  int iw[25];
};

void agent(tk &t) {
  static grid g;
  double par[10] = {0.1, 0.2, 0.3, 0, 1, 1, 1, 1, 1, 4};
  for (int i = 0; i < 25; ++i) {
    g.b[i] = 2;
    g.rho[i] = 0;
  }
  tk_ieq(t, "sdb.decay.rc", clab_sdbsim(5, 3, 0, par, g.b, g.rho, nullptr, 1, g.burg, g.mb, g.tot, g.w, g.iw),
         CLAB_OK);
  tk_eq(t, "sdb.decay", g.b[7], 2 * 0.8 * 0.8 * 0.8, 1e-15);
  tk_eq(t, "sdb.decay.mean", g.mb[2], 2 * 0.512, 1e-15);
  tk_eq(t, "sdb.decay.burg", g.burg[2], 0, 0);
  // certain burglary: five offenders at site 12 all burgle
  par[0] = 50;
  for (int i = 0; i < 25; ++i) {
    g.b[i] = 0;
    g.rho[i] = i == 12 ? 5 : 0;
  }
  tk_ieq(t, "sdb.burg.rc", clab_sdbsim(5, 1, 0, par, g.b, g.rho, nullptr, 1, g.burg, g.mb, g.tot, g.w, g.iw),
         CLAB_OK);
  tk_eq(t, "sdb.burg", g.burg[0], 5, 0);
  tk_eq(t, "sdb.burg.tot", g.tot[12], 5, 0);
  tk_eq(t, "sdb.burg.b", g.b[12], 5, 1e-15);
  tk_eq(t, "sdb.burg.left", g.rho[12], 0, 0);
  // no burglary: the offender at 12 moves to 7, 11, 13 or 17
  par[0] = 1e-12;
  for (int i = 0; i < 25; ++i) {
    g.b[i] = 0;
    g.rho[i] = i == 12 ? 1 : 0;
  }
  par[3] = 0;
  clab_sdbsim(5, 1, 0, par, g.b, g.rho, nullptr, 3, g.burg, g.mb, g.tot, g.w, g.iw);
  tk_eq(t, "sdb.move", g.rho[7] + g.rho[11] + g.rho[13] + g.rho[17], 1, 0);
  // arrivals: Poisson(gam dt) per site, seeded
  par[3] = 2;
  double r1 = 0, r2 = 0;
  clab_sdbsim(5, 1, 0, par, g.b, g.rho, nullptr, 9, g.burg, g.mb, g.tot, g.w, g.iw);
  for (int i = 0; i < 25; ++i) r1 += g.rho[i];
  tk_gt(t, "sdb.arrive", r1, 10);
  for (int i = 0; i < 25; ++i) {
    g.b[i] = 0;
    g.rho[i] = 0;
  }
  clab_sdbsim(5, 1, 0, par, g.b, g.rho, nullptr, 9, g.burg, g.mb, g.tot, g.w, g.iw);
  for (int i = 0; i < 25; ++i) r2 += g.rho[i];
  tk_gt(t, "sdb.arrive2", r2, 10);
  tk_ieq(t, "sdb.n", clab_sdbsim(2, 1, 0, par, g.b, g.rho, nullptr, 9, g.burg, g.mb, g.tot, g.w, g.iw),
         CLAB_EARG);
  par[1] = 2;
  tk_ieq(t, "sdb.om", clab_sdbsim(5, 1, 0, par, g.b, g.rho, nullptr, 9, g.burg, g.mb, g.tot, g.w, g.iw),
         CLAB_EDOM);
  par[1] = 0.2;
  g.b[3] = -1;
  tk_ieq(t, "sdb.neg", clab_sdbsim(5, 1, 0, par, g.b, g.rho, nullptr, 9, g.burg, g.mb, g.tot, g.w, g.iw),
         CLAB_EDOM);
}

void pde(tk &t) {
  static grid g;
  // a0, omega, eta, gam, theta, dt, dx, dif, eps, z
  double par[10] = {0.5, 1, 0.2, 0.3, 0, 0.01, 1, 1, 2, 4};
  const double bb = 2 * 1 * 0.3 / 1, rb = 0.3 / (0.5 + bb);
  for (int i = 0; i < 25; ++i) {
    g.b[i] = bb;
    g.rho[i] = rb;
  }
  tk_ieq(t, "pde.rc", clab_sdbsim(5, 5, 1, par, g.b, g.rho, nullptr, 0, g.burg, g.mb, g.tot, g.w, g.iw),
         CLAB_OK);
  double dev = 0;
  for (int i = 0; i < 25; ++i) dev = std::fmax(dev, std::fabs(g.b[i] - bb) + std::fabs(g.rho[i] - rb));
  tk_lt(t, "pde.steady", dev, 1e-13);
  tk_eq(t, "pde.burg", g.burg[0], 25 * rb * (0.5 + bb) * 0.01, 1e-13);
  // a suppressed site has no crime
  int sup[25] = {0};
  sup[6] = 1;
  clab_sdbsim(5, 1, 1, par, g.b, g.rho, sup, 0, g.burg, g.mb, g.tot, g.w, g.iw);
  tk_eq(t, "pde.sup", g.tot[6], 0, 0);
  tk_gt(t, "pde.sup2", g.tot[5], 0);
  // too long a step, and a field driven negative
  par[5] = 1.5;
  tk_ieq(t, "pde.dt", clab_sdbsim(5, 1, 1, par, g.b, g.rho, nullptr, 0, g.burg, g.mb, g.tot, g.w, g.iw),
         CLAB_EDOM);
  par[5] = 1;
  par[1] = 0.5;
  for (int i = 0; i < 25; ++i) {
    g.b[i] = 0;
    g.rho[i] = i == 12 ? 1 : 0;
  }
  tk_ieq(t, "pde.neg", clab_sdbsim(5, 1, 1, par, g.b, g.rho, nullptr, 0, g.burg, g.mb, g.tot, g.w, g.iw),
         CLAB_EDOM);
}

void thr(tk &t) {
  // a0 1, omega 1, eps 2, dif 1, gam 1, eta 0.01, z 1: B = 2, A = 3,
  // rho = 1/3; a = -1/3, b = 6, c0 = -1/3, c1 = 2/9, d = -3, d1 = 0.01, d2 = 1
  const double par[7] = {1, 1, 2, 1, 1, 0.01, 1}, k[3] = {0, 1, 5};
  double st[12], tr[3], de[3], gr[3];
  tk_ieq(t, "thr.rc", clab_sdbthr(par, k, 3, st, tr, de, gr), CLAB_OK);
  const double p = 0.01, q = -1.0 / 3 - 0.03 + 6 * 2.0 / 9, r = 3.0;
  tk_eq(t, "thr.p", st[0], p, 1e-15);
  tk_eq(t, "thr.q", st[1], q, 1e-14);
  tk_eq(t, "thr.r", st[2], r, 1e-14);
  tk_eq(t, "thr.lhs", st[6], 2, 1e-15);
  tk_eq(t, "thr.rhs", st[7], 1.03, 1e-15);
  tk_eq(t, "thr.unst", st[8], 1, 0);
  tk_eq(t, "thr.kc", st[9], q / (2 * p), 1e-12);
  tk_eq(t, "thr.det1", de[1], p - q + r, 1e-14);
  tk_eq(t, "thr.tr0", tr[0], -1.0 / 3 - 3, 1e-15);
  tk_gt(t, "thr.grow", st[11], 0);
  // complex eigenvalues at k = 0: the growth is half the trace
  tk_eq(t, "thr.gr0", gr[0], tr[0] / 2, 1e-15);
  // no reinforcement: stable
  const double p0[7] = {1, 1, 0, 1, 1, 0.01, 1};
  clab_sdbthr(p0, k, 0, st, tr, de, gr);
  tk_eq(t, "thr.eps0", st[8], 0, 0);
  tk_nan(t, "thr.eps0.k", st[9]);
  // no spread of attractiveness: P = 0, unstable when Q > 0
  const double pe[7] = {1, 1, 2, 1, 1, 0, 1};
  clab_sdbthr(pe, k, 0, st, tr, de, gr);
  tk_eq(t, "thr.p0", st[8], 1, 0);
  const double bad[7] = {0, 1, 2, 1, 1, 0, 1};
  tk_ieq(t, "thr.bad", clab_sdbthr(bad, k, 0, st, tr, de, gr), CLAB_EDOM);
  const double kn[1] = {-1};
  tk_ieq(t, "thr.k", clab_sdbthr(par, kn, 1, st, tr, de, gr), CLAB_EDOM);
}

}  // namespace

void t_sdb(tk &t, int big) {
  (void)big;
  agent(t);
  pde(t);
  thr(t);
}
