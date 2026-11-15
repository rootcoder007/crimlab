// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/sar.cpp: closed forms of search theory (Koopman 1956, 1957;
// Frost 1996) and the refusals of the kernels.
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "tchk.h"

void t_sar(tk &T, int) {
  const double pi = 3.14159265358979323846;
  // inverse cube half a sweep width off the track: 1 - exp(-1 / pi)
  double x = 0.5, p = 0;
  tk_ieq(T, "sar.lr", clab_sarlr(&x, 1, 1, 0, 1, 0, &p), CLAB_OK);
  tk_eq(T, "sar.lr1", p, 1 - std::exp(-1 / pi), 1e-15);
  tk_ieq(T, "sar.lrw", clab_sarlr(&x, 1, -1, 0, 1, 0, &p), CLAB_EARG);
  // random search 1 - exp(-C); inverse cube at C = 1 is erf(sqrt(pi) / 2)
  const double cv[2] = {1, 1};
  double pd[2];
  clab_sarpod(cv, 2, 0, 1, pd);
  tk_eq(T, "sar.rand", pd[0], 1 - std::exp(-1.0), 1e-15);
  clab_sarpod(cv, 1, 2, 1, pd);
  tk_eq(T, "sar.ic", pd[0], std::erf(std::sqrt(pi) / 2), 1e-14);
  tk_ieq(T, "sar.podm", clab_sarpod(cv, 1, 7, 1, pd), CLAB_EARG);
  // the inverse cube parallel sweep equals the closed form at W / S
  double ps = 0;
  clab_sarpsw(1, 1.3, 0, 1, 0, 2001, 60, &ps);
  tk_eq(T, "sar.psw", ps, std::erf(std::sqrt(pi) / 2 / 1.3), 1e-12);
  tk_ieq(T, "sar.pswS", clab_sarpsw(1, 0, 0, 1, 0, 11, 2, &ps), CLAB_EARG);
  // Koopman: u = 1 gives 1 - 2 / e
  double o[3];
  clab_sarcrc(1, 1, pi, o);
  tk_eq(T, "sar.crc", o[1], 1 - 2 / std::exp(1.0), 1e-15);
  tk_ieq(T, "sar.crcs", clab_sarcrc(0, 1, 1, o), CLAB_EARG);
  // Bayes on three areas: POS 0.55
  const double poc[3] = {0.5, 0.3, 0.2}, D[3] = {0.8, 0.5, 0};
  double pos, pa[3], up[3], cp[3];
  int cert;
  clab_sarscs(poc, 3, D, 1, &pos, pa, up, cp, &cert);
  tk_eq(T, "sar.scs", pos, 0.55, 1e-15);
  tk_eq(T, "sar.upd", up[2], 0.2 / 0.45, 1e-15);
  tk_ieq(T, "sar.scsn", clab_sarscs(poc, 0, D, 1, &pos, pa, up, cp, &cert), CLAB_EARG);
  // a point datum on one cell [-1, 1]^2: (2 Phi(1) - 1)^2
  const double e[2] = {-1, 1}, c0[2] = {0, 0};
  double g = 0;
  clab_sarmap(e, 2, e, 2, 0, c0, 1, c0, 64, nullptr, &g);
  tk_eq(T, "sar.map", g, std::pow(std::erf(1 / std::sqrt(2.0)), 2), 1e-15);
  tk_ieq(T, "sar.mapd", clab_sarmap(e, 2, e, 2, 9, c0, 1, c0, 64, nullptr, &g), CLAB_EARG);
  tk_ieq(T, "sar.area0", clab_sarmap(e, 2, e, 2, 2, c0, 1, c0, 64, nullptr, &g), CLAB_EDOM);
  // two equal cells share the effort
  const double q[2] = {0.5, 0.5}, a[2] = {1, 1};
  double z[2], cov[2], pod[2], ps2, lam;
  clab_saropt(q, a, 2, 1, 2, 1e-14, z, cov, pod, &ps2, &lam);
  tk_eq(T, "sar.opt", z[0], 1, 1e-10);
  tk_ieq(T, "sar.opt0", clab_saropt(q, a, 0, 1, 2, 1e-14, z, cov, pod, &ps2, &lam), CLAB_EARG);
  // Frost: an effort of 4 E^2 has the optimal factor 1.02
  const double z4 = 4;
  double sq[8];
  clab_sarsq(&z4, 1, 2, 1, sq);
  tk_eq(T, "sar.sq", std::round(sq[1] * 100) / 100, 1.02, 1e-12);
  tk_ieq(T, "sar.sq0", clab_sarsq(&z4, 0, 2, 1, sq), CLAB_EARG);
  // Gauss-Legendre weights add to 2; the figures are 36
  std::vector<double> t(16), w(16);
  clab_sargs(16, 1, t.data(), w.data());
  double sw = 0;
  for (double v : w) sw += v;
  tk_eq(T, "sar.gl", sw, 2, 1e-13);
  tk_ieq(T, "sar.gsn", clab_sargs(1, 1, t.data(), w.data()), CLAB_EARG);
  tk_ieq(T, "sar.fig", clab_sarfig(-1, nullptr), 36);
}
