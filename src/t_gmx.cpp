// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Self-test of clab_gmxlik (src/gmx.cpp): the gamma likelihood of peak
// heights by hand on one marker, a known contributor, one unknown summed
// over every genotype, coancestry, drop-in, stutter and degradation.
#include <cmath>

#include "crimlab.h"
#include "tchk.h"

namespace {

double ldg(double x, double a, double s) { return (a - 1) * std::log(x) - x / s - std::lgamma(a) - a * std::log(s); }
double lpg(double x, double a, double s) { return std::log(clab_pgamma(x / s, a, 1)); }

void run(tk &t) {
  const double mu = 1000, om = 0.25, at = 50, pc = 0.05, lam = 0.01;
  const double s = mu * om * om, a1 = 1 / (om * om);
  // one known contributor a/b, both peaks seen; position 2 is the rest
  const int np[1] = {3};
  const double p[3] = {0.2, 0.3, 0.5}, y[3] = {900, 1150, 0}, cn[3] = {0, 0, 0}, cl[1] = {0};
  const int kg[2] = {0, 1};
  double par[5] = {mu, om, 1, 0, 1}, ll = 0;
  tk_ieq(t, "gmx1", clab_gmxlik(1, np, p, y, nullptr, nullptr, 1, kg, 0, cn, cl, 0, at, pc, lam, 0, 0, par, &ll),
         CLAB_OK);
  tk_eq(t, "gmx2", ll, ldg(900, a1, s) + ldg(1150, a1, s) + std::log(1 - pc), 1e-13);
  // one unknown instead: a seen (900), b not seen; genotypes aa, ab, a2,
  // bb, b2, 22 by hand
  const double y2[3] = {900, 0, 0};
  tk_ieq(t, "gmx3", clab_gmxlik(1, np, p, y2, nullptr, nullptr, 0, nullptr, 1, cn, cl, 0, at, pc, lam, 0, 0, par, &ll),
         CLAB_OK);
  const double hom = ldg(900, 2 * a1, s) + std::log(1 - pc);
  const double het = ldg(900, a1, s) + lpg(at, a1, s) + std::log(1 - pc);
  const double din = std::log(pc * 0.2 * lam) - lam * (900 - at);
  const double e3 = std::exp(hom) * 0.04 + std::exp(het) * 2 * 0.2 * 0.3 + std::exp(het) * 2 * 0.2 * 0.5 +
                    std::exp(din + lpg(at, 2 * a1, s)) * (0.09 + 0.25) +
                    std::exp(din + 2 * lpg(at, a1, s)) * 2 * 0.3 * 0.5;
  tk_eq(t, "gmx4", ll, std::log(e3), 1e-13);
  // coancestry: the profile a/a conditions the unknown, theta 0.1; the
  // unknown must be a/a to explain one peak of height 2000 without drop-in
  const double y3[3] = {2000, 0, 0}, cn3[3] = {2, 0, 0}, cl3[1] = {2}, th = 0.1;
  clab_gmxlik(1, np, p, y3, nullptr, nullptr, 0, nullptr, 1, cn3, cl3, th, at, 0, lam, 0, 0, par, &ll);
  const double paa = (2 * th + (1 - th) * 0.2) * (3 * th + (1 - th) * 0.2) / ((1 + th) * (1 + 2 * th));
  const double pab = 2 * (2 * th + (1 - th) * 0.2) * ((1 - th) * 0.3) / ((1 + th) * (1 + 2 * th));
  const double pa2 = 2 * (2 * th + (1 - th) * 0.2) * ((1 - th) * 0.5) / ((1 + th) * (1 + 2 * th));
  const double e5 = paa * std::exp(ldg(2000, 2 * a1, s)) +
                    (pab + pa2) * std::exp(ldg(2000, a1, s) + lpg(at, a1, s));
  tk_eq(t, "gmx5", ll, std::log(e5), 1e-13);
  // two known contributors 0.7/0.3 with backward stutter xi and
  // degradation beta: positions 10 (size 150), 11 (154) and 9 (a stutter
  // position no one carries)
  const int np4[1] = {3}, st4[3] = {2, 0, -1}, kg4[4] = {0, 0, 0, 1};
  const double p4[3] = {0.4, 0.6, 0}, y4[3] = {1500, 300, 120}, z4[3] = {150, 154, 146};
  const double xi = 0.08, be = 0.7;
  double par4[6] = {mu, om, be, xi, 0.7, 0.3};
  tk_ieq(t, "gmx6", clab_gmxlik(1, np4, p4, y4, z4, st4, 2, kg4, 0, cn, cl, 0, at, pc, lam, 1, 1, par4, &ll),
         CLAB_OK);
  const double d10 = std::pow(be, 0.25), d11 = std::pow(be, 0.29);
  const double c10 = (0.7 * 2 + 0.3) * d10, c11 = 0.3 * d11;
  const double s10 = (1 - xi) * c10 + xi * c11, s11 = (1 - xi) * c11, s9 = xi * c10;
  tk_eq(t, "gmx7", ll,
        ldg(1500, s10 / (om * om), s) + ldg(300, s11 / (om * om), s) + ldg(120, s9 / (om * om), s) + std::log(1 - pc),
        1e-13);
  // refusals: a peak below the threshold, a proportion above 1
  const double y5[3] = {40, 0, 0};
  tk_ieq(t, "gmx8", clab_gmxlik(1, np, p, y5, nullptr, nullptr, 1, kg, 0, cn, cl, 0, at, pc, lam, 0, 0, par, &ll),
         CLAB_EDOM);
  double par6[5] = {mu, om, 1, 0, 1.5};
  tk_ieq(t, "gmx9", clab_gmxlik(1, np, p, y, nullptr, nullptr, 1, kg, 0, cn, cl, 0, at, pc, lam, 0, 0, par6, &ll),
         CLAB_EDOM);
}

}  // namespace

void t_gmx(tk &t, int) { run(t); }
