// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests of the journey-to-crime kernels (src/jtc.cpp).
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "cgeo.h"
#include "crimlab.h"
#include "tchk.h"

namespace {

// a message of the extended run (standalone, "ext"); silent by default
template <class... A>
void info(int big, const char *f, A... a) {
  if (big) std::printf(f, a...);
}

// Reference values from dev/ref/cpp/jtc.R (base R).
// R 4.5.0 
// 2.5 * dnorm(d, 4.2, 4.6)
const double rnorm[] = {
    0.14291170293506861, 0.14574222134047793, 0.15689290855517818,
    0.17021889060715442, 0.20250446773737507, 0.21681645673990912,
    0.1801511058496213, 0.042572803563679804, 3.2002955975362278e-08};
// 8.6 * dlnorm(d, 0.7, 0.9)
const double rlnorm[] = {
    0, 0.14596844358379713, 2.3008306593828585, 2.8171198969332711,
    1.4814419043110811, 0.65021690934028997, 0.20889345053906141,
    0.038963582819648603, 0.0014060015921874094};
// 1.89 * dexp(d, 0.06) / 0.06
const double rnexp[] = {
    1.8899999999999999, 1.8786939520619375, 1.8341420584066805,
    1.7799349684742301, 1.6267380754433591, 1.4689925549503082,
    1.2418184894504571, 0.89277278468051757, 0.31241489873879852};
// trips tr = ((i^2 31 + 7 i) mod 97) / 10, i = 0..59; t = 0, 0.5, ..., 10; h = 0.8
// normal kernel: sapply(t, function(s) sum(dnorm((s - tr) / h)) / (60 h))
const double kdn[] = {
    0.039552566407951718, 0.048067542743273323, 0.053812035621693725,
    0.063426296970558238, 0.076297572826433596, 0.08697097525670866,
    0.094803708078795304, 0.10182933681060218, 0.10577653910966357,
    0.10488396541790256, 0.10188260882695073, 0.097045495160413817,
    0.090546219492929569, 0.090094618188444262, 0.10183791440047169,
    0.12102387109099223, 0.13819893318318355, 0.14378385608694599,
    0.12817035740107871, 0.091607952721415586, 0.050012125957499792};
// reflected: adds dnorm((s + tr) / h)
const double kdnr[] = {
    0.079105132815903437, 0.073680738218569339, 0.065813843689874738,
    0.06735963535538475, 0.077187146871990422, 0.087109019818988734,
    0.094818366237732696, 0.10183040024507518, 0.1057765917719382,
    0.10488396719674369, 0.10188260886791052, 0.097045495161056358,
    0.090546219492936411, 0.090094618188444317, 0.10183791440047169,
    0.12102387109099223, 0.13819893318318355, 0.14378385608694599,
    0.12817035740107871, 0.091607952721415586, 0.050012125957499792};
// quartic kernel 15/16 (1 - u^2)^2
const double kdq[] = {
    0.068287849426269531, 0.064253807067871094, 0.023736953735351562,
    0.053658485412597649, 0.089273452758789049, 0.10317802429199219,
    0.076313018798828153, 0.10164260864257814, 0.14154911041259763,
    0.088758468627929674, 0.083799362182617243, 0.14283657073974609,
    0.071210861206054674, 0.052742958068847649, 0.09507656097412108,
    0.13858318328857422, 0.13369560241699222, 0.16872406005859378,
    0.19220829010009766, 0.086317062377929701, 0.014724731445312498};

const double PI = cgeo::pi;

void t_decay(tk &T, int /*big*/) {
  const double d[9] = {0, 0.1, 0.5, 1, 2.5, 4.2, 7, 12.5, 30};
  double f[9];
  const double pn[3] = {2.5, 4.2, 4.6};
  tk_true(T, "dec.ok", clab_jtdecay(d, 9, 2, pn, nullptr, 0, f) == CLAB_OK);
  for (int i = 0; i < 9; ++i) tk_eq(T, "dec.normal", f[i], rnorm[i], 1e-14);
  const double pl[3] = {8.6, 0.7, 0.9};
  clab_jtdecay(d, 9, 3, pl, nullptr, 0, f);
  for (int i = 0; i < 9; ++i) tk_eq(T, "dec.lognormal", f[i], rlnorm[i], 1e-14);
  const double pe[2] = {1.89, 0.06};
  clab_jtdecay(d, 9, 1, pe, nullptr, 0, f);
  for (int i = 0; i < 9; ++i) tk_eq(T, "dec.negexp", f[i], rnexp[i], 1e-14);
  // linear, cut at zero
  const double pa[2] = {1.9, -0.06};
  clab_jtdecay(d, 9, 0, pa, nullptr, 0, f);
  tk_eq(T, "dec.lin0", f[0], 1.9, 1e-15);
  tk_eq(T, "dec.lin7", f[6], 1.9 - 0.42, 1e-15);
  tk_eq(T, "dec.lin30", f[8], 1.9 - 1.8, 1e-14);
  const double far[1] = {40};
  double fb;
  clab_jtdecay(far, 1, 0, pa, nullptr, 0, &fb);
  tk_true(T, "dec.lin40", fb == 0.0);
  // truncated: rises linearly to P at dp, continuous, then decays
  const double pt[3] = {13.8, 0.4, 0.2};
  const double dt[5] = {0, 0.2, 0.4, 0.4 + 1e-12, 5.4};
  clab_jtdecay(dt, 5, 4, pt, nullptr, 0, f);
  tk_true(T, "dec.tr0", f[0] == 0.0);
  tk_eq(T, "dec.trmid", f[1], 6.9, 1e-14);
  tk_eq(T, "dec.trpeak", f[2], 13.8, 1e-14);
  tk_eq(T, "dec.trcont", f[3], 13.8, 1e-11);
  tk_eq(T, "dec.trtail", f[4], 13.8 * std::exp(-1.0), 1e-14);
  // table: exact at nodes, linear in between, 0 beyond
  const double tab[4] = {1, 3, 2, 0.5}, ptab[1] = {0.5};
  const double dd[6] = {0, 0.5, 0.75, 1.25, 1.5, 1.51};
  clab_jtdecay(dd, 6, 5, ptab, tab, 4, f);
  tk_eq(T, "dec.t0", f[0], 1, 1e-15);
  tk_eq(T, "dec.t1", f[1], 3, 1e-15);
  tk_eq(T, "dec.t15", f[2], 2.5, 1e-15);
  tk_eq(T, "dec.t25", f[3], 1.25, 1e-15);
  tk_eq(T, "dec.tend", f[4], 0.5, 1e-15);
  tk_true(T, "dec.tbeyond", f[5] == 0.0);
  // guards
  const double neg[1] = {-1};
  tk_true(T, "dec.negd", clab_jtdecay(neg, 1, 1, pe, nullptr, 0, f) == CLAB_EDOM);
  const double pbad[3] = {1, 0, 0};
  tk_true(T, "dec.sd0", clab_jtdecay(d, 9, 2, pbad, nullptr, 0, f) == CLAB_EDOM);
  tk_true(T, "dec.type", clab_jtdecay(d, 9, 6, pe, nullptr, 0, f) == CLAB_EARG);
  tk_true(T, "dec.tabnull", clab_jtdecay(d, 9, 5, ptab, nullptr, 4, f) == CLAB_EDOM);
}

void t_cal(tk &T, int big) {
  double tr[60], t[21], dn[21];
  for (int i = 0; i < 60; ++i) tr[i] = static_cast<double>((i * i * 31 + 7 * i) % 97) / 10;
  for (int j = 0; j < 21; ++j) t[j] = 0.5 * j;
  tk_true(T, "cal.ok", clab_jtcal(tr, nullptr, 60, 0, 0.8, 0, t, 21, dn) == CLAB_OK);
  for (int j = 0; j < 21; ++j) tk_eq(T, "cal.normal", dn[j], kdn[j], 1e-14);
  clab_jtcal(tr, nullptr, 60, 0, 0.8, 1, t, 21, dn);
  for (int j = 0; j < 21; ++j) tk_eq(T, "cal.reflect", dn[j], kdnr[j], 1e-14);
  clab_jtcal(tr, nullptr, 60, 1, 0.8, 0, t, 21, dn);
  for (int j = 0; j < 21; ++j) tk_eq(T, "cal.quartic", dn[j], kdq[j], 1e-14);
  // each kernel integrates to 1 (over [0, inf) when reflected): trapezoid
  // rule on a fine grid
  static double tt[40001], de[40001];
  for (int kern = 0; kern < 5; ++kern)
    for (int refl = 0; refl < 2; ++refl) {
      const double lo = refl ? 0.0 : -30.0, hi = 40.0;
      const int k = big ? 40001 : 4001;
      for (int j = 0; j < k; ++j) tt[j] = lo + (hi - lo) * j / (k - 1);
      clab_jtcal(tr, nullptr, 60, kern, 0.8, refl, tt, k, de);
      double s = 0;
      for (int j = 0; j < k - 1; ++j) s += 0.5 * (de[j] + de[j + 1]);
      s *= (hi - lo) / (k - 1);
      // trapezoid error: O(step^2) for the smooth kernels, O(step) for the
      // jumps of the uniform one
      const double tol = kern == 3 ? (big ? 2e-3 : 2e-2) : (big ? 2e-3 : 1e-3);
      tk_eq(T, "cal.mass", s, 1.0, kern < 2 || kern == 4 ? tol / 20 : tol);
    }
  // integer weights equal replicated trips
  double w[60], rep[180], a[21], b[21];
  int m = 0;
  for (int i = 0; i < 60; ++i) {
    w[i] = 1 + i % 3;
    for (int r = 0; r < w[i]; ++r) rep[m++] = tr[i];
  }
  clab_jtcal(tr, w, 60, 4, 0.5, 1, t, 21, a);
  clab_jtcal(rep, nullptr, m, 4, 0.5, 1, t, 21, b);
  for (int j = 0; j < 21; ++j) tk_eq(T, "cal.weights", a[j], b[j], 1e-14);
  // guards
  tk_true(T, "cal.h", clab_jtcal(tr, nullptr, 60, 0, 0, 0, t, 21, dn) == CLAB_EDOM);
  double tn[1] = {-1};
  tk_true(T, "cal.reflneg", clab_jtcal(tr, nullptr, 60, 0, 1, 1, tn, 1, dn) == CLAB_EDOM);
  tk_true(T, "cal.kern", clab_jtcal(tr, nullptr, 60, 5, 1, 0, t, 21, dn) == CLAB_EARG);
}

void t_surf(tk &T, int big) {
  static double raw[40000], prob[40000], raw2[40000], prob2[40000];
  double peak[3];
  // independent route: distances in the test, f from clab_jtdecay
  {
    const double x[4] = {2.2, 7.9, 5.1, 3.3}, y[4] = {1.4, 6.6, 9.0, 4.4};
    const double p[2] = {1.0, 0.4};
    tk_true(T, "srf.ok", clab_jtsurf(x, y, 4, 12, 10, 0, 0, 1, 1, 0, 0, 1, p, nullptr,
                             0, raw, prob, peak) == CLAB_OK);
    double s[120], tot = 0;
    for (int c = 0; c < 120; ++c) {
      const double cx = (c % 12) + 0.5, cy = (c / 12) + 0.5;
      double dd[4], f[4];
      for (int i = 0; i < 4; ++i) dd[i] = std::hypot(x[i] - cx, y[i] - cy);
      clab_jtdecay(dd, 4, 1, p, nullptr, 0, f);
      s[c] = f[0] + f[1] + f[2] + f[3];
      tot += s[c];
    }
    for (int c = 0; c < 120; ++c) {
      tk_eq(T, "srf.sum", raw[c], s[c], 1e-14);
      tk_eq(T, "srf.prob", prob[c], s[c] / tot, 1e-13);
    }
    // product mode against the direct product
    clab_jtsurf(x, y, 4, 12, 10, 0, 0, 1, 1, 0, 1, 1, p, nullptr, 0, raw2,
                prob2, peak);
    double q[120], qt = 0;
    for (int c = 0; c < 120; ++c) {
      const double cx = (c % 12) + 0.5, cy = (c / 12) + 0.5;
      q[c] = 1;
      for (int i = 0; i < 4; ++i)
        q[c] *= std::exp(-0.4 * std::hypot(x[i] - cx, y[i] - cy));
      qt += q[c];
    }
    for (int c = 0; c < 120; ++c) tk_eq(T, "srf.prod", prob2[c], q[c] / qt, 1e-12);
    // Manhattan
    clab_jtsurf(x, y, 4, 12, 10, 0, 0, 1, 1, 1, 0, 1, p, nullptr, 0, raw,
                prob, peak);
    const double cx = 3.5, cy = 2.5;
    double mm = 0;
    for (int i = 0; i < 4; ++i)
      mm += std::exp(-0.4 * (std::fabs(x[i] - cx) + std::fabs(y[i] - cy)));
    tk_eq(T, "srf.manhattan", raw[2 * 12 + 3], mm, 1e-14);
  }
  // incidents on a circle of radius 3 about a cell centre: the likelihood
  // surface peaks at the centre (it minimises the summed distance there,
  // the Weber point of the ring); the sum surface does not: exp(-d / 2) is
  // convex, so it peaks next to the ring
  {
    double x[8], y[8];
    for (int i = 0; i < 8; ++i) {
      x[i] = 5.5 + 3 * std::cos(2 * PI * i / 8);
      y[i] = 5.5 + 3 * std::sin(2 * PI * i / 8);
    }
    const double p[2] = {1.0, 0.5};
    clab_jtsurf(x, y, 8, 11, 11, 0, 0, 1, 1, 0, 1, 1, p, nullptr, 0, raw,
                prob, peak);
    tk_true(T, "srf.centre", peak[0] == 5 * 11 + 5 && peak[1] == 5.5 && peak[2] == 5.5);
    clab_jtsurf(x, y, 8, 11, 11, 0, 0, 1, 1, 0, 0, 1, p, nullptr, 0, raw,
                prob, peak);
    const double rr = std::hypot(peak[1] - 5.5, peak[2] - 5.5);
    tk_true(T, "srf.ring", rr > 2 && rr < 4);
  }
  // many incidents: the direct product underflows, the surface does not
  {
    cgeo::rng r(5);
    static double x[3000], y[3000];
    for (int i = 0; i < 3000; ++i) {
      x[i] = 10 + 4 * (r.unif() - 0.5);
      y[i] = 12 + 4 * (r.unif() - 0.5);
    }
    const double p[2] = {1.0, 1.0};
    tk_true(T, "srf.big", clab_jtsurf(x, y, 3000, 30, 30, 0, 0, 1, 1, 0, 1, 1, p,
                              nullptr, 0, raw, prob, peak) == CLAB_OK);
    double s = 0;
    bool fin = true;
    for (int c = 0; c < 900; ++c) {
      s += prob[c];
      fin = fin && std::isfinite(prob[c]);
    }
    tk_true(T, "srf.bigfin", fin);
    tk_eq(T, "srf.bigsum", s, 1, 1e-13);
    tk_true(T, "srf.bigpeak", std::fabs(peak[1] - 10) < 1.5 && std::fabs(peak[2] - 12) < 1.5);
  }
  // translation: grid and incidents moved by 1e5
  {
    const double x[3] = {1.2, 4.4, 2.0}, y[3] = {2.5, 1.0, 3.3};
    const double u[3] = {1.2 + 1e5, 4.4 + 1e5, 2.0 + 1e5},
                 v[3] = {2.5 - 1e5, 1.0 - 1e5, 3.3 - 1e5};
    const double p[3] = {1.0, 1.5, 1.0};
    clab_jtsurf(x, y, 3, 6, 5, 0, 0, 1, 1, 0, 0, 2, p, nullptr, 0, raw, prob,
                peak);
    clab_jtsurf(u, v, 3, 6, 5, 1e5, -1e5, 1, 1, 0, 0, 2, p, nullptr, 0, raw2,
                prob2, peak);
    for (int c = 0; c < 30; ++c) tk_eq(T, "srf.shift", prob2[c], prob[c], 1e-9);
  }
  // Calibration of the product surface: when the decay is the planar
  // density of an incident about the base and the base is uniform over the
  // cells, prob is the exact posterior, so a region holding posterior mass
  // q contains the base with probability E[q]. Bases on cell centres of a
  // 40 x 40 grid, 5 incidents each with planar density proportional to
  // exp(-r / 3) (radius Gamma(2, 3), uniform direction).
  {
    cgeo::rng r(77);
    const int G = 40, sims = big ? 300 : 100;
    const double s = 3.0, p[2] = {1.0, 1.0 / s};
    const double lev[2] = {0.5, 0.9};
    int hit[2] = {0, 0};
    double mass[2] = {0, 0}, hsp = 0, hss = 0;
    static int ord[1600];
    for (int k = 0; k < sims; ++k) {
      const int home = r.below(G * G);
      const double hx = home % G + 0.5, hy = home / G + 0.5;
      double x[5], y[5];
      for (int i = 0; i < 5; ++i) {
        double u1 = r.unif(), u2 = r.unif();
        while (u1 == 0) u1 = r.unif();
        while (u2 == 0) u2 = r.unif();
        const double rad = -s * (std::log(u1) + std::log(u2)), a = 2 * PI * r.unif();
        x[i] = hx + rad * std::cos(a);
        y[i] = hy + rad * std::sin(a);
      }
      clab_jtsurf(x, y, 5, G, G, 0, 0, 1, 1, 0, 1, 1, p, nullptr, 0, raw, prob,
                  peak);
      for (int c = 0; c < G * G; ++c) ord[c] = c;
      cgeo::hsort(ord, G * G, [](int a, int b) {
        return prob[a] == prob[b] ? a < b : prob[a] > prob[b];
      });
      for (int l = 0; l < 2; ++l) {
        double acc = 0;
        bool in = false;
        for (int c = 0; c < G * G && acc < lev[l]; ++c) {
          acc += prob[ord[c]];
          if (ord[c] == home) in = true;
        }
        mass[l] += acc;
        if (in) ++hit[l];
      }
      double acc[6];
      clab_jtacc(prob, G, G, 0, 0, 1, 1, hx, hy, acc);
      hsp += acc[2];
      clab_jtsurf(x, y, 5, G, G, 0, 0, 1, 1, 0, 0, 1, p, nullptr, 0, raw2,
                  prob2, peak);
      clab_jtacc(prob2, G, G, 0, 0, 1, 1, hx, hy, acc);
      hss += acc[2];
    }
    for (int l = 0; l < 2; ++l) {
      const double m = mass[l] / sims, c = static_cast<double>(hit[l]) / sims;
      info(big, "  product surface: region of mean mass %.3f covers the base "
                  "in %.3f of %d series\n",
                  m, c, sims);
      tk_true(T, "srf.calib", std::fabs(c - m) < 4 * std::sqrt(m * (1 - m) / sims));
    }
    info(big, "  mean hit score (share searched): product %.4f, sum %.4f\n",
                hsp / sims, hss / sims);
    tk_true(T, "srf.hits", hsp / sims < 0.1 && hss / sims < 0.1);
  }
  // guards
  const double x[1] = {1}, y[1] = {1}, p[2] = {1, -40};
  tk_true(T, "srf.zero", clab_jtsurf(x, y, 1, 3, 3, 10, 10, 1, 1, 0, 0, 0, p, nullptr,
                             0, raw, prob, peak) == CLAB_EDOM);
  tk_true(T, "srf.zerolog", clab_jtsurf(x, y, 1, 3, 3, 10, 10, 1, 1, 0, 1, 0, p,
                                nullptr, 0, raw, prob, peak) == CLAB_EDOM);
  const double pe[2] = {1, 1};
  tk_true(T, "srf.grid", clab_jtsurf(x, y, 1, 0, 3, 0, 0, 1, 1, 0, 0, 1, pe, nullptr,
                             0, raw, prob, peak) == CLAB_EARG);
  tk_true(T, "srf.dx", clab_jtsurf(x, y, 1, 3, 3, 0, 0, 0, 1, 0, 0, 1, pe, nullptr, 0,
                           raw, prob, peak) == CLAB_EARG);
  tk_true(T, "srf.m", clab_jtsurf(x, y, 0, 3, 3, 0, 0, 1, 1, 0, 0, 1, pe, nullptr, 0,
                          raw, prob, peak) == CLAB_EARG);
}

void t_acc(tk &T, int /*big*/) {
  // 3 x 2 grid of unit cells from (10, 20)
  const double pr[6] = {0.1, 0.3, 0.1, 0.3, 0.15, 0.05};
  double o[6];
  tk_true(T, "acc.ok", clab_jtacc(pr, 3, 2, 10, 20, 1, 1, 10.2, 20.7, o) == CLAB_OK);
  tk_true(T, "acc.cell", o[5] == 0 && o[0] == 0.1);
  // greater than 0.1: 0.3, 0.3, 0.15; at least 0.1: those and two 0.1s
  tk_eq(T, "acc.best", o[1], 4.0 / 6, 1e-15);
  tk_eq(T, "acc.worst", o[2], 5.0 / 6, 1e-15);
  // the first maximum is cell 1, centre (11.5, 20.5)
  tk_true(T, "acc.peak", o[4] == 1);
  tk_eq(T, "acc.err", o[3], std::hypot(11.5 - 10.2, 20.5 - 20.7), 1e-15);
  clab_jtacc(pr, 3, 2, 10, 20, 1, 1, 13, 22, o);  // the far corner: last cell
  tk_true(T, "acc.corner", o[5] == 5);
  tk_true(T, "acc.out", clab_jtacc(pr, 3, 2, 10, 20, 1, 1, 9.9, 20.5, o) == CLAB_EDOM);
}

void t_bjtc(tk &T, int /*big*/) {
  const double pj[3] = {0.2, 0.5, 0.3};
  const double od[6] = {5, 1, 0, 4, 3, 0};
  const int dest[3] = {0, 0, 1};
  double o[15];
  tk_true(T, "bj.ok", clab_bjtc(pj, 3, od, 2, dest, 3, o) == CLAB_OK);
  const double gen[3] = {6.0 / 13, 4.0 / 13, 3.0 / 13};
  const double cnd[3] = {11.0 / 21, 4.0 / 21, 6.0 / 21};
  double prd[3], pt = 0, rk[3], rt = 0;
  for (int i = 0; i < 3; ++i) {
    prd[i] = pj[i] * cnd[i];
    pt += prd[i];
  }
  for (int i = 0; i < 3; ++i) {
    prd[i] /= pt;
    rk[i] = prd[i] / gen[i];
    rt += rk[i];
  }
  for (int i = 0; i < 3; ++i) {
    tk_eq(T, "bj.jtc", o[i], pj[i], 1e-15);
    tk_eq(T, "bj.gen", o[3 + i], gen[i], 1e-15);
    tk_eq(T, "bj.cond", o[6 + i], cnd[i], 1e-15);
    tk_eq(T, "bj.prod", o[9 + i], prd[i], 1e-15);
    tk_eq(T, "bj.risk", o[12 + i], rk[i] / rt, 1e-15);
  }
  tk_eq(T, "bj.prodval", o[9], 2.2 / 6, 1e-15);
  // scaling the counts or the surface changes nothing
  double od7[6], pj3[3], o2[15];
  for (int i = 0; i < 6; ++i) od7[i] = 7 * od[i];
  for (int i = 0; i < 3; ++i) pj3[i] = 3 * pj[i];
  clab_bjtc(pj3, 3, od7, 2, dest, 3, o2);
  for (int i = 0; i < 15; ++i) tk_eq(T, "bj.scale", o2[i], o[i], 1e-15);
  // a zone where nobody lives: product and risk 0 there
  const double odz[8] = {5, 1, 0, 4, 3, 0, 0, 0};
  const double pjz[4] = {0.2, 0.4, 0.3, 0.1};
  double oz[20];
  clab_bjtc(pjz, 4, odz, 2, dest, 3, oz);
  tk_true(T, "bj.zero", oz[7] == 0 && oz[11] == 0 && oz[15] == 0 && oz[19] == 0);
  // no solved offence in the series' zones: undefined
  const double odn[6] = {5, 0, 0, 0, 3, 0};
  const int d1[1] = {1};
  tk_true(T, "bj.none", clab_bjtc(pj, 3, odn, 2, d1, 1, o) == CLAB_EDOM);
  const int dbad[1] = {2};
  tk_true(T, "bj.dest", clab_bjtc(pj, 3, od, 2, dbad, 1, o) == CLAB_EARG);
}

void famtest(tk &T, int big) {
  t_decay(T, big);
  t_cal(T, big);
  t_surf(T, big);
  t_acc(T, big);
  t_bjtc(T, big);
}

#ifdef CLAB_MAIN
double now() {
  return std::chrono::duration<double>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}

// a heap buffer for the timings (kept out of the program's static size)
template <class T>
struct hb {
  T *p;
  explicit hb(long n) : p(static_cast<T *>(std::malloc(sizeof(T) * n))) {
    if (p == nullptr) std::abort();
  }
  ~hb() { std::free(p); }
  hb(const hb &) = delete;
  hb &operator=(const hb &) = delete;
};

double med5(double *v) { return tk_q7(v, 5, 0.5); }

void fambench() {
  hb<double> raw_(1000000), prob_(1000000);
  double *raw = raw_.p, *prob = prob_.p;
  double x[20], y[20], peak[3];
  cgeo::rng r(3);
  const int g[3] = {32, 316, 1000};
  for (int G : g) {
    for (int i = 0; i < 20; ++i) {
      x[i] = G * (0.4 + 0.2 * r.unif());
      y[i] = G * (0.4 + 0.2 * r.unif());
    }
    const double p[3] = {13.8, 2.0, 0.2};
    double t[2][5];
    for (int k = 0; k < 5; ++k) {
      double a = now();
      clab_jtsurf(x, y, 20, G, G, 0, 0, 1, 1, 0, 0, 4, p, nullptr, 0, raw, prob,
                  peak);
      t[0][k] = now() - a;
      a = now();
      clab_jtsurf(x, y, 20, G, G, 0, 0, 1, 1, 0, 1, 4, p, nullptr, 0, raw, prob,
                  peak);
      t[1][k] = now() - a;
    }
    std::printf("bench jtc: 20 incidents, %d cells: sum %.5f product %.5f s "
                "(medians of 5)\n",
                G * G, med5(t[0]), med5(t[1]));
  }
}

#endif

}  // namespace

void t_jtc(tk &T, int big) { famtest(T, big); }
#ifdef CLAB_MAIN
void b_jtc() { fambench(); }
#endif
