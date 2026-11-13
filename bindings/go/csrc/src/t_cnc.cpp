// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for the concentration kernels (src/cnc.cpp): small cases
// by hand, closed forms of the zeta function, the Poisson generator's
// Gini against the exact mean difference, and the Hill estimators written
// out. Part of the self-test (selftest.cpp).
#include <cmath>

#include "crimlab.h"
#include "tchk.h"

namespace {

const double pi = 3.14159265358979323846;

void conc(tk &t) {
  // 0, 0, 1, 3, 6: Gini 0.6, of the positive units 1/3, zero share 0.4
  const double x[5] = {6, 0, 1, 0, 3}, fr[4] = {0.2, 0.5, 0, 1};
  double st[16], lp[6], lv[6], top[4], tn[4];
  int o[5];
  tk_ieq(t, "cconc.rc", clab_cconc(x, nullptr, 5, fr, 4, st, lp, lv, top, tn, o), CLAB_OK);
  tk_eq(t, "cconc.gini", st[10], 0.6, 1e-15);
  tk_eq(t, "cconc.gpos", st[11], 1.0 / 3, 1e-15);
  tk_eq(t, "cconc.ident", st[10], st[7] + (1 - st[7]) * st[11], 1e-15);
  tk_eq(t, "cconc.var", st[3], 6.5, 1e-15);
  tk_eq(t, "cconc.vpop", st[13], 5.2, 1e-15);
  tk_eq(t, "cconc.disp", st[4], 3.25, 1e-15);
  tk_eq(t, "cconc.ivar", st[5], 4.5, 1e-15);
  tk_eq(t, "cconc.znull", st[8], std::exp(-2.0), 1e-15);
  tk_eq(t, "cconc.gnull", st[12], std::exp(-2.0) + (1 - std::exp(-2.0)) / 3, 1e-15);
  tk_eq(t, "cconc.cvb", st[14], std::sqrt(2.6) / 2, 1e-15);
  tk_eq(t, "cconc.nullb", st[15], 0.5, 1e-15);
  const double lvw[6] = {0, 0, 0, 0.1, 0.4, 1}, lpw[6] = {0, 0.2, 0.4, 0.6, 0.8, 1};
  tk_veq(t, "cconc.lorv", lv, lvw, 6, 1e-15);
  tk_veq(t, "cconc.lorp", lp, lpw, 6, 1e-15);
  const double tw[4] = {0.6, 1, 0, 1}, nw[4] = {1, 3, 0, 5};
  tk_veq(t, "cconc.top", top, tw, 4, 1e-15);
  tk_veq(t, "cconc.topn", tn, nw, 4, 0);
  // the same units as weights: 0 twice, then 1, 3, 6
  const double xw[4] = {0, 1, 3, 6}, ww[4] = {2, 1, 1, 1};
  double s2[16], t2[4], n2[4];
  clab_cconc(xw, ww, 4, fr, 4, s2, lp, lv, t2, n2, o);
  tk_eq(t, "cconc.wgini", s2[10], 0.6, 1e-15);
  tk_eq(t, "cconc.wvar", s2[3], 6.5, 1e-15);
  tk_veq(t, "cconc.wtop", t2, tw, 4, 1e-15);
  // fractional weights: 1 with weight 1/2, 2 with 3/2, Gini 3/28
  const double xf[2] = {2, 1}, wf[2] = {1.5, 0.5};
  clab_cconc(xf, wf, 2, fr, 0, s2, lp, lv, t2, n2, o);
  tk_eq(t, "cconc.fgini", s2[10], 3.0 / 28, 1e-15);
  // one unit: no variance
  clab_cconc(x, nullptr, 1, fr, 0, s2, lp, lv, t2, n2, o);
  tk_nan(t, "cconc.w1var", s2[3]);
  const double xn[2] = {-1, 2}, x0[2] = {0, 0}, w0[2] = {1, 0}, fb[1] = {1.5};
  tk_ieq(t, "cconc.n", clab_cconc(x, nullptr, 0, fr, 0, st, lp, lv, top, tn, o), CLAB_EARG);
  tk_ieq(t, "cconc.neg", clab_cconc(xn, nullptr, 2, fr, 0, st, lp, lv, top, tn, o), CLAB_EDOM);
  tk_ieq(t, "cconc.w0", clab_cconc(xf, w0, 2, fr, 0, st, lp, lv, top, tn, o), CLAB_EDOM);
  tk_ieq(t, "cconc.fr", clab_cconc(xf, nullptr, 2, fb, 1, st, lp, lv, top, tn, o), CLAB_EDOM);
  tk_ieq(t, "cconc.zero", clab_cconc(x0, nullptr, 2, fr, 0, st, lp, lv, top, tn, o), CLAB_EDOM);
}

// E|X - Y| / (2 mu) for independent Poisson(mu) X and Y
double gpois(double mu) {
  double p[400], s = 0;
  for (int k = 0; k < 400; ++k) p[k] = std::exp(-mu + k * std::log(mu) - std::lgamma(k + 1.0));
  for (int j = 0; j < 400; ++j) {
    for (int k = 0; k < 400; ++k) s += std::fabs(j - k) * p[j] * p[k];
  }
  return s / (2 * mu);
}

void null(tk &t) {
  // one large sample: its Gini is the sample mean difference over twice
  // the sample mean, which sits within 1% of the population value
  static int wk[40000];
  double g[2];
  const double mus[2] = {1.0, 40.0};
  for (int k = 0; k < 2; ++k) {
    tk_ieq(t, "cnull.rc", clab_cnull(mus[k], 40000, 1, 7 + k, g, wk), CLAB_OK);
    const double want = gpois(mus[k]) * 39999.0 / 40000.0;
    tk_true(t, "cnull.gini", std::fabs(g[0] - want) < 0.01 * want);
  }
  // a mean so small that a draw has no events gives 0
  clab_cnull(1e-12, 5, 2, 1, g, wk);
  tk_eq(t, "cnull.none", g[1], 0, 0);
  // the same seed, the same draws
  double h[2];
  clab_cnull(3, 50, 2, 9, g, wk);
  clab_cnull(3, 50, 2, 9, h, wk);
  tk_veq(t, "cnull.seed", g, h, 2, 0);
  tk_ieq(t, "cnull.n", clab_cnull(3, 0, 1, 1, g, wk), CLAB_EARG);
  tk_ieq(t, "cnull.mu", clab_cnull(0, 5, 1, 1, g, wk), CLAB_EDOM);
  tk_ieq(t, "cnull.seedneg", clab_cnull(3, 5, 1, -1, g, wk), CLAB_EDOM);
}

void grow(tk &t) {
  const int code[6] = {0, 0, 1, 0, 2, 1};
  double d[6], ex[6], lo[6], hi[6], m, sl;
  int seen[3];
  tk_ieq(t, "cgrow.rc", clab_cgrow(code, 6, 3, d, &m, ex, lo, hi, &sl, seen), CLAB_OK);
  const double dw[6] = {1, 1, 2, 2, 3, 3};
  tk_veq(t, "cgrow.dist", d, dw, 6, 0);
  double f = 0;
  for (int j = 0; j < 6; ++j) f += m / (m + j);
  tk_eq(t, "cgrow.root", f, 3, 1e-13);
  tk_eq(t, "cgrow.exlast", ex[5], 3, 1e-13);
  tk_eq(t, "cgrow.lo", lo[2], m * std::log((m + 3) / m), 1e-13);
  tk_eq(t, "cgrow.hi", hi[2], 1 + m * std::log((m + 2) / m), 1e-13);
  // slope over i = 3..6 of log d on log i
  double mx = 0, my = 0, sxy = 0, sxx = 0;
  for (int i = 3; i <= 6; ++i) {
    mx += std::log(i) / 4;
    my += std::log(dw[i - 1]) / 4;
  }
  for (int i = 3; i <= 6; ++i) {
    sxy += (std::log(i) - mx) * (std::log(dw[i - 1]) - my);
    sxx += (std::log(i) - mx) * (std::log(i) - mx);
  }
  tk_eq(t, "cgrow.slope", sl, sxy / sxx, 1e-13);
  // all new: M infinite; all one place: M 0; three events: no slope
  const int nw[3] = {0, 1, 2}, one[3] = {0, 0, 0};
  clab_cgrow(nw, 3, 3, d, &m, ex, lo, hi, &sl, seen);
  tk_true(t, "cgrow.inf", std::isinf(m) && std::isnan(ex[0]) && std::isnan(sl));
  clab_cgrow(one, 3, 3, d, &m, ex, lo, hi, &sl, seen);
  tk_eq(t, "cgrow.zero", m, 0, 0);
  // a root beyond 1e8: 100000 events at 99999 places
  static int big[100000], bs[100000];
  static double bd[100000], be[100000], bl[100000], bh[100000];
  for (int i = 0; i < 100000; ++i) big[i] = i == 99999 ? 0 : i;
  clab_cgrow(big, 100000, 100000, bd, &m, be, bl, bh, &sl, bs);
  tk_gt(t, "cgrow.big", m, 1e8);
  tk_eq(t, "cgrow.bigex", be[99999], 99999, 1e-12);
  const int bad[2] = {0, 3};
  tk_ieq(t, "cgrow.n", clab_cgrow(code, 1, 3, d, &m, ex, lo, hi, &sl, seen), CLAB_EARG);
  tk_ieq(t, "cgrow.code", clab_cgrow(bad, 2, 3, d, &m, ex, lo, hi, &sl, seen), CLAB_EARG);
}

void scale(tk &t) {
  // x = 2 pop^0.8 exactly, and a unit with no events left out
  const double pop[5] = {10, 100, 1000, 10000, 50};
  double x[5], o[5];
  for (int i = 0; i < 4; ++i) x[i] = 2 * std::pow(pop[i], 0.8);
  x[4] = 0;
  tk_ieq(t, "cscale.rc", clab_urbsc(x, pop, 5, o), CLAB_OK);
  tk_eq(t, "cscale.beta", o[0], 0.8, 1e-14);
  tk_eq(t, "cscale.pref", o[1], 2, 1e-13);
  tk_true(t, "cscale.se", o[2] < 1e-7);
  tk_eq(t, "cscale.r2", o[3], 1, 1e-14);
  tk_eq(t, "cscale.used", o[4], 4, 0);
  const double p1[4] = {5, 5, 5, 5}, xc[4] = {1, 1, 1, 1}, xn[1] = {-1};
  clab_urbsc(xc, pop, 4, o);
  tk_nan(t, "cscale.r2nan", o[3]);
  tk_ieq(t, "cscale.few", clab_urbsc(x + 2, pop + 2, 3, o), CLAB_EDOM);
  tk_ieq(t, "cscale.flat", clab_urbsc(xc, p1, 4, o), CLAB_EDOM);
  tk_ieq(t, "cscale.neg", clab_urbsc(xn, pop, 1, o), CLAB_EDOM);
  tk_ieq(t, "cscale.n", clab_urbsc(x, pop, 0, o), CLAB_EARG);
}

void zeta(tk &t) {
  double z;
  clab_hzta(2, 1, &z);
  tk_eq(t, "hzta.2", z, pi * pi / 6, 1e-15);
  clab_hzta(4, 1, &z);
  tk_eq(t, "hzta.4", z, pi * pi * pi * pi / 90, 1e-15);
  clab_hzta(2, 0.5, &z);
  tk_eq(t, "hzta.half", z, pi * pi / 2, 1e-15);
  clab_hzta(3, 1, &z);
  tk_eq(t, "hzta.3", z, 1.2020569031595942854, 1e-15);
  // zeta(s, q) - zeta(s, q + 1) = q^-s, at a large s and a large q
  double a, b;
  clab_hzta(25, 1.5, &a);
  clab_hzta(25, 2.5, &b);
  tk_eq(t, "hzta.shift", a - b, std::pow(1.5, -25.0), 1e-13);
  clab_hzta(1.7, 300, &a);
  clab_hzta(1.7, 301, &b);
  tk_eq(t, "hzta.shiftq", a - b, std::pow(300.0, -1.7), 1e-12);
  // next to the pole: 1 / (s - 1) + Euler's gamma + gamma_1 (s - 1)
  // (Stieltjes constants gamma_1 and gamma_2; s - 1 is exact in doubles)
  const double e = 1.0001 - 1;
  clab_hzta(1 + e, 1, &z);
  tk_eq(t, "hzta.pole", z,
        1 / e + 0.5772156649015329 + 0.0728158454836767 * e - 0.0048451815964361 * e * e, 1e-14);
  tk_ieq(t, "hzta.s", clab_hzta(1, 1, &z), CLAB_EDOM);
  tk_ieq(t, "hzta.q", clab_hzta(2, 0, &z), CLAB_EDOM);
}

void hill(tk &t) {
  const double x[5] = {8, 1, 4, 2, 0.5};
  double o[5], w[8];
  // continuous Hill on 1, 2, 4, 8 above xmin = 1
  tk_ieq(t, "hilltl.rc", clab_hilltl(x, 5, 1, 0, 0, o, w), CLAB_OK);
  const double a = 1 + 4 / (6 * std::log(2.0));
  tk_eq(t, "hilltl.cont", o[0], a, 1e-15);
  tk_eq(t, "hilltl.se", o[1], (a - 1) / 2, 1e-15);
  tk_eq(t, "hilltl.ntl", o[3], 4, 0);
  // the KS distance against F(u) = 1 - u^(1 - a), both sides of each jump
  double d = 0;
  for (int i = 0; i < 4; ++i) {
    const double f = 1 - std::pow(std::pow(2.0, i), 1 - a);
    d = std::fmax(d, std::fmax(std::fabs((i + 1) / 4.0 - f), std::fabs(i / 4.0 - f)));
  }
  tk_eq(t, "hilltl.ks", o[2], d, 1e-15);
  // the discrete closed form: xmin - 1/2
  clab_hilltl(x, 5, 1, 1, 1, o, w);
  tk_eq(t, "hilltl.disc", o[0],
        1 + 4 / (std::log(2.0) + std::log(4.0) + std::log(8.0) + std::log(16.0)), 1e-15);
  // the discrete MLE: the score is 0 at the maximum
  const double c[8] = {1, 1, 1, 2, 2, 3, 5, 9};
  clab_hilltl(c, 8, 1, 1, 0, o, w);
  double ls = 0;
  for (double v : c) ls += std::log(v);
  auto nll = [&](double al) {
    double z;
    clab_hzta(al, 1, &z);
    return 8 * std::log(z) + al * ls;
  };
  const double h = 1e-5;
  tk_true(t, "hilltl.score", std::fabs(nll(o[0] + h) - nll(o[0] - h)) / (2 * h) < 1e-5);
  tk_true(t, "hilltl.mle", nll(o[0]) <= nll(o[0] + 1e-3) && nll(o[0]) <= nll(o[0] - 1e-3));
  tk_eq(t, "hilltl.bound", o[4], 0, 0);
  tk_true(t, "hilltl.mlese", o[1] > 0 && o[2] > 0 && o[2] < 1);
  // 999 values at the threshold and one above: the maximum is at 25
  static double y[1000], yw[1000];
  for (int i = 0; i < 1000; ++i) y[i] = i ? 1000 : 1001;
  clab_hilltl(y, 1000, 1000, 1, 0, o, yw);
  tk_eq(t, "hilltl.atend", o[4], 1, 0);
  const double eq[3] = {2, 2, 2}, neg[2] = {1, -1};
  tk_ieq(t, "hilltl.n", clab_hilltl(x, 0, 1, 0, 0, o, w), CLAB_EARG);
  tk_ieq(t, "hilltl.xmin", clab_hilltl(x, 5, 0, 0, 0, o, w), CLAB_EDOM);
  tk_ieq(t, "hilltl.pos", clab_hilltl(neg, 2, 1, 0, 0, o, w), CLAB_EDOM);
  tk_ieq(t, "hilltl.few", clab_hilltl(x, 5, 8, 0, 0, o, w), CLAB_EDOM);
  tk_ieq(t, "hilltl.equal", clab_hilltl(eq, 3, 1, 0, 0, o, w), CLAB_EDOM);
}

}  // namespace

void t_cnc(tk &t, int big) {
  (void)big;
  conc(t);
  null(t);
  grow(t);
  scale(t);
  zeta(t);
  hill(t);
}
