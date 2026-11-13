// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for the marked, multivariate and space-time Hawkes
// kernels (src/mhk.cpp): ETAS without marks is the Lomax Hawkes process
// of clab_hwkll, the multivariate process of one type the exponential
// one; the space-time log-likelihood and intensity by hand; every
// gradient against central differences; the spectral radius of a known
// matrix; fits that reach a vanishing gradient. Part of the self-test
// (selftest.cpp).
#include <cmath>

#include "crimlab.h"
#include "tchk.h"

namespace {

constexpr double pi = 3.141592653589793;

// sorted uniform times on [0, hor)
void unif(rng &r, double *t, int n, double hor) {
  for (int i = 0; i < n; ++i) t[i] = runi(r, 0, hor);
  for (int i = 1; i < n; ++i) {
    for (int j = i; j > 0 && t[j] < t[j - 1]; --j) {
      const double s = t[j];
      t[j] = t[j - 1];
      t[j - 1] = s;
    }
  }
}

void etas(tk &t) {
  rng r{7};
  double tm[60], dm[60], th[5], cov[25], st[4], g[5];
  unif(r, tm, 60, 100);
  for (int i = 0; i < 60; ++i) dm[i] = 0;
  // no marks: the Lomax process, shape a = p - 1, scale c, K = eta a c^a
  const double nu = 0.3, eta = 0.4, a = 1.7, c = 0.8;
  const double s0[5] = {nu, eta * a * std::pow(c, a), 0, c, a + 1}, h[4] = {std::log(nu), eta, a, c};
  tk_ieq(t, "etas.rc", clab_etasf(tm, dm, 60, 100, 0, s0, 0, th, cov, st, g), CLAB_OK);
  double ll;
  clab_hwkll(tm, 60, 100, 2, 0, 1, h, &ll, nullptr);
  tk_eq(t, "etas.lomax", st[0], ll, 1e-12);
  tk_eq(t, "etas.branch", st[3], eta, 1e-13);
  // the gradient with marks against central differences
  for (int i = 0; i < 60; ++i) dm[i] = runi(r, 0, 2);
  double p[5] = {0.3, 0.2, 0.6, 0.5, 1.4};
  clab_etasf(tm, dm, 60, 100, 1, p, 0, th, cov, st, g);
  for (int j = 0; j < 5; ++j) {
    const double s = p[j], e = 1e-6 * std::fmax(1, std::fabs(s));
    double hi, lo;
    p[j] = s + e;
    clab_etasf(tm, dm, 60, 100, 1, p, 0, th, cov, st, nullptr);
    hi = st[0];
    p[j] = s - e;
    clab_etasf(tm, dm, 60, 100, 1, p, 0, th, cov, st, nullptr);
    lo = st[0];
    p[j] = s;
    tk_eq(t, "etas.grad", g[j], (hi - lo) / (2 * e), 1e-6);
  }
  // fits: with and without marks, from the starts or a given one; the
  // marked fit on a clustered (Lomax) series
  static double ts[2000], wk[6000];
  static int gn[2000], iw[4000];
  static double dz[2000];
  int ns = 0;
  const double lp[4] = {std::log(0.3), 0.5, 0.6, 1};
  clab_hwksim(800, 2, 0, 1, lp, nullptr, 0, 5, 2000, ts, nullptr, nullptr, gn, &ns, wk, iw);
  for (int i = 0; i < ns; ++i) dz[i] = runi(r, 0, 1);
  tk_ieq(t, "etas.fit", clab_etasf(ts, dz, ns, 800, 1, nullptr, 500, th, cov, st, g), CLAB_OK);
  tk_ieq(t, "etas.fit.conv", static_cast<int>(st[2]), 1);
  // the gradient in the fit's own scale (log mu, log K, alpha, log c,
  // log(p - 1)) vanishes
  const double sc[5] = {th[0], th[1], 1, th[3], th[4] - 1};
  for (int j = 0; j < 5; ++j) tk_lt(t, "etas.fit.g", std::fabs(g[j] * sc[j]), 1e-3 * (1 + std::fabs(st[0])));
  tk_lt(t, "etas.fit.br", std::fabs(st[3] - 0.5), 0.15);
  tk_ieq(t, "etas.fit0", clab_etasf(tm, dm, 60, 100, 0, s0, 300, th, cov, st, g), CLAB_OK);
  tk_eq(t, "etas.fit0.a", th[2], 0, 0);
  tk_nan(t, "etas.fit0.cov", cov[2 + 5 * 2]);
  const double bad[5] = {0.3, 0.2, 0.6, 0.5, 0.9};
  tk_ieq(t, "etas.bad", clab_etasf(tm, dm, 60, 100, 1, bad, 10, th, cov, st, g), CLAB_EDOM);
  tk_ieq(t, "etas.bad0", clab_etasf(tm, dm, 60, 100, 1, bad, 0, th, cov, st, g), CLAB_EDOM);
  tk_ieq(t, "etas.arg", clab_etasf(tm, dm, 60, 50, 1, nullptr, 10, th, cov, st, g), CLAB_EARG);
}

void mv(tk &t) {
  rng r{8};
  static double tm[80], w[2 * 2 + 8 * 7 + 49];
  int ty[80], iw[2];
  double th[7], cov[49], st[4], g[7];
  unif(r, tm, 80, 60);
  // one type: the exponential Hawkes process
  for (int i = 0; i < 80; ++i) ty[i] = 0;
  const double s1[3] = {0.5, 0.35, 1.2}, h[3] = {std::log(0.5), 0.35, 1.2};
  tk_ieq(t, "mv.rc", clab_mvhwk(tm, ty, 80, 1, 60, 0, s1, 0, th, cov, st, g, w, iw), CLAB_OK);
  double ll;
  clab_hwkll(tm, 80, 60, 0, 0, 1, h, &ll, nullptr);
  tk_eq(t, "mv.one", st[0], ll, 1e-12);
  tk_eq(t, "mv.rad1", st[3], 0.35, 1e-12);
  // two types: the gradient, and the spectral radius of [0.2 0.3; 0.1 0.4]
  for (int i = 0; i < 80; ++i) ty[i] = i % 3 == 0;
  double p[7] = {0.4, 0.3, 0.2, 0.1, 0.3, 0.4, 0.9};
  clab_mvhwk(tm, ty, 80, 2, 60, 0, p, 0, th, cov, st, g, w, iw);
  tk_eq(t, "mv.rad", st[3], 0.5, 1e-12);
  for (int j = 0; j < 7; ++j) {
    const double s = p[j], e = 1e-6 * std::fmax(1, std::fabs(s));
    double hi, lo;
    p[j] = s + e;
    clab_mvhwk(tm, ty, 80, 2, 60, 0, p, 0, th, cov, st, nullptr, w, iw);
    hi = st[0];
    p[j] = s - e;
    clab_mvhwk(tm, ty, 80, 2, 60, 0, p, 0, th, cov, st, nullptr, w, iw);
    lo = st[0];
    p[j] = s;
    tk_eq(t, "mv.grad", g[j], (hi - lo) / (2 * e), 1e-6);
  }
  // a fit with beta given, and one from a start
  tk_ieq(t, "mv.fit", clab_mvhwk(tm, ty, 80, 2, 60, 1, nullptr, 300, th, cov, st, g, w, iw), CLAB_OK);
  for (int j = 0; j < 6; ++j) tk_lt(t, "mv.fit.g", std::fabs(g[j]) * th[j], 1e-2);
  tk_ieq(t, "mv.fit.s", clab_mvhwk(tm, ty, 80, 2, 60, 0, p, 50, th, cov, st, g, w, iw), CLAB_OK);
  const double bad[7] = {0.4, 0.3, 0.2, 0.1, 0.3, 0.4, 0};
  tk_ieq(t, "mv.bad", clab_mvhwk(tm, ty, 80, 2, 60, 0, bad, 0, th, cov, st, g, w, iw), CLAB_EDOM);
  tk_ieq(t, "mv.bads", clab_mvhwk(tm, ty, 80, 2, 60, 0, bad, 9, th, cov, st, g, w, iw), CLAB_EDOM);
  const double neg[7] = {0.4, 0.3, -0.2, 0.1, 0.3, 0.4, 1};
  tk_ieq(t, "mv.neg", clab_mvhwk(tm, ty, 80, 2, 60, 0, neg, 0, th, cov, st, g, w, iw), CLAB_EDOM);
  ty[3] = 2;
  tk_ieq(t, "mv.ty", clab_mvhwk(tm, ty, 80, 2, 60, 0, p, 0, th, cov, st, g, w, iw), CLAB_EARG);
}

void st(tk &t) {
  const double tm[2] = {1, 2}, x[2] = {0, 1}, y[2] = {0, 0};
  const double mu = 0.2, eta = 0.5, b = 1.3, s = 0.7, p[4] = {mu, eta, b, s};
  double th[4], cov[16], sv[3], g[4], ql[1];
  auto phi = [&](double r2) { return std::exp(-r2 / (2 * s * s)) / (2 * pi * s * s); };
  const double qt[1] = {2.5}, qx[1] = {0}, qy[1] = {0};
  tk_ieq(t, "st.rc", clab_sthwk(tm, x, y, 2, 3, 4, p, 0, th, cov, sv, g, qt, qx, qy, 1, ql), CLAB_OK);
  const double l2 = mu + eta * b * std::exp(-b) * phi(1);
  const double comp = mu * 12 + eta * (2 - std::exp(-2 * b) - std::exp(-b));
  tk_eq(t, "st.ll", sv[0], std::log(mu) + std::log(l2) - comp, 1e-14);
  tk_eq(t, "st.lam", ql[0], mu + eta * b * (std::exp(-1.5 * b) * phi(0) + std::exp(-0.5 * b) * phi(1)), 1e-14);
  // the gradient on 50 events
  rng r{9};
  double t5[50], x5[50], y5[50];
  unif(r, t5, 50, 40);
  for (int i = 0; i < 50; ++i) {
    x5[i] = runi(r, 0, 3);
    y5[i] = runi(r, 0, 2);
  }
  double q[4] = {0.1, 0.4, 0.9, 0.5};
  clab_sthwk(t5, x5, y5, 50, 40, 6, q, 0, th, cov, sv, g, nullptr, nullptr, nullptr, 0, nullptr);
  for (int j = 0; j < 4; ++j) {
    const double s0 = q[j], e = 1e-6 * std::fmax(1, std::fabs(s0));
    double hi, lo;
    q[j] = s0 + e;
    clab_sthwk(t5, x5, y5, 50, 40, 6, q, 0, th, cov, sv, nullptr, nullptr, nullptr, nullptr, 0, nullptr);
    hi = sv[0];
    q[j] = s0 - e;
    clab_sthwk(t5, x5, y5, 50, 40, 6, q, 0, th, cov, sv, nullptr, nullptr, nullptr, nullptr, 0, nullptr);
    lo = sv[0];
    q[j] = s0;
    tk_eq(t, "st.grad", g[j], (hi - lo) / (2 * e), 1e-6);
  }
  tk_ieq(t, "st.fit", clab_sthwk(t5, x5, y5, 50, 40, 6, nullptr, 300, th, cov, sv, g, nullptr, nullptr, nullptr, 0,
                                 nullptr),
         CLAB_OK);
  tk_ieq(t, "st.fits", clab_sthwk(t5, x5, y5, 50, 40, 6, q, 300, th, cov, sv, g, nullptr, nullptr, nullptr, 0,
                                  nullptr),
         CLAB_OK);
  const double bad[4] = {0.1, 1.4, 0.9, 0.5}, bz[4] = {0, 0.4, 0.9, 0.5};
  tk_ieq(t, "st.bad", clab_sthwk(t5, x5, y5, 50, 40, 6, bad, 30, th, cov, sv, g, nullptr, nullptr, nullptr, 0,
                                 nullptr),
         CLAB_EDOM);
  tk_ieq(t, "st.bz", clab_sthwk(t5, x5, y5, 50, 40, 6, bz, 0, th, cov, sv, g, nullptr, nullptr, nullptr, 0, nullptr),
         CLAB_EDOM);
  const double xn[2] = {0, HUGE_VAL};
  tk_ieq(t, "st.x", clab_sthwk(tm, xn, y, 2, 3, 4, p, 0, th, cov, sv, g, qt, qx, qy, 1, ql), CLAB_EDOM);
  tk_ieq(t, "st.arg", clab_sthwk(tm, x, y, 2, 3, 0, p, 0, th, cov, sv, g, qt, qx, qy, 1, ql), CLAB_EARG);
}

}  // namespace

void t_mhk(tk &t, int big) {
  (void)big;
  etas(t);
  mv(t);
  st(t);
}
