// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for the Hawkes kernels (src/hwk.cpp) and the numerical
// helpers (src/qn.cpp): the log-likelihood by hand on three events, the
// gradient against central differences, the Weibull kernel of shape 1
// against the exponential, recovery of known parameters from simulation,
// the residuals of a Poisson process, and the minimiser and matrix
// inverse on cases with known answers. Part of the self-test
// (selftest.cpp).
#include <cmath>

#include "crimlab.h"
#include "qn.h"
#include "tchk.h"

namespace {

double rosen(const double *u, double *g, void *) {
  const double a = 1 - u[0], b = u[1] - u[0] * u[0];
  g[0] = -2 * a - 400 * u[0] * b;
  g[1] = 200 * b;
  return a * a + 100 * b * b;
}

// finite only at 0: no step can be taken
double spike(const double *u, double *g, void *) {
  g[0] = 1;
  return u[0] == 0 ? 0 : std::nan("");
}

// a saddle: its Hessian is not positive definite
double saddle(const double *u, double *g, void *) {
  g[0] = 2 * u[0];
  g[1] = -2 * u[1];
  return u[0] * u[0] - u[1] * u[1];
}

double nanf(const double *, double *g, void *) {
  g[0] = 0;
  return std::nan("");
}

void helpers(tk &t) {
  static double w[3 * 199 * 199];
  double u[2] = {-1.2, 1}, fv;
  int it;
  tk_ieq(t, "qn.rosen.conv", qn::min(2, u, rosen, nullptr, 500, 1e-10, &fv, &it, w, 1.0), 1);
  tk_eq(t, "qn.rosen.x", u[0], 1, 1e-6);
  tk_eq(t, "qn.rosen.y", u[1], 1, 1e-6);
  double s[1] = {0};
  tk_ieq(t, "qn.spike", qn::min(1, s, spike, nullptr, 10, 1e-10, &fv, &it, w, 1.0), 0);
  tk_ieq(t, "qn.nan", qn::min(1, s, nanf, nullptr, 10, 1e-10, &fv, &it, w, 1.0), -1);
  double q[2] = {0, 0};
  tk_ieq(t, "qn.steps", qn::min(2, q, rosen, nullptr, 1, 1e-10, &fv, &it, w, 1.0), 0);
  // the inverse of a 3 x 3 positive definite matrix: A A^-1 = I
  const double a[9] = {4, 2, 0.5, 2, 3, 1, 0.5, 1, 2};
  double b[9];
  for (int i = 0; i < 9; ++i) b[i] = a[i];
  tk_ieq(t, "qn.inv", qn::spdinv(b, 3), 0);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      double c = 0;
      for (int k = 0; k < 3; ++k) c += a[i + 3 * k] * b[k + 3 * j];
      tk_eq(t, "qn.inv.id", c, i == j ? 1 : 0, 1e-12);
    }
  }
  // the covariance of a quadratic: the inverse of its Hessian
  double th[2] = {0.3, -0.2}, cv[4], cw[6];
  tk_ieq(t, "qn.hcov.saddle", qn::hcov(2, th, saddle, nullptr, cv, cw), 1);
  tk_nan(t, "qn.hcov.nan", cv[0]);
  // Kolmogorov-Smirnov: n = 1 has P(D < d) = 2 d - 1 for d in [1/2, 1]
  tk_eq(t, "ks.n1", qn::kspv(0.8, 1, true, w), 0.4, 1e-14);
  // n = 2, d in [1/2, 1]: P(D >= d) = 2 (1 - d)^2
  tk_eq(t, "ks.n2", qn::kspv(0.75, 2, true, w), 0.125, 1e-14);
  // the limiting distribution on both sides of 1: 1 - K(x), the values of
  // R 4.6.1's stats:::pkolmogorov_two_asymp(x / 10, 100, lower.tail = FALSE)
  tk_eq(t, "ks.lim.hi", qn::kspv(1.5 / 10, 100, true, w), 0.022217962616525123, 1e-13);
  tk_eq(t, "ks.lim.lo", qn::kspv(0.8 / 10, 100, false, w), 0.544142503082887474, 1e-13);
  double us[4] = {0.9, 0.1, 0.6, 0.35};
  tk_eq(t, "ks.d", qn::ksd(us, 4), 0.15, 1e-15);
}

// three events by hand: lambda at each and the compensator
void byhand(tk &t) {
  const double tm[3] = {1, 2.5, 3}, hor = 4;
  const double mu = 0.4, eta = 0.6, b = 1.3;
  const double th[3] = {std::log(mu), eta, b};
  double ll;
  tk_ieq(t, "hwk.ll.rc", clab_hwkll(tm, 3, hor, 0, 0, 1, th, &ll, nullptr), CLAB_OK);
  const double l2 = mu + eta * b * std::exp(-1.5 * b),
               l3 = mu + eta * b * (std::exp(-2 * b) + std::exp(-0.5 * b));
  const double comp = mu * hor + eta * ((1 - std::exp(-3 * b)) + (1 - std::exp(-1.5 * b)) +
                                        (1 - std::exp(-b)));
  const double want = std::log(mu) + std::log(l2) + std::log(l3) - comp;
  tk_eq(t, "hwk.ll.exp", ll, want, 1e-14);
  // Weibull of shape 1 and scale 1 / b is the exponential kernel
  const double tw[4] = {std::log(mu), eta, 1, 1 / b};
  double lw;
  tk_ieq(t, "hwk.ll.wrc", clab_hwkll(tm, 3, hor, 1, 0, 1, tw, &lw, nullptr), CLAB_OK);
  tk_eq(t, "hwk.ll.weib1", lw, want, 1e-13);
  // Lomax by hand
  const double a = 2.5, c = 0.7, tl[4] = {std::log(mu), eta, a, c};
  auto g = [&](double d) { return a / c * std::pow(1 + d / c, -a - 1); };
  auto gg = [&](double d) { return 1 - std::pow(1 + d / c, -a); };
  const double m2 = mu + eta * g(1.5), m3 = mu + eta * (g(2) + g(0.5));
  const double cl = mu * hor + eta * (gg(3) + gg(1.5) + gg(1));
  tk_ieq(t, "hwk.ll.lrc", clab_hwkll(tm, 3, hor, 2, 0, 1, tl, &lw, nullptr), CLAB_OK);
  tk_eq(t, "hwk.ll.lomax", lw, std::log(mu) + std::log(m2) + std::log(m3) - cl, 1e-13);
  // a tie: the second event at 1 sees only the background
  const double tt[3] = {1, 1, 3};
  tk_ieq(t, "hwk.tie.rc", clab_hwkll(tt, 3, hor, 0, 0, 1, th, &ll, nullptr), CLAB_OK);
  const double lt = mu + 2 * eta * b * std::exp(-2 * b);
  const double ct = mu * hor + eta * (2 * (1 - std::exp(-3 * b)) + (1 - std::exp(-b)));
  tk_eq(t, "hwk.tie", ll, 2 * std::log(mu) + std::log(lt) - ct, 1e-14);
  // an event at the end of the window adds no kernel mass
  const double te[2] = {1, 4};
  tk_ieq(t, "hwk.end.rc", clab_hwkll(te, 2, hor, 2, 0, 1, tl, &ll, nullptr), CLAB_OK);
  tk_eq(t, "hwk.end", ll, std::log(mu) + std::log(mu + eta * g(3)) - mu * hor - eta * gg(3),
        1e-14);
  // the refusals
  const double bad[3] = {2, 1, 3};
  tk_ieq(t, "hwk.ll.order", clab_hwkll(bad, 3, hor, 0, 0, 1, th, &ll, nullptr), CLAB_EARG);
  const double neg[3] = {std::log(mu), -0.1, b};
  tk_ieq(t, "hwk.ll.eta", clab_hwkll(tm, 3, hor, 0, 0, 1, neg, &ll, nullptr), CLAB_EDOM);
  const double inf[6] = {0, HUGE_VAL, 0, 0, 0.5, 1};
  tk_ieq(t, "hwk.ll.inf", clab_hwkll(tm, 3, hor, 0, 1, 1, inf, &ll, nullptr), CLAB_EDOM);
  const double kz[3] = {0, 0.5, 0};
  tk_ieq(t, "hwk.ll.k", clab_hwkll(tm, 3, hor, 0, 0, 1, kz, &ll, nullptr), CLAB_EDOM);
  tk_ieq(t, "hwk.ll.base", clab_hwkll(tm, 3, hor, 0, 2, 1, th, &ll, nullptr), CLAB_EARG);
}

// the analytic gradient against central differences, every kernel and
// baseline, on 40 events
void grad(tk &t) {
  rng r{42};
  double tm[40];
  for (int i = 0; i < 40; ++i) tm[i] = runi(r, 0, 50);
  for (int i = 1; i < 40; ++i) {
    for (int j = i; j > 0 && tm[j] < tm[j - 1]; --j) {
      const double s = tm[j];
      tm[j] = tm[j - 1];
      tm[j - 1] = s;
    }
  }
  for (int kind = 0; kind < 3; ++kind) {
    for (int base = 0; base < 2; ++base) {
      double th[7] = {std::log(0.3)}, g[7];
      int np = 1;
      if (base) {
        th[1] = 0.4;
        th[2] = -0.3;
        th[3] = 0.2;
        np = 4;
      }
      th[np++] = 0.45;
      th[np++] = kind == 0 ? 0.8 : (kind == 1 ? 0.7 : 2.2);
      if (kind) th[np++] = kind == 1 ? 1.9 : 1.4;
      double ll;
      tk_ieq(t, "hwk.grad.rc", clab_hwkll(tm, 40, 55, kind, base, 13, th, &ll, g), CLAB_OK);
      for (int j = 0; j < np; ++j) {
        double hi = 0, lo = 0;
        const double h = 1e-6 * std::fmax(1, std::fabs(th[j])), s = th[j];
        th[j] = s + h;
        clab_hwkll(tm, 40, 55, kind, base, 13, th, &hi, nullptr);
        th[j] = s - h;
        clab_hwkll(tm, 40, 55, kind, base, 13, th, &lo, nullptr);
        th[j] = s;
        tk_eq(t, "hwk.grad", g[j], (hi - lo) / (2 * h), 1e-6);
      }
    }
  }
}

// simulate, fit, and find the parameters again; then the residuals
void recover(tk &t) {
  static double tm[20000], x[1], y[1], w[3 * 20000];
  static int gen[20000], iw[2 * 20000];
  const double th[3] = {std::log(0.8), 0.5, 1.5};
  int n = 0;
  tk_ieq(t, "hwk.sim.rc",
         clab_hwksim(4000, 0, 0, 1, th, nullptr, 0, 7, 20000, tm, x, y, gen, &n, w, iw), CLAB_OK);
  // the expected count: mu hor / (1 - eta), less the offspring cut at hor
  tk_near(t, "hwk.sim.n", n, 0.8 * 4000 / 0.5, 0.05);
  double fit[3], cov[9], st[3];
  tk_ieq(t, "hwk.fit.rc", clab_hwkfit(tm, n, 4000, 0, 0, 1, nullptr, 500, fit, cov, st), CLAB_OK);
  tk_ieq(t, "hwk.fit.conv", static_cast<int>(st[2]), 1);
  for (int j = 0; j < 3; ++j) tk_lt(t, "hwk.fit.z", std::fabs(fit[j] - th[j]), 4 * std::sqrt(cov[4 * j]));
  // the gradient vanishes there
  double g[3], ll;
  clab_hwkll(tm, n, 4000, 0, 0, 1, fit, &ll, g);
  tk_eq(t, "hwk.fit.ll", st[0], ll, 1e-12);
  for (int j = 0; j < 3; ++j) tk_lt(t, "hwk.fit.g", std::fabs(g[j]), 1e-3);
  // a Lomax fit from the same events, from its own starts, and a given one
  double fl[4], cl[16];
  tk_ieq(t, "hwk.fitl.rc", clab_hwkfit(tm, 600, tm[599], 2, 0, 1, nullptr, 300, fl, cl, st),
         CLAB_OK);
  const double s0[4] = {0, 0.3, 3, 1};
  tk_ieq(t, "hwk.fits.rc", clab_hwkfit(tm, 600, tm[599], 1, 0, 1, s0, 300, fl, cl, st), CLAB_OK);
  const double sb[4] = {0, 1.3, 3, 1};
  tk_ieq(t, "hwk.fits.bad", clab_hwkfit(tm, 600, tm[599], 1, 0, 1, sb, 300, fl, cl, st),
         CLAB_EDOM);
  tk_ieq(t, "hwk.fit.n", clab_hwkfit(tm, 1, 4000, 0, 0, 1, nullptr, 300, fl, cl, st), CLAB_EARG);
  // with eta = 0 the gaps rescale to mu times the gaps
  double tau[20000], ks[2];
  const double t0[3] = {std::log(0.8), 0, 1.5};
  tk_ieq(t, "hwk.res.rc", clab_hwkres(tm, 50, 4000, 0, 0, 1, t0, tau, ks, w), CLAB_OK);
  tk_eq(t, "hwk.res.0", tau[0], 0.8 * tm[0], 1e-13);
  tk_eq(t, "hwk.res.1", tau[7], 0.8 * (tm[7] - tm[6]), 1e-12);
  tk_ieq(t, "hwk.res.bad", clab_hwkres(tm, 50, 1, 0, 0, 1, t0, tau, ks, w), CLAB_EARG);
  const double tb[3] = {0, -1, 1};
  tk_ieq(t, "hwk.res.dom", clab_hwkres(tm, 50, 4000, 0, 0, 1, tb, tau, ks, w), CLAB_EDOM);
  // at the truth the residuals are uniform: the test does not reject
  tk_ieq(t, "hwk.res.true", clab_hwkres(tm, n, 4000, 0, 0, 1, th, tau, ks, w), CLAB_OK);
  tk_gt(t, "hwk.res.p", ks[1], 0.01);
  double mean = 0;
  for (int i = 0; i < n; ++i) mean += tau[i] / n;
  tk_near(t, "hwk.res.mean", mean, 1, 0.05);
  // the cycle baseline: the integral by Simpson's rule against a fine
  // trapezoid, through the residual of eta = 0 over one gap
  const double tc[6] = {0.2, 0.1, 0.4, -0.3, 0, 1};
  const double one[1] = {37};
  tk_ieq(t, "hwk.res.cyc", clab_hwkres(one, 1, 40, 0, 1, 10, tc, tau, ks, w), CLAB_OK);
  double s = 0;
  const int m = 200000;
  for (int i = 0; i <= m; ++i) {
    const double u = 37.0 * i / m;
    const double v = std::exp(0.2 + 0.1 * u / 40 + 0.4 * std::sin(2 * 3.141592653589793 * u / 10) -
                              0.3 * std::cos(2 * 3.141592653589793 * u / 10));
    s += (i == 0 || i == m ? 0.5 : 1) * v * 37.0 / m;
  }
  tk_eq(t, "hwk.res.simpson", tau[0], s, 1e-7);
  // ties in a Weibull fit's residuals: the asymptotic p-value
  const double tw[5] = {0.5, 1, 1, 2, 3}, pw[4] = {0, 0.3, 1.2, 0.8};
  tk_ieq(t, "hwk.res.tie", clab_hwkres(tw, 5, 4, 1, 0, 1, pw, tau, ks, w), CLAB_OK);
  tk_gt(t, "hwk.res.tiep", ks[1], 0);
}

void sim(tk &t) {
  static double tm[5000], x[5000], y[5000], w[3 * 5000];
  static int gen[5000], iw[2 * 5000];
  int n = 0;
  // space and time: every event in the window, offspring after parents
  const double th[3] = {std::log(2.0), 0.6, 1}, win[4] = {0, 10, 0, 5};
  tk_ieq(t, "hwk.sim.st",
         clab_hwksim(100, 0, 0, 1, th, win, 0.3, 3, 5000, tm, x, y, gen, &n, w, iw), CLAB_OK);
  bool in = true, ord = true;
  int off = 0;
  for (int i = 0; i < n; ++i) {
    in = in && x[i] >= 0 && x[i] <= 10 && y[i] >= 0 && y[i] <= 5;
    ord = ord && (i == 0 || tm[i] >= tm[i - 1]);
    off += gen[i] > 0;
  }
  tk_true(t, "hwk.sim.in", in && ord);
  tk_gt(t, "hwk.sim.off", off, 0);
  // the same seed, the same events
  int n2 = 0;
  static double t2[5000];
  clab_hwksim(100, 0, 0, 1, th, win, 0.3, 3, 5000, t2, x, y, gen, &n2, w, iw);
  tk_ieq(t, "hwk.sim.seed", n2, n);
  tk_eq(t, "hwk.sim.same", t2[n - 1], tm[n - 1], 0);
  // Weibull and Lomax lags, the cycle baseline
  const double tw[7] = {0, 0.1, 0.5, 0.5, 0.5, 1.5, 2};
  tk_ieq(t, "hwk.sim.w", clab_hwksim(300, 1, 1, 50, tw, nullptr, 0, 9, 5000, tm, x, y, gen, &n, w, iw),
         CLAB_OK);
  tk_gt(t, "hwk.sim.wn", n, 100);
  const double tl[4] = {0, 0.5, 2, 1};
  tk_ieq(t, "hwk.sim.l", clab_hwksim(300, 2, 0, 1, tl, nullptr, 0, 9, 5000, tm, x, y, gen, &n, w, iw),
         CLAB_OK);
  // the mean count of a Lomax process: mu hor / (1 - eta), about
  tk_near(t, "hwk.sim.ln", n, 600, 0.2);
  // too many events
  tk_ieq(t, "hwk.sim.cap", clab_hwksim(300, 2, 0, 1, tl, nullptr, 0, 9, 10, tm, x, y, gen, &n, w, iw),
         CLAB_EDOM);
  const double big[3] = {std::log(50.0), 0.9, 1};
  tk_ieq(t, "hwk.sim.cap2",
         clab_hwksim(100, 0, 0, 1, big, win, 0.3, 3, 50, tm, x, y, gen, &n, w, iw), CLAB_EDOM);
  tk_ieq(t, "hwk.sim.arg", clab_hwksim(-1, 0, 0, 1, th, nullptr, 0, 9, 10, tm, x, y, gen, &n, w, iw),
         CLAB_EARG);
  const double bw[4] = {0, -1, 0, 1};
  tk_ieq(t, "hwk.sim.win", clab_hwksim(1, 0, 0, 1, th, bw, 1, 9, 10, tm, x, y, gen, &n, w, iw),
         CLAB_EDOM);
  const double bt[3] = {0, -1, 1};
  tk_ieq(t, "hwk.sim.par", clab_hwksim(1, 0, 0, 1, bt, nullptr, 0, 9, 10, tm, x, y, gen, &n, w, iw),
         CLAB_EDOM);
}

}  // namespace

void t_hwk(tk &t, int big) {
  (void)big;
  helpers(t);
  byhand(t);
  grad(t);
  recover(t);
  sim(t);
}
