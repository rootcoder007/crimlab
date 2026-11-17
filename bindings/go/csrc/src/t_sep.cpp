// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/sep.cpp: hand values, identities of the area structures,
// conservation in the SIR model, score equations of the fits, and the
// refusals.
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "tchk.h"

namespace {

constexpr double pi = 3.14159265358979323846;

// the CSR rows of a path of n areas
void path(int n, std::vector<int> &rp, std::vector<int> &ci, std::vector<double> &w) {
  rp.assign(1, 0);
  ci.clear();
  w.clear();
  for (int i = 0; i < n; ++i) {
    if (i > 0) ci.push_back(i - 1), w.push_back(1);
    if (i < n - 1) ci.push_back(i + 1), w.push_back(1);
    rp.push_back(static_cast<int>(ci.size()));
  }
}

// P(X <= y) of the negative binomial (size s, probability p)
double pnb(double y, double s, double p) {
  double f = std::pow(p, s), c = f;
  for (int k = 0; k < y; ++k) {
    f *= (k + s) / (k + 1) * (1 - p);
    c += f;
  }
  return c;
}

}  // namespace

void t_sep(tk &T, int big) {
  (void)big;
  // the Bayes detector at the last of eight weeks
  const double x[8] = {3, 5, 2, 4, 6, 3, 4, 12};
  const int at[2] = {8, 7};
  double ub[2];
  int al[2];
  tk_ieq(T, "bo.ok", clab_bayout(x, 8, 52, 0, 6, 1, 0.05, at, 1, ub, al), CLAB_OK);
  tk_true(T, "bo.q", pnb(ub[0], 24.5, 6.0 / 7) >= 0.95 && pnb(ub[0] - 1, 24.5, 6.0 / 7) < 0.95);
  tk_ieq(T, "bo.alarm", al[0], x[7] > ub[0]);
  const double xn[8] = {3, NAN, 2, 4, 6, 3, 4, NAN};
  tk_ieq(T, "bo.na", clab_bayout(xn, 8, 52, 0, 6, 1, 0.05, at, 1, ub, al), CLAB_OK);
  tk_ieq(T, "bo.naal", al[0], -1);
  tk_ieq(T, "bo.past", clab_bayout(x, 8, 2, 1, 1, 0, 0.05, at, 2, ub, al), CLAB_OK);
  tk_ieq(T, "bo.early", clab_bayout(x, 8, 52, 1, 6, 1, 0.05, at, 1, ub, al), CLAB_EDOM);
  const double xa[8] = {NAN, NAN, NAN, NAN, NAN, NAN, NAN, 1};
  tk_ieq(T, "bo.empty", clab_bayout(xa, 8, 52, 0, 6, 1, 0.05, at, 1, ub, al), CLAB_EDOM);
  tk_ieq(T, "bo.alpha", clab_bayout(x, 8, 52, 0, 6, 1, 1, at, 1, ub, al), CLAB_EARG);
  // a window summing to 4833: p^s is a subnormal double; R's
  // qnbinom(0.95, 4833.5, 6 / 7) is 856
  const double xb[7] = {4833, 0, 0, 0, 0, 0, 0};
  const int a7 = 7;
  tk_ieq(T, "bo.sub", clab_bayout(xb, 7, 52, 0, 6, 1, 0.05, &a7, 1, ub, al), CLAB_OK);
  tk_eq(T, "bo.subq", ub[0], 856, 0);
  // buffers: two sources, one receptor near both
  const double rx[2] = {0, 5}, ry[2] = {0, 5}, sx[2] = {0.5, 1}, sy[2] = {0.5, 0}, sw[2] = {2, 3};
  double cn[2], ws[2], nr[2];
  tk_ieq(T, "be.ok", clab_bufexp(rx, ry, 2, sx, sy, 2, 1, sw, 0, cn, ws, nr), CLAB_OK);
  tk_eq(T, "be.cnt", cn[0], 2, 0);
  tk_eq(T, "be.w", ws[0], 5, 0);
  tk_eq(T, "be.cnt1", cn[1], 0, 0);
  tk_eq(T, "be.near", nr[0], std::sqrt(0.5), 1e-15);
  tk_ieq(T, "be.ll", clab_bufexp(rx, ry, 2, sx, sy, 2, 1e5, nullptr, 1, cn, ws, nr), CLAB_OK);
  tk_eq(T, "be.llw", ws[0], 1, 0);
  const double bad[1] = {100};
  tk_ieq(T, "be.lat", clab_bufexp(rx, bad, 1, sx, sy, 1, 1, nullptr, 1, cn, ws, nr), CLAB_EDOM);
  const double inf[1] = {INFINITY};
  tk_ieq(T, "be.inf", clab_bufexp(inf, ry, 1, sx, sy, 1, 1, nullptr, 0, cn, ws, nr), CLAB_EDOM);
  tk_ieq(T, "be.r", clab_bufexp(rx, ry, 2, sx, sy, 2, -1, nullptr, 0, cn, ws, nr), CLAB_EARG);
  // kernels by hand
  double ke[1];
  const double kx[1] = {0.2}, ky[1] = {0.3}, qx[2] = {0.25, 0.6}, qy[2] = {0.35, 0.8};
  tk_ieq(T, "ke.q", clab_kerexp(kx, ky, 1, qx, qy, 2, 0.5, nullptr, 1, ke), CLAB_OK);
  const double u = std::sqrt(0.005) / 0.5;
  tk_eq(T, "ke.qv", ke[0], 3 / pi * (1 - u * u) * (1 - u * u) / 0.25, 1e-14);
  tk_ieq(T, "ke.g", clab_kerexp(kx, ky, 1, qx, qy, 1, 0.5, sw, 0, ke), CLAB_OK);
  tk_eq(T, "ke.gv", ke[0], 2 * std::exp(-u * u / 2) / (2 * pi) / 0.25, 1e-14);
  tk_ieq(T, "ke.h", clab_kerexp(kx, ky, 1, qx, qy, 1, 0, nullptr, 0, ke), CLAB_EDOM);
  tk_ieq(T, "ke.inf", clab_kerexp(inf, ky, 1, qx, qy, 1, 1, nullptr, 0, ke), CLAB_EDOM);
  tk_ieq(T, "ke.k", clab_kerexp(kx, ky, 1, qx, qy, 1, 1, nullptr, 2, ke), CLAB_EARG);
  // area structures on a path of four
  std::vector<int> rp, ci;
  std::vector<double> w;
  path(4, rp, ci, w);
  double q[16], g[16], cv[16], sc = 0;
  tk_ieq(T, "le.ok", clab_lerprc(rp.data(), ci.data(), w.data(), 4, 1, 2, q), CLAB_OK);
  tk_eq(T, "le.row", q[0] + q[4] + q[8] + q[12], 0, 1e-15);
  tk_eq(T, "le.d", q[5], 4, 0);
  tk_ieq(T, "le.rho", clab_lerprc(rp.data(), ci.data(), w.data(), 4, 2, 2, q), CLAB_EDOM);
  const int bci[2] = {7, 0};
  const int brp[5] = {0, 1, 2, 2, 2};
  tk_ieq(T, "le.ci", clab_lerprc(brp, bci, w.data(), 4, 1, 1, q), CLAB_EARG);
  tk_ieq(T, "bm.ok", clab_bym2(rp.data(), ci.data(), w.data(), 4, 2, 0.4, 1e-10, &sc, g, cv), CLAB_OK);
  // Q G Q = Q
  double Q[16] = {1, -1, 0, 0, -1, 2, -1, 0, 0, -1, 2, -1, 0, 0, -1, 1}, qg[16];
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) {
      double s = 0;
      for (int k = 0; k < 4; ++k) s += Q[i + 4 * k] * g[k + 4 * j];
      qg[i + 4 * j] = s;
    }
  double dev = 0;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) {
      double s = 0;
      for (int k = 0; k < 4; ++k) s += qg[i + 4 * k] * Q[k + 4 * j];
      dev = std::fmax(dev, std::fabs(s - Q[i + 4 * j]));
    }
  tk_lt(T, "bm.qgq", dev, 1e-12);
  double gm = 0;
  for (int i = 0; i < 4; ++i) gm += std::log(g[5 * i] / sc);
  tk_eq(T, "bm.gm", gm, 0, 1e-12);
  tk_eq(T, "bm.cov", cv[0], (0.6 + 0.4 * g[0] / sc) / 2, 1e-15);
  const int irp[5] = {0, 1, 2, 2, 2}, ici[2] = {1, 0};
  const double iw[2] = {1, 1};
  tk_ieq(T, "bm.iso", clab_bym2(irp, ici, iw, 4, 1, 0.5, 1e-10, &sc, g, cv), CLAB_EDOM);
  tk_ieq(T, "bm.tau", clab_bym2(rp.data(), ci.data(), w.data(), 4, 0, 0.5, 1e-10, &sc, g, cv), CLAB_EDOM);
  tk_ieq(T, "bm.n", clab_bym2(rp.data(), ci.data(), w.data(), 1, 1, 0.5, 1e-10, &sc, g, cv), CLAB_EARG);
  const double nw[6] = {1, NAN, 1, 1, 1, 1};
  tk_ieq(T, "bm.nan", clab_bym2(rp.data(), ci.data(), nw, 4, 1, 0.5, 1e-10, &sc, g, cv), CLAB_EARG);
  // ecological regression without covariates: the Poisson intercept is the
  // log of the ratio of observed to expected
  const double y[6] = {2, 0, 5, 3, 0, 7}, e[6] = {2.5, 1, 4, 3, 1.5, 5};
  const double xc[6] = {0.1, -0.3, 0.5, 0.2, -0.6, 0.9};
  double cf[2], zc[2], th = 0, fit[6], ll = 0;
  tk_ieq(T, "er.ok", clab_ecolrg(y, e, 6, nullptr, 0, 0, nullptr, 0, 1e-12, 500, cf, zc, &th, fit, &ll), CLAB_OK);
  tk_eq(T, "er.b0", cf[0], std::log(17.0 / 17.0), 1e-12);
  tk_nan(T, "er.th", th);
  tk_ieq(T, "er.nb", clab_ecolrg(y, e, 6, xc, 1, 1, nullptr, 0, 1e-12, 500, cf, zc, &th, fit, &ll), CLAB_OK);
  tk_gt(T, "er.thp", th, 0);
  tk_ieq(T, "er.zip", clab_ecolrg(y, e, 6, xc, 1, 2, nullptr, 0, 1e-12, 50, cf, zc, &th, fit, &ll), CLAB_OK);
  tk_true(T, "er.zll", std::isfinite(ll));
  const double yeq[6] = {1, 1, 1, 1, 1, 1}, eeq[6] = {1, 1, 1, 1, 1, 1};
  tk_ieq(T, "er.nbp", clab_ecolrg(yeq, eeq, 6, nullptr, 0, 1, nullptr, 0, 1e-12, 5, cf, zc, &th, fit, &ll), CLAB_OK);
  const double xs[12] = {1, 2, 3, 4, 5, 6, 2, 4, 6, 8, 10, 12};
  tk_ieq(T, "er.sing", clab_ecolrg(y, e, 6, xs, 2, 0, nullptr, 0, 1e-12, 500, cf, zc, &th, fit, &ll), CLAB_EDOM);
  tk_ieq(T, "er.nbs", clab_ecolrg(y, e, 6, xs, 2, 1, nullptr, 0, 1e-12, 500, cf, zc, &th, fit, &ll), CLAB_EDOM);
  tk_ieq(T, "er.zs", clab_ecolrg(y, e, 6, nullptr, 0, 2, xs, 2, 1e-12, 5, cf, zc, &th, fit, &ll), CLAB_EDOM);
  const double eb[6] = {1, 0, 1, 1, 1, 1};
  tk_ieq(T, "er.e", clab_ecolrg(y, eb, 6, nullptr, 0, 0, nullptr, 0, 1e-12, 500, cf, zc, &th, fit, &ll), CLAB_EDOM);
  tk_ieq(T, "er.fam", clab_ecolrg(y, e, 6, nullptr, 0, 3, nullptr, 0, 1e-12, 500, cf, zc, &th, fit, &ll), CLAB_EARG);
  // the attention bias
  const int d[4] = {0, 1, -1, 5};
  const double tb[3] = {1, 0.5, 0.25};
  double bo[4];
  tk_ieq(T, "sb.ok", clab_spbias(d, 2, 2, tb, 3, -10, bo), CLAB_OK);
  tk_eq(T, "sb.0", bo[0], 1, 0);
  tk_eq(T, "sb.un", bo[2], -10, 0);
  tk_eq(T, "sb.far", bo[3], 0.25, 0);
  const int dn[1] = {-2};
  tk_ieq(T, "sb.neg", clab_spbias(dn, 1, 1, tb, 3, -10, bo), CLAB_EDOM);
  tk_ieq(T, "sb.nb", clab_spbias(d, 2, 2, tb, 0, -10, bo), CLAB_EARG);
  // SIR: each patch keeps its population; isolated patches stay free
  const double i0[4] = {5, 0, 0, 0}, pop[4] = {1000, 800, 1200, 500};
  double S[44], I[44], R[44];
  tk_ieq(T, "si.ok", clab_spsir(i0, 4, pop, rp.data(), ci.data(), w.data(), 0.5, 0.2, 0.1, 1, 10, S, I, R), CLAB_OK);
  tk_eq(T, "si.cons", S[10 + 11 * 2] + I[10 + 11 * 2] + R[10 + 11 * 2], 1200, 1e-12);
  tk_gt(T, "si.reach", R[10 + 11 * 3], 0);
  tk_ieq(T, "si.iso", clab_spsir(i0, 4, pop, nullptr, nullptr, nullptr, 0.5, 0.2, 0.1, 1, 10, S, I, R), CLAB_OK);
  tk_eq(T, "si.free", R[10 + 11 * 3], 0, 0);
  tk_ieq(T, "si.k", clab_spsir(i0, 4, pop, nullptr, nullptr, nullptr, 0.5, 0.2, 2, 1, 10, S, I, R), CLAB_EDOM);
  const double i1[1] = {2000};
  tk_ieq(T, "si.pop", clab_spsir(i1, 1, pop, nullptr, nullptr, nullptr, 0.5, 0.2, 0, 1, 1, S, I, R), CLAB_EDOM);
  tk_ieq(T, "si.ci", clab_spsir(i0, 4, pop, brp, bci, w.data(), 0.5, 0.2, 0.1, 1, 1, S, I, R), CLAB_EARG);
  tk_ieq(T, "si.n", clab_spsir(i0, 0, pop, nullptr, nullptr, nullptr, 0.5, 0.2, 0.1, 1, 1, S, I, R), CLAB_EARG);
  // the frailty and cure models solve their score equations
  const int n = 30;
  std::vector<double> t(n), ev(n), xx(n), t2(36), e2(36), x2(36);
  std::vector<int> rg(n), r2(36);
  for (int i = 0; i < n; ++i) {
    rg[i] = i % 3;
    xx[i] = std::sin(i);
    t[i] = 0.5 + ((i * 7) % 11) / 5.0 * std::exp(-0.5 * xx[i] - 0.3 * rg[i]);
    ev[i] = i % 5 != 0;
  }
  path(3, rp, ci, w);
  double fth[5], fse[5], ff = 0;
  int it = 0;
  tk_ieq(T, "fr.ok", clab_sfrail(t.data(), ev.data(), n, xx.data(), 1, rg.data(), rp.data(), ci.data(), w.data(), 3, 1,
                                 1e-10, fth, fse, &ff, &it),
         CLAB_OK);
  double sc0 = 0;
  for (int i = 0; i < n; ++i) {
    const double H = std::pow(t[i], fth[0]) * std::exp(xx[i] * fth[1] + fth[2 + rg[i]]);
    sc0 += (ev[i] - H) * xx[i];
  }
  tk_lt(T, "fr.score", std::fabs(sc0), 1e-7);
  tk_true(T, "fr.se", fse[1] > 0);
  for (int i = 0; i < 36; ++i) {
    x2[i] = std::sin(i);
    t2[i] = i % 4 == 0 ? 5 : 0.2 + ((i * 7) % 11) / 8.0;
    e2[i] = !(i % 4 == 0 || i % 9 == 1);
    r2[i] = i % 3;
  }
  double cth[6], cse[6], cf2 = 0, cure[36];
  tk_ieq(T, "cu.ok", clab_sfcure(t2.data(), e2.data(), 36, x2.data(), 1, r2.data(), rp.data(), ci.data(), w.data(), 3, 1,
                                 1e-10, cth, cse, &cf2, cure, &it),
         CLAB_OK);
  tk_eq(T, "cu.cure", cure[0], std::exp(-std::exp(x2[0] * cth[2] + cth[3])), 1e-12);
  double sc1 = 0;
  for (int i = 0; i < 36; ++i) {
    const double tv = std::exp(x2[i] * cth[2] + cth[3 + r2[i]]), G = std::exp(cth[1] + cth[0] * std::log(t2[i]));
    sc1 += (e2[i] - tv * (1 - std::exp(-G))) * x2[i];
  }
  tk_lt(T, "cu.score", std::fabs(sc1), 1e-7);
  std::vector<double> tb2(t);
  tb2[0] = -1;
  tk_ieq(T, "fr.t", clab_sfrail(tb2.data(), ev.data(), n, xx.data(), 1, rg.data(), rp.data(), ci.data(), w.data(), 3, 1,
                                1e-10, fth, fse, &ff, &it),
         CLAB_EDOM);
  tk_ieq(T, "cu.t", clab_sfcure(tb2.data(), ev.data(), n, xx.data(), 1, rg.data(), rp.data(), ci.data(), w.data(), 3, 1,
                                1e-10, cth, cse, &cf2, cure, &it),
         CLAB_EDOM);
  std::vector<double> xnan(xx);
  xnan[3] = NAN;
  tk_ieq(T, "fr.x", clab_sfrail(t.data(), ev.data(), n, xnan.data(), 1, rg.data(), rp.data(), ci.data(), w.data(), 3, 1,
                                1e-10, fth, fse, &ff, &it),
         CLAB_EDOM);
  tk_ieq(T, "fr.w", clab_sfrail(t.data(), ev.data(), n, xx.data(), 1, rg.data(), rp.data(), ci.data(), nw, 3, 1, 1e-10,
                                fth, fse, &ff, &it),
         CLAB_EDOM);
  // a constant covariate in one area: the Hessian is singular
  const int one[1] = {0};
  const int orp[2] = {0, 0};
  std::vector<int> r0(n, 0);
  std::vector<double> ones(n, 1.0);
  tk_ieq(T, "fr.sing", clab_sfrail(t.data(), ev.data(), n, ones.data(), 1, r0.data(), orp, one, iw, 1, 0, 1e-10, fth, fse, &ff,
                                   &it),
         CLAB_EDOM);
  tk_ieq(T, "cu.sing", clab_sfcure(t.data(), ev.data(), n, ones.data(), 1, r0.data(), orp, one, iw, 1, 0, 1e-10, cth, cse, &cf2,
                                   cure, &it),
         CLAB_EDOM);
  // eigenvector filtering on a path of eight
  path(8, rp, ci, w);
  const double v[8] = {1, 2, 2.5, 4, 3.5, 3, 1.5, 1};
  int sel[8], ns = 0;
  double co[16], st[4], res[8], mi[7];
  for (int cr = 0; cr < 4; ++cr)
    for (int ca = 0; ca < 3; ++ca) {
      tk_ieq(T, "ev.ok", clab_evfilt(v, 8, nullptr, 0, rp.data(), ci.data(), w.data(), cr, ca, cr == 3 ? 0.01 : 0.1, -1,
                                     sel, &ns, co, st, res, mi),
             CLAB_OK);
    }
  const double x8[8] = {0.1, -0.3, 0.5, 0.2, -0.6, 0.9, 0.4, -0.2};
  tk_ieq(T, "ev.max", clab_evfilt(v, 8, x8, 1, rp.data(), ci.data(), w.data(), 0, 0, 1e-9, 1, sel, &ns, co, st, res, mi),
         CLAB_OK);
  tk_ieq(T, "ev.ns", ns, 1);
  tk_gt(T, "ev.mi", mi[0], mi[6]);
  const int erp[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  tk_ieq(T, "ev.now", clab_evfilt(v, 8, nullptr, 0, erp, ci.data(), w.data(), 0, 0, 0.1, -1, sel, &ns, co, st, res, mi),
         CLAB_EDOM);
  const double vn[8] = {1, NAN, 1, 1, 1, 1, 1, 1};
  tk_ieq(T, "ev.nan", clab_evfilt(vn, 8, nullptr, 0, rp.data(), ci.data(), w.data(), 0, 0, 0.1, -1, sel, &ns, co, st, res,
                                  mi),
         CLAB_EDOM);
  tk_ieq(T, "ev.crit", clab_evfilt(v, 8, nullptr, 0, rp.data(), ci.data(), w.data(), 4, 0, 0.1, -1, sel, &ns, co, st, res,
                                   mi),
         CLAB_EARG);
  const int brp8[9] = {0, 1, 1, 1, 1, 1, 1, 1, 1}, bci8[1] = {9};
  tk_ieq(T, "ev.w", clab_evfilt(v, 8, nullptr, 0, brp8, bci8, w.data(), 0, 0, 0.1, -1, sel, &ns, co, st, res, mi),
         CLAB_EARG);
  const double xd[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  tk_ieq(T, "ev.sing", clab_evfilt(v, 8, xd, 2, rp.data(), ci.data(), w.data(), 1, 0, 0.1, -1, sel, &ns, co, st, res,
                                   mi),
         CLAB_EDOM);
}
