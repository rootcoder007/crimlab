// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/apx.cpp and of clab_ppmlfe and clab_fesep (cnt.cpp): closed
// forms on tiny panels (the exponential GMM just identified by a binary
// instrument among them), the absorbed Poisson fit against the same model
// with the effects as columns (clab_cntreg), and the refusals.
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "tchk.h"

namespace {

// clab_ppmlfe with the work arrays it asks for.
int ppml(const double *y, const double *x, int n, int p, const int *f1, int n1, const int *f2,
         int n2, const int *cl, int ncl, double *b, double *vc, double *st) {
  std::vector<double> w(static_cast<std::size_t>(n) * (2 * p + 8) + 3 * p * p + p + 4 * (n1 + n2));
  std::vector<int> iw(static_cast<std::size_t>(n) + 1);
  return clab_ppmlfe(y, x, n, p, nullptr, nullptr, f1, n1, f2, n2, cl, ncl, nullptr, 0, 100,
                     1e-12, 1e-13, b, vc, nullptr, nullptr, st, w.data(), iw.data());
}

}  // namespace

void t_apx(tk &T, int) {
  // two periods, x 0 then 1, unit effects: the estimate is the log of the
  // ratio of the totals
  {
    const double y[4] = {3, 5, 2, 7}, x[4] = {0, 1, 0, 1};
    const int u[4] = {0, 0, 1, 1}, c[4] = {0, 0, 1, 1};
    double b = 0, vc = 0, st[6];
    tk_ieq(T, "apx.pp2", ppml(y, x, 4, 1, u, 2, nullptr, 0, c, 2, &b, &vc, st), CLAB_OK);
    tk_eq(T, "apx.pp2b", b, std::log(12.0 / 5.0), 1e-10);
    tk_ieq(T, "apx.pp2c", static_cast<int>(st[1]), 1);
  }
  // two effects absorbed equal the effects as columns (clab_cntreg), the
  // coefficient and its cluster sandwich (FWL; the G / (G - 1) factor apart)
  {
    const int nu = 4, nt = 5, n = nu * nt;
    rng r{20261115};
    std::vector<double> y(n), x(n);
    std::vector<int> u(n), t(n);
    for (int i = 0; i < n; ++i) {
      u[i] = i / nt;
      t[i] = i % nt;
      x[i] = runi(r, -1, 1);
      y[i] = rpois(r, std::exp(0.5 + 0.4 * x[i] + 0.2 * u[i] - 0.1 * t[i]));
    }
    double b = 0, vc = 0, st[6];
    tk_ieq(T, "apx.pfe", ppml(y.data(), x.data(), n, 1, u.data(), nu, t.data(), nt, u.data(), nu,
                              &b, &vc, st), CLAB_OK);
    const int p = 1 + 1 + (nu - 1) + (nt - 1);
    std::vector<double> xd(static_cast<std::size_t>(n) * p, 0.0);
    for (int i = 0; i < n; ++i) {
      xd[i] = 1;
      xd[n + i] = x[i];
      if (u[i] > 0) xd[(1 + u[i]) * n + i] = 1;
      if (t[i] > 0) xd[(1 + nu - 1 + t[i]) * n + i] = 1;
    }
    std::vector<double> bd(p), vb(p * p), vd(p * p), sd(8);
    std::vector<double> wk(static_cast<std::size_t>(n) * p + 6 * n + p + p * p + 2 * nu * p);
    std::vector<int> al(p), iw(p);
    clab_cntreg(y.data(), xd.data(), n, p, nullptr, u.data(), nu, 0, 100, 1e-14, bd.data(),
                vb.data(), vd.data(), sd.data(), al.data(), wk.data(), iw.data());
    tk_eq(T, "apx.pfeb", b, bd[1], 1e-8);
    tk_eq(T, "apx.pfev", vc, vd[1 * p + 1] * nu / (nu - 1.0), 1e-6);
    // a unit of zeros is dropped and changes nothing
    std::vector<int> keep(n);
    std::vector<double> sc(nu + nt);
    std::vector<double> y0 = y;
    for (int i = 0; i < nt; ++i) y0[i] = 0;
    tk_ieq(T, "apx.sep", clab_fesep(y0.data(), n, u.data(), nu, t.data(), nt, keep.data(), sc.data()), nt);
    tk_ieq(T, "apx.sepk", keep[0] + keep[nt], 1);
    // a regressor the effects absorb is refused
    std::vector<double> xa(n);
    for (int i = 0; i < n; ++i) xa[i] = u[i];
    tk_ieq(T, "apx.abs", ppml(y.data(), xa.data(), n, 1, u.data(), nu, t.data(), nt, u.data(), nu,
                              &b, &vc, st), CLAB_EDOM);
    tk_ieq(T, "apx.cl1", ppml(y.data(), x.data(), n, 1, u.data(), nu, nullptr, 0, u.data(), 1, &b,
                              &vc, st), CLAB_EARG);
  }
  // just identified with one group (an intercept): cov(z, y) / cov(z, x)
  {
    const double y[5] = {1, 3, 2, 6, 5}, x[5] = {0, 2, 1, 4, 2}, z[5] = {0, 1, 1, 2, 3};
    const int g[5] = {0, 0, 0, 0, 0}, c[5] = {0, 0, 1, 1, 1};
    double my = 0, mx = 0, mz = 0;
    for (int i = 0; i < 5; ++i) {
      my += y[i] / 5;
      mx += x[i] / 5;
      mz += z[i] / 5;
    }
    double szy = 0, szx = 0;
    for (int i = 0; i < 5; ++i) {
      szy += (z[i] - mz) * (y[i] - my);
      szx += (z[i] - mz) * (x[i] - mx);
    }
    double b = 0, vc = 0, f[4] = {0, 0, 0, 0}, st[3];
    tk_ieq(T, "apx.iv", clab_ivfe(y, x, 1, nullptr, 0, z, 1, 5, nullptr, g, 1, nullptr, 0, c, 2,
                                  nullptr, 0, 1e-13, &b, &vc, f, nullptr, nullptr, nullptr, nullptr, st), CLAB_OK);
    tk_eq(T, "apx.ivb", b, szy / szx, 1e-12);
    // one instrument: the effective F is the Wald F
    tk_eq(T, "apx.ivef", f[2], f[0], 1e-10);
    tk_eq(T, "apx.ivkp", f[1], f[0], 0);
    // one instrument: Montiel Olea and Pflueger's effective degrees of
    // freedom are 1 (equation 4 with a scalar W2)
    tk_eq(T, "apx.keff", f[3], 1, 1e-12);
    // the score variance by clusters equals the kernel's own sandwich
    double sco[5], brd = 0, v2 = 0;
    clab_ivfe(y, x, 1, nullptr, 0, z, 1, 5, nullptr, g, 1, nullptr, 0, c, 2, nullptr, 0, 1e-13, &b,
              &vc, f, nullptr, nullptr, sco, &brd, st);
    tk_ieq(T, "apx.vcx", clab_vcx(sco, 5, 1, &brd, 0, c, 2, nullptr, 0, nullptr, 0, 2.0, &v2), CLAB_OK);
    tk_eq(T, "apx.vcxe", v2, vc, 1e-12);
    // rows alone equal clusters of one row each; Newey-West with lag 0 too
    const int ci[5] = {0, 1, 2, 3, 4};
    const double tt[5] = {1, 2, 3, 4, 5};
    double vh = 0, vr = 0, vn = 0, vn1 = 0;
    clab_vcx(sco, 5, 1, &brd, 2, nullptr, 0, nullptr, 0, nullptr, 0, 1.0, &vh);
    clab_vcx(sco, 5, 1, &brd, 0, ci, 5, nullptr, 0, nullptr, 0, 1.0, &vr);
    tk_eq(T, "apx.hc", vh, vr, 1e-14);
    clab_vcx(sco, 5, 1, &brd, 3, g, 1, nullptr, 0, tt, 0, 1.0, &vn);
    tk_eq(T, "apx.nw0", vn, vh, 1e-14);
    // lag 1 by hand: sum s_i^2 + (1/2) 2 sum s_i s_{i+1}
    clab_vcx(sco, 5, 1, &brd, 3, g, 1, nullptr, 0, tt, 1, 1.0, &vn1);
    double mm = 0;
    for (int i = 0; i < 5; ++i) mm += sco[i] * sco[i] + (i < 4 ? sco[i] * sco[i + 1] : 0.0);
    tk_eq(T, "apx.nw1", vn1, brd * mm * brd, 1e-12);
    const double zc[5] = {1, 1, 1, 1, 1};
    tk_ieq(T, "apx.ivz", clab_ivfe(y, x, 1, nullptr, 0, zc, 1, 5, nullptr, g, 1, nullptr, 0, c, 2,
                                   nullptr, 0, 1e-13, &b, &vc, f, nullptr, nullptr, nullptr, nullptr, st), CLAB_EDOM);
    tk_ieq(T, "apx.ivq", clab_ivfe(y, x, 1, nullptr, 0, z, 0, 5, nullptr, g, 1, nullptr, 0, c, 2,
                                   nullptr, 0, 1e-13, &b, &vc, f, nullptr, nullptr, nullptr, nullptr, st), CLAB_EARG);
  }
  // the exponential GMM, just identified by a binary instrument w for a
  // binary d: with A_w, B_w the count totals of d = 0 and d = 1 among
  // the rows with instrument w, the multiplicative moments give exp(-b) =
  // (n A1 - n1 (A0 + A1)) / (n1 (B0 + B1) - n B1); with z = x (no
  // endogeneity) both residuals give the log ratio of the two means
  {
    const int n = 60;
    rng r{20261116};
    std::vector<double> x(2 * n), z(2 * n), y(n);
    std::vector<int> c(n);
    double a0 = 0, a1 = 0, bb0 = 0, bb1 = 0, n1 = 0, s0 = 0, s1 = 0, m0 = 0, m1 = 0;
    for (int i = 0; i < n; ++i) {
      const double w = runi(r) < 0.5 ? 1 : 0;
      const double d = runi(r) < 0.3 + 0.4 * w ? 1 : 0;
      y[i] = rpois(r, std::exp(0.5 + 0.7 * d));
      x[i] = 1;
      x[n + i] = d;
      z[i] = 1;
      z[n + i] = w;
      c[i] = i;
      n1 += w;
      if (w > 0) {
        (d > 0 ? bb1 : a1) += y[i];
      } else {
        (d > 0 ? bb0 : a0) += y[i];
      }
      (d > 0 ? s1 : s0) += y[i];
      (d > 0 ? m1 : m0) += 1;
    }
    const double t = (n * a1 - n1 * (a0 + a1)) / (n1 * (bb0 + bb1) - n * bb1);
    const double b0[2] = {std::log((s0 + s1) / n), 0};
    double b[2], v1[4], v2[4], st[6];
    tk_ieq(T, "apx.gm", clab_expgmm(y.data(), x.data(), 2, z.data(), 2, n, nullptr, b0, c.data(), n,
                                    0, 2, 100, 1e-12, b, v2, st), CLAB_OK);
    tk_eq(T, "apx.gmb", b[1], -std::log(t), 1e-9);
    tk_nan(T, "apx.gmj", st[0]);
    tk_ieq(T, "apx.gmc", static_cast<int>(st[3]), 1);
    // just identified: the one-step sandwich is the efficient variance
    clab_expgmm(y.data(), x.data(), 2, z.data(), 2, n, nullptr, b0, c.data(), n, 0, 1, 100, 1e-12,
                b, v1, st);
    tk_veq(T, "apx.gmv", v1, v2, 4, 1e-7);
    tk_ieq(T, "apx.gmx", clab_expgmm(y.data(), x.data(), 2, x.data(), 2, n, nullptr, b0, c.data(),
                                     n, 0, 2, 100, 1e-12, b, v1, st), CLAB_OK);
    tk_eq(T, "apx.gmxm", b[1], std::log((s1 / m1) / (s0 / m0)), 1e-9);
    clab_expgmm(y.data(), x.data(), 2, x.data(), 2, n, nullptr, b0, c.data(), n, 1, 3, 100, 1e-12,
                b, v1, st);
    tk_eq(T, "apx.gmxa", b[1], std::log((s1 / m1) / (s0 / m0)), 1e-9);
    tk_ieq(T, "apx.gmk", clab_expgmm(y.data(), x.data(), 2, z.data(), 1, n, nullptr, b0, c.data(),
                                     n, 0, 2, 100, 1e-12, b, v1, st), CLAB_EARG);
    // the cluster weight is singular with fewer clusters than moments
    std::vector<int> c1(n, 0);
    c1[0] = 1;
    std::vector<double> z3(3 * n);
    for (int i = 0; i < 2 * n; ++i) z3[i] = z[i];
    for (int i = 0; i < n; ++i) z3[2 * n + i] = runi(r);
    tk_ieq(T, "apx.gms", clab_expgmm(y.data(), x.data(), 2, z3.data(), 3, n, nullptr, b0, c1.data(),
                                     2, 0, 2, 100, 1e-12, b, v1, st), CLAB_EDOM);
  }
  // the over-identification tests: Sargan is n times the share of the
  // squared 2SLS residuals that the instruments explain (the residual's
  // regression on them, clab_ivfe without endogenous columns); just
  // identified it is 0, J is NaN and the two-step estimate is 2SLS
  {
    const int n = 50;
    rng r{20261117};
    std::vector<double> y(n), xe(n), zi(2 * n);
    std::vector<int> g(n), c(n);
    for (int i = 0; i < n; ++i) {
      g[i] = i % 5;
      c[i] = i % 10;
      zi[i] = runi(r, -1, 1);
      zi[n + i] = runi(r, -1, 1);
      const double e = runi(r, -1, 1);
      xe[i] = zi[i] - 0.5 * zi[n + i] + e + 0.1 * g[i];
      y[i] = 0.5 * xe[i] + e + runi(r, -1, 1);
    }
    double b = 0, b2 = 0, st[3], s3[3], f[4];
    std::vector<double> res(n);
    tk_ieq(T, "apx.oi", clab_ivoid(y.data(), xe.data(), 1, nullptr, 0, zi.data(), 2, n, g.data(), 5,
                                   nullptr, 0, c.data(), 10, nullptr, 0, 1e-13, &b2, st), CLAB_OK);
    clab_ivfe(y.data(), xe.data(), 1, nullptr, 0, zi.data(), 2, n, nullptr, g.data(), 5, nullptr, 0,
              nullptr, 0, nullptr, 0, 1e-13, &b, nullptr, f, nullptr, res.data(), nullptr, nullptr,
              s3);
    double uu = 0;
    for (int i = 0; i < n; ++i) uu += res[i] * res[i];
    double ba[2], sa[3];
    clab_ivfe(res.data(), nullptr, 0, zi.data(), 2, nullptr, 0, n, nullptr, g.data(), 5, nullptr, 0,
              nullptr, 0, nullptr, 0, 1e-13, ba, nullptr, nullptr, nullptr, nullptr, nullptr,
              nullptr, sa);
    tk_eq(T, "apx.oisg", st[0], n * (1 - sa[2] / uu), 1e-9);
    tk_ge(T, "apx.oij", st[1], 0);
    tk_eq(T, "apx.oidf", st[2], 1, 0);
    tk_ieq(T, "apx.oi1", clab_ivoid(y.data(), xe.data(), 1, nullptr, 0, zi.data(), 1, n, g.data(), 5,
                                    nullptr, 0, c.data(), 10, nullptr, 0, 1e-13, &b2, st), CLAB_OK);
    clab_ivfe(y.data(), xe.data(), 1, nullptr, 0, zi.data(), 1, n, nullptr, g.data(), 5, nullptr, 0,
              nullptr, 0, nullptr, 0, 1e-13, &b, nullptr, f, nullptr, nullptr, nullptr, nullptr, s3);
    tk_eq(T, "apx.oi1b", b2, b, 1e-10);
    tk_eq(T, "apx.oi1s", st[0], 0, 1e-9);
    tk_nan(T, "apx.oi1j", st[1]);
  }
  // the rk statistic of two endogenous columns does not change when they
  // are replaced by invertible combinations of themselves (Kleibergen and
  // Paap 2006, p. 103: the normalised estimator is invariant to such
  // transformations), and is the same in both rows
  {
    const int n = 60;
    rng r{20261118};
    std::vector<double> y(n), xe(2 * n), xa(2 * n), zi(3 * n);
    std::vector<int> g(n), c(n);
    for (int i = 0; i < n; ++i) {
      g[i] = i % 3;
      c[i] = i % 12;
      for (int a = 0; a < 3; ++a) zi[a * n + i] = runi(r, -1, 1);
      xe[i] = zi[i] + 0.5 * zi[n + i] + runi(r, -1, 1);
      xe[n + i] = 0.3 * zi[i] + zi[2 * n + i] + runi(r, -1, 1);
      xa[i] = 2 * xe[i] + xe[n + i];
      xa[n + i] = 0.5 * xe[i] + 3 * xe[n + i];
      y[i] = xe[i] - xe[n + i] + runi(r, -1, 1);
    }
    double b[2], f1[8], f2[8], st[3];
    tk_ieq(T, "apx.rk", clab_ivfe(y.data(), xe.data(), 2, nullptr, 0, zi.data(), 3, n, nullptr,
                                  g.data(), 3, nullptr, 0, c.data(), 12, nullptr, 0, 1e-13, b,
                                  nullptr, f1, nullptr, nullptr, nullptr, nullptr, st), CLAB_OK);
    clab_ivfe(y.data(), xa.data(), 2, nullptr, 0, zi.data(), 3, n, nullptr, g.data(), 3, nullptr, 0,
              c.data(), 12, nullptr, 0, 1e-13, b, nullptr, f2, nullptr, nullptr, nullptr, nullptr, st);
    tk_gt(T, "apx.rkp", f1[2], 0);
    tk_eq(T, "apx.rki", f2[2], f1[2], 1e-8);
    tk_eq(T, "apx.rk2", f1[3], f1[2], 0);
  }
  // runs: days 1, 2, 3 of a run, a gap breaks it, a new unit starts afresh
  {
    const int u[6] = {0, 0, 0, 0, 1, 1};
    const double t[6] = {1, 2, 3, 5, 3, 4}, e[6] = {1, 1, 1, 1, 0, 1};
    int run[6], ord[6];
    tk_ieq(T, "apx.rl", clab_runlen(u, t, e, 6, 1, run, ord), CLAB_OK);
    tk_ieq(T, "apx.rl3", run[2], 3);
    tk_ieq(T, "apx.rlg", run[3], 1);
    tk_ieq(T, "apx.rl0", run[4], 0);
    tk_ieq(T, "apx.rlu", run[5], 1);
    const double td[6] = {1, 1, 3, 5, 3, 4};
    tk_ieq(T, "apx.rld", clab_runlen(u, td, e, 6, 1, run, ord), CLAB_EDOM);
  }
  // quadrants from 0, the first left out; 360 is 0; by area
  {
    const double d[4] = {45, 135, 359.5, 360};
    double o[12];
    tk_ieq(T, "apx.wb", clab_windbin(d, 4, 90, 0, nullptr, 1, 0, 0, o), CLAB_OK);
    tk_eq(T, "apx.wb1", o[0 * 4 + 1], 1, 0);
    tk_eq(T, "apx.wb3", o[2 * 4 + 2], 1, 0);
    tk_eq(T, "apx.wb0", o[0 * 4 + 3] + o[1 * 4 + 3] + o[2 * 4 + 3], 0, 0);
    const int a[2] = {0, 1};
    double o2[12];
    clab_windbin(d, 2, 90, 0, a, 2, 0, 0, o2);
    tk_eq(T, "apx.wba", o2[(3 + 0) * 2 + 1], 1, 0);
    // the paper's edges: 90 is in the first quadrant, 91 in the second
    const double de[3] = {90, 91, 0};
    double oe[9];
    clab_windbin(de, 3, 90, 0, nullptr, 1, 0, 1, oe);
    tk_eq(T, "apx.wbe0", oe[0] + oe[3] + oe[6], 0, 0);
    tk_eq(T, "apx.wbe1", oe[1], 1, 0);
    tk_eq(T, "apx.wbe2", oe[2] + oe[5] + oe[8], 0, 0);
    tk_ieq(T, "apx.wbw", clab_windbin(d, 4, 70, 0, nullptr, 1, 0, 0, o), CLAB_EARG);
  }
  // the same unit a day later; a repeated day is refused
  {
    const int u[4] = {0, 0, 1, 1};
    const double t[4] = {1, 2, 1, 3};
    int idx[4], ord[4];
    tk_ieq(T, "apx.sh", clab_dshift(u, t, 4, 1, idx, ord), CLAB_OK);
    tk_ieq(T, "apx.sh0", idx[0], 1);
    tk_ieq(T, "apx.sh1", idx[1], -1);
    tk_ieq(T, "apx.sh2", idx[2], -1);
    const double td[4] = {1, 1, 1, 3};
    tk_ieq(T, "apx.shd", clab_dshift(u, td, 4, 1, idx, ord), CLAB_EDOM);
  }
  // a road running east: north of it is its left; a south wind blows north
  {
    const double ax[3] = {5, 5, 5}, ay[3] = {1, 1, -1}, wd[3] = {180, 90, 180};
    const double x0 = 0, y0 = 0, x1 = 10, y1 = 0;
    int s[3], sd[3], dn[3];
    tk_ieq(T, "apx.ud", clab_updwn(ax, ay, wd, 3, &x0, &y0, &x1, &y1, 1, 2, 60, 0, s, sd, dn, nullptr), CLAB_OK);
    tk_ieq(T, "apx.udd", dn[0], 1);
    tk_ieq(T, "apx.uds", sd[0], 1);
    tk_ieq(T, "apx.uda", dn[1], -1);
    tk_ieq(T, "apx.udu", dn[2], 0);
    const double far[1] = {50};
    tk_ieq(T, "apx.udf", clab_updwn(ax, far, wd, 1, &x0, &y0, &x1, &y1, 1, 2, 60, 0, s, sd, dn, nullptr), CLAB_OK);
    tk_ieq(T, "apx.udn", dn[0], -2);
    // near two points: dropped, or given to the nearer
    const double px0[2] = {0, 10}, py0[2] = {0, 0}, pn[2] = {NAN, NAN};
    const double qx[1] = {4}, qy[1] = {0}, qw[1] = {270};
    double qd[1];
    clab_updwn(qx, qy, qw, 1, px0, py0, pn, pn, 2, 20, 60, 0, s, sd, dn, qd);
    tk_ieq(T, "apx.tw0", dn[0], -3);
    clab_updwn(qx, qy, qw, 1, px0, py0, pn, pn, 2, 20, 60, 1, s, sd, dn, qd);
    tk_ieq(T, "apx.tw1", s[0], 0);
    tk_eq(T, "apx.twd", qd[0], 4, 1e-15);
    tk_ieq(T, "apx.twdn", dn[0], 1);
    tk_ieq(T, "apx.udt", clab_updwn(ax, ay, wd, 3, &x0, &y0, &x1, &y1, 1, 2, 95, 0, s, sd, dn, nullptr), CLAB_EARG);
  }
}
