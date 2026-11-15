// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/apx.cpp and of clab_ppmlfe and clab_fesep (cnt.cpp): closed
// forms on tiny panels, the absorbed Poisson fit against the same model
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
    double b = 0, vc = 0, f[3] = {0, 0, 0}, st[3];
    tk_ieq(T, "apx.iv", clab_ivfe(y, x, 1, nullptr, 0, z, 1, 5, nullptr, g, 1, nullptr, 0, c, 2,
                                  nullptr, 0, 1e-13, &b, &vc, f, nullptr, nullptr, nullptr, nullptr, st), CLAB_OK);
    tk_eq(T, "apx.ivb", b, szy / szx, 1e-12);
    // one instrument: the effective F is the Wald F
    tk_eq(T, "apx.ivef", f[2], f[0], 1e-10);
    tk_eq(T, "apx.ivkp", f[1], f[0], 0);
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
