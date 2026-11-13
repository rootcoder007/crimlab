// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/hot.cpp: the surface against the direct double sum and its
// integral, the edge corrections by hand, DBSCAN on a hand-made pattern,
// hot-cell persistence by hand, and the shrinkage factor's optimality.
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "tchk.h"

namespace {

using vd = std::vector<double>;
using vi = std::vector<int>;
const double PI = 3.14159265358979323846;

double ncdf(double z) { return 0.5 * std::erfc(-z / std::sqrt(2.0)); }

void famtest(tk &T, int big) {
  rng g = {909};
  // the grid surface: the direct sum, and its integral the total weight
  {
    const int n = big ? 300 : 60, na = 41, nb = 37;
    vd x(n), y(n), w(n), ax(na), ay(nb), z(na * nb);
    for (int i = 0; i < n; ++i) {
      x[i] = runi(g, 0, 10);
      y[i] = runi(g, 0, 6);
      w[i] = runi(g, 0.5, 2);
    }
    for (int a = 0; a < na; ++a) ax[a] = -6 + a * 22.0 / (na - 1);
    for (int b = 0; b < nb; ++b) ay[b] = -6 + b * 18.0 / (nb - 1);
    const double hx = 0.9, hy = 0.6;
    tk_ieq(T, "kd.ok",
           clab_kdsurf(x.data(), y.data(), w.data(), n, hx, hy, ax.data(), na, ay.data(), nb, 1,
                       nullptr, 0, 0, z.data()),
           CLAB_OK);
    double err = 0, tot = 0, sw = 0;
    for (int b = 0; b < nb; ++b)
      for (int a = 0; a < na; ++a) {
        double s = 0;
        for (int j = 0; j < n; ++j)
          s += w[j] * std::exp(-0.5 * ((ax[a] - x[j]) * (ax[a] - x[j]) / (hx * hx) +
                                       (ay[b] - y[j]) * (ay[b] - y[j]) / (hy * hy)));
        s /= 2 * PI * hx * hy;
        err = std::fmax(err, std::fabs(z[a + b * na] - s) / std::fmax(s, 1e-300));
        tot += z[a + b * na];
      }
    tk_lt(T, "kd.direct", err, 1e-13);
    for (int j = 0; j < n; ++j) sw += w[j];
    // the trapezoid rule is spectrally accurate for a Gaussian far inside
    tk_eq(T, "kd.mass", tot * (22.0 / (na - 1)) * (18.0 / (nb - 1)), sw, 1e-8);
  }
  // at points: leave one out, and the two edge corrections by hand
  {
    const vd x = {0.2, 0.4, 0.7}, y = {0.3, 0.5, 0.2};
    const double win[4] = {0, 1, 0, 1}, h = 0.2;
    vd o(3), u(3), d(3);
    clab_kdsurf(x.data(), y.data(), nullptr, 3, h, h, x.data(), 3, y.data(), 3, 0, win, 0, 1,
                o.data());
    clab_kdsurf(x.data(), y.data(), nullptr, 3, h, h, x.data(), 3, y.data(), 3, 0, win, 1, 1,
                u.data());
    clab_kdsurf(x.data(), y.data(), nullptr, 3, h, h, x.data(), 3, y.data(), 3, 0, win, 2, 1,
                d.data());
    auto k = [&](int a, int b) {
      return std::exp(-((x[a] - x[b]) * (x[a] - x[b]) + (y[a] - y[b]) * (y[a] - y[b])) /
                      (2 * h * h)) /
             (2 * PI * h * h);
    };
    auto e = [&](int a) {
      return (ncdf((1 - x[a]) / h) - ncdf(-x[a] / h)) * (ncdf((1 - y[a]) / h) - ncdf(-y[a] / h));
    };
    const vd wo = {k(0, 1) + k(0, 2), k(1, 0) + k(1, 2), k(2, 0) + k(2, 1)};
    const vd wu = {wo[0] / e(0), wo[1] / e(1), wo[2] / e(2)};
    const vd wd = {k(0, 1) / e(1) + k(0, 2) / e(2), k(1, 0) / e(0) + k(1, 2) / e(2),
                   k(2, 0) / e(0) + k(2, 1) / e(1)};
    tk_veq(T, "kd.loo", o.data(), wo.data(), 3, 1e-14);
    tk_veq(T, "kd.unif", u.data(), wu.data(), 3, 1e-14);
    tk_veq(T, "kd.diggle", d.data(), wd.data(), 3, 1e-14);
    tk_ieq(T, "kd.h",
           clab_kdsurf(x.data(), y.data(), nullptr, 3, 0, h, x.data(), 3, y.data(), 3, 0, win, 0, 0,
                       o.data()),
           CLAB_EDOM);
    tk_ieq(T, "kd.n",
           clab_kdsurf(x.data(), y.data(), nullptr, 0, h, h, x.data(), 3, y.data(), 3, 0, win, 0, 0,
                       o.data()),
           CLAB_EARG);
    tk_ieq(T, "kd.win",
           clab_kdsurf(x.data(), y.data(), nullptr, 3, h, h, x.data(), 3, y.data(), 3, 0, nullptr,
                       1, 0, o.data()),
           CLAB_EARG);
  }
  // DBSCAN: two groups, a border point and noise
  {
    const vd x = {0, 0.5, 1, 0.5, 1.9, 10, 10.4, 10.8, 30}, y = {0, 0, 0, 0.5, 0, 0, 0, 0.3, 0};
    vi cl(9), core(9);
    int k = 0;
    tk_ieq(T, "db.ok", clab_dbscn(x.data(), y.data(), 9, 1, 3, 1, cl.data(), core.data(), &k),
           CLAB_OK);
    const vi want = {0, 0, 0, 0, 0, 1, 1, 1, -1}, wc = {1, 1, 1, 1, 0, 1, 1, 1, 0};
    bool ok = k == 2;
    for (int i = 0; i < 9; ++i) ok = ok && cl[i] == want[i] && core[i] == wc[i];
    tk_true(T, "db.hand", ok);
    clab_dbscn(x.data(), y.data(), 9, 1, 3, 0, cl.data(), core.data(), &k);
    tk_true(T, "db.noborder", cl[4] == -1 && cl[2] == 0);
    tk_ieq(T, "db.eps", clab_dbscn(x.data(), y.data(), 9, 0, 3, 1, cl.data(), core.data(), &k),
           CLAB_EDOM);
  }
  // persistence by hand: three cells over three periods
  {
    // cell (0,0): 3, 2, 0; cell (1,0): 1, 2, 1; cell (5,5): 0, 0, 4
    vd x, y;
    vi p;
    auto add = [&](double a, double b, int q, int m) {
      for (int r = 0; r < m; ++r) {
        x.push_back(a);
        y.push_back(b);
        p.push_back(q);
      }
    };
    add(0.5, 0.5, 0, 3);
    add(0.5, 0.5, 1, 2);
    add(1.5, 0.5, 0, 1);
    add(1.5, 0.5, 1, 2);
    add(1.5, 0.5, 2, 1);
    add(5.5, 5.5, 2, 4);
    const int n = static_cast<int>(x.size());
    vi cid(n), hot(3 * n);
    vd ix(n), iy(n), per(27), cs(3 * n);
    int nc = 0;
    tk_ieq(T, "hp.ok",
           clab_hsper(x.data(), y.data(), p.data(), n, 3, 1, 0, 0, 0.34, cid.data(), ix.data(),
                      iy.data(), &nc, hot.data(), per.data(), cs.data()),
           CLAB_OK);
    tk_ieq(T, "hp.nc", nc, 3);
    // k = ceil(0.34 * 3) = 2: period 1 hot cells 0 and 1 (3, 1); period 2
    // a tie at 2 (cells 0 and 1, both hot); period 3 cells 2 and 1
    const vi wh = {1, 1, 0, 1, 1, 0, 0, 1, 1};
    bool ok = true;
    for (int q = 0; q < 9; ++q) ok = ok && hot[q] == wh[q];
    tk_true(T, "hp.hot", ok);
    const vd p1 = {4, 2, 1, 1.0 / 3, 2, 0, 0, 1, std::nan("")};
    tk_veq(T, "hp.per1", per.data(), p1.data(), 9, 1e-15);
    tk_eq(T, "hp.jac3", per[18 + 8], 1.0 / 3, 1e-15);
    const vd c1 = {4, 3, 3};
    tk_veq(T, "hp.cell", cs.data() + 3, c1.data(), 3, 0);
    tk_ieq(T, "hp.top",
           clab_hsper(x.data(), y.data(), p.data(), n, 3, 1, 0, 0, 0, cid.data(), ix.data(),
                      iy.data(), &nc, hot.data(), per.data(), cs.data()),
           CLAB_EDOM);
  }
  // shrinkage: with a noise law that holds exactly, the loss is (1 - B)^2
  // S_e + B^2 S_theta and B* minimises it
  {
    const vd th = {10, 20, 30, 40}, e = {3, -3, -3, 3};
    vd y(4), fit(6), s(4), f(4);
    for (int i = 0; i < 4; ++i) y[i] = th[i] + e[i];
    const double se = 36.0 / 4, ybar = 25;
    clab_ebshr(y.data(), nullptr, 4, se, std::nan(""), fit.data(), s.data(), f.data());
    tk_eq(T, "eb.mean", fit[0], ybar, 1e-15);
    // total scatter = (500 + 36) / 4; signal = 125
    tk_eq(T, "eb.sig", fit[3], 125, 1e-13);
    tk_eq(T, "eb.bstar", fit[4], 9.0 / 134, 1e-15);
    auto loss = [&](double b) {
      double l = 0;
      for (int i = 0; i < 4; ++i) {
        const double v = (1 - b) * y[i] + b * ybar - th[i];
        l += v * v;
      }
      return l;
    };
    bool mn = true;
    for (int q = 0; q <= 20; ++q) mn = mn && loss(fit[4]) <= loss(q / 20.0) + 1e-12;
    tk_true(T, "eb.min", mn);
    tk_eq(T, "eb.closed", loss(0.3), 0.49 * 36 + 0.09 * 500, 1e-12);
    tk_eq(T, "eb.fall", f[3], fit[5] * (y[3] - ybar), 1e-15);
    tk_ieq(T, "eb.b", clab_ebshr(y.data(), nullptr, 4, 1, 2, fit.data(), s.data(), f.data()),
           CLAB_EDOM);
  }
}

}  // namespace

void t_hot(tk &T, int big) { famtest(T, big); }
