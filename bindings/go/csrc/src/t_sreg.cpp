// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/sreg.cpp: Moran's I against its definition and its mean over
// every permutation, the bivariate statistic at y = x, the lag model's
// likelihood maximised against a grid, two-stage least squares on a
// noiseless model, and the effects against their closed form.
#include <algorithm>
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "tchk.h"

namespace {

using vd = std::vector<double>;
using vi = std::vector<int>;

// a ring of n units, each with its two neighbours, row-standardised
void ring(int n, vi &rp, vi &ci, vd &w) {
  rp.assign(n + 1, 0);
  ci.clear();
  w.clear();
  for (int i = 0; i < n; ++i) {
    int a = (i + n - 1) % n, b = (i + 1) % n;
    if (a > b) std::swap(a, b);
    ci.push_back(a);
    ci.push_back(b);
    w.push_back(0.5);
    w.push_back(0.5);
    rp[i + 1] = 2 * (i + 1);
  }
}

void famtest(tk &T, int big) {
  (void)big;
  rng g = {3131};
  vi rp, ci;
  vd w;
  // Moran's I by its definition; E[I] = -1 / (n - 1) over every order
  {
    const int n = 6;
    ring(n, rp, ci, w);
    vd x = {3, 1, 4, 1, 5, 9}, out(10);
    tk_ieq(T, "gm.ok",
           clab_gmoran(x.data(), n, 1, rp.data(), ci.data(), w.data(), 0, 1, out.data()), CLAB_OK);
    double m = 0;
    for (double v : x) m += v;
    m /= n;
    double num = 0, den = 0;
    for (int i = 0; i < n; ++i) {
      den += (x[i] - m) * (x[i] - m);
      num += 0.5 * (x[i] - m) * ((x[(i + 1) % n] - m) + (x[(i + n - 1) % n] - m));
    }
    tk_eq(T, "gm.def", out[0], num / den, 1e-14);
    tk_eq(T, "gm.s0", out[9], 6, 0);
    vi pm = {0, 1, 2, 3, 4, 5};
    double s = 0;
    int k = 0;
    vd xp(n), o(10);
    do {
      for (int i = 0; i < n; ++i) xp[i] = x[pm[i]];
      clab_gmoran(xp.data(), n, 1, rp.data(), ci.data(), w.data(), 0, 1, o.data());
      s += o[0];
      ++k;
    } while (std::next_permutation(pm.begin(), pm.end()));
    tk_eq(T, "gm.mean", s / k, -1.0 / (n - 1), 1e-13);
    tk_eq(T, "gm.e", out[1], -0.2, 1e-15);
    // the permutation p-value: in range and repeatable
    vd o1(10), o2(10);
    clab_gmoran(x.data(), n, 1, rp.data(), ci.data(), w.data(), 99, 3, o1.data());
    clab_gmoran(x.data(), n, 1, rp.data(), ci.data(), w.data(), 99, 3, o2.data());
    tk_true(T, "gm.perm", o1[8] == o2[8] && o1[8] >= 0.01 && o1[8] <= 1);
    // the bivariate statistic with y = x is Moran's I
    vd b(6);
    tk_ieq(T, "bv.ok",
           clab_bvmor(x.data(), x.data(), n, rp.data(), ci.data(), w.data(), 0, 1, b.data()),
           CLAB_OK);
    tk_eq(T, "bv.self", b[0], out[0], 1e-14);
    const vd c = {1, 1, 1, 1, 1, 1};
    tk_ieq(T, "bv.const",
           clab_bvmor(x.data(), c.data(), n, rp.data(), ci.data(), w.data(), 0, 1, b.data()),
           CLAB_EDOM);
    const int bad[5] = {0, 0, 1, 1, 1};
    tk_ieq(T, "gm.self", clab_gmoran(x.data(), 4, 1, bad, ci.data(), w.data(), 0, 1, out.data()),
           CLAB_EDOM);
  }
  // the lag model: the concentrated likelihood is largest at the estimate
  {
    const int n = 30, p = 2;
    ring(n, rp, ci, w);
    vd x(n * p), y(n), e(n);
    for (int i = 0; i < n; ++i) {
      x[i] = 1;
      x[n + i] = runi(g, -2, 2);
      e[i] = runi(g, -1, 1);
    }
    // y = (I - 0.5 W)^-1 (X b + e) by fixed-point iteration
    vd v(n), z(n, 0.0);
    for (int i = 0; i < n; ++i) v[i] = 1 + 2 * x[n + i] + e[i];
    for (int it = 0; it < 200; ++it) {
      for (int i = 0; i < n; ++i) z[i] = v[i] + 0.25 * (y[(i + 1) % n] + y[(i + n - 1) % n]);
      y = z;
    }
    vd cf(4 * p + 4), vc((p + 2) * (p + 2)), st(5);
    tk_ieq(T, "sar.ok",
           clab_sarml(y.data(), x.data(), n, p, rp.data(), ci.data(), w.data(), 0, 0, nullptr, 0,
                      cf.data(), vc.data(), st.data()),
           CLAB_OK);
    const double rho = cf[p];
    // the profile log-likelihood written out: OLS of y - r W y on X, and
    // log |I - r W| = sum log(1 - r cos(2 pi k / n)) on the ring
    auto prof = [&](double r) {
      vd yy(n);
      for (int i = 0; i < n; ++i) yy[i] = y[i] - r * 0.5 * (y[(i + 1) % n] + y[(i + n - 1) % n]);
      double sx = 0, sxx = 0, sy = 0, sxy = 0;
      for (int i = 0; i < n; ++i) {
        sx += x[n + i];
        sxx += x[n + i] * x[n + i];
        sy += yy[i];
        sxy += x[n + i] * yy[i];
      }
      const double b1 = (n * sxy - sx * sy) / (n * sxx - sx * sx), b0 = (sy - b1 * sx) / n;
      double sse = 0;
      for (int i = 0; i < n; ++i)
        sse += (yy[i] - b0 - b1 * x[n + i]) * (yy[i] - b0 - b1 * x[n + i]);
      double ld = 0;
      for (int k = 0; k < n; ++k)
        ld += std::log(1 - r * std::cos(2 * 3.14159265358979323846 * k / n));
      return -0.5 * n * std::log(sse / n) + ld;
    };
    bool mx = true;
    for (int q = -9; q <= 9; ++q) mx = mx && prof(rho) >= prof(rho + q * 1e-3);
    tk_true(T, "sar.max", mx);
    tk_lt(T, "sar.flat", std::fabs(prof(rho + 1e-5) - prof(rho - 1e-5)), 1e-8);
    // the full log-likelihood at the estimate equals the profile one
    tk_eq(T, "sar.ll", st[0], prof(rho) - 0.5 * n * (std::log(2 * 3.14159265358979323846) + 1),
          1e-10);
    // two-stage least squares with exact instruments recovers a noiseless
    // model
    vd y0(n, 0.0), h(n);
    for (int i = 0; i < n; ++i) v[i] = 1 + 2 * x[n + i];
    for (int it = 0; it < 300; ++it) {
      for (int i = 0; i < n; ++i) z[i] = v[i] + 0.25 * (y0[(i + 1) % n] + y0[(i + n - 1) % n]);
      y0 = z;
    }
    for (int i = 0; i < n; ++i) h[i] = 0.5 * (x[n + (i + 1) % n] + x[n + (i + n - 1) % n]);
    clab_sarml(y0.data(), x.data(), n, p, rp.data(), ci.data(), w.data(), 0, 1, h.data(), 1,
               cf.data(), vc.data(), st.data());
    const vd want = {1, 2, 0.5};
    tk_veq(T, "s2s.exact", cf.data(), want.data(), 3, 1e-9);
    tk_ieq(T, "sar.p",
           clab_sarml(y0.data(), x.data(), n, 0, rp.data(), ci.data(), w.data(), 0, 0, nullptr, 0,
                      cf.data(), vc.data(), st.data()),
           CLAB_EARG);
    tk_ieq(T, "s2s.err",
           clab_sarml(y0.data(), x.data(), n, p, rp.data(), ci.data(), w.data(), 1, 1, h.data(), 1,
                      cf.data(), vc.data(), st.data()),
           CLAB_EARG);
  }
  // effects on a ring: the total is (b + th) / (1 - rho), the direct the
  // mean diagonal of the inverse, by its series
  {
    const int n = 8;
    ring(n, rp, ci, w);
    const double rho = 0.4, b = 1.5, th = 0.5;
    vd out(6);
    tk_ieq(T, "sp.ok",
           clab_spill(n, rp.data(), ci.data(), w.data(), rho, &b, &th, 1, nullptr, out.data()),
           CLAB_OK);
    tk_eq(T, "sp.total", out[2], (b + th) / (1 - rho), 1e-14);
    // diagonal of (I - rho W)^-1 on the ring: mean over the eigenvalues
    double d1 = 0, d2 = 0;
    for (int k = 0; k < n; ++k) {
      const double l = std::cos(2 * 3.14159265358979323846 * k / n);
      d1 += 1 / (1 - rho * l);
      d2 += l / (1 - rho * l);
    }
    tk_eq(T, "sp.direct", out[0], (b * d1 + th * d2) / n, 1e-14);
    tk_nan(T, "sp.nose", out[3]);
    tk_ieq(T, "sp.sing",
           clab_spill(n, rp.data(), ci.data(), w.data(), 1, &b, &th, 1, nullptr, out.data()),
           CLAB_EDOM);
  }
}

}  // namespace

void t_sreg(tk &T, int big) { famtest(T, big); }
