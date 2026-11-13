// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/ppx.cpp: the isotropic weight at an edge and a corner, the
// cross-K of independent types against pi r^2, the segregation statistic
// by its definition, quadrat counts, Morisita's index and the segregation
// indices by hand.
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "tchk.h"

namespace {

using vd = std::vector<double>;
using vi = std::vector<int>;
const double PI = 3.14159265358979323846;

void famtest(tk &T, int big) {
  rng g = {5150};
  // the weight: a point on an edge sees half its circle, one at a corner a
  // quarter; one pair each, so K = A w / 1
  {
    const double win[4] = {0, 10, 0, 10}, r = 1;
    const vd x = {5, 5.5}, y = {0, 0}, xc = {0, 0.3}, yc = {0, 0.4};
    const vi m = {0, 1};
    double k = 0, seg[2];
    clab_crossk(x.data(), y.data(), m.data(), 2, 0, 1, win, &r, 1, 0, 0, 1, &k, seg);
    tk_eq(T, "ck.edge", k, 100 * 2, 1e-12);
    clab_crossk(xc.data(), yc.data(), m.data(), 2, 0, 1, win, &r, 1, 0, 0, 1, &k, seg);
    tk_eq(T, "ck.corner", k, 100 * 4, 1e-12);
    // a point inside, its circle crossing one side at distance 0.5: the
    // arc outside is 2 acos(0.5) = 2 pi / 3
    const vd xi = {5, 5}, yi = {0.5, 1.5};
    clab_crossk(xi.data(), yi.data(), m.data(), 2, 0, 1, win, &r, 1, 0, 0, 1, &k, seg);
    tk_eq(T, "ck.arc", k, 100 * 2 * PI / (2 * PI - 2 * PI / 3), 1e-12);
  }
  // independent types in a square: K near pi r^2, far from the edges or
  // not
  {
    const int n = big ? 4000 : 800;
    vd x(n), y(n);
    vi m(n);
    for (int a = 0; a < n; ++a) {
      x[a] = runi(g, 0, 1);
      y[a] = runi(g, 0, 1);
      m[a] = a % 2;
    }
    const double win[4] = {0, 1, 0, 1}, r[2] = {0.05, 0.1};
    double k[2], seg[2];
    tk_ieq(T, "ck.ok",
           clab_crossk(x.data(), y.data(), m.data(), n, 0, 1, win, r, 2, 0, 0, 1, k, seg), CLAB_OK);
    tk_lt(T, "ck.csr1", std::fabs(k[0] / (PI * 0.0025) - 1), big ? 0.05 : 0.12);
    tk_lt(T, "ck.csr2", std::fabs(k[1] / (PI * 0.01) - 1), big ? 0.03 : 0.08);
    double kii[2];
    clab_crossk(x.data(), y.data(), m.data(), n, 1, 1, win, r, 2, 0, 0, 1, kii, seg);
    tk_lt(T, "ck.same", std::fabs(kii[1] / (PI * 0.01) - 1), big ? 0.03 : 0.08);
  }
  // the segregation statistic by its definition on four points
  {
    const vd x = {0, 1, 0, 1}, y = {0, 0, 1, 1};
    const vi m = {0, 0, 1, 1};
    const double win[4] = {0, 1, 0, 1}, r = 0.5, sg = 1;
    double k, seg[2], seg2[2];
    clab_crossk(x.data(), y.data(), m.data(), 4, 0, 1, win, &r, 1, sg, 0, 1, &k, seg);
    // point 0: others at 1 (type 0), 1 (type 1) and sqrt 2 (type 1)
    const double a = std::exp(-0.5), b = std::exp(-1.0);
    const double p0 = a / (2 * a + b), p1 = (a + b) / (2 * a + b);
    tk_eq(T, "seg.def", seg[0], 4 * ((p0 - 0.5) * (p0 - 0.5) + (p1 - 0.5) * (p1 - 0.5)), 1e-15);
    tk_nan(T, "seg.nop", seg[1]);
    clab_crossk(x.data(), y.data(), m.data(), 4, 0, 1, win, &r, 1, sg, 19, 2, &k, seg);
    clab_crossk(x.data(), y.data(), m.data(), 4, 0, 1, win, &r, 1, sg, 19, 2, &k, seg2);
    tk_true(T, "seg.p", seg[1] == seg2[1] && seg[1] >= 0.05 && seg[1] <= 1);
    tk_ieq(T, "ck.out",
           clab_crossk(x.data(), y.data(), m.data(), 4, 0, 2, win, &r, 1, 0, 0, 1, &k, seg),
           CLAB_EDOM);
  }
  // quadrats by hand
  {
    const vd x = {0.1, 0.2, 0.9, 1.0, 0.6}, y = {0.1, 0.9, 0.9, 1.0, 0.2};
    const double win[4] = {0, 1, 0, 1};
    vd c(4), st(3);
    clab_qdrtst(x.data(), y.data(), 5, win, 2, 2, c.data(), st.data());
    // column 1: top 1, bottom 1; column 2: top 2, bottom 1
    const vd want = {1, 1, 2, 1};
    tk_veq(T, "qd.cnt", c.data(), want.data(), 4, 0);
    tk_eq(T, "qd.x2", st[0], (3 * 0.0625 + 0.5625) / 1.25, 1e-15);
    tk_eq(T, "qd.lam", st[2], 5, 0);
    tk_ieq(T, "qd.out", clab_qdrtst(x.data(), y.data(), 5, win, 1, 1, c.data(), st.data()),
           CLAB_EARG);
  }
  // Morisita: all in one quadrat gives n, an even spread 1 - (n - 1) / (s
  // - 1)
  {
    const vd c1 = {0, 0, 0, 10}, c2 = {5, 5, 5, 5};
    double o[5];
    clab_morsta(c1.data(), 4, 0.2, 9, o);
    tk_eq(T, "mo.clump", o[0], 4, 1e-15);
    tk_eq(T, "mo.chi", o[4], 4 * 9 + 4 - 10, 1e-14);
    clab_morsta(c2.data(), 4, 0.2, 9, o);
    tk_eq(T, "mo.even", o[0], 4 * (100 - 20) / (400.0 - 20), 1e-15);
    tk_ieq(T, "mo.one", clab_morsta(c1.data(), 1, 0.2, 9, o), CLAB_EARG);
  }
  // segregation indices: complete segregation 1, evenness 0
  {
    const vd seg = {10, 0, 0, 10}, ev = {5, 10, 5, 10};
    vd d(4), gr(8), dm(4);
    clab_dissix(seg.data(), 2, 2, 0.5, nullptr, nullptr, nullptr, d.data(), gr.data(), dm.data());
    tk_eq(T, "dx.d1", d[1], 1, 1e-15);
    const vd one = {1, 1, 1, 1};
    tk_veq(T, "dx.all1", gr.data(), one.data(), 4, 1e-15);
    clab_dissix(ev.data(), 2, 2, 0.5, nullptr, nullptr, nullptr, d.data(), gr.data(), dm.data());
    tk_eq(T, "dx.d0", d[1], 0, 0);
    tk_lt(T, "dx.h0", std::fabs(gr[1]) + std::fabs(gr[2]) + std::fabs(gr[3]), 1e-14);
    // Morrill on two neighbouring units of opposite shares: 1 - 1 = 0
    const int rp[3] = {0, 1, 2}, ci[2] = {1, 0};
    const double w[2] = {1, 1};
    clab_dissix(seg.data(), 2, 2, 0.5, rp, ci, w, d.data(), gr.data(), dm.data());
    tk_eq(T, "dx.morrill", dm[1], 0, 1e-15);
    tk_ieq(
        T, "dx.delta",
        clab_dissix(seg.data(), 2, 2, 1, nullptr, nullptr, nullptr, d.data(), gr.data(), dm.data()),
        CLAB_EDOM);
  }
}

}  // namespace

void t_ppx(tk &T, int big) { famtest(T, big); }
