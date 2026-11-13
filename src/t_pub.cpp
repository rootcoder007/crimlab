// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for the published-table kernels (src/pub.cpp) against
// values worked by hand, the adjusted p-values among them as p.adjust()
// defines them. Part of the self-test (selftest.cpp).
#include <cmath>
#include <limits>

#include "crimlab.h"
#include "tchk.h"

namespace {

const double nan = std::numeric_limits<double>::quiet_NaN();

void bands(tk &t) {
  const double lo[5] = {1, 6, nan, 61, nan}, hi[5] = {5, 15, nan, nan, 9};
  const int olo[5] = {0, 0, 0, 0, 1}, oup[5] = {0, 0, 0, 1, 0};
  double v[5], bl[5], bh[5];
  tk_ieq(t, "bandx.rc", clab_bandx(lo, hi, olo, oup, 5, 0, nan, 2, nan, v, bl, bh), CLAB_OK);
  tk_eq(t, "bandx.mid", v[1], 10.5, 0);
  tk_nan(t, "bandx.unread", v[2]);
  tk_eq(t, "bandx.openup", v[3], 91.5, 0);
  tk_eq(t, "bandx.openlo", v[4], 4.5, 0);
  tk_eq(t, "bandx.end", bh[3], 122, 0);
  clab_bandx(lo, hi, olo, oup, 5, 3, 100, 2, 2, v, bl, bh);
  tk_eq(t, "bandx.geo", v[0], std::sqrt(5.0), 1e-15);
  tk_eq(t, "bandx.cap", bh[3], 100, 0);
  tk_eq(t, "bandx.floor", bl[4], 2, 0);
  const double z0[1] = {0}, z4[1] = {4};
  const int no[1] = {0};
  clab_bandx(z0, z4, no, no, 1, 3, nan, 2, nan, v, bl, bh);
  tk_eq(t, "bandx.geo0", v[0], 2, 0);
  clab_bandx(lo, hi, olo, oup, 5, 1, nan, 2, nan, v, bl, bh);
  tk_eq(t, "bandx.lower", v[1], 6, 0);
  clab_bandx(lo, hi, olo, oup, 5, 2, nan, 2, nan, v, bl, bh);
  tk_eq(t, "bandx.upper", v[1], 15, 0);
  tk_ieq(t, "bandx.cap50", clab_bandx(lo, hi, olo, oup, 5, 0, 50, 2, nan, v, bl, bh), CLAB_EDOM);
  tk_ieq(t, "bandx.flr", clab_bandx(lo, hi, olo, oup, 5, 0, nan, 2, 20, v, bl, bh), CLAB_EDOM);
  tk_ieq(t, "bandx.rule", clab_bandx(lo, hi, olo, oup, 5, 4, nan, 2, nan, v, bl, bh), CLAB_EARG);
  tk_ieq(t, "bandx.fac", clab_bandx(lo, hi, olo, oup, 5, 0, nan, 0, nan, v, bl, bh), CLAB_EDOM);
}

void cells(tk &t) {
  const double val[4] = {100, 7, nan, nan}, lim[4] = {nan, nan, 5, nan};
  const int st[4] = {0, 0, 1, 2};
  double lo[4], hi[4];
  int so[4];
  tk_ieq(t, "pbbnd.rc", clab_pbbnd(val, st, lim, 4, 5, 0, nan, lo, hi, so), CLAB_OK);
  tk_eq(t, "pbbnd.lo", lo[0], 97.5, 0);
  tk_eq(t, "pbbnd.hi", hi[1], 9.5, 0);
  tk_eq(t, "pbbnd.sup", hi[2], 4, 0);
  tk_nan(t, "pbbnd.miss", lo[3]);
  tk_ieq(t, "pbbnd.st", so[0] * 1000 + so[1] * 100 + so[2] * 10 + so[3], 1123);
  clab_pbbnd(val, st, lim, 4, 5, 1, nan, lo, hi, so);
  tk_eq(t, "pbbnd.rand", lo[0], 96, 0);
  tk_eq(t, "pbbnd.floor0", lo[1], 3, 0);
  const double l2[4] = {nan, nan, nan, nan}, v3[1] = {3};
  clab_pbbnd(val, st, l2, 4, nan, 0, 10, lo, hi, so);
  tk_eq(t, "pbbnd.exact", hi[0], 100, 0);
  tk_eq(t, "pbbnd.suplim", hi[2], 9, 0);
  tk_ieq(t, "pbbnd.st0", so[0], 0);
  // rounding to 10 floors a small count's lower end at 0
  clab_pbbnd(v3, st, lim, 1, 10, 0, nan, lo, hi, so);
  tk_eq(t, "pbbnd.zero", lo[0], 0, 0);
  const double vn[1] = {-1};
  tk_ieq(t, "pbbnd.nolim", clab_pbbnd(val, st, l2, 4, nan, 0, nan, lo, hi, so), CLAB_EDOM);
  tk_ieq(t, "pbbnd.neg", clab_pbbnd(vn, st, lim, 1, nan, 0, nan, lo, hi, so), CLAB_EDOM);
  tk_ieq(t, "pbbnd.rnd", clab_pbbnd(val, st, lim, 4, -1, 0, nan, lo, hi, so), CLAB_EDOM);
  tk_ieq(t, "pbbnd.kind", clab_pbbnd(val, st, lim, 4, 5, 2, nan, lo, hi, so), CLAB_EARG);
  // the change between [10, 12] now and [8, 9] before
  const double nl[2] = {10, 10}, nh[2] = {12, 12}, pl[2] = {8, 0}, ph[2] = {9, 0};
  const double pop[2] = {1000, 1000};
  double o[16];
  tk_ieq(t, "pbchg.rc", clab_pbchg(nl, nh, pl, ph, 2, pop, pop, 1000, o), CLAB_OK);
  tk_eq(t, "pbchg.lo", o[0], 1, 0);
  tk_eq(t, "pbchg.hi", o[2], 4, 0);
  tk_eq(t, "pbchg.plo", o[4], 100 * (10.0 / 9 - 1), 1e-14);
  tk_eq(t, "pbchg.phi", o[6], 50, 1e-15);
  tk_nan(t, "pbchg.p0", o[5]);
  tk_nan(t, "pbchg.p0b", o[7]);
  tk_eq(t, "pbchg.rate", o[10], 12, 1e-15);
  tk_eq(t, "pbchg.prate", o[12], 8, 1e-15);
  clab_pbchg(nl, nh, pl, ph, 2, nullptr, nullptr, 1000, o);
  tk_nan(t, "pbchg.nopop", o[8]);
  const double p0[2] = {0, 1};
  tk_ieq(t, "pbchg.per", clab_pbchg(nl, nh, pl, ph, 2, pop, pop, 0, o), CLAB_EDOM);
  tk_ieq(t, "pbchg.pop", clab_pbchg(nl, nh, pl, ph, 2, p0, pop, 1, o), CLAB_EDOM);
  tk_ieq(t, "pbchg.n", clab_pbchg(nl, nh, pl, ph, -1, pop, pop, 1, o), CLAB_EARG);
}

void yoy(tk &t) {
  const double v[6] = {10, 12, 15, 5, 0, 7}, prd[6] = {2019, 2020, 2022, 2019, 2020, 2021};
  const int g[6] = {0, 0, 0, 1, 1, 1};
  double prev[6], chg[6], pct[6];
  tk_ieq(t, "yoypx.rc", clab_yoypx(v, prd, g, 6, 1, prev, chg, pct), CLAB_OK);
  tk_nan(t, "yoypx.first", prev[0]);
  tk_eq(t, "yoypx.prev", prev[1], 10, 0);
  tk_eq(t, "yoypx.pct", pct[1], 20, 1e-15);
  tk_nan(t, "yoypx.gap", prev[2]);
  tk_nan(t, "yoypx.group", prev[3]);
  tk_eq(t, "yoypx.chg", chg[4], -5, 0);
  tk_nan(t, "yoypx.zero", pct[5]);
  clab_yoypx(v, prd, g, 6, 2, prev, chg, pct);
  tk_eq(t, "yoypx.lag2", prev[2], 12, 0);
  const double pu[2] = {2020, 2019}, pi[1] = {std::numeric_limits<double>::infinity()};
  const int g0[2] = {0, 0};
  tk_ieq(t, "yoypx.sort", clab_yoypx(v, pu, g0, 2, 1, prev, chg, pct), CLAB_EARG);
  tk_ieq(t, "yoypx.lag", clab_yoypx(v, prd, g, 6, 0, prev, chg, pct), CLAB_EDOM);
  tk_ieq(t, "yoypx.inf", clab_yoypx(v, pi, g, 1, 1, prev, chg, pct), CLAB_EDOM);
  tk_ieq(t, "yoypx.n", clab_yoypx(v, prd, g, -1, 1, prev, chg, pct), CLAB_EARG);
}

void adj(tk &t) {
  const double p[4] = {0.01, 0.04, 0.03, 0.2};
  double q[4];
  int o[4];
  // as p.adjust(c(0.01, 0.04, 0.03, 0.2), method) gives them
  const double want[6][4] = {
      {0.01, 0.04, 0.03, 0.2},
      {0.04, 0.16, 0.12, 0.8},
      {0.04, 0.09, 0.09, 0.2},
      {0.04, 0.08, 0.08, 0.2},
      {0.04, 0.16 / 3, 0.16 / 3, 0.2},
      {0.04 * 25 / 12, 0.16 / 3 * 25 / 12, 0.16 / 3 * 25 / 12, 0.2 * 25 / 12}};
  for (int m = 0; m < 6; ++m) {
    tk_ieq(t, "pbadj.rc", clab_pbadj(p, 4, m, q, o), CLAB_OK);
    tk_veq(t, "pbadj.q", q, want[m], 4, 1e-15);
  }
  const double pb[1] = {1.5};
  tk_ieq(t, "pbadj.meth", clab_pbadj(p, 4, 6, q, o), CLAB_EARG);
  tk_ieq(t, "pbadj.p", clab_pbadj(pb, 1, 0, q, o), CLAB_EDOM);
}

}  // namespace

void t_pub(tk &t, int big) {
  (void)big;
  bands(t);
  cells(t);
  yoy(t);
  adj(t);
}
