// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for the downwind exposure kernel (dwn.cpp): bearings
// that are whole angles or arctangents of simple ratios, a wedge met only
// by its arc, the apex inside a polygon, the angle rule with and without
// the wrap about north, and the failure paths. Part of the self-test
// (selftest.cpp).
#include <cmath>

#include "crimlab.h"
#include "tchk.h"

namespace {

void sector(tk &tt) {
  // one source at the origin, the wind blowing to the north, 30 degrees
  // either side: (0, 100) at 0, (50, 100) at atan(0.5) = 26.57 degrees in;
  // (60, 100) at 30.96 and (100, 10) at 84.3 degrees out; the source's
  // own point in; with the radius 105, (50, 100) at 111.8 drops out
  const double sx[1] = {0}, sy[1] = {0}, to[1] = {0};
  const int sp[1] = {0};
  const double rx[5] = {0, 50, 60, 100, 0}, ry[5] = {100, 100, 100, 10, 0};
  double o[5];
  tk_ieq(tt, "dwn.sec.rc", clab_dwnexp(0, 1, sx, sy, sp, to, 5, rx, ry, nullptr, 1, 30, INFINITY,
                                       0, 1, o, nullptr),
         CLAB_OK);
  const double want[5] = {1, 1, 0, 0, 1};
  for (int i = 0; i < 5; ++i) tk_near(tt, "dwn.sec", o[i], want[i], 0);
  clab_dwnexp(0, 1, sx, sy, sp, to, 5, rx, ry, nullptr, 1, 30, 105, 0, 1, o, nullptr);
  tk_near(tt, "dwn.sec.rad", o[1], 0, 0);
  // on the sphere from (0, 0): (0, 10) due north, (10, 0) due east
  const double lx[2] = {0, 10}, ly[2] = {10, 0};
  clab_dwnexp(0, 1, sx, sy, sp, to, 2, lx, ly, nullptr, 1, 30, INFINITY, 1, 1, o, nullptr);
  tk_near(tt, "dwn.sec.ll0", o[0], 1, 0);
  tk_near(tt, "dwn.sec.ll1", o[1], 0, 0);
}

void wedge(tk &tt) {
  // the wedge north, 30 degrees, radius 100; four squares and a strip
  const double sx[1] = {0}, sy[1] = {0}, to[1] = {0};
  const int sp[1] = {0};
  // a square with vertices inside; one far to the east; one about the
  // apex; a strip y in [95, 99], x in [-80, 80]: no vertex within 100 (the
  // nearest is 124 away), no side reached (they end at y = 86.6), but the
  // arc crosses y = 95 at x = 31.2, 18.2 degrees off north; a strip just
  // above the arc, y in [100.5, 150]
  const double rx[20] = {-10, 10, 10, -10, 50, 70, 70, 50, -5, 5, 5, -5,
                         -80, 80, 80, -80, -10, 10, 10, -10};
  const double ry[20] = {40, 40, 60, 60, -5, -5, 5, 5, -5, -5, 5, 5,
                         95, 95, 99, 99, 100.5, 100.5, 150, 150};
  const int ro[6] = {0, 4, 8, 12, 16, 20};
  double o[5];
  tk_ieq(tt, "dwn.wdg.rc",
         clab_dwnexp(1, 1, sx, sy, sp, to, 5, rx, ry, ro, 1, 30, 100, 0, 1, o, nullptr),
         CLAB_OK);
  const double want[5] = {1, 0, 1, 1, 0};
  for (int i = 0; i < 5; ++i) tk_near(tt, "dwn.wdg", o[i], want[i], 0);
  // without the radius the last strip has vertices within the sector
  clab_dwnexp(1, 1, sx, sy, sp, to, 5, rx, ry, ro, 1, 30, INFINITY, 0, 1, o, nullptr);
  tk_near(tt, "dwn.wdg.inf", o[4], 1, 0);
  tk_ieq(tt, "dwn.wdg.ll",
         clab_dwnexp(1, 1, sx, sy, sp, to, 5, rx, ry, ro, 1, 30, 100, 1, 1, o, nullptr),
         CLAB_EARG);
  const int r2[2] = {0, 2};
  tk_ieq(tt, "dwn.wdg.ring",
         clab_dwnexp(1, 1, sx, sy, sp, to, 1, rx, ry, r2, 1, 30, 100, 0, 1, o, nullptr),
         CLAB_EARG);
}

void angle(tk &tt) {
  // sources at (0, 0) and (10, 0) blowing to 60 and 70 degrees: the mean
  // point (5, 0), the mean direction 65; (5, 10) lies at bearing 0, 65
  // off, (5, -10) at 180, 115 off; both sqrt(125) from each source
  const double sx[2] = {0, 10}, sy[2] = {0, 0}, to[2] = {60, 70};
  const int sp[2] = {0, 0};
  const double rx[2] = {5, 5}, ry[2] = {10, -10};
  double o[2], d[2];
  tk_ieq(tt, "dwn.ang.rc",
         clab_dwnexp(2, 2, sx, sy, sp, to, 2, rx, ry, nullptr, 1, 30, INFINITY, 0, 1, o, d),
         CLAB_OK);
  tk_near(tt, "dwn.ang.a0", o[0], 65, 1e-12);
  tk_near(tt, "dwn.ang.a1", o[1], 115, 1e-12);
  tk_near(tt, "dwn.ang.d", d[0], std::sqrt(125.0), 1e-12);
  // one source blowing to 10 degrees, a place at bearing 350 (tan 10
  // degrees to the west): 20 degrees about the circle, 340 plainly; a
  // second period without sources is NaN; a place on the source is left
  // out of its mean distance, so with no other source it is NaN
  const double s1[1] = {0}, t1[1] = {10};
  const int p1[1] = {0};
  const double qx[2] = {-100 * std::tan(10 * 3.14159265358979323846 / 180), 0}, qy[2] = {100, 0};
  double a[4], b[4];
  clab_dwnexp(2, 1, s1, s1, p1, t1, 2, qx, qy, nullptr, 2, 30, INFINITY, 0, 1, a, b);
  tk_near(tt, "dwn.ang.wrap", a[0], 20, 1e-9);
  tk_ieq(tt, "dwn.ang.empty", std::isnan(a[2]) ? 1 : 0, 1);
  tk_ieq(tt, "dwn.ang.self", std::isnan(b[1]) ? 1 : 0, 1);
  clab_dwnexp(2, 1, s1, s1, p1, t1, 2, qx, qy, nullptr, 2, 30, INFINITY, 0, 0, a, b);
  tk_near(tt, "dwn.ang.plain", a[0], 340, 1e-9);
}

void fail(tk &tt) {
  const double x[1] = {0}, nx[1] = {NAN};
  const int p[1] = {0}, bad[1] = {1};
  double o[1], d[1];
  tk_ieq(tt, "dwn.f.rule", clab_dwnexp(3, 1, x, x, p, x, 1, x, x, nullptr, 1, 30, 1, 0, 1, o, d),
         CLAB_EARG);
  tk_ieq(tt, "dwn.f.per", clab_dwnexp(0, 1, x, x, bad, x, 1, x, x, nullptr, 1, 30, 1, 0, 1, o, d),
         CLAB_EARG);
  tk_ieq(tt, "dwn.f.half", clab_dwnexp(0, 1, x, x, p, x, 1, x, x, nullptr, 1, 181, 1, 0, 1, o, d),
         CLAB_EARG);
  tk_ieq(tt, "dwn.f.rad", clab_dwnexp(0, 1, x, x, p, x, 1, x, x, nullptr, 1, 30, 0, 0, 1, o, d),
         CLAB_EARG);
  tk_ieq(tt, "dwn.f.nan", clab_dwnexp(0, 1, nx, x, p, x, 1, x, x, nullptr, 1, 30, 1, 0, 1, o, d),
         CLAB_EDOM);
  tk_ieq(tt, "dwn.f.out2",
         clab_dwnexp(2, 1, x, x, p, x, 1, x, x, nullptr, 1, 30, 1, 0, 1, o, nullptr), CLAB_EARG);
}

}  // namespace

void t_dwn(tk &t, int big) {
  (void)big;
  sector(t);
  wedge(t);
  angle(t);
  fail(t);
}
