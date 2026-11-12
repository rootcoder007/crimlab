// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of the inverse equal-area projection (src/geo2.cpp).
#include <cmath>

#include "cgeo.h"
#include "crimlab.h"
#include "tchk.h"

namespace {

const double R = 6371008.8;

double gr(int i, int a, int m) { return static_cast<double>((i * a) % m) / m; }

double dist(double lo1, double la1, double lo2, double la2) {
  double d = -1.0;
  clab_hvdist(&lo1, &la1, &lo2, &la2, 1, R, &d);
  return d;
}

void famtest(tk &T, int big) {
  // the round trip, centres and points spread over the globe, short of
  // the antipode (central angle below 170 degrees): back within ten
  // micrometres on the ground (SP6.0; the worst of the extended run is
  // about 1.1e-6 m, near 170 degrees)
  const int nc = big ? 200 : 40;
  double worst = 0.0;
  for (int c = 0; c < nc; ++c) {
    const double lon0 = gr(c, 7919, 1009) * 360 - 180;
    const double lat0 = gr(c, 104729, 1013) * 180 - 90;
    double lo[50], la[50], x[50], y[50], lo2[50], la2[50];
    int m = 0;
    for (int k = 0; k < 50; ++k) {
      const double a = gr(k + 50 * c, 31, 97) * 360 - 180;
      const double b = gr(k + 50 * c, 61, 101) * 180 - 90;
      if (dist(lon0, lat0, a, b) < R * 170 * cgeo::pi / 180) {
        lo[m] = a;
        la[m] = b;
        ++m;
      }
    }
    clab_eqproj(lo, la, m, lon0, lat0, R, x, y);
    tk_true(T, "eqinv.ok", clab_eqinv(x, y, m, lon0, lat0, R, lo2, la2) == CLAB_OK);
    for (int k = 0; k < m; ++k) {
      const double e = dist(lo[k], la[k], lo2[k], la2[k]);
      if (e > worst) worst = e;
    }
  }
  tk_lt(T, "eqinv.round", worst, 1e-5);
  // the poles and the centre itself (SP6.2); the longitude of a pole is
  // undefined, so only its latitude is checked
  {
    double lo[3] = {10, -170, 33}, la[3] = {90, -90, 45}, x[3], y[3], a[3], b[3];
    clab_eqproj(lo, la, 3, 33, 45, R, x, y);
    clab_eqinv(x, y, 3, 33, 45, R, a, b);
    tk_eq(T, "eqinv.npole", b[0], 90, 1e-12);
    tk_eq(T, "eqinv.spole", b[1], -90, 1e-12);
    tk_true(T, "eqinv.centre", a[2] == 33 && b[2] == 45);
    double z = 0;
    clab_eqinv(&z, &z, 1, 0, 90, R, a, b);
    tk_true(T, "eqinv.polecentre", b[0] == 90);
  }
  // across the date line, and the edge of the disc (the antipode, 2 R)
  {
    double lo = 179.5, la = -10, x, y, a, b;
    clab_eqproj(&lo, &la, 1, -179.5, -10, R, &x, &y);
    clab_eqinv(&x, &y, 1, -179.5, -10, R, &a, &b);
    tk_eq(T, "eqinv.dateline", a, 179.5, 1e-9);
    double ex = 2 * R, ey = 0;
    tk_true(T, "eqinv.edge", clab_eqinv(&ex, &ey, 1, 0, 0, R, &a, &b) == CLAB_OK);
    tk_eq(T, "eqinv.edgelat", b, 0, 1e-12);
    tk_eq(T, "eqinv.edgelon", std::fabs(a), 180, 1e-12);
  }
  // NaN passes through; aliasing the output with the input is allowed
  {
    double x[2] = {kutil::nan(), 1000}, y[2] = {0, -2000};
    tk_true(T, "eqinv.alias", clab_eqinv(x, y, 2, 5, 50, R, x, y) == CLAB_OK);
    tk_nan(T, "eqinv.nan", x[0]);
    tk_true(T, "eqinv.nan2", std::isnan(y[0]) && x[1] > 5 && y[1] < 50);
  }
  // guards
  {
    double x = 1, y = 1, a, b, big2 = 2.0001 * R, inf = HUGE_VAL;
    tk_ieq(T, "eqinv.n", clab_eqinv(&x, &y, -1, 0, 0, R, &a, &b), CLAB_EARG);
    tk_ieq(T, "eqinv.r", clab_eqinv(&x, &y, 1, 0, 0, 0, &a, &b), CLAB_EDOM);
    tk_ieq(T, "eqinv.lat0", clab_eqinv(&x, &y, 1, 0, 91, R, &a, &b), CLAB_EDOM);
    tk_ieq(T, "eqinv.far", clab_eqinv(&big2, &y, 1, 0, 0, R, &a, &b), CLAB_EDOM);
    tk_ieq(T, "eqinv.inf", clab_eqinv(&inf, &y, 1, 0, 0, R, &a, &b), CLAB_EDOM);
    tk_ieq(T, "eqinv.zero", clab_eqinv(&x, &y, 0, 0, 0, R, &a, &b), CLAB_OK);
  }
}

}  // namespace

void t_geo2(tk &T, int big) { famtest(T, big); }
