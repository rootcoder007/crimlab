// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests of the longitude/latitude kernels (src/geo.cpp).
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "cgeo.h"
#include "crimlab.h"
#include "tchk.h"

namespace {

const double R = 6371008.8;
const double PI = cgeo::pi;
double gr(int i, int a, int m) {
  return static_cast<double>((i * a) % m) / m;
}

double d1(double lo1, double la1, double lo2, double la2) {
  double d = -1.0;
  clab_hvdist(&lo1, &la1, &lo2, &la2, 1, R, &d);
  return d;
}

// independent route: unit vectors in R^3, angle = 2 asin(chord / 2)
double chord(double lo1, double la1, double lo2, double la2) {
  const double k = PI / 180.0;
  const double u[3] = {std::cos(la1 * k) * std::cos(lo1 * k),
                       std::cos(la1 * k) * std::sin(lo1 * k),
                       std::sin(la1 * k)};
  const double v[3] = {std::cos(la2 * k) * std::cos(lo2 * k),
                       std::cos(la2 * k) * std::sin(lo2 * k),
                       std::sin(la2 * k)};
  const double c = std::sqrt((u[0] - v[0]) * (u[0] - v[0]) +
                             (u[1] - v[1]) * (u[1] - v[1]) +
                             (u[2] - v[2]) * (u[2] - v[2]));
  return R * 2.0 * std::asin(0.5 * c);
}

void famtest(tk &T, int /*big*/) {
  // geosphere 1.5.20, distVincentySphere(p1, p2, r = 6371008.8), pairs
  // from dev/ref/cpp/cgr.R (R 4.5.0)
  const double gcd[50] = {
      555975.40116766398, 18876573.219166435, 11103759.137649018,
      17932774.037414499, 5158465.9279094785, 18839267.585578032,
      8992368.0914049074, 13991894.930613283, 15806057.435283799,
      9523083.1548771579, 17055588.394189652, 10369935.760833401,
      14386365.973060286, 16581627.140775314, 15124959.470544675,
      13950679.16382076,  16147723.831516633, 5312363.5927596251,
      4385577.3617087789, 3677530.8047800669, 15761230.268278403,
      12082713.516263362, 13217250.360789808, 2635069.0855921772,
      11542223.991292814, 11877426.703435924, 13844460.911190787,
      9742024.3571423199, 4115673.8853848297, 11874065.631591398,
      9153749.7541164216, 6846071.3613570295, 2728847.7032401212,
      7695859.9320466369, 8881073.3601547033, 12676537.943345433,
      2743683.7557668127, 11667866.399634887, 6197298.9226782043,
      11249438.675121836, 6254621.3378049908, 7677031.0256197471,
      7011424.1792598134, 9015640.7389638834, 8750626.0756016616,
      6188127.3200290427, 8605691.4249448199, 6795574.8388367491,
      10956607.354599977, 6021300.5682655675};
  double lo1[50], la1[50], lo2[50], la2[50], d[50];
  for (int k = 0; k < 50; ++k) {
    lo1[k] = gr(k, 7919, 1009) * 360 - 180;
    la1[k] = gr(k, 104729, 1013) * 180 - 90;
    lo2[k] = gr(k, 31, 97) * 360 - 180;
    la2[k] = gr(k, 61, 101) * 170 - 85;
  }
  tk_true(T, "geo.ok", clab_hvdist(lo1, la1, lo2, la2, 50, R, d) == CLAB_OK);
  for (int k = 0; k < 50; ++k) tk_eq(T, "geo.vincenty", d[k], gcd[k], 1e-13);
  // the R^3 chord route (away from antipodes, where it is well conditioned)
  for (int k = 0; k < 50; ++k)
    if (d[k] < 0.9 * PI * R)
      tk_eq(T, "geo.chord", d[k], chord(lo1[k], la1[k], lo2[k], la2[k]), 1e-12);

  // closed forms
  tk_eq(T, "geo.deg", d1(0, 0, 1, 0), R * PI / 180, 1e-15);
  tk_eq(T, "geo.antipode", d1(0, 0, 180, 0), R * PI, 1e-15);
  tk_eq(T, "geo.antipode2", d1(10, 20, -170, -20), R * PI, 1e-15);
  tk_eq(T, "geo.pole", d1(37, 90, -150, -90), R * PI, 1e-15);
  tk_eq(T, "geo.meridian", d1(5, -30, 5, 45), R * 75 * PI / 180, 1e-15);
  tk_true(T, "geo.same", d1(12.5, -33.25, 12.5, -33.25) == 0.0);
  tk_eq(T, "geo.wrap", d1(179.9, 0, -179.9, 0), R * 0.2 * PI / 180, 1e-11);
  // a nanodegree: relative accuracy kept far below a millimetre
  tk_eq(T, "geo.tiny", d1(0, 0, 0, 1e-9), R * 1e-9 * PI / 180, 1e-12);
  {
    const double a = -87.6, b = a + 1e-9, dl = b - a;  // dl exact
    tk_eq(T, "geo.tinylon", d1(a, 41.9, b, 41.9),
        R * dl * PI / 180 * std::cos(41.9 * PI / 180), 1e-9);
  }
  // nearly antipodal: pi R less a small arc
  tk_eq(T, "geo.nearanti", d1(0, 0, 180, 1e-7), R * (PI - 1e-7 * PI / 180),
      1e-15);
  tk_eq(T, "geo.sym", d1(lo1[3], la1[3], lo2[3], la2[3]),
      d1(lo2[3], la2[3], lo1[3], la1[3]), 1e-15);

  // equal-area projection: the centre goes to the origin, a point at
  // central angle c goes to radius 2 R sin(c / 2) whatever its azimuth
  {
    const double lon0 = -87.63, lat0 = 41.88;
    double lo[40], la[40], x[40], y[40];
    for (int k = 0; k < 40; ++k) {
      lo[k] = lon0 + (gr(k, 7919, 1009) - 0.5) * 30;
      la[k] = lat0 + (gr(k, 104729, 1013) - 0.5) * 30;
    }
    tk_true(T, "eqp.ok", clab_eqproj(lo, la, 40, lon0, lat0, R, x, y) == CLAB_OK);
    for (int k = 0; k < 40; ++k) {
      const double c = d1(lon0, lat0, lo[k], la[k]) / R;
      tk_eq(T, "eqp.radius", std::hypot(x[k], y[k]), 2 * R * std::sin(c / 2),
          1e-12);
    }
    double x0, y0;
    clab_eqproj(&lon0, &lat0, 1, lon0, lat0, R, &x0, &y0);
    tk_true(T, "eqp.centre", x0 == 0.0 && y0 == 0.0);
    // the area element |J| equals R^2 cos(lat) per square radian at a
    // spread of places (central differences, step 1e-4 degrees)
    for (int k = 0; k < 40; k += 4) {
      const double h = 1e-4;
      double a[4] = {lo[k] + h, lo[k] - h, lo[k], lo[k]};
      double b[4] = {la[k], la[k], la[k] + h, la[k] - h};
      double px[4], py[4];
      clab_eqproj(a, b, 4, lon0, lat0, R, px, py);
      const double rad = 2 * h * PI / 180;
      const double j = ((px[0] - px[1]) * (py[2] - py[3]) -
                        (py[0] - py[1]) * (px[2] - px[3])) /
                       (rad * rad);
      tk_eq(T, "eqp.area", j, R * R * std::cos(la[k] * PI / 180), 1e-6);
    }
    // north is up on the central meridian, east is +x on the parallel
    double lo2b[2] = {lon0, lon0 + 0.01}, la2b[2] = {lat0 + 0.01, lat0};
    double qx[2], qy[2];
    clab_eqproj(lo2b, la2b, 2, lon0, lat0, R, qx, qy);
    tk_true(T, "eqp.north", std::fabs(qx[0]) < 1e-9 && qy[0] > 0);
    tk_true(T, "eqp.east", qx[1] > 0);
  }

  // The unit contract. Three same-day incidents in lon/lat: A and B about
  // 84 m apart, C about 41 km away. Planar distances in degrees are all
  // below 1, so a "1000 m" link threshold applied to raw degrees joins all
  // three; projected first, only A and B are within 1000 m.
  {
    double lo[3] = {-87.6300, -87.6290, -87.1350};
    double la[3] = {41.8800, 41.8803, 41.8800};
    double x[3], y[3];
    clab_eqproj(lo, la, 3, -87.4, 41.88, R, x, y);
    const double ab = std::hypot(x[1] - x[0], y[1] - y[0]);
    const double ac = std::hypot(x[2] - x[0], y[2] - y[0]);
    tk_eq(T, "unit.ab", ab, d1(lo[0], la[0], lo[1], la[1]), 1e-5);
    tk_eq(T, "unit.ac", ac, d1(lo[0], la[0], lo[2], la[2]), 1e-5);
    tk_true(T, "unit.split", ab < 1000 && ac > 40000);
    tk_true(T, "unit.degrees", std::hypot(lo[2] - lo[0], la[2] - la[0]) < 1000);
  }

  // guards
  double z = 0, bad = 91, nanv = std::nan(""), dd, xx, yy;
  tk_true(T, "geo.lat", clab_hvdist(&z, &bad, &z, &z, 1, R, &dd) == CLAB_EDOM);
  tk_true(T, "geo.nan", clab_hvdist(&nanv, &z, &z, &z, 1, R, &dd) == CLAB_EDOM);
  tk_true(T, "geo.r", clab_hvdist(&z, &z, &z, &z, 1, 0.0, &dd) == CLAB_EDOM);
  tk_true(T, "geo.n", clab_hvdist(&z, &z, &z, &z, -1, R, &dd) == CLAB_EARG);
  tk_true(T, "geo.n0", clab_hvdist(&z, &z, &z, &z, 0, R, &dd) == CLAB_OK);
  double alo = 180, ala = 0;
  tk_true(T, "eqp.antipode",
     clab_eqproj(&alo, &ala, 1, 0, 0, R, &xx, &yy) == CLAB_EDOM);
  tk_true(T, "eqp.lat", clab_eqproj(&z, &bad, 1, 0, 0, R, &xx, &yy) == CLAB_EDOM);
  tk_true(T, "eqp.lat0", clab_eqproj(&z, &z, 1, 0, 95, R, &xx, &yy) == CLAB_EDOM);
}

#ifdef CLAB_MAIN
// a heap buffer for the timings (kept out of the program's static size)
template <class T>
struct hb {
  T *p;
  explicit hb(long n) : p(static_cast<T *>(std::malloc(sizeof(T) * n))) {
    if (p == nullptr) std::abort();
  }
  ~hb() { std::free(p); }
  hb(const hb &) = delete;
  hb &operator=(const hb &) = delete;
};

double med5(double *v) { return tk_q7(v, 5, 0.5); }

void fambench() {
  const int n = 1000000;
  hb<double> lo1_(n), la1_(n), lo2_(n), la2_(n), d_(n);

  double *lo1 = lo1_.p, *la1 = la1_.p, *lo2 = lo2_.p, *la2 = la2_.p, *d = d_.p;
  cgeo::rng r(7);
  for (int i = 0; i < n; ++i) {
    lo1[i] = r.unif() * 360 - 180;
    la1[i] = r.unif() * 180 - 90;
    lo2[i] = r.unif() * 360 - 180;
    la2[i] = r.unif() * 180 - 90;
  }
  double t[5];
  for (int k = 0; k < 5; ++k) {
    const auto a = std::chrono::steady_clock::now();
    clab_hvdist(lo1, la1, lo2, la2, n, R, d);
    t[k] = std::chrono::duration<double>(std::chrono::steady_clock::now() - a)
               .count();
  }
  std::printf("bench geo: hvdist n = 1e6 %.4f s (median of 5)\n", med5(t));
  for (int k = 0; k < 5; ++k) {
    const auto a = std::chrono::steady_clock::now();
    clab_eqproj(lo1, la1, n, 0.0, 0.0, R, lo2, la2);
    t[k] = std::chrono::duration<double>(std::chrono::steady_clock::now() - a)
               .count();
    for (int i = 0; i < n; ++i) la1[i] = 0.5 * la1[i];
  }
  std::printf("bench geo: eqproj n = 1e6 %.4f s (median of 5)\n", med5(t));
}

#endif

}  // namespace

void t_geo(tk &T, int big) { famtest(T, big); }
#ifdef CLAB_MAIN
void b_geo() { fambench(); }
#endif
