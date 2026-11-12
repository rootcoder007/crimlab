// SPDX-License-Identifier: AGPL-3.0-or-later
//
// libcrimlab, spatial description: longitude and latitude. Great-circle
// distances and a local equal-area projection, so that the planar kernels
// (centrography, distances, hot spots, journey to crime) receive metres
// and never degrees.
#include <cmath>

#include "cgeo.h"
#include "crimlab.h"

namespace {

constexpr double deg = cgeo::pi / 180.0;

bool lonlat(double lon, double lat) {
  return std::isfinite(lon) && lat >= -90.0 && lat <= 90.0;
}

// a longitude difference reduced exactly to [-180, 180] degrees
double dlon(double a, double b) {
  double d = std::fmod(a - b, 360.0);
  if (d > 180.0) d -= 360.0;
  if (d < -180.0) d += 360.0;
  return d;
}

}  // namespace

extern "C" {

int clab_hvdist(const double *lon1, const double *lat1, const double *lon2,
                const double *lat2, int n, double r, double *d) {
  if (n < 0) return CLAB_EARG;
  if (!(r > 0.0) || !std::isfinite(r)) return CLAB_EDOM;
  for (int i = 0; i < n; ++i)
    if (!lonlat(lon1[i], lat1[i]) || !lonlat(lon2[i], lat2[i]))
      return CLAB_EDOM;
  for (int i = 0; i < n; ++i) {
    const double p1 = lat1[i] * deg, p2 = lat2[i] * deg;
    const double dp = (lat2[i] - lat1[i]) * deg;
    const double l = dlon(lon2[i], lon1[i]) * deg;
    const double s1 = std::sin(p1), c1 = std::cos(p1), c2 = std::cos(p2);
    const double h = std::sin(0.5 * l), hh = 2.0 * h * h;  // 1 - cos l
    // the central angle as atan2(|u x v|, u . v), with
    //   c1 s2 - s1 c2 cos l = sin(p2 - p1) + s1 c2 (1 - cos l) and
    //   s1 s2 + c1 c2 cos l = cos(p2 - p1) - c1 c2 (1 - cos l),
    // so that nothing cancels: well conditioned from coincident points to
    // antipodes
    const double a = c2 * std::sin(l);
    const double b = std::sin(dp) + s1 * c2 * hh;
    const double num = std::hypot(a, b);
    const double den = std::cos(dp) - c1 * c2 * hh;
    d[i] = r * std::atan2(num, den);
  }
  return CLAB_OK;
}

int clab_eqproj(const double *lon, const double *lat, int n, double lon0,
                double lat0, double r, double *x, double *y) {
  if (n < 0) return CLAB_EARG;
  if (!(r > 0.0) || !std::isfinite(r) || !lonlat(lon0, lat0))
    return CLAB_EDOM;
  for (int i = 0; i < n; ++i)
    if (!lonlat(lon[i], lat[i])) return CLAB_EDOM;
  const double p0 = lat0 * deg, s0 = std::sin(p0), c0 = std::cos(p0);
  for (int i = 0; i < n; ++i) {
    const double p = lat[i] * deg, l = dlon(lon[i], lon0) * deg;
    const double dp = (lat[i] - lat0) * deg;
    const double c = std::cos(p);
    const double h = std::sin(0.5 * l), hh = 2.0 * h * h;  // 1 - cos l
    // Lambert azimuthal equal-area on the sphere (Snyder 1987,
    // doi:10.3133/pp1395, eqs. 24-2 to 24-4), with
    //   c0 s - s0 c cos l = sin(p - p0) + s0 c (1 - cos l) and
    //   1 + s0 s + c0 c cos l = 1 + cos(p - p0) - c0 c (1 - cos l)
    const double q = 1.0 + std::cos(dp) - c0 * c * hh;
    if (!(q > 1e-12)) return CLAB_EDOM;
    const double kp = std::sqrt(2.0 / q);
    x[i] = r * kp * c * std::sin(l);
    y[i] = r * kp * (std::sin(dp) + s0 * c * hh);
  }
  return CLAB_OK;
}

}  // extern "C"
