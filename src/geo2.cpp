// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, spatial description: the inverse of the equal-area
// projection of src/geo.cpp, so that centres computed in planar metres can
// be reported in longitude and latitude.
#include <cmath>

#include "cgeo.h"
#include "crimlab.h"
#include "kutil.h"

namespace {

constexpr double deg = cgeo::pi / 180.0;

bool lonlat(double lon, double lat) {
  return std::isfinite(lon) && lat >= -90.0 && lat <= 90.0;
}

}  // namespace

extern "C" {

int clab_eqinv(const double *x, const double *y, int n, double lon0,
               double lat0, double r, double *lon, double *lat) {
  if (n < 0) return CLAB_EARG;
  if (!kutil::pos(r) || !lonlat(lon0, lat0)) return CLAB_EDOM;
  for (int i = 0; i < n; ++i) {
    if (std::isnan(x[i]) || std::isnan(y[i])) continue;
    if (!std::isfinite(x[i]) || !std::isfinite(y[i])) return CLAB_EDOM;
    const double u = x[i] / (2.0 * r), v = y[i] / (2.0 * r);
    if (u * u + v * v > 1.0) return CLAB_EDOM;
  }
  const double p0 = lat0 * deg, s0 = std::sin(p0), c0 = std::cos(p0);
  for (int i = 0; i < n; ++i) {
    if (std::isnan(x[i]) || std::isnan(y[i])) {
      lon[i] = lat[i] = kutil::nan();
      continue;
    }
    // with s = rho / (2 r) and c the central angle (s = sin(c / 2)): cos c
    // = 1 - 2 s^2 and sin c (x, y) / rho = (x, y) sqrt(1 - s^2) / r, so the
    // point's unit vector, in the frame of the centre at longitude 0 (east,
    // north, up), needs no division by rho and no inverse sine; latitude
    // and longitude come from atan2, well conditioned everywhere (Snyder
    // 1987, doi:10.3133/pp1395, ch. 24, rewritten)
    const double u = x[i] / (2.0 * r), v = y[i] / (2.0 * r);
    const double s2 = u * u + v * v;
    const double q = std::sqrt(1.0 - s2), cc = 1.0 - 2.0 * s2;
    const double a = x[i] * q / r, b = y[i] * q / r;
    const double vx = cc * c0 - b * s0, vy = a, vz = cc * s0 + b * c0;
    const double la = std::atan2(vz, std::hypot(vx, vy)) / deg;
    double lo = lon0 + std::atan2(vy, vx) / deg;
    if (lo > 180.0) lo -= 360.0;
    if (lo < -180.0) lo += 360.0;
    lon[i] = lo;
    lat[i] = la;
  }
  return CLAB_OK;
}

}  // extern "C"
