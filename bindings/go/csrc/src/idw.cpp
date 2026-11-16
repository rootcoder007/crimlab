// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Inverse-distance-weighted interpolation (crimlab.h, idwint()): the
// value of monitors carried to target points as the mean of the k
// nearest within a radius, weighted by distance to the power -p, and the
// mean of the targets of each zone (the zonal statistic of a surface
// sampled at the zone's points). Planar or great-circle distances. Work:
// O(targets x monitors log k).
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

namespace {

using std::size_t;
constexpr double rdn = 3.14159265358979323846 / 180.0;

inline size_t z(int i) { return static_cast<size_t>(i); }

double dist(double x0, double y0, double x1, double y1, int ll) {
  if (!ll) return std::hypot(x1 - x0, y1 - y0);
  const double p0 = y0 * rdn, p1 = y1 * rdn;
  const double a = std::sin((p1 - p0) / 2), c = std::sin((x1 - x0) * rdn / 2);
  const double h = a * a + std::cos(p0) * std::cos(p1) * c * c;
  return 2.0 * 6371.0088 * std::asin(std::sqrt(std::min(1.0, h)));
}

// The value at one target: the monitors at its location (their mean),
// else the k nearest within rad by distance and index, weighted by d^-p.
double one(int nm, const double *mx, const double *my, const double *mv, double tx, double ty,
           double p, int k, double rad, int ll, std::vector<std::pair<double, int>> &dv,
           double &used) {
  dv.clear();
  double s0 = 0.0;
  int n0 = 0;
  for (int m = 0; m < nm; ++m) {
    if (!std::isfinite(mv[m])) continue;
    const double d = dist(mx[m], my[m], tx, ty, ll);
    if (d == 0.0) {
      s0 += mv[m];
      ++n0;
    } else if (d <= rad) {
      dv.emplace_back(d, m);
    }
  }
  if (n0 > 0) {
    used = n0;
    return s0 / n0;
  }
  const size_t kk = k > 0 ? std::min(z(k), dv.size()) : dv.size();
  std::partial_sort(dv.begin(), dv.begin() + static_cast<std::ptrdiff_t>(kk), dv.end());
  double sw = 0.0, sv = 0.0;
  for (size_t j = 0; j < kk; ++j) {
    const double w = std::pow(dv[j].first, -p);
    sw += w;
    sv += w * mv[dv[j].second];
  }
  used = static_cast<double>(kk);
  return kk > 0 ? sv / sw : NAN;
}

}  // namespace

extern "C" int clab_idwint(int nm, const double *mx, const double *my, const double *mv, int nt,
                           const double *tx, const double *ty, double p, int k, double rad,
                           int ll, const int *zone, int nz, double *out, double *used,
                           double *zout, double *zn) CLAB_TRY {
  if (nm < 1 || nt < 1 || !(p >= 0.0) || !std::isfinite(p) || k < 0 || !(rad > 0.0) ||
      (zone && nz < 1)) {
    return CLAB_EARG;
  }
  for (int m = 0; m < nm; ++m) {
    if (!std::isfinite(mx[m]) || !std::isfinite(my[m])) return CLAB_EDOM;
  }
  for (int t = 0; t < nt; ++t) {
    if (!std::isfinite(tx[t]) || !std::isfinite(ty[t])) return CLAB_EDOM;
    if (zone && (zone[t] < 0 || zone[t] >= nz)) return CLAB_EARG;
  }
  std::vector<std::pair<double, int>> dv;
  dv.reserve(z(nm));
  for (int t = 0; t < nt; ++t) {
    if (kutil::stop(static_cast<double>(nm) * 4.0)) return CLAB_EINTR;
    out[t] = one(nm, mx, my, mv, tx[t], ty[t], p, k, rad, ll, dv, used[t]);
  }
  if (zone) {
    std::fill(zout, zout + nz, 0.0);
    std::fill(zn, zn + nz, 0.0);
    for (int t = 0; t < nt; ++t) {
      if (std::isnan(out[t])) continue;
      zout[zone[t]] += out[t];
      zn[zone[t]] += 1.0;
    }
    for (int g = 0; g < nz; ++g) zout[g] = zn[g] > 0 ? zout[g] / zn[g] : NAN;
  }
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)
