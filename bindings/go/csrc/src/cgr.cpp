// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, spatial description: centrographic statistics of point
// patterns in planar coordinates (mean, median and minimum-distance
// centres, dispersion, the standard deviational ellipse, the convex hull).
// Coordinates must be projected (clab_eqproj turns longitude and latitude
// into planar metres); none of these kernels can tell degrees from metres.
#include <cmath>
#include <limits>

#include "cgeo.h"
#include "crimlab.h"
#include "kutil.h"

extern "C" {

int clab_cgstat(const double *x, const double *y, const double *w, int n,
                int k, double *out) {
  if (n < 1 || k < 0 || k > 2) return CLAB_EARG;
  double m[6];
  if (!cgeo::moments(x, y, w, n, m)) return CLAB_EDOM;
  const double W = m[0], den = W - k;
  if (!(den > 0.0)) return CLAB_EDOM;
  out[0] = W;
  out[1] = m[1];
  out[2] = m[2];
  out[3] = std::sqrt(m[3] / den);
  out[4] = std::sqrt(m[4] / den);
  out[5] = std::sqrt((m[3] + m[4]) / den);
  // geometric and harmonic means of each coordinate: defined only when
  // every coordinate carrying weight is positive
  bool px = true, py = true;
  kutil::ksum lx, ly, ix, iy;
  for (int i = 0; i < n; ++i) {
    const double wi = w ? w[i] : 1.0;
    if (wi == 0.0) continue;
    if (x[i] > 0.0) {
      lx.add(wi * std::log(x[i]));
      ix.add(wi / x[i]);
    } else {
      px = false;
    }
    if (y[i] > 0.0) {
      ly.add(wi * std::log(y[i]));
      iy.add(wi / y[i]);
    } else {
      py = false;
    }
  }
  const double nan = kutil::nan();
  out[6] = px ? std::exp(lx.val() / W) : nan;
  out[7] = py ? std::exp(ly.val() / W) : nan;
  out[8] = px ? W / ix.val() : nan;
  out[9] = py ? W / iy.val() : nan;
  return CLAB_OK;
}

int clab_cgmedn(const double *x, const double *y, int n, double *work,
                double *out) {
  if (n < 1) return CLAB_EARG;
  for (int i = 0; i < n; ++i)
    if (!std::isfinite(x[i]) || !std::isfinite(y[i])) return CLAB_EDOM;
  const double *c[2] = {x, y};
  for (int j = 0; j < 2; ++j) {
    for (int i = 0; i < n; ++i) work[i] = c[j][i];
    cgeo::hsort(work, n, [](double a, double b) { return a < b; });
    // the mean of the two middle values when n is even
    out[j] = (n % 2) ? work[n / 2]
                     : work[n / 2 - 1] + 0.5 * (work[n / 2] - work[n / 2 - 1]);
  }
  return CLAB_OK;
}

int clab_webpt(const double *x, const double *y, const double *w, int n,
               const double *start, double tol, int maxit, double *out) {
  if (n < 1 || maxit < 0) return CLAB_EARG;
  if (!(tol > 0.0) || !std::isfinite(tol)) return CLAB_EDOM;
  double m[6];
  if (!cgeo::moments(x, y, w, n, m)) return CLAB_EDOM;
  if (start != nullptr && (!std::isfinite(start[0]) || !std::isfinite(start[1])))
    return CLAB_EDOM;
  // work relative to the weighted mean centre; steps are measured against
  // the root mean squared distance to it
  const double ox = m[1], oy = m[2];
  double scale = std::sqrt((m[3] + m[4]) / m[0]);
  if (!(scale > 0.0)) scale = 1.0;
  double cx = start ? start[0] - ox : 0.0, cy = start ? start[1] - oy : 0.0;
  double r[4];
  cgeo::weber(x, y, w, n, ox, oy, scale, tol, maxit, &cx, &cy, r);
  out[0] = ox + cx;
  out[1] = oy + cy;
  out[2] = r[0];
  out[3] = r[1];
  out[4] = r[2];
  out[5] = r[3];
  return CLAB_OK;
}

int clab_sdell(const double *x, const double *y, const double *w, int n, int k,
               double c, double *out) {
  if (n < 1 || k < 0 || k > 2) return CLAB_EARG;
  if (!(c > 0.0) || !std::isfinite(c)) return CLAB_EDOM;
  double m[6];
  if (!cgeo::moments(x, y, w, n, m)) return CLAB_EDOM;
  const double den = m[0] - k;
  if (!(den > 0.0)) return CLAB_EDOM;
  double l1, l2, t;
  cgeo::eig2(m[3] / den, m[4] / den, m[5] / den, &l1, &l2, &t);
  out[0] = m[1];
  out[1] = m[2];
  out[2] = c * std::sqrt(l1);
  out[3] = c * std::sqrt(l2);
  out[4] = t;
  out[5] = cgeo::pi * out[2] * out[3];
  return CLAB_OK;
}

int clab_cvhull(const double *x, const double *y, int n, int *idx, int *hull,
                int *nh, double *out) {
  if (n < 1) return CLAB_EARG;
  for (int i = 0; i < n; ++i)
    if (!cgeo::okcoord(x[i]) || !cgeo::okcoord(y[i])) return CLAB_EDOM;
  for (int i = 0; i < n; ++i) idx[i] = i;
  cgeo::hsort(idx, n, [x, y](int a, int b) {
    if (x[a] != x[b]) return x[a] < x[b];
    if (y[a] != y[b]) return y[a] < y[b];
    return a < b;
  });
  // keep one of each set of coincident points: the smallest index
  int u = 0;
  for (int i = 0; i < n; ++i) {
    if (u > 0 && x[idx[i]] == x[idx[u - 1]] && y[idx[i]] == y[idx[u - 1]])
      continue;
    idx[u++] = idx[i];
  }
  // Andrew's monotone chain (Andrew 1979,
  // doi:10.1016/0020-0190(79)90072-3): lower then upper chain, a point
  // leaves the chain unless it makes a strict left turn
  int h = 0;
  auto turn = [x, y, hull](int a, int b, int p) {
    return cgeo::orient(x[hull[a]], y[hull[a]], x[hull[b]], y[hull[b]], x[p],
                        y[p]);
  };
  for (int i = 0; i < u; ++i) {
    while (h >= 2 && turn(h - 2, h - 1, idx[i]) <= 0) --h;
    hull[h++] = idx[i];
  }
  for (int i = u - 2, lo = h + 1; i >= 0; --i) {
    while (h >= lo && turn(h - 2, h - 1, idx[i]) <= 0) --h;
    hull[h++] = idx[i];
  }
  if (u > 1) --h;  // the last point repeats the first
  *nh = h;
  // shoelace area and perimeter, relative to the first vertex
  kutil::ksum a2, per;
  const double x0 = x[hull[0]], y0 = y[hull[0]];
  for (int i = 0; i < h; ++i) {
    const int p = hull[i], q = hull[(i + 1) % h];
    const double ux = x[p] - x0, uy = y[p] - y0;
    const double vx = x[q] - x0, vy = y[q] - y0;
    a2.add(ux * vy);
    a2.add(-uy * vx);
    per.add(std::hypot(x[q] - x[p], y[q] - y[p]));
  }
  out[0] = h >= 3 ? 0.5 * a2.val() : 0.0;
  out[1] = per.val();
  return CLAB_OK;
}

}  // extern "C"
