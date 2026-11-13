// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, place-based designs: the treated, comparison and excluded
// places around sites.
#include <cmath>

#include "crimlab.h"
#include "kutil.h"

namespace {

// the nearest of m points to (x, y): its squared distance
double near2(double x, double y, const double *qx, const double *qy, int m) {
  double b = HUGE_VAL;
  for (int k = 0; k < m; ++k) {
    const double d = (x - qx[k]) * (x - qx[k]) + (y - qy[k]) * (y - qy[k]);
    if (d < b) b = d;
  }
  return b;
}

}  // namespace

int clab_bufas(const double *x, const double *y, int n, const double *sx, const double *sy, int ns,
               const double *cx, const double *cy, int nc, double r, double buf, int *grp,
               double *ds, double *dc) {
  if (n < 1 || ns < 1 || nc < 1) return CLAB_EARG;
  if (!kutil::pos(r) || !kutil::nonneg(buf)) return CLAB_EDOM;
  for (int i = 0; i < n; ++i) {
    ds[i] = std::sqrt(near2(x[i], y[i], sx, sy, ns));
    dc[i] = std::sqrt(near2(x[i], y[i], cx, cy, nc));
    if (!std::isfinite(ds[i]) || !std::isfinite(dc[i])) return CLAB_EDOM;
    // treated within r of a site; comparison within r of a comparison point
    // and beyond buf of every site; excluded otherwise
    grp[i] = ds[i] <= r ? 0 : (dc[i] <= r && ds[i] > buf ? 1 : 2);
  }
  return CLAB_OK;
}
