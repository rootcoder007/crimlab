// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The layout of the maps (crimlab.h; chomap(), smmap(), ptmap() and
// prpmap()): class breaks and the class of each value, the bounding box
// of coordinates, the centres of point clusters and the radii and
// drawing order of proportional symbols. The R functions draw only what
// these compute. Work: O(n log n).
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

namespace {

using std::size_t;
inline size_t z(int i) { return static_cast<size_t>(i); }

// The box of the points with both coordinates finite, a range of one
// value widened by 0.5 on each side (as the plots' .prange() does); false
// when no point is finite.
bool box(const double *x, const double *y, int n, double *bb) {
  bool any = false;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(x[i]) || !std::isfinite(y[i])) continue;
    if (!any) {
      bb[0] = bb[1] = x[i];
      bb[2] = bb[3] = y[i];
      any = true;
    }
    bb[0] = std::min(bb[0], x[i]);
    bb[1] = std::max(bb[1], x[i]);
    bb[2] = std::min(bb[2], y[i]);
    bb[3] = std::max(bb[3], y[i]);
  }
  for (int a = 0; any && a < 4; a += 2) {
    if (bb[a + 1] - bb[a] == 0) {
      bb[a] -= 0.5;
      bb[a + 1] += 0.5;
    }
  }
  return any;
}

}  // namespace

extern "C" int clab_mapbrk(const double *v, int n, int k, int meth, double *br, int *nb, int *cls) CLAB_TRY {
  if (n < 1 || meth < 0 || meth > 2 || (meth < 2 && k < 1) || (meth == 2 && *nb < 2)) return CLAB_EARG;
  std::vector<double> f;
  for (int i = 0; i < n; ++i) {
    if (std::isfinite(v[i])) f.push_back(v[i]);
  }
  if (f.empty()) return CLAB_EDOM;
  const double lo = *std::min_element(f.begin(), f.end()), hi = *std::max_element(f.begin(), f.end());
  int m = *nb;
  if (meth == 2) {
    for (int j = 1; j < m; ++j) {
      if (!(br[j] > br[j - 1])) return CLAB_EARG;
    }
    if (br[0] > lo || br[m - 1] < hi) return CLAB_EDOM;
  } else {
    if (meth == 0) {
      // seq(lo, hi, length.out = k + 1)
      br[0] = lo;
      for (int j = 1; j < k; ++j) br[j] = lo + j * ((hi - lo) / k);
      br[k] = hi;
    } else {
      // quantile(v, seq(0, 1, length.out = k + 1), type = 7)
      std::vector<double> p(z(k) + 1, 0.0);
      for (int j = 1; j < k; ++j) p[z(j)] = j * (1.0 / k);
      p[z(k)] = 1;
      clab_q7(f.data(), static_cast<int>(f.size()), p.data(), k + 1, br);
    }
    // breaks that repeat merged; one value left is widened by 0.5
    m = 1;
    for (int j = 1; j <= k; ++j) {
      if (br[j] != br[m - 1]) br[m++] = br[j];
    }
    if (m == 1) {
      br[0] = lo - 0.5;
      br[1] = lo + 0.5;
      m = 2;
    }
  }
  // findInterval(v, br, rightmost.closed = TRUE, all.inside = TRUE)
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(v[i])) {
      cls[i] = 0;
      continue;
    }
    int b = static_cast<int>(std::upper_bound(br, br + m, v[i]) - br);
    if (v[i] == br[m - 1]) b = m - 1;
    cls[i] = std::min(std::max(b, 1), m - 1);
  }
  *nb = m;
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)

extern "C" int clab_mapbb(const double *x, const double *y, int n, double *bb) CLAB_TRY {
  if (n < 1) return CLAB_EARG;
  return box(x, y, n, bb) ? CLAB_OK : CLAB_EDOM;
}
CLAB_CATCH(CLAB_ENOMEM)

extern "C" int clab_ptcen(const double *x, const double *y, const int *cl, int n, int k, double *cx, double *cy,
                          int *cnt) CLAB_TRY {
  if (n < 1 || k < 0) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (cl[i] < 0 || cl[i] > k) return CLAB_EARG;
  }
  for (int j = 0; j < k; ++j) {
    cx[j] = 0.0;
    cy[j] = 0.0;
    cnt[j] = 0;
  }
  for (int i = 0; i < n; ++i) {
    if (cl[i] == 0) continue;
    const size_t j = z(cl[i] - 1);
    cx[j] += x[i];
    cy[j] += y[i];
    ++cnt[j];
  }
  for (int j = 0; j < k; ++j) {
    cx[j] = cnt[j] > 0 ? cx[j] / cnt[j] : NAN;
    cy[j] = cnt[j] > 0 ? cy[j] / cnt[j] : NAN;
  }
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)

extern "C" int clab_mapsym(const double *x, const double *y, const double *v, int n, double rmax, double *bb,
                           double *r, int *ord, int *nord, double *rm) CLAB_TRY {
  if (n < 1) return CLAB_EARG;
  if (!box(x, y, n, bb)) return CLAB_EDOM;
  double mx = -1.0;
  for (int i = 0; i < n; ++i) {
    if (std::isnan(v[i])) continue;
    if (!std::isfinite(v[i]) || v[i] < 0) return CLAB_EARG;
    mx = std::max(mx, v[i]);
  }
  if (mx < 0) return CLAB_EDOM;
  if (std::isnan(rmax)) rmax = std::max(bb[1] - bb[0], bb[3] - bb[2]) / 15;
  if (!(rmax > 0) || !std::isfinite(rmax)) return CLAB_EARG;
  std::vector<int> o;
  for (int i = 0; i < n; ++i) {
    if (std::isnan(v[i])) {
      r[i] = NAN;
      continue;
    }
    r[i] = mx > 0 ? rmax * std::sqrt(v[i] / mx) : 0.0;
    o.push_back(i);
  }
  // the largest first, so that the smaller are drawn over them; ties in
  // the order given
  std::stable_sort(o.begin(), o.end(), [v](int a, int b) { return v[a] > v[b]; });
  for (size_t q = 0; q < o.size(); ++q) ord[q] = o[q];
  *nord = static_cast<int>(o.size());
  bb[0] -= rmax;
  bb[1] += rmax;
  bb[2] -= rmax;
  bb[3] += rmax;
  *rm = rmax;
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)
