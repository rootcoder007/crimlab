// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, prospective hot-spotting (Bowers, Johnson and Pease 2004):
// the risk of each grid cell from recent nearby incidents, weighted down
// with distance and with time, its hot cells, and how many later
// incidents they caught.
#include <algorithm>
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

namespace {

// floor(v), clamped to [lo, hi] before it becomes an int
int ix(double v, int lo, int hi) {
  const double f = std::floor(v);
  return f < lo ? lo : f > hi ? hi : static_cast<int>(f);
}

}  // namespace

int clab_prohs(const double *px, const double *py, const double *t, int n, double at, const double *win,
               double cell, double bw, double tw, double tu, double top, double ahead, int nx, int ny, double *z,
               int *hot, double *ev) CLAB_TRY {
  if (n < 1 || nx < 1 || ny < 1 || static_cast<double>(nx) * ny > 1e8) return CLAB_EARG;
  if (!kutil::pos(cell) || !kutil::pos(bw) || !kutil::pos(tw) || !kutil::pos(tu) || !(top > 0 && top <= 1) ||
      !kutil::nonneg(ahead) || !std::isfinite(at)) {
    return CLAB_EDOM;
  }
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(px[i]) || !std::isfinite(py[i]) || !std::isfinite(t[i])) return CLAB_EDOM;
  }
  const std::size_t nc = static_cast<std::size_t>(nx) * ny;
  for (std::size_t k = 0; k < nc; ++k) z[k] = 0.0;
  // each incident of the window [at - tw, at) adds to the cells whose
  // centre lies within bw: 1 / (1 + d / cell) (1 + (at - t) / tu)
  if (bw / cell > 1e4) return CLAB_EDOM;
  const int rr = static_cast<int>(std::ceil(bw / cell)) + 1;
  for (int i = 0; i < n; ++i) {
    if (!(t[i] < at && t[i] >= at - tw)) continue;
    if (kutil::stop(static_cast<double>(rr) * rr * 4)) return CLAB_EINTR;
    const double wt = 1.0 / (1.0 + (at - t[i]) / tu);
    const int ci = ix((px[i] - win[0]) / cell, -rr - 1, nx + rr);
    const int cj = ix((py[i] - win[2]) / cell, -rr - 1, ny + rr);
    for (int a = std::max(0, ci - rr); a <= std::min(nx - 1, ci + rr); ++a) {
      const double cx = win[0] + (a + 0.5) * cell;
      for (int b = std::max(0, cj - rr); b <= std::min(ny - 1, cj + rr); ++b) {
        const double cy = win[2] + (b + 0.5) * cell;
        const double d = std::sqrt((cx - px[i]) * (cx - px[i]) + (cy - py[i]) * (cy - py[i]));
        if (d <= bw) z[static_cast<std::size_t>(a) + static_cast<std::size_t>(nx) * b] += wt / (1.0 + d / cell);
      }
    }
  }
  // the hot cells: the ceiling(top cells) of highest risk, of those with
  // any, ties by position
  std::vector<int> ord(nc);
  for (std::size_t k = 0; k < nc; ++k) ord[k] = static_cast<int>(k);
  std::stable_sort(ord.begin(), ord.end(), [&](int u, int v) { return z[u] > z[v]; });
  const std::size_t kh = static_cast<std::size_t>(std::ceil(top * static_cast<double>(nc) - 1e-9));
  for (std::size_t k = 0; k < nc; ++k) hot[k] = 0;
  std::size_t nh = 0;
  for (std::size_t q = 0; q < kh && q < nc; ++q) {
    if (z[ord[q]] > 0) {
      hot[ord[q]] = 1;
      ++nh;
    }
  }
  // ev: the incidents in [at, at + ahead) inside the window, those in hot
  // cells, the hit rate, the hot share of the area, and the predictive
  // accuracy index
  double ne = 0, nhit = 0;
  for (int i = 0; i < n; ++i) {
    if (!(t[i] >= at && t[i] < at + ahead)) continue;
    const int ci = ix((px[i] - win[0]) / cell, -1, nx);
    const int cj = ix((py[i] - win[2]) / cell, -1, ny);
    if (ci < 0 || ci >= nx || cj < 0 || cj >= ny) continue;
    ne += 1;
    nhit += hot[static_cast<std::size_t>(ci) + static_cast<std::size_t>(nx) * cj];
  }
  const double ash = static_cast<double>(nh) / static_cast<double>(nc);
  ev[0] = ne;
  ev[1] = nhit;
  ev[2] = ne > 0 ? nhit / ne : kutil::nan();
  ev[3] = ash;
  ev[4] = ne > 0 && ash > 0 ? ev[2] / ash : kutil::nan();
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)
