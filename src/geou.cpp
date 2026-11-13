// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, geographic utilities: point in polygon and the areal
// crosswalk of counts between two sets of units.
#include <algorithm>
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

namespace {

// (px, py) on the segment from a to b (closed)
bool onseg(double px, double py, double ax, double ay, double bx, double by) {
  const double cr = (bx - ax) * (py - ay) - (by - ay) * (px - ax);
  const double len = std::fabs(bx - ax) + std::fabs(by - ay);
  if (std::fabs(cr) > 1e-12 * len * (len + std::fabs(px - ax) + std::fabs(py - ay))) return false;
  return px >= std::min(ax, bx) && px <= std::max(ax, bx) && py >= std::min(ay, by) &&
         py <= std::max(ay, by);
}

}  // namespace

int clab_pinpl(const double *px, const double *py, int n, const double *vx, const double *vy,
               const int *pid, const int *rid, int nv, int np, int *reg, int *nreg) {
  if (n < 1 || nv < 3 || np < 1) return CLAB_EARG;
  for (int v = 0; v < nv; ++v) {
    if (pid[v] < 0 || pid[v] >= np ||
        (v > 0 && (pid[v] < pid[v - 1] || (pid[v] == pid[v - 1] && rid[v] < rid[v - 1])))) {
      return CLAB_EARG;
    }
    if (!std::isfinite(vx[v]) || !std::isfinite(vy[v])) return CLAB_EDOM;
  }
  // each polygon's vertices [s, e) and box
  std::vector<int> st(np + 1, 0);
  std::vector<double> box(4 * static_cast<std::size_t>(np));
  for (int v = 0; v < nv; ++v) ++st[pid[v] + 1];
  for (int q = 0; q < np; ++q) st[q + 1] += st[q];
  for (int q = 0; q < np; ++q) {
    double *b = box.data() + 4 * q;
    b[0] = b[2] = HUGE_VAL;
    b[1] = b[3] = -HUGE_VAL;
    for (int v = st[q]; v < st[q + 1]; ++v) {
      b[0] = std::min(b[0], vx[v]);
      b[1] = std::max(b[1], vx[v]);
      b[2] = std::min(b[2], vy[v]);
      b[3] = std::max(b[3], vy[v]);
    }
  }
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(px[i]) || !std::isfinite(py[i])) return CLAB_EDOM;
    reg[i] = -1;
    nreg[i] = 0;
    for (int q = 0; q < np; ++q) {
      const double *b = box.data() + 4 * q;
      if (px[i] < b[0] || px[i] > b[1] || py[i] < b[2] || py[i] > b[3]) continue;
      // even-odd over the polygon's rings, each closed from its last
      // vertex to its first; a point on an edge is inside
      bool in = false, edge = false;
      int r0 = st[q];
      for (int v = st[q]; v < st[q + 1] && !edge; ++v) {
        const bool last = v + 1 == st[q + 1] || rid[v + 1] != rid[v];
        const int w = last ? r0 : v + 1;
        if (onseg(px[i], py[i], vx[v], vy[v], vx[w], vy[w])) edge = true;
        if ((vy[v] > py[i]) != (vy[w] > py[i]) &&
            px[i] < vx[v] + (py[i] - vy[v]) * (vx[w] - vx[v]) / (vy[w] - vy[v])) {
          in = !in;
        }
        if (last) r0 = v + 1;
      }
      if (in || edge) {
        if (reg[i] < 0) reg[i] = q;
        ++nreg[i];
      }
    }
  }
  return CLAB_OK;
}

int clab_arwt(const double *val, int ns, const int *from, const int *to, const double *area, int np,
              int nt, double *out, double *share, int *main, double *lost) {
  if (ns < 1 || np < 1 || nt < 1) return CLAB_EARG;
  std::vector<double> sa(ns, 0.0);
  for (int q = 0; q < np; ++q) {
    if (from[q] < 0 || from[q] >= ns || to[q] < 0 || to[q] >= nt) return CLAB_EARG;
    if (!kutil::nonneg(area[q])) return CLAB_EDOM;
    sa[from[q]] += area[q];
  }
  for (int s = 0; s < ns; ++s) {
    if (!std::isfinite(val[s])) return CLAB_EDOM;
  }
  std::fill(out, out + nt, 0.0);
  std::vector<double> best(ns, -1.0);
  for (int s = 0; s < ns; ++s) main[s] = -1;
  for (int q = 0; q < np; ++q) {
    const int s = from[q];
    // a source whose pairs have no area splits evenly among them
    share[q] = sa[s] > 0 ? area[q] / sa[s] : kutil::nan();
    if (share[q] > best[s]) {
      best[s] = share[q];
      main[s] = to[q];
    }
  }
  std::vector<int> cnt(ns, 0);
  for (int q = 0; q < np; ++q) ++cnt[from[q]];
  for (int q = 0; q < np; ++q) {
    if (std::isnan(share[q])) {
      share[q] = 1.0 / cnt[from[q]];
      if (main[from[q]] < 0) main[from[q]] = to[q];
    }
    out[to[q]] += val[from[q]] * share[q];
  }
  double l = 0;
  for (int s = 0; s < ns; ++s) {
    if (cnt[s] == 0) l += val[s];
  }
  *lost = l;
  return CLAB_OK;
}
