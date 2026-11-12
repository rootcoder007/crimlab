// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, spatial weights: k nearest neighbours and distance bands of
// planar points, found on a uniform grid of cells, and the weights of a
// neighbour structure (binary or inverse distance, raw or row-standardised).
// Neighbour structures are compressed sparse rows: the neighbours of unit i
// are ci[rp[i] .. rp[i + 1]). The statistics of src/lsa.cpp read them.
#include <cmath>
#include <utility>

#include "cgeo.h"
#include "crimlab.h"
#include "kutil.h"

namespace {

double dist(const double *x, const double *y, int i, int j, int metric) {
  const double dx = x[i] - x[j], dy = y[i] - y[j];
  return metric ? std::fabs(dx) + std::fabs(dy) : std::hypot(dx, dy);
}

bool finite2(const double *x, const double *y, int n) {
  for (int i = 0; i < n; ++i)
    if (!std::isfinite(x[i]) || !std::isfinite(y[i])) return false;
  return true;
}

// the bounding box of the points, widened to a positive extent
void bbox(const double *x, const double *y, int n, double *b) {
  b[0] = b[1] = x[0];
  b[2] = b[3] = y[0];
  for (int i = 1; i < n; ++i) {
    b[0] = std::fmin(b[0], x[i]);
    b[1] = std::fmax(b[1], x[i]);
    b[2] = std::fmin(b[2], y[i]);
    b[3] = std::fmax(b[3], y[i]);
  }
  const double e = 1e-9 * (1.0 + std::fabs(b[1] - b[0]) + std::fabs(b[3] - b[2]));
  if (!(b[1] > b[0])) b[1] = b[0] + e;
  if (!(b[3] > b[2])) b[3] = b[2] + e;
}

// (d, j) before (e, k): the order of neighbours, nearest first, then by
// index
bool before(double d, int j, double e, int k) { return d < e || (d == e && j < k); }

// sift the worst candidate of a max-heap of m (distance, index) pairs down
void sift(double *hd, int *hj, int m, int at) {
  for (;;) {
    int w = at;
    const int l = 2 * at + 1, r = l + 1;
    if (l < m && before(hd[w], hj[w], hd[l], hj[l])) w = l;
    if (r < m && before(hd[w], hj[w], hd[r], hj[r])) w = r;
    if (w == at) return;
    std::swap(hd[w], hd[at]);
    std::swap(hj[w], hj[at]);
    at = w;
  }
}

}  // namespace

extern "C" {

int clab_swknb(const double *x, const double *y, int n, int k, int metric,
               int *nb, double *d, int *iwork) {
  if (n < 2 || k < 1 || k > n - 1 || metric < 0 || metric > 1) return CLAB_EARG;
  if (!finite2(x, y, n)) return CLAB_EDOM;
  double b[4];
  bbox(x, y, n, b);
  // cells holding about k + 1 points each on average
  const double cell = std::sqrt((b[1] - b[0]) * (b[3] - b[2]) * (k + 1.0) / n);
  cgeo::grid g;
  g.build(x, y, n, b[0], b[1], b[2], b[3], cell, iwork);
  const double wmin = std::fmin(g.wx, g.wy);
  for (int i = 0; i < n; ++i) {
    double *hd = d + static_cast<long>(i) * k;
    int *hj = nb + static_cast<long>(i) * k;
    int m = 0;
    const int ci = g.cx(x[i]), cj = g.cy(y[i]);
    // rings of cells around i's cell; a point of ring r + 1 is at least
    // r wmin away (Euclidean, and so also Manhattan)
    for (int r = 0;; ++r) {
      if (m == k && r > 0 && (r - 1) * wmin > hd[0]) break;
      if (r > g.gx && r > g.gy) break;
      for (int a = ci - r; a <= ci + r; ++a) {
        if (a < 0 || a >= g.gx) continue;
        for (int c = cj - r; c <= cj + r; ++c) {
          if (c < 0 || c >= g.gy) continue;
          if (a != ci - r && a != ci + r && c != cj - r && c != cj + r) continue;
          const int cell2 = c * g.gx + a;
          for (int q = g.start[cell2]; q < g.start[cell2 + 1]; ++q) {
            const int j = g.ord[q];
            if (j == i) continue;
            const double e = dist(x, y, i, j, metric);
            if (m < k) {
              // insert into the heap
              int at = m++;
              hd[at] = e;
              hj[at] = j;
              while (at > 0) {
                const int p = (at - 1) / 2;
                if (!before(hd[p], hj[p], hd[at], hj[at])) break;
                std::swap(hd[p], hd[at]);
                std::swap(hj[p], hj[at]);
                at = p;
              }
            } else if (before(e, j, hd[0], hj[0])) {
              hd[0] = e;
              hj[0] = j;
              sift(hd, hj, m, 0);
            }
          }
        }
      }
    }
    // heap to ascending order
    for (int e = m - 1; e > 0; --e) {
      std::swap(hd[0], hd[e]);
      std::swap(hj[0], hj[e]);
      sift(hd, hj, e, 0);
    }
  }
  return CLAB_OK;
}

int clab_swband(const double *x, const double *y, int n, double dmax,
                int metric, int cap, int *rp, int *ci, double *d, int *nnz,
                int *iwork) {
  if (n < 1 || cap < 0 || metric < 0 || metric > 1) return CLAB_EARG;
  if (!kutil::pos(dmax) || !finite2(x, y, n)) return CLAB_EDOM;
  double b[4];
  bbox(x, y, n, b);
  cgeo::grid g;
  g.build(x, y, n, b[0], b[1], b[2], b[3], dmax, iwork);
  long t = 0;
  rp[0] = 0;
  for (int i = 0; i < n; ++i) {
    const long row = t;
    // the cells the square of half-side dmax about the point touches (the
    // grid's cells may be narrower than dmax)
    const int a0 = g.cx(x[i] - dmax), a1 = g.cx(x[i] + dmax);
    const int c0 = g.cy(y[i] - dmax), c1 = g.cy(y[i] + dmax);
    for (int c = c0; c <= c1; ++c) {
      for (int a = a0; a <= a1; ++a) {
        const int cell = c * g.gx + a;
        for (int q = g.start[cell]; q < g.start[cell + 1]; ++q) {
          const int j = g.ord[q];
          if (j == i) continue;
          const double e = dist(x, y, i, j, metric);
          if (e <= dmax) {
            if (t < cap) {
              ci[t] = j;
              d[t] = e;
            }
            ++t;
          }
        }
      }
    }
    if (t > 2147483647L) return CLAB_EARG;
    // each row in increasing order of neighbour (insertion sort of the
    // stored part; rows are short)
    const long end = t < cap ? t : cap;
    for (long u = row + 1; u < end; ++u) {
      const int j = ci[u];
      const double e = d[u];
      long v = u - 1;
      while (v >= row && ci[v] > j) {
        ci[v + 1] = ci[v];
        d[v + 1] = d[v];
        --v;
      }
      ci[v + 1] = j;
      d[v + 1] = e;
    }
    rp[i + 1] = static_cast<int>(t);
  }
  *nnz = static_cast<int>(t);
  return CLAB_OK;
}

int clab_swwt(const int *rp, const double *d, int n, double pw, int rowstd,
              double *w) {
  if (n < 1 || rp[0] != 0 || rowstd < 0 || rowstd > 1) return CLAB_EARG;
  for (int i = 0; i < n; ++i)
    if (rp[i + 1] < rp[i]) return CLAB_EARG;
  if (!(pw >= 0.0) || !std::isfinite(pw)) return CLAB_EDOM;
  for (int i = 0; i < n; ++i) {
    double s = 0.0;
    for (int q = rp[i]; q < rp[i + 1]; ++q) {
      if (pw == 0.0) {
        w[q] = 1.0;
      } else {
        if (d == nullptr || !(d[q] > 0.0) || !std::isfinite(d[q])) return CLAB_EDOM;
        w[q] = std::pow(d[q], -pw);
      }
      s += w[q];
    }
    if (rowstd && s > 0.0)
      for (int q = rp[i]; q < rp[i + 1]; ++q) w[q] /= s;
  }
  return CLAB_OK;
}

}  // extern "C"
