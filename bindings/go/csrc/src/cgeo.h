// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Planar geometry shared by the spatial kernel families (centrography,
// distances, hot spots, journey to crime). Header only, nothing here
// allocates; every function reads and writes the caller's arrays.
#ifndef CLAB_CGEO_H
#define CLAB_CGEO_H

#include <cmath>
#include <cstdint>
#include <limits>

#include "kutil.h"

namespace cgeo {

constexpr double pi = 3.14159265358979323846;

// ---------------------------------------------------------------- sorting

// heapsort and the generator live in kutil.h, shared with the other
// kernel families
using kutil::hsort;
using kutil::rng;

// ---------------------------------------------------- exact orientation

// The sign of the exact value of ax*by - ax*cy - cx*by - ay*bx + ay*cx +
// cy*bx, the orientation of (a, b, c): +1 counter-clockwise, -1 clockwise,
// 0 collinear. A floating-point filter (Shewchuk 1997,
// doi:10.1007/PL00009321) answers the easy cases; the rest are decided on
// the exact sum of the six products, kept as a non-overlapping expansion
// (Shewchuk's grow-expansion), whose sign is that of its largest nonzero
// component. Exact for coordinates that are 0 or of magnitude in
// [1e-140, 1e150] (no product overflows, no lost low part underflows).
int orient(double ax, double ay, double bx, double by, double cx,
                  double cy);
// the coordinate range in which orient is exact
inline bool okcoord(double v) {
  const double a = std::fabs(v);
  return v == 0.0 || (a >= 1e-140 && a <= 1e150);
}

// ------------------------------------------------------------- moments

// Weighted centre and centred second moments of n points; w may be NULL
// (unit weights). Two passes: the mean, then sums of products of centred
// coordinates, each with a compensated sum, so a common offset of 1e7 in
// the coordinates costs nothing. Returns false when the weights are not
// finite and non-negative or sum to zero. m: W, mx, my, sxx, syy, sxy
// (sums, not divided).
bool moments(const double *x, const double *y, const double *w, int n,
                    double *m);

// Eigen-decomposition of the symmetric 2 x 2 matrix [[a, c], [c, b]]:
// l1 >= l2 and the angle t in [0, pi) of the eigenvector of l1,
// anticlockwise from the x axis. The smaller root is taken as det / l1 so
// that it does not cancel. A multiple of the identity has t = 0.
void eig2(double a, double b, double c, double *l1, double *l2,
                 double *t);

// ------------------------------------------------------- Weber point

// The point minimising sum w_i |p_i - c| (the centre of minimum distance,
// the spatial L1 median): Weiszfeld's iteration (Weiszfeld 1937,
// translated doi:10.1007/s10479-008-0352-z) with the modification of Vardi
// and Zhang (2000, doi:10.1073/pnas.97.4.1423) that keeps it well defined
// and monotone when an iterate lands on a data point. Coordinates are
// taken relative to (ox, oy); (cx, cy) holds the start on entry and the
// result on exit, both relative to (ox, oy).
//   With S = sum over p_i != c of w_i / d_i, R = sum of w_i (p_i - c) / d_i
//   and eta = the weight sitting exactly at c, the step is
//   c <- c + (1 - g) R / S, g = min(1, eta / |R|) (g = 0 when eta = 0).
// c is optimal if and only if |R| <= eta (the subgradient condition); the
// iteration stops there (code 2), when a step is at most tol * scale
// (code 1) or after maxit steps (code 0). out: objective at the result,
// max(0, |R| - eta) at the result, the number of steps, the code.
void weber(const double *x, const double *y, const double *w, int n,
                  double ox, double oy, double scale, double tol, int maxit,
                  double *cx, double *cy, double *out);


// ------------------------------------------------------------ point grid

// A uniform grid over the rectangle [x0, x1] x [y0, y1] holding n points
// by counting sort: the points of cell c are ord[start[c] .. start[c + 1]).
// The rectangle is split into ceil(width / c) by ceil(height / c) cells,
// c being `cell` grown by factors of 1.25 until there are at most n + 1
// cells, so a cell's side is at most c and may be shorter than `cell`: a
// search to distance r must reach from cx(x - r) to cx(x + r), not one
// cell either way. iwork needs n + (n + 2) ints. Points must lie in the
// rectangle.
struct grid {
  double x0, y0, wx, wy;
  int gx, gy;
  int *start, *ord;
  int cx(double x) const {
    int c = static_cast<int>((x - x0) / wx);
    return c < 0 ? 0 : (c >= gx ? gx - 1 : c);
  }
  int cy(double y) const {
    int c = static_cast<int>((y - y0) / wy);
    return c < 0 ? 0 : (c >= gy ? gy - 1 : c);
  }
  void build(const double *x, const double *y, int n, double ax0, double ax1,
             double ay0, double ay1, double cell, int *iwork);
};

// The k smallest squared distances from point i to the other points
// (coincident points count, at distance 0), ascending in best[0..k), by
// rings of cells around i's cell until no unvisited cell can be nearer.
// Requires k <= n - 1.
void knn(const grid &g, const double *x, const double *y, int i, int k,
                double *best);

}  // namespace cgeo

#endif
