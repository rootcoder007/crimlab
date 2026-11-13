// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, distance analysis: the nearest-neighbour index of Clark and
// Evans with Donnelly's edge correction, k-order nearest-neighbour
// distances against their expectations under complete spatial randomness
// (Thompson), and a Monte Carlo test in the study rectangle.
#include <cmath>
#include <limits>

#include "cgeo.h"
#include "crimlab.h"
#include "kutil.h"

namespace {

bool inwin(const double *x, const double *y, int n, const double *w) {
  if (!(std::isfinite(w[0]) && std::isfinite(w[1]) && std::isfinite(w[2]) &&
        std::isfinite(w[3]) && w[0] < w[1] && w[2] < w[3]))
    return false;
  if (std::fabs(w[0]) > 1e150 || std::fabs(w[1]) > 1e150 ||
      std::fabs(w[2]) > 1e150 || std::fabs(w[3]) > 1e150)
    return false;
  for (int i = 0; i < n; ++i)
    if (!(x[i] >= w[0] && x[i] <= w[1] && y[i] >= w[2] && y[i] <= w[3]))
      return false;
  return true;
}

// mean k-th nearest-neighbour distance for k = 1..kmax into md[0..kmax);
// optionally every distance into dist (row i, kmax per row)
void meannn(const double *x, const double *y, int n, const double *win,
            int kmax, int *iwork, double *best, double *md, double *dist) {
  cgeo::grid g;
  const double a = (win[1] - win[0]) * (win[3] - win[2]);
  g.build(x, y, n, win[0], win[1], win[2], win[3],
          std::sqrt(2.0 * a / n), iwork);
  // the kmax sums live in md as (sum, compensation) pairs: md has 2 kmax
  for (int k = 0; k < 2 * kmax; ++k) md[k] = 0.0;
  // points in cell order, for memory locality
  for (int t = 0; t < n; ++t) {
    const int i = g.ord[t];
    cgeo::knn(g, x, y, i, kmax, best);
    for (int k = 0; k < kmax; ++k) {
      const double d = std::sqrt(best[k]);
      if (dist) dist[static_cast<long>(i) * kmax + k] = d;
      kutil::ksum s{md[2 * k], md[2 * k + 1]};
      s.add(d);
      md[2 * k] = s.s;
      md[2 * k + 1] = s.c;
    }
  }
  for (int k = 0; k < kmax; ++k) md[k] = (md[2 * k] + md[2 * k + 1]) / n;
}

}  // namespace

extern "C" {

int clab_nnidx(const double *x, const double *y, int n, const double *win,
               int kmax, int nsim, unsigned int seed, double *dist,
               int *iwork, double *work, double *out) {
  if (n < 2 || kmax < 1 || kmax > n - 1 || nsim < 0) return CLAB_EARG;
  if (!inwin(x, y, n, win)) return CLAB_EDOM;
  const double w = win[1] - win[0], h = win[3] - win[2];
  const double area = w * h, perim = 2.0 * (w + h);
  if (!(area > 0.0) || !std::isfinite(area)) return CLAB_EDOM;
  double *best = work, *md = work + kmax, *sim = work + 3 * kmax;
  meannn(x, y, n, win, kmax, iwork, best, md, dist);
  // copy the observed means out of md before md is reused
  const double lam = n / area;
  double ck = 0.5;  // Thompson: E r_k = ck / sqrt(lam), c_{k+1} = c_k (2k+1)/(2k)
  for (int k = 0; k < kmax; ++k) {
    double *o = out + 6 * k;
    o[0] = md[k];
    o[1] = ck / std::sqrt(lam);
    o[2] = o[0] / o[1];
    o[3] = o[4] = kutil::nan();
    // the standard deviation of r_k under CSR: E r_k^2 = k / (pi lam)
    const double v = (k + 1) / (cgeo::pi * lam) - o[1] * o[1];
    o[5] = std::sqrt(v > 0.0 ? v : 0.0);
    ck *= (2.0 * (k + 1) + 1.0) / (2.0 * (k + 1));
  }
  // order 1: Donnelly's expectation and the Clark-Evans standard error
  double *e = out + 6 * kmax;
  const double ed = 0.5 / std::sqrt(lam) + (0.0514 + 0.0412 / std::sqrt(n)) *
                                               perim / n;
  const double se = std::sqrt((4.0 - cgeo::pi) / (4.0 * cgeo::pi * lam)) /
                    std::sqrt(static_cast<double>(n));
  e[0] = ed;
  e[1] = out[0] / ed;
  e[2] = (out[0] - out[1]) / se;
  e[3] = (out[0] - ed) / se;
  if (nsim > 0) {
    // Monte Carlo test: n uniform points in the rectangle, nsim times
    double *sx = sim, *sy = sim + n;
    int *lo = iwork + 2 * n + 2, *hi = lo + kmax;
    for (int k = 0; k < kmax; ++k) lo[k] = hi[k] = 0;
    cgeo::rng r(seed);
    for (int s = 0; s < nsim; ++s) {
      if (kutil::stop(16.0 * n * kmax)) return CLAB_EINTR;
      for (int i = 0; i < n; ++i) {
        sx[i] = win[0] + w * r.unif();
        sy[i] = win[2] + h * r.unif();
      }
      meannn(sx, sy, n, win, kmax, iwork, best, md, nullptr);
      for (int k = 0; k < kmax; ++k) {
        if (md[k] <= out[6 * k]) ++lo[k];
        if (md[k] >= out[6 * k]) ++hi[k];
      }
    }
    for (int k = 0; k < kmax; ++k) {
      out[6 * k + 3] = (1.0 + lo[k]) / (nsim + 1.0);
      out[6 * k + 4] = (1.0 + hi[k]) / (nsim + 1.0);
    }
  }
  return CLAB_OK;
}

}  // extern "C"
