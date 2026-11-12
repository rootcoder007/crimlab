// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, distance analysis: Ripley's K function in a rectangle with
// the uncorrected, border, modified border and isotropic estimators,
// Besag's L, and simulation envelopes with a global Monte Carlo test
// under complete spatial randomness.
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

bool radii(const double *r, int m) {
  for (int k = 0; k < m; ++k) {
    if (!std::isfinite(r[k]) || r[k] < 0.0) return false;
    if (k > 0 && !(r[k] > r[k - 1])) return false;
  }
  return true;
}

// The share of the circle of radius d about (px, py) that lies outside the
// rectangle, as an angle: each edge at distance e < d removes an arc of
// 2 acos(e / d); arcs of two adjacent edges overlap by acos(e1 / d) +
// acos(e2 / d) - pi / 2 when the corner lies inside the circle. Opposite
// edges' arcs never overlap (each is under pi wide), so this is exact.
double outside(double px, double py, const double *w, double d) {
  const double e[4] = {px - w[0], w[1] - px, py - w[2], w[3] - py};
  double a[4];
  double t = 0.0;
  for (int s = 0; s < 4; ++s) {
    a[s] = e[s] < d ? std::acos(e[s] / d) : 0.0;
    t += 2.0 * a[s];
  }
  const int cn[4][2] = {{0, 2}, {0, 3}, {1, 2}, {1, 3}};
  for (int c = 0; c < 4; ++c) {
    const double u = e[cn[c][0]], v = e[cn[c][1]];
    if (u * u + v * v < d * d) t -= a[cn[c][0]] + a[cn[c][1]] - 0.5 * cgeo::pi;
  }
  return t;
}

// first index k with r[k] >= d, or m
int bin(const double *r, int m, double d) {
  int lo = 0, hi = m;
  while (lo < hi) {
    const int mid = lo + (hi - lo) / 2;
    if (r[mid] >= d)
      hi = mid;
    else
      lo = mid + 1;
  }
  return lo;
}

// K at the m radii into k; work: 4 m doubles
void kfun(const double *x, const double *y, int n, const double *w,
          const double *r, int m, int corr, int *iwork, double *work,
          double *k) {
  const double A = (w[1] - w[0]) * (w[3] - w[2]);
  const double rmax = r[m - 1];
  double *acc = work, *loc = work + 2 * m, *den = work + 3 * m;
  for (int j = 0; j < 4 * m; ++j) acc[j] = 0.0;
  const double cell = rmax > 0.0 ? rmax : std::sqrt(A / n);
  cgeo::grid g;
  g.build(x, y, n, w[0], w[1], w[2], w[3], cell, iwork);
  // how many cells a search must reach on each side
  const int sx = static_cast<int>(std::ceil(rmax / g.wx));
  const int sy = static_cast<int>(std::ceil(rmax / g.wy));
  const double r2 = rmax * rmax;
  // centres in cell order, for memory locality
  for (int tt = 0; tt < n; ++tt) {
    const int i = g.ord[tt];
    const double px = x[i], py = y[i];
    const int ci = g.cx(px), cj = g.cy(py);
    const double b = std::fmin(std::fmin(px - w[0], w[1] - px),
                               std::fmin(py - w[2], w[3] - py));
    if (corr == 1 || corr == 2)
      for (int j = 0; j < m; ++j) loc[j] = 0.0;
    for (int bb = cj - sy; bb <= cj + sy; ++bb) {
      if (bb < 0 || bb >= g.gy) continue;
      for (int aa = ci - sx; aa <= ci + sx; ++aa) {
        if (aa < 0 || aa >= g.gx) continue;
        const int c = bb * g.gx + aa;
        for (int t = g.start[c]; t < g.start[c + 1]; ++t) {
          const int j = g.ord[t];
          if (j == i) continue;
          const double dx = x[j] - px, dy = y[j] - py;
          const double dd = dx * dx + dy * dy;
          if (dd > r2) continue;
          const double d = std::sqrt(dd);
          const int q = bin(r, m, d);
          if (q == m) continue;
          double wt = 1.0;
          if (corr == 3 && d > 0.0)
            wt = 2.0 * cgeo::pi / (2.0 * cgeo::pi - outside(px, py, w, d));
          if (corr == 1 || corr == 2) {
            loc[q] += 1.0;
          } else {
            kutil::ksum s{acc[2 * q], acc[2 * q + 1]};
            s.add(wt);
            acc[2 * q] = s.s;
            acc[2 * q + 1] = s.c;
          }
        }
      }
    }
    if (corr == 1 || corr == 2) {
      // this point's cumulative counts enter every radius it can see whole
      double run = 0.0;
      for (int j = 0; j < m && r[j] <= b; ++j) {
        run += loc[j];
        acc[2 * j] += run;  // integer counts: exact
        den[j] += 1.0;
      }
    }
  }
  // cumulate over radii (compensated) and scale
  kutil::ksum run;
  const double nn = static_cast<double>(n) * (n - 1);
  for (int j = 0; j < m; ++j) {
    double v;
    if (corr == 1 || corr == 2) {
      v = acc[2 * j];
    } else {
      run.add(acc[2 * j]);
      run.add(acc[2 * j + 1]);
      v = run.val();
    }
    switch (corr) {
      case 1:  // reduced sample, intensity n / A
        k[j] = den[j] > 0.0 ? v / den[j] * A / n : kutil::nan();
        break;
      case 2: {  // modified border: |W eroded by r|, NaN once it is empty
        const double ex = w[1] - w[0] - 2.0 * r[j], ey = w[3] - w[2] - 2.0 * r[j];
        const bool ok = ex > 0.0 && ey > 0.0;
        k[j] = ok ? v * A * A / (nn * (ex * ey)) : kutil::nan();
        break;
      }
      default:  // none, isotropic
        k[j] = v * A / nn;
    }
  }
}

}  // namespace

extern "C" {

int clab_ripk(const double *x, const double *y, int n, const double *win,
              const double *r, int m, int corr, int *iwork, double *work,
              double *k) {
  if (n < 2 || m < 1 || corr < 0 || corr > 3) return CLAB_EARG;
  if (!inwin(x, y, n, win) || !radii(r, m)) return CLAB_EDOM;
  kfun(x, y, n, win, r, m, corr, iwork, work, k);
  return CLAB_OK;
}

int clab_ripenv(const double *x, const double *y, int n, const double *win,
                const double *r, int m, int corr, int nsim, unsigned int seed,
                int *iwork, double *work, double *out) {
  if (n < 2 || m < 1 || corr < 0 || corr > 3 || nsim < 1) return CLAB_EARG;
  if (!inwin(x, y, n, win) || !radii(r, m)) return CLAB_EDOM;
  double *kw = work, *ks = work + 4 * m, *sx = work + 5 * m, *sy = sx + n;
  double *obs = out, *lo = out + m, *hi = out + 2 * m, *mean = out + 3 * m;
  // L(r) - r for the data
  kfun(x, y, n, win, r, m, corr, iwork, kw, ks);
  double tobs = 0.0;
  for (int j = 0; j < m; ++j) {
    obs[j] = std::sqrt(ks[j] / cgeo::pi) - r[j];
    lo[j] = std::numeric_limits<double>::infinity();
    hi[j] = -lo[j];
    mean[j] = 0.0;
    if (std::isfinite(obs[j]) && std::fabs(obs[j]) > tobs) tobs = std::fabs(obs[j]);
  }
  const double w = win[1] - win[0], h = win[3] - win[2];
  cgeo::rng rg(seed);
  int ge = 0;
  for (int s = 0; s < nsim; ++s) {
    for (int i = 0; i < n; ++i) {
      sx[i] = win[0] + w * rg.unif();
      sy[i] = win[2] + h * rg.unif();
    }
    kfun(sx, sy, n, win, r, m, corr, iwork, kw, ks);
    double t = 0.0;
    for (int j = 0; j < m; ++j) {
      const double l = std::sqrt(ks[j] / cgeo::pi) - r[j];
      if (l < lo[j]) lo[j] = l;
      if (l > hi[j]) hi[j] = l;
      mean[j] += l / nsim;
      if (std::isfinite(l) && std::fabs(l) > t) t = std::fabs(l);
    }
    if (t >= tobs) ++ge;
  }
  out[4 * m] = tobs;
  out[4 * m + 1] = (1.0 + ge) / (nsim + 1.0);
  return CLAB_OK;
}

}  // extern "C"
