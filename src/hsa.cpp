// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, the accuracy of a hot-spot map: the hit rate and predictive
// accuracy index of Chainey, Tompson and Uhlig (2008), the predictive
// efficiency index (Drawve 2016) and the recapture rate index of Levine
// (2008).
#include <algorithm>
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

namespace {

// The first cells of ord whose areas reach the share sh of the total;
// marks them in hot and returns their area.
double take(const std::vector<int> &ord, const double *area, double tot, double sh, std::vector<char> &hot) {
  double a = 0;
  for (char &h : hot) h = 0;
  for (int k : ord) {
    if (a >= sh * tot * (1 - 1e-12)) break;
    hot[static_cast<std::size_t>(k)] = 1;
    a += area[k];
  }
  return a;
}

}  // namespace

int clab_hsacc(const double *pred, const int *ishot, const double *obs, const double *area, const double *prev,
               int n, const double *top, int nt, double *res) {
  if (n < 1 || nt < 1) return CLAB_EARG;
  double A = 0, N = 0, Np = 0;
  for (int i = 0; i < n; ++i) {
    if (!kutil::nonneg(obs[i]) || !kutil::pos(area[i]) || (prev && !kutil::nonneg(prev[i]))) return CLAB_EDOM;
    if (!ishot && !std::isfinite(pred[i])) return CLAB_EDOM;
    A += area[i];
    N += obs[i];
    if (prev) Np += prev[i];
  }
  if (!(N > 0)) return CLAB_EDOM;
  std::vector<int> ord(static_cast<std::size_t>(n)), od(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) ord[static_cast<std::size_t>(i)] = od[static_cast<std::size_t>(i)] = i;
  if (!ishot) std::stable_sort(ord.begin(), ord.end(), [&](int u, int v) { return pred[u] > pred[v]; });
  std::stable_sort(od.begin(), od.end(), [&](int u, int v) { return obs[u] / area[u] > obs[v] / area[v]; });
  std::vector<char> hot(static_cast<std::size_t>(n)), best(static_cast<std::size_t>(n));
  // res: nt x 11, column-major: top, cells, area, area share, crimes in
  // the hot cells, all crimes, hit rate, PAI, the most crimes an area as
  // large could hold, PEI, RRI
  for (int q = 0; q < nt; ++q) {
    double a = 0;
    if (ishot) {
      for (int i = 0; i < n; ++i) {
        hot[static_cast<std::size_t>(i)] = ishot[i] != 0;
        if (ishot[i]) a += area[i];
      }
    } else {
      if (!(top[q] > 0 && top[q] <= 1)) return CLAB_EDOM;
      a = take(ord, area, A, top[q], hot);
    }
    double nh = 0, nc = 0, np = 0;
    for (int i = 0; i < n; ++i) {
      if (!hot[static_cast<std::size_t>(i)]) continue;
      nh += obs[i];
      nc += 1;
      if (prev) np += prev[i];
    }
    if (a > 0) {
      take(od, area, A, a / A, best);
    } else {
      std::fill(best.begin(), best.end(), 0);
    }
    double nm = 0;
    for (int i = 0; i < n; ++i) nm += best[static_cast<std::size_t>(i)] ? obs[i] : 0.0;
    const double hit = nh / N, sh = a / A;
    res[q] = ishot ? sh : top[q];
    res[nt + q] = nc;
    res[2 * nt + q] = a;
    res[3 * nt + q] = sh;
    res[4 * nt + q] = nh;
    res[5 * nt + q] = N;
    res[6 * nt + q] = hit;
    res[7 * nt + q] = sh > 0 ? hit / sh : kutil::nan();
    res[8 * nt + q] = nm;
    res[9 * nt + q] = nm > 0 ? nh / nm : kutil::nan();
    res[10 * nt + q] = prev && np > 0 && Np > 0 ? (nh / np) / (N / Np) : kutil::nan();
  }
  return CLAB_OK;
}
