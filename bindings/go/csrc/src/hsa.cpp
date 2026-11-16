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

// The first cells of ord whose areas reach the share sh of the total,
// as weights in hot (1 hot, 0 not); returns their area. The cells whose
// score ties with the last one taken and the first one left (a group the
// cut-off splits) are, by ties: 0 all hot, 1 none hot, 2 each hot by the
// share of the group's area that the cut-off still needs, so that the hot
// area is the target. Any of these makes the result one for every order
// of the rows; a split group needs a tie, so untied scores give the
// first cells whose areas reach the target whatever ties says.
double take(const std::vector<int> &ord, const double *pred, const double *area, double tot, double sh, int ties,
            std::vector<double> &hot) {
  const std::size_t n = ord.size();
  std::fill(hot.begin(), hot.end(), 0.0);
  double a = 0;
  std::size_t k = 0;
  for (; k < n; ++k) {
    if (a >= sh * tot * (1 - 1e-12)) break;
    hot[static_cast<std::size_t>(ord[k])] = 1;
    a += area[ord[k]];
  }
  if (!pred || k == 0 || k == n || pred[ord[k - 1]] != pred[ord[k]]) return a;
  const double s = pred[ord[k]];
  std::size_t b = k;  // the group is ord[b..e)
  while (b > 0 && pred[ord[b - 1]] == s) --b;
  std::size_t e = k;
  while (e < n && pred[ord[e]] == s) ++e;
  double ab = 0, ag = 0;
  for (std::size_t i = 0; i < b; ++i) ab += area[ord[i]];
  for (std::size_t i = b; i < e; ++i) ag += area[ord[i]];
  const double w = ties == 0 ? 1.0 : ties == 1 ? 0.0 : std::fmin(1.0, std::fmax(0.0, (sh * tot - ab) / ag));
  for (std::size_t i = b; i < e; ++i) hot[static_cast<std::size_t>(ord[i])] = w;
  return ab + w * ag;
}

}  // namespace

int clab_hsacc(const double *pred, const int *ishot, const double *obs, const double *area, const double *prev,
               int n, const double *top, int nt, int ties, double *res) CLAB_TRY {
  if (n < 1 || nt < 1 || ties < 0 || ties > 2) return CLAB_EARG;
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
  // the densest cells first; of equal density the larger first, so that
  // the order of the rows does not matter (cells of one density and one
  // area hold the same crimes)
  std::stable_sort(od.begin(), od.end(), [&](int u, int v) {
    const double du = obs[u] / area[u], dv = obs[v] / area[v];
    return du > dv || (du == dv && area[u] > area[v]);
  });
  std::vector<double> hot(static_cast<std::size_t>(n)), best(static_cast<std::size_t>(n));
  // res: nt x 11, column-major: top, cells, area, area share, crimes in
  // the hot cells, all crimes, hit rate, PAI, the most crimes an area as
  // large could hold, PEI, RRI
  for (int q = 0; q < nt; ++q) {
    double a = 0;
    if (ishot) {
      for (int i = 0; i < n; ++i) {
        hot[static_cast<std::size_t>(i)] = ishot[i] != 0 ? 1.0 : 0.0;
        if (ishot[i]) a += area[i];
      }
    } else {
      if (!(top[q] > 0 && top[q] <= 1)) return CLAB_EDOM;
      a = take(ord, pred, area, A, top[q], ties, hot);
    }
    double nh = 0, nc = 0, np = 0;
    for (int i = 0; i < n; ++i) {
      const double h = hot[static_cast<std::size_t>(i)];
      nh += h * obs[i];
      nc += h;
      if (prev) np += h * prev[i];
    }
    if (a > 0) {
      take(od, nullptr, area, A, a / A, 0, best);
    } else {
      std::fill(best.begin(), best.end(), 0.0);
    }
    double nm = 0;
    for (int i = 0; i < n; ++i) nm += best[static_cast<std::size_t>(i)] * obs[i];
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
} CLAB_CATCH(CLAB_ENOMEM)
