// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, geographic profiling and risk terrain: the circle
// hypothesis (clab_crchyp), the hit score of a profile (clab_hitscr) and
// risk terrain layers (clab_rtmlay).
#include <algorithm>
#include <cmath>

#include "crimlab.h"
#include "kutil.h"

namespace {

bool fin(const double *x, int n) {
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(x[i])) return false;
  }
  return true;
}

}  // namespace

int clab_crchyp(const double *x, const double *y, int n, const double *hx, const double *hy,
                int nh, double *st, int *pair, int *mar) {
  if (n < 2 || nh < 0) return CLAB_EARG;
  if (!fin(x, n) || !fin(y, n) || !fin(hx, nh) || !fin(hy, nh)) return CLAB_EDOM;
  // ponytail: O(n^2) pairs; a series has tens of sites, not thousands
  double best = -1.0;
  int bi = 0, bj = 1;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const double d = std::hypot(x[i] - x[j], y[i] - y[j]);
      if (d > best) {
        best = d;
        bi = i;
        bj = j;
      }
    }
  }
  const double cx = (x[bi] + x[bj]) / 2, cy = (y[bi] + y[bj]) / 2, r = best / 2;
  const double lim = r + 1e-12 * std::max(1.0, r);
  int in = 0;
  for (int i = 0; i < n; ++i) in += std::hypot(x[i] - cx, y[i] - cy) <= lim;
  for (int k = 0; k < nh; ++k) mar[k] = std::hypot(hx[k] - cx, hy[k] - cy) <= lim;
  st[0] = cx;
  st[1] = cy;
  st[2] = r;
  st[3] = static_cast<double>(in) / n;
  pair[0] = bi;
  pair[1] = bj;
  return CLAB_OK;
}

int clab_hitscr(const double *s, int n, const int *home, int nh, int *above, int *ties,
                double *hit, double *hitt) {
  if (n < 1 || nh < 1) return CLAB_EARG;
  if (!fin(s, n)) return CLAB_EDOM;
  for (int k = 0; k < nh; ++k) {
    if (home[k] < 0 || home[k] >= n) return CLAB_EARG;
    const double h = s[home[k]];
    int a = 0, t = 0;
    for (int i = 0; i < n; ++i) {
      a += s[i] > h;
      t += s[i] == h;
    }
    above[k] = a;
    ties[k] = t;
    hit[k] = 100.0 * a / n;
    hitt[k] = 100.0 * (a + t) / n;
  }
  return CLAB_OK;
}

int clab_rtmlay(const double *cx, const double *cy, int nc, const double *fx, const double *fy,
                int nf, double rad, int kind, double *out) {
  if (nc < 1 || nf < 1 || kind < 0 || kind > 2 || (kind == 2 && nc < 2)) return CLAB_EARG;
  if (!kutil::pos(rad) || !fin(cx, nc) || !fin(cy, nc) || !fin(fx, nf) || !fin(fy, nf)) {
    return CLAB_EDOM;
  }
  for (int i = 0; i < nc; ++i) {
    kutil::ksum s;
    for (int f = 0; f < nf; ++f) {
      const double d = std::hypot(cx[i] - fx[f], cy[i] - fy[f]);
      s.add(kind == 2 ? std::exp(-0.5 * (d / rad) * (d / rad)) : d <= rad ? 1.0 : 0.0);
    }
    out[i] = kind == 0 ? (s.val() > 0 ? 1.0 : 0.0) : s.val();
  }
  if (kind == 2) {
    kutil::ksum m;
    for (int i = 0; i < nc; ++i) m.add(out[i]);
    const double mu = m.val() / nc;
    kutil::ksum q;
    for (int i = 0; i < nc; ++i) q.add((out[i] - mu) * (out[i] - mu));
    const double sd = std::sqrt(q.val() / (nc - 1));
    for (int i = 0; i < nc; ++i) out[i] = sd > 0 ? (out[i] - mu) / sd : 0.0;
  }
  return CLAB_OK;
}
