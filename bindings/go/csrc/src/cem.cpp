// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Coarsened exact matching (crimlab.h, cemtch()): Iacus, King and Porro
// (2012). Rows with the same coarsened covariates form a stratum; strata
// holding both treated and control rows are kept, with weight 1 for the
// treated and (m_T^s / m_C^s) (m_C / m_T) for the controls (the weights
// of the estimate on the treated), and the multivariate imbalance L1 of
// Iacus, King and Porro (2011) is measured before and after on a second
// coarsening. Work: O(n k log n).
#include <cmath>
#include <cstddef>
#include <map>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

namespace {

using std::size_t;
inline size_t z(int i) { return static_cast<size_t>(i); }

// The cells of the n x k codes (column-major) numbered from 0 in order
// of first appearance.
int cells(const int *code, int n, int k, std::vector<int> &id) {
  std::map<std::vector<int>, int> h;
  std::vector<int> key(z(k));
  id.assign(z(n), 0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < k; ++j) key[z(j)] = code[z(j) * z(n) + z(i)];
    id[z(i)] = h.emplace(key, static_cast<int>(h.size())).first->second;
  }
  return static_cast<int>(h.size());
}

// L1 = sum over cells |f - g| / 2, f and g the treated and control
// shares of the weight w in each cell; NaN when a group has no weight.
double l1(const std::vector<int> &id, int nc, const int *tr, const double *w, int n) {
  std::vector<double> f(z(nc), 0.0), g(z(nc), 0.0);
  double sf = 0.0, sg = 0.0;
  for (int i = 0; i < n; ++i) {
    const double v = w ? w[i] : 1.0;
    if (tr[i]) {
      f[z(id[z(i)])] += v;
      sf += v;
    } else {
      g[z(id[z(i)])] += v;
      sg += v;
    }
  }
  if (!(sf > 0) || !(sg > 0)) return NAN;
  double s = 0.0;
  for (int c = 0; c < nc; ++c) s += std::fabs(f[z(c)] / sf - g[z(c)] / sg);
  return s / 2;
}

}  // namespace

extern "C" int clab_cemtch(const int *code, int n, int k, const int *tr, const int *lcode, int lk, int *strat,
                           double *w, double *st) CLAB_TRY {
  if (n < 1 || k < 0 || lk < 0) return CLAB_EARG;
  if (static_cast<double>(n) * (k > lk ? k : lk) > kutil::maxcell) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (tr[i] != 0 && tr[i] != 1) return CLAB_EARG;
  }
  std::vector<int> id;
  const int ns = cells(code, n, k, id);
  if (kutil::stop(static_cast<double>(n) * (k + 1))) return CLAB_EINTR;
  std::vector<int> nt(z(ns), 0), nc(z(ns), 0);
  for (int i = 0; i < n; ++i) (tr[i] ? nt : nc)[z(id[z(i)])] += 1;
  // the kept strata renumbered from 1 in order of first appearance
  std::vector<int> keep(z(ns), 0);
  int nk = 0;
  double mt = 0.0, mc = 0.0;
  for (int i = 0; i < n; ++i) {
    const int s = id[z(i)];
    if (nt[z(s)] > 0 && nc[z(s)] > 0 && keep[z(s)] == 0) keep[z(s)] = ++nk;
  }
  for (int i = 0; i < n; ++i) {
    const int s = id[z(i)];
    strat[i] = keep[z(s)] > 0 ? keep[z(s)] : 0;
    if (strat[i] > 0) (tr[i] ? mt : mc) += 1.0;
  }
  for (int i = 0; i < n; ++i) {
    const int s = id[z(i)];
    if (strat[i] == 0) {
      w[i] = 0.0;
    } else {
      w[i] = tr[i] ? 1.0 : (static_cast<double>(nt[z(s)]) / nc[z(s)]) * (mc / mt);
    }
  }
  std::vector<int> lid;
  const int nl = cells(lcode, n, lk, lid);
  st[0] = l1(lid, nl, tr, nullptr, n);
  st[1] = l1(lid, nl, tr, w, n);
  st[2] = nk;
  st[3] = mt;
  st[4] = mc;
  st[5] = ns;
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)
