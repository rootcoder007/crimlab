// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, near repeats: the Knox table of pairs of incidents by band of
// distance and band of time, its expectation under exchangeable times and
// a permutation test of each cell.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

namespace {

// increasing, positive and finite
bool edges(const double *e, int k) {
  for (int q = 0; q < k; ++q) {
    if (!kutil::pos(e[q])) return false;
    if (q > 0 && !(e[q] > e[q - 1])) return false;
  }
  return true;
}

// the band (lo, hi] holding v (0 in the first), or -1 above the last edge
int band(double v, const double *e, int k) {
  const int q = static_cast<int>(std::lower_bound(e, e + k, v) - e);
  return q < k ? q : -1;
}

// the observed count of each cell for the times t (cells: distance band
// fastest)
void cells(const std::vector<int> &pi, const std::vector<int> &pj, const std::vector<int> &sb,
           const double *t, const double *dt, int ks, int kt, double *obs) {
  for (int c = 0; c < ks * kt; ++c) obs[c] = 0;
  for (std::size_t q = 0; q < pi.size(); ++q) {
    const int b = band(std::fabs(t[pi[q]] - t[pj[q]]), dt, kt);
    if (b >= 0) obs[sb[q] + b * ks] += 1;
  }
}

}  // namespace

int clab_knxnr(const double *x, const double *y, const double *t, int n, const double *ds, int ks,
               const double *dt, int kt, int nsim, unsigned int seed, double *out) {
  if (n < 3 || ks < 1 || kt < 1 || nsim < 0) return CLAB_EARG;
  if (!edges(ds, ks) || !edges(dt, kt)) return CLAB_EDOM;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(x[i]) || !std::isfinite(y[i]) || !std::isfinite(t[i])) {
      return CLAB_EDOM;
    }
  }
  const int m = ks * kt;
  // the pairs within the largest distance, each once (i < j)
  std::vector<int> rp(n + 1), iw(2 * n + 2);
  int nnz = 0;
  clab_swband(x, y, n, ds[ks - 1], 0, 0, rp.data(), nullptr, nullptr, &nnz, iw.data());
  std::vector<int> ci(nnz > 0 ? nnz : 1);
  std::vector<double> d(nnz > 0 ? nnz : 1);
  clab_swband(x, y, n, ds[ks - 1], 0, nnz, rp.data(), ci.data(), d.data(), &nnz, iw.data());
  std::vector<int> pi, pj, sb;
  std::vector<double> scount(ks, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int q = rp[i]; q < rp[i + 1]; ++q) {
      if (ci[q] <= i) continue;
      const int s = band(d[q], ds, ks);
      pi.push_back(i);
      pj.push_back(ci[q]);
      sb.push_back(s);
      scount[s] += 1;
    }
  }
  // the pairs in each time band among all n (n - 1) / 2: a sliding window
  // over the sorted times for each edge
  std::vector<double> u(t, t + n);
  std::sort(u.begin(), u.end());
  std::vector<double> tcount(kt, 0.0);
  double prev = 0;
  for (int b = 0; b < kt; ++b) {
    std::int64_t c = 0;
    int j = 0;
    for (int i = 0; i < n; ++i) {
      if (j < i + 1) j = i + 1;
      while (j < n && u[j] - u[i] <= dt[b]) ++j;
      c += j - i - 1;
    }
    tcount[b] = static_cast<double>(c) - prev;
    prev = static_cast<double>(c);
  }
  const double np = 0.5 * n * (n - 1.0);
  double *obs = out, *ex = out + m, *p = out + 2 * m;
  cells(pi, pj, sb, t, dt, ks, kt, obs);
  for (int b = 0; b < kt; ++b) {
    for (int s = 0; s < ks; ++s) ex[s + b * ks] = scount[s] * tcount[b] / np;
  }
  out[3 * m] = static_cast<double>(pi.size());
  if (nsim == 0) {
    for (int c = 0; c < m; ++c) p[c] = kutil::nan();
    return CLAB_OK;
  }
  // times shuffled among the locations (Fisher-Yates), nsim times
  std::vector<double> tt(t, t + n), sim(m), ge(m, 0.0);
  kutil::rng g(seed);
  for (int r = 0; r < nsim; ++r) {
    for (int i = n - 1; i > 0; --i) std::swap(tt[i], tt[g.below(i + 1)]);
    cells(pi, pj, sb, tt.data(), dt, ks, kt, sim.data());
    for (int c = 0; c < m; ++c) ge[c] += sim[c] >= obs[c];
  }
  for (int c = 0; c < m; ++c) p[c] = (1 + ge[c]) / (nsim + 1.0);
  return CLAB_OK;
}
