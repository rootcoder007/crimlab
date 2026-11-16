// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The permutation test of a treatment effect with the assignment
// permuted within blocks (crimlab.h, blkprm()): every arrangement keeps
// each block's number of treated units, as a blocked design randomises
// them (Rosenbaum 2002, ch. 2; Good 2005, sec. 3.6). Both statistics are
// linear in the blocks' treated sums S_b, T = c0 + sum_b a_b S_b, and the
// observed value and every permuted one are made by that one formula, so
// an arrangement equal to the observed one gives the observed value bit
// for bit. Exact: every arrangement, block by block in lexicographic
// order of the treated positions; otherwise nperm draws of the package's
// generator. Work: O(nperm n) or O(arrangements x blocks).
#include <cmath>
#include <cstddef>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

namespace {

using std::size_t;
inline size_t z(int i) { return static_cast<size_t>(i); }

// The sums of y over every choice of t of the m positions of pos, in
// lexicographic order of the chosen positions, each summed in increasing
// order.
void sums(const double *y, const std::vector<int> &pos, int t, std::vector<double> &out) {
  const int m = static_cast<int>(pos.size());
  out.clear();
  std::vector<int> c(z(t));
  for (int i = 0; i < t; ++i) c[z(i)] = i;
  for (;;) {
    double s = 0;
    for (int i = 0; i < t; ++i) s += y[pos[z(c[z(i)])]];
    out.push_back(s);
    int i = t - 1;
    while (i >= 0 && c[z(i)] == m - t + i) --i;
    if (i < 0) return;
    ++c[z(i)];
    for (int j = i + 1; j < t; ++j) c[z(j)] = c[z(j - 1)] + 1;
  }
}

}  // namespace

extern "C" int clab_blkprm(const double *y, const int *d, const int *blk, int n, int nb, int stat, int alt,
                           int nperm, int exact, unsigned int seed, int ndist, double *st, double *dist) CLAB_TRY {
  if (n < 2 || nb < 1 || stat < 0 || stat > 1 || alt < 0 || alt > 2 || nperm < 1 || ndist < 0) return CLAB_EARG;
  std::vector<std::vector<int>> pos(z(nb));
  std::vector<int> tb(z(nb), 0);
  double n1 = 0, stot = 0, ymax = 0;
  for (int i = 0; i < n; ++i) {
    if (blk[i] < 0 || blk[i] >= nb || (d[i] != 0 && d[i] != 1) || !std::isfinite(y[i])) return CLAB_EARG;
    pos[z(blk[i])].push_back(i);
    tb[z(blk[i])] += d[i];
    n1 += d[i];
    stot += y[i];
    ymax = std::fmax(ymax, std::fabs(y[i]));
  }
  const double n0 = n - n1;
  if (n1 < 1 || n0 < 1) return CLAB_EDOM;
  // T = c0 + sum_b a_b S_b
  std::vector<double> a(z(nb)), sb(z(nb), 0.0), s1(z(nb), 0.0);
  for (int i = 0; i < n; ++i) {
    sb[z(blk[i])] += y[i];
    if (d[i]) s1[z(blk[i])] += y[i];
  }
  double c0 = 0;
  if (stat == 0) {
    c0 = -stot / n0;
    for (int b = 0; b < nb; ++b) a[z(b)] = 1 / n1 + 1 / n0;
  } else {
    for (int b = 0; b < nb; ++b) {
      const double m = static_cast<double>(pos[z(b)].size()), t = tb[z(b)];
      if (t < 1 || t > m - 1) return CLAB_EDOM;
      a[z(b)] = (m / n) * (1 / t + 1 / (m - t));
      c0 -= (m / n) * sb[z(b)] / (m - t);
    }
  }
  auto tval = [&](const double *s) {
    double v = c0;
    for (int b = 0; b < nb; ++b) v += a[z(b)] * s[b];
    return v;
  };
  // the treated sums of the observed arrangement, each in row order
  const double obs = tval(s1.data()), tol = 1e-9 * (std::fabs(obs) + ymax);
  auto ext = [&](double v) {
    if (alt == 0) return std::fabs(v) >= std::fabs(obs) - tol;
    if (alt == 1) return v >= obs - tol;
    return v <= obs + tol;
  };
  double nge = 0, cnt = 0;
  std::vector<double> cur(s1);
  if (exact) {
    std::vector<std::vector<double>> v(z(nb));
    double tot = 1;
    for (int b = 0; b < nb; ++b) {
      sums(y, pos[z(b)], tb[z(b)], v[z(b)]);
      tot *= static_cast<double>(v[z(b)].size());
      if (tot > ndist) return CLAB_EARG;
    }
    std::vector<size_t> at(z(nb), 0);
    for (int b = 0; b < nb; ++b) cur[z(b)] = v[z(b)][0];
    for (;;) {
      if (kutil::stop(nb)) return CLAB_EINTR;
      const double t = tval(cur.data());
      dist[static_cast<size_t>(cnt)] = t;
      cnt += 1;
      if (ext(t)) nge += 1;
      int b = nb - 1;
      while (b >= 0 && at[z(b)] + 1 == v[z(b)].size()) {
        at[z(b)] = 0;
        cur[z(b)] = v[z(b)][0];
        --b;
      }
      if (b < 0) break;
      cur[z(b)] = v[z(b)][++at[z(b)]];
    }
    st[1] = nge / cnt;
  } else {
    if (nperm > ndist) return CLAB_EARG;
    kutil::rng r(seed);
    std::vector<int> ix;
    for (int k = 0; k < nperm; ++k) {
      if (kutil::stop(n)) return CLAB_EINTR;
      for (int b = 0; b < nb; ++b) {
        // the first tb positions of a partial Fisher-Yates shuffle
        ix = pos[z(b)];
        const int m = static_cast<int>(ix.size());
        double s = 0;
        for (int i = 0; i < tb[z(b)]; ++i) {
          const int j = i + r.below(m - i);
          const int q = ix[z(j)];
          ix[z(j)] = ix[z(i)];
          ix[z(i)] = q;
          s += y[q];
        }
        cur[z(b)] = s;
      }
      const double t = tval(cur.data());
      dist[k] = t;
      if (ext(t)) nge += 1;
    }
    cnt = nperm;
    st[1] = (1 + nge) / (1 + cnt);
  }
  st[0] = obs;
  st[2] = cnt;
  st[3] = nge;
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)
