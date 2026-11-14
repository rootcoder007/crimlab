// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, forensics: the congruent matching profile segments (CMPS)
// score of two bullet land signatures (Chen et al. 2019; Ju and Hofmann
// 2022). The comparison signature is cut into basis segments; each finds
// its best registration positions on the reference by the peaks of its
// cross-correlation, at its own length and at doubled lengths about its
// centre; the position most segments agree on is the congruent one, and
// the segments that agree with it are counted.
#include <algorithm>
#include <climits>
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

namespace {

// Pearson correlation of a[i + t] and b[j + t] for t in [0, m) with j + t
// inside b, over the pairs with neither value missing (NaN); NaN with
// fewer than two pairs or no spread
double cor(const double *a, int i, const double *b, int nb, int j, int m) {
  const int t0 = std::max(0, -j), t1 = std::min(m, nb - j);
  double sa = 0, sb = 0;
  int k = 0;
  for (int t = t0; t < t1; ++t) {
    const double x = a[i + t], y = b[j + t];
    if (std::isnan(x) || std::isnan(y)) continue;
    sa += x;
    sb += y;
    ++k;
  }
  if (k < 2) return NAN;
  const double ma = sa / k, mb = sb / k;
  double sxy = 0, sxx = 0, syy = 0;
  for (int t = t0; t < t1; ++t) {
    const double x = a[i + t], y = b[j + t];
    if (std::isnan(x) || std::isnan(y)) continue;
    sxy += (x - ma) * (y - mb);
    sxx += (x - ma) * (x - ma);
    syy += (y - mb) * (y - mb);
  }
  if (!(sxx > 0) || !(syy > 0)) return NAN;
  return sxy / std::sqrt(sxx * syy);
}

// The lags (reference start minus segment start) of the k highest local
// maxima of the cross-correlation of a[i..i+m) with b, over the lags that
// keep at least ov points of the segment inside b; their values in val
int peaks(const double *a, int i, int m, const double *b, int nb, int k, int ov, std::vector<int> &pos,
          std::vector<double> &val) {
  pos.clear();
  val.clear();
  const int j0 = ov - m, nj = nb - ov - j0 + 1;
  if (nj < 1) return CLAB_OK;
  std::vector<double> c(static_cast<std::size_t>(nj));
  for (int q = 0; q < nj; ++q) c[static_cast<std::size_t>(q)] = cor(a, i, b, nb, j0 + q, m);
  if (kutil::stop(static_cast<double>(nj) * m)) return CLAB_EINTR;
  std::vector<int> pk;
  for (int q = 0; q < nj; ++q) {
    const double v = c[static_cast<std::size_t>(q)];
    if (std::isnan(v)) continue;
    const bool l = q == 0 || std::isnan(c[static_cast<std::size_t>(q - 1)]) || v > c[static_cast<std::size_t>(q - 1)];
    const bool r = q == nj - 1 || std::isnan(c[static_cast<std::size_t>(q + 1)]) || v >= c[static_cast<std::size_t>(q + 1)];
    if (l && r) pk.push_back(q);
  }
  std::stable_sort(pk.begin(), pk.end(), [&c](int x, int y) { return c[static_cast<std::size_t>(x)] > c[static_cast<std::size_t>(y)]; });
  for (std::size_t q = 0; q < pk.size() && static_cast<int>(q) < k; ++q) {
    pos.push_back(j0 + pk[q] - i);
    val.push_back(c[static_cast<std::size_t>(pk[q])]);
  }
  return CLAB_OK;
}

}  // namespace

extern "C" {

int clab_cmpsc(const double *a, int na, const double *b, int nb, int len, const int *npk, int nlev,
               int tx, int *off, int *nseg, int *score, int *pos, int *ccp, int *cmps) {
  if (na < 1 || nb < 1 || len < 2 || nlev < 1 || nlev > 8 || tx < 0) return CLAB_EARG;
  for (int l = 0; l < nlev; ++l) {
    if (npk[l] < 1) return CLAB_EARG;
  }
  for (int i = 0; i < na; ++i) {
    if (std::isinf(a[i])) return CLAB_EDOM;
  }
  for (int i = 0; i < nb; ++i) {
    if (std::isinf(b[i])) return CLAB_EDOM;
  }
  // the comparison from its first to its last value
  int f0 = 0, f1 = na - 1;
  while (f0 < na && std::isnan(a[f0])) ++f0;
  while (f1 >= f0 && std::isnan(a[f1])) --f1;
  if (f1 - f0 + 1 < len) return CLAB_EARG;
  const int ns = (f1 - f0 + 1) / len;
  *nseg = ns;
  *off = f0;
  // each segment's candidate positions: with one level, its top peaks;
  // with more, the peak of the last level (best first) that every other
  // level has a peak within tx of
  std::vector<std::vector<int>> cand(static_cast<std::size_t>(ns));
  std::vector<int> ps;
  std::vector<double> vs;
  std::vector<std::vector<int>> lv(static_cast<std::size_t>(nlev));
  std::vector<char> sk(static_cast<std::size_t>(nlev), 0);
  for (int s = 0; s < ns; ++s) {
    const int i0 = f0 + s * len;
    for (int l = 0; l < nlev; ++l) {
      // the segment at level l: len 2^l points about the basis segment's
      // centre, cut to the signature
      const double m = static_cast<double>(len) * std::pow(2.0, l);
      // a level longer than the signature is left out
      sk[static_cast<std::size_t>(l)] = m > static_cast<double>(f1 - f0 + 1);
      lv[static_cast<std::size_t>(l)].clear();
      if (sk[static_cast<std::size_t>(l)]) continue;
      const int ml = static_cast<int>(m);
      int lo = i0 + len / 2 - ml / 2;
      lo = std::max(f0, std::min(lo, f1 + 1 - ml));
      const int rc = peaks(a, lo, ml, b, nb, npk[l], ml, ps, vs);
      if (rc != CLAB_OK) return rc;
      // a lag d pairs a[t] with b[t + d], whatever the segment's length
      lv[static_cast<std::size_t>(l)] = ps;
    }
    std::vector<int> &cs = cand[static_cast<std::size_t>(s)];
    if (nlev == 1) {
      cs = lv[0];
      continue;
    }
    int last = nlev - 1;
    while (sk[static_cast<std::size_t>(last)]) --last;
    for (int p : lv[static_cast<std::size_t>(last)]) {
      bool ok = true;
      for (int l = 0; l < last && ok; ++l) {
        ok = std::any_of(lv[static_cast<std::size_t>(l)].begin(), lv[static_cast<std::size_t>(l)].end(),
                         [p, tx](int q) { return std::abs(q - p) <= tx; });
      }
      if (ok) {
        cs.push_back(p);
        break;
      }
    }
  }
  // votes: a candidate position gets one from each segment with a
  // candidate within tx; the most voted, the middle one of ties
  std::vector<int> all;
  for (const auto &cs : cand) all.insert(all.end(), cs.begin(), cs.end());
  std::sort(all.begin(), all.end());
  all.erase(std::unique(all.begin(), all.end()), all.end());
  int best = -1;
  std::vector<int> top;
  for (int p : all) {
    int v = 0;
    for (const auto &cs : cand) {
      v += std::any_of(cs.begin(), cs.end(), [p, tx](int q) { return std::abs(q - p) <= tx; });
    }
    if (v > best) {
      best = v;
      top.clear();
    }
    if (v == best) top.push_back(p);
  }
  *score = 0;
  *pos = top.empty() ? INT_MIN : top[(top.size() - 1) / 2];
  for (int s = 0; s < ns; ++s) {
    const std::vector<int> &cs = cand[static_cast<std::size_t>(s)];
    ccp[s] = INT_MIN;
    cmps[s] = 0;
    for (int q : cs) {
      if (ccp[s] == INT_MIN || (top.size() && std::abs(q - *pos) < std::abs(ccp[s] - *pos))) ccp[s] = q;
    }
    if (*pos != INT_MIN && ccp[s] != INT_MIN && std::abs(ccp[s] - *pos) <= tx) {
      cmps[s] = 1;
      ++*score;
    }
  }
  return CLAB_OK;
}

}  // extern "C"
