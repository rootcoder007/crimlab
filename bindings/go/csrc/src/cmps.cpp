// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, forensics: the congruent matching profile segments (CMPS)
// score of two bullet land signatures (Chen et al. 2019), with the rules
// of the open implementation of Ju and Hofmann (2022) for the details the
// paper leaves open: the comparison signature is cut, from its first to
// its last value, into basis segments; each is widened about its centre
// to the length of each level (cut at the ends of the signature); its
// cross-correlation with the reference is taken over every lag that leaves
// enough points in common (the smaller of 0.9 of the segment's points and
// a tenth of the reference's points, rounded); the highest local maxima
// are its peaks; a peak of the last level that every level confirms within
// the tolerance is a candidate; the position with the most segments voting
// (a candidate within the tolerance) is the congruent one, the middle of
// the positions with that many votes, each candidate counted as often as
// it occurs. With full, a lag is used only when the widened segment lies
// wholly inside the reference.
#include <algorithm>
#include <climits>
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

namespace {

// The correlation at each start j of y (n points) along x (the reference
// with its missing ends cut and padded with missing values), as the
// implementation computes it: running sums over the window, then the sums
// of the pairs with a value missing on one side taken out; NaN with fewer
// than mo complete pairs or no spread
void ccor(const std::vector<double> &x, const std::vector<double> &y, int mo, std::vector<double> &c) {
  const int nx = static_cast<int>(x.size()), ny = static_cast<int>(y.size());
  c.assign(static_cast<std::size_t>(std::max(0, nx - ny + 1)), NAN);
  double sx = 0, sx2 = 0, sy = 0, sy2 = 0;
  for (int i = 0; i < ny; ++i) {
    if (!std::isnan(y[static_cast<std::size_t>(i)])) {
      sy += y[static_cast<std::size_t>(i)];
      sy2 += y[static_cast<std::size_t>(i)] * y[static_cast<std::size_t>(i)];
    }
    if (i < nx && !std::isnan(x[static_cast<std::size_t>(i)])) {
      sx += x[static_cast<std::size_t>(i)];
      sx2 += x[static_cast<std::size_t>(i)] * x[static_cast<std::size_t>(i)];
    }
  }
  for (int j = 0; j < nx - ny + 1; ++j) {
    if (j != 0) {
      const double o = x[static_cast<std::size_t>(j - 1)], e = x[static_cast<std::size_t>(j - 1 + ny)];
      if (!std::isnan(o)) {
        sx -= o;
        sx2 -= o * o;
      }
      if (!std::isnan(e)) {
        sx += e;
        sx2 += e * e;
      }
    }
    double tx = sx, tx2 = sx2, ty = sy, ty2 = sy2, sxy = 0, a0 = 0, b0 = 0;
    int k = 0;
    bool va = false, vb = false;  // the pairs' values vary
    for (int i = 0; i < ny; ++i) {
      const double a = x[static_cast<std::size_t>(i + j)], b = y[static_cast<std::size_t>(i)];
      if (!std::isnan(a)) {
        if (!std::isnan(b)) {
          if (k == 0) {
            a0 = a;
            b0 = b;
          }
          va = va || a != a0;
          vb = vb || b != b0;
          ++k;
          sxy += a * b;
        } else {
          tx -= a;
          tx2 -= a * a;
        }
      } else if (!std::isnan(b)) {
        ty -= b;
        ty2 -= b * b;
      }
    }
    // a constant side has no correlation (cor() gives NA): the running
    // sums would leave rounding error in place of the 0 that says so
    if (k < mo || !va || !vb) continue;
    const double den = std::sqrt((k * tx2 - tx * tx) * (k * ty2 - ty * ty));
    if (den != 0) c[static_cast<std::size_t>(j)] = (sxy * k - tx * ty) / den;
  }
}

// The positions of the npk highest local maxima of c (missing values
// skipped; a value at least its neighbours; one value alone counts only
// when it is 1, as the implementation reads it), the highest first, ties in order; pos[k] is the lag of
// c[k] (NaN where it has none)
void peaks(const std::vector<double> &c, const std::vector<double> &pos, int npk, std::vector<double> &out) {
  out.clear();
  std::vector<int> ok;
  for (std::size_t k = 0; k < c.size(); ++k) {
    if (!std::isnan(c[k])) ok.push_back(static_cast<int>(k));
  }
  std::vector<int> pk;
  const std::size_t n = ok.size();
  auto v = [&c, &ok](std::size_t i) { return c[static_cast<std::size_t>(ok[i])]; };
  if (n == 1) {
    if (v(0) >= 1 && v(0) < 2) pk.push_back(ok[0]);
  } else if (n > 1) {
    for (std::size_t i = 0; i < n; ++i) {
      const bool l = i == 0 || v(i) >= v(i - 1), r = i == n - 1 || v(i) >= v(i + 1);
      if (l && r) pk.push_back(ok[i]);
    }
  }
  std::stable_sort(pk.begin(), pk.end(), [&c](int a, int b) { return c[static_cast<std::size_t>(a)] > c[static_cast<std::size_t>(b)]; });
  for (std::size_t q = 0; q < pk.size() && static_cast<int>(q) < npk; ++q) {
    const double p = pos[static_cast<std::size_t>(pk[q])];
    if (!std::isnan(p)) out.push_back(p);
  }
}

}  // namespace

extern "C" {

int clab_cmpsc(const double *a, int na, const double *b, int nb, int len, const int *npk, int nlev,
               int tx, int full, int *off, int *nseg, int *score, double *pos, double *ccp, int *cmps) CLAB_TRY {
  if (na < 1 || nb < 1 || len < 2 || nlev < 1 || nlev > 8 || tx < 0) return CLAB_EARG;
  for (int l = 0; l < nlev; ++l) {
    if (npk[l] < 1) return CLAB_EARG;
  }
  for (int i = 0; i < na; ++i) {
    if (std::isinf(a[i])) return CLAB_EDOM;
  }
  int b0 = 0, b1 = nb - 1, nbv = 0;
  for (int i = 0; i < nb; ++i) {
    if (std::isinf(b[i])) return CLAB_EDOM;
    nbv += !std::isnan(b[i]);
  }
  while (b0 < nb && std::isnan(b[b0])) ++b0;
  while (b1 >= b0 && std::isnan(b[b1])) --b1;
  int f0 = 0, f1 = na - 1;
  while (f0 < na && std::isnan(a[f0])) ++f0;
  while (f1 >= f0 && std::isnan(a[f1])) --f1;
  if (f1 - f0 + 1 < len) return CLAB_EARG;
  const int ns = (f1 - f0 + 1) / len;
  *nseg = ns;
  *off = f0;
  std::vector<std::vector<double>> cand(static_cast<std::size_t>(ns));
  std::vector<std::vector<double>> lv(static_cast<std::size_t>(nlev));
  std::vector<double> y, xx, c, lag;
  for (int s = 0; s < ns; ++s) {
    const int i0 = f0 + s * len;
    for (int l = 0; l < nlev; ++l) {
      const double ml = static_cast<double>(len) * std::pow(2.0, l);
      if (ml > 1e9) return CLAB_EARG;
      const int tc = static_cast<int>(ml) - len;
      const int lo = std::max(i0 - tc / 2, 0), hi = std::min(i0 + len - 1 + (tc + 1) / 2, na - 1);
      const int ny = hi - lo + 1;
      int nyv = 0;
      for (int i = lo; i <= hi; ++i) nyv += !std::isnan(a[i]);
      if (kutil::stop(static_cast<double>(ny) * nb)) return CLAB_EINTR;
      lv[static_cast<std::size_t>(l)].clear();
      c.clear();
      lag.clear();
      if (full) {
        // every start r of the segment wholly inside b; lag r - lo
        for (int r = 0; r + ny <= nb; ++r) {
          y.assign(a + lo, a + hi + 1);
          xx.assign(b + r, b + r + ny);
          std::vector<double> one;
          ccor(xx, y, 2, one);
          c.push_back(one[0]);
          lag.push_back(r - lo);
        }
      } else {
        if (ny > nb) return CLAB_EARG;
        const double mo = std::min(0.9 * nyv, std::nearbyint(0.1 * nbv));
        const int lagn = nb + 2 * (ny - static_cast<int>(std::ceil(mo))) - ny + 1;
        if (nyv > 0 && nbv > 0 && !(ny < mo) && lagn > 0) {
          int y0 = lo, y1 = hi;
          while (std::isnan(a[y0])) ++y0;
          while (std::isnan(a[y1])) --y1;
          const int nyy = y1 - y0 + 1, pad = static_cast<int>(nyy - mo);
          y.assign(a + y0, a + y1 + 1);
          xx.assign(static_cast<std::size_t>(pad), NAN);
          xx.insert(xx.end(), b + b0, b + b1 + 1);
          xx.insert(xx.end(), static_cast<std::size_t>(pad), NAN);
          std::vector<double> cc;
          ccor(xx, y, static_cast<int>(mo), cc);
          c.assign(static_cast<std::size_t>(b0 + (hi - y1)), NAN);
          c.insert(c.end(), cc.begin(), cc.end());
          c.insert(c.end(), static_cast<std::size_t>((nb - 1 - b1) + (y0 - lo)), NAN);
          // the lag of entry k, k - (ny - mo) less the segment's start;
          // entries past the lags have none
          for (std::size_t k = 0; k < c.size(); ++k) {
            lag.push_back(static_cast<int>(k) < lagn ? static_cast<double>(k) - (ny - mo) - lo : NAN);
          }
        }
      }
      peaks(c, lag, npk[l], lv[static_cast<std::size_t>(l)]);
    }
    std::vector<double> &cs = cand[static_cast<std::size_t>(s)];
    if (nlev == 1) {
      cs = lv[0];
      continue;
    }
    for (double p : lv[static_cast<std::size_t>(nlev - 1)]) {
      bool ok = true;
      for (int l = 0; l < nlev - 1 && ok; ++l) {
        ok = std::any_of(lv[static_cast<std::size_t>(l)].begin(), lv[static_cast<std::size_t>(l)].end(),
                         [p, tx](double q) { return std::fabs(q - p) <= tx; });
      }
      if (ok) cs.push_back(p);
    }
  }
  // votes for every candidate, in ascending order, repeats kept; the
  // middle (ceiling) of those with the most
  std::vector<double> all;
  for (const auto &cs : cand) all.insert(all.end(), cs.begin(), cs.end());
  std::stable_sort(all.begin(), all.end());
  int best = 0;
  std::vector<double> top;
  for (double p : all) {
    int v = 0;
    for (const auto &cs : cand) {
      v += std::any_of(cs.begin(), cs.end(), [p, tx](double q) { return std::fabs(q - p) <= tx; });
    }
    if (v > best) {
      best = v;
      top.clear();
    }
    if (v == best) top.push_back(p);
  }
  *score = 0;
  *pos = top.empty() ? NAN : top[(top.size() + 1) / 2 - 1];
  for (int s = 0; s < ns; ++s) {
    const std::vector<double> &cs = cand[static_cast<std::size_t>(s)];
    ccp[s] = NAN;
    cmps[s] = 0;
    for (double q : cs) {
      if (std::isnan(ccp[s]) || (!top.empty() && std::fabs(q - *pos) < std::fabs(ccp[s] - *pos))) ccp[s] = q;
    }
    if (!top.empty() && !std::isnan(ccp[s]) && std::fabs(ccp[s] - *pos) <= tx) {
      cmps[s] = 1;
      ++*score;
    }
  }
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

}  // extern "C"
