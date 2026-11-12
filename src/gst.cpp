// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, policing: local Getis-Ord G* hot-spot statistics on a regular
// grid, with an exposure (area inside the study window) per cell, an exact
// conditional test under complete spatial randomness, the classic statistic
// with its analytic moments and a conditional permutation test, and
// Benjamini-Hochberg q-values.
#include <cstdint>

#include "crimlab.h"
#include "kutil.h"

namespace {


// Benjamini-Hochberg q-values of the finite p[0..m) (NaN kept), using
// idx (m ints) for the order
void bhq(const double *p, int m, double *q, int *idx) {
  int k = 0;
  for (int i = 0; i < m; ++i) {
    q[i] = kutil::nan();
    if (!std::isnan(p[i])) idx[k++] = i;
  }
  kutil::hsort(idx, k, [p](int a, int b) { return p[a] < p[b] || (p[a] == p[b] && a < b); });
  double run = 1.0;
  for (int j = k - 1; j >= 0; --j) {
    const double v = p[idx[j]] * k / (j + 1.0);
    if (v < run) run = v;
    q[idx[j]] = run;
  }
}

}  // namespace

extern "C" {

int clab_gstr(const double *cnt, const double *ex, int nr, int nc, double rad,
              int nsim, unsigned int seed, double *out, double *work,
              int *iwork) {
  if (nr < 1 || nc < 1 || nsim < 0 || nr > 46340 || nc > 46340) return CLAB_EARG;
  if (!(rad >= 0.0) || rad > 1000.0) return CLAB_EDOM;
  const int nn = nr * nc;
  // neighbourhood offsets: cells whose centres lie within rad (self included)
  const int ri = static_cast<int>(std::floor(rad));
  int *off = iwork;  // 2 per offset, at most (2 ri + 1)^2 offsets
  int ko = 0;
  for (int dr = -ri; dr <= ri; ++dr)
    for (int dc = -ri; dc <= ri; ++dc)
      if (dr * dr + dc * dc <= rad * rad) {
        off[2 * ko] = dr;
        off[2 * ko + 1] = dc;
        ++ko;
      }
  int *in = off + 2 * ko;  // nn: index among the cells inside the window, or -1
  int *idx = in + nn;      // nn: a permutation of 0..m-1 (and scratch for q)
  int *pos = idx + nn;     // nn: pos[idx[t]] = t
  double sx = 0.0, sxx = 0.0, se = 0.0;
  int m = 0;
  for (int i = 0; i < nn; ++i) {
    if (!(ex[i] >= 0.0) || !std::isfinite(ex[i])) return CLAB_EDOM;
    if (ex[i] > 0.0) {
      if (!(cnt[i] >= 0.0) || cnt[i] != std::floor(cnt[i]) || cnt[i] > 9e15) return CLAB_EDOM;
      in[i] = m++;
      sx += cnt[i];
      sxx += cnt[i] * cnt[i];
      se += ex[i];
    } else {
      in[i] = -1;
    }
  }
  if (m < 2) return CLAB_EARG;
  double *gz = out, *pz = out + nn, *qz = out + 2 * nn, *pp = out + 3 * nn;
  double *ze = out + 4 * nn, *pe = out + 5 * nn, *qe = out + 6 * nn;
  const double xbar = sx / m;
  const double s = std::sqrt(sxx / m - xbar * xbar);
  double *vals = work;  // m: the counts of the cells inside, for permutation
  for (int i = 0; i < nn; ++i)
    if (in[i] >= 0) vals[in[i]] = cnt[i];
  kutil::rng r(seed);
  for (int j = 0; j < m; ++j) idx[j] = pos[j] = j;
  for (int i = 0; i < nn; ++i) {
    for (int c = 0; c < 7; ++c) out[c * nn + i] = kutil::nan();
    if (in[i] < 0) continue;
    const int ir = i % nr, ic = i / nr;
    double g = 0.0, w = 0.0, pe_sum = 0.0;
    for (int o = 0; o < ko; ++o) {
      const int rr = ir + off[2 * o], cc = ic + off[2 * o + 1];
      if (rr < 0 || rr >= nr || cc < 0 || cc >= nc) continue;
      const int j = rr + nr * cc;
      if (in[j] < 0) continue;
      g += cnt[j];
      w += 1.0;
      pe_sum += ex[j];
    }
    // classic G* (Ord and Getis 1995, binary weights, self included):
    // (G - W xbar) / (s sqrt((m W - W^2) / (m - 1)))
    const double den = s * std::sqrt((m * w - w * w) / (m - 1.0));
    if (den > 0.0) {
      gz[i] = (g - w * xbar) / den;
      pz[i] = std::erfc(std::fabs(gz[i]) / std::sqrt(2.0));
    }
    // exposure test: under complete spatial randomness, given the total n,
    // the points in the neighbourhood are Binomial(n, P), P its share of
    // the exposure; z and the exact two-sided p (twice the smaller tail)
    const double pr = pe_sum / se;
    const double mu = sx * pr, sd = std::sqrt(sx * pr * (1.0 - pr));
    if (sd > 0.0) ze[i] = (g - mu) / sd;
    {
      const double t = g >= mu ? kutil::btail(g, sx, pr, true) : kutil::btail(g, sx, pr, false);
      pe[i] = 2.0 * t < 1.0 ? 2.0 * t : 1.0;
    }
    // conditional permutation of the classic G*: the cell keeps its count,
    // its w - 1 neighbours draw w - 1 of the other m - 1 counts at random
    // (a partial Fisher-Yates shuffle of idx); two-sided on |G - W xbar|
    if (nsim > 0 && den > 0.0) {
      const int self = in[i], kk = static_cast<int>(w) - 1;
      const double obs = std::fabs(g - w * xbar);
      // the permutation idx is kept across draws; the cell's own index is
      // first moved to the last place, the draws use places 0..m-2
      auto swp = [&](int a, int b) {
        const int t = idx[a];
        idx[a] = idx[b];
        idx[b] = t;
        pos[idx[a]] = a;
        pos[idx[b]] = b;
      };
      swp(pos[self], m - 1);
      int ge = 0;
      for (int b = 0; b < nsim; ++b) {
        double gs = cnt[i];
        for (int t = 0; t < kk; ++t) {
          swp(t, t + static_cast<int>(r.below(m - 1 - t)));
          gs += vals[idx[t]];
        }
        if (std::fabs(gs - w * xbar) >= obs * (1.0 - 64 * std::numeric_limits<double>::epsilon())) ++ge;
      }
      pp[i] = (1.0 + ge) / (nsim + 1.0);
    }
  }
  bhq(pz, nn, qz, idx);
  bhq(pe, nn, qe, idx);
  return CLAB_OK;
}

}  // extern "C"
