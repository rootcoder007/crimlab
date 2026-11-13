// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, calibration audits of a risk tool: whether the order in
// which an algorithm ranks areas matches the order of what was then
// recorded there, and whether the disagreement follows the areas' group
// (clab_calaud: per-record data rolled up to areas, average ranks, the
// Spearman correlation with its AS 89 p-value, the mean rank gap of each
// group and its permutation test); and a risk score's spread across groups
// (clab_scrdsp: means, medians, gaps with their standard errors, the
// one-way analysis of variance and Welch's test).
#include <climits>
#include <cmath>

#include "crimlab.h"
#include "kutil.h"

namespace {

// Average ranks of x[0..n) with rank 1 the largest: tied values share the
// mean of the positions they occupy. idx: n ints of workspace. Returns
// whether any value is tied.
bool rdesc(const double *x, int n, int *idx, double *r) {
  for (int i = 0; i < n; ++i) idx[i] = i;
  kutil::hsort(idx, n, [x](int a, int b) { return x[a] == x[b] ? a < b : x[a] > x[b]; });
  bool tie = false;
  for (int i = 0; i < n;) {
    int j = i;
    while (j + 1 < n && x[idx[j + 1]] == x[idx[i]]) ++j;
    if (j > i) tie = true;
    const double v = 0.5 * (i + j) + 1.0;
    for (int k = i; k <= j; ++k) r[idx[k]] = v;
    i = j + 1;
  }
  return tie;
}

// The Pearson correlation of a and b, n > 1; NaN when either is constant.
double pear(const double *a, const double *b, int n) {
  kutil::ksum sa, sb;
  for (int i = 0; i < n; ++i) {
    sa.add(a[i]);
    sb.add(b[i]);
  }
  const double ma = sa.val() / n, mb = sb.val() / n;
  kutil::ksum xy, xx, yy;
  for (int i = 0; i < n; ++i) {
    const double u = a[i] - ma, v = b[i] - mb;
    xy.add(u * v);
    xx.add(u * u);
    yy.add(v * v);
  }
  if (!(xx.val() > 0.0) || !(yy.val() > 0.0)) return kutil::nan();
  // one square root of the product: rankings without ties have equal sums
  // of squares, and +-1 comes out exactly
  const double r = xy.val() / std::sqrt(xx.val() * yy.val());
  return std::fmax(-1.0, std::fmin(1.0, r));
}

// P(S >= is) (up) or P(S < is) for Spearman's S = sum (i - p_i)^2 under a
// uniformly random permutation p of n > 1 ranks, by Best and Roberts'
// algorithm AS 89 (1975, doi:10.2307/2347111): exact by enumerating the n!
// permutations for n <= 9 (AS 89 itself enumerates to 6), the Edgeworth
// series otherwise.
double prho(int n, double is, bool up) {
  if (is <= 0.0) return up ? 1.0 : 0.0;
  const double nn = n, n3 = nn * (nn * nn - 1.0) / 3.0;  // the largest S
  if (is > n3) return up ? 0.0 : 1.0;
  if (n <= 9) {
    // Heap's algorithm visits every permutation once
    int p[9], c[9];
    for (int i = 0; i < n; ++i) {
      p[i] = i;
      c[i] = 0;
    }
    auto ss = [&]() {
      double s = 0.0;
      for (int i = 0; i < n; ++i) s += static_cast<double>((i - p[i]) * (i - p[i]));
      return s;
    };
    double ge = ss() >= is ? 1.0 : 0.0, tot = 1.0;
    for (int i = 1; i < n;) {
      if (c[i] < i) {
        const int j = i % 2 == 0 ? 0 : c[i];
        const int t = p[j];
        p[j] = p[i];
        p[i] = t;
        ge += ss() >= is ? 1.0 : 0.0;
        tot += 1.0;
        ++c[i];
        i = 1;
      } else {
        c[i] = 0;
        ++i;
      }
    }
    return up ? ge / tot : (tot - ge) / tot;
  }
  // the Edgeworth series of AS 89, with its continuity correction
  const double c1 = 0.2274, c2 = 0.2531, c3 = 0.1745, c4 = 0.0758, c5 = 0.1033,
               c6 = 0.3932, c7 = 0.0879, c8 = 0.0151, c9 = 0.0072, c10 = 0.0831,
               c11 = 0.0131, c12 = 4.6e-4;
  const double b = 1.0 / nn;
  const double x = (6.0 * (is - 1.0) * b / (nn * nn - 1.0) - 1.0) * std::sqrt(1.0 / b - 1.0);
  double y = x * x;
  const double u =
      x * b *
      (c1 + b * (c2 + c3 * b) +
       y * (-c4 + b * (c5 + c6 * b) -
            y * b * (c7 + c8 * b - y * (c9 - c10 * b + y * b * (c11 - c12 * y)))));
  y = u / std::exp(0.5 * y);
  const double q = 0.5 * std::erfc((up ? x : -x) / std::sqrt(2.0));
  const double pv = (up ? y : -y) + q;
  return std::fmax(0.0, std::fmin(1.0, pv));
}

// The two-sided p-value of Spearman's statistic s (rounded to the even
// whole number it is) for n ranks without ties: twice the smaller tail.
double p2(int n, double s) {
  const double q = std::round(s), nn = n;
  const double p = q > nn * (nn * nn - 1.0) / 6.0 ? prho(n, q, true)
                                                  : prho(n, q + 2.0, false);
  return std::fmin(1.0, 2.0 * p);
}

}  // namespace

extern "C" {

int clab_calaud(const int *area, const double *risk, const double *outc,
                const int *grp, const double *pop, int n, int na, int g,
                int nperm, unsigned int seed, double *ar, double *gr,
                double *st, int *iwork, double *work) {
  if (n < 1 || na < 3 || na > n || g < 1 || nperm < 0) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (area[i] < 0 || area[i] >= na || grp[i] < 0 || grp[i] >= g) return CLAB_EARG;
    if (!std::isfinite(risk[i]) || !std::isfinite(outc[i])) return CLAB_EDOM;
    if (pop != nullptr && !kutil::pos(pop[i])) return CLAB_EDOM;
  }
  const double nan = kutil::nan();
  // ar[k * na + a], k: records, mean risk, outcome rate, group, risk rank,
  // outcome rank, rank gap
  auto A = [&](int k, int a) -> double & { return ar[k * static_cast<long long>(na) + a]; };
  // records sorted by area and group: the runs give each area's majority
  int *idx = iwork, *lab = iwork + n, *ix = lab + na;
  for (int i = 0; i < n; ++i) idx[i] = i;
  kutil::hsort(idx, n, [area, grp](int a, int b) {
    return area[a] == area[b] ? (grp[a] == grp[b] ? a < b : grp[a] < grp[b]) : area[a] < area[b];
  });
  for (int a = 0; a < na; ++a) A(0, a) = 0.0;
  for (int q = 0; q < n;) {
    const int a = area[idx[q]];
    kutil::ksum sr, so;
    int best = -1, bestn = 0, r = q;
    double pa = pop == nullptr ? 0.0 : pop[idx[q]];
    while (r < n && area[idx[r]] == a) {
      int t = r;
      while (t < n && area[idx[t]] == a && grp[idx[t]] == grp[idx[r]]) {
        const int i = idx[t];
        sr.add(risk[i]);
        so.add(outc[i]);
        if (pop != nullptr && pop[i] != pa) return CLAB_EDOM;  // varies within an area
        ++t;
      }
      // the most frequent group; on ties the first in group order
      if (t - r > bestn) {
        bestn = t - r;
        best = grp[idx[r]];
      }
      r = t;
    }
    const double m = r - q;
    A(0, a) = m;
    A(1, a) = sr.val() / m;
    A(2, a) = so.val() / (pop == nullptr ? m : pa);
    A(3, a) = best;
    lab[a] = best;
    q = r;
  }
  for (int a = 0; a < na; ++a) {
    if (!(A(0, a) > 0.0)) return CLAB_EARG;  // an area without records
  }
  const bool t1 = rdesc(ar + na, na, ix, ar + 4 * static_cast<long long>(na));
  const bool t2 = rdesc(ar + 2 * static_cast<long long>(na), na, ix,
                        ar + 5 * static_cast<long long>(na));
  // twice the gap is a whole number (ranks are multiples of one half), so
  // the group sums the permutation test compares are exact
  double *g2 = work, *sg = work + na, *sp = sg + g;
  for (int a = 0; a < na; ++a) {
    A(6, a) = A(5, a) - A(4, a);
    g2[a] = 2.0 * A(6, a);
  }
  const double rho = pear(ar + 4 * static_cast<long long>(na), ar + 5 * static_cast<long long>(na), na);
  const double nn = na, s = (nn * nn * nn - nn) * (1.0 - rho) / 6.0;
  st[0] = na;
  st[1] = rho;
  st[2] = s;
  st[3] = (t1 || t2 || na >= 1290 || std::isnan(rho)) ? nan : p2(na, s);
  st[4] = (t1 || t2) ? 1.0 : 0.0;
  // the group means, the between-group sum of squares of the gaps
  // T = sum_g n_g mean_g^2 (the gaps sum to 0) and their permutation tests
  for (int j = 0; j < g; ++j) {
    gr[j] = 0.0;
    sg[j] = 0.0;
  }
  for (int a = 0; a < na; ++a) {
    gr[lab[a]] += 1.0;
    sg[lab[a]] += g2[a];
  }
  auto tstat = [&](const double *x) {
    kutil::ksum t;
    for (int j = 0; j < g; ++j) {
      if (gr[j] > 0.0) t.add(x[j] * x[j] / gr[j]);
    }
    return 0.25 * t.val();
  };
  const double t0 = tstat(sg);
  for (int j = 0; j < g; ++j) {
    gr[g + j] = gr[j] > 0.0 ? 0.5 * sg[j] / gr[j] : nan;
    gr[2 * g + j] = nan;
  }
  st[5] = t0;
  st[6] = nan;
  if (nperm > 0) {
    kutil::rng gen(seed);
    double *cnt = sp + g, ct = 0.0;
    for (int j = 0; j < g; ++j) cnt[j] = 0.0;
    for (int a = 0; a < na; ++a) ix[a] = lab[a];
    for (int b = 0; b < nperm; ++b) {
      // Fisher-Yates on the area labels
      for (int a = na - 1; a > 0; --a) {
        const int k = gen.below(a + 1);
        const int t = ix[a];
        ix[a] = ix[k];
        ix[k] = t;
      }
      for (int j = 0; j < g; ++j) sp[j] = 0.0;
      for (int a = 0; a < na; ++a) sp[ix[a]] += g2[a];
      for (int j = 0; j < g; ++j) cnt[j] += std::fabs(sp[j]) >= std::fabs(sg[j]) ? 1.0 : 0.0;
      // T up to rounding: the same arrangement of sums gives the same T
      ct += tstat(sp) >= t0 * (1.0 - 1e-12) ? 1.0 : 0.0;
    }
    for (int j = 0; j < g; ++j) {
      if (gr[j] > 0.0) gr[2 * g + j] = (1.0 + cnt[j]) / (1.0 + nperm);
    }
    st[6] = (1.0 + ct) / (1.0 + nperm);
  }
  return CLAB_OK;
}

int clab_scrdsp(const double *x, const int *grp, int n, int g, int ref,
                double *gs, double *st, int *iwork) {
  if (n < 1 || g < 2 || ref < -1 || ref >= g) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (grp[i] < 0 || grp[i] >= g) return CLAB_EARG;
    if (!std::isfinite(x[i])) return CLAB_EDOM;
  }
  const double nan = kutil::nan();
  // gs[k * g + j], k: n, mean, median, sd, gap, standard error of the gap
  auto G = [&](int k, int j) -> double & { return gs[k * g + j]; };
  int *idx = iwork;
  for (int i = 0; i < n; ++i) idx[i] = i;
  kutil::hsort(idx, n, [x, grp](int a, int b) {
    return grp[a] == grp[b] ? (x[a] == x[b] ? a < b : x[a] < x[b]) : grp[a] < grp[b];
  });
  kutil::ksum all;
  for (int i = 0; i < n; ++i) all.add(x[i]);
  double mall = all.val() / n;
  {
    kutil::ksum c;
    for (int i = 0; i < n; ++i) c.add(x[i] - mall);
    mall += c.val() / n;
  }
  kutil::ksum ssb, ssw;
  int k = 0, q = 0;
  bool welch = true;
  for (int j = 0; j < g; ++j) {
    int r = q;
    while (r < n && grp[idx[r]] == j) ++r;
    const int m = r - q;
    if (m == 0) return CLAB_EDOM;  // a group without cases
    ++k;
    kutil::ksum s;
    for (int t = q; t < r; ++t) s.add(x[idx[t]]);
    double mu = s.val() / m;
    kutil::ksum c, v;
    for (int t = q; t < r; ++t) c.add(x[idx[t]] - mu);
    mu += c.val() / m;
    for (int t = q; t < r; ++t) v.add((x[idx[t]] - mu) * (x[idx[t]] - mu));
    G(0, j) = m;
    G(1, j) = mu;
    G(2, j) = m % 2 == 1 ? x[idx[q + m / 2]] : 0.5 * (x[idx[q + m / 2 - 1]] + x[idx[q + m / 2]]);
    G(3, j) = m > 1 ? std::sqrt(v.val() / (m - 1)) : nan;
    if (!(m > 1 && v.val() > 0.0)) welch = false;
    ssb.add(m * (mu - mall) * (mu - mall));
    ssw.add(v.val());
    q = r;
  }
  int b = ref;
  if (b < 0) {
    b = 0;
    for (int j = 1; j < g; ++j) {
      if (G(1, j) < G(1, b)) b = j;
    }
  }
  const double df1 = k - 1.0, df2 = n - static_cast<double>(k);
  const double msw = df2 > 0.0 ? ssw.val() / df2 : nan;
  double lo = G(1, 0), hi = G(1, 0);
  for (int j = 0; j < g; ++j) {
    G(4, j) = G(1, j) - G(1, b);
    G(5, j) = j == b ? nan : std::sqrt(msw * (1.0 / G(0, j) + 1.0 / G(0, b)));
    lo = std::fmin(lo, G(1, j));
    hi = std::fmax(hi, G(1, j));
  }
  // the analysis of variance F, equal variances assumed, and Welch's
  // (1951, doi:10.1093/biomet/38.3-4.330) for unequal ones
  st[0] = df2 > 0.0 ? (ssb.val() / df1) / msw : nan;
  st[1] = df1;
  st[2] = df2;
  st[3] = st[4] = st[5] = nan;
  if (welch) {
    kutil::ksum sw, swm;
    for (int j = 0; j < g; ++j) {
      const double w = G(0, j) / (G(3, j) * G(3, j));
      sw.add(w);
      swm.add(w * G(1, j));
    }
    const double m = swm.val() / sw.val();
    kutil::ksum a, t;
    for (int j = 0; j < g; ++j) {
      const double w = G(0, j) / (G(3, j) * G(3, j));
      a.add(w * (G(1, j) - m) * (G(1, j) - m));
      t.add((1.0 - w / sw.val()) * (1.0 - w / sw.val()) / (G(0, j) - 1.0));
    }
    const double tmp = t.val() / (df1 * (df1 + 2.0));  // k^2 - 1
    st[3] = a.val() / (df1 * (1.0 + 2.0 * (df1 - 1.0) * tmp));
    st[4] = df1;
    st[5] = 1.0 / (3.0 * tmp);
  }
  st[6] = hi - lo;
  st[7] = b;
  st[8] = msw;
  return CLAB_OK;
}

}  // extern "C"
