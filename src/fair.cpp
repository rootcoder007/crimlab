// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, fairness: group confusion-matrix metrics of a binary risk tool
// with Wilson score intervals (clab_fair); differences and ratios between
// groups with their intervals, the bias amplification score and a
// bootstrap (clab_fairdf); the same metrics cell by cell over time and
// place (clab_tmpaud); the Lipschitz audit of individual fairness
// (clab_lipaud); the bounds a noisy recorded base rate puts on the true
// one (clab_nsbnd), the false positive rate Chouldechova's identity implies
// (clab_impfpr) and the pairs a noisy ranking identifies (clab_rnkres).
#include <cmath>
#include <climits>
#include <cstdint>
#include <limits>

#include "crimlab.h"
#include "kutil.h"

namespace {

// Wilson score interval for x successes in m > 0 trials at normal quantile
// z > 0: the set of p with (x - m p)^2 <= z^2 m p (1 - p). With
// d = m + z^2, centre c = (x + z^2 / 2) / d and half-width
// h = (z / d) sqrt(x (m - x) / m + z^2 / 4), the ends are c -+ h. Since
// c^2 - h^2 = x^2 / (m d), the lower end is x^2 / (m d) / (c + h), with no
// cancellation, and by symmetry in x <-> m - x the upper end is 1 minus the
// same expression for m - x. x = 0 gives lo = 0 and x = m gives hi = 1
// exactly.
void wilson(double x, double m, double z, double *lo, double *hi) {
  if (!(m > 0.0)) {
    *lo = *hi = kutil::nan();
    return;
  }
  const double z2 = z * z;
  const double d = m + z2;
  const double y = m - x;
  const double h = (z / d) * std::sqrt(x * y / m + 0.25 * z2);
  const double cx = (x + 0.5 * z2) / d;
  const double cy = (y + 0.5 * z2) / d;
  *lo = x * x / (m * d) / (cx + h);
  *hi = 1.0 - y * y / (m * d) / (cy + h);
}

inline double ratio(double a, double b) {
  return b > 0.0 ? a / b : kutil::nan();
}

// Newcombe's hybrid score interval for p1 - p2 (his method 10, Newcombe
// 1998, doi:10.1002/(SICI)1097-0258(19980430)17:8<873::AID-SIM779>3.0.CO;2-I):
// the Wilson intervals (l1, u1), (l2, u2) combined as
// d - sqrt((p1 - l1)^2 + (u2 - p2)^2), d + sqrt((u1 - p1)^2 + (p2 - l2)^2).
// NaN when either denominator is 0.
void newc(double x1, double m1, double x2, double m2, double z, double *lo,
          double *hi) {
  if (!(m1 > 0.0) || !(m2 > 0.0)) {
    *lo = *hi = kutil::nan();
    return;
  }
  double l1, u1, l2, u2;
  wilson(x1, m1, z, &l1, &u1);
  wilson(x2, m2, z, &l2, &u2);
  const double p1 = x1 / m1, p2 = x2 / m2, d = p1 - p2;
  *lo = d - std::sqrt((p1 - l1) * (p1 - l1) + (u2 - p2) * (u2 - p2));
  *hi = d + std::sqrt((u1 - p1) * (u1 - p1) + (p2 - l2) * (p2 - l2));
}

// The log interval of Katz et al. (1978, doi:10.2307/2530610) for the ratio
// p1 / p2: exp(log(r) -+ z sqrt(1 / x1 - 1 / m1 + 1 / x2 - 1 / m2)); NaN when
// either count of successes is 0 (the logarithm is undefined there).
void katz(double x1, double m1, double x2, double m2, double z, double *lo,
          double *hi) {
  if (!(x1 > 0.0) || !(x2 > 0.0)) {
    *lo = *hi = kutil::nan();
    return;
  }
  const double r = (x1 / m1) / (x2 / m2);
  const double se = std::sqrt(1.0 / x1 - 1.0 / m1 + 1.0 / x2 - 1.0 / m2);
  *lo = r * std::exp(-z * se);
  *hi = r * std::exp(z * se);
}

// The Gini coefficient of x[0..n), x >= 0, as the sorted-rank formula
// 2 sum_i i x_(i) / (n sum x) - (n + 1) / n (the mean absolute difference
// over twice the mean, with n^2 pairs); 0 for n < 2 or an all-zero x.
// Sorts x.
double gini(double *x, int n) {
  if (n < 2) return 0.0;
  kutil::hsort(x, n, [](double a, double b) { return a < b; });
  kutil::ksum t, s;
  for (int i = 0; i < n; ++i) {
    t.add(x[i]);
    s.add((i + 1.0) * x[i]);
  }
  if (!(t.val() > 0.0)) return 0.0;
  return 2.0 * s.val() / (n * t.val()) - (n + 1.0) / n;
}

// The sample quantile of type 7 (R's default) of the sorted x[0..n), n > 0.
double q7(const double *x, int n, double p) {
  const double h = (n - 1) * p;
  const int i = static_cast<int>(std::floor(h));
  if (i >= n - 1) return x[n - 1];
  return x[i] + (h - i) * (x[i + 1] - x[i]);
}

// Of the values v[j], j != b: the one with the largest absolute value (the
// first on ties), NaN values skipped; NaN when none.
double worst(const double *v, int g, int b) {
  double best = kutil::nan();
  for (int j = 0; j < g; ++j) {
    if (j == b || std::isnan(v[j])) continue;
    if (std::isnan(best) || std::fabs(v[j]) > std::fabs(best)) best = v[j];
  }
  return best;
}

// The signed parity gap of largest size against b and the Gini coefficient
// of the g rates r (copied to w, g doubles): the bias amplification score
// is their product.
void bas(const double *r, int g, int b, double *w, double *dp, double *gi) {
  for (int j = 0; j < g; ++j) w[j] = r[j] - r[b];
  *dp = worst(w, g, b);
  for (int j = 0; j < g; ++j) w[j] = r[j];
  *gi = gini(w, g);
}

}  // namespace

extern "C" {

int clab_fair(const double *score, const int *y, const int *grp, int n, int g,
              double cut, int base, double z, double *cnt, double *met,
              double *rat, double *lo, double *hi) {
  if (n < 1 || g < 1 || base < 0 || base >= g) return CLAB_EARG;
  if (std::isnan(cut)) return CLAB_EDOM;
  if (lo != nullptr && hi != nullptr && !kutil::pos(z)) return CLAB_EDOM;
  for (int i = 0; i < n; ++i) {
    if (grp[i] < 0 || grp[i] >= g) return CLAB_EARG;
    if ((y[i] != 0 && y[i] != 1) || std::isnan(score[i])) return CLAB_EDOM;
  }
  for (int j = 0; j < 5 * g; ++j) cnt[j] = 0.0;
  // cnt[0..3][j]: tp, fp, tn, fn of group j; cnt[4][j]: its size. Counts
  // below 2^31 are exact in a double.
  for (int i = 0; i < n; ++i) {
    const int pr = score[i] > cut ? 1 : 0;
    const int k = y[i] == 1 ? (pr == 1 ? 0 : 3) : (pr == 1 ? 1 : 2);
    cnt[k * g + grp[i]] += 1.0;
    cnt[4 * g + grp[i]] += 1.0;
  }
  for (int j = 0; j < g; ++j) {
    const double tp = cnt[j], fp = cnt[g + j], tn = cnt[2 * g + j],
                 fn = cnt[3 * g + j], m = cnt[4 * g + j];
    // successes and trials of the eight proportions, in output order
    const double xs[8] = {tp + tn, tp, tn, fp, fn, tp, tn, tp + fp};
    const double ms[8] = {m, tp + fn, tn + fp, fp + tn, fn + tp,
                          tp + fp, tn + fn, m};
    for (int k = 0; k < 8; ++k) {
      met[k * g + j] = ratio(xs[k], ms[k]);
      if (lo != nullptr && hi != nullptr) {
        wilson(xs[k], ms[k], z, &lo[k * g + j], &hi[k * g + j]);
      }
    }
    met[8 * g + j] = tp + fp;
    // Matthews correlation: the numerator tp tn - fp fn in 64-bit integers
    // (each count < 2^31, each product < 2^62), so it is exact before its
    // one rounding; NaN when a margin is empty (the coefficient is 0 / 0).
    const std::int64_t a = static_cast<std::int64_t>(tp),
                       b = static_cast<std::int64_t>(tn),
                       c = static_cast<std::int64_t>(fp),
                       e = static_cast<std::int64_t>(fn);
    const double num = static_cast<double>(a * b - c * e);
    const double d1 = static_cast<double>((a + c) * (a + e));
    const double d2 = static_cast<double>((b + c) * (b + e));
    met[9 * g + j] = d1 > 0.0 && d2 > 0.0
                         ? num / (std::sqrt(d1) * std::sqrt(d2))
                         : kutil::nan();
  }
  // ratios to the base group: IEEE division (x / 0 is Inf, 0 / 0 is NaN)
  for (int k = 0; k < 10; ++k) {
    for (int j = 0; j < g; ++j) rat[k * g + j] = met[k * g + j] / met[k * g + base];
  }
  return CLAB_OK;
}

int clab_fairdf(const int *pred, const int *y, const int *grp, int n, int g,
                int base, double z, int nboot, unsigned int seed, double *cnt,
                double *met, double *sum, double *work) {
  if (n < 1 || g < 2 || base < -1 || base >= g || nboot < 0) return CLAB_EARG;
  if (!kutil::pos(z)) return CLAB_EDOM;
  for (int i = 0; i < n; ++i) {
    if (grp[i] < 0 || grp[i] >= g) return CLAB_EARG;
    if (pred[i] != 0 && pred[i] != 1) return CLAB_EDOM;
    if (y != nullptr && y[i] != 0 && y[i] != 1) return CLAB_EDOM;
  }
  const double nan = kutil::nan();
  // cnt[k * g + j], k: cases, predicted positive, true positives, of them
  // predicted positive, true negatives, of them predicted positive
  for (int k = 0; k < 6 * g; ++k) cnt[k] = 0.0;
  for (int i = 0; i < n; ++i) {
    const int j = grp[i];
    cnt[j] += 1.0;
    cnt[g + j] += pred[i];
    if (y != nullptr) {
      const int k = y[i] == 1 ? 2 : 4;
      cnt[k * g + j] += 1.0;
      cnt[(k + 1) * g + j] += pred[i];
    }
  }
  double *r = work;
  for (int j = 0; j < g; ++j) {
    if (!(cnt[j] > 0.0)) return CLAB_EDOM;  // a group without cases
    r[j] = cnt[g + j] / cnt[j];
  }
  int b = base;
  if (b < 0) {
    b = 0;
    for (int j = 1; j < g; ++j) {
      if (r[j] > r[b]) b = j;
    }
  }
  for (int k = 0; k < 23 * g; ++k) met[k] = nan;
  auto m = [&](int k, int j) -> double & { return met[k * g + j]; };
  for (int j = 0; j < g; ++j) {
    const double nj = cnt[j], xj = cnt[g + j], nb = cnt[b], xb = cnt[g + b];
    m(0, j) = nj;
    m(1, j) = xj;
    m(2, j) = r[j];
    wilson(xj, nj, z, &m(3, j), &m(4, j));
    m(5, j) = r[j] - r[b];
    m(8, j) = r[j] / r[b];
    if (j != b) {
      newc(xj, nj, xb, nb, z, &m(6, j), &m(7, j));
      katz(xj, nj, xb, nb, z, &m(9, j), &m(10, j));
    }
    if (y == nullptr) continue;
    const double pj = cnt[2 * g + j], tj = cnt[3 * g + j], qj = cnt[4 * g + j],
                 fj = cnt[5 * g + j], pb = cnt[2 * g + b], tb = cnt[3 * g + b],
                 qb = cnt[4 * g + b], fb = cnt[5 * g + b];
    m(11, j) = ratio(tj, pj);
    m(12, j) = ratio(fj, qj);
    m(13, j) = m(11, j) - ratio(tb, pb);
    m(16, j) = m(12, j) - ratio(fb, qb);
    m(19, j) = 0.5 * (m(13, j) + m(16, j));
    const double a = std::fabs(m(13, j)), c = std::fabs(m(16, j));
    m(22, j) = std::isnan(a) ? m(16, j)
               : (std::isnan(c) || a >= c ? m(13, j) : m(16, j));
    if (j == b) continue;
    newc(tj, pj, tb, pb, z, &m(14, j), &m(15, j));
    newc(fj, qj, fb, qb, z, &m(17, j), &m(18, j));
    // MOVER (Zou and Donner 2008, doi:10.1002/sim.3095): the interval of
    // a sum of independent estimates from the intervals of each
    const double d1 = m(13, j), d2 = m(16, j);
    m(20, j) = 0.5 * (d1 + d2 - std::sqrt((d1 - m(14, j)) * (d1 - m(14, j)) +
                                          (d2 - m(17, j)) * (d2 - m(17, j))));
    m(21, j) = 0.5 * (d1 + d2 + std::sqrt((m(15, j) - d1) * (m(15, j) - d1) +
                                          (m(18, j) - d2) * (m(18, j) - d2)));
  }
  // the summary: smallest finite ratio, the parity gap, Gini, BAS, their
  // bootstrap ends, the average odds difference and the equalized-odds gap
  // of largest size, the base group
  double dmin = nan;
  for (int j = 0; j < g; ++j) {
    if (j != b && std::isfinite(m(8, j)) && !(m(8, j) >= dmin)) dmin = m(8, j);
  }
  double *w = work + g;
  sum[0] = dmin;
  bas(r, g, b, w, &sum[1], &sum[2]);
  sum[3] = sum[1] * sum[2];
  for (int k = 4; k < 8; ++k) sum[k] = nan;
  sum[8] = y == nullptr ? nan : worst(met + 19 * g, g, b);
  sum[9] = y == nullptr ? nan : worst(met + 22 * g, g, b);
  sum[10] = b;
  if (nboot > 0) {
    // percentile bootstrap, resampling within each group: n_j draws with
    // replacement from a group's predictions count Binomial(n_j, r_j)
    // positives; the base group is the one chosen above
    kutil::rng gen(seed);
    double *rb = w + g, *gb = rb + g, *bb = gb + nboot;
    for (int t = 0; t < nboot; ++t) {
      for (int j = 0; j < g; ++j) {
        double c = 0.0;
        for (int i = 0; i < static_cast<int>(cnt[j]); ++i) c += gen.unif() < r[j] ? 1.0 : 0.0;
        rb[j] = c / cnt[j];
      }
      double dp, gi;
      bas(rb, g, b, w, &dp, &gi);
      gb[t] = gi;
      bb[t] = dp * gi;
    }
    auto less = [](double a, double c) { return a < c; };
    kutil::hsort(gb, nboot, less);
    kutil::hsort(bb, nboot, less);
    const double al = std::erfc(z / std::sqrt(2.0));
    sum[4] = q7(gb, nboot, 0.5 * al);
    sum[5] = q7(gb, nboot, 1.0 - 0.5 * al);
    sum[6] = q7(bb, nboot, 0.5 * al);
    sum[7] = q7(bb, nboot, 1.0 - 0.5 * al);
  }
  return CLAB_OK;
}

int clab_tmpaud(const int *pred, const int *grp, const int *cell, int n, int g,
                int nc, const int *city, int ncity, int base, double *cout,
                double *yout, double *sum, double *work) {
  if (n < 1 || g < 2 || nc < 1 || ncity < 1 || base < -1 || base >= g) return CLAB_EARG;
  for (int c = 0; c < nc; ++c) {
    if (city[c] < 0 || city[c] >= ncity) return CLAB_EARG;
  }
  for (int i = 0; i < n; ++i) {
    if (grp[i] < 0 || grp[i] >= g || cell[i] < 0 || cell[i] >= nc) return CLAB_EARG;
    if (pred[i] != 0 && pred[i] != 1) return CLAB_EDOM;
  }
  const double nan = kutil::nan();
  // work: cases and positives of group j in cell c at 2 (c g + j), then
  // the rates of one cell (g) and scratch (g)
  double *cg = work, *r = work + 2 * static_cast<long long>(nc) * g, *w = r + g;
  for (long long k = 0; k < 2LL * nc * g; ++k) cg[k] = 0.0;
  for (int i = 0; i < n; ++i) {
    const long long k = 2 * (static_cast<long long>(cell[i]) * g + grp[i]);
    cg[k] += 1.0;
    cg[k + 1] += pred[i];
  }
  int b = base;
  if (b < 0) {
    // the group with the highest pooled rate (the first on ties)
    double best = -1.0;
    for (int j = 0; j < g; ++j) {
      double nn = 0.0, pp = 0.0;
      for (int c = 0; c < nc; ++c) {
        nn += cg[2 * (static_cast<long long>(c) * g + j)];
        pp += cg[2 * (static_cast<long long>(c) * g + j) + 1];
      }
      const double rj = nn > 0.0 ? pp / nn : -1.0;
      if (rj > best) {
        best = rj;
        b = j;
      }
    }
  }
  for (int c = 0; c < nc; ++c) {
    const double *x = cg + 2 * static_cast<long long>(c) * g;
    double nn = 0.0;
    int k = 0;
    for (int j = 0; j < g; ++j) {
      nn += x[2 * j];
      if (x[2 * j] > 0.0) ++k;
    }
    double *o = cout + 6 * static_cast<long long>(c);
    o[0] = nn;
    o[1] = 0.0;
    for (int t = 2; t < 6; ++t) o[t] = nan;
    if (k < 2 || !(x[2 * b] > 0.0)) continue;  // skipped: no comparison
    o[1] = 1.0;
    // the present groups' rates, the base group first
    int m = 0, mb = 0;
    for (int j = 0; j < g; ++j) {
      if (!(x[2 * j] > 0.0)) continue;
      if (j == b) mb = m;
      r[m++] = x[2 * j + 1] / x[2 * j];
    }
    double dmin = nan;
    for (int j = 0; j < m; ++j) {
      const double q = r[j] / r[mb];
      if (j != mb && std::isfinite(q) && !(q >= dmin)) dmin = q;
    }
    o[2] = dmin;
    bas(r, m, mb, w, &o[3], &o[4]);
    o[5] = o[3] * o[4];
  }
  // by city: audited cells, means of the finite values, the range of the
  // finite ratios and how many exceed 1
  for (int y = 0; y < ncity; ++y) {
    double *o = yout + 9 * static_cast<long long>(y);
    kutil::ksum s[4];
    double k[4] = {0, 0, 0, 0}, lo = nan, hi = nan, gt = 0.0, na = 0.0;
    for (int c = 0; c < nc; ++c) {
      const double *e = cout + 6 * static_cast<long long>(c);
      if (city[c] != y || e[1] == 0.0) continue;
      na += 1.0;
      for (int t = 0; t < 4; ++t) {
        if (std::isfinite(e[2 + t])) {
          s[t].add(e[2 + t]);
          k[t] += 1.0;
        }
      }
      if (std::isfinite(e[2])) {
        if (!(e[2] >= lo)) lo = e[2];
        if (!(e[2] <= hi)) hi = e[2];
        gt += e[2] > 1.0 ? 1.0 : 0.0;
      }
    }
    o[0] = na;
    for (int t = 0; t < 4; ++t) o[1 + t] = k[t] > 0.0 ? s[t].val() / k[t] : nan;
    o[5] = lo;
    o[6] = hi;
    o[7] = hi - lo;
    o[8] = gt;
  }
  double wr = nan, mlo = nan, mhi = nan;
  int nm = 0;
  for (int y = 0; y < ncity; ++y) {
    const double *o = yout + 9 * static_cast<long long>(y);
    if (std::isfinite(o[7]) && !(o[7] <= wr)) wr = o[7];
    if (std::isfinite(o[1])) {
      ++nm;
      if (!(o[1] >= mlo)) mlo = o[1];
      if (!(o[1] <= mhi)) mhi = o[1];
    }
  }
  sum[0] = wr;
  sum[1] = nm >= 2 ? mhi - mlo : nan;
  sum[2] = b;
  return CLAB_OK;
}

int clab_lipaud(const double *h, int n, const int *pi, const int *pj,
                const double *d, int m, double lip, double *gap, double *rat,
                int *viol, double *sum) {
  if (n < 1 || m < 1) return CLAB_EARG;
  if (!kutil::nonneg(lip)) return CLAB_EDOM;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(h[i])) return CLAB_EDOM;
  }
  double lreq = 0.0, mg = -1.0, nv = 0.0;
  int kg = 0, kr = 0;
  for (int k = 0; k < m; ++k) {
    if (pi[k] < 0 || pi[k] >= n || pj[k] < 0 || pj[k] >= n) return CLAB_EARG;
    const double dk = d == nullptr ? 1.0 : d[k];
    if (!kutil::nonneg(dk)) return CLAB_EDOM;
    gap[k] = std::fabs(h[pi[k]] - h[pj[k]]);
    // a pair the metric cannot tell apart must be treated alike: no finite
    // constant allows a gap there
    rat[k] = dk > 0.0 ? gap[k] / dk
                      : (gap[k] > 0.0 ? std::numeric_limits<double>::infinity() : 0.0);
    viol[k] = gap[k] > lip * dk ? 1 : 0;
    nv += viol[k];
    if (gap[k] > mg) {
      mg = gap[k];
      kg = k;
    }
    if (rat[k] > lreq) {
      lreq = rat[k];
      kr = k;
    }
  }
  sum[0] = lreq;
  sum[1] = nv;
  sum[2] = nv / m;
  sum[3] = mg;
  sum[4] = kg;
  sum[5] = kr;
  return CLAB_OK;
}

int clab_nsbnd(const double *pobs, int n, double alo, double ahi, double blo,
               double bhi, double *lo, double *hi, int *ord, double *bd) {
  if (n < 1 || static_cast<double>(n) * (n - 1) / 2 > INT_MAX) return CLAB_EARG;
  if (!kutil::nonneg(alo) || !kutil::nonneg(blo) || !(ahi >= alo) ||
      !(bhi >= blo) || !(ahi + bhi < 1.0)) {
    return CLAB_EDOM;
  }
  for (int i = 0; i < n; ++i) {
    if (!(pobs[i] >= 0.0 && pobs[i] <= 1.0)) return CLAB_EDOM;
  }
  const double nan = kutil::nan();
  // p = (pobs - a) / (1 - a - b) falls in a and rises in b where p is in
  // [0, 1]; a rate no admissible noise reproduces (pobs < alo or
  // pobs > 1 - blo) has no bounds
  for (int i = 0; i < n; ++i) {
    const double q = pobs[i];
    if (q < alo || q > 1.0 - blo) {
      lo[i] = hi[i] = nan;
      continue;
    }
    const double a = std::fmin(ahi, q), bb = std::fmin(bhi, 1.0 - q);
    lo[i] = (q - a) / (1.0 - a - blo);
    hi[i] = (q - alo) / (1.0 - alo - bb);
  }
  int k = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j, ++k) {
      bd[k] = nan;
      if (std::isnan(lo[i]) || std::isnan(lo[j])) {
        ord[k] = 2;
        continue;
      }
      ord[k] = hi[i] < lo[j] ? 1 : (hi[j] < lo[i] ? -1 : 0);
      if (pobs[i] == pobs[j]) continue;
      // the largest upper end of the beta box (the alpha box held) at which
      // the lower recorded rate's upper bound stays below the other's lower
      // bound; capped at 1 - ahi, beyond which the box is not admissible
      const int l = pobs[i] < pobs[j] ? i : j, u = l == i ? j : i;
      if (!(lo[u] > 0.0)) continue;
      const double t = 1.0 - alo - (pobs[l] - alo) / lo[u];
      if (t > blo) bd[k] = std::fmin(t, 1.0 - ahi);
    }
  }
  return CLAB_OK;
}

int clab_impfpr(const double *p, const double *ppv, const double *fnr, int n,
                double *fpr) {
  if (n < 1) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (!(p[i] > 0.0 && p[i] < 1.0) || !(ppv[i] > 0.0 && ppv[i] <= 1.0) ||
        !(fnr[i] >= 0.0 && fnr[i] < 1.0)) {
      return CLAB_EDOM;
    }
    fpr[i] = p[i] / (1.0 - p[i]) * ((1.0 - ppv[i]) / ppv[i]) * (1.0 - fnr[i]);
  }
  return CLAB_OK;
}

int clab_rnkres(const double *est, const double *hw, int n, int *id,
                double *out) {
  if (n < 2 || static_cast<double>(n) * n > INT_MAX) return CLAB_EARG;
  double mh = 0.0;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(est[i]) || !kutil::nonneg(hw[i])) return CLAB_EDOM;
    mh = std::fmax(mh, hw[i]);
  }
  double nu = 0.0;
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      int &v = id[i + static_cast<long long>(n) * j];
      if (i == j) {
        v = -1;
        continue;
      }
      // the order of i and j survives every admissible noise exactly when
      // the estimates differ by more than the two half-widths together
      v = std::fabs(est[i] - est[j]) > hw[i] + hw[j] ? 1 : 0;
      if (i < j && v == 0) nu += 1.0;
    }
  }
  const double np = 0.5 * n * (n - 1.0);
  out[0] = np;
  out[1] = nu;
  out[2] = nu / np;
  out[3] = 2.0 * mh;
  return CLAB_OK;
}

}  // extern "C"
