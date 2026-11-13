// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, courts and sentencing: the judge-leniency instrument on
// observed data (clab_jdgiv) and on a known population (clab_jdgpop),
// monotonicity checks across judges (clab_jrate, clab_jslope, clab_jsub)
// and the shares sentenced at, below and above a mandatory minimum
// (clab_mmshr).
#include <cmath>
#include <initializer_list>

#include "crimlab.h"
#include "kutil.h"

namespace {

bool bin(double v) { return v == 0 || v == 1; }

}  // namespace

int clab_jdgiv(const double *z, const double *d, const double *y, const double *w, int n,
               const double *dsh, int nd, const double *deff, int ne, double *est, double *sens) {
  if (n < 2 || nd < 0 || ne < 0) return CLAB_EARG;
  kutil::ksum w1, w0, d1, d0, y1, y0;
  for (int i = 0; i < n; ++i) {
    const double wi = w ? w[i] : 1.0;
    if (!bin(z[i]) || !bin(d[i]) || !std::isfinite(y[i]) || !kutil::nonneg(wi)) return CLAB_EDOM;
    if (z[i] == 1) {
      w1.add(wi);
      d1.add(wi * d[i]);
      y1.add(wi * y[i]);
    } else {
      w0.add(wi);
      d0.add(wi * d[i]);
      y0.add(wi * y[i]);
    }
  }
  if (!(w1.val() > 0 && w0.val() > 0)) return CLAB_EDOM;
  for (int k = 0; k < nd; ++k) {
    if (!(dsh[k] >= 0 && dsh[k] <= 1)) return CLAB_EDOM;
  }
  for (int k = 0; k < ne; ++k) {
    if (!std::isfinite(deff[k])) return CLAB_EDOM;
  }
  const double fs = d1.val() / w1.val() - d0.val() / w0.val();
  const double itt = y1.val() / w1.val() - y0.val() / w0.val();
  const double wald = fs != 0 ? itt / fs : kutil::nan();
  const double al = d0.val() / w0.val(), nv = 1 - d1.val() / w1.val();
  est[0] = fs;
  est[1] = itt;
  est[2] = wald;
  est[3] = 1 - al - nv;
  est[4] = al;
  est[5] = nv;
  // with a defier share s and defier effect e, the compliers number fs +
  // s and their effect solves wald = (pc tc - s e) / fs
  for (int a = 0; a < nd; ++a) {
    for (int b = 0; b < ne; ++b) {
      const int k = a + nd * b;
      const double cs = est[3] + dsh[a];
      sens[k] = cs > 0 && !std::isnan(wald) ? (wald * (cs - dsh[a]) + dsh[a] * deff[b]) / cs
                                            : kutil::nan();
    }
  }
  return CLAB_OK;
}

int clab_jdgpop(const double *d0, const double *d1, const double *y0, const double *y1,
                const double *w, int n, double *out) {
  if (n < 1) return CLAB_EARG;
  kutil::ksum tw;
  for (int i = 0; i < n; ++i) {
    const double wi = w ? w[i] : 1.0;
    if (!bin(d0[i]) || !bin(d1[i]) || !std::isfinite(y0[i]) || !std::isfinite(y1[i]) ||
        !kutil::nonneg(wi)) {
      return CLAB_EDOM;
    }
    tw.add(wi);
  }
  if (!(tw.val() > 0)) return CLAB_EDOM;
  // types: 0 complier (0, 1), 1 defier (1, 0), 2 always, 3 never
  kutil::ksum sh[4], ef[4], itt, fs, ate;
  for (int i = 0; i < n; ++i) {
    const double wi = (w ? w[i] : 1.0) / tw.val();
    const int t = d0[i] == 0 && d1[i] == 1 ? 0 : d0[i] == 1 && d1[i] == 0 ? 1 : d0[i] == 1 ? 2 : 3;
    sh[t].add(wi);
    ef[t].add(wi * (y1[i] - y0[i]));
    itt.add(wi * ((d1[i] == 1 ? y1[i] : y0[i]) - (d0[i] == 1 ? y1[i] : y0[i])));
    fs.add(wi * (d1[i] - d0[i]));
    ate.add(wi * (y1[i] - y0[i]));
  }
  for (int t = 0; t < 4; ++t) {
    out[t] = sh[t].val();
    out[4 + t] = sh[t].val() > 0 ? ef[t].val() / sh[t].val() : kutil::nan();
  }
  out[8] = itt.val();
  out[9] = fs.val();
  out[10] = fs.val() != 0 ? itt.val() / fs.val() : kutil::nan();
  out[11] = ate.val();
  return CLAB_OK;
}

int clab_jrate(const int *jc, const double *d, const double *y, const double *w, int n, int nj,
               double *cnt, double *one, double *rate, double *outc, double *st) {
  if (n < 1 || nj < 1) return CLAB_EARG;
  for (int j = 0; j < nj; ++j) cnt[j] = one[j] = rate[j] = outc[j] = 0.0;
  // rate and outc first hold weighted sums; st[4] the weights' total
  double *ws = st + 4;
  for (int j = 0; j < nj; ++j) ws[j] = 0.0;
  for (int i = 0; i < n; ++i) {
    const double wi = w ? w[i] : 1.0;
    if (jc[i] < 0 || jc[i] >= nj) return CLAB_EARG;
    if (!bin(d[i]) || !kutil::nonneg(wi) || (y && !std::isfinite(y[i]))) return CLAB_EDOM;
    cnt[jc[i]] += 1;
    one[jc[i]] += d[i];
    ws[jc[i]] += wi;
    rate[jc[i]] += wi * d[i];
    if (y) outc[jc[i]] += wi * y[i];
  }
  kutil::ksum t1;
  double mn = INFINITY, mx = -INFINITY;
  for (int j = 0; j < nj; ++j) {
    if (!(ws[j] > 0)) return CLAB_EDOM;
    rate[j] /= ws[j];
    outc[j] = y ? outc[j] / ws[j] : kutil::nan();
    mn = std::fmin(mn, rate[j]);
    mx = std::fmax(mx, rate[j]);
    t1.add(one[j]);
  }
  // Pearson's chi-square of the decision against the judge, on the case
  // counts: expected cnt p ones and cnt (1 - p) zeros at the pooled p
  const double p = t1.val() / n;
  kutil::ksum x2;
  for (int j = 0; j < nj; ++j) {
    const double e1 = cnt[j] * p, e0 = cnt[j] * (1 - p);
    x2.add((one[j] - e1) * (one[j] - e1) / e1);
    x2.add((cnt[j] - one[j] - e0) * (cnt[j] - one[j] - e0) / e0);
  }
  st[0] = p > 0 && p < 1 && nj > 1 ? x2.val() : kutil::nan();
  st[1] = nj - 1;
  st[2] = mn;
  st[3] = mx;
  return CLAB_OK;
}

int clab_jslope(const double *rate, const double *outc, int nj, double lo, double hi, double *dp,
                double *dy, double *bnd, int *viol, double *late) {
  if (nj < 1) return CLAB_EARG;
  if (!(hi >= lo) || !std::isfinite(lo) || !std::isfinite(hi)) return CLAB_EDOM;
  long long k = 0;
  // every pair, the more lenient (lower rate) first: nested judges move
  // the outcome by at most (hi - lo) times the difference in rates
  for (int a = 0; a < nj; ++a) {
    for (int b = a + 1; b < nj; ++b, ++k) {
      dp[k] = rate[b] - rate[a];
      dy[k] = outc[b] - outc[a];
      bnd[k] = (hi - lo) * dp[k];
      viol[k] = std::fabs(dy[k]) > bnd[k] + 1e-12;
      late[k] = dp[k] > 0 ? dy[k] / dp[k] : kutil::nan();
    }
  }
  return CLAB_OK;
}

int clab_jsub(const int *jc, const int *gc, const double *d, int n, int nj, int ng, double minn,
              const double *jt, double *cn, double *cr, int *pg, int *ps, int *pl, double *pv,
              double *sl, int *sj, double *sc, int *npair) {
  if (n < 1 || nj < 1 || ng < 1) return CLAB_EARG;
  const long long nc = static_cast<long long>(nj) * ng;
  for (long long c = 0; c < nc; ++c) cn[c] = cr[c] = 0.0;
  for (int i = 0; i < n; ++i) {
    if (jc[i] < 0 || jc[i] >= nj || gc[i] < 0 || gc[i] >= ng) return CLAB_EARG;
    if (!bin(d[i])) return CLAB_EDOM;
    const long long c = jc[i] + static_cast<long long>(nj) * gc[i];
    cn[c] += 1;
    cr[c] += d[i];
  }
  for (long long c = 0; c < nc; ++c) cr[c] = cn[c] > 0 ? cr[c] / cn[c] : 0.0;
  long long k = 0;
  for (int g = 0; g < ng; ++g) {
    const double *cg = cn + static_cast<long long>(nj) * g,
                 *rg = cr + static_cast<long long>(nj) * g;
    // the slope of the group's rates on the overall ones, weighted by the
    // group's cases, over the judges with at least minn of them
    kutil::ksum sw, sx, sy;
    int m = 0;
    for (int j = 0; j < nj; ++j) {
      if (cg[j] < minn || cg[j] == 0) continue;
      ++m;
      sw.add(cg[j]);
      sx.add(cg[j] * jt[j]);
      sy.add(cg[j] * rg[j]);
    }
    sj[g] = m;
    sc[g] = sw.val();
    sl[g] = kutil::nan();
    if (m >= 2) {
      const double mx = sx.val() / sw.val(), my = sy.val() / sw.val();
      kutil::ksum sxy, sxx;
      for (int j = 0; j < nj; ++j) {
        if (cg[j] < minn || cg[j] == 0) continue;
        sxy.add(cg[j] * (jt[j] - mx) * (rg[j] - my));
        sxx.add(cg[j] * (jt[j] - mx) * (jt[j] - mx));
      }
      if (sxx.val() > 0) sl[g] = sxy.val() / sxx.val();
    }
    for (int a = 0; a < nj; ++a) {
      if (cg[a] < minn || cg[a] == 0) continue;
      for (int b = a + 1; b < nj; ++b) {
        if (cg[b] < minn || cg[b] == 0 || jt[a] == jt[b]) continue;
        // the stricter judge overall (the higher rate) first
        const int h = jt[a] < jt[b] ? b : a, l = h == a ? b : a;
        if (pg) {
          // the two-proportion z of the lenient judge's rate over the
          // stricter's in this group, pooled
          const double nh = cg[h], nl = cg[l], ph = rg[h], pq = rg[l];
          const double pp = (ph * nh + pq * nl) / (nh + nl);
          const double se = std::sqrt(pp * (1 - pp) * (1 / nh + 1 / nl));
          pg[k] = g;
          ps[k] = h;
          pl[k] = l;
          pv[k] = se > 0 ? (pq - ph) / se : 0.0;
        }
        ++k;
      }
    }
  }
  *npair = static_cast<int>(k);
  return CLAB_OK;
}

int clab_mmshr(const double *sent, const double *mn, const int *gc, int n, int ng, double *out) {
  if (n < 1 || ng < 0) return CLAB_EARG;
  // per group, then overall (row ng): n, at, below, above, the mean excess
  // of those above
  double *cnt = out, *at = cnt + ng + 1, *be = at + ng + 1, *ab = be + ng + 1, *mx = ab + ng + 1;
  for (int g = 0; g <= ng; ++g) cnt[g] = at[g] = be[g] = ab[g] = mx[g] = 0.0;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(sent[i]) || !std::isfinite(mn[i])) return CLAB_EDOM;
    const int g = gc ? gc[i] : 0;
    if (gc && (g < 0 || g >= ng)) return CLAB_EARG;
    const double e = sent[i] - mn[i];
    for (int r : {gc ? g : ng, ng}) {
      cnt[r] += 1;
      at[r] += e == 0;
      be[r] += e < 0;
      ab[r] += e > 0;
      if (e > 0) mx[r] += e;
      if (!gc) break;
    }
  }
  for (int g = 0; g <= ng; ++g) {
    mx[g] = ab[g] > 0 ? mx[g] / ab[g] : 0.0;
    if (cnt[g] > 0) {
      at[g] /= cnt[g];
      be[g] /= cnt[g];
      ab[g] /= cnt[g];
    }
  }
  return CLAB_OK;
}
