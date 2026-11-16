// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, eyewitness lineups: the receiver operating characteristic of
// identifications by confidence, its partial area, the diagnosticity ratio
// and the confidence-accuracy characteristic (Mickes, Flowe and Wixted
// 2012; Mickes 2015), and the fairness of a lineup from mock witnesses'
// choices (Malpass 1981; Tredoux 1998).
#include <algorithm>
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

namespace {

// The counts of one procedure on resampled or original witnesses: at each
// cut (confidence at least lev[l]) the hit and false alarm rates.
void rates(const int *tp, const int *pk, const double *cf, const int *ix, int m, const double *lev, int nl,
           int size, double *hr, double *far) {
  double ntp = 0, nta = 0;
  for (int l = 0; l < nl; ++l) hr[l] = far[l] = 0.0;
  for (int q = 0; q < m; ++q) {
    const int i = ix[q];
    if (tp[i]) ntp += 1;
    else nta += 1;
    if (pk[i] == 0) continue;
    const bool hit = tp[i] && pk[i] == 1, fa = !tp[i] && (pk[i] == 1 || (size > 0 && pk[i] == 2));
    if (!hit && !fa) continue;
    // the cuts at or below the confidence: lev is decreasing
    const int l0 = static_cast<int>(std::lower_bound(lev, lev + nl, cf[i], [](double a, double b) { return a > b; }) - lev);
    for (int l = l0; l < nl; ++l) {
      if (hit) hr[l] += 1;
      else far[l] += 1;
    }
  }
  const double dta = size > 0 ? nta * size : nta;
  for (int l = 0; l < nl; ++l) {
    hr[l] = ntp > 0 ? hr[l] / ntp : kutil::nan();
    far[l] = dta > 0 ? far[l] / dta : kutil::nan();
  }
}

// The area under (0, 0) and the points (far, hr) up to false alarm rate
// fmax. Beyond the last point the curve runs straight to (xe, 1), the
// point where every lineup counts as an identification (xe = 1, or 1 /
// size), as pROC joins an ROC to (1, 1); flat at 1 beyond xe.
double parea(const double *hr, const double *far, int nl, double fmax, double xe) {
  double a = 0.0, x0 = 0.0, y0 = 0.0;
  for (int l = 0; l < nl && x0 < fmax; ++l) {
    const double x1 = far[l], y1 = hr[l];
    if (x1 > fmax) {
      a += (fmax - x0) * (y0 + (y0 + (y1 - y0) * (fmax - x0) / (x1 - x0))) / 2;
      return a;
    }
    a += (x1 - x0) * (y0 + y1) / 2;
    x0 = x1;
    y0 = y1;
  }
  if (x0 < fmax && x0 < xe) {
    const double x1 = std::fmin(fmax, xe), y1 = y0 + (1 - y0) * (x1 - x0) / (xe - x0);
    a += (x1 - x0) * (y0 + y1) / 2;
    x0 = x1;
    y0 = y1;
  }
  if (x0 < fmax) a += (fmax - x0) * y0;
  return a;
}

}  // namespace

int clab_lnroc(const int *tp, const int *pk, const double *cf, const int *pr, int n, int np, int size,
               double fmax, int nboot, unsigned int seed, double *lev, int *nlev, double *roc, double *cac,
               double *pauc, double *bse, double *fused) CLAB_TRY {
  if (n < 2 || np < 1 || size < 0) return CLAB_EARG;
  if (nboot < 0 || static_cast<double>(nboot) * n > 4e12) return CLAB_EARG;
  if (static_cast<double>(nboot) * np > kutil::maxcell) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (pr[i] < 0 || pr[i] >= np || pk[i] < 0 || pk[i] > 2 || (tp[i] != 0 && tp[i] != 1)) return CLAB_EARG;
    if (pk[i] != 0 && !std::isfinite(cf[i])) return CLAB_EDOM;
  }
  // the cuts: the distinct confidences of the identifications, decreasing
  int nl = 0;
  for (int i = 0; i < n; ++i) {
    if (pk[i] != 0) lev[nl++] = cf[i];
  }
  std::sort(lev, lev + nl, [](double a, double b) { return a > b; });
  nl = static_cast<int>(std::unique(lev, lev + nl) - lev);
  *nlev = nl;
  if (nl < 1) return CLAB_EDOM;
  const std::size_t L = static_cast<std::size_t>(nl);
  std::vector<std::vector<int>> ix(static_cast<std::size_t>(np));
  for (int i = 0; i < n; ++i) ix[static_cast<std::size_t>(pr[i])].push_back(i);
  // roc: np nl rows (procedure-major) x 3: hr, far, dr; cac: np nl x 3:
  // the hit and false alarm rates at the level alone, and the accuracy
  double fm = fmax;
  double lo = HUGE_VAL;
  for (int p = 0; p < np; ++p) {
    const std::vector<int> &w = ix[static_cast<std::size_t>(p)];
    double *hr = roc + p * L, *fr = roc + np * L + p * L;
    int tps = 0, tas = 0;
    for (int i : w) (tp[i] ? tps : tas)++;
    if (tps < 1 || tas < 1) return CLAB_EDOM;
    rates(tp, pk, cf, w.data(), static_cast<int>(w.size()), lev, nl, size, hr, fr);
    for (std::size_t l = 0; l < L; ++l) {
      roc[2 * np * L + p * L + l] = hr[l] / fr[l];
      const double h1 = hr[l] - (l > 0 ? hr[l - 1] : 0.0), f1 = fr[l] - (l > 0 ? fr[l - 1] : 0.0);
      cac[p * L + l] = h1;
      cac[np * L + p * L + l] = f1;
      cac[2 * np * L + p * L + l] = h1 + f1 > 0 ? h1 / (h1 + f1) : kutil::nan();
    }
    lo = std::fmin(lo, fr[L - 1]);
  }
  if (!(fm > 0)) fm = lo;
  const double xe = size > 0 ? 1.0 / size : 1.0;
  *fused = fm;
  if (!(fm > 0)) return CLAB_EDOM;
  for (int p = 0; p < np; ++p) pauc[p] = parea(roc + p * L, roc + np * L + p * L, nl, fm, xe);
  // the bootstrap of each procedure's partial area, within the
  // target-present and target-absent lineups of each procedure
  for (int p = 0; p < np; ++p) bse[p] = kutil::nan();
  if (nboot > 0) {
    kutil::rng rg(seed);
    std::vector<double> bt(static_cast<std::size_t>(nboot) * np), hr(L), fr(L);
    std::vector<int> s;
    for (int b = 0; b < nboot; ++b) {
      if (kutil::stop(static_cast<double>(n) * 4)) return CLAB_EINTR;
      for (int p = 0; p < np; ++p) {
        const std::vector<int> &w = ix[static_cast<std::size_t>(p)];
        std::vector<int> t1, t0;
        for (int i : w) (tp[i] ? t1 : t0).push_back(i);
        s.clear();
        for (std::size_t k = 0; k < t1.size(); ++k) s.push_back(t1[static_cast<std::size_t>(rg.below(static_cast<int>(t1.size())))]);
        for (std::size_t k = 0; k < t0.size(); ++k) s.push_back(t0[static_cast<std::size_t>(rg.below(static_cast<int>(t0.size())))]);
        rates(tp, pk, cf, s.data(), static_cast<int>(s.size()), lev, nl, size, hr.data(), fr.data());
        bt[static_cast<std::size_t>(b) * np + p] = parea(hr.data(), fr.data(), nl, fm, xe);
      }
    }
    // bse: np x 2: the bootstrap standard error of each area and of its
    // difference from the first procedure's
    for (int p = 0; p < np; ++p) {
      for (int k = 0; k < 2; ++k) {
        if (k == 1 && p == 0) {
          bse[np] = kutil::nan();
          continue;
        }
        double m = 0, ss = 0;
        for (int b = 0; b < nboot; ++b) {
          const double v = bt[static_cast<std::size_t>(b) * np + p] - (k ? bt[static_cast<std::size_t>(b) * np] : 0.0);
          m += v;
        }
        m /= nboot;
        for (int b = 0; b < nboot; ++b) {
          const double v = bt[static_cast<std::size_t>(b) * np + p] - (k ? bt[static_cast<std::size_t>(b) * np] : 0.0);
          ss += (v - m) * (v - m);
        }
        bse[k * np + p] = nboot > 1 ? std::sqrt(ss / (nboot - 1)) : kutil::nan();
      }
    }
  }
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_lnsize(const int *ch, int n, int k, int susp, double z, double *cnt, double *out) CLAB_TRY {
  if (n < 1 || k < 2 || susp < 0 || susp >= k || !(z > 0)) return CLAB_EARG;
  for (int j = 0; j < k; ++j) cnt[j] = 0.0;
  int none = 0;
  for (int i = 0; i < n; ++i) {
    if (ch[i] == -1) {
      ++none;
      continue;
    }
    if (ch[i] < 0 || ch[i] >= k) return CLAB_EDOM;
    cnt[ch[i]] += 1;
  }
  // N all mock witnesses, M those who chose a member
  const double N = n, M = n - none;
  if (!(M > 0)) return CLAB_EDOM;
  const double e = M / k;
  double s2 = 0, s3 = 0, dev = 0;
  int ka = 0;
  for (int j = 0; j < k; ++j) {
    const double p = cnt[j] / M;
    s2 += p * p;
    s3 += p * p * p;
    dev += std::fabs(cnt[j] - e);
    ka += cnt[j] > 0;
  }
  // Malpass's adjustment: the null foils (members nobody chose) leave the
  // nominal size and the chance expectation
  const double ea = M / ka;
  double deva = 0;
  for (int j = 0; j < k; ++j) {
    if (cnt[j] > 0) deva += std::fabs(cnt[j] - ea);
  }
  const double ns = cnt[susp];
  // out: n, k, m, the suspect's share, the binomial P(X >= ns) of the m
  // choices under 1 / k, functional size, Malpass's effective size with
  // and without the adjustment, Tredoux's E' and its interval
  out[0] = N;
  out[1] = k;
  out[2] = M;
  out[3] = ns / N;
  out[4] = ns > 0 ? clab_pbinom(ns - 1, M, 1.0 / k, 0) : 1.0;
  out[5] = ns > 0 ? N / ns : HUGE_VAL;
  out[6] = ka - deva / (2 * ea);
  out[7] = k - dev / (2 * e);
  out[8] = 1 / s2;
  // the interval of I = 1 - sum p^2 (variance 4 (sum p^3 - (sum p^2)^2) /
  // m, Agresti and Agresti 1978) with its ends taken to E' = 1 / (1 - I)
  const double h = z * std::sqrt(std::fmax(4 * (s3 - s2 * s2) / M, 0.0));
  out[9] = 1 / (s2 + h);
  out[10] = s2 - h > 0 ? 1 / (s2 - h) : HUGE_VAL;
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)
