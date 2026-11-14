// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, custody and courts as stocks and flows: person-days,
// average daily population and length of stay (clab_stkflw), the
// distribution of stays (clab_lostay), Little's law on a docket
// (clab_bklog), case flow and the mix of outcomes over periods
// (clab_dktflw, clab_dkmix), the actuarial life table (clab_lftab) and
// the transition matrix between states (clab_trnmat).
#include <cmath>

#include "crimlab.h"
#include "kutil.h"

namespace {

// the change of v[i] against v[0] (base 0) or v[i - 1] (base 1), in percent
double chg(const double *v, int i, int base) {
  const double b = base == 0 ? v[0] : i > 0 ? v[i - 1] : kutil::nan();
  return b != 0 ? 100 * (v[i] / b - 1) : kutil::nan();
}

}  // namespace

int clab_stkflw(const double *pe, const double *dy, const double *ad, const double *al,
                const double *t, int n, int base, const double *expo, double per, const double *rel,
                double early, double *out) {
  if (n < 1 || base < 0 || base > 1) return CLAB_EARG;
  const int given = (pe != nullptr) + (dy != nullptr) + (ad != nullptr) + (al != nullptr);
  if (given != 2 || (dy && ad)) return CLAB_EARG;
  if (!kutil::pos(per) || !(std::isnan(early) || kutil::nonneg(early))) return CLAB_EDOM;
  double *p = out, *d = p + n, *a = d + n, *l = a + n;
  for (int i = 0; i < n; ++i) {
    if (!kutil::pos(t[i]) || (pe && !kutil::pos(pe[i])) || (dy && !kutil::nonneg(dy[i])) ||
        (ad && !kutil::nonneg(ad[i])) || (al && !kutil::pos(al[i])) ||
        (expo && !kutil::pos(expo[i])) || (rel && !kutil::nonneg(rel[i]))) {
      return CLAB_EDOM;
    }
    // days = adp t = people alos: any two of the four but days and adp
    d[i] = dy ? dy[i] : ad ? ad[i] * t[i] : pe[i] * al[i];
    p[i] = pe ? pe[i] : d[i] / al[i];
    a[i] = d[i] / t[i];
    l[i] = al ? al[i] : d[i] / p[i];
  }
  // the changes against the base: people, alos, days, adp; then the share
  // of the change in days (on the log scale) that people account for
  double *c = l + n, *sh = c + 4 * n;
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < 4; ++k) {
      const double *v = k == 0 ? p : k == 1 ? l : k == 2 ? d : a;
      c[k * n + i] = chg(v, i, base);
    }
    const double ld = std::log1p(c[2 * n + i] / 100);
    sh[i] = ld != 0 && std::isfinite(ld) ? std::log1p(c[i] / 100) / ld : kutil::nan();
  }
  // rates per per of the exposure, flow (people) and stock (adp), and their
  // changes; the stock an earlier release by early days removes
  double *fr = sh + n, *sr = fr + n, *fc = sr + n, *sc = fc + n, *cut = sc + n, *cp = cut + n;
  for (int i = 0; i < n; ++i) {
    fr[i] = expo ? per * p[i] / expo[i] : kutil::nan();
    sr[i] = expo ? per * a[i] / expo[i] : kutil::nan();
  }
  for (int i = 0; i < n; ++i) {
    fc[i] = expo ? chg(fr, i, base) : kutil::nan();
    sc[i] = expo ? chg(sr, i, base) : kutil::nan();
    cut[i] = rel && !std::isnan(early) ? rel[i] * early / t[i] : kutil::nan();
    cp[i] = a[i] > 0 ? 100 * cut[i] / a[i] : kutil::nan();
  }
  // the means over the periods
  double *mn = cp + n;
  for (int k = 0; k < 4; ++k) {
    kutil::ksum s;
    for (int i = 0; i < n; ++i) s.add(out[k * n + i]);
    mn[k] = s.val() / n;
  }
  return CLAB_OK;
}

int clab_lostay(const double *x, const double *w, int n, const double *br, int nb,
                const double *thr, int nt, const int *flag, double *st, double *bt, double *tt,
                double *fit, int *o) {
  if (n < 1 || nb < 0 || nt < 0) return CLAB_EARG;
  kutil::ksum sw, sx;
  for (int i = 0; i < n; ++i) {
    const double wi = w ? w[i] : 1.0;
    if (!kutil::nonneg(x[i]) || !kutil::nonneg(wi) || wi != std::floor(wi)) return CLAB_EDOM;
    sw.add(wi);
    sx.add(wi * x[i]);
    o[i] = i;
  }
  for (int k = 1; k < nb; ++k) {
    if (!(br[k] > br[k - 1])) return CLAB_EDOM;
  }
  const double nn = sw.val();
  if (!(nn > 0)) return CLAB_EDOM;
  kutil::hsort(o, n, [&](int a, int b) { return x[a] == x[b] ? a < b : x[a] < x[b]; });
  const double m = sx.val() / nn;
  kutil::ksum ss;
  for (int i = 0; i < n; ++i) ss.add((w ? w[i] : 1.0) * (x[i] - m) * (x[i] - m));
  const double sd = nn > 1 ? std::sqrt(ss.val() / (nn - 1)) : kutil::nan();
  // the type-7 quantile of the expanded values: the order statistics
  // floor(h) and ceil(h), h = (N - 1) p, from the cumulative weights
  auto ord = [&](double k) {
    double c = 0;
    for (int i = 0; i + 1 < n; ++i) {
      c += w ? w[o[i]] : 1.0;
      if (c > k) return x[o[i]];
    }
    return x[o[n - 1]];
  };
  auto q = [&](double p) {
    const double h = (nn - 1) * p, lo = std::floor(h);
    return ord(lo) + (h - lo) * (ord(std::ceil(h)) - ord(lo));
  };
  st[0] = nn;
  st[1] = sx.val();
  st[2] = m;
  st[3] = sd;
  st[4] = q(0.5);
  st[5] = q(0.25);
  st[6] = q(0.75);
  st[7] = x[o[n - 1]];
  st[8] = sd / std::sqrt(nn);
  // bands (br[k], br[k + 1]]: counts, shares and the flagged share
  for (int k = 0; k + 1 < nb; ++k) bt[k] = bt[nb - 1 + k] = 0.0;
  for (int i = 0; i < n; ++i) {
    const double wi = w ? w[i] : 1.0;
    int a = 0, b = nb - 1;
    if (nb < 2 || !(x[i] > br[0]) || x[i] > br[nb - 1]) continue;
    // the band k with br[k] < x <= br[k + 1]
    while (b - a > 1) {
      const int c = (a + b) / 2;
      if (x[i] > br[c])
        a = c;
      else
        b = c;
    }
    bt[a] += wi;
    if (flag && flag[i]) bt[nb - 1 + a] += wi;
  }
  for (int i = 0; i < nt; ++i) {
    kutil::ksum s;
    for (int j = 0; j < n; ++j) {
      if (x[j] <= thr[i]) s.add(w ? w[j] : 1.0);
    }
    tt[i] = s.val() / nn;
  }
  // the fits of the positive stays by maximum likelihood: lognormal,
  // exponential, Pareto from the smallest; AIC = 2 k - 2 log-likelihood
  kutil::ksum pw, pl, pl2, px;
  double xm = INFINITY;
  for (int i = 0; i < n; ++i) {
    if (!(x[i] > 0)) continue;
    const double wi = w ? w[i] : 1.0;
    pw.add(wi);
    pl.add(wi * std::log(x[i]));
    px.add(wi * x[i]);
    if (wi > 0) xm = std::fmin(xm, x[i]);
  }
  const double np = pw.val();
  for (int k = 0; k < 6; ++k) fit[k] = kutil::nan();
  if (np >= 2) {
    const double mu = pl.val() / np;
    for (int i = 0; i < n; ++i) {
      if (x[i] > 0) pl2.add((w ? w[i] : 1.0) * (std::log(x[i]) - mu) * (std::log(x[i]) - mu));
    }
    const double s2 = pl2.val() / np;
    const double pi = 3.14159265358979323846;
    // lognormal: -sum log x - n/2 log(2 pi s2) - n/2
    fit[0] = s2 > 0 ? 4 - 2 * (-pl.val() - np / 2 * std::log(2 * pi * s2) - np / 2) : kutil::nan();
    const double rate = np / px.val();
    fit[1] = 2 - 2 * (np * std::log(rate) - rate * px.val());
    const double lr = pl.val() - np * std::log(xm);
    if (lr > 0) {
      const double al = np / lr;
      fit[2] = 4 - 2 * (np * std::log(al) + np * al * std::log(xm) - (al + 1) * pl.val());
      fit[5] = al;
    }
    fit[3] = mu;
    fit[4] = std::sqrt(s2);
  }
  return CLAB_OK;
}

int clab_bklog(const double *a, const double *d, int n, double h, double *out) {
  if (n < 1) return CLAB_EARG;
  if (!kutil::pos(h)) return CLAB_EDOM;
  kutil::ksum occ, wt;
  int nc = 0, nd = 0;
  for (int i = 0; i < n; ++i) {
    if (!kutil::nonneg(a[i]) || a[i] > h) return CLAB_EDOM;
    // pending (no disposition, or one after the horizon) until the horizon
    const bool cen = std::isnan(d[i]) || d[i] > h;
    if (!cen && d[i] < a[i]) return CLAB_EDOM;
    if (cen) {
      ++nc;
      occ.add(h - a[i]);
    } else {
      ++nd;
      occ.add(d[i] - a[i]);
      wt.add(d[i] - a[i]);
    }
  }
  const double rate = n / h, w = nd ? wt.val() / nd : kutil::nan();
  out[0] = n;
  out[1] = nc;
  out[2] = rate;
  out[3] = w;
  out[4] = occ.val() / h;
  out[5] = occ.val();
  // L = lambda W holds exactly with every case disposed
  out[6] = nc == 0 ? occ.val() / h - rate * w : kutil::nan();
  return CLAB_OK;
}

int clab_dktflw(const double *rc, const double *ds, const double *pd, int n, int win, double *out) {
  if (n < 1 || win < 1) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (!kutil::nonneg(rc[i]) || !kutil::nonneg(ds[i]) || (pd && !kutil::nonneg(pd[i]))) {
      return CLAB_EDOM;
    }
  }
  // clearance, net, the rolling sums over win periods, the pending change,
  // the periods to clear the pending at the recent pace, the cumulative net
  double *cl = out, *nt = cl + n, *rr = nt + n, *rd = rr + n, *rcl = rd + n, *pc = rcl + n,
         *tc = pc + n, *cn = tc + n;
  kutil::ksum sr, sd, sn;
  for (int i = 0; i < n; ++i) {
    cl[i] = rc[i] > 0 ? ds[i] / rc[i] : kutil::nan();
    nt[i] = rc[i] - ds[i];
    sr.add(rc[i]);
    sd.add(ds[i]);
    if (i >= win) {
      sr.add(-rc[i - win]);
      sd.add(-ds[i - win]);
    }
    rr[i] = i >= win - 1 ? sr.val() : kutil::nan();
    rd[i] = i >= win - 1 ? sd.val() : kutil::nan();
    rcl[i] = rr[i] > 0 ? rd[i] / rr[i] : kutil::nan();
    sn.add(nt[i]);
    cn[i] = sn.val();
  }
  kutil::ksum sp;
  for (int i = 0; i < n; ++i) {
    pc[i] = pd && i > 0 ? pd[i] - pd[i - 1] : kutil::nan();
    // the mean disposals over the window ending at i (shorter at the start)
    sp.add(ds[i]);
    if (i >= win) sp.add(-ds[i - win]);
    const double pace = sp.val() / (i + 1 < win ? i + 1 : win);
    tc[i] = pd && pace > 0 ? pd[i] / pace : kutil::nan();
  }
  return CLAB_OK;
}

int clab_dkmix(const double *x, int nt, int k, const double *sc, double a, double *sh, double *lo,
               double *hi, double *st, double *tr) {
  if (nt < 1 || k < 1) return CLAB_EARG;
  if (!(a > 0)) return CLAB_EDOM;
  // x: nt periods by k outcomes, column-major; Goodman's simultaneous
  // intervals with a, the chi-square quantile at 1 - (1 - level) / k
  for (int t = 0; t < nt; ++t) {
    kutil::ksum s;
    for (int j = 0; j < k; ++j) {
      const double v = x[t + static_cast<long long>(nt) * j];
      if (!kutil::nonneg(v)) return CLAB_EDOM;
      s.add(v);
    }
    const double tot = s.val();
    if (!(tot > 0)) return CLAB_EDOM;
    for (int j = 0; j < k; ++j) {
      const long long c = t + static_cast<long long>(nt) * j;
      const double ni = x[c], hf = std::sqrt(a * (a + 4 * ni * (tot - ni) / tot));
      sh[c] = ni / tot;
      lo[c] = (a + 2 * ni - hf) / (2 * (tot + a));
      hi[c] = (a + 2 * ni + hf) / (2 * (tot + a));
    }
  }
  // the homogeneity chi-square over the outcomes that occur, Cramer's V
  kutil::ksum all;
  for (long long c = 0; c < static_cast<long long>(nt) * k; ++c) all.add(x[c]);
  const double tn = all.val();
  kutil::ksum x2;
  int kk = 0;
  for (int j = 0; j < k; ++j) {
    kutil::ksum s;
    for (int t = 0; t < nt; ++t) s.add(x[t + static_cast<long long>(nt) * j]);
    const double cs = s.val();
    if (!(cs > 0)) continue;
    ++kk;
    for (int t = 0; t < nt; ++t) {
      kutil::ksum r;
      for (int q = 0; q < k; ++q) r.add(x[t + static_cast<long long>(nt) * q]);
      const double e = r.val() * cs / tn, o = x[t + static_cast<long long>(nt) * j];
      x2.add((o - e) * (o - e) / e);
    }
  }
  st[0] = x2.val();
  st[1] = (nt - 1.0) * (kk - 1.0);
  st[2] =
      nt > 1 && kk > 1 ? std::sqrt(x2.val() / (tn * std::fmin(nt - 1.0, kk - 1.0))) : kutil::nan();
  // the Cochran-Armitage trend statistic of each outcome's share over the
  // scores (as prop.trend.test), and the change from first to last period
  for (int j = 0; j < k; ++j) {
    kutil::ksum sx, sn, sxs, ss, sss;
    for (int t = 0; t < nt; ++t) {
      kutil::ksum r;
      for (int q = 0; q < k; ++q) r.add(x[t + static_cast<long long>(nt) * q]);
      const double v = x[t + static_cast<long long>(nt) * j];
      sx.add(v);
      sn.add(r.val());
      sxs.add(v * sc[t]);
      ss.add(r.val() * sc[t]);
      sss.add(r.val() * sc[t] * sc[t]);
    }
    const double p = sx.val() / sn.val(), sb = ss.val() / sn.val();
    const double num = sxs.val() - p * ss.val();
    const double den = p * (1 - p) * (sss.val() - sn.val() * sb * sb);
    tr[j] = nt >= 3 && den > 0 ? num * num / den : kutil::nan();
    tr[k + j] = sh[nt - 1 + static_cast<long long>(nt) * j] - sh[static_cast<long long>(nt) * j];
  }
  return CLAB_OK;
}

int clab_lftab(const double *lo, const double *hi, const double *d, const double *w, int k,
               double z, int ci, double *out) {
  if (k < 1 || ci < 0 || ci > 2) return CLAB_EARG;
  kutil::ksum tot;
  for (int i = 0; i < k; ++i) {
    if (!std::isfinite(lo[i]) || !(hi[i] > lo[i]) || (i > 0 && lo[i] != hi[i - 1]))
      return CLAB_EDOM;
    if (i + 1 < k && !std::isfinite(hi[i])) return CLAB_EDOM;
    if (!kutil::nonneg(d[i]) || !kutil::nonneg(w[i])) return CLAB_EDOM;
    tot.add(d[i] + w[i]);
  }
  if (!(tot.val() > 0)) return CLAB_EDOM;
  double *en = out, *ar = en + k, *q = ar + k, *s = q + k, *se = s + k, *cl = se + k, *ch = cl + k,
         *hz = ch + k, *med = hz + k;
  double run = tot.val(), sv = 1, gw = 0;
  for (int i = 0; i < k; ++i) {
    // a censored case counts as half a case at risk in its band
    en[i] = run;
    ar[i] = run - w[i] / 2;
    q[i] = ar[i] > 0 ? std::fmin(d[i] / ar[i], 1.0) : kutil::nan();
    if (q[i] > 1 - 1e-9) q[i] = 1;
    sv *= std::isnan(q[i]) ? 1 : 1 - q[i];
    s[i] = sv;
    gw += ar[i] > 0 && q[i] < 1 ? d[i] / (ar[i] * (ar[i] - d[i])) : kutil::nan();
    se[i] = sv * std::sqrt(gw);
    const double g = std::sqrt(gw);
    double a, b;
    if (ci == 0) {
      a = std::fmax(sv - z * se[i], 0.0);
      b = std::fmin(sv + z * se[i], 1.0);
    } else if (ci == 1) {
      a = sv * std::exp(-z * g);
      b = std::fmin(sv * std::exp(z * g), 1.0);
    } else {
      const double u = g / std::log(sv);
      a = std::pow(sv, std::exp(-z * u));
      b = std::pow(sv, std::exp(z * u));
    }
    const bool ok = sv > 0 && sv < 1 && std::isfinite(gw);
    cl[i] = ok ? a : kutil::nan();
    ch[i] = ok ? b : kutil::nan();
    const double wd = hi[i] - lo[i];
    hz[i] =
        std::isfinite(wd) && ar[i] - d[i] / 2 > 0 ? d[i] / (wd * (ar[i] - d[i] / 2)) : kutil::nan();
    run -= d[i] + w[i];
  }
  // the median, linear within the band where survival crosses 1/2
  *med = INFINITY;
  for (int i = 0; i < k; ++i) {
    if (s[i] <= 0.5) {
      const double s0 = i > 0 ? s[i - 1] : 1.0;
      *med = std::isfinite(hi[i]) ? lo[i] + (hi[i] - lo[i]) * (s0 - 0.5) / (s0 - s[i]) : INFINITY;
      break;
    }
  }
  return CLAB_OK;
}

int clab_trnmat(const int *from, const int *to, int m, int k, double *cnt, double *p, double *pi,
                double *st, double *work) {
  if (m < 0 || k < 1) return CLAB_EARG;
  const long long kk = static_cast<long long>(k) * k;
  for (long long c = 0; c < kk; ++c) cnt[c] = p[c] = 0.0;
  for (int i = 0; i < m; ++i) {
    if (from[i] < 0 || from[i] >= k || to[i] < 0 || to[i] >= k) return CLAB_EARG;
    cnt[from[i] + static_cast<long long>(k) * to[i]] += 1;
  }
  bool exits = true;
  for (int a = 0; a < k; ++a) {
    kutil::ksum r;
    for (int b = 0; b < k; ++b) r.add(cnt[a + static_cast<long long>(k) * b]);
    if (r.val() > 0) {
      for (int b = 0; b < k; ++b) {
        p[a + static_cast<long long>(k) * b] = cnt[a + static_cast<long long>(k) * b] / r.val();
      }
    } else {
      exits = false;
    }
  }
  // the stationary distribution: pi (P - I) = 0 with sum pi = 1, by
  // Gaussian elimination with partial pivoting on (P - I)' with its last
  // equation replaced by the sum; NaN when a state has no exit or the
  // system is singular (no unique distribution)
  for (int a = 0; a < k; ++a) pi[a] = kutil::nan();
  if (exits && m > 0) {
    double *A = work, *b = work + kk;
    for (int r = 0; r < k; ++r) {
      for (int c = 0; c < k; ++c) {
        A[r + static_cast<long long>(k) * c] =
            r == k - 1 ? 1.0 : p[c + static_cast<long long>(k) * r] - (r == c ? 1.0 : 0.0);
      }
      b[r] = r == k - 1 ? 1.0 : 0.0;
    }
    bool sing = false;
    for (int c = 0; c < k && !sing; ++c) {
      int pv = c;
      for (int r = c + 1; r < k; ++r) {
        if (std::fabs(A[r + static_cast<long long>(k) * c]) >
            std::fabs(A[pv + static_cast<long long>(k) * c])) {
          pv = r;
        }
      }
      if (std::fabs(A[pv + static_cast<long long>(k) * c]) < 1e-12) {
        sing = true;
        break;
      }
      for (int q = 0; q < k; ++q) {
        const double t = A[c + static_cast<long long>(k) * q];
        A[c + static_cast<long long>(k) * q] = A[pv + static_cast<long long>(k) * q];
        A[pv + static_cast<long long>(k) * q] = t;
      }
      const double tb = b[c];
      b[c] = b[pv];
      b[pv] = tb;
      for (int r = c + 1; r < k; ++r) {
        const double f =
            A[r + static_cast<long long>(k) * c] / A[c + static_cast<long long>(k) * c];
        for (int q = c; q < k; ++q) {
          A[r + static_cast<long long>(k) * q] -= f * A[c + static_cast<long long>(k) * q];
        }
        b[r] -= f * b[c];
      }
    }
    if (!sing) {
      for (int r = k - 1; r >= 0; --r) {
        double s = b[r];
        for (int q = r + 1; q < k; ++q) s -= A[r + static_cast<long long>(k) * q] * pi[q];
        pi[r] = s / A[r + static_cast<long long>(k) * r];
      }
    }
  }
  // the share of transitions that stay in their state, and the Theil
  // index of the moves between states: sum p log(p q) over the q off-
  // diagonal cells with any (p their shares)
  kutil::ksum all, dg, off;
  for (int a = 0; a < k; ++a) {
    for (int b = 0; b < k; ++b) {
      const double v = cnt[a + static_cast<long long>(k) * b];
      all.add(v);
      if (a == b)
        dg.add(v);
      else
        off.add(v);
    }
  }
  st[0] = all.val();
  st[1] = all.val() > 0 ? dg.val() / all.val() : kutil::nan();
  st[2] = kutil::nan();
  if (off.val() > 0) {
    int q = 0;
    for (int a = 0; a < k; ++a) {
      for (int b = 0; b < k; ++b) q += a != b && cnt[a + static_cast<long long>(k) * b] > 0;
    }
    kutil::ksum th;
    for (int a = 0; a < k; ++a) {
      for (int b = 0; b < k; ++b) {
        const double v = cnt[a + static_cast<long long>(k) * b];
        if (a == b || v == 0) continue;
        const double s = v / off.val();
        th.add(s * std::log(s * q));
      }
    }
    st[2] = th.val();
  }
  return CLAB_OK;
}
