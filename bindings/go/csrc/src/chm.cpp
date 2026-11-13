// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, forensics: the Chumbley score for two striated tool-mark
// profiles, in the deterministic form of Hadler and Morris (2018): lowess
// detrending, a search for the best-matching pair of windows, and a rank
// test of same-shift against different-shift window correlations.
#include <cstdint>

#include "crimlab.h"
#include "kutil.h"

namespace {

// mean and root sum of squares of v[i..i+w), two passes
void wstat(const double *v, int i, int w, double *m, double *rs) {
  double s = 0.0;
  for (int k = 0; k < w; ++k) s += v[i + k];
  const double mu = s / w;
  double q = 0.0;
  for (int k = 0; k < w; ++k) q += (v[i + k] - mu) * (v[i + k] - mu);
  *m = mu;
  *rs = std::sqrt(q);
}

// Pearson correlation of a[i..i+w) and b[j..j+w), two passes
double wcor(const double *a, int i, const double *b, int j, int w) {
  double ma, ra, mb, rb;
  wstat(a, i, w, &ma, &ra);
  wstat(b, j, w, &mb, &rb);
  double s = 0.0;
  for (int k = 0; k < w; ++k) s += (a[i + k] - ma) * (b[j + k] - mb);
  return s / (ra * rb);
}

// R's round(): to the nearest integer, halves to even
inline int rint2(double x) {
  const double f = std::floor(x), d = x - f;
  if (d > 0.5) return static_cast<int>(f) + 1;
  if (d < 0.5) return static_cast<int>(f);
  return static_cast<int>(f) % 2 == 0 ? static_cast<int>(f) : static_cast<int>(f) + 1;
}

// trim to the central 98% (positions round(0.01 n) .. round(0.99 n),
// 1-based, a position 0 is dropped) and subtract the lowess smooth
// (f = coarse, 3 robustness steps, delta = 0.01 (m - 1)) of the profile
// against 1..m. out: the residuals; xw: m scratch (abscissae); lw: 4m.
int detrend(const double *v, int n, double coarse, double *out, double *xw, double *lw,
            int *m) {
  int s = rint2(0.01 * n), e = rint2(0.99 * n);
  if (s < 1) s = 1;
  const int len = e - s + 1;
  if (len < 1) return CLAB_EARG;
  for (int k = 0; k < len; ++k) {
    xw[k] = k + 1;
    out[k] = v[s - 1 + k];
  }
  const int rc = clab_lwss(xw, out, len, coarse, 3, 0.01 * (len - 1), lw + 3 * len, lw);
  if (rc != CLAB_OK) return rc;
  for (int k = 0; k < len; ++k) out[k] -= lw[3 * len + k];
  *m = len;
  return CLAB_OK;
}

}  // namespace

extern "C" {

int clab_chmb(const double *a, int na, const double *b, int nb, int wopt,
              int wval, double coarse, int *nsame, int *ndiff, double *u,
              double *p, double *work) {
  if (na < 1 || nb < 1 || wopt < 2 || wval < 2 || wval > wopt) return CLAB_EARG;
  if (!(coarse > 0.0) || !std::isfinite(coarse)) return CLAB_EDOM;
  for (int i = 0; i < na; ++i)
    if (!std::isfinite(a[i])) return CLAB_EDOM;
  for (int i = 0; i < nb; ++i)
    if (!std::isfinite(b[i])) return CLAB_EDOM;
  const int nmax = na > nb ? na : nb;
  double *y1 = work, *y2 = y1 + na, *xw = y2 + nb, *lw = xw + nmax;
  double *m1 = lw + 4 * nmax, *r1 = m1 + na, *m2 = r1 + na, *r2 = m2 + nb;
  double *vals = r2 + nb;  // 2 nmax validation correlations
  int n1, n2;
  int rc = detrend(a, na, coarse, y1, xw, lw, &n1);
  if (rc != CLAB_OK) return rc;
  rc = detrend(b, nb, coarse, y2, xw, lw, &n2);
  if (rc != CLAB_OK) return rc;
  // a trimmed profile shorter than the optimisation window has no window
  if (n1 < wopt || n2 < wopt) return CLAB_EARG;
  const int c1 = n1 - wopt + 1, c2 = n2 - wopt + 1;
  // a window whose spread is at the rounding level of its profile (root
  // sum of squares below 64 ulps of the profile's largest magnitude times
  // sqrt(wopt)) is flat: its correlation is noise, so it is left out
  double amax = 0.0, bmax = 0.0;
  for (int i = 0; i < na; ++i) amax = std::fabs(a[i]) > amax ? std::fabs(a[i]) : amax;
  for (int i = 0; i < nb; ++i) bmax = std::fabs(b[i]) > bmax ? std::fabs(b[i]) : bmax;
  const double e64 = 64.0 * std::numeric_limits<double>::epsilon() * std::sqrt(static_cast<double>(wopt));
  for (int i = 0; i < c1; ++i) {
    wstat(y1, i, wopt, &m1[i], &r1[i]);
    if (r1[i] <= e64 * amax) r1[i] = 0.0;
  }
  for (int j = 0; j < c2; ++j) {
    wstat(y2, j, wopt, &m2[j], &r2[j]);
    if (r2[j] <= e64 * bmax) r2[j] = 0.0;
  }
  // optimisation step: the pair (window j of mark 2, window i of mark 1)
  // of largest correlation; among equal values the first in column-major
  // order (smallest i, then smallest j). Windows of zero variance have no
  // correlation. All c1 c2 pairs, one lag j - i at a time with sliding
  // cross products (recomputed exactly every 128 steps): O(c1 c2).
  double best = -HUGE_VAL;
  int bi = -1, bj = -1;
  for (int lag = -(c1 - 1); lag <= c2 - 1; ++lag) {
    const int i0 = lag < 0 ? -lag : 0;
    const int i1 = c1 - 1 < c2 - 1 - lag ? c1 - 1 : c2 - 1 - lag;
    double s = 0.0;
    for (int i = i0; i <= i1; ++i) {
      const int j = i + lag;
      if ((i - i0) % 128 == 0) {
        s = 0.0;
        for (int k = 0; k < wopt; ++k) s += y1[i + k] * y2[j + k];
      } else {
        s += y1[i + wopt - 1] * y2[j + wopt - 1] - y1[i - 1] * y2[j - 1];
      }
      if (!(r1[i] > 0.0) || !(r2[j] > 0.0)) continue;
      const double c = (s - wopt * m1[i] * m2[j]) / (r1[i] * r2[j]);
      if (c > best || (c == best && (i < bi || (i == bi && j < bj)))) {
        best = c;
        bi = i;
        bj = j;
      }
    }
  }
  if (bi < 0) return CLAB_EDOM;  // every window is flat
  // validation step, with the field's 1-based indices: rows are windows of
  // mark 2, columns windows of mark 1, nr and nc their numbers
  const int nr = n2 - wval + 1, nc = n1 - wval + 1;
  const int r0 = bj + 1, col0 = bi + 1;
  int ns = 0, nd = 0;
  double *sv = vals, *dv = vals + nmax;
  // same shift: forward from the end of the best window, then backward
  for (int r = r0 + (wopt - wval), c = col0 + (wopt - wval); r + wval < nr && c + wval < nc;) {
    r += wval;
    c += wval;
    sv[ns++] = wcor(y2, r - 1, y1, c - 1, wval);
  }
  for (int r = r0, c = col0; r - wval > 0 && c - wval > 0;) {
    r -= wval;
    c -= wval;
    sv[ns++] = wcor(y2, r - 1, y1, c - 1, wval);
  }
  // different shift: along the anti-diagonal through the best pair
  for (int r = r0 + (wopt - wval), c = col0; r + wval < nr && c - wval > 0;) {
    r += wval;
    c -= wval;
    dv[nd++] = wcor(y2, r - 1, y1, c - 1, wval);
  }
  for (int r = r0, c = col0 + (wopt - wval); r - wval > 0 && c + wval < nc;) {
    r -= wval;
    c += wval;
    dv[nd++] = wcor(y2, r - 1, y1, c - 1, wval);
  }
  *nsame = ns;
  *ndiff = nd;
  if (ns == 0 || nd == 0) {
    *u = *p = kutil::nan();
    return CLAB_OK;
  }
  // rank-sum of the same-shift correlations among all (midranks), its
  // standardisation with the tie-corrected variance, upper normal tail
  const int nn = ns + nd;
  kutil::ksum t, r2s;
  for (int g = 0; g < nn; ++g) {
    const double x = g < ns ? sv[g] : dv[g - ns];
    double less = 0.0, eq = 0.0;
    for (int h = 0; h < nn; ++h) {
      const double v = h < ns ? sv[h] : dv[h - ns];
      less += v < x;
      eq += v == x;
    }
    const double rk = less + (eq + 1.0) / 2.0;
    if (g < ns) t.add(rk);
    r2s.add(rk * rk);
  }
  const double n = ns, m = nd, nt = nn;
  const double var = n * m / (nt * (nt - 1.0)) * r2s.val() - n * m * (nt + 1.0) * (nt + 1.0) / (4.0 * (nt - 1.0));
  const double z = (t.val() - n * (nt + 1.0) / 2.0) / std::sqrt(var);
  *u = z;
  *p = 0.5 * std::erfc(z / std::sqrt(2.0));
  return CLAB_OK;
}

}  // extern "C"
