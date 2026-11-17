// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Environmental crime (crimlab.h): the additive log-mean Divisia index
// (LMDI) decomposition of a change or a gap in a count over the factors of
// a Kaya identity, the variance inflation factors of a design, the
// Jarque-Bera test of residuals, Pesaran's CD test of cross-sectional
// dependence in a panel, and least squares with its diagnostics (olsdx).
// No I/O; only clab_olsdx() allocates, and returns CLAB_ENOMEM when that
// fails.
#include <cstdint>
#include <cmath>
#include <cstddef>
#include <vector>

#include "crimlab.h"
#include "kutil.h"
#include "pev.h"
#include "kutil.h"

namespace {

// the logarithmic mean (a - b) / (ln a - ln b), a when a = b; log1p keeps
// it exact when a and b are close
double lmean(double a, double b) {
  if (a == b) return a;
  return (a - b) / std::log1p((a - b) / b);
}

bool pos(double v) { return std::isfinite(v) && v > 0.0; }

}  // namespace

extern "C" int clab_lmdi(const double *y0, const double *y1, const double *f0,
                         const double *f1, int n, int k, double *out) CLAB_TRY {
  if (n < 1 || k < 1) return CLAB_EARG;
  const std::size_t sn = static_cast<std::size_t>(n), sk = static_cast<std::size_t>(k);
  for (std::size_t i = 0; i < sn; ++i) {
    if (!pos(y0[i]) || !pos(y1[i])) return CLAB_EDOM;
  }
  for (std::size_t q = 0; q < sn * sk; ++q) {
    if (!pos(f0[q]) || !pos(f1[q])) return CLAB_EDOM;
  }
  for (std::size_t i = 0; i < sn; ++i) {
    const double l = lmean(y1[i], y0[i]);
    for (std::size_t j = 0; j < sk; ++j) {
      out[j * sn + i] = l * std::log(f1[j * sn + i] / f0[j * sn + i]);
    }
  }
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

extern "C" int clab_vif(const double *x, int n, int p, double *out, double *work) CLAB_TRY {
  if (p < 1 || n < p + 2) return CLAB_EARG;
  const std::size_t sn = static_cast<std::size_t>(n), sp = static_cast<std::size_t>(p);
  double *d = work;         // the constant and the other columns: n p
  double *b = d + sn * sp;  // p
  double *xi = b + sp;      // p p
  double *r = xi + sp * sp;  // n
  double *wk = r + sn;      // p p + p
  for (std::size_t j = 0; j < sp; ++j) {
    const double *c = x + j * sn;
    double m = 0.0, tss = 0.0;
    for (std::size_t i = 0; i < sn; ++i) m += c[i];
    m /= n;
    for (std::size_t i = 0; i < sn; ++i) tss += (c[i] - m) * (c[i] - m);
    if (!(tss > 0.0)) {
      out[j] = NAN;
      continue;
    }
    if (p == 1) {
      out[j] = 1.0;
      continue;
    }
    for (std::size_t i = 0; i < sn; ++i) d[i] = 1.0;
    std::size_t col = 1;
    for (std::size_t h = 0; h < sp; ++h) {
      if (h == j) continue;
      for (std::size_t i = 0; i < sn; ++i) d[col * sn + i] = x[h * sn + i];
      ++col;
    }
    if (pev::ols(d, c, n, p, nullptr, b, xi, r, wk) != CLAB_OK) {
      out[j] = INFINITY;
      continue;
    }
    double rss = 0.0;
    for (std::size_t i = 0; i < sn; ++i) rss += r[i] * r[i];
    out[j] = rss > 0.0 ? tss / rss : INFINITY;
  }
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

extern "C" int clab_jb(const double *r, int n, double *out) CLAB_TRY {
  if (n < 4) return CLAB_EARG;
  double m = 0.0;
  for (int i = 0; i < n; ++i) m += r[i];
  m /= n;
  double m2 = 0.0, m3 = 0.0, m4 = 0.0;
  for (int i = 0; i < n; ++i) {
    const double e = r[i] - m, e2 = e * e;
    m2 += e2;
    m3 += e2 * e;
    m4 += e2 * e2;
  }
  m2 /= n;
  m3 /= n;
  m4 /= n;
  if (kutil::cnst(r, n) || !(m2 > 0.0)) return CLAB_EDOM;
  const double s = m3 / (m2 * std::sqrt(m2)), k = m4 / (m2 * m2);
  const double jb = n / 6.0 * (s * s + (k - 3.0) * (k - 3.0) / 4.0);
  out[0] = s;
  out[1] = k;
  out[2] = jb;
  out[3] = std::exp(-jb / 2.0);  // the upper tail of chi-squared on 2 df
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

extern "C" int clab_pcd(const double *e, int nu, int nt, double *out) CLAB_TRY {
  if (nu < 2 || nt < 2) return CLAB_EARG;
  const std::size_t st = static_cast<std::size_t>(nt);
  double cd = 0.0, sr = 0.0, sa = 0.0;
  std::int64_t m = 0;
  for (int a = 0; a < nu; ++a) {
    const double *ea = e + static_cast<std::size_t>(a) * st;
    for (int b = a + 1; b < nu; ++b) {
      const double *eb = e + static_cast<std::size_t>(b) * st;
      double ma = 0.0, mb = 0.0;
      int t = 0;
      for (std::size_t s = 0; s < st; ++s) {
        if (std::isnan(ea[s]) || std::isnan(eb[s])) continue;
        ma += ea[s];
        mb += eb[s];
        ++t;
      }
      if (t < 2) continue;
      ma /= t;
      mb /= t;
      double sab = 0.0, saa = 0.0, sbb = 0.0;
      for (std::size_t s = 0; s < st; ++s) {
        if (std::isnan(ea[s]) || std::isnan(eb[s])) continue;
        sab += (ea[s] - ma) * (eb[s] - mb);
        saa += (ea[s] - ma) * (ea[s] - ma);
        sbb += (eb[s] - mb) * (eb[s] - mb);
      }
      // a unit constant over the shared periods has no correlation (its
      // sum of squares about a rounded mean need not be 0)
      double fa = 0.0, fb = 0.0;
      bool va = false, vb = false, first = true;
      for (std::size_t s = 0; s < st; ++s) {
        if (std::isnan(ea[s]) || std::isnan(eb[s])) continue;
        if (first) {
          fa = ea[s];
          fb = eb[s];
          first = false;
        }
        va = va || ea[s] != fa;
        vb = vb || eb[s] != fb;
      }
      if (!va || !vb || !(saa > 0.0) || !(sbb > 0.0)) continue;
      const double rho = sab / std::sqrt(saa * sbb);
      cd += std::sqrt(static_cast<double>(t)) * rho;
      sr += rho;
      sa += std::fabs(rho);
      ++m;
    }
  }
  if (m == 0) return CLAB_EDOM;
  const double dm = static_cast<double>(m);
  out[0] = cd / std::sqrt(dm);
  out[1] = std::erfc(std::fabs(out[0]) / std::sqrt(2.0));
  out[2] = sr / dm;
  out[3] = sa / dm;
  out[4] = dm;
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

extern "C" int clab_olsdx(const double *y, const double *x, int n, int p, int icpt,
                          const double *at, int m, int pnew, double lv, double *co,
                          double *vc, double *cr, double *cs, double *pr,
                          double *st) CLAB_TRY {
  const int k = p + (icpt ? 1 : 0);
  if (p < 1 || n < k + 2 || m < 0 || (m > 0 && !at) || !(lv > 0.0 && lv < 1.0)) {
    return CLAB_EARG;
  }
  const std::size_t sn = static_cast<std::size_t>(n), sp = static_cast<std::size_t>(p),
                    sk = static_cast<std::size_t>(k), sm = static_cast<std::size_t>(m),
                    c0 = icpt ? 1 : 0;
  for (std::size_t i = 0; i < sn; ++i) {
    if (!std::isfinite(y[i])) return CLAB_EDOM;
  }
  for (std::size_t q = 0; q < sn * sp; ++q) {
    if (!std::isfinite(x[q])) return CLAB_EDOM;
  }
  for (std::size_t q = 0; q < sm * sp; ++q) {
    if (!std::isfinite(at[q])) return CLAB_EDOM;
  }
  // the design: the constant (with icpt), then the columns of x, centred
  // (with icpt; y too) and scaled by a power of 2 (exact): the fit is
  // equivariant to both, and the cross products of uncentred columns lose
  // the digits of their spread (x + 1e6 moved the slope by 3e-4, and x +
  // 1e7 was taken for collinear where lm() fits). b and xi are those of
  // the centred, scaled columns; ba and the variance vc are mapped back.
  std::vector<double> cx(sp, 0.0), sc(sp, 1.0), yc(sn);
  double cy = 0.0;
  for (std::size_t i = 0; i < sn && icpt; ++i) cy += (y[i] - cy) / static_cast<double>(i + 1);
  for (std::size_t j = 0; j < sp; ++j) {
    for (std::size_t i = 0; i < sn && icpt; ++i) cx[j] += (x[j * sn + i] - cx[j]) / static_cast<double>(i + 1);
    double big = 0.0;
    for (std::size_t i = 0; i < sn; ++i) big = std::fmax(big, std::fabs(x[j * sn + i] - cx[j]));
    int ex = 0;
    std::frexp(big, &ex);
    if (big > 0.0) sc[j] = std::ldexp(1.0, ex);
  }
  std::vector<double> d(sn * sk), b(sk), xi(sk * sk), r(sn), wk(sk * sk + sk);
  for (std::size_t i = 0; i < sn; ++i) yc[i] = y[i] - cy;
  for (std::size_t i = 0; i < sn && icpt; ++i) d[i] = 1.0;
  for (std::size_t j = 0; j < sp; ++j) {
    for (std::size_t i = 0; i < sn; ++i) d[(c0 + j) * sn + i] = (x[j * sn + i] - cx[j]) / sc[j];
  }
  if (pev::ols(d.data(), yc.data(), n, k, nullptr, b.data(), xi.data(), r.data(), wk.data()) !=
      CLAB_OK) {
    return CLAB_EDOM;
  }
  // back to the columns as given: beta = A b (+ cy in the constant), A
  // the identity but 1 / sc on the slopes and -cx / sc in the constant's
  // row; vc = A (s2 xi) A'
  std::vector<double> am(sk * sk, 0.0), ba(sk);
  for (std::size_t j = 0; j < sk; ++j) am[j * sk + j] = j < c0 ? 1.0 : 1.0 / sc[j - c0];
  for (std::size_t j = 0; j < sp && icpt; ++j) am[(c0 + j) * sk] = -cx[j] / sc[j];
  for (std::size_t a = 0; a < sk; ++a) {
    ba[a] = a < c0 ? cy : 0.0;
    for (std::size_t c = 0; c < sk; ++c) ba[a] += am[c * sk + a] * b[c];
  }
  double rss = 0.0, my = 0.0, tss = 0.0, sy = 0.0;
  for (std::size_t i = 0; i < sn; ++i) {
    rss += r[i] * r[i];
    my += y[i] / n;
  }
  for (std::size_t i = 0; i < sn; ++i) {
    sy += (y[i] - my) * (y[i] - my);
    tss += icpt ? (y[i] - my) * (y[i] - my) : y[i] * y[i];
  }
  sy = std::sqrt(sy / (n - 1));
  const double df2 = n - k, s2 = rss / df2, s = std::sqrt(s2);
  const double tq = clab_qt(1.0 - (1.0 - lv) / 2.0, df2, 1);
  for (std::size_t a = 0; a < sk; ++a) {
    for (std::size_t c = 0; c < sk; ++c) {
      double v = 0.0;
      for (std::size_t e = 0; e < sk; ++e) {
        for (std::size_t f = 0; f < sk; ++f) v += am[e * sk + a] * xi[f * sk + e] * am[f * sk + c];
      }
      vc[c * sk + a] = s2 * v;
    }
  }
  // the means and standard deviations of the columns of x
  std::vector<double> mx(sp, 0.0), sx(sp, 0.0);
  for (std::size_t j = 0; j < sp; ++j) {
    for (std::size_t i = 0; i < sn; ++i) mx[j] += x[j * sn + i] / n;
    for (std::size_t i = 0; i < sn; ++i) {
      sx[j] += (x[j * sn + i] - mx[j]) * (x[j * sn + i] - mx[j]);
    }
    sx[j] = std::sqrt(sx[j] / (n - 1));
  }
  std::vector<double> vf(sp), vw(sn * (sp + 1) + 2 * sp * (sp + 1));
  const bool vok = clab_vif(x, n, p, vf.data(), vw.data()) == CLAB_OK;
  // coefficients: est, se, t, p, lo, hi, beta, tol, vif (k x 9)
  for (std::size_t j = 0; j < sk; ++j) {
    const double e = ba[j], se = std::sqrt(vc[j * sk + j]), t = e / se;
    co[j] = e;
    co[sk + j] = se;
    co[2 * sk + j] = t;
    co[3 * sk + j] = 2.0 * clab_pt(-std::fabs(t), df2, 1);
    co[4 * sk + j] = e - tq * se;
    co[5 * sk + j] = e + tq * se;
    if (j < c0) {
      co[6 * sk + j] = co[7 * sk + j] = co[8 * sk + j] = NAN;
      continue;
    }
    const std::size_t h = j - c0;
    const double v = vok ? vf[h] : NAN;
    co[6 * sk + j] = e * sx[h] / sy;
    co[7 * sk + j] = 1.0 / v;
    co[8 * sk + j] = v;
  }
  // the correlations of the columns of x (p x p); NaN for a constant one
  for (std::size_t a = 0; a < sp; ++a) {
    for (std::size_t c = 0; c < sp; ++c) {
      double sab = 0.0;
      for (std::size_t i = 0; i < sn; ++i) {
        sab += (x[a * sn + i] - mx[a]) * (x[c * sn + i] - mx[c]);
      }
      cr[c * sp + a] = sab / (n - 1) / (sx[a] * sx[c]);
    }
  }
  // the rows: fit, resid, zres, sres, tres, hat, cook (n x 7)
  for (std::size_t i = 0; i < sn; ++i) {
    double h = 0.0;
    for (std::size_t a = 0; a < sk; ++a) {
      double u = 0.0;
      for (std::size_t c = 0; c < sk; ++c) u += xi[c * sk + a] * d[c * sn + i];
      h += d[a * sn + i] * u;
    }
    const double e = r[i], sr = e / (s * std::sqrt(1.0 - h));
    cs[i] = y[i] - e;
    cs[sn + i] = e;
    cs[2 * sn + i] = e / s;
    cs[3 * sn + i] = sr;
    cs[4 * sn + i] = sr * std::sqrt((df2 - 1.0) / (df2 - sr * sr));
    cs[5 * sn + i] = h;
    cs[6 * sn + i] = sr * sr * h / (k * (1.0 - h));
  }
  double dw = 0.0;
  for (std::size_t i = 1; i < sn; ++i) dw += (r[i] - r[i - 1]) * (r[i] - r[i - 1]);
  const double r2 = 1.0 - rss / tss, df1 = k - static_cast<double>(c0);
  const double f = (tss - rss) / df1 / s2;
  st[0] = std::sqrt(r2 > 0.0 ? r2 : 0.0);
  st[1] = r2;
  st[2] = 1.0 - (1.0 - r2) * (n - static_cast<double>(c0)) / df2;
  st[3] = s;
  st[4] = f;
  st[5] = df1;
  st[6] = df2;
  st[7] = clab_pf(f, df1, df2, 0);
  st[8] = dw / rss;
  st[9] = rss;
  // the points at (m x p): the fitted mean, its standard error (of a new
  // value with pnew) and interval (m x 4)
  std::vector<double> x0(sk);
  for (std::size_t a = 0; a < sm; ++a) {
    // in the centred, scaled columns, as fitted
    if (icpt) x0[0] = 1.0;
    for (std::size_t j = 0; j < sp; ++j) x0[c0 + j] = (at[j * sm + a] - cx[j]) / sc[j];
    double e = cy, v = 0.0;
    for (std::size_t j = 0; j < sk; ++j) {
      e += x0[j] * b[j];
      for (std::size_t c = 0; c < sk; ++c) v += x0[j] * xi[c * sk + j] * x0[c];
    }
    v = s2 * (v + (pnew ? 1.0 : 0.0));
    pr[a] = e;
    pr[sm + a] = std::sqrt(v);
    pr[2 * sm + a] = e - tq * pr[sm + a];
    pr[3 * sm + a] = e + tq * pr[sm + a];
  }
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)
