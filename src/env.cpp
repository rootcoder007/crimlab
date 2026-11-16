// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Environmental crime (crimlab.h): the additive log-mean Divisia index
// (LMDI) decomposition of a change or a gap in a count over the factors of
// a Kaya identity, the variance inflation factors of a design, the
// Jarque-Bera test of residuals, and Pesaran's CD test of cross-sectional
// dependence in a panel. No allocation, no exceptions, no I/O.
#include <cstdint>
#include <cmath>
#include <cstddef>

#include "crimlab.h"
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
  if (!(m2 > 0.0)) return CLAB_EDOM;
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
      if (!(saa > 0.0) || !(sbb > 0.0)) continue;
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
