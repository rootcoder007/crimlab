// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Falsification with negative-control outcomes (crimlab.h, falsctl()):
// the coefficient of the treatment d in the regression of each outcome on
// [1, x, d], by Frisch-Waugh-Lovell (xls.h), with its classical or HC1
// standard error, Student's t on the residual degrees of freedom, and
// optionally the permutation p-value of Kennedy (1995): the residuals of
// the outcome on [1, x] permuted against those of d, one permutation of
// the rows for every outcome at once. Outcomes that d cannot affect
// should show no association; one that does points to confounding or
// bias (Lipsitch, Tchetgen Tchetgen and Cohen 2010). The p-values of the
// controls are adjusted for their number (clab_pbadj).
#include <cmath>
#include <cstddef>
#include <vector>

#include "crimlab.h"
#include "kutil.h"
#include "xls.h"

extern "C" int clab_falsctl(const double *x, int n, int p, const double *d, const double *y, int m,
                            const int *ctl, int hc, double level, int adj, int nperm, unsigned int seed,
                            double *out) CLAB_TRY {
  using xls::z;
  if (n < 1 || p < 0 || m < 1 || nperm < 0 || adj < 0 || adj > 5 || !(level > 0 && level < 1)) {
    return CLAB_EARG;
  }
  xls::Proj pr;
  if (!pr.make(x, n, nullptr, p, nullptr)) return CLAB_EDOM;
  std::vector<double> rd(z(n)), ry(z(n) * z(m));
  pr.res(d, rd.data());
  const double dss = xls::css(d, nullptr, n);
  std::vector<double> b(z(m));
  for (int j = 0; j < m; ++j) {
    double *r = &ry[z(j) * z(n)];
    pr.res(y + z(j) * z(n), r);
    double v = 0, df = 0;
    if (!xls::coef(rd.data(), r, nullptr, n, p + 2, hc != 0, dss, b[z(j)], v, df)) return CLAB_EDOM;
    const double se = std::sqrt(v), t = b[z(j)] / se, q = clab_qt(0.5 + level / 2, df, 1);
    out[j] = b[z(j)];
    out[m + j] = se;
    out[2 * m + j] = t;
    out[3 * m + j] = 2 * clab_pt(-std::fabs(t), df, 1);
    out[4 * m + j] = b[z(j)] - q * se;
    out[5 * m + j] = b[z(j)] + q * se;
    out[6 * m + j] = kutil::nan();
    out[7 * m + j] = kutil::nan();
  }
  if (nperm > 0) {
    double srr = 0;
    for (double v : rd) srr += v * v;
    std::vector<double> k(z(m), 0.0);
    std::vector<int> ix(z(n));
    kutil::rng r(seed);
    for (int s = 0; s < nperm; ++s) {
      if (kutil::stop(static_cast<double>(n) * m)) return CLAB_EINTR;
      for (int i = 0; i < n; ++i) ix[z(i)] = i;
      for (int i = n - 1; i > 0; --i) {
        const int j = r.below(i + 1);
        const int t = ix[z(i)];
        ix[z(i)] = ix[z(j)];
        ix[z(j)] = t;
      }
      for (int j = 0; j < m; ++j) {
        const double *yr = &ry[z(j) * z(n)];
        double s1 = 0;
        for (int i = 0; i < n; ++i) s1 += rd[z(i)] * yr[ix[z(i)]];
        const double bs = s1 / srr, bo = b[z(j)];
        if (std::fabs(bs) >= std::fabs(bo) * (1 - 1e-9)) k[z(j)] += 1;
      }
    }
    for (int j = 0; j < m; ++j) out[6 * m + j] = (1 + k[z(j)]) / (1.0 + nperm);
  }
  // the controls' p-values (permutation ones when drawn), adjusted
  std::vector<double> pc, qc;
  std::vector<int> o;
  for (int j = 0; j < m; ++j) {
    if (ctl[j]) pc.push_back(out[(nperm > 0 ? 6 : 3) * m + j]);
  }
  qc.resize(pc.size());
  o.resize(pc.size());
  const int rc = clab_pbadj(pc.data(), static_cast<int>(pc.size()), adj, qc.data(), o.data());
  if (rc != CLAB_OK) return rc;
  for (int j = 0, c = 0; j < m; ++j) {
    if (ctl[j]) out[7 * m + j] = qc[z(c++)];
  }
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)
