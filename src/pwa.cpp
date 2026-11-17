// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The power of the one-way analysis of variance (crimlab.h, pwranv()):
// with k groups of n and Cohen's effect size f, the F statistic has df1 =
// k - 1, df2 = k (n - 1) and noncentrality lambda = k n f^2 (Cohen 1988,
// ch. 8), and the power is the chance that it passes the critical value
// of the central F at alpha. The upper tail of the noncentral F is the
// Poisson mixture of the upper tails of beta distributions (Johnson, Kotz
// and Balakrishnan 1995, ch. 30), clab_pnf() of src/dst.cpp, whose work
// does not grow with the noncentrality.
#include <cmath>

#include "crimlab.h"
#include "kutil.h"

namespace {

// One cell: k groups of n, effect f; out the 6 columns from df1.
void cell(int k, double n, double f, double alpha, double *o) {
  const double d1 = k - 1, d2 = k * (n - 1), l = k * n * f * f;
  const double fc = clab_qf(alpha, d1, d2, 0);
  o[0] = d1;
  o[1] = d2;
  o[2] = l;
  o[3] = fc;
  o[4] = clab_pnf(fc, d1, d2, l, 0);
}

}  // namespace

extern "C" int clab_pwranv(int k, int m, const double *n, const double *f, double alpha, double power,
                           double *out) CLAB_TRY {
  if (k < 2 || m < 0 || !(alpha > 0 && alpha < 1)) return CLAB_EARG;
  const double nmax = 1e7;
  for (int i = 0; i < m; ++i) {
    if (kutil::stop(1e3)) return CLAB_EINTR;
    if (!(f[i] >= 0) || !std::isfinite(f[i])) return CLAB_EARG;
    double o[5];
    double nn;
    if (n) {
      if (!(n[i] > 1) || !std::isfinite(n[i])) return CLAB_EARG;
      nn = n[i];
      cell(k, nn, f[i], alpha, o);
    } else {
      if (!(power > alpha && power < 1)) return CLAB_EARG;
      // the smallest whole n from 2 with power at least power: doubling,
      // then halving the bracket (the power grows with n when f > 0)
      cell(k, nmax, f[i], alpha, o);
      if (f[i] == 0 || o[4] < power) {
        nn = kutil::nan();
      } else {
        double a = 1, b = 2;
        for (;;) {
          cell(k, b, f[i], alpha, o);
          if (o[4] >= power || b >= nmax) break;
          a = b;
          b = std::fmin(2 * b, nmax);
        }
        while (b - a > 1) {
          const double c = std::floor((a + b) / 2);
          cell(k, c, f[i], alpha, o);
          if (o[4] >= power) {
            b = c;
          } else {
            a = c;
          }
        }
        nn = b;
        cell(k, nn, f[i], alpha, o);
      }
    }
    out[i] = nn;
    for (int j = 0; j < 5; ++j) out[(j + 1) * m + i] = std::isnan(nn) && j > 0 ? kutil::nan() : o[j];
    if (std::isnan(nn)) out[m + i] = k - 1;
  }
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)
