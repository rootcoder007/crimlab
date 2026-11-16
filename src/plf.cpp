// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Planted effects (crimlab.h, plntfx()): the outcome is shuffled, which
// removes whatever association it had with the treatment and the
// covariates, an effect of known size is added to the treated (y* =
// shuffled y + eff d), and the analysis (the coefficient of d in the least
// squares fit on [1, x, d], xls.h) is run again; over nrep shuffles its
// estimates should centre on the planted effect, and its intervals cover
// it at their level. The residual of y* on [1, x] is that of the shuffled
// y plus eff times that of d, so one projection serves every effect of a
// shuffle. One xorshift64* stream from the seed.
#include <cmath>
#include <cstddef>
#include <vector>

#include "crimlab.h"
#include "kutil.h"
#include "xls.h"

extern "C" int clab_plntfx(const double *x, int n, int p, const double *y, const double *d, const double *eff,
                           int ne, int nrep, int hc, double level, unsigned int seed, double *out) CLAB_TRY {
  using xls::z;
  if (n < 1 || p < 0 || ne < 1 || nrep < 2 || !(level > 0 && level < 1)) return CLAB_EARG;
  xls::Proj pr;
  if (!pr.make(x, n, nullptr, p, nullptr)) return CLAB_EDOM;
  std::vector<double> rd(z(n)), ys(y, y + n), rs(z(n)), ry(z(n));
  pr.res(d, rd.data());
  const double dss = xls::css(d, nullptr, n);
  {
    double b = 0, v = 0, df = 0;
    if (!xls::coef(rd.data(), rd.data(), nullptr, n, p + 2, hc != 0, dss, b, v, df)) return CLAB_EDOM;
  }
  // per effect: Welford's mean and sum of squares, the squared error,
  // the intervals covering eff, the tests rejecting 0
  std::vector<double> mn(z(ne), 0.0), m2(z(ne), 0.0), se2(z(ne), 0.0), cov(z(ne), 0.0), rej(z(ne), 0.0);
  kutil::rng r(seed);
  for (int s = 0; s < nrep; ++s) {
    if (kutil::stop(static_cast<double>(n) * (p + ne + 1))) return CLAB_EINTR;
    for (int i = n - 1; i > 0; --i) {
      const int j = r.below(i + 1);
      const double t = ys[z(i)];
      ys[z(i)] = ys[z(j)];
      ys[z(j)] = t;
    }
    pr.res(ys.data(), rs.data());
    for (int e = 0; e < ne; ++e) {
      for (int i = 0; i < n; ++i) ry[z(i)] = rs[z(i)] + eff[e] * rd[z(i)];
      double b = 0, v = 0, df = 0;
      xls::coef(rd.data(), ry.data(), nullptr, n, p + 2, hc != 0, dss, b, v, df);
      const double se = std::sqrt(v), q = clab_qt(0.5 + level / 2, df, 1);
      const double dl = b - mn[z(e)];
      mn[z(e)] += dl / (s + 1);
      m2[z(e)] += dl * (b - mn[z(e)]);
      se2[z(e)] += (b - eff[e]) * (b - eff[e]);
      if (std::fabs(b - eff[e]) <= q * se) cov[z(e)] += 1;
      if (std::fabs(b) > q * se) rej[z(e)] += 1;
    }
  }
  for (int e = 0; e < ne; ++e) {
    out[e] = mn[z(e)];
    out[ne + e] = std::sqrt(m2[z(e)] / (nrep - 1));
    out[2 * ne + e] = std::sqrt(se2[z(e)] / nrep);
    out[3 * ne + e] = cov[z(e)] / nrep;
    out[4 * ne + e] = rej[z(e)] / nrep;
  }
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)
