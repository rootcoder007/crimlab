// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The specification curve (crimlab.h, specrv()) of Simonsohn, Simmons and
// Nelson (2020): the coefficient of the treatment d in the least squares
// fit (xls.h) of each outcome on [1, x, the optional covariates in a
// subset, d], for every outcome and every subset, with the joint
// inference of their sec. 3: three summaries of the whole curve (the
// median estimate, the specifications significant in the dominant
// direction, and the sum of the t statistics over the square root of
// their number) are compared with their values on data where the null
// holds, made by permuting d (an experiment) or by drawing rows with
// replacement after removing each specification's estimate from its
// outcome (y - b d). One xorshift64* stream from the seed.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "crimlab.h"
#include "kutil.h"
#include "xls.h"

namespace {

using xls::z;

struct Fit {
  std::vector<double> b, se, df;
};

// The columns of subset s: the nx fixed ones, then the optional ones whose
// bit is set.
std::vector<int> cols(int nx, int m, int s) {
  std::vector<int> c;
  for (int j = 0; j < nx; ++j) c.push_back(j);
  for (int j = 0; j < m; ++j) {
    if ((s >> j) & 1) c.push_back(nx + j);
  }
  return c;
}

// The median, the specifications significant (p < 1 - level) with the
// sign sg, and the sum of the t statistics over sqrt(K).
void summ(const Fit &f, double level, double sg, double *o) {
  const std::size_t k = f.b.size();
  std::vector<double> v(f.b);
  std::sort(v.begin(), v.end());
  o[0] = k % 2 ? v[k / 2] : (v[k / 2 - 1] + v[k / 2]) / 2;
  double ns = 0, t = 0;
  for (std::size_t i = 0; i < k; ++i) {
    const double st = f.b[i] / f.se[i], p = 2 * clab_pt(-std::fabs(st), f.df[i], 1);
    if (p < 1 - level && f.b[i] * sg > 0) ns += 1;
    t += st;
  }
  o[1] = ns;
  o[2] = t / std::sqrt(static_cast<double>(k));
}

}  // namespace

extern "C" int clab_specrv(const double *x, int n, int nx, int m, const double *y, int ny, const double *d,
                           int hc, double level, int null, int nsim, unsigned int seed, double *spec, double *st,
                           double *dist) CLAB_TRY {
  if (n < 1 || nx < 0 || m < 0 || m > 10 || ny < 1 || null < 0 || null > 1 || nsim < 1 ||
      !(level > 0 && level < 1)) {
    return CLAB_EARG;
  }
  const int ns = 1 << m, K = ny * ns;
  const double cells = static_cast<double>(ns) * n * (nx + m + 1) + static_cast<double>(K) * n;
  if (cells > kutil::maxcell) return CLAB_EARG;
  // the observed fits, keeping each subset's projection and each
  // specification's outcome residuals for the permutations
  std::vector<xls::Proj> pr(z(ns));
  std::vector<double> ry(z(K) * z(n)), rd(z(n));
  Fit obs;
  const double dss = xls::css(d, nullptr, n);
  for (int s = 0; s < ns; ++s) {
    const std::vector<int> c = cols(nx, m, s);
    if (!pr[z(s)].make(x, n, c.data(), static_cast<int>(c.size()), nullptr)) return CLAB_EDOM;
    pr[z(s)].res(d, rd.data());
    for (int j = 0; j < ny; ++j) {
      double *r = &ry[(z(j) * z(ns) + z(s)) * z(n)];
      pr[z(s)].res(y + z(j) * z(n), r);
      double b = 0, v = 0, df = 0;
      if (!xls::coef(rd.data(), r, nullptr, n, static_cast<int>(c.size()) + 2, hc != 0, dss, b, v, df)) {
        return CLAB_EDOM;
      }
      obs.b.push_back(b);
      obs.se.push_back(std::sqrt(v));
      obs.df.push_back(df);
    }
  }
  // the specifications in output order: outcome, then subset
  std::vector<double> bo(z(K)), so(z(K)), fo(z(K));
  for (int s = 0, k = 0; s < ns; ++s) {
    for (int j = 0; j < ny; ++j, ++k) {
      const int q = j * ns + s;
      bo[z(q)] = obs.b[z(k)];
      so[z(q)] = obs.se[z(k)];
      fo[z(q)] = obs.df[z(k)];
    }
  }
  obs.b = bo;
  obs.se = so;
  obs.df = fo;
  for (int q = 0; q < K; ++q) {
    const double t = bo[z(q)] / so[z(q)], c = clab_qt(0.5 + level / 2, fo[z(q)], 1);
    spec[q] = bo[z(q)];
    spec[K + q] = so[z(q)];
    spec[2 * K + q] = t;
    spec[3 * K + q] = 2 * clab_pt(-std::fabs(t), fo[z(q)], 1);
    spec[4 * K + q] = bo[z(q)] - c * so[z(q)];
    spec[5 * K + q] = bo[z(q)] + c * so[z(q)];
  }
  double o0[3];
  {
    std::vector<double> v(bo);
    std::sort(v.begin(), v.end());
    const double med = K % 2 ? v[z(K / 2)] : (v[z(K / 2 - 1)] + v[z(K / 2)]) / 2;
    st[6] = med >= 0 ? 1 : -1;
  }
  const double sg = st[6];
  summ(obs, level, sg, o0);
  // the resamples
  kutil::rng r(seed);
  std::vector<double> dp(d, d + n), w(z(n)), rw(z(n));
  Fit f;
  f.b.resize(z(K));
  f.se.resize(z(K));
  f.df.resize(z(K));
  double ge[3] = {0, 0, 0}, nok = 0;
  for (int i = 0; i < nsim; ++i) {
    if (kutil::stop(static_cast<double>(n) * K * (nx + m + 2))) return CLAB_EINTR;
    bool ok = true;
    if (null == 0) {
      for (int a = n - 1; a > 0; --a) {
        const int b = r.below(a + 1);
        const double t = dp[z(a)];
        dp[z(a)] = dp[z(b)];
        dp[z(b)] = t;
      }
      for (int s = 0; s < ns && ok; ++s) {
        pr[z(s)].res(dp.data(), rd.data());
        const int k = pr[z(s)].k;
        for (int j = 0; j < ny && ok; ++j) {
          const int q = j * ns + s;
          double v = 0;
          ok = xls::coef(rd.data(), &ry[(z(j) * z(ns) + z(s)) * z(n)], nullptr, n, k + 2, hc != 0, dss,
                         f.b[z(q)], v, f.df[z(q)]);
          f.se[z(q)] = std::sqrt(v);
        }
      }
    } else {
      std::fill(w.begin(), w.end(), 0.0);
      for (int a = 0; a < n; ++a) w[z(r.below(n))] += 1;
      const double dw = xls::css(d, w.data(), n);
      xls::Proj pw;
      for (int s = 0; s < ns && ok; ++s) {
        const std::vector<int> c = cols(nx, m, s);
        ok = pw.make(x, n, c.data(), static_cast<int>(c.size()), w.data());
        if (!ok) break;
        pw.res(d, rd.data());
        for (int j = 0; j < ny && ok; ++j) {
          const int q = j * ns + s;
          pw.res(y + z(j) * z(n), rw.data());
          double b = 0, v = 0;
          ok = xls::coef(rd.data(), rw.data(), w.data(), n, static_cast<int>(c.size()) + 2, hc != 0, dw, b, v,
                         f.df[z(q)]);
          // the null imposed: the estimate on y - b_obs d is b - b_obs,
          // with the same residuals
          f.b[z(q)] = b - bo[z(q)];
          f.se[z(q)] = std::sqrt(v);
          ok = ok && f.se[z(q)] > 0;
        }
      }
    }
    if (!ok) {
      for (int a = 0; a < 3; ++a) dist[a * nsim + i] = kutil::nan();
      continue;
    }
    double o[3];
    summ(f, level, sg, o);
    nok += 1;
    for (int a = 0; a < 3; ++a) dist[a * nsim + i] = o[a];
    if (sg * o[0] >= sg * o0[0] - 1e-9 * std::fabs(o0[0])) ge[0] += 1;
    if (o[1] >= o0[1]) ge[1] += 1;
    if (sg * o[2] >= sg * o0[2] - 1e-9 * std::fabs(o0[2])) ge[2] += 1;
  }
  for (int a = 0; a < 3; ++a) {
    st[a] = o0[a];
    st[3 + a] = (1 + ge[a]) / (1 + nok);
  }
  st[7] = nok;
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)
