// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, internal benchmarking of officers (Ridgeway and MacDonald
// 2009): each officer's stops against the other officers' stops made in
// the same circumstances, by propensity weights.
#include <algorithm>
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

int clab_ibench(const double *y, const int *off, const double *x, int n, int p, int no, int minn,
                double *res, int *ok) {
  if (n < 2 || p < 0 || no < 2 || minn < 1) return CLAB_EARG;
  if (static_cast<long long>(n) * (p > 0 ? p : 1) > 2000000000LL) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (off[i] < 0 || off[i] >= no || !(y[i] == 0 || y[i] == 1)) return CLAB_EDOM;
  }
  const std::size_t un = static_cast<std::size_t>(n), up1 = static_cast<std::size_t>(p) + 1;
  std::vector<double> d(un), cf(up1), se(up1), w(2 * up1 * up1 + 4 * up1);
  std::vector<int> iw(un);
  // res: no x 7, column-major: stops, rate, benchmark, effective size of
  // the benchmark, z, p, q (Benjamini-Hochberg)
  for (std::size_t c = 0; c < 7 * static_cast<std::size_t>(no); ++c) res[c] = kutil::nan();
  for (int o = 0; o < no; ++o) {
    if (kutil::stop(static_cast<double>(n) * (p + 1) * (p + 1) * 10.0)) return CLAB_EINTR;
    ok[o] = 0;
    double nk = 0, s1 = 0;
    for (int i = 0; i < n; ++i) {
      d[static_cast<std::size_t>(i)] = off[i] == o ? 1.0 : 0.0;
      if (off[i] == o) {
        nk += 1;
        s1 += y[i];
      }
    }
    res[o] = nk;
    if (nk < minn || nk >= n) continue;
    res[no + o] = s1 / nk;
    double dev = 0.0;
    const int rc = clab_lgtreg(x, n, p, d.data(), 1e-10, 50, cf.data(), se.data(), &dev, w.data(), iw.data());
    if (rc == CLAB_EINTR) return rc;
    if (rc != CLAB_OK) continue;
    kutil::ksum sw, sw2, swy;
    for (int i = 0; i < n; ++i) {
      if (off[i] == o) continue;
      double eta = cf[0];
      for (int j = 0; j < p; ++j) eta += cf[static_cast<std::size_t>(j) + 1] * x[static_cast<std::size_t>(j) * un + i];
      const double wi = std::exp(eta);  // e / (1 - e)
      sw.add(wi);
      sw2.add(wi * wi);
      swy.add(wi * y[i]);
    }
    const double W = sw.val();
    if (!(W > 0) || !std::isfinite(W)) continue;
    const double b = swy.val() / W, ess = W * W / sw2.val(), r = res[no + o];
    const double v = r * (1 - r) / nk + b * (1 - b) / ess;
    res[2 * no + o] = b;
    res[3 * no + o] = ess;
    if (v > 0) {
      res[4 * no + o] = (r - b) / std::sqrt(v);
      res[5 * no + o] = std::erfc(std::fabs(res[4 * no + o]) / std::sqrt(2.0));
    }
    ok[o] = 1;
  }
  // Benjamini-Hochberg over the officers with a p-value
  std::vector<int> ix;
  for (int o = 0; o < no; ++o) {
    if (!std::isnan(res[5 * no + o])) ix.push_back(o);
  }
  const int m = static_cast<int>(ix.size());
  std::sort(ix.begin(), ix.end(), [&](int u, int v) {
    return res[5 * no + u] > res[5 * no + v] || (res[5 * no + u] == res[5 * no + v] && u < v);
  });
  double run = 1.0;
  for (int k = 0; k < m; ++k) {
    const int o = ix[static_cast<std::size_t>(k)];
    const double q = res[5 * no + o] * m / (m - k);
    run = std::fmin(run, q);
    res[6 * no + o] = run;
  }
  return CLAB_OK;
}
