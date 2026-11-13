// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Design-based variance by linearisation for stratified multistage
// samples: the variance of a total of influence values from the
// first-stage units (PSUs) within strata, with replacement or with the
// finite-population correction (Wolter 2007, ch. 6; the estimator of R's
// survey package, Lumley 2004).
#include <cmath>
#include <cstddef>
#include <vector>

#include "kutil.h"
#include "pev.h"

namespace pev {

double svar(const double *u, int n, const int *st, int nst, const int *ps, int nps,
            const double *fpc) {
  // the PSU totals, and each stratum's PSUs and their mean total
  std::vector<double> tot(static_cast<std::size_t>(nps), 0.0), sm(static_cast<std::size_t>(nst), 0.0),
      nh(static_cast<std::size_t>(nst), 0.0);
  std::vector<int> sp(static_cast<std::size_t>(nps), -1);
  for (int i = 0; i < n; ++i) {
    tot[static_cast<std::size_t>(ps[i])] += u[i];
    sp[static_cast<std::size_t>(ps[i])] = st[i];
  }
  for (int a = 0; a < nps; ++a) {
    sm[static_cast<std::size_t>(sp[static_cast<std::size_t>(a)])] += tot[static_cast<std::size_t>(a)];
    nh[static_cast<std::size_t>(sp[static_cast<std::size_t>(a)])] += 1.0;
  }
  kutil::ksum v;
  std::vector<double> ss(static_cast<std::size_t>(nst), 0.0);
  for (int a = 0; a < nps; ++a) {
    const std::size_t h = static_cast<std::size_t>(sp[static_cast<std::size_t>(a)]);
    const double dlt = tot[static_cast<std::size_t>(a)] - sm[h] / nh[h];
    ss[h] += dlt * dlt;
  }
  for (int h = 0; h < nst; ++h) {
    const std::size_t q = static_cast<std::size_t>(h);
    const double f = fpc ? nh[q] / fpc[q] : 0.0;
    v.add((1.0 - f) * nh[q] / (nh[q] - 1.0) * ss[q]);
  }
  return v.val();
}

}  // namespace pev
