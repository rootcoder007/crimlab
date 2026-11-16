// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Association rules between two categorical variables (crimlab.h,
// asrule()): the rules that the Apriori algorithm of Agrawal et al.
// (1996) finds when each transaction holds one level of each of two
// variables, as Wongchanapai et al. (2025) use it for binned exposure and
// offence rates: {} => {b} and {a} => {b} in both directions, with their
// support, confidence, coverage and lift, kept at a minimum support and
// confidence; a pair never observed together is no rule. Work: O(n +
// levels^2).
#include <cmath>
#include <cstddef>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

namespace {

inline std::size_t z(int i) { return static_cast<std::size_t>(i); }

}  // namespace

extern "C" int clab_asrule(int n, const int *xc, int nx, const int *yc, int ny, double minsup,
                           double minconf, int empty, int *nr, int *lhs, int *rhs,
                           double *out) CLAB_TRY {
  if (n < 1 || nx < 1 || ny < 1 || !(minsup >= 0.0 && minsup <= 1.0) ||
      !(minconf >= 0.0 && minconf <= 1.0)) {
    return CLAB_EARG;
  }
  if (static_cast<double>(nx) * ny > kutil::maxcell) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (xc[i] < 0 || xc[i] >= nx || yc[i] < 0 || yc[i] >= ny) return CLAB_EARG;
  }
  // items: x levels 0..nx-1, y levels nx..nx+ny-1; -1 the empty set
  const int ni = nx + ny;
  std::vector<double> cnt(z(ni), 0.0), joint(z(nx) * z(ny), 0.0);
  for (int i = 0; i < n; ++i) {
    cnt[z(xc[i])] += 1.0;
    cnt[z(nx + yc[i])] += 1.0;
    joint[z(yc[i]) * z(nx) + z(xc[i])] += 1.0;
  }
  const double dn = n;
  int r = 0;
  // out per rule: support, confidence, coverage, lift, count
  auto add = [&](int a, int b, double both, double ca) {
    const double sup = both / dn, conf = ca > 0.0 ? both / ca : NAN, cov = ca / dn;
    if (!(both > 0.0) || !(sup >= minsup) || !(conf >= minconf)) return;
    lhs[r] = a;
    rhs[r] = b;
    double *o = out + z(r) * 5;
    o[0] = sup;
    o[1] = conf;
    o[2] = cov;
    o[3] = conf / (cnt[z(b)] / dn);
    o[4] = both;
    ++r;
  };
  if (empty) {
    for (int b = 0; b < ni; ++b) add(-1, b, cnt[z(b)], dn);
  }
  for (int a = 0; a < nx; ++a) {
    for (int b = 0; b < ny; ++b) {
      const double j = joint[z(b) * z(nx) + z(a)];
      add(a, nx + b, j, cnt[z(a)]);
      add(nx + b, a, j, cnt[z(nx + b)]);
    }
  }
  *nr = r;
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)
