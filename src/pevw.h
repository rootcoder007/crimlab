// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The wild cluster bootstrap and the stepdown adjustment of wcb.cpp
// (wcbx(), wystep()).
#ifndef CLAB_PEVW_H
#define CLAB_PEVW_H

#include <vector>

namespace pev {

// A regression already within-transformed: y (n) on the p columns of x
// (n x p, the tested coefficient first), clusters g (n codes, ng of them),
// the unit and period codes u and t (nu, nt of them) of the two-way
// effects, and the CR1 factor adj. wsetup() fills the rest; false when
// X'X is singular.
struct wctx {
  int n = 0, p = 0, ng = 0, nu = 0, nt = 0;
  std::vector<double> x, y;
  std::vector<int> g, u, t;
  double adj = 1.0;
  std::vector<double> xi, rxi, h;
};
bool wsetup(wctx &c);

// The cluster-robust t statistic of the first coefficient against b0 for
// the outcome y; est (or NULL) gets the coefficient.
double tstat(const wctx &c, const double *y, double b0, double *est);

// The restricted wild cluster bootstrap p-value of b0: the share of the
// nd draws (rows of v, one weight per cluster) whose |t*| is at least |t|
// (ties 1) or above it (ties 0). Each draw's reweighted residuals are
// within-transformed again, as refitting the effects would.
double wpval(const wctx &c, double b0, const std::vector<double> &v, int nd, int ties);

// One end of the interval (dir -1 lower, 1 upper): the b0 where the
// p-value crosses alpha, by doubling steps of s0 then bisection; infinite
// when no b0 within 2^12 s0 is rejected.
double wend(const wctx &c, double bh, double s0, double alpha, int dir,
            const std::vector<double> &v, int nd, int ties);

// Westfall-Young stepdown (algorithm 2.8): obs (k), null (nd x k, column
// per hypothesis). Out: p (k), each hypothesis's share of null |draws| at
// least |obs|; pa (k), the adjusted p-values, monotone in |obs|.
void wystep(const double *obs, const double *null, int nd, int k, double *p, double *pa);

}  // namespace pev

#endif
