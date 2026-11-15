// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, displacement and diffusion of benefits: the weighted
// displacement quotient of Bowers and Johnson (2003) and the gross, net
// and total net effects as Guerette (2009, Appendix B) states them.
#include <cmath>

#include "crimlab.h"
#include "kutil.h"

int clab_wdqx(const double *a, const double *b, const double *c, double *out) {
  for (int k = 0; k < 2; ++k) {
    if (!kutil::nonneg(a[k]) || !kutil::nonneg(b[k]) || !kutil::pos(c[k])) return CLAB_EDOM;
  }
  // the control's counts as divisors: positive here, so c0 and c1 are
  // c[0] and c[1]; the selection keeps a divisor that a compiler can see
  // to be 0 (a caller's constant, after inlining) out of every division
  const double c0 = c[0] > 0 ? c[0] : 1.0, c1 = c[1] > 0 ? c[1] : 1.0;
  // the change in the treated area and in the buffer, each relative to the
  // control area's
  const double sm = a[1] / c1 - a[0] / c0, bdm = b[1] / c1 - b[0] / c0;
  const double g = c1 / c0;
  out[0] = sm;
  out[1] = bdm;
  // the quotient: NaN when sm is 0, the divisor 1 there for the same reason
  const double sd = sm != 0 ? sm : 1.0;
  out[2] = sm != 0 ? bdm / sd : kutil::nan();
  // the gross and net effects (Guerette 2009, equations 1 and 2)
  out[3] = a[0] - a[1];
  out[4] = a[0] / c0 - a[1] / c1;
  // the total net effect (equation 4) and the counts expected after in
  // the treated area and the buffer that it compares with the counts seen
  out[5] = a[0] * g;
  out[6] = b[0] * g;
  out[7] = (out[5] - a[1]) + (out[6] - b[1]);
  return CLAB_OK;
}
