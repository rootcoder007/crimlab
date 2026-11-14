// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, displacement and diffusion of benefits: the weighted
// displacement quotient of Bowers and Johnson (2003).
#include <cmath>

#include "crimlab.h"
#include "kutil.h"

int clab_wdqx(const double *a, const double *b, const double *c, double *out) {
  for (int k = 0; k < 2; ++k) {
    if (!kutil::nonneg(a[k]) || !kutil::nonneg(b[k]) || !kutil::pos(c[k])) return CLAB_EDOM;
  }
  // the change in the treated area and in the buffer, each relative to the
  // control area's
  const double sm = a[1] / c[1] - a[0] / c[0], bdm = b[1] / c[1] - b[0] / c[0];
  const double g = c[1] / c[0];
  out[0] = sm;
  out[1] = bdm;
  out[2] = sm != 0 ? bdm / sm : kutil::nan();
  out[3] = a[0] * g;
  out[4] = b[0] * g;
  out[5] = (out[3] - a[1]) + (out[4] - b[1]);
  return CLAB_OK;
}
