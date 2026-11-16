// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for psc.cpp: the credit at one and at one and a half
// days a day by hand, the sentence left, and the failure paths. Part of
// the self-test (selftest.cpp).
#include <cmath>

#include "crimlab.h"
#include "tchk.h"

namespace {

void psc(tk &tt) {
  // 10 days at 1.5 is 15, 45 days at 1.5 is 67.5; a sentence of 60 days
  // less 15 leaves 45, less 67.5 leaves nothing
  const double d[2] = {10, 45}, r[1] = {1.5}, s[1] = {60};
  double c[2], m[2];
  tk_ieq(tt, "psc.rc", clab_pscred(d, 2, r, 1, s, 1, c, m), CLAB_OK);
  tk_near(tt, "psc.c0", c[0], 15, 0);
  tk_near(tt, "psc.c1", c[1], 67.5, 0);
  tk_near(tt, "psc.m0", m[0], 45, 0);
  tk_near(tt, "psc.m1", m[1], 0, 0);
  const double r2[2] = {1, 1.5};
  clab_pscred(d, 2, r2, 2, nullptr, 0, c, m);
  tk_near(tt, "psc.r1", c[0], 10, 0);
  const double big[1] = {1.6}, neg[1] = {-1};
  tk_ieq(tt, "psc.ratio", clab_pscred(d, 2, big, 1, nullptr, 0, c, m), CLAB_EDOM);
  tk_ieq(tt, "psc.days", clab_pscred(neg, 1, r, 1, nullptr, 0, c, m), CLAB_EDOM);
  tk_ieq(tt, "psc.len", clab_pscred(d, 2, r2, 2, s, 3, c, m), CLAB_EARG);
}

}  // namespace

void t_ovr(tk &t, int big) {
  (void)big;
  psc(t);
}
