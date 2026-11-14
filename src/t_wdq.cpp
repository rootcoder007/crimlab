// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/wdq.cpp: the quotient by hand.
#include "crimlab.h"
#include "tchk.h"

void t_wdq(tk &T, int big) {
  (void)big;
  // treated 100 -> 60, buffer 50 -> 45, control 200 -> 180
  const double a[2] = {100, 60}, b[2] = {50, 45}, c[2] = {200, 180};
  double o[6];
  tk_ieq(T, "wd.ok", clab_wdqx(a, b, c, o), CLAB_OK);
  tk_eq(T, "wd.sm", o[0], 60.0 / 180 - 0.5, 1e-15);
  tk_eq(T, "wd.bdm", o[1], 0.25 - 0.25, 1e-15);
  tk_eq(T, "wd.wdq", o[2], 0, 1e-14);
  tk_eq(T, "wd.ea", o[3], 90, 1e-13);
  tk_eq(T, "wd.net", o[5], 30, 1e-13);
  const double b2[2] = {50, 60};
  tk_ieq(T, "wd.ok2", clab_wdqx(a, b2, c, o), CLAB_OK);
  tk_eq(T, "wd.disp", o[2], (60.0 / 180 - 0.25) / (60.0 / 180 - 0.5), 1e-14);
  const double c0[2] = {0, 1};
  tk_ieq(T, "wd.dom", clab_wdqx(a, b, c0, o), CLAB_EDOM);
  const double a0[2] = {10, 10}, c1[2] = {10, 10};
  tk_ieq(T, "wd.zero", clab_wdqx(a0, b, c1, o), CLAB_OK);
  tk_nan(T, "wd.nan", o[2]);
}
