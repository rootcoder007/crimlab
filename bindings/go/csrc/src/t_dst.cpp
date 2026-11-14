// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for the distributions (src/dst.cpp) against values of
// R's qnorm() and pnorm(). Part of the self-test (selftest.cpp).
#include <cmath>

#include "crimlab.h"
#include "tchk.h"

void t_dst(tk &t, int big) {
  (void)big;
  tk_eq(t, "qnorm.975", clab_qnorm(0.975), 1.959963984540054, 1e-15);
  tk_eq(t, "qnorm.3", clab_qnorm(0.3), -0.52440051270804067, 1e-15);
  tk_eq(t, "qnorm.995", clab_qnorm(0.995), 2.5758293035488999, 1e-15);
  tk_eq(t, "qnorm.tail", clab_qnorm(1e-10), -6.3613409024040557, 1e-15);
  tk_eq(t, "qnorm.far", clab_qnorm(1e-300), -37.04709629936120052, 1e-15);
  tk_eq(t, "qnorm.half", clab_qnorm(0.5), 0, 0);
  tk_true(t, "qnorm.ends", clab_qnorm(0) == -HUGE_VAL && clab_qnorm(1) == HUGE_VAL);
  tk_nan(t, "qnorm.out", clab_qnorm(1.5));
  tk_nan(t, "qnorm.nan", clab_qnorm(NAN));
  tk_eq(t, "pnorm.196", clab_pnorm(1.959963984540054), 0.975, 1e-15);
  tk_eq(t, "pnorm.tail", clab_pnorm(-10), 7.6198530241605269e-24, 1e-13);
}
