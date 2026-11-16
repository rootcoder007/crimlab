// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/gol.cpp: without covariates the intercepts are the logits
// of the shares above each threshold, the start itself; a category with
// no rows has no fit. Built with the self-test.
#include <cmath>

#include "crimlab.h"
#include "tchk.h"

void t_gol(tk &t, int big) {
  (void)big;
  const int y[4] = {0, 0, 1, 2};
  double est[2], cov[4], st[3];
  tk_ieq(t, "gol.rc", clab_gologt(y, 4, 3, nullptr, 0, nullptr, 100, est, cov, nullptr, st), CLAB_OK);
  tk_near(t, "gol.a0", est[0], 0, 1e-8);
  tk_near(t, "gol.a1", est[1], std::log(1.0 / 3), 1e-8);
  tk_near(t, "gol.ll", st[0], 2 * std::log(0.5) + 2 * std::log(0.25), 1e-12);
  tk_eq(t, "gol.conv", st[2], 1, 0);
  // the information of the first intercept alone: n F (1 - F) = 1 at F =
  // 1/2, less what the second shares with it
  tk_true(t, "gol.se", std::isfinite(cov[0]) && cov[0] > 0);
  const int z[4] = {0, 0, 2, 2};
  tk_ieq(t, "gol.empty", clab_gologt(z, 4, 3, nullptr, 0, nullptr, 100, est, cov, nullptr, st), CLAB_EDOM);
}
