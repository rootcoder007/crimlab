// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/spc.cpp: with the censored followed far beyond every
// failure the exponential cure model has a closed form, the share failing
// and the mean time of the failures.
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "tchk.h"

void t_spc(tk &T, int big) {
  (void)big;
  const std::vector<double> t = {1, 2, 3, 4, 1e6, 1e6, 1e6, 1e6}, ev = {1, 1, 1, 1, 0, 0, 0, 0};
  double par[5], cov[25], st[5];
  tk_ieq(T, "sc.ok", clab_spcure(t.data(), ev.data(), nullptr, 8, 0, 1, 500, par, cov, st), CLAB_OK);
  tk_near(T, "sc.g", par[0], 0, 1e-6);
  tk_near(T, "sc.e", par[1], std::log(2.5), 1e-6);
  tk_eq(T, "sc.conv", st[3], 1, 0);
  tk_near(T, "sc.pf", st[4], 0.5, 1e-6);
  // the information of the logit share: n delta (1 - delta) = 2
  tk_near(T, "sc.se", std::sqrt(cov[0]), std::sqrt(0.5), 1e-4);
  // Weibull and lognormal fits converge on the same data
  tk_ieq(T, "sc.w", clab_spcure(t.data(), ev.data(), nullptr, 8, 0, 0, 500, par, cov, st), CLAB_OK);
  tk_near(T, "sc.wg", par[0], 0, 1e-5);
  tk_ieq(T, "sc.ln", clab_spcure(t.data(), ev.data(), nullptr, 8, 0, 2, 500, par, cov, st), CLAB_OK);
  tk_near(T, "sc.lg", par[0], 0, 1e-5);
  tk_near(T, "sc.lm", par[1], (std::log(2.0) + std::log(3.0) + std::log(4.0)) / 4, 1e-5);
  const std::vector<double> e1(8, 1.0);
  tk_ieq(T, "sc.dom", clab_spcure(t.data(), e1.data(), nullptr, 8, 0, 1, 500, par, cov, st), CLAB_EDOM);
}
