// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/mtr.cpp: three sentences by hand.
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "tchk.h"

void t_mtr(tk &T, int big) {
  (void)big;
  // weighted shares 0.5, 0.3, 0.2; means 0.2, 0.4, 0.6
  const std::vector<double> y = {1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0};
  const std::vector<int> z = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2};
  const std::vector<double> w = {5, 5, 5, 5, 5, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2};
  // the middle sentence: below it share 0.5 (mean 0.2), above it 0.2
  // (mean 0.6); none, MTR, MTS, both
  const double want[4][2] = {{0.12, 0.82}, {0.1 + 0.12, 0.5 + 0.12 + 0.12}, {0.4 * 0.5, 0.4 * 0.8 + 0.2},
                             {0.1 + 0.4 * 0.5, 0.4 * 0.8 + 0.12}};
  double o[18];
  for (int as = 0; as < 4; ++as) {
    tk_ieq(T, "mt.ok", clab_mtrbnd(y.data(), z.data(), w.data(), 15, 3, 0, 1, as, o), CLAB_OK);
    tk_eq(T, "mt.lo", o[6 + 1], want[as][0], 1e-14);
    tk_eq(T, "mt.hi", o[9 + 1], want[as][1], 1e-14);
  }
  tk_eq(T, "mt.share", o[1], 0.3, 1e-15);
  tk_eq(T, "mt.mean", o[3 + 2], 0.6, 1e-15);
  tk_nan(T, "mt.ate0", o[12]);
  tk_ge(T, "mt.ate", o[12 + 1], 0);
  const std::vector<double> y2 = {2, 0};
  const std::vector<int> z2 = {0, 1};
  tk_ieq(T, "mt.dom", clab_mtrbnd(y2.data(), z2.data(), nullptr, 2, 2, 0, 1, 0, o), CLAB_EDOM);
}
