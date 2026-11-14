// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/prh.cpp: a two-cell grid by hand.
#include <vector>

#include "crimlab.h"
#include "tchk.h"

void t_prh(tk &T, int big) {
  (void)big;
  // one past incident at the centre of cell 0, a week before; two later
  // ones, in cell 0 and cell 1
  const std::vector<double> x = {0.5, 0.6, 1.6, 0.5}, y = {0.5, 0.4, 0.5, 0.5}, t = {0, 8, 9, -100};
  const double w[4] = {0, 2, 0, 1};
  double z[2], ev[5];
  int hot[2];
  tk_ieq(T, "ph.ok", clab_prohs(x.data(), y.data(), t.data(), 4, 7, w, 1, 2, 60, 7, 0.5, 7, 2, 1, z, hot, ev), CLAB_OK);
  tk_eq(T, "ph.z0", z[0], 0.5, 1e-15);
  tk_eq(T, "ph.z1", z[1], 0.25, 1e-15);
  tk_ieq(T, "ph.hot", hot[0] * 2 + hot[1], 2);
  tk_eq(T, "ph.n", ev[0], 2, 0);
  tk_eq(T, "ph.hit", ev[2], 0.5, 0);
  tk_eq(T, "ph.pai", ev[4], 1, 1e-15);
  // beyond the bandwidth nothing reaches cell 1
  tk_ieq(T, "ph.bw", clab_prohs(x.data(), y.data(), t.data(), 4, 7, w, 1, 0.9, 60, 7, 1, 7, 2, 1, z, hot, ev), CLAB_OK);
  tk_eq(T, "ph.bz", z[1], 0, 0);
  tk_ieq(T, "ph.bh", hot[1], 0);
  tk_ieq(T, "ph.dom", clab_prohs(x.data(), y.data(), t.data(), 4, 7, w, 0, 2, 60, 7, 0.5, 7, 2, 1, z, hot, ev), CLAB_EDOM);
}
