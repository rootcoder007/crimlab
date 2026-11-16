// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/hsa.cpp: ten cells by hand.
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "tchk.h"

void t_hsa(tk &T, int big) {
  (void)big;
  const std::vector<double> pr = {9, 8, 7, 6, 5, 4, 3, 2, 1, 0}, ob = {3, 0, 5, 1, 0, 0, 2, 9, 0, 0},
                            ar(10, 1.0), pv = {4, 2, 2, 0, 0, 0, 1, 1, 0, 0};
  const double top[2] = {0.2, 0.3};
  std::vector<double> r(22);
  tk_ieq(T, "ha.ok", clab_hsacc(pr.data(), nullptr, ob.data(), ar.data(), pv.data(), 10, top, 2, 0, r.data()), CLAB_OK);
  // top 20%: cells 0 and 1, 3 of 20 crimes on 2 of 10 cells
  tk_eq(T, "ha.cells", r[2], 2, 0);
  tk_eq(T, "ha.hit", r[12], 0.15, 1e-15);
  tk_eq(T, "ha.pai", r[14], 0.75, 1e-14);
  // the best two cells hold 9 + 5
  tk_eq(T, "ha.nmax", r[16], 14, 0);
  tk_eq(T, "ha.pei", r[18], 3.0 / 14, 1e-15);
  // recapture: 3 now of 6 before, against 20 of 10
  tk_eq(T, "ha.rri", r[20], (3.0 / 6) / (20.0 / 10), 1e-15);
  // top 30%: cells 0 to 2, 8 crimes
  tk_eq(T, "ha.n3", r[9], 8, 0);
  const std::vector<int> h = {0, 0, 0, 0, 0, 0, 0, 1, 0, 0};
  tk_ieq(T, "ha.lgl", clab_hsacc(nullptr, h.data(), ob.data(), ar.data(), nullptr, 10, top, 1, 0, r.data()), CLAB_OK);
  tk_eq(T, "ha.lpei", r[9], 1, 0);
  tk_nan(T, "ha.nrri", r[10]);
  const std::vector<double> z(10, 0.0);
  tk_ieq(T, "ha.dom", clab_hsacc(pr.data(), nullptr, z.data(), ar.data(), nullptr, 10, top, 1, 0, r.data()), CLAB_EDOM);
  tk_ieq(T, "ha.tie", clab_hsacc(pr.data(), nullptr, ob.data(), ar.data(), nullptr, 10, top, 1, 3, r.data()), CLAB_EARG);
  // four cells tie at the top, the cut-off (20%) takes two of their area:
  // all four, none, or each by half
  const std::vector<double> pt = {5, 5, 5, 5, 4, 3, 2, 1, 0, 0};
  const double t2[1] = {0.2};
  for (int k = 0; k < 3; ++k) {
    tk_ieq(T, "ha.tk", clab_hsacc(pt.data(), nullptr, ob.data(), ar.data(), nullptr, 10, t2, 1, k, r.data()), CLAB_OK);
    tk_eq(T, "ha.tcl", r[1], k == 0 ? 4 : k == 1 ? 0 : 2, 1e-15);
    tk_eq(T, "ha.tn", r[4], k == 0 ? 9 : k == 1 ? 0 : 4.5, 1e-15);
  }
}
