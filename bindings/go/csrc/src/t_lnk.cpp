// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/lnk.cpp: the pair measures and the areas by hand.
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "tchk.h"

void t_lnk(tk &T, int big) {
  (void)big;
  // two series of two crimes and one unsolved crime; behaviours in two
  // columns
  const std::vector<double> x = {0, 3, 10, 10, 50}, y = {0, 4, 0, 1, 50}, t = {0, 2, 0, 9, 1};
  const std::vector<int> mo = {1, 1, 0, 0, 1, 0, 1, 1, 1, 0}, ser = {0, 0, 1, 1, -1};
  std::vector<int> pi(10), pj(10), lk(10);
  std::vector<double> d(10), dt(10), j(10), pr(10);
  double cf[4], se[4], au[4], dev = 0;
  const int rc = clab_lnkpr(x.data(), y.data(), t.data(), mo.data(), 2, ser.data(), 5, pi.data(), pj.data(), d.data(),
                            dt.data(), j.data(), lk.data(), pr.data(), cf, se, au, &dev);
  // four linked-or-unlinked solved pairs this small separate perfectly
  tk_ieq(T, "lk.sep", rc, CLAB_EDOM);
  tk_eq(T, "lk.d", d[0], 5, 1e-15);
  tk_eq(T, "lk.dt", dt[0], 2, 0);
  // behaviours (column-major): crime 0 {a}, crime 1 {a, b}, crime 2 {b}:
  // Jaccard 1/2 for crimes 0 and 1, 0 for crimes 0 and 2
  tk_eq(T, "lk.j0", j[0], 0.5, 0);
  tk_eq(T, "lk.j1", j[1], 0, 0);
  tk_ieq(T, "lk.lk", lk[0] * 100 + lk[1] * 10 + (lk[3] + 1), 100);
  tk_ieq(T, "lk.pij", pi[9] * 10 + pj[9], 34);
  const std::vector<int> one = {0, 0, 0, 0, 0};
  tk_ieq(T, "lk.none",
         clab_lnkpr(x.data(), y.data(), t.data(), nullptr, 0, one.data(), 5, pi.data(), pj.data(), d.data(), dt.data(),
                    j.data(), lk.data(), pr.data(), cf, se, au, &dev),
         CLAB_EDOM);
}
