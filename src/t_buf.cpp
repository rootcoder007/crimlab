// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/buf.cpp: the buffer groups and distances by hand.
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "tchk.h"

namespace {

using vd = std::vector<double>;
using vi = std::vector<int>;

void famtest(tk &T, int big) {
  (void)big;
  // buffers: within r of a site; within r of a comparison point beyond buf
  // of every site; the rest
  {
    const vd x = {0.5, 10, 10, 4, 50}, y = {0, 0.5, 4.5, 0, 0};
    const double sx = 0, sy = 0, cx = 10, cy = 0;
    vi g(5);
    vd ds(5), dc(5);
    tk_ieq(T, "bf.ok",
           clab_bufas(x.data(), y.data(), 5, &sx, &sy, 1, &cx, &cy, 1, 1, 3, g.data(), ds.data(),
                      dc.data()),
           CLAB_OK);
    const vi wg = {0, 1, 2, 2, 2};
    bool ok = true;
    for (int i = 0; i < 5; ++i) ok = ok && g[i] == wg[i];
    tk_true(T, "bf.hand", ok);
    tk_eq(T, "bf.dc", dc[2], 4.5, 0);
    const double inf = HUGE_VAL;
    tk_ieq(T, "bf.inf",
           clab_bufas(x.data(), y.data(), 5, &inf, &sy, 1, &cx, &cy, 1, 1, 3, g.data(), ds.data(),
                      dc.data()),
           CLAB_EDOM);
    tk_ieq(T, "bf.r",
           clab_bufas(x.data(), y.data(), 5, &sx, &sy, 1, &cx, &cy, 1, 0, 3, g.data(), ds.data(),
                      dc.data()),
           CLAB_EDOM);
  }
}

}  // namespace

void t_buf(tk &T, int big) { famtest(T, big); }
