// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/geou.cpp: point in polygon with a hole, an edge, overlap and
// a polygon of two parts; the crosswalk by hand.
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "tchk.h"

namespace {

using vd = std::vector<double>;
using vi = std::vector<int>;

void famtest(tk &T, int big) {
  (void)big;
  // polygon 0: the square [0, 4]^2 less the hole [1, 2]^2; polygon 1: the
  // square [3, 5] x [0, 1] and, a second part, the triangle (6,0) (8,0)
  // (7,2)
  {
    const vd vx = {0, 4, 4, 0, 1, 2, 2, 1, 3, 5, 5, 3, 6, 8, 7};
    const vd vy = {0, 0, 4, 4, 1, 1, 2, 2, 0, 0, 1, 1, 0, 0, 2};
    const vi pid = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1};
    const vi rid = {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1};
    const vd px = {0.5, 1.5, 4, 3.5, 7, 9, 2, 7, 4.5};
    const vd py = {0.5, 1.5, 2, 0.5, 1, 1, 1.5, 1.99, 0.5};
    vi reg(9), nr(9);
    tk_ieq(T, "pip.ok",
           clab_pinpl(px.data(), py.data(), 9, vx.data(), vy.data(), pid.data(), rid.data(), 15, 2,
                      reg.data(), nr.data()),
           CLAB_OK);
    // inside; in the hole; on the outer edge; in both; in the second part;
    // outside; on the hole's edge; inside near the apex; in polygon 1 only
    const vi wr = {0, -1, 0, 0, 1, -1, 0, 1, 1}, wn = {1, 0, 1, 2, 1, 0, 1, 1, 1};
    bool ok = true;
    for (int i = 0; i < 9; ++i) ok = ok && reg[i] == wr[i] && nr[i] == wn[i];
    tk_true(T, "pip.hand", ok);
    const vi bad = {1, 0, 0};
    tk_ieq(T, "pip.order",
           clab_pinpl(px.data(), py.data(), 1, vx.data(), vy.data(), bad.data(), rid.data(), 3, 2,
                      reg.data(), nr.data()),
           CLAB_EARG);
  }
  // crosswalk: a split by area and an aggregation
  {
    const vd val = {100, 40, 60};
    const vi from = {0, 0, 1, 2}, to = {0, 1, 1, 1};
    const vd area = {3, 1, 5, 2};
    vd out(2), sh(4);
    vi mn(3);
    double lost = -1;
    tk_ieq(T, "aw.ok",
           clab_arwt(val.data(), 3, from.data(), to.data(), area.data(), 4, 2, out.data(),
                     sh.data(), mn.data(), &lost),
           CLAB_OK);
    const vd wo = {75, 25 + 40 + 60};
    tk_veq(T, "aw.out", out.data(), wo.data(), 2, 1e-15);
    tk_true(T, "aw.main", mn[0] == 0 && mn[1] == 1 && mn[2] == 1 && lost == 0);
    // pairs without area split evenly; a source without pairs is lost
    const vi f2 = {0, 0}, t2 = {0, 1};
    const vd a0 = {0, 0};
    clab_arwt(val.data(), 3, f2.data(), t2.data(), a0.data(), 2, 2, out.data(), sh.data(),
              mn.data(), &lost);
    const vd w2 = {50, 50};
    tk_veq(T, "aw.even", out.data(), w2.data(), 2, 0);
    tk_true(T, "aw.lost", lost == 100 && mn[0] == 0 && mn[1] == -1);
    const vd neg = {-1, 1};
    tk_ieq(T, "aw.neg",
           clab_arwt(val.data(), 3, f2.data(), t2.data(), neg.data(), 2, 2, out.data(), sh.data(),
                     mn.data(), &lost),
           CLAB_EDOM);
  }
}

}  // namespace

void t_geou(tk &T, int big) { famtest(T, big); }
