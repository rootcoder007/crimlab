// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/lnr.cpp: an eight-witness lineup study and a mock-witness
// lineup worked by hand.
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "tchk.h"

namespace {

using vd = std::vector<double>;
using vi = std::vector<int>;

void roc(tk &T) {
  // target present: suspect at 3, suspect at 1, filler at 2, none;
  // target absent: suspect at 2, filler at 3, none, none
  const vi tp = {1, 1, 1, 1, 0, 0, 0, 0}, pk = {1, 1, 2, 0, 1, 2, 0, 0}, pr(8, 0);
  const vd cf = {3, 1, 2, NAN, 2, 3, NAN, NAN};
  vd lev(8), r(24), c(24), pa(1), bs(2);
  int nl = 0;
  double fu = 0;
  tk_ieq(T, "ln.ok",
         clab_lnroc(tp.data(), pk.data(), cf.data(), pr.data(), 8, 1, 0, 0, 0, 0, lev.data(), &nl, r.data(),
                    c.data(), pa.data(), bs.data(), &fu),
         CLAB_OK);
  tk_ieq(T, "ln.nl", nl, 3);
  tk_eq(T, "ln.lev", lev[0], 3, 0);
  tk_eq(T, "ln.hr", r[2], 0.5, 1e-15);
  tk_eq(T, "ln.far", r[3 + 1], 0.25, 1e-15);
  tk_true(T, "ln.drinf", std::isinf(r[6]));
  tk_eq(T, "ln.dr", r[6 + 2], 2, 1e-15);
  tk_eq(T, "ln.pauc", pa[0], 0.0625, 1e-15);
  tk_eq(T, "ln.fused", fu, 0.25, 0);
  tk_eq(T, "ln.cac", c[6 + 1], 0, 0);
  tk_eq(T, "ln.cac1", c[6 + 2], 1, 0);
  // no designated innocent suspect: target-absent picks over the size
  tk_ieq(T, "ln.size",
         clab_lnroc(tp.data(), pk.data(), cf.data(), pr.data(), 8, 1, 3, 0.1, 0, 0, lev.data(), &nl, r.data(),
                    c.data(), pa.data(), bs.data(), &fu),
         CLAB_OK);
  tk_eq(T, "ln.sfar", r[3], 1.0 / 12, 1e-15);
  tk_eq(T, "ln.sfar2", r[4], 2.0 / 12, 1e-15);
  // up to 0.1: the segment from (0, 0) to (1/12, 0.25), then flat at 0.25
  // towards (1/6, 0.25): 1/96 + 1/240
  tk_eq(T, "ln.spauc", pa[0], 7.0 / 480.0, 1e-14);
  // the bootstrap gives a standard error
  tk_ieq(T, "ln.boot",
         clab_lnroc(tp.data(), pk.data(), cf.data(), pr.data(), 8, 1, 0, 0, 50, 7, lev.data(), &nl, r.data(),
                    c.data(), pa.data(), bs.data(), &fu),
         CLAB_OK);
  tk_ge(T, "ln.bse", bs[0], 0);
  const vi pk0(8, 0);
  tk_ieq(T, "ln.none",
         clab_lnroc(tp.data(), pk0.data(), cf.data(), pr.data(), 8, 1, 0, 0, 0, 0, lev.data(), &nl, r.data(),
                    c.data(), pa.data(), bs.data(), &fu),
         CLAB_EDOM);
}

void size(tk &T) {
  const vi ch = {0, 0, 0, 1, 2, 2};
  vd cnt(3), o(11);
  tk_ieq(T, "ls.ok", clab_lnsize(ch.data(), 6, 3, 0, 1.96, cnt.data(), o.data()), CLAB_OK);
  tk_eq(T, "ls.cnt", cnt[2], 2, 0);
  tk_eq(T, "ls.ps", o[3], 0.5, 1e-15);
  tk_eq(T, "ls.pbin", o[4], 233.0 / 729.0, 1e-13);
  tk_eq(T, "ls.fs", o[5], 2, 1e-15);
  tk_eq(T, "ls.mal", o[6], 2.5, 1e-15);
  tk_eq(T, "ls.malb", o[7], 2.5, 1e-15);
  tk_eq(T, "ls.etr", o[8], 36.0 / 14.0, 1e-14);
  tk_lt(T, "ls.lo", o[9], o[8]);
  // Tredoux (1996), Table 4.2, lineup g: 0 25 5 25 3 2, adjusted 2.83,
  // unadjusted 3.00
  vi g;
  const int gc[6] = {0, 25, 5, 25, 3, 2};
  for (int j = 0; j < 6; ++j) {
    for (int r = 0; r < gc[j]; ++r) g.push_back(j);
  }
  vd c6(6), o6(11);
  tk_ieq(T, "ls.gok", clab_lnsize(g.data(), 60, 6, 3, 1.96, c6.data(), o6.data()), CLAB_OK);
  tk_eq(T, "ls.ga", o6[6], 5 - 52.0 / 24, 1e-14);
  tk_eq(T, "ls.gb", o6[7], 3, 1e-14);
  // no choice (-1) counts for the share, not for the spread
  const vi nc = {0, 0, -1, 1};
  tk_ieq(T, "ls.nok", clab_lnsize(nc.data(), 4, 2, 0, 1.96, cnt.data(), o.data()), CLAB_OK);
  tk_eq(T, "ls.nps", o[3], 0.5, 1e-15);
  tk_eq(T, "ls.nm", o[2], 3, 0);
  tk_eq(T, "ls.nfs", o[5], 2, 1e-15);
  const vi none = {-1, -1};
  tk_ieq(T, "ls.nnone", clab_lnsize(none.data(), 2, 2, 0, 1.96, cnt.data(), o.data()), CLAB_EDOM);
  const vi bad = {0, 3};
  tk_ieq(T, "ls.dom", clab_lnsize(bad.data(), 2, 3, 0, 1.96, cnt.data(), o.data()), CLAB_EDOM);
}

}  // namespace

void t_lnr(tk &T, int big) {
  (void)big;
  roc(T);
  size(T);
}
