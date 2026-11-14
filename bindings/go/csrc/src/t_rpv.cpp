// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/rpv.cpp: a small record of incidents by hand.
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "tchk.h"

namespace {

using vd = std::vector<double>;
using vi = std::vector<int>;

void famtest(tk &T) {
  // target 0 three times (days 1, 4, 20), target 1 twice (5, 6), target 2
  // once
  const vi id = {0, 1, 0, 2, 0, 1};
  const vd t = {1, 5, 4, 9, 20, 6}, brk = {0, 7, 30};
  vd tab(21), st(10), gap(3), lag(6);
  int K = 0;
  tk_ieq(T, "rv.ok",
         clab_rptvic(id.data(), t.data(), 6, 3, 10, brk.data(), 3, tab.data(), &K, st.data(), gap.data(), lag.data()),
         CLAB_OK);
  tk_ieq(T, "rv.K", K, 3);
  tk_eq(T, "rv.k0", tab[0], 7, 0);
  tk_eq(T, "rv.k2", tab[2], 1, 0);
  tk_eq(T, "rv.inc3", tab[4 + 3], 3, 0);
  tk_eq(T, "rv.exp1", tab[8 + 1], 10 * 0.6 * std::exp(-0.6), 1e-14);
  tk_eq(T, "rv.rep", st[2], 3, 0);
  tk_eq(T, "rv.conc", st[4], 2, 0);
  tk_eq(T, "rv.prev", st[6], 0.3, 1e-15);
  tk_eq(T, "rv.erep", st[8], 10 * (1 - std::exp(-0.6) * 1.6), 1e-13);
  // gaps 3, 16 and 1
  tk_eq(T, "rv.gn", gap[0], 3, 0);
  tk_eq(T, "rv.gmed", gap[2], 3, 0);
  tk_eq(T, "rv.gmean", gap[1], 20.0 / 3, 1e-14);
  tk_eq(T, "rv.lag0", lag[0], 2, 0);
  tk_eq(T, "rv.cum", lag[4 + 1], 1, 0);
  tk_ieq(T, "rv.pop",
         clab_rptvic(id.data(), nullptr, 6, 3, 2, nullptr, 0, tab.data(), &K, st.data(), gap.data(), lag.data()),
         CLAB_EDOM);
  tk_ieq(T, "rv.nopop",
         clab_rptvic(id.data(), nullptr, 6, 3, NAN, nullptr, 0, tab.data(), &K, st.data(), gap.data(), lag.data()),
         CLAB_OK);
  tk_nan(T, "rv.nan", st[5]);
}

}  // namespace

void t_rpv(tk &T, int big) {
  (void)big;
  famtest(T);
}
