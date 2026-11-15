// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/ibn.cpp: the benchmark, z and false discovery rate by hand
// (no circumstances: the weights are equal), and a weighted case whose
// propensity is known in closed form (one indicator column).
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "tchk.h"

namespace {

using vd = std::vector<double>;
using vi = std::vector<int>;

void famtest(tk &T) {
  const vd y = {1, 0, 1, 1, 0, 0, 1, 0};
  const vi o = {0, 0, 0, 0, 1, 1, 1, 1};
  vd r(16), nl(3);
  vi ok(2);
  tk_ieq(T, "ib.ok", clab_ibench(y.data(), o.data(), nullptr, 8, 0, 2, 2, 0, 0, r.data(), ok.data(), nl.data()), CLAB_OK);
  tk_eq(T, "ib.rate", r[2], 0.75, 1e-15);
  tk_eq(T, "ib.bench", r[4], 0.25, 1e-12);
  tk_eq(T, "ib.ess", r[6], 4, 1e-12);
  const double z = 0.5 / std::sqrt(0.09375);
  tk_eq(T, "ib.z", r[8], z, 1e-12);
  tk_eq(T, "ib.zb", r[9], -z, 1e-12);
  tk_eq(T, "ib.q", r[12], std::erfc(z / std::sqrt(2.0)), 1e-12);
  // one 0/1 column: the propensity is the officer's share within each
  // level, so the weighted benchmark is the others' rate in each level
  // averaged over the officer's own stops by level
  const vd y2 = {1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0};
  const vd x2 = {1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0};
  const vi o2 = {0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1};
  vd r2(16);
  tk_ieq(T, "ib.x", clab_ibench(y2.data(), o2.data(), x2.data(), 12, 1, 2, 2, 0, 0, r2.data(), ok.data(), nl.data()), CLAB_OK);
  // officer 0: 3 stops with x = 1, 1 with x = 0; the others: rate 2/3 at
  // x = 1 (3 stops) and 2/5 at x = 0 (5 stops)
  tk_near(T, "ib.xb", r2[4], 0.75 * 2.0 / 3.0 + 0.25 * 0.4, 1e-8);
  vd r3(24);
  vi ok3(3);
  tk_ieq(T, "ib.y", clab_ibench(x2.data(), o.data(), nullptr, 8, 0, 3, 1, 0, 0, r3.data(), ok3.data(), nl.data()), CLAB_OK);
  tk_nan(T, "ib.empty", r3[2 * 3 + 2]);
  const vd bad = {2, 0};
  // doubly robust with one 0/1 column: the regression of y on the
  // indicator and x with these weights fits each level's weighted rate,
  // so the benchmark is again the others' rate by level averaged over the
  // officer's stops
  tk_ieq(T, "ib.dr", clab_ibench(y2.data(), o2.data(), x2.data(), 12, 1, 2, 2, 1, 0, r2.data(), ok.data(), nl.data()),
         CLAB_OK);
  tk_near(T, "ib.drb", r2[4], 0.75 * 2.0 / 3.0 + 0.25 * 0.4, 1e-6);
  tk_nan(T, "ib.nolf", r2[14]);
  tk_ieq(T, "ib.dom", clab_ibench(bad.data(), o.data(), nullptr, 2, 0, 2, 1, 0, 0, r.data(), ok.data(), nl.data()), CLAB_EDOM);
}

}  // namespace

void t_ibn(tk &T, int big) {
  (void)big;
  famtest(T);
}
