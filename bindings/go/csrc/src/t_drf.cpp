// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/drf.cpp: small cases by hand.
#include <cmath>

#include "crimlab.h"
#include "tchk.h"

void t_drf(tk &T, int big) {
  (void)big;
  const double x[4] = {1, 2, 3, 4}, y[4] = {3, 4, 5, 6};
  double r[4];
  clab_ks(x, 4, y, 4, r);
  tk_eq(T, "dr.ks", r[0], 0.5, 1e-15);
  tk_eq(T, "dr.neff", r[2], 2, 1e-15);
  tk_true(T, "dr.ksp", r[1] > 0 && r[1] <= 1);
  clab_ks(x, 4, x, 4, r);
  tk_eq(T, "dr.ks0", r[0], 0, 0);
  tk_eq(T, "dr.ksp0", r[1], 1, 0);
  const double p[2] = {0.5, 0.5}, q[2] = {1, 0};
  clab_psi(p, q, 2, 1e-6, r);
  tk_eq(T, "dr.psi", r[0], (0.5 - 1) * std::log(0.5) + (0.5 - 1e-6) * std::log(0.5 / 1e-6), 1e-14);
  tk_true(T, "dr.js", std::fabs(r[1] - (0.5 * 0.5 * std::log(0.5 / 0.25) + 0.5 * 0.5 * std::log(0.5 / 0.75) +
                                        0.5 * std::log(1 / 0.75))) < 1e-14);
  const double o[3] = {40, 60, 100}, e[3] = {50, 30, 20};
  clab_chig(o, e, 3, r);
  tk_eq(T, "dr.chig", r[0], 126, 1e-13);
  tk_eq(T, "dr.chidf", r[1], 2, 0);
  const double a[2] = {30, 20}, b[2] = {25, 25};
  clab_chih(a, b, 2, r);
  tk_eq(T, "dr.chih", r[0], 100.0 / 99, 1e-14);
  tk_eq(T, "dr.chimin", r[2], 22.5, 1e-14);
  double px[5], py[5];
  const int k = clab_psibin(x, 4, y, 4, 2, px, py);
  tk_ieq(T, "dr.bins", k, 2);
  tk_eq(T, "dr.bin0", px[0], 0.5, 0);
  tk_eq(T, "dr.bin1", py[1], 1, 0);
}
