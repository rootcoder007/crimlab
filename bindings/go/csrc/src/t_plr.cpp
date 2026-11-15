// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/plr.cpp: an outcome exactly linear in the treatment and the
// covariate gives the coefficient back with no error, by both fits; and the
// bootstrap summary against its definition.
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "tchk.h"

void t_plr(tk &T, int) {
  const int n = 20;
  std::vector<double> x(n), d(n), y(n), x2(2 * n);
  std::vector<int> rows(n), fold(n);
  for (int i = 0; i < n; ++i) {
    x[i] = std::sin(i + 1.0);
    d[i] = x[i] + std::cos(3.0 * i);
    y[i] = 3 * d[i] + 2 * x[i];
    x2[i] = x[i];
    x2[n + i] = 2 * x[i];  // collinear: the QR fit, the copy aliased
    rows[i] = i;
    fold[i] = i % 5 + 1;
  }
  double o[2];
  tk_ieq(T, "plr.ok", clab_plr(y.data(), d.data(), x.data(), n, 1, rows.data(), n, fold.data(), 5, o), CLAB_OK);
  tk_eq(T, "plr.ate", o[0], 3, 1e-12);
  tk_true(T, "plr.se", std::fabs(o[1]) < 1e-10);
  clab_plr(y.data(), d.data(), x2.data(), n, 2, rows.data(), n, fold.data(), 5, o);
  tk_eq(T, "plr.qr", o[0], 3, 1e-10);
  // the treatment equal to the covariate leaves nothing
  tk_ieq(T, "plr.dom", clab_plr(y.data(), x.data(), x.data(), n, 1, rows.data(), n, fold.data(), 5, o), CLAB_EDOM);
  const double v[4] = {4, 1, 3, 2}, pr[2] = {0.5, 0.025};
  double b[3];
  clab_bootci(v, 4, pr, 2, b);
  tk_eq(T, "plr.sd", b[0], std::sqrt(5.0 / 3), 1e-15);
  tk_eq(T, "plr.q50", b[1], 2.5, 0);
  tk_eq(T, "plr.q025", b[2], 1.075, 1e-15);
}
