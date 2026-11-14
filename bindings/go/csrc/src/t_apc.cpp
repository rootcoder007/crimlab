// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/apc.cpp: counts that equal the expected counts of effects
// orthogonal to the null vector give those effects back.
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "tchk.h"

void t_apc(tk &T, int big) {
  (void)big;
  const int na = 4, np = 3, nc = 6, p = 11;
  // effect-coded coefficients: intercept, age 1..3, period 1..2, cohort
  // 1..5
  std::vector<double> u = {-3, 0.4, 0.1, -0.2, 0.3, -0.1, -0.5, -0.2, 0.1, 0.25, 0.3};
  std::vector<double> B0(p, 0.0);
  for (int i = 0; i < na - 1; ++i) B0[static_cast<std::size_t>(1 + i)] = (i + 1) - 0.5 * (na + 1);
  for (int j = 0; j < np - 1; ++j) B0[static_cast<std::size_t>(na + j)] = -((j + 1) - 0.5 * (np + 1));
  for (int k = 0; k < nc - 1; ++k) B0[static_cast<std::size_t>(na + np - 1 + k)] = (k + 1) - 0.5 * (nc + 1);
  double bb = 0, bd = 0;
  for (int k = 0; k < p; ++k) {
    bb += B0[static_cast<std::size_t>(k)] * B0[static_cast<std::size_t>(k)];
    bd += B0[static_cast<std::size_t>(k)] * u[static_cast<std::size_t>(k)];
  }
  for (int k = 0; k < p; ++k) u[static_cast<std::size_t>(k)] -= bd / bb * B0[static_cast<std::size_t>(k)];
  // every level's effect (the last minus the sum of the others)
  auto lev = [&](int c0, int nl, int l) {
    if (l < nl - 1) return u[static_cast<std::size_t>(c0 + l)];
    double s = 0;
    for (int q = 0; q < nl - 1; ++q) s -= u[static_cast<std::size_t>(c0 + q)];
    return s;
  };
  std::vector<double> y(12), e(12, 1000.0);
  for (int j = 0; j < np; ++j) {
    for (int i = 0; i < na; ++i) {
      const double eta = u[0] + lev(1, na, i) + lev(na, np, j) + lev(na + np - 1, nc, j - i + na - 1);
      y[static_cast<std::size_t>(j * na + i)] = 1000.0 * std::exp(eta);
    }
  }
  std::vector<double> est(1 + na + np + nc), se(est.size()), st(5);
  tk_ieq(T, "ap.ok", clab_apcie(y.data(), e.data(), na, np, 100, est.data(), se.data(), st.data()), CLAB_OK);
  tk_near(T, "ap.int", est[0], u[0], 1e-8);
  tk_near(T, "ap.age4", est[4], lev(1, na, 3), 1e-8);
  tk_near(T, "ap.per3", est[1 + na + 2], lev(na, np, 2), 1e-8);
  tk_near(T, "ap.coh6", est[1 + na + np + 5], lev(na + np - 1, nc, 5), 1e-8);
  tk_near(T, "ap.dev", st[0], 0, 1e-8);
  tk_eq(T, "ap.df", st[1], 2, 0);
  tk_true(T, "ap.se", se[1] > 0);
  tk_ieq(T, "ap.small", clab_apcie(y.data(), e.data(), 2, 6, 100, est.data(), se.data(), st.data()), CLAB_EARG);
  e[3] = 0;
  tk_ieq(T, "ap.dom", clab_apcie(y.data(), e.data(), na, np, 100, est.data(), se.data(), st.data()), CLAB_EDOM);
}
