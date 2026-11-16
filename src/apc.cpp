// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, age-period-cohort models of rates: the Poisson log-linear
// model of a table of counts by age and period with effects of age, period
// and cohort (period - age), identified by the intrinsic estimator (Yang,
// Fu and Land 2004), the solution orthogonal to the null vector of the
// design.
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

int clab_apcie(const double *y, const double *expo, int na, int np, int maxit, double *est, double *se,
               double *st) CLAB_TRY {
  if (na < 3 || np < 3 || maxit < 1 || static_cast<double>(na) * np > 1e6) return CLAB_EARG;
  const int nc = na + np - 1, n = na * np, p = 1 + (na - 1) + (np - 1) + (nc - 1);
  const std::size_t un = static_cast<std::size_t>(n), up = static_cast<std::size_t>(p);
  for (int i = 0; i < n; ++i) {
    if (!kutil::nonneg(y[i]) || !kutil::pos(expo[i])) return CLAB_EDOM;
  }
  // the design with effect coding: cell (i, j) of the column-major table,
  // cohort k = j - i + na - 1 (0 the oldest)
  std::vector<double> X(un * up, 0.0), off(un);
  auto eff = [&](std::size_t row, int c0, int lev, int nl) {
    if (lev < nl - 1) {
      X[static_cast<std::size_t>(c0 + lev) * un + row] = 1.0;
    } else {
      for (int l = 0; l < nl - 1; ++l) X[static_cast<std::size_t>(c0 + l) * un + row] = -1.0;
    }
  };
  for (int j = 0; j < np; ++j) {
    for (int i = 0; i < na; ++i) {
      const std::size_t r = static_cast<std::size_t>(j) * na + i;
      X[r] = 1.0;
      eff(r, 1, i, na);
      eff(r, na, j, np);
      eff(r, na + np - 1, j - i + na - 1, nc);
      off[r] = std::log(expo[r]);
    }
  }
  std::vector<double> b(up), vb(up * up), vc(1), s8(8);
  std::vector<int> al(up);
  std::vector<double> w(un * up + 6 * un + up + up * up + 2 * un);
  std::vector<int> iw(up + 1 + un);
  const int rc = clab_cntreg(y, X.data(), n, p, off.data(), nullptr, 0, 0, maxit, 1e-10, b.data(), vb.data(), vc.data(),
                             s8.data(), al.data(), w.data(), iw.data());
  if (rc != CLAB_OK) return rc;
  // one column aliased (the null vector's last nonzero element); 0 there
  int nal = 0;
  for (int k = 0; k < p; ++k) {
    if (al[static_cast<std::size_t>(k)]) {
      ++nal;
      b[static_cast<std::size_t>(k)] = 0.0;
      for (int q = 0; q < p; ++q) {
        vb[static_cast<std::size_t>(k) * up + q] = 0.0;
        vb[static_cast<std::size_t>(q) * up + k] = 0.0;
      }
    }
  }
  if (nal != 1) return CLAB_EDOM;
  // the null vector of the effect-coded design and the projection on its
  // orthogonal complement
  std::vector<double> B0(up, 0.0);
  for (int i = 0; i < na - 1; ++i) B0[static_cast<std::size_t>(1 + i)] = (i + 1) - 0.5 * (na + 1);
  for (int j = 0; j < np - 1; ++j) B0[static_cast<std::size_t>(na + j)] = -((j + 1) - 0.5 * (np + 1));
  for (int k = 0; k < nc - 1; ++k) B0[static_cast<std::size_t>(na + np - 1 + k)] = (k + 1) - 0.5 * (nc + 1);
  double bb = 0, bd = 0;
  for (std::size_t k = 0; k < up; ++k) {
    bb += B0[k] * B0[k];
    bd += B0[k] * b[k];
  }
  for (std::size_t k = 0; k < up; ++k) b[k] -= bd / bb * B0[k];
  // V = P V* P with P = I - B0 B0' / bb
  std::vector<double> PV(up * up), V(up * up);
  for (std::size_t c = 0; c < up; ++c) {
    for (std::size_t r = 0; r < up; ++r) {
      double s = vb[c * up + r];
      for (std::size_t q = 0; q < up; ++q) s -= B0[r] * B0[q] / bb * vb[c * up + q];
      PV[c * up + r] = s;
    }
  }
  for (std::size_t c = 0; c < up; ++c) {
    for (std::size_t r = 0; r < up; ++r) {
      double s = PV[c * up + r];
      for (std::size_t q = 0; q < up; ++q) s -= PV[q * up + r] * B0[q] * B0[c] / bb;
      V[c * up + r] = s;
    }
  }
  // est, se: the intercept, then every level of age, period and cohort (the
  // last of each minus the sum of the others)
  est[0] = b[0];
  se[0] = std::sqrt(std::fmax(V[0], 0.0));
  int o = 1;
  const int c0s[3] = {1, na, na + np - 1}, nls[3] = {na, np, nc};
  std::vector<double> L(up);
  for (int f = 0; f < 3; ++f) {
    for (int l = 0; l < nls[f]; ++l, ++o) {
      for (double &e : L) e = 0.0;
      if (l < nls[f] - 1) {
        L[static_cast<std::size_t>(c0s[f] + l)] = 1.0;
      } else {
        for (int q = 0; q < nls[f] - 1; ++q) L[static_cast<std::size_t>(c0s[f] + q)] = -1.0;
      }
      double e = 0, v = 0;
      for (std::size_t r = 0; r < up; ++r) {
        e += L[r] * b[r];
        for (std::size_t c = 0; c < up; ++c) v += L[r] * V[c * up + r] * L[c];
      }
      est[o] = e;
      se[o] = std::sqrt(std::fmax(v, 0.0));
    }
  }
  // st: deviance, residual df, log-likelihood, steps, converged
  st[0] = s8[7];
  st[1] = n - (p - 1);
  st[2] = s8[0];
  st[3] = s8[4];
  st[4] = s8[5];
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)
