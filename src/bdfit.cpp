// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of dstfit() (src/dfit.cpp). The R wrapper (R/dfit.R)
// checks the arguments and draws the bootstrap resamples with R's own
// generator. Nothing calls back into the R evaluator.
#include <cmath>
#include <vector>

#include "bcol.h"
#include "kutil.h"
#include "rchk.h"

extern "C" {

// The estimate, or NULL when the fit fails.
SEXP C_clab_r_dfit(SEXP fam, SEXP x) {
  const int f = Rf_asInteger(fam), k = clab_dfnp(f);
  SEXP p = PROTECT(Rf_allocVector(REALSXP, k));
  const int rc = clab_dfit(f, REAL(x), static_cast<int>(clab::xlen(x, "x")), REAL(p));
  UNPROTECT(1);
  return rc == CLAB_OK ? p : R_NilValue;
}

// The log-likelihood, KS and AD statistics of a fit.
SEXP C_clab_r_dfgof(SEXP fam, SEXP x, SEXP par) {
  const int f = Rf_asInteger(fam), n = static_cast<int>(clab::xlen(x, "x"));
  SEXP o = PROTECT(Rf_allocVector(REALSXP, 3));
  REAL(o)[0] = clab_dfll(f, REAL(x), n, REAL(par));
  clab_dfgof(f, REAL(x), n, REAL(par), REAL(o) + 1);
  if (f >= 8) REAL(o)[2] = NA_REAL;
  UNPROTECT(1);
  return o;
}

// The percentile limits (a, 1 - a; type 7, failed fits left out) of the
// refits of family f on the resamples idx (n x B, from 1): a 2 x k matrix.
SEXP C_clab_r_dfboot(SEXP fam, SEXP x, SEXP idx, SEXP a) {
  const int f = Rf_asInteger(fam), k = clab_dfnp(f), n = Rf_nrows(idx), B = Rf_ncols(idx);
  std::vector<std::vector<double>> est(static_cast<size_t>(k));
  std::vector<double> xb(static_cast<size_t>(n)), p(static_cast<size_t>(k));
  for (int b = 0; b < B; ++b) {
    if (kutil::stop(50.0 * n)) clab::ichk(CLAB_EINTR);
    for (int i = 0; i < n; ++i) xb[static_cast<size_t>(i)] = REAL(x)[INTEGER(idx)[static_cast<R_xlen_t>(b) * n + i] - 1];
    if (clab_dfit(f, xb.data(), n, p.data()) != CLAB_OK) continue;
    for (int j = 0; j < k; ++j) est[static_cast<size_t>(j)].push_back(p[static_cast<size_t>(j)]);
  }
  SEXP o = PROTECT(Rf_allocMatrix(REALSXP, 2, k));
  const double pr[2] = {Rf_asReal(a), 1 - Rf_asReal(a)};
  for (int j = 0; j < k; ++j) {
    const std::vector<double> &e = est[static_cast<size_t>(j)];
    double q[3];
    clab_bootci(e.data(), static_cast<int>(e.size()), pr, 2, q);
    REAL(o)[2 * j] = e.empty() ? NA_REAL : q[1];
    REAL(o)[2 * j + 1] = e.empty() ? NA_REAL : q[2];
  }
  UNPROTECT(1);
  return o;
}

// Hill's tail index from the k largest of the sorted x: alpha and its
// asymptotic limits at conflvl.
SEXP C_clab_r_dfhill(SEXP x, SEXP k, SEXP lvl) {
  const int n = static_cast<int>(clab::xlen(x, "x")), kk = Rf_asInteger(k);
  long double s = 0;
  for (int i = 0; i < kk; ++i) s += std::log(REAL(x)[n - 1 - i] / REAL(x)[n - 1 - kk]);
  const double al = kk / static_cast<double>(s), z = clab_qnorm(1 - (1 - Rf_asReal(lvl)) / 2);
  SEXP o = PROTECT(Rf_allocVector(REALSXP, 3));
  REAL(o)[0] = al;
  REAL(o)[1] = al * (1 - z / std::sqrt(static_cast<double>(kk)));
  REAL(o)[2] = al * (1 + z / std::sqrt(static_cast<double>(kk)));
  UNPROTECT(1);
  return o;
}

}  // extern "C"
