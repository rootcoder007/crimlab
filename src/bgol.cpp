// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry point of gologt() (src/gol.cpp). The R wrapper (R/gol.R)
// checks the arguments and passes the categories as codes from 0, the
// covariates as a double matrix and the parallel ones as 0/1; this fits
// the model, and the proportional-odds model when a slope is free, and
// returns the estimates, their Wald intervals and the tests. Nothing calls
// back into the R evaluator.
#include <cmath>

#include "crimlab.h"
#include "rchk.h"

extern "C" double Rf_pnorm5(double, double, double, int, int);

extern "C" {

SEXP C_clab_r_gologt(SEXP y, SEXP M, SEXP x, SEXP par, SEXP maxit, SEXP level) {
  clab::who = "gologt";
  const int n = Rf_length(y), m = Rf_asInteger(M), p = Rf_ncols(x), mi = Rf_asInteger(maxit);
  const int J = m - 1;
  int np = J, free = 0;
  for (int k = 0; k < p; ++k) {
    np += INTEGER(par)[k] ? 1 : J;
    free += INTEGER(par)[k] ? 0 : 1;
  }
  const char *dom = "the likelihood could not be maximised from its start.";
  SEXP v[5];
  v[0] = PROTECT(Rf_allocMatrix(REALSXP, np, 6));
  v[1] = PROTECT(Rf_allocMatrix(REALSXP, np, np));
  v[2] = PROTECT(Rf_allocMatrix(REALSXP, 3, p));
  v[3] = PROTECT(Rf_allocVector(REALSXP, 3));
  v[4] = PROTECT(Rf_allocVector(REALSXP, 3));
  double *cf = REAL(v[0]), *st = REAL(v[3]), *lr = REAL(v[4]);
  int rc = clab_gologt(INTEGER(y), n, m, REAL(x), p, INTEGER(par), mi, cf, REAL(v[1]), REAL(v[2]), st);
  clab::ichk(rc);
  if (rc == CLAB_EDOM) clab::fail("%s", dom);
  clab::kchk(rc, nullptr);
  const double z = Rf_qnorm5(0.5 + 0.5 * Rf_asReal(level), 0.0, 1.0, 1, 0);
  for (int k = 0; k < np; ++k) {
    const double e = cf[k], s = std::sqrt(REAL(v[1])[static_cast<R_xlen_t>(k) * np + k]);
    cf[np + k] = s;
    cf[2 * np + k] = e / s;
    cf[3 * np + k] = 2 * Rf_pnorm5(-std::fabs(e / s), 0.0, 1.0, 1, 0);
    cf[4 * np + k] = e - z * s;
    cf[5 * np + k] = e + z * s;
  }
  lr[0] = lr[1] = lr[2] = NA_REAL;
  if (free > 0) {
    // the proportional-odds model: every slope held parallel
    const int n0 = J + p;
    // R's transient memory: nothing with a destructor when an error jumps
    int *all = reinterpret_cast<int *>(R_alloc(static_cast<std::size_t>(p), sizeof(int)));
    const std::size_t ne = static_cast<std::size_t>(n0) * (n0 + 1) + 3 * static_cast<std::size_t>(p) + 3;
    double *e0 = reinterpret_cast<double *>(R_alloc(ne, sizeof(double)));
    for (int k = 0; k < p; ++k) all[k] = 1;
    double *c0 = e0 + n0, *w0 = c0 + static_cast<std::size_t>(n0) * n0, *s0 = w0 + 3 * p;
    rc = clab_gologt(INTEGER(y), n, m, REAL(x), p, all, mi, e0, c0, w0, s0);
    clab::ichk(rc);
    if (rc == CLAB_EDOM) clab::fail("%s", dom);
    clab::kchk(rc, nullptr);
    lr[0] = std::fmax(2 * (st[0] - s0[0]), 0.0);
    lr[1] = static_cast<double>(free) * (J - 1);
    lr[2] = clab_pchisq(lr[0], lr[1], 0);
  }
  static const char *const nm[] = {"coef", "vcov", "wald", "fit", "lr"};
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
  for (int i = 0; i < 5; ++i) SET_VECTOR_ELT(out, i, v[i]);
  SEXP res = PROTECT(clab::nlist(out, nm));
  UNPROTECT(7);
  return res;
}

}  // extern "C"
