// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of exprspplr() (src/plr.cpp). The R wrapper
// (R/aexp.R) checks the arguments and draws the resamples and folds with
// R's own generator; these fit them. Nothing calls back into the R
// evaluator.
#include "bcol.h"
#include "rchk.h"

extern "C" {

// One fit per sample: rows (from 1) and folds concatenated, lens their
// lengths. Out: list(ate, se), NaN where the treatment had no variation.
SEXP C_clab_r_plr(SEXP y, SEXP d, SEXP x, SEXP rows, SEXP lens, SEXP folds) {
  clab::who = "exprspplr";
  const int n = static_cast<int>(clab::xlen(y, "y")), p = Rf_isNull(x) ? 0 : Rf_ncols(x);
  const int ns = static_cast<int>(clab::xlen(lens, "lens"));
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(r, 0, Rf_allocVector(REALSXP, ns));
  SET_VECTOR_ELT(r, 1, Rf_allocVector(REALSXP, ns));
  std::vector<int> rw(static_cast<size_t>(clab::xlen(rows, "rows")));
  for (R_xlen_t i = 0; i < clab::xlen(rows, "rows"); ++i) rw[static_cast<size_t>(i)] = INTEGER(rows)[i] - 1;
  R_xlen_t at = 0;
  for (int s = 0; s < ns; ++s) {
    const int m = INTEGER(lens)[s];
    double o[2];
    const int rc = clab_plr(REAL(y), REAL(d), p ? REAL(x) : nullptr, n, p, rw.data() + at, m,
                            INTEGER(folds) + at, 5, o);
    clab::ichk(rc);
    if (rc != CLAB_OK && rc != CLAB_EDOM) clab::fail("%s", clab_strerror(rc));
    REAL(VECTOR_ELT(r, 0))[s] = o[0];
    REAL(VECTOR_ELT(r, 1))[s] = o[1];
    at += m;
  }
  UNPROTECT(1);
  return r;
}

SEXP C_clab_r_bootci(SEXP v, SEXP pr) {
  SEXP o = PROTECT(Rf_allocVector(REALSXP, 1 + clab::xlen(pr, "pr")));
  clab_bootci(REAL(v), static_cast<int>(clab::xlen(v, "v")), REAL(pr), static_cast<int>(clab::xlen(pr, "pr")), REAL(o));
  UNPROTECT(1);
  return o;
}

}  // extern "C"
