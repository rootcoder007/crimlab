// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry point of the accuracy of a hot-spot map (hsacc,
// src/hsa.cpp). Nothing calls back into the R evaluator.
#include <cmath>
#include <vector>

#include "bcol.h"
#include "rchk.h"

using clab::fail;

extern "C" {

SEXP C_clab_r_hsacc(SEXP pred, SEXP obs, SEXP area, SEXP top, SEXP prev) {
  clab::who = "hsacc";
  int np = 0;
  SEXP ov = PROTECT(bj::nums(obs, "obs", -1, 0.0));
  ++np;
  const int n = bj::len(ov, "obs");
  const bool lg = TYPEOF(pred) == LGLSXP;
  std::vector<int> hot;
  SEXP pv = R_NilValue;
  if (lg) {
    clab::vec(pred, "pred");
    if (XLENGTH(pred) != n) fail("`pred` must have one value per cell, %d.", n);
    for (int i = 0; i < n; ++i) {
      if (LOGICAL(pred)[i] == NA_LOGICAL) fail("`pred` must not be NA.");
      hot.push_back(LOGICAL(pred)[i]);
    }
  } else {
    pv = PROTECT(bj::nums(pred, "pred", n));
    ++np;
  }
  SEXP av = R_NilValue;
  if (!Rf_isNull(area)) {
    av = PROTECT(bj::nums(area, "area", n, 0.0, true));
  } else {
    av = PROTECT(Rf_allocVector(REALSXP, n));
    for (int i = 0; i < n; ++i) REAL(av)[i] = 1.0;
  }
  ++np;
  SEXP tv = PROTECT(bj::nums(top, "top", -1, 0.0, true));
  ++np;
  const int nt = lg ? 1 : bj::len(tv, "top");
  for (R_xlen_t k = 0; k < XLENGTH(tv); ++k) {
    if (REAL(tv)[k] > 1) fail("`top` must be at most 1.");
  }
  SEXP rv = R_NilValue;
  if (!Rf_isNull(prev)) {
    rv = PROTECT(bj::nums(prev, "prev", n, 0.0));
    ++np;
  }
  double tot = 0;
  for (int i = 0; i < n; ++i) tot += REAL(ov)[i];
  if (!(tot > 0)) fail("`obs` has no crimes: every count is 0.");
  std::vector<double> res(11 * static_cast<std::size_t>(nt));
  bj::kc(clab_hsacc(lg ? nullptr : REAL(pv), lg ? hot.data() : nullptr, REAL(ov), REAL(av),
                    Rf_isNull(rv) ? nullptr : REAL(rv), n, REAL(tv), nt, res.data()),
         "the kernel refused the input.");
  static const char *const cn[] = {"top", "cells", "area", "ashare", "n", "N", "hit", "pai", "nmax", "pei", "rri"};
  bj::cols f(11, nt);
  for (int c = 0; c < 11; ++c) {
    double *d = f.dbl(cn[c]);
    for (int q = 0; q < nt; ++q) d[q] = bj::na(res[static_cast<std::size_t>(c) * nt + q]);
  }
  SEXP out = PROTECT(f.done());
  UNPROTECT(np + 2);
  return out;
}

}  // extern "C"
