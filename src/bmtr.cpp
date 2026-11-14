// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry point of the bounds under monotone treatment response and
// selection (mtrbnd, src/mtr.cpp). Nothing calls back into the R
// evaluator.
#include <cmath>
#include <vector>

#include "bcol.h"
#include "rchk.h"
#include "spin.h"

using clab::fail;

extern "C" {

SEXP C_clab_r_mtrbnd(SEXP y, SEXP z, SEXP ylo, SEXP yhi, SEXP assume, SEXP w) {
  clab::who = "mtrbnd";
  int np = 0;
  SEXP yv = PROTECT(bj::nums(y, "y"));
  ++np;
  const int n = bj::len(yv, "y");
  csp::lab g = bj::labs(z, "z", n);
  np += 2;
  if (g.g < 2) fail("`z` must hold two treatments at least.");
  const double lo = clab::one(ylo, "ylo", -HUGE_VAL, HUGE_VAL);
  const double hi = clab::one(yhi, "yhi", -HUGE_VAL, HUGE_VAL);
  if (!(lo < hi)) fail("`ylo` must be below `yhi`.");
  for (int i = 0; i < n; ++i) {
    if (REAL(yv)[i] < lo || REAL(yv)[i] > hi) fail("`y` must lie between `ylo` and `yhi`.");
  }
  static const char *const an[] = {"none", "mtr", "mts", "both"};
  const int as = csp::pick(assume, "assume", an, 4);
  SEXP wv = R_NilValue;
  if (!Rf_isNull(w)) {
    wv = PROTECT(bj::nums(w, "w", n, 0.0));
    ++np;
  }
  const int k = g.g;
  std::vector<double> o(6 * static_cast<std::size_t>(k));
  const int rc = clab_mtrbnd(REAL(yv), INTEGER(g.code), Rf_isNull(wv) ? nullptr : REAL(wv), n, k, lo, hi, as, o.data());
  if (rc != CLAB_OK) fail("every treatment needs people with a positive weight.");
  bj::cols f(7, k);
  f.put("z", g.lev);
  static const char *const cn[] = {"share", "mean", "lower", "upper", "atelo", "atehi"};
  for (int c = 0; c < 6; ++c) {
    double *d = f.dbl(cn[c]);
    for (int t = 0; t < k; ++t) d[t] = bj::na(o[static_cast<std::size_t>(c) * k + t]);
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, f.done());
  SET_VECTOR_ELT(out, 1, Rf_mkString(an[as]));
  const double r[2] = {lo, hi};
  static const char *const rn[] = {"ylo", "yhi"};
  SET_VECTOR_ELT(out, 2, clab::nvec(2, rn, r));
  static const char *const on[] = {"bnd", "assume", "range"};
  bj::named(out, on);
  Rf_classgets(out, PROTECT(Rf_mkString("mtrbnd")));
  UNPROTECT(np + 3);
  return out;
}

}  // extern "C"
