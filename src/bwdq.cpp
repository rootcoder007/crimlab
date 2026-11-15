// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry point of the weighted displacement quotient (wdqx,
// src/wdq.cpp). Nothing calls back into the R evaluator.
#include "bcol.h"
#include "rchk.h"

using clab::fail;

extern "C" {

SEXP C_clab_r_wdqx(SEXP a, SEXP b, SEXP c) {
  clab::who = "wdqx";
  SEXP av = PROTECT(bj::nums(a, "a", 2, 0.0));
  SEXP bv = PROTECT(bj::nums(b, "b", 2, 0.0));
  SEXP cv = PROTECT(bj::nums(c, "c", 2, 0.0, true));
  double o[8];
  bj::kc(clab_wdqx(REAL(av), REAL(bv), REAL(cv), o), "the kernel refused the input.");
  for (double &x : o) x = bj::na(x);
  static const char *const nm[] = {"sm", "bdm", "wdq", "ge", "ne", "ea", "eb", "tne"};
  SEXP out = PROTECT(clab::nvec(8, nm, o));
  UNPROTECT(4);
  return out;
}

}  // extern "C"
