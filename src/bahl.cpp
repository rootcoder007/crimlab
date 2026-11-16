// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the health side of air pollution (src/ahl.cpp):
// the concentration-response parameters of crfpm25() and crfno2(), and the
// burden of each area of polbrdarea(). The R wrappers check the arguments;
// these call the kernels and build the results. Nothing calls back into
// the R evaluator.
#include <cstddef>
#include <cstdint>
#include <vector>

#include "crimlab.h"
#include "rchk.h"

namespace {

void kc(int rc) {
  clab::ichk(rc);
  if (rc != CLAB_OK) clab::fail("%s", clab_strerror(rc));
}

}  // namespace

extern "C" {

// list(form, par (3), ci (2)) of pollutant pm (1 PM2.5, 0 NO2) and the
// outcome's zero-based place in its table.
SEXP C_clab_r_crfpar(SEXP pm, SEXP out) {
  clab::who = "crf";
  int form = 0;
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(r, 1, Rf_allocVector(REALSXP, 3));
  SET_VECTOR_ELT(r, 2, Rf_allocVector(REALSXP, 2));
  kc(clab_crfpar(Rf_asInteger(pm), Rf_asInteger(out), &form, REAL(VECTOR_ELT(r, 1)), REAL(VECTOR_ELT(r, 2))));
  SET_VECTOR_ELT(r, 0, Rf_ScalarInteger(form));
  UNPROTECT(1);
  return r;
}

// list(rr, paf, attrcases, basecases, bad, order): the areas worst first
// with order their rows in the table; bad the first row a check of
// polbrdn() refuses (NA when none), the other values then unset.
SEXP C_clab_r_pbarea(SEXP e, SEXP pop, SEXP br, SEXP key, SEXP rare, SEXP hirate) {
  clab::who = "polbrdarea";
  const int n = clab::xlen(e, "expcol");
  if (XLENGTH(pop) != n || XLENGTH(br) != n || XLENGTH(key) != 2) kc(CLAB_EARG);
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 6));
  for (int j = 0; j < 4; ++j) SET_VECTOR_ELT(r, j, Rf_allocVector(REALSXP, n));
  SET_VECTOR_ELT(r, 5, Rf_allocVector(INTSXP, n));
  std::vector<std::int64_t> o(static_cast<std::size_t>(n), 0);
  std::int64_t bad = 0;
  kc(clab_pbarea(static_cast<std::size_t>(n), REAL(e), REAL(pop), REAL(br), INTEGER(key)[0], INTEGER(key)[1],
                 Rf_asLogical(rare) == TRUE ? 1 : 0, Rf_asLogical(hirate) == TRUE ? 1 : 0, REAL(VECTOR_ELT(r, 0)),
                 REAL(VECTOR_ELT(r, 1)), REAL(VECTOR_ELT(r, 2)), REAL(VECTOR_ELT(r, 3)), o.data(), &bad));
  for (int i = 0; i < n; ++i) INTEGER(VECTOR_ELT(r, 5))[i] = static_cast<int>(o[static_cast<std::size_t>(i)]);
  SET_VECTOR_ELT(r, 4, Rf_ScalarInteger(bad > 0 ? static_cast<int>(bad) : NA_INTEGER));
  UNPROTECT(1);
  return r;
}

}  // extern "C"
