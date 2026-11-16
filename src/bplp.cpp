// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the plot data of R/plotm.R (src/plp.cpp). The
// plot methods check the object and draw; these compute what is drawn.
// Nothing calls back into the R evaluator.
#include <cstddef>

#include "crimlab.h"
#include "rchk.h"

extern "C" {

// list(fpr, tpr, auc) of plot.lnkpr(): the ROC curve from (0, 0) over the
// probabilities prob of the pairs, linked the logical flags.
SEXP C_clab_r_rocpts(SEXP prob, SEXP lk) {
  clab::who = "plot";
  const int n = clab::xlen(prob, "prob");
  if (TYPEOF(prob) != REALSXP || TYPEOF(lk) != LGLSXP || XLENGTH(lk) != n) clab::fail("%s", clab_strerror(CLAB_EARG));
  SEXP f = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(n) + 1));
  SEXP t = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(n) + 1));
  std::size_t k = 0;
  double auc = 0.0;
  int na[2] = {0, 0};
  const int rc = clab_rocpts(static_cast<std::size_t>(n), REAL(prob), LOGICAL(lk), REAL(f), REAL(t), &k, &auc, na);
  clab::ichk(rc);
  if (rc != CLAB_OK) clab::fail("%s", clab_strerror(rc));
  for (std::size_t i = 1; i < k; ++i) {
    if (na[0] != 0) REAL(f)[i] = NA_REAL;
    if (na[1] != 0) REAL(t)[i] = NA_REAL;
  }
  if ((na[0] != 0 || na[1] != 0) && k > 1) auc = NA_REAL;
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(r, 0, Rf_xlengthgets(f, static_cast<R_xlen_t>(k)));
  SET_VECTOR_ELT(r, 1, Rf_xlengthgets(t, static_cast<R_xlen_t>(k)));
  SET_VECTOR_ELT(r, 2, Rf_ScalarReal(auc));
  UNPROTECT(3);
  return r;
}

}  // extern "C"
