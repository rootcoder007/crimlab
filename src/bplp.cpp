// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the plot data of R/plotm.R (src/plp.cpp) and of
// the densities of R/plotx.R (src/dens.cpp). The plot methods check the
// object and draw; these compute what is drawn.
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

// list(x, y, bw, n, warn) of .dens() (R/plotx.R): stats::density() by
// clab_dens() (src/dens.cpp), the kernel of the command line's plotrdg()
// and plotvln(). w: NULL or one weight per value; ip: the rule, the kernel
// and the number of points; dp: bw, adjust, cut, from and to (from and to
// NA when not given). R's messages for the data it refuses.
SEXP C_clab_r_dens(SEXP x, SEXP w, SEXP ip, SEXP dp) {
  clab::who = "dens";
  const int n = clab::xlen(x, "x");
  if (TYPEOF(x) != REALSXP || TYPEOF(ip) != INTSXP || XLENGTH(ip) != 3 || TYPEOF(dp) != REALSXP ||
      XLENGTH(dp) != 5 || (!Rf_isNull(w) && (TYPEOF(w) != REALSXP || XLENGTH(w) != n))) {
    clab::fail("%s", clab_strerror(CLAB_EARG));
  }
  const int *pi = INTEGER(ip);
  const double *pd = REAL(dp);
  const R_xlen_t m = pi[2] > 0 ? pi[2] : 0;
  SEXP xo = PROTECT(Rf_allocVector(REALSXP, m));
  SEXP yo = PROTECT(Rf_allocVector(REALSXP, m));
  double st[2] = {0.0, 0.0};
  int why = 0;
  const int rc = clab_dens(REAL(x), Rf_isNull(w) ? nullptr : REAL(w), n, pi[0], pd[0], pd[1], pi[1], pi[2], pd[2],
                           !ISNAN(pd[3]), pd[3], !ISNAN(pd[4]), pd[4], REAL(xo), REAL(yo), st, &why);
  clab::ichk(rc);
  static const char *const msg[] = {"need at least 2 points to select a bandwidth automatically",
                                    "data are constant in bandwidth calculation",
                                    "sample is too sparse to find TD",
                                    "sample is too sparse to find alph2",
                                    "no solution in the specified range of bandwidths",
                                    "non-finite 'bw'",
                                    "'bw' is not positive.",
                                    "non-finite 'from'",
                                    "",
                                    "non-finite 'to'",
                                    "'xmin' not less than 'xmax'"};
  if (rc == CLAB_EDOM) clab::fail("%s", msg[why - 1]);
  if (rc != CLAB_OK) clab::fail("%s", clab_strerror(rc));
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 5));
  SET_VECTOR_ELT(r, 0, xo);
  SET_VECTOR_ELT(r, 1, yo);
  SET_VECTOR_ELT(r, 2, Rf_ScalarReal(st[0]));
  SET_VECTOR_ELT(r, 3, Rf_ScalarReal(st[1]));
  SET_VECTOR_ELT(r, 4, Rf_ScalarLogical(why == 9));
  static const char *const nm[] = {"x", "y", "bw", "n", "warn"};
  SEXP o = clab::nlist(r, nm);
  UNPROTECT(3);
  return o;
}

}  // extern "C"
