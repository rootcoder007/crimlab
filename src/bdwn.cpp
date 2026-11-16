// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of dwnexp() and tailwd() (src/dwn.cpp). The R wrapper (R/dwn.R)
// checks the arguments and passes doubles, the period codes from 0 and,
// for the wedge, the polygons' offsets; this calls the kernel and returns
// list(out, out2), out2 NULL unless the rule is the angle; tailwd() its
// vector. Nothing calls
// back into the R evaluator.
#include "crimlab.h"
#include "kutil.h"
#include "rchk.h"

extern "C" {

SEXP C_clab_r_dwnexp(SEXP rule, SEXP sx, SEXP sy, SEXP sp, SEXP sto, SEXP rx, SEXP ry, SEXP ro,
                     SEXP nr, SEXP np, SEXP half, SEXP rad, SEXP ll, SEXP wrap) {
  clab::who = "dwnexp";
  const int ru = Rf_asInteger(rule), ns = clab::xlen(sx, "sx");
  const int nrr = Rf_asInteger(nr), npp = Rf_asInteger(np);
  if (static_cast<double>(nrr) * npp > kutil::maxcell) {
    clab::fail("the places times the periods (%.0f) must be at most %.0f.",
               static_cast<double>(nrr) * npp, kutil::maxcell);
  }
  SEXP o1 = PROTECT(Rf_allocMatrix(REALSXP, nrr, npp));
  SEXP o2 = PROTECT(ru == 2 ? Rf_allocMatrix(REALSXP, nrr, npp) : R_NilValue);
  const int rc = clab_dwnexp(ru, ns, REAL(sx), REAL(sy), INTEGER(sp), REAL(sto), nrr, REAL(rx),
                             REAL(ry), Rf_isNull(ro) ? nullptr : INTEGER(ro), npp,
                             Rf_asReal(half), Rf_asReal(rad), Rf_asLogical(ll),
                             Rf_asLogical(wrap), REAL(o1), ru == 2 ? REAL(o2) : nullptr);
  clab::ichk(rc);
  if (rc != CLAB_OK) clab::fail("%s", clab_strerror(rc));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, o1);
  SET_VECTOR_ELT(out, 1, o2);
  UNPROTECT(3);
  return out;
}

SEXP C_clab_r_tailwd(SEXP dir, SEXP spd, SEXP course, SEXP grp, SEXP ng, SEXP agg) {
  clab::who = "tailwd";
  const int n = clab::xlen(dir, "dir"), g = Rf_isNull(grp) ? 0 : Rf_asInteger(ng);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, Rf_isNull(grp) ? n : g));
  const int rc = clab_tailwd(n, REAL(dir), REAL(spd), Rf_asReal(course),
                             Rf_isNull(grp) ? nullptr : INTEGER(grp), g, Rf_asInteger(agg),
                             REAL(out));
  clab::ichk(rc);
  if (rc != CLAB_OK) clab::fail("%s", clab_strerror(rc));
  UNPROTECT(1);
  return out;
}

}  // extern "C"
