// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the exposure-and-crime family (src/apx.cpp and
// clab_ppmlfe, clab_fesep of src/cnt.cpp). The R wrappers (R/apx.R) check
// the arguments and pass codes from 0, numeric matrices (with no columns
// for none) and NULL for what is left out; these call the kernels and
// build the results. Scratch comes from R_alloc, so an error leaks
// nothing. Nothing calls back into the R evaluator.
#include <cmath>

#include "bcol.h"
#include "rchk.h"

namespace {

int n32(SEXP x) { return clab::ilen(static_cast<double>(XLENGTH(x)), nullptr); }

const char *const absorbed =
    "the fit failed: a regressor or an instrument does not vary within the fixed effects (they "
    "absorb it), or the columns are collinear";

void kc(int rc) {
  clab::ichk(rc);
  if (rc == CLAB_EDOM) clab::fail("%s.", absorbed);
  if (rc != CLAB_OK) clab::fail("%s", clab_strerror(rc));
}

const int *icodes(SEXP x) { return Rf_isNull(x) ? nullptr : INTEGER(x); }
const double *dnul(SEXP x) { return Rf_isNull(x) ? nullptr : REAL(x); }

SEXP named(SEXP v, const char *const *nm, int k) {
  SEXP s = PROTECT(Rf_allocVector(STRSXP, k));
  for (int j = 0; j < k; ++j) SET_STRING_ELT(s, j, Rf_mkChar(nm[j]));
  Rf_setAttrib(v, R_NamesSymbol, s);
  UNPROTECT(1);
  return v;
}

}  // namespace

extern "C" {

SEXP C_clab_r_fesep(SEXP y, SEXP f1, SEXP n1, SEXP f2, SEXP n2) {
  clab::who = "ppmlfe";
  const int n = n32(y), a = Rf_asInteger(n1), b = Rf_isNull(f2) ? 0 : Rf_asInteger(n2);
  SEXP k = PROTECT(Rf_allocVector(LGLSXP, n));
  double *w = reinterpret_cast<double *>(R_alloc(static_cast<size_t>(a + b + 1), sizeof(double)));
  if (clab_fesep(REAL(y), n, INTEGER(f1), a, icodes(f2), b, LOGICAL(k), w) < 0) {
    clab::fail("%s", clab_strerror(CLAB_EARG));
  }
  UNPROTECT(1);
  return k;
}

SEXP C_clab_r_ppmlfe(SEXP y, SEXP x, SEXP off, SEXP f1, SEXP n1, SEXP f2, SEXP n2, SEXP cl,
                     SEXP ncl, SEXP cl2, SEXP ncl2, SEXP tol, SEXP maxit) {
  clab::who = "ppmlfe";
  const int n = n32(y), p = Rf_ncols(x), a = Rf_asInteger(n1);
  const int b = Rf_isNull(f2) ? 0 : Rf_asInteger(n2);
  const int g1 = Rf_isNull(cl) ? 0 : Rf_asInteger(ncl), g2 = Rf_isNull(cl2) ? 0 : Rf_asInteger(ncl2);
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(r, 0, Rf_allocVector(REALSXP, p));
  SET_VECTOR_ELT(r, 1, Rf_allocMatrix(REALSXP, p, p));
  SET_VECTOR_ELT(r, 2, Rf_allocVector(REALSXP, 6));
  const size_t sn = static_cast<size_t>(n), sp = static_cast<size_t>(p);
  double *w = reinterpret_cast<double *>(
      R_alloc(sn * (2 * sp + 8) + 3 * sp * sp + sp + 4 * static_cast<size_t>(a + b), sizeof(double)));
  int *iw = reinterpret_cast<int *>(
      R_alloc(sn + static_cast<size_t>(g1) * static_cast<size_t>(g2) + 1, sizeof(int)));
  kc(clab_ppmlfe(REAL(y), REAL(x), n, p, dnul(off), nullptr, INTEGER(f1), a, icodes(f2), b,
                 icodes(cl), g1, icodes(cl2), g2, Rf_asInteger(maxit), 1e-10, Rf_asReal(tol),
                 REAL(VECTOR_ELT(r, 0)), REAL(VECTOR_ELT(r, 1)), REAL(VECTOR_ELT(r, 2)), w, iw));
  static const char *const nm[] = {"beta", "vc", "st"};
  named(r, nm, 3);
  UNPROTECT(1);
  return r;
}

SEXP C_clab_r_ivfe(SEXP y, SEXP xe, SEXP xx, SEXP zi, SEXP f1, SEXP n1, SEXP f2, SEXP n2,
                   SEXP cl, SEXP ncl, SEXP cl2, SEXP ncl2, SEXP tol) {
  clab::who = "ivfe";
  const int n = n32(y), pe = Rf_ncols(xe), px = Rf_ncols(xx), q = Rf_ncols(zi), p = pe + px;
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 5));
  SET_VECTOR_ELT(r, 0, Rf_allocVector(REALSXP, p));
  SET_VECTOR_ELT(r, 1, Rf_allocMatrix(REALSXP, p, p));
  SET_VECTOR_ELT(r, 2, Rf_allocVector(REALSXP, pe));
  SET_VECTOR_ELT(r, 3, Rf_allocVector(REALSXP, 3));
  SET_VECTOR_ELT(r, 4, Rf_allocVector(REALSXP, n));
  kc(clab_ivfe(REAL(y), REAL(xe), pe, REAL(xx), px, REAL(zi), q, n, nullptr, INTEGER(f1),
               Rf_asInteger(n1), icodes(f2), Rf_isNull(f2) ? 0 : Rf_asInteger(n2), icodes(cl),
               Rf_isNull(cl) ? 0 : Rf_asInteger(ncl), icodes(cl2),
               Rf_isNull(cl2) ? 0 : Rf_asInteger(ncl2), Rf_asReal(tol), REAL(VECTOR_ELT(r, 0)),
               REAL(VECTOR_ELT(r, 1)), REAL(VECTOR_ELT(r, 2)), nullptr, REAL(VECTOR_ELT(r, 4)),
               REAL(VECTOR_ELT(r, 3))));
  static const char *const nm[] = {"beta", "vc", "fst", "st", "res"};
  named(r, nm, 5);
  UNPROTECT(1);
  return r;
}

SEXP C_clab_r_ppmlcf(SEXP y, SEXP xe, SEXP xx, SEXP zi, SEXP off, SEXP f1, SEXP n1, SEXP f2,
                     SEXP n2, SEXP cl, SEXP ncl, SEXP bm, SEXP tol, SEXP maxit) {
  clab::who = "ppmlcf";
  const int n = n32(y), pe = Rf_ncols(xe), px = Rf_ncols(xx), q = Rf_ncols(zi);
  const int k = 2 * pe + px, nb = Rf_ncols(bm);
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(r, 0, Rf_allocVector(REALSXP, k));
  SET_VECTOR_ELT(r, 1, Rf_allocVector(REALSXP, pe));
  SET_VECTOR_ELT(r, 2, Rf_allocMatrix(REALSXP, nb, k));
  SET_VECTOR_ELT(r, 3, Rf_allocVector(REALSXP, 3));
  kc(clab_ppmlcf(REAL(y), REAL(xe), pe, REAL(xx), px, REAL(zi), q, n, dnul(off), INTEGER(f1),
                 Rf_asInteger(n1), icodes(f2), Rf_isNull(f2) ? 0 : Rf_asInteger(n2), INTEGER(cl),
                 Rf_asInteger(ncl), INTEGER(bm), nb, Rf_asInteger(maxit), 1e-10, Rf_asReal(tol),
                 REAL(VECTOR_ELT(r, 0)), REAL(VECTOR_ELT(r, 1)), REAL(VECTOR_ELT(r, 2)),
                 REAL(VECTOR_ELT(r, 3))));
  static const char *const nm[] = {"beta", "fst", "draw", "st"};
  named(r, nm, 4);
  UNPROTECT(1);
  return r;
}

SEXP C_clab_r_windbin(SEXP dir, SEXP width, SEXP start, SEXP area, SEXP na, SEXP base) {
  clab::who = "windiv";
  const int n = n32(dir), a = Rf_asInteger(na);
  const int nb = static_cast<int>(std::nearbyint(360 / Rf_asReal(width)));
  SEXP r = PROTECT(Rf_allocMatrix(REALSXP, n, a * (nb - 1)));
  kc(clab_windbin(REAL(dir), n, Rf_asReal(width), Rf_asReal(start), icodes(area), a,
                  Rf_asInteger(base), REAL(r)));
  UNPROTECT(1);
  return r;
}

SEXP C_clab_r_dshift(SEXP u, SEXP t, SEXP k) {
  clab::who = "plcday";
  const int n = n32(t);
  SEXP r = PROTECT(Rf_allocVector(INTSXP, n));
  int *ord = reinterpret_cast<int *>(R_alloc(static_cast<size_t>(n), sizeof(int)));
  const int rc = clab_dshift(INTEGER(u), REAL(t), n, Rf_asReal(k), INTEGER(r), ord);
  clab::ichk(rc);
  if (rc == CLAB_EDOM) clab::fail("a unit has two rows on one day: give one row per unit and day.");
  if (rc != CLAB_OK) clab::fail("%s", clab_strerror(rc));
  for (int i = 0; i < n; ++i) INTEGER(r)[i] = INTEGER(r)[i] < 0 ? NA_INTEGER : INTEGER(r)[i] + 1;
  UNPROTECT(1);
  return r;
}

}  // extern "C"
