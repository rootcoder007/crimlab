// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of idwint(), gtwreg(), asrule(), ifesc() and
// cemtch() (src/idw.cpp, src/gtw.cpp, src/arl.cpp, src/ifx.cpp,
// src/cem.cpp). The R wrappers (R/spi.R) check the
// arguments and pass doubles, codes from 0 and NULL for what is left out;
// these call the kernels and return their arrays. Nothing calls back
// into the R evaluator.
#include "crimlab.h"
#include "kutil.h"
#include "rchk.h"

namespace {

void kc(int rc, const char *dom = nullptr) {
  clab::ichk(rc);
  if (rc == CLAB_EDOM && dom) clab::fail("%s", dom);
  if (rc != CLAB_OK) clab::fail("%s", clab_strerror(rc));
}

SEXP lst(int n, const char *const *nm, SEXP *v) {
  SEXP out = PROTECT(Rf_allocVector(VECSXP, n));
  SEXP s = PROTECT(Rf_allocVector(STRSXP, n));
  for (int i = 0; i < n; ++i) {
    SET_VECTOR_ELT(out, i, v[i]);
    SET_STRING_ELT(s, i, Rf_mkChar(nm[i]));
  }
  Rf_setAttrib(out, R_NamesSymbol, s);
  UNPROTECT(2);
  return out;
}

}  // namespace

extern "C" {

SEXP C_clab_r_idwint(SEXP mx, SEXP my, SEXP mv, SEXP tx, SEXP ty, SEXP p, SEXP k, SEXP rad,
                     SEXP ll, SEXP zone, SEXP nz) {
  clab::who = "idwint";
  const int nm = clab::xlen(mx, "mx"), nt = clab::xlen(tx, "tx");
  const int g = Rf_isNull(zone) ? 0 : Rf_asInteger(nz);
  SEXP v[4];
  v[0] = PROTECT(Rf_allocVector(REALSXP, nt));
  v[1] = PROTECT(Rf_allocVector(REALSXP, nt));
  v[2] = PROTECT(Rf_allocVector(REALSXP, g));
  v[3] = PROTECT(Rf_allocVector(REALSXP, g));
  kc(clab_idwint(nm, REAL(mx), REAL(my), REAL(mv), nt, REAL(tx), REAL(ty), Rf_asReal(p),
                 Rf_asInteger(k), Rf_asReal(rad), Rf_asLogical(ll),
                 Rf_isNull(zone) ? nullptr : INTEGER(zone), g, REAL(v[0]), REAL(v[1]),
                 g ? REAL(v[2]) : nullptr, g ? REAL(v[3]) : nullptr));
  static const char *const nms[] = {"val", "used", "zval", "zn"};
  SEXP out = lst(4, nms, v);
  UNPROTECT(4);
  return out;
}

SEXP C_clab_r_gtwreg(SEXP x, SEXP y, SEXP u, SEXP v, SEXP t, SEXP h, SEXP tau, SEXP ll) {
  clab::who = "gtwreg";
  const int n = clab::xlen(y, "y"), k = Rf_ncols(x);
  SEXP r[6];
  r[0] = PROTECT(Rf_allocMatrix(REALSXP, n, k));
  r[1] = PROTECT(Rf_allocMatrix(REALSXP, n, k));
  r[2] = PROTECT(Rf_allocVector(REALSXP, n));
  r[3] = PROTECT(Rf_allocVector(REALSXP, n));
  r[4] = PROTECT(Rf_allocVector(REALSXP, n));
  r[5] = PROTECT(Rf_allocVector(REALSXP, 7));
  kc(clab_gtwreg(n, k, REAL(x), REAL(y), REAL(u), REAL(v), Rf_isNull(t) ? nullptr : REAL(t),
                 Rf_asReal(h), Rf_asReal(tau), Rf_asLogical(ll), REAL(r[0]), REAL(r[1]),
                 REAL(r[2]), REAL(r[3]), REAL(r[4]), REAL(r[5])),
     "a local fit is singular: widen the bandwidth or drop a regressor.");
  static const char *const sn[] = {"rss", "trs", "trss", "sigma2", "cv", "aicc", "edf"};
  SEXP s = PROTECT(Rf_allocVector(STRSXP, 7));
  for (int i = 0; i < 7; ++i) SET_STRING_ELT(s, i, Rf_mkChar(sn[i]));
  Rf_setAttrib(r[5], R_NamesSymbol, s);
  static const char *const nms[] = {"coef", "se", "fit", "hat", "loo", "stat"};
  SEXP out = lst(6, nms, r);
  UNPROTECT(7);
  return out;
}

SEXP C_clab_r_asrule(SEXP xc, SEXP nx, SEXP yc, SEXP ny, SEXP minsup, SEXP minconf,
                     SEXP empty) {
  clab::who = "asrule";
  const int n = clab::xlen(xc, "x"), a = Rf_asInteger(nx), b = Rf_asInteger(ny);
  const double mx = (static_cast<double>(a) + b) + 2.0 * a * b;
  if (mx > kutil::maxcell) clab::fail("the levels of `x` and `y` are too many.");
  const int m = static_cast<int>(mx);
  int nr = 0;
  int *lh = reinterpret_cast<int *>(R_alloc(static_cast<size_t>(m), sizeof(int)));
  int *rh = reinterpret_cast<int *>(R_alloc(static_cast<size_t>(m), sizeof(int)));
  double *o = reinterpret_cast<double *>(R_alloc(static_cast<size_t>(m) * 5, sizeof(double)));
  kc(clab_asrule(n, INTEGER(xc), a, INTEGER(yc), b, Rf_asReal(minsup), Rf_asReal(minconf),
                 Rf_asLogical(empty), &nr, lh, rh, o));
  SEXP v[3];
  v[0] = PROTECT(Rf_allocVector(INTSXP, nr));
  v[1] = PROTECT(Rf_allocVector(INTSXP, nr));
  v[2] = PROTECT(Rf_allocMatrix(REALSXP, nr, 5));
  for (int i = 0; i < nr; ++i) {
    INTEGER(v[0])[i] = lh[i];
    INTEGER(v[1])[i] = rh[i];
    for (int j = 0; j < 5; ++j) REAL(v[2])[static_cast<size_t>(j) * nr + i] = o[i * 5 + j];
  }
  static const char *const nms[] = {"lhs", "rhs", "out"};
  SEXP out = lst(3, nms, v);
  UNPROTECT(3);
  return out;
}

SEXP C_clab_r_ifesc(SEXP y, SEXP x, SEXP a, SEXP r, SEXP rmax, SEXP maxit, SEXP tol, SEXP bm) {
  clab::who = "ifesc";
  const int nt = Rf_nrows(y), n = Rf_ncols(y);
  SEXP dm = Rf_getAttrib(x, R_DimSymbol);
  const int p = Rf_length(dm) == 3 ? INTEGER(dm)[2] : 0;
  const int rr = Rf_asInteger(r), rx = Rf_asInteger(rmax);
  const int nb = Rf_isNull(bm) ? 0 : Rf_ncols(bm);
  const int nf = rr > rx ? rr : rx;
  if (static_cast<double>(nt) * n > kutil::maxcell || static_cast<double>(nb) * (2.0 * nt + 1.0) > kutil::maxcell) {
    clab::fail("the panel or the bootstrap is too large.");
  }
  const int ne = 2 * nt;
  SEXP v[7];
  v[0] = PROTECT(Rf_allocVector(REALSXP, p));
  v[1] = PROTECT(Rf_allocMatrix(REALSXP, nt, n));
  v[2] = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(nt) * (nf > 0 ? nf : 0)));
  v[3] = PROTECT(Rf_allocVector(REALSXP, rx + 1));
  v[4] = PROTECT(Rf_allocVector(REALSXP, ne + 1));
  v[5] = PROTECT(Rf_allocMatrix(REALSXP, nb, ne + 1));
  v[6] = PROTECT(Rf_allocVector(REALSXP, 8));
  kc(clab_ifesc(REAL(y), REAL(x), p, INTEGER(a), n, nt, rr, rx, Rf_asInteger(maxit), Rf_asReal(tol),
                Rf_isNull(bm) ? nullptr : INTEGER(bm), nb, REAL(v[0]), REAL(v[1]), REAL(v[2]), REAL(v[3]),
                REAL(v[4]), REAL(v[5]), REAL(v[6])),
     "a treated unit has too few periods before adoption for the factors, or a fit is singular.");
  // the kernel's ne = nt - min(a) + max(a) < 2 nt: keep the first ne + 1
  const int k = static_cast<int>(REAL(v[6])[7]) + 1, r1 = static_cast<int>(REAL(v[6])[0]);
  SEXP e1 = PROTECT(Rf_allocVector(REALSXP, k));
  for (int i = 0; i < k; ++i) REAL(e1)[i] = REAL(v[4])[i];
  v[4] = e1;
  SEXP be = PROTECT(Rf_allocMatrix(REALSXP, nb, k));
  for (R_xlen_t i = 0; i < static_cast<R_xlen_t>(nb) * k; ++i) REAL(be)[i] = REAL(v[5])[i];
  v[5] = be;
  SEXP fc = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(nt) * r1));
  for (R_xlen_t i = 0; i < static_cast<R_xlen_t>(nt) * r1; ++i) REAL(fc)[i] = REAL(v[2])[i];
  v[2] = fc;
  static const char *const nms[] = {"beta", "y0", "fac", "cv", "ef", "bef", "st"};
  SEXP out = lst(7, nms, v);
  UNPROTECT(10);
  return out;
}

SEXP C_clab_r_cemtch(SEXP code, SEXP tr, SEXP lcode) {
  clab::who = "cemtch";
  const int n = Rf_nrows(code), k = Rf_ncols(code), lk = Rf_ncols(lcode);
  SEXP v[3];
  v[0] = PROTECT(Rf_allocVector(INTSXP, n));
  v[1] = PROTECT(Rf_allocVector(REALSXP, n));
  v[2] = PROTECT(Rf_allocVector(REALSXP, 6));
  kc(clab_cemtch(INTEGER(code), n, k, INTEGER(tr), INTEGER(lcode), lk, INTEGER(v[0]), REAL(v[1]), REAL(v[2])));
  static const char *const nms[] = {"strata", "w", "st"};
  SEXP out = lst(3, nms, v);
  UNPROTECT(3);
  return out;
}

}  // extern "C"
