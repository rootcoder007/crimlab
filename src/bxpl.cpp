// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of xpimp() and the explanation functions after it
// (src/xpl.cpp). The R wrappers (R/xpl.R) check the arguments, pass
// doubles and call the prediction function between the design and the
// summary; these call the kernels and return their arrays. Nothing calls
// back into the R evaluator.
#include <cmath>

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

unsigned int useed(SEXP s) { return static_cast<unsigned int>(std::fmod(Rf_asReal(s), 4294967296.0)); }

}  // namespace

extern "C" {

SEXP C_clab_r_xpperm(SEXP x, SEXP nrep, SEXP seed) {
  clab::who = "xpimp";
  const int n = Rf_nrows(x), p = Rf_ncols(x), k = Rf_asInteger(nrep);
  // the wrapper keeps n p nrep p within 2^26
  const int N = n * p * k;
  SEXP v[2];
  v[0] = PROTECT(Rf_allocMatrix(REALSXP, N, p));
  v[1] = PROTECT(Rf_allocVector(INTSXP, N));
  kc(clab_xpperm(REAL(x), n, p, k, useed(seed), REAL(v[0]), INTEGER(v[1])));
  int *pm = INTEGER(v[1]);
  for (int i = 0; i < N; ++i) pm[i] += 1;
  SEXP dm = PROTECT(Rf_allocVector(INTSXP, 3));
  INTEGER(dm)[0] = n;
  INTEGER(dm)[1] = k;
  INTEGER(dm)[2] = p;
  Rf_setAttrib(v[1], R_DimSymbol, dm);
  static const char *const nms[] = {"des", "perm"};
  SEXP out = lst(2, nms, v);
  UNPROTECT(3);
  return out;
}

SEXP C_clab_r_xpimp(SEXP y, SEXP f0, SEXP fp, SEXP p, SEXP nrep, SEXP loss, SEXP ratio) {
  clab::who = "xpimp";
  const int n = Rf_length(y), q = Rf_asInteger(p), k = Rf_asInteger(nrep);
  SEXP v[2];
  v[0] = PROTECT(Rf_allocMatrix(REALSXP, q, k));
  v[1] = PROTECT(Rf_allocVector(REALSXP, 1 + 2 * static_cast<R_xlen_t>(q)));
  kc(clab_xpimp(REAL(y), REAL(f0), REAL(fp), n, q, k, Rf_asInteger(loss), Rf_asLogical(ratio), REAL(v[0]),
                REAL(v[1])),
     "`ratio` needs a loss above 0 before permuting.");
  static const char *const nms[] = {"rep", "st"};
  SEXP out = lst(2, nms, v);
  UNPROTECT(2);
  return out;
}

SEXP C_clab_r_xpgrid(SEXP v, SEXP ng, SEXP type) {
  const int n = Rf_length(v), k = Rf_asInteger(ng);
  SEXP g = PROTECT(Rf_allocVector(REALSXP, k));
  int m = 0;
  kc(clab_xpgrid(REAL(v), n, k, Rf_asInteger(type), REAL(g), &m));
  SEXP out = PROTECT(Rf_lengthgets(g, m));
  UNPROTECT(2);
  return out;
}

SEXP C_clab_r_xpcpd(SEXP s, SEXP j, SEXP g) {
  const int ns = Rf_nrows(s), p = Rf_ncols(s), ng = Rf_length(g);
  // the wrapper keeps ns ng p within 2^26
  SEXP des = PROTECT(Rf_allocMatrix(REALSXP, ns * ng, p));
  kc(clab_xpcpd(REAL(s), ns, p, Rf_asInteger(j) - 1, REAL(g), ng, REAL(des)));
  UNPROTECT(1);
  return des;
}

SEXP C_clab_r_xppd(SEXP f, SEXP ns, SEXP ng) {
  const int k = Rf_asInteger(ng);
  SEXP pd = PROTECT(Rf_allocVector(REALSXP, k));
  kc(clab_xppd(REAL(f), Rf_asInteger(ns), k, REAL(pd)));
  UNPROTECT(1);
  return pd;
}

SEXP C_clab_r_xpaled(SEXP x, SEXP j, SEXP nbin) {
  clab::who = "xpale";
  const int n = Rf_nrows(x), p = Rf_ncols(x), K = Rf_asInteger(nbin);
  SEXP v[3];
  v[0] = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(K) + 1));
  v[1] = PROTECT(Rf_allocVector(INTSXP, n));
  v[2] = PROTECT(Rf_allocMatrix(REALSXP, 2 * n, p));
  int nz = 0;
  kc(clab_xpaled(REAL(x), n, p, Rf_asInteger(j) - 1, K, REAL(v[0]), &nz, INTEGER(v[1]), REAL(v[2])),
     "`v` names a column with one value: it has no local effects.");
  v[0] = Rf_lengthgets(v[0], nz);
  UNPROTECT(3);
  PROTECT(v[0]);
  PROTECT(v[1]);
  PROTECT(v[2]);
  static const char *const nms[] = {"z", "bin", "des"};
  SEXP out = lst(3, nms, v);
  UNPROTECT(3);
  return out;
}

SEXP C_clab_r_xpale(SEXP f, SEXP bin, SEXP nb) {
  clab::who = "xpale";
  const int n = Rf_length(bin), k = Rf_asInteger(nb);
  SEXP v[3];
  v[0] = PROTECT(Rf_allocVector(REALSXP, k));
  v[1] = PROTECT(Rf_allocVector(REALSXP, k));
  v[2] = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(k) + 1));
  kc(clab_xpale(REAL(f), n, INTEGER(bin), k, REAL(v[0]), REAL(v[1]), REAL(v[2])));
  static const char *const nms[] = {"eff", "cnt", "ale"};
  SEXP out = lst(3, nms, v);
  UNPROTECT(3);
  return out;
}

SEXP C_clab_r_xpshd(SEXP x, SEXP s, SEXP m, SEXP seed) {
  clab::who = "xpshap";
  const int n = Rf_nrows(x), ns = Rf_nrows(s), p = Rf_ncols(x), k = Rf_asInteger(m);
  // the wrapper keeps 2 ns p m p within 2^26
  SEXP des = PROTECT(Rf_allocMatrix(REALSXP, 2 * ns * p * k, p));
  kc(clab_xpshd(REAL(x), n, REAL(s), ns, p, k, useed(seed), REAL(des)));
  UNPROTECT(1);
  return des;
}

SEXP C_clab_r_xpsh(SEXP f, SEXP ns, SEXP p, SEXP m) {
  clab::who = "xpshap";
  const int a = Rf_asInteger(ns), q = Rf_asInteger(p);
  SEXP v[2];
  v[0] = PROTECT(Rf_allocMatrix(REALSXP, a, q));
  v[1] = PROTECT(Rf_allocMatrix(REALSXP, a, q));
  kc(clab_xpsh(REAL(f), a, q, Rf_asInteger(m), REAL(v[0]), REAL(v[1])));
  static const char *const nms[] = {"phi", "se"};
  SEXP out = lst(2, nms, v);
  UNPROTECT(2);
  return out;
}

}  // extern "C"
