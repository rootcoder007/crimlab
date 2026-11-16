// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the air pollution family (src/apol.cpp). The R
// wrappers (R/apol.R, R/ahlth.R) check the arguments; these call the
// kernels and build the results. Nothing calls back into the R evaluator.
#include <cmath>
#include <vector>

#include "bcol.h"
#include "rchk.h"

namespace {

int n32(SEXP x) { return clab::ilen(static_cast<double>(XLENGTH(x)), nullptr); }

void kc(int rc) {
  clab::ichk(rc);
  if (rc != CLAB_OK) clab::fail("%s", clab_strerror(rc));
}

SEXP vec(int n) { return Rf_allocVector(REALSXP, n); }

}  // namespace

extern "C" {

SEXP C_clab_r_pgsig(SEXP x, SEXP cls, SEXP urb) {
  clab::who = "pgsig";
  const int n = n32(x);
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(r, 0, vec(n));
  SET_VECTOR_ELT(r, 1, vec(n));
  kc(clab_pgsig(REAL(x), n, Rf_asInteger(cls), Rf_asInteger(urb), REAL(VECTOR_ELT(r, 0)), REAL(VECTOR_ELT(r, 1))));
  UNPROTECT(1);
  return r;
}

SEXP C_clab_r_briggs(SEXP x, SEXP u, SEXP diam, SEXP ev, SEXP ts, SEXP ta, SEXP cls, SEXP dz, SEXP g) {
  clab::who = "briggs";
  const int n = n32(x);
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(r, 0, vec(3));
  SET_VECTOR_ELT(r, 1, vec(n));
  kc(clab_briggs(REAL(x), n, Rf_asReal(u), Rf_asReal(diam), Rf_asReal(ev), Rf_asReal(ts), Rf_asReal(ta),
                 Rf_asInteger(cls), Rf_asReal(dz), Rf_asReal(g), REAL(VECTOR_ELT(r, 0)), REAL(VECTOR_ELT(r, 1))));
  UNPROTECT(1);
  return r;
}

SEXP C_clab_r_gplume(SEXP q, SEXP u, SEXP h, SEXP R, SEXP cls, SEXP urb, SEXP lid, SEXP nimg) {
  clab::who = "gplume";
  const int n = Rf_nrows(R);
  SEXP c = PROTECT(vec(n));
  kc(clab_gplume(Rf_asReal(q), Rf_asReal(u), Rf_asReal(h), REAL(R), n, Rf_asInteger(cls), Rf_asInteger(urb),
                 Rf_asReal(lid), Rf_asInteger(nimg), REAL(c)));
  UNPROTECT(1);
  return c;
}

SEXP C_clab_r_gspuff(SEXP m, SEXP u, SEXP h, SEXP R, SEXP t, SEXP cls, SEXP urb, SEXP lid, SEXP nimg) {
  clab::who = "gspuff";
  const int n = Rf_nrows(R);
  SEXP c = PROTECT(vec(n));
  kc(clab_gspuff(Rf_asReal(m), Rf_asReal(u), Rf_asReal(h), REAL(R), n, Rf_asReal(t), Rf_asInteger(cls),
                 Rf_asInteger(urb), Rf_asReal(lid), Rf_asInteger(nimg), REAL(c)));
  UNPROTECT(1);
  return c;
}

SEXP C_clab_r_advdif(SEXP c0, SEXP u, SEXP v, SEXP kx, SEXP ky, SEXP dx, SEXP dy, SEXP dt, SEXP ns, SEXP src) {
  clab::who = "advdif2d";
  const int nx = Rf_nrows(c0), ny = Rf_ncols(c0);
  SEXP cm = PROTECT(Rf_allocMatrix(REALSXP, nx, ny));
  for (R_xlen_t i = 0; i < clab::xlen(c0, "c0"); ++i) REAL(cm)[i] = REAL(c0)[i];
  kc(clab_advdif(REAL(cm), nx, ny, Rf_asReal(u), Rf_asReal(v), Rf_asReal(kx), Rf_asReal(ky), Rf_asReal(dx),
                 Rf_asReal(dy), Rf_asReal(dt), Rf_asReal(ns), Rf_isNull(src) ? nullptr : REAL(src)));
  UNPROTECT(1);
  return cm;
}

// The positions' means and, with a grid, the density on it.
SEXP C_clab_r_laggrid(SEXP xs, SEXP ys, SEXP g) {
  clab::who = "lagrprt";
  const int n = n32(xs);
  const bool hg = !Rf_isNull(g);
  SEXP r = PROTECT(Rf_allocVector(VECSXP, hg ? 2 : 1));
  SET_VECTOR_ELT(r, 0, vec(2));
  double *conc = nullptr;
  if (hg) {
    SET_VECTOR_ELT(r, 1, Rf_allocMatrix(REALSXP, static_cast<int>(REAL(g)[4]), static_cast<int>(REAL(g)[5])));
    conc = REAL(VECTOR_ELT(r, 1));
  }
  kc(clab_laggrid(REAL(xs), REAL(ys), n, hg ? REAL(g) : nullptr, conc, REAL(VECTOR_ELT(r, 0))));
  UNPROTECT(1);
  return r;
}

SEXP C_clab_r_lagr(SEXP n, SEXP x0, SEXP y0, SEXP u, SEXP v, SEXP kx, SEXP ky, SEXP dt, SEXP ns, SEXP seed,
                   SEXP g) {
  clab::who = "lagrprt";
  const int m = Rf_asInteger(n);
  SEXP xs = PROTECT(vec(m)), ys = PROTECT(vec(m));
  kc(clab_lagr(m, Rf_asReal(x0), Rf_asReal(y0), Rf_asReal(u), Rf_asReal(v), Rf_asReal(kx), Rf_asReal(ky),
               Rf_asReal(dt), Rf_asReal(ns), Rf_asReal(seed), REAL(xs), REAL(ys)));
  SEXP gr = PROTECT(C_clab_r_laggrid(xs, ys, g));
  const int k = Rf_length(gr);
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 2 + k));
  SET_VECTOR_ELT(r, 0, xs);
  SET_VECTOR_ELT(r, 1, ys);
  for (int i = 0; i < k; ++i) SET_VECTOR_ELT(r, 2 + i, VECTOR_ELT(gr, i));
  UNPROTECT(4);
  return r;
}

SEXP C_clab_r_lagstep(SEXP xs, SEXP ys, SEXP u, SEXP v, SEXP fx, SEXP fy, SEXP dt, SEXP zx, SEXP zy) {
  clab::who = "lagrprt";
  const int n = n32(xs);
  SEXP a = PROTECT(Rf_duplicate(xs)), b = PROTECT(Rf_duplicate(ys));
  kc(clab_lagstep(n, Rf_asReal(u), Rf_asReal(v), Rf_asReal(fx), Rf_asReal(fy), Rf_asReal(dt), REAL(zx), REAL(zy),
                  REAL(a), REAL(b)));
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(r, 0, a);
  SET_VECTOR_ELT(r, 1, b);
  UNPROTECT(3);
  return r;
}

SEXP C_clab_r_crf(SEXP z, SEXP refc, SEXP form, SEXP par) {
  clab::who = "crf";
  const int n = n32(z);
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(r, 0, vec(n));
  SET_VECTOR_ELT(r, 1, vec(n));
  kc(clab_crf(REAL(z), n, Rf_asReal(refc), Rf_asInteger(form), REAL(par), REAL(VECTOR_ELT(r, 0)),
              REAL(VECTOR_ELT(r, 1))));
  UNPROTECT(1);
  return r;
}

SEXP C_clab_r_levin(SEXP rr, SEXP p) {
  const int n = n32(rr);
  SEXP r = PROTECT(vec(n));
  for (int i = 0; i < n; ++i) REAL(r)[i] = clab_levin(REAL(rr)[i], Rf_asReal(p));
  UNPROTECT(1);
  return r;
}

SEXP C_clab_r_hif(SEXP dx, SEXP pop, SEXP rate, SEXP beta) {
  return Rf_ScalarReal(clab_hif(Rf_asReal(dx), Rf_asReal(pop), Rf_asReal(rate), Rf_asReal(beta)));
}

SEXP C_clab_r_concidx(SEXP h, SEXP inc) {
  clab::who = "expconcidx";
  const int n = n32(h);
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(r, 0, vec(3));
  SET_VECTOR_ELT(r, 1, vec(n));
  SET_VECTOR_ELT(r, 2, vec(n));
  kc(clab_concidx(REAL(h), REAL(inc), n, REAL(VECTOR_ELT(r, 0)), REAL(VECTOR_ELT(r, 1)), REAL(VECTOR_ELT(r, 2))));
  UNPROTECT(1);
  return r;
}

}  // extern "C"
