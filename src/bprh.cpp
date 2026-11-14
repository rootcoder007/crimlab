// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry point of prospective hot-spotting (prohs, src/prh.cpp): it
// checks the arguments, calls the kernel and builds the R result. Nothing
// calls back into the R evaluator.
#include <cmath>
#include <vector>

#include "bcol.h"
#include "rchk.h"
#include "spin.h"

using clab::fail;

extern "C" {

SEXP C_clab_r_prohs(SEXP x, SEXP y, SEXP t, SEXP at, SEXP win, SEXP cell, SEXP bw, SEXP tw, SEXP tunit,
                    SEXP top, SEXP ahead) {
  clab::who = "prohs";
  csp::pts p = csp::xy(x, y, R_NilValue, csp::planar);
  const int n = p.n;
  SEXP tv = PROTECT(bj::nums(t, "t", n));
  const double a0 = clab::one(at, "at", -HUGE_VAL, HUGE_VAL);
  double w[4];
  csp::win(win, w);
  const double ce = clab::one(cell, "cell", 0.0, HUGE_VAL);
  const double b = Rf_isNull(bw) ? 8 * ce : clab::one(bw, "bw", 0.0, HUGE_VAL);
  if (b / ce > 1e4) fail("`bw` must be at most 10000 cells.");
  const double tws = clab::one(tw, "tw", 0.0, HUGE_VAL);
  const double tu = clab::one(tunit, "tunit", 0.0, HUGE_VAL);
  const double tp = bj::inr(top, "top", 0.0, 1.0);
  if (!(tp > 0)) fail("`top` must be above 0.");
  const double ah = bj::inr(ahead, "ahead", 0.0, HUGE_VAL);
  const double gx = std::ceil((w[1] - w[0]) / ce), gy = std::ceil((w[3] - w[2]) / ce);
  if (gx * gy > 1e8) fail("the grid would have more than 1e8 cells: make `cell` larger.");
  const int nx = static_cast<int>(gx), ny = static_cast<int>(gy);
  SEXP z = PROTECT(Rf_allocMatrix(REALSXP, nx, ny));
  SEXP hot = PROTECT(Rf_allocMatrix(LGLSXP, nx, ny));
  double ev[5];
  const int rc = clab_prohs(p.x, p.y, REAL(tv), n, a0, w, ce, b, tws, tu, tp, ah, nx, ny, REAL(z), LOGICAL(hot), ev);
  clab::ichk(rc);
  bj::kc(rc, "the kernel refused the input.");
  SEXP gxv = PROTECT(Rf_allocVector(REALSXP, nx)), gyv = PROTECT(Rf_allocVector(REALSXP, ny));
  for (int i = 0; i < nx; ++i) REAL(gxv)[i] = w[0] + (i + 0.5) * ce;
  for (int j = 0; j < ny; ++j) REAL(gyv)[j] = w[2] + (j + 0.5) * ce;
  for (double &e : ev) e = bj::na(e);
  static const char *const en[] = {"n", "hits", "hit", "area", "pai"};
  const double par[7] = {a0, ce, b, tws, tu, tp, ah};
  static const char *const pn[] = {"at", "cell", "bw", "tw", "tunit", "top", "ahead"};
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 6));
  SET_VECTOR_ELT(out, 0, gxv);
  SET_VECTOR_ELT(out, 1, gyv);
  SET_VECTOR_ELT(out, 2, z);
  SET_VECTOR_ELT(out, 3, hot);
  SET_VECTOR_ELT(out, 4, clab::nvec(5, en, ev));
  SET_VECTOR_ELT(out, 5, clab::nvec(7, pn, par));
  static const char *const on[] = {"x", "y", "z", "hot", "eval", "par"};
  bj::named(out, on);
  Rf_classgets(out, PROTECT(Rf_mkString("prohs")));
  UNPROTECT(8);
  return out;
}

}  // extern "C"
