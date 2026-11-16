// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the map layouts (src/map.cpp) and of the kernels
// added beside them, each named above its entry. The R wrappers check the
// arguments and pass doubles, whole numbers and codes from 0; these call
// the kernels and return their arrays. Nothing calls back into the R
// evaluator.
#include <algorithm>

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

// src/map.cpp: the breaks (br) and the class of each value (cls, 0 for
// none); meth 0 equal widths, 1 quantiles, 2 the breaks br as given.
extern "C" {

SEXP C_clab_r_mapbrk(SEXP fn, SEXP v, SEXP k, SEXP meth, SEXP br) {
  clab::who = CHAR(STRING_ELT(fn, 0));
  const int n = Rf_length(v), kk = Rf_asInteger(k), m = Rf_asInteger(meth);
  int nb = m == 2 ? Rf_length(br) : kk + 1;
  SEXP b = PROTECT(Rf_allocVector(REALSXP, nb));
  if (m == 2) std::copy(REAL(br), REAL(br) + nb, REAL(b));
  SEXP v2[2];
  v2[1] = PROTECT(Rf_allocVector(INTSXP, n));
  kc(clab_mapbrk(REAL(v), n, kk, m, REAL(b), &nb, INTEGER(v2[1])), "`v` has no finite value.");
  v2[0] = PROTECT(Rf_allocVector(REALSXP, nb));
  std::copy(REAL(b), REAL(b) + nb, REAL(v2[0]));
  static const char *const nms[] = {"br", "cls"};
  SEXP out = lst(2, nms, v2);
  UNPROTECT(3);
  return out;
}

// src/map.cpp: the box of the points (xmin, xmax, ymin, ymax).
SEXP C_clab_r_mapbb(SEXP fn, SEXP x, SEXP y) {
  clab::who = CHAR(STRING_ELT(fn, 0));
  SEXP bb = PROTECT(Rf_allocVector(REALSXP, 4));
  kc(clab_mapbb(REAL(x), REAL(y), Rf_length(x), REAL(bb)), "no point has both coordinates.");
  UNPROTECT(1);
  return bb;
}

}  // extern "C"

// src/map.cpp: the centres (x, y) of the clusters 1 to k and their
// points (n); cl the cluster of each point, 0 for none.
extern "C" {

SEXP C_clab_r_ptcen(SEXP x, SEXP y, SEXP cl, SEXP k) {
  clab::who = "ptmap";
  const int n = Rf_length(x), kk = Rf_asInteger(k);
  SEXP v[3];
  v[0] = PROTECT(Rf_allocVector(REALSXP, kk));
  v[1] = PROTECT(Rf_allocVector(REALSXP, kk));
  v[2] = PROTECT(Rf_allocVector(INTSXP, kk));
  kc(clab_ptcen(REAL(x), REAL(y), INTEGER(cl), n, kk, REAL(v[0]), REAL(v[1]), INTEGER(v[2])));
  static const char *const nms[] = {"x", "y", "n"};
  SEXP out = lst(3, nms, v);
  UNPROTECT(3);
  return out;
}

}  // extern "C"

// src/map.cpp: the box widened by the largest radius (bb), the radii
// (r), the drawing order from 1 (ord) and the largest radius (rmax, from
// the box when NA).
extern "C" {

SEXP C_clab_r_mapsym(SEXP x, SEXP y, SEXP v, SEXP rmax) {
  clab::who = "prpmap";
  const int n = Rf_length(x);
  SEXP o[4];
  o[0] = PROTECT(Rf_allocVector(REALSXP, 4));
  o[1] = PROTECT(Rf_allocVector(REALSXP, n));
  SEXP od = PROTECT(Rf_allocVector(INTSXP, n));
  int no = 0;
  double rm = 0.0;
  kc(clab_mapsym(REAL(x), REAL(y), REAL(v), n, Rf_asReal(rmax), REAL(o[0]), REAL(o[1]), INTEGER(od), &no, &rm),
     "`v` must have a value.");
  o[2] = PROTECT(Rf_allocVector(INTSXP, no));
  for (int q = 0; q < no; ++q) INTEGER(o[2])[q] = INTEGER(od)[q] + 1;
  o[3] = PROTECT(Rf_ScalarReal(rm));
  static const char *const nms[] = {"bb", "r", "ord", "rmax"};
  SEXP out = lst(4, nms, o);
  UNPROTECT(5);
  return out;
}

}  // extern "C"

// src/eof.cpp: the column means (center), the EOFs (eofs, p x k), the
// series (pcs, m x k), the squared singular values (d2) and the total sum
// of squares (tot) of z (m x p); k 0 for every component.
extern "C" {

SEXP C_clab_r_eofst(SEXP z, SEXP k) {
  clab::who = "eofst";
  const int m = Rf_nrows(z), p = Rf_ncols(z), s = std::min(m, p);
  SEXP ctr = PROTECT(Rf_allocVector(REALSXP, p));
  SEXP ev = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(p) * s));
  SEXP pc = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(m) * s));
  SEXP d2 = PROTECT(Rf_allocVector(REALSXP, s));
  double tot = 0.0;
  int kk = 0;
  kc(clab_eofst(REAL(z), m, p, Rf_asInteger(k), REAL(ctr), REAL(ev), REAL(pc), REAL(d2), &tot, &kk),
     "`z` has no variation: every column is constant.");
  SEXP v[5];
  v[0] = ctr;
  v[1] = PROTECT(Rf_allocMatrix(REALSXP, p, kk));
  std::copy(REAL(ev), REAL(ev) + static_cast<R_xlen_t>(p) * kk, REAL(v[1]));
  v[2] = PROTECT(Rf_allocMatrix(REALSXP, m, kk));
  std::copy(REAL(pc), REAL(pc) + static_cast<R_xlen_t>(m) * kk, REAL(v[2]));
  v[3] = PROTECT(Rf_allocVector(REALSXP, kk));
  std::copy(REAL(d2), REAL(d2) + kk, REAL(v[3]));
  v[4] = PROTECT(Rf_ScalarReal(tot));
  static const char *const nms[] = {"center", "eofs", "pcs", "d2", "tot"};
  SEXP out = lst(5, nms, v);
  UNPROTECT(8);
  return out;
}

}  // extern "C"

// src/tps.cpp: the basis coefficients (a), the plane (b: const, x, y),
// the surface at the points (fit) and on the grid gx by gy (z).
extern "C" {

SEXP C_clab_r_tpsurf(SEXP x, SEXP y, SEXP z, SEXP lam, SEXP gx, SEXP gy) {
  clab::who = "tpsurf";
  const int n = Rf_length(x), nx = Rf_length(gx), ny = Rf_length(gy);
  SEXP v[4];
  v[0] = PROTECT(Rf_allocVector(REALSXP, n));
  v[1] = PROTECT(Rf_allocVector(REALSXP, 3));
  v[2] = PROTECT(Rf_allocVector(REALSXP, n));
  v[3] = PROTECT(Rf_allocMatrix(REALSXP, nx, ny));
  kc(clab_tpsurf(REAL(x), REAL(y), REAL(z), n, Rf_asReal(lam), REAL(gx), nx, REAL(gy), ny, REAL(v[0]), REAL(v[1]),
                 REAL(v[2]), REAL(v[3])),
     "the system is singular: the points lie on one line, or repeat with `lam` 0.");
  static const char *const nms[] = {"a", "b", "fit", "z"};
  SEXP out = lst(4, nms, v);
  UNPROTECT(4);
  return out;
}

}  // extern "C"

// src/dag.cpp: the colliders a -> c <- b (codes from 0) and whether each
// is shielded (sh), and the distinct edges (nedge); nn nodes.
extern "C" {

SEXP C_clab_r_dagcol(SEXP from, SEXP to, SEXP nn) {
  clab::who = "dagcol";
  const int ne = Rf_length(from), n = Rf_asInteger(nn);
  int nc = 0, nd = 0;
  kc(clab_dagcol(INTEGER(from), INTEGER(to), ne, n, 0, nullptr, nullptr, nullptr, nullptr, &nc, &nd),
     "the edges have a cycle; a DAG has none.");
  SEXP v[5];
  v[0] = PROTECT(Rf_allocVector(INTSXP, nc));
  v[1] = PROTECT(Rf_allocVector(INTSXP, nc));
  v[2] = PROTECT(Rf_allocVector(INTSXP, nc));
  v[3] = PROTECT(Rf_allocVector(LGLSXP, nc));
  kc(clab_dagcol(INTEGER(from), INTEGER(to), ne, n, nc, INTEGER(v[0]), INTEGER(v[1]), INTEGER(v[2]), LOGICAL(v[3]),
                 &nc, &nd));
  v[4] = PROTECT(Rf_ScalarInteger(nd));
  static const char *const nms[] = {"a", "c", "b", "sh", "nedge"};
  SEXP out = lst(5, nms, v);
  UNPROTECT(5);
  return out;
}

}  // extern "C"
