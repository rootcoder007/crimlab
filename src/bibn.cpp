// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry point of the internal benchmarking of officers (ibench,
// src/ibn.cpp): it checks the arguments, calls the kernel and builds the
// R result. Nothing calls back into the R evaluator.
#include <string>
#include <vector>

#include "bcol.h"
#include "rchk.h"

using clab::fail;

namespace {

// The stops' circumstances: NULL, a numeric matrix, a data frame or list
// of numeric columns, or one numeric vector; n rows, finite.
int xcols(SEXP x, int n, std::vector<double> &v) {
  v.clear();
  if (Rf_isNull(x)) return 0;
  if (Rf_isMatrix(x)) {
    if (!clab::isnum(x) || Rf_nrows(x) != n) fail("`x` must be a numeric matrix with one row per stop, %d.", n);
    SEXP xd = PROTECT(clab::dbl(x));
    v.assign(REAL(xd), REAL(xd) + XLENGTH(xd));
    UNPROTECT(1);
    for (double e : v) {
      if (!std::isfinite(e)) fail("`x` must hold finite numbers with no missing values.");
    }
    return Rf_ncols(x);
  }
  std::vector<SEXP> cs;
  if (TYPEOF(x) == VECSXP) {
    for (R_xlen_t j = 0; j < XLENGTH(x); ++j) cs.push_back(VECTOR_ELT(x, j));
  } else {
    cs.push_back(x);
  }
  for (SEXP c : cs) {
    SEXP cv = PROTECT(clab::num(c, "x", n));
    v.insert(v.end(), REAL(cv), REAL(cv) + n);
    UNPROTECT(1);
  }
  return static_cast<int>(cs.size());
}

}  // namespace

extern "C" {

SEXP C_clab_r_ibench(SEXP y, SEXP off, SEXP x, SEXP minn, SEXP fdr) {
  clab::who = "ibench";
  int np = 0;
  SEXP yv = PROTECT(csp::bin01(y, "y", XLENGTH(y)));
  ++np;
  const int n = bj::len(yv, "y");
  if (n < 2) fail("`y` must hold two stops at least.");
  csp::lab o = bj::labs(off, "off", n);
  np += 2;
  if (o.g < 2) fail("`off` must hold at least two officers.");
  std::vector<double> xv;
  const int p = xcols(x, n, xv);
  const int mn = clab::whole(minn, "minn", 1);
  const double fd = clab::one(fdr, "fdr", 0.0, 1.0);
  std::vector<double> yd(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) yd[static_cast<std::size_t>(i)] = INTEGER(yv)[i];
  const int no = o.g;
  std::vector<double> res(7 * static_cast<std::size_t>(no));
  std::vector<int> ok(static_cast<std::size_t>(no));
  const int rc = clab_ibench(yd.data(), INTEGER(o.code), p ? xv.data() : nullptr, n, p, no, mn, res.data(), ok.data());
  clab::ichk(rc);
  bj::kc(rc, "the kernel refused the input.");
  bj::cols t(10, no);
  t.put("off", o.lev);
  static const char *const cn[] = {"n", "rate", "bench", "ess", "z", "p", "q"};
  for (int k = 0; k < 7; ++k) {
    double *c = t.dbl(cn[k]);
    for (int j = 0; j < no; ++j) c[j] = bj::na(res[static_cast<std::size_t>(k) * no + j]);
  }
  int *fl = t.lgl("flag");
  for (int j = 0; j < no; ++j) {
    const double q = res[6 * static_cast<std::size_t>(no) + j];
    fl[j] = std::isnan(q) ? NA_LOGICAL : q < fd;
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  ++np;
  SET_VECTOR_ELT(out, 0, t.done());
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(fd));
  SET_VECTOR_ELT(out, 2, Rf_ScalarInteger(mn));
  SET_VECTOR_ELT(out, 3, Rf_ScalarInteger(p));
  static const char *const on[] = {"tab", "fdr", "minn", "nx"};
  bj::named(out, on);
  Rf_classgets(out, PROTECT(Rf_mkString("ibench")));
  UNPROTECT(np + 2);
  return out;
}

}  // extern "C"
