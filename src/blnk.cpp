// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry point of crime linkage (lnkpr, src/lnk.cpp): it checks the
// arguments, calls the kernel and builds the R result. Nothing calls back
// into the R evaluator.
#include <algorithm>
#include <cmath>
#include <cstring>
#include <map>
#include <string>
#include <vector>

#include "bcol.h"
#include "rchk.h"
#include "spin.h"

extern "C" double Rf_pnorm5(double, double, double, int, int);

using clab::fail;

namespace {

// 0/1 behaviours: a logical or numeric matrix, or a data frame or list of
// columns, n rows; column-major into v. The number of columns.
int mocols(SEXP mo, int n, std::vector<int> &v) {
  v.clear();
  if (Rf_isNull(mo)) return 0;
  if (Rf_isMatrix(mo)) {
    if (Rf_nrows(mo) != n || !(TYPEOF(mo) == LGLSXP || clab::isnum(mo))) {
      fail("`mo` must be a 0/1 or logical matrix with one row per crime, %d.", n);
    }
    SEXP d = PROTECT(TYPEOF(mo) == LGLSXP ? Rf_coerceVector(mo, REALSXP) : clab::dbl(mo));
    for (R_xlen_t k = 0; k < XLENGTH(d); ++k) {
      const double e = REAL(d)[k];
      if (!(e == 0 || e == 1)) fail("`mo` must hold 0 or 1 (or FALSE or TRUE), with no missing values.");
      v.push_back(static_cast<int>(e));
    }
    UNPROTECT(1);
    return Rf_ncols(mo);
  }
  if (TYPEOF(mo) != VECSXP) fail("`mo` must be a matrix, a data frame or a list of 0/1 columns.");
  for (R_xlen_t j = 0; j < XLENGTH(mo); ++j) {
    SEXP c = PROTECT(csp::bin01(VECTOR_ELT(mo, j), "mo", n));
    v.insert(v.end(), INTEGER(c), INTEGER(c) + n);
    UNPROTECT(1);
  }
  return static_cast<int>(XLENGTH(mo));
}

}  // namespace

extern "C" {

SEXP C_clab_r_lnkpr(SEXP x, SEXP y, SEXP t, SEXP series, SEXP mo) {
  clab::who = "lnkpr";
  csp::pts p = csp::xy(x, y, R_NilValue, csp::planar);
  const int n = p.n;
  if (n < 2) fail("`x` must hold two crimes at least.");
  if (n > 3000) fail("at most 3000 crimes: their pairs are kept in memory.");
  SEXP tv = PROTECT(bj::nums(t, "t", n));
  clab::unlab(series, "series");
  clab::vec(series, "series");
  if (clab::xlen(series, "series") != n) fail("`series` must have one label per crime, %d.", n);
  SEXP sl = PROTECT(clab::labs(series));
  if (Rf_isNull(sl)) fail("`series` must be labels (numbers, text or a factor), NA for an unsolved crime.");
  std::vector<int> ser(static_cast<std::size_t>(n), -1);
  {
    std::map<std::string, int> at;
    for (int i = 0; i < n; ++i) {
      if (STRING_ELT(sl, i) != NA_STRING) at.emplace(bj::str(sl, i), 0);
    }
    int g = 0;
    for (auto &e : at) e.second = g++;
    for (int i = 0; i < n; ++i) {
      if (STRING_ELT(sl, i) != NA_STRING) ser[static_cast<std::size_t>(i)] = at[bj::str(sl, i)];
    }
  }
  std::vector<int> mv;
  const int q = mocols(mo, n, mv);
  const std::size_t m = static_cast<std::size_t>(n) * (n - 1) / 2;
  bj::cols f(7, static_cast<int>(m));
  int *pi = f.itg("i"), *pj = f.itg("j");
  double *di = f.dbl("dist"), *dtv = f.dbl("dt"), *ja = f.dbl("jac");
  int *lk = f.lgl("linked");
  double *pr = f.dbl("prob");
  const int p3 = q > 0 ? 3 : 2;
  double cf[4], se[4], au[4], dev = 0.0;
  const int rc = clab_lnkpr(p.x, p.y, REAL(tv), q ? mv.data() : nullptr, q, ser.data(), n, pi, pj, di, dtv, ja, lk,
                            pr, cf, se, au, &dev);
  clab::ichk(rc);
  if (rc != CLAB_OK) {
    fail("the solved crimes need linked and unlinked pairs that do not separate perfectly.");
  }
  for (std::size_t k = 0; k < m; ++k) {
    ++pi[k];
    ++pj[k];
    if (lk[k] < 0) lk[k] = NA_LOGICAL;
    ja[k] = bj::na(ja[k]);
  }
  static const char *const tn[] = {"(Intercept)", "dist", "dt", "jac"};
  bj::cols c(5, p3 + 1);
  SEXP term = PROTECT(Rf_allocVector(STRSXP, p3 + 1));
  for (int j = 0; j <= p3; ++j) SET_STRING_ELT(term, j, Rf_mkChar(tn[j]));
  c.put("term", term);
  double *es = c.dbl("est"), *ss = c.dbl("se"), *zz = c.dbl("z"), *pv = c.dbl("p");
  for (int j = 0; j <= p3; ++j) {
    es[j] = cf[j];
    ss[j] = se[j];
    zz[j] = cf[j] / se[j];
    pv[j] = 2 * Rf_pnorm5(-std::fabs(zz[j]), 0.0, 1.0, 1, 0);
  }
  for (double &e : au) e = bj::na(e);
  static const char *const an[] = {"dist", "dt", "jac", "model"};
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
  SET_VECTOR_ELT(out, 0, c.done());
  SET_VECTOR_ELT(out, 1, clab::nvec(4, an, au));
  SET_VECTOR_ELT(out, 2, f.done());
  SET_VECTOR_ELT(out, 3, Rf_ScalarReal(dev));
  SET_VECTOR_ELT(out, 4, Rf_ScalarInteger(q));
  static const char *const on[] = {"coef", "auc", "pairs", "dev", "nmo"};
  bj::named(out, on);
  Rf_classgets(out, PROTECT(Rf_mkString("lnkpr")));
  UNPROTECT(8);
  return out;
}

}  // extern "C"
