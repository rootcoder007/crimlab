// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the learning and causal families: each checks its
// arguments, calls one kernel with scratch from R's transient memory and
// builds the R result. Nothing here calls back into the R evaluator.
#include <climits>
#include <cmath>
#include <cstring>

#include "crimlab.h"
#include "rchk.h"

using clab::fail;
using clab::kchk;
using clab::num;

namespace rl {

// x as a double matrix: n rows, p columns, every value finite

SEXP mat(SEXP x, const char *nm, int *n, int *p) {
  // is.matrix(x) && is.numeric(x) && all(is.finite(x)), at least one row
  SEXP d = Rf_getAttrib(x, R_DimSymbol);
  if (!clab::isnum(x) || TYPEOF(d) != INTSXP || XLENGTH(d) != 2 ||
      INTEGER_RO(d)[0] < 1) {
    fail("`%s` must be a numeric matrix with at least one row.", nm);
  }
  SEXP v = PROTECT(clab::dbl(x));
  const double *a = REAL_RO(v);
  for (R_xlen_t i = 0; i < XLENGTH(v); ++i) {
    if (!std::isfinite(a[i])) {
      fail("`%s` must be finite numbers with no missing values.", nm);
    }
  }
  *n = INTEGER_RO(d)[0];
  *p = INTEGER_RO(d)[1];
  UNPROTECT(1);
  return v;
}

int *folds(SEXP f, const char *nm, int n, int *k) {
  // whole numbers 1..K, length n, every one of them used, K >= 2
  SEXP v = PROTECT(num(f, nm, n, 1.0));
  const double *a = REAL_RO(v);
  double mx = 1.0;
  for (int i = 0; i < n; ++i) {
    if (a[i] != std::floor(a[i])) fail("`%s` must be whole numbers.", nm);
    mx = std::fmax(mx, a[i]);
  }
  if (mx < 2.0 || mx > n) fail("`%s` must number at least two folds.", nm);
  const int kk = static_cast<int>(mx);
  int *o = reinterpret_cast<int *>(R_alloc(static_cast<std::size_t>(n), sizeof(int)));
  char *seen = R_alloc(static_cast<std::size_t>(kk), 1);
  std::memset(seen, 0, static_cast<std::size_t>(kk));
  for (int i = 0; i < n; ++i) {
    o[i] = static_cast<int>(a[i]) - 1;
    seen[o[i]] = 1;
  }
  for (int j = 0; j < kk; ++j) {
    if (!seen[j]) fail("`%s` leaves fold %d empty.", nm, j + 1);
  }
  *k = kk;
  UNPROTECT(1);
  return o;
}

int pick(SEXP s, const char *nm, const char *const *ch, int k) {
  if (TYPEOF(s) == STRSXP && XLENGTH(s) == 1 && STRING_ELT(s, 0) != NA_STRING) {
    const char *c = CHAR(STRING_ELT(s, 0));
    for (int i = 0; i < k; ++i) {
      if (std::strcmp(c, ch[i]) == 0) return i;
    }
  }
  // the choices, quoted and joined
  char buf[256];
  std::size_t len = 0;
  for (int i = 0; i < k && len + 40 < sizeof buf; ++i) {
    len += static_cast<std::size_t>(std::snprintf(
        buf + len, sizeof buf - len, "%s\"%s\"", i == 0 ? "" : ", ", ch[i]));
  }
  fail("`%s` must be one of %s.", nm, buf);
}

// scratch in R's transient memory; R_alloc() refuses a block too large
double *dalloc(double n) {
  return reinterpret_cast<double *>(R_alloc(static_cast<std::size_t>(n), sizeof(double)));
}

int *ialloc(double n) {
  return reinterpret_cast<int *>(R_alloc(static_cast<std::size_t>(n), sizeof(int)));
}

void ychk(const double *y, int n, bool bin, const char *nm) {
  bool one = false, zero = false, diff = false;
  for (int i = 0; i < n; ++i) {
    if (bin && y[i] != 0.0 && y[i] != 1.0) {
      fail("`%s` must be 0 or 1 for a binary outcome.", nm);
    }
    one = one || y[i] == 1.0;
    zero = zero || y[i] == 0.0;
    diff = diff || y[i] != y[0];
  }
  if (!diff || (bin && !(one && zero))) fail("`%s` is constant.", nm);
}

// a matrix of n rows and p columns, column names those of x (if any)
SEXP rmat(int n, int p, SEXP x) {
  SEXP m = PROTECT(Rf_allocMatrix(REALSXP, n, p));
  SEXP dn = Rf_getAttrib(x, R_DimNamesSymbol);
  if (!Rf_isNull(dn) && !Rf_isNull(VECTOR_ELT(dn, 1))) {
    SEXP d = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(d, 0, VECTOR_ELT(dn, 1));
    Rf_setAttrib(m, R_DimNamesSymbol, d);
    UNPROTECT(1);
  }
  UNPROTECT(1);
  return m;
}

}  // namespace rl

extern "C" {

SEXP C_clab_r_enpath(SEXP x, SEXP y, SEXP fam, SEXP alpha, SEXP lambda,
                    SEXP nlam, SEXP lmr, SEXP fold, SEXP tol, SEXP maxit) {
  clab::who = "enpath";
  int n, p;
  SEXP xv = PROTECT(rl::mat(x, "x", &n, &p));
  if (n < 2 || p < 1) fail("`x` must have at least two rows and one column.");
  SEXP yv = PROTECT(num(y, "y", n));
  static const char *const fams[] = {"gaussian", "binomial"};
  const int f = rl::pick(fam, "fam", fams, 2);
  rl::ychk(REAL_RO(yv), n, f == 1, "y");
  const double a = clab::p01(alpha, "alpha");
  int nl;
  double r = 0.0;
  SEXP lam;
  if (Rf_isNull(lambda)) {
    nl = clab::whole(nlam, "nlam", 1);
    r = Rf_isNull(lmr) ? (n > p ? 1e-4 : 1e-2) : clab::one(lmr, "lmr", 0, 1);
    lam = PROTECT(Rf_allocVector(REALSXP, nl));
  } else {
    lam = PROTECT(clab::dcopy(num(lambda, "lambda", -1, 0.0)));
    nl = clab::ilen(static_cast<double>(XLENGTH(lam)), "lambda");
    for (int k = 1; k < nl; ++k) {
      if (REAL_RO(lam)[k] > REAL_RO(lam)[k - 1]) {
        fail("`lambda` must be decreasing.");
      }
    }
  }
  int nf = 0;
  const int *fd = Rf_isNull(fold) ? nullptr : rl::folds(fold, "fold", n, &nf);
  const double th = REAL_RO(num(tol, "tol", 1, 0.0, true))[0];
  const int mi = clab::whole(maxit, "maxit", 1);
  SEXP a0 = PROTECT(Rf_allocVector(REALSXP, nl));
  SEXP beta = PROTECT(rl::rmat(p, nl, x));
  SEXP cv = PROTECT(Rf_allocVector(REALSXP, 2 * static_cast<R_xlen_t>(nl)));
  int best[2];
  const double dp = p, dn = n;
  double *w = rl::dalloc(dp * dp + 6 * dp + 4 * dn + nl * (dp + 1) + nf * (nl + 1.0));
  int *iw = rl::ialloc(2 * dn + dp);
  const int rc = clab_enpath(REAL_RO(xv), n, p, REAL_RO(yv), f, a, REAL(lam), nl,
                            r, th, mi, fd, nf, REAL(a0), REAL(beta), REAL(cv),
                            best, w, iw);
  if (rc == CLAB_EDOM) {
    fail("no fit: no convergence within `maxit` passes, or `y` constant "
         "once a fold is held out.");
  }
  kchk(rc, nullptr);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 7));
  SET_VECTOR_ELT(out, 0, lam);
  SET_VECTOR_ELT(out, 1, a0);
  SET_VECTOR_ELT(out, 2, beta);
  if (fd != nullptr) {
    SEXP cm = PROTECT(Rf_allocVector(REALSXP, nl));
    SEXP cs = PROTECT(Rf_allocVector(REALSXP, nl));
    std::memcpy(REAL(cm), REAL_RO(cv), sizeof(double) * nl);
    std::memcpy(REAL(cs), REAL_RO(cv) + nl, sizeof(double) * nl);
    SET_VECTOR_ELT(out, 3, cm);
    SET_VECTOR_ELT(out, 4, cs);
    SET_VECTOR_ELT(out, 5, Rf_ScalarReal(REAL_RO(lam)[best[0]]));
    SET_VECTOR_ELT(out, 6, Rf_ScalarReal(REAL_RO(lam)[best[1]]));
    UNPROTECT(2);
  }
  static const char *const nm[] = {"lambda", "a0",   "beta", "cvm",
                                   "cvsd",   "lmin", "l1se"};
  UNPROTECT(7);
  return clab::nlist(out, nm);
}

SEXP C_clab_r_lgtreg(SEXP x, SEXP y, SEXP tol, SEXP maxit) {
  clab::who = "lgtreg";
  int n, p;
  SEXP xv = PROTECT(rl::mat(x, "x", &n, &p));
  SEXP yv = PROTECT(num(y, "y", n));
  rl::ychk(REAL_RO(yv), n, true, "y");
  if (n <= p + 1) fail("`x` needs more rows than columns plus one.");
  const double th = REAL_RO(num(tol, "tol", 1, 0.0, true))[0];
  const int mi = clab::whole(maxit, "maxit", 1);
  const double q = p + 1.0;
  double *w = rl::dalloc(2 * q * q + 4 * q);
  int *iw = rl::ialloc(n);
  SEXP cf = PROTECT(Rf_allocVector(REALSXP, p + 1));
  SEXP se = PROTECT(Rf_allocVector(REALSXP, p + 1));
  double dev;
  const int rc = clab_lgtreg(REAL_RO(xv), n, p, REAL_RO(yv), th, mi, REAL(cf),
                             REAL(se), &dev, w, iw);
  if (rc == CLAB_EDOM) {
    fail("no fit: collinear columns, or separated outcomes (no convergence "
         "within `maxit` steps).");
  }
  kchk(rc, nullptr);
  // names: "(Intercept)" and the column names of x, or x1, x2, ...
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, p + 1));
  SET_STRING_ELT(nms, 0, Rf_mkChar("(Intercept)"));
  SEXP dn = Rf_getAttrib(x, R_DimNamesSymbol);
  SEXP cn = Rf_isNull(dn) ? R_NilValue : VECTOR_ELT(dn, 1);
  char b[32];
  for (int j = 0; j < p; ++j) {
    if (!Rf_isNull(cn)) {
      SET_STRING_ELT(nms, j + 1, STRING_ELT(cn, j));
    } else {
      std::snprintf(b, sizeof b, "x%d", j + 1);
      SET_STRING_ELT(nms, j + 1, Rf_mkChar(b));
    }
  }
  Rf_setAttrib(cf, R_NamesSymbol, nms);
  Rf_setAttrib(se, R_NamesSymbol, nms);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, cf);
  SET_VECTOR_ELT(out, 1, se);
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(dev));
  static const char *const nm[] = {"coef", "se", "dev"};
  UNPROTECT(6);
  return clab::nlist(out, nm);
}

}  // extern "C"
