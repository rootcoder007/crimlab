// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of deterrence theory (detrsp, dtrid, inspg, insphs,
// colarr): each checks its arguments, calls the kernel and builds the R
// result. Nothing calls back into the R evaluator.
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

#include "bfx.h"
#include "rchk.h"

using clab::fail;

namespace {

// Three numbers, one per strategy.
SEXP three(SEXP x, const char *nm) { return bj::nums(x, nm, 3, 0.0); }

// Names "1", "2", ... or the values of x as R prints them.
SEXP vlab(const double *x, int n) {
  SEXP s = PROTECT(Rf_allocVector(STRSXP, n));
  for (int i = 0; i < n; ++i) {
    char b[48];
    clab::fmt(x[i], b, sizeof b);
    SET_STRING_ELT(s, i, Rf_mkChar(b));
  }
  UNPROTECT(1);
  return s;
}

// The numbers x (one, or n) as n of them. Not protected.
SEXP rec(SEXP x, const char *nm, int n) {
  const R_xlen_t k = XLENGTH(x);
  if (k != 1 && k != n) fail("`%s` must be one number or one per scenario, %d.", nm, n);
  SEXP o = Rf_allocVector(REALSXP, n);
  for (int i = 0; i < n; ++i) REAL(o)[i] = REAL(x)[k == 1 ? 0 : i];
  return o;
}

}  // namespace

extern "C" {

SEXP C_clab_r_detrsp(SEXP x, SEXP ben, SEXP sanc, SEXP p) {
  clab::who = "detrsp";
  SEXP xv = PROTECT(bj::nums(x, "x"));
  const int n = bj::len(xv, "x");
  SEXP bv = PROTECT(bj::nums(ben, "ben", n));
  SEXP sv = PROTECT(bj::nums(sanc, "sanc", n));
  SEXP pv = PROTECT(bj::nums(p, "p", -1, 0.0));
  const int np = bj::len(pv, "p");
  std::vector<int> o(static_cast<std::size_t>(n));
  std::iota(o.begin(), o.end(), 0);
  std::stable_sort(o.begin(), o.end(), [&](int a, int b) { return REAL(xv)[a] < REAL(xv)[b]; });
  std::vector<double> xs(static_cast<std::size_t>(n)), bs(xs), ss(xs);
  for (int i = 0; i < n; ++i) {
    const int k = o[static_cast<std::size_t>(i)];
    xs[static_cast<std::size_t>(i)] = REAL(xv)[k];
    bs[static_cast<std::size_t>(i)] = REAL(bv)[k];
    ss[static_cast<std::size_t>(i)] = REAL(sv)[k];
  }
  bj::cols f(3, np);
  f.put("p", pv);
  double *xo = f.dbl("xopt"), *v = f.dbl("value");
  bj::kc(clab_detrsp(xs.data(), bs.data(), ss.data(), n, REAL(pv), np, xo, v),
         "`x` must not repeat a level and `sanc` must increase with `x`.");
  SEXP out = f.done();
  UNPROTECT(5);
  return out;
}

SEXP C_clab_r_dtrid(SEXP p, SEXP s, SEXP c, SEXP y, SEXP lg, SEXP rho) {
  clab::who = "dtrid";
  SEXP pv = PROTECT(bj::nums(p, "p"));
  const int n = bj::len(pv, "p");
  SEXP sv = PROTECT(bj::nums(s, "s", n));
  SEXP cv = PROTECT(bj::nums(c, "c", n));
  SEXP yv = Rf_isNull(y) ? R_NilValue : bj::nums(y, "y", n);
  PROTECT(yv);
  clab::lgl(lg, "log", 1);
  const int l = LOGICAL(lg)[0];
  double rh = 0;
  if (!Rf_isNull(rho)) rh = bj::inr(rho, "rho", clab::ninf, R_PosInf);
  if (l) {
    for (int i = 0; i < n; ++i) {
      if (!(REAL(pv)[i] > 0 && REAL(sv)[i] > 0)) {
        fail("`p` and `s` must be positive with `log = TRUE`.");
      }
    }
  }
  int rank = 0, id[3], cst[3];
  double beta[4];
  std::vector<double> w(8 * static_cast<std::size_t>(n));
  if (clab_dtrid(REAL(pv), REAL(sv), REAL(cv), n, l, Rf_isNull(yv) ? nullptr : REAL(yv), &rank, id,
                 cst, beta, w.data()) != CLAB_OK) {
    fail("the design has rank %d of 4: a dimension held fixed is not identified, so `y` "
         "cannot be fitted.",
         rank);
  }
  static const char *const dn[] = {"certainty", "severity", "celerity"};
  SEXP iv = PROTECT(Rf_allocVector(LGLSXP, 3)), kv = PROTECT(Rf_allocVector(LGLSXP, 3));
  for (int j = 0; j < 3; ++j) {
    LOGICAL(iv)[j] = id[j];
    LOGICAL(kv)[j] = cst[j];
  }
  SEXP nm3 = PROTECT(bfx::strs(dn, 3));
  Rf_setAttrib(iv, R_NamesSymbol, nm3);
  Rf_setAttrib(kv, R_NamesSymbol, nm3);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 6));
  SET_VECTOR_ELT(out, 0, Rf_ScalarInteger(rank));
  SET_VECTOR_ELT(out, 1, iv);
  SET_VECTOR_ELT(out, 2, kv);
  if (!Rf_isNull(yv)) {
    static const char *const bn[] = {"(Intercept)", "p", "s", "c"};
    static const char *const ln[] = {"(Intercept)", "log_p", "log_s", "c"};
    SET_VECTOR_ELT(out, 3, clab::nvec(4, l ? ln : bn, beta));
    SET_VECTOR_ELT(out, 4, Rf_ScalarReal(l ? (beta[1] != 0 ? -beta[3] / beta[1] : NA_REAL)
                                           : NA_REAL));
  }
  if (!Rf_isNull(rho)) {
    SEXP ec = Rf_allocVector(REALSXP, n);
    SET_VECTOR_ELT(out, 5, ec);
    for (int i = 0; i < n; ++i) {
      REAL(ec)[i] = REAL(pv)[i] * REAL(sv)[i] * std::exp(-rh * REAL(cv)[i]);
    }
  }
  static const char *const on[] = {"rank", "ident", "const", "coef", "disc", "cost"};
  UNPROTECT(8);
  return bj::named(out, on);
}

SEXP C_clab_r_inspg(SEXP a, SEXP x, SEXP steps, SEXP dt) {
  clab::who = "inspg";
  const double av = bj::inr(a, "a", clab::ninf, R_PosInf);
  if (av == 0) fail("`a` must not be 0.");
  SEXP xv = PROTECT(three(x, "x"));
  const double *x0 = REAL(xv);
  if (std::fabs(x0[0] + x0[1] + x0[2] - 1) > 1e-10) fail("`x` must be three shares summing to 1.");
  const int ns = clab::whole(steps, "steps", 0);
  if (ns > 10000000) fail("`steps` must be at most 1e7.");
  const double h = clab::one(dt, "dt", 0.0, R_PosInf);
  const int m = ns + 1;
  SEXP tr = PROTECT(Rf_allocMatrix(REALSXP, m, 3));
  SEXP v = PROTECT(Rf_allocVector(REALSXP, m));
  double f0[3], dv = 0;
  clab_inspg(av, x0, ns, h, REAL(tr), REAL(v), f0, &dv);
  double drift = 0;
  for (int k = 0; k < m; ++k) drift = std::max(drift, std::fabs(REAL(v)[k] - REAL(v)[0]));
  static const char *const xn[] = {"x1", "x2", "x3"};
  bfx::dnm(tr, R_NilValue, PROTECT(bfx::strs(xn, 3)));
  const double third[3] = {1.0 / 3, 1.0 / 3, 1.0 / 3};
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 6));
  SET_VECTOR_ELT(out, 0, tr);
  SET_VECTOR_ELT(out, 1, v);
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(drift));
  SET_VECTOR_ELT(out, 3, Rf_ScalarReal(dv));
  SET_VECTOR_ELT(out, 4, clab::nvec(3, xn, f0));
  SET_VECTOR_ELT(out, 5, clab::nvec(3, xn, third));
  static const char *const on[] = {"traj", "inv", "drift", "dvdt", "field", "fixed"};
  UNPROTECT(5);
  return bj::named(out, on);
}

SEXP C_clab_r_insphs(SEXP temp, SEXP cost, SEXP steps, SEXP rate, SEXP x0) {
  clab::who = "insphs";
  SEXP tv = PROTECT(bj::nums(temp, "temp"));
  SEXP gv = PROTECT(bj::nums(cost, "cost"));
  const int nt = bj::len(tv, "temp"), ng = bj::len(gv, "cost");
  if (static_cast<double>(nt) * ng > 1e7) fail("`temp` and `cost` make more than 1e7 cells.");
  const int ns = clab::whole(steps, "steps", 0);
  const double r = clab::one(rate, "rate", 0.0, R_PosInf);
  SEXP xv = PROTECT(three(x0, "x0"));
  if (!(REAL(xv)[0] + REAL(xv)[1] + REAL(xv)[2] > 0)) fail("`x0` must have a positive sum.");
  const R_xlen_t cells = static_cast<R_xlen_t>(nt) * ng;
  std::vector<double> o(3 * static_cast<std::size_t>(cells));
  clab_insphs(REAL(tv), nt, REAL(gv), ng, ns, r, REAL(xv), o.data());
  SEXP rn = PROTECT(vlab(REAL(tv), nt)), cn = PROTECT(vlab(REAL(gv), ng));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
  // the kernel's order is C, P, O; the result's crime (P), coop (C), insp (O)
  static const int at[3] = {1, 0, 2};
  for (int k = 0; k < 3; ++k) {
    SEXP m = bfx::dmat(nt, ng, o.data() + at[k] * cells);
    SET_VECTOR_ELT(out, k, m);
    bfx::dnm(m, rn, cn);
  }
  SET_VECTOR_ELT(out, 3, tv);
  SET_VECTOR_ELT(out, 4, gv);
  static const char *const on[] = {"crime", "coop", "insp", "temp", "cost"};
  UNPROTECT(6);
  return bj::named(out, on);
}

SEXP C_clab_r_colarr(SEXP a, SEXP e, SEXP bg) {
  clab::who = "colarr";
  SEXP a0 = PROTECT(bj::nums(a, "a")), e0 = PROTECT(bj::nums(e, "e"));
  SEXP b0 = PROTECT(bj::nums(bg, "bg"));
  const int n = std::max({bj::len(a0, "a"), bj::len(e0, "e"), bj::len(b0, "bg")});
  SEXP av = PROTECT(rec(a0, "a", n)), ev = PROTECT(rec(e0, "e", n)), bv = PROTECT(rec(b0, "bg", n));
  bj::cols f(9, n);
  f.put("a", av);
  f.put("e", ev);
  f.put("bg", bv);
  double *orr = f.dbl("or"), *lor = f.dbl("lor");
  std::vector<double> c(4 * static_cast<std::size_t>(n));
  bj::kc(clab_colarr(REAL(av), REAL(ev), REAL(bv), n, c.data(), orr),
         "`a` and `e` must be strictly between 0 and 1, `bg` from 0 to 1.");
  static const char *const cn[] = {"c11", "c10", "c01", "c00"};
  for (int k = 0; k < 4; ++k) {
    double *d = f.dbl(cn[k]);
    std::copy(c.begin() + k * n, c.begin() + (k + 1) * n, d);
  }
  for (int i = 0; i < n; ++i) lor[i] = std::log(orr[i]);
  SEXP out = f.done();
  UNPROTECT(7);
  return out;
}

}  // extern "C"
