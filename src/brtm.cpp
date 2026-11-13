// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of geographic profiling and risk terrain (crchyp,
// hitscr, rtmfit): each checks its arguments, calls the kernels and
// builds the R result. Nothing calls back into the R evaluator.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "bfx.h"
#include "rchk.h"

extern "C" double Rf_pnorm5(double, double, double, int, int);

using clab::fail;

namespace {

// The layers of rtmfit given as columns: a matrix or vector of n rows, or
// a list of n numbers each. Columns into L (n by k), names into nm.
int cols(SEXP lay, int n, std::vector<double> &L, std::vector<std::string> &nm) {
  int k = 0;
  if (TYPEOF(lay) == VECSXP) {
    k = static_cast<int>(XLENGTH(lay));
    if (k < 1) fail("`lay` must hold at least one layer.");
    SEXP ln = Rf_getAttrib(lay, R_NamesSymbol);
    for (int j = 0; j < k; ++j) {
      SEXP v = PROTECT(bj::nums(VECTOR_ELT(lay, j), "lay", n));
      L.insert(L.end(), REAL(v), REAL(v) + n);
      UNPROTECT(1);
      nm.emplace_back(Rf_isNull(ln) || !*bj::str(ln, j) ? "" : bj::str(ln, j));
    }
  } else {
    int r = 0;
    SEXP v = PROTECT(bfx::mat(lay, "lay", -1, &r));
    if (r != n) fail("`lay` must have one row per cell, %d.", n);
    k = static_cast<int>(XLENGTH(v) / n);
    L.assign(REAL(v), REAL(v) + XLENGTH(v));
    SEXP dn = Rf_getAttrib(lay, R_DimNamesSymbol);
    SEXP cn = Rf_isNull(dn) ? R_NilValue : VECTOR_ELT(dn, 1);
    for (int j = 0; j < k; ++j) nm.emplace_back(Rf_isNull(cn) ? "" : bj::str(cn, j));
    UNPROTECT(1);
  }
  return k;
}

// The layers of rtmfit built from features: cell (n by 2), lay a list of
// feature matrices (or one), each two columns.
int feats(SEXP lay, SEXP cell, SEXP rad, SEXP kind, int n, std::vector<double> &L,
          std::vector<std::string> &nm) {
  int nr = 0;
  SEXP c = PROTECT(bfx::mat(cell, "cell", 2, &nr));
  if (nr != n) fail("`cell` must have one row per count, %d.", n);
  const double r = clab::one(rad, "rad", 0.0, R_PosInf);
  static const char *const ko[] = {"prox", "dens", "kde"};
  const int kd = csp::pick(kind, "kind", ko, 3);
  if (kd == 2 && n < 2) fail("`kind = \"kde\"` needs at least two cells.");
  const bool one = TYPEOF(lay) != VECSXP;
  const int k = one ? 1 : static_cast<int>(XLENGTH(lay));
  if (k < 1) fail("`lay` must hold at least one layer.");
  SEXP ln = one ? R_NilValue : Rf_getAttrib(lay, R_NamesSymbol);
  L.resize(static_cast<std::size_t>(n) * k);
  for (int j = 0; j < k; ++j) {
    int nf = 0;
    SEXP f = PROTECT(bfx::mat(one ? lay : VECTOR_ELT(lay, j), "lay", 2, &nf));
    clab_rtmlay(REAL(c), REAL(c) + n, n, REAL(f), REAL(f) + nf, nf, r, kd,
                &L[static_cast<std::size_t>(n) * j]);
    UNPROTECT(1);
    nm.emplace_back(Rf_isNull(ln) || !*bj::str(ln, j) ? "" : bj::str(ln, j));
  }
  UNPROTECT(1);
  return k;
}

}  // namespace

extern "C" {

SEXP C_clab_r_crchyp(SEXP x, SEXP y, SEXP home) {
  clab::who = "crchyp";
  SEXP xv = PROTECT(bj::nums(x, "x"));
  const int n = bj::len(xv, "x");
  SEXP yv = PROTECT(bj::nums(y, "y", n));
  if (n < 2) fail("`x` and `y` must give at least two sites.");
  int nh = 0;
  SEXP hv = R_NilValue;
  if (!Rf_isNull(home)) {
    SEXP d = Rf_getAttrib(home, R_DimSymbol);
    if (Rf_isNull(d) && clab::isnum(home) && XLENGTH(home) == 2) {
      hv = PROTECT(bj::nums(home, "home", 2));
      nh = 1;
    } else {
      hv = PROTECT(bfx::mat(home, "home", 2, &nh));
    }
  } else {
    PROTECT(hv);
  }
  double st[4];
  int pair[2];
  std::vector<int> mar(static_cast<std::size_t>(nh) + 1);
  const double *hx = nh ? REAL(hv) : nullptr;
  clab_crchyp(REAL(xv), REAL(yv), n, hx, nh ? hx + nh : nullptr, nh, st, pair, mar.data());
  static const char *const cn[] = {"x", "y"};
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
  SET_VECTOR_ELT(out, 0, clab::nvec(2, cn, st));
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(st[2]));
  SEXP pr = Rf_allocVector(INTSXP, 2);
  SET_VECTOR_ELT(out, 2, pr);
  INTEGER(pr)[0] = pair[0] + 1;
  INTEGER(pr)[1] = pair[1] + 1;
  SET_VECTOR_ELT(out, 3, Rf_ScalarReal(st[3]));
  if (nh) {
    SEXP m = Rf_allocVector(LGLSXP, nh);
    SET_VECTOR_ELT(out, 4, m);
    for (int k = 0; k < nh; ++k) LOGICAL(m)[k] = mar[static_cast<std::size_t>(k)];
  }
  static const char *const nm[] = {"centre", "r", "pair", "share", "mar"};
  UNPROTECT(4);
  return bj::named(out, nm);
}

SEXP C_clab_r_hitscr(SEXP s, SEXP home) {
  clab::who = "hitscr";
  SEXP sv = PROTECT(bj::nums(s, "s"));
  const int n = bj::len(sv, "s");
  SEXP hv = PROTECT(bj::nums(home, "home", -1, 1.0));
  const int nh = bj::len(hv, "home");
  std::vector<int> h(static_cast<std::size_t>(nh));
  for (int k = 0; k < nh; ++k) {
    const double v = REAL(hv)[k];
    if (v != std::floor(v) || v > n) fail("`home` must be cell numbers from 1 to %d.", n);
    h[static_cast<std::size_t>(k)] = static_cast<int>(v) - 1;
  }
  bj::cols f(6, nh);
  int *hm = f.itg("home");
  for (int k = 0; k < nh; ++k) hm[k] = h[static_cast<std::size_t>(k)] + 1;
  int *ab = f.itg("above"), *ti = f.itg("ties"), *rk = f.itg("rank");
  double *hit = f.dbl("hit"), *hitt = f.dbl("hitt");
  clab_hitscr(REAL(sv), n, h.data(), nh, ab, ti, hit, hitt);
  for (int k = 0; k < nh; ++k) rk[k] = ab[k] + 1;
  SEXP out = f.done();
  UNPROTECT(3);
  return out;
}

SEXP C_clab_r_rtmfit(SEXP y, SEXP lay, SEXP cell, SEXP rad, SEXP kind, SEXP level) {
  clab::who = "rtmfit";
  SEXP yv = PROTECT(bj::nums(y, "y", -1, 0.0));
  const int n = bj::len(yv, "y");
  bool any = false;
  for (int i = 0; i < n; ++i) any = any || REAL(yv)[i] > 0;
  if (!any) fail("`y` must hold at least one count above 0.");
  const double lv = clab::one(level, "level", 0.0, 1.0);
  std::vector<double> L;
  std::vector<std::string> nm;
  const int k = Rf_isNull(cell) ? cols(lay, n, L, nm) : feats(lay, cell, rad, kind, n, L, nm);
  for (int j = 0; j < k; ++j) {
    if (nm[static_cast<std::size_t>(j)].empty()) {
      char b[32];
      std::snprintf(b, sizeof b, "layer%d", j + 1);
      nm[static_cast<std::size_t>(j)] = b;
    }
  }
  const int p = k + 1;
  std::vector<double> x(static_cast<std::size_t>(n) * p, 1.0);
  std::copy(L.begin(), L.end(), x.begin() + n);
  std::vector<double> beta(static_cast<std::size_t>(p)), vb(static_cast<std::size_t>(p) * p),
      vc(1), st(8),
      w(static_cast<std::size_t>(n) * p + 6 * static_cast<std::size_t>(n) + p +
        static_cast<std::size_t>(p) * p);
  std::vector<int> al(static_cast<std::size_t>(p)), iw(static_cast<std::size_t>(p));
  bj::kc(clab_cntreg(REAL(yv), x.data(), n, p, nullptr, nullptr, 0, 0, 100, 1e-14, beta.data(),
                     vb.data(), vc.data(), st.data(), al.data(), w.data(), iw.data()),
         "the Poisson fit failed: the deviance left the finite range.");
  const double q = Rf_qnorm5(0.5 + lv / 2, 0.0, 1.0, 1, 0);
  bj::cols f(8, p);
  SEXP tn = PROTECT(Rf_allocVector(STRSXP, p));
  SET_STRING_ELT(tn, 0, Rf_mkChar("(Intercept)"));
  for (int j = 0; j < k; ++j) {
    SET_STRING_ELT(tn, j + 1, Rf_mkCharCE(nm[static_cast<std::size_t>(j)].c_str(), CE_UTF8));
  }
  f.put("term", tn);
  double *est = f.dbl("est"), *se = f.dbl("se"), *z = f.dbl("z"), *pv = f.dbl("p"),
         *rrv = f.dbl("rrv"), *lo = f.dbl("lo"), *hi = f.dbl("hi");
  for (int j = 0; j < p; ++j) {
    const double b = beta[static_cast<std::size_t>(j)];
    const double s = std::sqrt(vb[static_cast<std::size_t>(j) * p + j]);
    est[j] = bj::na(b);
    se[j] = bj::na(s);
    z[j] = bj::na(b / s);
    pv[j] = bj::na(2 * Rf_pnorm5(-std::fabs(b / s), 0.0, 1.0, 1, 0));
    rrv[j] = bj::na(std::exp(b));
    lo[j] = bj::na(std::exp(b - q * s));
    hi[j] = bj::na(std::exp(b + q * s));
  }
  bj::cols g(4, n);
  double *rr = g.dbl("rr"), *sc = g.dbl("score"), *fit = g.dbl("fitted");
  g.put("y", yv);
  double mn = R_PosInf;
  for (int i = 0; i < n; ++i) {
    double e = 0;
    for (int j = 1; j < p; ++j) {
      if (!al[static_cast<std::size_t>(j)]) {
        e += L[static_cast<std::size_t>(j - 1) * n + i] * beta[static_cast<std::size_t>(j)];
      }
    }
    rr[i] = std::exp(e);
    fit[i] = std::exp(e + beta[0]);
    mn = std::min(mn, rr[i]);
  }
  for (int i = 0; i < n; ++i) sc[i] = rr[i] / mn;
  SEXP lm = PROTECT(bfx::dmat(n, k, L.data()));
  SEXP ln = PROTECT(Rf_allocVector(STRSXP, k));
  for (int j = 0; j < k; ++j) SET_STRING_ELT(ln, j, STRING_ELT(tn, j + 1));
  bfx::dnm(lm, R_NilValue, ln);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 7));
  SET_VECTOR_ELT(out, 0, f.done());
  SET_VECTOR_ELT(out, 1, g.done());
  SET_VECTOR_ELT(out, 2, lm);
  SET_VECTOR_ELT(out, 3, Rf_ScalarReal(st[7]));
  SET_VECTOR_ELT(out, 4, Rf_ScalarReal(st[1]));
  SET_VECTOR_ELT(out, 5, Rf_ScalarLogical(st[5] == 1));
  SET_VECTOR_ELT(out, 6, Rf_ScalarReal(lv));
  static const char *const on[] = {"coef", "cells", "lay", "dev", "aic", "conv", "level"};
  UNPROTECT(7);
  return bj::named(out, on);
}

}  // extern "C"
