// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the area functions (bndxw, synexp): each checks
// its arguments, calls the kernels and builds the R result. Nothing calls
// back into the R evaluator.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <unordered_map>
#include <vector>

#include "bfx.h"
#include "rchk.h"

using clab::fail;

extern "C" {

SEXP C_clab_r_bndxw(SEXP from, SEXP to, SEXP x, SEXP unit, SEXP w) {
  clab::who = "bndxw";
  const R_xlen_t m0 = XLENGTH(from);
  if (m0 < 1) fail("`from` must hold at least one link.");
  const int m = clab::ilen(static_cast<double>(m0), "from");
  csp::lab f = bj::labs(from, "from", m);
  csp::lab g = bj::labs(to, "to", m);
  std::vector<double> wv(static_cast<std::size_t>(m), 1.0);
  if (!Rf_isNull(w)) {
    SEXP d = PROTECT(bj::nums(w, "w", m, 0.0));
    std::copy(REAL(d), REAL(d) + m, wv.begin());
    UNPROTECT(1);
  }
  SEXP xv = PROTECT(bj::nums(x, "x"));
  const int n = bj::len(xv, "x");
  if (XLENGTH(unit) != n) fail("`unit` must have one label per value, %d.", n);
  SEXP ul = PROTECT(clab::labs(unit));
  if (Rf_isNull(ul) || clab::anyna(ul)) fail("`unit` must be labels with no missing values.");
  std::vector<int> u(static_cast<std::size_t>(n));
  double lost = 0;
  int nlost = 0;
  SEXP miss = PROTECT(Rf_allocVector(STRSXP, n));
  {
    std::unordered_map<std::string, int> at;
    for (int s = 0; s < f.g; ++s) at.emplace(bj::str(f.lev, s), s);
    std::unordered_map<std::string, int> seen;
    for (int i = 0; i < n; ++i) {
      auto it = at.find(bj::str(ul, i));
      u[static_cast<std::size_t>(i)] = it == at.end() ? -1 : it->second;
      if (it == at.end()) {
        lost += REAL(xv)[i];
        if (seen.emplace(bj::str(ul, i), 1).second) SET_STRING_ELT(miss, nlost++, STRING_ELT(ul, i));
      }
    }
  }
  std::vector<double> ws(static_cast<std::size_t>(f.g)), bw(ws);
  std::vector<int> best(static_cast<std::size_t>(f.g));
  bj::cols tb(2, g.g);
  tb.put("to", g.lev);
  double *out = tb.dbl("value");
  clab_bndxw(REAL(xv), u.data(), n, INTEGER(f.code), INTEGER(g.code), wv.data(), m, f.g, g.g, out,
             ws.data(), best.data(), bw.data());
  bj::cols mp(4, f.g);
  mp.put("from", f.lev);
  SEXP bt = PROTECT(Rf_allocVector(STRSXP, f.g));
  mp.put("to", bt);
  double *sh = mp.dbl("share"), *sw = mp.dbl("wsum");
  for (int s = 0; s < f.g; ++s) {
    const std::size_t k = static_cast<std::size_t>(s);
    SET_STRING_ELT(bt, s, best[k] < 0 ? NA_STRING : STRING_ELT(g.lev, best[k]));
    sh[s] = bj::na(bw[k]);
    sw[s] = ws[k];
  }
  SEXP mv = PROTECT(Rf_lengthgets(miss, nlost));
  SEXP out2 = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out2, 0, tb.done());
  SET_VECTOR_ELT(out2, 1, mp.done());
  SET_VECTOR_ELT(out2, 2, Rf_ScalarReal(lost));
  SET_VECTOR_ELT(out2, 3, mv);
  static const char *const nm[] = {"tab", "map", "lost", "unmatched"};
  UNPROTECT(12);
  return bj::named(out2, nm);
}

SEXP C_clab_r_synexp(SEXP y, SEXP x, SEXP ax, SEXP pop) {
  clab::who = "synexp";
  int n = 0, na = 0;
  SEXP xv = PROTECT(bfx::mat(x, "x", -1, &n));
  const int p = static_cast<int>(XLENGTH(xv) / n);
  SEXP yv = PROTECT(bj::nums(y, "y", n));
  int n0 = 0, n1 = 0;
  for (int i = 0; i < n; ++i) {
    const double v = REAL(yv)[i];
    if (v != 0 && v != 1) fail("`y` must be 0 or 1.");
    n0 += v == 0;
    n1 += v == 1;
  }
  if (!n0 || !n1) fail("`y` must hold both 0 and 1.");
  if (n <= p + 1) fail("`x` needs more rows than columns plus one.");
  SEXP av = PROTECT(bfx::mat(ax, "ax", p, &na));
  SEXP pv = PROTECT(bj::nums(pop, "pop", na, 0.0));
  const std::size_t q = static_cast<std::size_t>(p) + 1;
  std::vector<double> cf(q), se(q), wk(2 * q * q + 4 * q);
  std::vector<int> iw(static_cast<std::size_t>(n));
  double dev = 0;
  bj::kc(clab_lgtreg(REAL(xv), n, p, REAL(yv), 1e-10, 50, cf.data(), se.data(), &dev, wk.data(),
                     iw.data()),
         "no fit: collinear columns of `x`, or separated outcomes.");
  bj::cols c(3, p + 1);
  SEXP tn = PROTECT(Rf_allocVector(STRSXP, p + 1));
  c.put("term", tn);
  SET_STRING_ELT(tn, 0, Rf_mkChar("(Intercept)"));
  SEXP dn = Rf_getAttrib(x, R_DimNamesSymbol);
  SEXP cn = Rf_isNull(dn) ? R_NilValue : VECTOR_ELT(dn, 1);
  for (int j = 0; j < p; ++j) {
    char b[24];
    std::snprintf(b, sizeof b, "x%d", j + 1);
    SET_STRING_ELT(tn, j + 1, Rf_isNull(cn) ? Rf_mkChar(b) : STRING_ELT(cn, j));
  }
  double *est = c.dbl("est"), *sev = c.dbl("se");
  std::copy(cf.begin(), cf.end(), est);
  std::copy(se.begin(), se.end(), sev);
  bj::cols a(3, na);
  double *rate = a.dbl("rate");
  a.put("pop", pv);
  double *ex = a.dbl("expo");
  clab_synexp(cf.data(), REAL(av), na, p, REAL(pv), rate, ex);
  double tot = 0;
  for (int i = 0; i < na; ++i) tot += ex[i];
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out, 0, c.done());
  SET_VECTOR_ELT(out, 1, a.done());
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(tot));
  SET_VECTOR_ELT(out, 3, Rf_ScalarReal(dev));
  static const char *const nm[] = {"coef", "areas", "total", "dev"};
  UNPROTECT(8);
  return bj::named(out, nm);
}

}  // extern "C"
