// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of taphonomy (clrilr, preslr, tpdmc, tpbhm): each
// checks its arguments, calls the kernels and builds the R result.
// Nothing calls back into the R evaluator.
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

#include "bfx.h"
#include "rchk.h"

extern "C" double Rf_pnorm5(double, double, double, int, int);

using clab::fail;

namespace {

// One finite number or n of them (pos: positive), as n numbers. Not
// protected.
SEXP oneorn(SEXP v, const char *nm, int n, bool pos) {
  SEXP d = PROTECT(pos ? bj::nums(v, nm, -1, 0.0, true) : clab::num(v, nm));
  const R_xlen_t k = XLENGTH(d);
  if (k != 1 && k != n) fail("`%s` must be one number or one per value, %d.", nm, n);
  SEXP o = Rf_allocVector(REALSXP, n);
  for (int i = 0; i < n; ++i) REAL(o)[i] = REAL(d)[k == 1 ? 0 : i];
  UNPROTECT(1);
  return o;
}

}  // namespace

extern "C" {

SEXP C_clab_r_clrilr(SEXP x, SEXP type, SEXP ps) {
  clab::who = "clrilr";
  int n = 0;
  SEXP xv = PROTECT(bfx::mat(x, "x", -1, &n, 0.0));
  const int d = static_cast<int>(XLENGTH(xv) / n);
  if (d < 2) fail("`x` must have at least two parts (columns).");
  static const char *const to[] = {"clr", "ilr"};
  const int ty = csp::pick(type, "type", to, 2);
  const double pc = bj::inr(ps, "ps", 0.0, R_PosInf);
  const int dc = ty == 0 ? d : d - 1;
  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, n, dc));
  bj::kc(clab_clrilr(REAL(xv), n, d, pc, ty, REAL(out)),
         "a part is 0 after the pseudocount: give `ps` above 0.");
  SEXP dn = Rf_getAttrib(x, R_DimNamesSymbol);
  SEXP rn = Rf_isNull(dn) ? R_NilValue : VECTOR_ELT(dn, 0);
  SEXP cn = Rf_isNull(dn) ? R_NilValue : VECTOR_ELT(dn, 1);
  if (ty == 1) {
    cn = PROTECT(Rf_allocVector(STRSXP, dc));
    for (int j = 0; j < dc; ++j) {
      char b[24];
      std::snprintf(b, sizeof b, "ilr%d", j + 1);
      SET_STRING_ELT(cn, j, Rf_mkChar(b));
    }
  } else {
    PROTECT(cn);
  }
  bfx::dnm(out, rn, cn);
  UNPROTECT(3);
  return out;
}

SEXP C_clab_r_preslr(SEXP x, SEXP m1, SEXP s1, SEXP m2, SEXP s2) {
  clab::who = "preslr";
  SEXP xv = PROTECT(bj::nums(x, "x"));
  const int n = bj::len(xv, "x");
  SEXP a = PROTECT(oneorn(m1, "m1", n, false)), b = PROTECT(oneorn(s1, "s1", n, true));
  SEXP c = PROTECT(oneorn(m2, "m2", n, false)), d = PROTECT(oneorn(s2, "s2", n, true));
  double o[5];
  int band = 0, side = 0;
  clab_preslr(REAL(xv), n, REAL(a), REAL(b), n, REAL(c), REAL(d), n, o, &band, &side);
  static const char *const bw[] = {"no support either way", "weak support", "moderate support",
                                   "moderately strong support", "strong support",
                                   "very strong support", "extremely strong support"};
  char vb[64];
  std::snprintf(vb, sizeof vb, "%s%s", bw[band], side == 1 ? " for H1" : side == 2 ? " for H2" : "");
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 6));
  for (int j = 0; j < 5; ++j) SET_VECTOR_ELT(out, j, Rf_ScalarReal(o[j]));
  SET_VECTOR_ELT(out, 5, Rf_mkString(vb));
  static const char *const nm[] = {"ll1", "ll2", "llr", "lr", "l10", "verbal"};
  UNPROTECT(6);
  return bj::named(out, nm);
}

SEXP C_clab_r_tpdmc(SEXP pres, SEXP decay, SEXP mum, SEXP states, SEXP start, SEXP steps,
                    SEXP seed) {
  clab::who = "tpdmc";
  const double pr = clab::p01(pres, "pres"), mr = clab::p01(mum, "mum");
  const double dr = bj::inr(decay, "decay", 0.0, 1.0);
  if (!(dr > 0)) fail("`decay` must be above 0.");
  if (TYPEOF(states) != STRSXP || XLENGTH(states) < 1 || clab::anyna(states)) {
    fail("`states` must name the transient stages: labels, at least one.");
  }
  const int k = static_cast<int>(XLENGTH(states)), m = k + 2;
  for (int i = 0; i < k; ++i) {
    for (int j = 0; j < i; ++j) {
      if (clab::seql(STRING_ELT(states, i), STRING_ELT(states, j))) {
        fail("`states` must name each stage once.");
      }
    }
  }
  int s0 = -1;
  if (TYPEOF(start) == STRSXP && XLENGTH(start) == 1) {
    for (int i = 0; i < k; ++i) {
      if (clab::seql(STRING_ELT(start, 0), STRING_ELT(states, i))) s0 = i;
    }
  } else if (clab::isnum(start) && XLENGTH(start) == 1) {
    const double v = clab::val1(start);
    if (v == std::floor(v) && v >= 1 && v <= k) s0 = static_cast<int>(v) - 1;
  }
  if (s0 < 0) fail("`start` must be a transient stage, by name or number from 1 to %d.", k);
  const int ns = clab::whole(steps, "steps", 0);
  SEXP pm = PROTECT(Rf_allocMatrix(REALSXP, m, m));
  SEXP nm = PROTECT(Rf_allocMatrix(REALSXP, k, k));
  SEXP bm = PROTECT(Rf_allocMatrix(REALSXP, k, 2));
  SEXP tv = PROTECT(Rf_allocVector(REALSXP, k));
  bj::kc(clab_tpdmc(pr, dr, mr, k, REAL(pm), REAL(nm), REAL(bm), REAL(tv)),
         "nothing leaves a stage: with `pres` 1, `mum` must be above 0.");
  SEXP all = PROTECT(Rf_allocVector(STRSXP, m));
  for (int i = 0; i < k; ++i) SET_STRING_ELT(all, i, STRING_ELT(states, i));
  SET_STRING_ELT(all, k, Rf_mkChar("skeletal"));
  SET_STRING_ELT(all, k + 1, Rf_mkChar("mummified"));
  SEXP ab = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(ab, 0, Rf_mkChar("skeletal"));
  SET_STRING_ELT(ab, 1, Rf_mkChar("mummified"));
  bfx::dnm(pm, all, all);
  bfx::dnm(nm, states, states);
  bfx::dnm(bm, states, ab);
  Rf_setAttrib(tv, R_NamesSymbol, states);
  SEXP path = R_NilValue;
  if (ns > 0) {
    const unsigned int sd = csp::seed(seed);
    std::vector<int> p(static_cast<std::size_t>(ns) + 1);
    int len = 0;
    clab_tpdsim(REAL(pm), m, k, s0, ns, sd, p.data(), &len);
    path = Rf_allocVector(STRSXP, len);
    for (int i = 0; i < len; ++i) SET_STRING_ELT(path, i, STRING_ELT(all, p[static_cast<std::size_t>(i)]));
  }
  PROTECT(path);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 6));
  SET_VECTOR_ELT(out, 0, pm);
  SET_VECTOR_ELT(out, 1, nm);
  SET_VECTOR_ELT(out, 2, bm);
  SET_VECTOR_ELT(out, 3, tv);
  SET_VECTOR_ELT(out, 4, Rf_ScalarInteger(s0 + 1));
  SET_VECTOR_ELT(out, 5, path);
  static const char *const on[] = {"P", "N", "absorb", "time", "start", "path"};
  UNPROTECT(8);
  return bj::named(out, on);
}

SEXP C_clab_r_tpbhm(SEXP y, SEXP x, SEXP grp, SEXP pm, SEXP psd, SEXP level) {
  clab::who = "tpbhm";
  SEXP yv = PROTECT(bj::nums(y, "y"));
  const int n = bj::len(yv, "y");
  int nr = 0;
  SEXP xv = PROTECT(bfx::mat(x, "x", -1, &nr));
  if (nr != n) fail("`x` must have one row per case, %d.", n);
  const int q = static_cast<int>(XLENGTH(xv) / n), p = q + 1;
  if (n < 2) fail("`y` must hold at least two cases.");
  SEXP mv = PROTECT(clab::num(pm, "pm"));
  SEXP sv = PROTECT(bj::nums(psd, "psd", -1, 0.0, true));
  if ((XLENGTH(mv) != 1 && XLENGTH(mv) != p) || (XLENGTH(sv) != 1 && XLENGTH(sv) != p)) {
    fail("`pm` and `psd` must be one number or one per term, %d.", p);
  }
  const double lv = clab::one(level, "level", 0.0, 1.0);
  std::vector<double> m0(static_cast<std::size_t>(p)), s0(static_cast<std::size_t>(p));
  for (int j = 0; j < p; ++j) {
    m0[static_cast<std::size_t>(j)] = REAL(mv)[XLENGTH(mv) == 1 ? 0 : j];
    s0[static_cast<std::size_t>(j)] = REAL(sv)[XLENGTH(sv) == 1 ? 0 : j];
  }
  int ng = 0;
  const int *gc = nullptr;
  SEXP glev = R_NilValue;
  if (!Rf_isNull(grp)) {
    csp::lab g = bj::labs(grp, "grp", n);
    ng = g.g;
    gc = INTEGER(g.code);
    glev = g.lev;
    if (ng < 2) fail("`grp` must hold at least two groups.");
  } else {
    PROTECT(glev);
    PROTECT(glev);
  }
  std::vector<double> mu(static_cast<std::size_t>(p)), sd(static_cast<std::size_t>(p)), st(2),
      gm(static_cast<std::size_t>(ng) + 1), lam(static_cast<std::size_t>(ng) + 1),
      w(2 * static_cast<std::size_t>(p) * p + p);
  std::vector<int> nj(static_cast<std::size_t>(ng) + 1);
  bj::kc(clab_tpbhm(REAL(yv), REAL(xv), n, q, m0.data(), s0.data(), gc, ng, mu.data(), sd.data(),
                    st.data(), gm.data(), lam.data(), nj.data(), w.data()),
         "`x` has dependent columns (or fewer cases than terms), or `y` does not vary.");
  const double zq = Rf_qnorm5(0.5 + lv / 2, 0.0, 1.0, 1, 0);
  SEXP tn = PROTECT(Rf_allocVector(STRSXP, p));
  SET_STRING_ELT(tn, 0, Rf_mkChar("(Intercept)"));
  SEXP dn = Rf_getAttrib(x, R_DimNamesSymbol);
  SEXP cn = Rf_isNull(dn) ? R_NilValue : VECTOR_ELT(dn, 1);
  for (int j = 1; j < p; ++j) {
    char b[24];
    std::snprintf(b, sizeof b, "x%d", j);
    SET_STRING_ELT(tn, j, Rf_isNull(cn) ? Rf_mkChar(b) : STRING_ELT(cn, j - 1));
  }
  bj::cols f(7, p);
  f.put("term", tn);
  double *pmu = f.dbl("mean"), *psd2 = f.dbl("sd"), *lo = f.dbl("lo"), *hi = f.dbl("hi"),
         *pp = f.dbl("ppos"), *pr = f.dbl("prior");
  for (int j = 0; j < p; ++j) {
    const double a = mu[static_cast<std::size_t>(j)], s = sd[static_cast<std::size_t>(j)];
    pmu[j] = a;
    psd2[j] = s;
    lo[j] = a - zq * s;
    hi[j] = a + zq * s;
    pp[j] = Rf_pnorm5(a / s, 0.0, 1.0, 1, 0);
    pr[j] = s0[static_cast<std::size_t>(j)];
  }
  SEXP fit = PROTECT(Rf_allocVector(REALSXP, n));
  for (int i = 0; i < n; ++i) {
    double v = mu[0];
    for (int j = 1; j < p; ++j) v += REAL(xv)[static_cast<R_xlen_t>(j - 1) * n + i] * mu[static_cast<std::size_t>(j)];
    REAL(fit)[i] = v;
  }
  SEXP gr = R_NilValue;
  if (ng) {
    bj::cols h(5, ng);
    h.put("group", glev);
    int *hn = h.itg("n");
    double *raw = h.dbl("raw"), *sh = h.dbl("shrink"), *po = h.dbl("pooled");
    for (int j = 0; j < ng; ++j) {
      const std::size_t u = static_cast<std::size_t>(j);
      hn[j] = nj[u];
      raw[j] = gm[u];
      sh[j] = lam[u];
      po[j] = lam[u] * gm[u];
    }
    gr = h.done();
    UNPROTECT(1);
  }
  PROTECT(gr);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 6));
  SET_VECTOR_ELT(out, 0, f.done());
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(std::sqrt(st[0])));
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(bj::na(st[1])));
  SET_VECTOR_ELT(out, 3, gr);
  SET_VECTOR_ELT(out, 4, fit);
  SET_VECTOR_ELT(out, 5, Rf_ScalarReal(lv));
  static const char *const on[] = {"coef", "sigma", "tau2", "groups", "fitted", "level"};
  UNPROTECT(11);
  return bj::named(out, on);
}

}  // extern "C"
