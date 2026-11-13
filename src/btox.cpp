// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of forensic toxicology (txcal, lodimp, pmrrat,
// etohlr, mtxrel): each checks its arguments, calls the kernel and builds
// the R result. Nothing calls back into the R evaluator.
#include <cmath>
#include <vector>

#include "bfx.h"
#include "rchk.h"

using clab::fail;

namespace {

// A factor of the codes v (n; -1 for NA) over the k labels lev.
SEXP fac(const int *v, int n, const char *const *lev, int k) {
  SEXP f = PROTECT(Rf_allocVector(INTSXP, n));
  for (int i = 0; i < n; ++i) INTEGER(f)[i] = v[i] < 0 ? NA_INTEGER : v[i] + 1;
  Rf_setAttrib(f, R_LevelsSymbol, bfx::strs(lev, k));
  Rf_setAttrib(f, R_ClassSymbol, Rf_mkString("factor"));
  UNPROTECT(1);
  return f;
}

}  // namespace

extern "C" {

SEXP C_clab_r_txcal(SEXP conc, SEXP resp, SEXP w, SEXP unk, SEXP k) {
  clab::who = "txcal";
  SEXP xv = PROTECT(bj::nums(conc, "conc", -1, 0.0, true));
  const int n = bj::len(xv, "conc");
  SEXP yv = PROTECT(bj::nums(resp, "resp", n));
  if (n < 3) fail("`conc` and `resp` must give at least three calibrators.");
  SEXP kv = PROTECT(bj::nums(k, "k", 2, 0.0, true));
  std::vector<double> wt(static_cast<std::size_t>(n));
  const char *wl = "numeric";
  if (TYPEOF(w) == STRSXP) {
    static const char *const wo[] = {"1/x^2", "1/x", "none"};
    const int m = csp::pick(w, "w", wo, 3);
    wl = wo[m];
    for (int i = 0; i < n; ++i) {
      const double x = REAL(xv)[i];
      wt[static_cast<std::size_t>(i)] = m == 0 ? 1 / (x * x) : m == 1 ? 1 / x : 1.0;
    }
  } else {
    SEXP wv = PROTECT(bj::nums(w, "w", n, 0.0, true));
    std::copy(REAL(wv), REAL(wv) + n, wt.begin());
    UNPROTECT(1);
  }
  int nu = 0;
  SEXP uv = R_NilValue;
  if (!Rf_isNull(unk)) {
    uv = bfx::nna(unk, "unk", -1);
    nu = bj::len(uv, "unk");
  }
  PROTECT(uv);
  double st[6];
  std::vector<double> xh(static_cast<std::size_t>(nu) + 1);
  std::vector<int> fl(static_cast<std::size_t>(nu) + 1);
  bj::kc(clab_txcal(REAL(xv), REAL(yv), wt.data(), n, nu ? REAL(uv) : nullptr, nu, REAL(kv)[0],
                    REAL(kv)[1], st, xh.data(), fl.data()),
         "the calibrators give no line: `conc` does not vary or the slope is 0.");
  static const char *const fn[] = {"int", "slope", "r2", "sigma", "lod", "loq"};
  for (double &v : st) v = bj::na(v);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, clab::nvec(6, fn, st));
  if (nu) {
    bj::cols f(3, nu);
    f.put("resp", uv);
    double *c = f.dbl("conc");
    for (int j = 0; j < nu; ++j) c[j] = bj::na(xh[static_cast<std::size_t>(j)]);
    static const char *const lv[] = {"below_lod", "below_loq", "quantifiable"};
    f.put("flag", fac(fl.data(), nu, lv, 3));
    SET_VECTOR_ELT(out, 1, f.done());
    UNPROTECT(1);
  }
  SET_VECTOR_ELT(out, 2, Rf_mkString(wl));
  static const char *const nm[] = {"fit", "unk", "w"};
  UNPROTECT(5);
  return bj::named(out, nm);
}

SEXP C_clab_r_lodimp(SEXP x, SEXP lod, SEXP meth) {
  clab::who = "lodimp";
  SEXP xv = PROTECT(bfx::nna(x, "x", -1));
  const int n = bj::len(xv, "x");
  SEXP lv = PROTECT(bj::nums(lod, "lod", -1, 0.0, true));
  const int nl = bj::len(lv, "lod");
  if (nl != 1 && nl != n) fail("`lod` must be one limit, or one per value, %d.", n);
  static const char *const mo[] = {"half", "sqrt2", "lod"};
  const int m = csp::pick(meth, "meth", mo, 3);
  bj::cols f(4, n);
  f.put("x", xv);
  double *l = f.dbl("lod");
  for (int i = 0; i < n; ++i) l[i] = REAL(lv)[nl == 1 ? 0 : i];
  double *o = f.dbl("imp"), st[2];
  int *c = f.lgl("cens");
  clab_lodimp(REAL(xv), n, REAL(lv), nl, m, o, c, st);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, f.done());
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(bj::na(st[0])));
  SET_VECTOR_ELT(out, 2, Rf_ScalarInteger(static_cast<int>(st[1])));
  static const char *const nm[] = {"tab", "frac", "nmiss"};
  UNPROTECT(4);
  return bj::named(out, nm);
}

SEXP C_clab_r_pmrrat(SEXP cen, SEXP per, SEXP cut) {
  clab::who = "pmrrat";
  SEXP cv = PROTECT(bj::nums(cen, "cen", -1, 0.0, true));
  const int n = bj::len(cv, "cen");
  SEXP pv = PROTECT(bj::nums(per, "per", n, 0.0, true));
  SEXP kv = PROTECT(bj::nums(cut, "cut", 2, 0.0, true));
  if (REAL(kv)[1] < REAL(kv)[0]) fail("`cut` must be two increasing ratios.");
  bj::cols f(4, n);
  f.put("cen", cv);
  f.put("per", pv);
  double *r = f.dbl("ratio");
  std::vector<int> k(static_cast<std::size_t>(n));
  clab_pmrrat(REAL(cv), REAL(pv), n, REAL(kv)[0], REAL(kv)[1], r, k.data());
  static const char *const lv[] = {"minimal", "modest", "significant"};
  f.put("redist", fac(k.data(), n, lv, 3));
  SEXP out = f.done();
  UNPROTECT(4);
  return out;
}

SEXP C_clab_r_etohlr(SEXP eth, SEXP prop, SEXP but, SEXP etg, SEXP thr) {
  clab::who = "etohlr";
  SEXP ev = PROTECT(bj::nums(eth, "eth", -1, 0.0));
  const int n = bj::len(ev, "eth");
  SEXP pv = PROTECT(bfx::nna(prop, "prop", n));
  SEXP bv = PROTECT(bfx::nna(but, "but", n));
  SEXP gv = PROTECT(bfx::nna(etg, "etg", n));
  const double th = bj::inr(thr, "thr", 0.0, R_PosInf);
  bj::cols f(6, n);
  f.put("eth", ev);
  f.put("prop", pv);
  f.put("but", bv);
  f.put("etg", gv);
  double *hi = f.dbl("higher");
  std::vector<int> v(static_cast<std::size_t>(n));
  clab_etohlr(REAL(ev), REAL(pv), REAL(bv), REAL(gv), n, th, hi, v.data());
  for (int i = 0; i < n; ++i) hi[i] = bj::na(hi[i]);
  static const char *const lv[] = {"antemortem", "postmortem", "indeterminate"};
  f.put("verdict", fac(v.data(), n, lv, 3));
  SEXP out = f.done();
  UNPROTECT(5);
  return out;
}

SEXP C_clab_r_mtxrel(SEXP base, SEXP dil, SEXP days, SEXP stage, SEXP dmax, SEXP dcap,
                     SEXP smax, SEXP scap) {
  clab::who = "mtxrel";
  SEXP bv = PROTECT(bj::nums(base, "base", -1, 0.0));
  const int n = bj::len(bv, "base");
  clab::lgl(dil, "dil", n);
  const double dy = bj::inr(days, "days", 0.0, R_PosInf), sg = bj::inr(stage, "stage", 0.0, R_PosInf);
  const double dm = clab::one(dmax, "dmax", 0.0, R_PosInf), sm = clab::one(smax, "smax", 0.0, R_PosInf);
  const double dc = clab::p01(dcap, "dcap"), sc = clab::p01(scap, "scap");
  std::vector<double> rel(static_cast<std::size_t>(n));
  std::vector<int> ord(static_cast<std::size_t>(n));
  clab_mtxrel(REAL(bv), LOGICAL(dil), n, dy, sg, dm, dc, sm, sc, rel.data(), ord.data());
  SEXP nm = Rf_getAttrib(base, R_NamesSymbol);
  SEXP lab = PROTECT(Rf_isNull(nm) ? clab::seqs(n) : nm);
  bj::cols f(5, n);
  SEXP mx = PROTECT(Rf_allocVector(STRSXP, n));
  f.put("matrix", mx);
  double *b = f.dbl("base"), *r = f.dbl("rel");
  int *d = f.lgl("dil"), *rk = f.itg("rank");
  for (int j = 0; j < n; ++j) {
    const int i = ord[static_cast<std::size_t>(j)];
    SET_STRING_ELT(mx, j, STRING_ELT(lab, i));
    b[j] = REAL(bv)[i];
    r[j] = rel[static_cast<std::size_t>(i)];
    d[j] = LOGICAL(dil)[i];
    rk[j] = j + 1;
  }
  SEXP out = f.done();
  UNPROTECT(4);
  return out;
}

}  // extern "C"
