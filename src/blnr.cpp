// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the eyewitness lineup methods (lnroc, lnsize;
// src/lnr.cpp): each checks its arguments, calls the kernel and builds
// the R result. Nothing calls back into the R evaluator.
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstring>
#include <vector>

#include "bcol.h"
#include "kutil.h"
#include "rchk.h"

extern "C" double Rf_pnorm5(double, double, double, int, int);

using clab::fail;

namespace {

// The procedure label of each row of a procedure-major table of np x nl.
SEXP prep(SEXP lev, int np, int nl) {
  SEXP s = PROTECT(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(np) * nl));
  for (int p = 0; p < np; ++p) {
    for (int l = 0; l < nl; ++l) SET_STRING_ELT(s, static_cast<R_xlen_t>(p) * nl + l, STRING_ELT(lev, p));
  }
  UNPROTECT(1);
  return s;
}

}  // namespace

extern "C" {

SEXP C_clab_r_lnroc(SEXP tp, SEXP pick, SEXP conf, SEXP proc, SEXP size, SEXP far, SEXP nboot, SEXP seed,
                    SEXP level) {
  clab::who = "lnroc";
  int np0 = 0;
  SEXP tv = PROTECT(csp::bin01(tp, "tp", XLENGTH(tp)));
  ++np0;
  const int n = bj::len(tv, "tp");
  if (n < 2) fail("`tp` must hold two witnesses at least.");
  csp::lab pl = bj::labs(pick, "pick", n);
  np0 += 2;
  std::vector<int> pk(static_cast<std::size_t>(n));
  {
    std::vector<int> map(static_cast<std::size_t>(pl.g));
    for (int j = 0; j < pl.g; ++j) {
      const char *s = bj::str(pl.lev, j);
      map[static_cast<std::size_t>(j)] =
          !std::strcmp(s, "none") ? 0 : !std::strcmp(s, "suspect") ? 1 : !std::strcmp(s, "filler") ? 2 : -1;
      if (map[static_cast<std::size_t>(j)] < 0) {
        fail("`pick` must be \"suspect\", \"filler\" or \"none\" for every witness, not %s.", s);
      }
    }
    for (int i = 0; i < n; ++i) pk[static_cast<std::size_t>(i)] = map[static_cast<std::size_t>(INTEGER(pl.code)[i])];
  }
  clab::unlab(conf, "conf");
  clab::vec(conf, "conf");
  // all NA (no identification has a confidence) comes in as logical
  bool nal = TYPEOF(conf) == LGLSXP;
  for (R_xlen_t i = 0; nal && i < XLENGTH(conf); ++i) nal = LOGICAL(conf)[i] == NA_LOGICAL;
  if (!(clab::isnum(conf) || nal) || XLENGTH(conf) != n) fail("`conf` must be numbers, one per witness, %d.", n);
  SEXP cv = PROTECT(nal ? Rf_coerceVector(conf, REALSXP) : clab::dbl(conf));
  ++np0;
  for (int i = 0; i < n; ++i) {
    if (pk[static_cast<std::size_t>(i)] != 0 && !std::isfinite(REAL(cv)[i])) {
      fail("`conf` must be a finite number for every identification.");
    }
  }
  SEXP zero = PROTECT(Rf_allocVector(INTSXP, n));
  ++np0;
  for (int i = 0; i < n; ++i) INTEGER(zero)[i] = 0;
  const int *pc = INTEGER(zero);
  int np = 1;
  SEXP plev = R_NilValue;
  if (!Rf_isNull(proc)) {
    csp::lab pr = bj::labs(proc, "proc", n);
    np0 += 2;
    pc = INTEGER(pr.code);
    np = pr.g;
    plev = pr.lev;
  }
  const int sz = Rf_isNull(size) ? 0 : clab::whole(size, "size", 2);
  const double fm = Rf_isNull(far) ? 0.0 : bj::inr(far, "far", 0.0, 1.0);
  if (!Rf_isNull(far) && !(fm > 0)) fail("`far` must be one number above 0, at most 1.");
  const int nb = clab::whole(nboot, "nboot", 0);
  // the bootstrap areas are kept, one per replicate and procedure
  if (static_cast<double>(nb) * np > kutil::maxcell) fail("`nboot` times the procedures must be at most %.0f.", kutil::maxcell);
  const unsigned int sd = nb > 0 ? csp::seed(seed) : 0u;
  const double lv = clab::one(level, "level", 0.0, 1.0);
  std::vector<int> ti(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) ti[static_cast<std::size_t>(i)] = INTEGER(tv)[i];
  {
    std::vector<int> t1(static_cast<std::size_t>(np), 0), t0(static_cast<std::size_t>(np), 0);
    for (int i = 0; i < n; ++i) (ti[static_cast<std::size_t>(i)] ? t1 : t0)[static_cast<std::size_t>(pc[i])]++;
    for (int p = 0; p < np; ++p) {
      if (!t1[static_cast<std::size_t>(p)] || !t0[static_cast<std::size_t>(p)]) {
        fail("every procedure needs target-present and target-absent lineups.");
      }
    }
    bool any = false;
    for (int v : pk) any = any || v != 0;
    if (!any) fail("no witness made an identification.");
  }
  std::vector<double> lev(static_cast<std::size_t>(n));
  int nl = 0;
  std::vector<double> roc(3 * static_cast<std::size_t>(np) * n), cac(3 * static_cast<std::size_t>(np) * n);
  std::vector<double> pa(static_cast<std::size_t>(np)), bse(2 * static_cast<std::size_t>(np));
  double fu = 0.0;
  const int rc = clab_lnroc(ti.data(), pk.data(), REAL(cv), pc, n, np, sz, fm, nb, sd, lev.data(), &nl,
                            roc.data(), cac.data(), pa.data(), bse.data(), &fu);
  clab::ichk(rc);
  bj::kc(rc, "the false alarm rates are all 0: no partial area.");
  const std::size_t un = static_cast<std::size_t>(np);
  const int m = np * nl;
  const bool hp = !Rf_isNull(plev);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 6));
  ++np0;
  for (int t = 0; t < 2; ++t) {
    bj::cols f(6, m);
    if (hp) f.put("proc", prep(plev, np, nl));
    double *cf = f.dbl("conf");
    for (int k = 0; k < m; ++k) cf[k] = lev[static_cast<std::size_t>(k % nl)];
    const double *src = t == 0 ? roc.data() : cac.data();
    static const char *const rn[] = {"hr", "far", "dr"}, *const cn[] = {"hr", "far", "acc"};
    for (int q = 0; q < 3; ++q) {
      double *c = f.dbl(t == 0 ? rn[q] : cn[q]);
      for (int k = 0; k < m; ++k) c[k] = bj::na(src[static_cast<std::size_t>(q) * m + k]);
    }
    SET_VECTOR_ELT(out, t, f.done());
  }
  const double z = Rf_qnorm5(0.5 + 0.5 * lv, 0.0, 1.0, 1, 0);
  bj::cols s(8, np);
  if (hp) s.put("proc", plev);
  double *a = s.dbl("pauc"), *se = s.dbl("se"), *lo = s.dbl("lo"), *hi = s.dbl("hi");
  double *df = s.dbl("diff"), *ds = s.dbl("dse"), *pv = s.dbl("p");
  for (int p = 0; p < np; ++p) {
    a[p] = pa[static_cast<std::size_t>(p)];
    se[p] = bj::na(bse[static_cast<std::size_t>(p)]);
    lo[p] = bj::na(a[p] - z * bse[static_cast<std::size_t>(p)]);
    hi[p] = bj::na(a[p] + z * bse[static_cast<std::size_t>(p)]);
    df[p] = a[p] - pa[0];
    ds[p] = bj::na(bse[un + p]);
    const double zz = df[p] / bse[un + p];
    pv[p] = std::isnan(zz) ? NA_REAL : 2 * Rf_pnorm5(-std::fabs(zz), 0.0, 1.0, 1, 0);
  }
  SET_VECTOR_ELT(out, 2, s.done());
  SET_VECTOR_ELT(out, 3, Rf_ScalarReal(fu));
  SET_VECTOR_ELT(out, 4, sz ? Rf_ScalarInteger(sz) : R_NilValue);
  SET_VECTOR_ELT(out, 5, Rf_ScalarInteger(nb));
  static const char *const on[] = {"roc", "cac", "pauc", "far", "size", "nboot"};
  bj::named(out, on);
  Rf_classgets(out, PROTECT(Rf_mkString("lnroc")));
  UNPROTECT(np0 + 4);
  return out;
}

SEXP C_clab_r_lnsize(SEXP choice, SEXP k, SEXP susp, SEXP level) {
  clab::who = "lnsize";
  SEXP cv = PROTECT(bj::nums(choice, "choice", -1, 0.0));
  const int n = bj::len(cv, "choice");
  int mx = 0;
  for (int i = 0; i < n; ++i) {
    const double v = REAL(cv)[i];
    if (v != std::floor(v) || v > INT_MAX) {
      fail("`choice` must be whole numbers from 0: lineup positions, 0 for no choice.");
    }
    mx = std::max(mx, static_cast<int>(v));
  }
  if (mx == 0) fail("`choice` must name a lineup member at least once.");
  const int kk = Rf_isNull(k) ? std::max(mx, 2) : clab::whole(k, "k", 2);
  if (kk > kutil::maxcell) fail("`k` must be at most %.0f.", kutil::maxcell);
  if (mx > kk) fail("`choice` names position %d of a lineup of %d.", mx, kk);
  const int sp = clab::whole(susp, "susp", 1);
  if (sp > kk) fail("`susp` must be a position from 1 to %d.", kk);
  const double lv = clab::one(level, "level", 0.0, 1.0);
  std::vector<int> ch(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) ch[static_cast<std::size_t>(i)] = static_cast<int>(REAL(cv)[i]) - 1;
  SEXP cnt = PROTECT(Rf_allocVector(REALSXP, kk));
  double o[11];
  bj::kc(clab_lnsize(ch.data(), n, kk, sp - 1, Rf_qnorm5(0.5 + 0.5 * lv, 0.0, 1.0, 1, 0), REAL(cnt), o),
         "the kernel refused the input.");
  static const char *const sn[] = {"n",      "k",   "nm",    "psusp", "pbin", "fsize",
                                   "esize", "esizeb", "etr", "etrlo", "etrhi"};
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, cnt);
  SET_VECTOR_ELT(out, 1, clab::nvec(11, sn, o));
  SET_VECTOR_ELT(out, 2, Rf_ScalarInteger(sp));
  static const char *const on[] = {"count", "stat", "susp"};
  bj::named(out, on);
  Rf_classgets(out, PROTECT(Rf_mkString("lnsize")));
  UNPROTECT(4);
  return out;
}

}  // extern "C"
