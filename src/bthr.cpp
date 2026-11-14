// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry point of the tests of discrimination in searches (thrtst,
// src/thr.cpp): it checks the arguments, calls the kernels and builds the
// R result. Nothing calls back into the R evaluator.
#include <cmath>
#include <vector>

#include "bcol.h"
#include "rchk.h"

extern "C" double Rf_pnorm5(double, double, double, int, int);

using clab::fail;

namespace {

SEXP labsel(SEXP lev, const std::vector<int> &code) {
  SEXP s = PROTECT(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(code.size())));
  for (std::size_t i = 0; i < code.size(); ++i) {
    SET_STRING_ELT(s, static_cast<R_xlen_t>(i), STRING_ELT(lev, code[i]));
  }
  UNPROTECT(1);
  return s;
}

}  // namespace

extern "C" {

SEXP C_clab_r_thrtst(SEXP stops, SEXP srch, SEXP hits, SEXP grp, SEXP area, SEXP off, SEXP ref,
                     SEXP level, SEXP maxit) {
  clab::who = "thrtst";
  int np = 0;
  SEXP sv = PROTECT(bj::nums(stops, "stops", -1, 0.0));
  ++np;
  const int n = bj::len(sv, "stops");
  SEXP kv = PROTECT(bj::nums(srch, "srch", n, 0.0));
  SEXP hv = PROTECT(bj::nums(hits, "hits", n, 0.0));
  np += 2;
  for (int i = 0; i < n; ++i) {
    if (REAL(kv)[i] > REAL(sv)[i] || REAL(hv)[i] > REAL(kv)[i]) {
      fail("`srch` must be at most `stops`, and `hits` at most `srch`, in every row.");
    }
  }
  csp::lab g = bj::labs(grp, "grp", n);
  np += 2;
  if (g.g < 2) fail("`grp` must hold at least two groups.");
  SEXP zero = PROTECT(Rf_allocVector(INTSXP, n));
  ++np;
  for (int i = 0; i < n; ++i) INTEGER(zero)[i] = 0;
  const int *ac = INTEGER(zero);
  int na = 1;
  SEXP alev = R_NilValue;
  if (!Rf_isNull(area)) {
    csp::lab a = bj::labs(area, "area", n);
    np += 2;
    ac = INTEGER(a.code);
    na = a.g;
    alev = a.lev;
  }
  const int rf = Rf_isNull(ref) ? 1 : clab::ref(ref, g.lev, g.g);
  const double lv = clab::one(level, "level", 0.0, 1.0);
  const int mi = clab::whole(maxit, "maxit", 1);
  const int ng = g.g;
  std::vector<double> tk(static_cast<std::size_t>(ng), 0.0);
  for (int i = 0; i < n; ++i) tk[static_cast<std::size_t>(INTEGER(g.code)[i])] += REAL(kv)[i];
  for (int j = 0; j < ng; ++j) {
    if (!(tk[static_cast<std::size_t>(j)] > 0)) {
      fail("every group needs searches: %s has none.", bj::str(g.lev, j));
    }
  }
  // the threshold model only where it is identified: more rates than
  // parameters
  std::vector<char> seen(static_cast<std::size_t>(ng) * na, 0);
  int nc = 0;
  for (int i = 0; i < n; ++i) {
    char &s = seen[static_cast<std::size_t>(ac[i]) * ng + INTEGER(g.code)[i]];
    if (REAL(sv)[i] > 0 && !s) {
      s = 1;
      ++nc;
    }
  }
  const bool fit = 2 * nc > 3 * (ng + na - 1);
  const std::size_t un = static_cast<std::size_t>(ng), cn = static_cast<std::size_t>(ng) * na;
  std::vector<double> gs(6 * un), hz(4 * un), chi(3), thr(un), tse(3 * un), cfit(5 * cn), lrt(6);
  std::vector<int> conv(4, 0);
  const int rc = clab_thrtst(REAL(sv), REAL(kv), REAL(hv), INTEGER(g.code), ac, n, ng, na, rf - 1,
                             fit ? mi : 0, gs.data(), hz.data(), chi.data(), thr.data(), tse.data(),
                             cfit.data(), lrt.data(), conv.data());
  clab::ichk(rc);
  bj::kc(rc, "the threshold model could not be fitted to these counts.");
  const double z = Rf_qnorm5(0.5 + 0.5 * lv, 0.0, 1.0, 1, 0);
  bj::cols h(10, ng);
  h.put("group", g.lev);
  static const char *const hn[] = {"stops", "srch", "hits", "srate", "hrate"};
  for (int q = 0; q < 5; ++q) {
    double *c = h.dbl(hn[q]);
    for (int j = 0; j < ng; ++j) c[j] = bj::na(gs[static_cast<std::size_t>(q) * un + j]);
  }
  static const char *const zn[] = {"diff", "se", "z", "p"};
  for (int q = 0; q < 4; ++q) {
    double *c = h.dbl(zn[q]);
    for (int j = 0; j < ng; ++j) c[j] = bj::na(hz[static_cast<std::size_t>(q) * un + j]);
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 8));
  ++np;
  SET_VECTOR_ELT(out, 0, h.done());
  static const char *const kn[] = {"chisq", "df", "p"};
  SET_VECTOR_ELT(out, 1, clab::nvec(3, kn, chi.data()));
  if (fit) {
    bj::cols t(9, ng);
    t.put("group", g.lev);
    double *tv = t.dbl("thr"), *se = t.dbl("se"), *lo = t.dbl("lo"), *hi = t.dbl("hi");
    double *df = t.dbl("diff"), *ds = t.dbl("dse"), *dz = t.dbl("z"), *dp = t.dbl("p");
    for (int j = 0; j < ng; ++j) {
      tv[j] = thr[static_cast<std::size_t>(j)];
      se[j] = bj::na(tse[static_cast<std::size_t>(j)]);
      lo[j] = bj::na(tv[j] - z * tse[static_cast<std::size_t>(j)]);
      hi[j] = bj::na(tv[j] + z * tse[static_cast<std::size_t>(j)]);
      df[j] = tse[un + j];
      ds[j] = bj::na(tse[2 * un + j]);
      const double zz = df[j] / tse[2 * un + j];
      dz[j] = j == rf - 1 ? NA_REAL : bj::na(zz);
      dp[j] = j == rf - 1 || std::isnan(zz) ? NA_REAL : 2 * Rf_pnorm5(-std::fabs(zz), 0.0, 1.0, 1, 0);
    }
    SET_VECTOR_ELT(out, 2, t.done());
    static const char *const ln[] = {"stat", "df", "p", "llfull", "llsame", "cells"};
    SET_VECTOR_ELT(out, 3, clab::nvec(6, ln, lrt.data()));
    bj::cols f(3, 2);
    SEXP mn = PROTECT(Rf_allocVector(STRSXP, 2));
    ++np;
    SET_STRING_ELT(mn, 0, Rf_mkChar("full"));
    SET_STRING_ELT(mn, 1, Rf_mkChar("same"));
    f.put("model", mn);
    int *cv = f.lgl("conv"), *it = f.itg("iter");
    for (int r = 0; r < 2; ++r) {
      cv[r] = conv[static_cast<std::size_t>(r)] == 1;
      it[r] = conv[static_cast<std::size_t>(2 + r)];
    }
    SET_VECTOR_ELT(out, 4, f.done());
    std::vector<int> cg, ca;
    for (int q = 0; q < na; ++q) {
      for (int j = 0; j < ng; ++j) {
        if (seen[static_cast<std::size_t>(q) * ng + j]) {
          cg.push_back(j);
          ca.push_back(q);
        }
      }
    }
    const int m = static_cast<int>(cg.size());
    bj::cols e(7, m);
    e.put("group", labsel(g.lev, cg));
    if (!Rf_isNull(alev)) e.put("area", labsel(alev, ca));
    static const char *const en[] = {"phi", "delta", "thr", "srate", "hrate"};
    for (int q = 0; q < 5; ++q) {
      double *c = e.dbl(en[q]);
      for (int k = 0; k < m; ++k) {
        c[k] = bj::na(cfit[static_cast<std::size_t>(q) * cn + static_cast<std::size_t>(ca[k]) * ng + cg[k]]);
      }
    }
    SET_VECTOR_ELT(out, 5, e.done());
  }
  if (!Rf_isNull(off)) {
    csp::lab o = bj::labs(off, "off", n);
    np += 2;
    if (o.g < 2) fail("`off` must hold at least two officer groups.");
    const int no = o.g, npair = no * (no - 1);
    std::vector<double> rate(5 * static_cast<std::size_t>(no) * ng), res(3 * static_cast<std::size_t>(npair));
    const int r2 = clab_thrrk(REAL(sv), REAL(kv), REAL(hv), INTEGER(o.code), INTEGER(g.code), n, no, ng, z,
                              rate.data(), res.data());
    clab::ichk(r2);
    bj::kc(r2, "the kernel refused the input.");
    std::vector<int> a1, a2;
    for (int x = 0; x < no; ++x) {
      for (int y = x + 1; y < no; ++y) {
        for (int ms = 0; ms < 2; ++ms) {
          a1.push_back(x);
          a2.push_back(y);
        }
      }
    }
    bj::cols r(6, npair);
    r.put("off1", labsel(o.lev, a1));
    r.put("off2", labsel(o.lev, a2));
    SEXP me = PROTECT(Rf_allocVector(STRSXP, npair));
    ++np;
    for (int k = 0; k < npair; ++k) SET_STRING_ELT(me, k, Rf_mkChar(k % 2 ? "hit" : "search"));
    r.put("meas", me);
    double *lo = r.dbl("zmin"), *hi = r.dbl("zmax");
    int *rj = r.lgl("rej");
    for (int k = 0; k < npair; ++k) {
      lo[k] = bj::na(res[static_cast<std::size_t>(k)]);
      hi[k] = bj::na(res[static_cast<std::size_t>(npair + k)]);
      rj[k] = res[static_cast<std::size_t>(2 * npair + k)] > 0;
    }
    SET_VECTOR_ELT(out, 6, r.done());
  }
  SEXP rl = PROTECT(Rf_allocVector(STRSXP, 1));
  ++np;
  SET_STRING_ELT(rl, 0, STRING_ELT(g.lev, rf - 1));
  SET_VECTOR_ELT(out, 7, rl);
  static const char *const on[] = {"hit", "kpt", "thr", "lrt", "fit", "cells", "rank", "ref"};
  bj::named(out, on);
  Rf_classgets(out, PROTECT(Rf_mkString("thrtst")));
  ++np;
  const int ncol = 1 + (fit ? 3 : 0) + (Rf_isNull(off) ? 0 : 1);
  UNPROTECT(np + ncol);
  return out;
}

}  // extern "C"
