// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the courts and sentencing family (jdgiv, jmono,
// mmshr): each checks its arguments, calls the kernels and builds the R
// result. Normal tails and chi-square tails come from R's maths library.
// Nothing calls back into the R evaluator.
#include <algorithm>
#include <cmath>
#include <vector>

#include "bcol.h"
#include "rchk.h"

extern "C" double Rf_pnorm5(double, double, double, int, int);
extern "C" double Rf_pchisq(double, double, int, int);

using clab::fail;

namespace {

// The decisions as doubles 0 and 1.
SEXP dec(SEXP x, const char *nm, R_xlen_t n) {
  clab::unlab(x, nm);
  clab::vec(x, nm);
  SEXP i = PROTECT(csp::bin01(x, nm, n));
  SEXP d = PROTECT(Rf_allocVector(REALSXP, n));
  for (R_xlen_t k = 0; k < n; ++k) REAL(d)[k] = INTEGER(i)[k];
  UNPROTECT(2);
  return d;
}

// One column of a two-column matrix or data frame, as doubles.
SEXP col(SEXP x, int j, const char *nm) {
  if (TYPEOF(x) == VECSXP) {
    if (XLENGTH(x) != 2) fail("`%s` must have two columns.", nm);
    return VECTOR_ELT(x, j);
  }
  SEXP dm = Rf_getAttrib(x, R_DimSymbol);
  if (TYPEOF(dm) != INTSXP || XLENGTH(dm) != 2 || INTEGER(dm)[1] != 2) {
    fail("`%s` must have two columns.", nm);
  }
  if (TYPEOF(x) != REALSXP && TYPEOF(x) != LGLSXP && TYPEOF(x) != INTSXP) {
    fail("`%s` must hold numbers.", nm);
  }
  const int n = INTEGER(dm)[0];
  SEXP v = PROTECT(Rf_allocVector(TYPEOF(x), n));
  for (int i = 0; i < n; ++i) {
    if (TYPEOF(x) == REALSXP)
      REAL(v)[i] = REAL(x)[i + static_cast<R_xlen_t>(n) * j];
    else if (TYPEOF(x) == LGLSXP)
      LOGICAL(v)[i] = LOGICAL(x)[i + static_cast<R_xlen_t>(n) * j];
    else
      INTEGER(v)[i] = INTEGER(x)[i + static_cast<R_xlen_t>(n) * j];
  }
  UNPROTECT(1);
  return v;
}

}  // namespace

extern "C" {

SEXP C_clab_r_jdgiv(SEXP z, SEXP d, SEXP y, SEXP w, SEXP dsh, SEXP deff) {
  clab::who = "jdgiv";
  int np = 0;
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  ++np;
  if (Rf_isNull(z)) {
    // a known population: d = (d0, d1), y = (y0, y1)
    SEXP d0 = PROTECT(col(d, 0, "d")), d1 = PROTECT(col(d, 1, "d"));
    SEXP y0 = PROTECT(col(y, 0, "y")), y1 = PROTECT(col(y, 1, "y"));
    np += 4;
    const int n = bj::len(d0, "d");
    SEXP a = PROTECT(dec(d0, "d", n)), b = PROTECT(dec(d1, "d", n));
    SEXP u = PROTECT(bj::nums(y0, "y", n)), v = PROTECT(bj::nums(y1, "y", n));
    np += 4;
    SEXP wv = R_NilValue;
    if (!Rf_isNull(w)) {
      wv = PROTECT(bj::nums(w, "w", n, 0.0));
      ++np;
    }
    double o[12];
    if (clab_jdgpop(REAL(a), REAL(b), REAL(u), REAL(v), Rf_isNull(wv) ? nullptr : REAL(wv), n, o) !=
        CLAB_OK) {
      fail("`w` must have a positive total.");
    }
    static const char *const sn[] = {"complier", "defier", "always", "never"};
    static const char *const en[] = {"itt", "fs", "wald", "ate"};
    SET_VECTOR_ELT(out, 0, Rf_mkString("pop"));
    SET_VECTOR_ELT(out, 1, clab::nvec(4, en, o + 8));
    SET_VECTOR_ELT(out, 2, clab::nvec(4, sn, o));
    SEXP ef = PROTECT(clab::nvec(4, sn, o + 4));
    ++np;
    for (int k = 0; k < 4; ++k) REAL(ef)[k] = bj::na(REAL(ef)[k]);
    SET_VECTOR_ELT(out, 3, ef);
    static const char *const nm[] = {"type", "est", "shares", "effects"};
    UNPROTECT(np);
    return bj::named(out, nm);
  }
  const int n = bj::len(z, "z");
  SEXP zd = PROTECT(dec(z, "z", n)), dd = PROTECT(dec(d, "d", n));
  SEXP yv = PROTECT(bj::nums(y, "y", n));
  np += 3;
  SEXP wv = R_NilValue;
  if (!Rf_isNull(w)) {
    wv = PROTECT(bj::nums(w, "w", n, 0.0));
    ++np;
  }
  SEXP sv = PROTECT(bj::nums(dsh, "dsh", -1, 0.0)), ev = PROTECT(bj::nums(deff, "deff"));
  np += 2;
  const int nd = bj::len(sv, "dsh"), ne = bj::len(ev, "deff");
  for (int k = 0; k < nd; ++k) {
    if (REAL(sv)[k] > 1) fail("`dsh` must be shares from 0 to 1.");
  }
  double est[6];
  SEXP sens = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(nd) * ne));
  ++np;
  if (clab_jdgiv(REAL(zd), REAL(dd), REAL(yv), Rf_isNull(wv) ? nullptr : REAL(wv), n, REAL(sv), nd,
                 REAL(ev), ne, est, REAL(sens)) != CLAB_OK) {
    fail("`z` must take both values, each with a positive weight.");
  }
  static const char *const en[] = {"fs", "itt", "wald"};
  static const char *const sn[] = {"complier", "always", "never"};
  SEXP e = PROTECT(clab::nvec(3, en, est));
  ++np;
  REAL(e)[2] = bj::na(REAL(e)[2]);
  bj::cols g(4, nd * ne);
  ++np;
  double *gs = g.dbl("dsh"), *ge = g.dbl("deff"), *gc = g.dbl("cshare"), *gt = g.dbl("ceff");
  for (int b = 0; b < ne; ++b) {
    for (int a = 0; a < nd; ++a) {
      const int k = a + nd * b;
      gs[k] = REAL(sv)[a];
      ge[k] = REAL(ev)[b];
      gc[k] = est[3] + REAL(sv)[a];
      gt[k] = bj::na(REAL(sens)[k]);
    }
  }
  SET_VECTOR_ELT(out, 0, Rf_mkString("obs"));
  SET_VECTOR_ELT(out, 1, e);
  SET_VECTOR_ELT(out, 2, clab::nvec(3, sn, est + 3));
  SET_VECTOR_ELT(out, 3, g.done());
  static const char *const nm[] = {"type", "est", "shares", "sens"};
  UNPROTECT(np);
  return bj::named(out, nm);
}

SEXP C_clab_r_jmono(SEXP judge, SEXP d, SEXP y, SEXP grp, SEXP w, SEXP lo, SEXP hi, SEXP minn,
                    SEXP adj, SEXP alpha) {
  clab::who = "jmono";
  int np = 0;
  const int n = bj::len(d, "d");
  SEXP dd = PROTECT(dec(d, "d", n));
  ++np;
  csp::lab j = bj::labs(judge, "judge", n);
  np += 2;
  const int nj = j.g;
  SEXP yv = R_NilValue, wv = R_NilValue;
  if (!Rf_isNull(y)) {
    yv = PROTECT(bj::nums(y, "y", n));
    ++np;
  }
  if (!Rf_isNull(w)) {
    wv = PROTECT(bj::nums(w, "w", n, 0.0));
    ++np;
  }
  const double mn = bj::inr(minn, "minn", 0.0, R_PosInf);
  static const char *const am[] = {"none", "bonferroni", "holm", "hochberg", "BH", "BY"};
  const int ad = csp::pick(adj, "adj", am, 6);
  const double al = clab::one(alpha, "alpha", 0.0, 1.0);
  SEXP jw = PROTECT(Rf_allocVector(REALSXP, 4 * static_cast<R_xlen_t>(nj) + 4 + nj));
  ++np;
  double *cnt = REAL(jw), *one = cnt + nj, *rate = one + nj, *outc = rate + nj, *st = outc + nj;
  if (clab_jrate(INTEGER(j.code), REAL(dd), Rf_isNull(yv) ? nullptr : REAL(yv),
                 Rf_isNull(wv) ? nullptr : REAL(wv), n, nj, cnt, one, rate, outc, st) != CLAB_OK) {
    fail("every judge needs cases with a positive total weight.");
  }
  // the judges, most lenient (lowest rate) first
  std::vector<int> o(static_cast<std::size_t>(nj));
  for (int k = 0; k < nj; ++k) o[static_cast<std::size_t>(k)] = k;
  std::stable_sort(o.begin(), o.end(), [&](int a, int b) { return rate[a] < rate[b]; });
  bj::cols f(5, nj);
  ++np;
  SEXP jl = PROTECT(Rf_allocVector(STRSXP, nj));
  ++np;
  for (int k = 0; k < nj; ++k)
    SET_STRING_ELT(jl, k, STRING_ELT(j.lev, o[static_cast<std::size_t>(k)]));
  f.put("judge", jl);
  double *fn = f.dbl("n"), *fr = f.dbl("rate"), *fo = Rf_isNull(yv) ? nullptr : f.dbl("outcome");
  SEXP so = PROTECT(Rf_allocVector(REALSXP, 2 * static_cast<R_xlen_t>(nj)));
  ++np;
  double *sr = REAL(so), *sy = sr + nj;
  for (int k = 0; k < nj; ++k) {
    const int q = o[static_cast<std::size_t>(k)];
    fn[k] = cnt[q];
    fr[k] = sr[k] = rate[q];
    sy[k] = outc[q];
    if (fo) fo[k] = outc[q];
  }
  o.clear();
  o.shrink_to_fit();
  const double x2 = st[0], df = st[1];
  const double hv[5] = {x2, df, std::isnan(x2) ? NA_REAL : Rf_pchisq(x2, df, 0, 0), st[2], st[3]};
  static const char *const hn[] = {"chisq", "df", "p", "min", "max"};
  SEXP hom = PROTECT(clab::nvec(5, hn, hv));
  ++np;
  REAL(hom)[0] = bj::na(REAL(hom)[0]);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
  ++np;
  SET_VECTOR_ELT(out, 0, f.done());
  SET_VECTOR_ELT(out, 1, hom);
  if (!Rf_isNull(yv)) {
    double l = INFINITY, h = -INFINITY;
    for (int i = 0; i < n; ++i) {
      l = std::fmin(l, REAL(yv)[i]);
      h = std::fmax(h, REAL(yv)[i]);
    }
    if (!Rf_isNull(lo)) {
      const double v = clab::one(lo, "lo", R_NegInf, R_PosInf);
      if (v > l) fail("`lo` must be at most the smallest outcome.");
      l = v;
    }
    if (!Rf_isNull(hi)) {
      const double v = clab::one(hi, "hi", R_NegInf, R_PosInf);
      if (v < h) fail("`hi` must be at least the largest outcome.");
      h = v;
    }
    const long long npr = static_cast<long long>(nj) * (nj - 1) / 2;
    if (npr > 50000000LL) fail("`judge` has too many judges for every pair (more than 1e4).");
    const int m = static_cast<int>(npr);
    bj::cols p(7, m);
    ++np;
    SEXP ja = PROTECT(Rf_allocVector(STRSXP, m)), jb = PROTECT(Rf_allocVector(STRSXP, m));
    np += 2;
    for (int a = 0, k = 0; a < nj; ++a) {
      for (int b = a + 1; b < nj; ++b, ++k) {
        SET_STRING_ELT(ja, k, STRING_ELT(jl, a));
        SET_STRING_ELT(jb, k, STRING_ELT(jl, b));
      }
    }
    p.put("j", ja);
    p.put("k", jb);
    double *dp = p.dbl("dp"), *dy = p.dbl("dy"), *bn = p.dbl("bound");
    int *vi = p.lgl("viol");
    double *lt = p.dbl("late");
    clab_jslope(sr, sy, nj, l, h, dp, dy, bn, vi, lt);
    for (int k = 0; k < m; ++k) lt[k] = bj::na(lt[k]);
    SET_VECTOR_ELT(out, 2, p.done());
  }
  if (!Rf_isNull(grp)) {
    csp::lab g = bj::labs(grp, "grp", n);
    np += 2;
    if (static_cast<long long>(nj) * g.g > 100000000LL)
      fail("`judge` and `grp` make more than 1e8 cells.");
    SEXP cw = PROTECT(Rf_allocVector(REALSXP, 2 * static_cast<R_xlen_t>(nj) * g.g + 2 * g.g));
    ++np;
    double *cn = REAL(cw), *cr = cn + static_cast<R_xlen_t>(nj) * g.g,
           *sl = cr + static_cast<R_xlen_t>(nj) * g.g, *sc = sl + g.g;
    SEXP sjv = PROTECT(Rf_allocVector(INTSXP, g.g));
    ++np;
    int npair = 0;
    clab_jsub(INTEGER(j.code), INTEGER(g.code), REAL(dd), n, nj, g.g, mn, rate, cn, cr, nullptr,
              nullptr, nullptr, nullptr, sl, INTEGER(sjv), sc, &npair);
    SEXP iw = PROTECT(Rf_allocVector(INTSXP, 3 * static_cast<R_xlen_t>(npair)));
    SEXP zw = PROTECT(Rf_allocVector(REALSXP, 2 * static_cast<R_xlen_t>(npair)));
    np += 2;
    int *pg = INTEGER(iw), *ps = pg + npair, *pl = ps + npair;
    double *pz = REAL(zw), *pa = pz + npair;
    clab_jsub(INTEGER(j.code), INTEGER(g.code), REAL(dd), n, nj, g.g, mn, rate, cn, cr, pg, ps, pl,
              pz, sl, INTEGER(sjv), sc, &npair);
    bj::cols c(14, npair);
    ++np;
    SEXP cg = PROTECT(Rf_allocVector(STRSXP, npair)), cs = PROTECT(Rf_allocVector(STRSXP, npair)),
         cl = PROTECT(Rf_allocVector(STRSXP, npair));
    np += 3;
    double *rs = nullptr, *rl, *gs, *gl, *ns, *nl, *zz, *pv, *pj;
    int *rv, *sg;
    c.put("group", cg);
    c.put("stricter", cs);
    c.put("lenient", cl);
    rs = c.dbl("rates");
    rl = c.dbl("ratel");
    gs = c.dbl("grates");
    gl = c.dbl("gratel");
    ns = c.dbl("ns");
    nl = c.dbl("nl");
    zz = c.dbl("z");
    pv = c.dbl("p");
    pj = c.dbl("padj");
    rv = c.lgl("reversed");
    sg = c.lgl("signif");
    SEXP ow = PROTECT(Rf_allocVector(INTSXP, npair));
    ++np;
    for (int k = 0; k < npair; ++k) {
      const R_xlen_t ch = ps[k] + static_cast<R_xlen_t>(nj) * pg[k];
      const R_xlen_t cl2 = pl[k] + static_cast<R_xlen_t>(nj) * pg[k];
      SET_STRING_ELT(cg, k, STRING_ELT(g.lev, pg[k]));
      SET_STRING_ELT(cs, k, STRING_ELT(j.lev, ps[k]));
      SET_STRING_ELT(cl, k, STRING_ELT(j.lev, pl[k]));
      rs[k] = rate[ps[k]];
      rl[k] = rate[pl[k]];
      gs[k] = cr[ch];
      gl[k] = cr[cl2];
      ns[k] = cn[ch];
      nl[k] = cn[cl2];
      zz[k] = pz[k];
      pv[k] = Rf_pnorm5(-pz[k], 0.0, 1.0, 1, 0);
      rv[k] = gs[k] < gl[k];
    }
    clab_pbadj(pv, npair, ad, pa, INTEGER(ow));
    int nrev = 0, nsig = 0;
    for (int k = 0; k < npair; ++k) {
      pj[k] = pa[k];
      sg[k] = pa[k] < al;
      nrev += rv[k];
      nsig += sg[k];
    }
    // the pairs and judges in a significant reversal
    std::vector<long long> sp;
    std::vector<int> sj;
    for (int k = 0; k < npair; ++k) {
      if (!sg[k]) continue;
      sp.push_back(static_cast<long long>(ps[k]) * nj + pl[k]);
      sj.push_back(ps[k]);
      sj.push_back(pl[k]);
    }
    std::sort(sp.begin(), sp.end());
    std::sort(sj.begin(), sj.end());
    const double sm[6] = {static_cast<double>(npair),
                          static_cast<double>(nrev),
                          npair ? static_cast<double>(nrev) / npair : NA_REAL,
                          static_cast<double>(nsig),
                          static_cast<double>(std::unique(sp.begin(), sp.end()) - sp.begin()),
                          static_cast<double>(std::unique(sj.begin(), sj.end()) - sj.begin())};
    static const char *const smn[] = {"pairs", "reversed", "share", "signif", "refuted", "judges"};
    // the slopes of the groups with two judges or more
    int ns2 = 0;
    for (int q = 0; q < g.g; ++q) ns2 += INTEGER(sjv)[q] >= 2;
    bj::cols b(4, ns2);
    ++np;
    SEXP bg = PROTECT(Rf_allocVector(STRSXP, ns2));
    ++np;
    b.put("group", bg);
    double *bj2 = b.dbl("judges"), *bc = b.dbl("cases"), *bs = b.dbl("slope");
    for (int q = 0, k = 0; q < g.g; ++q) {
      if (INTEGER(sjv)[q] < 2) continue;
      SET_STRING_ELT(bg, k, STRING_ELT(g.lev, q));
      bj2[k] = INTEGER(sjv)[q];
      bc[k] = sc[q];
      bs[k] = bj::na(sl[q]);
      ++k;
    }
    SEXP sub = PROTECT(Rf_allocVector(VECSXP, 3));
    ++np;
    SET_VECTOR_ELT(sub, 0, c.done());
    SET_VECTOR_ELT(sub, 1, clab::nvec(6, smn, sm));
    SET_VECTOR_ELT(sub, 2, b.done());
    static const char *const subn[] = {"pairs", "summary", "slopes"};
    SET_VECTOR_ELT(out, 3, bj::named(sub, subn));
  }
  SET_VECTOR_ELT(out, 4, Rf_ScalarReal(al));
  static const char *const nm[] = {"judges", "homog", "slope", "sub", "alpha"};
  UNPROTECT(np);
  return bj::named(out, nm);
}

SEXP C_clab_r_mmshr(SEXP sent, SEXP mn, SEXP grp) {
  clab::who = "mmshr";
  int np = 0;
  SEXP sv = PROTECT(bj::nums(sent, "sent"));
  ++np;
  const int n = bj::len(sv, "sent");
  SEXP mv = PROTECT(bj::nums(mn, "min"));
  ++np;
  if (XLENGTH(mv) != 1 && XLENGTH(mv) != n)
    fail("`min` must be one number or one per case, %d.", n);
  SEXP mr = PROTECT(Rf_allocVector(REALSXP, n));
  ++np;
  for (int i = 0; i < n; ++i) REAL(mr)[i] = REAL(mv)[XLENGTH(mv) == 1 ? 0 : i];
  int ng = 0;
  const int *gc = nullptr;
  SEXP lev = R_NilValue;
  if (!Rf_isNull(grp)) {
    csp::lab g = bj::labs(grp, "grp", n);
    np += 2;
    ng = g.g;
    gc = INTEGER(g.code);
    lev = g.lev;
  }
  SEXP ov = PROTECT(Rf_allocVector(REALSXP, 5 * (static_cast<R_xlen_t>(ng) + 1)));
  ++np;
  clab_mmshr(REAL(sv), REAL(mr), gc, n, ng, REAL(ov));
  // the groups, then the overall row; without groups the overall row only
  const int r0 = gc ? 0 : ng, nr = gc ? ng + 1 : 1;
  bj::cols f(6, nr);
  ++np;
  SEXP gl = PROTECT(Rf_allocVector(STRSXP, nr));
  ++np;
  for (int k = 0; k < nr; ++k) {
    SET_STRING_ELT(gl, k, k < nr - 1 ? STRING_ELT(lev, k) : Rf_mkChar("all"));
  }
  f.put("group", gl);
  static const char *const cn[] = {"n", "at", "below", "above", "mexcess"};
  for (int c = 0; c < 5; ++c) {
    double *d = f.dbl(cn[c]);
    for (int k = 0; k < nr; ++k) d[k] = REAL(ov)[c * (ng + 1) + r0 + k];
  }
  SEXP out = PROTECT(f.done());
  ++np;
  UNPROTECT(np);
  return out;
}

}  // extern "C"
