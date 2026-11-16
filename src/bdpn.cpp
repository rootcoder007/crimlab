// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the dynamic and heterogeneous panels: difference
// and system GMM (gmmpnl), Wooldridge's test of serial correlation
// (wsctst) and the common correlated effects estimators (ccepnl). Each
// reads and checks its arguments, hands them to the shared work of
// dpnd.cpp (which clab_jrun() calls too) and builds the R result. Nothing
// calls back into the R evaluator.
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "crimlab.h"
#include "bmat.h"
#include "dpnd.h"
#include "pevr.h"
#include "rchk.h"
#include "spin.h"

using clab::fail;
using pevr::cols;
using pevr::nlist;
using pevr::z;

namespace {

// The panel of dpnd from pevr's.
dpnd::pnl plain(const pevr::pan &pn) {
  dpnd::pnl p;
  p.n = pn.n;
  p.y = pn.y;
  p.u = pn.u;
  p.t = pn.t;
  p.nu = pn.nu;
  p.nt = pn.nt;
  p.tv = pn.tv;
  for (int a = 0; a < pn.nu; ++a) p.ulev.emplace_back(Rf_translateCharUTF8(STRING_ELT(pn.ulev, a)));
  return p;
}

// An interrupted kernel unwinds to the entry point (glue.cpp); a refusal
// is raised in the function's name.
void done(const std::string &msg, int code) {
  if (code == CLAB_EINTR) throw clab::Intr{};
  if (!msg.empty()) fail("%s", msg.c_str());
}

// Columns of x named by s (character), each once.
std::vector<int> cnames(SEXP s, const char *an, const bmat::mat &xm) {
  if (TYPEOF(s) != STRSXP) fail("`%s` must be names of columns of `x`.", an);
  std::vector<int> o;
  for (R_xlen_t i = 0; i < XLENGTH(s); ++i) {
    if (STRING_ELT(s, i) == NA_STRING) fail("`%s` must not hold NA.", an);
    const std::string v = Rf_translateCharUTF8(STRING_ELT(s, i));
    const auto it = std::find(xm.nm.begin(), xm.nm.end(), v);
    if (it == xm.nm.end()) fail("`%s` names %s, which is not a column of `x`.", an, v.c_str());
    const int j = static_cast<int>(it - xm.nm.begin());
    if (std::find(o.begin(), o.end(), j) != o.end()) fail("`%s` names %s twice.", an, v.c_str());
    o.push_back(j);
  }
  return o;
}

// Lag limits, one value or one per GMM-style variable: whole numbers from
// lo, or Inf (big) when inf.
std::vector<int> lags(SEXP s, const char *an, int nv, int lo, bool inf) {
  if (!clab::isnum(s) || (XLENGTH(s) != 1 && XLENGTH(s) != nv)) {
    fail("`%s` must be one number or %d, one per GMM-style variable.", an, nv);
  }
  SEXP d = PROTECT(clab::dbl(s));
  std::vector<int> o;
  for (int i = 0; i < nv; ++i) {
    const double v = REAL(d)[XLENGTH(s) == 1 ? 0 : i];
    if (inf && v == R_PosInf) {
      o.push_back(1 << 30);
      continue;
    }
    if (!std::isfinite(v) || v != std::floor(v) || v < lo || v > 1e6) {
      fail("`%s` must hold whole numbers from %d%s.", an, lo, inf ? ", or Inf" : "");
    }
    o.push_back(static_cast<int>(v));
  }
  UNPROTECT(1);
  return o;
}

// The options of gmmpnl() into sp, in the order clab_jrun() reads them.
void gopts(dpnd::spec &sp, const bmat::mat &xm, SEXP ylag, SEXP gmm, SEXP iv, SEXP lmin,
           SEXP lmax, SEXP coll, SEXP type, SEXP trans, SEXP hmat, SEXP cons, SEXP tdum) {
  static const char *const tyo[] = {"system", "diff"};
  static const char *const tro[] = {"fd", "fod"};
  static const char *const hmo[] = {"ab", "dpd", "id"};
  sp.ylag = clab::whole(ylag, "ylag", 0);
  if (sp.ylag > 1000) fail("`ylag` must be at most 1000.");
  if (sp.ylag == 0 && xm.p == 0) fail("the model needs `ylag` above 0 or columns in `x`.");
  sp.sys = csp::pick(type, "type", tyo, 2) == 0;
  sp.fod = csp::pick(trans, "trans", tro, 2);
  sp.hm = csp::pick(hmat, "hmat", hmo, 3);
  sp.coll = csp::flag(coll, "collapse");
  sp.cons = csp::flag(cons, "cons");
  sp.tdum = csp::flag(tdum, "tdum");
  const std::vector<int> gx = Rf_isNull(gmm) ? std::vector<int>() : cnames(gmm, "gmm", xm);
  if (sp.ylag > 0) sp.gv.push_back(-1);
  sp.gv.insert(sp.gv.end(), gx.begin(), gx.end());
  if (Rf_isNull(iv)) {
    for (int j = 0; j < xm.p; ++j) {
      if (std::find(gx.begin(), gx.end(), j) == gx.end()) sp.iv.push_back(j);
    }
  } else {
    sp.iv = cnames(iv, "iv", xm);
    for (int j : sp.iv) {
      if (std::find(gx.begin(), gx.end(), j) != gx.end()) {
        fail("`iv` and `gmm` both name %s.", xm.nm[z(j)].c_str());
      }
    }
  }
  const int nv = static_cast<int>(sp.gv.size());
  if (nv > 0) {
    sp.lo = lags(lmin, "lagmin", nv, 1, false);
    sp.hi = lags(lmax, "lagmax", nv, 1, true);
    for (int i = 0; i < nv; ++i) {
      if (sp.hi[z(i)] < sp.lo[z(i)]) fail("`lagmax` must be at least `lagmin`.");
    }
  }
}

SEXP strs(const std::vector<std::string> &v) {
  SEXP s = PROTECT(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(v.size())));
  for (size_t i = 0; i < v.size(); ++i) {
    SET_STRING_ELT(s, static_cast<R_xlen_t>(i), Rf_mkCharCE(v[i].c_str(), CE_UTF8));
  }
  UNPROTECT(1);
  return s;
}

SEXP rmat(const std::vector<double> &v, int r, int c) {
  SEXP m = PROTECT(Rf_allocMatrix(REALSXP, r, c));
  std::copy(v.begin(), v.end(), REAL(m));
  UNPROTECT(1);
  return m;
}

double *dcol(cols &c, const char *nm, const std::vector<double> &v) {
  const int n = static_cast<int>(v.size());
  double *p = c.dbl(nm, n);
  for (int i = 0; i < n; ++i) p[i] = pevr::na(v[z(i)]);
  return p;
}
SEXP scol(cols &c, const char *nm, const std::vector<std::string> &v) {
  const int n = static_cast<int>(v.size());
  SEXP s = c.str(nm, n);
  for (int i = 0; i < n; ++i) SET_STRING_ELT(s, i, Rf_mkCharCE(v[z(i)].c_str(), CE_UTF8));
  return s;
}

// A coefficient table: term, est, se, z (or t), p, lo, hi.
SEXP ctab(const dpnd::tab &t, bool tst) {
  cols c(7);
  scol(c, "term", t.name);
  dcol(c, "est", t.est);
  dcol(c, "se", t.se);
  dcol(c, tst ? "t" : "z", t.stat);
  dcol(c, "p", t.p);
  dcol(c, "lo", t.lo);
  dcol(c, "hi", t.hi);
  return c.frame();
}

}  // namespace

extern "C" {

//' @srrstats {G2.6} gmmpnl(), wsctst() and ccepnl() read the outcome, units
//'   and periods as plain vectors of one length; regressors, instruments
//'   and observed factors as a numeric matrix with one row per row.
//' @noRd
SEXP C_clab_r_gmmpnl(SEXP y, SEXP x, SEXP unit, SEXP time, SEXP ylag, SEXP gmm, SEXP iv,
                     SEXP zz, SEXP lmin, SEXP lmax, SEXP coll, SEXP type, SEXP trans, SEXP steps,
                     SEXP robust, SEXP wc, SEXP small, SEXP hmat, SEXP cons, SEXP tdum, SEXP ar,
                     SEXP level) {
  clab::who = "gmmpnl";
  int np = 0;
  pevr::pan pn = pevr::panel(y, unit, time, clab::ninf, &np);
  const bmat::mat xm = Rf_isNull(x) ? bmat::mat() : bmat::readm(x, "x", "x", pn.n, clab::ninf);
  const bmat::mat zx = Rf_isNull(zz) ? bmat::mat() : bmat::readm(zz, "z", "z", pn.n, clab::ninf);
  dpnd::spec sp;
  gopts(sp, xm, ylag, gmm, iv, lmin, lmax, coll, type, trans, hmat, cons, tdum);
  sp.steps = clab::whole(steps, "steps", 1);
  if (sp.steps > 2) fail("`steps` must be 1 or 2.");
  sp.rob = csp::flag(robust, "robust");
  sp.wc = Rf_isNull(wc) ? -1 : csp::flag(wc, "wc");
  done(dpnd::gcheck(sp), CLAB_OK);
  sp.small = csp::flag(small, "small");
  sp.nar = clab::whole(ar, "ar", 0);
  if (sp.nar > 100) fail("`ar` must be at most 100.");
  sp.level = clab::one(level, "level", 0.0, 1.0);
  dpnd::gout o;
  int code = CLAB_OK;
  const std::string msg = dpnd::gfit(plain(pn), xm, zx, sp, o, code);
  done(msg, code);
  const int k = static_cast<int>(o.k);
  SEXP coef = PROTECT(ctab(o.coef, o.tstat));
  cols tc(4);
  scol(tc, "test", o.tests.name);
  dcol(tc, "stat", o.tests.stat);
  dcol(tc, "df", o.tests.df);
  dcol(tc, "p", o.tests.p);
  SEXP te = PROTECT(tc.frame());
  cols wc5(5);
  scol(wc5, "test", o.wald.name);
  dcol(wc5, "stat", o.wald.stat);
  dcol(wc5, "df1", o.wald.df);
  dcol(wc5, "df2", o.wald.df2);
  dcol(wc5, "p", o.wald.p);
  SEXP wd = PROTECT(wc5.frame());
  cols xw5(5);
  scol(xw5, "test", o.xwald.name);
  dcol(xw5, "stat", o.xwald.stat);
  dcol(xw5, "df1", o.xwald.df);
  dcol(xw5, "df2", o.xwald.df2);
  dcol(xw5, "p", o.xwald.p);
  SEXP xwd = PROTECT(xw5.frame());
  SEXP rs = PROTECT(Rf_allocVector(REALSXP, pn.n));
  for (int i = 0; i < pn.n; ++i) REAL(rs)[i] = pevr::na(o.resid[z(i)]);
  cols fs(9);
  dcol(fs, "n", {o.n});
  dcol(fs, "groups", {o.groups});
  dcol(fs, "inst", {o.inst});
  dcol(fs, "omin", {o.omin});
  dcol(fs, "oavg", {o.oavg});
  dcol(fs, "omax", {o.omax});
  dcol(fs, "steps", {o.steps});
  dcol(fs, "k", {o.k});
  dcol(fs, "sigma2", {o.sigma2});
  SEXP fit = PROTECT(fs.frame());
  SEXP vm = PROTECT(rmat(o.vcov, k, k));
  SEXP dn = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dn, 0, strs(o.coef.name));
  SET_VECTOR_ELT(dn, 1, strs(o.coef.name));
  Rf_setAttrib(vm, R_DimNamesSymbol, dn);
  SEXP in = PROTECT(strs(o.inst_names)), dr = PROTECT(strs(o.dropped));
  SEXP du = PROTECT(strs(o.dummies)), nt = PROTECT(strs(o.note));
  SEXP ty = PROTECT(Rf_mkString(o.type.c_str())), tr = PROTECT(Rf_mkString(o.trans.c_str()));
  SEXP lvs = PROTECT(Rf_ScalarReal(sp.level));
  SEXP rb = PROTECT(Rf_ScalarLogical(sp.rob)), sm = PROTECT(Rf_ScalarLogical(sp.small));
  np += 17;
  SEXP out = PROTECT(nlist({{"coef", coef}, {"tests", te}, {"wald", wd}, {"xwald", xwd},
                            {"fit", fit}, {"vcov", vm}, {"resid", rs}, {"inst", in},
                            {"dropped", dr}, {"dummies", du}, {"note", nt}, {"type", ty},
                            {"trans", tr}, {"robust", rb}, {"small", sm}, {"level", lvs}}));
  Rf_classgets(out, PROTECT(Rf_mkString("gmmpnl")));
  UNPROTECT(np + 2);
  return out;
}

SEXP C_clab_r_wsctst(SEXP y, SEXP x, SEXP unit, SEXP time, SEXP cons, SEXP ssc) {
  clab::who = "wsctst";
  static const char *const sso[] = {"cr1", "none"};
  int np = 0;
  pevr::pan pn = pevr::panel(y, unit, time, clab::ninf, &np);
  const bmat::mat xm = Rf_isNull(x) ? bmat::mat() : bmat::readm(x, "x", "x", pn.n, clab::ninf);
  const int c = csp::flag(cons, "cons");
  const int raw = csp::pick(ssc, "ssc", sso, 2);
  dpnd::wout o;
  int code = CLAB_OK;
  const std::string msg = dpnd::wfit(plain(pn), xm, c, raw, o, code);
  done(msg, code);
  static const char *const sn[] = {"F", "df1", "df2", "p", "coef", "se", "pairs", "clusters",
                                   "ndiff"};
  SEXP sv = PROTECT(clab::nvec(9, sn, o.stat));
  SEXP co = R_NilValue;
  if (xm.p > 0) {
    cols cf(3);
    scol(cf, "term", xm.nm);
    dcol(cf, "est", o.b);
    dcol(cf, "se", o.se);
    co = cf.frame();
  }
  PROTECT(co);
  SEXP sv2 = PROTECT(Rf_mkString(sso[raw]));
  np += 3;
  SEXP out = PROTECT(nlist({{"stat", sv}, {"coef", co}, {"ssc", sv2}}));
  Rf_classgets(out, PROTECT(Rf_mkString("wsctst")));
  UNPROTECT(np + 2);
  return out;
}

SEXP C_clab_r_ccepnl(SEXP y, SEXP x, SEXP unit, SEXP time, SEXP fac, SEXP cons, SEXP w,
                     SEXP avg, SEXP ylag, SEXP lags, SEXP bias, SEXP pool, SEXP pvar,
                     SEXP level) {
  clab::who = "ccepnl";
  static const char *const bso[] = {"none", "jack", "rma"};
  static const char *const pvo[] = {"pesaran", "bsw"};
  int np = 0;
  pevr::pan pn = pevr::panel(y, unit, time, clab::ninf, &np);
  const int n = pn.n;
  const bmat::mat xm = bmat::readm(x, "x", "x", n, clab::ninf);
  const bmat::mat fm = Rf_isNull(fac) ? bmat::mat() : bmat::readm(fac, "fac", "f", n, clab::ninf);
  const int cn = csp::flag(cons, "cons");
  std::vector<std::string> an;
  if (!Rf_isNull(avg)) {
    if (TYPEOF(avg) != STRSXP) fail("`avg` must be names: y or columns of `x`.");
    for (R_xlen_t i = 0; i < XLENGTH(avg); ++i) {
      if (STRING_ELT(avg, i) == NA_STRING) fail("`avg` must not hold NA.");
      an.emplace_back(Rf_translateCharUTF8(STRING_ELT(avg, i)));
    }
  }
  dpnd::copt co;
  co.ylag = clab::whole(ylag, "ylag", 0);
  if (co.ylag > 100) fail("`ylag` must be at most 100.");
  if (!Rf_isNull(lags)) {
    co.lags = clab::whole(lags, "lags", 0);
    if (co.lags > 1000) fail("`lags` must be at most 1000.");
  }
  co.bias = csp::pick(bias, "bias", bso, 3);
  co.pool = Rf_isNull(pool) ? -1 : csp::flag(pool, "pool");
  co.pvar = csp::pick(pvar, "pvar", pvo, 2);
  const double lv = clab::one(level, "level", 0.0, 1.0);
  std::vector<double> wv;
  if (!Rf_isNull(w)) {
    clab::vec(w, "w");
    SEXP ws = PROTECT(clab::num(w, "w", n, 0.0, true));
    ++np;
    wv.assign(REAL(ws), REAL(ws) + n);
  }
  dpnd::cres o;
  int code = CLAB_OK;
  const std::string msg =
      dpnd::cfit(plain(pn), xm, fm, cn, wv, Rf_isNull(avg), an, co, lv, o, code);
  done(msg, code);
  const int k = static_cast<int>(o.k), nu = pn.nu;
  SEXP mg = PROTECT(ctab(o.mg, false));
  SEXP pl = PROTECT(o.pool ? ctab(o.pooled, false) : R_NilValue);
  SEXP m0 = PROTECT(o.bias ? ctab(o.mg0, false) : R_NilValue);
  SEXP bm = PROTECT(rmat(o.bi, nu, k));
  SEXP dn = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dn, 0, pn.ulev);
  SET_VECTOR_ELT(dn, 1, strs(o.nm));
  Rf_setAttrib(bm, R_DimNamesSymbol, dn);
  SEXP vm1 = PROTECT(rmat(o.vmg, k, k));
  SEXP vm2 = PROTECT(o.pool ? rmat(o.vp, k, k) : R_NilValue);
  SEXP vm0 = PROTECT(o.bias ? rmat(o.vmg0, k, k) : R_NilValue);
  static const char *const hn[] = {"stat", "df", "p"};
  static const char *const fn[] = {"er", "thr", "avg"};
  double hv[3], fv[3];
  for (int j = 0; j < 3; ++j) {
    hv[j] = pevr::na(o.hb[j]);
    fv[j] = pevr::na(o.nf[j]);
  }
  SEXP hb = PROTECT(o.pool ? clab::nvec(3, hn, hv) : R_NilValue);
  SEXP nf = PROTECT(clab::nvec(3, fn, fv));
  cols fs(7);
  dcol(fs, "units", {o.units});
  dcol(fs, "periods", {o.periods});
  dcol(fs, "k", {o.k});
  dcol(fs, "h", {o.h});
  dcol(fs, "n", {o.n});
  dcol(fs, "lags", {o.lags});
  dcol(fs, "ylag", {o.ylag});
  SEXP fit = PROTECT(fs.frame());
  SEXP wr = PROTECT(Rf_allocVector(REALSXP, nu));
  std::copy(o.w.begin(), o.w.end(), REAL(wr));
  Rf_setAttrib(wr, R_NamesSymbol, pn.ulev);
  SEXP lvs = PROTECT(Rf_ScalarReal(lv));
  SEXP bs = PROTECT(Rf_mkString(bso[co.bias])), pvs = PROTECT(Rf_mkString(pvo[co.pvar]));
  np += 15;
  SEXP out = PROTECT(nlist({{"mg", mg}, {"pooled", pl}, {"mg0", m0}, {"bi", bm}, {"vmg", vm1},
                            {"vp", vm2}, {"vmg0", vm0}, {"hb", hb}, {"nfac", nf}, {"fit", fit},
                            {"w", wr}, {"bias", bs}, {"pvar", pvs}, {"level", lvs}}));
  Rf_classgets(out, PROTECT(Rf_mkString("ccepnl")));
  UNPROTECT(np + 2);
  return out;
}

}  // extern "C"
