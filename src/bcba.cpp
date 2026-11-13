// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the design-based causal estimators (cbatt), the
// balance and overlap table (smdbal) and the E-value (evlx): check the
// arguments, call cba.cpp and build the R result. Nothing calls back into
// the R evaluator.
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "kutil.h"
#include "pevc.h"
#include "pevr.h"
#include "rchk.h"
#include "spin.h"

extern "C" double Rf_qnorm5(double, double, double, int, int);

using clab::fail;
using pevr::cols;
using pevr::z;

namespace {

// The covariates: a numeric matrix, a data frame or list of numeric
// columns, or one numeric vector; n rows, finite.
void xmat(SEXP x, int n, std::vector<double> &v, std::vector<std::string> &nm) {
  std::vector<SEXP> cs;
  SEXP names = R_NilValue;
  if (TYPEOF(x) == VECSXP) {
    for (R_xlen_t j = 0; j < XLENGTH(x); ++j) cs.push_back(VECTOR_ELT(x, j));
    names = Rf_getAttrib(x, R_NamesSymbol);
  } else if (Rf_isMatrix(x)) {
    if (!clab::isnum(x) || Rf_nrows(x) != n) fail("`x` must be a numeric matrix with one row per unit, %d.", n);
    SEXP dn = Rf_getAttrib(x, R_DimNamesSymbol);
    SEXP cn = Rf_isNull(dn) ? R_NilValue : VECTOR_ELT(dn, 1);
    SEXP xd = PROTECT(clab::dbl(x));
    v.assign(REAL(xd), REAL(xd) + XLENGTH(xd));
    UNPROTECT(1);
    for (double e : v) {
      if (!std::isfinite(e)) fail("`x` must hold finite numbers with no missing values.");
    }
    nm.clear();
    for (int j = 0; j < Rf_ncols(x); ++j) {
      const bool has = !Rf_isNull(cn) && STRING_ELT(cn, j) != NA_STRING && CHAR(STRING_ELT(cn, j))[0] != '\0';
      nm.push_back(has ? Rf_translateCharUTF8(STRING_ELT(cn, j)) : "x" + std::to_string(j + 1));
    }
    if (nm.empty()) fail("`x` must hold one covariate at least.");
    return;
  } else {
    cs.push_back(x);
  }
  if (cs.empty()) fail("`x` must hold one covariate at least.");
  v.clear();
  nm.clear();
  for (std::size_t j = 0; j < cs.size(); ++j) {
    SEXP c = PROTECT(clab::num(cs[j], "x", n));
    v.insert(v.end(), REAL(c), REAL(c) + n);
    UNPROTECT(1);
    const bool has = !Rf_isNull(names) && STRING_ELT(names, static_cast<R_xlen_t>(j)) != NA_STRING &&
                     CHAR(STRING_ELT(names, static_cast<R_xlen_t>(j)))[0] != '\0';
    nm.push_back(has ? Rf_translateCharUTF8(STRING_ELT(names, static_cast<R_xlen_t>(j)))
                     : "x" + std::to_string(j + 1));
  }
}

// The treatment as 0.0 and 1.0, with units in both arms.
std::vector<double> arms(SEXP d, int n) {
  SEXP dv = PROTECT(csp::bin01(d, "d", n));
  std::vector<double> o(INTEGER(dv), INTEGER(dv) + n);
  UNPROTECT(1);
  const double n1 = std::count(o.begin(), o.end(), 1.0);
  if (n1 < 2 || n - n1 < 2) fail("`d` must have two treated and two untreated units at least.");
  return o;
}

// Given propensities: n numbers strictly between 0 and 1.
std::vector<double> given(SEXP ps, int n) {
  SEXP p = PROTECT(clab::num(ps, "ps", n));
  std::vector<double> e(REAL(p), REAL(p) + n);
  UNPROTECT(1);
  for (double v : e) {
    if (!(v > 0 && v < 1)) fail("`ps` must be propensities strictly between 0 and 1.");
  }
  return e;
}

const char *const sepmsg =
    "the logistic propensity has no finite fit (a covariate separates the arms); give `ps`.";

}  // namespace

extern "C" {

//' @srrstats {G2.6} cbatt() and smdbal() read the outcome, treatment,
//'   weights and propensities as plain vectors, the covariates as a matrix,
//'   a data frame or one vector.
//' @noRd
SEXP C_clab_r_cbatt(SEXP y, SEXP d, SEXP x, SEXP est, SEXP tgt, SEXP ps, SEXP k, SEXP calip,
                    SEXP nsub, SEXP reps, SEXP seed, SEXP level) {
  clab::who = "cbatt";
  static const char *const eo[] = {"ipw", "gcomp", "sub", "psm", "ebal"};
  static const char *const to[] = {"ate", "att", "atc"};
  const int me = csp::pick(est, "est", eo, 5), tg = csp::pick(tgt, "tgt", to, 3);
  const int kk = clab::whole(k, "k", 1), ns = clab::whole(nsub, "nsub", 2);
  const int nr = clab::whole(reps, "reps", 0);
  const double lv = clab::one(level, "level", 0.0, 1.0);
  clab::vec(y, "y");
  SEXP yv = PROTECT(clab::num(y, "y"));
  const int n = clab::ilen(static_cast<double>(XLENGTH(yv)), "y");
  std::vector<double> dv = arms(d, n), xv;
  std::vector<std::string> nm;
  xmat(x, n, xv, nm);
  if (me == 4 && tg != 1) fail("entropy balancing (`est = \"ebal\"`) estimates the ATT only.");
  if (me == 3 && tg == 0) fail("matching (`est = \"psm\"`) estimates the ATT or the ATC.");
  pev::cdat c{REAL(yv), dv.data(), xv.data(), n, static_cast<int>(nm.size())};
  std::vector<int> all(z(n));
  for (int i = 0; i < n; ++i) all[z(i)] = i;
  // the propensity, for the estimators that use it
  const bool need = me == 0 || me == 2 || me == 3, own = need && Rf_isNull(ps);
  std::vector<double> e(z(n), kutil::nan());
  if (own) {
    if (!pev::pscore(c, all, e)) fail("%s", sepmsg);
  } else if (need) {
    e = given(ps, n);
  }
  std::vector<double> w;
  // one estimate on the units idx (a bootstrap sample has repeats)
  auto one = [&](const std::vector<int> &idx, std::vector<double> &ee) {
    if (own && &ee != &e && !pev::pscore(c, idx, ee)) return kutil::nan();
    switch (me) {
      case 0: return pev::ipw(c, idx, ee, tg);
      case 1: return pev::gcomp(c, idx, tg);
      case 2: return pev::subcl(c, idx, ee, tg, ns);
      default: {
        std::vector<double> ww;
        return pev::ebal(c, idx, ww) ? pev::ebatt(c, idx, ww) : kutil::nan();
      }
    }
  };
  double b = NA_REAL, se = NA_REAL, med = NA_REAL;
  SEXP mt = R_NilValue;
  int np = 1;
  if (me == 3) {
    double cal = 0.0;
    if (!Rf_isNull(calip)) {
      double m = 0.0, s = 0.0;
      for (double v : e) m += v / n;
      for (double v : e) s += (v - m) * (v - m) / (n - 1);
      cal = clab::one(calip, "calip", 0.0, R_PosInf) * std::sqrt(s);
    }
    std::vector<int> ft, fc;
    std::vector<double> dif;
    const int m = pev::psm(c, e, tg, kk, cal, ft, fc, dif);
    if (m < 2) fail("fewer than two units were matched within the caliper.");
    double s = 0.0, ss = 0.0;
    for (double v : dif) s += v / m;
    for (double v : dif) ss += (v - s) * (v - s) / (m - 1);
    b = s;
    se = std::sqrt(ss / m);
    std::vector<double> sd = dif;
    std::sort(sd.begin(), sd.end());
    med = m % 2 ? sd[z(m / 2)] : 0.5 * (sd[z(m / 2 - 1)] + sd[z(m / 2)]);
    const int np2 = static_cast<int>(ft.size());
    cols cm(2);
    double *a = cm.dbl("focal", np2), *bb = cm.dbl("match", np2);
    for (int q = 0; q < np2; ++q) {
      a[q] = ft[z(q)] + 1;
      bb[q] = fc[z(q)] + 1;
    }
    mt = PROTECT(cm.frame());
    ++np;
  } else {
    b = one(all, e);
    if (std::isnan(b)) {
      fail("the estimate is undefined: an outcome model failed, a subclass lacks an arm, or "
           "the balancing weights do not exist.");
    }
    if (me == 0 || me == 4) {
      w.assign(z(n), 1.0);
      if (me == 4) {
        pev::ebal(c, all, w);
      } else {
        for (int i = 0; i < n; ++i) {
          const double p = e[z(i)];
          w[z(i)] = dv[z(i)] > 0.5 ? (tg == 0 ? 1 / p : tg == 1 ? 1 : (1 - p) / p)
                                   : (tg == 0 ? 1 / (1 - p) : tg == 1 ? p / (1 - p) : 1);
        }
      }
    }
    if (nr > 0) {
      kutil::rng g(csp::seed(seed));
      kutil::ksum s1, s2;
      double kq = 0;
      std::vector<int> idx(z(n));
      std::vector<double> eb = e;
      for (int r = 0; r < nr; ++r) {
        for (int i = 0; i < n; ++i) idx[z(i)] = g.below(n);
        const double v = one(idx, eb);
        if (std::isnan(v)) continue;
        s1.add(v);
        s2.add(v * v);
        kq += 1;
      }
      if (kq > 1) {
        const double m = s1.val() / kq;
        se = std::sqrt(std::max(0.0, (s2.val() - kq * m * m) / (kq - 1)));
      }
    }
  }
  const double zq = Rf_qnorm5(0.5 + 0.5 * lv, 0.0, 1.0, 1, 0);
  double n1 = 0;
  for (double v : dv) n1 += v;
  cols r(7);
  double *re = r.dbl("est", 1), *rs = r.dbl("se", 1), *rl = r.dbl("lo", 1), *rh = r.dbl("hi", 1),
         *rm = r.dbl("med", 1), *r1 = r.dbl("n1", 1), *r0 = r.dbl("n0", 1);
  re[0] = b;
  rs[0] = se;
  rl[0] = ISNAN(se) ? NA_REAL : b - zq * se;
  rh[0] = ISNAN(se) ? NA_REAL : b + zq * se;
  rm[0] = med;
  r1[0] = n1;
  r0[0] = n - n1;
  SEXP res = PROTECT(r.frame());
  SEXP pv = PROTECT(Rf_allocVector(REALSXP, n));
  std::copy(e.begin(), e.end(), REAL(pv));
  SEXP wv = R_NilValue;
  if (!w.empty()) {
    wv = PROTECT(Rf_allocVector(REALSXP, n));
    std::copy(w.begin(), w.end(), REAL(wv));
    ++np;
  }
  SEXP e1 = PROTECT(Rf_mkString(eo[me])), t1 = PROTECT(Rf_mkString(to[tg])),
       l1 = PROTECT(Rf_ScalarReal(lv));
  np += 5;
  SEXP out = PROTECT(pevr::nlist({{"res", res}, {"ps", pv}, {"w", wv}, {"match", mt},
                                  {"est", e1}, {"tgt", t1}, {"level", l1}}));
  Rf_classgets(out, PROTECT(Rf_mkString("cbatt")));
  UNPROTECT(np + 2);
  return out;
}

SEXP C_clab_r_smdbal(SEXP x, SEXP d, SEXP w, SEXP ps, SEXP thr) {
  clab::who = "smdbal";
  const double th = clab::one(thr, "thr", 0.0, R_PosInf);
  const int n = clab::ilen(static_cast<double>(Rf_xlength(d)), "d");
  std::vector<double> dv = arms(d, n), xv, wv;
  std::vector<std::string> nm;
  xmat(x, n, xv, nm);
  if (!Rf_isNull(w)) {
    SEXP wx = PROTECT(clab::num(w, "w", n, 0.0));
    wv.assign(REAL(wx), REAL(wx) + n);
    UNPROTECT(1);
  }
  const int p = static_cast<int>(nm.size());
  pev::cdat c{dv.data(), dv.data(), xv.data(), n, p};
  std::vector<double> e;
  if (Rf_isNull(ps)) {
    std::vector<int> all(z(n));
    for (int i = 0; i < n; ++i) all[z(i)] = i;
    if (!pev::pscore(c, all, e)) fail("%s", sepmsg);
  } else {
    e = given(ps, n);
  }
  cols t(wv.empty() ? 5 : 9);
  SEXP tv = t.str("var", p);
  double *m1 = t.dbl("m1", p), *m0 = t.dbl("m0", p), *sm = t.dbl("smd", p), *vr = t.dbl("vr", p);
  double *a1 = nullptr, *a0 = nullptr, *as = nullptr, *av = nullptr;
  if (!wv.empty()) {
    a1 = t.dbl("m1w", p);
    a0 = t.dbl("m0w", p);
    as = t.dbl("smdw", p);
    av = t.dbl("vrw", p);
  }
  bool bal = true;
  for (int j = 0; j < p; ++j) {
    double o[4];
    SET_STRING_ELT(tv, j, Rf_mkCharCE(nm[z(j)].c_str(), CE_UTF8));
    pev::smd(xv.data() + z(j) * z(n), n, dv.data(), nullptr, o);
    m1[j] = o[0];
    m0[j] = o[1];
    sm[j] = o[2];
    vr[j] = o[3];
    if (!wv.empty()) {
      pev::smd(xv.data() + z(j) * z(n), n, dv.data(), wv.data(), o);
      a1[j] = o[0];
      a0[j] = o[1];
      as[j] = o[2];
      av[j] = o[3];
    }
    bal = bal && std::fabs(wv.empty() ? sm[j] : as[j]) < th;
  }
  SEXP tab = PROTECT(t.frame());
  // overlap: the propensity's range in each arm, the share of units outside
  // the common range, the largest inverse probability weight
  double lo[2] = {1, 1}, hi[2] = {0, 0}, wmax = 0;
  for (int i = 0; i < n; ++i) {
    const int g = dv[z(i)] > 0.5 ? 1 : 0;
    lo[g] = std::min(lo[g], e[z(i)]);
    hi[g] = std::max(hi[g], e[z(i)]);
    wmax = std::max(wmax, g ? 1 / e[z(i)] : 1 / (1 - e[z(i)]));
  }
  const double cl = std::max(lo[0], lo[1]), ch = std::min(hi[0], hi[1]);
  double out = 0;
  for (double v : e) out += (v < cl || v > ch) ? 1.0 / n : 0.0;
  cols ov(6);
  double *o1 = ov.dbl("min1", 1), *o2 = ov.dbl("max1", 1), *o3 = ov.dbl("min0", 1),
         *o4 = ov.dbl("max0", 1), *o5 = ov.dbl("out", 1), *o6 = ov.dbl("wmax", 1);
  o1[0] = lo[1];
  o2[0] = hi[1];
  o3[0] = lo[0];
  o4[0] = hi[0];
  o5[0] = out;
  o6[0] = wmax;
  SEXP ovl = PROTECT(ov.frame());
  SEXP b1 = PROTECT(Rf_ScalarLogical(bal)), b2 = PROTECT(Rf_ScalarLogical(out == 0)),
       t1 = PROTECT(Rf_ScalarReal(th));
  SEXP res = PROTECT(pevr::nlist({{"tab", tab}, {"overlap", ovl}, {"bal", b1}, {"ovl", b2},
                                  {"thr", t1}}));
  Rf_classgets(res, PROTECT(Rf_mkString("smdbal")));
  UNPROTECT(7);
  return res;
}

//' @srrstats {G2.6} evlx() reads one estimate and its optional limits.
//' @noRd
SEXP C_clab_r_evlx(SEXP est, SEXP lo, SEXP hi, SEXP type, SEXP rare, SEXP sd) {
  clab::who = "evlx";
  static const char *const to[] = {"rr", "or", "hr", "smd"};
  const int ty = csp::pick(type, "type", to, 4);
  const bool rr = csp::flag(rare, "rare");
  const double lw = ty == 3 ? clab::ninf : 0.0;
  auto get = [&](SEXP v, const char *nm) {
    if (Rf_isNull(v)) return kutil::nan();
    SEXP q = PROTECT(clab::num(v, nm, 1, lw, ty != 3));
    const double r = REAL(q)[0];
    UNPROTECT(1);
    return r;
  };
  const double e = get(est, "est"), l = get(lo, "lo"), h = get(hi, "hi");
  if (std::isnan(e)) fail("`est` must be one estimate.");
  if (!std::isnan(l) && !std::isnan(h) && !(l <= e && e <= h)) fail("`lo` and `hi` must surround `est`.");
  double s = kutil::nan();
  if (!Rf_isNull(sd)) {
    if (ty != 3) fail("`sd` is the standard error of a standardised difference (`type = \"smd\"`).");
    s = clab::one(sd, "sd", 0.0, R_PosInf);
  }
  double o[3];
  pev::evalue(e, l, h, ty, rr, s, o);
  static const char *const on[] = {"rr", "evalue", "evlim"};
  for (double &v : o) v = pevr::na(v);
  return clab::nvec(3, on, o);
}

}  // extern "C"
