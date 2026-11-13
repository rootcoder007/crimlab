// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry point of design-based survey estimates (svymx): checks the
// design and the analysis variable, forms each domain's estimate and its
// influence values, and takes their linearisation variance (svy.cpp).
// Nothing calls back into the R evaluator.
#include <algorithm>
#include <cmath>
#include <vector>

#include "pev.h"
#include "pevr.h"
#include "rchk.h"
#include "spin.h"

extern "C" double Rf_qt(double, double, int, int);
extern "C" double Rf_pt(double, double, int, int);

using clab::fail;
using pevr::cols;
using pevr::z;

extern "C" {

//' @srrstats {G2.6} svymx() reads the analysis variable, weights, strata,
//'   first-stage units, population counts and domains as plain vectors of
//'   one length.
//' @noRd
SEXP C_clab_r_svymx(SEXP y, SEXP w, SEXP strata, SEXP psu, SEXP fpc, SEXP stat, SEXP den,
                    SEXP by, SEXP ref, SEXP level) {
  clab::who = "svymx";
  static const char *const so[] = {"mean", "total", "prop", "ratio", "var"};
  const int sk = csp::pick(stat, "stat", so, 5);
  const double lv = clab::one(level, "level", 0.0, 1.0);
  clab::vec(y, "y");
  SEXP yv = PROTECT(clab::num(y, "y"));
  int np = 1;
  const int n = clab::ilen(static_cast<double>(XLENGTH(yv)), "y");
  if (n < 2) fail("`y` must have two values at least.");
  const double *yy = REAL(yv);
  if (sk == 2) {
    for (int i = 0; i < n; ++i) {
      if (yy[i] != 0 && yy[i] != 1) fail("`stat = \"prop\"` needs `y` of 0 and 1.");
    }
  }
  std::vector<double> wv(z(n), 1.0), xv;
  if (!Rf_isNull(w)) {
    SEXP wx = PROTECT(clab::num(w, "w", n, 0.0, true));
    ++np;
    wv.assign(REAL(wx), REAL(wx) + n);
  }
  if (sk == 3) {
    if (Rf_isNull(den)) fail("`stat = \"ratio\"` needs `den`.");
    SEXP dx = PROTECT(clab::num(den, "den", n));
    ++np;
    xv.assign(REAL(dx), REAL(dx) + n);
  } else if (!Rf_isNull(den)) {
    fail("`den` is the denominator of `stat = \"ratio\"`.");
  }
  // strata and first-stage units (nested: a unit label counts within its
  // stratum)
  std::vector<int> st(z(n), 0), ps(z(n));
  int nst = 1;
  if (!Rf_isNull(strata)) {
    csp::lab l = csp::codes(strata, "strata", n);
    st.assign(INTEGER(l.code), INTEGER(l.code) + n);
    nst = l.g;
    UNPROTECT(2);
  }
  std::vector<int> pu(z(n));
  int npu = n;
  if (Rf_isNull(psu)) {
    for (int i = 0; i < n; ++i) pu[z(i)] = i;
  } else {
    csp::lab l = csp::codes(psu, "psu", n);
    pu.assign(INTEGER(l.code), INTEGER(l.code) + n);
    npu = l.g;
    UNPROTECT(2);
  }
  std::vector<int> iw(z(nst) * z(npu));
  const int nps = pev::pairs(st.data(), pu.data(), n, nst, npu, ps.data(), iw.data());
  std::vector<double> nh(z(nst), 0.0);
  {
    std::vector<int> sp(z(nps), -1);
    for (int i = 0; i < n; ++i) sp[z(ps[z(i)])] = st[z(i)];
    for (int a : sp) nh[z(a)] += 1.0;
  }
  for (double v : nh) {
    if (v < 2) fail("a stratum has one first-stage unit: its variance cannot be estimated.");
  }
  std::vector<double> fv;
  if (!Rf_isNull(fpc)) {
    SEXP fx = PROTECT(clab::num(fpc, "fpc", n, 0.0, true));
    ++np;
    fv.assign(z(nst), -1.0);
    for (int i = 0; i < n; ++i) {
      double &f = fv[z(st[z(i)])];
      if (f >= 0 && f != REAL(fx)[i]) fail("`fpc` must be constant within each stratum.");
      f = REAL(fx)[i];
    }
    for (int h = 0; h < nst; ++h) {
      if (fv[z(h)] < nh[z(h)]) fail("`fpc` must be at least the first-stage units sampled in its stratum.");
    }
  }
  // domains
  std::vector<int> dm(z(n), 0);
  int nd = 1;
  SEXP dlev = R_NilValue;
  if (!Rf_isNull(by)) {
    csp::lab l = csp::codes(by, "by", n);
    dm.assign(INTEGER(l.code), INTEGER(l.code) + n);
    nd = l.g;
    dlev = l.lev;
    np += 2;
  }
  const int r0 = Rf_isNull(ref) || Rf_isNull(by) ? 0 : clab::ref(ref, dlev, nd, "ref", "domain") - 1;
  const double df = nps - nst, tq = Rf_qt(0.5 + 0.5 * lv, df, 1, 0);
  std::vector<double> th(z(nd)), u(z(nd) * z(n), 0.0), cnt(z(nd), 0.0);
  for (int d = 0; d < nd; ++d) {
    double sw = 0, swy = 0, swx = 0, k = 0;
    for (int i = 0; i < n; ++i) {
      if (dm[z(i)] != d) continue;
      sw += wv[z(i)];
      swy += wv[z(i)] * yy[i];
      if (sk == 3) swx += wv[z(i)] * xv[z(i)];
      k += 1;
    }
    cnt[z(d)] = k;
    double *ud = u.data() + z(d) * z(n);
    double t;
    if (sk == 1) {
      t = swy;
      for (int i = 0; i < n; ++i) ud[i] = dm[z(i)] == d ? wv[z(i)] * yy[i] : 0.0;
    } else if (sk == 3) {
      if (swx == 0) fail("`den` sums to 0 in a domain: the ratio is undefined.");
      t = swy / swx;
      for (int i = 0; i < n; ++i) ud[i] = dm[z(i)] == d ? wv[z(i)] * (yy[i] - t * xv[z(i)]) / swx : 0.0;
    } else {
      // a mean; the variance is the mean of (y - ybar)^2 k / (k - 1)
      const double m = swy / sw;
      std::vector<double> zz(z(n));
      for (int i = 0; i < n; ++i) zz[z(i)] = sk == 4 ? (yy[i] - m) * (yy[i] - m) * k / (k - 1) : yy[i];
      if (sk == 4 && k < 2) fail("a domain has one unit: its variance is undefined.");
      t = 0;
      for (int i = 0; i < n; ++i) t += dm[z(i)] == d ? wv[z(i)] * zz[z(i)] / sw : 0.0;
      for (int i = 0; i < n; ++i) ud[i] = dm[z(i)] == d ? wv[z(i)] * (zz[z(i)] - t) / sw : 0.0;
    }
    th[z(d)] = t;
  }
  const double *fp = fv.empty() ? nullptr : fv.data();
  cols c(8);
  SEXP cd = c.str("dom", nd);
  double *ce = c.dbl("est", nd), *cs = c.dbl("se", nd), *cl = c.dbl("lo", nd),
         *ch = c.dbl("hi", nd), *cv = c.dbl("cv", nd), *cn = c.dbl("n", nd);
  double *cdf = c.dbl("df", nd);
  for (int d = 0; d < nd; ++d) {
    const double se = std::sqrt(pev::svar(u.data() + z(d) * z(n), n, st.data(), nst, ps.data(), nps, fp));
    SET_STRING_ELT(cd, d, Rf_isNull(dlev) ? Rf_mkChar("all") : STRING_ELT(dlev, d));
    ce[d] = th[z(d)];
    cs[d] = se;
    cl[d] = th[z(d)] - tq * se;
    ch[d] = th[z(d)] + tq * se;
    cv[d] = th[z(d)] != 0 ? 100 * se / std::fabs(th[z(d)]) : NA_REAL;
    cn[d] = cnt[z(d)];
    cdf[d] = df;
  }
  SEXP est = PROTECT(c.frame());
  ++np;
  SEXP cmp = R_NilValue;
  if (nd > 1) {
    cols k(5);
    SEXP kd = k.str("dom", nd - 1);
    double *kf = k.dbl("diff", nd - 1), *ks = k.dbl("se", nd - 1), *kt = k.dbl("stat", nd - 1),
           *kp = k.dbl("p", nd - 1);
    std::vector<double> ud(z(n));
    for (int d = 0, q = 0; d < nd; ++d) {
      if (d == r0) continue;
      for (int i = 0; i < n; ++i) ud[z(i)] = u[z(d) * z(n) + z(i)] - u[z(r0) * z(n) + z(i)];
      const double se = std::sqrt(pev::svar(ud.data(), n, st.data(), nst, ps.data(), nps, fp));
      SET_STRING_ELT(kd, q, STRING_ELT(dlev, d));
      kf[q] = th[z(d)] - th[z(r0)];
      ks[q] = se;
      kt[q] = kf[q] / se;
      kp[q] = 2 * Rf_pt(-std::fabs(kt[q]), df, 1, 0);
      ++q;
    }
    cmp = PROTECT(k.frame());
    ++np;
  }
  // the design: Kish's effective size and design effect of the weights
  double sw = 0, sw2 = 0;
  for (double v : wv) {
    sw += v;
    sw2 += v * v;
  }
  cols g(7);
  double *g1 = g.dbl("n", 1), *g2 = g.dbl("strata", 1), *g3 = g.dbl("psu", 1),
         *g4 = g.dbl("sumw", 1), *g5 = g.dbl("df", 1), *g6 = g.dbl("effn", 1),
         *g7 = g.dbl("deffw", 1);
  g1[0] = n;
  g2[0] = nst;
  g3[0] = nps;
  g4[0] = sw;
  g5[0] = df;
  g6[0] = sw * sw / sw2;
  g7[0] = n * sw2 / (sw * sw);
  SEXP des = PROTECT(g.frame());
  SEXP s1 = PROTECT(Rf_mkString(so[sk])), l1 = PROTECT(Rf_ScalarReal(lv));
  np += 3;
  SEXP out = PROTECT(pevr::nlist({{"est", est}, {"cmp", cmp}, {"des", des}, {"stat", s1},
                                  {"level", l1}}));
  Rf_classgets(out, PROTECT(Rf_mkString("svymx")));
  UNPROTECT(np + 2);
  return out;
}

}  // extern "C"
