// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the synthetic control (synthc) and synthetic
// difference in differences (sdidx): check the arguments, hold the
// balanced panel wide, solve the weights (syn.cpp), run the placebos or
// the jackknife, and build the R result. Nothing calls back into the R
// evaluator.
#include <algorithm>
#include <cmath>
#include <vector>

#include "kutil.h"
#include "pev.h"
#include "pevr.h"
#include "rchk.h"
#include "spin.h"

using clab::fail;
using pevr::cols;
using pevr::z;

namespace {

// The balanced panel wide ([s nu + i]), the adoption period of the
// treated units (a code, all the same), and the treated and comparison
// units.
struct wide {
  std::vector<double> yw;
  std::vector<int> trd, ctl;
  int t0 = 0;
};

wide widen(const pevr::pan &pn, SEXP first) {
  if (pn.n != pn.nu * pn.nt) fail("the panel must be balanced: every unit in every period.");
  wide w;
  w.yw.resize(z(pn.n));
  for (int i = 0; i < pn.n; ++i) w.yw[z(pn.t[z(i)]) * z(pn.nu) + z(pn.u[z(i)])] = pn.y[i];
  const std::vector<double> fv = pevr::firstcol(first, pn);
  double g = R_PosInf;
  for (int a = 0; a < pn.nu; ++a) {
    if (!std::isfinite(fv[z(a)])) {
      w.ctl.push_back(a);
      continue;
    }
    if (std::isfinite(g) && fv[z(a)] != g) fail("the treated units must share one adoption period.");
    g = fv[z(a)];
    w.trd.push_back(a);
  }
  if (w.trd.empty()) fail("no unit is treated: `first` is NA, 0 or Inf throughout.");
  if (w.ctl.size() < 2) fail("the comparison needs two units never treated at least.");
  const auto it = std::find(pn.tv.begin(), pn.tv.end(), g);
  if (it == pn.tv.begin() || it == pn.tv.end()) {
    fail("the adoption period must be a period of `time` after the first.");
  }
  w.t0 = static_cast<int>(it - pn.tv.begin());
  return w;
}

// The synthetic control of unit i from the units d: the weights and the
// gap in every period. false when the weights fail.
bool synth(const wide &w, int nu, int nt, int i, const std::vector<int> &d, double *om,
           double *gap) {
  const int n0 = static_cast<int>(d.size());
  std::vector<double> a(z(w.t0) * z(n0)), b(z(w.t0)), g(z(n0) * z(n0));
  for (int s = 0; s < w.t0; ++s) {
    b[z(s)] = w.yw[z(s) * z(nu) + z(i)];
    for (int q = 0; q < n0; ++q) a[z(q) * z(w.t0) + z(s)] = w.yw[z(s) * z(nu) + z(d[z(q)])];
  }
  pev::gram(a.data(), b.data(), w.t0, n0, 0, 0.0, g.data());
  if (pev::mnp(g.data(), n0, om, 10000) < 0) return false;
  for (int s = 0; s < nt; ++s) {
    double v = w.yw[z(s) * z(nu) + z(i)];
    for (int q = 0; q < n0; ++q) v -= om[q] * w.yw[z(s) * z(nu) + z(d[z(q)])];
    gap[s] = v;
  }
  return true;
}

// The root mean squared gap before (pre) or from (post) period t0.
double rmspe(const double *gap, int t0, int nt, bool post) {
  double s = 0.0;
  const int lo = post ? t0 : 0, hi = post ? nt : t0;
  for (int q = lo; q < hi; ++q) s += gap[q] * gap[q];
  return std::sqrt(s / (hi - lo));
}

const char *const wfail = "the weights could not be solved (the units' paths are degenerate).";

}  // namespace

extern "C" {

//' @srrstats {G2.6} synthc() and sdidx() read the outcome, units, periods
//'   and adoption periods as plain vectors of one length.
//' @noRd
SEXP C_clab_r_synthc(SEXP y, SEXP unit, SEXP time, SEXP first, SEXP plac) {
  clab::who = "synthc";
  const bool dp = csp::flag(plac, "plac");
  int np = 0;
  pevr::pan pn = pevr::panel(y, unit, time, clab::ninf, &np);
  const wide w = widen(pn, first);
  if (w.trd.size() != 1) fail("one unit must be treated: `first` is finite for %d.", static_cast<int>(w.trd.size()));
  const int nu = pn.nu, nt = pn.nt, n0 = static_cast<int>(w.ctl.size()), i1 = w.trd[0];
  std::vector<double> om(z(n0)), gap(z(nt));
  if (!synth(w, nu, nt, i1, w.ctl, om.data(), gap.data())) fail("%s", wfail);
  const double pre = rmspe(gap.data(), w.t0, nt, false), post = rmspe(gap.data(), w.t0, nt, true);
  const double ratio = post / pre;
  // placebos: each comparison unit in turn, the others its donors
  SEXP pl = R_NilValue;
  double pv = NA_REAL;
  if (dp) {
    cols c(2);
    SEXP cu = c.str("unit", n0);
    double *cr = c.dbl("ratio", n0);
    std::vector<double> o2(z(n0 - 1)), g2(z(nt));
    int ge = 0;
    for (int q = 0; q < n0; ++q) {
      std::vector<int> d;
      for (int r = 0; r < n0; ++r) {
        if (r != q) d.push_back(w.ctl[z(r)]);
      }
      const bool ok = synth(w, nu, nt, w.ctl[z(q)], d, o2.data(), g2.data());
      const double rq = ok ? rmspe(g2.data(), w.t0, nt, true) / rmspe(g2.data(), w.t0, nt, false) : NA_REAL;
      SET_STRING_ELT(cu, q, STRING_ELT(pn.ulev, w.ctl[z(q)]));
      cr[q] = rq;
      ge += !ISNAN(rq) && rq >= ratio ? 1 : 0;
    }
    pl = PROTECT(c.frame());
    ++np;
    pv = (1.0 + ge) / (n0 + 1.0);
  }
  cols cw(2);
  SEXP wu = cw.str("unit", n0);
  double *ww = cw.dbl("w", n0);
  for (int q = 0; q < n0; ++q) {
    SET_STRING_ELT(wu, q, STRING_ELT(pn.ulev, w.ctl[z(q)]));
    ww[q] = om[z(q)];
  }
  SEXP wt = PROTECT(cw.frame());
  cols cp(4);
  double *pt = cp.dbl("time", nt), *py = cp.dbl("treated", nt), *ps = cp.dbl("synth", nt),
         *pg = cp.dbl("gap", nt);
  for (int s = 0; s < nt; ++s) {
    pt[s] = pn.tv[z(s)];
    py[s] = w.yw[z(s) * z(nu) + z(i1)];
    ps[s] = py[s] - gap[z(s)];
    pg[s] = gap[z(s)];
  }
  SEXP path = PROTECT(cp.frame());
  cols cf(4);
  double *f1 = cf.dbl("rmspepre", 1), *f2 = cf.dbl("rmspepost", 1), *f3 = cf.dbl("ratio", 1),
         *f4 = cf.dbl("p", 1);
  f1[0] = pre;
  f2[0] = post;
  f3[0] = pevr::na(ratio);
  f4[0] = pv;
  SEXP fit = PROTECT(cf.frame());
  np += 3;
  SEXP tu = PROTECT(Rf_ScalarString(STRING_ELT(pn.ulev, i1)));
  ++np;
  SEXP out = PROTECT(pevr::nlist({{"weights", wt}, {"path", path}, {"fit", fit}, {"plac", pl},
                                  {"unit", tu}}));
  ++np;
  Rf_classgets(out, PROTECT(Rf_mkString("synthc")));
  UNPROTECT(np + 1);
  return out;
}

SEXP C_clab_r_sdidx(SEXP y, SEXP unit, SEXP time, SEXP first, SEXP reps, SEXP seed,
                    SEXP level) {
  clab::who = "sdidx";
  const int nr = clab::whole(reps, "reps", 2);
  const double lv = clab::one(level, "level", 0.0, 1.0);
  int np = 0;
  pevr::pan pn = pevr::panel(y, unit, time, clab::ninf, &np);
  const wide w = widen(pn, first);
  const int nu = pn.nu, nt = pn.nt, n0 = static_cast<int>(w.ctl.size());
  const int n1 = static_cast<int>(w.trd.size());
  std::vector<double> om(z(n0)), la(z(w.t0));
  const double tau = pev::sdtau(w.yw.data(), nu, nt, w.ctl, w.trd, w.t0, true, om.data(), la.data());
  if (std::isnan(tau)) fail("%s", wfail);
  double se = NA_REAL;
  const char *meth;
  if (n1 == 1) {
    // placebo: n1 comparison units at random treated as if adopting
    meth = "placebo";
    kutil::rng g(csp::seed(seed));
    std::vector<double> o2(z(n0 - 1)), l2(z(w.t0));
    kutil::ksum s1, s2;
    double k = 0;
    for (int r = 0; r < nr; ++r) {
      std::vector<int> c = w.ctl;
      const int pick = g.below(n0);
      const std::vector<int> t1 = {c[z(pick)]};
      c.erase(c.begin() + pick);
      const double v = pev::sdtau(w.yw.data(), nu, nt, c, t1, w.t0, true, o2.data(), l2.data());
      if (std::isnan(v)) continue;
      s1.add(v);
      s2.add(v * v);
      k += 1;
    }
    if (k > 1) {
      const double m = s1.val() / k;
      se = std::sqrt((nr - 1.0) / nr * std::max(0.0, (s2.val() - k * m * m) / (k - 1)));
    }
  } else {
    // the jackknife over units with the weights held
    meth = "jackknife";
    std::vector<double> u;
    for (int a = 0; a < n0 + n1; ++a) {
      std::vector<int> c, t;
      std::vector<double> o2;
      double sw = 0.0;
      for (int q = 0; q < n0; ++q) {
        if (q == a) continue;
        c.push_back(w.ctl[z(q)]);
        o2.push_back(om[z(q)]);
        sw += om[z(q)];
      }
      for (int q = 0; q < n1; ++q) {
        if (q + n0 != a) t.push_back(w.trd[z(q)]);
      }
      if (!(sw > 0)) {
        u.clear();
        break;
      }
      for (double &v : o2) v /= sw;
      u.push_back(pev::sdtau(w.yw.data(), nu, nt, c, t, w.t0, false, o2.data(), la.data()));
    }
    if (!u.empty()) {
      const double n = static_cast<double>(u.size());
      double m = 0.0, ss = 0.0;
      for (double v : u) m += v / n;
      for (double v : u) ss += (v - m) * (v - m);
      se = std::sqrt((n - 1.0) / n * ss);
    }
  }
  const double zq = Rf_qnorm5(0.5 + 0.5 * lv, 0.0, 1.0, 1, 0);
  cols ce(4);
  double *e1 = ce.dbl("est", 1), *e2 = ce.dbl("se", 1), *e3 = ce.dbl("lo", 1), *e4 = ce.dbl("hi", 1);
  e1[0] = tau;
  e2[0] = se;
  e3[0] = ISNAN(se) ? NA_REAL : tau - zq * se;
  e4[0] = ISNAN(se) ? NA_REAL : tau + zq * se;
  SEXP est = PROTECT(ce.frame());
  cols cu(2);
  SEXP uu = cu.str("unit", n0);
  double *uw = cu.dbl("w", n0);
  for (int q = 0; q < n0; ++q) {
    SET_STRING_ELT(uu, q, STRING_ELT(pn.ulev, w.ctl[z(q)]));
    uw[q] = om[z(q)];
  }
  SEXP unitw = PROTECT(cu.frame());
  cols ct(2);
  double *tt = ct.dbl("time", w.t0), *tw = ct.dbl("w", w.t0);
  for (int s = 0; s < w.t0; ++s) {
    tt[s] = pn.tv[z(s)];
    tw[s] = la[z(s)];
  }
  SEXP timew = PROTECT(ct.frame());
  cols cp(3);
  double *pt = cp.dbl("time", nt), *py = cp.dbl("treated", nt), *ps = cp.dbl("synth", nt);
  for (int s = 0; s < nt; ++s) {
    pt[s] = pn.tv[z(s)];
    py[s] = 0.0;
    for (int i : w.trd) py[s] += w.yw[z(s) * z(nu) + z(i)] / n1;
    ps[s] = 0.0;
    for (int q = 0; q < n0; ++q) ps[s] += om[z(q)] * w.yw[z(s) * z(nu) + z(w.ctl[z(q)])];
  }
  SEXP path = PROTECT(cp.frame());
  SEXP m1 = PROTECT(Rf_mkString(meth)), l1 = PROTECT(Rf_ScalarReal(lv));
  np += 6;
  SEXP out = PROTECT(pevr::nlist({{"est", est}, {"meth", m1}, {"unitw", unitw}, {"timew", timew},
                                  {"path", path}, {"level", l1}}));
  ++np;
  Rf_classgets(out, PROTECT(Rf_mkString("sdidx")));
  UNPROTECT(np + 1);
  return out;
}

}  // extern "C"
