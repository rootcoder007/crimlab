// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry point of the staggered-adoption estimators (stgdid): checks
// the arguments, holds the balanced panel wide, runs the estimator of
// stg.cpp and its cluster bootstrap, and builds the R result. Nothing
// calls back into the R evaluator.
#include <algorithm>
#include <cmath>
#include <vector>

#include "pevr.h"
#include "pevs.h"
#include "rchk.h"
#include "spin.h"

using clab::fail;
using pevr::cols;
using pevr::z;

extern "C" {

//' @srrstats {G2.6} stgdid() reads the outcome, units, periods and
//'   treatment as plain vectors of one length; covariates as for didfe().
//' @noRd
SEXP C_clab_r_stgdid(SEXP y, SEXP unit, SEXP time, SEXP first, SEXP trt, SEXP meth,
                     SEXP ctrl, SEXP adj, SEXP est, SEXP clus, SEXP reps, SEXP seed,
                     SEXP level) {
  clab::who = "stgdid";
  static const char *const mo[] = {"csa", "wool", "dcdh"};
  static const char *const co[] = {"notyet", "never"};
  static const char *const eo[] = {"dr", "ipw", "reg"};
  pev::stg s{};
  s.meth = csp::pick(meth, "meth", mo, 3);
  s.ctrl = csp::pick(ctrl, "ctrl", co, 2);
  s.est = csp::pick(est, "est", eo, 3);
  const int nr = clab::whole(reps, "reps", 0);
  const double lv = clab::one(level, "level", 0.0, 1.0);
  int np = 0;
  pevr::pan pn = pevr::panel(y, unit, time, clab::ninf, &np);
  const int n = pn.n, nu = pn.nu, nt = pn.nt;
  if (n != nu * nt) fail("the panel must be balanced: every unit in every period.");
  if (Rf_isNull(first) == Rf_isNull(trt)) fail("give one of `first` and `trt`.");
  // first treated periods, and the treatment of every unit and period
  std::vector<int> g(z(nu), nt), dw(z(nu) * z(nt), 0);
  if (!Rf_isNull(trt)) {
    SEXP dv = PROTECT(csp::bin01(trt, "trt", n));
    ++np;
    for (int i = 0; i < n; ++i) dw[z(pn.t[z(i)]) * z(nu) + z(pn.u[z(i)])] = INTEGER(dv)[i];
    if (s.meth != 2) g = pevr::firsts(pn, INTEGER(dv));
  } else {
    std::vector<double> fv = pevr::firstcol(first, pn);
    for (int a = 0; a < nu; ++a) {
      const double v = fv[z(a)];
      if (v > pn.tv.back()) continue;
      if (v <= pn.tv.front()) {
        g[z(a)] = 0;
      } else {
        const auto it = std::lower_bound(pn.tv.begin(), pn.tv.end(), v);
        if (*it != v) fail("`first` must be one of the periods of `time` (or NA, 0, Inf).");
        g[z(a)] = static_cast<int>(it - pn.tv.begin());
      }
      for (int b = g[z(a)]; b < nt; ++b) dw[z(b) * z(nu) + z(a)] = 1;
    }
  }
  // covariates, wide
  std::vector<pevr::term> x;
  pevr::adjcols(adj, n, x);
  s.k = static_cast<int>(x.size());
  if (s.k > 0 && s.meth != 0) fail("`adj` is taken by `meth = \"csa\"` only.");
  std::vector<double> xw(z(s.k) * z(nt) * z(nu));
  for (int j = 0; j < s.k; ++j) {
    for (int i = 0; i < n; ++i) {
      xw[(z(j) * z(nt) + z(pn.t[z(i)])) * z(nu) + z(pn.u[z(i)])] = x[z(j)].x[z(i)];
    }
  }
  std::vector<double> yw(z(n));
  for (int i = 0; i < n; ++i) yw[z(pn.t[z(i)]) * z(nu) + z(pn.u[z(i)])] = pn.y[i];
  // the clusters of the bootstrap, one per unit
  std::vector<int> cl(z(nu));
  int ncl = nu;
  if (Rf_isNull(clus)) {
    for (int a = 0; a < nu; ++a) cl[z(a)] = a;
  } else {
    csp::lab l = csp::codes(clus, "clus", n);
    std::vector<int> seen(z(nu), -1);
    for (int i = 0; i < n; ++i) {
      const int a = pn.u[z(i)], c = INTEGER(l.code)[i];
      if (seen[z(a)] >= 0 && seen[z(a)] != c) fail("`clus` must be constant within each unit.");
      seen[z(a)] = c;
      cl[z(a)] = c;
    }
    ncl = l.g;
    UNPROTECT(2);
    if (ncl < 2) fail("`clus` has one cluster: the bootstrap needs two or more.");
  }
  const std::uint64_t sd = nr > 0 ? csp::seed(seed) : 0u;
  s.yw = yw.data();
  s.dw = dw.data();
  s.xw = xw.data();
  s.g = g.data();
  s.nu = nu;
  s.nt = nt;
  // units treated from the first period have no period before: left out
  std::vector<int> idx;
  for (int a = 0; a < nu; ++a) {
    if (s.meth == 2 || g[z(a)] > 0) idx.push_back(a);
  }
  const double zq = Rf_qnorm5(0.5 + 0.5 * lv, 0.0, 1.0, 1, 0);
  SEXP att = R_NilValue, dyn = R_NilValue, ov;
  if (s.meth == 2) {
    double v[4], se[4] = {NA_REAL, NA_REAL, NA_REAL, NA_REAL};
    std::vector<int> ev;
    pev::stat(s, idx, ev, v);
    if (std::isnan(v[0])) fail("no unit switches treatment against a unit that did not.");
    if (nr > 0) pev::boot(s, cl.data(), ncl, nr, sd, ev, 4, se);
    cols c(6);
    SEXP tm = c.str("term", 2);
    double *e = c.dbl("est", 2), *sv = c.dbl("se", 2), *lo = c.dbl("lo", 2), *hi = c.dbl("hi", 2),
           *ns = c.dbl("n", 2);
    for (int q = 0; q < 2; ++q) {
      SET_STRING_ELT(tm, q, Rf_mkChar(q == 0 ? "effect" : "placebo"));
      e[q] = pevr::na(v[q]);
      sv[q] = pevr::na(se[q]);
      lo[q] = pevr::na(v[q] - zq * se[q]);
      hi[q] = pevr::na(v[q] + zq * se[q]);
      ns[q] = v[q + 2];
    }
    ov = PROTECT(c.frame());
    ++np;
  } else {
    pev::cells(g.data(), nu, nt, s.cg, s.ct, s.meth == 1 ? 1 : 0);
    if (s.cg.empty()) fail("no unit is first treated after the first period.");
    std::vector<int> ev;
    for (std::size_t c = 0; c < s.cg.size(); ++c) ev.push_back(s.ct[c] - s.cg[c]);
    std::sort(ev.begin(), ev.end());
    ev.erase(std::unique(ev.begin(), ev.end()), ev.end());
    const int nc = static_cast<int>(s.cg.size()), ne = static_cast<int>(ev.size());
    const int ns = nc + 1 + ne;
    std::vector<double> v(z(ns)), se(z(ns), NA_REAL);
    pev::stat(s, idx, ev, v.data());
    if (nr > 0) pev::boot(s, cl.data(), ncl, nr, sd, ev, ns, se.data());
    std::vector<double> ng(z(nt), 0.0);
    for (int a : idx) {
      if (g[z(a)] < nt) ng[z(g[z(a)])] += 1.0;
    }
    cols c(8);
    double *cc = c.dbl("cohort", nc), *cp = c.dbl("period", nc), *ce = c.dbl("evtm", nc),
           *e = c.dbl("est", nc), *sv = c.dbl("se", nc), *lo = c.dbl("lo", nc),
           *hi = c.dbl("hi", nc), *nn = c.dbl("ncoh", nc);
    for (int q = 0; q < nc; ++q) {
      cc[q] = pn.tv[z(s.cg[z(q)])];
      cp[q] = pn.tv[z(s.ct[z(q)])];
      ce[q] = s.ct[z(q)] - s.cg[z(q)];
      e[q] = pevr::na(v[z(q)]);
      sv[q] = pevr::na(se[z(q)]);
      lo[q] = pevr::na(v[z(q)] - zq * se[z(q)]);
      hi[q] = pevr::na(v[z(q)] + zq * se[z(q)]);
      nn[q] = ng[z(s.cg[z(q)])];
    }
    att = PROTECT(c.frame());
    cols o(5);
    SEXP tm = o.str("term", 1);
    SET_STRING_ELT(tm, 0, Rf_mkChar("att"));
    double *oe = o.dbl("est", 1), *os = o.dbl("se", 1), *ol = o.dbl("lo", 1), *oh = o.dbl("hi", 1);
    oe[0] = pevr::na(v[z(nc)]);
    os[0] = pevr::na(se[z(nc)]);
    ol[0] = pevr::na(v[z(nc)] - zq * se[z(nc)]);
    oh[0] = pevr::na(v[z(nc)] + zq * se[z(nc)]);
    ov = PROTECT(o.frame());
    cols d(5);
    double *de = d.dbl("evtm", ne), *dv = d.dbl("est", ne), *ds = d.dbl("se", ne),
           *dl = d.dbl("lo", ne), *dh = d.dbl("hi", ne);
    for (int q = 0; q < ne; ++q) {
      const std::size_t a = z(nc + 1 + q);
      de[q] = ev[z(q)];
      dv[q] = pevr::na(v[a]);
      ds[q] = pevr::na(se[a]);
      dl[q] = pevr::na(v[a] - zq * se[a]);
      dh[q] = pevr::na(v[a] + zq * se[a]);
    }
    dyn = PROTECT(d.frame());
    np += 3;
  }
  SEXP m1 = PROTECT(Rf_mkString(mo[s.meth])), c1 = PROTECT(Rf_mkString(co[s.ctrl])),
       r1 = PROTECT(Rf_ScalarInteger(nr)), l1 = PROTECT(Rf_ScalarReal(lv));
  np += 4;
  SEXP out = PROTECT(pevr::nlist({{"att", att}, {"overall", ov}, {"dynamic", dyn},
                                  {"meth", m1}, {"ctrl", c1}, {"reps", r1}, {"level", l1}}));
  ++np;
  Rf_classgets(out, PROTECT(Rf_mkString("stgdid")));
  UNPROTECT(np + 1);
  return out;
}

}  // extern "C"
