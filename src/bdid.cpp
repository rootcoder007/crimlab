// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the difference-in-differences families: the
// two-way fixed-effects regression (didfe) and the event study (evtst).
// Each checks its arguments, builds the kernel's arrays, calls the kernels
// (did.cpp, and cnt.cpp for counts) and builds the R result. Distribution
// functions come from R's maths library. Nothing calls back into the R
// evaluator.
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

#include "crimlab.h"
#include "pev.h"
#include "pevr.h"
#include "rchk.h"
#include "spin.h"

extern "C" double Rf_pt(double, double, int, int);
extern "C" double Rf_qt(double, double, int, int);
extern "C" double Rf_pf(double, double, double, int, int);
extern "C" double Rf_pnorm5(double, double, double, int, int);

using clab::fail;
using pevr::cols;
using pevr::nlist;
using pevr::z;

namespace {

const double qnan = std::numeric_limits<double>::quiet_NaN();

// One fit of the two-way model on the rows idx.
struct fit {
  bool ok = false;
  std::vector<double> b, v;  // p, p x p
  double g = 0, kk = 0, n = 0, theta = qnan;
  int nu = 0, nt = 0, drop = 0;
};

// The model of didfe() and evtst(): the panel, the regressors, the
// clusters and the options.
struct model {
  const double *y;
  const pevr::pan *pn;
  const std::vector<pevr::term> *x;
  const int *cl, *cl2;
  const double *w;    // NULL or n weights
  const double *off;  // NULL or n log exposures
  int fe, fam;        // fam: 0 least squares, 1 Poisson, 2 negative binomial
};

// Codes from 0, in order of first appearance, of src over the rows idx.
int rec(const int *src, const std::vector<int> &idx, std::vector<int> &out) {
  std::vector<int> map;
  out.resize(idx.size());
  int k = 0;
  for (std::size_t q = 0; q < idx.size(); ++q) {
    const int s = src[idx[q]];
    if (s >= static_cast<int>(map.size())) map.resize(z(s) + 1, -1);
    if (map[z(s)] < 0) map[z(s)] = k++;
    out[q] = map[z(s)];
  }
  return k;
}

// fixest's "nested" count of the fixed effects (did.cpp has the same for
// least squares): an effect nested in a cluster costs nothing.
bool nests(const std::vector<int> &f, int nf, const std::vector<int> &c) {
  std::vector<int> seen(z(nf), -1);
  for (std::size_t i = 0; i < f.size(); ++i) {
    if (seen[z(f[i])] < 0) seen[z(f[i])] = c[i];
    if (seen[z(f[i])] != c[i]) return false;
  }
  return true;
}

fit run(const model &m, std::vector<int> idx) {
  fit r;
  const int p = static_cast<int>(m.x->size());
  if (m.fam > 0) {
    // counts: units and periods with no count carry nothing; drop them
    // until none is left
    for (;;) {
      std::vector<double> su(z(m.pn->nu), 0.0), st(z(m.pn->nt), 0.0);
      for (int i : idx) {
        su[z(m.pn->u[z(i)])] += m.y[i];
        st[z(m.pn->t[z(i)])] += m.y[i];
      }
      std::vector<int> keep;
      for (int i : idx) {
        if ((su[z(m.pn->u[z(i)])] > 0 || !m.fe) && st[z(m.pn->t[z(i)])] > 0) keep.push_back(i);
      }
      if (keep.size() == idx.size()) break;
      idx.swap(keep);
    }
  }
  const int n = static_cast<int>(idx.size());
  std::vector<int> u, t, c1, c2;
  const int nu = rec(m.pn->u.data(), idx, u), nt = rec(m.pn->t.data(), idx, t);
  const int ncl = rec(m.cl, idx, c1), ncl2 = m.cl2 ? rec(m.cl2, idx, c2) : 0;
  r.nu = nu;
  r.nt = nt;
  r.n = n;
  if (n < 2 || ncl < 2 || (m.cl2 && ncl2 < 2)) return r;
  std::vector<double> y(z(n)), w, x(z(n) * z(p));
  for (int q = 0; q < n; ++q) y[z(q)] = m.y[idx[z(q)]];
  if (m.w) {
    w.resize(z(n));
    for (int q = 0; q < n; ++q) w[z(q)] = m.w[idx[z(q)]];
  }
  for (int j = 0; j < p; ++j) {
    for (int q = 0; q < n; ++q) x[z(j) * z(n) + z(q)] = (*m.x)[z(j)].x[z(idx[z(q)])];
  }
  r.b.assign(z(p), qnan);
  r.v.assign(z(p) * z(p), qnan);
  if (m.fam == 0) {
    const std::size_t sp = z(p), sn = z(n);
    std::vector<double> work(sn * (sp + 2) + sp * sp + 2 * z(nu + nt) + (sn + 2 * sp) * sp + sp);
    std::vector<int> iwork(sn + std::max<std::size_t>({z(nu), z(nt), z(ncl) * z(ncl2)}));
    double st[5];
    const int rc = clab_didfe(y.data(), x.data(), n, p, u.data(), nu, t.data(), nt, m.fe,
                              m.w ? w.data() : nullptr, c1.data(), ncl,
                              m.cl2 ? c2.data() : nullptr, ncl2, 1e-10, 10000, r.b.data(), r.v.data(),
                              st, work.data(), iwork.data());
    if (rc != CLAB_OK) return r;
    r.g = st[1];
    r.kk = st[2];
    r.ok = true;
    return r;
  }
  // counts: the dummies of the effects after an intercept and the
  // regressors, the fit and its sandwich by cnt.cpp
  const int pd = 1 + p + (m.fe ? nu - 1 : 0) + nt - 1;
  std::vector<double> xd(z(n) * z(pd), 0.0), off;
  for (int q = 0; q < n; ++q) xd[z(q)] = 1.0;
  std::copy(x.begin(), x.end(), xd.begin() + n);
  int col = 1 + p;
  if (m.fe) {
    for (int a = 1; a < nu; ++a, ++col) {
      for (int q = 0; q < n; ++q) xd[z(col) * z(n) + z(q)] = u[z(q)] == a ? 1.0 : 0.0;
    }
  }
  for (int a = 1; a < nt; ++a, ++col) {
    for (int q = 0; q < n; ++q) xd[z(col) * z(n) + z(q)] = t[z(q)] == a ? 1.0 : 0.0;
  }
  if (m.off) {
    off.resize(z(n));
    for (int q = 0; q < n; ++q) off[z(q)] = m.off[idx[z(q)]];
  }
  std::vector<int> c12(z(n)), pw(z(ncl) * z(std::max(ncl2, 1)));
  const int n12 = m.cl2 ? pev::pairs(c1.data(), c2.data(), n, ncl, ncl2, c12.data(), pw.data()) : 0;
  const int big = std::max({ncl, ncl2, n12});
  const std::size_t sp = z(pd), sn = z(n);
  std::vector<double> work(sn * sp + 8 * sn + sp + sp * sp + 2 * z(big) * sp), beta(sp),
      vb(sp * sp), vc(sp * sp), vs(sp * sp, 0.0), st(8);
  std::vector<int> alias(sp), iwork(sp + z(big) + 1 + sn);
  const int *cls[3] = {c1.data(), c2.data(), c12.data()};
  const int ncs[3] = {ncl, ncl2, n12};
  const double sg[3] = {1.0, 1.0, -1.0};
  for (int k = 0; k < (m.cl2 ? 3 : 1); ++k) {
    const int rc = clab_cntreg(y.data(), xd.data(), n, pd, m.off ? off.data() : nullptr, cls[k],
                               ncs[k], m.fam == 2 ? 1 : 0, 25, 1e-8, beta.data(), vb.data(),
                               vc.data(), st.data(), alias.data(), work.data(), iwork.data());
    if (rc != CLAB_OK) return r;
    for (int j = 1; j <= p; ++j) {
      if (alias[z(j)]) return r;
    }
    for (std::size_t q = 0; q < sp * sp; ++q) vs[q] += sg[k] * vc[q];
  }
  r.g = m.cl2 ? std::min(ncl, ncl2) : ncl;
  r.kk = p + 1.0;
  if (m.fe && !(nests(u, nu, c1) || (m.cl2 && nests(u, nu, c2)))) r.kk += nu - 1.0;
  if (!(nests(t, nt, c1) || (m.cl2 && nests(t, nt, c2)))) r.kk += nt - 1.0;
  const double c = r.g / (r.g - 1) * (n - 1.0) / (n - r.kk);
  for (int j = 0; j < p; ++j) {
    r.b[z(j)] = beta[z(j) + 1];
    for (int k = 0; k < p; ++k) r.v[z(k) * z(p) + z(j)] = c * vs[z(k + 1) * sp + z(j) + 1];
  }
  r.theta = m.fam == 2 ? st[2] : qnan;
  r.drop = m.pn->nu - nu;
  r.ok = true;
  return r;
}

// The coefficient table of a fit: est, se, lo, hi, stat, p; with t on G -
// 1 degrees of freedom (least squares) or the normal (counts); for counts
// also the rate ratios.
SEXP ctab(const fit &f, const std::vector<pevr::term> &x, int p0, int p, bool cnt, double lv) {
  const int k = p - p0;
  cols c(cnt ? 10 : 7);
  SEXP tm = c.str("term", k);
  double *e = c.dbl("est", k), *s = c.dbl("se", k), *lo = c.dbl("lo", k), *hi = c.dbl("hi", k),
         *st = c.dbl("stat", k), *pv = c.dbl("p", k);
  double *rr = nullptr, *rl = nullptr, *rh = nullptr;
  if (cnt) {
    rr = c.dbl("rr", k);
    rl = c.dbl("rrlo", k);
    rh = c.dbl("rrhi", k);
  }
  const double df = f.g - 1;
  const double q = cnt ? Rf_qnorm5(0.5 + 0.5 * lv, 0.0, 1.0, 1, 0) : Rf_qt(0.5 + 0.5 * lv, df, 1, 0);
  for (int j = 0; j < k; ++j) {
    const int a = j + p0;
    SET_STRING_ELT(tm, j, Rf_mkCharCE(x[z(a)].nm.c_str(), CE_UTF8));
    const double b = f.b[z(a)], v = f.v[z(a) * z(p) + z(a)];
    const double se = v >= 0 ? std::sqrt(v) : qnan;
    e[j] = b;
    s[j] = pevr::na(se);
    lo[j] = pevr::na(b - q * se);
    hi[j] = pevr::na(b + q * se);
    st[j] = pevr::na(b / se);
    pv[j] = pevr::na(cnt ? 2 * Rf_pnorm5(-std::fabs(b / se), 0.0, 1.0, 1, 0)
                         : 2 * Rf_pt(-std::fabs(b / se), df, 1, 0));
    if (cnt) {
      rr[j] = std::exp(b);
      rl[j] = pevr::na(std::exp(b - q * se));
      rh[j] = pevr::na(std::exp(b + q * se));
    }
  }
  return c.frame();
}

// The clusters: NULL (the units), one vector of labels, or a list (data
// frame) of two such vectors for two-way clustering.
int clusters(SEXP clus, const pevr::pan &pn, std::vector<int> &c1, std::vector<int> &c2,
             int *ncl, int *ncl2) {
  const int n = pn.n;
  *ncl2 = 0;
  if (Rf_isNull(clus)) {
    c1 = pn.u;
    *ncl = pn.nu;
    return 0;
  }
  SEXP a = clus, b = R_NilValue;
  if (TYPEOF(clus) == VECSXP) {
    if (XLENGTH(clus) < 1 || XLENGTH(clus) > 2) {
      fail("`clus` must be one vector of cluster labels or a list of two.");
    }
    a = VECTOR_ELT(clus, 0);
    if (XLENGTH(clus) == 2) b = VECTOR_ELT(clus, 1);
  }
  csp::lab l = csp::codes(a, "clus", n);
  c1.assign(INTEGER(l.code), INTEGER(l.code) + n);
  *ncl = l.g;
  UNPROTECT(2);
  if (!Rf_isNull(b)) {
    csp::lab m = csp::codes(b, "clus", n);
    c2.assign(INTEGER(m.code), INTEGER(m.code) + n);
    *ncl2 = m.g;
    UNPROTECT(2);
  }
  if (*ncl < 2 || (!Rf_isNull(b) && *ncl2 < 2)) {
    fail("`clus` has one cluster: a cluster-robust variance needs two or more.");
  }
  if (*ncl < 10 || (!Rf_isNull(b) && *ncl2 < 10)) {
    Rf_warningcall(R_NilValue,
                   "%s: fewer than ten clusters; cluster-robust errors are unreliable with so "
                   "few (wcbx() is made for this case).",
                   clab::who);
  }
  return Rf_isNull(b) ? 1 : 2;
}

// The pre-period mean of the comparison units: the units never treated
// in the periods before the first treatment, else every untreated row
// before it.
double ctlmean(const pevr::pan &pn, const int *d, const double *w) {
  std::vector<char> ever(z(pn.nu), 0);
  int first = pn.nt;
  for (int i = 0; i < pn.n; ++i) {
    if (d[i]) {
      ever[z(pn.u[z(i)])] = 1;
      first = std::min(first, pn.t[z(i)]);
    }
  }
  for (int pass = 0; pass < 2; ++pass) {
    double s = 0.0, sw = 0.0;
    for (int i = 0; i < pn.n; ++i) {
      const bool in = pn.t[z(i)] < first && (pass == 0 ? !ever[z(pn.u[z(i)])] : d[i] == 0);
      if (!in) continue;
      const double wi = w ? w[i] : 1.0;
      s += wi * pn.y[i];
      sw += wi;
    }
    if (sw > 0) return s / sw;
  }
  return NA_REAL;
}

const char *const absorbed =
    "the fit failed: a regressor is constant within units and periods (the fixed effects "
    "absorb it), or the regressors are collinear.";

}  // namespace

extern "C" {

//' @srrstats {G2.6} didfe() and evtst() read the outcome, units, periods,
//'   treatment and weights as plain vectors of one length; adjustment
//'   columns come as a data frame, a list, a matrix or one vector.
//' @srrstats {G2.4c} Units and clusters are read as labels; periods as
//'   numbers, sorted by value.
//' @noRd
SEXP C_clab_r_didfe(SEXP y, SEXP unit, SEXP time, SEXP trt, SEXP adj, SEXP clus, SEXP w,
                    SEXP fe, SEXP fam, SEXP expo, SEXP loo, SEXP bacon, SEXP mtch, SEXP exact,
                    SEXP calip, SEXP level) {
  clab::who = "didfe";
  static const char *const feo[] = {"both", "time"};
  static const char *const fmo[] = {"ols", "pois", "nb"};
  const int fx = csp::pick(fe, "fe", feo, 2) == 0 ? 1 : 0;
  const int fm = csp::pick(fam, "fam", fmo, 3);
  const double lv = clab::one(level, "level", 0.0, 1.0);
  const bool dloo = csp::flag(loo, "loo"), dbac = csp::flag(bacon, "bacon");
  const int ratio = clab::whole(mtch, "mtch", 0);
  int np = 0;
  pevr::pan pn = pevr::panel(y, unit, time, fm > 0 ? 0.0 : clab::ninf, &np);
  const int n = pn.n;
  SEXP dv = PROTECT(csp::bin01(trt, "trt", n));
  ++np;
  const int *d = INTEGER(dv);
  int nd = 0;
  for (int i = 0; i < n; ++i) nd += d[i];
  if (nd == 0 || nd == n) fail("`trt` must be 1 in some rows and 0 in others.");
  std::vector<double> wv;
  if (!Rf_isNull(w)) {
    SEXP wx = PROTECT(clab::num(w, "w", n, 0.0));
    ++np;
    wv.assign(REAL(wx), REAL(wx) + n);
  }
  std::vector<double> off;
  if (!Rf_isNull(expo)) {
    if (fm == 0) fail("`expo` is an exposure for counts: give it with `fam = \"pois\"` or \"nb\".");
    SEXP e = PROTECT(clab::num(expo, "expo", n, 0.0, true));
    ++np;
    off.resize(z(n));
    for (int i = 0; i < n; ++i) off[z(i)] = std::log(REAL(e)[i]);
  }
  std::vector<int> c1, c2;
  int ncl = 0, ncl2 = 0;
  const int kc = clusters(clus, pn, c1, c2, &ncl, &ncl2);
  std::vector<pevr::term> x;
  {
    pevr::term tr;
    tr.nm = "trt";
    tr.x.assign(d, d + n);
    x.push_back(tr);
  }
  pevr::adjcols(adj, n, x);
  const int p = static_cast<int>(x.size());
  if ((dbac || ratio > 0) && fm != 0) fail("`bacon` and `mtch` work with `fam = \"ols\"`.");
  if (dbac && (p > 1 || !wv.empty() || !fx)) {
    fail("`bacon = TRUE` decomposes the plain two-way estimate: no `adj`, no `w`, `fe = \"both\"`.");
  }
  if (!Rf_isNull(exact) && ratio == 0) fail("`exact` groups the matching: give `mtch` too.");
  if (!Rf_isNull(calip) && ratio == 0) fail("`calip` limits the matching: give `mtch` too.");
  // each unit's first treated period (nt for never), when needed
  std::vector<int> g;
  if (dbac || ratio > 0) {
    if (n != pn.nu * pn.nt) fail("`bacon` and `mtch` need a balanced panel: every unit in every period.");
    g = pevr::firsts(pn, d);
  }
  model m{pn.y, &pn, &x, c1.data(), kc == 2 ? c2.data() : nullptr,
          wv.empty() ? nullptr : wv.data(), off.empty() ? nullptr : off.data(), fx, fm};
  // matching on the outcomes before the first adoption
  SEXP mt = R_NilValue;
  std::vector<double> mw;
  if (ratio > 0) {
    std::vector<int> tr(z(pn.nu));
    int g0 = pn.nt;
    for (int a = 0; a < pn.nu; ++a) {
      tr[z(a)] = g[z(a)] < pn.nt ? 1 : 0;
      if (tr[z(a)]) g0 = std::min(g0, g[z(a)]);
    }
    if (std::count(tr.begin(), tr.end(), 1) == pn.nu) fail("`mtch` needs units never treated to match to.");
    if (g0 == 0) fail("`mtch` needs a period before the first adoption to match on.");
    std::vector<int> ex(z(pn.nu), 0);
    if (!Rf_isNull(exact)) {
      csp::lab l = csp::codes(exact, "exact", n);
      std::vector<int> seen(z(pn.nu), -1);
      for (int i = 0; i < n; ++i) {
        const int a = pn.u[z(i)], e = INTEGER(l.code)[i];
        if (seen[z(a)] >= 0 && seen[z(a)] != e) fail("`exact` must be constant within each unit.");
        seen[z(a)] = e;
        ex[z(a)] = e;
      }
      UNPROTECT(2);
    }
    double cal = 0.0;
    if (!Rf_isNull(calip)) cal = clab::one(calip, "calip", 0.0, R_PosInf);
    std::vector<double> f(z(pn.nu) * z(g0));
    for (int i = 0; i < n; ++i) {
      if (pn.t[z(i)] < g0) f[z(pn.t[z(i)]) * z(pn.nu) + z(pn.u[z(i)])] = pn.y[i];
    }
    const int n1 = static_cast<int>(std::count(tr.begin(), tr.end(), 1));
    std::vector<int> ma(z(n1) * z(ratio)), mc(z(n1) * z(ratio));
    std::vector<double> md(z(n1) * z(ratio)), wu(z(pn.nu)), work(z(g0) * z(g0) + 3 * z(g0));
    const int nm = pev::mtch(f.data(), pn.nu, g0, tr.data(), ex.data(), ratio, cal, ma.data(),
                             mc.data(), md.data(), wu.data(), work.data());
    if (nm < 0) {
      fail("the pre-period outcomes are collinear across units: too many periods before the "
           "first adoption for the units to match on.");
    }
    if (nm == 0) fail("no treated unit has a comparison unit within its `exact` group and the caliper.");
    cols c(3);
    SEXP ct = c.str("treated", nm), cc = c.str("control", nm);
    double *cd = c.dbl("dist", nm);
    for (int q = 0; q < nm; ++q) {
      SET_STRING_ELT(ct, q, STRING_ELT(pn.ulev, ma[z(q)]));
      SET_STRING_ELT(cc, q, STRING_ELT(pn.ulev, mc[z(q)]));
      cd[q] = md[z(q)];
    }
    mt = PROTECT(c.frame());
    ++np;
    mw.resize(z(n));
    for (int i = 0; i < n; ++i) mw[z(i)] = wu[z(pn.u[z(i)])] * (wv.empty() ? 1.0 : wv[z(i)]);
    m.w = mw.data();
  }
  std::vector<int> all;
  for (int i = 0; i < n; ++i) {
    if (!m.w || m.w[i] > 0) all.push_back(i);
  }
  const fit f = run(m, all);
  if (!f.ok) fail("%s", absorbed);
  SEXP co = PROTECT(ctab(f, x, 0, p, fm > 0, lv));
  ++np;
  // leave one treated unit out at a time
  SEXP lo = R_NilValue;
  if (dloo) {
    std::vector<char> ever(z(pn.nu), 0);
    for (int i : all) {
      if (d[i]) ever[z(pn.u[z(i)])] = 1;
    }
    std::vector<int> us;
    for (int a = 0; a < pn.nu; ++a) {
      if (ever[z(a)]) us.push_back(a);
    }
    const int k = static_cast<int>(us.size());
    cols c(3);
    SEXP cu = c.str("unit", k);
    double *ce = c.dbl("est", k), *cs = c.dbl("se", k);
    for (int q = 0; q < k; ++q) {
      std::vector<int> idx;
      for (int i : all) {
        if (pn.u[z(i)] != us[z(q)]) idx.push_back(i);
      }
      const fit h = run(m, idx);
      SET_STRING_ELT(cu, q, STRING_ELT(pn.ulev, us[z(q)]));
      ce[q] = h.ok ? h.b[0] : NA_REAL;
      cs[q] = h.ok ? pevr::na(std::sqrt(h.v[0])) : NA_REAL;
    }
    lo = PROTECT(c.frame());
    ++np;
  }
  // the Goodman-Bacon decomposition
  SEXP bc = R_NilValue;
  if (dbac) {
    std::vector<double> yw(z(n));
    for (int i = 0; i < n; ++i) yw[z(pn.t[z(i)]) * z(pn.nu) + z(pn.u[z(i)])] = pn.y[i];
    const int mx = pn.nt * pn.nt + pn.nt;
    std::vector<int> ty(z(mx)), tr(z(mx)), cg(z(mx)), iw(2 * z(pn.nu));
    std::vector<double> es(z(mx)), wt(z(mx)), work(z(pn.nt + 1) * z(pn.nt + 1) + z(pn.nt) + 1);
    const int r = pev::bacon(yw.data(), pn.nu, pn.nt, g.data(), ty.data(), tr.data(), cg.data(),
                             es.data(), wt.data(), work.data(), iw.data());
    static const char *const tyn[] = {"treated vs never", "earlier vs later", "later vs earlier"};
    cols c(5);
    SEXP ct = c.str("type", r);
    double *cr = c.dbl("treated", r), *cc = c.dbl("control", r), *ce = c.dbl("est", r),
           *cw = c.dbl("weight", r);
    for (int q = 0; q < r; ++q) {
      SET_STRING_ELT(ct, q, Rf_mkChar(tyn[ty[z(q)]]));
      cr[q] = pn.tv[z(tr[z(q)])];
      cc[q] = cg[z(q)] < pn.nt ? pn.tv[z(cg[z(q)])] : R_PosInf;
      ce[q] = es[z(q)];
      cw[q] = wt[z(q)];
    }
    bc = PROTECT(c.frame());
    ++np;
  }
  const double cm = ctlmean(pn, d, m.w);
  const double re = fm == 0 ? f.b[0] / cm : std::exp(f.b[0]) - 1.0;
  cols fs(9);
  double *fn = fs.dbl("n", 1), *fu = fs.dbl("units", 1), *ft = fs.dbl("periods", 1),
         *fg = fs.dbl("ncl", 1), *fk = fs.dbl("k", 1), *fc = fs.dbl("ctlmean", 1),
         *fr = fs.dbl("releff", 1), *fth = fs.dbl("theta", 1), *fd = fs.dbl("dropped", 1);
  fn[0] = f.n;
  fu[0] = f.nu;
  ft[0] = f.nt;
  fg[0] = f.g;
  fk[0] = f.kk;
  fc[0] = cm;
  fr[0] = pevr::na(re);
  fth[0] = pevr::na(f.theta);
  fd[0] = f.drop;
  SEXP info = PROTECT(fs.frame());
  SEXP vc = PROTECT(Rf_allocMatrix(REALSXP, p, p));
  np += 2;
  for (int q = 0; q < p * p; ++q) REAL(vc)[q] = f.v[z(q)];
  SEXP fam1 = PROTECT(Rf_mkString(fmo[fm])), lvs = PROTECT(Rf_ScalarReal(lv));
  np += 2;
  SEXP out = PROTECT(nlist({{"coef", co}, {"fit", info}, {"vcov", vc}, {"loo", lo},
                            {"bacon", bc}, {"match", mt}, {"fam", fam1}, {"level", lvs}}));
  ++np;
  Rf_classgets(out, PROTECT(Rf_mkString("didfe")));
  UNPROTECT(np + 1);
  return out;
}

SEXP C_clab_r_evtst(SEXP y, SEXP unit, SEXP time, SEXP first, SEXP leads, SEXP lags,
                    SEXP ref, SEXP bin, SEXP adj, SEXP clus, SEXP w, SEXP ma, SEXP norm,
                    SEXP level) {
  clab::who = "evtst";
  const int ld = clab::whole(leads, "leads", 1), lg = clab::whole(lags, "lags", 0);
  const double rf = clab::val1(ref);
  if (!(rf == std::floor(rf) && rf >= -ld && rf <= lg)) {
    fail("`ref` must be one whole number from -leads to lags.");
  }
  const int r0 = static_cast<int>(rf);
  const bool bn = csp::flag(bin, "bin"), nm = csp::flag(norm, "norm");
  const int mw = clab::whole(ma, "ma", 1);
  const double lv = clab::one(level, "level", 0.0, 1.0);
  int np = 0;
  pevr::pan pn = pevr::panel(y, unit, time, clab::ninf, &np);
  const int n = pn.n;
  for (double v : pn.tv) {
    if (v != std::floor(v)) fail("`time` must be whole numbers (years, or months counted from a start).");
  }
  std::vector<double> gv = pevr::firstcol(first, pn);
  bool any = false;
  for (double v : gv) any = any || std::isfinite(v);
  if (!any) fail("no unit is ever treated: `first` is NA, 0 or Inf throughout.");
  std::vector<double> wv;
  if (!Rf_isNull(w)) {
    SEXP wx = PROTECT(clab::num(w, "w", n, 0.0));
    ++np;
    wv.assign(REAL(wx), REAL(wx) + n);
  }
  std::vector<int> c1, c2;
  int ncl = 0, ncl2 = 0;
  const int kc = clusters(clus, pn, c1, c2, &ncl, &ncl2);
  // event time of each row, binned at the ends or out of the window
  std::vector<double> k(z(n));
  std::vector<char> keep(z(n), 1);
  for (int i = 0; i < n; ++i) {
    const double gi = gv[z(pn.u[z(i)])];
    k[z(i)] = std::isfinite(gi) ? pn.tv[z(pn.t[z(i)])] - gi : R_PosInf;
    if (!std::isfinite(k[z(i)])) continue;
    if (bn) {
      k[z(i)] = std::min(std::max(k[z(i)], static_cast<double>(-ld)), static_cast<double>(lg));
    } else if (k[z(i)] < -ld || k[z(i)] > lg) {
      keep[z(i)] = 0;
    }
  }
  for (int i = 0; i < n; ++i) {
    if (!wv.empty() && !(wv[z(i)] > 0)) keep[z(i)] = 0;
  }
  std::vector<int> ks;
  for (int e = -ld; e <= lg; ++e) {
    if (e != r0) ks.push_back(e);
  }
  std::vector<pevr::term> x;
  for (int e : ks) {
    pevr::term tm;
    tm.nm = std::to_string(e);
    tm.x.resize(z(n));
    bool seen = false;
    for (int i = 0; i < n; ++i) {
      tm.x[z(i)] = k[z(i)] == e ? 1.0 : 0.0;
      seen = seen || (keep[z(i)] && k[z(i)] == e);
    }
    if (!seen) fail("no row has event time %d: narrow `leads` or `lags`.", e);
    x.push_back(tm);
  }
  const int q = static_cast<int>(x.size());
  pevr::adjcols(adj, n, x);
  const int p = static_cast<int>(x.size());
  model m{pn.y, &pn, &x, c1.data(), kc == 2 ? c2.data() : nullptr,
          wv.empty() ? nullptr : wv.data(), nullptr, 1, 0};
  std::vector<int> idx;
  for (int i = 0; i < n; ++i) {
    if (keep[z(i)]) idx.push_back(i);
  }
  const fit f = run(m, idx);
  if (!f.ok) fail("%s", absorbed);
  // the table with the reference period as a row of zeros
  const int ne = q + 1;
  const double df = f.g - 1, tq = Rf_qt(0.5 + 0.5 * lv, df, 1, 0);
  cols c(7);
  double *ev = c.dbl("evtm", ne), *e = c.dbl("est", ne), *s = c.dbl("se", ne),
         *lo = c.dbl("lo", ne), *hi = c.dbl("hi", ne), *st = c.dbl("stat", ne),
         *pv = c.dbl("p", ne);
  for (int a = 0, j = 0; a < ne; ++a) {
    const int et = -ld + a;
    ev[a] = et;
    if (et == r0) {
      e[a] = 0.0;
      s[a] = lo[a] = hi[a] = st[a] = pv[a] = NA_REAL;
      continue;
    }
    const double b = f.b[z(j)], vv = f.v[z(j) * z(p) + z(j)];
    const double se = vv >= 0 ? std::sqrt(vv) : qnan;
    e[a] = b;
    s[a] = pevr::na(se);
    lo[a] = pevr::na(b - tq * se);
    hi[a] = pevr::na(b + tq * se);
    st[a] = pevr::na(b / se);
    pv[a] = pevr::na(2 * Rf_pt(-std::fabs(b / se), df, 1, 0));
    ++j;
  }
  SEXP tab = PROTECT(c.frame());
  ++np;
  // the joint test of the leads: F = W / q on q and G - 1 degrees of freedom
  std::vector<int> pre;
  for (int j = 0; j < q; ++j) {
    if (ks[z(j)] < 0) pre.push_back(j);
  }
  const int nq = static_cast<int>(pre.size());
  double fst = NA_REAL, fp = NA_REAL;
  if (nq > 0) {
    std::vector<double> vp(z(nq) * z(nq)), vi(z(nq) * z(nq));
    for (int a = 0; a < nq; ++a) {
      for (int b = 0; b < nq; ++b) vp[z(b) * z(nq) + z(a)] = f.v[z(pre[z(b)]) * z(p) + z(pre[z(a)])];
    }
    if (pev::chol(vp.data(), nq)) {
      pev::cinv(vp.data(), nq, vi.data());
      double wd = 0.0;
      for (int a = 0; a < nq; ++a) {
        for (int b = 0; b < nq; ++b) {
          wd += f.b[z(pre[z(a)])] * vi[z(b) * z(nq) + z(a)] * f.b[z(pre[z(b)])];
        }
      }
      fst = wd / nq;
      fp = Rf_pf(fst, nq, df, 0, 0);
    }
  }
  cols pt(4);
  double *pf = pt.dbl("f", 1), *p1 = pt.dbl("df1", 1), *p2 = pt.dbl("df2", 1), *pp = pt.dbl("p", 1);
  pf[0] = fst;
  p1[0] = nq;
  p2[0] = df;
  pp[0] = fp;
  SEXP pret = PROTECT(pt.frame());
  ++np;
  SEXP vc = PROTECT(Rf_allocMatrix(REALSXP, q, q));
  ++np;
  for (int a = 0; a < q; ++a) {
    for (int b = 0; b < q; ++b) REAL(vc)[b * q + a] = f.v[z(b) * z(p) + z(a)];
  }
  // raw trends of each adoption cohort around its adoption
  std::vector<double> cv;
  for (double v : gv) {
    if (std::isfinite(v) && std::find(cv.begin(), cv.end(), v) == cv.end()) cv.push_back(v);
  }
  std::sort(cv.begin(), cv.end());
  std::vector<double> rc, re, rv;
  for (double cg : cv) {
    std::vector<double> sm(z(pn.nt), 0.0), cnt(z(pn.nt), 0.0);
    for (int i = 0; i < n; ++i) {
      if (gv[z(pn.u[z(i)])] == cg) {
        sm[z(pn.t[z(i)])] += pn.y[i];
        cnt[z(pn.t[z(i)])] += 1.0;
      }
    }
    std::vector<double> tt, mv;
    for (int a = 0; a < pn.nt; ++a) {
      if (cnt[z(a)] > 0) {
        tt.push_back(pn.tv[z(a)]);
        mv.push_back(sm[z(a)] / cnt[z(a)]);
      }
    }
    const int nc = static_cast<int>(tt.size());
    std::vector<double> av(z(nc), NA_REAL);
    for (int a = mw - 1; a < nc; ++a) {
      double s2 = 0.0;
      for (int b = a - mw + 1; b <= a; ++b) s2 += mv[z(b)];
      av[z(a)] = s2 / mw;
    }
    double base = NA_REAL;
    for (int a = 0; a < nc; ++a) {
      if (tt[z(a)] - cg == -1) base = av[z(a)];
    }
    for (int a = 0; a < nc; ++a) {
      const double et = tt[z(a)] - cg;
      if (et < -ld || et > lg) continue;
      rc.push_back(cg);
      re.push_back(et);
      double v = av[z(a)];
      if (nm) v = ISNAN(base) || base == 0 || ISNAN(v) ? NA_REAL : v / base;
      rv.push_back(v);
    }
  }
  const int nr = static_cast<int>(rc.size());
  cols ch(3);
  std::copy(rc.begin(), rc.end(), ch.dbl("cohort", nr));
  std::copy(re.begin(), re.end(), ch.dbl("evtm", nr));
  std::copy(rv.begin(), rv.end(), ch.dbl("value", nr));
  SEXP coh = PROTECT(ch.frame());
  SEXP nn = PROTECT(Rf_ScalarReal(f.n)), gg = PROTECT(Rf_ScalarReal(f.g)),
       lvs = PROTECT(Rf_ScalarReal(lv));
  np += 4;
  SEXP out = PROTECT(nlist({{"est", tab}, {"pre", pret}, {"vcov", vc}, {"coh", coh},
                            {"n", nn}, {"ncl", gg}, {"level", lvs}}));
  ++np;
  Rf_classgets(out, PROTECT(Rf_mkString("evtst")));
  UNPROTECT(np + 1);
  return out;
}

}  // extern "C"
