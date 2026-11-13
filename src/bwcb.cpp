// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of inference for a difference in differences with
// few clusters (wcbx: the wild cluster bootstrap, permutation of the
// adoption dates across clusters, and simulated power) and of the
// Westfall-Young stepdown (wystep). Each checks its arguments, runs the
// kernels (wcb.cpp, did.cpp) and builds the R result. Nothing calls back
// into the R evaluator.
#include <algorithm>
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "kutil.h"
#include "pev.h"
#include "pevr.h"
#include "pevw.h"
#include "rchk.h"
#include "spin.h"

extern "C" double Rf_pt(double, double, int, int);

using clab::fail;
using pevr::cols;
using pevr::z;

namespace {

const char *const absorbed =
    "the fit failed: a regressor is constant within units and periods (the fixed effects "
    "absorb it), or the regressors are collinear.";

// One two-way fit by clab_didfe: the first coefficient and its CR1 t
// statistic; false when the fit fails.
bool fit1(const std::vector<double> &y, const std::vector<double> &x, int p,
          const pevr::pan &pn, const std::vector<int> &cl, int ncl, double *b, double *t) {
  const int n = pn.n;
  std::vector<double> beta(z(p)), vc(z(p) * z(p)), st(5),
      work(z(n) * z(p + 2) + z(p) * z(p) + 2 * z(pn.nu + pn.nt) + z(n + 2 * p) * z(p) + z(p));
  std::vector<int> iwork(z(n) + z(std::max(pn.nu, pn.nt)));
  if (clab_didfe(y.data(), x.data(), n, p, pn.u.data(), pn.nu, pn.t.data(), pn.nt, 1, nullptr,
                 cl.data(), ncl, nullptr, 0, 1e-10, 10000, beta.data(), vc.data(), st.data(),
                 work.data(), iwork.data()) != CLAB_OK) {
    return false;
  }
  *b = beta[0];
  *t = beta[0] / std::sqrt(vc[0]);
  return true;
}

}  // namespace

extern "C" {

//' @srrstats {G2.6} wcbx() reads the outcome, units, periods, treatment
//'   and clusters as plain vectors of one length (several outcomes as the
//'   columns of a matrix or data frame).
//' @noRd
SEXP C_clab_r_wcbx(SEXP y, SEXP unit, SEXP time, SEXP trt, SEXP adj, SEXP clus, SEXP meth,
                   SEXP reps, SEXP wtype, SEXP ties, SEXP stat, SEXP ntrt, SEXP adopt, SEXP eff,
                   SEXP alpha, SEXP power, SEXP seed, SEXP level) {
  clab::who = "wcbx";
  static const char *const mo[] = {"wild", "perm", "power"};
  static const char *const wo[] = {"rad", "webb"};
  static const char *const to[] = {"count", "strict"};
  static const char *const so[] = {"t", "est"};
  const int mt = csp::pick(meth, "meth", mo, 3), wt = csp::pick(wtype, "wtype", wo, 2);
  const int tie = csp::pick(ties, "ties", to, 2) == 0 ? 1 : 0;
  const int sk = csp::pick(stat, "stat", so, 2);
  const int nr = Rf_isNull(reps) ? (mt == 2 ? 200 : 9999) : clab::whole(reps, "reps", 9);
  const double lv = clab::one(level, "level", 0.0, 1.0);
  int np = 0;
  // the outcomes: one vector, or the columns of a matrix or data frame
  std::vector<SEXP> ys;
  if (TYPEOF(y) == VECSXP) {
    for (R_xlen_t j = 0; j < XLENGTH(y); ++j) ys.push_back(VECTOR_ELT(y, j));
  } else if (Rf_isMatrix(y) && clab::isnum(y) && Rf_ncols(y) > 1) {
    const int nc = Rf_ncols(y), nrw = Rf_nrows(y);
    SEXP yd = PROTECT(clab::dbl(y));
    ++np;
    for (int j = 0; j < nc; ++j) {
      SEXP c = PROTECT(Rf_allocVector(REALSXP, nrw));
      ++np;
      std::copy(REAL(yd) + z(j) * z(nrw), REAL(yd) + z(j + 1) * z(nrw), REAL(c));
      ys.push_back(c);
    }
  } else {
    ys.push_back(y);
  }
  if (ys.empty()) fail("`y` has no outcome.");
  if (ys.size() > 1 && mt != 1) fail("several outcomes are taken by `meth = \"perm\"` only.");
  pevr::pan pn = pevr::panel(ys[0], unit, time, clab::ninf, &np);
  const int n = pn.n, k = static_cast<int>(ys.size());
  std::vector<std::vector<double>> yv(z(k));
  for (int j = 0; j < k; ++j) {
    SEXP v = PROTECT(clab::num(ys[z(j)], "y", n));
    ++np;
    yv[z(j)].assign(REAL(v), REAL(v) + n);
  }
  // clusters, one-way
  std::vector<int> cl(pn.u);
  int ncl = pn.nu;
  if (!Rf_isNull(clus)) {
    csp::lab l = csp::codes(clus, "clus", n);
    cl.assign(INTEGER(l.code), INTEGER(l.code) + n);
    ncl = l.g;
    UNPROTECT(2);
    if (ncl < 2) fail("`clus` has one cluster: the inference needs two or more.");
  }
  std::vector<pevr::term> x;
  pevr::adjcols(adj, n, x);
  const int pa = static_cast<int>(x.size());
  const int p = pa + 1;
  std::vector<double> xa(z(n) * z(p), 0.0);
  for (int j = 0; j < pa; ++j) std::copy(x[z(j)].x.begin(), x[z(j)].x.end(), xa.begin() + static_cast<std::ptrdiff_t>(j + 1) * n);
  SEXP res;
  int nd = 0;
  bool en = false;
  double mde = NA_REAL, size = NA_REAL;
  if (mt == 0 || mt == 1) {
    if (Rf_isNull(trt)) fail("`trt` is needed for the wild bootstrap and the permutations.");
    SEXP dv = PROTECT(csp::bin01(trt, "trt", n));
    ++np;
    for (int i = 0; i < n; ++i) xa[z(i)] = INTEGER(dv)[i];
    if (mt == 0) {
      // the regression within-transformed once; K by fixest's nested count
      pev::wctx c;
      c.n = n;
      c.p = p;
      c.x = xa;
      c.y = yv[0];
      c.g = cl;
      c.ng = ncl;
      c.u = pn.u;
      c.t = pn.t;
      c.nu = pn.nu;
      c.nt = pn.nt;
      std::vector<double> wk(2 * z(pn.nu + pn.nt));
      pev::within(c.x.data(), n, p, pn.u.data(), pn.nu, pn.t.data(), pn.nt, 1, nullptr, 1e-12,
                  10000, wk.data());
      pev::within(c.y.data(), n, 1, pn.u.data(), pn.nu, pn.t.data(), pn.nt, 1, nullptr, 1e-12,
                  10000, wk.data());
      auto nested = [&](const std::vector<int> &f, int nf) {
        std::vector<int> seen(z(nf), -1);
        for (int i = 0; i < n; ++i) {
          if (seen[z(f[z(i)])] < 0) seen[z(f[z(i)])] = cl[z(i)];
          if (seen[z(f[z(i)])] != cl[z(i)]) return false;
        }
        return true;
      };
      double kk = p + 1.0;
      if (!nested(pn.u, pn.nu)) kk += pn.nu - 1.0;
      if (!nested(pn.t, pn.nt)) kk += pn.nt - 1.0;
      c.adj = static_cast<double>(ncl) / (ncl - 1) * (n - 1.0) / (n - kk);
      if (!pev::wsetup(c)) fail("%s", absorbed);
      double bh = 0.0;
      const double th = pev::tstat(c, c.y.data(), 0.0, &bh);
      // the draws: every sign pattern when there are few enough
      en = wt == 0 && ncl <= 30 && std::ldexp(1.0, ncl) <= nr;
      std::vector<double> v;
      if (en) {
        nd = 1 << ncl;
        v.resize(z(nd) * z(ncl));
        for (int d = 0; d < nd; ++d) {
          for (int g = 0; g < ncl; ++g) v[z(d) * z(ncl) + z(g)] = (d >> g) & 1 ? 1.0 : -1.0;
        }
      } else {
        static const double webb[6] = {-1.224744871391589, -1.0, -0.7071067811865476,
                                       0.7071067811865476, 1.0, 1.224744871391589};
        kutil::rng g(csp::seed(seed));
        nd = nr;
        v.resize(z(nd) * z(ncl));
        for (double &e : v) e = wt == 0 ? (g.below(2) ? 1.0 : -1.0) : webb[g.below(6)];
      }
      const double pv = pev::wpval(c, 0.0, v, nd, tie);
      const double al = 1.0 - lv, s0 = std::fabs(bh / th);
      const bool open = tie && en && 2.0 / nd > al;
      const double lo = open ? R_NegInf : pev::wend(c, bh, s0, al, -1, v, nd, tie);
      const double hi = open ? R_PosInf : pev::wend(c, bh, s0, al, 1, v, nd, tie);
      cols r(5);
      double *re = r.dbl("est", 1), *rt = r.dbl("stat", 1), *rp = r.dbl("p", 1),
             *rl = r.dbl("lo", 1), *rh = r.dbl("hi", 1);
      re[0] = bh;
      rt[0] = th;
      rp[0] = pv;
      rl[0] = lo;
      rh[0] = hi;
      res = PROTECT(r.frame());
      ++np;
    } else {
      // permutation of the clusters' adoption periods
      const std::vector<int> g = pevr::firsts(pn, INTEGER(dv));
      std::vector<int> ga(z(ncl), -1);
      for (int i = 0; i < n; ++i) {
        const int a = g[z(pn.u[z(i)])], c = cl[z(i)];
        if (ga[z(c)] >= 0 && ga[z(c)] != a) fail("the units of a cluster must be treated from the same period.");
        ga[z(c)] = a;
      }
      std::vector<double> ob(z(k)), oe(z(k));
      for (int j = 0; j < k; ++j) {
        double b, t;
        if (!fit1(yv[z(j)], xa, p, pn, cl, ncl, &b, &t)) fail("%s", absorbed);
        oe[z(j)] = b;
        ob[z(j)] = sk == 0 ? t : b;
      }
      // every distinct assignment when there are at most reps of them
      std::vector<int> s = ga;
      std::sort(s.begin(), s.end());
      double lc = std::lgamma(ncl + 1.0);
      for (std::size_t a = 0; a < s.size();) {
        std::size_t b = a;
        while (b < s.size() && s[b] == s[a]) ++b;
        lc -= std::lgamma(static_cast<double>(b - a) + 1.0);
        a = b;
      }
      en = std::exp(lc) <= nr + 0.5;
      std::vector<std::vector<int>> as;
      if (en) {
        do as.push_back(s); while (std::next_permutation(s.begin(), s.end()));
      } else {
        kutil::rng rg(csp::seed(seed));
        for (int r = 0; r < nr; ++r) {
          std::vector<int> q = ga;
          for (int a = ncl - 1; a > 0; --a) std::swap(q[z(a)], q[z(rg.below(a + 1))]);
          as.push_back(q);
        }
      }
      nd = static_cast<int>(as.size());
      std::vector<double> nul(z(nd) * z(k), NA_REAL), xp = xa;
      std::vector<double> cnt(z(k), 0.0), ok(z(k), 0.0);
      for (int d = 0; d < nd; ++d) {
        for (int i = 0; i < n; ++i) xp[z(i)] = pn.t[z(i)] >= as[z(d)][z(cl[z(i)])] ? 1.0 : 0.0;
        for (int j = 0; j < k; ++j) {
          double b, t;
          if (!fit1(yv[z(j)], xp, p, pn, cl, ncl, &b, &t)) continue;
          const double sv = sk == 0 ? t : b;
          nul[z(j) * z(nd) + z(d)] = sv;
          ok[z(j)] += 1.0;
          cnt[z(j)] += tie ? std::fabs(sv) >= std::fabs(ob[z(j)]) - 1e-12
                           : std::fabs(sv) > std::fabs(ob[z(j)]) + 1e-12;
        }
      }
      // the stepdown over the draws where every outcome was estimated
      std::vector<double> full, pr(z(k)), pst(z(k), NA_REAL);
      int nf = 0;
      for (int d = 0; d < nd; ++d) {
        bool all = true;
        for (int j = 0; j < k; ++j) all = all && !ISNAN(nul[z(j) * z(nd) + z(d)]);
        nf += all;
      }
      full.resize(z(nf) * z(k));
      for (int d = 0, q = 0; d < nd; ++d) {
        bool all = true;
        for (int j = 0; j < k; ++j) all = all && !ISNAN(nul[z(j) * z(nd) + z(d)]);
        if (!all) continue;
        for (int j = 0; j < k; ++j) full[z(j) * z(nf) + z(q)] = nul[z(j) * z(nd) + z(d)];
        ++q;
      }
      if (k > 1 && nf > 0) pev::wystep(ob.data(), full.data(), nf, k, pr.data(), pst.data());
      cols r(5);
      SEXP ro = r.str("outcome", k);
      double *re = r.dbl("est", k), *rt = r.dbl("stat", k), *rp = r.dbl("p", k),
             *rs = r.dbl("pstep", k);
      SEXP yn = TYPEOF(y) == VECSXP ? Rf_getAttrib(y, R_NamesSymbol) : R_NilValue;
      for (int j = 0; j < k; ++j) {
        const bool nm = !Rf_isNull(yn) && STRING_ELT(yn, j) != NA_STRING;
        SET_STRING_ELT(ro, j, nm ? STRING_ELT(yn, j) : Rf_mkChar(("y" + std::to_string(j + 1)).c_str()));
        re[j] = oe[z(j)];
        rt[j] = ob[z(j)];
        rp[j] = en ? cnt[z(j)] / ok[z(j)] : (cnt[z(j)] + 1.0) / (ok[z(j)] + 1.0);
        rs[j] = k > 1 ? pst[z(j)] : rp[j];
      }
      res = PROTECT(r.frame());
      ++np;
    }
  } else {
    // power: ntrt clusters drawn, adoption periods drawn from adopt, the
    // effect added from adoption on, the CR1 t test at alpha
    const int nt1 = clab::whole(ntrt, "ntrt", 1);
    if (nt1 >= ncl) fail("`ntrt` must leave some clusters untreated: below %d.", ncl);
    if (Rf_isNull(adopt)) fail("`adopt` must give the adoption periods to draw from.");
    SEXP av = PROTECT(clab::num(adopt, "adopt"));
    ++np;
    std::vector<int> ac;
    for (R_xlen_t q = 0; q < XLENGTH(av); ++q) {
      const auto it = std::find(pn.tv.begin(), pn.tv.end(), REAL(av)[q]);
      if (it == pn.tv.end() || it == pn.tv.begin()) fail("`adopt` must be periods of `time` after the first.");
      ac.push_back(static_cast<int>(it - pn.tv.begin()));
    }
    const double al = clab::one(alpha, "alpha", 0.0, 1.0), pw = clab::one(power, "power", 0.0, 1.0);
    std::vector<double> ev;
    if (Rf_isNull(eff)) {
      const double sd = std::sqrt(std::max(0.0, [&] {
        double m = 0, s = 0;
        for (double v : yv[0]) m += v / n;
        for (double v : yv[0]) s += (v - m) * (v - m) / (n - 1);
        return s;
      }()));
      for (int q = 0; q < 5; ++q) ev.push_back(0.5 * q * sd);
    } else {
      SEXP e = PROTECT(clab::num(eff, "eff"));
      ++np;
      ev.assign(REAL(e), REAL(e) + XLENGTH(e));
      std::sort(ev.begin(), ev.end());
    }
    const int ne = static_cast<int>(ev.size());
    std::vector<double> rej(z(ne), 0.0), tot(z(ne), 0.0), yy(z(n)), xp = xa;
    kutil::rng rg(csp::seed(seed));
    for (int r = 0; r < nr; ++r) {
      std::vector<int> q(z(ncl)), ga(z(ncl), pn.nt);
      for (int a = 0; a < ncl; ++a) q[z(a)] = a;
      for (int a = 0; a < nt1; ++a) {
        std::swap(q[z(a)], q[z(a + rg.below(ncl - a))]);
        ga[z(q[z(a)])] = ac[z(rg.below(static_cast<int>(ac.size())))];
      }
      for (int i = 0; i < n; ++i) xp[z(i)] = pn.t[z(i)] >= ga[z(cl[z(i)])] ? 1.0 : 0.0;
      for (int e = 0; e < ne; ++e) {
        for (int i = 0; i < n; ++i) yy[z(i)] = yv[0][z(i)] + ev[z(e)] * xp[z(i)];
        double b, t;
        if (!fit1(yy, xp, p, pn, cl, ncl, &b, &t)) continue;
        tot[z(e)] += 1.0;
        rej[z(e)] += 2 * Rf_pt(-std::fabs(t), ncl - 1.0, 1, 0) < al ? 1.0 : 0.0;
      }
    }
    cols r(2);
    double *re = r.dbl("eff", ne), *rp = r.dbl("power", ne);
    for (int e = 0; e < ne; ++e) {
      re[e] = ev[z(e)];
      rp[e] = tot[z(e)] > 0 ? rej[z(e)] / tot[z(e)] : NA_REAL;
      if (ev[z(e)] == 0) size = rp[e];
    }
    for (int e = 1; e < ne && ISNAN(mde); ++e) {
      if (rp[e - 1] < pw && rp[e] >= pw) {
        mde = re[e - 1] + (pw - rp[e - 1]) / (rp[e] - rp[e - 1]) * (re[e] - re[e - 1]);
      }
    }
    if (ne > 0 && rp[0] >= pw) mde = re[0];
    nd = nr;
    res = PROTECT(r.frame());
    ++np;
  }
  SEXP m1 = PROTECT(Rf_mkString(mo[mt])), g1 = PROTECT(Rf_ScalarInteger(ncl)),
       d1 = PROTECT(Rf_ScalarInteger(nd)), e1 = PROTECT(Rf_ScalarLogical(en ? 1 : 0)),
       mm = PROTECT(Rf_ScalarReal(mde)), sz = PROTECT(Rf_ScalarReal(size)),
       l1 = PROTECT(Rf_ScalarReal(lv));
  np += 7;
  SEXP out = PROTECT(pevr::nlist({{"res", res}, {"meth", m1}, {"ncl", g1}, {"draws", d1},
                                  {"enum", e1}, {"mde", mm}, {"size", sz}, {"level", l1}}));
  ++np;
  Rf_classgets(out, PROTECT(Rf_mkString("wcbx")));
  UNPROTECT(np + 1);
  return out;
}

//' @srrstats {G2.6} wystep() reads the observed statistics as a vector and
//'   the null draws as a numeric matrix, one column per statistic.
//' @noRd
SEXP C_clab_r_wystep(SEXP obs, SEXP null) {
  clab::who = "wystep";
  clab::vec(obs, "obs");
  SEXP ob = PROTECT(clab::num(obs, "obs"));
  const int k = clab::ilen(static_cast<double>(XLENGTH(ob)), "obs");
  if (!Rf_isMatrix(null) || !clab::isnum(null) || Rf_ncols(null) != k) {
    fail("`null` must be a numeric matrix with one column per statistic of `obs`, %d.", k);
  }
  SEXP nv = PROTECT(clab::dbl(null));
  const int nr = Rf_nrows(null);
  std::vector<double> full;
  int nd = 0;
  for (int d = 0; d < nr; ++d) {
    bool ok = true;
    for (int j = 0; j < k; ++j) ok = ok && std::isfinite(REAL(nv)[z(j) * z(nr) + z(d)]);
    nd += ok;
  }
  if (nd == 0) fail("`null` has no draw with every statistic finite.");
  full.resize(z(nd) * z(k));
  for (int d = 0, q = 0; d < nr; ++d) {
    bool ok = true;
    for (int j = 0; j < k; ++j) ok = ok && std::isfinite(REAL(nv)[z(j) * z(nr) + z(d)]);
    if (!ok) continue;
    for (int j = 0; j < k; ++j) full[z(j) * z(nd) + z(q)] = REAL(nv)[z(j) * z(nr) + z(d)];
    ++q;
  }
  cols c(3);
  double *cs = c.dbl("stat", k), *cp = c.dbl("p", k), *ca = c.dbl("padj", k);
  std::copy(REAL(ob), REAL(ob) + k, cs);
  pev::wystep(REAL(ob), full.data(), nd, k, cp, ca);
  SEXP out = c.frame();
  UNPROTECT(2);
  return out;
}

}  // extern "C"
