// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Panels of units over periods (src/pevr.cpp, src/bdid.cpp): the panel
// reader, two-way fixed effects difference in differences and the event
// study.

inline std::size_t zz(int i) { return static_cast<std::size_t>(i); }

// pevr::panel(): y, unit and time, one row per unit and period.
struct Pan {
  int n = 0;
  std::vector<double> y;
  std::vector<int> u, t;
  int nu = 0, nt = 0;
  std::vector<std::string> ulev;
  std::vector<double> tv;
};
bool isvec(Ctx &c, const Val &o, const char *nm) {
  if (o.t == Val::ARR) {
    for (const Val &e : o.a) {
      if (e.t == Val::ARR || e.t == Val::OBJ) return notvec(c, nm);
    }
  }
  return true;
}
bool panel(Ctx &c, double lo, Pan &p, const Val *y = nullptr) {
  if (!y) y = &c["y"];
  CK(isvec(c, *y, "y"));
  Dv yv;
  CK(num(c, *y, "y", yv, -1, lo));
  p.n = static_cast<int>(yv.n());
  if (p.n < 4) return c.fail("`y` must have at least four rows (two units in two periods).");
  p.y = yv.v;
  Lab l;
  CK(codes(c, c["unit"], "unit", yv.n(), l));
  p.u = l.code;
  p.nu = l.g;
  p.ulev = l.lev;
  CK(isvec(c, c["time"], "time"));
  Dv tv;
  CK(num(c, c["time"], "time", tv, p.n));
  p.tv = tv.v;
  std::sort(p.tv.begin(), p.tv.end());
  p.tv.erase(std::unique(p.tv.begin(), p.tv.end()), p.tv.end());
  p.nt = static_cast<int>(p.tv.size());
  p.t.resize(zz(p.n));
  for (std::size_t i = 0; i < zz(p.n); ++i) {
    p.t[i] = static_cast<int>(std::lower_bound(p.tv.begin(), p.tv.end(), tv[i]) - p.tv.begin());
  }
  if (p.nu < 2 || p.nt < 2) return c.fail("the panel needs two units and two periods at least.");
  std::vector<char> seen(zz(p.nu) * zz(p.nt), 0);
  for (std::size_t i = 0; i < zz(p.n); ++i) {
    char &s = seen[zz(p.u[i]) * zz(p.nt) + zz(p.t[i])];
    if (s) return c.fail("`unit` and `time` repeat a pair: one row per unit and period.");
    s = 1;
  }
  return true;
}

// pevr::firsts(): each unit's first treated period (nt for never).
bool firsts(Ctx &c, const Pan &pn, const int *d, std::vector<int> &g) {
  std::vector<int> dw(zz(pn.nu) * zz(pn.nt), -1);
  g.assign(zz(pn.nu), pn.nt);
  for (std::size_t i = 0; i < zz(pn.n); ++i) dw[zz(pn.u[i]) * zz(pn.nt) + zz(pn.t[i])] = d[i];
  for (std::size_t a = 0; a < zz(pn.nu); ++a) {
    for (int s = 0; s < pn.nt; ++s) {
      const int v = dw[a * zz(pn.nt) + zz(s)];
      if (v == 1 && g[a] == pn.nt) g[a] = s;
      if (v == 0 && g[a] < pn.nt) return c.fail("`trt` must stay 1 once a unit is treated (unit %s).", pn.ulev[a].c_str());
    }
  }
  return true;
}

// pevr::firstcol(): `first` per row, the period each unit is first
// treated, Inf for never (NA, 0 or Inf).
bool firstcol(Ctx &c, const Pan &pn, std::vector<double> &g) {
  const Val &f = c["first"];
  // all nulls is logical NA in R, not numbers
  bool any = f.t == Val::NUM;
  for (const Val &e : f.a) any = any || e.t != Val::NUL;
  if (vlen(f) != zz(pn.n) || !isnum(f) || !any) {
    return c.fail("`first` must be numbers, one per row: the period each unit is first treated.");
  }
  const std::vector<double> fv = dbl(f);
  g.assign(zz(pn.nu), 0.0);
  std::vector<char> seen(zz(pn.nu), 0);
  for (std::size_t i = 0; i < zz(pn.n); ++i) {
    double v = fv[i];
    if (std::isnan(v) || v == 0 || v == Inf) v = Inf;
    if (v == -Inf) return c.fail("`first` must be a period, or NA, 0 or Inf for never.");
    const std::size_t a = zz(pn.u[i]);
    if (seen[a] && !(g[a] == v)) return c.fail("`first` must be constant within each unit.");
    seen[a] = 1;
    g[a] = v;
  }
  return true;
}

// clusters(): the units, or one vector of labels, or a list of two.
int pclus(Ctx &c, const Pan &pn, std::vector<int> &c1, std::vector<int> &c2, int &ncl, int &ncl2) {
  ncl2 = 0;
  if (!c.has("clus")) {
    c1 = pn.u;
    ncl = pn.nu;
    return 0;
  }
  const Val &cl = c["clus"];
  const Val *a = &cl, *b = nullptr;
  if (cl.t == Val::OBJ || (cl.t == Val::ARR && !cl.a.empty() && cl.a[0].t == Val::ARR)) {
    if (cl.a.empty() || cl.a.size() > 2) {
      c.fail("`clus` must be one vector of cluster labels or a list of two.");
      return -1;
    }
    a = &cl.a[0];
    if (cl.a.size() == 2) b = &cl.a[1];
  }
  Lab l;
  if (!codes(c, *a, "clus", zz(pn.n), l)) return -1;
  c1 = l.code;
  ncl = l.g;
  if (b) {
    Lab m;
    if (!codes(c, *b, "clus", zz(pn.n), m)) return -1;
    c2 = m.code;
    ncl2 = m.g;
  }
  if (ncl < 2 || (b && ncl2 < 2)) {
    c.fail("`clus` has one cluster: a cluster-robust variance needs two or more.");
    return -1;
  }
  if (ncl < 10 || (b && ncl2 < 10)) {
    c.warn.push_back(std::string(c.fn) +
                     ": fewer than ten clusters; cluster-robust errors are unreliable with so few (wcbx() is made "
                     "for this case).");
  }
  return b ? 2 : 1;
}

struct Pfit {
  bool ok = false;
  std::vector<double> b, v;
  double g = 0, kk = 0, n = 0, theta = NaN;
  int nu = 0, nt = 0, drop = 0;
};
struct Pmod {
  const double *y;
  const Pan *pn;
  const std::vector<Term> *x;
  const int *cl, *cl2;
  const double *w, *off;
  int fe = 0, fam = 0;
};

int recode(const int *src, const std::vector<int> &idx, std::vector<int> &out) {
  std::vector<int> map;
  out.resize(idx.size());
  int k = 0;
  for (std::size_t q = 0; q < idx.size(); ++q) {
    const int s = src[idx[q]];
    if (s >= static_cast<int>(map.size())) map.resize(zz(s) + 1, -1);
    if (map[zz(s)] < 0) map[zz(s)] = k++;
    out[q] = map[zz(s)];
  }
  return k;
}
bool nests(const std::vector<int> &f, int nf, const std::vector<int> &cl) {
  std::vector<int> seen(zz(nf), -1);
  for (std::size_t i = 0; i < f.size(); ++i) {
    if (seen[zz(f[i])] < 0) seen[zz(f[i])] = cl[i];
    if (seen[zz(f[i])] != cl[i]) return false;
  }
  return true;
}

// run() of src/bdid.cpp: the fit on the rows idx.
Pfit prun(const Pmod &m, std::vector<int> idx) {
  Pfit r;
  const int p = static_cast<int>(m.x->size());
  if (m.fam > 0) {
    for (;;) {
      std::vector<double> su(zz(m.pn->nu), 0.0), st(zz(m.pn->nt), 0.0);
      for (int i : idx) {
        su[zz(m.pn->u[zz(i)])] += m.y[i];
        st[zz(m.pn->t[zz(i)])] += m.y[i];
      }
      std::vector<int> keep;
      for (int i : idx) {
        if ((su[zz(m.pn->u[zz(i)])] > 0 || !m.fe) && st[zz(m.pn->t[zz(i)])] > 0) keep.push_back(i);
      }
      if (keep.size() == idx.size()) break;
      idx.swap(keep);
    }
  }
  const int n = static_cast<int>(idx.size());
  std::vector<int> u, t, c1, c2;
  const int nu = recode(m.pn->u.data(), idx, u), nt = recode(m.pn->t.data(), idx, t);
  const int ncl = recode(m.cl, idx, c1), ncl2 = m.cl2 ? recode(m.cl2, idx, c2) : 0;
  r.nu = nu;
  r.nt = nt;
  r.n = n;
  if (n < 2 || ncl < 2 || (m.cl2 && ncl2 < 2)) return r;
  std::vector<double> y(zz(n)), w, x(zz(n) * zz(p));
  for (int q = 0; q < n; ++q) y[zz(q)] = m.y[idx[zz(q)]];
  if (m.w) {
    w.resize(zz(n));
    for (int q = 0; q < n; ++q) w[zz(q)] = m.w[idx[zz(q)]];
  }
  for (int j = 0; j < p; ++j) {
    for (int q = 0; q < n; ++q) x[zz(j) * zz(n) + zz(q)] = (*m.x)[zz(j)].x[zz(idx[zz(q)])];
  }
  r.b.assign(zz(p), NaN);
  r.v.assign(zz(p) * zz(p), NaN);
  if (m.fam == 0) {
    const std::size_t sp = zz(p), sn = zz(n);
    std::vector<double> work(sn * (sp + 2) + sp * sp + 2 * zz(nu + nt) + (sn + 2 * sp) * sp + sp);
    std::vector<int> iwork(sn + std::max<std::size_t>({zz(nu), zz(nt), zz(ncl) * zz(ncl2)}));
    double st[5] = {0};
    if (clab_didfe(y.data(), x.data(), n, p, u.data(), nu, t.data(), nt, m.fe, m.w ? w.data() : nullptr, c1.data(), ncl,
                   m.cl2 ? c2.data() : nullptr, ncl2, 1e-10, 10000, r.b.data(), r.v.data(), st, work.data(),
                   iwork.data()) != CLAB_OK) {
      return r;
    }
    r.g = st[1];
    r.kk = st[2];
    r.ok = true;
    return r;
  }
  const int pd = 1 + p + (m.fe ? nu - 1 : 0) + nt - 1;
  std::vector<double> xd(zz(n) * zz(pd), 0.0), off;
  for (int q = 0; q < n; ++q) xd[zz(q)] = 1.0;
  std::copy(x.begin(), x.end(), xd.begin() + n);
  int cc = 1 + p;
  if (m.fe) {
    for (int a = 1; a < nu; ++a, ++cc) {
      for (int q = 0; q < n; ++q) xd[zz(cc) * zz(n) + zz(q)] = u[zz(q)] == a ? 1.0 : 0.0;
    }
  }
  for (int a = 1; a < nt; ++a, ++cc) {
    for (int q = 0; q < n; ++q) xd[zz(cc) * zz(n) + zz(q)] = t[zz(q)] == a ? 1.0 : 0.0;
  }
  if (m.off) {
    off.resize(zz(n));
    for (int q = 0; q < n; ++q) off[zz(q)] = m.off[idx[zz(q)]];
  }
  std::vector<int> c12(zz(n)), pw(zz(ncl) * zz(std::max(ncl2, 1)));
  const int n12 = m.cl2 ? pev::pairs(c1.data(), c2.data(), n, ncl, ncl2, c12.data(), pw.data()) : 0;
  const int big = std::max({ncl, ncl2, n12});
  const std::size_t sp = zz(pd), sn = zz(n);
  std::vector<double> work(sn * sp + 8 * sn + sp + sp * sp + 2 * zz(big) * sp), beta(sp), vb(sp * sp), vc(sp * sp),
      vs(sp * sp, 0.0), st(8);
  std::vector<int> alias(sp), iwork(sp + zz(big) + 1 + sn);
  const int *cls[3] = {c1.data(), c2.data(), c12.data()};
  const int ncs[3] = {ncl, ncl2, n12};
  const double sg[3] = {1.0, 1.0, -1.0};
  for (int k = 0; k < (m.cl2 ? 3 : 1); ++k) {
    if (clab_cntreg(y.data(), xd.data(), n, pd, m.off ? off.data() : nullptr, cls[k], ncs[k], m.fam == 2 ? 1 : 0, 25,
                    1e-8, beta.data(), vb.data(), vc.data(), st.data(), alias.data(), work.data(),
                    iwork.data()) != CLAB_OK) {
      return r;
    }
    for (int j = 1; j <= p; ++j) {
      if (alias[zz(j)]) return r;
    }
    for (std::size_t q = 0; q < sp * sp; ++q) vs[q] += sg[k] * vc[q];
  }
  r.g = m.cl2 ? std::min(ncl, ncl2) : ncl;
  r.kk = p + 1.0;
  if (m.fe && !(nests(u, nu, c1) || (m.cl2 && nests(u, nu, c2)))) r.kk += nu - 1.0;
  if (!(nests(t, nt, c1) || (m.cl2 && nests(t, nt, c2)))) r.kk += nt - 1.0;
  const double cf = r.g / (r.g - 1) * (n - 1.0) / (n - r.kk);
  for (int j = 0; j < p; ++j) {
    r.b[zz(j)] = beta[zz(j) + 1];
    for (int k = 0; k < p; ++k) r.v[zz(k) * zz(p) + zz(j)] = cf * vs[zz(k + 1) * sp + zz(j) + 1];
  }
  r.theta = m.fam == 2 ? st[2] : NaN;
  r.drop = m.pn->nu - nu;
  r.ok = true;
  return r;
}

// ctab(): the coefficients p0 to p - 1 of a fit.
Val pctab(const Pfit &f, const std::vector<Term> &x, int p0, int p, bool cnt, double lv) {
  const double df = f.g - 1;
  const double q = cnt ? clab_qnorm(0.5 + 0.5 * lv) : clab_qt(0.5 + 0.5 * lv, df, 1);
  std::vector<std::string> tm;
  std::vector<double> e, s, lo, hi, st, pv, rr, rl, rh;
  for (int a = p0; a < p; ++a) {
    tm.push_back(x[zz(a)].nm);
    const double b = f.b[zz(a)], v = f.v[zz(a) * zz(p) + zz(a)];
    const double se = v >= 0 ? std::sqrt(v) : NaN;
    e.push_back(b);
    s.push_back(na(se));
    lo.push_back(na(b - q * se));
    hi.push_back(na(b + q * se));
    st.push_back(na(b / se));
    pv.push_back(na(cnt ? 2 * clab_pnorm(-std::fabs(b / se)) : 2 * clab_pt(-std::fabs(b / se), df, 1)));
    rr.push_back(std::exp(b));
    rl.push_back(na(std::exp(b - q * se)));
    rh.push_back(na(std::exp(b + q * se)));
  }
  Frame fr;
  fr.str("term", tm);
  fr.num("est", e);
  fr.num("se", s);
  fr.num("lo", lo);
  fr.num("hi", hi);
  fr.num("stat", st);
  fr.num("p", pv);
  if (cnt) {
    fr.num("rr", rr);
    fr.num("rrlo", rl);
    fr.num("rrhi", rh);
  }
  return fr.o;
}

double ctlmean(const Pan &pn, const int *d, const double *w) {
  std::vector<char> ever(zz(pn.nu), 0);
  int first = pn.nt;
  for (std::size_t i = 0; i < zz(pn.n); ++i) {
    if (d[i]) {
      ever[zz(pn.u[i])] = 1;
      first = std::min(first, pn.t[i]);
    }
  }
  for (int pass = 0; pass < 2; ++pass) {
    double s = 0.0, sw = 0.0;
    for (std::size_t i = 0; i < zz(pn.n); ++i) {
      const bool in = pn.t[i] < first && (pass == 0 ? !ever[zz(pn.u[i])] : d[i] == 0);
      if (!in) continue;
      const double wi = w ? w[i] : 1.0;
      s += wi * pn.y[i];
      sw += wi;
    }
    if (sw > 0) return s / sw;
  }
  return NA;
}

const char *const absorbed =
    "the fit failed: a regressor is constant within units and periods (the fixed effects absorb it), or the "
    "regressors are collinear.";

// A p x p matrix as rows.
Val pmat(const std::vector<double> &v, int p) { return vmat(v, p, p); }

bool didfe(Ctx &c, Val &out) {
  static const char *const feo[] = {"both", "time"};
  static const char *const fmo[] = {"ols", "pois", "nb"};
  int fe0 = 0, fm = 0, ratio = 0;
  CK(pickd(c, "fe", feo, 2, 0, fe0));
  const int fx = fe0 == 0 ? 1 : 0;
  CK(pickd(c, "fam", fmo, 3, 0, fm));
  double lv = 0.0;
  CK(oned(c, "level", 0.95, 0.0, 1.0, lv));
  bool dloo = false, dbac = false;
  CK(flagd(c, "loo", false, dloo));
  CK(flagd(c, "bacon", false, dbac));
  CK(wholed(c, "mtch", 0, 0, ratio));
  Pan pn;
  CK(panel(c, fm > 0 ? 0.0 : -HUGE_VAL, pn));
  const int n = pn.n;
  std::vector<int> dv;
  CK(bin01(c, c["trt"], "trt", zz(n), dv));
  const int *d = dv.data();
  int nd = 0;
  for (int e : dv) nd += e;
  if (nd == 0 || nd == n) return c.fail("`trt` must be 1 in some rows and 0 in others.");
  std::vector<double> wv;
  if (c.has("w")) {
    Dv wx;
    CK(num(c, c["w"], "w", wx, n, 0.0));
    wv = wx.v;
  }
  std::vector<double> off;
  if (c.has("expo")) {
    if (fm == 0) return c.fail("`expo` is an exposure for counts: give it with `fam = \"pois\"` or \"nb\".");
    Dv e;
    CK(num(c, c["expo"], "expo", e, n, 0.0, true));
    for (double v : e.v) off.push_back(std::log(v));
  }
  std::vector<int> c1, c2;
  int ncl = 0, ncl2 = 0;
  const int kc = pclus(c, pn, c1, c2, ncl, ncl2);
  if (kc < 0) return false;
  std::vector<Term> x(1);
  x[0].nm = "trt";
  x[0].x.assign(d, d + n);
  CK(adjcols(c, c["adj"], zz(n), x, "row of the panel"));
  const int p = static_cast<int>(x.size());
  if ((dbac || ratio > 0) && fm != 0) return c.fail("`bacon` and `mtch` work with `fam = \"ols\"`.");
  if (dbac && (p > 1 || !wv.empty() || !fx)) {
    return c.fail("`bacon = TRUE` decomposes the plain two-way estimate: no `adj`, no `w`, `fe = \"both\"`.");
  }
  if (c.has("exact") && ratio == 0) return c.fail("`exact` groups the matching: give `mtch` too.");
  if (c.has("calip") && ratio == 0) return c.fail("`calip` limits the matching: give `mtch` too.");
  std::vector<int> g;
  if (dbac || ratio > 0) {
    if (n != pn.nu * pn.nt) return c.fail("`bacon` and `mtch` need a balanced panel: every unit in every period.");
    CK(firsts(c, pn, d, g));
  }
  Pmod m{pn.y.data(), &pn, &x, c1.data(), kc == 2 ? c2.data() : nullptr, wv.empty() ? nullptr : wv.data(),
         off.empty() ? nullptr : off.data(), fx, fm};
  Val mt;
  std::vector<double> mw;
  if (ratio > 0) {
    std::vector<int> tr(zz(pn.nu));
    int g0 = pn.nt;
    for (std::size_t a = 0; a < zz(pn.nu); ++a) {
      tr[a] = g[a] < pn.nt ? 1 : 0;
      if (tr[a]) g0 = std::min(g0, g[a]);
    }
    if (std::count(tr.begin(), tr.end(), 1) == pn.nu) return c.fail("`mtch` needs units never treated to match to.");
    if (g0 == 0) return c.fail("`mtch` needs a period before the first adoption to match on.");
    std::vector<int> ex(zz(pn.nu), 0);
    if (c.has("exact")) {
      Lab l;
      CK(codes(c, c["exact"], "exact", zz(n), l));
      std::vector<int> seen(zz(pn.nu), -1);
      for (std::size_t i = 0; i < zz(n); ++i) {
        const std::size_t a = zz(pn.u[i]);
        const int e = l.code[i];
        if (seen[a] >= 0 && seen[a] != e) return c.fail("`exact` must be constant within each unit.");
        seen[a] = e;
        ex[a] = e;
      }
    }
    double cal = 0.0;
    if (c.has("calip")) CK(one(c, c["calip"], "calip", 0.0, Inf, cal));
    std::vector<double> f(zz(pn.nu) * zz(g0));
    for (std::size_t i = 0; i < zz(n); ++i) {
      if (pn.t[i] < g0) f[zz(pn.t[i]) * zz(pn.nu) + zz(pn.u[i])] = pn.y[i];
    }
    const int n1 = static_cast<int>(std::count(tr.begin(), tr.end(), 1));
    std::vector<int> ma(zz(n1) * zz(ratio)), mc(zz(n1) * zz(ratio));
    std::vector<double> md(zz(n1) * zz(ratio)), wu(zz(pn.nu)), work(zz(g0) * zz(g0) + 3 * zz(g0));
    const int nm = pev::mtch(f.data(), pn.nu, g0, tr.data(), ex.data(), ratio, cal, ma.data(), mc.data(), md.data(),
                             wu.data(), work.data());
    if (nm < 0) {
      return c.fail("the pre-period outcomes are collinear across units: too many periods before the first adoption "
                    "for the units to match on.");
    }
    if (nm == 0) return c.fail("no treated unit has a comparison unit within its `exact` group and the caliper.");
    std::vector<std::string> ct, cc;
    for (int q = 0; q < nm; ++q) {
      ct.push_back(pn.ulev[zz(ma[zz(q)])]);
      cc.push_back(pn.ulev[zz(mc[zz(q)])]);
    }
    Frame fr;
    fr.str("treated", ct);
    fr.str("control", cc);
    fr.num("dist", md.data(), zz(nm));
    mt = fr.o;
    mw.resize(zz(n));
    for (std::size_t i = 0; i < zz(n); ++i) mw[i] = wu[zz(pn.u[i])] * (wv.empty() ? 1.0 : wv[i]);
    m.w = mw.data();
  }
  std::vector<int> all;
  for (int i = 0; i < n; ++i) {
    if (!m.w || m.w[i] > 0) all.push_back(i);
  }
  const Pfit f = prun(m, all);
  if (!f.ok) return c.fail("%s", absorbed);
  out = vobj();
  put(out, "coef", pctab(f, x, 0, p, fm > 0, lv));
  Val lo;
  if (dloo) {
    std::vector<char> ever(zz(pn.nu), 0);
    for (int i : all) {
      if (d[i]) ever[zz(pn.u[zz(i)])] = 1;
    }
    std::vector<std::string> cu;
    std::vector<double> ce, cs;
    for (int a = 0; a < pn.nu; ++a) {
      if (!ever[zz(a)]) continue;
      std::vector<int> idx;
      for (int i : all) {
        if (pn.u[zz(i)] != a) idx.push_back(i);
      }
      const Pfit h = prun(m, idx);
      cu.push_back(pn.ulev[zz(a)]);
      ce.push_back(h.ok ? h.b[0] : NA);
      cs.push_back(h.ok ? na(std::sqrt(h.v[0])) : NA);
    }
    Frame fr;
    fr.str("unit", cu);
    fr.num("est", ce);
    fr.num("se", cs);
    lo = fr.o;
  }
  Val bc;
  if (dbac) {
    std::vector<double> yw(zz(n));
    for (std::size_t i = 0; i < zz(n); ++i) yw[zz(pn.t[i]) * zz(pn.nu) + zz(pn.u[i])] = pn.y[i];
    const int mx = pn.nt * pn.nt + pn.nt;
    std::vector<int> ty(zz(mx)), tr(zz(mx)), cg(zz(mx)), iw(2 * zz(pn.nu));
    std::vector<double> es(zz(mx)), wt(zz(mx)), work(zz(pn.nt + 1) * zz(pn.nt + 1) + zz(pn.nt) + 1);
    const int r = pev::bacon(yw.data(), pn.nu, pn.nt, g.data(), ty.data(), tr.data(), cg.data(), es.data(), wt.data(),
                             work.data(), iw.data());
    static const char *const tyn[] = {"treated vs never", "earlier vs later", "later vs earlier"};
    std::vector<std::string> ct;
    std::vector<double> cr, cc;
    for (int q = 0; q < r; ++q) {
      ct.push_back(tyn[ty[zz(q)]]);
      cr.push_back(pn.tv[zz(tr[zz(q)])]);
      cc.push_back(cg[zz(q)] < pn.nt ? pn.tv[zz(cg[zz(q)])] : Inf);
    }
    Frame fr;
    fr.str("type", ct);
    fr.num("treated", cr);
    fr.num("control", cc);
    fr.num("est", es.data(), zz(r));
    fr.num("weight", wt.data(), zz(r));
    bc = fr.o;
  }
  const double cm = ctlmean(pn, d, m.w);
  const double re = fm == 0 ? f.b[0] / cm : std::exp(f.b[0]) - 1.0;
  Frame fs;
  col1(fs, "n", f.n);
  col1(fs, "units", f.nu);
  col1(fs, "periods", f.nt);
  col1(fs, "ncl", f.g);
  col1(fs, "k", f.kk);
  {
    // the control mean as R keeps it: NA stays NA
    fs.num("ctlmean", &cm, 1);
  }
  col1(fs, "releff", re);
  col1(fs, "theta", f.theta);
  col1(fs, "dropped", f.drop);
  put(out, "fit", fs.o);
  put(out, "vcov", pmat(f.v, p));
  put(out, "loo", lo);
  put(out, "bacon", bc);
  put(out, "match", mt);
  put(out, "fam", vstr(fmo[fm]));
  put(out, "level", vnum(lv));
  return true;
}

bool evtst(Ctx &c, Val &out) {
  int ld = 0, lg = 0, mw = 0;
  CK(wholed(c, "leads", 4, 1, ld));
  CK(wholed(c, "lags", 4, 0, lg));
  const double rf = val1d(c, "ref", -1);
  if (!(rf == std::floor(rf) && rf >= -ld && rf <= lg)) return c.fail("`ref` must be one whole number from -leads to lags.");
  const int r0 = static_cast<int>(rf);
  bool bn = true, nm = true;
  CK(flagd(c, "bin", true, bn));
  CK(flagd(c, "norm", true, nm));
  CK(wholed(c, "ma", 1, 1, mw));
  double lv = 0.0;
  CK(oned(c, "level", 0.95, 0.0, 1.0, lv));
  Pan pn;
  CK(panel(c, -HUGE_VAL, pn));
  const int n = pn.n;
  for (double v : pn.tv) {
    if (v != std::floor(v)) return c.fail("`time` must be whole numbers (years, or months counted from a start).");
  }
  std::vector<double> gv;
  CK(firstcol(c, pn, gv));
  bool any = false;
  for (double v : gv) any = any || std::isfinite(v);
  if (!any) return c.fail("no unit is ever treated: `first` is NA, 0 or Inf throughout.");
  std::vector<double> wv;
  if (c.has("w")) {
    Dv wx;
    CK(num(c, c["w"], "w", wx, n, 0.0));
    wv = wx.v;
  }
  std::vector<int> c1, c2;
  int ncl = 0, ncl2 = 0;
  const int kc = pclus(c, pn, c1, c2, ncl, ncl2);
  if (kc < 0) return false;
  std::vector<double> k(zz(n));
  std::vector<char> keep(zz(n), 1);
  for (std::size_t i = 0; i < zz(n); ++i) {
    const double gi = gv[zz(pn.u[i])];
    k[i] = std::isfinite(gi) ? pn.tv[zz(pn.t[i])] - gi : Inf;
    if (!std::isfinite(k[i])) continue;
    if (bn) k[i] = std::min(std::max(k[i], static_cast<double>(-ld)), static_cast<double>(lg));
    else if (k[i] < -ld || k[i] > lg) keep[i] = 0;
  }
  for (std::size_t i = 0; i < zz(n); ++i) {
    if (!wv.empty() && !(wv[i] > 0)) keep[i] = 0;
  }
  std::vector<int> ks;
  for (int e = -ld; e <= lg; ++e) {
    if (e != r0) ks.push_back(e);
  }
  std::vector<Term> x;
  for (int e : ks) {
    Term tm;
    tm.nm = std::to_string(e);
    tm.x.resize(zz(n));
    bool seen = false;
    for (std::size_t i = 0; i < zz(n); ++i) {
      tm.x[i] = k[i] == e ? 1.0 : 0.0;
      seen = seen || (keep[i] && k[i] == e);
    }
    if (!seen) return c.fail("no row has event time %d: narrow `leads` or `lags`.", e);
    x.push_back(tm);
  }
  const int q = static_cast<int>(x.size());
  CK(adjcols(c, c["adj"], zz(n), x, "row of the panel"));
  const int p = static_cast<int>(x.size());
  Pmod m{pn.y.data(), &pn, &x, c1.data(), kc == 2 ? c2.data() : nullptr, wv.empty() ? nullptr : wv.data(), nullptr, 1, 0};
  std::vector<int> idx;
  for (int i = 0; i < n; ++i) {
    if (keep[zz(i)]) idx.push_back(i);
  }
  const Pfit f = prun(m, idx);
  if (!f.ok) return c.fail("%s", absorbed);
  const int ne = q + 1;
  const double df = f.g - 1, tq = clab_qt(0.5 + 0.5 * lv, df, 1);
  std::vector<double> ev, e, s, lo, hi, st, pv;
  for (int a = 0, j = 0; a < ne; ++a) {
    const int et = -ld + a;
    ev.push_back(et);
    if (et == r0) {
      e.push_back(0.0);
      s.push_back(NA);
      lo.push_back(NA);
      hi.push_back(NA);
      st.push_back(NA);
      pv.push_back(NA);
      continue;
    }
    const double b = f.b[zz(j)], vv = f.v[zz(j) * zz(p) + zz(j)];
    const double se = vv >= 0 ? std::sqrt(vv) : NaN;
    e.push_back(b);
    s.push_back(na(se));
    lo.push_back(na(b - tq * se));
    hi.push_back(na(b + tq * se));
    st.push_back(na(b / se));
    pv.push_back(na(2 * clab_pt(-std::fabs(b / se), df, 1)));
    ++j;
  }
  Frame tab;
  tab.num("evtm", ev);
  tab.num("est", e);
  tab.num("se", s);
  tab.num("lo", lo);
  tab.num("hi", hi);
  tab.num("stat", st);
  tab.num("p", pv);
  std::vector<int> pre;
  for (int j = 0; j < q; ++j) {
    if (ks[zz(j)] < 0) pre.push_back(j);
  }
  const int nq = static_cast<int>(pre.size());
  double fst = NA, fp = NA;
  if (nq > 0) {
    std::vector<double> vp(zz(nq) * zz(nq)), vi(zz(nq) * zz(nq));
    for (int a = 0; a < nq; ++a) {
      for (int b = 0; b < nq; ++b) vp[zz(b) * zz(nq) + zz(a)] = f.v[zz(pre[zz(b)]) * zz(p) + zz(pre[zz(a)])];
    }
    if (pev::chol(vp.data(), nq)) {
      pev::cinv(vp.data(), nq, vi.data());
      double wd = 0.0;
      for (int a = 0; a < nq; ++a) {
        for (int b = 0; b < nq; ++b) wd += f.b[zz(pre[zz(a)])] * vi[zz(b) * zz(nq) + zz(a)] * f.b[zz(pre[zz(b)])];
      }
      fst = wd / nq;
      fp = clab_pf(fst, nq, df, 0);
    }
  }
  Frame pt;
  pt.num("f", &fst, 1);
  const double nqd = nq;
  pt.num("df1", &nqd, 1);
  pt.num("df2", &df, 1);
  pt.num("p", &fp, 1);
  std::vector<double> vc(zz(q) * zz(q));
  for (int a = 0; a < q; ++a) {
    for (int b = 0; b < q; ++b) vc[zz(b) * zz(q) + zz(a)] = f.v[zz(b) * zz(p) + zz(a)];
  }
  std::vector<double> cv;
  for (double v : gv) {
    if (std::isfinite(v) && std::find(cv.begin(), cv.end(), v) == cv.end()) cv.push_back(v);
  }
  std::sort(cv.begin(), cv.end());
  std::vector<double> rc, re, rv;
  for (double cg : cv) {
    std::vector<double> sm(zz(pn.nt), 0.0), cnt(zz(pn.nt), 0.0);
    for (std::size_t i = 0; i < zz(n); ++i) {
      if (gv[zz(pn.u[i])] == cg) {
        sm[zz(pn.t[i])] += pn.y[i];
        cnt[zz(pn.t[i])] += 1.0;
      }
    }
    std::vector<double> tt, mv;
    for (int a = 0; a < pn.nt; ++a) {
      if (cnt[zz(a)] > 0) {
        tt.push_back(pn.tv[zz(a)]);
        mv.push_back(sm[zz(a)] / cnt[zz(a)]);
      }
    }
    const int nc = static_cast<int>(tt.size());
    std::vector<double> av(zz(nc), NA);
    for (int a = mw - 1; a < nc; ++a) {
      double s2 = 0.0;
      for (int b = a - mw + 1; b <= a; ++b) s2 += mv[zz(b)];
      av[zz(a)] = s2 / mw;
    }
    double base = NA;
    for (int a = 0; a < nc; ++a) {
      if (tt[zz(a)] - cg == -1) base = av[zz(a)];
    }
    for (int a = 0; a < nc; ++a) {
      const double et = tt[zz(a)] - cg;
      if (et < -ld || et > lg) continue;
      rc.push_back(cg);
      re.push_back(et);
      double v = av[zz(a)];
      if (nm) v = std::isnan(base) || base == 0 || std::isnan(v) ? NA : v / base;
      rv.push_back(v);
    }
  }
  Frame ch;
  ch.num("cohort", rc);
  ch.num("evtm", re);
  ch.num("value", rv);
  Val sens, brk;
  if (c.has("mbar")) {
    Dv mb;
    CK(num(c, c["mbar"], "mbar", mb, -1, 0.0));
    const int nb = static_cast<int>(mb.n());
    int nd = 0;
    CK(wholed(c, "draws", 5000, 100, nd));
    if (nq == 0) return c.fail("`mbar` needs one lead at least besides the reference period.");
    std::vector<int> tg;
    if (!c.has("target")) {
      for (int j = 0; j < q; ++j) {
        if (ks[zz(j)] >= 0) tg.push_back(j);
      }
    } else {
      const double tv = val1(c["target"]);
      for (int j = 0; j < q; ++j) {
        if (ks[zz(j)] >= 0 && ks[zz(j)] == tv) tg.push_back(j);
      }
    }
    if (tg.empty()) return c.fail("`target` must be NULL or an event time from 0 estimated.");
    double kk = 0.0;
    for (int j : tg) kk += (ks[zz(j)] + 1.0) / static_cast<double>(tg.size());
    auto parts = [&](const double *b, double *cen, double *dm) {
      double cc = 0.0;
      for (int j : tg) cc += b[j] / static_cast<double>(tg.size());
      std::vector<std::pair<int, double>> pc;
      for (int j : pre) pc.emplace_back(ks[zz(j)], b[j]);
      if (r0 < 0) pc.emplace_back(r0, 0.0);
      std::sort(pc.begin(), pc.end());
      double dd = 0.0;
      for (std::size_t a = 1; a < pc.size(); ++a) dd = std::max(dd, std::fabs(pc[a].second - pc[a - 1].second));
      *cen = cc;
      *dm = dd;
    };
    double c0 = 0.0, d0 = 0.0;
    parts(f.b.data(), &c0, &d0);
    std::vector<double> lv2(zz(q) * zz(q)), cs(zz(nd)), ds(zz(nd)), bs(zz(q)), eps(zz(q));
    for (int a = 0; a < q; ++a) {
      for (int b = 0; b < q; ++b) lv2[zz(b) * zz(q) + zz(a)] = f.v[zz(b) * zz(p) + zz(a)] + (a == b ? 1e-12 : 0.0);
    }
    pev::chol(lv2.data(), q);
    unsigned int sd = 0;
    CK(seed(c, sd));
    kutil::rng rg(sd);
    for (int dr = 0; dr < nd; ++dr) {
      for (int a = 0; a < q; a += 2) {
        const double u1 = 1.0 - rg.unif(), u2 = rg.unif();
        const double rr = std::sqrt(-2.0 * std::log(u1));
        eps[zz(a)] = rr * std::cos(6.283185307179586 * u2);
        if (a + 1 < q) eps[zz(a + 1)] = rr * std::sin(6.283185307179586 * u2);
      }
      for (int a = 0; a < q; ++a) {
        double v = f.b[zz(a)];
        for (int b = 0; b <= a; ++b) v += lv2[zz(b) * zz(q) + zz(a)] * eps[zz(b)];
        bs[zz(a)] = v;
      }
      parts(bs.data(), &cs[zz(dr)], &ds[zz(dr)]);
    }
    const double al = 1.0 - lv;
    auto q7 = [](std::vector<double> xs, double pr) {
      std::sort(xs.begin(), xs.end());
      const double h = (static_cast<double>(xs.size()) - 1) * pr;
      const std::size_t l = static_cast<std::size_t>(std::floor(h));
      const std::size_t u = std::min(l + 1, xs.size() - 1);
      return xs[l] + (h - static_cast<double>(l)) * (xs[u] - xs[l]);
    };
    auto ends = [&](double mm, double *l, double *h) {
      std::vector<double> a(zz(nd)), b(zz(nd));
      for (int dr = 0; dr < nd; ++dr) {
        a[zz(dr)] = cs[zz(dr)] - mm * ds[zz(dr)] * kk;
        b[zz(dr)] = cs[zz(dr)] + mm * ds[zz(dr)] * kk;
      }
      *l = q7(a, al / 2);
      *h = q7(b, 1 - al / 2);
    };
    std::vector<double> sl(zz(nb)), su(zz(nb)), slo(zz(nb)), shi(zz(nb));
    for (int a = 0; a < nb; ++a) {
      const double mm = mb[zz(a)];
      sl[zz(a)] = c0 - mm * d0 * kk;
      su[zz(a)] = c0 + mm * d0 * kk;
      ends(mm, &slo[zz(a)], &shi[zz(a)]);
    }
    Frame sc;
    sc.num("mbar", mb.v);
    sc.num("lob", sl);
    sc.num("upb", su);
    sc.num("lo", slo);
    sc.num("hi", shi);
    sens = sc.o;
    auto covers = [&](double mm) {
      double l = 0.0, h = 0.0;
      ends(mm, &l, &h);
      return l <= 0 && h >= 0;
    };
    double bd = covers(0) ? 0.0 : Inf;
    if (bd > 0 && covers(10)) {
      double l = 0, h = 10;
      for (int it = 0; it < 40; ++it) {
        const double mid = 0.5 * (l + h);
        (covers(mid) ? h : l) = mid;
      }
      bd = h;
    }
    brk = vnum(bd);
  }
  out = vobj();
  put(out, "est", tab.o);
  put(out, "pre", pt.o);
  put(out, "vcov", vmat(vc, q, q));
  put(out, "coh", ch.o);
  put(out, "sens", sens);
  put(out, "brk", brk);
  put(out, "n", vnum(f.n));
  put(out, "ncl", vnum(f.g));
  put(out, "level", vnum(lv));
  return true;
}

// One estimate of a one-row data frame.
inline void erow(Frame &f, double e, double se, double zq) {
  col1(f, "est", e);
  col1(f, "se", se);
  col1(f, "lo", e - zq * se);
  col1(f, "hi", e + zq * se);
}

bool stgdid(Ctx &c, Val &out) {
  static const char *const mo[] = {"csa", "wool", "dcdh"};
  static const char *const co[] = {"notyet", "never"};
  static const char *const eo[] = {"dr", "ipw", "reg"};
  pev::stg s{};
  CK(pickd(c, "meth", mo, 3, 0, s.meth));
  CK(pickd(c, "ctrl", co, 2, 0, s.ctrl));
  CK(pickd(c, "est", eo, 3, 0, s.est));
  int nr = 0;
  CK(wholed(c, "reps", 499, 0, nr));
  double lv = 0.0;
  CK(oned(c, "level", 0.95, 0.0, 1.0, lv));
  Pan pn;
  CK(panel(c, -HUGE_VAL, pn));
  const int n = pn.n, nu = pn.nu, nt = pn.nt;
  if (n != nu * nt) return c.fail("the panel must be balanced: every unit in every period.");
  if (c.has("first") == c.has("trt")) return c.fail("give one of `first` and `trt`.");
  std::vector<int> g(zz(nu), nt), dw(zz(nu) * zz(nt), 0);
  if (c.has("trt")) {
    std::vector<int> dv;
    CK(bin01(c, c["trt"], "trt", zz(n), dv));
    for (std::size_t i = 0; i < zz(n); ++i) dw[zz(pn.t[i]) * zz(nu) + zz(pn.u[i])] = dv[i];
    if (s.meth != 2) CK(firsts(c, pn, dv.data(), g));
  } else {
    std::vector<double> fv;
    CK(firstcol(c, pn, fv));
    for (std::size_t a = 0; a < zz(nu); ++a) {
      const double v = fv[a];
      if (v > pn.tv.back()) continue;
      if (v <= pn.tv.front()) {
        g[a] = 0;
      } else {
        const auto it = std::lower_bound(pn.tv.begin(), pn.tv.end(), v);
        if (*it != v) return c.fail("`first` must be one of the periods of `time` (or NA, 0, Inf).");
        g[a] = static_cast<int>(it - pn.tv.begin());
      }
      for (int b = g[a]; b < nt; ++b) dw[zz(b) * zz(nu) + a] = 1;
    }
  }
  std::vector<Term> x;
  CK(adjcols(c, c["adj"], zz(n), x, "row of the panel"));
  s.k = static_cast<int>(x.size());
  if (s.k > 0 && s.meth != 0) return c.fail("`adj` is taken by `meth = \"csa\"` only.");
  std::vector<double> xw(zz(s.k) * zz(nt) * zz(nu));
  for (int j = 0; j < s.k; ++j) {
    for (std::size_t i = 0; i < zz(n); ++i) xw[(zz(j) * zz(nt) + zz(pn.t[i])) * zz(nu) + zz(pn.u[i])] = x[zz(j)].x[i];
  }
  std::vector<double> yw(zz(n));
  for (std::size_t i = 0; i < zz(n); ++i) yw[zz(pn.t[i]) * zz(nu) + zz(pn.u[i])] = pn.y[i];
  std::vector<int> cl(zz(nu));
  int ncl = nu;
  if (!c.has("clus")) {
    for (int a = 0; a < nu; ++a) cl[zz(a)] = a;
  } else {
    Lab l;
    CK(codes(c, c["clus"], "clus", zz(n), l));
    std::vector<int> seen(zz(nu), -1);
    for (std::size_t i = 0; i < zz(n); ++i) {
      const std::size_t a = zz(pn.u[i]);
      const int e = l.code[i];
      if (seen[a] >= 0 && seen[a] != e) return c.fail("`clus` must be constant within each unit.");
      seen[a] = e;
      cl[a] = e;
    }
    ncl = l.g;
    if (ncl < 2) return c.fail("`clus` has one cluster: the bootstrap needs two or more.");
  }
  unsigned int sd0 = 0;
  if (nr > 0) CK(seed(c, sd0));
  const std::uint64_t sd = sd0;
  s.yw = yw.data();
  s.dw = dw.data();
  s.xw = xw.data();
  s.g = g.data();
  s.nu = nu;
  s.nt = nt;
  std::vector<int> idx;
  for (int a = 0; a < nu; ++a) {
    if (s.meth == 2 || g[zz(a)] > 0) idx.push_back(a);
  }
  const double zq = clab_qnorm(0.5 + 0.5 * lv);
  Val att, dyn, ov;
  if (s.meth == 2) {
    double v[4] = {0}, se[4] = {NA, NA, NA, NA};
    std::vector<int> ev;
    pev::stat(s, idx, ev, v);
    if (std::isnan(v[0])) return c.fail("no unit switches treatment against a unit that did not.");
    if (nr > 0) pev::boot(s, cl.data(), ncl, nr, sd, ev, 4, se);
    Frame f;
    f.str("term", {"effect", "placebo"});
    std::vector<double> e(2), sv(2), lo(2), hi(2), ns(2);
    for (int q = 0; q < 2; ++q) {
      e[zz(q)] = na(v[q]);
      sv[zz(q)] = na(se[q]);
      lo[zz(q)] = na(v[q] - zq * se[q]);
      hi[zz(q)] = na(v[q] + zq * se[q]);
      ns[zz(q)] = v[q + 2];
    }
    f.num("est", e);
    f.num("se", sv);
    f.num("lo", lo);
    f.num("hi", hi);
    f.num("n", ns);
    ov = f.o;
  } else {
    pev::cells(g.data(), nu, nt, s.cg, s.ct, s.meth == 1 ? 1 : 0);
    if (s.cg.empty()) return c.fail("no unit is first treated after the first period.");
    std::vector<int> ev;
    for (std::size_t q = 0; q < s.cg.size(); ++q) ev.push_back(s.ct[q] - s.cg[q]);
    std::sort(ev.begin(), ev.end());
    ev.erase(std::unique(ev.begin(), ev.end()), ev.end());
    const int nc = static_cast<int>(s.cg.size()), ne = static_cast<int>(ev.size());
    const int ns = nc + 1 + ne;
    std::vector<double> v(zz(ns)), se(zz(ns), NA);
    pev::stat(s, idx, ev, v.data());
    if (nr > 0) pev::boot(s, cl.data(), ncl, nr, sd, ev, ns, se.data());
    std::vector<double> ng(zz(nt), 0.0);
    for (int a : idx) {
      if (g[zz(a)] < nt) ng[zz(g[zz(a)])] += 1.0;
    }
    std::vector<double> cc, cp, ce, e, sv, lo, hi, nn;
    for (std::size_t q = 0; q < zz(nc); ++q) {
      cc.push_back(pn.tv[zz(s.cg[q])]);
      cp.push_back(pn.tv[zz(s.ct[q])]);
      ce.push_back(s.ct[q] - s.cg[q]);
      e.push_back(na(v[q]));
      sv.push_back(na(se[q]));
      lo.push_back(na(v[q] - zq * se[q]));
      hi.push_back(na(v[q] + zq * se[q]));
      nn.push_back(ng[zz(s.cg[q])]);
    }
    Frame f;
    f.num("cohort", cc);
    f.num("period", cp);
    f.num("evtm", ce);
    f.num("est", e);
    f.num("se", sv);
    f.num("lo", lo);
    f.num("hi", hi);
    f.num("ncoh", nn);
    att = f.o;
    Frame o;
    o.str("term", {"att"});
    erow(o, v[zz(nc)], se[zz(nc)], zq);
    ov = o.o;
    std::vector<double> de, dv, ds, dl, dh;
    for (int q = 0; q < ne; ++q) {
      const std::size_t a = zz(nc + 1 + q);
      de.push_back(ev[zz(q)]);
      dv.push_back(na(v[a]));
      ds.push_back(na(se[a]));
      dl.push_back(na(v[a] - zq * se[a]));
      dh.push_back(na(v[a] + zq * se[a]));
    }
    Frame d;
    d.num("evtm", de);
    d.num("est", dv);
    d.num("se", ds);
    d.num("lo", dl);
    d.num("hi", dh);
    dyn = d.o;
  }
  out = vobj();
  put(out, "att", att);
  put(out, "overall", ov);
  put(out, "dynamic", dyn);
  put(out, "meth", vstr(mo[s.meth]));
  put(out, "ctrl", vstr(co[s.ctrl]));
  put(out, "reps", vnum(nr));
  put(out, "level", vnum(lv));
  return true;
}

// widen(): a balanced panel wide, its treated units (one adoption) and
// its comparison units.
struct Wide {
  std::vector<double> yw;
  std::vector<int> trd, ctl;
  int t0 = 0;
};
bool widen(Ctx &c, const Pan &pn, Wide &w) {
  if (pn.n != pn.nu * pn.nt) return c.fail("the panel must be balanced: every unit in every period.");
  w.yw.resize(zz(pn.n));
  for (std::size_t i = 0; i < zz(pn.n); ++i) w.yw[zz(pn.t[i]) * zz(pn.nu) + zz(pn.u[i])] = pn.y[i];
  std::vector<double> fv;
  CK(firstcol(c, pn, fv));
  double g = Inf;
  for (int a = 0; a < pn.nu; ++a) {
    if (!std::isfinite(fv[zz(a)])) {
      w.ctl.push_back(a);
      continue;
    }
    if (std::isfinite(g) && fv[zz(a)] != g) return c.fail("the treated units must share one adoption period.");
    g = fv[zz(a)];
    w.trd.push_back(a);
  }
  if (w.trd.empty()) return c.fail("no unit is treated: `first` is NA, 0 or Inf throughout.");
  if (w.ctl.size() < 2) return c.fail("the comparison needs two units never treated at least.");
  const auto it = std::find(pn.tv.begin(), pn.tv.end(), g);
  if (it == pn.tv.begin() || it == pn.tv.end()) return c.fail("the adoption period must be a period of `time` after the first.");
  w.t0 = static_cast<int>(it - pn.tv.begin());
  return true;
}

bool synth(const Wide &w, int nu, int nt, int i, const std::vector<int> &d, double *om, double *gap) {
  const int n0 = static_cast<int>(d.size());
  std::vector<double> a(zz(w.t0) * zz(n0)), b(zz(w.t0)), g(zz(n0) * zz(n0));
  for (int s = 0; s < w.t0; ++s) {
    b[zz(s)] = w.yw[zz(s) * zz(nu) + zz(i)];
    for (int q = 0; q < n0; ++q) a[zz(q) * zz(w.t0) + zz(s)] = w.yw[zz(s) * zz(nu) + zz(d[zz(q)])];
  }
  pev::gram(a.data(), b.data(), w.t0, n0, 0, 0.0, g.data());
  if (pev::mnp(g.data(), n0, om, 10000) < 0) return false;
  for (int s = 0; s < nt; ++s) {
    double v = w.yw[zz(s) * zz(nu) + zz(i)];
    for (int q = 0; q < n0; ++q) v -= om[q] * w.yw[zz(s) * zz(nu) + zz(d[zz(q)])];
    gap[s] = v;
  }
  return true;
}
double rmspe(const double *gap, int t0, int nt, bool post) {
  double s = 0.0;
  const int lo = post ? t0 : 0, hi = post ? nt : t0;
  for (int q = lo; q < hi; ++q) s += gap[q] * gap[q];
  return std::sqrt(s / (hi - lo));
}
const char *const wfail = "the weights could not be solved (the units' paths are degenerate).";

bool synthc(Ctx &c, Val &out) {
  bool dp = true;
  CK(flagd(c, "plac", true, dp));
  Pan pn;
  CK(panel(c, -HUGE_VAL, pn));
  Wide w;
  CK(widen(c, pn, w));
  if (w.trd.size() != 1) return c.fail("one unit must be treated: `first` is finite for %d.", static_cast<int>(w.trd.size()));
  const int nu = pn.nu, nt = pn.nt, n0 = static_cast<int>(w.ctl.size()), i1 = w.trd[0];
  std::vector<double> om(zz(n0)), gap(zz(nt));
  if (!synth(w, nu, nt, i1, w.ctl, om.data(), gap.data())) return c.fail("%s", wfail);
  const double pre = rmspe(gap.data(), w.t0, nt, false), post = rmspe(gap.data(), w.t0, nt, true);
  const double ratio = post / pre;
  Val pl;
  double pv = NA;
  if (dp) {
    std::vector<std::string> cu;
    std::vector<double> cr;
    std::vector<double> o2(zz(n0 - 1)), g2(zz(nt));
    int ge = 0;
    for (int q = 0; q < n0; ++q) {
      std::vector<int> d;
      for (int r = 0; r < n0; ++r) {
        if (r != q) d.push_back(w.ctl[zz(r)]);
      }
      const bool ok = synth(w, nu, nt, w.ctl[zz(q)], d, o2.data(), g2.data());
      const double rq = ok ? rmspe(g2.data(), w.t0, nt, true) / rmspe(g2.data(), w.t0, nt, false) : NA;
      cu.push_back(pn.ulev[zz(w.ctl[zz(q)])]);
      cr.push_back(rq);
      ge += !std::isnan(rq) && rq >= ratio ? 1 : 0;
    }
    Frame f;
    f.str("unit", cu);
    f.num("ratio", cr);
    pl = f.o;
    pv = (1.0 + ge) / (n0 + 1.0);
  }
  std::vector<std::string> wu;
  for (int q = 0; q < n0; ++q) wu.push_back(pn.ulev[zz(w.ctl[zz(q)])]);
  Frame wt;
  wt.str("unit", wu);
  wt.num("w", om);
  std::vector<double> pt(zz(nt)), py(zz(nt)), ps(zz(nt));
  for (std::size_t s = 0; s < zz(nt); ++s) {
    pt[s] = pn.tv[s];
    py[s] = w.yw[s * zz(nu) + zz(i1)];
    ps[s] = py[s] - gap[s];
  }
  Frame path;
  path.num("time", pt);
  path.num("treated", py);
  path.num("synth", ps);
  path.num("gap", gap);
  Frame fit;
  fit.num("rmspepre", &pre, 1);
  fit.num("rmspepost", &post, 1);
  col1(fit, "ratio", ratio);
  fit.num("p", &pv, 1);
  out = vobj();
  put(out, "weights", wt.o);
  put(out, "path", path.o);
  put(out, "fit", fit.o);
  put(out, "plac", pl);
  put(out, "unit", vstr(pn.ulev[zz(i1)]));
  return true;
}

bool sdidx(Ctx &c, Val &out) {
  int nr = 0;
  CK(wholed(c, "reps", 200, 2, nr));
  double lv = 0.0;
  CK(oned(c, "level", 0.95, 0.0, 1.0, lv));
  Pan pn;
  CK(panel(c, -HUGE_VAL, pn));
  Wide w;
  CK(widen(c, pn, w));
  const int nu = pn.nu, nt = pn.nt, n0 = static_cast<int>(w.ctl.size());
  const int n1 = static_cast<int>(w.trd.size());
  std::vector<double> om(zz(n0)), la(zz(w.t0));
  const double tau = pev::sdtau(w.yw.data(), nu, nt, w.ctl, w.trd, w.t0, true, om.data(), la.data());
  if (std::isnan(tau)) return c.fail("%s", wfail);
  double se = NA;
  const char *meth = "placebo";
  if (n1 == 1) {
    unsigned int sd = 0;
    CK(seed(c, sd));
    kutil::rng g(sd);
    std::vector<double> o2(zz(n0 - 1)), l2(zz(w.t0));
    kutil::ksum s1, s2;
    double k = 0;
    for (int r = 0; r < nr; ++r) {
      std::vector<int> cc = w.ctl;
      const int pk = g.below(n0);
      const std::vector<int> t1 = {cc[zz(pk)]};
      cc.erase(cc.begin() + pk);
      const double v = pev::sdtau(w.yw.data(), nu, nt, cc, t1, w.t0, true, o2.data(), l2.data());
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
    meth = "jackknife";
    std::vector<double> u;
    for (int a = 0; a < n0 + n1; ++a) {
      std::vector<int> cc, t;
      std::vector<double> o2;
      double sw = 0.0;
      for (int q = 0; q < n0; ++q) {
        if (q == a) continue;
        cc.push_back(w.ctl[zz(q)]);
        o2.push_back(om[zz(q)]);
        sw += om[zz(q)];
      }
      for (int q = 0; q < n1; ++q) {
        if (q + n0 != a) t.push_back(w.trd[zz(q)]);
      }
      if (!(sw > 0)) {
        u.clear();
        break;
      }
      for (double &v : o2) v /= sw;
      u.push_back(pev::sdtau(w.yw.data(), nu, nt, cc, t, w.t0, false, o2.data(), la.data()));
    }
    if (!u.empty()) {
      const double nn = static_cast<double>(u.size());
      double m = 0.0, ss = 0.0;
      for (double v : u) m += v / nn;
      for (double v : u) ss += (v - m) * (v - m);
      se = std::sqrt((nn - 1.0) / nn * ss);
    }
  }
  const double zq = clab_qnorm(0.5 + 0.5 * lv);
  Frame est;
  const double lo = std::isnan(se) ? NA : tau - zq * se, hi = std::isnan(se) ? NA : tau + zq * se;
  est.num("est", &tau, 1);
  est.num("se", &se, 1);
  est.num("lo", &lo, 1);
  est.num("hi", &hi, 1);
  std::vector<std::string> uu;
  for (int q = 0; q < n0; ++q) uu.push_back(pn.ulev[zz(w.ctl[zz(q)])]);
  Frame unitw;
  unitw.str("unit", uu);
  unitw.num("w", om);
  Frame timew;
  timew.num("time", pn.tv.data(), zz(w.t0));
  timew.num("w", la);
  std::vector<double> pt(zz(nt)), py(zz(nt), 0.0), ps(zz(nt), 0.0);
  for (std::size_t s = 0; s < zz(nt); ++s) {
    pt[s] = pn.tv[s];
    for (int i : w.trd) py[s] += w.yw[s * zz(nu) + zz(i)] / n1;
    for (int q = 0; q < n0; ++q) ps[s] += om[zz(q)] * w.yw[s * zz(nu) + zz(w.ctl[zz(q)])];
  }
  Frame path;
  path.num("time", pt);
  path.num("treated", py);
  path.num("synth", ps);
  out = vobj();
  put(out, "est", est.o);
  put(out, "meth", vstr(meth));
  put(out, "unitw", unitw.o);
  put(out, "timew", timew.o);
  put(out, "path", path.o);
  put(out, "level", vnum(lv));
  return true;
}
