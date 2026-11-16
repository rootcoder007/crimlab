// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// idwint(), gtwreg() and asrule() (src/idw.cpp, src/gtw.cpp, src/arl.cpp)
// for clab_jrun(): the R wrappers of R/spi.R in C++, with their checks in
// their order and words. Labels are numbers when every label is a number,
// else text, sorted as numbers or byte by byte (R sorts text by locale).

namespace spi {

inline bool numvec(Ctx &c, const char *nm, Nv &d) { return cknum(c, c[nm], nm, d, false); }

inline bool radius(Ctx &c, double &rad) {
  rad = HUGE_VAL;
  if (!c.has("rad")) return true;
  const std::vector<double> r = ax::nums(c["rad"]);
  if (!ax::numv(c["rad"]) || r.size() != 1 || std::isnan(r[0]) || !(r[0] > 0)) {
    return ckstop(c, "`rad` must be a positive number or Inf");
  }
  rad = r[0];
  return true;
}

inline Val nav(const std::vector<double> &v) {
  Val a = varr();
  for (double x : v) a.a.push_back(std::isnan(x) ? Val() : vnum(x));
  return a;
}

}  // namespace spi

bool idwint(Ctx &c, Val &out) {
  Nv mx, my, tx, ty, p, k;
  CK(spi::numvec(c, "mx", mx));
  CK(spi::numvec(c, "my", my));
  const std::vector<double> mv = ax::nums(c["mv"]);
  if (!ax::numv(c["mv"]) || my.n() != mx.n() || mv.size() != mx.n()) {
    return ckstop(c, "`mx`, `my` and `mv` must be numeric vectors of one length");
  }
  CK(spi::numvec(c, "tx", tx));
  CK(spi::numvec(c, "ty", ty));
  if (ty.n() != tx.n()) return ckstop(c, "`tx` and `ty` must have one length");
  if (c.has("p")) {
    CK(cknum(c, c["p"], "p", p, true, 0));
  } else {
    p.v.push_back(2);
  }
  if (c.has("k")) {
    CK(cknum(c, c["k"], "k", k, true, 0, 2147483647, false, true));
  } else {
    k.v.push_back(0);
  }
  double rad = 0;
  CK(spi::radius(c, rad));
  bool ll = false;
  CK(dw::flg(c, "lonlat", false, ll));
  const std::size_t nm = mx.n(), nt = tx.n();
  if (nm > 2147483647u || nt > 2147483647u) return ckstop(c, "too many monitors or targets");
  dw::Lab z;
  const bool hz = c.has("zone");
  if (hz && !dw::labs(c["zone"], nt, true, z)) {
    return ckstop(c, "`zone` must give the zone of each target, with no NA");
  }
  std::vector<double> o(nt), u(nt), zo(hz ? z.nlev : 0), zn(hz ? z.nlev : 0);
  const int rc = clab_idwint(static_cast<int>(nm), mx.p(), my.p(), mv.data(), static_cast<int>(nt), tx.p(), ty.p(),
                             p.v[0], static_cast<int>(k.v[0]), rad, ll ? 1 : 0, hz ? z.code.data() : nullptr,
                             hz ? static_cast<int>(z.nlev) : 0, o.data(), u.data(), hz ? zo.data() : nullptr,
                             hz ? zn.data() : nullptr);
  if (rc != CLAB_OK) return c.fail("%s", clab_strerror(rc));
  Frame f;
  if (!hz) {
    f.add("val", spi::nav(o));
    f.num("used", u);
  } else {
    f.add("zone", dw::labv(z, 1, z.nlev, false));
    f.add("val", spi::nav(zo));
    f.num("n", zn);
  }
  out = f.o;
  return true;
}

namespace spi {

// One fit of gtwreg() at (h, tau).
struct Gfit {
  std::vector<double> b, se, fit, hat, loo, st;
};

inline int gfit(const std::vector<double> &x, int k, const Nv &y, const Nv &u, const Nv &v, const Nv *t, double h,
                double tau, bool ll, Gfit &g) {
  const std::size_t n = y.n();
  g.b.assign(n * zz(k), 0.0);
  g.se.assign(n * zz(k), 0.0);
  g.fit.assign(n, 0.0);
  g.hat.assign(n, 0.0);
  g.loo.assign(n, 0.0);
  g.st.assign(7, 0.0);
  return clab_gtwreg(static_cast<int>(n), k, x.data(), y.p(), u.p(), v.p(), t ? t->p() : nullptr, h, tau, ll ? 1 : 0,
                     g.b.data(), g.se.data(), g.fit.data(), g.hat.data(), g.loo.data(), g.st.data());
}

}  // namespace spi

bool gtwreg(Ctx &c, Val &out) {
  Nv y, u, v, t, h, tau;
  CK(cknum(c, c["y"], "y", y, false));
  const std::size_t n = y.n();
  if (n > 2147483647u) return ckstop(c, "`y` has too many values");
  bool cons = true, ll = false;
  CK(dw::flg(c, "cons", true, cons));
  dpnd::xm xm;
  if (c.has("x")) CK(env::rdm(c, c["x"], "x", "x", static_cast<int>(n), -HUGE_VAL, xm));
  const int k = xm.p + (cons ? 1 : 0);
  if (k == 0) return ckstop(c, "the model needs `x` or `cons = TRUE`");
  std::vector<double> x;
  std::vector<std::string> cn;
  if (cons) {
    x.assign(n, 1.0);
    cn.push_back("(const)");
  }
  x.insert(x.end(), xm.v.begin(), xm.v.end());
  cn.insert(cn.end(), xm.nm.begin(), xm.nm.end());
  CK(cknum(c, c["u"], "u", u, false));
  CK(cknum(c, c["v"], "v", v, false));
  if (u.n() != n || v.n() != n) return ckstop(c, "`u` and `v` must have one value per value of `y`");
  const bool ht = c.has("t");
  if (ht) {
    CK(cknum(c, c["t"], "t", t, false));
    if (t.n() != n) return ckstop(c, "`t` must have one value per value of `y`");
  }
  if (n <= zz(k) + 1) return ckstop(c, "the model needs more observations than coefficients plus one");
  CK(cknum(c, c["h"], "h", h, false, -HUGE_VAL, HUGE_VAL, true));
  if (c.has("tau")) {
    CK(cknum(c, c["tau"], "tau", tau, false, 0));
  } else {
    tau.v.push_back(0);
  }
  int crit = 0;
  CK(dw::pick(c, "crit", {"cv", "aicc"}, crit));
  CK(dw::flg(c, "lonlat", false, ll));
  // the grid: h fastest, as expand.grid()
  std::vector<double> gh, gt, gc, ga;
  spi::Gfit g;
  for (double tt : tau.v) {
    for (double hh : h.v) {
      const int rc = spi::gfit(x, k, y, u, v, ht ? &t : nullptr, hh, tt, ll, g);
      if (rc == CLAB_EDOM) return c.fail("a local fit is singular: widen the bandwidth or drop a regressor.");
      if (rc != CLAB_OK) return c.fail("%s", clab_strerror(rc));
      gh.push_back(hh);
      gt.push_back(tt);
      gc.push_back(g.st[4]);
      ga.push_back(g.st[5]);
    }
  }
  const std::vector<double> &cr = crit == 0 ? gc : ga;
  std::size_t b = cr.size();
  for (std::size_t i = 0; i < cr.size(); ++i) {
    if (!std::isnan(cr[i]) && (b == cr.size() || cr[i] < cr[b])) b = i;
  }
  if (b == cr.size()) return ckstop(c, "no bandwidth gives a finite `%s`", crit == 0 ? "cv" : "aicc");
  const int rc = spi::gfit(x, k, y, u, v, ht ? &t : nullptr, gh[b], gt[b], ll, g);
  if (rc != CLAB_OK) return c.fail("%s", clab_strerror(rc));
  std::vector<double> res(n);
  long double my = 0;
  for (std::size_t i = 0; i < n; ++i) my += y.v[i];
  const double ym = static_cast<double>(my / static_cast<long double>(n));
  double tss = 0;
  for (std::size_t i = 0; i < n; ++i) {
    res[i] = y.v[i] - g.fit[i];
    tss += (y.v[i] - ym) * (y.v[i] - ym);
  }
  out = vobj();
  put(out, "coef", vmat(g.b, static_cast<int>(n), k));
  put(out, "se", vmat(g.se, static_cast<int>(n), k));
  put(out, "fit", vnums(g.fit));
  put(out, "resid", vnums(res));
  put(out, "hat", vnums(g.hat));
  put(out, "loo", spi::nav(g.loo));
  static const char *const sn[] = {"h", "tau", "rss", "trs", "trss", "sigma2", "cv", "aicc", "edf", "r2"};
  const double sv[] = {gh[b], gt[b], g.st[0], g.st[1], g.st[2], g.st[3], g.st[4], g.st[5], g.st[6], 1 - g.st[0] / tss};
  put(out, "stat", nvec(10, sn, sv));
  Frame gf;
  gf.num("h", gh);
  gf.num("tau", gt);
  gf.add("cv", spi::nav(gc));
  gf.add("aicc", spi::nav(ga));
  put(out, "grid", gf.o);
  return true;
}

bool asrule(Ctx &c, Val &out) {
  const Val x = ax::asarr(c["x"]), y = ax::asarr(c["y"]);
  bool na = x.a.size() != y.a.size() || x.a.empty();
  for (const Val &e : x.a) na = na || ax::vna(e);
  for (const Val &e : y.a) na = na || ax::vna(e);
  if (na) return ckstop(c, "`x` and `y` must have one length and no NA");
  Nv ms, mc;
  if (c.has("minsup")) {
    CK(cknum(c, c["minsup"], "minsup", ms, true, 0, 1));
  } else {
    ms.v.push_back(0.1);
  }
  if (c.has("minconf")) {
    CK(cknum(c, c["minconf"], "minconf", mc, true, 0, 1));
  } else {
    mc.v.push_back(0.2);
  }
  bool empty = true;
  CK(dw::flg(c, "empty", true, empty));
  std::vector<std::string> nm = {"x", "y"};
  if (c.has("nm")) {
    const Val a = ax::asarr(c["nm"]);
    bool ok = a.a.size() == 2;
    for (const Val &e : a.a) ok = ok && e.t == Val::STR;
    if (!ok) return ckstop(c, "`nm` must be two names");
    nm = {a.a[0].s, a.a[1].s};
  }
  const std::size_t n = x.a.size();
  if (n > 2147483647u) return ckstop(c, "`x` has too many values");
  dw::Lab lx, ly;
  dw::labs(c["x"], n, true, lx);
  dw::labs(c["y"], n, true, ly);
  const double mxr = static_cast<double>(lx.nlev + ly.nlev) + 2.0 * static_cast<double>(lx.nlev) * static_cast<double>(ly.nlev);
  if (mxr > kutil::maxcell) return ckstop(c, "the levels of `x` and `y` are too many");
  const std::size_t m = static_cast<std::size_t>(mxr);
  std::vector<int> lh(m), rh(m);
  std::vector<double> o(m * 5);
  int nr = 0;
  const int rc = clab_asrule(static_cast<int>(n), lx.code.data(), static_cast<int>(lx.nlev), ly.code.data(),
                             static_cast<int>(ly.nlev), ms.v[0], mc.v[0], empty ? 1 : 0, &nr, lh.data(), rh.data(),
                             o.data());
  if (rc != CLAB_OK) return c.fail("%s", clab_strerror(rc));
  auto lab = [&](const dw::Lab &l, std::size_t i) {
    if (!l.num) return l.sv[i];
    return rfmt(l.nv[i]);
  };
  std::vector<std::string> it;
  for (std::size_t i = 0; i < lx.nlev; ++i) it.push_back(nm[0] + "=" + lab(lx, i));
  for (std::size_t i = 0; i < ly.nlev; ++i) it.push_back(nm[1] + "=" + lab(ly, i));
  std::vector<std::string> L, R;
  std::vector<double> col[5];
  for (int r = 0; r < nr; ++r) {
    L.push_back(lh[zz(r)] < 0 ? std::string() : it[zz(lh[zz(r)])]);
    R.push_back(it[zz(rh[zz(r)])]);
    for (int j = 0; j < 5; ++j) col[j].push_back(o[zz(r) * 5 + zz(j)]);
  }
  Frame f;
  f.str("lhs", L);
  f.str("rhs", R);
  static const char *const cnm[] = {"support", "confidence", "coverage", "lift", "count"};
  for (int j = 0; j < 5; ++j) f.num(cnm[j], col[j]);
  out = f.o;
  return true;
}

namespace ax {

// ifesc() (src/ifx.cpp): R/ifx.R in C++, its checks in their order and
// words. Units in the order they first appear, times sorted by value.
inline Val ifrow(const std::vector<double> &b, const std::vector<double> &se, double lv) {
  const double q = clab_qnorm(0.5 + lv / 2);
  std::vector<double> st, p, lo, hi;
  for (std::size_t j = 0; j < b.size(); ++j) {
    const bool na = std::isnan(se[j]);
    st.push_back(na ? NA : b[j] / se[j]);
    p.push_back(na ? NA : 2 * clab_pnorm(-std::fabs(b[j] / se[j])));
    lo.push_back(na ? NA : b[j] - q * se[j]);
    hi.push_back(na ? NA : b[j] + q * se[j]);
  }
  Frame f;
  f.num("est", b);
  f.num("se", se);
  f.num("stat", st);
  f.num("p", p);
  f.num("lo", lo);
  f.num("hi", hi);
  return f.o;
}

struct Ipan {
  std::vector<double> y, x, tt;
  std::vector<int> a;
  Cols units;
  int n = 0, nt = 0, p = 0;
};

inline bool ifpan(Ctx &c, const Val &d, Cols &x, Ipan &o) {
  Cols y, tr, un, tm;
  CK(acol(c, d, "y", {1}, false, y));
  CK(acol(c, d, "d", {1}, false, tr));
  CK(acol(c, d, "unit", {1}, false, un));
  CK(acol(c, d, "time", {1}, false, tm));
  CK(acol(c, d, "x", {}, true, x));
  std::vector<double> yv, xm, dv;
  CK(amat(c, d, y, yv));
  CK(amat(c, d, x, xm));
  CK(amat(c, d, tr, dv));
  for (double v : dv) {
    if (v != 0 && v != 1) return ckstop(c, "column `%s` must hold 0 or 1", tr[0].c_str());
  }
  const Val &tc = *fc(d, tm[0]);
  if (!isnumc(tc)) return ckstop(c, "column `%s` must hold numbers or dates", tm[0].c_str());
  const std::vector<double> tv = colv(tc);
  const std::size_t m = nr(d);
  std::map<std::string, int> um;
  std::vector<int> iu(m), it(m);
  for (std::size_t i = 0; i < m; ++i) {
    const std::string l = vlab((*fc(d, un[0])).a[i]);
    const auto e = um.emplace(l, static_cast<int>(o.units.size()));
    if (e.second) o.units.push_back(l);
    iu[i] = e.first->second;
  }
  o.tt = tv;
  std::sort(o.tt.begin(), o.tt.end());
  o.tt.erase(std::unique(o.tt.begin(), o.tt.end()), o.tt.end());
  o.n = static_cast<int>(o.units.size());
  o.nt = static_cast<int>(o.tt.size());
  o.p = static_cast<int>(x.size());
  const std::size_t sn = zz(o.n), st = zz(o.nt);
  std::vector<char> seen(sn * st, 0);
  bool bal = m == sn * st;
  for (std::size_t i = 0; i < m && bal; ++i) {
    it[i] = static_cast<int>(std::lower_bound(o.tt.begin(), o.tt.end(), tv[i]) - o.tt.begin());
    char &s = seen[zz(iu[i]) * st + zz(it[i])];
    bal = !s;
    s = 1;
  }
  if (!bal) return ckstop(c, "the panel must be balanced: each unit once at each time");
  o.y.assign(sn * st, 0.0);
  o.x.assign(sn * st * zz(o.p), 0.0);
  std::vector<double> dm(sn * st);
  for (std::size_t i = 0; i < m; ++i) {
    const std::size_t k = zz(iu[i]) * st + zz(it[i]);
    o.y[k] = yv[i];
    dm[k] = dv[i];
    for (std::size_t j = 0; j < zz(o.p); ++j) o.x[j * sn * st + k] = xm[j * m + i];
  }
  bool mono = true;
  o.a.assign(sn, o.nt);
  for (std::size_t u = 0; u < sn; ++u) {
    for (std::size_t t = 0; t < st; ++t) {
      if (dm[u * st + t] == 1 && o.a[u] == o.nt) o.a[u] = static_cast<int>(t);
      if (t > 0 && dm[u * st + t] < dm[u * st + t - 1]) mono = false;
    }
  }
  if (!mono) return ckstop(c, "the treatment must stay 1 once it turns 1");
  int nco = 0;
  for (int v : o.a) {
    if (v == 0) {
      return ckstop(c, "a treated unit is treated at the first time: it has no periods before adoption");
    }
    nco += v == o.nt ? 1 : 0;
  }
  if (nco == o.n) return ckstop(c, "no unit is treated");
  if (nco < 2) return ckstop(c, "the panel needs two or more never-treated units");
  return true;
}

bool ifesc(Ctx &c, Val &out) {
  static const std::vector<const char *> ses = {"boot", "none"};
  int se = 0;
  CK(marg(c, "se", ses, se));
  Val d;
  CK(dat(c, c["data"], "data", d));
  Cols x;
  Ipan pn;
  CK(ifpan(c, d, x, pn));
  int r = -1;
  if (c.has("r") && c["r"].t != Val::NUL) {
    Nv q;
    CK(cknum(c, c["r"], "r", q, true, 0, Inf, false, true));
    r = q[0] > INT_MAX ? INT_MAX : static_cast<int>(q[0]);
  }
  Nv rm, tl, mx;
  CK(cknumd(c, "rmax", 5, rm, 0, 50, false, true));
  double lv = 0;
  CK(alevel(c, lv));
  CK(cknumd(c, "tol", 1e-8, tl, -Inf, Inf, true));
  CK(cknumd(c, "maxit", 1000, mx, 1, Inf, false, true));
  const int rmax = static_cast<int>(rm[0]), maxit = mx[0] > INT_MAX ? INT_MAX : static_cast<int>(mx[0]);
  int nb = 0;
  std::vector<int> bm;
  if (se == 0) {
    bool hs = false;
    double sd = 0;
    CK(bargs(c, nb, hs, sd));
    std::uint32_t s0 = 0;
    CK(adseed(c, hs, sd, s0));
    Adraw ad(s0, pn.n);
    std::vector<int> b1;
    for (std::size_t b0 = 0; b0 < zz(nb); b0 += ad.blk) {
      const std::size_t k = std::min(ad.blk, zz(nb) - b0);
      ad.next(k, b1);
      bm.insert(bm.end(), b1.begin(), b1.end());
    }
  }
  const std::size_t sn = zz(pn.n), st = zz(pn.nt), mxe = 2 * st + 1;
  if (static_cast<double>(nb) * static_cast<double>(mxe) > kutil::maxcell) {
    return ckstop(c, "the panel or the bootstrap is too large");
  }
  std::vector<double> beta(zz(pn.p)), y0(sn * st), fac(st * zz(std::max(r, rmax))), cv(zz(rmax) + 1), ef(mxe),
      bef(zz(nb) * mxe), s8(8);
  const int rc = clab_ifesc(pd(pn.y), pd(pn.x), pn.p, pn.a.data(), pn.n, pn.nt, r, rmax, maxit, tl[0],
                            bm.empty() ? nullptr : bm.data(), nb, beta.data(), y0.data(), fac.data(), cv.data(),
                            ef.data(), bef.data(), s8.data());
  if (rc == CLAB_EDOM) {
    return ckstop(c, "a treated unit has too few periods before adoption for the factors, or a fit is singular");
  }
  if (rc != CLAB_OK) return kerr(c, "ifesc", rc);
  const int rr = static_cast<int>(s8[0]), ne = static_cast<int>(s8[7]), emin = static_cast<int>(s8[6]);
  const std::size_t k = zz(ne) + 1;
  std::vector<double> sev(k, NA);
  Val draws;
  const int nf = static_cast<int>(s8[5]);
  if (se == 0) {
    if (nf > 0) {
      wrn(c, std::to_string(nf) + " of " + std::to_string(nb) +
                 " bootstrap draws could not be fitted and are left out.");
    }
    if (nb - nf < 2) return ckstop(c, "fewer than two bootstrap draws could be fitted");
    sev = csd(bef, zz(nb), static_cast<int>(k));
    draws = vmat(bef.data(), nb, static_cast<int>(k));
  }
  // eff: event times with their treated units counted
  std::vector<double> e, cnt(zz(ne), 0.0), b2(ef.begin() + 1, ef.begin() + static_cast<std::ptrdiff_t>(k)),
      s2(sev.begin() + 1, sev.end());
  for (int i = 0; i < ne; ++i) e.push_back(emin + i);
  Cols cu;
  std::vector<double> ct, ce, cy, c0, cf;
  for (std::size_t u = 0; u < sn; ++u) {
    if (pn.a[u] == pn.nt) continue;
    for (std::size_t t = 0; t < st; ++t) {
      const double ev = static_cast<double>(t) - pn.a[u];
      cnt[zz(static_cast<int>(ev) - emin)] += 1;
      cu.push_back(pn.units[u]);
      ct.push_back(pn.tt[t]);
      ce.push_back(ev);
      cy.push_back(pn.y[u * st + t]);
      c0.push_back(y0[u * st + t]);
      cf.push_back(pn.y[u * st + t] - y0[u * st + t]);
    }
  }
  Val eff = vobj();
  put(eff, "e", vnums(e));
  put(eff, "n", vnums(cnt));
  const Val er = ifrow(b2, s2, lv);
  for (std::size_t j = 0; j < er.k.size(); ++j) put(eff, er.k[j].c_str(), er.a[j]);
  Frame fc1;
  fc1.str("unit", cu);
  fc1.num("time", ct);
  fc1.num("e", ce);
  fc1.num("y", cy);
  fc1.num("y0", c0);
  fc1.num("eff", cf);
  Val cvv;
  if (r < 0) {
    std::vector<double> cr, cm;
    for (int i = 0; i <= rmax; ++i) {
      if (std::isnan(cv[zz(i)])) continue;
      cr.push_back(i);
      cm.push_back(cv[zz(i)]);
    }
    Frame f;
    f.num("r", cr);
    f.num("mspe", cm);
    cvv = f.o;
  }
  Val cvg;
  cvg.t = Val::BOOL;
  cvg.b = s8[2] == 1;
  out = vobj();
  put(out, "att", ifrow({ef[0]}, {sev[0]}, lv));
  put(out, "eff", eff);
  put(out, "r", vnum(rr));
  put(out, "cv", cvv);
  put(out, "beta", nvec(x, beta.data()));
  put(out, "cf", fc1.o);
  put(out, "factors", vmat(fac.data(), pn.nt, rr));
  put(out, "draws", draws);
  put(out, "ntr", vnum(s8[3]));
  put(out, "nco", vnum(s8[4]));
  put(out, "nt", vnum(pn.nt));
  put(out, "iter", vnum(s8[1]));
  put(out, "conv", cvg);
  put(out, "nboot", vnum(nb));
  put(out, "nfail", vnum(se == 0 ? nf : 0));
  put(out, "level", vnum(lv));
  return true;
}

// lrcce() (R/lrc.R): the object a data frame of unit estimates, one
// column per term; x an array of terms or an object of named arrays.
bool lrcce(Ctx &c, Val &out) {
  const Val &o = c["object"];
  bool okf = o.t == Val::OBJ && !o.k.empty();
  std::size_t n = 0;
  for (std::size_t j = 0; okf && j < o.a.size(); ++j) {
    okf = o.a[j].t == Val::ARR && numv(o.a[j]) && (j == 0 || o.a[j].a.size() == n);
    if (okf) n = o.a[j].a.size();
  }
  if (!okf) return ckstop(c, "`object` must be a ccepnl fit or a matrix of unit estimates with named columns");
  if (n < 2) return ckstop(c, "`object` needs two or more units");
  auto colof = [&](const std::string &k) -> const Val * {
    for (std::size_t j = 0; j < o.k.size(); ++j) {
      if (o.k[j] == k) return &o.a[j];
    }
    return nullptr;
  };
  auto strs = [](const Val &v, Cols &s) {
    if (v.t == Val::STR) {
      s = {v.s};
      return true;
    }
    if (v.t != Val::ARR) return false;
    for (const Val &e : v.a) {
      if (e.t != Val::STR) return false;
      s.push_back(e.s);
    }
    return true;
  };
  std::vector<Cols> xs;
  Cols nm;
  const Val &xv = c["x"];
  bool ok = true;
  if (xv.t == Val::OBJ) {
    for (std::size_t j = 0; j < xv.k.size() && ok; ++j) {
      Cols s;
      ok = strs(xv.a[j], s);
      xs.push_back(s);
      nm.push_back(xv.k[j]);
    }
  } else {
    Cols s;
    ok = strs(xv, s) && !s.empty();
    xs.push_back(s);
    nm.push_back(ok ? s[0] : "");
  }
  Cols yl;
  if (c.has("ylag") && c["ylag"].t != Val::NUL) {
    ok = ok && strs(c["ylag"], yl);
  } else {
    for (const std::string &k : o.k) {
      std::size_t q = 1;
      const bool m = k.size() > 3 && k[0] == 'L' && k.compare(k.size() - 2, 2, ".y") == 0;
      while (m && q < k.size() - 2 && std::isdigit(static_cast<unsigned char>(k[q]))) ++q;
      if (m && q == k.size() - 2 && q > 1) yl.push_back(k);
    }
  }
  std::vector<Cols> all = xs;
  all.push_back(yl);
  for (const Cols &s : all) {
    for (const std::string &k : s) ok = ok && colof(k) != nullptr;
  }
  if (!ok) return ckstop(c, "`x` and `ylag` must name columns of the unit estimates");
  if (yl.empty()) return ckstop(c, "the fit has no lags of `y`");
  double lv = 0;
  CK(alevel(c, lv));
  auto rsum = [&](const Cols &s, std::size_t i) {
    double v = 0.0;
    for (const std::string &k : s) v += colof(k)->a[i].n;
    return v;
  };
  const std::size_t m = xs.size() + 1;
  std::vector<double> u(n * m);
  for (std::size_t i = 0; i < n; ++i) {
    const double phi = 1 - rsum(yl, i);
    for (std::size_t j = 0; j < xs.size(); ++j) u[j * n + i] = rsum(xs[j], i) / phi;
    u[xs.size() * n + i] = phi;
  }
  for (double v : u) {
    if (!std::isfinite(v)) return ckstop(c, "a unit's lags of `y` sum to 1");
  }
  nm.push_back("phi");
  std::vector<double> mu(m), se(m);
  for (std::size_t j = 0; j < m; ++j) {
    long double s = 0.0, q = 0.0;
    for (std::size_t i = 0; i < n; ++i) s += u[j * n + i];
    mu[j] = static_cast<double>(s / static_cast<long double>(n));
    for (std::size_t i = 0; i < n; ++i) q += (u[j * n + i] - mu[j]) * (u[j * n + i] - mu[j]);
    se[j] = std::sqrt(static_cast<double>(q) / (static_cast<double>(n) * (static_cast<double>(n) - 1)));
  }
  const Val tb = atab(nm, mu, se, lv, NaN, false);
  out = vobj();
  for (std::size_t j = 0; j < tb.k.size(); ++j) put(out, tb.k[j] == "stat" ? "z" : tb.k[j].c_str(), tb.a[j]);
  put(out, "units", vmat(u.data(), static_cast<int>(n), static_cast<int>(m)));
  return true;
}

// cemtch() (src/cem.cpp): R/cem.R in C++. Numeric covariates cut into
// equal-width bins as seq(min, max, length.out = k + 1) makes them, k from
// cut or the rule (0 Sturges, 1 Scott); other columns matched on their
// labels.
inline bool cembr(Ctx &c, const std::vector<double> &z, const Val *k, int rule, const std::string &v,
                  const char *what, std::vector<double> &br) {
  const double lo = *std::min_element(z.begin(), z.end()), hi = *std::max_element(z.begin(), z.end());
  const std::string nm = std::string(what) + "$" + v;
  double nb = 0;
  const bool one = !k || k->t == Val::NUM || (k->t == Val::ARR && k->a.size() == 1);
  if (!k) {
    const double n = static_cast<double>(z.size());
    if (rule == 0) {
      nb = std::ceil(std::log2(n) + 1);
    } else {
      long double m = 0, q = 0;
      for (double x : z) m += x;
      m /= n;
      for (double x : z) q += (x - m) * (x - m);
      const double h = 3.5 * std::sqrt(static_cast<double>(q / (n - 1))) * std::pow(n, -1.0 / 3);
      nb = h > 0 ? std::max(1.0, std::ceil((hi - lo) / h)) : 1.0;
    }
  } else if (one) {
    const Val &e = k->t == Val::ARR ? k->a[0] : *k;
    nb = e.t == Val::NUM ? e.n : NaN;
    if (!std::isfinite(nb) || nb < 1 || nb != std::round(nb)) {
      return ckstop(c, "`%s` must be a whole number of bins from 1", nm.c_str());
    }
  }
  br.clear();
  if (one) {
    if (lo == hi) {
      br = {lo - 0.5, lo + 0.5};
      return true;
    }
    const double n1 = nb, del = hi - lo;
    br.push_back(lo);
    for (int j = 1; j < static_cast<int>(nb); ++j) br.push_back(lo + j * (del / n1));
    br.push_back(hi);
    return true;
  }
  bool ok = k->t == Val::ARR;
  for (std::size_t j = 0; ok && j < k->a.size(); ++j) {
    ok = k->a[j].t == Val::NUM && !std::isnan(k->a[j].n) && (j == 0 || k->a[j].n > k->a[j - 1].n);
    if (ok) br.push_back(k->a[j].n);
  }
  if (!ok || br.front() > lo || br.back() < hi) {
    return ckstop(c, "`%s` must be increasing breaks covering the column", nm.c_str());
  }
  return true;
}

inline bool cemcode(Ctx &c, const Val &d, const Cols &x, const char *what, int rule, std::vector<int> &code,
                    Val &cuts) {
  const Val &cv = c.has(what) ? c[what] : nul;
  if (cv.t != Val::NUL) {
    bool ok = cv.t == Val::OBJ && !cv.k.empty();
    for (const std::string &k : cv.k) ok = ok && has(x, k);
    if (!ok) return ckstop(c, "`%s` must be a list named by columns of `x`", what);
  }
  const std::size_t n = nr(d);
  code.assign(n * x.size(), 0);
  // list(): an empty array until a numeric column adds its breaks
  cuts = varr();
  for (std::size_t j = 0; j < x.size(); ++j) {
    const Val &col = *fc(d, x[j]);
    if (!isnumc(col)) {
      std::map<std::string, int> h;
      for (std::size_t i = 0; i < n; ++i) {
        code[j * n + i] = h.emplace(vlab(col.a[i]), static_cast<int>(h.size()) + 1).first->second;
      }
      continue;
    }
    const std::vector<double> z = colv(col);
    for (double v : z) {
      if (!std::isfinite(v)) return ckstop(c, "column `%s` must be finite", x[j].c_str());
    }
    const Val *k = nullptr;
    for (std::size_t q = 0; q < cv.k.size(); ++q) {
      if (cv.k[q] == x[j]) k = &cv.a[q];
    }
    std::vector<double> br;
    CK(cembr(c, z, k, rule, x[j], what, br));
    put(cuts, x[j].c_str(), vnums(br));
    // findInterval(rightmost.closed, all.inside)
    const int nb = static_cast<int>(br.size());
    for (std::size_t i = 0; i < n; ++i) {
      int b = static_cast<int>(std::upper_bound(br.begin(), br.end(), z[i]) - br.begin());
      if (z[i] == br.back()) b = nb - 1;
      code[j * n + i] = std::min(std::max(b, 1), nb - 1);
    }
  }
  return true;
}

bool cemtch(Ctx &c, Val &out) {
  Val d;
  CK(dat(c, c["data"], "data", d));
  Cols dc, x;
  CK(acol(c, d, "d", {1}, false, dc));
  CK(acol(c, d, "x", {}, false, x));
  std::vector<double> tv;
  CK(amat(c, d, dc, tv));
  for (double v : tv) {
    if (v != 0 && v != 1) return ckstop(c, "column `%s` must hold 0 or 1", dc[0].c_str());
  }
  std::vector<int> a, b;
  Val cuts, lc;
  CK(cemcode(c, d, x, "cut", 0, a, cuts));
  CK(cemcode(c, d, x, "lcut", 1, b, lc));
  const std::size_t n = nr(d);
  std::vector<int> tr(n), s(n);
  for (std::size_t i = 0; i < n; ++i) tr[i] = static_cast<int>(tv[i]);
  std::vector<double> w(n), st(6);
  const int k = static_cast<int>(x.size());
  const int rc = clab_cemtch(a.data(), static_cast<int>(n), k, tr.data(), b.data(), k, s.data(), w.data(), st.data());
  if (rc != CLAB_OK) return kerr(c, "cemtch", rc);
  Val keep = varr(), sv = varr();
  double tab[6] = {0, 0, 0, 0, 0, 0};
  for (std::size_t i = 0; i < n; ++i) {
    Val e;
    e.t = Val::BOOL;
    e.b = s[i] > 0;
    keep.a.push_back(e);
    sv.a.push_back(vnum(s[i] > 0 ? s[i] : NA));
    const int g = tr[i];
    tab[g * 3] += 1;
    tab[g * 3 + (s[i] > 0 ? 1 : 2)] += 1;
  }
  Val l1 = vobj(), ns = vobj();
  put(l1, "before", vnum(st[0]));
  put(l1, "after", vnum(st[1]));
  put(ns, "kept", vnum(st[2]));
  put(ns, "all", vnum(st[5]));
  out = vobj();
  put(out, "keep", keep);
  put(out, "strata", sv);
  put(out, "w", vnums(w));
  put(out, "tab", vmat(tab, 3, 2));
  put(out, "l1", l1);
  put(out, "nstrata", ns);
  put(out, "cut", cuts);
  return true;
}

}  // namespace ax
