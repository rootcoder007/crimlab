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
