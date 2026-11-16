// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Environmental crime and dynamic panels for clab_jrun(): lmdix(),
// fepnl() and pcdtst() (src/benv.cpp) with their checks in their order
// and words, and gmmpnl(), wsctst() and ccepnl(), whose work after the
// arguments is shared with the R binding (src/dpnd.cpp). A matrix
// argument is an array of rows or an object of numeric columns (a data
// frame, which R turns into a matrix first); its results are written as
// dev/clref.R writes R's.

namespace env {

// pevr::readm() after as.matrix(): the columns of x, n rows, finite and at
// least lo, named by the object's keys or pre1, pre2, ...
bool rdm(Ctx &c, const Val &o, const char *an, const char *pre, int n, double lo, dpnd::xm &m) {
  m = dpnd::xm();
  const char *bad = "`%s` must be a numeric matrix or vector, one row per row of `y`.";
  std::vector<const Val *> cl;
  if (o.t == Val::OBJ) {
    for (std::size_t j = 0; j < o.a.size(); ++j) {
      if (!isnum(o.a[j])) return c.fail(bad, an);
      cl.push_back(&o.a[j]);
      m.nm.push_back(o.k[j]);
    }
    for (const Val *e : cl) {
      if (vlen(*e) != zz(n)) return c.fail("`%s` must have one row per row of `y`, %d.", an, n);
    }
  }
  Val flat = varr();
  if (o.t == Val::OBJ) {
    m.p = static_cast<int>(cl.size());
    if (m.p < 1) return c.fail("`%s` must have one row per row of `y`, %d.", an, n);
    for (const Val *e : cl) {
      const Val a = ax::asarr(*e);
      flat.a.insert(flat.a.end(), a.a.begin(), a.a.end());
    }
  } else if (o.t == Val::ARR && !o.a.empty() && o.a[0].t == Val::ARR) {
    const std::size_t p = o.a[0].a.size();
    for (const Val &r : o.a) {
      if (r.t != Val::ARR || r.a.size() != p || !isnum(r)) return c.fail(bad, an);
    }
    if (o.a.size() != zz(n) || p < 1) return c.fail("`%s` must have one row per row of `y`, %d.", an, n);
    m.p = static_cast<int>(p);
    for (std::size_t j = 0; j < p; ++j) {
      for (const Val &r : o.a) flat.a.push_back(r.a[j]);
    }
  } else {
    if (!isnum(o)) return c.fail(bad, an);
    flat = ax::asarr(o);
    if (flat.a.size() != zz(n)) return c.fail("`%s` must have one row per row of `y`, %d.", an, n);
    m.p = 1;
  }
  if (m.nm.empty()) {
    for (int j = 0; j < m.p; ++j) m.nm.push_back(pre + std::to_string(j + 1));
  }
  Dv d;
  CK(num(c, flat, an, d, static_cast<long>(n) * m.p, lo));
  m.v = d.v;
  return true;
}

// The panel of dpnd from the one of jdid.inc.
dpnd::pnl plain(const Pan &p) {
  dpnd::pnl q;
  q.n = p.n;
  q.y = p.y.data();
  q.u = p.u;
  q.t = p.t;
  q.nu = p.nu;
  q.nt = p.nt;
  q.tv = p.tv;
  q.ulev = p.ulev;
  return q;
}

// Names: one string or an array of strings, none missing (character(0)
// is an empty array); false, with the message, otherwise.
bool names(Ctx &c, const Val &o, const char *an, const char *what, std::vector<std::string> &s) {
  s.clear();
  if (o.t == Val::STR) {
    s.push_back(o.s);
    return true;
  }
  if (o.t != Val::ARR) return c.fail("`%s` must be %s.", an, what);
  for (const Val &e : o.a) {
    if (e.t == Val::NUL) return c.fail("`%s` must not hold NA.", an);
    if (e.t != Val::STR) return c.fail("`%s` must be %s.", an, what);
    s.push_back(e.s);
  }
  return true;
}

// Columns of x named, each once (bdpn.cpp, cnames()).
bool cnames(Ctx &c, const Val &o, const char *an, const dpnd::xm &x, std::vector<int> &out) {
  std::vector<std::string> s;
  CK(names(c, o, an, "names of columns of `x`", s));
  out.clear();
  for (const std::string &v : s) {
    const auto it = std::find(x.nm.begin(), x.nm.end(), v);
    if (it == x.nm.end()) return c.fail("`%s` names %s, which is not a column of `x`.", an, v.c_str());
    const int j = static_cast<int>(it - x.nm.begin());
    if (std::find(out.begin(), out.end(), j) != out.end()) return c.fail("`%s` names %s twice.", an, v.c_str());
    out.push_back(j);
  }
  return true;
}

// Lag limits (bdpn.cpp, lags()).
bool lags(Ctx &c, const Val &o, const char *an, int nv, int lo, bool inf, std::vector<int> &out) {
  Dv d;
  const Got g = readd(o, d);
  if (g != NUM || (d.n() != 1 && d.n() != zz(nv))) {
    return c.fail("`%s` must be one number or %d, one per GMM-style variable.", an, nv);
  }
  out.clear();
  for (int i = 0; i < nv; ++i) {
    const double v = d[d.n() == 1 ? 0 : zz(i)];
    if (inf && v == Inf) {
      out.push_back(1 << 30);
      continue;
    }
    if (!std::isfinite(v) || v != std::floor(v) || v < lo || v > 1e6) {
      return c.fail("`%s` must hold whole numbers from %d%s.", an, lo, inf ? ", or Inf" : "");
    }
    out.push_back(static_cast<int>(v));
  }
  return true;
}

Val strs1(const std::vector<std::string> &v) { return unbox(vstrs(v)); }

Val nas(const std::vector<double> &v) {
  Val a = varr();
  for (double x : v) a.a.push_back(vnum(na(x)));
  return a;
}

Val ctab(const dpnd::tab &t, bool tst) {
  Frame f;
  f.str("term", t.name);
  f.add("est", nas(t.est));
  f.add("se", nas(t.se));
  f.add(tst ? "t" : "z", nas(t.stat));
  f.add("p", nas(t.p));
  f.add("lo", nas(t.lo));
  f.add("hi", nas(t.hi));
  return f.o;
}

bool done(Ctx &c, const std::string &msg, int code) {
  if (code == CLAB_EINTR) return c.fail("%s", clab_strerror(CLAB_EINTR));
  if (!msg.empty()) return c.fail("%s", msg.c_str());
  return true;
}

}  // namespace env

bool gmmpnl(Ctx &c, Val &out) {
  Pan pn;
  CK(panel(c, -HUGE_VAL, pn));
  dpnd::xm xm, zx;
  if (c.has("x")) CK(env::rdm(c, c["x"], "x", "x", pn.n, -HUGE_VAL, xm));
  if (c.has("z")) CK(env::rdm(c, c["z"], "z", "z", pn.n, -HUGE_VAL, zx));
  dpnd::spec sp;
  static const char *const tyo[] = {"system", "diff"};
  static const char *const tro[] = {"fd", "fod"};
  static const char *const hmo[] = {"ab", "dpd", "id"};
  CK(wholed(c, "ylag", 1, 0, sp.ylag));
  if (sp.ylag > 1000) return c.fail("`ylag` must be at most 1000.");
  if (sp.ylag == 0 && xm.p == 0) return c.fail("the model needs `ylag` above 0 or columns in `x`.");
  int ty = 0;
  CK(pickd(c, "type", tyo, 2, 0, ty));
  sp.sys = ty == 0;
  CK(pickd(c, "trans", tro, 2, 0, sp.fod));
  CK(pickd(c, "hmat", hmo, 3, 0, sp.hm));
  bool b = false;
  CK(flagd(c, "collapse", false, b));
  sp.coll = b;
  CK(flagd(c, "cons", true, b));
  sp.cons = b;
  CK(flagd(c, "tdum", true, b));
  sp.tdum = b;
  std::vector<int> gx;
  if (c.has("gmm")) CK(env::cnames(c, c["gmm"], "gmm", xm, gx));
  if (sp.ylag > 0) sp.gv.push_back(-1);
  sp.gv.insert(sp.gv.end(), gx.begin(), gx.end());
  if (!c.has("iv")) {
    for (int j = 0; j < xm.p; ++j) {
      if (std::find(gx.begin(), gx.end(), j) == gx.end()) sp.iv.push_back(j);
    }
  } else {
    CK(env::cnames(c, c["iv"], "iv", xm, sp.iv));
    for (int j : sp.iv) {
      if (std::find(gx.begin(), gx.end(), j) != gx.end()) {
        return c.fail("`iv` and `gmm` both name %s.", xm.nm[zz(j)].c_str());
      }
    }
  }
  const int nv = static_cast<int>(sp.gv.size());
  if (nv > 0) {
    Val d1, d2;
    CK(env::lags(c, arg(c, "lagmin", d1, 2), "lagmin", nv, 1, false, sp.lo));
    CK(env::lags(c, arg(c, "lagmax", d2, Inf), "lagmax", nv, 1, true, sp.hi));
    for (int i = 0; i < nv; ++i) {
      if (sp.hi[zz(i)] < sp.lo[zz(i)]) return c.fail("`lagmax` must be at least `lagmin`.");
    }
  }
  CK(wholed(c, "steps", 1, 1, sp.steps));
  if (sp.steps > 2) return c.fail("`steps` must be 1 or 2.");
  CK(flagd(c, "robust", true, b));
  sp.rob = b;
  if (c.has("wc")) {
    CK(flag(c, c["wc"], "wc", b));
    sp.wc = b;
  }
  CK(env::done(c, dpnd::gcheck(sp), CLAB_OK));
  CK(flagd(c, "small", false, b));
  sp.small = b;
  CK(wholed(c, "ar", 2, 0, sp.nar));
  if (sp.nar > 100) return c.fail("`ar` must be at most 100.");
  CK(oned(c, "level", 0.95, 0.0, 1.0, sp.level));
  dpnd::gout o;
  int code = CLAB_OK;
  const std::string msg = dpnd::gfit(env::plain(pn), xm, zx, sp, o, code);
  CK(env::done(c, msg, code));
  const int k = static_cast<int>(o.k);
  out = vobj();
  put(out, "coef", env::ctab(o.coef, o.tstat));
  Frame te;
  te.str("test", o.tests.name);
  te.add("stat", env::nas(o.tests.stat));
  te.add("df", env::nas(o.tests.df));
  te.add("p", env::nas(o.tests.p));
  put(out, "tests", te.o);
  Frame wd;
  wd.str("test", o.wald.name);
  wd.add("stat", env::nas(o.wald.stat));
  wd.add("df1", env::nas(o.wald.df));
  wd.add("df2", env::nas(o.wald.df2));
  wd.add("p", env::nas(o.wald.p));
  put(out, "wald", wd.o);
  Frame fs;
  col1(fs, "n", o.n);
  col1(fs, "groups", o.groups);
  col1(fs, "inst", o.inst);
  col1(fs, "omin", o.omin);
  col1(fs, "oavg", o.oavg);
  col1(fs, "omax", o.omax);
  col1(fs, "steps", o.steps);
  col1(fs, "k", o.k);
  col1(fs, "sigma2", o.sigma2);
  put(out, "fit", fs.o);
  put(out, "vcov", vmat(o.vcov, k, k));
  put(out, "inst", env::strs1(o.inst_names));
  put(out, "dropped", env::strs1(o.dropped));
  put(out, "dummies", env::strs1(o.dummies));
  put(out, "note", env::strs1(o.note));
  put(out, "type", vstr(o.type));
  put(out, "trans", vstr(o.trans));
  put(out, "robust", vbool(sp.rob != 0));
  put(out, "small", vbool(sp.small != 0));
  put(out, "level", vnum(sp.level));
  return true;
}

bool wsctst(Ctx &c, Val &out) {
  static const char *const sso[] = {"cr1", "none"};
  Pan pn;
  CK(panel(c, -HUGE_VAL, pn));
  dpnd::xm xm;
  if (c.has("x")) CK(env::rdm(c, c["x"], "x", "x", pn.n, -HUGE_VAL, xm));
  bool cn = false;
  CK(flagd(c, "cons", false, cn));
  int raw = 0;
  CK(pickd(c, "ssc", sso, 2, 0, raw));
  dpnd::wout o;
  int code = CLAB_OK;
  const std::string msg = dpnd::wfit(env::plain(pn), xm, cn, raw, o, code);
  CK(env::done(c, msg, code));
  static const char *const sn[] = {"F", "df1", "df2", "p", "coef", "se", "pairs", "clusters", "ndiff"};
  out = vobj();
  put(out, "stat", nvec(9, sn, o.stat));
  Val co;
  if (xm.p > 0) {
    Frame f;
    f.str("term", xm.nm);
    f.add("est", env::nas(o.b));
    f.add("se", env::nas(o.se));
    co = f.o;
  }
  put(out, "coef", co);
  put(out, "ssc", vstr(sso[raw]));
  return true;
}

bool ccepnl(Ctx &c, Val &out) {
  Pan pn;
  CK(panel(c, -HUGE_VAL, pn));
  const int n = pn.n;
  dpnd::xm xm, fm;
  CK(env::rdm(c, c["x"], "x", "x", n, -HUGE_VAL, xm));
  if (c.has("fac")) CK(env::rdm(c, c["fac"], "fac", "f", n, -HUGE_VAL, fm));
  bool cn = true;
  CK(flagd(c, "cons", true, cn));
  std::vector<std::string> an;
  if (c.has("avg")) CK(env::names(c, c["avg"], "avg", "names: y or columns of `x`", an));
  double lv = 0.0;
  CK(oned(c, "level", 0.95, 0.0, 1.0, lv));
  std::vector<double> wv;
  if (c.has("w")) {
    CK(isvec(c, c["w"], "w"));
    Dv w;
    CK(num(c, c["w"], "w", w, n, 0.0, true));
    wv = w.v;
  }
  dpnd::cres o;
  int code = CLAB_OK;
  const std::string msg = dpnd::cfit(env::plain(pn), xm, fm, cn, wv, !c.has("avg"), an, lv, o, code);
  CK(env::done(c, msg, code));
  const int k = xm.p, nu = pn.nu;
  out = vobj();
  put(out, "mg", env::ctab(o.mg, false));
  put(out, "pooled", env::ctab(o.pooled, false));
  put(out, "bi", vmat(o.bi, nu, k));
  put(out, "vmg", vmat(o.vmg, k, k));
  put(out, "vp", vmat(o.vp, k, k));
  Frame fs;
  col1(fs, "units", o.units);
  col1(fs, "periods", o.periods);
  col1(fs, "k", o.k);
  col1(fs, "h", o.h);
  col1(fs, "n", o.n);
  put(out, "fit", fs.o);
  Val w = vobj();
  for (int a = 0; a < nu; ++a) put(w, pn.ulev[zz(a)].c_str(), vnum(o.w[zz(a)]));
  put(out, "w", w);
  put(out, "level", vnum(lv));
  return true;
}

// pcdtst() (src/benv.cpp).
bool pcdtst(Ctx &c, Val &out) {
  Pan pn;
  CK(panel(c, -HUGE_VAL, pn, &c["x"]));
  if (static_cast<double>(pn.nu) * pn.nt > 67108864.0) {
    return c.fail("the units by periods (%.0f cells) are too many.", static_cast<double>(pn.nu) * pn.nt);
  }
  std::vector<double> e(zz(pn.nu) * zz(pn.nt), NaN);
  for (std::size_t i = 0; i < zz(pn.n); ++i) e[zz(pn.u[i]) * zz(pn.nt) + zz(pn.t[i])] = pn.y[i];
  double o[7];
  if (clab_pcd(e.data(), pn.nu, pn.nt, o) != CLAB_OK) {
    return c.fail("no pair of units shares two periods with variation in both: CD is undefined.");
  }
  o[5] = pn.nu;
  o[6] = pn.nt;
  static const char *const nm[] = {"cd", "p", "rho", "absrho", "pairs", "units", "periods"};
  out = nvec(7, nm, o);
  return true;
}

namespace env {

// benv.cpp, periods(): NULL for def (0 the first and last, 1 all), else
// periods of the panel, increasing; their codes.
bool periods(Ctx &c, const Pan &pn, int def, int need, std::vector<int> &out) {
  out.clear();
  if (!c.has("at")) {
    if (def == 0) {
      out.push_back(0);
      out.push_back(pn.nt - 1);
    } else {
      for (int s = 0; s < pn.nt; ++s) out.push_back(s);
    }
    return true;
  }
  Dv a;
  CK(num(c, c["at"], "at", a));
  for (std::size_t q = 0; q < a.n(); ++q) {
    int s = -1;
    for (int r = 0; r < pn.nt; ++r) {
      if (pn.tv[zz(r)] == a[q]) s = r;
    }
    if (s < 0) return c.fail("`at` must hold periods of `time`.");
    if (!out.empty() && s <= out.back()) return c.fail("`at` must be increasing.");
    out.push_back(s);
  }
  if (static_cast<int>(out.size()) < need) return c.fail("`at` must hold at least %d periods.", need);
  return true;
}

// The factors of row values y (1) and q (m): y / q1, q1 / q2, ..., qm.
void factors(double y, const double *q, int m, std::size_t step, double *f) {
  f[0] = y / q[0];
  for (int j = 1; j < m; ++j) f[zz(j)] = q[zz(j - 1) * step] / q[zz(j) * step];
  f[zz(m)] = q[zz(m - 1) * step];
}

struct dec {
  std::vector<int> unit, s0, s1;
  std::vector<double> y0, y1, dy, f0, f1, c;
};

inline bool good(double v) { return std::isfinite(v) && v > 0; }

void decompose(dec &d, int k) {
  const int r = static_cast<int>(d.unit.size());
  d.c.assign(zz(r) * zz(k), NA);
  std::vector<int> ok;
  for (int i = 0; i < r; ++i) {
    bool g = good(d.y0[zz(i)]) && good(d.y1[zz(i)]);
    for (int j = 0; j < k; ++j) g = g && good(d.f0[zz(i) * zz(k) + zz(j)]) && good(d.f1[zz(i) * zz(k) + zz(j)]);
    if (g) ok.push_back(i);
  }
  const int m = static_cast<int>(ok.size());
  if (m == 0) return;
  std::vector<double> y0(zz(m)), y1(zz(m)), a(zz(m) * zz(k)), b(zz(m) * zz(k)), cc(zz(m) * zz(k));
  for (int h = 0; h < m; ++h) {
    const std::size_t i = zz(ok[zz(h)]);
    y0[zz(h)] = d.y0[i];
    y1[zz(h)] = d.y1[i];
    for (int j = 0; j < k; ++j) {
      a[zz(j) * zz(m) + zz(h)] = d.f0[i * zz(k) + zz(j)];
      b[zz(j) * zz(m) + zz(h)] = d.f1[i * zz(k) + zz(j)];
    }
  }
  clab_lmdi(y0.data(), y1.data(), a.data(), b.data(), m, k, cc.data());
  for (int h = 0; h < m; ++h) {
    for (int j = 0; j < k; ++j) d.c[zz(j) * zz(r) + zz(ok[zz(h)])] = cc[zz(j) * zz(m) + zz(h)];
  }
}

Val decframe(const dec &d, const Pan &pn, const std::vector<std::string> &fn, bool sp) {
  const int r = static_cast<int>(d.unit.size()), k = static_cast<int>(fn.size());
  std::vector<std::string> u;
  std::vector<double> a, b, v0, v1, dy;
  for (int i = 0; i < r; ++i) {
    u.push_back(pn.ulev[zz(d.unit[zz(i)])]);
    a.push_back(pn.tv[zz(d.s0[zz(i)])]);
    b.push_back(pn.tv[zz(d.s1[zz(i)])]);
    v0.push_back(d.y0[zz(i)]);
    v1.push_back(d.y1[zz(i)]);
    dy.push_back(d.dy[zz(i)]);
  }
  Frame f;
  f.str("unit", u);
  f.num(sp ? "time" : "from", a);
  if (!sp) f.num("to", b);
  f.num(sp ? "ref" : "y0", v0);
  f.num(sp ? "y" : "y1", v1);
  f.num("dy", dy);
  for (int j = 0; j < k; ++j) f.num(fn[zz(j)].c_str(), d.c.data() + zz(j) * zz(r), zz(r));
  return f.o;
}

// benv.cpp, refer(): the reference region of a period.
void refer(const Pan &pn, const dpnd::xm &q, const std::vector<int> &rows, int gm, double zr, double *ry,
           double *rf) {
  const int m = q.p, n = pn.n;
  std::vector<double> a(zz(m) + 1, 0.0);
  for (int i : rows) {
    for (int j = 0; j <= m; ++j) {
      double v = j == 0 ? pn.y[zz(i)] : q.v[zz(j - 1) * zz(n) + zz(i)];
      if (gm) v = std::log(v > 0 ? v : zr);
      a[zz(j)] += v / static_cast<double>(rows.size());
    }
  }
  for (int j = 0; j <= m; ++j) {
    if (gm) a[zz(j)] = std::exp(a[zz(j)]);
    if (!gm && a[zz(j)] <= 0) a[zz(j)] = zr;
  }
  *ry = a[0];
  factors(a[0], a.data() + 1, m, 1, rf);
}

}  // namespace env

// lmdix() (src/benv.cpp).
bool lmdix(Ctx &c, Val &out) {
  static const char *const tyo[] = {"time", "space"};
  static const char *const rfo[] = {"mean", "gmean"};
  int sp = 0, gm = 0;
  CK(pickd(c, "type", tyo, 2, 0, sp));
  CK(pickd(c, "ref", rfo, 2, 0, gm));
  double zr = 0.0;
  CK(oned(c, "zero", 0.0, -1.0, Inf, zr));
  if (zr < 0) return c.fail("`zero` must be 0 (leave zeros out) or a positive number.");
  Pan pn;
  CK(panel(c, 0.0, pn));
  const int n = pn.n;
  dpnd::xm qm;
  CK(env::rdm(c, c["q"], "q", "q", n, 0.0, qm));
  const int m = qm.p, k = m + 1;
  std::vector<std::string> fn;
  if (c.has("fnm")) {
    const Val &f = c["fnm"];
    const Val a = ax::asarr(f);
    if (a.a.size() != zz(k)) return c.fail("`fnm` must be %d names, one per factor.", k);
    for (const Val &e : a.a) {
      if (e.t == Val::NUL) return c.fail("`fnm` must not hold NA.");
      if (e.t != Val::STR) return c.fail("`fnm` must be %d names, one per factor.", k);
      fn.push_back(e.s);
    }
  } else {
    fn.push_back("y/" + qm.nm[0]);
    for (int j = 1; j < m; ++j) fn.push_back(qm.nm[zz(j - 1)] + "/" + qm.nm[zz(j)]);
    fn.push_back(qm.nm[zz(m - 1)]);
  }
  std::vector<double> yz(zz(n)), fz(zz(n) * zz(k)), qz(qm.v);
  for (double &v : qz) v = v > 0 ? v : zr;
  for (int i = 0; i < n; ++i) {
    yz[zz(i)] = pn.y[zz(i)] > 0 ? pn.y[zz(i)] : zr;
    env::factors(yz[zz(i)], qz.data() + i, m, zz(n), fz.data() + zz(i) * zz(k));
  }
  const double cells = static_cast<double>(pn.nu) * pn.nt;
  if (cells > 67108864.0) return c.fail("the units by periods (%.0f cells) are too many.", cells);
  std::vector<int> row(zz(pn.nu) * zz(pn.nt), -1);
  for (int i = 0; i < n; ++i) row[zz(pn.u[zz(i)]) * zz(pn.nt) + zz(pn.t[zz(i)])] = i;
  std::vector<int> ps;
  CK(env::periods(c, pn, sp, sp ? 1 : 2, ps));
  env::dec d;
  Val rf;
  auto add = [&](int a, int s0, int s1, double y0, double y1, double dy, const double *f0, const double *f1) {
    d.unit.push_back(a);
    d.s0.push_back(s0);
    d.s1.push_back(s1);
    d.y0.push_back(y0);
    d.y1.push_back(y1);
    d.dy.push_back(dy);
    d.f0.insert(d.f0.end(), f0, f0 + k);
    d.f1.insert(d.f1.end(), f1, f1 + k);
  };
  if (!sp) {
    for (std::size_t h = 1; h < ps.size(); ++h) {
      for (int a = 0; a < pn.nu; ++a) {
        const int i0 = row[zz(a) * zz(pn.nt) + zz(ps[h - 1])], i1 = row[zz(a) * zz(pn.nt) + zz(ps[h])];
        if (i0 < 0 || i1 < 0) continue;
        add(a, ps[h - 1], ps[h], pn.y[zz(i0)], pn.y[zz(i1)], pn.y[zz(i1)] - pn.y[zz(i0)], fz.data() + zz(i0) * zz(k),
            fz.data() + zz(i1) * zz(k));
      }
    }
  } else {
    const int nr = static_cast<int>(ps.size());
    std::vector<double> rt(zz(nr)), ry(zz(nr)), rc(zz(nr) * zz(k)), f(zz(k));
    for (int h = 0; h < nr; ++h) {
      std::vector<int> rows;
      for (int a = 0; a < pn.nu; ++a) {
        const int i = row[zz(a) * zz(pn.nt) + zz(ps[zz(h)])];
        if (i >= 0) rows.push_back(i);
      }
      double yr = 0;
      env::refer(pn, qm, rows, gm, zr, &yr, f.data());
      rt[zz(h)] = pn.tv[zz(ps[zz(h)])];
      ry[zz(h)] = yr;
      for (int j = 0; j < k; ++j) rc[zz(j) * zz(nr) + zz(h)] = f[zz(j)];
      for (int i : rows) {
        add(pn.u[zz(i)], ps[zz(h)], ps[zz(h)], yr, pn.y[zz(i)], pn.y[zz(i)] - yr, f.data(), fz.data() + zz(i) * zz(k));
      }
    }
    Frame fr;
    fr.num("time", rt);
    fr.num("y", ry);
    for (int j = 0; j < k; ++j) fr.num(fn[zz(j)].c_str(), rc.data() + zz(j) * zz(nr), zz(nr));
    rf = fr.o;
  }
  std::vector<double> y0(d.y0), y1(d.y1);
  for (double &v : d.y0) v = v > 0 ? v : zr;
  for (double &v : d.y1) v = v > 0 ? v : zr;
  env::decompose(d, k);
  d.y0.swap(y0);
  d.y1.swap(y1);
  // the factors of every row as given (zeros kept)
  std::vector<std::string> fu;
  std::vector<double> ft, fy, fv(zz(n) * zz(k)), f(zz(k));
  for (int i = 0; i < n; ++i) {
    fu.push_back(pn.ulev[zz(pn.u[zz(i)])]);
    ft.push_back(pn.tv[zz(pn.t[zz(i)])]);
    fy.push_back(pn.y[zz(i)]);
    env::factors(pn.y[zz(i)], qm.v.data() + i, m, zz(n), f.data());
    for (int j = 0; j < k; ++j) fv[zz(j) * zz(n) + zz(i)] = na(f[zz(j)]);
  }
  Frame fa;
  fa.str("unit", fu);
  fa.num("time", ft);
  fa.num("y", fy);
  for (int j = 0; j < k; ++j) fa.num(fn[zz(j)].c_str(), fv.data() + zz(j) * zz(n), zz(n));
  out = vobj();
  put(out, "dec", env::decframe(d, pn, fn, sp != 0));
  put(out, "fac", fa.o);
  put(out, "ref", rf);
  put(out, "factors", vstrs(fn));
  put(out, "type", vstr(tyo[sp]));
  put(out, "refkind", vstr(rfo[gm]));
  return true;
}

namespace env {

double sdev(const double *v, int n) {
  double m = 0, s = 0;
  for (int i = 0; i < n; ++i) m += v[i] / n;
  for (int i = 0; i < n; ++i) s += (v[i] - m) * (v[i] - m);
  return std::sqrt(s / (n - 1));
}

double mean(const double *v, int n) {
  double m = 0;
  for (int i = 0; i < n; ++i) m += v[i] / n;
  return m;
}

}  // namespace env

// fepnl() (src/benv.cpp).
bool fepnl(Ctx &c, Val &out) {
  static const char *const sso[] = {"cr1", "none"};
  int raw = 0;
  CK(pickd(c, "ssc", sso, 2, 0, raw));
  double lv = 0.0;
  CK(oned(c, "level", 0.95, 0.0, 1.0, lv));
  Pan pn;
  CK(panel(c, -HUGE_VAL, pn));
  const int n = pn.n;
  dpnd::xm xm;
  CK(env::rdm(c, c["x"], "x", "x", n, -HUGE_VAL, xm));
  int ia = -1, ib = -1;
  if (c.has("int")) {
    const Val a = ax::asarr(c["int"]);
    if (a.a.size() != 2 || a.a[0].t != Val::STR || a.a[1].t != Val::STR) {
      return c.fail("`int` must be two names of columns of `x`.");
    }
    for (int j = 0; j < 2; ++j) {
      const int h = find(xm.nm, a.a[zz(j)].s.c_str());
      if (h < 0) return c.fail("`int` names %s, which is not a column of `x`.", a.a[zz(j)].s.c_str());
      (j ? ib : ia) = h;
    }
    if (ia == ib) return c.fail("`int` must name two different columns of `x`.");
    for (int i = 0; i < n; ++i) xm.v.push_back(xm.v[zz(ia) * zz(n) + zz(i)] * xm.v[zz(ib) * zz(n) + zz(i)]);
    xm.nm.push_back(xm.nm[zz(ia)] + ":" + xm.nm[zz(ib)]);
    ++xm.p;
  } else if (c.has("at")) {
    return c.fail("`at` gives moderator values for the simple slopes: give `int` too.");
  }
  const int p = xm.p;
  std::vector<int> cl;
  int ncl = pn.nu;
  if (!c.has("clus")) {
    cl = pn.u;
  } else {
    Lab l;
    CK(codes(c, c["clus"], "clus", zz(n), l));
    cl = l.code;
    ncl = l.g;
    if (ncl < 2) return c.fail("`clus` has one cluster: a cluster-robust variance needs two or more.");
  }
  if (ncl < 10) {
    c.warn.push_back(std::string(c.fn) +
                     ": fewer than ten clusters; cluster-robust errors are unreliable with so few (wcbx() is made "
                     "for this case).");
  }
  const std::size_t sp = zz(p), sn = zz(n);
  std::vector<double> b(sp), vc(sp * sp),
      work(sn * (sp + 2) + sp * sp + 2 * zz(pn.nu + pn.nt) + (sn + 2 * sp) * sp + sp);
  std::vector<int> iwork(sn + zz(std::max(pn.nu, std::max(pn.nt, ncl))));
  double st[5];
  if (n <= p + pn.nu + pn.nt ||
      clab_didfe(pn.y.data(), xm.v.data(), n, p, pn.u.data(), pn.nu, pn.t.data(), pn.nt, 1, nullptr, cl.data(), ncl,
                 nullptr, 0, 1e-10, 10000, b.data(), vc.data(), st, work.data(), iwork.data()) != CLAB_OK) {
    return c.fail("a column of `x` is absorbed by the unit and period effects, the columns are collinear, or there "
                  "are too few rows for the effects.");
  }
  if (raw) {
    for (double &v : vc) v /= st[3];
  }
  const double *yw = work.data(), *res = work.data() + sn * (sp + 1);
  std::vector<double> vif(sp), wv(sn * (sp + 1) + 2 * sp * (sp + 1));
  if (clab_vif(xm.v.data(), n, p, vif.data(), wv.data()) != CLAB_OK) vif.assign(sp, NA);
  double jb[4];
  const bool jok = clab_jb(res, n, jb) == CLAB_OK;
  const double df = st[1] - 1, tq = clab_qt(1 - (1 - lv) / 2, df, 1);
  const double sy = env::sdev(pn.y.data(), n);
  std::vector<double> ce(sp), cs(sp), cl0(sp), ch(sp), cz(sp), cp(sp), cb(sp), cv(sp);
  for (std::size_t j = 0; j < sp; ++j) {
    ce[j] = b[j];
    cs[j] = std::sqrt(vc[j * sp + j]);
    cl0[j] = ce[j] - tq * cs[j];
    ch[j] = ce[j] + tq * cs[j];
    cz[j] = ce[j] / cs[j];
    cp[j] = 2 * clab_pt(-std::fabs(cz[j]), df, 1);
    cb[j] = ce[j] * env::sdev(xm.v.data() + j * sn, n) / sy;
    cv[j] = na(vif[j]);
  }
  Frame co;
  co.str("term", xm.nm);
  co.num("est", ce);
  co.num("se", cs);
  co.num("lo", cl0);
  co.num("hi", ch);
  co.num("t", cz);
  co.num("p", cp);
  co.num("beta", cb);
  co.num("vif", cv);
  Val sl;
  if (ia >= 0) {
    const double *mv = xm.v.data() + zz(ib) * sn;
    std::vector<double> mm;
    if (!c.has("at")) {
      const double mu = env::mean(mv, n), s = env::sdev(mv, n);
      mm = {mu - s, mu, mu + s};
    } else {
      Dv a;
      CK(num(c, c["at"], "at", a));
      mm = a.v;
    }
    const std::size_t a = zz(ia), cc = sp - 1;
    std::vector<double> oa, oe, os, ol, oh, ot, op;
    for (double x : mm) {
      const double e = b[a] + b[cc] * x;
      const double v = vc[a * sp + a] + x * x * vc[cc * sp + cc] + 2 * x * vc[a * sp + cc];
      const double s = std::sqrt(v);
      oa.push_back(x);
      oe.push_back(e);
      os.push_back(s);
      ol.push_back(e - tq * s);
      oh.push_back(e + tq * s);
      ot.push_back(e / s);
      op.push_back(2 * clab_pt(-std::fabs(e / s), df, 1));
    }
    Frame f;
    f.num("at", oa);
    f.num("slope", oe);
    f.num("se", os);
    f.num("lo", ol);
    f.num("hi", oh);
    f.num("t", ot);
    f.num("p", op);
    sl = f.o;
  }
  double tss = 0, tw = 0;
  const double my = env::mean(pn.y.data(), n);
  for (int i = 0; i < n; ++i) {
    tss += (pn.y[zz(i)] - my) * (pn.y[zz(i)] - my);
    tw += yw[i] * yw[i];
  }
  Frame fs;
  col1(fs, "n", n);
  col1(fs, "units", pn.nu);
  col1(fs, "periods", pn.nt);
  col1(fs, "ncl", st[1]);
  col1(fs, "k", st[2]);
  col1(fs, "r2", 1 - st[4] / tss);
  col1(fs, "r2w", 1 - st[4] / tw);
  col1(fs, "rss", st[4]);
  static const char *const jn[] = {"skew", "kurt", "jb", "p"};
  const double na4[4] = {NA, NA, NA, NA};
  out = vobj();
  put(out, "coef", co.o);
  put(out, "slopes", sl);
  put(out, "fit", fs.o);
  put(out, "jb", nvec(4, jn, jok ? jb : na4));
  put(out, "vcov", vmat(vc, p, p));
  put(out, "resid", vnums(res, sn));
  put(out, "ssc", vstr(sso[raw]));
  put(out, "level", vnum(lv));
  return true;
}
