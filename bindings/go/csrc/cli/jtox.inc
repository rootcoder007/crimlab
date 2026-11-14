// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Forensic toxicology (src/btox.cpp): calibration lines, values below the
// limit of detection, postmortem redistribution, ethanol after death and
// the reliability of specimen matrices.

// One NA, the default of an argument that R gives as NA.
Val nas() {
  Val v = varr();
  v.a.push_back(Val());
  return v;
}
const Val &argna(Ctx &c, const char *nm, Val &d, double x) {
  if (c.has(nm)) return c[nm];
  d = std::isnan(x) ? nas() : vnum(x);
  return d;
}

bool txcal(Ctx &c, Val &out) {
  Dv xv, yv, kv;
  CK(num(c, c["conc"], "conc", xv, -1, 0.0, true));
  const std::size_t n = xv.n();
  CK(num(c, c["resp"], "resp", yv, static_cast<long>(n)));
  if (n < 3) return c.fail("`conc` and `resp` must give at least three calibrators.");
  Val dk = varr();
  dk.a = {vnum(3.3), vnum(10)};
  CK(num(c, c.has("k") ? c["k"] : dk, "k", kv, 2, 0.0, true));
  std::vector<double> wt(n);
  const char *wl = "numeric";
  const Val &w = c["w"];
  if (!c.has("w") || w.t == Val::STR) {
    static const char *const wo[] = {"1/x^2", "1/x", "none"};
    int m;
    CK(pickd(c, "w", wo, 3, 0, m));
    wl = wo[m];
    for (std::size_t i = 0; i < n; ++i) {
      const double x = xv[i];
      wt[i] = m == 0 ? 1 / (x * x) : m == 1 ? 1 / x : 1.0;
    }
  } else {
    Dv wv;
    CK(num(c, w, "w", wv, static_cast<long>(n), 0.0, true));
    wt = wv.v;
  }
  std::vector<double> uv;
  if (c.has("unk")) CK(nna(c, c["unk"], "unk", -1, uv));
  const int nu = static_cast<int>(uv.size());
  double st[6];
  std::vector<double> xh(uv.size() + 1);
  std::vector<int> fl(uv.size() + 1);
  CK(kc(c, clab_txcal(xv.p(), yv.p(), wt.data(), static_cast<int>(n), nu ? uv.data() : nullptr, nu, kv[0], kv[1], st,
                      xh.data(), fl.data()),
        "the calibrators give no line: `conc` does not vary or the slope is 0."));
  static const char *const fn[] = {"int", "slope", "r2", "sigma", "lod", "loq"};
  for (double &v : st) v = na(v);
  out = vobj();
  put(out, "fit", nvec(6, fn, st));
  if (nu) {
    Frame f;
    f.num("resp", uv);
    for (int j = 0; j < nu; ++j) xh[static_cast<std::size_t>(j)] = na(xh[static_cast<std::size_t>(j)]);
    f.num("conc", xh.data(), uv.size());
    static const char *const lv[] = {"below_lod", "below_loq", "quantifiable"};
    f.add("flag", fac(fl.data(), uv.size(), lv));
    put(out, "unk", f.o);
  } else {
    put(out, "unk", Val());
  }
  put(out, "w", vstr(wl));
  return true;
}

bool lodimp(Ctx &c, Val &out) {
  std::vector<double> xv;
  CK(nna(c, c["x"], "x", -1, xv));
  const std::size_t n = xv.size();
  Dv lv;
  CK(num(c, c["lod"], "lod", lv, -1, 0.0, true));
  if (lv.n() != 1 && lv.n() != n) return c.fail("`lod` must be one limit, or one per value, %d.", static_cast<int>(n));
  static const char *const mo[] = {"half", "sqrt2", "lod"};
  int m;
  CK(pickd(c, "meth", mo, 3, 0, m));
  std::vector<double> l(n), o(n);
  std::vector<int> cs(n);
  for (std::size_t i = 0; i < n; ++i) l[i] = lv[lv.n() == 1 ? 0 : i];
  double st[2];
  clab_lodimp(xv.data(), static_cast<int>(n), lv.p(), static_cast<int>(lv.n()), m, o.data(), cs.data(), st);
  Frame f;
  f.num("x", xv);
  f.num("lod", l);
  f.num("imp", o);
  f.lgl("cens", cs.data(), n);
  out = vobj();
  put(out, "tab", f.o);
  put(out, "frac", vnum(na(st[0])));
  put(out, "nmiss", vnum(st[1]));
  return true;
}

bool pmrrat(Ctx &c, Val &out) {
  Dv cv, pv, kv;
  CK(num(c, c["cen"], "cen", cv, -1, 0.0, true));
  CK(num(c, c["per"], "per", pv, static_cast<long>(cv.n()), 0.0, true));
  Val dk = varr();
  dk.a = {vnum(1), vnum(2)};
  CK(num(c, c.has("cut") ? c["cut"] : dk, "cut", kv, 2, 0.0, true));
  if (kv[1] < kv[0]) return c.fail("`cut` must be two increasing ratios.");
  const std::size_t n = cv.n();
  std::vector<double> r(n);
  std::vector<int> k(n);
  clab_pmrrat(cv.p(), pv.p(), static_cast<int>(n), kv[0], kv[1], r.data(), k.data());
  static const char *const lv[] = {"minimal", "modest", "significant"};
  Frame f;
  f.num("cen", cv.v);
  f.num("per", pv.v);
  f.num("ratio", r);
  f.add("redist", fac(k.data(), n, lv));
  out = f.o;
  return true;
}

bool etohlr(Ctx &c, Val &out) {
  Dv ev;
  CK(num(c, c["eth"], "eth", ev, -1, 0.0));
  const long n = static_cast<long>(ev.n());
  std::vector<double> pv, bv, gv;
  Val d1, d2, d3;
  CK(nna(c, argna(c, "prop", d1, 0), "prop", n, pv));
  CK(nna(c, argna(c, "but", d2, 0), "but", n, bv));
  CK(nna(c, argna(c, "etg", d3, NaN), "etg", n, gv));
  double th;
  CK(binrd(c, "thr", 0, 0.0, HUGE_VAL, th));
  std::vector<double> hi(ev.n());
  std::vector<int> v(ev.n());
  clab_etohlr(ev.p(), pv.data(), bv.data(), gv.data(), static_cast<int>(n), th, hi.data(), v.data());
  for (double &e : hi) e = na(e);
  static const char *const lv[] = {"antemortem", "postmortem", "indeterminate"};
  Frame f;
  f.num("eth", ev.v);
  f.num("prop", pv);
  f.num("but", bv);
  f.num("etg", gv);
  f.num("higher", hi);
  f.add("verdict", fac(v.data(), ev.n(), lv));
  out = f.o;
  return true;
}

bool mtxrel(Ctx &c, Val &out) {
  Dv bv;
  CK(num(c, c["base"], "base", bv, -1, 0.0));
  const std::size_t n = bv.n();
  std::vector<int> dl;
  CK(lgl(c, c["dil"], "dil", static_cast<long>(n), dl));
  double dy, sg, dm, sm, dc, sc;
  CK(binrd(c, "days", 0, 0.0, HUGE_VAL, dy));
  CK(binrd(c, "stage", 0, 0.0, HUGE_VAL, sg));
  CK(oned(c, "dmax", 30, 0.0, HUGE_VAL, dm));
  CK(oned(c, "smax", 10, 0.0, HUGE_VAL, sm));
  CK(p01d(c, "dcap", 0.8, dc));
  CK(p01d(c, "scap", 0.6, sc));
  std::vector<double> rel(n), b(n), r(n);
  std::vector<int> ord(n), d(n), rk(n);
  clab_mtxrel(bv.p(), dl.data(), static_cast<int>(n), dy, sg, dm, dc, sm, sc, rel.data(), ord.data());
  const std::vector<std::string> lab = nmsof(bv.nm, n);
  std::vector<std::string> mx(n);
  for (std::size_t j = 0; j < n; ++j) {
    const std::size_t i = static_cast<std::size_t>(ord[j]);
    mx[j] = lab[i];
    b[j] = bv[i];
    r[j] = rel[i];
    d[j] = dl[i];
    rk[j] = static_cast<int>(j) + 1;
  }
  Frame f;
  f.str("matrix", mx);
  f.num("base", b);
  f.num("rel", r);
  f.lgl("dil", d.data(), n);
  f.ints("rank", rk.data(), n);
  out = f.o;
  return true;
}
