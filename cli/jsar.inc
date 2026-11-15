// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Search theory for search and rescue (src/sar.cpp), with the checks of
// R/sar.R and R/ckx.R in their words.

// A numeric argument as .cknum() reads it: its values in R's order and,
// for a matrix, its dimensions (nr < 0 for a vector).
struct Nv {
  std::vector<double> v;
  int nr = -1, nc = -1;
  std::size_t n() const { return v.size(); }
  double operator[](std::size_t i) const { return v[i]; }
  const double *p() const { return v.data(); }
};

// R's class of a JSON value that is not numbers: logical, character, list.
inline const char *jcls(const Val &o, bool &lg) {
  lg = false;
  if (o.t == Val::BOOL) {
    lg = true;
    return "logical";
  }
  if (o.t == Val::STR) {
    double x = 0;
    return snum(o, x) ? nullptr : "character";
  }
  if (o.t == Val::ARR || o.t == Val::OBJ) {
    bool anyb = false, anys = false, anyl = false, alln = !o.a.empty();
    double x = 0;
    for (const Val &e : o.a) {
      if (e.t == Val::BOOL) anyb = true;
      else if (e.t == Val::STR && !snum(e, x)) anys = true;
      else if ((e.t == Val::ARR || e.t == Val::OBJ) && o.t == Val::OBJ) anyl = true;
      if (e.t != Val::NUL) alln = false;
    }
    if (anyl) return "list";
    if (anys) return "character";
    if (anyb) {
      lg = true;
      return "logical";
    }
    (void)alln;
  }
  return nullptr;
}

inline bool ckstop(Ctx &c, const char *f, ...) {
  std::va_list ap;
  va_start(ap, f);
  std::vector<char> b(1024);
  std::vsnprintf(b.data(), b.size(), f, ap);
  va_end(ap);
  return c.fail("%s.", b.data());
}

inline bool cknum(Ctx &c, const Val &o, const char *w, Nv &d, bool scalar = true, double mn = -HUGE_VAL,
                  double mx = HUGE_VAL, bool positive = false, bool whole = false) {
  d = Nv();
  if (o.t == Val::NUL) return ckstop(c, "`%s` must be numeric, got NULL", w);
  bool lg = false;
  const char *cl = jcls(o, lg);
  if (cl != nullptr) {
    // a logical vector of NA only is read as numbers (all missing)
    bool allna = lg && o.t != Val::BOOL;
    for (const Val &e : o.a) allna = allna && e.t == Val::NUL;
    if (lg && !allna) return ckstop(c, "`%s` must be numeric, not logical", w);
    if (!lg) return ckstop(c, "`%s` must be numeric, got %s", w, cl);
  }
  if (o.t == Val::ARR && !o.a.empty() && o.a[0].t == Val::ARR) {
    int dim[2];
    for (int &k : dim) k = -1;
    if (!shape(o, 2, dim, 0)) return ckstop(c, "`%s` must be numeric, got list", w);
    d.v.assign(static_cast<std::size_t>(dim[0]) * static_cast<std::size_t>(dim[1]), 0.0);
    bool bad = false;
    fill(o, 2, dim, 0, 0, 1, d.v, bad, true);
    d.nr = dim[0];
    d.nc = dim[1];
  } else {
    Dv t;
    readd(o, t);
    d.v = t.v;
  }
  if (d.v.empty()) return ckstop(c, "`%s` must have at least one value", w);
  if (scalar && d.v.size() != 1) return ckstop(c, "`%s` must be a single number", w);
  for (double x : d.v) {
    if (!std::isfinite(x)) return ckstop(c, "`%s` must be finite, with no missing values", w);
  }
  for (double x : d.v) {
    if (positive && x <= 0) return ckstop(c, "`%s` must be positive", w);
  }
  for (double x : d.v) {
    if (x < mn || x > mx) {
      if (std::isfinite(mx)) {
        return ckstop(c, "`%s` must be between %s and %s", w, rfmt(mn).c_str(), rfmt(mx).c_str());
      }
      return ckstop(c, "`%s` must be at least %s", w, rfmt(mn).c_str());
    }
  }
  for (double x : d.v) {
    if (whole && x != std::nearbyint(x)) return ckstop(c, "`%s` must be a whole number", w);
  }
  return true;
}

inline bool cknumd(Ctx &c, const char *nm, double def, Nv &d, double mn = -HUGE_VAL, double mx = HUGE_VAL,
                   bool positive = false, bool whole = false) {
  if (!c.has(nm)) {
    d = Nv();
    d.v.push_back(def);
    return true;
  }
  return cknum(c, c[nm], nm, d, true, mn, mx, positive, whole);
}

// deparse() of a short value, for the refusal of .ckarg().
inline std::string jdep(const Val &o) {
  if (o.t == Val::STR) return "\"" + o.s + "\"";
  if (o.t == Val::NUM) return rfmt(o.n);
  if (o.t == Val::BOOL) return o.b ? "TRUE" : "FALSE";
  if (o.t == Val::ARR) {
    std::string s = "c(";
    for (std::size_t i = 0; i < o.a.size(); ++i) s += (i ? ", " : "") + jdep(o.a[i]);
    return s + ")";
  }
  return "NULL";
}

// .ckarg(): one of the choices, matched partially; left out, null or every
// choice in order picks the first.
inline bool ckarg(Ctx &c, const char *nm, const std::vector<std::string> &ch, int &k) {
  k = 0;
  const Val &o = c[nm];
  if (o.t == Val::NUL) return true;
  if (o.t == Val::ARR && o.a.size() == ch.size()) {
    bool all = true;
    for (std::size_t i = 0; i < ch.size(); ++i) all = all && o.a[i].t == Val::STR && o.a[i].s == ch[i];
    if (all) return true;
  }
  const Val *s = &o;
  if (o.t == Val::ARR && o.a.size() == 1) s = &o.a[0];
  int hit = -1, part = 0;
  if (s->t == Val::STR && !s->s.empty()) {
    for (std::size_t i = 0; i < ch.size(); ++i) {
      if (ch[i] == s->s) hit = static_cast<int>(i);
    }
    if (hit < 0) {
      for (std::size_t i = 0; i < ch.size(); ++i) {
        if (ch[i].compare(0, s->s.size(), s->s) == 0) {
          ++part;
          hit = static_cast<int>(i);
        }
      }
      if (part != 1) hit = -1;
    }
  }
  if (hit < 0) {
    std::string sh;
    for (std::size_t i = 0; i < ch.size(); ++i) sh += (i ? ", \"" : "\"") + ch[i] + "\"";
    return ckstop(c, "`%s` must be one of %s, not %s", nm, sh.c_str(), jdep(o).c_str());
  }
  k = hit;
  return true;
}

inline Val nvval(const Nv &d) { return d.nr >= 0 ? vmat(d.v, d.nr, d.nc) : vauto(d.v); }

inline bool sarincr(Ctx &c, const Val &o, const char *nm, Nv &d) {
  if (o.t == Val::ARR && !o.a.empty() && o.a[0].t == Val::ARR) {
    return ckstop(c, "`%s` must be a vector of edges, not a matrix", nm);
  }
  CK(cknum(c, o, nm, d, false));
  bool up = d.n() >= 2;
  for (std::size_t i = 1; up && i < d.n(); ++i) up = d[i] > d[i - 1];
  if (!up) return ckstop(c, "`%s` must hold at least two strictly increasing values", nm);
  return true;
}

// R's format() of a sum of probabilities (7 significant digits).
inline std::string fmt7(double x) {
  char b[64];
  std::snprintf(b, sizeof b, "%.7g", x);
  return b;
}

inline bool sarpoc(Ctx &c, const Val &o, Nv &d) {
  CK(cknum(c, o, "poc", d, false, 0, 1));
  double s = 0;
  bool any = false;
  for (double x : d.v) {
    s += x;
    any = any || x > 0;
  }
  if (s > 1 + 1e-9) return ckstop(c, "`poc` must sum to at most 1, got %s", fmt7(s).c_str());
  if (!any) return ckstop(c, "`poc` must have at least one cell with a positive probability");
  return true;
}

// A datum map given as the object datumpmap() returns: one with poc, xedg
// and yedg.
inline bool ismap(const Val &o) { return o.t == Val::OBJ && col(o, "poc").t != Val::NUL && col(o, "xedg").t != Val::NUL && col(o, "yedg").t != Val::NUL; }

inline std::vector<double> mapedg(const Val &o, const char *k) {
  Dv d;
  readd(col(o, k), d);
  return d.v;
}

// outer(diff(yedg), diff(xedg)): a map's cell areas, column-major.
inline std::vector<double> mapar(const Val &m) {
  const std::vector<double> xe = mapedg(m, "xedg"), ye = mapedg(m, "yedg");
  std::vector<double> a;
  for (std::size_t j = 0; j + 1 < xe.size(); ++j) {
    for (std::size_t i = 0; i + 1 < ye.size(); ++i) a.push_back((ye[i + 1] - ye[i]) * (xe[j + 1] - xe[j]));
  }
  return a;
}

const std::vector<std::string> sarlrm = {"invcube", "definite", "mbeta"};

bool ltrlrng(Ctx &c, Val &out) {
  int md = 0;
  CK(ckarg(c, "model", sarlrm, md));
  Nv x, W, m, sd;
  CK(cknum(c, c["x"], "x", x, false));
  CK(cknum(c, c["W"], "W", W, true, -HUGE_VAL, HUGE_VAL, true));
  CK(cknumd(c, "m", 1, m, 0, 1, true));
  CK(cknumd(c, "navsd", 0, sd, 0));
  std::vector<double> p(x.n());
  CK(kchk(c, clab_sarlr(x.p(), static_cast<int>(x.n()), W[0], md, m[0], sd[0], p.data())));
  out = vauto(p);
  return true;
}

bool swpwdth(Ctx &c, Val &out) {
  Nv x, p;
  CK(sarincr(c, c["x"], "x", x));
  CK(cknum(c, c["p"], "p", p, false, 0, 1));
  if (p.n() != x.n()) {
    return ckstop(c, "`p` must have one value per `x` (%d), got %d", static_cast<int>(x.n()),
                  static_cast<int>(p.n()));
  }
  double w = 0;
  CK(kchk(c, clab_sarsw(x.p(), p.p(), static_cast<int>(x.n()), &w)));
  out = vnum(w);
  return true;
}

bool srchpod(Ctx &c, Val &out) {
  int md = 0;
  CK(ckarg(c, "model", {"random", "definite", "invcube", "mbeta"}, md));
  const bool hw = c.has("W"), he = c.has("effort"), ha = c.has("area");
  std::vector<double> cv;
  if (c.has("coverage")) {
    if (hw || he || ha) return ckstop(c, "give either `coverage` or `W`, `effort` and `area`, not both");
    Nv d;
    CK(cknum(c, c["coverage"], "coverage", d, false, 0));
    cv = d.v;
  } else {
    if (!(hw && he && ha)) return ckstop(c, "give `coverage`, or all three of `W`, `effort` and `area`");
    Nv W, e, a;
    CK(cknum(c, c["W"], "W", W, true, -HUGE_VAL, HUGE_VAL, true));
    CK(cknum(c, c["effort"], "effort", e, false, 0));
    if (ismap(c["area"])) a.v = mapar(c["area"]);
    else CK(cknum(c, c["area"], "area", a, false, -HUGE_VAL, HUGE_VAL, true));
    if (e.n() != a.n() && std::min(e.n(), a.n()) > 1) {
      return ckstop(c, "`effort` and `area` must have the same length, or one of them length 1");
    }
    const std::size_t n = std::max(e.n(), a.n());
    for (std::size_t i = 0; i < n; ++i) cv.push_back(W[0] * e[i % e.n()] / a[i % a.n()]);
  }
  Nv m;
  CK(cknumd(c, "m", 1, m, 0, 1, true));
  std::vector<double> p(cv.size());
  CK(kchk(c, clab_sarpod(cv.data(), static_cast<int>(cv.size()), md, m[0], p.data())));
  out = vauto(p);
  return true;
}

bool pswppod(Ctx &c, Val &out) {
  int md = 0;
  CK(ckarg(c, "model", sarlrm, md));
  Nv W, S, m, sd, np, nt;
  CK(cknum(c, c["W"], "W", W, true, -HUGE_VAL, HUGE_VAL, true));
  CK(cknum(c, c["S"], "S", S, true, -HUGE_VAL, 1e12, true));
  CK(cknumd(c, "m", 1, m, 0, 1, true));
  CK(cknumd(c, "navsd", 0, sd, 0));
  CK(cknumd(c, "npts", 2001, np, 3, 1e5, false, true));
  CK(cknumd(c, "ntrk", 60, nt, 1, 1e4, false, true));
  double p = 0;
  CK(kchk(c, clab_sarpsw(W[0], S[0], md, m[0], sd[0], static_cast<int>(np[0]), static_cast<int>(nt[0]), &p)));
  out = vnum(p);
  return true;
}

bool srchscs(Ctx &c, Val &out) {
  const Val &pm = c["poc"];
  const bool map = ismap(pm);
  Val pv = map ? col(pm, "poc") : pm;
  Val dv = c["pod"];
  if (dv.t == Val::OBJ && col(dv, "pod").t != Val::NUL) dv = col(dv, "pod");
  Nv poc, pod;
  CK(sarpoc(c, pv, poc));
  // a POD matrix laid out like the map is one search, cell by cell
  if (map && poc.nr >= 0 && dv.t == Val::ARR && !dv.a.empty() && dv.a[0].t == Val::ARR) {
    int dim[2] = {-1, -1};
    if (shape(dv, 2, dim, 0) && dim[0] == poc.nr && dim[1] == poc.nc) {
      std::vector<double> v(static_cast<std::size_t>(dim[0]) * static_cast<std::size_t>(dim[1]));
      bool bad = false;
      fill(dv, 2, dim, 0, 0, 1, v, bad, true);
      dv = vnums(v);
    }
  }
  CK(cknum(c, dv, "pod", pod, false, 0, 1));
  const int ns = pod.nr >= 0 ? pod.nr : 1, nd = pod.nr >= 0 ? pod.nc : static_cast<int>(pod.n());
  if (nd != static_cast<int>(poc.n())) {
    return ckstop(c, "`pod` must have one value (or column) per area of `poc` (%d), got %d",
                  static_cast<int>(poc.n()), nd);
  }
  const std::size_t nar = poc.n();
  std::vector<double> pa(nar), up(nar), cp(nar);
  double pos = 0;
  int cert = 0;
  CK(kchk(c, clab_sarscs(poc.p(), static_cast<int>(nar), pod.p(), ns, &pos, pa.data(), up.data(), cp.data(), &cert)));
  if (cert) {
    c.warn.push_back("srchscs: the search is certain to find the object (POS = 1), so the updated containment is "
                     "undefined and returned as NA.");
  }
  for (double &v : up) v = na(v);
  out = vobj();
  put(out, "pos", vnum(pos));
  put(out, "posarea", vauto(pa));
  put(out, "pocupd", vauto(up));
  put(out, "cumpod", vauto(cp));
  put(out, "poc", vauto(poc.v));
  put(out, "nsrch", vnum(ns));
  if (map) {
    Val m = pm;
    double s = 0;
    for (double v : up) s += v;
    for (std::size_t i = 0; i < m.k.size(); ++i) {
      if (m.k[i] == "poc") m.a[i] = vmat(up, poc.nr, poc.nc);
      if (m.k[i] == "outside") m.a[i] = vnum(std::isnan(s) ? NA : 1 - s);
    }
    put(out, "map", m);
  }
  return true;
}

bool datumpmap(Ctx &c, Val &out) {
  int dt = 0;
  static const std::vector<std::string> ch = {"point", "line", "area", "corridor", "general"};
  CK(ckarg(c, "datum", ch, dt));
  Nv xe, ye, cen, sg, end, nd;
  CK(sarincr(c, c["xedg"], "xedg", xe));
  CK(sarincr(c, c["yedg"], "yedg", ye));
  if (c.has("center")) {
    CK(cknum(c, c["center"], "center", cen, false));
  } else {
    cen.v = {0, 0};
  }
  if (cen.n() != 2) return ckstop(c, "`center` must be two numbers, x and y");
  CK(cknumd(c, "sigma", 1, sg, -HUGE_VAL, HUGE_VAL, true));
  const bool ne = dt >= 1 && dt <= 3;
  if (ne) {
    if (!c.has("end")) return ckstop(c, "a %s datum needs `end`, two numbers", ch[static_cast<std::size_t>(dt)].c_str());
    CK(cknum(c, c["end"], "end", end, false));
    if (end.n() != 2) return ckstop(c, "`end` must be two numbers, x and y");
  }
  CK(cknumd(c, "nodes", 64, nd, 2, 1000, false, true));
  const int nx = static_cast<int>(xe.n()) - 1, ny = static_cast<int>(ye.n()) - 1;
  if (dt == 3) {
    if (cen[1] != end[1] && cen[0] != end[0]) return ckstop(c, "a corridor datum must run along x or along y");
    if (cen[0] == end[0] && cen[1] == end[1]) return ckstop(c, "a corridor datum needs `center` and `end` to differ");
  }
  if (dt == 2 && (cen[0] == end[0] || cen[1] == end[1])) {
    return ckstop(c, "an area datum needs `center` and `end` to differ in x and in y");
  }
  Nv w;
  if (dt == 4) {
    if (!c.has("wts")) return ckstop(c, "a general datum needs `wts`, a matrix of cell weights");
    CK(cknum(c, c["wts"], "wts", w, false, 0));
    if (w.nr != ny || w.nc != nx) {
      return ckstop(c, "`wts` must be a matrix with one row per y cell (%d) and one column per x cell (%d)", ny, nx);
    }
    bool any = false;
    for (double v : w.v) any = any || v > 0;
    if (!any) return ckstop(c, "`wts` must have at least one positive weight");
  }
  const double e0[2] = {0, 0};
  std::vector<double> g(static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny));
  CK(kchk(c, clab_sarmap(xe.p(), static_cast<int>(xe.n()), ye.p(), static_cast<int>(ye.n()), dt, cen.p(), sg[0],
                         ne ? end.p() : e0, static_cast<int>(nd[0]), dt == 4 ? w.p() : nullptr, g.data())));
  double s = 0;
  for (double v : g) s += v;
  out = vobj();
  put(out, "poc", vmat(g, ny, nx));
  put(out, "outside", vnum(1 - s));
  put(out, "xedg", vauto(xe.v));
  put(out, "yedg", vauto(ye.v));
  put(out, "datum", vstr(ch[static_cast<std::size_t>(dt)]));
  put(out, "center", vauto(cen.v));
  put(out, "end", ne ? vauto(end.v) : c["end"]);
  put(out, "sigma", vnum(sg[0]));
  return true;
}

bool optcrcsrch(Ctx &c, Val &out) {
  Nv s, W, e;
  CK(cknum(c, c["sigma"], "sigma", s, true, -HUGE_VAL, HUGE_VAL, true));
  CK(cknum(c, c["W"], "W", W, true, -HUGE_VAL, HUGE_VAL, true));
  CK(cknum(c, c["effort"], "effort", e, true, 0));
  double r[3];
  CK(kchk(c, clab_sarcrc(s[0], W[0], e[0], r)));
  out = vobj();
  static const char *const nm[] = {"radius", "pos", "cntrdns"};
  for (int k = 0; k < 3; ++k) put(out, nm[k], vnum(r[k]));
  put(out, "sigma", vnum(s[0]));
  put(out, "W", vnum(W[0]));
  put(out, "total", vnum(e[0]));
  return true;
}

bool opteffort(Ctx &c, Val &out) {
  const Val &pm = c["poc"];
  const bool map = ismap(pm);
  Nv poc, ar;
  CK(sarpoc(c, map ? col(pm, "poc") : pm, poc));
  if (c.has("area")) {
    CK(cknum(c, c["area"], "area", ar, false, -HUGE_VAL, HUGE_VAL, true));
  } else if (map) {
    ar.v = mapar(pm);
  } else {
    return ckstop(c, "`area` must be given, one positive number per cell");
  }
  if (ar.n() != poc.n()) {
    return ckstop(c, "`area` must have one value per cell of `poc` (%d), got %d", static_cast<int>(poc.n()),
                  static_cast<int>(ar.n()));
  }
  Nv W, e, tol;
  CK(cknum(c, c["W"], "W", W, true, -HUGE_VAL, HUGE_VAL, true));
  CK(cknum(c, c["effort"], "effort", e, true, 0));
  CK(cknumd(c, "tol", 1e-14, tol, -HUGE_VAL, HUGE_VAL, true));
  const std::size_t n = poc.n();
  std::vector<double> z(n), cv(n), pd(n);
  double pos = 0, lam = 0;
  CK(kchk(c, clab_saropt(poc.p(), ar.p(), static_cast<int>(n), W[0], e[0], tol[0], z.data(), cv.data(), pd.data(),
                         &pos, &lam)));
  const int nr = poc.nr >= 0 ? poc.nr : ar.nr, nc = poc.nr >= 0 ? poc.nc : ar.nc;
  auto shp = [&](const std::vector<double> &v) { return nr >= 0 ? vmat(v, nr, nc) : vauto(v); };
  out = vobj();
  put(out, "effort", shp(z));
  put(out, "coverage", shp(cv));
  put(out, "pod", shp(pd));
  put(out, "pos", vnum(pos));
  put(out, "lambda", vnum(lam));
  put(out, "poc", nvval(poc));
  put(out, "W", vnum(W[0]));
  put(out, "total", vnum(e[0]));
  if (map) {
    put(out, "xedg", col(pm, "xedg"));
    put(out, "yedg", col(pm, "yedg"));
  }
  return true;
}

bool optsqr(Ctx &c, Val &out) {
  int md = 0;
  CK(ckarg(c, "model", {"invcube", "random", "definite", "mbeta"}, md));
  static const int code[] = {2, 0, 1, 3};
  Nv e, m, E;
  CK(cknum(c, c["effort"], "effort", e, false, -HUGE_VAL, HUGE_VAL, true));
  CK(cknumd(c, "m", 1, m, 0, 1, true));
  CK(cknumd(c, "E", 1, E, -HUGE_VAL, HUGE_VAL, true));
  const std::size_t k = e.n();
  std::vector<double> z(k), o(8 * k);
  for (std::size_t i = 0; i < k; ++i) z[i] = e[i] / (E[0] * E[0]);
  CK(kchk(c, clab_sarsq(z.data(), static_cast<int>(k), code[md], m[0], o.data())));
  std::vector<int> sn;
  std::vector<std::vector<double>> cl(7);
  std::vector<double> side;
  for (std::size_t i = 0; i < k; ++i) {
    sn.push_back(static_cast<int>(i) + 1);
    for (std::size_t j = 0; j < 7; ++j) cl[j].push_back(o[8 * i + j]);
    side.push_back(2 * o[8 * i + 1] * E[0]);
  }
  Frame f;
  f.ints("search", sn.data(), k);
  f.num("effort", e.v);
  f.num("releff", z);
  f.num("cumeff", cl[0]);
  f.num("factor", cl[1]);
  f.num("side", side);
  static const char *const nm[] = {"poc", "coverage", "pod", "pos", "cumpos"};
  for (int j = 0; j < 5; ++j) f.num(nm[j], cl[static_cast<std::size_t>(j) + 2]);
  out = f.o;
  return true;
}

bool srchfig(Ctx &c, Val &out) {
  const int n = clab_sarfig(-1, nullptr);
  std::vector<std::vector<std::string>> t(4);
  const char *row[4];
  for (int i = 0; i < n; ++i) {
    clab_sarfig(i, row);
    for (int k = 0; k < 4; ++k) t[static_cast<std::size_t>(k)].push_back(row[k]);
  }
  int pick = -1;
  if (c.has("fig")) {
    std::string fg;
    if (c.has("num")) {
      Nv a, b;
      CK(cknum(c, c["fig"], "fig", a, true, -HUGE_VAL, HUGE_VAL, false, true));
      CK(cknum(c, c["num"], "num", b, true, -HUGE_VAL, HUGE_VAL, false, true));
      fg = rfmt(a[0]) + "-" + rfmt(b[0]);
    } else {
      const Val &f = c["fig"];
      const Val *s = f.t == Val::ARR && f.a.size() == 1 ? &f.a[0] : &f;
      if (s->t != Val::STR || s->s.empty()) return ckstop(c, "`fig` must be a single non-empty string");
      fg = s->s;
    }
    for (int i = 0; i < n; ++i) {
      if (t[0][static_cast<std::size_t>(i)] == fg) pick = i;
    }
    if (pick < 0) {
      std::string all;
      for (int i = 0; i < n; ++i) all += (i ? ", " : "") + t[0][static_cast<std::size_t>(i)];
      return ckstop(c, "`fig` must be one of the book's figures, %s; got \"%s\"", all.c_str(), fg.c_str());
    }
  }
  static const char *const nm[] = {"fig", "caption", "how", "yours"};
  Frame f;
  for (int k = 0; k < 4; ++k) {
    std::vector<std::string> v = t[static_cast<std::size_t>(k)];
    if (pick >= 0) v = {v[static_cast<std::size_t>(pick)]};
    f.str(nm[k], v);
  }
  out = f.o;
  return true;
}
