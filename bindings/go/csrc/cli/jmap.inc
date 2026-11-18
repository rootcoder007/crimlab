// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The maps of R/map.R and the functions added beside them, for
// clab_jrun(), with the checks of R in their order and words. A map draws
// nothing here: it returns what the R function returns invisibly, the
// layout of src/map.cpp. The drawing arguments (main, col, pal and those
// of `...`) are not taken, so the colours are those of the default
// palettes. An area is an array of rows [x, y], a row of nulls between
// rings; an object of areas names them.

namespace mp {

// One area: rows of two numbers, both or neither missing, none infinite,
// at least three points.
inline bool ring(const Val &a, std::vector<double> &x, std::vector<double> &y) {
  if (a.t != Val::ARR) return false;
  int nf = 0;
  for (const Val &r : a.a) {
    if (r.t != Val::ARR || r.a.size() != 2 || !isnum(r)) return false;
    const double u = ax::nv(r.a[0]), v = ax::nv(r.a[1]);
    if (std::isnan(u) != std::isnan(v) || std::isinf(u) || std::isinf(v)) return false;
    x.push_back(u);
    y.push_back(v);
    nf += std::isnan(u) ? 0 : 1;
  }
  return nf >= 3;
}

// .mpoly(): the areas, all their coordinates (a NaN between areas) and
// their labels, the names or the positions from 1.
struct Poly {
  std::vector<double> ax, ay;
  Val lab;
  std::size_t n = 0;
};
inline bool poly(Ctx &c, Poly &p) {
  const Val &o = c["polys"];
  if ((o.t != Val::ARR && o.t != Val::OBJ) || o.a.empty()) {
    return c.fail("`polys` must be a non-empty list of areas.");
  }
  p.lab = varr();
  for (std::size_t i = 0; i < o.a.size(); ++i) {
    std::vector<double> x, y;
    if (!ring(o.a[i], x, y)) {
      return c.fail("`polys[[%d]]` must be a two-column numeric matrix of at least three points (rows of NA between "
                    "rings).",
                    sat<int>(i + 1));
    }
    p.ax.insert(p.ax.end(), x.begin(), x.end());
    p.ay.insert(p.ay.end(), y.begin(), y.end());
    p.ax.push_back(NaN);
    p.ay.push_back(NaN);
    p.lab.a.push_back(o.t == Val::OBJ ? vstr(o.k[i]) : vnum(static_cast<double>(i + 1)));
  }
  p.n = o.a.size();
  return true;
}

// .mbrk(): the breaks and the class of each value (0 for none).
inline bool brk(Ctx &c, const std::vector<double> &v, std::vector<double> &br, std::vector<int> &cls) {
  std::string meth = "quantile";
  if (c.has("meth")) {
    const Val &m = c["meth"];
    const Val &e = m.t == Val::ARR && m.a.size() == 1 ? m.a[0] : m;
    meth = e.t == Val::STR ? e.s : std::string();
  }
  if (meth != "quantile" && meth != "equal") return c.fail("`meth` must be \"quantile\" or \"equal\".");
  double lo = HUGE_VAL, hi = -HUGE_VAL;
  for (double x : v) {
    if (!std::isfinite(x)) continue;
    lo = std::min(lo, x);
    hi = std::max(hi, x);
  }
  if (!(lo <= hi)) return c.fail("`v` has no finite value.");
  const bool given = c.has("brk");
  const std::vector<double> b = given ? plt::vals(c["brk"]) : std::vector<double>{5.0};
  const bool one = b.size() == 1;
  bool ok = !given || isnum(c["brk"]);
  if (one) {
    ok = ok && std::isfinite(b[0]) && b[0] >= 1 && b[0] <= 9 && b[0] == std::round(b[0]);
  } else {
    ok = ok && b.size() >= 2;
    for (std::size_t j = 0; ok && j < b.size(); ++j) ok = std::isfinite(b[j]) && (j == 0 || b[j] > b[j - 1]);
    ok = ok && b.front() <= lo && b.back() >= hi;
  }
  if (!ok) return c.fail("`brk` must be a whole number of classes from 1 to 9, or increasing breaks covering `v`.");
  const int k = one ? sat<int>(b[0]) : 0, mt = one ? (meth == "equal" ? 0 : 1) : 2;
  int nb = one ? k + 1 : sat<int>(b.size());
  br.assign(sat<std::size_t>(nb), 0.0);
  if (!one) br = b;
  cls.assign(v.size(), 0);
  CK(kchk(c, clab_mapbrk(v.data(), sat<int>(v.size()), k, mt, br.data(), &nb, cls.data())));
  br.resize(sat<std::size_t>(nb));
  return true;
}

// blpal(name, n) of R/plot.R.
inline ax::Cols pal(const char *nm, int n) {
  const plt::Pal *p = &plt::pals[0];
  for (const plt::Pal &q : plt::pals) {
    if (std::strcmp(q.name, nm) == 0) p = &q;
  }
  ax::Cols r;
  if (std::strcmp(p->name, "okabe") == 0 || std::strcmp(p->name, "grey") == 0) {
    for (int i = 0; i < n; ++i) r.push_back(p->col[sat<std::size_t>(i) % p->col.size()]);
    return r;
  }
  if (n == 1) return {p->col[0]};
  return plt::ramp(p->col, n);
}

inline bool bbox(Ctx &c, const std::vector<double> &x, const std::vector<double> &y, std::vector<double> &bb) {
  bb.assign(4, 0.0);
  return kchk(c, clab_mapbb(x.data(), y.data(), sat<int>(x.size()), bb.data()));
}

// The columns cls (NA for none) and col (null for none) of a frame.
inline void cls(Frame &f, const std::vector<int> &k, const ax::Cols &pl) {
  std::vector<double> cv;
  Val col = varr();
  for (int v : k) {
    cv.push_back(v > 0 ? static_cast<double>(v) : NA);
    col.a.push_back(v > 0 ? vstr(pl[sat<std::size_t>(v - 1)]) : Val());
  }
  f.num("cls", cv);
  f.add("col", col);
}

}  // namespace mp

// chomap()
bool chomap(Ctx &c, Val &out) {
  mp::Poly p;
  CK(mp::poly(c, p));
  const Val &vo = c["v"];
  const std::vector<double> v = plt::vals(vo);
  if (!isnum(vo) || v.size() != p.n) return c.fail("`v` must be numbers, one per area.");
  std::vector<double> br, bb;
  std::vector<int> k;
  CK(mp::brk(c, v, br, k));
  CK(mp::bbox(c, p.ax, p.ay, bb));
  Frame f;
  f.add("area", p.lab);
  f.num("v", v);
  mp::cls(f, k, mp::pal("sequential", sat<int>(br.size()) - 1));
  out = vobj();
  put(out, "areas", f.o);
  put(out, "breaks", vnums(br));
  put(out, "bbox", vnums(bb));
  return true;
}

// smmap(): v an array of rows (an area each, periods 1, 2, ...) or an
// object of columns (a period each, named).
bool smmap(Ctx &c, Val &out) {
  mp::Poly p;
  CK(mp::poly(c, p));
  const Val &vo = c["v"];
  std::vector<double> v;
  ax::Cols prd;
  bool ok = false;
  if (vo.t == Val::OBJ && !vo.a.empty()) {
    ok = true;
    for (std::size_t j = 0; ok && j < vo.a.size(); ++j) {
      const Val &col = vo.a[j];
      ok = col.t == Val::ARR && col.a.size() == p.n && ax::isnumc(col);
      if (!ok) break;
      const std::vector<double> cv = plt::vals(col);
      v.insert(v.end(), cv.begin(), cv.end());
      prd.push_back(vo.k[j]);
    }
  } else if (vo.t == Val::ARR && vo.a.size() == p.n) {
    const std::size_t np = vo.a[0].t == Val::ARR ? vo.a[0].a.size() : 0;
    ok = np > 0;
    for (const Val &r : vo.a) ok = ok && r.t == Val::ARR && r.a.size() == np && isnum(r);
    if (ok) {
      v.assign(np * p.n, 0.0);
      for (std::size_t i = 0; i < p.n; ++i) {
        for (std::size_t j = 0; j < np; ++j) v[j * p.n + i] = ax::nv(vo.a[i].a[j]);
      }
      for (std::size_t j = 0; j < np; ++j) prd.push_back(std::to_string(j + 1));
    }
  }
  if (!ok) return c.fail("`v` must be a numeric matrix, one row per area and one column per period.");
  std::vector<double> br, bb;
  std::vector<int> k;
  CK(mp::brk(c, v, br, k));
  CK(mp::bbox(c, p.ax, p.ay, bb));
  Frame f;
  Val lab = varr();
  ax::Cols pv;
  for (const std::string &s : prd) {
    lab.a.insert(lab.a.end(), p.lab.a.begin(), p.lab.a.end());
    for (std::size_t i = 0; i < p.n; ++i) pv.push_back(s);
  }
  f.add("area", lab);
  f.str("prd", pv);
  f.num("v", v);
  mp::cls(f, k, mp::pal("sequential", sat<int>(br.size()) - 1));
  out = vobj();
  put(out, "areas", f.o);
  put(out, "breaks", vnums(br));
  put(out, "bbox", vnums(bb));
  return true;
}

namespace mp {

// .mxy(): x and y, or x alone an array of rows of two or an object of
// two numeric columns.
inline bool xy(Ctx &c, std::vector<double> &x, std::vector<double> &y) {
  const Val &xo = c["x"], &yo = c["y"];
  if (yo.t == Val::NUL) {
    bool ok = false;
    if (xo.t == Val::ARR && !xo.a.empty()) {
      ok = true;
      for (const Val &r : xo.a) {
        ok = ok && r.t == Val::ARR && r.a.size() == 2 && isnum(r);
        if (!ok) break;
        x.push_back(ax::nv(r.a[0]));
        y.push_back(ax::nv(r.a[1]));
      }
    } else if (xo.t == Val::OBJ && xo.a.size() == 2) {
      ok = xo.a[0].t == Val::ARR && xo.a[1].t == Val::ARR && xo.a[0].a.size() == xo.a[1].a.size() &&
           ax::isnumc(xo.a[0]) && ax::isnumc(xo.a[1]);
      if (ok) {
        x = plt::vals(xo.a[0]);
        y = plt::vals(xo.a[1]);
      }
    }
    if (!ok) {
      return c.fail("`x` must be a two-column matrix or data frame of coordinates, or numbers with `y`.");
    }
  } else {
    x = plt::vals(xo);
    y = plt::vals(yo);
    if (!isnum(xo) || !isnum(yo) || x.size() != y.size()) {
      return c.fail("`x` and `y` must be numbers of the same length.");
    }
  }
  bool fin = !x.empty();
  for (std::size_t i = 0; fin && i < x.size(); ++i) fin = std::isfinite(x[i]) && std::isfinite(y[i]);
  if (!fin) return c.fail("the coordinates must be finite, at least one point.");
  return true;
}

}  // namespace mp

// ptmap(): cl an array, or a dbscn result (an object with cl).
bool ptmap(Ctx &c, Val &out) {
  std::vector<double> x, y;
  CK(mp::xy(c, x, y));
  const std::size_t n = x.size();
  const bool given = c.has("cl");
  std::vector<int> cl(n, 0);
  if (given) {
    const Val *co = &c["cl"];
    if (co->t == Val::OBJ) {
      const Val *e = ax::fc(*co, "cl");
      if (e != nullptr) co = e;
    }
    const std::vector<double> v = plt::vals(*co);
    bool ok = co->t != Val::OBJ && isnum(*co) && v.size() == n;
    for (std::size_t i = 0; ok && i < n; ++i) {
      ok = std::isnan(v[i]) || (v[i] >= 1 && v[i] == std::round(v[i]) && v[i] <= INT_MAX);
      if (ok && !std::isnan(v[i])) cl[i] = sat<int>(v[i]);
    }
    if (!ok) return c.fail("`cl` must be a dbscn result or a cluster per point: whole numbers from 1, NA for noise.");
  }
  int k = 0;
  for (int v : cl) k = std::max(k, v);
  const ax::Cols pl = mp::pal("okabe", std::max(k, 1));
  static const int pch[8] = {16, 17, 15, 18, 1, 2, 0, 5};
  std::vector<double> cv(n, NA), ph(n, 16);
  Val col = varr();
  for (std::size_t i = 0; i < n; ++i) {
    if (!given) {
      col.a.push_back(vstr(pl[0]));
    } else if (cl[i] == 0) {
      col.a.push_back(vstr("#969696"));
      ph[i] = 4;
    } else {
      cv[i] = cl[i];
      col.a.push_back(vstr(pl[sat<std::size_t>(cl[i] - 1)]));
      ph[i] = pch[(cl[i] - 1) % 8];
    }
  }
  std::vector<double> bb;
  CK(mp::bbox(c, x, y, bb));
  const std::size_t kk = sat<std::size_t>(k);
  std::vector<double> cx(kk), cy(kk), id(kk);
  std::vector<int> cn(kk);
  CK(kchk(c, clab_ptcen(x.data(), y.data(), cl.data(), sat<int>(n), k, cx.data(), cy.data(), cn.data())));
  for (std::size_t j = 0; j < kk; ++j) id[j] = static_cast<double>(j + 1);
  Frame f, g;
  f.num("x", x);
  f.num("y", y);
  f.num("cl", cv);
  f.add("col", col);
  f.num("pch", ph);
  g.num("cl", id);
  g.num("x", cx);
  g.num("y", cy);
  g.ints("n", cn.data(), cn.size());
  out = vobj();
  put(out, "pts", f.o);
  put(out, "cent", g.o);
  put(out, "bbox", vnums(bb));
  return true;
}

// prpmap()
bool prpmap(Ctx &c, Val &out) {
  std::vector<double> x, y;
  CK(mp::xy(c, x, y));
  const std::size_t n = x.size();
  const Val &vo = c["v"];
  const std::vector<double> v = plt::vals(vo);
  bool ok = isnum(vo) && v.size() == n, any = false;
  for (std::size_t i = 0; ok && i < n; ++i) {
    if (std::isnan(v[i])) continue;
    any = true;
    ok = std::isfinite(v[i]) && v[i] >= 0;
  }
  if (!ok || !any) return c.fail("`v` must be numbers from 0, one per point, not all missing.");
  double rmax = NaN;
  if (c.has("rmax")) {
    const std::vector<double> r = plt::vals(c["rmax"]);
    if (!isnum(c["rmax"]) || r.size() != 1 || !std::isfinite(r[0]) || r[0] < std::numeric_limits<double>::min()) {
      return c.fail("`rmax` must be one positive number.");
    }
    rmax = r[0];
  }
  std::vector<double> bb(4), r(n);
  std::vector<int> od(n);
  int no = 0;
  double rm = 0.0;
  CK(kchk(c, clab_mapsym(x.data(), y.data(), v.data(), sat<int>(n), rmax, bb.data(), r.data(), od.data(),
                         &no, &rm)));
  std::vector<double> dr(n, NA);
  for (int q = 0; q < no; ++q) dr[sat<std::size_t>(od[sat<std::size_t>(q)])] = q + 1;
  for (std::size_t i = 0; i < n; ++i) {
    if (std::isnan(v[i])) r[i] = NA;
  }
  Frame f;
  f.num("x", x);
  f.num("y", y);
  f.num("v", v);
  f.num("r", r);
  f.num("ord", dr);
  out = vobj();
  put(out, "pts", f.o);
  put(out, "rmax", vnum(rm));
  put(out, "bbox", vnums(bb));
  return true;
}

// eofst(): z an array of rows (periods) or an object of columns (places,
// named).
bool eofst(Ctx &c, Val &out) {
  const Val &zo = c["z"];
  std::size_t m = 0, p = 0;
  std::vector<double> a;
  ax::Cols cn;
  bool ok = false;
  if (zo.t == Val::OBJ && !zo.a.empty()) {
    ok = true;
    m = zo.a[0].t == Val::ARR ? zo.a[0].a.size() : 0;
    p = zo.a.size();
    for (std::size_t j = 0; ok && j < p; ++j) {
      ok = zo.a[j].t == Val::ARR && zo.a[j].a.size() == m && ax::isnumc(zo.a[j]);
      if (!ok) break;
      const std::vector<double> v = plt::vals(zo.a[j]);
      a.insert(a.end(), v.begin(), v.end());
      cn.push_back(zo.k[j]);
    }
  } else if (zo.t == Val::ARR && !zo.a.empty()) {
    m = zo.a.size();
    p = zo.a[0].t == Val::ARR ? zo.a[0].a.size() : 0;
    ok = true;
    for (const Val &r : zo.a) ok = ok && r.t == Val::ARR && r.a.size() == p && isnum(r);
    if (ok) {
      a.assign(m * p, 0.0);
      for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < p; ++j) a[j * m + i] = ax::nv(zo.a[i].a[j]);
      }
    }
  }
  if (!ok || m < 2 || p < 1) return c.fail("`z` must be a numeric matrix of at least two rows: periods by places.");
  for (double v : a) {
    if (!std::isfinite(v)) return c.fail("`z` must be finite.");
  }
  const std::size_t s = std::min(m, p);
  int k = 0;
  if (c.has("k")) {
    const std::vector<double> kv = plt::vals(c["k"]);
    if (!isnum(c["k"]) || kv.size() != 1 || !std::isfinite(kv[0]) || kv[0] < 1 || kv[0] != std::round(kv[0])) {
      return c.fail("`k` must be NULL or one whole number from 1.");
    }
    k = sat<int>(std::min(kv[0], static_cast<double>(s)));
  }
  std::vector<double> ctr(p), ev(p * s), pc(m * s), d2(s);
  double tot = 0.0;
  int kk = 0;
  const int rc = clab_eofst(a.data(), sat<int>(m), sat<int>(p), k, ctr.data(), ev.data(), pc.data(),
                            d2.data(), &tot, &kk);
  if (rc == CLAB_EDOM) return c.fail("`z` has no variation: every column is constant.");
  CK(kchk(c, rc));
  std::vector<double> ex, sd;
  for (int j = 0; j < kk; ++j) {
    const double d = d2[sat<std::size_t>(j)];
    ex.push_back(d / tot);
    sd.push_back(std::sqrt(d) / std::sqrt(static_cast<double>(m) - 1));
  }
  Val ce = vobj();
  if (cn.empty()) {
    ce = unbox(vnums(ctr));
  } else {
    for (std::size_t j = 0; j < p; ++j) put(ce, cn[j].c_str(), vnum(ctr[j]));
  }
  out = vobj();
  put(out, "eofs", vmat(ev.data(), sat<int>(p), kk));
  put(out, "pcs", vmat(pc.data(), sat<int>(m), kk));
  put(out, "explained", unbox(vnums(ex)));
  put(out, "sdev", unbox(vnums(sd)));
  put(out, "center", ce);
  return true;
}

// tpsurf()
bool tpsurf(Ctx &c, Val &out) {
  const Val &xo = c["x"], &yo = c["y"], &zo = c["z"];
  const std::vector<double> x = plt::vals(xo), y = plt::vals(yo), z = plt::vals(zo);
  if (!isnum(xo) || !isnum(yo) || !isnum(zo) || y.size() != x.size() || z.size() != x.size()) {
    return c.fail("`x`, `y` and `z` must be numbers of the same length.");
  }
  bool fin = x.size() >= 3;
  for (std::size_t i = 0; fin && i < x.size(); ++i) {
    fin = std::isfinite(x[i]) && std::isfinite(y[i]) && std::isfinite(z[i]);
  }
  if (!fin) return c.fail("`x`, `y` and `z` must be finite, at least three points.");
  double lam = 0.0;
  if (c.has("lam")) {
    const std::vector<double> l = plt::vals(c["lam"]);
    if (!isnum(c["lam"]) || l.size() != 1 || !std::isfinite(l[0]) || l[0] < 0) {
      return c.fail("`lam` must be one number from 0.");
    }
    lam = l[0];
  }
  double g[2] = {50, 50};
  const char *gn[2] = {"nx", "ny"};
  for (int q = 0; q < 2; ++q) {
    if (!c.has(gn[q])) {
      g[q] = q == 1 ? g[0] : 50;
      continue;
    }
    const std::vector<double> v = plt::vals(c[gn[q]]);
    if (!isnum(c[gn[q]]) || v.size() != 1 || !std::isfinite(v[0]) || v[0] < 2 || v[0] > 1000 ||
        v[0] != std::round(v[0])) {
      return c.fail("`nx` and `ny` must be whole numbers from 2 to 1000.");
    }
    g[q] = v[0];
  }
  if (x.size() > 4000) return c.fail("at most 4000 points.");
  const int n = sat<int>(x.size()), nx = sat<int>(g[0]), ny = sat<int>(g[1]);
  const auto rx = std::minmax_element(x.begin(), x.end()), ry = std::minmax_element(y.begin(), y.end());
  const std::vector<double> gx = ax::seqn(*rx.first, *rx.second, nx), gy = ax::seqn(*ry.first, *ry.second, ny);
  std::vector<double> a(x.size()), b(3), fit(x.size()), gz(sat<std::size_t>(nx) * sat<std::size_t>(ny));
  const int rc = clab_tpsurf(x.data(), y.data(), z.data(), n, lam, gx.data(), nx, gy.data(), ny, a.data(), b.data(),
                             fit.data(), gz.data());
  if (rc == CLAB_EDOM) return c.fail("the system is singular: the points lie on one line, or repeat with `lam` 0.");
  CK(kchk(c, rc));
  std::vector<double> res(x.size());
  for (std::size_t i = 0; i < x.size(); ++i) res[i] = z[i] - fit[i];
  Val be = vobj();
  put(be, "const", vnum(b[0]));
  put(be, "x", vnum(b[1]));
  put(be, "y", vnum(b[2]));
  out = vobj();
  put(out, "x", vnums(gx));
  put(out, "y", vnums(gy));
  put(out, "z", vmat(gz.data(), nx, ny));
  put(out, "a", vnums(a));
  put(out, "beta", be);
  put(out, "fitted", vnums(fit));
  put(out, "resid", vnums(res));
  put(out, "lam", vnum(lam));
  return true;
}

// dagcol()
bool dagcol(Ctx &c, Val &out) {
  const Val &fo = c["from"], &to = c["to"];
  const std::vector<Val> f = plt::elems(fo), t = plt::elems(to);
  bool ok = fo.t != Val::OBJ && to.t != Val::OBJ && !f.empty() && f.size() == t.size();
  for (std::size_t i = 0; ok && i < f.size(); ++i) {
    ok = f[i].t != Val::ARR && f[i].t != Val::OBJ && t[i].t != Val::ARR && t[i].t != Val::OBJ;
  }
  if (!ok) return c.fail("`from` and `to` must be labels of the same length, one edge each, at least one edge.");
  for (std::size_t i = 0; i < f.size(); ++i) {
    if (ax::vna(f[i]) || ax::vna(t[i])) return c.fail("`from` and `to` must have no missing values.");
  }
  ax::Cols nodes;
  std::vector<int> fc, tc;
  const auto code = [&nodes](const std::string &s) {
    const auto it = std::find(nodes.begin(), nodes.end(), s);
    if (it != nodes.end()) return sat<int>(it - nodes.begin());
    nodes.push_back(s);
    return sat<int>(nodes.size()) - 1;
  };
  for (std::size_t i = 0; i < f.size(); ++i) {
    fc.push_back(code(ax::vlab(f[i])));
    tc.push_back(code(ax::vlab(t[i])));
  }
  const int ne = sat<int>(f.size()), nn = sat<int>(nodes.size());
  int nc = 0, nd = 0;
  int rc = clab_dagcol(fc.data(), tc.data(), ne, nn, 0, nullptr, nullptr, nullptr, nullptr, &nc, &nd);
  if (rc == CLAB_EDOM) return c.fail("the edges have a cycle; a DAG has none.");
  CK(kchk(c, rc));
  const std::size_t k = sat<std::size_t>(nc);
  std::vector<int> a(k), cc(k), b(k), sh(k);
  CK(kchk(c, clab_dagcol(fc.data(), tc.data(), ne, nn, nc, a.data(), cc.data(), b.data(), sh.data(), &nc, &nd)));
  ax::Cols la, lc, lb;
  for (std::size_t i = 0; i < k; ++i) {
    la.push_back(nodes[sat<std::size_t>(a[i])]);
    lc.push_back(nodes[sat<std::size_t>(cc[i])]);
    lb.push_back(nodes[sat<std::size_t>(b[i])]);
  }
  Val q;
  if (c.has("triple")) {
    const Val &tr = c["triple"];
    const std::vector<Val> e = plt::elems(tr);
    bool tok = tr.t != Val::OBJ && e.size() == 3;
    std::string s[3];
    for (std::size_t i = 0; tok && i < 3; ++i) {
      tok = e[i].t != Val::ARR && e[i].t != Val::OBJ && !ax::vna(e[i]);
      if (tok) s[i] = ax::vlab(e[i]);
      tok = tok && std::find(nodes.begin(), nodes.end(), s[i]) != nodes.end();
    }
    if (!tok) return c.fail("`triple` must be three labels of nodes: a, c, b.");
    bool hit = false, shd = false;
    for (std::size_t i = 0; i < k; ++i) {
      const bool h = lc[i] == s[1] && ((la[i] == s[0] && lb[i] == s[2]) || (la[i] == s[2] && lb[i] == s[0]));
      hit = hit || h;
      shd = shd || (h && sh[i] != 0);
    }
    q = vobj();
    put(q, "collider", vbool(hit));
    put(q, "shielded", vbool(shd));
  }
  Frame fr;
  fr.str("a", la);
  fr.str("c", lc);
  fr.str("b", lb);
  fr.lgl("shielded", sh.data(), sh.size());
  out = vobj();
  put(out, "col", fr.o);
  put(out, "nodes", vstrs(nodes));
  put(out, "nedge", vnum(nd));
  put(out, "query", q);
  return true;
}
