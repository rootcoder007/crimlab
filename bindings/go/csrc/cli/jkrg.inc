// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Geostatistics (src/bkrg.cpp): the semivariogram and variogram models,
// point and block kriging, cross-validation, trend surfaces, the spatial
// linear and mixed models, and the spatial data object of R/krg.R.

// zmat() of src/bkh.h: a design argument, column-major, p columns.
bool zmat(Ctx &c, const char *nm, int n, bool one, int &p, std::vector<double> &v) {
  const Val o = asmat(c[nm]);
  if (o.t == Val::NUL) {
    p = one ? 1 : 0;
    v.assign(one ? zz(n) : 0, 1.0);
    return true;
  }
  if (!isnum(o) && !(o.t == Val::ARR && !o.a.empty() && o.a[0].t == Val::ARR)) {
    return c.fail("`%s` must be a numeric matrix or vector.", nm);
  }
  if (o.t == Val::ARR && !o.a.empty() && o.a[0].t == Val::ARR) {
    int dim[2] = {-1, -1};
    if (!shape(o, 2, dim, 0)) return c.fail("`%s` must be a numeric matrix or vector.", nm);
    if (dim[0] != n) return c.fail("`%s` must have %d rows.", nm, n);
    int nr = 0;
    if (!nmat(c, o, nm, nr, p, v, true)) return c.fail("`%s` must be a numeric matrix or vector.", nm);
  } else {
    v = dbl(o);
    if (static_cast<int>(v.size()) != n) return c.fail("`%s` must have %d rows.", nm, n);
    p = 1;
  }
  for (double e : v) {
    if (!std::isfinite(e)) return c.fail("`%s` must hold finite numbers.", nm);
  }
  return true;
}

const char *const vmodels[] = {"spherical", "exponential", "gaussian"};

bool xyv(Ctx &c, Pts &p, Dv &v) {
  CK(xy(c, c["x"], c["y"], PLANAR, p));
  return num(c, c["v"], "v", v, p.n);
}

bool vgemp(Ctx &c, Val &out) {
  Pts p;
  CK(xy(c, c["x"], c["y"], PLANAR, p));
  if (p.n < 2) return c.fail("`x` must hold at least 2 points.");
  Dv v;
  CK(num(c, c["v"], "v", v, p.n));
  int nb = 0;
  CK(wholed(c, "nbins", 15, 1, nb));
  if (nb > kutil::maxcell) return c.fail("`nbins` must be at most %.0f.", kutil::maxcell);
  double cut = NaN;
  if (c.has("cutoff")) CK(one(c, c["cutoff"], "cutoff", 0.0, Inf, cut));
  std::vector<double> d(zz(nb)), g(zz(nb)), np(zz(nb));
  int nk = 0;
  CK(kc(c, clab_vgemp(p.x.data(), p.y.data(), p.n, v.p(), nb, cut, d.data(), g.data(), np.data(), &nk),
        "the points all coincide."));
  std::vector<int> ni(zz(nk));
  for (int k = 0; k < nk; ++k) ni[zz(k)] = static_cast<int>(np[zz(k)]);
  Frame f;
  f.num("dist", d.data(), zz(nk));
  f.num("gamma", g.data(), zz(nk));
  f.ints("np", ni.data(), zz(nk));
  out = f.o;
  return true;
}

bool vgfit(Ctx &c, Val &out) {
  Pts p;
  CK(xy(c, c["x"], c["y"], PLANAR, p));
  if (p.n < 3) return c.fail("`x` must hold at least 3 points.");
  Dv v;
  CK(num(c, c["v"], "v", v, p.n));
  int md = 1;
  CK(pickd(c, "model", vmodels, 3, 1, md));
  double o[6];
  CK(kc(c, clab_vgfit(p.x.data(), p.y.data(), p.n, v.p(), md, o),
        "`v` must vary, and the points must not all coincide."));
  out = vobj();
  put(out, "model", vstr(vmodels[md]));
  static const char *const nm[] = {"nugget", "psill", "range", "loglik"};
  for (int k = 0; k < 4; ++k) put(out, nm[k], vnum(o[k]));
  put(out, "conv", vbool(o[4] != 0));
  put(out, "nfn", vnum(o[5]));
  return true;
}

bool okrge(Ctx &c, Val &out) {
  Pts p, a;
  Dv v;
  CK(xyv(c, p, v));
  CK(xy(c, c["at"], nul, PLANAR, a));
  int md = 1;
  double par[3];
  const Val &g = c["vgm"];
  if (g.t == Val::NUL) {
    if (p.n < 2) return c.fail("`x` must hold at least 2 points to fit a variogram.");
    CK(kc(c, clab_vgwls(p.x.data(), p.y.data(), p.n, v.p(), 1, par),
          "no pair of points is close enough for a variogram."));
  } else {
    if (g.t != Val::OBJ) return c.fail("`vgm` must be a list with `model`, `nugget`, `psill` and `range`.");
    CK(pick(c, col(g, "model"), "vgm$model", vmodels, 3, md));
    Dv a1, a2, a3;
    CK(num(c, col(g, "nugget"), "vgm$nugget", a1, 1, 0.0));
    CK(num(c, col(g, "psill"), "vgm$psill", a2, 1, 0.0));
    CK(num(c, col(g, "range"), "vgm$range", a3, 1, 0.0, true));
    par[0] = a1[0];
    par[1] = a2[0];
    par[2] = a3[0];
  }
  std::vector<double> pr(zz(a.n)), vr(zz(a.n));
  CK(kc(c, clab_okrge(p.x.data(), p.y.data(), p.n, v.p(), a.x.data(), a.y.data(), a.n, md, par[0], par[1],
                      par[2], pr.data(), vr.data()),
        "the kriging system is singular: two points coincide without a nugget, or the variogram is flat."));
  Frame f;
  f.num("pred", pr);
  f.num("var", vr);
  out = f.o;
  return true;
}

// the covariance parameters of the block methods
bool covpar(Ctx &c, double *cp) {
  Val d;
  Dv a, b;
  CK(num(c, arg(c, "nugget", d, 0), "nugget", a, 1, 0.0));
  Val e;
  CK(num(c, arg(c, "sill", e, 1), "sill", b, 1));
  cp[0] = a[0];
  cp[1] = b[0];
  CK(oned(c, "range", 1, 0.0, Inf, cp[2]));
  if (cp[1] < cp[0]) return c.fail("`sill` must be at least `nugget`.");
  return true;
}

bool bkrig(Ctx &c, Val &out) {
  Pts p;
  Dv v;
  CK(xyv(c, p, v));
  int q = 0;
  CK(wholed(c, "nq", 25, 1, q));
  double cp[3];
  CK(covpar(c, cp));
  const Val &bl = c["blocks"];
  if (bl.t != Val::ARR || bl.a.empty()) return c.fail("`blocks` must be a list of blocks.");
  const int nb = static_cast<int>(bl.a.size());
  const int k = static_cast<int>(std::nearbyint(std::sqrt(static_cast<double>(q))));
  std::vector<int> off(zz(nb) + 1, 0);
  std::vector<std::vector<double>> m(zz(nb));
  std::vector<int> nr(zz(nb));
  double tot = 0;  // the points of all blocks, counted before any is made
  for (int b = 0; b < nb; ++b) {
    const Val &e = bl.a[zz(b)];
    int dim[2] = {-1, -1};
    if (!(e.t == Val::ARR && !e.a.empty() && e.a[0].t == Val::ARR && shape(e, 2, dim, 0) && dim[1] == 2 &&
          isnum(e.a[0]))) {
      return c.fail("`blocks[[%d]]` must be a numeric matrix of two columns.", b + 1);
    }
    nr[zz(b)] = dim[0];
    m[zz(b)].assign(zz(dim[0]) * 2, 0.0);
    bool bad = false;
    fill(e, 2, dim, 0, 0, 1, m[zz(b)], bad, true);
    tot += dim[0] == 2 ? static_cast<double>(k) * k : dim[0];
    if (tot > kutil::maxcell) {
      return c.fail("the points of `blocks`, `nq` for each box, must be at most %.0f.", kutil::maxcell);
    }
    off[zz(b) + 1] = off[zz(b)] + (dim[0] == 2 ? k * k : dim[0]);
  }
  std::vector<double> bx(zz(off[zz(nb)])), by(bx);
  for (int b = 0; b < nb; ++b) {
    const std::vector<double> &mm = m[zz(b)];
    const int r = nr[zz(b)];
    for (double e : mm) {
      if (!std::isfinite(e)) return c.fail("`blocks[[%d]]` must hold finite numbers.", b + 1);
    }
    const std::size_t o = zz(off[zz(b)]);
    if (r == 2) {
      const double lo[2] = {mm[0], mm[2]}, hi[2] = {mm[1], mm[3]};
      clab_bkpts(lo, hi, q, bx.data() + o, by.data() + o);
    } else {
      for (int i = 0; i < r; ++i) {
        bx[o + zz(i)] = mm[zz(i)];
        by[o + zz(i)] = mm[zz(r + i)];
      }
    }
  }
  std::vector<double> es(zz(nb)), se(zz(nb));
  CK(kc(c, clab_bkrig(p.x.data(), p.y.data(), p.n, v.p(), bx.data(), by.data(), off.data(), nb, cp[0], cp[1],
                      cp[2], es.data(), se.data()),
        "the kriging system is singular: two points coincide without a nugget."));
  Frame f;
  f.num("est", es);
  f.num("se", se);
  out = f.o;
  return true;
}

bool krigcv(Ctx &c, Val &out) {
  Pts p;
  CK(xy(c, c["x"], c["y"], PLANAR, p));
  if (p.n < 3) return c.fail("`x` must hold at least 3 points.");
  Dv v;
  CK(num(c, c["v"], "v", v, p.n));
  double cp[3];
  CK(covpar(c, cp));
  std::vector<double> r(zz(p.n));
  CK(kc(c, clab_krigcv(p.x.data(), p.y.data(), p.n, v.p(), cp[0], cp[1], cp[2], r.data()),
        "a kriging system is singular: two points coincide without a nugget."));
  double s2 = 0, s1 = 0;
  for (double e : r) {
    s2 += e * e;
    s1 += std::fabs(e);
  }
  out = vobj();
  put(out, "mspe", vnum(s2 / p.n));
  put(out, "rmspe", vnum(std::sqrt(s2 / p.n)));
  put(out, "mae", vnum(s1 / p.n));
  put(out, "resid", vauto(r));
  return true;
}

bool trsurf(Ctx &c, Val &out) {
  Pts p;
  Dv v;
  CK(xyv(c, p, v));
  int o = 0;
  CK(wholed(c, "order", 2, 0, o));
  if (o > 3) return c.fail("`order` must be 0, 1, 2 or 3.");
  static const int np[] = {1, 3, 6, 10};
  const int k = np[o];
  if (p.n < k) return c.fail("order %d needs at least %d points.", o, k);
  static const char *const cn[] = {"1", "x", "y", "x2", "y2", "xy", "x3", "y3", "x2y", "xy2"};
  double b[10], s[10], r2 = 0;
  CK(kc(c, clab_trsurf(p.x.data(), p.y.data(), p.n, v.p(), o, b, s, &r2),
        "the monomials are collinear at these points."));
  out = vobj();
  put(out, "coef", nvec(k, cn, b));
  put(out, "se", nvec(k, cn, s));
  put(out, "r2", vnum(r2));
  put(out, "order", vnum(o));
  put(out, "n", vnum(p.n));
  return true;
}

// the spatial linear models: sglmx (kind 0) and smixm (kind 1)
bool slm(Ctx &c, Val &out, int kind) {
  Pts p;
  CK(xy(c, c["x"], c["y"], PLANAR, p));
  const int n = p.n;
  Dv v;
  CK(num(c, c["v"], "v", v, n));
  int k = 0;
  std::vector<double> z;
  CK(zmat(c, "z", n, true, k, z));
  if (n <= k) return c.fail("`v` must have more values than `z` has columns.");
  if (!kutil::span(p.x.data(), p.y.data(), n)) return c.fail("the points are too far apart: their squared distances overflow a double; rescale `x` and `y`.");
  std::vector<double> b(zz(k)), s(zz(k));
  double o[4];
  const int rc = kind == 0 ? clab_sglmx(z.data(), n, k, v.p(), p.x.data(), p.y.data(), b.data(), s.data(), o)
                           : clab_smixm(z.data(), n, k, v.p(), p.x.data(), p.y.data(), b.data(), s.data(), o);
  CK(kc(c, rc, "the columns of `z` are collinear."));
  std::vector<std::string> cn;
  for (int j = 0; j < k; ++j) cn.push_back(c.has("z") ? "z" + std::to_string(j + 1) : "1");
  out = vobj();
  put(out, "coef", nvec(cn, b.data()));
  put(out, "se", nvec(cn, s.data()));
  if (kind == 0) {
    put(out, "sigma2", vnum(o[0]));
    put(out, "phi", vnum(o[1]));
    put(out, "tau2", vnum(o[2]));
  } else {
    put(out, "sigma2", vnum(o[0]));
    put(out, "tau2", vnum(o[1]));
    put(out, "phi", vnum(o[2]));
    put(out, "conv", vbool(o[3] != 0));
  }
  put(out, "n", vnum(n));
  return true;
}
bool sglmx(Ctx &c, Val &out) { return slm(c, out, 0); }
bool smixm(Ctx &c, Val &out) { return slm(c, out, 1); }

bool spobj(Ctx &c, Val &out) {
  const Val &d = c["data"];
  if (!isdf(d)) return c.fail("`data` must be a data frame.");
  std::vector<std::string> cs = {"lon", "lat"};
  if (c.has("coords")) {
    const Val &k = c["coords"];
    cs.clear();
    if (k.t == Val::ARR) {
      for (const Val &e : k.a) cs.push_back(e.t == Val::STR ? e.s : std::string(1, '\0'));
    }
    if (k.t != Val::ARR || cs.size() != 2) return c.fail("`coords` must name two columns of `data`.");
  }
  for (const std::string &s : cs) {
    if (col(d, s.c_str()).t == Val::NUL) return c.fail("`coords` must name two columns of `data`.");
  }
  double crs = 4326;
  if (c.has("crs")) {
    crs = val1(c["crs"]);
    if (std::isnan(crs) || crs != std::floor(crs) || crs < 1) {
      return c.fail("`crs` must be one EPSG code, a whole number.");
    }
  }
  const Val &cx = col(d, cs[0].c_str()), &cy = col(d, cs[1].c_str());
  if (!isnum(cx) || !isnum(cy)) return c.fail("the coordinates must be finite numbers.");
  const std::vector<double> x = dbl(cx), y = dbl(cy);
  for (std::size_t i = 0; i < x.size(); ++i) {
    if (!std::isfinite(x[i]) || !std::isfinite(y[i])) return c.fail("the coordinates must be finite numbers.");
  }
  const std::size_t n = x.size();
  std::vector<double> m(n * 2);
  for (std::size_t i = 0; i < n; ++i) {
    m[i] = x[i];
    m[n + i] = y[i];
  }
  Val rest = vobj();
  for (std::size_t j = 0; j < d.k.size(); ++j) {
    if (d.k[j] != cs[0] && d.k[j] != cs[1]) put(rest, d.k[j].c_str(), d.a[j]);
  }
  out = vobj();
  put(out, "coords", vmat(m, static_cast<int>(n), 2));
  put(out, "data", rest);
  put(out, "cnames", vstrs(cs));
  put(out, "crs", vnum(crs));
  put(out, "units", vstr(crs == 4326 ? "degrees" : "metres"));
  put(out, "n", vnum(static_cast<double>(n)));
  return true;
}

// the coordinates (column-major) and EPSG code of an spobj() object
bool sobj(Ctx &c, std::vector<double> &m, int &n, double &crs) {
  const Val &s = c["x"];
  const Val &co = col(s, "coords");
  if (s.t != Val::OBJ || co.t == Val::NUL) return c.fail("`x` must be the object spobj() returns.");
  int dim[2] = {-1, -1};
  if (!(co.t == Val::ARR && !co.a.empty() && co.a[0].t == Val::ARR && shape(co, 2, dim, 0) && dim[1] == 2)) {
    return c.fail("`x$coords` must be a matrix of two columns.");
  }
  n = dim[0];
  CK(nmat(c, co, "x$coords", dim[0], dim[1], m, true));
  crs = val1(col(s, "crs"));
  return true;
}

bool spxy(Ctx &c, Val &out) {
  const Val &s = c["x"];
  if (s.t != Val::OBJ || col(s, "coords").t == Val::NUL) {
    return c.fail("`x` must be the object spobj() returns.");
  }
  out = col(s, "coords");
  return true;
}

bool crsof(Ctx &c, Val &out) {
  const Val &s = c["x"];
  if (s.t != Val::OBJ || col(s, "coords").t == Val::NUL) {
    return c.fail("`x` must be the object spobj() returns.");
  }
  const double crs = val1(col(s, "crs"));
  out = vobj();
  put(out, "epsg", col(s, "crs"));
  put(out, "wkt", vstr(crs == 4326   ? "GEOGCRS[\"WGS 84\", EPSG:4326]"
                       : crs == 3857 ? "PROJCRS[\"WGS 84 / Pseudo-Mercator\", EPSG:3857]"
                                     : "EPSG:" + rfmt(crs)));
  return true;
}

bool lfspec(Ctx &c, Val &out) {
  std::vector<double> m;
  int n = 0;
  double crs = 0;
  CK(sobj(c, m, n, crs));
  if (crs != 4326 && crs != 3857) return c.fail("only EPSG:4326 and EPSG:3857 can be shown on a web map.");
  std::vector<double> ll(m.size());
  if (crs == 3857) {
    CK(kc(c, clab_wmerc(m.data(), m.data() + n, n, 1, ll.data(), ll.data() + n),
          "the coordinates must be finite."));
  } else {
    ll = m;
  }
  const double ctr[2] = {nmean(ll.data(), zz(n)), nmean(ll.data() + n, zz(n))};
  static const char *const cn[] = {"lon", "lat"};
  out = vobj();
  put(out, "provider", vstr("OpenStreetMap"));
  put(out, "interactive", vbool(true));
  put(out, "center", nvec(2, cn, ctr));
  put(out, "zoom", vnum(8));
  put(out, "points", vmat(ll, n, 2));
  return true;
}

bool spsamp(Ctx &c, Val &out) {
  std::vector<double> m;
  int n = 0;
  double crs = 0;
  CK(sobj(c, m, n, crs));
  int k = 0;
  CK(whole(c, c["size"], "size", 0, k));
  if (k > n) return c.fail("`size` must be at most the number of points, %d.", n);
  bool bd = true;
  CK(flagd(c, "dens", true, bd));
  unsigned int sd = 0;
  CK(seed(c, sd));
  std::vector<double> pr(zz(n), 1.0 / n);
  if (bd) {
    if (n < 2) return c.fail("`x` must hold at least 2 points to sample by density.");
    CK(kc(c, clab_spdens(m.data(), m.data() + n, n, crs == 4326, pr.data()), "the points all coincide."));
  }
  std::vector<int> idx(zz(k));
  CK(kc(c, clab_spsamp(bd ? pr.data() : nullptr, n, k, sd, idx.data()), "the kernel refused the input."));
  out = vints(idx.data(), zz(k));
  return true;
}
