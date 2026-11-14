// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Fairness and the audits of risk tools and predictive policing: group
// metrics, disparities over time and place, the Lipschitz audit, base
// rates behind noisy records, noisy rankings, patrol detection, simulated
// records and the Oaxaca-Blinder decomposition.

// A frame of k columns of g rows from a block of k * g values, after a
// label column.
Val gframe(const char *first, const std::vector<std::string> &lab, const char *const *nm, int k, int g,
           const double *src) {
  Frame f;
  f.str(first, lab);
  for (int t = 0; t < k; ++t) f.num(nm[t], src + static_cast<std::size_t>(t) * g, static_cast<std::size_t>(g));
  return f.o;
}

// glab() and b01(): labels and 0/1 values given as a vector.
bool glab(Ctx &c, const Val &o, const char *nm, std::size_t n, Lab &L) { return codes(c, o, nm, n, L); }

bool fair(Ctx &c, Val &out) {
  Args a(c, {"score", "y", "grp"});
  CK(a.byname(c));
  Dv s;
  CK(num(c, a[0], "score", s));
  const std::size_t n = s.n();
  std::vector<int> y;
  CK(bin01(c, a[1], "y", n, y));
  Lab g;
  CK(codes(c, a[2], "grp", n, g));
  const double ct = val1d(c, "cut", 0.5);
  if (std::isnan(ct)) return c.fail("`cut` must be one number.");
  Val d1;
  int b = 0;
  CK(rref(c, arg(c, "base", d1, 1), g.lev, "ref", "group", b));
  double z = 0.0;
  CK(zconf(c, "conf", 0.95, z));
  const int G = g.g;
  std::vector<double> w(static_cast<std::size_t>(41) * G);
  double *cnt = w.data(), *met = cnt + 5 * G, *rat = met + 10 * G, *lo = rat + 10 * G, *hi = lo + 8 * G;
  if (clab_fair(s.p(), y.data(), g.code.data(), static_cast<int>(n), G, ct, b - 1, z, cnt, met, rat, lo, hi) !=
      CLAB_OK) {
    return c.fail("the kernel refused the input.");
  }
  static const char *const cn[] = {"tp", "fp", "tn", "fn", "n"};
  static const char *const mn[] = {"acc", "tpr", "tnr", "fpr", "fnr", "ppv", "npv", "ppos", "npos", "mcc"};
  const double *src[5] = {cnt, met, rat, lo, hi};
  const int wd[5] = {5, 10, 10, 8, 8};
  static const char *const nm[] = {"cnt", "met", "rat", "lo", "hi"};
  out = vobj();
  for (int t = 0; t < 5; ++t) put(out, nm[t], gframe("group", g.lev, t == 0 ? cn : mn, wd[t], G, src[t]));
  return true;
}

bool fairdf(Ctx &c, Val &out) {
  Args a(c, {"pred", "grp", "y"});
  CK(a.byname(c));
  const std::size_t n = vlen(a[0]);
  if (n < 1) return c.fail("`pred` must hold at least one case.");
  std::vector<int> pv, yv;
  CK(bin01(c, a[0], "pred", n, pv));
  Lab g;
  CK(glab(c, a[1], "grp", n, g));
  if (g.g < 2) return c.fail("`grp` must hold at least two groups.");
  const bool hy = a[2].t != Val::NUL;
  if (hy) CK(bin01(c, a[2], "y", n, yv));
  int b = 0;
  if (c.has("base")) CK(rref(c, c["base"], g.lev, "base", "group", b));
  double z = 0.0;
  CK(zconf(c, "conf", 0.95, z));
  int B = 0;
  CK(wholed(c, "nboot", 0, 0, B));
  unsigned int s = 0;
  if (B > 0) CK(seed(c, s));
  const int G = g.g;
  std::vector<double> w(static_cast<std::size_t>(40) * G + 11 + 2 * static_cast<std::size_t>(B));
  double *cnt = w.data(), *met = cnt + 6 * G, *sum = met + 23 * G, *wk = sum + 11;
  if (clab_fairdf(pv.data(), hy ? yv.data() : nullptr, g.code.data(), static_cast<int>(n), G, b - 1, z, B, s, cnt,
                  met, sum, wk) != CLAB_OK) {
    return c.fail("the kernel refused the input.");
  }
  static const char *const mn[] = {"n",    "npos", "ppos",   "lo",     "hi",   "dp",     "dplo",   "dphi",
                                   "di",   "dilo", "dihi",   "tpr",    "fpr",  "dtpr",   "dtprlo", "dtprhi",
                                   "dfpr", "dfprlo", "dfprhi", "aod", "aodlo", "aodhi", "eo"};
  static const char *const sn[] = {"di", "dp", "gini", "bas", "ginilo", "ginihi", "baslo", "bashi", "aod", "eo"};
  out = vobj();
  put(out, "tab", gframe("group", g.lev, mn, hy ? 23 : 11, G, met));
  put(out, "sum", nvec(hy ? 10 : 8, sn, sum));
  put(out, "base", vstr(g.lev[static_cast<std::size_t>(sum[10])]));
  return true;
}

bool tmpaud(Ctx &c, Val &out) {
  Args a(c, {"pred", "grp", "period", "city"});
  CK(a.byname(c));
  const std::size_t n = vlen(a[0]);
  if (n < 1) return c.fail("`pred` must hold at least one case.");
  std::vector<int> pv;
  CK(bin01(c, a[0], "pred", n, pv));
  Lab g, per, ct;
  CK(glab(c, a[1], "grp", n, g));
  if (g.g < 2) return c.fail("`grp` must hold at least two groups.");
  CK(glab(c, a[2], "period", n, per));
  if (a[3].t == Val::NUL) {
    ct.code.assign(n, 0);
    ct.lev = {"all"};
    ct.g = 1;
  } else {
    CK(glab(c, a[3], "city", n, ct));
  }
  std::vector<long long> key(n);
  for (std::size_t i = 0; i < n; ++i) key[i] = static_cast<long long>(ct.code[i]) * per.g + per.code[i];
  std::vector<long long> u(key);
  std::sort(u.begin(), u.end());
  u.erase(std::unique(u.begin(), u.end()), u.end());
  const int nc = static_cast<int>(u.size());
  std::vector<int> cell(n), cy(static_cast<std::size_t>(nc));
  for (std::size_t i = 0; i < n; ++i) {
    cell[i] = static_cast<int>(std::lower_bound(u.begin(), u.end(), key[i]) - u.begin());
  }
  for (int k = 0; k < nc; ++k) cy[static_cast<std::size_t>(k)] = static_cast<int>(u[static_cast<std::size_t>(k)] / per.g);
  int b = 0;
  if (c.has("base")) CK(rref(c, c["base"], g.lev, "base", "group", b));
  const int ncity = ct.g;
  std::vector<double> w(6 * static_cast<std::size_t>(nc) + 9 * static_cast<std::size_t>(ncity) + 3 +
                        2 * static_cast<std::size_t>(nc) * g.g + 2 * static_cast<std::size_t>(g.g));
  double *co = w.data(), *yo = co + 6 * static_cast<std::size_t>(nc), *sum = yo + 9 * ncity, *wk = sum + 3;
  if (clab_tmpaud(pv.data(), g.code.data(), cell.data(), static_cast<int>(n), g.g, nc, cy.data(), ncity, b - 1, co,
                  yo, sum, wk) != CLAB_OK) {
    return c.fail("the kernel refused the input.");
  }
  Frame fc;
  std::vector<std::string> cn, pn;
  std::vector<int> cnn(static_cast<std::size_t>(nc)), ok(static_cast<std::size_t>(nc));
  std::vector<double> v[4];
  for (int k = 0; k < nc; ++k) {
    const long long kk = u[static_cast<std::size_t>(k)];
    cn.push_back(ct.lev[static_cast<std::size_t>(kk / per.g)]);
    pn.push_back(per.lev[static_cast<std::size_t>(kk % per.g)]);
    cnn[static_cast<std::size_t>(k)] = static_cast<int>(co[6 * k]);
    ok[static_cast<std::size_t>(k)] = co[6 * k + 1] > 0.0;
    for (int t = 0; t < 4; ++t) v[t].push_back(co[6 * k + 2 + t]);
  }
  fc.str("city", cn);
  fc.str("period", pn);
  fc.ints("n", cnn.data(), cnn.size());
  fc.lgl("ok", ok.data(), ok.size());
  static const char *const vn[] = {"di", "dp", "gini", "bas"};
  for (int t = 0; t < 4; ++t) fc.num(vn[t], v[t]);
  Frame fy;
  fy.str("city", ct.lev);
  std::vector<int> ncl, gt;
  std::vector<double> yv[7];
  for (int k = 0; k < ncity; ++k) {
    ncl.push_back(static_cast<int>(yo[9 * k]));
    for (int t = 0; t < 7; ++t) yv[t].push_back(yo[9 * k + 1 + t]);
    gt.push_back(static_cast<int>(yo[9 * k + 8]));
  }
  fy.ints("ncell", ncl.data(), ncl.size());
  static const char *const yn[] = {"di", "dp", "gini", "bas", "dimin", "dimax", "dirng"};
  for (int t = 0; t < 7; ++t) fy.num(yn[t], yv[t]);
  fy.ints("ngt1", gt.data(), gt.size());
  static const char *const sn[] = {"rng", "spread"};
  out = vobj();
  put(out, "cells", fc.o);
  put(out, "cities", fy.o);
  put(out, "sum", nvec(2, sn, sum));
  put(out, "base", vstr(g.lev[static_cast<std::size_t>(sum[2])]));
  return true;
}

bool lipaud(Ctx &c, Val &out) {
  Dv hv, lv;
  CK(num(c, c["h"], "h", hv));
  const int n = static_cast<int>(hv.n());
  Val d1;
  CK(num(c, arg(c, "lip", d1, 1), "lip", lv, 1, 0.0));
  int m = 0;
  std::vector<int> pi, pj;
  std::vector<double> dv;
  const bool hd = c.has("d");
  if (!c.has("pairs")) {
    if (!hd) {
      return c.fail("give `d`, the distances: a dist object or a square matrix, or `pairs` with one distance each.");
    }
    std::vector<double> dd;
    int dim[2] = {-1, -1};
    if (!shape(c["d"], 2, dim, 0) || dim[0] != n || dim[1] != n) {
      return c.fail("`d` must be a dist object or a square matrix of distances between the %d cases.", n);
    }
    dd.assign(static_cast<std::size_t>(n) * static_cast<std::size_t>(n), 0.0);
    bool bad = false;
    fill(c["d"], 2, dim, 0, 0, 1, dd, bad, true);
    if (n < 2) return c.fail("`h` must hold at least two cases.");
    if (0.5 * n * (n - 1.0) > INT_MAX) return c.fail("`h` holds too many cases for every pair.");
    m = n * (n - 1) / 2;
    for (int j = 0; j < n; ++j) {
      for (int i = j + 1; i < n; ++i) {
        pi.push_back(j);
        pj.push_back(i);
        const double e = dd[static_cast<std::size_t>(i) + static_cast<std::size_t>(n) * j];
        if (!(e == dd[static_cast<std::size_t>(j) + static_cast<std::size_t>(n) * i])) {
          return c.fail("`d` must be symmetric.");
        }
        if (!(e >= 0.0) || !std::isfinite(e)) return c.fail("`d` must hold finite distances, at least 0.");
        dv.push_back(e);
      }
    }
  } else {
    int r = 0, k2 = 0;
    std::vector<double> pm;
    CK(nmat(c, c["pairs"], "pairs", r, k2, pm));
    if (k2 != 2 || r < 1) return c.fail("`pairs` must be a matrix of two columns, one row per pair.");
    m = r;
    pi.assign(static_cast<std::size_t>(m), 0);
    pj.assign(static_cast<std::size_t>(m), 0);
    for (int k = 0; k < 2 * m; ++k) {
      const double e = pm[static_cast<std::size_t>(k)];
      if (e != std::floor(e) || e < 1 || e > n) return c.fail("`pairs` must hold case numbers from 1 to %d.", n);
      (k < m ? pi : pj)[static_cast<std::size_t>(k % m)] = static_cast<int>(e) - 1;
    }
    if (hd) {
      Dv d;
      CK(num(c, c["d"], "d", d, m, 0.0));
      dv = d.v;
    }
  }
  std::vector<double> gap(static_cast<std::size_t>(m)), rat(static_cast<std::size_t>(m));
  std::vector<int> vi(static_cast<std::size_t>(m));
  double sum[6];
  if (clab_lipaud(hv.p(), n, pi.data(), pj.data(), dv.empty() ? nullptr : dv.data(), m, lv[0], gap.data(), rat.data(),
                  vi.data(), sum) != CLAB_OK) {
    return c.fail("the kernel refused the input.");
  }
  std::vector<int> ci(static_cast<std::size_t>(m)), cj(static_cast<std::size_t>(m));
  std::vector<double> cd(static_cast<std::size_t>(m));
  for (std::size_t k = 0; k < static_cast<std::size_t>(m); ++k) {
    ci[k] = pi[k] + 1;
    cj[k] = pj[k] + 1;
    cd[k] = dv.empty() ? 1.0 : dv[k];
  }
  Frame f;
  f.ints("i", ci.data(), ci.size());
  f.ints("j", cj.data(), cj.size());
  f.num("d", cd);
  f.num("gap", gap);
  f.num("ratio", rat);
  f.lgl("viol", vi.data(), vi.size());
  const std::size_t kr = static_cast<std::size_t>(sum[5]);
  out = vobj();
  put(out, "lreq", vnum(sum[0]));
  put(out, "nviol", vnum(sum[1]));
  put(out, "rate", vnum(sum[2]));
  put(out, "ok", vbool(sum[1] == 0.0));
  put(out, "maxgap", vnum(sum[3]));
  const int wp[2] = {pi[kr] + 1, pj[kr] + 1};
  put(out, "worst", vints(wp, 2));
  put(out, "pairs", f.o);
  return true;
}

// box(): one number a for [0, a], or two, lo <= hi, all from 0.
bool box(Ctx &c, const char *nm, double &lo, double &hi) {
  Dv v;
  CK(num(c, c[nm], nm, v, -1, 0.0));
  if (v.n() > 2) return c.fail("`%s` must be one number (the box [0, %s]) or two (its ends).", nm, nm);
  lo = v.n() == 1 ? 0.0 : v[0];
  hi = v[v.n() - 1];
  if (lo > hi) return c.fail("`%s` must give its lower end first.", nm);
  return true;
}

bool nsbnd(Ctx &c, Val &out) {
  Dv pv;
  CK(num(c, c["pobs"], "pobs", pv, -1, 0.0));
  const int n = static_cast<int>(pv.n());
  for (double x : pv.v) {
    if (x > 1.0) return c.fail("`pobs` must be rates from 0 to 1.");
  }
  double alo = 0.0, ahi = 0.0, blo = 0.0, bhi = 0.0;
  CK(box(c, "a", alo, ahi));
  CK(box(c, "b", blo, bhi));
  if (!(ahi + bhi < 1.0)) return c.fail("the upper ends of `a` and `b` must add to less than 1.");
  if (0.5 * n * (n - 1.0) > INT_MAX) return c.fail("`pobs` holds too many rates for every pair.");
  const int m = n * (n - 1) / 2;
  const std::vector<std::string> lab = nmsof(pv.nm, pv.n());
  std::vector<double> lo(pv.n()), hi(pv.n()), wd(pv.n()), bd(static_cast<std::size_t>(m));
  std::vector<int> ord(static_cast<std::size_t>(m));
  if (clab_nsbnd(pv.p(), n, alo, ahi, blo, bhi, lo.data(), hi.data(), ord.data(), bd.data()) != CLAB_OK) {
    return c.fail("the kernel refused the input.");
  }
  for (std::size_t i = 0; i < pv.n(); ++i) wd[i] = hi[i] - lo[i];
  Frame f;
  f.str("rate", lab);
  f.num("pobs", pv.v);
  f.num("lo", lo);
  f.num("hi", hi);
  f.num("width", wd);
  Val li = varr(), lj = varr(), os = varr();
  static const char *const ow[] = {"higher", "undecided", "lower"};
  std::size_t k = 0;
  for (std::size_t i = 0; i < pv.n(); ++i) {
    for (std::size_t j = i + 1; j < pv.n(); ++j, ++k) {
      li.a.push_back(vstr(lab[i]));
      lj.a.push_back(vstr(lab[j]));
      const int o = ord[k];
      os.a.push_back(o == 2 ? Val() : vstr(ow[o + 1]));
    }
  }
  Frame e;
  e.add("i", li);
  e.add("j", lj);
  e.add("ord", os);
  e.num("bd", bd);
  out = vobj();
  put(out, "bnd", f.o);
  put(out, "pairs", e.o);
  return true;
}

bool rnkres(Ctx &c, Val &out) {
  Dv ev, hv;
  CK(num(c, c["est"], "est", ev));
  const int n = static_cast<int>(ev.n());
  if (n < 2) return c.fail("`est` must hold at least two estimates.");
  if (static_cast<double>(n) * n > INT_MAX) return c.fail("`est` holds too many estimates.");
  CK(num(c, c["hw"], "hw", hv, -1, 0.0));
  if (hv.n() != 1 && hv.n() != ev.n()) return c.fail("`hw` must be one half-width or one per estimate.");
  std::vector<double> hh(ev.n());
  for (std::size_t i = 0; i < ev.n(); ++i) hh[i] = hv[hv.n() == 1 ? 0 : i];
  std::vector<int> id(ev.n() * ev.n());
  double o[4];
  if (clab_rnkres(ev.p(), hh.data(), n, id.data(), o) != CLAB_OK) return c.fail("the kernel refused the input.");
  out = vobj();
  put(out, "id", vlmat(id.data(), n, n));
  put(out, "npair", vnum(o[0]));
  put(out, "nunid", vnum(o[1]));
  put(out, "share", vnum(o[2]));
  put(out, "res", vnum(o[3]));
  return true;
}

// pts2(): planar points, a two-column matrix (an array of rows) or an
// object with columns x and y (or two columns); n of them into px, py.
bool pts2(Ctx &c, const Val &x, const char *nm, std::vector<double> &px, std::vector<double> &py, bool empty) {
  px.clear();
  py.clear();
  if (x.t == Val::NUL && empty) return true;
  Val cx, cy;
  int dim[2] = {-1, -1};
  if (x.t == Val::ARR && shape(x, 2, dim, 0) && dim[1] == 2) {
    cx = varr();
    cy = varr();
    for (const Val &r : x.a) {
      cx.a.push_back(r.a[0]);
      cy.a.push_back(r.a[1]);
    }
  } else if (x.t == Val::OBJ) {
    cx = col(x, "x");
    cy = col(x, "y");
    if ((cx.t == Val::NUL || cy.t == Val::NUL) && x.a.size() == 2) {
      cx = x.a[0];
      cy = x.a[1];
    }
    if (cx.t == Val::NUL || cy.t == Val::NUL) return c.fail("`%s` must have columns x and y.", nm);
  } else {
    return c.fail("`%s` must be a two-column numeric matrix or a data frame with columns x and y.", nm);
  }
  if (!isnum(cx) || !isnum(cy) || vlen(cx) != vlen(cy)) {
    return c.fail("`%s` must hold numeric coordinates, x and y of one length.", nm);
  }
  if (vlen(cx) == 0) {
    if (!empty) return c.fail("`%s` must hold at least one point.", nm);
    return true;
  }
  px = dbl(cx);
  py = dbl(cy);
  for (std::size_t i = 0; i < px.size(); ++i) {
    if (!std::isfinite(px[i]) || !std::isfinite(py[i])) return c.fail("`%s` must hold finite coordinates.", nm);
  }
  return true;
}

bool nordet(Ctx &c, Val &out) {
  std::vector<double> cx, cy, ox, oy;
  CK(pts2(c, c["x"], "x", cx, cy, false));
  CK(pts2(c, c["off"], "off", ox, oy, true));
  Dv r;
  CK(num(c, c["radius"], "radius", r, 1, 0.0, true));
  const double pd = val1d(c, "p", 0.85);
  if (!(pd > 0.0 && pd <= 1.0)) return c.fail("`p` must be one number in (0, 1].");
  bool dr = false;
  CK(flagd(c, "draw", false, dr));
  unsigned int s = 0;
  if (dr) CK(seed(c, s));
  const std::size_t n = cx.size(), m = ox.size();
  std::vector<int> k(n), det(n);
  std::vector<double> pr(n);
  if (clab_nordet(cx.data(), cy.data(), static_cast<int>(n), m ? ox.data() : nullptr, m ? oy.data() : nullptr,
                  static_cast<int>(m), r[0], pd, dr ? 1 : 0, s, k.data(), pr.data(), dr ? det.data() : nullptr) !=
      CLAB_OK) {
    return c.fail("the kernel refused the input.");
  }
  Frame f;
  f.ints("k", k.data(), n);
  f.num("prob", pr);
  if (dr) f.ints("det", det.data(), n);
  out = f.o;
  return true;
}

bool fairsim(Ctx &c, Val &out) {
  int N = 0;
  CK(wholed(c, "n", 2000, 1, N));
  std::vector<std::string> gl = {"A", "B"};
  if (c.has("groups")) {
    if (!alabs(c["groups"], gl) || gl.size() < 2) {
      return c.fail("`groups` must be two or more labels with no missing values.");
    }
  }
  const int G = static_cast<int>(gl.size());
  for (int i = 0; i < G; ++i) {
    for (int j = 0; j < i; ++j) {
      if (gl[static_cast<std::size_t>(i)] == gl[static_cast<std::size_t>(j)]) {
        return c.fail("`groups` must not repeat a label.");
      }
    }
  }
  Dv pv;
  if (c.has("props")) {
    CK(num(c, c["props"], "props", pv, G, 0.0));
    double t = 0.0;
    for (double x : pv.v) t += x;
    if (!(t > 0.0)) return c.fail("`props` must not all be 0.");
  }
  int na = 0;
  CK(wholed(c, "nareas", 20, 1, na));
  if (na < G) return c.fail("`nareas` must be at least the number of groups, %d.", G);
  double b0 = 0.0;
  CK(p01d(c, "base", 0.3, b0));
  const double bs = val1d(c, "bias", 0.5);
  if (!(bs >= -1.0 && bs <= 1.0)) return c.fail("`bias` must be one number from -1 to 1.");
  unsigned int s;
  CK(seed(c, s));
  const std::size_t n = static_cast<std::size_t>(N);
  std::vector<int> ai(n), gi(n), tr(n), dt(n);
  std::vector<double> rk(n);
  if (clab_fairsim(N, G, pv.n() ? pv.p() : nullptr, na, b0, bs, s, ai.data(), gi.data(), tr.data(), dt.data(),
                   rk.data()) != CLAB_OK) {
    return c.fail("the kernel refused the input.");
  }
  for (int &e : ai) ++e;
  Frame f;
  f.ints("area", ai.data(), n);
  f.add("group", vstrs(gl, gi));
  f.ints("truth", tr.data(), n);
  f.ints("det", dt.data(), n);
  f.num("risk", rk);
  out = f.o;
  return true;
}

bool oaxbd(Ctx &c, Val &out) {
  Args a(c, {"y", "grp"});
  CK(a.byname(c));
  Dv yv;
  CK(num(c, a[0], "y", yv));
  const int n = static_cast<int>(yv.n());
  Lab g;
  CK(glab(c, a[1], "grp", yv.n(), g));
  if (g.g != 2) return c.fail("`grp` must hold exactly two groups.");
  Val d1;
  int r = 0;
  CK(rref(c, arg(c, "ref", d1, 1), g.lev, "ref", "group", r));
  --r;
  double z = 0.0;
  CK(zconf(c, "conf", 0.95, z));
  static const char *const fmsg =
      "`x` must be numeric covariates, one row per case: code a factor as indicator columns (how the unexplained "
      "part splits across them then depends on the category left out).";
  int p = 0;
  std::vector<double> xv;
  std::vector<std::string> tn;
  if (c.has("x")) {
    const Val &x = c["x"];
    int dim[2] = {-1, -1};
    if (x.t == Val::OBJ) {
      p = static_cast<int>(x.a.size());
      for (const Val &e : x.a) {
        if (!isnum(e) || vlen(e) != yv.n() || e.t == Val::OBJ) return c.fail("%s", fmsg);
        const std::vector<double> v = dbl(e);
        xv.insert(xv.end(), v.begin(), v.end());
      }
      tn = x.k;
    } else if (isnum(x)) {
      if (vlen(x) != yv.n()) return c.fail("%s", fmsg);
      p = 1;
      xv = dbl(x);
    } else if (x.t == Val::ARR && shape(x, 2, dim, 0)) {
      if (dim[0] != n) return c.fail("%s", fmsg);
      int r2 = 0, c2 = 0;
      CK(nmat(c, x, "x", r2, c2, xv, true));
      p = c2;
    } else {
      return c.fail("%s", fmsg);
    }
    for (double e : xv) {
      if (!std::isfinite(e)) return c.fail("`x` must hold finite numbers with no missing values.");
    }
  }
  const int q = p + 1;
  int nn[2] = {0, 0};
  for (int e : g.code) ++nn[e];
  if (nn[0] <= q || nn[1] <= q) return c.fail("each group needs more cases than coefficients, %d.", q);
  std::vector<std::string> lab = {"const"};
  for (int j = 0; j < p; ++j) {
    if (!tn.empty() && !tn[static_cast<std::size_t>(j)].empty()) {
      lab.push_back(tn[static_cast<std::size_t>(j)]);
    } else {
      char b[24];
      std::snprintf(b, sizeof b, p == 1 && tn.empty() ? "x" : "x%d", j + 1);
      lab.push_back(b);
    }
  }
  const std::size_t Q = static_cast<std::size_t>(q), N = yv.n();
  std::vector<double> w(2 * Q * 3 + 14 + 5 * Q + Q * (N + 5 * Q + 1) + N);
  double *xm = w.data(), *bb = xm + 2 * q, *bse = bb + 2 * q, *dec = bse + 2 * q, *var = dec + 14, *wk = var + 5 * q;
  const int rc = clab_oaxbd(yv.p(), xv.empty() ? nullptr : xv.data(), g.code.data(), n, p, r, xm, bb, bse, dec, var, wk);
  if (rc == CLAB_EDOM) return c.fail("a covariate is constant, or a combination of the others, within a group.");
  if (rc != CLAB_OK) return c.fail("the kernel refused the input.");
  Val gl = vobj();
  put(gl, "a", vstr(g.lev[static_cast<std::size_t>(1 - r)]));
  put(gl, "b", vstr(g.lev[static_cast<std::size_t>(r)]));
  Frame fm, fc, fv, fd;
  fm.str("term", lab);
  fm.num("a", xm, Q);
  fm.num("b", xm + q, Q);
  fc.str("term", lab);
  fc.num("a", bb, Q);
  fc.num("ase", bse, Q);
  fc.num("b", bb + q, Q);
  fc.num("bse", bse + q, Q);
  fv.str("term", lab);
  static const char *const vn[] = {"exp", "expse", "unexp", "unexpse", "shift"};
  for (int t = 0; t < 5; ++t) fv.num(vn[t], var + static_cast<std::size_t>(t) * Q, Q);
  static const char *const pn[] = {"gap", "exp", "unexp", "expa", "unexpa", "endow", "coef", "inter"};
  const int ix[8] = {2, 4, 6, 8, 10, 4, 10, 12};
  std::vector<double> de(8), ds(8), dl(8), dh(8);
  std::vector<std::string> ps;
  for (int k = 0; k < 8; ++k) {
    ps.push_back(pn[k]);
    de[static_cast<std::size_t>(k)] = dec[ix[k]];
    ds[static_cast<std::size_t>(k)] = dec[ix[k] + 1];
    dl[static_cast<std::size_t>(k)] = de[static_cast<std::size_t>(k)] - z * ds[static_cast<std::size_t>(k)];
    dh[static_cast<std::size_t>(k)] = de[static_cast<std::size_t>(k)] + z * ds[static_cast<std::size_t>(k)];
  }
  fd.str("part", ps);
  fd.num("est", de);
  fd.num("se", ds);
  fd.num("lo", dl);
  fd.num("hi", dh);
  out = vobj();
  put(out, "groups", gl);
  put(out, "means", fm.o);
  put(out, "coef", fc.o);
  put(out, "dec", fd.o);
  put(out, "by", fv.o);
  return true;
}
