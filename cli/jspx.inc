// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Spatial methods of src/bspx.cpp: near repeats and space-time scans,
// kernel density, DBSCAN, hot-spot persistence, empirical-Bayes
// shrinkage, Moran's I and the spatial regressions with their spillovers,
// cross-type K, quadrat tests, Morisita's index, segregation indices,
// points in polygons, areal crosswalks and buffer assignment.

// edges(): increasing positive numbers.
bool edges(Ctx &c, const char *nm, Dv &v) {
  CK(num(c, c[nm], nm, v, -1, 0.0, true));
  for (std::size_t q = 1; q < v.n(); ++q) {
    if (!(v[q] > v[q - 1])) return c.fail("`%s` must increase.", nm);
  }
  return true;
}

// mnum(): the values of a numeric matrix, finite (from 0 when lo is 0).
bool mnum(Ctx &c, const Val &x, const char *nm, double lo, int &nr, int &nc, std::vector<double> &v) {
  int dim[2] = {-1, -1};
  if (!(x.t == Val::ARR && !x.a.empty() && x.a[0].t == Val::ARR && shape(x, 2, dim, 0))) {
    return c.fail("`%s` must be numeric.", nm);
  }
  nr = dim[0];
  nc = dim[1];
  v.assign(zz(nr) * zz(nc), 0.0);
  bool bad = false;
  fill(x, 2, dim, 0, 0, 1, v, bad, true);
  for (double e : v) {
    if (!std::isfinite(e) || e < lo) {
      return c.fail("`%s` must hold finite numbers%s with no missing values.", nm, lo == 0 ? " from 0" : "");
    }
  }
  return true;
}

// A matrix given as rows, or a data frame (as.matrix(): its names): rows,
// columns, values column-major and column names.
bool dfmat(Ctx &c, const Val &x, const char *nm, double lo, int &nr, int &nc, std::vector<double> &v,
           std::vector<std::string> &cn, const char *pre) {
  cn.clear();
  if (isdf(x)) {
    nc = static_cast<int>(x.a.size());
    nr = static_cast<int>(vlen(x.a[0]));
    v.clear();
    for (std::size_t j = 0; j < x.a.size(); ++j) {
      const Val &e = x.a[j];
      if (!isnum(e) || static_cast<int>(vlen(e)) != nr) return c.fail("`%s` must be numeric.", nm);
      const std::vector<double> d = dbl(e);
      for (double q : d) {
        if (!std::isfinite(q) || q < lo) {
          return c.fail("`%s` must hold finite numbers%s with no missing values.", nm, lo == 0 ? " from 0" : "");
        }
      }
      v.insert(v.end(), d.begin(), d.end());
      cn.push_back(x.k[j]);
    }
    return true;
  }
  CK(mnum(c, x, nm, lo, nr, nc, v));
  for (int j = 0; j < nc; ++j) cn.push_back(pre + std::to_string(j + 1));
  return true;
}
inline bool ismat2(const Val &x) {
  int dim[2] = {-1, -1};
  return isdf(x) || (x.t == Val::ARR && !x.a.empty() && x.a[0].t == Val::ARR && shape(x, 2, dim, 0));
}

// wread() of src/bspx.cpp: an spdep listw (neighbours and weights), an
// n x n matrix, or the list spwt() returns.
bool wreadx(Ctx &c, const Val &w, int n, Wcsr &o) {
  std::vector<int> ei, ej;
  std::vector<double> ew;
  if (w.t == Val::OBJ && col(w, "neighbours").t != Val::NUL) {
    const Val &nb = col(w, "neighbours"), &wt = col(w, "weights");
    if (nb.t != Val::ARR || wt.t != Val::ARR || static_cast<int>(nb.a.size()) != n ||
        static_cast<int>(wt.a.size()) != n) {
      return c.fail("`w` must have one row of neighbours and weights per unit.");
    }
    for (int i = 0; i < n; ++i) {
      const Val &a = nb.a[zz(i)], &b = wt.a[zz(i)];
      if (!isnum(a)) return c.fail("`w` must have integer neighbour lists.");
      const std::vector<double> av = dbl(a);
      for (double e : av) {
        if (e != std::floor(e)) return c.fail("`w` must have integer neighbour lists.");
      }
      const std::size_t k = av.size();
      if (k == 1 && av[0] == 0) continue;
      if (!isnum(b) || vlen(b) != k) return c.fail("`w` must have one weight per neighbour.");
      const std::vector<double> bv = dbl(b);
      for (std::size_t q = 0; q < k; ++q) {
        const double j = av[q];
        if (j < 1 || j > n) return c.fail("`w` names a neighbour outside 1 to %d.", n);
        ei.push_back(i);
        ej.push_back(static_cast<int>(j) - 1);
        ew.push_back(bv[q]);
      }
    }
  } else if (w.t == Val::ARR && !w.a.empty() && w.a[0].t == Val::ARR) {
    int dim[2] = {-1, -1};
    if (!shape(w, 2, dim, 0) || dim[0] != n || dim[1] != n) return c.fail("`w` must be a %d x %d matrix.", n, n);
    int nr = 0, nc = 0;
    std::vector<double> m;
    CK(mnum(c, w, "w", 0.0, nr, nc, m));
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        const double v = m[zz(i) + zz(j) * zz(n)];
        if (v == 0) continue;
        if (i == j) return c.fail("a unit cannot be its own neighbour in `w`.");
        ei.push_back(i);
        ej.push_back(j);
        ew.push_back(v);
      }
    }
  } else {
    const Val &nn = col(w, "n");
    if (nn.t == Val::NUL) return c.fail("`w` must be the list spwt() returns, an spdep listw or a matrix.");
    if (val1(nn) != n) return c.fail("`w` is for %g units, not %d.", val1(nn), n);
    Dv i, j, v;
    CK(num(c, col(w, "i"), "w$i", i, -1, 1.0));
    CK(num(c, col(w, "j"), "w$j", j, static_cast<long>(i.n()), 1.0));
    CK(num(c, col(w, "w"), "w$w", v, static_cast<long>(i.n()), 0.0));
    for (std::size_t q = 0; q < i.n(); ++q) {
      const double a = i[q], b = j[q];
      if (a > n || b > n || a != std::floor(a) || b != std::floor(b)) {
        return c.fail("`w$i` and `w$j` must be unit numbers from 1 to %d.", n);
      }
      ei.push_back(static_cast<int>(a) - 1);
      ej.push_back(static_cast<int>(b) - 1);
      ew.push_back(v[q]);
    }
  }
  for (double e : ew) {
    if (!(e >= 0) || !std::isfinite(e)) return c.fail("`w` must hold finite weights from 0.");
  }
  return wbuild(c, n, ei, ej, ew, o);
}

// The seed when nsim draws are asked for.
bool nseed(Ctx &c, int ns, unsigned int &sd) {
  sd = 0;
  return ns <= 0 || seed(c, sd);
}

bool knxnr(Ctx &c, Val &out) {
  Pts p;
  CK(xy(c, PROJ, p));
  const int n = p.n;
  if (n < 3) return c.fail("at least 3 incidents are needed.");
  Dv tv, sv, bv;
  CK(num(c, c["t"], "t", tv, n));
  CK(edges(c, "ds", sv));
  CK(edges(c, "dt", bv));
  const int ks = static_cast<int>(sv.n()), kt = static_cast<int>(bv.n());
  int ns = 0;
  CK(wholed(c, "nsim", 999, 0, ns));
  unsigned int sd = 0;
  CK(nseed(c, ns, sd));
  const std::size_t m = zz(ks) * zz(kt);
  std::vector<double> o(3 * m + 1);
  if (clab_knxnr(p.x.data(), p.y.data(), tv.p(), n, sv.p(), ks, bv.p(), kt, ns, sd, o.data()) != CLAB_OK) {
    return c.fail("the kernel refused the input.");
  }
  std::vector<double> col0[9];
  for (auto &v : col0) v.resize(m);
  const double *ob = o.data(), *ex = ob + m, *pp = ob + 2 * m;
  for (int b = 0; b < kt; ++b) {
    for (int s = 0; s < ks; ++s) {
      const std::size_t q = zz(s) + zz(b) * zz(ks);
      col0[0][q] = s ? sv[zz(s - 1)] : 0.0;
      col0[1][q] = sv[zz(s)];
      col0[2][q] = b ? bv[zz(b - 1)] : 0.0;
      col0[3][q] = bv[zz(b)];
      col0[4][q] = ob[q];
      col0[5][q] = ex[q];
      col0[6][q] = ex[q] > 0 ? ob[q] / ex[q] : NaN;
      col0[7][q] = pp[q];
      col0[8][q] = ob[q] > 0 ? clab_ppois(ob[q] - 1, ex[q], 0) : 1.0;
    }
  }
  static const char *const tn[] = {"dlo", "dhi", "tlo", "thi", "obs", "exp", "ratio", "p", "ppois"};
  Frame f;
  for (int k = 0; k < 9; ++k) f.num(tn[k], col0[k]);
  out = vobj();
  put(out, "tab", f.o);
  put(out, "n", vnum(n));
  put(out, "pairs", vnum(o[3 * m]));
  put(out, "nsim", vnum(ns));
  return true;
}

// nrd(): the normal reference bandwidth of MASS::bandwidth.nrd() / 4.
double nrd(const std::vector<double> &v) {
  const int n = static_cast<int>(v.size());
  std::vector<double> u(v);
  std::sort(u.begin(), u.end());
  auto q = [&](double pr) {
    const double h = (n - 1) * pr;
    const int lo = static_cast<int>(std::floor(h));
    return lo + 1 < n ? u[zz(lo)] + (h - lo) * (u[zz(lo + 1)] - u[zz(lo)]) : u[zz(lo)];
  };
  double m = 0, ss = 0;
  for (double e : v) m += e;
  m /= n;
  for (double e : v) ss += (e - m) * (e - m);
  const double sd = std::sqrt(ss / (n - 1));
  return 1.06 * std::min(sd, (q(0.75) - q(0.25)) / 1.34) * std::pow(n, -0.2);
}

bool kdsurf(Ctx &c, Val &out) {
  Pts p;
  CK(xy(c, PLANAR, p));
  const int n = p.n;
  if (n < 2) return c.fail("at least 2 points are needed.");
  Dv wv;
  if (c.has("w")) CK(num(c, c["w"], "w", wv, n, 0.0));
  double hh[2] = {0, 0};
  if (!c.has("h")) {
    hh[0] = nrd(p.x);
    hh[1] = nrd(p.y);
    if (!(hh[0] > 0 && hh[1] > 0)) return c.fail("the points do not spread along both axes: give `h`.");
  } else {
    Dv hv;
    CK(num(c, c["h"], "h", hv, -1, 0.0, true));
    if (hv.n() != 1 && hv.n() != 2) return c.fail("`h` must be one or two positive numbers.");
    hh[0] = hv[0];
    hh[1] = hv[hv.n() - 1];
  }
  static const char *const eo[] = {"none", "uniform", "diggle"};
  int ed = 0;
  CK(pickd(c, "edge", eo, 3, 0, ed));
  double wn[4] = {0};
  const bool hw = c.has("win");
  if (hw) CK(win(c, c["win"], wn));
  if (ed > 0 && !hw) return c.fail("an edge correction needs `win`.");
  const double *pw = wv.n() ? wv.p() : nullptr;
  static const char *const xyn[] = {"x", "y"};
  out = vobj();
  if (c.has("at")) {
    const Val &at = c["at"];
    const bool own = at.t == Val::STR || (at.t == Val::ARR && !at.a.empty() && at.a[0].t == Val::STR);
    if (own && !(at.t == Val::STR ? at.s == "points" : at.a.size() == 1 && at.a[0].s == "points")) {
      return c.fail("`at` must be NULL, \"points\" or coordinates.");
    }
    Pts q;
    if (own) q = p;
    else CK(xy(c, at, Val(), PLANAR, q));
    std::vector<double> z(zz(q.n));
    if (clab_kdsurf(p.x.data(), p.y.data(), pw, n, hh[0], hh[1], q.x.data(), q.n, q.y.data(), q.n, 0, hw ? wn : nullptr,
                    ed, own ? 1 : 0, z.data()) != CLAB_OK) {
      return c.fail("a point lies where its window holds none of its kernel.");
    }
    put(out, "x", vauto(q.x));
    put(out, "y", vauto(q.y));
    put(out, "z", vauto(z));
    put(out, "h", nvec(2, xyn, hh));
    put(out, "regions", Val());
    return true;
  }
  int ax = 0, ay = 0;
  CK(wholed(c, "nx", 100, 2, ax));
  Val dny;
  CK(whole(c, arg(c, "ny", dny, ax), "ny", 2, ay));
  if (static_cast<double>(ax) * ay > 1e8) return c.fail("`nx` times `ny` must be at most 1e8.");
  Dv cv;
  {
    Val dc = varr();
    dc.a = {vnum(0.5), vnum(0.9)};
    CK(num(c, c.has("cover") ? c["cover"] : dc, "cover", cv, -1, 0.0, true));
  }
  for (double e : cv.v) {
    if (!(e < 1)) return c.fail("`cover` must hold shares between 0 and 1.");
  }
  double b[4] = {0};
  if (hw) {
    std::copy(wn, wn + 4, b);
  } else {
    b[0] = *std::min_element(p.x.begin(), p.x.end()) - 3 * hh[0];
    b[1] = *std::max_element(p.x.begin(), p.x.end()) + 3 * hh[0];
    b[2] = *std::min_element(p.y.begin(), p.y.end()) - 3 * hh[1];
    b[3] = *std::max_element(p.y.begin(), p.y.end()) + 3 * hh[1];
  }
  std::vector<double> gx(zz(ax)), gy(zz(ay));
  for (int a = 0; a < ax; ++a) gx[zz(a)] = b[0] + a * (b[1] - b[0]) / (ax - 1);
  for (int a = 0; a < ay; ++a) gy[zz(a)] = b[2] + a * (b[3] - b[2]) / (ay - 1);
  const std::size_t m = zz(ax) * zz(ay);
  std::vector<double> z(m);
  if (clab_kdsurf(p.x.data(), p.y.data(), pw, n, hh[0], hh[1], gx.data(), ax, gy.data(), ay, 1, hw ? wn : nullptr, ed, 0,
                  z.data()) != CLAB_OK) {
    return c.fail("a point lies where its window holds none of its kernel.");
  }
  std::vector<double> zs(z);
  std::sort(zs.begin(), zs.end(), [](double a, double d) { return a > d; });
  double tot = 0;
  for (double e : zs) tot += e;
  const double ca = (b[1] - b[0]) / (ax - 1) * (b[3] - b[2]) / (ay - 1);
  std::vector<double> rc[4];
  for (double e : cv.v) {
    double run = 0;
    std::size_t q = 0;
    while (q < m && run < e * tot) run += zs[q++];
    rc[0].push_back(e);
    rc[1].push_back(zs[q > 0 ? q - 1 : 0]);
    rc[2].push_back(static_cast<double>(q) * ca);
    rc[3].push_back(static_cast<double>(q) / static_cast<double>(m));
  }
  Frame reg;
  static const char *const rn[] = {"cover", "level", "area", "share"};
  for (int k = 0; k < 4; ++k) reg.num(rn[k], rc[k]);
  put(out, "x", vauto(gx));
  put(out, "y", vauto(gy));
  put(out, "z", vmat(z, ax, ay));
  put(out, "h", nvec(2, xyn, hh));
  put(out, "regions", reg.o);
  return true;
}

bool dbscn(Ctx &c, Val &out) {
  Pts p;
  CK(xy(c, PROJ, p));
  const int n = p.n;
  double e = 0.0;
  CK(one(c, c["eps"], "eps", 0.0, Inf, e));
  int mp = 0;
  CK(wholed(c, "minpts", 5, 1, mp));
  bool bd = true;
  CK(flagd(c, "border", true, bd));
  std::vector<int> cl(zz(n)), co(zz(n));
  int k = 0;
  if (clab_dbscn(p.x.data(), p.y.data(), n, e, mp, bd ? 1 : 0, cl.data(), co.data(), &k) != CLAB_OK) {
    return c.fail("the kernel refused the input.");
  }
  std::vector<int> sz(zz(k), 0);
  for (int &v : cl) {
    if (v >= 0) {
      ++sz[zz(v)];
      ++v;
    } else {
      v = INT_MIN;
    }
  }
  out = vobj();
  put(out, "cl", unbox(vints(cl.data(), cl.size())));
  put(out, "core", unbox(vlgls(co.data(), co.size())));
  put(out, "sizes", unbox(vints(sz.data(), sz.size())));
  return true;
}

bool hsper(Ctx &c, Val &out) {
  Pts p;
  CK(xy(c, PROJ, p));
  const int n = p.n;
  Lab g;
  CK(codes(c, c["prd"], "prd", zz(n), g));
  double cs = 0.0;
  CK(one(c, c["cell"], "cell", 0.0, Inf, cs));
  const double tp = val1d(c, "top", 0.05);
  if (!(tp > 0 && tp <= 1)) return c.fail("`top` must be one share above 0, at most 1.");
  double o[2] = {0, 0};
  if (!c.has("org")) {
    o[0] = *std::min_element(p.x.begin(), p.x.end());
    o[1] = *std::min_element(p.y.begin(), p.y.end());
  } else {
    Dv ov;
    CK(num(c, c["org"], "org", ov, 2));
    o[0] = ov[0];
    o[1] = ov[1];
  }
  const int np = g.g;
  std::vector<int> cid(zz(n)), hot(zz(n) * zz(np));
  std::vector<double> ix(zz(n)), iy(zz(n)), per(9 * zz(np)), cv(3 * zz(n));
  int nc = 0;
  if (clab_hsper(p.x.data(), p.y.data(), g.code.data(), n, np, cs, o[0], o[1], tp, cid.data(), ix.data(), iy.data(), &nc,
                 hot.data(), per.data(), cv.data()) != CLAB_OK) {
    return c.fail("the cells are too many for the grid: choose a larger `cell` or an `org` nearer the points.");
  }
  static const char *const pn[] = {"n", "active", "nhalf", "shalf", "hot", "tiedin", "tiedout", "hotshr", "jacc"};
  Frame pf;
  pf.str("period", g.lev);
  for (int j = 0; j < 9; ++j) {
    std::vector<double> v(zz(np));
    for (int q = 0; q < np; ++q) v[zz(q)] = per[9 * zz(q) + zz(j)];
    pf.num(pn[j], v);
  }
  const std::size_t sc = zz(nc);
  std::vector<double> cc[5];
  for (auto &v : cc) v.resize(sc);
  for (std::size_t j = 0; j < sc; ++j) {
    cc[0][j] = o[0] + (ix[j] + 0.5) * cs;
    cc[1][j] = o[1] + (iy[j] + 0.5) * cs;
    for (int q = 0; q < 3; ++q) cc[2 + q][j] = cv[3 * j + zz(q)];
  }
  back(p, cc[0].data(), cc[1].data(), nc);
  static const char *const cn[] = {"x", "y", "n", "nhot", "run"};
  Frame cf;
  for (int k = 0; k < 5; ++k) cf.num(cn[k], cc[k]);
  std::vector<int> hm(sc * zz(np));
  for (int q = 0; q < np; ++q) {
    for (std::size_t j = 0; j < sc; ++j) hm[j + zz(q) * sc] = hot[j + zz(q) * sc] != 0;
  }
  for (int &v : cid) ++v;
  out = vobj();
  put(out, "periods", pf.o);
  put(out, "cells", cf.o);
  put(out, "hot", vlmat(hm.data(), nc, np));
  put(out, "cell", unbox(vints(cid.data(), cid.size())));
  put(out, "k", vnum(std::max(1, static_cast<int>(std::ceil(tp * nc)))));
  return true;
}

bool ebshr(Ctx &c, Val &out) {
  Dv yv;
  CK(num(c, c["y"], "y", yv));
  const int n = static_cast<int>(yv.n());
  if (n < 2) return c.fail("at least 2 places are needed.");
  Dv wv;
  if (c.has("w")) CK(num(c, c["w"], "w", wv, n, 0.0));
  const double *pw = wv.n() ? wv.p() : nullptr;
  double nv = 0.0;
  if (!c.has("v")) {
    double sw = 0, sy = 0;
    for (int i = 0; i < n; ++i) {
      sw += pw ? pw[i] : 1.0;
      sy += (pw ? pw[i] : 1.0) * yv[zz(i)];
    }
    if (!(sw > 0)) return c.fail("`w` must have a positive sum.");
    nv = sy / sw;
    if (!(nv >= 0)) return c.fail("the mean of `y` is negative: give the noise variance `v`.");
  } else {
    Dv vv;
    CK(num(c, c["v"], "v", vv, 1, 0.0));
    nv = vv[0];
  }
  double bb = NaN;
  if (c.has("b")) CK(p01(c, c["b"], "b", bb));
  double fit[6] = {0};
  std::vector<double> sh(zz(n)), fa(zz(n));
  if (clab_ebshr(yv.p(), pw, n, nv, bb, fit, sh.data(), fa.data()) != CLAB_OK) return c.fail("`w` must have a positive sum.");
  static const char *const fn[] = {"mean", "total", "noise", "signal", "bstar", "b"};
  Frame est;
  est.num("y", yv.v);
  est.num("shrunk", sh);
  est.num("fall", fa);
  out = vobj();
  put(out, "fit", nvec(6, fn, fit));
  put(out, "est", est.o);
  return true;
}

bool kscan(Ctx &c, Val &out) {
  Pts p;
  CK(xy(c, PROJ, p));
  const int n = p.n;
  if (n < 3) return c.fail("at least 3 incidents are needed.");
  Dv tv, rv, wv;
  CK(num(c, c["t"], "t", tv, n));
  CK(edges(c, "r", rv));
  CK(num(c, c["tw"], "tw", wv, -1, 0.0, true));
  int mc = 0, tp = 0, ns = 0;
  CK(wholed(c, "minc", 2, 1, mc));
  CK(wholed(c, "top", 1, 1, tp));
  CK(wholed(c, "nsim", 999, 0, ns));
  unsigned int sd = 0;
  CK(nseed(c, ns, sd));
  const int nr = static_cast<int>(rv.n()), nw = static_cast<int>(wv.n());
  const int tq = tp < n ? tp : n;
  std::vector<double> o(10 * zz(tq));
  int nc = 0;
  if (clab_kscan(p.x.data(), p.y.data(), tv.p(), n, rv.p(), nr, wv.p(), nw, mc, ns, sd, tq, o.data(), &nc) != CLAB_OK) {
    return c.fail("the kernel refused the input.");
  }
  std::vector<double> cc[12];
  for (auto &v : cc) v.resize(zz(nc));
  std::vector<int> mb(zz(n), INT_MIN);
  for (int k = nc - 1; k >= 0; --k) {
    const double *v = o.data() + 10 * zz(k);
    const std::size_t ci = static_cast<std::size_t>(v[0]), sk = zz(k);
    cc[0][sk] = p.x[ci];
    cc[1][sk] = p.y[ci];
    for (int j = 0; j < 7; ++j) cc[2 + j][sk] = v[1 + j];
    cc[9][sk] = (v[4] / v[5]) / ((n - v[4]) / (n - v[5]));
    cc[10][sk] = v[8];
    cc[11][sk] = v[9];
    for (int i = 0; i < n; ++i) {
      const double dx = p.x[zz(i)] - p.x[ci], dy = p.y[zz(i)] - p.y[ci];
      if (dx * dx + dy * dy <= v[1] * v[1] && tv[zz(i)] >= v[2] && tv[zz(i)] < v[3]) mb[zz(i)] = k + 1;
    }
  }
  back(p, cc[0].data(), cc[1].data(), nc);
  static const char *const cn[] = {"x", "y", "r", "start", "end", "n", "exp", "ncirc", "nwin", "rr", "llr", "p"};
  Frame cl;
  for (int k = 0; k < 12; ++k) cl.num(cn[k], cc[k]);
  out = vobj();
  put(out, "cl", cl.o);
  put(out, "memb", unbox(vints(mb.data(), mb.size())));
  put(out, "n", vnum(n));
  put(out, "nsim", vnum(ns));
  return true;
}

bool gmoran(Ctx &c, Val &out) {
  const Val &x = c["x"];
  const bool mat = ismat2(x);
  int n = 0, k = 1;
  std::vector<double> xv;
  std::vector<std::string> cn;
  if (mat) {
    CK(dfmat(c, x, "x", -HUGE_VAL, n, k, xv, cn, "x"));
  } else {
    Dv d;
    CK(num(c, x, "x", d));
    xv = d.v;
    n = static_cast<int>(d.n());
    cn = {"x"};
  }
  if (n < 4) return c.fail("at least 4 units are needed.");
  Wcsr w;
  CK(wreadx(c, c["w"], n, w));
  int ns = 0;
  CK(wholed(c, "nsim", 0, 0, ns));
  unsigned int sd = 0;
  CK(nseed(c, ns, sd));
  std::vector<double> o(10 * zz(k));
  if (clab_gmoran(xv.data(), n, k, w.rp.data(), w.ci.data(), w.w.data(), ns, sd, o.data()) != CLAB_OK) {
    return c.fail("`w` has no positive weight: Moran's I is undefined.");
  }
  static const char *const nm[] = {"i", "e", "vn", "vr", "zn", "zr", "pn", "pr", "pperm", "s0"};
  Frame f;
  f.str("var", cn);
  for (int j = 0; j < 10; ++j) {
    std::vector<double> v(zz(k));
    for (int q = 0; q < k; ++q) v[zz(q)] = o[10 * zz(q) + zz(j)];
    f.num(nm[j], v);
  }
  out = f.o;
  return true;
}

bool bvmor(Ctx &c, Val &out) {
  Dv xv, yv;
  CK(num(c, c["x"], "x", xv));
  const int n = static_cast<int>(xv.n());
  if (n < 3) return c.fail("at least 3 units are needed.");
  CK(num(c, c["y"], "y", yv, n));
  Wcsr w;
  CK(wreadx(c, c["w"], n, w));
  int ns = 0;
  CK(wholed(c, "nsim", 999, 0, ns));
  unsigned int sd = 0;
  CK(nseed(c, ns, sd));
  double o[6] = {0};
  if (clab_bvmor(xv.p(), yv.p(), n, w.rp.data(), w.ci.data(), w.w.data(), ns, sd, o) != CLAB_OK) {
    return c.fail("`x` or `y` is constant, or `w` has no positive weight: the statistic is undefined.");
  }
  static const char *const nm[] = {"i", "mean", "sd", "z", "p", "ptwo"};
  out = nvec(6, nm, o);
  return true;
}

bool sarml(Ctx &c, Val &out) {
  Dv yv;
  CK(num(c, c["y"], "y", yv));
  const int n = static_cast<int>(yv.n());
  const std::size_t sn = yv.n();
  if (n < 4) return c.fail("at least 4 units are needed.");
  if (n > 5000) return c.fail("at most 5000 units: the likelihood takes a dense determinant.");
  const Val &x = c["x"];
  int nr = 0, px = 0;
  std::vector<double> xv;
  std::vector<std::string> xn;
  if (ismat2(x)) {
    if (isdf(x) ? static_cast<int>(vlen(x.a[0])) != n : static_cast<int>(x.a.size()) != n) {
      return c.fail("`x` must have one row per value of `y`.");
    }
    CK(dfmat(c, x, "x", -HUGE_VAL, nr, px, xv, xn, "x"));
  } else {
    if (static_cast<int>(vlen(x)) != n) return c.fail("`x` must have one row per value of `y`.");
    if (!isnum(x)) return c.fail("`x` must be numeric.");
    xv = dbl(x);
    for (double e : xv) {
      if (!std::isfinite(e)) return c.fail("`x` must hold finite numbers with no missing values.");
    }
    px = 1;
    xn = {"x"};
  }
  Wcsr w;
  CK(wreadx(c, c["w"], n, w));
  static const char *const mo[] = {"lag", "error", "sac"};
  static const char *const me[] = {"ml", "2sls"};
  int md = 0, mt = 0;
  CK(pickd(c, "model", mo, 3, 0, md));
  CK(pickd(c, "method", me, 2, 0, mt));
  bool db = false, ic = true;
  CK(flagd(c, "durbin", false, db));
  CK(flagd(c, "int", true, ic));
  if (mt == 1 && md != 0) return c.fail("two-stage least squares fits the lag model only.");
  std::vector<int> var;
  for (int j = 0; j < px; ++j) {
    const double *v = xv.data() + zz(j) * sn;
    bool same = true;
    for (int i = 1; i < n && same; ++i) same = v[i] == v[0];
    if (!same) var.push_back(j);
  }
  const int nv = static_cast<int>(var.size());
  auto lag = [&](const double *v, double *o) {
    for (int i = 0; i < n; ++i) {
      double s = 0;
      for (int q = w.rp[zz(i)]; q < w.rp[zz(i) + 1]; ++q) s += w.w[zz(q)] * v[w.ci[zz(q)]];
      o[i] = s;
    }
  };
  const int p = (ic ? 1 : 0) + px + (db ? nv : 0);
  if (p >= n) return c.fail("`x` has too many columns for the number of units.");
  std::vector<double> X(sn * zz(p));
  std::vector<std::string> tn;
  int cc = 0;
  if (ic) {
    std::fill(X.begin(), X.begin() + n, 1.0);
    tn.push_back("(Intercept)");
    ++cc;
  }
  for (int j = 0; j < px; ++j) {
    std::copy(xv.begin() + static_cast<std::ptrdiff_t>(zz(j) * sn), xv.begin() + static_cast<std::ptrdiff_t>(zz(j + 1) * sn),
              X.begin() + static_cast<std::ptrdiff_t>(zz(cc) * sn));
    tn.push_back(xn[zz(j)]);
    ++cc;
  }
  if (db) {
    for (int j : var) {
      lag(xv.data() + zz(j) * sn, X.data() + zz(cc) * sn);
      tn.push_back("lag." + xn[zz(j)]);
      ++cc;
    }
  }
  tn.push_back("rho");
  tn.push_back("lambda");
  const int ni = mt == 1 ? (db ? nv : 2 * nv) : 0;
  std::vector<double> H(sn * zz(ni > 0 ? ni : 1));
  if (mt == 1) {
    if (nv == 0) return c.fail("two-stage least squares needs a column of `x` that varies.");
    std::vector<double> t(sn);
    int h = 0;
    for (int j : var) {
      lag(xv.data() + zz(j) * sn, t.data());
      if (!db) std::copy(t.begin(), t.end(), H.begin() + static_cast<std::ptrdiff_t>(zz(h++) * sn));
      lag(t.data(), H.data() + zz(h++) * sn);
    }
  }
  const int k = p + 2;
  std::vector<double> cf(4 * zz(p) + 4), vc(zz(k) * zz(k));
  double st[5] = {0};
  if (clab_sarml(yv.p(), X.data(), n, p, w.rp.data(), w.ci.data(), w.w.data(), md, mt, mt == 1 ? H.data() : nullptr, ni,
                 cf.data(), vc.data(), st) != CLAB_OK) {
    return c.fail("the model cannot be fitted: `x` is collinear (after filtering), or `w` has no positive weight.");
  }
  std::vector<std::string> term;
  std::vector<double> c0[6];
  for (int j = 0; j < k; ++j) {
    if ((j == p && md == 1) || (j == p + 1 && md == 0)) continue;
    const std::size_t sj = zz(j), sk = zz(k);
    term.push_back(tn[sj]);
    c0[0].push_back(cf[sj]);
    c0[1].push_back(cf[sk + sj]);
    c0[2].push_back(cf[sj] / cf[sk + sj]);
    c0[3].push_back(2 * clab_pnorm(-std::fabs(c0[2].back())));
    c0[4].push_back(j < p ? cf[2 * zz(p) + 4 + sj] : NaN);
    c0[5].push_back(j < p ? cf[3 * zz(p) + 4 + sj] : NaN);
  }
  static const char *const cn[] = {"est", "se", "z", "p", "ols", "olsse"};
  Frame tab;
  tab.str("term", term);
  for (int j = 0; j < 6; ++j) tab.num(cn[j], c0[j]);
  static const char *const fn[] = {"loglik", "aic", "s2", "llols", "n"};
  const double fv[5] = {st[0], st[2], st[1], st[3], static_cast<double>(n)};
  out = vobj();
  put(out, "coef", tab.o);
  put(out, "fit", nvec(5, fn, fv));
  put(out, "vcov", vmat(vc, k, k));
  put(out, "model", vstr(mo[md]));
  put(out, "method", vstr(me[mt]));
  put(out, "durbin", vbool(db));
  put(out, "nlag", vnum(db ? nv : 0));
  put(out, "w", c["w"]);
  return true;
}

// spill(): from a sarml() fit (its JSON), or from w, rho, b and th.
bool spill(Ctx &c, Val &out) {
  Val w, rho, b, th, vc;
  std::vector<std::string> nm;
  if (c.has("fit")) {
    const Val &f = c["fit"];
    const Val &md = col(f, "model");
    if (f.t != Val::OBJ || md.t != Val::STR || (md.s != "lag" && md.s != "sac")) {
      return failas(c, "spill: `fit` must be a sarml() fit of a lag or combined model; in an error model the "
                       "coefficients are the effects.");
    }
    const Val &cf = col(f, "coef"), &v = col(f, "vcov");
    std::vector<std::string> tm;
    if (!alabs(col(cf, "term"), tm)) return c.fail("`fit` is not a sarml() fit.");
    const std::vector<double> est = dbl(col(cf, "est"));
    int dim[2] = {-1, -1};
    if (est.size() != tm.size() || v.t != Val::ARR || !shape(v, 2, dim, 0) || dim[0] != dim[1] || dim[0] < 2) {
      return c.fail("`fit` is not a sarml() fit.");
    }
    const int pv = dim[0] - 2;
    if (pv < 0 || static_cast<std::size_t>(pv) > tm.size()) return c.fail("`fit` is not a sarml() fit.");
    // the names of the vcov rows: the terms, then rho and lambda
    std::vector<std::string> rn(tm.begin(), tm.begin() + pv);
    rn.push_back("rho");
    rn.push_back("lambda");
    const int nl = static_cast<int>(val1(col(f, "nlag")));
    std::vector<std::string> lg, xs;
    if (nl > 0) lg.assign(rn.begin() + (pv - nl), rn.begin() + pv);
    for (int j = 0; j < pv; ++j) {
      const std::string &s = rn[zz(j)];
      if (s != "(Intercept)" && std::find(lg.begin(), lg.end(), s) == lg.end() &&
          std::find(xs.begin(), xs.end(), s) == xs.end()) {
        xs.push_back(s);
      }
    }
    if (nl > 0) {
      xs.clear();
      for (const std::string &s : lg) xs.push_back(s.compare(0, 4, "lag.") == 0 ? s.substr(4) : s);
    }
    auto at = [&](const std::string &s) {
      for (std::size_t q = 0; q < tm.size(); ++q) {
        if (tm[q] == s) return est[q];
      }
      return NA;
    };
    std::vector<std::string> id = {"rho"};
    id.insert(id.end(), xs.begin(), xs.end());
    id.insert(id.end(), lg.begin(), lg.end());
    std::vector<int> ix;
    for (const std::string &s : id) {
      const auto it = std::find(rn.begin(), rn.end(), s);
      if (it == rn.end()) return c.fail("`fit` is not a sarml() fit.");
      ix.push_back(static_cast<int>(it - rn.begin()));
    }
    rho = vnum(at("rho"));
    b = varr();
    for (const std::string &s : xs) b.a.push_back(vnum(at(s)));
    if (nl > 0) {
      th = varr();
      for (const std::string &s : lg) th.a.push_back(vnum(at(s)));
    }
    vc = varr();
    for (int r : ix) {
      Val row = varr();
      for (int q : ix) row.a.push_back(v.a[zz(r)].a[zz(q)]);
      vc.a.push_back(row);
    }
    w = col(f, "w");
    nm = xs;
  } else {
    if (!c.has("w") || !c.has("rho") || !c.has("b")) return failas(c, "spill: give a sarml() fit, or `w`, `rho` and `b`.");
    w = c["w"];
    rho = c["rho"];
    b = c["b"];
    th = c["th"];
  }
  Dv bv, tv;
  CK(num(c, b, "b", bv));
  const int k = static_cast<int>(bv.n());
  if (nm.empty()) {
    nm = bv.nm;
    if (nm.empty()) nm = nmsof({}, zz(k));
    if (bv.nm.empty()) {
      for (int j = 0; j < k; ++j) nm[zz(j)] = "x" + std::to_string(j + 1);
    }
  }
  if (th.t != Val::NUL) CK(num(c, th, "th", tv, k));
  double r = 0.0;
  CK(one(c, rho, "rho", -HUGE_VAL, HUGE_VAL, r));
  int n = 0;
  if (w.t == Val::OBJ && col(w, "neighbours").t != Val::NUL) {
    n = static_cast<int>(vlen(col(w, "neighbours")));
  } else if (w.t == Val::ARR && !w.a.empty() && w.a[0].t == Val::ARR) {
    n = static_cast<int>(w.a.size());
  } else if (w.t == Val::OBJ && col(w, "n").t != Val::NUL) {
    n = static_cast<int>(val1(col(w, "n")));
  } else {
    return c.fail("`w` must be the list spwt() returns, an spdep listw or a matrix.");
  }
  if (n < 2) return c.fail("`w` must be for at least 2 units.");
  if (n > 5000) return c.fail("at most 5000 units: the effects take a dense inverse.");
  Wcsr wc;
  CK(wreadx(c, w, n, wc));
  std::vector<double> vv;
  if (vc.t != Val::NUL) {
    int a = 0, d = 0;
    CK(mnum(c, vc, "vcov", -HUGE_VAL, a, d, vv));
  }
  std::vector<double> o(6 * zz(k));
  if (clab_spill(n, wc.rp.data(), wc.ci.data(), wc.w.data(), r, bv.p(), tv.n() ? tv.p() : nullptr, k,
                 vv.empty() ? nullptr : vv.data(), o.data()) != CLAB_OK) {
    return c.fail("I - rho W is singular: no spillover multiplier.");
  }
  static const char *const cn[] = {"direct", "indirect", "total", "sed", "sei", "set"};
  Frame f;
  f.str("term", nm);
  for (int j = 0; j < 6; ++j) {
    std::vector<double> v(zz(k));
    for (int q = 0; q < k; ++q) v[zz(q)] = o[6 * zz(q) + zz(j)];
    f.num(cn[j], v);
  }
  out = f.o;
  return true;
}

bool crossk(Ctx &c, Val &out) {
  Pts p;
  CK(xy(c, PLANAR, p));
  const int n = p.n;
  if (n < 2) return c.fail("at least 2 points are needed.");
  Lab g;
  CK(codes(c, c["marks"], "marks", zz(n), g));
  int ti = 0, tj = 0;
  CK(rref(c, c["i"], g.lev, "i", "type", ti));
  CK(rref(c, c["j"], g.lev, "j", "type", tj));
  double wn[4] = {0};
  CK(win(c, c["win"], wn));
  for (int a = 0; a < n; ++a) {
    if (!inside(wn, p.x[zz(a)], p.y[zz(a)])) return c.fail("every point must lie inside `win`.");
  }
  Dv rv;
  CK(num(c, c["r"], "r", rv, -1, 0.0));
  for (std::size_t q = 1; q < rv.n(); ++q) {
    if (!(rv[q] > rv[q - 1])) return c.fail("`r` must increase.");
  }
  const int nr = static_cast<int>(rv.n());
  const bool sg = c.has("sigma");
  double sd = 0.0;
  if (sg) CK(one(c, c["sigma"], "sigma", 0.0, HUGE_VAL, sd));
  int ns = 0;
  CK(wholed(c, "nsim", 99, 0, ns));
  unsigned int se = 0;
  if (sg && ns > 0) CK(seed(c, se));
  std::vector<double> kv(zz(nr)), l(zz(nr)), th(zz(nr));
  double seg[2] = {0};
  if (clab_crossk(p.x.data(), p.y.data(), g.code.data(), n, ti - 1, tj - 1, wn, rv.p(), nr, sd, ns, se, kv.data(), seg) !=
      CLAB_OK) {
    return c.fail("type `i` needs a point and type `j` another, or the kernel weights of a point vanish: widen "
                  "`sigma`.");
  }
  for (std::size_t q = 0; q < zz(nr); ++q) {
    l[q] = std::sqrt(kv[q] / 3.14159265358979323846);
    th[q] = 3.14159265358979323846 * rv[q] * rv[q];
  }
  Frame f;
  f.num("r", rv.v);
  f.num("k", kv);
  f.num("l", l);
  f.num("theo", th);
  Val sv;
  if (sg) {
    static const char *const sn[] = {"t", "p", "nsim"};
    const double v[3] = {seg[0], seg[1], static_cast<double>(ns)};
    sv = nvec(3, sn, v);
  }
  out = vobj();
  put(out, "k", f.o);
  put(out, "seg", sv);
  put(out, "i", vstr(g.lev[zz(ti - 1)]));
  put(out, "j", vstr(g.lev[zz(tj - 1)]));
  return true;
}

bool qdrtst(Ctx &c, Val &out) {
  Pts p;
  CK(xy(c, PLANAR, p));
  const int n = p.n;
  double wn[4] = {0};
  CK(win(c, c["win"], wn));
  for (int a = 0; a < n; ++a) {
    if (!inside(wn, p.x[zz(a)], p.y[zz(a)])) return c.fail("every point must lie inside `win`.");
  }
  int ax = 0, ay = 0;
  CK(whole(c, c["nx"], "nx", 1, ax));
  Val d1;
  CK(whole(c, arg(c, "ny", d1, ax), "ny", 1, ay));
  if (static_cast<double>(ax) * ay < 2 || static_cast<double>(ax) * ay > 1e7) {
    return c.fail("`nx` times `ny` must be from 2 to 1e7.");
  }
  std::vector<double> cn(zz(ay) * zz(ax));
  double st[3] = {0};
  if (clab_qdrtst(p.x.data(), p.y.data(), n, wn, ax, ay, cn.data(), st) != CLAB_OK) {
    return c.fail("the kernel refused the input.");
  }
  const double up = clab_pchisq(st[0], st[1], 0);
  static const char *const sn[] = {"x2", "df", "p", "ptwo", "lambda"};
  const double v[5] = {st[0], st[1], up, std::fmin(1.0, 2 * std::fmin(up, 1 - up)), st[2]};
  out = vobj();
  put(out, "counts", vmat(cn, ay, ax));
  put(out, "stat", nvec(5, sn, v));
  return true;
}

bool morsta(Ctx &c, Val &out) {
  const Val &x = c["counts"];
  std::vector<double> cv;
  if (ismat2(x)) {
    int a = 0, b = 0;
    std::vector<std::string> nm;
    CK(dfmat(c, x, "counts", 0.0, a, b, cv, nm, "x"));
  } else {
    if (!isnum(x)) return c.fail("`counts` must be numeric.");
    cv = dbl(x);
    for (double e : cv) {
      if (!std::isfinite(e) || e < 0) return c.fail("`counts` must hold finite numbers from 0 with no missing values.");
    }
  }
  const int n = static_cast<int>(cv.size());
  if (n < 2) return c.fail("at least 2 quadrats are needed.");
  double a = 0.0;
  CK(oned(c, "crit", 0.05, 0.0, 1.0, a));
  double o[5] = {0};
  if (clab_morsta(cv.data(), n, clab_qchisq(1 - a / 2, n - 1, 0), clab_qchisq(a / 2, n - 1, 0), o) != CLAB_OK) {
    return c.fail("the counts must add up to more than 1.");
  }
  static const char *const nm[] = {"imor", "mclu", "muni", "imst", "chi", "p"};
  const double v[6] = {o[0], o[1], o[2], o[3], o[4], clab_pchisq(o[4], n - 1, 0)};
  out = nvec(6, nm, v);
  return true;
}

bool dissix(Ctx &c, Val &out) {
  const Val &x = c["x"];
  if (!ismat2(x)) return c.fail("`x` must be a matrix or data frame of counts, units by groups.");
  int n = 0, g = 0;
  std::vector<double> xv;
  std::vector<std::string> gn;
  {
    // the dimensions first, then the values (the order of R's checks)
    int dim[2] = {-1, -1};
    g = isdf(x) ? static_cast<int>(x.a.size()) : (shape(x, 2, dim, 0) ? dim[1] : 0);
  }
  if (g < 2) return c.fail("`x` must have at least 2 groups (columns).");
  CK(dfmat(c, x, "x", 0.0, n, g, xv, gn, "g"));
  double dl = 0.0;
  CK(oned(c, "delta", 0.5, 0.0, 1.0, dl));
  const bool hw = c.has("w");
  Wcsr w;
  if (hw) CK(wreadx(c, c["w"], n, w));
  std::vector<double> d(zz(g) * zz(g)), dm(hw ? zz(g) * zz(g) : 0), gr(4 * zz(g));
  if (clab_dissix(xv.data(), n, g, dl, hw ? w.rp.data() : nullptr, hw ? w.ci.data() : nullptr, hw ? w.w.data() : nullptr,
                  d.data(), gr.data(), hw ? dm.data() : nullptr) != CLAB_OK) {
    return c.fail("every group needs a positive total, and `w` a positive weight.");
  }
  static const char *const cn[] = {"is", "h", "atk", "gini"};
  Frame gf;
  gf.str("group", gn);
  for (int j = 0; j < 4; ++j) {
    std::vector<double> v(zz(g));
    for (int k = 0; k < g; ++k) v[zz(k)] = gr[4 * zz(k) + zz(j)];
    gf.num(cn[j], v);
  }
  out = vobj();
  put(out, "d", vmat(d, g, g));
  put(out, "grp", gf.o);
  put(out, "dm", hw ? vmat(dm, g, g) : Val());
  put(out, "delta", vnum(dl));
  return true;
}

bool pinpl(Ctx &c, Val &out) {
  Pts p;
  CK(xy(c, PLANAR, p));
  const int n = p.n;
  const Val &poly = c["poly"];
  const char *pm = "`poly` must be a data frame of x, y and id (and ring).";
  if (poly.t != Val::OBJ) return c.fail("%s", pm);
  Val px = col(poly, "x"), py = col(poly, "y");
  if (px.t == Val::NUL) px = col(poly, "lon");
  if (py.t == Val::NUL) py = col(poly, "lat");
  const Val &pi = col(poly, "id"), &pr = col(poly, "ring");
  if (px.t == Val::NUL || py.t == Val::NUL || pi.t == Val::NUL) return c.fail("%s", pm);
  Dv vx, vy;
  CK(num(c, px, "poly$x", vx));
  const int nv = static_cast<int>(vx.n());
  if (nv < 3) return c.fail("`poly` needs at least 3 vertices.");
  CK(num(c, py, "poly$y", vy, nv));
  Lab id, rr;
  CK(codes(c, pi, "poly$id", zz(nv), id));
  if (pr.t == Val::NUL) rr.code.assign(zz(nv), 0);
  else CK(codes(c, pr, "poly$ring", zz(nv), rr));
  std::vector<int> o(zz(nv));
  for (int v = 0; v < nv; ++v) o[zz(v)] = v;
  const int *ic = id.code.data(), *rc = rr.code.data();
  std::stable_sort(o.begin(), o.end(), [&](int a, int b) { return ic[a] < ic[b] || (ic[a] == ic[b] && rc[a] < rc[b]); });
  std::vector<double> sx(zz(nv)), sy(zz(nv));
  std::vector<int> sp(zz(nv)), sr(zz(nv));
  for (std::size_t v = 0; v < zz(nv); ++v) {
    sx[v] = vx[zz(o[v])];
    sy[v] = vy[zz(o[v])];
    sp[v] = ic[o[v]];
    sr[v] = rc[o[v]];
  }
  const int np = id.g;
  Dv pv;
  const bool hp = c.has("pop");
  if (hp) CK(num(c, c["pop"], "pop", pv, np, 0.0));
  std::vector<int> reg(zz(n)), nreg(zz(n));
  if (clab_pinpl(p.x.data(), p.y.data(), n, sx.data(), sy.data(), sp.data(), sr.data(), nv, np, reg.data(), nreg.data()) !=
      CLAB_OK) {
    return c.fail("the kernel refused the input.");
  }
  Val lab = varr();
  std::vector<double> units(zz(np), 0.0);
  for (int g : reg) {
    lab.a.push_back(g >= 0 ? vstr(id.lev[zz(g)]) : Val());
    if (g >= 0) units[zz(g)] += 1;
  }
  double tot = 0, cov = 0, with = 0, nu = 0;
  for (int q = 0; q < np; ++q) {
    const double pp = hp ? pv[zz(q)] : NaN;
    tot += pp;
    if (units[zz(q)] > 0) {
      cov += pp;
      with += 1;
    }
    nu += units[zz(q)];
  }
  std::vector<double> cpop(zz(np)), cs(zz(np));
  std::vector<int> ch(zz(np));
  for (std::size_t q = 0; q < zz(np); ++q) {
    cpop[q] = hp ? pv[q] : NaN;
    ch[q] = units[q] > 0;
    cs[q] = tot > 0 ? cpop[q] / tot : NaN;
  }
  Frame rf;
  rf.str("region", id.lev);
  rf.num("units", units);
  rf.num("pop", cpop);
  rf.lgl("hasunit", ch.data(), ch.size());
  rf.num("popshr", cs);
  static const char *const cn[] = {"regions", "withunit", "units", "totpop", "covpop", "uncpop", "covshr"};
  const double cv[7] = {static_cast<double>(np), with, nu, tot, cov, tot - cov, tot > 0 ? cov / tot : NaN};
  out = vobj();
  put(out, "reg", unbox(lab));
  put(out, "nreg", unbox(vints(nreg.data(), nreg.size())));
  put(out, "regs", rf.o);
  put(out, "cov", nvec(7, cn, cv));
  return true;
}

bool arwt(Ctx &c, Val &out) {
  Dv xv;
  CK(num(c, c["x"], "x", xv));
  const int ns = static_cast<int>(xv.n());
  Val idv;
  if (c.has("id")) {
    idv = c["id"];
  } else if (!xv.nm.empty()) {
    idv = vstrs(xv.nm);
  } else {
    return failas(c, "arwt: `id` must name each source unit (or `x` have names).");
  }
  Lab src;
  CK(codes(c, idv, "id", zz(ns), src));
  if (src.g != ns) return c.fail("`id` must name each source unit once.");
  const int np = static_cast<int>(vlen(c["from"]));
  if (np < 1) return c.fail("`from` must hold at least one pair.");
  Lab tg, fr;
  CK(codes(c, c["to"], "to", zz(np), tg));
  CK(codes(c, c["from"], "from", zz(np), fr));
  std::vector<int> fi(zz(np));
  for (std::size_t q = 0; q < zz(np); ++q) {
    const std::string &lb = fr.lev[zz(fr.code[q])];
    int k = -1;
    for (int s = 0; s < src.g && k < 0; ++s) {
      if (src.lev[zz(s)] == lb) k = s;
    }
    if (k < 0) return c.fail("`from` names a unit not in `id`: %s.", lb.c_str());
    fi[q] = k;
  }
  std::vector<double> val(zz(ns));
  for (std::size_t s = 0; s < zz(ns); ++s) val[zz(src.code[s])] = xv[s];
  std::vector<double> ar(zz(np), 1.0);
  if (c.has("area")) {
    Dv av;
    CK(num(c, c["area"], "area", av, np, 0.0));
    ar = av.v;
  }
  Dv tv;
  {
    Val dt = varr();
    dt.a = {vnum(10), vnum(20)};
    CK(num(c, c.has("top") ? c["top"] : dt, "top", tv, -1, 1.0));
  }
  const int nt = tg.g;
  std::vector<double> ov(zz(nt)), sh(zz(np));
  std::vector<int> mn(zz(ns));
  double lost = 0;
  if (clab_arwt(val.data(), ns, fi.data(), tg.code.data(), ar.data(), np, nt, ov.data(), sh.data(), mn.data(), &lost) !=
      CLAB_OK) {
    return c.fail("the kernel refused the input.");
  }
  Frame of;
  of.str("to", tg.lev);
  of.num("value", ov);
  std::vector<double> best(zz(ns), -1.0);
  for (std::size_t q = 0; q < zz(np); ++q) best[zz(fi[q])] = std::max(best[zz(fi[q])], sh[q]);
  bool exact = true;
  std::vector<int> cnt(zz(ns), 0);
  for (int f : fi) ++cnt[zz(f)];
  Val mt = varr();
  std::vector<double> ms(zz(ns));
  for (std::size_t s = 0; s < zz(ns); ++s) {
    mt.a.push_back(mn[s] >= 0 ? vstr(tg.lev[zz(mn[s])]) : Val());
    ms[s] = mn[s] >= 0 ? best[s] : NaN;
    exact = exact && cnt[s] <= 1;
  }
  Frame mf;
  mf.str("from", src.lev);
  mf.add("to", mt);
  mf.num("share", ms);
  std::vector<double> v(ov);
  std::sort(v.begin(), v.end());
  double sum = 0, cum = 0, scum = 0;
  for (double e : v) sum += e;
  for (double e : v) {
    cum += e;
    scum += cum;
  }
  std::vector<std::string> cnn = {"gini"};
  std::vector<double> cc = {sum > 0 ? (nt + 1 - 2 * scum / sum) / nt : 0.0};
  for (double t : tv.v) {
    const int m = static_cast<int>(std::min<double>(nt, std::floor(t)));
    double s = 0;
    for (int q = nt - m; q < nt; ++q) s += v[zz(q)];
    cc.push_back(sum > 0 ? s / sum : NaN);
    char buf[48];
    std::snprintf(buf, sizeof buf, "top%g", std::floor(t));
    cnn.push_back(buf);
  }
  out = vobj();
  put(out, "out", of.o);
  put(out, "map", mf.o);
  put(out, "exact", vbool(exact));
  put(out, "lost", vnum(lost));
  put(out, "conc", nvec(cnn, cc.data()));
  return true;
}

bool bufas(Ctx &c, Val &out) {
  Pts p;
  CK(xy(c, PROJ, p));
  const int n = p.n;
  Pts s, q;
  CK(xy(c, c["site"], Val(), p.ll ? RAW : PLANAR, s));
  CK(xy(c, c["comp"], Val(), p.ll ? RAW : PLANAR, q));
  std::vector<double> sx = s.x, sy = s.y, cx = q.x, cy = q.y;
  CK(fwd(c, p, sx.data(), sy.data(), s.n));
  CK(fwd(c, p, cx.data(), cy.data(), q.n));
  double rr = 0.0;
  CK(one(c, c["r"], "r", 0.0, HUGE_VAL, rr));
  Dv bv;
  Val d1;
  CK(num(c, arg(c, "buf", d1, 3 * rr), "buf", bv, 1, 0.0));
  const double bf = bv[0];
  std::vector<int> gi(zz(n));
  std::vector<double> ds(zz(n)), dc(zz(n));
  if (clab_bufas(p.x.data(), p.y.data(), n, sx.data(), sy.data(), s.n, cx.data(), cy.data(), q.n, rr, bf, gi.data(),
                 ds.data(), dc.data()) != CLAB_OK) {
    return c.fail("the kernel refused the input.");
  }
  static const char *const lv[] = {"treated", "comparison", "excluded"};
  double cnt[3] = {0, 0, 0};
  for (int g : gi) cnt[g] += 1;
  Frame pf;
  pf.add("grp", fac(gi.data(), gi.size(), lv));
  pf.num("dsite", ds);
  pf.num("dcomp", dc);
  const double par[2] = {rr, bf};
  static const char *const qn[] = {"r", "buf"};
  out = vobj();
  put(out, "pts", pf.o);
  put(out, "n", nvec(3, lv, cnt));
  put(out, "par", nvec(2, qn, par));
  return true;
}
