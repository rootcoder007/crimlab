// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Spatial functions (src/bsp.cpp): geography, centrography, distance
// analysis, hot spots, space-time tests and the journey to crime.

// posn(), in(): a positive finite number; one number from lo to hi.
bool posn(Ctx &c, const char *nm, double def, double &v) { return oned(c, nm, def, 0.0, Inf, v); }
bool inr(Ctx &c, const Val &x, const char *nm, double lo, double hi, double &v) {
  v = val1(x);
  if (!(v >= lo && v <= hi)) return c.fail("`%s` must be one number from %g to %g.", nm, lo, hi);
  return true;
}
double wsum(const std::vector<double> &w, int n) {
  if (w.empty()) return n;
  double s = 0;
  for (double e : w) s += e;
  return s;
}
// kdiv(): the divisor offset k in {0, 1, 2}, below the total weight.
bool kdiv(Ctx &c, double W, int &v) {
  CK(wholed(c, "k", 1, 0, v));
  if (v > 2) return c.fail("`k` must be 0, 1 or 2.");
  if (!(W > v)) return c.fail("the total weight must exceed `k`.");
  return true;
}
// A kernel's code: its own words when msg is null.
bool kc(Ctx &c, int rc, const char *msg = nullptr) {
  if (rc == CLAB_OK) return true;
  return c.fail("%s", msg ? msg : clab_strerror(rc));
}
bool metric(Ctx &c, int &me) {
  static const char *const mo[] = {"euclid", "manhattan"};
  return pickd(c, "metric", mo, 2, 0, me);
}

bool hvdist(Ctx &c, Val &out) {
  double rad = 0.0;
  CK(posn(c, "r", 6371008.8, rad));
  Pts a, b;
  CK(xy(c, c["from"], nul, RAW, a));
  CK(xy(c, c["to"], nul, RAW, b));
  if (a.n != b.n && a.n != 1 && b.n != 1) return c.fail("`from` and `to` must have the same number of points, or one.");
  const int n = std::max(a.n, b.n);
  std::vector<double> q(4 * static_cast<std::size_t>(n)), d(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    const std::size_t ia = static_cast<std::size_t>(a.n == 1 ? 0 : i), ib = static_cast<std::size_t>(b.n == 1 ? 0 : i);
    q[static_cast<std::size_t>(i)] = a.x[ia];
    q[static_cast<std::size_t>(n + i)] = a.y[ia];
    q[static_cast<std::size_t>(2 * n + i)] = b.x[ib];
    q[static_cast<std::size_t>(3 * n + i)] = b.y[ib];
  }
  CK(kc(c, clab_hvdist(q.data(), q.data() + n, q.data() + 2 * n, q.data() + 3 * n, n, rad, d.data())));
  out = vauto(d);
  return true;
}

bool eqproj(Ctx &c, Val &out) {
  double rad = 0.0;
  CK(posn(c, "r", 6371008.8, rad));
  Pts p;
  CK(xy(c, RAW, p));
  double c0 = 0.0, c1 = 0.0;
  if (c.has("lon0") != c.has("lat0")) return c.fail("give both `lon0` and `lat0`, or neither.");
  if (!c.has("lon0")) {
    c0 = p.lon0;
    c1 = p.lat0;
  } else {
    CK(inr(c, c["lon0"], "lon0", -180.0, 180.0, c0));
    CK(inr(c, c["lat0"], "lat0", -90.0, 90.0, c1));
  }
  std::vector<double> x(static_cast<std::size_t>(p.n)), y(static_cast<std::size_t>(p.n));
  CK(kc(c, clab_eqproj(p.x.data(), p.y.data(), p.n, c0, c1, rad, x.data(), y.data()),
        "a point lies at the antipode of the centre."));
  Frame f;
  f.num("x", x);
  f.num("y", y);
  out = f.o;
  return true;
}

bool eqinv(Ctx &c, Val &out) {
  double rad = 0.0;
  CK(posn(c, "r", 6371008.8, rad));
  Pts p;
  CK(xy(c, PLANAR, p));
  double c0 = 0.0, c1 = 0.0;
  CK(inr(c, c["lon0"], "lon0", -180.0, 180.0, c0));
  CK(inr(c, c["lat0"], "lat0", -90.0, 90.0, c1));
  std::vector<double> lo(static_cast<std::size_t>(p.n)), la(static_cast<std::size_t>(p.n));
  CK(kc(c, clab_eqinv(p.x.data(), p.y.data(), p.n, c0, c1, rad, lo.data(), la.data()),
        "a point lies farther than 2 r from the centre."));
  Frame f;
  f.num("lon", lo);
  f.num("lat", la);
  out = f.o;
  return true;
}

// x, y and w matched by name, then the points and the weights.
bool xyw(Ctx &c, Mode m, Pts &p, std::vector<double> &w) {
  Args a(c, {"x", "y", "w"});
  CK(a.byname(c));
  CK(xy(c, a[0], a[1], m, p));
  return pwts(c, a[2], static_cast<std::size_t>(p.n), w);
}

bool cgstat(Ctx &c, Val &out) {
  Pts p;
  std::vector<double> w;
  CK(xyw(c, PROJ, p, w));
  int kk = 0;
  CK(kdiv(c, wsum(w, p.n), kk));
  double o[10];
  CK(kc(c, clab_cgstat(p.x.data(), p.y.data(), ptr(w), p.n, kk, o)));
  if (p.ll) {
    back(p, o + 1, o + 2, 1);
    for (int j = 6; j < 10; ++j) o[j] = NaN;
  }
  static const char *const nm[] = {"w", "mx", "my", "sx", "sy", "sdd", "gx", "gy", "hx", "hy"};
  out = nvec(10, nm, o);
  return true;
}

bool cgmedn(Ctx &c, Val &out) {
  Pts p;
  CK(xy(c, PROJ, p));
  std::vector<double> wk(static_cast<std::size_t>(p.n));
  double o[2];
  CK(kc(c, clab_cgmedn(p.x.data(), p.y.data(), p.n, wk.data(), o)));
  back(p, o, o + 1, 1);
  static const char *const nm[] = {"x", "y"};
  out = nvec(2, nm, o);
  return true;
}

bool webpt(Ctx &c, Val &out) {
  Pts p;
  std::vector<double> w;
  CK(xyw(c, PROJ, p, w));
  double s[2];
  const bool has = c.has("start");
  if (has) {
    Dv sv;
    CK(num(c, c["start"], "start", sv, 2));
    s[0] = sv[0];
    s[1] = sv[1];
    CK(fwd(c, p, s, s + 1, 1));
  }
  double tl = 0.0;
  int mi = 0;
  CK(posn(c, "tol", 1e-10, tl));
  CK(wholed(c, "maxit", 1000, 0, mi));
  double o[6];
  CK(kc(c, clab_webpt(p.x.data(), p.y.data(), ptr(w), p.n, has ? s : nullptr, tl, mi, o)));
  back(p, o, o + 1, 1);
  static const char *const nm[] = {"x", "y", "obj", "gap", "iter", "code"};
  out = nvec(6, nm, o);
  return true;
}

bool sdell(Ctx &c, Val &out) {
  Pts p;
  std::vector<double> w;
  CK(xyw(c, PROJ, p, w));
  int kk = 0;
  CK(kdiv(c, wsum(w, p.n), kk));
  double cc = 0.0;
  CK(posn(c, "c", 1, cc));
  double o[6];
  CK(kc(c, clab_sdell(p.x.data(), p.y.data(), ptr(w), p.n, kk, cc, o)));
  back(p, o, o + 1, 1);
  static const char *const nm[] = {"mx", "my", "a", "b", "angle", "area"};
  out = nvec(6, nm, o);
  return true;
}

bool cvhull(Ctx &c, Val &out) {
  Pts p;
  CK(xy(c, PROJ, p));
  std::vector<int> iw(static_cast<std::size_t>(p.n) + 1), h(2 * static_cast<std::size_t>(p.n) + 1);
  int nh = 0;
  double o[2];
  CK(kc(c, clab_cvhull(p.x.data(), p.y.data(), p.n, iw.data(), h.data(), &nh, o),
        "coordinates must be 0 or of magnitude from 1e-140 to 1e150."));
  const std::vector<int> id = plus1(h.data(), nh);
  out = vobj();
  put(out, "idx", unbox(vints(id.data(), id.size())));
  put(out, "area", vnum(o[0]));
  put(out, "perim", vnum(o[1]));
  return true;
}

// Planar points inside the rectangle win, at least two of them.
bool inwin(Ctx &c, Pts &p, double *w) {
  CK(xy(c, PLANAR, p));
  CK(win(c, c["win"], w));
  if (p.n < 2) return c.fail("at least 2 points are needed.");
  for (int i = 0; i < p.n; ++i) {
    if (!inside(w, p.x[static_cast<std::size_t>(i)], p.y[static_cast<std::size_t>(i)])) {
      return c.fail("every point must lie inside `win`.");
    }
  }
  return true;
}
// Radii from 0, strictly increasing.
bool radii(Ctx &c, Dv &v) {
  CK(num(c, c["r"], "r", v, -1, 0.0));
  for (std::size_t j = 1; j < v.n(); ++j) {
    if (!(v[j] > v[j - 1])) return c.fail("`r` must be strictly increasing.");
  }
  return true;
}
const char *const kcor[] = {"none", "border", "modif", "iso"};

bool nnidx(Ctx &c, Val &out) {
  double w[4];
  Pts p;
  CK(inwin(c, p, w));
  const int n = p.n;
  int km = 0, ns = 0, me = 0;
  CK(wholed(c, "kmax", 1, 1, km));
  if (km > n - 1) return c.fail("`kmax` must be below the number of points.");
  CK(wholed(c, "nsim", 0, 0, ns));
  bool dd = false;
  CK(flagd(c, "dist", false, dd));
  CK(metric(c, me));
  unsigned int sd = 0;
  if (ns > 0) CK(seed(c, sd));
  const std::size_t nk = static_cast<std::size_t>(n) * static_cast<std::size_t>(km);
  std::vector<int> iw(2 * static_cast<std::size_t>(n) + 2 + 2 * static_cast<std::size_t>(km) + (me ? nk : 0));
  std::vector<double> wk(3 * static_cast<std::size_t>(km) + 2 * static_cast<std::size_t>(n) + (me ? nk : 0));
  std::vector<double> dv(dd ? nk : 0), o(6 * static_cast<std::size_t>(km) + 4);
  CK(kc(c, (me ? clab_nnidm : clab_nnidx)(p.x.data(), p.y.data(), n, w, km, ns, sd, dd ? dv.data() : nullptr, iw.data(),
                                         wk.data(), o.data())));
  static const char *const on[] = {"mean", "exp", "ratio", "pclu", "preg", "sd"};
  Frame f;
  std::vector<double> kk(static_cast<std::size_t>(km));
  for (int k = 0; k < km; ++k) kk[static_cast<std::size_t>(k)] = k + 1;
  f.num("k", kk);
  for (int j = 0; j < 6; ++j) {
    std::vector<double> col(static_cast<std::size_t>(km));
    for (int k = 0; k < km; ++k) col[static_cast<std::size_t>(k)] = o[static_cast<std::size_t>(6 * k + j)];
    f.num(on[j], col);
  }
  static const char *const cn[] = {"donnelly", "ratio", "z", "zc"};
  out = vobj();
  put(out, "order", f.o);
  put(out, "ce", nvec(4, cn, o.data() + 6 * km));
  if (dd) {
    // the kernel's rows are points; the matrix is n by kmax
    std::vector<double> dm(nk);
    for (int i = 0; i < n; ++i) {
      for (int k = 0; k < km; ++k) {
        dm[static_cast<std::size_t>(i) + static_cast<std::size_t>(n) * k] = dv[static_cast<std::size_t>(i) * km + k];
      }
    }
    put(out, "dist", vmat(dm, n, km));
  } else {
    put(out, "dist", Val());
  }
  return true;
}

bool ripk(Ctx &c, Val &out) {
  double w[4];
  Pts p;
  CK(inwin(c, p, w));
  Dv rv;
  CK(radii(c, rv));
  const int m = static_cast<int>(rv.n());
  int co = 0, me = 0;
  CK(pickd(c, "corr", kcor, 4, 3, co));
  CK(metric(c, me));
  if (me == 1 && co == 3) {
    return c.fail("the isotropic correction has no Manhattan form: use corr = \"border\", \"modif\" or \"none\".");
  }
  std::vector<int> iw(2 * static_cast<std::size_t>(p.n) + 2);
  std::vector<double> wk(4 * rv.n()), k(rv.n()), l(rv.n());
  CK(kc(c, clab_ripkm(p.x.data(), p.y.data(), p.n, w, rv.p(), m, co, me, iw.data(), wk.data(), k.data())));
  for (std::size_t j = 0; j < rv.n(); ++j) l[j] = std::sqrt(k[j] / (me ? 2.0 : 3.14159265358979323846));
  Frame f;
  f.num("r", rv.v);
  f.num("k", k);
  f.num("l", l);
  out = f.o;
  return true;
}

bool ripenv(Ctx &c, Val &out) {
  double w[4];
  Pts p;
  CK(inwin(c, p, w));
  Dv rv;
  CK(radii(c, rv));
  const int m = static_cast<int>(rv.n());
  int co = 0, ns = 0;
  CK(pickd(c, "corr", kcor, 4, 3, co));
  CK(wholed(c, "nsim", 99, 1, ns));
  unsigned int sd;
  CK(seed(c, sd));
  std::vector<int> iw(2 * static_cast<std::size_t>(p.n) + 2);
  std::vector<double> wk(5 * rv.n() + 2 * static_cast<std::size_t>(p.n)), o(4 * rv.n() + 2);
  CK(kc(c, clab_ripenv(p.x.data(), p.y.data(), p.n, w, rv.p(), m, co, ns, sd, iw.data(), wk.data(), o.data())));
  static const char *const en[] = {"obs", "lo", "hi", "mean"};
  Frame f;
  f.num("r", rv.v);
  for (int b = 0; b < 4; ++b) f.num(en[b], o.data() + static_cast<std::size_t>(b) * rv.n(), rv.n());
  static const char *const tn[] = {"t", "p"};
  out = vobj();
  put(out, "env", f.o);
  put(out, "test", nvec(2, tn, o.data() + 4 * rv.n()));
  return true;
}

// zcount(): zone counts, whole numbers from 0, at least one positive.
bool zcount(Ctx &c, int n, std::vector<double> &v) {
  Dv d;
  CK(num(c, c["count"], "count", d, n, 0.0));
  bool any = false;
  for (double e : d.v) {
    if (e != std::floor(e)) return c.fail("`count` must be whole numbers, one per zone.");
    any = any || e > 0;
  }
  if (!any) return c.fail("`count` must have a positive sum.");
  v = d.v;
  return true;
}

// Zones as points: each zone's centroid repeated count times; first[i]
// the copy that stands for zone i, -1 for an empty zone.
struct Zones {
  std::vector<double> x, y;
  std::vector<int> first;
};
bool expand(Ctx &c, const Pts &p, Zones &z) {
  std::vector<double> cv;
  CK(zcount(c, p.n, cv));
  double tot = 0;
  for (double e : cv) tot += e;
  if (tot > 1e8) return c.fail("`count` must sum to at most 1e8.");
  z.first.assign(static_cast<std::size_t>(p.n), -1);
  int t = 0;
  for (std::size_t i = 0; i < static_cast<std::size_t>(p.n); ++i) {
    const int k = static_cast<int>(cv[i]);
    z.first[i] = k > 0 ? t : -1;
    for (int q = 0; q < k; ++q, ++t) {
      z.x.push_back(p.x[i]);
      z.y.push_back(p.y[i]);
    }
  }
  return true;
}
// Per-zone group numbers from the copies' (NA for empty zones), cols
// columns of nc copies each.
std::vector<int> zgroups(const Zones &z, const std::vector<int> &cl, int nz, int cols) {
  std::vector<int> v(static_cast<std::size_t>(nz) * static_cast<std::size_t>(cols));
  const std::size_t nc = z.x.size();
  for (int j = 0; j < cols; ++j) {
    for (int i = 0; i < nz; ++i) {
      const int f = z.first[static_cast<std::size_t>(i)];
      v[static_cast<std::size_t>(i) + static_cast<std::size_t>(nz) * j] =
          f < 0 ? INT_MIN : cl[static_cast<std::size_t>(f) + nc * static_cast<std::size_t>(j)];
    }
  }
  return v;
}
// Group numbers from 0 (or -1) to R's from 1 (or NA).
std::vector<int> grpno(const int *cl, std::size_t n) {
  std::vector<int> v(n);
  for (std::size_t i = 0; i < n; ++i) v[i] = cl[i] >= 0 ? cl[i] + 1 : INT_MIN;
  return v;
}
// A column-major integer matrix as an array of rows (INT_MIN is null).
Val vimat(const std::vector<int> &x, int nr, int nc) {
  Val v = varr();
  for (int i = 0; i < nr; ++i) {
    Val r = varr();
    for (int j = 0; j < nc; ++j) {
      const int e = x[static_cast<std::size_t>(i) + static_cast<std::size_t>(nr) * j];
      r.a.push_back(e == INT_MIN ? Val() : vnum(e));
    }
    v.a.push_back(r);
  }
  return v;
}
// The cluster table of nnhclu and stacel: nc rows of w values, columns nm;
// columns cx, cx + 1 and cx + 2, cx + 3 hold centres.
Val ctab(const Pts &p, const std::vector<double> &cs, int nc, int w, const char *const *nm, int cx) {
  std::vector<std::vector<double>> col(static_cast<std::size_t>(w), std::vector<double>(static_cast<std::size_t>(nc)));
  for (int i = 0; i < nc; ++i) {
    for (int j = 0; j < w; ++j) {
      col[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)] = cs[static_cast<std::size_t>(i) * w + j];
    }
  }
  back(p, col[static_cast<std::size_t>(cx)].data(), col[static_cast<std::size_t>(cx + 1)].data(), nc);
  back(p, col[static_cast<std::size_t>(cx + 2)].data(), col[static_cast<std::size_t>(cx + 3)].data(), nc);
  Frame f;
  for (int j = 0; j < w; ++j) f.num(nm[j], col[static_cast<std::size_t>(j)]);
  return f.o;
}

bool fzmode(Ctx &c, Val &out) {
  Pts p;
  CK(xy(c, PROJ, p));
  double rad = 0.0;
  CK(inr(c, c["r"], "r", 0.0, std::numeric_limits<double>::max(), rad));
  int me = 0;
  CK(metric(c, me));
  const std::size_t n = static_cast<std::size_t>(p.n);
  std::vector<int> iw(2 * n + 2);
  if (!c.has("count") && me == 0) {
    std::vector<int> o(n);
    CK(kc(c, clab_fzmode(p.x.data(), p.y.data(), p.n, rad, iw.data(), o.data())));
    out = unbox(vints(o.data(), n));
    return true;
  }
  std::vector<double> cv;
  if (c.has("count")) CK(zcount(c, p.n, cv));
  const double dm = rad > 0 ? rad : std::numeric_limits<double>::min();
  std::vector<int> rp(n + 1);
  int nnz = 0;
  CK(kc(c, clab_swband(p.x.data(), p.y.data(), p.n, dm, me, 0, rp.data(), nullptr, nullptr, &nnz, iw.data())));
  std::vector<int> ci(static_cast<std::size_t>(nnz) + 1);
  std::vector<double> dv(static_cast<std::size_t>(nnz) + 1);
  CK(kc(c, clab_swband(p.x.data(), p.y.data(), p.n, dm, me, nnz, rp.data(), ci.data(), dv.data(), &nnz, iw.data())));
  std::vector<double> o(n);
  for (std::size_t i = 0; i < n; ++i) {
    double s = cv.empty() ? 1.0 : cv[i];
    for (int q = rp[i]; q < rp[i + 1]; ++q) s += cv.empty() ? 1.0 : cv[static_cast<std::size_t>(ci[static_cast<std::size_t>(q)])];
    o[i] = s;
  }
  out = vauto(o);
  return true;
}

bool kmpart(Ctx &c, Val &out) {
  Pts p;
  CK(xy(c, PROJ, p));
  const bool zon = c.has("count");
  Zones zn;
  if (zon) CK(expand(c, p, zn));
  Pts pz = p;
  if (zon) {
    pz.x = zn.x;
    pz.y = zn.y;
    pz.n = static_cast<int>(zn.x.size());
  }
  const int n = pz.n;
  int kk = 0;
  CK(whole(c, c["k"], "k", 1, kk));
  if (kk > n) return c.fail("`k` must be at most the number of points.");
  std::vector<double> cen(2 * static_cast<std::size_t>(kk));
  int in0 = 1;
  double par[2] = {0, 0};
  const Val &ini = c["init"];
  if (!c.has("init") || ini.t == Val::STR) {
    static const char *const io[] = {"pp", "cells"};
    int k2 = 0;
    CK(pickd(c, "init", io, 2, 0, k2));
    in0 = 1 + k2;
  } else {
    Pts q;
    CK(xy(c, ini, nul, pz.ll ? RAW : PLANAR, q));
    if (q.n != kk) return c.fail("`init` must give `k` starting centres.");
    std::copy(q.x.begin(), q.x.end(), cen.begin());
    std::copy(q.y.begin(), q.y.end(), cen.begin() + kk);
    CK(fwd(c, pz, cen.data(), cen.data() + kk, kk));
    in0 = 0;
  }
  if (in0 == 2) {
    Dv cv;
    CK(num(c, c["cell"], "cell", cv, 2, 0.0));
    par[0] = cv[0];
    par[1] = cv[1];
    if (!(par[0] > 0)) return c.fail("`cell` must give a positive cell side.");
  }
  static const char *const uo[] = {"mean", "weber"};
  int up = 0, mi = 0;
  CK(pickd(c, "upd", uo, 2, 0, up));
  CK(wholed(c, "maxit", 100, 0, mi));
  unsigned int sd = 0;
  if (in0 == 1) CK(seed(c, sd));
  const std::size_t N = static_cast<std::size_t>(n);
  std::vector<int> cl(N), iw(2 * N);
  std::vector<double> wk(2 * N), o(4 + 3 * static_cast<std::size_t>(kk));
  CK(kc(c, clab_kmpart(pz.x.data(), pz.y.data(), n, kk, in0, par, sd, up, mi, cen.data(), cl.data(), iw.data(),
                       wk.data(), o.data()),
        in0 == 1 ? "fewer distinct points than `k`." : "fewer than `k` cells hold points far enough apart."));
  const std::size_t K = static_cast<std::size_t>(kk);
  std::vector<double> cx(cen.begin(), cen.begin() + kk), cy(cen.begin() + kk, cen.end()), cn3[3];
  for (std::size_t j = 0; j < K; ++j) {
    for (int b = 0; b < 3; ++b) cn3[b].push_back(o[4 + 3 * j + static_cast<std::size_t>(b)]);
  }
  back(p, cx.data(), cy.data(), kk);
  Frame f;
  f.num("x", cx);
  f.num("y", cy);
  f.num("n", cn3[0]);
  f.num("ss", cn3[1]);
  f.num("sd", cn3[2]);
  static const char *const fn[] = {"ss", "sd", "iter", "conv"};
  const std::vector<int> gr = grpno(cl.data(), N);
  const std::vector<int> g2 = zon ? zgroups(zn, gr, p.n, 1) : gr;
  out = vobj();
  put(out, "cl", unbox(vints(g2.data(), g2.size())));
  put(out, "cen", f.o);
  put(out, "fit", nvec(4, fn, o.data()));
  return true;
}

bool nnhclu(Ctx &c, Val &out) {
  Pts p;
  CK(xy(c, PROJ, p));
  const bool zon = c.has("count");
  Zones zn;
  if (zon) CK(expand(c, p, zn));
  Pts pz = p;
  if (zon) {
    pz.x = zn.x;
    pz.y = zn.y;
    pz.n = static_cast<int>(zn.x.size());
  }
  const int n = pz.n;
  if (n < 2) return c.fail("at least 2 points are needed.");
  if (c.has("area") == c.has("fixd")) {
    return c.fail("give one of `area` (for the random-expectation threshold) and `fixd` (a fixed threshold).");
  }
  double ar = 0, fd = 0, pv, cc;
  if (c.has("area")) CK(one(c, c["area"], "area", 0.0, Inf, ar));
  if (c.has("fixd")) CK(one(c, c["fixd"], "fixd", 0.0, Inf, fd));
  CK(oned(c, "p", 0.5, 0.0, 1.0, pv));
  const double z = qnorm(pv);
  int mp = 0, mh = 0, mo = 0, kk = 0;
  CK(wholed(c, "minpts", 10, 2, mp));
  Val d1;
  CK(whole(c, c.has("minhi") ? c["minhi"] : (d1 = vnum(mp)), "minhi", 2, mh));
  CK(wholed(c, "maxord", 3, 1, mo));
  if (mo > 20) return c.fail("`maxord` must be at most 20.");
  CK(wholed(c, "k", 1, 0, kk));
  if (kk > 2) return c.fail("`k` must be 0, 1 or 2.");
  CK(posn(c, "c", 1, cc));
  const std::size_t N = static_cast<std::size_t>(n);
  std::vector<int> iw(6 * N + 2), cl(N * static_cast<std::size_t>(mo));
  std::vector<double> wk(4 * N), cs(11 * N);
  int nc = 0;
  CK(kc(c, clab_nnhclu(pz.x.data(), pz.y.data(), n, ar, fd, z, mp, mh, mo, kk, cc, iw.data(), wk.data(), cl.data(), &nc,
                       cs.data())));
  for (int &q : cl) q = q >= 0 ? q + 1 : INT_MIN;
  static const char *const cn[] = {"order", "points", "units", "mx", "my", "wx", "wy", "a", "b", "angle", "thresh"};
  const std::vector<int> g2 = zon ? zgroups(zn, cl, p.n, mo) : cl;
  out = vobj();
  put(out, "cl", vimat(g2, zon ? p.n : n, mo));
  put(out, "clusters", ctab(pz, cs, nc, 11, cn, 3));
  return true;
}

bool stacel(Ctx &c, Val &out) {
  Pts p;
  CK(xy(c, PLANAR, p));
  const bool zon = c.has("count");
  Zones zn;
  if (zon) CK(expand(c, p, zn));
  Pts pz = p;
  if (zon) {
    pz.x = zn.x;
    pz.y = zn.y;
    pz.n = static_cast<int>(zn.x.size());
  }
  double w[4];
  CK(win(c, c["win"], w));
  int gx = 0, gy = 0, tp = 0, mp = 0, kk = 0;
  CK(wholed(c, "nx", 20, 2, gx));
  CK(wholed(c, "ny", 20, 2, gy));
  if (gx > 10000 || gy > 10000) return c.fail("`nx` and `ny` must be at most 10000.");
  double rad = 0.0, cc = 0.0;
  CK(one(c, c["r"], "r", 0.0, Inf, rad));
  CK(wholed(c, "top", 25, 1, tp));
  CK(wholed(c, "minpts", 5, 2, mp));
  CK(wholed(c, "k", 1, 0, kk));
  if (kk > 2) return c.fail("`k` must be 0, 1 or 2.");
  CK(posn(c, "c", 1, cc));
  const std::size_t n = static_cast<std::size_t>(pz.n), nn = static_cast<std::size_t>(gx) * static_cast<std::size_t>(gy);
  std::vector<int> iw(3 * n + 2 + nn + static_cast<std::size_t>(tp) + std::max(n, 2 * static_cast<std::size_t>(tp))),
      cl(n + 1);
  std::vector<double> wk(nn + 2 * n), cs(8 * n + 1);
  int nc = 0;
  CK(kc(c, clab_stacel(pz.x.data(), pz.y.data(), pz.n, w, gx, gy, rad, tp, mp, kk, cc, iw.data(), wk.data(), cl.data(),
                       &nc, cs.data())));
  static const char *const cn[] = {"points", "mx", "my", "wx", "wy", "a", "b", "angle"};
  const std::vector<int> gr = grpno(cl.data(), n);
  const std::vector<int> g2 = zon ? zgroups(zn, gr, p.n, 1) : gr;
  out = vobj();
  put(out, "cl", unbox(vints(g2.data(), g2.size())));
  put(out, "clusters", ctab(pz, cs, nc, 8, cn, 1));
  return true;
}

// x, y and t matched by name, then the points and the times.
bool xyt(Ctx &c, Pts &p, Dv &tv) {
  Args a(c, {"x", "y", "t"});
  CK(a.byname(c));
  CK(xy(c, a[0], a[1], PROJ, p));
  return num(c, a[2], "t", tv, p.n);
}

bool mantst(Ctx &c, Val &out) {
  Pts p;
  Dv tv;
  CK(xyt(c, p, tv));
  if (p.n < 3 || p.n > 46340) return c.fail("the number of events must be from 3 to 46340.");
  Val d1, d2;
  double a = 0.0, b = 0.0;
  CK(inr(c, arg(c, "cs", d1, 0), "cs", 0.0, std::numeric_limits<double>::max(), a));
  CK(inr(c, arg(c, "ct", d2, 0), "ct", 0.0, std::numeric_limits<double>::max(), b));
  int ns = 0, me = 0;
  CK(wholed(c, "nsim", 999, 0, ns));
  CK(metric(c, me));
  unsigned int sd = 0;
  if (ns > 0) CK(seed(c, sd));
  std::vector<int> iw(static_cast<std::size_t>(p.n));
  double o[5];
  CK(kc(c, clab_mantm(p.x.data(), p.y.data(), tv.p(), p.n, me, a, b, ns, sd, iw.data(), o),
        "the distances or the time intervals are all equal: r is undefined."));
  static const char *const nm[] = {"r", "pairs", "pup", "plo", "p2"};
  out = nvec(5, nm, o);
  return true;
}

bool stmavg(Ctx &c, Val &out) {
  Pts p;
  Dv tv;
  CK(xyt(c, p, tv));
  int sp = 0;
  CK(wholed(c, "span", 5, 1, sp));
  if (sp % 2 == 0) return c.fail("`span` must be odd.");
  const std::size_t n = static_cast<std::size_t>(p.n);
  std::vector<int> iw(n);
  std::vector<double> o(4 * n), mx(n), my(n), k(n);
  CK(kc(c, clab_stmavg(p.x.data(), p.y.data(), tv.p(), p.n, sp, iw.data(), o.data())));
  for (std::size_t q = 0; q < n; ++q) {
    const double *r = o.data() + 4 * q;
    const std::size_t i = static_cast<std::size_t>(r[0]);
    mx[i] = r[1];
    my[i] = r[2];
    k[i] = r[3];
  }
  back(p, mx.data(), my.data(), p.n);
  Frame f;
  f.num("mx", mx);
  f.num("my", my);
  f.num("n", k);
  out = f.o;
  return true;
}

// nummat(): finite numbers from 0.
bool nummat(Ctx &c, const std::vector<double> &v, const char *nm) {
  for (double e : v) {
    if (!(e >= 0) || !std::isfinite(e)) return c.fail("`%s` must be finite numbers from 0 with no missing values.", nm);
  }
  return true;
}

// A grid of nx by ny cells over win.
struct Grid {
  int nx = 0, ny = 0;
  double x0 = 0.0, y0 = 0.0, dx = 0.0, dy = 0.0;
};
bool mkgrid(Ctx &c, int nx, int ny, Grid &g) {
  double w[4];
  CK(win(c, c["win"], w));
  g = Grid{nx, ny, w[0], w[2], 0, 0};
  if (static_cast<double>(g.nx) * g.ny > 1e9) return c.fail("`nx` times `ny` must be at most 1e9 cells.");
  g.dx = (w[1] - w[0]) / g.nx;
  g.dy = (w[3] - w[2]) / g.ny;
  return true;
}

// A decay function: its type, parameters and table (type "table").
struct Decay {
  int type = 0;
  Dv par, tab;
};
bool mkdecay(Ctx &c, Decay &d) {
  static const char *const dtype[] = {"linear", "negexp", "normal", "lognormal", "trunc", "table"};
  static const int dlen[] = {2, 2, 3, 3, 3, 1};
  CK(pick(c, c["type"], "type", dtype, 6, d.type));
  CK(num(c, c["par"], "par", d.par, dlen[d.type]));
  if (d.type == 5) {
    CK(num(c, c["tab"], "tab", d.tab, -1, 0.0));
    if (d.tab.n() < 2) return c.fail("`tab` must hold at least 2 values.");
  }
  return true;
}

bool jtdecay(Ctx &c, Val &out) {
  Dv dv;
  CK(num(c, c["d"], "d", dv, -1, 0.0));
  Decay f;
  CK(mkdecay(c, f));
  std::vector<double> o(dv.n());
  CK(kc(c, clab_jtdecay(dv.p(), static_cast<int>(dv.n()), f.type, f.par.p(), f.type == 5 ? f.tab.p() : nullptr,
                        f.type == 5 ? static_cast<int>(f.tab.n()) : 0, o.data()),
        "`par` is outside the domain of this decay type."));
  out = vauto(o);
  return true;
}

bool jtsurf(Ctx &c, Val &out) {
  Pts p;
  CK(xy(c, PLANAR, p));
  int nx = 0, ny = 0;
  CK(whole(c, c["nx"], "nx", 1, nx));
  CK(whole(c, c["ny"], "ny", 1, ny));
  Grid g;
  CK(mkgrid(c, nx, ny, g));
  Decay f;
  CK(mkdecay(c, f));
  int me = 0, md = 0;
  CK(metric(c, me));
  static const char *const so[] = {"sum", "lik"};
  CK(pickd(c, "mode", so, 2, 0, md));
  const std::size_t m = static_cast<std::size_t>(g.nx) * static_cast<std::size_t>(g.ny);
  std::vector<double> raw(m), prob(m);
  double pk[3];
  CK(kc(c, clab_jtsurf(p.x.data(), p.y.data(), p.n, g.nx, g.ny, g.x0, g.y0, g.dx, g.dy, me, md, f.type, f.par.p(),
                       f.type == 5 ? f.tab.p() : nullptr, f.type == 5 ? static_cast<int>(f.tab.n()) : 0, raw.data(),
                       prob.data(), pk),
        "the surface is 0 everywhere: no cell is within reach of the decay function, or `par` is outside its domain."));
  pk[0] += 1;
  static const char *const pn[] = {"cell", "x", "y"};
  out = vobj();
  put(out, "raw", vmat(raw, g.nx, g.ny));
  put(out, "prob", vmat(prob, g.nx, g.ny));
  put(out, "peak", nvec(3, pn, pk));
  return true;
}

bool jtcal(Ctx &c, Val &out) {
  Args a(c, {"d", "w"});
  CK(a.byname(c));
  Dv dv, tv;
  CK(num(c, a[0], "d", dv, -1, 0.0));
  bool rf = false;
  CK(flagd(c, "refl", true, rf));
  CK(num(c, c["t"], "t", tv, -1, rf ? 0.0 : -HUGE_VAL));
  double bw = 0.0;
  CK(one(c, c["h"], "h", 0.0, Inf, bw));
  std::vector<double> w;
  CK(pwts(c, a[1], dv.n(), w));
  static const char *const ko[] = {"normal", "quartic", "triangular", "uniform", "laplace"};
  int kk = 0;
  CK(pickd(c, "kern", ko, 5, 0, kk));
  std::vector<double> o(tv.n());
  CK(kc(c, clab_jtcal(dv.p(), ptr(w), static_cast<int>(dv.n()), kk, bw, rf ? 1 : 0, tv.p(), static_cast<int>(tv.n()),
                      o.data())));
  out = vauto(o);
  return true;
}

bool jtacc(Ctx &c, Val &out) {
  int dim[2] = {-1, -1};
  if (!shape(c["prob"], 2, dim, 0)) {
    return c.fail("`prob` must be a matrix of cell values (nx by ny), as jtsurf() returns.");
  }
  std::vector<double> pv;
  int nr = 0, nc = 0;
  bool bad = false;
  pv.assign(static_cast<std::size_t>(dim[0]) * static_cast<std::size_t>(dim[1]), 0.0);
  fill(c["prob"], 2, dim, 0, 0, 1, pv, bad, true);
  nr = dim[0];
  nc = dim[1];
  if (pv.empty()) return c.fail("`prob` must be finite numbers from 0 with no missing values.");
  CK(nummat(c, pv, "prob"));
  Grid g;
  CK(mkgrid(c, nr, nc, g));
  Dv bv;
  CK(num(c, c["base"], "base", bv, 2));
  const double hx = bv[0], hy = bv[1];
  const double w[4] = {g.x0, g.x0 + g.dx * g.nx, g.y0, g.y0 + g.dy * g.ny};
  if (!inside(w, hx, hy)) return c.fail("`base` must lie inside `win`.");
  double o[6];
  CK(kc(c, clab_jtacc(pv.data(), g.nx, g.ny, g.x0, g.y0, g.dx, g.dy, hx, hy, o)));
  o[4] += 1;
  o[5] += 1;
  static const char *const nm[] = {"p", "hitbest", "hitworst", "err", "peak", "cell"};
  out = nvec(6, nm, o);
  return true;
}

bool bjtc(Ctx &c, Val &out) {
  Dv pv;
  CK(num(c, gvec(c["pj"]), "pj", pv, -1, 0.0));
  const int no = static_cast<int>(pv.n());
  int dim[2] = {-1, -1};
  if (!shape(c["od"], 2, dim, 0) || dim[0] != no) {
    return c.fail("`od` must be a matrix with one row per origin zone (the length of `pj`).");
  }
  const int nd = dim[1];
  std::vector<double> ov(static_cast<std::size_t>(no) * static_cast<std::size_t>(nd));
  bool bad = false;
  fill(c["od"], 2, dim, 0, 0, 1, ov, bad, true);
  if (ov.empty()) return c.fail("`od` must be finite numbers from 0 with no missing values.");
  CK(nummat(c, ov, "od"));
  std::vector<double> tr(ov.size());
  for (int o = 0; o < no; ++o) {
    for (int j = 0; j < nd; ++j) {
      tr[static_cast<std::size_t>(o) * nd + j] = ov[static_cast<std::size_t>(o) + static_cast<std::size_t>(no) * j];
    }
  }
  Dv dv;
  CK(num(c, c["dest"], "dest", dv, -1, 1.0));
  std::vector<int> di(dv.n());
  for (std::size_t k = 0; k < dv.n(); ++k) {
    const double v = dv[k];
    if (v != std::floor(v) || v > nd) return c.fail("`dest` must be whole numbers from 1 to %d, the columns of `od`.", nd);
    di[k] = static_cast<int>(v) - 1;
  }
  std::vector<double> o(5 * static_cast<std::size_t>(no));
  CK(kc(c, clab_bjtc(pv.p(), no, tr.data(), nd, di.data(), static_cast<int>(dv.n()), o.data()),
        "a block sums to 0: no journey-to-crime value, no solved offence, or no solved offence in the series' zones."));
  static const char *const nm[] = {"pj", "po", "pod", "prod", "risk"};
  Frame f;
  for (int b = 0; b < 5; ++b) f.num(nm[b], o.data() + static_cast<std::size_t>(b) * no, static_cast<std::size_t>(no));
  out = f.o;
  return true;
}
