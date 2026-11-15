// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Spatial epidemiology (src/bsep.cpp): outbreak detection, exposure, the
// area structures, ecological regression, the graph attention bias, the
// metapopulation SIR model, the frailty and cure models and eigenvector
// filtering.

bool bayout(Ctx &c, Val &out) {
  const Val &xo = c["x"];
  Dv x;
  const Got g = readd(xo, x);
  if (g == NOTNUM || x.n() < 1) return c.fail("`x` must be a numeric vector of counts.");
  if (g == MATRIX) return notvec(c, "x");
  const int n = static_cast<int>(x.n());
  for (double e : x.v) {
    if (!std::isnan(e) && (!(e >= 0) || !std::isfinite(e))) return c.fail("`x` must be counts from 0, or NA.");
  }
  int fq = 0, bb = 0, ww = 0;
  CK(wholed(c, "freq", 52, 1, fq));
  CK(wholed(c, "b", 0, 0, bb));
  CK(wholed(c, "w", 6, 0, ww));
  bool ay = true;
  CK(flagd(c, "acty", true, ay));
  double al = 0;
  CK(oned(c, "alpha", 0.05, 0.0, 1.0, al));
  const double need = static_cast<double>(bb) * fq + ww + 1;
  std::vector<int> t;
  if (!c.has("at")) {
    if (n < need) return c.fail("`x` has %d values; the reference window needs at least %.0f.", n, need);
    for (int q = static_cast<int>(need); q <= n; ++q) t.push_back(q);
  } else {
    Dv a;
    CK(num(c, c["at"], "at", a, -1, 1.0));
    for (double e : a.v) {
      if (e != std::floor(e) || e > n) return c.fail("`at` must be times from 1 to %d.", n);
      if (e < need) return c.fail("time %.0f is too early for the reference window, which needs %.0f.", e, need);
      t.push_back(static_cast<int>(e));
    }
  }
  const int m = static_cast<int>(t.size());
  std::vector<double> ub(zz(m));
  std::vector<int> al2(zz(m));
  CK(kc(c, clab_bayout(x.p(), n, fq, bb, ww, ay, al, t.data(), m, ub.data(), al2.data()),
        "a reference window holds no observed count."));
  Frame f;
  f.ints("t", t.data(), zz(m));
  f.num("ub", ub);
  f.lgl("alarm", al2.data(), zz(m));
  out = f.o;
  return true;
}

// receptors and sources, planar or longitude and latitude
bool rsrc(Ctx &c, bool ll, Pts &r, Pts &s) {
  CK(xy(c, c["x"], c["y"], ll ? RAW : PLANAR, r));
  return xy(c, c["src"], nul, ll ? RAW : PLANAR, s);
}

bool swts(Ctx &c, int m, std::vector<double> &w) {
  w.clear();
  if (!c.has("w")) return true;
  Dv d;
  CK(num(c, c["w"], "w", d, m, 0.0));
  w = d.v;
  return true;
}

bool bufexp(Ctx &c, Val &out) {
  bool ll = false;
  CK(flagd(c, "lonlat", false, ll));
  Pts a, s;
  CK(rsrc(c, ll, a, s));
  Dv r;
  CK(num(c, c["r"], "r", r, 1, 0.0));
  std::vector<double> w;
  CK(swts(c, s.n, w));
  std::vector<double> cn(zz(a.n)), ws(cn), nr(cn);
  CK(kc(c, clab_bufexp(a.x.data(), a.y.data(), a.n, s.x.data(), s.y.data(), s.n, r[0], ptr(w), ll, cn.data(),
                       ws.data(), nr.data()),
        "the kernel refused the input."));
  Frame f;
  f.num("count", cn);
  f.num("wsum", ws);
  f.num("near", nr);
  out = f.o;
  return true;
}

bool kerexp(Ctx &c, Val &out) {
  Pts a, s;
  CK(rsrc(c, false, a, s));
  double h = 0;
  CK(one(c, c["h"], "h", 0.0, Inf, h));
  std::vector<double> w;
  CK(swts(c, s.n, w));
  static const char *const ko[] = {"gaussian", "quartic"};
  int k = 0;
  CK(pickd(c, "kern", ko, 2, 0, k));
  std::vector<double> o(zz(a.n));
  CK(kc(c, clab_kerexp(a.x.data(), a.y.data(), a.n, s.x.data(), s.y.data(), s.n, h, ptr(w), k, o.data()),
        "the kernel refused the input."));
  out = vauto(o);
  return true;
}

// nareas() of src/bsep.cpp: the areas the weights w are for
bool nareas(Ctx &c, int &n) {
  const Val &w = c["w"];
  int dim[2] = {-1, -1};
  if (w.t == Val::ARR && !w.a.empty() && w.a[0].t == Val::ARR && shape(w, 2, dim, 0)) {
    n = dim[0];
    return true;
  }
  if (w.t == Val::OBJ) {
    const Val &nb = col(w, "neighbours");
    if (nb.t == Val::ARR && col(w, "weights").t != Val::NUL) {
      n = static_cast<int>(nb.a.size());
      return true;
    }
    const double v = val1(col(w, "n"));
    if (v >= 1 && v <= 1e9 && v == std::floor(v)) {
      n = static_cast<int>(v);
      return true;
    }
  }
  return c.fail("`w` must be the list spwt() returns, an spdep listw or a matrix.");
}

bool lerprc(Ctx &c, Val &out) {
  int n = 0;
  CK(nareas(c, n));
  Wcsr w;
  CK(wreadx(c, c["w"], n, w));
  double r = 0, t = 0;
  CK(p01(c, c["rho"], "rho", r));
  CK(oned(c, "tau", 1, 0.0, Inf, t));
  std::vector<double> q(zz(n) * zz(n));
  CK(kc(c, clab_lerprc(w.rp.data(), w.ci.data(), ptr(w.w), n, r, t, q.data()), "the kernel refused the input."));
  out = vmat(q, n, n);
  return true;
}

bool bym2sc(Ctx &c, Val &out) {
  int n = 0;
  CK(nareas(c, n));
  if (n < 2) return c.fail("`w` must link at least 2 areas.");
  Wcsr w;
  CK(wreadx(c, c["w"], n, w));
  double t = 0, ph = 0, tl = 0;
  CK(oned(c, "tau", 1, 0.0, Inf, t));
  CK(p01d(c, "phi", 0.5, ph));
  CK(oned(c, "tol", 1e-10, 0.0, 1.0, tl));
  std::vector<double> g(zz(n) * zz(n)), cv(g);
  double s = 0;
  CK(kc(c, clab_bym2(w.rp.data(), w.ci.data(), ptr(w.w), n, t, ph, tl, &s, g.data(), cv.data()),
        "every area needs a neighbour."));
  out = vobj();
  put(out, "scale", vnum(s));
  put(out, "ginv", vmat(g, n, n));
  put(out, "cov", vmat(cv, n, n));
  return true;
}

bool ecolrg(Ctx &c, Val &out) {
  Dv y, e;
  CK(num(c, c["y"], "y", y, -1, 0.0));
  const int n = static_cast<int>(y.n());
  CK(num(c, c["e"], "e", e, n, 0.0, true));
  int p = 0, q = 0;
  std::vector<double> x, z;
  CK(zmat(c, "x", n, false, p, x));
  static const char *const fo[] = {"poisson", "negbin", "zip"};
  int fam = 0;
  CK(pickd(c, "family", fo, 3, 0, fam));
  CK(zmat(c, "z", n, false, q, z));
  double tl = 0;
  CK(oned(c, "tol", 1e-12, 0.0, 1.0, tl));
  int mi = 0;
  CK(wholed(c, "maxit", 500, 1, mi));
  if (n <= p + 1) return c.fail("`y` must have more values than the model has coefficients.");
  std::vector<double> cf(zz(p) + 1), zc(zz(q) + 1), fit(zz(n));
  double th = 0, ll = 0;
  CK(kc(c, clab_ecolrg(y.p(), e.p(), n, ptr(x), p, fam, ptr(z), q, tl, mi, cf.data(), zc.data(), &th, fit.data(),
                       &ll),
        "the design is singular, or the fit diverged."));
  out = vobj();
  put(out, "coef", vauto(cf));
  if (fam == 1) put(out, "theta", vnum(th));
  if (fam == 2) put(out, "zcoef", vauto(zc));
  put(out, "fitted", vauto(fit));
  put(out, "loglik", vnum(ll));
  return true;
}

bool spbias(Ctx &c, Val &out) {
  const Val &d = c["d"];
  int dim[2] = {-1, -1};
  if (!(d.t == Val::ARR && !d.a.empty() && d.a[0].t == Val::ARR && shape(d, 2, dim, 0) && isnum(d.a[0]))) {
    return c.fail("`d` must be a numeric matrix of hop distances.");
  }
  std::vector<double> dv(zz(dim[0]) * zz(dim[1]));
  bool bad = false;
  fill(d, 2, dim, 0, 0, 1, dv, bad, true);
  std::vector<int> di(dv.size());
  for (std::size_t q = 0; q < dv.size(); ++q) {
    const double e = dv[q];
    if (!std::isfinite(e) || e != std::floor(e) || e < -1 || e > 2147483647.0) {
      return c.fail("`d` must hold whole numbers from 0, or -1 for an unreachable pair.");
    }
    di[q] = static_cast<int>(e);
  }
  Dv b, u;
  CK(num(c, c["b"], "b", b));
  Val dd;
  CK(num(c, arg(c, "ub", dd, -10), "ub", u, 1));
  std::vector<double> o(dv.size());
  CK(kc(c, clab_spbias(di.data(), dim[0], dim[1], b.p(), static_cast<int>(b.n()), u[0], o.data()),
        "the kernel refused the input."));
  out = vobj();
  put(out, "bias", vmat(o, dim[0], dim[1]));
  put(out, "ub", vnum(u[0]));
  return true;
}

bool spsir(Ctx &c, Val &out) {
  Dv iv, pv;
  CK(num(c, c["i0"], "i0", iv, -1, 0.0));
  const int n = static_cast<int>(iv.n());
  if (c.has("pop")) {
    CK(num(c, c["pop"], "pop", pv, n, 0.0, true));
  } else {
    pv.v.assign(zz(n), 1000.0);
  }
  for (int i = 0; i < n; ++i) {
    if (iv[zz(i)] > pv[zz(i)]) return c.fail("`i0` must not exceed `pop`.");
  }
  Wcsr w;
  const bool hw = c.has("w");
  if (hw) CK(wreadx(c, c["w"], n, w));
  Dv bv, gv;
  Val d1, d2;
  CK(num(c, arg(c, "beta", d1, 0.5), "beta", bv, 1, 0.0));
  CK(num(c, arg(c, "gamma", d2, 0.2), "gamma", gv, 1, 0.0));
  double k = 0, tm = 0, h = 0;
  CK(p01d(c, "kappa", 0.1, k));
  CK(oned(c, "tmax", 100, 0.0, Inf, tm));
  CK(oned(c, "dt", 0.1, 0.0, Inf, h));
  const double st = std::nearbyint(tm / h);
  if (st > 1e6) return c.fail("`tmax` / `dt` must be at most 1e6 steps.");
  const int steps = static_cast<int>(st);
  const std::size_t T = zz(steps) + 1;
  std::vector<double> S(T * zz(n)), I(S), R(S);
  CK(kc(c, clab_spsir(iv.p(), n, pv.p(), hw ? w.rp.data() : nullptr, hw ? w.ci.data() : nullptr,
                      hw ? ptr(w.w) : nullptr, bv[0], gv[0], k, h, steps, S.data(), I.data(), R.data()),
        "the kernel refused the input."));
  std::vector<double> t(T);
  double rs = 0, ns = 0;
  for (std::size_t q = 0; q < T; ++q) t[q] = static_cast<double>(q) * h;
  for (int i = 0; i < n; ++i) {
    rs += R[zz(steps) + T * zz(i)];
    ns += pv[zz(i)];
  }
  out = vobj();
  put(out, "attack", vnum(rs / ns));
  put(out, "r0", vnum(bv[0] / gv[0]));
  put(out, "t", vauto(t));
  put(out, "S", vmat(S, static_cast<int>(T), n));
  put(out, "I", vmat(I, static_cast<int>(T), n));
  put(out, "R", vmat(R, static_cast<int>(T), n));
  return true;
}

// the frailty (kind 0) and cure (kind 1) models
bool surv(Ctx &c, Val &out, int kind) {
  Dv t;
  CK(num(c, c["time"], "time", t, -1, 0.0, true));
  const int n = static_cast<int>(t.n());
  std::vector<int> ev;
  CK(bin01(c, c["event"], "event", zz(n), ev));
  std::vector<double> d(ev.begin(), ev.end());
  int p = 0;
  std::vector<double> x;
  CK(zmat(c, "x", n, false, p, x));
  int nr = 0;
  CK(nareas(c, nr));
  Wcsr w;
  CK(wreadx(c, c["w"], nr, w));
  Dv rg;
  CK(num(c, c["reg"], "reg", rg, n, 1.0));
  std::vector<int> reg(zz(n));
  for (int i = 0; i < n; ++i) {
    const double v = rg[zz(i)];
    if (v != std::floor(v) || v > nr) return c.fail("`reg` must be area numbers from 1 to %d, the rows of `w`.", nr);
    reg[zz(i)] = static_cast<int>(v) - 1;
  }
  Dv tu;
  Val d1;
  CK(num(c, arg(c, "tau", d1, 1), "tau", tu, 1, 0.0));
  double tl = 0;
  CK(oned(c, "tol", 1e-10, 0.0, 1.0, tl));
  const int h = kind == 0 ? 1 : 2, k = h + p + nr;
  std::vector<double> th(zz(k)), se(zz(k)), cure(zz(n));
  double f = 0;
  int it = 0;
  const int rc = kind == 0 ? clab_sfrail(t.p(), d.data(), n, ptr(x), p, reg.data(), w.rp.data(), w.ci.data(),
                                         ptr(w.w), nr, tu[0], tl, th.data(), se.data(), &f, &it)
                           : clab_sfcure(t.p(), d.data(), n, ptr(x), p, reg.data(), w.rp.data(), w.ci.data(),
                                         ptr(w.w), nr, tu[0], tl, th.data(), se.data(), &f, cure.data(), &it);
  CK(kc(c, rc, "the Hessian is singular: the model is not identified by these data."));
  const std::vector<double> cf(th.begin() + h, th.begin() + h + p), u(th.begin() + h + p, th.end());
  out = vobj();
  put(out, "pll", vnum(f));
  put(out, "shape", vnum(th[0]));
  if (kind == 1) put(out, "lscale", vnum(th[1]));
  put(out, "coef", p == 0 ? varr() : vauto(cf));
  put(out, kind == 0 ? "frail" : "effect", vauto(u));
  put(out, "se", vauto(se));
  if (kind == 1) put(out, "cure", vauto(cure));
  put(out, "iter", vnum(it));
  put(out, "tau", vnum(tu[0]));
  return true;
}
bool sfrail(Ctx &c, Val &out) { return surv(c, out, 0); }
bool sfcure(Ctx &c, Val &out) { return surv(c, out, 1); }

bool evfilt(Ctx &c, Val &out) {
  Dv y;
  CK(num(c, c["v"], "v", y));
  const int n = static_cast<int>(y.n());
  if (n < 3) return c.fail("`v` must have at least 3 values.");
  Wcsr w;
  CK(wreadx(c, c["w"], n, w));
  int p = 0;
  std::vector<double> x;
  CK(zmat(c, "x", n, false, p, x));
  static const char *const mo[] = {"moran", "aic", "press", "r2"};
  static const char *const co[] = {"positive", "negative", "all"};
  int cr = 0, ca = 0;
  CK(pickd(c, "method", mo, 4, 0, cr));
  CK(pickd(c, "cand", co, 3, 0, ca));
  double tl = 0;
  CK(oned(c, "tol", 0.1, 0.0, Inf, tl));
  int mv = -1;
  if (c.has("maxv")) CK(whole(c, c["maxv"], "maxv", 0, mv));
  std::vector<int> sel(zz(n));
  int ns = 0;
  std::vector<double> cf(zz(p + n)), res(zz(n)), mi(zz(n - 1));
  double st[4];
  CK(kc(c, clab_evfilt(y.p(), n, ptr(x), p, w.rp.data(), w.ci.data(), ptr(w.w), cr, ca, tl, mv, sel.data(), &ns,
                       cf.data(), st, res.data(), mi.data()),
        "the weights have no links, or the design is singular."));
  cf.resize(zz(p + 1 + ns));
  sel.resize(zz(ns));
  out = vobj();
  put(out, "sel", ns == 1 ? vnum(sel[0]) : vints(sel));
  put(out, "coef", vauto(cf));
  static const char *const nm[] = {"imoran", "aic", "r2", "fr2"};
  for (int j = 0; j < 4; ++j) put(out, nm[j], vnum(st[j]));
  put(out, "resid", vauto(res));
  put(out, "mi", vauto(mi));
  return true;
}
