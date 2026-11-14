// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Offender profiling and risk terrain (src/brtm.cpp): the circle
// hypothesis, hit scores, risk terrain models; and the life course
// (src/blcr.cpp): the age-crime curve as a mixture, cumulative risk and
// pipelines of cases.

// rtmfit()'s layers as columns: an object of n-value columns (a data frame
// or list) or a matrix of n rows.
bool lcols(Ctx &c, const Val &lay, int n, std::vector<double> &L, std::vector<std::string> &nm) {
  if (lay.t == Val::OBJ) {
    if (lay.a.empty()) return c.fail("`lay` must hold at least one layer.");
    for (std::size_t j = 0; j < lay.a.size(); ++j) {
      Dv v;
      CK(num(c, lay.a[j], "lay", v, n));
      L.insert(L.end(), v.v.begin(), v.v.end());
      nm.push_back(lay.k[j]);
    }
    return true;
  }
  int r = 0;
  std::vector<double> v;
  CK(bmat(c, lay, "lay", -1, r, v));
  if (r != n) return c.fail("`lay` must have one row per cell, %d.", n);
  L = v;
  for (std::size_t j = 0; j < v.size() / static_cast<std::size_t>(n); ++j) nm.push_back("");
  return true;
}
// rtmfit()'s layers from features: one matrix of points, or an object of them.
bool lfeats(Ctx &c, const Val &lay, int n, std::vector<double> &L, std::vector<std::string> &nm) {
  int nr = 0;
  std::vector<double> cv;
  CK(bmat(c, asmat(c["cell"]), "cell", 2, nr, cv));
  if (nr != n) return c.fail("`cell` must have one row per count, %d.", n);
  double r;
  CK(one(c, c["rad"], "rad", 0.0, HUGE_VAL, r));
  static const char *const ko[] = {"prox", "dens", "kde"};
  int kd;
  CK(pickd(c, "kind", ko, 3, 0, kd));
  if (kd == 2 && n < 2) return c.fail("`kind = \"kde\"` needs at least two cells.");
  const bool one1 = lay.t != Val::OBJ;
  const std::size_t k = one1 ? 1 : lay.a.size();
  if (k < 1) return c.fail("`lay` must hold at least one layer.");
  L.resize(static_cast<std::size_t>(n) * k);
  for (std::size_t j = 0; j < k; ++j) {
    int nf = 0;
    std::vector<double> f;
    CK(bmat(c, asmat(one1 ? lay : lay.a[j]), "lay", 2, nf, f));
    clab_rtmlay(cv.data(), cv.data() + n, n, f.data(), f.data() + nf, nf, r, kd, &L[static_cast<std::size_t>(n) * j]);
    nm.push_back(one1 ? "" : lay.k[j]);
  }
  return true;
}

bool crchyp(Ctx &c, Val &out) {
  Dv xv, yv;
  CK(num(c, c["x"], "x", xv));
  const int n = static_cast<int>(xv.n());
  CK(num(c, c["y"], "y", yv, n));
  if (n < 2) return c.fail("`x` and `y` must give at least two sites.");
  int nh = 0;
  std::vector<double> hv;
  if (c.has("home")) {
    const Val &h = c["home"];
    if (h.t == Val::ARR && isnum(h) && h.a.size() == 2) {
      Dv d;
      CK(num(c, h, "home", d, 2));
      hv = d.v;
      nh = 1;
    } else {
      CK(bmat(c, asmat(h), "home", 2, nh, hv));
    }
  }
  double st[4];
  int pair[2];
  std::vector<int> mar(static_cast<std::size_t>(nh) + 1);
  const double *hx = nh ? hv.data() : nullptr;
  clab_crchyp(xv.p(), yv.p(), n, hx, nh ? hx + nh : nullptr, nh, st, pair, mar.data());
  static const char *const cn[] = {"x", "y"};
  const int pr[2] = {pair[0] + 1, pair[1] + 1};
  out = vobj();
  put(out, "centre", nvec(2, cn, st));
  put(out, "r", vnum(st[2]));
  put(out, "pair", vints(pr, 2));
  put(out, "share", vnum(st[3]));
  put(out, "mar", nh ? unbox(vlgls(mar.data(), static_cast<std::size_t>(nh))) : Val());
  return true;
}

bool hitscr(Ctx &c, Val &out) {
  Dv sv, hv;
  CK(num(c, c["s"], "s", sv));
  const int n = static_cast<int>(sv.n());
  CK(num(c, c["home"], "home", hv, -1, 1.0));
  const std::size_t nh = hv.n();
  std::vector<int> h(nh), hm(nh), ab(nh), ti(nh), rk(nh);
  for (std::size_t k = 0; k < nh; ++k) {
    const double v = hv[k];
    if (v != std::floor(v) || v > n) return c.fail("`home` must be cell numbers from 1 to %d.", n);
    h[k] = static_cast<int>(v) - 1;
    hm[k] = h[k] + 1;
  }
  std::vector<double> hit(nh), hitt(nh);
  clab_hitscr(sv.p(), n, h.data(), static_cast<int>(nh), ab.data(), ti.data(), hit.data(), hitt.data());
  for (std::size_t k = 0; k < nh; ++k) rk[k] = ab[k] + 1;
  Frame f;
  f.ints("home", hm.data(), nh);
  f.ints("above", ab.data(), nh);
  f.ints("ties", ti.data(), nh);
  f.ints("rank", rk.data(), nh);
  f.num("hit", hit);
  f.num("hitt", hitt);
  out = f.o;
  return true;
}

bool rtmfit(Ctx &c, Val &out) {
  Dv yv;
  CK(num(c, c["y"], "y", yv, -1, 0.0));
  const int n = static_cast<int>(yv.n());
  bool any = false;
  for (double e : yv.v) any = any || e > 0;
  if (!any) return c.fail("`y` must hold at least one count above 0.");
  double lv;
  CK(oned(c, "level", 0.95, 0.0, 1.0, lv));
  std::vector<double> L;
  std::vector<std::string> nm;
  if (c.has("cell")) {
    CK(lfeats(c, c["lay"], n, L, nm));
  } else {
    CK(lcols(c, c["lay"], n, L, nm));
  }
  const int k = static_cast<int>(nm.size());
  for (int j = 0; j < k; ++j) {
    if (nm[static_cast<std::size_t>(j)].empty()) nm[static_cast<std::size_t>(j)] = "layer" + std::to_string(j + 1);
  }
  const int p = k + 1;
  const std::size_t N = static_cast<std::size_t>(n), P = static_cast<std::size_t>(p);
  std::vector<double> x(N * P, 1.0);
  std::copy(L.begin(), L.end(), x.begin() + n);
  std::vector<double> beta(P), vb(P * P), vc(1), st(8), w(N * P + 6 * N + P + P * P);
  std::vector<int> al(P), iw(P);
  CK(kc(c, clab_cntreg(yv.p(), x.data(), n, p, nullptr, nullptr, 0, 0, 100, 1e-14, beta.data(), vb.data(), vc.data(),
                       st.data(), al.data(), w.data(), iw.data()),
        "the Poisson fit failed: the deviance left the finite range."));
  const double q = qnorm(0.5 + lv / 2);
  std::vector<std::string> tn = {"(Intercept)"};
  tn.insert(tn.end(), nm.begin(), nm.end());
  std::vector<double> est(P), se(P), z(P), pv(P), rrv(P), lo(P), hi(P);
  for (std::size_t j = 0; j < P; ++j) {
    const double b = beta[j], s = std::sqrt(vb[j * P + j]);
    est[j] = na(b);
    se[j] = na(s);
    z[j] = na(b / s);
    pv[j] = na(2 * clab_pnorm(-std::fabs(b / s)));
    rrv[j] = na(std::exp(b));
    lo[j] = na(std::exp(b - q * s));
    hi[j] = na(std::exp(b + q * s));
  }
  std::vector<double> rr(N), sc(N), fit(N);
  double mn = HUGE_VAL;
  for (std::size_t i = 0; i < N; ++i) {
    double e = 0;
    for (std::size_t j = 1; j < P; ++j) {
      if (!al[j]) e += L[(j - 1) * N + i] * beta[j];
    }
    rr[i] = std::exp(e);
    fit[i] = std::exp(e + beta[0]);
    mn = std::min(mn, rr[i]);
  }
  for (std::size_t i = 0; i < N; ++i) sc[i] = rr[i] / mn;
  Frame f, g;
  f.str("term", tn);
  f.num("est", est);
  f.num("se", se);
  f.num("z", z);
  f.num("p", pv);
  f.num("rrv", rrv);
  f.num("lo", lo);
  f.num("hi", hi);
  g.num("rr", rr);
  g.num("score", sc);
  g.num("fitted", fit);
  g.num("y", yv.v);
  out = vobj();
  put(out, "coef", f.o);
  put(out, "cells", g.o);
  put(out, "lay", vmat(L, n, k));
  put(out, "dev", vnum(st[7]));
  put(out, "aic", vnum(st[1]));
  put(out, "conv", vbool(st[5] == 1));
  put(out, "level", vnum(lv));
  return true;
}

// The local mat(): an object of columns or a matrix of finite numbers;
// column-major into v, column names into cn.
bool lcmat(Ctx &c, const Val &x, const char *nm, int &nr, int &nc, std::vector<double> &v, std::vector<std::string> &cn) {
  cn.clear();
  v.clear();
  if (x.t == Val::OBJ) {
    nc = static_cast<int>(x.a.size());
    if (nc < 1) return c.fail("`%s` must have at least one column.", nm);
    nr = static_cast<int>(vlen(x.a[0]));
    for (const Val &e : x.a) {
      Dv d;
      CK(num(c, e, nm, d, nr));
      v.insert(v.end(), d.v.begin(), d.v.end());
    }
    cn = x.k;
    return true;
  }
  int dim[2] = {-1, -1};
  if (x.t != Val::ARR || x.a.empty() || !shape(x, 2, dim, 0)) {
    return c.fail("`%s` must be a numeric matrix or a data frame of numbers.", nm);
  }
  nr = dim[0];
  nc = dim[1];
  v.assign(static_cast<std::size_t>(nr) * static_cast<std::size_t>(nc), 0.0);
  bool bad = false;
  fill(x, 2, dim, 0, 0, 1, v, bad, true);
  for (double e : v) {
    if (!std::isfinite(e)) return c.fail("`%s` must hold finite numbers.", nm);
  }
  return true;
}

bool agemix(Ctx &c, Val &out) {
  if (!c.has("curves") && !c.has("panel")) return c.fail("give `curves` and `shares`, or `panel`.");
  Val crv, cov, pan;
  if (c.has("curves")) {
    int na2, nt;
    std::vector<double> cv;
    std::vector<std::string> cn;
    CK(lcmat(c, c["curves"], "curves", na2, nt, cv, cn));
    if (na2 < 1) return c.fail("`curves` must have at least one age.");
    if (nt > 30) return c.fail("`curves` may have at most 30 types.");
    Dv w, av;
    CK(num(c, c["shares"], "shares", w, nt, 0.0));
    if (c.has("ages")) {
      CK(num(c, c["ages"], "ages", av, na2));
    } else {
      for (int a = 0; a < na2; ++a) av.v.push_back(a + 1);
    }
    const std::size_t A = static_cast<std::size_t>(na2);
    std::vector<double> ag(A), cm(A * A);
    if (clab_agemix(cv.data(), na2, nt, w.p(), ag.data(), cm.data()) != CLAB_OK) return c.fail("`shares` must sum to 1.");
    Frame f;
    f.num("age", av.v);
    f.num("agg", ag);
    for (int k = 0; k < nt; ++k) {
      const std::string nk = !cn.empty() && !cn[static_cast<std::size_t>(k)].empty() ? cn[static_cast<std::size_t>(k)]
                                                                                    : "type" + std::to_string(k + 1);
      f.num(nk.c_str(), cv.data() + A * static_cast<std::size_t>(k), A);
    }
    crv = f.o;
    cov = vmat(cm, na2, na2);
  }
  if (c.has("panel")) {
    int n, k;
    std::vector<double> pv;
    std::vector<std::string> cn;
    CK(lcmat(c, c["panel"], "panel", n, k, pv, cn));
    if (k != 2) return c.fail("`panel` must have two columns, the outcomes at two ages.");
    if (n < 3) return c.fail("`panel` must have at least three individuals.");
    int nps;
    CK(wholed(c, "nperm", 2000, 0, nps));
    unsigned int sd = 0;
    if (nps > 0) CK(seed(c, sd));
    std::vector<double> wk(static_cast<std::size_t>(n));
    double o[2];
    clab_agepan(pv.data(), pv.data() + n, n, nps, sd, o, wk.data());
    static const char *const pn[] = {"cov", "p", "n", "nperm"};
    const double pr[4] = {o[0], nps > 0 ? o[1] : NA, static_cast<double>(n), static_cast<double>(nps)};
    pan = nvec(4, pn, pr);
  }
  out = vobj();
  put(out, "curve", crv);
  put(out, "cov", cov);
  put(out, "panel", pan);
  return true;
}

bool cumrsk(Ctx &c, Val &out) {
  Dv av, ev, pv, rv, dv;
  CK(num(c, c["age"], "age", av));
  const long n = static_cast<long>(av.n());
  for (std::size_t i = 1; i < av.n(); ++i) {
    if (!(av[i] > av[i - 1])) return c.fail("`age` must be increasing.");
  }
  CK(num(c, c["ev"], "ev", ev, n, 0.0));
  CK(num(c, c["pop"], "pop", pv, n, 0.0, true));
  if (c.has("atrisk")) CK(num(c, c["atrisk"], "atrisk", rv, n, 0.0, true));
  if (c.has("dth")) {
    CK(num(c, c["dth"], "dth", dv, n, 0.0));
    for (double e : dv.v) {
      if (e >= 1) return c.fail("`dth` must be death rates from 0 to below 1.");
    }
  }
  double lv;
  CK(oned(c, "level", 0.95, 0.0, 1.0, lv));
  const double z = qnorm(0.5 + lv / 2);
  const std::size_t N = av.n();
  std::vector<double> h(N), nr(N), cr(N), lo(N), hi(N), vw(N);
  if (clab_cumrsk(av.p(), ev.p(), pv.p(), rv.n() ? rv.p() : nullptr, dv.n() ? dv.p() : nullptr, static_cast<int>(n), h.data(),
                  nr.data(), cr.data(), vw.data()) != CLAB_OK) {
    return c.fail("more first events (and deaths) than people at risk at some age.");
  }
  for (std::size_t i = 0; i < N; ++i) {
    const double s = 1 - cr[i], sv = std::sqrt(vw[i]);
    lo[i] = std::isnan(sv) ? NA : 1 - s * std::exp(z * sv);
    hi[i] = std::isnan(sv) ? NA : 1 - s * std::exp(-z * sv);
  }
  Frame f;
  f.num("age", av.v);
  f.num("hazard", h);
  f.num("atrisk", nr);
  f.num("cumrisk", cr);
  f.num("lo", lo);
  f.num("hi", hi);
  out = f.o;
  return true;
}

bool jpipe(Ctx &c, Val &out) {
  Dv nv, dv;
  CK(num(c, c["num"], "num", nv, -1, 0.0));
  const std::size_t n = nv.n();
  CK(num(c, c["den"], "den", dv, static_cast<long>(n), 0.0));
  int w, lg;
  CK(wholed(c, "win", 12, 1, w));
  CK(wholed(c, "lag", 0, 0, lg));
  Lab g;
  const bool hg = c.has("grp");
  if (hg) CK(codes(c, c["grp"], "grp", n, g));
  std::vector<int> o(n);
  for (std::size_t i = 0; i < n; ++i) o[i] = static_cast<int>(i);
  if (hg) {
    std::stable_sort(o.begin(), o.end(), [&](int a, int b) {
      return g.code[static_cast<std::size_t>(a)] < g.code[static_cast<std::size_t>(b)];
    });
  }
  std::vector<double> sn(n), sd(n), on(n), od(n), rt(n);
  for (std::size_t i = 0; i < n; ++i) {
    sn[i] = nv[static_cast<std::size_t>(o[i])];
    sd[i] = dv[static_cast<std::size_t>(o[i])];
  }
  std::vector<std::string> gl(n);
  std::vector<int> ix(n);
  auto code = [&](std::size_t i) { return hg ? g.code[static_cast<std::size_t>(o[i])] : 0; };
  for (std::size_t a = 0; a < n;) {
    std::size_t b = a;
    const int ca = code(a);
    while (b < n && code(b) == ca) ++b;
    clab_jpipe(sn.data() + a, sd.data() + a, static_cast<int>(b - a), w, lg, on.data() + a, od.data() + a, rt.data() + a);
    for (std::size_t i = a; i < b; ++i) {
      gl[i] = hg ? g.lev[static_cast<std::size_t>(ca)] : "all";
      ix[i] = static_cast<int>(i - a) + 1;
      on[i] = na(on[i]);
      od[i] = na(od[i]);
      rt[i] = na(rt[i]);
    }
    a = b;
  }
  Frame f;
  f.str("group", gl);
  f.ints("index", ix.data(), n);
  f.num("num", on);
  f.num("den", od);
  f.num("rate", rt);
  out = f.o;
  return true;
}
