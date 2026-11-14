// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Deterrence (src/bdtr.cpp): the response to sanctions, identification of
// certainty, severity and celerity, the inspection game, collider bias;
// and crime severity (src/bsev.cpp): severity-weighted counts and
// z-score composites.

bool detrsp(Ctx &c, Val &out) {
  Dv xv, bv, sv, pv;
  CK(num(c, c["x"], "x", xv));
  const long n = static_cast<long>(xv.n());
  CK(num(c, c["ben"], "ben", bv, n));
  CK(num(c, c["sanc"], "sanc", sv, n));
  CK(num(c, c["p"], "p", pv, -1, 0.0));
  std::vector<int> o(xv.n());
  for (std::size_t i = 0; i < o.size(); ++i) o[i] = static_cast<int>(i);
  std::stable_sort(o.begin(), o.end(), [&](int a, int b) {
    return xv[static_cast<std::size_t>(a)] < xv[static_cast<std::size_t>(b)];
  });
  std::vector<double> xs(xv.n()), bs(xv.n()), ss(xv.n()), xo(pv.n()), v(pv.n());
  for (std::size_t i = 0; i < xv.n(); ++i) {
    const std::size_t k = static_cast<std::size_t>(o[i]);
    xs[i] = xv[k];
    bs[i] = bv[k];
    ss[i] = sv[k];
  }
  CK(kc(c, clab_detrsp(xs.data(), bs.data(), ss.data(), static_cast<int>(n), pv.p(), static_cast<int>(pv.n()), xo.data(),
                       v.data()),
        "`x` must not repeat a level and `sanc` must increase with `x`."));
  Frame f;
  f.num("p", pv.v);
  f.num("xopt", xo);
  f.num("value", v);
  out = f.o;
  return true;
}

bool dtrid(Ctx &c, Val &out) {
  Dv pv, sv, cv, yv;
  CK(num(c, c["p"], "p", pv));
  const long n = static_cast<long>(pv.n());
  CK(num(c, c["s"], "s", sv, n));
  CK(num(c, c["c"], "c", cv, n));
  if (c.has("y")) CK(num(c, c["y"], "y", yv, n));
  std::vector<int> lg = {0};
  if (c.has("log")) CK(lgl(c, c["log"], "log", 1, lg));
  const int l = lg[0];
  double rh = 0;
  if (c.has("rho")) CK(binr(c, c["rho"], "rho", -HUGE_VAL, HUGE_VAL, rh));
  if (l) {
    for (std::size_t i = 0; i < pv.n(); ++i) {
      if (!(pv[i] > 0 && sv[i] > 0)) return c.fail("`p` and `s` must be positive with `log = TRUE`.");
    }
  }
  int rank = 0, id[3], cst[3];
  double beta[4];
  std::vector<double> w(8 * pv.n());
  if (clab_dtrid(pv.p(), sv.p(), cv.p(), static_cast<int>(n), l, yv.n() ? yv.p() : nullptr, &rank, id, cst, beta,
                 w.data()) != CLAB_OK) {
    return c.fail("the design has rank %d of 4: a dimension held fixed is not identified, so `y` cannot be fitted.",
                  rank);
  }
  static const char *const dn[] = {"certainty", "severity", "celerity"};
  Val iv = vobj(), kv = vobj();
  for (int j = 0; j < 3; ++j) {
    put(iv, dn[j], vbool(id[j] != 0));
    put(kv, dn[j], vbool(cst[j] != 0));
  }
  out = vobj();
  put(out, "rank", vnum(rank));
  put(out, "ident", iv);
  put(out, "const", kv);
  if (yv.n()) {
    static const char *const bn[] = {"(Intercept)", "p", "s", "c"};
    static const char *const ln[] = {"(Intercept)", "log_p", "log_s", "c"};
    put(out, "coef", nvec(4, l ? ln : bn, beta));
    put(out, "disc", vnum(l ? (beta[1] != 0 ? -beta[3] / beta[1] : NA) : NA));
  } else {
    put(out, "coef", Val());
    put(out, "disc", Val());
  }
  if (c.has("rho")) {
    std::vector<double> ec(pv.n());
    for (std::size_t i = 0; i < pv.n(); ++i) ec[i] = pv[i] * sv[i] * std::exp(-rh * cv[i]);
    put(out, "cost", vauto(ec));
  } else {
    put(out, "cost", Val());
  }
  return true;
}

bool inspg(Ctx &c, Val &out) {
  double av = 0.0;
  CK(binr(c, c["a"], "a", -HUGE_VAL, HUGE_VAL, av));
  if (av == 0) return c.fail("`a` must not be 0.");
  Dv xv;
  CK(num(c, c["x"], "x", xv, 3, 0.0));
  if (std::fabs(xv[0] + xv[1] + xv[2] - 1) > 1e-10) return c.fail("`x` must be three shares summing to 1.");
  int ns = 0;
  CK(wholed(c, "steps", 1000, 0, ns));
  if (ns > 10000000) return c.fail("`steps` must be at most 1e7.");
  double h = 0.0;
  CK(oned(c, "dt", 0.01, 0.0, HUGE_VAL, h));
  const int m = ns + 1;
  std::vector<double> tr(3 * static_cast<std::size_t>(m)), v(static_cast<std::size_t>(m));
  double f0[3], dv = 0;
  clab_inspg(av, xv.p(), ns, h, tr.data(), v.data(), f0, &dv);
  double drift = 0;
  for (double e : v) drift = std::max(drift, std::fabs(e - v[0]));
  static const char *const xn[] = {"x1", "x2", "x3"};
  const double third[3] = {1.0 / 3, 1.0 / 3, 1.0 / 3};
  out = vobj();
  put(out, "traj", vmat(tr, m, 3));
  put(out, "inv", vauto(v));
  put(out, "drift", vnum(drift));
  put(out, "dvdt", vnum(dv));
  put(out, "field", nvec(3, xn, f0));
  put(out, "fixed", nvec(3, xn, third));
  return true;
}

// seq(from, to, length.out = n), as R computes it.
std::vector<double> rseq(double from, double to, int n) {
  std::vector<double> v(static_cast<std::size_t>(n));
  const double by = (to - from) / (n - 1);
  for (int i = 0; i < n; ++i) v[static_cast<std::size_t>(i)] = from + i * by;
  if (n > 1) v[static_cast<std::size_t>(n) - 1] = to;
  return v;
}

bool insphs(Ctx &c, Val &out) {
  Dv tv, gv, xv;
  Val dt = vnums(rseq(0.05, 1.8, 20)), dg = vnums(rseq(0.05, 1.2, 20)), dx = varr();
  dx.a = {vnum(0.34), vnum(0.33), vnum(0.33)};
  CK(num(c, c.has("temp") ? c["temp"] : dt, "temp", tv));
  CK(num(c, c.has("cost") ? c["cost"] : dg, "cost", gv));
  const int nt = static_cast<int>(tv.n()), ng = static_cast<int>(gv.n());
  if (static_cast<double>(nt) * ng > 1e7) return c.fail("`temp` and `cost` make more than 1e7 cells.");
  int ns = 0;
  double r = 0.0;
  CK(wholed(c, "steps", 600, 0, ns));
  CK(oned(c, "rate", 0.05, 0.0, HUGE_VAL, r));
  CK(num(c, c.has("x0") ? c["x0"] : dx, "x0", xv, 3, 0.0));
  if (!(xv[0] + xv[1] + xv[2] > 0)) return c.fail("`x0` must have a positive sum.");
  const std::size_t cells = static_cast<std::size_t>(nt) * static_cast<std::size_t>(ng);
  std::vector<double> o(3 * cells);
  clab_insphs(tv.p(), nt, gv.p(), ng, ns, r, xv.p(), o.data());
  static const int at[3] = {1, 0, 2};
  static const char *const on[] = {"crime", "coop", "insp"};
  out = vobj();
  for (int k = 0; k < 3; ++k) put(out, on[k], vmat(o.data() + at[k] * cells, nt, ng));
  put(out, "temp", vauto(tv.v));
  put(out, "cost", vauto(gv.v));
  return true;
}

bool colarr(Ctx &c, Val &out) {
  Dv a0, e0, b0;
  CK(num(c, c["a"], "a", a0));
  CK(num(c, c["e"], "e", e0));
  CK(num(c, c["bg"], "bg", b0));
  const std::size_t n = std::max({a0.n(), e0.n(), b0.n()});
  std::vector<double> v[3];
  const Dv *src[3] = {&a0, &e0, &b0};
  static const char *const nm[] = {"a", "e", "bg"};
  for (int t = 0; t < 3; ++t) {
    const std::size_t k = src[t]->n();
    if (k != 1 && k != n) return c.fail("`%s` must be one number or one per scenario, %d.", nm[t], static_cast<int>(n));
    for (std::size_t i = 0; i < n; ++i) v[t].push_back((*src[t])[k == 1 ? 0 : i]);
  }
  std::vector<double> cl(4 * n), orr(n), lor(n);
  CK(kc(c, clab_colarr(v[0].data(), v[1].data(), v[2].data(), static_cast<int>(n), cl.data(), orr.data()),
        "`a` and `e` must be strictly between 0 and 1, `bg` from 0 to 1."));
  for (std::size_t i = 0; i < n; ++i) lor[i] = std::log(orr[i]);
  Frame f;
  for (int t = 0; t < 3; ++t) f.num(nm[t], v[t]);
  f.num("or", orr);
  f.num("lor", lor);
  static const char *const cn[] = {"c11", "c10", "c01", "c00"};
  for (int k = 0; k < 4; ++k) f.num(cn[k], cl.data() + static_cast<std::size_t>(k) * n, n);
  out = f.o;
  return true;
}

bool crsev(Ctx &c, Val &out) {
  Dv xv;
  CK(num(c, c["x"], "x", xv, -1, 0.0));
  const std::size_t n = xv.n();
  Lab o;
  CK(codes(c, c["off"], "off", n, o));
  std::vector<double> wv, sv;
  CK(perlev(c, c["w"], "w", o.lev, wv));
  if (c.has("sent")) {
    CK(perlev(c, c["sent"], "sent", o.lev, sv));
    for (std::size_t j = 0; j < wv.size(); ++j) wv[j] *= sv[j];
  }
  std::vector<int> zero(n, 0);
  const int *uc = zero.data(), *tc = zero.data();
  int nt = 1;
  Lab u, t;
  const bool hu = c.has("unit"), ht = c.has("time");
  if (hu) {
    CK(codes(c, c["unit"], "unit", n, u));
    uc = u.code.data();
  }
  if (ht) {
    CK(codes(c, c["time"], "time", n, t));
    tc = t.code.data();
    nt = t.g;
  }
  std::vector<long long> k(n);
  for (std::size_t i = 0; i < n; ++i) k[i] = static_cast<long long>(uc[i]) * nt + tc[i];
  std::vector<long long> s(k);
  std::sort(s.begin(), s.end());
  s.erase(std::unique(s.begin(), s.end()), s.end());
  const int nk = static_cast<int>(s.size());
  std::vector<int> key(n);
  for (std::size_t i = 0; i < n; ++i) key[i] = static_cast<int>(std::lower_bound(s.begin(), s.end(), k[i]) - s.begin());
  std::vector<double> pv;
  if (c.has("pop")) {
    Dv pr;
    CK(num(c, c["pop"], "pop", pr, static_cast<long>(n), 0.0, true));
    pv.assign(static_cast<std::size_t>(nk), -1.0);
    for (std::size_t i = 0; i < n; ++i) {
      double &e = pv[static_cast<std::size_t>(key[i])];
      if (e >= 0 && e != pr[i]) return c.fail("`pop` must be the same on every row of a unit and period.");
      e = pr[i];
    }
  }
  double pe = 0.0;
  CK(oned(c, "per", 1e5, 0.0, HUGE_VAL, pe));
  int b = 0;
  if (c.has("base")) CK(whole(c, c["base"], "base", 1, b));
  --b;
  if (b >= nk) return c.fail("`base` must be a row of the result, from 1 to %d.", nk);
  const std::size_t K = static_cast<std::size_t>(nk);
  std::vector<double> cnt(K), sev(K), rate(K), srate(K), idx(K);
  CK(kc(c, clab_crsev(xv.p(), o.code.data(), key.data(), static_cast<int>(n), o.g, nk, wv.data(), ptr(pv), pe, b,
                      cnt.data(), sev.data(), rate.data(), srate.data(), idx.data()),
        "the kernel refused the input."));
  Frame f;
  std::vector<std::string> ul, tl;
  for (std::size_t j = 0; j < K; ++j) {
    if (hu) ul.push_back(u.lev[static_cast<std::size_t>(s[j] / nt)]);
    if (ht) tl.push_back(t.lev[static_cast<std::size_t>(s[j] % nt)]);
  }
  if (hu) f.str("unit", ul);
  if (ht) f.str("time", tl);
  f.num("count", cnt);
  f.num("sev", sev);
  if (!pv.empty()) {
    f.num("pop", pv);
    f.num("rate", rate);
    f.num("srate", srate);
  }
  if (b >= 0) {
    for (double &e : idx) e = na(e);
    f.num("idx", idx);
  }
  Frame g;
  g.str("off", o.lev);
  g.num("w", wv);
  out = vobj();
  put(out, "tab", f.o);
  put(out, "w", g.o);
  put(out, "per", vnum(pe));
  put(out, "base", b < 0 ? Val() : vnum(b + 1));
  return true;
}

bool zcomp(Ctx &c, Val &out) {
  Dv xv;
  CK(num(c, c["x"], "x", xv, -1, 0.0));
  const std::size_t n = xv.n();
  Lab u, ct;
  CK(codes(c, c["unit"], "unit", n, u));
  CK(codes(c, c["cat"], "cat", n, ct));
  if (u.g < 2) return c.fail("`unit` must hold at least two units.");
  const long long cells = static_cast<long long>(u.g) * ct.g;
  if (cells > 100000000LL) return c.fail("`unit` and `cat` make more than 1e8 cells.");
  std::vector<double> xm(static_cast<std::size_t>(cells), 0.0), zm(static_cast<std::size_t>(cells));
  for (std::size_t i = 0; i < n; ++i) {
    xm[static_cast<std::size_t>(u.code[i]) + static_cast<std::size_t>(u.g) * static_cast<std::size_t>(ct.code[i])] += xv[i];
  }
  std::vector<double> wv(static_cast<std::size_t>(ct.g), 1.0);
  if (c.has("w")) CK(perlev(c, c["w"], "w", ct.lev, wv));
  std::vector<double> comp(static_cast<std::size_t>(u.g)), mu(static_cast<std::size_t>(ct.g)), sd(static_cast<std::size_t>(ct.g));
  clab_zcomp(xm.data(), u.g, ct.g, wv.data(), zm.data(), comp.data(), mu.data(), sd.data());
  std::vector<int> ord(static_cast<std::size_t>(u.g)), rk(static_cast<std::size_t>(u.g));
  for (int j = 0; j < u.g; ++j) ord[static_cast<std::size_t>(j)] = j;
  std::stable_sort(ord.begin(), ord.end(), [&](int a, int b) {
    return comp[static_cast<std::size_t>(a)] > comp[static_cast<std::size_t>(b)];
  });
  for (std::size_t j = 0; j < ord.size(); ++j) {
    const std::size_t k = static_cast<std::size_t>(ord[j]);
    rk[k] = j > 0 && comp[static_cast<std::size_t>(ord[j - 1])] == comp[k] ? rk[static_cast<std::size_t>(ord[j - 1])]
                                                                           : static_cast<int>(j) + 1;
  }
  Frame f, g;
  f.str("unit", u.lev);
  f.num("comp", comp);
  f.ints("rank", rk.data(), rk.size());
  g.str("cat", ct.lev);
  g.num("mean", mu);
  g.num("sd", sd);
  g.num("w", wv);
  out = vobj();
  put(out, "tab", f.o);
  put(out, "cat", g.o);
  put(out, "z", vmat(zm, u.g, ct.g));
  return true;
}
