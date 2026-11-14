// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Crosswalks (src/bxw.cpp): counts moved between boundaries and synthetic
// exposure; and taphonomy (src/btph.cpp): log-ratio transforms, the
// likelihood ratio of a post-mortem interval, the decomposition chain and
// the Bayesian hierarchical model of decomposition.

// The column names of a data frame given for a matrix (as.matrix() keeps
// them); none for a matrix given as rows.
std::vector<std::string> mcols(const Val &x) {
  return asmat(x).t == Val::ARR && x.t == Val::OBJ ? x.k : std::vector<std::string>();
}

bool bndxw(Ctx &c, Val &out) {
  const std::size_t m = vlen(c["from"]);
  if (m < 1) return c.fail("`from` must hold at least one link.");
  Lab f, g;
  CK(codes(c, c["from"], "from", m, f));
  CK(codes(c, c["to"], "to", m, g));
  std::vector<double> wv(m, 1.0);
  if (c.has("w")) {
    Dv d;
    CK(num(c, c["w"], "w", d, static_cast<long>(m), 0.0));
    wv = d.v;
  }
  Dv xv;
  CK(num(c, c["x"], "x", xv));
  const std::size_t n = xv.n();
  if (vlen(c["unit"]) != n) return c.fail("`unit` must have one label per value, %d.", static_cast<int>(n));
  std::vector<std::string> ul;
  if (!alabs(c["unit"], ul)) return c.fail("`unit` must be labels with no missing values.");
  std::vector<int> u(n);
  double lost = 0;
  std::vector<std::string> miss;
  for (std::size_t i = 0; i < n; ++i) {
    const auto it = std::find(f.lev.begin(), f.lev.end(), ul[i]);
    u[i] = it == f.lev.end() ? -1 : static_cast<int>(it - f.lev.begin());
    if (u[i] < 0) {
      lost += xv[i];
      if (std::find(miss.begin(), miss.end(), ul[i]) == miss.end()) miss.push_back(ul[i]);
    }
  }
  const std::size_t F = static_cast<std::size_t>(f.g);
  std::vector<double> ws(F), bw(F), o(static_cast<std::size_t>(g.g));
  std::vector<int> best(F);
  clab_bndxw(xv.p(), u.data(), static_cast<int>(n), f.code.data(), g.code.data(), wv.data(), static_cast<int>(m), f.g, g.g,
             o.data(), ws.data(), best.data(), bw.data());
  Frame tb, mp;
  tb.str("to", g.lev);
  tb.num("value", o);
  Val bt = varr();
  for (std::size_t s = 0; s < F; ++s) {
    bt.a.push_back(best[s] < 0 ? Val() : vstr(g.lev[static_cast<std::size_t>(best[s])]));
    bw[s] = na(bw[s]);
  }
  mp.str("from", f.lev);
  mp.add("to", bt);
  mp.num("share", bw);
  mp.num("wsum", ws);
  out = vobj();
  put(out, "tab", tb.o);
  put(out, "map", mp.o);
  put(out, "lost", vnum(lost));
  put(out, "unmatched", unbox(vstrs(miss)));
  return true;
}

bool synexp(Ctx &c, Val &out) {
  int n = 0, na2 = 0;
  std::vector<double> xv, av;
  const Val x = asmat(c["x"]);
  CK(bmat(c, x, "x", -1, n, xv));
  const int p = static_cast<int>(xv.size() / static_cast<std::size_t>(n));
  Dv yv;
  CK(num(c, c["y"], "y", yv, n));
  int n0 = 0, n1 = 0;
  for (double v : yv.v) {
    if (v != 0 && v != 1) return c.fail("`y` must be 0 or 1.");
    n0 += v == 0;
    n1 += v == 1;
  }
  if (!n0 || !n1) return c.fail("`y` must hold both 0 and 1.");
  if (n <= p + 1) return c.fail("`x` needs more rows than columns plus one.");
  CK(bmat(c, asmat(c["ax"]), "ax", p, na2, av));
  Dv pv;
  CK(num(c, c["pop"], "pop", pv, na2, 0.0));
  const std::size_t q = static_cast<std::size_t>(p) + 1;
  std::vector<double> cf(q), se(q), wk(2 * q * q + 4 * q);
  std::vector<int> iw(static_cast<std::size_t>(n));
  double dev = 0;
  CK(kc(c, clab_lgtreg(xv.data(), n, p, yv.p(), 1e-10, 50, cf.data(), se.data(), &dev, wk.data(), iw.data()),
        "no fit: collinear columns of `x`, or separated outcomes."));
  const std::vector<std::string> cn = mcols(c["x"]);
  std::vector<std::string> tn = {"(Intercept)"};
  for (int j = 0; j < p; ++j) tn.push_back(cn.empty() ? "x" + std::to_string(j + 1) : cn[static_cast<std::size_t>(j)]);
  std::vector<double> rate(static_cast<std::size_t>(na2)), ex(static_cast<std::size_t>(na2));
  clab_synexp(cf.data(), av.data(), na2, p, pv.p(), rate.data(), ex.data());
  double tot = 0;
  for (double e : ex) tot += e;
  Frame cc, a;
  cc.str("term", tn);
  cc.num("est", cf);
  cc.num("se", se);
  a.num("rate", rate);
  a.num("pop", pv.v);
  a.num("expo", ex);
  out = vobj();
  put(out, "coef", cc.o);
  put(out, "areas", a.o);
  put(out, "total", vnum(tot));
  put(out, "dev", vnum(dev));
  return true;
}

bool clrilr(Ctx &c, Val &out) {
  int n = 0;
  std::vector<double> xv;
  CK(bmat(c, asmat(c["x"]), "x", -1, n, xv, 0.0));
  const int d = static_cast<int>(xv.size() / static_cast<std::size_t>(n));
  if (d < 2) return c.fail("`x` must have at least two parts (columns).");
  static const char *const to[] = {"clr", "ilr"};
  int ty = 0;
  CK(pickd(c, "type", to, 2, 0, ty));
  double pc = 0.0;
  CK(binrd(c, "ps", 1e-6, 0.0, HUGE_VAL, pc));
  const int dc = ty == 0 ? d : d - 1;
  std::vector<double> o(static_cast<std::size_t>(n) * static_cast<std::size_t>(dc));
  CK(kc(c, clab_clrilr(xv.data(), n, d, pc, ty, o.data()), "a part is 0 after the pseudocount: give `ps` above 0."));
  out = vmat(o, n, dc);
  return true;
}

// oneorn(): one number or one per value, recycled to n.
bool oneorn(Ctx &c, const char *nm, std::size_t n, bool pos, std::vector<double> &o) {
  Dv d;
  if (pos) {
    CK(num(c, c[nm], nm, d, -1, 0.0, true));
  } else {
    CK(num(c, c[nm], nm, d));
  }
  if (d.n() != 1 && d.n() != n) return c.fail("`%s` must be one number or one per value, %d.", nm, static_cast<int>(n));
  o.assign(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) o[i] = d[d.n() == 1 ? 0 : i];
  return true;
}

bool preslr(Ctx &c, Val &out) {
  Dv xv;
  CK(num(c, c["x"], "x", xv));
  const std::size_t n = xv.n();
  std::vector<double> a, b, cc, d;
  CK(oneorn(c, "m1", n, false, a));
  CK(oneorn(c, "s1", n, true, b));
  CK(oneorn(c, "m2", n, false, cc));
  CK(oneorn(c, "s2", n, true, d));
  double o[5];
  int band = 0, side = 0;
  const int N = static_cast<int>(n);
  clab_preslr(xv.p(), N, a.data(), b.data(), N, cc.data(), d.data(), N, o, &band, &side);
  static const char *const bw[] = {"no support either way", "weak support", "moderate support",
                                   "moderately strong support", "strong support", "very strong support",
                                   "extremely strong support"};
  const std::string vb = std::string(bw[band]) + (side == 1 ? " for H1" : side == 2 ? " for H2" : "");
  static const char *const nm[] = {"ll1", "ll2", "llr", "lr", "l10"};
  out = vobj();
  for (int j = 0; j < 5; ++j) put(out, nm[j], vnum(o[j]));
  put(out, "verbal", vstr(vb));
  return true;
}

bool tpdmc(Ctx &c, Val &out) {
  double pr = 0.0, mr = 0.0, dr = 0.0;
  CK(p01d(c, "pres", 0, pr));
  CK(p01d(c, "mum", 0.5, mr));
  CK(binrd(c, "decay", 0.5, 0.0, 1.0, dr));
  if (!(dr > 0)) return c.fail("`decay` must be above 0.");
  std::vector<std::string> st = {"fresh", "bloat", "active", "advanced"};
  if (c.has("states")) {
    const Val &s = c["states"];
    bool ok = s.t == Val::STR || s.t == Val::ARR;
    st.clear();
    for (const Val &e : s.t == Val::STR ? std::vector<Val>{s} : s.a) {
      ok = ok && e.t == Val::STR;
      if (ok) st.push_back(e.s);
    }
    if (!ok || st.empty()) return c.fail("`states` must name the transient stages: labels, at least one.");
  }
  const int k = static_cast<int>(st.size()), m = k + 2;
  for (int i = 0; i < k; ++i) {
    for (int j = 0; j < i; ++j) {
      if (st[static_cast<std::size_t>(i)] == st[static_cast<std::size_t>(j)]) return c.fail("`states` must name each stage once.");
    }
  }
  int s0 = -1;
  const Val &sv = c.has("start") ? c["start"] : vnum(1);
  const Val *e = sv.t == Val::ARR && sv.a.size() == 1 ? &sv.a[0] : &sv;
  if (e->t == Val::STR) {
    for (int i = 0; i < k; ++i) {
      if (e->s == st[static_cast<std::size_t>(i)]) s0 = i;
    }
  } else {
    const double v = val1(sv);
    if (v == std::floor(v) && v >= 1 && v <= k) s0 = static_cast<int>(v) - 1;
  }
  if (s0 < 0) return c.fail("`start` must be a transient stage, by name or number from 1 to %d.", k);
  int ns = 0;
  CK(wholed(c, "steps", 0, 0, ns));
  const std::size_t M = static_cast<std::size_t>(m), K = static_cast<std::size_t>(k);
  std::vector<double> pm(M * M), nm(K * K), bm(K * 2), tv(K);
  CK(kc(c, clab_tpdmc(pr, dr, mr, k, pm.data(), nm.data(), bm.data(), tv.data()),
        "nothing leaves a stage: with `pres` 1, `mum` must be above 0."));
  std::vector<std::string> all = st;
  all.push_back("skeletal");
  all.push_back("mummified");
  Val path;
  if (ns > 0) {
    unsigned int sd = 0;
    CK(seed(c, sd));
    std::vector<int> p(static_cast<std::size_t>(ns) + 1);
    int len = 0;
    clab_tpdsim(pm.data(), m, k, s0, ns, sd, p.data(), &len);
    std::vector<std::string> ps;
    for (int i = 0; i < len; ++i) ps.push_back(all[static_cast<std::size_t>(p[static_cast<std::size_t>(i)])]);
    path = unbox(vstrs(ps));
  }
  out = vobj();
  put(out, "P", vmat(pm, m, m));
  put(out, "N", vmat(nm, k, k));
  put(out, "absorb", vmat(bm, k, 2));
  put(out, "time", nvec(st, tv.data()));
  put(out, "start", vnum(s0 + 1));
  put(out, "path", path);
  return true;
}

bool tpbhm(Ctx &c, Val &out) {
  Dv yv;
  CK(num(c, c["y"], "y", yv));
  const int n = static_cast<int>(yv.n());
  int nr = 0;
  std::vector<double> xv;
  CK(bmat(c, asmat(c["x"]), "x", -1, nr, xv));
  if (nr != n) return c.fail("`x` must have one row per case, %d.", n);
  const int q = static_cast<int>(xv.size() / static_cast<std::size_t>(n)), p = q + 1;
  if (n < 2) return c.fail("`y` must hold at least two cases.");
  Dv mv, sv;
  Val d1, d2;
  CK(num(c, arg(c, "pm", d1, 0), "pm", mv));
  CK(num(c, arg(c, "psd", d2, 10), "psd", sv, -1, 0.0, true));
  if ((mv.n() != 1 && mv.n() != static_cast<std::size_t>(p)) || (sv.n() != 1 && sv.n() != static_cast<std::size_t>(p))) {
    return c.fail("`pm` and `psd` must be one number or one per term, %d.", p);
  }
  double lv = 0.0;
  CK(oned(c, "level", 0.95, 0.0, 1.0, lv));
  const std::size_t P = static_cast<std::size_t>(p);
  std::vector<double> m0(P), s0(P);
  for (std::size_t j = 0; j < P; ++j) {
    m0[j] = mv[mv.n() == 1 ? 0 : j];
    s0[j] = sv[sv.n() == 1 ? 0 : j];
  }
  Lab g;
  int ng = 0;
  if (c.has("grp")) {
    CK(codes(c, c["grp"], "grp", yv.n(), g));
    ng = g.g;
    if (ng < 2) return c.fail("`grp` must hold at least two groups.");
  }
  const std::size_t G = static_cast<std::size_t>(ng) + 1;
  std::vector<double> mu(P), sd(P), st(2), gm(G), lam(G), w(2 * P * P + P);
  std::vector<int> nj(G);
  CK(kc(c, clab_tpbhm(yv.p(), xv.data(), n, q, m0.data(), s0.data(), ng ? g.code.data() : nullptr, ng, mu.data(), sd.data(),
                      st.data(), gm.data(), lam.data(), nj.data(), w.data()),
        "`x` has dependent columns (or fewer cases than terms), or `y` does not vary."));
  const double zq = qnorm(0.5 + lv / 2);
  const std::vector<std::string> cn = mcols(c["x"]);
  std::vector<std::string> tn = {"(Intercept)"};
  for (int j = 1; j < p; ++j) tn.push_back(cn.empty() ? "x" + std::to_string(j) : cn[static_cast<std::size_t>(j) - 1]);
  std::vector<double> lo(P), hi(P), pp(P);
  for (std::size_t j = 0; j < P; ++j) {
    lo[j] = mu[j] - zq * sd[j];
    hi[j] = mu[j] + zq * sd[j];
    pp[j] = clab_pnorm(mu[j] / sd[j]);
  }
  std::vector<double> fit(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    double v = mu[0];
    for (int j = 1; j < p; ++j) v += xv[static_cast<std::size_t>(j - 1) * static_cast<std::size_t>(n) + static_cast<std::size_t>(i)] * mu[static_cast<std::size_t>(j)];
    fit[static_cast<std::size_t>(i)] = v;
  }
  Frame f;
  f.str("term", tn);
  f.num("mean", mu);
  f.num("sd", sd);
  f.num("lo", lo);
  f.num("hi", hi);
  f.num("ppos", pp);
  f.num("prior", s0);
  Val gr;
  if (ng) {
    std::vector<double> po(static_cast<std::size_t>(ng));
    for (std::size_t j = 0; j < po.size(); ++j) po[j] = lam[j] * gm[j];
    Frame h;
    h.str("group", g.lev);
    h.ints("n", nj.data(), po.size());
    h.num("raw", gm.data(), po.size());
    h.num("shrink", lam.data(), po.size());
    h.num("pooled", po);
    gr = h.o;
  }
  out = vobj();
  put(out, "coef", f.o);
  put(out, "sigma", vnum(std::sqrt(st[0])));
  put(out, "tau2", vnum(na(st[1])));
  put(out, "groups", gr);
  put(out, "fitted", vauto(fit));
  put(out, "level", vnum(lv));
  return true;
}
