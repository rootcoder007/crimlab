// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Concentration of crime (src/bcnc.cpp): concentration measures and the
// tail of a count distribution; and the dark figure (src/bdkf.cpp):
// capture-recapture and the bounds on unrecorded crime.

// The units, total and Gini of each group.
Val bygini(const std::vector<double> &x, const double *w, const std::vector<int> &code, const Lab &g) {
  const std::size_t n = x.size(), G = static_cast<std::size_t>(g.g);
  std::vector<double> un(G), tot(G), gi(G);
  std::vector<double> xs(n + 1), ws(n + 1), lp(n + 2), lv(n + 2);
  std::vector<int> o(n + 1);
  double st[16];
  for (std::size_t k = 0; k < G; ++k) {
    std::size_t m = 0;
    for (std::size_t i = 0; i < n; ++i) {
      if (static_cast<std::size_t>(code[i]) != k) continue;
      xs[m] = x[i];
      ws[m++] = w ? w[i] : 1.0;
    }
    if (clab_cconc(xs.data(), ws.data(), static_cast<int>(m), nullptr, 0, st, lp.data(), lv.data(), nullptr, nullptr,
                   o.data()) == CLAB_OK) {
      un[k] = st[0];
      tot[k] = st[1];
      gi[k] = st[10];
    } else {
      double tw = 0;
      for (std::size_t j = 0; j < m; ++j) tw += ws[j];
      un[k] = tw;
      tot[k] = 0;
      gi[k] = NA;
    }
  }
  Frame f;
  f.str("group", g.lev);
  f.num("units", un);
  f.num("total", tot);
  f.num("gini", gi);
  return f.o;
}

bool cconc(Ctx &c, Val &out) {
  const Val &x = c["x"];
  const bool labs = !isnum(x);
  std::vector<double> xv, wv;
  Lab p;
  std::size_t n, ne = 0;
  if (labs) {
    ne = vlen(x);
    if (ne < 2) return c.fail("`x` must hold at least two events when it gives their places.");
    if (c.has("w")) return c.fail("`w` goes with counts: with places of events, each event counts once.");
    if (c.has("pop")) return c.fail("`pop` goes with counts per unit, not with places of events.");
    CK(codes(c, x, "x", ne, p));
    n = static_cast<std::size_t>(p.g);
    xv.assign(n, 0.0);
    for (int e : p.code) xv[static_cast<std::size_t>(e)] += 1;
  } else {
    Dv d;
    CK(num(c, x, "x", d, -1, 0.0));
    xv = d.v;
    n = xv.size();
    if (c.has("w")) {
      Dv w;
      CK(num(c, c["w"], "w", w, static_cast<long>(n), 0.0, true));
      wv = w.v;
    }
  }
  Dv fr;
  Val df = vnums(std::vector<double>{0.01, 0.05, 0.1, 0.25, 0.5});
  CK(num(c, c.has("fracs") ? c["fracs"] : df, "fracs", fr, -1, 0.0));
  for (double e : fr.v) {
    if (e > 1) return c.fail("`fracs` must be shares from 0 to 1.");
  }
  int ns = 0;
  CK(wholed(c, "nsim", 2000, 0, ns));
  const std::size_t nf = fr.n();
  std::vector<double> st(16), lp(n + 1), lv(n + 1), ts(nf), tn(nf);
  std::vector<int> ow(n + 1);
  if (clab_cconc(xv.data(), ptr(wv), static_cast<int>(n), fr.p(), static_cast<int>(nf), st.data(), lp.data(), lv.data(),
                 ts.data(), tn.data(), ow.data()) != CLAB_OK) {
    return c.fail("`x` has no events: every count is 0.");
  }
  static const char *const sn[] = {"units", "total", "mean", "var", "disp", "ivar", "isd", "zero",
                                   "znull", "zgap", "gini", "gpos", "gnull", "vpop", "cvb", "nullb"};
  std::vector<double> sv(16);
  for (std::size_t k = 0; k < 16; ++k) sv[k] = na(st[k]);
  Val nres;
  if (ns > 0) {
    const double un = st[0];
    if (std::fabs(un - std::round(un)) > 1e-9 || un > 1e8) {
      return c.fail("`nsim` needs whole-number weights (and at most 1e8 units).");
    }
    unsigned int sd;
    CK(seed(c, sd));
    const int nu = static_cast<int>(std::round(un));
    std::vector<double> g(static_cast<std::size_t>(ns));
    std::vector<int> wk(static_cast<std::size_t>(nu) + 1);
    CK(kc(c, clab_cnull(st[2], nu, ns, sd, g.data(), wk.data()), "the mean count is above 1e8: no null to simulate."));
    const double gobs = st[10];
    int ge = 0;
    double sum = 0;
    for (double e : g) {
      ge += e >= gobs - 1e-12;
      sum += e;
    }
    std::sort(g.begin(), g.end());
    const double h = 0.95 * (ns - 1);
    const int lo = static_cast<int>(std::floor(h));
    const double q = g[static_cast<std::size_t>(lo)] +
                     (h - lo) * (g[static_cast<std::size_t>(std::min(lo + 1, ns - 1))] - g[static_cast<std::size_t>(lo)]);
    const double nv[5] = {static_cast<double>(ns), sum / ns, q, (1.0 + ge) / (ns + 1.0), gobs - sum / ns};
    static const char *const nn[] = {"nsim", "mean", "q95", "p", "excess"};
    nres = nvec(5, nn, nv);
  }
  Val byf;
  if (c.has("by")) {
    Lab g;
    if (labs) {
      CK(codes(c, c["by"], "by", ne, g));
      std::vector<long long> key(ne);
      for (std::size_t i = 0; i < ne; ++i) key[i] = static_cast<long long>(g.code[i]) * static_cast<long long>(n) + p.code[i];
      std::sort(key.begin(), key.end());
      std::vector<double> cx;
      std::vector<int> cg;
      for (std::size_t i = 0; i < ne; ++i) {
        if (i == 0 || key[i] != key[i - 1]) {
          cx.push_back(0);
          cg.push_back(static_cast<int>(key[i] / static_cast<long long>(n)));
        }
        cx.back() += 1;
      }
      byf = bygini(cx, nullptr, cg, g);
    } else {
      CK(codes(c, c["by"], "by", n, g));
      byf = bygini(xv, ptr(wv), g.code, g);
    }
  }
  Val scl;
  if (c.has("pop")) {
    Dv pv;
    CK(num(c, c["pop"], "pop", pv, static_cast<long>(n), 0.0, true));
    double sc[5];
    if (clab_urbsc(xv.data(), pv.p(), static_cast<int>(n), sc) != CLAB_OK) {
      for (int k = 0; k < 4; ++k) sc[k] = NA;
    }
    static const char *const cn[] = {"beta", "pref", "se", "r2", "used"};
    scl = nvec(5, cn, sc);
  }
  Val grw;
  if (labs) {
    std::vector<int> ix(ne), sw(n + 1);
    std::vector<double> di(ne), ex(ne), glo(ne), ghi(ne);
    double mh = 0.0, sl = 0.0;
    clab_cgrow(p.code.data(), static_cast<int>(ne), static_cast<int>(n), di.data(), &mh, ex.data(), glo.data(), ghi.data(),
               &sl, sw.data());
    for (std::size_t i = 0; i < ne; ++i) {
      ix[i] = static_cast<int>(i) + 1;
      ex[i] = na(ex[i]);
      glo[i] = na(glo[i]);
      ghi[i] = na(ghi[i]);
    }
    Frame gc;
    gc.ints("i", ix.data(), ne);
    gc.num("dist", di);
    gc.num("expd", ex);
    gc.num("lo", glo);
    gc.num("hi", ghi);
    grw = vobj();
    put(grw, "mhat", vnum(mh));
    put(grw, "slope", vnum(na(sl)));
    put(grw, "curve", gc.o);
  }
  Frame tc, lc;
  tc.num("frac", fr.v);
  tc.num("units", tn);
  tc.num("share", ts);
  lc.num("pop", lp);
  lc.num("share", lv);
  out = vobj();
  put(out, "stat", nvec(16, sn, sv.data()));
  put(out, "top", tc.o);
  put(out, "lorenz", lc.o);
  put(out, "null", nres);
  put(out, "by", byf);
  put(out, "scale", scl);
  put(out, "grow", grw);
  put(out, "places", labs ? unbox(vstrs(p.lev)) : Val());
  return true;
}

bool hilltl(Ctx &c, Val &out) {
  Dv xv;
  CK(num(c, c["x"], "x", xv, -1, 0.0));
  std::vector<double> pv;
  for (double e : xv.v) {
    if (e > 0) pv.push_back(e);
  }
  if (pv.empty()) return c.fail("`x` has no values above 0 to fit a tail to.");
  const int n = static_cast<int>(pv.size());
  bool dis = false, apx = false;
  CK(flagd(c, "discrete", true, dis));
  CK(flagd(c, "approx", false, apx));
  int mt = 0;
  CK(wholed(c, "mintl", 3, 1, mt));
  mt = std::max(2, mt);
  double xm = 0.0;
  if (!c.has("xmin")) {
    std::vector<double> s(pv);
    std::sort(s.begin(), s.end(), [](double a, double b) { return a > b; });
    xm = n >= 50 ? s[49] : s[static_cast<std::size_t>(n) - 1];
  } else {
    CK(one(c, c["xmin"], "xmin", 0.0, HUGE_VAL, xm));
  }
  std::vector<double> wk(pv.size());
  double o[5];
  const int rc = clab_hilltl(pv.data(), n, xm, dis, apx, o, wk.data());
  const int ntl = static_cast<int>(o[3]);
  const char *meth;
  if (ntl < mt) {
    meth = "too few tail observations";
  } else if (rc != CLAB_OK) {
    meth = "not estimable: every tail value is the same";
  } else if (dis && !apx) {
    meth = o[4] > 0 ? "discrete maximum likelihood (alpha at an end of the search)" : "discrete maximum likelihood";
  } else {
    meth = dis ? "continuity-corrected Hill estimator (approximate)" : "continuous Hill estimator";
  }
  const bool fit = ntl >= mt && rc == CLAB_OK;
  out = vobj();
  put(out, "alpha", vnum(fit ? o[0] : NA));
  put(out, "se", vnum(fit ? na(o[1]) : NA));
  put(out, "xmin", vnum(xm));
  put(out, "ntl", vnum(ntl));
  put(out, "ks", vnum(fit ? o[2] : NA));
  put(out, "reliable", vbool(fit && ntl >= 50 && o[4] == 0));
  put(out, "method", vstr(meth));
  return true;
}

bool crcap(Ctx &c, Val &out) {
  Dv xv;
  CK(num(c, c["x"], "x", xv, -1, 0.0));
  const int k = static_cast<int>(xv.n());
  if (k != 3 && k != 7) return c.fail("`x` must be c(n1, n2, m) for two lists or the seven cells of three.");
  double kp = 0.0;
  CK(binrd(c, "kappa", 1, 1.0, HUGE_VAL, kp));
  double cl[7];
  static const char *const two[] = {"n1", "n2", "m"};
  static const char *const three[] = {"100", "010", "001", "110", "101", "011", "111"};
  const char *const *nm = k == 3 ? two : three;
  const bool need = k == 7;
  if (xv.nm.empty()) {
    if (need) return c.fail("`x` must name its seven cells 100, 010, 001, 110, 101, 011 and 111.");
    for (int j = 0; j < k; ++j) cl[j] = xv[static_cast<std::size_t>(j)];
  } else {
    for (int j = 0; j < k; ++j) {
      int hit = -1;
      for (int i = 0; i < k; ++i) {
        if (xv.nm[static_cast<std::size_t>(i)] == nm[j]) {
          if (hit >= 0) return c.fail("`x` names %s twice.", nm[j]);
          hit = i;
        }
      }
      if (hit < 0) {
        if (need) return c.fail("`x` must name its seven cells 100, 010, 001, 110, 101, 011 and 111.");
        return c.fail("`x` must be named n1, n2 and m, or not named at all.");
      }
      cl[j] = xv[static_cast<std::size_t>(hit)];
    }
  }
  Dv iv;
  if (c.has("int")) {
    if (k == 3) return c.fail("`int` goes with three lists.");
    CK(num(c, c["int"], "int", iv));
  }
  const int ni = static_cast<int>(iv.n());
  double st[6], pw[18];
  std::vector<double> sv(2 * iv.n() + 1);
  if (clab_crcap(cl, k, kp, iv.n() ? iv.p() : nullptr, ni, st, pw, sv.data()) != CLAB_OK) {
    return c.fail("%s", k == 3 ? "`x` must have an overlap m above 0 and no larger than either list."
                               : "the three cells of two lists only (110, 101, 011) must be above 0.");
  }
  out = vobj();
  if (k == 3) {
    static const char *const en[] = {"petersen", "lower", "upper", "chapman", "se", "floor"};
    put(out, "type", vstr("two"));
    put(out, "est", nvec(6, en, st));
    put(out, "pairs", Val());
    put(out, "sens", Val());
    return true;
  }
  static const char *const en[] = {"observed", "missing", "total"};
  put(out, "type", vstr("three"));
  put(out, "est", nvec(3, en, st));
  Frame p;
  p.str("lists", {"1-2", "1-3", "2-3"});
  static const char *const pn[] = {"n1", "n2", "m", "floor", "petersen", "chapman"};
  for (int q = 0; q < 6; ++q) p.num(pn[q], pw + 3 * q, 3);
  put(out, "pairs", p.o);
  if (ni > 0) {
    Frame s;
    s.num("int", iv.v);
    s.num("missing", sv.data(), iv.n());
    s.num("total", sv.data() + ni, iv.n());
    put(out, "sens", s.o);
  } else {
    put(out, "sens", Val());
  }
  return true;
}

bool dkfig(Ctx &c, Val &out) {
  if (!c.has("v") && !c.has("va") && !c.has("k")) {
    return c.fail("give `v` and `r` (survey and recorded rates), `va`, `ra` and `rb`, or `k`.");
  }
  out = vobj();
  Val bnd, cmp, hier;
  if (c.has("v")) {
    Dv vv, r0;
    CK(num(c, c["v"], "v", vv, -1, 0.0));
    const std::size_t n = vv.n();
    if (!c.has("r")) return c.fail("`r` must give the recorded rate with `v`.");
    CK(num(c, c["r"], "r", r0, -1, 0.0));
    if (r0.n() != 1 && r0.n() != n) return c.fail("`r` must have length 1 or the length of `v`, %d.", static_cast<int>(n));
    for (std::size_t i = 0; i < n; ++i) {
      if (vv[i] > 1) return c.fail("`v` must be rates from 0 to 1.");
      if (r0[r0.n() == 1 ? 0 : i] > 1) return c.fail("`r` must be rates from 0 to 1.");
    }
    double a = 0.0, b = 0.0;
    CK(binrd(c, "amax", 0, 0.0, 1.0, a));
    CK(binrd(c, "bmax", 0, 0.0, 1.0, b));
    if (a + b >= 1) return c.fail("`amax + bmax` must be below 1.");
    double th = NaN;
    if (c.has("thr")) {
      CK(one(c, c["thr"], "thr", 0.0, HUGE_VAL, th));
      for (double e : vv.v) {
        if (!(e > 0 && e < th)) return c.fail("`thr` must lie above every `v`, and each `v` above 0.");
      }
    }
    std::vector<double> rr(n), lo(n), hi(n), dlo(n), dhi(n), ru(n), bk(n);
    for (std::size_t i = 0; i < n; ++i) rr[i] = r0[r0.n() == 1 ? 0 : i];
    clab_dkfig(vv.p(), rr.data(), static_cast<int>(n), a, b, th, lo.data(), hi.data(), bk.data());
    for (std::size_t i = 0; i < n; ++i) {
      dlo[i] = lo[i] - rr[i];
      dhi[i] = hi[i] - rr[i];
      ru[i] = rr[i] > 0 ? hi[i] / rr[i] : NA;
    }
    Frame f;
    f.num("v", vv.v);
    f.num("r", rr);
    f.num("lo", lo);
    f.num("hi", hi);
    f.num("dlo", dlo);
    f.num("dhi", dhi);
    f.num("ratio", ru);
    if (c.has("thr")) f.num("brk", bk);
    bnd = f.o;
  }
  if (c.has("va")) {
    if (!c.has("ra") || !c.has("rb")) return c.fail("`ra` and `rb` must come with `va`.");
    double a = 0.0, b = 0.0, cc = 0.0, l = 0.0;
    CK(one(c, c["va"], "va", 0.0, HUGE_VAL, a));
    CK(one(c, c["ra"], "ra", 0.0, HUGE_VAL, b));
    CK(binr(c, c["rb"], "rb", 0.0, HUGE_VAL, cc));
    CK(binrd(c, "lam", 1, 0.0, 1.0, l));
    if (l == 0) return c.fail("`lam` must be above 0.");
    if (b > a) return c.fail("`ra` must not exceed `va`: a record is an offence.");
    double o[7];
    clab_dkcmp(a, b, cc, l, o);
    static const char *const cn[] = {"rate", "rateup", "totlo", "tothi", "dklo", "dkhi", "width"};
    cmp = nvec(7, cn, o);
  }
  if (c.has("k")) {
    Dv kv;
    CK(num(c, c["k"], "k", kv, -1, 1.0));
    double o[6];
    if (clab_dkhier(kv.p(), static_cast<int>(kv.n()), o) != CLAB_OK) {
      return c.fail("`k` must be whole numbers of offences from 1.");
    }
    static const char *const hn[] = {"incid", "offen", "ratio", "extra", "lower", "upper"};
    hier = nvec(6, hn, o);
  }
  put(out, "bnd", bnd);
  put(out, "cmp", cmp);
  put(out, "hier", hier);
  return true;
}
