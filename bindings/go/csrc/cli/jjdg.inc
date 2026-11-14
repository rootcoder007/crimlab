// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Courts (src/bjdg.cpp): judge designs as instruments, and the share of
// sentences at a mandatory minimum.

// dec(): 0 or 1 decisions as numbers.
bool dec01(Ctx &c, const Val &x, const char *nm, std::size_t n, std::vector<double> &d) {
  std::vector<int> b;
  CK(bin01(c, x, nm, n, b));
  d.assign(b.begin(), b.end());
  return true;
}
// One of the two columns of x: an object of two columns, or a matrix of
// two columns given as rows.
bool col2(Ctx &c, const Val &x, int j, const char *nm, Val &v) {
  if (x.t == Val::OBJ) {
    if (x.a.size() != 2) return c.fail("`%s` must have two columns.", nm);
    v = x.a[static_cast<std::size_t>(j)];
    return true;
  }
  int dim[2] = {-1, -1};
  if (x.t != Val::ARR || x.a.empty() || !shape(x, 2, dim, 0) || dim[1] != 2) return c.fail("`%s` must have two columns.", nm);
  v = varr();
  for (const Val &r : x.a) {
    if (r.a[static_cast<std::size_t>(j)].t == Val::STR) return c.fail("`%s` must hold numbers.", nm);
    v.a.push_back(r.a[static_cast<std::size_t>(j)]);
  }
  return true;
}

bool jdgiv(Ctx &c, Val &out) {
  out = vobj();
  if (!c.has("z")) {
    Val d0, d1, y0, y1;
    CK(col2(c, c["d"], 0, "d", d0));
    CK(col2(c, c["d"], 1, "d", d1));
    CK(col2(c, c["y"], 0, "y", y0));
    CK(col2(c, c["y"], 1, "y", y1));
    const std::size_t n = vlen(d0);
    std::vector<double> a, b;
    CK(dec01(c, d0, "d", n, a));
    CK(dec01(c, d1, "d", n, b));
    Dv u, v, wv;
    CK(num(c, y0, "y", u, static_cast<long>(n)));
    CK(num(c, y1, "y", v, static_cast<long>(n)));
    if (c.has("w")) CK(num(c, c["w"], "w", wv, static_cast<long>(n), 0.0));
    double o[12];
    if (clab_jdgpop(a.data(), b.data(), u.p(), v.p(), wv.n() ? wv.p() : nullptr, static_cast<int>(n), o) != CLAB_OK) {
      return c.fail("`w` must have a positive total.");
    }
    static const char *const sn[] = {"complier", "defier", "always", "never"};
    static const char *const en[] = {"itt", "fs", "wald", "ate"};
    for (int k = 4; k < 8; ++k) o[k] = na(o[k]);
    put(out, "type", vstr("pop"));
    put(out, "est", nvec(4, en, o + 8));
    put(out, "shares", nvec(4, sn, o));
    put(out, "effects", nvec(4, sn, o + 4));
    return true;
  }
  const std::size_t n = vlen(c["z"]);
  std::vector<double> zd, dd;
  CK(dec01(c, c["z"], "z", n, zd));
  CK(dec01(c, c["d"], "d", n, dd));
  Dv yv, wv, sv, ev;
  CK(num(c, c["y"], "y", yv, static_cast<long>(n)));
  if (c.has("w")) CK(num(c, c["w"], "w", wv, static_cast<long>(n), 0.0));
  Val ds = vnums(std::vector<double>{0, 0.05, 0.1}), d1;
  CK(num(c, c.has("dsh") ? c["dsh"] : ds, "dsh", sv, -1, 0.0));
  CK(num(c, arg(c, "deff", d1, 0), "deff", ev));
  const std::size_t nd = sv.n(), ne = ev.n();
  for (double e : sv.v) {
    if (e > 1) return c.fail("`dsh` must be shares from 0 to 1.");
  }
  double est[6];
  std::vector<double> sens(nd * ne);
  if (clab_jdgiv(zd.data(), dd.data(), yv.p(), wv.n() ? wv.p() : nullptr, static_cast<int>(n), sv.p(), static_cast<int>(nd),
                 ev.p(), static_cast<int>(ne), est, sens.data()) != CLAB_OK) {
    return c.fail("`z` must take both values, each with a positive weight.");
  }
  static const char *const en[] = {"fs", "itt", "wald"};
  static const char *const sn[] = {"complier", "always", "never"};
  double e3[3] = {est[0], est[1], na(est[2])};
  std::vector<double> gs, ge, gc, gt;
  for (std::size_t b = 0; b < ne; ++b) {
    for (std::size_t a = 0; a < nd; ++a) {
      gs.push_back(sv[a]);
      ge.push_back(ev[b]);
      gc.push_back(est[3] + sv[a]);
      gt.push_back(na(sens[a + nd * b]));
    }
  }
  Frame g;
  g.num("dsh", gs);
  g.num("deff", ge);
  g.num("cshare", gc);
  g.num("ceff", gt);
  put(out, "type", vstr("obs"));
  put(out, "est", nvec(3, en, e3));
  put(out, "shares", nvec(3, sn, est + 3));
  put(out, "sens", g.o);
  return true;
}

bool mmshr(Ctx &c, Val &out) {
  Dv sv, mv;
  CK(num(c, c["sent"], "sent", sv));
  const std::size_t n = sv.n();
  CK(num(c, c["min"], "min", mv));
  if (mv.n() != 1 && mv.n() != n) return c.fail("`min` must be one number or one per case, %d.", static_cast<int>(n));
  std::vector<double> mr(n);
  for (std::size_t i = 0; i < n; ++i) mr[i] = mv[mv.n() == 1 ? 0 : i];
  Lab g;
  const bool hg = c.has("grp");
  if (hg) CK(codes(c, c["grp"], "grp", n, g));
  const int ng = hg ? g.g : 0;
  std::vector<double> ov(5 * (static_cast<std::size_t>(ng) + 1));
  clab_mmshr(sv.p(), mr.data(), hg ? g.code.data() : nullptr, static_cast<int>(n), ng, ov.data());
  const int r0 = hg ? 0 : ng, nr = hg ? ng + 1 : 1;
  std::vector<std::string> gl;
  for (int k = 0; k < nr; ++k) gl.push_back(k < nr - 1 ? g.lev[static_cast<std::size_t>(k)] : "all");
  Frame f;
  f.str("group", gl);
  static const char *const cn[] = {"n", "at", "below", "above", "mexcess"};
  for (int q = 0; q < 5; ++q) f.num(cn[q], ov.data() + q * (ng + 1) + r0, static_cast<std::size_t>(nr));
  out = f.o;
  return true;
}

const char *const padjs[] = {"none", "bonferroni", "holm", "hochberg", "BH", "BY"};

bool jmono(Ctx &c, Val &out) {
  const std::size_t n = vlen(c["d"]);
  const int ni = static_cast<int>(n);
  std::vector<double> dd;
  CK(dec01(c, c["d"], "d", n, dd));
  Lab j;
  CK(codes(c, c["judge"], "judge", n, j));
  const int nj = j.g;
  const std::size_t sj = static_cast<std::size_t>(nj);
  Dv yv, wv;
  const bool hy = c.has("y");
  if (hy) CK(num(c, c["y"], "y", yv, static_cast<long>(n)));
  if (c.has("w")) CK(num(c, c["w"], "w", wv, static_cast<long>(n), 0.0));
  double mn = 0.0, al = 0.0;
  CK(binrd(c, "minn", 30, 0.0, Inf, mn));
  int ad = 0;
  CK(pickd(c, "adj", padjs, 6, 2, ad));
  CK(oned(c, "alpha", 0.05, 0.0, 1.0, al));
  std::vector<double> jw(4 * sj + 4 + sj);
  double *cnt = jw.data(), *one1 = cnt + nj, *rate = one1 + nj, *outc = rate + nj, *st = outc + nj;
  if (clab_jrate(j.code.data(), dd.data(), hy ? yv.p() : nullptr, wv.n() ? wv.p() : nullptr, ni, nj, cnt, one1, rate,
                 outc, st) != CLAB_OK) {
    return c.fail("every judge needs cases with a positive total weight.");
  }
  std::vector<int> o(sj);
  for (int k = 0; k < nj; ++k) o[static_cast<std::size_t>(k)] = k;
  std::stable_sort(o.begin(), o.end(), [&](int a, int b) { return rate[a] < rate[b]; });
  std::vector<std::string> jl(sj);
  std::vector<double> fn(sj), sr(sj), sy(sj);
  for (std::size_t k = 0; k < sj; ++k) {
    const int q = o[k];
    jl[k] = j.lev[static_cast<std::size_t>(q)];
    fn[k] = cnt[q];
    sr[k] = rate[q];
    sy[k] = outc[q];
  }
  Frame f;
  f.str("judge", jl);
  f.num("n", fn);
  f.num("rate", sr);
  if (hy) f.num("outcome", sy);
  const double x2 = st[0], df = st[1];
  const double hv[5] = {na(x2), df, std::isnan(x2) ? NA : clab_pchisq(x2, df, 0), st[2], st[3]};
  static const char *const hn[] = {"chisq", "df", "p", "min", "max"};
  out = vobj();
  put(out, "judges", f.o);
  put(out, "homog", nvec(5, hn, hv));
  Val slope, sub;
  if (hy) {
    double l = Inf, h = -Inf;
    for (double e : yv.v) {
      l = std::fmin(l, e);
      h = std::fmax(h, e);
    }
    if (c.has("lo")) {
      double v = 0.0;
      CK(one(c, c["lo"], "lo", -Inf, Inf, v));
      if (v > l) return c.fail("`lo` must be at most the smallest outcome.");
      l = v;
    }
    if (c.has("hi")) {
      double v = 0.0;
      CK(one(c, c["hi"], "hi", -Inf, Inf, v));
      if (v < h) return c.fail("`hi` must be at least the largest outcome.");
      h = v;
    }
    const long long npr = static_cast<long long>(nj) * (nj - 1) / 2;
    if (npr > 50000000LL) return c.fail("`judge` has too many judges for every pair (more than 1e4).");
    const std::size_t m = static_cast<std::size_t>(npr);
    std::vector<std::string> ja, jb;
    for (std::size_t a = 0; a < sj; ++a) {
      for (std::size_t b = a + 1; b < sj; ++b) {
        ja.push_back(jl[a]);
        jb.push_back(jl[b]);
      }
    }
    std::vector<double> dp(m), dy(m), bn(m), lt(m);
    std::vector<int> vi(m);
    clab_jslope(sr.data(), sy.data(), nj, l, h, dp.data(), dy.data(), bn.data(), vi.data(), lt.data());
    for (double &e : lt) e = na(e);
    Frame p;
    p.str("j", ja);
    p.str("k", jb);
    p.num("dp", dp);
    p.num("dy", dy);
    p.num("bound", bn);
    p.lgl("viol", vi.data(), m);
    p.num("late", lt);
    slope = p.o;
  }
  if (c.has("grp")) {
    Lab g;
    CK(codes(c, c["grp"], "grp", n, g));
    if (static_cast<long long>(nj) * g.g > 100000000LL) return c.fail("`judge` and `grp` make more than 1e8 cells.");
    const std::size_t sg = static_cast<std::size_t>(g.g);
    std::vector<double> cw(2 * sj * sg + 2 * sg);
    double *cn = cw.data(), *cr = cn + sj * sg, *sl = cr + sj * sg, *sc = sl + sg;
    std::vector<int> sjv(sg);
    int npair = 0;
    clab_jsub(j.code.data(), g.code.data(), dd.data(), ni, nj, g.g, mn, rate, cn, cr, nullptr, nullptr, nullptr,
              nullptr, sl, sjv.data(), sc, &npair);
    const std::size_t np = static_cast<std::size_t>(npair);
    std::vector<int> iw(3 * np + 1), ow(np + 1);
    std::vector<double> zw(2 * np + 1);
    int *pg = iw.data(), *ps = pg + np, *pl = ps + np;
    double *pz = zw.data(), *pa = pz + np;
    clab_jsub(j.code.data(), g.code.data(), dd.data(), ni, nj, g.g, mn, rate, cn, cr, pg, ps, pl, pz, sl, sjv.data(),
              sc, &npair);
    std::vector<std::string> cg(np), cs(np), cl(np);
    std::vector<double> rs(np), rl(np), gs(np), gl(np), ns(np), nl(np), zz(np), pv(np), pj(np);
    std::vector<int> rv(np), sgf(np);
    for (std::size_t k = 0; k < np; ++k) {
      const std::size_t ch = static_cast<std::size_t>(ps[k]) + sj * static_cast<std::size_t>(pg[k]);
      const std::size_t c2 = static_cast<std::size_t>(pl[k]) + sj * static_cast<std::size_t>(pg[k]);
      cg[k] = g.lev[static_cast<std::size_t>(pg[k])];
      cs[k] = j.lev[static_cast<std::size_t>(ps[k])];
      cl[k] = j.lev[static_cast<std::size_t>(pl[k])];
      rs[k] = rate[ps[k]];
      rl[k] = rate[pl[k]];
      gs[k] = cr[ch];
      gl[k] = cr[c2];
      ns[k] = cn[ch];
      nl[k] = cn[c2];
      zz[k] = pz[k];
      pv[k] = clab_pnorm(-pz[k]);
      rv[k] = gs[k] < gl[k];
    }
    clab_pbadj(pv.data(), npair, ad, pa, ow.data());
    int nrev = 0, nsig = 0;
    for (std::size_t k = 0; k < np; ++k) {
      pj[k] = pa[k];
      sgf[k] = pa[k] < al;
      nrev += rv[k];
      nsig += sgf[k];
    }
    std::vector<long long> spr;
    std::vector<int> sjj;
    for (std::size_t k = 0; k < np; ++k) {
      if (!sgf[k]) continue;
      spr.push_back(static_cast<long long>(ps[k]) * nj + pl[k]);
      sjj.push_back(ps[k]);
      sjj.push_back(pl[k]);
    }
    std::sort(spr.begin(), spr.end());
    std::sort(sjj.begin(), sjj.end());
    const double sm[6] = {static_cast<double>(npair),
                          static_cast<double>(nrev),
                          npair ? static_cast<double>(nrev) / npair : NA,
                          static_cast<double>(nsig),
                          static_cast<double>(std::unique(spr.begin(), spr.end()) - spr.begin()),
                          static_cast<double>(std::unique(sjj.begin(), sjj.end()) - sjj.begin())};
    static const char *const smn[] = {"pairs", "reversed", "share", "signif", "refuted", "judges"};
    Frame pf;
    pf.str("group", cg);
    pf.str("stricter", cs);
    pf.str("lenient", cl);
    pf.num("rates", rs);
    pf.num("ratel", rl);
    pf.num("grates", gs);
    pf.num("gratel", gl);
    pf.num("ns", ns);
    pf.num("nl", nl);
    pf.num("z", zz);
    pf.num("p", pv);
    pf.num("padj", pj);
    pf.lgl("reversed", rv.data(), np);
    pf.lgl("signif", sgf.data(), np);
    std::vector<std::string> bg;
    std::vector<double> bjn, bc, bs;
    for (std::size_t q = 0; q < sg; ++q) {
      if (sjv[q] < 2) continue;
      bg.push_back(g.lev[q]);
      bjn.push_back(sjv[q]);
      bc.push_back(sc[q]);
      bs.push_back(na(sl[q]));
    }
    Frame bf;
    bf.str("group", bg);
    bf.num("judges", bjn);
    bf.num("cases", bc);
    bf.num("slope", bs);
    sub = vobj();
    put(sub, "pairs", pf.o);
    put(sub, "summary", nvec(6, smn, sm));
    put(sub, "slopes", bf.o);
  }
  put(out, "slope", slope);
  put(out, "sub", sub);
  put(out, "alpha", vnum(al));
  return true;
}
