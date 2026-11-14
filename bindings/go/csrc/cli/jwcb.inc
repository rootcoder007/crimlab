// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Inference with few clusters (src/bwcb.cpp): the wild cluster bootstrap,
// randomisation of adoption periods, design power, and the Westfall-Young
// stepdown.

// One two-way fixed-effects fit: the first coefficient and its t.
bool fit1(const std::vector<double> &y, const std::vector<double> &x, int p, const Pan &pn, const std::vector<int> &cl,
          int ncl, double *b, double *t) {
  const int n = pn.n;
  std::vector<double> beta(zz(p)), vc(zz(p) * zz(p)), st(5),
      work(zz(n) * zz(p + 2) + zz(p) * zz(p) + 2 * zz(pn.nu + pn.nt) + zz(n + 2 * p) * zz(p) + zz(p));
  std::vector<int> iwork(zz(n) + zz(std::max(pn.nu, pn.nt)));
  if (clab_didfe(y.data(), x.data(), n, p, pn.u.data(), pn.nu, pn.t.data(), pn.nt, 1, nullptr, cl.data(), ncl, nullptr,
                 0, 1e-10, 10000, beta.data(), vc.data(), st.data(), work.data(), iwork.data()) != CLAB_OK) {
    return false;
  }
  *b = beta[0];
  *t = beta[0] / std::sqrt(vc[0]);
  return true;
}

bool wcbx(Ctx &c, Val &out) {
  static const char *const mo[] = {"wild", "perm", "power"};
  static const char *const wo[] = {"rad", "webb"};
  static const char *const to[] = {"count", "strict"};
  static const char *const so[] = {"t", "est"};
  int mt = 0, wt = 0, ti = 0, sk = 0, nr = 0;
  CK(pickd(c, "meth", mo, 3, 0, mt));
  CK(pickd(c, "wtype", wo, 2, 0, wt));
  CK(pickd(c, "ties", to, 2, 0, ti));
  const int tie = ti == 0 ? 1 : 0;
  CK(pickd(c, "stat", so, 2, 0, sk));
  if (c.has("reps")) CK(whole(c, c["reps"], "reps", 9, nr));
  else nr = mt == 2 ? 200 : 9999;
  double lv = 0.0;
  CK(oned(c, "level", 0.95, 0.0, 1.0, lv));
  const Val &y = c["y"];
  std::vector<Val> ys;
  std::vector<std::string> ynm;
  int dim[2] = {-1, -1};
  if (y.t == Val::OBJ && !y.a.empty() && y.a[0].t == Val::ARR) {
    ys = y.a;
    ynm = y.k;
  } else if (y.t == Val::ARR && !y.a.empty() && y.a[0].t == Val::ARR && shape(y, 2, dim, 0) && dim[1] > 1) {
    for (int j = 0; j < dim[1]; ++j) {
      Val v = varr();
      for (const Val &r : y.a) v.a.push_back(r.a[zz(j)]);
      ys.push_back(v);
    }
  } else {
    ys.push_back(y);
  }
  if (ys.size() > 1 && mt != 1) return c.fail("several outcomes are taken by `meth = \"perm\"` only.");
  Pan pn;
  CK(panel(c, -HUGE_VAL, pn, &ys[0]));
  const int n = pn.n, k = static_cast<int>(ys.size());
  std::vector<std::vector<double>> yv(zz(k));
  for (int j = 0; j < k; ++j) {
    Dv v;
    CK(num(c, ys[zz(j)], "y", v, n));
    yv[zz(j)] = v.v;
  }
  std::vector<int> cl(pn.u);
  int ncl = pn.nu;
  if (c.has("clus")) {
    Lab l;
    CK(codes(c, c["clus"], "clus", zz(n), l));
    cl = l.code;
    ncl = l.g;
    if (ncl < 2) return c.fail("`clus` has one cluster: the inference needs two or more.");
  }
  std::vector<Term> x;
  CK(adjcols(c, c["adj"], zz(n), x, "row of the panel"));
  const int pa = static_cast<int>(x.size());
  const int p = pa + 1;
  std::vector<double> xa(zz(n) * zz(p), 0.0);
  for (int j = 0; j < pa; ++j) std::copy(x[zz(j)].x.begin(), x[zz(j)].x.end(), xa.begin() + static_cast<std::ptrdiff_t>(j + 1) * n);
  Val res;
  int nd = 0;
  bool en = false;
  double mde = NA, size = NA;
  if (mt == 0 || mt == 1) {
    if (!c.has("trt")) return c.fail("`trt` is needed for the wild bootstrap and the permutations.");
    std::vector<int> dv;
    CK(bin01(c, c["trt"], "trt", zz(n), dv));
    for (std::size_t i = 0; i < zz(n); ++i) xa[i] = dv[i];
    if (mt == 0) {
      pev::wctx w;
      w.n = n;
      w.p = p;
      w.x = xa;
      w.y = yv[0];
      w.g = cl;
      w.ng = ncl;
      w.u = pn.u;
      w.t = pn.t;
      w.nu = pn.nu;
      w.nt = pn.nt;
      std::vector<double> wk(2 * zz(pn.nu + pn.nt));
      pev::within(w.x.data(), n, p, pn.u.data(), pn.nu, pn.t.data(), pn.nt, 1, nullptr, 1e-12, 10000, wk.data());
      pev::within(w.y.data(), n, 1, pn.u.data(), pn.nu, pn.t.data(), pn.nt, 1, nullptr, 1e-12, 10000, wk.data());
      auto nested = [&](const std::vector<int> &f, int nf) {
        std::vector<int> seen(zz(nf), -1);
        for (std::size_t i = 0; i < zz(n); ++i) {
          if (seen[zz(f[i])] < 0) seen[zz(f[i])] = cl[i];
          if (seen[zz(f[i])] != cl[i]) return false;
        }
        return true;
      };
      double kk = p + 1.0;
      if (!nested(pn.u, pn.nu)) kk += pn.nu - 1.0;
      if (!nested(pn.t, pn.nt)) kk += pn.nt - 1.0;
      w.adj = static_cast<double>(ncl) / (ncl - 1) * (n - 1.0) / (n - kk);
      if (!pev::wsetup(w)) return c.fail("%s", absorbed);
      double bh = 0.0;
      const double th = pev::tstat(w, w.y.data(), 0.0, &bh);
      en = wt == 0 && ncl <= 30 && std::ldexp(1.0, ncl) <= nr;
      std::vector<double> v;
      if (en) {
        nd = 1 << ncl;
        v.resize(zz(nd) * zz(ncl));
        for (int d = 0; d < nd; ++d) {
          for (int g = 0; g < ncl; ++g) v[zz(d) * zz(ncl) + zz(g)] = (d >> g) & 1 ? 1.0 : -1.0;
        }
      } else {
        static const double webb[6] = {-1.224744871391589, -1.0, -0.7071067811865476, 0.7071067811865476, 1.0,
                                       1.224744871391589};
        unsigned int sd = 0;
        CK(seed(c, sd));
        kutil::rng g(sd);
        nd = nr;
        v.resize(zz(nd) * zz(ncl));
        for (double &e : v) e = wt == 0 ? (g.below(2) ? 1.0 : -1.0) : webb[g.below(6)];
      }
      const double pv = pev::wpval(w, 0.0, v, nd, tie);
      const double al = 1.0 - lv, s0 = std::fabs(bh / th);
      const bool open = tie && en && 2.0 / nd > al;
      const double lo = open ? -Inf : pev::wend(w, bh, s0, al, -1, v, nd, tie);
      const double hi = open ? Inf : pev::wend(w, bh, s0, al, 1, v, nd, tie);
      Frame r;
      r.num("est", &bh, 1);
      r.num("stat", &th, 1);
      r.num("p", &pv, 1);
      r.num("lo", &lo, 1);
      r.num("hi", &hi, 1);
      res = r.o;
    } else {
      std::vector<int> g;
      CK(firsts(c, pn, dv.data(), g));
      std::vector<int> ga(zz(ncl), -1);
      for (std::size_t i = 0; i < zz(n); ++i) {
        const int a = g[zz(pn.u[i])], q = cl[i];
        if (ga[zz(q)] >= 0 && ga[zz(q)] != a) return c.fail("the units of a cluster must be treated from the same period.");
        ga[zz(q)] = a;
      }
      std::vector<double> ob(zz(k)), oe(zz(k));
      for (int j = 0; j < k; ++j) {
        double b = 0.0, t = 0.0;
        if (!fit1(yv[zz(j)], xa, p, pn, cl, ncl, &b, &t)) return c.fail("%s", absorbed);
        oe[zz(j)] = b;
        ob[zz(j)] = sk == 0 ? t : b;
      }
      std::vector<int> s = ga;
      std::sort(s.begin(), s.end());
      double lc = std::lgamma(ncl + 1.0);
      for (std::size_t a = 0; a < s.size();) {
        std::size_t b = a;
        while (b < s.size() && s[b] == s[a]) ++b;
        lc -= std::lgamma(static_cast<double>(b - a) + 1.0);
        a = b;
      }
      en = std::exp(lc) <= nr + 0.5;
      std::vector<std::vector<int>> as;
      if (en) {
        do as.push_back(s);
        while (std::next_permutation(s.begin(), s.end()));
      } else {
        unsigned int sd = 0;
        CK(seed(c, sd));
        kutil::rng rg(sd);
        for (int r = 0; r < nr; ++r) {
          std::vector<int> q = ga;
          for (int a = ncl - 1; a > 0; --a) std::swap(q[zz(a)], q[zz(rg.below(a + 1))]);
          as.push_back(q);
        }
      }
      nd = static_cast<int>(as.size());
      std::vector<double> nv0(zz(nd) * zz(k), NA), xp = xa;
      std::vector<double> cnt(zz(k), 0.0), ok(zz(k), 0.0);
      for (int d = 0; d < nd; ++d) {
        for (std::size_t i = 0; i < zz(n); ++i) xp[i] = pn.t[i] >= as[zz(d)][zz(cl[i])] ? 1.0 : 0.0;
        for (int j = 0; j < k; ++j) {
          double b = 0.0, t = 0.0;
          if (!fit1(yv[zz(j)], xp, p, pn, cl, ncl, &b, &t)) continue;
          const double sv = sk == 0 ? t : b;
          nv0[zz(j) * zz(nd) + zz(d)] = sv;
          ok[zz(j)] += 1.0;
          cnt[zz(j)] += tie ? std::fabs(sv) >= std::fabs(ob[zz(j)]) - 1e-12 : std::fabs(sv) > std::fabs(ob[zz(j)]) + 1e-12;
        }
      }
      std::vector<double> full, pr(zz(k)), pst(zz(k), NA);
      int nf = 0;
      for (int d = 0; d < nd; ++d) {
        bool all = true;
        for (int j = 0; j < k; ++j) all = all && !std::isnan(nv0[zz(j) * zz(nd) + zz(d)]);
        nf += all;
      }
      full.resize(zz(nf) * zz(k));
      for (int d = 0, q = 0; d < nd; ++d) {
        bool all = true;
        for (int j = 0; j < k; ++j) all = all && !std::isnan(nv0[zz(j) * zz(nd) + zz(d)]);
        if (!all) continue;
        for (int j = 0; j < k; ++j) full[zz(j) * zz(nf) + zz(q)] = nv0[zz(j) * zz(nd) + zz(d)];
        ++q;
      }
      if (k > 1 && nf > 0) pev::wystep(ob.data(), full.data(), nf, k, pr.data(), pst.data());
      std::vector<std::string> ro;
      std::vector<double> rp, rs;
      for (int j = 0; j < k; ++j) {
        ro.push_back(zz(j) < ynm.size() ? ynm[zz(j)] : "y" + std::to_string(j + 1));
        rp.push_back(en ? cnt[zz(j)] / ok[zz(j)] : (cnt[zz(j)] + 1.0) / (ok[zz(j)] + 1.0));
        rs.push_back(k > 1 ? pst[zz(j)] : rp.back());
      }
      Frame r;
      r.str("outcome", ro);
      r.num("est", oe);
      r.num("stat", ob);
      r.num("p", rp);
      r.num("pstep", rs);
      res = r.o;
    }
  } else {
    int nt1 = 0;
    CK(wholed(c, "ntrt", 1, 1, nt1));
    if (nt1 >= ncl) return c.fail("`ntrt` must leave some clusters untreated: below %d.", ncl);
    if (!c.has("adopt")) return c.fail("`adopt` must give the adoption periods to draw from.");
    Dv av;
    CK(num(c, c["adopt"], "adopt", av));
    std::vector<int> ac;
    for (double v : av.v) {
      const auto it = std::find(pn.tv.begin(), pn.tv.end(), v);
      if (it == pn.tv.end() || it == pn.tv.begin()) return c.fail("`adopt` must be periods of `time` after the first.");
      ac.push_back(static_cast<int>(it - pn.tv.begin()));
    }
    double al = 0.0, pw = 0.0;
    CK(oned(c, "alpha", 0.05, 0.0, 1.0, al));
    CK(oned(c, "power", 0.8, 0.0, 1.0, pw));
    std::vector<double> ev;
    if (!c.has("eff")) {
      double m = 0, s = 0;
      for (double v : yv[0]) m += v / n;
      for (double v : yv[0]) s += (v - m) * (v - m) / (n - 1);
      const double sd = std::sqrt(std::max(0.0, s));
      for (int q = 0; q < 5; ++q) ev.push_back(0.5 * q * sd);
    } else {
      Dv e;
      CK(num(c, c["eff"], "eff", e));
      ev = e.v;
      std::sort(ev.begin(), ev.end());
    }
    const int ne = static_cast<int>(ev.size());
    std::vector<double> rej(zz(ne), 0.0), tot(zz(ne), 0.0), yy(zz(n)), xp = xa;
    unsigned int sd = 0;
    CK(seed(c, sd));
    kutil::rng rg(sd);
    for (int r = 0; r < nr; ++r) {
      std::vector<int> q(zz(ncl)), ga(zz(ncl), pn.nt);
      for (int a = 0; a < ncl; ++a) q[zz(a)] = a;
      for (int a = 0; a < nt1; ++a) {
        std::swap(q[zz(a)], q[zz(a + rg.below(ncl - a))]);
        ga[zz(q[zz(a)])] = ac[zz(rg.below(static_cast<int>(ac.size())))];
      }
      for (std::size_t i = 0; i < zz(n); ++i) xp[i] = pn.t[i] >= ga[zz(cl[i])] ? 1.0 : 0.0;
      for (int e = 0; e < ne; ++e) {
        for (std::size_t i = 0; i < zz(n); ++i) yy[i] = yv[0][i] + ev[zz(e)] * xp[i];
        double b = 0.0, t = 0.0;
        if (!fit1(yy, xp, p, pn, cl, ncl, &b, &t)) continue;
        tot[zz(e)] += 1.0;
        rej[zz(e)] += 2 * clab_pt(-std::fabs(t), ncl - 1.0, 1) < al ? 1.0 : 0.0;
      }
    }
    std::vector<double> rp(zz(ne));
    for (int e = 0; e < ne; ++e) {
      rp[zz(e)] = tot[zz(e)] > 0 ? rej[zz(e)] / tot[zz(e)] : NA;
      if (ev[zz(e)] == 0) size = rp[zz(e)];
    }
    for (int e = 1; e < ne && std::isnan(mde); ++e) {
      if (rp[zz(e - 1)] < pw && rp[zz(e)] >= pw) {
        mde = ev[zz(e - 1)] + (pw - rp[zz(e - 1)]) / (rp[zz(e)] - rp[zz(e - 1)]) * (ev[zz(e)] - ev[zz(e - 1)]);
      }
    }
    if (ne > 0 && rp[0] >= pw) mde = ev[0];
    nd = nr;
    Frame r;
    r.num("eff", ev);
    r.num("power", rp);
    res = r.o;
  }
  out = vobj();
  put(out, "res", res);
  put(out, "meth", vstr(mo[mt]));
  put(out, "ncl", vnum(ncl));
  put(out, "draws", vnum(nd));
  put(out, "enum", vbool(en));
  put(out, "mde", vnum(mde));
  put(out, "size", vnum(size));
  put(out, "level", vnum(lv));
  return true;
}

bool wystep(Ctx &c, Val &out) {
  CK(isvec(c, c["obs"], "obs"));
  Dv ob;
  CK(num(c, c["obs"], "obs", ob));
  const int k = static_cast<int>(ob.n());
  const Val &nl = c["null"];
  int dim[2] = {-1, -1};
  if (!(nl.t == Val::ARR && !nl.a.empty() && nl.a[0].t == Val::ARR && shape(nl, 2, dim, 0)) || dim[1] != k) {
    return c.fail("`null` must be a numeric matrix with one column per statistic of `obs`, %d.", k);
  }
  const int nr = dim[0];
  std::vector<double> nv(zz(nr) * zz(k));
  bool bad = false;
  fill(nl, 2, dim, 0, 0, 1, nv, bad, true);
  std::vector<double> full;
  int nd = 0;
  for (int d = 0; d < nr; ++d) {
    bool ok = true;
    for (int j = 0; j < k; ++j) ok = ok && std::isfinite(nv[zz(j) * zz(nr) + zz(d)]);
    nd += ok;
  }
  if (nd == 0) return c.fail("`null` has no draw with every statistic finite.");
  full.resize(zz(nd) * zz(k));
  for (int d = 0, q = 0; d < nr; ++d) {
    bool ok = true;
    for (int j = 0; j < k; ++j) ok = ok && std::isfinite(nv[zz(j) * zz(nr) + zz(d)]);
    if (!ok) continue;
    for (int j = 0; j < k; ++j) full[zz(j) * zz(nd) + zz(q)] = nv[zz(j) * zz(nr) + zz(d)];
    ++q;
  }
  std::vector<double> cp(zz(k)), ca(zz(k));
  pev::wystep(ob.p(), full.data(), nd, k, cp.data(), ca.data());
  Frame f;
  f.num("stat", ob.v);
  f.num("p", cp);
  f.num("padj", ca);
  out = f.o;
  return true;
}
