// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Counts, rates and tables (src/bstat.cpp): count regression, rates by
// group, direct and indirect standardisation, and contingency tables.

// level: the level given (default 0.95), checked as zq() checks it.
bool lvl(Ctx &c, double &z, double &lv) {
  CK(zconf(c, "level", 0.95, z));
  lv = val1d(c, "level", 0.95);
  return true;
}

// One column of a count regression's design: values and term label.
struct Term {
  std::vector<double> x;
  std::string nm;
};

// adjcols(): an object of columns, a numeric matrix (its columns) or one
// vector; labels become indicators of every level but the first.
bool adjcols(Ctx &c, const Val &adj, std::size_t n, std::vector<Term> &out, const char *per = "count") {
  if (adj.t == Val::NUL) return true;
  std::vector<std::pair<std::string, const Val *>> cs;
  int dim[2] = {-1, -1};
  if (adj.t == Val::OBJ) {
    for (std::size_t j = 0; j < adj.a.size(); ++j) {
      cs.emplace_back(adj.k[j].empty() ? "adj" + std::to_string(j + 1) : adj.k[j], &adj.a[j]);
    }
  } else if (adj.t == Val::ARR && !adj.a.empty() && adj.a[0].t == Val::ARR && shape(adj, 2, dim, 0)) {
    const int nc = dim[1];
    if (static_cast<std::size_t>(dim[0]) != n) return c.fail("`adj` must have one row per %s, %d.", per, static_cast<int>(n));
    std::vector<double> v(n * static_cast<std::size_t>(nc));
    bool bad = false;
    fill(adj, 2, dim, 0, 0, 1, v, bad, true);
    for (int j = 0; j < nc; ++j) {
      Term t;
      t.nm = "adj" + std::to_string(j + 1);
      t.x.assign(v.begin() + static_cast<std::ptrdiff_t>(j * n), v.begin() + static_cast<std::ptrdiff_t>((j + 1) * n));
      for (double e : t.x) {
        if (!std::isfinite(e)) return c.fail("`adj` must hold finite numbers with no missing values.");
      }
      out.push_back(t);
    }
    return true;
  } else {
    cs.emplace_back("adj", &adj);
  }
  for (auto &p : cs) {
    const Val &x = *p.second;
    const std::string an = "adj " + p.first;
    const bool vec = x.t != Val::OBJ && !(x.t == Val::ARR && !x.a.empty() && x.a[0].t == Val::ARR);
    if (!vec || vlen(x) != n) {
      return c.fail("`adj` column %s must be a vector with one value per %s, %d.", p.first.c_str(), per,
                    static_cast<int>(n));
    }
    if (isnum(x)) {
      Dv v;
      CK(num(c, x, an.c_str(), v, static_cast<long>(n)));
      Term t;
      t.nm = p.first;
      t.x = v.v;
      out.push_back(t);
      continue;
    }
    Lab l;
    CK(codes(c, x, an.c_str(), n, l));
    for (int k = 1; k < l.g; ++k) {
      Term t;
      t.nm = p.first + "=" + l.lev[static_cast<std::size_t>(k)];
      t.x.resize(n);
      for (std::size_t i = 0; i < n; ++i) t.x[i] = l.code[i] == k ? 1.0 : 0.0;
      out.push_back(t);
    }
  }
  return true;
}

bool cntreg(Ctx &c, Val &out) {
  Args a(c, {"y", "grp", "expo", "clus"});
  CK(a.byname(c));
  Dv yv;
  CK(num(c, a[0], "y", yv, -1, 0.0));
  const std::size_t sn = yv.n();
  const int n = static_cast<int>(sn);
  double tot = 0.0;
  for (double e : yv.v) tot += e;
  if (!(tot > 0)) return c.fail("`y` has no events: every count is 0.");
  if (a[1].t == Val::NUL) return c.fail("`grp` must give the group of each count.");
  Lab g;
  CK(codes(c, a[1], "grp", sn, g));
  if (g.g < 2) return c.fail("`grp` has one group: a rate ratio needs two.");
  int r0 = 1;
  CK(ref(c, c["ref"], g.lev, r0));
  --r0;
  double z = 0.0, lv = 0.0;
  CK(lvl(c, z, lv));
  bool dogee = false, adjc = false;
  CK(flagd(c, "gee", false, dogee));
  CK(flagd(c, "cadj", true, adjc));
  std::vector<double> off;
  if (a[2].t != Val::NUL) {
    Dv e;
    CK(num(c, a[2], "expo", e, n, 0.0, true));
    off.resize(sn);
    for (std::size_t i = 0; i < sn; ++i) off[i] = std::log(e[i]);
  }
  int ncl = 0;
  Lab cl;
  if (a[3].t != Val::NUL) {
    CK(codes(c, a[3], "clus", sn, cl));
    if (cl.g < 2) return c.fail("`clus` has one cluster: a cluster-robust variance needs two or more.");
    ncl = cl.g;
  }
  const int *clp = ncl ? cl.code.data() : nullptr;
  if (dogee && !clp) return c.fail("`gee = TRUE` needs `clus`.");
  std::vector<Term> ts(1);
  ts[0].nm = "intercept";
  ts[0].x.assign(sn, 1.0);
  std::vector<int> gt;
  for (int k = 0; k < g.g; ++k) {
    if (k == r0) continue;
    Term t;
    t.nm = g.lev[static_cast<std::size_t>(k)];
    t.x.resize(sn);
    for (std::size_t i = 0; i < sn; ++i) t.x[i] = g.code[i] == k ? 1.0 : 0.0;
    ts.push_back(t);
    gt.push_back(k);
  }
  CK(adjcols(c, c["adj"], sn, ts));
  const int p = static_cast<int>(ts.size());
  const std::size_t sp = static_cast<std::size_t>(p);
  std::vector<double> x(sn * sp);
  for (std::size_t j = 0; j < sp; ++j) std::copy(ts[j].x.begin(), ts[j].x.end(), x.begin() + static_cast<std::ptrdiff_t>(j * sn));
  const std::size_t sc = static_cast<std::size_t>(ncl);
  std::vector<double> work(sn * sp + 6 * sn + sp + sp * sp + 2 * sc * sp + 2 * sn);
  std::vector<int> iwork(sp + sc + 1 + sn);
  std::vector<double> beta(3 * sp, NaN), vb(2 * sp * sp, NaN), vc(3 * sp * sp, NaN), st(3 * 8, NaN);
  std::vector<int> alias(2 * sp, 1);
  const double *of = off.empty() ? nullptr : off.data();
  bool ok[3] = {false, false, false};
  for (int f = 0; f < 2; ++f) {
    const std::size_t fo = static_cast<std::size_t>(f);
    const int rc = clab_cntreg(yv.p(), x.data(), n, p, of, clp, ncl, f, 25, 1e-8, &beta[fo * sp], &vb[fo * sp * sp],
                               &vc[fo * sp * sp], &st[fo * 8], &alias[fo * sp], work.data(), iwork.data());
    if (f == 0 && rc != CLAB_OK) {
      return c.fail("the fit failed: the deviance or a coefficient left the finite range (a group or cell with no "
                    "counts can do this).");
    }
    if (rc == CLAB_OK) ok[f] = true;
    else c.warn.push_back("cntreg: the negative binomial fit failed (theta has no finite estimate); its rows are NA.");
  }
  if (dogee) {
    std::copy(beta.begin(), beta.begin() + static_cast<std::ptrdiff_t>(sp), beta.begin() + static_cast<std::ptrdiff_t>(2 * sp));
    if (clab_geeex(yv.p(), x.data(), n, p, of, clp, ncl, alias.data(), 100, 1e-10, &beta[2 * sp], &vc[2 * sp * sp],
                   &st[16], work.data(), iwork.data()) == CLAB_OK) {
      ok[2] = true;
    } else {
      c.warn.push_back("cntreg: the GEE fit failed (the exchangeable working correlation is not positive definite, "
                       "or the scoring diverged); its rows are NA.");
    }
  }
  const double fc = adjc && ncl > 1 ? static_cast<double>(ncl) / (ncl - 1) : 1.0;
  struct Row {
    const char *lab;
    int f = 0, v = 0;
  };
  std::vector<Row> rows = {{"pois", 0, 0}, {"nb", 1, 0}};
  if (clp) {
    rows.push_back({"poiscl", 0, 1});
    rows.push_back({"nbcl", 1, 1});
  }
  if (dogee) rows.push_back({"gee", 2, 1});
  std::vector<double> gsum(static_cast<std::size_t>(g.g), 0.0);
  for (std::size_t i = 0; i < sn; ++i) gsum[static_cast<std::size_t>(g.code[i])] += yv[i];
  auto est = [&](const Row &w, std::size_t j, double &b, double &s) {
    const std::size_t fo = static_cast<std::size_t>(w.f);
    b = ok[w.f] ? beta[fo * sp + j] : NaN;
    const double v = !ok[w.f] ? NaN
                     : w.v == 0 ? vb[fo * sp * sp + j * sp + j]
                                : vc[fo * sp * sp + j * sp + j] * (w.f == 2 ? 1.0 : fc);
    s = std::sqrt(v);
  };
  std::vector<std::string> cfit, cterm, ifit, iterm;
  std::vector<double> ce, cs, cz, cp, ii, il, ih, ip;
  std::vector<int> is;
  for (const Row &w : rows) {
    for (std::size_t j = 0; j < sp; ++j) {
      double b = 0.0, s = 0.0;
      est(w, j, b, s);
      cfit.push_back(w.lab);
      cterm.push_back(ts[j].nm);
      ce.push_back(na(b));
      cs.push_back(na(s));
      cz.push_back(na(b / s));
      cp.push_back(std::isnan(b / s) ? NA : 2 * clab_pnorm(-std::fabs(b / s)));
    }
    for (std::size_t k = 0; k < gt.size(); ++k) {
      double b = 0.0, s = 0.0;
      est(w, k + 1, b, s);
      ifit.push_back(w.lab);
      iterm.push_back(ts[k + 1].nm);
      const bool fin = !std::isnan(b) && !std::isnan(s);
      ii.push_back(na(std::exp(b)));
      il.push_back(fin ? std::exp(b - z * s) : NA);
      ih.push_back(fin ? std::exp(b + z * s) : NA);
      ip.push_back(fin ? 2 * clab_pnorm(-std::fabs(b / s)) : NA);
      is.push_back(gsum[static_cast<std::size_t>(gt[k])] == 0 || gsum[static_cast<std::size_t>(r0)] == 0 ? 1 : 0);
    }
  }
  Frame fi, fco, ff;
  fi.str("fit", ifit);
  fi.str("term", iterm);
  fi.num("irr", ii);
  fi.num("lo", il);
  fi.num("hi", ih);
  fi.num("p", ip);
  fi.lgl("sep", is.data(), is.size());
  fco.str("fit", cfit);
  fco.str("term", cterm);
  fco.num("est", ce);
  fco.num("se", cs);
  fco.num("z", cz);
  fco.num("p", cp);
  const int nf = dogee ? 3 : 2;
  static const char *const fl[] = {"pois", "nb", "gee"};
  std::vector<std::string> fn;
  std::vector<double> fll, fa, ft, fst, fal, fit, frk;
  std::vector<int> fcv;
  for (int f = 0; f < nf; ++f) {
    const double *s = &st[static_cast<std::size_t>(f) * 8];
    fn.push_back(fl[f]);
    auto nv = [&](double v) { return ok[f] && !std::isnan(v) ? v : NA; };
    if (f < 2) {
      fll.push_back(nv(s[0]));
      fa.push_back(nv(s[1]));
      ft.push_back(nv(s[2]));
      fst.push_back(nv(s[3]));
      fal.push_back(NA);
      fit.push_back(nv(s[4]));
      fcv.push_back(ok[f] ? (s[5] != 0 ? 1 : 0) : -1);
      frk.push_back(nv(s[6]));
    } else {
      fll.push_back(NA);
      fa.push_back(NA);
      ft.push_back(NA);
      fst.push_back(NA);
      fal.push_back(nv(s[0]));
      fit.push_back(nv(s[2]));
      fcv.push_back(ok[f] ? (s[3] != 0 ? 1 : 0) : -1);
      frk.push_back(nv(st[6]));
    }
  }
  ff.str("fit", fn);
  ff.num("loglik", fll);
  ff.num("aic", fa);
  ff.num("theta", ft);
  ff.num("setheta", fst);
  ff.num("alpha", fal);
  ff.num("iter", fit);
  ff.lgl("conv", fcv.data(), fcv.size());
  ff.num("rank", frk);
  out = vobj();
  put(out, "irr", fi.o);
  put(out, "coef", fco.o);
  put(out, "fit", ff.o);
  put(out, "ref", vstr(g.lev[static_cast<std::size_t>(r0)]));
  put(out, "n", vnum(n));
  put(out, "ncl", vnum(ncl));
  put(out, "level", vnum(lv));
  return true;
}

// Clopper-Pearson limits of x of n, and the exact Poisson limits of a
// count x, at tail a / 2 each side.
double cplo(double x, double n, double a) { return x > 0 ? clab_qbeta(a / 2, x, n - x + 1, 1) : 0.0; }
double cphi(double x, double n, double a) { return x < n ? clab_qbeta(1 - a / 2, x + 1, n - x, 1) : 1.0; }
double gmlo(double x, double a) { return x > 0 ? clab_qgamma(a / 2, x, 1) : 0.0; }
double gmhi(double x, double a) { return clab_qgamma(1 - a / 2, x + 1, 1); }

// groups(): sums of x and n by grp, or each row its own group labelled by
// the names of x or 1, 2, ...
struct Grpd {
  std::vector<double> x, n;
  std::vector<int> code;
  std::vector<std::string> lev;
  int g = 0;
};
bool groups(Ctx &c, const Val &grp, const Dv &xv, const double *nv, Grpd &r) {
  const std::size_t m = xv.n();
  if (grp.t == Val::NUL) {
    r.lev = nmsof(xv.nm, m);
    r.g = static_cast<int>(m);
    r.code.resize(m);
    for (std::size_t i = 0; i < m; ++i) r.code[i] = static_cast<int>(i);
  } else {
    Lab l;
    CK(codes(c, grp, "grp", m, l));
    r.lev = l.lev;
    r.g = l.g;
    r.code = l.code;
  }
  r.x.assign(static_cast<std::size_t>(r.g), 0.0);
  r.n.assign(static_cast<std::size_t>(r.g), 0.0);
  for (std::size_t i = 0; i < m; ++i) {
    const std::size_t k = static_cast<std::size_t>(r.code[i]);
    r.x[k] += xv[i];
    r.n[k] += nv ? nv[i] : 1.0;
  }
  return true;
}

// perstr(): one value per stratum, named ones matched by label.
bool perstr(Ctx &c, const Val &v, const char *nm, const std::vector<std::string> &lev, std::vector<double> &out) {
  Dv d;
  CK(num(c, v, nm, d, -1, 0.0));
  const std::size_t ns = lev.size();
  out.assign(ns, 0.0);
  if (d.nm.empty()) {
    if (d.n() != ns) return c.fail("`%s` must have one value per stratum, %d.", nm, static_cast<int>(ns));
    out = d.v;
    return true;
  }
  for (std::size_t s = 0; s < ns; ++s) {
    int hit = -1;
    for (std::size_t j = 0; j < d.n(); ++j) {
      if (d.nm[j] == lev[s]) {
        if (hit >= 0) return c.fail("`%s` names stratum %s twice.", nm, lev[s].c_str());
        hit = static_cast<int>(j);
      }
    }
    if (hit < 0) return c.fail("`%s` has no value for stratum %s.", nm, lev[s].c_str());
    out[s] = d[static_cast<std::size_t>(hit)];
  }
  return true;
}

const char *const rtype[] = {"rate", "prop"};

bool grprr(Ctx &c, Val &out) {
  Args ar(c, {"x", "n", "grp"});
  CK(ar.byname(c));
  Dv xv;
  CK(num(c, ar[0], "x", xv, -1, 0.0));
  const std::size_t m = xv.n();
  int prop = 0;
  CK(pickd(c, "type", rtype, 2, 0, prop));
  Dv nn;
  const double *nv = nullptr;
  if (ar[1].t != Val::NUL) {
    CK(num(c, ar[1], "n", nn, static_cast<long>(m), 0.0, true));
    nv = nn.p();
  }
  double pr = 0.0, z = 0.0, lv = 0.0;
  CK(oned(c, "per", 1, 0.0, Inf, pr));
  CK(lvl(c, z, lv));
  const double a = 1 - lv;
  Grpd g;
  CK(groups(c, ar[2], xv, nv, g));
  for (int k = 0; k < g.g; ++k) {
    if (prop && g.x[static_cast<std::size_t>(k)] > g.n[static_cast<std::size_t>(k)]) {
      return c.fail("group %s has more events than trials: `x` must not exceed `n` for proportions.",
                    g.lev[static_cast<std::size_t>(k)].c_str());
    }
  }
  if (!prop && !nv) return c.fail("`n` must give the exposure of each count for `type = \"rate\"`.");
  int r0 = 1;
  CK(ref(c, c["ref"], g.lev, r0));
  --r0;
  const int k = g.g;
  const std::size_t sk = static_cast<std::size_t>(k);
  std::vector<double> o(10 * sk);
  double *rt = &o[0], *lo = &o[sk], *hi = &o[2 * sk], *rp = &o[3 * sk], *ra = &o[4 * sk], *rl = &o[5 * sk],
         *rh = &o[6 * sk], *df = &o[7 * sk], *dl = &o[8 * sk], *dh = &o[9 * sk];
  if (clab_grprr(g.x.data(), g.n.data(), k, r0, prop, z, rt, lo, hi, rp, ra, rl, rh, df, dl, dh) != CLAB_OK) {
    return c.fail("the kernel refused the input.");
  }
  std::vector<double> cr(sk), cl(sk), ch(sk), el(sk), eh(sk), cp(sk);
  for (std::size_t i = 0; i < sk; ++i) {
    const double xi = g.x[i], ni = g.n[i];
    cr[i] = pr * rt[i];
    cl[i] = na(pr * lo[i]);
    ch[i] = na(pr * hi[i]);
    el[i] = pr * (prop ? cplo(xi, ni, a) : gmlo(xi, a) / ni);
    eh[i] = pr * (prop ? cphi(xi, ni, a) : gmhi(xi, a) / ni);
    cp[i] = na(rp[i]);
  }
  Frame rs;
  rs.str("group", g.lev);
  rs.num("x", g.x);
  rs.num("n", g.n);
  rs.num("rate", cr);
  rs.num("lo", cl);
  rs.num("hi", ch);
  rs.num("elo", el);
  rs.num("ehi", eh);
  rs.num("rep", cp);
  std::vector<std::string> cg;
  std::vector<double> mr, ml, mh, me, mf, md, mdl, mdh;
  const double x0 = g.x[static_cast<std::size_t>(r0)], n0 = g.n[static_cast<std::size_t>(r0)];
  for (int i = 0; i < k; ++i) {
    if (i == r0) continue;
    const std::size_t si = static_cast<std::size_t>(i);
    const double x1 = g.x[si], n1 = g.n[si];
    cg.push_back(g.lev[si]);
    mr.push_back(na(ra[i]));
    ml.push_back(na(rl[i]));
    mh.push_back(na(rh[i]));
    if (!prop && x1 + x0 > 0) {
      const double pl = cplo(x1, x1 + x0, a), pu = cphi(x1, x1 + x0, a);
      me.push_back(pl * n0 / ((1 - pl) * n1));
      mf.push_back(pu < 1 ? pu * n0 / ((1 - pu) * n1) : Inf);
    } else {
      me.push_back(NA);
      mf.push_back(NA);
    }
    md.push_back(pr * df[i]);
    mdl.push_back(na(pr * dl[i]));
    mdh.push_back(na(pr * dh[i]));
  }
  Frame cm;
  cm.str("group", cg);
  cm.num("ratio", mr);
  cm.num("lo", ml);
  cm.num("hi", mh);
  cm.num("elo", me);
  cm.num("ehi", mf);
  cm.num("diff", md);
  cm.num("dlo", mdl);
  cm.num("dhi", mdh);
  out = vobj();
  put(out, "rates", rs.o);
  put(out, "cmp", cm.o);
  put(out, "ref", vstr(g.lev[static_cast<std::size_t>(r0)]));
  put(out, "type", vstr(rtype[prop]));
  put(out, "per", vnum(pr));
  put(out, "level", vnum(lv));
  return true;
}

bool dastd(Ctx &c, Val &out) {
  Args ar(c, {"x", "n", "strat", "grp"});
  CK(ar.byname(c));
  Dv xv, nv;
  CK(num(c, ar[0], "x", xv, -1, 0.0));
  const std::size_t m = xv.n();
  CK(num(c, ar[1], "n", nv, static_cast<long>(m), 0.0, true));
  Lab s, g;
  CK(codes(c, ar[2], "strat", m, s));
  if (ar[3].t != Val::NUL) {
    CK(codes(c, ar[3], "grp", m, g));
  } else {
    g.code.assign(m, 0);
    g.lev = {"all"};
    g.g = 1;
  }
  std::vector<double> w;
  CK(perstr(c, c["std"], "std", s.lev, w));
  int r0 = 1;
  CK(ref(c, c["ref"], g.lev, r0));
  --r0;
  double pr = 0.0, z = 0.0, lv = 0.0;
  CK(oned(c, "per", 1, 0.0, Inf, pr));
  CK(lvl(c, z, lv));
  const double a = 1 - lv;
  const int k = g.g;
  const std::size_t sk = static_cast<std::size_t>(k);
  std::vector<double> o(10 * sk);
  std::vector<int> seen(sk * static_cast<std::size_t>(s.g));
  double *rt = &o[0], *vr = &o[sk], *wm = &o[2 * sk], *ev = &o[3 * sk], *ra = &o[4 * sk], *rl = &o[5 * sk],
         *rh = &o[6 * sk], *df = &o[7 * sk], *dl = &o[8 * sk], *dh = &o[9 * sk];
  const int rc = clab_dastd(xv.p(), nv.p(), g.code.data(), s.code.data(), static_cast<int>(m), k, s.g, w.data(), r0,
                            z, rt, vr, wm, ev, ra, rl, rh, df, dl, dh, seen.data());
  if (rc == CLAB_EARG) return c.fail("each group must have one row per stratum.");
  if (rc != CLAB_OK) {
    return c.fail("every group needs a row for each stratum of positive weight, and `std` a positive sum.");
  }
  std::vector<double> cx(sk), cr(sk), cl(sk), ch(sk), cv(sk);
  for (std::size_t i = 0; i < sk; ++i) {
    const double y = rt[i], v = vr[i], wx = wm[i];
    cx[i] = ev[i];
    cr[i] = pr * y;
    cl[i] = y > 0 ? pr * v / (2 * y) * clab_qchisq(a / 2, 2 * y * y / v, 1) : 0.0;
    ch[i] = pr * (v + wx * wx) / (2 * (y + wx)) * clab_qchisq(1 - a / 2, 2 * (y + wx) * (y + wx) / (v + wx * wx), 1);
    cv[i] = pr * pr * v;
  }
  Frame rs;
  rs.str("group", g.lev);
  rs.num("x", cx);
  rs.num("rate", cr);
  rs.num("lo", cl);
  rs.num("hi", ch);
  rs.num("var", cv);
  std::vector<std::string> cg;
  std::vector<double> mr, ml, mh, md, mdl, mdh;
  for (int i = 0; i < k; ++i) {
    if (i == r0) continue;
    cg.push_back(g.lev[static_cast<std::size_t>(i)]);
    mr.push_back(na(ra[i]));
    ml.push_back(na(rl[i]));
    mh.push_back(na(rh[i]));
    md.push_back(pr * df[i]);
    mdl.push_back(pr * dl[i]);
    mdh.push_back(pr * dh[i]);
  }
  Frame cm;
  cm.str("group", cg);
  cm.num("ratio", mr);
  cm.num("lo", ml);
  cm.num("hi", mh);
  cm.num("diff", md);
  cm.num("dlo", mdl);
  cm.num("dhi", mdh);
  out = vobj();
  put(out, "rates", rs.o);
  put(out, "cmp", cm.o);
  put(out, "ref", vstr(g.lev[static_cast<std::size_t>(r0)]));
  put(out, "per", vnum(pr));
  put(out, "level", vnum(lv));
  return true;
}

bool sirex(Ctx &c, Val &out) {
  Args ar(c, {"x", "n", "strat", "grp"});
  CK(ar.byname(c));
  Dv xv, nv;
  CK(num(c, ar[0], "x", xv, -1, 0.0));
  const std::size_t m = xv.n();
  CK(num(c, ar[1], "n", nv, static_cast<long>(m), 0.0, true));
  Lab s, g;
  CK(codes(c, ar[2], "strat", m, s));
  if (ar[3].t == Val::NUL) return c.fail("`grp` must give the group (area) of each row.");
  CK(codes(c, ar[3], "grp", m, g));
  std::vector<double> sr;
  if (c.has("rates")) CK(perstr(c, c["rates"], "rates", s.lev, sr));
  double z = 0.0, lv = 0.0;
  CK(lvl(c, z, lv));
  const double a = 1 - lv;
  const int k = g.g;
  const std::size_t sk = static_cast<std::size_t>(k);
  std::vector<double> o(5 * sk), wk(2 * static_cast<std::size_t>(s.g));
  double prior[2];
  if (clab_sirex(xv.p(), nv.p(), g.code.data(), s.code.data(), static_cast<int>(m), k, s.g,
                 sr.empty() ? nullptr : sr.data(), &o[0], &o[sk], &o[2 * sk], &o[3 * sk], &o[4 * sk], prior,
                 wk.data()) != CLAB_OK) {
    return c.fail("every group needs a positive expected count: a stratum rate of 0 or a group with no exposure "
                  "leaves none.");
  }
  std::vector<double> cl(sk), ch(sk), fl(sk), fh(sk);
  std::vector<int> cf(sk);
  for (std::size_t i = 0; i < sk; ++i) {
    const double ob = o[i], e = o[sk + i], r = o[2 * sk + i];
    cl[i] = gmlo(ob, a) / e;
    ch[i] = gmhi(ob, a) / e;
    fl[i] = clab_qpois(a / 2, e, 1) / e;
    fh[i] = clab_qpois(1 - a / 2, e, 1) / e;
    cf[i] = r < fl[i] || r > fh[i] ? 1 : 0;
  }
  Frame cs;
  cs.str("group", g.lev);
  cs.num("obs", &o[0], sk);
  cs.num("expd", &o[sk], sk);
  cs.num("sir", &o[2 * sk], sk);
  cs.num("lo", cl);
  cs.num("hi", ch);
  cs.num("eb", &o[3 * sk], sk);
  cs.num("shr", &o[4 * sk], sk);
  cs.num("flo", fl);
  cs.num("fhi", fh);
  cs.lgl("out", cf.data(), sk);
  static const char *const pn[] = {"mean", "var"};
  out = vobj();
  put(out, "sir", cs.o);
  put(out, "prior", nvec(2, pn, prior));
  put(out, "level", vnum(lv));
  return true;
}

// The noncentrality at which pnchisq(st, df, ncp) is q (0 when ncp = 0
// leaves it below q), by bisection.
double ncpq(double st, double df, double q) {
  if (clab_pchisq(st, df, 1) <= q) return 0.0;
  double a = 0, b = std::fmax(1.0, st);
  while (clab_pnchisq(st, df, b, 1) > q) b *= 2;
  for (int it = 0; it < 200 && b - a > 1e-12 * b; ++it) {
    const double m = (a + b) / 2;
    if (clab_pnchisq(st, df, m, 1) > q) a = m;
    else b = m;
  }
  return (a + b) / 2;
}

const char *const aggs[] = {"n", "sum", "max", "min", "mean", "dist"};

bool ctabx(Ctx &c, Val &out) {
  bool yates = false, ex = false;
  CK(flagd(c, "correct", false, yates));
  CK(flagd(c, "exact", false, ex));
  double z = 0.0, conf = 0.0;
  CK(lvl(c, z, conf));
  int r = 0, cc = 0, ag = 0;
  std::vector<double> t;
  std::vector<std::string> rl, cl;
  if (!c.has("col")) {
    const Val &row = c["row"];
    int dim[2] = {-1, -1};
    if (!(row.t == Val::ARR && !row.a.empty() && row.a[0].t == Val::ARR && shape(row, 2, dim, 0))) {
      return c.fail("`row` must be the row labels (with `col`) or a numeric matrix or 2-d table.");
    }
    if (c.has("w")) return c.fail("`w` goes with `row` and `col` labels, not with a table.");
    r = dim[0];
    cc = dim[1];
    t.assign(static_cast<std::size_t>(r) * cc, 0.0);
    bool bad = false;
    fill(row, 2, dim, 0, 0, 1, t, bad, true);
    for (double e : t) {
      if (!(e >= 0) || !std::isfinite(e)) return c.fail("`row` must hold counts: finite numbers from 0.");
    }
    rl = nmsof({}, static_cast<std::size_t>(r));
    cl = nmsof({}, static_cast<std::size_t>(cc));
    CK(pickd(c, "agg", aggs, 6, 0, ag));
    if (ag > 1) return c.fail("`agg` applies to `row` and `col` labels, not to a table.");
  } else {
    Args a(c, {"row", "col", "w"});
    CK(a.byname(c));
    const std::size_t n = vlen(a[0]);
    if (n < 1) return c.fail("`row` must have at least one label.");
    Lab la, lb;
    CK(codes(c, a[0], "row", n, la));
    CK(codes(c, a[1], "col", n, lb));
    r = la.g;
    cc = lb.g;
    rl = la.lev;
    cl = lb.lev;
    CK(pickd(c, "agg", aggs, 6, a[2].t == Val::NUL ? 0 : 1, ag));
    const std::size_t rc = static_cast<std::size_t>(r) * static_cast<std::size_t>(cc);
    t.assign(rc, 0.0);
    std::vector<double> cnt(rc, 0.0);
    auto at = [&](std::size_t i) {
      return static_cast<std::size_t>(lb.code[i]) * static_cast<std::size_t>(r) + static_cast<std::size_t>(la.code[i]);
    };
    if (ag == 0) {
      if (a[2].t != Val::NUL) return c.fail("`agg = \"n\"` counts rows: leave `w` out.");
      for (std::size_t i = 0; i < n; ++i) t[at(i)] += 1;
    } else if (a[2].t == Val::NUL) {
      return c.fail("`agg = \"%s\"` needs `w`, the values to aggregate.", aggs[ag]);
    } else if (ag == 5) {
      Lab d;
      CK(codes(c, a[2], "w", n, d));
      std::vector<std::vector<int>> seen(rc);
      for (std::size_t i = 0; i < n; ++i) seen[at(i)].push_back(d.code[i]);
      for (std::size_t q = 0; q < rc; ++q) {
        std::sort(seen[q].begin(), seen[q].end());
        t[q] = static_cast<double>(std::unique(seen[q].begin(), seen[q].end()) - seen[q].begin());
      }
    } else {
      Dv wv;
      CK(num(c, a[2], "w", wv, static_cast<long>(n), ag == 1 ? 0.0 : -HUGE_VAL));
      for (std::size_t i = 0; i < n; ++i) {
        const std::size_t q = at(i);
        const double v = wv[i];
        if (ag == 1 || ag == 4) t[q] += v;
        else if (cnt[q] == 0) t[q] = v;
        else if (ag == 2) t[q] = std::fmax(t[q], v);
        else t[q] = std::fmin(t[q], v);
        cnt[q] += 1;
      }
      if (ag == 4) {
        for (std::size_t q = 0; q < rc; ++q) t[q] = cnt[q] > 0 ? t[q] / cnt[q] : 0.0;
      }
    }
  }
  const std::size_t sr = static_cast<std::size_t>(r);
  std::vector<int> kr, kc;
  for (int i = 0; i < r; ++i) {
    double s = 0;
    for (int j = 0; j < cc; ++j) s += t[static_cast<std::size_t>(j) * sr + static_cast<std::size_t>(i)];
    if (s > 0) kr.push_back(i);
  }
  for (int j = 0; j < cc; ++j) {
    double s = 0;
    for (int i = 0; i < r; ++i) s += t[static_cast<std::size_t>(j) * sr + static_cast<std::size_t>(i)];
    if (s > 0) kc.push_back(j);
  }
  const int r2 = static_cast<int>(kr.size()), c2 = static_cast<int>(kc.size());
  const std::size_t s2 = static_cast<std::size_t>(r2);
  const bool counts = ag == 0 || ag == 1 || ag == 5;
  std::vector<double> u(s2 * static_cast<std::size_t>(c2)), e(u.size()), rs(u.size()), st(9, NaN);
  std::vector<int> dg(s2, -1);
  for (std::size_t j = 0; j < static_cast<std::size_t>(c2); ++j) {
    for (std::size_t i = 0; i < s2; ++i) u[j * s2 + i] = t[static_cast<std::size_t>(kc[j]) * sr + static_cast<std::size_t>(kr[i])];
  }
  for (std::size_t i = 0; i < s2; ++i) {
    for (int j = 0; j < c2; ++j) {
      if (rl[static_cast<std::size_t>(kr[i])] == cl[static_cast<std::size_t>(kc[static_cast<std::size_t>(j)])]) dg[i] = j;
    }
  }
  bool any = false;
  for (int v : dg) any = any || v >= 0;
  double fs[5] = {NaN, NaN, NaN, NaN, NaN};
  const bool ok = counts && r2 >= 1 && c2 >= 1;
  if (ok) {
    if (clab_ctabx(u.data(), r2, c2, yates ? 1 : 0, any ? dg.data() : nullptr, e.data(), rs.data(), st.data()) !=
        CLAB_OK) {
      return c.fail("the kernel refused the table.");
    }
    if (ex && r2 >= 2 && c2 >= 2) {
      std::vector<double> rem(static_cast<std::size_t>(r2 + c2));
      if (clab_fshx(u.data(), r2, c2, 1e7, conf, fs, rem.data()) != CLAB_OK) {
        for (double v : u) {
          if (v != std::floor(v)) return c.fail("the exact test needs whole counts.");
        }
        return c.fail("the exact test would visit more than 1e7 tables: leave `exact` FALSE.");
      }
    }
  } else if (ex) {
    return c.fail("the exact test needs counts: `agg` \"n\", \"sum\" or \"dist\".");
  }
  const bool two = ok && r2 >= 2 && c2 >= 2;
  for (double &v : e) v = ok ? v : NA;
  for (double &v : rs) v = ok ? v : NA;
  Frame sc;
  auto one1 = [&](const char *nm, double v) {
    const double x = std::isnan(v) ? NA : v;
    sc.num(nm, &x, 1);
  };
  const double df = st[2], ch = st[1], nn = st[0];
  const int kk = r2 < c2 ? r2 : c2;
  one1("n", st[0]);
  one1("chisq", two ? ch : NaN);
  one1("df", two ? df : NaN);
  one1("p", two ? clab_pchisq(ch, df, 0) : NaN);
  one1("g2", two ? st[3] : NaN);
  one1("g2p", two ? clab_pchisq(st[3], df, 0) : NaN);
  one1("v", two ? st[4] : NaN);
  const double chu = two ? (yates && r2 == 2 && c2 == 2 ? st[4] * st[4] * nn * (kk - 1) : ch) : NaN;
  one1("vlo", two ? std::fmin(1.0, std::sqrt(ncpq(chu, df, 1 - (1 - conf) / 2) / (nn * (kk - 1)))) : NaN);
  one1("vhi", two ? std::fmin(1.0, std::sqrt(ncpq(chu, df, (1 - conf) / 2) / (nn * (kk - 1)))) : NaN);
  one1("vbc", two ? st[5] : NaN);
  one1("cc", two ? st[6] : NaN);
  one1("phi", st[7]);
  one1("diag", st[8]);
  one1("fp", fs[0]);
  one1("or", fs[1]);
  one1("orlo", fs[2]);
  one1("orhi", fs[3]);
  out = vobj();
  put(out, "tab", vmat(t, r, cc));
  put(out, "expd", vmat(e, r2, c2));
  put(out, "res", vmat(rs, r2, c2));
  put(out, "stat", sc.o);
  put(out, "drop", vnum((r - r2) + (cc - c2)));
  put(out, "agg", vstr(aggs[ag]));
  put(out, "correct", vbool(yates));
  put(out, "level", vnum(conf));
  return true;
}
