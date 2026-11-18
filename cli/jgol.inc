// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The generalised ordered logit (src/gol.cpp, src/bgol.cpp): R/gol.R in
// C++, and the table of its slopes by threshold.

namespace gol {

// .golev(): the levels of the outcome column a, lowest first, as text.
inline bool lev(Ctx &c, const Val &a, const std::string &y, std::vector<std::string> &lv) {
  lv.clear();
  if (c.has("lev")) {
    const Val l = ax::asarr(c["lev"]);
    bool ok = !l.a.empty();
    for (const Val &e : l.a) ok = ok && (e.t == Val::STR || e.t == Val::NUM);
    for (const Val &e : l.a) lv.push_back(ax::vlab(e));
    for (std::size_t i = 0; ok && i < lv.size(); ++i) {
      for (std::size_t j = 0; j < i; ++j) ok = ok && lv[i] != lv[j];
    }
    if (!ok) return ckstop(c, "`lev` must be distinct levels with no missing values");
  } else if (ax::isnumc(a)) {
    std::vector<double> v = ax::colv(a);
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    for (double e : v) lv.push_back(ax::vlab(vnum(e)));
  } else {
    return ckstop(c, "`lev` must give the order of the levels of column `%s`, text", y.c_str());
  }
  std::vector<std::string> bad, miss;
  for (const Val &e : a.a) {
    const std::string s = ax::vlab(e);
    if (std::find(lv.begin(), lv.end(), s) == lv.end() && std::find(bad.begin(), bad.end(), s) == bad.end()) {
      bad.push_back(s);
    }
  }
  if (!bad.empty()) {
    if (bad.size() > 5) bad.resize(5);
    return ckstop(c, "column `%s` has values not in its levels: %s", y.c_str(), ax::join(bad, ", ").c_str());
  }
  for (const std::string &l : lv) {
    bool in = false;
    for (const Val &e : a.a) in = in || ax::vlab(e) == l;
    if (!in) miss.push_back(l);
  }
  if (!miss.empty()) {
    return ckstop(c, "column `%s` has no rows at level(s) %s", y.c_str(), ax::join(miss, ", ").c_str());
  }
  if (lv.size() < 3) return ckstop(c, "column `%s` must have at least three levels", y.c_str());
  return true;
}

}  // namespace gol

bool gologt(Ctx &c, Val &out) {
  Val d;
  CK(ax::dat(c, c["data"], "data", d));
  ax::Cols yc, x;
  CK(ax::acol(c, d, "y", {1}, false, yc));
  CK(ax::acol(c, d, "x", {}, false, x));
  std::vector<double> xm;
  CK(ax::amat(c, d, x, xm));
  std::vector<int> pv(x.size(), 0);
  if (c.has("par")) {
    const Val pa = ax::asarr(c["par"]);
    bool ok = !pa.a.empty();
    std::vector<std::string> seen;
    for (const Val &e : pa.a) {
      ok = ok && e.t == Val::STR && ax::has(x, e.s) && !ax::has(seen, e.s);
      if (ok) seen.push_back(e.s);
    }
    if (!ok) return ckstop(c, "`par` must name columns of `x`, each once");
    for (std::size_t k = 0; k < x.size(); ++k) pv[k] = ax::has(seen, x[k]) ? 1 : 0;
  }
  const Val &ya = *ax::fc(d, yc[0]);
  std::vector<std::string> lv;
  CK(gol::lev(c, ya, yc[0], lv));
  const int n = sat<int>(ax::nr(d)), M = sat<int>(lv.size()), J = M - 1;
  const int p = sat<int>(x.size());
  std::vector<int> yv(sat<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    const std::string s = ax::vlab(ya.a[sat<std::size_t>(i)]);
    yv[sat<std::size_t>(i)] = sat<int>(std::find(lv.begin(), lv.end(), s) - lv.begin());
  }
  Nv lvl, mi;
  CK(cknumd(c, "level", 0.95, lvl));
  if (lvl[0] <= 0 || lvl[0] >= 1) return ckstop(c, "`level` must be strictly between 0 and 1");
  CK(cknumd(c, "maxit", 500, mi, 1, HUGE_VAL, false, true));
  const int mx = mi[0] > 2147483647.0 ? 2147483647 : sat<int>(mi[0]);
  int np = J, fr = 0;
  for (int k = 0; k < p; ++k) {
    np += pv[sat<std::size_t>(k)] ? 1 : J;
    fr += pv[sat<std::size_t>(k)] ? 0 : 1;
  }
  const std::size_t snp = sat<std::size_t>(np);
  std::vector<double> est(snp), cov(snp * snp), wd(3 * sat<std::size_t>(p)), st(3);
  const char *dom = "the likelihood could not be maximised from its start";
  int rc = clab_gologt(yv.data(), n, M, xm.data(), p, pv.data(), mx, est.data(), cov.data(), wd.data(), st.data());
  if (rc == CLAB_EDOM) return ckstop(c, "%s", dom);
  CK(kchk(c, rc));
  double lr[3] = {NA, NA, NA};
  if (fr > 0) {
    const std::size_t n0 = sat<std::size_t>(J + p);
    std::vector<int> all(sat<std::size_t>(p), 1);
    std::vector<double> e0(n0), c0(n0 * n0), w0(3 * sat<std::size_t>(p)), s0(3);
    rc = clab_gologt(yv.data(), n, M, xm.data(), p, all.data(), mx, e0.data(), c0.data(), w0.data(), s0.data());
    if (rc == CLAB_EDOM) return ckstop(c, "%s", dom);
    CK(kchk(c, rc));
    lr[0] = std::fmax(2 * (st[0] - s0[0]), 0.0);
    lr[1] = static_cast<double>(fr) * (J - 1);
    lr[2] = clab_pchisq(lr[0], lr[1], 0);
  }
  const double z = qnorm(0.5 + 0.5 * lvl[0]);
  std::vector<std::string> term, thr, lb;
  for (int j = 0; j < J; ++j) lb.push_back(lv[sat<std::size_t>(j)] + "|" + lv[sat<std::size_t>(j + 1)]);
  for (int j = 0; j < J; ++j) {
    term.push_back("(Intercept)");
    thr.push_back(lb[sat<std::size_t>(j)]);
  }
  for (int k = 0; k < p; ++k) {
    for (int j = 0; j < (pv[sat<std::size_t>(k)] ? 1 : J); ++j) {
      term.push_back(x[sat<std::size_t>(k)]);
      thr.push_back(pv[sat<std::size_t>(k)] ? "all" : lb[sat<std::size_t>(j)]);
    }
  }
  std::vector<double> se(snp), zv(snp), pr(snp), lo(snp), hi(snp);
  for (std::size_t k = 0; k < snp; ++k) {
    const double e = est[k], s = std::sqrt(cov[k * snp + k]);
    se[k] = s;
    zv[k] = e / s;
    pr[k] = 2 * clab_pnorm(-std::fabs(e / s));
    lo[k] = e - z * s;
    hi[k] = e + z * s;
  }
  Frame f, w;
  f.str("term", term);
  f.str("thr", thr);
  f.num("est", est);
  f.num("se", se);
  f.num("z", zv);
  f.num("p", pr);
  f.num("lo", lo);
  f.num("hi", hi);
  std::vector<double> ws, wf, wp;
  for (int k = 0; k < p; ++k) {
    ws.push_back(wd[3 * sat<std::size_t>(k)]);
    wf.push_back(wd[3 * sat<std::size_t>(k) + 1]);
    wp.push_back(wd[3 * sat<std::size_t>(k) + 2]);
  }
  w.str("term", x);
  w.num("stat", ws);
  w.num("df", wf);
  w.num("p", wp);
  static const char *const ln[] = {"stat", "df", "p"}, *const fn[] = {"loglik", "aic", "iter", "conv"};
  const double ft[4] = {st[0], -2 * st[0] + 2 * np, st[1], st[2]};
  out = vobj();
  put(out, "coef", f.o);
  put(out, "wald", w.o);
  put(out, "lr", nvec(3, ln, lr));
  put(out, "fit", nvec(4, fn, ft));
  put(out, "vcov", vmat(cov, np, np));
  put(out, "lev", vstrs(lv));
  put(out, "n", vnum(n));
  return true;
}

bool gocoef(Ctx &c, Val &out) {
  const Val &o = c["object"];
  const Val &cf = o.t == Val::OBJ ? col(o, "coef") : nul, &wd = o.t == Val::OBJ ? col(o, "wald") : nul;
  const Val &tv = cf.t == Val::OBJ ? col(cf, "term") : nul, &hv = cf.t == Val::OBJ ? col(cf, "thr") : nul;
  const Val &ev = cf.t == Val::OBJ ? col(cf, "est") : nul, &sv = cf.t == Val::OBJ ? col(cf, "se") : nul;
  const Val &wt = wd.t == Val::OBJ ? col(wd, "term") : nul;
  const std::size_t n = vlen(tv);
  bool ok = n > 0 && tv.t == Val::ARR && hv.t == Val::ARR && vlen(hv) == n && vlen(ev) == n && vlen(sv) == n &&
            wt.t == Val::ARR && isnum(ev) && isnum(sv);
  for (std::size_t i = 0; ok && i < n; ++i) ok = tv.a[i].t == Val::STR && hv.a[i].t == Val::STR;
  for (std::size_t j = 0; ok && j < wd.a.size(); ++j) ok = wd.a[j].t == Val::ARR && wd.a[j].a.size() == wt.a.size();
  if (!ok) return ckstop(c, "`object` must be a result of gologt()");
  const std::vector<double> e = dbl(ev), s = dbl(sv);
  std::vector<std::string> tm, lb;
  for (std::size_t i = 0; i < n; ++i) {
    if (!ax::has(tm, tv.a[i].s)) tm.push_back(tv.a[i].s);
    if (tv.a[i].s == "(Intercept)") lb.push_back(hv.a[i].s);
  }
  std::vector<std::string> v = tm;
  if (c.has("v")) {
    const Val a = ax::asarr(c["v"]);
    v.clear();
    bool good = !a.a.empty();
    for (const Val &x : a.a) {
      good = good && x.t == Val::STR && ax::has(tm, x.s);
      if (good) v.push_back(x.s);
    }
    if (!good) return ckstop(c, "`v` must name terms of the fit: %s", ax::join(tm, ", ").c_str());
  }
  const std::size_t J = lb.size(), K = v.size();
  std::vector<double> me(J * K), ms(J * K);
  for (std::size_t k = 0; k < K; ++k) {
    std::vector<std::size_t> r;
    for (std::size_t i = 0; i < n; ++i) {
      if (tv.a[i].s == v[k]) r.push_back(i);
    }
    for (std::size_t j = 0; j < J; ++j) {
      const std::size_t i = r.size() == 1 ? r[0] : r[j < r.size() ? j : r.size() - 1];
      me[k * J + j] = e[i];
      ms[k * J + j] = s[i];
    }
  }
  // the rows of wald for the covariates of v
  std::vector<std::size_t> keep;
  for (std::size_t i = 0; i < wt.a.size(); ++i) {
    if (wt.a[i].t == Val::STR && ax::has(v, wt.a[i].s)) keep.push_back(i);
  }
  Val w = vobj();
  for (std::size_t j = 0; j < wd.k.size(); ++j) {
    Val cv = varr();
    for (std::size_t i : keep) cv.a.push_back(wd.a[j].a[i]);
    put(w, wd.k[j].c_str(), cv);
  }
  out = vobj();
  put(out, "est", vmat(me.data(), sat<int>(J), sat<int>(K)));
  put(out, "se", vmat(ms.data(), sat<int>(J), sat<int>(K)));
  put(out, "wald", w);
  return true;
}
