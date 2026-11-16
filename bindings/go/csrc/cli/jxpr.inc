// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The experiments, the trials and the survey tests (src/rnd.cpp and the
// kernels after it) for clab_jrun(): the R wrappers in C++, with their
// checks in their order and words. In R a seed left out is drawn from R's
// generator; outside R there is none, so a function that draws needs
// `seed`.

namespace ax {

// .xpseed(): one whole number from 0 to 2147483647.
inline bool xseed(Ctx &c, unsigned int &s) {
  if (!c.has("seed")) return c.fail("`seed` is needed outside R: one whole number from 0.");
  Nv q;
  CK(cknum(c, c["seed"], "seed", q, true, 0, 2147483647.0, false, true));
  s = static_cast<unsigned int>(q[0]);
  return true;
}

// Text values, each a string (is.character() with no NA).
inline bool strs(const Val &o, Cols &s) {
  s.clear();
  const Val a = asarr(o);
  for (const Val &e : a.a) {
    if (e.t != Val::STR) return false;
    s.push_back(e.s);
  }
  return true;
}

// One string.
inline bool str1(Ctx &c, const char *nm, const char *def, std::string &s) {
  s = def;
  if (!c.has(nm)) return true;
  Cols v;
  if (!strs(c[nm], v) || v.size() != 1) return ckstop(c, "`%s` must be one string", nm);
  s = v[0];
  return true;
}

// rndalc() (src/rnd.cpp): R/rndalc.R in C++. strata is an object from
// stratum to places.
bool rndalc(Ctx &c, Val &out) {
  Cols arms = {"A", "B"};
  if (c.has("arms")) {
    bool ok = strs(c["arms"], arms) && arms.size() >= 2;
    for (std::size_t i = 0; ok && i < arms.size(); ++i) {
      ok = !arms[i].empty();
      for (std::size_t j = 0; ok && j < i; ++j) ok = arms[i] != arms[j];
    }
    if (!ok) return ckstop(c, "`arms` must name two or more arms, each once");
  }
  const int k = static_cast<int>(arms.size());
  Nv ra;
  if (c.has("ratio")) {
    CK(cknum(c, c["ratio"], "ratio", ra, false, 1, 1e6, false, true));
  } else {
    ra.v.assign(arms.size(), 1.0);
  }
  if (ra.n() != arms.size()) return ckstop(c, "`ratio` must give one number per arm");
  double sr = 0;
  for (double v : ra.v) sr += v;
  Nv bl;
  if (c.has("block")) {
    CK(cknum(c, c["block"], "block", bl, false, 0, 1e6, false, true));
  } else {
    bl.v = {2 * sr, 3 * sr};
  }
  const bool simple = bl.n() == 1 && bl[0] == 0;
  for (double v : bl.v) {
    if (!simple && (v == 0 || rmod(v, sr) != 0)) {
      return ckstop(c, "`block` must be 0, or sizes that are multiples of sum(ratio) (%s)", rfmt(sr).c_str());
    }
  }
  std::vector<double> cnt;
  Cols nm;
  Nv q;
  if (!c.has("strata")) {
    CK(cknum(c, c["n"], "n", q, true, 1, 1e7, false, true));
    cnt = q.v;
  } else {
    CK(cknum(c, c["strata"], "strata", q, false, 0, Inf, false, true));
    cnt = q.v;
    const Val &so = c["strata"];
    bool ok = so.t == Val::OBJ;
    if (ok) nm = so.k;
    for (std::size_t i = 0; ok && i < nm.size(); ++i) {
      ok = !nm[i].empty();
      for (std::size_t j = 0; ok && j < i; ++j) ok = nm[i] != nm[j];
    }
    if (!ok) return ckstop(c, "`strata` must be named, each stratum once");
    double tot = 0;
    for (double v : cnt) tot += v;
    if (tot < 1 || tot > 1e7) return ckstop(c, "`strata` must hold from 1 to 10000000 places in all");
    if (c.has("n")) {
      CK(cknum(c, c["n"], "n", q, true, 1, 1e7, false, true));
      if (q[0] != tot) return ckstop(c, "`n` must equal the places of `strata` (%s)", rfmt(tot).c_str());
    }
  }
  unsigned int sd = 0;
  CK(xseed(c, sd));
  std::string pre;
  CK(str1(c, "prefix", "P", pre));
  std::vector<int> ci, ri, bi;
  int n = 0;
  for (double v : cnt) {
    ci.push_back(static_cast<int>(v));
    n += static_cast<int>(v);
  }
  for (double v : ra.v) ri.push_back(static_cast<int>(v));
  if (!simple) {
    for (double v : bl.v) bi.push_back(static_cast<int>(v));
  }
  std::vector<int> arm(zz(n)), bno(zz(n));
  const int rc = clab_rndalc(static_cast<int>(ci.size()), ci.data(), k, ri.data(), bi.empty() ? nullptr : bi.data(),
                             static_cast<int>(bi.size()), sd, arm.data(), bno.data());
  if (rc != CLAB_OK) return kerr(c, "rndalc", rc);
  const int w = static_cast<int>(std::to_string(n).size());
  Cols id, st, av;
  Val bv = varr();
  std::vector<double> tab(ci.size() * zz(k), 0.0);
  for (int i = 0, s = 0, used = 0; i < n; ++i) {
    while (used == ci[zz(s)]) {
      ++s;
      used = 0;
    }
    ++used;
    std::string d = std::to_string(i + 1);
    id.push_back(pre + std::string(zz(w) - d.size(), '0') + d);
    if (!nm.empty()) st.push_back(nm[zz(s)]);
    bv.a.push_back(vnum(bno[zz(i)] > 0 ? bno[zz(i)] : NA));
    av.push_back(arms[zz(arm[zz(i)])]);
    tab[zz(s) + ci.size() * zz(arm[zz(i)])] += 1;
  }
  Frame f;
  f.str("id", id);
  if (!nm.empty()) f.str("stratum", st);
  f.add("block", bv);
  f.str("arm", av);
  out = vobj();
  put(out, "alloc", f.o);
  put(out, "tab", vmat(tab, static_cast<int>(ci.size()), k));
  put(out, "seed", vnum(sd));
  return true;
}

// choose() of R for whole n >= k >= 0: the product for k < 30, else from
// the log gamma function, rounded to a whole number.
inline double rchoose(double n, double k) {
  if (n - k < k) k = n - k;
  if (k < 30) {
    double r = k < 1 ? 1 : n;
    for (int j = 2; j <= k; ++j) r *= (n - j + 1) / j;
    return std::nearbyint(r);
  }
  return std::nearbyint(std::exp(std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1)));
}

// blkprm() (src/bpm.cpp): R/blkprm.R in C++. Blocks are coded in order
// of first appearance.
bool blkprm(Ctx &c, Val &out) {
  Val d;
  CK(dat(c, c["data"], "data", d));
  Cols yc, dc, bc;
  CK(acol(c, d, "y", {1}, false, yc));
  CK(acol(c, d, "d", {1}, false, dc));
  CK(acol(c, d, "block", {1}, false, bc));
  int st = 0, al = 0;
  CK(marg(c, "stat", {"diff", "strat"}, st));
  CK(marg(c, "alt", {"two.sided", "greater", "less"}, al));
  std::vector<double> yv, tv;
  CK(amat(c, d, yc, yv));
  CK(amat(c, d, dc, tv));
  double n1 = 0;
  for (double v : tv) {
    if (v != 0 && v != 1) return ckstop(c, "column `%s` must hold 0 or 1", dc[0].c_str());
    n1 += v;
  }
  const std::size_t n = nr(d);
  if (n1 == 0 || n1 == static_cast<double>(n)) {
    return ckstop(c, "column `%s` must hold treated and control units", dc[0].c_str());
  }
  const Code b = acode(*fc(d, bc[0]));
  std::vector<double> m(zz(b.n), 0.0), t1(zz(b.n), 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    m[zz(b.c[i])] += 1;
    t1[zz(b.c[i])] += tv[i];
  }
  if (st == 1) {
    for (int k = 0; k < b.n; ++k) {
      if (t1[zz(k)] == 0 || t1[zz(k)] == m[zz(k)]) {
        return ckstop(c, "with `stat = \"strat\"` every block needs treated and control units");
      }
    }
  }
  bool ex = false;
  CK(flagd(c, "exact", false, ex));
  double np = 0;
  unsigned int sd = 0;
  if (ex) {
    np = 1;
    for (int k = 0; k < b.n; ++k) np *= rchoose(m[zz(k)], t1[zz(k)]);
    if (np > 1e6) {
      return ckstop(c, "`exact` enumerates at most 1000000 arrangements; these blocks have %s", rfmt(np).c_str());
    }
  } else {
    Nv q;
    CK(cknumd(c, "nperm", 2000, q, 1, 1e7, false, true));
    np = q[0];
    CK(xseed(c, sd));
  }
  std::vector<int> di(n);
  for (std::size_t i = 0; i < n; ++i) di[i] = static_cast<int>(tv[i]);
  const int nd = static_cast<int>(np);
  std::vector<double> s4(4), dist(zz(nd));
  const int rc = clab_blkprm(yv.data(), di.data(), b.c.data(), static_cast<int>(n), b.n, st, al, nd, ex ? 1 : 0, sd,
                             nd, s4.data(), dist.data());
  if (rc != CLAB_OK) return kerr(c, "blkprm", rc);
  static const char *const sn[] = {"diff", "strat"}, *const an[] = {"two.sided", "greater", "less"};
  out = vobj();
  put(out, "est", vnum(s4[0]));
  put(out, "p", vnum(s4[1]));
  put(out, "nperm", vnum(s4[2]));
  put(out, "nge", vnum(s4[3]));
  put(out, "nblock", vnum(b.n));
  put(out, "stat", vstr(sn[st]));
  put(out, "alt", vstr(an[al]));
  put(out, "exact", vbool(ex));
  put(out, "dist", vnums(dist));
  put(out, "seed", vnum(ex ? NA : static_cast<double>(sd)));
  return true;
}

// pwranv() (src/pwa.cpp): R/pwranv.R in C++, the grid of n and f with n
// varying fastest.
bool pwranv(Ctx &c, Val &out) {
  Nv kv, fv, av, nv, pv;
  CK(cknum(c, c["k"], "k", kv, true, 2, 1e6, false, true));
  CK(cknum(c, c["f"], "f", fv, false, 0));
  CK(cknumd(c, "alpha", 0.05, av));
  const double alpha = av[0];
  if (alpha <= 0 || alpha >= 1) return ckstop(c, "`alpha` must be one number strictly between 0 and 1");
  if (c.has("n") == c.has("power")) return ckstop(c, "give one of `n` and `power`");
  std::vector<double> n, f;
  double tg = NaN;
  if (c.has("n")) {
    CK(cknum(c, c["n"], "n", nv, false, 1, 1e7));
    for (double v : nv.v) {
      if (v <= 1) return ckstop(c, "`n` must be above 1");
    }
    for (double b : fv.v) {
      for (double a : nv.v) {
        n.push_back(a);
        f.push_back(b);
      }
    }
  } else {
    CK(cknum(c, c["power"], "power", pv));
    tg = pv[0];
    if (tg <= alpha || tg >= 1) return ckstop(c, "`power` must be above `alpha` and below 1");
    f = fv.v;
  }
  const int m = static_cast<int>(f.size());
  std::vector<double> o(zz(m) * 6);
  const int rc = clab_pwranv(static_cast<int>(kv[0]), m, n.empty() ? nullptr : n.data(), f.data(), alpha, tg, o.data());
  if (rc != CLAB_OK) return kerr(c, "pwranv", rc);
  for (double &v : o) {
    if (std::isnan(v)) v = NA;
  }
  Frame fr;
  fr.num("k", std::vector<double>(zz(m), kv[0]));
  fr.num("n", o.data(), zz(m));
  fr.num("f", f);
  static const char *const nm[] = {"df1", "df2", "ncp", "fcrit", "power"};
  for (int j = 0; j < 5; ++j) fr.num(nm[j], o.data() + zz(j + 1) * zz(m), zz(m));
  out = fr.o;
  return true;
}

// .xpopen(): one number strictly between 0 and 1, def when left out.
inline bool xopen(Ctx &c, const char *nm, double def, double &v) {
  Nv q;
  CK(cknumd(c, nm, def, q));
  v = q[0];
  if (v <= 0 || v >= 1) return ckstop(c, "`%s` must be one number strictly between 0 and 1", nm);
  return true;
}

// pwrsim() (src/pws.cpp): R/pwrsim.R in C++, the grid of n and eff with
// n varying fastest.
bool pwrsim(Ctx &c, Val &out) {
  Nv nv, ev, q;
  CK(cknum(c, c["n"], "n", nv, false, 2, 1e6, false, true));
  CK(cknum(c, c["eff"], "eff", ev, false));
  int fam = 0, test = 0, alt = 0;
  CK(marg(c, "fam", {"normal", "binary", "poisson"}, fam));
  CK(marg(c, "test", {"welch", "perm"}, test));
  CK(marg(c, "alt", {"two.sided", "greater", "less"}, alt));
  if (fam == 0) {
    CK(cknumd(c, "sd", 1, q, -Inf, Inf, true));
  } else if (fam == 1) {
    CK(cknumd(c, "p0", 0.5, q, 0, 1));
  } else {
    CK(cknumd(c, "mu0", 1, q, -Inf, Inf, true));
  }
  const double par = q[0];
  for (double e : ev.v) {
    if (fam == 1 && (par + e < 0 || par + e > 1)) return ckstop(c, "`p0 + eff` must be from 0 to 1");
    if (fam == 2 && par + e < 0) return ckstop(c, "`mu0 + eff` must be at least 0");
  }
  CK(cknumd(c, "ratio", 1, q, -Inf, Inf, true));
  const double ratio = q[0];
  double alpha = 0, level = 0;
  CK(xopen(c, "alpha", 0.05, alpha));
  CK(cknumd(c, "nsim", 1000, q, 1, 1e6, false, true));
  const int nsim = static_cast<int>(q[0]);
  int nperm = 1;
  if (test == 1) {
    CK(cknumd(c, "nperm", 199, q, 1, 1e5, false, true));
    if (1 / (q[0] + 1) > alpha) {
      return ckstop(c, "`nperm` is too few for `alpha`: the least p-value, 1 / (nperm + 1), is above it");
    }
    nperm = static_cast<int>(q[0]);
  }
  CK(xopen(c, "level", 0.95, level));
  unsigned int sd = 0;
  CK(xseed(c, sd));
  std::vector<int> n0, n1;
  std::vector<double> gn, gn1, ge;
  for (double e : ev.v) {
    for (double a : nv.v) {
      const double b = std::nearbyint(ratio * a);
      if (b < 2 || b > 1e6) return ckstop(c, "`ratio * n` must round to treated arms of 2 to 1000000");
      n0.push_back(static_cast<int>(a));
      n1.push_back(static_cast<int>(b));
      gn.push_back(a);
      gn1.push_back(b);
      ge.push_back(e);
    }
  }
  const int m = static_cast<int>(ge.size());
  std::vector<double> o(zz(m) * 5);
  const int rc = clab_pwrsim(fam, par, m, n0.data(), n1.data(), ge.data(), test, alt, alpha, nsim, nperm, level,
                             sd, o.data());
  if (rc != CLAB_OK) return kerr(c, "pwrsim", rc);
  static const char *const fn[] = {"normal", "binary", "poisson"}, *const tn[] = {"welch", "perm"},
                           *const an[] = {"two.sided", "greater", "less"};
  Frame f;
  f.num("n", gn);
  f.num("n1", gn1);
  f.num("eff", ge);
  static const char *const nm[] = {"power", "se", "lo", "hi", "rej"};
  for (int j = 0; j < 5; ++j) f.num(nm[j], o.data() + zz(j) * zz(m), zz(m));
  out = vobj();
  put(out, "tab", f.o);
  put(out, "fam", vstr(fn[fam]));
  put(out, "test", vstr(tn[test]));
  put(out, "alt", vstr(an[alt]));
  put(out, "alpha", vnum(alpha));
  put(out, "nsim", vnum(nsim));
  put(out, "seed", vnum(sd));
  return true;
}

// falsctl() (src/fls.cpp): R/falsctl.R in C++.
bool falsctl(Ctx &c, Val &out) {
  Val d;
  CK(dat(c, c["data"], "data", d));
  Cols yc, dc, nc, xc;
  CK(acol(c, d, "y", {1}, true, yc));
  CK(acol(c, d, "d", {1}, false, dc));
  CK(acol(c, d, "nc", {}, false, nc));
  CK(acol(c, d, "x", {}, true, xc));
  bool clash = has(xc, dc[0]) || has(yc, dc[0]) || (!yc.empty() && has(xc, yc[0]));
  for (const std::string &s : nc) clash = clash || has(yc, s) || s == dc[0] || has(xc, s);
  if (clash) return ckstop(c, "`y`, `d`, `nc` and `x` must name different columns");
  int se = 0, adj = 0;
  CK(marg(c, "se", {"hc1", "iid"}, se));
  CK(marg(c, "adj", {"holm", "bonferroni", "hochberg", "BH", "BY", "none"}, adj));
  Cols oc = yc;
  oc.insert(oc.end(), nc.begin(), nc.end());
  std::vector<double> ym, dv, xm;
  CK(amat(c, d, oc, ym));
  CK(amat(c, d, dc, dv));
  CK(amat(c, d, xc, xm));
  const std::size_t n = nr(d);
  if (n <= xc.size() + 2) return ckstop(c, "`data` needs more rows than covariates plus two");
  bool vary = false;
  for (double v : dv) vary = vary || v != dv[0];
  if (!vary) return ckstop(c, "column `%s` must vary", dc[0].c_str());
  Nv q;
  CK(cknumd(c, "nperm", 0, q, 0, 1e6, false, true));
  const int np = static_cast<int>(q[0]);
  double alpha = 0, level = 0;
  CK(xopen(c, "alpha", 0.05, alpha));
  CK(xopen(c, "level", 0.95, level));
  unsigned int sd = 0;
  if (np > 0) CK(xseed(c, sd));
  std::vector<int> ctl;
  if (!yc.empty()) ctl.push_back(0);
  ctl.insert(ctl.end(), nc.size(), 1);
  static const int code[] = {2, 1, 3, 4, 5, 0};
  const int m = static_cast<int>(oc.size());
  std::vector<double> o(zz(m) * 8);
  const int rc = clab_falsctl(pd(xm), static_cast<int>(n), static_cast<int>(xc.size()), dv.data(), ym.data(), m,
                              ctl.data(), se == 0 ? 1 : 0, level, code[adj], np, sd, o.data());
  if (rc == CLAB_EDOM) return c.fail("the covariates are collinear, or `d` does not vary once they are taken out.");
  if (rc != CLAB_OK) return kerr(c, "falsctl", rc);
  for (std::size_t i = 6 * zz(m); i < o.size(); ++i) {
    if (std::isnan(o[i])) o[i] = NA;
  }
  static const char *const an[] = {"holm", "bonferroni", "hochberg", "BH", "BY", "none"};
  Frame f;
  f.str("outcome", oc);
  Cols role;
  for (int v : ctl) role.push_back(v ? "control" : "primary");
  f.str("role", role);
  static const char *const nm[] = {"est", "se", "stat", "p", "lo", "hi", "pperm", "padj"};
  double nflag = 0;
  for (int j = 0; j < 8; ++j) f.num(nm[j], o.data() + zz(j) * zz(m), zz(m));
  for (int j = 0; j < m; ++j) {
    const double v = o[7 * zz(m) + zz(j)];
    if (!std::isnan(v) && v <= alpha) nflag += 1;
  }
  out = vobj();
  put(out, "tab", f.o);
  put(out, "nflag", vnum(nflag));
  put(out, "ncontrol", vnum(static_cast<double>(nc.size())));
  put(out, "se", vstr(se == 0 ? "hc1" : "iid"));
  put(out, "adj", vstr(an[adj]));
  put(out, "alpha", vnum(alpha));
  put(out, "level", vnum(level));
  put(out, "nperm", vnum(np));
  put(out, "seed", vnum(np > 0 ? static_cast<double>(sd) : NA));
  return true;
}

// plntfx() (src/plf.cpp): R/plntfx.R in C++.
bool plntfx(Ctx &c, Val &out) {
  Val d;
  CK(dat(c, c["data"], "data", d));
  Cols yc, dc, xc;
  CK(acol(c, d, "y", {1}, false, yc));
  CK(acol(c, d, "d", {1}, false, dc));
  CK(acol(c, d, "x", {}, true, xc));
  if (dc[0] == yc[0] || has(xc, yc[0]) || has(xc, dc[0])) {
    return ckstop(c, "`y`, `d` and `x` must name different columns");
  }
  int se = 0;
  CK(marg(c, "se", {"hc1", "iid"}, se));
  std::vector<double> yv, dv, xm;
  CK(amat(c, d, yc, yv));
  CK(amat(c, d, dc, dv));
  CK(amat(c, d, xc, xm));
  const std::size_t n = nr(d);
  if (n <= xc.size() + 2) return ckstop(c, "`data` needs more rows than covariates plus two");
  bool vary = false;
  for (double v : dv) vary = vary || v != dv[0];
  if (!vary) return ckstop(c, "column `%s` must vary", dc[0].c_str());
  Nv ev, q, tv;
  CK(cknum(c, c["eff"], "eff", ev, false));
  CK(cknumd(c, "nrep", 200, q, 2, 1e6, false, true));
  const double nrep = q[0];
  double level = 0;
  CK(xopen(c, "level", 0.95, level));
  if (c.has("tol")) CK(cknum(c, c["tol"], "tol", tv, true, 0));
  unsigned int sd = 0;
  CK(xseed(c, sd));
  const int m = static_cast<int>(ev.n());
  std::vector<double> o(zz(m) * 5);
  const int rc = clab_plntfx(pd(xm), static_cast<int>(n), static_cast<int>(xc.size()), yv.data(), dv.data(),
                             ev.p(), m, static_cast<int>(nrep), se == 0 ? 1 : 0, level, sd, o.data());
  if (rc == CLAB_EDOM) return c.fail("the covariates are collinear, or `d` does not vary once they are taken out.");
  if (rc != CLAB_OK) return kerr(c, "plntfx", rc);
  std::vector<double> bias, tl;
  Val ok = varr();
  for (int e = 0; e < m; ++e) {
    bias.push_back(o[zz(e)] - ev[zz(e)]);
    tl.push_back(tv.n() ? tv[0] : 3 * o[zz(m) + zz(e)] / std::sqrt(nrep));
    ok.a.push_back(vbool(std::fabs(o[zz(e)] - ev[zz(e)]) <= tl.back()));
  }
  Frame f;
  f.num("eff", ev.v);
  f.num("mean", o.data(), zz(m));
  f.num("bias", bias);
  f.num("sd", o.data() + zz(m), zz(m));
  f.num("rmse", o.data() + 2 * zz(m), zz(m));
  f.num("cover", o.data() + 3 * zz(m), zz(m));
  f.num("rej", o.data() + 4 * zz(m), zz(m));
  f.num("tol", tl);
  f.add("ok", ok);
  out = vobj();
  put(out, "tab", f.o);
  put(out, "nrep", vnum(nrep));
  put(out, "se", vstr(se == 0 ? "hc1" : "iid"));
  put(out, "level", vnum(level));
  put(out, "seed", vnum(sd));
  return true;
}

// specrv() (src/spr.cpp): R/specrv.R in C++.
bool specrv(Ctx &c, Val &out) {
  Val d;
  CK(dat(c, c["data"], "data", d));
  Cols yc, dc, xc, oc;
  CK(acol(c, d, "y", {}, false, yc));
  CK(acol(c, d, "d", {1}, false, dc));
  CK(acol(c, d, "x", {}, true, xc));
  CK(acol(c, d, "opt", {}, true, oc));
  Cols all = yc;
  all.insert(all.end(), dc.begin(), dc.end());
  all.insert(all.end(), xc.begin(), xc.end());
  all.insert(all.end(), oc.begin(), oc.end());
  for (std::size_t i = 0; i < all.size(); ++i) {
    for (std::size_t j = 0; j < i; ++j) {
      if (all[i] == all[j]) return ckstop(c, "`y`, `d`, `x` and `opt` must name different columns");
    }
  }
  if (oc.size() > 10) return ckstop(c, "`opt` may name at most ten columns");
  int se = 0, nl = 0;
  CK(marg(c, "se", {"hc1", "iid"}, se));
  CK(marg(c, "null", {"perm", "boot"}, nl));
  Cols xo = xc;
  xo.insert(xo.end(), oc.begin(), oc.end());
  std::vector<double> ym, dv, xm;
  CK(amat(c, d, yc, ym));
  CK(amat(c, d, dc, dv));
  CK(amat(c, d, xo, xm));
  const std::size_t n = nr(d);
  if (n <= xo.size() + 2) return ckstop(c, "`data` needs more rows than covariates plus two");
  bool vary = false;
  for (double v : dv) vary = vary || v != dv[0];
  if (!vary) return ckstop(c, "column `%s` must vary", dc[0].c_str());
  const int m = static_cast<int>(oc.size()), ns = 1 << m, ny = static_cast<int>(yc.size()), K = ny * ns;
  const double dn = static_cast<double>(n);
  if (ns * dn * (static_cast<double>(xo.size()) + 1) + ns * static_cast<double>(ny) * dn > 67108864.0) {
    return ckstop(c, "too many specifications for these rows");
  }
  double level = 0;
  CK(xopen(c, "level", 0.95, level));
  Nv q;
  CK(cknumd(c, "nsim", 500, q, 1, 1e5, false, true));
  const int nsim = static_cast<int>(q[0]);
  unsigned int sd = 0;
  CK(xseed(c, sd));
  std::vector<double> sp(zz(K) * 6), st(8), dist(zz(nsim) * 3);
  const int rc = clab_specrv(pd(xm), static_cast<int>(n), static_cast<int>(xc.size()), m, ym.data(), ny, dv.data(),
                             se == 0 ? 1 : 0, level, nl, nsim, sd, sp.data(), st.data(), dist.data());
  if (rc == CLAB_EDOM) return c.fail("the covariates are collinear, or `d` does not vary once they are taken out.");
  if (rc != CLAB_OK) return kerr(c, "specrv", rc);
  Cols ctrl, outc;
  for (int b = 0; b < ns; ++b) {
    std::string s;
    for (int i = 0; i < m; ++i) {
      if ((b >> i) & 1) s += (s.empty() ? "" : " + ") + oc[zz(i)];
    }
    ctrl.push_back(s.empty() ? "none" : s);
  }
  Cols cl;
  Val sig = varr();
  for (int j = 0; j < ny; ++j) {
    for (int b = 0; b < ns; ++b) {
      outc.push_back(yc[zz(j)]);
      cl.push_back(ctrl[zz(b)]);
    }
  }
  for (int k = 0; k < K; ++k) sig.a.push_back(vbool(sp[3 * zz(K) + zz(k)] < 1 - level));
  Frame f;
  f.str("outcome", outc);
  f.str("ctrl", cl);
  static const char *const nm[] = {"est", "se", "stat", "p", "lo", "hi"};
  for (int j = 0; j < 6; ++j) f.num(nm[j], sp.data() + zz(j) * zz(K), zz(K));
  f.add("sig", sig);
  Frame t;
  t.str("stat", {"median", "nsig", "z"});
  t.num("obs", st.data(), 3);
  t.num("p", st.data() + 3, 3);
  for (double &v : dist) {
    if (std::isnan(v)) v = NA;
  }
  out = vobj();
  put(out, "spec", f.o);
  put(out, "test", t.o);
  put(out, "sign", vnum(st[6]));
  put(out, "dist", vmat(dist, nsim, 3));
  put(out, "nok", vnum(st[7]));
  put(out, "null", vstr(nl == 0 ? "perm" : "boot"));
  put(out, "nsim", vnum(nsim));
  put(out, "seed", vnum(sd));
  return true;
}

// gensw() (src/gsw.cpp): R/gensw.R in C++.
bool gensw(Ctx &c, Val &out) {
  Val d;
  CK(dat(c, c["data"], "data", d));
  Cols sc, xc;
  std::vector<double> sv, xm;
  CK(acol(c, d, "s", {1}, false, sc));
  CK(acol(c, d, "x", {}, false, xc));
  if (has(xc, sc[0])) return ckstop(c, "`s` and `x` must name different columns");
  int tg = 0;
  CK(marg(c, "target", {"generalize", "transport"}, tg));
  CK(amat(c, d, sc, sv));
  CK(amat(c, d, xc, xm));
  double n1 = 0;
  for (double v : sv) {
    if (v != 0 && v != 1) return ckstop(c, "column `%s` must hold 0 or 1", sc[0].c_str());
    n1 += v;
  }
  const std::size_t n = nr(d);
  if (n1 == 0 || n1 == static_cast<double>(n)) {
    return ckstop(c, "column `%s` must mark units in and out of the trial", sc[0].c_str());
  }
  Nv q;
  CK(cknumd(c, "trim", 0, q, 0));
  const double trim = q[0];
  if (trim >= 0.5) return ckstop(c, "`trim` must be below 0.5");
  CK(cknumd(c, "bins", 10, q, 1, 1e4, false, true));
  const int p = static_cast<int>(xc.size());
  std::vector<int> si(n);
  for (std::size_t i = 0; i < n; ++i) si[i] = static_cast<int>(sv[i]);
  std::vector<double> ps(n), w(n), cf(zz(p) + 1), st(6), bal(zz(p) * 5);
  const int rc = clab_gensw(xm.data(), static_cast<int>(n), p, si.data(), tg, trim, static_cast<int>(q[0]),
                            ps.data(), w.data(), cf.data(), st.data(), bal.data());
  if (rc == CLAB_EDOM) {
    return c.fail("the sampling model cannot be fitted: the covariates separate the trial from the rest, or are "
                  "collinear.");
  }
  if (rc != CLAB_OK) return kerr(c, "gensw", rc);
  Cols cn = {"(Intercept)"};
  cn.insert(cn.end(), xc.begin(), xc.end());
  static const char *const dn[] = {"ntrial", "ntarget", "smd", "dmean", "gindex", "ess"};
  Frame b;
  b.str("var", xc);
  static const char *const bn[] = {"target", "trial", "wtrial", "smd", "wsmd"};
  for (int j = 0; j < 5; ++j) b.num(bn[j], bal.data() + zz(j) * zz(p), zz(p));
  out = vobj();
  put(out, "ps", vnums(ps));
  put(out, "w", vnums(w));
  put(out, "coef", nvec(cn, cf.data()));
  put(out, "diag", nvec(6, dn, st.data()));
  put(out, "bal", b.o);
  put(out, "target", vstr(tg == 0 ? "generalize" : "transport"));
  return true;
}

// .trcol(): the name checks of .axcols() for one column, its values not
// read.
inline bool tcol(Ctx &c, const Val &d, const char *nm, std::string &v) {
  const Nm m = nms(c, {nm});
  if (m.null) return ckstop(c, "`%s` must name columns of `data`, got NULL", nm);
  if (!m.chr || m.na || m.v.empty()) return ckstop(c, "`%s` must name columns of `data`", nm);
  if (m.v.size() != 1) return ckstop(c, "`%s` must name 1 column(s)", nm);
  if (fc(d, m.v[0]) == nullptr) return ckstop(c, "`%s` names columns not in `data`: %s", nm, m.v[0].c_str());
  v = m.v[0];
  return true;
}

// .trval(): numbers, finite where keep, 0 elsewhere.
inline bool tval(Ctx &c, const Val &d, const std::string &v, const std::vector<double> &keep,
                 std::vector<double> &z) {
  const Val &a = *fc(d, v);
  const char *cl = ccls(a);
  if (cl[0] == 'c') return ckstop(c, "column `%s` must be numeric, got %s", v.c_str(), cl);
  z = colv(a);
  for (std::size_t i = 0; i < z.size(); ++i) {
    if (keep[i] == 1 && !std::isfinite(z[i])) return ckstop(c, "column `%s` must be finite in the trial", v.c_str());
    if (keep[i] != 1) z[i] = 0;
  }
  return true;
}

// trnsp() (src/trs.cpp): R/trnsp.R in C++.
bool trnsp(Ctx &c, Val &out) {
  Val d;
  CK(dat(c, c["data"], "data", d));
  Cols sc, xc;
  CK(acol(c, d, "s", {1}, false, sc));
  CK(acol(c, d, "x", {}, false, xc));
  std::string yn, dn;
  CK(tcol(c, d, "y", yn));
  CK(tcol(c, d, "d", dn));
  Cols all = {yn, dn, sc[0]};
  all.insert(all.end(), xc.begin(), xc.end());
  for (std::size_t i = 0; i < all.size(); ++i) {
    for (std::size_t j = 0; j < i; ++j) {
      if (all[i] == all[j]) return ckstop(c, "`y`, `d`, `s` and `x` must name different columns");
    }
  }
  int tg = 0;
  CK(marg(c, "target", {"transport", "generalize"}, tg));
  std::vector<double> sv, xm, yv, dv;
  CK(amat(c, d, sc, sv));
  CK(amat(c, d, xc, xm));
  double n1 = 0;
  for (double v : sv) {
    if (v != 0 && v != 1) return ckstop(c, "column `%s` must hold 0 or 1", sc[0].c_str());
    n1 += v;
  }
  const std::size_t n = nr(d);
  if (n1 == 0 || n1 == static_cast<double>(n)) {
    return ckstop(c, "column `%s` must mark units in and out of the trial", sc[0].c_str());
  }
  CK(tval(c, d, yn, sv, yv));
  CK(tval(c, d, dn, sv, dv));
  double a1 = 0;
  for (std::size_t i = 0; i < n; ++i) {
    if (sv[i] != 1) continue;
    if (dv[i] != 0 && dv[i] != 1) return ckstop(c, "column `%s` must hold 0 or 1 in the trial", dn.c_str());
    a1 += dv[i];
  }
  if (a1 == 0 || a1 == n1) return ckstop(c, "column `%s` must hold both arms in the trial", dn.c_str());
  Nv q;
  CK(cknumd(c, "trim", 0, q, 0));
  const double trim = q[0];
  if (trim >= 0.5) return ckstop(c, "`trim` must be below 0.5");
  CK(cknumd(c, "nboot", 200, q, 0, 1e5, false, true));
  const int nb = static_cast<int>(q[0]);
  double lv = 0;
  CK(xopen(c, "level", 0.95, lv));
  unsigned int sd = 0;
  if (nb > 0) CK(xseed(c, sd));
  std::vector<int> si(n), di(n);
  for (std::size_t i = 0; i < n; ++i) {
    si[i] = static_cast<int>(sv[i]);
    di[i] = static_cast<int>(dv[i]);
  }
  std::vector<double> est(2), se(2), dr(zz(nb) * 2), ps(n), w(n), nf(1);
  const int rc = clab_trnsp(xm.data(), static_cast<int>(n), static_cast<int>(xc.size()), yv.data(), di.data(),
                            si.data(), tg == 0 ? 1 : 0, trim, nb, sd, est.data(), se.data(), dr.data(), ps.data(),
                            w.data(), nf.data());
  if (rc == CLAB_EDOM) {
    return c.fail("a model cannot be fitted: an arm of the trial has too few rows for the covariates, or the "
                  "covariates separate the trial from the rest.");
  }
  if (rc != CLAB_OK) return kerr(c, "trnsp", rc);
  const double z = clab_qnorm(0.5 + lv / 2);
  std::vector<double> lo(2), hi(2);
  for (int a = 0; a < 2; ++a) {
    if (std::isnan(se[zz(a)])) se[zz(a)] = NA;
    lo[zz(a)] = std::isnan(se[zz(a)]) ? NA : est[zz(a)] - z * se[zz(a)];
    hi[zz(a)] = std::isnan(se[zz(a)]) ? NA : est[zz(a)] + z * se[zz(a)];
  }
  for (double &v : dr) {
    if (std::isnan(v)) v = NA;
  }
  Frame f;
  f.str("method", {"om", "iow"});
  f.num("est", est);
  f.num("se", se);
  f.num("lo", lo);
  f.num("hi", hi);
  out = vobj();
  put(out, "tab", f.o);
  put(out, "ps", vnums(ps));
  put(out, "w", vnums(w));
  put(out, "draws", vmat(dr, nb, 2));
  put(out, "nfail", vnum(nf[0]));
  put(out, "target", vstr(tg == 0 ? "transport" : "generalize"));
  put(out, "nboot", vnum(nb));
  put(out, "level", vnum(lv));
  put(out, "seed", vnum(nb > 0 ? static_cast<double>(sd) : NA));
  return true;
}

// rcttml() (src/rct.cpp): R/rcttml.R in C++.
bool rcttml(Ctx &c, Val &out) {
  Val d;
  CK(dat(c, c["data"], "data", d));
  Cols yc, dc, xc;
  CK(acol(c, d, "y", {1}, false, yc));
  CK(acol(c, d, "d", {1}, false, dc));
  CK(acol(c, d, "x", {}, true, xc));
  if (dc[0] == yc[0] || has(xc, yc[0]) || has(xc, dc[0])) {
    return ckstop(c, "`y`, `d` and `x` must name different columns");
  }
  int fam = 0;
  CK(marg(c, "fam", {"gaussian", "binomial"}, fam));
  std::vector<double> yv, dv, xm;
  CK(amat(c, d, yc, yv));
  CK(amat(c, d, dc, dv));
  CK(amat(c, d, xc, xm));
  double n1 = 0;
  for (double v : dv) {
    if (v != 0 && v != 1) return ckstop(c, "column `%s` must hold 0 or 1", dc[0].c_str());
    n1 += v;
  }
  const std::size_t n = nr(d);
  if (n1 < 2 || static_cast<double>(n) - n1 < 2) {
    return ckstop(c, "column `%s` must have two participants or more in each arm", dc[0].c_str());
  }
  if (fam == 1) {
    for (double v : yv) {
      if (v != 0 && v != 1) return ckstop(c, "with `fam = \"binomial\"`, column `%s` must hold 0 or 1", yc[0].c_str());
    }
  }
  bool inter = false, gx = false;
  CK(flagd(c, "inter", false, inter));
  CK(flagd(c, "gx", false, gx));
  double g0 = NaN;
  if (c.has("g")) {
    CK(xopen(c, "g", 0, g0));
    if (gx) return ckstop(c, "give `g` or `gx`, not both");
  }
  if (gx && xc.empty()) return ckstop(c, "`gx` needs covariates `x`");
  Nv q;
  CK(cknumd(c, "gbound", 0.01, q, 0));
  const double gb = q[0];
  if (gb >= 0.5) return ckstop(c, "`gbound` must be below 0.5");
  double lo = 0, hi = 1;
  if (fam == 0) {
    lo = *std::min_element(yv.begin(), yv.end());
    hi = *std::max_element(yv.begin(), yv.end());
    const double ylo = lo, yhi = hi;
    if (c.has("bounds")) {
      CK(cknum(c, c["bounds"], "bounds", q, false));
      if (q.n() != 2 || q[0] >= q[1] || q[0] > ylo || q[1] < yhi) {
        return ckstop(c, "`bounds` must be two increasing numbers around `y`");
      }
      lo = q[0];
      hi = q[1];
    } else if (lo >= hi) {
      return ckstop(c, "`bounds` must be two increasing numbers around `y`");
    }
  }
  CK(cknumd(c, "qbound", 0.005, q, -Inf, Inf, true));
  const double qb = q[0];
  if (qb >= 0.5) return ckstop(c, "`qbound` must be below 0.5");
  double lv = 0;
  CK(xopen(c, "level", 0.95, lv));
  std::vector<int> di(n);
  for (std::size_t i = 0; i < n; ++i) di[i] = static_cast<int>(dv[i]);
  std::vector<double> o(6), aux(3);
  const int rc = clab_rcttml(pd(xm), static_cast<int>(n), static_cast<int>(xc.size()), yv.data(), di.data(), fam,
                             inter, gx, g0, gb, lo, hi, qb, o.data(), aux.data());
  if (rc == CLAB_EDOM) {
    return c.fail("the outcome or treatment model cannot be fitted: the covariates are collinear, or separate the "
                  "outcome or the arms.");
  }
  if (rc != CLAB_OK) return kerr(c, "rcttml", rc);
  const double z = clab_qnorm(0.5 + lv / 2);
  std::vector<double> est(o.begin(), o.begin() + 3), se(o.begin() + 3, o.end()), lw(3), up(3), pv(3);
  for (int k = 0; k < 3; ++k) {
    if (std::isnan(se[zz(k)])) {
      se[zz(k)] = lw[zz(k)] = up[zz(k)] = pv[zz(k)] = NA;
      continue;
    }
    lw[zz(k)] = est[zz(k)] - z * se[zz(k)];
    up[zz(k)] = est[zz(k)] + z * se[zz(k)];
    pv[zz(k)] = 2 * clab_pnorm(-std::fabs(est[zz(k)] / se[zz(k)]));
  }
  Frame f;
  f.str("method", {"tmle", "gcomp", "unadj"});
  f.num("est", est);
  f.num("se", se);
  f.num("lo", lw);
  f.num("hi", up);
  f.num("p", pv);
  out = vobj();
  put(out, "tab", f.o);
  put(out, "eps", vnum(aux[0]));
  put(out, "g", vnums(aux.data() + 1, 2));
  put(out, "fam", vstr(fam == 0 ? "gaussian" : "binomial"));
  put(out, "bounds", vnums(std::vector<double>{lo, hi}));
  put(out, "n", vnum(static_cast<double>(n)));
  put(out, "level", vnum(lv));
  return true;
}

// rscchi() (src/rsc.cpp): src/bxpr.cpp in C++, the design read as
// svymx() reads it.
bool rscchi(Ctx &c, Val &out) {
  const std::size_t sn = vlen(c["x"]);
  const int n = static_cast<int>(sn);
  if (sn > static_cast<std::size_t>(INT_MAX)) return c.fail("`x` is too long for crimlab: more than %d values.", INT_MAX);
  if (n < 2) return c.fail("`x` must have two values at least.");
  Lab lx, ly;
  CK(codes(c, c["x"], "x", sn, lx));
  CK(codes(c, c["y"], "y", sn, ly));
  if (lx.g < 2 || ly.g < 2) return c.fail("`x` and `y` must each hold two labels or more.");
  std::vector<double> wv(sn, 1.0);
  if (c.has("w")) {
    Dv wx;
    CK(num(c, c["w"], "w", wx, n, 0.0, true));
    wv = wx.v;
  }
  std::vector<int> st(sn, 0), ps(sn);
  int nst = 1;
  if (c.has("strata")) {
    Lab l;
    CK(codes(c, c["strata"], "strata", sn, l));
    st = l.code;
    nst = l.g;
  }
  std::vector<int> pu(sn);
  int npu = n;
  if (!c.has("psu")) {
    for (int i = 0; i < n; ++i) pu[zz(i)] = i;
  } else {
    Lab l;
    CK(codes(c, c["psu"], "psu", sn, l));
    pu = l.code;
    npu = l.g;
  }
  std::vector<int> iw(zz(nst) * zz(npu));
  const int nps = pev::pairs(st.data(), pu.data(), n, nst, npu, ps.data(), iw.data());
  std::vector<double> nh(zz(nst), 0.0);
  {
    std::vector<int> sp(zz(nps), -1);
    for (std::size_t i = 0; i < sn; ++i) sp[zz(ps[i])] = st[i];
    for (int a : sp) nh[zz(a)] += 1.0;
  }
  for (double v : nh) {
    if (v < 2) return c.fail("a stratum has one first-stage unit: its variance cannot be estimated.");
  }
  std::vector<double> fv;
  if (c.has("fpc")) {
    Dv fx;
    CK(num(c, c["fpc"], "fpc", fx, n, 0.0, true));
    fv.assign(zz(nst), -1.0);
    for (std::size_t i = 0; i < sn; ++i) {
      double &f = fv[zz(st[i])];
      if (f >= 0 && f != fx[i]) return c.fail("`fpc` must be constant within each stratum.");
      f = fx[i];
    }
    for (int h = 0; h < nst; ++h) {
      if (fv[zz(h)] < nh[zz(h)]) return c.fail("`fpc` must be at least the first-stage units sampled in its stratum.");
    }
  }
  const int K = (lx.g - 1) * (ly.g - 1);
  std::vector<double> tab(zz(lx.g) * zz(ly.g)), test(16), dlt(zz(K)), aux(3);
  const int rc = clab_rscchi(lx.code.data(), ly.code.data(), lx.g, ly.g, n, wv.data(), st.data(), nst, ps.data(),
                             nps, fv.empty() ? nullptr : fv.data(), tab.data(), test.data(), dlt.data(), aux.data());
  if (rc == CLAB_EDOM) return c.fail("the design effects cannot be found: a margin of the table, or its variance, is 0.");
  if (rc != CLAB_OK) return kerr(c, "rscchi", rc);
  for (double &v : test) {
    if (std::isnan(v)) v = NA;
  }
  Frame f;
  f.str("test", {"pearson", "rs1", "rs2", "f"});
  f.num("stat", test.data(), 4);
  f.num("df1", test.data() + 4, 4);
  f.num("df2", test.data() + 8, 4);
  f.num("p", test.data() + 12, 4);
  const double dv[4] = {static_cast<double>(n), static_cast<double>(nst), static_cast<double>(nps), aux[2]};
  static const char *const dn[] = {"n", "strata", "psu", "df"};
  out = vobj();
  put(out, "tab", vmat(tab, lx.g, ly.g));
  put(out, "test", f.o);
  put(out, "deff", vnums(dlt));
  put(out, "dbar", vnum(aux[0]));
  put(out, "a2", vnum(aux[1]));
  put(out, "des", nvec(4, dn, dv));
  return true;
}

// mipmm() (src/mip.cpp): R/mipmm.R in C++. A missing value is null.
bool mipmm(Ctx &c, Val &out) {
  Val d;
  CK(dat(c, c["data"], "data", d));
  const Nm ym = nms(c, {"y"});
  if (ym.null) return ckstop(c, "`y` must name columns of `data`, got NULL");
  if (!ym.chr || ym.na || ym.v.empty()) return ckstop(c, "`y` must name columns of `data`");
  Cols miss;
  for (const std::string &s : ym.v) {
    if (fc(d, s) == nullptr && !has(miss, s)) miss.push_back(s);
  }
  if (!miss.empty()) return ckstop(c, "`y` names columns not in `data`: %s", join(miss, ", ").c_str());
  const Cols &yc = ym.v;
  for (std::size_t i = 0; i < yc.size(); ++i) {
    for (std::size_t j = 0; j < i; ++j) {
      if (yc[i] == yc[j]) return ckstop(c, "`y` names a column twice");
    }
  }
  Cols xc;
  CK(acol(c, d, "x", {}, true, xc));
  for (const std::string &s : yc) {
    if (has(xc, s)) return ckstop(c, "`y` and `x` must name different columns");
  }
  const std::size_t n = nr(d);
  std::vector<double> xm;
  for (const std::string &s : yc) {
    const Val &a = *fc(d, s);
    const char *cl = ccls(a);
    if (cl[0] == 'c') return ckstop(c, "column `%s` must be numeric, got %s", s.c_str(), cl);
    for (const Val &e : a.a) {
      if (e.t == Val::NUL) {
        xm.push_back(NaN);
        continue;
      }
      const double v = nv(e);
      if (!std::isfinite(v)) return ckstop(c, "column `%s` must be finite where it is observed", s.c_str());
      xm.push_back(v);
    }
  }
  std::vector<double> xv;
  CK(amat(c, d, xc, xv));
  xm.insert(xm.end(), xv.begin(), xv.end());
  Nv q;
  CK(cknumd(c, "m", 5, q, 2, 1000, false, true));
  const int m = static_cast<int>(q[0]);
  CK(cknumd(c, "k", 5, q, 1, 1e4, false, true));
  const int k = static_cast<int>(q[0]);
  CK(cknumd(c, "maxit", 5, q, 1, 1000, false, true));
  const int mx = static_cast<int>(q[0]);
  CK(cknumd(c, "ridge", 1e-5, q, 0));
  const double ridge = q[0];
  Cols a = yc;
  a.insert(a.end(), xc.begin(), xc.end());
  int fy = -1;
  if (c.has("fy")) {
    Cols v;
    if (!strs(c["fy"], v) || v.size() != 1 || !has(a, v[0])) {
      return ckstop(c, "`fy` must name one column of `y` or `x`");
    }
    fy = static_cast<int>(std::find(a.begin(), a.end(), v[0]) - a.begin());
  }
  Cols fxn;
  std::vector<int> fx;
  if (c.has("fx")) {
    bool ok = fy >= 0 && strs(c["fx"], fxn);
    for (std::size_t i = 0; ok && i < fxn.size(); ++i) {
      ok = has(a, fxn[i]) && fxn[i] != a[zz(fy)];
      for (std::size_t j = 0; ok && j < i; ++j) ok = fxn[i] != fxn[j];
      if (ok) fx.push_back(static_cast<int>(std::find(a.begin(), a.end(), fxn[i]) - a.begin()));
    }
    if (!ok) return ckstop(c, "`fx` must name columns of `y` or `x` other than `fy`, each once");
  }
  double lv = 0;
  CK(xopen(c, "level", 0.95, lv));
  unsigned int sd = 0;
  CK(xseed(c, sd));
  const int nt = static_cast<int>(yc.size()), qn = static_cast<int>(a.size());
  std::vector<int> tg(zz(nt));
  std::vector<double> nmis(zz(nt), 0.0);
  std::size_t tot = 0;
  for (int j = 0; j < nt; ++j) {
    tg[zz(j)] = j;
    for (std::size_t i = 0; i < n; ++i) nmis[zz(j)] += std::isnan(xm[zz(j) * n + i]) ? 1 : 0;
    tot += static_cast<std::size_t>(nmis[zz(j)]);
  }
  const int np = fy >= 0 ? static_cast<int>(fx.size()) + 1 : nt;
  std::vector<double> imp(tot * zz(m)), pool(zz(np) * 8);
  const int rc = clab_mipmm(xm.data(), static_cast<int>(n), qn, tg.data(), nt, m, k, mx, ridge, sd, fy,
                            fx.empty() ? nullptr : fx.data(), static_cast<int>(fx.size()), lv, imp.data(),
                            pool.data());
  if (rc == CLAB_EDOM) {
    return c.fail("a column of `y` has too few observed values for its predictors, or they are collinear.");
  }
  if (rc != CLAB_OK) return kerr(c, "mipmm", rc);
  Val iv = vobj();
  std::size_t at = 0;
  for (int j = 0; j < nt; ++j) {
    const int r = static_cast<int>(nmis[zz(j)]);
    put(iv, yc[zz(j)].c_str(), vmat(imp.data() + at, r, m));
    at += zz(r) * zz(m);
  }
  Cols term;
  if (fy < 0) {
    term = yc;
  } else {
    term.push_back("(Intercept)");
    term.insert(term.end(), fxn.begin(), fxn.end());
  }
  Frame f;
  f.str("term", term);
  static const char *const pn[] = {"est", "se", "df", "lo", "hi", "p", "riv", "fmi"};
  for (int j = 0; j < 8; ++j) f.num(pn[j], pool.data() + zz(j) * zz(np), zz(np));
  out = vobj();
  put(out, "imp", iv);
  put(out, "pool", f.o);
  put(out, "nmis", nvec(yc, nmis.data()));
  put(out, "m", vnum(m));
  put(out, "k", vnum(k));
  put(out, "maxit", vnum(mx));
  put(out, "seed", vnum(sd));
  return true;
}

}  // namespace ax
