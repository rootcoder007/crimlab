// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Events in time (src/btm.cpp): the Hawkes fit and its marked,
// multivariate and spatio-temporal kin, trends in counts, CUSUM
// surveillance, cross-correlation, survival of recurrences, the
// hypercube queue and the burglary model.

// picks(): each text of a character vector among the options.
bool picks(Ctx &c, const char *nm, const char *const *opt, int n, int def, std::vector<int> &out) {
  out.clear();
  if (!c.has(nm)) {
    out.push_back(def);
    return true;
  }
  const Val &v = c[nm];
  std::vector<Val> e;
  if (v.t == Val::STR) e.push_back(v);
  else if (v.t == Val::ARR && !v.a.empty()) e = v.a;
  bool ok = !e.empty();
  for (const Val &x : e) ok = ok && x.t == Val::STR;
  if (!ok) e.assign(1, Val());
  for (const Val &x : e) {
    int k = 0;
    CK(pick(c, x, nm, opt, n, k));
    out.push_back(k);
  }
  return true;
}

std::vector<std::string> hpnames(int kind, int base) {
  std::vector<std::string> s;
  if (base) s = {"a0", "a1", "a2", "a3"};
  else s = {"a0"};
  s.emplace_back("eta");
  if (kind) {
    s.emplace_back("shape");
    s.emplace_back("scale");
  } else {
    s.emplace_back("beta");
  }
  return s;
}

// evtimes() on given values: finite, sorted, in [0, hor], at least nmin.
bool evtv(Ctx &c, const std::vector<double> &t, double hor, int nmin) {
  const int n = static_cast<int>(t.size());
  if (n < nmin) return c.fail("`times` must hold at least %d events, not %d.", nmin, n);
  for (std::size_t i = 0; i < t.size(); ++i) {
    if (t[i] < 0 || t[i] > hor) return c.fail("`times` must lie between 0 and `hor`.");
    if (i > 0 && t[i] < t[i - 1]) return c.fail("`times` must be sorted from earliest to latest.");
  }
  return true;
}

void tiewarn(Ctx &c, const std::vector<double> &t, const char *fn) {
  for (std::size_t i = 1; i < t.size(); ++i) {
    if (t[i] == t[i - 1]) {
      c.warn.push_back(std::string(fn) +
                       ": some times are tied; events at one instant do not excite each other. Times recorded to a "
                       "resolution (days, say) are better spread across it first (`res`).");
      return;
    }
  }
}

// A covariance matrix as rows, NaN as NA.
Val covm(const double *cov, int np) {
  std::vector<double> v(cov, cov + static_cast<std::size_t>(np) * np);
  for (double &e : v) e = na(e);
  return vmat(v, np, np);
}
// par, se and cov of a fit, named by nm.
void pse(Val &out, const std::vector<std::string> &nm, const double *th, const double *cov) {
  const int np = static_cast<int>(nm.size());
  std::vector<double> p(th, th + np), s(static_cast<std::size_t>(np));
  for (int j = 0; j < np; ++j) s[static_cast<std::size_t>(j)] = na(std::sqrt(cov[j * np + j]));
  for (double &e : p) e = na(e);
  put(out, "par", nvec(nm, p.data()));
  put(out, "se", nvec(nm, s.data()));
  put(out, "cov", covm(cov, np));
}

bool hwkfit(Ctx &c, Val &out) {
  double h = 0.0, pr = 0.0;
  CK(pos1(c, "hor", NaN, h));
  CK(pos1(c, "per", 365.25, pr));
  int mx = 0;
  CK(wholed(c, "maxit", 500, 1, mx));
  std::vector<int> ks, bs;
  CK(picks(c, "kernel", hkinds, 3, 0, ks));
  CK(picks(c, "base", hbases, 2, 0, bs));
  std::vector<double> t;
  bool jit = false;
  if (c.has("res")) {
    double r = 0.0;
    CK(pos1(c, "res", NaN, r));
    Dv tv;
    CK(num(c, c["times"], "times", tv));
    t = tv.v;
    unsigned int s = 0;
    CK(seed(c, s));
    kutil::rng g(s);
    for (double &e : t) {
      if (e < 0 || e > h) return c.fail("`times` must lie between 0 and `hor`.");
      e += g.unif() * std::fmin(r, h - e);
    }
    kutil::hsort(t.data(), static_cast<int>(t.size()), [](double a, double b) { return a < b; });
    CK(evtv(c, t, h, 10));
    jit = true;
  } else {
    CK(evtimes(c, h, 10, t));
    tiewarn(c, t, "hwkfit");
  }
  const int n = static_cast<int>(t.size());
  const int nc = static_cast<int>(ks.size() * bs.size());
  std::vector<double> st0;
  if (c.has("start")) {
    if (nc > 1) return c.fail("`start` needs one kernel and one baseline.");
    CK(parv(c, "start", ks[0], bs[0], st0));
    if (!(st0[static_cast<std::size_t>(hnpar(ks[0], bs[0]) - (ks[0] ? 3 : 2))] < 1)) {
      return c.fail("`start`: eta must be below 1.");
    }
  }
  std::vector<std::string> ck, cb;
  std::vector<double> cl, ckp, ca, cbi, ce;
  std::vector<int> cv;
  int bi = -1, bk = 0, bb = 0;
  std::vector<double> bth, bcov, bst;
  for (int q = 0; q < nc; ++q) {
    const int kd = ks[static_cast<std::size_t>(q) / bs.size()], bd = bs[static_cast<std::size_t>(q) % bs.size()];
    const int np = hnpar(kd, bd);
    std::vector<double> th(static_cast<std::size_t>(np)), cov(static_cast<std::size_t>(np * np)), st(3);
    if (clab_hwkfit(t.data(), n, h, kd, bd, pr, st0.empty() ? nullptr : st0.data(), mx, th.data(), cov.data(),
                    st.data()) != CLAB_OK) {
      return c.fail("the %s fit failed: the likelihood is not finite at any start.", hkinds[kd]);
    }
    ck.push_back(hkinds[kd]);
    cb.push_back(hbases[bd]);
    cl.push_back(st[0]);
    ckp.push_back(np);
    ca.push_back(2 * np - 2 * st[0]);
    cbi.push_back(np * std::log(static_cast<double>(n)) - 2 * st[0]);
    ce.push_back(th[static_cast<std::size_t>(np - (kd ? 3 : 2))]);
    cv.push_back(st[2] > 0);
    if (bi < 0 || ca.back() < ca[static_cast<std::size_t>(bi)]) {
      bi = q;
      bth = th;
      bcov = cov;
      bst = st;
      bk = kd;
      bb = bd;
    }
  }
  const int np = hnpar(bk, bb);
  const double eta = bth[static_cast<std::size_t>(np - (bk ? 3 : 2))];
  const double br[4] = {eta, eta < 1 ? 1 / (1 - eta) : Inf, bb == 0 && eta < 1 ? std::exp(bth[0]) / (1 - eta) : NA, eta};
  static const char *const bn[] = {"eta", "size", "rate", "share"};
  if (eta > 0.99) {
    char m[256];
    std::snprintf(m, sizeof m,
                  "hwkfit: the branching ratio is %.3f, near 1: the fit is near-critical and probably mistakes a trend "
                  "or a cycle for contagion (try base = \"cycle\").",
                  eta);
    c.warn.push_back(m);
  }
  Frame cmp;
  cmp.str("kernel", ck);
  cmp.str("base", cb);
  cmp.num("loglik", cl);
  cmp.num("k", ckp);
  cmp.num("aic", ca);
  cmp.num("bic", cbi);
  cmp.num("eta", ce);
  cmp.lgl("conv", cv.data(), cv.size());
  out = vobj();
  pse(out, hpnames(bk, bb), bth.data(), bcov.data());
  put(out, "loglik", vnum(bst[0]));
  put(out, "aic", vnum(ca[static_cast<std::size_t>(bi)]));
  put(out, "bic", vnum(cbi[static_cast<std::size_t>(bi)]));
  put(out, "branch", nvec(4, bn, br));
  put(out, "conv", vbool(bst[2] > 0));
  put(out, "iter", vnum(static_cast<int>(bst[1])));
  put(out, "kernel", vstr(hkinds[bk]));
  put(out, "base", vstr(hbases[bb]));
  put(out, "per", vnum(pr));
  put(out, "n", vnum(n));
  put(out, "hor", vnum(h));
  put(out, "jit", vbool(jit));
  std::vector<double> tt = t;
  for (double &e : tt) e = na(e);
  put(out, "times", vauto(tt));
  put(out, "cmp", cmp.o);
  return true;
}

// periods(): 1..n, or n finite numbers, strictly increasing.
bool periods(Ctx &c, std::size_t n, std::vector<double> &v) {
  v.resize(n);
  if (!c.has("x")) {
    for (std::size_t i = 0; i < n; ++i) v[i] = static_cast<double>(i + 1);
    return true;
  }
  Dv xv;
  CK(num(c, c["x"], "x", xv, static_cast<long>(n)));
  for (std::size_t i = 0; i < n; ++i) {
    v[i] = xv[i];
    if (i > 0 && !(v[i] > v[i - 1])) return c.fail("`x` must be strictly increasing.");
  }
  return true;
}
// One value as a one-row data frame column.
inline void col1(Frame &f, const char *nm, double v) {
  const double x = na(v);
  f.num(nm, &x, 1);
}
inline void lgl1(Frame &f, const char *nm, double v) {
  const int b = v != 0;
  f.lgl(nm, &b, 1);
}

const char *const talts[] = {"both", "up", "down"};

bool trndrr(Ctx &c, Val &out) {
  Dv yv;
  CK(num(c, c["y"], "y", yv, -1, 0.0));
  const std::size_t sn = yv.n();
  const int n = static_cast<int>(sn);
  if (n < 3) return c.fail("`y` must hold at least 3 periods, not %d.", n);
  const double *yy = yv.p();
  std::vector<double> xv;
  CK(periods(c, sn, xv));
  int a = 0;
  CK(pickd(c, "alt", talts, 3, 0, a));
  double z = 0.0, lv = 0.0;
  CK(zconf(c, "level", 0.95, z));
  lv = val1d(c, "level", 0.95);
  std::vector<double> off(sn, 0.0);
  if (c.has("expo")) {
    Dv e;
    CK(num(c, c["expo"], "expo", e, n, 0.0, true));
    for (std::size_t i = 0; i < sn; ++i) off[i] = std::log(e[i]);
  }
  double xm = 0, tot = 0;
  for (std::size_t i = 0; i < sn; ++i) {
    xm += xv[i] / n;
    tot += yy[i];
  }
  if (!(tot > 0)) return c.fail("`y` has no events: every count is 0.");
  int lo = n, hi = -1;
  for (int i = 0; i < n; ++i) {
    if (yy[i] > 0) {
      lo = std::min(lo, i);
      hi = std::max(hi, i);
    }
  }
  std::vector<double> d(2 * sn), mu(sn, NA);
  for (std::size_t i = 0; i < sn; ++i) {
    d[i] = 1;
    d[sn + i] = xv[i] - xm;
  }
  double rr[8];
  for (double &v : rr) v = NA;
  rr[7] = 0;
  if (lo == hi && (lo == 0 || hi == n - 1)) {
    rr[0] = lo == 0 ? 0.0 : Inf;
    c.warn.push_back(std::string("trndrr: every event falls in the ") + (lo == 0 ? "first" : "last") +
                     " period: the trend has no finite estimate.");
  } else {
    double beta[2] = {0}, vb[4] = {0}, vc[4] = {0}, st[8] = {0};
    int al[2] = {0}, iw[2] = {0};
    std::vector<double> w(2 * sn + 6 * sn + 2 + 4);
    if (clab_cntreg(yy, d.data(), n, 2, off.data(), nullptr, 0, 0, 25, 1e-8, beta, vb, vc, st, al, w.data(), iw) !=
        CLAB_OK) {
      return c.fail("the Poisson fit failed: the deviance left the finite range.");
    }
    double pe = 0;
    for (std::size_t i = 0; i < sn; ++i) {
      const double m = std::exp(beta[0] + beta[1] * d[sn + i] + off[i]);
      mu[i] = m;
      pe += (yy[i] - m) * (yy[i] - m) / m;
    }
    const double disp = pe / (n - 2), b = beta[1];
    double se = std::sqrt(vb[3]), q = z, pv = 0.0;
    const bool over = disp > 1.5;
    if (over) {
      se *= std::sqrt(disp);
      q = clab_qt(0.5 + 0.5 * lv, n - 2, 1);
      pv = 2 * clab_pt(-std::fabs(b / se), n - 2, 1);
    } else {
      pv = 2 * clab_pnorm(-std::fabs(b / se));
    }
    const double r0[8] = {std::exp(b), std::exp(b - q * se), std::exp(b + q * se), pv, b, se, disp, over ? 1.0 : 0.0};
    for (int j = 0; j < 8; ++j) rr[j] = r0[j];
  }
  Frame rrf;
  static const char *const rn[] = {"rr", "lo", "hi", "p", "b", "se", "disp"};
  for (int j = 0; j < 7; ++j) col1(rrf, rn[j], rr[j]);
  lgl1(rrf, "over", rr[7]);
  double st[10] = {0};
  std::vector<double> w(sn * (sn + 1) / 2);
  std::vector<int> iw(sn);
  clab_trndrr(xv.data(), yy, n, a, z, st, w.data(), iw.data());
  Frame mk, sen;
  col1(mk, "s", st[0]);
  col1(mk, "var", st[1]);
  col1(mk, "tau", st[2]);
  col1(mk, "p", st[3]);
  lgl1(mk, "exact", n <= 8);
  col1(sen, "slope", st[4]);
  col1(sen, "int", st[5]);
  col1(sen, "lo", st[6]);
  col1(sen, "hi", st[7]);
  lgl1(sen, "clamp", st[8] + st[9]);
  std::vector<double> chg(sn, NA), pct(sn, NA);
  for (std::size_t i = 1; i < sn; ++i) {
    chg[i] = yy[i] - yy[i - 1];
    pct[i] = yy[i - 1] > 0 ? 100 * (yy[i] - yy[i - 1]) / yy[i - 1] : NA;
  }
  Frame yoy;
  yoy.num("x", xv);
  yoy.num("y", yv.v);
  yoy.num("chg", chg);
  yoy.num("pct", pct);
  Val sv;
  if (c.has("seas")) {
    Lab l;
    CK(codes(c, c["seas"], "seas", sn, l));
    if (l.g < 2) return c.fail("`seas` has one season: the test needs two or more.");
    const std::size_t sg = static_cast<std::size_t>(l.g);
    std::vector<double> nn(sg, 0.0), ex(sg, tot / l.g);
    for (std::size_t i = 0; i < sn; ++i) nn[static_cast<std::size_t>(l.code[i])] += yy[i];
    double x2 = 0;
    for (std::size_t g = 0; g < sg; ++g) x2 += (nn[g] - ex[g]) * (nn[g] - ex[g]) / ex[g];
    Frame tab;
    tab.str("season", l.lev);
    tab.num("n", nn);
    tab.num("expd", ex);
    sv = vobj();
    put(sv, "tab", tab.o);
    put(sv, "stat", vnum(x2));
    put(sv, "df", vnum(l.g - 1));
    put(sv, "p", vnum(clab_pchisq(x2, l.g - 1, 0)));
  }
  for (double &e : mu) e = na(e);
  out = vobj();
  put(out, "rr", rrf.o);
  put(out, "mk", mk.o);
  put(out, "sen", sen.o);
  put(out, "yoy", yoy.o);
  put(out, "seas", sv);
  put(out, "fitted", vauto(mu));
  put(out, "n", vnum(n));
  put(out, "alt", vstr(talts[a]));
  put(out, "level", vnum(lv));
  return true;
}

bool cusrv(Ctx &c, Val &out) {
  static const char *const tr[] = {"std", "rossi", "anscombe", "none"};
  Dv yv;
  CK(num(c, c["y"], "y", yv, -1, 0.0));
  const std::size_t sn = yv.n();
  const int n = static_cast<int>(sn);
  Dv kd;
  Val d1;
  CK(num(c, arg(c, "k", d1, 1.04), "k", kd, 1));
  double hh = 0.0;
  CK(pos1(c, "h", 2.26, hh));
  int t = 0;
  CK(pickd(c, "trans", tr, 4, 0, t));
  bool rs = false;
  CK(flagd(c, "reset", false, rs));
  int s0 = 0;
  std::vector<double> m(sn);
  if (!c.has("expd")) {
    if (!c.has("start")) return c.fail("give `expd`, or `start`, the first period watched.");
    CK(whole(c, c["start"], "start", 2, s0));
    --s0;
    if (s0 >= n) return c.fail("`start` must be at most the number of counts, %d.", n);
    double b = 0;
    for (int i = 0; i < s0; ++i) b += yv[static_cast<std::size_t>(i)] / s0;
    if (!(b > 0)) return c.fail("the counts before `start` are all 0: no expected count.");
    for (double &v : m) v = b;
  } else {
    Dv e;
    CK(num(c, c["expd"], "expd", e, n, 0.0, true));
    m = e.v;
  }
  const std::size_t nw = sn - static_cast<std::size_t>(s0);
  std::vector<double> z(nw), s(nw);
  std::vector<int> al(nw);
  clab_cusrv(yv.p() + s0, m.data() + s0, static_cast<int>(nw), kd[0], hh, t, rs ? 1 : 0, z.data(), s.data(), al.data());
  std::vector<int> ix(nw);
  for (std::size_t i = 0; i < nw; ++i) ix[i] = s0 + static_cast<int>(i) + 1;
  std::vector<double> yy(yv.v.begin() + s0, yv.v.end()), mm(m.begin() + s0, m.end());
  for (double &e : z) e = na(e);
  for (double &e : s) e = na(e);
  Frame f;
  f.ints("i", ix.data(), nw);
  f.num("y", yy);
  f.num("expd", mm);
  f.num("z", z);
  f.num("s", s);
  f.lgl("alarm", al.data(), nw);
  out = f.o;
  return true;
}

bool xcser(Ctx &c, Val &out) {
  Dv xv, yv;
  CK(num(c, c["x"], "x", xv));
  CK(num(c, c["y"], "y", yv));
  std::vector<double> a, b;
  Val common;
  if (!xv.nm.empty() && !yv.nm.empty()) {
    std::vector<std::string> cm;
    for (std::size_t i = 0; i < xv.n(); ++i) {
      for (std::size_t j = 0; j < yv.n(); ++j) {
        if (xv.nm[i] == yv.nm[j]) {
          a.push_back(xv[i]);
          b.push_back(yv[j]);
          cm.push_back(xv.nm[i]);
          break;
        }
      }
    }
    common = unbox(vstrs(cm));
  } else {
    if (xv.n() != yv.n()) return c.fail("`x` and `y` must have one length, or names to match them by.");
    a = xv.v;
    b = yv.v;
  }
  bool dif = false;
  CK(flagd(c, "dif", false, dif));
  if (dif) {
    for (std::size_t i = a.empty() ? 0 : a.size() - 1; i > 0; --i) {
      a[i] -= a[i - 1];
      b[i] -= b[i - 1];
    }
    if (!a.empty()) {
      a.erase(a.begin());
      b.erase(b.begin());
    }
  }
  const int n = static_cast<int>(a.size());
  if (n < 3) return c.fail("the series have %d periods in common: at least 3 are needed.", n);
  double z = 0.0;
  CK(zconf(c, "level", 0.95, z));
  int lag = std::max(0, std::min(n - 1, static_cast<int>(std::floor(10 * std::log10(n / 2.0)))));
  if (c.has("lagmax")) CK(whole(c, c["lagmax"], "lagmax", 0, lag));
  if (lag >= n) return c.fail("`lagmax` must be below the number of periods, %d.", n);
  std::vector<double> r(static_cast<std::size_t>(2 * lag + 1));
  if (clab_xcser(a.data(), b.data(), n, lag, r.data()) != CLAB_OK) return c.fail("a series is constant: no correlation.");
  const double r0 = r[static_cast<std::size_t>(lag)];
  const double fz = std::atanh(r0), sz = n > 3 ? 1 / std::sqrt(n - 3.0) : NA;
  const double tt = r0 * std::sqrt((n - 2) / (1 - r0 * r0));
  const double cc[4] = {r0, n > 3 ? std::tanh(fz - z * sz) : NA, n > 3 ? std::tanh(fz + z * sz) : NA,
                        std::fabs(r0) < 1 ? 2 * clab_pt(-std::fabs(tt), n - 2, 1) : 0.0};
  std::vector<int> lg;
  for (int k = -lag; k <= lag; ++k) lg.push_back(k);
  for (double &e : r) e = na(e);
  Frame cor, ccf;
  col1(cor, "r", cc[0]);
  col1(cor, "lo", cc[1]);
  col1(cor, "hi", cc[2]);
  col1(cor, "p", cc[3]);
  const int nn = n;
  cor.ints("n", &nn, 1);
  ccf.ints("lag", lg.data(), lg.size());
  ccf.num("r", r);
  out = vobj();
  put(out, "cor", cor.o);
  put(out, "ccf", ccf.o);
  put(out, "common", common);
  put(out, "level", vnum(val1d(c, "level", 0.95)));
  return true;
}

bool kmrec(Ctx &c, Val &out) {
  Dv tv;
  CK(num(c, c["time"], "time", tv, -1, 0.0));
  const std::size_t s0 = tv.n();
  const int n0 = static_cast<int>(s0);
  double z = 0.0;
  CK(zconf(c, "level", 0.95, z));
  std::vector<double> t;
  std::vector<int> ev, gr;
  int g = 1;
  Lab gl;
  const bool hg = c.has("grp");
  if (hg) {
    CK(codes(c, c["grp"], "grp", s0, gl));
    g = gl.g;
  }
  if (!c.has("unit")) {
    t = tv.v;
    if (!c.has("event")) {
      ev.assign(s0, 1);
    } else {
      CK(bin01(c, c["event"], "event", s0, ev));
    }
    for (std::size_t i = 0; i < s0; ++i) gr.push_back(hg ? gl.code[i] : 0);
  } else {
    if (c.has("event")) return c.fail("`event` is not read with `unit`: the gaps are the events.");
    Lab u;
    CK(codes(c, c["unit"], "unit", s0, u));
    double e = Inf;
    if (c.has("end")) {
      Dv ed;
      CK(num(c, c["end"], "end", ed, 1));
      e = ed[0];
      for (double v : tv.v) {
        if (v > e) return c.fail("`end` must be at or after every time.");
      }
    }
    std::vector<int> ix(s0);
    for (int i = 0; i < n0; ++i) ix[static_cast<std::size_t>(i)] = i;
    const int *uc = u.code.data();
    const double *tp = tv.p();
    kutil::hsort(ix.data(), n0, [&](int a, int b) { return uc[a] == uc[b] ? (tp[a] == tp[b] ? a < b : tp[a] < tp[b]) : uc[a] < uc[b]; });
    for (int i = 0; i < n0; ++i) {
      const int a = ix[static_cast<std::size_t>(i)];
      const int gg = hg ? gl.code[static_cast<std::size_t>(a)] : 0;
      if (i > 0 && uc[ix[static_cast<std::size_t>(i - 1)]] == uc[a]) {
        const int b = ix[static_cast<std::size_t>(i - 1)];
        if (hg && gl.code[static_cast<std::size_t>(b)] != gg) return c.fail("`grp` must be the same for every time of a unit.");
        t.push_back(tp[a] - tp[b]);
        ev.push_back(1);
        gr.push_back(gg);
      }
      const bool last = i == n0 - 1 || uc[ix[static_cast<std::size_t>(i + 1)]] != uc[a];
      if (last && std::isfinite(e)) {
        t.push_back(e - tp[a]);
        ev.push_back(0);
        gr.push_back(gg);
      }
    }
    if (t.empty()) return c.fail("no unit has two times and `end` is not given: there are no gaps.");
  }
  const std::vector<std::string> lev = hg ? gl.lev : std::vector<std::string>{"all"};
  const int n = static_cast<int>(t.size());
  const std::size_t sn = t.size(), sg = static_cast<std::size_t>(g);
  std::vector<int> go(sn), iw(sn);
  std::vector<double> tt(sn), nr(sn), nd(sn), nc(sn), s(sn), se(sn), lo(sn), hi(sn), med(2 * sg), w(3 * sg + sg * sg);
  double lr[2] = {0};
  int m = 0;
  clab_kmrec(t.data(), ev.data(), gr.data(), n, g, z, &m, go.data(), tt.data(), nr.data(), nd.data(), nc.data(),
             s.data(), se.data(), lo.data(), hi.data(), med.data(), lr, w.data(), iw.data());
  const std::size_t sm = static_cast<std::size_t>(m);
  std::vector<std::string> gk(sm);
  for (std::size_t i = 0; i < sm; ++i) gk[i] = lev[static_cast<std::size_t>(go[i])];
  auto nv = [](const std::vector<double> &v, std::size_t k) {
    std::vector<double> r(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(k));
    for (double &e : r) e = na(e);
    return r;
  };
  Frame km;
  km.str("grp", gk);
  km.num("time", nv(tt, sm));
  km.num("nrisk", nv(nr, sm));
  km.num("nevent", nv(nd, sm));
  km.num("ncens", nv(nc, sm));
  km.num("surv", nv(s, sm));
  km.num("se", nv(se, sm));
  km.num("lo", nv(lo, sm));
  km.num("hi", nv(hi, sm));
  std::vector<double> gn(sg, 0.0), ge(sg, 0.0);
  for (std::size_t i = 0; i < sn; ++i) {
    gn[static_cast<std::size_t>(gr[i])] += 1;
    ge[static_cast<std::size_t>(gr[i])] += ev[i];
  }
  Frame md;
  md.str("grp", lev);
  md.num("n", gn);
  md.num("events", ge);
  md.num("median", nv(med, sg));
  Val lrv;
  if (g > 1) {
    const double pv = std::isnan(lr[0]) ? NA : clab_pchisq(lr[0], lr[1], 0);
    lrv = vobj();
    put(lrv, "stat", vnum(na(lr[0])));
    put(lrv, "df", vnum(g - 1));
    put(lrv, "p", vnum(pv));
  }
  out = vobj();
  put(out, "km", km.o);
  put(out, "med", md.o);
  put(out, "lr", lrv);
  put(out, "n", vnum(n));
  put(out, "level", vnum(val1d(c, "level", 0.95)));
  return true;
}

// numat(): a numeric matrix's values, finite, at least lo, column-major.
bool numat(Ctx &c, const Val &x, const char *nm, double lo, int &nr, int &nc, std::vector<double> &v) {
  int dim[2] = {-1, -1};
  if (!(x.t == Val::ARR && !x.a.empty() && x.a[0].t == Val::ARR && shape(x, 2, dim, 0))) return c.fail("`%s` must hold numbers.", nm);
  nr = dim[0];
  nc = dim[1];
  v.assign(static_cast<std::size_t>(nr) * static_cast<std::size_t>(nc), 0.0);
  bool bad = false;
  fill(x, 2, dim, 0, 0, 1, v, bad, true);
  for (double e : v) {
    if (!std::isfinite(e) || e < lo) return c.fail("`%s` must hold finite numbers from %g.", nm, lo);
  }
  return true;
}
inline bool ismat(const Val &x, int *dim) {
  dim[0] = dim[1] = -1;
  return x.t == Val::ARR && !x.a.empty() && x.a[0].t == Val::ARR && shape(x, 2, dim, 0);
}

bool hypq(Ctx &c, Val &out) {
  static const char *const qs[] = {"loss", "wait"};
  Dv lv;
  CK(num(c, c["lam"], "lam", lv, -1, 0.0, true));
  const int m = static_cast<int>(lv.n());
  double u = 0.0;
  CK(pos1(c, "mu", 1, u));
  int w = 0;
  CK(pickd(c, "queue", qs, 2, 0, w));
  const Val &pf = c["pref"];
  int n = -1;
  std::vector<int> pr;
  int dim[2];
  if (pf.t == Val::OBJ || (pf.t == Val::ARR && !pf.a.empty() && pf.a[0].t == Val::ARR && !ismat(pf, dim))) {
    if (static_cast<int>(pf.a.size()) != m) return c.fail("`pref` must hold one order per atom, %d.", m);
    n = static_cast<int>(vlen(pf.a[0]));
    pr.assign(static_cast<std::size_t>(m * n), 0);
    for (int j = 0; j < m; ++j) {
      const Val &v = pf.a[static_cast<std::size_t>(j)];
      if (static_cast<int>(vlen(v)) != n) return c.fail("`pref`: every order must list all %d units.", n);
      Dv d;
      CK(num(c, v, "pref", d));
      for (int r = 0; r < n; ++r) pr[static_cast<std::size_t>(j + r * m)] = static_cast<int>(d[static_cast<std::size_t>(r)]) - 1;
    }
  } else {
    if (!ismat(pf, dim) || dim[0] != m) return c.fail("`pref` must be a matrix with one row per atom, %d, or a list.", m);
    int nr = 0;
    std::vector<double> d;
    CK(numat(c, pf, "pref", -HUGE_VAL, nr, n, d));
    pr.resize(d.size());
    for (std::size_t i = 0; i < d.size(); ++i) pr[i] = static_cast<int>(d[i]) - 1;
  }
  for (int e : pr) {
    if (e < 0 || e >= n) return c.fail("`pref` must hold unit numbers from 1 to %d.", n);
  }
  if (n > 10) return c.fail("at most 10 units: the model has 2^n states.");
  const std::size_t s = static_cast<std::size_t>(1) << n;
  std::vector<double> pst(s), wl(static_cast<std::size_t>(n)), dp(static_cast<std::size_t>(n * m)), work(s * s + s);
  double st[3] = {0};
  if (clab_hypq(lv.p(), pr.data(), m, n, u, w, pst.data(), wl.data(), dp.data(), st, work.data()) != CLAB_OK) {
    return c.fail(w ? "every order must list each unit once, and the calls must be fewer than the units can serve "
                      "(rho below 1) for a queue."
                    : "every order must list each unit once.");
  }
  std::vector<std::string> bs(s);
  for (std::size_t b = 0; b < s; ++b) {
    for (int v = 0; v < n; ++v) bs[b] += ((b >> v) & 1u) ? '1' : '0';
  }
  double big = 0;
  for (double e : lv.v) big += e;
  Val trv;
  if (c.has("tt")) {
    const Val &tv0 = c["tt"];
    if (!ismat(tv0, dim) || dim[0] != n || dim[1] != m) {
      return c.fail("`tt` must be a matrix of travel times, one row per unit and one column per atom.");
    }
    int nr = 0, nc2 = 0;
    std::vector<double> tv;
    CK(numat(c, tv0, "tt", 0.0, nr, nc2, tv));
    std::vector<double> at(static_cast<std::size_t>(m));
    double nm2 = 0, den = 0;
    for (int j = 0; j < m; ++j) {
      double a = 0, b = 0;
      for (int v = 0; v < n; ++v) {
        a += dp[static_cast<std::size_t>(v + j * n)] * tv[static_cast<std::size_t>(v + j * n)];
        b += dp[static_cast<std::size_t>(v + j * n)];
      }
      at[static_cast<std::size_t>(j)] = na(a / b);
      nm2 += a;
      den += b;
    }
    trv = vobj();
    put(trv, "mean", vnum(nm2 / den));
    put(trv, "atom", vauto(at));
  }
  for (double &e : pst) e = na(e);
  for (double &e : wl) e = na(e);
  Frame sf;
  sf.str("busy", bs);
  sf.num("p", pst);
  out = vobj();
  put(out, "state", sf.o);
  put(out, "wl", vauto(wl));
  put(out, "disp", vmat(dp, n, m));
  put(out, "queue", vstr(qs[w]));
  put(out, "pfull", vnum(st[0]));
  put(out, "lq", vnum(st[1]));
  put(out, "wq", vnum(st[1] / big));
  put(out, "rho", vnum(st[2]));
  put(out, "travel", trv);
  return true;
}

// field(): an n x n matrix of numbers from 0, or fill.
bool field(Ctx &c, const char *nm, double fl, int nn, std::vector<double> &v) {
  const std::size_t s = static_cast<std::size_t>(nn) * static_cast<std::size_t>(nn);
  v.assign(s, fl);
  if (!c.has(nm)) return true;
  int dim[2];
  if (!ismat(c[nm], dim) || dim[0] != nn || dim[1] != nn) return c.fail("`%s` must be an n x n matrix.", nm);
  int nr = 0, nc = 0;
  return numat(c, c[nm], nm, 0.0, nr, nc, v);
}
// One number from lo (num(x, nm, 1, lo)).
bool num1(Ctx &c, const char *nm, double def, double lo, double &v) {
  Dv d;
  Val d1;
  CK(num(c, arg(c, nm, d1, def), nm, d, 1, lo));
  v = d[0];
  return true;
}

bool sdbsim(Ctx &c, Val &out) {
  static const char *const md[] = {"agent", "pde"};
  int nn = 0, ns = 0, mo = 0;
  CK(whole(c, c["n"], "n", 3, nn));
  CK(whole(c, c["steps"], "steps", 0, ns));
  CK(pickd(c, "mode", md, 2, 0, mo));
  if (nn > 2000) return c.fail("`n` must be at most 2000.");
  double par[10] = {0};
  CK(pos1(c, "a0", 1.0 / 30, par[0]));
  CK(num1(c, "omega", 1.0 / 15, 0.0, par[1]));
  CK(p01d(c, "eta", 0.03, par[2]));
  CK(num1(c, "gam", 0.002, 0.0, par[3]));
  CK(num1(c, "theta", 1, 0.0, par[4]));
  CK(pos1(c, "dt", mo ? 0.01 : 1, par[5]));
  CK(pos1(c, "dx", 1, par[6]));
  CK(pos1(c, "dif", 1, par[7]));
  CK(num1(c, "eps", 1, 0.0, par[8]));
  CK(pos1(c, "z", 4, par[9]));
  if (par[1] * par[5] > 1) return c.fail("`omega` times `dt` must be at most 1.");
  const std::size_t s = static_cast<std::size_t>(nn) * static_cast<std::size_t>(nn);
  std::vector<double> bv, rv;
  CK(field(c, "b", 0.0, nn, bv));
  CK(field(c, "rho", mo ? par[3] / par[0] : 0.0, nn, rv));
  if (mo == 0) {
    for (double v : rv) {
      if (v != std::floor(v) || v > 1e6) return c.fail("`rho` must hold whole numbers of offenders, at most 1e6 a site.");
    }
  }
  std::vector<int> sp;
  if (c.has("sup")) {
    int dim[2];
    if (!ismat(c["sup"], dim) || dim[0] != nn || dim[1] != nn) return c.fail("`sup` must be an n x n matrix.");
    sp.assign(s, 0);
    for (int i = 0; i < nn; ++i) {
      for (int j = 0; j < nn; ++j) {
        const Val &e = c["sup"].a[static_cast<std::size_t>(i)].a[static_cast<std::size_t>(j)];
        if (e.t != Val::BOOL) return c.fail("`sup` must be TRUE or FALSE at each site, with no missing values.");
        sp[static_cast<std::size_t>(i) + static_cast<std::size_t>(nn) * static_cast<std::size_t>(j)] = e.b;
      }
    }
  }
  unsigned int sd = 0;
  CK(seed(c, sd));
  std::vector<double> burg(static_cast<std::size_t>(ns)), mb(static_cast<std::size_t>(ns)), tot(s), w(3 * s);
  std::vector<int> iw(s);
  if (clab_sdbsim(nn, ns, mo, par, bv.data(), rv.data(), sp.empty() ? nullptr : sp.data(), sd, burg.data(), mb.data(),
                  tot.data(), w.data(), iw.data()) != CLAB_OK) {
    return c.fail("the explicit scheme is unstable: lower `dt` below z dx^2 / (4 dif max(1, eta)) = %g.",
                  par[9] * par[6] * par[6] / (4 * par[7] * std::fmax(1.0, par[2])));
  }
  for (double &e : burg) e = na(e);
  for (double &e : mb) e = na(e);
  out = vobj();
  put(out, "b", vmat(bv, nn, nn));
  put(out, "rho", vmat(rv, nn, nn));
  put(out, "tot", vmat(tot, nn, nn));
  put(out, "burg", vauto(burg));
  put(out, "meanb", vauto(mb));
  put(out, "mode", vstr(md[mo]));
  return true;
}

bool sdbthr(Ctx &c, Val &out) {
  double par[7] = {0};
  CK(pos1(c, "a0", NaN, par[0]));
  CK(pos1(c, "omega", NaN, par[1]));
  {
    Dv d;
    CK(num(c, c["eps"], "eps", d, 1, 0.0));
    par[2] = d[0];
  }
  CK(pos1(c, "dif", NaN, par[3]));
  CK(pos1(c, "gam", NaN, par[4]));
  {
    Dv d;
    CK(num(c, c["eta"], "eta", d, 1, 0.0));
    par[5] = d[0];
  }
  CK(pos1(c, "z", 1, par[6]));
  std::vector<double> kv;
  if (c.has("k")) {
    Dv d;
    CK(num(c, c["k"], "k", d, -1, 0.0));
    kv = d.v;
  }
  const int nk = static_cast<int>(kv.size());
  double st[12] = {0};
  std::vector<double> tr(kv.size()), de(kv.size()), gr(kv.size());
  clab_sdbthr(par, kv.data(), nk, st, tr.data(), de.data(), gr.data());
  static const char *const sn[] = {"b", "a", "rho"};
  static const char *const tn[] = {"lhs", "rhs"};
  Val disp;
  if (nk) {
    for (double &e : tr) e = na(e);
    for (double &e : de) e = na(e);
    for (double &e : gr) e = na(e);
    Frame f;
    f.num("k", kv);
    f.num("trace", tr);
    f.num("det", de);
    f.num("growth", gr);
    disp = f.o;
  }
  out = vobj();
  put(out, "steady", nvec(3, sn, st + 3));
  put(out, "p", vnum(st[0]));
  put(out, "q", vnum(st[1]));
  put(out, "r", vnum(st[2]));
  put(out, "thr", nvec(2, tn, st + 6));
  put(out, "unstable", vbool(st[8] != 0));
  put(out, "kcrit", vnum(na(st[9])));
  put(out, "detcrit", vnum(na(st[10])));
  put(out, "grcrit", vnum(na(st[11])));
  put(out, "disp", disp);
  return true;
}

bool etasf(Ctx &c, Val &out) {
  double h = 0.0;
  CK(pos1(c, "hor", NaN, h));
  std::vector<double> t;
  CK(evtimes(c, h, 10, t));
  tiewarn(c, t, "etasf");
  const int n = static_cast<int>(t.size());
  int mx = 0;
  CK(wholed(c, "maxit", 500, 1, mx));
  std::vector<double> dm(t.size(), 0.0);
  bool af = false;
  if (c.has("mag")) {
    Dv mv;
    CK(num(c, c["mag"], "mag", mv, n));
    double lo = mv[0], hi = lo;
    for (double e : mv.v) {
      lo = std::fmin(lo, e);
      hi = std::fmax(hi, e);
    }
    double r = lo;
    if (c.has("m0")) {
      Dv d;
      CK(num(c, c["m0"], "m0", d, 1));
      r = d[0];
    }
    for (std::size_t i = 0; i < t.size(); ++i) dm[i] = mv[i] - r;
    af = hi - lo > 1e-12;
  }
  std::vector<double> s0;
  if (c.has("start")) {
    Dv d;
    CK(num(c, c["start"], "start", d, 5));
    s0 = d.v;
  }
  double th[5] = {0}, cov[25] = {0}, st[4] = {0};
  if (clab_etasf(t.data(), dm.data(), n, h, af ? 1 : 0, s0.empty() ? nullptr : s0.data(), mx, th, cov, st, nullptr) !=
      CLAB_OK) {
    return c.fail("`start` must have mu, K and c above 0 and p above 1, finite alpha.");
  }
  const int k = af ? 5 : 4;
  out = vobj();
  pse(out, {"mu", "k", "alpha", "c", "p"}, th, cov);
  put(out, "loglik", vnum(st[0]));
  put(out, "aic", vnum(2 * k - 2 * st[0]));
  put(out, "branch", vnum(st[3]));
  put(out, "conv", vbool(st[2] > 0));
  put(out, "iter", vnum(static_cast<int>(st[1])));
  put(out, "afree", vbool(af));
  put(out, "n", vnum(n));
  put(out, "hor", vnum(h));
  return true;
}

bool mvhwk(Ctx &c, Val &out) {
  double h = 0.0;
  CK(pos1(c, "hor", NaN, h));
  std::vector<double> t;
  CK(evtimes(c, h, 10, t));
  tiewarn(c, t, "mvhwk");
  const int n = static_cast<int>(t.size());
  int mx = 0;
  CK(wholed(c, "maxit", 500, 1, mx));
  Lab l;
  CK(codes(c, c["type"], "type", t.size(), l));
  const int k = l.g;
  if (k > 10) return c.fail("`type` has %d types: at most 10.", k);
  double b = 0.0;
  if (c.has("beta")) CK(pos1(c, "beta", NaN, b));
  const int np = k + k * k + (b > 0 ? 0 : 1);
  std::vector<double> s0;
  if (c.has("start")) {
    Dv d;
    CK(num(c, c["start"], "start", d, np));
    s0 = d.v;
  }
  const std::size_t sp = static_cast<std::size_t>(np), sk = static_cast<std::size_t>(k);
  std::vector<double> th(sp), cov(sp * sp), st(4), w(2 * sk + 8 * sp + sp * sp);
  std::vector<int> iw(sk);
  if (clab_mvhwk(t.data(), l.code.data(), n, k, h, b, s0.empty() ? nullptr : s0.data(), mx, th.data(), cov.data(),
                 st.data(), nullptr, w.data(), iw.data()) != CLAB_OK) {
    return c.fail("`start` must hold positive numbers.");
  }
  std::vector<std::string> nm;
  for (int q = 0; q < k; ++q) nm.push_back("mu." + l.lev[static_cast<std::size_t>(q)]);
  for (int j = 0; j < k; ++j) {
    for (int q = 0; q < k; ++q) nm.push_back("a." + l.lev[static_cast<std::size_t>(q)] + "." + l.lev[static_cast<std::size_t>(j)]);
  }
  if (b <= 0) nm.emplace_back("beta");
  out = vobj();
  pse(out, nm, th.data(), cov.data());
  std::vector<double> mu(th.begin(), th.begin() + k);
  for (double &e : mu) e = na(e);
  put(out, "mu", nvec(l.lev, mu.data()));
  put(out, "a", vmat(th.data() + k, k, k));
  put(out, "beta", vnum(b > 0 ? b : th[sp - 1]));
  put(out, "rad", vnum(st[3]));
  put(out, "loglik", vnum(st[0]));
  put(out, "aic", vnum(2 * np - 2 * st[0]));
  put(out, "conv", vbool(st[2] > 0));
  put(out, "iter", vnum(static_cast<int>(st[1])));
  put(out, "n", vnum(n));
  put(out, "hor", vnum(h));
  return true;
}

bool sthwk(Ctx &c, Val &out) {
  double h = 0.0;
  CK(pos1(c, "hor", NaN, h));
  std::vector<double> t;
  CK(evtimes(c, h, 10, t));
  tiewarn(c, t, "sthwk");
  const int n = static_cast<int>(t.size());
  int mx = 0;
  CK(wholed(c, "maxit", 500, 1, mx));
  Dv xv, yv;
  CK(num(c, c["x"], "x", xv, n));
  CK(num(c, c["y"], "y", yv, n));
  double w[4] = {0};
  if (!c.has("win")) {
    w[0] = w[1] = xv[0];
    w[2] = w[3] = yv[0];
    for (std::size_t i = 0; i < xv.n(); ++i) {
      w[0] = std::fmin(w[0], xv[i]);
      w[1] = std::fmax(w[1], xv[i]);
      w[2] = std::fmin(w[2], yv[i]);
      w[3] = std::fmax(w[3], yv[i]);
    }
    if (!(w[0] < w[1] && w[2] < w[3])) return c.fail("the events lie on a line: give `win`.");
  } else {
    CK(win(c, c["win"], w));
    for (std::size_t i = 0; i < xv.n(); ++i) {
      if (!inside(w, xv[i], yv[i])) return c.fail("every event must lie in `win`.");
    }
  }
  const double area = (w[1] - w[0]) * (w[3] - w[2]);
  std::vector<double> s0;
  const bool hp = c.has("par");
  if (hp || c.has("start")) {
    Dv d;
    CK(num(c, c[hp ? "par" : "start"], hp ? "par" : "start", d, 4, 0.0));
    s0 = d.v;
  }
  std::vector<double> qt, qx, qy;
  const bool ha = c.has("at");
  if (ha) {
    const Val &at = c["at"];
    const char *msg = "`at` must be a data frame or list with columns t, x and y.";
    if (at.t != Val::OBJ) return c.fail("%s", msg);
    const char *want[3] = {"t", "x", "y"};
    std::vector<double> *dst[3] = {&qt, &qx, &qy};
    for (int q = 0; q < 3; ++q) {
      const Val &cl = col(at, want[q]);
      if (cl.t == Val::NUL) return c.fail("%s", msg);
      Dv d;
      CK(num(c, cl, "at", d));
      *dst[q] = d.v;
    }
    if (qx.size() != qt.size() || qy.size() != qt.size()) return c.fail("`at`'s columns must have one length.");
  }
  double th[4] = {0}, cov[16] = {0}, st[3] = {0};
  std::vector<double> ql(qt.size());
  if (clab_sthwk(t.data(), xv.p(), yv.p(), n, h, area, s0.empty() ? nullptr : s0.data(), hp ? 0 : mx, th, cov, st,
                 nullptr, qt.data(), qx.data(), qy.data(), static_cast<int>(qt.size()), ql.data()) != CLAB_OK) {
    return c.fail("`par` or `start`: mu, beta and sigma must be above 0 (and eta below 1 for a start).");
  }
  static const char *const wn[] = {"xmin", "xmax", "ymin", "ymax"};
  out = vobj();
  pse(out, {"mu", "eta", "beta", "sigma"}, th, cov);
  put(out, "loglik", vnum(st[0]));
  put(out, "aic", vnum(8 - 2 * st[0]));
  put(out, "conv", vbool(st[2] > 0));
  put(out, "iter", vnum(static_cast<int>(st[1])));
  put(out, "fit", vbool(!hp));
  put(out, "n", vnum(n));
  put(out, "hor", vnum(h));
  put(out, "win", nvec(4, wn, w));
  put(out, "area", vnum(area));
  for (double &e : ql) e = na(e);
  put(out, "lam", ha ? vauto(ql) : Val());
  return true;
}
