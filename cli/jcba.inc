// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Covariate balance (src/bcba.cpp): standardised differences and overlap of
// two arms, and the E-value of an estimate.

// xmat(): covariates, an object of columns, a matrix of n rows or one
// vector; column-major into v, their names into nm.
bool cxmat(Ctx &c, const Val &x, std::size_t n, std::vector<double> &v, std::vector<std::string> &nm) {
  v.clear();
  nm.clear();
  int dim[2] = {-1, -1};
  if (x.t == Val::ARR && !x.a.empty() && x.a[0].t == Val::ARR) {
    if (!shape(x, 2, dim, 0) || !isnum(x.a[0]) || static_cast<std::size_t>(dim[0]) != n) {
      return c.fail("`x` must be a numeric matrix with one row per unit, %d.", static_cast<int>(n));
    }
    v.assign(n * static_cast<std::size_t>(dim[1]), 0.0);
    bool bad = false;
    fill(x, 2, dim, 0, 0, 1, v, bad, true);
    for (double e : v) {
      if (!std::isfinite(e)) return c.fail("`x` must hold finite numbers with no missing values.");
    }
    for (int j = 0; j < dim[1]; ++j) nm.push_back("x" + std::to_string(j + 1));
    if (nm.empty()) return c.fail("`x` must hold one covariate at least.");
    return true;
  }
  std::vector<const Val *> cs;
  std::vector<std::string> names;
  if (x.t == Val::OBJ) {
    for (const Val &e : x.a) cs.push_back(&e);
    names = x.k;
  } else {
    cs.push_back(&x);
  }
  if (cs.empty()) return c.fail("`x` must hold one covariate at least.");
  for (std::size_t j = 0; j < cs.size(); ++j) {
    Dv d;
    CK(num(c, *cs[j], "x", d, static_cast<long>(n)));
    v.insert(v.end(), d.v.begin(), d.v.end());
    nm.push_back(j < names.size() && !names[j].empty() ? names[j] : "x" + std::to_string(j + 1));
  }
  return true;
}

bool smdbal(Ctx &c, Val &out) {
  double th = 0.0;
  CK(oned(c, "thr", 0.1, 0.0, HUGE_VAL, th));
  const std::size_t n = vlen(c["d"]);
  std::vector<int> d01;
  CK(bin01(c, c["d"], "d", n, d01));
  std::vector<double> dv(d01.begin(), d01.end());
  const double n1 = static_cast<double>(std::count(dv.begin(), dv.end(), 1.0));
  if (n1 < 2 || static_cast<double>(n) - n1 < 2) return c.fail("`d` must have two treated and two untreated units at least.");
  std::vector<double> xv, wv;
  std::vector<std::string> nm;
  CK(cxmat(c, c["x"], n, xv, nm));
  if (c.has("w")) {
    Dv w;
    CK(num(c, c["w"], "w", w, static_cast<long>(n), 0.0));
    wv = w.v;
  }
  const int N = static_cast<int>(n), p = static_cast<int>(nm.size());
  pev::cdat cd{dv.data(), dv.data(), xv.data(), N, p};
  std::vector<double> e;
  if (!c.has("ps")) {
    std::vector<int> all(n);
    for (std::size_t i = 0; i < n; ++i) all[i] = static_cast<int>(i);
    if (!pev::pscore(cd, all, e)) {
      return c.fail("the logistic propensity has no finite fit (a covariate separates the arms); give `ps`.");
    }
  } else {
    Dv pp;
    CK(num(c, c["ps"], "ps", pp, static_cast<long>(n)));
    for (double x : pp.v) {
      if (!(x > 0 && x < 1)) return c.fail("`ps` must be propensities strictly between 0 and 1.");
    }
    e = pp.v;
  }
  std::vector<double> m1, m0, sm, vr, a1, a0, as, av;
  bool bal = true;
  for (int j = 0; j < p; ++j) {
    double o[4];
    const double *xj = xv.data() + static_cast<std::size_t>(j) * n;
    pev::smd(xj, N, dv.data(), nullptr, o);
    m1.push_back(o[0]);
    m0.push_back(o[1]);
    sm.push_back(o[2]);
    vr.push_back(o[3]);
    if (!wv.empty()) {
      pev::smd(xj, N, dv.data(), wv.data(), o);
      a1.push_back(o[0]);
      a0.push_back(o[1]);
      as.push_back(o[2]);
      av.push_back(o[3]);
    }
    bal = bal && std::fabs(wv.empty() ? sm.back() : as.back()) < th;
  }
  Frame t;
  t.str("var", nm);
  t.num("m1", m1);
  t.num("m0", m0);
  t.num("smd", sm);
  t.num("vr", vr);
  if (!wv.empty()) {
    t.num("m1w", a1);
    t.num("m0w", a0);
    t.num("smdw", as);
    t.num("vrw", av);
  }
  double lo[2] = {1, 1}, hi[2] = {0, 0}, wmax = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const int g = dv[i] > 0.5 ? 1 : 0;
    lo[g] = std::min(lo[g], e[i]);
    hi[g] = std::max(hi[g], e[i]);
    wmax = std::max(wmax, g ? 1 / e[i] : 1 / (1 - e[i]));
  }
  const double cl = std::max(lo[0], lo[1]), ch = std::min(hi[0], hi[1]);
  double ou = 0;
  for (double v : e) ou += (v < cl || v > ch) ? 1.0 / static_cast<double>(n) : 0.0;
  Frame ov;
  const double ovv[6] = {lo[1], hi[1], lo[0], hi[0], ou, wmax};
  static const char *const on[] = {"min1", "max1", "min0", "max0", "out", "wmax"};
  for (int k = 0; k < 6; ++k) ov.num(on[k], ovv + k, 1);
  out = vobj();
  put(out, "tab", t.o);
  put(out, "overlap", ov.o);
  put(out, "bal", vbool(bal));
  put(out, "ovl", vbool(ou == 0));
  put(out, "thr", vnum(th));
  return true;
}

bool evlx(Ctx &c, Val &out) {
  static const char *const to[] = {"rr", "or", "hr", "smd"};
  int ty = 0;
  CK(pickd(c, "type", to, 4, 0, ty));
  bool rr = false;
  CK(flagd(c, "rare", false, rr));
  const double lw = ty == 3 ? -HUGE_VAL : 0.0;
  double v3[3];
  static const char *const an[] = {"est", "lo", "hi"};
  for (int k = 0; k < 3; ++k) {
    v3[k] = NaN;
    if (!c.has(an[k])) continue;
    Dv q;
    CK(num(c, c[an[k]], an[k], q, 1, lw, ty != 3));
    v3[k] = q[0];
  }
  const double e = v3[0], l = v3[1], h = v3[2];
  if (std::isnan(e)) return c.fail("`est` must be one estimate.");
  if (!std::isnan(l) && !std::isnan(h) && !(l <= e && e <= h)) return c.fail("`lo` and `hi` must surround `est`.");
  double s = NaN;
  if (c.has("sd")) {
    if (ty != 3) return c.fail("`sd` is the standard error of a standardised difference (`type = \"smd\"`).");
    CK(one(c, c["sd"], "sd", 0.0, HUGE_VAL, s));
  }
  double o[3];
  pev::evalue(e, l, h, ty, rr, s, o);
  for (double &v : o) v = na(v);
  static const char *const on[] = {"rr", "evalue", "evlim"};
  out = nvec(3, on, o);
  return true;
}

bool cbatt(Ctx &c, Val &out) {
  static const char *const eo[] = {"ipw", "gcomp", "sub", "psm", "ebal"};
  static const char *const to[] = {"ate", "att", "atc"};
  int me = 0, tg = 0, kk = 0, ns = 0, nr = 0;
  CK(pickd(c, "est", eo, 5, 0, me));
  CK(pickd(c, "tgt", to, 3, 0, tg));
  CK(wholed(c, "k", 1, 1, kk));
  CK(wholed(c, "nsub", 5, 2, ns));
  CK(wholed(c, "reps", 499, 0, nr));
  double lv = 0.0;
  CK(oned(c, "level", 0.95, 0.0, 1.0, lv));
  const Val &y = c["y"];
  if (y.t == Val::ARR) {
    for (const Val &e : y.a) {
      if (e.t == Val::ARR || e.t == Val::OBJ) return notvec(c, "y");
    }
  }
  Dv yv;
  CK(num(c, y, "y", yv));
  const std::size_t sn = yv.n();
  const int n = static_cast<int>(sn);
  std::vector<int> di;
  CK(bin01(c, c["d"], "d", sn, di));
  std::vector<double> dv(di.begin(), di.end());
  const double cn1 = static_cast<double>(std::count(dv.begin(), dv.end(), 1.0));
  if (cn1 < 2 || n - cn1 < 2) return c.fail("`d` must have two treated and two untreated units at least.");
  std::vector<double> xv;
  std::vector<std::string> nm;
  CK(cxmat(c, c["x"], sn, xv, nm));
  if (me == 4 && tg != 1) return c.fail("entropy balancing (`est = \"ebal\"`) estimates the ATT only.");
  if (me == 3 && tg == 0) return c.fail("matching (`est = \"psm\"`) estimates the ATT or the ATC.");
  pev::cdat cd{yv.p(), dv.data(), xv.data(), n, static_cast<int>(nm.size())};
  std::vector<int> all(sn);
  for (int i = 0; i < n; ++i) all[static_cast<std::size_t>(i)] = i;
  const bool need = me == 0 || me == 2 || me == 3, own = need && !c.has("ps");
  std::vector<double> e(sn, NaN);
  if (own) {
    if (!pev::pscore(cd, all, e)) {
      return c.fail("the logistic propensity has no finite fit (a covariate separates the arms); give `ps`.");
    }
  } else if (need) {
    Dv p;
    CK(num(c, c["ps"], "ps", p, n));
    e = p.v;
    for (double v : e) {
      if (!(v > 0 && v < 1)) return c.fail("`ps` must be propensities strictly between 0 and 1.");
    }
  }
  std::vector<double> w;
  auto one1 = [&](const std::vector<int> &idx, std::vector<double> &ee) {
    if (own && &ee != &e && !pev::pscore(cd, idx, ee)) return NaN;
    switch (me) {
      case 0: return pev::ipw(cd, idx, ee, tg);
      case 1: return pev::gcomp(cd, idx, tg);
      case 2: return pev::subcl(cd, idx, ee, tg, ns);
      default: {
        std::vector<double> ww;
        return pev::ebal(cd, idx, ww) ? pev::ebatt(cd, idx, ww) : NaN;
      }
    }
  };
  double b = NA, se = NA, med = NA;
  Val mt;
  if (me == 3) {
    double cal = 0.0;
    if (c.has("calip")) {
      double m = 0.0, s = 0.0, k1 = 0.0;
      for (double v : e) m += v / n;
      for (double v : e) s += (v - m) * (v - m) / (n - 1);
      CK(one(c, c["calip"], "calip", 0.0, Inf, k1));
      cal = k1 * std::sqrt(s);
    }
    std::vector<int> ft, fc;
    std::vector<double> dif;
    const int m = pev::psm(cd, e, tg, kk, cal, ft, fc, dif);
    if (m < 2) return c.fail("fewer than two units were matched within the caliper.");
    double s = 0.0, ss = 0.0;
    for (double v : dif) s += v / m;
    for (double v : dif) ss += (v - s) * (v - s) / (m - 1);
    b = s;
    se = std::sqrt(ss / m);
    std::vector<double> sd = dif;
    std::sort(sd.begin(), sd.end());
    med = m % 2 ? sd[static_cast<std::size_t>(m / 2)] : 0.5 * (sd[static_cast<std::size_t>(m / 2 - 1)] + sd[static_cast<std::size_t>(m / 2)]);
    std::vector<double> a, bb;
    for (std::size_t q = 0; q < ft.size(); ++q) {
      a.push_back(ft[q] + 1);
      bb.push_back(fc[q] + 1);
    }
    Frame cm;
    cm.num("focal", a);
    cm.num("match", bb);
    mt = cm.o;
  } else {
    b = one1(all, e);
    if (std::isnan(b)) {
      return c.fail("the estimate is undefined: an outcome model failed, a subclass lacks an arm, or the balancing "
                    "weights do not exist.");
    }
    if (me == 0 || me == 4) {
      w.assign(sn, 1.0);
      if (me == 4) {
        pev::ebal(cd, all, w);
      } else {
        for (std::size_t i = 0; i < sn; ++i) {
          const double p = e[i];
          w[i] = dv[i] > 0.5 ? (tg == 0 ? 1 / p : tg == 1 ? 1 : (1 - p) / p) : (tg == 0 ? 1 / (1 - p) : tg == 1 ? p / (1 - p) : 1);
        }
      }
    }
    if (nr > 0) {
      unsigned int sd = 0;
      CK(seed(c, sd));
      kutil::rng g(sd);
      kutil::ksum s1, s2;
      double kq = 0;
      std::vector<int> idx(sn);
      std::vector<double> eb = e;
      for (int r = 0; r < nr; ++r) {
        for (std::size_t i = 0; i < sn; ++i) idx[i] = g.below(n);
        const double v = one1(idx, eb);
        if (std::isnan(v)) continue;
        s1.add(v);
        s2.add(v * v);
        kq += 1;
      }
      if (kq > 1) {
        const double m = s1.val() / kq;
        se = std::sqrt(std::max(0.0, (s2.val() - kq * m * m) / (kq - 1)));
      }
    }
  }
  const double zq = clab_qnorm(0.5 + 0.5 * lv);
  const double lo = std::isnan(se) ? NA : b - zq * se, hi = std::isnan(se) ? NA : b + zq * se, n0 = n - cn1;
  Frame r;
  r.num("est", &b, 1);
  r.num("se", &se, 1);
  r.num("lo", &lo, 1);
  r.num("hi", &hi, 1);
  r.num("med", &med, 1);
  r.num("n1", &cn1, 1);
  r.num("n0", &n0, 1);
  out = vobj();
  put(out, "res", r.o);
  put(out, "ps", vauto(e));
  put(out, "w", w.empty() ? Val() : vauto(w));
  put(out, "match", mt);
  put(out, "est", vstr(eo[me]));
  put(out, "tgt", vstr(to[tg]));
  put(out, "level", vnum(lv));
  return true;
}
