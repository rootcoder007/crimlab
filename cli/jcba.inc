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
