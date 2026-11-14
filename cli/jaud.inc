// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Audits of risk scores (src/bdos.cpp): calibration across areas and the
// spread of scores between groups.

// A data frame given as JSON: an object whose columns are arrays.
inline bool isdf(const Val &o) { return o.t == Val::OBJ && !o.a.empty() && o.a[0].t == Val::ARR; }

bool calaud(Ctx &c, Val &out) {
  Val uv[4] = {c["area"], c["risk"], c["outcome"], c["group"]};
  Val pop = c["population"];
  if (isdf(uv[0])) {
    if (c.has("risk") || c.has("outcome") || c.has("group") || c.has("population")) {
      return c.fail("give `area` as a data frame of the columns, or every argument as a vector, not both.");
    }
    const Val df = uv[0];
    uv[0] = col(df, "area");
    uv[1] = col(df, "risk");
    uv[2] = col(df, "outcome");
    uv[3] = col(df, "group");
    pop = col(df, "population");
    if (uv[0].t == Val::NUL || uv[1].t == Val::NUL || uv[2].t == Val::NUL || uv[3].t == Val::NUL) {
      return c.fail("`area` as a data frame needs columns area, risk, outcome and group.");
    }
  }
  if (uv[1].t == Val::NUL || uv[2].t == Val::NUL || uv[3].t == Val::NUL) {
    return c.fail("`risk`, `outcome` and `group` are needed.");
  }
  CK(vbyname(c, {{"area", &uv[0]}, {"risk", &uv[1]}, {"outcome", &uv[2]}, {"group", &uv[3]}}));
  const std::size_t sn = vlen(uv[0]);
  const int n = static_cast<int>(sn);
  if (n < 1) return c.fail("`area` must hold at least one record.");
  Lab a, g;
  CK(codes(c, uv[0], "area", sn, a));
  Dv rv, ov;
  CK(num(c, uv[1], "risk", rv, n));
  CK(num(c, uv[2], "outcome", ov, n));
  CK(codes(c, uv[3], "group", sn, g));
  if (a.g < 3) return c.fail("at least 3 areas are needed.");
  std::vector<double> pp;
  if (pop.t != Val::NUL) {
    Dv pv;
    CK(num(c, pop, "population", pv, -1, 0.0, true));
    pp.resize(sn);
    if (!pv.nm.empty() && pv.n() == static_cast<std::size_t>(a.g)) {
      std::vector<int> at(static_cast<std::size_t>(a.g), -1);
      for (std::size_t k = 0; k < pv.n(); ++k) {
        for (std::size_t t = 0; t < at.size(); ++t) {
          if (pv.nm[k] == a.lev[t]) {
            if (at[t] >= 0) return c.fail("`population` names an area twice.");
            at[t] = static_cast<int>(k);
          }
        }
      }
      for (std::size_t t = 0; t < at.size(); ++t) {
        if (at[t] < 0) return c.fail("`population` has no value for area %s.", a.lev[t].c_str());
      }
      for (std::size_t i = 0; i < sn; ++i) pp[i] = pv[static_cast<std::size_t>(at[static_cast<std::size_t>(a.code[i])])];
    } else if (pv.n() == sn) {
      pp = pv.v;
    } else {
      return c.fail("`population` must be one number per record, or one per area named by area.");
    }
  }
  int B = 0;
  CK(wholed(c, "nperm", 999, 0, B));
  unsigned int s = 0;
  if (B > 0) CK(seed(c, s));
  const int na = a.g, G = g.g;
  const std::size_t sa = static_cast<std::size_t>(na), sg = static_cast<std::size_t>(G);
  std::vector<int> iw(sn + 2 * sa);
  std::vector<double> w(7 * sa + 3 * sg + 7 + sa + 3 * sg);
  double *ar = w.data(), *gr = ar + 7 * sa, *st = gr + 3 * sg, *wk = st + 7;
  const int rc = clab_calaud(a.code.data(), rv.p(), ov.p(), g.code.data(), pp.empty() ? nullptr : pp.data(), n, na, G,
                             B, s, ar, gr, st, iw.data(), wk);
  if (rc == CLAB_EDOM) return c.fail("`population` must take one value within each area.");
  if (rc != CLAB_OK) return c.fail("the kernel refused the input.");
  std::vector<std::string> ag(sa);
  for (std::size_t k = 0; k < sa; ++k) ag[k] = g.lev[static_cast<std::size_t>(ar[3 * sa + k])];
  Frame fa;
  fa.str("area", a.lev);
  fa.num("n", ar, sa);
  fa.num("risk", ar + sa, sa);
  fa.num("rate", ar + 2 * sa, sa);
  fa.str("group", ag);
  fa.num("riskrk", ar + 4 * sa, sa);
  fa.num("raterk", ar + 5 * sa, sa);
  fa.num("gap", ar + 6 * sa, sa);
  Frame fg;
  fg.str("group", g.lev);
  fg.num("n", gr, sg);
  fg.num("gap", gr + sg, sg);
  fg.num("p", gr + 2 * sg, sg);
  // the large-sample p-value of cor.test(): t on n - 2 degrees of freedom
  const double rho = st[1], df = na - 2.0;
  const double tt = rho * std::sqrt(df / ((1.0 - rho) * (1.0 + rho)));
  const double sp[5] = {rho, st[2], st[3], tt, std::isnan(tt) ? NaN : std::fmin(1.0, 2.0 * clab_pt(-std::fabs(tt), df, 1))};
  static const char *const spn[] = {"rho", "s", "pex", "t", "pt"};
  static const char *const bn[] = {"ss", "p"};
  out = vobj();
  put(out, "areas", fa.o);
  put(out, "groups", fg.o);
  put(out, "spear", nvec(5, spn, sp));
  put(out, "ties", vbool(st[4] > 0.0));
  put(out, "between", nvec(2, bn, st + 5));
  put(out, "nperm", vnum(B));
  return true;
}

bool scrdsp(Ctx &c, Val &out) {
  Args a(c, {"score", "group"});
  CK(a.byname(c));
  Dv xv;
  CK(num(c, a[0], "score", xv));
  const std::size_t sn = xv.n();
  const int n = static_cast<int>(sn);
  Lab g;
  CK(codes(c, a[1], "group", sn, g));
  if (g.g < 2) return c.fail("`group` must hold at least two groups.");
  int r = -1;
  if (c.has("ref")) {
    CK(ref(c, c["ref"], g.lev, r));
    --r;
  }
  double cf = 0.0;
  CK(oned(c, "conf", 0.95, 0.0, 1.0, cf));
  const int G = g.g;
  const std::size_t sg = static_cast<std::size_t>(G);
  std::vector<int> iw(sn);
  std::vector<double> w(6 * sg + 9);
  double *gs = w.data(), *st = gs + 6 * sg;
  if (clab_scrdsp(xv.p(), g.code.data(), n, G, r, gs, st, iw.data()) != CLAB_OK) {
    return c.fail("the kernel refused the input.");
  }
  const double q = st[2] > 0.0 ? clab_qt(0.5 + 0.5 * cf, st[2], 1) : NaN;
  Frame f;
  f.str("group", g.lev);
  f.num("n", gs, sg);
  static const char *const cn[] = {"mean", "median", "sd", "gap", "se"};
  for (int t = 0; t < 5; ++t) f.num(cn[t], gs + static_cast<std::size_t>(t + 1) * sg, sg);
  std::vector<double> lo(sg), hi(sg);
  for (std::size_t j = 0; j < sg; ++j) {
    lo[j] = gs[4 * sg + j] - q * gs[5 * sg + j];
    hi[j] = gs[4 * sg + j] + q * gs[5 * sg + j];
  }
  f.num("lo", lo);
  f.num("hi", hi);
  auto pf = [](double F, double d1, double d2) { return std::isnan(F) ? NaN : clab_pf(F, d1, d2, 0); };
  const double an[4] = {st[0], st[1], st[2], pf(st[0], st[1], st[2])};
  const double wl[4] = {st[3], st[4], st[5], pf(st[3], st[4], st[5])};
  static const char *const fn[] = {"f", "df1", "df2", "p"};
  out = vobj();
  put(out, "tab", f.o);
  put(out, "anova", nvec(4, fn, an));
  put(out, "welch", nvec(4, fn, wl));
  put(out, "spread", vnum(st[6]));
  put(out, "ref", vstr(g.lev[static_cast<std::size_t>(st[7])]));
  return true;
}
