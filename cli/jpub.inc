// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Published figures (src/bpub.cpp): changes from one period to the next
// with their exact intervals and the envelope of rounded counts.

// opt(): left out (NaN), or one number from lo up.
bool optd(Ctx &c, const char *nm, double lo, double &v) {
  v = NaN;
  if (!c.has(nm)) return true;
  return binr(c, c[nm], nm, lo, Inf, v);
}

// The exact two-sided binomial p-value of x of n at p, as binom.test().
double binp(double x, double n, double p) {
  const double d = clab_dbinom(x, n, p), m = n * p, re = 1 + 1e-7;
  if (x == m) return 1.0;
  if (x < m) {
    double a = std::ceil(m), b = n + 1;
    while (a < b) {
      const double h = std::floor((a + b) / 2);
      if (clab_dbinom(h, n, p) <= d * re) b = h;
      else a = h + 1;
    }
    return std::fmin(1.0, clab_pbinom(x, n, p, 1) + (a <= n ? clab_pbinom(a - 1, n, p, 0) : 0));
  }
  double a = -1, b = std::floor(m);
  while (a < b) {
    const double h = std::ceil((a + b) / 2);
    if (clab_dbinom(h, n, p) <= d * re) a = h;
    else b = h - 1;
  }
  return std::fmin(1.0, (a >= 0 ? clab_pbinom(a, n, p, 1) : 0) + clab_pbinom(x - 1, n, p, 0));
}

const char *const rkinds[] = {"nearest", "random"};

bool yoypx(Ctx &c, Val &out) {
  static const char *const un[] = {"count", "continuous", "percent"};
  int u = 0;
  CK(pickd(c, "units", un, 3, 0, u));
  Dv xv, pv;
  CK(num(c, c["x"], "x", xv, -1, u == 0 ? 0.0 : -HUGE_VAL));
  const std::size_t sn = xv.n();
  const int n = static_cast<int>(sn);
  CK(num(c, c["prd"], "prd", pv, n));
  double lg = 0.0, mb = 0.0, lv = 0.0, rn = 0.0, al = 0.0;
  CK(oned(c, "lag", 1, 0.0, Inf, lg));
  CK(binrd(c, "minbase", 0, 0.0, Inf, mb));
  CK(oned(c, "level", 0.95, 0.0, 1.0, lv));
  CK(optd(c, "rnd", 0.0, rn));
  if (rn == 0) return c.fail("`rnd` must be above 0.");
  if (u != 0 && !std::isnan(rn)) return c.fail("`rnd` goes with counts (`units = \"count\"`).");
  int kd = 0, ad = 0;
  CK(pickd(c, "kind", rkinds, 2, 0, kd));
  CK(pickd(c, "adj", padjs, 6, 4, ad));
  CK(oned(c, "alpha", 0.05, 0.0, 1.0, al));
  Lab g;
  const bool hb = c.has("by");
  if (hb) CK(codes(c, c["by"], "by", sn, g));
  std::vector<int> o(sn);
  for (int i = 0; i < n; ++i) o[static_cast<std::size_t>(i)] = i;
  const int *gi = hb ? g.code.data() : nullptr;
  const double *pp = pv.p();
  std::stable_sort(o.begin(), o.end(), [&](int a, int b) {
    const int ga = gi ? gi[a] : 0, gb = gi ? gi[b] : 0;
    return ga < gb || (ga == gb && pp[a] < pp[b]);
  });
  std::vector<double> cx, cp;
  std::vector<int> cgr, cn;
  for (int i = 0; i < n; ++i) {
    const int a = o[static_cast<std::size_t>(i)];
    const int b = i ? o[static_cast<std::size_t>(i - 1)] : -1;
    if (i == 0 || (gi ? gi[a] != gi[b] : false) || pp[a] != pp[b]) {
      cx.push_back(0);
      cp.push_back(pp[a]);
      cgr.push_back(gi ? gi[a] : 0);
      cn.push_back(0);
    }
    cx.back() += xv[static_cast<std::size_t>(a)];
    ++cn.back();
  }
  const std::size_t sm = cx.size();
  const int m = static_cast<int>(sm);
  if (u != 0) {
    for (std::size_t j = 0; j < sm; ++j) cx[j] /= cn[j];
  }
  std::vector<double> prev(sm), chg(sm), pct(sm);
  const std::vector<double> val = cx;
  clab_yoypx(cx.data(), cp.data(), cgr.data(), m, lg, prev.data(), chg.data(), pct.data());
  if (u == 2) pct = chg;
  std::vector<int> fc(sm, -1);
  static const char *const fls[] = {"no earlier period", "earlier period is 0", "small base"};
  for (std::size_t j = 0; j < sm; ++j) {
    if (std::isnan(prev[j])) fc[j] = 0;
    else if (u != 2 && prev[j] == 0) fc[j] = 1;
    else if (u != 2 && prev[j] < mb) fc[j] = 2;
  }
  for (std::size_t j = 0; j < sm; ++j) {
    prev[j] = na(prev[j]);
    chg[j] = na(chg[j]);
    pct[j] = na(pct[j]);
  }
  Frame f;
  if (hb) {
    std::vector<std::string> gl(sm);
    for (std::size_t j = 0; j < sm; ++j) gl[j] = g.lev[static_cast<std::size_t>(cgr[j])];
    f.str("group", gl);
  }
  f.num("prd", cp);
  f.num("value", val);
  f.num("prev", prev);
  f.num("change", chg);
  f.num(u == 2 ? "ppchg" : "pctchg", pct);
  Val fv = varr();
  for (int k : fc) fv.a.push_back(k < 0 ? Val() : vstr(fls[k]));
  f.add("flag", fv);
  if (u == 0) {
    const double a = 1 - lv;
    std::vector<double> clo(sm, NA), chi(sm, NA), pv2(sm, NA), pa(sm, NA), pk, qk(sm);
    std::vector<int> sg(sm, -1), oi(sm + 1);
    for (std::size_t j = 0; j < sm; ++j) {
      const double v = val[j], b = prev[j];
      if (std::isnan(b) || v + b == 0) continue;
      const double l = v == 0 ? 0 : clab_qbeta(a / 2, v, b + 1, 1);
      const double h = b == 0 ? 1 : clab_qbeta(1 - a / 2, v + 1, b, 1);
      clo[j] = 100 * (l / (1 - l) - 1);
      chi[j] = h >= 1 ? Inf : 100 * (h / (1 - h) - 1);
      pv2[j] = binp(std::round(v), std::round(v + b), 0.5);
      pk.push_back(pv2[j]);
    }
    clab_pbadj(pk.data(), static_cast<int>(pk.size()), ad, qk.data(), oi.data());
    for (std::size_t j = 0, i = 0; j < sm; ++j) {
      if (std::isnan(pv2[j])) continue;
      pa[j] = qk[i++];
      sg[j] = pa[j] < al;
    }
    f.num("pctlo", clo);
    f.num("pctup", chi);
    f.num("p", pv2);
    f.num("padj", pa);
    f.lgl("signif", sg.data(), sm);
    if (!std::isnan(rn)) {
      std::vector<double> bw(12 * sm);
      std::vector<int> bs(3 * sm, 0);
      double *nlo = bw.data(), *nhi = nlo + m, *plo = nhi + m, *phi = plo + m, *eo = phi + m;
      int *zst = bs.data(), *s1 = zst + m, *s2 = s1 + m;
      std::vector<double> lim(sm, NaN), pz(sm);
      for (std::size_t j = 0; j < sm; ++j) pz[j] = std::isnan(prev[j]) ? NaN : prev[j];
      std::vector<double> vv = val;
      clab_pbbnd(vv.data(), zst, lim.data(), m, rn, kd, NaN, nlo, nhi, s1);
      clab_pbbnd(pz.data(), zst, lim.data(), m, rn, kd, NaN, plo, phi, s2);
      clab_pbchg(nlo, nhi, plo, phi, m, nullptr, nullptr, 1, eo);
      std::vector<double> el(sm), eh(sm), ql(sm), qh(sm), cl(sm), ch(sm);
      for (std::size_t j = 0; j < sm; ++j) {
        el[j] = na(eo[j]);
        eh[j] = na(eo[sm + j]);
        ql[j] = na(eo[2 * sm + j]);
        qh[j] = na(eo[3 * sm + j]);
        cl[j] = std::isnan(clo[j]) || std::isnan(ql[j]) ? NA : std::fmin(clo[j], ql[j]);
        ch[j] = std::isnan(chi[j]) || std::isnan(qh[j]) ? NA : std::fmax(chi[j], qh[j]);
      }
      f.num("envlo", el);
      f.num("envup", eh);
      f.num("envplo", ql);
      f.num("envpup", qh);
      f.num("cpctlo", cl);
      f.num("cpctup", ch);
    }
  }
  out = f.o;
  return true;
}
