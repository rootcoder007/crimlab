// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Aoristic weights, digit tests (Benford) and constant reporters of a
// reporting panel.

// inst(): numbers of seconds (date-times in R); null is NA, none infinite.
bool inst(Ctx &c, const Val &x, const char *nm, std::vector<double> &v) {
  if (!isnum(x)) return c.fail("`%s` must be numbers of seconds or date-times (POSIXct).", nm);
  int dim[2] = {-1, -1};
  if (x.t == Val::ARR && shape(x, 2, dim, 0)) return c.fail("`%s` must be numbers of seconds or date-times (POSIXct).", nm);
  v = dbl(x);
  for (double e : v) {
    if (std::isinf(e)) return c.fail("`%s` must not be infinite.", nm);
  }
  return true;
}

bool aorwt(Ctx &c, Val &out) {
  std::vector<double> f, t, org, ttv;
  CK(inst(c, c["from"], "from", f));
  const int n = static_cast<int>(f.size());
  CK(inst(c, c["to"], "to", t));
  if (t.size() != f.size()) return c.fail("`to` must have the length of `from`.");
  double w = 0.0;
  CK(oned(c, "bw", 3600, 0.0, Inf, w));
  int k = 0;
  CK(wholed(c, "nb", 168, 1, k));
  Val d1;
  CK(inst(c, arg(c, "origin", d1, 259200), "origin", org));
  if (org.size() != 1 || std::isnan(org[0])) return c.fail("`origin` must be one instant.");
  Val d2;
  Dv o0;
  CK(num(c, arg(c, "off", d2, 0), "off", o0, 1));
  if (c.has("tt") != c.has("toff")) return c.fail("give both `tt` and `toff`, or neither.");
  Dv tov;
  if (c.has("tt")) {
    CK(inst(c, c["tt"], "tt", ttv));
    CK(num(c, c["toff"], "toff", tov, static_cast<long>(ttv.size())));
  }
  const int nt = static_cast<int>(ttv.size());
  const double l = val1d(c, "lo", 120), h = val1d(c, "hi", 604800);
  if (!(l >= 0) || !(h > l)) return c.fail("`lo` and `hi` must be numbers with 0 <= lo < hi.");
  bool ea = false;
  CK(flagd(c, "each", false, ea));
  std::vector<double> wm(ea ? static_cast<std::size_t>(n) * k : 0), tot(static_cast<std::size_t>(k));
  int cnt[5];
  if (clab_aorwt(f.data(), t.data(), n, ttv.data(), tov.v.data(), nt, o0[0], org[0], w, k, l, h,
                 ea ? wm.data() : nullptr, tot.data(), cnt) != CLAB_OK) {
    return c.fail("`tt` must be strictly increasing, with finite offsets.");
  }
  double cd[5];
  for (int j = 0; j < 5; ++j) cd[j] = cnt[j];
  static const char *const cn[] = {"left", "noend", "reversed", "short", "long"};
  out = vobj();
  put(out, "tot", vauto(tot));
  put(out, "w", ea ? vmat(wm, n, k) : Val());
  put(out, "cnt", nvec(5, cn, cd));
  return true;
}

bool bdig(Ctx &c, Val &out) {
  if (!isnum(c["x"]) || c["x"].t == Val::NUL) return c.fail("`x` must be numbers.");
  const std::vector<double> v = dbl(c["x"]);
  const int n = static_cast<int>(v.size());
  int d = 0;
  CK(wholed(c, "nd", 1, 1, d));
  if (d > 9) return c.fail("`nd` must be from 1 to 9.");
  const double dv = val1d(c, "dp", -1);
  if (!(dv == std::floor(dv) && dv >= -1 && dv <= 22)) return c.fail("`dp` must be a whole number from -1 to 22.");
  static const char *const so[] = {"pos", "neg", "both"};
  int s = 0;
  CK(pickd(c, "sgn", so, 3, 0, s));
  std::vector<int> dig(v.size());
  int cnt[4];
  if (clab_bdig(v.data(), n, d, static_cast<int>(dv), s == 0 ? 1 : (s == 1 ? -1 : 0), dig.data(), cnt) != CLAB_OK) {
    return c.fail("the kernel refused the input.");
  }
  double cd[4];
  for (int j = 0; j < 4; ++j) cd[j] = cnt[j];
  static const char *const cn[] = {"zero", "nonfinite", "sign", "range"};
  out = vobj();
  put(out, "dig", unbox(vints(dig.data(), dig.size())));
  put(out, "cnt", nvec(4, cn, cd));
  return true;
}

bool bnft(Ctx &c, Val &out) {
  Dv v;
  CK(num(c, c["dig"], "dig", v, -1, 0.0));
  int d = 0;
  CK(wholed(c, "nd", 1, 1, d));
  if (d > 6) return c.fail("`nd` must be from 1 to 6.");
  const int lo = static_cast<int>(std::pow(10.0, d - 1));
  std::vector<int> iv(v.n());
  for (std::size_t i = 0; i < v.n(); ++i) {
    const double e = v[i];
    if (e != std::floor(e) || (e != 0 && (e < lo || e >= 10.0 * lo))) {
      return c.fail("`dig` must hold %d-digit numbers (or 0, skipped), as bdig() gives.", d);
    }
    iv[i] = static_cast<int>(e);
  }
  int ns = 0;
  CK(wholed(c, "nsim", 0, 0, ns));
  unsigned int s = 0;
  if (ns > 0) CK(seed(c, s));
  const int k = 9 * lo;
  std::vector<double> obs(static_cast<std::size_t>(k)), wk(static_cast<std::size_t>(k));
  double o[6];
  if (clab_bnft(iv.data(), static_cast<int>(v.n()), d, ns, s, obs.data(), wk.data(), o) != CLAB_OK) {
    return c.fail("the kernel refused the input.");
  }
  std::vector<double> dg(static_cast<std::size_t>(k)), ex(static_cast<std::size_t>(k));
  for (int j = 0; j < k; ++j) {
    dg[static_cast<std::size_t>(j)] = lo + j;
    ex[static_cast<std::size_t>(j)] = o[0] * std::log10(1.0 + 1.0 / (lo + j));
  }
  Frame f;
  f.num("digit", dg);
  f.num("count", obs);
  f.num("exp", ex);
  static const char *const sn[] = {"m", "chisq", "df", "p", "pmc", "mad"};
  out = vobj();
  put(out, "obs", f.o);
  put(out, "stat", nvec(6, sn, o));
  return true;
}


bool rptc(Ctx &c, Val &out) {
  Args a(c, {"unit", "year", "rep"});
  CK(a.byname(c));
  Dv yv, rv;
  CK(num(c, a[1], "year", yv));
  const std::size_t n = yv.n();
  Lab u;
  CK(codes(c, a[0], "unit", n, u));
  CK(num(c, a[2], "rep", rv, static_cast<long>(n)));
  int a0 = 0, a1 = 0, mr = 0;
  CK(whole(c, c["y0"], "y0", 0, a0));
  CK(whole(c, c["y1"], "y1", a0, a1));
  CK(wholed(c, "minrep", 12, 0, mr));
  if (static_cast<double>(a1 - a0 + 1) * u.g > 1e9) return c.fail("too many units times years.");
  std::vector<int> yi(n), ri(n);
  for (std::size_t i = 0; i < n; ++i) {
    const double e = yv[i], r = rv[i];
    if (e != std::floor(e) || std::fabs(e) > 1e9 || r != std::floor(r) || std::fabs(r) > 1e9) {
      return c.fail("`year` and `rep` must be whole numbers.");
    }
    yi[i] = static_cast<int>(e);
    ri[i] = static_cast<int>(r);
  }
  std::vector<int> fl(static_cast<std::size_t>(u.g)), wk(static_cast<std::size_t>(u.g) * static_cast<std::size_t>(a1 - a0 + 1));
  int nc = 0;
  if (clab_rptc(u.code.data(), yi.data(), ri.data(), static_cast<int>(n), u.g, a0, a1, mr, fl.data(), &nc,
                wk.data()) != CLAB_OK) {
    return c.fail("the kernel refused the input.");
  }
  Frame f;
  f.str("unit", u.lev);
  f.lgl("const", fl.data(), fl.size());
  out = f.o;
  return true;
}
