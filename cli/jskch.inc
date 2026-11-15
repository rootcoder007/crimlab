// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// One-pass summaries (src/skch.cpp), with the checks of R/skch.R in their
// words. An accumulator is the object onlsmry() returns: {"state": [5]}.

// .sknum(): numbers, text that reads as numbers, or missing values.
inline bool sknum(Ctx &c, const Val &o, const char *w, std::vector<double> &v) {
  v.clear();
  if (o.t == Val::NUL) return true;
  const std::vector<Val> one = {o};
  const std::vector<Val> &es = o.t == Val::ARR ? o.a : one;
  bool anyb = false;
  std::string firstb;
  for (const Val &e : es) {
    double x = 0;
    if (e.t == Val::NUM) v.push_back(e.n);
    else if (e.t == Val::NUL) v.push_back(NA);
    else if (snum(e, x)) v.push_back(x);
    else if (e.t == Val::STR) {
      char *end = nullptr;
      std::string s = e.s;
      const double d = std::strtod(s.c_str(), &end);
      bool blank = true;
      for (char ch : s) blank = blank && std::isspace(static_cast<unsigned char>(ch));
      if (s == "NA" || blank) v.push_back(NA);
      else if (end != nullptr && *end == '\0' && end != s.c_str()) v.push_back(d);
      else return ckstop(c, "`%s` must be numeric, not text such as \"%s\"", w, s.c_str());
    } else if (e.t == Val::BOOL) {
      if (!anyb) firstb = e.b ? "TRUE" : "FALSE";
      anyb = true;
      v.push_back(e.b);
    } else {
      return ckstop(c, "`%s` must be numeric, got list", w);
    }
  }
  if (anyb) return ckstop(c, "`%s` must be a number, not %s", w, firstb.c_str());
  return true;
}

inline bool skacc(Ctx &c, const char *nm, std::vector<double> &st, const char *msg) {
  const Val &o = c[nm];
  Dv d;
  if (o.t != Val::OBJ || readd(col(o, "state"), d) != NUM || d.n() != 5) return ckstop(c, "%s", msg);
  st = d.v;
  return true;
}

inline Val skout(const std::vector<double> &st) {
  Val o = vobj();
  put(o, "state", vnums(st));
  return o;
}

inline bool skupd(Ctx &c, std::vector<double> &st, const char *nm) {
  std::vector<double> x;
  CK(sknum(c, c[nm], "x", x));
  std::vector<double> y;
  for (double e : x) {
    if (std::isnan(e)) continue;
    if (std::isinf(e)) return ckstop(c, "`x` has infinite values, which no moment can hold; set them to NA first");
    y.push_back(e);
  }
  if (y.empty()) return true;
  double b[5], o[5];
  clab_momacc(y.data(), static_cast<double>(y.size()), b);
  clab_mommrg(st.data(), b, o);
  st.assign(o, o + 5);
  return true;
}

bool onlsmry(Ctx &c, Val &out) {
  std::vector<double> st(5, 0.0);
  if (c.has("x")) CK(skupd(c, st, "x"));
  out = skout(st);
  return true;
}

bool smryupd(Ctx &c, Val &out) {
  std::vector<double> st;
  CK(skacc(c, "acc", st, "`acc` must come from onlsmry()"));
  CK(skupd(c, st, "x"));
  out = skout(st);
  return true;
}

bool smrymrg(Ctx &c, Val &out) {
  std::vector<double> a, b;
  CK(skacc(c, "a", a, "both arguments must come from onlsmry()"));
  CK(skacc(c, "b", b, "both arguments must come from onlsmry()"));
  double o[5];
  clab_mommrg(a.data(), b.data(), o);
  out = skout(std::vector<double>(o, o + 5));
  return true;
}

bool smrystat(Ctx &c, Val &out) {
  std::vector<double> st;
  CK(skacc(c, "acc", st, "`acc` must come from onlsmry()"));
  double o[6];
  clab_momstat(st.data(), o);
  static const char *const nm[] = {"n", "mean", "variance", "sd", "skewness", "kurtosis"};
  out = nvec(6, nm, o);
  return true;
}

inline bool rsvarg(Ctx &c, const char *nm, double &v) {
  std::vector<double> x;
  CK(sknum(c, c[nm], nm, x));
  if (x.size() != 1 || std::isnan(x[0]) || !std::isfinite(x[0]) || x[0] < 0 || x[0] != std::floor(x[0])) {
    return ckstop(c, "`%s` must be a single non-negative whole number", nm);
  }
  if (x[0] > INT_MAX) return ckstop(c, "`%s` must be at most 2147483647", nm);
  v = x[0];
  return true;
}

inline bool rsvrun(Ctx &c, double n, std::vector<double> &idx) {
  double k = 0, seed = 42;
  CK(rsvarg(c, "k", k));
  if (c.has("seed")) {
    std::vector<double> s;
    CK(sknum(c, c["seed"], "seed", s));
    if (s.size() != 1 || std::isnan(s[0]) || s[0] < 0 || s[0] != std::floor(s[0]) || s[0] >= 9007199254740992.0) {
      return ckstop(c, "`seed` must be a whole number in [0, 2^53)");
    }
    seed = s[0];
  }
  idx.assign(static_cast<std::size_t>(std::min(n, k)), 0.0);
  return kchk(c, clab_rsv(n, k, seed, idx.data()));
}

bool rsvidx(Ctx &c, Val &out) {
  double n = 0;
  CK(rsvarg(c, "n", n));
  std::vector<double> idx;
  CK(rsvrun(c, n, idx));
  out = idx.size() == 1 ? vnum(idx[0]) : vnums(idx);
  return true;
}

bool rsvsmp(Ctx &c, Val &out) {
  const Val &x = c["x"];
  const std::size_t n = x.t == Val::ARR ? x.a.size() : x.t == Val::NUL ? 0 : 1;
  std::vector<double> idx;
  Ctx r{"rsvidx", c.in, std::string(), {}};
  if (!rsvrun(r, static_cast<double>(n), idx)) {
    c.err = r.err;
    return false;
  }
  Val o = varr();
  for (double i : idx) o.a.push_back(x.t == Val::ARR ? x.a[static_cast<std::size_t>(i) - 1] : x);
  out = o.a.size() == 1 ? o.a[0] : o;
  return true;
}

inline bool skreg(Ctx &c, const char *nm, std::vector<int> &r) {
  Dv d;
  const Val &o = c[nm];
  bool lg = false;
  if (jcls(o, lg) != nullptr && !lg) return ckstop(c, "`%s` must be numeric, got %s", nm, jcls(o, lg));
  readd(o, d);
  r.clear();
  for (double v : d.v) {
    if (!std::isfinite(v) || v < 0 || v > 64 || v != std::floor(v)) {
      return ckstop(c, "`%s` must be non-negative integer registers", nm);
    }
    r.push_back(static_cast<int>(v));
  }
  return true;
}

bool dstnskch(Ctx &c, Val &out) {
  double p = 14;
  if (c.has("p")) {
    std::vector<double> v;
    CK(sknum(c, c["p"], "p", v));
    if (v.size() != 1 || std::isnan(v[0]) || v[0] < 4 || v[0] > 20 || v[0] != std::nearbyint(v[0])) {
      return ckstop(c, "`p` must be a single integer between 4 and 20");
    }
    p = v[0];
  }
  const int pp = static_cast<int>(p);
  std::vector<int> reg(static_cast<std::size_t>(1) << pp, 0);
  if (c.has("regstr")) {
    std::vector<int> r;
    CK(skreg(c, "regstr", r));
    if (r.size() != reg.size()) return ckstop(c, "`regstr` must have 2^p entries, from a sketch with the same `p`");
    reg = r;
  }
  const Val &x = c["x"];
  const std::vector<Val> one = {x};
  for (const Val &e : x.t == Val::ARR ? x.a : one) {
    std::string s;
    if (e.t == Val::NUL || !lab(e, s)) continue;
    clab_hlladd(s.c_str(), s.size(), pp, reg.data());
  }
  out = vints(reg);
  return true;
}

bool dstncnt(Ctx &c, Val &out) {
  std::vector<int> r;
  if (c.has("regstr")) CK(skreg(c, "regstr", r));
  out = vnum(clab_hllcnt(r.data(), r.size()));
  return true;
}

bool skchmrg(Ctx &c, Val &out) {
  std::vector<int> a, b;
  CK(skreg(c, "a", a));
  CK(skreg(c, "b", b));
  if (a.size() != b.size()) return ckstop(c, "both sketches must have the same number of registers");
  for (std::size_t i = 0; i < a.size(); ++i) a[i] = std::max(a[i], b[i]);
  out = vints(a);
  return true;
}
