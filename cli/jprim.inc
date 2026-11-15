// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The summary statistics (R/prim.R, src/prim.cpp) for clab_jrun(), with the
// checks and words of the R functions.

// .pnum(): numbers, true and false, null (NA); left out: no values.
bool pnum(Ctx &c, const char *nm, std::vector<double> &v) {
  v.clear();
  const Val &o = c[nm];
  if (o.t == Val::NUL) return true;
  auto one1 = [&v](const Val &e) {
    double x = 0.0;
    if (e.t == Val::NUM) v.push_back(e.n);
    else if (e.t == Val::BOOL) v.push_back(e.b ? 1.0 : 0.0);
    else if (e.t == Val::NUL) v.push_back(NA);
    else if (snum(e, x)) v.push_back(x);
    else return false;
    return true;
  };
  bool ok = true;
  if (o.t == Val::ARR || o.t == Val::OBJ) {
    for (const Val &e : o.a) ok = ok && one1(e);
  } else {
    ok = one1(o);
  }
  return ok || c.fail("`%s` must be numbers (or TRUE and FALSE).", nm);
}

// .pone(): one finite number from lo (to hi), whole when asked; the
// default when left out.
bool pone(Ctx &c, const char *nm, double def, double lo, double hi, bool whole, double &x) {
  x = c.has(nm) ? val1(c[nm]) : def;
  if (c.has(nm) && c[nm].t == Val::BOOL) x = NaN;
  if (std::isfinite(x) && x >= lo && x <= hi && (!whole || x == std::floor(x))) return true;
  const std::string to = std::isfinite(hi) && !whole ? " to " + rfmt(hi) : "";
  return c.fail("`%s` must be %s from %s%s.", nm, whole ? "one whole number" : "one number", rfmt(lo).c_str(),
                to.c_str());
}

bool psame(Ctx &c, const std::vector<double> &x, const std::vector<double> &y, const char *a, const char *b) {
  return x.size() == y.size() || c.fail("`%s` and `%s` must have the same length.", a, b);
}

bool pseed(Ctx &c, double def, double &s) {
  s = c.has("seed") ? val1(c["seed"]) : def;
  if (c.has("seed") && c["seed"].t == Val::BOOL) s = NaN;
  return (s >= 0 && s == std::floor(s) && s < 9007199254740992.0) ||
         c.fail("`seed` must be one whole number from 0, below 2^53.");
}

int ilen(const std::vector<double> &v) { return static_cast<int>(v.size()); }

bool cmn(Ctx &c, Val &out) {
  std::vector<double> x;
  CK(pnum(c, "x", x));
  out = vnum(clab_mean(x.data(), ilen(x)));
  return true;
}

bool cvari(Ctx &c, Val &out) {
  std::vector<double> x;
  CK(pnum(c, "x", x));
  out = vnum(clab_var(x.data(), ilen(x), 1));
  return true;
}

bool ccorr(Ctx &c, Val &out) {
  std::vector<double> x, y;
  CK(pnum(c, "x", x));
  CK(pnum(c, "y", y));
  CK(psame(c, x, y, "x", "y"));
  out = vnum(clab_cor(x.data(), y.data(), ilen(x)));
  return true;
}

bool csdev(Ctx &c, Val &out) {
  std::vector<double> x;
  CK(pnum(c, "x", x));
  double d = 0;
  CK(pone(c, "ddof", 1, 0, INT_MAX, true, d));
  out = vnum(std::sqrt(clab_var(x.data(), ilen(x), static_cast<int>(d))));
  return true;
}

bool cdistm(Ctx &c, Val &out) {
  std::vector<double> a, b;
  CK(pnum(c, "a", a));
  CK(pnum(c, "b", b));
  CK(psame(c, a, b, "a", "b"));
  out = vnum(clab_edist(a.data(), b.data(), ilen(a)));
  return true;
}

bool pdens(Ctx &c, Val &out, int lg) {
  std::vector<double> x;
  CK(pnum(c, "x", x));
  const double m = c.has("mean") ? val1(c["mean"]) : 0, s = c.has("sd") ? val1(c["sd"]) : 1;
  if (!std::isfinite(m)) return c.fail("`mean` must be one finite number.");
  if (!std::isfinite(s) || s <= 0) return c.fail("`sd` must be one positive finite number.");
  std::vector<double> o(x.size());
  clab_dnorm(x.data(), ilen(x), m, s, lg, o.data());
  out = unbox(vnums(o));
  return true;
}
bool cnormpdf(Ctx &c, Val &out) { return pdens(c, out, 0); }
bool cnormlpdf(Ctx &c, Val &out) { return pdens(c, out, 1); }

bool cmom(Ctx &c, Val &out) {
  std::vector<double> x;
  CK(pnum(c, "x", x));
  double m[4];
  clab_moms(x.data(), ilen(x), m);
  static const char *const nm[] = {"mean", "variance", "skewness", "kurtosis"};
  out = nvec(4, nm, m);
  return true;
}

// formatC(100 p, format = "fg", digits = 7) followed by %
std::string pctnm(double p) {
  char b[64];
  std::snprintf(b, sizeof b, "%.7g", 100 * p);
  if (std::strchr(b, 'e') != nullptr) {
    std::snprintf(b, sizeof b, "%.15f", 100 * p);
    std::string s = b;
    while (s.back() == '0') s.pop_back();
    if (s.back() == '.') s.pop_back();
    return s + "%";
  }
  return std::string(b) + "%";
}

bool cqnt(Ctx &c, Val &out) {
  std::vector<double> x, p;
  CK(pnum(c, "x", x));
  if (c.has("probs")) {
    CK(pnum(c, "probs", p));
  } else {
    p = {0, 0.25, 0.5, 0.75, 1};
  }
  if (p.empty()) return c.fail("`probs` is empty.");
  for (double v : p) {
    if (std::isnan(v) || v < 0 || v > 1) return c.fail("`probs` must lie in [0, 1].");
  }
  std::vector<double> q(p.size());
  clab_quant(x.data(), ilen(x), p.data(), ilen(p), q.data());
  out = vobj();
  for (std::size_t k = 0; k < p.size(); ++k) put(out, pctnm(p[k]).c_str(), vnum(q[k]));
  return true;
}

bool cmed(Ctx &c, Val &out) {
  std::vector<double> x;
  CK(pnum(c, "x", x));
  out = vnum(clab_med(x.data(), ilen(x)));
  return true;
}

bool cmad(Ctx &c, Val &out) {
  std::vector<double> x;
  CK(pnum(c, "x", x));
  double k = 0;
  CK(pone(c, "const", 1.4826, 0, Inf, false, k));
  out = vnum(clab_mad(x.data(), ilen(x), k));
  return true;
}

void pq2(const std::vector<double> &x, double *q) {
  const double p[2] = {0.25, 0.75};
  clab_quant(x.data(), ilen(x), p, 2, q);
}

bool ciqr(Ctx &c, Val &out) {
  std::vector<double> x;
  CK(pnum(c, "x", x));
  double q[2];
  pq2(x, q);
  out = vnum(q[1] - q[0]);
  return true;
}

bool ctkyfnc(Ctx &c, Val &out) {
  std::vector<double> x;
  CK(pnum(c, "x", x));
  double k = 0, q[2];
  CK(pone(c, "k", 1.5, 0, Inf, false, k));
  pq2(x, q);
  const double r[2] = {q[0] - k * (q[1] - q[0]), q[1] + k * (q[1] - q[0])};
  static const char *const nm[] = {"lower", "upper"};
  out = nvec(2, nm, r);
  return true;
}

bool ptrim(Ctx &c, Val &out, bool w) {
  std::vector<double> x;
  CK(pnum(c, "x", x));
  double t = 0;
  CK(pone(c, "trim", 0.1, 0, 0.5, false, t));
  out = vnum(w ? clab_wins(x.data(), ilen(x), t) : clab_trim(x.data(), ilen(x), t));
  return true;
}
bool ctrmmean(Ctx &c, Val &out) { return ptrim(c, out, false); }
bool cwinsmean(Ctx &c, Val &out) { return ptrim(c, out, true); }

bool cwtd(Ctx &c, Val &out) {
  std::vector<double> x, w;
  CK(pnum(c, "x", x));
  CK(pnum(c, "w", w));
  CK(psame(c, x, w, "x", "w"));
  bool bad = false, pos = false;
  for (double v : w) {
    bad = bad || v < 0 || std::isinf(v);
    pos = pos || v > 0;
  }
  if (bad || !pos) return c.fail("`w` must be non-negative and finite, not all zero.");
  double r[2];
  clab_wtd(x.data(), w.data(), ilen(x), r);
  static const char *const nm[] = {"mean", "variance"};
  out = nvec(2, nm, r);
  return true;
}

bool cspr(Ctx &c, Val &out) {
  std::vector<double> x, y;
  CK(pnum(c, "x", x));
  CK(pnum(c, "y", y));
  CK(psame(c, x, y, "x", "y"));
  out = vnum(clab_spear(x.data(), y.data(), ilen(x)));
  return true;
}

bool cmidrk(Ctx &c, Val &out) {
  std::vector<double> x, k, r;
  CK(pnum(c, "x", x));
  for (double v : x) {
    if (!std::isnan(v)) k.push_back(v);
  }
  r.resize(k.size());
  clab_midrk(k.data(), ilen(k), r.data());
  std::vector<double> o(x.size(), NA);
  for (std::size_t i = 0, j = 0; i < x.size(); ++i) {
    if (!std::isnan(x[i])) o[i] = r[j++];
  }
  out = unbox(vnums(o));
  return true;
}

// x: a matrix (an array of rows) or a data frame (an object of columns).
bool ccovm(Ctx &c, Val &out) {
  const Val &x = c["x"];
  std::vector<double> v;
  std::vector<std::string> cn;
  int nr = 0, nc = 0;
  if (x.t == Val::OBJ) {
    for (const Val &e : x.a) {
      if (!isnum(e)) return c.fail("every column of `x` must be numeric.");
    }
    cn = x.k;
    nc = static_cast<int>(x.a.size());
    nr = nc ? static_cast<int>(vlen(x.a[0])) : 0;
    for (const Val &e : x.a) {
      if (static_cast<int>(vlen(e)) != nr) return c.fail("`x` must be a numeric matrix or data frame.");
      const std::vector<double> d = dbl(e);
      v.insert(v.end(), d.begin(), d.end());
    }
  } else if (x.t != Val::ARR || x.a.empty() || x.a[0].t != Val::ARR || !nmat(c, x, "x", nr, nc, v, true)) {
    c.err.clear();
    return c.fail("`x` must be a numeric matrix or data frame.");
  }
  std::vector<double> o(static_cast<std::size_t>(nc) * static_cast<std::size_t>(nc));
  clab_covm(v.data(), nr, nc, o.data());
  out = varr();
  for (int i = 0; i < nc; ++i) {
    Val r = varr();
    for (int j = 0; j < nc; ++j) r.a.push_back(vnum(o[static_cast<std::size_t>(j) * nc + i]));
    out.a.push_back(r);
  }
  return true;
}

bool cbootmean(Ctx &c, Val &out) {
  std::vector<double> x;
  CK(pnum(c, "x", x));
  double b = 0, s = 0;
  CK(pone(c, "B", 1000, 1, INT_MAX, true, b));
  CK(pseed(c, 42, s));
  std::vector<double> o(static_cast<std::size_t>(b));
  CK(kchk(c, clab_bmean(x.data(), ilen(x), static_cast<int>(b), s, o.data())));
  out = unbox(vnums(o));
  return true;
}

bool cipw(Ctx &c, Val &out) {
  std::vector<double> t, p;
  CK(pnum(c, "treat", t));
  CK(pnum(c, "ps", p));
  CK(psame(c, t, p, "treat", "ps"));
  double lo = 0, hi = 0;
  CK(pone(c, "trimlo", 0.01, 0, 1, false, lo));
  CK(pone(c, "trimhi", 0.99, 0, 1, false, hi));
  if (!(lo > 0 && hi < 1 && lo < hi)) return c.fail("need 0 < `trimlo` < `trimhi` < 1.");
  for (double v : t) {
    if (!std::isnan(v) && v != 0 && v != 1) return c.fail("`treat` must be coded 0/1.");
  }
  for (double v : p) {
    if (v < 0 || v > 1) return c.fail("`ps` must lie in [0, 1].");
  }
  std::vector<double> o(t.size(), NA);
  for (std::size_t i = 0; i < t.size(); ++i) {
    if (!std::isnan(t[i]) && !std::isnan(p[i])) clab_ipw(&t[i], &p[i], 1, lo, hi, &o[i]);
  }
  out = unbox(vnums(o));
  return true;
}

bool cgamcdf(Ctx &c, Val &out) {
  const double a = val1(c["shape"]);
  if (c["shape"].t == Val::BOOL || std::isnan(a) || a <= 0) return c.fail("`shape` must be one positive number.");
  std::vector<double> x;
  CK(pnum(c, "x", x));
  for (double &v : x) {
    if (!std::isnan(v)) v = clab_pgamma(v, a, 1);
  }
  out = unbox(vnums(x));
  return true;
}

bool cunif(Ctx &c, Val &out) {
  double n = 0, s = 0;
  CK(pone(c, "n", 0, 0, INT_MAX, true, n));
  CK(pseed(c, 0, s));
  std::vector<double> o(static_cast<std::size_t>(n));
  clab_unif(s, static_cast<int>(n), o.data());
  out = unbox(vnums(o));
  return true;
}
