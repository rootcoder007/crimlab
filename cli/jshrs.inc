// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// shrs() (R/shrs.R) for clab_jrun(): a data frame (an object of columns,
// with count naming one) or a vector of counts (an object gives the label
// column). Groups are sorted as aggregate() sorts them, the first by
// column varying fastest; text is ordered by its bytes.

// One column of a data frame as values that can be compared and written
// back: numbers (NaN missing) or text (miss).
struct Cv {
  bool num = true;
  std::vector<double> d;
  std::vector<std::string> s;
  std::vector<char> miss;
  std::size_t size() const { return num ? d.size() : s.size(); }
  bool na(std::size_t i) const { return num ? std::isnan(d[i]) : miss[i] != 0; }
  Val get(std::size_t i) const {
    if (na(i)) return Val();
    return num ? vnum(d[i]) : vstr(s[i]);
  }
};

bool cvof(const Val &o, Cv &v) {
  v = Cv();
  if (isnum(o)) {
    v.d = dbl(o);
    return true;
  }
  v.num = false;
  const std::vector<Val> one = {o};
  const std::vector<Val> &es = o.t == Val::ARR ? o.a : one;
  for (const Val &e : es) {
    std::string t;
    if (e.t == Val::NUL) {
      v.s.emplace_back();
      v.miss.push_back(1);
    } else if (lab(e, t)) {
      v.s.push_back(t);
      v.miss.push_back(0);
    } else {
      return false;
    }
  }
  return true;
}

int cvcmp(const Cv &v, std::size_t i, std::size_t j) {
  if (v.num) return v.d[i] < v.d[j] ? -1 : v.d[i] > v.d[j] ? 1 : 0;
  return std::strcmp(v.s[i].c_str(), v.s[j].c_str());
}

bool shfin(Ctx &c, std::vector<double> &cnt, double den, Val &out, Val keys) {
  double conf = 0.95;
  if (c.has("conflvl")) conf = val1(c["conflvl"]);
  if (std::isnan(conf) || conf <= 0 || conf >= 1) return c.fail("`conflvl` must be one number strictly between 0 and 1.");
  for (double x : cnt) {
    if (x < 0) return c.fail("`count` must be non-negative.");
  }
  if (std::isnan(den) || den <= 0) return c.fail("the total must be positive.");
  for (double x : cnt) {
    if (x > den) return c.fail("a count is larger than the total, so these are not shares of it.");
  }
  const std::size_t k = cnt.size();
  std::vector<double> sh(k, NA), lo(k, NA), hi(k, NA);
  for (std::size_t i = 0; i < k; ++i) {
    if (!std::isnan(cnt[i])) clab_shrs(&cnt[i], 1, den, conf, &sh[i], &lo[i], &hi[i]);
  }
  out = keys.t == Val::OBJ ? keys : vobj();
  put(out, "count", vnums(cnt));
  put(out, "total", vnums(std::vector<double>(k, den)));
  put(out, "share", vnums(sh));
  put(out, "lower", vnums(lo));
  put(out, "upper", vnums(hi));
  return true;
}

bool shtot(Ctx &c, const char *also, double &t) {
  t = val1(c["total"]);
  if (!std::isfinite(t) || t < 0) {
    return c.fail("`total` must be NULL, %sone finite non-negative number.", also);
  }
  return true;
}

bool shrs(Ctx &c, Val &out) {
  const Val &x = c["x"];
  bool df = false;
  for (const Val &e : x.a) df = df || (x.t == Val::OBJ && (e.t == Val::ARR || e.t == Val::OBJ));
  for (const char *a : {"count", "by"}) {
    if (!df && c.has(a)) return c.fail("unused argument(s): %s.", a);
  }
  if (!df) {
    Dv d;
    bool fin = readd(x, d) == NUM && x.t != Val::BOOL;
    bool all = true;
    for (double v : d.v) {
      fin = fin && !std::isinf(v);
      all = all && std::isnan(v);
    }
    if (!fin || all) return c.fail("`x` must be finite counts (a numeric vector) or a data frame of them.");
    std::vector<double> cnt = dbl(x);
    double den = 0;
    if (c.has("total")) {
      CK(shtot(c, "", den));
    } else {
      for (double v : cnt) den += std::isnan(v) ? 0 : v;
    }
    Val keys;
    if (x.t == Val::OBJ) {
      keys = vobj();
      put(keys, "label", vstrs(x.k));
    }
    return shfin(c, cnt, den, out, keys);
  }
  const Val &cn = c["count"];
  if (cn.t != Val::STR) return c.fail("`count` must name the column of counts.");
  const Val &cv = col(x, cn.s.c_str());
  if (cv.t == Val::NUL) return c.fail("`count` column not found: %s.", cn.s.c_str());
  if (!isnum(cv)) return c.fail("`%s` must be numeric.", cn.s.c_str());
  const std::vector<double> y = dbl(cv);
  const std::size_t n = y.size();
  std::vector<std::string> by;
  if (c.has("by") && !labs(c["by"], by)) return c.fail("`by` column not found: %s.", "by");
  std::vector<std::string> miss;
  for (const std::string &b : by) {
    if (col(x, b.c_str()).t == Val::NUL) miss.push_back(b);
  }
  if (!miss.empty()) {
    std::string m;
    for (const std::string &b : miss) m += (m.empty() ? "" : ", ") + b;
    return c.fail("`by` column not found: %s.", m.c_str());
  }
  std::vector<Cv> g(by.size());
  for (std::size_t j = 0; j < by.size(); ++j) {
    if (!cvof(col(x, by[j].c_str()), g[j]) || g[j].size() != n) return c.fail("`by` column not found: %s.", by[j].c_str());
  }
  std::vector<double> cnt;
  Val keys;
  if (by.empty()) {
    double s = 0;
    for (double v : y) s += std::isnan(v) ? 0 : v;
    cnt.push_back(s);
  } else {
    std::vector<std::size_t> r;
    for (std::size_t i = 0; i < n; ++i) {
      bool ok = true;
      for (const Cv &v : g) ok = ok && !v.na(i);
      if (ok) r.push_back(i);
    }
    auto less = [&g](std::size_t a, std::size_t b) {
      for (std::size_t j = g.size(); j-- > 0;) {
        const int k = cvcmp(g[j], a, b);
        if (k) return k < 0;
      }
      return false;
    };
    std::stable_sort(r.begin(), r.end(), less);
    std::vector<std::size_t> first;
    for (std::size_t q = 0; q < r.size(); ++q) {
      if (q == 0 || less(r[q - 1], r[q])) {
        first.push_back(r[q]);
        cnt.push_back(0);
      }
      if (!std::isnan(y[r[q]])) cnt.back() += y[r[q]];
    }
    keys = vobj();
    for (std::size_t j = 0; j < by.size(); ++j) {
      Val a = varr();
      for (std::size_t i : first) a.a.push_back(g[j].get(i));
      put(keys, by[j].c_str(), a);
    }
  }
  double den = 0;
  if (!c.has("total")) {
    for (double v : cnt) den += v;
  } else if (c["total"].t == Val::STR) {
    const Val &t = col(x, c["total"].s.c_str());
    if (t.t == Val::NUL) return c.fail("`total` column not found: %s.", c["total"].s.c_str());
    for (double v : dbl(t)) den += std::isnan(v) ? 0 : v;
  } else {
    CK(shtot(c, "a column name or ", den));
  }
  return shfin(c, cnt, den, out, keys);
}
