// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Region-map checks (R/rgnmap.R) for clab_jrun(): maps are data frames
// (objects of columns); text is ordered by its bytes.

// A column as.character(): text, numbers as R writes them, null missing.
bool rmtxt(const Val &v, std::vector<std::string> &s, std::vector<char> &na) {
  s.clear();
  na.clear();
  const std::vector<Val> one = {v};
  for (const Val &e : v.t == Val::ARR ? v.a : one) {
    std::string t;
    if (e.t == Val::NUL) {
      s.emplace_back("NA");
      na.push_back(1);
    } else if (lab(e, t)) {
      s.push_back(t);
      na.push_back(0);
    } else {
      return false;
    }
  }
  return true;
}

std::string rmtrim(const std::string &s) {
  const char *ws = " \t\n\r";
  const std::size_t a = s.find_first_not_of(ws);
  if (a == std::string::npos) return "";
  return s.substr(a, s.find_last_not_of(ws) - a + 1);
}

bool rmcol(Ctx &c, const char *arg, const char *nmarg, std::vector<std::string> &s, std::vector<char> &na,
           const Val **raw = nullptr) {
  const Val &d = c[arg];
  if (d.t != Val::OBJ) return c.fail("`%s` must be a data frame.", arg);
  std::string nm;
  if (c[nmarg].t != Val::STR) {
    std::vector<std::string> v;
    labs(c[nmarg], v);
    for (std::size_t i = 0; i < v.size(); ++i) nm += (i ? ", " : "") + v[i];
    return c.fail("`%s` has no column `%s`.", arg, nm.c_str());
  }
  nm = c[nmarg].s;
  const Val &v = col(d, nm.c_str());
  if (v.t == Val::NUL || !rmtxt(v, s, na)) return c.fail("`%s` has no column `%s`.", arg, nm.c_str());
  // every column of a data frame has one length: the callers index one
  // column by the rows of another
  for (const Val &e : d.a) {
    if ((e.t == Val::ARR ? e.a.size() : 1) != s.size()) return c.fail("the columns of `%s` differ in length.", arg);
  }
  if (raw) *raw = &v;
  return true;
}

bool rgnmapintg(Ctx &c, Val &out) {
  std::vector<std::string> u, r;
  std::vector<char> un, rn;
  CK(rmcol(c, "map", "unit", u, un));
  CK(rmcol(c, "map", "reg", r, rn));
  auto blank = [](const std::string &s, char na) { return na || rmtrim(s).empty(); };
  double dup = 0, gap = 0, unk = 0;
  std::vector<std::string> seen;
  for (std::size_t i = 0; i < u.size(); ++i) {
    if (blank(u[i], un[i])) continue;
    if (std::find(seen.begin(), seen.end(), u[i]) != seen.end()) dup += 1;
    else seen.push_back(u[i]);
  }
  for (std::size_t i = 0; i < r.size(); ++i) gap += blank(r[i], rn[i]);
  const bool hr = c.has("regs");
  if (hr) {
    std::vector<std::string> g;
    std::vector<char> gn;
    rmtxt(c["regs"], g, gn);
    for (std::string &s : g) s = rmtrim(s);
    for (std::size_t i = 0; i < r.size(); ++i) {
      if (!blank(r[i], rn[i]) && std::find(g.begin(), g.end(), rmtrim(r[i])) == g.end()) unk += 1;
    }
  }
  std::vector<std::string> chk = {"units assigned more than one region", "units assigned no region"};
  std::vector<double> obs = {dup, gap};
  if (hr) {
    chk.push_back("region codes not in the reference geography");
    obs.push_back(unk);
  }
  Frame f;
  f.str("check", chk);
  f.num("observed", obs);
  f.num("expected", std::vector<double>(obs.size(), 0.0));
  std::vector<int> p;
  for (double v : obs) p.push_back(v == 0);
  f.lgl("pass", p.data(), p.size());
  out = f.o;
  return true;
}

// all.equal(a, b) of two numbers at its default tolerance
bool aeq(double a, double b) {
  if (std::isnan(a) || std::isnan(b)) return std::isnan(a) && std::isnan(b);
  const double tol = 1.5e-8;
  double xy = std::fabs(a - b);
  if (std::isfinite(a) && std::fabs(a) > tol) xy /= std::fabs(a);
  return !(xy > tol);
}

bool rgnmapcmp(Ctx &c, Val &out) {
  std::vector<std::string> pu, ou;
  std::vector<char> pn, on;
  CK(rmcol(c, "pub", "unit", pu, pn));
  CK(rmcol(c, "obs", "unit", ou, on));
  for (std::string &s : pu) s = rmtrim(s);
  for (std::string &s : ou) s = rmtrim(s);
  const Val &P = c["pub"], &O = c["obs"];
  std::vector<std::string> both, cols;
  for (const std::string &k : P.k) {
    if (std::find(O.k.begin(), O.k.end(), k) != O.k.end()) both.push_back(k);
  }
  if (c.has("cols")) {
    labs(c["cols"], cols);
  } else {
    for (const std::string &k : both) {
      if (k != c["unit"].s) cols.push_back(k);
    }
  }
  std::vector<std::string> miss;
  for (const std::string &k : cols) {
    if (std::find(both.begin(), both.end(), k) == both.end()) miss.push_back(k);
  }
  if (!miss.empty()) {
    std::string m;
    for (const std::string &k : miss) m += (m.empty() ? "" : ", ") + k;
    return c.fail("both frames need the compared columns; missing: %s.", m.c_str());
  }
  auto uniq = [](const std::vector<std::string> &v) {
    std::vector<std::string> o;
    for (const std::string &s : v) {
      if (std::find(o.begin(), o.end(), s) == o.end()) o.push_back(s);
    }
    return o;
  };
  auto in = [](const std::vector<std::string> &v, const std::string &s) {
    return std::find(v.begin(), v.end(), s) != v.end();
  };
  std::vector<std::string> onlyp, onlyo, shared, all = pu;
  for (const std::string &s : uniq(pu)) (in(ou, s) ? shared : onlyp).push_back(s);
  for (const std::string &s : uniq(ou)) {
    if (!in(pu, s)) onlyo.push_back(s);
  }
  all.insert(all.end(), ou.begin(), ou.end());
  std::vector<std::string> cn = {"rows"}, first;
  std::vector<double> cells = {static_cast<double>(uniq(all).size())},
                      mis = {static_cast<double>(onlyp.size() + onlyo.size())};
  first.push_back(!onlyp.empty() ? "published only: " + onlyp[0]
                                 : !onlyo.empty() ? "recomputed only: " + onlyo[0] : std::string());
  auto pos = [](const std::vector<std::string> &v, const std::string &s) {
    return static_cast<std::size_t>(std::find(v.begin(), v.end(), s) - v.begin());
  };
  for (const std::string &k : cols) {
    const Val &a = col(P, k.c_str()), &b = col(O, k.c_str());
    std::vector<std::string> as, bs;
    std::vector<char> an, bn;
    rmtxt(a, as, an);
    rmtxt(b, bs, bn);
    const bool num = isnum(a) && isnum(b);
    const std::vector<double> ad = num ? dbl(a) : std::vector<double>(), bd = num ? dbl(b) : std::vector<double>();
    std::size_t nb = 0;
    std::string f1;
    for (const std::string &s : shared) {
      const std::size_t i = pos(pu, s), j = pos(ou, s);
      bool same;
      if (num) {
        same = aeq(ad[i], bd[j]);
      } else if (an[i] || bn[j]) {
        same = an[i] && bn[j];
      } else {
        same = rmtrim(as[i]) == rmtrim(bs[j]);
      }
      if (!same) {
        if (nb == 0) f1 = s + ": published " + as[i] + ", recomputed " + bs[j];
        ++nb;
      }
    }
    cn.push_back(k);
    cells.push_back(static_cast<double>(shared.size()));
    mis.push_back(static_cast<double>(nb));
    first.push_back(f1);
  }
  Frame f;
  f.str("column", cn);
  f.num("cells", cells);
  f.num("mismatched", mis);
  f.str("first", first);
  out = f.o;
  return true;
}

bool rgnmapalt(Ctx &c, Val &out) {
  std::vector<std::string> u, r;
  std::vector<char> un, rn;
  CK(rmcol(c, "map", "unit", u, un));
  CK(rmcol(c, "map", "reg", r, rn));
  for (std::string &s : u) s = rmtrim(s);
  for (std::string &s : r) s = rmtrim(s);
  std::vector<std::string> rk, rv;
  const Val &ro = c["route"];
  bool isdf = ro.t == Val::OBJ && !ro.a.empty();
  for (const Val &e : ro.a) isdf = isdf && e.t == Val::ARR;
  if (isdf) {
    std::vector<std::string> a, b;
    std::vector<char> an, bn;
    CK(rmcol(c, "route", "reg", a, an));
    CK(rmcol(c, "route", "unit", b, bn));
    for (std::size_t i = 0; i < a.size(); ++i) {
      if (an[i] || rmtrim(a[i]).empty()) continue;
      rk.push_back(bn[i] ? "NA" : rmtrim(b[i]));
      rv.push_back(rmtrim(a[i]));
    }
  } else if (ro.t == Val::OBJ) {
    for (std::size_t i = 0; i < ro.a.size(); ++i) {
      std::string t;
      if (!lab(ro.a[i], t) || rmtrim(t).empty()) continue;
      rk.push_back(rmtrim(ro.k[i]));
      rv.push_back(rmtrim(t));
    }
  } else {
    return c.fail("`route` must be named by unit, or be a data frame.");
  }
  std::vector<std::string> known;
  std::vector<char> kn;
  rmtxt(c["known"], known, kn);
  if (!c.has("known")) known.clear();
  for (std::string &s : known) s = rmtrim(s);
  std::vector<std::size_t> bad;
  std::vector<std::string> sec;
  for (std::size_t i = 0; i < u.size(); ++i) {
    const auto it = std::find(rk.begin(), rk.end(), u[i]);
    if (it == rk.end()) continue;
    const std::string &s = rv[static_cast<std::size_t>(it - rk.begin())];
    if (s != r[i]) {
      bad.push_back(i);
      sec.push_back(s);
    }
  }
  std::vector<std::size_t> o(bad.size());
  std::vector<int> kw(bad.size());
  for (std::size_t q = 0; q < bad.size(); ++q) {
    o[q] = q;
    kw[q] = std::find(known.begin(), known.end(), u[bad[q]]) != known.end();
  }
  std::stable_sort(o.begin(), o.end(), [&](std::size_t a, std::size_t b) {
    if (kw[a] != kw[b]) return kw[a] < kw[b];
    return std::strcmp(u[bad[a]].c_str(), u[bad[b]].c_str()) < 0;
  });
  std::vector<std::string> fu, fp, fs;
  std::vector<int> fk;
  for (std::size_t q : o) {
    fu.push_back(u[bad[q]]);
    fp.push_back(r[bad[q]]);
    fs.push_back(sec[q]);
    fk.push_back(kw[q]);
  }
  Frame f;
  f.str("unit", fu);
  f.str("primary", fp);
  f.str("second", fs);
  f.lgl("known", fk.data(), fk.size());
  out = f.o;
  return true;
}
