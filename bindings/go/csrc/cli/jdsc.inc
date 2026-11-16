// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Data description (src/dsc.cpp), with the checks of R/dsc.R in their
// words. A data frame is an object of columns; a column of numbers and
// nulls is numeric, one with text character, one of true and false
// logical. Text is sorted in byte order, where R sorts it in the
// session's collation: tables of mixed-case labels can order their ties
// differently.

// A column: its name, its R type (0 numeric, 1 character, 2 logical) and
// its values (NA: null).
struct Dcol {
  std::string nm;
  int ty = 0;
  std::vector<double> x;       // numeric; logical as 0/1
  std::vector<std::string> s;  // as.character()
  std::vector<char> na;
  const Val *v = nullptr;
};

inline bool dframe(Ctx &c, const Val &d, std::vector<Dcol> &cs, std::size_t &n, const char *msg) {
  cs.clear();
  n = 0;
  if (d.t != Val::OBJ) return ckstop(c, "%s", msg);
  for (const Val &v : d.a) n = std::max(n, vlen(v));
  for (std::size_t j = 0; j < d.a.size(); ++j) {
    Dcol k;
    k.nm = d.k[j];
    k.v = &d.a[j];
    const Val &v = d.a[j];
    bool str = false, lg = false, num = false;
    for (std::size_t i = 0; i < vlen(v); ++i) {
      const Val &e = v.t == Val::ARR ? v.a[i] : v;
      double z = 0;
      if (e.t == Val::STR && !snum(e, z)) str = true;
      else if (e.t == Val::BOOL) lg = true;
      else if (e.t != Val::NUL) num = true;
    }
    k.ty = str ? 1 : lg && !num ? 2 : 0;
    for (std::size_t i = 0; i < n; ++i) {
      const Val &e = v.t == Val::ARR ? (i < v.a.size() ? v.a[i] : nul) : v;
      double z = NA;
      std::string t = "NA";
      bool isn = e.t == Val::NUL;
      if (e.t == Val::NUM) {
        z = e.n;
        t = rfmt(z);
      } else if (e.t == Val::BOOL) {
        z = e.b;
        t = e.b ? "TRUE" : "FALSE";
      } else if (e.t == Val::STR) {
        if (snum(e, z) && k.ty != 1) {
          t = e.s;
          isn = std::isnan(z);
        } else {
          t = e.s;
          z = NA;
        }
      }
      if (k.ty == 0 && e.t == Val::NUM) isn = std::isnan(e.n);
      k.x.push_back(isn ? NA : z);
      k.s.push_back(isn ? "NA" : t);
      k.na.push_back(isn);
    }
    cs.push_back(k);
  }
  return true;
}

// .dfnm(): named columns, no name twice.
inline bool dfnm(Ctx &c, const std::vector<Dcol> &cs) {
  for (std::size_t j = 0; j < cs.size(); ++j) {
    if (cs[j].nm.empty()) {
      return ckstop(c, "every column of `data` must have a name (names(data) has an empty one)");
    }
    for (std::size_t k = 0; k < j; ++k) {
      if (cs[k].nm == cs[j].nm) return ckstop(c, "`data` has two columns named %s: rename one", cs[j].nm.c_str());
    }
  }
  return true;
}

// .ckpos(): one positive whole number (def when left out).
inline bool ckposj(Ctx &c, const char *nm, int def, int &v) {
  v = def;
  if (!c.has(nm)) return true;
  const Val &o = c[nm];
  const Val *e = o.t == Val::ARR && o.a.size() == 1 ? &o.a[0] : &o;
  double x = NaN;
  if (e->t == Val::BOOL) return ckstop(c, "`%s` must be a whole number, not %s", nm, e->b ? "TRUE" : "FALSE");
  if (e->t == Val::NUM) x = e->n;
  else if (e->t == Val::STR) {
    char *end = nullptr;
    x = std::strtod(e->s.c_str(), &end);
    if (end == e->s.c_str() || *end != 0) {
      return ckstop(c, "`%s` must be numeric, not text such as \"%s\"", nm, e->s.c_str());
    }
  }
  if (std::isfinite(x) && x != std::floor(x)) {
    return ckstop(c, "`%s` must be a whole number, not %s", nm, rfmt(x).c_str());
  }
  if (o.t == Val::ARR && o.a.size() != 1) x = NaN;
  if (!(x >= 1) || x > INT_MAX) return ckstop(c, "`%s` must be a single positive integer", nm);
  v = static_cast<int>(x);
  return true;
}

inline bool flagj(Ctx &c, const char *nm, bool def, bool &b) {
  b = def;
  if (!c.has(nm)) return true;
  const Val &o = c[nm];
  const Val *e = o.t == Val::ARR && o.a.size() == 1 ? &o.a[0] : &o;
  if (e->t != Val::BOOL) return ckstop(c, "`%s` must be TRUE or FALSE", nm);
  b = e->b;
  return true;
}

inline bool uniok() {
  const char *a = std::getenv("CRIMLAB_ASCII_ONLY");
  return a == nullptr || *a == 0;
}

// .spark(): the sketch of finite values x.
inline std::string spark(std::vector<double> x, int bins) {
  x.erase(std::remove_if(x.begin(), x.end(), [](double v) { return !std::isfinite(v); }), x.end());
  if (x.empty()) return std::string(static_cast<std::size_t>(bins), ' ');
  const bool u = uniok();
  if (*std::min_element(x.begin(), x.end()) == *std::max_element(x.begin(), x.end())) {
    return std::string(static_cast<std::size_t>(bins / 2), ' ') + (u ? "█" : "#") +
           std::string(static_cast<std::size_t>(bins - bins / 2 - 1), ' ');
  }
  std::vector<int> cnt(static_cast<std::size_t>(bins));
  clab_hbin(x.data(), static_cast<int>(x.size()), bins, cnt.data());
  static const char *const uni[] = {"▁", "▂", "▃", "▄", "▅", "▆", "▇", "█"};
  static const char *const asc[] = {".", ":", "-", "=", "#"};
  const int nl = u ? 8 : 5;
  const int mx = *std::max_element(cnt.begin(), cnt.end());
  std::string r;
  for (int k : cnt) {
    if (k == 0) {
      r += " ";
      continue;
    }
    const int i = std::max(1, static_cast<int>(std::ceil(static_cast<double>(k) / mx * nl)));
    r += u ? uni[i - 1] : asc[i - 1];
  }
  return r;
}

bool inlhist(Ctx &c, Val &out) {
  std::vector<double> x;
  const Val &o = c["x"];
  bool lg = false;
  const char *cl = o.t == Val::NUL ? nullptr : jcls(o, lg);
  if (cl != nullptr && !lg) return ckstop(c, "`x` must be numeric, got %s", cl);
  if (o.t != Val::NUL) {
    for (std::size_t i = 0; i < vlen(o); ++i) {
      const Val &e = o.t == Val::ARR ? o.a[i] : o;
      double z = NaN;
      if (e.t == Val::NUM) z = e.n;
      else if (e.t == Val::BOOL) z = e.b;
      else snum(e, z);
      x.push_back(z);
    }
  }
  int bins = 10;
  CK(ckposj(c, "bins", 10, bins));
  if (bins > kutil::maxcell) return c.fail("`bins` must be at most %.0f.", kutil::maxcell);
  out = vstr(spark(x, bins));
  return true;
}

// A frequency table from labels (NA: missing) as R's table() and sort() give it.
inline void freq(const std::vector<std::string> &s, const std::vector<char> &na, bool srt, int maxlvl, bool fold,
                 std::vector<std::string> &lv, std::vector<int> &lev, std::vector<int> &cnt) {
  lv.clear();
  for (std::size_t i = 0; i < s.size(); ++i) {
    if (!na[i]) lv.push_back(s[i]);
  }
  std::sort(lv.begin(), lv.end());
  lv.erase(std::unique(lv.begin(), lv.end()), lv.end());
  std::vector<int> code;
  for (std::size_t i = 0; i < s.size(); ++i) {
    if (!na[i]) code.push_back(static_cast<int>(std::lower_bound(lv.begin(), lv.end(), s[i]) - lv.begin()) + 1);
  }
  lev.assign(lv.size() + 1, 0);
  cnt.assign(lv.size() + 1, 0);
  const int k = clab_freq(code.data(), static_cast<int>(code.size()), static_cast<int>(lv.size()), srt, maxlvl, fold,
                          lev.data(), cnt.data());
  lev.resize(static_cast<std::size_t>(k));
  cnt.resize(static_cast<std::size_t>(k));
}

bool freqtbl(Ctx &c, Val &out) {
  bool srt = true;
  CK(flagj(c, "sort", true, srt));
  const Val &d = c["data"];
  std::vector<std::string> s;
  std::vector<char> na;
  if (d.t == Val::OBJ) {
    if (!c.has("col")) return ckstop(c, "`col` is required when `data` is a data frame");
    std::string col;
    CK(colnm(c, "col", "", col));
    std::vector<Dcol> cs;
    std::size_t n = 0;
    CK(dframe(c, d, cs, n, ""));
    const Dcol *k = nullptr;
    for (const Dcol &e : cs) {
      if (e.nm == col) k = &e;
    }
    if (k == nullptr) return ckstop(c, "no such column: %s", col.c_str());
    s = k->s;
    na = k->na;
  } else if (d.t != Val::NUL) {
    for (std::size_t i = 0; i < vlen(d); ++i) {
      const Val &e = d.t == Val::ARR ? d.a[i] : d;
      std::string t;
      na.push_back(!lab(e, t));
      s.push_back(t);
    }
  }
  int maxlvl = 25;
  CK(ckposj(c, "maxlvl", 25, maxlvl));
  std::vector<std::string> lv;
  std::vector<int> lev, cnt;
  freq(s, na, srt, maxlvl, true, lv, lev, cnt);
  int nna = 0;
  for (char e : na) nna += e;
  const double all = static_cast<double>(s.size()), valid = all - nna;
  Val v = varr(), nn = varr(), pc = varr(), pv = varr();
  auto row = [&](const Val &lab1, int k, bool missing) {
    v.a.push_back(lab1);
    nn.a.push_back(vnum(k));
    pc.a.push_back(vnum(100 * static_cast<double>(k) / all));
    pv.a.push_back(missing || !(valid > 0) ? Val() : vnum(100 * static_cast<double>(k) / valid));
  };
  if (!s.empty()) {
    for (std::size_t i = 0; i < lev.size(); ++i) {
      row(vstr(lev[i] > static_cast<int>(lv.size()) ? "(other)" : lv[static_cast<std::size_t>(lev[i] - 1)]), cnt[i],
          false);
    }
    if (nna > 0) row(Val(), nna, true);
  }
  out = vobj();
  put(out, "value", v);
  put(out, "n", nn);
  put(out, "pct", pc);
  put(out, "pctvalid", pv);
  return true;
}

inline bool methj(Ctx &c, int &m) {
  CK(ckarg(c, "meth", {"spearman", "pearson"}, m));
  return true;
}

// The numeric columns of data for the correlation functions.
inline bool corin(Ctx &c, std::vector<Dcol> &cs, std::vector<const Dcol *> &num, std::size_t &n) {
  CK(dframe(c, c["data"], cs, n, "`data` must be a data frame"));
  CK(dfnm(c, cs));
  for (const Dcol &k : cs) {
    if (k.ty == 0) num.push_back(&k);
  }
  if (num.size() < 2) return ckstop(c, "need at least two numeric columns to correlate");
  return true;
}

inline void corrun(const std::vector<const Dcol *> &num, std::size_t n, int meth, int minprs, std::vector<double> &r,
                   std::vector<int> &np) {
  const int p = static_cast<int>(num.size());
  std::vector<double> x(n * static_cast<std::size_t>(p));
  for (int j = 0; j < p; ++j) {
    for (std::size_t i = 0; i < n; ++i) x[static_cast<std::size_t>(j) * n + i] = num[static_cast<std::size_t>(j)]->x[i];
  }
  const std::size_t k = static_cast<std::size_t>(p) * (p - 1) / 2;
  r.assign(k, 0.0);
  np.assign(k, 0);
  clab_pcor(x.data(), static_cast<int>(n), p, meth, minprs, r.data(), np.data());
  for (std::size_t i = 0; i < k; ++i) {
    if (np[i] < minprs) r[i] = NA;
  }
}

bool cortbl(Ctx &c, Val &out) {
  int meth = 0;
  CK(methj(c, meth));
  std::vector<Dcol> cs;
  std::vector<const Dcol *> num;
  std::size_t n = 0;
  if (c["data"].t != Val::OBJ) return ckstop(c, "`data` must be a data frame");
  CK(corin(c, cs, num, n));
  int mp = 3;
  if (c.has("minprs")) {
    const double v = val1(c["minprs"]);
    if (std::isfinite(v) && v != std::floor(v)) {
      return ckstop(c, "`minprs` must be a whole number, not %s", rfmt(v).c_str());
    }
    if (!(v >= 2) || v > INT_MAX) return ckstop(c, "`minprs` must be a single integer of at least 2");
    mp = static_cast<int>(v);
  }
  std::vector<double> r;
  std::vector<int> np;
  corrun(num, n, meth, mp, r, np);
  Val x = varr(), y = varr();
  for (std::size_t i = 0; i < num.size(); ++i) {
    for (std::size_t j = i + 1; j < num.size(); ++j) {
      x.a.push_back(vstr(num[i]->nm));
      y.a.push_back(vstr(num[j]->nm));
    }
  }
  out = vobj();
  put(out, "x", x);
  put(out, "y", y);
  put(out, "correlation", vnums(r));
  put(out, "nprs", vints(np));
  return true;
}

bool topcors(Ctx &c, Val &out) {
  int meth = 0;
  CK(methj(c, meth));
  if (c["data"].t != Val::OBJ) return ckstop(c, "`data` must be a data frame");
  std::vector<Dcol> cs;
  std::vector<const Dcol *> num;
  std::size_t n = 0;
  CK(corin(c, cs, num, n));
  int top = 10;
  CK(ckposj(c, "n", 10, top));
  double mn = 0;
  if (c.has("minabs")) {
    mn = val1(c["minabs"]);
    if (!(mn >= 0 && mn <= 1)) return ckstop(c, "`minabs` must be a single value in [0, 1]");
  }
  std::vector<double> r;
  std::vector<int> np;
  corrun(num, n, meth, 3, r, np);
  std::vector<std::size_t> o;
  std::vector<std::pair<std::size_t, std::size_t>> pr;
  for (std::size_t i = 0; i < num.size(); ++i) {
    for (std::size_t j = i + 1; j < num.size(); ++j) pr.push_back({i, j});
  }
  for (std::size_t k = 0; k < r.size(); ++k) {
    if (!std::isnan(r[k]) && std::fabs(r[k]) >= mn) o.push_back(k);
  }
  std::stable_sort(o.begin(), o.end(), [&](std::size_t a, std::size_t b) { return std::fabs(r[a]) > std::fabs(r[b]); });
  if (o.size() > static_cast<std::size_t>(top)) o.resize(static_cast<std::size_t>(top));
  Val x = varr(), y = varr(), cr = varr(), ab = varr();
  for (std::size_t k : o) {
    x.a.push_back(vstr(num[pr[k].first]->nm));
    y.a.push_back(vstr(num[pr[k].second]->nm));
    cr.a.push_back(vnum(r[k]));
    ab.a.push_back(vnum(std::fabs(r[k])));
  }
  out = vobj();
  put(out, "x", x);
  put(out, "y", y);
  put(out, "correlation", cr);
  put(out, "abscor", ab);
  return true;
}

// The rows keep of the columns cs of a frame, as an object of columns.
inline Val subframe(const std::vector<Dcol> &cs, const std::vector<std::size_t> &rows,
                    const std::vector<char> &keepc) {
  Val o = vobj();
  for (std::size_t j = 0; j < cs.size(); ++j) {
    if (!keepc[j]) continue;
    Val a = varr();
    const Val &v = *cs[j].v;
    for (std::size_t i : rows) a.a.push_back(v.t == Val::ARR ? (i < v.a.size() ? v.a[i] : Val()) : v);
    put(o, cs[j].nm.c_str(), a);
  }
  return o;
}

bool drpempty(Ctx &c, Val &out) {
  std::vector<Dcol> cs;
  std::size_t n = 0;
  CK(dframe(c, c["data"], cs, n, "`data` must be a data frame"));
  CK(dfnm(c, cs));
  bool rows = true, colsw = true;
  if (c.has("which")) {
    const Val &w = c["which"];
    std::vector<std::string> ws;
    if (!labs(w, ws) || ws.empty()) return ckstop(c, "`which` must be one or more of \"rows\", \"cols\", not %s", jdep(w).c_str());
    rows = colsw = false;
    for (const std::string &s : ws) {
      const bool r = !s.empty() && std::string("rows").compare(0, s.size(), s) == 0;
      const bool k = !s.empty() && std::string("cols").compare(0, s.size(), s) == 0;
      if (!r && !k) return ckstop(c, "`which` must be one or more of \"rows\", \"cols\", not %s", jdep(w).c_str());
      rows = rows || r;
      colsw = colsw || k;
    }
  }
  bool quiet = false;
  CK(flagj(c, "quiet", false, quiet));
  std::vector<char> keepc(cs.size(), 1);
  if (colsw) {
    for (std::size_t j = 0; j < cs.size(); ++j) {
      bool all = true;
      for (char e : cs[j].na) all = all && e;
      if (all) keepc[j] = 0;
    }
  }
  std::vector<std::size_t> rr;
  bool anyc = false;
  for (char k : keepc) anyc = anyc || k;
  for (std::size_t i = 0; i < n; ++i) {
    bool all = true;
    for (std::size_t j = 0; j < cs.size(); ++j) {
      if (keepc[j]) all = all && cs[j].na[i];
    }
    if (!(rows && anyc && all)) rr.push_back(i);
  }
  out = subframe(cs, rr, keepc);
  return true;
}

bool drpconst(Ctx &c, Val &out) {
  std::vector<Dcol> cs;
  std::size_t n = 0;
  CK(dframe(c, c["data"], cs, n, "`data` must be a data frame"));
  CK(dfnm(c, cs));
  bool quiet = false, naas = false;
  CK(flagj(c, "quiet", false, quiet));
  if (c.has("naasval")) naas = c["naasval"].t == Val::BOOL && c["naasval"].b;
  std::vector<char> keepc(cs.size(), 1);
  for (std::size_t j = 0; j < cs.size(); ++j) {
    std::vector<std::string> u;
    for (std::size_t i = 0; i < n; ++i) {
      if (cs[j].na[i] && !naas) continue;
      const std::string k = cs[j].na[i] ? std::string("\x01NA") : cs[j].s[i];
      if (std::find(u.begin(), u.end(), k) == u.end()) u.push_back(k);
      if (u.size() > 1) break;
    }
    keepc[j] = u.size() > 1;
  }
  std::vector<std::size_t> rr;
  for (std::size_t i = 0; i < n; ++i) rr.push_back(i);
  out = subframe(cs, rr, keepc);
  return true;
}

// iconv(, to = "ASCII//TRANSLIT") for the Latin-1 letters; "?" for others.
inline std::string translit(const std::string &s) {
  static const char *const lat = "AAAAAAACEEEEIIIIDNOOOOOxOUUUUYTsaaaaaaaceeeeiiiidnooooo/ouuuuyty";
  std::string r;
  for (std::size_t i = 0; i < s.size(); ++i) {
    const unsigned char ch = static_cast<unsigned char>(s[i]);
    if (ch < 0x80) {
      r += static_cast<char>(ch);
      continue;
    }
    std::size_t len = ch >= 0xF0 ? 4 : ch >= 0xE0 ? 3 : 2;
    if (ch == 0xC3 && i + 1 < s.size()) {
      const int k = (static_cast<unsigned char>(s[i + 1]) & 0x3F) | 0x40;
      const char t = lat[k - 0x40];
      if (k == 0x46) r += "AE";
      else if (k == 0x66) r += "ae";
      else if (k == 0x5F) r += "ss";
      else r += t;
    } else {
      r += '?';
    }
    i += len - 1;
  }
  return r;
}

bool clncolnms(Ctx &c, Val &out) {
  static const std::vector<std::string> cases = {"snake", "lower_camel", "upper_camel", "screaming_snake", "none"};
  int cs = 0;
  CK(ckarg(c, "case", cases, cs));
  std::string sep = "_";
  if (c.has("sep")) {
    const Val &o = c["sep"];
    const Val *e = o.t == Val::ARR && o.a.size() == 1 ? &o.a[0] : &o;
    if (e->t != Val::STR) return ckstop(c, "`sep` must be a single string");
    sep = e->s;
  }
  bool quiet = false;
  CK(flagj(c, "quiet", false, quiet));
  const Val &d = c["data"];
  std::vector<std::string> nms;
  if (d.t == Val::OBJ) {
    nms = d.k;
  } else {
    if (!labs(d, nms) || d.t == Val::NUM || (d.t == Val::ARR && !d.a.empty() && d.a[0].t != Val::STR)) {
      return ckstop(c, "`data` must be a data frame or column names, without NA");
    }
  }
  if (nms.empty()) {
    out = d.t == Val::OBJ ? d : varr();
    return true;
  }
  auto low = [](std::string s) {
    for (char &ch : s) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    return s;
  };
  auto up = [](std::string s) {
    for (char &ch : s) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
    return s;
  };
  std::vector<std::string> res;
  for (const std::string &nm : nms) {
    std::string t = translit(nm), u;
    for (char ch : t) {
      if (ch == '%') u += " pct ";
      else if (ch == '&') u += " and ";
      else if (ch == '#') u += " n ";
      else u += ch;
    }
    std::string w;
    for (std::size_t i = 0; i < u.size(); ++i) {
      const unsigned char a = static_cast<unsigned char>(u[i]);
      w += u[i];
      if (i + 1 < u.size() && ((a >= 'a' && a <= 'z') || (a >= '0' && a <= '9')) && u[i + 1] >= 'A' &&
          u[i + 1] <= 'Z') {
        w += ' ';
      }
    }
    std::vector<std::string> words(1);
    for (char ch : w) {
      if (std::isalnum(static_cast<unsigned char>(ch)) && (ch & 0x80) == 0) words.back() += ch;
      else if (!words.back().empty()) words.emplace_back();
    }
    if (words.back().empty()) words.pop_back();
    if (words.empty()) words.push_back("x");
    std::string o;
    for (std::size_t i = 0; i < words.size(); ++i) {
      const std::string &x = words[i];
      std::string capw = up(x.substr(0, 1)) + low(x.substr(1));
      if (cs == 0) o += (i ? sep : "") + low(x);
      else if (cs == 3) o += (i ? sep : "") + up(x);
      else if (cs == 4) o += (i ? sep : "") + x;
      else if (cs == 1) o += i ? capw : low(x);
      else o += capw;
    }
    if (!o.empty() && std::isdigit(static_cast<unsigned char>(o[0]))) o = "x" + o;
    res.push_back(o);
  }
  std::vector<std::string> fin = res;
  for (std::size_t i = 0; i < res.size(); ++i) {
    int seen = 0;
    for (std::size_t k = 0; k < i; ++k) seen += res[k] == res[i];
    if (seen) fin[i] = res[i] + sep + std::to_string(seen + 1);
  }
  if (d.t == Val::OBJ) {
    out = d;
    out.k = fin;
  } else {
    out = unbox(vstrs(fin));
  }
  return true;
}

// The numeric columns of a frame as a column-major matrix (NA as NaN).
inline std::vector<double> nmat(const std::vector<const Dcol *> &num, std::size_t n) {
  std::vector<double> x(n * num.size());
  for (std::size_t j = 0; j < num.size(); ++j) {
    for (std::size_t i = 0; i < n; ++i) x[j * n + i] = num[j]->na[i] ? NaN : num[j]->x[i];
  }
  return x;
}

bool mahaoutl(Ctx &c, Val &out) {
  std::vector<Dcol> cs;
  std::size_t n = 0;
  CK(dframe(c, c["data"], cs, n, "`data` must be a data frame or a numeric matrix"));
  bool rob = true;
  CK(flagj(c, "robust", true, rob));
  std::vector<const Dcol *> num;
  for (const Dcol &k : cs) {
    if (k.ty == 0) num.push_back(&k);
  }
  if (num.empty()) return ckstop(c, "`data` has no numeric columns");
  double alpha = 0.001;
  if (c.has("alpha")) {
    alpha = val1(c["alpha"]);
    if (!(alpha > 0 && alpha < 1)) return ckstop(c, "`alpha` must be a single value strictly inside (0, 1)");
  }
  const int p = static_cast<int>(num.size());
  const std::vector<double> x = nmat(num, n);
  std::vector<int> ok;
  for (std::size_t i = 0; i < n; ++i) {
    bool all = true;
    for (int j = 0; j < p; ++j) all = all && !std::isnan(x[static_cast<std::size_t>(j) * n + i]);
    if (all) ok.push_back(static_cast<int>(i));
  }
  const int m = static_cast<int>(ok.size());
  if (m < p + 1) return ckstop(c, "need more complete rows than columns to estimate a covariance");
  std::vector<double> xc(static_cast<std::size_t>(m) * p), d2(static_cast<std::size_t>(m));
  for (int j = 0; j < p; ++j) {
    for (int i = 0; i < m; ++i) {
      xc[static_cast<std::size_t>(j) * m + i] = x[static_cast<std::size_t>(j) * n + static_cast<std::size_t>(ok[i])];
    }
  }
  if (clab_maha(xc.data(), m, p, rob, d2.data()) != CLAB_OK) {
    return ckstop(c,
                  "the columns are collinear, so no Mahalanobis distance is defined. Drop the redundant column(s) "
                  "-- drpconst() and topcors() will find them");
  }
  struct R {
    int row;
    double d, pv, lp;
  };
  std::vector<R> rs;
  for (int i = 0; i < m; ++i) {
    const double v = d2[static_cast<std::size_t>(i)];
    rs.push_back({ok[i] + 1, std::sqrt(v), clab_pchisq(v, p, 0), clab_lpchisq(v, p)});
  }
  std::stable_sort(rs.begin(), rs.end(), [](const R &a, const R &b) { return a.d > b.d; });
  Val row = varr(), dist = varr(), pv = varr(), lp = varr(), ol = varr();
  for (const R &r : rs) {
    row.a.push_back(vnum(r.row));
    dist.a.push_back(vnum(r.d));
    pv.a.push_back(vnum(r.pv));
    lp.a.push_back(vnum(r.lp));
    ol.a.push_back(vbool(r.pv < alpha));
  }
  std::size_t k = 0;
  for (std::size_t i = 0; i < n; ++i) {
    if (k < ok.size() && ok[k] == static_cast<int>(i)) {
      ++k;
      continue;
    }
    row.a.push_back(vnum(static_cast<double>(i) + 1));
    dist.a.push_back(Val());
    pv.a.push_back(Val());
    lp.a.push_back(Val());
    ol.a.push_back(Val());
  }
  out = vobj();
  put(out, "row", row);
  put(out, "distance", dist);
  put(out, "p_value", pv);
  put(out, "logpvalue", lp);
  put(out, "outlier", ol);
  return true;
}

// The non-empty frame of the missingness functions.
inline bool mframe(Ctx &c, std::vector<Dcol> &cs, std::size_t &n) {
  CK(dframe(c, c["data"], cs, n, "`data` must be a data frame"));
  CK(dfnm(c, cs));
  if (cs.empty()) return ckstop(c, "`data` has no columns");
  return true;
}

inline std::vector<int> namat(const std::vector<Dcol> &cs, std::size_t n) {
  std::vector<int> na(n * cs.size());
  for (std::size_t j = 0; j < cs.size(); ++j) {
    for (std::size_t i = 0; i < n; ++i) na[j * n + i] = cs[j].na[i];
  }
  return na;
}

bool missruns(Ctx &c, Val &out) {
  std::vector<Dcol> cs;
  std::size_t n = 0;
  CK(mframe(c, cs, n));
  int mr = 2;
  CK(ckposj(c, "minrun", 2, mr));
  struct R {
    std::string col;
    int st, en;
  };
  std::vector<R> rs;
  std::vector<int> st(n + 1), en(n + 1);
  for (const Dcol &k : cs) {
    std::vector<int> na(k.na.begin(), k.na.end());
    const int m = clab_nrun(na.data(), static_cast<int>(n), mr, st.data(), en.data());
    for (int i = 0; i < m; ++i) rs.push_back({k.nm, st[static_cast<std::size_t>(i)], en[static_cast<std::size_t>(i)]});
  }
  std::stable_sort(rs.begin(), rs.end(), [](const R &a, const R &b) {
    const int la = a.en - a.st, lb = b.en - b.st;
    return la != lb ? la > lb : a.col < b.col;
  });
  Val cl = varr(), s = varr(), e = varr(), l = varr();
  for (const R &r : rs) {
    cl.a.push_back(vstr(r.col));
    s.a.push_back(vnum(r.st));
    e.a.push_back(vnum(r.en));
    l.a.push_back(vnum(r.en - r.st + 1));
  }
  out = vobj();
  put(out, "column", cl);
  put(out, "start", s);
  put(out, "end", e);
  put(out, "length", l);
  return true;
}

bool missgrid(Ctx &c, Val &out) {
  std::vector<Dcol> cs;
  std::size_t n = 0;
  CK(mframe(c, cs, n));
  int h = 20, w = 12;
  CK(ckposj(c, "height", 20, h));
  if (c.has("width")) {
    const double v = val1(c["width"]);
    if (std::isfinite(v) && v != std::floor(v)) {
      return ckstop(c, "`width` must be a whole number, not %s", rfmt(v).c_str());
    }
    if (!(v >= 1 && v <= 10000)) return ckstop(c, "`width` must be one whole number from 1 to 10000");
    w = static_cast<int>(v);
  }
  const bool u = uniok();
  const std::string hr = u ? "─" : "-";
  auto rep = [](const std::string &s, int k) {
    std::string r;
    for (int i = 0; i < k; ++i) r += s;
    return r;
  };
  const std::string top = rep(hr, 2) + " Missingness map " + rep(hr, 66 - 19), bot = rep(hr, 66);
  std::vector<std::string> lines = {top};
  if (n == 0) {
    lines.push_back("  (no rows)");
    lines.push_back(bot);
    out = vstrs(lines);
    return true;
  }
  // names cut to w characters (UTF-8 aware), written down
  std::vector<std::vector<std::string>> ch;
  std::size_t mx = 0;
  for (const Dcol &k : cs) {
    std::vector<std::string> cc;
    for (std::size_t i = 0; i < k.nm.size() && static_cast<int>(cc.size()) < w;) {
      const unsigned char b = static_cast<unsigned char>(k.nm[i]);
      const std::size_t len = b >= 0xF0 ? 4 : b >= 0xE0 ? 3 : b >= 0xC0 ? 2 : 1;
      cc.push_back(k.nm.substr(i, len));
      i += len;
    }
    mx = std::max(mx, cc.size());
    ch.push_back(cc);
  }
  for (std::size_t r = 0; r < mx; ++r) {
    std::string l = "         ";
    for (const auto &cc : ch) l += r < cc.size() ? cc[r] : " ";
    lines.push_back(l);
  }
  static const char *const su[] = {" ", "░", "▒", "▓", "█"};
  static const char *const sa[] = {" ", ".", "-", "+", "#"};
  const std::vector<int> na = namat(cs, n);
  const int b = std::min(h, static_cast<int>(n)), p = static_cast<int>(cs.size());
  std::vector<int> st(static_cast<std::size_t>(b) + 1), sh(static_cast<std::size_t>(b + 1) * p);
  const int k = clab_mgrid(na.data(), static_cast<int>(n), p, b, st.data(), sh.data());
  for (int i = 0; i < k; ++i) {
    char lb[32];
    std::snprintf(lb, sizeof lb, "  %6d ", st[static_cast<std::size_t>(i)]);
    std::string l = lb;
    for (int j = 0; j < p; ++j) l += (u ? su : sa)[sh[static_cast<std::size_t>(j) * k + i]];
    lines.push_back(l);
  }
  lines.push_back("");
  lines.push_back(std::string("  legend: '") + (u ? su : sa)[0] + "' none  '" + (u ? su : sa)[4] + "' all missing");
  lines.push_back(bot);
  out = vstrs(lines);
  return true;
}

bool mcartst(Ctx &c, Val &out) {
  int maxit = 500;
  if (c.has("maxit")) {
    const double v = val1(c["maxit"]);
    if (std::isfinite(v) && v != std::floor(v)) {
      return ckstop(c, "`maxit` must be a whole number, not %s", rfmt(v).c_str());
    }
    if (!(v >= 1 && v <= 1e6)) return ckstop(c, "`maxit` must be one whole number from 1 to 1e+06");
    maxit = static_cast<int>(v);
  }
  double tol = 1e-7;
  if (c.has("tol")) {
    tol = val1(c["tol"]);
    if (!(tol >= 1e-12 && tol <= 1)) return ckstop(c, "`tol` must be one number from 1e-12 to 1");
  }
  std::vector<Dcol> cs;
  std::size_t n = 0;
  CK(dframe(c, c["data"], cs, n, "`data` must be a data frame or a numeric matrix"));
  std::vector<const Dcol *> num;
  std::string drop;
  for (const Dcol &k : cs) {
    if (k.ty == 0) num.push_back(&k);
    else drop += (drop.empty() ? "" : ", ") + k.nm;
  }
  if (!drop.empty()) c.warn.push_back(std::string(c.fn) + ": dropping non-numeric column(s): " + drop + ".");
  if (num.size() < 2) return ckstop(c, "Little's test needs at least two numeric columns");
  // drop wholly missing columns, then rows
  std::vector<const Dcol *> keep;
  for (const Dcol *k : num) {
    bool all = true;
    for (char e : k->na) all = all && e;
    if (!all) keep.push_back(k);
  }
  std::vector<std::size_t> rows;
  for (std::size_t i = 0; i < n; ++i) {
    bool all = true;
    for (const Dcol *k : keep) all = all && k->na[i];
    if (!all) rows.push_back(i);
  }
  const std::size_t m = rows.size(), p = keep.size();
  std::vector<double> x(m * p);
  bool inf = false;
  for (std::size_t j = 0; j < p; ++j) {
    for (std::size_t i = 0; i < m; ++i) {
      const double v = keep[j]->na[rows[i]] ? NaN : keep[j]->x[rows[i]];
      inf = inf || std::isinf(v);
      x[j * m + i] = v;
    }
  }
  if (inf) return ckstop(c, "`data` has infinite values, which Little's test cannot use; set them to NA first");
  if (p < 2 || m < 2) return ckstop(c, "too few complete columns or rows remain for Little's test");
  double o[8];
  std::vector<double> mu(p), sig(p * p);
  const int rc = clab_mcar(x.data(), static_cast<int>(m), static_cast<int>(p), maxit, tol, o, mu.data(), sig.data());
  if (rc != CLAB_OK) {
    if (o[0] == 1) {
      return ckstop(c,
                    "the columns are collinear (or one is constant), so the multivariate normal likelihood Little's "
                    "test is built on is degenerate. Drop the redundant column(s) first -- topcors() and drpconst() "
                    "will find them");
    }
    return ckstop(c, "the EM estimates did not converge to a usable covariance: the columns may be collinear, or a "
                     "pattern may have too few rows");
  }
  const double df = o[1];
  std::string note;
  if (df <= 0) {
    note = "df is " + rfmt(df) +
           ": with a single missingness pattern there is nothing to test. A complete table trivially satisfies MCAR.";
  }
  if (o[7] > 0) {
    note += (note.empty() ? "" : " ") + std::to_string(static_cast<int>(o[7])) +
            " pattern(s) were skipped: their observed variables were collinear.";
  }
  out = vobj();
  put(out, "statistic", vnum(o[0]));
  put(out, "df", vnum(df));
  put(out, "logpvalue", df > 0 ? vnum(o[3]) : Val());
  put(out, "p_value", df > 0 ? vnum(o[2]) : Val());
  put(out, "npat", vnum(o[4]));
  put(out, "nused", vnum(o[5]));
  put(out, "nvars", vnum(static_cast<double>(p)));
  put(out, "iterations", vnum(o[6]));
  Val mv = vobj();
  for (std::size_t j = 0; j < p; ++j) put(mv, keep[j]->nm.c_str(), vnum(mu[j]));
  put(out, "mu", mv);
  Val sm = varr();
  for (std::size_t i = 0; i < p; ++i) {
    Val r = varr();
    for (std::size_t j = 0; j < p; ++j) r.a.push_back(vnum(sig[j * p + i]));
    sm.a.push_back(r);
  }
  put(out, "sigma", sm);
  put(out, "note", note.empty() ? Val() : vstr(note));
  return true;
}

bool misspat(Ctx &c, Val &out) {
  std::vector<Dcol> cs;
  std::size_t n = 0;
  CK(mframe(c, cs, n));
  int mp = 20;
  CK(ckposj(c, "maxpat", 20, mp));
  std::vector<std::string> code(n);
  std::vector<char> na(n, 0);
  for (std::size_t i = 0; i < n; ++i) {
    for (const Dcol &k : cs) code[i] += k.na[i] ? 'X' : '.';
  }
  std::vector<std::string> lv;
  std::vector<int> lev, cnt;
  if (n > 0) freq(code, na, true, mp, false, lv, lev, cnt);
  Val pt = varr(), nr = varr(), pc = varr(), nm = varr(), cl = varr();
  for (std::size_t i = 0; i < lev.size(); ++i) {
    const std::string &s = lv[static_cast<std::size_t>(lev[i] - 1)];
    pt.a.push_back(vstr(s));
    nr.a.push_back(vnum(cnt[i]));
    pc.a.push_back(vnum(100 * static_cast<double>(cnt[i]) / static_cast<double>(n)));
    int k = 0;
    std::string names;
    for (std::size_t j = 0; j < s.size(); ++j) {
      if (s[j] != 'X') continue;
      ++k;
      names += (names.empty() ? "" : ", ") + cs[j].nm;
    }
    nm.a.push_back(vnum(k));
    cl.a.push_back(vstr(names));
  }
  out = vobj();
  put(out, "pattern", pt);
  put(out, "n_rows", nr);
  put(out, "pctrows", pc);
  put(out, "n_missing", nm);
  put(out, "columns", cl);
  return true;
}

bool duprows(Ctx &c, Val &out) {
  std::vector<Dcol> cs;
  std::size_t n = 0;
  CK(dframe(c, c["data"], cs, n, "`data` must be a data frame"));
  bool quiet = false;
  CK(flagj(c, "quiet", false, quiet));
  std::vector<std::size_t> pos;
  if (n > 0) {
    if (c.has("cols")) {
      std::vector<std::string> cc;
      labs(c["cols"], cc);
      std::string miss;
      for (const std::string &s : cc) {
        bool f = false;
        for (const Dcol &k : cs) f = f || k.nm == s;
        if (!f && miss.find(s) == std::string::npos) miss += (miss.empty() ? "" : ", ") + s;
      }
      if (!miss.empty()) return ckstop(c, "no such column(s): %s", miss.c_str());
      for (const std::string &s : cc) {
        int seen = 0;
        std::size_t at = 0;
        for (std::size_t j = 0; j < cs.size(); ++j) {
          if (cs[j].nm == s) {
            if (!seen) at = j;
            ++seen;
          }
        }
        if (seen > 1) {
          return ckstop(c, "`data` has two columns named %s: give the columns by position, or rename one", s.c_str());
        }
        pos.push_back(at);
      }
    } else {
      for (std::size_t j = 0; j < cs.size(); ++j) pos.push_back(j);
    }
  }
  std::vector<std::string> key(n);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t k = 0; k < pos.size(); ++k) key[i] += (k ? "\r" : "") + cs[pos[k]].s[i];
  }
  std::vector<std::string> lv = key;
  std::sort(lv.begin(), lv.end());
  lv.erase(std::unique(lv.begin(), lv.end()), lv.end());
  std::vector<int> code(n), rows(n), ct(n);
  for (std::size_t i = 0; i < n; ++i) {
    code[i] = static_cast<int>(std::lower_bound(lv.begin(), lv.end(), key[i]) - lv.begin()) + 1;
  }
  const int k = n ? clab_dups(code.data(), static_cast<int>(n), static_cast<int>(lv.size()), rows.data(), ct.data())
                  : 0;
  std::vector<std::size_t> rr;
  for (int i = 0; i < k; ++i) rr.push_back(static_cast<std::size_t>(rows[static_cast<std::size_t>(i)] - 1));
  std::vector<char> all(cs.size(), 1);
  out = subframe(cs, rr, all);
  Val dc = varr();
  for (int i = 0; i < k; ++i) dc.a.push_back(vnum(ct[static_cast<std::size_t>(i)]));
  put(out, "dupe_count", dc);
  return true;
}

bool misssmry(Ctx &c, Val &out) {
  std::vector<Dcol> cs;
  std::size_t n = 0;
  CK(dframe(c, c["data"], cs, n, "`data` must be a data frame"));
  if (cs.empty()) return ckstop(c, "`data` has no columns");
  const std::vector<int> na = namat(cs, n);
  double o[8];
  clab_msum(na.data(), static_cast<int>(n), static_cast<int>(cs.size()), o);
  static const char *const nm[] = {"n_rows",          "n_cols",           "n_missing",          "pct_missing",
                                   "n_complete_rows", "pct_complete_rows", "n_cols_any_missing", "n_cols_all_missing"};
  out = nvec(8, nm, o);
  return true;
}

// format(x, trim = TRUE) of a numeric vector: 7 significant digits, the
// same decimals for all.
inline std::vector<std::string> rformat(const std::vector<double> &x) {
  int dec = 0;
  for (double v : x) {
    for (int d = 0; d < 15; ++d) {
      char b[64];
      std::snprintf(b, sizeof b, "%.*f", d, v);
      char s7[64];
      std::snprintf(s7, sizeof s7, "%.7g", v);
      if (std::strtod(b, nullptr) == std::strtod(s7, nullptr)) {
        dec = std::max(dec, d);
        break;
      }
    }
  }
  std::vector<std::string> r;
  for (double v : x) {
    char b[64];
    std::snprintf(b, sizeof b, "%.*f", dec, v);
    r.push_back(b);
  }
  return r;
}

bool profcols(Ctx &c, Val &out) {
  std::vector<Dcol> cs;
  std::size_t n = 0;
  CK(dframe(c, c["data"], cs, n, "`data` must be a data frame"));
  CK(dfnm(c, cs));
  if (cs.empty()) return ckstop(c, "`data` has no columns to profile");
  std::vector<double> q = {0.25, 0.5, 0.75};
  if (c.has("qnts")) {
    Dv d;
    if (readd(c["qnts"], d) != NUM) return ckstop(c, "`qnts` must be probabilities in [0, 1]");
    q = d.v;
  } else if (c.in != nullptr) {
    for (std::size_t i = 0; i < c.in->k.size(); ++i) {
      if (c.in->k[i] == "qnts") q.clear();  // qnts = NULL: none
    }
  }
  for (double v : q) {
    if (!(v >= 0 && v <= 1)) return ckstop(c, "`qnts` must be probabilities in [0, 1]");
  }
  static const char *const nm[] = {"column", "type", "n", "n_missing", "pct_missing", "n_distinct", "nzero",
                                   "nneg", "ninf", "mean", "sd", "median", "mad", "min", "max", "skewness",
                                   "n_outliers", "hist", "top"};
  std::vector<Val> col(19, varr());
  std::vector<Val> qc(q.size(), varr());
  for (const Dcol &k : cs) {
    int nmiss = 0;
    for (char e : k.na) nmiss += e;
    col[0].a.push_back(vstr(k.nm));
    col[1].a.push_back(vstr(k.ty == 0 ? "numeric" : k.ty == 1 ? "character" : "logical"));
    col[2].a.push_back(vnum(static_cast<double>(n)));
    col[3].a.push_back(vnum(nmiss));
    col[4].a.push_back(vnum(100 * static_cast<double>(nmiss) / std::max<double>(1, static_cast<double>(n))));
    std::vector<Val> st(14, Val());
    st[0] = vnum(0);
    std::vector<double> fin;
    if (k.ty != 0) {
      std::vector<std::string> lv;
      std::vector<int> lev, cnt;
      freq(k.s, k.na, true, INT_MAX, false, lv, lev, cnt);
      st[0] = vnum(static_cast<double>(lv.size()));
      if (!lev.empty()) st[13] = vstr(lv[static_cast<std::size_t>(lev[0] - 1)]);
    } else {
      std::vector<double> ok;
      for (std::size_t i = 0; i < n; ++i) {
        if (!k.na[i]) ok.push_back(k.x[i]);
      }
      if (!ok.empty()) {
        std::vector<double> u = ok;
        std::sort(u.begin(), u.end());
        u.erase(std::unique(u.begin(), u.end()), u.end());
        int z = 0, neg = 0, inf = 0;
        for (double v : ok) {
          z += v == 0;
          neg += v < 0;
          inf += std::isinf(v);
          if (std::isfinite(v)) fin.push_back(v);
        }
        st[0] = vnum(static_cast<double>(u.size()));
        st[1] = vnum(z);
        st[2] = vnum(neg);
        st[3] = vnum(inf);
        if (!fin.empty()) {
          double s[7];
          clab_pstat(fin.data(), static_cast<int>(fin.size()), s);
          st[4] = vnum(s[0]);
          st[5] = vnum(s[1]);
          st[6] = vnum(s[3]);
          st[7] = vnum(s[4]);
          st[8] = vnum(*std::min_element(fin.begin(), fin.end()));
          st[9] = vnum(*std::max_element(fin.begin(), fin.end()));
          st[10] = vnum(s[2]);
          int ol = 0;
          for (double v : fin) ol += v < s[5] || v > s[6];
          st[11] = vnum(ol);
          st[12] = vstr(spark(fin, 10));
        }
      }
    }
    for (int j = 0; j < 14; ++j) col[static_cast<std::size_t>(5 + j)].a.push_back(st[static_cast<std::size_t>(j)]);
    std::vector<double> qv(q.size(), NA);
    if (!fin.empty()) clab_q7(fin.data(), static_cast<int>(fin.size()), q.data(), static_cast<int>(q.size()), qv.data());
    for (std::size_t j = 0; j < q.size(); ++j) qc[j].a.push_back(vnum(qv[j]));
  }
  out = vobj();
  for (int j = 0; j < 19; ++j) put(out, nm[j], col[static_cast<std::size_t>(j)]);
  std::vector<double> pq;
  for (double v : q) pq.push_back(100 * v);
  const std::vector<std::string> qn = rformat(pq);
  for (std::size_t j = 0; j < q.size(); ++j) {
    const std::string k = "q" + qn[j];
    put(out, k.c_str(), qc[j]);
  }
  return true;
}

bool dstntst(Ctx &c, Val &out) {
  std::vector<Dcol> cs;
  std::size_t n = 0;
  CK(dframe(c, c["data"], cs, n, "`data` must be a data frame"));
  auto colj = [&](const char *what, const Dcol *&k) {
    const Val &o = c[what];
    const Val *e = o.t == Val::ARR && o.a.size() == 1 ? &o.a[0] : &o;
    if (e->t != Val::STR || e->s.find_first_not_of(" \t\n") == std::string::npos) {
      return ckstop(c, "`%s` must be a single non-empty string", what);
    }
    k = nullptr;
    for (const Dcol &d : cs) {
      if (d.nm == e->s && k == nullptr) k = &d;
    }
    if (k == nullptr) return ckstop(c, "column `%s` not found in `data`", e->s.c_str());
    return true;
  };
  const Dcol *g = nullptr, *es = nullptr, *lo = nullptr, *hi = nullptr, *ql = nullptr;
  CK(colj("grp", g));
  CK(colj("est", es));
  CK(colj("lower", lo));
  CK(colj("upper", hi));
  if (c.has("qual")) CK(colj("qual", ql));
  double lvl = 0.95;
  if (c.has("conflvl")) {
    lvl = val1(c["conflvl"]);
    if (!(lvl > 0 && lvl < 1)) return ckstop(c, "`conflvl` must be a single number strictly between 0 and 1");
  }
  const double z = clab_qnorm(1 - (1 - lvl) / 2);
  auto trim = [](const std::string &s) {
    const std::size_t a = s.find_first_not_of(" \t\n\r\f\v");
    if (a == std::string::npos) return std::string();
    return s.substr(a, s.find_last_not_of(" \t\n\r\f\v") - a + 1);
  };
  std::vector<std::string> gl;
  for (std::size_t i = 0; i < n; ++i) gl.push_back(g->na[i] ? "NA" : trim(g->s[i]));
  std::vector<std::string> pooled;
  const bool hasp = c.has("pooled");
  if (hasp) {
    if (!labs(c["pooled"], pooled)) return ckstop(c, "`pooled` must be the labels of the rows that pool other groups");
    std::string miss;
    for (const std::string &p : pooled) {
      if (std::find(gl.begin(), gl.end(), trim(p)) == gl.end()) miss += (miss.empty() ? "\"" : ", \"") + trim(p) + "\"";
    }
    if (!miss.empty()) return ckstop(c, "`pooled` names %s, not a group of `data`", miss.c_str());
  }
  for (const Dcol *k : {es, lo, hi}) {
    if (k->ty == 1) {
      for (std::size_t i = 0; i < n; ++i) {
        if (!k->na[i]) {
          return ckstop(c, "`%s` must be numeric, not text such as \"%s\"", k == es ? "est" : k == lo ? "lower" : "upper",
                        k->s[i].c_str());
        }
      }
    }
  }
  std::vector<double> y(n), s(n);
  std::string bad;
  for (std::size_t i = 0; i < n; ++i) {
    y[i] = es->na[i] ? NaN : es->x[i];
    const double l = lo->na[i] ? NaN : lo->x[i], h = hi->na[i] ? NaN : hi->x[i];
    if (std::isfinite(l) && std::isfinite(h) && std::isfinite(y[i]) && (l > y[i] || h < y[i])) {
      if (bad.find(gl[i]) == std::string::npos) bad += (bad.empty() ? "" : ", ") + gl[i];
      y[i] = NaN;
    }
    s[i] = (h - l) / (2 * z);
  }
  if (!bad.empty()) {
    c.warn.push_back(std::string(c.fn) + ": the interval does not contain the estimate for: " + bad +
                     "; not compared.");
  }
  std::vector<std::string> q(n, "ok");
  if (ql != nullptr) {
    for (std::size_t i = 0; i < n; ++i) q[i] = ql->na[i] ? "ok" : trim(ql->s[i]);
  }
  std::string odd;
  for (std::string &f : q) {
    if (f.empty()) f = "ok";
    if (f != "ok" && f != "E" && f != "F" && f != "x" && f != ".." && f != "...") {
      if (odd.find("\"" + f + "\"") == std::string::npos) odd += (odd.empty() ? "\"" : ", \"") + f + "\"";
      f = "ok";
    }
  }
  if (!odd.empty()) c.warn.push_back(std::string(c.fn) + ": quality flag(s) " + odd + " not recognised: treated as ok.");
  std::vector<double> yu, su;
  Val gu = varr(), left = varr();
  std::vector<std::string> used;
  for (std::size_t i = 0; i < n; ++i) {
    bool pan;
    if (hasp) {
      pan = false;
      for (const std::string &p : pooled) pan = pan || trim(p) == gl[i];
    } else {
      std::string t = gl[i];
      for (char &ch : t) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
      pan = t.rfind("total", 0) == 0 || t.rfind("indigenous identity", 0) == 0 ||
            t.rfind("aboriginal identity", 0) == 0 || t.rfind("indigenous population", 0) == 0;
    }
    const bool sup = q[i] == "F" || q[i] == "x" || q[i] == ".." || q[i] == "...";
    if (!pan && !sup && std::isfinite(y[i]) && std::isfinite(s[i]) && s[i] > 0) {
      yu.push_back(y[i]);
      su.push_back(s[i]);
      used.push_back(gl[i]);
      gu.a.push_back(vstr(gl[i]));
    } else {
      left.a.push_back(vstr(gl[i]));
    }
  }
  if (yu.size() < 2) return ckstop(c, "the test needs two groups at least with an estimate and an interval");
  for (std::size_t i = 0; i < used.size(); ++i) {
    for (std::size_t k = 0; k < i; ++k) {
      if (used[i] == used[k]) return ckstop(c, "a group appears twice");
    }
  }
  double r[6];
  clab_cochq(yu.data(), su.data(), static_cast<int>(yu.size()), r);
  Val grp = vobj();
  put(grp, "group", gu);
  put(grp, "estimate", vnums(yu));
  put(grp, "se", vnums(su));
  out = vobj();
  put(out, "groups", grp);
  put(out, "pooled", vnum(r[0]));
  put(out, "poose", vnum(r[1]));
  put(out, "q", vnum(r[2]));
  put(out, "df", vnum(static_cast<double>(yu.size()) - 1));
  put(out, "p_value", vnum(r[4]));
  put(out, "i2", vnum(r[5]));
  put(out, "leftout", unbox(left));
  put(out, "conflvl", vnum(lvl));
  return true;
}
