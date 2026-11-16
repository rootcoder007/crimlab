// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Data description (src/dsc.cpp, src/dscd.cpp), with the checks of R/dsc.R in their
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

// .spark(): the sketch of finite values x (clab_spark()).
inline std::string spark(const std::vector<double> &x, int bins) {
  std::vector<char> b(static_cast<std::size_t>(bins) * 3 + 1, 0);
  std::int64_t len = 0;
  clab_spark(x.data(), static_cast<std::int64_t>(x.size()), bins, uniok(), b.data(), &len);
  return std::string(b.data(), static_cast<std::size_t>(len));
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
  std::vector<std::pair<std::size_t, std::size_t>> pr;
  for (std::size_t i = 0; i < num.size(); ++i) {
    for (std::size_t j = i + 1; j < num.size(); ++j) pr.push_back({i, j});
  }
  std::vector<int> idx(r.size() + 1, 0);
  std::int64_t m = 0;
  clab_topc(r.data(), static_cast<std::int64_t>(r.size()), mn, top, idx.data(), &m);
  Val x = varr(), y = varr(), cr = varr(), ab = varr();
  for (std::int64_t t = 0; t < m; ++t) {
    const std::size_t k = static_cast<std::size_t>(idx[static_cast<std::size_t>(t)] - 1);
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
  const std::vector<std::string> fin = dscd::cnames(nms, cs, sep);
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
  std::vector<int> rw(n, 0), ob(n, 0);
  std::vector<double> d(n, 0.0), pv(n, 0.0), lp(n, 0.0);
  int why = 0;
  clab_mhtab(x.data(), static_cast<std::int64_t>(n), p, rob, alpha, rw.data(), d.data(), pv.data(), lp.data(),
             ob.data(), &why);
  if (why == 1) return ckstop(c, "need more complete rows than columns to estimate a covariance");
  if (why == 2) {
    return ckstop(c,
                  "the columns are collinear, so no Mahalanobis distance is defined. Drop the redundant column(s) "
                  "-- drpconst() and topcors() will find them");
  }
  Val row = varr(), dist = varr(), pval = varr(), lpv = varr(), ol = varr();
  for (std::size_t i = 0; i < n; ++i) {
    const bool na = ob[i] < 0;
    row.a.push_back(vnum(rw[i]));
    dist.a.push_back(na ? Val() : vnum(d[i]));
    pval.a.push_back(na ? Val() : vnum(pv[i]));
    lpv.a.push_back(na ? Val() : vnum(lp[i]));
    ol.a.push_back(na ? Val() : vbool(ob[i] == 1));
  }
  out = vobj();
  put(out, "row", row);
  put(out, "distance", dist);
  put(out, "p_value", pval);
  put(out, "logpvalue", lpv);
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
  // ties in byte order of the names
  std::vector<std::size_t> o(cs.size(), 0);
  for (std::size_t j = 0; j < cs.size(); ++j) o[j] = j;
  std::sort(o.begin(), o.end(), [&](std::size_t a, std::size_t b) { return cs[a].nm < cs[b].nm; });
  std::vector<int> rank(cs.size() + 1, 0);
  for (std::size_t j = 0; j < o.size(); ++j) rank[o[j]] = static_cast<int>(j) + 1;
  const std::vector<int> na = namat(cs, n);
  const std::size_t cap = cs.size() * (n / 2 + 1);
  std::vector<int> col(cap, 0), st(cap, 0), en(cap, 0);
  std::int64_t k = 0;
  clab_mruns(na.data(), static_cast<std::int64_t>(n), static_cast<int>(cs.size()), mr, rank.data(), col.data(),
             st.data(), en.data(), &k);
  Val cl = varr(), s = varr(), e = varr(), l = varr();
  for (std::size_t i = 0; i < static_cast<std::size_t>(k); ++i) {
    cl.a.push_back(vstr(cs[static_cast<std::size_t>(col[i] - 1)].nm));
    s.a.push_back(vnum(st[i]));
    e.a.push_back(vnum(en[i]));
    l.a.push_back(vnum(en[i] - st[i] + 1));
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
  std::vector<std::string> nm;
  for (const Dcol &k : cs) nm.push_back(k.nm);
  const std::vector<int> na = namat(cs, n);
  out = vstrs(dscd::mgrid(na.data(), n, cs.size(), nm, h, w, uniok()));
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
  const std::vector<double> xa = nmat(num, n);
  std::vector<double> x(xa.size() + 1, 0.0);
  std::vector<int> kc(num.size() + 1, 0);
  std::int64_t mm = 0;
  int kk = 0, why = 0;
  clab_mcdrop(xa.data(), static_cast<std::int64_t>(n), static_cast<int>(num.size()), x.data(), &mm, kc.data(), &kk,
              &why);
  if (why == 1) return ckstop(c, "`data` has infinite values, which Little's test cannot use; set them to NA first");
  if (why == 2) return ckstop(c, "too few complete columns or rows remain for Little's test");
  const std::size_t m = static_cast<std::size_t>(mm), p = static_cast<std::size_t>(kk);
  std::vector<const Dcol *> keep;
  for (std::size_t j = 0; j < p; ++j) keep.push_back(num[static_cast<std::size_t>(kc[j] - 1)]);
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
  const std::string note = dscd::mcnote(o);
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
  std::vector<std::string> nm;
  for (const Dcol &k : cs) nm.push_back(k.nm);
  const std::vector<int> na = namat(cs, n);
  const dscd::pat r = dscd::mpat(na.data(), n, cs.size(), nm, mp);
  Val pt = varr(), nr = varr(), pc = varr(), nmv = varr(), cl = varr();
  for (std::size_t i = 0; i < r.code.size(); ++i) {
    pt.a.push_back(vstr(r.code[i]));
    nr.a.push_back(vnum(r.rows[i]));
    pc.a.push_back(vnum(r.pct[i]));
    nmv.a.push_back(vnum(r.nmiss[i]));
    cl.a.push_back(vstr(r.cols[i]));
  }
  out = vobj();
  put(out, "pattern", pt);
  put(out, "n_rows", nr);
  put(out, "pctrows", pc);
  put(out, "n_missing", nmv);
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
    std::vector<double> qv(q.size() + 1, NA);
    if (k.ty != 0) {
      std::vector<std::string> lv;
      std::vector<int> lev, cnt;
      freq(k.s, k.na, true, INT_MAX, false, lv, lev, cnt);
      st[0] = vnum(static_cast<double>(lv.size()));
      if (!lev.empty()) st[13] = vstr(lv[static_cast<std::size_t>(lev[0] - 1)]);
    } else {
      std::vector<double> x(n, 0.0);
      for (std::size_t i = 0; i < n; ++i) x[i] = k.na[i] ? NaN : k.x[i];
      double v[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
      clab_prof(x.data(), static_cast<std::int64_t>(n), q.data(), static_cast<int>(q.size()), v, qv.data());
      for (int j = 0; j < 12; ++j) {
        if (!std::isnan(v[j]) || j == 0) st[static_cast<std::size_t>(j)] = vnum(v[j]);
      }
      if (!std::isnan(v[8])) st[12] = vstr(spark(x, 10));
    }
    for (int j = 0; j < 14; ++j) col[static_cast<std::size_t>(5 + j)].a.push_back(st[static_cast<std::size_t>(j)]);
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
  std::vector<double> y(n), l(n), h(n);
  std::vector<std::string> q, pl;
  std::vector<char> gna(n, 0), qna;
  for (std::size_t i = 0; i < n; ++i) {
    y[i] = es->na[i] ? NaN : es->x[i];
    l[i] = lo->na[i] ? NaN : lo->x[i];
    h[i] = hi->na[i] ? NaN : hi->x[i];
    gna[i] = g->na[i];
    if (ql != nullptr) {
      q.push_back(ql->na[i] ? std::string() : trim(ql->s[i]));
      qna.push_back(ql->na[i]);
    }
  }
  for (const std::string &p : pooled) pl.push_back(trim(p));
  const dscd::sel sl = dscd::dtsel(gl, gna, y.data(), l.data(), h.data(), q, qna, hasp, pl, z);
  auto join = [](const std::vector<std::string> &v, const char *qt) {
    std::string o;
    for (const std::string &e : v) o += (o.empty() ? "" : ", ") + std::string(qt) + e + qt;
    return o;
  };
  std::vector<std::string> bad = sl.bad;
  for (std::size_t i = 0; i < bad.size(); ++i) {
    if (sl.bna[i]) bad[i] = "NA";
  }
  if (!bad.empty()) {
    c.warn.push_back(std::string(c.fn) + ": the interval does not contain the estimate for: " + join(bad, "") +
                     "; not compared.");
  }
  if (!sl.odd.empty()) {
    c.warn.push_back(std::string(c.fn) + ": quality flag(s) " + join(sl.odd, "\"") + " not recognised: treated as ok.");
  }
  std::vector<double> yu, su;
  Val gu = varr(), left = varr();
  std::vector<std::string> used;
  for (std::size_t i = 0; i < n; ++i) {
    if (sl.use[i]) {
      yu.push_back(sl.y[i]);
      su.push_back(sl.s[i]);
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
