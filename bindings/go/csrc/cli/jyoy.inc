// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The change-table labels, summary and text writers (R/yoyw.R) for
// clab_jrun(): a table is a data frame (an object of columns), numbers
// are written by src/ynum.h as in R. With `file` the text is written
// there and the path returned; without, the text is returned. The HTML
// and PDF writers are R's alone (yoyhtml(), yoypdf()).

struct Yc {
  int k = 0;  // 0 numbers, 1 text, 2 true and false
  std::vector<double> d;
  std::vector<std::string> s;
  std::vector<char> na;
  std::size_t n() const { return na.size(); }
};

bool ycol(const Val &v, Yc &c) {
  c = Yc();
  const std::vector<Val> one = {v};
  const std::vector<Val> &es = v.t == Val::ARR ? v.a : one;
  bool num = true, lgl = true;
  for (const Val &e : es) {
    double x = 0;
    num = num && (e.t == Val::NUM || e.t == Val::NUL || snum(e, x));
    lgl = lgl && (e.t == Val::BOOL || e.t == Val::NUL);
  }
  c.k = lgl && !es.empty() && !num ? 2 : num ? 0 : 1;
  if (lgl && num) c.k = 2;
  for (const Val &e : es) {
    std::string t;
    double x = NA;
    const bool miss = e.t == Val::NUL;
    if (e.t == Val::NUM) x = e.n;
    else snum(e, x);
    if (!miss && !lab(e, t)) return false;
    c.s.push_back(t);
    c.d.push_back(e.t == Val::BOOL ? (e.b ? 1 : 0) : x);
    c.na.push_back(miss);
  }
  return true;
}

struct Yt {
  std::vector<std::string> nm;
  std::vector<Yc> c;
  std::size_t n = 0;
  const Yc *get(const char *k) const {
    for (std::size_t j = 0; j < nm.size(); ++j) {
      if (nm[j] == k) return &c[j];
    }
    return nullptr;
  }
};

bool ytab(Ctx &c, Yt &t, bool change) {
  const Val &x = c["x"];
  bool ok = x.t == Val::OBJ;
  t = Yt();
  for (std::size_t j = 0; ok && j < x.a.size(); ++j) {
    Yc v;
    ok = ycol(x.a[j], v) && (j == 0 || v.n() == t.n);
    t.n = v.n();
    t.nm.push_back(x.k[j]);
    t.c.push_back(v);
  }
  if (change) {
    for (const char *k : {"prd", "value", "prev", "change"}) ok = ok && t.get(k) && t.get(k)->k == 0;
    ok = ok && (t.get("pctchg") || t.get("ppchg"));
    if (!ok) {
      return c.fail("`x` must be a change table: a data frame with the columns prd, value, prev, change and "
                    "pctchg or ppchg, as yoypx() gives.");
    }
  }
  return ok || c.fail("`x` must be a data frame.");
}

bool yfile(Ctx &c, std::string &f) {
  f.clear();
  if (!c.has("file")) return true;
  if (c["file"].t != Val::STR || c["file"].s.empty()) return c.fail("`file` must be NULL or one file path.");
  f = c["file"].s;
  return true;
}

bool ydone(Ctx &c, const std::string &txt, const std::string &f, Val &out) {
  if (f.empty()) {
    out = vstr(txt);
    return true;
  }
  std::FILE *h = std::fopen(f.c_str(), "wb");
  if (h == nullptr) return c.fail("cannot open file '%s'.", f.c_str());
  const bool ok = std::fwrite(txt.data(), 1, txt.size(), h) == txt.size();
  if (std::fclose(h) != 0 || !ok) return c.fail("cannot write file '%s'.", f.c_str());
  out = vstr(f);
  return true;
}

bool ytxt1(Ctx &c, const char *nm, const char *def, std::string &s) {
  if (!c.has(nm)) {
    s = def;
    return true;
  }
  if (c[nm].t != Val::STR) return c.fail("`%s` must be one string.", nm);
  s = c[nm].s;
  return true;
}

bool yflag(Ctx &c, const char *nm, bool def, bool &b) {
  b = def;
  if (!c.has(nm)) return true;
  if (c[nm].t != Val::BOOL) return c.fail("`%s` must be TRUE or FALSE.", nm);
  b = c[nm].b;
  return true;
}

bool ydig(Ctx &c, const char *nm, int def, int &d, bool null) {
  d = def;
  if (!c.has(nm)) return true;
  const double v = val1(c[nm]);
  if (std::isnan(v) || v < 0 || v > 15 || v != std::floor(v)) {
    return c.fail("`%s` must be one whole number from 0 to 15.", nm);
  }
  (void)null;
  d = static_cast<int>(v);
  return true;
}

bool ymeta(Ctx &c, std::string &m) {
  m.clear();
  if (!c.has("meta")) return true;
  std::vector<std::string> v;
  const Val &o = c["meta"];
  bool ok = o.t == Val::STR || o.t == Val::ARR;
  for (const Val &e : o.a) ok = ok && e.t == Val::STR;
  if (!ok || !labs(o, v)) return c.fail("`meta` must be NULL or text.");
  for (std::size_t i = 0; i < v.size(); ++i) m += (i ? "; " : "") + v[i];
  std::string r;
  for (std::size_t i = 0; i < m.size(); ++i) {
    if (m[i] == '\r' || m[i] == '\n') {
      if (r.empty() || r.back() != ' ' || (i > 0 && m[i - 1] != '\r' && m[i - 1] != '\n')) r += ' ';
    } else {
      r += m[i];
    }
  }
  m = r;
  return true;
}

// round(x, d) as R rounds: to the nearest of the two candidates
double rround(double x, int d) {
  if (!std::isfinite(x)) return x;
  const double p = std::pow(10.0, d);
  return std::nearbyint(x * p) / p;
}

void ypctround(Yt &t, int d) {
  for (const char *k : {"pctchg", "ppchg", "pctlo", "pctup"}) {
    for (std::size_t j = 0; j < t.nm.size(); ++j) {
      if (t.nm[j] == k) {
        for (double &v : t.c[j].d) v = rround(v, d);
      }
    }
  }
}

std::string ycell(const Yc &v, std::size_t i) {
  if (v.k == 0) return ynum::all(v.d[i]);
  return v.s[i];
}

bool ydelim(Ctx &c, Val &out, char sep) {
  std::string f, na, com, m;
  CK(yfile(c, f));
  CK(ytxt1(c, "na", "", na));
  CK(ytxt1(c, "comment", "#", com));
  Yt t;
  CK(ytab(c, t, false));
  if (c.has("digits")) {
    int d = 0;
    CK(ydig(c, "digits", 0, d, true));
    ypctround(t, d);
  }
  CK(ymeta(c, m));
  auto cell = [sep, &na](const std::string &s, bool miss, bool txt) {
    if (miss) return na;
    std::string o = s;
    if (sep != ',') {
      for (char &ch : o) {
        if (ch == '\t' || ch == '\r' || ch == '\n') ch = ' ';
      }
      return o;
    }
    if (o.find_first_of(",\"\r\n") != std::string::npos || (txt && o.empty())) {
      std::string q = "\"";
      for (char ch : o) q += ch == '"' ? std::string("\"\"") : std::string(1, ch);
      return q + "\"";
    }
    return o;
  };
  std::string txt;
  if (!m.empty()) txt += com + " " + m + "\n";
  for (std::size_t j = 0; j < t.nm.size(); ++j) txt += (j ? std::string(1, sep) : "") + cell(t.nm[j], false, true);
  txt += "\n";
  for (std::size_t i = 0; i < t.n; ++i) {
    for (std::size_t j = 0; j < t.nm.size(); ++j) {
      const Yc &v = t.c[j];
      txt += (j ? std::string(1, sep) : "") + cell(ycell(v, i), v.na[i] != 0, v.k == 1);
    }
    txt += "\n";
  }
  return ydone(c, txt, f, out);
}
bool yoycsv(Ctx &c, Val &out) { return ydelim(c, out, ','); }
bool yoytsv(Ctx &c, Val &out) { return ydelim(c, out, '\t'); }

std::string yjs1s(const std::string &s) {
  std::string o = "\"";
  for (unsigned char ch : s) {
    if (ch == '\\') o += "\\\\";
    else if (ch == '"') o += "\\\"";
    else if (ch == '\n') o += "\\n";
    else if (ch == '\r') o += "\\r";
    else if (ch == '\t') o += "\\t";
    else if (ch < 32) {
      char b[8];
      std::snprintf(b, sizeof b, "\\u%04x", ch);
      o += b;
    } else {
      o += static_cast<char>(ch);
    }
  }
  return o + "\"";
}

std::string yjs1(const Yc &v, std::size_t i) {
  if (v.na[i]) return "null";
  if (v.k == 2) return v.d[i] != 0 ? "true" : "false";
  if (v.k == 0) return std::isfinite(v.d[i]) ? ynum::all(v.d[i]) : "\"" + ynum::all(v.d[i]) + "\"";
  return yjs1s(v.s[i]);
}

std::string yjsrows(const Yt &t, bool pretty, int ind) {
  if (t.n == 0) return "[]";
  std::string o = pretty ? "[\n" : "[";
  for (std::size_t i = 0; i < t.n; ++i) {
    if (i) o += pretty ? ",\n" : ",";
    if (pretty) o += std::string(static_cast<std::size_t>(ind) + 2, ' ') + "{\n";
    else o += "{";
    for (std::size_t j = 0; j < t.nm.size(); ++j) {
      if (j) o += pretty ? ",\n" : ",";
      if (pretty) o += std::string(static_cast<std::size_t>(ind) + 4, ' ');
      o += yjs1s(t.nm[j]) + (pretty ? ": " : ":") + yjs1(t.c[j], i);
    }
    o += pretty ? "\n" + std::string(static_cast<std::size_t>(ind) + 2, ' ') + "}" : "}";
  }
  return o + (pretty ? "\n" + std::string(static_cast<std::size_t>(ind), ' ') + "]" : "]");
}

// yoysmry() of the table t as a table of the same form
bool ysmry(Ctx &c, const Yt &t, double lag, Yt &s) {
  const Yc *g = t.get("group"), *p = t.get("prd"), *v = t.get("value");
  const Yc *ch = t.get("ppchg") ? t.get("ppchg") : t.get("pctchg");
  std::vector<std::string> keys;
  for (std::size_t i = 0; i < t.n; ++i) {
    const std::string k = g ? (g->na[i] ? "NA" : g->s[i]) : "";
    if (std::find(keys.begin(), keys.end(), k) == keys.end()) keys.push_back(k);
  }
  s = Yt();
  const char *cn[] = {"from", "to", "first", "last", "totalpct", "cagrpct", "periods", "up", "down", "flat"};
  if (g) {
    s.nm.push_back("group");
    Yc a;
    a.k = 1;
    s.c.push_back(a);
  }
  for (const char *k : cn) {
    s.nm.push_back(k);
    s.c.push_back(Yc());
  }
  const std::size_t o = g ? 1 : 0;
  for (const std::string &k : keys) {
    std::vector<std::size_t> ix;
    for (std::size_t i = 0; i < t.n; ++i) {
      if (!g || (g->na[i] ? "NA" : g->s[i]) == k) ix.push_back(i);
    }
    std::stable_sort(ix.begin(), ix.end(), [p](std::size_t a, std::size_t b) { return p->d[a] < p->d[b]; });
    std::vector<double> vv, cc;
    std::vector<std::size_t> fin;
    for (std::size_t i : ix) {
      vv.push_back(v->na[i] ? NaN : v->d[i]);
      cc.push_back(ch->na[i] ? NaN : ch->d[i]);
      if (!v->na[i] && !std::isnan(v->d[i])) fin.push_back(i);
    }
    double r[8];
    clab_ysum(vv.data(), cc.data(), static_cast<int>(vv.size()), lag, r);
    if (g) {
      s.c[0].s.push_back(k);
      s.c[0].d.push_back(NaN);
      s.c[0].na.push_back(0);
    }
    const double fr = fin.empty() ? NA : p->d[fin.front()], to = fin.empty() ? NA : p->d[fin.back()];
    const double vals[10] = {fr, to, r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]};
    for (std::size_t q = 0; q < 10; ++q) {
      s.c[o + q].d.push_back(vals[q]);
      s.c[o + q].s.emplace_back();
      s.c[o + q].na.push_back(std::isnan(vals[q]));
    }
  }
  s.n = keys.size();
  (void)c;
  return true;
}

bool ylag(Ctx &c, double &lag) {
  lag = c.has("lag") ? val1(c["lag"]) : 1;
  return (std::isfinite(lag) && lag > 0) || c.fail("`lag` must be one positive number.");
}

bool yoysmry(Ctx &c, Val &out) {
  Yt t, s;
  CK(ytab(c, t, true));
  double lag = 1;
  CK(ylag(c, lag));
  CK(ysmry(c, t, lag, s));
  out = vobj();
  for (std::size_t j = 0; j < s.nm.size(); ++j) {
    Val a = varr();
    for (std::size_t i = 0; i < s.n; ++i) {
      a.a.push_back(s.c[j].k == 1 ? vstr(s.c[j].s[i]) : s.c[j].na[i] ? Val() : vnum(s.c[j].d[i]));
    }
    put(out, s.nm[j].c_str(), a);
  }
  return true;
}

bool yoyjs(Ctx &c, Val &out) {
  std::string f, m;
  CK(yfile(c, f));
  bool pretty = true;
  CK(yflag(c, "pretty", true, pretty));
  Yt t;
  CK(ytab(c, t, false));
  if (c.has("digits")) {
    int d = 0;
    CK(ydig(c, "digits", 0, d, true));
    ypctround(t, d);
  }
  CK(ymeta(c, m));
  const bool chg = t.get("prd") && t.get("value") && t.get("prev") && t.get("change") &&
                   (t.get("pctchg") || t.get("ppchg"));
  std::vector<std::pair<std::string, std::string>> kv;
  if (!m.empty()) kv.emplace_back("meta", yjs1s(m));
  if (chg) {
    Yt s;
    CK(ysmry(c, t, 1, s));
    kv.emplace_back("summary", yjsrows(s, pretty, 2));
  }
  kv.emplace_back("rows", yjsrows(t, pretty, 2));
  std::string txt = pretty ? "{\n" : "{";
  for (std::size_t q = 0; q < kv.size(); ++q) {
    if (q) txt += pretty ? ",\n" : ",";
    txt += (pretty ? "  " : "") + yjs1s(kv[q].first) + (pretty ? ": " : ":") + kv[q].second;
  }
  txt += pretty ? "\n}\n" : "}\n";
  return ydone(c, txt, f, out);
}

bool fylbl(Ctx &c, Val &out) {
  std::string sep;
  CK(ytxt1(c, "sep", "/", sep));
  bool shrt = true;
  CK(yflag(c, "short", true, shrt));
  const Val &e = c["endyear"];
  std::vector<std::string> s;
  std::vector<char> na;
  bool num = isnum(e), txt = e.t == Val::STR || e.t == Val::ARR;
  for (const Val &x : e.a) txt = txt && (x.t == Val::STR || x.t == Val::NUL);
  if (!num && !txt) return c.fail("`endyear` must be numeric or character years.");
  if (num) {
    for (double y : dbl(e)) {
      if (!std::isnan(y) && y < 1) return c.fail("`endyear` must be positive years.");
    }
    for (double y : dbl(e)) {
      if (!std::isnan(y) && y != std::round(y)) return c.fail("`endyear` must be whole years.");
    }
  }
  Val o = varr();
  const std::vector<Val> one = {e};
  for (const Val &x : e.t == Val::ARR ? e.a : one) {
    long y = 0;
    bool ok = false;
    if (x.t == Val::NUM) {
      y = static_cast<long>(x.n);
      ok = true;
    } else if (x.t == Val::STR) {
      char *end = nullptr;
      const double v = std::strtod(x.s.c_str(), &end);
      ok = end != x.s.c_str() && *end == '\0' && std::isfinite(v);
      y = static_cast<long>(v);
    }
    if (!ok) {
      o.a.push_back(Val());
      continue;
    }
    char b[64];
    if (shrt) std::snprintf(b, sizeof b, "%ld%s%02ld", y - 1, sep.c_str(), y % 100);
    else std::snprintf(b, sizeof b, "%ld%s%ld", y - 1, sep.c_str(), y);
    o.a.push_back(vstr(b));
  }
  out = unbox(o);
  return true;
}

// yoylbl() with labels a named vector (an object) or one per period
bool yoylbl(Ctx &c, Val &out) {
  Yt t;
  CK(ytab(c, t, true));
  const Val &l = c["labels"];
  const Yc &p = *t.get("prd");
  std::vector<std::string> lab1(t.n);
  for (std::size_t i = 0; i < t.n; ++i) lab1[i] = p.na[i] ? "NA" : rfmt(p.d[i]);
  std::vector<std::string> out1 = lab1;
  if (l.t == Val::OBJ) {
    for (std::size_t i = 0; i < t.n; ++i) {
      for (std::size_t q = 0; q < l.k.size(); ++q) {
        std::string s;
        if (l.k[q] == lab1[i] && lab(l.a[q], s)) {
          out1[i] = s;
          break;
        }
      }
    }
  } else {
    std::vector<std::string> v;
    if (!labs(l, v)) return c.fail("`labels` must be a function of the periods or a vector of labels.");
    std::vector<std::string> u;
    for (const std::string &s : lab1) {
      if (std::find(u.begin(), u.end(), s) == u.end()) u.push_back(s);
    }
    if (v.size() != u.size()) {
      return c.fail("`labels` has %d entries for %d distinct periods; name them, or pass a function.",
                    static_cast<int>(v.size()), static_cast<int>(u.size()));
    }
    for (std::size_t i = 0; i < t.n; ++i) {
      out1[i] = v[static_cast<std::size_t>(std::find(u.begin(), u.end(), lab1[i]) - u.begin())];
    }
  }
  out = c["x"];
  Val a = varr();
  for (const std::string &s : out1) a.a.push_back(vstr(s));
  bool had = false;
  for (std::size_t j = 0; j < out.k.size(); ++j) {
    if (out.k[j] == "perlbl") {
      out.a[j] = a;
      had = true;
    }
  }
  if (!had) put(out, "perlbl", a);
  return true;
}

// Markdown's own characters, line breaks and HTML escaped, as .ymdesc();
// control bytes as \ooo or their names
std::string ymdesc(const std::string &s0) {
  std::string s;
  for (std::size_t i = 0; i < s0.size(); ++i) {
    const char ch = s0[i];
    if (ch == '\r' || ch == '\n') {
      if (s.empty() || (s0[i - 1] != '\r' && s0[i - 1] != '\n')) s += ' ';
      continue;
    }
    if (ch == '\\' || ch == '`' || ch == '*') s += '\\';
    s += ch;
  }
  std::string o;
  for (unsigned char ch : s) {
    if (ch < 32 || ch == 127) {
      static const char *const nm[] = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, "\\a",
                                       "\\b",   "\\t",   nullptr, nullptr, "\\f",   "\\r",   nullptr, nullptr};
      char b[8];
      if (ch < 16 && nm[ch]) o += nm[ch];
      else if (ch == 27) o += "\\033";
      else {
        std::snprintf(b, sizeof b, "\\%03o", ch);
        o += b;
      }
    } else if (ch == '&') o += "&amp;";
    else if (ch == '<') o += "&lt;";
    else if (ch == '>') o += "&gt;";
    else if (ch == '|') o += "\\|";
    else o += static_cast<char>(ch);
  }
  return o;
}

std::size_t uwidth(const std::string &s) {
  std::size_t n = 0;
  for (unsigned char ch : s) n += (ch & 0xC0) != 0x80;
  return n;
}

bool yoymd(Ctx &c, Val &out) {
  std::string f, m;
  CK(yfile(c, f));
  Yt t;
  CK(ytab(c, t, true));
  int dg = 1;
  CK(ydig(c, "digits", 1, dg, false));
  bool align = true;
  CK(yflag(c, "align", true, align));
  CK(ymeta(c, m));
  const Yc &val = *t.get("value");
  int vd = 0;
  if (!t.get("pctlo")) {
    for (std::size_t i = 0; i < t.n; ++i) {
      if (val.na[i] || !std::isfinite(val.d[i])) continue;
      int k = 0;
      while (k < 6 && std::fabs(val.d[i] - rround(val.d[i], k)) >= 1e-9) ++k;
      vd = std::max(vd, k);
    }
  }
  const bool pp = t.get("ppchg") != nullptr;
  const Yc &ch = *t.get(pp ? "ppchg" : "pctchg");
  const std::string dash = "\xe2\x80\x94";
  auto num = [&](const Yc &v, std::size_t i) { return v.na[i] || std::isnan(v.d[i]) ? dash : ynum::fix(v.d[i], vd); };
  auto pct = [&](double v, int d) { return std::isnan(v) ? dash : ynum::pct(v, d, pp ? "pp" : "%"); };
  std::vector<std::string> hd;
  std::vector<std::vector<std::string>> cols;
  const Yc *g = t.get("group"), *pl = t.get("perlbl"), *lo = t.get("pctlo"), *hi = t.get("pctup"),
           *fl = t.get("flag");
  auto addc = [&](const char *h) {
    hd.push_back(h);
    cols.emplace_back();
  };
  if (g) addc("group");
  addc("prd");
  addc("value");
  addc("previous");
  addc("change");
  addc(pp ? "points" : "change %");
  if (lo && hi) addc("interval");
  addc("note");
  for (std::size_t i = 0; i < t.n; ++i) {
    std::size_t q = 0;
    if (g) cols[q++].push_back(g->na[i] ? "NA" : g->s[i]);
    const Yc &p = *t.get("prd");
    cols[q++].push_back(pl ? (pl->na[i] ? "NA" : pl->s[i]) : (p.na[i] ? "NA" : rfmt(p.d[i])));
    cols[q++].push_back(num(val, i));
    cols[q++].push_back(num(*t.get("prev"), i));
    cols[q++].push_back(num(*t.get("change"), i));
    cols[q++].push_back(pct(ch.na[i] ? NaN : ch.d[i], dg));
    if (lo && hi) {
      cols[q++].push_back(lo->na[i] ? "" : "[" + pct(lo->d[i], 0) + ", " + pct(hi->d[i], 0) + "]");
    }
    cols[q++].push_back(fl && !fl->na[i] ? fl->s[i] : "");
  }
  for (auto &cv : cols) {
    for (std::string &s : cv) s = ymdesc(s);
  }
  for (std::string &s : hd) s = ymdesc(s);
  const std::size_t nc = hd.size(), nby = g ? 1 : 0;
  std::vector<std::size_t> w(nc, 0);
  if (align) {
    for (std::size_t j = 0; j < nc; ++j) {
      w[j] = std::max<std::size_t>(3, uwidth(hd[j]));
      for (const std::string &s : cols[j]) w[j] = std::max(w[j], uwidth(s));
    }
  }
  auto right = [&](std::size_t j) { return j > nby && hd[j] != "note"; };
  auto row = [&](const std::vector<std::string> &v) {
    std::string o = "| ";
    for (std::size_t j = 0; j < nc; ++j) {
      const std::size_t gap = w[j] > uwidth(v[j]) ? w[j] - uwidth(v[j]) : 0;
      if (j) o += " | ";
      o += right(j) ? std::string(gap, ' ') + v[j] : v[j] + std::string(gap, ' ');
    }
    return o + " |";
  };
  std::string txt = row(hd) + "\n| ";
  for (std::size_t j = 0; j < nc; ++j) {
    const std::size_t n = std::max<std::size_t>(3, w[j]);
    if (j) txt += " | ";
    txt += right(j) ? std::string(n - 1, '-') + ":" : std::string(n, '-');
  }
  txt += " |\n";
  for (std::size_t i = 0; i < t.n; ++i) {
    std::vector<std::string> v;
    for (std::size_t j = 0; j < nc; ++j) v.push_back(cols[j][i]);
    txt += row(v) + "\n";
  }
  if (!m.empty()) txt += "\n_" + ymdesc(m) + "_\n";
  return ydone(c, txt, f, out);
}

bool yoywr(Ctx &c, Val &out) {
  if (c["file"].t != Val::STR) return c.fail("`file` must be one file path.");
  std::string fm;
  CK(ytxt1(c, "format", "auto", fm));
  for (char &ch : fm) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
  const std::string &f = c["file"].s;
  if (fm == "auto") {
    std::string b = f.substr(f.find_last_of("/\\") == std::string::npos ? 0 : f.find_last_of("/\\") + 1);
    std::string e = b.find('.') == std::string::npos ? b : b.substr(b.find_last_of('.') + 1);
    for (char &ch : e) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    static const char *const ex[][2] = {{"html", "html"}, {"htm", "html"}, {"pdf", "pdf"},      {"csv", "csv"},
                                        {"tsv", "tsv"},   {"tab", "tsv"},  {"json", "json"},    {"md", "markdown"},
                                        {"markdown", "markdown"}};
    fm.clear();
    for (const auto &p : ex) {
      if (e == p[0]) fm = p[1];
    }
    if (fm.empty()) return c.fail("cannot tell the format from \"%s\"; pass `format`.", b.c_str());
  }
  if (fm == "html" || fm == "pdf") {
    return c.fail("the %s writer is R's alone: call yoy%s() in R.", fm.c_str(), fm == "html" ? "html" : "pdf");
  }
  bool ok;
  if (fm == "csv") ok = yoycsv(c, out);
  else if (fm == "tsv") ok = yoytsv(c, out);
  else if (fm == "json") ok = yoyjs(c, out);
  else if (fm == "markdown") ok = yoymd(c, out);
  else return c.fail("unknown format: %s.", fm.c_str());
  if (ok) out = vstr(f);
  return ok;
}
