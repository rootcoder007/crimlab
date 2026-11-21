// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// cramerv() (R/crv.R), drfttab() and drftcal() (R/drtab.R), and
// valrules() with its ready-made rules (R/rul.R) for clab_jrun(). The
// Monte Carlo p-value of a sparse table of labels is drawn by R's own
// generator, so drfttab() leaves such a column untested here and says so
// in its note, as drfthomog() refuses that branch. A rule that holds an R
// function (rul()) runs in R only. A pattern (rulerx()) is read as
// ECMAScript here and as Perl's in R; the common syntax is the same. A
// matcher that gives up on a pattern (libc++ and MSVC bound their work)
// is reported as a rule that could not be evaluated.

namespace vr {

inline const Val &vget(const Val &o, const char *k) {
  if (o.t == Val::OBJ) {
    for (std::size_t i = 0; i < o.k.size(); ++i) {
      if (o.k[i] == k) return o.a[i];
    }
  }
  return nul;
}

// One value: itself, or the one item of an array of one.
inline const Val &one(const Val &o) { return o.t == Val::ARR && o.a.size() == 1 ? o.a[0] : o; }

// One number; NaN when it is anything else.
inline double num1(const Val &o, double def) {
  if (o.t == Val::NUL) return def;
  const Val &e = one(o);
  return e.t == Val::NUM ? e.n : NaN;
}

bool cramerv(Ctx &c, Val &out) {
  bool bc = true;
  CK(flagj(c, "biascor", true, bc));
  const Val &t = c["tbl"];
  bool ok = t.t == Val::ARR && !t.a.empty();
  for (const Val &r : t.a) {
    ok = ok && r.t == Val::ARR;
    for (const Val &e : r.a) {
      double z = 0;
      ok = ok && (e.t == Val::NUM || e.t == Val::NUL || snum(e, z));
    }
  }
  if (!ok) return ckstop(c, "`tbl` must be a two-way table of counts (a numeric matrix or a 2-d table)");
  Val in = vobj();
  put(in, "row", t);
  if (c.has("level")) put(in, "level", c["level"]);
  Ctx s{c.fn, &in, std::string(), {}};
  Val x;
  if (!ctabx(s, x)) {
    c.err = s.err;
    return false;
  }
  for (const std::string &w : s.warn) c.warn.push_back(w);
  const Val &st = vget(x, "stat");
  auto f = [&](const char *k) {
    const Val &v = one(vget(st, k));
    return v.t == Val::NUM ? v.n : NA;
  };
  double mn = NA;
  bool any = false;
  for (const Val &r : vget(x, "expd").a) {
    for (const Val &e : r.a) {
      const double v = e.t == Val::NUM ? e.n : NA;
      if (!any) mn = v;
      else if (std::isnan(v) || std::isnan(mn)) mn = NA;
      else mn = std::fmin(mn, v);
      any = true;
    }
  }
  const double r[10] = {bc ? f("vbc") : f("v"), f("v"), f("vlo"), f("vhi"), f("chisq"), f("df"), f("p"),
                        f("n"),                 mn,     num1(vget(x, "drop"), NA)};
  static const char *const nm[] = {"v", "vraw", "lo", "hi", "chisq", "df", "p", "n", "minexp", "drop"};
  out = nvec(10, nm, r);
  return true;
}

// drfttab(): one column.
struct Dtrow {
  std::string type, test, note;
  double st = NA, p = NA, psi = NA;
  int nx = 0, ny = 0, flag = -1;
};
struct Dtarg {
  double alpha = 0.01, psithr = 0.25, psiminn = 100;
};

inline bool dtargs(Ctx &c, Dtarg &a) {
  a.alpha = num1(c["alpha"], 0.01);
  if (!(a.alpha > 0 && a.alpha < 1)) return ckstop(c, "`alpha` must be one number in (0, 1)");
  a.psithr = num1(c["psithr"], 0.25);
  if (!(a.psithr > 0) || std::isinf(a.psithr)) return ckstop(c, "`psithr` must be one positive number");
  a.psiminn = num1(c["psiminn"], 100);
  if (!(a.psiminn >= 1) || std::isinf(a.psiminn) || a.psiminn != std::floor(a.psiminn)) {
    return ckstop(c, "`psiminn` must be a whole number of at least 1");
  }
  return true;
}

// Numbers (0), labels (1), or nothing but missing values (2).
inline int dtkind(const Dcol &k) {
  if (k.ty != 0) return 1;
  for (char v : k.na) {
    if (!v) return 0;
  }
  return 2;
}

// Text tabulated, the labels sorted by their bytes.
inline void dttab(const std::vector<std::string> &s, std::vector<std::string> &lv, std::vector<double> &n) {
  lv = s;
  std::sort(lv.begin(), lv.end());
  lv.erase(std::unique(lv.begin(), lv.end()), lv.end());
  n.assign(lv.size(), 0);
  for (const std::string &v : s) n[sat<std::size_t>(std::lower_bound(lv.begin(), lv.end(), v) - lv.begin())] += 1;
}

// .dtcol(): rows ra of the reference column a against rows rb of the
// current column b.
bool dtcol(Ctx &c, const Dcol &a, const std::vector<std::size_t> &ra, const Dcol &b, const std::vector<std::size_t> &rb,
           const Dtarg &g, Dtrow &r) {
  r = Dtrow();
  int ka = dtkind(a), kb = dtkind(b);
  if (ka == 2) ka = kb == 2 ? 1 : kb;
  if (kb == 2) kb = ka;
  if (ka != kb) {
    r.type = "other";
    r.test = "none";
    r.nx = sat<int>(ra.size());
    r.ny = sat<int>(rb.size());
    r.note = "not numbers on both sides nor labels on both sides";
    return true;
  }
  if (ka == 0) {
    std::vector<double> x, y;
    for (std::size_t i : ra) {
      if (!a.na[i] && std::isfinite(a.x[i])) x.push_back(a.x[i]);
    }
    for (std::size_t i : rb) {
      if (!b.na[i] && std::isfinite(b.x[i])) y.push_back(b.x[i]);
    }
    r.type = "numeric";
    r.test = "ks";
    r.nx = sat<int>(x.size());
    r.ny = sat<int>(y.size());
    if (x.empty() || y.empty()) {
      r.note = "no values on one side";
      return true;
    }
    double k[3];
    clab_ks(x.data(), r.nx, y.data(), r.ny, k);
    r.st = k[0];
    r.p = k[1];
    const double bn = c.has("bins") ? val1(c["bins"]) : 10;
    if (std::isnan(bn) || bn < 2 || bn != std::floor(bn)) return c.fail("`bins` must be one whole number of at least 2.");
    if (bn > 1e6) return c.fail("`bins` must be at most 1,000,000.");
    std::vector<double> p(sat<std::size_t>(bn) + 2), q(p.size());
    const int nb = clab_psibin(x.data(), r.nx, y.data(), r.ny, sat<int>(bn), p.data(), q.data());
    double ps[2];
    clab_psi(p.data(), q.data(), nb, 1e-6, ps);
    r.psi = ps[0];
  } else {
    std::vector<std::string> sx, sy, lx, ly;
    std::vector<double> nx, ny, x, y;
    for (std::size_t i : ra) {
      if (!a.na[i]) sx.push_back(a.s[i]);
    }
    for (std::size_t i : rb) {
      if (!b.na[i]) sy.push_back(b.s[i]);
    }
    r.type = "label";
    r.test = "homogeneity";
    r.nx = sat<int>(sx.size());
    r.ny = sat<int>(sy.size());
    if (sx.empty() || sy.empty()) {
      r.note = "no values on one side";
      return true;
    }
    dttab(sx, lx, nx);
    dttab(sy, ly, ny);
    std::vector<std::string> lv = unionv(lx, ly);
    dralign(lv, lx, nx, x);
    dralign(lv, ly, ny, y);
    double h[4];
    clab_chih(x.data(), y.data(), sat<int>(lv.size()), h);
    r.st = h[0];
    const double n = static_cast<double>(sx.size() + sy.size());
    if (h[1] < 1) {
      r.note = "inapplicable: one category";
      return true;
    }
    if (n / h[3] < 2) {
      r.note = "inapplicable: " + rfmt(h[3]) + " categories for " + rfmt(n) +
               " observations (fewer than two per category)";
      return true;
    }
    if (h[2] < 5) {
      r.note = "small expected counts: the Monte Carlo p-value is drawn by R's random number generator "
               "(call drfttab() in R)";
      return true;
    }
    r.p = clab_pchisq(h[0], h[1], 0);
  }
  r.flag = (r.p < g.alpha || (r.psi > g.psithr && std::min(r.nx, r.ny) >= g.psiminn)) ? 1 : 0;
  return true;
}

// The columns of a table, each name once (its first column).
inline void dtnames(const std::vector<Dcol> &cs, std::vector<std::size_t> &at) {
  at.clear();
  for (std::size_t j = 0; j < cs.size(); ++j) {
    bool seen = false;
    for (std::size_t q : at) seen = seen || cs[q].nm == cs[j].nm;
    if (!seen) at.push_back(j);
  }
}

bool drfttab(Ctx &c, Val &out) {
  if (c["ref"].t != Val::OBJ || c["cur"].t != Val::OBJ) return ckstop(c, "`ref` and `cur` must be data frames");
  std::vector<Dcol> a, b;
  std::size_t na = 0, nb = 0;
  CK(dframe(c, c["ref"], a, na, "`ref` and `cur` must be data frames"));
  CK(dframe(c, c["cur"], b, nb, "`ref` and `cur` must be data frames"));
  Dtarg g;
  CK(dtargs(c, g));
  std::vector<std::size_t> ia, ib, ra(na), rb(nb);
  dtnames(a, ia);
  dtnames(b, ib);
  for (std::size_t i = 0; i < na; ++i) ra[i] = i;
  for (std::size_t i = 0; i < nb; ++i) rb[i] = i;
  std::vector<std::string> nm, ty, te, no;
  std::vector<double> st, p, psi;
  std::vector<int> nx, ny, fl;
  for (std::size_t j : ia) {
    for (std::size_t q : ib) {
      if (a[j].nm != b[q].nm) continue;
      Dtrow r;
      CK(dtcol(c, a[j], ra, b[q], rb, g, r));
      nm.push_back(a[j].nm);
      ty.push_back(r.type);
      te.push_back(r.test);
      no.push_back(r.note);
      st.push_back(r.st);
      p.push_back(r.p);
      psi.push_back(r.psi);
      nx.push_back(r.nx);
      ny.push_back(r.ny);
      fl.push_back(r.flag);
    }
  }
  if (nm.empty()) return ckstop(c, "`ref` and `cur` share no column");
  Frame f;
  f.str("column", nm);
  f.str("type", ty);
  f.str("test", te);
  f.num("statistic", st);
  f.num("p", p);
  f.num("psi", psi);
  f.ints("nref", nx.data(), nx.size());
  f.ints("ncur", ny.data(), ny.size());
  f.lgl("flag", fl.data(), fl.size());
  f.str("note", no);
  out = f.o;
  return true;
}

bool drftcal(Ctx &c, Val &out) {
  static const char *const msg = "`data` must be a data frame with at least four rows";
  if (c["data"].t != Val::OBJ) return ckstop(c, "%s", msg);
  std::vector<Dcol> d;
  std::size_t m = 0;
  CK(dframe(c, c["data"], d, m, msg));
  if (m < 4 || m > sat<std::size_t>(INT_MAX)) return ckstop(c, "%s", msg);
  const double n = num1(c["n"], 50);
  if (!(n >= 1) || n > 1e6 || n != std::floor(n)) return ckstop(c, "`n` must be a whole number from 1 to 1000000");
  double s = 0;
  CK(pseed(c, 1, s));
  Dtarg g;
  CK(dtargs(c, g));
  std::vector<std::size_t> at;
  dtnames(d, at);
  if (at.empty()) return ckstop(c, "`ref` and `cur` share no column");
  // each split has a seed of its own, drawn from the seed given
  std::vector<double> hit(at.size(), 0), u(m), each(sat<std::size_t>(n));
  clab_unif(s, sat<int>(n), each.data());
  std::vector<std::size_t> o(m), ra, rb;
  double any = 0;
  const int reps = sat<int>(n);
  for (int i = 0; i < reps; ++i) {
    clab_unif(std::floor(each[sat<std::size_t>(i)] * 4503599627370496.0), sat<int>(m), u.data());
    for (std::size_t q = 0; q < m; ++q) o[q] = q;
    std::stable_sort(o.begin(), o.end(), [&](std::size_t x, std::size_t y) { return u[x] < u[y]; });
    ra.assign(o.begin(), o.begin() + sat<std::ptrdiff_t>(m / 2));
    rb.assign(o.begin() + sat<std::ptrdiff_t>(m / 2), o.end());
    std::sort(ra.begin(), ra.end());
    std::sort(rb.begin(), rb.end());
    bool one1 = false;
    for (std::size_t j = 0; j < at.size(); ++j) {
      Dtrow r;
      CK(dtcol(c, d[at[j]], ra, d[at[j]], rb, g, r));
      if (r.flag == 1) {
        hit[j] += 1;
        one1 = true;
      }
    }
    if (one1) any += 1;
  }
  std::vector<std::string> nm;
  for (std::size_t j = 0; j < at.size(); ++j) {
    nm.push_back(d[at[j]].nm);
    hit[j] /= n;
  }
  Frame f;
  f.str("column", nm);
  f.num("far", hit);
  out = vobj();
  put(out, "columns", f.o);
  put(out, "anyflag", vnum(any / n));
  put(out, "alpha", vnum(g.alpha));
  put(out, "n", vnum(n));
  put(out, "alphafam", vnum(g.alpha / static_cast<double>(at.size())));
  return true;
}

// The ready-made rules. A rule is an object: rule (its kind), name, col
// (when it reads one column), sev, msg, then what it was given.
enum { RIN, RBTW, RNUL, RUNQ, RRX, RINC, RMAD, RCMP, RDST, RCNT, RKINDS };
const char *const rlkind[RKINDS] = {"inset", "btw", "null", "uniq", "rx", "incr", "mad", "cmpl", "dstn", "colcnt"};
const char *const rlfn[RKINDS] = {"ruleinset", "rulebtw",  "rulenull", "ruleuniq", "rulerx",
                                  "ruleincr",  "rulemad",  "rulecmpl", "ruledstn", "rulecolcnt"};
// The longest pattern and the longest value a pattern is tried on, in
// bytes: the matcher recurses with the length of what it reads.
const std::size_t RXPAT = 256, RXVAL = 4096;

inline bool rlstr(Ctx &c, const char *nm, std::string &s) {
  const Val &e = one(c[nm]);
  if (e.t != Val::STR || e.s.empty()) return ckstop(c, "`%s` must be one non-empty string", nm);
  s = e.s;
  return true;
}

inline bool rlone(Ctx &c, const char *nm, double &x) {
  const Val &e = one(c[nm]);
  if (e.t != Val::NUM || !std::isfinite(e.n)) return ckstop(c, "`%s` must be one number", nm);
  x = e.n;
  return true;
}

// A pattern compiled: 0 done, 1 not a pattern, 2 this build has none.
inline int rlrx(const std::string &pat, std::regex &re) {
#if defined(__cpp_exceptions) || defined(_CPPUNWIND)
  try {
    auto fl = std::regex::ECMAScript | std::regex::nosubs;
#ifdef __GLIBCXX__
    // the matcher of libstdc++ backtracks without a bound unless told to
    // run every state side by side: time then grows with the value times
    // the pattern, never exponentially (it has no back-references so)
    fl |= std::regex_constants::__polynomial;
#endif
    re = std::regex(pat, fl);
    return 0;
  } catch (const std::regex_error &) {
    return 1;
  }
#else
  (void)pat;
  (void)re;
  return 2;
#endif
}

bool mkrule(Ctx &c, int k, Val &out) {
  std::string col, name, msg;
  if (k <= RMAD) CK(rlstr(c, "col", col));
  static const char *const sv[] = {"warning", "fatal"};
  int si = 0;
  CK(pickd(c, "sev", sv, 2, 0, si));
  bool napass = true;
  const bool hasna = k == RIN || k == RBTW || k == RRX || k == RMAD;
  if (hasna) CK(flagj(c, "napass", true, napass));
  Val more = vobj();
  if (k == RIN) {
    const Val &o = c["set"];
    Val set = varr();
    bool ok = o.t != Val::NUL && o.t != Val::OBJ;
    for (std::size_t i = 0; ok && i < vlen(o); ++i) {
      const Val &e = o.t == Val::ARR ? o.a[i] : o;
      double z = 0;
      if (e.t == Val::STR) set.a.push_back(snum(e, z) ? vstr(rfmt(z)) : e);
      else if (e.t == Val::NUM) set.a.push_back(std::isnan(e.n) ? Val() : vstr(rfmt(e.n)));
      else if (e.t == Val::BOOL) set.a.push_back(vstr(e.b ? "TRUE" : "FALSE"));
      else if (e.t == Val::NUL) set.a.push_back(Val());
      else ok = false;
    }
    if (!ok || set.a.empty()) return ckstop(c, "`set` must be a vector of at least one allowed value");
    name = col + "_in_set";
    msg = "Column '" + col + "' has values outside the allowed set";
    put(more, "set", unbox(set));
  } else if (k == RBTW) {
    double lo = 0, hi = 0;
    CK(rlone(c, "lo", lo));
    CK(rlone(c, "hi", hi));
    if (lo > hi) return ckstop(c, "`lo` must not be above `hi`");
    name = col + "_between";
    msg = "Column '" + col + "' has values outside [" + rfmt(lo) + ", " + rfmt(hi) + "]";
    put(more, "lo", vnum(lo));
    put(more, "hi", vnum(hi));
  } else if (k == RNUL) {
    name = col + "_not_null";
    msg = "Column '" + col + "' has missing values";
  } else if (k == RUNQ) {
    name = col + "_unique";
    msg = "Column '" + col + "' has duplicate values";
  } else if (k == RRX) {
    std::string pat;
    CK(rlstr(c, "pat", pat));
    if (pat.size() > RXPAT) return ckstop(c, "`pat` must be at most %d bytes", sat<int>(RXPAT));
    std::regex re;
    const int e = rlrx(pat, re);
    if (e == 2) return ckstop(c, "this build has no regular expressions");
    if (e) return ckstop(c, "`pat` is not a valid regular expression");
    name = col + "_regex";
    msg = "Column '" + col + "' has values not matching /" + pat + "/";
    put(more, "pat", vstr(pat));
  } else if (k == RINC) {
    bool strict = false;
    CK(flagj(c, "strict", false, strict));
    name = col + "_increasing";
    msg = "Column '" + col + "' is not " + (strict ? "strictly increasing" : "non-decreasing");
    put(more, "strict", vbool(strict));
  } else if (k == RMAD) {
    double n = 3;
    if (c.has("n")) CK(rlone(c, "n", n));
    if (n <= 0) return ckstop(c, "`n` must be positive");
    name = col + "_within_mads";
    msg = "Column '" + col + "' has values more than " + rfmt(n) + " MADs from its median";
    put(more, "n", vnum(n));
  } else if (k == RCMP) {
    name = "complete_rows";
    msg = "Some rows have missing values";
  } else if (k == RDST) {
    name = "distinct_rows";
    msg = "Some rows are duplicated";
    if (c.has("cols")) {
      const Val &o = c["cols"];
      Val cols = varr();
      bool ok = o.t == Val::STR || o.t == Val::ARR;
      for (std::size_t i = 0; ok && i < vlen(o); ++i) {
        const Val &e = o.t == Val::ARR ? o.a[i] : o;
        ok = e.t == Val::STR && !e.s.empty();
        if (ok) cols.a.push_back(e);
      }
      if (!ok || cols.a.empty()) return ckstop(c, "`cols` must be NULL or a character vector of column names");
      msg += " on (";
      for (std::size_t i = 0; i < cols.a.size(); ++i) msg += (i ? ", " : "") + cols.a[i].s;
      msg += ")";
      put(more, "cols", unbox(cols));
    }
  } else {
    double n = 0;
    CK(rlone(c, "n", n));
    if (n < 0 || n != std::floor(n)) return ckstop(c, "`n` must be a whole number of at least 0");
    name = "col_count";
    char b[400];
    std::snprintf(b, sizeof b, "Expected %.0f columns", n);
    msg = b;
    put(more, "n", vnum(n));
  }
  out = vobj();
  put(out, "rule", vstr(rlkind[k]));
  put(out, "name", vstr(name));
  if (k <= RMAD) put(out, "col", vstr(col));
  put(out, "sev", vstr(sv[si]));
  put(out, "msg", vstr(msg));
  for (std::size_t i = 0; i < more.k.size(); ++i) {
    put(out, more.k[i].c_str(), more.a[i]);
    if (hasna && i + 1 == more.k.size()) put(out, "napass", vbool(napass));
  }
  return true;
}

bool ruleinset(Ctx &c, Val &out) { return mkrule(c, RIN, out); }
bool rulebtw(Ctx &c, Val &out) { return mkrule(c, RBTW, out); }
bool rulenull(Ctx &c, Val &out) { return mkrule(c, RNUL, out); }
bool ruleuniq(Ctx &c, Val &out) { return mkrule(c, RUNQ, out); }
bool rulerx(Ctx &c, Val &out) { return mkrule(c, RRX, out); }
bool ruleincr(Ctx &c, Val &out) { return mkrule(c, RINC, out); }
bool rulemad(Ctx &c, Val &out) { return mkrule(c, RMAD, out); }
bool rulecmpl(Ctx &c, Val &out) { return mkrule(c, RCMP, out); }
bool ruledstn(Ctx &c, Val &out) { return mkrule(c, RDST, out); }
bool rulecolcnt(Ctx &c, Val &out) { return mkrule(c, RCNT, out); }

// One cell as a key of equality: a missing value, text, or a number
// (0 and -0 the same).
inline std::string rlkey(const Dcol &k, std::size_t i) {
  if (k.na[i]) return std::string(1, '\x01');
  if (k.ty == 1) return "s" + k.s[i];
  char b[40];
  std::snprintf(b, sizeof b, "n%.17g", k.x[i] == 0 ? 0.0 : k.x[i]);
  return b;
}

inline double rlmed(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  const std::size_t n = v.size();
  return n % 2 ? v[n / 2] : (v[n / 2 - 1] + v[n / 2]) / 2;
}

// The rows that repeat an earlier key.
inline void rldup(const std::vector<std::string> &key, std::vector<int> &ok) {
  std::vector<std::size_t> o(key.size());
  for (std::size_t i = 0; i < o.size(); ++i) o[i] = i;
  std::stable_sort(o.begin(), o.end(), [&](std::size_t a, std::size_t b) { return key[a] < key[b]; });
  ok.assign(key.size(), 1);
  for (std::size_t i = 1; i < o.size(); ++i) {
    if (key[o[i]] == key[o[i - 1]]) ok[o[i]] = 0;
  }
}

// .rlrun(): rule r (from mkrule(), kind k) on the table. what: 0 it
// passes or does not apply, 1 an answer per row in ok, 2 it fails as a
// whole, 3 it could not be evaluated (why).
int rlrun(const Val &r, int k, const std::vector<Dcol> &cs, std::size_t n, std::vector<int> &ok, std::string &why) {
  ok.clear();
  const Dcol *v = nullptr;
  if (k <= RMAD) {
    const std::string &col = vget(r, "col").s;
    for (const Dcol &d : cs) {
      if (v == nullptr && d.nm == col) v = &d;
    }
    if (v == nullptr) return 0;
  }
  const bool napass = vget(r, "napass").b;
  const bool numrule = k == RBTW || k == RINC || k == RMAD;
  if (numrule && v->ty == 1) {
    why = "the column is not numeric";
    return 3;
  }
  if (k == RIN) {
    const Val &set = vget(r, "set");
    for (std::size_t i = 0; i < n; ++i) {
      bool in = false;
      for (std::size_t j = 0; !in && j < vlen(set); ++j) {
        const Val &e = set.t == Val::ARR ? set.a[j] : set;
        in = e.t == Val::STR && e.s == v->s[i];
      }
      ok.push_back(v->na[i] ? napass : in);
    }
    return 1;
  }
  if (k == RBTW) {
    const double lo = vget(r, "lo").n, hi = vget(r, "hi").n;
    for (std::size_t i = 0; i < n; ++i) ok.push_back(v->na[i] ? napass : v->x[i] >= lo && v->x[i] <= hi);
    return 1;
  }
  if (k == RNUL) {
    for (std::size_t i = 0; i < n; ++i) ok.push_back(!v->na[i]);
    return 1;
  }
  if (k == RUNQ) {
    std::vector<std::string> key;
    for (std::size_t i = 0; i < n; ++i) key.push_back(rlkey(*v, i));
    rldup(key, ok);
    return 1;
  }
  if (k == RRX) {
    std::regex re;
    if (rlrx(vget(r, "pat").s, re) != 0) {
      why = "the pattern could not be read";
      return 3;
    }
    for (std::size_t i = 0; i < n; ++i) {
      if (v->na[i]) {
        ok.push_back(napass);
        continue;
      }
      if (v->s[i].size() > RXVAL) {
        why = "a value is longer than " + rfmt(static_cast<double>(RXVAL)) + " bytes";
        return 3;
      }
#if defined(__cpp_exceptions) || defined(_CPPUNWIND)
      try {
        ok.push_back(std::regex_search(v->s[i], re));
      } catch (const std::regex_error &) {
        why = "the pattern is too complex for a value";
        return 3;
      }
#else
      ok.push_back(std::regex_search(v->s[i], re));
#endif
    }
    return 1;
  }
  if (k == RINC) {
    const bool strict = vget(r, "strict").b;
    bool first = true, good = true;
    double last = 0;
    for (std::size_t i = 0; i < n; ++i) {
      if (v->na[i]) continue;
      if (!first) good = good && (strict ? v->x[i] > last : v->x[i] >= last);
      last = v->x[i];
      first = false;
    }
    return good ? 0 : 2;
  }
  if (k == RMAD) {
    const double lim = vget(r, "n").n;
    std::vector<double> x;
    for (std::size_t i = 0; i < n; ++i) {
      if (!v->na[i]) x.push_back(v->x[i]);
    }
    if (x.empty()) {
      ok.assign(n, napass);
      return 1;
    }
    const double med = rlmed(x);
    for (double &e : x) e = std::fabs(e - med);
    const double dev = 1.4826 * rlmed(x);
    for (std::size_t i = 0; i < n; ++i) {
      if (v->na[i]) ok.push_back(napass);
      else ok.push_back(dev > 0 ? std::fabs(v->x[i] - med) <= lim * dev : v->x[i] == med);
    }
    return 1;
  }
  if (k == RCMP) {
    if (cs.empty()) return 0;
    for (std::size_t i = 0; i < n; ++i) {
      bool full = true;
      for (const Dcol &d : cs) full = full && !d.na[i];
      ok.push_back(full);
    }
    return 1;
  }
  if (k == RDST) {
    std::vector<const Dcol *> have;
    const Val &cols = vget(r, "cols");
    if (cols.t == Val::NUL) {
      std::vector<std::size_t> at;
      dtnames(cs, at);
      for (std::size_t j : at) have.push_back(&cs[j]);
    } else {
      for (std::size_t j = 0; j < vlen(cols); ++j) {
        const std::string &s = (cols.t == Val::ARR ? cols.a[j] : cols).s;
        const Dcol *f = nullptr;
        for (const Dcol &d : cs) {
          if (f == nullptr && d.nm == s) f = &d;
        }
        if (f != nullptr && std::find(have.begin(), have.end(), f) == have.end()) have.push_back(f);
      }
    }
    if (have.empty()) return 0;
    std::vector<std::string> key(n);
    for (std::size_t i = 0; i < n; ++i) {
      for (const Dcol *d : have) key[i] += rlkey(*d, i) + '\x1f';
    }
    rldup(key, ok);
    return 1;
  }
  return static_cast<double>(cs.size()) == vget(r, "n").n ? 0 : 2;
}

bool valrules(Ctx &c, Val &out) {
  static const char *const bad = "`rules` must be a rule or a list of rules, from rul() or the rule functions";
  std::vector<Dcol> cs;
  std::size_t n = 0;
  CK(dframe(c, c["data"], cs, n, "`data` must be a data frame"));
  const Val &rs = c["rules"];
  auto isr = [](const Val &v) { return v.t == Val::OBJ && vget(v, "rule").t == Val::STR; };
  std::vector<const Val *> rl;
  if (isr(rs)) {
    rl.push_back(&rs);
  } else if (rs.t == Val::ARR) {
    for (const Val &e : rs.a) {
      if (!isr(e)) return ckstop(c, "%s", bad);
      rl.push_back(&e);
    }
  } else {
    return ckstop(c, "%s", bad);
  }
  out = varr();
  for (const Val *r : rl) {
    const std::string &kind = vget(*r, "rule").s;
    if (kind == "fn") return ckstop(c, "a rule made by rul() holds an R function, and runs in R only");
    int k = 0;
    while (k < RKINDS && kind != rlkind[k]) ++k;
    if (k == RKINDS) return ckstop(c, "%s", bad);
    Ctx s{rlfn[k], r, std::string(), {}};
    Val spec;
    if (!mkrule(s, k, spec)) {
      c.err = s.err;
      return false;
    }
    std::vector<int> ok;
    std::string why;
    const int what = rlrun(spec, k, cs, n, ok, why);
    const std::string &msg = vget(spec, "msg").s;
    Val iss = vobj();
    put(iss, "severity", vget(spec, "sev"));
    if (what == 1) {
      std::vector<int> rows;
      int nf = 0;
      for (std::size_t i = 0; i < ok.size(); ++i) {
        if (ok[i]) continue;
        if (++nf <= 10) rows.push_back(sat<int>(i) + 1);
      }
      if (nf == 0) continue;
      put(iss, "message", vstr(msg + " (" + rfmt(nf) + " of " + rfmt(static_cast<double>(n)) + " rows)"));
      put(iss, "nfail", vnum(nf));
      put(iss, "rows", unbox(vints(rows)));
    } else if (what == 2) {
      put(iss, "message", vstr(msg));
    } else if (what == 3) {
      put(iss, "message", vstr(msg + " (the rule could not be evaluated: " + why + ")"));
    } else {
      continue;
    }
    const std::string &name = vget(spec, "name").s;
    bool had = false;
    for (std::size_t i = 0; i < out.k.size(); ++i) {
      if (out.k[i] == name) {
        out.a[i] = iss;
        had = true;
      }
    }
    if (!had) put(out, name.c_str(), iss);
  }
  return true;
}

}  // namespace vr
