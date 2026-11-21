// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// cramerv() (R/crv.R), drfttab() and drftcal() (R/drtab.R), valrules()
// with its ready-made rules (R/rul.R), and the alcohol and drug use
// functions (R/adx.R) for clab_jrun(). The
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

// ebac(), stddrk(), hvydrk(), lrdg(), auditc(), dast10(), cudit() and
// canast() (R/adx.R).

// Numbers given as one number or an array (null: missing); absent: none.
inline bool adnum(Ctx &c, const char *nm, std::vector<double> &v) {
  v.clear();
  const Val &o = c[nm];
  bool ok = o.t == Val::NUM || o.t == Val::ARR;
  for (std::size_t i = 0; ok && i < vlen(o); ++i) {
    const Val &e = o.t == Val::ARR ? o.a[i] : o;
    double z = 0;
    if (e.t == Val::NUM) v.push_back(e.n);
    else if (e.t == Val::NUL) v.push_back(NA);
    else if (snum(e, z)) v.push_back(z);
    else ok = false;
  }
  if (!ok || v.empty()) return ckstop(c, "`%s` must be numeric", nm);
  return true;
}

// One finite number, or NaN.
inline double adone(const Val &o) {
  const Val &e = one(o);
  return e.t == Val::NUM && std::isfinite(e.n) ? e.n : NaN;
}

inline double adat(const std::vector<double> &v, std::size_t i) { return v.size() == 1 ? v[0] : v[i]; }

bool ebac(Ctx &c, Val &out) {
  static const char *const un[] = {"kg", "lb"};
  int u = 0;
  CK(pickd(c, "unit", un, 2, 0, u));
  std::vector<double> d, w, h(1, 0.0), r;
  CK(adnum(c, "drinks", d));
  CK(adnum(c, "weight", w));
  if (c.has("hours")) CK(adnum(c, "hours", h));
  for (double x : d) {
    if (x < 0) return ckstop(c, "`drinks` must be numbers from 0");
  }
  for (double x : w) {
    if (x <= 0) return ckstop(c, "`weight` must be positive numbers");
  }
  for (double x : h) {
    if (x < 0) return ckstop(c, "`hours` must be numbers from 0");
  }
  std::vector<double> ht;
  if (c.has("height")) {
    CK(adnum(c, "height", ht));
    for (double x : ht) {
      if (x <= 0) return ckstop(c, "`height` must be positive numbers");
    }
  }
  if (!c.has("r")) {
    const Val &s = c["sex"];
    std::vector<int> male;
    bool ok = s.t == Val::NUL || s.t == Val::STR || (s.t == Val::ARR && !s.a.empty());
    if (s.t == Val::NUL) male.push_back(1);
    for (std::size_t i = 0; ok && s.t != Val::NUL && i < vlen(s); ++i) {
      const Val &e = s.t == Val::ARR ? s.a[i] : s;
      ok = e.t == Val::STR && (e.s == "male" || e.s == "female");
      if (ok) male.push_back(e.s == "male");
    }
    if (!ok) return ckstop(c, "`sex` must be \"male\" or \"female\"");
    if (ht.empty()) {
      for (int k : male) r.push_back(k ? 0.73 : 0.66);
    } else {
      // the factor of Seidl, Jensen and Alt (2000) from weight and height
      const std::size_t q = std::max(std::max(male.size(), w.size()), ht.size());
      for (std::size_t k : {male.size(), w.size(), ht.size()}) {
        if (k != 1 && k != q) return ckstop(c, "the arguments must have one length, or length 1");
      }
      for (std::size_t i = 0; i < q; ++i) {
        const double kg = u == 1 ? adat(w, i) * 0.45359237 : adat(w, i), cm = adat(ht, i);
        const double f = (male.size() == 1 ? male[0] : male[i]) ? 0.31608 - 0.004821 * kg + 0.004632 * cm
                                                                 : 0.31223 - 0.006446 * kg + 0.004466 * cm;
        if (f <= 0) return ckstop(c, "the weight and height give a factor that is not positive");
        r.push_back(std::isnan(f) ? NA : f);
      }
    }
  } else {
    CK(adnum(c, "r", r));
    for (double x : r) {
      if (x <= 0) return ckstop(c, "`r` must be positive numbers");
    }
  }
  const double beta = c.has("beta") ? adone(c["beta"]) : 0.015;
  if (!(beta >= 0)) return ckstop(c, "`beta` must be one number from 0");
  const double g = c.has("grams") ? adone(c["grams"]) : 14;
  if (!(g > 0)) return ckstop(c, "`grams` must be one positive number");
  const std::size_t n = std::max(std::max(d.size(), w.size()), std::max(h.size(), r.size()));
  for (std::size_t k : {d.size(), w.size(), h.size(), r.size()}) {
    if (k != 1 && k != n) return ckstop(c, "the arguments must have one length, or length 1");
  }
  std::vector<double> o(n);
  for (std::size_t i = 0; i < n; ++i) {
    const double kg = u == 1 ? adat(w, i) * 0.45359237 : adat(w, i);
    const double v = g * adat(d, i) / (10 * adat(r, i) * kg) - beta * adat(h, i);
    o[i] = std::isnan(v) ? NA : std::fmax(0.0, v);
  }
  out = unbox(vnums(o));
  return true;
}

bool stddrk(Ctx &c, Val &out) {
  static const char *const sd[] = {"ca", "us", "uk", "au", "who"};
  static const double gs[] = {13.45, 14, 8, 10, 10};
  int s = 0;
  CK(pickd(c, "std", sd, 5, 0, s));
  std::vector<double> v, a;
  CK(adnum(c, "vol", v));
  CK(adnum(c, "abv", a));
  for (double x : v) {
    if (x < 0) return ckstop(c, "`vol` must be numbers from 0");
  }
  for (double x : a) {
    if (x < 0 || x > 100) return ckstop(c, "`abv` must be per cent, from 0 to 100");
  }
  const double g = c.has("grams") ? adone(c["grams"]) : gs[s];
  if (!(g > 0)) return ckstop(c, "`grams` must be one positive number");
  const std::size_t n = std::max(v.size(), a.size());
  if ((v.size() != 1 && v.size() != n) || (a.size() != 1 && a.size() != n)) {
    return ckstop(c, "the arguments must have one length, or length 1");
  }
  std::vector<double> o(n);
  for (std::size_t i = 0; i < n; ++i) {
    const double x = adat(v, i) * adat(a, i) / 100 * 0.789 / g;
    o[i] = std::isnan(x) ? NA : x;
  }
  out = unbox(vnums(o));
  return true;
}

// .adxmat(): the answers, a matrix (an array of rows) or a data frame (an
// object of columns) of numbers, as rows of k columns; k is k1 or k2.
bool adxmat(Ctx &c, std::size_t k1, std::size_t k2, const char *what, std::vector<std::vector<double>> &m,
            std::size_t &k) {
  const Val &x = c["x"];
  m.clear();
  k = 0;
  bool ok = true;
  auto cell = [&](const Val &e, double &z) {
    if (e.t == Val::NUM) z = e.n;
    else if (e.t == Val::NUL) z = NA;
    else if (!snum(e, z)) ok = false;
  };
  if (x.t == Val::OBJ) {
    k = x.a.size();
    std::size_t n = 0;
    for (const Val &col : x.a) {
      ok = ok && (col.t == Val::ARR || col.t == Val::NUM || col.t == Val::NUL);
      n = std::max(n, vlen(col));
    }
    m.assign(n, std::vector<double>(k, NA));
    for (std::size_t j = 0; ok && j < k; ++j) {
      const Val &col = x.a[j];
      if (vlen(col) != n) ok = false;
      for (std::size_t i = 0; ok && i < n; ++i) cell(col.t == Val::ARR ? col.a[i] : col, m[i][j]);
    }
  } else if (x.t == Val::ARR && !x.a.empty() && x.a[0].t == Val::ARR) {
    k = x.a[0].a.size();
    for (const Val &row : x.a) {
      ok = ok && row.t == Val::ARR && row.a.size() == k;
      if (!ok) break;
      std::vector<double> r(k, NA);
      for (std::size_t j = 0; ok && j < k; ++j) cell(row.a[j], r[j]);
      m.push_back(r);
    }
  } else {
    ok = false;
  }
  if (!ok || (k != k1 && k != k2)) return ckstop(c, "`x` must be a numeric matrix or data frame of %s", what);
  return true;
}

// .adxrng(): items scored 0 to 4, the items t1 and t2 (from 1; 0: none)
// 0, 2 or 4.
bool adxrng(Ctx &c, const std::vector<std::vector<double>> &m, std::size_t k, std::size_t t1, std::size_t t2) {
  for (std::size_t j = 0; j < k; ++j) {
    const bool two = j + 1 == t1 || j + 1 == t2;
    for (const std::vector<double> &r : m) {
      const double v = r[j];
      if (std::isnan(v)) continue;
      if (two ? !(v == 0 || v == 2 || v == 4) : !(v == 0 || v == 1 || v == 2 || v == 3 || v == 4)) {
        return two ? ckstop(c, "item %d must be scored 0, 2 or 4", sat<int>(j) + 1)
                   : ckstop(c, "item %d must be scored 0 to 4", sat<int>(j) + 1);
      }
    }
  }
  return true;
}

inline double adsum(const std::vector<double> &r, std::size_t from, std::size_t to) {
  double s = 0;
  for (std::size_t j = from; j < to; ++j) {
    if (std::isnan(r[j])) return NA;
    s += r[j];
  }
  return s;
}

// A column of labels, a missing one null.
inline Val adlab(const std::vector<double> &t, const double *cut, const char *const *lab, int n) {
  Val v = varr();
  for (double x : t) {
    if (std::isnan(x)) {
      v.a.push_back(Val());
      continue;
    }
    int k = 0;
    while (k + 1 < n && x >= cut[k + 1]) ++k;
    v.a.push_back(vstr(lab[k]));
  }
  return v;
}

bool auditc(Ctx &c, Val &out) {
  std::vector<std::vector<double>> m;
  std::size_t k = 0;
  CK(adxmat(c, 10, 3, "10 columns (the AUDIT) or 3 (the AUDIT-C)", m, k));
  CK(adxrng(c, m, k, 9, 10));
  const std::size_t n = m.size();
  std::vector<int> male;
  const Val &s = c["sex"];
  if (s.t != Val::NUL) {
    bool ok = s.t == Val::STR || s.t == Val::ARR;
    for (std::size_t i = 0; ok && i < vlen(s); ++i) {
      const Val &e = s.t == Val::ARR ? s.a[i] : s;
      ok = e.t == Val::STR && (e.s == "male" || e.s == "female");
      if (ok) male.push_back(e.s == "male");
    }
    if (!ok || (male.size() != 1 && male.size() != n)) {
      return ckstop(c, "`sex` must be \"male\" or \"female\", one or one per row");
    }
  }
  std::vector<double> total, c3;
  std::vector<int> cpos;
  for (std::size_t i = 0; i < n; ++i) {
    c3.push_back(adsum(m[i], 0, 3));
    total.push_back(k == 10 ? adsum(m[i], 0, 10) : NA);
    if (male.empty() || std::isnan(c3[i])) cpos.push_back(-1);
    else cpos.push_back(c3[i] >= ((male.size() == 1 ? male[0] : male[i]) ? 4 : 3));
  }
  static const double cut[] = {0, 8, 16, 20};
  static const char *const lab[] = {"low", "hazardous", "harmful", "dependence"};
  Frame f;
  f.num("total", total);
  f.num("auditc", c3);
  f.add("zone", adlab(total, cut, lab, 4));
  f.lgl("cpos", cpos.data(), cpos.size());
  out = f.o;
  return true;
}

bool dast10(Ctx &c, Val &out) {
  std::vector<std::vector<double>> m;
  std::size_t k = 0;
  CK(adxmat(c, 10, 10, "10 columns (the DAST-10)", m, k));
  for (std::size_t j = 0; j < k; ++j) {
    for (const std::vector<double> &r : m) {
      if (!std::isnan(r[j]) && r[j] != 0 && r[j] != 1) {
        return ckstop(c, "item %d must be 0 (no) or 1 (yes)", sat<int>(j) + 1);
      }
    }
  }
  std::vector<double> total;
  for (std::vector<double> &r : m) {
    r[2] = 1 - r[2];
    total.push_back(adsum(r, 0, 10));
  }
  static const double cut[] = {0, 1, 3, 6, 9};
  static const char *const lab[] = {"none", "low", "moderate", "substantial", "severe"};
  Frame f;
  f.num("total", total);
  f.add("level", adlab(total, cut, lab, 5));
  out = f.o;
  return true;
}

bool cudit(Ctx &c, Val &out) {
  std::vector<std::vector<double>> m;
  std::size_t k = 0;
  CK(adxmat(c, 8, 8, "8 columns (the CUDIT-R)", m, k));
  CK(adxrng(c, m, k, 8, 0));
  std::vector<double> total;
  for (const std::vector<double> &r : m) total.push_back(adsum(r, 0, 8));
  static const double cut[] = {0, 8, 12};
  static const char *const lab[] = {"low", "hazardous", "disorder"};
  Frame f;
  f.num("total", total);
  f.add("level", adlab(total, cut, lab, 3));
  out = f.o;
  return true;
}

bool canast(Ctx &c, Val &out) {
  static const double w[6][5] = {{0, 2, 3, 4, 6}, {0, 3, 4, 5, 6}, {0, 4, 5, 6, 7},
                                 {0, 5, 6, 7, 8}, {0, 3, 6, 0, 0}, {0, 3, 6, 0, 0}};
  std::vector<std::vector<double>> m;
  std::size_t k = 0;
  CK(adxmat(c, 6, 6, "6 columns (the ASSIST items on cannabis)", m, k));
  for (std::size_t j = 0; j < 6; ++j) {
    const int top = j < 4 ? 4 : 2;
    for (std::vector<double> &r : m) {
      const double v = r[j];
      if (std::isnan(v)) continue;
      if (v != std::floor(v) || v < 0 || v > top) {
        return ckstop(c, "item %d must be answered 0 to %d", sat<int>(j) + 1, top);
      }
    }
  }
  std::vector<double> total;
  for (std::vector<double> &r : m) {
    for (std::size_t j = 0; j < 6; ++j) {
      if (!std::isnan(r[j])) r[j] = w[j][sat<std::size_t>(r[j])];
    }
    total.push_back(adsum(r, 0, 6));
  }
  static const double cut[] = {0, 4, 27};
  static const char *const lab[] = {"low", "moderate", "high"};
  Frame f;
  f.num("total", total);
  f.add("level", adlab(total, cut, lab, 3));
  out = f.o;
  return true;
}

// .adxflag(): drinks against the threshold of each sex, at or above it
// (atl) or above it.
bool adxflag(Ctx &c, Val &out, double fdef, double mdef, bool atl) {
  std::vector<double> d;
  CK(adnum(c, "drinks", d));
  for (double x : d) {
    if (x < 0) return ckstop(c, "`drinks` must be numbers from 0");
  }
  const Val &s = c["sex"];
  std::vector<int> male;
  bool ok = s.t == Val::STR || (s.t == Val::ARR && !s.a.empty());
  for (std::size_t i = 0; ok && i < vlen(s); ++i) {
    const Val &e = s.t == Val::ARR ? s.a[i] : s;
    ok = e.t == Val::STR && (e.s == "male" || e.s == "female");
    if (ok) male.push_back(e.s == "male");
  }
  if (!ok) return ckstop(c, "`sex` must be \"male\" or \"female\"");
  const double f = c.has("f") ? adone(c["f"]) : fdef, m = c.has("m") ? adone(c["m"]) : mdef;
  if (!(f >= 0)) return ckstop(c, "`f` must be one number from 0");
  if (!(m >= 0)) return ckstop(c, "`m` must be one number from 0");
  const std::size_t n = std::max(d.size(), male.size());
  if ((d.size() != 1 && d.size() != n) || (male.size() != 1 && male.size() != n)) {
    return ckstop(c, "the arguments must have one length, or length 1");
  }
  Val v = varr();
  for (std::size_t i = 0; i < n; ++i) {
    const double x = adat(d, i), t = (male.size() == 1 ? male[0] : male[i]) ? m : f;
    if (std::isnan(x)) v.a.push_back(Val());
    else v.a.push_back(vbool(atl ? x >= t : x > t));
  }
  out = unbox(v);
  return true;
}

bool hvydrk(Ctx &c, Val &out) { return adxflag(c, out, 4, 5, true); }
bool lrdg(Ctx &c, Val &out) { return adxflag(c, out, 2, 3, false); }

// svypci() and svyqual() (R/svyci.R).

bool svypci(Ctx &c, Val &out) {
  static const char *const ms[] = {"beta", "wilson", "logit", "wald"};
  int m = 0;
  CK(pickd(c, "method", ms, 4, 0, m));
  bool trunc = true;
  CK(flagj(c, "trunc", true, trunc));
  std::vector<double> est, se, n, df;
  CK(adnum(c, "est", est));
  CK(adnum(c, "se", se));
  CK(adnum(c, "n", n));
  if (c.has("df")) {
    CK(adnum(c, "df", df));
  } else {
    for (double v : n) df.push_back(std::isnan(v) ? NA : v - 1);
  }
  for (double v : est) {
    if (v < 0 || v > 1) return ckstop(c, "`est` must be proportions, from 0 to 1");
  }
  for (double v : se) {
    if (v < 0) return ckstop(c, "`se` must be numbers from 0");
  }
  for (double v : n) {
    if (v < 2 || v != std::floor(v)) return ckstop(c, "`n` must be whole numbers of at least 2");
  }
  for (double v : df) {
    if (v < 1) return ckstop(c, "`df` must be at least 1");
  }
  const double level = c.has("level") ? adone(c["level"]) : 0.95;
  if (!(level > 0 && level < 1)) return ckstop(c, "`level` must be one number in (0, 1)");
  std::size_t k = 0;
  for (std::size_t q : {est.size(), se.size(), n.size(), df.size()}) k = std::max(k, q);
  for (std::size_t q : {est.size(), se.size(), n.size(), df.size()}) {
    if (q != 1 && q != k) return ckstop(c, "the arguments must have one length, or length 1");
  }
  const double a = 1 - level, z = clab_qnorm(1 - a / 2);
  std::vector<double> lo(k, NA), hi(k, NA), ne(k, NA), ev(k, NA);
  for (std::size_t i = 0; i < k; ++i) {
    const double p = adat(est, i), s = adat(se, i), nn = adat(n, i), d = adat(df, i);
    ev[i] = std::isnan(p) ? NA : p;
    if (std::isnan(p) || std::isnan(s) || std::isnan(nn) || std::isnan(d)) continue;
    const double tq = clab_qt(1 - a / 2, d, 1);
    double e = nn;
    if (p > 0 && p < 1 && s > 0) {
      const double r = clab_qt(1 - a / 2, nn - 1, 1) / tq;
      e = p * (1 - p) / (s * s) * r * r;
      if (trunc) e = std::fmin(nn, e);
    }
    ne[i] = e;
    const double x = p * e;
    if (m == 0) {
      lo[i] = x > 0 ? clab_qbeta(a / 2, x, e - x + 1, 1) : 0;
      hi[i] = x < e ? clab_qbeta(1 - a / 2, x + 1, e - x, 1) : 1;
    } else if (m == 1) {
      const double mid = p + z * z / (2 * e), half = z * std::sqrt(p * (1 - p) / e + z * z / (4 * e * e));
      lo[i] = (mid - half) / (1 + z * z / e);
      hi[i] = (mid + half) / (1 + z * z / e);
    } else if (m == 2) {
      if (p > 0 && p < 1) {
        const double l = std::log(p / (1 - p)), h = tq * s / (p * (1 - p));
        lo[i] = 1 / (1 + std::exp(-(l - h)));
        hi[i] = 1 / (1 + std::exp(-(l + h)));
      }
    } else {
      lo[i] = p - tq * s;
      hi[i] = p + tq * s;
    }
  }
  Frame f;
  f.num("est", ev);
  f.num("lo", lo);
  f.num("hi", hi);
  f.num("neff", ne);
  out = f.o;
  return true;
}

// Numbers from 0 (null: missing), or the refusal.
inline bool adpos(Ctx &c, const char *nm, std::vector<double> &v) {
  CK(adnum(c, nm, v));
  for (double x : v) {
    if (x < 0) return ckstop(c, "`%s` must be numbers from 0", nm);
  }
  return true;
}

bool svyqual(Ctx &c, Val &out) {
  std::vector<double> cv, n, d;
  CK(adpos(c, "cv", cv));
  CK(adpos(c, "n", n));
  const bool hasd = c.has("d");
  if (hasd) CK(adpos(c, "d", d));
  static const char *const tn[] = {"cvlo", "cvhi", "nmin", "dmin"};
  double t[4] = {16.6, 33.3, 10, 50};
  for (int i = 0; i < 4; ++i) {
    if (c.has(tn[i])) t[i] = adone(c[tn[i]]);
    if (!(t[i] >= 0)) return ckstop(c, "`%s` must be one number from 0", tn[i]);
  }
  if (t[0] > t[1]) return ckstop(c, "`cvlo` must not be above `cvhi`");
  std::size_t k = std::max(cv.size(), n.size());
  if (hasd) k = std::max(k, d.size());
  if ((cv.size() != 1 && cv.size() != k) || (n.size() != 1 && n.size() != k) ||
      (hasd && d.size() != 1 && d.size() != k)) {
    return ckstop(c, "the arguments must have one length, or length 1");
  }
  Val v = varr();
  for (std::size_t i = 0; i < k; ++i) {
    const double x = adat(cv, i), nn = adat(n, i), dd = hasd ? adat(d, i) : 0;
    // few: 1 too few, 0 enough, -1 not known
    int few = std::isnan(nn) ? -1 : nn < t[2];
    if (hasd) few = few == 1 || (!std::isnan(dd) && dd < t[3]) ? 1 : (few == -1 || std::isnan(dd) ? -1 : 0);
    if (few == 1 || (!std::isnan(x) && x > t[1])) v.a.push_back(vstr("unacceptable"));
    else if (few == -1 || std::isnan(x)) v.a.push_back(Val());
    else v.a.push_back(vstr(x <= t[0] ? "acceptable" : "marginal"));
  }
  out = unbox(v);
  return true;
}

}  // namespace vr
