// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Checks (src/bind.cpp): a reported odds ratio, risk ratio or risk
// difference against the table it came from (relchk), and the audit of
// categorical variables before a model (catchk).

// tlab(): one side's labels, none missing or repeated, or "1", "2", ...
bool tlab(Ctx &c, const std::vector<std::string> &lab, const std::vector<int> &na, std::size_t n, const char *side,
          std::vector<std::string> &out) {
  if (lab.size() != n) {
    out = grp({}, n);
    return true;
  }
  for (std::size_t i = 0; i < n; ++i) {
    if (na[i]) return c.fail("`x` has a missing %s label.", side);
    for (std::size_t j = 0; j < i; ++j) {
      if (lab[j] == lab[i]) return c.fail("`x` names a %s twice: %s.", side, lab[i].c_str());
    }
  }
  out = lab;
  return true;
}

// The names of an object, with null values marked missing (a named
// character vector); none for anything else.
void onames(const Val &o, std::vector<std::string> &s, std::vector<int> &na) {
  s.clear();
  na.clear();
  if (o.t != Val::OBJ) return;
  s = o.k;
  na.assign(s.size(), 0);
}

// ctab(): relchk()'s table as counts, column-major, r by c, with labels.
bool ctab(Ctx &c, const Val &x, std::vector<double> &p, int &r, int &k, std::vector<std::string> &rn,
          std::vector<std::string> &cn) {
  static const char *const tforms =
      "a matrix, a 2-d table, a data frame of numeric columns or a list of numeric vectors of one length";
  std::size_t nr = 0, nc = 0;
  std::vector<std::string> rl, cl;
  std::vector<int> rna, cna;
  const bool lst = x.t == Val::OBJ;
  int dim[2] = {-1, -1};
  if (lst) {
    nc = x.a.size();
    for (std::size_t j = 0; j < nc; ++j) {
      const Val &e = x.a[j];
      int d2[2] = {-1, -1};
      if (!isnum(e) || (e.t == Val::ARR && shape(e, 2, d2, 0)) || (j > 0 && vlen(e) != nr)) {
        return c.fail("`x` must be counts: %s.", tforms);
      }
      if (j == 0) nr = vlen(e);
    }
    onames(x, cl, cna);
    if (nc > 0) onames(x.a[0], rl, rna);
  } else if (x.t == Val::ARR && shape(x, 2, dim, 0)) {
    nr = static_cast<std::size_t>(dim[0]);
    nc = static_cast<std::size_t>(dim[1]);
  } else {
    return c.fail("`x` must be counts: %s.", tforms);
  }
  if (nr < 2 || nc < 2) return c.fail("`x` must have at least two rows and two columns.");
  if (3.0 * static_cast<double>(nr) * static_cast<double>(nc) * (static_cast<double>(nr) + static_cast<double>(nc) - 2) >
      INT_MAX) {
    return c.fail("`x` is too large: its readings would outnumber the largest integer.");
  }
  p.assign(nr * nc, 0.0);
  for (std::size_t j = 0; j < nc; ++j) {
    for (std::size_t i = 0; i < nr; ++i) {
      const Val &e = lst ? x.a[j].a[i] : x.a[i].a[j];
      p[i + nr * j] = e.t == Val::NUM ? e.n : NA;
    }
  }
  for (double e : p) {
    if (std::isnan(e)) return c.fail("`x` has a missing count: the table must be complete.");
    if (!std::isfinite(e)) return c.fail("`x` must be finite counts.");
    if (e < 0) return c.fail("`x` has a negative count.");
  }
  CK(tlab(c, rl, rna, nr, "row", rn));
  CK(tlab(c, cl, cna, nc, "column", cn));
  for (std::size_t i = 0; i < nr; ++i) {
    double s = 0;
    for (std::size_t j = 0; j < nc; ++j) s += p[i + nr * j];
    if (s == 0) return c.fail("`x` has no count in row %s: a group with no units has no contrast.", rn[i].c_str());
  }
  for (std::size_t j = 0; j < nc; ++j) {
    double s = 0;
    for (std::size_t i = 0; i < nr; ++i) s += p[i + nr * j];
    if (s == 0) {
      return c.fail("`x` has no count in column %s: an outcome no unit has gives no contrast.", cn[j].c_str());
    }
  }
  r = static_cast<int>(nr);
  k = static_cast<int>(nc);
  return true;
}

// The fewest decimals, 0 to 15, at which v prints as itself; -1 if none.
int ndec(double v) {
  char b[400];
  for (int d = 0; d <= 15; ++d) {
    std::snprintf(b, sizeof b, "%.*f", d, v);
    if (std::strtod(b, nullptr) == v) return d;
  }
  return -1;
}

bool relchk(Ctx &c, Val &out) {
  static const char *const kinds[] = {"or", "rr", "rd"};
  int r, k;
  std::vector<double> tb;
  std::vector<std::string> rn, cn;
  CK(ctab(c, c["x"], tb, r, k, rn, cn));
  Dv vv;
  CK(num(c, c["value"], "value", vv, 1));
  const double val = vv[0];
  int q = 0;
  if (c.has("type")) {
    const Val &t = c["type"];
    q = -1;
    const Val *e = t.t == Val::ARR && t.a.size() == 1 ? &t.a[0] : &t;
    for (int i = 0; i < 3 && e->t == Val::STR; ++i) {
      if (e->s == kinds[i]) q = i;
    }
    if (q < 0) return c.fail("`type` must be \"or\", \"rr\" or \"rd\".");
  }
  if (c.has("digits") && c.has("tol")) return c.fail("give `digits` or `tol`, not both.");
  int dg = INT_MIN;
  double tl = 0;
  if (c.has("tol")) {
    Dv t;
    CK(num(c, c["tol"], "tol", t, 1, 0.0));
    tl = t[0];
  } else {
    int d = ndec(val);
    if (c.has("digits")) CK(whole(c, c["digits"], "digits", 0, d));
    if (d > 15) return c.fail("`digits` must be at most 15.");
    if (d >= 0) {
      dg = d;
      tl = 0.5 * std::pow(10.0, -d);
    }
  }
  int gi = 0, hi = 0, ei = 0;
  if (c.has("group")) CK(rref(c, c["group"], rn, "group", "row", gi));
  if (c.has("ref")) CK(rref(c, c["ref"], rn, "ref", "row", hi));
  if (c.has("event")) CK(rref(c, c["event"], cn, "event", "column", ei));
  --gi;
  --hi;
  --ei;
  if (gi >= 0 && gi == hi) return c.fail("`group` and `ref` must name two different rows.");
  const std::size_t m = 3 * static_cast<std::size_t>(r) * static_cast<std::size_t>(k) * static_cast<std::size_t>(r + k - 2);
  std::vector<int> gr(m), g(m), h(m), e(m), kd(m), mt(m);
  std::vector<double> v(m);
  int nm = 0;
  CK(kchk(c, clab_relchk(tb.data(), r, k, val, tl, gr.data(), g.data(), h.data(), e.data(), kd.data(), v.data(), mt.data(),
                         &nm)));
  static const char *const side[] = {"rows", "cols"};
  Val c0 = varr(), c1 = varr(), c2 = varr(), c3 = varr(), c4 = varr(), c7 = varr();
  int hold = 0;
  std::vector<int> cl(m);
  for (std::size_t j = 0; j < m; ++j) {
    const std::vector<std::string> &gl = gr[j] == 0 ? rn : cn, &el = gr[j] == 0 ? cn : rn;
    c0.a.push_back(vstr(side[gr[j]]));
    c1.a.push_back(vstr(kinds[kd[j]]));
    c2.a.push_back(vstr(gl[static_cast<std::size_t>(g[j])]));
    c3.a.push_back(vstr(gl[static_cast<std::size_t>(h[j])]));
    c4.a.push_back(vstr(el[static_cast<std::size_t>(e[j])]));
    cl[j] = gr[j] == 0 && kd[j] == q && (gi < 0 || g[j] == gi) && (hi < 0 || h[j] == hi) && (ei < 0 || e[j] == ei);
    if (cl[j] && mt[j] != 0) ++hold;
  }
  Frame rd;
  rd.add("grp", c0);
  rd.add("type", c1);
  rd.add("group", c2);
  rd.add("ref", c3);
  rd.add("event", c4);
  rd.num("value", v);
  rd.lgl("match", mt.data(), m);
  rd.lgl("claim", cl.data(), m);
  Val cv = vobj();
  put(cv, "group", gi < 0 ? Val() : vstr(rn[static_cast<std::size_t>(gi)]));
  put(cv, "ref", hi < 0 ? Val() : vstr(rn[static_cast<std::size_t>(hi)]));
  put(cv, "event", ei < 0 ? Val() : vstr(cn[static_cast<std::size_t>(ei)]));
  out = vobj();
  put(out, "verdict", vstr(hold > 0 ? "holds" : (nm > 0 ? "other" : "none")));
  put(out, "ok", vbool(hold > 0));
  put(out, "value", vnum(val));
  put(out, "tol", vnum(tl));
  put(out, "digits", dg == INT_MIN ? Val() : vnum(dg));
  put(out, "type", vstr(kinds[q]));
  put(out, "claim", cv);
  put(out, "reads", rd.o);
  return true;
}

// lfold(): white space trimmed and each inner run one space, ASCII letters
// in lower case.
std::string lfold(const std::string &s) {
  std::string o;
  bool sp = false;
  for (const char ch0 : s) {
    const unsigned char ch = static_cast<unsigned char>(ch0);
    if (ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' || ch == '\f' || ch == '\v') {
      sp = !o.empty();
      continue;
    }
    if (sp) o += ' ';
    sp = false;
    o += ch >= 'A' && ch <= 'Z' ? static_cast<char>(ch - 'A' + 'a') : static_cast<char>(ch);
  }
  return o;
}

// What catchk() finds in one variable.
struct Var {
  const char *type = "";
  bool cat = true, bin = false, codes = false;
  int n = 0, k = 0, dref = -1;
  const char *why = "";
  std::vector<std::string> lev;
  std::vector<int> cnt, twin;
  int s[4] = {0, 0, 0, 0};
};

// The R type of a JSON vector: logical (true, false, null), double
// (numbers and null), character (text, or text and numbers); a nested
// value is a list.
enum Jt { JLGL, JDBL, JSTR, JLIST };
Jt jtype(const Val &x) {
  if (x.t == Val::BOOL || x.t == Val::NUL) return JLGL;
  if (x.t == Val::NUM) return JDBL;
  if (x.t == Val::STR) return JSTR;
  bool num = false, str = false, lgl = false;
  for (const Val &e : x.a) {
    if (e.t == Val::ARR || e.t == Val::OBJ) return JLIST;
    num = num || e.t == Val::NUM;
    str = str || e.t == Val::STR;
    lgl = lgl || e.t == Val::BOOL;
  }
  if (str) return JSTR;
  if (num) return JDBL;
  (void)lgl;
  return JLGL;
}
std::vector<const Val *> elems(const Val &x) {
  std::vector<const Val *> v;
  if (x.t == Val::ARR || x.t == Val::OBJ) {
    for (const Val &e : x.a) v.push_back(&e);
  } else {
    v.push_back(&x);
  }
  return v;
}

bool vaud(Ctx &c, const Val &x, int nmin, Var &v) {
  const std::vector<const Val *> el = elems(x);
  if (el.size() > static_cast<std::size_t>(INT_MAX)) return c.fail("`x` is too long for crimlab.");
  v.n = static_cast<int>(el.size());
  const Jt t = jtype(x);
  int dim[2] = {-1, -1};
  if (t == JLIST) {
    v.type = x.t == Val::ARR && shape(x, 2, dim, 0) ? "matrix" : "list";
    v.cat = false;
  } else {
    static const char *const tn[] = {"logical", "double", "character"};
    v.type = tn[t];
  }
  if (v.cat && t == JDBL) {
    for (const Val *e : el) {
      if (e->t == Val::NUM && e->n != std::floor(e->n)) v.cat = false;
    }
  }
  if (!v.cat) {
    v.why = "not categorical";
    int na = 0;
    if (t != JLIST) {
      for (const Val *e : el) na += e->t == Val::NUL;
    }
    v.s[0] = na;
    return true;
  }
  std::vector<int> code(el.size(), -1);
  const bool num = t == JDBL;
  bool b01 = num;
  int seen = 0;
  for (const Val *e : el) {
    if (!b01) break;
    if (e->t != Val::NUM) continue;
    ++seen;
    b01 = e->n == 0 || e->n == 1;
  }
  v.bin = t == JLGL || (b01 && seen > 0);
  if (v.bin) {
    v.lev = {t == JLGL ? "FALSE" : "0", t == JLGL ? "TRUE" : "1"};
    v.k = 2;
    for (std::size_t i = 0; i < el.size(); ++i) {
      const Val *e = el[i];
      code[i] = e->t == Val::BOOL ? (e->b ? 1 : 0) : e->t == Val::NUM ? static_cast<int>(e->n) : -1;
    }
  } else {
    // the distinct values as text, by first appearance
    std::vector<std::string> s(el.size());
    std::vector<int> na(el.size(), 0);
    for (std::size_t i = 0; i < el.size(); ++i) {
      if (!lab(*el[i], s[i])) na[i] = 1;
    }
    for (std::size_t i = 0; i < el.size(); ++i) {
      if (na[i]) continue;
      std::size_t j = 0;
      for (; j < v.lev.size() && v.lev[j] != s[i]; ++j) {
      }
      if (j == v.lev.size()) v.lev.push_back(s[i]);
      code[i] = static_cast<int>(j);
    }
    v.k = static_cast<int>(v.lev.size());
  }
  std::vector<std::string> fl;
  for (const std::string &e : v.lev) fl.push_back(lfold(e));
  int nd;
  const std::vector<int> fk = fcode(fl, nd, nullptr);
  v.cnt.assign(static_cast<std::size_t>(v.k) + 1, 0);
  v.twin.assign(static_cast<std::size_t>(v.k) + 1, 0);
  std::vector<int> iw(static_cast<std::size_t>(v.k) + 1);
  CK(kchk(c, clab_catchk(code.data(), v.n, v.k, fk.data(), nmin, v.cnt.data(), v.twin.data(), v.s, iw.data())));
  v.codes = num && !v.bin && v.k >= 1 && v.k <= 10;
  if (t == JLGL) {
    v.dref = 0;
    v.why = "logical: TRUE is the treatment";
  } else if (v.bin) {
    v.dref = 0;
    v.why = "0 and 1: 1 is the treatment";
  } else {
    v.why = "no order: state it";
  }
  return true;
}

// catchk()'s verdict on a stated reference w for v.
const char *refst(const Var &v, const std::string &w, int nmin) {
  for (int j = 0; j < v.k; ++j) {
    if (v.lev[static_cast<std::size_t>(j)] != w) continue;
    if (v.cnt[static_cast<std::size_t>(j)] == 0) return "empty";
    if (v.cnt[static_cast<std::size_t>(j)] < nmin) return "low";
    return j == v.dref ? "default" : "other";
  }
  const std::string fw = lfold(w);
  for (const std::string &e : v.lev) {
    if (lfold(e) == fw) return "case";
  }
  return "absent";
}

bool catchk(Ctx &c, Val &out) {
  int mn;
  CK(wholed(c, "nmin", 5, 1, mn));
  const Val &x = c["x"];
  std::vector<Val> vars;
  std::vector<std::string> nms;
  if (x.t == Val::OBJ) {
    vars = x.a;
    for (std::size_t i = 0; i < x.k.size(); ++i) nms.push_back(x.k[i].empty() ? rfmt(static_cast<double>(i + 1)) : x.k[i]);
  } else {
    int dim[2] = {-1, -1};
    if (x.t == Val::ARR && shape(x, 2, dim, 0) && dim[0] > 1 && dim[1] > 1) return notvec(c, "x");
    if (x.t == Val::ARR) {
      for (const Val &e : x.a) {
        if (e.t == Val::ARR || e.t == Val::OBJ) {
          return c.fail("`x` must be a vector, a data frame or a list of vectors.");
        }
      }
    }
    vars.push_back(x);
    nms.push_back("x");
  }
  const std::size_t nv = vars.size();
  if (nv == 0) return c.fail("`x` has no variables.");
  for (std::size_t i = 0; i < nv; ++i) {
    for (std::size_t j = 0; j < i; ++j) {
      if (nms[i] == nms[j]) return c.fail("`x` names a variable twice: %s.", nms[i].c_str());
    }
  }
  std::vector<std::string> want(nv);
  std::vector<int> has(nv, 0);
  if (c.has("ref")) {
    const Val &rf = c["ref"];
    const std::vector<const Val *> re = elems(rf);
    if (re.empty()) return c.fail("`ref` must name at least one level.");
    for (std::size_t j = 0; j < re.size(); ++j) {
      std::string w;
      const Val *e = re[j]->t == Val::ARR && re[j]->a.size() == 1 ? &re[j]->a[0] : re[j];
      if (!lab(*e, w)) return c.fail("`ref` must be levels with no missing values, one per variable.");
      if (rf.t != Val::OBJ) {
        if (re.size() != 1) return c.fail("`ref` must be one level, or levels named by variable.");
        for (std::size_t i = 0; i < nv; ++i) {
          want[i] = w;
          has[i] = 1;
        }
        continue;
      }
      const auto it = std::find(nms.begin(), nms.end(), rf.k[j]);
      if (it == nms.end()) return c.fail("`ref` names no variable: %s.", rf.k[j].c_str());
      const std::size_t i = static_cast<std::size_t>(it - nms.begin());
      if (has[i]) return c.fail("`ref` names a variable twice: %s.", rf.k[j].c_str());
      want[i] = w;
      has[i] = 1;
    }
  }
  static const char *const vnm[] = {"var", "type", "n", "nna", "nlev", "unused", "low", "twin", "bin",
                                    "codes", "nolab", "dlab", "ref", "why", "want", "state", "ok"};
  std::vector<Val> col(17, varr());
  Val lv[7];
  for (Val &e : lv) e = varr();
  bool all = true;
  for (std::size_t i = 0; i < nv; ++i) {
    Var v;
    CK(vaud(c, vars[i], mn, v));
    const char *st = !has[i] ? nullptr : (v.cat ? refst(v, want[i], mn) : "absent");
    const std::string sst = st == nullptr ? "" : st;
    const bool good = st == nullptr || sst == "default" || sst == "other";
    const bool ok = v.cat ? v.s[0] == 0 && v.s[1] == 0 && v.s[2] == 0 && v.s[3] == 0 && !v.codes && v.k - v.s[1] >= 2 &&
                                good
                          : v.s[0] == 0 && good;
    all = all && ok;
    col[0].a.push_back(vstr(nms[i]));
    col[1].a.push_back(vstr(v.type));
    col[2].a.push_back(vnum(v.n));
    col[3].a.push_back(vnum(v.s[0]));
    col[4].a.push_back(v.cat ? vnum(v.k) : Val());
    col[5].a.push_back(v.cat ? vnum(v.s[1]) : Val());
    col[6].a.push_back(v.cat ? vnum(v.s[2]) : Val());
    col[7].a.push_back(v.cat ? vnum(v.s[3]) : Val());
    col[8].a.push_back(vbool(v.cat && v.bin));
    col[9].a.push_back(vbool(v.cat && v.codes));
    col[10].a.push_back(vnum(0));
    col[11].a.push_back(vnum(0));
    col[12].a.push_back(v.dref < 0 ? Val() : vstr(v.lev[static_cast<std::size_t>(v.dref)]));
    col[13].a.push_back(vstr(v.why));
    col[14].a.push_back(has[i] ? vstr(want[i]) : Val());
    col[15].a.push_back(st == nullptr ? Val() : vstr(st));
    col[16].a.push_back(vbool(ok));
    for (std::size_t j = 0; j < v.lev.size(); ++j) {
      lv[0].a.push_back(vstr(nms[i]));
      lv[1].a.push_back(vstr(v.lev[j]));
      lv[2].a.push_back(Val());
      lv[3].a.push_back(vnum(v.cnt[j]));
      lv[4].a.push_back(vbool(v.cnt[j] > 0 && v.cnt[j] < mn));
      lv[5].a.push_back(v.twin[j] < 0 ? Val() : vstr(v.lev[static_cast<std::size_t>(v.twin[j])]));
      lv[6].a.push_back(vbool(static_cast<int>(j) == v.dref));
    }
  }
  Frame vd, ld;
  for (int j = 0; j < 17; ++j) vd.add(vnm[j], col[static_cast<std::size_t>(j)]);
  static const char *const lnm[] = {"var", "level", "label", "n", "low", "twin", "ref"};
  for (int j = 0; j < 7; ++j) ld.add(lnm[j], lv[j]);
  out = vobj();
  put(out, "ok", vbool(all));
  put(out, "nmin", vnum(mn));
  put(out, "vars", vd.o);
  put(out, "levels", ld.o);
  return true;
}
