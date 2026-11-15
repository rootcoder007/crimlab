// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The Indigenous health, published-table, census and Calls to Action
// functions and polex() for clab_jrun(): the JSON arguments made values
// of igh::V (src/igh.h), the function run, and its result written as
// dev/clref.R writes the R result (a data frame an object of columns, a
// named vector an object, one unnamed value that value; attributes are
// not written).

namespace ij {

igh::V fromj(const Val &v) {
  igh::V o;
  o.js = true;
  switch (v.t) {
    case Val::NUL:
      return o;
    case Val::BOOL:
      o.t = igh::V::LGL;
      o.i.push_back(v.b);
      return o;
    case Val::NUM:
      o.t = igh::V::DBL;
      o.d.push_back(v.n);
      return o;
    case Val::STR:
      o.t = igh::V::STR;
      o.s.push_back(v.s);
      o.na.push_back(0);
      return o;
    case Val::OBJ:
      o.t = igh::V::LIST;
      o.nm = v.k;
      for (const Val &e : v.a) o.a.push_back(fromj(e));
      return o;
    default:
      break;
  }
  bool nest = false, str = false, nm = false;
  for (const Val &e : v.a) {
    nest = nest || e.t == Val::ARR || e.t == Val::OBJ;
    str = str || e.t == Val::STR;
    nm = nm || e.t == Val::NUM;
  }
  if (nest) {
    o.t = igh::V::LIST;
    for (const Val &e : v.a) o.a.push_back(fromj(e));
    return o;
  }
  if (str) {
    o.t = igh::V::STR;
    for (const Val &e : v.a) {
      o.na.push_back(e.t == Val::NUL);
      o.s.push_back(e.t == Val::STR ? e.s : e.t == Val::NUM ? rfmt(e.n) : e.t == Val::BOOL ? (e.b ? "TRUE" : "FALSE")
                                                                                          : "");
    }
  } else if (nm) {
    o.t = igh::V::DBL;
    for (const Val &e : v.a) o.d.push_back(e.t == Val::NUM ? e.n : e.t == Val::BOOL ? (e.b ? 1.0 : 0.0) : NA);
  } else {
    o.t = igh::V::LGL;
    for (const Val &e : v.a) o.i.push_back(e.t == Val::BOOL ? (e.b ? 1 : 0) : INT_MIN);
  }
  return o;
}

Val one(const igh::V &v, std::size_t j) {
  switch (v.t) {
    case igh::V::DBL:
      return vnum(v.d[j]);
    case igh::V::INT:
      return v.i[j] == INT_MIN ? Val() : vnum(v.i[j]);
    case igh::V::LGL:
      return v.i[j] == INT_MIN ? Val() : vbool(v.i[j] != 0);
    case igh::V::FAC:
      return v.i[j] == INT_MIN ? Val() : vstr(v.s[static_cast<std::size_t>(v.i[j] - 1)]);
    case igh::V::STR:
      return !v.na.empty() && v.na[j] ? Val() : vstr(v.s[j]);
    default:
      return Val();
  }
}

Val arr(const igh::V &v) {
  Val a = varr();
  for (std::size_t j = 0; j < v.len(); ++j) a.a.push_back(one(v, j));
  return a;
}

Val toj(const igh::V &v) {
  if (v.t == igh::V::NUL) return Val();
  if (v.t == igh::V::LIST) {
    Val o;
    if (v.df) {
      o.t = Val::OBJ;
      for (std::size_t j = 0; j < v.a.size(); ++j) put(o, v.nm[j].c_str(), arr(v.a[j]));
      return o;
    }
    if (!v.nm.empty()) {
      o.t = Val::OBJ;
      for (std::size_t j = 0; j < v.a.size(); ++j) put(o, v.nm[j].c_str(), toj(v.a[j]));
      return o;
    }
    o = varr();
    for (const igh::V &e : v.a) o.a.push_back(toj(e));
    return o;
  }
  if (v.nr >= 0 && v.nc >= 0) {
    Val o = varr();
    for (int r = 0; r < v.nr; ++r) {
      Val row = varr();
      for (int k = 0; k < v.nc; ++k) row.a.push_back(one(v, static_cast<std::size_t>(k * v.nr + r)));
      o.a.push_back(row);
    }
    return o;
  }
  if (!v.nm.empty()) {
    Val o;
    o.t = Val::OBJ;
    for (std::size_t j = 0; j < v.len(); ++j) put(o, v.nm[j].c_str(), one(v, j));
    return o;
  }
  if (v.len() == 1) return one(v, 0);
  return arr(v);
}

bool run(Ctx &c, const char *fn, Val &out) {
  std::vector<std::string> an;
  std::vector<igh::V> av;
  for (std::size_t j = 0; j < c.in->k.size(); ++j) {
    an.push_back(c.in->k[j]);
    av.push_back(fromj(c.in->a[j]));
  }
  igh::V r;
  std::string err;
  std::vector<std::string> warn;
  if (igh::call(fn, an, av, r, err, warn) != CLAB_OK) return c.fail("%s", err.c_str());
  for (const std::string &w : warn) c.warn.push_back(w);
  out = toj(r);
  return true;
}

}  // namespace ij

bool sctidy(Ctx &c, Val &out) { return ij::run(c, "sctidy", out); }
bool qltyrcd(Ctx &c, Val &out) { return ij::run(c, "qltyrcd", out); }
bool limfrom(Ctx &c, Val &out) { return ij::run(c, "limfrom", out); }
bool indgap(Ctx &c, Val &out) { return ij::run(c, "indgap", out); }
bool genexpo(Ctx &c, Val &out) { return ij::run(c, "genexpo", out); }
bool exptrnd(Ctx &c, Val &out) { return ij::run(c, "exptrnd", out); }
bool trccmpl(Ctx &c, Val &out) { return ij::run(c, "trccmpl", out); }
bool trcprog(Ctx &c, Val &out) { return ij::run(c, "trcprog", out); }
bool mhcont(Ctx &c, Val &out) { return ij::run(c, "mhcont", out); }
bool cwbscr(Ctx &c, Val &out) { return ij::run(c, "cwbscr", out); }
bool connassoc(Ctx &c, Val &out) { return ij::run(c, "connassoc", out); }
bool healchg(Ctx &c, Val &out) { return ij::run(c, "healchg", out); }
bool commidx(Ctx &c, Val &out) { return ij::run(c, "commidx", out); }
bool frmaud(Ctx &c, Val &out) { return ij::run(c, "frmaud", out); }
bool scests(Ctx &c, Val &out) { return ij::run(c, "scests", out); }
bool ordcmp(Ctx &c, Val &out) { return ij::run(c, "ordcmp", out); }
bool stdprev(Ctx &c, Val &out) { return ij::run(c, "stdprev", out); }
bool welprof(Ctx &c, Val &out) { return ij::run(c, "welprof", out); }
bool holbal(Ctx &c, Val &out) { return ij::run(c, "holbal", out); }
bool dualcont(Ctx &c, Val &out) { return ij::run(c, "dualcont", out); }
bool cnsprof(Ctx &c, Val &out) { return ij::run(c, "cnsprof", out); }
bool profshr(Ctx &c, Val &out) { return ij::run(c, "profshr", out); }
bool cnsstd(Ctx &c, Val &out) { return ij::run(c, "cnsstd", out); }
bool cnscwb(Ctx &c, Val &out) { return ij::run(c, "cnscwb", out); }
bool trcind(Ctx &c, Val &out) { return ij::run(c, "trcind", out); }
bool trchlthgap(Ctx &c, Val &out) { return ij::run(c, "trchlthgap", out); }
bool polex(Ctx &c, Val &out) { return ij::run(c, "polex", out); }
