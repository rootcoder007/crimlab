// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Sentencing arithmetic, synthetic example data, published tables and
// the coding of police oversight reports (src/psc.cpp, src/sfx.cpp,
// src/ovr.cpp) for clab_jrun(): the R wrappers of R/psc.R, R/sfx.R,
// R/tabrpt.R and R/ovr.R in C++, with their checks in their order and
// words, and their results written as dev/clref.R writes R's.

namespace ovr {

// .pscone(): one value, or one per value of `days`.
inline bool pone(Ctx &c, const Nv &x, std::size_t n, const char *w) {
  if (x.n() != 1 && x.n() != n) return ckstop(c, "`%s` must have one value or one per value of `days`", w);
  return true;
}

bool pscred(Ctx &c, Val &out) {
  Nv d, r, s;
  CK(cknum(c, c["days"], "days", d, false, 0));
  if (c.has("ratio")) {
    CK(cknum(c, c["ratio"], "ratio", r, false, 0, 1.5));
  } else {
    r.v.push_back(1.5);
  }
  const std::size_t n = d.n();
  CK(pone(c, r, n, "ratio"));
  const bool hs = c.has("sentence");
  if (hs) {
    CK(cknum(c, c["sentence"], "sentence", s, false, 0));
    CK(pone(c, s, n, "sentence"));
  }
  std::vector<double> cr(n), rm(n), rr(n), ss(n);
  const int rc = clab_pscred(d.p(), static_cast<int64_t>(n), r.p(), static_cast<int64_t>(r.n()),
                             hs ? s.p() : nullptr, hs ? static_cast<int64_t>(s.n()) : 0, cr.data(), rm.data());
  if (rc != CLAB_OK) return kchk(c, rc);
  for (std::size_t i = 0; i < n; ++i) {
    rr[i] = r[r.n() == 1 ? 0 : i];
    if (hs) ss[i] = s[s.n() == 1 ? 0 : i];
  }
  out = vobj();
  put(out, "days", vnums(d.v));
  put(out, "ratio", vnums(rr));
  put(out, "credit", vnums(cr));
  if (hs) {
    put(out, "sentence", vnums(ss));
    put(out, "remain", vnums(rm));
  }
  return true;
}

// The presets of synfix() (src/sfx.cpp), column by column.
struct Pcol {
  std::string kind, name;
  int type, lo, hi;
  std::string lev;
};
inline std::vector<Pcol> ptab() {
  std::vector<Pcol> t;
  const char *k, *nm, *lv;
  int ty, lo, hi;
  for (int j = 0; clab_sfxcol(j, &k, &nm, &ty, &lo, &hi, &lv) == CLAB_OK; ++j) t.push_back({k, nm, ty, lo, hi, lv});
  return t;
}

// n rows drawn from the values of each column (vals[j], as JSON values).
inline bool sdraw(Ctx &c, const std::vector<std::string> &nm, const std::vector<std::vector<Val>> &vals, double n,
                  double seed, Val &out) {
  const std::size_t k = nm.size(), nn = static_cast<std::size_t>(n);
  std::vector<double> dn(k);
  for (std::size_t j = 0; j < k; ++j) dn[j] = static_cast<double>(vals[j].size());
  std::vector<int> code(nn * k);
  const int rc = clab_synfix(dn.data(), static_cast<int>(k), static_cast<int64_t>(nn), static_cast<uint32_t>(seed),
                             code.data());
  if (rc != CLAB_OK) return kchk(c, rc);
  out = vobj();
  for (std::size_t j = 0; j < k; ++j) {
    Val a = varr();
    a.a.reserve(nn);
    for (std::size_t i = 0; i < nn; ++i) a.a.push_back(vals[j][static_cast<std::size_t>(code[j * nn + i])]);
    put(out, nm[j].c_str(), std::move(a));
  }
  return true;
}

// .sfxpre(): the values of the columns of preset kind.
inline void spre(const std::vector<Pcol> &t, const std::string &kind, double n, std::vector<std::string> &nm,
                 std::vector<std::vector<Val>> &vals) {
  nm.clear();
  vals.clear();
  for (const Pcol &p : t) {
    if (p.kind != kind) continue;
    nm.push_back(p.name);
    std::vector<Val> v;
    if (p.type == 0) {
      std::size_t a = 0;
      for (;;) {
        const std::size_t b = p.lev.find('|', a);
        v.push_back(vstr(p.lev.substr(a, b == std::string::npos ? std::string::npos : b - a)));
        if (b == std::string::npos) break;
        a = b + 1;
      }
    } else if (p.type == 1) {
      for (int x = p.lo; x <= p.hi; ++x) v.push_back(vnum(x));
    } else {
      const double m = std::max(2.0, std::floor(n / 4));
      char b[32];
      for (double x = 1; x <= m; ++x) {
        std::snprintf(b, sizeof b, "syn-%05.0f", x);
        v.push_back(vstr(b));
      }
    }
    vals.push_back(std::move(v));
  }
}

bool synfix(Ctx &c, Val &out) {
  Nv n, sd;
  CK(cknumd(c, "n", 100, n, 1, 1e7, false, true));
  CK(cknumd(c, "seed", 1, sd, 0, 2147483647.0, false, true));
  const std::vector<Pcol> t = ptab();
  const Val &k = c["kind"];
  std::vector<std::string> nm;
  std::vector<std::vector<Val>> vals;
  if (k.t == Val::NUL) {
    static const char *const ty[] = {"labels", "range", "id"};
    Val a = varr(), b = varr(), d = varr(), e = varr();
    for (const Pcol &p : t) {
      a.a.push_back(vstr(p.kind));
      b.a.push_back(vstr(p.name));
      d.a.push_back(vstr(ty[p.type]));
      e.a.push_back(vstr(p.type == 1 ? std::to_string(p.lo) + "|" + std::to_string(p.hi) : p.lev));
    }
    out = vobj();
    put(out, "kind", std::move(a));
    put(out, "column", std::move(b));
    put(out, "type", std::move(d));
    put(out, "values", std::move(e));
    return true;
  }
  if (k.t == Val::OBJ) {
    // .sfxlist()
    bool ok = !k.k.empty();
    for (std::size_t j = 0; ok && j < k.k.size(); ++j) {
      ok = !k.k[j].empty() && std::find(k.k.begin(), k.k.begin() + static_cast<std::ptrdiff_t>(j), k.k[j]) ==
                                  k.k.begin() + static_cast<std::ptrdiff_t>(j);
    }
    if (!ok) return ckstop(c, "a list `kind` must name each column once");
    for (std::size_t j = 0; j < k.k.size(); ++j) {
      const Val &v = k.a[j];
      std::vector<Val> e;
      if (v.t == Val::ARR) e = v.a;
      else e.push_back(v);
      bool good = !e.empty() && e[0].t != Val::NUL && e[0].t != Val::ARR && e[0].t != Val::OBJ;
      for (const Val &x : e) good = good && x.t == e[0].t;
      if (!good) {
        return ckstop(c, "column `%s` of `kind` must be a vector of the values it takes, with no missing values",
                      k.k[j].c_str());
      }
      nm.push_back(k.k[j]);
      vals.push_back(std::move(e));
    }
    return sdraw(c, nm, vals, n[0], sd[0], out);
  }
  const Val *s = k.t == Val::ARR && k.a.size() == 1 ? &k.a[0] : &k;
  if (s->t != Val::STR) {
    return ckstop(c, "`kind` must be the name of a preset, a named list of the values of each column, or NULL");
  }
  std::vector<std::string> fam;
  for (const Pcol &p : t) {
    const std::string f = p.kind.substr(0, p.kind.find('-'));
    if (f == s->s && (fam.empty() || fam.back() != p.kind)) fam.push_back(p.kind);
  }
  if (!fam.empty()) {
    out = vobj();
    for (std::size_t i = 0; i < fam.size(); ++i) {
      spre(t, fam[i], n[0], nm, vals);
      Val f;
      CK(sdraw(c, nm, vals, n[0], std::fmod(sd[0] + static_cast<double>(i + 1), 2147483648.0), f));
      put(out, fam[i].substr(fam[i].find('-') + 1).c_str(), std::move(f));
    }
    return true;
  }
  spre(t, s->s, n[0], nm, vals);
  if (nm.empty()) return ckstop(c, "`kind` names no preset: \"%s\"; synfix() lists them", s->s.c_str());
  return sdraw(c, nm, vals, n[0], sd[0], out);
}

// .ovrin(): text, an array of text and nulls (null missing in miss).
inline bool otxt(Ctx &c, const char *nm, std::vector<std::string> &s, std::vector<char> &miss) {
  s.clear();
  miss.clear();
  const Val &o = c[nm];
  std::vector<const Val *> e;
  if (o.t == Val::ARR) {
    for (const Val &v : o.a) e.push_back(&v);
  } else {
    e.push_back(&o);
  }
  for (const Val *v : e) {
    if (v->t != Val::STR && v->t != Val::NUL) {
      return ckstop(c, "`%s` must be text (a character vector), not %s", nm, ax::vcls(o));
    }
    s.push_back(v->t == Val::STR ? v->s : std::string());
    miss.push_back(v->t == Val::NUL);
  }
  return true;
}

// .ovrhtml(): null (left out), true or false; -1, 1 or 0.
inline bool ohtml(Ctx &c, int &h) {
  const Val &o = c["html"];
  const Val *v = o.t == Val::ARR && o.a.size() == 1 ? &o.a[0] : &o;
  if (v->t == Val::NUL && o.t != Val::ARR) {
    h = -1;
    return true;
  }
  if (v->t != Val::BOOL && v->t != Val::NUL) return ckstop(c, "`html` must be NA, TRUE or FALSE");
  h = v->t == Val::NUL ? -1 : v->b;
  return true;
}

// The text of s as the rules read it (src/ovr.cpp).
inline std::string otext(const std::string &s, int h) {
  std::string o(s.size() + 1, '\0');
  int64_t len = 0;
  clab_ovrtxt(s.data(), static_cast<int64_t>(s.size()), h, &o[0], &len);
  o.resize(static_cast<std::size_t>(len));
  return o;
}

// Text and nulls: one value is that value.
inline Val ostrs(const std::vector<std::string> &s, const std::vector<char> &miss) {
  Val v = varr();
  for (std::size_t i = 0; i < s.size(); ++i) v.a.push_back(miss[i] ? Val() : vstr(s[i]));
  return unbox(v);
}

bool ovrtxt(Ctx &c, Val &out) {
  std::vector<std::string> s;
  std::vector<char> m;
  CK(otxt(c, "x", s, m));
  int h = -1;
  CK(ohtml(c, h));
  for (std::size_t i = 0; i < s.size(); ++i) {
    if (!m[i]) s[i] = otext(s[i], h);
  }
  out = ostrs(s, m);
  return true;
}

bool ovrdte(Ctx &c, Val &out) {
  std::vector<std::string> s;
  std::vector<char> m;
  CK(otxt(c, "x", s, m));
  Val v = varr();
  for (std::size_t i = 0; i < s.size(); ++i) {
    char b[16];
    if (!m[i]) clab_ovrdte(s[i].data(), static_cast<int64_t>(s[i].size()), b);
    v.a.push_back(m[i] || !b[0] ? Val() : vstr(b));
  }
  out = unbox(v);
  return true;
}

}  // namespace ovr
