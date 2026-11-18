// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Categorical integrity (R/catg.R) for clab_jrun(), with the checks and
// words of R. JSON carries no value labels: declbld() and xfrvrf() take
// the code book as an argument. Text is ordered by its bytes, as R orders
// it under LC_COLLATE=C.

struct Labs {
  std::vector<std::string> s;
  std::vector<char> na;
  bool num = false;
  std::vector<double> d;  // the numbers, when num
  std::size_t size() const { return s.size(); }
};

// .catx(): labels (text, numbers, true and false), null missing.
bool catx(Ctx &c, const char *nm, Labs &L) {
  L = Labs();
  const Val &o = c[nm];
  const std::string m = std::string("`") + nm + "` must be a vector of labels.";
  if (o.t == Val::NUL) return c.fail("%s", m.c_str());
  L.num = isnum(o);
  if (L.num) L.d = dbl(o);
  const std::vector<Val> one = {o};
  for (const Val &e : (o.t == Val::ARR || o.t == Val::OBJ) ? o.a : one) {
    std::string t;
    if (e.t == Val::NUL) {
      L.s.emplace_back();
      L.na.push_back(1);
    } else if (lab(e, t)) {
      L.s.push_back(t);
      L.na.push_back(0);
    } else {
      return c.fail("%s", m.c_str());
    }
  }
  return true;
}

// A named vector of labels: an object of text (or numbers).
struct Map {
  std::vector<std::string> k, v;
  const std::string *get(const std::string &x) const {
    for (std::size_t i = 0; i < k.size(); ++i) {
      if (k[i] == x) return &v[i];
    }
    return nullptr;
  }
};

bool catmap(Ctx &c, const char *nm, Map &m, bool text = true) {
  m = Map();
  const Val &o = c[nm];
  bool ok = o.t == Val::OBJ && !o.a.empty();
  for (std::size_t i = 0; ok && i < o.a.size(); ++i) {
    std::string t;
    ok = !o.k[i].empty() && (text ? o.a[i].t == Val::STR : lab(o.a[i], t));
    if (ok) {
      m.k.push_back(o.k[i]);
      m.v.push_back(text ? o.a[i].s : t);
    }
  }
  return ok || c.fail("`%s` must be a named character vector (old = new).", nm);
}

std::string cq(const std::string &s) { return "'" + s + "'"; }
std::string cqs(const std::vector<std::string> &v) {
  std::string o;
  for (const std::string &s : v) o += (o.empty() ? "" : ", ") + cq(s);
  return o;
}
bool has(const std::vector<std::string> &v, const std::string &s) {
  return std::find(v.begin(), v.end(), s) != v.end();
}
// unique non-missing labels, in the order of their first appearance
std::vector<std::string> seen(const Labs &L) {
  std::vector<std::string> u;
  for (std::size_t i = 0; i < L.size(); ++i) {
    if (!L.na[i] && !has(u, L.s[i])) u.push_back(L.s[i]);
  }
  return u;
}

Val lout(const Labs &L) {
  Val v = varr();
  for (std::size_t i = 0; i < L.size(); ++i) v.a.push_back(L.na[i] ? Val() : vstr(L.s[i]));
  return unbox(v);
}

// .recode()
bool recode(Ctx &c, Labs &L, const Map &m, const std::vector<std::string> &keep) {
  std::vector<std::string> un;
  for (const std::string &s : seen(L)) {
    if (!m.get(s) && !has(keep, s)) un.push_back(s);
  }
  if (!un.empty()) {
    return c.fail("values with no mapping: %s. Every observed category must be mapped by name (or listed in `keep`).",
                  cqs(un).c_str());
  }
  for (std::size_t i = 0; i < L.size(); ++i) {
    const std::string *t = L.na[i] ? nullptr : m.get(L.s[i]);
    if (t) L.s[i] = *t;
  }
  L.num = false;
  return true;
}

bool keeps(Ctx &c, std::vector<std::string> &k) {
  k.clear();
  return !c.has("keep") || labs(c["keep"], k) || c.fail("`keep` must be labels.");
}

bool grdrcd(Ctx &c, Val &out) {
  Labs L;
  Map m;
  std::vector<std::string> k;
  CK(catx(c, "x", L));
  CK(catmap(c, "map", m));
  CK(keeps(c, k));
  CK(recode(c, L, m, k));
  out = lout(L);
  return true;
}

bool deccodes(Ctx &c, Val &out) {
  Labs L;
  Map m;
  CK(catx(c, "x", L));
  CK(catmap(c, "vlabs", m));
  const bool kn = !c.has("keepna") || (c["keepna"].t == Val::BOOL && c["keepna"].b);
  for (char n : L.na) {
    if (n && !kn) return c.fail("`x` has missing codes and `keepna` is FALSE.");
  }
  CK(recode(c, L, m, {}));
  out = lout(L);
  return true;
}

bool relbl(Ctx &c, Val &out) {
  Labs L;
  if (!catx(c, "x", L)) {
    c.err.clear();
    return c.fail("`x` must be a vector of labels (text, a factor or numbers).");
  }
  const Val &mo = c["map"];
  bool named = mo.t == Val::OBJ && !mo.a.empty();
  for (const std::string &k : mo.k) named = named && !k.empty();
  if (!named) return c.fail("`map` must be named (old label = new label); labels are never assigned by position.");
  Map m;
  CK(catmap(c, "map", m, false));
  std::vector<std::string> k;
  CK(keeps(c, k));
  // the sorted distinct values: numbers by value, text by its bytes
  std::vector<std::size_t> ix;
  for (std::size_t i = 0; i < L.size(); ++i) {
    if (!L.na[i]) ix.push_back(i);
  }
  std::stable_sort(ix.begin(), ix.end(), [&L](std::size_t a, std::size_t b) {
    return L.num ? L.d[a] < L.d[b] : std::strcmp(L.s[a].c_str(), L.s[b].c_str()) < 0;
  });
  std::vector<std::string> old, un;
  for (std::size_t i : ix) {
    if (!has(old, L.s[i])) old.push_back(L.s[i]);
  }
  for (const std::string &s : old) {
    if (!m.get(s) && !has(k, s)) un.push_back(s);
  }
  if (!un.empty()) {
    return c.fail("levels with no mapping: %s. Name every level or list it in `keep`.", cqs(un).c_str());
  }
  CK(recode(c, L, m, k));
  out = lout(L);
  return true;
}

bool grdlvls(Ctx &c, Val &out) {
  Labs L;
  CK(catx(c, "x", L));
  std::vector<std::string> lv;
  const Val &lo = c["lvls"];
  bool ok = (lo.t == Val::STR || lo.t == Val::ARR) && labs(lo, lv) && !lv.empty();
  for (const Val &e : lo.a) ok = ok && e.t == Val::STR;
  if (!ok) return c.fail("`lvls` must be a non-empty character vector.");
  std::vector<std::string> stray;
  for (const std::string &s : seen(L)) {
    if (!has(lv, s)) stray.push_back(s);
  }
  if (!stray.empty()) return c.fail("values outside the declared levels: %s.", cqs(stray).c_str());
  std::string r;
  if (c.has("ref") && (!lab(c["ref"], r) || r != lv[0])) {
    return c.fail("the declared reference %s is not the first level (%s); put the reference first in `lvls`.",
                  cq(r).c_str(), cq(lv[0]).c_str());
  }
  out = lout(L);
  return true;
}

bool declbld(Ctx &c, Val &out) {
  Labs L;
  if (!catx(c, "x", L)) {
    c.err.clear();
    return c.fail("`x` must be a vector of codes.");
  }
  if (!c.has("vlabs")) {
    return c.fail("`x` carries no `labels` attribute and `vlabs` was not given; the codes alone are not the categories.");
  }
  Map m;
  CK(catmap(c, "vlabs", m));
  CK(recode(c, L, m, {}));
  out = lout(L);
  return true;
}

bool grdbin(Ctx &c, Val &out) {
  const Val &x = c["x"];
  std::string col;
  lab(c["col"], col);
  if (x.t == Val::NUL || x.t == Val::OBJ) return c.fail("column `%s` must be an atomic vector.", col.c_str());
  Labs L;
  if (!catx(c, "x", L)) {
    c.err.clear();
    return c.fail("column `%s` must be an atomic vector.", col.c_str());
  }
  bool txt = false;
  for (const Val &e : (x.t == Val::ARR ? x.a : std::vector<Val>{x})) txt = txt || e.t == Val::STR;
  if (txt) {
    std::vector<std::string> lv;
    for (std::size_t i = 0; i < L.size() && lv.size() < 4; ++i) {
      const std::string s = L.na[i] ? "NA" : L.s[i];
      if (!has(lv, s)) lv.push_back(s);
    }
    return c.fail("column %s is categorical (%s...). It is not coerced: as.numeric() on a factor returns level "
                  "indices, not the data. Encode it first, for example as.integer(x == \"treated\").",
                  cq(col).c_str(), cqs(lv).c_str());
  }
  std::vector<std::string> u;
  bool bad = false;
  for (std::size_t i = 0; i < L.size(); ++i) {
    if (L.na[i]) continue;
    const std::string s = L.num ? rfmt(L.d[i]) : L.s[i];
    const bool b01 = L.num ? (L.d[i] == 0 || L.d[i] == 1) : (s == "TRUE" || s == "FALSE");
    bad = bad || !b01;
    if (!has(u, s)) u.push_back(s);
  }
  if (bad) {
    std::string o;
    for (std::size_t i = 0; i < u.size() && i < 5; ++i) o += (i ? ", " : "") + u[i];
    return c.fail("column %s must be 0 or 1 (saw: %s).", cq(col).c_str(), o.c_str());
  }
  out = vbool(true);
  return true;
}

bool vrfrcd(Ctx &c, Val &out) {
  Labs o, r;
  CK(catx(c, "orig", o));
  CK(catx(c, "recd", r));
  if (o.size() != r.size()) {
    return c.fail("length mismatch (%d vs %d): rows were lost or duplicated during the recode.",
                  sat<int>(o.size()), sat<int>(r.size()));
  }
  if (o.na != r.na) return c.fail("missingness changed during the recode (values silently became NA, or NAs were filled).");
  Map d;
  if (c["decl"].t == Val::OBJ) catmap(c, "decl", d, false);
  c.err.clear();
  std::vector<std::string> ro, rc;
  for (std::size_t i = 0; i < o.size(); ++i) {
    if (o.na[i]) continue;
    if (!has(ro, o.s[i])) ro.push_back(o.s[i]);
    if (!has(rc, r.s[i])) rc.push_back(r.s[i]);
  }
  auto by = [](const std::string &a, const std::string &b) { return std::strcmp(a.c_str(), b.c_str()) < 0; };
  std::sort(ro.begin(), ro.end(), by);
  std::sort(rc.begin(), rc.end(), by);
  std::vector<double> f(ro.size() * rc.size(), 0.0);
  auto at = [](const std::vector<std::string> &v, const std::string &s) {
    return sat<std::size_t>(std::find(v.begin(), v.end(), s) - v.begin());
  };
  for (std::size_t i = 0; i < o.size(); ++i) {
    if (!o.na[i]) f[at(ro, o.s[i]) + ro.size() * at(rc, r.s[i])] += 1;
  }
  std::vector<std::string> fan;
  for (std::size_t a = 0; a < ro.size(); ++a) {
    int k = 0;
    for (std::size_t b = 0; b < rc.size(); ++b) k += f[a + ro.size() * b] > 0;
    if (k > 1) fan.push_back(ro[a]);
  }
  if (!fan.empty()) {
    return c.fail("original category mapped to MULTIPLE new categories: %s. "
                  "The recode is not a function of the category label.", cqs(fan).c_str());
  }
  for (std::size_t a = 0; a < ro.size(); ++a) {
    std::string got;
    for (std::size_t b = 0; b < rc.size(); ++b) {
      if (f[a + ro.size() * b] > 0) got = rc[b];
    }
    const std::string *w = d.get(ro[a]);
    const std::string want = w ? *w : ro[a];
    if (got != want) {
      return c.fail("%s was mapped to %s but the declared mapping says %s. "
                  "THIS is how groups get swapped; fix the recode before any model runs.", cq(ro[a]).c_str(), cq(got).c_str(),
                    cq(want).c_str());
    }
  }
  Val oo = varr(), rr = varr(), ff = varr();
  for (std::size_t b = 0; b < rc.size(); ++b) {
    for (std::size_t a = 0; a < ro.size(); ++a) {
      oo.a.push_back(vstr(ro[a]));
      rr.a.push_back(vstr(rc[b]));
      ff.a.push_back(vnum(f[a + ro.size() * b]));
    }
  }
  out = vobj();
  put(out, "original", oo);
  put(out, "recoded", rr);
  put(out, "Freq", ff);
  return true;
}

// The first permutation (in the order of R's recursion: each first
// element in turn, then the permutations of the rest) of the labels under
// which the observed counts, relabelled, match the published ones; pre
// holds it on success. Relabelled, label j has the count of the position
// whose new label is j.
bool pfind(std::vector<std::size_t> &pre, const std::vector<std::size_t> &rest, const std::vector<double> &obs,
           const std::vector<double> &pub, double tol) {
  const std::size_t k = obs.size();
  if (rest.empty()) {
    bool same = true, fit = true;
    for (std::size_t j = 0; j < k; ++j) same = same && pre[j] == j;
    for (std::size_t j = 0; j < k && fit; ++j) {
      const std::size_t q = sat<std::size_t>(std::find(pre.begin(), pre.end(), j) - pre.begin());
      fit = std::fabs(obs[q] - pub[j]) <= tol;
    }
    return !same && fit;
  }
  for (std::size_t i = 0; i < rest.size(); ++i) {
    std::vector<std::size_t> r2 = rest;
    pre.push_back(rest[i]);
    r2.erase(r2.begin() + sat<std::ptrdiff_t>(i));
    if (pfind(pre, r2, obs, pub, tol)) return true;
    pre.pop_back();
  }
  return false;
}

// vrfmrg() on the labels L; strict stops at a failure (fail() true)
bool mrg(Ctx &c, const Labs &L, Val &out, bool &failed) {
  failed = false;
  bool strict = !c.has("strict") || (c["strict"].t == Val::BOOL && c["strict"].b);
  double tol = 0;
  if (c.has("tol")) tol = val1(c["tol"]);
  if (!std::isfinite(tol) || tol < 0) return c.fail("`tol` must be one number from 0.");
  const Val &p = c[c.has("pub") ? "pub" : "soucnts"];
  bool ok = p.t == Val::OBJ && !p.a.empty() && isnum(p);
  for (std::size_t i = 0; ok && i < p.k.size(); ++i) ok = !p.k[i].empty() && p.a[i].t == Val::NUM;
  if (!ok) return c.fail("`pub` must be a named numeric vector.");
  const std::vector<std::string> &labs = p.k;
  const std::size_t k = labs.size();
  std::vector<double> obs(k, 0), pub(k);
  for (std::size_t j = 0; j < k; ++j) {
    pub[j] = p.a[j].n;
    for (std::size_t i = 0; i < L.size(); ++i) obs[j] += !L.na[i] && L.s[i] == labs[j];
  }
  std::vector<std::string> extra;
  for (const std::string &s : seen(L)) {
    if (!has(labs, s)) extra.push_back(s);
  }
  Val cv = vobj(), pv = vobj();
  for (std::size_t j = 0; j < k; ++j) {
    put(cv, labs[j].c_str(), vnum(obs[j]));
    put(pv, labs[j].c_str(), vnum(pub[j]));
  }
  out = vobj();
  put(out, "counts", cv);
  put(out, "published", pv);
  std::string msg;
  bool good = false;
  Val perm;
  if (!extra.empty()) {
    msg = "vrfmrg: labels present in the data but not in the published counts: " + cqs(extra) + ".";
  } else {
    good = true;
    for (std::size_t j = 0; j < k; ++j) good = good && std::fabs(obs[j] - pub[j]) <= tol;
    if (!good) {
      std::vector<std::string> hit;
      if (k > 1 && k <= 7) {
        std::vector<std::size_t> pre, all(k);
        for (std::size_t j = 0; j < k; ++j) all[j] = j;
        if (pfind(pre, all, obs, pub, tol)) {
          for (std::size_t j = 0; j < k; ++j) hit.push_back(labs[pre[j]]);
        }
      }
      std::string hint;
      if (!hit.empty()) {
        perm = vobj();
        std::string mv;
        for (std::size_t j = 0; j < k; ++j) {
          put(perm, labs[j].c_str(), vstr(hit[j]));
          if (hit[j] != labs[j]) mv += (mv.empty() ? "" : ", ") + labs[j] + " -> " + hit[j];
        }
        hint = " The observed counts match the published ones if the labels are permuted (" + mv + "): the labels are attached to the wrong groups.";
      }
      msg = "vrfmrg: recoded counts do not match the published counts. ";
      for (std::size_t j = 0; j < k; ++j) {
        msg += (j ? "; " : "") + labs[j] + ": observed " + rfmt(obs[j]) + ", published " + rfmt(pub[j]);
      }
      msg += "." + hint;
    }
  }
  put(out, "ok", vbool(good));
  put(out, "permutation", perm);
  put(out, "message", msg.empty() ? Val() : vstr(msg));
  if (!good && strict) {
    failed = true;
    c.err = msg;
    return false;
  }
  return true;
}

bool vrfmrg(Ctx &c, Val &out) {
  Labs L;
  CK(catx(c, "x", L));
  bool f = false;
  return mrg(c, L, out, f);
}

bool xfrvrf(Ctx &c, Val &out) {
  Labs L;
  CK(catx(c, "imported", L));
  const Val &cb = c["codebook"];
  if (c.has("codebook") && cb.t != Val::OBJ) return c.fail("`codebook` must be a named vector (code = label), or NULL.");
  if (c.has("codebook")) {
    Map m;
    if (!catmap(c, "codebook", m, false)) {
      c.err.clear();
      return c.fail("`codebook` must be a named vector (code = label), or NULL.");
    }
    // a refusal of the decoding or of the counts names xfrvrf, as R's
    // .clas() names it
    if (!recode(c, L, m, {})) return false;
  }
  bool f = false;
  Val m;
  if (!mrg(c, L, m, f)) {
    if (f) c.err = "xfrvrf: " + c.err.substr(c.err.find(": ") + 2);
    return false;
  }
  const bool ok = col(m, "ok").b;
  out = vobj();
  put(out, "ok", vbool(ok));
  put(out, "decoded", lout(L));
  put(out, "marginals", m);
  put(out, "codebook", Val());
  Val rs = varr();
  if (!ok) rs = vstr(col(m, "message").s);
  put(out, "reasons", ok ? Val() : rs);
  return true;
}
