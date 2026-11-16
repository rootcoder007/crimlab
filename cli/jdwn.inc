// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// dwnexp() (src/dwn.cpp) for clab_jrun(): the R wrapper of R/dwn.R in
// C++, with its checks in their order and words. Periods and places are
// labels: numbers when every label is a number, else text, sorted as
// numbers or byte by byte.

namespace dw {

using Cols = std::vector<std::string>;

// The labels of an argument, their sorted unique values (numbers or
// text) and each value's code from 0; false when one is missing.
struct Lab {
  bool num = true;
  std::vector<double> nv;
  Cols sv;
  std::vector<int> code;
  std::size_t nlev = 0;
};

inline bool labs(const Val &o, std::size_t n, bool sorted, Lab &l) {
  const Val a = ax::asarr(o);
  if (a.a.size() != n) return false;
  for (const Val &e : a.a) {
    if (ax::vna(e)) return false;
    l.num = l.num && (e.t == Val::NUM || (e.t == Val::BOOL));
  }
  l.code.assign(n, 0);
  if (l.num) {
    std::vector<double> v;
    for (const Val &e : a.a) v.push_back(ax::nv(e));
    l.nv = v;
    if (sorted) std::sort(l.nv.begin(), l.nv.end());
    std::vector<double> u;
    for (double x : l.nv) {
      if (std::find(u.begin(), u.end(), x) == u.end()) u.push_back(x);
    }
    l.nv = u;
    for (std::size_t i = 0; i < n; ++i) {
      l.code[i] = static_cast<int>(std::find(u.begin(), u.end(), v[i]) - u.begin());
    }
    l.nlev = u.size();
  } else {
    Cols v;
    for (const Val &e : a.a) v.push_back(ax::vlab(e));
    Cols u;
    for (const std::string &s : v) {
      if (std::find(u.begin(), u.end(), s) == u.end()) u.push_back(s);
    }
    if (sorted) std::sort(u.begin(), u.end());
    l.sv = u;
    for (std::size_t i = 0; i < n; ++i) {
      l.code[i] = static_cast<int>(std::find(u.begin(), u.end(), v[i]) - u.begin());
    }
    l.nlev = u.size();
  }
  return true;
}

// The labels, each repeated: one level per place (ip) or per period.
inline Val labv(const Lab &l, std::size_t np, std::size_t ni, bool per) {
  Val v = varr();
  for (std::size_t p = 0; p < np; ++p) {
    for (std::size_t i = 0; i < ni; ++i) {
      const std::size_t k = per ? p : i;
      v.a.push_back(l.num ? vnum(l.nv[k]) : vstr(l.sv[k]));
    }
  }
  return v;
}

// match.arg() of rule and dwt, in the package's words.
inline bool pick(Ctx &c, const char *nm, const std::vector<const char *> &ch, int &k) {
  k = 0;
  if (!c.has(nm)) return true;
  const Val &o = c[nm];
  const Val *s = o.t == Val::ARR && o.a.size() == 1 ? &o.a[0] : &o;
  int hit = -1;
  if (s->t == Val::STR && !s->s.empty()) {
    for (std::size_t i = 0; hit < 0 && i < ch.size(); ++i) {
      if (std::strncmp(ch[i], s->s.c_str(), s->s.size()) == 0) hit = static_cast<int>(i);
    }
  }
  if (hit >= 0) {
    k = hit;
    return true;
  }
  std::string q;
  for (std::size_t i = 0; i < ch.size(); ++i) q += std::string(i ? ", \"" : "\"") + ch[i] + "\"";
  const std::string got = s->t == Val::STR ? "\"" + s->s + "\"" : "a value that is not text";
  return ckstop(c, "`%s` must be one of %s, not %s", nm, q.c_str(), got.c_str());
}

inline bool flg(Ctx &c, const char *nm, bool def, bool &b) {
  b = def;
  if (!c.has(nm)) return true;
  if (!ax::tf(c[nm], b)) return ckstop(c, "`%s` must be TRUE or FALSE", nm);
  return true;
}

}  // namespace dw

bool dwnexp(Ctx &c, Val &out) {
  int ru = 0, wt = 0;
  CK(dw::pick(c, "rule", {"sector", "wedge", "angle"}, ru));
  Nv sx, sy, sd, rx, ry, hf;
  CK(cknum(c, c["sx"], "sx", sx, false));
  CK(cknum(c, c["sy"], "sy", sy, false));
  CK(cknum(c, c["sdir"], "sdir", sd, false));
  const std::size_t ns = sx.n();
  if (sy.n() != ns || sd.n() != ns) return ckstop(c, "`sx`, `sy` and `sdir` must have one length");
  dw::Lab pl;
  if (c.has("sper") && !dw::labs(c["sper"], ns, true, pl)) {
    return ckstop(c, "`sper` must give the period of each source, with no NA");
  }
  CK(cknum(c, c["rx"], "rx", rx, false));
  CK(cknum(c, c["ry"], "ry", ry, false));
  const std::size_t nv = rx.n();
  if (ry.n() != nv) return ckstop(c, "`rx` and `ry` must have one length");
  dw::Lab il;
  if (c.has("rid") && !dw::labs(c["rid"], nv, false, il)) {
    return ckstop(c, "`rid` must give the place of each row of `rx`, with no NA");
  }
  if (c.has("half")) {
    CK(cknum(c, c["half"], "half", hf, true, 0, 180));
  } else {
    hf.v.push_back(30);
  }
  double rad = HUGE_VAL;
  if (c.has("rad")) {
    const std::vector<double> r = ax::nums(c["rad"]);
    if (!ax::numv(c["rad"]) || r.size() != 1 || std::isnan(r[0]) || !(r[0] > 0)) {
      return ckstop(c, "`rad` must be a positive number or Inf");
    }
    rad = r[0];
  }
  if (!c.has("sper")) {
    pl.nv = {1};
    pl.code.assign(ns, 0);
    pl.nlev = 1;
  }
  std::vector<double> px = rx.v, py = ry.v;
  std::vector<int> ro;
  std::size_t npl = nv;
  if (ru != 1) {
    if (!c.has("rid")) {
      il.nv.clear();
      for (std::size_t i = 0; i < nv; ++i) il.nv.push_back(static_cast<double>(i + 1));
      il.nlev = nv;
    } else if (il.nlev != nv) {
      return ckstop(c, "`rid` must name each place once");
    }
  } else {
    if (!c.has("rid")) return ckstop(c, "`rule = \"wedge\"` needs `rid`, the polygons");
    std::vector<int> cnt(il.nlev, 0);
    for (int k : il.code) ++cnt[zz(k)];
    for (int k : cnt) {
      if (k < 3) return ckstop(c, "each polygon in `rid` needs three vertices or more");
    }
    ro.assign(1, 0);
    px.clear();
    py.clear();
    for (std::size_t g = 0; g < il.nlev; ++g) {
      for (std::size_t i = 0; i < nv; ++i) {
        if (zz(il.code[i]) != g) continue;
        px.push_back(rx.v[i]);
        py.push_back(ry.v[i]);
      }
      ro.push_back(static_cast<int>(px.size()));
    }
    npl = il.nlev;
  }
  bool from = true, ll = false, wrap = true;
  CK(dw::flg(c, "from", true, from));
  CK(dw::flg(c, "lonlat", false, ll));
  CK(dw::flg(c, "wrap", true, wrap));
  CK(dw::pick(c, "dwt", {"mult", "div", "none"}, wt));
  if (ru == 1 && ll) return ckstop(c, "`rule = \"wedge\"` needs planar coordinates (`lonlat = FALSE`)");
  if (ll) {
    bool ok = true;
    for (double v : sy.v) ok = ok && std::fabs(v) <= 90;
    for (double v : ry.v) ok = ok && std::fabs(v) <= 90;
    if (!ok) return ckstop(c, "with `lonlat = TRUE`, `sy` and `ry` must be latitudes, -90 to 90");
  }
  const std::size_t np = pl.nlev;
  if (static_cast<double>(npl) * static_cast<double>(np) > kutil::maxcell || ns > 2147483647u ||
      px.size() > 2147483647u) {
    return ckstop(c, "the places times the periods must be at most %.0f", kutil::maxcell);
  }
  std::vector<double> to(ns), o1(npl * np), o2(ru == 2 ? npl * np : 0);
  for (std::size_t i = 0; i < ns; ++i) to[i] = from ? ax::rmod(sd.v[i] + 180, 360) : sd.v[i];
  const int rc = clab_dwnexp(ru, static_cast<int>(ns), sx.p(), sy.p(), pl.code.data(), to.data(),
                             static_cast<int>(npl), px.data(), py.data(), ro.empty() ? nullptr : ro.data(),
                             static_cast<int>(np), hf.v[0], rad, ll ? 1 : 0, wrap ? 1 : 0, o1.data(),
                             ru == 2 ? o2.data() : nullptr);
  if (rc != CLAB_OK) return c.fail("%s", clab_strerror(rc));
  Frame f;
  f.add("id", dw::labv(il, np, npl, false));
  f.add("per", dw::labv(pl, np, npl, true));
  if (ru == 2) {
    std::vector<double> zv(o1.size());
    for (std::size_t i = 0; i < zv.size(); ++i) {
      zv[i] = wt == 0 ? o1[i] * o2[i] : (wt == 1 ? o1[i] / o2[i] : o1[i]);
    }
    f.add("angle", env::nas(o1));
    f.add("dist", env::nas(o2));
    f.add("z", env::nas(zv));
  } else {
    std::vector<int> an(o1.size());
    for (std::size_t i = 0; i < an.size(); ++i) an[i] = o1[i] > 0 ? 1 : 0;
    f.num("n", o1);
    f.ints("any", an.data(), an.size());
  }
  out = f.o;
  return true;
}

bool tailwd(Ctx &c, Val &out) {
  const Val &d = c["dir"], &s = c["spd"];
  const std::vector<double> dv = ax::nums(d), sv = ax::nums(s);
  if (!ax::numv(d) || !ax::numv(s) || dv.size() != sv.size() || dv.empty()) {
    return ckstop(c, "`dir` and `spd` must be numeric vectors of one length");
  }
  for (double v : sv) {
    if (v < 0) return ckstop(c, "`spd` must hold wind speeds, 0 or more (or NA)");
  }
  Nv co;
  if (c.has("course")) {
    CK(cknum(c, c["course"], "course", co));
  } else {
    co.v.push_back(250);
  }
  int agg = 0;
  CK(dw::pick(c, "agg", {"max", "mean"}, agg));
  const std::size_t n = dv.size();
  if (n > 2147483647u) return ckstop(c, "`dir` must have at most 2147483647 values");
  if (!c.has("by")) {
    std::vector<double> o(n);
    const int rc = clab_tailwd(static_cast<int>(n), dv.data(), sv.data(), co.v[0], nullptr, 0, 0, o.data());
    if (rc != CLAB_OK) return c.fail("%s", clab_strerror(rc));
    out = env::nas(o);
    return true;
  }
  dw::Lab g;
  if (!dw::labs(c["by"], n, true, g)) return ckstop(c, "`by` must give the group of each reading, with no NA");
  std::vector<double> o(g.nlev);
  const int rc = clab_tailwd(static_cast<int>(n), dv.data(), sv.data(), co.v[0], g.code.data(),
                             static_cast<int>(g.nlev), agg, o.data());
  if (rc != CLAB_OK) return c.fail("%s", clab_strerror(rc));
  Frame f;
  f.add("by", dw::labv(g, 1, g.nlev, false));
  f.add("tw", env::nas(o));
  out = f.o;
  return true;
}
