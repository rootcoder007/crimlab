// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The 2021 Census Indigenous Population Profile read as published
// (cnsprof), shares with intervals that carry the long-form sample and
// random rounding (profshr), an Indigenous standard population (cnsstd),
// a Community Well-Being score from the profile (cnscwb); the indicators
// of the Calls to Action 19 and 55 (trcind) and their health gaps
// (trchlthgap); the simulated policy panel of the examples (polex); and
// the table of the functions with call() (igh.h).
#include <algorithm>
#include <climits>
#include <cmath>
#include <map>
#include <set>

#include "crimlab.h"
#include "igk.h"

namespace igh {

V f_sctidy(const Cx &);
V f_qltyrcd(const Cx &);
V f_limfrom(const Cx &);
V f_indgap(const Cx &);
V f_genexpo(const Cx &);
V f_exptrnd(const Cx &);
V f_trccmpl(const Cx &);
V f_trcprog(const Cx &);
V f_mhcont(const Cx &);
V f_cwbscr(const Cx &);
V f_connassoc(const Cx &);
V f_healchg(const Cx &);
V f_commidx(const Cx &);
V f_frmaud(const Cx &);
V f_scests(const Cx &);
V f_ordcmp(const Cx &);
V f_stdprev(const Cx &);
V f_welprof(const Cx &);
V f_holbal(const Cx &);
V f_dualcont(const Cx &);

namespace {

const int NAI = INT_MIN;
const std::vector<std::string> qlev = {"ok", "E", "F", "x", "..", "..."};

bool in(const std::string &s, const std::vector<std::string> &v) { return std::find(v.begin(), v.end(), s) != v.end(); }

// a column of the profile as labels, "NA" for missing (as paste() writes it)
std::vector<std::string> cl(const V &x, const char *nm) {
  std::vector<std::string> o(nrow(x), "NA");
  if (const V *c = col(x, nm)) {
    const Lab l = lab(*c);
    for (std::size_t j = 0; j < l.n(); ++j) {
      if (!l.na[j]) o[j] = l.s[j];
    }
  }
  return o;
}

std::vector<double> cd(const V &x, const char *nm) {
  const V *c = col(x, nm);
  return c ? num(*c, nm) : std::vector<double>(nrow(x), na());
}

std::vector<int> ci(const V &x, const char *nm) {
  const std::vector<double> d = cd(x, nm);
  if (bad()) return {};
  std::vector<int> o(d.size());
  for (std::size_t j = 0; j < d.size(); ++j) o[j] = std::isnan(d[j]) ? NAI : static_cast<int>(d[j]);
  return o;
}

// x has every column of a cnsprof() result
bool hascols(const V &x) {
  for (const char *n : {"geo", "dguid", "residence", "gender", "age", "group", "chaid", "characteristic", "parid",
                        "rootid", "sample", "value", "quality"}) {
    if (!col(x, n)) return false;
  }
  return true;
}

// x is a cnsprof() result (or, from JSON, a frame with its columns). The
// class is taken on trust here, so that the other arguments are checked
// first; a function that reads a column without a further check calls
// hascols() before it does.
bool isprof(const V &x) {
  if (inh(x, "cnsprof")) return true;
  return x.js && isdf(x) && hascols(x);
}

// .rmbl_ipp_groups
std::string grouplab(const std::string &b) {
  static const std::map<std::string, std::string> m = {
      {"TOTAL_IND_IDENT", "Total"},
      {"IND_IDENT", "Indigenous identity"},
      {"SINGLE_IND_RESP", "Single Indigenous responses"},
      {"FIRST_NATIONS", "First Nations"},
      {"M\xc3\x89TIS", "M\xc3\xa9tis"},
      {"INUK", "Inuit"},
      {"MULTI_IND_RESP", "Multiple Indigenous responses"},
      {"IND_RESP_NIE", "Indigenous responses not included elsewhere"},
      {"NON_IND", "Non-Indigenous"},
      {"TOTAL_IND_HHLD", "Total households"},
      {"IND_HHLD", "Indigenous households"},
      {"NON_IND_HHLD", "Non-Indigenous households"}};
  const auto f = m.find(b);
  return f == m.end() ? b : f->second;
}

// sub("\\s*\\([0-9]+\\)$", "", nm)
std::string nonum(const std::string &nm) {
  if (nm.empty() || nm.back() != ')') return nm;
  std::size_t b = nm.size() - 1;
  while (b > 0 && std::isdigit(static_cast<unsigned char>(nm[b - 1]))) --b;
  if (b == nm.size() - 1 || b == 0 || nm[b - 1] != '(') return nm;
  --b;
  while (b > 0 && std::isspace(static_cast<unsigned char>(nm[b - 1]))) --b;
  return nm.substr(0, b);
}

// a column of x as labels, or def where x has no such column
std::vector<std::string> orv(const V &x, const char *c0, const char *def) {
  const std::size_t n = nrow(x);
  std::vector<std::string> o(n, def);
  if (const V *v = col(x, c0)) {
    const Lab l = lab(*v);
    for (std::size_t j = 0; j < n; ++j) o[j] = l.na[j] ? "NA" : l.s[j];
  }
  return o;
}

typedef std::vector<std::size_t> Rows;

// the rows of a profile whose characteristic is name, at age ag (any age
// when ag is null), among the rows tot marks (all when tot is null) and the
// topics root (any when null)
Rows sel(const std::vector<std::string> &ch, const std::vector<std::string> &age, const std::vector<int> &rootid,
         const std::vector<char> *tot, const std::string &name, const char *ag, const std::vector<int> *root) {
  Rows o;
  for (std::size_t j = 0; j < ch.size(); ++j) {
    if ((tot && !(*tot)[j]) || ch[j] != name || (ag && age[j] != ag)) continue;
    if (root && std::find(root->begin(), root->end(), rootid[j]) == root->end()) continue;
    o.push_back(j);
  }
  return o;
}

}  // namespace

V f_cnsprof(const Cx &c) {
  const V &x = c["x"];
  nofn(x, "x");
  if (bad()) return {};
  if (!isdf(x)) return stop("`x` must be a data frame read from the profile");
  for (const char *n : {"GEO_NAME", "CHARACTERISTIC_ID", "CHARACTERISTIC_NAME"}) {
    if (!col(x, n)) {
      return stop("`x` is not an Indigenous Population Profile table (no GEO_NAME, CHARACTERISTIC_ID, CHARACTERISTIC_NAME)");
    }
  }
  std::vector<std::size_t> sym, vcol;
  std::vector<std::string> glab;
  for (std::size_t j = 0; j < x.nm.size(); ++j) {
    if (x.nm[j].compare(0, 7, "SYMBOL_") != 0) continue;
    const std::string b = x.nm[j].substr(7);
    for (std::size_t k = 0; k < x.nm.size(); ++k) {
      if (x.nm[k].compare(0, 7, "SYMBOL_") != 0 && nonum(x.nm[k]) == b) {
        sym.push_back(j);
        vcol.push_back(k);
        glab.push_back(grouplab(b));
        break;
      }
    }
  }
  if (sym.empty()) return stop("`x` has no group columns with their SYMBOL_ columns");
  const std::size_t n = nrow(x);
  const std::vector<int> cid = ci(x, "CHARACTERISTIC_ID");
  if (bad()) return {};
  const Lab raw = lab(*col(x, "CHARACTERISTIC_NAME"));
  // .rmbl_ipp_tree(): each characteristic once, in id order, with its level
  // (two spaces of indent each), parent, root and sample
  std::vector<int> id;
  std::vector<std::string> nm;
  {
    std::set<int> seen;
    std::vector<std::pair<int, std::string>> u;
    for (std::size_t j = 0; j < n; ++j) {
      if (seen.insert(cid[j]).second) u.push_back({cid[j], raw.na[j] ? "NA" : raw.s[j]});
    }
    std::stable_sort(u.begin(), u.end(), [](const auto &a, const auto &b) {
      if (a.first == NAI || b.first == NAI) return a.first != NAI && b.first == NAI;
      return a.first < b.first;
    });
    for (const auto &p : u) {
      id.push_back(p.first);
      nm.push_back(p.second);
    }
  }
  const std::size_t m = id.size();
  std::vector<int> level(m), parent(m, NAI), root(m, NAI);
  std::vector<std::string> name(m), sample(m);
  {
    std::vector<int> stack;
    int last = NAI;
    for (std::size_t i = 0; i < m; ++i) {
      std::size_t ws = 0;
      while (ws < nm[i].size() && (nm[i][ws] == ' ' || nm[i][ws] == '\t' || nm[i][ws] == '\n' || nm[i][ws] == '\r')) {
        ++ws;
      }
      level[i] = static_cast<int>(ws / 2);
      std::string s = nm[i];
      for (std::size_t p = s.find("\xc2\xa0"); p != std::string::npos; p = s.find("\xc2\xa0")) s.replace(p, 2, " ");
      name[i] = trim(s);
      const int lv = level[i];
      if (static_cast<int>(stack.size()) > lv) stack.resize(static_cast<std::size_t>(lv));
      parent[i] = lv > 0 && static_cast<int>(stack.size()) >= lv ? stack[static_cast<std::size_t>(lv - 1)] : NAI;
      root[i] = lv > 0 && !stack.empty() ? stack[0] : id[i];
      const bool tot = name[i].compare(0, 5, "Total") == 0;
      if (lv == 0 && tot) last = id[i];
      if (lv == 0 && !tot && last != NAI) root[i] = last;
      if (static_cast<int>(stack.size()) < lv + 1) stack.resize(static_cast<std::size_t>(lv + 1), NAI);
      stack[static_cast<std::size_t>(lv)] = id[i];
    }
    for (std::size_t i = 0; i < m; ++i) {
      const auto f = std::find(id.begin(), id.end(), root[i]);
      const std::string rn = f == id.end() ? "" : name[static_cast<std::size_t>(f - id.begin())];
      sample[i] = rn.find("25% sample") != std::string::npos ? "25%" : "100%";
    }
  }
  const std::size_t k = vcol.size();
  std::vector<std::string> geo, dguid, geolvl, res, gen, age, chn, smp, grp, q;
  std::vector<char> dgn, gln;
  std::vector<int> chaid, lvl, par, rt;
  std::vector<double> val, tsf, tlf;
  const Lab lgeo = lab(*col(x, "GEO_NAME"));
  const V *dg = col(x, "DGUID"), *gl = col(x, "GEO_LEVEL");
  const Lab ldg = dg ? lab(*dg) : Lab{std::vector<std::string>(n), std::vector<char>(n, 1)};
  const Lab lgl = gl ? lab(*gl) : Lab{std::vector<std::string>(n), std::vector<char>(n, 1)};
  const std::vector<std::string> vres = orv(x, "RESIDENCE", "Total"), vgen = orv(x, "GENDER", "Total"),
                                 vage = orv(x, "AGE_GROUP", "Total");
  std::vector<double> vsf(n, na()), vlf(n, na());
  for (int t = 0; t < 2; ++t) {
    if (const V *v = col(x, t ? "TNR_LF" : "TNR_SF")) {
      const Lab l = lab(*v);
      for (std::size_t j = 0; j < n; ++j) {
        double y = na();
        if (!l.na[j] && !tnum(l.s[j], y)) y = na();
        (t ? vlf : vsf)[j] = y;
      }
    }
  }
  std::vector<char> geon(lgeo.na);
  for (std::size_t g = 0; g < k; ++g) {
    const Lab ls = lab(x.a[sym[g]]), lv = lab(x.a[vcol[g]]);
    for (std::size_t j = 0; j < n; ++j) {
      geo.push_back(lgeo.s[j]);
      dguid.push_back(ldg.s[j]);
      dgn.push_back(ldg.na[j]);
      geolvl.push_back(lgl.s[j]);
      gln.push_back(lgl.na[j]);
      res.push_back(vres[j]);
      gen.push_back(vgen[j]);
      age.push_back(vage[j]);
      const std::size_t at = static_cast<std::size_t>(std::find(id.begin(), id.end(), cid[j]) - id.begin());
      chaid.push_back(cid[j]);
      chn.push_back(name[at]);
      lvl.push_back(level[at]);
      par.push_back(parent[at]);
      rt.push_back(root[at]);
      smp.push_back(sample[at]);
      grp.push_back(glab[g]);
      std::string s = ls.na[j] ? "" : trim(ls.s[j]);
      std::string qq = s.empty() ? "ok" : s;
      if (qq == "rE") qq = "E";
      if (!in(qq, qlev)) qq = "ok";
      double y = na();
      if (!lv.na[j] && !tnum(lv.s[j], y)) y = na();
      val.push_back(gone(qq) ? na() : y);
      q.push_back(qq);
      tsf.push_back(vsf[j]);
      tlf.push_back(vlf[j]);
    }
  }
  std::vector<char> gn;
  for (std::size_t g = 0; g < k; ++g) gn.insert(gn.end(), geon.begin(), geon.end());
  V out = frame({"geo", "dguid", "geolvl", "residence", "gender", "age", "chaid", "characteristic", "level", "parid",
                 "rootid", "sample", "group", "value", "quality", "tnrsf", "tnrlf"},
                {svna(geo, gn), svna(dguid, dgn), svna(geolvl, gln), sv(res), sv(gen), sv(age), iv(chaid), sv(chn),
                 iv(lvl), iv(par), iv(rt), sv(smp), sv(grp), dv(val), fac(q, qlev), dv(tsf), dv(tlf)});
  out.cls = {"cnsprof", "data.frame"};
  return out;
}

namespace {

// .rmbl_profshr_pick(): one characteristic by id or name
int chpick(const V &x, const V &v, const char *what) {
  const bool numv = (v.t == V::DBL || v.t == V::INT) && v.len() == 1;
  const bool strv = v.t == V::STR && v.s.size() == 1 && (v.na.empty() || !v.na[0]);
  double dv0 = numv ? (v.t == V::DBL ? v.d[0] : (v.i[0] == NAI ? na() : v.i[0])) : 0;
  if (!(numv || strv) || (numv && std::isnan(dv0))) return stop("`%s` must be one characteristic id or name", what);
  const std::vector<int> chaid = ci(x, "chaid");
  if (bad()) return 0;
  const std::vector<std::string> ch = cl(x, "characteristic");
  std::vector<int> ids;
  std::string shown;
  if (numv) {
    shown = rch(dv0);
    if (std::isfinite(dv0) && dv0 == std::round(dv0) && std::fabs(dv0) <= INT_MAX) ids.push_back(static_cast<int>(dv0));
    else ids.push_back(NAI);
  } else {
    shown = v.s[0];
    std::set<int> s;
    for (std::size_t j = 0; j < ch.size(); ++j) {
      if (ch[j] == v.s[0] && s.insert(chaid[j]).second) ids.push_back(chaid[j]);
    }
  }
  const bool there = ids.size() == 1 && std::find(chaid.begin(), chaid.end(), ids[0]) != chaid.end() && ids[0] != NAI;
  if (!there) {
    return stop("`%s` (%s) names %s characteristic of `x`", what, shown.c_str(), ids.size() > 1 ? "more than one" : "no");
  }
  return ids[0];
}

}  // namespace

V f_profshr(const Cx &c) {
  const V &x = c["x"];
  if (!isprof(x)) return stop("`x` must come from cnsprof()");
  const double smp = scalar(c["smpfrac"], "smpfrac", 0.001, 1, false, 0.25);
  const double deff = scalar(c["deff"], "deff", 0.1, 100, false, 1);
  const double rb = scalar(c["roubase"], "roubase", 1, 1000, true, 5);
  const double cf = conf(c["conflvl"], "conflvl", 0.95);
  if (bad()) return {};
  const int nid = chpick(x, c["num"], "num");
  if (bad()) return {};
  const std::vector<int> chaid = ci(x, "chaid"), parid = ci(x, "parid");
  if (bad()) return {};
  int did;
  if (c.has("den")) {
    did = chpick(x, c["den"], "den");
    if (bad()) return {};
  } else {
    std::vector<int> ps;
    for (std::size_t j = 0; j < chaid.size(); ++j) {
      if (chaid[j] == nid && std::find(ps.begin(), ps.end(), parid[j]) == ps.end()) ps.push_back(parid[j]);
    }
    if (ps.size() != 1 || ps[0] == NAI) return stop("the numerator has no parent: give `den`");
    did = ps[0];
  }
  const std::vector<std::string> dguid = cl(x, "dguid"), geo = cl(x, "geo"), res = cl(x, "residence"),
                                 gen = cl(x, "gender"), age = cl(x, "age"), grp = cl(x, "group"),
                                 smpv = cl(x, "sample"), ch = cl(x, "characteristic"), ql = cl(x, "quality");
  const std::vector<double> val = cd(x, "value");
  if (bad()) return {};
  if (!hascols(x)) return stop("`x` must come from cnsprof()");
  const Lab rl = lab(*col(x, "residence"));
  std::vector<std::size_t> a, b;
  std::map<std::string, std::size_t> kb;
  auto key = [&](std::size_t j) {
    return dguid[j] + "\r" + geo[j] + "\r" + res[j] + "\r" + gen[j] + "\r" + age[j] + "\r" + grp[j];
  };
  for (std::size_t j = 0; j < chaid.size(); ++j) {
    if (chaid[j] == nid) a.push_back(j);
    if (chaid[j] == did) {
      b.push_back(j);
      if (!kb.count(key(j))) kb[key(j)] = j;
    }
  }
  const double z = qn(cf), cr = (rb * rb - 1) / 6, NA = na();
  std::vector<double> nu, de, es, lo, hi;
  std::vector<std::string> qu;
  std::vector<char> qn2;
  for (std::size_t j : a) {
    const auto f = kb.find(key(j));
    const double num1 = val[j], den = f == kb.end() ? NA : val[f->second];
    const bool ok = std::isfinite(num1) && std::isfinite(den) && den > 0 && num1 >= 0 && num1 <= den;
    const double p = ok ? num1 / den : NA;
    const bool off = rl.na[j] || res[j] != "On reserve";
    const double fr = smpv[j] == "25%" && off ? smp : 1;
    const double v = ok ? deff * (1 - fr) * p * (1 - p) / (fr * den) + cr * (1 + p * p) / (den * den) : NA;
    const double ne = ok && v > 0 ? std::max(p * (1 - p), 1 / den) / v : NA;
    const double centre = (p + z * z / (2 * ne)) / (1 + z * z / ne);
    const double half = z / (1 + z * z / ne) * std::sqrt(p * (1 - p) / ne + z * z / (4 * ne * ne));
    nu.push_back(num1);
    de.push_back(den);
    es.push_back(100 * p);
    const double l0 = centre - half, h0 = centre + half;
    lo.push_back(std::isnan(l0) ? 100 * l0 : 100 * std::max(0.0, l0));
    hi.push_back(std::isnan(h0) ? 100 * h0 : 100 * std::min(1.0, h0));
    if (ql[j] == "ok") {
      qu.push_back(f == kb.end() ? "" : ql[f->second]);
      qn2.push_back(f == kb.end());
    } else {
      qu.push_back(ql[j]);
      qn2.push_back(0);
    }
  }
  auto pickl = [&](const char *n) {
    const Lab l = lab(*col(x, n));
    std::vector<std::string> o;
    std::vector<char> m;
    for (std::size_t j : a) {
      o.push_back(l.s[j]);
      m.push_back(l.na[j]);
    }
    return svna(o, m);
  };
  V out = frame({"geo", "dguid", "residence", "gender", "age", "group", "numerator", "denominator", "estimate", "lower",
                 "upper", "quality"},
                {pickl("geo"), pickl("dguid"), pickl("residence"), pickl("gender"), pickl("age"), pickl("group"),
                 dv(nu), dv(de), dv(es), dv(lo), dv(hi), svna(qu, qn2)});
  std::vector<std::string> cn, cdn;
  for (std::size_t j : a) {
    if (!in(ch[j], cn)) cn.push_back(ch[j]);
  }
  for (std::size_t j : b) {
    if (!in(ch[j], cdn)) cdn.push_back(ch[j]);
  }
  std::vector<std::string> both = cn, bn(cn.size(), "numerator");
  both.insert(both.end(), cdn.begin(), cdn.end());
  bn.insert(bn.end(), cdn.size(), "denominator");
  // one name per id: the profile's tree names each characteristic once
  out.at.push_back({"characteristics", named(sv(both), bn)});
  out.at.push_back({"conflvl", dv1(cf)});
  out.cls = {"profshr", "data.frame"};
  return out;
}

V f_cnsstd(const Cx &c) {
  const V &x = c["x"];
  if (!isprof(x)) return stop("`x` must come from cnsprof()");
  const std::string grp = !dflt(c["grp"]) ? str1(c["grp"], "grp") : "Indigenous identity";
  const std::string geo = !dflt(c["geo"]) ? str1(c["geo"], "geo") : "Canada";
  const V &ag = c["ages"];
  nofn(ag, "ages");
  if (bad()) return {};
  std::vector<std::string> ages = {"0 to 4 years",   "5 to 9 years",   "10 to 14 years", "15 to 19 years",
                                   "20 to 24 years", "25 to 34 years", "35 to 44 years", "45 to 54 years",
                                   "55 to 64 years", "65 to 74 years", "75 years and over"};
  if (!dflt(ag)) {
    if (!isstr(ag)) return stop("`ages` must name age groups");
    ages = ag.s;
  }
  const std::string resid = !dflt(c["resid"]) ? str1(c["resid"], "resid") : "Total - Residence on or off reserve";
  if (bad()) return {};
  const std::vector<int> chaid = ci(x, "chaid");
  if (bad()) return {};
  int first = INT_MAX;
  for (int v : chaid) {
    if (v != NAI) first = std::min(first, v);
  }
  const std::vector<std::string> gv = cl(x, "group"), geov = cl(x, "geo"), dg = cl(x, "dguid"), gen = cl(x, "gender"),
                                 res = cl(x, "residence"), age = cl(x, "age");
  if (!hascols(x)) return stop("`x` must come from cnsprof()");
  const Lab dgl = lab(*col(x, "dguid"));
  const std::vector<double> val = cd(x, "value");
  if (bad()) return {};
  std::vector<std::size_t> d;
  std::set<std::string> places;
  for (std::size_t j = 0; j < chaid.size(); ++j) {
    if (chaid[j] == first && gv[j] == grp && (geov[j] == geo || dg[j] == geo) &&
        (gen[j] == "Total - Gender" || gen[j] == "Total") && (res[j] == resid || res[j] == "Total")) {
      d.push_back(j);
      places.insert(dgl.na[j] ? std::string("\x01NA") : dg[j]);
    }
  }
  if (d.empty()) return stop("no population for group \"%s\" in \"%s\"", grp.c_str(), geo.c_str());
  if (places.size() > 1) {
    return stop("\"%s\" names %d places: give its DGUID as `geo`", geo.c_str(), static_cast<int>(places.size()));
  }
  std::vector<double> pop;
  std::vector<std::string> miss;
  for (const std::string &a : ages) {
    double v = na();
    for (std::size_t j : d) {
      if (age[j] == a) {
        v = val[j];
        break;
      }
    }
    if (std::isnan(v)) miss.push_back(a);
    pop.push_back(v);
  }
  if (!miss.empty()) return stop("no population for: %s", tostr(miss).c_str());
  const double s = rsum(pop);
  std::vector<double> w;
  for (double v : pop) w.push_back(v / s);
  return frame({"age", "population", "weight"}, {sv(ages), dv(pop), dv(w)});
}

V f_cnscwb(const Cx &c) {
  const V &x = c["x"], &dwell = c["dwell"];
  if (!isprof(x)) return stop("`x` must come from cnsprof()");
  if (dwell.t != V::NUL && !isprof(dwell)) return stop("`dwell` must come from cnsprof()");
  const V &gv = c["groups"];
  nofn(gv, "groups");
  if (bad()) return {};
  std::vector<std::string> groups = {"First Nations", "M\xc3\xa9tis", "Inuit", "Non-Indigenous"};
  if (!dflt(gv)) {
    bool ok = isstr(gv);
    for (const std::string &s : gv.s) ok = ok && !s.empty();
    if (!ok) return stop("`groups` must name one group or more");
    groups = gv.s;
  }
  const std::string resid = !dflt(c["resid"]) ? str1(c["resid"], "resid") : "Total - Residence on or off reserve";
  const double minpop = scalar(c["minpop"], "minpop", 0, HUGE_VAL, false, 250);
  if (bad()) return {};
  const std::vector<std::string> gen = cl(x, "gender"), res = cl(x, "residence"), age = cl(x, "age"),
                                 ch = cl(x, "characteristic"), dg = cl(x, "dguid"), gr = cl(x, "group"),
                                 geo = cl(x, "geo");
  const std::vector<int> chaid = ci(x, "chaid"), rootid = ci(x, "rootid");
  const std::vector<double> val = cd(x, "value");
  if (bad()) return {};
  const std::size_t n = gen.size();
  std::vector<char> tot(n);
  int pop_id = INT_MAX;
  for (std::size_t j = 0; j < n; ++j) {
    tot[j] = (gen[j] == "Total - Gender" || gen[j] == "Total") && (res[j] == resid || res[j] == "Total");
    if (chaid[j] != NAI) pop_id = std::min(pop_id, chaid[j]);
  }
  auto selt = [&](const std::string &name, const char *ag, const std::vector<int> *root) {
    return sel(ch, age, rootid, &tot, name, ag, root);
  };
  Rows edu, pop;
  std::vector<int> eid;
  for (std::size_t j = 0; j < n; ++j) {
    if (tot[j] && age[j] == "25 years and over" && ch[j].compare(0, 27, "Total - Highest certificate") == 0) {
      edu.push_back(j);
      if (std::find(eid.begin(), eid.end(), chaid[j]) == eid.end()) eid.push_back(chaid[j]);
    }
    if (tot[j] && chaid[j] == pop_id && age[j] == "Total - Age groups") pop.push_back(j);
  }
  const Rows rec = selt("Number of total income recipients aged 15 years and over in private households in 2020 - "
                       "25% sample data",
                       "Total - Age groups", nullptr);
  const Rows avg = selt("Average total income in 2020 among recipients ($)", "Total - Age groups", nullptr);
  const Rows none = selt("No certificate, diploma or degree", "25 years and over", &eid);
  const Rows ba = selt("Bachelor's degree or higher", "25 years and over", &eid);
  const Rows part = selt("Participation rate", "25 to 64 years", nullptr);
  const Rows unemp = selt("Unemployment rate", "25 to 64 years", nullptr);
  auto look = [&](const Rows &r, const std::vector<std::string> &d0, const std::vector<std::string> &g0,
                  const std::vector<double> &v0, const std::string &g, const std::string &grp) {
    for (std::size_t j : r) {
      if (d0[j] == g && g0[j] == grp) return v0[j];
    }
    return na();
  };
  std::vector<std::string> geos;
  for (const std::string &s : dg) {
    if (!in(s, geos)) geos.push_back(s);
  }
  std::vector<std::string> og, od, ogr;
  std::vector<double> inc, hs, uni, popn, pa, em;
  for (const std::string &g : geos) {
    for (const std::string &grp : groups) {
      og.push_back(geo[static_cast<std::size_t>(std::find(dg.begin(), dg.end(), g) - dg.begin())]);
      od.push_back(g);
      ogr.push_back(grp);
      inc.push_back(look(rec, dg, gr, val, g, grp) * look(avg, dg, gr, val, g, grp) / look(pop, dg, gr, val, g, grp));
      const double et = look(edu, dg, gr, val, g, grp);
      hs.push_back(100 * (1 - look(none, dg, gr, val, g, grp) / et));
      uni.push_back(100 * look(ba, dg, gr, val, g, grp) / et);
      popn.push_back(look(pop, dg, gr, val, g, grp));
      pa.push_back(look(part, dg, gr, val, g, grp));
      em.push_back(100 - look(unemp, dg, gr, val, g, grp));
    }
  }
  const std::size_t m = og.size();
  std::vector<double> notcro(m, na()), nomaj(m, na());
  if (dwell.t != V::NUL) {
    const std::vector<std::string> hch = cl(dwell, "characteristic"), hdg = cl(dwell, "dguid"),
                                   hgr = cl(dwell, "group");
    const std::vector<double> hv = cd(dwell, "value");
    if (bad()) return {};
    const std::vector<int> none0(hch.size(), 0);
    auto hh = [&](const std::string &name) { return sel(hch, hch, none0, nullptr, name, nullptr, nullptr); };
    const Rows ct = hh("Total - Private households by number of persons per room - 25% sample data");
    const Rows co = hh("One person or fewer per room");
    const Rows rt = hh("Total - Occupied private dwellings by dwelling condition - 25% sample data");
    const Rows rm = hh("Major repairs needed");
    for (std::size_t i = 0; i < m; ++i) {
      const std::string hg = ogr[i] == "Non-Indigenous" ? "Non-Indigenous households" : "Indigenous households";
      notcro[i] = 100 * look(co, hdg, hgr, hv, od[i], hg) / look(ct, hdg, hgr, hv, od[i], hg);
      nomaj[i] = 100 * (1 - look(rm, hdg, hgr, hv, od[i], hg) / look(rt, hdg, hgr, hv, od[i], hg));
    }
  }
  const double NA = na();
  std::vector<double> si(m, NA), se(m, NA), sh(m, NA), sl(m, NA), sc(m, NA);
  std::vector<std::size_t> ok;
  for (std::size_t i = 0; i < m; ++i) {
    if (std::isfinite(inc[i]) && inc[i] > 0 && std::isfinite(popn[i]) && popn[i] >= minpop) ok.push_back(i);
  }
  if (!ok.empty()) {
    auto pm = [](double v) { return std::isnan(v) ? v : std::min(100.0, v); };
    std::vector<double> a1, a2, a3, a4, a5, a6, a7;
    for (std::size_t i : ok) {
      a1.push_back(inc[i]);
      a2.push_back(pm(hs[i]));
      a3.push_back(pm(uni[i]));
      a4.push_back(std::isnan(notcro[i]) ? 50 : notcro[i]);
      a5.push_back(std::isnan(nomaj[i]) ? 50 : nomaj[i]);
      a6.push_back(pa[i]);
      a7.push_back(em[i]);
    }
    std::vector<std::string> an = {"incpercap", "highsch20p", "uni25p", "notcro", "nomajrep", "partic", "empl"};
    std::vector<V> av = {dv(a1), dv(a2), dv(a3), dv(a4), dv(a5), dv(a6), dv(a7)};
    if (c.has("bounds")) {
      an.push_back("bounds");
      av.push_back(c["bounds"]);
    }
    const Cx c2{&an, &av, c.warn};
    const V s = f_cwbscr(c2);
    if (bad()) return {};
    for (std::size_t k = 0; k < ok.size(); ++k) {
      const std::size_t i = ok[k];
      si[i] = s.a[0].d[k];
      se[i] = s.a[1].d[k];
      sh[i] = std::isnan(notcro[i]) || std::isnan(nomaj[i]) ? NA : s.a[2].d[k];
      sl[i] = s.a[3].d[k];
      const long double t = static_cast<long double>(si[i]) + se[i] + sh[i] + sl[i];
      sc[i] = static_cast<double>(t / 4);
    }
  }
  V out = frame({"geo", "dguid", "group", "income", "education", "housing", "labour", "cwb", "incpercap", "highsch25p",
                 "uni25p", "population", "participation", "employment", "nomajrep", "notcro"},
                {sv(og), sv(od), sv(ogr), dv(si), dv(se), dv(sh), dv(sl), dv(sc), dv(inc), dv(hs), dv(uni), dv(popn),
                 dv(pa), dv(em), dv(nomaj), dv(notcro)});
  out.cls = {"cnscwb", "data.frame"};
  return out;
}

namespace {

V trctab() {
  struct R {
    const char *call, *ind, *tab, *meas, *bet, *note;
  };
  static const R r[] = {
      {"19", "infant mortality", "17-10-0161", "Infant mortality rate", "lower",
       "by Indigenous identity, Canada; Inuit regions in 13-10-0100"},
      {"19", "maternal health", "13-10-0064", "Low birth weight; Pre-term births", "lower",
       "birth outcomes only, Canada and Inuit regions"},
      {"19", "suicide", "41-10-0021", "Seriously contemplated suicide in the past 12 months", "lower",
       "suicidal thoughts; deaths by cause in Inuit regions in 13-10-0157"},
      {"19", "mental health", "13-10-0924", "Perceived mental health, very good or excellent", "higher", ""},
      {"19", "mental health", "13-10-0924", "Mood disorder", "lower", ""},
      {"19", "mental health", "13-10-0924", "Anxiety disorder", "lower", ""},
      {"19", "mental health", "13-10-0924", "Perceived life stress, most days quite a bit or extremely stressful",
       "lower", ""},
      {"19", "addictions", "13-10-0924", "Heavy drinking", "lower", "a proxy"},
      {"19", "addictions", "13-10-0924", "Current smoker, daily or occasional", "lower", "a proxy"},
      {"19", "life expectancy", "17-10-0160", "Life expectancy at birth", "higher",
       "by Indigenous identity, Canada and regions; Inuit regions in 13-10-0403"},
      {"19", "birth rates", "census profile", "Population aged 0 to 4", "",
       "no direction: a measure of population, not of health"},
      {"19", "infant and child health issues", "13-10-0064", "Small for gestational age", "lower",
       "Canada and Inuit regions"},
      {"19", "chronic diseases", "13-10-0924", "At least one underlying chronic health condition", "lower", ""},
      {"19", "chronic diseases", "13-10-0924", "Diabetes", "lower", ""},
      {"19", "chronic diseases", "13-10-0924", "High blood pressure", "lower", ""},
      {"19", "chronic diseases", "13-10-0924", "Arthritis", "lower", ""},
      {"19", "illness and injury incidence", "13-10-0157", "Mortality", "lower",
       "deaths by cause, injuries among them, Canada and Inuit regions"},
      {"19", "availability of appropriate health services", "13-10-0924", "Has a regular healthcare provider", "higher",
       "a proxy: having a provider, not its cultural safety"},
      {"55", "children in care", "census profile", "Foster children", "lower",
       "foster children among children in census families (profshr)"},
      {"55", "educational attainment", "census profile", "No certificate, diploma or degree", "lower",
       "aged 25 and over"},
      {"55", "educational attainment", "census profile", "Bachelor's degree or higher", "higher", "aged 25 and over"},
      {"55", "income attainment", "census profile",
       "Prevalence of low income based on the Low-income measure, after tax (LIM-AT) (%)", "lower", ""},
      {"55", "income attainment", "census profile", "Average total income in 2020 among recipients ($)", "higher", ""}};
  std::vector<std::string> a, b, t, m, e, no;
  for (const R &x : r) {
    a.push_back(x.call);
    b.push_back(x.ind);
    t.push_back(x.tab);
    m.push_back(x.meas);
    e.push_back(x.bet);
    no.push_back(x.note);
  }
  return frame({"call", "indicator", "table", "measure", "better", "note"}, {sv(a), sv(b), sv(t), sv(m), sv(e), sv(no)});
}

// the last of sort(unique(p)): numbers by value, text byte by byte
V latest(const V &p, const std::vector<std::size_t> &r) {
  if (p.t == V::DBL || p.t == V::INT) {
    const std::vector<double> d = num(p, "period");
    if (bad()) return {};
    double best = -HUGE_VAL;
    bool any = false;
    for (std::size_t j : r) {
      if (!std::isnan(d[j])) {
        best = std::max(best, d[j]);
        any = true;
      }
    }
    if (!any) return V();
    V o = p.t == V::INT ? iv({static_cast<int>(best)}) : dv1(best);
    return o;
  }
  const Lab l = lab(p);
  std::string best;
  bool any = false;
  for (std::size_t j : r) {
    if (l.na[j]) continue;
    if (!any || l.s[j] > best) best = l.s[j];
    any = true;
  }
  return any ? sv1(best) : V();
}

bool same(const V &p, std::size_t j, const V &pr) {
  if (pr.t == V::NUL) return false;
  const Lab a = lab(rows(p, {j})), b = lab(pr);
  return !a.na[0] && a.s[0] == b.s[0];
}

}  // namespace

V f_trcind(const Cx &) { return trctab(); }

namespace {

// .trc_gap1(): one measure's gaps against the reference
bool gap1(const Cx &c, const V &d0, const std::vector<std::size_t> &r0, const std::string &measure, const V &period,
          const std::string &indicator, const std::string &better, const std::string &ref, std::vector<V> &rows1) {
  std::vector<std::size_t> r;
  const Lab g = lab(*col(d0, "group"));
  for (std::size_t j : r0) {
    if (!g.na[j]) r.push_back(j);
  }
  bool hasref = false;
  std::set<std::string> seen;
  bool dup = false;
  for (std::size_t j : r) {
    hasref = hasref || g.s[j] == ref;
    dup = dup || !seen.insert(g.s[j]).second;
  }
  if (!hasref) return false;
  const std::string pn = lab(period).s[0];
  if (dup) {
    c.w(sf("%s (%s): more than one row per group, skipped; filter the table to one row each", measure.c_str(),
           pn.c_str()));
    return false;
  }
  V d = rows(d0, r);
  d.df = true;
  std::vector<std::string> an = {"data", "grp", "est", "lower", "upper", "ref", "qual"};
  std::vector<V> av = {d, sv1("group"), sv1("estimate"), sv1("lower"), sv1("upper"), sv1(ref), sv1("quality")};
  std::vector<std::string> quiet;
  const Cx c2{&an, &av, &quiet};
  const V gp = f_indgap(c2);
  if (bad()) return false;
  const bool low = better == "lower";
  const std::size_t n = gp.a[0].len();
  std::vector<double> ra = gp.a[3].d, rl = gp.a[4].d, ru = gp.a[5].d, dr(n), lo(n), hi(n);
  for (std::size_t j = 0; j < n; ++j) {
    dr[j] = low ? ra[j] : 1 / ra[j];
    lo[j] = low ? rl[j] : 1 / ru[j];
    hi[j] = low ? ru[j] : 1 / rl[j];
  }
  V per = period;
  // assign(n, x[0]) reads x[0] after freeing it (libc++): a copy first
  if (per.t == V::INT) per.i = std::vector<int>(n, per.i[0]);
  else if (per.t == V::DBL) per.d = std::vector<double>(n, per.d[0]);
  else per = sv(std::vector<std::string>(n, pn));
  per.nm.clear();
  rows1.push_back(frame({"indicator", "measure", "period", "group", "estimate", "refest", "difference", "disratio",
                         "lower", "upper", "usewithcau"},
                        {sv(std::vector<std::string>(n, indicator)), sv(std::vector<std::string>(n, measure)), per,
                         gp.a[0], gp.a[1], gp.a[2], gp.a[6], dv(dr), dv(lo), dv(hi), gp.a[9]}));
  return true;
}

std::vector<std::size_t> allages(const V &d, const std::vector<std::string> &agetot, std::size_t which) {
  std::vector<std::size_t> r;
  const std::size_t n = nrow(d);
  const V *a = col(d, "age");
  const Lab l = a ? lab(*a) : Lab{};
  for (std::size_t j = 0; j < n; ++j) {
    if (!a || (!l.na[j] && l.s[j] == agetot[which])) r.push_back(j);
  }
  return r;
}

// not bind: a call on a std::vector also finds std::bind
V vbind(const std::vector<V> &p) {
  V o = p[0];
  for (std::size_t k = 1; k < p.size(); ++k) {
    for (std::size_t j = 0; j < o.a.size(); ++j) {
      V &a = o.a[j];
      const V &b = p[k].a[j];
      if (a.t != b.t) {
        const Lab la = lab(a), lb = lab(b);
        V s = svna(la.s, la.na);
        s.s.insert(s.s.end(), lb.s.begin(), lb.s.end());
        s.na.insert(s.na.end(), lb.na.begin(), lb.na.end());
        a = s;
        continue;
      }
      a.d.insert(a.d.end(), b.d.begin(), b.d.end());
      a.i.insert(a.i.end(), b.i.begin(), b.i.end());
      if (a.t == V::STR) {
        a.s.insert(a.s.end(), b.s.begin(), b.s.end());
        a.na.insert(a.na.end(), b.na.begin(), b.na.end());
      }
    }
  }
  o.nr = static_cast<int>(o.a[0].len());
  return o;
}

}  // namespace

V f_trchlthgap(const Cx &c) {
  const std::string ref = !dflt(c["ref"]) ? str1(c["ref"], "ref") : "Non-Indigenous";
  if (bad()) return {};
  std::vector<std::string> agetot = {"Total, 18 years and over", "Total, 15 years and over"};
  if (c.has("agetot")) {
    const V &a = c["agetot"];
    if (!isstr(a) || a.s.size() > 2) return stop("`agetot` must be one or two age labels (health, then suicide)");
    agetot = {a.s[0], a.s[a.s.size() - 1]};
  }
  static const char *const tn[] = {"health", "lifeexp", "infmort", "suiide"};
  for (const char *t : tn) {
    const V &v = c[t];
    if (v.t == V::NUL) continue;
    if (!isdf(v) || !col(v, "group") || !col(v, "estimate")) {
      return stop("`%s` must be a data frame of indicators (columns group, estimate, lower, upper, quality, and period, "
           "indicator or measure, age): the tables must be data frames or NULL",
           t);
    }
  }
  const V ti0 = trctab();
  std::vector<std::size_t> ti;
  for (std::size_t j = 0; j < ti0.a[0].len(); ++j) {
    if (ti0.a[0].s[j] == "19" && !ti0.a[4].s[j].empty()) ti.push_back(j);
  }
  auto tf = [&](int k, std::size_t j) { return ti0.a[static_cast<std::size_t>(k)].s[j]; };
  std::vector<V> pieces;
  const V &health = c["health"];
  if (health.t != V::NUL) {
    const std::vector<std::size_t> h = allages(health, agetot, 0);
    const V *pc = col(health, "period");
    const V *ic = col(health, "indicator");
    const V pr = pc ? latest(*pc, h) : V();
    if (bad()) return {};
    for (std::size_t j : ti) {
      if (tf(2, j) != "13-10-0924") continue;
      std::vector<std::size_t> r;
      const Lab il = ic ? lab(*ic) : Lab{};
      for (std::size_t k : h) {
        if (pc && same(*pc, k, pr) && ic && !il.na[k] && il.s[k] == tf(3, j)) r.push_back(k);
      }
      gap1(c, health, r, tf(3, j), pr.t == V::NUL ? sv1("NA") : pr, tf(1, j), tf(4, j), ref, pieces);
      if (bad()) return {};
    }
  }
  const char *sn[] = {"lifeexp", "infmort"};
  const char *st[] = {"17-10-0160", "17-10-0161"};
  for (int s = 0; s < 2; ++s) {
    const V &d = c[sn[s]];
    if (d.t == V::NUL) continue;
    std::size_t i = 0;
    for (std::size_t j : ti) {
      if (tf(2, j) == st[s]) i = j;
    }
    const std::vector<double> e = num(*col(d, "estimate"), "estimate");
    if (bad()) return {};
    std::vector<std::size_t> okr, all;
    for (std::size_t j = 0; j < e.size(); ++j) {
      all.push_back(j);
      if (std::isfinite(e[j])) okr.push_back(j);
    }
    const V *pc = col(d, "period");
    const V pr = pc ? latest(*pc, okr) : V();
    if (bad()) return {};
    std::vector<std::size_t> r;
    for (std::size_t j : all) {
      if (pc && same(*pc, j, pr)) r.push_back(j);
    }
    gap1(c, d, r, tf(3, i), pr.t == V::NUL ? sv1("NA") : pr, tf(1, i), tf(4, i), ref, pieces);
    if (bad()) return {};
  }
  const V &su = c["suiide"];
  if (su.t != V::NUL) {
    std::size_t i = 0;
    for (std::size_t j : ti) {
      if (tf(2, j) == "41-10-0021") i = j;
    }
    std::vector<std::size_t> r0 = allages(su, agetot, 1), r;
    const V *mc = col(su, "measure");
    const Lab ml = mc ? lab(*mc) : Lab{};
    for (std::size_t j : r0) {
      if (mc && !ml.na[j] && ml.s[j] == tf(3, i)) r.push_back(j);
    }
    V pr = sv1("2017");
    const V *pc = col(su, "period");
    if (pc) {
      const V lp = latest(*pc, r);
      if (bad()) return {};
      if (lp.t != V::NUL) {
        pr = lp;
        std::vector<std::size_t> r2;
        for (std::size_t j : r) {
          if (same(*pc, j, pr)) r2.push_back(j);
        }
        r = r2;
      }
    }
    const Lab gl = lab(*col(su, "group"));
    bool hasref = false;
    for (std::size_t j : r) hasref = hasref || (!gl.na[j] && gl.s[j] == ref);
    if (hasref) gap1(c, su, r, tf(3, i), pr, tf(1, i), tf(4, i), ref, pieces);
    if (bad()) return {};
  }
  if (pieces.empty()) return stop("no table given measures a Call 19 indicator against \"%s\"", ref.c_str());
  V out = vbind(pieces);
  const std::vector<double> &lo = out.a[8].d;
  std::vector<int> worse;
  for (double v : lo) worse.push_back(std::isfinite(v) && v > 1);
  out.nm.push_back("worse");
  out.a.push_back(lv(worse));
  out.cls = {"trchlthgap", "data.frame"};
  return out;
}

V f_polex(const Cx &c) {
  std::uint32_t sd = 20240101;
  seedarg(c["seed"], sd);
  if (bad()) return {};
  Rng g(sd);
  const int nd = 20, nm = 36;
  std::vector<double> pop(nd), a(nd), m(nm);
  for (double &v : pop) v = rround(40000 + (250000 - 40000) * g.unif(), 0);
  for (double &v : a) v = 30 + 5 * g.norm();
  for (double &v : m) v = 1.5 * g.norm();
  const double stage[4] = {13, 19, 25, na()};
  const std::size_t n = static_cast<std::size_t>(nd * nm);
  std::vector<std::string> dist, reg;
  std::vector<int> month, crimes;
  std::vector<double> refm, rf, pp, lp, lam, cr, ar;
  for (int i = 1; i <= nd; ++i) {
    for (int t = 1; t <= nm; ++t) {
      char b[16];
      std::snprintf(b, sizeof b, "D%02d", i);
      dist.push_back(b);
      reg.push_back("R" + std::to_string((i - 1) % 4 + 1));
      month.push_back(t);
      const double s = stage[(i - 1) % 4];
      refm.push_back(s);
      const double r = !std::isnan(s) && t >= s ? 1 : 0;
      rf.push_back(r);
      pp.push_back(pop[static_cast<std::size_t>(i - 1)]);
      lp.push_back(std::log(pop[static_cast<std::size_t>(i - 1)]));
      lam.push_back(std::max(1.0, a[static_cast<std::size_t>(i - 1)] + m[static_cast<std::size_t>(t - 1)] + 0.05 * t - 3 * r));
    }
  }
  for (std::size_t j = 0; j < n; ++j) crimes.push_back(static_cast<int>(g.pois(lam[j] * pp[j] / 1e5)));
  for (std::size_t j = 0; j < n; ++j) cr.push_back(rround(lam[j] + g.norm(), 3));
  for (std::size_t j = 0; j < n; ++j) ar.push_back(rround(0.4 * lam[j] + 1.2 * rf[j] + 0.8 * g.norm(), 3));
  return frame({"district", "region", "month", "refmonth", "reformed", "population", "logpop", "crimes", "crmrate",
                "arrrate"},
               {sv(dist), sv(reg), iv(month), dv(refm), dv(rf), dv(pp), dv(lp), iv(crimes), dv(cr), dv(ar)});
}

const Fn fns[] = {
    {"sctidy", "x"},
    {"qltyrcd", "x suppr? caution?"},
    {"limfrom", "est se? cv? moe? moelev? conflvl?"},
    {"indgap", "data grp est lower upper ref? qual? conflvl? pooled?"},
    {"genexpo", "parent gparent partial?"},
    {"exptrnd", "cases n dose? conflvl?"},
    {"trccmpl", ""},
    {"trcprog", "done? yearend? total? start? conflvl?"},
    {"mhcont", "data grp ind est lower upper hope? belong? meaning? purp? qual? conflvl?"},
    {"cwbscr", "incpercap highsch20p uni25p notcro nomajrep partic empl bounds?"},
    {"connassoc", "data place x y xlo? xup? ylo? yup? reps? seed? conflvl?"},
    {"healchg", "data grp ind prd est lower upper from to well? qual? conflvl? distress?"},
    {"commidx", "data grp ind est lower upper w qual? conflvl?"},
    {"frmaud", "x"},
    {"scests", "x"},
    {"ordcmp", "data grp cat est lower upper lvls ref qual? conflvl? pooled?"},
    {"stdprev", "est lower upper std grp? qual? conflvl? scale?"},
    {"welprof", "data grp ind est lower upper well distress domain? qual?"},
    {"holbal", "data grp domain est lower upper w? qual? reps? seed? conflvl?"},
    {"dualcont", "x well distress conflvl?"},
    {"cnsprof", "x"},
    {"profshr", "x num den? smpfrac? deff? conflvl? roubase?"},
    {"cnsstd", "x grp? geo? ages? resid?"},
    {"cnscwb", "x dwell? groups? bounds? resid? minpop?"},
    {"trcind", ""},
    {"trchlthgap", "health? lifeexp? infmort? suiide? ref? agetot?"},
    {"polex", "seed?"},
};
const int nfns = static_cast<int>(sizeof fns / sizeof fns[0]);

int call(const char *fn, const std::vector<std::string> &an, const std::vector<V> &av, V &out, std::string &err,
         std::vector<std::string> &warn) {
  typedef V (*F)(const Cx &);
  static const F fp[] = {f_sctidy,  f_qltyrcd, f_limfrom, f_indgap,  f_genexpo,  f_exptrnd,  f_trccmpl,
                         f_trcprog, f_mhcont,  f_cwbscr,  f_connassoc, f_healchg, f_commidx,  f_frmaud,
                         f_scests,  f_ordcmp,  f_stdprev, f_welprof, f_holbal,   f_dualcont, f_cnsprof,
                         f_profshr, f_cnsstd,  f_cnscwb,  f_trcind,  f_trchlthgap, f_polex};
  static_assert(sizeof fp / sizeof fp[0] == sizeof fns / sizeof fns[0], "one function per name");
  const Cx c{&an, &av, &warn};
  for (int k = 0; k < nfns; ++k) {
    if (std::string(fns[k].name) != fn) continue;
    // the refusal of this call (see stop()); the one of a call further up
    // the stack, if any, comes back after it
    Stop st;
    Stop *const up = setst(&st);
    // .rmbl_no_fn(): no argument may be a function or an environment
    for (std::size_t j = 0; j < av.size() && !st.on; ++j) nofn(av[j], an[j].c_str());
    V r;
    if (!st.on) r = fp[k](c);
    setst(up);
    if (st.on) {
      err = st.msg;
      return CLAB_EARG;
    }
    // every warning starts with the function's name, as the refusals do
    const std::string pre = std::string(fn) + ": ";
    for (std::string &w : warn) {
      if (w.compare(0, pre.size(), pre) != 0) w = pre + w;
    }
    out = std::move(r);
    return CLAB_OK;
  }
  err = std::string("no function ") + fn;
  return CLAB_EARG;
}

}  // namespace igh
