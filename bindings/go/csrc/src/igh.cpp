// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Statistics Canada's published tables read as published (quality flags
// kept, suppressed cells never imputed), compared group by group, the
// generations of residential school attendance and their dose-response,
// and the pace of the Truth and Reconciliation Commission's Calls to
// Action: sctidy(), qltyrcd(), limfrom(), indgap(), genexpo(), exptrnd(),
// trccmpl() and trcprog() (igh.h).
#include <algorithm>
#include <climits>
#include <cmath>
#include <map>
#include <set>

#include "crimlab.h"
#include "igk.h"

namespace igh {

namespace {

const std::vector<std::string> qlev = {"ok", "E", "F", "x", "..", "..."};

// A name ending in "[n]" after its colon: a column of one member of a
// census table's dimension.
bool membercol(const std::string &nm) {
  const std::size_t c = nm.find(':');
  if (c == std::string::npos) return false;
  std::size_t e = nm.size();
  while (e > c + 1 && (nm[e - 1] == ' ' || nm[e - 1] == '\t' || nm[e - 1] == '\n' || nm[e - 1] == '\r')) --e;
  if (e == 0 || nm[e - 1] != ']') return false;
  std::size_t b = e - 1;
  while (b > c + 1 && std::isdigit(static_cast<unsigned char>(nm[b - 1]))) --b;
  return b < e - 1 && b > c + 1 && nm[b - 1] == '[';
}

// the member label: after the first colon, less its "[n]" and the space
// before it
std::string memberlab(const std::string &nm) {
  std::string s = nm.substr(nm.find(':') + 1);
  std::size_t e = s.size();
  while (e > 0 && std::isspace(static_cast<unsigned char>(s[e - 1]))) --e;
  std::size_t b = e - 1;
  while (b > 0 && std::isdigit(static_cast<unsigned char>(s[b - 1]))) --b;
  --b;
  while (b > 0 && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
  return trim(s.substr(0, b));
}

void setcol(V &df, const std::string &nm, V v) {
  const int j = coli(df, nm);
  if (j >= 0) {
    df.a[static_cast<std::size_t>(j)] = std::move(v);
  } else {
    df.nm.push_back(nm);
    df.a.push_back(std::move(v));
  }
}

// .rmbl_statcan_long(): a census table's wide layout made long, or false
bool sclong(V &x) {
  const std::vector<std::string> &nm = x.nm;
  std::vector<std::size_t> vc;
  for (std::size_t j = 0; j < nm.size(); ++j) {
    if (membercol(nm[j])) vc.push_back(j);
  }
  if (vc.empty()) return false;
  std::vector<std::string> dims;
  for (std::size_t j : vc) dims.push_back(trim(nm[j].substr(0, nm[j].find(':'))));
  dims = uniq(dims);
  if (dims.size() != 1) stop("`x` spreads more than one dimension across its columns");
  std::vector<char> hs(vc.size()), drop(nm.size(), 0);
  for (std::size_t k = 0; k < vc.size(); ++k) {
    const std::size_t sj = vc[k] + 1;
    hs[k] = sj < nm.size() && nm[sj].compare(0, 6, "Symbol") == 0;
    drop[vc[k]] = 1;
    if (hs[k]) drop[sj] = 1;
  }
  const std::size_t n = nrow(x);
  std::vector<std::size_t> rr;
  for (std::size_t k = 0; k < vc.size(); ++k) {
    for (std::size_t r = 0; r < n; ++r) rr.push_back(r);
  }
  V out;
  out.t = V::LIST;
  out.df = true;
  out.cls = {"data.frame"};
  for (std::size_t j = 0; j < nm.size(); ++j) {
    if (drop[j]) continue;
    out.nm.push_back(nm[j]);
    out.a.push_back(rows(x.a[j], rr));
  }
  std::vector<std::string> mem, val, st;
  std::vector<char> vna, sna;
  for (std::size_t k = 0; k < vc.size(); ++k) {
    const std::string m = memberlab(nm[vc[k]]);
    const Lab lv = lab(x.a[vc[k]]);
    for (std::size_t r = 0; r < n; ++r) {
      mem.push_back(m);
      val.push_back(lv.s[r]);
      vna.push_back(lv.na[r]);
    }
    if (hs[k]) {
      const Lab ls = lab(x.a[vc[k] + 1]);
      for (std::size_t r = 0; r < n; ++r) {
        st.push_back(ls.s[r]);
        sna.push_back(ls.na[r]);
      }
    } else {
      for (std::size_t r = 0; r < n; ++r) {
        st.push_back("");
        sna.push_back(0);
      }
    }
  }
  setcol(out, dims[0], sv(mem));
  setcol(out, "VALUE", svna(val, vna));
  setcol(out, "STATUS", svna(st, sna));
  out.nr = static_cast<int>(rr.size());
  x = out;
  return true;
}

// .rmbl_statcan_members(): a label that names two members of a dimension
// gets the member number of its coordinate
void scmembers(V &x) {
  int cc = coli(x, "COORDINATE");
  if (cc < 0) cc = coli(x, "Coordinate");
  const int gi = coli(x, "GEO");
  if (cc < 0 || gi < 0) return;
  static const std::set<std::string> meta = {"DGUID", "UOM", "UOM_ID", "SCALAR_FACTOR", "SCALAR_ID", "VECTOR"};
  std::vector<int> dims;
  // seq(match("GEO"), match(cc) - 1), downwards when the coordinate comes first
  const int to = cc - 1;
  for (int j = gi;; j += (to >= gi ? 1 : -1)) {
    if (!meta.count(x.nm[static_cast<std::size_t>(j)])) dims.push_back(j);
    if (j == to) break;
  }
  const std::size_t k = dims.size();
  if (k < 2) return;
  const Lab co = lab(x.a[static_cast<std::size_t>(cc)]);
  const std::size_t n = co.n();
  std::vector<std::vector<std::string>> ids(n);
  for (std::size_t r = 0; r < n; ++r) {
    if (co.na[r]) return;
    std::string p;
    for (char ch : co.s[r]) {
      if (ch == '.') {
        ids[r].push_back(p);
        p.clear();
      } else {
        p += ch;
      }
    }
    // strsplit() drops an empty last piece
    if (!p.empty()) ids[r].push_back(p);
    if (ids[r].size() < k) return;
  }
  for (std::size_t j = 1; j < k; ++j) {
    V &cv = x.a[static_cast<std::size_t>(dims[j])];
    const Lab lb = lab(cv);
    std::map<std::string, std::set<std::string>> seen;
    for (std::size_t r = 0; r < n; ++r) {
      if (!lb.na[r]) seen[lb.s[r]].insert(ids[r][j]);
    }
    bool any = false;
    std::vector<std::string> s = lb.s;
    for (std::size_t r = 0; r < n; ++r) {
      if (!lb.na[r] && seen[lb.s[r]].size() > 1) {
        s[r] = lb.s[r] + " [" + ids[r][j] + "]";
        any = true;
      }
    }
    if (any) cv = svna(s, lb.na);
  }
}

}  // namespace

V f_sctidy(const Cx &c) {
  V x = c["x"];
  nofn(x, "x");
  if (!isdf(x)) stop("`x` must be a data frame read from a Statistics Canada table");
  if (coli(x, "VALUE") < 0 && !sclong(x)) {
    stop("`x` has no VALUE column, nor census columns: it is not a Statistics Canada table");
  }
  scmembers(x);
  const std::size_t n = nrow(x);
  Lab st;
  if (const V *s = col(x, "STATUS")) {
    st = lab(*s);
  } else {
    st.s.assign(n, "");
    st.na.assign(n, 0);
  }
  std::vector<std::string> q(n);
  for (std::size_t r = 0; r < n; ++r) {
    const std::string t = st.na[r] ? "" : trim(st.s[r]);
    q[r] = t.empty() ? "ok" : t;
    if (std::find(qlev.begin(), qlev.end(), q[r]) == qlev.end()) q[r] = "ok";
  }
  const Lab vl = lab(*col(x, "VALUE"));
  std::vector<double> v(n);
  for (std::size_t r = 0; r < n; ++r) {
    V one = vl.na[r] ? svna({""}, {1}) : sv1(vl.s[r]);
    double y = na();
    if (!vl.na[r]) {
      try {
        y = num(one, "VALUE")[0];
      } catch (const Stop &) {
        y = na();
      }
    }
    v[r] = gone(q[r]) ? na() : y;
  }
  setcol(x, "VALUE", dv(v));
  setcol(x, "quality", fac(q, qlev));
  x.df = true;
  x.nr = static_cast<int>(n);
  x.rn.clear();
  x.cls = {"sctidy", "data.frame"};
  x.js = false;
  return x;
}

V f_qltyrcd(const Cx &c) {
  const V &x = c["x"], &sp = c["suppr"], &ca = c["caution"];
  nofn(x, "x");
  nofn(sp, "suppr");
  nofn(ca, "caution");
  if (x.t == V::NUL || x.t == V::LIST || x.odd) stop("`x` must be a vector of flags");
  bool allna = true;
  for (std::size_t j = 0; j < x.len(); ++j) {
    const Lab l = lab(rows(x, {j}));
    allna = allna && l.na[0];
  }
  if (x.t != V::STR && x.t != V::FAC && !allna) {
    const char *k = x.t == V::DBL ? "numeric" : x.t == V::INT ? "integer" : "logical";
    stop("`x` must be a vector of flags (character), not %s", k);
  }
  const V spv = dflt(sp) ? sv({"F", "x", "..", "..."}) : sp;
  const V cav = dflt(ca) ? sv1("E") : ca;
  for (int k = 0; k < 2; ++k) {
    const V &v = k ? cav : spv;
    bool ok = v.t == V::STR;
    for (char m : v.na) ok = ok && !m;
    if (!ok) stop("`%s` must be flags (character)", k ? "caution" : "suppr");
  }
  for (const std::string &a : spv.s) {
    if (std::find(cav.s.begin(), cav.s.end(), a) != cav.s.end()) {
      stop("a flag cannot be both suppressed and caution");
    }
  }
  std::set<std::string> ts, tc;
  for (const std::string &a : spv.s) ts.insert(trim(a));
  for (const std::string &a : cav.s) tc.insert(trim(a));
  const Lab f = lab(x);
  std::vector<std::string> out(f.n(), "ok");
  std::vector<std::string> odd;
  for (std::size_t j = 0; j < f.n(); ++j) {
    if (f.na[j]) continue;
    const std::string t = trim(f.s[j]);
    if (tc.count(t)) out[j] = "E";
    if (ts.count(t)) out[j] = "x";
    if (!t.empty() && !ts.count(t) && !tc.count(t) && t != "ok" &&
        std::find(odd.begin(), odd.end(), t) == odd.end()) {
      odd.push_back(t);
    }
  }
  if (!odd.empty()) {
    c.w(sf("qltyrcd: flag(s) %s are in neither `suppr` nor `caution`: treated as ok", quoted(odd).c_str()));
  }
  return sv(out);
}

V f_limfrom(const Cx &c) {
  const std::vector<double> est = num(c["est"], "est");
  const bool g[3] = {c.has("se"), c.has("cv"), c.has("moe")};
  if (g[0] + g[1] + g[2] != 1) stop("give one of `se`, `cv` or `moe`");
  const double cl = conf(c["conflvl"], "conflvl", 0.95);
  const double ml = conf(c["moelev"], "moelev", 0.90);
  std::vector<double> s;
  if (g[0]) {
    s = num(c["se"], "se");
  } else if (g[1]) {
    const std::vector<double> cv = num(c["cv"], "cv");
    const std::size_t n = est.empty() || cv.empty() ? 0 : std::max(est.size(), cv.size());
    for (std::size_t j = 0; j < n; ++j) s.push_back(est[j % est.size()] * cv[j % cv.size()] / 100);
  } else {
    const std::vector<double> m = num(c["moe"], "moe");
    const double zm = qn(ml);
    for (double x : m) s.push_back(x / zm);
  }
  if (s.size() != est.size() && s.size() != 1) stop("the standard errors must be one per estimate");
  for (double x : s) {
    if (x < 0) stop("standard errors, CVs and margins must not be negative");
  }
  const double z = qn(cl);
  std::vector<double> lo, hi, se;
  for (std::size_t j = 0; j < est.size(); ++j) {
    const double sj = s[j % s.size()];
    lo.push_back(est[j] - z * sj);
    hi.push_back(est[j] + z * sj);
    se.push_back(sj);
  }
  V out = frame({"estimate", "lower", "upper", "se"}, {dv(est), dv(lo), dv(hi), dv(se)});
  out.cls = {"limfrom", "data.frame"};
  return out;
}

namespace {

// the column arguments of indgap(), ordcmp(): their names, in data
std::string gapcol(const V &data, const V &v, const char *what, const char *fn) {
  const std::string s = str1(v, what);
  (void)fn;
  if (!col(data, s)) stop("column `%s` not found in `data`", s.c_str());
  return s;
}

}  // namespace

// One comparison of published estimates with a reference: the frame of
// indgap(), also used by trchlthgap().
V gapframe(const Cx &c, const Lab &g, std::vector<double> est, const std::vector<double> &lo,
           const std::vector<double> &hi, const std::vector<std::string> &q, const std::string &ref,
           const std::string &grp, double cl, const V &pooled) {
  std::vector<std::string> gs;
  for (std::size_t j = 0; j < g.n(); ++j) gs.push_back(g.na[j] ? std::string("\x01NA") : g.s[j]);
  {
    std::set<std::string> u(gs.begin(), gs.end());
    if (u.size() != gs.size()) {
      stop("a group appears twice: filter `data` to one period, age group, sex and place");
    }
  }
  const std::string tr = trim(ref);
  std::vector<std::size_t> rr;
  for (std::size_t j = 0; j < g.n(); ++j) {
    if (!g.na[j] && g.s[j] == tr) rr.push_back(j);
  }
  if (rr.size() != 1) {
    stop("`ref` must name exactly one row of `%s` (\"%s\" names %d)", grp.c_str(), ref.c_str(),
         static_cast<int>(rr.size()));
  }
  const std::size_t r = rr[0];
  if (gone(q[r]) || !std::isfinite(est[r])) {
    stop("the reference estimate is suppressed or missing: there is nothing to compare with");
  }
  std::vector<char> bad(est.size());
  for (std::size_t j = 0; j < est.size(); ++j) {
    bad[j] = std::isfinite(lo[j]) && std::isfinite(hi[j]) && std::isfinite(est[j]) && (lo[j] > est[j] || hi[j] < est[j]);
  }
  if (bad[r]) stop("the reference group's published interval does not contain its estimate");
  std::vector<std::string> bl;
  for (std::size_t j = 0; j < est.size(); ++j) {
    if (bad[j]) {
      bl.push_back(g.na[j] ? "NA" : g.s[j]);
      est[j] = na();
    }
  }
  if (!bl.empty()) {
    c.w(sf("indgap: the interval does not contain the estimate for: %s; not compared", tostr(bl).c_str()));
  }
  const double z = qn(cl);
  const double NA = na();
  auto sel = [&](std::size_t j) {
    if (std::isnan(lo[j]) || std::isnan(hi[j])) return NA;
    return lo[j] > 0 && hi[j] > 0 ? (std::log(hi[j]) - std::log(lo[j])) / (2 * z) : NA;
  };
  std::vector<std::string> grs;
  std::vector<char> gna;
  std::vector<double> e, re, ra, rl, ru, di, dl, du;
  std::vector<int> uc, pg;
  std::vector<std::string> gu;
  for (std::size_t j = 0; j < g.n(); ++j) {
    if (j == r) continue;
    grs.push_back(g.s[j]);
    gna.push_back(g.na[j]);
    gu.push_back(g.na[j] ? "" : g.s[j]);
    const bool ok = !gone(q[j]) && std::isfinite(est[j]);
    const double dif = est[j] - est[r];
    const double sej = (hi[j] - lo[j]) / (2 * z), ser = (hi[r] - lo[r]) / (2 * z);
    const double sd = std::sqrt(sej * sej + ser * ser);
    const double ratio = est[j] / est[r];
    const double a = sel(j), b = sel(r);
    const double sr = std::sqrt(a * a + b * b);
    e.push_back(est[j]);
    re.push_back(est[r]);
    ra.push_back(ok ? ratio : NA);
    const bool rp = ok && !std::isnan(ratio) && ratio > 0;
    rl.push_back(rp ? std::exp(std::log(ratio) - z * sr) : NA);
    ru.push_back(rp ? std::exp(std::log(ratio) + z * sr) : NA);
    di.push_back(ok ? dif : NA);
    dl.push_back(ok ? dif - z * sd : NA);
    du.push_back(ok ? dif + z * sd : NA);
    uc.push_back(q[j] == "E" || q[r] == "E");
  }
  pg = pooledrows(gu, pooled);
  for (std::size_t j = 0; j < gna.size(); ++j) {
    if (gna[j] && pooled.t == V::NUL) pg[j] = 0;
    if (gna[j] && pooled.t != V::NUL) pg[j] = 0;
  }
  V out = frame({"group", "estimate", "refest", "ratio", "ratiolo", "ratioup", "difference", "diflo", "difup",
                 "usewithcau", "poogrp"},
                {svna(grs, gna), dv(e), dv(re), dv(ra), dv(rl), dv(ru), dv(di), dv(dl), dv(du), lv(uc), lv(pg)});
  out.at.push_back({"reference", sv1(ref)});
  out.at.push_back({"conflvl", dv1(cl)});
  out.cls = {"indgap", "data.frame"};
  return out;
}

V f_indgap(const Cx &c) {
  const V &data = c["data"];
  if (!isdf(data)) stop("`data` must be a data frame");
  const std::string grp = gapcol(data, c["grp"], "grp", "indgap");
  const std::string est = gapcol(data, c["est"], "est", "indgap");
  const std::string lwc = gapcol(data, c["lower"], "lower", "indgap");
  const std::string upc = gapcol(data, c["upper"], "upper", "indgap");
  std::string qual;
  if (c.has("qual")) qual = gapcol(data, c["qual"], "qual", "indgap");
  const std::string ref = !dflt(c["ref"]) ? str1(c["ref"], "ref") : "Non-Indigenous";
  const V &pooled = c["pooled"];
  {
    const std::string tr = trim(ref);
    const bool tot = lower(tr).compare(0, 5, "total") == 0;
    if (tot || pooledrows({ref}, pooled)[0]) {
      c.w("the reference looks like a total that contains the groups: the intervals assume independent "
          "samples and are not right for it; compare with the groups' complement (for example "
          "non-Indigenous people)");
    }
  }
  const double cl = conf(c["conflvl"], "conflvl", 0.95);
  const Lab g = tlab(*col(data, grp));
  {
    std::vector<std::string> gs;
    for (std::size_t j = 0; j < g.n(); ++j) gs.push_back(g.na[j] ? "\x01" : g.s[j]);
    pooledcheck(pooled, gs);
  }
  const std::vector<double> e = num(*col(data, est), "est");
  const std::vector<double> lo = num(*col(data, lwc), "lower");
  const std::vector<double> hi = num(*col(data, upc), "upper");
  const std::vector<std::string> q = quality(c, data, qual, g.n());
  return gapframe(c, g, e, lo, hi, q, ref, grp, cl, pooled);
}

V f_genexpo(const Cx &c) {
  const std::string partial = pick(c["partial"], "partial", {"na", "lowerbnd"});
  const V &p = c["parent"], &gp = c["gparent"];
  for (int k = 0; k < 2; ++k) {
    const V &v = k ? gp : p;
    if (v.t != V::LGL) stop("`%s` must be TRUE, FALSE or NA for each person", k ? "gparent" : "parent");
  }
  if (p.len() != gp.len()) stop("`parent` and `gparent` must be the same length");
  std::vector<int> out(p.len());
  for (std::size_t j = 0; j < out.size(); ++j) {
    const int a = p.i[j], b = gp.i[j];
    int o = (a == 1) + (b == 1);
    const bool en = a == INT_MIN || b == INT_MIN;
    if (en && o == 0) o = INT_MIN;
    if (partial == "na" && en && o == 1) o = INT_MIN;
    out[j] = o;
  }
  return iv(out);
}

V f_exptrnd(const Cx &c) {
  const std::vector<double> a = num(c["cases"], "cases");
  const std::vector<double> n = num(c["n"], "n");
  std::vector<double> dose;
  if (!dflt(c["dose"])) {
    dose = num(c["dose"], "dose");
  } else {
    for (std::size_t j = 0; j < a.size(); ++j) dose.push_back(static_cast<double>(j));
  }
  const double cl = conf(c["conflvl"], "conflvl", 0.95);
  const std::size_t k = a.size();
  if (k < 2 || n.size() != k || dose.size() != k) {
    stop("`cases`, `n` and `dose` must have one entry per level, two levels at least");
  }
  bool bad = false, a0 = true, an = true;
  for (std::size_t j = 0; j < k; ++j) {
    bad = bad || std::isnan(a[j]) || std::isnan(n[j]) || std::isnan(dose[j]) || a[j] < 0 || n[j] <= 0 ||
          a[j] > n[j] || std::fabs(a[j] - std::round(a[j])) > 1e-8 || std::fabs(n[j] - std::round(n[j])) > 1e-8;
  }
  if (bad) stop("`cases` and `n` must be whole numbers with 0 <= cases <= n and n > 0");
  for (std::size_t j = 1; j < k; ++j) {
    if (!(dose[j] > dose[j - 1])) stop("`dose` must be strictly increasing");
  }
  for (std::size_t j = 0; j < k; ++j) {
    a0 = a0 && a[j] == 0;
    an = an && a[j] == n[j];
  }
  if (a0 || an) {
    stop("no level has a case (or every person at every level is one): there is no trend to test");
  }
  const double z = qn(cl);
  double adj = 0;
  for (std::size_t j = 0; j < k; ++j) {
    if (a[j] == 0 || n[j] - a[j] == 0) adj = 0.5;
  }
  std::vector<double> pr, orr, lo, hi;
  const double b1 = n[0] - a[0];
  for (std::size_t j = 0; j < k; ++j) {
    const double b = n[j] - a[j];
    const double o = ((a[j] + adj) / (b + adj)) / ((a[0] + adj) / (b1 + adj));
    const double se = j == 0 ? 0 : std::sqrt(1 / (a[j] + adj) + 1 / (b + adj) + 1 / (a[0] + adj) + 1 / (b1 + adj));
    pr.push_back(a[j] / n[j]);
    orr.push_back(o);
    lo.push_back(j == 0 ? 1 : std::exp(std::log(o) - z * se));
    hi.push_back(j == 0 ? 1 : std::exp(std::log(o) + z * se));
  }
  V out = frame({"dose", "cases", "n", "proportion", "orratio", "lower", "upper"},
                {dv(dose), dv(a), dv(n), dv(pr), dv(orr), dv(lo), dv(hi)});
  // prop.trend.test(): the score's sum of squares in the weighted fit of
  // the proportions on the scores
  double sa = 0, sn = 0;
  for (std::size_t j = 0; j < k; ++j) {
    sa += a[j];
    sn += n[j];
  }
  const double p = sa / sn;
  double sw = 0, ss = 0, sf = 0;
  std::vector<double> w(k), f(k);
  for (std::size_t j = 0; j < k; ++j) {
    w[j] = n[j] / p / (1 - p);
    f[j] = a[j] / n[j];
    sw += w[j];
    ss += w[j] * dose[j];
    sf += w[j] * f[j];
  }
  const double mb = ss / sw, fb = sf / sw;
  double sxy = 0, sxx = 0;
  for (std::size_t j = 0; j < k; ++j) {
    sxy += w[j] * (dose[j] - mb) * (f[j] - fb);
    sxx += w[j] * (dose[j] - mb) * (dose[j] - mb);
  }
  const double chi = sxy * sxy / sxx;
  std::string dn = "a out of n ,\n using scores:";
  for (double d : dose) dn += " " + rch(d);
  V tr = nlist({"statistic", "parameter", "p.value", "method", "data.name"},
               {named(dv1(chi), {"X-squared"}), named(dv1(1), {"df"}), dv1(clab_pchisq(chi, 1, 0)),
                sv1("Chi-squared Test for Trend in Proportions"), sv1(dn)});
  tr.cls = {"htest"};
  out.at.push_back({"trend", tr});
  out.cls = {"exptrnd", "data.frame"};
  return out;
}

V f_trccmpl(const Cx &) {
  return frame({"period", "yearend", "completed", "cumulative", "source"},
               {sv({"2015-2019", "2020", "2021", "2022"}), iv({2019, 2020, 2021, 2022}), iv({8, 0, 3, 2}),
                iv({8, 8, 11, 13}),
                sv({"Jewell and Mosby (2022): 13 complete in all", "Jewell and Mosby (2021)",
                    "Jewell and Mosby (2021)", "Jewell and Mosby (2022)"})});
}

V f_trcprog(const Cx &c) {
  const std::vector<double> done = !dflt(c["done"]) ? num(c["done"], "done", true) : std::vector<double>{8, 8, 11, 13};
  const std::vector<double> ye = !dflt(c["yearend"]) ? num(c["yearend"], "yearend")
                                                  : std::vector<double>{2019, 2020, 2021, 2022};
  const double total = scalar(c["total"], "total", 1, 1e6, true, 94);
  const double start = scalar(c["start"], "start", -HUGE_VAL, HUGE_VAL, false, 2015.96);
  const double cl = conf(c["conflvl"], "conflvl", 0.95);
  bool nas = false;
  for (double x : done) nas = nas || std::isnan(x);
  for (double x : ye) nas = nas || std::isnan(x);
  if (done.empty() || done.size() != ye.size() || nas) {
    stop("`done` and `yearend` must be one count per year, with no NA");
  }
  for (std::size_t j = 1; j < ye.size(); ++j) {
    if (!(ye[j] > ye[j - 1])) stop("`yearend` must be increasing");
  }
  for (std::size_t j = 0; j < done.size(); ++j) {
    if ((j && done[j] < done[j - 1]) || done[j] < 0 || done[j] > total) {
      stop("`done` must be cumulative counts between 0 and `total`");
    }
  }
  const double c1 = done.back(), t1 = ye.back() + 1, span = t1 - start;
  if (span <= 0) stop("the last count must come after `start`");
  const double a = 1 - cl;
  const double lo = c1 > 0 ? clab_qgamma(a / 2, c1, 1) / span : 0;
  const double hi = clab_qgamma(1 - a / 2, c1 + 1, 1) / span;
  const double rem = total - c1;
  auto by = [&](double r) { return rem == 0 ? t1 : r > 0 ? t1 + rem / r : HUGE_VAL; };
  V out = frame({"completed", "remaining", "years", "rate", "ratelo", "rateup", "cmplby", "cmplbyearly", "cmplbylate"},
                {dv1(c1), dv1(rem), dv1(span), dv1(c1 / span), dv1(lo), dv1(hi), dv1(by(c1 / span)), dv1(by(hi)),
                 dv1(by(lo))});
  std::vector<double> yr{start}, cm{0};
  for (std::size_t j = 0; j < ye.size(); ++j) {
    yr.push_back(ye[j] + 1);
    cm.push_back(done[j]);
  }
  out.at.push_back({"history", frame({"year", "completed"}, {dv(yr), dv(cm)})});
  out.at.push_back({"total", dv1(total)});
  out.cls = {"trcprog", "data.frame"};
  return out;
}

}  // namespace igh
