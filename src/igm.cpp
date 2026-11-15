// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// First Nations, Metis and Inuit mental health and wellbeing from
// Statistics Canada's published tables: the estimates and their limits
// taken out of a table's statistics column (scests), rated scales compared
// as whole distributions (ordcmp), survey percentages standardised by age
// (stdprev), wellness beside distress (welprof, dualcont) and the balance
// of wellness across domains (holbal). igh.h.
#include <algorithm>
#include <climits>
#include <cmath>
#include <map>
#include <regex>
#include <cstring>
#include <tuple>
#include <set>

#include "crimlab.h"
#include "igk.h"

namespace igh {

V f_sctidy(const Cx &c);

namespace {

bool in(const std::string &s, const std::vector<std::string> &v) { return std::find(v.begin(), v.end(), s) != v.end(); }

std::vector<std::string> tl(const V &v) {
  const Lab l = tlab(v);
  std::vector<std::string> o(l.n());
  for (std::size_t j = 0; j < l.n(); ++j) o[j] = l.na[j] ? "NA" : l.s[j];
  return o;
}

bool istext(const V &v) { return v.t == V::STR || v.t == V::FAC; }

// the row of indicator ind for each group, or SIZE_MAX
std::vector<std::size_t> rowsof(const Lab &il, const std::vector<std::string> &gs,
                                const std::vector<std::string> &groups, const std::string &ind) {
  std::vector<std::size_t> o;
  for (const std::string &gi : groups) {
    std::size_t at = SIZE_MAX;
    for (std::size_t j = 0; j < gs.size() && at == SIZE_MAX; ++j) {
      if (!il.na[j] && il.s[j] == ind && gs[j] == gi) at = j;
    }
    o.push_back(at);
  }
  return o;
}

// a column argument checked in the words of the R package
std::string ccol(const V &data, const V &v, const char *what) {
  const std::string s = str1(v, what);
  if (bad()) return {};
  if (!col(data, s)) return stop("column `%s` not found in `data`", s.c_str());
  return s;
}

std::string gsubws(const std::string &s) {
  std::string o;
  bool sp = false;
  for (char ch : s) {
    if (std::isspace(static_cast<unsigned char>(ch))) {
      if (!sp) o += ' ';
      sp = true;
    } else {
      o += ch;
      sp = false;
    }
  }
  return o;
}

// "^(Statistics|Characteristics)\b"
bool statname(const std::string &n) {
  for (const char *w : {"Statistics", "Characteristics"}) {
    const std::size_t m = std::strlen(w);
    if (n.compare(0, m, w) == 0) {
      if (n.size() == m) return true;
      const unsigned char ch = static_cast<unsigned char>(n[m]);
      if (!(std::isalnum(ch) || ch == '_')) return true;
    }
  }
  return false;
}

bool hasstatwords(const V &v) {
  const Lab l = lab(v);
  for (std::size_t j = 0; j < l.n(); ++j) {
    if (l.na[j]) continue;
    const std::string s = lower(l.s[j]);
    if (s.find("confidence interval") != std::string::npos || s.find("coefficient of variation") != std::string::npos ||
        s.find("standard error") != std::string::npos) {
      return true;
    }
  }
  return false;
}

}  // namespace

V f_scests(const Cx &c) {
  V x = c["x"];
  nofn(x, "x");
  if (bad()) return {};
  if (!inh(x, "sctidy")) {
    std::vector<std::string> an{"x"};
    std::vector<V> av{x};
    Cx c2{&an, &av, c.warn};
    x = f_sctidy(c2);
    if (bad()) return {};
  }
  static const std::set<std::string> meta = {"UOM",    "UOM_ID", "SCALAR_FACTOR", "SCALAR_ID", "VECTOR",   "COORDINATE",
                                             "Coordinate", "VALUE", "STATUS",     "SYMBOL",    "TERMINATED", "DECIMALS",
                                             "quality"};
  std::vector<std::string> cols;
  for (const std::string &n : x.nm) {
    if (!meta.count(n)) cols.push_back(n);
  }
  const std::size_t n = nrow(x);
  std::vector<std::string> uom(n, "");
  if (const V *u = col(x, "UOM")) {
    const Lab l = tlab(*u);
    for (std::size_t r = 0; r < n; ++r) uom[r] = l.na[r] ? "NA" : l.s[r];
  }
  std::vector<char> st(cols.size());
  int ns = 0;
  for (std::size_t j = 0; j < cols.size(); ++j) {
    const V &v = *col(x, cols[j]);
    st[j] = istext(v) && (statname(cols[j]) || hasstatwords(v));
    ns += st[j];
  }
  if (ns > 1) {
    std::vector<std::size_t> named;
    std::size_t first = SIZE_MAX;
    for (std::size_t j = 0; j < cols.size(); ++j) {
      if (!st[j]) continue;
      if (first == SIZE_MAX) first = j;
      if (cols[j].compare(0, 10, "Statistics") == 0) named.push_back(j);
    }
    const std::size_t keep = named.size() == 1 ? named[0] : first;
    for (std::size_t j = 0; j < cols.size(); ++j) st[j] = j == keep;
    ns = 1;
  }
  // .sce_trim(): the text columns trimmed
  for (const std::string &cn : cols) {
    V &v = x.a[static_cast<std::size_t>(coli(x, cn))];
    if (istext(v)) {
      const Lab l = tlab(v);
      v = svna(l.s, l.na);
    }
  }
  const std::vector<double> value = col(x, "VALUE")->d;
  const V qual = *col(x, "quality");
  if (ns == 0) {
    std::vector<std::size_t> all(n);
    for (std::size_t r = 0; r < n; ++r) all[r] = r;
    V out;
    out.t = V::LIST;
    out.df = true;
    for (const std::string &cn : cols) {
      out.nm.push_back(cn);
      out.a.push_back(*col(x, cn));
    }
    const std::vector<double> nas(n, na());
    out.nm.insert(out.nm.end(), {"measure", "estimate", "upper", "lower", "cisou", "quality"});
    out.a.push_back(sv(uom));
    out.a.push_back(dv(value));
    out.a.push_back(dv(nas));
    out.a.push_back(dv(nas));
    out.a.push_back(svna(std::vector<std::string>(n), std::vector<char>(n, 1)));
    out.a.push_back(qual);
    out.nr = static_cast<int>(n);
    out.cls = {"sctidy", "data.frame"};
    out.at.push_back({"conflvl", dv1(0.95)});
    out.at.push_back({"unmatched", sv({})});
    return out;
  }
  std::string sc;
  for (std::size_t j = 0; j < cols.size(); ++j) {
    if (st[j]) sc = cols[j];
  }
  const Lab sl = lab(*col(x, sc));
  std::vector<std::string> lb(n), low(n);
  for (std::size_t r = 0; r < n; ++r) {
    lb[r] = sl.na[r] ? "NA" : gsubws(trim(sl.s[r]));
    low[r] = lower(lb[r]);
  }
  // .sce_roles()
  std::vector<std::string> role(n, "estimate"), target(n);
  std::vector<char> tna(n, 1);
  const std::vector<std::pair<std::string, std::vector<std::string>>> pats = {
      {"lower",
       {"^low [0-9.]+% confidence interval,? ?(.*)$", "^confidence interval lower bound,? ?(?:of )?(.*)$",
        "^[0-9.]+% confidence interval \\(lower limit\\),? ?(.*)$", "^[0-9.]+% confidence interval lower bound,? ?(.*)$"}},
      {"upper",
       {"^high [0-9.]+% confidence interval,? ?(.*)$", "^confidence interval upper bound,? ?(?:of )?(.*)$",
        "^[0-9.]+% confidence interval \\(upper limit\\),? ?(.*)$", "^[0-9.]+% confidence interval upper bound,? ?(.*)$"}},
      {"cv", {"^coefficient of variation,? ?(?:for |of )?(.*)$"}},
      {"se", {"^standard error,? ?(?:for |of )?(.*)$"}}};
  for (const auto &pr : pats) {
    for (const std::string &pt : pr.second) {
      const std::regex re(pt);
      for (std::size_t r = 0; r < n; ++r) {
        std::smatch m;
        if (role[r] == "estimate" && std::regex_search(low[r], m, re)) {
          role[r] = pr.first;
          target[r] = trim(m[1].str());
          tna[r] = 0;
        }
      }
    }
  }
  // .sce_conf()
  double conf = 0.95;
  {
    const std::regex re("[0-9.]+(?=% confidence)");
    std::vector<std::string> lv;
    for (std::size_t r = 0; r < n; ++r) {
      std::smatch m;
      if (std::regex_search(low[r], m, re)) lv.push_back(m[0].str());
    }
    lv = uniq(lv);
    if (lv.size() == 1) {
      const std::vector<double> cv = num(sv1(lv[0]), "conf");
      if (bad()) return {};
      conf = cv[0] / 100;
    }
  }
  std::vector<std::string> dims;
  for (const std::string &cn : cols) {
    if (cn != sc) dims.push_back(cn);
  }
  std::vector<std::string> key(n);
  {
    std::vector<Lab> dl;
    for (const std::string &d : dims) dl.push_back(lab(*col(x, d)));
    for (std::size_t r = 0; r < n; ++r) {
      std::string k;
      for (std::size_t j = 0; j < dl.size(); ++j) k += (j ? "\r" : "") + (dl[j].na[r] ? std::string("NA") : dl[j].s[r]);
      key[r] = k;
    }
  }
  std::vector<std::size_t> e;
  std::map<std::string, std::size_t> ekey;
  for (std::size_t r = 0; r < n; ++r) {
    if (role[r] != "estimate") continue;
    const std::string k = key[r] + "\r" + lb[r];
    if (ekey.count(k)) return stop("`x` has two estimates for one cell and measure of the table");
    ekey[k] = e.size();
    e.push_back(r);
  }
  std::vector<std::string> ests;
  for (std::size_t r : e) ests.push_back(lb[r]);
  ests = uniq(ests);
  // .sce_limit_fit(): the share of a limit's cells on the right side of
  // the estimates of measure cn
  auto fit = [&](const std::string &cn, const std::vector<std::size_t> &i, const std::vector<std::size_t> &j,
                 const std::string &r) {
    std::map<std::string, std::size_t> kj;
    for (std::size_t t : j) {
      if (!kj.count(key[t])) kj[key[t]] = t;
    }
    double hits = 0, cnt = 0;
    for (std::size_t t : i) {
      const auto f = ekey.find(key[t] + "\r" + cn);
      const double ev = f == ekey.end() ? na() : value[e[f->second]];
      const double bv = value[t];
      bool ok = std::isfinite(ev) && std::isfinite(bv);
      bool ft = r == "lower" ? bv <= ev + 1e-9 : bv >= ev - 1e-9;
      if (!j.empty()) {
        const auto g = kj.find(key[t]);
        const double ov = g == kj.end() ? na() : value[g->second];
        ok = ok && std::isfinite(ov);
        ft = ft && (r == "lower" ? ov >= ev - 1e-9 : ov <= ev + 1e-9);
      }
      if (ok) {
        cnt += 1;
        hits += ft;
      }
    }
    return cnt > 0 ? hits / cnt : 0.0;
  };
  // .sce_bounds(): the measure each precision label belongs to
  std::vector<std::string> bls, unmatched;
  for (std::size_t r = 0; r < n; ++r) {
    if (role[r] != "estimate") bls.push_back(lb[r]);
  }
  bls = uniq(bls);
  std::vector<std::pair<std::string, std::pair<std::string, std::string>>> bound;
  for (const std::string &bl : bls) {
    std::vector<std::size_t> i;
    for (std::size_t r = 0; r < n; ++r) {
      if (lb[r] == bl) i.push_back(r);
    }
    const std::string r = role[i[0]], t = target[i[0]];
    std::vector<std::string> cand = ests;
    if (!t.empty()) {
      std::vector<std::string> ex, fx;
      for (const std::string &s : ests) {
        if (lower(s) == t) ex.push_back(s);
        if (lower(s).find(t) != std::string::npos) fx.push_back(s);
      }
      cand = ex.empty() ? fx : ex;
    }
    bool anyu = false;
    for (const std::string &u : uom) anyu = anyu || !u.empty();
    if (cand.size() > 1 && anyu) {
      std::set<std::string> u;
      for (std::size_t k : i) u.insert(uom[k]);
      std::vector<std::string> byu;
      for (const std::string &cn : cand) {
        bool hit = false;
        for (std::size_t k = 0; k < n; ++k) hit = hit || (lb[k] == cn && u.count(uom[k]));
        if (hit) byu.push_back(cn);
      }
      if (!byu.empty()) cand = byu;
    }
    if (cand.size() > 1 && (r == "lower" || r == "upper")) {
      const std::string orole = r == "lower" ? "upper" : "lower";
      std::vector<std::string> other;
      for (std::size_t k = 0; k < n; ++k) {
        if (role[k] == orole && !tna[k] && target[k] == t) other.push_back(lb[k]);
      }
      other = uniq(other);
      std::vector<std::size_t> j;
      if (other.size() == 1) {
        for (std::size_t k = 0; k < n; ++k) {
          if (lb[k] == other[0]) j.push_back(k);
        }
      }
      std::vector<double> fits;
      double mx = -HUGE_VAL;
      for (const std::string &cn : cand) {
        fits.push_back(fit(cn, i, j, r));
        mx = std::max(mx, fits.back());
      }
      std::vector<std::string> keep;
      for (std::size_t k = 0; k < cand.size(); ++k) {
        if (fits[k] >= 0.95 && fits[k] == mx) keep.push_back(cand[k]);
      }
      cand = keep;
    }
    if (cand.size() == 1) bound.push_back({bl, {r, cand[0]}});
    else unmatched.push_back(bl);
  }
  // .sce_out() and .sce_intervals()
  const std::size_t ne = e.size();
  const double NA = na();
  std::vector<double> est(ne), lo(ne, NA), hi(ne, NA), cv(ne, NA), se(ne, NA);
  for (std::size_t k = 0; k < ne; ++k) est[k] = value[e[k]];
  for (const auto &b : bound) {
    for (std::size_t r = 0; r < n; ++r) {
      if (lb[r] != b.first) continue;
      const auto f = ekey.find(key[r] + "\r" + b.second.second);
      if (f == ekey.end()) continue;
      const std::string &ro = b.second.first;
      if (ro == "lower") lo[f->second] = value[r];
      else if (ro == "upper") hi[f->second] = value[r];
      else if (ro == "cv") cv[f->second] = value[r];
      else se[f->second] = value[r];
    }
  }
  const double z = qn(conf);
  std::vector<std::string> cs(ne);
  std::vector<char> csn(ne, 1);
  for (std::size_t k = 0; k < ne; ++k) {
    const bool pub = std::isfinite(lo[k]) && std::isfinite(hi[k]);
    if (pub) {
      cs[k] = "published";
      csn[k] = 0;
      if (std::isfinite(est[k]) && (lo[k] > est[k] || hi[k] < est[k])) cs[k] = "inconsistent";
    } else if (std::isfinite(cv[k]) && std::isfinite(est[k])) {
      const double sg = std::sqrt(std::log(1 + (cv[k] / 100) * (cv[k] / 100)));
      lo[k] = est[k] * std::exp(-z * sg);
      hi[k] = est[k] * std::exp(z * sg);
      cs[k] = "cv";
      csn[k] = 0;
    } else if (std::isfinite(se[k]) && std::isfinite(est[k])) {
      lo[k] = est[k] - z * se[k];
      hi[k] = est[k] + z * se[k];
      cs[k] = "se";
      csn[k] = 0;
    }
    if (!std::isfinite(est[k])) {
      lo[k] = hi[k] = NA;
      csn[k] = 1;
      cs[k].clear();
    }
  }
  V out;
  out.t = V::LIST;
  out.df = true;
  for (const std::string &d : dims) {
    out.nm.push_back(d);
    out.a.push_back(rows(*col(x, d), e));
  }
  out.nm.insert(out.nm.end(), {"measure", "estimate", "lower", "upper", "cisou", "quality"});
  std::vector<std::string> ms;
  for (std::size_t r : e) ms.push_back(lb[r]);
  out.a.push_back(sv(ms));
  out.a.push_back(dv(est));
  out.a.push_back(dv(lo));
  out.a.push_back(dv(hi));
  out.a.push_back(svna(cs, csn));
  out.a.push_back(rows(qual, e));
  out.nr = static_cast<int>(ne);
  out.cls = {"sctidy", "data.frame"};
  out.at.push_back({"conflvl", dv1(conf)});
  out.at.push_back({"unmatched", sv(unmatched)});
  return out;
}

V f_ordcmp(const Cx &c) {
  const V &data = c["data"];
  if (!isdf(data)) return stop("`data` must be a data frame");
  std::string cn[5];
  static const char *const an[] = {"grp", "cat", "est", "lower", "upper"};
  for (int k = 0; k < 5; ++k) cn[k] = ccol(data, c[an[k]], an[k]);
  std::string qn1;
  if (c.has("qual")) qn1 = ccol(data, c["qual"], "qual");
  const V &lv0 = c["lvls"];
  nofn(lv0, "lvls");
  if (bad()) return {};
  if (!isstr(lv0) || lv0.s.size() < 2 || uniq(lv0.s).size() != lv0.s.size()) {
    return stop("`lvls` must be two or more distinct categories, from the lowest to the highest");
  }
  const std::string ref = str1(c["ref"], "ref");
  const double cl = conf(c["conflvl"], "conflvl", 0.95);
  if (bad()) return {};
  const double z = qn(cl);
  const std::vector<std::string> g = tl(*col(data, cn[0]));
  const V &pooled = c["pooled"];
  pooledcheck(pooled, g);
  if (bad()) return {};
  const std::vector<std::string> cat = tl(*col(data, cn[1]));
  std::vector<std::string> lvls;
  for (const std::string &s : lv0.s) lvls.push_back(trim(s));
  const std::vector<double> est = num(*col(data, cn[2]), "est");
  const std::vector<double> lo = num(*col(data, cn[3]), "lower");
  const std::vector<double> hi = num(*col(data, cn[4]), "upper");
  if (bad()) return {};
  const std::vector<std::string> q = quality(c, data, qn1, g.size());
  if (!in(ref, g)) return stop("\"%s\" is not a group of `data`", ref.c_str());
  bool anyc = false;
  for (const std::string &s : cat) anyc = anyc || in(s, lvls);
  if (!anyc) return stop("no row of `data` has a category in `lvls`");
  const std::vector<std::string> groups = uniq(g);
  const std::size_t k = lvls.size(), ng = groups.size();
  const double NA = na();
  std::vector<double> sh(ng * k, NA), neff(ng, NA);
  std::vector<int> caut(ng, 0);
  for (std::size_t gi = 0; gi < ng; ++gi) {
    std::vector<std::size_t> i;
    std::set<std::string> seen;
    for (std::size_t j = 0; j < g.size(); ++j) {
      if (g[j] == groups[gi] && in(cat[j], lvls)) {
        if (!seen.insert(cat[j]).second) return stop("group \"%s\" has a category twice", groups[gi].c_str());
        i.push_back(j);
      }
    }
    std::vector<std::size_t> ii;
    bool miss = false;
    for (const std::string &lv : lvls) {
      std::size_t at = SIZE_MAX;
      for (std::size_t j : i) {
        if (cat[j] == lv) at = j;
      }
      if (at == SIZE_MAX) miss = true;
      ii.push_back(at);
    }
    if (miss) continue;
    bool bad = false;
    for (std::size_t j : ii) bad = bad || gone(q[j]) || !std::isfinite(est[j]);
    if (bad) continue;
    bool out = false;
    for (std::size_t j : ii) out = out || (lo[j] > est[j]) || (hi[j] < est[j]);
    if (out) {
      c.w(sf("an interval of \"%s\" does not contain its estimate; not compared", groups[gi].c_str()));
      continue;
    }
    long double tot = 0;
    for (std::size_t j : ii) tot += est[j];
    const double t = static_cast<double>(tot);
    std::vector<double> ne;
    bool ce = false;
    for (std::size_t m = 0; m < k; ++m) {
      const std::size_t j = ii[m];
      const double p = est[j] / t;
      const double s = (hi[j] - lo[j]) / (2 * z) / t;
      sh[gi * k + m] = p;
      if (std::isfinite(s) && s > 0 && p > 0 && p < 1) ne.push_back(p * (1 - p) / (s * s));
      ce = ce || q[j] == "E";
    }
    neff[gi] = ne.empty() ? NA : median(ne);
    caut[gi] = ce;
  }
  const std::size_t ri = static_cast<std::size_t>(std::find(groups.begin(), groups.end(), ref) - groups.begin());
  for (std::size_t m = 0; m < k; ++m) {
    if (std::isnan(sh[ri * k + m])) {
      return stop("the reference group's distribution is suppressed or incomplete: there is nothing to compare with");
    }
  }
  std::vector<double> q0(k), g0(k);
  {
    long double cs = 0;
    for (std::size_t m = 0; m < k; ++m) {
      q0[m] = sh[ri * k + m];
      const double below = static_cast<double>(cs);
      g0[m] = below + q0[m] / 2;
      cs += q0[m];
    }
  }
  const double n0 = neff[ri];
  std::vector<std::string> og;
  std::vector<double> op, ol, ou, od, on;
  std::vector<int> oc;
  for (std::size_t gi = 0; gi < ng; ++gi) {
    if (gi == ri) continue;
    og.push_back(groups[gi]);
    on.push_back(neff[gi]);
    oc.push_back(caut[gi] || caut[ri]);
    std::vector<double> p1(k);
    bool nap = false;
    for (std::size_t m = 0; m < k; ++m) {
      p1[m] = sh[gi * k + m];
      nap = nap || std::isnan(p1[m]);
    }
    if (nap) {
      op.push_back(NA);
      ol.push_back(NA);
      ou.push_back(NA);
      od.push_back(NA);
      continue;
    }
    long double a = 0, s1 = 0, s0 = 0;
    for (std::size_t m = 0; m < k; ++m) a += static_cast<long double>(p1[m] * g0[m]);
    const double ad = static_cast<double>(a);
    std::vector<double> g1(k);
    {
      long double cs = 0;
      std::vector<double> ab(k);
      // the share strictly above each category
      for (std::size_t m = k; m-- > 0;) {
        ab[m] = static_cast<double>(cs);
        cs += p1[m];
      }
      for (std::size_t m = 0; m < k; ++m) g1[m] = ab[m] + p1[m] / 2;
    }
    for (std::size_t m = 0; m < k; ++m) {
      s1 += static_cast<long double>(p1[m] * (g0[m] * g0[m]));
      s0 += static_cast<long double>(q0[m] * (g1[m] * g1[m]));
    }
    const double v = (static_cast<double>(s1) - ad * ad) / neff[gi] + (static_cast<double>(s0) - ad * ad) / n0;
    op.push_back(ad);
    od.push_back(2 * ad - 1);
    if (!std::isfinite(v) || ad <= 0 || ad >= 1) {
      ol.push_back(NA);
      ou.push_back(NA);
      continue;
    }
    const double sl = std::sqrt(v) / (ad * (1 - ad));
    ol.push_back(plogis(qlogis(ad) - z * sl));
    ou.push_back(plogis(qlogis(ad) + z * sl));
  }
  const std::vector<int> pg = pooledrows(og, pooled);
  if (bad()) return {};
  V out = frame({"group", "reference", "probbet", "lower", "upper", "delta", "neff", "neffref", "usewithcau", "poogrp"},
                {sv(og), sv(std::vector<std::string>(og.size(), ref)), dv(op), dv(ol), dv(ou), dv(od), dv(on),
                 dv(std::vector<double>(og.size(), n0)), lv(oc), lv(pg)});
  V m = dv(std::vector<double>(ng * k));
  for (std::size_t gi = 0; gi < ng; ++gi) {
    for (std::size_t j = 0; j < k; ++j) m.d[j * ng + gi] = sh[gi * k + j];
  }
  m.nr = static_cast<int>(ng);
  m.nc = static_cast<int>(k);
  m.rn = groups;
  m.cn = lvls;
  out.at.push_back({"shares", m});
  out.at.push_back({"conflvl", dv1(cl)});
  out.cls = {"ordcmp", "data.frame"};
  return out;
}

V f_stdprev(const Cx &c) {
  std::vector<double> est = num(c["est"], "est");
  const std::vector<double> lo = num(c["lower"], "lower"), hi = num(c["upper"], "upper");
  const std::vector<double> std = num(c["std"], "std");
  const double cl = conf(c["conflvl"], "conflvl", 0.95);
  double scale = 100;
  if (!dflt(c["scale"])) {
    const V &s = c["scale"];
    const bool inf = (s.t == V::DBL && s.d.size() == 1 && s.d[0] == HUGE_VAL) ||
                     (s.t == V::STR && s.s.size() == 1 && s.s[0] == "Inf");
    scale = inf ? HUGE_VAL : scalar(s, "scale", 0, HUGE_VAL, false, 100);
  }
  if (bad()) return {};
  const std::size_t k = std.size();
  bool bad = k == 0;
  long double ss = 0;
  for (double x : std) {
    bad = bad || std::isnan(x) || x < 0;
    ss += x;
  }
  if (bad || ss <= 0) return stop("`std` must be non-negative counts, one per age group, not all zero");
  if (lo.size() != est.size() || hi.size() != est.size()) return stop("`est`, `lower` and `upper` must be the same length");
  const std::size_t n = est.size();
  const V &gv = c["grp"], &qv = c["qual"];
  nofn(gv, "grp");
  nofn(qv, "qual");
  if (igh::bad()) return {};
  std::vector<std::string> grp(n, "all");
  if (gv.t != V::NUL) {
    const Lab l = lab(gv);
    if (l.n() != n) return stop("`grp` must have one label per row");
    for (std::size_t j = 0; j < n; ++j) grp[j] = l.na[j] ? "NA" : l.s[j];
  }
  Lab ql;
  if (qv.t == V::NUL) {
    ql.s.assign(n, "ok");
    ql.na.assign(n, 0);
  } else {
    ql = lab(qv);
  }
  if (ql.n() != n) return stop("`qual` must have one flag per row");
  const std::vector<std::string> q = qflags(c, ql);
  int nb = 0;
  for (std::size_t j = 0; j < n; ++j) {
    if (std::isfinite(lo[j]) && std::isfinite(hi[j]) && std::isfinite(est[j]) && (lo[j] > est[j] || hi[j] < est[j])) {
      ++nb;
      est[j] = na();
    }
  }
  if (nb) c.w(sf("%d interval(s) do not contain their estimate; their groups are not standardised", nb));
  const double z = qn(cl);
  std::vector<double> w(k);
  for (std::size_t j = 0; j < k; ++j) w[j] = std[j] / static_cast<double>(ss);
  std::vector<std::string> og;
  std::vector<double> oe, ol, ou, os;
  std::vector<int> oc;
  const double NA = na();
  for (const std::string &gi : uniq(grp)) {
    std::vector<std::size_t> i;
    for (std::size_t j = 0; j < n; ++j) {
      if (grp[j] == gi) i.push_back(j);
    }
    if (i.size() != k) {
      return stop("group \"%s\" has %d age groups; `std` has %d", gi.c_str(), static_cast<int>(i.size()), static_cast<int>(k));
    }
    og.push_back(gi);
    bool gn = false, ce = false;
    for (std::size_t j : i) {
      gn = gn || gone(q[j]) || !std::isfinite(est[j]) || !std::isfinite(lo[j]) || !std::isfinite(hi[j]);
      ce = ce || q[j] == "E";
    }
    if (gn) {
      oe.push_back(NA);
      ol.push_back(NA);
      ou.push_back(NA);
      os.push_back(NA);
      oc.push_back(0);
      continue;
    }
    long double y = 0, v = 0;
    for (std::size_t m = 0; m < k; ++m) {
      const std::size_t j = i[m];
      y += static_cast<long double>(w[m] * est[j]);
      const double s = (hi[j] - lo[j]) / (2 * z);
      v += static_cast<long double>(w[m] * w[m] * (s * s));
    }
    const double yd = static_cast<double>(y), se = std::sqrt(static_cast<double>(v));
    oe.push_back(yd);
    ol.push_back(std::max(0.0, yd - z * se));
    ou.push_back(std::min(scale, yd + z * se));
    os.push_back(se);
    oc.push_back(ce);
  }
  V out = frame({"group", "estimate", "lower", "upper", "se", "usewithcau"},
                {sv(og), dv(oe), dv(ol), dv(ou), dv(os), lv(oc)});
  out.at.push_back({"conflvl", dv1(cl)});
  out.cls = {"stdprev", "data.frame"};
  return out;
}

V f_welprof(const Cx &c) {
  const V &data = c["data"];
  if (!isdf(data)) return stop("`data` must be a data frame");
  std::string cn[5];
  static const char *const an[] = {"grp", "ind", "est", "lower", "upper"};
  for (int k = 0; k < 5; ++k) cn[k] = ccol(data, c[an[k]], an[k]);
  std::string qn1;
  if (c.has("qual")) qn1 = ccol(data, c["qual"], "qual");
  const V &wv = c["well"], &dv0 = c["distress"], &dom = c["domain"];
  nofn(wv, "well");
  nofn(dv0, "distress");
  nofn(dom, "domain");
  if (bad()) return {};
  const std::vector<std::string> ind = tl(*col(data, cn[1]));
  if (!isstr(wv)) return stop("`well` must name indicators of `data`");
  if (!isstr(dv0)) return stop("`distress` must name indicators of `data`");
  std::vector<std::string> well, dis;
  for (const std::string &s : wv.s) well.push_back(trim(s));
  for (const std::string &s : dv0.s) dis.push_back(trim(s));
  for (const std::string &s : well) {
    if (in(s, dis)) return stop("an indicator cannot be both one of wellness and one of distress");
  }
  std::vector<std::string> both = well, miss;
  both.insert(both.end(), dis.begin(), dis.end());
  for (const std::string &s : uniq(both)) {
    if (!in(s, ind)) miss.push_back(s);
  }
  if (!miss.empty()) return stop("not an indicator of `data`: %s", tostr(miss).c_str());
  V dm = dom;
  if (dm.t == V::LIST && dm.js) {
    V d;
    d.t = V::STR;
    bool ok = true;
    for (const V &e : dm.a) {
      ok = ok && e.t == V::STR && e.s.size() == 1;
      if (ok) {
        d.s.push_back(e.s[0]);
        d.na.push_back(e.na.empty() ? 0 : e.na[0]);
      }
    }
    if (ok) {
      d.nm = dm.nm;
      dm = d;
    }
  }
  if (dm.t != V::NUL && (dm.t != V::STR || dm.nm.empty())) {
    return stop("`domain` must be a named character vector: indicator = domain");
  }
  std::vector<double> est = num(*col(data, cn[2]), "est");
  if (bad()) return {};
  const std::vector<std::string> q = quality(c, data, qn1, ind.size());
  const std::size_t n = ind.size();
  std::vector<std::string> domv(n);
  std::vector<char> domn(n, 1);
  if (dm.t != V::NUL) {
    for (std::size_t j = 0; j < n; ++j) {
      for (std::size_t k = 0; k < dm.nm.size(); ++k) {
        if (dm.nm[k] == ind[j]) {
          domv[j] = dm.s[k];
          domn[j] = dm.na.empty() ? 0 : dm.na[k];
          break;
        }
      }
    }
  }
  const std::vector<double> lr = num(*col(data, cn[3]), "lower"), hr = num(*col(data, cn[4]), "upper");
  if (bad()) return {};
  std::vector<double> lo(n), hi(n);
  std::vector<char> keep(n);
  const std::vector<std::string> g = tl(*col(data, cn[0]));
  std::vector<std::string> labs(n);
  for (std::size_t j = 0; j < n; ++j) {
    const bool gn = gone(q[j]);
    if (gn) est[j] = na();
    lo[j] = gn ? na() : lr[j];
    hi[j] = gn ? na() : hr[j];
    keep[j] = in(ind[j], both);
    labs[j] = g[j] + ": " + ind[j];
  }
  contains(c, "welprof", est, lo, hi, labs, &keep);
  std::vector<std::size_t> kr;
  for (std::size_t j = 0; j < n; ++j) {
    if (keep[j]) kr.push_back(j);
  }
  std::vector<std::string> ug;
  for (std::size_t j : kr) ug.push_back(g[j]);
  ug = uniq(ug);
  auto pos = [](const std::vector<std::string> &v, const std::string &s) {
    return static_cast<int>(std::find(v.begin(), v.end(), s) - v.begin());
  };
  std::stable_sort(kr.begin(), kr.end(), [&](std::size_t a, std::size_t b) {
    const auto ka = std::make_tuple(!in(ind[a], well), pos(both, ind[a]), pos(ug, g[a]));
    const auto kb = std::make_tuple(!in(ind[b], well), pos(both, ind[b]), pos(ug, g[b]));
    return ka < kb;
  });
  std::vector<std::string> og, oi, os, od;
  std::vector<char> odn;
  std::vector<double> oe, ol, ou;
  std::vector<int> oc;
  for (std::size_t j : kr) {
    og.push_back(g[j]);
    oi.push_back(ind[j]);
    os.push_back(in(ind[j], well) ? "wellness" : "distress");
    od.push_back(domv[j]);
    odn.push_back(domn[j]);
    oe.push_back(est[j]);
    ol.push_back(lo[j]);
    ou.push_back(hi[j]);
    oc.push_back(q[j] == "E");
  }
  V out = frame({"group", "indicator", "side", "domain", "estimate", "lower", "upper", "usewithcau"},
                {sv(og), sv(oi), sv(os), svna(od, odn), dv(oe), dv(ol), dv(ou), lv(oc)});
  out.cls = {"welprof", "data.frame"};
  return out;
}

V f_holbal(const Cx &c) {
  const V &data = c["data"];
  if (!isdf(data)) return stop("`data` must be a data frame");
  const std::string gc = ccol(data, c["grp"], "grp"), dc = ccol(data, c["domain"], "domain"),
                    ec = ccol(data, c["est"], "est"), lc = ccol(data, c["lower"], "lower"),
                    uc = ccol(data, c["upper"], "upper");
  std::string qn1;
  if (c.has("qual")) qn1 = ccol(data, c["qual"], "qual");
  const int reps = static_cast<int>(scalar(c["reps"], "reps", 100, 1e6, true, 10000));
  const double cl = conf(c["conflvl"], "conflvl", 0.95);
  if (bad()) return {};
  const double z = qn(cl);
  const std::vector<std::string> g = tl(*col(data, gc));
  const Lab dml = tlab(*col(data, dc));
  std::vector<double> raw = num(*col(data, ec), "est");
  const std::vector<double> lr = num(*col(data, lc), "lower"), hr = num(*col(data, uc), "upper");
  if (bad()) return {};
  {
    bool anyf = false, all1 = true;
    for (double x : raw) {
      if (std::isfinite(x)) {
        anyf = true;
        all1 = all1 && x <= 1;
      }
    }
    if (anyf && all1) return stop("the estimates are all at or below 1: are they proportions? give percentages (0 to 100)");
  }
  const std::size_t n = g.size();
  const std::vector<std::string> q = quality(c, data, qn1, n);
  std::vector<char> indom(n);
  std::vector<std::string> dm(n), labs(n);
  for (std::size_t j = 0; j < n; ++j) {
    if (gone(q[j])) raw[j] = na();
    indom[j] = !dml.na[j] && !dml.s[j].empty();
    dm[j] = dml.na[j] ? "\x01NA" : dml.s[j];
    labs[j] = g[j] + ": " + (dml.na[j] ? "NA" : dml.s[j]);
  }
  contains(c, "holbal", raw, lr, hr, labs, &indom);
  std::vector<double> est(n), se(n);
  for (std::size_t j = 0; j < n; ++j) {
    est[j] = raw[j] / 100;
    se[j] = (hr[j] - lr[j]) / (2 * z) / 100;
  }
  std::vector<std::string> doms;
  for (std::size_t j = 0; j < n; ++j) {
    if (indom[j]) doms.push_back(dm[j]);
  }
  doms = uniq(doms);
  const std::size_t k = doms.size();
  if (k < 2) return stop("balance needs two domains at least");
  V wv = c["w"];
  std::vector<double> w(k, 1);
  if (wv.t != V::NUL) {
    nofn(wv, "w");
    if (bad()) return {};
    if (wv.t == V::LIST && wv.js) {
      V d;
      d.t = V::DBL;
      for (const V &e : wv.a) {
        if (e.t == V::DBL && e.d.size() == 1) d.d.push_back(e.d[0]);
      }
      if (d.d.size() == wv.a.size()) {
        d.nm = wv.nm;
        wv = d;
      }
    }
    bool bad = !(wv.t == V::DBL || wv.t == V::INT) || wv.nm.empty();
    std::vector<double> x;
    if (!bad) {
      x = num(wv, "w");
      if (igh::bad()) return {};
      std::set<std::string> a(wv.nm.begin(), wv.nm.end()), b(doms.begin(), doms.end());
      bad = a != b;
      double s = 0;
      for (double v : x) {
        bad = bad || !std::isfinite(v) || v < 0;
        s += v;
      }
      bad = bad || !(s > 0);
    }
    if (bad) return stop("`w` must be non-negative numbers named by the domains (%s), not all zero", tostr(doms).c_str());
    for (std::size_t m = 0; m < k; ++m) {
      for (std::size_t j = 0; j < wv.nm.size(); ++j) {
        if (wv.nm[j] == doms[m]) {
          w[m] = x[j];
          break;
        }
      }
    }
  }
  {
    long double s = 0;
    for (double v : w) s += v;
    for (double &v : w) v /= static_cast<double>(s);
  }
  const std::vector<std::string> groups = uniq(g);
  const std::size_t ng = groups.size();
  const double NA = na();
  std::vector<double> lev(ng * k, NA), sem(ng * k, NA);
  std::vector<int> caut(ng, 0);
  for (std::size_t gi = 0; gi < ng; ++gi) {
    for (std::size_t m = 0; m < k; ++m) {
      std::vector<double> ev, sv2;
      bool nae = false;
      for (std::size_t j = 0; j < n; ++j) {
        if (g[j] == groups[gi] && dm[j] == doms[m]) {
          ev.push_back(est[j]);
          sv2.push_back(se[j] * se[j]);
          nae = nae || std::isnan(est[j]);
        }
      }
      if (ev.empty() || nae) continue;
      lev[gi * k + m] = rmean(ev);
      sem[gi * k + m] = std::sqrt(rsum(sv2)) / static_cast<double>(ev.size());
    }
    for (std::size_t j = 0; j < n; ++j) {
      if (g[j] == groups[gi] && indom[j] && q[j] == "E") caut[gi] = 1;
    }
  }
  for (double v : lev) {
    if (!std::isnan(v) && (v <= 0 || v >= 1)) return stop("domain levels must be percentages strictly between 0 and 100");
  }
  std::uint32_t sd = 1;
  if (!seedarg(c["seed"], sd)) sd = 1;
  if (bad()) return {};
  Rng rg(sd);
  const double a = (1 - cl) / 2;
  std::vector<double> ol(ng), oll(ng), olu(ng), ob(ng), obl(ng), obu(ng);
  std::vector<std::string> ow(ng);
  std::vector<char> own(ng, 0);
  for (std::size_t gi = 0; gi < ng; ++gi) {
    std::vector<double> p(k), s(k);
    bool nap = false;
    for (std::size_t m = 0; m < k; ++m) {
      p[m] = lev[gi * k + m];
      s[m] = sem[gi * k + m];
      nap = nap || std::isnan(p[m]) || std::isnan(s[m]);
    }
    if (nap) {
      ol[gi] = oll[gi] = olu[gi] = ob[gi] = obl[gi] = obu[gi] = NA;
      own[gi] = 1;
      continue;
    }
    std::vector<double> z0(static_cast<std::size_t>(reps) * k);
    for (double &v : z0) v = rg.norm();
    std::vector<double> bb(static_cast<std::size_t>(reps));
    for (int b = 0; b < reps; ++b) {
      double lw = 0, aw = 0;
      for (std::size_t m = 0; m < k; ++m) {
        const double sl = s[m] / (p[m] * (1 - p[m]));
        const double d = plogis(qlogis(p[m]) + z0[m * static_cast<std::size_t>(reps) + static_cast<std::size_t>(b)] * sl);
        lw += std::log(d) * w[m];
        aw += d * w[m];
      }
      bb[static_cast<std::size_t>(b)] = std::exp(lw) / aw;
    }
    long double lv = 0, lse = 0, lg = 0;
    for (std::size_t m = 0; m < k; ++m) {
      lv += static_cast<long double>(w[m] * p[m]);
      lse += static_cast<long double>(w[m] * w[m] * (s[m] * s[m]));
      lg += static_cast<long double>(w[m] * std::log(p[m]));
    }
    const double lvd = static_cast<double>(lv), lsd = std::sqrt(static_cast<double>(lse));
    const double b0 = std::exp(static_cast<double>(lg)) / lvd;
    const double bias = rmean(bb) - b0;
    ol[gi] = 100 * lvd;
    oll[gi] = 100 * std::max(0.0, lvd - z * lsd);
    olu[gi] = 100 * std::min(1.0, lvd + z * lsd);
    ob[gi] = b0;
    obl[gi] = std::min(b0, std::max(0.0, quant7(bb, a) - bias));
    obu[gi] = std::max(b0, std::min(1.0, quant7(bb, 1 - a) - bias));
    std::size_t wm = 0;
    for (std::size_t m = 1; m < k; ++m) {
      if (p[m] < p[wm]) wm = m;
    }
    ow[gi] = doms[wm];
  }
  V out = frame({"group", "level", "lvllo", "lvlup", "balance", "ballo", "balup", "weadom", "usewithcau"},
                {sv(groups), dv(ol), dv(oll), dv(olu), dv(ob), dv(obl), dv(obu), svna(ow, own), lv(caut)});
  V m = dv(std::vector<double>(ng * k));
  for (std::size_t gi = 0; gi < ng; ++gi) {
    for (std::size_t j = 0; j < k; ++j) m.d[j * ng + gi] = 100 * lev[gi * k + j];
  }
  m.nr = static_cast<int>(ng);
  m.nc = static_cast<int>(k);
  m.rn = groups;
  m.cn = doms;
  out.at.push_back({"domains", m});
  out.at.push_back({"weights", named(dv(w), doms)});
  out.at.push_back({"conflvl", dv1(cl)});
  out.cls = {"holbal", "data.frame"};
  return out;
}

V f_dualcont(const Cx &c) {
  const V &x = c["x"];
  static const char *const need[] = {"group", "indicator", "side", "estimate", "lower", "upper", "usewithcau"};
  bool from = inh(x, "welprof");
  if (!from && x.js && isdf(x)) {
    from = true;
    for (const char *n : need) from = from && col(x, n);
  }
  if (!from) return stop("`x` must come from welprof()");
  const std::string well = str1(c["well"], "well"), dis = str1(c["distress"], "distress");
  if (bad()) return {};
  const Lab gl = lab(*col(x, "group")), il = lab(*col(x, "indicator")), sl = lab(*col(x, "side"));
  const std::size_t n = gl.n();
  bool hw = false, hd = false;
  for (std::size_t j = 0; j < n; ++j) {
    hw = hw || (!il.na[j] && !sl.na[j] && il.s[j] == well && sl.s[j] == "wellness");
    hd = hd || (!il.na[j] && !sl.na[j] && il.s[j] == dis && sl.s[j] == "distress");
  }
  if (!hw) return stop("\"%s\" is not a wellness indicator of `x`", well.c_str());
  if (!hd) return stop("\"%s\" is not a distress indicator of `x`", dis.c_str());
  const double cl = conf(c["conflvl"], "conflvl", 0.95);
  if (bad()) return {};
  const double z = qn(cl);
  std::vector<std::string> gs(n);
  for (std::size_t j = 0; j < n; ++j) gs[j] = gl.na[j] ? "\x01NA" : gl.s[j];
  const std::vector<std::string> groups = uniq(gs);
  const std::vector<double> est = num(*col(x, "estimate"), "estimate"), lo = num(*col(x, "lower"), "lower"),
                            hi = num(*col(x, "upper"), "upper");
  if (bad()) return {};
  const V &ucv = *col(x, "usewithcau");
  const double NA = na();
  // the row of indicator ind for each group, or none
  const std::vector<std::size_t> pw = rowsof(il, gs, groups, well), pd = rowsof(il, gs, groups, dis);
  auto get = [&](const std::vector<double> &v, std::size_t j) { return j == SIZE_MAX ? NA : v[j]; };
  auto slog = [&](std::size_t j) {
    const double l = get(lo, j), h = get(hi, j);
    if (std::isnan(l) || std::isnan(h)) return NA;
    return l > 0 && h > 0 ? (std::log(h) - std::log(l)) / (2 * z) : NA;
  };
  auto caut = [&](std::size_t j) {
    if (j == SIZE_MAX) return false;
    const int v = ucv.t == V::LGL ? ucv.i[j] : INT_MIN;
    return v != INT_MIN && v != 0;
  };
  std::vector<std::string> og;
  std::vector<char> ogn;
  std::vector<double> w, wl, wu, d, dl, du, r, rl, ru;
  std::vector<int> uc;
  for (std::size_t k = 0; k < groups.size(); ++k) {
    og.push_back(groups[k] == "\x01NA" ? "" : groups[k]);
    ogn.push_back(groups[k] == "\x01NA");
    const double a = get(est, pw[k]), b = get(est, pd[k]);
    w.push_back(a);
    wl.push_back(get(lo, pw[k]));
    wu.push_back(get(hi, pw[k]));
    d.push_back(b);
    dl.push_back(get(lo, pd[k]));
    du.push_back(get(hi, pd[k]));
    const double ra = a / b;
    const double s1 = slog(pw[k]), s2 = slog(pd[k]);
    const double sr = std::sqrt(s1 * s1 + s2 * s2);
    const bool ok = std::isfinite(ra) && ra > 0;
    r.push_back(ok ? ra : NA);
    rl.push_back(ok ? std::exp(std::log(ra) - z * sr) : NA);
    ru.push_back(ok ? std::exp(std::log(ra) + z * sr) : NA);
    uc.push_back(caut(pw[k]) || caut(pd[k]));
  }
  V out = frame({"group", "wellness", "wello", "welup", "distress", "distrlo", "distrup", "strratio", "ratiolo",
                 "ratioup", "usewithcau"},
                {svna(og, ogn), dv(w), dv(wl), dv(wu), dv(d), dv(dl), dv(du), dv(r), dv(rl), dv(ru), lv(uc)});
  out.at.push_back({"labels", named(sv({well, dis}), {"wellness", "distress"})});
  out.cls = {"dualcont", "data.frame"};
  return out;
}

}  // namespace igh
