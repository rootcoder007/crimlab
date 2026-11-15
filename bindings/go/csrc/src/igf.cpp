// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// First Nations, Inuit and Metis frameworks for wellness put to work on
// published tables: the First Nations Mental Wellness Continuum's four
// directions (mhcont), Indigenous Services Canada's Community Well-Being
// index (cwbscr), connection across places (connassoc), healing as change
// over time (healchg), a community-weighted index of protective and risk
// factors (commidx), and an audit of how deficit-framed a set of results
// is (frmaud). igh.h.
#include <algorithm>
#include <climits>
#include <cmath>
#include <map>
#include <set>

#include "crimlab.h"
#include "igk.h"

namespace igh {

namespace {

// .rmbl_ind_cols(): data a data frame, each argument the name of one of
// its columns; the optional ones may be NULL
void indcols(const Cx &c, const V &data, std::initializer_list<const char *> req,
             std::initializer_list<const char *> opt) {
  if (!isdf(data)) stop("`data` must be a data frame");
  for (const char *nm : req) {
    if (!c.has(nm)) stop("`%s` must name a column of `data`", nm);
    colnm(data, c[nm], nm);
  }
  for (const char *nm : opt) {
    if (c.has(nm)) colnm(data, c[nm], nm);
  }
}

std::string cs(const Cx &c, const char *nm) { return c.has(nm) ? str1(c[nm], nm) : std::string(); }

// group or indicator labels, trimmed; NA as "NA"
std::vector<std::string> tl(const V &v) {
  const Lab l = tlab(v);
  std::vector<std::string> o(l.n());
  for (std::size_t j = 0; j < l.n(); ++j) o[j] = l.na[j] ? "NA" : l.s[j];
  return o;
}

bool in(const std::string &s, const std::vector<std::string> &v) { return std::find(v.begin(), v.end(), s) != v.end(); }

}  // namespace

V f_mhcont(const Cx &c) {
  const V &data = c["data"];
  indcols(c, data, {"grp", "ind", "est", "lower", "upper"}, {"qual"});
  const double cl = conf(c["conflvl"], "conflvl", 0.95);
  const double z = qn(cl);
  static const char *const an[] = {"hope", "belong", "meaning", "purp"};
  static const char *const dn[] = {"hope", "belonging", "meaning", "purpose"};
  for (const char *a : an) nofn(c[a], a);
  std::vector<std::vector<std::string>> dirs(4);
  std::vector<char> isnull(4);
  for (int k = 0; k < 4; ++k) {
    const V &v = c[an[k]];
    isnull[static_cast<std::size_t>(k)] = v.t == V::NUL;
    if (v.t == V::NUL) continue;
    if (!isstr(v)) stop("`%s` must name indicators of `data`, or be NULL", an[k]);
    for (const std::string &s : v.s) dirs[static_cast<std::size_t>(k)].push_back(trim(s));
  }
  if (isnull[0] && isnull[1] && isnull[2] && isnull[3]) stop("name a proxy for one direction at least");
  const std::vector<std::string> ind = tl(*col(data, str1(c["ind"], "ind")));
  std::vector<std::string> all, miss;
  for (const auto &d : dirs) all.insert(all.end(), d.begin(), d.end());
  for (const std::string &s : uniq(all)) {
    if (!in(s, ind)) miss.push_back(s);
  }
  if (!miss.empty()) stop("not an indicator of `data`: %s", tostr(miss).c_str());
  const std::vector<std::string> g = tl(*col(data, str1(c["grp"], "grp")));
  std::vector<double> est = num(*col(data, str1(c["est"], "est")), "est");
  const std::vector<double> lo = num(*col(data, str1(c["lower"], "lower")), "lower");
  const std::vector<double> hi = num(*col(data, str1(c["upper"], "upper")), "upper");
  std::vector<char> named(ind.size());
  std::vector<std::string> labs(ind.size());
  for (std::size_t j = 0; j < ind.size(); ++j) {
    named[j] = in(ind[j], all);
    labs[j] = g[j] + ": " + ind[j];
  }
  contains(c, "mhcont", est, lo, hi, labs, &named);
  std::vector<double> se(est.size());
  for (std::size_t j = 0; j < se.size(); ++j) se[j] = (hi[j] - lo[j]) / (2 * z);
  const std::vector<std::string> q = quality(c, data, cs(c, "qual"), g.size());
  for (std::size_t j = 0; j < est.size(); ++j) {
    if (gone(q[j])) est[j] = na();
  }
  std::map<std::string, int> used;
  for (const auto &d : dirs) {
    for (const std::string &s : uniq(d)) used[s]++;
  }
  std::vector<std::string> og, oi, od;
  std::vector<double> oe, ol, ou;
  std::vector<int> oc, osh;
  const double NA = na();
  for (const std::string &gi : uniq(g)) {
    for (int k = 0; k < 4; ++k) {
      og.push_back(gi);
      oi.push_back(dn[k]);
      if (isnull[static_cast<std::size_t>(k)]) {
        od.push_back("no published proxy");
        oe.push_back(NA);
        ol.push_back(NA);
        ou.push_back(NA);
        oc.push_back(0);
        osh.push_back(0);
        continue;
      }
      const std::vector<std::string> &p = dirs[static_cast<std::size_t>(k)];
      std::vector<std::size_t> ii;
      for (std::size_t j = 0; j < g.size(); ++j) {
        if (g[j] == gi && in(ind[j], p)) ii.push_back(j);
      }
      bool ok = ii.size() == p.size();
      std::vector<double> ev, sv2;
      bool caut = false;
      for (std::size_t j : ii) {
        ok = ok && !std::isnan(est[j]);
        ev.push_back(est[j]);
        sv2.push_back(se[j] * se[j]);
        caut = caut || q[j] == "E";
      }
      const double m = ok ? rmean(ev) : NA;
      const double s = ok ? std::sqrt(rsum(sv2)) / static_cast<double>(ii.size()) : NA;
      std::string dom;
      bool sh = false;
      for (std::size_t j = 0; j < p.size(); ++j) {
        dom += (j ? "; " : "") + p[j];
        sh = sh || used[p[j]] > 1;
      }
      od.push_back(dom);
      oe.push_back(m);
      ol.push_back(m - z * s);
      ou.push_back(m + z * s);
      oc.push_back(caut);
      osh.push_back(sh);
    }
  }
  V out = frame({"group", "indicator", "side", "domain", "estimate", "lower", "upper", "usewithcau", "shared", "proxy"},
                {sv(og), sv(oi), sv(std::vector<std::string>(og.size(), "wellness")), sv(od), dv(oe), dv(ol), dv(ou),
                 lv(oc), lv(osh), lv(std::vector<int>(og.size(), 1))});
  out.at.push_back({"conflvl", dv1(cl)});
  out.cls = {"mhcont", "welprof", "data.frame"};
  return out;
}

V f_cwbscr(const Cx &c) {
  static const char *const an[] = {"incpercap", "highsch20p", "uni25p", "notcro", "nomajrep", "partic", "empl"};
  static const char *const nn[] = {"incpercap", "highsch20p", "uni25p", "notcro", "nomajrep", "participation",
                                   "employment"};
  std::vector<std::vector<double>> a(7);
  for (int k = 0; k < 7; ++k) a[static_cast<std::size_t>(k)] = num(c[an[k]], nn[k]);
  std::size_t n = 0;
  for (const auto &v : a) n = std::max(n, v.size());
  bool bad = n == 0;
  for (const auto &v : a) bad = bad || (v.size() != n && v.size() != 1);
  if (bad) stop("every argument must have one value per community (or one value for all)");
  for (int k = 1; k < 7; ++k) {
    for (double x : a[static_cast<std::size_t>(k)]) {
      if (x < 0 || x > 100) stop("`%s` must be percentages from 0 to 100", nn[k]);
    }
  }
  for (double x : a[0]) {
    if (x <= 0) stop("`incpercap` must be positive");
  }
  const std::vector<double> b = !dflt(c["bounds"]) ? num(c["bounds"], "bounds") : std::vector<double>{2650, 75000};
  if (b.size() != 2 || std::isnan(b[0]) || std::isnan(b[1]) || b[0] <= 0 || b[1] <= b[0]) {
    stop("`bounds` must be two increasing positive incomes");
  }
  auto at = [&](int k, std::size_t j) { return a[static_cast<std::size_t>(k)][j % a[static_cast<std::size_t>(k)].size()]; };
  std::vector<double> inc(n), edu(n), hou(n), lab2(n), cwb(n);
  for (std::size_t j = 0; j < n; ++j) {
    double v = 100 * (std::log(at(0, j)) - std::log(b[0])) / (std::log(b[1]) - std::log(b[0]));
    if (!std::isnan(v)) v = std::min(100.0, std::max(0.0, v));
    inc[j] = v;
    edu[j] = 2.0 / 3 * at(1, j) + 1.0 / 3 * at(2, j);
    hou[j] = (at(3, j) + at(4, j)) / 2;
    lab2[j] = (at(5, j) + at(6, j)) / 2;
    const long double s = static_cast<long double>(inc[j]) + edu[j] + hou[j] + lab2[j];
    cwb[j] = static_cast<double>(s / 4);
  }
  V out = frame({"income", "education", "housing", "labour", "cwb"}, {dv(inc), dv(edu), dv(hou), dv(lab2), dv(cwb)});
  out.cls = {"cwbscr", "data.frame"};
  return out;
}

V f_connassoc(const Cx &c) {
  const V &data = c["data"];
  indcols(c, data, {"place", "x", "y"}, {"xlo", "xup", "ylo", "yup"});
  const int nl = c.has("xlo") + c.has("xup") + c.has("ylo") + c.has("yup");
  if (nl % 4) stop("give all four limit columns, or none");
  const int reps = static_cast<int>(scalar(c["reps"], "reps", 100, 1e6, true, 10000));
  const double cl = conf(c["conflvl"], "conflvl", 0.95);
  const double z = qn(cl);
  const std::string xn = str1(c["x"], "x"), yn = str1(c["y"], "y");
  const std::vector<double> xv = num(*col(data, xn), "x"), yv = num(*col(data, yn), "y");
  std::vector<char> ok(xv.size());
  int nok = 0;
  for (std::size_t j = 0; j < ok.size(); ++j) {
    ok[j] = std::isfinite(xv[j]) && std::isfinite(yv[j]);
    nok += ok[j];
  }
  if (nok < 3) stop("the association needs three places at least with both indicators");
  const Lab pl = lab(*col(data, str1(c["place"], "place")));
  const bool lim = nl == 4;
  std::vector<double> xse(xv.size()), yse(xv.size());
  if (lim) {
    const std::vector<double> xu = num(*col(data, str1(c["xup"], "xup")), "xup");
    const std::vector<double> xl = num(*col(data, str1(c["xlo"], "xlo")), "xlo");
    const std::vector<double> yu = num(*col(data, str1(c["yup"], "yup")), "yup");
    const std::vector<double> yl = num(*col(data, str1(c["ylo"], "ylo")), "ylo");
    for (std::size_t j = 0; j < ok.size(); ++j) {
      xse[j] = (xu[j] - xl[j]) / (2 * z);
      yse[j] = (yu[j] - yl[j]) / (2 * z);
      ok[j] = ok[j] && std::isfinite(xse[j]) && std::isfinite(yse[j]);
    }
  }
  std::vector<std::size_t> keep;
  std::vector<std::string> rn;
  for (std::size_t j = 0; j < ok.size(); ++j) {
    if (ok[j]) {
      keep.push_back(j);
      rn.push_back(std::to_string(j + 1));
    }
  }
  std::vector<double> x, y, xs, ys;
  std::vector<std::string> ps;
  std::vector<char> pna;
  for (std::size_t j : keep) {
    x.push_back(xv[j]);
    y.push_back(yv[j]);
    xs.push_back(xse[j]);
    ys.push_back(yse[j]);
    ps.push_back(pl.s[j]);
    pna.push_back(pl.na[j]);
  }
  const std::size_t n = x.size();
  const std::vector<double> rx = rank(x), ry = rank(y);
  const double r0 = cor(rx, ry);
  std::uint32_t sd = 1;
  if (!seedarg(c["seed"], sd)) sd = 1;
  Rng g(sd);
  const bool exact = n <= 8;
  double hits = 0, p;
  const double NA = na();
  if (exact && std::isnan(r0)) {
    p = NA;
  } else if (exact) {
    std::vector<std::size_t> pm(n);
    for (std::size_t j = 0; j < n; ++j) pm[j] = j;
    double cnt = 0;
    std::vector<double> yp(n);
    do {
      for (std::size_t j = 0; j < n; ++j) yp[j] = ry[pm[j]];
      const double r = cor(rx, yp);
      if (std::fabs(r) >= std::fabs(r0) - 1e-12) hits += 1;
      cnt += 1;
    } while (std::next_permutation(pm.begin(), pm.end()));
    p = hits / cnt;
  } else {
    std::vector<int> pool(n);
    std::vector<double> yp(n);
    bool nan = false;
    for (int b = 0; b < reps; ++b) {
      for (std::size_t j = 0; j < n; ++j) pool[j] = static_cast<int>(j);
      int left = static_cast<int>(n);
      for (std::size_t j = 0; j < n; ++j) {
        const int k = static_cast<int>(g.index(left));
        yp[j] = ry[static_cast<std::size_t>(pool[static_cast<std::size_t>(k)])];
        pool[static_cast<std::size_t>(k)] = pool[static_cast<std::size_t>(--left)];
      }
      const double r = cor(rx, yp);
      if (std::isnan(r)) nan = true;
      else if (std::fabs(r) >= std::fabs(r0) - 1e-12) hits += 1;
    }
    p = nan || std::isnan(r0) ? NA : (hits + 1) / (reps + 1);
  }
  double lo = NA, hi = NA;
  if (lim) {
    std::vector<double> draws(static_cast<std::size_t>(reps));
    std::vector<double> dx(n), dy(n);
    for (int b = 0; b < reps; ++b) {
      // cor() reads its second argument first: the y draws come first
      for (std::size_t j = 0; j < n; ++j) dy[j] = ys[j] == 0 ? y[j] : y[j] + ys[j] * g.norm();
      for (std::size_t j = 0; j < n; ++j) dx[j] = xs[j] == 0 ? x[j] : x[j] + xs[j] * g.norm();
      draws[static_cast<std::size_t>(b)] = cor(rank(dx), rank(dy));
    }
    const double a = (1 - cl) / 2;
    lo = quant7(draws, a);
    hi = quant7(draws, 1 - a);
  }
  std::vector<std::string> dnm{"place", "x", "y"};
  std::vector<V> dc{svna(ps, pna), dv(x), dv(y)};
  if (lim) {
    dnm.push_back("xse");
    dnm.push_back("yse");
    dc.push_back(dv(xs));
    dc.push_back(dv(ys));
  }
  V d = frame(dnm, dc);
  if (keep.size() != ok.size()) d.rn = rn;
  V out = nlist({"data", "rho", "p_value", "exact", "rholo", "rhoup", "n", "labels", "conflvl"},
                {d, dv1(r0), dv1(p), lv({exact}), dv1(lo), dv1(hi), iv({static_cast<int>(n)}),
                 named(sv({xn, yn}), {"x", "y"}), dv1(cl)});
  out.cls = {"connassoc"};
  return out;
}

V f_healchg(const Cx &c) {
  const V &data = c["data"];
  indcols(c, data, {"grp", "ind", "prd", "est", "lower", "upper"}, {"qual"});
  const std::string from = str1(c["from"], "from"), to = str1(c["to"], "to");
  const V &well = c["well"], &dis = c["distress"];
  nofn(well, "well");
  nofn(dis, "distress");
  const double cl = conf(c["conflvl"], "conflvl", 0.95);
  const double z = qn(cl);
  const Lab prl = lab(*col(data, str1(c["prd"], "prd")));
  std::vector<std::string> pr(prl.n());
  for (std::size_t j = 0; j < pr.size(); ++j) pr[j] = prl.na[j] ? "\x01NA" : prl.s[j];
  for (const std::string &v : {from, to}) {
    if (!in(v, pr)) stop("\"%s\" is not a period of `data`", v.c_str());
  }
  const std::vector<std::string> g = tl(*col(data, str1(c["grp"], "grp")));
  const std::vector<std::string> ind = tl(*col(data, str1(c["ind"], "ind")));
  std::vector<double> est = num(*col(data, str1(c["est"], "est")), "est");
  const std::vector<double> lo = num(*col(data, str1(c["lower"], "lower")), "lower");
  const std::vector<double> hi = num(*col(data, str1(c["upper"], "upper")), "upper");
  std::vector<std::string> labs(g.size());
  for (std::size_t j = 0; j < g.size(); ++j) labs[j] = g[j] + ": " + ind[j] + " (" + (prl.na[j] ? "NA" : pr[j]) + ")";
  contains(c, "healchg", est, lo, hi, labs, nullptr);
  std::vector<double> se(est.size());
  for (std::size_t j = 0; j < se.size(); ++j) se[j] = (hi[j] - lo[j]) / (2 * z);
  const std::vector<std::string> q = quality(c, data, cs(c, "qual"), g.size());
  for (std::size_t j = 0; j < est.size(); ++j) {
    if (gone(q[j])) est[j] = na();
  }
  std::vector<std::string> key(g.size());
  for (std::size_t j = 0; j < g.size(); ++j) key[j] = g[j] + "\r" + ind[j];
  std::vector<std::size_t> k1, k2;
  for (std::size_t j = 0; j < pr.size(); ++j) {
    if (pr[j] == from) k1.push_back(j);
    if (pr[j] == to) k2.push_back(j);
  }
  std::map<std::string, std::size_t> m1, m2;
  bool dup = false;
  std::vector<std::string> o1;
  for (std::size_t j : k1) {
    dup = dup || m1.count(key[j]);
    if (!m1.count(key[j])) {
      m1[key[j]] = j;
      o1.push_back(key[j]);
    }
  }
  for (std::size_t j : k2) {
    dup = dup || m2.count(key[j]);
    if (!m2.count(key[j])) m2[key[j]] = j;
  }
  if (dup) stop("a group and indicator appear twice in one period: filter `data` to one age group, sex and place");
  std::vector<std::size_t> i1, i2;
  for (const std::string &k : o1) {
    if (m2.count(k)) {
      i1.push_back(m1[k]);
      i2.push_back(m2[k]);
    }
  }
  if (i1.empty()) stop("no group and indicator is in both periods");
  std::vector<std::string> wl, dl;
  const bool sides = well.t != V::NUL || dis.t != V::NUL;
  if (well.t != V::NUL) {
    const Lab l = lab(well);
    for (std::size_t j = 0; j < l.n(); ++j) wl.push_back(l.na[j] ? "\x01" : trim(l.s[j]));
  }
  if (dis.t != V::NUL) {
    const Lab l = lab(dis);
    for (std::size_t j = 0; j < l.n(); ++j) dl.push_back(l.na[j] ? "\x01" : trim(l.s[j]));
  }
  struct Row {
    std::string g, i, side;
    bool sna;
    double f, t, ch, lo, hi;
    int uc;
  };
  std::vector<Row> rows;
  for (std::size_t k = 0; k < i1.size(); ++k) {
    const std::size_t a = i1[k], b = i2[k];
    Row r;
    r.g = g[a];
    r.i = ind[a];
    r.sna = !sides;
    if (sides) r.side = in(ind[a], wl) ? "wellness" : dis.t == V::NUL ? "distress" : in(ind[a], dl) ? "distress" : "other";
    r.f = est[a];
    r.t = est[b];
    r.ch = est[b] - est[a];
    const double s = std::sqrt(se[a] * se[a] + se[b] * se[b]);
    r.lo = r.ch - z * s;
    r.hi = r.ch + z * s;
    r.uc = q[a] == "E" || q[b] == "E";
    rows.push_back(r);
  }
  std::vector<std::string> ui;
  for (const Row &r : rows) ui.push_back(r.i);
  ui = uniq(ui);
  std::vector<std::size_t> o(rows.size());
  for (std::size_t j = 0; j < o.size(); ++j) o[j] = j;
  auto rk = [&](std::size_t j) {
    const int a = rows[j].sna || rows[j].side != "wellness";
    const int b = static_cast<int>(std::find(ui.begin(), ui.end(), rows[j].i) - ui.begin());
    return std::make_pair(a, b);
  };
  std::stable_sort(o.begin(), o.end(), [&](std::size_t a, std::size_t b) { return rk(a) < rk(b); });
  std::vector<std::string> og, oi, os;
  std::vector<char> osn;
  std::vector<double> of, ot, oc, ol, oh;
  std::vector<int> ou;
  for (std::size_t j : o) {
    const Row &r = rows[j];
    og.push_back(r.g);
    oi.push_back(r.i);
    os.push_back(r.side);
    osn.push_back(r.sna);
    of.push_back(r.f);
    ot.push_back(r.t);
    oc.push_back(r.ch);
    ol.push_back(r.lo);
    oh.push_back(r.hi);
    ou.push_back(r.uc);
  }
  V out = frame({"group", "indicator", "side", "from", "to", "change", "lower", "upper", "usewithcau"},
                {sv(og), sv(oi), svna(os, osn), dv(of), dv(ot), dv(oc), dv(ol), dv(oh), lv(ou)});
  out.at.push_back({"periods", named(sv({from, to}), {"from", "to"})});
  out.cls = {"healchg", "data.frame"};
  return out;
}

V f_commidx(const Cx &c) {
  const V &data = c["data"];
  indcols(c, data, {"grp", "ind", "est", "lower", "upper"}, {"qual"});
  V w = c["w"];
  nofn(w, "w");
  // a JSON object of numbers is a named vector
  if (w.t == V::LIST && w.js) {
    V d;
    d.t = V::DBL;
    bool okj = true;
    for (const V &e : w.a) {
      okj = okj && e.t == V::DBL && e.d.size() == 1;
      if (okj) d.d.push_back(e.d[0]);
    }
    if (okj) {
      d.nm = w.nm;
      w = d;
    }
  }
  bool bad = !(w.t == V::DBL || w.t == V::INT) || w.nm.empty();
  std::vector<double> wv;
  if (!bad) {
    wv = num(w, "w");
    bool all0 = true;
    for (double x : wv) {
      bad = bad || std::isnan(x);
      all0 = all0 && x == 0;
    }
    bad = bad || all0 || uniq(w.nm).size() != w.nm.size();
  }
  if (bad) {
    stop("`w` must be numbers named by indicator (positive protective, negative risk), set by the community: "
         "the function does not choose them");
  }
  const double cl = conf(c["conflvl"], "conflvl", 0.95);
  const double z = qn(cl);
  const std::vector<std::string> ind = tl(*col(data, str1(c["ind"], "ind")));
  std::vector<std::string> miss;
  for (const std::string &s : w.nm) {
    if (!in(s, ind)) miss.push_back(s);
  }
  if (!miss.empty()) stop("not an indicator of `data`: %s", tostr(uniq(miss)).c_str());
  const std::vector<std::string> g = tl(*col(data, str1(c["grp"], "grp")));
  std::vector<double> est = num(*col(data, str1(c["est"], "est")), "est");
  const std::vector<double> lo = num(*col(data, str1(c["lower"], "lower")), "lower");
  const std::vector<double> hi = num(*col(data, str1(c["upper"], "upper")), "upper");
  std::vector<char> used(ind.size());
  std::vector<std::string> labs(ind.size());
  for (std::size_t j = 0; j < ind.size(); ++j) {
    used[j] = in(ind[j], w.nm);
    labs[j] = g[j] + ": " + ind[j];
  }
  contains(c, "commidx", est, lo, hi, labs, &used);
  std::vector<double> se(est.size());
  for (std::size_t j = 0; j < se.size(); ++j) se[j] = (hi[j] - lo[j]) / (2 * z);
  const std::vector<std::string> q = quality(c, data, cs(c, "qual"), g.size());
  for (std::size_t j = 0; j < est.size(); ++j) {
    if (gone(q[j])) est[j] = na();
  }
  std::vector<double> aw;
  for (double x : wv) aw.push_back(std::fabs(x));
  const double tot = rsum(aw);
  std::vector<std::string> og;
  std::vector<double> oi, ol, ou;
  std::vector<int> oc;
  const double NA = na();
  for (const std::string &gi : uniq(g)) {
    std::vector<std::size_t> ii;
    std::set<std::string> seen;
    for (std::size_t j = 0; j < g.size(); ++j) {
      if (g[j] == gi && used[j]) {
        if (!seen.insert(ind[j]).second) stop("group \"%s\" has an indicator twice", gi.c_str());
        ii.push_back(j);
      }
    }
    og.push_back(gi);
    std::vector<double> num1, var1;
    bool ok = true, caut = false;
    for (std::size_t k = 0; k < w.nm.size(); ++k) {
      std::size_t at = SIZE_MAX;
      for (std::size_t j : ii) {
        if (ind[j] == w.nm[k]) at = j;
      }
      if (at == SIZE_MAX || std::isnan(est[at])) {
        ok = false;
        break;
      }
      num1.push_back(wv[k] * est[at]);
      var1.push_back(wv[k] * wv[k] * se[at] * se[at]);
      caut = caut || q[at] == "E";
    }
    if (!ok) {
      oi.push_back(NA);
      ol.push_back(NA);
      ou.push_back(NA);
      oc.push_back(0);
      continue;
    }
    const double v = rsum(num1) / tot, s = std::sqrt(rsum(var1)) / tot;
    oi.push_back(v);
    ol.push_back(v - z * s);
    ou.push_back(v + z * s);
    oc.push_back(caut);
  }
  V out = frame({"group", "index", "lower", "upper", "usewithcau"}, {sv(og), dv(oi), dv(ol), dv(ou), lv(oc)});
  V wa = w;
  wa.js = false;
  out.at.push_back({"weights", wa});
  out.at.push_back({"conflvl", dv1(cl)});
  out.cls = {"commidx", "data.frame"};
  return out;
}

V f_frmaud(const Cx &c) {
  const V &x = c["x"];
  const bool from = inh(x, "welprof") || (x.js && isdf(x) && col(x, "group") && col(x, "side"));
  if (!from) stop("`x` must come from welprof()");
  const Lab g = lab(*col(x, "group")), s = lab(*col(x, "side"));
  bool any = false;
  for (std::size_t j = 0; j < s.n(); ++j) any = any || (!s.na[j] && (s.s[j] == "wellness" || s.s[j] == "distress"));
  if (!any) stop("`x` has no wellness or distress indicator: there is no framing to count");
  std::vector<std::string> gs(g.n());
  for (std::size_t j = 0; j < g.n(); ++j) gs[j] = g.na[j] ? "\x01NA" : g.s[j];
  const std::vector<std::string> gr = uniq(gs);
  std::vector<int> wc, dc;
  std::vector<double> sh;
  std::vector<std::string> df;
  int tw = 0, td = 0;
  for (const std::string &gi : gr) {
    int a = 0, b = 0;
    for (std::size_t j = 0; j < gs.size(); ++j) {
      if (gs[j] != gi || s.na[j]) continue;
      a += s.s[j] == "wellness";
      b += s.s[j] == "distress";
    }
    wc.push_back(a);
    dc.push_back(b);
    tw += a;
    td += b;
    const int t = a + b;
    sh.push_back(t > 0 ? static_cast<double>(b) / std::max(t, 1) : na());
    if (a == 0) df.push_back(gi);
  }
  std::vector<std::string> gl;
  std::vector<char> gn;
  for (const std::string &gi : gr) {
    gl.push_back(gi == "\x01NA" ? "" : gi);
    gn.push_back(gi == "\x01NA");
  }
  std::vector<char> dn;
  for (std::string &d : df) {
    dn.push_back(d == "\x01NA");
    if (dn.back()) d.clear();
  }
  V bg = frame({"group", "wellness", "distress", "shrdef"}, {svna(gl, gn), iv(wc), iv(dc), dv(sh)});
  V out = nlist({"bygrp", "shrdef", "defonly"},
                {bg, dv1(static_cast<double>(td) / (tw + td)), svna(df, dn)});
  out.cls = {"frmaud"};
  return out;
}

}  // namespace igh
