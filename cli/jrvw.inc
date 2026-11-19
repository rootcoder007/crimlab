// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Reviews, entitlements and decision makers (R/rvw.R): whether a review
// due by a set day was held, on how many days of a stay an entitlement was
// missed, and whether decision makers decide alike. The checks, their
// order and their messages are those of the R functions.

namespace rvw {

// numbers with NA for null; false when the argument is not numeric
inline bool nums(const Val &o, std::vector<double> &v) {
  if (!ax::numv(o)) return false;
  v = plt::vals(o);
  return true;
}
// numbers, or nothing but missing values
inline bool numsna(const Val &o, std::vector<double> &v) {
  if (nums(o, v)) return true;
  const std::vector<Val> e = plt::elems(o);
  for (const Val &x : e) {
    if (x.t != Val::NUL) return false;
  }
  v.assign(e.size(), NA);
  return true;
}
inline bool anyna(const std::vector<double> &v) {
  for (double x : v) {
    if (std::isnan(x)) return true;
  }
  return false;
}
// .sdnum(): numbers from 0, none missing
inline bool from0(Ctx &c, const Val &o, const char *nm, std::vector<double> &v) {
  bool ok = nums(o, v) && !v.empty() && !anyna(v);
  for (std::size_t i = 0; ok && i < v.size(); ++i) ok = v[i] >= 0;
  if (!ok) return c.fail("`%s` must be numbers from 0 with no missing value.", nm);
  return true;
}
// .sdpct(): percentages from 0 to 100, none missing
inline bool pct(Ctx &c, const Val &o, const char *nm, std::vector<double> &v) {
  bool ok = nums(o, v) && !v.empty() && !anyna(v);
  for (std::size_t i = 0; ok && i < v.size(); ++i) ok = v[i] >= 0 && v[i] <= 100;
  if (!ok) return c.fail("`%s` must be percentages from 0 to 100 with no missing value.", nm);
  return true;
}
// .sdcuts(): increasing numbers above 0 and below `top`
inline bool cuts(Ctx &c, const Val &o, double top, std::vector<double> &v) {
  bool ok = nums(o, v) && !v.empty() && !anyna(v);
  for (std::size_t i = 0; ok && i < v.size(); ++i) ok = v[i] > 0 && v[i] < top && (i == 0 || v[i] > v[i - 1]);
  if (!ok) {
    return std::isfinite(top) ? c.fail("`cuts` must be increasing numbers above 0 and below %.15g.", top)
                              : c.fail("`cuts` must be increasing numbers above 0.");
  }
  return true;
}
// .sdlev(): one number between 0 and 1
inline bool level(Ctx &c, double &lv) {
  const Val &o = c["level"];
  lv = 0.95;
  if (o.t == Val::NUL) return true;
  std::vector<double> v;
  if (!nums(o, v) || v.size() != 1 || !std::isfinite(v[0]) || v[0] <= 0 || v[0] >= 1) {
    return c.fail("`level` must be one number between 0 and 1.");
  }
  lv = v[0];
  return true;
}
// labels, one per row and none missing; false otherwise
inline bool labels(const Val &o, std::size_t n, ax::Cols &l) {
  const std::vector<Val> e = plt::elems(o);
  if (e.size() != n) return false;
  for (const Val &x : e) {
    if (ax::vna(x) || x.t == Val::ARR || x.t == Val::OBJ) return false;
    l.push_back(ax::vlab(x));
  }
  return true;
}
// .sdgrp(): the groups in the order they appear, then "all"; code[i] the
// group of row i
struct Grp {
  ax::Cols lev;
  std::vector<std::size_t> code;
  bool in(std::size_t k, std::size_t i) const { return k + 1 == lev.size() || code[i] == k; }
};
inline bool groups(Ctx &c, std::size_t n, Grp &g) {
  const Val &o = c["by"];
  g.code.assign(n, 0);
  if (o.t != Val::NUL) {
    ax::Cols l;
    if (!labels(o, n, l)) return c.fail("`by` must have one label per stay and no missing value.");
    for (std::size_t i = 0; i < n; ++i) {
      std::size_t k = 0;
      while (k < g.lev.size() && g.lev[k] != l[i]) ++k;
      if (k == g.lev.size()) g.lev.push_back(l[i]);
      g.code[i] = k;
    }
  }
  g.lev.push_back("all");
  return true;
}
inline std::string num15(double x) {
  char b[40];
  std::snprintf(b, sizeof b, "%.15g", x);
  return b;
}
// .sdband(): the names of the bands the cuts part
inline ax::Cols bands(const std::vector<double> &ct, bool top) {
  ax::Cols r;
  r.push_back((top ? "over 0 to " : "up to ") + num15(ct[0]));
  for (std::size_t j = 1; j < ct.size(); ++j) r.push_back("over " + num15(ct[j - 1]) + " to " + num15(ct[j]));
  r.push_back("over " + num15(ct.back()) + (top ? " to under 100" : ""));
  return r;
}
inline std::size_t above(const std::vector<double> &ct, double x) {
  std::size_t k = 0;
  for (double v : ct) k += x > v;
  return k;
}
// .sdcp(): k of n, its share and exact limits (NA of no cases)
inline void cp(double k, double n, double a, double *o) {
  o[0] = n > 0 ? k / n : NA;
  o[1] = n > 0 ? cplo(k, n, a) : NA;
  o[2] = n > 0 ? cphi(k, n, a) : NA;
}
// .dmchi(): Pearson's chi-square that k of n is one rate
inline void chi(const std::vector<double> &n, const std::vector<double> &k, double *o) {
  const std::size_t m = n.size();
  long double tot = 0.0L, tk = 0.0L;
  for (std::size_t j = 0; j < m; ++j) {
    tot += n[j];
    tk += k[j];
  }
  o[0] = static_cast<double>(m);
  o[1] = o[2] = o[3] = o[4] = NA;
  if (m < 2) return;
  const double p = static_cast<double>(tk) / static_cast<double>(tot);
  if (p == 0 || p == 1) return;
  long double st = 0.0L;
  double low = 0;
  for (int side = 0; side < 2; ++side) {
    for (std::size_t j = 0; j < m; ++j) {
      const double e = side == 0 ? n[j] * p : n[j] * (1 - p), ob = side == 0 ? k[j] : n[j] - k[j];
      st += (ob - e) * (ob - e) / e;
      low += e < 5;
    }
  }
  o[1] = static_cast<double>(st);
  o[2] = static_cast<double>(m) - 1;
  o[3] = clab_pchisq(o[1], o[2], 0);
  o[4] = low;
}

}  // namespace rvw

// rvwtim()
bool rvwtim(Ctx &c, Val &out) {
  std::vector<double> days, first, dv, ct;
  CK(rvw::from0(c, c["days"], "days", days));
  const std::size_t n = days.size();
  bool ok = rvw::numsna(c["first"], first) && first.size() == n;
  for (std::size_t i = 0; ok && i < n; ++i) ok = std::isnan(first[i]) || first[i] >= 0;
  if (!ok) return c.fail("`first` must be numbers from 0 (NA: never reviewed), one per stay.");
  if (!rvw::nums(c["due"], dv) || dv.size() != 1 || !std::isfinite(dv[0]) || dv[0] <= 0) {
    return c.fail("`due` must be one number above 0.");
  }
  const double due = dv[0];
  if (c["cuts"].t == Val::NUL) ct.assign(1, due);
  else CK(rvw::cuts(c, c["cuts"], Inf, ct));
  rvw::Grp g;
  CK(rvw::groups(c, n, g));
  double lv = 0.0;
  CK(rvw::level(c, lv));
  const double a = 1 - lv;
  const ax::Cols lab = rvw::bands(ct, false);
  const std::size_t nb = lab.size();
  ax::Cols cg, cb, sg;
  std::vector<double> cn, cr, cv, cs, cl, ch, sn, sr, sp, spn, ss, sl, sh, slate, smed;
  for (std::size_t k = 0; k < g.lev.size(); ++k) {
    std::vector<double> nn(nb, 0.0), rv(nb, 0.0), fr;
    double tn = 0, tr = 0, past = 0, pn = 0, late = 0;
    for (std::size_t i = 0; i < n; ++i) {
      if (!g.in(k, i)) continue;
      const std::size_t b = rvw::above(ct, days[i]);
      const bool r = !std::isnan(first[i]);
      nn[b] += 1;
      rv[b] += r;
      tn += 1;
      tr += r;
      past += days[i] > due;
      pn += days[i] > due && !r;
      late += r && first[i] > due;
      if (r) fr.push_back(first[i]);
    }
    for (std::size_t b = 0; b < nb; ++b) {
      double o[3];
      rvw::cp(nn[b] - rv[b], nn[b], a, o);
      cg.push_back(g.lev[k]);
      cb.push_back(lab[b]);
      cn.push_back(nn[b]);
      cr.push_back(rv[b]);
      cv.push_back(nn[b] - rv[b]);
      cs.push_back(o[0]);
      cl.push_back(o[1]);
      ch.push_back(o[2]);
    }
    double o[3];
    rvw::cp(pn, past, a, o);
    std::sort(fr.begin(), fr.end());
    const std::size_t m = fr.size();
    sg.push_back(g.lev[k]);
    sn.push_back(tn);
    sr.push_back(tr);
    sp.push_back(past);
    spn.push_back(pn);
    ss.push_back(o[0]);
    sl.push_back(o[1]);
    sh.push_back(o[2]);
    slate.push_back(late);
    smed.push_back(m == 0 ? NA : m % 2 ? fr[m / 2] : (fr[m / 2 - 1] + fr[m / 2]) / 2);
  }
  Frame f, s;
  f.str("group", cg);
  f.str("band", cb);
  f.num("n", cn);
  f.num("rev", cr);
  f.num("never", cv);
  f.num("share", cs);
  f.num("lo", cl);
  f.num("hi", ch);
  s.str("group", sg);
  s.num("n", sn);
  s.num("rev", sr);
  s.num("past", sp);
  s.num("pastnever", spn);
  s.num("share", ss);
  s.num("lo", sl);
  s.num("hi", sh);
  s.num("late", slate);
  s.num("med", smed);
  out = vobj();
  put(out, "cov", f.o);
  put(out, "sum", s.o);
  put(out, "due", vnum(due));
  put(out, "level", vnum(lv));
  return true;
}

// entcmp()
bool entcmp(Ctx &c, Val &out) {
  std::vector<double> miss, days, ct, hours, ref;
  CK(rvw::pct(c, c["miss"], "miss", miss));
  const std::size_t n = miss.size();
  const bool hd = c["days"].t != Val::NUL, hh = c["hours"].t != Val::NUL, hr = c["ref"].t != Val::NUL;
  if (hd) {
    CK(rvw::from0(c, c["days"], "days", days));
    if (days.size() != n) return c.fail("`days` must have one value per stay.");
  }
  if (c["cuts"].t == Val::NUL) ct = {20, 50};
  else CK(rvw::cuts(c, c["cuts"], 100, ct));
  rvw::Grp g;
  CK(rvw::groups(c, n, g));
  double lv = 0.0;
  CK(rvw::level(c, lv));
  const double a = 1 - lv;
  ax::Cols lab = rvw::bands(ct, true);
  lab.insert(lab.begin(), "0");
  lab.push_back("100");
  const std::size_t nb = lab.size();
  ax::Cols tg, tb, sg;
  std::vector<double> tn, ts, sn, s0, s1, ss, sl, sh, shalf, sday;
  for (std::size_t k = 0; k < g.lev.size(); ++k) {
    std::vector<double> nn(nb, 0.0);
    double tot = 0, none = 0, all = 0, half = 0;
    long double md = 0.0L, dd = 0.0L;
    for (std::size_t i = 0; i < n; ++i) {
      if (!g.in(k, i)) continue;
      const std::size_t b = miss[i] == 0 ? 0 : miss[i] == 100 ? nb - 1 : 1 + rvw::above(ct, miss[i]);
      nn[b] += 1;
      tot += 1;
      none += miss[i] == 0;
      all += miss[i] == 100;
      half += miss[i] >= 50;
      if (hd) {
        md += miss[i] * days[i];
        dd += days[i];
      }
    }
    for (std::size_t b = 0; b < nb; ++b) {
      tg.push_back(g.lev[k]);
      tb.push_back(lab[b]);
      tn.push_back(nn[b]);
      ts.push_back(nn[b] / tot);
    }
    double o[3];
    rvw::cp(all, tot, a, o);
    sg.push_back(g.lev[k]);
    sn.push_back(tot);
    s0.push_back(none);
    s1.push_back(all);
    ss.push_back(o[0]);
    sl.push_back(o[1]);
    sh.push_back(o[2]);
    shalf.push_back(half);
    sday.push_back(hd ? static_cast<double>(md) / (100 * static_cast<double>(dd)) : NA);
  }
  Frame t, s;
  t.str("group", tg);
  t.str("band", tb);
  t.num("n", tn);
  t.num("share", ts);
  s.str("group", sg);
  s.num("n", sn);
  s.num("none", s0);
  s.num("all", s1);
  s.num("share", ss);
  s.num("lo", sl);
  s.num("hi", sh);
  s.num("half", shalf);
  s.num("dayshr", sday);
  out = vobj();
  put(out, "tab", t.o);
  put(out, "sum", s.o);
  if (hh) {
    bool ok = rvw::numsna(c["hours"], hours) && hours.size() == n;
    for (std::size_t i = 0; ok && i < n; ++i) ok = std::isnan(hours[i]) || hours[i] >= 0;
    if (!ok) return c.fail("`hours` must be numbers from 0 (NA: none missed), one per stay.");
    const std::vector<double> hc = {0.5, 1, 2, 3};
    std::vector<double> nn(5, 0.0), sh2;
    double tot = 0;
    for (std::size_t i = 0; i < n; ++i) {
      if (!(miss[i] > 0) || std::isnan(hours[i])) continue;
      nn[rvw::above(hc, hours[i])] += 1;
      tot += 1;
    }
    for (double v : nn) sh2.push_back(v / tot);
    Frame h;
    h.str("band", rvw::bands(hc, false));
    h.num("n", nn);
    h.num("share", sh2);
    put(out, "hrs", h.o);
  }
  if (hr) {
    CK(rvw::pct(c, c["ref"], "ref", ref));
    if (ref.size() != n) return c.fail("`ref` must have one value per stay.");
    ax::Cols rg;
    std::vector<double> rn, rm, rt;
    for (std::size_t k = 0; k < g.lev.size(); ++k) {
      double tot = 0, most = 0, net = 0;
      for (std::size_t i = 0; i < n; ++i) {
        if (!g.in(k, i)) continue;
        tot += 1;
        most += ref[i] > 50;
        net += miss[i] > ref[i];
      }
      rg.push_back(g.lev[k]);
      rn.push_back(tot);
      rm.push_back(most);
      rt.push_back(net);
    }
    Frame r;
    r.str("group", rg);
    r.num("n", rn);
    r.num("most", rm);
    r.num("net", rt);
    put(out, "ref", r.o);
  }
  put(out, "level", vnum(lv));
  return true;
}

// dmvar()
bool dmvar(Ctx &c, Val &out) {
  ax::Cols dec, dm, lev;
  const std::size_t n = plt::len(c["dec"]);
  if (n == 0 || !rvw::labels(c["dec"], n, dec)) return c.fail("`dec` must be labels with no missing value.");
  if (!rvw::labels(c["dm"], n, dm)) return c.fail("`dm` must have one label per case and no missing value.");
  const std::vector<Val> oe = plt::elems(c["out"]);
  const bool one = oe.size() == 1 && !ax::vna(oe[0]) && oe[0].t != Val::ARR && oe[0].t != Val::OBJ;
  const std::string ov = one ? ax::vlab(oe[0]) : "";
  if (!one || !ax::has(dec, ov)) return c.fail("`out` must be one of the decisions in `dec`.");
  std::vector<std::size_t> code(n);
  for (std::size_t i = 0; i < n; ++i) {
    std::size_t k = 0;
    while (k < lev.size() && lev[k] != dm[i]) ++k;
    if (k == lev.size()) lev.push_back(dm[i]);
    code[i] = k;
  }
  if (lev.size() < 2) return c.fail("`dm` must name two or more decision makers.");
  double nmin = 5;
  if (c["nmin"].t != Val::NUL) {
    std::vector<double> v;
    if (!rvw::nums(c["nmin"], v) || v.size() != 1 || !std::isfinite(v[0]) || v[0] < 1 || v[0] != std::round(v[0])) {
      return c.fail("`nmin` must be one whole number from 1.");
    }
    nmin = v[0];
  }
  double lv = 0.0;
  CK(rvw::level(c, lv));
  const double a = 1 - lv;
  const std::size_t m = lev.size();
  std::vector<double> nn(m, 0.0), kk(m, 0.0), rate, lo, hi, share;
  double tk = 0;
  for (std::size_t i = 0; i < n; ++i) {
    nn[code[i]] += 1;
    if (dec[i] == ov) {
      kk[code[i]] += 1;
      tk += 1;
    }
  }
  for (std::size_t j = 0; j < m; ++j) {
    double o[3];
    rvw::cp(kk[j], nn[j], a, o);
    rate.push_back(o[0]);
    lo.push_back(o[1]);
    hi.push_back(o[2]);
    share.push_back(kk[j] / tk);
  }
  Frame t, ts, d;
  t.str("dm", lev);
  t.num("n", nn);
  t.num("k", kk);
  t.num("rate", rate);
  t.num("lo", lo);
  t.num("hi", hi);
  t.num("share", share);
  std::vector<double> bn, bk;
  for (std::size_t j = 0; j < m; ++j) {
    if (nn[j] >= nmin) {
      bn.push_back(nn[j]);
      bk.push_back(kk[j]);
    }
  }
  double r1[5], r2[5];
  rvw::chi(nn, kk, r1);
  rvw::chi(bn, bk, r2);
  ts.str("set", {"all", "n >= nmin"});
  ts.num("m", {r1[0], r2[0]});
  ts.num("stat", {r1[1], r2[1]});
  ts.num("df", {r1[2], r2[2]});
  ts.num("p", {r1[3], r2[3]});
  ts.num("low", {r1[4], r2[4]});
  std::vector<double> ds, df, dp;
  for (std::size_t i = 0; i < m; ++i) {
    std::vector<double> xn, xk;
    for (std::size_t j = 0; j < m; ++j) {
      if (j != i) {
        xn.push_back(nn[j]);
        xk.push_back(kk[j]);
      }
    }
    double r[5];
    rvw::chi(xn, xk, r);
    ds.push_back(r[1]);
    df.push_back(r[2]);
    dp.push_back(r[3]);
  }
  d.str("dm", lev);
  d.num("stat", ds);
  d.num("df", df);
  d.num("p", dp);
  out = vobj();
  put(out, "tab", t.o);
  put(out, "test", ts.o);
  put(out, "drop", d.o);
  put(out, "out", vstr(ov));
  put(out, "level", vnum(lv));
  return true;
}
