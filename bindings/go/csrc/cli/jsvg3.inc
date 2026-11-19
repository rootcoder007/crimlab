// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// plotsvg(), third part: the plots of the Indigenous indicators
// (R/ploti.R), of the exposure models (R/apxplt.R), of the maps of points
// and of the forensic curves. Included by jsvg.inc inside namespace sv.

// ---- the Indigenous indicators: tables with one row per group -----------

// .pstar(): a label, starred where the source says to use it with caution.
inline ax::Cols star(const Val &d, const char *lab = "group") {
  ax::Cols l = sv_(d, lab);
  const std::vector<bool> f = bv(d, "usewithcau");
  for (std::size_t i = 0; i < l.size() && i < f.size(); ++i) {
    if (f[i]) l[i] += " *";
  }
  return l;
}
inline bool rindgap(In s) {
  CK(need(s, "indgap", {"group", "ratio"}));
  Dots d = tab(s.d, "group", "ratio", "ratiolo", "ratioup", "Ratio to the reference group",
               "Ratio (log scale); 1 is no gap", 1);
  d.lab = star(s.d);
  const std::string rf = str(get(s.a, "ref"));
  if (!rf.empty()) d.main = "Ratio to " + rf;
  const std::vector<bool> po = bv(s.d, "poogrp");
  if (std::find(po.begin(), po.end(), true) != po.end()) {
    for (bool b : po) d.grp.push_back(b ? "pooled total" : "distinct group");
  }
  d.logx = true;
  push(s, dots(s.f, d));
  return true;
}
// Two or more columns of one table, each a group of the plot.
inline Dots cols2(const Val &t, const ax::Cols &lab, std::initializer_list<const char *> est,
                  std::initializer_list<const char *> lo_, std::initializer_list<const char *> hi_,
                  const ax::Cols &names, const char *main, const char *xlab) {
  Dots d;
  d.main = main;
  d.xlab = xlab;
  std::size_t k = 0;
  for (const char *e : est) {
    const std::vector<double> v = nv(t, e);
    const std::vector<double> l = lo_.size() ? nv(t, *(lo_.begin() + k)) : std::vector<double>(v.size(), NA);
    const std::vector<double> h = hi_.size() ? nv(t, *(hi_.begin() + k)) : std::vector<double>(v.size(), NA);
    for (std::size_t i = 0; i < v.size(); ++i) {
      d.lab.push_back(lab[i]);
      d.est.push_back(v[i]);
      d.lo.push_back(l[i]);
      d.hi.push_back(h[i]);
      d.grp.push_back(names[k]);
    }
    ++k;
  }
  return d;
}
inline bool rdualcont(In s) {
  CK(need(s, "dualcont", {"group", "wellness", "distress"}));
  push(s, dots(s.f, cols2(s.d, star(s.d), {"wellness", "distress"}, {"wello", "distrlo"}, {"welup", "distrup"},
                          {"wellness", "distress"}, "Wellness beside distress", "Estimate")));
  return true;
}
inline bool rcommidx(In s) {
  CK(need(s, "commidx", {"group", "index"}));
  Dots d = tab(s.d, "group", "index", "lower", "upper", "Community-weighted index", "Index, under the weights given");
  d.lab = star(s.d);
  push(s, dots(s.f, d));
  return true;
}
inline bool rstdprev(In s) {
  CK(need(s, "stdprev", {"group", "estimate"}));
  Dots d = tab(s.d, "group", "estimate", "lower", "upper", "Age-standardised prevalence", "Prevalence");
  d.lab = star(s.d);
  push(s, dots(s.f, d));
  return true;
}
inline bool rordcmp(In s) {
  CK(need(s, "ordcmp", {"group", "probbet"}));
  const ax::Cols rf = sv_(s.d, "reference");
  Dots d = tab(s.d, "group", "probbet", "lower", "upper", "", "Probability; 0.5 is no difference", 0.5);
  d.main = "Chance of a higher rating than " + (rf.empty() ? std::string("the reference") : rf[0]);
  d.lab = star(s.d);
  push(s, dots(s.f, d));
  return true;
}
// The rows of the profile and change tables: indicator (side), by group.
inline Dots profile(const Val &t, const char *est, const char *main, const char *xlab, double r) {
  const ax::Cols ind = sv_(t, "indicator"), side = sv_(t, "side"), g = star(t);
  const std::vector<double> e = nv(t, est), l = nv(t, "lower"), h = nv(t, "upper");
  Dots d;
  d.main = main;
  d.xlab = xlab;
  d.ref = r;
  for (std::size_t i = 0; i < e.size(); ++i) {
    if (!std::isfinite(e[i])) continue;
    d.lab.push_back(ind[i] + " (" + side[i] + ")");
    d.est.push_back(e[i]);
    d.lo.push_back(l[i]);
    d.hi.push_back(h[i]);
    d.grp.push_back(g[i]);
  }
  return d;
}
inline bool rwelprof(In s) {
  CK(need(s, "welprof", {"group", "indicator", "estimate"}));
  push(s, dots(s.f, profile(s.d, "estimate", "Wellness profile", "Estimate", NA)));
  return true;
}
inline bool rhealchg(In s) {
  CK(need(s, "healchg", {"group", "indicator", "change"}));
  push(s, dots(s.f, profile(s.d, "change", "Change between the two periods",
                            "Change in the estimate's units; 0 is none", 0)));
  return true;
}
inline bool rholbal(In s) {
  CK(need(s, "holbal", {"group", "level", "balance"}));
  Dots a = tab(s.d, "group", "level", "lvllo", "lvlup", "Level across the domains", "Mean level");
  Dots b = tab(s.d, "group", "balance", "ballo", "balup", "Balance of the domains", "Balance; 1 is even", 1);
  a.lab = b.lab = star(s.d);
  push(s, dots(s.f, a));
  push(s, dots(s.f, b));
  return true;
}
inline bool rtrchlthgap(In s) {
  CK(need(s, "trchlthgap", {"group", "disratio"}));
  Dots d = tab(s.d, "group", "disratio", "lower", "upper", "Health gaps (Call to Action 19)",
               "Disparity ratio (log scale); above 1 is worse", 1);
  d.lab = star(s.d);
  const ax::Cols ind = sv_(s.d, "indicator"), pe = sv_(s.d, "period");
  for (std::size_t i = 0; i < ind.size(); ++i) d.grp.push_back(ind[i] + " " + pe[i]);
  d.logx = true;
  push(s, dots(s.f, d));
  return true;
}
// .pcwb(): the four components of a Community Well-Being score and the score.
inline bool cwb(In s, const ax::Cols &lab, const char *main) {
  push(s, dots(s.f, cols2(s.d, lab, {"income", "education", "housing", "labour", "cwb"}, {}, {},
                          {"income", "education", "housing", "labour force", "score"}, main,
                          "0 to 100; the score is the mean of the four")));
  return true;
}
inline bool rcnscwb(In s) {
  CK(need(s, "cnscwb", {"geo", "group", "income"}));
  const ax::Cols geo = sv_(s.d, "geo"), g = sv_(s.d, "group");
  ax::Cols lab;
  for (std::size_t i = 0; i < geo.size(); ++i) lab.push_back(geo[i] + ": " + g[i]);
  return cwb(s, lab, "Community Well-Being (census)");
}
inline bool rcwbscr(In s) {
  CK(need(s, "cwbscr", {"income", "cwb"}));
  return cwb(s, seqlab("community", nv(s.d, "cwb").size()), "Community Well-Being scores");
}

// ---- the exposure models ------------------------------------------------

// .axci(): est - z se and est + z se at the level of the fit.
inline void levci(Dots &d, const std::vector<double> &se, double level) { zci(d, se, clab_qnorm(0.5 + level / 2)); }
// .axcoef() after .axterm(): the coefficients of a fit, those named by
// `terms` when it is given; as rate ratios on a log axis with `rr`.
inline bool coefs(In s, bool rr, const char *main) {
  const Val &t = get(s.d, "coef");
  const ax::Cols all = sv_(t, "term"), want = sv_(*s.c.in, "terms");
  for (const std::string &w : want) {
    if (!ax::has(all, w)) {
      return s.c.fail("`terms` must be NULL or names of terms of the fit: %s", ax::join(all, ", ").c_str());
    }
  }
  const Dots a = tab(t, "term", "est", "lo", "hi", main, rr ? "Rate ratio (log scale)" : "Coefficient", rr ? 1 : 0);
  Dots d = a;
  d.lab.clear();
  d.est.clear();
  d.lo.clear();
  d.hi.clear();
  for (std::size_t i = 0; i < a.lab.size(); ++i) {
    if (!want.empty() && !ax::has(want, a.lab[i])) continue;
    d.lab.push_back(a.lab[i]);
    d.est.push_back(rr ? std::exp(a.est[i]) : a.est[i]);
    d.lo.push_back(rr ? std::exp(a.lo[i]) : a.lo[i]);
    d.hi.push_back(rr ? std::exp(a.hi[i]) : a.hi[i]);
  }
  d.logx = rr;
  push(s, dots(s.f, d));
  return true;
}
// .axfst(): the strength of the first stage, four statistics per exposure.
inline void first(In s) {
  push(s, dots(s.f, cols2(get(s.d, "first"), sv_(get(s.d, "first"), "endo"), {"wald", "kp", "eff", "cv"}, {}, {},
                          {"Wald F", "Kleibergen-Paap F", "effective F", "5% critical value"},
                          "First-stage strength", "F statistic")));
}
inline bool rppmlfe(In s) {
  CK(need(s, "ppmlfe", {"coef"}));
  return coefs(s, true, "Poisson fixed-effects estimates");
}
inline bool rexpgmm(In s) {
  CK(need(s, "expgmm", {"coef"}));
  return coefs(s, true, "Exponential-mean GMM estimates");
}
inline bool rivfe(In s) {
  CK(need(s, "ivfe", {"coef", "first"}));
  const bool iv = !sv_(get(s.d, "first"), "endo").empty() && has(s.d, "fstab");
  CK(coefs(s, false, iv ? "Two-stage least-squares estimates" : "Least-squares estimates"));
  if (!iv) return true;
  const Val &fb = get(s.d, "fstab");
  Dots d = tab(fb, "inst", "est", nullptr, nullptr, "First-stage coefficients", "Change in the exposure", 0);
  const ax::Cols en = sv_(fb, "endo");
  for (std::size_t i = 0; i < d.lab.size(); ++i) d.lab[i] += " (" + en[i] + ")";
  levci(d, nv(fb, "se"), num1(get(s.d, "level")));
  push(s, dots(s.f, d));
  first(s);
  return true;
}
inline bool rppmlcf(In s) {
  CK(need(s, "ppmlcf", {"coef", "first"}));
  CK(coefs(s, true, str(get(s.d, "how")) == "plug" ? "Poisson plug-in estimates"
                                                    : "Poisson control-function estimates"));
  first(s);
  return true;
}
// .axprof(): estimates along a quantity, joined by a dotted line, with
// their intervals; the marked ones as squares in the second colour.
struct Prof {
  std::vector<double> x, est, lo, hi, mark;
  std::string main, xlab, ylab;
  ax::Cols key, lab;
};
inline Pan prof(const Fig &f, const Prof &d) {
  const ax::Cols cl = cols(f, 2);
  Pan p = frame(d.main.c_str(), d.xlab.c_str(), d.ylab.c_str());
  rng(d.x, p.x0, p.x1);
  std::vector<double> all = d.est;
  all.insert(all.end(), d.lo.begin(), d.lo.end());
  all.insert(all.end(), d.hi.begin(), d.hi.end());
  all.push_back(0);
  rng(all, p.y0, p.y1);
  if (!d.lab.empty()) labs(p.ax, d.lab, d.x, true);
  const double wx = p.x1 - p.x0;
  p.it.push_back(ref(p.x0 - wx, 0, p.x1 + wx, 0, 2));
  std::vector<std::size_t> o(d.x.size());
  for (std::size_t i = 0; i < o.size(); ++i) o[i] = i;
  std::stable_sort(o.begin(), o.end(), [&](std::size_t a, std::size_t b) { return d.x[a] < d.x[b]; });
  std::vector<double> lx, ly;
  for (std::size_t i : o) {
    lx.push_back(d.x[i]);
    ly.push_back(d.est[i]);
  }
  p.it.push_back(line(lx, ly, cl[0], 3, 1.2));
  bool any = false;
  for (int m = 0; m < 2; ++m) {
    It sg, pt;
    sg.kind = SEG;
    sg.lw = 1.6;
    pt.kind = PTS;
    pt.pch = m ? 3 : 0;
    pt.r = 3.8;
    sg.col = pt.col = cl[sat<std::size_t>(m)];
    for (std::size_t i = 0; i < d.x.size(); ++i) {
      const bool mk = std::find(d.mark.begin(), d.mark.end(), d.x[i]) != d.mark.end();
      if (mk != (m == 1)) continue;
      any = any || mk;
      pt.x.push_back(d.x[i]);
      pt.y.push_back(d.est[i]);
      if (i < d.lo.size() && std::isfinite(d.lo[i]) && std::isfinite(d.hi[i])) {
        sg.x.insert(sg.x.end(), {d.x[i], d.x[i]});
        sg.y.insert(sg.y.end(), {d.lo[i], d.hi[i]});
      }
    }
    p.it.push_back(sg);
    p.it.push_back(pt);
  }
  if (any && d.key.size() == 2) {
    p.key = {{d.key[0], cl[0], 0, 0}, {d.key[1], cl[1], 0, 3}};
    p.corner = 1;
  }
  return p;
}
inline Prof proftab(const Val &t, const char *x) {
  Prof d;
  d.x = nv(t, x);
  d.est = nv(t, "est");
  d.lo = nv(t, "lo");
  d.hi = nv(t, "hi");
  return d;
}
inline bool rupdwn(In s) {
  CK(need(s, "updwn", {"coef", "class"}));
  const Val &cf = get(s.d, "coef");
  const ax::Cols tm = sv_(cf, "term");
  if (tm.size() > 1) {
    Prof d = proftab(cf, "est");
    d.x.clear();
    // .axring(): the midpoint of the ring a term names, "name:[lo,hi)"
    for (const std::string &t : tm) {
      const std::size_t a = t.find_first_of("[("), c = t.find(',', a == std::string::npos ? 0 : a);
      d.x.push_back(a == std::string::npos || c == std::string::npos
                        ? NA
                        : (std::strtod(t.c_str() + a + 1, nullptr) + std::strtod(t.c_str() + c + 1, nullptr)) / 2);
    }
    d.main = "Downwind against upwind";
    d.xlab = "Distance from the source (ring midpoint)";
    d.ylab = "Downwind difference";
    push(s, prof(s.f, d));
  } else {
    push(s, dots(s.f, tab(cf, "term", "est", "lo", "hi", "Downwind against upwind", "Downwind difference", 0)));
  }
  if (has(s.d, "placebo")) {
    Prof d = proftab(get(s.d, "placebo"), "off");
    d.main = "Placebo segments moved sideways";
    d.xlab = "Offset";
    d.ylab = "Downwind difference";
    d.mark = {0};
    d.key = {"moved", "the source itself"};
    push(s, prof(s.f, d));
  }
  const Val &k = get(s.d, "class");
  push(s, bars(s.f, ax::nums(k), k.k, "Rows by classification", "Rows"));
  return true;
}
inline bool rrunfx(In s) {
  CK(need(s, "runfx", {"coef"}));
  const Val &cf = get(s.d, "coef");
  Prof d = proftab(cf, "est");
  d.x.clear();
  for (const std::string &t : sv_(cf, "term")) {
    // .axrun(): the day a term names, "run3" or "run5p" (5 and beyond)
    const double k = std::strtod(t.c_str() + (t.compare(0, 3, "run") == 0 ? 3 : 0), nullptr);
    d.x.push_back(k);
    d.lab.push_back(gfmt(k, 7) + (!t.empty() && t.back() == 'p' ? "+" : ""));
  }
  d.main = "Effect by day of a run of exposed days";
  d.xlab = "Day of the run";
  d.ylab = "Coefficient";
  push(s, prof(s.f, d));
  return true;
}
inline bool rdosebin(In s) {
  CK(need(s, "dosebin", {"tab", "test"}));
  const Val &t = get(s.d, "tab"), &tt = get(s.d, "test");
  Prof d = proftab(t, "mid");
  const std::vector<double> se = nv(t, "se"), pv = nv(tt, "p");
  for (std::size_t i = 0; i < d.x.size(); ++i) {
    if (!std::isfinite(se[i]) && d.est[i] == 0) d.mark.push_back(d.x[i]);
  }
  const ax::Cols nm = sv_(tt, "test");
  std::string tests;
  for (std::size_t i = 0; i < nm.size(); ++i) tests += (i ? "; " : "") + nm[i] + " p " + gfmt(pv[i], 2);
  d.main = "Effect of each bin against the reference";
  d.xlab = "Exposure, median of each bin (" + tests + ")";
  d.ylab = "Coefficient";
  d.key = {"bin", "reference bin"};
  push(s, prof(s.f, d));
  return true;
}
inline bool rjntx(In s) {
  CK(need(s, "jntx", {"cond", "surf"}));
  const Val &cd = get(s.d, "cond"), &sf = get(s.d, "surf");
  Dots d = tab(cd, "x2", "est", "lo", "hi", "Effect of x1 at levels of x2", "Rate ratio of x1 (log scale)", 1);
  const std::vector<double> p2 = nv(cd, "p2"), x2 = nv(cd, "x2");
  for (std::size_t i = 0; i < d.lab.size(); ++i) {
    d.lab[i] = "x2 at " + gfmt(100 * p2[i], 7) + "% (" + gfmt(x2[i], 3) + ")";
    d.est[i] = std::exp(d.est[i]);
    d.lo[i] = std::exp(d.lo[i]);
    d.hi[i] = std::exp(d.hi[i]);
  }
  d.logx = true;
  push(s, dots(s.f, d));
  // .axsurf(): the log rate ratio over the grid of the two exposures
  const std::vector<double> a = nv(sf, "x1"), b = nv(sf, "x2"), rr = nv(sf, "rr");
  std::vector<double> u = a, v = b;
  std::sort(u.begin(), u.end());
  u.erase(std::unique(u.begin(), u.end()), u.end());
  std::sort(v.begin(), v.end());
  v.erase(std::unique(v.begin(), v.end()), v.end());
  std::vector<std::vector<double>> z(u.size(), std::vector<double>(v.size(), NA));
  for (std::size_t i = 0; i < rr.size(); ++i) {
    const std::size_t iu = sat<std::size_t>(std::lower_bound(u.begin(), u.end(), a[i]) - u.begin());
    const std::size_t iv = sat<std::size_t>(std::lower_bound(v.begin(), v.end(), b[i]) - v.begin());
    z[iu][iv] = std::log(rr[i]);
  }
  push(s, surf(s.f, u, v, z, "Rate ratio against the reference point", "First exposure (x1)",
               "Second exposure (x2)", true));
  return true;
}

// ---- points on a plane, and the forensic curves --------------------------

// ptmap(): the points in the colours and shapes the function gives them,
// the centre of each cluster marked with its number.
inline bool rptmap(In s) {
  CK(need(s, "ptmap", {"pts"}));
  const Val &t = get(s.d, "pts");
  const std::vector<double> x = nv(t, "x"), y = nv(t, "y");
  const ax::Cols col = sv_(t, "col");
  Pan p = frame("Points", "x", "y");
  rng(x, p.x0, p.x1);
  rng(y, p.y0, p.y1);
  It pt;
  pt.kind = PTS;
  pt.x = x;
  pt.y = y;
  pt.cols = col;
  p.it.push_back(pt);
  if (has(s.d, "cent")) {
    const Val &c = get(s.d, "cent");
    const std::vector<double> cx = nv(c, "x"), cy = nv(c, "y"), cl = nv(c, "cl");
    for (std::size_t i = 0; i < cx.size(); ++i) {
      if (std::isfinite(cx[i]) && std::isfinite(cy[i])) p.it.push_back(text(cx[i], cy[i], gfmt(cl[i], 7), 1));
    }
  }
  push(s, p);
  return true;
}
// plot() of dbscn(): the points it was given, by cluster (0 is noise).
inline bool rdbscn(In s) {
  CK(need(s, "dbscn", {"cl"}));
  const Val &xo = get(s.a, "x");
  std::vector<double> x, y;
  if (xo.t == Val::OBJ) {
    x = nv(xo, "x");
    y = nv(xo, "y");
  } else if (xo.t == Val::ARR && !xo.a.empty() && xo.a[0].t == Val::ARR) {
    // a matrix of two columns, written as rows
    for (const std::vector<double> &r : mat(xo)) {
      x.push_back(r.size() > 0 ? r[0] : NA);
      y.push_back(r.size() > 1 ? r[1] : NA);
    }
  } else {
    x = ax::nums(xo);
    y = nv(s.a, "y");
  }
  std::vector<double> cl = nv(s.d, "cl");
  for (double &v : cl) v = std::isfinite(v) ? v : 0;  // noise has no cluster
  if (x.size() != cl.size() || y.size() != cl.size()) {
    return s.c.fail("the points of `x` and `y` could not be matched to the clusters.");
  }
  std::vector<double> lv = cl;
  std::sort(lv.begin(), lv.end());
  lv.erase(std::unique(lv.begin(), lv.end()), lv.end());
  const ax::Cols c0 = cols(s.f, sat<int>(lv.size()));
  Pan p = frame("DBSCAN clusters", "x", "y");
  rng(x, p.x0, p.x1);
  rng(y, p.y0, p.y1);
  for (std::size_t k = 0; k < lv.size(); ++k) {
    std::vector<double> qx, qy;
    for (std::size_t i = 0; i < cl.size(); ++i) {
      if (cl[i] == lv[k]) {
        qx.push_back(x[i]);
        qy.push_back(y[i]);
      }
    }
    const bool noise = lv[k] == 0;
    p.it.push_back(points(qx, qy, noise ? "#969696" : c0[k], noise ? 2 : pch(k)));
    p.key.push_back({noise ? "noise" : "cluster " + gfmt(lv[k], 7), noise ? "#969696" : c0[k], 0, noise ? 2 : pch(k)});
  }
  p.corner = 1;
  push(s, p);
  return true;
}
inline bool rtipplt(In s) {
  CK(need(s, "tipplt", {"same", "diff"}));
  const ax::Cols cl = two(s.f);
  Pan p = frame("Tippett plot", "log10 likelihood ratio", "share of ratios at or above");
  std::vector<double> all;
  int k = 0;
  for (const char *key : {"same", "diff"}) {
    const Val &t = get(s.d, key);
    It l = line(nv(t, "log10lr"), nv(t, "share"), cl[sat<std::size_t>(k)], k + 1, 2);
    l.kind = STEP;
    all.insert(all.end(), l.x.begin(), l.x.end());
    p.it.push_back(l);
    p.key.push_back({k ? "different source" : "same source", cl[sat<std::size_t>(k)], k + 1, -1});
    ++k;
  }
  all.push_back(0);
  rng(all, p.x0, p.x1);
  p.it.push_back(ref(0, -1, 0, 2, 3));
  p.corner = 2;
  push(s, p);
  return true;
}
inline bool receplt(In s) {
  CK(need(s, "eceplt", {"lo", "ece"}));
  const std::vector<double> x = nv(s.d, "lo");
  Grp d = series(get(s.d, "lo"), {"reported", "calibrated", "always 1"},
                 {nv(s.d, "ece"), nv(s.d, "min"), nv(s.d, "ref")});
  titled(d, "Empirical cross-entropy", "prior log10 odds", "empirical cross-entropy (bits)");
  d.type = 'l';
  if (s.f.palette == "grey") d.col = {"#000000", "#000000", "#000000"};
  d.corner = 0;
  push(s, grp(s.f, d));
  return true;
}

// ---- regression discontinuity (R/rdplt.R) --------------------------------

// .rdwls(): the weighted least-squares line y = a + b u; false when the
// weights rest on fewer than two distinct values of u.
inline bool wls(const std::vector<double> &u, const std::vector<double> &y, const std::vector<double> &w, double &a,
                double &b) {
  double s0 = 0, s1 = 0, s2 = 0, t0 = 0, t1 = 0;
  for (std::size_t i = 0; i < u.size(); ++i) {
    s0 += w[i];
    s1 += w[i] * u[i];
    s2 += w[i] * u[i] * u[i];
    t0 += w[i] * y[i];
    t1 += w[i] * u[i] * y[i];
  }
  const double dn = s0 * s2 - s1 * s1;
  if (!(dn > 1e-12 * s0 * s2)) return false;
  a = (s2 * t0 - s1 * t1) / dn;
  b = (s0 * t1 - s1 * t0) / dn;
  return true;
}
// .rddraw(): the bins as points (rings below the cutoff, discs from it on)
// and the fitted line of each side.
struct Rd {
  std::vector<double> bx[2], by[2], fx[2], fy[2];
};
inline void rddraw(In s, const Rd &d, double cut, const std::string &main, const std::string &xlab,
                   const char *ylab) {
  const ax::Cols cl = cols(s.f, 2);
  Pan p = frame(main.c_str(), xlab.c_str(), ylab);
  std::vector<double> ax_, ay_;
  for (int k = 0; k < 2; ++k) {
    ax_.insert(ax_.end(), d.bx[k].begin(), d.bx[k].end());
    ax_.insert(ax_.end(), d.fx[k].begin(), d.fx[k].end());
    ay_.insert(ay_.end(), d.by[k].begin(), d.by[k].end());
    ay_.insert(ay_.end(), d.fy[k].begin(), d.fy[k].end());
  }
  rng(ax_, p.x0, p.x1);
  rng(ay_, p.y0, p.y1);
  const double wy = p.y1 - p.y0;
  p.it.push_back(ref(cut, p.y0 - wy, cut, p.y1 + wy, 2));
  for (int k = 0; k < 2; ++k) {
    const std::string &c0 = cl[sat<std::size_t>(k)];
    p.it.push_back(points(d.bx[k], d.by[k], c0, k ? 0 : 2));
    p.it.push_back(line(d.fx[k], d.fy[k], c0, k + 1, 2));
    p.key.push_back({k ? "from the cutoff on" : "below the cutoff", c0, k + 1, k ? 0 : 2});
  }
  p.corner = 0;
  push(s, p);
}
// plot() of rdloc(): the outcome's means in 20 bins on each side and the
// local-linear line of each side over the rows the kernel weighs.
inline bool rrdloc(In s) {
  CK(need(s, "rdloc", {"coef", "h", "cut"}));
  Val dt;
  CK(ax::dat(s.c, get(s.a, "data"), "data", dt));
  const std::string yn = has(s.a, "y") ? str(get(s.a, "y")) : "y", rn = has(s.a, "run") ? str(get(s.a, "run")) : "run";
  const std::vector<double> y0 = nv(dt, yn.c_str()), r0 = nv(dt, rn.c_str());
  const double cut = num1(get(s.d, "cut")), h = num1(get(s.d, "h"));
  const std::string kern = str(get(s.d, "kern"));
  std::vector<double> u[2], y[2];
  for (std::size_t i = 0; i < y0.size() && i < r0.size(); ++i) {
    if (!std::isfinite(y0[i]) || !std::isfinite(r0[i])) continue;
    const int k = r0[i] - cut < 0 ? 0 : 1;
    u[k].push_back(r0[i] - cut);
    y[k].push_back(y0[i]);
  }
  if (u[0].empty() || u[1].empty()) return s.c.fail("the rows of `data` could not be matched to the fit.");
  Rd d;
  const int nb = 20;
  for (int k = 0; k < 2; ++k) {
    // .rdbin(): equal bins over the side, the empty ones left out
    const double e0 = k ? 0 : lo(u[0]), e1 = k ? hi(u[1]) : 0, wd = (e1 - e0) / nb;
    std::vector<double> sum(nb, 0.0), cnt(nb, 0.0);
    for (std::size_t i = 0; i < u[k].size(); ++i) {
      const int g = wd > 0 ? std::min(nb - 1, std::max(0, sat<int>(std::floor((u[k][i] - e0) / wd)))) : 0;
      sum[sat<std::size_t>(g)] += y[k][i];
      cnt[sat<std::size_t>(g)] += 1;
    }
    for (int g = 0; g < nb; ++g) {
      if (cnt[sat<std::size_t>(g)] == 0) continue;
      d.bx[k].push_back(cut + (wd > 0 ? e0 + (g + 0.5) * wd : e0));
      d.by[k].push_back(sum[sat<std::size_t>(g)] / cnt[sat<std::size_t>(g)]);
    }
    // .rdkw(): the kernel's weights at u / h
    std::vector<double> ku, ky, kw;
    for (std::size_t i = 0; i < u[k].size(); ++i) {
      const double v = u[k][i] / h;
      const double w = kern == "uni" ? (std::fabs(v) <= 1 ? 0.5 : 0)
                                     : kern == "epa" ? std::max(0.75 * (1 - v * v), 0.0) : std::max(1 - std::fabs(v), 0.0);
      if (w > 0) {
        ku.push_back(u[k][i]);
        ky.push_back(y[k][i]);
        kw.push_back(w);
      }
    }
    double a = NA, b = NA;
    if (!ku.empty() && wls(ku, ky, kw, a, b)) {
      const double x0 = k ? 0 : lo(ku), x1 = k ? hi(ku) : 0;
      d.fx[k] = {cut + x0, cut + x1};
      d.fy[k] = {a + b * x0, a + b * x1};
    }
  }
  const Val &cf = get(s.d, "coef");
  rddraw(s, d, cut, "Outcome around the cutoff", rn + " (estimate " + gfmt(nv(cf, "est")[0], 3) + ", p " +
                                                    gfmt(nv(cf, "p")[0], 2) + ")", yn.c_str());
  return true;
}
// plot() of rddens(): McCrary's first-step histogram and the local-linear
// smooth of its heights on each side, at 40 points up to the cutoff.
inline bool rrddens(In s) {
  CK(need(s, "rddens", {"tab", "b", "h", "cut"}));
  const std::vector<double> r = nv(s.a, "r");
  const double cut = num1(get(s.d, "cut")), b = num1(get(s.d, "b")), h = num1(get(s.d, "h"));
  if (r.empty() || !(b > 0) || !(h > 0)) return s.c.fail("the values of `r` could not be matched to the fit.");
  double jlo = HUGE_VAL, jhi = -HUGE_VAL;
  for (double v : r) {
    jlo = std::min(jlo, std::floor((v - cut) / b));
    jhi = std::max(jhi, std::floor((v - cut) / b));
  }
  const std::size_t m = sat<std::size_t>(jhi - jlo + 1);
  if (m > 100000) return s.c.fail("the bin width `b` gives more than 100000 bins to draw.");
  std::vector<double> ht(m, 0.0);
  for (double v : r) ht[sat<std::size_t>(std::floor((v - cut) / b) - jlo)] += 1 / (static_cast<double>(r.size()) * b);
  Rd d;
  for (std::size_t i = 0; i < m; ++i) {
    const double x = (jlo + static_cast<double>(i)) * b + b / 2 + cut;
    const int k = x < cut ? 0 : 1;
    d.bx[k].push_back(x);
    d.by[k].push_back(ht[i]);
  }
  const int ng = 40;
  for (int k = 0; k < 2; ++k) {
    if (d.bx[k].empty()) continue;
    const double g0 = k ? 0 : lo(d.bx[0]) - cut, g1 = k ? hi(d.bx[1]) - cut : 0;
    for (int g = 0; g < ng; ++g) {
      const double at = g0 + (g1 - g0) * g / (ng - 1);
      std::vector<double> u, w;
      for (double x : d.bx[k]) {
        u.push_back(x - cut - at);
        w.push_back(std::max(1 - std::fabs((x - cut - at) / h), 0.0));
      }
      double a = NA, sl = NA;
      if (wls(u, d.by[k], w, a, sl) && std::isfinite(a)) {
        d.fx[k].push_back(at + cut);
        d.fy[k].push_back(a);
      }
    }
  }
  const Val &t = get(s.d, "tab");
  rddraw(s, d, cut, "Density at the cutoff", "Running variable (log difference " + gfmt(nv(t, "theta")[0], 3) +
                                                 ", p " + gfmt(nv(t, "p")[0], 2) + ")",
         "Density: bin heights and fits");
  return true;
}
