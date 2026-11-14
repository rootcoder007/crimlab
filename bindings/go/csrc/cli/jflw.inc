// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Stocks and flows of cases and people (src/bflw.cpp, src/bpub.cpp, and
// the band readers of R/bandx.R): stays and their bands, court dockets,
// stock and flow, backlogs, life tables, transition matrices, published
// bands and suppressed or rounded cells.

// A refusal whose words name another function (an R helper's stop()).
bool failas(Ctx &c, const std::string &m) {
  c.err = m;
  return false;
}

// plab(): period labels as given, or 1, 2, ...
bool plab(Ctx &c, std::size_t n, std::vector<std::string> &s) {
  if (!c.has("prd")) {
    s = nmsof({}, n);
    return true;
  }
  const Val &p = c["prd"];
  if (p.t == Val::ARR) {
    for (const Val &e : p.a) {
      if (e.t == Val::ARR || e.t == Val::OBJ) return notvec(c, "prd");
    }
  }
  if (vlen(p) != n) return c.fail("`prd` must have one label per period, %d.", static_cast<int>(n));
  if (!alabs(p, s)) return c.fail("`prd` must be labels: numbers, text or a factor.");
  return true;
}

// rci(): the exact interval of the ratio a / b from the conditional
// binomial of a among a + b.
void rci(double a, double b, double p, double &lo, double &hi) {
  if (std::isnan(a) || std::isnan(b) || a + b == 0) {
    lo = hi = NA;
    return;
  }
  const double l = a == 0 ? 0 : clab_qbeta(p, a, b + 1, 1);
  const double h = b == 0 ? 1 : clab_qbeta(1 - p, a + 1, b, 1);
  lo = l / (1 - l);
  hi = h >= 1 ? Inf : h / (1 - h);
}

// One number, or one per period, recycled to n.
bool rec1(Ctx &c, const Val &v, const char *nm, std::size_t n, std::vector<double> &out) {
  Dv d;
  CK(num(c, v, nm, d, -1, 0.0, true));
  if (d.n() != 1 && d.n() != n) return c.fail("`%s` must be one number or one per period, %d.", nm, static_cast<int>(n));
  out.assign(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) out[i] = d[d.n() == 1 ? 0 : i];
  return true;
}

bool lostay(Ctx &c, Val &out) {
  Dv xv;
  CK(num(c, c["days"], "days", xv, -1, 0.0));
  const std::size_t sn = xv.n();
  const int n = static_cast<int>(sn);
  Dv wv;
  if (c.has("w")) {
    CK(num(c, c["w"], "w", wv, n, 0.0));
    double tot = 0;
    for (double e : wv.v) {
      if (e != std::floor(e)) return c.fail("`w` must be whole-number counts.");
      tot += e;
    }
    if (!(tot > 0)) return c.fail("`w` must have a positive total.");
  }
  std::vector<double> bv;
  int nb = 0;
  if (c.has("breaks")) {
    const Val &b = c["breaks"];
    if (!isnum(b) || vlen(b) < 2) return c.fail("`breaks` must be at least two band edges.");
    bv = dbl(b);
    nb = static_cast<int>(bv.size());
    for (int k = 0; k < nb; ++k) {
      const std::size_t sk = static_cast<std::size_t>(k);
      if (std::isnan(bv[sk]) || (k > 0 && !(bv[sk] > bv[sk - 1]))) {
        return c.fail("`breaks` must be increasing band edges with no missing values.");
      }
    }
  }
  Dv tv;
  int nt = 0;
  if (c.has("thr")) {
    CK(num(c, c["thr"], "thr", tv, -1, 0.0));
    nt = static_cast<int>(tv.n());
  }
  std::vector<int> fv;
  if (c.has("flag")) CK(lgl(c, c["flag"], "flag", n, fv));
  double lv = 0.0;
  CK(oned(c, "level", 0.95, 0.0, 1.0, lv));
  double st[9] = {0}, fit[6] = {0};
  std::vector<double> bw(2 * static_cast<std::size_t>(nb) + static_cast<std::size_t>(nt) + 2);
  std::vector<int> ow(sn);
  double *bt = bw.data(), *tt = bt + 2 * nb;
  clab_lostay(xv.p(), wv.n() ? wv.p() : nullptr, n, bv.empty() ? nullptr : bv.data(), nb, nt ? tv.p() : nullptr, nt,
              fv.empty() ? nullptr : fv.data(), st, bt, tt, fit, ow.data());
  const double nn = st[0];
  const double q = nn > 1 ? clab_qt(0.5 + lv / 2, nn - 1, 1) : NA;
  const double sv[11] = {nn,    st[1], st[2], na(st[3]), st[4], st[5], st[6], st[7], na(st[8]),
                         nn > 1 ? st[2] - q * st[8] : NA, nn > 1 ? st[2] + q * st[8] : NA};
  static const char *const sn2[] = {"n", "days", "mean", "sd", "median", "q1", "q3", "max", "se", "lo", "hi"};
  out = vobj();
  put(out, "stay", nvec(11, sn2, sv));
  Val bands, thr;
  if (nb >= 2) {
    const std::size_t k = static_cast<std::size_t>(nb - 1);
    std::vector<double> gl(k), gh(k), gn(k), gs(k), gf(k);
    for (std::size_t j = 0; j < k; ++j) {
      gl[j] = bv[j];
      gh[j] = bv[j + 1];
      gn[j] = bt[j];
      gs[j] = bt[j] / nn;
      gf[j] = bt[j] > 0 ? bt[k + j] / bt[j] : NA;
    }
    Frame g;
    g.num("lo", gl);
    g.num("hi", gh);
    g.num("n", gn);
    g.num("share", gs);
    if (!fv.empty()) g.num("flagged", gf);
    bands = g.o;
  }
  if (nt > 0) {
    std::vector<double> at(tt, tt + nt), ab(static_cast<std::size_t>(nt));
    for (std::size_t i = 0; i < ab.size(); ++i) ab[i] = 1 - at[i];
    Frame h;
    h.num("thr", tv.v);
    h.num("atmost", at);
    h.num("above", ab);
    thr = h.o;
  }
  put(out, "bands", bands);
  put(out, "thr", thr);
  static const char *const mo[] = {"lognormal", "exponential", "pareto"};
  std::vector<double> ma(3);
  std::vector<int> mb(3);
  int best = -1;
  for (int k = 0; k < 3; ++k) {
    ma[static_cast<std::size_t>(k)] = na(fit[k]);
    if (!std::isnan(fit[k]) && (best < 0 || fit[k] < fit[best])) best = k;
  }
  for (int k = 0; k < 3; ++k) mb[static_cast<std::size_t>(k)] = k == best;
  Frame m;
  m.str("model", {mo[0], mo[1], mo[2]});
  m.num("aic", ma);
  m.lgl("best", mb.data(), 3);
  put(out, "fit", m.o);
  static const char *const pn[] = {"meanlog", "sdlog", "alpha"};
  const double pv[3] = {na(fit[3]), na(fit[4]), na(fit[5])};
  put(out, "par", nvec(3, pn, pv));
  return true;
}

bool bklog(Ctx &c, Val &out) {
  Dv av;
  CK(num(c, c["arrive"], "arrive", av, -1, 0.0));
  const std::size_t sn = av.n();
  const int n = static_cast<int>(sn);
  const Val &dp = c["dispose"];
  if (dp.t == Val::ARR) {
    for (const Val &e : dp.a) {
      if (e.t == Val::ARR || e.t == Val::OBJ) return notvec(c, "dispose");
    }
  }
  if (!isnum(dp) || vlen(dp) != sn) {
    return c.fail("`dispose` must be numbers (NA for a pending case), one per case, %d.", n);
  }
  std::vector<double> dv = dbl(dp);
  double h = 0;
  for (std::size_t i = 0; i < sn; ++i) {
    const double d = dv[i];
    if (std::isinf(d)) return c.fail("`dispose` must be finite or NA.");
    if (!std::isnan(d) && d < av[i]) return c.fail("a disposition precedes its arrival (case %d).", static_cast<int>(i) + 1);
    h = std::fmax(h, std::fmax(av[i], std::isnan(d) ? 0 : d));
  }
  if (c.has("horizon")) CK(one(c, c["horizon"], "horizon", 0.0, Inf, h));
  double o[7] = {0};
  if (clab_bklog(av.p(), dv.data(), n, h, o) != CLAB_OK) return c.fail("every arrival must lie inside the horizon, after 0.");
  for (double &e : o) e = na(e);
  static const char *const on[] = {"cases", "pending", "rate", "wait", "backlog", "occup", "resid"};
  out = vobj();
  put(out, "est", nvec(7, on, o));
  put(out, "horizon", vnum(h));
  Val need;
  if (c.has("target")) {
    double tg = 0.0;
    CK(binr(c, c["target"], "target", 0.0, Inf, tg));
    need = vnum(tg / o[2]);
  }
  put(out, "need", need);
  return true;
}

bool dktflw(Ctx &c, Val &out) {
  if (!c.has("recv") && !c.has("mix")) return c.fail("give `recv` and `disp`, or `mix`.");
  int w = 0;
  CK(wholed(c, "win", 12, 1, w));
  double lv = 0.0;
  CK(oned(c, "level", 0.95, 0.0, 1.0, lv));
  const double a = 1 - lv;
  out = vobj();
  Val flow, mixo;
  if (c.has("recv")) {
    Dv rv, dv, pv;
    CK(num(c, c["recv"], "recv", rv, -1, 0.0));
    const std::size_t sn = rv.n();
    const int n = static_cast<int>(sn);
    if (!c.has("disp")) return c.fail("`disp` must come with `recv`.");
    CK(num(c, c["disp"], "disp", dv, n, 0.0));
    const bool hp = c.has("pend");
    if (hp) CK(num(c, c["pend"], "pend", pv, n, 0.0));
    std::vector<std::string> pl;
    CK(plab(c, sn, pl));
    std::vector<double> ow(8 * sn);
    const double *o = ow.data();
    clab_dktflw(rv.p(), dv.p(), hp ? pv.p() : nullptr, n, w, ow.data());
    Frame f;
    f.str("period", pl);
    f.num("recv", rv.v);
    f.num("disp", dv.v);
    std::vector<double> cl(sn), cl0(sn), ch0(sn), nt(sn), rc(sn), rl(sn), rh(sn);
    for (std::size_t i = 0; i < sn; ++i) {
      cl[i] = na(o[i]);
      rci(dv[i], rv[i], a / 2, cl0[i], ch0[i]);
      nt[i] = o[sn + i];
      rc[i] = na(o[4 * sn + i]);
      rci(o[3 * sn + i], o[2 * sn + i], a / 2, rl[i], rh[i]);
    }
    f.num("clear", cl);
    f.num("lo", cl0);
    f.num("hi", ch0);
    f.num("net", nt);
    f.num("rclear", rc);
    f.num("rlo", rl);
    f.num("rhi", rh);
    if (hp) {
      f.num("pend", pv.v);
      std::vector<double> pc(sn), tc(sn);
      for (std::size_t i = 0; i < sn; ++i) {
        pc[i] = na(o[5 * sn + i]);
        tc[i] = na(o[6 * sn + i]);
      }
      f.num("pchg", pc);
      f.num("toclear", tc);
    } else {
      f.num("cumnet", o + 7 * sn, sn);
    }
    if (c.has("span")) {
      if (!hp) return c.fail("`span` needs `pend`.");
      std::vector<double> sp;
      CK(rec1(c, c["span"], "span", sn, sp));
      std::vector<double> fr(sn), im(sn);
      for (std::size_t i = 0; i < sn; ++i) {
        fr[i] = rv[i] / sp[i];
        im[i] = fr[i] > 0 ? pv[i] / fr[i] : NA;
      }
      f.num("frate", fr);
      f.num("impl", im);
      if (c.has("mtime")) {
        Dv mt;
        CK(num(c, c["mtime"], "mtime", mt, n, 0.0, true));
        std::vector<double> ra(sn), lb(sn);
        for (std::size_t i = 0; i < sn; ++i) {
          ra[i] = std::isnan(im[i]) ? NA : im[i] / mt[i];
          lb[i] = fr[i] * mt[i];
        }
        f.num("ratio", ra);
        f.num("ltlbk", lb);
      }
    } else if (c.has("mtime")) {
      return c.fail("`mtime` needs `span` and `pend`.");
    }
    if (c.has("shares")) {
      const Val &sh = c["shares"];
      int dim[2] = {-1, -1};
      if (!(sh.t == Val::ARR && !sh.a.empty() && sh.a[0].t == Val::ARR && shape(sh, 2, dim, 0)) || dim[0] != n) {
        return c.fail("`shares` must be a numeric matrix with one row per period, %d.", n);
      }
      const int k = dim[1];
      Dv fl;
      CK(num(c, c["floors"], "floors", fl, k, 0.0));
      std::vector<double> sv(sn * static_cast<std::size_t>(k));
      bool bad = false;
      fill(sh, 2, dim, 0, 0, 1, sv, bad, true);
      for (std::size_t j = 1; j < static_cast<std::size_t>(k); ++j) {
        if (fl[j] < fl[j - 1]) return c.fail("`floors` must be increasing.");
      }
      std::vector<double> bf(sn);
      for (std::size_t i = 0; i < sn; ++i) {
        double s = 0, sw = 0;
        for (std::size_t j = 0; j < static_cast<std::size_t>(k); ++j) {
          const double v = sv[i + sn * j];
          if (!(v >= 0) || !std::isfinite(v)) return c.fail("`shares` must be numbers from 0.");
          s += v * fl[j];
          sw += v;
        }
        if (std::fabs(sw - 1) > 0.01) return c.fail("each row of `shares` must sum to 1 (within rounding).");
        bf[i] = s / sw;
      }
      f.num("floorbd", bf);
    }
    flow = f.o;
  }
  if (c.has("mix")) {
    const Val &mx = c["mix"];
    int dim[2] = {-1, -1};
    if (!(mx.t == Val::ARR && !mx.a.empty() && mx.a[0].t == Val::ARR && shape(mx, 2, dim, 0))) {
      return c.fail("`mix` must be a numeric matrix of periods by outcomes.");
    }
    const int nt = dim[0], k = dim[1];
    if (nt < 1 || k < 1) return c.fail("`mix` must have a period and an outcome.");
    const std::size_t st = static_cast<std::size_t>(nt), sk = static_cast<std::size_t>(k);
    std::vector<double> xv(st * sk);
    bool bad = false;
    fill(mx, 2, dim, 0, 0, 1, xv, bad, true);
    Dv sc;
    if (c.has("scores")) {
      CK(num(c, c["scores"], "scores", sc, nt));
    } else {
      for (int t = 0; t < nt; ++t) sc.v.push_back(t + 1);
    }
    const std::size_t cc = st * sk;
    std::vector<double> mw(3 * cc + 3 + 2 * sk);
    double *sh = mw.data(), *lo = sh + cc, *hi = lo + cc, *s3 = hi + cc, *tr = s3 + 3;
    if (clab_dkmix(xv.data(), nt, k, sc.p(), clab_qchisq(1 - a / k, 1.0, 1), sh, lo, hi, s3, tr) != CLAB_OK) {
      return c.fail("`mix` must hold counts from 0, every period with a case.");
    }
    std::vector<std::string> fp(cc), fo(cc), go(sk);
    for (std::size_t j = 0; j < sk; ++j) {
      go[j] = "outcome" + std::to_string(j + 1);
      for (std::size_t t = 0; t < st; ++t) {
        fp[t + st * j] = std::to_string(t + 1);
        fo[t + st * j] = go[j];
      }
    }
    Frame f;
    f.str("period", fp);
    f.str("outcome", fo);
    f.num("count", xv);
    f.num("share", sh, cc);
    f.num("lo", lo, cc);
    f.num("hi", hi, cc);
    std::vector<double> gs(sk), gp(sk), gh(sk), gc(sk), pw(2 * sk + 1);
    std::vector<int> iw(sk + 1);
    int m = 0;
    for (std::size_t j = 0; j < sk; ++j) {
      gs[j] = na(tr[j]);
      gp[j] = std::isnan(gs[j]) ? NA : clab_pchisq(gs[j], 1.0, 0);
      gc[j] = tr[sk + j];
      if (!std::isnan(gp[j])) pw[static_cast<std::size_t>(m++)] = gp[j];
    }
    clab_pbadj(pw.data(), m, 2, pw.data() + sk, iw.data());
    for (std::size_t j = 0, i = 0; j < sk; ++j) gh[j] = std::isnan(gp[j]) ? NA : pw[sk + i++];
    Frame g;
    g.str("outcome", go);
    g.num("stat", gs);
    g.num("p", gp);
    g.num("holm", gh);
    g.num("change", gc);
    const double hv[4] = {s3[0], s3[1], nt > 1 ? clab_pchisq(s3[0], s3[1], 0) : NA, na(s3[2])};
    static const char *const hn[] = {"chisq", "df", "p", "cramerv"};
    mixo = vobj();
    put(mixo, "shares", f.o);
    put(mixo, "homog", nvec(4, hn, hv));
    put(mixo, "trend", g.o);
  }
  put(out, "flow", flow);
  put(out, "mix", mixo);
  put(out, "level", vnum(lv));
  return true;
}

// A date "YYYY-MM-DD" or "YYYY/MM/DD" (as.Date()'s formats; the first
// that reads the first value is used for all) as days from 1970-01-01.
bool rdate(const std::string &s, char sep, double &d) {
  const char *p = s.c_str();
  while (*p == ' ') ++p;
  char *e = nullptr;
  const long y = std::strtol(p, &e, 10);
  if (e == p || *e != sep) return false;
  p = e + 1;
  const long m = std::strtol(p, &e, 10);
  if (e == p || *e != sep || m < 1 || m > 12) return false;
  p = e + 1;
  const long dd = std::strtol(p, &e, 10);
  if (e == p || dd < 1) return false;
  static const int ml[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  const bool leap = (y % 4 == 0 && y % 100 != 0) || y % 400 == 0;
  if (dd > ml[m - 1] + (m == 2 && leap)) return false;
  // days from civil (H. Hinnant's algorithm)
  const long yy = y - (m <= 2);
  const long era = (yy >= 0 ? yy : yy - 399) / 400;
  const long yoe = yy - era * 400;
  const long doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + dd - 1;
  const long doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
  d = static_cast<double>(era * 146097 + doe - 719468);
  return true;
}
bool rdates(Ctx &c, const char *nm, std::vector<double> &out) {
  const Val &v = c[nm];
  const std::string bad = "stkflw: `from` and `to` must both be dates or date strings.";
  if (v.t == Val::NUL || isnum(v)) return failas(c, bad);
  std::vector<std::string> s;
  const bool arr = v.t == Val::ARR || v.t == Val::OBJ;
  for (std::size_t i = 0; i < vlen(v); ++i) {
    const Val &e = arr ? v.a[i] : v;
    if (e.t == Val::STR) s.push_back(e.s);
    else if (e.t == Val::NUL) s.push_back(std::string(1, '\0'));
    else return failas(c, bad);
  }
  char sep = '-';
  double d = 0.0;
  for (const std::string &t : s) {
    if (!t.empty() && t[0] == '\0') continue;
    sep = rdate(t, '-', d) ? '-' : '/';
    break;
  }
  out.clear();
  for (const std::string &t : s) {
    if ((!t.empty() && t[0] == '\0') || !rdate(t, sep, d)) {
      return failas(c, std::string("stkflw: `") + nm + "` holds something that is not a date.");
    }
    out.push_back(d);
  }
  return true;
}

bool stkflw(Ctx &c, Val &out) {
  Val tdays;
  if (c.has("from") || c.has("to")) {
    std::vector<double> a, b;
    CK(rdates(c, "from", a));
    CK(rdates(c, "to", b));
    bool bad = a.size() != b.size();
    for (std::size_t i = 0; !bad && i < a.size(); ++i) bad = b[i] < a[i];
    if (bad) return failas(c, "stkflw: `to` must not precede `from`, one per period.");
    tdays = varr();
    for (std::size_t i = 0; i < a.size(); ++i) tdays.a.push_back(vnum(b[i] - a[i] + 1));
  }
  const bool hpe = c.has("people"), hd = c.has("days"), ha = c.has("adp"), hl = c.has("alos");
  if (hpe + hd + ha + hl != 2 || (hd && ha)) {
    return c.fail("give two of `people`, `days`, `adp` and `alos` (not `days` with `adp`).");
  }
  Dv av;
  if (ha) {
    const Val &ad = c["adp"];
    const bool lst = (ad.t == Val::OBJ || ad.t == Val::ARR) && !ad.a.empty() && ad.a[0].t == Val::ARR;
    if (lst) {
      for (const Val &e : ad.a) {
        Dv d;
        CK(num(c, e, "adp", d, -1, 0.0));
        av.v.push_back(nmean(d.p(), d.n()));
      }
    } else {
      CK(num(c, ad, "adp", av, -1, 0.0));
    }
  }
  const std::size_t sn = hpe ? vlen(c["people"]) : hd ? vlen(c["days"]) : av.n();
  const int n = static_cast<int>(sn);
  Dv pv, dv, lv;
  if (hpe) CK(num(c, c["people"], "people", pv, n, 0.0, true));
  if (hd) CK(num(c, c["days"], "days", dv, n, 0.0));
  if (hl) CK(num(c, c["alos"], "alos", lv, n, 0.0, true));
  if (ha && av.n() != sn) return c.fail("`adp` must have one value per period, %d.", n);
  Val d365;
  std::vector<double> tv, xv, rv;
  CK(rec1(c, tdays.t != Val::NUL ? tdays : arg(c, "t", d365, 365), "t", sn, tv));
  if (c.has("expo")) {
    CK(rec1(c, c["expo"], "expo", sn, xv));
    // the words of R's message name periods
  }
  double pe = 0.0;
  CK(oned(c, "per", 1e5, 0.0, Inf, pe));
  static const char *const bases[] = {"first", "previous"};
  int b = 0;
  CK(pickd(c, "base", bases, 2, 0, b));
  if (c.has("rel") != c.has("early")) return c.fail("`rel` and `early` go together.");
  double ea = NaN;
  if (c.has("rel")) {
    Dv r;
    CK(num(c, c["rel"], "rel", r, n, 0.0));
    rv = r.v;
    CK(binr(c, c["early"], "early", 0.0, Inf, ea));
  }
  std::vector<std::string> pl;
  CK(plab(c, sn, pl));
  std::vector<double> ow(16 * sn + 4);
  const double *o = ow.data();
  clab_stkflw(hpe ? pv.p() : nullptr, hd ? dv.p() : nullptr, ha ? av.p() : nullptr, hl ? lv.p() : nullptr, tv.data(), n,
              b, xv.empty() ? nullptr : xv.data(), pe, rv.empty() ? nullptr : rv.data(), ea, ow.data());
  Frame f;
  f.str("period", pl);
  f.num("t", tv);
  static const char *const cn[] = {"people", "days",     "adp",     "alos",  "peochg", "aloschg", "dayschg", "adpchg",
                                   "pshare", "flowrate", "stkrate", "frchg", "srchg",  "cut",     "cutpct"};
  for (int k = 0; k < 15; ++k) {
    if (k >= 9 && k <= 12 && xv.empty()) continue;
    if (k >= 13 && rv.empty()) continue;
    std::vector<double> d(sn);
    for (std::size_t i = 0; i < sn; ++i) d[i] = na(o[static_cast<std::size_t>(k) * sn + i]);
    f.num(cn[k], d);
  }
  static const char *const mn[] = {"people", "days", "adp", "alos"};
  out = vobj();
  put(out, "flow", f.o);
  put(out, "mean", nvec(4, mn, o + 15 * sn));
  put(out, "base", vstr(bases[b]));
  put(out, "per", vnum(xv.empty() ? NA : pe));
  return true;
}

bool trnmat(Ctx &c, Val &out) {
  const std::size_t sn = vlen(c["state"]);
  const int n = static_cast<int>(sn);
  Lab s, p;
  CK(codes(c, c["state"], "state", sn, s));
  CK(codes(c, c["id"], "id", sn, p));
  Dv tv;
  const bool ht = c.has("time");
  if (ht) CK(num(c, c["time"], "time", tv, n));
  std::vector<int> o(sn);
  for (int i = 0; i < n; ++i) o[static_cast<std::size_t>(i)] = i;
  const int *pc = p.code.data();
  const double *tt = ht ? tv.p() : nullptr;
  std::stable_sort(o.begin(), o.end(), [&](int a, int b) { return pc[a] < pc[b] || (pc[a] == pc[b] && tt && tt[a] < tt[b]); });
  std::vector<int> fr, to;
  for (std::size_t i = 1; i < sn; ++i) {
    const int a = o[i - 1], b = o[i];
    if (pc[a] != pc[b]) continue;
    fr.push_back(s.code[static_cast<std::size_t>(a)]);
    to.push_back(s.code[static_cast<std::size_t>(b)]);
  }
  const int k = s.g;
  if (static_cast<long long>(k) * k > 25000000LL) return c.fail("`state` has more than 5000 states.");
  const std::size_t sk = static_cast<std::size_t>(k);
  std::vector<double> cnt(sk * sk), pm(sk * sk), pv(sk), wk(sk * (sk + 1));
  double st[3] = {0};
  clab_trnmat(fr.data(), to.data(), static_cast<int>(fr.size()), k, cnt.data(), pm.data(), pv.data(), st, wk.data());
  for (double &e : pv) e = na(e);
  for (double &e : st) e = na(e);
  static const char *const smn[] = {"moves", "stay", "theil"};
  out = vobj();
  put(out, "count", vmat(cnt, k, k));
  put(out, "prob", vmat(pm, k, k));
  put(out, "pi", nvec(s.lev, pv.data()));
  put(out, "summ", nvec(3, smn, st));
  return true;
}

// .bands(): band labels read into their ends.
struct Bands {
  std::vector<std::string> lab;
  std::vector<int> isna;
  std::vector<double> lo, hi;
  std::vector<int> olo, oup;
};

// The range "a to b": its ends unless a word or a plus opens it or the
// ends are reversed.
void brange(const std::smatch &m, double *r) {
  static const std::regex cue("\\b(over|above|more|plus|under|less|below|fewer|greater|than|gt|lt)\\b");
  const double a = std::strtod(m[2].str().c_str(), nullptr), b = std::strtod(m[3].str().c_str(), nullptr);
  const std::string pre = m[1].str(), post = m[4].str();
  if (!std::regex_search(pre, cue) && !std::regex_search(post, cue) && post.find('+') == std::string::npos && a <= b) {
    r[0] = a;
    r[1] = b;
  }
}

// .blab(): one normalised label: lower, upper, open below, open above.
void blab(const std::string &z, bool clsup, bool intsca, double *r) {
  r[0] = r[1] = NA;
  r[2] = r[3] = 0;
  if (z.empty()) return;
  static const std::string nm = "(-?[0-9]*\\.?[0-9]+)";
  static const std::regex rng("^([a-z]+ )?" + nm + " ?(?:to|-|through) ?" + nm + "([^0-9]*)$");
  static const std::regex one1("-?[0-9]*\\.?[0-9]+");
  std::smatch m;
  if (std::regex_match(z, m, rng)) {
    brange(m, r);
    return;
  }
  auto it = std::sregex_iterator(z.begin(), z.end(), one1);
  std::vector<std::string> k;
  for (; it != std::sregex_iterator(); ++it) k.push_back(it->str());
  if (k.size() != 1) return;
  const double v = std::strtod(k[0].c_str(), nullptr);
  static const std::regex dec("[0-9]\\.[0-9]");
  const double st = intsca && !std::regex_search(z, dec) ? 1 : 0;
  static const std::regex up1("(\\+\\s*$|and over|and above|or more|or over|plus\\s*$)");
  static const std::regex up2("(greater than|more than|over|above|\\bgt\\b)");
  static const std::regex dn1("(or less|or fewer|and under|and below|or under)");
  static const std::regex dn2("(less than|fewer than|under|below|^<|\\blt\\b)");
  static const std::regex bare("^-?[0-9]*\\.?[0-9]+$");
  if (std::regex_search(z, up1)) {
    r[0] = clsup ? v : v + st;
    r[3] = 1;
  } else if (std::regex_search(z, up2)) {
    r[0] = v + st;
    r[3] = 1;
  } else if (std::regex_search(z, dn1)) {
    r[1] = v;
    r[2] = 1;
  } else if (std::regex_search(z, dn2)) {
    r[1] = v - st;
    r[2] = 1;
  } else if (std::regex_search(z, bare)) {
    r[0] = r[1] = v;
  }
}

// The labels normalised as .bands() does: lower case (ASCII), trimmed,
// dashes and the minus sign as "-", runs of spaces as one, and the commas
// of thousands dropped.
std::string bnorm(const std::string &lab) {
  std::string s;
  for (std::size_t i = 0; i < lab.size(); ++i) {
    const unsigned char ch = static_cast<unsigned char>(lab[i]);
    // U+2010 to U+2015 (E2 80 90 to 95) and U+2212 (E2 88 92)
    if (ch == 0xE2 && i + 2 < lab.size()) {
      const unsigned char b1 = static_cast<unsigned char>(lab[i + 1]), b2 = static_cast<unsigned char>(lab[i + 2]);
      if ((b1 == 0x80 && b2 >= 0x90 && b2 <= 0x95) || (b1 == 0x88 && b2 == 0x92)) {
        s += '-';
        i += 2;
        continue;
      }
    }
    s += static_cast<char>(std::tolower(ch));
  }
  auto ws = [](char ch) { return ch == ' ' || ch == '\t' || ch == '\r' || ch == '\n' || ch == '\v' || ch == '\f'; };
  std::size_t a = 0, e = s.size();
  while (a < e && (s[a] == ' ' || s[a] == '\t' || s[a] == '\r' || s[a] == '\n')) ++a;
  while (e > a && (s[e - 1] == ' ' || s[e - 1] == '\t' || s[e - 1] == '\r' || s[e - 1] == '\n')) --e;
  std::string t;
  for (std::size_t i = a; i < e; ++i) {
    if (ws(s[i])) {
      if (t.empty() || t.back() != ' ') t += ' ';
    } else {
      t += s[i];
    }
  }
  std::string u;
  auto dig = [&t](std::size_t i) { return i < t.size() && t[i] >= '0' && t[i] <= '9'; };
  for (std::size_t i = 0; i < t.size(); ++i) {
    if (t[i] == ',' && i > 0 && dig(i - 1) && dig(i + 1) && dig(i + 2) && dig(i + 3) && !dig(i + 4)) continue;
    u += t[i];
  }
  return u;
}

bool bxok(Ctx &c, const Val &x) {
  const char *m = "bandx: `x` must be band labels: text or a factor.";
  const bool arr = x.t == Val::ARR || x.t == Val::OBJ;
  const std::size_t n = vlen(x);
  bool any = false;
  for (std::size_t i = 0; i < n; ++i) {
    const Val &e = arr ? x.a[i] : x;
    if (e.t == Val::STR) any = true;
    else if (e.t != Val::NUL) return failas(c, m);
  }
  return (n && any) || failas(c, m);
}

bool bands(Ctx &c, const Val &x, bool clsup, bool intsca, Bands &b) {
  CK(bxok(c, x));
  const bool arr = x.t == Val::ARR || x.t == Val::OBJ;
  const std::size_t n = vlen(x);
  b = Bands();
  for (std::size_t i = 0; i < n; ++i) {
    const Val &e = arr ? x.a[i] : x;
    double r[4];
    const bool isn = e.t == Val::NUL;
    if (isn) {
      r[0] = r[1] = NA;
      r[2] = r[3] = 0;
    } else {
      blab(bnorm(e.s), clsup, intsca, r);
    }
    b.lab.push_back(isn ? std::string() : e.s);
    b.isna.push_back(isn);
    b.lo.push_back(r[0]);
    b.hi.push_back(r[1]);
    b.olo.push_back(r[2] == 1);
    b.oup.push_back(r[3] == 1);
  }
  return true;
}
Val blabs(const Bands &b) {
  Val v = varr();
  for (std::size_t i = 0; i < b.lab.size(); ++i) v.a.push_back(b.isna[i] ? Val() : vstr(b.lab[i]));
  return v;
}
// toString() of the labels where bad.
std::string blist(const Bands &b, const std::vector<int> &bad) {
  std::string s;
  for (std::size_t i = 0; i < bad.size(); ++i) {
    if (!bad[i]) continue;
    if (!s.empty()) s += ", ";
    s += b.isna[i] ? "NA" : b.lab[i];
  }
  return s;
}

bool bflag(Ctx &c, const char *nm, bool &v) {
  if (!c.has(nm)) {
    v = true;
    return true;
  }
  const Val &o = c[nm];
  const Val *e = o.t == Val::ARR && o.a.size() == 1 ? &o.a[0] : &o;
  if (e->t != Val::BOOL) return failas(c, std::string("bandx: `") + nm + "` must be TRUE or FALSE.");
  v = e->b;
  return true;
}

bool bandx(Ctx &c, Val &out) {
  bool clsup = true, intsca = true;
  CK(bxok(c, c["x"]));
  CK(bflag(c, "clsup", clsup));
  CK(bflag(c, "intsca", intsca));
  Bands b;
  CK(bands(c, c["x"], clsup, intsca, b));
  const std::size_t sn = b.lab.size();
  const int n = static_cast<int>(sn);
  static const char *const rules[] = {"midpoint", "lower", "upper", "geometric"};
  int r = 0;
  CK(pickd(c, "rule", rules, 4, 0, r));
  double cp = NaN, fc = 0.0, fl = NaN;
  if (c.has("cap")) CK(binr(c, c["cap"], "cap", -Inf, Inf, cp));
  CK(oned(c, "fac", 2, 0.0, Inf, fc));
  if (c.has("floor")) CK(binr(c, c["floor"], "floor", -Inf, Inf, fl));
  std::vector<double> v(sn), bl(sn), bh(sn);
  if (clab_bandx(b.lo.data(), b.hi.data(), b.olo.data(), b.oup.data(), n, r, cp, fc, fl, v.data(), bl.data(),
                 bh.data()) != CLAB_OK) {
    return c.fail("the open band's assumed end (`cap` or `floor`) lies on the wrong side of its other end.");
  }
  for (double &e : v) e = na(e);
  Frame f;
  f.add("label", blabs(b));
  f.num("lower", b.lo);
  f.num("upper", b.hi);
  f.lgl("openlo", b.olo.data(), sn);
  f.lgl("openup", b.oup.data(), sn);
  f.num("value", v);
  f.num("lo", bl);
  f.num("hi", bh);
  Val vals;
  if (c.has("cnts")) {
    const Val &cn = c["cnts"];
    std::vector<double> cv = dbl(cn);
    bool ok = isnum(cn) && cv.size() == sn;
    for (double e : cv) ok = ok && !std::isnan(e) && e >= 0;
    if (!ok) return failas(c, "bandx: `cnts` must be numbers from 0, one per band, " + std::to_string(sn) + ".");
    std::vector<int> lost(sn);
    bool anyl = false;
    for (std::size_t i = 0; i < sn; ++i) anyl |= (lost[i] = std::isnan(v[i]) && cv[i] > 0) != 0;
    bool drop = false;
    if (c.has("drop")) {
      const Val &d = c["drop"];
      drop = (d.t == Val::BOOL && d.b) || (d.t == Val::ARR && d.a.size() == 1 && d.a[0].t == Val::BOOL && d.a[0].b);
    }
    if (anyl && !drop) {
      return failas(c, "bandx: these bands cannot be read: " + blist(b, lost) + " (drop = TRUE leaves out their units).");
    }
    std::vector<double> e;
    for (std::size_t i = 0; i < sn; ++i) {
      if (std::isnan(v[i])) continue;
      const double k = std::nearbyint(cv[i]);
      for (double q = 0; q < k; ++q) e.push_back(v[i]);
    }
    vals = vauto(e);
  }
  out = vobj();
  put(out, "bands", f.o);
  put(out, "vals", vals);
  put(out, "sens", Val());
  put(out, "span", Val());
  put(out, "relspan", Val());
  return true;
}

// .cells(): published cells read: the value (NA for none), the status (0
// a number, 1 suppressed, 2 missing) and a suppressed cell's "<k" limit.
bool cells(Ctx &c, const Val &x, const std::vector<std::string> &sup, const char *nm, std::vector<double> &v,
           std::vector<int> &st, std::vector<double> &lim) {
  const std::size_t n = vlen(x);
  const bool arr = x.t == Val::ARR || x.t == Val::OBJ;
  bool text = false, okt = true;
  for (std::size_t i = 0; i < n; ++i) {
    const Val &e = arr ? x.a[i] : x;
    if (e.t == Val::STR) text = true;
    else if (e.t != Val::NUM && e.t != Val::NUL) okt = false;
  }
  if (!n || !okt) return failas(c, std::string("pbbnd: `") + nm + "` must be numbers or text cells.");
  v.assign(n, NA);
  st.assign(n, 0);
  lim.assign(n, NA);
  if (!text) {
    v = dbl(x);
  } else {
    static const std::regex lt("^<\\s*([0-9]+)$");
    std::vector<std::string> bad;
    for (std::size_t i = 0; i < n; ++i) {
      const Val &e = arr ? x.a[i] : x;
      if (e.t == Val::NUL) continue;
      std::string s = e.t == Val::STR ? e.s : rfmt(e.n);
      std::size_t a = 0, z = s.size();
      while (a < z && (s[a] == ' ' || s[a] == '\t' || s[a] == '\r' || s[a] == '\n')) ++a;
      while (z > a && (s[z - 1] == ' ' || s[z - 1] == '\t' || s[z - 1] == '\r' || s[z - 1] == '\n')) --z;
      s = s.substr(a, z - a);
      std::smatch m;
      const bool islt = std::regex_match(s, m, lt);
      if (islt) lim[i] = std::strtod(m[1].str().c_str(), nullptr);
      if (islt || std::find(sup.begin(), sup.end(), s) != sup.end()) st[i] = 1;
      if (st[i] != 0 || s.empty()) continue;
      std::string t;
      for (char ch : s) {
        if (ch != ',') t += ch;
      }
      char *end = nullptr;
      const double d = std::strtod(t.c_str(), &end);
      while (end && (*end == ' ' || *end == '\t')) ++end;
      if (t.empty() || end == t.c_str() || *end != '\0') {
        if (std::find(bad.begin(), bad.end(), s) == bad.end()) bad.push_back(s);
      } else {
        v[i] = d;
      }
    }
    if (!bad.empty()) {
      std::string l;
      for (const std::string &s : bad) l += (l.empty() ? "" : ", ") + s;
      return failas(c, "pbbnd: cannot read these cells as numbers or suppression marks: " + l + ".");
    }
  }
  for (std::size_t i = 0; i < n; ++i) {
    if (!text && std::isinf(v[i])) return failas(c, std::string("pbbnd: `") + nm + "` must be finite.");
    if (st[i] == 0 && std::isnan(v[i])) st[i] = 2;
  }
  return true;
}

bool pbbnd(Ctx &c, Val &out) {
  std::vector<std::string> sup = {"x", "X", "s", "..", "F", "suppressed", "n/a"};
  if (c.has("supmar")) {
    const Val &s = c["supmar"];
    sup.clear();
    const bool arr = s.t == Val::ARR;
    for (std::size_t i = 0; i < vlen(s); ++i) {
      const Val &e = arr ? s.a[i] : s;
      if (e.t != Val::STR) return failas(c, "pbbnd: `supmar` must be text markers.");
      sup.push_back(e.s);
    }
  }
  std::vector<double> val, lim, pvv, plim;
  std::vector<int> st, pst;
  CK(cells(c, c["x"], sup, "x", val, st, lim));
  const bool hp = c.has("prev");
  if (hp) CK(cells(c, c["prev"], sup, "prev", pvv, pst, plim));
  const std::size_t sn = val.size();
  const int n = static_cast<int>(sn);
  double rn = NaN, sp = NaN;
  CK(optd(c, "rnd", 0.0, rn));
  CK(optd(c, "suplim", 0.0, sp));
  if (rn == 0) return c.fail("`rnd` must be above 0.");
  int kd = 0;
  CK(pickd(c, "kind", rkinds, 2, 0, kd));
  for (double e : val) {
    if (e < 0) return c.fail("`x` must hold counts from 0.");
  }
  std::vector<double> lo(sn), hi(sn);
  std::vector<int> so(sn);
  if (clab_pbbnd(val.data(), st.data(), lim.data(), n, rn, kd, sp, lo.data(), hi.data(), so.data()) != CLAB_OK) {
    return c.fail("a suppressed cell has no limit: give `suplim` or write the cell as \"<k\".");
  }
  static const char *const sn2[] = {"exact", "rounded", "suppressed", "missing"};
  std::vector<std::string> ss(sn);
  for (std::size_t i = 0; i < sn; ++i) {
    ss[i] = sn2[so[i]];
    lo[i] = na(lo[i]);
    hi[i] = na(hi[i]);
  }
  Frame f;
  f.num("value", val);
  f.num("lower", lo);
  f.num("upper", hi);
  f.str("status", ss);
  Val chg;
  if (hp) {
    if (pvv.size() != sn) return c.fail("`prev` must have one cell per cell of `x`, %d.", n);
    for (double e : pvv) {
      if (e < 0) return c.fail("`prev` must hold counts from 0.");
    }
    std::vector<double> w(10 * sn);
    std::vector<int> so2(sn);
    double *plo = w.data(), *phi = plo + sn, *o = phi + sn;
    if (clab_pbbnd(pvv.data(), pst.data(), plim.data(), n, rn, kd, sp, plo, phi, so2.data()) != CLAB_OK) {
      return c.fail("a suppressed cell of `prev` has no limit: give `suplim` or write it as \"<k\".");
    }
    Dv pv, qv;
    const bool hpo = c.has("pop"), hqo = c.has("prepop");
    if (hpo) {
      CK(num(c, c["pop"], "pop", pv, -1, 0.0, true));
      if (pv.n() != 1 && pv.n() != sn) return c.fail("`pop` must be one number or one per cell, %d.", n);
    }
    if (hqo) {
      CK(num(c, c["prepop"], "prepop", qv, -1, 0.0, true));
      if (qv.n() != 1 && qv.n() != sn) return c.fail("`prepop` must be one number or one per cell, %d.", n);
    }
    std::vector<double> pr(2 * sn);
    for (std::size_t i = 0; i < sn; ++i) {
      if (hpo) pr[i] = pv[pv.n() == 1 ? 0 : i];
      if (hqo) pr[sn + i] = qv[qv.n() == 1 ? 0 : i];
    }
    double pe = 0.0;
    CK(oned(c, "per", 1000, 0.0, Inf, pe));
    clab_pbchg(lo.data(), hi.data(), plo, phi, n, hpo ? pr.data() : nullptr, hqo ? pr.data() + sn : nullptr, pe, o);
    static const char *const cn[] = {"chglo", "chgup", "pctlo", "pctup", "ratelo", "rateup", "preratelo", "prerateup"};
    const int nc = !hpo ? 4 : !hqo ? 6 : 8;
    Frame g;
    for (int k = 0; k < nc; ++k) {
      std::vector<double> d(sn);
      for (std::size_t i = 0; i < sn; ++i) d[i] = na(o[static_cast<std::size_t>(k) * sn + i]);
      g.num(cn[k], d);
    }
    chg = g.o;
  }
  out = vobj();
  put(out, "bnd", f.o);
  put(out, "chg", chg);
  return true;
}

bool lftab(Ctx &c, Val &out) {
  const Val &br = c["breaks"];
  std::vector<double> lo, hi;
  const bool txt = !isnum(br) && (br.t == Val::STR || br.t == Val::ARR);
  if (txt) {
    Bands b;
    CK(bands(c, br, true, false, b));
    std::vector<int> bad(b.lab.size());
    bool any = false;
    for (std::size_t i = 0; i < b.lab.size(); ++i) {
      lo.push_back(b.olo[i] ? 0.0 : b.lo[i]);
      hi.push_back(b.oup[i] ? Inf : b.hi[i]);
      bad[i] = std::isnan(lo[i]) || std::isnan(hi[i]);
      any = any || bad[i];
    }
    if (any) return failas(c, "lftab: these bands cannot be read: " + blist(b, bad) + ".");
  } else {
    std::vector<double> e = dbl(br);
    bool ok = isnum(br) && e.size() >= 2;
    for (double v : e) ok = ok && !std::isnan(v);
    if (!ok) return failas(c, "lftab: `breaks` must be at least two band edges.");
    lo.assign(e.begin(), e.end() - 1);
    hi.assign(e.begin() + 1, e.end());
  }
  const std::size_t sk = lo.size();
  const int k = static_cast<int>(sk);
  Dv dv;
  CK(num(c, c["events"], "events", dv, k, 0.0));
  std::vector<double> cv;
  {
    Dv c0;
    Val d0;
    CK(num(c, arg(c, "cens", d0, 0), "cens", c0, -1, 0.0));
    if (c0.n() != 1 && c0.n() != sk) return c.fail("`cens` must be one number or one per band, %d.", k);
    for (std::size_t i = 0; i < sk; ++i) cv.push_back(c0[c0.n() == 1 ? 0 : i]);
  }
  double lv = 0.0;
  CK(oned(c, "level", 0.95, 0.0, 1.0, lv));
  static const char *const cis[] = {"plain", "log", "log-log"};
  int ci = 0;
  CK(pickd(c, "ci", cis, 3, 2, ci));
  std::vector<double> ow(8 * sk + 1);
  if (clab_lftab(lo.data(), hi.data(), dv.p(), cv.data(), k, clab_qnorm(0.5 + lv / 2), ci, ow.data()) != CLAB_OK) {
    return c.fail("the bands must be contiguous, only the last open, with cases in them.");
  }
  const double *o = ow.data();
  static const char *const cn[] = {"entering", "atrisk", "q", "surv", "se", "lo", "hi", "hazard"};
  Frame f;
  f.num("lower", lo);
  f.num("upper", hi);
  f.num(cn[0], o, sk);
  f.num(cn[1], o + sk, sk);
  f.num("events", dv.v);
  f.num("cens", cv);
  for (std::size_t j = 2; j < 8; ++j) {
    std::vector<double> d(sk);
    for (std::size_t i = 0; i < sk; ++i) d[i] = na(o[j * sk + i]);
    f.num(cn[j], d);
  }
  out = vobj();
  put(out, "tab", f.o);
  put(out, "median", vnum(o[8 * sk]));
  put(out, "ci", vstr(cis[ci]));
  return true;
}
