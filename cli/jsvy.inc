// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Survey estimates (src/bsvy.cpp): means, totals, proportions, ratios and
// variances under strata, clusters and weights, by domain.

bool svymx(Ctx &c, Val &out) {
  static const char *const so[] = {"mean", "total", "prop", "ratio", "var"};
  int sk = 0;
  CK(pickd(c, "stat", so, 5, 0, sk));
  double lv = 0.0;
  CK(oned(c, "level", 0.95, 0.0, 1.0, lv));
  CK(isvec(c, c["y"], "y"));
  Dv yv;
  CK(num(c, c["y"], "y", yv));
  const int n = static_cast<int>(yv.n());
  const std::size_t sn = yv.n();
  if (n < 2) return c.fail("`y` must have two values at least.");
  const double *yy = yv.p();
  if (sk == 2) {
    for (double e : yv.v) {
      if (e != 0 && e != 1) return c.fail("`stat = \"prop\"` needs `y` of 0 and 1.");
    }
  }
  std::vector<double> wv(sn, 1.0), xv;
  if (c.has("w")) {
    Dv wx;
    CK(num(c, c["w"], "w", wx, n, 0.0, true));
    wv = wx.v;
  }
  if (sk == 3) {
    if (!c.has("den")) return c.fail("`stat = \"ratio\"` needs `den`.");
    Dv dx;
    CK(num(c, c["den"], "den", dx, n));
    xv = dx.v;
  } else if (c.has("den")) {
    return c.fail("`den` is the denominator of `stat = \"ratio\"`.");
  }
  std::vector<int> st(sn, 0), ps(sn);
  int nst = 1;
  if (c.has("strata")) {
    Lab l;
    CK(codes(c, c["strata"], "strata", sn, l));
    st = l.code;
    nst = l.g;
  }
  std::vector<int> pu(sn);
  int npu = n;
  if (!c.has("psu")) {
    for (int i = 0; i < n; ++i) pu[zz(i)] = i;
  } else {
    Lab l;
    CK(codes(c, c["psu"], "psu", sn, l));
    pu = l.code;
    npu = l.g;
  }
  std::vector<int> iw(zz(nst) * zz(npu));
  const int nps = pev::pairs(st.data(), pu.data(), n, nst, npu, ps.data(), iw.data());
  std::vector<double> nh(zz(nst), 0.0);
  {
    std::vector<int> sp(zz(nps), -1);
    for (std::size_t i = 0; i < sn; ++i) sp[zz(ps[i])] = st[i];
    for (int a : sp) nh[zz(a)] += 1.0;
  }
  for (double v : nh) {
    if (v < 2) return c.fail("a stratum has one first-stage unit: its variance cannot be estimated.");
  }
  std::vector<double> fv;
  if (c.has("fpc")) {
    Dv fx;
    CK(num(c, c["fpc"], "fpc", fx, n, 0.0, true));
    fv.assign(zz(nst), -1.0);
    for (std::size_t i = 0; i < sn; ++i) {
      double &f = fv[zz(st[i])];
      if (f >= 0 && f != fx[i]) return c.fail("`fpc` must be constant within each stratum.");
      f = fx[i];
    }
    for (int h = 0; h < nst; ++h) {
      if (fv[zz(h)] < nh[zz(h)]) return c.fail("`fpc` must be at least the first-stage units sampled in its stratum.");
    }
  }
  std::vector<int> dm(sn, 0);
  int nd = 1;
  std::vector<std::string> dlev;
  const bool hb = c.has("by");
  if (hb) {
    Lab l;
    CK(codes(c, c["by"], "by", sn, l));
    dm = l.code;
    nd = l.g;
    dlev = l.lev;
  }
  int r0 = 1;
  if (c.has("ref") && hb) CK(rref(c, c["ref"], dlev, "ref", "domain", r0));
  --r0;
  const double df = nps - nst, tq = clab_qt(0.5 + 0.5 * lv, df, 1);
  std::vector<double> th(zz(nd)), u(zz(nd) * sn, 0.0), cnt(zz(nd), 0.0);
  for (int d = 0; d < nd; ++d) {
    double sw = 0, swy = 0, swx = 0, k = 0;
    for (std::size_t i = 0; i < sn; ++i) {
      if (dm[i] != d) continue;
      sw += wv[i];
      swy += wv[i] * yy[i];
      if (sk == 3) swx += wv[i] * xv[i];
      k += 1;
    }
    cnt[zz(d)] = k;
    double *ud = u.data() + zz(d) * sn;
    double t = 0.0;
    if (sk == 1) {
      t = swy;
      for (std::size_t i = 0; i < sn; ++i) ud[i] = dm[i] == d ? wv[i] * yy[i] : 0.0;
    } else if (sk == 3) {
      if (swx == 0) return c.fail("`den` sums to 0 in a domain: the ratio is undefined.");
      t = swy / swx;
      for (std::size_t i = 0; i < sn; ++i) ud[i] = dm[i] == d ? wv[i] * (yy[i] - t * xv[i]) / swx : 0.0;
    } else {
      const double m = swy / sw;
      std::vector<double> zv(sn);
      for (std::size_t i = 0; i < sn; ++i) zv[i] = sk == 4 ? (yy[i] - m) * (yy[i] - m) * k / (k - 1) : yy[i];
      if (sk == 4 && k < 2) return c.fail("a domain has one unit: its variance is undefined.");
      for (std::size_t i = 0; i < sn; ++i) t += dm[i] == d ? wv[i] * zv[i] / sw : 0.0;
      for (std::size_t i = 0; i < sn; ++i) ud[i] = dm[i] == d ? wv[i] * (zv[i] - t) / sw : 0.0;
    }
    th[zz(d)] = t;
  }
  const double *fp = fv.empty() ? nullptr : fv.data();
  std::vector<std::string> cd;
  std::vector<double> cs, cl, ch, cv, cdf;
  for (int d = 0; d < nd; ++d) {
    const double se = std::sqrt(pev::svar(u.data() + zz(d) * sn, n, st.data(), nst, ps.data(), nps, fp));
    cd.push_back(hb ? dlev[zz(d)] : "all");
    cs.push_back(se);
    cl.push_back(th[zz(d)] - tq * se);
    ch.push_back(th[zz(d)] + tq * se);
    cv.push_back(th[zz(d)] != 0 ? 100 * se / std::fabs(th[zz(d)]) : NA);
    cdf.push_back(df);
  }
  Frame est;
  est.str("dom", cd);
  est.num("est", th);
  est.num("se", cs);
  est.num("lo", cl);
  est.num("hi", ch);
  est.num("cv", cv);
  est.num("n", cnt);
  est.num("df", cdf);
  Val cmp;
  if (nd > 1) {
    std::vector<std::string> kd;
    std::vector<double> kf, ks, kt, kp, ud(sn);
    for (int d = 0; d < nd; ++d) {
      if (d == r0) continue;
      for (std::size_t i = 0; i < sn; ++i) ud[i] = u[zz(d) * sn + i] - u[zz(r0) * sn + i];
      const double se = std::sqrt(pev::svar(ud.data(), n, st.data(), nst, ps.data(), nps, fp));
      kd.push_back(dlev[zz(d)]);
      kf.push_back(th[zz(d)] - th[zz(r0)]);
      ks.push_back(se);
      kt.push_back(kf.back() / se);
      kp.push_back(2 * clab_pt(-std::fabs(kt.back()), df, 1));
    }
    Frame k;
    k.str("dom", kd);
    k.num("diff", kf);
    k.num("se", ks);
    k.num("stat", kt);
    k.num("p", kp);
    cmp = k.o;
  }
  double sw = 0, sw2 = 0;
  for (double v : wv) {
    sw += v;
    sw2 += v * v;
  }
  Frame g;
  const double gv[7] = {static_cast<double>(n), static_cast<double>(nst), static_cast<double>(nps), sw, df,
                        sw * sw / sw2, n * sw2 / (sw * sw)};
  static const char *const gn[] = {"n", "strata", "psu", "sumw", "df", "effn", "deffw"};
  for (int j = 0; j < 7; ++j) g.num(gn[j], gv + j, 1);
  out = vobj();
  put(out, "est", est.o);
  put(out, "cmp", cmp);
  put(out, "des", g.o);
  put(out, "stat", vstr(so[sk]));
  put(out, "level", vnum(lv));
  return true;
}
