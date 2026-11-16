// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Eyewitness lineups (src/blnr.cpp): ROC curves and lineup fairness.

bool lnroc(Ctx &c, Val &out) {
  const std::size_t n = vlen(c["tp"]);
  std::vector<int> tv;
  CK(bin01(c, c["tp"], "tp", n, tv));
  if (n < 2) return c.fail("`tp` must hold two witnesses at least.");
  Lab pl;
  CK(codes(c, c["pick"], "pick", n, pl));
  std::vector<int> pk(n), map;
  for (const std::string &s : pl.lev) {
    const int v = s == "none" ? 0 : s == "suspect" ? 1 : s == "filler" ? 2 : -1;
    if (v < 0) {
      return c.fail("`pick` must be \"suspect\", \"filler\" or \"none\" for every witness, not %s.", s.c_str());
    }
    map.push_back(v);
  }
  for (std::size_t i = 0; i < n; ++i) pk[i] = map[zz(pl.code[i])];
  Dv cv;
  if (readd(c["conf"], cv) != NUM || cv.n() != n) {
    return c.fail("`conf` must be numbers, one per witness, %d.", static_cast<int>(n));
  }
  for (std::size_t i = 0; i < n; ++i) {
    if (pk[i] != 0 && !std::isfinite(cv[i])) return c.fail("`conf` must be a finite number for every identification.");
  }
  Lab pr;
  int np = 1;
  const bool hp = c.has("proc");
  if (hp) {
    CK(codes(c, c["proc"], "proc", n, pr));
    np = pr.g;
  } else {
    pr.code.assign(n, 0);
  }
  int sz = 0;
  if (c.has("size")) CK(whole(c, c["size"], "size", 2, sz));
  double fm = 0.0;
  if (c.has("far")) {
    CK(binr(c, c["far"], "far", 0.0, 1.0, fm));
    if (!(fm > 0)) return c.fail("`far` must be one number above 0, at most 1.");
  }
  int nb = 0;
  CK(wholed(c, "nboot", 0, 0, nb));
  unsigned int sd = 0;
  if (nb > 0) CK(seed(c, sd));
  double lv = 0.0;
  CK(oned(c, "level", 0.95, 0.0, 1.0, lv));
  {
    std::vector<int> t1(zz(np), 0), t0(zz(np), 0);
    for (std::size_t i = 0; i < n; ++i) (tv[i] ? t1 : t0)[zz(pr.code[i])]++;
    for (int p = 0; p < np; ++p) {
      if (!t1[zz(p)] || !t0[zz(p)]) return c.fail("every procedure needs target-present and target-absent lineups.");
    }
    bool any = false;
    for (int v : pk) any = any || v != 0;
    if (!any) return c.fail("no witness made an identification.");
  }
  std::vector<double> lev(n), roc(3 * zz(np) * n), cac(3 * zz(np) * n), pa(zz(np)), bse(2 * zz(np));
  int nl = 0;
  double fu = 0.0;
  const int rc = clab_lnroc(tv.data(), pk.data(), cv.p(), pr.code.data(), static_cast<int>(n), np, sz, fm, nb, sd,
                            lev.data(), &nl, roc.data(), cac.data(), pa.data(), bse.data(), &fu);
  if (rc == CLAB_EINTR) return kchk(c, rc);
  if (rc != CLAB_OK) return c.fail("the false alarm rates are all 0: no partial area.");
  const std::size_t m = zz(np) * zz(nl), un = zz(np);
  out = vobj();
  for (int t = 0; t < 2; ++t) {
    Frame f;
    if (hp) {
      std::vector<int> pc;
      for (int p = 0; p < np; ++p) {
        for (int l = 0; l < nl; ++l) pc.push_back(p);
      }
      f.add("proc", vstrs(pr.lev, pc));
    }
    std::vector<double> cf;
    for (std::size_t k = 0; k < m; ++k) cf.push_back(lev[k % zz(nl)]);
    f.num("conf", cf);
    const double *src = t == 0 ? roc.data() : cac.data();
    static const char *const rn[] = {"hr", "far", "dr"}, *const cn[] = {"hr", "far", "acc"};
    for (int q = 0; q < 3; ++q) f.num(t == 0 ? rn[q] : cn[q], nav(src + zz(q) * m, m));
    put(out, t == 0 ? "roc" : "cac", f.o);
  }
  const double z = qnorm(0.5 + 0.5 * lv);
  Frame s;
  if (hp) s.str("proc", pr.lev);
  std::vector<double> se(un), lo(un), hi(un), df(un), ds(un), pv(un);
  for (std::size_t p = 0; p < un; ++p) {
    se[p] = na(bse[p]);
    lo[p] = na(pa[p] - z * bse[p]);
    hi[p] = na(pa[p] + z * bse[p]);
    df[p] = pa[p] - pa[0];
    ds[p] = na(bse[un + p]);
    const double zv = df[p] / bse[un + p];
    pv[p] = std::isnan(zv) ? NA : 2 * clab_pnorm(-std::fabs(zv));
  }
  s.num("pauc", pa);
  s.num("se", se);
  s.num("lo", lo);
  s.num("hi", hi);
  s.num("diff", df);
  s.num("dse", ds);
  s.num("p", pv);
  put(out, "pauc", s.o);
  put(out, "far", vnum(fu));
  put(out, "size", sz ? vnum(sz) : Val());
  put(out, "nboot", vnum(nb));
  return true;
}

bool lnsize(Ctx &c, Val &out) {
  Dv cv;
  CK(num(c, c["choice"], "choice", cv, -1, 0.0));
  int mx = 0;
  for (double v : cv.v) {
    if (v != std::floor(v) || v > INT_MAX) {
      return c.fail("`choice` must be whole numbers from 0: lineup positions, 0 for no choice.");
    }
    mx = std::max(mx, static_cast<int>(v));
  }
  if (mx == 0) return c.fail("`choice` must name a lineup member at least once.");
  int kk = std::max(mx, 2);
  if (c.has("k")) CK(whole(c, c["k"], "k", 2, kk));
  if (kk > kutil::maxcell) return c.fail("`k` must be at most %.0f.", kutil::maxcell);
  if (mx > kk) return c.fail("`choice` names position %d of a lineup of %d.", mx, kk);
  int sp = 1;
  CK(wholed(c, "susp", 1, 1, sp));
  if (sp > kk) return c.fail("`susp` must be a position from 1 to %d.", kk);
  double lv = 0.0;
  CK(oned(c, "level", 0.95, 0.0, 1.0, lv));
  std::vector<int> ch;
  for (double v : cv.v) ch.push_back(static_cast<int>(v) - 1);
  std::vector<double> cnt(zz(kk));
  double o[11];
  if (clab_lnsize(ch.data(), static_cast<int>(ch.size()), kk, sp - 1, qnorm(0.5 + 0.5 * lv), cnt.data(), o) != CLAB_OK) {
    return c.fail("the kernel refused the input.");
  }
  static const char *const sn[] = {"n",      "k",   "nm",    "psusp", "pbin", "fsize",
                                   "esize", "esizeb", "etr", "etrlo", "etrhi"};
  out = vobj();
  put(out, "count", vnums(cnt));
  put(out, "stat", nvec(11, sn, o));
  put(out, "susp", vnum(sp));
  return true;
}
