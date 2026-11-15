// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of discrimination in searches (src/bthr.cpp): hit rates, the
// threshold model and the rank order of officer groups.

// NaN as R's NA (bj::na()), for a vector.
inline std::vector<double> nav(const double *x, std::size_t n) {
  std::vector<double> v(x, x + n);
  for (double &e : v) e = na(e);
  return v;
}

bool thrtst(Ctx &c, Val &out) {
  Dv sv, kv, hv;
  CK(num(c, c["stops"], "stops", sv, -1, 0.0));
  const std::size_t n = sv.n();
  CK(num(c, c["srch"], "srch", kv, static_cast<long>(n), 0.0));
  CK(num(c, c["hits"], "hits", hv, static_cast<long>(n), 0.0));
  for (std::size_t i = 0; i < n; ++i) {
    if (kv[i] > sv[i] || hv[i] > kv[i]) {
      return c.fail("`srch` must be at most `stops`, and `hits` at most `srch`, in every row.");
    }
  }
  Lab g;
  CK(codes(c, c["grp"], "grp", n, g));
  if (g.g < 2) return c.fail("`grp` must hold at least two groups.");
  Lab a;
  int nar = 1;
  if (c.has("area")) {
    CK(codes(c, c["area"], "area", n, a));
    nar = a.g;
  } else {
    a.code.assign(n, 0);
  }
  int rf = 1;
  CK(ref(c, c["ref"], g.lev, rf));
  double lv = 0.0;
  CK(oned(c, "level", 0.95, 0.0, 1.0, lv));
  int mi = 0;
  CK(wholed(c, "maxit", 2000, 1, mi));
  static const char *const rko[] = {"pair", "joint"};
  int rk = 0;
  if (c.has("rank")) CK(pick(c, c["rank"], "rank", rko, 2, rk));
  const int ng = g.g, in = static_cast<int>(n);
  std::vector<double> tk(zz(ng), 0.0);
  for (std::size_t i = 0; i < n; ++i) tk[zz(g.code[i])] += kv[i];
  for (int j = 0; j < ng; ++j) {
    if (!(tk[zz(j)] > 0)) return c.fail("every group needs searches: %s has none.", g.lev[zz(j)].c_str());
  }
  std::vector<char> seen(zz(ng) * zz(nar), 0);
  int nc = 0;
  for (std::size_t i = 0; i < n; ++i) {
    char &s = seen[zz(a.code[i]) * zz(ng) + zz(g.code[i])];
    if (sv[i] > 0 && !s) {
      s = 1;
      ++nc;
    }
  }
  const bool fit = 2 * nc > 3 * (ng + nar - 1);
  const std::size_t un = zz(ng), cn = zz(ng) * zz(nar);
  std::vector<double> gs(6 * un), hz(4 * un), chi(3), thr(un), tse(3 * un), cfit(5 * cn), lrt(6);
  std::vector<int> conv(4, 0);
  const int rc = clab_thrtst(sv.p(), kv.p(), hv.p(), g.code.data(), a.code.data(), in, ng, nar, rf - 1, fit ? mi : 0,
                             gs.data(), hz.data(), chi.data(), thr.data(), tse.data(), cfit.data(), lrt.data(),
                             conv.data());
  if (rc == CLAB_EINTR) return kchk(c, rc);
  if (rc != CLAB_OK) return c.fail("the threshold model could not be fitted to these counts.");
  const double z = qnorm(0.5 + 0.5 * lv);
  Frame h;
  h.str("group", g.lev);
  static const char *const hn[] = {"stops", "srch", "hits", "srate", "hrate"};
  for (int q = 0; q < 5; ++q) h.num(hn[q], nav(&gs[zz(q) * un], un));
  static const char *const zn[] = {"diff", "se", "z", "p"};
  for (int q = 0; q < 4; ++q) h.num(zn[q], nav(&hz[zz(q) * un], un));
  out = vobj();
  put(out, "hit", h.o);
  static const char *const kn[] = {"chisq", "df", "p"};
  put(out, "kpt", nvec(3, kn, chi.data()));
  if (fit) {
    Frame t;
    t.str("group", g.lev);
    std::vector<double> se(un), lo(un), hi(un), df(un), ds(un), dz(un), dp(un);
    for (std::size_t j = 0; j < un; ++j) {
      se[j] = na(tse[j]);
      lo[j] = na(thr[j] - z * tse[j]);
      hi[j] = na(thr[j] + z * tse[j]);
      df[j] = tse[un + j];
      ds[j] = na(tse[2 * un + j]);
      const double zv = df[j] / tse[2 * un + j];
      const bool r0 = static_cast<int>(j) == rf - 1 || std::isnan(zv);
      dz[j] = r0 ? NA : zv;
      dp[j] = r0 ? NA : 2 * clab_pnorm(-std::fabs(zv));
    }
    t.num("thr", thr);
    t.num("se", se);
    t.num("lo", lo);
    t.num("hi", hi);
    t.num("diff", df);
    t.num("dse", ds);
    t.num("z", dz);
    t.num("p", dp);
    put(out, "thr", t.o);
    static const char *const ln[] = {"stat", "df", "p", "llfull", "llsame", "cells"};
    put(out, "lrt", nvec(6, ln, lrt.data()));
    Frame f;
    f.str("model", {"full", "same"});
    const int cv[2] = {conv[0] == 1, conv[1] == 1};
    f.lgl("conv", cv, 2);
    f.ints("iter", &conv[2], 2);
    put(out, "fit", f.o);
    std::vector<int> cg, ca;
    for (int q = 0; q < nar; ++q) {
      for (int j = 0; j < ng; ++j) {
        if (seen[zz(q) * zz(ng) + zz(j)]) {
          cg.push_back(j);
          ca.push_back(q);
        }
      }
    }
    Frame e;
    e.add("group", vstrs(g.lev, cg));
    if (c.has("area")) e.add("area", vstrs(a.lev, ca));
    static const char *const en[] = {"phi", "delta", "thr", "srate", "hrate"};
    for (int q = 0; q < 5; ++q) {
      std::vector<double> v;
      for (std::size_t k = 0; k < cg.size(); ++k) v.push_back(na(cfit[zz(q) * cn + zz(ca[k]) * zz(ng) + zz(cg[k])]));
      e.num(en[q], v);
    }
    put(out, "cells", e.o);
  } else {
    put(out, "thr", Val());
    put(out, "lrt", Val());
    put(out, "fit", Val());
    put(out, "cells", Val());
  }
  Val monov;
  if (c.has("off")) {
    Lab o;
    CK(codes(c, c["off"], "off", n, o));
    if (o.g < 2) return c.fail("`off` must hold at least two officer groups.");
    const int no = o.g, npair = no * (no - 1);
    std::vector<double> rate(5 * zz(no) * zz(ng)), res(3 * zz(npair)), mono(6 * zz(ng)), jt(2 * (1 + zz(no)));
    CK(kchk(c, clab_thrrk(sv.p(), kv.p(), hv.p(), o.code.data(), g.code.data(), in, no, ng, z, qnorm(lv),
                          rate.data(), res.data(), mono.data(), jt.data())));
    Frame mo;
    std::vector<int> mg;
    std::vector<std::string> mm;
    for (int j = 0; j < ng; ++j) {
      mg.insert(mg.end(), {j, j});
      mm.insert(mm.end(), {"search", "hit"});
    }
    mo.add("grp", vstrs(g.lev, mg));
    mo.str("meas", mm);
    static const char *const mn2[] = {"chisq", "df", "p"};
    for (int q = 0; q < 3; ++q) mo.num(mn2[q], nav(&mono[zz(q) * 2 * zz(ng)], 2 * zz(ng)));
    monov = mo.o;
    if (rk == 1) {
      Frame r;
      std::vector<int> jo, rj;
      std::vector<std::string> jm;
      std::vector<double> rr;
      for (int ms = 0; ms < 2; ++ms) {
        for (int x = 0; x < no; ++x) {
          jo.push_back(x);
          jm.push_back(ms ? "hit" : "search");
          rr.push_back(na(jt[2 * (1 + zz(x)) + zz(ms)]));
          rj.push_back(jt[zz(ms)] > 0);
        }
      }
      r.str("meas", jm);
      r.add("off", vstrs(o.lev, jo));
      r.num("rank", rr);
      r.lgl("rej", rj.data(), rj.size());
      put(out, "rank", r.o);
    } else {
    std::vector<int> a1, a2, rj;
    std::vector<std::string> me;
    for (int x = 0; x < no; ++x) {
      for (int y = x + 1; y < no; ++y) {
        for (int ms = 0; ms < 2; ++ms) {
          a1.push_back(x);
          a2.push_back(y);
          me.push_back(ms ? "hit" : "search");
        }
      }
    }
    for (int k = 0; k < npair; ++k) rj.push_back(res[zz(2 * npair + k)] > 0);
    Frame r;
    r.add("off1", vstrs(o.lev, a1));
    r.add("off2", vstrs(o.lev, a2));
    r.str("meas", me);
    r.num("zmin", nav(&res[0], zz(npair)));
    r.num("zmax", nav(&res[zz(npair)], zz(npair)));
    r.lgl("rej", rj.data(), zz(npair));
    put(out, "rank", r.o);
    }
  } else {
    put(out, "rank", Val());
  }
  put(out, "ref", vstr(g.lev[zz(rf - 1)]));
  put(out, "mono", monov);
  return true;
}
