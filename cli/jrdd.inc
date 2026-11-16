// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Regression discontinuity (src/rdd.cpp) for clab_jrun(): rdloc() and
// rddens() of R/rdloc.R, with their checks in their order and words. The
// local fits are clab_rdfit, which the R function calls too. The data
// R keeps for the plots (the attribute plt) are not part of the result.

namespace rdd {

const char *const kerns[] = {"tri", "uni", "epa"};

// .rdfit(): one fit, or false with R's message.
bool fit(Ctx &c, const std::vector<double> &y, const std::vector<double> &r, const std::vector<double> *t,
         double h, int kern, int p, const std::vector<int> *cl, bool g, double *st) {
  const int rc = clab_rdfit(y.data(), r.data(), t ? t->data() : nullptr, static_cast<int>(y.size()), h, kern, p,
                            cl ? cl->data() : nullptr, g ? 1 : 0, st);
  if (rc != CLAB_OK) return ax::kerr(c, "rdloc", rc);
  if (st[6] == 1) {
    return ckstop(c, "fewer than %d rows on a side within the bandwidth %g: widen it", p + 2, h);
  }
  if (st[6] == 2) return ckstop(c, "one cluster within the bandwidth");
  return true;
}

}  // namespace rdd

bool rdloc(Ctx &c, Val &out) {
  using namespace ax;
  int bw = 0, kern = 0, ci = 0, se = 0, ssc = 0;
  CK(marg(c, "bw", {"ik", "user"}, bw));
  CK(marg(c, "kern", {"tri", "uni", "epa"}, kern));
  CK(marg(c, "ci", {"robust", "conventional"}, ci));
  CK(marg(c, "se", {"hc1", "cluster"}, se));
  CK(marg(c, "ssc", {"G", "none"}, ssc));
  Val d;
  CK(dat(c, c["data"], "data", d));
  Cols yr, tr;
  CK(acols(c, d, nms(c, {"y", "run"}), "y, run", {2}, false, yr));
  CK(acol(c, d, "trt", {1}, true, tr));
  std::vector<double> yv, rv, tv;
  CK(amat(c, d, {yr[0]}, yv));
  CK(amat(c, d, {yr[1]}, rv));
  if (!tr.empty()) CK(amat(c, d, tr, tv));
  Nv q;
  CK(cknumd(c, "cut", 0, q));
  const double cut = q[0];
  double lv = 0;
  CK(alevel(c, lv));
  bool lo = false, hi = false;
  for (double v : rv) {
    lo = lo || v < cut;
    hi = hi || v >= cut;
  }
  if (!lo || !hi) return ckstop(c, "`cut` must have rows on both sides");
  std::vector<int> cl;
  if (se == 1) {
    if (!c.has("clus")) return ckstop(c, "`se = \"cluster\"` needs `clus`");
    Cols k;
    CK(acol(c, d, "clus", {1}, false, k));
    cl = acode(*fc(d, k[0])).c;
  }
  double h = 0;
  Val ik;
  if (bw == 1) {
    if (!c.has("h")) return ckstop(c, "`bw = \"user\"` needs `h`");
    CK(cknum(c, c["h"], "h", q, true, -HUGE_VAL, HUGE_VAL, true));
    h = q[0];
  } else {
    if (c.has("h")) return ckstop(c, "give `h` with `bw = \"user\"`");
    double st[10];
    const int rc = clab_ikbw(yv.data(), rv.data(), static_cast<int>(yv.size()), cut, kern, &h, st);
    if (rc == CLAB_EDOM) {
      return c.fail("the bandwidth of Imbens and Kalyanaraman cannot be computed: too few rows near the cutoff or a "
                    "third derivative of 0; give `h` with `bw = \"user\"`.");
    }
    if (rc != CLAB_OK) return c.fail("%s", clab_strerror(rc));
    static const char *const nm[] = {"h1", "f", "vl", "vr", "m3", "h2l", "h2r", "m2l", "m2r", "ck"};
    ik = nvec(10, nm, st);
  }
  std::vector<double> rc0(rv.size());
  for (std::size_t i = 0; i < rv.size(); ++i) rc0[i] = rv[i] - cut;
  const std::vector<int> *pc = se == 1 ? &cl : nullptr;
  double s1[7], s2[7];
  CK(rdd::fit(c, yv, rc0, tr.empty() ? nullptr : &tv, h, kern, 1, pc, ssc == 0, s1));
  CK(rdd::fit(c, yv, rc0, tr.empty() ? nullptr : &tv, h, kern, 2, pc, ssc == 0, s2));
  Cols tm = {"conventional", "robust"};
  std::vector<double> b = {s1[0], s2[0]}, sv = {s1[1], s2[1]};
  if (ci == 0) {
    std::swap(tm[0], tm[1]);
    std::swap(b[0], b[1]);
    std::swap(sv[0], sv[1]);
  }
  Val tab = atab(tm, b, sv, lv, NaN, false);
  tab.k[0] = "method";
  Val pl;
  if (c.has("plac")) {
    CK(cknum(c, c["plac"], "plac", q, false));
    for (double v : q.v) {
      if (v == cut) return ckstop(c, "`plac` must not hold `cut` itself");
    }
    std::vector<double> pe, ps, pp, pn;
    for (double p0 : q.v) {
      std::vector<double> ys, rs;
      std::vector<int> cs;
      for (std::size_t i = 0; i < rv.size(); ++i) {
        if (p0 < cut ? rv[i] < cut : rv[i] >= cut) {
          ys.push_back(yv[i]);
          rs.push_back(rv[i] - p0);
          if (se == 1) cs.push_back(cl[i]);
        }
      }
      double s[7];
      CK(rdd::fit(c, ys, rs, nullptr, h, kern, 1, se == 1 ? &cs : nullptr, ssc == 0, s));
      pe.push_back(s[0]);
      ps.push_back(s[1]);
      pp.push_back(2 * clab_pnorm(-std::fabs(s[0] / s[1])));
      pn.push_back(s[2]);
    }
    Frame f;
    f.num("cut", q.v);
    f.num("est", pe);
    f.num("se", ps);
    f.num("p", pp);
    f.num("n", pn);
    pl = f.o;
  }
  out = vobj();
  put(out, "coef", tab);
  put(out, "h", vnum(h));
  put(out, "n", vnum(s1[2]));
  put(out, "nl", vnum(s1[3]));
  put(out, "nr", vnum(s1[4]));
  put(out, "ncl", vnum(s1[5]));
  put(out, "ik", ik);
  put(out, "placebo", pl);
  put(out, "cut", vnum(cut));
  put(out, "kern", vstr(rdd::kerns[kern]));
  put(out, "fuzzy", vbool(!tr.empty()));
  put(out, "se", vstr(se == 1 ? "cluster" : "hc1"));
  put(out, "level", vnum(lv));
  return true;
}

bool rddens(Ctx &c, Val &out) {
  using namespace ax;
  const Val &o = c["r"];
  std::vector<double> r = nums(o);
  bool ok = numv(o);
  for (double v : r) ok = ok && std::isfinite(v);
  if (!ok) return ckstop(c, "`r` must hold finite numbers");
  Nv q;
  CK(cknumd(c, "cut", 0, q));
  const double cut = q[0];
  double b = 0, h = 0;
  if (c.has("b")) {
    CK(cknum(c, c["b"], "b", q, true, -HUGE_VAL, HUGE_VAL, true));
    b = q[0];
  }
  if (c.has("h")) {
    CK(cknum(c, c["h"], "h", q, true, -HUGE_VAL, HUGE_VAL, true));
    h = q[0];
  }
  double lv = 0;
  CK(alevel(c, lv));
  double s[7];
  const int rc = clab_mccr(r.data(), static_cast<int>(r.size()), cut, b, h, s);
  if (rc == CLAB_EDOM) {
    return c.fail("the density cannot be estimated: the cutoff is outside the data, or a side has too few bins (give "
                  "`b` or `h`).");
  }
  if (rc != CLAB_OK) return c.fail("%s", clab_strerror(rc));
  Val tab = atab({"theta"}, {s[0]}, {s[1]}, lv, NaN, false);
  tab.k.erase(tab.k.begin());
  tab.a.erase(tab.a.begin());
  tab.k[0] = "theta";
  put(tab, "fl", vnums(&s[2], 1));
  put(tab, "fr", vnums(&s[3], 1));
  out = vobj();
  put(out, "tab", tab);
  put(out, "b", vnum(s[4]));
  put(out, "h", vnum(s[5]));
  put(out, "nbin", vnum(s[6]));
  put(out, "n", vnum(static_cast<double>(r.size())));
  put(out, "cut", vnum(cut));
  put(out, "level", vnum(lv));
  return true;
}
