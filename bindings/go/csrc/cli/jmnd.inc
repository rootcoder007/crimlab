// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Long stays in segregation (src/bmnd.cpp): the stays that reach a
// threshold of days, alone or by the hours-and-contact rule, and their
// shares by group.

// word(): "yes" or "no" (any case, spaces around), 1 or 0; -1 else.
int yword(const std::string &s) {
  std::size_t i = 0, e = s.size();
  while (i < e && s[i] == ' ') ++i;
  while (e > i && s[e - 1] == ' ') --e;
  std::string w;
  for (std::size_t k = i; k < e; ++k) w += static_cast<char>(std::tolower(static_cast<unsigned char>(s[k])));
  if (w.find(' ') != std::string::npos) return -1;
  return w == "yes" ? 1 : w == "no" ? 0 : -1;
}
// yes(): TRUE or FALSE, 0 or 1, or "yes" or "no", one per stay.
bool yes(Ctx &c, const Val &x, const char *nm, std::size_t n, int *out) {
  if (vlen(x) != n) return c.fail("`%s` must have one value per stay, %d.", nm, static_cast<int>(n));
  const bool arr = x.t == Val::ARR || x.t == Val::OBJ;
  for (std::size_t i = 0; i < n; ++i) {
    const Val &e = arr ? x.a[i] : x;
    int v = -1;
    if (e.t == Val::BOOL) v = e.b;
    else if (e.t == Val::NUM) v = e.n != 0;
    else if (e.t == Val::STR) v = yword(e.s);
    if (v < 0) return c.fail("`%s` must be TRUE or FALSE, 0 or 1, or \"yes\" or \"no\", with no missing values.", nm);
    out[i] = v;
  }
  return true;
}

bool mndla(Ctx &c, Val &out) {
  Dv dv;
  CK(num(c, c["days"], "days", dv, -1, 0.0));
  const std::size_t sn = dv.n();
  const int n = static_cast<int>(sn);
  double th = 0.0, lv = 0.0;
  CK(binrd(c, "thr", 15, 0.0, Inf, th));
  CK(oned(c, "level", 0.95, 0.0, 1.0, lv));
  Dv hv, mv;
  const bool hh = c.has("hours");
  if (hh != c.has("miss")) return c.fail("`hours` and `miss` go together.");
  if (hh) {
    CK(num(c, c["hours"], "hours", hv, n, 0.0));
    CK(num(c, c["miss"], "miss", mv, n, 0.0));
    for (double e : mv.v) {
      if (e > 100) return c.fail("`miss` must be percentages from 0 to 100.");
    }
  }
  std::vector<int> iw(4 * sn);
  int *nal = iw.data(), *el = nal + n, *ex = el + n, *cls = ex + n;
  int am = 0;
  CK(wholed(c, "amin", 2, 1, am));
  bool hasal = false;
  if (c.has("alert")) {
    if (hh) return c.fail("`alert` goes with the duration rule, not with `hours`.");
    hasal = true;
    const Val &al = c["alert"];
    int dim[2] = {-1, -1};
    const bool m = al.t == Val::ARR && !al.a.empty() && al.a[0].t == Val::ARR && shape(al, 2, dim, 0);
    const int k = al.t == Val::OBJ ? static_cast<int>(al.a.size()) : m ? dim[1] : 1;
    for (int j = 0; j < k; ++j) {
      if (m) {
        if (dim[0] != n) return c.fail("`alert` must have one row per stay, %d.", n);
        Val s = varr();
        for (const Val &r : al.a) s.a.push_back(r.a[static_cast<std::size_t>(j)]);
        CK(yes(c, s, "alert", sn, cls));
      } else {
        CK(yes(c, al.t == Val::OBJ ? al.a[static_cast<std::size_t>(j)] : al, "alert", sn, cls));
      }
      for (int i = 0; i < n; ++i) nal[i] += cls[i];
    }
  }
  const bool he = c.has("elig"), hx = c.has("excl");
  if (he) CK(yes(c, c["elig"], "elig", sn, el));
  if (hx) CK(yes(c, c["excl"], "excl", sn, ex));
  clab_mndla(dv.p(), hh ? hv.p() : nullptr, hh ? mv.p() : nullptr, hasal ? nal : nullptr, am, he ? el : nullptr,
             hx ? ex : nullptr, n, th, cls);
  std::vector<int> gcode(sn, 0);
  int ng = 1;
  std::vector<std::string> glev;
  const bool hb = c.has("by");
  if (hb) {
    Lab g;
    CK(codes(c, c["by"], "by", sn, g));
    ng = g.g;
    glev = g.lev;
    gcode = g.code;
  }
  Lab pid;
  if (c.has("id")) CK(codes(c, c["id"], "id", sn, pid));
  const int *idc = c.has("id") ? pid.code.data() : nullptr;
  Dv wv;
  if (c.has("w")) CK(num(c, c["w"], "w", wv, n, 0.0));
  std::vector<double> ow(7 * (static_cast<std::size_t>(ng) + 1));
  std::vector<int> oi(2 * sn);
  clab_mndtab(cls, idc, gcode.data(), dv.p(), wv.n() ? wv.p() : nullptr, he ? el : nullptr, n, ng, th, ow.data(),
              oi.data());
  const int nd = idc ? 3 : 1, g0 = hb ? 0 : 1;
  static const char *const dn[] = {"stays", "people", "cumul"};
  const double a = 1 - lv;
  std::vector<std::string> gs, ds;
  std::vector<double> fn, fk, fs, fl, fh, fsol, ftor;
  for (int g = g0; g <= ng; ++g) {
    const double *o = ow.data() + 7 * (g == ng ? ng : g);
    for (int d = 0; d < nd; ++d) {
      gs.push_back(g == ng ? "all" : glev[static_cast<std::size_t>(g)]);
      ds.push_back(dn[d]);
      const double nn = d == 0 ? o[0] : o[4], kk = d == 0 ? o[1] : d == 1 ? o[5] : o[6];
      fn.push_back(nn);
      fk.push_back(kk);
      fs.push_back(nn > 0 ? kk / nn : NA);
      const bool whole = nn > 0 && nn == std::floor(nn) && kk == std::floor(kk);
      fl.push_back(!whole ? NA : kk > 0 ? clab_qbeta(a / 2, kk, nn - kk + 1, 1) : 0.0);
      fh.push_back(!whole ? NA : kk < nn ? clab_qbeta(1 - a / 2, kk + 1, nn - kk, 1) : 1.0);
      fsol.push_back(d == 0 ? o[2] : NA);
      ftor.push_back(d == 0 ? o[3] : NA);
    }
  }
  Frame f;
  f.str("group", gs);
  f.str("den", ds);
  f.num("n", fn);
  f.num("k", fk);
  f.num("share", fs);
  f.num("lo", fl);
  f.num("hi", fh);
  if (hh) {
    f.num("sol", fsol);
    f.num("tor", ftor);
  }
  out = vobj();
  put(out, "tab", f.o);
  put(out, "cls", unbox(vints(cls, sn)));
  put(out, "alerts", hasal ? unbox(vints(nal, sn)) : Val());
  put(out, "rule", vstr(hh ? "rule44" : "duration"));
  return true;
}
