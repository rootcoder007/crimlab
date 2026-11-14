// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Courts (src/bjdg.cpp): judge designs as instruments, and the share of
// sentences at a mandatory minimum.

// dec(): 0 or 1 decisions as numbers.
bool dec01(Ctx &c, const Val &x, const char *nm, std::size_t n, std::vector<double> &d) {
  std::vector<int> b;
  CK(bin01(c, x, nm, n, b));
  d.assign(b.begin(), b.end());
  return true;
}
// One of the two columns of x: an object of two columns, or a matrix of
// two columns given as rows.
bool col2(Ctx &c, const Val &x, int j, const char *nm, Val &v) {
  if (x.t == Val::OBJ) {
    if (x.a.size() != 2) return c.fail("`%s` must have two columns.", nm);
    v = x.a[static_cast<std::size_t>(j)];
    return true;
  }
  int dim[2] = {-1, -1};
  if (x.t != Val::ARR || x.a.empty() || !shape(x, 2, dim, 0) || dim[1] != 2) return c.fail("`%s` must have two columns.", nm);
  v = varr();
  for (const Val &r : x.a) {
    if (r.a[static_cast<std::size_t>(j)].t == Val::STR) return c.fail("`%s` must hold numbers.", nm);
    v.a.push_back(r.a[static_cast<std::size_t>(j)]);
  }
  return true;
}

bool jdgiv(Ctx &c, Val &out) {
  out = vobj();
  if (!c.has("z")) {
    Val d0, d1, y0, y1;
    CK(col2(c, c["d"], 0, "d", d0));
    CK(col2(c, c["d"], 1, "d", d1));
    CK(col2(c, c["y"], 0, "y", y0));
    CK(col2(c, c["y"], 1, "y", y1));
    const std::size_t n = vlen(d0);
    std::vector<double> a, b;
    CK(dec01(c, d0, "d", n, a));
    CK(dec01(c, d1, "d", n, b));
    Dv u, v, wv;
    CK(num(c, y0, "y", u, static_cast<long>(n)));
    CK(num(c, y1, "y", v, static_cast<long>(n)));
    if (c.has("w")) CK(num(c, c["w"], "w", wv, static_cast<long>(n), 0.0));
    double o[12];
    if (clab_jdgpop(a.data(), b.data(), u.p(), v.p(), wv.n() ? wv.p() : nullptr, static_cast<int>(n), o) != CLAB_OK) {
      return c.fail("`w` must have a positive total.");
    }
    static const char *const sn[] = {"complier", "defier", "always", "never"};
    static const char *const en[] = {"itt", "fs", "wald", "ate"};
    for (int k = 4; k < 8; ++k) o[k] = na(o[k]);
    put(out, "type", vstr("pop"));
    put(out, "est", nvec(4, en, o + 8));
    put(out, "shares", nvec(4, sn, o));
    put(out, "effects", nvec(4, sn, o + 4));
    return true;
  }
  const std::size_t n = vlen(c["z"]);
  std::vector<double> zd, dd;
  CK(dec01(c, c["z"], "z", n, zd));
  CK(dec01(c, c["d"], "d", n, dd));
  Dv yv, wv, sv, ev;
  CK(num(c, c["y"], "y", yv, static_cast<long>(n)));
  if (c.has("w")) CK(num(c, c["w"], "w", wv, static_cast<long>(n), 0.0));
  Val ds = vnums(std::vector<double>{0, 0.05, 0.1}), d1;
  CK(num(c, c.has("dsh") ? c["dsh"] : ds, "dsh", sv, -1, 0.0));
  CK(num(c, arg(c, "deff", d1, 0), "deff", ev));
  const std::size_t nd = sv.n(), ne = ev.n();
  for (double e : sv.v) {
    if (e > 1) return c.fail("`dsh` must be shares from 0 to 1.");
  }
  double est[6];
  std::vector<double> sens(nd * ne);
  if (clab_jdgiv(zd.data(), dd.data(), yv.p(), wv.n() ? wv.p() : nullptr, static_cast<int>(n), sv.p(), static_cast<int>(nd),
                 ev.p(), static_cast<int>(ne), est, sens.data()) != CLAB_OK) {
    return c.fail("`z` must take both values, each with a positive weight.");
  }
  static const char *const en[] = {"fs", "itt", "wald"};
  static const char *const sn[] = {"complier", "always", "never"};
  double e3[3] = {est[0], est[1], na(est[2])};
  std::vector<double> gs, ge, gc, gt;
  for (std::size_t b = 0; b < ne; ++b) {
    for (std::size_t a = 0; a < nd; ++a) {
      gs.push_back(sv[a]);
      ge.push_back(ev[b]);
      gc.push_back(est[3] + sv[a]);
      gt.push_back(na(sens[a + nd * b]));
    }
  }
  Frame g;
  g.num("dsh", gs);
  g.num("deff", ge);
  g.num("cshare", gc);
  g.num("ceff", gt);
  put(out, "type", vstr("obs"));
  put(out, "est", nvec(3, en, e3));
  put(out, "shares", nvec(3, sn, est + 3));
  put(out, "sens", g.o);
  return true;
}

bool mmshr(Ctx &c, Val &out) {
  Dv sv, mv;
  CK(num(c, c["sent"], "sent", sv));
  const std::size_t n = sv.n();
  CK(num(c, c["min"], "min", mv));
  if (mv.n() != 1 && mv.n() != n) return c.fail("`min` must be one number or one per case, %d.", static_cast<int>(n));
  std::vector<double> mr(n);
  for (std::size_t i = 0; i < n; ++i) mr[i] = mv[mv.n() == 1 ? 0 : i];
  Lab g;
  const bool hg = c.has("grp");
  if (hg) CK(codes(c, c["grp"], "grp", n, g));
  const int ng = hg ? g.g : 0;
  std::vector<double> ov(5 * (static_cast<std::size_t>(ng) + 1));
  clab_mmshr(sv.p(), mr.data(), hg ? g.code.data() : nullptr, static_cast<int>(n), ng, ov.data());
  const int r0 = hg ? 0 : ng, nr = hg ? ng + 1 : 1;
  std::vector<std::string> gl;
  for (int k = 0; k < nr; ++k) gl.push_back(k < nr - 1 ? g.lev[static_cast<std::size_t>(k)] : "all");
  Frame f;
  f.str("group", gl);
  static const char *const cn[] = {"n", "at", "below", "above", "mexcess"};
  for (int q = 0; q < 5; ++q) f.num(cn[q], ov.data() + q * (ng + 1) + r0, static_cast<std::size_t>(nr));
  out = f.o;
  return true;
}
