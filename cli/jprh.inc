// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Prospective hot-spotting (src/bprh.cpp).

bool prohs(Ctx &c, Val &out) {
  Pts p;
  CK(xy(c, PLANAR, p));
  const int n = p.n;
  Dv tv;
  CK(num(c, c["t"], "t", tv, n));
  double a0 = 0.0;
  CK(one(c, c["at"], "at", -Inf, Inf, a0));
  double w[4];
  CK(win(c, c["win"], w));
  double ce = 0.0;
  CK(one(c, c["cell"], "cell", 0.0, Inf, ce));
  double b = 8 * ce;
  if (c.has("bw")) CK(one(c, c["bw"], "bw", 0.0, Inf, b));
  if (b / ce > 1e4) return c.fail("`bw` must be at most 10000 cells.");
  double tws = 0.0, tu = 0.0, tp = 0.0, ah = 0.0;
  CK(oned(c, "tw", 60, 0.0, Inf, tws));
  CK(oned(c, "tunit", 7, 0.0, Inf, tu));
  CK(binrd(c, "top", 0.1, 0.0, 1.0, tp));
  if (!(tp > 0)) return c.fail("`top` must be above 0.");
  CK(binrd(c, "ahead", 7, 0.0, Inf, ah));
  const double gx = std::ceil((w[1] - w[0]) / ce), gy = std::ceil((w[3] - w[2]) / ce);
  if (gx * gy > 1e8) return c.fail("the grid would have more than 1e8 cells: make `cell` larger.");
  const int nx = sat<int>(gx), ny = sat<int>(gy);
  std::vector<double> z(zz(nx) * zz(ny));
  std::vector<int> hot(z.size());
  double ev[5];
  const int rc = clab_prohs(p.x.data(), p.y.data(), tv.p(), n, a0, w, ce, b, tws, tu, tp, ah, nx, ny, z.data(),
                            hot.data(), ev);
  CK(kchk(c, rc));
  std::vector<double> gxv, gyv;
  for (int i = 0; i < nx; ++i) gxv.push_back(w[0] + (i + 0.5) * ce);
  for (int j = 0; j < ny; ++j) gyv.push_back(w[2] + (j + 0.5) * ce);
  for (double &e : ev) e = na(e);
  static const char *const en[] = {"n", "hits", "hit", "area", "pai"};
  const double par[7] = {a0, ce, b, tws, tu, tp, ah};
  static const char *const pn[] = {"at", "cell", "bw", "tw", "tunit", "top", "ahead"};
  out = vobj();
  put(out, "x", vauto(gxv));
  put(out, "y", vauto(gyv));
  put(out, "z", vmat(z, nx, ny));
  put(out, "hot", vlmat(hot.data(), nx, ny));
  put(out, "eval", nvec(5, en, ev));
  put(out, "par", nvec(7, pn, par));
  return true;
}
