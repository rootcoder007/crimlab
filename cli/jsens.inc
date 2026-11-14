// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Rosenbaum sensitivity bounds for matched pairs.

bool rswp(Ctx &c, Val &out) {
  Dv v, g;
  CK(num(c, c["d"], "d", v));
  Val d1;
  CK(num(c, arg(c, "gam", d1, 1), "gam", g, -1, 1.0));
  double al;
  CK(oned(c, "alpha", 0.05, 0.0, 0.5, al));
  const double z = qnorm(1.0 - al);
  const std::size_t ng = g.n();
  std::vector<double> plo(ng), pup(ng), wk(v.n());
  std::vector<int> iw(v.n());
  double gc = 0, ts = 0;
  if (clab_rswp(v.p(), static_cast<int>(v.n()), g.p(), static_cast<int>(ng), z, plo.data(), pup.data(), &gc, &ts,
                wk.data(), iw.data()) != CLAB_OK) {
    return c.fail("the kernel refused the input.");
  }
  Frame f;
  f.num("gam", g.v);
  f.num("plo", plo);
  f.num("pup", pup);
  out = vobj();
  put(out, "p", f.o);
  put(out, "gcrit", vnum(gc));
  put(out, "t", vnum(ts));
  return true;
}

bool rshl(Ctx &c, Val &out) {
  Dv v, g;
  CK(num(c, c["d"], "d", v));
  Val d1;
  CK(num(c, arg(c, "gam", d1, 1), "gam", g, -1, 1.0));
  std::vector<double> lo(g.n()), hi(g.n()), wk(v.n());
  if (clab_rshl(v.p(), static_cast<int>(v.n()), g.p(), static_cast<int>(g.n()), lo.data(), hi.data(), wk.data()) !=
      CLAB_OK) {
    return c.fail("`d` must hold at most 2^26 values.");
  }
  Frame f;
  f.num("gam", g.v);
  f.num("lo", lo);
  f.num("hi", hi);
  out = f.o;
  return true;
}
