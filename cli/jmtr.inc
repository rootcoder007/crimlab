// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Bounds under monotone treatment response and selection (src/bmtr.cpp).

bool mtrbnd(Ctx &c, Val &out) {
  Dv yv;
  CK(num(c, c["y"], "y", yv));
  const std::size_t n = yv.n();
  Lab g;
  CK(codes(c, c["z"], "z", n, g));
  if (g.g < 2) return c.fail("`z` must hold two treatments at least.");
  double lo = 0.0, hi = 0.0;
  CK(oned(c, "ylo", 0, -Inf, Inf, lo));
  CK(oned(c, "yhi", 1, -Inf, Inf, hi));
  if (!(lo < hi)) return c.fail("`ylo` must be below `yhi`.");
  for (double x : yv.v) {
    if (x < lo || x > hi) return c.fail("`y` must lie between `ylo` and `yhi`.");
  }
  static const char *const an[] = {"none", "mtr", "mts", "both"};
  int as = 0;
  CK(pickd(c, "assume", an, 4, 3, as));
  Dv wv;
  const bool hw = c.has("w");
  if (hw) CK(num(c, c["w"], "w", wv, static_cast<long>(n), 0.0));
  const int k = g.g;
  std::vector<double> o(6 * zz(k));
  if (clab_mtrbnd(yv.p(), g.code.data(), hw ? wv.p() : nullptr, static_cast<int>(n), k, lo, hi, as, o.data()) !=
      CLAB_OK) {
    return c.fail("every treatment needs people with a positive weight.");
  }
  Frame f;
  f.str("z", g.lev);
  static const char *const cn[] = {"share", "mean", "lower", "upper", "atelo", "atehi"};
  for (int q = 0; q < 6; ++q) f.num(cn[q], nav(&o[zz(q) * zz(k)], zz(k)));
  const double r[2] = {lo, hi};
  static const char *const rn[] = {"ylo", "yhi"};
  out = vobj();
  put(out, "bnd", f.o);
  put(out, "assume", vstr(an[as]));
  put(out, "range", nvec(2, rn, r));
  return true;
}
