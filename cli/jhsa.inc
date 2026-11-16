// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The accuracy of a hot-spot map (src/bhsa.cpp).

bool hsacc(Ctx &c, Val &out) {
  Dv ov;
  CK(num(c, c["obs"], "obs", ov, -1, 0.0));
  const std::size_t n = ov.n();
  const Val &pr = c["pred"];
  bool lg = pr.t == Val::BOOL;
  if (pr.t == Val::ARR) {
    lg = !pr.a.empty();
    for (const Val &e : pr.a) lg = lg && (e.t == Val::BOOL || e.t == Val::NUL);
  }
  std::vector<int> hot;
  Dv pv;
  if (lg) {
    if (vlen(pr) != n) return c.fail("`pred` must have one value per cell, %d.", static_cast<int>(n));
    std::vector<int> b;
    readb(pr, b);
    for (int v : b) {
      if (v < 0) return c.fail("`pred` must not be NA.");
    }
    hot = b;
  } else {
    CK(num(c, pr, "pred", pv, static_cast<std::int64_t>(n)));
  }
  Dv av;
  if (c.has("area")) {
    CK(num(c, c["area"], "area", av, static_cast<std::int64_t>(n), 0.0, true));
  } else {
    ones(av, n);
  }
  Dv tv;
  Val d1;
  CK(num(c, arg(c, "top", d1, 0.1), "top", tv, -1, 0.0, true));
  for (double x : tv.v) {
    if (x > 1) return c.fail("`top` must be at most 1.");
  }
  const int nt = lg ? 1 : static_cast<int>(tv.n());
  Dv rv;
  const bool hp = c.has("prev");
  if (hp) CK(num(c, c["prev"], "prev", rv, static_cast<std::int64_t>(n), 0.0));
  double tot = 0;
  for (double x : ov.v) tot += x;
  if (!(tot > 0)) return c.fail("`obs` has no crimes: every count is 0.");
  std::vector<double> res(11 * zz(nt));
  if (clab_hsacc(lg ? nullptr : pv.p(), lg ? hot.data() : nullptr, ov.p(), av.p(), hp ? rv.p() : nullptr,
                 static_cast<int>(n), tv.p(), nt, res.data()) != CLAB_OK) {
    return c.fail("the kernel refused the input.");
  }
  static const char *const cn[] = {"top", "cells", "area", "ashare", "n", "N", "hit", "pai", "nmax", "pei", "rri"};
  Frame f;
  for (int k = 0; k < 11; ++k) f.num(cn[k], nav(&res[zz(k) * zz(nt)], zz(nt)));
  out = f.o;
  return true;
}
