// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Repeat victimisation (src/brpv.cpp).

bool rptvic(Ctx &c, Val &out) {
  const std::size_t xn = vlen(c["id"]);
  if (xn < 1) return c.fail("`id` must hold one incident at least.");
  Lab v;
  CK(codes(c, c["id"], "id", xn, v));
  const int n = sat<int>(xn);
  Dv tv;
  const bool ht = c.has("t");
  if (ht) CK(num(c, c["t"], "t", tv, n));
  double pp = NaN;
  if (c.has("pop")) {
    CK(binr(c, c["pop"], "pop", 1.0, Inf, pp));
    if (pp < v.g) return c.fail("`pop` must be at least the number of victims, %d.", v.g);
  }
  Dv bv;
  int nb = 0;
  if (c.has("brk")) {
    if (!ht) return c.fail("`brk` needs the times `t`.");
    CK(num(c, c["brk"], "brk", bv));
    nb = sat<int>(bv.n());
    if (nb < 2) return c.fail("`brk` must hold two breaks at least.");
    for (int j = 1; j < nb; ++j) {
      if (!(bv[zz(j)] > bv[zz(j - 1)])) return c.fail("`brk` must increase.");
    }
  }
  std::vector<double> tab(3 * (xn + 1)), st(10), gap(3), lag(3 * zz(nb > 1 ? nb - 1 : 0) + 1);
  int K = 0;
  if (clab_rptvic(v.code.data(), ht ? tv.p() : nullptr, n, v.g, pp, nb ? bv.p() : nullptr, nb, tab.data(), &K,
                  st.data(), gap.data(), lag.data()) != CLAB_OK) {
    return c.fail("the kernel refused the input.");
  }
  const bool hp = !std::isnan(pp);
  const int r0 = hp ? 0 : 1, nr = K + 1 - r0;
  const std::size_t R = zz(K) + 1;
  std::vector<int> kk;
  std::vector<double> vi, in, ex;
  for (int r = 0; r < nr; ++r) {
    const std::size_t q = zz(r + r0);
    kk.push_back(r + r0);
    vi.push_back(tab[q]);
    in.push_back(tab[R + q]);
    ex.push_back(tab[2 * R + q]);
  }
  Frame f;
  f.ints("k", kk.data(), kk.size());
  f.num("victims", vi);
  f.num("incidents", in);
  if (hp) f.num("expected", ex);
  out = vobj();
  put(out, "tab", f.o);
  static const char *const sn[] = {"incidents", "victims", "repeats", "rshare", "conc",
                                   "incid", "prev", "rvict", "erep", "ratio"};
  put(out, "stat", nvec(10, sn, nav(st.data(), 10).data()));
  if (ht) {
    static const char *const gn[] = {"n", "mean", "median"};
    put(out, "gap", nvec(3, gn, nav(gap.data(), 3).data()));
  } else {
    put(out, "gap", Val());
  }
  if (nb > 1) {
    const std::size_t m = zz(nb - 1);
    Frame l;
    l.num("lo", bv.p(), m);
    l.num("hi", bv.p() + 1, m);
    l.num("n", lag.data(), m);
    l.num("share", nav(lag.data() + m, m));
    l.num("cum", nav(lag.data() + 2 * m, m));
    put(out, "lag", l.o);
  } else {
    put(out, "lag", Val());
  }
  return true;
}
