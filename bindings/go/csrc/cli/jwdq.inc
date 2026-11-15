// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The weighted displacement quotient (src/bwdq.cpp).

bool wdqx(Ctx &c, Val &out) {
  Dv a, b, d;
  CK(num(c, c["a"], "a", a, 2, 0.0));
  CK(num(c, c["b"], "b", b, 2, 0.0));
  CK(num(c, c["c"], "c", d, 2, 0.0, true));
  double o[8];
  if (clab_wdqx(a.p(), b.p(), d.p(), o) != CLAB_OK) return c.fail("the kernel refused the input.");
  for (double &x : o) x = na(x);
  static const char *const nm[] = {"sm", "bdm", "wdq", "ge", "ne", "ea", "eb", "tne"};
  out = nvec(8, nm, o);
  return true;
}
