// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The weighted displacement quotient (src/bwdq.cpp).

bool wdqx(Ctx &c, Val &out) {
  Dv a, b, d;
  CK(num(c, c["a"], "a", a, 2, 0.0));
  CK(num(c, c["b"], "b", b, 2, 0.0));
  CK(num(c, c["c"], "c", d, 2, 0.0, true));
  double o[6];
  if (clab_wdqx(a.p(), b.p(), d.p(), o) != CLAB_OK) return c.fail("the kernel refused the input.");
  for (double &x : o) x = na(x);
  static const char *const nm[] = {"sm", "bdm", "wdq", "ea", "eb", "net"};
  out = nvec(6, nm, o);
  return true;
}
