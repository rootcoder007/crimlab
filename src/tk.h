// SPDX-License-Identifier: AGPL-3.0-or-later
//
// Minimal checks for the standalone kernel tests (src/t_*.cpp, compiled only
// with -DCLAB_MAIN). TK(id, cond) records a check; TKN(id, got, want, tol)
// passes when |got - want| <= tol * max(1, |want|), when both are NaN, or
// when got == want (equal infinities). A failure prints the id, the observed
// and the expected value. Each test file registers its family with TK_FAMILY.
#ifndef CLAB_TK_H
#define CLAB_TK_H

#include <cmath>
#include <cstdio>

namespace tk {

struct ent {
  const char *name;
  void (*test)();
  void (*bench)();
};
inline ent *tab() {
  static ent t[32];
  return t;
}
inline int &ntab() {
  static int n = 0;
  return n;
}
inline int &nrun() {
  static int n = 0;
  return n;
}
inline int &nfail() {
  static int n = 0;
  return n;
}
struct reg {
  reg(const char *nm, void (*t)(), void (*b)()) {
    if (ntab() < 32) {
      tab()[ntab()++] = ent{nm, t, b};
    } else {
      ++nfail();
      std::printf("FAIL tk: more than 32 families\n");
    }
  }
};
inline void ok(const char *id, bool c, const char *txt) {
  ++nrun();
  if (!c) {
    ++nfail();
    std::printf("FAIL %s: %s\n", id, txt);
  }
}
inline void near(const char *id, double got, double want, double tol) {
  ++nrun();
  const double a = std::fabs(want);
  const bool c = (std::isnan(got) && std::isnan(want)) || got == want ||
                 std::fabs(got - want) <= tol * (a > 1.0 ? a : 1.0);
  if (!c) {
    ++nfail();
    std::printf("FAIL %s: got %.17g want %.17g tol %g\n", id, got, want, tol);
  }
}

}  // namespace tk

#define TK(id, cond) tk::ok((id), static_cast<bool>(cond), #cond)
#define TKN(id, got, want, tol) tk::near((id), (got), (want), (tol))
#define TK_FAMILY(name, test, bench) \
  static const tk::reg tk_reg_##test(name, test, bench)

#endif
