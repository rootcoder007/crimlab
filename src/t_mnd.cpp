// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for the Mandela-Rules kernels (src/mnd.cpp) against
// classes and counts worked by hand. Part of the self-test (selftest.cpp).
#include "crimlab.h"
#include "tchk.h"

namespace {

void cls(tk &t) {
  const double d[4] = {5, 16, 20, 0};
  const int nal[4] = {2, 0, 0, 1}, ex[4] = {0, 1, 0, 0}, el[4] = {1, 1, 0, 1};
  int c[4];
  tk_ieq(t, "mndla.rc", clab_mndla(d, nullptr, nullptr, nullptr, 2, nullptr, nullptr, 4, 15, c),
         CLAB_OK);
  tk_ieq(t, "mndla.dur", c[0] * 1000 + c[1] * 100 + c[2] * 10 + c[3], 110);
  clab_mndla(d, nullptr, nullptr, nal, 2, nullptr, nullptr, 4, 15, c);
  tk_ieq(t, "mndla.alert", c[0] * 1000 + c[1] * 100 + c[2] * 10 + c[3], 1110);
  clab_mndla(d, nullptr, nullptr, nal, 2, nullptr, ex, 4, 15, c);
  tk_ieq(t, "mndla.excl", c[0] * 1000 + c[1] * 100 + c[2] * 10 + c[3], 1010);
  clab_mndla(d, nullptr, nullptr, nal, 2, el, ex, 4, 15, c);
  tk_ieq(t, "mndla.elig", c[0] * 1000 + c[1] * 100 + c[2] * 10 + c[3], 1000);
  // Rules 43 and 44: at most 2 hours out, every day missed
  const double dh[3] = {10, 20, 20}, h[3] = {1, 2, 3}, m[3] = {100, 100, 100},
               m9[3] = {99, 100, 100};
  clab_mndla(dh, h, m, nullptr, 2, nullptr, nullptr, 3, 15, c);
  tk_ieq(t, "mndla.rule44", c[0] * 100 + c[1] * 10 + c[2], 120);
  clab_mndla(dh, h, m9, nullptr, 2, nullptr, nullptr, 3, 15, c);
  tk_ieq(t, "mndla.miss", c[0], 0);
  const double dn[1] = {-1}, hn[1] = {-1}, mb[1] = {101};
  tk_ieq(t, "mndla.n", clab_mndla(d, nullptr, nullptr, nullptr, 2, nullptr, nullptr, -1, 15, c),
         CLAB_EARG);
  tk_ieq(t, "mndla.thr", clab_mndla(d, nullptr, nullptr, nullptr, 2, nullptr, nullptr, 4, -1, c),
         CLAB_EDOM);
  tk_ieq(t, "mndla.days", clab_mndla(dn, nullptr, nullptr, nullptr, 2, nullptr, nullptr, 1, 15, c),
         CLAB_EDOM);
  tk_ieq(t, "mndla.hours", clab_mndla(dh, hn, m, nullptr, 2, nullptr, nullptr, 1, 15, c),
         CLAB_EDOM);
  tk_ieq(t, "mndla.missb", clab_mndla(dh, h, mb, nullptr, 2, nullptr, nullptr, 1, 15, c),
         CLAB_EDOM);
}

void tab(tk &t) {
  const int c[4] = {1, 0, 1, 2}, g[4] = {0, 0, 1, 1}, id[4] = {0, 0, 1, 2};
  const double d[4] = {16, 3, 20, 10};
  double o[21];
  int wk[8];
  tk_ieq(t, "mndtab.rc", clab_mndtab(c, id, g, d, nullptr, nullptr, 4, 2, 15, o, wk), CLAB_OK);
  // per group: stays, qualifying, class 1, class 2, people, any, cumulative
  const double want[21] = {2, 1, 1, 0, 1, 1, 1, 2, 2, 1, 1, 2, 2, 1, 4, 3, 2, 1, 3, 3, 2};
  tk_veq(t, "mndtab.out", o, want, 21, 0);
  const int el[4] = {0, 1, 0, 1};
  clab_mndtab(c, id, g, d, nullptr, el, 4, 2, 15, o, wk);
  tk_eq(t, "mndtab.elig", o[20], 1, 0);
  const double w[4] = {2, 1, 1, 1};
  clab_mndtab(c, nullptr, g, d, w, nullptr, 4, 2, 15, o, wk);
  tk_true(t, "mndtab.w", o[0] == 3 && o[1] == 2 && o[4] == 0);
  const int gb[1] = {2}, cb[1] = {3};
  const double wn[1] = {-1}, dn[1] = {-1};
  tk_ieq(t, "mndtab.n", clab_mndtab(c, id, g, d, nullptr, nullptr, 0, 2, 15, o, wk), CLAB_EARG);
  tk_ieq(t, "mndtab.g", clab_mndtab(c, id, gb, d, nullptr, nullptr, 1, 2, 15, o, wk), CLAB_EARG);
  tk_ieq(t, "mndtab.c", clab_mndtab(cb, id, g, d, nullptr, nullptr, 1, 2, 15, o, wk), CLAB_EARG);
  tk_ieq(t, "mndtab.w0", clab_mndtab(c, id, g, d, wn, nullptr, 1, 2, 15, o, wk), CLAB_EDOM);
  tk_ieq(t, "mndtab.d", clab_mndtab(c, id, g, dn, nullptr, nullptr, 1, 2, 15, o, wk), CLAB_EDOM);
}

}  // namespace

void t_mnd(tk &t, int big) {
  (void)big;
  cls(t);
  tab(t);
}
