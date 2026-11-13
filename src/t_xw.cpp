// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for clab_bndxw and clab_synexp (src/xw.cpp) against
// sums done by hand. Part of the self-test (selftest.cpp).
#include <cmath>

#include "crimlab.h"
#include "tchk.h"

namespace {

void bndxw(tk &t) {
  // sources 0, 1, 2; 0 and 1 merge into target 0, 2 splits 0.25 / 0.75
  const double x[4] = {3, 4, 5, 8};
  const int u[4] = {0, 1, 1, 2}, from[4] = {0, 1, 2, 2}, to[4] = {0, 0, 0, 1};
  const double w[4] = {1, 1, 0.25, 0.75};
  double out[2], ws[3], bw[3];
  int best[3];
  tk_ieq(t, "bndxw.rc", clab_bndxw(x, u, 4, from, to, w, 4, 3, 2, out, ws, best, bw), CLAB_OK);
  tk_eq(t, "bndxw.t0", out[0], 3 + 9 + 2, 0);
  tk_eq(t, "bndxw.t1", out[1], 6, 0);
  tk_eq(t, "bndxw.ws", ws[2], 1, 0);
  tk_ieq(t, "bndxw.best", best[2], 1);
  tk_eq(t, "bndxw.bw", bw[2], 0.75, 0);
  // a record without a source is skipped; a source without a link
  const int u2[4] = {-1, 1, 1, 1}, from2[2] = {1, 1}, to2[2] = {0, 1};
  const double w2[2] = {0.5, 0.5};
  clab_bndxw(x, u2, 4, from2, to2, w2, 2, 3, 2, out, ws, best, bw);
  tk_eq(t, "bndxw.skip", out[0], 8.5, 0);
  tk_ieq(t, "bndxw.none", best[0], -1);
  tk_nan(t, "bndxw.nonew", bw[0]);
  tk_ieq(t, "bndxw.tie", best[1], 0);
  const int bad[1] = {3}, badu[1] = {-2};
  const double wn[1] = {-1}, xn[1] = {NAN};
  tk_ieq(t, "bndxw.n", clab_bndxw(x, u, -1, from, to, w, 4, 3, 2, out, ws, best, bw), CLAB_EARG);
  tk_ieq(t, "bndxw.m", clab_bndxw(x, u, 4, from, to, w, 0, 3, 2, out, ws, best, bw), CLAB_EARG);
  tk_ieq(t, "bndxw.ns", clab_bndxw(x, u, 4, from, to, w, 4, 0, 2, out, ws, best, bw), CLAB_EARG);
  tk_ieq(t, "bndxw.nt", clab_bndxw(x, u, 4, from, to, w, 4, 3, 0, out, ws, best, bw), CLAB_EARG);
  tk_ieq(t, "bndxw.from", clab_bndxw(x, u, 4, bad, to, w, 1, 3, 2, out, ws, best, bw), CLAB_EARG);
  tk_ieq(t, "bndxw.to", clab_bndxw(x, u, 4, from, bad, w, 1, 3, 2, out, ws, best, bw), CLAB_EARG);
  tk_ieq(t, "bndxw.w", clab_bndxw(x, u, 4, from, to, wn, 1, 3, 2, out, ws, best, bw), CLAB_EDOM);
  tk_ieq(t, "bndxw.u", clab_bndxw(x, badu, 1, from, to, w, 4, 3, 2, out, ws, best, bw), CLAB_EARG);
  tk_ieq(t, "bndxw.u3", clab_bndxw(x, bad, 1, from, to, w, 4, 3, 2, out, ws, best, bw), CLAB_EARG);
  tk_ieq(t, "bndxw.x", clab_bndxw(xn, u, 1, from, to, w, 4, 3, 2, out, ws, best, bw), CLAB_EDOM);
}

void synexp(tk &t) {
  // eta = 0 and log 3: rates 1/2 and 3/4
  const double cf[2] = {0, 1}, ax[2] = {0, std::log(3.0)}, pop[2] = {100, 40};
  double r[2], e[2];
  tk_ieq(t, "synexp.rc", clab_synexp(cf, ax, 2, 1, pop, r, e), CLAB_OK);
  tk_eq(t, "synexp.r0", r[0], 0.5, 0);
  tk_eq(t, "synexp.r1", r[1], 0.75, 1e-15);
  tk_eq(t, "synexp.e1", e[1], 30, 1e-14);
  const double cn[2] = {NAN, 1}, an[2] = {0, NAN}, pn[2] = {-1, 1};
  tk_ieq(t, "synexp.na", clab_synexp(cf, ax, 0, 1, pop, r, e), CLAB_EARG);
  tk_ieq(t, "synexp.p", clab_synexp(cf, ax, 2, -1, pop, r, e), CLAB_EARG);
  tk_ieq(t, "synexp.coef", clab_synexp(cn, ax, 2, 1, pop, r, e), CLAB_EDOM);
  tk_ieq(t, "synexp.ax", clab_synexp(cf, an, 2, 1, pop, r, e), CLAB_EDOM);
  tk_ieq(t, "synexp.pop", clab_synexp(cf, ax, 2, 1, pn, r, e), CLAB_EDOM);
}

}  // namespace

void t_xw(tk &t, int big) {
  (void)big;
  bndxw(t);
  synexp(t);
}
