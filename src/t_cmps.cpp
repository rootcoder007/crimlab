// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Self-test of clab_cmpsc (src/cmps.cpp): a signature against a shifted
// copy of itself finds the shift and counts every segment that can reach
// it; against an unrelated signature few segments agree; refusals.
#include <climits>
#include <cmath>
#include <cstdint>

#include "crimlab.h"
#include "tchk.h"

namespace {

double a[800], b[800], c[800];

// a walk with noise from MINSTD seed s, n points
void walk(double *x, int n, std::int64_t s) {
  double w = 0;
  for (int i = 0; i < n; ++i) {
    s = (48271 * s) % 2147483647;
    w += static_cast<double>(s) / 2147483647.0 - 0.5;
    x[i] = w + std::sin(i / 9.0);
  }
}

void run(tk &t) {
  walk(a, 800, 7);
  walk(c, 760, 99);
  for (int i = 0; i < 760; ++i) b[i] = a[i + 20];
  const int npk[3] = {5, 3, 1};
  int off = 0, ns = 0, sc = 0, cm[16];
  double pos = 0, ccp[16];
  tk_ieq(t, "cmps1", clab_cmpsc(a, 600, b, 760, 50, npk, 3, 25, 0, &off, &ns, &sc, &pos, ccp, cm), CLAB_OK);
  tk_ieq(t, "cmps2", ns, 12);
  // a[t] pairs with b[t - 20]
  tk_eq(t, "cmps3", pos, -20, 0);
  // segment 0 cannot reach lag -20, where fewer of its points than the
  // overlap needed (45) fall inside b; the others can
  tk_ieq(t, "cmps4", sc, 11);
  tk_eq(t, "cmps5", ccp[5], -20, 0);
  tk_ieq(t, "cmps6", cm[0], 0);
  // one level: the top five peaks of each segment
  const int one[1] = {5};
  clab_cmpsc(a, 600, b, 760, 50, one, 1, 25, 0, &off, &ns, &sc, &pos, ccp, cm);
  tk_eq(t, "cmps7", pos, -20, 0);
  tk_ge(t, "cmps8", sc, 11);
  // an unrelated signature
  clab_cmpsc(a, 600, c, 760, 50, npk, 3, 25, 0, &off, &ns, &sc, &pos, ccp, cm);
  tk_le(t, "cmps9", sc, 4);
  // missing values in the comparison leave a segment without candidates
  for (int i = 100; i < 150; ++i) a[i] = NAN;
  clab_cmpsc(a, 600, b, 760, 50, npk, 3, 25, 0, &off, &ns, &sc, &pos, ccp, cm);
  tk_nan(t, "cmps10", ccp[2]);
  tk_eq(t, "cmps11", pos, -20, 0);
  // refusals
  tk_ieq(t, "cmps12", clab_cmpsc(a, 40, b, 760, 50, npk, 3, 25, 0, &off, &ns, &sc, &pos, ccp, cm), CLAB_EARG);
  b[3] = HUGE_VAL;
  tk_ieq(t, "cmps13", clab_cmpsc(a, 600, b, 760, 50, npk, 3, 25, 0, &off, &ns, &sc, &pos, ccp, cm), CLAB_EDOM);
}

}  // namespace

void t_cmps(tk &t, int) { run(t); }
