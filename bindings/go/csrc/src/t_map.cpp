// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for map.cpp and the kernels beside it (eof.cpp,
// tps.cpp, dag.cpp): class breaks, boxes, cluster centres, symbol radii,
// an EOF, a thin-plate plane and the colliders of a small graph, each
// worked by hand, and the failure paths. Part of the self-test
// (selftest.cpp).
#include <cmath>

#include "crimlab.h"
#include "tchk.h"

namespace {

void brk(tk &tt) {
  // 1 to 5 and a missing value in two classes of equal width: breaks 1,
  // 3, 5; 1 and 2 in class 1, 3 (on a break) and 5 (the top) in class 2,
  // the missing value in none
  const double v[6] = {1, 2, 3, 4, 5, NAN};
  double br[6];
  int nb = 0, cls[6];
  tk_ieq(tt, "map.brk.rc", clab_mapbrk(v, 6, 2, 0, br, &nb, cls), CLAB_OK);
  tk_ieq(tt, "map.brk.nb", nb, 3);
  tk_near(tt, "map.brk.mid", br[1], 3, 0);
  tk_ieq(tt, "map.brk.cls", cls[0] * 10000 + cls[1] * 1000 + cls[2] * 100 + cls[4] * 10 + cls[5], 11220);
  // quantiles of 1 1 1 2 at 0, 1/2, 1: 1, 1, 2, merged to 1, 2
  const double q[4] = {1, 1, 1, 2};
  clab_mapbrk(q, 4, 2, 1, br, &nb, cls);
  tk_ieq(tt, "map.brk.qnb", nb, 2);
  tk_near(tt, "map.brk.qhi", br[1], 2, 0);
  tk_ieq(tt, "map.brk.qcls", cls[3], 1);
  // one value: one class from 0.5 below to 0.5 above
  const double c[2] = {4, 4};
  clab_mapbrk(c, 2, 3, 0, br, &nb, cls);
  tk_near(tt, "map.brk.lo", br[0], 3.5, 0);
  tk_near(tt, "map.brk.hi", br[1], 4.5, 0);
  // given breaks that do not cover the values
  br[0] = 2;
  br[1] = 3;
  nb = 2;
  tk_ieq(tt, "map.brk.cover", clab_mapbrk(v, 6, 0, 2, br, &nb, cls), CLAB_EDOM);
  const double na[1] = {NAN};
  tk_ieq(tt, "map.brk.none", clab_mapbrk(na, 1, 2, 0, br, &nb, cls), CLAB_EDOM);
  // the box of the finite points (0, 1) and (2, 1): y widened to 0.5, 1.5
  const double x[3] = {0, 2, NAN}, y[3] = {1, 1, 5};
  double bb[4];
  tk_ieq(tt, "map.bb.rc", clab_mapbb(x, y, 3, bb), CLAB_OK);
  tk_near(tt, "map.bb.x", bb[1] - bb[0], 2, 0);
  tk_near(tt, "map.bb.y0", bb[2], 0.5, 0);
  tk_near(tt, "map.bb.y1", bb[3], 1.5, 0);
}

void cen(tk &tt) {
  // clusters 1 = {(0, 0), (2, 0)} and 2 = {(4, 3)}, the last point noise;
  // a third cluster with no point
  const double x[4] = {0, 2, 4, 10}, y[4] = {0, 0, 3, 1};
  const int cl[4] = {1, 1, 2, 0};
  double cx[3], cy[3];
  int cn[3];
  tk_ieq(tt, "map.cen.rc", clab_ptcen(x, y, cl, 4, 3, cx, cy, cn), CLAB_OK);
  tk_near(tt, "map.cen.x1", cx[0], 1, 0);
  tk_near(tt, "map.cen.y2", cy[1], 3, 0);
  tk_ieq(tt, "map.cen.n", cn[0] * 100 + cn[1] * 10 + cn[2], 210);
  tk_ieq(tt, "map.cen.nan", std::isnan(cx[2]) ? 1 : 0, 1);
  const int bad[4] = {1, 4, 0, 0};
  tk_ieq(tt, "map.cen.bad", clab_ptcen(x, y, bad, 4, 3, cx, cy, cn), CLAB_EARG);
}

void sym(tk &tt) {
  // values 4 and 1 (and none) at (0, 0), (3, 0), (0, 4): the box is 3 by
  // 4, so rmax = 4 / 15; radii rmax and rmax / 2; the first drawn first
  const double x[3] = {0, 3, 0}, y[3] = {0, 0, 4}, v[3] = {4, 1, NAN};
  double bb[4], r[3], rm = 0;
  int od[3], no = 0;
  tk_ieq(tt, "map.sym.rc", clab_mapsym(x, y, v, 3, NAN, bb, r, od, &no, &rm), CLAB_OK);
  tk_near(tt, "map.sym.rmax", rm, 4.0 / 15, 1e-15);
  tk_near(tt, "map.sym.r1", r[1], 2.0 / 15, 1e-15);
  tk_ieq(tt, "map.sym.ord", no * 100 + od[0] * 10 + od[1], 201);
  tk_ieq(tt, "map.sym.none", std::isnan(r[2]) ? 1 : 0, 1);
  tk_near(tt, "map.sym.bb", bb[3], 4 + 4.0 / 15, 1e-15);
  const double neg[3] = {1, -1, 2};
  tk_ieq(tt, "map.sym.neg", clab_mapsym(x, y, neg, 3, 1, bb, r, od, &no, &rm), CLAB_EARG);
}

void eof(tk &tt) {
  // columns (1, 3) and (2, 2): centred (-1, 1) and (0, 0); A'A = diag(2,
  // 0), one component: the EOF (1, 0), the series (-1, 1), d^2 = 2
  const double z[4] = {1, 3, 2, 2};
  double ctr[2], ev[4], pc[4], d2[2], tot = 0;
  int nk = 0;
  tk_ieq(tt, "map.eof.rc", clab_eofst(z, 2, 2, 0, ctr, ev, pc, d2, &tot, &nk), CLAB_OK);
  tk_ieq(tt, "map.eof.nk", nk, 1);
  tk_near(tt, "map.eof.ctr", ctr[0], 2, 0);
  tk_near(tt, "map.eof.v", ev[0], 1, 1e-15);
  tk_near(tt, "map.eof.pc", pc[0], -1, 1e-15);
  tk_near(tt, "map.eof.d2", d2[0], 2, 1e-14);
  tk_near(tt, "map.eof.tot", tot, 2, 0);
  // a wide matrix (2 periods, 3 places), through A A': rows (0, 0, 0)
  // and (2, 4, 4) centre to -/+ (1, 2, 2); d^2 = 2 * 9 = 18 and the EOF
  // (1, 2, 2) / 3, its largest loading the first of the two 2/3s
  const double w[6] = {0, 2, 0, 4, 0, 4};
  double ew[6], pw[4];
  tk_ieq(tt, "map.eof.wrc", clab_eofst(w, 2, 3, 0, ctr, ew, pw, d2, &tot, &nk), CLAB_OK);
  tk_near(tt, "map.eof.wd2", d2[0], 18, 1e-13);
  tk_near(tt, "map.eof.wv", ew[1], 2.0 / 3, 1e-14);
  tk_near(tt, "map.eof.wpc", pw[1], 3, 1e-13);
  const double one[4] = {5, 5, 1, 1};
  tk_ieq(tt, "map.eof.flat", clab_eofst(one, 2, 2, 0, ctr, ev, pc, d2, &tot, &nk), CLAB_EDOM);
}

void tps(tk &tt) {
  // three points fix a plane, 1 + x + 2 y, with no bending (a = 0: T'a =
  // 0 has only that solution); at (1, 1) the surface is 4
  const double x[3] = {0, 1, 0}, y[3] = {0, 0, 1}, z[3] = {1, 2, 3}, gx[2] = {0, 1}, gy[2] = {0, 1};
  double a[3], b[3], fit[3], gz[4];
  tk_ieq(tt, "map.tps.rc", clab_tpsurf(x, y, z, 3, 0, gx, 2, gy, 2, a, b, fit, gz), CLAB_OK);
  tk_near(tt, "map.tps.b0", b[0], 1, 1e-12);
  tk_near(tt, "map.tps.b2", b[2], 2, 1e-12);
  tk_near(tt, "map.tps.a", a[1], 0, 1e-12);
  tk_near(tt, "map.tps.fit", fit[2], 3, 1e-12);
  tk_near(tt, "map.tps.gz", gz[3], 4, 1e-12);
  const double lx[3] = {0, 1, 2};
  tk_ieq(tt, "map.tps.line", clab_tpsurf(lx, lx, z, 3, 0, gx, 2, gy, 2, a, b, fit, gz), CLAB_EDOM);
}

void dag(tk &tt) {
  // 0 -> 2, 1 -> 2, 0 -> 1, 3 -> 2: node 2 has causes 0, 1, 3: (0, 1)
  // shielded by 0 -> 1, (0, 3) and (1, 3) not
  const int f[4] = {0, 1, 0, 3}, t[4] = {2, 2, 1, 2};
  int a[3], c[3], b[3], s[3], nc = 0, nd = 0;
  tk_ieq(tt, "map.dag.rc", clab_dagcol(f, t, 4, 4, 3, a, c, b, s, &nc, &nd), CLAB_OK);
  tk_ieq(tt, "map.dag.nc", nc * 10 + nd, 34);
  tk_ieq(tt, "map.dag.c", c[0] * 100 + c[1] * 10 + c[2], 222);
  tk_ieq(tt, "map.dag.ab", a[1] * 10 + b[1], 3);
  tk_ieq(tt, "map.dag.sh", s[0] * 100 + s[1] * 10 + s[2], 100);
  const int cf[2] = {0, 1}, ct[2] = {1, 0};
  tk_ieq(tt, "map.dag.cyc", clab_dagcol(cf, ct, 2, 2, 0, nullptr, nullptr, nullptr, nullptr, &nc, &nd), CLAB_EDOM);
  tk_ieq(tt, "map.dag.self", clab_dagcol(cf, cf, 2, 2, 0, nullptr, nullptr, nullptr, nullptr, &nc, &nd), CLAB_EDOM);
}

}  // namespace

void t_map(tk &t, int big) {
  (void)big;
  brk(t);
  cen(t);
  sym(t);
  eof(t);
  tps(t);
  dag(t);
}
