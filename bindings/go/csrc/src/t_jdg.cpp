// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for the courts and sentencing kernels (src/jdg.cpp)
// against values worked by hand. Part of the self-test (selftest.cpp).
#include <cmath>
#include <limits>

#include "crimlab.h"
#include "tchk.h"

namespace {

const double inf = std::numeric_limits<double>::infinity();

void iv(tk &t) {
  const double z[8] = {1, 1, 1, 1, 0, 0, 0, 0}, d[8] = {1, 1, 1, 0, 1, 0, 0, 0};
  const double y[8] = {5, 4, 3, 1, 2, 1, 0, 1}, ds[2] = {0, 0.1}, de[2] = {0, 2};
  double est[6], sens[4];
  tk_ieq(t, "jdgiv.rc", clab_jdgiv(z, d, y, nullptr, 8, ds, 2, de, 2, est, sens), CLAB_OK);
  tk_eq(t, "jdgiv.fs", est[0], 0.5, 1e-15);
  tk_eq(t, "jdgiv.itt", est[1], 2.25, 1e-15);
  tk_eq(t, "jdgiv.wald", est[2], 4.5, 1e-15);
  tk_eq(t, "jdgiv.c", est[3], 0.5, 1e-15);
  tk_eq(t, "jdgiv.a", est[4], 0.25, 1e-15);
  tk_eq(t, "jdgiv.s0", sens[0], 4.5, 1e-15);
  tk_eq(t, "jdgiv.s1", sens[1], 2.25 / 0.6, 1e-14);
  tk_eq(t, "jdgiv.s3", sens[3], 2.45 / 0.6, 1e-14);
  const double w[8] = {2, 1, 1, 1, 1, 1, 1, 1};
  clab_jdgiv(z, d, y, w, 8, ds, 2, de, 2, est, sens);
  tk_eq(t, "jdgiv.wfs", est[0], 0.55, 1e-15);
  // no first stage: no ratio
  const double d1[8] = {1, 1, 1, 1, 1, 1, 1, 1};
  clab_jdgiv(z, d1, y, nullptr, 8, ds, 2, de, 2, est, sens);
  tk_true(t, "jdgiv.nofs", std::isnan(est[2]) && std::isnan(sens[0]));
  const double z2[2] = {2, 0}, z1[2] = {1, 1}, bad[1] = {1.5}, binf[1] = {inf};
  const double wn[8] = {-1, 1, 1, 1, 1, 1, 1, 1};
  tk_ieq(t, "jdgiv.n", clab_jdgiv(z, d, y, nullptr, 1, ds, 2, de, 2, est, sens), CLAB_EARG);
  tk_ieq(t, "jdgiv.z", clab_jdgiv(z2, d, y, nullptr, 2, ds, 2, de, 2, est, sens), CLAB_EDOM);
  tk_ieq(t, "jdgiv.arm", clab_jdgiv(z1, d, y, nullptr, 2, ds, 2, de, 2, est, sens), CLAB_EDOM);
  tk_ieq(t, "jdgiv.dsh", clab_jdgiv(z, d, y, nullptr, 8, bad, 1, de, 2, est, sens), CLAB_EDOM);
  tk_ieq(t, "jdgiv.deff", clab_jdgiv(z, d, y, nullptr, 8, ds, 2, binf, 1, est, sens), CLAB_EDOM);
  tk_ieq(t, "jdgiv.w", clab_jdgiv(z, d, y, wn, 8, ds, 2, de, 2, est, sens), CLAB_EDOM);
  // a known population: complier, defier, always, never
  const double a0[4] = {0, 1, 1, 0}, a1[4] = {1, 0, 1, 0}, y0[4] = {0, 0, 1, 2},
               y1[4] = {3, 1, 1, 2};
  double o[12];
  tk_ieq(t, "jdgpop.rc", clab_jdgpop(a0, a1, y0, y1, nullptr, 4, o), CLAB_OK);
  tk_eq(t, "jdgpop.share", o[1], 0.25, 0);
  tk_eq(t, "jdgpop.eff", o[4], 3, 0);
  tk_eq(t, "jdgpop.itt", o[8], 0.5, 1e-15);
  tk_eq(t, "jdgpop.fs", o[9], 0, 0);
  tk_nan(t, "jdgpop.wald", o[10]);
  tk_eq(t, "jdgpop.ate", o[11], 1, 1e-15);
  // no defiers: the Wald ratio is the compliers' effect
  const double b0[4] = {0, 0, 1, 0}, b1[4] = {1, 1, 1, 0};
  clab_jdgpop(b0, b1, y0, y1, nullptr, 4, o);
  tk_eq(t, "jdgpop.late", o[10], o[4], 1e-15);
  tk_eq(t, "jdgpop.late2", o[10], 2, 1e-15);
  tk_nan(t, "jdgpop.nodef", o[5]);
  const double w0[4] = {0, 0, 0, 0}, bb[4] = {2, 0, 0, 0};
  tk_ieq(t, "jdgpop.n", clab_jdgpop(a0, a1, y0, y1, nullptr, 0, o), CLAB_EARG);
  tk_ieq(t, "jdgpop.d", clab_jdgpop(bb, a1, y0, y1, nullptr, 4, o), CLAB_EDOM);
  tk_ieq(t, "jdgpop.w", clab_jdgpop(a0, a1, y0, y1, w0, 4, o), CLAB_EDOM);
}

void mono(tk &t) {
  const int jc[6] = {0, 0, 0, 1, 1, 2};
  const double d[6] = {1, 0, 1, 0, 0, 1}, y[6] = {1, 2, 3, 4, 5, 6};
  double cnt[3], one[3], rate[3], outc[3], st[7];
  tk_ieq(t, "jrate.rc", clab_jrate(jc, d, y, nullptr, 6, 3, cnt, one, rate, outc, st), CLAB_OK);
  tk_eq(t, "jrate.rate", rate[0], 2.0 / 3, 1e-15);
  tk_eq(t, "jrate.outc", outc[1], 4.5, 0);
  tk_eq(t, "jrate.one", one[0], 2, 0);
  tk_eq(t, "jrate.x2", st[0], 1.0 / 3 + 2 + 1, 1e-14);
  tk_eq(t, "jrate.df", st[1], 2, 0);
  tk_eq(t, "jrate.range", st[3] - st[2], 1, 0);
  clab_jrate(jc, d, nullptr, nullptr, 6, 3, cnt, one, rate, outc, st);
  tk_nan(t, "jrate.noy", outc[0]);
  const double d1[6] = {1, 1, 1, 1, 1, 1}, dd[1] = {2}, w[6] = {1, 1, 1, 1, 1, 0};
  const int jb[1] = {5};
  clab_jrate(jc, d1, nullptr, nullptr, 6, 3, cnt, one, rate, outc, st);
  tk_nan(t, "jrate.same", st[0]);
  tk_ieq(t, "jrate.n", clab_jrate(jc, d, y, nullptr, 0, 3, cnt, one, rate, outc, st), CLAB_EARG);
  tk_ieq(t, "jrate.jc", clab_jrate(jb, d, y, nullptr, 1, 3, cnt, one, rate, outc, st), CLAB_EARG);
  tk_ieq(t, "jrate.d", clab_jrate(jc, dd, y, nullptr, 1, 3, cnt, one, rate, outc, st), CLAB_EDOM);
  tk_ieq(t, "jrate.w", clab_jrate(jc, d, y, w, 6, 3, cnt, one, rate, outc, st), CLAB_EDOM);
  // the slope test: three judges, outcomes in [0, 2]
  const double r[3] = {0, 0.5, 1}, oc[3] = {0, 1, 1.5}, ov[3] = {0, 2, 1.5}, rz[2] = {0, 0};
  double dp[3], dy[3], bn[3], lt[3];
  int vi[3];
  tk_ieq(t, "jslope.rc", clab_jslope(r, oc, 3, 0, 2, dp, dy, bn, vi, lt), CLAB_OK);
  tk_eq(t, "jslope.dp", dp[1], 1, 0);
  tk_eq(t, "jslope.bnd", bn[2], 1, 0);
  tk_eq(t, "jslope.late", lt[0], 2, 0);
  tk_ieq(t, "jslope.ok", vi[0] + vi[1] + vi[2], 0);
  clab_jslope(r, ov, 3, 0, 2, dp, dy, bn, vi, lt);
  tk_ieq(t, "jslope.viol", vi[0], 1);
  clab_jslope(rz, oc, 2, 0, 2, dp, dy, bn, vi, lt);
  tk_nan(t, "jslope.nan", lt[0]);
  tk_ieq(t, "jslope.n", clab_jslope(r, oc, 0, 0, 2, dp, dy, bn, vi, lt), CLAB_EARG);
  tk_ieq(t, "jslope.hl", clab_jslope(r, oc, 3, 2, 0, dp, dy, bn, vi, lt), CLAB_EDOM);
}

void sub(tk &t) {
  const int jc[8] = {0, 0, 0, 0, 1, 1, 1, 1}, gc[8] = {0, 0, 1, 1, 0, 0, 1, 1};
  const double d[8] = {1, 1, 1, 0, 1, 0, 0, 0}, jt[2] = {0.75, 0.25};
  double cn[4], cr[4], pv[2], sl[2], sc[2];
  int pg[2], ps[2], pl[2], sj[2], np = -1;
  tk_ieq(t, "jsub.count",
         clab_jsub(jc, gc, d, 8, 2, 2, 1, jt, cn, cr, nullptr, nullptr, nullptr, nullptr, sl, sj,
                   sc, &np),
         CLAB_OK);
  tk_ieq(t, "jsub.np", np, 2);
  clab_jsub(jc, gc, d, 8, 2, 2, 1, jt, cn, cr, pg, ps, pl, pv, sl, sj, sc, &np);
  tk_ieq(t, "jsub.strict", ps[0], 0);
  tk_ieq(t, "jsub.group", pg[1], 1);
  tk_eq(t, "jsub.z", pv[0], -0.5 / std::sqrt(0.75 * 0.25), 1e-14);
  tk_eq(t, "jsub.z1", pv[1], -0.5 / std::sqrt(0.25 * 0.75), 1e-14);
  tk_eq(t, "jsub.slope", sl[0], 1, 1e-14);
  tk_eq(t, "jsub.cases", sc[1], 4, 0);
  tk_eq(t, "jsub.rate", cr[1], 0.5, 0);
  // too few cases: no judges, no pairs, no slope
  clab_jsub(jc, gc, d, 8, 2, 2, 3, jt, cn, cr, pg, ps, pl, pv, sl, sj, sc, &np);
  tk_true(t, "jsub.minn", np == 0 && sj[0] == 0 && std::isnan(sl[0]));
  // equal overall rates: no pair; equal group rates: z = 0
  const double je[2] = {0.5, 0.5}, d1[8] = {1, 1, 0, 0, 1, 1, 0, 0};
  clab_jsub(jc, gc, d, 8, 2, 2, 1, je, cn, cr, pg, ps, pl, pv, sl, sj, sc, &np);
  tk_true(t, "jsub.tie", np == 0 && std::isnan(sl[0]));
  clab_jsub(jc, gc, d1, 8, 2, 2, 1, jt, cn, cr, pg, ps, pl, pv, sl, sj, sc, &np);
  tk_eq(t, "jsub.se0", pv[0], 0, 0);
  const int jb[1] = {4};
  const double db[1] = {3};
  tk_ieq(t, "jsub.n", clab_jsub(jc, gc, d, 0, 2, 2, 1, jt, cn, cr, pg, ps, pl, pv, sl, sj, sc, &np),
         CLAB_EARG);
  tk_ieq(t, "jsub.jc",
         clab_jsub(jb, gc, d, 1, 2, 2, 1, jt, cn, cr, pg, ps, pl, pv, sl, sj, sc, &np), CLAB_EARG);
  tk_ieq(t, "jsub.d",
         clab_jsub(jc, gc, db, 1, 2, 2, 1, jt, cn, cr, pg, ps, pl, pv, sl, sj, sc, &np), CLAB_EDOM);
}

void mm(tk &t) {
  const double s[4] = {10, 12, 8, 10}, m[4] = {10, 10, 10, 10};
  const int g[4] = {0, 0, 1, 1};
  double o[15];
  tk_ieq(t, "mmshr.rc", clab_mmshr(s, m, g, 4, 2, o), CLAB_OK);
  // n, at, below, above, mean excess; rows: group 0, group 1, overall
  const double want[15] = {2, 2, 4, 0.5, 0.5, 0.5, 0, 0.5, 0.25, 0.5, 0, 0.25, 2, 0, 2};
  tk_veq(t, "mmshr.out", o, want, 15, 1e-15);
  double p[5];
  clab_mmshr(s, m, nullptr, 4, 0, p);
  const double wp[5] = {4, 0.5, 0.25, 0.25, 2};
  tk_veq(t, "mmshr.all", p, wp, 5, 1e-15);
  const double si[1] = {inf};
  const int gb[1] = {3};
  tk_ieq(t, "mmshr.n", clab_mmshr(s, m, g, 0, 2, o), CLAB_EARG);
  tk_ieq(t, "mmshr.inf", clab_mmshr(si, m, nullptr, 1, 0, o), CLAB_EDOM);
  tk_ieq(t, "mmshr.g", clab_mmshr(s, m, gb, 1, 2, o), CLAB_EARG);
}

}  // namespace

void t_jdg(tk &t, int big) {
  (void)big;
  iv(t);
  mono(t);
  sub(t);
  mm(t);
}
