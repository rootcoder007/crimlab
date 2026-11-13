// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for the stock, flow and court kernels (src/flw.cpp)
// against values worked by hand. Part of the self-test (selftest.cpp).
#include <cmath>
#include <limits>

#include "crimlab.h"
#include "tchk.h"

namespace {

const double nan = std::numeric_limits<double>::quiet_NaN();
const double inf = std::numeric_limits<double>::infinity();

void stock(tk &t) {
  const double pe[2] = {100, 120}, dy[2] = {3650, 4380}, tt[2] = {365, 365};
  const double ex[2] = {1000, 1000}, rel[2] = {50, 50};
  double o[36];
  tk_ieq(t, "stkflw.rc", clab_stkflw(pe, dy, nullptr, nullptr, tt, 2, 0, ex, 100, rel, 182.5, o),
         CLAB_OK);
  tk_eq(t, "stkflw.adp", o[4], 10, 0);
  tk_eq(t, "stkflw.alos", o[7], 36.5, 0);
  tk_eq(t, "stkflw.peochg", o[9], 20, 1e-14);
  tk_eq(t, "stkflw.aloschg", o[11], 0, 0);
  tk_eq(t, "stkflw.dayschg", o[13], 20, 1e-14);
  tk_eq(t, "stkflw.pshare", o[17], 1, 1e-14);
  tk_nan(t, "stkflw.pshare0", o[16]);
  tk_eq(t, "stkflw.flowrate", o[19], 12, 1e-14);
  tk_eq(t, "stkflw.stkrate", o[20], 1, 0);
  tk_eq(t, "stkflw.frchg", o[23], 20, 1e-13);
  tk_eq(t, "stkflw.cut", o[26], 25, 0);
  tk_eq(t, "stkflw.cutpct", o[29], 100 * 25 / 12.0, 1e-14);
  tk_eq(t, "stkflw.mean", o[30], 110, 0);
  // the other pairs
  const double ad[1] = {10}, al[1] = {36.5}, p1[1] = {100}, d1[1] = {3650}, t1[1] = {365};
  clab_stkflw(nullptr, nullptr, ad, al, t1, 1, 0, nullptr, 1, nullptr, nan, o);
  tk_eq(t, "stkflw.admits", o[0], 100, 1e-14);
  clab_stkflw(p1, nullptr, nullptr, al, t1, 1, 0, nullptr, 1, nullptr, nan, o);
  tk_eq(t, "stkflw.pal", o[1], 3650, 0);
  clab_stkflw(p1, nullptr, ad, nullptr, t1, 1, 0, nullptr, 1, nullptr, nan, o);
  tk_eq(t, "stkflw.pad", o[3], 36.5, 0);
  clab_stkflw(nullptr, d1, nullptr, al, t1, 1, 0, nullptr, 1, nullptr, nan, o);
  tk_eq(t, "stkflw.dal", o[0], 100, 0);
  tk_true(t, "stkflw.norate", std::isnan(o[11]) && std::isnan(o[14]));
  clab_stkflw(pe, dy, nullptr, nullptr, tt, 2, 1, nullptr, 1, nullptr, nan, o);
  tk_nan(t, "stkflw.prev0", o[8]);
  tk_eq(t, "stkflw.prev", o[9], 20, 1e-14);
  const double z[1] = {0};
  tk_ieq(t, "stkflw.n",
         clab_stkflw(pe, dy, nullptr, nullptr, tt, 0, 0, nullptr, 1, nullptr, nan, o), CLAB_EARG);
  tk_ieq(t, "stkflw.base",
         clab_stkflw(pe, dy, nullptr, nullptr, tt, 2, 2, nullptr, 1, nullptr, nan, o), CLAB_EARG);
  tk_ieq(t, "stkflw.dadp",
         clab_stkflw(nullptr, d1, ad, nullptr, t1, 1, 0, nullptr, 1, nullptr, nan, o), CLAB_EARG);
  tk_ieq(t, "stkflw.three", clab_stkflw(p1, d1, nullptr, al, t1, 1, 0, nullptr, 1, nullptr, nan, o),
         CLAB_EARG);
  tk_ieq(t, "stkflw.per",
         clab_stkflw(pe, dy, nullptr, nullptr, tt, 2, 0, nullptr, 0, nullptr, nan, o), CLAB_EDOM);
  tk_ieq(t, "stkflw.early", clab_stkflw(pe, dy, nullptr, nullptr, tt, 2, 0, nullptr, 1, rel, -1, o),
         CLAB_EDOM);
  tk_ieq(t, "stkflw.t", clab_stkflw(p1, d1, nullptr, nullptr, z, 1, 0, nullptr, 1, nullptr, nan, o),
         CLAB_EDOM);
}

void stay(tk &t) {
  const double x[5] = {10, 2, 3, 4, 1}, br[4] = {0, 2, 5, inf}, th[1] = {3};
  const int fl[5] = {1, 1, 0, 0, 1};
  double st[9], bt[6], tt[1], fit[6];
  int o[6];
  tk_ieq(t, "lostay.rc", clab_lostay(x, nullptr, 5, br, 4, th, 1, fl, st, bt, tt, fit, o), CLAB_OK);
  tk_eq(t, "lostay.mean", st[2], 4, 0);
  tk_eq(t, "lostay.sd", st[3], std::sqrt(12.5), 1e-15);
  tk_eq(t, "lostay.median", st[4], 3, 0);
  tk_eq(t, "lostay.q1", st[5], 2, 0);
  tk_eq(t, "lostay.max", st[7], 10, 0);
  const double bw[6] = {2, 2, 1, 2, 0, 1};
  tk_veq(t, "lostay.bands", bt, bw, 6, 0);
  tk_eq(t, "lostay.thr", tt[0], 0.6, 1e-15);
  // the exponential's AIC at rate 1 / 4, the lognormal's and Pareto's
  tk_eq(t, "lostay.exp", fit[1], 2 - 2 * (5 * std::log(0.25) - 5), 1e-14);
  double sl = 0, s2 = 0;
  for (double v : x) sl += std::log(v);
  for (double v : x) s2 += (std::log(v) - sl / 5) * (std::log(v) - sl / 5) / 5;
  const double pi = 3.14159265358979323846;
  tk_eq(t, "lostay.lnorm", fit[0], 4 - 2 * (-sl - 2.5 * std::log(2 * pi * s2) - 2.5), 1e-14);
  tk_eq(t, "lostay.pareto", fit[5], 5 / sl, 1e-15);
  tk_eq(t, "lostay.parAIC", fit[2], 4 - 2 * (5 * std::log(5 / sl) - (5 / sl + 1) * sl), 1e-14);
  // frequencies: 10 counted twice
  const double w[5] = {2, 1, 1, 1, 1};
  clab_lostay(x, w, 5, br, 0, th, 0, nullptr, st, bt, tt, fit, o);
  tk_eq(t, "lostay.wmed", st[4], 3.5, 0);
  tk_eq(t, "lostay.wq1", st[5], 2.25, 1e-15);
  tk_eq(t, "lostay.wn", st[0], 6, 0);
  // one positive stay: no fit; equal stays: no lognormal or Pareto
  const double x1[2] = {0, 5}, xe[3] = {4, 4, 4};
  clab_lostay(x1, nullptr, 2, br, 0, th, 0, nullptr, st, bt, tt, fit, o);
  tk_nan(t, "lostay.nofit", fit[1]);
  clab_lostay(xe, nullptr, 3, br, 0, th, 0, nullptr, st, bt, tt, fit, o);
  tk_true(t, "lostay.flat", std::isnan(fit[0]) && std::isnan(fit[2]) && !std::isnan(fit[1]));
  clab_lostay(x, nullptr, 1, br, 0, th, 0, nullptr, st, bt, tt, fit, o);
  tk_nan(t, "lostay.sd1", st[3]);
  const double xn[1] = {-1}, wf[1] = {1.5}, w0[1] = {0}, bb[2] = {3, 1};
  tk_ieq(t, "lostay.n", clab_lostay(x, nullptr, 0, br, 0, th, 0, nullptr, st, bt, tt, fit, o),
         CLAB_EARG);
  tk_ieq(t, "lostay.x", clab_lostay(xn, nullptr, 1, br, 0, th, 0, nullptr, st, bt, tt, fit, o),
         CLAB_EDOM);
  tk_ieq(t, "lostay.w", clab_lostay(x, wf, 1, br, 0, th, 0, nullptr, st, bt, tt, fit, o),
         CLAB_EDOM);
  tk_ieq(t, "lostay.w0", clab_lostay(x, w0, 1, br, 0, th, 0, nullptr, st, bt, tt, fit, o),
         CLAB_EDOM);
  tk_ieq(t, "lostay.br", clab_lostay(x, nullptr, 1, bb, 2, th, 0, nullptr, st, bt, tt, fit, o),
         CLAB_EDOM);
}

void docket(tk &t) {
  const double a[5] = {0, 2, 3, 7, 9}, d[5] = {5, 4, 12, nan, 15};
  double o[7];
  tk_ieq(t, "bklog.rc", clab_bklog(a, d, 5, 20, o), CLAB_OK);
  tk_eq(t, "bklog.backlog", o[4], 1.75, 1e-15);
  tk_eq(t, "bklog.wait", o[3], 5.5, 0);
  tk_eq(t, "bklog.rate", o[2], 0.25, 0);
  tk_nan(t, "bklog.resid", o[6]);
  const double a2[2] = {0, 1}, d2[2] = {2, 4}, d3[2] = {2, 25};
  clab_bklog(a2, d2, 2, 4, o);
  tk_eq(t, "bklog.little", o[6], 0, 1e-15);
  clab_bklog(a2, d3, 2, 20, o);
  tk_eq(t, "bklog.after", o[1], 1, 0);
  const double ab[1] = {30}, db[1] = {1}, a3[1] = {5};
  tk_ieq(t, "bklog.n", clab_bklog(a, d, 0, 20, o), CLAB_EARG);
  tk_ieq(t, "bklog.h", clab_bklog(a, d, 5, 0, o), CLAB_EDOM);
  tk_ieq(t, "bklog.a", clab_bklog(ab, db, 1, 20, o), CLAB_EDOM);
  tk_ieq(t, "bklog.d", clab_bklog(a3, db, 1, 20, o), CLAB_EDOM);
  // case flow
  const double rc[3] = {10, 20, 0}, ds[3] = {5, 25, 5}, pd[3] = {10, 5, 0};
  double f[24];
  tk_ieq(t, "dktflw.rc", clab_dktflw(rc, ds, pd, 3, 2, f), CLAB_OK);
  tk_eq(t, "dktflw.clear", f[1], 1.25, 0);
  tk_nan(t, "dktflw.noclear", f[2]);
  tk_eq(t, "dktflw.net", f[4], -5, 0);
  tk_nan(t, "dktflw.roll0", f[6]);
  tk_eq(t, "dktflw.rr", f[8], 20, 0);
  tk_eq(t, "dktflw.rclear", f[14], 1.5, 0);
  tk_eq(t, "dktflw.pchg", f[16], -5, 0);
  tk_eq(t, "dktflw.toclear0", f[18], 2, 0);
  tk_eq(t, "dktflw.toclear", f[19], 5.0 / 15, 1e-15);
  tk_eq(t, "dktflw.cumnet", f[23], -5, 0);
  clab_dktflw(rc, ds, nullptr, 3, 2, f);
  tk_true(t, "dktflw.nopend", std::isnan(f[16]) && std::isnan(f[19]));
  const double rn[1] = {-1};
  tk_ieq(t, "dktflw.n", clab_dktflw(rc, ds, pd, 0, 2, f), CLAB_EARG);
  tk_ieq(t, "dktflw.win", clab_dktflw(rc, ds, pd, 3, 0, f), CLAB_EARG);
  tk_ieq(t, "dktflw.neg", clab_dktflw(rn, ds, pd, 1, 1, f), CLAB_EDOM);
}

void mix(tk &t) {
  // two periods, two outcomes: 10 and 30, then 20 and 20
  const double x[4] = {10, 20, 30, 20}, sc[2] = {1, 2};
  double sh[4], lo[4], hi[4], st[3], tr[4];
  const double a = 3.841458820694124;
  tk_ieq(t, "dkmix.rc", clab_dkmix(x, 2, 2, sc, a, sh, lo, hi, st, tr), CLAB_OK);
  tk_eq(t, "dkmix.share", sh[0], 0.25, 0);
  const double hf = std::sqrt(a * (a + 4 * 10 * 30 / 40.0));
  tk_eq(t, "dkmix.lo", lo[0], (a + 20 - hf) / (2 * (40 + a)), 1e-15);
  tk_eq(t, "dkmix.chisq", st[0], 16.0 / 3, 1e-14);
  tk_eq(t, "dkmix.df", st[1], 1, 0);
  tk_eq(t, "dkmix.v", st[2], std::sqrt(16.0 / 3 / 80), 1e-14);
  tk_nan(t, "dkmix.notrend", tr[0]);
  tk_eq(t, "dkmix.change", tr[2], 0.25, 0);
  // three periods: the trend statistic; an outcome never seen drops out
  const double y[9] = {10, 20, 30, 30, 20, 10, 0, 0, 0}, s3[3] = {1, 2, 3};
  double sh3[9], lo3[9], hi3[9], tr3[6];
  clab_dkmix(y, 3, 3, s3, a, sh3, lo3, hi3, st, tr3);
  // sum x s - p sum n s over p (1 - p)(sum n s^2 - (sum n s)^2 / sum n)
  const double num = (10 + 40 + 90) - 0.5 * 240, den = 0.25 * (40 * 14 - 240.0 * 240 / 120);
  tk_eq(t, "dkmix.trend", tr3[0], num * num / den, 1e-13);
  tk_nan(t, "dkmix.never", tr3[2]);
  tk_eq(t, "dkmix.df3", st[1], 2, 0);
  const double xn[1] = {-1}, x0[2] = {0, 0};
  tk_ieq(t, "dkmix.n", clab_dkmix(x, 0, 2, sc, a, sh, lo, hi, st, tr), CLAB_EARG);
  tk_ieq(t, "dkmix.a", clab_dkmix(x, 2, 2, sc, 0, sh, lo, hi, st, tr), CLAB_EDOM);
  tk_ieq(t, "dkmix.neg", clab_dkmix(xn, 1, 1, sc, a, sh, lo, hi, st, tr), CLAB_EDOM);
  tk_ieq(t, "dkmix.empty", clab_dkmix(x0, 1, 2, sc, a, sh, lo, hi, st, tr), CLAB_EDOM);
}

void life(tk &t) {
  const double lo[3] = {0, 3, 6}, hi[3] = {3, 6, inf}, d[3] = {40, 30, 20}, w[3] = {2, 2, 4};
  const double z = 1.959963984540054;
  double o[25];
  tk_ieq(t, "lftab.rc", clab_lftab(lo, hi, d, w, 3, z, 2, o), CLAB_OK);
  tk_eq(t, "lftab.enter", o[1], 56, 0);
  tk_eq(t, "lftab.atrisk", o[3], 97, 0);
  const double s0 = 57.0 / 97, s1 = s0 * 25 / 55;
  tk_eq(t, "lftab.surv", o[10], s1, 1e-15);
  const double g0 = 40.0 / (97 * 57);
  tk_eq(t, "lftab.se", o[12], s0 * std::sqrt(g0), 1e-15);
  const double u = std::sqrt(g0) / std::log(s0);
  tk_eq(t, "lftab.loglog", o[15], std::pow(s0, std::exp(-z * u)), 1e-14);
  tk_eq(t, "lftab.hazard", o[21], 40 / (3 * (97 - 20.0)), 1e-15);
  tk_nan(t, "lftab.openhz", o[23]);
  tk_eq(t, "lftab.median", o[24], 3 + 3 * (s0 - 0.5) / (s0 - s1), 1e-14);
  clab_lftab(lo, hi, d, w, 3, z, 0, o);
  tk_eq(t, "lftab.plain", o[15], s0 - z * s0 * std::sqrt(g0), 1e-14);
  clab_lftab(lo, hi, d, w, 3, z, 1, o);
  tk_eq(t, "lftab.log", o[15], s0 * std::exp(-z * std::sqrt(g0)), 1e-14);
  // every case ends in the first band: survival 0, no interval
  const double da[1] = {5}, wa[1] = {0}, h1[1] = {2}, l1[1] = {0};
  clab_lftab(l1, h1, da, wa, 1, z, 2, o);
  tk_true(t, "lftab.all", o[3] == 0 && std::isnan(o[5]) && o[8] == 1);
  // survival above one half throughout, or crossing in the open band
  const double dl[3] = {1, 1, 90};
  const double dz[2] = {1, 0};
  clab_lftab(lo, hi, dz, w, 2, z, 2, o);
  tk_true(t, "lftab.never", std::isinf(o[16]));
  clab_lftab(lo, hi, dl, w, 3, z, 2, o);
  tk_true(t, "lftab.open", std::isinf(o[24]));
  const double hb[3] = {3, 7, inf}, ho[3] = {3, inf, 9}, dn[3] = {-1, 0, 0}, z3[3] = {0, 0, 0};
  tk_ieq(t, "lftab.k", clab_lftab(lo, hi, d, w, 0, z, 2, o), CLAB_EARG);
  tk_ieq(t, "lftab.ci", clab_lftab(lo, hi, d, w, 3, z, 3, o), CLAB_EARG);
  tk_ieq(t, "lftab.gap", clab_lftab(lo, hb, d, w, 3, z, 2, o), CLAB_EDOM);
  tk_ieq(t, "lftab.openmid", clab_lftab(lo, ho, d, w, 3, z, 2, o), CLAB_EDOM);
  tk_ieq(t, "lftab.neg", clab_lftab(lo, hi, dn, w, 3, z, 2, o), CLAB_EDOM);
  tk_ieq(t, "lftab.none", clab_lftab(lo, hi, z3, z3, 3, z, 2, o), CLAB_EDOM);
  tk_ieq(t, "lftab.order", clab_lftab(hi, lo, d, w, 1, z, 2, o), CLAB_EDOM);
}

void trans(tk &t) {
  const int fr[4] = {0, 0, 1, 2}, to[4] = {0, 1, 2, 0};
  double c[9], p[9], pi[3], st[3], wk[12];
  tk_ieq(t, "trnmat.rc", clab_trnmat(fr, to, 4, 3, c, p, pi, st, wk), CLAB_OK);
  tk_eq(t, "trnmat.p", p[3], 0.5, 0);
  tk_eq(t, "trnmat.p2", p[2], 1, 0);
  const double pw[3] = {0.5, 0.25, 0.25};
  tk_veq(t, "trnmat.pi", pi, pw, 3, 1e-15);
  tk_eq(t, "trnmat.stay", st[1], 0.25, 0);
  tk_eq(t, "trnmat.theil", st[2], 0, 1e-15);
  // a state never left; two closed classes; no moves
  const int f1[1] = {0}, t1[1] = {1}, f2[2] = {0, 1};
  clab_trnmat(f1, t1, 1, 2, c, p, pi, st, wk);
  tk_nan(t, "trnmat.absorb", pi[0]);
  clab_trnmat(f2, f2, 2, 2, c, p, pi, st, wk);
  tk_true(t, "trnmat.closed", std::isnan(pi[0]) && std::isnan(st[2]) && st[1] == 1);
  clab_trnmat(f2, f2, 0, 2, c, p, pi, st, wk);
  tk_nan(t, "trnmat.none", st[1]);
  const int fb[1] = {5};
  tk_ieq(t, "trnmat.k", clab_trnmat(fr, to, 4, 0, c, p, pi, st, wk), CLAB_EARG);
  tk_ieq(t, "trnmat.from", clab_trnmat(fb, to, 1, 3, c, p, pi, st, wk), CLAB_EARG);
}

}  // namespace

void t_flw(tk &t, int big) {
  (void)big;
  stock(t);
  stay(t);
  docket(t);
  mix(t);
  life(t);
  trans(t);
}
