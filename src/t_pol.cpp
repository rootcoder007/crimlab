// SPDX-License-Identifier: AGPL-3.0-or-later
//
// Self-test of the policing kernels: clearance decomposition and Simpson's
// paradox, exposure bounds, misclassification, regression to the mean.
#include <cmath>

#include "crimlab.h"
#include "tchk.h"

namespace {

struct clr {
  double ra[8], rb[8], wa[8], wb[8], re[8], me[8], tot[5];
  int cnt[3];
  int code;
};

clr clrsmp(const double *ca, const double *na, const double *cb,
           const double *nb, int n) {
  clr r{};
  r.code = clab_clrsmp(ca, na, cb, nb, n, r.ra, r.rb, r.wa, r.wb, r.re, r.me,
                       r.tot, r.cnt);
  return r;
}

void cnt(tk &t, const char *id, const clr &r, int up, int down, int kind) {
  tk_ieq(t, id, r.cnt[0], up);
  tk_ieq(t, id, r.cnt[1], down);
  tk_ieq(t, id, r.cnt[2], kind);
}

void clearance(tk &t) {
  // the toy Simpson case by hand
  const double ca[2] = {90, 20}, na[2] = {100, 50}, cb[2] = {19, 60},
               nb[2] = {20, 140};
  clr r = clrsmp(ca, na, cb, nb, 2);
  const double e6[2] = {0.9, 0.4}, e7[2] = {0.95, 60.0 / 140},
               e8[2] = {2.0 / 3, 1.0 / 3}, e9[2] = {1.0 / 8, 7.0 / 8},
               e10[2] = {110.0 / 150, 79.0 / 160};
  tk_veq(t, "pol6", r.ra, e6, 2, 1e-15);
  tk_veq(t, "pol7", r.rb, e7, 2, 1e-15);
  tk_veq(t, "pol8", r.wa, e8, 2, 1e-15);
  tk_veq(t, "pol9", r.wb, e9, 2, 1e-15);
  tk_veq(t, "pol10", r.tot, e10, 2, 1e-15);
  // rate terms (wa + wb) / 2 (rb - ra)
  const double e14[2] = {(2.0 / 3 + 1.0 / 8) / 2 * 0.05,
                         (1.0 / 3 + 7.0 / 8) / 2 * (60.0 / 140 - 0.4)};
  tk_veq(t, "pol14", r.re, e14, 2, 1e-12);
  // mix terms centred on the mean of the aggregates
  const double cen = (110.0 / 150 + 79.0 / 160) / 2;
  double e19[2];
  for (int k = 0; k < 2; ++k) e19[k] = ((e6[k] + e7[k]) / 2 - cen) * (e9[k] - e8[k]);
  tk_veq(t, "pol19", r.me, e19, 2, 1e-12);
  tk_eq(t, "pol20", r.tot[2], 79.0 / 160 - 110.0 / 150, 1e-15);
  tk_ieq(t, "pol21", r.cnt[0], 2);
  tk_ieq(t, "pol22", r.cnt[1], 0);
  tk_ieq(t, "pol23", r.cnt[2], 2);
  // every type up, so the fall is all mix (Crimlab.Clearance.change_ge_mix)
  tk_gt(t, "pol25", r.tot[3], 0);
  tk_lt(t, "pol26", r.tot[4], r.tot[2] + 1e-15);
  // every rate down, aggregate up
  r = clrsmp(cb, nb, ca, na, 2);
  cnt(t, "pol59", r, 0, 2, 2);
  // one part up, one down, the aggregate and the rate effect both fall
  const double c80[2] = {80, 10}, n100[2] = {100, 100},
               c62[2] = {85 * 0.5, 0.05 * 150}, n62[2] = {50, 150};
  r = clrsmp(c80, n100, c62, n62, 2);
  tk_eq(t, "pol65", r.tot[3],
        (0.5 + 0.25) / 2 * (0.85 - 0.8) + (0.5 + 0.75) / 2 * (0.05 - 0.1),
        1e-12);
  cnt(t, "pol66", r, 1, 1, 0);
  const double c68[2] = {85, 9.6};
  r = clrsmp(c80, n100, c68, n100, 2);
  tk_ieq(t, "pol69", r.cnt[2], 0);
  tk_eq(t, "pol70", r.tot[4], 0, 1e-15);
  // aggregate falls, rate effect rises, not every part up: kind 1
  const double c72[2] = {45, 7.2};
  r = clrsmp(c80, n100, c72, n62, 2);
  tk_gt(t, "pol73", r.tot[3], 0);
  tk_lt(t, "pol74", r.tot[2], 0);
  cnt(t, "pol75", r, 1, 1, 1);
  // no change in the counts: nothing moves
  const double c5[2] = {5, 6}, n10[2] = {10, 10};
  r = clrsmp(c5, n10, c5, n10, 2);
  const double z4[4] = {0, 0, 0, 0}, g79[4] = {r.tot[2], r.tot[3], r.tot[4],
                                               static_cast<double>(r.cnt[2])};
  tk_veq(t, "pol79", g79, z4, 4, 1e-15);
  // rates above one (clearances of earlier years' offences) are allowed
  const double c12 = 12, n10s = 10, c5s = 5;
  r = clrsmp(&c12, &n10s, &c5s, &n10s, 1);
  tk_eq(t, "pol82", r.tot[0], 1.2, tk_tol);
}

void clrrand(tk &t) {
  rng g = {11};
  for (int k = 0; k < 200; ++k) {
    const int n = rsam(g, 2, 6);
    double da[6], db[6], na[6], nb[6];
    for (int i = 0; i < n; ++i) {
      da[i] = runi(g, 1, 500);
      db[i] = runi(g, 1, 500);
      na[i] = da[i] * runi(g);
      nb[i] = db[i] * runi(g);
    }
    const clr r = clrsmp(na, da, nb, db, n);
    double sna = 0, sda = 0, snb = 0, sdb = 0, sre = 0, unc = 0, agg = 0;
    bool noneup = true;
    for (int i = 0; i < n; ++i) {
      sna += na[i];
      sda += da[i];
      snb += nb[i];
      sdb += db[i];
      sre += r.re[i];
      noneup = noneup && nb[i] / db[i] >= na[i] / da[i];
    }
    for (int i = 0; i < n; ++i) {
      unc += (na[i] / da[i] + nb[i] / db[i]) / 2 * (db[i] / sdb - da[i] / sda);
      agg += da[i] / sda * na[i] / da[i];
    }
    const double chg = snb / sdb - sna / sda;
    tk_eq(t, "pol41", r.tot[2], chg, 1e-12);
    tk_eq(t, "pol42", r.tot[3] + r.tot[4], chg, 1e-12);
    tk_eq(t, "pol43", r.tot[3], sre, 1e-12);
    // the uncentred mix effect has the same total (Crimlab.Clearance.kitagawa)
    tk_eq(t, "pol46", r.tot[4], unc, 1e-12);
    // Crimlab.Clearance.agg_eq_weighted
    tk_eq(t, "pol48", r.tot[0], agg, 1e-12);
    // Crimlab.Clearance.change_ge_mix whenever no rate falls
    if (noneup) tk_ge(t, "pol51", r.tot[2] - r.tot[4], -1e-12);
  }
}

struct eb {
  double rate[3], lo[3], hi[3], ratio[3], rlo[3], rhi[3];
  int sign[3];
};

eb expbnd(const double *y, const double *m, int n, int ref, double gam) {
  eb r{};
  clab_expbnd(y, m, n, ref, gam, r.rate, r.lo, r.hi, r.ratio, r.rlo, r.rhi,
              r.sign);
  return r;
}

void exposure(tk &t) {
  const double y[2] = {100, 300}, m[2] = {1000, 1000};
  eb r = expbnd(y, m, 2, 0, 1.5);
  const double e87[2] = {0.1, 0.3}, e88[2] = {0.1 / 1.5, 0.3 / 1.5},
               e89[2] = {0.1 * 1.5, 0.3 * 1.5}, e90[2] = {1, 3},
               e91[2] = {1 / 2.25, 3 / 2.25}, e92[2] = {1 * 2.25, 3 * 2.25};
  tk_veq(t, "pol87", r.rate, e87, 2, tk_tol);
  tk_veq(t, "pol88", r.lo, e88, 2, 1e-15);
  tk_veq(t, "pol89", r.hi, e89, 2, 1e-15);
  tk_veq(t, "pol90", r.ratio, e90, 2, 1e-15);
  tk_veq(t, "pol91", r.rlo, e91, 2, 1e-15);
  tk_veq(t, "pol92", r.rhi, e92, 2, 1e-15);
  tk_ieq(t, "pol93", r.sign[0], 0);
  tk_ieq(t, "pol93", r.sign[1], 1);
  // the second group as reference; a group below the band
  r = expbnd(y, m, 2, 1, 1.5);
  tk_ieq(t, "pol96", r.sign[0], -1);
  tk_ieq(t, "pol96", r.sign[1], 0);
  // the band closes at gamma = 1, and a wide band leaves the sign open
  const double y23[2] = {2, 3}, y13[2] = {1, 3}, m11[2] = {1, 1};
  r = expbnd(y23, m11, 2, 0, 1);
  tk_ieq(t, "pol98", r.sign[0], 0);
  tk_ieq(t, "pol98", r.sign[1], 1);
  r = expbnd(y13, m11, 2, 0, 2);
  tk_ieq(t, "pol99", r.sign[0], 0);
  tk_ieq(t, "pol99", r.sign[1], 0);
  // random exposures inside the band give ratios inside [rlo, rhi]; the
  // ends come from putting the two groups at opposite ends of their bands
  const double yy[3] = {40, 90, 15}, mm[3] = {200, 250, 90}, gm = 1.7;
  const eb b = expbnd(yy, mm, 3, 1, gm);
  rng g = {3};
  for (int k = 0; k < 500; ++k) {
    double e[3];
    for (int i = 0; i < 3; ++i) {
      e[i] = mm[i] * std::exp(runi(g, -std::log(gm), std::log(gm)));
    }
    for (int i = 0; i < 3; ++i) {
      const double tr = (yy[i] / e[i]) / (yy[1] / e[1]), tr1 = yy[i] / e[i];
      tk_ge(t, "pol110", tr, b.rlo[i] * (1 - 1e-12));
      tk_le(t, "pol111", tr, b.rhi[i] * (1 + 1e-12));
      tk_ge(t, "pol113", tr1, b.lo[i] * (1 - 1e-12));
      tk_le(t, "pol114", tr1, b.hi[i] * (1 + 1e-12));
    }
  }
  const double lo1 = (yy[0] / (gm * mm[0])) / (yy[1] / (mm[1] / gm));
  const double hi1 = (yy[0] / (mm[0] / gm)) / (yy[1] / (gm * mm[1]));
  tk_eq(t, "pol118", b.rlo[0], lo1, 1e-14);
  tk_eq(t, "pol119", b.rhi[0], hi1, 1e-14);
}

// The measured difference from the 2 x 2 table, cell by cell.
double cells(double r1, double r0, double p, double se, double sp) {
  const double pa[2] = {p, 1 - p}, rec[2] = {se, 1 - sp}, risk[2] = {r1, r0};
  double n1 = 0, d1 = 0, n0 = 0, d0 = 0;
  for (int a = 0; a < 2; ++a) {
    n1 += pa[a] * rec[a] * risk[a];
    d1 += pa[a] * rec[a];
    n0 += pa[a] * (1 - rec[a]) * risk[a];
    d0 += pa[a] * (1 - rec[a]);
  }
  return n1 / d1 - n0 / d0;
}

void misclass(tk &t) {
  const double r1s[3] = {0, 0.3, 0.9}, r0s[2] = {0.1, 0.5},
               ps[3] = {0.05, 0.4, 0.9}, ses[4] = {0, 0.3, 0.8, 1},
               sps[3] = {0.2, 0.7, 1};
  for (double sp : sps) {
    for (double se : ses) {
      for (double p : ps) {
        for (double r0 : r0s) {
          for (double r1 : r1s) {
            const double d1 = p * se + (1 - p) * (1 - sp);
            if (d1 <= 0 || d1 >= 1) continue;
            double x[5];
            clab_misrd(r1, r0, p, se, sp, x);
            tk_eq(t, "pol143", x[1], cells(r1, r0, p, se, sp), 1e-12);
            // Crimlab.Misclass.binary_attenuation
            tk_eq(t, "pol145", x[1], x[2] * (r1 - r0), 1e-12);
            tk_eq(t, "pol147", x[3] + x[4], 1, 1e-15);
            if (se + sp >= 1) {
              // never reversed (binary_sign) and never larger
              tk_ge(t, "pol150", x[2], -1e-15);
              tk_le(t, "pol151", x[2], 1 + 1e-12);
            } else {
              tk_lt(t, "pol153", x[2], 0);
              // two checks either way: where se + sp is 1 to rounding the
              // branch taken differs between x87 and SSE2, the count not
              tk_le(t, "pol151", x[2], 1 + 1e-12);
            }
          }
        }
      }
    }
  }
  // perfect measurement changes nothing
  double x[5];
  clab_misrd(0.3, 0.1, 0.4, 1, 1, x);
  tk_eq(t, "pol157", x[1], 0.2, 1e-15);
}

void regmean(tk &t) {
  const double x1[6] = {10, 2, 7, 4, 9, 10}, x2[6] = {6, 3, 8, 5, 2, 6},
               w1[6] = {1, 1, 1, 1, 1, 1}, w2[5] = {2, 1, 1, 1, 1};
  double r[7], rw[7], rr[7];
  clab_rtmhs(x1, x2, w1, 5, 6, r);
  tk_eq(t, "pol164", r[0], 3, 1e-15);
  tk_eq(t, "pol165", r[1], ((6.0 - 10) + (8.0 - 7) + (2.0 - 9)) / 3, 1e-15);
  tk_eq(t, "pol166", r[2], 7.0 - 8, 1e-15);
  // selection on x1 > cut in the stacked data and its period-swapped copy
  double s = 0, ns = 0, l = 0, nl = 0;
  for (int i = 0; i < 10; ++i) {
    const double a = i < 5 ? x1[i] : x2[i - 5], b = i < 5 ? x2[i] : x1[i - 5];
    if (a > 6) {
      s += b - a;
      ++ns;
    }
    if (a < 6) {
      l += b - a;
      ++nl;
    }
  }
  tk_eq(t, "pol170", r[3], s / ns, 1e-15);
  tk_eq(t, "pol171", r[6], l / nl, 1e-15);
  tk_eq(t, "pol172", r[5], ((3.0 - 2) + (5.0 - 4)) / 2, 1e-15);
  tk_eq(t, "pol173", r[4], r[1] - r[3], 1e-15);
  // weights: a place of weight 2 counts as two copies
  clab_rtmhs(x1, x2, w2, 5, 6, rw);
  clab_rtmhs(x1, x2, w1, 6, 6, rr);
  tk_veq(t, "pol177", rw + 1, rr + 1, 6, 1e-14);
  // nothing above the cut in period 2, nothing below in either
  const double a179[2] = {5, 1}, b179[2] = {1, 1}, a180[2] = {5, 6},
               b180[2] = {7, 8};
  clab_rtmhs(a179, b179, w1, 2, 4, r);
  tk_nan(t, "pol179", r[2]);
  clab_rtmhs(a180, b180, w1, 2, 4, r);
  tk_nan(t, "pol180", r[5]);
}

void rtmrand(tk &t) {
  // Crimlab.RegMean signs hold on any data
  rng g = {19};
  // the one large buffer of the self-test: 2,000 places of the simulation
  static double x1[2000], x2[2000], srt[2000];
  double w[60];
  for (int k = 0; k < 300; ++k) {
    const int n = rsam(g, 5, 60);
    const double sc = runi(g, 0.5, 3);
    const int sh = rsam(g, -3, 6);
    for (int i = 0; i < n; ++i) {
      x1[i] = rpois(g, 5) * sc;
      x2[i] = rpois(g, 5) + sh;
      w[i] = runi(g);
      srt[i] = x1[i];
    }
    const double cut = tk_q7(srt, n, runi(g, 0.1, 0.8));
    double sw = 0;
    for (int i = 0; i < n; ++i) sw += x1[i] > cut ? w[i] : 0;
    if (sw <= 0) continue;
    double r[7];
    clab_rtmhs(x1, x2, w, n, cut, r);
    tk_le(t, "pol193", r[3], 1e-12);
    if (!std::isnan(r[6])) tk_ge(t, "pol194", r[6], -1e-12);
  }
  // under exchangeable periods with no treatment the selected places fall:
  // 400 replications, each the mean over 2,000 places; the mean change is
  // several standard errors below 0
  g = rng{7};
  double ch[400], m = 0, v = 0;
  for (double &c : ch) {
    for (int i = 0; i < 2000; ++i) {
      const double mu = rexp(g, 0.5) + rexp(g, 0.5);  // Gamma(2, rate 0.5)
      x1[i] = rpois(g, mu);
      x2[i] = rpois(g, mu);
      srt[i] = x1[i];
    }
    const double cut = tk_q7(srt, 2000, 0.9);
    for (double &v : srt) v = 1;  // now the weights
    double r[7];
    clab_rtmhs(x1, x2, srt, 2000, cut, r);
    c = r[1];
    m += c / 400;
  }
  for (double c : ch) v += (c - m) * (c - m) / 399;
  tk_lt(t, "pol206", m + 5 * std::sqrt(v) / std::sqrt(400.0), 0);
}

void guards(tk &t) {
  double o[16];
  int c3[3];
  const double m1 = -1, one = 1, z = 0;
  tk_ieq(t, "crim73", clab_clrsmp(&m1, &one, &one, &one, 1, o, o, o, o, o, o,
                                  o, c3), CLAB_EDOM);
  tk_ieq(t, "crim74", clab_clrsmp(&one, &z, &one, &one, 1, o, o, o, o, o, o,
                                  o, c3), CLAB_EDOM);
  tk_ieq(t, "crim75", clab_misrd(1.5, 0.2, 0.5, 0.9, 0.9, o), CLAB_EDOM);
  tk_ieq(t, "crim76", clab_misrd(0.5, 0.2, 1, 0.9, 0.9, o), CLAB_EDOM);
  tk_ieq(t, "crim77", clab_misrd(0.5, 0.2, 0.5, 0.9, -0.1, o), CLAB_EDOM);
  const double xn[2] = {1, NAN}, x12[2] = {1, 2}, w11[2] = {1, 1},
               w1m[2] = {1, -1};
  tk_ieq(t, "crim78", clab_rtmhs(xn, x12, w11, 2, 0.5, o), CLAB_EDOM);
  tk_ieq(t, "crim79", clab_rtmhs(x12, x12, w1m, 2, 0.5, o), CLAB_EDOM);
}

}  // namespace

void t_pol(tk &t) {
  clearance(t);
  clrrand(t);
  exposure(t);
  misclass(t);
  regmean(t);
  rtmrand(t);
  guards(t);
}
