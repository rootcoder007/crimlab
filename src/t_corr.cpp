// SPDX-License-Identifier: AGPL-3.0-or-later
//
// Self-test of the corrections kernels: incapacitation, careers, selective
// labels, spell lengths, sentence-effect bounds, hazards after depletion.
#include <algorithm>
#include <cmath>

#include "crimlab.h"
#include "tchk.h"

namespace {

struct inc {
  double rate[64], prev[64], marg[64], agg[3];
};

inc incap(const double *lam, const double *q, const double *s, int n,
          const double *share) {
  inc r{};
  clab_incap(lam, q, s, n, share, r.rate, r.prev, r.marg, r.agg);
  return r;
}

void incap(tk &t) {
  const double lam[3] = {2, 10, 0.5}, q[3] = {0.1, 0.3, 1}, s[3] = {1, 2.5, 0};
  const inc r = incap(lam, q, s, 3, nullptr);
  // one cycle: 1 / q crimes in a free spell of 1 / (lam q) years, then s
  // years inside; the share prevented is the time inside over the cycle
  double cyc[3], pr[3], s1[3];
  for (int i = 0; i < 3; ++i) {
    cyc[i] = (1 / q[i]) / (1 / (lam[i] * q[i]) + s[i]);
    pr[i] = s[i] / (1 / (lam[i] * q[i]) + s[i]);
    s1[i] = s[i] + 1;
  }
  tk_veq(t, "corr9", r.rate, cyc, 3, 1e-14);
  tk_veq(t, "corr11", r.prev, pr, 3, 1e-14);
  // one more year, by difference
  const inc r1 = incap(lam, q, s1, 3, nullptr);
  double dm[3];
  for (int i = 0; i < 3; ++i) dm[i] = r.rate[i] - r1.rate[i];
  tk_veq(t, "corr14", r.marg, dm, 3, 1e-12);
  tk_eq(t, "corr16", r.prev[2], 0, 1e-15);
  for (int i = 0; i < 3; ++i) tk_lt(t, "corr17", r.prev[i], 1);
  // the documented monotone claims on a grid
  const double lams[4] = {0.2, 1, 5, 40}, qs[3] = {0.01, 0.2, 1};
  for (double l : lams) {
    for (double qq : qs) {
      double lv[61], qv[61], sv[61];
      for (int i = 0; i < 61; ++i) {
        lv[i] = l;
        qv[i] = qq;
        sv[i] = 0.5 * i;
      }
      const inc g = incap(lv, qv, sv, 61, nullptr);
      for (int i = 0; i < 60; ++i) {
        tk_lt(t, "corr24", g.rate[i + 1], g.rate[i]);
        // the marginal year buys less each time (tested, not a theorem)
        tk_lt(t, "corr26", g.marg[i + 1], g.marg[i]);
      }
      for (int i = 0; i < 61; ++i) tk_gt(t, "corr27", g.marg[i], 0);
    }
  }
  // surer sentences lower the rate
  const double l4[5] = {4, 4, 4, 4, 4}, q5[5] = {0.1, 0.2, 0.4, 0.8, 1},
               s2[5] = {2, 2, 2, 2, 2};
  const inc a = incap(l4, q5, s2, 5, nullptr);
  for (int i = 0; i < 4; ++i) tk_lt(t, "corr31", a.rate[i + 1], a.rate[i]);
  // a uniform sentence prevents a larger share of the high-rate group
  const double lu[4] = {0.5, 2, 10, 50}, qu[4] = {0.2, 0.2, 0.2, 0.2},
               su[4] = {1, 1, 1, 1};
  const inc b = incap(lu, qu, su, 4, nullptr);
  for (int i = 0; i < 3; ++i) tk_gt(t, "corr33", b.prev[i + 1], b.prev[i]);
  // aggregate over groups
  const double lg[2] = {2, 10}, qg[2] = {0.1, 0.1}, sg[2] = {1, 1},
               sh[2] = {0.8, 0.2};
  const inc c = incap(lg, qg, sg, 2, sh);
  const double fr = 0.8 * 2 + 0.2 * 10, in = 0.8 * 2 / 1.2 + 0.2 * 10 / 2.0;
  tk_eq(t, "corr41", c.agg[0], fr, 1e-15);
  tk_eq(t, "corr42", c.agg[1], in, 1e-15);
  tk_eq(t, "corr43", c.agg[2], 1 - in / fr, 1e-15);
}

void career(tk &t) {
  const double lam[10] = {12, 10, 8, 6, 5, 4, 3, 2, 2, 1};
  double o[5];
  clab_incapr(lam, 10, 2, 3, 0.25, o);
  tk_eq(t, "corr50", o[0], 8 + 6 + 5, 1e-15);
  tk_eq(t, "corr51", o[1], 0.75 * 19, 1e-15);
  tk_eq(t, "corr52", o[2], 3 * 8, 1e-15);
  tk_eq(t, "corr53", o[3], 3 * 4, 1e-15);
  tk_eq(t, "corr54", o[4], 6 + 5 + 4, 1e-15);
  // the shortest path that holds the sentence still has a later value
  const double l3[3] = {5, 4, 3}, e57[5] = {9, 9, 10, 6, 7};
  clab_incapr(l3, 3, 0, 2, 0, o);
  tk_veq(t, "corr57", o, e57, 5, 1e-15);
  // constant rate: the bounds meet
  const double l6[6] = {3, 3, 3, 3, 3, 3}, e60[5] = {12, 0, 12, 12, 12};
  clab_incapr(l6, 6, 1, 4, 1, o);
  tk_veq(t, "corr60", o, e60, 5, 1e-15);
  rng g = {2};
  for (int k = 0; k < 200; ++k) {
    double lv[20];
    const int len = rsam(g, 3, 20);
    for (int i = 0; i < len; ++i) lv[i] = rexp(g, 1);
    std::sort(lv, lv + len, [](double x, double y) { return x > y; });
    const int s = rsam(g, 1, len - 1);
    const int t0 = rsam(g, 0, len - 1 - s);
    clab_incapr(lv, len, t0, s, runi(g), o);
    tk_le(t, "corr68", o[3], o[0] * (1 + 1e-12));
    tk_le(t, "corr69", o[0], o[2] * (1 + 1e-12));
    tk_le(t, "corr70", o[4], o[0] * (1 + 1e-12));
    tk_le(t, "corr71", o[1], o[2] * (1 + 1e-12));
  }
}

void sellab(tk &t) {
  // y of the detained is not read; 0 here as the R layer passes it
  double y[8] = {1, 0, 1, 0, 0, 0, 0, 0};
  const int rel[8] = {1, 1, 1, 1, 0, 0, 0, 1}, rule[8] = {1, 1, 0, 1, 1, 1, 0, 0};
  const double w[8] = {1, 2, 1, 1, 0.5, 1.5, 1, 1};
  double o[6];
  clab_sellab(y, rel, rule, w, 8, o);
  tk_eq(t, "corr81", o[0], (1.0 + 1) / (1 + 2 + 1 + 1 + 1), 1e-15);
  // enumerate the outcomes of the two detained defendants the rule releases
  double lo = 1, hi = 0;
  for (int f = 0; f < 4; ++f) {
    double yy[8] = {1, 0, 1, 0, 0, 0, 0, 0};
    yy[4] = f & 1;
    yy[5] = (f >> 1) & 1;
    double num = 0, den = 0;
    for (int i = 0; i < 8; ++i) {
      num += rule[i] ? w[i] * yy[i] : 0;
      den += rule[i] ? w[i] : 0;
    }
    lo = std::min(lo, num / den);
    hi = std::max(hi, num / den);
  }
  tk_eq(t, "corr88", o[2], lo, 1e-15);
  tk_eq(t, "corr89", o[3], hi, 1e-15);
  tk_eq(t, "corr90", o[1], 2.0 / 6, 1e-15);
  // Crimlab.Select.unobserved_width
  tk_eq(t, "corr92", o[4], o[1], 1e-15);
  tk_eq(t, "corr93", o[5], 1.0 / 4, 1e-15);
  // a rule inside the judge's releases is identified
  const int in[8] = {1, 1, 0, 0, 0, 0, 0, 1};
  clab_sellab(y, rel, in, w + 0, 8, o);
  const double one8[8] = {1, 1, 1, 1, 1, 1, 1, 1};
  clab_sellab(y, rel, in, one8, 8, o);
  tk_eq(t, "corr96", o[4], 0, 1e-15);
  tk_eq(t, "corr97", o[2], 1.0 / 3, 1e-15);
  tk_eq(t, "corr98", o[5], o[2], 1e-15);
  // the naive score sits at an end of the interval when the observed part
  // never fails or always fails
  const int r4[4] = {1, 1, 0, 0}, u4[4] = {1, 0, 1, 0};
  const double y0[4] = {0, 0, 0, 0}, y1[4] = {1, 1, 0, 0};
  clab_sellab(y0, r4, u4, one8, 4, o);
  const double e104[3] = {0, 0.5, 0}, e106[3] = {0.5, 1, 1};
  double g3[3] = {o[2], o[3], o[5]};
  tk_veq(t, "corr104", g3, e104, 3, tk_tol);
  clab_sellab(y1, r4, u4, one8, 4, o);
  g3[0] = o[2];
  g3[1] = o[3];
  g3[2] = o[5];
  tk_veq(t, "corr106", g3, e106, 3, tk_tol);
  // a rule that releases only detained defendants
  const int r2[2] = {1, 0}, u2[2] = {0, 1};
  const double y2[2] = {1, 0};
  clab_sellab(y2, r2, u2, one8, 2, o);
  const double e109[2] = {0, 1};
  tk_veq(t, "corr109", o + 2, e109, 2, 1e-15);
  tk_nan(t, "corr110", o[5]);
}

void selrand(tk &t) {
  // Crimlab.Select.naive_between: the naive score never leaves the interval
  rng g = {9};
  for (int k = 0; k < 300; ++k) {
    const int n = rsam(g, 2, 40);
    int rel[40], rule[40];
    double y[40], w[40], o[6];
    bool anyr = false, anyu = false;
    for (int i = 0; i < n; ++i) {
      rel[i] = runi(g) < 0.6;
      rule[i] = runi(g) < 0.6;
      y[i] = rel[i] && runi(g) < 0.4 ? 1 : 0;
      w[i] = runi(g);
      anyr = anyr || rel[i];
      anyu = anyu || rule[i];
    }
    if (!anyr || !anyu) continue;
    clab_sellab(y, rel, rule, w, n, o);
    if (std::isnan(o[5])) continue;
    tk_ge(t, "corr122", o[5], o[2] - 1e-12);
    tk_le(t, "corr123", o[5], o[3] + 1e-12);
  }
}

void losbnd(tk &t) {
  const double done[5] = {30, 45, 60, 90, 120}, open[3] = {100, 150, 200};
  double o[4];
  clab_losbnd(done, 5, open, 3, o);
  tk_eq(t, "corr129", o[0], 69, 1e-15);
  tk_eq(t, "corr130", o[1], 150, 1e-15);
  tk_eq(t, "corr131", o[2], (345.0 + 450) / 8, 1e-15);
  // Crimlab.Censor.lower_bound_sub
  tk_eq(t, "corr133", o[3], o[2] - o[0], 1e-12);
  rng g = {8};
  for (int k = 0; k < 200; ++k) {
    double tv[30], av[30];
    const int n = rsam(g, 1, 30), m = rsam(g, 1, 30);
    for (int i = 0; i < n; ++i) tv[i] = rexp(g, 1.0 / 50);
    for (int i = 0; i < m; ++i) av[i] = rexp(g, 1.0 / 80);
    clab_losbnd(tv, n, av, m, o);
    // any completion u >= a has a mean at least the bound
    double tm = 0;
    for (int i = 0; i < n; ++i) tm += tv[i];
    for (int i = 0; i < m; ++i) tm += av[i] + rexp(g, 1.0 / 40);
    tm /= n + m;
    tk_ge(t, "corr142", tm, o[2] - 1e-12);
    tk_ge(t, "corr143", tm - o[0], o[3] - 1e-12);
  }
  // open spells younger than the completed mean: the bound lies below it
  const double d2[2] = {100, 200}, a1 = 10;
  clab_losbnd(d2, 2, &a1, 1, o);
  tk_lt(t, "corr146", o[3], 0);
}

void sntbnd(tk &t) {
  // z: "a" is 0, "b" (the treatment, last in sorted order) is 1
  const double y[6] = {1, 0, 1, 1, 0, 0}, w[6] = {1, 2, 1, 0.5, 1, 1.5};
  const int z[6] = {0, 0, 1, 1, 1, 0}, za[6] = {1, 1, 0, 0, 0, 1};
  double o[11], s[11];
  clab_sntbnd(y, z, w, 6, o);
  double ww = 0;
  for (double x : w) ww += x;
  const double e156[2] = {4.5 / ww, 2.5 / ww};
  tk_veq(t, "corr156", o + 2, e156, 2, 1e-15);
  // brute force: fill each person's missing outcome with 0 or 1
  double lo[3] = {1, 1, 1}, hi[3] = {-1, -1, -1};
  for (int f = 0; f < 64; ++f) {
    double e[3] = {0, 0, 0};
    for (int i = 0; i < 6; ++i) {
      const double fi = (f >> i) & 1;
      const double y1 = z[i] == 1 ? y[i] : fi, y0 = z[i] == 0 ? y[i] : fi;
      e[0] += w[i] * y0 / ww;
      e[1] += w[i] * y1 / ww;
      e[2] += w[i] * (y1 - y0) / ww;
    }
    for (int j = 0; j < 3; ++j) {
      lo[j] = std::min(lo[j], e[j]);
      hi[j] = std::max(hi[j], e[j]);
    }
  }
  const double gl[3] = {o[4], o[6], o[8]}, gh[3] = {o[5], o[7], o[9]};
  tk_veq(t, "corr164", gl, lo, 3, 1e-14);
  tk_veq(t, "corr165", gh, hi, 3, 1e-14);
  // Crimlab.Sentence.ate_width_one, ate_contains_zero
  tk_eq(t, "corr167", o[9] - o[8], 1, 1e-15);
  tk_le(t, "corr168", o[8], 0);
  tk_ge(t, "corr169", o[9], 0);
  // the other value as treatment flips the effect
  clab_sntbnd(y, za, w, 6, s);
  tk_eq(t, "corr175", s[8], -o[9], 1e-15);
  // an arm with zero weight leaves its observed rate undefined
  const double y2[2] = {1, 0}, w2[2] = {1, 0};
  const int z2[2] = {0, 1};
  clab_sntbnd(y2, z2, w2, 2, o);
  tk_nan(t, "corr178", o[10]);
}

void hzdsel(tk &t) {
  const double ah[2] = {0.5, 0.8}, al[2] = {0.9, 0.9};
  double o[5];
  clab_hzdsel(0.5, 0.5, 0.1, ah, al, o);
  tk_eq(t, "corr183", o[0], 5.0 / 14, 1e-15);
  tk_eq(t, "corr184", o[1], 8.0 / 17, 1e-15);
  tk_eq(t, "corr185", o[2], 17.0 / 70, 1e-15);
  tk_eq(t, "corr186", o[3], 49.0 / 170, 1e-15);
  tk_eq(t, "corr187", o[4], 343.0 / 289, 1e-14);
  // a cohort of 10,000: survivors and their reconvictions
  rng g = {6};
  for (int k = 0; k < 100; ++k) {
    const double s = runi(g), h = runi(g, 0.2, 1), l = runi(g, 0, h);
    double sh[2] = {runi(g, 0.05, 1), runi(g, 0.05, 1)};
    const double sl[2] = {runi(g, 0.05, 1), runi(g, 0.05, 1)};
    clab_hzdsel(s, h, l, sh, sl, o);
    double hz[2];
    for (int j = 0; j < 2; ++j) {
      const double hi = 1e4 * s * sh[j], lo = 1e4 * (1 - s) * sl[j];
      hz[j] = (hi * h + lo * l) / (hi + lo);
    }
    tk_veq(t, "corr200", o + 2, hz, 2, 1e-12);
    // Crimlab.Hazard.hr2_gt_one_of_depletion, same low-risk survival
    const bool differ = sh[0] != sh[1];
    std::sort(sh, sh + 2);
    const double s1[2] = {sl[0], sl[0]};
    clab_hzdsel(s, h, l, sh, s1, o);
    if (differ) tk_gt(t, "corr203", o[4], 1);
  }
}

void guards(tk &t) {
  double o[11];
  const double z = 0, h = 0.5, one = 1;
  tk_ieq(t, "crim80", clab_incap(&z, &h, &one, 1, nullptr, o, o, o, o),
         CLAB_EDOM);
  const double two = 2, q15 = 1.5, m1 = -1;
  tk_ieq(t, "crim81", clab_incap(&two, &q15, &one, 1, nullptr, o, o, o, o),
         CLAB_EDOM);
  tk_ieq(t, "crim82", clab_incap(&two, &h, &m1, 1, nullptr, o, o, o, o),
         CLAB_EDOM);
  const double y05[2] = {0, 0.5}, y01[2] = {0, 1}, w11[2] = {1, 1},
               w1m[2] = {1, -1};
  const int z01[2] = {0, 1}, z02[2] = {0, 2};
  tk_ieq(t, "crim83", clab_sntbnd(y05, z01, w11, 2, o), CLAB_EDOM);
  tk_ieq(t, "crim84", clab_sntbnd(y01, z02, w11, 2, o), CLAB_EDOM);
  tk_ieq(t, "crim85", clab_sntbnd(y01, z01, w1m, 2, o), CLAB_EDOM);
  const double a90[2] = {0.9, 0}, b98[2] = {0.9, 0.8}, b12[2] = {1.2, 0.8};
  tk_ieq(t, "crim86", clab_hzdsel(0.5, 0.4, 0.2, a90, b98, o), CLAB_EDOM);
  tk_ieq(t, "crim88", clab_hzdsel(0.5, 0.4, 0.2, b98, b12, o), CLAB_EDOM);
}

}  // namespace

void t_corr(tk &t) {
  incap(t);
  career(t);
  sellab(t);
  selrand(t);
  losbnd(t);
  sntbnd(t);
  hzdsel(t);
  guards(t);
}
