// SPDX-License-Identifier: AGPL-3.0-or-later
//
// Self-test of the disparity kernels: veil of darkness, the risks behind an
// odds ratio, the outcome test, benchmark decomposition, dyadic rates. A
// check's id names the test it came from (file and line, e.g. disp3).
#include <cmath>
#include <cstdint>
#include <limits>

#include "crimlab.h"
#include "tchk.h"

namespace {

void vod(tk &t) {
  const double day[2] = {60, 40}, night[2] = {50, 50};
  double sd[2], sn[2], ratio[2], adj[2];
  int sign[2];
  tk_ieq(t, "disp3", clab_vod(day, night, 2, 0, nullptr, sd, sn, ratio, sign,
                              nullptr), CLAB_OK);
  const double r3[2] = {1, (40.0 / 60) / (50.0 / 50)};
  tk_veq(t, "disp3", ratio, r3, 2, 1e-10);
  const double s4[2] = {0.6, 0.4};
  tk_veq(t, "disp4", sd, s4, 2, tk_tol);
  tk_ieq(t, "disp5", sign[0], 1);
  tk_ieq(t, "disp5", sign[1], -1);
  // an exposure shift equal to the ratio explains all of it
  const double shift[2] = {1, 2.0 / 3};
  clab_vod(day, night, 2, 0, shift, sd, sn, ratio, sign, adj);
  const double a9[2] = {1, 1};
  tk_veq(t, "disp9", adj, a9, 2, 1e-10);
  // the second group as reference moves the baseline
  clab_vod(day, night, 2, 1, nullptr, sd, sn, ratio, sign, nullptr);
  const double r11[2] = {1.5, 1};
  tk_veq(t, "disp11", ratio, r11, 2, tk_tol);
  const double d15[2] = {5, 5}, n15[2] = {7, 7};
  clab_vod(d15, n15, 2, 0, nullptr, sd, sn, ratio, sign, nullptr);
  tk_ieq(t, "disp15", sign[0], 0);
  tk_ieq(t, "disp15", sign[1], 0);
  // a group never seen in daylight gives a ratio of 0
  const double d52[2] = {50, 0}, n52[2] = {40, 3};
  tk_ieq(t, "edge52", clab_vod(d52, n52, 2, 0, nullptr, sd, sn, ratio, sign,
                               nullptr), CLAB_OK);
  const double r52[2] = {1, 0};
  tk_veq(t, "edge52", ratio, r52, 2, 1e-10);
  const double d53[2] = {1, 1e6}, n53[2] = {1e6, 1};
  clab_vod(d53, n53, 2, 0, nullptr, sd, sn, ratio, sign, nullptr);
  const double r53[2] = {1, 1e12};
  tk_veq(t, "edge53", ratio, r53, 2, 1e-10);
}

void orrisk(tk &t) {
  const double ors[5] = {0.05, 0.5, 1, 3, 400}, qs[3] = {0.01, 0.4, 0.99},
               ss[3] = {0.02, 0.5, 0.9};
  for (double orr : ors) {
    for (double q : qs) {
      for (double s : ss) {
        double a = -1, b = -1;
        tk_ieq(t, "disp22", clab_orrisk(orr, q, s, &a, &b), CLAB_OK);
        tk_gt(t, "disp25", a, 0);
        tk_lt(t, "disp26", a, 1);
        tk_gt(t, "disp27", b, 0);
        tk_lt(t, "disp28", b, 1);
        tk_eq(t, "disp29", a * (1 - b) / (b * (1 - a)), orr, 1e-9);
        tk_eq(t, "disp30", s * a + (1 - s) * b, q, 1e-12);
      }
    }
  }
  // a tiny odds ratio: OR - 1 must not round to -1
  const double tiny[4] = {1e-8, 1e-15, 1e-17, 1e-300},
               q2[3] = {0.01, 0.5, 0.99}, s2[3] = {0.02, 0.5, 0.98};
  for (double orr : tiny) {
    for (double q : q2) {
      for (double s : s2) {
        double a = -1, b = -1;
        tk_ieq(t, "edge78", clab_orrisk(orr, q, s, &a, &b), CLAB_OK);
        tk_gt(t, "edge81", a, 0);
        tk_lt(t, "edge82", a, 1);
        tk_gt(t, "edge83", b, 0);
        tk_le(t, "edge85", b, 1);
        tk_eq(t, "edge86", s * a + (1 - s) * b, q, 1e-12);
      }
    }
  }
  // the mirror image: the reciprocal with the groups exchanged
  double la, lb, ha, hb;
  clab_orrisk(1e-6, 0.3, 0.2, &la, &lb);
  clab_orrisk(1e6, 0.3, 0.8, &ha, &hb);
  tk_eq(t, "edge93", la, hb, 1e-14);
  tk_eq(t, "edge94", lb, ha, 1e-14);
  tk_ieq(t, "edge96", clab_orrisk(1e-320, 0.5, 0.5, &la, &lb), CLAB_EDOM);
}

void outtst(tk &t) {
  const double q[3] = {0.1, 0.5, 0.9}, w1[3] = {1, 1, 1}, w2[3] = {1, 3, 1};
  double o[2];
  clab_outtst(q, w1, 3, 0.3, o, o + 1);
  const double e39[2] = {2.0 / 3, 0.7};
  tk_veq(t, "disp39", o, e39, 2, tk_tol);
  clab_outtst(q, w2, 3, 0.3, o, o + 1);
  const double e40[2] = {0.8, (1.5 + 0.9) / 4};
  tk_veq(t, "disp40", o, e40, 2, tk_tol);
  const double q42[2] = {0.1, 0.2};
  clab_outtst(q42, w1, 2, 0.5, o, o + 1);
  tk_nan(t, "disp42", o[1]);
  // same threshold, different hit rates
  const double qa[2] = {0.2, 0.9}, qb[2] = {0.2, 0.6};
  double h1, h2;
  clab_outtst(qa, w1, 2, 0.5, o, &h1);
  clab_outtst(qb, w1, 2, 0.5, o, &h2);
  tk_gt(t, "disp44", h1, h2);
  // a stop exactly at the threshold is not searched: q > t
  const double qt[2] = {0.3, 0.5};
  clab_outtst(qt, w1, 2, 0.3, o, o + 1);
  tk_eq(t, "outtst.t", o[0], 0.5, 1e-15);
  tk_eq(t, "outtst.t", o[1], 0.5, 1e-15);
}

void dspbm(tk &t) {
  const double pop[2] = {9000, 1000}, con[2] = {900, 300}, frc[2] = {9, 6};
  double cd[2], fd[2], rd[2];
  clab_dspbm(pop, con, frc, 2, 0, cd, fd, rd);
  const double e50[2] = {1, 3}, e51[2] = {1, 2};
  tk_veq(t, "disp50", cd, e50, 2, 1e-10);
  tk_veq(t, "disp51", fd, e51, 2, 1e-10);
  const double e52[2] = {cd[0] * fd[0], cd[1] * fd[1]};
  tk_veq(t, "disp52", rd, e52, 2, 1e-10);
}

void dyadrt(tk &t) {
  // counts proportional to the pairs exposed: random mixing; dyad k has
  // offender k % 3 and victim k / 3, R's column-major order
  const double pop[3] = {800, 150, 50};
  double cnt[9], ro[9], rp[9], no[9], ra[9];
  int off[9], vic[9];
  for (int k = 0; k < 9; ++k) {
    off[k] = k % 3;
    vic[k] = k / 3;
    cnt[k] = pop[off[k]] * pop[vic[k]] * 0.002;
  }
  clab_dyadrt(cnt, off, vic, 9, pop, 3, ro, rp, no, ra);
  const double one[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
  tk_veq(t, "disp58", ra, one, 9, 1e-10);
  tk_veq(t, "disp59", ro, no, 9, 1e-10);
  // per offender rises with the victim group, per pair does not
  tk_gt(t, "disp61", ro[0], ro[6]);
  tk_eq(t, "disp62", rp[0], rp[6], 1e-10);
}

// The kernels check their own input (test-crimlab.R): a caller that skips
// the R layer gets an error code.
void guards(tk &t) {
  double o[8];
  int si[2];
  const double one2[2] = {1, 1}, m1[2] = {-1, 1}, z1[2] = {1, 0},
               z0[2] = {0, 1};
  tk_ieq(t, "crim9", clab_vod(one2, one2, 2, 2, nullptr, o, o, o, si, o),
         CLAB_EARG);
  tk_ieq(t, "crim10", clab_vod(one2, one2, 2, -1, nullptr, o, o, o, si, o),
         CLAB_EARG);
  tk_ieq(t, "crim11", clab_vod(m1, one2, 2, 1, nullptr, o, o, o, si, o),
         CLAB_EDOM);
  tk_ieq(t, "crim12", clab_vod(one2, z1, 2, 0, nullptr, o, o, o, si, o),
         CLAB_EDOM);
  tk_ieq(t, "crim13", clab_vod(one2, one2, 2, 0, z1, o, o, o, si, o),
         CLAB_EDOM);
  tk_ieq(t, "crim14", clab_vod(z0, one2, 2, 0, nullptr, o, o, o, si, o),
         CLAB_EDOM);
  tk_ieq(t, "crim15", clab_orrisk(0, 0.5, 0.5, o, o + 1), CLAB_EDOM);
  tk_ieq(t, "crim16", clab_orrisk(2, 1, 0.5, o, o + 1), CLAB_EDOM);
  tk_ieq(t, "crim17", clab_orrisk(2, 0.5, 0, o, o + 1), CLAB_EDOM);
  const double h = 0.5, w1 = 1, wm = -1, w0 = 0, big = 1.5;
  tk_ieq(t, "crim18", clab_outtst(&h, &w1, 0, 0.5, o, o + 1), CLAB_EARG);
  tk_ieq(t, "crim19", clab_outtst(&h, &w1, 1, INFINITY, o, o + 1), CLAB_EDOM);
  tk_ieq(t, "crim20", clab_outtst(&big, &w1, 1, 0.5, o, o + 1), CLAB_EDOM);
  tk_ieq(t, "crim21", clab_outtst(&h, &wm, 1, 0.5, o, o + 1), CLAB_EDOM);
  tk_ieq(t, "crim22", clab_outtst(&h, &w0, 1, 0.5, o, o + 1), CLAB_EDOM);
  tk_ieq(t, "crim23", clab_dspbm(&w1, &w1, &w1, 1, 1, o, o, o), CLAB_EARG);
  tk_ieq(t, "crim24", clab_dspbm(&w1, &w0, &w1, 1, 0, o, o, o), CLAB_EDOM);
  const int zero = 0, onei = 1, mone = -1;
  tk_ieq(t, "crim26", clab_dyadrt(&w1, &zero, &zero, 0, &w1, 1, o, o, o, o),
         CLAB_EARG);
  tk_ieq(t, "crim28", clab_dyadrt(&w1, &zero, &zero, 1, &w1, 0, o, o, o, o),
         CLAB_EARG);
  tk_ieq(t, "crim29", clab_dyadrt(&w1, &zero, &zero, 1, &w0, 1, o, o, o, o),
         CLAB_EDOM);
  tk_ieq(t, "crim30", clab_dyadrt(&w1, &onei, &zero, 1, &w1, 1, o, o, o, o),
         CLAB_EARG);
  tk_ieq(t, "crim31", clab_dyadrt(&w1, &zero, &mone, 1, &w1, 1, o, o, o, o),
         CLAB_EARG);
  tk_ieq(t, "crim32", clab_dyadrt(&wm, &zero, &zero, 1, &w1, 1, o, o, o, o),
         CLAB_EDOM);
}

// Stops simulated with a known visibility factor: group b (exposure 0.3
// against 0.7) is stopped 1.5 times as often in daylight, so a share
// 0.45 / 1.15 of daylight stops and 0.3 of night stops are of group b, n
// of each. Returns the estimated factor.
double vodsim(rng &g, int n) {
  double day[2] = {0, 0}, night[2] = {0, 0};
  for (int i = 0; i < n; ++i) {
    day[runi(g) < 0.45 / 1.15 ? 1 : 0] += 1;
    night[runi(g) < 0.3 ? 1 : 0] += 1;
  }
  double sd[2], sn[2], r[2];
  int si[2];
  clab_vod(day, night, 2, 0, nullptr, sd, sn, r, si, nullptr);
  return r[1];
}

//' @srrstats {G5.6} vodsim() simulates stops with a known visibility
//'   factor of 1.5 and clab_vod() recovers it.
//' @srrstats {G5.6a} The recovered factor is checked to a relative
//'   tolerance of 2 per cent, about three standard errors.
//' @srrstats {G5.6b} The recovery is checked under five seeds.
//' @srrstats {G5.7} The mean error over 20 seeds is smaller with 100 times
//'   more stops.
void recover(tk &t) {
  for (std::uint64_t s = 1; s <= 5; ++s) {
    rng g = {s};
    tk_eq(t, "recover", vodsim(g, 200000), 1.5, 0.02);
  }
  double small = 0, large = 0;
  for (std::uint64_t s = 1; s <= 20; ++s) {
    rng g = {100 + s};
    small += std::fabs(vodsim(g, 1000) - 1.5);
    large += std::fabs(vodsim(g, 100000) - 1.5);
  }
  tk_lt(t, "recover", large, small);
}

// x with every element moved by one part in 2^52, alternately up and down.
void nudge(const double *x, double *y, int n) {
  const double e = std::numeric_limits<double>::epsilon();
  for (int i = 0; i < n; ++i) y[i] = x[i] * (i % 2 ? 1 - e : 1 + e);
}

//' @srrstats {G5.9a} Inputs moved by .Machine$double.eps leave the results
//'   of clab_vod(), clab_cllr() and clab_clrsmp() unchanged to 1e-12.
void noise(tk &t) {
  const double day[2] = {600, 400}, night[2] = {700, 300};
  double dn[2], a[2], b[2], r0[2], r1[2];
  int si[2];
  nudge(day, dn, 2);
  clab_vod(day, night, 2, 0, nullptr, a, b, r0, si, nullptr);
  clab_vod(dn, night, 2, 0, nullptr, a, b, r1, si, nullptr);
  tk_veq(t, "noise", r1, r0, 2, 1e-12);
  const double s[4] = {50, 200, 8, 0.5}, f[4] = {0.01, 0.2, 1.5, 0.001};
  double s2[4], f2[4], c0 = 0, c1 = 0;
  nudge(s, s2, 4);
  nudge(f, f2, 4);
  clab_cllr(s, 4, f, 4, &c0);
  clab_cllr(s2, 4, f2, 4, &c1);
  tk_eq(t, "noise", c1, c0, 1e-12);
  const double ca[2] = {90, 20}, ra[2] = {100, 50}, cb[2] = {19, 60},
               rb[2] = {20, 140};
  double ca2[2], w[6][2], t0[5], t1[5];
  int cnt[3];
  nudge(ca, ca2, 2);
  clab_clrsmp(ca, ra, cb, rb, 2, w[0], w[1], w[2], w[3], w[4], w[5], t0, cnt);
  clab_clrsmp(ca2, ra, cb, rb, 2, w[0], w[1], w[2], w[3], w[4], w[5], t1, cnt);
  tk_veq(t, "noise", t1, t0, 5, 1e-12);
}

}  // namespace

void t_disp(tk &t) {
  vod(t);
  orrisk(t);
  outtst(t);
  dspbm(t);
  dyadrt(t);
  guards(t);
  recover(t);
  noise(t);
}
