// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for clab_grprr, clab_dastd and clab_sirex (src/rat.cpp)
// against the formulas written out here. Part of the self-test
// (selftest.cpp).
#include <cmath>

#include "crimlab.h"
#include "tchk.h"

namespace {

const double z = 1.959963984540054;

void rat(tk &t) {
  double rt[3], lo[3], hi[3], rp[3], ra[3], rl[3], rh[3], df[3], dl[3], dh[3];
  // 1. proportions 15 / 60 and 40 / 80 against the second
  {
    const double x[2] = {15, 40}, n[2] = {60, 80};
    tk_ieq(t, "rat.rc", clab_grprr(x, n, 2, 1, 1, z, rt, lo, hi, rp, ra, rl, rh, df, dl, dh),
           CLAB_OK);
    // Wilson by hand
    const double p = 0.25, m = 60, c = (p + z * z / (2 * m)) / (1 + z * z / m);
    const double h = z * std::sqrt(p * (1 - p) / m + z * z / (4 * m * m)) / (1 + z * z / m);
    tk_eq(t, "rat.wil.lo", lo[0], c - h, 1e-14);
    tk_eq(t, "rat.wil.hi", hi[0], c + h, 1e-14);
    tk_eq(t, "rat.ratio", ra[0], 0.5, 1e-15);
    const double s = std::sqrt(1 / 15.0 - 1 / 60.0 + 1 / 40.0 - 1 / 80.0);
    tk_eq(t, "rat.katz", rl[0], 0.5 * std::exp(-z * s), 1e-14);
    tk_eq(t, "rat.newc", dl[0], -0.25 - std::sqrt((p - lo[0]) * (p - lo[0]) +
                                                  (hi[1] - 0.5) * (hi[1] - 0.5)),
          1e-14);
    tk_nan(t, "rat.refnan", rl[1]);
    tk_eq(t, "rat.rep", rp[0], (15.0 / 55) / (60.0 / 140), 1e-14);
  }
  // 2. rates with a zero: ratio 0 and Inf, no log-normal interval
  {
    const double x[3] = {0, 10, 5}, n[3] = {100, 200, 50};
    tk_ieq(t, "rat.rc2", clab_grprr(x, n, 3, 0, 0, z, rt, lo, hi, rp, ra, rl, rh, df, dl, dh),
           CLAB_OK);
    tk_true(t, "rat.inf", std::isinf(ra[1]) && std::isnan(rl[1]) && std::isnan(lo[0]));
    tk_eq(t, "rat.wald", dh[2], 0.1 + z * std::sqrt(5 / 2500.0), 1e-14);
    tk_eq(t, "rat.lnrate", hi[1], 0.05 * std::exp(z / std::sqrt(10.0)), 1e-14);
    tk_ieq(t, "rat.zeroref", clab_grprr(x, n, 3, 0, 0, z, rt, lo, hi, rp, ra, rl, rh, df, dl,
                                        dh),
           CLAB_OK);
  }
  // 3. refusals
  {
    const double x[2] = {5, 2}, n[2] = {3, 4};
    tk_ieq(t, "rat.xgtn", clab_grprr(x, n, 2, 0, 1, z, rt, lo, hi, rp, ra, rl, rh, df, dl, dh),
           CLAB_EDOM);
    const double xn[2] = {-1, 2};
    tk_ieq(t, "rat.neg", clab_grprr(xn, n, 2, 0, 0, z, rt, lo, hi, rp, ra, rl, rh, df, dl, dh),
           CLAB_EDOM);
    tk_ieq(t, "rat.ref", clab_grprr(x, n, 2, 2, 0, z, rt, lo, hi, rp, ra, rl, rh, df, dl, dh),
           CLAB_EARG);
  }
  // 4. direct standardisation: two strata, weights 3 : 1
  {
    const double x[4] = {10, 4, 6, 9}, n[4] = {100, 50, 120, 60}, w[2] = {3, 1};
    const int g[4] = {0, 0, 1, 1}, s[4] = {0, 1, 0, 1};
    double r[2], v[2], wm[2], ev[2];
    int seen[4];
    tk_ieq(t, "dst.rc", clab_dastd(x, n, g, s, 4, 2, 2, w, 0, z, r, v, wm, ev, ra, rl, rh, df,
                                   dl, dh, seen),
           CLAB_OK);
    tk_eq(t, "dst.r0", r[0], 0.75 * 0.1 + 0.25 * 0.08, 1e-15);
    tk_eq(t, "dst.v0", v[0], 0.75 * 0.75 * 10 / 1e4 + 0.0625 * 4 / 2500, 1e-15);
    tk_eq(t, "dst.wm", wm[1], 0.75 / 120, 1e-15);
    tk_eq(t, "dst.ratio", ra[1], (0.75 * 0.05 + 0.25 * 0.15) / r[0], 1e-14);
    // a standardised rate lies between its stratum rates
    tk_true(t, "dst.convex", r[1] >= 0.05 && r[1] <= 0.15);
    const int s2[4] = {0, 0, 0, 1};
    tk_ieq(t, "dst.dup", clab_dastd(x, n, g, s2, 4, 2, 2, w, 0, z, r, v, wm, ev, ra, rl, rh,
                                    df, dl, dh, seen),
           CLAB_EARG);
    const double w0[2] = {0, 0};
    tk_ieq(t, "dst.w0", clab_dastd(x, n, g, s, 4, 2, 2, w0, 0, z, r, v, wm, ev, ra, rl, rh,
                                   df, dl, dh, seen),
           CLAB_EDOM);
  }
  // 5. indirect standardisation and Marshall's shrinkage
  {
    const double x[6] = {12, 3, 60, 20, 4, 1}, n[6] = {4000, 2000, 8000, 3000, 1500, 900};
    const int g[6] = {0, 0, 1, 1, 2, 2}, s[6] = {0, 1, 0, 1, 0, 1};
    double o[3], e[3], sr[3], eb[3], sh[3], pr[2], wk[4];
    tk_ieq(t, "sir.rc", clab_sirex(x, n, g, s, 6, 3, 2, nullptr, o, e, sr, eb, sh, pr, wk),
           CLAB_OK);
    const double r0 = 76.0 / 13500, r1 = 24.0 / 5900;
    tk_eq(t, "sir.e0", e[0], 4000 * r0 + 2000 * r1, 1e-14);
    tk_eq(t, "sir.mean", pr[0], 1.0, 1e-14);  // pooled rates: obs total = exp total
    double s2 = 0;
    for (int i = 0; i < 3; ++i) s2 += e[i] * (sr[i] - 1) * (sr[i] - 1);
    s2 = s2 / 100 - 1 / (100.0 / 3);
    // these areas vary more than Poisson noise: the prior has a variance
    tk_gt(t, "sir.var", s2, 0);
    tk_eq(t, "sir.eb", eb[0], (o[0] + 1 / s2) / (e[0] + 1 / s2), 1e-12);
    const double rz[2] = {0, 0};
    tk_ieq(t, "sir.zero", clab_sirex(x, n, g, s, 6, 3, 2, rz, o, e, sr, eb, sh, pr, wk),
           CLAB_EDOM);
  }
}

}  // namespace

void t_rat(tk &t, int big) {
  (void)big;
  rat(t);
}
