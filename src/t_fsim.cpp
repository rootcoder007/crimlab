// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for clab_nordet and clab_fairsim (src/fsim.cpp). Part of
// the self-test (selftest.cpp).

#include <cmath>

#include "crimlab.h"
#include "tchk.h"

namespace {

int rec[100000], ar[100000], gr[100000], tr[100000], dt[100000];
double rk[100000];

void run(tk &tt, int big) {
  // 1. noisy-OR by hand: incidents (0, 0) and (10, 10); officers (3, 4)
  // at distance 5 (counted: the range is inclusive), (0, 1) and (20, 20)
  {
    const double cx[2] = {0, 10}, cy[2] = {0, 10}, ox[3] = {3, 0, 20}, oy[3] = {4, 1, 20};
    int k[2], det[2];
    double p[2];
    tk_true(tt, "nordet.rc", clab_nordet(cx, cy, 2, ox, oy, 3, 5, 0.5, 0, 0, k, p, det) == CLAB_OK);
    tk_ieq(tt, "nordet.k0", k[0], 2);
    tk_ieq(tt, "nordet.k1", k[1], 0);
    tk_near(tt, "nordet.p0", p[0], 0.75, 1e-16);
    tk_near(tt, "nordet.p1", p[1], 0, 0);
    clab_nordet(cx, cy, 2, ox, oy, 3, 5, 1, 1, 3, k, p, det);
    tk_true(tt, "nordet.certain", p[0] == 1 && det[0] == 1 && det[1] == 0);
    clab_nordet(cx, cy, 2, ox, oy, 0, 5, 0.5, 0, 0, k, p, det);
    tk_true(tt, "nordet.none", k[0] == 0 && p[0] == 0);
    tk_true(tt, "nordet.guard", clab_nordet(cx, cy, 2, ox, oy, 3, 0, 0.5, 0, 0, k, p, det) == CLAB_EDOM);
    tk_true(tt, "nordet.guard.p", clab_nordet(cx, cy, 2, ox, oy, 3, 1, 0, 0, 0, k, p, det) == CLAB_EDOM);
  }
  // 2. the simulator: segregated areas, rates and moments within four
  // standard errors of their rules
  {
    const int n = big ? 100000 : 40000;
    const double pr[3] = {2, 1, 1};
    tk_true(tt, "fairsim.rc",
            clab_fairsim(n, 3, pr, 11, 0.4, 0.3, 5, ar, gr, tr, dt, rk) == CLAB_OK);
    int seg = 0;
    double ng[3] = {0, 0, 0}, dg[3] = {0, 0, 0}, tg[3] = {0, 0, 0}, sg[3] = {0, 0, 0},
           s2[3] = {0, 0, 0}, lo = 500, hi = 0;
    for (int i = 0; i < n; ++i) {
      seg += ar[i] % 3 == gr[i] && ar[i] >= 0 && ar[i] < 11 ? 0 : 1;
      ng[gr[i]] += 1;
      dg[gr[i]] += dt[i];
      tg[gr[i]] += tr[i];
      sg[gr[i]] += rk[i];
      s2[gr[i]] += rk[i] * rk[i];
      lo = std::fmin(lo, rk[i]);
      hi = std::fmax(hi, rk[i]);
    }
    tk_ieq(tt, "fairsim.seg", seg, 0);
    auto prop = [&](const char *id, double x, double m, double p) {
      tk_lt(tt, id, std::fabs(x / m - p), 4 * std::sqrt(p * (1 - p) / m));
    };
    prop("fairsim.share", ng[0], n, 0.5);
    prop("fairsim.truth", tg[1], ng[1], 0.4);
    prop("fairsim.det0", dg[0], ng[0], 0.4);
    prop("fairsim.det2", dg[2], ng[2], 0.4 * 0.7);
    const double m2 = sg[2] / ng[2], v2 = s2[2] / ng[2] - m2 * m2;
    tk_lt(tt, "fairsim.risk", std::fabs(m2 - 280), 4 * 40 / std::sqrt(ng[2]));
    tk_lt(tt, "fairsim.sd", std::fabs(std::sqrt(v2) - 40), 1.5);
    tk_true(tt, "fairsim.clamp", lo >= 0 && hi <= 500);
    // the same seed, the same records
    clab_fairsim(1000, 3, pr, 11, 0.4, 0.3, 5, rec, gr, tr, dt, rk);
    int same = 1;
    for (int i = 0; i < 1000; ++i) same &= rec[i] == ar[i];
    tk_true(tt, "fairsim.seed", same == 1);
    // a recording rate above 1 is 1
    clab_fairsim(2000, 2, nullptr, 4, 0.9, -1, 1, ar, gr, tr, dt, rk);
    int all = 1;
    for (int i = 0; i < 2000; ++i) all &= gr[i] == 0 || dt[i] == 1;
    tk_true(tt, "fairsim.clamp.det", all == 1);
    tk_true(tt, "fairsim.guard", clab_fairsim(10, 3, pr, 2, 0.4, 0.3, 5, ar, gr, tr, dt, rk) == CLAB_EARG);
    tk_true(tt, "fairsim.guard.bias", clab_fairsim(10, 2, nullptr, 2, 0.4, 2, 5, ar, gr, tr, dt, rk) == CLAB_EDOM);
    const double z0[2] = {0, 0};
    tk_true(tt, "fairsim.guard.props", clab_fairsim(10, 2, z0, 2, 0.4, 0, 5, ar, gr, tr, dt, rk) == CLAB_EDOM);
  }
}

}  // namespace

void t_fsim(tk &t, int big) { run(t, big); }
