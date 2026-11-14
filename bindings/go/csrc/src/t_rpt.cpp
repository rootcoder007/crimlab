// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for clab_rptc (src/rpt.cpp). Part of the self-test (selftest.cpp).

#include <cstdint>

#include "crimlab.h"
#include "tchk.h"

namespace {

struct xrng {
  std::uint64_t s;
  double unif() {
    s ^= s >> 12;
    s ^= s << 25;
    s ^= s >> 27;
    return static_cast<double>((s * 2685821657736338717ULL) >> 11) * 0x1p-53;
  }
};

// the panel of dev/ref/cpp/rpt.R: units 0..39, years 2000..2009,
// months missing (u * 7 + y * 3 + (u * y) % 5) % 13, no record when
// (u + 2 y) % 17 == 0; rep = 12 - missing
int panel(int *u, int *y, int *rep) {
  int n = 0;
  for (int a = 0; a < 40; ++a)
    for (int b = 0; b < 10; ++b) {
      if ((a + 2 * b) % 17 == 0) continue;
      u[n] = a;
      y[n] = 2000 + b;
      rep[n] = 12 - (a * 7 + b * 3 + (a * b) % 5) % 13;
      ++n;
    }
  return n;
}

void run(tk &tt) {
  int flag[40], work[400], nc = 0;
  // 1. A and B report in 2000 and 2002 and miss all of 2001: neither is a
  // constant reporter.
  {
    const int u[6] = {0, 1, 0, 1, 0, 1}, y[6] = {2000, 2000, 2001, 2001, 2002, 2002},
              r[6] = {12, 12, 0, 0, 12, 12};
    tk_true(tt, "rpt.gap.rc", clab_rptc(u, y, r, 6, 2, 2000, 2002, 12, flag, &nc, work) == CLAB_OK);
    tk_true(tt, "rpt.gap", nc == 0 && flag[0] == 0 && flag[1] == 0);
    // A passes 2000, B 2001, both 2002: neither reported every year
    const int u2[5] = {0, 1, 1, 0, 1}, y2[5] = {2000, 2000, 2001, 2002, 2002},
              r2[5] = {12, 0, 12, 12, 12};
    clab_rptc(u2, y2, r2, 5, 2, 2000, 2002, 12, flag, &nc, work);
    tk_true(tt, "rpt.reset", nc == 0);
  }
  // 2. the deterministic panel at four thresholds: no unit at m = 3 and
  // m = 6, {14} at m = 9 and {2, 4, ..., 38} at m = 12 (reference sets from
  // crimeutils 0.5.1 at m = 9 and m = 12, dev/ref/cpp/rpt.R)
  {
    static int u[400], y[400], r[400];
    const int n = panel(u, y, r);
    const int want12[16] = {2, 4, 6, 8, 10, 12, 14, 19, 21, 23, 25, 27, 29, 31, 36, 38};
    clab_rptc(u, y, r, n, 40, 2000, 2009, 12 - 3, flag, &nc, work);
    tk_true(tt, "rpt.p3", nc == 0);
    clab_rptc(u, y, r, n, 40, 2000, 2009, 12 - 6, flag, &nc, work);
    tk_true(tt, "rpt.p6", nc == 0);
    clab_rptc(u, y, r, n, 40, 2000, 2009, 12 - 9, flag, &nc, work);
    tk_true(tt, "rpt.p9", nc == 1 && flag[14] == 1);
    clab_rptc(u, y, r, n, 40, 2000, 2009, 0, flag, &nc, work);
    tk_true(tt, "rpt.p12.n", nc == 16);
    for (int k = 0; k < 16; ++k) tk_true(tt, "rpt.p12", flag[want12[k]] == 1);
  }
  // 3. brute force on 15 random panels with duplicates and out-of-range years:
  // a unit is constant iff for each year some record qualifies; and the
  // result is monotone in minrep
  {
    xrng g{12};
    static int u[3000], y[3000], r[3000];
    int bad = 0, mono = 0;
    for (int rep = 0; rep < 15; ++rep) {
      const int n = 1 + static_cast<int>(g.unif() * 1000);
      for (int i = 0; i < n; ++i) {
        u[i] = static_cast<int>(g.unif() * 40);
        y[i] = 1998 + static_cast<int>(g.unif() * 14);
        r[i] = static_cast<int>(g.unif() * 13);
      }
      int prev[40] = {0};
      for (int m = 12; m >= 0; --m) {
        clab_rptc(u, y, r, n, 40, 2000, 2009, m, flag, &nc, work);
        for (int a = 0; a < 40; ++a) {
          int all = 1;
          for (int yy = 2000; yy <= 2009; ++yy) {
            int any = 0;
            for (int i = 0; i < n; ++i) any |= u[i] == a && y[i] == yy && r[i] >= m;
            all &= any;
          }
          bad += flag[a] != all;
          mono += prev[a] && !flag[a];
          prev[a] = flag[a];
        }
      }
    }
    tk_true(tt, "rpt.brute", bad == 0);
    tk_true(tt, "rpt.monotone", mono == 0);
  }
  // 4. guards and empty input
  {
    const int u[1] = {3}, y[1] = {2000}, r[1] = {12};
    tk_true(tt, "rpt.guard.unit", clab_rptc(u, y, r, 1, 2, 2000, 2001, 12, flag, &nc, work) == CLAB_EARG);
    tk_true(tt, "rpt.guard.years", clab_rptc(u, y, r, 1, 4, 2001, 2000, 12, flag, &nc, work) == CLAB_EARG);
    tk_true(tt, "rpt.empty", clab_rptc(u, y, r, 0, 4, 2000, 2001, 12, flag, &nc, work) == CLAB_OK && nc == 0);
  }
}


}  // namespace

void t_rpt(tk &t, int) { run(t); }

#ifdef CLAB_MAIN
#include <chrono>
#include <cstdio>

// 20,000 agencies x 30 years (crimeutils::ucr_constant_reporter_oris
// 0.465 s median of 3, dev/ref/cpp/rpt.R)
void b_rpt() {
  const int na = 20000, ny = 30, n = na * ny;
  static int u[n], y[n], r[n], flag[na], work[n];
  xrng g{1};
  for (int i = 0; i < n; ++i) {
    u[i] = i / ny;
    y[i] = 1991 + i % ny;
    int miss = 0;
    for (int k = 0; k < 12; ++k) miss += g.unif() < 0.02;
    r[i] = 12 - miss;
  }
  double t[5];
  int nc = 0;
  for (int k = 0; k < 5; ++k) {
    const auto a = std::chrono::steady_clock::now();
    clab_rptc(u, y, r, n, na, 1991, 2020, 12, flag, &nc, work);
    t[k] = std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count();
  }
  std::printf("BENCH clab_rptc 20000 x 30: median %.6f s (k=5), %d constant\n", tk_med5(t), nc);
}
#endif
