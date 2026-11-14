// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for clab_chmb (src/chm.cpp). Built only with
// the self-test (selftest.cpp).

#include <cmath>
#include <cstdint>

#include "crimlab.h"
#include "tchk.h"

namespace {

// toolmaRk 0.0.2, R 4.5.0, dev/ref/cpp/chm.R section 2:
// chumbley_non_random(matrix(m1), matrix(m2), 300, 30, 0.25) (same source),
// (m1, m3, 300, 30, 0.25) (different sources) and (m1, m2) at the
// defaults 500, 50, 0.25: n_same, n_diff, U, p. The marks are built by
// marks() below.
const double ch12[4] = {29, 22, 6.066934341798361, 6.5187411024680841e-10};
const double ch13[4] = {33, 5, -2.0080258544570162, 0.97767972783597323};
const double ch12d[4] = {12, 12, 4.156921938165306, 1.6128207281251861e-05}; // defaults 500, 50, 0.25

double m1[1500], m2[1500], m3[1500], work[60000];

// MINSTD from 4242: u[0..7799]; walks w1 = cumsum(u[0..2599] - 1/2), w3 =
// cumsum(u[2600..5199] - 1/2), noise e = (u[5200..7799] - 1/2) / 5;
// m1 = w1[0..1499] + e[0..1499], m2 = w1[300..1799] + e[1000..2499],
// m3 = w3[0..1499] + e[0..1499]
void marks() {
  static double u[7800], w1[2600], w3[2600];
  std::int64_t s = 4242;
  for (int i = 0; i < 7800; ++i) {
    s = (48271 * s) % 2147483647;
    u[i] = static_cast<double>(s) / 2147483647.0;
  }
  double a = 0, b = 0;
  for (int i = 0; i < 2600; ++i) {
    a += u[i] - 0.5;
    b += u[2600 + i] - 0.5;
    w1[i] = a;
    w3[i] = b;
  }
  for (int i = 0; i < 1500; ++i) {
    const double e1 = (u[5200 + i] - 0.5) * 0.2, e2 = (u[5200 + 1000 + i] - 0.5) * 0.2;
    m1[i] = w1[i] + e1;
    m2[i] = w1[300 + i] + e2;
    m3[i] = w3[i] + e1;
  }
}

void run(tk &tt) {
  marks();
  int ns, nd;
  double u, p;
  // 1. parity with toolmaRk (same source, different sources, defaults).
  // The reference p-value is 1 - pnorm(U); here p is erfc(U / sqrt 2) / 2,
  // so tiny p-values are compared relative to 1e-6.
  tk_true(tt, "chm.rc", clab_chmb(m1, 1500, m2, 1500, 300, 30, 0.25, &ns, &nd, &u, &p, work) == CLAB_OK);
  tk_true(tt, "chm.12.n", ns == ch12[0] && nd == ch12[1]);
  tk_near(tt, "chm.12.U", u, ch12[2], 1e-12);
  tk_near(tt, "chm.12.p", p / ch12[3], 1, 1e-6);
  clab_chmb(m1, 1500, m3, 1500, 300, 30, 0.25, &ns, &nd, &u, &p, work);
  tk_true(tt, "chm.13.n", ns == ch13[0] && nd == ch13[1]);
  tk_near(tt, "chm.13.U", u, ch13[2], 1e-12);
  tk_near(tt, "chm.13.p", p, ch13[3], 1e-12);
  clab_chmb(m1, 1500, m2, 1500, 500, 50, 0.25, &ns, &nd, &u, &p, work);
  tk_true(tt, "chm.12d.n", ns == ch12d[0] && nd == ch12d[1]);
  tk_near(tt, "chm.12d.U", u, ch12d[2], 1e-12);
  tk_near(tt, "chm.12d.p", p / ch12d[3], 1, 1e-9);
  // 2. a bound from the closed form of the rank sum: its largest value,
  // with every same-shift correlation above every different-shift one, is
  // n m + n (n + 1) / 2, so U <= (n m / 2) / sqrt(n m (n + m + 1) / 12)
  // (ties only lower it)
  {
    clab_chmb(m1, 1500, m2, 1500, 300, 30, 0.25, &ns, &nd, &u, &p, work);
    const double n = ns, m = nd, umax = (n * m / 2) / std::sqrt(n * m * (n + m + 1) / 12);
    tk_true(tt, "chm.bound", u <= umax + 1e-12);
  }
  // 3. defined behaviour: a 300-point profile with the 500-point window is
  // CLAB_EARG, a missing value CLAB_EDOM.
  tk_true(tt, "chm.short", clab_chmb(m1, 300, m2, 300, 500, 50, 0.25, &ns, &nd, &u, &p, work) == CLAB_EARG);
  {
    static double mn[1500];
    for (int i = 0; i < 1500; ++i) mn[i] = m1[i];
    mn[699] = std::nan("");
    tk_true(tt, "chm.nan", clab_chmb(mn, 1500, m2, 1500, 300, 30, 0.25, &ns, &nd, &u, &p, work) == CLAB_EDOM);
    for (int i = 0; i < 1500; ++i) mn[i] = 1.0;
    tk_true(tt, "chm.flat", clab_chmb(mn, 1500, m2, 1500, 300, 30, 0.25, &ns, &nd, &u, &p, work) == CLAB_EDOM);
  }
  // 3b. R's round() takes halves to even: 0.01 n = 12.5 for n = 1,250
  // trims from position 12 (not 13); and a pair too short for any
  // validation window on one side gives NaN with zero counts
  {
    tk_true(tt, "chm.half", clab_chmb(m1, 1250, m2, 1250, 300, 30, 0.25, &ns, &nd, &u, &p, work) == CLAB_OK &&
                                std::isfinite(u));
    tk_true(tt, "chm.empty", clab_chmb(m1, 330, m2, 330, 300, 30, 0.25, &ns, &nd, &u, &p, work) == CLAB_OK &&
                                 (ns == 0 || nd == 0) && std::isnan(u) && std::isnan(p));
  }
  // 4. symmetry of the detrending: adding a smooth trend to a mark leaves
  // the score almost unchanged (lowess removes it); exact for a constant
  {
    static double mc[1500];
    for (int i = 0; i < 1500; ++i) mc[i] = m2[i] + 7.0;
    double u2, p2;
    int a, b;
    clab_chmb(m1, 1500, m2, 1500, 300, 30, 0.25, &ns, &nd, &u, &p, work);
    clab_chmb(m1, 1500, mc, 1500, 300, 30, 0.25, &a, &b, &u2, &p2, work);
    tk_true(tt, "chm.shift", a == ns && b == nd && std::fabs(u2 - u) < 1e-9);
  }
}


}  // namespace

void t_chm(tk &t, int) { run(t); }

#ifdef CLAB_MAIN
#include <chrono>
#include <cstdio>

// timing: 2,000- and 9,500-point profiles at the defaults (toolmaRk took
// 0.423 s and 5.267 s, medians on zeus)
void b_chm() {
  static double a[9500], b[9500], w[250000];
  std::int64_t s = 7;
  double c = 0;
  static double walk[10000];
  for (int i = 0; i < 10000; ++i) {
    s = (48271 * s) % 2147483647;
    c += static_cast<double>(s) / 2147483647.0 - 0.5;
    walk[i] = c;
  }
  const int len[2] = {2000, 9500}, sh[2] = {100, 500};
  for (int q = 0; q < 2; ++q) {
    for (int i = 0; i < len[q]; ++i) {
      a[i] = walk[i];
      b[i] = walk[i + sh[q]];
    }
    double t[5], u, p;
    int ns, nd;
    for (int k = 0; k < 5; ++k) {
      const auto t0 = std::chrono::steady_clock::now();
      clab_chmb(a, len[q], b, len[q], 500, 50, 0.25, &ns, &nd, &u, &p, w);
      t[k] = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    }
    std::printf("BENCH clab_chmb %d-point profiles: median %.4f s (k=5)\n", len[q], tk_med5(t));
  }
}
#endif
