// SPDX-License-Identifier: AGPL-3.0-or-later
//
// Standalone tests for clab_skp (src/sk.cpp). Built only with -DCLAB_MAIN.
#ifdef CLAB_MAIN

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>

#include "crimlab.h"
#include "tk.h"

namespace {

// SK4FGA 0.1.1, R 4.5.0, dev/doss/SK4FGA.R: partition(ri)$groups for the 24
// values of ri() below (its fixed variance 1.6e-9, alpha 0.05)
const int ref_ri[24] = {2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1};
// partition((ri - 1.5) * 1000)$groups: 24 singletons
const int ref_resc[24] = {8, 11, 23, 2, 9, 12, 24, 10, 15, 13, 19, 17, 16, 14, 20, 18, 5, 21, 3, 6, 1, 22, 4, 7};

void ri(double *x) {
  const double mu[3] = {1.5180, 1.5186, 1.5200};
  for (int j = 0; j < 24; ++j) x[j] = mu[j / 8] + ((j * 7) % 9 - 4) * 1e-5;
}

struct rng {
  std::uint64_t s;
  double unif() {
    s ^= s >> 12;
    s ^= s << 25;
    s ^= s >> 27;
    return static_cast<double>((s * 2685821657736338717ULL) >> 11) * 0x1p-53;
  }
  double norm() {
    const double a = 1.0 - unif(), b = unif();
    return std::sqrt(-2.0 * std::log(a)) * std::cos(6.283185307179586 * b);
  }
};

// the same partition: i, j together in a iff together in b
bool same(const int *a, const int *b, int n) {
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      if ((a[i] == a[j]) != (b[i] == b[j])) return false;
  return true;
}

void t_sk() {
  static int iw[3 * 400], g[400], g2[400];
  static double w[2 * 400], x[400];
  int ng = 0;
  // 1. parity with SK4FGA on refractive indices at its own variance
  ri(x);
  TK("sk.ri.rc", clab_skp(x, 24, 1.6e-9, 0.05, 999, 1, g, &ng, iw, w) == CLAB_OK);
  TK("sk.ri.ng", ng == 3);
  TK("sk.ri.same", same(g, ref_ri, 24));
  for (int j = 0; j < 24; ++j) TK("sk.ri.order", g[j] == 1 + j / 8);
  // 2. disagreement (dossier csls/d/SK4FGA.md section 3, reproduced in
  // dev/doss/SK4FGA.R): rescaled data with the rescaled variance give the
  // same three groups; SK4FGA 0.1.1 keeps 1.6e-9 and returns 24 singletons
  for (int j = 0; j < 24; ++j) x[j] = (x[j] - 1.5) * 1000;
  clab_skp(x, 24, 1.6e-9 * 1e6, 0.05, 999, 1, g2, &ng, iw, w);
  TK("sk.resc", ng == 3 && same(g, g2, 24));
  TK("sk.resc.field", !same(g2, ref_resc, 24));
  // 3. homogeneous sample: one group
  for (int j = 0; j < 20; ++j) x[j] = 1.5183 + ((j * 7) % 9 - 4) * 1e-5;
  clab_skp(x, 20, 1.6e-9, 0.05, 999, 1, g, &ng, iw, w);
  TK("sk.homog", ng == 1);
  // 4. k = 2 against the exact null: lambda = c (x1 - x2)^2 / (2 s2) and
  // (x1 - x2)^2 / (2 s2) ~ chi-square(1), so p = erfc(|x1 - x2| / (2 sd)).
  // d = 3.6428 sd gives p = 0.0100 (split at 0.05); d = 1.8124 sd gives
  // p = 0.2000 (no split). With 999 simulations the Monte Carlo p is
  // within 0.013 (4 sd) of either.
  {
    TKN("sk.k2.p01", std::erfc(3.6428 / 2), 0.0100, 1e-4);
    TKN("sk.k2.p20", std::erfc(1.8124 / 2), 0.2000, 1e-4);
    double y[2] = {0, 3.6428};
    clab_skp(y, 2, 1.0, 0.05, 999, 3, g, &ng, iw, w);
    TK("sk.k2.split", ng == 2);
    y[1] = 1.8124;
    clab_skp(y, 2, 1.0, 0.05, 999, 3, g, &ng, iw, w);
    TK("sk.k2.nosplit", ng == 1);
  }
  // 5. invariance: a + b x with variance b^2 s2 gives the same groups
  {
    rng r{4};
    for (int j = 0; j < 60; ++j) x[j] = (j % 3) * 10.0 + r.norm();
    clab_skp(x, 60, 1.0, 0.05, 499, 7, g, &ng, iw, w);
    const int ng1 = ng;
    for (int j = 0; j < 60; ++j) x[j] = -3.0 + 0.25 * x[j];
    clab_skp(x, 60, 0.0625, 0.05, 499, 7, g2, &ng, iw, w);
    TK("sk.inv", ng == ng1 && same(g, g2, 60));
    // clusters 10 sd apart are never merged (an extreme value may still be
    // split off inside one: with this seed, one value of the middle cluster)
    bool pure = ng1 >= 3;
    for (int i = 0; i < 60; ++i)
      for (int j = 0; j < 60; ++j) pure = pure && (g[i] != g[j] || i % 3 == j % 3);
    TK("sk.inv.pure", pure);
  }
  // 6. size of the first split under the null: 400 samples of 50 from
  // N(0, 1) with s2 = 1, 199 simulations each; the rate at 0.05 has sd
  // 0.011, tolerance 4 sd
  {
    rng r{9};
    int rej = 0;
    for (int rep = 0; rep < 400; ++rep) {
      for (int j = 0; j < 50; ++j) x[j] = r.norm();
      clab_skp(x, 50, 1.0, 0.05, 199, static_cast<unsigned>(rep), g, &ng, iw, w);
      rej += ng > 1;
    }
    TK("sk.size", rej >= 0.006 * 400 && rej <= 0.094 * 400);
  }
  // 7. guards: NaN is refused (no silent drop that would misalign labels)
  {
    double y[3] = {1, std::nan(""), 2};
    TK("sk.guard.nan", clab_skp(y, 3, 1, 0.05, 9, 1, g, &ng, iw, w) == CLAB_EDOM);
    TK("sk.guard.s2", clab_skp(y, 1, 0, 0.05, 9, 1, g, &ng, iw, w) == CLAB_EDOM);
    TK("sk.guard.a", clab_skp(y, 1, 1, 1.0, 9, 1, g, &ng, iw, w) == CLAB_EDOM);
    TK("sk.guard.n", clab_skp(y, 0, 1, 0.05, 9, 1, g, &ng, iw, w) == CLAB_EARG);
    TK("sk.n1", clab_skp(y, 1, 1, 0.05, 9, 1, g, &ng, iw, w) == CLAB_OK && ng == 1 && g[0] == 1);
  }
}

// the dossier's standard case: n = 5,000 values N(1.518, 1e-4^2) at the
// fixed variance 1.6e-9 (SK4FGA::partition 0.294 s median of 3, 6 groups)
void b_sk() {
  const int n = 5000;
  static double x[n], w[2 * n];
  static int g[n], iw[3 * n];
  rng r{2};
  for (int i = 0; i < n; ++i) x[i] = 1.518 + 1e-4 * r.norm();
  double t[5];
  int ng = 0;
  for (int ns = 199; ns <= 999; ns += 800) {
    for (int k = 0; k < 5; ++k) {
      const auto a = std::chrono::steady_clock::now();
      clab_skp(x, n, 1.6e-9, 0.05, ns, 1, g, &ng, iw, w);
      t[k] = std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count();
    }
    std::printf("BENCH clab_skp n=5000 nsim=%d: median %.6f s (k=5), %d groups\n", ns, tk::med5(t), ng);
  }
}

}  // namespace

TK_FAMILY("sk", t_sk, b_sk);

#endif
