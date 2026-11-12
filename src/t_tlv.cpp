// SPDX-License-Identifier: AGPL-3.0-or-later
//
// Standalone tests for clab_tlvc and clab_tlvlr (src/tlv.cpp). Part of the
// self-test (selftest.cpp).

#include <cmath>
#include <cstdint>

#include "crimlab.h"
#include "tchk.h"

namespace {

// Reference values from dev/doss/comparison.R (comparison 1.0.8, R 4.5.0).
// gU, gC, gmu: two.level.components(glass, c(7, 8, 9), 1) (glass K, Ca, Fe);
// gctl, gr1, gr2: two.level.comparison.items of glass rows 1:6, 7:12 and
// 19:24 (6 replicates each); glr: two.level.normal.LR(ctl, r1, Z) and
// (ctl, r2, Z). gfar: natural-log LR of (gs40, gs40) and (gp20, gm20) from
// the 2p-variate joint density in base R (comparison returns NaN for both).
// bU, bC, bmu: two.level.components on the balanced MINSTD sample (30 x 4);
// uU, uC: the same on the unbalanced sample (n_i = 2, 3, 5, 8 repeating);
// sU, sC: the pooled-within and ANOVA estimates of that sample in base R.
const double gU[9] = {0.39896219988105996, 0.0030193434087121131, -0.027026314339015214, 0.0030193434087121131, 0.024478855009848442, 0.0026876797352272763, -0.027026314339015214, 0.0026876797352272763, 0.21975470811590886};
const double gC[9] = {2.2721186922751548, -0.5024275724461893, -0.28341571893493966, -0.5024275724461893, 0.85597421139504692, 0.16968221937235764, -0.28341571893493966, 0.16968221937235764, 1.8734637411683064};
const double gmu[3] = {-3.1066036666666665, -1.13828, -4.7975165833333335};
const double gctl[3] = {-1.7326999999999999, -1.1209666666666667, -5.1186666666666669};
const double gr1[3] = {-1.69025, -1.0625166666666666, -4.6171166666666661};
const double gr2[3] = {-3.1130833333333334, -1.1454833333333334, -5.5865333333333336};
const double glr[2] = {51.142431056211358, 0.0289990757442239};
const double gfar[2] = {3262.3198454486483, -114466.36982375706};
const double gs40[3] = {58.561494642935962, 35.886583079730691, 49.631145989353357};
const double gp20[3] = {28.41439732146798, 17.382808206532012, 22.256239661343347};
const double gm20[3] = {-31.879797321467983, -19.624741539865347, -32.493572994676683};
const double bU[9] = {0.00087950588744766748, -3.4264752085989355e-05, -0.00016449772197268156, -3.4264752085989355e-05, 0.00084538221210306861, -7.406400636017366e-05, -0.00016449772197268156, -7.406400636017366e-05, 0.00095989403454484516};
const double bC[9] = {0.23701590060677694, 0.048684372364825597, 0.16627847958109138, 0.048684372364825597, 0.23827424115560622, 0.16113471709988367, 0.16627847958109138, 0.16113471709988367, 0.42693765808822753};
const double bmu[3] = {0.53802389427617059, 0.622597592264304, 0.75781491213259977};
const double uU[9] = {0.001205862433681829, -9.382354722479472e-05, -0.00018144133802358169, -9.382354722479472e-05, 0.0010080883770283437, 1.6153555542244327e-05, -0.00018144133802358169, 1.6153555542244327e-05, 0.0011968601462841091};
const double uC[9] = {0.2376928737417188, 0.07353947616016944, 0.18937107453608676, 0.07353947616016944, 0.2129592942245247, 0.13237864895464579, 0.18937107453608676, 0.13237864895464579, 0.43532267506947486};
const double sU[9] = {0.00089448904616081534, -4.1955214802495819e-05, -8.0635826245155956e-05, -4.1955214802495819e-05, 0.00078905299226618138, -8.2856124695310251e-05, -8.0635826245155956e-05, -8.2856124695310251e-05, 0.00088611526329112884};
const double sC[9] = {0.24000621071588024, 0.074220936528557074, 0.19113347225086691, 0.074220936528557074, 0.21501805699611853, 0.13364981943910245, 0.19113347225086691, 0.13364981943910245, 0.43949944418779197};

struct xrng {
  std::uint64_t s;
  double unif() {
    s ^= s >> 12;
    s ^= s << 25;
    s ^= s >> 27;
    return static_cast<double>((s * 2685821657736338717ULL) >> 11) * 0x1p-53;
  }
  double norm() {  // Box-Muller, one draw per call
    const double a = 1.0 - unif(), b = unif();
    return std::sqrt(-2.0 * std::log(a)) * std::cos(6.283185307179586 * b);
  }
};

// the sample of dev/doss/comparison.R: MINSTD from seed 12345; m item
// effects (column-major m x 3, element e scaled by {1, 2, 0.5}[e % 3]), then
// nrep[i] rows per item of effect + 0.1 (u - 0.5). x is n x 3 column-major.
int mksample(int m, const int *nrep, double *x, int *item) {
  std::int64_t s = 12345;
  auto u = [&s]() {
    s = (48271 * s) % 2147483647;
    return static_cast<double>(s) / 2147483647.0;
  };
  double eff[3 * 64];
  const double mult[3] = {1, 2, 0.5};
  for (int e = 0; e < 3 * m; ++e) eff[e] = u() * mult[e % 3];
  int n = 0;
  for (int i = 0; i < m; ++i) n += nrep[i];
  int r = 0;
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < nrep[i]; ++j, ++r) {
      for (int k = 0; k < 3; ++k) x[r + n * k] = eff[i + m * k] + 0.1 * (u() - 0.5);
      item[r] = i;
    }
  }
  return n;
}

// log density of the 2p-variate normal of (a, b) given the model: the
// covariance [[U/na + C, sC], [sC, U/nb + C]] with s = 1 (same source) or
// 0; a dense Cholesky of the 2p x 2p matrix, a route independent of the
// difference / weighted-mean factorisation the kernel uses.
double ljoint(const double *u, const double *c, const double *mu, int p,
              const double *a, double na, const double *b, double nb, int same) {
  double s[36], r[6], v[6];
  const int q = 2 * p;
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j < p; ++j) {
      s[i + q * j] = u[i + p * j] / na + c[i + p * j];
      s[i + p + q * (j + p)] = u[i + p * j] / nb + c[i + p * j];
      s[i + p + q * j] = s[i + q * (j + p)] = same ? c[i + p * j] : 0.0;
    }
    r[i] = a[i] - mu[i];
    r[i + p] = b[i] - mu[i];
  }
  double ld = 0.0, qf = 0.0;
  for (int j = 0; j < q; ++j) {
    double d = s[j + q * j];
    for (int k = 0; k < j; ++k) d -= s[j + q * k] * s[j + q * k];
    d = std::sqrt(d);
    s[j + q * j] = d;
    ld += 2.0 * std::log(d);
    for (int i = j + 1; i < q; ++i) {
      double e = s[i + q * j];
      for (int k = 0; k < j; ++k) e -= s[i + q * k] * s[j + q * k];
      s[i + q * j] = e / d;
    }
  }
  for (int i = 0; i < q; ++i) {
    double e = r[i];
    for (int k = 0; k < i; ++k) e -= s[i + q * k] * v[k];
    v[i] = e / s[i + q * i];
    qf += v[i] * v[i];
  }
  return -0.5 * (ld + qf);  // the (2 pi)^p factors cancel in a ratio
}

void run(tk &tt) {
  double out[2], work[9 + 6];
  // 1. parity with comparison 1.0.8 on its glass example
  tk_true(tt, "tlv.glass.rc", clab_tlvlr(gU, gC, gmu, 3, gctl, 6, gr1, 6, out, work) == CLAB_OK);
  tk_near(tt, "tlv.glass.same", out[0], std::log(glr[0]), 1e-12);
  tk_near(tt, "tlv.glass.log10", out[1], std::log10(glr[0]), 1e-12);
  clab_tlvlr(gU, gC, gmu, 3, gctl, 6, gr2, 6, out, work);
  tk_near(tt, "tlv.glass.diff", out[0], std::log(glr[1]), 1e-12);
  // 2. the joint-density route on the same inputs
  tk_near(tt, "tlv.joint.same", out[0],
      ljoint(gU, gC, gmu, 3, gctl, 6, gr2, 6, 1) - ljoint(gU, gC, gmu, 3, gctl, 6, gr2, 6, 0),
      1e-12);
  // 3. items 40 between-item SDs from the mean (same source) and 40 SDs
  // apart: the log LR is finite and equals the joint-density value from
  // base R (dev/doss/comparison.R).
  clab_tlvlr(gU, gC, gmu, 3, gs40, 6, gs40, 6, out, work);
  tk_near(tt, "tlv.far.same", out[0], gfar[0], 1e-12);
  tk_near(tt, "tlv.far.same.c", out[0],
      ljoint(gU, gC, gmu, 3, gs40, 6, gs40, 6, 1) - ljoint(gU, gC, gmu, 3, gs40, 6, gs40, 6, 0),
      1e-12);
  clab_tlvlr(gU, gC, gmu, 3, gp20, 6, gm20, 6, out, work);
  tk_near(tt, "tlv.far.apart", out[0], gfar[1], 1e-12);
  // 4. identities: a common shift of the items and the mean leaves the LR
  // unchanged; exchanging the items (with their counts) too
  {
    double a[3], b[3], m2[3], o2[2];
    for (int k = 0; k < 3; ++k) {
      a[k] = gctl[k] + 0.75;
      b[k] = gr1[k] + 0.75;
      m2[k] = gmu[k] + 0.75;
    }
    clab_tlvlr(gU, gC, gmu, 3, gctl, 3, gr1, 7, out, work);
    clab_tlvlr(gU, gC, m2, 3, a, 3, b, 7, o2, work);
    tk_near(tt, "tlv.shift", o2[0], out[0], 1e-12);
    clab_tlvlr(gU, gC, gmu, 3, gr1, 7, gctl, 3, o2, work);
    tk_near(tt, "tlv.swap", o2[0], out[0], 1e-13);
  }
  // 5. univariate closed form: U = 1, C = 4, mu = 0, n = 1 and 1, x = 1, 3:
  // ln LR = -1/2 [4/2 + ln 2 + 4/4.5 + ln 4.5] + 1/2 [1/5 + 9/5 + 2 ln 5]
  {
    const double u = 1, c = 4, mu = 0, xa = 1, xb = 3;
    clab_tlvlr(&u, &c, &mu, 1, &xa, 1, &xb, 1, out, work);
    const double want = -0.5 * (2.0 + std::log(2.0) + 4.0 / 4.5 + std::log(4.5)) +
                        0.5 * (2.0 + 2.0 * std::log(5.0));
    tk_near(tt, "tlv.uni", out[0], want, 1e-14);
  }
  // 6. components: parity on the balanced sample, the ANOVA estimates on
  // the unbalanced one
  {
    static double x[200 * 3], w[30 * 3 + 60], mu[3], u[9], c[9];
    static int item[200], nrep[30];
    for (int i = 0; i < 30; ++i) nrep[i] = 4;
    int n = mksample(30, nrep, x, item);
    tk_true(tt, "tlv.comp.rc", clab_tlvc(x, n, 3, item, 30, mu, u, c, w) == CLAB_OK);
    for (int k = 0; k < 9; ++k) {
      tk_near(tt, "tlv.comp.C", c[k], bC[k], 1e-12);
    }
    for (int k = 0; k < 3; ++k) tk_near(tt, "tlv.comp.mu", mu[k], bmu[k], 1e-14);
    // U is of order 1e-3: compared relative to its own entries
    for (int k = 0; k < 9; ++k) tk_near(tt, "tlv.comp.U", u[k] / bU[k], 1, 1e-11);
    const int pat[4] = {2, 3, 5, 8};
    for (int i = 0; i < 30; ++i) nrep[i] = pat[i % 4];
    n = mksample(30, nrep, x, item);
    clab_tlvc(x, n, 3, item, 30, mu, u, c, w);
    for (int k = 0; k < 9; ++k) {
      tk_near(tt, "tlv.unbal.U", u[k] / sU[k], 1, 1e-11);
      tk_near(tt, "tlv.unbal.C", c[k], sC[k], 1e-12);
    }
    // the size-weighted within estimate (uU, uC) is a different estimator:
    // they must not coincide on unbalanced data
    tk_true(tt, "tlv.unbal.differs", std::fabs(uU[0] / u[0] - 1) > 0.2);
    // its between-item C differs too, through U (about 1% here)
    tk_true(tt, "tlv.unbal.C.differs", std::fabs(uC[0] / c[0] - 1) > 0.005);
  }
  // 7. seeded simulation: unbiasedness on the unbalanced design, true
  // U = I, C = 4 I, 30 items, 400 data sets. One U_kk estimate has
  // variance 2 / (N - m) = 2 / 105, so the mean of 400 has sd 0.0069;
  // C_kk ~ (S_B - (m-1) U) / (N - sum n_i^2 / N) has sd about 1.1 per set,
  // 0.055 for the mean. Tolerances are 4 sd: 0.03 and 0.25.
  {
    static double x[200 * 3], w[30 * 3 + 60];
    static int item[200];
    int nrep[30], n = 0;
    const int pat[4] = {2, 3, 5, 8};
    for (int i = 0; i < 30; ++i) n += nrep[i] = pat[i % 4];
    xrng r{2026};
    double su[3] = {0, 0, 0}, sc[3] = {0, 0, 0};
    for (int rep = 0; rep < 400; ++rep) {
      int row = 0;
      for (int i = 0; i < 30; ++i) {
        double th[3];
        for (int k = 0; k < 3; ++k) th[k] = 2.0 * r.norm();
        for (int j = 0; j < nrep[i]; ++j, ++row) {
          for (int k = 0; k < 3; ++k) x[row + n * k] = th[k] + r.norm();
          item[row] = i;
        }
      }
      double mu[3], u[9], c[9];
      clab_tlvc(x, n, 3, item, 30, mu, u, c, w);
      for (int k = 0; k < 3; ++k) {
        su[k] += u[4 * k] / 400;
        sc[k] += c[4 * k] / 400;
      }
    }
    for (int k = 0; k < 3; ++k) {
      tk_near(tt, "tlv.sim.U", su[k], 1.0, 0.03);
      tk_true(tt, "tlv.sim.C", std::fabs(sc[k] - 4.0) < 0.25);
    }
  }
  // 8. guards and a non-positive-definite C
  {
    const double u = 1, c = -2, mu = 0, xa = 0, xb = 0;
    tk_true(tt, "tlv.guard.pd", clab_tlvlr(&u, &c, &mu, 1, &xa, 1, &xb, 1, out, work) == CLAB_EDOM);
    const double c2 = 1, um = -1, nanx = std::nan("");
    tk_true(tt, "tlv.guard.nan", clab_tlvlr(&u, &c2, &mu, 1, &nanx, 1, &xb, 1, out, work) == CLAB_EDOM);
    tk_true(tt, "tlv.guard.U", clab_tlvlr(&um, &c2, &mu, 1, &xa, 1, &xb, 1, out, work) == CLAB_EDOM);
    tk_true(tt, "tlv.guard.n", clab_tlvlr(&u, &c2, &mu, 1, &xa, 0, &xb, 1, out, work) == CLAB_EDOM);
    tk_true(tt, "tlv.guard.p", clab_tlvlr(&u, &c2, &mu, 0, &xa, 1, &xb, 1, out, work) == CLAB_EARG);
    double x[3] = {1, 2, 3}, w[10], m3[1], u3[1], c3[1];
    int it[3] = {0, 1, 1}, it2[3] = {0, 2, 2};
    tk_true(tt, "tlv.guard.nm", clab_tlvc(x, 2, 1, it, 2, m3, u3, c3, w) == CLAB_EARG);
    tk_true(tt, "tlv.guard.empty", clab_tlvc(x, 3, 1, it2, 3, m3, u3, c3, w) == CLAB_EARG);
    x[1] = std::nan("");
    tk_true(tt, "tlv.guard.nan", clab_tlvc(x, 3, 1, it, 2, m3, u3, c3, w) == CLAB_EDOM);
  }
}


}  // namespace

void t_tlv(tk &t, int) { run(t); }

#ifdef CLAB_MAIN
#include <chrono>
#include <cstdio>

// timing: components of 200 items x 12 replicates x 3
// variables (comparison::two.level.components on glass took 2.423 s,
// median of 5 on zeus), and one LR (0.00019 s per call)
void b_tlv() {
  static double x[2400 * 3], w[200 * 3 + 400];
  static int item[2400];
  xrng r{99};
  for (int i = 0; i < 2400; ++i) {
    item[i] = i / 12;
    for (int k = 0; k < 3; ++k) x[i + 2400 * k] = r.norm() + 3.0 * (item[i] % 7);
  }
  double mu[3], u[9], c[9], t[5], t2[5], out[2], wk[15];
  for (int k = 0; k < 5; ++k) {
    auto a = std::chrono::steady_clock::now();
    clab_tlvc(x, 2400, 3, item, 200, mu, u, c, w);
    t[k] = std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count();
    a = std::chrono::steady_clock::now();
    for (int j = 0; j < 1000; ++j) clab_tlvlr(gU, gC, gmu, 3, gctl, 6, gr1, 6, out, wk);
    t2[k] = std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count() / 1000;
  }
  std::printf("BENCH clab_tlvc 200 x 12 x 3: median %.6f s (k=5)\n", tk_med5(t));
  std::printf("BENCH clab_tlvlr p=3: %.9f s per call (median of 5 x 1000)\n", tk_med5(t2));
}
#endif
