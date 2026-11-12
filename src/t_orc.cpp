// SPDX-License-Identifier: AGPL-3.0-or-later
//
// Parity of crimlab's existing forensic kernels (src/core.cpp: clab_gtmp,
// clab_kinidx, clab_patidx, clab_mixp, clab_cllr, clab_lrcal) with
// independent implementations in other packages. Built only with
// -DCLAB_MAIN.
// tk-deps: core fair
#ifdef CLAB_MAIN

#include <cmath>
#include <cstdio>

#include "crimlab.h"
#include "tk.h"

namespace {

// Reference values, dev/doss/oracles.R, R 4.5.0. One locus, alleles A..F
// with frequencies FR, no mutation.
// ref_pi: Familias 2.6.4, FamiliasPosterior(list(af is father, unrelated),
//   ...)$likelihoods ratio for the trios TRIO below (mother, child, AF);
// ref_ki: Familias 2.6.4, likelihood ratios of full sibs, half sibs and
//   parent-child against unrelated for the pairs PAIR below (column-major:
//   the 7 full-sib values, then half sibs, then parent-child);
// ref_mp: Familias 2.6.4, P(G, G | two unrelated, kinship theta) /
//   P(G | theta) for the genotypes GT below (fastest) and theta = 0, 0.01,
//   0.03, 0.1;
// ref_nd1..3: numberofalleles 1.0.1, pr_number_of_distinct_alleles(2 x,
//   FR) for x = 1, 2, 3 contributors, counts 1..6 (NA written as 0);
// ref_dna1..3: DNAtools 0.2.5, Pnm_locus(m = x, theta = 0, alleleProbs =
//   FR), the same distribution by its compiled routine;
// ref_fm, ref_fm_low: fairmodels 1.2.2, fairness_check(explain(fixed
//   scores p), protected = g, privileged = "a", cutoff = 0.5 or 0.001)
//   $groups_data$s: ACC TPR TNR FPR FNR PPV NPV STP for groups a b c, on
//   the deterministic sample of dev/doss/fairness.R;
// ref_ece: comparison 1.0.8, calc.ece(lss, lds) at prior 0.5: $ece and
//   $ece.cal, with lss = exp(seq(-1, 3, length.out = 15)) and
//   lds = exp(seq(-4, 1.2, length.out = 18)).
const double ref_pi[8] = {3.3333333333333335, 6.6666666666666661, 6.6666666666666661, 2.5, 6.666666666666667, 0, 2.5, 4.9999999999999991};
const double ref_ki[21] = {28.999999999999993, 5.2500000000000009, 0.25, 1.9166666666666665, 7.6880165289256199, 6.4107142857142856, 1.0833333333333335, 4.2499999999999991, 5.4999999999999991, 0.50000000000000011, 2.1666666666666665, 2.7727272727272725, 2.1964285714285716, 1.3333333333333333, 7.4999999999999991, 9.9999999999999982, 0, 3.3333333333333335, 4.545454545454545, 3.3928571428571428, 1.6666666666666667};
const double ref_mp[20] = {0.0025000000000000005, 0.010000000000000002, 0.084000000000000005, 0.078400000000000011, 0.088000000000000009, 0.0053632789749563196, 0.012590759075907591, 0.088373519704911668, 0.088623412929528272, 0.09198679867986799, 0.013763738779996336, 0.018262502289796666, 0.096960615497343805, 0.10982465653050007, 0.099874702326433418, 0.064034090909090915, 0.041742424242424248, 0.12533333333333338, 0.18901818181818186, 0.12642424242424244};
const double ref_nd1[6] = {0.20180000000000003, 0.79820000000000002, 0, 0, 0, 0};
const double ref_nd2[6] = {0.010701620000000002, 0.22765838000000005, 0.55368000000000006, 0.20795999999999992, 0, 0};
const double ref_nd3[6] = {0.0006716764580000003, 0.047767833542000017, 0.33173658000000017, 0.46752951000000009, 0.14564159999999993, 0.0066527999999999449};
const double ref_ece[2] = {0.63034589985460643, 0.48597032561095654};
const double ref_dna1[2] = {0.20180000000000003, 0.79820000000000002};
const double ref_dna2[4] = {0.010701620000000002, 0.22765838000000005, 0.55368000000000006, 0.20795999999999992};
const double ref_dna3[6] = {0.0006716764580000003, 0.047767833542000017, 0.33173658000000017, 0.46752951000000009, 0.14564159999999993, 0.0066527999999999726};
const double ref_fm[24] = {0.64071856287425155, 0.65765765765765771, 0.64564564564564564, 0.63888888888888884, 0.65492957746478875, 0.63380281690140849, 0.64210526315789473, 0.65968586387434558, 0.65445026178010468, 0.35789473684210527, 0.34031413612565448, 0.34554973821989526, 0.3611111111111111, 0.34507042253521125, 0.36619718309859156, 0.57499999999999996, 0.58860759493670889, 0.57692307692307687, 0.70114942528735635, 0.71999999999999997, 0.70621468926553677, 0.47904191616766467, 0.47447447447447449, 0.46846846846846846};
const double ref_fm_low[24] = {0.43712574850299402, 0.43243243243243246, 0.42942942942942941, 1, 1, 1, 0.010526315789473684, 0.010471204188481676, 0.005235602094240838, 0.98947368421052628, 0.98952879581151831, 0.99476439790575921, 0, 0, 0, 0.43373493975903615, 0.42900302114803623, 0.42771084337349397, 1, 1, 1, 0.99401197604790414, 0.99399399399399402, 0.99699699699699695};

const double FR[6] = {0.05, 0.1, 0.15, 0.2, 0.22, 0.28};
// trios (mother, child, alleged father), alleles A..F as 0..5
const int TRIO[8][6] = {{0, 1, 0, 2, 2, 3}, {0, 1, 0, 1, 1, 1}, {0, 1, 0, 1, 0, 1}, {2, 2, 2, 3, 3, 4},
                        {0, 1, 1, 2, 2, 2}, {4, 5, 4, 5, 0, 1}, {3, 3, 3, 3, 3, 5}, {1, 2, 1, 3, 3, 3}};
const int PAIR[7][4] = {{0, 1, 0, 1}, {0, 0, 0, 1}, {0, 1, 2, 3}, {2, 3, 2, 2},
                        {4, 4, 4, 4}, {1, 5, 5, 1}, {0, 2, 2, 4}};
const int GT[5][2] = {{0, 0}, {0, 1}, {2, 5}, {5, 5}, {3, 4}};

void t_orc() {
  // 1. paternity index against Familias
  for (int k = 0; k < 8; ++k) {
    const int *t = TRIO[k];
    const int c[2] = {t[2], t[3]}, m[2] = {t[0], t[1]}, af[2] = {t[4], t[5]};
    const double f[2] = {FR[c[0]], FR[c[1]]};
    double idx;
    TK("orc.pi.rc", clab_patidx(c, m, af, f, 1, &idx) == CLAB_OK);
    TKN("orc.pi.familias", idx, ref_pi[k], 1e-13);
  }
  // 2. kinship index (full sibs, half sibs, parent-child) against Familias
  {
    const double kap[3][3] = {{0.25, 0.5, 0.25}, {0.5, 0.5, 0}, {0, 1, 0}};
    for (int r = 0; r < 3; ++r)
      for (int k = 0; k < 7; ++k) {
        const int g1[2] = {PAIR[k][0], PAIR[k][1]}, g2[2] = {PAIR[k][2], PAIR[k][3]};
        const double f[2] = {FR[g1[0]], FR[g1[1]]};
        double idx;
        TK("orc.ki.rc", clab_kinidx(g1, g2, f, 1, kap[r], &idx) == CLAB_OK);
        TKN("orc.ki.familias", idx, ref_ki[7 * r + k], 1e-13);
      }
  }
  // 3. match probability with the coancestry correction (NRC II 4.10)
  // against Familias' joint over single likelihood under kinship theta
  {
    const double th[4] = {0, 0.01, 0.03, 0.1};
    for (int j = 0; j < 4; ++j)
      for (int k = 0; k < 5; ++k) {
        const double p = FR[GT[k][0]], q = FR[GT[k][1]];
        const int hom = GT[k][0] == GT[k][1];
        double loc, l10;
        TK("orc.mp.rc", clab_gtmp(&p, &q, &hom, 1, th[j], &loc, &l10) == CLAB_OK);
        TKN("orc.mp.familias", loc, ref_mp[5 * j + k], 1e-13);
      }
  }
  // 4. mixture probabilities: summing clab_mixp (the chance that x unknown
  // contributors carry exactly the alleles of S) over all sets S of k of
  // the 6 alleles gives the chance of k distinct alleles, which
  // numberofalleles and DNAtools compute by their own routines
  {
    const double *ref[3] = {ref_nd1, ref_nd2, ref_nd3};
    for (int x = 1; x <= 3; ++x) {
      double pk[7] = {0, 0, 0, 0, 0, 0, 0};
      for (int s = 1; s < 64; ++s) {
        double p[6];
        int need[6], n = 0;
        for (int j = 0; j < 6; ++j)
          if ((s >> j) & 1) {
            p[n] = FR[j];
            need[n] = 1;
            ++n;
          }
        double v;
        clab_mixp(p, need, n, x, &v);
        pk[n] += v;
      }
      const double *dna = x == 1 ? ref_dna1 : (x == 2 ? ref_dna2 : ref_dna3);
      for (int k = 1; k <= 6; ++k) {
        TKN("orc.nd", pk[k], ref[x - 1][k - 1], 1e-13);
        if (k <= 2 * x) TKN("orc.dnatools", pk[k], dna[k - 1], 1e-13);
      }
    }
  }
  // 5. Cllr and the minimum Cllr against the empirical cross-entropy at
  // prior 0.5 (base 2), raw and after pool-adjacent-violators calibration
  {
    double ss[15], ds[18], lr[33], post[33], work[66], c, cmin;
    int same[33];
    for (int i = 0; i < 15; ++i) ss[i] = std::exp(-1 + 4.0 * i / 14);
    for (int i = 0; i < 18; ++i) ds[i] = std::exp(-4 + 5.2 * i / 17);
    TK("orc.cllr.rc", clab_cllr(ss, 15, ds, 18, &c) == CLAB_OK);
    TKN("orc.cllr.ece", c, ref_ece[0], 1e-13);
    // merge in ascending order
    int a = 0, b = 0;
    for (int k = 0; k < 33; ++k) {
      const bool takes = b >= 18 || (a < 15 && ss[a] < ds[b]);
      lr[k] = takes ? ss[a++] : ds[b++];
      same[k] = takes;
    }
    TK("orc.lrcal.rc", clab_lrcal(lr, same, 33, post, work, &cmin) == CLAB_OK);
    std::printf("  orc: cllrmin %.17g, comparison ece.cal %.17g\n", cmin, ref_ece[1]);
    TKN("orc.lrcal.ece", cmin, ref_ece[1], 1e-12);
  }
  // 6. clab_fair against fairmodels (a second oracle for it, one that is
  // right where a nearly all-positive cutoff breaks the other)
  {
    static double p[1000];
    static int y[1000], g[1000];
    for (int i = 0; i < 1000; ++i) {
      y[i] = ((i * 13) % 7) < 3 ? 1 : 0;
      g[i] = i % 3;
      p[i] = (((i * 37) % 101) + 30.0 * y[i]) / 131.0;
    }
    double cnt[15], met[30], rat[30];
    for (int c = 0; c < 2; ++c) {
      const double *ref = c == 0 ? ref_fm : ref_fm_low;
      clab_fair(p, y, g, 1000, 3, c == 0 ? 0.5 : 0.001, 0, 1.96, cnt, met, rat, nullptr, nullptr);
      for (int k = 0; k < 8; ++k)
        for (int j = 0; j < 3; ++j) TKN("orc.fairmodels", met[3 * k + j], ref[3 * k + j], 1e-14);
    }
  }
}

}  // namespace

TK_FAMILY("orc", t_orc, nullptr);

#endif
