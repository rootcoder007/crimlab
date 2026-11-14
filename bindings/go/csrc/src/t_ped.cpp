// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Self-test of clab_mutmat and clab_pedlik (src/ped.cpp): likelihoods by
// hand, and against brute force over every allele each person may carry
// (a pedigree with a loop among them), with and without mutation.
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "tchk.h"

namespace {

// the likelihood by summing over every ordered allele pair of every
// person: founders from p, children from a parent's allele chosen with
// probability 1/2 then mutated (m: NULL for none); typed persons must
// carry their unordered genotype
double brute(int n, const int *fa, const int *mo, const int *g, int na, const double *p,
             const double *m) {
  const int nv = 2 * n;
  std::vector<int> a(static_cast<std::size_t>(nv), 0);
  double tot = 0.0;
  for (;;) {
    double pr = 1.0;
    for (int i = 0; i < n && pr > 0; ++i) {
      const int x = a[static_cast<std::size_t>(2 * i)], y = a[static_cast<std::size_t>(2 * i + 1)];
      if (g[2 * i] >= 0 && !((x == g[2 * i] && y == g[2 * i + 1]) || (y == g[2 * i] && x == g[2 * i + 1]))) {
        pr = 0;
        break;
      }
      if (fa[i] < 0) {
        pr *= p[x] * p[y];
        continue;
      }
      for (int s = 0; s < 2; ++s) {
        const int par = s == 0 ? fa[i] : mo[i], c = s == 0 ? x : y;
        const int u = a[static_cast<std::size_t>(2 * par)], w = a[static_cast<std::size_t>(2 * par + 1)];
        const double tu = m ? m[u * na + c] : (u == c), tw = m ? m[w * na + c] : (w == c);
        pr *= 0.5 * (tu + tw);
      }
    }
    tot += pr;
    int k = nv - 1;
    while (k >= 0 && ++a[static_cast<std::size_t>(k)] == na) a[static_cast<std::size_t>(k--)] = 0;
    if (k < 0) break;
  }
  return tot;
}

void run(tk &t) {
  // stepwise matrix, alleles 8 9 10 12, rate 0.01, range 0.5 (the values
  // of pedmut 0.9.1's mutationMatrix("stepwise", rate2 = 0))
  const double lab[4] = {8, 9, 10, 12};
  double m[16];
  tk_ieq(t, "ped1", clab_mutmat(2, 4, lab, 0.01, 0.5, m), CLAB_OK);
  const double r0[4] = {0.99, 0.006153846153846154, 0.003076923076923077, 0.0007692307692307692};
  tk_veq(t, "ped2", m, r0, 4, 1e-15);
  tk_eq(t, "ped3", m[3 * 4 + 2], 0.01 * 0.5 / (0.5 + 0.25 + 0.125), 1e-15);
  tk_ieq(t, "ped4", clab_mutmat(1, 4, nullptr, 0.03, 0.5, m), CLAB_OK);
  tk_eq(t, "ped5", m[1], 0.01, 1e-15);
  tk_eq(t, "ped6", m[5], 0.97, 1e-15);
  tk_ieq(t, "ped7", clab_mutmat(2, 4, lab, 0.01, 1.5, m), CLAB_EDOM);
  // father F 0/1 and child 0/2, mother untyped: LR against unrelated is
  // 1 / (4 p0)
  const double p[3] = {0.1, 0.3, 0.6};
  const int fa[3] = {-1, -1, 0}, mo[3] = {-1, -1, 1};
  const int g[6] = {0, 1, -1, -1, 0, 2};
  double l1 = 0, l2 = 0;
  tk_ieq(t, "ped8", clab_pedlik(3, fa, mo, g, 3, p, nullptr, nullptr, &l1), CLAB_OK);
  const int fu[2] = {-1, -1}, mu[2] = {-1, -1}, gu[4] = {0, 1, 0, 2};
  clab_pedlik(2, fu, mu, gu, 3, p, nullptr, nullptr, &l2);
  tk_eq(t, "ped9", std::exp(l1 - l2), 1 / (4 * 0.1), 1e-13);
  tk_eq(t, "ped10", std::exp(l2), 2 * 0.1 * 0.3 * 2 * 0.1 * 0.6, 1e-13);
  // an exclusion: F 0/0, child 2/2. None: impossible; equal mutation at
  // rate r: the father's allele mutates to 2 (r / 2), the untyped mother
  // passes 2 (p2 (1 - r) + (1 - p2) r / 2), against p2^2 unrelated
  const int gx[6] = {0, 0, -1, -1, 2, 2}, gxu[4] = {0, 0, 2, 2};
  clab_pedlik(3, fa, mo, gx, 3, p, nullptr, nullptr, &l1);
  tk_true(t, "ped11", std::isinf(l1) && l1 < 0);
  double me[9];
  clab_mutmat(1, 3, nullptr, 0.002, 0.5, me);
  clab_pedlik(3, fa, mo, gx, 3, p, me, me, &l1);
  clab_pedlik(2, fu, mu, gxu, 3, p, me, me, &l2);
  tk_eq(t, "ped12", std::exp(l1 - l2), 0.001 * (0.6 * 0.998 + 0.4 * 0.001) / 0.36, 1e-12);
  // a child of half siblings (a loop), against brute force with and
  // without mutation: 0 (a man), 1 and 2 founders; 3, the son of 0 and 1,
  // and 4, the daughter of 0 and 2, half siblings; 5 a founder; 6 the
  // child of 3 and 4; 7 the child of 3 and 5; 1, 4, 6 and 7 typed
  const int lfa[8] = {-1, -1, -1, 0, 0, -1, 3, 3}, lmo[8] = {-1, -1, -1, 1, 2, -1, 4, 5};
  const int lg[16] = {-1, -1, 0, 1, -1, -1, -1, -1, 1, 1, -1, -1, 0, 1, 1, 1};
  const double p2[2] = {0.35, 0.65};
  double mb[4];
  clab_mutmat(1, 2, nullptr, 0.05, 0.5, mb);
  tk_ieq(t, "ped13", clab_pedlik(8, lfa, lmo, lg, 2, p2, nullptr, nullptr, &l1), CLAB_OK);
  tk_eq(t, "ped14", std::exp(l1), brute(8, lfa, lmo, lg, 2, p2, nullptr), 1e-12);
  clab_pedlik(8, lfa, lmo, lg, 2, p2, mb, mb, &l1);
  tk_eq(t, "ped15", std::exp(l1), brute(8, lfa, lmo, lg, 2, p2, mb), 1e-12);
  // full siblings, three alleles, stepwise mutation, against brute force
  const int sfa[4] = {-1, -1, 0, 0}, smo[4] = {-1, -1, 1, 1};
  const int sg[8] = {-1, -1, -1, -1, 0, 1, 1, 2};
  const double l3[3] = {14, 15, 16};
  double ms[9];
  clab_mutmat(2, 3, l3, 0.01, 0.3, ms);
  clab_pedlik(4, sfa, smo, sg, 3, p, ms, ms, &l1);
  tk_eq(t, "ped16", std::exp(l1), brute(4, sfa, smo, sg, 3, p, ms), 1e-12);
  // full siblings without mutation: the kinship index of siblings,
  // (1 + p0 / ... ) by the IBD formula 1/4 + 1/2 k1-term + 1/4 k2-term;
  // genotypes 0/1 and 1/2 share one allele: LR = 1/4 + 1/2 / (4 p1)
  clab_pedlik(4, sfa, smo, sg, 3, p, nullptr, nullptr, &l1);
  const int gs[4] = {0, 1, 1, 2};
  clab_pedlik(2, fu, mu, gs, 3, p, nullptr, nullptr, &l2);
  tk_eq(t, "ped17", std::exp(l1 - l2), 0.25 + 0.5 / (4 * 0.3), 1e-13);
  // refusals: a cycle, one parent, frequencies that do not sum to 1
  const int cfa[3] = {-1, 2, 1}, cmo[3] = {-1, 0, 0};
  const int cg[6] = {0, 0, -1, -1, -1, -1};
  tk_ieq(t, "ped18", clab_pedlik(3, cfa, cmo, cg, 3, p, nullptr, nullptr, &l1), CLAB_EARG);
  const int ofa[3] = {-1, -1, 0}, omo[3] = {-1, -1, -1};
  tk_ieq(t, "ped19", clab_pedlik(3, ofa, omo, cg, 3, p, nullptr, nullptr, &l1), CLAB_EARG);
  const double pb[3] = {0.1, 0.3, 0.5};
  tk_ieq(t, "ped20", clab_pedlik(3, fa, mo, g, 3, pb, nullptr, nullptr, &l1), CLAB_EDOM);
}

}  // namespace

void t_ped(tk &t, int) { run(t); }
