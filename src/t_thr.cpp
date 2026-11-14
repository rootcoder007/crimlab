// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/thr.cpp: hit rates, the chi-square and the rank-order z by
// hand, and the threshold model recovering the thresholds of counts made
// exactly from its own formulas.
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "tchk.h"

namespace {

using vd = std::vector<double>;
using vi = std::vector<int>;

double up(double z) { return 0.5 * std::erfc(z / std::sqrt(2.0)); }
double lg(double x) { return 1.0 / (1.0 + std::exp(-x)); }

void hand(tk &T) {
  const vd st = {100, 100}, sr = {20, 40}, ht = {8, 10};
  const vi g = {0, 1}, a = {0, 0};
  vd gs(12), hz(8), chi(3), thr(2), tse(6), cf(10), lrt(6);
  vi cv(4);
  tk_ieq(T, "thr.ok",
         clab_thrtst(st.data(), sr.data(), ht.data(), g.data(), a.data(), 2, 2, 1, 0, 0, gs.data(),
                     hz.data(), chi.data(), thr.data(), tse.data(), cf.data(), lrt.data(), cv.data()),
         CLAB_OK);
  tk_eq(T, "thr.hr", gs[9], 0.25, 1e-15);
  tk_eq(T, "thr.sr", gs[6], 0.2, 1e-15);
  tk_eq(T, "thr.diff", hz[1], -0.15, 1e-14);
  tk_eq(T, "thr.se", hz[3], std::sqrt(0.4 * 0.6 / 20 + 0.25 * 0.75 / 40), 1e-14);
  tk_eq(T, "thr.chi", chi[0], 10.0 / 7.0, 1e-14);
  tk_eq(T, "thr.df", chi[1], 1, 0);
  tk_nan(T, "thr.refz", hz[4]);
  // refusals: a search count above the stops; one group with no searches
  const vd sr2 = {120, 40}, sr3 = {0, 40}, ht3 = {0, 10};
  tk_ieq(T, "thr.dom",
         clab_thrtst(st.data(), sr2.data(), ht.data(), g.data(), a.data(), 2, 2, 1, 0, 0, gs.data(),
                     hz.data(), chi.data(), thr.data(), tse.data(), cf.data(), lrt.data(), cv.data()),
         CLAB_EDOM);
  tk_ieq(T, "thr.none",
         clab_thrtst(st.data(), sr3.data(), ht3.data(), g.data(), a.data(), 2, 2, 1, 0, 0, gs.data(),
                     hz.data(), chi.data(), thr.data(), tse.data(), cf.data(), lrt.data(), cv.data()),
         CLAB_EDOM);
  // the model needs more rates than parameters
  tk_ieq(T, "thr.ident",
         clab_thrtst(st.data(), sr.data(), ht.data(), g.data(), a.data(), 2, 2, 1, 0, 100, gs.data(),
                     hz.data(), chi.data(), thr.data(), tse.data(), cf.data(), lrt.data(), cv.data()),
         CLAB_EDOM);
  // rank order: z of -2 and 3.65 for the search rates of two officer groups
  const vd s4 = {100, 100, 100, 100}, k4 = {10, 20, 30, 10}, h4 = {5, 5, 5, 5};
  const vi o4 = {0, 1, 0, 1}, g4 = {0, 0, 1, 1};
  vd rate(20), res(6);
  tk_ieq(T, "rk.ok",
         clab_thrrk(s4.data(), k4.data(), h4.data(), o4.data(), g4.data(), 4, 2, 2, 1.959963984540054,
                    rate.data(), res.data()),
         CLAB_OK);
  tk_eq(T, "rk.zmin", res[0], -2.0, 1e-14);
  tk_eq(T, "rk.zmax", res[2], 0.2 / std::sqrt(0.003), 1e-14);
  tk_eq(T, "rk.rej", res[4], 1, 0);
  tk_eq(T, "rk.hit", res[5], 0, 0);
}

// counts made from the model's own search and hit rates: the fit gives
// back the thresholds, and with the same thresholds no evidence of a
// difference
void model(tk &T, double dt) {
  const int ng = 2, na = 4;
  const double ag[2] = {-1.2, -0.6}, bg[2] = {0.4, 0.2}, cg[2] = {-1.0, -1.0 + dt};
  const double aa[4] = {0, 0.3, -0.4, 0.2}, ba[4] = {0, 0.15, -0.1, 0.05}, ca[4] = {0, 0.2, -0.3, 0.1};
  vd st, sr, ht;
  vi g, a;
  double want[2] = {0, 0};
  for (int q = 0; q < na; ++q) {
    for (int j = 0; j < ng; ++j) {
      const double ph = lg(ag[j] + aa[q]), dl = std::exp(bg[j] + ba[q]), gg = cg[j] + ca[q];
      const double x = (gg - (ag[j] + aa[q])) / dl + 0.5 * dl;
      const double A = ph * up(x - dl), B = (1 - ph) * up(x);
      st.push_back(20000);
      sr.push_back(20000 * (A + B));
      ht.push_back(20000 * A);
      g.push_back(j);
      a.push_back(q);
      want[j] += lg(gg) / na;
    }
  }
  vd gs(12), hz(8), chi(3), thr(2), tse(6), cf(40), lrt(6);
  vi cv(4);
  tk_ieq(T, "thm.ok",
         clab_thrtst(st.data(), sr.data(), ht.data(), g.data(), a.data(), 8, ng, na, 0, 500, gs.data(),
                     hz.data(), chi.data(), thr.data(), tse.data(), cf.data(), lrt.data(), cv.data()),
         CLAB_OK);
  tk_ieq(T, "thm.conv", cv[0], 1);
  tk_near(T, "thm.t0", thr[0], want[0], 1e-6);
  tk_near(T, "thm.t1", thr[1], want[1], 1e-6);
  tk_near(T, "thm.srate", cf[3 * 8 + 5], sr[5] / st[5], 1e-6);
  if (dt == 0) {
    tk_lt(T, "thm.lrt0", lrt[0], 1e-4);
  } else {
    tk_ge(T, "thm.lrt1", lrt[0], 0);
    tk_ge(T, "thm.ll", lrt[3], lrt[4]);
    tk_gt(T, "thm.se", tse[2 + 1], 0);
  }
  tk_eq(T, "thm.cells", lrt[5], 8, 0);
}

}  // namespace

void t_thr(tk &T, int big) {
  (void)big;
  hand(T);
  model(T, 0.0);
  model(T, 0.5);
}
