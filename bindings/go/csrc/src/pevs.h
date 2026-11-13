// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The staggered-adoption estimators of stg.cpp (stgdid()). Units are
// rows of a balanced panel held wide: the value of unit i in period s is
// at [s nu + i]. A sample of units is a list of indices, with repeats in
// a bootstrap sample.
#ifndef CLAB_PEVS_H
#define CLAB_PEVS_H

#include <cstdint>
#include <vector>

namespace pev {

struct stg {
  const double *yw;   // outcomes, nu x nt
  const int *dw;      // treatments 0 or 1, nu x nt (meth 2)
  const double *xw;   // covariates, (k nt) x nu: column j of period s at [(j nt + s) nu + i]
  const int *g;       // first treated period of each unit (code; nt for never)
  int nu, nt, k;
  int meth;           // 0 Callaway-Sant'Anna, 1 Wooldridge, 2 de Chaisemartin-D'Haultfoeuille
  int ctrl;           // 0 not yet treated, 1 never treated
  int est;            // with covariates: 0 doubly robust, 1 IPW, 2 regression
  std::vector<int> cg, ct;  // the group-time cells: cohort and period codes
};

// The cells (cohort g from 1 to nt - 1 that occurs, period s): every
// period from 1 (post = 0) or the periods from g on (post = 1).
void cells(const int *g, int nu, int nt, std::vector<int> &cg, std::vector<int> &ct, int post);

// The group-time effects of the cells over the units idx: the cohort's
// change from the base period (g - 1 after adoption, s - 1 before) less
// the comparison units' (never treated, or not yet treated by max(s,
// base)), or with covariates the estimator est. NaN for a cell without
// cohort or comparison units.
void csa(const stg &s, const std::vector<int> &idx, double *att);

// The post-adoption cells as coefficients of the two-way regression on
// one indicator per cell; without never-treated units the periods from
// the last adoption on are left out. NaN for a cell not estimated.
void wool(const stg &s, const std::vector<int> &idx, double *att);

// DID_M and its placebo (out[0], out[1]) and their switchers (out[2],
// out[3]); NaN where no unit switches against a stable one.
void dcdh(const stg &s, const std::vector<int> &idx, double *out);

// The overall effect (agg[0]: the post cells weighted by cohort size)
// and the effect at each event time ev (sorted; agg[1..]).
void aggr(const stg &s, const std::vector<int> &idx, const double *att, std::vector<int> &ev,
          double *agg);

// Every statistic of the method over the units idx into v: the cells,
// the overall effect, the event times (methods 0 and 1); or the four
// numbers of dcdh().
void stat(const stg &s, const std::vector<int> &idx, std::vector<int> &ev, double *v);

// The cluster bootstrap: reps draws of ncl clusters with replacement (cl:
// the cluster of each unit), each statistic's standard deviation over the
// draws where it is defined into se (ns).
void boot(const stg &s, const int *cl, int ncl, int reps, std::uint64_t seed,
          std::vector<int> &ev, int ns, double *se);

}  // namespace pev

#endif
