// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The design-based causal estimators of cba.cpp (cbatt(), smdbal(),
// evlx()). A sample of units is a list of indices into the data, with
// repeats in a bootstrap sample.
#ifndef CLAB_PEVC_H
#define CLAB_PEVC_H

#include <vector>

namespace pev {

// Outcome y, treatment d (0 or 1, as doubles) and covariates x (n x p,
// column-major) of n units.
struct cdat {
  const double *y, *d, *x;
  int n, p;
};

// The logistic propensity of the units idx (intercept and x), into e (n;
// NaN outside idx); false when the fit fails (separation).
bool pscore(const cdat &c, const std::vector<int> &idx, std::vector<double> &e);

// Hajek-normalised inverse probability weighting: tgt 0 ATE, 1 ATT, 2 ATC.
double ipw(const cdat &c, const std::vector<int> &idx, const std::vector<double> &e, int tgt);

// g-computation: least squares of y on x in each arm, the mean difference
// of the two predictions over the target units; NaN when a fit fails.
double gcomp(const cdat &c, const std::vector<int> &idx, int tgt);

// Subclassification at the ns quantiles of the propensity among the
// target units: the subclasses' differences of means weighted by their
// target units; NaN when a subclass with target units lacks an arm.
double subcl(const cdat &c, const std::vector<int> &idx, const std::vector<double> &e, int tgt,
             int ns);

// Greedy k:1 nearest-neighbour matching without replacement on the
// propensity, the focal arm (treated for tgt 1, controls for tgt 2) in
// decreasing propensity, within cal (0: none). Out: the pairs (mt focal,
// mc match) and, per matched focal unit, the treated less the control
// (mean of its matches). Returns the matched focal units.
int psm(const cdat &c, const std::vector<double> &e, int tgt, int k, double cal,
        std::vector<int> &mt, std::vector<int> &mc, std::vector<double> &dif);

// Entropy balancing for the ATT (Hainmueller 2012): w (n), 1 for the
// treated and, for the controls of idx, weights summing to the number of
// treated whose weighted means equal the treated means of x. false when
// the dual's Hessian is singular or Newton does not converge.
bool ebal(const cdat &c, const std::vector<int> &idx, std::vector<double> &w);

// The ATT with the weights of ebal().
double ebatt(const cdat &c, const std::vector<int> &idx, const std::vector<double> &w);

// One covariate's balance: out (4) the treated and control means
// (weighted by w, or not when NULL), the standardised mean difference
// over the unweighted pooled standard deviation, and the variance ratio.
void smd(const double *x, int n, const double *d, const double *w, double *out);

// The E-value (VanderWeele and Ding 2017) of an estimate of type 0 (risk
// ratio), 1 (odds ratio), 2 (hazard ratio) or 3 (standardised difference,
// with sd its standard error) and its limits lo, hi (NaN for none); rare:
// the outcome is rare (an OR or HR is then a RR). out (3): the risk ratio,
// its E-value, the E-value of the limit nearer the null (1 when the
// interval holds the null; NaN without limits).
void evalue(double est, double lo, double hi, int type, int rare, double sd, double *out);

}  // namespace pev

#endif
