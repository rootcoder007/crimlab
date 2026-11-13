// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The linear algebra of the policy-evaluation families (did.cpp, stg.cpp,
// syn.cpp, wcb.cpp): the within transformation of a panel by alternating
// projections, weighted least squares by Cholesky, and cluster-robust
// variances, one-way and two-way. No allocation, no exceptions: scratch
// comes from the caller. Matrices are column-major.
#ifndef CLAB_PEV_H
#define CLAB_PEV_H

#include <vector>

namespace pev {

// Subtracts from each of the k columns of v (n rows) its weighted means
// by the n unit codes u (0..nu-1) and the period codes t (0..nt-1),
// alternately, until the largest unit or period mean left is below tol
// times the column's largest absolute value (two = 1), or its period
// means only, once (two = 0). w: NULL or n positive weights. work: nu +
// nt doubles. Returns the sweeps of the slowest column, or -1 when maxit
// sweeps did not converge.
int within(double *v, int n, int k, const int *u, int nu, const int *t, int nt, int two,
           const double *w, double tol, int maxit, double *work);

// In place: the lower Cholesky factor of the p x p matrix a. false when a
// pivot is not above 1e-12 times the largest diagonal element (a is then
// singular to working precision).
bool chol(double *a, int p);

// From the factor l (chol()): the inverse of l l' into ai (p x p).
void cinv(const double *l, int p, double *ai);

// Weighted least squares of y on the p columns of x (n rows; w NULL or
// n weights): beta (p), xi = (X'WX)^-1 (p x p), r = y - x beta (n).
// work: p p + p doubles. CLAB_EDOM when X'WX is singular.
int ols(const double *x, const double *y, int n, int p, const double *w, double *beta,
        double *xi, double *r, double *work);

// The meat of a cluster sandwich: m (p x p) = sum over the clusters g of
// s_g s_g', s_g = sum over the rows i of g of w_i r_i x_i. cl: n codes
// 0..ncl-1. work: ncl p doubles.
void meat(const double *x, const double *r, const double *w, int n, int p, const int *cl,
          int ncl, double *m, double *work);

// The sandwich xi m xi into v (p x p), times c. work: p p doubles.
void sand(const double *xi, const double *m, int p, double c, double *v, double *work);

// The cluster-robust (CR1) variance of an ols() fit into v (p x p):
// one-way by cl (cl2 NULL), or two-way (Cameron, Gelbach and Miller 2011)
// by cl and cl2 with their intersection c12 (n12 codes): V1 + V2 - V12,
// every term with the factor G / (G - 1) (n - 1) / (n - K) of G, the
// smaller number of clusters. work: (max cluster count) p + 2 p p doubles.
// Returns G.
int crv(const double *x, const double *r, const double *w, int n, int p, const double *xi,
        const int *cl, int ncl, const int *cl2, int ncl2, const int *c12, int n12, double kk,
        double *v, double *work);

// Codes from 0 for pairs (a[i], b[i]), a < na: into c (n), in order of a
// na + b; returns the number of distinct pairs. iwork: na nb ints, or
// fewer when nb is 0 (then c is a).
int pairs(const int *a, const int *b, int n, int na, int nb, int *c, int *iwork);

// The Goodman-Bacon (2021) decomposition of the two-way fixed-effects
// estimate on a balanced panel without weights or covariates: yw (nu x
// nt), the outcome of unit i in period s at yw[s nu + i]; g, the period
// (0..nt-1) each unit is first treated, nt or more for never (treatment
// stays on). Out, one row per two-group comparison: ty (0 treated against
// never, 1 earlier against later, 2 later against earlier), tr and co (the
// adoption periods of the treated and the comparison group, nt for
// never), est (the 2 x 2 estimate) and wt (its weight; they sum to 1).
// Rows: at most k + k (k - 1), k the adoption periods. work: k + 1 + (k +
// 1) nt doubles; iwork: 2 nu ints. Returns the number of rows.
int bacon(const double *yw, int nu, int nt, const int *g, int *ty, int *tr, int *co,
          double *est, double *wt, double *work, int *iwork);

// Nearest-neighbour matching with replacement of the treated units (trt 1)
// to comparison units (trt 0) with the same exact code ex, ratio per
// treated unit, by the Mahalanobis distance of their features f (nu x k)
// with the pooled within-group covariance; with calip > 0, only matches
// within calip standard deviations of all the distances computed. Out:
// the matches (mt, mc, md: treated, comparison, distance; at most n1
// ratio), and the unit weights wu (nu): 1 for a matched treated unit, the
// sum of 1 / (its treated unit's matches) for a comparison unit, 0
// otherwise. work: k k + 3 k doubles. Returns the number of matches, or
// -1 when the covariance is singular.
int mtch(const double *f, int nu, int k, const int *trt, const int *ex, int ratio,
         double calip, int *mt, int *mc, double *md, double *wu, double *work);

// The Gram matrix g (k x k) of the weight problem min (1/m) ||A w + w0 -
// b||^2 + z2 ||w||^2: a (m x k), b (m); icpt 1 centres the columns of A -
// b over the m rows (the free intercept w0), 0 leaves them.
void gram(const double *a, const double *b, int m, int k, int icpt, double z2, double *g);

// Wolfe's minimum-norm point: w (k) on the simplex minimising w'G w for
// the Gram matrix g (k x k). Returns the major cycles taken, or -1 when
// the corral's system is singular or maxit cycles did not end.
int mnp(const double *g, int k, double *w, int maxit);

// The synthetic difference-in-differences estimate (Arkhangelsky et al.
// 2021) of the treated units trd against the comparison units ctl, from
// the wide panel yw (nu x nt, [s nu + i]) adopted at period t0. With fit,
// the unit weights om (one per ctl) and the time weights la (t0) are
// solved for (ridge (n1 t1)^(1/4) sigma and 1e-6 sigma, sigma the standard
// deviation of the comparison units' first differences before t0);
// otherwise the given ones are used. NaN when a weight problem fails.
double sdtau(const double *yw, int nu, int nt, const std::vector<int> &ctl,
             const std::vector<int> &trd, int t0, bool fit, double *om, double *la);

}  // namespace pev

#endif
