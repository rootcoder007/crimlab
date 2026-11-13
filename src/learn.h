// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The learners of libcrimlab as C++ for the other kernel families (the
// cross-fitting of dml.cpp): elastic net paths and generalised linear
// models by Newton steps. Each fits on a subset of the rows of a
// column-major matrix, given by an index list, so that a fold is never
// copied. No allocation, no exceptions: scratch comes from the caller.
#ifndef CLAB_LEARN_H
#define CLAB_LEARN_H

namespace lrn {

// A design: x column-major, n rows (its leading dimension), p columns;
// y the n outcomes; the rows idx[0..m) are the ones used.
struct Dat {
  const double *x;
  int n, p;
  const double *y;
  const int *idx;
  int m;
};

// Scratch for path(): doubles; and p ints.
long pwork(int n, int p);

// The largest lambda of the path (all coefficients 0 there) for family fam
// (0 gaussian, 1 binomial) and mixing alpha; work: pwork() doubles.
// Returns CLAB_OK or CLAB_EDOM (a constant outcome, a binomial outcome
// that is not 0 or 1).
int lmax(const Dat &d, int fam, double alpha, double *lam, double *work);

// The elastic net path at lam[0..nl) (decreasing) by cyclical coordinate
// descent, columns standardised over the rows used. Out: a0[nl] and
// beta[p nl] (column k for lam[k]) on the original scale. thresh: the
// convergence threshold on the largest weighted squared change of a
// standardised coefficient in a pass; maxit: the most passes in all.
int path(const Dat &d, int fam, double alpha, const double *lam, int nl,
         double thresh, int maxit, double *a0, double *beta, double *work,
         int *iwork);

// The path at lam (or, when 0 < lmr < 1, at nl lambdas from lmax() down
// to lmr times it, written to lam) and, when fold is given (fold[i] in [0,
// nf) for the row d.idx[i]), K-fold cross-validation: cv (2 nl), the mean
// held-out loss and its standard error; best (2), the indices of the
// minimum and of the largest lambda within one standard error. work:
// cvwork() doubles; iwork: d.m + p ints.
long cvwork(int n, int p, int nl, int nf);
int cvnet(const Dat &d, int fam, double alpha, double *lam, int nl,
          double lmr, double thresh, int maxit, const int *fold, int nf,
          double *a0, double *beta, double *cv, int *best, double *work,
          int *iwork);

// A generalised linear model with intercept by Newton steps on the full
// system (gaussian: one least-squares solve; binomial: logistic
// regression). Out: coef[p + 1] (intercept first), se[p + 1] from the
// inverse information, dev (the deviance). work: gwork() doubles. Returns
// CLAB_EDOM for a singular system or no convergence within maxit steps.
long gwork(int p);
int glm(const Dat &d, int fam, double tol, int maxit, double *coef,
        double *se, double *dev, double *work);

// The linear predictor a0 + x[i, ] b for the rows idx[0..m).
void lpred(const double *x, int n, int p, const int *idx, int m, double a0,
           const double *b, double *eta);

}  // namespace lrn

#endif
