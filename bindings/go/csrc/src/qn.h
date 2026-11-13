// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Numerical helpers of the point-process and series kernels (hwk.cpp,
// trn.cpp, kmr.cpp): a quasi-Newton minimiser, the inverse of a symmetric
// positive definite matrix, the covariance from a numerical Hessian of an
// analytic gradient, and the Kolmogorov-Smirnov p-value. No allocation:
// the caller passes the work arrays.
#ifndef CLAB_QN_H
#define CLAB_QN_H

namespace qn {

// An objective: f(u), and its gradient into g, for p parameters; a value
// that is not finite marks a point outside the domain.
typedef double (*fn)(const double *u, double *g, void *ctx);

// Minimises f from u (p values, overwritten with the minimum) by BFGS on
// the inverse Hessian with a backtracking (Armijo) line search, at most
// maxit steps, until the largest gradient element is below gtol (1 +
// |f|). Out: fv, the value at u; iter, the steps taken. Returns 1 when
// converged, 0 when the steps ran out or the line search found no lower
// point, -1 when f is not finite at the start. work: p^2 + 5 p doubles.
int min(int p, double *u, fn f, void *ctx, int maxit, double gtol, double *fv, int *iter,
        double *work);

// The inverse of the symmetric positive definite p x p matrix a, in place
// (column-major, both triangles). Returns 0, or 1 when a is not positive
// definite (a is then left in a partial state).
int spdinv(double *a, int p);

// The covariance at th (p values) as the inverse of the Hessian of a
// negative log-likelihood, the Hessian by central differences of the
// analytic gradient gr (steps 1e-5 max(1, |th_k|)), symmetrised. th is
// restored on return. Returns 0, or 1 when the Hessian is not positive
// definite or a gradient is not finite (cov is then all NaN). work: 3 p
// doubles.
int hcov(int p, double *th, fn gr, void *ctx, double *cov, double *work);

// The p-value of the two-sided one-sample Kolmogorov-Smirnov statistic d
// on n points, as R's ks.test() gives it: exact (Marsaglia, Tsang and
// Wang 2003) when exact and n < 100, otherwise from the limiting
// distribution of sqrt(n) d (tolerance 1e-6). work: 3 m^2 doubles with
// m = 2 floor(n d) + 1 when exact (at most 3 (2 n + 1)^2).
double kspv(double d, int n, bool exact, double *work);

// The statistic of that test for the n values u, sorted in place, against
// the uniform on [0, 1].
double ksd(double *u, int n);

}  // namespace qn

#endif
