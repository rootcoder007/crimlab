// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Dense linear algebra and the two minimisers of the geostatistics and
// spatial epidemiology families (krg.cpp, sep.cpp): the symmetric
// eigendecomposition by cyclic Jacobi rotations, LU and Cholesky solves,
// and the Nelder-Mead and Brent minimisers in the form R's optim() and
// optimize() use, so that a fit stops where R's would. Matrices are
// column-major. Nothing here allocates except the eigendecomposition's
// and the simplex's own small arrays.
#ifndef CLAB_GLA_H
#define CLAB_GLA_H

#include <cfloat>
#include <cmath>
#include <cstddef>
#include <vector>

#include "kutil.h"

namespace gla {

inline std::size_t at(int i, int j, int n) {
  return static_cast<std::size_t>(j) * static_cast<std::size_t>(n) + static_cast<std::size_t>(i);
}

// The eigenvalues of the symmetric n x n matrix a (destroyed) in
// decreasing order into val, the unit eigenvectors as the columns of vec
// (n x n), by cyclic Jacobi rotations (Golub and Van Loan 2013, sec. 8.5)
// until the off-diagonal mass is below 1e-30 of the total. False when a
// is not finite or the sweeps do not converge.
bool jacobi(double *a, int n, double *val, double *vec);

// LU with partial pivoting of the n x n matrix a (overwritten) and the
// solution of a z = b for the m columns of b (n x m, overwritten), as
// LAPACK's dgesv: false when a pivot is 0 or a value not finite. piv: n.
bool lusolve(double *a, int n, double *b, int m, int *piv);

// In place: the lower Cholesky factor l of the symmetric n x n matrix a
// (its lower triangle read), as R's chol(): false when a pivot is not
// positive. The upper triangle is set to 0.
bool chol(double *a, int n);
// l z = b and l' z = b for the m columns of b (n x m), in place.
void fsolve(const double *l, int n, double *b, int m);
void bsolve(const double *l, int n, double *b, int m);

// The Nelder-Mead simplex of R's optim() (nmmin, after Nash 1990, with
// reflection 1, contraction 0.5 and expansion 2): the start x (n values)
// is replaced by the minimum found; a non-finite value of f counts as
// 1e35. Stops when the spread of the simplex is at most reltol (|f| +
// reltol), or after maxit evaluations (return 1), or when a shrink does
// not reduce the simplex (return 10); 0 on convergence; 2, x unchanged,
// when the start or its value is not finite. *fx: the value at
// x; *nf: the evaluations.
template <class F>
int nmmin(int n, double *x, double *fx, F f, double reltol, int maxit, int *nf) {
  const double big = 1.0e35;
  const int n1 = n + 1, C = n + 2;
  // P: n + 1 rows (the parameters, then f), n + 2 columns
  std::vector<double> P(static_cast<std::size_t>(n1) * C);
  auto p = [&](int i, int j) -> double & { return P[at(i, j, n1)]; };
  std::vector<double> b(x, x + n);
  double v = f(b.data());
  int cnt = 1, fail = 0;
  bool fin = true;
  for (int i = 0; i < n; ++i) fin = fin && std::isfinite(b[i]);
  // an infinite start has no simplex: b + step is b however large the step
  if (!fin || !std::isfinite(v)) {
    *fx = v;
    *nf = cnt;
    return 2;
  }
  const double conv = reltol * (std::fabs(v) + reltol);
  p(n1 - 1, 0) = v;
  for (int i = 0; i < n; ++i) p(i, 0) = b[i];
  int L = 1, H = 1;
  double size = 0, step = 0;
  for (int i = 0; i < n; ++i) {
    if (0.1 * std::fabs(b[i]) > step) step = 0.1 * std::fabs(b[i]);
  }
  if (step == 0.0) step = 0.1;
  for (int j = 2; j <= n1; ++j) {
    for (int i = 0; i < n; ++i) p(i, j - 1) = b[i];
    double tr = step;
    while (p(j - 2, j - 1) == b[j - 2]) {
      p(j - 2, j - 1) = b[j - 2] + tr;
      tr *= 10;
    }
    size += tr;
  }
  double old = size;
  bool calc = true;
  do {
    if (kutil::stop(1e4)) break;
    if (calc) {
      for (int j = 0; j < n1; ++j) {
        if (j + 1 == L) continue;
        for (int i = 0; i < n; ++i) b[i] = p(i, j);
        v = f(b.data());
        if (!std::isfinite(v)) v = big;
        ++cnt;
        p(n1 - 1, j) = v;
      }
      calc = false;
    }
    double VL = p(n1 - 1, L - 1), VH = VL;
    H = L;
    for (int j = 1; j <= n1; ++j) {
      if (j == L) continue;
      v = p(n1 - 1, j - 1);
      if (v < VL) {
        L = j;
        VL = v;
      }
      if (v > VH) {
        H = j;
        VH = v;
      }
    }
    if (VH <= VL + conv) break;
    for (int i = 0; i < n; ++i) {
      double t = -p(i, H - 1);
      for (int j = 0; j < n1; ++j) t += p(i, j);
      p(i, C - 1) = t / n;
    }
    for (int i = 0; i < n; ++i) b[i] = 2.0 * p(i, C - 1) - p(i, H - 1);
    v = f(b.data());
    if (!std::isfinite(v)) v = big;
    ++cnt;
    const double VR = v;
    if (VR < VL) {
      p(n1 - 1, C - 1) = v;
      for (int i = 0; i < n; ++i) {
        const double e = 2.0 * b[i] - p(i, C - 1);
        p(i, C - 1) = b[i];
        b[i] = e;
      }
      v = f(b.data());
      if (!std::isfinite(v)) v = big;
      ++cnt;
      if (v < VR) {
        for (int i = 0; i < n; ++i) p(i, H - 1) = b[i];
        p(n1 - 1, H - 1) = v;
      } else {
        for (int i = 0; i < n; ++i) p(i, H - 1) = p(i, C - 1);
        p(n1 - 1, H - 1) = VR;
      }
    } else {
      if (VR < VH) {
        for (int i = 0; i < n; ++i) p(i, H - 1) = b[i];
        p(n1 - 1, H - 1) = VR;
      }
      for (int i = 0; i < n; ++i) b[i] = 0.5 * p(i, H - 1) + 0.5 * p(i, C - 1);
      v = f(b.data());
      if (!std::isfinite(v)) v = big;
      ++cnt;
      if (v < p(n1 - 1, H - 1)) {
        for (int i = 0; i < n; ++i) p(i, H - 1) = b[i];
        p(n1 - 1, H - 1) = v;
      } else if (VR >= VH) {
        calc = true;
        size = 0.0;
        for (int j = 0; j < n1; ++j) {
          if (j + 1 == L) continue;
          for (int i = 0; i < n; ++i) {
            p(i, j) = 0.5 * (p(i, j) - p(i, L - 1)) + p(i, L - 1);
            size += std::fabs(p(i, j) - p(i, L - 1));
          }
        }
        if (size < old) {
          old = size;
        } else {
          fail = 10;
          break;
        }
      }
    }
  } while (cnt <= maxit);
  *fx = p(n1 - 1, L - 1);
  for (int i = 0; i < n; ++i) x[i] = p(i, L - 1);
  if (cnt > maxit) fail = 1;
  *nf = cnt;
  return fail;
}

// Brent's minimiser of f on [lo, hi] as R's optimize() runs it (Brent
// 1973, ch. 5; R's Brent_fmin): golden-section and parabolic steps until
// the bracket is within 2 (sqrt(eps) |x| + tol / 3) of the best point.
template <class F>
double fmin(F f, double lo, double hi, double tol) {
  const double c = (3. - std::sqrt(5.)) * .5, eps = std::sqrt(DBL_EPSILON), tol3 = tol / 3.;
  double a = lo, b = hi, v = a + c * (b - a), w = v, x = v, d = 0., e = 0.;
  double fx = f(x), fv = fx, fw = fx;
  for (int it = 0; it < 10000; ++it) {  // Brent's steps: a cap, never reached
    if (kutil::stop(1e4)) break;
    const double xm = (a + b) * .5, tol1 = eps * std::fabs(x) + tol3, t2 = tol1 * 2.;
    if (std::fabs(x - xm) <= t2 - (b - a) * .5) break;
    double p = 0., q = 0., r = 0.;
    if (std::fabs(e) > tol1) {
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = (q - r) * 2.;
      if (q > 0.)
        p = -p;
      else
        q = -q;
      r = e;
      e = d;
    }
    double u;
    if (std::fabs(p) >= std::fabs(q * .5 * r) || p <= q * (a - x) || p >= q * (b - x)) {
      e = x < xm ? b - x : a - x;
      d = c * e;
    } else {
      d = p / q;
      u = x + d;
      if (u - a < t2 || b - u < t2) {
        d = tol1;
        if (x >= xm) d = -d;
      }
    }
    if (std::fabs(d) >= tol1)
      u = x + d;
    else if (d > 0.)
      u = x + tol1;
    else
      u = x - tol1;
    const double fu = f(u);
    if (fu <= fx) {
      if (u < x)
        b = x;
      else
        a = x;
      v = w;
      w = x;
      x = u;
      fv = fw;
      fw = fx;
      fx = fu;
    } else {
      if (u < x)
        a = u;
      else
        b = u;
      if (fu <= fw || w == x) {
        v = w;
        fv = fw;
        w = u;
        fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v = u;
        fv = fu;
      }
    }
  }
  return x;
}

}  // namespace gla

#endif
