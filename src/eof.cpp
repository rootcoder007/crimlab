// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Empirical orthogonal functions of a space-time matrix (crimlab.h,
// eofst()): Lorenz (1956); Hannachi, Jolliffe and Stephenson (2007). The
// columns of z (m periods by p places) are centred and the centred
// matrix A = U D V' decomposed through the symmetric eigendecomposition
// (gla::jacobi) of the smaller of A'A and A A'. The EOFs are the columns
// of V and the principal component series A V, each signed so that its
// largest loading in absolute value is positive. Work: O(m p s + s^3), s
// = min(m, p).
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "crimlab.h"
#include "gla.h"
#include "kutil.h"

extern "C" int clab_eofst(const double *zm, int m, int p, int k, double *ctr, double *eof, double *pc, double *d2,
                          double *tot, int *nk) CLAB_TRY {
  using std::size_t;
  if (m < 2 || p < 1 || k < 0) return CLAB_EARG;
  const int s = std::min(m, p);
  if (static_cast<double>(m) * p > kutil::maxcell || static_cast<double>(s) * s > kutil::maxcell) return CLAB_EARG;
  const size_t M = static_cast<size_t>(m), P = static_cast<size_t>(p), S = static_cast<size_t>(s);
  std::vector<double> a(M * P);
  double t = 0.0;
  for (size_t j = 0; j < P; ++j) {
    double c = 0.0;
    for (size_t i = 0; i < M; ++i) c += zm[j * M + i];
    c /= static_cast<double>(m);
    ctr[j] = c;
    for (size_t i = 0; i < M; ++i) {
      const double e = zm[j * M + i] - c;
      a[j * M + i] = e;
      t += e * e;
    }
  }
  *tot = t;
  *nk = 0;
  if (!std::isfinite(t)) return CLAB_EARG;
  if (!(t > 0)) return CLAB_EDOM;
  if (kutil::stop(static_cast<double>(m) * p * s)) return CLAB_EINTR;
  // the cross-product of the shorter side: A'A when p <= m, else A A'
  const bool tall = p <= m;
  std::vector<double> g(S * S), val(S), vec(S * S);
  for (size_t u = 0; u < S; ++u) {
    for (size_t w = 0; w <= u; ++w) {
      double q = 0.0;
      if (tall) {
        for (size_t i = 0; i < M; ++i) q += a[u * M + i] * a[w * M + i];
      } else {
        for (size_t j = 0; j < P; ++j) q += a[j * M + u] * a[j * M + w];
      }
      g[w * S + u] = q;
      g[u * S + w] = q;
    }
  }
  if (!gla::jacobi(g.data(), s, val.data(), vec.data()) || !(val[0] > 0)) return CLAB_EDOM;
  // the components whose squared singular value is above 1e-10 of the
  // largest (the eigenvalues carry an error of a few eps times the
  // largest, so a rank-deficient direction stays below it)
  int r = 0;
  while (r < s && val[static_cast<size_t>(r)] > 1e-10 * val[0]) ++r;
  const int kk = k > 0 ? std::min(k, r) : r;
  for (int ci = 0; ci < kk; ++ci) {
    const size_t c = static_cast<size_t>(ci);
    double *v = eof + c * P;
    if (tall) {
      for (size_t j = 0; j < P; ++j) v[j] = vec[c * S + j];
    } else {
      // v = A' u / d
      const double d = std::sqrt(val[c]);
      for (size_t j = 0; j < P; ++j) {
        double q = 0.0;
        for (size_t i = 0; i < M; ++i) q += a[j * M + i] * vec[c * S + i];
        v[j] = q / d;
      }
    }
    size_t b = 0;
    for (size_t j = 1; j < P; ++j) {
      if (std::fabs(v[j]) > std::fabs(v[b])) b = j;
    }
    if (v[b] < 0) {
      for (size_t j = 0; j < P; ++j) v[j] = -v[j];
    }
    for (size_t i = 0; i < M; ++i) {
      double q = 0.0;
      for (size_t j = 0; j < P; ++j) q += a[j * M + i] * v[j];
      pc[c * M + i] = q;
    }
    d2[c] = val[c];
  }
  *nk = kk;
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)
