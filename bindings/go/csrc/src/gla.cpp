// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, the dense linear algebra of gla.h: the cyclic Jacobi
// eigendecomposition, LU with partial pivoting and Cholesky.
#include "gla.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace gla {

bool jacobi(double *a, int n, double *val, double *vec) {
  for (std::size_t q = 0; q < static_cast<std::size_t>(n) * n; ++q) {
    if (!std::isfinite(a[q])) return false;
  }
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) vec[at(i, j, n)] = i == j ? 1.0 : 0.0;
  bool done = false;
  for (int sw = 0; sw < 100 && !done; ++sw) {
    if (kutil::stop(static_cast<double>(n) * n * n)) return false;
    double off = 0, tot = 0;
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < n; ++i) {
        const double s = a[at(i, j, n)] * a[at(i, j, n)];
        tot += s;
        if (i != j) off += s;
      }
    if (off <= 1e-30 * tot) {
      done = true;
      break;
    }
    for (int p = 0; p < n - 1; ++p)
      for (int q = p + 1; q < n; ++q) {
        const double apq = a[at(p, q, n)];
        if (apq == 0) continue;
        const double app = a[at(p, p, n)], aqq = a[at(q, q, n)];
        // the rotation that zeroes a[p, q] (the smaller angle)
        const double th = (aqq - app) / (2 * apq);
        const double t = (th >= 0 ? 1.0 : -1.0) / (std::fabs(th) + std::sqrt(th * th + 1));
        const double c = 1 / std::sqrt(t * t + 1), s = t * c;
        for (int k = 0; k < n; ++k) {
          const double akp = a[at(k, p, n)], akq = a[at(k, q, n)];
          a[at(k, p, n)] = c * akp - s * akq;
          a[at(k, q, n)] = s * akp + c * akq;
        }
        for (int k = 0; k < n; ++k) {
          const double apk = a[at(p, k, n)], aqk = a[at(q, k, n)];
          a[at(p, k, n)] = c * apk - s * aqk;
          a[at(q, k, n)] = s * apk + c * aqk;
        }
        a[at(p, q, n)] = a[at(q, p, n)] = 0;
        for (int k = 0; k < n; ++k) {
          const double vkp = vec[at(k, p, n)], vkq = vec[at(k, q, n)];
          vec[at(k, p, n)] = c * vkp - s * vkq;
          vec[at(k, q, n)] = s * vkp + c * vkq;
        }
      }
  }
  if (!done) return false;
  // decreasing order, by selection (n is small: these are dense methods)
  std::vector<double> col(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) val[i] = a[at(i, i, n)];
  for (int i = 0; i < n; ++i) {
    int m = i;
    for (int j = i + 1; j < n; ++j) {
      if (val[j] > val[m]) m = j;
    }
    if (m == i) continue;
    std::swap(val[i], val[m]);
    for (int k = 0; k < n; ++k) std::swap(vec[at(k, i, n)], vec[at(k, m, n)]);
  }
  return true;
}

bool lusolve(double *a, int n, double *b, int m, int *piv) {
  for (int k = 0; k < n; ++k) {
    int r = k;
    for (int i = k + 1; i < n; ++i) {
      if (std::fabs(a[at(i, k, n)]) > std::fabs(a[at(r, k, n)])) r = i;
    }
    const double pv = a[at(r, k, n)];
    if (pv == 0 || !std::isfinite(pv)) return false;
    piv[k] = r;
    if (r != k) {
      for (int j = 0; j < n; ++j) std::swap(a[at(k, j, n)], a[at(r, j, n)]);
      for (int j = 0; j < m; ++j) std::swap(b[at(k, j, n)], b[at(r, j, n)]);
    }
    for (int i = k + 1; i < n; ++i) {
      const double l = a[at(i, k, n)] / pv;
      a[at(i, k, n)] = l;
      if (l == 0) continue;
      for (int j = k + 1; j < n; ++j) a[at(i, j, n)] -= l * a[at(k, j, n)];
      for (int j = 0; j < m; ++j) b[at(i, j, n)] -= l * b[at(k, j, n)];
    }
  }
  for (int j = 0; j < m; ++j)
    for (int i = n - 1; i >= 0; --i) {
      double s = b[at(i, j, n)];
      for (int k = i + 1; k < n; ++k) s -= a[at(i, k, n)] * b[at(k, j, n)];
      b[at(i, j, n)] = s / a[at(i, i, n)];
    }
  for (std::size_t q = 0; q < static_cast<std::size_t>(n) * m; ++q) {
    if (!std::isfinite(b[q])) return false;
  }
  return true;
}

bool chol(double *a, int n) {
  for (int j = 0; j < n; ++j) {
    double d = a[at(j, j, n)];
    for (int k = 0; k < j; ++k) d -= a[at(j, k, n)] * a[at(j, k, n)];
    if (!(d > 0)) return false;
    const double l = std::sqrt(d);
    a[at(j, j, n)] = l;
    for (int i = j + 1; i < n; ++i) {
      double s = a[at(i, j, n)];
      for (int k = 0; k < j; ++k) s -= a[at(i, k, n)] * a[at(j, k, n)];
      a[at(i, j, n)] = s / l;
    }
    for (int i = 0; i < j; ++i) a[at(i, j, n)] = 0;
  }
  return true;
}

void fsolve(const double *l, int n, double *b, int m) {
  for (int j = 0; j < m; ++j)
    for (int i = 0; i < n; ++i) {
      double s = b[at(i, j, n)];
      for (int k = 0; k < i; ++k) s -= l[at(i, k, n)] * b[at(k, j, n)];
      b[at(i, j, n)] = s / l[at(i, i, n)];
    }
}

void bsolve(const double *l, int n, double *b, int m) {
  for (int j = 0; j < m; ++j)
    for (int i = n - 1; i >= 0; --i) {
      double s = b[at(i, j, n)];
      for (int k = i + 1; k < n; ++k) s -= l[at(k, i, n)] * b[at(k, j, n)];
      b[at(i, j, n)] = s / l[at(i, i, n)];
    }
}

}  // namespace gla
