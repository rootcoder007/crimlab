// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Thin-plate spline surface in the plane (crimlab.h, tpsurf()): Duchon
// (1977); Wahba (1990, ch. 2). f(s) = sum_i a_i phi(|s - s_i|) + b_0 +
// b_1 x + b_2 y, phi(r) = r^2 log r, the coefficients solving [K + lam I,
// T; T', 0] [a; b] = [z; 0], T the rows (1, x_i, y_i); lam 0
// interpolates. Solved by LU (gla::lusolve). Work: O(n^3 + n nx ny).
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "crimlab.h"
#include "gla.h"
#include "kutil.h"

namespace {

// phi of the distance between (dx, dy) and 0: r^2 log r = r^2 log(r^2) / 2
inline double phi(double dx, double dy) {
  const double r2 = dx * dx + dy * dy;
  return r2 > 0 ? 0.5 * r2 * std::log(r2) : 0.0;
}

}  // namespace

extern "C" int clab_tpsurf(const double *x, const double *y, const double *z, int n, double lam, const double *gx,
                           int nx, const double *gy, int ny, double *a, double *b, double *fit, double *gz) CLAB_TRY {
  using std::size_t;
  if (n < 3 || nx < 1 || ny < 1 || !(lam >= 0) || !std::isfinite(lam)) return CLAB_EARG;
  const int q = n + 3;
  if (static_cast<double>(q) * q > kutil::maxcell || static_cast<double>(nx) * ny > kutil::maxcell) return CLAB_EARG;
  // points on one line leave the plane undetermined: the largest cross
  // product from the first point, against the farthest point, is at most
  // 1e-12 of its squared distance
  double d2 = 0.0, fx = 0.0, fy = 0.0;
  for (int i = 1; i < n; ++i) {
    const double dx = x[i] - x[0], dy = y[i] - y[0], d = dx * dx + dy * dy;
    if (d > d2) {
      d2 = d;
      fx = dx;
      fy = dy;
    }
  }
  double cr = 0.0;
  for (int i = 1; i < n; ++i) cr = std::max(cr, std::fabs(fx * (y[i] - y[0]) - fy * (x[i] - x[0])));
  if (!(d2 > 0) || cr <= 1e-12 * d2) return CLAB_EDOM;
  const size_t Q = static_cast<size_t>(q), N = static_cast<size_t>(n);
  std::vector<double> m(Q * Q, 0.0), r(Q, 0.0);
  std::vector<int> piv(Q);
  for (size_t j = 0; j < N; ++j) {
    for (size_t i = 0; i < N; ++i) m[j * Q + i] = phi(x[i] - x[j], y[i] - y[j]) + (i == j ? lam : 0.0);
    m[N * Q + j] = 1.0;
    m[(N + 1) * Q + j] = x[j];
    m[(N + 2) * Q + j] = y[j];
    m[j * Q + N] = 1.0;
    m[j * Q + N + 1] = x[j];
    m[j * Q + N + 2] = y[j];
    r[j] = z[j];
  }
  if (kutil::stop(static_cast<double>(q) * q * q)) return CLAB_EINTR;
  if (!gla::lusolve(m.data(), q, r.data(), 1, piv.data())) return CLAB_EDOM;
  for (size_t i = 0; i < Q; ++i) {
    if (!std::isfinite(r[i])) return CLAB_EDOM;
  }
  std::copy(r.begin(), r.begin() + n, a);
  std::copy(r.begin() + n, r.end(), b);
  const auto f = [&](double u, double v) {
    double s = b[0] + b[1] * u + b[2] * v;
    for (int j = 0; j < n; ++j) s += a[j] * phi(u - x[j], v - y[j]);
    return s;
  };
  for (int i = 0; i < n; ++i) fit[i] = f(x[i], y[i]);
  if (kutil::stop(static_cast<double>(nx) * ny * n)) return CLAB_EINTR;
  const size_t NX = static_cast<size_t>(nx);
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) gz[static_cast<size_t>(j) * NX + static_cast<size_t>(i)] = f(gx[i], gy[j]);
  }
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)
