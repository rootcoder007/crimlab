// SPDX-License-Identifier: AGPL-3.0-or-later
//
// libcrimlab, space-time: Mantel's correlation between the spatial and
// the temporal distances of all pairs of events with a permutation test,
// and the spatial-temporal moving average of a series of events.
#include <cmath>
#include <limits>

#include "cgeo.h"
#include "crimlab.h"
#include "kutil.h"

namespace {

// the pair dissimilarity: the distance, or 1 / (distance + c) when c > 0
inline double tr(double d, double c) { return c > 0.0 ? 1.0 / (d + c) : d; }

// sum over pairs i < j of X_ij Y_{p(i) p(j)}
double cross(const double *x, const double *y, const double *t, const int *p,
             int n, double cs, double ct) {
  kutil::ksum s;
  for (int i = 1; i < n; ++i) {
    const double ti = t[p[i]];
    for (int j = 0; j < i; ++j) {
      const double dx = tr(std::hypot(x[i] - x[j], y[i] - y[j]), cs);
      const double dt = tr(std::fabs(ti - t[p[j]]), ct);
      s.add(dx * dt);
    }
  }
  return s.val();
}

}  // namespace

extern "C" {

int clab_mantst(const double *x, const double *y, const double *t, int n,
                double cs, double ct, int nsim, unsigned int seed, int *iwork,
                double *out) {
  if (n < 3 || nsim < 0 || n > 46340) return CLAB_EARG;
  if (!kutil::nonneg(cs) || !kutil::nonneg(ct)) return CLAB_EDOM;
  for (int i = 0; i < n; ++i)
    if (!std::isfinite(x[i]) || !std::isfinite(y[i]) || !std::isfinite(t[i]))
      return CLAB_EDOM;
  // means and standard deviations of the two pair dissimilarities do not
  // change under a permutation of the times: two passes, once
  const double np = 0.5 * n * (n - 1.0);
  kutil::ksum sx, st;
  for (int i = 1; i < n; ++i)
    for (int j = 0; j < i; ++j) {
      sx.add(tr(std::hypot(x[i] - x[j], y[i] - y[j]), cs));
      st.add(tr(std::fabs(t[i] - t[j]), ct));
    }
  const double mx = sx.val() / np, mt = st.val() / np;
  kutil::ksum vx, vt;
  for (int i = 1; i < n; ++i)
    for (int j = 0; j < i; ++j) {
      const double a = tr(std::hypot(x[i] - x[j], y[i] - y[j]), cs) - mx;
      const double b = tr(std::fabs(t[i] - t[j]), ct) - mt;
      vx.add(a * a);
      vt.add(b * b);
    }
  const double sdx = std::sqrt(vx.val() / np), sdt = std::sqrt(vt.val() / np);
  if (!(sdx > 0.0) || !(sdt > 0.0)) return CLAB_EDOM;  // a constant: r undefined
  int *p = iwork;
  for (int i = 0; i < n; ++i) p[i] = i;
  auto rof = [&](double c) { return (c / np - mx * mt) / (sdx * sdt); };
  const double r0 = rof(cross(x, y, t, p, n, cs, ct));
  out[0] = r0;
  out[1] = np;
  out[2] = out[3] = out[4] = kutil::nan();
  if (nsim > 0) {
    cgeo::rng rg(seed);
    int ge = 0, le = 0, ab = 0;
    // a tolerance for ties of r with the observed value (equal sums may
    // round apart by a few ulps)
    const double eps = 1e-12 * (1.0 + std::fabs(r0));
    for (int s = 0; s < nsim; ++s) {
      for (int i = n - 1; i > 0; --i) {  // Fisher-Yates
        const int j = rg.below(i + 1);
        const int tmp = p[i];
        p[i] = p[j];
        p[j] = tmp;
      }
      const double r = rof(cross(x, y, t, p, n, cs, ct));
      if (r >= r0 - eps) ++ge;
      if (r <= r0 + eps) ++le;
      if (std::fabs(r) >= std::fabs(r0) - eps) ++ab;
    }
    out[2] = (1.0 + ge) / (nsim + 1.0);
    out[3] = (1.0 + le) / (nsim + 1.0);
    out[4] = (1.0 + ab) / (nsim + 1.0);
  }
  return CLAB_OK;
}

int clab_stmavg(const double *x, const double *y, const double *t, int n,
                int span, int *iwork, double *out) {
  if (n < 1 || span < 1 || span % 2 == 0) return CLAB_EARG;
  for (int i = 0; i < n; ++i)
    if (!std::isfinite(x[i]) || !std::isfinite(y[i]) || !std::isfinite(t[i]))
      return CLAB_EDOM;
  int *o = iwork;
  for (int i = 0; i < n; ++i) o[i] = i;
  cgeo::hsort(o, n, [t](int a, int b) { return t[a] < t[b] || (t[a] == t[b] && a < b); });
  const int h = span / 2;
  for (int q = 0; q < n; ++q) {
    const int a = q - h < 0 ? 0 : q - h, b = q + h > n - 1 ? n - 1 : q + h;
    kutil::ksum sx, sy;
    for (int k = a; k <= b; ++k) {
      sx.add(x[o[k]]);
      sy.add(y[o[k]]);
    }
    const int m = b - a + 1;
    double *r = out + 4 * q;
    r[0] = o[q];
    r[1] = sx.val() / m;
    r[2] = sy.val() / m;
    r[3] = m;
  }
  return CLAB_OK;
}

}  // extern "C"
