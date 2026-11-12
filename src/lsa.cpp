// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, spatial autocorrelation on a weights structure (compressed
// sparse rows, src/swt.cpp): Geary's C and the global G of Getis and Ord
// with their moments, local Moran with conditional permutations, and the
// local G and G* of Getis and Ord. Weights rows must hold distinct
// neighbours in increasing order, none of them the unit itself.
#include <cmath>

#include "cgeo.h"
#include "crimlab.h"
#include "kutil.h"

namespace {

// checks of a weights structure: rows non-decreasing, columns in range,
// increasing within a row, no unit its own neighbour, finite weights >= 0
int wcheck(int n, const int *rp, const int *ci, const double *w) {
  if (n < 1 || rp[0] != 0) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (rp[i + 1] < rp[i]) return CLAB_EARG;
    for (int q = rp[i]; q < rp[i + 1]; ++q) {
      if (ci[q] < 0 || ci[q] >= n) return CLAB_EARG;
      if (ci[q] == i || (q > rp[i] && ci[q] <= ci[q - 1])) return CLAB_EDOM;
      if (!kutil::nonneg(w[q])) return CLAB_EDOM;
    }
  }
  return CLAB_OK;
}

// w_ji, 0 when j has no edge to i (binary search in row j)
double wat(const int *rp, const int *ci, const double *w, int j, int i) {
  int lo = rp[j], hi = rp[j + 1];
  while (lo < hi) {
    const int m = lo + (hi - lo) / 2;
    if (ci[m] < i) lo = m + 1; else hi = m;
  }
  return lo < rp[j + 1] && ci[lo] == i ? w[lo] : 0.0;
}

// S0, S1, S2 of Cliff and Ord (1981); cs: n doubles of scratch
void sums(int n, const int *rp, const int *ci, const double *w, double *cs,
          double *s) {
  kutil::ksum s0, s1, s2;
  for (int i = 0; i < n; ++i) cs[i] = 0.0;
  for (int i = 0; i < n; ++i) {
    for (int q = rp[i]; q < rp[i + 1]; ++q) {
      s0.add(w[q]);
      s1.add(w[q] * w[q] + w[q] * wat(rp, ci, w, ci[q], i));
      cs[ci[q]] += w[q];
    }
  }
  for (int i = 0; i < n; ++i) {
    double r = 0.0;
    for (int q = rp[i]; q < rp[i + 1]; ++q) r += w[q];
    s2.add((r + cs[i]) * (r + cs[i]));
  }
  s[0] = s0.val();
  s[1] = s1.val();
  s[2] = s2.val();
}

double upper(double z) { return 0.5 * std::erfc(z / std::sqrt(2.0)); }

// Geary's numerator sum_ij w_ij (v_i - v_j)^2
double gnum(const double *v, int n, const int *rp, const int *ci, const double *w) {
  kutil::ksum s;
  for (int i = 0; i < n; ++i)
    for (int q = rp[i]; q < rp[i + 1]; ++q) {
      const double d = v[i] - v[ci[q]];
      s.add(w[q] * d * d);
    }
  return s.val();
}

// sum_ij w_ij v_i v_j
double gcross(const double *v, int n, const int *rp, const int *ci, const double *w) {
  kutil::ksum s;
  for (int i = 0; i < n; ++i)
    for (int q = rp[i]; q < rp[i + 1]; ++q) s.add(w[q] * v[i] * v[ci[q]]);
  return s.val();
}

void shuffle(double *v, int n, kutil::rng &r) {
  for (int i = n - 1; i > 0; --i) {
    const int j = r.below(i + 1);
    const double t = v[i];
    v[i] = v[j];
    v[j] = t;
  }
}

}  // namespace

extern "C" {

int clab_gearc(const double *x, int n, const int *rp, const int *ci,
               const double *w, int nsim, unsigned int seed, double *work,
               double *out) {
  const int rc = wcheck(n, rp, ci, w);
  if (rc != CLAB_OK) return rc;
  if (n < 4 || nsim < 0) return CLAB_EARG;
  for (int i = 0; i < n; ++i)
    if (!std::isfinite(x[i])) return CLAB_EDOM;
  double s[3];
  sums(n, rp, ci, w, work, s);
  kutil::ksum mx;
  for (int i = 0; i < n; ++i) mx.add(x[i]);
  const double m = mx.val() / n;
  kutil::ksum q2, q4;
  for (int i = 0; i < n; ++i) {
    const double z = x[i] - m;
    q2.add(z * z);
    q4.add(z * z * z * z);
  }
  const double m2 = q2.val(), dn = n;
  if (!(m2 > 0.0) || !(s[0] > 0.0)) return CLAB_EDOM;
  const double b2 = dn * q4.val() / (m2 * m2), S0 = s[0], S1 = s[1], S2 = s[2];
  const double c = (dn - 1.0) * gnum(x, n, rp, ci, w) / (2.0 * S0 * m2);
  const double vn = ((2.0 * S1 + S2) * (dn - 1.0) - 4.0 * S0 * S0) /
                    (2.0 * (dn + 1.0) * S0 * S0);
  const double vr =
      ((dn - 1.0) * S1 * (dn * dn - 3.0 * dn + 3.0 - (dn - 1.0) * b2) -
       0.25 * (dn - 1.0) * S2 * (dn * dn + 3.0 * dn - 6.0 - (dn * dn - dn + 2.0) * b2) +
       S0 * S0 * (dn * dn - 3.0 - (dn - 1.0) * (dn - 1.0) * b2)) /
      (dn * (dn - 2.0) * (dn - 3.0) * S0 * S0);
  out[0] = c;
  out[1] = vn;
  out[2] = vr;
  out[3] = (1.0 - c) / std::sqrt(vn);
  out[4] = (1.0 - c) / std::sqrt(vr);
  out[5] = upper(out[3]);
  out[6] = upper(out[4]);
  out[7] = kutil::nan();
  if (nsim > 0) {
    double *v = work;
    for (int i = 0; i < n; ++i) v[i] = x[i];
    kutil::rng r(seed);
    const double num = c * 2.0 * S0 * m2 / (dn - 1.0);
    const double eps = 1e-12 * (1.0 + num);
    int le = 0;
    for (int k = 0; k < nsim; ++k) {
      shuffle(v, n, r);
      if (gnum(v, n, rp, ci, w) <= num + eps) ++le;
    }
    out[7] = (1.0 + le) / (nsim + 1.0);
  }
  out[8] = S0;
  return CLAB_OK;
}

int clab_gorg(const double *x, int n, const int *rp, const int *ci,
              const double *w, int nsim, unsigned int seed, double *work,
              double *out) {
  const int rc = wcheck(n, rp, ci, w);
  if (rc != CLAB_OK) return rc;
  if (n < 4 || nsim < 0) return CLAB_EARG;
  kutil::ksum a1, a2, a3, a4;
  for (int i = 0; i < n; ++i) {
    if (!kutil::nonneg(x[i])) return CLAB_EDOM;
    a1.add(x[i]);
    a2.add(x[i] * x[i]);
    a3.add(x[i] * x[i] * x[i]);
    a4.add(x[i] * x[i] * x[i] * x[i]);
  }
  double s[3];
  sums(n, rp, ci, w, work, s);
  const double m1 = a1.val(), m2 = a2.val(), m3 = a3.val(), m4 = a4.val();
  const double den = m1 * m1 - m2, dn = n, W = s[0], S1 = s[1], S2 = s[2];
  if (!(den > 0.0) || !(W > 0.0)) return CLAB_EDOM;
  const double g = gcross(x, n, rp, ci, w) / den;
  const double e = W / (dn * (dn - 1.0));
  const double b0 = (dn * dn - 3.0 * dn + 3.0) * S1 - dn * S2 + 3.0 * W * W;
  const double b1 = -((dn * dn - dn) * S1 - 2.0 * dn * S2 + 6.0 * W * W);
  const double b2 = -(2.0 * dn * S1 - (dn + 3.0) * S2 + 6.0 * W * W);
  const double b3 = 4.0 * (dn - 1.0) * S1 - 2.0 * (dn + 1.0) * S2 + 8.0 * W * W;
  const double b4 = S1 - S2 + W * W;
  const double n4 = dn * (dn - 1.0) * (dn - 2.0) * (dn - 3.0);
  const double eg2 = (b0 * m2 * m2 + b1 * m4 + b2 * m1 * m1 * m2 + b3 * m1 * m3 +
                      b4 * m1 * m1 * m1 * m1) /
                     (den * den * n4);
  const double v = eg2 - e * e;
  out[0] = g;
  out[1] = e;
  out[2] = v;
  out[3] = (g - e) / std::sqrt(v);
  out[4] = upper(out[3]);
  out[5] = kutil::nan();
  if (nsim > 0) {
    double *u = work;
    for (int i = 0; i < n; ++i) u[i] = x[i];
    kutil::rng r(seed);
    const double num = g * den, eps = 1e-12 * (1.0 + num);
    int ge = 0;
    for (int k = 0; k < nsim; ++k) {
      shuffle(u, n, r);
      if (gcross(u, n, rp, ci, w) >= num - eps) ++ge;
    }
    out[5] = (1.0 + ge) / (nsim + 1.0);
  }
  return CLAB_OK;
}

int clab_lmorn(const double *x, int n, const int *rp, const int *ci,
               const double *w, int nsim, unsigned int seed, int *iwork,
               double *work, double *out) {
  const int rc = wcheck(n, rp, ci, w);
  if (rc != CLAB_OK) return rc;
  if (n < 3 || nsim < 0) return CLAB_EARG;
  for (int i = 0; i < n; ++i)
    if (!std::isfinite(x[i])) return CLAB_EDOM;
  kutil::ksum mx;
  for (int i = 0; i < n; ++i) mx.add(x[i]);
  const double mean = mx.val() / n, dn = n;
  double *z = work;
  kutil::ksum q2, q4;
  for (int i = 0; i < n; ++i) {
    z[i] = x[i] - mean;
    q2.add(z[i] * z[i]);
    q4.add(z[i] * z[i] * z[i] * z[i]);
  }
  const double m2 = q2.val() / dn;
  if (!(m2 > 0.0)) return CLAB_EDOM;
  const double b2 = (q4.val() / dn) / (m2 * m2);
  double *I = out, *E = out + n, *V = out + 2 * n, *Z = out + 3 * n,
         *P = out + 4 * n, *PP = out + 5 * n, *Q = out + 6 * n;
  for (int i = 0; i < n; ++i) {
    double lag = 0.0, wi = 0.0, wi2 = 0.0;
    for (int q = rp[i]; q < rp[i + 1]; ++q) {
      lag += w[q] * z[ci[q]];
      wi += w[q];
      wi2 += w[q] * w[q];
    }
    I[i] = z[i] / m2 * lag;
    E[i] = -wi / (dn - 1.0);
    V[i] = wi2 * (dn - b2) / (dn - 1.0) +
           (wi * wi - wi2) * (2.0 * b2 - dn) / ((dn - 1.0) * (dn - 2.0)) -
           E[i] * E[i];
    Z[i] = (I[i] - E[i]) / std::sqrt(V[i]);
    P[i] = 2.0 * upper(std::fabs(Z[i]));
    PP[i] = kutil::nan();
    // quadrant: 1 high-high, 2 low-low, 3 low-high, 4 high-low, 0 neither
    Q[i] = z[i] > 0 && lag > 0 ? 1 : z[i] < 0 && lag < 0 ? 2 : z[i] < 0 && lag > 0 ? 3
           : z[i] > 0 && lag < 0 ? 4 : 0;
  }
  if (nsim > 0) {
    // conditional permutation (Anselin 1995): unit i keeps its value; its
    // neighbours' values are drawn without replacement from the others
    int *pool = iwork;
    kutil::rng r(seed);
    for (int i = 0; i < n; ++i) {
      const int k = rp[i + 1] - rp[i];
      int m = 0;
      for (int j = 0; j < n; ++j)
        if (j != i) pool[m++] = j;
      int ge = 0, le = 0;
      const double eps = 1e-12 * (1.0 + std::fabs(I[i]));
      for (int s = 0; s < nsim; ++s) {
        double lag = 0.0;
        for (int t = 0; t < k; ++t) {
          const int u = t + r.below(m - t);
          const int tmp = pool[t];
          pool[t] = pool[u];
          pool[u] = tmp;
          lag += w[rp[i] + t] * z[pool[t]];
        }
        const double v = z[i] / m2 * lag;
        if (v >= I[i] - eps) ++ge;
        if (v <= I[i] + eps) ++le;
      }
      PP[i] = (1.0 + (ge < le ? ge : le)) / (nsim + 1.0);
    }
  }
  return CLAB_OK;
}

int clab_locgi(const double *x, int n, const int *rp, const int *ci,
               const double *w, int star, double *out) {
  const int rc = wcheck(n, rp, ci, w);
  if (rc != CLAB_OK) return rc;
  if (n < 3 || star < 0 || star > 1) return CLAB_EARG;
  kutil::ksum a1, a2;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(x[i])) return CLAB_EDOM;
    a1.add(x[i]);
    a2.add(x[i] * x[i]);
  }
  const double t1 = a1.val(), t2 = a2.val(), dn = n;
  for (int i = 0; i < n; ++i) {
    double lag = 0.0, wi = 0.0, s1 = 0.0;
    for (int q = rp[i]; q < rp[i + 1]; ++q) {
      lag += w[q] * x[ci[q]];
      wi += w[q];
      s1 += w[q] * w[q];
    }
    double z;
    if (star) {
      // the unit in its own neighbourhood with weight 1
      lag += x[i];
      wi += 1.0;
      s1 += 1.0;
      const double xb = t1 / dn, s = std::sqrt(t2 / dn - xb * xb);
      z = (lag - wi * xb) / (s * std::sqrt((dn * s1 - wi * wi) / (dn - 1.0)));
    } else {
      const double xb = (t1 - x[i]) / (dn - 1.0);
      const double s = std::sqrt((t2 - x[i] * x[i]) / (dn - 1.0) - xb * xb);
      z = (lag - wi * xb) / (s * std::sqrt(((dn - 1.0) * s1 - wi * wi) / (dn - 2.0)));
    }
    out[i] = z;
  }
  return CLAB_OK;
}

int clab_dirmn(const double *a, const double *w, int n, double *out) {
  if (n < 1) return CLAB_EARG;
  kutil::ksum sc, ss, sw;
  for (int i = 0; i < n; ++i) {
    const double wi = w ? w[i] : 1.0;
    if (!std::isfinite(a[i]) || !kutil::nonneg(wi)) return CLAB_EDOM;
    sc.add(wi * std::cos(a[i]));
    ss.add(wi * std::sin(a[i]));
    sw.add(wi);
  }
  const double W = sw.val();
  if (!(W > 0.0)) return CLAB_EDOM;
  const double c = sc.val() / W, s = ss.val() / W;
  double r = std::hypot(c, s);
  if (r > 1.0) r = 1.0;  // rounding
  double t = std::atan2(s, c);
  if (t < 0.0) t += 2.0 * cgeo::pi;
  // a resultant at the level of rounding has no direction
  out[0] = r > 8.0 * 2.220446049250313e-16 ? t : kutil::nan();
  out[1] = r;
  out[2] = 1.0 - r;
  out[3] = std::sqrt(-2.0 * std::log(r));
  out[4] = W;
  return CLAB_OK;
}

}  // extern "C"
