// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, the dark figure of crime: capture-recapture from two or
// three lists (clab_crcap) and the bounds on a true victimisation rate
// and on the offences behind a recorded count (clab_dkfig).
#include <cmath>

#include "crimlab.h"
#include "kutil.h"

int clab_crcap(const double *x, int k, double kappa, const double *it, int ni, double *st,
               double *pw, double *sens) {
  if (k != 3 && k != 7) return CLAB_EARG;
  if (!(kappa >= 1) || !std::isfinite(kappa)) return CLAB_EDOM;
  for (int j = 0; j < k; ++j) {
    if (!kutil::nonneg(x[j])) return CLAB_EDOM;
  }
  for (int j = 0; j < ni; ++j) {
    if (!std::isfinite(it[j])) return CLAB_EDOM;
  }
  if (k == 3) {
    const double n1 = x[0], n2 = x[1], m = x[2];
    if (!(m > 0) || m > n1 || m > n2) return CLAB_EDOM;
    const double p = n1 * n2 / m;
    st[0] = p;
    st[1] = p / kappa;
    st[2] = kappa * p;
    st[3] = (n1 + 1) * (n2 + 1) / (m + 1) - 1;
    // Seber's variance of the Chapman estimate
    st[4] = std::sqrt((n1 + 1) * (n2 + 1) * (n1 - m) * (n2 - m) / ((m + 1) * (m + 1) * (m + 2)));
    st[5] = n1 + n2 - m;
    return CLAB_OK;
  }
  // the cells 100, 010, 001, 110, 101, 011, 111
  if (!(x[3] > 0 && x[4] > 0 && x[5] > 0)) return CLAB_EDOM;
  kutil::ksum o;
  for (int j = 0; j < 7; ++j) o.add(x[j]);
  const double miss = x[0] * x[1] * x[2] * x[6] / (x[3] * x[4] * x[5]);
  st[0] = o.val();
  st[1] = miss;
  st[2] = o.val() + miss;
  // the three pairs of lists, 1-2, 1-3, 2-3: each list's total, the
  // overlap, the distinct floor, Petersen and Chapman
  const int in[3][4] = {{0, 3, 4, 6}, {1, 3, 5, 6}, {2, 4, 5, 6}};
  const int ov[3][2] = {{3, 6}, {4, 6}, {5, 6}};
  const int pr[3][2] = {{0, 1}, {0, 2}, {1, 2}};
  for (int p = 0; p < 3; ++p) {
    const int *a = in[pr[p][0]], *b = in[pr[p][1]];
    const double n1 = x[a[0]] + x[a[1]] + x[a[2]] + x[a[3]];
    const double n2 = x[b[0]] + x[b[1]] + x[b[2]] + x[b[3]];
    const double m = x[ov[p][0]] + x[ov[p][1]];
    pw[p] = n1;
    pw[3 + p] = n2;
    pw[6 + p] = m;
    pw[9 + p] = n1 + n2 - m;
    pw[12 + p] = n1 * n2 / m;
    pw[15 + p] = (n1 + 1) * (n2 + 1) / (m + 1) - 1;
  }
  // a three-way log interaction t scales the missing cell by exp(-t)
  for (int j = 0; j < ni; ++j) {
    sens[j] = miss * std::exp(-it[j]);
    sens[ni + j] = o.val() + sens[j];
  }
  return CLAB_OK;
}

int clab_dkfig(const double *v, const double *r, int n, double amax, double bmax, double thr,
               double *lo, double *hi, double *brk) {
  if (n < 0) return CLAB_EARG;
  if (!kutil::nonneg(amax) || !kutil::nonneg(bmax) || !(amax + bmax < 1)) return CLAB_EDOM;
  for (int i = 0; i < n; ++i) {
    if (!(v[i] >= 0 && v[i] <= 1) || !(r[i] >= 0 && r[i] <= 1)) return CLAB_EDOM;
    if (!std::isnan(thr) && !(thr > v[i] && v[i] > 0)) return CLAB_EDOM;
  }
  for (int i = 0; i < n; ++i) {
    // the true rate: at least the recorded one and what the false
    // positives leave, at most what the false negatives can hide
    lo[i] = std::fmax(r[i], (v[i] - amax) / (1 - amax));
    hi[i] = std::fmax(lo[i], std::fmin(1.0, v[i] / (1 - bmax)));
    brk[i] = std::isnan(thr) ? kutil::nan() : 1 - v[i] / thr;
  }
  return CLAB_OK;
}

int clab_dkcmp(double va, double ra, double rb, double lam, double *out) {
  if (!kutil::pos(va) || !kutil::pos(ra) || !kutil::nonneg(rb)) return CLAB_EDOM;
  if (ra > va || !(lam > 0 && lam <= 1)) return CLAB_EDOM;
  // offences of type A are recorded at rate ra / va; type B, which the
  // survey does not cover, at no less than lam times that
  const double ka = ra / va, k = lam * ka;
  out[0] = ka;
  out[1] = (ra + rb) / (va + rb);
  out[2] = va + rb;
  out[3] = va + rb / k;
  out[4] = va - ra;
  out[5] = va - ra + rb * (1 / k - 1);
  out[6] = rb * (1 - k) / k;
  return CLAB_OK;
}

int clab_dkhier(const double *k, int n, double *out) {
  if (n < 1) return CLAB_EARG;
  kutil::ksum s;
  double mx = 0;
  for (int i = 0; i < n; ++i) {
    if (!(k[i] >= 1) || k[i] != std::floor(k[i]) || !std::isfinite(k[i])) return CLAB_EDOM;
    s.add(k[i]);
    mx = std::fmax(mx, k[i]);
  }
  out[0] = n;
  out[1] = s.val();
  out[2] = s.val() / n;
  out[3] = s.val() / n - 1;
  out[4] = n;
  out[5] = mx * n;
  return CLAB_OK;
}
