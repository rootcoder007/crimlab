// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, rates: rates and proportions by group with their ratios and
// differences against a reference group (clab_grprr), directly
// standardised rates (clab_dastd) and indirect standardisation with
// empirical Bayes shrinkage (clab_sirex). The exact intervals need the
// quantiles of the gamma and beta distributions; the binding computes
// them from the counts these kernels return.
#include <cmath>
#include <limits>

#include "crimlab.h"

namespace {

constexpr double qnan = std::numeric_limits<double>::quiet_NaN();
constexpr double inf = std::numeric_limits<double>::infinity();

// The Wilson score interval of x events in n trials (Wilson 1927).
void wilson(double x, double n, double z, double *lo, double *hi) {
  const double p = x / n, z2 = z * z, d = 1 + z2 / n;
  const double c = (p + z2 / (2 * n)) / d;
  const double h = z * std::sqrt(p * (1 - p) / n + z2 / (4 * n * n)) / d;
  *lo = std::fmax(0.0, c - h);
  *hi = std::fmin(1.0, c + h);
}

bool cnts(const double *x, const double *n, int k, bool prop) {
  for (int i = 0; i < k; ++i) {
    if (!(x[i] >= 0) || !std::isfinite(x[i]) || !(n[i] > 0) || !std::isfinite(n[i])) {
      return false;
    }
    if (prop && x[i] > n[i]) return false;
  }
  return true;
}

}  // namespace

extern "C" {

int clab_grprr(const double *x, const double *n, int k, int ref, int prop, double z,
               double *rate, double *lo, double *hi, double *rep, double *ratio,
               double *rlo, double *rhi, double *diff, double *dlo, double *dhi) {
  if (k < 1 || ref < 0 || ref >= k) return CLAB_EARG;
  if (!(z > 0) || !std::isfinite(z) || !cnts(x, n, k, prop != 0)) return CLAB_EDOM;
  double xs = 0, ns = 0;
  for (int i = 0; i < k; ++i) {
    xs += x[i];
    ns += n[i];
  }
  for (int i = 0; i < k; ++i) {
    rate[i] = x[i] / n[i];
    if (prop) {
      wilson(x[i], n[i], z, &lo[i], &hi[i]);
    } else if (x[i] > 0) {
      // the log-normal interval of a Poisson rate
      const double s = z / std::sqrt(x[i]);
      lo[i] = rate[i] * std::exp(-s);
      hi[i] = rate[i] * std::exp(s);
    } else {
      lo[i] = hi[i] = qnan;
    }
    // the group's share of the events over its share of the denominator
    rep[i] = xs > 0 ? (x[i] / xs) / (n[i] / ns) : qnan;
  }
  const double x0 = x[ref], n0 = n[ref], r0 = rate[ref];
  for (int i = 0; i < k; ++i) {
    const double x1 = x[i], n1 = n[i], r1 = rate[i];
    ratio[i] = r0 > 0 ? r1 / r0 : (r1 > 0 ? inf : qnan);
    diff[i] = r1 - r0;
    // log-normal intervals of the ratio: Katz et al. (1978) for risks,
    // 1 / x1 + 1 / x0 for Poisson rates; none with a zero count
    if (i == ref || !(x1 > 0) || !(x0 > 0)) {
      rlo[i] = rhi[i] = qnan;
    } else {
      const double v = prop ? 1 / x1 - 1 / n1 + 1 / x0 - 1 / n0 : 1 / x1 + 1 / x0;
      const double s = z * std::sqrt(v);
      rlo[i] = ratio[i] * std::exp(-s);
      rhi[i] = ratio[i] * std::exp(s);
    }
    if (i == ref) {
      dlo[i] = dhi[i] = qnan;
    } else if (prop) {
      // Newcombe's (1998) hybrid score interval from the Wilson limits
      const double a = std::fmax(r1 - lo[i], 0.0), b = std::fmax(hi[i] - r1, 0.0);
      const double c = std::fmax(r0 - lo[ref], 0.0), d = std::fmax(hi[ref] - r0, 0.0);
      dlo[i] = diff[i] - std::sqrt(a * a + d * d);
      dhi[i] = diff[i] + std::sqrt(b * b + c * c);
    } else {
      const double s = z * std::sqrt(x1 / (n1 * n1) + x0 / (n0 * n0));
      dlo[i] = diff[i] - s;
      dhi[i] = diff[i] + s;
    }
  }
  return CLAB_OK;
}

int clab_dastd(const double *x, const double *n, const int *grp, const int *str, int m,
               int ng, int ns, const double *w, int ref, double z, double *rate,
               double *var, double *wmax, double *ev, double *ratio, double *rlo,
               double *rhi, double *diff, double *dlo, double *dhi, int *seen) {
  if (m < 1 || ng < 1 || ns < 1 || ref < 0 || ref >= ng) return CLAB_EARG;
  if (!(z > 0) || !std::isfinite(z) || !cnts(x, n, m, false)) return CLAB_EDOM;
  double ws = 0;
  for (int s = 0; s < ns; ++s) {
    if (!(w[s] >= 0) || !std::isfinite(w[s])) return CLAB_EDOM;
    ws += w[s];
  }
  if (!(ws > 0)) return CLAB_EDOM;
  for (int g = 0; g < ng; ++g) rate[g] = var[g] = wmax[g] = ev[g] = 0;
  for (int q = 0; q < ng * ns; ++q) seen[q] = 0;
  for (int i = 0; i < m; ++i) {
    const int g = grp[i], s = str[i];
    if (g < 0 || g >= ng || s < 0 || s >= ns) return CLAB_EARG;
    // one row per group and stratum
    if (seen[g * ns + s]++) return CLAB_EARG;
    const double wi = w[s] / ws / n[i];
    rate[g] += wi * x[i];
    var[g] += wi * wi * x[i];
    wmax[g] = std::fmax(wmax[g], wi);
    ev[g] += x[i];
  }
  // every group needs every stratum that carries weight
  for (int g = 0; g < ng; ++g) {
    for (int s = 0; s < ns; ++s) {
      if (w[s] > 0 && !seen[g * ns + s]) return CLAB_EDOM;
    }
  }
  const double y0 = rate[ref];
  for (int g = 0; g < ng; ++g) {
    const double y1 = rate[g];
    ratio[g] = y0 > 0 ? y1 / y0 : (y1 > 0 ? inf : qnan);
    diff[g] = y1 - y0;
    if (g == ref) {
      rlo[g] = rhi[g] = dlo[g] = dhi[g] = qnan;
      continue;
    }
    if (y1 > 0 && y0 > 0) {
      const double s = z * std::sqrt(var[g] / (y1 * y1) + var[ref] / (y0 * y0));
      rlo[g] = ratio[g] * std::exp(-s);
      rhi[g] = ratio[g] * std::exp(s);
    } else {
      rlo[g] = rhi[g] = qnan;
    }
    const double s = z * std::sqrt(var[g] + var[ref]);
    dlo[g] = diff[g] - s;
    dhi[g] = diff[g] + s;
  }
  return CLAB_OK;
}

int clab_sirex(const double *x, const double *n, const int *grp, const int *str, int m,
               int ng, int ns, const double *sr, double *obs, double *expd, double *sir,
               double *eb, double *shr, double *prior, double *work) {
  if (m < 1 || ng < 1 || ns < 1) return CLAB_EARG;
  if (!cnts(x, n, m, false)) return CLAB_EDOM;
  double *num = work, *den = work + ns;
  for (int s = 0; s < ns; ++s) num[s] = den[s] = 0;
  for (int g = 0; g < ng; ++g) obs[g] = expd[g] = 0;
  for (int i = 0; i < m; ++i) {
    if (grp[i] < 0 || grp[i] >= ng || str[i] < 0 || str[i] >= ns) return CLAB_EARG;
    num[str[i]] += x[i];
    den[str[i]] += n[i];
  }
  // the stratum rates: given, or pooled over every group
  for (int s = 0; s < ns; ++s) {
    const double r = sr ? sr[s] : num[s] / den[s];
    if (!(r >= 0) || !std::isfinite(r)) return CLAB_EDOM;
    num[s] = r;
  }
  for (int i = 0; i < m; ++i) {
    obs[grp[i]] += x[i];
    expd[grp[i]] += n[i] * num[str[i]];
  }
  double os = 0, es = 0;
  for (int g = 0; g < ng; ++g) {
    if (!(expd[g] > 0)) return CLAB_EDOM;
    sir[g] = obs[g] / expd[g];
    os += obs[g];
    es += expd[g];
  }
  // Marshall's (1991) moment estimates of a gamma prior on the relative
  // risks: mean mm, variance s2 net of the Poisson noise
  const double mm = os / es;
  double s2 = 0;
  for (int g = 0; g < ng; ++g) s2 += expd[g] * (sir[g] - mm) * (sir[g] - mm);
  s2 = s2 / es - mm / (es / ng);
  if (!(s2 > 0)) {
    for (int g = 0; g < ng; ++g) {
      eb[g] = mm;
      shr[g] = 1;
    }
    prior[0] = mm;
    prior[1] = 0;
    return CLAB_OK;
  }
  const double al = mm / s2, nu = mm * mm / s2;
  for (int g = 0; g < ng; ++g) {
    eb[g] = (obs[g] + nu) / (expd[g] + al);
    shr[g] = al / (expd[g] + al);
  }
  prior[0] = mm;
  prior[1] = s2;
  return CLAB_OK;
}

}  // extern "C"
