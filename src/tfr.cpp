// SPDX-License-Identifier: AGPL-3.0-or-later
//
// libcrimlab, forensics: the distribution of the number of glass fragments
// recovered under the transfer, persistence and recovery model of Curran
// et al. (1998). Conditional on the per-case random effects the count is
// Poisson; the transfer intensity (normal transfer amount, gamma distance)
// is integrated by deterministic quadrature, the time since the offence is
// summed exactly, and only the hourly loss rates and the recovery rate,
// whose dimension grows with time, are averaged by Monte Carlo.
#include <cstdint>

#include "crimlab.h"
#include "kutil.h"

namespace {

// 8-point Gauss-Legendre rule on [-1, 1]
const double GX[4] = {0.1834346424956498, 0.5255324099163290, 0.7966664774136267,
                      0.9602898564975363};
const double GW[4] = {0.3626837833783620, 0.3137066458778873, 0.2223810344533745,
                      0.1012285362903763};

// composite 8-point rule on [a, b] with np panels: node k, weight k
void glnode(double a, double b, int np, int k, double *x, double *w) {
  const int panel = k / 8, j = k % 8;
  const double h = (b - a) / np, c = a + (panel + 0.5) * h;
  const double t = j < 4 ? -GX[3 - j] : GX[j - 4];
  const double v = j < 4 ? GW[3 - j] : GW[j - 4];
  *x = c + 0.5 * h * t;
  *w = 0.5 * h * v;
}

struct rng {
  std::uint64_t s;
  explicit rng(std::uint64_t seed) {
    std::uint64_t z = seed + 0x9E3779B97F4A7C15ULL;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    s = (z ^ (z >> 31)) | 1ULL;
  }
  double unif() {
    s ^= s >> 12;
    s ^= s << 25;
    s ^= s >> 27;
    return static_cast<double>((s * 2685821657736338717ULL) >> 11) * 0x1p-53;
  }
};

enum { LAM, LSD, DIS, Q, L0, U0, LS0, US0, LJ, UJ, LSJ, USJ, LR, UR, LT, UT, RR, NPAR };

const double TEPS = 1e-16;  // negative-binomial tail left out

// the time distribution in hours: ph[h] for h = 1..hmax (t = 0 and t = 1
// both give one hour of loss); returns hmax, or 0 on a bad parameter, or
// with ph == NULL only the length. *cut: the probability left out.
int hours(const double *p, int tdist, double *ph, double *cut) {
  const double lt = p[LT], ut = p[UT];
  *cut = 0.0;
  if (lt == ut) {
    const int t = static_cast<int>(lt);
    const int h = t < 1 ? 1 : t;
    if (ph != nullptr) {
      for (int k = 1; k < h; ++k) ph[k] = 0.0;
      ph[h] = 1.0;
    }
    return h;
  }
  const int t0 = static_cast<int>(std::floor(lt)), t1 = static_cast<int>(std::ceil(ut));
  const double size = 0.5 * (lt + ut), r = p[RR];
  if (tdist == 0) {
    // NB(size, r) on 0, 1, 2, ...: P(0) = r^size, P(k) = P(k-1) (k - 1 + size)(1 - r) / k
    double pk = std::exp(size * std::log(r)), cum = 0.0;
    int k = 0;
    const double mean = size * (1.0 - r) / r;
    for (;; ++k) {
      if (ph != nullptr) {
        const int h = k < 1 ? 1 : k;
        ph[h] += pk;
      }
      cum += pk;
      // past the mean the ratio q = P(k+1)/P(k) = (k + size)(1 - r)/(k + 1)
      // falls, so the tail beyond k is at most P(k) q / (1 - q)
      const double q = (k + size) * (1.0 - r) / (k + 1);
      if (k >= mean && q < 1.0 && pk * q / (1.0 - q) <= TEPS) {
        *cut = pk * q / (1.0 - q);
        break;
      }
      if (k > 100000) {
        *cut = 1.0 - cum > 0.0 ? 1.0 - cum : 0.0;
        break;
      }
      pk *= q;
    }
    return k < 1 ? 1 : k;
  }
  // cnegbin and uniform: times t0..t1
  double tot = 0.0;
  const int hmax = t1 < 1 ? 1 : t1;
  if (ph != nullptr)
    for (int h = 1; h <= hmax; ++h) ph[h] = 0.0;
  double pk = std::exp(size * std::log(r));
  for (int k = 0; k <= t1; ++k) {
    if (k >= t0) {
      const double v = tdist == 1 ? pk : 1.0;
      tot += v;
      if (ph != nullptr) ph[k < 1 ? 1 : k] += v;
    }
    pk *= (k + size) * (1.0 - r) / (k + 1);
  }
  if (ph != nullptr)
    for (int h = 1; h <= hmax; ++h) ph[h] /= tot;
  return hmax;
}

int check(const double *p, int deffect, int tdist, int ymax) {
  for (int k = 0; k < NPAR; ++k)
    if (!std::isfinite(p[k])) return CLAB_EDOM;
  if (ymax < 0 || tdist < 0 || tdist > 2 || deffect < 0 || deffect > 1) return CLAB_EARG;
  if (!(p[LAM] >= 0.0) || !(p[LSD] >= 0.0)) return CLAB_EDOM;
  if (deffect && !(p[DIS] > 0.0)) return CLAB_EDOM;
  if (!(p[Q] >= 0.0 && p[Q] <= 1.0)) return CLAB_EDOM;
  for (int k = L0; k <= LR; k += 2)
    if (!(p[k] >= 0.0 && p[k] <= p[k + 1] && p[k + 1] <= 1.0)) return CLAB_EDOM;
  if (!(p[LT] >= 0.0 && p[LT] <= p[UT])) return CLAB_EDOM;
  if (p[LT] == p[UT] && p[LT] != std::floor(p[LT])) return CLAB_EDOM;
  if (p[LT] != p[UT] && !(p[RR] > 0.0 && p[RR] <= 1.0)) return CLAB_EDOM;
  return CLAB_OK;
}

// density of log|X| at z, X ~ N(mu, sd^2), sd > 0
inline double flx(double z, double mu, double sd) {
  const double x = std::exp(z), z1 = (x - mu) / sd, z2 = (x + mu) / sd;
  return x * (std::exp(-0.5 * z1 * z1) + std::exp(-0.5 * z2 * z2)) / (sd * 2.5066282746310002);
}

const double LCUT = -39.14394658089878;  // log(1e-17): below, e^l s < 1e-17

// The transfer intensity M = |X| e^{1 - D/d} (X ~ N(lam, lsd^2), D ~
// Gamma(d, 1); M = |X| without the distance effect) as weighted nodes
// (m, w). With lsd > 0, L = log M has a smooth density f_L, tabulated on
// the grid l = LCUT, LCUT + hl, ..., up to log(|lam| + 12 lsd) + 1, and the
// nodes carry the trapezoid weights hl f_L(l) (the trapezoid rule converges
// geometrically for smooth integrands that vanish at the ends; below LCUT
// the Poisson law is a point mass at 0 to 1e-17, handled by the caller's
// expm1 form). With the distance effect f_L(l) = int f_log|X|(l - 1 + D/d)
// f_D(D) dD, taken with D = w^{1/d} (which removes the D^{d-1} singularity:
// f_D dD = e^{-D} dw / Gamma(d + 1)) on w in [0, Dmax^d], Dmax = d + 60 +
// 10 sqrt(d), by a composite 8-point Gauss-Legendre rule with nq panels.
// With lsd = 0 the nodes are that w rule itself (or one node). Returns the
// number of nodes; with m == NULL only counts.
int nodes(const double *p, int deffect, double hl, int nq, double *m, double *w) {
  const double mu = p[LAM], sd = p[LSD], d = p[DIS];
  const double dmax = d + 60.0 + 10.0 * std::sqrt(d);
  const double wmax = std::pow(dmax, d), lgd1 = std::lgamma(d + 1.0);
  if (sd == 0.0) {
    if (!deffect) {
      if (m != nullptr) {
        m[0] = std::fabs(mu);
        w[0] = 1.0;
      }
      return 1;
    }
    if (m != nullptr) {
      for (int k = 0; k < 8 * nq; ++k) {
        double ww, wt;
        glnode(0.0, wmax, nq, k, &ww, &wt);
        const double dd = std::pow(ww, 1.0 / d);
        m[k] = std::fabs(mu) * std::exp(1.0 - dd / d);
        w[k] = wt * std::exp(-dd - lgd1);
      }
    }
    return 8 * nq;
  }
  const double lhi = std::log(std::fabs(mu) + 12.0 * sd) + (deffect ? 1.0 : 0.0);
  const int g = static_cast<int>(std::ceil((lhi - LCUT) / hl)) + 1;
  if (m == nullptr) return g;
  for (int k = 0; k < g; ++k) {
    const double l = LCUT + k * hl;
    double f;
    if (!deffect) {
      f = flx(l, mu, sd);
    } else {
      kutil::ksum s;
      for (int j = 0; j < 8 * nq; ++j) {
        double ww, wt;
        glnode(0.0, wmax, nq, j, &ww, &wt);
        const double dd = std::pow(ww, 1.0 / d);
        s.add(wt * flx(l - 1.0 + dd / d, mu, sd) * std::exp(-dd - lgd1));
      }
      f = s.val();
    }
    m[k] = std::exp(l);
    w[k] = hl * f * (k == 0 || k == g - 1 ? 0.5 : 1.0);
  }
  return g;
}

}  // namespace

extern "C" {

int clab_tfrw(const double *par, int deffect, int tdist, int ymax, double hl,
              int nq, int *len) {
  const int rc = check(par, deffect, tdist, ymax);
  if (rc != CLAB_OK) return rc;
  if (nq < 1 || nq > 100000) return CLAB_EARG;
  if (!(hl >= 1e-4 && hl <= 1.0)) return CLAB_EDOM;
  double cut;
  const int hmax = hours(par, tdist, nullptr, &cut);
  const long long k = nodes(par, deffect, hl, nq, nullptr, nullptr);
  const long long n = 2 * k + (hmax + 1) + 2 * (hmax + 1) + 3LL * (ymax + 1);
  if (n > 2147483647LL) return CLAB_EARG;
  *len = static_cast<int>(n);
  return CLAB_OK;
}

int clab_tfrp(const double *par, int deffect, int tdist, int ymax, int nsim,
              unsigned int seed, double hl, int nq, double *prob, double *se,
              double *cut, double *work) {
  int len;
  int rc = clab_tfrw(par, deffect, tdist, ymax, hl, nq, &len);
  if (rc != CLAB_OK) return rc;
  if (nsim < 1) return CLAB_EARG;
  double tcut;
  const int hmax = hours(par, tdist, nullptr, &tcut);
  const int kk = nodes(par, deffect, hl, nq, nullptr, nullptr);
  double *m = work;
  double *w = m + kk, *ph = w + kk, *ah = ph + (hmax + 1), *bh = ah + (hmax + 1);
  double *acc = bh + (hmax + 1), *sum = acc + (ymax + 1), *sq = sum + (ymax + 1);
  for (int h = 0; h <= hmax; ++h) ph[h] = 0.0;
  hours(par, tdist, ph, &tcut);
  int k = nodes(par, deffect, hl, nq, m, w);
  double wsum = 0.0;
  for (int j = 0; j < k; ++j) wsum += w[j];
  // nodes of negligible weight (below 1e-20 of the total) are left out
  {
    int k2 = 0;
    for (int j = 0; j < k; ++j) {
      if (w[j] < 1e-20 * wsum) continue;
      m[k2] = m[j];
      w[k2] = w[j];
      ++k2;
    }
    k = k2;
  }
  for (int y = 0; y <= ymax; ++y) sum[y] = sq[y] = 0.0;
  // fixed rates and recovery: one replicate is exact
  bool fixed = true;
  for (int j = L0; j <= LR; j += 2) fixed = fixed && par[j] == par[j + 1];
  const int reps = fixed ? 1 : nsim;
  rng r(seed);
  const double qq = par[Q];
  for (int rep = 0; rep < reps; ++rep) {
    // draw order: recovery, hour-1 rates (low, high), then hours 2..hmax
    const double rec = par[LR] + (par[UR] - par[LR]) * r.unif();
    double a = 1.0, b = 1.0;
    for (int h = 1; h <= hmax; ++h) {
      const double la = h == 1 ? par[L0] : par[LJ], ua = h == 1 ? par[U0] : par[UJ];
      const double lb = h == 1 ? par[LS0] : par[LSJ], ub = h == 1 ? par[US0] : par[USJ];
      a *= 1.0 - (la + (ua - la) * r.unif());
      b *= 1.0 - (lb + (ub - lb) * r.unif());
      ah[h] = a;
      bh[h] = b;
    }
    for (int y = 0; y <= ymax; ++y) acc[y] = 0.0;
    // P(0) = 1 + E[e^-V - 1] (expm1: exact for small V, and the part of
    // the intensity left below the grid adds nothing to it); P(y > 0) =
    // E[e^-V V^y / y!]
    double z0 = 0.0;
    for (int h = 1; h <= hmax; ++h) {
      if (!(ph[h] > 0.0)) continue;
      const double s = rec * ((1.0 - qq) * ah[h] + qq * bh[h]);
      for (int j = 0; j < k; ++j) {
        const double v = m[j] * s, wt = ph[h] * w[j];
        z0 += wt * std::expm1(-v);
        if (v < 600.0) {
          double pk = std::exp(-v);
          for (int y = 1; y <= ymax; ++y) {
            pk *= v / y;
            // past the mode the terms only fall: stop once they are below
            // the smallest normal double (their sum is then < 1e-290)
            if (y > v && pk < 1e-300) break;
            acc[y] += wt * pk;
          }
        } else {
          const double lv = std::log(v);
          double lf = 0.0;
          for (int y = 1; y <= ymax; ++y) {
            lf += std::log(static_cast<double>(y));
            acc[y] += wt * std::exp(-v + y * lv - lf);
          }
        }
      }
    }
    acc[0] = 1.0 + z0;
    for (int y = 0; y <= ymax; ++y) {
      sum[y] += acc[y];
      sq[y] += acc[y] * acc[y];
    }
  }
  for (int y = 0; y <= ymax; ++y) {
    prob[y] = sum[y] / reps;
    if (reps > 1) {
      const double var = (sq[y] - reps * prob[y] * prob[y]) / (reps - 1);
      se[y] = var > 0.0 ? std::sqrt(var / reps) : 0.0;
    } else {
      se[y] = 0.0;
    }
  }
  cut[0] = tcut;          // time distribution left out
  cut[1] = wsum;  // the intensity mass on the nodes (below the grid: zeros)
  return CLAB_OK;
}

}  // extern "C"
