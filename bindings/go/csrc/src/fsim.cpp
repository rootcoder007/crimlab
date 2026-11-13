// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, simulation for the predictive-policing audits (after Barman
// and Barman, arXiv:2603.18987): the noisy-OR chance that patrols detect an
// incident (clab_nordet) and records of incidents whose recording and
// predicted risk differ by group by a known amount (clab_fairsim).
#include <cmath>

#include "crimlab.h"
#include "kutil.h"

namespace {

// A standard normal by the Box-Muller transform (cosine half), from two
// uniforms; 1 - u keeps the logarithm's argument in (0, 1].
double rnorm(kutil::rng &g) {
  const double u = 1.0 - g.unif(), v = g.unif();
  return std::sqrt(-2.0 * std::log(u)) * std::cos(6.283185307179586 * v);
}

}  // namespace

extern "C" {

int clab_nordet(const double *cx, const double *cy, int n, const double *ox,
                const double *oy, int m, double radius, double p, int draw,
                unsigned int seed, int *k, double *prob, int *det) {
  if (n < 1 || m < 0) return CLAB_EARG;
  if (!kutil::pos(radius) || !(p > 0.0 && p <= 1.0)) return CLAB_EDOM;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(cx[i]) || !std::isfinite(cy[i])) return CLAB_EDOM;
  }
  for (int j = 0; j < m; ++j) {
    if (!std::isfinite(ox[j]) || !std::isfinite(oy[j])) return CLAB_EDOM;
  }
  // ponytail: every pair, O(n m); a grid of cells of side radius would
  // make it O(n + m) per occupied cell when the counts grow large
  const double lq = p < 1.0 ? std::log1p(-p) : 0.0;
  for (int i = 0; i < n; ++i) {
    int c = 0;
    for (int j = 0; j < m; ++j) {
      const double dx = cx[i] - ox[j], dy = cy[i] - oy[j];
      c += std::sqrt(dx * dx + dy * dy) <= radius ? 1 : 0;
    }
    k[i] = c;
    // 1 - (1 - p)^k: each officer in range detects independently
    prob[i] = c == 0 ? 0.0 : (p < 1.0 ? -std::expm1(c * lq) : 1.0);
  }
  if (draw != 0) {
    kutil::rng g(seed);
    for (int i = 0; i < n; ++i) det[i] = g.unif() < prob[i] ? 1 : 0;
  }
  return CLAB_OK;
}

int clab_fairsim(int n, int g, const double *props, int nareas, double base,
                 double bias, unsigned int seed, int *area, int *grp,
                 int *truth, int *det, double *risk) {
  if (n < 1 || g < 2 || nareas < g) return CLAB_EARG;
  if (!(base >= 0.0 && base <= 1.0) || !(bias >= -1.0 && bias <= 1.0)) return CLAB_EDOM;
  double tot = 0.0;
  if (props != nullptr) {
    for (int j = 0; j < g; ++j) {
      if (!kutil::nonneg(props[j])) return CLAB_EDOM;
      tot += props[j];
    }
    if (!(tot > 0.0) || !std::isfinite(tot)) return CLAB_EDOM;
  }
  const double dr = std::fmin(1.0, std::fmax(0.0, base * (1.0 - bias)));
  kutil::rng gen(seed);
  for (int i = 0; i < n; ++i) {
    // the group, by its share
    int j = 0;
    if (props == nullptr) {
      j = gen.below(g);
    } else {
      double u = gen.unif() * tot, c = props[0];
      while (j < g - 1 && u >= c) c += props[++j];
    }
    grp[i] = j;
    // areas are segregated: area a holds group a mod g only
    const int na = (nareas - 1 - j) / g + 1;
    area[i] = j + g * gen.below(na);
    // the true outcome does not depend on the group; its recording does,
    // and so does the predicted risk: the injected disparity
    truth[i] = gen.unif() < base ? 1 : 0;
    det[i] = gen.unif() < (j == 0 ? base : dr) ? 1 : 0;
    const double s = 250.0 + (j == 0 ? 0.0 : 100.0 * bias) + 40.0 * rnorm(gen);
    risk[i] = std::fmin(500.0, std::fmax(0.0, s));
  }
  return CLAB_OK;
}

}  // extern "C"
