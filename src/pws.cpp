// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Power by simulation (crimlab.h, pwrsim()): a two-arm trial is drawn
// nsim times for each design (control and treated sizes, effect) and
// tested; the power is the share of rejections at alpha, with its Monte
// Carlo standard error and Wilson's interval (Morris, White and Crowther
// 2019). Outcomes are drawn by inversion from uniforms of the package's
// generator: normal (qnorm), binary (u < p) or Poisson (qpois). The test
// is Welch's t or the permutation test of the difference in means. Each
// design has its own stream, made from the seed and the design's number,
// so that a design's power does not depend on the others asked for.
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

namespace {

using std::size_t;

// A uniform strictly inside (0, 1).
double uo(kutil::rng &r) { return (static_cast<double>(r.next() >> 11) + 0.5) * 0x1.0p-53; }

double draw(kutil::rng &r, int fam, double par, double eff, bool trt) {
  const double u = uo(r);
  if (fam == 0) return (trt ? eff : 0.0) + par * clab_qnorm(u);
  if (fam == 1) return u < par + (trt ? eff : 0.0) ? 1.0 : 0.0;
  return clab_qpois(u, par + (trt ? eff : 0.0), 1);
}

// Welch's test of y[0..n1) (treated) against y[n1..n1 + n0); with no
// spread in either arm, p is 0 when the means differ in the direction
// tested and 1 otherwise.
double welch(const double *y, int n1, int n0, int alt) {
  double m1 = 0, m0 = 0, v1 = 0, v0 = 0;
  for (int i = 0; i < n1; ++i) m1 += y[i];
  for (int i = 0; i < n0; ++i) m0 += y[n1 + i];
  m1 /= n1;
  m0 /= n0;
  for (int i = 0; i < n1; ++i) v1 += (y[i] - m1) * (y[i] - m1);
  for (int i = 0; i < n0; ++i) v0 += (y[n1 + i] - m0) * (y[n1 + i] - m0);
  const double a = v1 / (n1 - 1) / n1, b = v0 / (n0 - 1) / n0, se = std::sqrt(a + b), d = m1 - m0;
  if (!(se > 0)) {
    if (alt == 0) return d != 0 ? 0.0 : 1.0;
    return (alt == 1 ? d > 0 : d < 0) ? 0.0 : 1.0;
  }
  const double t = d / se, df = (a + b) * (a + b) / (a * a / (n1 - 1) + b * b / (n0 - 1));
  if (alt == 0) return 2 * clab_pt(-std::fabs(t), df, 1);
  return alt == 1 ? clab_pt(t, df, 0) : clab_pt(t, df, 1);
}

// The permutation test of the difference in means, nperm relabellings.
double perm(kutil::rng &r, std::vector<double> &y, int n1, int n0, int alt, int nperm) {
  const int n = n1 + n0;
  double s = 0, s1 = 0, ym = 0;
  for (int i = 0; i < n; ++i) {
    s += y[static_cast<size_t>(i)];
    ym = std::fmax(ym, std::fabs(y[static_cast<size_t>(i)]));
  }
  for (int i = 0; i < n1; ++i) s1 += y[static_cast<size_t>(i)];
  auto tv = [&](double t1) { return t1 / n1 - (s - t1) / n0; };
  const double obs = tv(s1), tol = 1e-9 * (std::fabs(obs) + ym);
  double k = 0;
  for (int q = 0; q < nperm; ++q) {
    double t1 = 0;
    for (int i = 0; i < n1; ++i) {
      const int j = i + r.below(n - i);
      const double v = y[static_cast<size_t>(j)];
      y[static_cast<size_t>(j)] = y[static_cast<size_t>(i)];
      y[static_cast<size_t>(i)] = v;
      t1 += v;
    }
    const double t = tv(t1);
    if (alt == 0 ? std::fabs(t) >= std::fabs(obs) - tol : alt == 1 ? t >= obs - tol : t <= obs + tol) k += 1;
  }
  return (1 + k) / (1 + nperm);
}

}  // namespace

extern "C" int clab_pwrsim(int fam, double par, int m, const int *n0, const int *n1, const double *eff,
                           int test, int alt, double alpha, int nsim, int nperm, double level,
                           unsigned int seed, double *out) CLAB_TRY {
  if (fam < 0 || fam > 2 || test < 0 || test > 1 || alt < 0 || alt > 2 || m < 0 || nsim < 1 ||
      (test == 1 && nperm < 1) || !(alpha > 0 && alpha < 1) || !(level > 0 && level < 1)) {
    return CLAB_EARG;
  }
  const double z = clab_qnorm(0.5 + level / 2), z2 = z * z, ns = nsim;
  std::vector<double> y;
  for (int c = 0; c < m; ++c) {
    if (n0[c] < 2 || n1[c] < 2) return CLAB_EARG;
    const double p1 = par + eff[c];
    if (fam == 0 && !(par > 0)) return CLAB_EARG;
    if (fam == 1 && !(par >= 0 && par <= 1 && p1 >= 0 && p1 <= 1)) return CLAB_EDOM;
    if (fam == 2 && !(par > 0 && p1 >= 0)) return CLAB_EDOM;
    kutil::rng r((static_cast<std::uint64_t>(seed) << 32) | static_cast<std::uint64_t>(c));
    y.assign(static_cast<size_t>(n1[c]) + static_cast<size_t>(n0[c]), 0.0);
    double rej = 0;
    for (int s = 0; s < nsim; ++s) {
      if (kutil::stop(static_cast<double>(y.size()) * (test == 1 ? nperm + 1.0 : 1.0))) return CLAB_EINTR;
      for (int i = 0; i < n1[c]; ++i) y[static_cast<size_t>(i)] = draw(r, fam, par, eff[c], true);
      for (int i = 0; i < n0[c]; ++i) y[static_cast<size_t>(n1[c] + i)] = draw(r, fam, par, eff[c], false);
      const double p = test == 0 ? welch(y.data(), n1[c], n0[c], alt) : perm(r, y, n1[c], n0[c], alt, nperm);
      if (p <= alpha) rej += 1;
    }
    const double pw = rej / ns, se = std::sqrt(pw * (1 - pw) / ns);
    const double mid = (pw + z2 / (2 * ns)) / (1 + z2 / ns),
                 hw = z * std::sqrt(pw * (1 - pw) / ns + z2 / (4 * ns * ns)) / (1 + z2 / ns);
    out[c] = pw;
    out[m + c] = se;
    out[2 * m + c] = mid - hw;
    out[3 * m + c] = mid + hw;
    out[4 * m + c] = rej;
  }
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)
