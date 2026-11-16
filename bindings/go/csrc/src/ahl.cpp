// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The burden of air pollution on health, shared by the R wrappers and the
// command line: the published concentration-response parameters by
// pollutant and outcome (Chen and Hoek 2020; Burnett et al. 2014; Huangfu
// and Atkinson 2020; Khreis et al. 2017) and the burden of each area of a
// table (Levin's attributable fraction of the whole population exposed).
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

namespace {

using std::size_t;

// A value a scalar check of polbrdn() takes: finite and not negative.
bool okv(double v) { return std::isfinite(v) && v >= 0; }

// A baseline rate per person-year inside the band (0, or 1e-5 to 0.1,
// widened to 1e-11 by rare and to 1 by hirate).
bool okrate(double r, int rare, int hirate) {
  if (!okv(r) || r > 1) return false;
  if (r == 0) return true;
  if (r < 1e-11) return false;
  if (r < 1e-5 && rare == 0) return false;
  return !(r > 0.1 && hirate == 0);
}

}  // namespace

int clab_crfpar(int pm, int out, int *form, double *par, double *ci) CLAB_TRY {
  if (form == nullptr || par == nullptr || ci == nullptr) return CLAB_EARG;
  if (out < 0 || out > 2) return CLAB_EARG;
  const double na = std::nan("");
  *form = 0;
  par[0] = 0.0;
  par[1] = 0.0;
  par[2] = 0.0;
  ci[0] = na;
  ci[1] = na;
  if (pm != 0 && out == 0) {
    par[0] = std::log(1.08);
    ci[0] = std::log(1.06);
    ci[1] = std::log(1.09);
  } else if (pm != 0) {
    // the integrated exposure-response function: alpha, gamma, delta
    *form = 1;
    const bool ihd = out == 1;
    par[0] = ihd ? 1.91 : 1.46;
    par[1] = ihd ? 0.14 : 0.13;
    par[2] = ihd ? 0.49 : 0.61;
  } else if (out == 0) {
    par[0] = std::log(1.02);
    ci[0] = std::log(1.01);
    ci[1] = std::log(1.04);
  } else if (out == 1) {
    par[0] = std::log(1.03);
    ci[0] = std::log(1.00);
    ci[1] = std::log(1.05);
  } else {
    // asthma onset in children: the RR per 4 ug/m3 put per 10
    par[0] = std::log(1.05) * 10 / 4;
    ci[0] = std::log(1.02) * 10 / 4;
    ci[1] = std::log(1.07) * 10 / 4;
  }
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_pbarea(std::size_t n, const double *e, const double *pop, const double *br, int pm, int out,
                int rare, int hirate, double *rr, double *paf, double *att, double *base, std::int64_t *ord,
                std::int64_t *bad) CLAB_TRY {
  if (bad == nullptr) return CLAB_EARG;
  *bad = 0;
  if (n > 0 && (e == nullptr || pop == nullptr || br == nullptr || rr == nullptr || paf == nullptr ||
                att == nullptr || base == nullptr || ord == nullptr)) {
    return CLAB_EARG;
  }
  int form = 0;
  double par[3] = {0.0, 0.0, 0.0};
  double ci[2] = {0.0, 0.0};
  const int rc = clab_crfpar(pm, out, &form, par, ci);
  if (rc != CLAB_OK) return rc;
  const double refc = pm != 0 ? 5.8 : 10.0;
  for (size_t i = 0; i < n; ++i) {
    if (!okv(e[i]) || !okrate(br[i], rare, hirate) || !okv(pop[i])) {
      *bad = static_cast<std::int64_t>(i) + 1;
      return CLAB_OK;
    }
  }
  std::vector<double> a(n, 0.0), b(n, 0.0), c(n, 0.0), d(n, 0.0);
  for (size_t i = 0; i < n; ++i) {
    double r = 0.0;
    double lr = 0.0;
    const int k = clab_crf(&e[i], 1, refc, form, par, &r, &lr);
    if (k != CLAB_OK) return k;
    const double bs = br[i] * pop[i];
    a[i] = r;
    b[i] = clab_levin(r, 1.0);
    c[i] = b[i] * bs;
    d[i] = bs;
  }
  // worst first; ties keep their order and NaN goes last, as R's order()
  std::vector<size_t> o(n, 0);
  std::iota(o.begin(), o.end(), size_t{0});
  std::stable_sort(o.begin(), o.end(), [&c](size_t x, size_t y) {
    if (std::isnan(c[x])) return false;
    return std::isnan(c[y]) || c[x] > c[y];
  });
  for (size_t i = 0; i < n; ++i) {
    const size_t j = o[i];
    rr[i] = a[j];
    paf[i] = b[j];
    att[i] = c[j];
    base[i] = d[j];
    ord[i] = static_cast<std::int64_t>(j) + 1;
  }
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)
