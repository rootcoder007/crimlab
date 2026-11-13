// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, the life course and the justice pipeline: the age-crime
// curve of a mixture of latent types (clab_agemix) and its panel test
// (clab_agepan), the cumulative risk of a first event by age
// (clab_cumrsk), and rates between stages of the pipeline with a lag
// (clab_jpipe).
#include <cmath>
#include <cstdint>

#include "crimlab.h"
#include "kutil.h"

int clab_agemix(const double *cv, int na, int nt, const double *w, double *agg, double *cov) {
  if (na < 1 || nt < 1) return CLAB_EARG;
  kutil::ksum sw;
  for (int k = 0; k < nt; ++k) {
    if (!kutil::nonneg(w[k])) return CLAB_EDOM;
    sw.add(w[k]);
  }
  if (std::fabs(sw.val() - 1) > 1e-8) return CLAB_EDOM;
  for (long long c = 0; c < static_cast<long long>(na) * nt; ++c) {
    if (!std::isfinite(cv[c])) return CLAB_EDOM;
  }
  // the aggregate curve, then the covariance across types of the curves
  // at two ages: what a single-type model forces to 0
  for (int a = 0; a < na; ++a) {
    kutil::ksum s;
    for (int k = 0; k < nt; ++k) s.add(w[k] * cv[a + static_cast<long long>(na) * k]);
    agg[a] = s.val();
  }
  for (int a = 0; a < na; ++a) {
    for (int b = a; b < na; ++b) {
      kutil::ksum s;
      for (int k = 0; k < nt; ++k) {
        s.add(w[k] * (cv[a + static_cast<long long>(na) * k] - agg[a]) *
              (cv[b + static_cast<long long>(na) * k] - agg[b]));
      }
      cov[a + static_cast<long long>(na) * b] = cov[b + static_cast<long long>(na) * a] = s.val();
    }
  }
  return CLAB_OK;
}

int clab_agepan(const double *ya, const double *yb, int n, int nperm, double seed, double *out,
                double *work) {
  if (n < 3 || nperm < 0) return CLAB_EARG;
  if (!(seed >= 0)) return CLAB_EDOM;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(ya[i]) || !std::isfinite(yb[i])) return CLAB_EDOM;
  }
  auto cv = [&](const double *b) {
    kutil::ksum sab, sa, sb;
    for (int i = 0; i < n; ++i) {
      sab.add(ya[i] * b[i]);
      sa.add(ya[i]);
      sb.add(b[i]);
    }
    return sab.val() / n - (sa.val() / n) * (sb.val() / n);
  };
  const double obs = cv(yb);
  for (int i = 0; i < n; ++i) work[i] = yb[i];
  kutil::rng r(static_cast<std::uint64_t>(seed));
  int ge = 0;
  for (int p = 0; p < nperm; ++p) {
    // Fisher-Yates
    for (int i = n - 1; i > 0; --i) {
      const int j = r.below(i + 1);
      const double t = work[i];
      work[i] = work[j];
      work[j] = t;
    }
    ge += std::fabs(cv(work)) >= std::fabs(obs) - 1e-12;
  }
  out[0] = obs;
  out[1] = (1.0 + ge) / (nperm + 1.0);
  return CLAB_OK;
}

int clab_cumrsk(const double *age, const double *e, const double *pop, const double *atrisk,
                const double *dth, int n, double *h, double *nr, double *cr, double *v) {
  if (n < 1) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(age[i]) || !kutil::nonneg(e[i]) || !kutil::pos(pop[i])) return CLAB_EDOM;
    if (i > 0 && !(age[i] > age[i - 1])) return CLAB_EDOM;
    if (atrisk && !kutil::pos(atrisk[i])) return CLAB_EDOM;
    if (dth && !(dth[i] >= 0 && dth[i] < 1)) return CLAB_EDOM;
  }
  double surv = 1, l = 1, acc = 0, vs = 0;
  for (int i = 0; i < n; ++i) {
    // those at risk: given, or the population less those already with a
    // first event
    nr[i] = atrisk ? atrisk[i] : pop[i] * (1 - (i > 0 ? cr[i - 1] : 0));
    h[i] = e[i] / nr[i];
    if (h[i] > 1) return CLAB_EDOM;
    if (!dth) {
      surv *= 1 - h[i];
      cr[i] = 1 - surv;
      // Greenwood's variance of log survival
      vs += h[i] < 1 ? h[i] / (nr[i] * (1 - h[i])) : INFINITY;
      v[i] = vs;
    } else {
      // with deaths competing: the first events of the cohort still alive
      // and without one
      acc += l * h[i];
      l *= 1 - h[i] - dth[i];
      if (l < 0) return CLAB_EDOM;
      cr[i] = acc;
      v[i] = kutil::nan();
    }
  }
  return CLAB_OK;
}

int clab_jpipe(const double *num, const double *den, int n, int win, int lag, double *onum,
               double *oden, double *rate) {
  if (n < 0 || win < 1 || lag < 0) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (!kutil::nonneg(num[i]) || !kutil::nonneg(den[i])) return CLAB_EDOM;
  }
  // window sums ending at e (numerator) and at e - lag (denominator), each
  // kept by adding the new term and taking off the one leaving
  kutil::ksum sn, sd;
  for (int e = 0; e < n; ++e) {
    sn.add(num[e]);
    if (e >= win) sn.add(-num[e - win]);
    onum[e] = e >= win - 1 ? sn.val() : kutil::nan();
    const int d = e - lag;
    if (d >= 0) {
      sd.add(den[d]);
      if (d >= win) sd.add(-den[d - win]);
    }
    oden[e] = d >= win - 1 ? sd.val() : kutil::nan();
    rate[e] = oden[e] > 0 ? onum[e] / oden[e] : kutil::nan();
  }
  return CLAB_OK;
}
