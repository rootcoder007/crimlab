// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, policing: Andresen's spatial point pattern test with the
// multinomial bootstrap streamed over replicates: memory linear in the
// number of units, not units times replicates.
#include <cstdint>

#include "crimlab.h"
#include "kutil.h"

namespace {


// Walker's alias table (Vose's construction) for the probabilities
// c[i] / tot over the m units; pr, al: m each; st: m ints of scratch
void alias(const double *c, int m, double tot, double *pr, int *al, int *st) {
  int ns = 0, nl = m;  // small ones from the front of st, large from the back
  for (int i = 0; i < m; ++i) {
    pr[i] = c[i] * m / tot;
    al[i] = i;
    if (pr[i] < 1.0) st[ns++] = i; else st[--nl] = i;
  }
  while (ns > 0 && nl < m) {
    const int s = st[--ns], l = st[nl];
    al[s] = l;
    pr[l] -= 1.0 - pr[s];
    if (pr[l] < 1.0) {
      ++nl;
      st[ns++] = l;
    }
  }
  // what is left is 1 up to rounding
  for (int i = 0; i < ns; ++i) pr[st[i]] = 1.0;
  for (int i = nl; i < m; ++i) pr[st[i]] = 1.0;
}

inline int draw(kutil::rng &r, const double *pr, const int *al, int m) {
  const double u = r.unif() * m;
  int i = static_cast<int>(u);
  if (i >= m) i = m - 1;
  return (u - i) < pr[i] ? i : al[i];
}

}  // namespace

extern "C" {

int clab_sptst(const double *a, const double *b, int n, int nboot,
               unsigned int seed, double alpha, double scale, double *out,
               double *glob, double *work, int *iwork) {
  if (n < 1 || nboot < 1) return CLAB_EARG;
  if (!(alpha > 0.0 && alpha < 1.0) || !kutil::pos(scale)) return CLAB_EDOM;
  // units with a missing count (NaN) are left out entirely and counted;
  // others must be whole numbers >= 0
  int nna = 0;
  double na = 0.0, nb = 0.0;
  for (int i = 0; i < n; ++i) {
    if (std::isnan(a[i]) || std::isnan(b[i])) {
      ++nna;
      continue;
    }
    if (!(a[i] >= 0.0) || !(b[i] >= 0.0) || a[i] != std::floor(a[i]) || b[i] != std::floor(b[i]) ||
        a[i] > 2e9 || b[i] > 2e9)
      return CLAB_EDOM;
    na += a[i];
    nb += b[i];
  }
  double *diff = out, *pv = out + n, *sig = out + 2 * n;
  double *ca = work, *cb = ca + n, *mean = cb + n, *ge = mean + n;
  double *pra = ge + n, *prb = pra + n;
  int *ala = iwork, *alb = ala + n, *st = alb + n, *ix = st + n;
  // the units in use, compacted: ix[k] = unit
  int m = 0;
  for (int i = 0; i < n; ++i) {
    diff[i] = pv[i] = sig[i] = kutil::nan();
    if (!(std::isnan(a[i]) || std::isnan(b[i]))) ix[m++] = i;
  }
  int nzero = 0;
  for (int k = 0; k < m; ++k) {
    const int i = ix[k];
    diff[i] = (nb > 0.0 ? b[i] / nb * scale : 0.0) - (na > 0.0 ? a[i] / na * scale : 0.0);
    nzero += a[i] == 0.0 && b[i] == 0.0;
  }
  const int nelig = m - nzero;
  glob[0] = nna;
  glob[1] = nzero;
  glob[2] = nelig;
  for (int g = 3; g < 8; ++g) glob[g] = kutil::nan();
  if (na == 0.0 || nb == 0.0 || nelig == 0) return CLAB_OK;  // no inference
  // alias tables over the units in use
  for (int k = 0; k < m; ++k) ca[k] = a[ix[k]];
  alias(ca, m, na, pra, ala, st);
  for (int k = 0; k < m; ++k) cb[k] = b[ix[k]];
  alias(cb, m, nb, prb, alb, st);
  const long long da = static_cast<long long>(na), db = static_cast<long long>(nb);
  // two passes over the same replicates (the generator is restarted from
  // the seed): the first for the replicate mean of each difference, the
  // second to count centred replicates at least as far from 0 as observed
  for (int pass = 0; pass < 2; ++pass) {
    kutil::rng r(seed);
    for (int k = 0; k < m; ++k) {
      if (pass == 0) mean[k] = 0.0;
      ge[k] = 0.0;
    }
    for (int rep = 0; rep < nboot; ++rep) {
      if (kutil::stop(static_cast<double>(da + db) + m)) return CLAB_EINTR;
      for (int k = 0; k < m; ++k) ca[k] = cb[k] = 0.0;
      for (long long t = 0; t < da; ++t) ca[draw(r, pra, ala, m)] += 1.0;
      for (long long t = 0; t < db; ++t) cb[draw(r, prb, alb, m)] += 1.0;
      for (int k = 0; k < m; ++k) {
        const double d = cb[k] / nb * scale - ca[k] / na * scale;
        if (pass == 0) {
          mean[k] += d;
        } else {
          const int i = ix[k];
          if (std::fabs(d - mean[k]) >= std::fabs(diff[i])) ge[k] += 1.0;
        }
      }
    }
    if (pass == 0)
      for (int k = 0; k < m; ++k) mean[k] /= nboot;
  }
  int nsig = 0;
  for (int k = 0; k < m; ++k) {
    const int i = ix[k];
    if (a[i] == 0.0 && b[i] == 0.0) continue;
    pv[i] = (1.0 + ge[k]) / (nboot + 1.0);
    sig[i] = pv[i] <= alpha ? (diff[i] > 0.0 ? 1.0 : (diff[i] < 0.0 ? -1.0 : 0.0)) : 0.0;
    nsig += pv[i] <= alpha;
  }
  // the S index over the eligible units, and the sameness test: under
  // sameness the flags are Binomial(nelig, alpha); its exact upper tail
  glob[3] = nsig;
  glob[4] = 1.0 - static_cast<double>(nsig) / nelig;
  const double se = std::sqrt(alpha * (1.0 - alpha) / nelig);
  glob[5] = (glob[4] - (1.0 - alpha)) / se;
  const double mu = nelig * alpha;
  glob[6] = nsig == 0 ? 1.0
            : (nsig >= mu ? kutil::btail(nsig, nelig, alpha, true)
                          : 1.0 - kutil::btail(nsig - 1.0, nelig, alpha, false));
  glob[7] = nboot;
  return CLAB_OK;
}

}  // extern "C"
