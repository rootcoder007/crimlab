// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for clab_crsev and clab_zcomp (src/sev.cpp) against
// sums done by hand. Part of the self-test (selftest.cpp).
#include <cmath>
#include <limits>

#include "crimlab.h"
#include "tchk.h"

namespace {

void crsev(tk &t) {
  // two cells, two offences of weights 1 and 50, populations 1000
  const double x[4] = {10, 2, 8, 3}, w[2] = {1, 50}, pop[2] = {1000, 1000};
  const int off[4] = {0, 1, 0, 1}, key[4] = {0, 0, 1, 1};
  double cnt[2], sev[2], rate[2], sr[2], idx[2];
  tk_ieq(t, "crsev.rc", clab_crsev(x, off, key, 4, 2, 2, w, pop, 1e5, 0, cnt, sev, rate, sr, idx),
         CLAB_OK);
  tk_eq(t, "crsev.cnt", cnt[1], 11, 0);
  tk_eq(t, "crsev.sev", sev[0], 110, 0);
  tk_eq(t, "crsev.sev2", sev[1], 158, 0);
  tk_eq(t, "crsev.rate", rate[0], 1200, 1e-15);
  tk_eq(t, "crsev.srate", sr[1], 15800, 1e-15);
  tk_eq(t, "crsev.idx0", idx[0], 100, 0);
  tk_eq(t, "crsev.idx1", idx[1], 15800.0 / 110, 1e-15);
  // no population: the index is on the weighted sums; no base, no index
  clab_crsev(x, off, key, 4, 2, 2, w, nullptr, 1e5, 1, cnt, sev, rate, sr, idx);
  tk_nan(t, "crsev.norate", rate[0]);
  tk_eq(t, "crsev.idxsum", idx[0], 11000.0 / 158, 1e-15);
  clab_crsev(x, off, key, 4, 2, 2, w, nullptr, 1e5, -1, cnt, sev, rate, sr, idx);
  tk_nan(t, "crsev.nobase", idx[1]);
  const double x0[4] = {0, 0, 8, 3};
  clab_crsev(x0, off, key, 4, 2, 2, w, nullptr, 1e5, 0, cnt, sev, rate, sr, idx);
  tk_nan(t, "crsev.zerobase", idx[1]);
  // refusals
  const double wn[2] = {-1, 1}, xn[4] = {-1, 0, 0, 0}, p0[2] = {0, 1};
  const int ob[4] = {0, 2, 0, 0}, kb[4] = {0, 0, 2, 0};
  tk_ieq(t, "crsev.n", clab_crsev(x, off, key, 0, 2, 2, w, pop, 1, 0, cnt, sev, rate, sr, idx),
         CLAB_EARG);
  tk_ieq(t, "crsev.base", clab_crsev(x, off, key, 4, 2, 2, w, pop, 1, 2, cnt, sev, rate, sr, idx),
         CLAB_EARG);
  tk_ieq(t, "crsev.per", clab_crsev(x, off, key, 4, 2, 2, w, pop, 0, 0, cnt, sev, rate, sr, idx),
         CLAB_EDOM);
  tk_ieq(t, "crsev.w", clab_crsev(x, off, key, 4, 2, 2, wn, pop, 1, 0, cnt, sev, rate, sr, idx),
         CLAB_EDOM);
  tk_ieq(t, "crsev.pop", clab_crsev(x, off, key, 4, 2, 2, w, p0, 1, 0, cnt, sev, rate, sr, idx),
         CLAB_EDOM);
  tk_ieq(t, "crsev.off", clab_crsev(x, ob, key, 4, 2, 2, w, pop, 1, 0, cnt, sev, rate, sr, idx),
         CLAB_EARG);
  tk_ieq(t, "crsev.key", clab_crsev(x, off, kb, 4, 2, 2, w, pop, 1, 0, cnt, sev, rate, sr, idx),
         CLAB_EARG);
  tk_ieq(t, "crsev.x", clab_crsev(xn, off, key, 4, 2, 2, w, pop, 1, 0, cnt, sev, rate, sr, idx),
         CLAB_EDOM);
}

void zcomp(tk &t) {
  // category 1: 1, 2, 3 (mean 2, sd 1); category 2 does not vary
  const double x[6] = {1, 2, 3, 5, 5, 5}, w[2] = {2, 1};
  double z[6], comp[3], mu[2], sd[2];
  tk_ieq(t, "zcomp.rc", clab_zcomp(x, 3, 2, w, z, comp, mu, sd), CLAB_OK);
  tk_eq(t, "zcomp.mu", mu[1], 5, 0);
  tk_eq(t, "zcomp.sd", sd[0], 1, 1e-15);
  tk_eq(t, "zcomp.sd1", sd[1], 1, 0);
  tk_eq(t, "zcomp.z", z[2], 1, 1e-15);
  tk_eq(t, "zcomp.z0", z[4], 0, 0);
  tk_eq(t, "zcomp.c0", comp[0], -2, 1e-15);
  tk_eq(t, "zcomp.c2", comp[2], 2, 1e-15);
  const double wn[2] = {std::numeric_limits<double>::quiet_NaN(), 1};
  const double xi[6] = {1, 2, std::numeric_limits<double>::infinity(), 5, 5, 5};
  tk_ieq(t, "zcomp.nu", clab_zcomp(x, 1, 2, w, z, comp, mu, sd), CLAB_EARG);
  tk_ieq(t, "zcomp.w", clab_zcomp(x, 3, 2, wn, z, comp, mu, sd), CLAB_EDOM);
  tk_ieq(t, "zcomp.x", clab_zcomp(xi, 3, 2, w, z, comp, mu, sd), CLAB_EDOM);
}

}  // namespace

void t_sev(tk &t, int big) {
  (void)big;
  crsev(t);
  zcomp(t);
}
