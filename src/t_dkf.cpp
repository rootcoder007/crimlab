// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for the dark-figure kernels (src/dkf.cpp) against
// values worked by hand. Part of the self-test (selftest.cpp).
#include <cmath>
#include <limits>

#include "crimlab.h"
#include "tchk.h"

namespace {

const double nan = std::numeric_limits<double>::quiet_NaN();

void cap(tk &t) {
  const double x[3] = {500, 400, 80};
  double st[6], pw[18], sens[4];
  tk_ieq(t, "crcap.rc", clab_crcap(x, 3, 1.5, nullptr, 0, st, pw, sens), CLAB_OK);
  tk_eq(t, "crcap.pet", st[0], 2500, 0);
  tk_eq(t, "crcap.lo", st[1], 2500 / 1.5, 1e-15);
  tk_eq(t, "crcap.hi", st[2], 3750, 1e-15);
  tk_eq(t, "crcap.chap", st[3], 501.0 * 401 / 81 - 1, 1e-15);
  tk_eq(t, "crcap.se", st[4], std::sqrt(501.0 * 401 * 420 * 320 / (81.0 * 81 * 82)), 1e-15);
  tk_eq(t, "crcap.floor", st[5], 820, 0);
  // three lists: the missing cell 120 90 60 8 / (30 20 15) = 576
  const double c[7] = {120, 90, 60, 30, 20, 15, 8}, it[2] = {0, std::log(2.0)};
  tk_ieq(t, "crcap.rc3", clab_crcap(c, 7, 1, it, 2, st, pw, sens), CLAB_OK);
  tk_eq(t, "crcap.obs", st[0], 343, 0);
  tk_eq(t, "crcap.miss", st[1], 576, 1e-15);
  tk_eq(t, "crcap.tot", st[2], 919, 1e-15);
  tk_eq(t, "crcap.n1", pw[0], 178, 0);
  tk_eq(t, "crcap.n2", pw[3], 143, 0);
  tk_eq(t, "crcap.m", pw[6], 38, 0);
  tk_eq(t, "crcap.pair", pw[12], 178.0 * 143 / 38, 1e-15);
  // the pair 2-3: list 3 holds 001, 101, 011, 111
  tk_eq(t, "crcap.n3", pw[5], 60 + 20 + 15 + 8, 0);
  tk_eq(t, "crcap.sens", sens[1], 288, 1e-15);
  tk_eq(t, "crcap.senstot", sens[3], 343 + 288, 1e-15);
  const double neg[3] = {-1, 4, 1}, m0[3] = {5, 5, 0}, mbig[3] = {5, 5, 6};
  const double c0[7] = {1, 1, 1, 0, 1, 1, 1}, itn[1] = {nan};
  tk_ieq(t, "crcap.k", clab_crcap(x, 5, 1, nullptr, 0, st, pw, sens), CLAB_EARG);
  tk_ieq(t, "crcap.kappa", clab_crcap(x, 3, 0.5, nullptr, 0, st, pw, sens), CLAB_EDOM);
  tk_ieq(t, "crcap.neg", clab_crcap(neg, 3, 1, nullptr, 0, st, pw, sens), CLAB_EDOM);
  tk_ieq(t, "crcap.int", clab_crcap(c, 7, 1, itn, 1, st, pw, sens), CLAB_EDOM);
  tk_ieq(t, "crcap.m0", clab_crcap(m0, 3, 1, nullptr, 0, st, pw, sens), CLAB_EDOM);
  tk_ieq(t, "crcap.mbig", clab_crcap(mbig, 3, 1, nullptr, 0, st, pw, sens), CLAB_EDOM);
  tk_ieq(t, "crcap.pair0", clab_crcap(c0, 7, 1, nullptr, 0, st, pw, sens), CLAB_EDOM);
}

void fig(tk &t) {
  const double v[2] = {0.12, 0.02}, r[2] = {0.03, 0.05};
  double lo[2], hi[2], bk[2];
  tk_ieq(t, "dkfig.rc", clab_dkfig(v, r, 2, 0.02, 0.25, 0.15, lo, hi, bk), CLAB_OK);
  tk_eq(t, "dkfig.lo", lo[0], 0.1 / 0.98, 1e-15);
  tk_eq(t, "dkfig.hi", hi[0], 0.16, 1e-15);
  tk_eq(t, "dkfig.brk", bk[0], 0.2, 1e-14);
  // the recorded rate binds both ends
  tk_eq(t, "dkfig.clip", lo[1], 0.05, 0);
  tk_eq(t, "dkfig.clip2", hi[1], 0.05, 0);
  clab_dkfig(v, r, 2, 0, 0, nan, lo, hi, bk);
  tk_nan(t, "dkfig.nobrk", bk[0]);
  const double v2[1] = {1.5}, v0[1] = {0};
  tk_ieq(t, "dkfig.n", clab_dkfig(v, r, -1, 0, 0, nan, lo, hi, bk), CLAB_EARG);
  tk_ieq(t, "dkfig.amax", clab_dkfig(v, r, 2, -1, 0, nan, lo, hi, bk), CLAB_EDOM);
  tk_ieq(t, "dkfig.sum", clab_dkfig(v, r, 2, 0.5, 0.5, nan, lo, hi, bk), CLAB_EDOM);
  tk_ieq(t, "dkfig.v", clab_dkfig(v2, r, 1, 0, 0, nan, lo, hi, bk), CLAB_EDOM);
  tk_ieq(t, "dkfig.thr", clab_dkfig(v, r, 2, 0, 0, 0.1, lo, hi, bk), CLAB_EDOM);
  tk_ieq(t, "dkfig.thr0", clab_dkfig(v0, r, 1, 0, 0, 0.1, lo, hi, bk), CLAB_EDOM);
  // a survey of type A only: k = 0.5 * 300 / 900 = 1/6
  double o[7];
  tk_ieq(t, "dkcmp.rc", clab_dkcmp(900, 300, 450, 0.5, o), CLAB_OK);
  const double ow[7] = {1.0 / 3, 750.0 / 1350, 1350, 3600, 600, 2850, 2250};
  tk_veq(t, "dkcmp.out", o, ow, 7, 1e-14);
  tk_ieq(t, "dkcmp.va", clab_dkcmp(0, 0, 1, 1, o), CLAB_EDOM);
  tk_ieq(t, "dkcmp.ra", clab_dkcmp(10, 20, 1, 1, o), CLAB_EDOM);
  tk_ieq(t, "dkcmp.lam", clab_dkcmp(10, 2, 1, 0, o), CLAB_EDOM);
  tk_ieq(t, "dkcmp.lam1", clab_dkcmp(10, 2, 1, 1.5, o), CLAB_EDOM);
  tk_ieq(t, "dkcmp.rb", clab_dkcmp(10, 2, -1, 1, o), CLAB_EDOM);
  const double k[5] = {1, 1, 2, 1, 3}, kf[1] = {1.5}, k0[1] = {0},
               ki[1] = {std::numeric_limits<double>::infinity()};
  tk_ieq(t, "dkhier.rc", clab_dkhier(k, 5, o), CLAB_OK);
  const double hw[6] = {5, 8, 1.6, 0.6, 5, 15};
  tk_veq(t, "dkhier.out", o, hw, 6, 1e-15);
  tk_ieq(t, "dkhier.n", clab_dkhier(k, 0, o), CLAB_EARG);
  tk_ieq(t, "dkhier.frac", clab_dkhier(kf, 1, o), CLAB_EDOM);
  tk_ieq(t, "dkhier.zero", clab_dkhier(k0, 1, o), CLAB_EDOM);
  tk_ieq(t, "dkhier.inf", clab_dkhier(ki, 1, o), CLAB_EDOM);
}

}  // namespace

void t_dkf(tk &t, int big) {
  (void)big;
  cap(t);
  fig(t);
}
