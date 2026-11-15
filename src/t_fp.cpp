// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/fp.cpp: the table lookups, and the energy of runs worked by
// hand from the Green Algorithms formula and codecarbon's fallbacks.
#include <cmath>
#include <cstring>

#include "crimlab.h"
#include "tchk.h"

void t_fp(tk &T, int) {
  tk_ieq(T, "fp.rows", clab_cin(), 378);
  const char *loc, *iso3, *name, *src;
  double g;
  int yr;
  const int i = clab_cimatch("  ca-on ");
  tk_true(T, "fp.match", i >= 0);
  if (i >= 0) {
    clab_cirow(i, &loc, &iso3, &name, &g, &yr, &src);
    tk_true(T, "fp.loc", std::strcmp(loc, "CA-ON") == 0 && yr == 2024);
    tk_eq(T, "fp.g", g, 90.97, 0);
  }
  tk_ieq(T, "fp.iso3", clab_cimatch("CAN"), clab_cimatch("Canada"));
  tk_ieq(T, "fp.none", clab_cimatch("Atlantis"), -1);
  const char *t = clab_tzterr("posix/Europe/Berlin");
  tk_true(T, "fp.tz", t != nullptr && std::strcmp(t, "DE") == 0);
  tk_true(T, "fp.tz0", clab_tzterr("Etc/UTC") == nullptr);
  tk_true(T, "fp.proxy", std::strcmp(clab_terrproxy("MC"), "FR") == 0 && clab_terrproxy("AQ")[0] == 0);
  // one hour, half of 4 cores busy, 100 W, 10 GB at 0.3725 W/GB, PUE 1.5:
  // (100 x 0.5 + 3.725) W for 1 h is 0.053725 kWh, times 1.5
  double in[16] = {3600, 7200, 4, 0, 0, NAN, NAN, NAN, NAN, 100, 0, 10, 0.3725, 1.5, 100, 0};
  double out[10], dout[2];
  clab_fpe(in, 0, nullptr, nullptr, nullptr, 0, nullptr, nullptr, nullptr, out, dout);
  tk_eq(T, "fp.u", out[1], 0.5, 0);
  tk_eq(T, "fp.tot", out[8], 0.053725 * 1.5, 1e-15);
  tk_eq(T, "fp.co2", out[9], out[8] * 100, 0);
  // the unknown CPU: 42.5 W whatever the utilisation; a device at 200 W
  // and one with 0.01 kWh measured; a CPU row is not counted again
  in[9] = NAN;
  const double dw[3] = {200, NAN, 50}, de[3] = {NAN, 0.01, NAN};
  const int dk[3] = {0, 0, 1};
  double d3[3];
  clab_fpe(in, 0, nullptr, nullptr, nullptr, 3, dw, de, dk, out, d3);
  tk_true(T, "fp.flag", (static_cast<int>(out[4]) & 8) != 0);
  tk_eq(T, "fp.dev", d3[0], 0.2, 1e-15);
  tk_true(T, "fp.devcpu", std::isnan(d3[2]));
  tk_eq(T, "fp.dtot", out[8], (0.0425 + 0.003725 + 0.21) * 1.5, 1e-15);
  // RAPL: a counter that wrapped through its range of 1000 uJ
  in[15] = 1;
  const double r0[1] = {900}, r1[1] = {100}, mx[1] = {1000};
  clab_fpe(in, 1, r0, r1, mx, 0, nullptr, nullptr, nullptr, out, dout);
  tk_eq(T, "fp.rapl", out[5], 200 * 1e-6 / 3.6e6 * 1.5, 1e-30);
  const double nx[1] = {NAN};
  clab_fpe(in, 1, r0, r1, nx, 0, nullptr, nullptr, nullptr, out, dout);
  tk_eq(T, "fp.wrap", out[3], 1, 0);
  // the machine mode: 30 of 100 ticks idle is 70% busy
  in[15] = 0;
  in[3] = 1;
  in[5] = 10;
  in[6] = 100;
  in[7] = 40;
  in[8] = 200;
  clab_fpe(in, 0, nullptr, nullptr, nullptr, 0, nullptr, nullptr, nullptr, out, dout);
  tk_eq(T, "fp.mach", out[1], 0.7, 1e-15);
  double e[2];
  clab_fpeqv(1750, e);
  tk_eq(T, "fp.km", e[0], 10, 0);
}
