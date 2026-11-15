// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of src/apol.cpp: the dispersion formulas and the health side by
// hand, and the refusals of the kernels.
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "tchk.h"

void t_apol(tk &T, int) {
  const double pi = 3.14159265358979323846;
  // class D rural at 1 km: 0.08 x / sqrt(1 + 0.1), 0.06 x / sqrt(1 + 1.5)
  const double x = 1000;
  double sy, sz;
  tk_ieq(T, "apol.sig", clab_pgsig(&x, 1, 3, 0, &sy, &sz), CLAB_OK);
  tk_eq(T, "apol.sy", sy, 80 / std::sqrt(1.1), 1e-14);
  tk_eq(T, "apol.sz", sz, 60 / std::sqrt(2.5), 1e-14);
  tk_ieq(T, "apol.sigc", clab_pgsig(&x, 1, 6, 0, &sy, &sz), CLAB_EARG);
  // the plume on the centreline at ground level from a ground release:
  // q / (pi u sy sz)
  const double r[3] = {1000, 0, 0};
  double c = 0;
  clab_gplume(1, 2, 0, r, 1, 3, 0, NAN, 3, &c);
  tk_eq(T, "apol.gpl", c, 1 / (pi * 2 * sy * sz), 1e-14);
  const double up[3] = {-5, 0, 0};
  clab_gplume(1, 2, 0, up, 1, 3, 0, NAN, 3, &c);
  tk_eq(T, "apol.up", c, 0, 0);
  tk_ieq(T, "apol.gplu", clab_gplume(1, 0, 0, r, 1, 3, 0, NAN, 3, &c), CLAB_EARG);
  tk_ieq(T, "apol.gspu", clab_gspuff(1, 0, 0, r, 1, 1, 3, 0, NAN, 3, &c), CLAB_EARG);
  // Briggs: buoyancy flux g v d^2 dT / (4 Ts)
  double o[3], rise;
  clab_briggs(&x, 1, 4, 1, 12, 420, 285, 3, NAN, 9.80616, o, &rise);
  tk_eq(T, "apol.flux", o[0], 9.80616 * 12 * 135 / (4 * 420), 1e-14);
  tk_ieq(T, "apol.brgu", clab_briggs(&x, 1, 0, 1, 12, 420, 285, 3, NAN, 9.8, o, &rise), CLAB_EARG);
  // a field without wind or diffusion keeps its interior and loses its ring
  std::vector<double> f(9, 1.0);
  clab_advdif(f.data(), 3, 3, 0, 0, 0, 0, 1, 1, 1, 1, nullptr);
  tk_eq(T, "apol.adv", f[4], 1, 0);
  tk_eq(T, "apol.ring", f[0], 0, 0);
  tk_ieq(T, "apol.advd", clab_advdif(f.data(), 3, 3, 0, 0, 0, 0, 0, 1, 1, 1, nullptr), CLAB_EARG);
  // particles without diffusion move with the wind
  double xs[2], ys[2], m[2];
  clab_lagr(2, 0, 0, 1, 2, 0, 0, 1, 3, 0, xs, ys);
  tk_eq(T, "apol.lag", xs[1], 3, 1e-15);
  clab_laggrid(xs, ys, 2, nullptr, nullptr, m);
  tk_eq(T, "apol.lagm", m[1], 6, 1e-15);
  tk_ieq(T, "apol.lagn", clab_lagr(0, 0, 0, 1, 2, 0, 0, 1, 3, 0, xs, ys), CLAB_EARG);
  tk_ieq(T, "apol.lgs", clab_lagstep(-1, 0, 0, 0, 0, 1, xs, ys, xs, ys), CLAB_EARG);
  tk_ieq(T, "apol.lgg", clab_laggrid(xs, ys, 0, nullptr, nullptr, m), CLAB_EARG);
  // the health side: no excess risk at the reference; Levin with all
  // exposed is 1 - 1/RR
  const double z = 5, b = std::log(1.08);
  double rr, lr;
  clab_crf(&z, 1, 5.8, 0, &b, &rr, &lr);
  tk_eq(T, "apol.crf", rr, 1, 0);
  tk_ieq(T, "apol.crff", clab_crf(&z, 1, 5.8, 3, &b, &rr, &lr), CLAB_EARG);
  tk_eq(T, "apol.levin", clab_levin(2, 1), 0.5, 1e-15);
  tk_eq(T, "apol.hif", clab_hif(0, 1e6, 0.008, 0.01), 0, 0);
  const double h[4] = {20, 20, 20, 20}, inc[4] = {1, 2, 3, 4};
  double ci[3], cp[4], ce[4];
  clab_concidx(h, inc, 4, ci, cp, ce);
  tk_eq(T, "apol.ci0", ci[0], 0, 1e-15);
  tk_ieq(T, "apol.cin", clab_concidx(h, inc, 1, ci, cp, ce), CLAB_EARG);
}
