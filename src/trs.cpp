// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Transporting a trial's average effect (crimlab.h, trnsp()): from the
// trial's rows (s = 1) to a target, the rows outside it (transport) or
// all rows (generalise), by two estimators (Dahabreh et al. 2020): the
// outcome model (least squares of y on [1, x] in each arm of the trial,
// xls.h, averaged over the target of the difference of the two
// predictions) and inverse odds weighting (the difference of the arms'
// means weighted by the sampling weights of clab_gensw(), normalised
// within each arm; Westreich et al. 2017). The standard errors are those
// of the bootstrap of the rows within the trial and within the rest,
// drawn from the package's generator.
#include <cmath>
#include <cstddef>
#include <vector>

#include "crimlab.h"
#include "kutil.h"
#include "xls.h"

namespace {

using xls::z;

// Both estimates on one data set; false when a fit fails or an arm is
// empty. y and d are not read outside the trial.
bool est1(const double *x, int n, int p, const double *y, const int *d, const int *s, int target, double trim,
          double *om, double *iow, double *ps, double *w) {
  std::vector<double> yz(z(n), 0.0), wa(z(n)), r(z(n)), b(z(p) + 1), b1(z(p) + 1), cf(z(p) + 1), st(6),
      bal(z(p) * 5);
  for (int i = 0; i < n; ++i) yz[z(i)] = s[i] ? y[i] : 0.0;
  for (int a = 0; a < 2; ++a) {
    for (int i = 0; i < n; ++i) wa[z(i)] = s[i] && d[i] == a ? 1.0 : 0.0;
    xls::Proj pr;
    if (!pr.make(x, n, nullptr, p, wa.data())) return false;
    pr.res(yz.data(), r.data(), a ? b1.data() : b.data());
  }
  double sm = 0, nt = 0;
  for (int i = 0; i < n; ++i) {
    if (target == 1 && s[i]) continue;
    double f = b1[0] - b[0];
    for (int j = 0; j < p; ++j) f += (b1[z(j) + 1] - b[z(j) + 1]) * x[z(j) * z(n) + z(i)];
    sm += f;
    nt += 1;
  }
  *om = sm / nt;
  if (clab_gensw(x, n, p, s, target, trim, 1, ps, w, cf.data(), st.data(), bal.data()) != CLAB_OK) return false;
  double sw[2] = {0, 0}, swy[2] = {0, 0};
  for (int i = 0; i < n; ++i) {
    if (!s[i]) continue;
    sw[d[i]] += w[i];
    swy[d[i]] += w[i] * y[i];
  }
  if (!(sw[0] > 0 && sw[1] > 0)) return false;
  *iow = swy[1] / sw[1] - swy[0] / sw[0];
  return std::isfinite(*om) && std::isfinite(*iow);
}

}  // namespace

extern "C" int clab_trnsp(const double *x, int n, int p, const double *y, const int *d, const int *s, int target,
                          double trim, int nboot, unsigned int seed, double *est, double *se, double *draws,
                          double *ps, double *w, double *nf) CLAB_TRY {
  if (n < 2 || p < 0 || target < 0 || target > 1 || nboot < 0 || !(trim >= 0 && trim < 0.5)) return CLAB_EARG;
  std::vector<int> in, out;
  for (int i = 0; i < n; ++i) {
    if (s[i] != 0 && s[i] != 1) return CLAB_EARG;
    if (s[i] && d[i] != 0 && d[i] != 1) return CLAB_EARG;
    (s[i] ? in : out).push_back(i);
  }
  if (in.empty() || out.empty()) return CLAB_EDOM;
  if (!est1(x, n, p, y, d, s, target, trim, &est[0], &est[1], ps, w)) return CLAB_EDOM;
  // the bootstrap: rows drawn with replacement within the trial and
  // within the rest, the data set rebuilt
  kutil::rng r(seed);
  std::vector<double> xb(z(n) * z(p)), yb(z(n)), pb(z(n)), wb(z(n));
  std::vector<int> db(z(n)), sb(z(n));
  double ok = 0, fail = 0, m[2] = {0, 0}, q[2] = {0, 0};
  for (int k = 0; k < nboot; ++k) {
    if (kutil::stop(static_cast<double>(n) * (p + 1) * (p + 4))) return CLAB_EINTR;
    for (int i = 0; i < n; ++i) {
      const std::vector<int> &g = i < static_cast<int>(in.size()) ? in : out;
      const int src = g[z(r.below(static_cast<int>(g.size())))];
      for (int j = 0; j < p; ++j) xb[z(j) * z(n) + z(i)] = x[z(j) * z(n) + z(src)];
      yb[z(i)] = s[src] ? y[src] : 0.0;
      db[z(i)] = s[src] ? d[src] : 0;
      sb[z(i)] = s[src];
    }
    double e[2];
    if (!est1(xb.data(), n, p, yb.data(), db.data(), sb.data(), target, trim, &e[0], &e[1], pb.data(), wb.data())) {
      draws[k] = draws[nboot + k] = kutil::nan();
      fail += 1;
      continue;
    }
    ok += 1;
    for (int a = 0; a < 2; ++a) {
      draws[a * nboot + k] = e[a];
      // Welford's mean and sum of squares
      const double dl = e[a] - m[a];
      m[a] += dl / ok;
      q[a] += dl * (e[a] - m[a]);
    }
  }
  for (int a = 0; a < 2; ++a) se[a] = ok > 1 ? std::sqrt(q[a] / (ok - 1)) : kutil::nan();
  *nf = fail;
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)
