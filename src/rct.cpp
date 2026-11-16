// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Targeted maximum likelihood for a randomised trial with covariate
// adjustment (crimlab.h, rcttml()): Moore and van der Laan (2009). The
// outcome regression Q(d, x) (least squares, xls.h, or logistic,
// clab_lgtreg) on [1, d, x] (and d x with inter) is fluctuated along the
// clever covariate H = d / g - (1 - d) / (1 - g), g the chance of
// treatment (the trial's share, a known value, or a logistic regression
// on x, bounded to [gbound, 1 - gbound]), by a logistic regression with
// offset logit(Q) (Newton's method, halving a step that lowers the
// likelihood); a continuous outcome is first put on [0, 1] by its bounds
// and Q bounded to [qbound, 1 - qbound] (Gruber and van der Laan 2010).
// The average effect is the mean of Q*(1, x) - Q*(0, x), its standard
// error that of its influence curve. The unfluctuated plug-in (the
// g-computation) and the difference in means come with it.
#include <cmath>
#include <cstddef>
#include <vector>

#include "crimlab.h"
#include "kutil.h"
#include "xls.h"

namespace {

using xls::z;

inline double expit(double v) { return 1 / (1 + std::exp(-v)); }
inline double logit(double p) { return std::log(p / (1 - p)); }

// Logistic coefficients of y on [1, x] (q columns): clab_lgtreg.
int lgt(const double *x, int n, int q, const double *y, std::vector<double> &cf) {
  cf.assign(z(q) + 1, 0.0);
  std::vector<double> se(z(q) + 1), wk(2 * (z(q) + 1) * (z(q) + 1) + 4 * (z(q) + 1));
  std::vector<int> iw(z(n));
  double dev;
  return clab_lgtreg(x, n, q, y, 1e-10, 50, cf.data(), se.data(), &dev, wk.data(), iw.data());
}

}  // namespace

extern "C" int clab_rcttml(const double *x, int n, int p, const double *y, const int *d, int fam, int inter, int gx,
                           double g0, double gbound, double lo, double hi, double qbound, double *out,
                           double *aux) CLAB_TRY {
  if (n < 2 || p < 0 || fam < 0 || fam > 1 || !(gbound >= 0 && gbound < 0.5) || !(qbound > 0 && qbound < 0.5) ||
      (gx && p < 1) || (fam == 0 && !(hi > lo))) {
    return CLAB_EARG;
  }
  double n1 = 0, s1 = 0, s0 = 0;
  for (int i = 0; i < n; ++i) {
    if (d[i] != 0 && d[i] != 1) return CLAB_EARG;
    if (fam == 1 && y[i] != 0 && y[i] != 1) return CLAB_EARG;
    n1 += d[i];
    (d[i] ? s1 : s0) += y[i];
  }
  const double n0 = n - n1;
  if (n1 < 2 || n0 < 2) return CLAB_EDOM;
  // the difference in means and its unpooled standard error
  {
    const double m1 = s1 / n1, m0 = s0 / n0;
    double v1 = 0, v0 = 0;
    for (int i = 0; i < n; ++i) {
      const double e = y[i] - (d[i] ? m1 : m0);
      (d[i] ? v1 : v0) += e * e;
    }
    out[2] = m1 - m0;
    out[5] = std::sqrt(v1 / (n1 - 1) / n1 + v0 / (n0 - 1) / n0);
  }
  // the design of Q: d, x, then d x
  const int q = 1 + p + (inter ? p : 0);
  std::vector<double> zm(z(n) * z(q)), dv(z(n)), cf(z(q) + 1);
  for (int i = 0; i < n; ++i) {
    dv[z(i)] = d[i];
    zm[z(i)] = d[i];
    for (int j = 0; j < p; ++j) {
      const double v = x[z(j) * z(n) + z(i)];
      zm[(z(j) + 1) * z(n) + z(i)] = v;
      if (inter) zm[(z(p) + z(j) + 1) * z(n) + z(i)] = d[i] * v;
    }
  }
  if (fam == 0) {
    xls::Proj pr;
    if (!pr.make(zm.data(), n, nullptr, q, nullptr)) return CLAB_EDOM;
    std::vector<double> r(z(n));
    pr.res(y, r.data(), cf.data());
  } else {
    const int rc = lgt(zm.data(), n, q, y, cf);
    if (rc != CLAB_OK) return rc;
  }
  // Q(1, x) and Q(0, x) on the scale of the fit
  std::vector<double> q1(z(n)), q0(z(n)), g(z(n));
  for (int i = 0; i < n; ++i) {
    double e0 = cf[0], e1 = cf[0] + cf[1];
    for (int j = 0; j < p; ++j) {
      const double v = x[z(j) * z(n) + z(i)];
      e0 += cf[z(j) + 2] * v;
      e1 += cf[z(j) + 2] * v + (inter ? cf[z(p) + z(j) + 2] * v : 0.0);
    }
    q1[z(i)] = fam == 0 ? e1 : expit(e1);
    q0[z(i)] = fam == 0 ? e0 : expit(e0);
  }
  double gc = 0;
  for (int i = 0; i < n; ++i) gc += q1[z(i)] - q0[z(i)];
  out[1] = gc / n;
  out[4] = kutil::nan();
  // g: a logistic regression on x, a known value or the share treated
  if (gx) {
    std::vector<double> gf;
    const int rc = lgt(x, n, p, dv.data(), gf);
    if (rc != CLAB_OK) return rc;
    for (int i = 0; i < n; ++i) {
      double e = gf[0];
      for (int j = 0; j < p; ++j) e += gf[z(j) + 1] * x[z(j) * z(n) + z(i)];
      g[z(i)] = expit(e);
    }
  } else {
    const double gv = std::isnan(g0) ? n1 / n : g0;
    for (int i = 0; i < n; ++i) g[z(i)] = gv;
  }
  double gmin = 1, gmax = 0;
  for (double &v : g) {
    v = std::fmin(1 - gbound, std::fmax(gbound, v));
    gmin = std::fmin(gmin, v);
    gmax = std::fmax(gmax, v);
  }
  // the bounded scale: y and Q on [0, 1], Q kept in [qbound, 1 - qbound]
  const double sc = fam == 0 ? hi - lo : 1, o0 = fam == 0 ? lo : 0;
  std::vector<double> ys(z(n)), off(z(n)), hh(z(n));
  for (int i = 0; i < n; ++i) {
    ys[z(i)] = (y[i] - o0) / sc;
    q1[z(i)] = std::fmin(1 - qbound, std::fmax(qbound, (q1[z(i)] - o0) / sc));
    q0[z(i)] = std::fmin(1 - qbound, std::fmax(qbound, (q0[z(i)] - o0) / sc));
    off[z(i)] = logit(d[i] ? q1[z(i)] : q0[z(i)]);
    hh[z(i)] = d[i] ? 1 / g[z(i)] : -1 / (1 - g[z(i)]);
  }
  // the fluctuation: Newton's method on eps, halving a step that lowers
  // the log-likelihood
  auto ll = [&](double e) {
    double s = 0;
    for (int i = 0; i < n; ++i) {
      const double pr = expit(off[z(i)] + e * hh[z(i)]);
      s += ys[z(i)] * std::log(pr) + (1 - ys[z(i)]) * std::log1p(-pr);
    }
    return s;
  };
  double eps = 0, l0 = ll(0);
  bool conv = false;
  for (int it = 0; it < 200 && !conv; ++it) {
    if (kutil::stop(n)) return CLAB_EINTR;
    double u = 0, info = 0;
    for (int i = 0; i < n; ++i) {
      const double pr = expit(off[z(i)] + eps * hh[z(i)]);
      u += hh[z(i)] * (ys[z(i)] - pr);
      info += hh[z(i)] * hh[z(i)] * pr * (1 - pr);
    }
    if (!(info > 0)) return CLAB_EDOM;
    double st = u / info, e1 = eps + st, l1 = ll(e1);
    for (int h = 0; h < 60 && !(l1 >= l0); ++h) {
      st /= 2;
      e1 = eps + st;
      l1 = ll(e1);
    }
    conv = std::fabs(st) <= 1e-12 * (1 + std::fabs(eps));
    eps = e1;
    l0 = l1;
  }
  if (!conv || !std::isfinite(eps)) return CLAB_EDOM;
  // the targeted Q, the effect and its influence curve
  std::vector<double> ic(z(n));
  double ps = 0;
  for (int i = 0; i < n; ++i) {
    q1[z(i)] = expit(logit(q1[z(i)]) + eps / g[z(i)]);
    q0[z(i)] = expit(logit(q0[z(i)]) - eps / (1 - g[z(i)]));
    ps += q1[z(i)] - q0[z(i)];
  }
  ps /= n;
  double im = 0;
  for (int i = 0; i < n; ++i) {
    const double qd = d[i] ? q1[z(i)] : q0[z(i)];
    ic[z(i)] = hh[z(i)] * (ys[z(i)] - qd) + q1[z(i)] - q0[z(i)] - ps;
    im += ic[z(i)];
  }
  im /= n;
  double v = 0;
  for (int i = 0; i < n; ++i) v += (ic[z(i)] - im) * (ic[z(i)] - im);
  out[0] = ps * sc;
  out[3] = std::sqrt(v / (n - 1) / n) * sc;
  aux[0] = eps;
  aux[1] = gmin;
  aux[2] = gmax;
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)
