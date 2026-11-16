// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, the split-population (cure) model of recidivism (Schmidt
// and Witte 1989): a logistic share who would ever fail and, among them,
// a Weibull, exponential or lognormal time to failure, by maximum
// likelihood with the analytic gradient.
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "kutil.h"
#include "qn.h"

namespace {

constexpr double ir2pi = 0.39894228040143267794;

struct cm {
  const double *t, *ev, *x;
  int n, p, dist;  // dist: 0 Weibull, 1 exponential, 2 lognormal
};

double xb(const cm &m, const double *b, int i) {
  double e = b[0];
  for (int j = 0; j < m.p; ++j) e += b[j + 1] * m.x[static_cast<std::size_t>(j) * m.n + i];
  return e;
}

// The negative log-likelihood at u = (gamma, beta, log shape) and its
// gradient.
double cnll(const double *u, double *g, void *ctx) {
  const cm &m = *static_cast<const cm *>(ctx);
  const int q = m.p + 1, np = 2 * q + (m.dist != 1);
  if (g != nullptr) {
    for (int k = 0; k < np; ++k) g[k] = 0.0;
  }
  const double ls = m.dist == 1 ? 0.0 : u[2 * q], sh = std::exp(ls);
  kutil::ksum s;
  for (int i = 0; i < m.n; ++i) {
    const double dl = 1.0 / (1.0 + std::exp(-xb(m, u, i))), eta = xb(m, u + q, i);
    const double lt = std::log(m.t[i]);
    double lf = 0.0, S = 0.0, fe = 0.0, fs = 0.0, Se = 0.0, Ss = 0.0;
    if (m.dist == 2) {
      const double w = (lt - eta) / sh, ph = ir2pi * std::exp(-0.5 * w * w);
      lf = -lt - ls - 0.5 * w * w - 0.91893853320467274178;
      S = 0.5 * std::erfc(w / 1.41421356237309504880);
      fe = w / sh;
      fs = -1.0 + w * w;
      Se = ph / sh;
      Ss = ph * w;
    } else {
      const double lz = sh * (lt - eta), z = std::exp(lz);
      lf = ls + (sh - 1) * lt - sh * eta - z;
      S = std::exp(-z);
      fe = -sh + sh * z;
      fs = 1 + lz - z * lz;
      Se = S * sh * z;
      Ss = -S * z * lz;
    }
    double li = 0.0, dg = 0.0, de = 0.0, ds = 0.0;
    if (m.ev[i] > 0) {
      li = std::log(dl) + lf;
      dg = 1 - dl;
      de = fe;
      ds = fs;
    } else {
      const double L = 1 - dl + dl * S;
      li = std::log(L);
      dg = dl * (1 - dl) * (S - 1) / L;
      de = dl * Se / L;
      ds = dl * Ss / L;
    }
    if (!std::isfinite(li)) return HUGE_VAL;
    s.add(li);
    if (g == nullptr) continue;
    g[0] -= dg;
    g[q] -= de;
    for (int j = 0; j < m.p; ++j) {
      const double xv = m.x[static_cast<std::size_t>(j) * m.n + i];
      g[j + 1] -= dg * xv;
      g[q + j + 1] -= de * xv;
    }
    if (m.dist != 1) g[2 * q] -= ds;
  }
  const double v = -s.val();
  return std::isfinite(v) ? v : HUGE_VAL;
}

}  // namespace

int clab_spcure(const double *t, const double *ev, const double *x, int n, int p, int dist, int maxit, double *par,
                double *cov, double *st) CLAB_TRY {
  if (n < 2 || p < 0 || dist < 0 || dist > 2 || maxit < 1) return CLAB_EARG;
  double ne = 0, sl = 0, sl2 = 0;
  for (int i = 0; i < n; ++i) {
    if (!kutil::pos(t[i]) || !(ev[i] == 0 || ev[i] == 1)) return CLAB_EDOM;
    if (ev[i] > 0) {
      ne += 1;
      sl += std::log(t[i]);
      sl2 += std::log(t[i]) * std::log(t[i]);
    }
  }
  if (ne < 1 || ne >= n) return CLAB_EDOM;
  const int q = p + 1, np = 2 * q + (dist != 1);
  cm m{t, ev, x, n, p, dist};
  // start: the share failing, the mean log time of the failures, shape 1
  // (the standard deviation of their log times for the lognormal)
  for (int k = 0; k < np; ++k) par[k] = 0.0;
  const double sh = std::fmin(0.95, ne / n + 0.05);
  par[0] = std::log(sh / (1 - sh));
  const double ml = sl / ne;
  par[q] = ml;
  if (dist == 2) par[2 * q] = 0.5 * std::log(std::fmax(sl2 / ne - ml * ml, 1e-4));
  std::vector<double> w(static_cast<std::size_t>(np) * np + 5 * static_cast<std::size_t>(np));
  double fv = 0.0;
  int it = 0;
  const int c = qn::min(np, par, cnll, &m, maxit, 1e-8, &fv, &it, w.data(), static_cast<double>(n) * (p + 4));
  if (kutil::stop(0.0)) return CLAB_EINTR;
  if (c < 0 || !std::isfinite(fv)) return CLAB_EDOM;
  std::vector<double> hw(3 * static_cast<std::size_t>(np));
  qn::hcov(np, par, cnll, &m, cov, hw.data());
  // st: log-likelihood, AIC, steps, converged, the mean share failing
  double md = 0;
  for (int i = 0; i < n; ++i) md += 1.0 / (1.0 + std::exp(-xb(m, par, i)));
  st[0] = -fv;
  st[1] = 2 * fv + 2 * np;
  st[2] = it;
  st[3] = c;
  st[4] = md / n;
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)
