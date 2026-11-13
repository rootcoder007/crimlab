// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, Hawkes processes with marks, types and places: the temporal
// ETAS model with the modified Omori kernel (clab_etasf), the multivariate
// Hawkes process with exponential kernels (clab_mvhwk) and the
// space-time Hawkes process with a Gaussian spread (clab_sthwk). Each
// gives the log-likelihood with its analytic gradient at given
// parameters (maxit 0) or fits by BFGS (qn.cpp). No allocation: the
// caller passes the work arrays.
#include <cmath>
#include <cstddef>

#include "crimlab.h"
#include "kutil.h"
#include "qn.h"

namespace {

constexpr double pi = 3.14159265358979323846;

bool timesok(const double *t, int n, double hor) {
  if (n < 1 || !kutil::pos(hor)) return false;
  for (int i = 0; i < n; ++i) {
    if (!(t[i] >= 0 && t[i] <= hor) || (i > 0 && t[i] < t[i - 1])) return false;
  }
  return true;
}

// ------------------------------------------------------------ ETAS

struct em {
  const double *t, *dm;
  int n;
  double hor;
  bool af;  // alpha free
};

// -log L at th = (mu, K, alpha, c, p) and its gradient; NaN outside
double enll(const em &m, const double *th, double *g) {
  const double mu = th[0], k = th[1], al = th[2], c = th[3], p = th[4];
  if (!kutil::pos(mu) || !kutil::nonneg(k) || !std::isfinite(al) || !kutil::pos(c) || !(p > 1) ||
      !std::isfinite(p)) {
    return kutil::nan();
  }
  double gr[5] = {0, 0, 0, 0, 0}, ll = 0;
  int lo = 0;
  for (int i = 0; i < m.n; ++i) {
    if (i > 0 && m.t[i] > m.t[i - 1]) lo = i;
    double s = 0, sa = 0, sc = 0, sp = 0;
    for (int j = 0; j < lo; ++j) {
      const double w = m.t[i] - m.t[j] + c, e = std::exp(al * m.dm[j]), h = e * std::pow(w, -p);
      s += h;
      sa += m.dm[j] * h;
      sc += -p * h / w;
      sp += -std::log(w) * h;
    }
    const double lam = mu + k * s;
    ll += std::log(lam);
    gr[0] -= 1 / lam;
    gr[1] -= s / lam;
    gr[2] -= k * sa / lam;
    gr[3] -= k * sc / lam;
    gr[4] -= k * sp / lam;
  }
  const double q = p - 1, cq = std::pow(c, -q), lc = std::log(c);
  double comp = mu * m.hor;
  gr[0] += m.hor;
  for (int j = 0; j < m.n; ++j) {
    const double w = m.hor - m.t[j] + c, wq = std::pow(w, -q), e = std::exp(al * m.dm[j]);
    const double f = (cq - wq) / q;
    comp += k * e * f;
    gr[1] += e * f;
    gr[2] += k * m.dm[j] * e * f;
    gr[3] += k * e * (-cq / c + wq / w);
    gr[4] += k * e * ((std::log(w) * wq - lc * cq) / q - f / q);
  }
  if (!m.af) gr[2] = 0;
  for (int j = 0; j < 5; ++j) g[j] = gr[j];
  return comp - ll;
}

// u: log mu, log K, alpha, log c, log(p - 1)
void eth(const double *u, double *th) {
  th[0] = std::exp(u[0]);
  th[1] = std::exp(u[1]);
  th[2] = u[2];
  th[3] = std::exp(u[3]);
  th[4] = 1 + std::exp(u[4]);
}

double efu(const double *u, double *g, void *ctx) {
  const em &m = *static_cast<const em *>(ctx);
  double th[5];
  eth(u, th);
  const double v = enll(m, th, g);
  g[0] *= th[0];
  g[1] *= th[1];
  g[3] *= th[3];
  g[4] *= th[4] - 1;
  return v;
}

double efth(const double *th, double *g, void *ctx) { return enll(*static_cast<const em *>(ctx), th, g); }

// the same with alpha fixed at 0: (mu, K, c, p)
double efth4(const double *t4, double *g4, void *ctx) {
  const double th[5] = {t4[0], t4[1], 0, t4[2], t4[3]};
  double g[5];
  const double v = enll(*static_cast<const em *>(ctx), th, g);
  g4[0] = g[0];
  g4[1] = g[1];
  g4[2] = g[3];
  g4[3] = g[4];
  return v;
}

// ------------------------------------------------------------ multivariate

struct mm {
  const double *t;
  const int *ty;
  int n, k;
  double hor, beta;  // beta <= 0: estimated (the last parameter)
  double *r, *q;     // 2 k scratch each
  int *pend;
};

double mnll(const mm &m, const double *th, double *g) {
  const int k = m.k, kk = k * k, np = k + kk + (m.beta > 0 ? 0 : 1);
  const double b = m.beta > 0 ? m.beta : th[k + kk];
  if (!kutil::pos(b)) return kutil::nan();
  for (int j = 0; j < k + kk; ++j) {
    if (!kutil::nonneg(th[j]) || (j < k && !(th[j] > 0))) return kutil::nan();
  }
  const double *mu = th, *a = th + k;
  for (int j = 0; j < np; ++j) g[j] = 0;
  for (int l = 0; l < k; ++l) {
    m.r[l] = m.q[l] = 0;
    m.pend[l] = 0;
  }
  double ll = 0, tl = 0, gb = 0;
  for (int i = 0; i < m.n; ++i) {
    const double ti = m.t[i];
    if (ti > tl) {
      const double dt = ti - tl, e = std::exp(-b * dt);
      for (int l = 0; l < k; ++l) {
        m.q[l] = e * (m.q[l] + dt * (m.r[l] + m.pend[l]));
        m.r[l] = e * (m.r[l] + m.pend[l]);
        m.pend[l] = 0;
      }
    }
    const int c = m.ty[i];
    double lam = mu[c];
    for (int l = 0; l < k; ++l) lam += b * a[c + l * k] * m.r[l];
    ll += std::log(lam);
    g[c] -= 1 / lam;
    for (int l = 0; l < k; ++l) {
      g[k + c + l * k] -= b * m.r[l] / lam;
      gb -= a[c + l * k] * (m.r[l] - b * m.q[l]) / lam;
    }
    ++m.pend[c];
    tl = ti;
  }
  double comp = 0;
  for (int c = 0; c < k; ++c) {
    comp += mu[c] * m.hor;
    g[c] += m.hor;
  }
  for (int j = 0; j < m.n; ++j) {
    const double d = m.hor - m.t[j], gg = -std::expm1(-b * d), dg = d * std::exp(-b * d);
    const int l = m.ty[j];
    for (int c = 0; c < k; ++c) {
      comp += a[c + l * k] * gg;
      g[k + c + l * k] += gg;
      gb += a[c + l * k] * dg;
    }
  }
  if (m.beta <= 0) g[k + kk] = gb;
  return comp - ll;
}

double mfu(const double *u, double *g, void *ctx) {
  const mm &m = *static_cast<const mm *>(ctx);
  const int np = m.k + m.k * m.k + (m.beta > 0 ? 0 : 1);
  double *th = m.r + m.k;  // the second half of r's scratch: np values
  for (int j = 0; j < np; ++j) th[j] = std::exp(u[j]);
  const double v = mnll(m, th, g);
  for (int j = 0; j < np; ++j) g[j] *= std::exp(u[j]);
  return v;
}

double mfth(const double *th, double *g, void *ctx) { return mnll(*static_cast<const mm *>(ctx), th, g); }

// ------------------------------------------------------------ space-time

struct sm {
  const double *t, *x, *y;
  int n;
  double hor, area;
};

double snll(const sm &m, const double *th, double *g) {
  const double mu = th[0], eta = th[1], b = th[2], s = th[3];
  if (!kutil::pos(mu) || !kutil::nonneg(eta) || !kutil::pos(b) || !kutil::pos(s)) return kutil::nan();
  const double s2 = s * s, nrm = 1 / (2 * pi * s2);
  double gr[4] = {0, 0, 0, 0}, ll = 0;
  int lo = 0;
  for (int i = 0; i < m.n; ++i) {
    if (i > 0 && m.t[i] > m.t[i - 1]) lo = i;
    double v = 0, vb = 0, vs = 0;
    for (int j = 0; j < lo; ++j) {
      const double d = m.t[i] - m.t[j], dx = m.x[i] - m.x[j], dy = m.y[i] - m.y[j], r2 = dx * dx + dy * dy;
      const double e = std::exp(-b * d), ph = nrm * std::exp(-r2 / (2 * s2));
      v += b * e * ph;
      vb += e * (1 - b * d) * ph;
      vs += b * e * ph * (r2 / (s2 * s) - 2 / s);
    }
    const double lam = mu + eta * v;
    ll += std::log(lam);
    gr[0] -= 1 / lam;
    gr[1] -= v / lam;
    gr[2] -= eta * vb / lam;
    gr[3] -= eta * vs / lam;
  }
  double comp = mu * m.hor * m.area;
  gr[0] += m.hor * m.area;
  for (int j = 0; j < m.n; ++j) {
    const double d = m.hor - m.t[j];
    comp += eta * -std::expm1(-b * d);
    gr[1] += -std::expm1(-b * d);
    gr[2] += eta * d * std::exp(-b * d);
  }
  for (int j = 0; j < 4; ++j) g[j] = gr[j];
  return comp - ll;
}

// u: log mu, logit eta, log beta, log sigma
void sth(const double *u, double *th) {
  th[0] = std::exp(u[0]);
  th[1] = 1 / (1 + std::exp(-u[1]));
  th[2] = std::exp(u[2]);
  th[3] = std::exp(u[3]);
}

double sfu(const double *u, double *g, void *ctx) {
  const sm &m = *static_cast<const sm *>(ctx);
  double th[4];
  sth(u, th);
  const double v = snll(m, th, g);
  g[0] *= th[0];
  g[1] *= th[1] * (1 - th[1]);
  g[2] *= th[2];
  g[3] *= th[3];
  return v;
}

double sfth(const double *th, double *g, void *ctx) { return snll(*static_cast<const sm *>(ctx), th, g); }

// a fit from u (overwritten): keeps the better in best (theta via conv)
struct best {
  double u[128] = {};
  double v = kutil::nan();
  int iter = 0, conv = 0;
};

void try1(int np, double *u, qn::fn f, void *ctx, int maxit, double *work, best &b) {
  double v;
  int it;
  const int c = qn::min(np, u, f, ctx, maxit, 1e-7, &v, &it, work);
  if (c >= 0 && std::isfinite(v) && (std::isnan(b.v) || v < b.v)) {
    for (int j = 0; j < np; ++j) b.u[j] = u[j];
    b.v = v;
    b.iter = it;
    b.conv = c;
  }
}

}  // namespace

extern "C" {

int clab_etasf(const double *t, const double *dm, int n, double hor, int afree,
               const double *start, int maxit, double *th, double *cov, double *st,
               double *grad) {
  if (!timesok(t, n, hor) || maxit < 0 || (maxit == 0 && !start)) return CLAB_EARG;
  em m{t, dm, n, hor, afree != 0};
  double g[5], work[5 * 5 + 5 * 5], u[5];
  best b;
  if (maxit == 0) {
    for (int j = 0; j < 5; ++j) th[j] = start[j];
    if (!afree) th[2] = 0;
    b.v = enll(m, th, g);
    if (std::isnan(b.v)) return CLAB_EDOM;
    b.iter = 0;
    b.conv = 0;
  } else {
    // two starts, or the one given
    const double s0[2][5] = {{n / (2 * hor), 0.2, 0.5, 0.01, 1.1}, {n / (2 * hor), 0.05, 0.2, 0.1, 1.5}};
    for (int r = 0; r < (start ? 1 : 2); ++r) {
      const double *s = start ? start : s0[r];
      if (!(s[0] > 0 && s[1] > 0 && s[3] > 0 && s[4] > 1 && std::isfinite(s[2]))) return CLAB_EDOM;
      u[0] = std::log(s[0]);
      u[1] = std::log(s[1]);
      u[2] = afree ? s[2] : 0;
      u[3] = std::log(s[3]);
      u[4] = std::log(s[4] - 1);
      try1(5, u, efu, &m, maxit, work, b);
    }
    if (std::isnan(b.v)) return CLAB_EDOM;
    eth(b.u, th);
    if (!afree) th[2] = 0;
  }
  enll(m, th, g);
  if (grad) {
    for (int j = 0; j < 5; ++j) grad[j] = -g[j];
  }
  if (afree) {
    qn::hcov(5, th, efth, &m, cov, work);
  } else {
    // alpha fixed at 0: the covariance of the other four, NaN for alpha
    double t4[4] = {th[0], th[1], th[3], th[4]}, c4[16];
    qn::hcov(4, t4, efth4, &m, c4, work);
    const int ix[4] = {0, 1, 3, 4};
    for (int i = 0; i < 25; ++i) cov[i] = kutil::nan();
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) cov[ix[i] + 5 * ix[j]] = c4[i + 4 * j];
    }
  }
  double me = 0;
  for (int j = 0; j < n; ++j) me += std::exp(th[2] * dm[j]) / n;
  st[0] = -b.v;
  st[1] = b.iter;
  st[2] = b.conv;
  st[3] = th[1] * me * std::pow(th[3], 1 - th[4]) / (th[4] - 1);
  return CLAB_OK;
}

int clab_mvhwk(const double *t, const int *ty, int n, int k, double hor, double beta,
               const double *start, int maxit, double *th, double *cov, double *st,
               double *grad, double *work, int *iwork) {
  if (!timesok(t, n, hor) || k < 1 || k > 10 || maxit < 0 || (maxit == 0 && !start)) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (ty[i] < 0 || ty[i] >= k) return CLAB_EARG;
  }
  const int kk = k * k, np = k + kk + (beta > 0 ? 0 : 1);
  double *r = work, *q = r + k + np, *u = q + k, *g = u + np, *w = g + np;  // w: np^2 + 5 np
  mm m{t, ty, n, k, hor, beta, r, q, iwork};
  best b;
  if (maxit == 0) {
    for (int j = 0; j < np; ++j) th[j] = start[j];
    b.v = mnll(m, th, g);
    if (std::isnan(b.v)) return CLAB_EDOM;
    b.iter = b.conv = 0;
  } else {
    // from start, or from mu_k = n_k / (2 hor), every A 0.1 and beta =
    // n / hor
    for (int j = 0; j < np; ++j) u[j] = 0;
    for (int i = 0; i < n; ++i) u[ty[i]] += 1;
    for (int j = 0; j < np; ++j) {
      const double s = start ? start[j] : (j < k ? std::fmax(u[j], 1.0) / (2 * hor) : (j < k + kk ? 0.1 : n / hor));
      if (!(s > 0)) return CLAB_EDOM;
      u[j] = std::log(s);
    }
    try1(np, u, mfu, &m, maxit, w, b);
    if (std::isnan(b.v)) return CLAB_EDOM;
    for (int j = 0; j < np; ++j) th[j] = std::exp(b.u[j]);
    mnll(m, th, g);
  }
  if (grad) {
    for (int j = 0; j < np; ++j) grad[j] = -g[j];
  }
  qn::hcov(np, th, mfth, &m, cov, w);
  // the spectral radius of A: power iteration on A + I (non-negative)
  double *v = u, *v2 = g, lam = 0;
  for (int c = 0; c < k; ++c) v[c] = 1;
  for (int it = 0; it < 100000; ++it) {
    double s = 0;
    for (int c = 0; c < k; ++c) {
      v2[c] = v[c];
      for (int l = 0; l < k; ++l) v2[c] += th[k + c + l * k] * v[l];
      s = std::fmax(s, v2[c]);
    }
    const double old = lam;
    lam = s;
    for (int c = 0; c < k; ++c) v[c] = v2[c] / s;
    if (std::fabs(lam - old) <= 1e-15 * lam) break;
  }
  st[0] = -b.v;
  st[1] = b.iter;
  st[2] = b.conv;
  st[3] = lam - 1;
  return CLAB_OK;
}

int clab_sthwk(const double *t, const double *x, const double *y, int n, double hor,
               double area, const double *start, int maxit, double *th, double *cov,
               double *st, double *grad, const double *qt, const double *qx, const double *qy,
               int nq, double *qlam) {
  if (!timesok(t, n, hor) || !kutil::pos(area) || maxit < 0 || (maxit == 0 && !start) || nq < 0) {
    return CLAB_EARG;
  }
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(x[i]) || !std::isfinite(y[i])) return CLAB_EDOM;
  }
  sm m{t, x, y, n, hor, area};
  double g[4], work[4 * 4 + 5 * 4], u[4];
  best b;
  if (maxit == 0) {
    for (int j = 0; j < 4; ++j) th[j] = start[j];
    b.v = snll(m, th, g);
    if (std::isnan(b.v)) return CLAB_EDOM;
    b.iter = b.conv = 0;
  } else {
    const double rt = n / hor, sd = std::sqrt(area);
    const double s0[2][4] = {{rt / (2 * area), 0.4, rt, 0.05 * sd}, {rt / (2 * area), 0.1, 0.25 * rt, 0.2 * sd}};
    for (int r = 0; r < (start ? 1 : 2); ++r) {
      const double *s = start ? start : s0[r];
      if (!(s[0] > 0 && s[1] > 0 && s[1] < 1 && s[2] > 0 && s[3] > 0)) return CLAB_EDOM;
      u[0] = std::log(s[0]);
      u[1] = std::log(s[1] / (1 - s[1]));
      u[2] = std::log(s[2]);
      u[3] = std::log(s[3]);
      try1(4, u, sfu, &m, maxit, work, b);
    }
    if (std::isnan(b.v)) return CLAB_EDOM;
    sth(b.u, th);
    snll(m, th, g);
  }
  if (grad) {
    for (int j = 0; j < 4; ++j) grad[j] = -g[j];
  }
  qn::hcov(4, th, sfth, &m, cov, work);
  st[0] = -b.v;
  st[1] = b.iter;
  st[2] = b.conv;
  // the conditional intensity at the query points, from strictly earlier
  // events
  const double s2 = th[3] * th[3];
  for (int i = 0; i < nq; ++i) {
    double v = th[0];
    for (int j = 0; j < n && t[j] < qt[i]; ++j) {
      const double dx = qx[i] - x[j], dy = qy[i] - y[j];
      v += th[1] * th[2] * std::exp(-th[2] * (qt[i] - t[j])) * std::exp(-(dx * dx + dy * dy) / (2 * s2)) /
           (2 * pi * s2);
    }
    qlam[i] = v;
  }
  return CLAB_OK;
}

}  // extern "C"
