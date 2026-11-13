// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, the temporal Hawkes process: the log-likelihood and its
// analytic gradient (clab_hwkll), the maximum-likelihood fit (clab_hwkfit),
// the time-rescaling residuals (clab_hwkres) and simulation by the
// branching construction, in time or in space and time (clab_hwksim).
// Kernels: exponential (Ozaki's O(n) recursion), Weibull and Lomax (the
// direct O(n^2) sum); baselines: constant, or a log-linear trend with one
// cycle. No allocation: the caller passes the work arrays.
#include <cmath>
#include <cstddef>

#include "crimlab.h"
#include "kutil.h"
#include "qn.h"

namespace {

constexpr double pi = 3.14159265358979323846;

struct hm {
  const double *t;
  int n;
  double hor;
  int kind;  // 0 exponential, 1 Weibull, 2 Lomax
  int base;  // 0 constant, 1 trend and cycle
  double per;
};

int nbase(int base) { return base ? 4 : 1; }
int nkern(int kind) { return kind ? 2 : 1; }

// the baseline at t and the terms of its log: 1, t / hor, sin, cos
double nu(const hm &m, const double *a, double t, double *f) {
  f[0] = 1;
  if (!m.base) return std::exp(a[0]);
  const double w = 2 * pi / m.per;
  f[1] = t / m.hor;
  f[2] = std::sin(w * t);
  f[3] = std::cos(w * t);
  return std::exp(a[0] + a[1] * f[1] + a[2] * f[2] + a[3] * f[3]);
}

// The integral of the baseline over [lo, hi] and its derivatives in the
// baseline terms: exact for the constant one, otherwise Simpson's rule on
// steps of at most per / 32 (an even number of them, at least 2).
double nuint(const hm &m, const double *a, double lo, double hi, double *da) {
  const int nb = nbase(m.base);
  double f[4];
  if (!m.base) {
    const double v = std::exp(a[0]) * (hi - lo);
    da[0] = v;
    return v;
  }
  int k = 2 * static_cast<int>(std::ceil((hi - lo) / (m.per / 16)));
  if (k < 2) k = 2;
  const double h = (hi - lo) / k;
  double s = 0;
  for (int j = 0; j < nb; ++j) da[j] = 0;
  for (int i = 0; i <= k; ++i) {
    const double w = (i == 0 || i == k) ? 1 : (i % 2 ? 4 : 2);
    const double v = w * nu(m, a, lo + i * h, f);
    s += v;
    for (int j = 0; j < nb; ++j) da[j] += v * f[j];
  }
  for (int j = 0; j < nb; ++j) da[j] *= h / 3;
  return s * h / 3;
}

// The kernel g at lag d > 0 (log), the distribution function G at d >= 0,
// and their derivatives in the two kernel parameters p0, p1 (Weibull:
// shape, scale; Lomax: shape, scale).
void kern(int kind, double p0, double p1, double d, double *lg, double *dlg, double *gg,
          double *dg) {
  if (kind == 1) {
    const double l = std::log(d / p1), z = std::exp(p0 * l), ez = std::exp(-z);
    *lg = std::log(p0) - std::log(p1) + (p0 - 1) * l - z;
    dlg[0] = 1 / p0 + l - z * l;
    dlg[1] = p0 * (z - 1) / p1;
    *gg = -std::expm1(-z);
    // past z = 700 the mass is all in: e^-z z underflows to 0
    dg[0] = z > 700 ? 0 : ez * z * l;
    dg[1] = z > 700 ? 0 : -ez * p0 * z / p1;
  } else {
    const double lw = std::log1p(d / p1), wa = std::exp(-p0 * lw);
    *lg = std::log(p0) - std::log(p1) - (p0 + 1) * lw;
    dlg[0] = 1 / p0 - lw;
    dlg[1] = -1 / p1 + (p0 + 1) * d / (p1 * p1 * (1 + d / p1));
    *gg = -std::expm1(-p0 * lw);
    dg[0] = wa * lw;
    dg[1] = -p0 * wa / (1 + d / p1) * d / (p1 * p1);
  }
}

bool okpar(const hm &m, const double *th) {
  const int nb = nbase(m.base), nk = nkern(m.kind);
  for (int j = 0; j < nb; ++j) {
    if (!std::isfinite(th[j])) return false;
  }
  if (!kutil::nonneg(th[nb])) return false;
  for (int j = 0; j < nk; ++j) {
    if (!kutil::pos(th[nb + 1 + j])) return false;
  }
  return true;
}

// The negative log-likelihood at th and, when g is not null, its gradient.
// NaN outside the parameter space.
double nll(const hm &m, const double *th, double *g) {
  if (!okpar(m, th)) return kutil::nan();
  const int nb = nbase(m.base), nk = nkern(m.kind), np = nb + 1 + nk;
  const double *a = th, eta = th[nb], p0 = th[nb + 1], p1 = nk > 1 ? th[nb + 2] : 0;
  double gr[7] = {0, 0, 0, 0, 0, 0, 0}, f[4], da[4];
  double ll = 0;
  // the exponential recursion over strictly earlier events: A = sum
  // e^(-b d), B = sum d e^(-b d); pend events at the last distinct time
  double ea = 0, eb = 0, tl = 0;
  int pend = 0, lo = 0;
  for (int i = 0; i < m.n; ++i) {
    const double ti = m.t[i];
    if (i > 0 && ti > m.t[i - 1]) lo = i;
    double s = 0, ds[2] = {0, 0};
    if (m.kind == 0) {
      if (i > 0 && ti > tl) {
        const double dt = ti - tl, e = std::exp(-p0 * dt);
        eb = e * (eb + dt * (ea + pend));
        ea = e * (ea + pend);
        pend = 0;
      }
      s = p0 * ea;
      ds[0] = ea - p0 * eb;
      ++pend;
      tl = ti;
    } else {
      for (int j = 0; j < lo; ++j) {
        double lg, dl[2], gg, dg[2];
        kern(m.kind, p0, p1, ti - m.t[j], &lg, dl, &gg, dg);
        const double v = std::exp(lg);
        if (v > 0) {
          s += v;
          ds[0] += v * dl[0];
          ds[1] += v * dl[1];
        }
      }
    }
    const double v = nu(m, a, ti, f), lam = v + eta * s;
    ll += std::log(lam);
    for (int j = 0; j < nb; ++j) gr[j] -= v * f[j] / lam;
    gr[nb] -= s / lam;
    for (int j = 0; j < nk; ++j) gr[nb + 1 + j] -= eta * ds[j] / lam;
  }
  // the compensator: the baseline's integral and each event's kernel mass
  // inside the window
  double comp = nuint(m, a, 0, m.hor, da);
  for (int j = 0; j < nb; ++j) gr[j] += da[j];
  for (int j = 0; j < m.n; ++j) {
    const double d = m.hor - m.t[j];
    double gg = 0, dg[2] = {0, 0};
    if (m.kind == 0) {
      const double e = std::exp(-p0 * d);
      gg = -std::expm1(-p0 * d);
      dg[0] = d * e;
    } else if (d > 0) {
      double lg, dl[2];
      kern(m.kind, p0, p1, d, &lg, dl, &gg, dg);
    }
    comp += eta * gg;
    gr[nb] += gg;
    for (int k = 0; k < nk; ++k) gr[nb + 1 + k] += eta * dg[k];
  }
  if (g) {
    for (int j = 0; j < np; ++j) g[j] = gr[j];
  }
  return comp - ll;
}

// the objective of the fit in unconstrained u: the baseline terms as they
// are, eta = 1 / (1 + e^-u), kernel parameters e^u
void totheta(const hm &m, const double *u, double *th) {
  const int nb = nbase(m.base), nk = nkern(m.kind);
  for (int j = 0; j < nb; ++j) th[j] = u[j];
  th[nb] = 1 / (1 + std::exp(-u[nb]));
  for (int j = 0; j < nk; ++j) th[nb + 1 + j] = std::exp(u[nb + 1 + j]);
}

double fu(const double *u, double *g, void *ctx) {
  const hm &m = *static_cast<const hm *>(ctx);
  const int nb = nbase(m.base), np = nb + 1 + nkern(m.kind);
  double th[7];
  totheta(m, u, th);
  const double v = nll(m, th, g);
  for (int j = nb + 1; j < np; ++j) g[j] *= th[j];
  g[nb] *= th[nb] * (1 - th[nb]);
  return v;
}

double fth(const double *th, double *g, void *ctx) {
  return nll(*static_cast<const hm *>(ctx), th, g);
}

// One fit from the start th (in theta); u, the result in theta, nll and
// steps into best when it improves on it.
struct run {
  double th[7];
  double v;
  int iter, conv;
};

void fit1(const hm &m, const double *th0, int maxit, run &best) {
  const int nb = nbase(m.base), nk = nkern(m.kind), np = nb + 1 + nk;
  double u[7], work[7 * 7 + 5 * 7];
  for (int j = 0; j < nb; ++j) u[j] = th0[j];
  u[nb] = std::log(th0[nb] / (1 - th0[nb]));
  for (int j = 0; j < nk; ++j) u[nb + 1 + j] = std::log(th0[nb + 1 + j]);
  double v;
  int it;
  const int c = qn::min(np, u, fu, const_cast<hm *>(&m), maxit, 1e-7, &v, &it, work);
  if (c >= 0 && std::isfinite(v) && (std::isnan(best.v) || v < best.v)) {
    totheta(m, u, best.th);
    best.v = v;
    best.iter = it;
    best.conv = c;
  }
}

bool sorted(const double *t, int n, double hor) {
  for (int i = 0; i < n; ++i) {
    if (!(t[i] >= 0 && t[i] <= hor) || (i > 0 && t[i] < t[i - 1])) return false;
  }
  return true;
}

bool okmodel(const double *t, int n, double hor, int kind, int base, double per) {
  return n >= 1 && kutil::pos(hor) && kind >= 0 && kind <= 2 && base >= 0 && base <= 1 &&
         (!base || kutil::pos(per)) && sorted(t, n, hor);
}

// a standard normal pair from two uniforms (Box and Muller 1958)
void bm(kutil::rng &r, double *z0, double *z1) {
  const double a = std::sqrt(-2 * std::log(1 - r.unif())), b = 2 * pi * r.unif();
  *z0 = a * std::cos(b);
  *z1 = a * std::sin(b);
}

}  // namespace

extern "C" {

int clab_hwkll(const double *t, int n, double hor, int kind, int base, double per,
               const double *th, double *ll, double *grad) {
  if (!okmodel(t, n, hor, kind, base, per)) return CLAB_EARG;
  hm m{t, n, hor, kind, base, per};
  if (!okpar(m, th)) return CLAB_EDOM;
  const double v = nll(m, th, grad);
  *ll = -v;
  if (grad) {
    for (int j = 0; j < nbase(base) + 1 + nkern(kind); ++j) grad[j] = -grad[j];
  }
  return CLAB_OK;
}

int clab_hwkfit(const double *t, int n, double hor, int kind, int base, double per,
                const double *start, int maxit, double *th, double *cov, double *st) {
  if (!okmodel(t, n, hor, kind, base, per) || n < 2 || maxit < 1) return CLAB_EARG;
  hm m{t, n, hor, kind, base, per};
  const int nb = nbase(base), np = nb + 1 + nkern(kind);
  run best;
  best.v = kutil::nan();
  if (start) {
    if (!okpar(m, start) || !(start[nb] > 0 && start[nb] < 1)) return CLAB_EDOM;
    fit1(m, start, maxit, best);
  } else {
    const double mg = hor / n;
    double s[7] = {std::log(0.6 * n / hor), 0, 0, 0, 0, 0, 0};
    // the exponential fit from two starts; for the others it is a third
    // start, the shape-1 member of the family
    hm me = m;
    me.kind = 0;
    run ex;
    ex.v = kutil::nan();
    s[nb] = 0.4;
    s[nb + 1] = 1 / mg;
    fit1(me, s, maxit, ex);
    s[nb] = 0.1;
    s[nb + 1] = 0.25 / mg;
    fit1(me, s, maxit, ex);
    if (kind == 0) {
      best = ex;
    } else {
      s[nb] = 0.4;
      s[nb + 1] = kind == 1 ? 1.5 : 2.5;
      s[nb + 2] = kind == 1 ? 1.2 * mg : 5 * mg;
      fit1(m, s, maxit, best);
      if (!std::isnan(ex.v)) {
        for (int j = 0; j <= nb; ++j) s[j] = ex.th[j];
        s[nb] = std::fmin(std::fmax(s[nb], 1e-6), 1 - 1e-6);
        s[nb + 1] = kind == 1 ? 1 : 10;
        s[nb + 2] = kind == 1 ? 1 / ex.th[nb + 1] : 10 / ex.th[nb + 1];
        fit1(m, s, maxit, best);
      }
    }
  }
  if (std::isnan(best.v)) return CLAB_EDOM;
  for (int j = 0; j < np; ++j) th[j] = best.th[j];
  double work[21];
  qn::hcov(np, th, fth, &m, cov, work);
  st[0] = -best.v;
  st[1] = best.iter;
  st[2] = best.conv;
  return CLAB_OK;
}

int clab_hwkres(const double *t, int n, double hor, int kind, int base, double per,
                const double *th, double *tau, double *ks, double *work) {
  if (!okmodel(t, n, hor, kind, base, per)) return CLAB_EARG;
  hm m{t, n, hor, kind, base, per};
  if (!okpar(m, th)) return CLAB_EDOM;
  const int nb = nbase(base);
  const double eta = th[nb], p0 = th[nb + 1], p1 = kind ? th[nb + 2] : 0;
  double da[4], cum = 0, prev = 0, ea = 0, tl = 0;
  int pend = 0, lo = 0;
  bool tie = false;
  for (int i = 0; i < n; ++i) {
    const double ti = t[i];
    if (i > 0 && ti > t[i - 1]) lo = i;
    tie = tie || (i > 0 && ti == t[i - 1]);
    cum += nuint(m, th, i ? t[i - 1] : 0, ti, da);
    double s = 0;
    if (kind == 0) {
      if (i > 0 && ti > tl) {
        ea = std::exp(-p0 * (ti - tl)) * (ea + pend);
        pend = 0;
      }
      s = lo - ea;
      ++pend;
      tl = ti;
    } else {
      for (int j = 0; j < lo; ++j) {
        double lg, dl[2], gg, dg[2];
        kern(kind, p0, p1, ti - t[j], &lg, dl, &gg, dg);
        s += gg;
      }
    }
    const double big = cum + eta * s;
    tau[i] = big - prev;
    prev = big;
  }
  for (int i = 0; i < n; ++i) work[i] = -std::expm1(-tau[i]);
  ks[0] = qn::ksd(work, n);
  ks[1] = qn::kspv(ks[0], n, !tie, work);
  return CLAB_OK;
}

int clab_hwksim(double hor, int kind, int base, double per, const double *th,
                const double *win, double sig, unsigned int seed, int cap, double *t,
                double *x, double *y, int *gen, int *n, double *work, int *iwork) {
  if (!kutil::pos(hor) || kind < 0 || kind > 2 || base < 0 || base > 1 ||
      (base && !kutil::pos(per)) || cap < 1) {
    return CLAB_EARG;
  }
  hm m{nullptr, 0, hor, kind, base, per};
  if (!okpar(m, th)) return CLAB_EDOM;
  if (win && !(win[0] < win[1] && win[2] < win[3] && std::isfinite(win[0]) &&
               std::isfinite(win[1]) && std::isfinite(win[2]) && std::isfinite(win[3]) &&
               kutil::pos(sig))) {
    return CLAB_EDOM;
  }
  const int nb = nbase(base);
  const double eta = th[nb], p0 = th[nb + 1], p1 = kind ? th[nb + 2] : 0;
  kutil::rng r(seed);
  int k = 0;
  double f[4];
  // the immigrants: thinning against the largest the baseline can be
  const double top = base ? std::exp(th[0] + std::fmax(0.0, th[1]) + std::hypot(th[2], th[3]))
                          : std::exp(th[0]);
  for (double s = -std::log(1 - r.unif()) / top; s < hor; s += -std::log(1 - r.unif()) / top) {
    if (base && r.unif() * top > nu(m, th, s, f)) continue;
    if (k == cap) return CLAB_EDOM;
    t[k] = s;
    gen[k] = 0;
    if (win) {
      x[k] = win[0] + (win[1] - win[0]) * r.unif();
      y[k] = win[2] + (win[3] - win[2]) * r.unif();
    }
    ++k;
  }
  // each event's offspring: Poisson(eta) of them, at lags drawn from the
  // kernel by inversion, kept when inside the window
  for (int i = 0; i < k; ++i) {
    int c = 0;
    double pr = std::exp(-eta), cdf = pr;
    for (const double u = r.unif(); u > cdf && c < 100000; cdf += pr) pr *= eta / ++c;
    for (int j = 0; j < c; ++j) {
      const double u = 1 - r.unif();
      const double lag = kind == 0   ? -std::log(u) / p0
                         : kind == 1 ? p1 * std::pow(-std::log(u), 1 / p0)
                                     : p1 * std::expm1(-std::log(u) / p0);
      double cx = 0, cy = 0;
      if (win) {
        double z0, z1;
        bm(r, &z0, &z1);
        cx = x[i] + sig * z0;
        cy = y[i] + sig * z1;
      }
      if (t[i] + lag >= hor || (win && !(cx >= win[0] && cx <= win[1] && cy >= win[2] && cy <= win[3]))) {
        continue;
      }
      if (k == cap) return CLAB_EDOM;
      t[k] = t[i] + lag;
      gen[k] = gen[i] + 1;
      if (win) {
        x[k] = cx;
        y[k] = cy;
      }
      ++k;
    }
  }
  // in time order
  for (int i = 0; i < k; ++i) iwork[i] = i;
  kutil::hsort(iwork, k, [t](int a, int b) { return t[a] < t[b] || (t[a] == t[b] && a < b); });
  double *wt = work, *wx = work + cap, *wy = work + 2 * static_cast<std::size_t>(cap);
  int *wg = iwork + cap;
  for (int i = 0; i < k; ++i) {
    wt[i] = t[iwork[i]];
    wg[i] = gen[iwork[i]];
    if (win) {
      wx[i] = x[iwork[i]];
      wy[i] = y[iwork[i]];
    }
  }
  for (int i = 0; i < k; ++i) {
    t[i] = wt[i];
    gen[i] = wg[i];
    if (win) {
      x[i] = wx[i];
      y[i] = wy[i];
    }
  }
  *n = k;
  return CLAB_OK;
}

}  // extern "C"
