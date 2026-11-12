// SPDX-License-Identifier: AGPL-3.0-or-later
//
// libcrimlab, courts and corrections: group-based trajectory models for
// counts (Nagin), Poisson or zero-inflated Poisson within each group, fitted
// by EM with a damped Newton M-step.
#include <cstdint>

#include "crimlab.h"
#include "kutil.h"

namespace {

// log(1 + e^a) and log(e^a + e^b) without overflow
inline double softplus(double a) {
  return a > 0.0 ? a + std::log1p(std::exp(-a)) : std::log1p(std::exp(a));
}
inline double logistic(double a) {
  return a >= 0.0 ? 1.0 / (1.0 + std::exp(-a)) : std::exp(a) / (1.0 + std::exp(a));
}

struct dat {
  const double *y, *off, *x, *z;
  int n, nt, p, q;
  // y missing: NaN or negative
  bool miss(int i, int t) const {
    const double v = y[i + static_cast<std::int64_t>(n) * t];
    return !(v >= 0.0);
  }
  double yv(int i, int t) const { return y[i + static_cast<std::int64_t>(n) * t]; }
  double ov(int i, int t) const {
    return off == nullptr ? 0.0 : off[i + static_cast<std::int64_t>(n) * t];
  }
};

// The log-likelihood of one cell, with eta = off + x'beta, lam = e^eta and
// a = z'gamma (zero inflation rho = logistic(a); q = 0: Poisson):
//   y = 0:  log(e^a + e^-lam) - log(1 + e^a)
//   y > 0:  y eta - lam - log y! - log(1 + e^a)
// and, when d is given, its derivatives d[0] = dl/deta, d[1] = dl/da,
// d[2] = d2l/deta2, d[3] = d2l/deta da, d[4] = d2l/da2. With
// s = logistic(a + lam) (the chance that a zero is structural):
//   y = 0: dl/da = s - rho, dl/deta = -(1 - s) lam,
//          d2l/deta2 = s (1 - s) lam^2 - (1 - s) lam,
//          d2l/deta da = s (1 - s) lam, d2l/da2 = s (1 - s) - rho (1 - rho)
//   y > 0: dl/da = -rho, dl/deta = y - lam, d2l/deta2 = -lam,
//          d2l/da2 = -rho (1 - rho)
double cell(double y, double eta, double a, bool zip, double *d, bool lg = true) {
  const double lam = std::exp(eta);
  const double lyf = lg && y > 1.0 ? std::lgamma(y + 1.0) : 0.0;  // log y!
  if (!zip) {
    if (d != nullptr) {
      d[0] = y - lam;
      d[2] = -lam;
      d[1] = d[3] = d[4] = 0.0;
    }
    return y * eta - lam - lyf;
  }
  const double rho = logistic(a), sp = softplus(a);
  if (y == 0.0) {
    const double s = logistic(a + lam);
    if (d != nullptr) {
      d[0] = -(1.0 - s) * lam;
      d[1] = s - rho;
      d[2] = s * (1.0 - s) * lam * lam - (1.0 - s) * lam;
      d[3] = s * (1.0 - s) * lam;
      d[4] = s * (1.0 - s) - rho * (1.0 - rho);
    }
    const double hi = a > -lam ? a : -lam, lo = a > -lam ? -lam : a;
    return hi + std::log1p(std::exp(lo - hi)) - sp;
  }
  if (d != nullptr) {
    d[0] = y - lam;
    d[1] = -rho;
    d[2] = -lam;
    d[3] = 0.0;
    d[4] = -rho * (1.0 - rho);
  }
  return y * eta - lam - lyf - sp;
}

// the linear predictors of every period under th = (beta, gamma), without
// the offsets: le[t] = x_t' beta, la[t] = z_t' gamma
void lin(const dat &D, const double *th, double *le, double *la) {
  for (int t = 0; t < D.nt; ++t) {
    double e = 0.0, b = 0.0;
    for (int j = 0; j < D.p; ++j) e += D.x[t + static_cast<std::int64_t>(D.nt) * j] * th[j];
    for (int j = 0; j < D.q; ++j) b += D.z[t + static_cast<std::int64_t>(D.nt) * j] * th[D.p + j];
    le[t] = e;
    la[t] = b;
  }
}

// E-step: w (n x ng) posterior membership, returns the log-likelihood;
// work: 2 nt
double estep(const dat &D, const double *par, const double *prob, int ng, double *w,
             double *work) {
  const int m = D.p + D.q;
  double *le = work, *la = work + D.nt;
  for (int k = 0; k < ng; ++k) {
    lin(D, par + k * m, le, la);
    const double lp = std::log(prob[k]);
    for (int i = 0; i < D.n; ++i) {
      kutil::ksum s;
      for (int t = 0; t < D.nt; ++t) {
        if (D.miss(i, t)) continue;
        s.add(cell(D.yv(i, t), D.ov(i, t) + le[t], la[t], D.q > 0, nullptr));
      }
      w[i + static_cast<std::int64_t>(D.n) * k] = lp + s.val();
    }
  }
  kutil::ksum ll;
  for (int i = 0; i < D.n; ++i) {
    double mx = -HUGE_VAL;
    for (int k = 0; k < ng; ++k) {
      const double v = w[i + static_cast<std::int64_t>(D.n) * k];
      if (v > mx) mx = v;
    }
    double s = 0.0;
    for (int k = 0; k < ng; ++k) {
      double &v = w[i + static_cast<std::int64_t>(D.n) * k];
      v = std::exp(v - mx);
      s += v;
    }
    for (int k = 0; k < ng; ++k) w[i + static_cast<std::int64_t>(D.n) * k] /= s;
    ll.add(mx + std::log(s));
  }
  return ll.val();
}

// Q_k(th) = sum_i w_ik sum_t l_it(th) (log y! left out: it does not depend
// on th), with gradient g and Hessian h (m x m, full) when g is given. The
// design rows depend on t only, so the cell derivatives are first summed
// over units for each period (acc: 5 nt), then expanded: O(n nt + nt m^2).
// work: 2 nt + 5 nt.
double qfun(const dat &D, const double *wk, const double *th, double *g, double *h,
            double *work) {
  const int m = D.p + D.q, nt = D.nt;
  double *le = work, *la = work + nt, *acc = work + 2 * nt;
  lin(D, th, le, la);
  if (g != nullptr)
    for (int j = 0; j < 5 * nt; ++j) acc[j] = 0.0;
  kutil::ksum s;
  double d[5];
  for (int i = 0; i < D.n; ++i) {
    const double wi = wk[i];
    if (!(wi > 0.0)) continue;
    for (int t = 0; t < nt; ++t) {
      if (D.miss(i, t)) continue;
      s.add(wi * cell(D.yv(i, t), D.ov(i, t) + le[t], la[t], D.q > 0,
                      g != nullptr ? d : nullptr, false));
      if (g == nullptr) continue;
      for (int c = 0; c < 5; ++c) acc[c * nt + t] += wi * d[c];
    }
  }
  if (g == nullptr) return s.val();
  // design value j of period t: x_tj for j < p, z_t(j-p) after
  auto v = [&D](int t, int j) {
    return j < D.p ? D.x[t + static_cast<std::int64_t>(D.nt) * j]
                   : D.z[t + static_cast<std::int64_t>(D.nt) * (j - D.p)];
  };
  for (int j = 0; j < m; ++j) {
    double gj = 0.0;
    for (int t = 0; t < nt; ++t) gj += acc[(j < D.p ? 0 : 1) * nt + t] * v(t, j);
    g[j] = gj;
    for (int l = 0; l <= j; ++l) {
      const int c = j < D.p && l < D.p ? 2 : (j >= D.p && l >= D.p ? 4 : 3);
      double hj = 0.0;
      for (int t = 0; t < nt; ++t) hj += acc[c * nt + t] * v(t, j) * v(t, l);
      h[j + m * l] = h[l + m * j] = hj;
    }
  }
  return s.val();
}

// Solve (-H + mu I) dx = g by Cholesky in a (m x m scratch); false when not
// positive definite
bool dsolve(const double *h, const double *g, double mu, int m, double *a, double *dx) {
  for (int j = 0; j < m; ++j)
    for (int i = j; i < m; ++i) a[i + m * j] = -h[i + m * j] + (i == j ? mu : 0.0);
  for (int j = 0; j < m; ++j) {
    double dd = a[j + m * j];
    for (int k = 0; k < j; ++k) dd -= a[j + m * k] * a[j + m * k];
    if (!(dd > 0.0) || !std::isfinite(dd)) return false;
    dd = std::sqrt(dd);
    a[j + m * j] = dd;
    for (int i = j + 1; i < m; ++i) {
      double e = a[i + m * j];
      for (int k = 0; k < j; ++k) e -= a[i + m * k] * a[j + m * k];
      a[i + m * j] = e / dd;
    }
  }
  for (int i = 0; i < m; ++i) {
    double e = g[i];
    for (int k = 0; k < i; ++k) e -= a[i + m * k] * dx[k];
    dx[i] = e / a[i + m * i];
  }
  for (int i = m - 1; i >= 0; --i) {
    double e = dx[i];
    for (int k = i + 1; k < m; ++k) e -= a[k + m * i] * dx[k];
    dx[i] = e / a[i + m * i];
  }
  return true;
}

// M-step for one group: up to 4 damped Newton steps of ascent on Q_k (a
// generalised EM: every accepted step raises Q_k, so the log-likelihood
// never decreases, and a fixed point has a zero gradient of Q_k). work:
// 2 m^2 + 3 m + 7 nt.
void mstep(const dat &D, const double *wk, double *th, double *work) {
  const int m = D.p + D.q;
  double *g = work, *h = g + m, *a = h + m * m, *dx = a + m * m, *tr = dx + m;
  double *qw = tr + m;  // 7 nt for qfun
  double mu = 0.0;
  for (int it = 0; it < 4; ++it) {
    const double q0 = qfun(D, wk, th, g, h, qw);
    double gmax = 0.0, hmax = 0.0;
    for (int j = 0; j < m; ++j) {
      gmax = std::fabs(g[j]) > gmax ? std::fabs(g[j]) : gmax;
      hmax = std::fabs(h[j + m * j]) > hmax ? std::fabs(h[j + m * j]) : hmax;
    }
    if (gmax <= 1e-10 * (1.0 + std::fabs(q0))) return;
    bool ok = false;
    for (int tries = 0; tries < 60 && !ok; ++tries) {
      if (dsolve(h, g, mu, m, a, dx)) {
        for (int j = 0; j < m; ++j) tr[j] = th[j] + dx[j];
        const double q1 = qfun(D, wk, tr, nullptr, nullptr, qw);
        if (std::isfinite(q1) && q1 >= q0) {
          for (int j = 0; j < m; ++j) th[j] = tr[j];
          ok = true;
          mu *= 0.1;
          if (mu < 1e-12 * (1.0 + hmax)) mu = 0.0;
          if (q1 - q0 <= 1e-14 * (1.0 + std::fabs(q0))) return;
          break;
        }
      }
      mu = mu > 0.0 ? 10.0 * mu : 1e-8 * (1.0 + hmax);
    }
    if (!ok) return;
  }
}

int check(const dat &D, int ng) {
  if (D.n < 1 || D.nt < 1 || D.p < 1 || D.q < 0 || ng < 1) return CLAB_EARG;
  for (int i = 0; i < D.n; ++i)
    for (int t = 0; t < D.nt; ++t) {
      const double v = D.yv(i, t);
      if (std::isinf(v) || (v >= 0.0 && v != std::floor(v))) return CLAB_EDOM;
      if (!std::isfinite(D.ov(i, t))) return CLAB_EDOM;
    }
  for (int t = 0; t < D.nt; ++t) {
    for (int j = 0; j < D.p; ++j)
      if (!std::isfinite(D.x[t + static_cast<std::int64_t>(D.nt) * j])) return CLAB_EDOM;
    for (int j = 0; j < D.q; ++j)
      if (!std::isfinite(D.z[t + static_cast<std::int64_t>(D.nt) * j])) return CLAB_EDOM;
  }
  return CLAB_OK;
}

}  // namespace

extern "C" {

int clab_zipll(const double *y, int n, int nt, const double *off,
               const double *x, int p, const double *z, int q, int ng,
               const double *par, const double *prob, double *w, double *ll,
               double *work) {
  const dat D{y, off, x, z, n, nt, p, q};
  const int rc = check(D, ng);
  if (rc != CLAB_OK) return rc;
  double s = 0.0;
  for (int k = 0; k < ng; ++k) {
    if (!(prob[k] > 0.0)) return CLAB_EDOM;
    s += prob[k];
  }
  if (std::fabs(s - 1.0) > 1e-9) return CLAB_EDOM;
  for (int j = 0; j < ng * (p + q); ++j)
    if (!std::isfinite(par[j])) return CLAB_EDOM;
  *ll = estep(D, par, prob, ng, w, work);
  return CLAB_OK;
}

int clab_zipem(const double *y, int n, int nt, const double *off,
               const double *x, int p, const double *z, int q, int ng,
               double *par, double *prob, double *w, double *ll, int maxit,
               double tol, int *iter, double *work) {
  const dat D{y, off, x, z, n, nt, p, q};
  if (maxit < 1 || !(tol > 0.0)) return CLAB_EARG;
  int rc = clab_zipll(y, n, nt, off, x, p, z, q, ng, par, prob, w, ll, work);
  if (rc != CLAB_OK) return rc;
  const int m = p + q;
  double cur = *ll;
  int it = 0;
  for (; it < maxit; ++it) {
    for (int k = 0; k < ng; ++k) {
      double s = 0.0;
      for (int i = 0; i < n; ++i) s += w[i + static_cast<std::int64_t>(n) * k];
      prob[k] = s / n;
      if (prob[k] > 0.0) mstep(D, w + static_cast<std::int64_t>(n) * k, par + k * m, work);
    }
    // a group that lost every unit keeps a tiny share so logs stay finite
    for (int k = 0; k < ng; ++k)
      if (!(prob[k] > 1e-300)) prob[k] = 1e-300;
    const double nl = estep(D, par, prob, ng, w, work);
    const double gain = nl - cur;
    cur = nl;
    if (gain <= tol * (1.0 + std::fabs(nl))) {
      ++it;
      break;
    }
  }
  *ll = cur;
  *iter = it;
  return CLAB_OK;
}

}  // extern "C"
