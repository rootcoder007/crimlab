// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, counts over time and across categories: ARIMA and seasonal
// ARIMA forecasts by exact maximum likelihood (the Kalman filter of
// Gardner, Harvey and Phillips 1980 on the state-space form of Harvey
// 1989, conditional sum of squares for the start, BFGS as R's optim()
// runs it), the Ornstein-Uhlenbeck (Langevin) fit and simulation of a
// count series and the Fokker-Planck evolution of its density, the
// Lotka-Volterra rates of a yearly series, and the correlation of counts
// per unit between categories.
#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>
#include <limits>
#include <map>

#include "kutil.h"
#include "svk.h"

namespace svk {

namespace {

const double NaN = std::numeric_limits<double>::quiet_NaN();

std::size_t z(int i) { return static_cast<std::size_t>(i); }

double mean(const std::vector<double> &v) {
  double s = 0;
  for (double e : v) s += e;
  return s / static_cast<double>(v.size());
}

double median(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  const std::size_t n = v.size();
  return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

// --- ARIMA ---------------------------------------------------------------

struct Arma {
  int mp, mq, msp, msq, ns, d, sd;
  int np() const { return mp + mq + msp + msq; }
};

void partrans(int p, const double *raw, double *nw) {
  std::vector<double> work(z(p));
  for (int j = 0; j < p; ++j) work[z(j)] = nw[j] = std::tanh(raw[j]);
  for (int j = 1; j < p; ++j) {
    const double a = nw[j];
    for (int k = 0; k < j; ++k) work[z(k)] -= a * nw[j - k - 1];
    for (int k = 0; k < j; ++k) nw[k] = work[z(k)];
  }
}

void invpartrans(int p, const double *phi, double *nw) {
  std::vector<double> work(z(p));
  for (int j = 0; j < p; ++j) work[z(j)] = nw[j] = phi[j];
  for (int j = p - 1; j > 0; --j) {
    const double a = nw[j];
    for (int k = 0; k < j; ++k) work[z(k)] = (nw[k] + a * nw[j - k - 1]) / (1 - a * a);
    for (int k = 0; k < j; ++k) nw[k] = work[z(k)];
  }
  for (int j = 0; j < p; ++j) nw[j] = std::atanh(nw[j]);
}

void transpars(const std::vector<double> &in, const Arma &a, bool trans, std::vector<double> &phi,
               std::vector<double> &theta) {
  const int p = a.mp + a.ns * a.msp, q = a.mq + a.ns * a.msq;
  std::vector<double> par = in;
  if (trans) {
    if (a.mp > 0) partrans(a.mp, in.data(), par.data());
    const int v = a.mp + a.mq;
    if (a.msp > 0) partrans(a.msp, in.data() + v, par.data() + v);
  }
  phi.assign(z(p), 0.0);
  theta.assign(z(q), 0.0);
  for (int i = 0; i < a.mp; ++i) phi[z(i)] = par[z(i)];
  for (int i = 0; i < a.mq; ++i) theta[z(i)] = par[z(i + a.mp)];
  if (a.ns > 0) {
    for (int j = 0; j < a.msp; ++j) {
      phi[z((j + 1) * a.ns - 1)] += par[z(j + a.mp + a.mq)];
      for (int i = 0; i < a.mp; ++i) phi[z((j + 1) * a.ns + i)] -= par[z(i)] * par[z(j + a.mp + a.mq)];
    }
    for (int j = 0; j < a.msq; ++j) {
      theta[z((j + 1) * a.ns - 1)] += par[z(j + a.mp + a.mq + a.msp)];
      for (int i = 0; i < a.mq; ++i) theta[z((j + 1) * a.ns + i)] += par[z(i + a.mp)] * par[z(j + a.mp + a.mq + a.msp)];
    }
  }
}

std::vector<double> undopars(const std::vector<double> &in, const Arma &a) {
  std::vector<double> r = in;
  if (a.mp > 0) partrans(a.mp, in.data(), r.data());
  const int v = a.mp + a.mq;
  if (a.msp > 0) partrans(a.msp, in.data() + v, r.data() + v);
  return r;
}

std::vector<double> invtrans(const std::vector<double> &in, const Arma &a) {
  std::vector<double> r = in;
  if (a.mp > 0) invpartrans(a.mp, in.data(), r.data());
  const int v = a.mp + a.mq;
  if (a.msp > 0) invpartrans(a.msp, in.data() + v, r.data() + v);
  return r;
}

double css(const std::vector<double> &y, const Arma &a, const std::vector<double> &phi,
           const std::vector<double> &theta, int ncond) {
  const int n = static_cast<int>(y.size()), p = static_cast<int>(phi.size()), q = static_cast<int>(theta.size());
  std::vector<double> w = y, res(z(n));
  for (int i = 0; i < a.d; ++i) {
    for (int l = n - 1; l > 0; --l) w[z(l)] -= w[z(l - 1)];
  }
  for (int i = 0; i < a.sd; ++i) {
    for (int l = n - 1; l >= a.ns; --l) w[z(l)] -= w[z(l - a.ns)];
  }
  double ssq = 0;
  int nu = 0;
  for (int l = ncond; l < n; ++l) {
    double t = w[z(l)];
    for (int j = 0; j < p; ++j) t -= phi[z(j)] * w[z(l - j - 1)];
    for (int j = 0; j < std::min(l - ncond, q); ++j) t -= theta[z(j)] * res[z(l - j - 1)];
    res[z(l)] = t;
    if (!std::isnan(t)) {
      ++nu;
      ssq += t * t;
    }
  }
  return ssq / nu;
}

// Roots of the polynomial c[0] + c[1] x + ... (Aberth's iteration).
std::vector<std::complex<double>> roots(std::vector<double> c) {
  while (!c.empty() && c.back() == 0.0) c.pop_back();
  const int n = static_cast<int>(c.size()) - 1;
  std::vector<std::complex<double>> r(z(std::max(n, 0)));
  if (n < 1) return r;
  const std::complex<double> one(1, 0);
  double rad = 0;
  for (int i = 0; i < n; ++i) rad = std::max(rad, std::pow(std::fabs(c[z(i)] / c[z(n)]), 1.0 / (n - i)));
  rad = 2 * rad + 1e-3;
  for (int i = 0; i < n; ++i) r[z(i)] = std::polar(rad, 2 * 3.141592653589793 * (i + 0.25) / n);
  auto ev = [&](std::complex<double> x, std::complex<double> &d) {
    std::complex<double> p(c[z(n)], 0);
    d = 0;
    for (int i = n - 1; i >= 0; --i) {
      d = d * x + p;
      p = p * x + c[z(i)];
    }
    return p;
  };
  for (int it = 0; it < 500; ++it) {
    double mx = 0;
    for (int i = 0; i < n; ++i) {
      std::complex<double> d;
      const std::complex<double> p = ev(r[z(i)], d);
      const std::complex<double> ratio = p / d;
      std::complex<double> s(0, 0);
      for (int j = 0; j < n; ++j) {
        if (j != i) s += one / (r[z(i)] - r[z(j)]);
      }
      const std::complex<double> w = ratio / (one - ratio * s);
      r[z(i)] -= w;
      mx = std::max(mx, std::abs(w) / std::max(1.0, std::abs(r[z(i)])));
    }
    if (mx < 1e-15) break;
  }
  return r;
}

bool archeck(const std::vector<double> &ar) {
  int p = 0;
  for (int i = 0; i < static_cast<int>(ar.size()); ++i) {
    if (ar[z(i)] != 0) p = i + 1;
  }
  if (!p) return true;
  std::vector<double> c(z(p + 1));
  c[0] = 1;
  for (int i = 0; i < p; ++i) c[z(i + 1)] = -ar[z(i)];
  for (const auto &r : roots(c)) {
    if (!(std::abs(r) > 1)) return false;
  }
  return true;
}

std::vector<double> mainvert(const std::vector<double> &ma) {
  const int q = static_cast<int>(ma.size());
  int q0 = 0;
  for (int i = 0; i < q; ++i) {
    if (ma[z(i)] != 0) q0 = i + 1;
  }
  if (!q0) return ma;
  std::vector<double> c(z(q0 + 1));
  c[0] = 1;
  for (int i = 0; i < q0; ++i) c[z(i + 1)] = ma[z(i)];
  std::vector<std::complex<double>> rt = roots(c);
  bool any = false;
  for (const auto &r : rt) any = any || std::abs(r) < 1;
  if (!any) return ma;
  std::vector<double> o(z(q), 0.0);
  if (q0 == 1) {
    o[0] = 1 / ma[0];
    return o;
  }
  for (auto &r : rt) {
    if (std::abs(r) < 1) r = 1.0 / r;
  }
  std::vector<std::complex<double>> x(1, 1.0);
  for (const auto &r : rt) {
    std::vector<std::complex<double>> nx(x.size() + 1, 0.0);
    for (std::size_t i = 0; i < x.size(); ++i) {
      nx[i] += x[i];
      nx[i + 1] -= x[i] / r;
    }
    x = nx;
  }
  for (int i = 0; i < q0; ++i) o[z(i)] = x[z(i + 1)].real();
  return o;
}

struct Model {
  std::vector<double> phi, theta, delta, Zv, a, P, T, V, Pn;
  int r = 0, rd = 0;
};

// The stationary covariance of the ARMA part of the state: the solution of
// P = T P T' + R R' (the Q0 of Gardner, Harvey and Phillips 1980).
bool q0(const std::vector<double> &phi, const std::vector<double> &theta, int r, std::vector<double> &Q) {
  const int p = static_cast<int>(phi.size()), q = static_cast<int>(theta.size()), n = r * r;
  std::vector<double> R(z(r), 0.0);
  R[0] = 1;
  for (int i = 1; i < r; ++i) R[z(i)] = i - 1 < q ? theta[z(i - 1)] : 0.0;
  auto T = [&](int i, int j) { return (j == 0 && i < p ? phi[z(i)] : 0.0) + (j == i + 1 ? 1.0 : 0.0); };
  Mat A(n, n), b(n, 1), x;
  for (int i = 0; i < r; ++i) {
    for (int j = 0; j < r; ++j) {
      const int row = i + r * j;
      b(row, 0) = R[z(i)] * R[z(j)];
      for (int k = 0; k < r; ++k) {
        for (int l = 0; l < r; ++l) A(row, k + r * l) -= T(i, k) * T(j, l);
      }
      A(row, row) += 1;
    }
  }
  if (!solve(A, b, x)) return false;
  Q.assign(x.v.begin(), x.v.end());
  for (int i = 0; i < r; ++i) {
    for (int j = 0; j < i; ++j) Q[z(i + r * j)] = Q[z(j + r * i)] = 0.5 * (Q[z(i + r * j)] + Q[z(j + r * i)]);
  }
  return true;
}

bool makearima(const std::vector<double> &phi, std::vector<double> theta, const std::vector<double> &delta,
               double kappa, Model &m) {
  const int p = static_cast<int>(phi.size()), q = static_cast<int>(theta.size());
  const int r = std::max(p, q + 1), d = static_cast<int>(delta.size()), rd = r + d;
  m.r = r;
  m.rd = rd;
  m.phi = phi;
  if (q < r - 1) theta.resize(z(r - 1), 0.0);
  m.theta = theta;
  m.delta = delta;
  m.Zv.assign(z(rd), 0.0);
  m.Zv[0] = 1;
  for (int i = 0; i < d; ++i) m.Zv[z(r + i)] = delta[z(i)];
  m.T.assign(z(rd) * z(rd), 0.0);
  for (int i = 0; i < p; ++i) m.T[z(i)] = phi[z(i)];
  for (int i = 1; i < r; ++i) m.T[z(i - 1 + rd * i)] = 1;
  if (d > 0) {
    for (int j = 0; j < rd; ++j) m.T[z(r + rd * j)] = m.Zv[z(j)];
    for (int i = 1; i < d; ++i) m.T[z(r + i + rd * (r + i - 1))] = 1;
  }
  std::vector<double> R(z(rd), 0.0);
  R[0] = 1;
  for (int i = 0; i < static_cast<int>(theta.size()); ++i) R[z(i + 1)] = theta[z(i)];
  m.V.assign(z(rd) * z(rd), 0.0);
  for (int i = 0; i < rd; ++i) {
    for (int j = 0; j < rd; ++j) m.V[z(i + rd * j)] = R[z(i)] * R[z(j)];
  }
  m.a.assign(z(rd), 0.0);
  m.P.assign(z(rd) * z(rd), 0.0);
  m.Pn.assign(z(rd) * z(rd), 0.0);
  if (r > 1) {
    std::vector<double> Q;
    if (!q0(phi, theta, r, Q)) return false;
    for (int i = 0; i < r; ++i) {
      for (int j = 0; j < r; ++j) m.Pn[z(i + rd * j)] = Q[z(i + r * j)];
    }
  } else {
    m.Pn[0] = p > 0 ? 1 / (1 - phi[0] * phi[0]) : 1;
  }
  for (int i = 0; i < d; ++i) m.Pn[z(r + i + rd * (r + i))] = kappa;
  return true;
}

// The Kalman filter of the exact likelihood: ssq, sum of log gains, n.
void like(const std::vector<double> &y, Model &m, double *out) {
  const int n = static_cast<int>(y.size()), rd = m.rd, r = m.r, d = rd - r;
  const int p = static_cast<int>(m.phi.size()), q = static_cast<int>(m.theta.size());
  const double *phi = m.phi.data(), *theta = m.theta.data(), *delta = m.delta.data();
  double *a = m.a.data(), *P = m.P.data(), *Pnew = m.Pn.data();
  std::vector<double> anew(z(rd)), M(z(rd)), mm(z(rd) * z(rd));
  double sumlog = 0, ssq = 0;
  int nu = 0;
  for (int l = 0; l < n; ++l) {
    for (int i = 0; i < r; ++i) {
      double t = i < r - 1 ? a[i + 1] : 0.0;
      if (i < p) t += phi[i] * a[0];
      anew[z(i)] = t;
    }
    if (d > 0) {
      for (int i = r + 1; i < rd; ++i) anew[z(i)] = a[i - 1];
      double t = a[0];
      for (int i = 0; i < d; ++i) t += delta[i] * a[r + i];
      anew[z(r)] = t;
    }
    if (l > 0) {
      if (d == 0) {
        for (int i = 0; i < r; ++i) {
          const double vi = i == 0 ? 1.0 : i - 1 < q ? theta[i - 1] : 0.0;
          for (int j = 0; j < r; ++j) {
            double t = j == 0 ? vi : j - 1 < q ? vi * theta[j - 1] : 0.0;
            if (i < p && j < p) t += phi[i] * phi[j] * P[0];
            if (i < r - 1 && j < r - 1) t += P[i + 1 + r * (j + 1)];
            if (i < p && j < r - 1) t += phi[i] * P[j + 1];
            if (j < p && i < r - 1) t += phi[j] * P[i + 1];
            Pnew[i + r * j] = t;
          }
        }
      } else {
        for (int i = 0; i < r; ++i) {
          for (int j = 0; j < rd; ++j) {
            double t = 0;
            if (i < p) t += phi[i] * P[rd * j];
            if (i < r - 1) t += P[i + 1 + rd * j];
            mm[z(i + rd * j)] = t;
          }
        }
        for (int j = 0; j < rd; ++j) {
          double t = P[rd * j];
          for (int k = 0; k < d; ++k) t += delta[k] * P[r + k + rd * j];
          mm[z(r + rd * j)] = t;
        }
        for (int i = 1; i < d; ++i) {
          for (int j = 0; j < rd; ++j) mm[z(r + i + rd * j)] = P[r + i - 1 + rd * j];
        }
        for (int i = 0; i < r; ++i) {
          for (int j = 0; j < rd; ++j) {
            double t = 0;
            if (i < p) t += phi[i] * mm[z(j)];
            if (i < r - 1) t += mm[z(rd * (i + 1) + j)];
            Pnew[j + rd * i] = t;
          }
        }
        for (int j = 0; j < rd; ++j) {
          double t = mm[z(j)];
          for (int k = 0; k < d; ++k) t += delta[k] * mm[z(rd * (r + k) + j)];
          Pnew[j + rd * r] = t;
        }
        for (int i = 1; i < d; ++i) {
          for (int j = 0; j < rd; ++j) Pnew[j + rd * (r + i)] = mm[z(rd * (r + i - 1) + j)];
        }
        for (int i = 0; i <= q; ++i) {
          const double vi = i == 0 ? 1.0 : theta[i - 1];
          for (int j = 0; j <= q; ++j) Pnew[i + rd * j] += vi * (j == 0 ? 1.0 : theta[j - 1]);
        }
      }
    }
    if (!std::isnan(y[z(l)])) {
      double res = y[z(l)] - anew[0];
      for (int i = 0; i < d; ++i) res -= delta[i] * anew[z(r + i)];
      double gain = M[0] = Pnew[0];
      for (int i = 0; i < d; ++i) gain += delta[i] * Pnew[r + i];
      for (int i = 0; i < rd; ++i) {
        double t = Pnew[i];
        for (int k = 0; k < d; ++k) t += Pnew[i + rd * (r + k)] * delta[k];
        M[z(i)] = t;
      }
      for (int i = 0; i < d; ++i) gain += delta[i] * M[z(r + i)];
      if (gain < 1e4) {
        ++nu;
        ssq += res * res / gain;
        sumlog += std::log(gain);
      }
      for (int i = 0; i < rd; ++i) a[i] = anew[z(i)] + M[z(i)] * res / gain;
      for (int i = 0; i < rd; ++i) {
        for (int j = 0; j < rd; ++j) P[i + rd * j] = Pnew[i + rd * j] - M[z(i)] * M[z(j)] / gain;
      }
    } else {
      for (int i = 0; i < rd; ++i) a[i] = anew[z(i)];
      for (int i = 0; i < rd * rd; ++i) P[i] = Pnew[i];
    }
  }
  out[0] = ssq;
  out[1] = sumlog;
  out[2] = nu;
}

// False when interrupted (the outputs are then incomplete).
bool fore(int h, Model m, std::vector<double> &fc, std::vector<double> &se) {
  const int p = m.rd;
  std::vector<double> anew(z(p)), Pnew(z(p) * z(p)), mm(z(p) * z(p));
  fc.assign(z(h), 0.0);
  se.assign(z(h), 0.0);
  for (int l = 0; l < h; ++l) {
    if (kutil::stop(3.0 * p * p * p)) return false;
    double f = 0;
    for (int i = 0; i < p; ++i) {
      double t = 0;
      for (int k = 0; k < p; ++k) t += m.T[z(i + p * k)] * m.a[z(k)];
      anew[z(i)] = t;
      f += t * m.Zv[z(i)];
    }
    m.a = anew;
    fc[z(l)] = f;
    for (int i = 0; i < p; ++i) {
      for (int j = 0; j < p; ++j) {
        double t = 0;
        for (int k = 0; k < p; ++k) t += m.T[z(i + p * k)] * m.P[z(k + p * j)];
        mm[z(i + p * j)] = t;
      }
    }
    for (int i = 0; i < p; ++i) {
      for (int j = 0; j < p; ++j) {
        double t = m.V[z(i + p * j)];
        for (int k = 0; k < p; ++k) t += mm[z(i + p * k)] * m.T[z(j + p * k)];
        Pnew[z(i + p * j)] = t;
      }
    }
    double t = 0;
    for (int i = 0; i < p; ++i) {
      for (int j = 0; j < p; ++j) {
        m.P[z(i + j * p)] = Pnew[z(i + j * p)];
        t += m.Zv[z(i)] * m.Zv[z(j)] * m.P[z(i + j * p)];
      }
    }
    se[z(l)] = t;
  }
  return true;
}

// R's vmmin (BFGS) with the central differences of optim() (ndeps 1e-3,
// parscale 1), reltol sqrt(epsilon), maxit 100. False when a difference
// is not finite.
using Fn = std::function<double(const std::vector<double> &)>;

bool bfgs(const Fn &fn, std::vector<double> &b, double &fmin, int maxit, bool &conv, std::string &err) {
  const int n = static_cast<int>(b.size());
  const double reltol = 1.490116119384765625e-08, abstol = -HUGE_VAL;
  auto grad = [&](const std::vector<double> &p, std::vector<double> &g) {
    std::vector<double> x = p;
    g.assign(z(n), 0.0);
    for (int i = 0; i < n; ++i) {
      x[z(i)] = p[z(i)] + 1e-3;
      const double v1 = fn(x);
      x[z(i)] = p[z(i)] - 1e-3;
      const double v2 = fn(x);
      g[z(i)] = (v1 - v2) / (2 * 1e-3);
      if (!std::isfinite(g[z(i)])) {
        err = msg("non-finite finite-difference value [%d]", i + 1);
        return false;
      }
      x[z(i)] = p[z(i)];
    }
    return true;
  };
  conv = true;
  if (maxit <= 0) {
    fmin = fn(b);
    return true;
  }
  double f = fn(b);
  if (!std::isfinite(f)) {
    err = "initial value in 'vmmin' is not finite";
    return false;
  }
  fmin = f;
  std::vector<double> g, t(z(n)), X(z(n)), c(z(n));
  std::vector<std::vector<double>> B(z(n), std::vector<double>(z(n)));
  if (!grad(b, g)) return false;
  int gradcount = 1, iter = 1, ilast = gradcount, count = 0;
  do {
    if (ilast == gradcount) {
      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) B[z(i)][z(j)] = 0.0;
        B[z(i)][z(i)] = 1.0;
      }
    }
    for (int i = 0; i < n; ++i) {
      X[z(i)] = b[z(i)];
      c[z(i)] = g[z(i)];
    }
    double gradproj = 0;
    for (int i = 0; i < n; ++i) {
      double s = 0;
      for (int j = 0; j <= i; ++j) s -= B[z(i)][z(j)] * g[z(j)];
      for (int j = i + 1; j < n; ++j) s -= B[z(j)][z(i)] * g[z(j)];
      t[z(i)] = s;
      gradproj += s * g[z(i)];
    }
    if (gradproj < 0.0) {
      double step = 1.0;
      bool acc = false;
      do {
        count = 0;
        for (int i = 0; i < n; ++i) {
          b[z(i)] = X[z(i)] + step * t[z(i)];
          if (10.0 + X[z(i)] == 10.0 + b[z(i)]) ++count;
        }
        if (count < n) {
          f = fn(b);
          acc = std::isfinite(f) && f <= fmin + gradproj * step * 0.0001;
          if (!acc) step *= 0.2;
        }
      } while (!(count == n || acc));
      const bool enough = f > abstol && std::fabs(f - fmin) > reltol * (std::fabs(fmin) + reltol);
      if (!enough) {
        count = n;
        fmin = f;
      }
      if (count < n) {
        fmin = f;
        if (!grad(b, g)) return false;
        ++gradcount;
        ++iter;
        double D1 = 0;
        for (int i = 0; i < n; ++i) {
          t[z(i)] = step * t[z(i)];
          c[z(i)] = g[z(i)] - c[z(i)];
          D1 += t[z(i)] * c[z(i)];
        }
        if (D1 > 0) {
          double D2 = 0;
          for (int i = 0; i < n; ++i) {
            double s = 0;
            for (int j = 0; j <= i; ++j) s += B[z(i)][z(j)] * c[z(j)];
            for (int j = i + 1; j < n; ++j) s += B[z(j)][z(i)] * c[z(j)];
            X[z(i)] = s;
            D2 += s * c[z(i)];
          }
          D2 = 1.0 + D2 / D1;
          for (int i = 0; i < n; ++i) {
            for (int j = 0; j <= i; ++j) {
              B[z(i)][z(j)] += (D2 * t[z(i)] * t[z(j)] - X[z(i)] * t[z(j)] - t[z(i)] * X[z(j)]) / D1;
            }
          }
        } else {
          ilast = gradcount;
        }
      } else {
        if (ilast < gradcount) {
          count = 0;
          ilast = gradcount;
        }
      }
    } else {
      count = 0;
      if (ilast == gradcount) {
        count = n;
      } else {
        ilast = gradcount;
      }
    }
    if (iter >= maxit) break;
    if (gradcount - ilast > 2 * n) ilast = gradcount;
  } while (count != n || ilast != gradcount);
  conv = iter < maxit;
  return true;
}

struct Fit {
  std::vector<double> coef;
  double sigma2 = 0, loglik = 0, aic = 0;
  int nused = 0;
  Model mod;
};

// stats::arima(x, order, seasonal) by its default "CSS-ML", without a mean
// (the series is differenced at least once).
bool arima(const std::vector<double> &x, const Arma &a, Fit &ft, std::string &err) {
  std::vector<double> Dl(1, 1.0);
  auto conv = [](const std::vector<double> &u, const std::vector<double> &v) {
    std::vector<double> o(u.size() + v.size() - 1, 0.0);
    for (std::size_t i = 0; i < u.size(); ++i) {
      for (std::size_t j = 0; j < v.size(); ++j) o[i + j] += u[i] * v[j];
    }
    return o;
  };
  for (int i = 0; i < a.d; ++i) Dl = conv(Dl, {1.0, -1.0});
  for (int i = 0; i < a.sd; ++i) {
    std::vector<double> s(z(a.ns + 1), 0.0);
    s[0] = 1;
    s[z(a.ns)] = -1;
    Dl = conv(Dl, s);
  }
  std::vector<double> delta(Dl.begin() + 1, Dl.end());
  for (double &e : delta) e = -e;
  int nobs = 0;
  for (double e : x) nobs += !std::isnan(e);
  const int nused = nobs - static_cast<int>(delta.size());
  const int narma = a.np();
  bool anyna = nobs < static_cast<int>(x.size());
  int ncond = a.d + a.sd * a.ns + a.mp + a.ns * a.msp;
  std::vector<double> init(z(narma), 0.0);
  if (!anyna) {
    auto armacss = [&](const std::vector<double> &p) {
      std::vector<double> phi, theta;
      transpars(p, a, false, phi, theta);
      return 0.5 * std::log(css(x, a, phi, theta, ncond));
    };
    std::vector<double> b = init;
    double fm;
    bool cv;
    if (narma > 0) {
      if (!bfgs(armacss, b, fm, 100, cv, err)) return false;
      if (cv) init = b;
    }
    if (a.mp > 0 && !archeck(std::vector<double>(init.begin(), init.begin() + a.mp))) {
      err = "non-stationary AR part from CSS";
      return false;
    }
    if (a.msp > 0 &&
        !archeck(std::vector<double>(init.begin() + a.mp + a.mq, init.begin() + a.mp + a.mq + a.msp))) {
      err = "non-stationary seasonal AR part from CSS";
      return false;
    }
  }
  init = invtrans(init, a);
  auto mai = [&](std::vector<double> &v) {
    if (a.mq > 0) {
      const std::vector<double> m = mainvert(std::vector<double>(v.begin() + a.mp, v.begin() + a.mp + a.mq));
      std::copy(m.begin(), m.end(), v.begin() + a.mp);
    }
    if (a.msq > 0) {
      const int o = a.mp + a.mq + a.msp;
      const std::vector<double> m = mainvert(std::vector<double>(v.begin() + o, v.begin() + o + a.msq));
      std::copy(m.begin(), m.end(), v.begin() + o);
    }
  };
  mai(init);
  auto armafn = [&](const std::vector<double> &p) {
    std::vector<double> phi, theta;
    transpars(p, a, true, phi, theta);
    Model m;
    // upARIMA: theta is not padded, the filter reads its length only
    if (!makearima(phi, theta, delta, 1e6, m)) return std::numeric_limits<double>::max();
    m.theta = theta;
    double r[3];
    like(x, m, r);
    const double s2 = r[0] / r[2];
    return 0.5 * (std::log(s2) + r[1] / r[2]);
  };
  std::vector<double> coef = init;
  double val = 0;
  bool cv = true;
  if (narma > 0) {
    if (!bfgs(armafn, coef, val, 100, cv, err)) return false;
    std::vector<double> c2 = coef;
    mai(c2);
    if (c2 != coef) {
      coef = c2;
      val = armafn(coef);
    }
  } else {
    val = armafn(coef);
  }
  coef = undopars(coef, a);
  std::vector<double> phi, theta;
  transpars(coef, a, false, phi, theta);
  // an AR(1) has no system for q0() to find singular: its state variance
  // 1 / (1 - phi^2) is refused where 1 - phi^2 is below the 1e-14 that
  // solve() asks of a pivot (so at and past |phi| = 1 too)
  const bool unit = phi.size() == 1 && theta.empty() && !(1 - phi[0] * phi[0] > 1e-14);
  if (unit || !makearima(phi, theta, delta, 1e6, ft.mod)) {
    err = "the fitted ARMA part is not stationary.";
    return false;
  }
  double r[3];
  like(x, ft.mod, r);
  ft.sigma2 = r[0] / nused;
  ft.coef = coef;
  ft.nused = nused;
  const double v2 = 2 * nused * val + nused + nused * std::log(2 * 3.141592653589793238462643383280);
  ft.loglik = -0.5 * v2;
  ft.aic = v2 + 2 * narma + 2;
  return true;
}

bool getord(const Node &a, const char *nm, int n, std::vector<int> &o, std::string &err) {
  std::vector<double> v;
  if (!getv(a, nm, v, err)) return false;
  if (static_cast<int>(v.size()) != n) {
    err = msg("`%s` must have %d whole numbers from 0.", nm, n);
    return false;
  }
  o.clear();
  for (double e : v) {
    if (e < 0 || e != std::floor(e) || e > 1000) {
      err = msg("`%s` must have %d whole numbers from 0.", nm, n);
      return false;
    }
    o.push_back(static_cast<int>(e));
  }
  return true;
}

}  // namespace

bool k_arimfc(const Node &a, Node &out, std::string &err) {
  std::vector<double> y;
  int h;
  if (!getv(a, "y", y, err) || !geti(a, "h", h, err, 1) || !cap(h, kutil::maxcell, "`h`", err)) return false;
  if (y.size() < 24) {
    err = msg("`y` must have 24 values at least, got %d.", static_cast<int>(y.size()));
    return false;
  }
  Fit f;
  if (!arima(y, {1, 1, 0, 0, 12, 1, 0}, f, err)) return false;
  std::vector<double> fc, se;
  if (!fore(h, f.mod, fc, se)) {
    err.clear();
    return false;
  }
  for (double &e : se) e = std::sqrt(e * f.sigma2);
  const double k = static_cast<double>(f.coef.size()) + 1;
  out.add("fc", vec(fc));
  out.add("se", vec(se));
  out.add("coef", vec(f.coef));
  out.add("sigma2", num(f.sigma2));
  out.add("aic", num(f.aic));
  out.add("bic", num(f.aic + k * (std::log(static_cast<double>(y.size())) - 2)));
  out.add("n", ival(static_cast<double>(y.size())));
  return true;
}

bool k_sarfc(const Node &a, Node &out, std::string &err) {
  std::vector<double> y;
  std::vector<int> o, s;
  int h;
  if (!getv(a, "y", y, err) || !geti(a, "h", h, err, 1) || !getord(a, "order", 3, o, err) ||
      !getord(a, "seasonal", 4, s, err)) {
    return false;
  }
  if (s[3] < 1) {
    err = "the period, the fourth value of `seasonal`, must be at least 1.";
    return false;
  }
  if (o[1] + s[1] < 1) {
    err = "the model must difference the series at least once: `order[2]` or `seasonal[2]` above 0.";
    return false;
  }
  const int n = static_cast<int>(y.size()), tn = n - h;
  if (n < 36) {
    err = msg("`y` must have 36 values at least, got %d.", n);
    return false;
  }
  if (tn < 24) {
    err = msg("only %d training values are left before the last `h`; 24 are needed.", tn);
    return false;
  }
  std::vector<double> tr(y.begin(), y.begin() + tn), te(y.begin() + tn, y.end());
  Fit f;
  if (!arima(tr, {o[0], o[2], s[0], s[2], s[3], o[1], s[1]}, f, err)) return false;
  std::vector<double> fc, se;
  if (!fore(h, f.mod, fc, se)) {
    err.clear();
    return false;
  }
  double ap = 0, sq = 0;
  for (int i = 0; i < h; ++i) {
    const double e = te[z(i)] - fc[z(i)];
    ap += std::fabs(e / std::max(te[z(i)], 1.0));
    sq += e * e;
  }
  const double k = static_cast<double>(f.coef.size()) + 1;
  out.add("fc", vec(fc));
  out.add("actual", vec(te));
  out.add("coef", vec(f.coef));
  out.add("aic", num(f.aic));
  out.add("bic", num(f.aic + k * (std::log(static_cast<double>(tn)) - 2)));
  out.add("mape", num(ap / h * 100));
  out.add("rmse", num(std::sqrt(sq / h)));
  out.add("ntrain", ival(tn));
  return true;
}

// --- correlation of counts per unit ------------------------------------------

bool k_catcor(const Node &a, Node &out, std::string &err) {
  const Node *u = a.at("units");
  if (!u || u->k != Node::LIST || u->l.size() < 2 || u->nm.size() != u->l.size()) {
    err = "`units` must be a named list of two categories at least, each the unit of every event.";
    return false;
  }
  std::vector<std::size_t> ord(u->l.size());
  for (std::size_t i = 0; i < ord.size(); ++i) ord[i] = i;
  std::stable_sort(ord.begin(), ord.end(), [&](std::size_t p, std::size_t q) { return u->nm[p] < u->nm[q]; });
  std::map<std::string, std::vector<double>> tab;
  const int K = static_cast<int>(ord.size());
  std::vector<std::string> cats;
  for (int c = 0; c < K; ++c) {
    const Node &e = u->l[ord[z(c)]];
    cats.push_back(u->nm[ord[z(c)]]);
    if (e.k != Node::STR) {
      err = "`units` must hold character vectors.";
      return false;
    }
    for (const std::string &s : e.s) {
      std::vector<double> &r = tab[s];
      r.resize(z(K), 0.0);
      r[z(c)] += 1;
    }
  }
  const int n = static_cast<int>(tab.size());
  if (n < 2) {
    err = "`units` must name two units at least.";
    return false;
  }
  Mat P(n, K);
  Node un;
  un.k = Node::STR;
  int i = 0;
  for (auto &kv : tab) {
    kv.second.resize(z(K), 0.0);
    for (int c = 0; c < K; ++c) P(i, c) = kv.second[z(c)];
    un.s.push_back(kv.first);
    ++i;
  }
  Mat R(K, K);
  for (int p = 0; p < K; ++p) {
    for (int q = 0; q < K; ++q) {
      double mp = 0, mq = 0;
      for (int t = 0; t < n; ++t) {
        mp += P(t, p);
        mq += P(t, q);
      }
      mp /= n;
      mq /= n;
      double sxy = 0, sxx = 0, syy = 0;
      for (int t = 0; t < n; ++t) {
        sxy += (P(t, p) - mp) * (P(t, q) - mq);
        sxx += (P(t, p) - mp) * (P(t, p) - mp);
        syy += (P(t, q) - mq) * (P(t, q) - mq);
      }
      R(p, q) = sxx > 0 && syy > 0 ? std::min(1.0, std::max(-1.0, sxy / std::sqrt(sxx * syy))) : NaN;
    }
  }
  Node cn;
  cn.k = Node::STR;
  cn.s = cats;
  out.add("r", mnode(R));
  out.add("counts", mnode(P));
  out.add("units", un);
  out.add("cats", cn);
  return true;
}

// --- Ornstein-Uhlenbeck fit, simulation and density ---------------------------

namespace {

struct Ou {
  double theta, mu, sigma;
};

// Least squares of the first differences on the lagged values:
// dx = a + b x, theta = -b (at least 1e-6), mu = -a / b, sigma the
// standard deviation of the residuals.
Ou oufit(const std::vector<double> &x) {
  const int n = static_cast<int>(x.size()) - 1;
  std::vector<double> dx(z(n)), xl(z(n));
  for (int i = 0; i < n; ++i) {
    dx[z(i)] = x[z(i + 1)] - x[z(i)];
    xl[z(i)] = x[z(i)];
  }
  const double mx = mean(xl), my = mean(dx);
  double sxy = 0, sxx = 0;
  for (int i = 0; i < n; ++i) {
    sxy += (xl[z(i)] - mx) * (dx[z(i)] - my);
    sxx += (xl[z(i)] - mx) * (xl[z(i)] - mx);
  }
  const double b = sxx > 0 ? sxy / sxx : NaN, a0 = my - (std::isnan(b) ? 0.0 : b) * mx;
  Ou o;
  o.theta = std::max(1e-6, std::isnan(b) ? 1e-6 : -b);
  o.mu = !std::isnan(b) && b != 0 ? -a0 / b : mean(x);
  std::vector<double> r(z(n));
  for (int i = 0; i < n; ++i) r[z(i)] = dx[z(i)] - (a0 + (std::isnan(b) ? 0.0 : b) * xl[z(i)]);
  o.sigma = sdv(r.data(), n);
  if (!std::isfinite(o.sigma)) o.sigma = 0.0;
  return o;
}

bool series(const Node &a, std::vector<double> &y, std::string &err) {
  if (!getv(a, "y", y, err)) return false;
  if (y.size() < 3) {
    err = "`y` must have three values at least.";
    return false;
  }
  return true;
}

}  // namespace

bool k_ousim(const Node &a, Node &out, std::string &err) {
  std::vector<double> y;
  int np, seed;
  double days, dt;
  if (!series(a, y, err) || !geti(a, "npath", np, err, 1) || !getd(a, "days", days, err, 0, 1e7) ||
      !getd(a, "dt", dt, err, 1e-9, HUGE_VAL) || !geti(a, "seed", seed, err, -2147483647)) {
    return false;
  }
  const Ou f = oufit(y);
  const int ns = std::max(1, static_cast<int>(days / dt));
  if (static_cast<double>(ns) * np > 1e8) {
    err = "`npath` times the steps must be at most 1e8.";
    return false;
  }
  rng g(seed);
  Mat P(np, ns);
  for (int i = 0; i < np; ++i) P(i, 0) = y.back();
  const double sq = std::sqrt(dt);
  for (int k = 1; k < ns; ++k) {
    if (kutil::stop(np)) {
      err.clear();
      return false;
    }
    for (int i = 0; i < np; ++i) P(i, k) = P(i, k - 1) + f.theta * (f.mu - P(i, k - 1)) * dt + f.sigma * sq * g.norm();
  }
  std::vector<double> fin(P.v.end() - np, P.v.end());
  out.add("theta", num(f.theta));
  out.add("mu", num(f.mu));
  out.add("sigma", num(f.sigma));
  out.add("halflife", num(std::log(2.0) / std::max(f.theta, 1e-6)));
  out.add("q", vec({quant(fin, 0.05), quant(fin, 0.5), quant(fin, 0.95)}));
  out.add("paths", mnode(P));
  out.add("n", ival(static_cast<double>(y.size())));
  return true;
}

bool k_fpdens(const Node &a, Node &out, std::string &err) {
  std::vector<double> y;
  int ng, ns;
  if (!series(a, y, err) || !geti(a, "ngrid", ng, err, 3) || !geti(a, "nstep", ns, err, 0) ||
      !cap(ng, kutil::maxcell, "`ngrid`", err)) {
    return false;
  }
  const Ou f = oufit(y);
  double mx = y[0];
  for (double e : y) mx = std::max(mx, e);
  const double xmax = mx * 1.5 + 1;
  std::vector<double> gr(z(ng));
  for (int i = 0; i < ng; ++i) gr[z(i)] = 0 + i * ((xmax - 0) / (ng - 1));
  gr[z(ng - 1)] = xmax;
  const double dx = gr[1] - gr[0], dts = 0.05;
  std::vector<double> p(z(ng)), dr(z(ng));
  const double last = y.back();
  for (int i = 0; i < ng; ++i) {
    p[z(i)] = std::exp(-((gr[z(i)] - last) * (gr[z(i)] - last)) / (2 * (f.sigma + 1) * (f.sigma + 1)));
    dr[z(i)] = f.theta * (f.mu - gr[z(i)]);
  }
  double s0 = 0;
  for (double e : p) s0 += e;
  s0 *= dx;
  if (s0 > 0) {
    for (double &e : p) e /= s0;
  }
  const double dc = 0.5 * f.sigma * f.sigma;
  for (int st = 0; st < ns; ++st) {
    if (kutil::stop(10.0 * ng)) {
      err.clear();
      return false;
    }
    std::vector<double> np(z(ng));
    for (int i = 0; i < ng; ++i) {
      const double l = p[z((i + 1) % ng)], r = p[z((i + ng - 1) % ng)];
      const double d2 = (l - 2 * p[z(i)] + r) / (dx * dx), d1 = (l - r) / (2 * dx);
      np[z(i)] = p[z(i)] + dts * (-dr[z(i)] * d1 + dc * d2);
    }
    np[0] = np[1];
    np[z(ng - 1)] = np[z(ng - 2)];
    double s = 0;
    for (double &e : np) {
      e = std::max(e, 0.0);
      s += e;
    }
    s *= dx;
    if (s > 0) {
      for (double &e : np) e /= s;
    }
    p = np;
  }
  const int pk = static_cast<int>(std::max_element(p.begin(), p.end()) - p.begin());
  out.add("theta", num(f.theta));
  out.add("mu", num(f.mu));
  out.add("sigma", num(f.sigma));
  out.add("grid", vec(gr));
  out.add("density", vec(p));
  out.add("peak", num(gr[z(pk)]));
  out.add("statvar", num(f.theta > 0 ? f.sigma * f.sigma / (2 * f.theta) : NaN));
  return true;
}

// --- Lotka-Volterra rates -------------------------------------------------------

bool k_lvfit(const Node &a, Node &out, std::string &err) {
  std::vector<double> c, yr;
  if (!getv(a, "y", c, err) || !getv(a, "year", yr, err, false, true)) return false;
  const int n = static_cast<int>(c.size());
  if (n < 5) {
    err = msg("`y` must have five years at least, got %d.", n);
    return false;
  }
  for (double e : c) {
    if (e < 0) {
      err = "`y` must hold counts from 0.";
      return false;
    }
  }
  if (yr.empty()) {
    for (int i = 0; i < n; ++i) yr.push_back(i + 1);
  } else if (static_cast<int>(yr.size()) != n) {
    err = msg("`year` must have one value per count, %d.", n);
    return false;
  }
  std::vector<double> x(z(n)), y(z(n));
  for (int i = 0; i < n; ++i) x[z(i)] = std::max(c[z(i)], 1.0);
  const double f3 = 1.0 / 3;
  for (int i = 0; i < n; ++i) {
    if (i < 2) {
      y[z(i)] = x[z(i)];
      continue;
    }
    double s = 0;
    for (int j = 0; j < 3; ++j) s += f3 * x[z(i - j)];
    y[z(i)] = s;
  }
  std::vector<double> lx(z(n - 1)), ly(z(n - 1));
  for (int i = 0; i + 1 < n; ++i) {
    lx[z(i)] = std::log(x[z(i + 1)]) - std::log(x[z(i)]);
    ly[z(i)] = std::log(std::max(y[z(i + 1)], 1.0)) - std::log(std::max(y[z(i)], 1.0));
  }
  const double al = median(lx) + 0.05, ga = -median(ly) + 0.05;
  const double T = 2 * 3.141592653589793238462643383280 / std::sqrt(std::max(al * ga, 1e-6));
  out.add("alpha", num(al));
  out.add("beta", num(al / std::max(median(y), 1.0)));
  out.add("gamma", num(ga));
  out.add("delta", num(ga / std::max(median(x), 1.0)));
  out.add("period", num(T));
  out.add("year", vec(yr));
  out.add("prey", vec(x));
  out.add("pred", vec(y));
  return true;
}

}  // namespace svk
