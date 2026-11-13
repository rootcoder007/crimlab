// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, hot-spot dynamics: the Short-D'Orsogna-Brantingham model of
// residential burglary, its agent version on a periodic lattice and its
// reaction-diffusion limit (clab_sdbsim), and the linear stability of the
// uniform state of the latter (clab_sdbthr). No allocation: the caller
// passes the work arrays.
#include <cmath>
#include <cstddef>

#include "crimlab.h"
#include "kutil.h"

namespace {

// the four neighbours of site i on the n x n torus (column-major)
void nbrs(int i, int n, int *nb) {
  const int r = i % n, c = i / n;
  nb[0] = (r + n - 1) % n + c * n;
  nb[1] = (r + 1) % n + c * n;
  nb[2] = r + ((c + n - 1) % n) * n;
  nb[3] = r + ((c + 1) % n) * n;
}

int pois(kutil::rng &g, double mu) {
  int k = 0;
  double p = std::exp(-mu), c = p;
  for (const double u = g.unif(); u > c && k < 100000; c += p) p *= mu / ++k;
  return k;
}

}  // namespace

extern "C" {

int clab_sdbsim(int n, int steps, int mode, const double *par, double *b, double *rho,
                const int *sup, unsigned int seed, double *burg, double *meanb, double *tot,
                double *work, int *iwork) {
  if (n < 3 || steps < 0 || mode < 0 || mode > 1) return CLAB_EARG;
  const double a0 = par[0], om = par[1], eta = par[2], gam = par[3], th = par[4], dt = par[5],
               dx = par[6], dif = par[7], eps = par[8], z = par[9];
  if (!kutil::pos(a0) || !kutil::nonneg(om) || !(eta >= 0 && eta <= 1) || !kutil::nonneg(gam) ||
      !kutil::nonneg(th) || !kutil::pos(dt) || !kutil::pos(dx) || !kutil::pos(dif) ||
      !kutil::nonneg(eps) || !kutil::pos(z) || om * dt > 1) {
    return CLAB_EDOM;
  }
  const int nn = n * n;
  for (int i = 0; i < nn; ++i) {
    if (!kutil::nonneg(b[i]) || !kutil::nonneg(rho[i])) return CLAB_EDOM;
    tot[i] = 0;
  }
  double *nb2 = work, *nr = work + nn, *a = work + 2 * static_cast<std::size_t>(nn);
  int nb[4];
  if (mode == 0) {
    // agents: rho holds whole numbers of offenders per site
    kutil::rng g(seed);
    int *e = iwork;
    for (int t = 0; t < steps; ++t) {
      if (kutil::stop(8.0 * nn)) return CLAB_EINTR;
      for (int i = 0; i < nn; ++i) {
        a[i] = a0 + b[i];
        nr[i] = 0;
        e[i] = 0;
      }
      double sb = 0;
      for (int i = 0; i < nn; ++i) {
        const int c = static_cast<int>(rho[i]);
        const double pb = -std::expm1(-a[i] * dt);
        nbrs(i, n, nb);
        const double w = a[nb[0]] + a[nb[1]] + a[nb[2]] + a[nb[3]];
        for (int k = 0; k < c; ++k) {
          if (g.unif() < pb) {
            ++e[i];
            continue;
          }
          double u = g.unif() * w;
          int j = 0;
          while (j < 3 && u >= a[nb[j]]) u -= a[nb[j++]];
          nr[nb[j]] += 1;
        }
      }
      int ne = 0;
      for (int i = 0; i < nn; ++i) {
        nbrs(i, n, nb);
        nb2[i] = ((1 - eta) * b[i] + eta / 4 * (b[nb[0]] + b[nb[1]] + b[nb[2]] + b[nb[3]])) * (1 - om * dt) +
                 th * e[i];
        ne += e[i];
        tot[i] += e[i];
      }
      for (int i = 0; i < nn; ++i) {
        b[i] = nb2[i];
        sb += b[i];
        rho[i] = nr[i] + pois(g, gam * dt);
      }
      burg[t] = ne;
      meanb[t] = sb / nn;
    }
    return CLAB_OK;
  }
  // the reaction-diffusion limit by explicit Euler steps; the step must
  // keep the diffusion stable
  if (dt * 4 * dif * std::fmax(1.0, eta) / z > dx * dx) return CLAB_EDOM;
  const double h2 = dx * dx;
  for (int t = 0; t < steps; ++t) {
    if (kutil::stop(8.0 * nn)) return CLAB_EINTR;
    double sb = 0, sc = 0;
    for (int i = 0; i < nn; ++i) a[i] = a0 + b[i];
    for (int i = 0; i < nn; ++i) {
      nbrs(i, n, nb);
      double lb = -4 * b[i], lr = -4 * rho[i], adv = 0;
      const double q = rho[i] / a[i];
      for (int k = 0; k < 4; ++k) {
        lb += b[nb[k]];
        lr += rho[nb[k]];
        adv += (q + rho[nb[k]] / a[nb[k]]) / 2 * (a[nb[k]] - a[i]);
      }
      const double cr = sup && sup[i] ? 0.0 : rho[i] * a[i];
      nb2[i] = b[i] + dt * (eta * dif / z * lb / h2 - om * b[i] + eps * dif * cr);
      nr[i] = rho[i] + dt * (dif / z * (lr - 2 * adv) / h2 - cr + gam);
      sc += cr;
      tot[i] += cr * dt;
    }
    for (int i = 0; i < nn; ++i) {
      if (!(nb2[i] >= 0 && nr[i] >= 0)) return CLAB_EDOM;
      b[i] = nb2[i];
      rho[i] = nr[i];
      sb += b[i];
    }
    burg[t] = sc * dt;
    meanb[t] = sb / nn;
  }
  return CLAB_OK;
}

int clab_sdbthr(const double *par, const double *k, int nk, double *st, double *tr, double *det,
                double *gr) {
  const double a0 = par[0], om = par[1], eps = par[2], dif = par[3], gam = par[4], eta = par[5],
               z = par[6];
  if (!kutil::pos(a0) || !kutil::pos(om) || !kutil::nonneg(eps) || !kutil::pos(dif) ||
      !kutil::pos(gam) || !kutil::nonneg(eta) || !kutil::pos(z) || nk < 0) {
    return CLAB_EDOM;
  }
  // the uniform steady state and the linearisation at wavenumber-squared k:
  // [[a - d1 k, bb], [c0 + c1 k, d - d2 k]]
  const double bb0 = eps * dif * gam / om, ab = a0 + bb0, rb = gam / ab;
  const double a = -om + eps * dif * rb, bb = eps * dif * ab, c0 = -rb, c1 = 2 * dif * rb / (z * ab),
               d = -ab, d1 = eta * dif / z, d2 = dif / z;
  const double p = d1 * d2, q = a * d2 + d * d1 + bb * c1, r = a * d - bb * c0;
  auto trk = [&](double x) { return a - d1 * x + d - d2 * x; };
  auto dtk = [&](double x) { return (a - d1 * x) * (d - d2 * x) - bb * (c0 + c1 * x); };
  auto grk = [&](double x) {
    const double t = trk(x), dd = dtk(x);
    return (t + std::sqrt(std::fmax(t * t - 4 * dd, 0.0))) / 2;
  };
  st[0] = p;
  st[1] = q;
  st[2] = r;
  st[3] = bb0;
  st[4] = ab;
  st[5] = rb;
  st[6] = 3 * eps * dif * rb;
  st[7] = om + eta * ab;
  st[8] = p > 0 ? (q > 0 && q * q > 4 * p * r) : q > 0;
  st[9] = p > 0 && q > 0 ? q / (2 * p) : kutil::nan();
  st[10] = std::isnan(st[9]) ? kutil::nan() : dtk(st[9]);
  st[11] = std::isnan(st[9]) ? kutil::nan() : grk(st[9]);
  for (int i = 0; i < nk; ++i) {
    if (!kutil::nonneg(k[i])) return CLAB_EDOM;
    tr[i] = trk(k[i]);
    det[i] = dtk(k[i]);
    gr[i] = grk(k[i]);
  }
  return CLAB_OK;
}

}  // extern "C"
