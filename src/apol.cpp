// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Air pollution: Pasquill-Gifford dispersion coefficients (Briggs 1973),
// Briggs plume rise (Briggs 1975; EPA ISC3), the Gaussian plume and puff
// with ground and lid reflection (Turner 1970; Seinfeld and Pandis 2016),
// an explicit Eulerian advection-diffusion solver, a Lagrangian random walk
// (Thomson 1987), and the health side: concentration-response functions
// (Chen and Hoek 2020; Huangfu and Atkinson 2020; Khreis et al. 2017;
// Burnett et al. 2014), Levin's attributable fraction, the BenMAP
// health-impact function and the concentration index of Wagstaff, Paci and
// van Doorslaer (1991).
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

namespace {

using std::size_t;
const double pi = 3.14159265358979323846;

// The Briggs (1973) coefficients by class A-F: a_y, a_z, b_z, e_z.
const double rural[6][4] = {{0.22, 0.20, 0, 0},           {0.16, 0.12, 0, 0},
                            {0.11, 0.08, 0.0002, -0.5},   {0.08, 0.06, 0.0015, -0.5},
                            {0.06, 0.03, 0.0003, -1},     {0.04, 0.016, 0.0003, -1}};
const double urban[6][4] = {{0.32, 0.24, 0.001, 0.5},     {0.32, 0.24, 0.001, 0.5},
                            {0.22, 0.20, 0, 0},           {0.16, 0.14, 0.0003, -0.5},
                            {0.11, 0.08, 0.0015, -0.5},   {0.11, 0.08, 0.0015, -0.5}};

void sig(double x, int cls, int urb, double &sy, double &sz) {
  const double *c = urb ? urban[cls] : rural[cls];
  const double cy = urb ? 0.0004 : 0.0001;
  sy = c[0] * x * std::pow(1 + cy * x, -0.5);
  sz = c[1] * x * std::pow(1 + c[2] * x, c[3]);
}

double sq(double v) { return v * v; }

// The vertical term: ground reflection, and nimg image pairs under a lid.
double vert(double z, double h, double sz, double lid, int nimg) {
  double t = std::exp(-sq(z - h) / (2 * sz * sz)) + std::exp(-sq(z + h) / (2 * sz * sz));
  if (!std::isnan(lid)) {
    for (int j = 1; j <= nimg; ++j) {
      for (int s = 1; s >= -1; s -= 2) {
        const double off = 2 * s * j * lid;
        t = t + std::exp(-sq(z - h + off) / (2 * sz * sz)) + std::exp(-sq(z + h + off) / (2 * sz * sz));
      }
    }
  }
  return t;
}

// R's x %% y for doubles, in extended precision.
double rmod(double x, double y) {
  const double q = x / y;
  const long double t = static_cast<long double>(x) - std::floor(q) * static_cast<long double>(y);
  return static_cast<double>(t - std::floor(t / y) * y);
}

// The uniform stream of Steele, Lea and Flood (2014), splitmix64.
void smix(std::uint64_t seed, size_t n, double *out) {
  std::uint64_t x = seed;
  for (size_t i = 0; i < n; ++i) {
    x += 0x9E3779B97F4A7C15ULL;
    std::uint64_t z = x;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    z ^= z >> 31;
    out[i] = static_cast<double>(z >> 11) * (1.0 / 9007199254740992.0);
  }
}

}  // namespace

extern "C" {

int clab_pgsig(const double *x, int n, int cls, int urb, double *sy, double *sz) {
  if (n < 0 || cls < 0 || cls > 5) return CLAB_EARG;
  for (int i = 0; i < n; ++i) sig(x[i], cls, urb, sy[i], sz[i]);
  return CLAB_OK;
}

int clab_briggs(const double *x, int n, double u, double diam, double ev, double ts, double ta, int cls,
                double dthdz, double g, double *out, double *rise) {
  if (n < 0 || !(u > 0) || cls < 0 || cls > 5 || !(ts > 0) || !(ta > 0)) return CLAB_EARG;
  const double F = g * ev * (diam * diam) * (ts - ta) / (4 * ts);
  double fin, xf;
  if (cls >= 4) {
    const double lapse = std::isnan(dthdz) ? (cls == 4 ? 0.020 : 0.035) : dthdz;
    const double s = g * lapse / ta;
    fin = 2.6 * std::pow(F / (u * s), 1.0 / 3);
    xf = 2.0715 * u / std::sqrt(s);
  } else if (F < 55) {
    fin = 21.425 * std::pow(F, 0.75) / u;
    xf = 49 * std::pow(F, 5.0 / 8);
  } else {
    fin = 38.71 * std::pow(F, 0.6) / u;
    xf = 119 * std::pow(F, 0.4);
  }
  out[0] = F;
  out[1] = fin;
  out[2] = xf;
  for (int i = 0; i < n; ++i) {
    const double tr = 1.6 * std::pow(F, 1.0 / 3) * std::pow(x[i], 2.0 / 3) / u;
    rise[i] = x[i] >= xf ? fin : (tr < fin ? tr : fin);
  }
  return CLAB_OK;
}

int clab_gplume(double q, double u, double h, const double *r, int n, int cls, int urb, double lid, int nimg,
                double *c) {
  if (n < 0 || !(u > 0) || cls < 0 || cls > 5 || nimg < 0) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (kutil::stop(1e3)) return CLAB_EINTR;
    const double x = r[i], y = r[static_cast<size_t>(n) + i], z = r[2 * static_cast<size_t>(n) + i];
    if (x <= 0) {
      c[i] = 0;
      continue;
    }
    double sy, sz;
    sig(x, cls, urb, sy, sz);
    c[i] = q / (2 * pi * u * sy * sz) * std::exp(-(y * y) / (2 * sy * sy)) * vert(z, h, sz, lid, nimg);
  }
  return CLAB_OK;
}

int clab_gspuff(double mass, double u, double h, const double *r, int n, double t, int cls, int urb, double lid,
                int nimg, double *c) {
  if (n < 0 || !(u > 0) || !(t > 0) || cls < 0 || cls > 5 || nimg < 0) return CLAB_EARG;
  double sy, sz;
  sig(u * t, cls, urb, sy, sz);
  const double k = mass / (std::pow(2 * pi, 1.5) * sy * sy * sz);
  for (int i = 0; i < n; ++i) {
    if (kutil::stop(1e3)) return CLAB_EINTR;
    const double x = r[i], y = r[static_cast<size_t>(n) + i], z = r[2 * static_cast<size_t>(n) + i];
    c[i] = k * std::exp(-sq(x - u * t) / (2 * sy * sy)) * std::exp(-(y * y) / (2 * sy * sy)) *
           vert(z, h, sz, lid, nimg);
  }
  return CLAB_OK;
}

int clab_advdif(double *cm, int nx, int ny, double u, double v, double kx, double ky, double dx, double dy, double dt,
                double nsteps, const double *src) {
  if (nx < 1 || ny < 1 || !(dx > 0) || !(dy > 0) || !(dt > 0) || !(nsteps >= 0)) return CLAB_EARG;
  const double ax = u * dt / dx, ay = v * dt / dy, dxn = kx * dt / (dx * dx), dyn = ky * dt / (dy * dy);
  const size_t m = static_cast<size_t>(nx) * ny;
  std::vector<double> nw(m);
  auto at = [&](int i, int j) { return cm[static_cast<size_t>(j) * nx + i]; };
  for (double s = 0; s < nsteps; ++s) {
    if (kutil::stop(static_cast<double>(m))) return CLAB_EINTR;
    std::fill(nw.begin(), nw.end(), 0.0);
    for (int i = 1; i + 1 < nx; ++i) {
      for (int j = 1; j + 1 < ny; ++j) {
        const double c = at(i, j);
        const double adx = u >= 0 ? ax * (c - at(i - 1, j)) : ax * (at(i + 1, j) - c);
        const double ady = v >= 0 ? ay * (c - at(i, j - 1)) : ay * (at(i, j + 1) - c);
        const double dif = dxn * (at(i + 1, j) - 2 * c + at(i - 1, j)) + dyn * (at(i, j + 1) - 2 * c + at(i, j - 1));
        const double so = src != nullptr ? src[static_cast<size_t>(j) * nx + i] * dt : 0.0;
        nw[static_cast<size_t>(j) * nx + i] = c - adx - ady + dif + so;
      }
    }
    std::copy(nw.begin(), nw.end(), cm);
  }
  return CLAB_OK;
}

int clab_lagr(int n, double x0, double y0, double u, double v, double kx, double ky, double dt, double nsteps,
              double seed, double *xs, double *ys) {
  if (n < 1 || !(dt > 0) || !(nsteps >= 0) || kx < 0 || ky < 0) return CLAB_EARG;
  const double fx = std::sqrt(2 * kx * dt), fy = std::sqrt(2 * ky * dt), eps = 2.220446049250313e-16;
  for (int i = 0; i < n; ++i) {
    xs[i] = x0;
    ys[i] = y0;
  }
  std::vector<double> z(static_cast<size_t>(n));
  for (double k = 0; k < nsteps; ++k) {
    if (kutil::stop(2.0 * n)) return CLAB_EINTR;
    for (int d = 0; d < 2; ++d) {
      const double s = rmod(seed * 65536 + (2 * k + d), 2147483647);
      smix(static_cast<std::uint64_t>(static_cast<std::int64_t>(s)), static_cast<size_t>(n), z.data());
      double *p = d == 0 ? xs : ys;
      const double step = (d == 0 ? u : v) * dt, f = d == 0 ? fx : fy;
      for (int i = 0; i < n; ++i) {
        const double w = std::min(std::max(z[static_cast<size_t>(i)], eps), 1 - eps);
        p[i] = p[i] + step + f * clab_qnorm(w);
      }
    }
  }
  return CLAB_OK;
}

int clab_lagstep(int n, double u, double v, double fx, double fy, double dt, const double *zx, const double *zy,
                 double *xs, double *ys) {
  if (n < 0) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    xs[i] = xs[i] + u * dt + fx * zx[i];
    ys[i] = ys[i] + v * dt + fy * zy[i];
  }
  return CLAB_OK;
}

int clab_laggrid(const double *xs, const double *ys, int n, const double *g, double *conc, double *mean) {
  if (n < 1) return CLAB_EARG;
  long double sx = 0, sy = 0;
  for (int i = 0; i < n; ++i) {
    sx += xs[i];
    sy += ys[i];
  }
  mean[0] = static_cast<double>(sx / n);
  mean[1] = static_cast<double>(sy / n);
  if (g == nullptr) return CLAB_OK;
  const int gx = static_cast<int>(g[4]), gy = static_cast<int>(g[5]);
  const double wx = (g[1] - g[0]) / g[4], wy = (g[3] - g[2]) / g[5];
  std::fill(conc, conc + static_cast<size_t>(gx) * gy, 0.0);
  for (int p = 0; p < n; ++p) {
    const double ii = std::floor((xs[p] - g[0]) / wx), jj = std::floor((ys[p] - g[2]) / wy);
    if (ii >= 0 && ii < g[4] && jj >= 0 && jj < g[5]) {
      conc[static_cast<size_t>(jj) * gx + static_cast<size_t>(ii)] += 1 / (n * wx * wy);
    }
  }
  return CLAB_OK;
}

// Concentration-response: form 0 the log-linear exp(beta e / 10), 1 the IER
// 1 + a (1 - exp(-g e^d)), e the excess over refc. rr, lrr (n each).
int clab_crf(const double *z, int n, double refc, int form, const double *par, double *rr, double *lrr) {
  if (n < 0 || form < 0 || form > 1) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    const double e = z[i] - refc > 0 ? z[i] - refc : 0.0;
    if (form == 0) {
      lrr[i] = par[0] * e / 10;
      rr[i] = std::exp(lrr[i]);
    } else {
      rr[i] = 1 + par[0] * (1 - std::exp(-par[1] * std::pow(e, par[2])));
      lrr[i] = std::log(rr[i]);
    }
  }
  return CLAB_OK;
}

double clab_levin(double rr, double p) { return p * (rr - 1) / (1 + p * (rr - 1)); }

double clab_hif(double dx, double pop, double rate, double beta) { return rate * pop * (1 - std::exp(-beta * dx)); }

int clab_concidx(const double *h, const double *inc, int n, double *out, double *cpop, double *cexp) {
  if (n < 2) return CLAB_EARG;
  std::vector<int> o(static_cast<size_t>(n));
  std::iota(o.begin(), o.end(), 0);
  std::stable_sort(o.begin(), o.end(), [&](int a, int b) { return inc[a] < inc[b]; });
  // mid-ranks of the incomes
  std::vector<double> rk(static_cast<size_t>(n));
  for (int i = 0; i < n;) {
    int j = i;
    while (j + 1 < n && inc[o[static_cast<size_t>(j + 1)]] == inc[o[static_cast<size_t>(i)]]) ++j;
    const double r = (i + j) / 2.0 + 1;
    for (int k = i; k <= j; ++k) rk[static_cast<size_t>(o[static_cast<size_t>(k)])] = r;
    i = j + 1;
  }
  long double sh = 0, sr = 0;
  std::vector<double> R(static_cast<size_t>(n));
  for (int i = 0; i < n; ++i) {
    R[static_cast<size_t>(i)] = (rk[static_cast<size_t>(i)] - 0.5) / n;
    sh += h[i];
    sr += R[static_cast<size_t>(i)];
  }
  // R's mean(): the sum over n, then one correction pass
  double mh = static_cast<double>(sh / n), mr = static_cast<double>(sr / n);
  long double th = 0, tr = 0;
  for (int i = 0; i < n; ++i) {
    th += h[i] - mh;
    tr += R[static_cast<size_t>(i)] - mr;
  }
  mh = static_cast<double>(mh + th / n);
  mr = static_cast<double>(mr + tr / n);
  long double cv = 0;
  for (int i = 0; i < n; ++i) cv += (h[i] - mh) * (R[static_cast<size_t>(i)] - mr);
  const double cov = static_cast<double>(cv) / n;
  out[0] = 2 * cov / mh;
  out[1] = cov;
  out[2] = mh;
  // the concentration curve, poorest first, ties broken by exposure
  std::stable_sort(o.begin(), o.end(), [&](int a, int b) { return inc[a] < inc[b] || (inc[a] == inc[b] && h[a] < h[b]); });
  long double run = 0;
  for (int i = 0; i < n; ++i) {
    run += h[o[static_cast<size_t>(i)]];
    cpop[i] = static_cast<double>(i + 1) / n;
    cexp[i] = static_cast<double>(run) / static_cast<double>(sh);
  }
  return CLAB_OK;
}

}  // extern "C"
