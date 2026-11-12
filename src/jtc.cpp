// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, journey to crime: distance-decay functions, the surface they
// give over a grid of candidate bases for a series of incidents, kernel
// calibration of the decay from known trips, accuracy measures, and the
// Bayesian combination with an origin-destination matrix.
#include <cmath>
#include <limits>

#include "cgeo.h"
#include "crimlab.h"
#include "kutil.h"

namespace {

const double ninf = -std::numeric_limits<double>::infinity();
const double ls2pi = 0.91893853320467274178;  // log(sqrt(2 pi))

// validity of the parameters of a decay function
bool parok(int type, const double *p, const double *tab, int nt) {
  switch (type) {
    case 0:
      return std::isfinite(p[0]) && std::isfinite(p[1]);
    case 1:
      return kutil::pos(p[0]) && kutil::nonneg(p[1]);
    case 2:
    case 3:
      return kutil::pos(p[0]) && std::isfinite(p[1]) && kutil::pos(p[2]);
    case 4:
      return kutil::pos(p[0]) && kutil::pos(p[1]) && kutil::nonneg(p[2]);
    default:  // 5, a table (the type is checked by the callers)
      if (!kutil::pos(p[0]) || tab == nullptr || nt < 2) return false;
      for (int k = 0; k < nt; ++k)
        if (!kutil::nonneg(tab[k])) return false;
      return true;
  }
}

// log f(d) for d >= 0 (-inf where f is 0)
double lgdec(int type, const double *p, const double *tab, int nt, double d) {
  switch (type) {
    case 0: {
      const double v = p[0] + p[1] * d;
      return v > 0.0 ? std::log(v) : ninf;
    }
    case 1:
      return std::log(p[0]) - p[1] * d;
    case 2: {
      const double z = (d - p[1]) / p[2];
      return std::log(p[0]) - std::log(p[2]) - ls2pi - 0.5 * z * z;
    }
    case 3: {
      if (d == 0.0) return ninf;
      const double z = (std::log(d) - p[1]) / p[2];
      return std::log(p[0]) - std::log(d) - std::log(p[2]) - ls2pi -
             0.5 * z * z;
    }
    case 4:
      if (d <= p[1]) return d > 0.0 ? std::log(p[0] * (d / p[1])) : ninf;
      return std::log(p[0]) - p[2] * (d - p[1]);
    default: {  // table, linear between nodes k * delta, 0 beyond the last
      const double u = d / p[0];
      if (u > nt - 1) return ninf;
      int k = static_cast<int>(u);
      if (k >= nt - 1) k = nt - 2;
      const double t = u - k;
      const double v = tab[k] + t * (tab[k + 1] - tab[k]);
      return v > 0.0 ? std::log(v) : ninf;
    }
  }
}

bool gridok(int nx, int ny, double x0, double y0, double dx, double dy) {
  return nx >= 1 && ny >= 1 && nx <= 1000000000 / ny && std::isfinite(x0) &&
         std::isfinite(y0) && kutil::pos(dx) && kutil::pos(dy);
}

}  // namespace

extern "C" {

int clab_jtdecay(const double *d, int n, int type, const double *par,
                 const double *tab, int nt, double *f) {
  if (n < 0 || type < 0 || type > 5) return CLAB_EARG;
  if (!parok(type, par, tab, nt)) return CLAB_EDOM;
  for (int i = 0; i < n; ++i)
    if (!kutil::nonneg(d[i])) return CLAB_EDOM;
  for (int i = 0; i < n; ++i) f[i] = std::exp(lgdec(type, par, tab, nt, d[i]));
  return CLAB_OK;
}

int clab_jtsurf(const double *x, const double *y, int m, int nx, int ny,
                double x0, double y0, double dx, double dy, int metric,
                int mode, int type, const double *par, const double *tab,
                int nt, double *raw, double *prob, double *peak) {
  if (m < 1 || metric < 0 || metric > 1 || mode < 0 || mode > 1 || type < 0 ||
      type > 5)
    return CLAB_EARG;
  if (!gridok(nx, ny, x0, y0, dx, dy)) return CLAB_EARG;
  if (!parok(type, par, tab, nt)) return CLAB_EDOM;
  for (int i = 0; i < m; ++i)
    if (!std::isfinite(x[i]) || !std::isfinite(y[i])) return CLAB_EDOM;
  const int nc = nx * ny;
  for (int c = 0; c < nc; ++c) {
    const double cx = x0 + (c % nx + 0.5) * dx, cy = y0 + (c / nx + 0.5) * dy;
    kutil::ksum s;
    bool dead = false;
    for (int i = 0; i < m; ++i) {
      const double ax = std::fabs(x[i] - cx), ay = std::fabs(y[i] - cy);
      const double d = metric == 0 ? std::hypot(ax, ay) : ax + ay;
      const double l = lgdec(type, par, tab, nt, d);
      if (mode == 0) {
        s.add(std::exp(l));
      } else if (l == ninf) {
        dead = true;
      } else {
        s.add(l);
      }
    }
    raw[c] = (mode == 1 && dead) ? ninf : s.val();
  }
  // normalise: sum of the surface, or of exp(log-likelihood - max)
  int arg = 0;
  for (int c = 1; c < nc; ++c)
    if (raw[c] > raw[arg]) arg = c;
  const double top = raw[arg];
  if (mode == 0 ? !(top > 0.0) : top == ninf) return CLAB_EDOM;
  if (!std::isfinite(top)) return CLAB_EDOM;
  kutil::ksum tot;
  for (int c = 0; c < nc; ++c) {
    prob[c] = mode == 0 ? raw[c] : std::exp(raw[c] - top);
    tot.add(prob[c]);
  }
  const double t = tot.val();
  for (int c = 0; c < nc; ++c) prob[c] /= t;
  peak[0] = arg;
  peak[1] = x0 + (arg % nx + 0.5) * dx;
  peak[2] = y0 + (arg / nx + 0.5) * dy;
  return CLAB_OK;
}

int clab_jtcal(const double *d, const double *w, int m, int kern, double h,
               int refl, const double *t, int k, double *dens) {
  if (m < 1 || k < 0 || kern < 0 || kern > 4 || refl < 0 || refl > 1)
    return CLAB_EARG;
  if (!kutil::pos(h)) return CLAB_EDOM;
  kutil::ksum sw;
  for (int i = 0; i < m; ++i) {
    if (!kutil::nonneg(d[i])) return CLAB_EDOM;
    const double wi = w ? w[i] : 1.0;
    if (!kutil::nonneg(wi)) return CLAB_EDOM;
    sw.add(wi);
  }
  const double W = sw.val();
  if (!(W > 0.0) || !std::isfinite(W)) return CLAB_EDOM;
  for (int j = 0; j < k; ++j)
    if (!std::isfinite(t[j]) || (refl && t[j] < 0.0)) return CLAB_EDOM;
  auto K = [kern](double u) {
    const double a = std::fabs(u);
    switch (kern) {
      case 0:
        return std::exp(-0.5 * u * u - ls2pi);
      case 1:
        return a < 1.0 ? 0.9375 * (1.0 - u * u) * (1.0 - u * u) : 0.0;
      case 2:
        return a < 1.0 ? 1.0 - a : 0.0;
      case 3:
        return a <= 1.0 ? 0.5 : 0.0;
      default:
        return 0.5 * std::exp(-a);
    }
  };
  for (int j = 0; j < k; ++j) {
    kutil::ksum s;
    for (int i = 0; i < m; ++i) {
      const double wi = w ? w[i] : 1.0;
      if (wi == 0.0) continue;
      double v = K((t[j] - d[i]) / h);
      if (refl) v += K((t[j] + d[i]) / h);
      s.add(wi * v);
    }
    dens[j] = s.val() / (h * W);
  }
  return CLAB_OK;
}

int clab_jtacc(const double *prob, int nx, int ny, double x0, double y0,
               double dx, double dy, double hx, double hy, double *out) {
  if (!gridok(nx, ny, x0, y0, dx, dy)) return CLAB_EARG;
  if (!std::isfinite(hx) || !std::isfinite(hy)) return CLAB_EDOM;
  const double fx = (hx - x0) / dx, fy = (hy - y0) / dy;
  if (!(fx >= 0.0 && fx <= nx && fy >= 0.0 && fy <= ny)) return CLAB_EDOM;
  const int ix = fx >= nx ? nx - 1 : static_cast<int>(fx);
  const int iy = fy >= ny ? ny - 1 : static_cast<int>(fy);
  const int nc = nx * ny, home = iy * nx + ix;
  for (int c = 0; c < nc; ++c)
    if (!kutil::nonneg(prob[c])) return CLAB_EDOM;
  const double ph = prob[home];
  int gt = 0, ge = 0, arg = 0;
  for (int c = 0; c < nc; ++c) {
    if (prob[c] > ph) ++gt;
    if (prob[c] >= ph) ++ge;
    if (prob[c] > prob[arg]) arg = c;
  }
  out[0] = ph;
  out[1] = (gt + 1.0) / nc;
  out[2] = static_cast<double>(ge) / nc;
  out[3] = std::hypot(x0 + (arg % nx + 0.5) * dx - hx,
                      y0 + (arg / nx + 0.5) * dy - hy);
  out[4] = arg;
  out[5] = home;
  return CLAB_OK;
}

int clab_bjtc(const double *pj, int no, const double *od, int nd,
              const int *dest, int m, double *out) {
  if (no < 1 || nd < 1 || m < 1 || no > 1000000000 / nd) return CLAB_EARG;
  for (int k = 0; k < m; ++k)
    if (dest[k] < 0 || dest[k] >= nd) return CLAB_EARG;
  for (int o = 0; o < no; ++o) {
    if (!kutil::nonneg(pj[o])) return CLAB_EDOM;
    for (int j = 0; j < nd; ++j)
      if (!kutil::nonneg(od[static_cast<long>(o) * nd + j])) return CLAB_EDOM;
  }
  double *pjt = out, *gen = out + no, *cnd = out + 2 * no, *prd = out + 3 * no,
         *rsk = out + 4 * no;
  for (int o = 0; o < no; ++o) {
    kutil::ksum g, c;
    const double *row = od + static_cast<long>(o) * nd;
    for (int j = 0; j < nd; ++j) g.add(row[j]);
    for (int k = 0; k < m; ++k) c.add(row[dest[k]]);  // repeats count
    pjt[o] = pj[o];
    gen[o] = g.val();
    cnd[o] = c.val();
  }
  auto norm = [no](double *v) {
    kutil::ksum s;
    for (int o = 0; o < no; ++o) s.add(v[o]);
    const double t = s.val();
    if (!(t > 0.0) || !std::isfinite(t)) return false;
    for (int o = 0; o < no; ++o) v[o] /= t;
    return true;
  };
  if (!norm(pjt) || !norm(gen) || !norm(cnd)) return CLAB_EDOM;
  for (int o = 0; o < no; ++o) {
    prd[o] = pjt[o] * cnd[o];
    // where no offender lives the product is 0 too (cnd <= gen): risk 0
    rsk[o] = gen[o] > 0.0 ? prd[o] / gen[o] : 0.0;
  }
  if (!norm(prd) || !norm(rsk)) return CLAB_EDOM;
  return CLAB_OK;
}

}  // extern "C"
