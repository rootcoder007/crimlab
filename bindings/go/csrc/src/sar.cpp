// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Search theory for search and rescue planning (Koopman 1946, 1956, 1957;
// Frost 1996; Stone 1975): lateral range curves, coverage and the
// probability of detection, parallel sweeps with navigation error, the
// probability of success with Bayesian updating of containment, datum
// probability maps, the optimal allocation of effort and Frost's optimal
// search square. Nothing here allocates on the caller's behalf.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

namespace {

using std::size_t;
const double pi = 3.14159265358979323846;

double pn(double x) { return clab_pnorm(x); }

// The normal density with standard deviation sd at x, split for accuracy in
// the tails (Loader's device: exp(-x^2/2) as a product of two exponentials
// of an exactly representable part and the rest).
double dn(double x, double sd) {
  const double u = std::fabs(x / sd);
  if (u >= 2 * std::sqrt(1.7976931348623157e308)) return 0;
  const double k = 0.398942280401432677939946059934;  // 1 / sqrt(2 pi)
  if (u < 5) return k * std::exp(-0.5 * u * u) / sd;
  const double u1 = std::ldexp(std::nearbyint(std::ldexp(u, 16)), -16), u2 = u - u1;
  return k / sd * (std::exp(-0.5 * u1 * u1) * std::exp((-0.5 * u2 - u1) * u2));
}

// The lateral range curve: 0 inverse cube, 1 definite range, 2 M-Beta.
double lr(double x, double W, int model, double m) {
  const double a = std::fabs(x);
  if (model == 1) return a <= W / 2 ? 1.0 : 0.0;
  if (model == 2) return a <= W / (2 * m) ? m : 0.0;
  return a == 0 ? 1.0 : 1 - std::exp(-(W * W) / (4 * pi * (a * a)));
}

// The eigenvalues of the symmetric tridiagonal matrix with zero diagonal
// and off-diagonal b (n - 1), with the first component of each unit
// eigenvector, by the implicit QL method (Bowdler, Martin, Reinsch and
// Wilkinson 1968, the procedure tql2), sorted ascending: the Golub-Welsch
// nodes and weights of a Gauss rule.
bool gw(int n, const std::vector<double> &b, std::vector<double> &d, std::vector<double> &z) {
  d.assign(static_cast<size_t>(n), 0.0);
  std::vector<double> e(static_cast<size_t>(n), 0.0);
  for (int i = 0; i + 1 < n; ++i) e[static_cast<size_t>(i)] = b[static_cast<size_t>(i)];
  // the full eigenvector matrix: the rotations mix columns, row by row
  std::vector<double> v(static_cast<size_t>(n) * n, 0.0);
  for (int i = 0; i < n; ++i) v[static_cast<size_t>(i) * n + i] = 1;
  for (int l = 0; l < n; ++l) {
    int it = 0, mm;
    do {
      for (mm = l; mm < n - 1; ++mm) {
        const double dd = std::fabs(d[mm]) + std::fabs(d[mm + 1]);
        if (std::fabs(e[mm]) <= 1e-17 * dd || dd + std::fabs(e[mm]) == dd) break;
      }
      if (mm != l) {
        if (++it > 60) return false;
        double g = (d[l + 1] - d[l]) / (2 * e[l]);
        double r = std::hypot(g, 1.0);
        g = d[mm] - d[l] + e[l] / (g + (g >= 0 ? std::fabs(r) : -std::fabs(r)));
        double s = 1, c = 1, p = 0;
        int i;
        for (i = mm - 1; i >= l; --i) {
          double f = s * e[i];
          const double bb = c * e[i];
          r = std::hypot(f, g);
          e[i + 1] = r;
          if (r == 0) {
            d[i + 1] -= p;
            e[mm] = 0;
            break;
          }
          s = f / r;
          c = g / r;
          g = d[i + 1] - p;
          r = (d[i] - g) * s + 2 * c * bb;
          p = s * r;
          d[i + 1] = g + p;
          g = c * r - bb;
          for (int k = 0; k < n; ++k) {
            f = v[static_cast<size_t>(i + 1) * n + k];
            v[static_cast<size_t>(i + 1) * n + k] = s * v[static_cast<size_t>(i) * n + k] + c * f;
            v[static_cast<size_t>(i) * n + k] = c * v[static_cast<size_t>(i) * n + k] - s * f;
          }
        }
        if (r == 0 && i >= l) continue;
        d[l] -= p;
        e[l] = g;
        e[mm] = 0;
      }
    } while (mm != l);
  }
  // sort ascending, carrying the first components
  std::vector<int> o(static_cast<size_t>(n));
  for (int i = 0; i < n; ++i) o[static_cast<size_t>(i)] = i;
  for (int i = 1; i < n; ++i) {
    for (int j = i; j > 0 && d[static_cast<size_t>(o[j])] < d[static_cast<size_t>(o[j - 1])]; --j) {
      std::swap(o[j], o[j - 1]);
    }
  }
  std::vector<double> ds(static_cast<size_t>(n));
  z.assign(static_cast<size_t>(n), 0.0);
  for (int i = 0; i < n; ++i) {
    ds[static_cast<size_t>(i)] = d[static_cast<size_t>(o[i])];
    z[static_cast<size_t>(i)] = v[static_cast<size_t>(o[i]) * n];
  }
  d = ds;
  return true;
}

// Brent's minimiser on [a, b] to tol (Brent 1973, chapter 5, procedure
// localmin), with the relative tolerance sqrt(DBL_EPSILON) of the book.
template <class F>
double bmin(F f, double a, double b, double tol) {
  const double c = (3 - std::sqrt(5.0)) * 0.5, eps = std::sqrt(2.220446049250313e-16);
  double v = a + c * (b - a), w = v, x = v, d = 0, e = 0;
  double fx = f(x), fv = fx, fw = fx;
  const double t3 = tol / 3;
  for (int it = 0; it < 1000; ++it) {
    const double xm = (a + b) * 0.5, t1 = eps * std::fabs(x) + t3, t2 = t1 * 2;
    if (std::fabs(x - xm) <= t2 - (b - a) * 0.5) break;
    double p = 0, q = 0, r = 0;
    if (std::fabs(e) > t1) {
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = (q - r) * 2;
      if (q > 0) p = -p;
      else q = -q;
      r = e;
      e = d;
    }
    if (std::fabs(p) >= std::fabs(q * 0.5 * r) || p <= q * (a - x) || p >= q * (b - x)) {
      e = x < xm ? b - x : a - x;
      d = c * e;
    } else {
      d = p / q;
      const double u = x + d;
      if (u - a < t2 || b - u < t2) d = x >= xm ? -t1 : t1;
    }
    const double u = std::fabs(d) >= t1 ? x + d : (d > 0 ? x + t1 : x - t1);
    const double fu = f(u);
    if (fu <= fx) {
      if (u < x) b = x;
      else a = x;
      v = w;
      w = x;
      x = u;
      fv = fw;
      fw = fx;
      fx = fu;
    } else {
      if (u < x) a = u;
      else b = u;
      if (fu <= fw || w == x) {
        v = w;
        fv = fw;
        w = u;
        fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v = u;
        fv = fu;
      }
    }
  }
  return x;
}

double pod1(double cv, int model, double m) {
  switch (model) {
    case 0: return 1 - std::exp(-cv);
    case 1: return cv < 1 ? cv : 1.0;
    case 2: return 2 * pn(std::sqrt(pi) * cv / 2 * std::sqrt(2.0)) - 1;
    default: return cv < m ? cv : 1 - std::pow(1 - m, cv / m);
  }
}

// One cell's mass of a circular normal about (px, py).
void cell(const double *xe, int nx, const double *ye, int ny, double s, double px, double py, double wt,
          double *poc) {
  for (int j = 0; j < nx; ++j) {
    const double dx = pn((xe[j + 1] - px) / s) - pn((xe[j] - px) / s);
    for (int i = 0; i < ny; ++i) {
      const double dy = pn((ye[i + 1] - py) / s) - pn((ye[i] - py) / s);
      poc[static_cast<size_t>(j) * ny + i] += wt * (dy * dx);
    }
  }
}

// The share of [e[i], e[i + 1]] inside [a, b], over b - a.
double along(const double *e, int i, double a, double b) {
  const double hi = e[i + 1] < b ? e[i + 1] : b, lo = e[i] > a ? e[i] : a;
  const double o = hi - lo;
  return (o > 0 ? o : 0.0) / (b - a);
}

}  // namespace

extern "C" {

int clab_sargs(int n, int kind, double *t, double *w) CLAB_TRY {
  if (n < 2 || n > 1000 || t == nullptr || w == nullptr) return CLAB_EARG;
  std::vector<double> b(static_cast<size_t>(n - 1)), d, z;
  for (int i = 1; i < n; ++i) {
    const double di = i;
    b[static_cast<size_t>(i - 1)] = kind == 0 ? std::sqrt(di / 2) : di / std::sqrt(4 * di * di - 1);
  }
  if (!gw(n, b, d, z)) return CLAB_EDOM;
  for (int i = 0; i < n; ++i) {
    t[i] = d[static_cast<size_t>(i)];
    w[i] = (kind == 0 ? std::sqrt(pi) : 2.0) * (z[static_cast<size_t>(i)] * z[static_cast<size_t>(i)]);
  }
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_sarlr(const double *x, int n, double W, int model, double m, double navsd, double *p) CLAB_TRY {
  if (n < 0 || !(W > 0) || model < 0 || model > 2 || !(m > 0 && m <= 1) || !(navsd >= 0)) return CLAB_EARG;
  if (navsd == 0) {
    for (int i = 0; i < n; ++i) p[i] = lr(x[i], W, model, m);
    return CLAB_OK;
  }
  if (model != 0) {
    const double h = model == 2 ? m : 1.0, r = W / (2 * h);
    for (int i = 0; i < n; ++i) p[i] = h * (pn((x[i] + r) / navsd) - pn((x[i] - r) / navsd));
    return CLAB_OK;
  }
  // the trapezoid rule on 641 points over 8 standard deviations either side
  const int k = 641;
  std::vector<double> t(k), w(k);
  for (int j = 0; j < k; ++j) t[static_cast<size_t>(j)] = (j == k - 1 ? 8.0 : -8 + j * (16.0 / (k - 1))) * navsd;
  const double h = t[1] - t[0];
  for (int j = 0; j < k; ++j) w[static_cast<size_t>(j)] = dn(t[static_cast<size_t>(j)], navsd) * h;
  for (int i = 0; i < n; ++i) {
    double s = 0;
    for (int j = 0; j < k; ++j) s += lr(x[i] - t[static_cast<size_t>(j)], W, model, m) * w[static_cast<size_t>(j)];
    p[i] = s;
  }
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_sarsw(const double *x, const double *p, int n, double *w) CLAB_TRY {
  if (n < 2) return CLAB_EARG;
  double s = 0;
  for (int i = 0; i + 1 < n; ++i) s += (x[i + 1] - x[i]) * (p[i + 1] + p[i]) / 2;
  *w = s;
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_sarpod(const double *cv, int n, int model, double m, double *pod) CLAB_TRY {
  if (n < 0 || model < 0 || model > 3 || !(m > 0 && m <= 1)) return CLAB_EARG;
  for (int i = 0; i < n; ++i) pod[i] = pod1(cv[i], model, m);
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_sarpsw(double W, double S, int model, double m, double navsd, int npts, int ntrk, double *pod) CLAB_TRY {
  if (!(W > 0) || !(S > 0) || model < 0 || model > 2 || !(m > 0 && m <= 1) || !(navsd >= 0) || npts < 3 ||
      ntrk < 1) {
    return CLAB_EARG;
  }
  if (npts % 2 == 0) ++npts;
  std::vector<double> gt, gwt;
  if (navsd > 0) {
    gt.resize(40);
    gwt.resize(40);
    const int rc = clab_sargs(40, 0, gt.data(), gwt.data());
    if (rc != CLAB_OK) return rc;
  }
  const double sq = std::sqrt(2.0) * navsd, spi = std::sqrt(pi);
  auto lrn = [&](double d) {
    if (gt.empty()) return lr(d, W, model, m);
    double s = 0;
    for (size_t j = 0; j < gt.size(); ++j) s += lr(d + sq * gt[j], W, model, m) * gwt[j];
    return s / spi;
  };
  const double h = S / (npts - 1);
  std::vector<double> vals(static_cast<size_t>(npts));
  for (int i = 0; i < npts; ++i) {
    if (kutil::stop(static_cast<double>(2 * ntrk + 1) * (gt.empty() ? 1 : 40))) return CLAB_EINTR;
    const double x = i * h;
    double miss = 1, near = 0;
    for (int k = -ntrk; k <= ntrk; ++k) {
      const double dd = x - k * S;
      miss *= 1 - lrn(dd);
      if (dd != 0) near += 1 / (dd * dd);
    }
    if (model == 0 && miss > 0) {
      const double sn = std::sin(pi * x / S);
      const double full = sn != 0 ? (pi / S) * (pi / S) / (sn * sn) : near;
      const double ex = full - near > 0 ? full - near : 0.0;
      miss *= std::exp(-(W * W) / (4 * pi) * ex);
    }
    vals[static_cast<size_t>(i)] = 1 - miss;
  }
  double so = 0, se = 0;
  for (int i = 1; i <= npts - 2; i += 2) so += vals[static_cast<size_t>(i)];
  for (int i = 2; i <= npts - 3; i += 2) se += vals[static_cast<size_t>(i)];
  *pod = (vals[0] + vals[static_cast<size_t>(npts - 1)] + 4 * so + 2 * se) * h / 3 / S;
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_sarscs(const double *poc, int na, const double *D, int ns, double *pos, double *posarea, double *pocupd,
                double *cumpod, int *certain) CLAB_TRY {
  if (na < 1 || ns < 1) return CLAB_EARG;
  *certain = 0;
  for (int j = 0; j < na; ++j) pocupd[j] = poc[j];
  for (int k = 0; k < ns; ++k) {
    double p = 0;
    for (int j = 0; j < na; ++j) p += pocupd[j] * D[static_cast<size_t>(j) * ns + k];
    if (1 - p <= 0) {
      *certain = 1;
      for (int j = 0; j < na; ++j) pocupd[j] = kutil::nan();
      break;
    }
    for (int j = 0; j < na; ++j) pocupd[j] = pocupd[j] * (1 - D[static_cast<size_t>(j) * ns + k]) / (1 - p);
  }
  double t = 0;
  for (int j = 0; j < na; ++j) {
    double q = 1;
    for (int k = 0; k < ns; ++k) q *= 1 - D[static_cast<size_t>(j) * ns + k];
    cumpod[j] = 1 - q;
    posarea[j] = poc[j] * cumpod[j];
    t += posarea[j];
  }
  *pos = t;
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_sarmap(const double *xe, int nxe, const double *ye, int nye, int datum, const double *cen, double sigma,
                const double *end, int nodes, const double *wts, double *poc) CLAB_TRY {
  if (nxe < 2 || nye < 2 || datum < 0 || datum > 4) return CLAB_EARG;
  const int nx = nxe - 1, ny = nye - 1;
  const size_t nc = static_cast<size_t>(nx) * ny;
  for (size_t i = 0; i < nc; ++i) poc[i] = 0;
  switch (datum) {
    case 0: cell(xe, nx, ye, ny, sigma, cen[0], cen[1], 1.0, poc); break;
    case 1: {
      if (nodes < 2 || nodes > 1000) return CLAB_EARG;
      std::vector<double> t(static_cast<size_t>(nodes)), w(static_cast<size_t>(nodes));
      const int rc = clab_sargs(nodes, 1, t.data(), w.data());
      if (rc != CLAB_OK) return rc;
      std::vector<double> one(nc);
      for (int k = 0; k < nodes; ++k) {
        const double f = (t[static_cast<size_t>(k)] + 1) / 2;
        for (size_t i = 0; i < nc; ++i) one[i] = 0;
        cell(xe, nx, ye, ny, sigma, cen[0] + f * (end[0] - cen[0]), cen[1] + f * (end[1] - cen[1]), 1.0,
             one.data());
        for (size_t i = 0; i < nc; ++i) poc[i] += w[static_cast<size_t>(k)] / 2 * one[i];
      }
      break;
    }
    case 2: {
      const double x0 = cen[0] < end[0] ? cen[0] : end[0], x1 = cen[0] < end[0] ? end[0] : cen[0];
      const double y0 = cen[1] < end[1] ? cen[1] : end[1], y1 = cen[1] < end[1] ? end[1] : cen[1];
      if (!(x1 - x0 > 0) || !(y1 - y0 > 0)) return CLAB_EDOM;
      for (int j = 0; j < nx; ++j) {
        const double hx = xe[j + 1] < x1 ? xe[j + 1] : x1, lx = xe[j] > x0 ? xe[j] : x0;
        const double ox = hx - lx > 0 ? hx - lx : 0.0;
        for (int i = 0; i < ny; ++i) {
          const double hy = ye[i + 1] < y1 ? ye[i + 1] : y1, ly = ye[i] > y0 ? ye[i] : y0;
          const double oy = hy - ly > 0 ? hy - ly : 0.0;
          poc[static_cast<size_t>(j) * ny + i] = oy * ox / ((x1 - x0) * (y1 - y0));
        }
      }
      break;
    }
    case 3: {
      const bool hz = cen[1] == end[1];
      if (!hz && cen[0] != end[0]) return CLAB_EDOM;
      if (cen[0] == end[0] && cen[1] == end[1]) return CLAB_EDOM;
      if (hz) {
        const double a = cen[0] < end[0] ? cen[0] : end[0], b = cen[0] < end[0] ? end[0] : cen[0];
        for (int j = 0; j < nx; ++j) {
          const double al = along(xe, j, a, b);
          for (int i = 0; i < ny; ++i) {
            const double ac = pn((ye[i + 1] - cen[1]) / sigma) - pn((ye[i] - cen[1]) / sigma);
            poc[static_cast<size_t>(j) * ny + i] = ac * al;
          }
        }
      } else {
        const double a = cen[1] < end[1] ? cen[1] : end[1], b = cen[1] < end[1] ? end[1] : cen[1];
        for (int j = 0; j < nx; ++j) {
          const double ac = pn((xe[j + 1] - cen[0]) / sigma) - pn((xe[j] - cen[0]) / sigma);
          for (int i = 0; i < ny; ++i) poc[static_cast<size_t>(j) * ny + i] = along(ye, i, a, b) * ac;
        }
      }
      break;
    }
    default: {
      double s = 0;
      for (size_t i = 0; i < nc; ++i) {
        if (!(wts[i] >= 0) || !std::isfinite(wts[i])) return CLAB_EDOM;
        s += wts[i];
      }
      if (!(s > 0)) return CLAB_EDOM;
      for (size_t i = 0; i < nc; ++i) poc[i] = wts[i] / s;
    }
  }
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_sarcrc(double sigma, double W, double effort, double *out) CLAB_TRY {
  if (!(sigma > 0) || !(W > 0) || !(effort >= 0)) return CLAB_EARG;
  const double R = std::pow(4 * (sigma * sigma) * W * effort / pi, 0.25);
  const double u = R * R / (2 * (sigma * sigma));
  out[0] = R;
  out[1] = 1 - (1 + u) * std::exp(-u);
  out[2] = R * R / (2 * (sigma * sigma) * W);
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_saropt(const double *poc, const double *area, int n, double W, double effort, double tol, double *z,
                double *cov, double *pod, double *pos, double *lam) CLAB_TRY {
  if (n < 1 || !(W > 0) || !(effort >= 0) || !(tol > 0)) return CLAB_EARG;
  std::vector<double> r(static_cast<size_t>(n));
  double rmin = HUGE_VAL, rmax = -HUGE_VAL;
  bool any = false;
  for (int i = 0; i < n; ++i) {
    if (!(area[i] > 0) || !(poc[i] >= 0)) return CLAB_EDOM;
    r[static_cast<size_t>(i)] = poc[i] * W / area[i];
    // a density that underflows to 0 (or overflows) leaves no root to
    // bracket: the search below would step down for ever
    if (!std::isfinite(r[static_cast<size_t>(i)]) || !std::isfinite(area[i] / W) ||
        (poc[i] > 0 && !(r[static_cast<size_t>(i)] > 0))) {
      return CLAB_EDOM;
    }
    if (poc[i] > 0 && r[static_cast<size_t>(i)] < rmin) rmin = r[static_cast<size_t>(i)];
    if (r[static_cast<size_t>(i)] > rmax) rmax = r[static_cast<size_t>(i)];
    any = any || poc[i] > 0;
  }
  if (!any) return CLAB_EDOM;
  auto total = [&](double l) {
    double s = 0;
    for (int i = 0; i < n; ++i) {
      const double q = r[static_cast<size_t>(i)] / l;
      if (q > 1) s += area[i] / W * std::log(q);
    }
    return s;
  };
  double lo = std::log(rmin) - 50;
  while (total(std::exp(lo)) <= effort) lo -= 50;
  double hi = std::log(rmax);
  for (int it = 0; it < 400; ++it) {
    const double mid = (lo + hi) / 2;
    if (total(std::exp(mid)) > effort) lo = mid;
    else hi = mid;
    if (hi - lo < tol) break;
  }
  const double l = std::exp((lo + hi) / 2);
  double s = 0;
  for (int i = 0; i < n; ++i) {
    const double q = r[static_cast<size_t>(i)] / l;
    z[i] = q > 1 ? area[i] / W * std::log(q) : 0.0;
    cov[i] = W * z[i] / area[i];
    pod[i] = 1 - std::exp(-cov[i]);
    s += poc[i] * pod[i];
  }
  *pos = s;
  *lam = l;
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_sarsq(const double *z, int k, int model, double m, double *out) CLAB_TRY {
  if (k < 1 || model < 0 || model > 3 || !(m > 0 && m <= 1)) return CLAB_EARG;
  const double e0 = std::sqrt(2 * std::log(2.0));
  std::vector<double> a, pd;
  // the unnormalised mass inside a centred square of half side h after the
  // squares a searched with PODs pd: shells between the half sides, each
  // scaled by the chance every earlier search covering it missed
  auto mass = [&](double h) {
    std::vector<double> lo = a;
    lo.push_back(0);
    std::sort(lo.begin(), lo.end());
    lo.erase(std::unique(lo.begin(), lo.end()), lo.end());
    double s = 0;
    for (size_t q = 0; q < lo.size(); ++q) {
      const double hi = q + 1 < lo.size() ? lo[q + 1] : HUGE_VAL;
      double miss = 1;
      for (size_t j = 0; j < a.size(); ++j) {
        if (a[j] >= hi) miss *= 1 - pd[j];
      }
      const double bh = 2 * pn(h < hi ? h : hi) - 1, bl = 2 * pn(h < lo[q] ? h : lo[q]) - 1;
      s += (bh * bh - bl * bl) * miss;
    }
    return s;
  };
  double left = 1, cum = 0;
  for (int s = 0; s < k; ++s) {
    if (!(z[s] > 0)) return CLAB_EDOM;
    const double zs = z[s];
    auto fn = [&](double f) { return mass(f * e0) * pod1(zs / (4 * (f * f)), model, m); };
    const int ng = 241;
    const double l0 = std::log(0.02), l1 = std::log(10.0), by = (l1 - l0) / (ng - 1);
    std::vector<double> fs(ng);
    int best = 0;
    double bv = -HUGE_VAL;
    for (int i = 0; i < ng; ++i) {
      fs[static_cast<size_t>(i)] = std::exp(i == ng - 1 ? l1 : l0 + i * by);
      const double v = fn(fs[static_cast<size_t>(i)]);
      if (v > bv) {
        bv = v;
        best = i;
      }
    }
    const double f = bmin([&](double u) { return -fn(u); }, fs[static_cast<size_t>(best > 0 ? best - 1 : 0)],
                          fs[static_cast<size_t>(best < ng - 1 ? best + 1 : ng - 1)], 1e-10);
    const double ms = mass(f * e0), cv = zs / (4 * (f * f)), pd1 = pod1(cv, model, m);
    cum += zs;
    double *o = out + static_cast<size_t>(s) * 8;
    o[0] = cum;
    o[1] = f;
    o[2] = ms / left;
    o[3] = cv;
    o[4] = pd1;
    o[5] = ms * pd1 / left;
    o[6] = 1 - (left - ms * pd1);
    o[7] = zs;
    left -= ms * pd1;
    a.push_back(f * e0);
    pd.push_back(pd1);
  }
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

}  // extern "C"

namespace {

// The figures of Frost (1996): number, caption, how drawn, and one line of
// R that computes a version of it from the reader's own numbers.
const char *const figs[][4] = {
    {"3-1", "Probability density distribution for a point datum with top view", "model",
     "plot(datumpmap(seq(-4, 4, by = 0.25), seq(-4, 4, by = 0.25), center = c(0, 0), sigma = 1.5), view = \"both\")"},
    {"3-2", "Circle of radius one standard deviation and its circumscribed square", "model",
     "sum(datumpmap(c(-1.5, 1.5), c(-1.5, 1.5), sigma = 1.5)$poc)"},
    {"3-3", "Circle of radius E and its circumscribed square", "model",
     "sum(datumpmap(c(-1, 1) * 1.177 * 1.5, c(-1, 1) * 1.177 * 1.5, sigma = 1.5)$poc)"},
    {"3-4", "Probability map for a point datum", "model",
     "plot(datumpmap(c(-3, -1, 1, 3) * 1.77, c(-3, -1, 1, 3) * 1.77, sigma = 1.5), view = \"cells\")"},
    {"3-5", "Circle of radius 1.1 E and its circumscribed square", "model",
     "sum(datumpmap(c(-1, 1) * 1.295 * 1.5, c(-1, 1) * 1.295 * 1.5, sigma = 1.5)$poc)"},
    {"3-6", "Line datum: circles at the end points boxed in", "schematic",
     "plot(datumpmap(seq(-4, 4, by = 0.25), seq(-4, 4, by = 0.25), \"line\", center = c(0, -2), end = c(0, 2), "
     "sigma = 1), view = \"both\")"},
    {"3-7", "Probability density distribution for a line datum connecting two point datums", "approximation",
     "plot(datumpmap(seq(-4, 4, by = 0.25), seq(-4, 4, by = 0.25), \"line\", center = c(0, -2), end = c(0, 2), "
     "sigma = 1), view = \"both\")"},
    {"3-8", "Probability density distribution for a line datum with top view", "model",
     "plot(datumpmap(seq(-4, 4, by = 0.25), seq(-4, 4, by = 0.25), \"corridor\", center = c(0, -2), end = c(0, 2), "
     "sigma = 1), view = \"both\")"},
    {"3-9", "Strip of width 2 E centred on the datum line", "model",
     "plot(datumpmap(c(-5, 5), c(-1.18, 1.18), \"corridor\", center = c(-4, 0), end = c(4, 0), sigma = 1), "
     "view = \"cells\")"},
    {"3-10", "Example of a completed probability map for a line datum", "model",
     "plot(datumpmap(0:10, c(-3, -1.8, -0.6, 0.6, 1.8, 3), \"corridor\", center = c(0, 0), end = c(10, 0), "
     "sigma = 1), view = \"cells\")"},
    {"3-11", "Uniform probability density distribution with top view", "model",
     "plot(datumpmap(seq(-4, 4, by = 0.25), seq(-4, 4, by = 0.25), \"area\", center = c(-2, -1.5), end = c(2, 1.5)), "
     "view = \"both\")"},
    {"3-12", "Probability map for an area datum", "model",
     "plot(datumpmap(0:5, 0:4, \"area\", center = c(0, 0), end = c(5, 4)), view = \"cells\")"},
    {"3-13", "Generalized probability density distribution with top view", "approximation",
     "plot(datumpmap(0:5, 0:4, \"general\", wts = matrix(c(1:10, 10:1), 4, 5)), view = \"both\")"},
    {"3-14", "Probability map for a generalized datum", "model",
     "plot(datumpmap(0:5, 0:4, \"general\", wts = matrix(c(1:10, 10:1), 4, 5)), view = \"cells\")"},
    {"4-1", "Lateral range curves: definite range, M-Beta, inverse cube and another sensor", "model",
     "plot(ltrlrng(seq(-3, 3, by = 0.05), W = 0.8, model = \"mbeta\", m = 0.5))"},
    {"4-2", "Definite range lateral range curve", "model",
     "plot(ltrlrng(seq(-3, 3, by = 0.05), W = 0.8, model = \"definite\"))"},
    {"4-3", "The area swept by a definite range sensor", "schematic",
     "plot(ltrlrng(seq(-3, 3, by = 0.05), W = 0.8, model = \"definite\"), equiv = TRUE)"},
    {"4-4", "M-Beta lateral range curve, M = 0.5", "model",
     "plot(ltrlrng(seq(-3, 3, by = 0.05), W = 0.8, model = \"mbeta\", m = 0.5), equiv = TRUE)"},
    {"4-5", "M-Beta lateral range curve, M = 0.25", "model",
     "plot(ltrlrng(seq(-3, 3, by = 0.05), W = 0.8, model = \"mbeta\", m = 0.25), equiv = TRUE)"},
    {"4-6", "The geometry of the inverse cube law", "schematic",
     "plot(ltrlrng(seq(-3, 3, by = 0.05), W = 0.8, model = \"invcube\"))"},
    {"4-7", "Inverse cube lateral range curve and its sweep width", "model",
     "plot(ltrlrng(seq(-3, 3, by = 0.05), W = 0.8, model = \"invcube\"), equiv = TRUE)"},
    {"5-1", "Random search", "model", "plot(srchpod(seq(0, 3, by = 0.1), model = \"random\"))"},
    {"5-2", "Parallel sweep search", "schematic", "pswppod(W = 0.8, S = 1, model = \"invcube\")"},
    {"5-3", "Definite range curves of two tracks 2 W apart", "model",
     "plot(ltrlrng(seq(-3, 3, by = 0.05), W = 0.8, model = \"definite\"), spacing = 1.6)"},
    {"5-4", "Definite range curves of two tracks W apart", "model",
     "plot(ltrlrng(seq(-3, 3, by = 0.05), W = 0.8, model = \"definite\"), spacing = 0.8)"},
    {"5-5", "Definite range and random search", "model",
     "plot(srchpod(seq(0, 3, by = 0.1), \"definite\"), models = c(\"definite\", \"random\"))"},
    {"5-6", "M-Beta (M = 0.5) curves of two tracks 2 W apart", "model",
     "plot(ltrlrng(seq(-3, 3, by = 0.05), W = 0.8, model = \"mbeta\", m = 0.5), spacing = 1.6)"},
    {"5-7", "M-Beta (M = 0.25 and M = 0.5)", "model",
     "plot(srchpod(seq(0, 3, by = 0.1), \"mbeta\", m = 0.5), models = c(\"mbeta\", \"definite\", \"random\"), "
     "m = c(0.5, 0.25))"},
    {"5-8", "Inverse cube", "model",
     "plot(srchpod(seq(0, 3, by = 0.1), \"invcube\"), models = c(\"invcube\", \"definite\", \"random\"))"},
    {"5-9", "Lateral range curves under ideal and poor search conditions", "model",
     "plot(ltrlrng(seq(-3, 3, by = 0.05), W = 0.8, model = \"invcube\", navsd = 0.3), spacing = 1)"},
    {"6-1", "A simple non-uniform distribution", "model", "optsqr(effort = 2)"},
    {"6-2", "The trial search square of half the area", "model", "optsqr(effort = 2)$side"},
    {"6-3", "Optimal search factors for point datums (effort to 18)", "model",
     "with(do.call(rbind, lapply(seq(0.5, 18, by = 0.5), optsqr)), plot(factor, cumeff, type = \"l\"))"},
    {"6-4", "Optimal search factors for point datums (effort to 100)", "model",
     "with(do.call(rbind, lapply(seq(2, 100, by = 2), optsqr)), plot(factor, cumeff, type = \"l\"))"},
    {"6-5", "Cumulative probability of success for optimal searches", "model",
     "with(do.call(rbind, lapply(seq(0.5, 25, by = 0.5), optsqr)), plot(cumeff, 100 * cumpos, type = \"l\"))"},
    {"6-6", "Optimal search factors for 1-5 sequential searches with random relative efforts", "approximation",
     "optsqr(c(1, 0.5, 2, 1.5))"},
};

}  // namespace

extern "C" int clab_sarfig(int i, const char **out) CLAB_TRY {
  const int n = static_cast<int>(sizeof figs / sizeof figs[0]);
  if (i < 0 || out == nullptr) return n;
  if (i >= n) return -1;
  for (int k = 0; k < 4; ++k) out[k] = figs[i][k];
  return n;
} CLAB_CATCH(CLAB_ENOMEM)
