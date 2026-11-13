// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, spatial econometrics: global and bivariate Moran's I with
// their inference, the spatial lag, error and combined models by maximum
// likelihood (and the lag model by two-stage least squares), and the
// direct, indirect and total effects of a lag model.
#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

namespace {

using vd = std::vector<double>;

// the upper normal tail
double upper(double z) { return 0.5 * std::erfc(z / std::sqrt(2.0)); }

// valid weights: rows of distinct neighbours, none the unit itself, finite
// weights from 0
bool wok(int n, const int *rp, const int *ci, const double *w) {
  if (rp[0] != 0) return false;
  for (int i = 0; i < n; ++i) {
    if (rp[i + 1] < rp[i]) return false;
    for (int q = rp[i]; q < rp[i + 1]; ++q) {
      if (ci[q] < 0 || ci[q] >= n || ci[q] == i || !kutil::nonneg(w[q])) return false;
    }
  }
  return true;
}

// S0, S1, S2 of Cliff and Ord (1981): S1 = sum_ij (w_ij + w_ji)^2 / 2 over
// ordered pairs, the reverse weight found by binary search in the sorted
// pair keys
void sums(int n, const int *rp, const int *ci, const double *w, double *s) {
  std::vector<std::pair<long long, double>> e;
  vd rs(n, 0.0), cs(n, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int q = rp[i]; q < rp[i + 1]; ++q) {
      e.emplace_back(static_cast<long long>(i) * n + ci[q], w[q]);
      rs[i] += w[q];
      cs[ci[q]] += w[q];
    }
  }
  std::sort(e.begin(), e.end());
  auto find = [&](long long key) {
    auto it = std::lower_bound(e.begin(), e.end(), std::make_pair(key, -HUGE_VAL));
    return it != e.end() && it->first == key ? it->second : -1.0;
  };
  kutil::ksum s0, s1, s2;
  for (int i = 0; i < n; ++i) {
    s0.add(rs[i]);
    s2.add((rs[i] + cs[i]) * (rs[i] + cs[i]));
  }
  for (const auto &v : e) {
    const long long i = v.first / n, j = v.first % n;
    const double r = find(j * n + i);
    const double a = v.second + (r >= 0 ? r : 0.0);
    // an ordered pair whose reverse is absent stands for both orders
    s1.add(r >= 0 ? a * a : 2 * a * a);
  }
  s[0] = s0.val();
  s[1] = 0.5 * s1.val();
  s[2] = s2.val();
}

// sum_ij w_ij a_i b_j
double cross(int n, const int *rp, const int *ci, const double *w, const double *a,
             const double *b) {
  kutil::ksum s;
  for (int i = 0; i < n; ++i) {
    double l = 0;
    for (int q = rp[i]; q < rp[i + 1]; ++q) l += w[q] * b[ci[q]];
    s.add(a[i] * l);
  }
  return s.val();
}

// ------------------------------------------------------- dense algebra

// LU with partial pivoting of the n x n matrix a (row-major, in place);
// piv: n ints. The log of |det| into *ld; false when singular (a pivot
// within n ulps of the largest entry).
bool lu(double *a, int n, int *piv, double *ld) {
  double l = 0, big = 0;
  for (std::size_t q = 0; q < static_cast<std::size_t>(n) * n; ++q)
    big = std::max(big, std::fabs(a[q]));
  // a pivot at rounding level of the largest entry counts as zero
  const double tiny = big * n * 2.220446049250313e-16;
  for (int k = 0; k < n; ++k) {
    int p = k;
    for (int i = k + 1; i < n; ++i) {
      if (std::fabs(a[static_cast<std::size_t>(i) * n + k]) >
          std::fabs(a[static_cast<std::size_t>(p) * n + k]))
        p = i;
    }
    piv[k] = p;
    if (p != k) {
      for (int j = 0; j < n; ++j) {
        std::swap(a[static_cast<std::size_t>(k) * n + j], a[static_cast<std::size_t>(p) * n + j]);
      }
    }
    const double d = a[static_cast<std::size_t>(k) * n + k];
    if (!(std::fabs(d) > tiny)) return false;
    l += std::log(std::fabs(d));
    for (int i = k + 1; i < n; ++i) {
      double *ri = a + static_cast<std::size_t>(i) * n;
      const double f = ri[k] / d;
      ri[k] = f;
      if (f == 0) continue;
      const double *rk = a + static_cast<std::size_t>(k) * n;
      for (int j = k + 1; j < n; ++j) ri[j] -= f * rk[j];
    }
  }
  *ld = l;
  return true;
}

// solve with the factors of lu() for m right-hand sides b (n x m, row-major)
void lus(const double *a, int n, const int *piv, double *b, int m) {
  for (int k = 0; k < n; ++k) {
    if (piv[k] != k) {
      for (int j = 0; j < m; ++j) {
        std::swap(b[static_cast<std::size_t>(k) * m + j],
                  b[static_cast<std::size_t>(piv[k]) * m + j]);
      }
    }
  }
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < i; ++k) {
      const double f = a[static_cast<std::size_t>(i) * n + k];
      if (f == 0) continue;
      for (int j = 0; j < m; ++j)
        b[static_cast<std::size_t>(i) * m + j] -= f * b[static_cast<std::size_t>(k) * m + j];
    }
  }
  for (int i = n - 1; i >= 0; --i) {
    for (int k = i + 1; k < n; ++k) {
      const double f = a[static_cast<std::size_t>(i) * n + k];
      if (f == 0) continue;
      for (int j = 0; j < m; ++j)
        b[static_cast<std::size_t>(i) * m + j] -= f * b[static_cast<std::size_t>(k) * m + j];
    }
    const double d = a[static_cast<std::size_t>(i) * n + i];
    for (int j = 0; j < m; ++j) b[static_cast<std::size_t>(i) * m + j] /= d;
  }
}

// the inverse of the symmetric positive definite p x p matrix a, in place
// (Gauss-Jordan through lu); false when singular
bool inv(double *a, int p) {
  vd f(a, a + static_cast<std::size_t>(p) * p), b(static_cast<std::size_t>(p) * p, 0.0);
  std::vector<int> piv(p);
  double ld;
  if (!lu(f.data(), p, piv.data(), &ld)) return false;
  for (int i = 0; i < p; ++i) b[static_cast<std::size_t>(i) * p + i] = 1;
  lus(f.data(), p, piv.data(), b.data(), p);
  std::copy(b.begin(), b.end(), a);
  return true;
}

// least squares of the columns y (n x m) on x (n x p, column-major) by
// Householder QR: coefficients into b (p x m, column-major); false when x
// has rank below p
bool ls(const double *x, int n, int p, const double *y, int m, double *b) {
  vd a(x, x + static_cast<std::size_t>(n) * p), r(y, y + static_cast<std::size_t>(n) * m);
  double big = 0;
  for (std::size_t q = 0; q < a.size(); ++q) big = std::max(big, std::fabs(a[q]));
  for (int k = 0; k < p; ++k) {
    double *c = a.data() + static_cast<std::size_t>(k) * n;
    double nr = 0;
    for (int i = k; i < n; ++i) nr += c[i] * c[i];
    nr = std::sqrt(nr);
    if (!(nr > 1e-11 * big * std::sqrt(static_cast<double>(n)))) return false;
    const double al = c[k] > 0 ? -nr : nr;
    c[k] -= al;
    double vv = 0;
    for (int i = k; i < n; ++i) vv += c[i] * c[i];
    auto refl = [&](double *v) {
      double s = 0;
      for (int i = k; i < n; ++i) s += c[i] * v[i];
      s = 2 * s / vv;
      for (int i = k; i < n; ++i) v[i] -= s * c[i];
    };
    for (int j = k + 1; j < p; ++j) refl(a.data() + static_cast<std::size_t>(j) * n);
    for (int j = 0; j < m; ++j) refl(r.data() + static_cast<std::size_t>(j) * n);
    c[k] = al;  // the diagonal of R (the reflector is no longer needed)
  }
  for (int j = 0; j < m; ++j) {
    double *bj = b + static_cast<std::size_t>(j) * p;
    const double *rj = r.data() + static_cast<std::size_t>(j) * n;
    for (int i = p - 1; i >= 0; --i) {
      double s = rj[i];
      for (int k = i + 1; k < p; ++k) s -= a[static_cast<std::size_t>(k) * n + i] * bj[k];
      bj[i] = s / a[static_cast<std::size_t>(i) * n + i];
    }
  }
  return true;
}

// Brent's minimiser of f on [lo, hi] (Brent 1973, ch. 5), to tol
template <class F>
double fmin(F f, double lo, double hi, double tol) {
  const double c = 0.5 * (3 - std::sqrt(5.0));
  double a = lo, b = hi, v = a + c * (b - a), w = v, x = v, e = 0, d = 0;
  double fx = f(x), fv = fx, fw = fx;
  for (int it = 0; it < 500; ++it) {
    if (kutil::stop(1e5)) break;
    const double xm = 0.5 * (a + b), t1 = 1e-12 * std::fabs(x) + tol / 3, t2 = 2 * t1;
    if (std::fabs(x - xm) <= t2 - 0.5 * (b - a)) break;
    double p = 0, q = 0, r = 0;
    if (std::fabs(e) > t1) {
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = 2 * (q - r);
      if (q > 0) p = -p;
      q = std::fabs(q);
      r = e;
      e = d;
    }
    if (std::fabs(p) >= std::fabs(0.5 * q * r) || p <= q * (a - x) || p >= q * (b - x)) {
      e = x < xm ? b - x : a - x;
      d = c * e;
    } else {
      d = p / q;
      const double u = x + d;
      if (u - a < t2 || b - u < t2) d = x < xm ? t1 : -t1;
    }
    const double u = x + (std::fabs(d) >= t1 ? d : (d > 0 ? t1 : -t1));
    const double fu = f(u);
    if (fu <= fx) {
      if (u < x)
        b = x;
      else
        a = x;
      v = w;
      fv = fw;
      w = x;
      fw = fx;
      x = u;
      fx = fu;
    } else {
      if (u < x)
        a = u;
      else
        b = u;
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

// the model's pieces: dense W (row-major), X (column-major), y
struct mdl {
  int n, p;
  vd W, X, y, Wy;
  // log |I - a W|, or -inf when singular
  double ldet(double a) const {
    vd m(static_cast<std::size_t>(n) * n);
    for (std::size_t q = 0; q < m.size(); ++q) m[q] = -a * W[q];
    for (int i = 0; i < n; ++i) m[static_cast<std::size_t>(i) * n + i] += 1;
    std::vector<int> piv(n);
    double l;
    return lu(m.data(), n, piv.data(), &l) ? l : -HUGE_VAL;
  }
  // (I - a W) v
  void filt(double a, const double *v, double *o) const {
    for (int i = 0; i < n; ++i) {
      double s = 0;
      const double *r = W.data() + static_cast<std::size_t>(i) * n;
      for (int j = 0; j < n; ++j) s += r[j] * v[j];
      o[i] = v[i] - a * s;
    }
  }
  // for the lag rho and the error lambda: beta (p) and the residual sum of
  // squares of B (A y - X beta), A = I - rho W, B = I - lam W; false when
  // the filtered X loses rank
  bool fit(double rho, double lam, double *beta, double *sse) const {
    vd ay(n), by(n), bx(static_cast<std::size_t>(n) * p);
    filt(rho, y.data(), ay.data());
    filt(lam, ay.data(), by.data());
    for (int j = 0; j < p; ++j) {
      filt(lam, X.data() + static_cast<std::size_t>(j) * n,
           bx.data() + static_cast<std::size_t>(j) * n);
    }
    if (!ls(bx.data(), n, p, by.data(), 1, beta)) return false;
    kutil::ksum s;
    for (int i = 0; i < n; ++i) {
      double f = by[i];
      for (int j = 0; j < p; ++j) f -= bx[static_cast<std::size_t>(j) * n + i] * beta[j];
      s.add(f * f);
    }
    *sse = s.val();
    return true;
  }
  // the full log-likelihood at (beta, rho, lam, s2)
  double ll(const double *beta, double rho, double lam, double s2) const {
    vd ay(n), r(n), br(n);
    filt(rho, y.data(), ay.data());
    for (int i = 0; i < n; ++i) {
      double f = ay[i];
      for (int j = 0; j < p; ++j) f -= X[static_cast<std::size_t>(j) * n + i] * beta[j];
      r[i] = f;
    }
    filt(lam, r.data(), br.data());
    double s = 0;
    for (int i = 0; i < n; ++i) s += br[i] * br[i];
    return -0.5 * n * std::log(2 * 3.14159265358979323846 * s2) - s / (2 * s2) + ldet(rho) +
           ldet(lam);
  }
};

// S = (I - a W)^-1 (row-major); false when singular
bool sinv(const vd &W, int n, double a, vd &S) {
  vd m(static_cast<std::size_t>(n) * n);
  for (std::size_t q = 0; q < m.size(); ++q) m[q] = -a * W[q];
  for (int i = 0; i < n; ++i) m[static_cast<std::size_t>(i) * n + i] += 1;
  std::vector<int> piv(n);
  double l;
  if (!lu(m.data(), n, piv.data(), &l)) return false;
  S.assign(static_cast<std::size_t>(n) * n, 0.0);
  for (int i = 0; i < n; ++i) S[static_cast<std::size_t>(i) * n + i] = 1;
  lus(m.data(), n, piv.data(), S.data(), n);
  return true;
}

// the product of two n x n row-major matrices
vd mul(const vd &a, const vd &b, int n) {
  vd c(static_cast<std::size_t>(n) * n, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < n; ++k) {
      const double f = a[static_cast<std::size_t>(i) * n + k];
      if (f == 0) continue;
      const double *rb = b.data() + static_cast<std::size_t>(k) * n;
      double *rc = c.data() + static_cast<std::size_t>(i) * n;
      for (int j = 0; j < n; ++j) rc[j] += f * rb[j];
    }
  }
  return c;
}

double tr(const vd &a, int n) {
  double s = 0;
  for (int i = 0; i < n; ++i) s += a[static_cast<std::size_t>(i) * n + i];
  return s;
}

// tr(a' a) = sum of squares, and tr(a a)
double trtt(const vd &a, int n) {
  double s = 0;
  for (double v : a) s += v * v;
  (void)n;
  return s;
}
double trsq(const vd &a, int n) {
  double s = 0;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      s += a[static_cast<std::size_t>(i) * n + j] * a[static_cast<std::size_t>(j) * n + i];
  return s;
}

vd dense(int n, const int *rp, const int *ci, const double *w) {
  vd W(static_cast<std::size_t>(n) * n, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int q = rp[i]; q < rp[i + 1]; ++q) W[static_cast<std::size_t>(i) * n + ci[q]] = w[q];
  }
  return W;
}

}  // namespace

int clab_gmoran(const double *x, int n, int k, const int *rp, const int *ci, const double *w,
                int nsim, unsigned int seed, double *out) {
  if (n < 4 || k < 1 || nsim < 0) return CLAB_EARG;
  if (!wok(n, rp, ci, w)) return CLAB_EDOM;
  double s[3];
  sums(n, rp, ci, w, s);
  if (!(s[0] > 0)) return CLAB_EDOM;
  const double N = n, e = -1.0 / (n - 1), s0 = s[0], s1 = s[1], s2 = s[2];
  vd z(n), zp(n);
  for (int c = 0; c < k; ++c) {
    const double *v = x + static_cast<std::size_t>(c) * n;
    double m = 0;
    for (int i = 0; i < n; ++i) {
      if (!std::isfinite(v[i])) return CLAB_EDOM;
      m += v[i];
    }
    m /= n;
    double m2 = 0, m4 = 0;
    for (int i = 0; i < n; ++i) {
      z[i] = v[i] - m;
      m2 += z[i] * z[i];
      m4 += z[i] * z[i] * z[i] * z[i];
    }
    double *o = out + static_cast<std::size_t>(c) * 10;
    if (!(m2 > 0)) {
      for (int q = 0; q < 10; ++q) o[q] = kutil::nan();
      continue;
    }
    const double I = N / s0 * cross(n, rp, ci, w, z.data(), z.data()) / m2;
    const double b2 = N * m4 / (m2 * m2);
    const double vn = (N * N * s1 - N * s2 + 3 * s0 * s0) / ((N * N - 1) * s0 * s0) - e * e;
    const double vr = (N * ((N * N - 3 * N + 3) * s1 - N * s2 + 3 * s0 * s0) -
                       b2 * ((N * N - N) * s1 - 2 * N * s2 + 6 * s0 * s0)) /
                          ((N - 1) * (N - 2) * (N - 3) * s0 * s0) -
                      e * e;
    o[0] = I;
    o[1] = e;
    o[2] = vn;
    o[3] = vr;
    o[4] = (I - e) / std::sqrt(vn);
    o[5] = vr > 0 ? (I - e) / std::sqrt(vr) : kutil::nan();
    o[6] = upper(o[4]);
    o[7] = vr > 0 ? upper(o[5]) : kutil::nan();
    o[9] = s0;
    o[8] = kutil::nan();
    if (nsim > 0) {
      kutil::rng g(seed + static_cast<unsigned int>(c));
      zp = z;
      double ge = 0;
      for (int r = 0; r < nsim; ++r) {
        if (kutil::stop(n + rp[n])) return CLAB_EINTR;
        for (int i = n - 1; i > 0; --i) std::swap(zp[i], zp[g.below(i + 1)]);
        const double Ir = N / s0 * cross(n, rp, ci, w, zp.data(), zp.data()) / m2;
        ge += Ir >= I;
      }
      o[8] = (1 + ge) / (nsim + 1.0);
    }
  }
  return CLAB_OK;
}

int clab_bvmor(const double *x, const double *y, int n, const int *rp, const int *ci,
               const double *w, int nsim, unsigned int seed, double *out) {
  if (n < 3 || nsim < 0) return CLAB_EARG;
  if (!wok(n, rp, ci, w)) return CLAB_EDOM;
  double s[3];
  sums(n, rp, ci, w, s);
  double mx = 0, my = 0;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(x[i]) || !std::isfinite(y[i])) return CLAB_EDOM;
    mx += x[i];
    my += y[i];
  }
  mx /= n;
  my /= n;
  vd zx(n), zy(n);
  double sx = 0, sy = 0;
  for (int i = 0; i < n; ++i) {
    zx[i] = x[i] - mx;
    zy[i] = y[i] - my;
    sx += zx[i] * zx[i];
    sy += zy[i] * zy[i];
  }
  if (!(s[0] > 0) || !(sx > 0) || !(sy > 0)) return CLAB_EDOM;
  const double f = n / s[0] / std::sqrt(sx * sy);
  const double I = f * cross(n, rp, ci, w, zx.data(), zy.data());
  out[0] = I;
  for (int q = 1; q < 6; ++q) out[q] = kutil::nan();
  if (nsim > 0) {
    // y's values shuffled over the units, x held
    kutil::rng g(seed);
    double sm = 0, ss = 0, ge = 0, le = 0;
    for (int r = 0; r < nsim; ++r) {
      if (kutil::stop(n + rp[n])) return CLAB_EINTR;
      for (int i = n - 1; i > 0; --i) std::swap(zy[i], zy[g.below(i + 1)]);
      const double v = f * cross(n, rp, ci, w, zx.data(), zy.data());
      sm += v;
      ss += v * v;
      ge += v >= I;
      le += v <= I;
    }
    const double mean = sm / nsim;
    const double sd =
        nsim > 1 ? std::sqrt(std::max(0.0, (ss - nsim * mean * mean) / (nsim - 1))) : kutil::nan();
    out[1] = mean;
    out[2] = sd;
    out[3] = sd > 0 ? (I - mean) / sd : kutil::nan();
    out[4] = (1 + ge) / (nsim + 1.0);
    out[5] = (1 + std::min(ge, le)) / (nsim + 1.0);
  }
  return CLAB_OK;
}

int clab_sarml(const double *y, const double *x, int n, int p, const int *rp, const int *ci,
               const double *w, int model, int method, const double *inst, int ni, double *coef,
               double *vc, double *st) {
  if (n < 3 || p < 1 || p >= n || model < 0 || model > 2 || method < 0 || method > 1) {
    return CLAB_EARG;
  }
  if (method == 1 && (model != 0 || inst == nullptr || ni < 1)) return CLAB_EARG;
  if (!wok(n, rp, ci, w)) return CLAB_EDOM;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(y[i])) return CLAB_EDOM;
  }
  for (long q = 0; q < static_cast<long>(n) * p; ++q) {
    if (!std::isfinite(x[q])) return CLAB_EDOM;
  }
  mdl m;
  m.n = n;
  m.p = p;
  m.W = dense(n, rp, ci, w);
  m.X.assign(x, x + static_cast<std::size_t>(n) * p);
  m.y.assign(y, y + n);
  m.Wy.assign(n, 0.0);
  double rmax = 0;
  for (int i = 0; i < n; ++i) {
    double rs = 0;
    for (int q = rp[i]; q < rp[i + 1]; ++q) {
      m.Wy[i] += w[q] * y[ci[q]];
      rs += w[q];
    }
    rmax = std::max(rmax, rs);
  }
  if (!(rmax > 0)) return CLAB_EDOM;
  // OLS, the comparison: coefficients and their usual standard errors
  vd bo(p), xtx(static_cast<std::size_t>(p) * p);
  if (!ls(x, n, p, y, 1, bo.data())) return CLAB_EDOM;
  double sso = 0;
  for (int i = 0; i < n; ++i) {
    double f = y[i];
    for (int j = 0; j < p; ++j) f -= x[static_cast<std::size_t>(j) * n + i] * bo[j];
    sso += f * f;
  }
  for (int a = 0; a < p; ++a)
    for (int b = 0; b < p; ++b) {
      double s = 0;
      for (int i = 0; i < n; ++i)
        s += x[static_cast<std::size_t>(a) * n + i] * x[static_cast<std::size_t>(b) * n + i];
      xtx[static_cast<std::size_t>(a) * p + b] = s;
    }
  if (!inv(xtx.data(), p)) return CLAB_EDOM;
  for (int j = 0; j < p; ++j) {
    coef[2 * p + 4 + j] = bo[j];
    coef[3 * p + 4 + j] = std::sqrt(sso / (n - p) * xtx[static_cast<std::size_t>(j) * p + j]);
  }
  // the parameters: beta (p), rho, lambda; vc is (p + 2) square, row-major
  const int k = p + 2;
  std::fill(vc, vc + static_cast<std::size_t>(k) * k, kutil::nan());
  vd beta(p);
  double rho = 0, lam = 0, s2 = 0;
  // the search interval: the spectral radius of W is at most its largest
  // row sum
  const double lo = -1 / rmax + 1e-9, hi = 1 / rmax - 1e-9;
  const double tol = 1e-10;
  if (method == 1) {
    // two-stage least squares with the instruments [X, inst]
    const int q = p + ni;
    vd H(static_cast<std::size_t>(n) * q), Z(static_cast<std::size_t>(n) * (p + 1)),
        G(static_cast<std::size_t>(q) * (p + 1)), Zh(static_cast<std::size_t>(n) * (p + 1));
    std::copy(x, x + static_cast<std::size_t>(n) * p, H.begin());
    std::copy(inst, inst + static_cast<std::size_t>(n) * ni,
              H.begin() + static_cast<std::size_t>(n) * p);
    std::copy(m.Wy.begin(), m.Wy.end(), Z.begin());
    std::copy(x, x + static_cast<std::size_t>(n) * p, Z.begin() + n);
    if (!ls(H.data(), n, q, Z.data(), p + 1, G.data())) return CLAB_EDOM;
    for (int j = 0; j <= p; ++j)
      for (int i = 0; i < n; ++i) {
        double s = 0;
        for (int a = 0; a < q; ++a)
          s += H[static_cast<std::size_t>(a) * n + i] * G[static_cast<std::size_t>(j) * q + a];
        Zh[static_cast<std::size_t>(j) * n + i] = s;
      }
    vd d(p + 1);
    if (!ls(Zh.data(), n, p + 1, y, 1, d.data())) return CLAB_EDOM;
    rho = d[0];
    for (int j = 0; j < p; ++j) beta[j] = d[j + 1];
    double ss = 0;
    for (int i = 0; i < n; ++i) {
      double f = y[i] - rho * m.Wy[i];
      for (int j = 0; j < p; ++j) f -= x[static_cast<std::size_t>(j) * n + i] * beta[j];
      ss += f * f;
    }
    s2 = ss / n;
    vd zz(static_cast<std::size_t>(p + 1) * (p + 1));
    for (int a = 0; a <= p; ++a)
      for (int b = 0; b <= p; ++b) {
        double s = 0;
        for (int i = 0; i < n; ++i)
          s += Zh[static_cast<std::size_t>(a) * n + i] * Zh[static_cast<std::size_t>(b) * n + i];
        zz[static_cast<std::size_t>(a) * (p + 1) + b] = s;
      }
    if (!inv(zz.data(), p + 1)) return CLAB_EDOM;
    // order (beta, rho)
    for (int a = 0; a <= p; ++a)
      for (int b = 0; b <= p; ++b) {
        const int ia = a == 0 ? p : a - 1, ib = b == 0 ? p : b - 1;
        vc[static_cast<std::size_t>(ia) * k + ib] =
            ss / (n - p - 1) * zz[static_cast<std::size_t>(a) * (p + 1) + b];
      }
    st[0] = kutil::nan();
  } else if (model == 0 || model == 1) {
    const bool lag = model == 0;
    auto nll = [&](double a) {
      double sse;
      vd b(p);
      if (!(lag ? m.fit(a, 0, b.data(), &sse) : m.fit(0, a, b.data(), &sse))) return HUGE_VAL;
      return 0.5 * n * std::log(sse / n) - m.ldet(a);
    };
    const double a = fmin(nll, lo, hi, tol);
    double sse;
    if (lag) {
      rho = a;
      m.fit(rho, 0, beta.data(), &sse);
    } else {
      lam = a;
      m.fit(0, lam, beta.data(), &sse);
    }
    s2 = sse / n;
    // the information matrix (Anselin 1988, ch. 6) in (beta, a, s2)
    vd S;
    if (!sinv(m.W, n, a, S)) return CLAB_EDOM;
    const vd WA = mul(m.W, S, n);
    const double t1 = tr(WA, n), t2 = trsq(WA, n), t3 = trtt(WA, n);
    const int kk = p + 2;
    vd I(static_cast<std::size_t>(kk) * kk, 0.0);
    if (lag) {
      // X'X / s2, X' WA X beta / s2, and for rho the traces plus the
      // squared length of WA X beta
      vd xb(n, 0.0), wxb(n, 0.0);
      for (int i = 0; i < n; ++i)
        for (int j = 0; j < p; ++j) xb[i] += x[static_cast<std::size_t>(j) * n + i] * beta[j];
      for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) wxb[i] += WA[static_cast<std::size_t>(i) * n + j] * xb[j];
      for (int a1 = 0; a1 < p; ++a1) {
        for (int b1 = 0; b1 < p; ++b1) {
          double s = 0;
          for (int i = 0; i < n; ++i)
            s += x[static_cast<std::size_t>(a1) * n + i] * x[static_cast<std::size_t>(b1) * n + i];
          I[static_cast<std::size_t>(a1) * kk + b1] = s / s2;
        }
        double s = 0;
        for (int i = 0; i < n; ++i) s += x[static_cast<std::size_t>(a1) * n + i] * wxb[i];
        I[static_cast<std::size_t>(a1) * kk + p] = I[static_cast<std::size_t>(p) * kk + a1] =
            s / s2;
      }
      double q = 0;
      for (int i = 0; i < n; ++i) q += wxb[i] * wxb[i];
      I[static_cast<std::size_t>(p) * kk + p] = t2 + t3 + q / s2;
    } else {
      // beta from (I - lam W) X, separate from (lam, s2)
      vd bx(static_cast<std::size_t>(n) * p);
      for (int j = 0; j < p; ++j)
        m.filt(lam, x + static_cast<std::size_t>(j) * n,
               bx.data() + static_cast<std::size_t>(j) * n);
      for (int a1 = 0; a1 < p; ++a1)
        for (int b1 = 0; b1 < p; ++b1) {
          double s = 0;
          for (int i = 0; i < n; ++i)
            s +=
                bx[static_cast<std::size_t>(a1) * n + i] * bx[static_cast<std::size_t>(b1) * n + i];
          I[static_cast<std::size_t>(a1) * kk + b1] = s / s2;
        }
      I[static_cast<std::size_t>(p) * kk + p] = t2 + t3;
    }
    I[static_cast<std::size_t>(p) * kk + p + 1] = I[static_cast<std::size_t>(p + 1) * kk + p] =
        t1 / s2;
    I[static_cast<std::size_t>(p + 1) * kk + p + 1] = n / (2 * s2 * s2);
    if (!inv(I.data(), kk)) return CLAB_EDOM;
    // into (beta, rho, lambda)
    const int at = lag ? p : p + 1;
    for (int a1 = 0; a1 <= p; ++a1)
      for (int b1 = 0; b1 <= p; ++b1) {
        const int ia = a1 < p ? a1 : at, ib = b1 < p ? b1 : at;
        vc[static_cast<std::size_t>(ia) * k + ib] = I[static_cast<std::size_t>(a1) * kk + b1];
      }
  } else {
    // both: Nelder-Mead on (rho, lambda) inside the box, then the
    // covariance from a finite-difference Hessian of the full likelihood
    auto nll = [&](const double *v) {
      if (!(v[0] > lo && v[0] < hi && v[1] > lo && v[1] < hi)) return HUGE_VAL;
      double sse;
      vd b(p);
      if (!m.fit(v[0], v[1], b.data(), &sse)) return HUGE_VAL;
      return 0.5 * n * std::log(sse / n) - m.ldet(v[0]) - m.ldet(v[1]);
    };
    double sx[3][2] = {{0, 0}, {0.2 * hi, 0}, {0, 0.2 * hi}}, fv[3];
    for (int r = 0; r < 3; ++r) fv[r] = nll(sx[r]);
    for (int it = 0; it < 2000; ++it) {
      if (kutil::stop(static_cast<double>(n) * n)) break;
      int o[3] = {0, 1, 2};
      std::sort(o, o + 3, [&](int a1, int b1) { return fv[a1] < fv[b1]; });
      const int bst = o[0], mid = o[1], wst = o[2];
      if (std::fabs(sx[wst][0] - sx[bst][0]) + std::fabs(sx[wst][1] - sx[bst][1]) < 1e-11 &&
          std::fabs(fv[wst] - fv[bst]) < 1e-12 * (1 + std::fabs(fv[bst]))) {
        break;
      }
      double c[2], rf[2];
      for (int j = 0; j < 2; ++j) {
        c[j] = 0.5 * (sx[bst][j] + sx[mid][j]);
        rf[j] = c[j] + (c[j] - sx[wst][j]);
      }
      const double fr = nll(rf);
      if (fr < fv[bst]) {
        double ex[2] = {c[0] + 2 * (rf[0] - c[0]), c[1] + 2 * (rf[1] - c[1])};
        const double fe = nll(ex);
        const double *nw = fe < fr ? ex : rf;
        sx[wst][0] = nw[0];
        sx[wst][1] = nw[1];
        fv[wst] = std::min(fe, fr);
      } else if (fr < fv[mid]) {
        sx[wst][0] = rf[0];
        sx[wst][1] = rf[1];
        fv[wst] = fr;
      } else {
        double ct[2] = {c[0] + 0.5 * (sx[wst][0] - c[0]), c[1] + 0.5 * (sx[wst][1] - c[1])};
        const double fc = nll(ct);
        if (fc < fv[wst]) {
          sx[wst][0] = ct[0];
          sx[wst][1] = ct[1];
          fv[wst] = fc;
        } else {
          for (int r = 0; r < 3; ++r) {
            if (r == bst) continue;
            sx[r][0] = 0.5 * (sx[r][0] + sx[bst][0]);
            sx[r][1] = 0.5 * (sx[r][1] + sx[bst][1]);
            fv[r] = nll(sx[r]);
          }
        }
      }
    }
    int bst = 0;
    for (int r = 1; r < 3; ++r) {
      if (fv[r] < fv[bst]) bst = r;
    }
    rho = sx[bst][0];
    lam = sx[bst][1];
    double sse;
    if (!m.fit(rho, lam, beta.data(), &sse)) return CLAB_EDOM;
    s2 = sse / n;
    // the information matrix in (beta, rho, lambda, s2): with SA = A^-1,
    // SB = B^-1, WA = W SA, WB = W SB and C = B WA SB (the expectations of
    // the second derivatives of the log-likelihood)
    vd SA, SB;
    if (!sinv(m.W, n, rho, SA) || !sinv(m.W, n, lam, SB)) return CLAB_EDOM;
    const vd WA = mul(m.W, SA, n), WB = mul(m.W, SB, n);
    vd BWA = mul(m.W, WA, n);
    for (std::size_t q = 0; q < BWA.size(); ++q) BWA[q] = WA[q] - lam * BWA[q];
    const vd C = mul(BWA, SB, n);
    double cr = 0, cs = 0;
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {
        cr += WB[static_cast<std::size_t>(i) * n + j] * C[static_cast<std::size_t>(i) * n + j];
        cs += WB[static_cast<std::size_t>(i) * n + j] * C[static_cast<std::size_t>(j) * n + i];
      }
    vd xb(n, 0.0), g(n, 0.0), bx(static_cast<std::size_t>(n) * p);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < p; ++j) xb[i] += x[static_cast<std::size_t>(j) * n + i] * beta[j];
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) g[i] += BWA[static_cast<std::size_t>(i) * n + j] * xb[j];
    for (int j = 0; j < p; ++j)
      m.filt(lam, x + static_cast<std::size_t>(j) * n, bx.data() + static_cast<std::size_t>(j) * n);
    const int kk = p + 3;
    vd I(static_cast<std::size_t>(kk) * kk, 0.0);
    auto at = [&](int r, int c) -> double & { return I[static_cast<std::size_t>(r) * kk + c]; };
    for (int a1 = 0; a1 < p; ++a1) {
      for (int b1 = 0; b1 < p; ++b1) {
        double sm = 0;
        for (int i = 0; i < n; ++i)
          sm += bx[static_cast<std::size_t>(a1) * n + i] * bx[static_cast<std::size_t>(b1) * n + i];
        at(a1, b1) = sm / s2;
      }
      double sm = 0;
      for (int i = 0; i < n; ++i) sm += bx[static_cast<std::size_t>(a1) * n + i] * g[i];
      at(a1, p) = at(p, a1) = sm / s2;
    }
    double q = 0;
    for (int i = 0; i < n; ++i) q += g[i] * g[i];
    at(p, p) = trsq(WA, n) + trtt(C, n) + q / s2;
    at(p, p + 1) = at(p + 1, p) = cr + cs;
    at(p + 1, p + 1) = trsq(WB, n) + trtt(WB, n);
    at(p, p + 2) = at(p + 2, p) = tr(WA, n) / s2;
    at(p + 1, p + 2) = at(p + 2, p + 1) = tr(WB, n) / s2;
    at(p + 2, p + 2) = n / (2 * s2 * s2);
    if (!inv(I.data(), kk)) return CLAB_EDOM;
    for (int a1 = 0; a1 < k; ++a1)
      for (int b1 = 0; b1 < k; ++b1)
        vc[static_cast<std::size_t>(a1) * k + b1] = I[static_cast<std::size_t>(a1) * kk + b1];
  }
  for (int j = 0; j < p; ++j) coef[j] = beta[j];
  coef[p] = model == 1 ? kutil::nan() : rho;
  coef[p + 1] = model == 0 ? kutil::nan() : lam;
  for (int j = 0; j < k; ++j) {
    const double v = vc[static_cast<std::size_t>(j) * k + j];
    coef[k + j] = v >= 0 ? std::sqrt(v) : kutil::nan();
  }
  st[1] = s2;
  if (method == 0) {
    st[0] = m.ll(beta.data(), rho, lam, s2);
  }
  // parameters: beta, rho or lambda (or both), s2
  const int np = p + (model == 2 ? 2 : 1) + (method == 0 ? 1 : 0);
  st[2] = method == 0 ? -2 * st[0] + 2 * np : kutil::nan();
  st[3] = -0.5 * n * std::log(2 * 3.14159265358979323846 * sso / n) - 0.5 * n;
  st[4] = rmax;
  return CLAB_OK;
}

int clab_spill(int n, const int *rp, const int *ci, const double *w, double rho, const double *b,
               const double *th, int k, const double *vc, double *out) {
  if (n < 2 || k < 1) return CLAB_EARG;
  if (!wok(n, rp, ci, w) || !std::isfinite(rho)) return CLAB_EDOM;
  for (int j = 0; j < k; ++j) {
    if (!std::isfinite(b[j]) || (th != nullptr && !std::isfinite(th[j]))) return CLAB_EDOM;
  }
  const vd W = dense(n, rp, ci, w);
  vd S;
  if (!sinv(W, n, rho, S)) return CLAB_EDOM;
  const vd SW = mul(S, W, n), SWS = mul(SW, S, n), SWSW = mul(SWS, W, n);
  auto tot = [&](const vd &a) {
    kutil::ksum s;
    for (double v : a) s.add(v);
    return s.val() / n;
  };
  // per unit of beta and theta: the mean diagonal (direct) and the mean
  // row sum (total), and their derivatives in rho (dS / drho = S W S)
  const double d1 = tr(S, n) / n, d2 = tr(SW, n) / n, t1 = tot(S), t2 = tot(SW);
  const double dd1 = tr(SWS, n) / n, dd2 = tr(SWSW, n) / n, dt1 = tot(SWS), dt2 = tot(SWSW);
  const int m = th == nullptr ? 1 + k : 1 + 2 * k;
  for (int j = 0; j < k; ++j) {
    const double tj = th == nullptr ? 0.0 : th[j];
    const double dir = b[j] * d1 + tj * d2, all = b[j] * t1 + tj * t2;
    double *o = out + 6 * j;
    o[0] = dir;
    o[1] = all - dir;
    o[2] = all;
    o[3] = o[4] = o[5] = kutil::nan();
    if (vc == nullptr) continue;
    // the delta method in (rho, beta_j, theta_j)
    const double gd[3] = {b[j] * dd1 + tj * dd2, d1, d2};
    const double gt[3] = {b[j] * dt1 + tj * dt2, t1, t2};
    double gi[3];
    for (int q = 0; q < 3; ++q) gi[q] = gt[q] - gd[q];
    const int ix[3] = {0, 1 + j, th == nullptr ? -1 : 1 + k + j};
    auto qf = [&](const double *g) {
      double s = 0;
      for (int a = 0; a < 3; ++a)
        for (int c = 0; c < 3; ++c) {
          if (ix[a] < 0 || ix[c] < 0) continue;
          s += g[a] * g[c] * vc[static_cast<std::size_t>(ix[a]) * m + ix[c]];
        }
      return s;
    };
    o[3] = std::sqrt(std::max(0.0, qf(gd)));
    o[4] = std::sqrt(std::max(0.0, qf(gi)));
    o[5] = std::sqrt(std::max(0.0, qf(gt)));
  }
  return CLAB_OK;
}
