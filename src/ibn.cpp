// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, internal benchmarking of officers (Ridgeway and MacDonald
// 2009): each officer's stops against the other officers' stops made in
// the same circumstances, by propensity weights, with the doubly robust
// regression adjustment, and the local false discovery rate of Efron
// (2004) over the officers.
#include <algorithm>
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

namespace {

constexpr double ir2pi = 0.39894228040143267794;

double lgs(double x) { return 1.0 / (1.0 + std::exp(-x)); }
double pn(double z) { return 0.5 * std::erfc(-z / 1.41421356237309504880); }
double dn(double z) { return ir2pi * std::exp(-0.5 * z * z); }

// Cholesky of the k x k symmetric positive definite a (column-major) in
// place, lower; false when not positive definite.
bool chol(double *a, int k) {
  for (int j = 0; j < k; ++j) {
    double s = a[j * k + j];
    for (int m = 0; m < j; ++m) s -= a[m * k + j] * a[m * k + j];
    if (!(s > 0) || !std::isfinite(s)) return false;
    const double d = std::sqrt(s);
    a[j * k + j] = d;
    for (int i = j + 1; i < k; ++i) {
      double t = a[j * k + i];
      for (int m = 0; m < j; ++m) t -= a[m * k + i] * a[m * k + j];
      a[j * k + i] = t / d;
    }
  }
  return true;
}

// Solves L L' x = b in place, L from chol().
void csolve(const double *L, int k, double *b) {
  for (int i = 0; i < k; ++i) {
    double t = b[i];
    for (int m = 0; m < i; ++m) t -= L[m * k + i] * b[m];
    b[i] = t / L[i * k + i];
  }
  for (int i = k - 1; i >= 0; --i) {
    double t = b[i];
    for (int m = i + 1; m < k; ++m) t -= L[i * k + m] * b[m];
    b[i] = t / L[i * k + i];
  }
}

// The weighted logistic regression of y on the n x k design z
// (column-major) with weights w, by Newton steps from 0; coef (k); the
// Cholesky factor of the information at the solution in L (k x k).
// false when it does not converge in 100 steps or the information is
// singular.
bool wlogit(const double *z, const double *y, const double *w, int n, int k, double *coef, double *L) {
  std::vector<double> g(static_cast<std::size_t>(k)), eta(static_cast<std::size_t>(n));
  for (int j = 0; j < k; ++j) coef[j] = 0.0;
  for (int it = 0; it < 100; ++it) {
    for (std::size_t c = 0; c < static_cast<std::size_t>(k) * static_cast<std::size_t>(k); ++c) L[c] = 0.0;
    for (int j = 0; j < k; ++j) g[static_cast<std::size_t>(j)] = 0.0;
    for (int i = 0; i < n; ++i) {
      double e = 0;
      for (int j = 0; j < k; ++j) e += z[static_cast<std::size_t>(j) * n + i] * coef[j];
      const double mu = lgs(e), v = w[i] * mu * (1 - mu), r = w[i] * (y[i] - mu);
      for (int j = 0; j < k; ++j) {
        const double zj = z[static_cast<std::size_t>(j) * n + i];
        g[static_cast<std::size_t>(j)] += r * zj;
        for (int m = j; m < k; ++m) L[j * k + m] += v * zj * z[static_cast<std::size_t>(m) * n + i];
      }
    }
    if (!chol(L, k)) return false;
    std::vector<double> d(g);
    csolve(L, k, d.data());
    double mx = 0;
    for (int j = 0; j < k; ++j) {
      coef[j] += d[static_cast<std::size_t>(j)];
      mx = std::fmax(mx, std::fabs(d[static_cast<std::size_t>(j)]) / (1 + std::fabs(coef[j])));
      if (!std::isfinite(coef[j])) return false;
    }
    if (mx < 1e-10) {
      // the information at the solution
      for (std::size_t c = 0; c < static_cast<std::size_t>(k) * static_cast<std::size_t>(k); ++c) L[c] = 0.0;
      for (int i = 0; i < n; ++i) {
        double e = 0;
        for (int j = 0; j < k; ++j) e += z[static_cast<std::size_t>(j) * n + i] * coef[j];
        const double mu = lgs(e), v = w[i] * mu * (1 - mu);
        for (int j = 0; j < k; ++j) {
          for (int m = j; m < k; ++m) {
            L[j * k + m] += v * z[static_cast<std::size_t>(j) * n + i] * z[static_cast<std::size_t>(m) * n + i];
          }
        }
      }
      return chol(L, k);
    }
  }
  return false;
}

// The quantile of type 7 of sorted x (n >= 1).
double q7(const std::vector<double> &x, double p) {
  const double h = (static_cast<double>(x.size()) - 1) * p;
  const std::size_t lo = static_cast<std::size_t>(std::floor(h));
  const std::size_t hi = lo + 1 < x.size() ? lo + 1 : lo;
  return x[lo] + (h - static_cast<double>(lo)) * (x[hi] - x[lo]);
}

// The maximum likelihood fit of the null N(d, s) to the z's in xm -+ xw,
// as locmle() of the package locfdr: false when it fails.
bool locmle(const std::vector<double> &z, double xm, double xw, double &d, double &s, double &p0) {
  const double ao = xm - xw, bo = xm + xw;
  double y1 = 0, y2 = 0;
  int n0 = 0;
  for (double v : z) {
    if (v >= ao && v <= bo) {
      y1 += v;
      y2 += v * v;
      ++n0;
    }
  }
  if (n0 < 1) return false;
  y1 /= n0;
  y2 /= n0;
  const double that = static_cast<double>(n0) / static_cast<double>(z.size());
  double b1 = 0, b2 = 0;
  bool ok = false;
  for (int j = 1; j <= 35; ++j) {
    const double be1 = d / (s * s), be2 = -1 / (2 * s * s);
    const double aa = (ao - d) / s, bb = (bo - d) / s;
    const double H0 = pn(bb) - pn(aa), fa = dn(aa), fb = dn(bb);
    const double H1 = fa - fb, H2 = H0 + aa * fa - bb * fb;
    const double H3 = (2 + aa * aa) * fa - (2 + bb * bb) * fb;
    const double H4 = 3 * H0 + (3 * aa + aa * aa * aa) * fa - (3 * bb + bb * bb * bb) * fb;
    const double H[5] = {H0, H1, H2, H3, H4};
    // E[k] = sum_i choose(k, i) s^k r^(k - i) H_i / H0, r = d / s
    const double r = d / s;
    double E[5];
    for (int k = 0; k < 5; ++k) {
      double v = 0, c = 1;
      for (int i = 0; i <= k; ++i) {
        v += c * std::pow(r, k - i) * H[i];
        c = c * (k - i) / (i + 1);
      }
      E[k] = std::pow(s, k) * v / H0;
    }
    const double m1 = E[1], m2 = E[2];
    const double v11 = E[2] - m1 * m1, v12 = E[3] - m1 * m2, v22 = E[4] - m2 * m2;
    const double det = v11 * v22 - v12 * v12;
    const double r1 = y1 - m1, r2 = y2 - m2;
    double s1 = (v22 * r1 - v12 * r2) / det, s2 = (v11 * r2 - v12 * r1) / det;
    const double f = 1 + 1.0 / (static_cast<double>(j) * j);
    b1 = be1 + s1 / f;
    b2 = be2 + s2 / f;
    if (b2 > 0) {
      b1 = be1 + 0.1 * s1 / f;
      b2 = be2 + 0.1 * s2 / f;
    }
    if (std::isnan(b2) || b2 >= 0) {
      ok = false;
      break;
    }
    ok = true;
    d = -b1 / (2 * b2);
    s = 1 / std::sqrt(-2 * b2);
    if (std::sqrt((b1 - be1) * (b1 - be1) + (b2 - be2) * (b2 - be2)) < 1e-5) break;
  }
  if (!ok || !std::isfinite(d) || !(s > 0)) return false;
  const double aa = (ao - d) / s, bb = (bo - d) / s;
  p0 = that / (pn(bb) - pn(aa));
  return std::isfinite(p0);
}

// The local false discovery rate of each of the m z's by the method and
// defaults of locfdr() (Efron 2004; package locfdr 1.1-8: 120 breaks,
// natural spline of 7 df, the empirical null by maximum likelihood);
// nul: its delta, sigma and p0. false when the null cannot be fitted.
bool lfdr(const std::vector<double> &zz, double *fdr, double *nul) {
  const int N = static_cast<int>(zz.size());
  if (N < 3) return false;
  std::vector<double> sz(zz);
  std::sort(sz.begin(), sz.end());
  const double lo = sz.front(), up = sz.back();
  if (!(up > lo)) return false;
  const int bre = 120, K = bre - 1;
  std::vector<double> br(static_cast<std::size_t>(bre)), x(static_cast<std::size_t>(K)), y(static_cast<std::size_t>(K), 0.0);
  const double by = (up - lo) / (bre - 1);
  for (int i = 0; i < bre - 1; ++i) br[static_cast<std::size_t>(i)] = lo + i * by;
  br[static_cast<std::size_t>(bre - 1)] = up;
  // hist(): right-closed bins, the first closed, breaks moved by a fuzz
  std::vector<double> dif(static_cast<std::size_t>(K));
  for (int i = 0; i < K; ++i) dif[static_cast<std::size_t>(i)] = br[static_cast<std::size_t>(i + 1)] - br[static_cast<std::size_t>(i)];
  std::sort(dif.begin(), dif.end());
  const double fz = 1e-7 * (K % 2 ? dif[static_cast<std::size_t>(K / 2)]
                                  : 0.5 * (dif[static_cast<std::size_t>(K / 2 - 1)] + dif[static_cast<std::size_t>(K / 2)]));
  std::vector<double> fb(br);
  fb[0] -= fz;
  for (int i = 1; i < bre; ++i) fb[static_cast<std::size_t>(i)] += fz;
  for (double v : zz) {
    int j = static_cast<int>(std::upper_bound(fb.begin(), fb.end(), v) - fb.begin()) - 1;
    // (fb[j], fb[j + 1]]: a value equal to a break belongs below it
    if (j >= 0 && j < bre && v == fb[static_cast<std::size_t>(j)]) --j;
    if (j < 0 || j >= K) continue;
    y[static_cast<std::size_t>(j)] += 1;
  }
  for (int i = 0; i < K; ++i) x[static_cast<std::size_t>(i)] = (br[static_cast<std::size_t>(i)] + br[static_cast<std::size_t>(i + 1)]) / 2;
  // the natural cubic splines with knots at the quantiles 1/7..6/7 of x
  // and boundary knots at its ends (the space of ns(x, df = 7) and an
  // intercept), in the truncated power basis on x scaled to [0, 1]
  std::vector<double> kn(8);
  kn[0] = x.front();
  kn[7] = x.back();
  for (int k = 1; k <= 6; ++k) kn[static_cast<std::size_t>(k)] = q7(x, k / 7.0);
  const double x0 = x.front(), xs = x.back() - x.front();
  const int nb = 8;
  std::vector<double> B(static_cast<std::size_t>(K) * nb);
  auto tp = [](double v) { return v > 0 ? v * v * v : 0.0; };
  for (int i = 0; i < K; ++i) {
    const double u = (x[static_cast<std::size_t>(i)] - x0) / xs;
    std::vector<double> t(8);
    for (int k = 0; k < 8; ++k) t[static_cast<std::size_t>(k)] = (kn[static_cast<std::size_t>(k)] - x0) / xs;
    auto dk = [&](int k) { return (tp(u - t[static_cast<std::size_t>(k)]) - tp(u - t[7])) / (t[7] - t[static_cast<std::size_t>(k)]); };
    B[static_cast<std::size_t>(i)] = 1.0;
    B[static_cast<std::size_t>(K) + i] = u;
    for (int k = 0; k < 6; ++k) B[static_cast<std::size_t>(k + 2) * K + i] = dk(k) - dk(6);
  }
  // the Poisson regression of the counts on the basis, by Newton steps
  // from log(y + 0.1) (glm()'s start)
  std::vector<double> cf(static_cast<std::size_t>(nb), 0.0), L(static_cast<std::size_t>(nb) * nb), f(static_cast<std::size_t>(K));
  {
    // the start: least squares of log(y + 0.1) on the basis
    std::vector<double> g(static_cast<std::size_t>(nb), 0.0);
    for (double &v : L) v = 0.0;
    for (int i = 0; i < K; ++i) {
      const double t = std::log(y[static_cast<std::size_t>(i)] + 0.1);
      for (int j = 0; j < nb; ++j) {
        g[static_cast<std::size_t>(j)] += B[static_cast<std::size_t>(j) * K + i] * t;
        for (int m = j; m < nb; ++m) {
          L[static_cast<std::size_t>(j) * nb + m] += B[static_cast<std::size_t>(j) * K + i] * B[static_cast<std::size_t>(m) * K + i];
        }
      }
    }
    if (!chol(L.data(), nb)) return false;
    csolve(L.data(), nb, g.data());
    cf = g;
  }
  bool conv = false;
  for (int it = 0; it < 200 && !conv; ++it) {
    std::vector<double> g(static_cast<std::size_t>(nb), 0.0);
    for (double &v : L) v = 0.0;
    for (int i = 0; i < K; ++i) {
      double e = 0;
      for (int j = 0; j < nb; ++j) e += B[static_cast<std::size_t>(j) * K + i] * cf[static_cast<std::size_t>(j)];
      const double mu = std::exp(e), r = y[static_cast<std::size_t>(i)] - mu;
      for (int j = 0; j < nb; ++j) {
        const double bj = B[static_cast<std::size_t>(j) * K + i];
        g[static_cast<std::size_t>(j)] += r * bj;
        for (int m = j; m < nb; ++m) L[static_cast<std::size_t>(j) * nb + m] += mu * bj * B[static_cast<std::size_t>(m) * K + i];
      }
    }
    if (!chol(L.data(), nb)) return false;
    csolve(L.data(), nb, g.data());
    double mx = 0;
    for (int j = 0; j < nb; ++j) {
      cf[static_cast<std::size_t>(j)] += g[static_cast<std::size_t>(j)];
      mx = std::fmax(mx, std::fabs(g[static_cast<std::size_t>(j)]) / (1 + std::fabs(cf[static_cast<std::size_t>(j)])));
    }
    conv = mx < 1e-11;
  }
  if (!conv) return false;
  double sf = 0;
  for (int i = 0; i < K; ++i) {
    double e = 0;
    for (int j = 0; j < nb; ++j) e += B[static_cast<std::size_t>(j) * K + i] * cf[static_cast<std::size_t>(j)];
    f[static_cast<std::size_t>(i)] = std::exp(e);
    sf += f[static_cast<std::size_t>(i)];
  }
  // the empirical null: two passes of locmle(), the first from N(0, 1)
  // over the median -+ b times a robust scale, the second over the first
  // estimate's delta -+ b sigma
  const double med = N % 2 ? sz[static_cast<std::size_t>(N / 2)] : 0.5 * (sz[static_cast<std::size_t>(N / 2 - 1)] + sz[static_cast<std::size_t>(N / 2)]);
  const double sc = (q7(sz, 0.75) - q7(sz, 0.25)) / (2 * 0.67448975019608174);
  const double b = 4.3 * std::exp(-0.26 * std::log10(static_cast<double>(N)));
  double d = 0, s = 1, p0 = 0;
  if (!locmle(zz, med, b * sc, d, s, p0)) return false;
  const double d1 = d, s1 = s;
  if (!locmle(zz, d1, b * s1, d, s, p0)) return false;
  nul[0] = d;
  nul[1] = s;
  nul[2] = p0;
  std::vector<double> f0(static_cast<std::size_t>(K)), fd(static_cast<std::size_t>(K));
  double s0 = 0;
  for (int i = 0; i < K; ++i) {
    f0[static_cast<std::size_t>(i)] = dn((x[static_cast<std::size_t>(i)] - d) / s) / s;
    s0 += f0[static_cast<std::size_t>(i)];
  }
  for (int i = 0; i < K; ++i) {
    fd[static_cast<std::size_t>(i)] = std::fmin(p0 * (sf * f0[static_cast<std::size_t>(i)] / s0) / f[static_cast<std::size_t>(i)], 1.0);
  }
  // 1 between the outermost bins at 1 on either side of delta, and within
  // delta -+ sigma
  double xl = d, xh = d;
  bool any = false;
  for (int i = 0; i < K; ++i) {
    const double xi = x[static_cast<std::size_t>(i)];
    if (xi <= d && fd[static_cast<std::size_t>(i)] == 1 && (!any || xi < xl)) {
      xl = xi;
      any = true;
    }
  }
  any = false;
  for (int i = 0; i < K; ++i) {
    const double xi = x[static_cast<std::size_t>(i)];
    if (xi >= d && fd[static_cast<std::size_t>(i)] == 1 && (!any || xi > xh)) {
      xh = xi;
      any = true;
    }
  }
  for (int i = 0; i < K; ++i) {
    const double xi = x[static_cast<std::size_t>(i)];
    if ((xi >= xl && xi <= xh) || (xi >= d - s && xi <= d + s)) fd[static_cast<std::size_t>(i)] = 1;
  }
  // linear interpolation at each z, constant beyond the end bins
  for (int q = 0; q < N; ++q) {
    const double v = zz[static_cast<std::size_t>(q)];
    if (v <= x.front()) {
      fdr[q] = fd.front();
    } else if (v >= x.back()) {
      fdr[q] = fd.back();
    } else {
      const int j = static_cast<int>(std::upper_bound(x.begin(), x.end(), v) - x.begin()) - 1;
      const double t = (v - x[static_cast<std::size_t>(j)]) / (x[static_cast<std::size_t>(j + 1)] - x[static_cast<std::size_t>(j)]);
      fdr[q] = fd[static_cast<std::size_t>(j)] + t * (fd[static_cast<std::size_t>(j + 1)] - fd[static_cast<std::size_t>(j)]);
    }
  }
  return true;
}

}  // namespace

int clab_ibench(const double *y, const int *off, const double *x, int n, int p, int no, int minn, int dr,
                int local, double *res, int *ok, double *nul) {
  if (n < 2 || p < 0 || no < 2 || minn < 1) return CLAB_EARG;
  if (static_cast<long long>(n) * (p > 0 ? p : 1) > 2000000000LL) return CLAB_EARG;
  // the regressions' (p + 2) x (p + 2) work arrays, bounded as the
  // binding says (at most 8190 columns)
  if (p > 8190) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (off[i] < 0 || off[i] >= no || !(y[i] == 0 || y[i] == 1)) return CLAB_EDOM;
  }
  const std::size_t un = static_cast<std::size_t>(n), up1 = static_cast<std::size_t>(p) + 1;
  std::vector<double> d(un), cf(up1), se(up1), w(2 * up1 * up1 + 4 * up1);
  std::vector<int> iw(un);
  // the doubly robust fit: the design (1, officer, x), weights
  const int kd = p + 2;
  std::vector<double> Z(dr ? un * static_cast<std::size_t>(kd) : 0), wt(dr ? un : 0), bd(static_cast<std::size_t>(kd)),
      Ld(static_cast<std::size_t>(kd) * kd);
  // res: no x 8, column-major: stops, rate, benchmark, effective size of
  // the benchmark, z, p, q (Benjamini-Hochberg), local fdr
  for (std::size_t c = 0; c < 8 * static_cast<std::size_t>(no); ++c) res[c] = kutil::nan();
  for (int k = 0; k < 3; ++k) nul[k] = kutil::nan();
  for (int o = 0; o < no; ++o) {
    if (kutil::stop(static_cast<double>(n) * (p + 2) * (p + 2) * (dr ? 40.0 : 10.0))) return CLAB_EINTR;
    ok[o] = 0;
    double nk = 0, s1 = 0;
    for (int i = 0; i < n; ++i) {
      d[static_cast<std::size_t>(i)] = off[i] == o ? 1.0 : 0.0;
      if (off[i] == o) {
        nk += 1;
        s1 += y[i];
      }
    }
    res[o] = nk;
    if (nk < minn || nk >= n) continue;
    res[no + o] = s1 / nk;
    double dev = 0.0;
    const int rc = clab_lgtreg(x, n, p, d.data(), 1e-10, 50, cf.data(), se.data(), &dev, w.data(), iw.data());
    if (rc == CLAB_EINTR) return rc;
    if (rc != CLAB_OK) continue;
    kutil::ksum sw, sw2, swy;
    for (int i = 0; i < n; ++i) {
      double eta = cf[0];
      for (int j = 0; j < p; ++j) eta += cf[static_cast<std::size_t>(j) + 1] * x[static_cast<std::size_t>(j) * un + i];
      const double wi = std::exp(eta);  // e / (1 - e)
      if (dr) wt[static_cast<std::size_t>(i)] = off[i] == o ? 1.0 : wi;
      if (off[i] == o) continue;
      sw.add(wi);
      sw2.add(wi * wi);
      swy.add(wi * y[i]);
    }
    const double W = sw.val();
    if (!(W > 0) || !std::isfinite(W)) continue;
    const double b = swy.val() / W, ess = W * W / sw2.val(), r = res[no + o];
    res[3 * no + o] = ess;
    if (!dr) {
      const double v = r * (1 - r) / nk + b * (1 - b) / ess;
      res[2 * no + o] = b;
      if (v > 0) {
        res[4 * no + o] = (r - b) / std::sqrt(v);
        res[5 * no + o] = std::erfc(std::fabs(res[4 * no + o]) / std::sqrt(2.0));
      }
      ok[o] = 1;
      continue;
    }
    // doubly robust: the weighted logistic regression of y on the
    // officer's indicator and x over all stops (the officer's weighted 1,
    // the others e / (1 - e)); the benchmark is the mean over the
    // officer's stops of the fitted chance with the indicator at 0, and
    // the difference from the officer's rate has the sandwich variance
    for (int i = 0; i < n; ++i) {
      Z[static_cast<std::size_t>(i)] = 1.0;
      Z[un + static_cast<std::size_t>(i)] = d[static_cast<std::size_t>(i)];
      for (int j = 0; j < p; ++j) Z[static_cast<std::size_t>(j + 2) * un + i] = x[static_cast<std::size_t>(j) * un + i];
    }
    if (!wlogit(Z.data(), y, wt.data(), n, kd, bd.data(), Ld.data())) continue;
    std::vector<double> gc(static_cast<std::size_t>(kd), 0.0), M(static_cast<std::size_t>(kd) * kd, 0.0);
    double m1 = 0, m0 = 0;
    for (int i = 0; i < n; ++i) {
      double e0 = bd[0];
      for (int j = 0; j < p; ++j) e0 += bd[static_cast<std::size_t>(j + 2)] * x[static_cast<std::size_t>(j) * un + i];
      const double e = e0 + bd[1] * d[static_cast<std::size_t>(i)], mu = lgs(e);
      const double rr = wt[static_cast<std::size_t>(i)] * (y[i] - mu);
      for (int j = 0; j < kd; ++j) {
        for (int m = 0; m < kd; ++m) {
          M[static_cast<std::size_t>(j) * kd + m] += rr * rr * Z[static_cast<std::size_t>(j) * un + i] * Z[static_cast<std::size_t>(m) * un + i];
        }
      }
      if (off[i] != o) continue;
      const double u1 = lgs(e0 + bd[1]), u0 = lgs(e0);
      m1 += u1;
      m0 += u0;
      // the gradient of mean(u1 - u0) over the officer's stops
      gc[0] += u1 * (1 - u1) - u0 * (1 - u0);
      gc[1] += u1 * (1 - u1);
      for (int j = 0; j < p; ++j) {
        gc[static_cast<std::size_t>(j + 2)] += (u1 * (1 - u1) - u0 * (1 - u0)) * x[static_cast<std::size_t>(j) * un + i];
      }
    }
    for (double &v : gc) v /= nk;
    // v = g' H^-1 M H^-1 g
    std::vector<double> hg(gc);
    csolve(Ld.data(), kd, hg.data());
    double v = 0;
    for (int j = 0; j < kd; ++j) {
      for (int m = 0; m < kd; ++m) v += hg[static_cast<std::size_t>(j)] * M[static_cast<std::size_t>(j) * kd + m] * hg[static_cast<std::size_t>(m)];
    }
    const double te = (m1 - m0) / nk;
    res[2 * no + o] = m0 / nk;
    if (v > 0) {
      res[4 * no + o] = te / std::sqrt(v);
      res[5 * no + o] = std::erfc(std::fabs(res[4 * no + o]) / std::sqrt(2.0));
    }
    ok[o] = 1;
  }
  // Benjamini-Hochberg over the officers with a p-value
  std::vector<int> ix;
  for (int o = 0; o < no; ++o) {
    if (!std::isnan(res[5 * no + o])) ix.push_back(o);
  }
  const int m = static_cast<int>(ix.size());
  std::sort(ix.begin(), ix.end(), [&](int u, int v) {
    return res[5 * no + u] > res[5 * no + v] || (res[5 * no + u] == res[5 * no + v] && u < v);
  });
  double run = 1.0;
  for (int k = 0; k < m; ++k) {
    const int o = ix[static_cast<std::size_t>(k)];
    const double q = res[5 * no + o] * m / (m - k);
    run = std::fmin(run, q);
    res[6 * no + o] = run;
  }
  if (local) {
    std::vector<double> zs;
    std::vector<int> zo;
    for (int o = 0; o < no; ++o) {
      if (!std::isnan(res[4 * no + o])) {
        zs.push_back(res[4 * no + o]);
        zo.push_back(o);
      }
    }
    std::vector<double> fd(zs.size());
    if (!lfdr(zs, fd.data(), nul)) return CLAB_EDOM;
    for (std::size_t k = 0; k < zs.size(); ++k) res[7 * no + zo[k]] = fd[k];
  }
  return CLAB_OK;
}
