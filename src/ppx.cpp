// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, point patterns and segregation: the cross-type K function
// with Ripley's isotropic correction and the kernel segregation test, the
// quadrat count test, Morisita's index, and the dissimilarity, entropy,
// Atkinson and Gini indices of residential segregation.
#include <algorithm>
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

namespace {

constexpr double pi = 3.14159265358979323846;

// Ripley's isotropic weight in the rectangle w of the circle of radius r
// about (x, y): the circle's length over its length inside
double ripw(double x, double y, double r, const double *w) {
  const double d[4] = {x - w[0], y - w[2], w[1] - x, w[3] - y};  // left, down, right, up
  double a[4], out = 0;
  for (int k = 0; k < 4; ++k) {
    a[k] = d[k] < r ? std::acos(d[k] / r) : 0.0;
    out += 2 * a[k];
  }
  // arcs beyond two adjacent sides overlap past the corner
  for (int k = 0; k < 4; ++k) out -= std::max(0.0, a[k] + a[(k + 1) % 4] - pi / 2);
  return 2 * pi / (2 * pi - out);
}

}  // namespace

int clab_crossk(const double *x, const double *y, const int *mk, int n, int ti, int tj,
                const double *win, const double *r, int nr, double sigma, int nsim,
                unsigned int seed, double *k, double *seg) {
  if (n < 2 || nr < 1 || nsim < 0) return CLAB_EARG;
  if (!(win[0] < win[1] && win[2] < win[3])) return CLAB_EDOM;
  int ni = 0, nj = 0, nt = 0;
  for (int a = 0; a < n; ++a) {
    if (!(x[a] >= win[0] && x[a] <= win[1] && y[a] >= win[2] && y[a] <= win[3])) return CLAB_EDOM;
    if (mk[a] < 0) return CLAB_EARG;
    nt = std::max(nt, mk[a] + 1);
    ni += mk[a] == ti;
    nj += mk[a] == tj;
  }
  for (int q = 0; q < nr; ++q) {
    if (!kutil::nonneg(r[q]) || (q > 0 && !(r[q] > r[q - 1]))) return CLAB_EDOM;
  }
  if (ni < 1 || nj < 1 || (ti == tj && ni < 2)) return CLAB_EDOM;
  // K_ij(r) = A / (n_i n_j) sum over pairs (a of type i, b of type j, a !=
  // b) with d <= r of the isotropic weight at a
  const double area = (win[1] - win[0]) * (win[3] - win[2]);
  std::vector<double> s(nr, 0.0);
  for (int a = 0; a < n; ++a) {
    if (mk[a] != ti) continue;
    for (int b = 0; b < n; ++b) {
      if (b == a || mk[b] != tj) continue;
      const double d = std::sqrt((x[a] - x[b]) * (x[a] - x[b]) + (y[a] - y[b]) * (y[a] - y[b]));
      if (d > r[nr - 1]) continue;
      const double e = d > 0 ? ripw(x[a], y[a], d, win) : 1.0;
      for (int q = nr - 1; q >= 0 && d <= r[q]; --q) s[q] += e;
    }
  }
  const double den = static_cast<double>(ni) * (ti == tj ? ni - 1 : nj);
  for (int q = 0; q < nr; ++q) k[q] = area * s[q] / den;
  if (!(sigma > 0)) return CLAB_OK;
  // the segregation test (Diggle, Zheng and Durr 2005): leave-one-out
  // kernel shares of each type at each point against the overall shares
  std::vector<double> kk(static_cast<std::size_t>(n) * n), den2(n, 0.0), bar(nt, 0.0);
  for (int a = 0; a < n; ++a) {
    bar[mk[a]] += 1.0 / n;
    for (int b = 0; b < n; ++b) {
      const double d2 = (x[a] - x[b]) * (x[a] - x[b]) + (y[a] - y[b]) * (y[a] - y[b]);
      const double v = a == b ? 0.0 : std::exp(-d2 / (2 * sigma * sigma));
      kk[static_cast<std::size_t>(a) * n + b] = v;
      den2[a] += v;
    }
    if (!(den2[a] > 0)) return CLAB_EDOM;
  }
  std::vector<double> pm(nt);
  auto stat = [&](const int *lab) {
    double t = 0;
    for (int a = 0; a < n; ++a) {
      std::fill(pm.begin(), pm.end(), 0.0);
      const double *row = kk.data() + static_cast<std::size_t>(a) * n;
      for (int b = 0; b < n; ++b) pm[lab[b]] += row[b];
      for (int m = 0; m < nt; ++m) {
        const double e = pm[m] / den2[a] - bar[m];
        t += e * e;
      }
    }
    return t;
  };
  const double t0 = stat(mk);
  std::vector<int> lab(mk, mk + n);
  kutil::rng g(seed);
  double ge = 0;
  for (int q = 0; q < nsim; ++q) {
    for (int a = n - 1; a > 0; --a) std::swap(lab[a], lab[g.below(a + 1)]);
    ge += stat(lab.data()) >= t0;
  }
  seg[0] = t0;
  seg[1] = nsim > 0 ? (1 + ge) / (nsim + 1.0) : kutil::nan();
  return CLAB_OK;
}

int clab_qdrtst(const double *x, const double *y, int n, const double *win, int nx, int ny,
                double *cnt, double *st) {
  if (n < 1 || nx < 1 || ny < 1 || nx * ny < 2) return CLAB_EARG;
  if (!(win[0] < win[1] && win[2] < win[3])) return CLAB_EDOM;
  std::fill(cnt, cnt + static_cast<std::size_t>(nx) * ny, 0.0);
  for (int a = 0; a < n; ++a) {
    if (!(x[a] >= win[0] && x[a] <= win[1] && y[a] >= win[2] && y[a] <= win[3])) return CLAB_EDOM;
    const int i =
        std::min(nx - 1, static_cast<int>(std::floor((x[a] - win[0]) / (win[1] - win[0]) * nx)));
    const int j =
        std::min(ny - 1, static_cast<int>(std::floor((y[a] - win[2]) / (win[3] - win[2]) * ny)));
    // rows from the top (the highest y), columns from the left
    cnt[(ny - 1 - j) + static_cast<std::size_t>(i) * ny] += 1;
  }
  const double e = static_cast<double>(n) / (nx * ny);
  kutil::ksum s;
  for (std::size_t q = 0; q < static_cast<std::size_t>(nx) * ny; ++q)
    s.add((cnt[q] - e) * (cnt[q] - e) / e);
  st[0] = s.val();
  st[1] = nx * ny - 1;
  st[2] = n / ((win[1] - win[0]) * (win[3] - win[2]));
  return CLAB_OK;
}

int clab_morsta(const double *x, int n, double chiu, double chil, double *out) {
  if (n < 2) return CLAB_EARG;
  double s = 0, s2 = 0;
  for (int i = 0; i < n; ++i) {
    if (!kutil::nonneg(x[i])) return CLAB_EDOM;
    s += x[i];
    s2 += x[i] * x[i];
  }
  if (!(s > 1)) return CLAB_EDOM;
  // Morisita's index, the uniform and clumped indices at the chi-square
  // quantiles chiu (upper) and chil (lower), and Smith-Gill's standardised
  // index
  const double id = n * (s2 - s) / (s * s - s);
  const double mu = (chiu - n + s) / (s - 1), mc = (chil - n + s) / (s - 1);
  double st = id;
  if (id >= mc && mc > 1) {
    st = 0.5 + 0.5 * (id - mc) / (n - mc);
  } else if (mc > id && id >= 1) {
    st = 0.5 * (id - 1) / (mc - 1);
  } else if (1 > id && id > mu) {
    st = -0.5 * (id - 1) / (mu - 1);
  } else if (1 > mu && mu > id) {
    st = -0.5 + 0.5 * (id - mu) / mu;
  }
  out[0] = id;
  out[1] = mc;
  out[2] = mu;
  out[3] = st;
  out[4] = id * (s - 1) + n - s;
  return CLAB_OK;
}

int clab_dissix(const double *x, int n, int g, double delta, const int *rp, const int *ci,
                const double *w, double *d, double *grp, double *dm) {
  if (n < 1 || g < 2) return CLAB_EARG;
  if (!(delta > 0 && delta < 1)) return CLAB_EDOM;
  std::vector<double> t(n, 0.0), tot(g, 0.0);
  for (int k = 0; k < g; ++k) {
    for (int i = 0; i < n; ++i) {
      const double v = x[i + static_cast<std::size_t>(k) * n];
      if (!kutil::nonneg(v)) return CLAB_EDOM;
      t[i] += v;
      tot[k] += v;
    }
  }
  double T = 0;
  for (int k = 0; k < g; ++k) {
    if (!(tot[k] > 0)) return CLAB_EDOM;
    T += tot[k];
  }
  auto X = [&](int i, int k) { return x[i + static_cast<std::size_t>(k) * n]; };
  // Duncan's D between each pair of groups
  for (int a = 0; a < g; ++a) {
    d[a + static_cast<std::size_t>(a) * g] = 0;
    for (int b = a + 1; b < g; ++b) {
      kutil::ksum s;
      for (int i = 0; i < n; ++i) s.add(std::fabs(X(i, a) / tot[a] - X(i, b) / tot[b]));
      d[a + static_cast<std::size_t>(b) * g] = d[b + static_cast<std::size_t>(a) * g] =
          0.5 * s.val();
    }
  }
  auto ent = [](double p) {
    return (p > 0 ? -p * std::log(p) : 0.0) + (p < 1 ? -(1 - p) * std::log(1 - p) : 0.0);
  };
  for (int k = 0; k < g; ++k) {
    const double P = tot[k] / T;
    kutil::ksum is, h, at, gi;
    for (int i = 0; i < n; ++i) {
      if (!(t[i] > 0)) continue;
      const double p = X(i, k) / t[i];
      is.add(std::fabs(X(i, k) / tot[k] - (t[i] - X(i, k)) / (T - tot[k])));
      h.add(t[i] * (ent(P) - ent(p)));
      at.add(std::pow(1 - p, 1 - delta) * std::pow(p, delta) * t[i] / (P * T));
      for (int j = 0; j < n; ++j) {
        if (t[j] > 0) gi.add(t[i] * t[j] * std::fabs(p - X(j, k) / t[j]));
      }
    }
    double *o = grp + 4 * k;
    o[0] = 0.5 * is.val();
    o[1] = P < 1 ? h.val() / (ent(P) * T) : kutil::nan();
    o[2] = P < 1 ? 1 - P / (1 - P) * std::pow(at.val(), 1 / (1 - delta)) : kutil::nan();
    o[3] = P < 1 ? gi.val() / (2 * T * T * P * (1 - P)) : kutil::nan();
  }
  if (rp == nullptr) return CLAB_OK;
  // Morrill's adjustment: D less the weighted mean absolute difference of
  // the pair's shares between neighbouring units
  for (int a = 0; a < g; ++a) {
    dm[a + static_cast<std::size_t>(a) * g] = 0;
    for (int b = a + 1; b < g; ++b) {
      kutil::ksum s, c;
      auto z = [&](int i) {
        const double u = X(i, a) + X(i, b);
        return u > 0 ? X(i, a) / u : 0.0;
      };
      for (int i = 0; i < n; ++i) {
        if (!(t[i] > 0)) continue;
        for (int q = rp[i]; q < rp[i + 1]; ++q) {
          if (!(t[ci[q]] > 0)) continue;
          s.add(w[q] * std::fabs(z(i) - z(ci[q])));
          c.add(w[q]);
        }
      }
      if (!(c.val() > 0)) return CLAB_EDOM;
      dm[a + static_cast<std::size_t>(b) * g] = dm[b + static_cast<std::size_t>(a) * g] =
          d[a + static_cast<std::size_t>(b) * g] - s.val() / c.val();
    }
  }
  return CLAB_OK;
}
