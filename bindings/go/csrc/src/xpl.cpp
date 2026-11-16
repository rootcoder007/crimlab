// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, the explanation of a prediction function (crimlab.h, xpimp()
// and the functions after it): the rows at which the function is asked to
// predict, and the summaries of its answers there. The function itself is
// the caller's (an R function, or the learners of suplrn()): each method
// has a design kernel, which writes the rows, and a summary kernel, which
// reads the predictions at them. Matrices are column-major. No R here.
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

namespace {

using i64 = std::int64_t;

// Fisher-Yates: pm[0..n) a uniform permutation of 0..n-1 drawn from g.
void shuf(int *pm, int n, kutil::rng &g) {
  for (int i = 0; i < n; ++i) pm[i] = i;
  for (int i = n - 1; i > 0; --i) {
    const int k = g.below(i + 1), t = pm[i];
    pm[i] = pm[k];
    pm[k] = t;
  }
}

// The mean of v[0..n) and its standard deviation (divisor n - 1; NaN for
// one value).
void msd(const double *v, int n, double *m, double *s) {
  kutil::ksum a;
  for (int i = 0; i < n; ++i) a.add(v[i]);
  *m = a.val() / n;
  if (n < 2) {
    *s = kutil::nan();
    return;
  }
  kutil::ksum q;
  for (int i = 0; i < n; ++i) q.add((v[i] - *m) * (v[i] - *m));
  *s = std::sqrt(q.val() / (n - 1));
}

bool fin(const double *v, i64 n) {
  for (i64 i = 0; i < n; ++i) {
    if (!std::isfinite(v[i])) return false;
  }
  return true;
}

}  // namespace

extern "C" {

// Permutation importance (Breiman 2001; Fisher, Rudin and Dominici 2019):
// for column j and repetition r the block b = j nrep + r of n rows is x
// with column j shuffled.
int clab_xpperm(const double *x, int n, int p, int nrep, unsigned int seed, double *des, int *perm) CLAB_TRY {
  if (n < 1 || p < 1 || nrep < 1) return CLAB_EARG;
  if (static_cast<double>(n) * p * nrep * p > kutil::maxcell) return CLAB_EARG;
  const i64 N = static_cast<i64>(n) * p * nrep;
  kutil::rng g(seed);
  for (int j = 0; j < p; ++j) {
    for (int r = 0; r < nrep; ++r) {
      const i64 r0 = (static_cast<i64>(j) * nrep + r) * n;
      int *pm = perm + r0;
      shuf(pm, n, g);
      for (int c = 0; c < p; ++c) {
        const double *xc = x + static_cast<i64>(c) * n;
        double *dc = des + static_cast<i64>(c) * N + r0;
        for (int i = 0; i < n; ++i) dc[i] = c == j ? xc[pm[i]] : xc[i];
      }
      if (kutil::stop(static_cast<double>(n) * p)) return CLAB_EINTR;
    }
  }
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)

int clab_xpimp(const double *y, const double *f0, const double *fp, int n, int p, int nrep, int loss, int ratio,
               double *rep, double *st) CLAB_TRY {
  if (n < 1 || p < 1 || nrep < 1 || loss < 0 || loss > 1 || ratio < 0 || ratio > 1) return CLAB_EARG;
  if (static_cast<double>(n) * p * nrep > kutil::maxcell) return CLAB_EARG;
  const i64 N = static_cast<i64>(n) * p * nrep;
  if (!fin(y, n) || !fin(f0, n) || !fin(fp, N)) return CLAB_EDOM;
  auto err = [&](const double *f) {
    kutil::ksum s;
    for (int i = 0; i < n; ++i) {
      const double e = y[i] - f[i];
      s.add(loss == 0 ? e * e : std::fabs(e));
    }
    return s.val() / n;
  };
  const double b = err(f0);
  if (ratio && !(b > 0)) return CLAB_EDOM;
  st[0] = b;
  for (int j = 0; j < p; ++j) {
    for (int r = 0; r < nrep; ++r) {
      const double e = err(fp + (static_cast<i64>(j) * nrep + r) * n);
      rep[j + static_cast<i64>(p) * r] = ratio ? e / b : e - b;
    }
  }
  std::vector<double> v(static_cast<std::size_t>(nrep));
  for (int j = 0; j < p; ++j) {
    for (int r = 0; r < nrep; ++r) v[static_cast<std::size_t>(r)] = rep[j + static_cast<i64>(p) * r];
    msd(v.data(), nrep, st + 1 + j, st + 1 + p + j);
  }
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)

// The grid of a column v[0..n): type 0, ng points from its minimum to its
// maximum as seq(length.out = ng) spaces them; type 1, its quantiles at
// seq(0, 1, length.out = ng) as quantile(type = 7) interpolates them.
// Repeated values are dropped: g holds *m increasing values.
int clab_xpgrid(const double *v, int n, int ng, int type, double *g, int *m) CLAB_TRY {
  if (n < 1 || ng < 2 || type < 0 || type > 1) return CLAB_EARG;
  if (!fin(v, n)) return CLAB_EDOM;
  std::vector<double> s(v, v + n);
  std::sort(s.begin(), s.end());
  const double lo = s.front(), hi = s.back();
  for (int k = 0; k < ng; ++k) {
    if (type == 0) {
      g[k] = k == 0 ? lo : (k == ng - 1 ? hi : lo + kutil::rd(k * ((hi - lo) / (ng - 1))));
      continue;
    }
    // index 1 + (n - 1) p from 1, as quantile() forms it
    const double pk = k == ng - 1 ? 1.0 : k * (1.0 / (ng - 1));
    const double ix = 1 + kutil::rd((n - 1) * pk), a = std::floor(ix);
    const std::size_t la = static_cast<std::size_t>(a) - 1, hb = static_cast<std::size_t>(std::ceil(ix)) - 1;
    const double h = ix - a;
    g[k] = ix > a && s[hb] != s[la] ? kutil::rd((1 - h) * s[la]) + kutil::rd(h * s[hb]) : s[la];
  }
  int u = 1;
  for (int k = 1; k < ng; ++k) {
    if (g[k] != g[u - 1]) g[u++] = g[k];
  }
  *m = u;
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)

// Partial dependence (Friedman 2001) and ceteris paribus profiles: for
// each grid value g[k] and each row i of s (ns rows, p columns), row
// k ns + i of des is s[i, ] with column j set to g[k].
int clab_xpcpd(const double *s, int ns, int p, int j, const double *g, int ng, double *des) CLAB_TRY {
  if (ns < 1 || p < 1 || j < 0 || j >= p || ng < 1) return CLAB_EARG;
  if (static_cast<double>(ns) * ng * p > kutil::maxcell) return CLAB_EARG;
  const i64 N = static_cast<i64>(ns) * ng;
  for (int c = 0; c < p; ++c) {
    const double *sc = s + static_cast<i64>(c) * ns;
    double *dc = des + static_cast<i64>(c) * N;
    for (int k = 0; k < ng; ++k) {
      for (int i = 0; i < ns; ++i) dc[static_cast<i64>(k) * ns + i] = c == j ? g[k] : sc[i];
    }
  }
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)

// pd[k]: the mean of the predictions f over the block k of ns rows.
int clab_xppd(const double *f, int ns, int ng, double *pd) CLAB_TRY {
  if (ns < 1 || ng < 1) return CLAB_EARG;
  if (static_cast<double>(ns) * ng > kutil::maxcell) return CLAB_EARG;
  if (!fin(f, static_cast<i64>(ns) * ng)) return CLAB_EDOM;
  for (int k = 0; k < ng; ++k) {
    kutil::ksum a;
    for (int i = 0; i < ns; ++i) a.add(f[static_cast<i64>(k) * ns + i]);
    pd[k] = a.val() / ns;
  }
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)

// Accumulated local effects (Apley and Zhu 2020, sec. 3, as their ALEPlot
// computes them): the edges of column j of x (n rows, p columns) are z_0,
// its minimum, and its quantiles of type 1 at seq(1/K, 1, length.out = K),
// repeated values dropped (*nz of them, at most K + 1); bin[i], from 0, is
// the interval (z_b, z_b+1] of row i, the first closed below. Rows i and
// n + i of des (2 n rows) are row i of x with column j set to the upper and
// to the lower edge of its interval. CLAB_EDOM when the column has one
// value.
int clab_xpaled(const double *x, int n, int p, int j, int K, double *zz, int *nz, int *bin,
                double *des) CLAB_TRY {
  if (n < 1 || p < 1 || j < 0 || j >= p || K < 1) return CLAB_EARG;
  if (2.0 * n * p > kutil::maxcell) return CLAB_EARG;
  const double *v = x + static_cast<i64>(j) * n;
  if (!fin(v, n)) return CLAB_EDOM;
  std::vector<double> s(v, v + n);
  std::sort(s.begin(), s.end());
  zz[0] = s.front();
  // quantile(type = 1): the order statistic j, or j + 1 when n p passes j,
  // j = floor(n p + 4 eps)
  const double a = 1.0 / K, fz = 4 * DBL_EPSILON;
  for (int k = 1; k <= K; ++k) {
    const double pk = k == K ? 1.0 : (k == 1 ? a : a + kutil::rd((k - 1) * ((1.0 - a) / (K - 1))));
    const double np = n * pk, jj = std::floor(np + fz);
    const int o = std::min(std::max(static_cast<int>(jj) + (np > jj ? 1 : 0), 1), n);
    zz[k] = s[static_cast<std::size_t>(o - 1)];
  }
  int u = 1;
  for (int k = 1; k <= K; ++k) {
    if (zz[k] != zz[u - 1]) zz[u++] = zz[k];
  }
  *nz = u;
  if (u < 2) return CLAB_EDOM;
  for (int i = 0; i < n; ++i) bin[i] = static_cast<int>(std::lower_bound(zz + 1, zz + u, v[i]) - (zz + 1));
  const i64 N = 2 * static_cast<i64>(n);
  for (int c = 0; c < p; ++c) {
    const double *xc = x + static_cast<i64>(c) * n;
    double *dc = des + static_cast<i64>(c) * N;
    for (int i = 0; i < n; ++i) {
      dc[i] = c == j ? zz[bin[i] + 1] : xc[i];
      dc[n + i] = c == j ? zz[bin[i]] : xc[i];
    }
  }
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)

// eff[b]: the mean over the rows of interval b of f[i] - f[n + i]; cnt[b]
// their number; ale (nb + 1 values): 0 and the running sums of eff at the
// edges, less the mean over the rows of the midpoints of their intervals
// (the centring of ALEPlot).
int clab_xpale(const double *f, int n, const int *bin, int nb, double *eff, double *cnt, double *ale) CLAB_TRY {
  if (n < 1 || nb < 1) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (bin[i] < 0 || bin[i] >= nb) return CLAB_EARG;
  }
  if (!fin(f, 2 * static_cast<i64>(n))) return CLAB_EDOM;
  std::vector<kutil::ksum> s(static_cast<std::size_t>(nb));
  for (int b = 0; b < nb; ++b) cnt[b] = 0.0;
  for (int i = 0; i < n; ++i) {
    s[static_cast<std::size_t>(bin[i])].add(f[i] - f[n + i]);
    cnt[bin[i]] += 1.0;
  }
  ale[0] = 0.0;
  for (int b = 0; b < nb; ++b) {
    eff[b] = cnt[b] > 0 ? s[static_cast<std::size_t>(b)].val() / cnt[b] : 0.0;
    ale[b + 1] = ale[b] + eff[b];
  }
  kutil::ksum c;
  for (int b = 0; b < nb; ++b) c.add((ale[b] + ale[b + 1]) / 2 * cnt[b]);
  const double m = c.val() / n;
  for (int b = 0; b <= nb; ++b) ale[b] -= m;
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)

// Sampling Shapley values (Strumbelj and Kononenko 2014, Algorithm 1): for
// instance a of s (ns rows, p columns), feature j and sample t < m, a
// permutation of the features and a row w of x (n rows) are drawn; row q =
// (a p + j) m + t of des takes the instance's values on the features up to
// and including j in the permutation and w's on the others, and row Q + q
// the same with j taken from w (Q = ns p m).
int clab_xpshd(const double *x, int n, const double *s, int ns, int p, int m, unsigned int seed,
               double *des) CLAB_TRY {
  if (n < 1 || ns < 1 || p < 1 || m < 1) return CLAB_EARG;
  if (2.0 * ns * p * m * p > kutil::maxcell) return CLAB_EARG;
  const i64 Q = static_cast<i64>(ns) * p * m, N = 2 * Q;
  std::vector<int> pm(static_cast<std::size_t>(p)), pos(static_cast<std::size_t>(p));
  kutil::rng g(seed);
  for (int a = 0; a < ns; ++a) {
    for (int j = 0; j < p; ++j) {
      for (int t = 0; t < m; ++t) {
        shuf(pm.data(), p, g);
        const int w = g.below(n);
        for (int k = 0; k < p; ++k) pos[static_cast<std::size_t>(pm[static_cast<std::size_t>(k)])] = k;
        const i64 q = (static_cast<i64>(a) * p + j) * m + t;
        const int pj = pos[static_cast<std::size_t>(j)];
        for (int c = 0; c < p; ++c) {
          const double xi = s[static_cast<i64>(c) * ns + a], wi = x[static_cast<i64>(c) * n + w];
          const bool before = pos[static_cast<std::size_t>(c)] < pj;
          des[static_cast<i64>(c) * N + q] = before || c == j ? xi : wi;
          des[static_cast<i64>(c) * N + Q + q] = before ? xi : wi;
        }
      }
      if (kutil::stop(static_cast<double>(m) * p)) return CLAB_EINTR;
    }
  }
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)

// phi[a + ns j]: the mean over the m samples of f[q] - f[Q + q]; se: their
// standard deviation over sqrt(m) (NaN for one sample).
int clab_xpsh(const double *f, int ns, int p, int m, double *phi, double *se) CLAB_TRY {
  if (ns < 1 || p < 1 || m < 1) return CLAB_EARG;
  if (2.0 * ns * p * m > kutil::maxcell) return CLAB_EARG;
  const i64 Q = static_cast<i64>(ns) * p * m;
  if (!fin(f, 2 * Q)) return CLAB_EDOM;
  std::vector<double> v(static_cast<std::size_t>(m));
  for (int a = 0; a < ns; ++a) {
    for (int j = 0; j < p; ++j) {
      const i64 q0 = (static_cast<i64>(a) * p + j) * m, o = a + static_cast<i64>(ns) * j;
      for (int t = 0; t < m; ++t) v[static_cast<std::size_t>(t)] = f[q0 + t] - f[Q + q0 + t];
      msd(v.data(), m, phi + o, se + o);
      se[o] /= std::sqrt(static_cast<double>(m));
    }
  }
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)

}  // extern "C"
