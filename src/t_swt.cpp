// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of the spatial weights (src/swt.cpp), the autocorrelation
// statistics (src/lsa.cpp) and the directional mean: brute force, dense
// re-derivations of each formula written independently of the kernels,
// identities and seeded simulations of the permutation tests.
#include <cmath>
#include <vector>

#include "cgeo.h"
#include "crimlab.h"
#include "tchk.h"

namespace {

using vd = std::vector<double>;
using vi = std::vector<int>;

double dm(const vd &x, const vd &y, int i, int j, int m) {
  const double a = x[i] - x[j], b = y[i] - y[j];
  return m ? std::fabs(a) + std::fabs(b) : std::hypot(a, b);
}

// points: uniform, then a clustered set with exact duplicates
void pts(rng &g, int n, vd &x, vd &y, bool clu) {
  x.resize(n);
  y.resize(n);
  for (int i = 0; i < n; ++i) {
    if (clu && i % 3 == 0 && i > 0) {
      x[i] = x[i - 1];
      y[i] = y[i - 1];
    } else if (clu) {
      x[i] = (i % 4) * 10 + runi(g);
      y[i] = (i % 5) * 10 + runi(g);
    } else {
      x[i] = runi(g, 0, 100);
      y[i] = runi(g, 0, 50);
    }
  }
}

// a dense weights matrix from rows
vd dense(int n, const vi &rp, const vi &ci, const vd &w) {
  vd W(static_cast<std::size_t>(n) * n, 0.0);
  for (int i = 0; i < n; ++i)
    for (int q = rp[i]; q < rp[i + 1]; ++q) W[i * n + ci[q]] = w[q];
  return W;
}

void famtest(tk &T, int big) {
  rng g = {41};
  // k nearest neighbours and bands against brute force
  for (int rep = 0; rep < (big ? 40 : 8); ++rep) {
    const int n = 30 + 17 * rep, k = 1 + rep % 6, met = rep % 2;
    vd x, y;
    pts(g, n, x, y, rep % 4 >= 2);
    vi nb(static_cast<std::size_t>(n) * k), iw(2 * n + 2);
    vd d(static_cast<std::size_t>(n) * k);
    tk_true(T, "swknb.ok", clab_swknb(x.data(), y.data(), n, k, met, nb.data(),
                                      d.data(), iw.data()) == CLAB_OK);
    bool same = true;
    for (int i = 0; i < n; ++i) {
      // brute force: sort all others by (distance, index), take k
      vi o;
      for (int j = 0; j < n; ++j)
        if (j != i) o.push_back(j);
      for (std::size_t a = 1; a < o.size(); ++a)
        for (std::size_t b = a; b > 0; --b) {
          const double e1 = dm(x, y, i, o[b - 1], met), e2 = dm(x, y, i, o[b], met);
          if (e2 < e1 || (e2 == e1 && o[b] < o[b - 1])) std::swap(o[b - 1], o[b]);
          else break;
        }
      for (int q = 0; q < k; ++q) {
        if (nb[i * k + q] != o[q] || d[i * k + q] != dm(x, y, i, o[q], met)) same = false;
      }
    }
    tk_true(T, "swknb.brute", same);
    const double r = 5.0 + rep;
    vi rp(n + 1);
    int nnz = -1;
    clab_swband(x.data(), y.data(), n, r, met, 0, rp.data(), nullptr, nullptr, &nnz,
                iw.data());
    vi ci(nnz + 1);
    vd dd(nnz + 1);
    int nn2 = -1;
    tk_true(T, "swband.ok", clab_swband(x.data(), y.data(), n, r, met, nnz, rp.data(),
                                        ci.data(), dd.data(), &nn2, iw.data()) == CLAB_OK);
    bool ok = nn2 == nnz;
    int t = 0;
    for (int i = 0; i < n && ok; ++i) {
      ok = ok && rp[i] == t;
      for (int j = 0; j < n; ++j) {
        if (j == i || dm(x, y, i, j, met) > r) continue;
        ok = ok && ci[t] == j && dd[t] == dm(x, y, i, j, met);
        ++t;
      }
    }
    tk_true(T, "swband.brute", ok && t == nnz);
  }
  // weights: binary, inverse distance, row-standardised; an empty row
  {
    const int rp[4] = {0, 2, 3, 3};
    const double d[3] = {2.0, 4.0, 0.5};
    double w[3];
    clab_swwt(rp, d, 3, 1.0, 1, w);
    tk_eq(T, "swwt.idw", w[0], (1 / 2.0) / (1 / 2.0 + 1 / 4.0), 1e-15);
    tk_eq(T, "swwt.idw2", w[2], 1.0, 1e-15);
    clab_swwt(rp, nullptr, 3, 0.0, 0, w);
    tk_true(T, "swwt.bin", w[0] == 1 && w[1] == 1 && w[2] == 1);
    clab_swwt(rp, d, 3, 2.0, 0, w);
    tk_eq(T, "swwt.pow", w[1], 1.0 / 16, 1e-15);
    const double z[3] = {0.0, 1, 1};
    tk_ieq(T, "swwt.zero", clab_swwt(rp, z, 3, 1.0, 0, w), CLAB_EDOM);
    const int bad[4] = {0, 2, 1, 3};
    tk_ieq(T, "swwt.rows", clab_swwt(bad, d, 3, 0.0, 0, w), CLAB_EARG);
    tk_ieq(T, "swwt.pw", clab_swwt(rp, d, 3, -1.0, 0, w), CLAB_EDOM);
  }
  // Geary's C, global G, local Moran, local G: dense re-derivations on
  // random weights (k nearest neighbours, row-standardised or inverse
  // distance)
  for (int rep = 0; rep < (big ? 30 : 6); ++rep) {
    const int n = 12 + 5 * rep, k = 2 + rep % 4;
    vd x, y, v(n);
    pts(g, n, x, y, false);
    for (int i = 0; i < n; ++i) v[i] = rexp(g, 0.2) + (i % 3);
    vi nb(static_cast<std::size_t>(n) * k), iw(2 * n + 2);
    vd d(static_cast<std::size_t>(n) * k);
    clab_swknb(x.data(), y.data(), n, k, 0, nb.data(), d.data(), iw.data());
    const std::size_t nk = static_cast<std::size_t>(n) * k;
    vi rp(n + 1), ci(nk);
    vd dd(nk), w(nk);
    for (int i = 0; i <= n; ++i) rp[i] = i * k;
    // rows sorted by neighbour
    for (int i = 0; i < n; ++i) {
      const std::ptrdiff_t r0 = static_cast<std::ptrdiff_t>(i) * k;
      vi o(nb.begin() + r0, nb.begin() + r0 + k);
      vd e(d.begin() + r0, d.begin() + r0 + k);
      for (int a = 0; a < k; ++a)
        for (int b = a + 1; b < k; ++b)
          if (o[b] < o[a]) {
            std::swap(o[a], o[b]);
            std::swap(e[a], e[b]);
          }
      for (int a = 0; a < k; ++a) {
        ci[i * k + a] = o[a];
        dd[i * k + a] = e[a];
      }
    }
    clab_swwt(rp.data(), dd.data(), n, rep % 2 ? 1.0 : 0.0, rep % 3 == 0, w.data());
    const vd W = dense(n, rp, ci, w);
    double S0 = 0, S1 = 0, S2 = 0, m = 0, m2 = 0, m4 = 0;
    for (int i = 0; i < n; ++i) m += v[i] / n;
    for (int i = 0; i < n; ++i) {
      double r = 0, c = 0;
      for (int j = 0; j < n; ++j) {
        S0 += W[i * n + j];
        S1 += 0.5 * (W[i * n + j] + W[j * n + i]) * (W[i * n + j] + W[j * n + i]);
        r += W[i * n + j];
        c += W[j * n + i];
      }
      S2 += (r + c) * (r + c);
      m2 += (v[i] - m) * (v[i] - m);
      m4 += std::pow(v[i] - m, 4);
    }
    double num = 0, cross = 0, all = 0;
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {
        num += W[i * n + j] * (v[i] - v[j]) * (v[i] - v[j]);
        cross += W[i * n + j] * v[i] * v[j];
        if (i != j) all += v[i] * v[j];
      }
    const double N = n, b2 = N * m4 / (m2 * m2);
    double o[9];
    vd wk(n);
    clab_gearc(v.data(), n, rp.data(), ci.data(), w.data(), 0, 1, wk.data(), o);
    tk_eq(T, "gearc.c", o[0], (N - 1) * num / (2 * S0 * m2), 1e-12);
    tk_eq(T, "gearc.vn", o[1], ((2 * S1 + S2) * (N - 1) - 4 * S0 * S0) /
                                   (2 * (N + 1) * S0 * S0), 1e-12);
    const double vr = ((N - 1) * S1 * (N * N - 3 * N + 3 - (N - 1) * b2) -
                       0.25 * (N - 1) * S2 * (N * N + 3 * N - 6 - (N * N - N + 2) * b2) +
                       S0 * S0 * (N * N - 3 - (N - 1) * (N - 1) * b2)) /
                      (N * (N - 2) * (N - 3) * S0 * S0);
    tk_eq(T, "gearc.vr", o[2], vr, 1e-12);
    tk_nan(T, "gearc.noperm", o[7]);
    clab_gorg(v.data(), n, rp.data(), ci.data(), w.data(), 0, 1, wk.data(), o);
    tk_eq(T, "gorg.g", o[0], cross / all, 1e-12);
    tk_eq(T, "gorg.e", o[1], S0 / (N * (N - 1)), 1e-12);
    // local Moran: sum_i I_i = S0 times global Moran's I
    vd lo(7 * n);
    clab_lmorn(v.data(), n, rp.data(), ci.data(), w.data(), 0, 1, iw.data(), wk.data(),
               lo.data());
    double si = 0, moran = 0;
    for (int i = 0; i < n; ++i) si += lo[i];
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) moran += W[i * n + j] * (v[i] - m) * (v[j] - m);
    moran *= N / (S0 * m2);
    tk_eq(T, "lmorn.sum", si, S0 * moran, 1e-10);
    // local G* and G by the formulas written out
    vd gs(n), gx(n);
    clab_locgi(v.data(), n, rp.data(), ci.data(), w.data(), 1, gs.data());
    clab_locgi(v.data(), n, rp.data(), ci.data(), w.data(), 0, gx.data());
    double t1 = 0, t2 = 0;
    for (int i = 0; i < n; ++i) {
      t1 += v[i];
      t2 += v[i] * v[i];
    }
    bool ok = true;
    for (int i = 0; i < n; ++i) {
      double lag = v[i], wi = 1, s1 = 1, lag0 = 0, w0 = 0, s10 = 0;
      for (int j = 0; j < n; ++j) {
        lag += W[i * n + j] * v[j];
        wi += W[i * n + j];
        s1 += W[i * n + j] * W[i * n + j];
        lag0 += W[i * n + j] * v[j];
        w0 += W[i * n + j];
        s10 += W[i * n + j] * W[i * n + j];
      }
      const double xb = t1 / N, s = std::sqrt(t2 / N - xb * xb);
      const double zs = (lag - wi * xb) / (s * std::sqrt((N * s1 - wi * wi) / (N - 1)));
      const double xb0 = (t1 - v[i]) / (N - 1),
                   sd0 = std::sqrt((t2 - v[i] * v[i]) / (N - 1) - xb0 * xb0);
      const double z0 = (lag0 - w0 * xb0) / (sd0 * std::sqrt(((N - 1) * s10 - w0 * w0) / (N - 2)));
      ok = ok && std::fabs(gs[i] - zs) < 1e-10 * (1 + std::fabs(zs)) &&
           std::fabs(gx[i] - z0) < 1e-10 * (1 + std::fabs(z0));
    }
    tk_true(T, "locgi.formula", ok);
  }
  // permutation tests: their size under exchangeable values (each p <=
  // 0.05 with chance 0.05 when nsim = 19), and seeds reproduce
  {
    const int n = 25, k = 4, sets = big ? 400 : 80;
    int rg = 0, rgo = 0, rl = 0;
    vd x, y;
    pts(g, n, x, y, false);
    vi nb(n * k), iw(2 * n + 2), rp(n + 1), ci(n * k);
    vd d(n * k), w(n * k), v(n), wk(n), lo(7 * n);
    clab_swknb(x.data(), y.data(), n, k, 0, nb.data(), d.data(), iw.data());
    for (int i = 0; i <= n; ++i) rp[i] = i * k;
    for (int i = 0; i < n; ++i) {
      vi o(nb.begin() + i * k, nb.begin() + (i + 1) * k);
      for (int a = 0; a < k; ++a)
        for (int b = a + 1; b < k; ++b)
          if (o[b] < o[a]) std::swap(o[a], o[b]);
      for (int a = 0; a < k; ++a) ci[i * k + a] = o[a];
    }
    clab_swwt(rp.data(), nullptr, n, 0.0, 1, w.data());
    for (int s = 0; s < sets; ++s) {
      for (int i = 0; i < n; ++i) v[i] = rexp(g, 1.0);
      double o[9];
      clab_gearc(v.data(), n, rp.data(), ci.data(), w.data(), 19, 7 + s, wk.data(), o);
      rg += o[7] <= 0.05 ? 1 : 0;
      clab_gorg(v.data(), n, rp.data(), ci.data(), w.data(), 19, 7 + s, wk.data(), o);
      rgo += o[5] <= 0.05 ? 1 : 0;
      clab_lmorn(v.data(), n, rp.data(), ci.data(), w.data(), 19, 7 + s, iw.data(),
                 wk.data(), lo.data());
      rl += lo[5 * n] <= 0.05 ? 1 : 0;
    }
    // 4 binomial standard deviations about 0.05 (the local test is two
    // one-sided tests folded: about 0.1)
    const double sd = 4 * std::sqrt(0.05 * 0.95 / sets);
    tk_lt(T, "gearc.size", std::fabs(rg / static_cast<double>(sets) - 0.05), sd);
    tk_lt(T, "gorg.size", std::fabs(rgo / static_cast<double>(sets) - 0.05), sd);
    tk_lt(T, "lmorn.size", std::fabs(rl / static_cast<double>(sets) - 0.1),
          4 * std::sqrt(0.1 * 0.9 / sets));
    double a[9], b[9];
    clab_gearc(v.data(), n, rp.data(), ci.data(), w.data(), 49, 3, wk.data(), a);
    clab_gearc(v.data(), n, rp.data(), ci.data(), w.data(), 49, 3, wk.data(), b);
    tk_true(T, "gearc.seed", a[7] == b[7]);
    // a smooth surface: strong positive autocorrelation, C well below 1
    for (int i = 0; i < n; ++i) v[i] = x[i] + y[i];
    clab_gearc(v.data(), n, rp.data(), ci.data(), w.data(), 99, 3, wk.data(), a);
    tk_true(T, "gearc.power", a[0] < 0.5 && a[7] == 0.01);
  }
  // directional mean: closed forms
  {
    const double pi = cgeo::pi;
    double a[4] = {0, pi / 2, 0, 0}, o[5];
    clab_dirmn(a, nullptr, 2, o);
    tk_eq(T, "dirmn.mean", o[0], pi / 4, 1e-15);
    tk_eq(T, "dirmn.r", o[1], std::cos(pi / 4), 1e-15);
    tk_eq(T, "dirmn.sd", o[3], std::sqrt(-2 * std::log(std::cos(pi / 4))), 1e-15);
    const double q[4] = {0, pi / 2, pi, 3 * pi / 2};
    clab_dirmn(q, nullptr, 4, o);
    tk_true(T, "dirmn.zero", std::isnan(o[0]) && o[1] < 1e-15 && o[2] > 1 - 1e-15);
    const double wv[2] = {3, 1}, r2[2] = {-pi / 2, 0};
    clab_dirmn(r2, wv, 2, o);
    tk_eq(T, "dirmn.w", o[0], 2 * pi - std::atan2(3, 1), 1e-14);
    tk_eq(T, "dirmn.W", o[4], 4, 0);
    const double one = 1;
    clab_dirmn(&one, nullptr, 1, o);
    tk_true(T, "dirmn.one", o[1] == 1 && o[2] == 0 && o[3] == 0);
  }
  // guards
  {
    double x[3] = {0, 1, 2}, y[3] = {0, 0, 0}, d[3], o[9];
    int nb[3], iw[8], rp[4] = {0, 1, 2, 3}, ci[3] = {1, 2, 0}, nnz;
    double w[3] = {1, 1, 1};
    tk_ieq(T, "swknb.k", clab_swknb(x, y, 3, 3, 0, nb, d, iw), CLAB_EARG);
    double nanv[3] = {0, NAN, 1};
    tk_ieq(T, "swknb.nan", clab_swknb(nanv, y, 3, 1, 0, nb, d, iw), CLAB_EDOM);
    tk_ieq(T, "swband.d", clab_swband(x, y, 3, 0, 0, 0, rp, ci, d, &nnz, iw), CLAB_EDOM);
    tk_ieq(T, "swband.m", clab_swband(x, y, 3, 1, 2, 0, rp, ci, d, &nnz, iw), CLAB_EARG);
    tk_ieq(T, "gearc.n", clab_gearc(x, 3, rp, ci, w, 0, 1, d, o), CLAB_EARG);
    const int self[4] = {0, 1, 2, 3}, cs[3] = {0, 2, 0};
    tk_ieq(T, "lsa.self", clab_locgi(x, 3, self, cs, w, 0, o), CLAB_EDOM);
    const int oob[3] = {1, 5, 0};
    tk_ieq(T, "lsa.range", clab_locgi(x, 3, self, oob, w, 0, o), CLAB_EARG);
    double neg[3] = {1, -1, 1};
    tk_ieq(T, "lsa.negw", clab_locgi(x, 3, self, ci, neg, 0, o), CLAB_EDOM);
    tk_ieq(T, "gorg.neg", clab_gorg(neg, 3, rp, ci, w, 0, 1, d, o), CLAB_EARG);
    double cst[3] = {2, 2, 2};
    tk_ieq(T, "lmorn.const", clab_lmorn(cst, 3, rp, ci, w, 0, 1, iw, d, o), CLAB_EDOM);
    tk_ieq(T, "locgi.star", clab_locgi(x, 3, rp, ci, w, 2, o), CLAB_EARG);
    double z[2] = {0, 0};
    tk_ieq(T, "dirmn.w0", clab_dirmn(x, z, 2, o), CLAB_EDOM);
    tk_ieq(T, "dirmn.n", clab_dirmn(x, nullptr, 0, o), CLAB_EARG);
  }
}

}  // namespace

void t_swt(tk &T, int big) { famtest(T, big); }
