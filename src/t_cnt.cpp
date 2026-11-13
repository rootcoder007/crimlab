// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests for clab_cntreg and clab_geeex (src/cnt.cpp): closed
// forms of the group-indicator model, the maximum of the likelihood, the
// sandwich by hand and the estimating equations solved by elimination.
// Part of the self-test (selftest.cpp).
#include <cmath>

#include "crimlab.h"
#include "tchk.h"

namespace {

constexpr int nmax = 40, pmax = 4;

struct out {
  double beta[pmax], vb[pmax * pmax], vc[pmax * pmax], st[8];
  int alias[pmax];
  double work[nmax * pmax + 6 * nmax + pmax + pmax * pmax + 2 * nmax * pmax];
  int iwork[pmax + nmax + 1 + nmax];
};

int fit(const double *y, const double *x, int n, int p, const double *off, const int *cl,
        int ncl, int nb, out &o) {
  return clab_cntreg(y, x, n, p, off, cl, ncl, nb, 100, 1e-13, o.beta, o.vb, o.vc, o.st,
                     o.alias, o.work, o.iwork);
}

double llnb(const double *y, const double *mu, int n, double th) {
  double s = 0;
  for (int i = 0; i < n; ++i) {
    s += std::lgamma(th + y[i]) - std::lgamma(th) - std::lgamma(y[i] + 1) +
         th * std::log(th / (th + mu[i])) + y[i] * std::log(mu[i] / (th + mu[i]));
  }
  return s;
}

void cnt(tk &t) {
  static out o;
  // 1. two groups with exposure: the Poisson rate ratio is the ratio of
  // the pooled rates, its variance 1 / y0 + 1 / y1
  const double y[6] = {3, 5, 10, 2, 0, 0}, e[4] = {2, 3, 4, 1};
  double x[24], off[4];
  for (int i = 0; i < 4; ++i) {
    x[i] = 1;
    x[4 + i] = i >= 2 ? 1 : 0;
    off[i] = std::log(e[i]);
  }
  tk_ieq(t, "cnt.rc", fit(y, x, 4, 2, off, nullptr, 0, 0, o), CLAB_OK);
  tk_eq(t, "cnt.b0", o.beta[0], std::log(8.0 / 5), 1e-12);
  tk_eq(t, "cnt.b1", o.beta[1], std::log(1.5), 1e-12);
  tk_eq(t, "cnt.v1", o.vb[3], 1.0 / 8 + 1.0 / 12, 1e-10);
  tk_eq(t, "cnt.v0", o.vb[0], 1.0 / 8, 1e-10);
  tk_ieq(t, "cnt.rank", static_cast<int>(o.st[6]), 2);
  tk_ieq(t, "cnt.conv", static_cast<int>(o.st[5]), 1);
  // the log-likelihood at the pooled rates
  {
    double ll = 0;
    for (int i = 0; i < 4; ++i) {
      const double mu = e[i] * (i >= 2 ? 12.0 / 5 : 8.0 / 5);
      ll += y[i] * std::log(mu) - mu - std::lgamma(y[i] + 1);
    }
    tk_eq(t, "cnt.ll", o.st[0], ll, 1e-12);
    tk_eq(t, "cnt.aic", o.st[1], -2 * ll + 4, 1e-12);
  }
  // 2. an aliased copy of the group column: dropped, nothing else moves
  for (int i = 0; i < 4; ++i) x[8 + i] = x[4 + i];
  tk_ieq(t, "cnt.alias.rc", fit(y, x, 4, 3, off, nullptr, 0, 0, o), CLAB_OK);
  tk_true(t, "cnt.alias", o.alias[2] == 1 && o.alias[1] == 0 && std::isnan(o.beta[2]) &&
                              std::isnan(o.vb[8]));
  tk_eq(t, "cnt.alias.b1", o.beta[1], std::log(1.5), 1e-12);
  // 3. negative binomial on two groups of 6: the fitted means are the
  // group means for any theta, the variance of log(mean) is (1 / mu +
  // 1 / theta) / n, and theta maximises the likelihood
  {
    const double yn[12] = {0, 3, 9, 1, 14, 2, 5, 22, 1, 8, 30, 4};
    double xn[24];
    for (int i = 0; i < 12; ++i) {
      xn[i] = 1;
      xn[12 + i] = i >= 6 ? 1 : 0;
    }
    tk_ieq(t, "cnt.nb.rc", fit(yn, xn, 12, 2, nullptr, nullptr, 0, 1, o), CLAB_OK);
    const double m0 = 29.0 / 6, m1 = 70.0 / 6, th = o.st[2];
    tk_eq(t, "cnt.nb.b1", o.beta[1], std::log(m1 / m0), 1e-9);
    tk_eq(t, "cnt.nb.v1", o.vb[3], (1 / m0 + 1 / th + 1 / m1 + 1 / th) / 6, 1e-6);
    double mu[12];
    for (int i = 0; i < 12; ++i) mu[i] = i >= 6 ? m1 : m0;
    const double l0 = llnb(yn, mu, 12, th);
    tk_eq(t, "cnt.nb.ll", o.st[0], l0, 1e-9);
    tk_gt(t, "cnt.nb.max1", l0, llnb(yn, mu, 12, th * 1.001));
    tk_gt(t, "cnt.nb.max2", l0, llnb(yn, mu, 12, th / 1.001));
    tk_eq(t, "cnt.nb.aic", o.st[1], -2 * l0 + 6, 1e-9);
    tk_gt(t, "cnt.nb.se", o.st[3], 0);
  }
  // 4. one cluster per count: the sandwich (HC0) by hand, 2 x 2
  {
    int cl[4] = {0, 1, 2, 3};
    tk_ieq(t, "cnt.cl.rc", fit(y, x, 4, 2, off, cl, 4, 0, o), CLAB_OK);
    double mu[4], a = 0, b = 0, c = 0, m00 = 0, m01 = 0, m11 = 0;
    for (int i = 0; i < 4; ++i) {
      mu[i] = std::exp(o.beta[0] + o.beta[1] * x[4 + i] + off[i]);
      a += mu[i];
      b += mu[i] * x[4 + i];
      c += mu[i] * x[4 + i];
      const double r2 = (y[i] - mu[i]) * (y[i] - mu[i]);
      m00 += r2;
      m01 += r2 * x[4 + i];
      m11 += r2 * x[4 + i];
    }
    const double det = a * c - b * b, i01 = -b / det, i11 = a / det;
    const double v11 = i01 * (m00 * i01 + m01 * i11) + i11 * (m01 * i01 + m11 * i11);
    tk_eq(t, "cnt.cl.v11", o.vc[3], v11, 1e-10);
    // and the GEE with clusters of one: no pairs, alpha 0, the same fit
    static double wk[2 * 4 + 2 * 4 + 4 + 2 * 4 * 2];
    static int iw[2 + 4 + 1 + 4];
    double bg[2] = {o.beta[0], o.beta[1]}, vg[4], sg[4];
    const int al[2] = {0, 0};
    bg[0] += 0.2;  // a start away from the solution
    tk_ieq(t, "gee.rc1", clab_geeex(y, x, 4, 2, off, cl, 4, al, 100, 1e-13, bg, vg, sg, wk, iw),
           CLAB_OK);
    tk_eq(t, "gee.alpha0", sg[0], 0, 0);
    tk_eq(t, "gee.b1", bg[1], std::log(1.5), 1e-10);
    tk_eq(t, "gee.vc", vg[3], v11, 1e-8);
  }
  // 5. GEE on clusters of 3 and 4: the estimating equations sum_i D_i'
  // V_i^-1 (y_i - mu_i), V_i solved by elimination, vanish at the estimate
  {
    const double yg[14] = {2, 5, 3, 9, 4, 6, 1, 0, 2, 7, 12, 3, 5, 8};
    const int cg[14] = {0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3};
    double xg[28];
    for (int i = 0; i < 14; ++i) {
      xg[i] = 1;
      xg[14 + i] = (i % 3 == 0) ? 1 : 0;
    }
    tk_ieq(t, "gee.glm", fit(yg, xg, 14, 2, nullptr, cg, 4, 0, o), CLAB_OK);
    static double wk[2 * 14 + 8 + 4 + 2 * 4 * 2];
    static int iw[2 + 4 + 1 + 14];
    double bg[2] = {o.beta[0], o.beta[1]}, vg[4], sg[4];
    const int al[2] = {0, 0};
    tk_ieq(t, "gee.rc2",
           clab_geeex(yg, xg, 14, 2, nullptr, cg, 4, al, 100, 1e-13, bg, vg, sg, wk, iw),
           CLAB_OK);
    const double alpha = sg[0];
    double u0 = 0, u1 = 0, ss = 0, num = 0, cnt2 = 0;
    for (int g = 0, s = 0; g < 4; ++g) {
      int m = 0;
      while (s + m < 14 && cg[s + m] == g) ++m;
      double v[16], rhs[4], mu[4];
      for (int i = 0; i < m; ++i) {
        mu[i] = std::exp(bg[0] + bg[1] * xg[14 + s + i]);
        rhs[i] = yg[s + i] - mu[i];
        const double r = rhs[i] / std::sqrt(mu[i]);
        ss += r * r;
        for (int j = 0; j < i; ++j) num += r * (yg[s + j] - mu[j]) / std::sqrt(mu[j]);
      }
      cnt2 += m * (m - 1) / 2.0;
      for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) v[i * 4 + j] = std::sqrt(mu[i] * mu[j]) * (i == j ? 1 : alpha);
      }
      // Gaussian elimination: v w = rhs
      for (int c = 0; c < m; ++c) {
        for (int r = c + 1; r < m; ++r) {
          const double f = v[r * 4 + c] / v[c * 4 + c];
          for (int j = c; j < m; ++j) v[r * 4 + j] -= f * v[c * 4 + j];
          rhs[r] -= f * rhs[c];
        }
      }
      for (int c = m - 1; c >= 0; --c) {
        for (int j = c + 1; j < m; ++j) rhs[c] -= v[c * 4 + j] * rhs[j];
        rhs[c] /= v[c * 4 + c];
      }
      for (int i = 0; i < m; ++i) {
        u0 += mu[i] * rhs[i];
        u1 += mu[i] * xg[14 + s + i] * rhs[i];
      }
      s += m;
    }
    tk_eq(t, "gee.u0", u0, 0, 1e-9);
    tk_eq(t, "gee.u1", u1, 0, 1e-9);
    tk_eq(t, "gee.alpha", alpha, num / (cnt2 * ss / 14), 1e-9);
    tk_ieq(t, "gee.conv", static_cast<int>(sg[3]), 1);
    tk_gt(t, "gee.var", vg[3], 0);
  }
  // 6. one step is not enough: no convergence, reported as such
  tk_ieq(t, "cnt.maxit1", clab_cntreg(y, x, 4, 2, off, nullptr, 0, 0, 1, 1e-13, o.beta, o.vb,
                                      o.vc, o.st, o.alias, o.work, o.iwork),
         CLAB_OK);
  tk_ieq(t, "cnt.noconv", static_cast<int>(o.st[5]), 0);
  // 7. refusals
  tk_ieq(t, "cnt.zero", fit(y + 4, x, 2, 1, nullptr, nullptr, 0, 0, o), CLAB_EDOM);
  const double yneg[2] = {-1, 2};
  tk_ieq(t, "cnt.neg", fit(yneg, x, 2, 1, nullptr, nullptr, 0, 0, o), CLAB_EDOM);
  tk_ieq(t, "cnt.maxit", clab_cntreg(y, x, 4, 2, nullptr, nullptr, 0, 0, 0, 1e-8, o.beta, o.vb,
                                     o.vc, o.st, o.alias, o.work, o.iwork),
         CLAB_EARG);
  {
    const int cl[4] = {0, 1, 2, 9};
    tk_ieq(t, "cnt.cl.range", fit(y, x, 4, 2, off, cl, 4, 0, o), CLAB_EDOM);
    const int al[2] = {0, 0};
    double bg[2] = {0, 0}, vg[4], sg[4];
    tk_ieq(t, "gee.nocl", clab_geeex(y, x, 4, 2, off, nullptr, 4, al, 10, 1e-8, bg, vg, sg,
                                     o.work, o.iwork),
           CLAB_EDOM);
    const int all[2] = {1, 1};
    const int c2[4] = {0, 0, 1, 1};
    tk_ieq(t, "gee.none", clab_geeex(y, x, 4, 2, off, c2, 2, all, 10, 1e-8, bg, vg, sg,
                                     o.work, o.iwork),
           CLAB_EARG);
  }
}

}  // namespace

void t_cnt(tk &t, int big) {
  (void)big;
  cnt(t);
}
