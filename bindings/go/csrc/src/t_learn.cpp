// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Tests of the learners (src/learn.cpp): the ridge closed form, the
// optimality (Karush-Kuhn-Tucker) conditions of every fit on the path,
// cross-validation recomputed fold by fold, the score equations of the
// logistic regression, and every refusal. Built with the self-test.
#include <cmath>
#include <cstring>
#include <initializer_list>

#include "crimlab.h"
#include "learn.h"
#include "tchk.h"

namespace {

const int N = 120, P = 5, NL = 12;
double X[N * P], Y[N], YB[N];
double A0[NL], B[P * NL], LAM[NL], CV[2 * NL];
double W[P * P + 6 * P + 5 * N + N * P + NL * (P + 1) + 5 * (NL + 1)];
int IW[2 * N + P], BEST[2], NFIT, SCR = 1;

double rnorm(rng &g) {
  return std::sqrt(-2.0 * std::log(runi(g))) * std::cos(6.283185307179586 * runi(g));
}

// column j of X standardised over rows 0..m: mean, population sd
void cstat(const double *x, int n, int m, int j, double *mu, double *sd) {
  double s = 0.0, v = 0.0;
  for (int i = 0; i < m; ++i) s += x[i + j * n];
  *mu = s / m;
  for (int i = 0; i < m; ++i) v += (x[i + j * n] - *mu) * (x[i + j * n] - *mu);
  *sd = std::sqrt(v / m);
}

// The largest violation of the optimality conditions at lambda l: the
// gradient of the smooth part, in standardised units, against the
// subgradient of the penalty; and the intercept's score.
double kkt(const double *x, int n, int p, const double *y, int fam,
           double alpha, double l, double a0, const double *b) {
  double ym = 0.0, ys = 0.0;
  for (int i = 0; i < n; ++i) ym += y[i] / n;
  for (int i = 0; i < n; ++i) ys += (y[i] - ym) * (y[i] - ym) / n;
  ys = std::sqrt(ys);
  const double sc = fam == 0 ? ys : 1.0;  // the gaussian works on y / ys
  double r[N], worst = 0.0, s0 = 0.0;
  for (int i = 0; i < n; ++i) {
    double e = a0;
    for (int j = 0; j < p; ++j) e += x[i + j * n] * b[j];
    r[i] = fam == 0 ? (y[i] - e) / sc : y[i] - 1.0 / (1.0 + std::exp(-e));
    s0 += r[i] / n;
  }
  worst = std::fabs(s0);
  for (int j = 0; j < p; ++j) {
    double mu, sd, g = 0.0;
    cstat(x, n, n, j, &mu, &sd);
    for (int i = 0; i < n; ++i) g += (x[i + j * n] - mu) / sd * r[i] / n;
    const double bs = b[j] * sd / sc, l1 = l / sc * alpha, l2 = l / sc * (1 - alpha);
    const double v = bs != 0.0 ? std::fabs(g - l2 * bs - l1 * (bs > 0 ? 1 : -1))
                               : std::fmax(std::fabs(g) - l1, 0.0);
    worst = std::fmax(worst, v);
  }
  return worst;
}

// Gaussian elimination with partial pivoting: solve a x = b, q x q.
void gsolve(double *a, double *b, int q) {
  for (int c = 0; c < q; ++c) {
    int r = c;
    for (int i = c + 1; i < q; ++i) {
      if (std::fabs(a[i + c * q]) > std::fabs(a[r + c * q])) r = i;
    }
    for (int k = 0; k < q; ++k) {
      const double t = a[c + k * q];
      a[c + k * q] = a[r + k * q];
      a[r + k * q] = t;
    }
    const double t = b[c];
    b[c] = b[r];
    b[r] = t;
    for (int i = c + 1; i < q; ++i) {
      const double f = a[i + c * q] / a[c + c * q];
      for (int k = c; k < q; ++k) a[i + k * q] -= f * a[c + k * q];
      b[i] -= f * b[c];
    }
  }
  for (int c = q - 1; c >= 0; --c) {
    for (int k = c + 1; k < q; ++k) b[c] -= a[c + k * q] * b[k];
    b[c] /= a[c + c * q];
  }
}

int fit(int fam, double alpha, int nl, double lmr, const int *fold, int nf,
        double thresh = 1e-14, int maxit = 100000, const double *x = X,
        int n = N, int p = P, const double *y = nullptr) {
  if (y == nullptr) y = fam == 0 ? Y : YB;
  return clab_enpath(x, n, p, y, fam, alpha, LAM, nl, lmr, thresh, maxit, SCR, fold,
                    nf, A0, B, CV, BEST, &NFIT, W, IW);
}

void t_net(tk &t) {
  rng g = {11};
  const double beta[P] = {1.0, -0.5, 0.0, 0.25, 0.0};
  for (int i = 0; i < N; ++i) {
    double e = 0.0;
    for (int j = 0; j < P; ++j) {
      X[i + j * N] = (j + 1) * rnorm(g) + 2.0 * j;
      e += beta[j] * (X[i + j * N] - 2.0 * j) / (j + 1);
    }
    Y[i] = 3.0 + e + rnorm(g);
    YB[i] = runi(g) < 1.0 / (1.0 + std::exp(-e)) ? 1.0 : 0.0;
  }
  // the ridge closed form: (X~'X~ / n + (lam / ys) I) b = X~'y~ / n
  {
    LAM[0] = 0.4;
    tk_ieq(t, "net.ridge.rc", fit(0, 0.0, 1, 0.0, nullptr, 0, 1e-26), CLAB_OK);
    double mu[P], sd[P], a[P * P], r[P], ym = 0.0, ys = 0.0;
    for (int i = 0; i < N; ++i) ym += Y[i] / N;
    for (int i = 0; i < N; ++i) ys += (Y[i] - ym) * (Y[i] - ym) / N;
    ys = std::sqrt(ys);
    for (int j = 0; j < P; ++j) cstat(X, N, N, j, mu + j, sd + j);
    for (int j = 0; j < P; ++j) {
      r[j] = 0.0;
      for (int i = 0; i < N; ++i) r[j] += (X[i + j * N] - mu[j]) / sd[j] * (Y[i] - ym) / ys / N;
      for (int k = 0; k < P; ++k) {
        double s = 0.0;
        for (int i = 0; i < N; ++i) {
          s += (X[i + j * N] - mu[j]) / sd[j] * (X[i + k * N] - mu[k]) / sd[k] / N;
        }
        a[j + k * P] = s + (j == k ? 0.4 / ys : 0.0);
      }
    }
    gsolve(a, r, P);
    double want[P], a0 = ym;
    for (int j = 0; j < P; ++j) {
      want[j] = r[j] * ys / sd[j];
      a0 -= mu[j] * want[j];
    }
    tk_veq(t, "net.ridge.b", B, want, P, 1e-10);
    tk_eq(t, "net.ridge.a0", A0[0], a0, 1e-10);
  }
  // the path from the data: the first lambda zeroes every coefficient,
  // the next does not; log-even spacing down to lmr times the first
  for (int fam = 0; fam < 2; ++fam) {
    for (double alpha : {1.0, 0.5, 0.0}) {
      const int rc = fit(fam, alpha, NL, 0.01, nullptr, 0);
      tk_ieq(t, "net.path.rc", rc, CLAB_OK);
      double z0 = 0.0, z1 = 0.0, worst = 0.0;
      for (int j = 0; j < P; ++j) {
        z0 += std::fabs(B[j]);
        z1 += std::fabs(B[P + j]);
      }
      if (alpha > 0.0) tk_true(t, "net.path.zero", z0 == 0.0 && z1 > 0.0);
      tk_eq(t, "net.path.lmr", LAM[NL - 1] / LAM[0], 0.01, 1e-12);
      for (int k = 0; k < NFIT; ++k) {
        worst = std::fmax(worst, kkt(X, N, P, fam == 0 ? Y : YB, fam, alpha,
                                     LAM[k], A0[k], B + k * P));
      }
      tk_lt(t, "net.path.kkt", worst, 1e-6);
    }
  }
  // screening by the strong rule leaves the binomial solution unchanged
  {
    fit(1, 0.7, NL, 0.01, nullptr, 0, 1e-16);
    double bs[P * NL], as[NL];
    const int ns = NFIT;
    std::memcpy(bs, B, sizeof bs);
    std::memcpy(as, A0, sizeof as);
    SCR = 0;
    fit(1, 0.7, NL, 0.01, nullptr, 0, 1e-16);
    SCR = 1;
    tk_ieq(t, "net.scr.n", NFIT, ns);
    tk_veq(t, "net.scr.b", B, bs, P * ns, 1e-7);
    tk_veq(t, "net.scr.a0", A0, as, ns, 1e-7);
  }
  // the largest lambda of the binomial path is attained exactly: just
  // below it one coefficient leaves zero
  {
    fit(1, 1.0, 1, 0.5, nullptr, 0);
    const double lm = LAM[0];
    LAM[0] = lm * (1 - 1e-6);
    fit(1, 1.0, 1, 0.0, nullptr, 0);
    double z = 0.0;
    for (int j = 0; j < P; ++j) z += std::fabs(B[j]);
    tk_gt(t, "net.lmax.edge", z, 0.0);
  }
  // cross-validation, recomputed by fitting each training set on its own
  int fold[N];
  for (int i = 0; i < N; ++i) fold[i] = (i * 7) % 4;
  for (int fam = 0; fam < 2; ++fam) {
    const double *y = fam == 0 ? Y : YB;
    tk_ieq(t, "net.cv.rc", fit(fam, 0.8, NL, 0.05, fold, 4), CLAB_OK);
    const int nl = NFIT;
    double lam[NL], cm[NL], cs[NL], fsum[4 * NL], fcnt[4] = {0, 0, 0, 0};
    std::memcpy(lam, LAM, sizeof lam);
    std::memcpy(cm, CV, sizeof cm);
    std::memcpy(cs, CV + nl, sizeof cs);
    const int b0 = BEST[0], b1 = BEST[1];
    for (int f = 0; f < 4; ++f) {
      static double xt[N * P], yt[N];
      int m = 0;
      for (int i = 0; i < N; ++i) {
        if (fold[i] == f) continue;
        for (int j = 0; j < P; ++j) xt[m + j * N] = X[i + j * N];
        yt[m++] = y[i];
      }
      // xt has leading dimension N; pack it to m rows
      static double xp[N * P];
      for (int j = 0; j < P; ++j) {
        for (int i = 0; i < m; ++i) xp[i + j * m] = xt[i + j * N];
      }
      std::memcpy(LAM, lam, sizeof lam);
      fit(fam, 0.8, nl, 0.0, nullptr, 0, 1e-14, 100000, xp, m, P, yt);
      for (int k = 0; k < nl; ++k) {
        double s = 0.0;
        for (int i = 0; i < N; ++i) {
          if (fold[i] != f) continue;
          double e = A0[k];
          for (int j = 0; j < P; ++j) e += X[i + j * N] * B[k * P + j];
          if (fam == 0) {
            s += (y[i] - e) * (y[i] - e);
          } else {
            double pr = 1.0 / (1.0 + std::exp(-e));
            pr = std::fmin(std::fmax(pr, 1e-5), 1 - 1e-5);
            s += -2.0 * (y[i] == 1.0 ? std::log(pr) : std::log(1 - pr));
          }
        }
        fsum[f * NL + k] = s;
      }
      fcnt[f] = N - m;
    }
    double wm[NL], ws[NL];
    int im = 0;
    for (int k = 0; k < nl; ++k) {
      wm[k] = (fsum[k] + fsum[NL + k] + fsum[2 * NL + k] + fsum[3 * NL + k]) / N;
      double v = 0.0;
      for (int f = 0; f < 4; ++f) {
        v += fcnt[f] * (fsum[f * NL + k] / fcnt[f] - wm[k]) * (fsum[f * NL + k] / fcnt[f] - wm[k]);
      }
      ws[k] = std::sqrt(v / N / 3);
      if (wm[k] < wm[im]) im = k;
    }
    int i1 = 0;
    while (wm[i1] > wm[im] + ws[im]) ++i1;
    tk_veq(t, "net.cv.m", cm, wm, nl, 1e-9);
    tk_veq(t, "net.cv.s", cs, ws, nl, 1e-7);
    tk_ieq(t, "net.cv.min", b0, im);
    tk_ieq(t, "net.cv.1se", b1, i1);
    tk_true(t, "net.cv.order", b1 <= b0);
  }
  // a constant column keeps a zero coefficient and does not disturb the
  // others; a monotone shift of a column moves only the intercept
  {
    static double xc[N * P], yc[N];
    std::memcpy(xc, X, sizeof X);
    for (int i = 0; i < N; ++i) xc[i + 2 * N] = 7.0;
    fit(0, 1.0, NL, 0.01, nullptr, 0, 1e-14, 100000, xc);
    bool zero = true;
    for (int k = 0; k < NFIT; ++k) zero = zero && B[k * P + 2] == 0.0;
    tk_true(t, "net.const", zero);
    fit(1, 1.0, NL, 0.01, nullptr, 0, 1e-14, 100000, xc);
    zero = true;
    for (int k = 0; k < NFIT; ++k) zero = zero && B[k * P + 2] == 0.0;
    tk_true(t, "net.bconst", zero);
    std::memcpy(xc, X, sizeof X);
    for (int i = 0; i < N; ++i) xc[i] += 100.0;
    double b1[P * NL], a1[NL];
    fit(0, 0.5, NL, 0.01, nullptr, 0);
    const int n1 = NFIT;
    std::memcpy(b1, B, sizeof b1);
    std::memcpy(a1, A0, sizeof a1);
    fit(0, 0.5, NL, 0.01, nullptr, 0, 1e-14, 100000, xc);
    tk_ieq(t, "net.shift.n", NFIT, n1);
    tk_veq(t, "net.shift.b", B, b1, P * n1, 1e-9);
    tk_eq(t, "net.shift.a0", A0[n1 - 1] + 100.0 * B[n1 * P - P], a1[n1 - 1], 1e-9);
    // the early stop of a path chosen here (glmnet's rule): a gaussian
    // outcome almost exactly linear, and a binary one at penalties down to
    // 1e-9 of the largest
    for (int i = 0; i < N; ++i) yc[i] = X[i] + 1e-3 * Y[i];
    fit(0, 1.0, NL, 1e-4, nullptr, 0, 1e-14, 100000, X, N, P, yc);
    tk_true(t, "net.stop.g", NFIT >= 5 && NFIT < NL);
    fit(1, 1.0, NL, 1e-9, nullptr, 0, 1e-14, 1000000, X, N, P, YB);
    tk_true(t, "net.stop.b", NFIT >= 5 && NFIT < NL);
    // and a constant outcome or a non-binary one has no path
    for (int i = 0; i < N; ++i) yc[i] = 2.0;
    tk_ieq(t, "net.yconst", fit(0, 1.0, NL, 0.01, nullptr, 0, 1e-7, 100, X, N, P, yc), CLAB_EDOM);
    tk_ieq(t, "net.yconst2", fit(0, 1.0, 1, 0.0, nullptr, 0, 1e-7, 100, X, N, P, yc), CLAB_EDOM);
    tk_ieq(t, "net.ybin", fit(1, 1.0, NL, 0.01, nullptr, 0, 1e-7, 100, X, N, P, Y), CLAB_EDOM);
    tk_ieq(t, "net.ybin2", fit(1, 1.0, 1, 0.0, nullptr, 0, 1e-7, 100, X, N, P, Y), CLAB_EDOM);
    for (int i = 0; i < N; ++i) yc[i] = 1.0;
    tk_ieq(t, "net.ybconst", fit(1, 1.0, 1, 0.0, nullptr, 0, 1e-7, 100, X, N, P, yc), CLAB_EDOM);
    // a training set with a constant outcome: rows 0..9 hold the only ones
    for (int i = 0; i < N; ++i) yc[i] = i < 10 ? 1.0 : 0.0;
    int fd[N];
    for (int i = 0; i < N; ++i) fd[i] = i < 10 ? 0 : 1;
    tk_ieq(t, "net.foldconst", fit(1, 1.0, 3, 0.1, fd, 2, 1e-7, 100000, X, N, P, yc), CLAB_EDOM);
  }
  // refusals
  int fd[N];
  for (int i = 0; i < N; ++i) fd[i] = i % 3;
  tk_ieq(t, "net.e.n", fit(0, 1, NL, 0.01, nullptr, 0, 1e-7, 100, X, 1), CLAB_EARG);
  tk_ieq(t, "net.e.p", fit(0, 1, NL, 0.01, nullptr, 0, 1e-7, 100, X, N, 0), CLAB_EARG);
  tk_ieq(t, "net.e.nl", fit(0, 1, 0, 0.01, nullptr, 0), CLAB_EARG);
  tk_ieq(t, "net.e.it", fit(0, 1, NL, 0.01, nullptr, 0, 1e-7, 0), CLAB_EARG);
  tk_ieq(t, "net.e.fam", fit(2, 1, NL, 0.01, nullptr, 0, 1e-7, 100, X, N, P, Y), CLAB_EARG);
  tk_ieq(t, "net.e.nf", fit(0, 1, NL, 0.01, fd, 1), CLAB_EARG);
  tk_ieq(t, "net.e.nf2", fit(0, 1, NL, 0.01, fd, N + 1), CLAB_EARG);
  tk_ieq(t, "net.e.fid", fit(0, 1, NL, 0.01, fd, 2), CLAB_EARG);
  tk_ieq(t, "net.e.empty", fit(0, 1, NL, 0.01, fd, 4), CLAB_EARG);
  tk_ieq(t, "net.e.alpha", fit(0, 1.5, NL, 0.01, nullptr, 0), CLAB_EDOM);
  tk_ieq(t, "net.e.alpha2", fit(0, -0.1, NL, 0.01, nullptr, 0), CLAB_EDOM);
  tk_ieq(t, "net.e.th", fit(0, 1, NL, 0.01, nullptr, 0, 0.0), CLAB_EDOM);
  tk_ieq(t, "net.e.lmr", fit(0, 1, NL, 1.0, nullptr, 0), CLAB_EDOM);
  tk_ieq(t, "net.e.it2", fit(0, 1, NL, 0.01, nullptr, 0, 1e-14, 1), CLAB_EDOM);
  tk_ieq(t, "net.e.it3", fit(1, 1, NL, 0.01, nullptr, 0, 1e-14, 2), CLAB_EDOM);
  tk_ieq(t, "net.e.it4", fit(0, 1, NL, 0.01, fd, 3, 1e-14, 30), CLAB_EDOM);
  LAM[0] = 1.0;
  LAM[1] = 2.0;
  tk_ieq(t, "net.e.dec", fit(0, 1, 2, 0.0, nullptr, 0), CLAB_EDOM);
  LAM[0] = -1.0;
  tk_ieq(t, "net.e.neg", fit(0, 1, 1, 0.0, nullptr, 0), CLAB_EDOM);
  LAM[0] = NAN;
  tk_ieq(t, "net.e.nan", fit(0, 1, 1, 0.0, nullptr, 0), CLAB_EDOM);
  {
    static double xn[N * P], yn[N];
    std::memcpy(xn, X, sizeof X);
    std::memcpy(yn, Y, sizeof Y);
    xn[7] = NAN;
    tk_ieq(t, "net.e.x", fit(0, 1, NL, 0.01, nullptr, 0, 1e-7, 100, xn), CLAB_EDOM);
    yn[3] = INFINITY;
    tk_ieq(t, "net.e.y", fit(0, 1, NL, 0.01, nullptr, 0, 1e-7, 100, X, N, P, yn), CLAB_EDOM);
  }
}

void t_glm(tk &t) {
  // logistic regression: the score equations hold at the fit, and the
  // standard errors are the root diagonal of the inverse information
  const int q = P + 1;
  double cf[q], se[q], dev, w[2 * q * q + 4 * q];
  int iw[N];
  tk_ieq(t, "glm.rc", clab_lgtreg(X, N, P, YB, 1e-12, 50, cf, se, &dev, w, iw), CLAB_OK);
  double u[q] = {0}, info[q * q] = {0}, d = 0.0;
  for (int i = 0; i < N; ++i) {
    double z[q] = {1.0}, e = cf[0];
    for (int j = 0; j < P; ++j) {
      z[j + 1] = X[i + j * N];
      e += cf[j + 1] * z[j + 1];
    }
    const double pr = 1.0 / (1.0 + std::exp(-e));
    d -= 2.0 * (YB[i] == 1.0 ? std::log(pr) : std::log(1 - pr));
    for (int c = 0; c < q; ++c) {
      u[c] += z[c] * (YB[i] - pr);
      for (int k = 0; k < q; ++k) info[c + k * q] += pr * (1 - pr) * z[c] * z[k];
    }
  }
  double umax = 0.0;
  for (int c = 0; c < q; ++c) umax = std::fmax(umax, std::fabs(u[c]));
  tk_lt(t, "glm.score", umax, 1e-9);
  tk_eq(t, "glm.dev", dev, d, 1e-12);
  for (int c = 0; c < q; ++c) {
    double a[q * q], e[q] = {0};
    std::memcpy(a, info, sizeof a);
    e[c] = 1.0;
    gsolve(a, e, q);
    tk_eq(t, "glm.se", se[c], std::sqrt(e[c]), 1e-9);
  }
  // the gaussian family: least squares, the residual variance in the
  // standard errors; checked against the normal equations
  int idx[N];
  for (int i = 0; i < N; ++i) idx[i] = i;
  const lrn::Dat dd = {X, N, P, Y, idx, N};
  tk_ieq(t, "glm.g.rc", lrn::glm(dd, 0, 1e-12, 5, cf, se, &dev, w), CLAB_OK);
  double xty[q] = {0}, xtx[q * q] = {0};
  for (int i = 0; i < N; ++i) {
    double z[q] = {1.0};
    for (int j = 0; j < P; ++j) z[j + 1] = X[i + j * N];
    for (int c = 0; c < q; ++c) {
      xty[c] += z[c] * Y[i];
      for (int k = 0; k < q; ++k) xtx[c + k * q] += z[c] * z[k];
    }
  }
  double a[q * q];
  std::memcpy(a, xtx, sizeof a);
  gsolve(a, xty, q);
  tk_veq(t, "glm.g.coef", cf, xty, q, 1e-10);
  double rss = 0.0;
  for (int i = 0; i < N; ++i) {
    double e = xty[0];
    for (int j = 0; j < P; ++j) e += xty[j + 1] * X[i + j * N];
    rss += (Y[i] - e) * (Y[i] - e);
  }
  tk_eq(t, "glm.g.dev", dev, rss, 1e-10);
  double e1[q] = {0, 1};
  std::memcpy(a, xtx, sizeof a);
  gsolve(a, e1, q);
  tk_eq(t, "glm.g.se", se[1], std::sqrt(rss / (N - q) * e1[1]), 1e-9);
  // prediction
  double eta[3];
  const double b[P] = {1, 0, 2, 0, 0};
  lrn::lpred(X, N, P, idx, 3, 0.5, b, eta);
  tk_eq(t, "glm.pred", eta[2], 0.5 + X[2] + 2 * X[2 + 2 * N], 1e-15);
  // refusals: separation, collinearity, too few rows, bad arguments
  double ys[N], xc[N * 2];
  for (int i = 0; i < N; ++i) {
    ys[i] = X[i] > 0.0 ? 1.0 : 0.0;
    xc[i] = X[i];
    xc[i + N] = 2.0 * X[i];
  }
  tk_ieq(t, "glm.e.sep", clab_lgtreg(X, N, 1, ys, 1e-12, 50, cf, se, &dev, w, iw), CLAB_EDOM);
  tk_ieq(t, "glm.e.col", clab_lgtreg(xc, N, 2, YB, 1e-12, 50, cf, se, &dev, w, iw), CLAB_EDOM);
  tk_ieq(t, "glm.e.m", clab_lgtreg(X, 2, 1, YB, 1e-12, 50, cf, se, &dev, w, iw), CLAB_EDOM);
  tk_ieq(t, "glm.e.y", clab_lgtreg(X, N, 1, Y, 1e-12, 50, cf, se, &dev, w, iw), CLAB_EDOM);
  for (int i = 0; i < N; ++i) ys[i] = 0.0;
  tk_ieq(t, "glm.e.yc", clab_lgtreg(X, N, 1, ys, 1e-12, 50, cf, se, &dev, w, iw), CLAB_EDOM);
  tk_ieq(t, "glm.e.tol", clab_lgtreg(X, N, 1, YB, 0.0, 50, cf, se, &dev, w, iw), CLAB_EDOM);
  tk_ieq(t, "glm.e.n", clab_lgtreg(X, 0, 1, YB, 1e-9, 50, cf, se, &dev, w, iw), CLAB_EARG);
  tk_ieq(t, "glm.e.it", clab_lgtreg(X, N, 1, YB, 1e-9, 0, cf, se, &dev, w, iw), CLAB_EARG);
  tk_ieq(t, "glm.e.it1", clab_lgtreg(X, N, 1, YB, 1e-9, 1, cf, se, &dev, w, iw), CLAB_EDOM);
  xc[3] = NAN;
  tk_ieq(t, "glm.e.x", clab_lgtreg(xc, N, 1, YB, 1e-9, 50, cf, se, &dev, w, iw), CLAB_EDOM);
}

}  // namespace

void t_learn(tk &t, int) {
  // the scratch the sizes promise is what the tests hand over
  tk_true(t, "net.size", lrn::cvwork(N, P, NL, 5) <= static_cast<long>(sizeof W / sizeof W[0]) &&
                             lrn::gwork(P) == 2 * (P + 1) * (P + 1) + 4 * (P + 1));
  t_net(t);
  t_glm(t);
}
