// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, causal: double machine learning by cross-fitting
// (Chernozhukov et al. 2018, doi:10.1111/ectj.12097). One engine fits
// every nuisance function fold by fold with the learners of learn.cpp and
// writes the out-of-fold predictions; the estimators are then closed forms
// of those predictions: the linear Neyman-orthogonal scores of the
// partially linear models, the interactive models, difference in
// differences and sample selection, the augmented inverse probability
// weighting estimator and targeted maximum likelihood. Nothing is
// allocated; no R here.
#include <cmath>
#include <initializer_list>

#include "crimlab.h"
#include "kutil.h"
#include "learn.h"

namespace {

inline double expit(double e) { return 1.0 / (1.0 + std::exp(-e)); }
inline double clip(double v, double lo, double hi) {
  return v < lo ? lo : (v > hi ? hi : v);
}

// Scratch of one learner: doubles and ints, for a design of n rows and p
// columns.
long lw(const clab_lrn &l, int n, int p) {
  if (l.kind == 0) return 2L * (p + 1) + lrn::gwork(p) + n;
  if (l.kind == 2) return lrn::fwork(n, p, n) + n;
  return static_cast<long>(l.nlam) * (p + 4) + lrn::cvwork(n, p, l.nlam, l.nfold) + n;
}
long liw(const clab_lrn &l, int n, int p) {
  if (l.kind == 0) return 0;
  if (l.kind == 2) return lrn::fiwork(n, p, n);
  return 3L * n + p;
}

// Fit learner l to y over the training rows tr[0..mt) and write its
// predictions at the test rows te[0..me) to out[te[i]]: the mean of y, or
// the probability for fam 1. A training outcome that is constant is its
// own prediction.
int nfit(const clab_lrn &l, int fam, const double *x, int n, int p,
         const double *y, const int *tr, int mt, const int *te, int me,
         unsigned int seed, double *out, double *w, int *iw) {
  if (mt < 1) return CLAB_EDOM;
  bool cst = true;
  for (int i = 1; i < mt && cst; ++i) cst = y[tr[i]] == y[tr[0]];
  if (cst) {
    for (int i = 0; i < me; ++i) out[te[i]] = y[tr[0]];
    return CLAB_OK;
  }
  const lrn::Dat d = {x, n, p, y, tr, mt};
  double *eta = w;
  if (l.kind == 2) {
    // the forest's mean is already a probability for a 0/1 target; the
    // smallest node left unsplit defaults to ranger's 5 and 10
    const int mn = l.minnode > 0 ? l.minnode : (fam == 1 ? 10 : 5);
    const int mtry = l.mtry > 0 ? (l.mtry < p ? l.mtry : p)
                                : static_cast<int>(std::fmin(std::ceil(std::sqrt(p) + 20.0), p));
    const lrn::Forest f = {l.ntree, mtry, mn, l.hon ? 0.5 : 1.0, l.hon ? 0 : 1, seed, l.hon,
                           l.hon ? 0.05 : 0.0};
    // fails only when interrupted: ntree is checked by clab_dmlx(), mtry
    // is in range and the training set has two rows or more (one row is a
    // constant)
    const int rc = lrn::forest(x, n, p, y, tr, mt, x, n, te, me, f, eta, nullptr, w + n, iw);
    if (rc != CLAB_OK) return rc;
    for (int i = 0; i < me; ++i) out[te[i]] = eta[i];
    return CLAB_OK;
  }
  if (l.kind == 0) {
    double *cf = eta + n, *se = cf + p + 1, dev;
    const int rc = lrn::glm(d, fam, 1e-10, 50, cf, se, &dev, se + p + 1);
    if (rc != CLAB_OK) return rc;
    lrn::lpred(x, n, p, te, me, cf[0], cf + 1, eta);
  } else {
    // the inner folds: a seeded shuffle of the training rows dealt round
    // robin, so each fold differs in size by at most one
    const int nl = l.nlam, k = l.nfold;
    int *f = iw, *pm = f + mt, *ciw = pm + mt;
    kutil::rng g(seed);
    for (int i = 0; i < mt; ++i) pm[i] = i;
    for (int i = mt - 1; i > 0; --i) {
      const int j = g.below(i + 1), t = pm[i];
      pm[i] = pm[j];
      pm[j] = t;
    }
    for (int i = 0; i < mt; ++i) f[pm[i]] = i % k;
    double *lam = eta + n, *a0 = lam + nl, *cv = a0 + nl, *beta = cv + 2 * nl;
    double *cw = beta + static_cast<long>(p) * nl;
    int best[2];
    const double lmr = mt > p ? 1e-4 : 1e-2;
    int nfit;
    const int rc = lrn::cvnet(d, fam, l.alpha, lam, nl, lmr, l.tol, 100000, l.scr, f, k,
                              a0, beta, cv, best, &nfit, cw, ciw);
    if (rc != CLAB_OK) return rc;
    lrn::lpred(x, n, p, te, me, a0[best[0]], beta + static_cast<long>(p) * best[0], eta);
  }
  for (int i = 0; i < me; ++i) out[te[i]] = fam == 1 ? expit(eta[i]) : eta[i];
  return CLAB_OK;
}

// The nuisances of each model: what is learned (the outcome y, the
// treatment d, the instrument z, the selection s, or a working target),
// with which family and learner, on which training rows.
enum Tgt { TY, TD, TZ, TS, TW };
struct Nu {
  Tgt tgt;
  int fam;    // 0 regression, 1 probability
  int out;    // 1: the outcome learner, 0: the treatment learner
  int cond;   // the training rows: see keep()
  int aug;    // 1: the design is x with d appended
};
// training-row conditions: 0 all; 1 d == 0; 2 d == 1; 3 z == 0; 4 z == 1;
// 5..8 (d, t) = (0, 0), (0, 1), (1, 0), (1, 1); 9 d == 0 and s == 1;
// 10 d == 1 and s == 1
bool keep(int c, int i, const double *d, const double *z, const double *s) {
  switch (c) {
    case 1: return d[i] == 0.0;
    case 2: return d[i] == 1.0;
    case 3: return z[i] == 0.0;
    case 4: return z[i] == 1.0;
    case 5: return d[i] == 0.0 && s[i] == 0.0;
    case 6: return d[i] == 0.0 && s[i] == 1.0;
    case 7: return d[i] == 1.0 && s[i] == 0.0;
    case 8: return d[i] == 1.0 && s[i] == 1.0;
    case 9: return d[i] == 0.0 && s[i] == 1.0;
    case 10: return d[i] == 1.0 && s[i] == 1.0;
    default: return true;
  }
}

const Nu kplr[] = {{TY, 0, 1, 0, 0}, {TD, 0, 0, 0, 0}, {TW, 0, 1, 0, 0}};
const Nu kpliv[] = {{TY, 0, 1, 0, 0}, {TZ, 0, 0, 0, 0}, {TD, 0, 0, 0, 0},
                    {TW, 0, 1, 0, 0}};
const Nu kirm[] = {{TY, 0, 1, 1, 0}, {TY, 0, 1, 2, 0}, {TD, 1, 0, 0, 0}};
const Nu kiivm[] = {{TY, 0, 1, 3, 0}, {TY, 0, 1, 4, 0}, {TZ, 1, 0, 0, 0},
                    {TD, 1, 0, 3, 0}, {TD, 1, 0, 4, 0}};
const Nu kdid[] = {{TY, 0, 1, 1, 0}, {TD, 1, 0, 0, 0}};
const Nu kdidcs[] = {{TY, 0, 1, 5, 0}, {TY, 0, 1, 6, 0}, {TY, 0, 1, 7, 0},
                     {TY, 0, 1, 8, 0}, {TD, 1, 0, 0, 0}};
const Nu kssm[] = {{TS, 1, 0, 0, 1}, {TD, 1, 0, 0, 0}, {TY, 0, 1, 9, 0},
                   {TY, 0, 1, 10, 0}};

// the experimental scores of difference in differences need no propensity
const Nu kdidx[] = {{TY, 0, 1, 1, 0}, {TY, 0, 1, 2, 0}};

const Nu *spec(int model, int sc) {
  switch (model) {
    case 0: return kplr;
    case 1: return kpliv;
    case 6: return sc == 1 ? kdidx : kdid;
    case 7: return kdidcs;
    case 8: return kssm;
    case 3: return kiivm;
    default: return kirm;  // irm, aipw, tmle
  }
}

// -sum(psi_b) / sum(psi_a): the root of the mean of a linear score
double root(const double *pa, const double *pb, int n) {
  double sa = 0.0, sb = 0.0;
  for (int i = 0; i < n; ++i) {
    sa += pa[i];
    sb += pb[i];
  }
  return -sb / sa;
}

// The score elements psi_a, psi_b of one repetition from its predictions
// (pr: k columns of n), and theta and its variance as DoubleML computes
// them (dml2): theta = -mean(psi_b) / mean(psi_a), var = mean(psi^2) /
// mean(psi_a)^2 / n; clustered (one-way, fold constant within cluster) as
// Chiang et al. (2022, doi:10.1080/07350015.2021.1895815). The variance of
// aipw and tmle is var(IC) / n with the n - 1 divisor, as the AIPW and
// tmle packages compute it.
// Two-way clusters (Chiang et al. 2022; DoubleML's two-way form): rows
// carry clusters c1 = clu[i], c2 = clu[n + i] and the folds of their
// clusters, f1 = fold[i], f2 = fold[n + i]; the test set of block (k, l)
// is f1 = k and f2 = l. With I_k, J_l the clusters of each kind in fold k
// or l, theta = -sum_b sum psi_b / (|I_k| |J_l|) over the same with psi_a;
// gamma sums, block by block, min(|I_k|, |J_l|) / (|I_k| |J_l|)^2 times
// the squared sums of psi by c1 over rows with f2 = l and by c2 over rows
// with f1 = k; var = gamma / K^2 / (j / K^2)^2 / min(C1, C2).
int twoway(const double *pa, const double *pb, int n, const int *fold, int nf,
           const int *clu, int nc, int nc2, double *th, double *var, double *w) {
  const int *c2 = clu + n, *f2 = fold + n;
  double *k1 = w, *k2 = k1 + nf, *g1 = k2 + nf, *g2 = g1 + nc, *s1 = g2 + nc2,
         *s2 = s1 + static_cast<long>(nc) * nf;
  for (int k = 0; k < nf; ++k) k1[k] = k2[k] = 0.0;
  for (int c = 0; c < nc; ++c) g1[c] = -1.0;
  for (int c = 0; c < nc2; ++c) g2[c] = -1.0;
  for (int i = 0; i < n; ++i) {
    if (g1[clu[i]] < 0.0) {
      g1[clu[i]] = fold[i];
      k1[fold[i]] += 1.0;
    }
    if (g2[c2[i]] < 0.0) {
      g2[c2[i]] = f2[i];
      k2[f2[i]] += 1.0;
    }
  }
  double sa = 0.0, sb = 0.0;
  for (int i = 0; i < n; ++i) {
    const double q = k1[fold[i]] * k2[f2[i]];
    sa += pa[i] / q;
    sb += pb[i] / q;
  }
  *th = -sb / sa;
  if (!std::isfinite(*th)) return CLAB_EDOM;
  for (long c = 0; c < static_cast<long>(nc) * nf; ++c) s1[c] = 0.0;
  for (long c = 0; c < static_cast<long>(nc2) * nf; ++c) s2[c] = 0.0;
  for (int i = 0; i < n; ++i) {
    const double ps = pa[i] * *th + pb[i];
    s1[static_cast<long>(clu[i]) * nf + f2[i]] += ps;
    s2[static_cast<long>(c2[i]) * nf + fold[i]] += ps;
  }
  double gam = 0.0;
  for (int c = 0; c < nc; ++c) {
    const int k = static_cast<int>(g1[c]);
    for (int l = 0; l < nf; ++l) {
      const double q = k1[k] * k2[l], v = s1[static_cast<long>(c) * nf + l];
      gam += std::fmin(k1[k], k2[l]) / (q * q) * v * v;
    }
  }
  for (int c = 0; c < nc2; ++c) {
    const int l = static_cast<int>(g2[c]);
    for (int k = 0; k < nf; ++k) {
      const double q = k1[k] * k2[l], v = s2[static_cast<long>(c) * nf + k];
      gam += std::fmin(k1[k], k2[l]) / (q * q) * v * v;
    }
  }
  const double kk = static_cast<double>(nf) * nf;
  gam /= kk;
  const double j = sa / kk;
  *var = gam / (j * j) / std::fmin(nc, nc2);
  return CLAB_OK;
}

int score(int model, int sc, const double *y, const double *d,
          const double *z, const double *s, const double *pr, int n,
          const int *fold, int nf, const int *clu, int nc, int nc2,
          double trim, double *pa, double *pb, double *th, double *var,
          double *w) {
  const double *p0 = pr, *p1 = pr + n, *p2 = pr + 2L * n, *p3 = pr + 3L * n,
               *p4 = pr + 4L * n;
  switch (model) {
    case 0:
      for (int i = 0; i < n; ++i) {
        const double v = d[i] - p1[i];
        pa[i] = sc == 0 ? -v * v : -v * d[i];
        pb[i] = sc == 0 ? v * (y[i] - p0[i]) : v * (y[i] - p2[i]);
      }
      break;
    case 1:
      for (int i = 0; i < n; ++i) {
        const double v = z[i] - p1[i];
        pa[i] = sc == 0 ? -(d[i] - p2[i]) * v : -d[i] * v;
        pb[i] = sc == 0 ? v * (y[i] - p0[i]) : v * (y[i] - p3[i]);
      }
      break;
    case 2:
    case 4: {
      // ATTE: the treated share of each test fold
      double *ph = w, *cnt = w + nf;
      for (int k = 0; k < nf; ++k) ph[k] = cnt[k] = 0.0;
      for (int i = 0; i < n; ++i) {
        ph[fold[i]] += d[i];
        cnt[fold[i]] += 1.0;
      }
      for (int i = 0; i < n; ++i) {
        const double m = clip(p2[i], trim, 1.0 - trim);
        if (sc == 0) {
          pa[i] = -1.0;
          pb[i] = p1[i] - p0[i] + d[i] * (y[i] - p1[i]) / m -
                  (1.0 - d[i]) * (y[i] - p0[i]) / (1.0 - m);
        } else {
          const double q = ph[fold[i]] / cnt[fold[i]];
          pa[i] = -d[i] / q;
          pb[i] = d[i] * (y[i] - p0[i]) / q -
                  m * (1.0 - d[i]) * (y[i] - p0[i]) / (q * (1.0 - m));
        }
      }
      break;
    }
    case 3:
      for (int i = 0; i < n; ++i) {
        const double m = clip(p2[i], trim, 1.0 - trim);
        pb[i] = p1[i] - p0[i] + z[i] * (y[i] - p1[i]) / m -
                (1.0 - z[i]) * (y[i] - p0[i]) / (1.0 - m);
        pa[i] = -(p4[i] - p3[i] + z[i] * (d[i] - p4[i]) / m -
                  (1.0 - z[i]) * (d[i] - p3[i]) / (1.0 - m));
      }
      break;
    case 6: {
      if (sc == 1) {
        // the experimental score (treatment at random): p0 = g0, p1 = g1
        double md = 0.0;
        for (int i = 0; i < n; ++i) md += d[i] / n;
        if (!(md > 0.0 && md < 1.0)) return CLAB_EDOM;
        for (int i = 0; i < n; ++i) {
          const double q = d[i] / md;
          pa[i] = -1.0;
          pb[i] = (q - 1.0) * p0[i] + (1.0 - q) * p1[i] +
                  (q - (1.0 - d[i]) / (1.0 - md)) * (y[i] - p0[i]);
        }
        break;
      }
      // in-sample normalised weights (DoubleML for Python's default)
      double md = 0.0, mw = 0.0;
      for (int i = 0; i < n; ++i) {
        const double m = clip(p1[i], trim, 1.0 - trim);
        md += d[i] / n;
        mw += (1.0 - d[i]) * m / (1.0 - m) / n;
      }
      if (!(md > 0.0 && mw > 0.0)) return CLAB_EDOM;
      for (int i = 0; i < n; ++i) {
        const double m = clip(p1[i], trim, 1.0 - trim);
        pa[i] = -d[i] / md;
        pb[i] = (d[i] / md - (1.0 - d[i]) * m / (1.0 - m) / mw) * (y[i] - p0[i]);
      }
      break;
    }
    case 7: {
      // t is passed in s; p0..p3: g(d, t) at (0,0), (0,1), (1,0), (1,1)
      if (sc == 1) {
        // experimental: cell shares in place of propensities
        double c[4] = {0.0, 0.0, 0.0, 0.0};
        for (int i = 0; i < n; ++i) c[static_cast<int>(2 * d[i] + s[i])] += 1.0 / n;
        if (!(c[0] > 0.0 && c[1] > 0.0 && c[2] > 0.0 && c[3] > 0.0)) return CLAB_EDOM;
        for (int i = 0; i < n; ++i) {
          const int cell = static_cast<int>(2 * d[i] + s[i]);
          const double *g = cell == 0 ? p0 : (cell == 1 ? p1 : (cell == 2 ? p2 : p3));
          const double sg = cell == 0 || cell == 3 ? 1.0 : -1.0;
          pa[i] = -1.0;
          pb[i] = p3[i] - p2[i] - p1[i] + p0[i] + sg * (y[i] - g[i]) / c[cell];
        }
        break;
      }
      double md = 0.0, m11 = 0.0, m10 = 0.0, m01 = 0.0, m00 = 0.0;
      for (int i = 0; i < n; ++i) {
        const double m = clip(p4[i], trim, 1.0 - trim), r = m / (1.0 - m);
        md += d[i];
        m11 += d[i] * s[i];
        m10 += d[i] * (1.0 - s[i]);
        m01 += (1.0 - d[i]) * s[i] * r;
        m00 += (1.0 - d[i]) * (1.0 - s[i]) * r;
      }
      if (!(md > 0.0 && m11 > 0.0 && m10 > 0.0 && m01 > 0.0 && m00 > 0.0))
        return CLAB_EDOM;
      md /= n;
      m11 /= n;
      m10 /= n;
      m01 /= n;
      m00 /= n;
      for (int i = 0; i < n; ++i) {
        const double m = clip(p4[i], trim, 1.0 - trim), r = m / (1.0 - m);
        const double wa = d[i] / md;
        pa[i] = -wa;
        pb[i] = wa * (p3[i] - p2[i] - p1[i] + p0[i]) +
                d[i] * s[i] / m11 * (y[i] - p3[i]) -
                d[i] * (1.0 - s[i]) / m10 * (y[i] - p2[i]) +
                (1.0 - d[i]) * (1.0 - s[i]) * r / m00 * (y[i] - p0[i]) -
                (1.0 - d[i]) * s[i] * r / m01 * (y[i] - p1[i]);
      }
      break;
    }
    case 8:
      // p0 = P(s = 1 | d, x), p1 = P(d = 1 | x), p2, p3: E[y | d, s = 1, x]
      for (int i = 0; i < n; ++i) {
        const double m = clip(p1[i], trim, 1.0 - trim);
        pa[i] = -1.0;
        const double a = s[i] == 1.0 ? (y[i] - p3[i]) : 0.0;
        const double b = s[i] == 1.0 ? (y[i] - p2[i]) : 0.0;
        pb[i] = d[i] * s[i] * a / (m * p0[i]) + p3[i] -
                (1.0 - d[i]) * s[i] * b / ((1.0 - m) * p0[i]) - p2[i];
      }
      break;
    default: {
      // tmle: the fluctuation of the tmle package. Y scaled to [0, 1] by
      // its range; Q to the range widened by 1 per cent, then scaled and
      // kept in [5e-4, 1 - 5e-4]; g bounded below by trim; one logistic
      // fluctuation per arm (the two clever covariates have disjoint
      // support), solved by Newton's method.
      double a = y[0], b = y[0];
      for (int i = 1; i < n; ++i) {
        a = std::fmin(a, y[i]);
        b = std::fmax(b, y[i]);
      }
      const double ab = b - a, qlo = a - 0.01 * std::fabs(a), qhi = b + 0.01 * std::fabs(b);
      auto lq = [&](double q) {
        const double u = clip((clip(q, qlo, qhi) - a) / ab, 5e-4, 1.0 - 5e-4);
        return std::log(u / (1.0 - u));
      };
      double eps[2] = {0.0, 0.0};
      for (int arm = 0; arm < 2; ++arm) {
        const double *q = arm == 0 ? p0 : p1;
        // the score in eps is decreasing; Newton from 0, 100 steps at most
        for (int it = 0; it < 100; ++it) {
          double u = 0.0, h = 0.0;
          for (int i = 0; i < n; ++i) {
            if (d[i] != arm) continue;
            const double g = std::fmax(arm == 1 ? p2[i] : 1.0 - p2[i], trim);
            const double c = 1.0 / g, e = expit(lq(q[i]) + eps[arm] * c);
            u += c * ((y[i] - a) / ab - e);
            h += c * c * e * (1.0 - e);
          }
          const double st = u / h;
          eps[arm] += st;
          if (std::fabs(st) <= 1e-13 * (1.0 + std::fabs(eps[arm]))) break;
        }
      }
      double m1 = 0.0, m0 = 0.0;
      for (int i = 0; i < n; ++i) {
        const double g1 = std::fmax(p2[i], trim), g0 = std::fmax(1.0 - p2[i], trim);
        const double q1 = a + ab * expit(lq(p1[i]) + eps[1] / g1);
        const double q0 = a + ab * expit(lq(p0[i]) + eps[0] / g0);
        const double qa = d[i] == 1.0 ? q1 : q0;
        pa[i] = -1.0;
        // the influence function, plus the estimate (added below)
        pb[i] = (d[i] / g1 - (1.0 - d[i]) / g0) * (y[i] - qa) + q1 - q0;
        m1 += q1 / n;
        m0 += q0 / n;
      }
      *th = m1 - m0;
      double v = 0.0, mu = 0.0;
      for (int i = 0; i < n; ++i) mu += (pb[i] - *th) / n;
      for (int i = 0; i < n; ++i) v += (pb[i] - *th - mu) * (pb[i] - *th - mu);
      // psi_b = IC + theta, so that theta = -mean(psi_b) / mean(psi_a)
      // would give mean(IC) + theta: the IC has mean 0 only up to the
      // fluctuation's tolerance; theta is the plug-in mean, as tmle
      for (int i = 0; i < n; ++i) pb[i] -= mu;
      *var = v / (n - 1) / n;
      return CLAB_OK;
    }
  }
  if (clu == nullptr || model > 1) {
    *th = root(pa, pb, n);
    if (!std::isfinite(*th)) return CLAB_EDOM;
    if (model == 4) {
      double v = 0.0;
      for (int i = 0; i < n; ++i) v += (pb[i] - *th) * (pb[i] - *th);
      *var = v / (n - 1) / n;
      return CLAB_OK;
    }
    double sa = 0.0, sp = 0.0;
    for (int i = 0; i < n; ++i) {
      const double ps = pa[i] * *th + pb[i];
      sa += pa[i];
      sp += ps * ps;
    }
    sa /= n;
    *var = sp / n / (sa * sa) / n;
    return CLAB_OK;
  }
  if (nc2 > 0) return twoway(pa, pb, n, fold, nf, clu, nc, nc2, th, var, w);
  // one-way clusters: each fold's clusters counted, sums scaled by them
  double *ck = w, *cs = w + nf;
  for (int k = 0; k < nf; ++k) ck[k] = 0.0;
  for (int c = 0; c < nc; ++c) cs[c] = -1.0;
  for (int i = 0; i < n; ++i) {
    if (cs[clu[i]] < 0.0) {
      cs[clu[i]] = fold[i];
      ck[fold[i]] += 1.0;
    }
  }
  double sa = 0.0, sb = 0.0;
  for (int i = 0; i < n; ++i) {
    sa += pa[i] / ck[fold[i]];
    sb += pb[i] / ck[fold[i]];
  }
  *th = -sb / sa;
  if (!std::isfinite(*th)) return CLAB_EDOM;
  // gamma: per-cluster sums of psi, squared, over the fold's cluster count
  for (int c = 0; c < nc; ++c) cs[c] = 0.0;
  for (int i = 0; i < n; ++i) cs[clu[i]] += pa[i] * *th + pb[i];
  double gam = 0.0;
  // the fold of each cluster, found again from its first row
  double *cf = cs + nc;
  for (int c = 0; c < nc; ++c) cf[c] = -1.0;
  for (int i = 0; i < n; ++i) {
    if (cf[clu[i]] < 0.0) cf[clu[i]] = fold[i];
  }
  for (int c = 0; c < nc; ++c) {
    gam += cs[c] * cs[c] / ck[static_cast<int>(cf[c])];
  }
  gam /= nf;
  const double j = sa / nf;
  *var = gam / (j * j) / nc;
  return CLAB_OK;
}

}  // namespace

extern "C" {

long clab_dmlxsz(int model, int n, int p, int nf, int nrep,
                const clab_lrn *ly, const clab_lrn *ld, long *niw) {
  const int pa = model == 8 ? p + 1 : p;
  long a = lw(*ly, n, pa), b = lw(*ld, n, pa);
  long ia = liw(*ly, n, pa), ib = liw(*ld, n, pa);
  *niw = 2L * n + (ia > ib ? ia : ib);
  const long m = (a > b ? a : b) + 2L * n + nf + (model == 8 ? static_cast<long>(n) * (p + 1) : 0) +
                 2L * nf + 2L * n * (nf + 1);
  return m > nrep ? m : nrep;
}

int clab_dmlx(int model, int sc, const double *x, int n, int p,
             const double *y, const double *d, const double *z,
             const double *s, const int *fold, int nf, int nrep,
             const int *clu, int nc, int nc2, const clab_lrn *ly,
             const clab_lrn *ld, double trim, int fit, double *pred,
             double *est, double *psi, double *work, int *iwork) {
  if (model < 0 || model > 8 || sc < 0 || sc > 1 || n < 2 || p < 1 ||
      nf < 2 || nf > n || nrep < 1 || (clu != nullptr && (model > 1 || nc < 2)) ||
      (nc2 != 0 && (clu == nullptr || nc2 < 2 || nrep != 1)))
    return CLAB_EARG;
  const int nw = nc2 > 0 ? 2 : 1;  // two-way: fold and clu hold two columns
  const int k = clab_dmlxk(model, sc);
  for (long i = 0; i < static_cast<long>(n) * nrep * nw; ++i) {
    if (fold[i] < 0 || fold[i] >= nf) return CLAB_EARG;
  }
  if (clu != nullptr) {
    for (long i = 0; i < static_cast<long>(n) * nw; ++i) {
      if (clu[i] < 0 || clu[i] >= (i < n ? nc : nc2)) return CLAB_EARG;
    }
  }
  for (const clab_lrn *l : {ly, ld}) {
    if (l->kind < 0 || l->kind > 2 || (l->kind == 1 && (l->nlam < 1 || l->nfold < 2)) ||
        (l->kind == 2 && l->ntree < 1))
      return CLAB_EARG;
    if (!(l->alpha >= 0.0 && l->alpha <= 1.0) || (l->kind == 1 && !(l->tol > 0.0)))
      return CLAB_EDOM;
  }
  if (!(trim >= 0.0 && trim < 0.5)) return CLAB_EDOM;
  const Nu *nu = spec(model, sc);
  int *tr = iwork, *te = tr + n, *liwk = te + n;
  double *tw = work + 2L * n + nf, *xa = nullptr;
  double *tgt = work;  // the working target of the IV-type score
  if (model == 8) {
    xa = tw;
    tw = xa + static_cast<long>(n) * (p + 1);
    for (long i = 0; i < static_cast<long>(n) * p; ++i) xa[i] = x[i];
    for (int i = 0; i < n; ++i) xa[static_cast<long>(n) * p + i] = d[i];
  }
  double *th = est + 2, *se = est + 2 + nrep;
  for (int r = 0; r < nrep; ++r) {
    const int *fr = fold + static_cast<long>(r) * n;
    double *pr = pred + static_cast<long>(r) * k * n;
    for (int j = 0; fit && j < k; ++j) {
      const Nu &u = nu[j];
      const clab_lrn &l = u.out ? *ly : *ld;
      const double *tv = u.tgt == TY ? y : u.tgt == TD ? d : u.tgt == TZ ? z : u.tgt == TS ? s : tgt;
      if (u.tgt == TW) {
        // y - theta0 d, theta0 from the partialling-out score
        const double *pd = pr + (model == 0 ? 1L : 2L) * n;
        const double *pz = pr + n;
        double sa = 0.0, sb = 0.0;
        for (int i = 0; i < n; ++i) {
          const double v = model == 0 ? d[i] - pd[i] : z[i] - pz[i];
          const double w = model == 0 ? v : d[i] - pd[i];
          sa += v * w;
          sb += v * (y[i] - pr[i]);
        }
        const double t0 = sb / sa;
        for (int i = 0; i < n; ++i) tgt[i] = y[i] - t0 * d[i];
      }
      // one block per fold; two-way, one per pair of folds, whose training
      // rows share neither cluster with its test rows
      for (int f = 0; f < nf * (nw == 2 ? nf : 1); ++f) {
        const int f1 = nw == 2 ? f / nf : f, f2 = f % nf;
        int mt = 0, me = 0;
        for (int i = 0; i < n; ++i) {
          const bool in = fr[i] == f1 && (nw == 1 || fr[n + i] == f2);
          if (in) {
            te[me++] = i;
          } else if ((nw == 1 || (fr[i] != f1 && fr[n + i] != f2)) &&
                     keep(u.cond, i, d, z, s)) {
            tr[mt++] = i;
          }
        }
        const unsigned int sd = l.seed + 7919U * static_cast<unsigned int>(r * 64 + j * 8 + f);
        const int rc = nfit(l, u.fam, u.aug ? xa : x, n, u.aug ? p + 1 : p, tv, tr, mt,
                            te, me, sd, pr + static_cast<long>(j) * n, tw, liwk);
        if (rc != CLAB_OK) return rc;
      }
    }
    double *pa = psi + static_cast<long>(r) * n, *pb = psi + static_cast<long>(nrep + r) * n;
    double v;
    const int rc = score(model, sc, y, d, z, s, pr, n, fr, nf, clu, nc, nc2, trim, pa,
                         pb, th + r, &v, work + n);
    if (rc != CLAB_OK) return rc;
    se[r] = std::sqrt(v);
  }
  // the median over repetitions, and the spread of the repetitions added
  // to the variance (DoubleML's aggregation): n, or the cluster count
  const double sf = clu == nullptr ? n : (nc2 > 0 ? std::fmin(nc, nc2) : nc);
  double *t = work;
  for (int r = 0; r < nrep; ++r) t[r] = th[r];
  kutil::hsort(t, nrep, [](double a, double b) { return a < b; });
  est[0] = nrep % 2 ? t[nrep / 2] : 0.5 * (t[nrep / 2 - 1] + t[nrep / 2]);
  for (int r = 0; r < nrep; ++r) {
    t[r] = sf * se[r] * se[r] + (th[r] - est[0]) * (th[r] - est[0]);
  }
  kutil::hsort(t, nrep, [](double a, double b) { return a < b; });
  const double md = nrep % 2 ? t[nrep / 2] : 0.5 * (t[nrep / 2 - 1] + t[nrep / 2]);
  est[1] = std::sqrt(md / sf);
  return CLAB_OK;
}

int clab_dmlxk(int model, int sc) {
  static const int k[] = {2, 3, 3, 5, 3, 3, 2, 5, 4};
  if (model < 0 || model > 8) return 0;
  // IV-type adds g; the experimental cross-sections drop m
  return k[model] + (model <= 1 && sc == 1 ? 1 : 0) - (model == 7 && sc == 1 ? 1 : 0);
}

}  // extern "C"
