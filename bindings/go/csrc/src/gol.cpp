// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, the generalised ordered logit (crimlab.h, gologt()): Peterson
// and Harrell (1990), Williams (2006). Category y in 0..M-1; for each
// threshold j < M - 1, P(y > j) = F(a_j + x b_j) with F the logistic
// distribution function, the slopes of a covariate either one per
// threshold or one for all (held parallel). Fitted by maximum likelihood
// (quasi-Newton steps with the analytic gradient), with the covariance
// from the Hessian and, for each covariate not held parallel, the Wald
// test that its slopes are equal across the thresholds.
#include <cmath>
#include <cstdint>
#include <vector>

#include "crimlab.h"
#include "kutil.h"
#include "qn.h"

namespace {

using i64 = std::int64_t;

struct gm {
  const int *y;
  const double *x;
  int n, M, p;
  const int *off;  // the first parameter of each covariate's slopes
  const int *par;  // 1: one slope for all the thresholds
};

double expit(double e) { return 1.0 / (1.0 + std::exp(-e)); }

// The negative log-likelihood at u and its gradient into g (when not null);
// HUGE_VAL where a category's probability is not positive.
double nll(const double *u, double *g, void *ctx) {
  const gm &m = *static_cast<const gm *>(ctx);
  const int J = m.M - 1;
  const int np = m.off[m.p];
  if (g != nullptr) {
    for (int k = 0; k < np; ++k) g[k] = 0.0;
  }
  kutil::ksum s;
  for (int i = 0; i < m.n; ++i) {
    const int c = m.y[i];
    // eta at the thresholds c - 1 and c, those that bound category c
    double eta[2] = {0.0, 0.0};
    const int jj[2] = {c - 1, c};
    for (int h = 0; h < 2; ++h) {
      const int j = jj[h];
      if (j < 0 || j >= J) continue;
      double e = u[j];
      for (int k = 0; k < m.p; ++k) e += u[m.off[k] + (m.par[k] ? 0 : j)] * m.x[static_cast<i64>(k) * m.n + i];
      eta[h] = e;
    }
    const double Fa = c == 0 ? 1.0 : expit(eta[0]), Fb = c == J ? 0.0 : expit(eta[1]);
    const double pr = Fa - Fb;
    if (!(pr > 0.0)) return HUGE_VAL;
    s.add(std::log(pr));
    if (g == nullptr) continue;
    const double da = c == 0 ? 0.0 : Fa * (1.0 - Fa) / pr, db = c == J ? 0.0 : -Fb * (1.0 - Fb) / pr;
    for (int h = 0; h < 2; ++h) {
      const int j = jj[h];
      const double dv = h == 0 ? da : db;
      if (j < 0 || j >= J || dv == 0.0) continue;
      g[j] -= dv;
      for (int k = 0; k < m.p; ++k) g[m.off[k] + (m.par[k] ? 0 : j)] -= dv * m.x[static_cast<i64>(k) * m.n + i];
    }
  }
  const double v = -s.val();
  return std::isfinite(v) ? v : HUGE_VAL;
}

}  // namespace

int clab_gologt(const int *y, int n, int M, const double *x, int p, const int *par, int maxit, double *est,
                double *cov, double *wd, double *st) CLAB_TRY {
  if (n < 2 || M < 3 || p < 0 || maxit < 1) return CLAB_EARG;
  const int J = M - 1;
  std::vector<int> off(static_cast<std::size_t>(p) + 1), cnt(static_cast<std::size_t>(M), 0);
  off[0] = J;
  for (int k = 0; k < p; ++k) {
    if (par[k] != 0 && par[k] != 1) return CLAB_EARG;
    off[static_cast<std::size_t>(k) + 1] = off[static_cast<std::size_t>(k)] + (par[k] ? 1 : J);
  }
  for (int i = 0; i < n; ++i) {
    if (y[i] < 0 || y[i] >= M) return CLAB_EARG;
    cnt[static_cast<std::size_t>(y[i])] += 1;
  }
  for (int c = 0; c < M; ++c) {
    if (cnt[static_cast<std::size_t>(c)] == 0) return CLAB_EDOM;
  }
  for (i64 q = 0; q < static_cast<i64>(n) * p; ++q) {
    if (!std::isfinite(x[q])) return CLAB_EDOM;
  }
  const int np = off[static_cast<std::size_t>(p)];
  // start: the intercepts at the logits of the shares above each
  // threshold, every slope 0
  for (int k = 0; k < np; ++k) est[k] = 0.0;
  double above = n;
  for (int j = 0; j < J; ++j) {
    above -= cnt[static_cast<std::size_t>(j)];
    est[j] = std::log(above / (n - above));
  }
  gm m{y, x, n, M, p, off.data(), par};
  std::vector<double> w(static_cast<std::size_t>(np) * np + 5 * static_cast<std::size_t>(np));
  double fv = 0.0;
  int it = 0;
  const int c = qn::min(np, est, nll, &m, maxit, 1e-8, &fv, &it, w.data(), static_cast<double>(n) * (p + 2));
  if (kutil::stop(0.0)) return CLAB_EINTR;
  if (c < 0 || !std::isfinite(fv)) return CLAB_EDOM;
  std::vector<double> hw(3 * static_cast<std::size_t>(np));
  qn::hcov(np, est, nll, &m, cov, hw.data());
  // the Wald test that a covariate's J slopes are equal: the J - 1
  // differences from the first, D b, against D V D'
  const int q = J - 1;
  std::vector<double> a(static_cast<std::size_t>(q) * q), b(static_cast<std::size_t>(q));
  for (int k = 0; k < p; ++k) {
    double *o = wd + 3 * k;
    if (par[k]) {
      o[0] = o[1] = o[2] = kutil::nan();
      continue;
    }
    const int f = off[static_cast<std::size_t>(k)];
    auto V = [&](int r, int s) { return cov[static_cast<i64>(f + s) * np + f + r]; };
    for (int r = 0; r < q; ++r) {
      b[static_cast<std::size_t>(r)] = est[f] - est[f + r + 1];
      for (int s = 0; s < q; ++s) {
        a[static_cast<std::size_t>(s) * q + r] = V(0, 0) - V(0, s + 1) - V(r + 1, 0) + V(r + 1, s + 1);
      }
    }
    if (!std::isfinite(V(0, 0)) || qn::spdinv(a.data(), q) != 0) {
      o[0] = o[2] = kutil::nan();
      o[1] = q;
      continue;
    }
    double w2 = 0.0;
    for (int r = 0; r < q; ++r) {
      for (int s = 0; s < q; ++s) w2 += b[static_cast<std::size_t>(r)] * a[static_cast<std::size_t>(s) * q + r] *
                                         b[static_cast<std::size_t>(s)];
    }
    o[0] = w2;
    o[1] = q;
    o[2] = clab_pchisq(w2, q, 0);
  }
  st[0] = -fv;
  st[1] = it;
  st[2] = c;
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)
