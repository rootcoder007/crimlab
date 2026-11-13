// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Causal effects of a 0/1 treatment from (y, d, x) by design-based
// estimators: the logistic propensity (the learners of learn.h),
// Hajek-normalised inverse probability weighting, g-computation by a
// least-squares outcome model in each arm, subclassification on the
// propensity, greedy nearest-neighbour matching on the propensity with a
// caliper, and entropy balancing for the effect on the treated
// (Hainmueller 2012); standardised mean differences; and the E-value of
// VanderWeele and Ding (2017). Units enter through a list of indices, so a
// bootstrap sample is a list with repeats. Scratch is held in vectors.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <vector>

#include "kutil.h"
#include "learn.h"
#include "pev.h"
#include "pevc.h"

namespace pev {

namespace {

inline std::size_t z(int i) { return static_cast<std::size_t>(i); }
const double qnan = std::numeric_limits<double>::quiet_NaN();

}  // namespace

bool pscore(const cdat &c, const std::vector<int> &idx, std::vector<double> &e) {
  std::vector<double> cf(z(c.p) + 1), se(z(c.p) + 1), work(z(lrn::gwork(c.p))), eta(idx.size());
  double dev = 0.0;
  lrn::Dat d{c.x, c.n, c.p, c.d, idx.data(), static_cast<int>(idx.size())};
  if (lrn::glm(d, 1, 1e-10, 50, cf.data(), se.data(), &dev, work.data()) != 0) return false;
  lrn::lpred(c.x, c.n, c.p, idx.data(), static_cast<int>(idx.size()), cf[0], cf.data() + 1,
             eta.data());
  e.assign(z(c.n), qnan);
  for (std::size_t q = 0; q < idx.size(); ++q) e[z(idx[q])] = 1.0 / (1.0 + std::exp(-eta[q]));
  return true;
}

double ipw(const cdat &c, const std::vector<int> &idx, const std::vector<double> &e, int tgt) {
  // tgt: 0 ATE, 1 ATT, 2 ATC; Hajek-normalised in each arm
  kutil::ksum s1, w1, s0, w0;
  for (int i : idx) {
    const double p = e[z(i)];
    if (c.d[i] > 0.5) {
      const double w = tgt == 0 ? 1.0 / p : tgt == 1 ? 1.0 : (1.0 - p) / p;
      s1.add(w * c.y[i]);
      w1.add(w);
    } else {
      const double w = tgt == 0 ? 1.0 / (1.0 - p) : tgt == 1 ? p / (1.0 - p) : 1.0;
      s0.add(w * c.y[i]);
      w0.add(w);
    }
  }
  return s1.val() / w1.val() - s0.val() / w0.val();
}

double gcomp(const cdat &c, const std::vector<int> &idx, int tgt) {
  // the least-squares outcome model of each arm, averaged over the target
  std::vector<int> a[2];
  for (int i : idx) a[c.d[i] > 0.5 ? 1 : 0].push_back(i);
  std::vector<double> m[2];
  std::vector<double> cf(z(c.p) + 1), se(z(c.p) + 1), work(z(lrn::gwork(c.p))),
      eta(idx.size());
  for (int g = 0; g < 2; ++g) {
    double dev = 0.0;
    lrn::Dat d{c.x, c.n, c.p, c.y, a[g].data(), static_cast<int>(a[g].size())};
    if (lrn::glm(d, 0, 1e-10, 50, cf.data(), se.data(), &dev, work.data()) != 0) return qnan;
    lrn::lpred(c.x, c.n, c.p, idx.data(), static_cast<int>(idx.size()), cf[0], cf.data() + 1,
               eta.data());
    m[g] = eta;
  }
  kutil::ksum s;
  double k = 0;
  for (std::size_t q = 0; q < idx.size(); ++q) {
    const bool t = c.d[idx[q]] > 0.5;
    if ((tgt == 1 && !t) || (tgt == 2 && t)) continue;
    s.add(m[1][q] - m[0][q]);
    k += 1;
  }
  return s.val() / k;
}

double subcl(const cdat &c, const std::vector<int> &idx, const std::vector<double> &e, int tgt,
             int ns) {
  // cut points: quantiles (type 7) of the propensity over the target units
  std::vector<double> v;
  for (int i : idx) {
    const bool t = c.d[i] > 0.5;
    if (tgt == 0 || (tgt == 1 && t) || (tgt == 2 && !t)) v.push_back(e[z(i)]);
  }
  std::sort(v.begin(), v.end());
  std::vector<double> cut(z(ns - 1));
  for (int k = 1; k < ns; ++k) {
    const double h = (static_cast<double>(v.size()) - 1) * k / ns;
    const std::size_t lo = static_cast<std::size_t>(std::floor(h));
    const std::size_t hi = std::min(lo + 1, v.size() - 1);
    cut[z(k - 1)] = v[lo] + (h - static_cast<double>(lo)) * (v[hi] - v[lo]);
  }
  std::vector<double> sy(2 * z(ns), 0.0), nn(2 * z(ns), 0.0);
  for (int i : idx) {
    const int s = static_cast<int>(std::lower_bound(cut.begin(), cut.end(), e[z(i)]) - cut.begin());
    const int g = c.d[i] > 0.5 ? 1 : 0;
    sy[z(2 * s + g)] += c.y[i];
    nn[z(2 * s + g)] += 1;
  }
  kutil::ksum est;
  double tot = 0;
  for (int s = 0; s < ns; ++s) {
    const double n0 = nn[z(2 * s)], n1 = nn[z(2 * s + 1)];
    const double w = tgt == 0 ? n0 + n1 : tgt == 1 ? n1 : n0;
    if (w == 0) continue;
    if (n0 == 0 || n1 == 0) return qnan;
    est.add(w * (sy[z(2 * s + 1)] / n1 - sy[z(2 * s)] / n0));
    tot += w;
  }
  return est.val() / tot;
}

int psm(const cdat &c, const std::vector<double> &e, int tgt, int k, double cal,
        std::vector<int> &mt, std::vector<int> &mc, std::vector<double> &dif) {
  // greedy, without replacement, the focal units (treated for the ATT,
  // controls for the ATC) in decreasing propensity, each taking its k
  // nearest free units of the other arm within cal (absolute propensity;
  // 0 for none)
  const double foc = tgt == 2 ? 0.0 : 1.0;
  std::vector<int> f, o;
  for (int i = 0; i < c.n; ++i) ((c.d[i] > 0.5) == (foc > 0.5) ? f : o).push_back(i);
  std::stable_sort(f.begin(), f.end(), [&](int a, int b) { return e[z(a)] > e[z(b)]; });
  std::vector<char> used(z(c.n), 0);
  mt.clear();
  mc.clear();
  dif.clear();
  for (int i : f) {
    std::vector<int> got;
    for (int r = 0; r < k; ++r) {
      int best = -1;
      for (int j : o) {
        if (used[z(j)]) continue;
        const double dd = std::fabs(e[z(i)] - e[z(j)]);
        if (cal > 0 && dd > cal) continue;
        if (best < 0 || dd < std::fabs(e[z(i)] - e[z(best)])) best = j;
      }
      if (best < 0) break;
      used[z(best)] = 1;
      got.push_back(best);
    }
    if (got.empty()) continue;
    double m = 0.0;
    for (int j : got) {
      mt.push_back(i);
      mc.push_back(j);
      m += c.y[j] / static_cast<double>(got.size());
    }
    // the focal unit's outcome less its matches' mean, treated minus control
    dif.push_back(foc > 0.5 ? c.y[i] - m : m - c.y[i]);
  }
  return static_cast<int>(dif.size());
}

bool ebal(const cdat &c, const std::vector<int> &idx, std::vector<double> &w) {
  // control weights exp(l'(x - xbar1)) / sum, matching the treated means:
  // Newton on the convex dual log sum exp(l'(x - xbar1))
  const int p = c.p;
  std::vector<double> m1(z(p), 0.0);
  std::vector<int> ctl;
  double n1 = 0;
  for (int i : idx) {
    if (c.d[i] > 0.5) {
      for (int j = 0; j < p; ++j) m1[z(j)] += c.x[z(j) * z(c.n) + z(i)];
      n1 += 1;
    } else {
      ctl.push_back(i);
    }
  }
  for (double &v : m1) v /= n1;
  std::vector<double> l(z(p), 0.0), g(z(p)), h(z(p) * z(p)), hi(z(p) * z(p)), ww(ctl.size());
  auto eval = [&](const std::vector<double> &lam, double *obj) {
    double mx = -HUGE_VAL;
    for (std::size_t q = 0; q < ctl.size(); ++q) {
      double s = 0.0;
      for (int j = 0; j < p; ++j) s += lam[z(j)] * (c.x[z(j) * z(c.n) + z(ctl[q])] - m1[z(j)]);
      ww[q] = s;
      mx = std::max(mx, s);
    }
    double t = 0.0;
    for (double &v : ww) {
      v = std::exp(v - mx);
      t += v;
    }
    for (double &v : ww) v /= t;
    *obj = mx + std::log(t);
  };
  double f = 0.0;
  eval(l, &f);
  for (int it = 0; it < 200; ++it) {
    // gradient: the weighted means less the treated means; Hessian: their
    // weighted covariance
    std::fill(g.begin(), g.end(), 0.0);
    std::fill(h.begin(), h.end(), 0.0);
    for (std::size_t q = 0; q < ctl.size(); ++q) {
      for (int j = 0; j < p; ++j) g[z(j)] += ww[q] * (c.x[z(j) * z(c.n) + z(ctl[q])] - m1[z(j)]);
    }
    double gm = 0.0;
    for (double v : g) gm = std::max(gm, std::fabs(v));
    if (gm < 1e-10) break;
    for (std::size_t q = 0; q < ctl.size(); ++q) {
      for (int j = 0; j < p; ++j) {
        const double a = c.x[z(j) * z(c.n) + z(ctl[q])] - m1[z(j)] - g[z(j)];
        for (int r = 0; r < p; ++r) {
          h[z(r) * z(p) + z(j)] += ww[q] * a * (c.x[z(r) * z(c.n) + z(ctl[q])] - m1[z(r)] - g[z(r)]);
        }
      }
    }
    if (!chol(h.data(), p)) return false;
    cinv(h.data(), p, hi.data());
    std::vector<double> st(z(p), 0.0), l2(z(p));
    for (int j = 0; j < p; ++j) {
      for (int r = 0; r < p; ++r) st[z(j)] -= hi[z(r) * z(p) + z(j)] * g[z(r)];
    }
    // halve the step until the dual falls
    double a = 1.0, f2 = 0.0;
    for (int b = 0; b < 60; ++b, a *= 0.5) {
      for (int j = 0; j < p; ++j) l2[z(j)] = l[z(j)] + a * st[z(j)];
      eval(l2, &f2);
      if (f2 <= f) break;
    }
    l = l2;
    f = f2;
    if (it == 199) return false;
  }
  eval(l, &f);
  w.assign(z(c.n), 0.0);
  for (std::size_t q = 0; q < ctl.size(); ++q) w[z(ctl[q])] += ww[q] * n1;
  for (int i : idx) {
    if (c.d[i] > 0.5) w[z(i)] = 1.0;
  }
  return true;
}

double ebatt(const cdat &c, const std::vector<int> &idx, const std::vector<double> &w) {
  kutil::ksum s1, s0, w0;
  double n1 = 0;
  for (int i : idx) {
    if (c.d[i] > 0.5) {
      s1.add(c.y[i]);
      n1 += 1;
    }
  }
  // a unit drawn several times carries its weight each time
  std::vector<double> cnt(z(c.n), 0.0);
  for (int i : idx) cnt[z(i)] += 1;
  for (int i = 0; i < c.n; ++i) {
    if (c.d[i] < 0.5 && cnt[z(i)] > 0) {
      s0.add(w[z(i)] * c.y[i]);
      w0.add(w[z(i)]);
    }
  }
  return s1.val() / n1 - s0.val() / w0.val();
}

void smd(const double *x, int n, const double *d, const double *w, double *out) {
  // weighted means (w NULL: none) over the unweighted pooled standard
  // deviation sqrt((s1^2 + s0^2) / 2), as cobalt's default; out: m1, m0,
  // smd, the variance ratio s1^2 / s0^2 (weighted)
  double sw[2] = {0, 0}, sm[2] = {0, 0}, u[2] = {0, 0}, uu[2] = {0, 0}, k[2] = {0, 0};
  for (int i = 0; i < n; ++i) {
    const int g = d[i] > 0.5 ? 1 : 0;
    const double wi = w ? w[i] : 1.0;
    sw[g] += wi;
    sm[g] += wi * x[i];
    u[g] += x[i];
    k[g] += 1;
  }
  const double m[2] = {sm[0] / sw[0], sm[1] / sw[1]}, um[2] = {u[0] / k[0], u[1] / k[1]};
  double vw[2] = {0, 0}, v1[2] = {0, 0};
  for (int i = 0; i < n; ++i) {
    const int g = d[i] > 0.5 ? 1 : 0;
    const double wi = w ? w[i] : 1.0;
    uu[g] += (x[i] - um[g]) * (x[i] - um[g]);
    vw[g] += wi * (x[i] - m[g]) * (x[i] - m[g]);
    v1[g] += wi * wi;
  }
  const double s1 = uu[1] / (k[1] - 1), s0 = uu[0] / (k[0] - 1);
  // weighted variances with the unbiased correction of cobalt's col_w_sd
  const double q1 = vw[1] * sw[1] / (sw[1] * sw[1] - v1[1]),
               q0 = vw[0] * sw[0] / (sw[0] * sw[0] - v1[0]);
  out[0] = m[1];
  out[1] = m[0];
  out[2] = (m[1] - m[0]) / std::sqrt((s1 + s0) / 2);
  out[3] = q1 / q0;
}

void evalue(double est, double lo, double hi, int type, int rare, double sd, double *out) {
  // to a risk ratio: type 0 RR, 1 OR, 2 HR, 3 a standardised difference
  // (exp(0.91 d); sd the standard error of d for the limits)
  auto rr = [&](double v) {
    if (type == 1 && !rare) return std::sqrt(v);
    if (type == 2 && !rare) {
      return (1.0 - std::pow(0.5, std::sqrt(v))) / (1.0 - std::pow(0.5, std::sqrt(1.0 / v)));
    }
    if (type == 3) return std::exp(0.91 * v);
    return v;
  };
  auto ev = [](double r) {
    if (r < 1) r = 1.0 / r;
    return r + std::sqrt(r * (r - 1.0));
  };
  const double r = rr(est);
  double l = std::isnan(lo) ? qnan : rr(lo), u = std::isnan(hi) ? qnan : rr(hi);
  if (type == 3 && !std::isnan(sd)) {
    l = std::exp(0.91 * est - 1.78 * sd);
    u = std::exp(0.91 * est + 1.78 * sd);
  }
  out[0] = r;
  out[1] = ev(r);
  // the limit nearer the null; 1 when the interval holds it
  double lim = qnan;
  if (r >= 1 && !std::isnan(l)) lim = l <= 1 ? 1.0 : ev(l);
  if (r < 1 && !std::isnan(u)) lim = u >= 1 ? 1.0 : ev(u);
  out[2] = lim;
}

}  // namespace pev
