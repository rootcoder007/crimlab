// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Staggered adoption on a balanced panel: the group-time effects of
// Callaway and Sant'Anna (2021), with the doubly robust, inverse
// probability weighting and outcome regression estimators of Sant'Anna
// and Zhao (2020) for covariates; Wooldridge's (2021) regression on
// cohort-period cells; de Chaisemartin and D'Haultfoeuille's (2020)
// estimator for treatments that switch on and off; and the cluster
// bootstrap of any of them. Units enter through a list of indices, so a
// bootstrap sample is a list with repeats. Scratch is held in vectors.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "crimlab.h"
#include "kutil.h"
#include "learn.h"
#include "pev.h"
#include "pevs.h"

namespace pev {

namespace {

inline std::size_t z(int i) { return static_cast<std::size_t>(i); }
const double qnan = std::numeric_limits<double>::quiet_NaN();

// The doubly robust, IPW or regression estimate of one group-time cell:
// the units idx[0..m) with dd 1 (the cohort) or 0 (the comparison),
// outcome change dy, covariates x (nu x k at the base period).
double drcell(const double *x, int nu, int k, const double *dy, const double *dd,
              const std::vector<int> &all, const std::vector<int> &ctl, int est) {
  std::vector<double> mh(z(nu), 0.0), ps(z(nu), 0.0);
  std::vector<double> work(z(lrn::gwork(k)));
  std::vector<double> cf(z(k) + 1), se(z(k) + 1), eta(all.size());
  double dev = 0.0;
  if (est != 1) {
    // the outcome regression among the comparison units
    lrn::Dat d{x, nu, k, dy, ctl.data(), static_cast<int>(ctl.size())};
    if (lrn::glm(d, 0, 1e-10, 50, cf.data(), se.data(), &dev, work.data()) != CLAB_OK) return qnan;
    lrn::lpred(x, nu, k, all.data(), static_cast<int>(all.size()), cf[0], cf.data() + 1, eta.data());
    for (std::size_t q = 0; q < all.size(); ++q) mh[z(all[q])] = eta[q];
  }
  if (est != 2) {
    lrn::Dat d{x, nu, k, dd, all.data(), static_cast<int>(all.size())};
    if (lrn::glm(d, 1, 1e-10, 50, cf.data(), se.data(), &dev, work.data()) != CLAB_OK) return qnan;
    lrn::lpred(x, nu, k, all.data(), static_cast<int>(all.size()), cf[0], cf.data() + 1, eta.data());
    for (std::size_t q = 0; q < all.size(); ++q) {
      ps[z(all[q])] = std::min(1.0 / (1.0 + std::exp(-eta[q])), 1.0 - 1e-6);
    }
  }
  kutil::ksum s1, w1, s0, w0;
  for (int i : all) {
    const double r = dy[i] - mh[z(i)];
    if (dd[i] > 0.5) {
      s1.add(r);
      w1.add(1.0);
    } else if (est != 2) {
      // comparison units whose propensity is 0.995 or more are trimmed
      const double w = ps[z(i)] < 0.995 ? ps[z(i)] / (1.0 - ps[z(i)]) : 0.0;
      s0.add(w * r);
      w0.add(w);
    }
  }
  const double a = s1.val() / w1.val();
  if (est == 2) return a;
  return w0.val() > 0 ? a - s0.val() / w0.val() : qnan;
}

}  // namespace

void cells(const int *g, int nu, int nt, std::vector<int> &cg, std::vector<int> &ct, int post) {
  std::vector<char> has(z(nt), 0);
  for (int i = 0; i < nu; ++i) {
    if (g[i] >= 1 && g[i] < nt) has[z(g[i])] = 1;
  }
  for (int a = 1; a < nt; ++a) {
    if (!has[z(a)]) continue;
    for (int s = post ? a : 1; s < nt; ++s) {
      cg.push_back(a);
      ct.push_back(s);
    }
  }
}

void csa(const stg &s, const std::vector<int> &idx, double *att) {
  // multiplicities and the indicator vectors
  const int nc = static_cast<int>(s.cg.size());
  std::vector<double> dy(z(s.nu)), dd(z(s.nu)), xb(s.k > 0 ? z(s.nu) * z(s.k) : 0);
  for (int c = 0; c < nc; ++c) {
    const int gc = s.cg[z(c)], tc = s.ct[z(c)];
    const int b = tc >= gc ? gc - 1 : tc - 1;
    const int hi = std::max(tc, b);
    std::vector<int> all, ctl;
    for (int i : idx) {
      const int gi = s.g[i];
      const bool coh = gi == gc;
      const bool cmp = s.ctrl == 1 ? gi >= s.nt : (gi > hi && gi != gc);
      if (!coh && !cmp) continue;
      all.push_back(i);
      if (cmp) ctl.push_back(i);
    }
    int n1 = 0;
    for (int i : all) n1 += s.g[i] == gc;
    if (n1 == 0 || ctl.empty()) {
      att[c] = qnan;
      continue;
    }
    for (int i : all) {
      dy[z(i)] = s.yw[z(tc) * z(s.nu) + z(i)] - s.yw[z(b) * z(s.nu) + z(i)];
      dd[z(i)] = s.g[i] == gc ? 1.0 : 0.0;
    }
    if (s.k == 0) {
      kutil::ksum a1, a0;
      for (int i : all) (dd[z(i)] > 0.5 ? a1 : a0).add(dy[z(i)]);
      att[c] = a1.val() / n1 - a0.val() / static_cast<double>(ctl.size());
      continue;
    }
    for (int j = 0; j < s.k; ++j) {
      for (int i = 0; i < s.nu; ++i) {
        xb[z(j) * z(s.nu) + z(i)] = s.xw[(z(j) * z(s.nt) + z(b)) * z(s.nu) + z(i)];
      }
    }
    att[c] = drcell(xb.data(), s.nu, s.k, dy.data(), dd.data(), all, ctl, s.est);
  }
}

void wool(const stg &s, const std::vector<int> &idx, double *att) {
  // the rows of the sampled units (a unit drawn twice is two units), the
  // periods with a comparison, and one indicator per post cell
  const int nc = static_cast<int>(s.cg.size());
  bool never = false;
  int glast = 0;
  for (int i : idx) {
    never = never || s.g[i] >= s.nt;
    if (s.g[i] < s.nt) glast = std::max(glast, s.g[i]);
  }
  const int tmax = never ? s.nt : glast;
  std::vector<double> y, x;
  std::vector<int> u, t;
  int nu = 0;
  for (int i : idx) {
    for (int a = 0; a < tmax; ++a) {
      y.push_back(s.yw[z(a) * z(s.nu) + z(i)]);
      u.push_back(nu);
      t.push_back(a);
    }
    ++nu;
  }
  const int n = static_cast<int>(y.size());
  // the cells that occur in the sample and have a period kept
  std::vector<int> col(z(nc), -1);
  int p = 0;
  for (int c = 0; c < nc; ++c) {
    bool in = s.ct[z(c)] < tmax;
    bool seen = false;
    for (int i : idx) seen = seen || s.g[i] == s.cg[z(c)];
    if (in && seen) col[z(c)] = p++;
    att[c] = qnan;
  }
  if (p == 0) return;
  x.assign(z(n) * z(p), 0.0);
  int r = 0;
  for (int i : idx) {
    for (int a = 0; a < tmax; ++a, ++r) {
      for (int c = 0; c < nc; ++c) {
        if (col[z(c)] >= 0 && s.g[i] == s.cg[z(c)] && a == s.ct[z(c)]) {
          x[z(col[z(c)]) * z(n) + z(r)] = 1.0;
        }
      }
    }
  }
  std::vector<double> beta(z(p)), vc(z(p) * z(p)), st(5),
      work(z(n) * z(p + 2) + z(p) * z(p) + 2 * z(nu + tmax) + z(n + 2 * p) * z(p) + z(p));
  std::vector<int> iwork(z(n) + z(std::max(nu, tmax)));
  // the cells are distinct and every period kept has comparison rows, so
  // the fit does not fail; NaN stays if it did
  const int rc = clab_didfe(y.data(), x.data(), n, p, u.data(), nu, t.data(), tmax, 1, nullptr,
                            u.data(), nu, nullptr, 0, 1e-10, 10000, beta.data(), vc.data(),
                            st.data(), work.data(), iwork.data());
  for (int c = 0; c < nc; ++c) {
    if (rc == CLAB_OK && col[z(c)] >= 0) att[c] = beta[z(col[z(c)])];
  }
}

void dcdh(const stg &s, const std::vector<int> &idx, double *out) {
  // out: the effect, the placebo, the switchers of each
  kutil::ksum num, den, pnum, pden;
  for (int a = 1; a < s.nt; ++a) {
    // joiners (0 to 1) against stayers at 0; leavers (1 to 0) against
    // stayers at 1; the placebo one period earlier for units stable then
    for (int pl = 0; pl < 2; ++pl) {
      if (pl == 1 && a < 2) continue;
      for (int dir = 0; dir < 2; ++dir) {
        kutil::ksum ms, mc;
        double ns = 0, nc = 0;
        for (int i : idx) {
          const int d0 = s.dw[z(a - 1) * z(s.nu) + z(i)], d1 = s.dw[z(a) * z(s.nu) + z(i)];
          if (d0 != dir) continue;
          if (pl == 1 && s.dw[z(a - 2) * z(s.nu) + z(i)] != dir) continue;
          const int lo = pl == 1 ? a - 2 : a - 1, hi = pl == 1 ? a - 1 : a;
          const double dy = s.yw[z(hi) * z(s.nu) + z(i)] - s.yw[z(lo) * z(s.nu) + z(i)];
          if (d1 != dir) {
            ms.add(dy);
            ns += 1;
          } else {
            mc.add(dy);
            nc += 1;
          }
        }
        if (ns == 0 || nc == 0) continue;
        double did = ms.val() / ns - mc.val() / nc;
        if (dir == 1) did = -did;
        (pl ? pnum : num).add(ns * did);
        (pl ? pden : den).add(ns);
      }
    }
  }
  out[0] = den.val() > 0 ? num.val() / den.val() : qnan;
  out[1] = pden.val() > 0 ? pnum.val() / pden.val() : qnan;
  out[2] = den.val();
  out[3] = pden.val();
}

void aggr(const stg &s, const std::vector<int> &idx, const double *att, std::vector<int> &ev,
          double *agg) {
  // cohort sizes in the sample
  std::vector<double> ng(z(s.nt) + 1, 0.0);
  for (int i : idx) {
    if (s.g[i] < s.nt) ng[z(s.g[i])] += 1.0;
  }
  const int nc = static_cast<int>(s.cg.size()), ne = static_cast<int>(ev.size());
  kutil::ksum on, od;
  std::vector<double> en(z(ne), 0.0), ed(z(ne), 0.0);
  for (int c = 0; c < nc; ++c) {
    if (std::isnan(att[c])) continue;
    const double w = ng[z(s.cg[z(c)])];
    const int e = s.ct[z(c)] - s.cg[z(c)];
    if (e >= 0) {
      on.add(w * att[c]);
      od.add(w);
    }
    const int q = static_cast<int>(std::lower_bound(ev.begin(), ev.end(), e) - ev.begin());
    en[z(q)] += w * att[c];
    ed[z(q)] += w;
  }
  agg[0] = od.val() > 0 ? on.val() / od.val() : qnan;
  for (int q = 0; q < ne; ++q) agg[q + 1] = ed[z(q)] > 0 ? en[z(q)] / ed[z(q)] : qnan;
}

void stat(const stg &s, const std::vector<int> &idx, std::vector<int> &ev, double *v) {
  const int nc = static_cast<int>(s.cg.size());
  if (s.meth == 2) {
    dcdh(s, idx, v);
    return;
  }
  if (s.meth == 0) {
    csa(s, idx, v);
  } else {
    wool(s, idx, v);
  }
  aggr(s, idx, v, ev, v + nc);
}

void boot(const stg &s, const int *cl, int ncl, int reps, std::uint64_t seed,
          std::vector<int> &ev, int ns, double *se) {
  // the units of each cluster
  std::vector<std::vector<int>> mem(z(ncl));
  for (int i = 0; i < s.nu; ++i) mem[z(cl[i])].push_back(i);
  kutil::rng g(seed);
  std::vector<kutil::ksum> s1(z(ns)), s2(z(ns));
  std::vector<double> nn(z(ns), 0.0), v(z(ns));
  std::vector<int> idx;
  for (int r = 0; r < reps; ++r) {
    if (kutil::stop(8.0 * s.nu)) return;
    idx.clear();
    for (int c = 0; c < ncl; ++c) {
      const auto &m = mem[z(g.below(ncl))];
      idx.insert(idx.end(), m.begin(), m.end());
    }
    stat(s, idx, ev, v.data());
    for (int q = 0; q < ns; ++q) {
      if (std::isnan(v[z(q)])) continue;
      s1[z(q)].add(v[z(q)]);
      s2[z(q)].add(v[z(q)] * v[z(q)]);
      nn[z(q)] += 1.0;
    }
  }
  for (int q = 0; q < ns; ++q) {
    const double k = nn[z(q)];
    if (k < 2) {
      se[q] = qnan;
      continue;
    }
    const double m = s1[z(q)].val() / k;
    const double var = (s2[z(q)].val() - k * m * m) / (k - 1);
    se[q] = std::sqrt(std::max(var, 0.0));
  }
}

}  // namespace pev
