// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, the space-time permutation scan statistic: cylinders of
// circles about the incidents and windows of time, the Poisson
// log-likelihood ratio of each against the expectation from the margins,
// and a permutation test of the most likely clusters.
#include <algorithm>
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

namespace {

// the log-likelihood ratio of c incidents in a cylinder expecting e, out
// of N (0 unless c > e; 0 log 0 = 0)
double llr(double c, double e, double N) {
  if (!(c > e)) return 0.0;
  const double out = N - c;
  return c * std::log(c / e) + (out > 0 ? out * std::log(out / (N - e)) : 0.0);
}

struct cyl {
  double llr;
  int ctr, rad, a, b, w;  // centre, radius index, window ranks [a, b), length
  double c, e, ns, nt;
};

struct scan {
  int n, nr, nw, nu, minc;
  const double *x, *y, *r, *tw;
  std::vector<double> u;               // the distinct times, sorted
  std::vector<int> ctr;                // the distinct locations' first incident
  std::vector<std::vector<int>> memb;  // per centre: incidents by distance
  std::vector<std::vector<int>> cut;   // per centre: members within r[k]
  std::vector<int> wend;               // per window length and start: end rank
  std::vector<double> cg;              // incidents before each rank

  // the best cylinder of each centre and radius for the time ranks tr, or
  // only the largest ratio (all = false)
  double run(const int *tr, bool all, std::vector<cyl> *best) {
    std::vector<double> cnt(nu), pre(nu + 1);
    double top = 0;
    for (std::size_t c = 0; c < ctr.size(); ++c) {
      std::fill(cnt.begin(), cnt.end(), 0.0);
      int got = 0;
      for (int k = 0; k < nr; ++k) {
        const int m = cut[c][k];
        if (m == got && k > 0) {
          if (all) {
            cyl z = (*best)[(c * nr) + k - 1];
            z.rad = k;
            (*best)[c * nr + k] = z;
          }
          continue;
        }
        for (; got < m; ++got) cnt[tr[memb[c][got]]] += 1;
        pre[0] = 0;
        for (int q = 0; q < nu; ++q) pre[q + 1] = pre[q] + cnt[q];
        cyl z = {0, static_cast<int>(c), k, 0, 0, 0, 0, 0, static_cast<double>(m), 0};
        for (int w = 0; w < nw; ++w) {
          const int *we = wend.data() + static_cast<std::size_t>(w) * nu;
          for (int a = 0; a < nu; ++a) {
            const int b = we[a];
            const double ci = pre[b] - pre[a];
            if (ci < minc) continue;
            const double nt = cg[b] - cg[a];
            const double e = m * nt / n;
            const double l = llr(ci, e, n);
            if (l > z.llr) {
              z.llr = l;
              z.a = a;
              z.b = b;
              z.w = w;
              z.c = ci;
              z.e = e;
              z.nt = nt;
            }
          }
        }
        if (z.llr > top) top = z.llr;
        if (all) (*best)[c * nr + k] = z;
      }
    }
    return top;
  }
};

double d2(const double *x, const double *y, int i, int j) {
  return (x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]);
}

}  // namespace

int clab_kscan(const double *x, const double *y, const double *t, int n, const double *r, int nr,
               const double *tw, int nw, int minc, int nsim, unsigned int seed, int top,
               double *out, int *ncl) {
  if (n < 3 || nr < 1 || nw < 1 || minc < 1 || nsim < 0 || top < 1) return CLAB_EARG;
  for (int k = 0; k < nr; ++k) {
    if (!kutil::pos(r[k]) || (k > 0 && !(r[k] > r[k - 1]))) return CLAB_EDOM;
  }
  for (int w = 0; w < nw; ++w) {
    if (!kutil::pos(tw[w])) return CLAB_EDOM;
  }
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(x[i]) || !std::isfinite(y[i]) || !std::isfinite(t[i])) return CLAB_EDOM;
  }
  scan s;
  s.n = n;
  s.nr = nr;
  s.nw = nw;
  s.minc = minc;
  s.x = x;
  s.y = y;
  s.r = r;
  s.tw = tw;
  s.u.assign(t, t + n);
  std::sort(s.u.begin(), s.u.end());
  s.u.erase(std::unique(s.u.begin(), s.u.end()), s.u.end());
  s.nu = static_cast<int>(s.u.size());
  std::vector<int> tr(n);
  for (int i = 0; i < n; ++i) {
    tr[i] = static_cast<int>(std::lower_bound(s.u.begin(), s.u.end(), t[i]) - s.u.begin());
  }
  s.cg.assign(s.nu + 1, 0.0);
  for (int i = 0; i < n; ++i) s.cg[tr[i] + 1] += 1;
  for (int q = 0; q < s.nu; ++q) s.cg[q + 1] += s.cg[q];
  // windows [u_a, u_a + tw) as ranks [a, b)
  s.wend.resize(static_cast<std::size_t>(nw) * s.nu);
  for (int w = 0; w < nw; ++w) {
    for (int a = 0; a < s.nu; ++a) {
      s.wend[static_cast<std::size_t>(w) * s.nu + a] =
          static_cast<int>(std::lower_bound(s.u.begin(), s.u.end(), s.u[a] + tw[w]) - s.u.begin());
    }
  }
  // centres: each distinct location once, its incidents within the largest
  // radius ordered by distance
  std::vector<int> o(n);
  for (int i = 0; i < n; ++i) o[i] = i;
  std::sort(o.begin(), o.end(), [&](int a, int b) {
    return x[a] == x[b] ? (y[a] == y[b] ? a < b : y[a] < y[b]) : x[a] < x[b];
  });
  for (int q = 0; q < n; ++q) {
    if (q == 0 || x[o[q]] != x[o[q - 1]] || y[o[q]] != y[o[q - 1]]) s.ctr.push_back(o[q]);
  }
  std::sort(s.ctr.begin(), s.ctr.end());
  const double rmax2 = r[nr - 1] * r[nr - 1];
  s.memb.resize(s.ctr.size());
  s.cut.resize(s.ctr.size());
  for (std::size_t c = 0; c < s.ctr.size(); ++c) {
    const int i = s.ctr[c];
    std::vector<int> &m = s.memb[c];
    for (int j = 0; j < n; ++j) {
      if (d2(x, y, i, j) <= rmax2) m.push_back(j);
    }
    std::sort(m.begin(), m.end(), [&](int a, int b) {
      const double da = d2(x, y, i, a), db = d2(x, y, i, b);
      return da < db || (da == db && a < b);
    });
    s.cut[c].resize(nr);
    std::size_t q = 0;
    for (int k = 0; k < nr; ++k) {
      while (q < m.size() && d2(x, y, i, m[q]) <= r[k] * r[k]) ++q;
      s.cut[c][k] = static_cast<int>(q);
    }
  }
  std::vector<cyl> best(s.ctr.size() * nr);
  s.run(tr.data(), true, &best);
  // the most likely cluster and the secondary ones whose circles overlap
  // none chosen before
  std::sort(best.begin(), best.end(), [](const cyl &a, const cyl &b) {
    return a.llr > b.llr ||
           (a.llr == b.llr && (a.ctr < b.ctr || (a.ctr == b.ctr && a.rad < b.rad)));
  });
  std::vector<cyl> pick;
  for (const cyl &z : best) {
    if (static_cast<int>(pick.size()) >= top || !(z.llr > 0)) break;
    bool far = true;
    for (const cyl &q : pick) {
      const double sr = r[z.rad] + r[q.rad];
      far = far && d2(x, y, s.ctr[z.ctr], s.ctr[q.ctr]) > sr * sr;
    }
    if (far) pick.push_back(z);
  }
  // the null distribution of the largest ratio: times shuffled
  std::vector<double> mx(nsim);
  kutil::rng g(seed);
  for (int k = 0; k < nsim; ++k) {
    for (int i = n - 1; i > 0; --i) std::swap(tr[i], tr[g.below(i + 1)]);
    mx[k] = s.run(tr.data(), false, nullptr);
  }
  *ncl = static_cast<int>(pick.size());
  for (std::size_t q = 0; q < pick.size(); ++q) {
    const cyl &z = pick[q];
    double ge = 0;
    for (int k = 0; k < nsim; ++k) ge += mx[k] >= z.llr;
    double *v = out + 10 * q;
    v[0] = s.ctr[z.ctr];
    v[1] = r[z.rad];
    v[2] = s.u[z.a];
    v[3] = s.u[z.a] + tw[z.w];
    v[4] = z.c;
    v[5] = z.e;
    v[6] = z.ns;
    v[7] = z.nt;
    v[8] = z.llr;
    v[9] = nsim > 0 ? (1 + ge) / (nsim + 1.0) : kutil::nan();
  }
  return CLAB_OK;
}
