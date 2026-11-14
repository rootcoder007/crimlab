// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, repeat victimisation: incidence, prevalence and
// concentration (Trickett, Osborn, Seymour and Pease 1992), the victims by
// number of victimisations against a Poisson count, and the time between
// a target's successive victimisations (Johnson 2008).
#include <algorithm>
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

int clab_rptvic(const int *id, const double *t, int n, int nid, double pop, const double *brk, int nb,
                double *tab, int *kmax, double *st, double *gap, double *lag) {
  if (n < 1 || nid < 1 || nid > n || (t == nullptr && nb > 0) || nb < 0 || nb == 1) return CLAB_EARG;
  if (std::isfinite(pop) && pop < nid) return CLAB_EDOM;
  std::vector<int> k(static_cast<std::size_t>(nid), 0);
  for (int i = 0; i < n; ++i) {
    if (id[i] < 0 || id[i] >= nid) return CLAB_EARG;
    ++k[static_cast<std::size_t>(id[i])];
  }
  int K = 0;
  for (int v : k) K = std::max(K, v);
  *kmax = K;
  // tab: (K + 1) x 3, rows k = 0..K: victims with exactly k, their
  // incidents, and the victims a Poisson count with the same mean over pop
  // targets would give (NaN without pop; row 0 needs pop)
  const bool hp = std::isfinite(pop);
  const double lam = hp ? n / pop : kutil::nan();
  const std::size_t R = static_cast<std::size_t>(K) + 1;
  for (std::size_t r = 0; r < 3 * R; ++r) tab[r] = 0.0;
  for (int v : k) tab[v] += 1;
  tab[0] = hp ? pop - nid : kutil::nan();
  double lp = -lam;  // log of the Poisson probability of k
  for (std::size_t r = 0; r < R; ++r) {
    tab[R + r] = static_cast<double>(r) * tab[r];
    if (r > 0) lp += std::log(lam) - std::log(static_cast<double>(r));
    tab[2 * R + r] = hp ? pop * std::exp(lp) : kutil::nan();
  }
  // st: incidents, victims, repeats, repeat share, concentration,
  // incidence, prevalence, victims with 2 or more, the Poisson's expected
  // number of them, and their ratio
  const double N = n, V = nid;
  double v2 = 0;
  for (int v : k) v2 += v >= 2;
  st[0] = N;
  st[1] = V;
  st[2] = N - V;
  st[3] = (N - V) / N;
  st[4] = N / V;
  st[5] = hp ? N / pop : kutil::nan();
  st[6] = hp ? V / pop : kutil::nan();
  st[7] = v2;
  st[8] = hp ? pop * -std::expm1(-lam) - pop * lam * std::exp(-lam) : kutil::nan();
  st[9] = hp ? v2 / st[8] : kutil::nan();
  // the gaps between a target's successive incidents
  gap[0] = 0;
  gap[1] = gap[2] = kutil::nan();
  if (t == nullptr) return CLAB_OK;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(t[i])) return CLAB_EDOM;
  }
  std::vector<int> o(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) o[static_cast<std::size_t>(i)] = i;
  std::sort(o.begin(), o.end(), [&](int a, int b) {
    return id[a] < id[b] || (id[a] == id[b] && (t[a] < t[b] || (t[a] == t[b] && a < b)));
  });
  std::vector<double> g;
  for (int q = 1; q < n; ++q) {
    const int a = o[static_cast<std::size_t>(q) - 1], b = o[static_cast<std::size_t>(q)];
    if (id[a] == id[b]) g.push_back(t[b] - t[a]);
  }
  const std::size_t m = g.size();
  gap[0] = static_cast<double>(m);
  if (m > 0) {
    std::sort(g.begin(), g.end());
    kutil::ksum s;
    for (double x : g) s.add(x);
    gap[1] = s.val() / static_cast<double>(m);
    gap[2] = m % 2 ? g[m / 2] : 0.5 * (g[m / 2 - 1] + g[m / 2]);
  }
  // lag: (nb - 1) x 3: gaps in [brk[j], brk[j + 1]), their share of all
  // gaps and the share up to brk[j + 1]
  for (int j = 0; j + 1 < nb; ++j) {
    if (!(brk[j] < brk[j + 1])) return CLAB_EDOM;
  }
  double cum = 0;
  for (int j = 0; j + 1 < nb; ++j) {
    const double c = static_cast<double>(std::lower_bound(g.begin(), g.end(), brk[j + 1]) -
                                         std::lower_bound(g.begin(), g.end(), brk[j]));
    cum += c;
    lag[j] = c;
    lag[nb - 1 + j] = m > 0 ? c / static_cast<double>(m) : kutil::nan();
    lag[2 * (nb - 1) + j] = m > 0 ? cum / static_cast<double>(m) : kutil::nan();
  }
  return CLAB_OK;
}
