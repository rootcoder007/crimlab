// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, more spatial methods: the nearest-neighbour index under
// Manhattan distance, head-bang smoothing, the correlated random walk of a
// series of events, and risk-adjusted nearest-neighbour clustering.
#include <cmath>
#include <limits>

#include "cgeo.h"
#include "crimlab.h"
#include "kutil.h"

namespace {

constexpr double pi = cgeo::pi;

bool finite2(const double *x, const double *y, int n) {
  for (int i = 0; i < n; ++i)
    if (!std::isfinite(x[i]) || !std::isfinite(y[i])) return false;
  return true;
}

bool inwin(const double *x, const double *y, int n, const double *w) {
  if (!(w[0] < w[1] && w[2] < w[3])) return false;
  for (int i = 0; i < n; ++i)
    if (!(x[i] >= w[0] && x[i] <= w[1] && y[i] >= w[2] && y[i] <= w[3])) return false;
  return true;
}

// the mean of each order's k-th Manhattan neighbour distance into md
void meanl1(const double *x, const double *y, int n, int kmax, int *nb,
            double *d, int *iwork, double *md) {
  clab_swknb(x, y, n, kmax, 1, nb, d, iwork);
  for (int k = 0; k < kmax; ++k) {
    kutil::ksum s;
    for (int i = 0; i < n; ++i) s.add(d[static_cast<long>(i) * kmax + k]);
    md[k] = s.val() / n;
  }
}

// the median of v[0..m) (sorted in place); the mean of the middle two for
// even m
double median(double *v, int m) {
  kutil::hsort(v, m, [](double a, double b) { return a < b; });
  return m % 2 ? v[m / 2] : v[m / 2 - 1] + 0.5 * (v[m / 2] - v[m / 2 - 1]);
}

// the weighted median of v[0..m) with weights u: the smallest value whose
// cumulative weight reaches half the total (o: m ints of scratch)
double wmedian(const double *v, const double *u, int m, int *o) {
  for (int q = 0; q < m; ++q) o[q] = q;
  kutil::hsort(o, m, [v](int a, int b) { return v[a] < v[b] || (v[a] == v[b] && a < b); });
  double tot = 0.0;
  for (int q = 0; q < m; ++q) tot += u[q];
  int q = 0;
  double run = u[o[0]];
  while (q + 1 < m && run < 0.5 * tot) run += u[o[++q]];
  return v[o[q]];
}

double med3(double a, double b, double c) {
  return std::fmax(std::fmin(a, b), std::fmin(std::fmax(a, b), c));
}

}  // namespace

extern "C" {

int clab_nnidm(const double *x, const double *y, int n, const double *win,
               int kmax, int nsim, unsigned int seed, double *dist,
               int *iwork, double *work, double *out) {
  if (n < 2 || kmax < 1 || kmax > n - 1 || nsim < 0) return CLAB_EARG;
  if (!finite2(x, y, n) || !inwin(x, y, n, win)) return CLAB_EDOM;
  const double w = win[1] - win[0], h = win[3] - win[2], area = w * h;
  if (!(area > 0.0) || !std::isfinite(area)) return CLAB_EDOM;
  int *nb = iwork + 2 * n + 2;
  double *d = work, *md = work + static_cast<long>(n) * kmax, *sx = md + kmax,
         *sy = sx + n;
  meanl1(x, y, n, kmax, nb, d, iwork, md);
  if (dist != nullptr)
    for (long q = 0; q < static_cast<long>(n) * kmax; ++q) dist[q] = d[q];
  const double lam = n / area;
  // a Manhattan ball of radius r has area 2 r^2, so under complete spatial
  // randomness (no edge) R_k^2 is Gamma(k, 2 lam): E R_k = Gamma(k + 1/2) /
  // (Gamma(k) sqrt(2 lam)) and E R_k^2 = k / (2 lam)
  double gk = std::sqrt(pi) / 2.0;  // Gamma(k + 1/2) / Gamma(k) at k = 1
  for (int k = 0; k < kmax; ++k) {
    double *o = out + 6 * k;
    o[0] = md[k];
    o[1] = gk / std::sqrt(2.0 * lam);
    o[2] = o[0] / o[1];
    o[3] = o[4] = kutil::nan();
    const double v = (k + 1) / (2.0 * lam) - o[1] * o[1];
    o[5] = std::sqrt(v > 0.0 ? v : 0.0);
    gk *= (k + 1.5) / (k + 1.0);
  }
  // order 1: no edge-corrected expectation is published for Manhattan
  // distance (NaN); the z statistic about the uncorrected expectation uses
  // Var R_1 = (4 - pi) / (8 lam)
  double *e = out + 6 * kmax;
  e[0] = e[1] = e[3] = kutil::nan();
  e[2] = (out[0] - out[1]) / std::sqrt((4.0 - pi) / (8.0 * lam * n));
  if (nsim > 0) {
    int *lo = nb + static_cast<long>(n) * kmax, *hi = lo + kmax;
    for (int k = 0; k < kmax; ++k) lo[k] = hi[k] = 0;
    double *ms = sy + n;
    cgeo::rng r(seed);
    for (int s = 0; s < nsim; ++s) {
      if (kutil::stop(16.0 * n * kmax)) return CLAB_EINTR;
      for (int i = 0; i < n; ++i) {
        sx[i] = win[0] + w * r.unif();
        sy[i] = win[2] + h * r.unif();
      }
      meanl1(sx, sy, n, kmax, nb, d, iwork, ms);
      for (int k = 0; k < kmax; ++k) {
        if (ms[k] <= out[6 * k]) ++lo[k];
        if (ms[k] >= out[6 * k]) ++hi[k];
      }
    }
    for (int k = 0; k < kmax; ++k) {
      out[6 * k + 3] = (1.0 + lo[k]) / (nsim + 1.0);
      out[6 * k + 4] = (1.0 + hi[k]) / (nsim + 1.0);
    }
  }
  return CLAB_OK;
}

int clab_hdbang(const double *x, const double *y, const double *v,
                const double *w, int n, int nn, double ang, int ntrip,
                int iter, int *iwork, double *work, double *out) {
  if (n < 3 || nn < 2 || nn > n - 1 || iter < 0 || ntrip < 0) return CLAB_EARG;
  if (!finite2(x, y, n) || !(ang >= 90.0 && ang <= 180.0)) return CLAB_EDOM;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(v[i])) return CLAB_EDOM;
    if (w != nullptr && !(w[i] > 0.0 && std::isfinite(w[i]))) return CLAB_EDOM;
  }
  const int T = nn * (nn - 1) / 2;
  int *nb = iwork + 2 * n + 2, *tj = nb + n * nn, *tk = tj + T, *ord = tk + T;
  double *d = work, *cur = d + n * nn, *lo = cur + n, *hi = lo + T, *wl = hi + T,
         *wh = wl + T, *key = wh + T;
  clab_swknb(x, y, n, nn, 0, nb, d, iwork);
  const double cmax = std::cos(ang * pi / 180.0) + 1e-12;
  for (int i = 0; i < n; ++i) out[i] = v[i];
  for (int it = 0; it < iter; ++it) {
    if (kutil::stop(static_cast<double>(n) * (T + nn))) return CLAB_EINTR;
    for (int i = 0; i < n; ++i) cur[i] = out[i];
    for (int i = 0; i < n; ++i) {
      // the triples (j, i, k) of i's neighbours whose angle at i is at
      // least ang, i.e. with i roughly between j and k
      int m = 0;
      const int *ni = nb + i * nn;
      for (int a = 0; a < nn; ++a) {
        const double ax = x[ni[a]] - x[i], ay = y[ni[a]] - y[i];
        const double la = std::hypot(ax, ay);
        if (!(la > 0.0)) continue;
        for (int b = a + 1; b < nn; ++b) {
          const double bx = x[ni[b]] - x[i], by = y[ni[b]] - y[i];
          const double lb = std::hypot(bx, by);
          if (!(lb > 0.0)) continue;
          if ((ax * bx + ay * by) / (la * lb) > cmax) continue;
          tj[m] = ni[a];
          tk[m] = ni[b];
          key[m] = la + lb;
          ++m;
        }
      }
      if (m == 0) continue;
      // keep the ntrip triples whose ends are nearest (all when 0)
      for (int q = 0; q < m; ++q) ord[q] = q;
      kutil::hsort(ord, m, [key](int p, int q) { return key[p] < key[q] || (key[p] == key[q] && p < q); });
      const int use = ntrip > 0 && ntrip < m ? ntrip : m;
      for (int q = 0; q < use; ++q) {
        const int j = tj[ord[q]], k = tk[ord[q]];
        const bool jl = cur[j] <= cur[k];
        lo[q] = jl ? cur[j] : cur[k];
        hi[q] = jl ? cur[k] : cur[j];
        if (w != nullptr) {
          wl[q] = jl ? w[j] : w[k];
          wh[q] = jl ? w[k] : w[j];
        }
      }
      double L, H;
      if (w == nullptr) {
        L = median(lo, use);
        H = median(hi, use);
      } else {
        L = wmedian(lo, wl, use, ord);
        H = wmedian(hi, wh, use, ord);
      }
      out[i] = med3(L, cur[i], H);
    }
  }
  return CLAB_OK;
}

int clab_crwfit(const double *x, const double *y, int n, double *step,
                double *turn, double *out) {
  if (n < 3) return CLAB_EARG;
  if (!finite2(x, y, n)) return CLAB_EDOM;
  for (int i = 0; i + 1 < n; ++i) {
    step[i] = std::hypot(x[i + 1] - x[i], y[i + 1] - y[i]);
    if (!(step[i] > 0.0)) return CLAB_EDOM;  // a repeat has no heading
  }
  kutil::ksum sl, sl2, sc, ss;
  for (int i = 0; i + 1 < n; ++i) {
    sl.add(step[i]);
    sl2.add(step[i] * step[i]);
  }
  for (int i = 0; i + 2 < n; ++i) {
    const double h0 = std::atan2(y[i + 1] - y[i], x[i + 1] - x[i]);
    const double h1 = std::atan2(y[i + 2] - y[i + 1], x[i + 2] - x[i + 1]);
    double t = h1 - h0;
    if (t > pi) t -= 2.0 * pi;
    if (t <= -pi) t += 2.0 * pi;
    turn[i] = t;
    sc.add(std::cos(t));
    ss.add(std::sin(t));
  }
  out[0] = sl.val() / (n - 1);
  out[1] = sl2.val() / (n - 1);
  out[2] = sc.val() / (n - 2);
  out[3] = ss.val() / (n - 2);
  out[4] = std::atan2(y[n - 1] - y[n - 2], x[n - 1] - x[n - 2]);
  return CLAB_OK;
}

int clab_crwmsd(const double *x, const double *y, int n, double ml,
                double ml2, double c, int h, double *obs, double *msd) {
  if (h < 1 || n < 0) return CLAB_EARG;
  for (int k = 1; k <= h; ++k) {
    // the observed mean squared displacement over k moves, over every start
    if (k < n) {
      kutil::ksum s;
      for (int i = 0; i + k < n; ++i) {
        const double dx = x[i + k] - x[i], dy = y[i + k] - y[i];
        s.add(dx * dx + dy * dy);
      }
      obs[k - 1] = s.val() / (n - k);
    } else {
      obs[k - 1] = kutil::nan();
    }
  }
  if (!kutil::nonneg(ml) || !kutil::nonneg(ml2) || !(c >= -1.0 && c <= 1.0)) return CLAB_EDOM;
  for (int k = 1; k <= h; ++k) {
    // Kareiva and Shigesada (1983), symmetric turning: E R_k^2 = k E l^2 +
    // 2 (E l)^2 c / (1 - c) (k - (1 - c^k) / (1 - c)); the limit c -> 1 is
    // k E l^2 + (E l)^2 k (k - 1), a straight line
    double v;
    if (c == 1.0) {
      v = k * ml2 + ml * ml * k * (k - 1.0);
    } else {
      v = k * ml2 + 2.0 * ml * ml * c / (1.0 - c) * (k - (1.0 - std::pow(c, k)) / (1.0 - c));
    }
    msd[k - 1] = v;
  }
  return CLAB_OK;
}

int clab_crwsim(const double *step, int ns, const double *turn, int nt,
                double x0, double y0, double head, int h, int nsim,
                unsigned int seed, double *px, double *py) {
  if (ns < 1 || nt < 1 || h < 1 || nsim < 1) return CLAB_EARG;
  if (!std::isfinite(x0) || !std::isfinite(y0) || !std::isfinite(head)) return CLAB_EDOM;
  for (int i = 0; i < ns; ++i)
    if (!kutil::nonneg(step[i])) return CLAB_EDOM;
  for (int i = 0; i < nt; ++i)
    if (!std::isfinite(turn[i])) return CLAB_EDOM;
  cgeo::rng r(seed);
  for (int s = 0; s < nsim; ++s) {
    if (kutil::stop(4.0 * h)) return CLAB_EINTR;
    double cx = x0, cy = y0, hd = head;
    for (int k = 0; k < h; ++k) {
      // a turn and a step drawn independently from the observed ones
      hd += turn[r.below(nt)];
      const double l = step[r.below(ns)];
      cx += l * std::cos(hd);
      cy += l * std::sin(hd);
      px[static_cast<long>(k) * nsim + s] = cx;
      py[static_cast<long>(k) * nsim + s] = cy;
    }
  }
  return CLAB_OK;
}

int clab_rnnhcl(const double *x, const double *y, int n, const double *bx,
                const double *by, const double *bw, int nb, double bwid,
                double z, int minpts, int k, double c, int *iwork,
                double *work, double *thr, int *cl, int *ncl, double *cs) {
  if (n < 2 || nb < 1 || minpts < 2 || k < 0 || k > 2) return CLAB_EARG;
  if (!finite2(x, y, n) || !finite2(bx, by, nb) || !kutil::pos(bwid) ||
      !std::isfinite(z) || !kutil::pos(c))
    return CLAB_EDOM;
  kutil::ksum sb;
  for (int q = 0; q < nb; ++q) {
    if (!kutil::nonneg(bw[q])) return CLAB_EDOM;
    sb.add(bw[q]);
  }
  const double B = sb.val();
  if (!(B > 0.0)) return CLAB_EDOM;
  // the incidents' expected intensity under the baseline: n times the
  // baseline's normal-kernel density over its total, and from it each
  // point's random-expectation threshold
  const double f = 0.5 + z * std::sqrt((4.0 - pi) / (4.0 * pi * n));
  if (!(f > 0.0)) return CLAB_EDOM;
  double *lam = work, *mx = work + n, *my = mx + n;
  for (int i = 0; i < n; ++i) {
    kutil::ksum s;
    for (int q = 0; q < nb; ++q) {
      const double dx = x[i] - bx[q], dy = y[i] - by[q];
      s.add(bw[q] * std::exp(-(dx * dx + dy * dy) / (2.0 * bwid * bwid)));
    }
    lam[i] = n * s.val() / (2.0 * pi * bwid * bwid * B);
    thr[i] = lam[i] > 0.0 ? f / std::sqrt(lam[i]) : std::numeric_limits<double>::infinity();
  }
  // links: d_ij below the threshold at the mean of the two intensities
  // ponytail: O(n^2) pair scan; a grid needs a bound on the largest threshold
  int *deg = iwork, *ord = iwork + n, *asg = iwork + 2 * n, *mem = iwork + 3 * n;
  for (int i = 0; i < n; ++i) deg[i] = 0;
  auto linked = [&](int i, int j) {
    const double l = 0.5 * (lam[i] + lam[j]);
    if (!(l > 0.0)) return true;
    return std::hypot(x[i] - x[j], y[i] - y[j]) < f / std::sqrt(l);
  };
  for (int i = 0; i < n; ++i)
    for (int j = i + 1; j < n; ++j)
      if (linked(i, j)) {
        ++deg[i];
        ++deg[j];
      }
  for (int i = 0; i < n; ++i) {
    ord[i] = i;
    asg[i] = -1;
    cl[i] = -1;
  }
  kutil::hsort(ord, n, [deg](int a, int b) { return deg[a] > deg[b] || (deg[a] == deg[b] && a < b); });
  // groups: by decreasing links, each unclustered seed takes every
  // unclustered unit linked to it
  int ng = 0;
  for (int t = 0; t < n; ++t) {
    const int i = ord[t];
    if (deg[i] == 0 || asg[i] >= 0) continue;
    asg[i] = ng;
    for (int j = 0; j < n; ++j)
      if (j != i && asg[j] < 0 && linked(i, j)) asg[j] = ng;
    ++ng;
  }
  // keep groups of at least minpts, numbered by size (then first seed)
  int *size = ord;  // reuse: ng <= n
  for (int g = 0; g < ng; ++g) size[g] = 0;
  for (int i = 0; i < n; ++i)
    if (asg[i] >= 0) ++size[asg[i]];
  int *gid = deg;  // reuse
  for (int g = 0; g < ng; ++g) gid[g] = g;
  kutil::hsort(gid, ng, [size](int a, int b) { return size[a] > size[b] || (size[a] == size[b] && a < b); });
  int kept = 0;
  for (int q = 0; q < ng; ++q) {
    const int g = gid[q];
    if (size[g] < minpts) break;
    int m = 0;
    for (int i = 0; i < n; ++i)
      if (asg[i] == g) {
        mem[m] = i;
        mx[m] = x[i];
        my[m] = y[i];
        cl[i] = kept;
        ++m;
      }
    double *o = cs + 8 * kept;
    double st[10], wb[6], el[6];
    o[0] = m;
    clab_cgstat(mx, my, nullptr, m, 0, st);
    o[1] = st[1];
    o[2] = st[2];
    clab_webpt(mx, my, nullptr, m, nullptr, 1e-10, 1000, wb);
    o[3] = wb[0];
    o[4] = wb[1];
    if (clab_sdell(mx, my, nullptr, m, k, c, el) == CLAB_OK) {
      o[5] = el[2];
      o[6] = el[3];
      o[7] = el[4];
    } else {
      o[5] = o[6] = o[7] = kutil::nan();
    }
    ++kept;
  }
  *ncl = kept;
  return CLAB_OK;
}

}  // extern "C"
