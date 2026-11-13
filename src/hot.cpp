// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, hot spots: the kernel intensity surface with edge
// corrections, DBSCAN clusters, the persistence of hot grid cells over
// periods, and the empirical-Bayes shrinkage of hot-spot counts.
#include <algorithm>
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

namespace {

constexpr double pi = 3.14159265358979323846;

// the standard normal distribution function, from erfc (no cancellation in
// either tail)
double phi(double z) { return 0.5 * std::erfc(-z / std::sqrt(2.0)); }

// the share of the product normal kernel centred at (u, v) inside win
double emass(double u, double v, double hx, double hy, const double *win) {
  return (phi((win[1] - u) / hx) - phi((win[0] - u) / hx)) *
         (phi((win[3] - v) / hy) - phi((win[2] - v) / hy));
}

}  // namespace

int clab_kdsurf(const double *x, const double *y, const double *w, int n, double hx, double hy,
                const double *ax, int na, const double *ay, int nb, int grid, const double *win,
                int edge, int loo, double *out) {
  if (n < 1 || na < 1 || nb < 0 || edge < 0 || edge > 2 || grid < 0 || grid > 1) {
    return CLAB_EARG;
  }
  if (!grid && nb != na) return CLAB_EARG;
  if (loo && (grid || na != n)) return CLAB_EARG;
  if (edge > 0 && win == nullptr) return CLAB_EARG;
  if (!kutil::pos(hx) || !kutil::pos(hy)) return CLAB_EDOM;
  if (edge > 0 && !(win[0] < win[1] && win[2] < win[3])) return CLAB_EDOM;
  for (int j = 0; j < n; ++j) {
    if (!std::isfinite(x[j]) || !std::isfinite(y[j])) return CLAB_EDOM;
    if (w != nullptr && !kutil::nonneg(w[j])) return CLAB_EDOM;
  }
  // each point's mass: its weight, over its own edge mass for Diggle's
  // (Jones's) correction
  std::vector<double> m(n);
  for (int j = 0; j < n; ++j) {
    m[j] = w == nullptr ? 1.0 : w[j];
    if (edge == 2) {
      const double e = emass(x[j], y[j], hx, hy, win);
      if (!(e > 0)) return CLAB_EDOM;
      m[j] /= e;
    }
  }
  const double c = 1.0 / (2 * pi * hx * hy);
  if (grid) {
    // separable: z[a, b] = sum_j m_j kx(a, j) ky(b, j)
    std::vector<double> ky(static_cast<std::size_t>(nb) * n), kx(n);
    for (int b = 0; b < nb; ++b) {
      for (int j = 0; j < n; ++j) {
        const double d = (ay[b] - y[j]) / hy;
        ky[static_cast<std::size_t>(b) * n + j] = std::exp(-0.5 * d * d);
      }
    }
    for (int a = 0; a < na; ++a) {
      if (kutil::stop(static_cast<double>(n) * nb)) return CLAB_EINTR;
      for (int j = 0; j < n; ++j) {
        const double d = (ax[a] - x[j]) / hx;
        kx[j] = m[j] * std::exp(-0.5 * d * d);
      }
      for (int b = 0; b < nb; ++b) {
        const double *k = ky.data() + static_cast<std::size_t>(b) * n;
        kutil::ksum s;
        for (int j = 0; j < n; ++j) s.add(kx[j] * k[j]);
        out[a + static_cast<std::size_t>(b) * na] = c * s.val();
      }
    }
  } else {
    for (int a = 0; a < na; ++a) {
      if (kutil::stop(n)) return CLAB_EINTR;
      kutil::ksum s;
      for (int j = 0; j < n; ++j) {
        if (loo && j == a) continue;
        const double dx = (ax[a] - x[j]) / hx, dy = (ay[a] - y[j]) / hy;
        s.add(m[j] * std::exp(-0.5 * (dx * dx + dy * dy)));
      }
      out[a] = c * s.val();
    }
  }
  if (edge == 1) {
    for (int a = 0; a < na; ++a) {
      for (int b = 0; b < (grid ? nb : 1); ++b) {
        const double u = ax[a], v = grid ? ay[b] : ay[a];
        const double e = emass(u, v, hx, hy, win);
        const std::size_t q = a + static_cast<std::size_t>(b) * na;
        out[q] = e > 0 ? out[q] / e : kutil::nan();
      }
    }
  }
  return CLAB_OK;
}

int clab_dbscn(const double *x, const double *y, int n, double eps, int minpts, int border, int *cl,
               int *core, int *ncl) {
  if (n < 1 || minpts < 1 || border < 0 || border > 1) return CLAB_EARG;
  if (!kutil::pos(eps)) return CLAB_EDOM;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(x[i]) || !std::isfinite(y[i])) return CLAB_EDOM;
  }
  std::vector<int> rp(n + 1), iw(2 * n + 2);
  int nnz = 0;
  clab_swband(x, y, n, eps, 0, 0, rp.data(), nullptr, nullptr, &nnz, iw.data());
  std::vector<int> ci(nnz > 0 ? nnz : 1);
  std::vector<double> d(nnz > 0 ? nnz : 1);
  clab_swband(x, y, n, eps, 0, nnz, rp.data(), ci.data(), d.data(), &nnz, iw.data());
  // a core point has at least minpts points within eps, itself included
  for (int i = 0; i < n; ++i) {
    core[i] = rp[i + 1] - rp[i] + 1 >= minpts;
    cl[i] = -1;
  }
  // clusters grown from the core points in index order, breadth first;
  // a border point joins the first cluster that reaches it
  std::vector<int> queue;
  int k = 0;
  for (int i = 0; i < n; ++i) {
    if (!core[i] || cl[i] >= 0) continue;
    cl[i] = k;
    queue.assign(1, i);
    for (std::size_t h = 0; h < queue.size(); ++h) {
      const int a = queue[h];
      for (int q = rp[a]; q < rp[a + 1]; ++q) {
        const int b = ci[q];
        if (cl[b] >= 0) continue;
        if (core[b]) {
          cl[b] = k;
          queue.push_back(b);
        } else if (border) {
          cl[b] = k;
        }
      }
    }
    ++k;
  }
  *ncl = k;
  return CLAB_OK;
}

int clab_hsper(const double *x, const double *y, const int *prd, int n, int np, double cell,
               double x0, double y0, double top, int *cid, double *ix, double *iy, int *nc,
               int *hot, double *per, double *cs) {
  if (n < 1 || np < 1) return CLAB_EARG;
  if (!kutil::pos(cell) || !(top > 0 && top <= 1) || !std::isfinite(x0) || !std::isfinite(y0)) {
    return CLAB_EDOM;
  }
  std::vector<double> gx(n), gy(n);
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(x[i]) || !std::isfinite(y[i])) return CLAB_EDOM;
    if (prd[i] < 0 || prd[i] >= np) return CLAB_EARG;
    gx[i] = std::floor((x[i] - x0) / cell);
    gy[i] = std::floor((y[i] - y0) / cell);
    if (!(std::fabs(gx[i]) < 1e15 && std::fabs(gy[i]) < 1e15)) return CLAB_EDOM;
  }
  // cells in order of column, then row
  std::vector<int> o(n);
  for (int i = 0; i < n; ++i) o[i] = i;
  std::sort(o.begin(), o.end(), [&](int a, int b) {
    return gx[a] == gx[b] ? (gy[a] == gy[b] ? a < b : gy[a] < gy[b]) : gx[a] < gx[b];
  });
  int k = -1;
  for (int q = 0; q < n; ++q) {
    const int i = o[q];
    if (q == 0 || gx[i] != gx[o[q - 1]] || gy[i] != gy[o[q - 1]]) {
      ++k;
      ix[k] = gx[i];
      iy[k] = gy[i];
    }
    cid[i] = k;
  }
  const int m = k + 1;
  *nc = m;
  std::vector<double> cnt(static_cast<std::size_t>(m) * np, 0.0);
  for (int i = 0; i < n; ++i) cnt[cid[i] + static_cast<std::size_t>(prd[i]) * m] += 1;
  const int kh = std::max(1, static_cast<int>(std::ceil(top * m)));
  std::vector<double> srt(m);
  for (int p = 0; p < np; ++p) {
    const double *c = cnt.data() + static_cast<std::size_t>(p) * m;
    int *h = hot + static_cast<std::size_t>(p) * m;
    double tot = 0;
    int act = 0;
    for (int j = 0; j < m; ++j) {
      srt[j] = c[j];
      tot += c[j];
      act += c[j] > 0;
    }
    std::sort(srt.begin(), srt.end(), [](double a, double b) { return a > b; });
    const double cut = srt[kh - 1];
    int above = 0, tied = 0;
    for (int j = 0; j < m; ++j) {
      above += c[j] > cut;
      tied += c[j] == cut && c[j] > 0;
    }
    // a tie that would make more than twice kh cells hot singles none out:
    // the tied cells are left out
    const bool drop = above + tied > 2 * kh;
    int nh = 0;
    double inhot = 0;
    for (int j = 0; j < m; ++j) {
      h[j] = c[j] > cut || (!drop && c[j] == cut && c[j] > 0);
      nh += h[j];
      if (h[j]) inhot += c[j];
    }
    // the fewest cells holding half the incidents
    double run = 0;
    int half = 0;
    while (tot > 0 && run < 0.5 * tot) run += srt[half++];
    double *r = per + static_cast<std::size_t>(p) * 9;
    r[0] = tot;
    r[1] = act;
    r[2] = tot > 0 ? half : kutil::nan();
    r[3] = tot > 0 ? static_cast<double>(half) / m : kutil::nan();
    r[4] = nh;
    r[5] = std::max(0, nh - kh);
    r[6] = drop ? tied : 0;
    r[7] = tot > 0 ? inhot / tot : kutil::nan();
    // Jaccard overlap with the period before; unknown when either hot set
    // is empty
    r[8] = kutil::nan();
    if (p > 0) {
      const int *g = hot + static_cast<std::size_t>(p - 1) * m;
      int un = 0, in = 0, ng = 0;
      for (int j = 0; j < m; ++j) {
        un += h[j] || g[j];
        in += h[j] && g[j];
        ng += g[j];
      }
      if (nh > 0 && ng > 0) r[8] = static_cast<double>(in) / un;
    }
  }
  // per cell: incidents, periods hot, longest run of successive hot periods
  for (int j = 0; j < m; ++j) {
    double tot = 0;
    int nh = 0, len = 0, best = 0;
    for (int p = 0; p < np; ++p) {
      tot += cnt[j + static_cast<std::size_t>(p) * m];
      if (hot[j + static_cast<std::size_t>(p) * m]) {
        ++nh;
        best = std::max(best, ++len);
      } else {
        len = 0;
      }
    }
    cs[3 * j] = tot;
    cs[3 * j + 1] = nh;
    cs[3 * j + 2] = best;
  }
  return CLAB_OK;
}

int clab_ebshr(const double *y, const double *w, int n, double v, double b, double *fit,
               double *shr, double *fall) {
  if (n < 2) return CLAB_EARG;
  if (!kutil::nonneg(v) || !(std::isnan(b) || (b >= 0 && b <= 1))) return CLAB_EDOM;
  kutil::ksum sw, sy;
  for (int i = 0; i < n; ++i) {
    const double wi = w == nullptr ? 1.0 : w[i];
    if (!std::isfinite(y[i]) || !kutil::nonneg(wi)) return CLAB_EDOM;
    sw.add(wi);
    sy.add(wi * y[i]);
  }
  const double W = sw.val();
  if (!(W > 0)) return CLAB_EDOM;
  const double ybar = sy.val() / W;
  kutil::ksum st;
  for (int i = 0; i < n; ++i) {
    const double wi = w == nullptr ? 1.0 : w[i];
    st.add(wi * (y[i] - ybar) * (y[i] - ybar));
  }
  const double tot = st.val() / W;
  const double sig = std::max(tot - v, 0.0);
  // B* = S_e / (S_e + S_theta), the minimiser of (1 - B)^2 S_e + B^2 S_theta
  const double bs = v + sig > 0 ? v / (v + sig) : 0.0;
  const double bb = std::isnan(b) ? bs : b;
  for (int i = 0; i < n; ++i) {
    shr[i] = (1 - bb) * y[i] + bb * ybar;
    fall[i] = bb * (y[i] - ybar);
  }
  fit[0] = ybar;
  fit[1] = tot;
  fit[2] = v;
  fit[3] = sig;
  fit[4] = bs;
  fit[5] = bb;
  return CLAB_OK;
}
