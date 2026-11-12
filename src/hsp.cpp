// SPDX-License-Identifier: AGPL-3.0-or-later
//
// libcrimlab, hot spots of points: the fuzzy mode, K-means partitioning,
// nearest-neighbour hierarchical clustering and search-circle hot-spot
// areas, each with the clusters' centres and standard deviational ellipses.
#include <cmath>
#include <limits>

#include "cgeo.h"
#include "crimlab.h"
#include "kutil.h"

namespace {

bool finite2(const double *x, const double *y, int n) {
  for (int i = 0; i < n; ++i)
    if (!std::isfinite(x[i]) || !std::isfinite(y[i]) ||
        std::fabs(x[i]) > 1e150 || std::fabs(y[i]) > 1e150)
      return false;
  return true;
}

void bbox(const double *x, const double *y, int n, double *b) {
  b[0] = b[1] = x[0];
  b[2] = b[3] = y[0];
  for (int i = 1; i < n; ++i) {
    if (x[i] < b[0]) b[0] = x[i];
    if (x[i] > b[1]) b[1] = x[i];
    if (y[i] < b[2]) b[2] = y[i];
    if (y[i] > b[3]) b[3] = y[i];
  }
}

// visit every point j with |p_j - (px, py)| <= r (closed) or < r (open)
template <class F>
void within(const cgeo::grid &g, const double *x, const double *y, double px,
            double py, double r, bool closed, F f) {
  const int a0 = g.cx(px - r), a1 = g.cx(px + r);
  const int b0 = g.cy(py - r), b1 = g.cy(py + r);
  const double r2 = r * r;
  for (int b = b0; b <= b1; ++b)
    for (int a = a0; a <= a1; ++a) {
      const int c = b * g.gx + a;
      for (int t = g.start[c]; t < g.start[c + 1]; ++t) {
        const int j = g.ord[t];
        const double dx = x[j] - px, dy = y[j] - py;
        const double d = dx * dx + dy * dy;
        if (closed ? d <= r2 : d < r2) f(j);
      }
    }
}

// centre (mean and Weber point) and ellipse of the points idx[0..m) of
// (x, y): row of 8: count, mean x, mean y, Weber x, Weber y, a, b, theta.
// work: 2m doubles.
void summ(const double *x, const double *y, const int *idx, int m, int k,
          double c, double *work, double *row) {
  double *gx = work, *gy = work + m;
  for (int t = 0; t < m; ++t) {
    gx[t] = x[idx[t]];
    gy[t] = y[idx[t]];
  }
  double mo[6] = {0, 0, 0, 0, 0, 0};
  cgeo::moments(gx, gy, nullptr, m, mo);  // m >= 1 finite points: fills mo
  row[0] = m;
  row[1] = mo[1];
  row[2] = mo[2];
  double sc = std::sqrt((mo[3] + mo[4]) / m), wx = 0.0, wy = 0.0, wo[4];
  if (!(sc > 0.0)) sc = 1.0;
  cgeo::weber(gx, gy, nullptr, m, mo[1], mo[2], sc, 1e-12, 1000, &wx, &wy, wo);
  row[3] = mo[1] + wx;
  row[4] = mo[2] + wy;
  const double den = m - k;
  if (den > 0.0) {
    double l1, l2, th;
    cgeo::eig2(mo[3] / den, mo[4] / den, mo[5] / den, &l1, &l2, &th);
    row[5] = c * std::sqrt(l1);
    row[6] = c * std::sqrt(l2);
    row[7] = th;
  } else {
    row[5] = row[6] = row[7] = kutil::nan();
  }
}

}  // namespace

extern "C" {

int clab_fzmode(const double *x, const double *y, int n, double r, int *iwork,
                int *cnt) {
  if (n < 1) return CLAB_EARG;
  if (!kutil::nonneg(r) || !finite2(x, y, n)) return CLAB_EDOM;
  double b[4];
  bbox(x, y, n, b);
  const double side = std::fmax(b[1] - b[0], b[3] - b[2]);
  cgeo::grid g;
  g.build(x, y, n, b[0], b[1], b[2], b[3],
          std::fmax(r, side / std::sqrt(static_cast<double>(n))), iwork);
  for (int t = 0; t < n; ++t) {  // in cell order, for memory locality
    const int i = g.ord[t];
    int c = 0;
    within(g, x, y, x[i], y[i], r, true, [&c](int) { ++c; });
    cnt[i] = c;
  }
  return CLAB_OK;
}

int clab_kmpart(const double *x, const double *y, int n, int k, int init,
                const double *par, unsigned int seed, int upd, int maxit,
                double *cen, int *cl, int *iwork, double *work, double *out) {
  if (n < 1 || k < 1 || k > n || init < 0 || init > 2 || upd < 0 ||
      upd > 1 || maxit < 0)
    return CLAB_EARG;
  if (!finite2(x, y, n)) return CLAB_EDOM;
  double *cx = cen, *cy = cen + k;
  if (init == 0) {
    if (!finite2(cx, cy, k)) return CLAB_EDOM;
  } else if (init == 1) {
    // k-means++ (Arthur and Vassilvitskii 2007): the first centre uniform,
    // each next with probability proportional to D^2
    cgeo::rng r(seed);
    double *d2 = work;
    int f = r.below(n);
    cx[0] = x[f];
    cy[0] = y[f];
    for (int i = 0; i < n; ++i) {
      const double dx = x[i] - cx[0], dy = y[i] - cy[0];
      d2[i] = dx * dx + dy * dy;
    }
    for (int j = 1; j < k; ++j) {
      kutil::ksum s;
      for (int i = 0; i < n; ++i) s.add(d2[i]);
      const double tot = s.val();
      int pick = n - 1;
      if (tot > 0.0) {
        const double u = r.unif() * tot;
        double acc = 0.0;
        for (int i = 0; i < n; ++i) {
          acc += d2[i];
          if (acc > u && d2[i] > 0.0) {
            pick = i;
            break;
          }
        }
        while (d2[pick] == 0.0) --pick;  // rounding at the end of the walk
      } else {
        return CLAB_EDOM;  // fewer than k distinct points
      }
      cx[j] = x[pick];
      cy[j] = y[pick];
      for (int i = 0; i < n; ++i) {
        const double dx = x[i] - cx[j], dy = y[i] - cy[j];
        const double v = dx * dx + dy * dy;
        if (v < d2[i]) d2[i] = v;
      }
    }
  } else {
    // separated grid seeds: cells of side par[0] over the bounding box,
    // taken by decreasing count (then cell number), a cell's centre
    // accepted when at least par[1] from every accepted seed
    const double side = par[0], sep = par[1];
    if (!kutil::pos(side) || !kutil::nonneg(sep)) return CLAB_EDOM;
    double b[4];
    bbox(x, y, n, b);
    const double fx = std::floor((b[1] - b[0]) / side) + 1.0;
    const double fy = std::floor((b[3] - b[2]) / side) + 1.0;
    if (fx * fy > 1e15) return CLAB_EDOM;
    auto cell = [&](int i) {
      return std::floor((x[i] - b[0]) / side) +
             fx * std::floor((y[i] - b[2]) / side);
    };
    int *ord = iwork, *run = iwork + n;
    for (int i = 0; i < n; ++i) ord[i] = i;
    cgeo::hsort(ord, n, [&](int a, int c) {
      const double u = cell(a), v = cell(c);
      return u < v || (u == v && a < c);
    });
    int nr = 0;
    for (int t = 0; t < n; ++t)
      if (t == 0 || cell(ord[t]) != cell(ord[t - 1])) run[nr++] = t;
    // run[q] start of the q-th occupied cell; sort cells by count desc
    double *cntv = work;
    for (int q = 0; q < nr; ++q)
      cntv[q] = (q + 1 < nr ? run[q + 1] : n) - run[q];
    double *cc = work + n;  // cell code of each run
    for (int q = 0; q < nr; ++q) cc[q] = cell(ord[run[q]]);
    for (int q = 0; q < nr; ++q) run[q] = q;
    cgeo::hsort(run, nr, [&](int a, int c) {
      return cntv[a] > cntv[c] || (cntv[a] == cntv[c] && cc[a] < cc[c]);
    });
    int got = 0;
    for (int t = 0; t < nr && got < k; ++t) {
      const double code = cc[run[t]];
      const double iy = std::floor(code / fx), ix = code - fx * iy;
      const double sx = b[0] + (ix + 0.5) * side, sy = b[2] + (iy + 0.5) * side;
      bool ok = true;
      for (int j = 0; j < got && ok; ++j)
        ok = std::hypot(sx - cx[j], sy - cy[j]) >= sep;
      if (ok) {
        cx[got] = sx;
        cy[got] = sy;
        ++got;
      }
    }
    if (got < k) return CLAB_EDOM;
  }
  // Lloyd's iterations (Lloyd 1982, doi:10.1109/TIT.1982.1056489): assign
  // to the nearest centre (lowest index on ties), then move each centre to
  // its points' mean (upd 0) or Weber point (upd 1); stop when no point
  // changes cluster. An emptied cluster keeps its centre.
  for (int i = 0; i < n; ++i) cl[i] = -1;
  int it = 0, conv = 0;
  int *mem = iwork;
  for (;;) {
    int changed = 0;
    for (int i = 0; i < n; ++i) {
      int bj = 0;
      double bd = std::numeric_limits<double>::infinity();
      for (int j = 0; j < k; ++j) {
        const double dx = x[i] - cx[j], dy = y[i] - cy[j];
        const double d = dx * dx + dy * dy;
        if (d < bd) {
          bd = d;
          bj = j;
        }
      }
      if (cl[i] != bj) {
        cl[i] = bj;
        ++changed;
      }
    }
    if (changed == 0) {
      conv = 1;
      break;
    }
    if (it >= maxit) break;
    ++it;
    for (int j = 0; j < k; ++j) {
      int m = 0;
      for (int i = 0; i < n; ++i)
        if (cl[i] == j) mem[m++] = i;
      if (m == 0) continue;
      double row[8];
      if (upd == 0) {
        kutil::ksum sx, sy;
        for (int t = 0; t < m; ++t) {
          sx.add(x[mem[t]]);
          sy.add(y[mem[t]]);
        }
        cx[j] = sx.val() / m;
        cy[j] = sy.val() / m;
      } else {
        summ(x, y, mem, m, 0, 1.0, work, row);
        cx[j] = row[3];
        cy[j] = row[4];
      }
    }
  }
  // summaries about the final centres
  kutil::ksum tw, td;
  for (int j = 0; j < k; ++j) out[4 + 3 * j] = out[5 + 3 * j] = out[6 + 3 * j] = 0;
  for (int j = 0; j < k; ++j) {
    kutil::ksum w, s;
    int m = 0;
    for (int i = 0; i < n; ++i) {
      if (cl[i] != j) continue;
      const double dx = x[i] - cx[j], dy = y[i] - cy[j];
      w.add(dx * dx + dy * dy);
      s.add(std::hypot(dx, dy));
      ++m;
    }
    out[4 + 3 * j] = m;
    out[5 + 3 * j] = w.val();
    out[6 + 3 * j] = s.val();
    tw.add(w.val());
    td.add(s.val());
  }
  out[0] = tw.val();
  out[1] = td.val();
  out[2] = it;
  out[3] = conv;
  return CLAB_OK;
}

int clab_nnhclu(const double *x, const double *y, int n, double area,
                double fixd, double z, int minpts, int minhi, int maxord,
                int k, double c, int *iwork, double *work, int *cl, int *ncl,
                double *cs) {
  if (n < 2 || minpts < 2 || minhi < 2 || maxord < 1 || k < 0 || k > 2)
    return CLAB_EARG;
  if (!finite2(x, y, n) || !kutil::nonneg(fixd) || !std::isfinite(z) ||
      !kutil::pos(c))
    return CLAB_EDOM;
  if (fixd == 0.0 && !kutil::pos(area)) return CLAB_EDOM;
  // units of the current order: ux, uy (m of them); unit u of order o
  // corresponds to cluster number base + u of order o - 1 (order 1: points)
  double *ux = work, *uy = work + n, *sw = work + 2 * n;  // sw: 2n scratch
  int *deg = iwork, *ord = iwork + n, *asg = iwork + 2 * n,
      *mem = iwork + 3 * n, *giw = iwork + 4 * n;  // giw: 2n + 2
  for (int i = 0; i < n; ++i) {
    ux[i] = x[i];
    uy[i] = y[i];
  }
  for (int t = 0; t < n * maxord; ++t) cl[t] = -1;
  int m = n, nc = 0, prevbase = 0;
  for (int o = 0; o < maxord; ++o) {
    if (o > 0 && m < 4) break;  // too few units to cluster further
    const double thr =
        fixd > 0.0 ? fixd
                   : 0.5 * std::sqrt(area / m) +
                         z * std::sqrt((4.0 - cgeo::pi) * area /
                                       (4.0 * cgeo::pi)) /
                             m;
    if (!(thr > 0.0)) break;
    double b[4];
    bbox(ux, uy, m, b);
    cgeo::grid g;
    g.build(ux, uy, m, b[0], b[1], b[2], b[3], thr, giw);
    // neighbours closer than the threshold
    for (int t = 0; t < m; ++t) {  // in cell order, for memory locality
      const int i = g.ord[t];
      int d = 0;
      within(g, ux, uy, ux[i], uy[i], thr, false, [&d, i](int j) {
        if (j != i) ++d;
      });
      deg[i] = d;
      ord[i] = i;
      asg[i] = -1;
    }
    cgeo::hsort(ord, m, [deg](int a, int bb) {
      return deg[a] > deg[bb] || (deg[a] == deg[bb] && a < bb);
    });
    // seeds in that order take every unassigned unit within the threshold
    int nfound = 0;
    for (int t = 0; t < m && deg[ord[t]] > 0; ++t) {
      const int s = ord[t];
      if (asg[s] >= 0) continue;
      const int id = nfound++;
      asg[s] = id;
      within(g, ux, uy, ux[s], uy[s], thr, false, [&asg, id](int j) {
        if (asg[j] < 0) asg[j] = id;
      });
    }
    // keep the clusters with enough members, numbered by size (desc), then
    // by order of their seed
    const int need = o == 0 ? minpts : minhi;
    int *sz = deg;  // reuse
    for (int q = 0; q < nfound; ++q) sz[q] = 0;
    for (int i = 0; i < m; ++i)
      if (asg[i] >= 0) ++sz[asg[i]];
    for (int q = 0; q < nfound; ++q) ord[q] = q;
    cgeo::hsort(ord, nfound, [sz](int a, int bb) {
      return sz[a] > sz[bb] || (sz[a] == sz[bb] && a < bb);
    });
    int kept = 0;
    int *newid = mem;  // mem doubles as the map old id -> new id
    for (int q = 0; q < nfound; ++q) newid[q] = -1;
    for (int t = 0; t < nfound; ++t)
      if (sz[ord[t]] >= need) newid[ord[t]] = kept++;
    if (kept == 0) break;
    const int base = nc;
    // membership of the original points at this order
    for (int i = 0; i < n; ++i) {
      int u = o == 0 ? i : cl[(o - 1) * n + i];
      if (o > 0) u = u >= 0 ? u - prevbase : -1;
      const int a = (u >= 0 && asg[u] >= 0) ? newid[asg[u]] : -1;
      cl[o * n + i] = a >= 0 ? base + a : -1;
    }
    // summaries: of the member units (points at order 1, lower-order
    // centres above); Weber centres become the next order's units
    for (int q = 0; q < kept; ++q) {
      double *row = cs + 11 * (base + q);
      int cnt = 0, pts = 0;
      int *list = ord;  // reuse
      for (int u = 0; u < m; ++u)
        if (asg[u] >= 0 && newid[asg[u]] == q) list[cnt++] = u;
      for (int i = 0; i < n; ++i)
        if (cl[o * n + i] == base + q) ++pts;
      summ(ux, uy, list, cnt, k, c, sw, row + 2);
      row[0] = o + 1;
      row[1] = pts;
      row[10] = thr;
    }
    nc += kept;
    for (int q = 0; q < kept; ++q) {
      ux[q] = cs[11 * (base + q) + 5];
      uy[q] = cs[11 * (base + q) + 6];
    }
    prevbase = base;
    m = kept;
  }
  *ncl = nc;
  return CLAB_OK;
}

int clab_stacel(const double *x, const double *y, int n, const double *win,
                int nx, int ny, double r, int top, int minpts, int k, double c,
                int *iwork, double *work, int *cl, int *ncl, double *cs) {
  if (n < 1 || nx < 2 || ny < 2 || nx > 10000 || ny > 10000 || top < 1 ||
      minpts < 2 || k < 0 || k > 2)
    return CLAB_EARG;
  if (!finite2(x, y, n) || !kutil::pos(r) || !kutil::pos(c)) return CLAB_EDOM;
  if (!(std::isfinite(win[0]) && std::isfinite(win[1]) &&
        std::isfinite(win[2]) && std::isfinite(win[3]) && win[0] < win[1] &&
        win[2] < win[3]))
    return CLAB_EDOM;
  const int nn = nx * ny;
  // iwork: grid 2n + 2, nodes nn, union-find top, first n, lists
  // max(n, 2 top); work: nn counts, 2n scratch
  int *giw = iwork, *nodes = iwork + 2 * n + 2, *par = nodes + nn,
      *first = par + top, *lst = first + n;
  double *cntv = work, *sw = work + nn;
  double b[4];
  bbox(x, y, n, b);
  cgeo::grid g;
  g.build(x, y, n, b[0], b[1], b[2], b[3],
          std::fmax(r, std::fmax(b[1] - b[0], b[3] - b[2]) /
                           std::sqrt(static_cast<double>(n))),
          giw);
  auto nodex = [&](int q) {
    return win[0] + (win[1] - win[0]) * (q % nx) / (nx - 1.0);
  };
  auto nodey = [&](int q) {
    return win[2] + (win[3] - win[2]) * (q / nx) / (ny - 1.0);
  };
  for (int q = 0; q < nn; ++q) {
    int cc = 0;
    within(g, x, y, nodex(q), nodey(q), r, true, [&cc](int) { ++cc; });
    cntv[q] = cc;
    nodes[q] = q;
  }
  cgeo::hsort(nodes, nn, [cntv](int a, int bb) {
    return cntv[a] > cntv[bb] || (cntv[a] == cntv[bb] && a < bb);
  });
  int nk = 0;
  while (nk < top && nk < nn && cntv[nodes[nk]] >= 2) ++nk;
  // circles sharing a point are merged (union-find, smaller rank as root)
  for (int t = 0; t < nk; ++t) par[t] = t;
  auto root = [par](int a) {
    while (par[a] != a) a = par[a] = par[par[a]];
    return a;
  };
  for (int i = 0; i < n; ++i) first[i] = -1;
  for (int t = 0; t < nk; ++t) {
    const int q = nodes[t];
    within(g, x, y, nodex(q), nodey(q), r, true, [&](int j) {
      if (first[j] < 0) {
        first[j] = t;
      } else {
        const int a = root(first[j]), bb = root(t);
        // the better-ranked circle (smaller rank) stays the root
        par[a < bb ? bb : a] = a < bb ? a : bb;
      }
    });
  }
  // a cluster is the set of points in any circle of a component; count
  // them by the component's root (the best-ranked circle)
  int *sz = lst;  // nk ints
  for (int t = 0; t < nk; ++t) sz[t] = 0;
  for (int i = 0; i < n; ++i) {
    cl[i] = first[i] >= 0 ? root(first[i]) : -1;
    if (cl[i] >= 0) ++sz[cl[i]];
  }
  int *rk = lst + nk;  // nk ints: roots by size desc, then rank
  int nr = 0;
  for (int t = 0; t < nk; ++t)
    if (root(t) == t && sz[t] >= minpts) rk[nr++] = t;
  cgeo::hsort(rk, nr, [sz](int a, int bb) {
    return sz[a] > sz[bb] || (sz[a] == sz[bb] && a < bb);
  });
  // map roots to cluster numbers (the node counts are no longer needed)
  double *rid = cntv;
  for (int t = 0; t < nk; ++t) rid[t] = -1;
  for (int q = 0; q < nr; ++q) rid[rk[q]] = q;
  for (int i = 0; i < n; ++i) cl[i] = cl[i] >= 0 ? static_cast<int>(rid[cl[i]]) : -1;
  for (int q = 0; q < nr; ++q) {
    int m = 0;
    int *mem = lst;  // sizes no longer needed
    for (int i = 0; i < n; ++i)
      if (cl[i] == q) mem[m++] = i;
    summ(x, y, mem, m, k, c, sw, cs + 8 * q);
  }
  *ncl = nr;
  return CLAB_OK;
}

}  // extern "C"
