// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, learning: random forests (Breiman 2001,
// doi:10.1023/A:1010933404324) for regression and for the probability of
// a 0/1 outcome (the mean of the outcome in a leaf; the variance split of a
// 0/1 outcome is the Gini split). Each column is cut once into at most 256
// bins at its own order statistics, so a node's best split comes from one
// pass over its rows into a histogram, and the forest depends on the
// columns only through their ranks: a strictly increasing transform of a
// column leaves every tree unchanged. Trees are grown and used one at a
// time, so memory does not grow with the number of trees.
#include <algorithm>
#include <cmath>

#include "crimlab.h"
#include "kutil.h"
#include "learn.h"

namespace {

const int NB = 256;

// The cut points of column j over the training rows: with at most NB
// distinct values, the midpoints between neighbours (every split a sort
// would find); otherwise NB - 1 order statistics. Returns their number.
int cuts(const double *xj, const int *tr, int mt, double *buf, double *cut) {
  for (int i = 0; i < mt; ++i) buf[i] = xj[tr[i]];
  std::sort(buf, buf + mt);
  int k = 0;
  for (int i = 1; i < mt; ++i) {
    if (buf[i] != buf[k]) buf[++k] = buf[i];
  }
  const int u = k + 1;  // distinct values, sorted, in buf[0..u)
  if (u <= NB) {
    for (int i = 0; i + 1 < u; ++i) cut[i] = buf[i] + 0.5 * (buf[i + 1] - buf[i]);
    return u - 1;
  }
  int c = 0;
  for (int b = 1; b < NB; ++b) {
    const double v = buf[static_cast<long>(b) * u / NB];
    if (c == 0 || v > cut[c - 1]) cut[c++] = v;
  }
  return c;
}

// the bin of v: the number of cut points below it
inline unsigned char bin(double v, const double *cut, int c) {
  return static_cast<unsigned char>(std::lower_bound(cut, cut + c, v) - cut);
}

struct Node {
  int j, b, l, r;  // split: column j, bins <= b to l; leaf: j = -1
  double v;        // leaf value
};

}  // namespace

namespace lrn {

long fwork(int n, int p, int nt) {
  // cut points, a sort buffer, the histogram, the sums of predictions,
  // the nodes (three doubles each, at most 2n + 1)
  return static_cast<long>(p) * (NB - 1) + n + 2L * NB + nt + 2L * n +
         5L * (2L * n + 1);
}
long fiwork(int n, int p, int nt) {
  // bins (bytes, four to an int), the bag, its counts, the stack, columns
  return (static_cast<long>(n + nt) * p + 3) / 4 + 3L * n + 1 + p;
}

int forest(const double *x, int n, int p, const double *y, const int *tr,
           int mt, const double *xp, int np, const int *te, int me,
           const Forest &f, double *out, double *oob, double *work,
           int *iwork) {
  if (mt < 2 || f.ntree < 1 || f.mtry < 1 || f.mtry > p || f.minnode < 1 ||
      !(f.frac > 0.0 && f.frac <= 1.0) || !(f.imb >= 0.0 && f.imb < 0.5))
    return CLAB_EARG;
  const int m = static_cast<int>(std::lround(f.frac * mt));
  if (m < (f.hon ? 2 : 1)) return CLAB_EARG;
  // honest trees (Wager and Athey 2018): the first half of a shuffled bag
  // places the splits, the second half the leaf means
  const int mg = f.hon ? m / 2 : m;
  double *cut = work, *buf = cut + static_cast<long>(p) * (NB - 1);
  double *hs = buf + n, *hc = hs + NB, *acc = hc + NB, *oa = acc + me,
         *on = oa + n;
  // 2m + 1 nodes in the doubles that follow (a node is 24 bytes)
  static_assert(sizeof(Node) == 3 * sizeof(double), "a node is three doubles");
  Node *node = reinterpret_cast<Node *>(on + n);
  double *es = on + n + 3L * (2L * n + 1), *ec = es + 2L * n + 1;
  unsigned char *bx = reinterpret_cast<unsigned char *>(iwork);
  int *bag = iwork + (static_cast<long>(n + me) * p + 3) / 4, *cnt = bag + n;
  int *stk = cnt + n, *col = stk + n + 1;
  // bins: training rows first (rows of x), then the rows to predict
  const long ld = static_cast<long>(n) + me;
  for (int j = 0; j < p; ++j) {
    const double *xj = x + static_cast<long>(j) * n;
    double *cj = cut + static_cast<long>(j) * (NB - 1);
    const int c = cuts(xj, tr, mt, buf, cj);
    unsigned char *bj = bx + j * ld;
    for (int i = 0; i < n; ++i) bj[i] = bin(xj[i], cj, c);
    const double *pj = xp + static_cast<long>(j) * np;
    for (int i = 0; i < me; ++i) bj[n + i] = bin(pj[te[i]], cj, c);
  }
  for (int i = 0; i < me; ++i) acc[i] = 0.0;
  if (oob != nullptr) {
    for (int i = 0; i < n; ++i) oa[i] = on[i] = 0.0;
  }
  for (int t = 0; t < f.ntree; ++t) {
    kutil::rng g(static_cast<std::uint64_t>(f.seed) * 1000003ULL + t);
    // the bag: m draws with replacement, or the first m of a shuffle
    for (int i = 0; i < n; ++i) cnt[i] = 0;
    if (f.repl) {
      for (int k = 0; k < m; ++k) bag[k] = tr[g.below(mt)];
    } else {
      for (int k = 0; k < mt; ++k) stk[k] = k;
      for (int k = 0; k < m; ++k) {
        const int r = k + g.below(mt - k), s = stk[k];
        stk[k] = stk[r];
        stk[r] = s;
        bag[k] = tr[stk[k]];
      }
    }
    for (int k = 0; k < m; ++k) ++cnt[bag[k]];
    if (f.hon) {
      for (int k = m - 1; k > 0; --k) {
        const int r = g.below(k + 1), s = bag[k];
        bag[k] = bag[r];
        bag[r] = s;
      }
    }
    // grow: a node is the range [lo, hi) of bag, kept on a stack
    int nn = 1, sp = 0;
    node[0].l = 0;
    node[0].r = mg;
    stk[sp++] = 0;
    while (sp > 0) {
      const int id = stk[--sp], lo = node[id].l, hi = node[id].r, sz = hi - lo;
      double s = 0.0;
      bool pure = true;
      for (int k = lo; k < hi; ++k) {
        s += y[bag[k]];
        pure = pure && y[bag[k]] == y[bag[lo]];
      }
      node[id].j = -1;
      node[id].v = s / sz;
      if (sz <= f.minnode || pure) continue;
      // mtry columns without replacement
      for (int j = 0; j < p; ++j) col[j] = j;
      double best = s * s / sz;
      int bj = -1, bb = 0;
      for (int k = 0; k < f.mtry; ++k) {
        const int r = k + g.below(p - k), c = col[k];
        col[k] = col[r];
        col[r] = c;
        const int j = col[k];
        const unsigned char *bv = bx + j * ld;
        for (int b = 0; b < NB; ++b) hs[b] = hc[b] = 0.0;
        for (int q = lo; q < hi; ++q) {
          hs[bv[bag[q]]] += y[bag[q]];
          hc[bv[bag[q]]] += 1.0;
        }
        // each child holds at least imb of the node (and one row)
        const double cmin = std::fmax(1.0, std::ceil(f.imb * sz));
        double sl = 0.0, nl = 0.0;
        for (int b = 0; b + 1 < NB; ++b) {
          sl += hs[b];
          nl += hc[b];
          if (nl < cmin) continue;
          if (sz - nl < cmin) break;
          const double sr = s - sl, nr = sz - nl;
          const double v = sl * sl / nl + sr * sr / nr;
          if (v > best * (1.0 + 1e-12)) {
            best = v;
            bj = j;
            bb = b;
          }
        }
      }
      if (bj < 0) continue;
      // partition the range: bins <= bb first
      const unsigned char *bv = bx + bj * ld;
      int a = lo, z = hi - 1;
      while (a <= z) {
        if (bv[bag[a]] <= bb) {
          ++a;
        } else {
          const int u = bag[a];
          bag[a] = bag[z];
          bag[z--] = u;
        }
      }
      node[id].j = bj;
      node[id].b = bb;
      const int l = nn++, r = nn++;
      node[l].l = lo;
      node[l].r = a;
      node[r].l = a;
      node[r].r = hi;
      node[id].l = l;
      node[id].r = r;
      stk[sp++] = r;
      stk[sp++] = l;
    }
    // honest: the estimation half's sums and counts at every node on
    // its rows' paths; a leaf without estimation rows takes the mean of
    // its deepest ancestor that has some (the root always has)
    if (f.hon) {
      for (int k = 0; k < nn; ++k) es[k] = ec[k] = 0.0;
      for (int q = mg; q < m; ++q) {
        const int i = bag[q];
        int k = 0;
        for (;;) {
          es[k] += y[i];
          ec[k] += 1.0;
          if (node[k].j < 0) break;
          k = bx[node[k].j * ld + i] <= node[k].b ? node[k].l : node[k].r;
        }
      }
    }
    // the value at the leaf of the row at column offset o of the bins
    auto val = [&](long o) {
      int k = 0, h = 0;
      for (;;) {
        if (f.hon && ec[k] > 0.0) h = k;
        if (node[k].j < 0) break;
        k = bx[node[k].j * ld + o] <= node[k].b ? node[k].l : node[k].r;
      }
      return f.hon ? es[h] / ec[h] : node[k].v;
    };
    // predict: the rows asked for, and the training rows out of the bag
    for (int i = 0; i < me; ++i) acc[i] += val(n + i);
    if (oob != nullptr) {
      for (int q = 0; q < mt; ++q) {
        const int i = tr[q];
        if (cnt[i] > 0) continue;
        oa[i] += val(i);
        on[i] += 1.0;
      }
    }
  }
  for (int i = 0; i < me; ++i) out[i] = acc[i] / f.ntree;
  if (oob != nullptr) {
    for (int q = 0; q < mt; ++q) {
      const int i = tr[q];
      oob[q] = on[i] > 0.0 ? oa[i] / on[i] : kutil::nan();
    }
  }
  return CLAB_OK;
}

}  // namespace lrn

extern "C" {

int clab_rfor(const double *x, int n, int p, const double *y,
              const double *xt, int nt, int ntree, int mtry, int minnode,
              double frac, int repl, int hon, double imb, unsigned int seed,
              double *pred, double *oob, double *work, int *iwork) {
  if (n < 2 || p < 1 || nt < 0) return CLAB_EARG;
  for (long k = 0; k < static_cast<long>(n) * p; ++k) {
    if (!std::isfinite(x[k])) return CLAB_EDOM;
  }
  for (long k = 0; k < static_cast<long>(nt) * p; ++k) {
    if (!std::isfinite(xt[k])) return CLAB_EDOM;
  }
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(y[i])) return CLAB_EDOM;
  }
  // rows: 0..n-1 of x to train on, 0..nt-1 of xt to predict; the index
  // lists sit at the front of iwork, the forest's scratch after them
  int *tr = iwork, *te = tr + n;
  for (int i = 0; i < n; ++i) tr[i] = i;
  for (int i = 0; i < nt; ++i) te[i] = i;
  const lrn::Forest f = {ntree, mtry, minnode, frac, repl, seed, hon, imb};
  return lrn::forest(x, n, p, y, tr, n, xt, nt, te, nt, f, pred, oob, work,
                     te + nt);
}

}  // extern "C"
