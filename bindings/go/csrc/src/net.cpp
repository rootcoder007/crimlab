// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, networks: co-offending graphs and the measures of criminal
// networks on a graph in compressed sparse rows (CSR). Node i's
// neighbours are adj[off[i]] .. adj[off[i + 1] - 1]; an undirected graph
// holds each edge in both rows. No allocation: the caller owns every
// array, the work arrays included.
#include <climits>
#include <cmath>
#include <cstdint>

#include "crimlab.h"
#include "kutil.h"

namespace {

// A graph the kernels can read: n >= 1 nodes, offsets that never fall and
// every neighbour a node.
bool csrok(int n, const int *off, const int *adj) {
  if (n < 1 || off[0] != 0) return false;
  for (int i = 0; i < n; ++i) {
    if (off[i + 1] < off[i]) return false;
  }
  for (int p = 0; p < off[n]; ++p) {
    if (adj[p] < 0 || adj[p] >= n) return false;
  }
  return true;
}

// Weights, when given, positive and finite.
bool wtok(const double *wt, int nnz) {
  if (wt == nullptr) return true;
  for (int p = 0; p < nnz; ++p) {
    if (!kutil::pos(wt[p])) return false;
  }
  return true;
}

inline double wof(const double *wt, int p) { return wt ? wt[p] : 1.0; }

// Breadth-first search from s over one or two CSR graphs (the second may be
// null): every node reached gets comp[v] = id; q is a queue of n ints.
void bfs(int s, int id, const int *off, const int *adj, const int *ioff,
         const int *iadj, int *comp, int *q) {
  int head = 0, tail = 0;
  comp[s] = id;
  q[tail++] = s;
  while (head < tail) {
    const int v = q[head++];
    for (int p = off[v]; p < off[v + 1]; ++p) {
      if (comp[adj[p]] < 0) {
        comp[adj[p]] = id;
        q[tail++] = adj[p];
      }
    }
    if (ioff == nullptr) continue;
    for (int p = ioff[v]; p < ioff[v + 1]; ++p) {
      if (comp[iadj[p]] >= 0) continue;
      comp[iadj[p]] = id;
      q[tail++] = iadj[p];
    }
  }
}

// Components numbered 0, 1, ... in the order of their first node; returns
// how many.
int comps(int n, const int *off, const int *adj, const int *ioff,
          const int *iadj, int *comp, int *q) {
  for (int i = 0; i < n; ++i) comp[i] = -1;
  int nc = 0;
  for (int s = 0; s < n; ++s) {
    if (comp[s] < 0) bfs(s, nc++, off, adj, ioff, iadj, comp, q);
  }
  return nc;
}

// A mix of the two endpoints of a pair, for the co-offence hash table.
inline std::uint64_t phash(int u, int v) {
  std::uint64_t z = (static_cast<std::uint64_t>(static_cast<unsigned>(u)) << 32) |
                    static_cast<unsigned>(v);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

// Rows of a co-offence table grouped by incident, stably: rows of incident
// c are ord[bo[c]] .. ord[bo[c + 1] - 1], in input order.
void bucket(const int *inc, int k, int ni, int *bo, int *ord) {
  for (int c = 0; c <= ni; ++c) bo[c] = 0;
  for (int r = 0; r < k; ++r) ++bo[inc[r] + 1];
  for (int c = 0; c < ni; ++c) bo[c + 1] += bo[c];
  for (int r = 0; r < k; ++r) ord[bo[inc[r]]++] = r;
  for (int c = ni; c > 0; --c) bo[c] = bo[c - 1];
  bo[0] = 0;
}

bool coofok(const int *inc, const int *off, int k, int ni, int no) {
  if (k < 0 || ni < 0 || no < 0 || ni == INT_MAX) return false;
  for (int r = 0; r < k; ++r) {
    if (inc[r] < 0 || inc[r] >= ni || off[r] < 0 || off[r] >= no) return false;
  }
  return true;
}

}  // namespace

extern "C" {

int clab_ntcsr(const int *a, const int *b, const double *w, int m, int n,
               int dir, int *off, int *adj, double *wt, int *iwork,
               int *nloop, int *ndup) {
  if (n < 1 || m < 0 || (!dir && m > INT_MAX / 2)) return CLAB_EARG;
  for (int k = 0; k < m; ++k) {
    if (a[k] < 0 || a[k] >= n || b[k] < 0 || b[k] >= n) return CLAB_EARG;
    if (w != nullptr && !kutil::pos(w[k])) return CLAB_EDOM;
  }
  // row lengths with duplicates, then their starts
  for (int i = 0; i <= n; ++i) off[i] = 0;
  int loops = 0;
  for (int k = 0; k < m; ++k) {
    if (a[k] == b[k]) {
      ++loops;
      continue;
    }
    ++off[a[k] + 1];
    if (!dir) ++off[b[k] + 1];
  }
  for (int i = 0; i < n; ++i) off[i + 1] += off[i];
  // each row in input order
  for (int i = 0; i < n; ++i) iwork[i] = off[i];
  for (int k = 0; k < m; ++k) {
    if (a[k] == b[k]) continue;
    const double wk = w ? w[k] : 1.0;
    int p = iwork[a[k]]++;
    adj[p] = b[k];
    wt[p] = wk;
    if (dir) continue;
    p = iwork[b[k]]++;
    adj[p] = a[k];
    wt[p] = wk;
  }
  // merge repeated neighbours into their first entry, compacting in place:
  // iwork[j] is where j last went, which is in row i when it is at least
  // the row's new start
  for (int i = 0; i < n; ++i) iwork[i] = -1;
  int q = 0, dup = 0;
  for (int i = 0; i < n; ++i) {
    const int s = off[i], e = off[i + 1];
    off[i] = q;
    for (int p = s; p < e; ++p) {
      const int j = adj[p];
      if (iwork[j] >= off[i]) {
        wt[iwork[j]] += wt[p];
        ++dup;
      } else {
        iwork[j] = q;
        adj[q] = j;
        wt[q] = wt[p];
        ++q;
      }
    }
  }
  off[n] = q;
  *nloop = loops;
  *ndup = dir ? dup : dup / 2;
  return CLAB_OK;
}

int clab_ntdeg(int n, const int *off, const double *wt, int *deg,
               double *str) {
  if (n < 1 || off[0] != 0) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (off[i + 1] < off[i]) return CLAB_EARG;
  }
  if (!wtok(wt, off[n])) return CLAB_EDOM;
  for (int i = 0; i < n; ++i) {
    deg[i] = off[i + 1] - off[i];
    kutil::ksum s;
    for (int p = off[i]; p < off[i + 1]; ++p) s.add(wof(wt, p));
    str[i] = s.val();
  }
  return CLAB_OK;
}

//' @srrstats {G3.0} The network kernels compare path lengths as whole
//'   numbers and path counts exactly (Brandes); the iterative ones stop on
//'   a stated tolerance (eigenvector, Fiedler vector, Louvain moves).
//' @noRd
int clab_ntbtw(int n, const int *off, const int *adj, int dir, int norm,
               double *bc, double *work, int *iwork) {
  if (!csrok(n, off, adj)) return CLAB_EARG;
  double *sig = work, *del = work + n;
  int *dist = iwork, *ord = iwork + n;
  for (int i = 0; i < n; ++i) {
    bc[i] = 0.0;
    dist[i] = -1;
    sig[i] = 0.0;
    del[i] = 0.0;
  }
  for (int s = 0; s < n; ++s) {
    // shortest-path counts by breadth-first search
    int head = 0, tail = 0;
    dist[s] = 0;
    sig[s] = 1.0;
    ord[tail++] = s;
    while (head < tail) {
      const int v = ord[head++];
      for (int p = off[v]; p < off[v + 1]; ++p) {
        const int w = adj[p];
        if (dist[w] < 0) {
          dist[w] = dist[v] + 1;
          ord[tail++] = w;
        }
        if (dist[w] == dist[v] + 1) sig[w] += sig[v];
      }
    }
    // dependencies back from the farthest nodes: v's successors w are
    // the neighbours one step farther, done before v
    for (int t = tail - 1; t >= 0; --t) {
      const int v = ord[t];
      for (int p = off[v]; p < off[v + 1]; ++p) {
        const int w = adj[p];
        if (dist[w] == dist[v] + 1) del[v] += sig[v] / sig[w] * (1.0 + del[w]);
      }
      if (v != s) bc[v] += del[v];
    }
    for (int t = 0; t < tail; ++t) {
      dist[ord[t]] = -1;
      sig[ord[t]] = 0.0;
      del[ord[t]] = 0.0;
    }
  }
  // an undirected pair was counted from both ends
  // (the normalised scores divide by the ordered pairs of other nodes,
  // (n - 1)(n - 2), or by the unordered ones when undirected)
  const double nn = static_cast<double>(n);
  const double f = dir ? 1.0 : 0.5;
  const double sc = !norm ? f : n < 3 ? kutil::nan() : 1.0 / ((nn - 1) * (nn - 2));
  for (int i = 0; i < n; ++i) bc[i] *= sc;
  return CLAB_OK;
}

int clab_ntcls(int n, const int *off, const int *adj, double *cls,
               double *ncls, double *hrm, int *reach, int *iwork) {
  if (!csrok(n, off, adj)) return CLAB_EARG;
  int *dist = iwork, *q = iwork + n;
  for (int i = 0; i < n; ++i) dist[i] = -1;
  for (int s = 0; s < n; ++s) {
    int head = 0, tail = 0;
    dist[s] = 0;
    q[tail++] = s;
    // distances are whole numbers: their sum is exact in a double; the
    // harmonic sum adds count / d once per distance d
    double sum = 0.0, h = 0.0;
    int lev = 1, cnt = 0;
    while (head < tail) {
      const int v = q[head++];
      if (dist[v] > lev) {
        h += cnt / static_cast<double>(lev);
        lev = dist[v];
        cnt = 0;
      }
      if (v != s) {
        sum += dist[v];
        ++cnt;
      }
      for (int p = off[v]; p < off[v + 1]; ++p) {
        if (dist[adj[p]] < 0) {
          dist[adj[p]] = dist[v] + 1;
          q[tail++] = adj[p];
        }
      }
    }
    h += cnt / static_cast<double>(lev);
    const int r = tail - 1;
    reach[s] = r;
    cls[s] = r > 0 ? 1.0 / sum : kutil::nan();
    ncls[s] = r > 0 ? r / sum : kutil::nan();
    hrm[s] = h;
    for (int t = 0; t < tail; ++t) dist[q[t]] = -1;
  }
  return CLAB_OK;
}

int clab_ntevc(int n, const int *off, const int *adj, const double *wt,
               double tol, int maxit, double *x, double *val, int *iter,
               int *conv, double *work, int *iwork) {
  if (!csrok(n, off, adj) || maxit < 1) return CLAB_EARG;
  if (!wtok(wt, off[n]) || !(tol > 0 && tol < 1)) return CLAB_EDOM;
  int *comp = iwork;
  const int nc = comps(n, off, adj, nullptr, nullptr, comp, iwork + n);
  double *y = work, *cm = work + n, *lam = work + 2 * n, *den = work + 3 * n;
  // the shift c, the largest weight, makes A + c I's leading eigenvalue
  // strictly dominant on a connected bipartite graph too
  double c = 0.0;
  for (int p = 0; p < off[n]; ++p) c = wof(wt, p) > c ? wof(wt, p) : c;
  if (c == 0.0) c = 1.0;
  for (int i = 0; i < n; ++i) x[i] = 1.0;
  int it = 0;
  bool ok = false;
  while (!ok && it < maxit) {
    if (kutil::stop(n + off[n])) return CLAB_EINTR;
    ++it;
    for (int k = 0; k < nc; ++k) cm[k] = 0.0;
    for (int i = 0; i < n; ++i) {
      double s = c * x[i];
      for (int p = off[i]; p < off[i + 1]; ++p) s += wof(wt, p) * x[adj[p]];
      y[i] = s;
      cm[comp[i]] = s > cm[comp[i]] ? s : cm[comp[i]];
    }
    // each component scaled to a maximum of 1, its own power iteration
    double d = 0.0;
    for (int i = 0; i < n; ++i) {
      const double v = y[i] / cm[comp[i]];
      const double e = std::fabs(v - x[i]);
      d = e > d ? e : d;
      x[i] = v;
    }
    ok = d < tol;
  }
  // each component's eigenvalue by its Rayleigh quotient
  for (int k = 0; k < nc; ++k) {
    lam[k] = 0.0;
    den[k] = 0.0;
  }
  for (int i = 0; i < n; ++i) {
    double s = 0.0;
    for (int p = off[i]; p < off[i + 1]; ++p) s += wof(wt, p) * x[adj[p]];
    lam[comp[i]] += x[i] * s;
    den[comp[i]] += x[i] * x[i];
  }
  double lmax = 0.0;
  for (int k = 0; k < nc; ++k) {
    lam[k] /= den[k];
    lmax = lam[k] > lmax ? lam[k] : lmax;
  }
  // the eigenvector of the largest eigenvalue is 0 on every component
  // whose own leading eigenvalue is smaller
  for (int i = 0; i < n; ++i) {
    if (lam[comp[i]] < lmax * (1.0 - 1e-10)) x[i] = 0.0;
  }
  *val = lmax;
  *iter = it;
  *conv = ok;
  return CLAB_OK;
}

int clab_ntcore(int n, const int *off, const int *adj, int *core,
                int *iwork) {
  if (!csrok(n, off, adj)) return CLAB_EARG;
  // Batagelj and Zaversnik: nodes kept sorted by current degree in vert,
  // bin[d] the first position of degree d
  int *bin = iwork, *pos = iwork + n, *vert = iwork + 2 * n;
  int md = 0;
  for (int v = 0; v < n; ++v) {
    core[v] = off[v + 1] - off[v];
    md = core[v] > md ? core[v] : md;
  }
  for (int d = 0; d <= md; ++d) bin[d] = 0;
  for (int v = 0; v < n; ++v) ++bin[core[v]];
  for (int d = 0, s = 0; d <= md; ++d) {
    const int k = bin[d];
    bin[d] = s;
    s += k;
  }
  for (int v = 0; v < n; ++v) {
    pos[v] = bin[core[v]]++;
    vert[pos[v]] = v;
  }
  for (int d = md; d > 0; --d) bin[d] = bin[d - 1];
  bin[0] = 0;
  for (int i = 0; i < n; ++i) {
    const int v = vert[i];
    for (int p = off[v]; p < off[v + 1]; ++p) {
      const int u = adj[p];
      if (core[u] <= core[v]) continue;
      const int du = core[u], pu = pos[u], pw = bin[du], w = vert[pw];
      if (u != w) {
        pos[u] = pw;
        vert[pu] = w;
        pos[w] = pu;
        vert[pw] = u;
      }
      ++bin[du];
      --core[u];
    }
  }
  return CLAB_OK;
}

int clab_ntwcc(int n, const int *off, const int *adj, const int *ioff,
               const int *iadj, int *comp, int *nc, int *iwork) {
  if (!csrok(n, off, adj)) return CLAB_EARG;
  if (ioff != nullptr && !csrok(n, ioff, iadj)) return CLAB_EARG;
  *nc = comps(n, off, adj, ioff, iadj, comp, iwork);
  return CLAB_OK;
}

int clab_ntscc(int n, const int *off, const int *adj, int *comp, int *nc,
               int *iwork) {
  if (!csrok(n, off, adj)) return CLAB_EARG;
  // Tarjan's algorithm with an explicit call stack (cs, ep: the node and
  // its next edge) so that a long path cannot overflow the C stack
  int *idx = iwork, *low = iwork + n, *stk = iwork + 2 * n,
      *on = iwork + 3 * n, *cs = iwork + 4 * n, *ep = iwork + 5 * n;
  for (int i = 0; i < n; ++i) {
    idx[i] = -1;
    on[i] = 0;
  }
  int cnt = 0, sp = 0, k = 0;
  for (int s = 0; s < n; ++s) {
    if (idx[s] >= 0) continue;
    int top = 0;
    idx[s] = low[s] = cnt++;
    stk[sp++] = s;
    on[s] = 1;
    cs[top] = s;
    ep[top++] = off[s];
    while (top > 0) {
      const int v = cs[top - 1];
      if (ep[top - 1] < off[v + 1]) {
        const int w = adj[ep[top - 1]++];
        if (idx[w] < 0) {
          idx[w] = low[w] = cnt++;
          stk[sp++] = w;
          on[w] = 1;
          cs[top] = w;
          ep[top++] = off[w];
        } else if (on[w] && idx[w] < low[v]) {
          low[v] = idx[w];
        }
        continue;
      }
      if (low[v] == idx[v]) {
        int w;
        do {
          w = stk[--sp];
          on[w] = 0;
          comp[w] = k;
        } while (w != v);
        ++k;
      }
      if (--top > 0 && low[v] < low[cs[top - 1]]) low[cs[top - 1]] = low[v];
    }
  }
  // renumbered in the order of each component's first node
  for (int c = 0; c < k; ++c) idx[c] = -1;
  int m = 0;
  for (int i = 0; i < n; ++i) {
    if (idx[comp[i]] < 0) idx[comp[i]] = m++;
    comp[i] = idx[comp[i]];
  }
  *nc = k;
  return CLAB_OK;
}

int clab_ntmod(int n, const int *off, const int *adj, const double *wt,
               const int *memb, double res, double *q, double *work) {
  if (!csrok(n, off, adj)) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (memb[i] < 0 || memb[i] >= n) return CLAB_EARG;
  }
  if (!wtok(wt, off[n]) || !kutil::nonneg(res)) return CLAB_EDOM;
  double *in = work, *tot = work + n;
  for (int c = 0; c < n; ++c) in[c] = tot[c] = 0.0;
  kutil::ksum m2;
  for (int i = 0; i < n; ++i) {
    for (int p = off[i]; p < off[i + 1]; ++p) {
      const double w = wof(wt, p);
      tot[memb[i]] += w;
      if (memb[adj[p]] == memb[i]) in[memb[i]] += w;
      m2.add(w);
    }
  }
  const double t = m2.val();
  if (t == 0.0) {
    *q = kutil::nan();
    return CLAB_OK;
  }
  kutil::ksum s;
  for (int c = 0; c < n; ++c) {
    if (tot[c] == 0.0) continue;
    const double f = tot[c] / t;
    s.add(in[c] / t);
    s.add(-res * f * f);
  }
  *q = s.val();
  return CLAB_OK;
}

int clab_ntlvnw(int n, int nnz, double *len, double *ilen) {
  if (n < 1 || nnz < 0) return CLAB_EARG;
  *len = 2.0 * nnz + 3.0 * n;
  *ilen = 2.0 * nnz + 9.0 * n + 3.0;
  return CLAB_OK;
}

int clab_ntlvn(int n, const int *off, const int *adj, const double *wt,
               double res, unsigned int seed, int useseed, int *memb,
               double *q, int *lev, double *work, int *iwork) {
  if (!csrok(n, off, adj)) return CLAB_EARG;
  if (!wtok(wt, off[n]) || !kutil::pos(res)) return CLAB_EDOM;
  const int nnz = off[n];
  // two graphs in turn: the level read and the aggregate written
  double *wa = work, *wb = work + nnz, *k = work + 2 * nnz,
         *tot = k + n, *nw = tot + n;
  int *oa = iwork, *ob = oa + n + 1, *aa = ob + n + 1, *ab = aa + nnz,
      *comm = ab + nnz, *ord = comm + n, *tch = ord + n, *mark = tch + n,
      *map = mark + n, *bord = map + n, *bo = bord + n;
  kutil::ksum m2s;
  for (int i = 0; i <= n; ++i) oa[i] = off[i];
  for (int p = 0; p < nnz; ++p) {
    aa[p] = adj[p];
    wa[p] = wof(wt, p);
    m2s.add(wa[p]);
  }
  const double m2 = m2s.val();
  for (int i = 0; i < n; ++i) memb[i] = i;
  *lev = 0;
  kutil::rng g(seed);
  int cn = n;
  while (m2 > 0.0) {
    for (int i = 0; i < cn; ++i) {
      double s = 0.0;
      for (int p = oa[i]; p < oa[i + 1]; ++p) s += wa[p];
      k[i] = s;
      tot[i] = s;
      comm[i] = i;
      ord[i] = i;
      mark[i] = 0;
    }
    // the visiting order: node order, or a permutation drawn from the seed
    if (useseed) {
      for (int i = cn - 1; i > 0; --i) {
        const int j = g.below(i + 1), t = ord[i];
        ord[i] = ord[j];
        ord[j] = t;
      }
    }
    bool any = false;
    int moves;
    do {
      if (kutil::stop(cn + oa[cn])) return CLAB_EINTR;
      moves = 0;
      for (int t = 0; t < cn; ++t) {
        const int i = ord[t], ci = comm[i];
        int nt = 0;
        for (int p = oa[i]; p < oa[i + 1]; ++p) {
          const int j = aa[p];
          if (j == i) continue;
          const int c = comm[j];
          if (mark[c] != i + 1) {
            mark[c] = i + 1;
            nw[c] = 0.0;
            tch[nt++] = c;
          }
          nw[c] += wa[p];
        }
        tot[ci] -= k[i];
        const double own = (mark[ci] == i + 1 ? nw[ci] : 0.0) -
                           res * tot[ci] * k[i] / m2;
        int best = ci;
        double bg = own;
        for (int u = 0; u < nt; ++u) {
          const int c = tch[u];
          if (c == ci) continue;
          const double gc = nw[c] - res * tot[c] * k[i] / m2;
          if (best == ci ? gc > own + 1e-12 * k[i] : gc > bg) {
            best = c;
            bg = gc;
          }
        }
        // marks of this node are spent: a later node i' never equals i
        tot[best] += k[i];
        if (best != ci) {
          comm[i] = best;
          ++moves;
        }
      }
      any = any || moves > 0;
      for (int i = 0; i < cn; ++i) mark[i] = 0;
    } while (moves > 0);
    if (!any) break;
    ++*lev;
    // communities renumbered by their first node
    for (int c = 0; c < cn; ++c) map[c] = -1;
    int nc = 0;
    for (int i = 0; i < cn; ++i) {
      if (map[comm[i]] < 0) map[comm[i]] = nc++;
      comm[i] = map[comm[i]];
    }
    for (int v = 0; v < n; ++v) memb[v] = comm[memb[v]];
    // the aggregate: one node per community, rows gathered over its members
    // in node order, neighbours in order of first meeting
    for (int c = 0; c <= nc; ++c) bo[c] = 0;
    for (int i = 0; i < cn; ++i) ++bo[comm[i] + 1];
    for (int c = 0; c < nc; ++c) bo[c + 1] += bo[c];
    for (int i = 0; i < cn; ++i) bord[bo[comm[i]]++] = i;
    for (int c = nc; c > 0; --c) bo[c] = bo[c - 1];
    bo[0] = 0;
    int pq = 0;
    ob[0] = 0;
    for (int c = 0; c < nc; ++c) {
      int nt = 0;
      for (int u = bo[c]; u < bo[c + 1]; ++u) {
        const int i = bord[u];
        for (int p = oa[i]; p < oa[i + 1]; ++p) {
          const int d = comm[aa[p]];
          if (mark[d] != c + 1) {
            mark[d] = c + 1;
            nw[d] = 0.0;
            tch[nt++] = d;
          }
          nw[d] += wa[p];
        }
      }
      for (int u = 0; u < nt; ++u) {
        ab[pq] = tch[u];
        wb[pq++] = nw[tch[u]];
      }
      ob[c + 1] = pq;
    }
    int *ti = oa;
    oa = ob;
    ob = ti;
    ti = aa;
    aa = ab;
    ab = ti;
    double *td = wa;
    wa = wb;
    wb = td;
    cn = nc;
  }
  return clab_ntmod(n, off, adj, wt, memb, res, q, work + 2 * nnz);
}

int clab_ntbrk(int n, const int *off, const int *adj, const int *ioff,
               const int *iadj, const int *grp, double *cnt, int *iwork) {
  if (!csrok(n, off, adj) || !csrok(n, ioff, iadj)) return CLAB_EARG;
  int *stamp = iwork;
  for (int i = 0; i < n; ++i) stamp[i] = -1;
  for (int i = 0; i < 5 * n; ++i) cnt[i] = 0.0;
  // every two-path i -> v -> j with i != j and no tie i -> j, by the
  // groups of the three (Gould and Fernandez 1989)
  for (int v = 0; v < n; ++v) {
    for (int p = ioff[v]; p < ioff[v + 1]; ++p) {
      const int i = iadj[p];
      for (int r = off[i]; r < off[i + 1]; ++r) stamp[adj[r]] = p;
      for (int r = off[v]; r < off[v + 1]; ++r) {
        const int j = adj[r];
        if (j == i || stamp[j] == p) continue;
        const int gi = grp[i], gv = grp[v], gj = grp[j];
        const int role = gi == gv && gv == gj ? 0    // coordinator
                         : gi == gj           ? 1    // itinerant
                         : gi == gv           ? 2    // representative
                         : gv == gj           ? 3    // gatekeeper
                                              : 4;   // liaison
        cnt[role * n + v] += 1.0;
      }
    }
  }
  return CLAB_OK;
}

int clab_ntchg(int n, const int *off, const int *adj, const double *wt,
               double tol, int maxit, double *out, int *set, int *iter,
               int *conv, double *work, int *iwork) {
  if (!csrok(n, off, adj) || maxit < 1) return CLAB_EARG;
  if (!wtok(wt, off[n]) || !(tol > 0 && tol < 1) || n < 2) return CLAB_EDOM;
  int one = 0;
  for (int i = 0; i < n; ++i) set[i] = -1;
  bfs(0, 0, off, adj, nullptr, nullptr, set, iwork);
  for (int i = 0; i < n; ++i) one += set[i] == 0;
  if (one < n) return CLAB_EDOM;
  double *d = work, *s = work + n, *x = work + 2 * n, *y = work + 3 * n;
  double ss = 0.0;
  for (int i = 0; i < n; ++i) {
    kutil::ksum a;
    for (int p = off[i]; p < off[i + 1]; ++p) a.add(wof(wt, p));
    d[i] = a.val();
    s[i] = std::sqrt(d[i]);
    ss += d[i];
  }
  const double sn = std::sqrt(ss);
  // x against u = s / |s|, the eigenvector of eigenvalue 0, then unit
  auto deflate = [&](double *v) {
    double dot = 0.0, nn = 0.0;
    for (int i = 0; i < n; ++i) dot += v[i] * s[i] / sn;
    for (int i = 0; i < n; ++i) {
      v[i] -= dot * s[i] / sn;
      nn += v[i] * v[i];
    }
    nn = std::sqrt(nn);
    for (int i = 0; i < n; ++i) v[i] /= nn;
  };
  kutil::rng g(0x5EEDu);
  for (int i = 0; i < n; ++i) x[i] = g.unif() - 0.5;
  deflate(x);
  // power iteration on M = (I + D^-1/2 A D^-1/2) / 2, whose spectrum is in
  // [0, 1]; on the complement of u its top eigenvalue is 1 - lambda2 / 2.
  // Stop when the residual |M x - mu x| falls below tol.
  int it = 0;
  bool ok = false;
  while (!ok && it < maxit) {
    if (kutil::stop(n + off[n])) return CLAB_EINTR;
    ++it;
    double mu = 0.0;
    for (int i = 0; i < n; ++i) {
      double a = 0.0;
      for (int p = off[i]; p < off[i + 1]; ++p) {
        a += wof(wt, p) * x[adj[p]] / s[adj[p]];
      }
      y[i] = 0.5 * (x[i] + a / s[i]);
      mu += x[i] * y[i];
    }
    double r = 0.0, dot = 0.0;
    for (int i = 0; i < n; ++i) dot += y[i] * s[i] / sn;
    for (int i = 0; i < n; ++i) {
      const double e = y[i] - dot * s[i] / sn - mu * x[i];
      r += e * e;
    }
    ok = std::sqrt(r) < tol;
    for (int i = 0; i < n; ++i) x[i] = y[i];
    deflate(x);
  }
  // lambda2 as the Rayleigh quotient of the normalised Laplacian, a sum of
  // squares, so that a small value keeps its relative accuracy
  kutil::ksum lq;
  for (int i = 0; i < n; ++i) {
    for (int p = off[i]; p < off[i + 1]; ++p) {
      const double e = x[i] / s[i] - x[adj[p]] / s[adj[p]];
      lq.add(0.5 * wof(wt, p) * e * e);
    }
  }
  const double lam = lq.val();
  // the sweep: prefixes of the nodes in order of x / s (ties by node)
  int *o = iwork;
  for (int i = 0; i < n; ++i) {
    o[i] = i;
    y[i] = x[i] / s[i];
    set[i] = 0;
  }
  kutil::hsort(o, n, [&](int a, int b) {
    return y[a] < y[b] || (y[a] == y[b] && a < b);
  });
  double vol = 0.0, cut = 0.0, best = kutil::nan();
  int bk = 0;
  for (int t = 0; t < n - 1; ++t) {
    const int v = o[t];
    double in = 0.0;
    for (int p = off[v]; p < off[v + 1]; ++p) {
      if (set[adj[p]]) in += wof(wt, p);
    }
    set[v] = 1;
    vol += d[v];
    cut += d[v] - 2.0 * in;
    const double m = vol < ss - vol ? vol : ss - vol;
    const double phi = cut / m;
    if (!(phi >= best)) {
      best = phi;
      bk = t + 1;
    }
  }
  for (int i = 0; i < n; ++i) set[i] = 0;
  for (int t = 0; t < bk; ++t) set[o[t]] = 1;
  out[0] = lam;
  out[1] = lam / 2.0;
  out[2] = std::sqrt(2.0 * lam);
  out[3] = best;
  *iter = it;
  *conv = ok;
  return CLAB_OK;
}

int clab_coofw(const int *inc, const int *off, int k, int ni, int no,
               double *npair, double *ilen, int *iwork) {
  if (!coofok(inc, off, k, ni, no)) return CLAB_EARG;
  int *bo = iwork, *ord = iwork + ni + 1, *stamp = ord + k;
  bucket(inc, k, ni, bo, ord);
  for (int j = 0; j < no; ++j) stamp[j] = -1;
  double np = 0.0;
  for (int c = 0; c < ni; ++c) {
    double m = 0.0;
    for (int u = bo[c]; u < bo[c + 1]; ++u) {
      const int j = off[ord[u]];
      if (stamp[j] == c) continue;
      stamp[j] = c;
      m += 1.0;
    }
    np += m * (m - 1.0) / 2.0;
  }
  // the pairs must fit clab_coofg's table (at most INT_MAX / 4 of them)
  if (np > INT_MAX / 4) return CLAB_EARG;
  // the hash table of clab_coofg: a power of two at least twice the pairs
  double ts = 2.0;
  while (ts < 2.0 * np) ts *= 2.0;
  *npair = np;
  *ilen = (ni + 1.0) + 2.0 * k + no + ts;
  return CLAB_OK;
}

int clab_coofg(const int *inc, const int *off, int k, int ni, int no,
               int cap, int *pa, int *pb, double *pw, int *np, int *iwork) {
  if (!coofok(inc, off, k, ni, no) || cap < 0 || cap > INT_MAX / 4) {
    return CLAB_EARG;
  }
  int *bo = iwork, *ord = iwork + ni + 1, *mem = ord + k, *stamp = mem + k,
      *tab = stamp + no;
  int ts = 2;
  while (ts < 2 * cap) ts *= 2;
  for (int u = 0; u < ts; ++u) tab[u] = -1;
  for (int j = 0; j < no; ++j) stamp[j] = -1;
  bucket(inc, k, ni, bo, ord);
  int npr = 0;
  for (int c = 0; c < ni; ++c) {
    // the incident's distinct offenders, in the order of the rows
    int m = 0;
    for (int u = bo[c]; u < bo[c + 1]; ++u) {
      const int j = off[ord[u]];
      if (stamp[j] == c) continue;
      stamp[j] = c;
      mem[m++] = j;
    }
    for (int x = 0; x < m; ++x) {
      for (int y = x + 1; y < m; ++y) {
        const int a = mem[x] < mem[y] ? mem[x] : mem[y];
        const int b = mem[x] < mem[y] ? mem[y] : mem[x];
        std::uint64_t h = phash(a, b) & static_cast<std::uint64_t>(ts - 1);
        while (tab[h] >= 0 && !(pa[tab[h]] == a && pb[tab[h]] == b)) {
          h = (h + 1) & static_cast<std::uint64_t>(ts - 1);
        }
        if (tab[h] >= 0) {
          pw[tab[h]] += 1.0;
          continue;
        }
        if (npr == cap) return CLAB_EARG;
        tab[h] = npr;
        pa[npr] = a;
        pb[npr] = b;
        pw[npr++] = 1.0;
      }
    }
  }
  *np = npr;
  return CLAB_OK;
}

}  // extern "C"
