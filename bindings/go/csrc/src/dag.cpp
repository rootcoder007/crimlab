// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Colliders of a directed acyclic graph (crimlab.h, dagcol()): every
// triple a -> c <- b of two causes of one node, shielded when a and b
// are adjacent (an edge either way), else an unshielded collider, a
// v-structure (Pearl 2009; Verma and Pearl 1991). A graph with a cycle is
// refused (Kahn's algorithm). Work: O(e log e + the sum over the nodes of
// their in-degree squared).
#include <algorithm>
#include <cstddef>
#include <set>
#include <utility>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

extern "C" int clab_dagcol(const int *from, const int *to, int ne, int nn, int cap, int *ca, int *cc, int *cb,
                           int *sh, int *nc, int *nd) CLAB_TRY {
  using std::size_t;
  if (ne < 1 || nn < 1 || cap < 0) return CLAB_EARG;
  std::set<std::pair<int, int>> e;
  for (int i = 0; i < ne; ++i) {
    if (from[i] < 0 || from[i] >= nn || to[i] < 0 || to[i] >= nn) return CLAB_EARG;
    if (from[i] == to[i]) return CLAB_EDOM;
    e.emplace(from[i], to[i]);
  }
  const size_t N = static_cast<size_t>(nn);
  std::vector<std::vector<int>> par(N), ch(N);
  std::vector<int> indeg(N, 0);
  // the set is ordered by cause, so each node's causes come in order
  for (const auto &p : e) {
    par[static_cast<size_t>(p.second)].push_back(p.first);
    ch[static_cast<size_t>(p.first)].push_back(p.second);
    ++indeg[static_cast<size_t>(p.second)];
  }
  std::vector<int> q;
  for (int v = 0; v < nn; ++v) {
    if (indeg[static_cast<size_t>(v)] == 0) q.push_back(v);
  }
  for (size_t h = 0; h < q.size(); ++h) {
    for (int w : ch[static_cast<size_t>(q[h])]) {
      if (--indeg[static_cast<size_t>(w)] == 0) q.push_back(w);
    }
  }
  if (q.size() != N) return CLAB_EDOM;
  double tot = 0.0;
  for (const auto &ps : par) tot += 0.5 * static_cast<double>(ps.size()) * (static_cast<double>(ps.size()) - 1);
  if (tot > kutil::maxcell) return CLAB_EARG;
  int k = 0;
  for (int c = 0; c < nn; ++c) {
    const std::vector<int> &ps = par[static_cast<size_t>(c)];
    for (size_t i = 0; i < ps.size(); ++i) {
      for (size_t j = i + 1; j < ps.size(); ++j) {
        if (k < cap) {
          ca[k] = ps[i];
          cc[k] = c;
          cb[k] = ps[j];
          sh[k] = e.count({ps[i], ps[j]}) + e.count({ps[j], ps[i]}) > 0 ? 1 : 0;
        }
        ++k;
      }
    }
  }
  *nc = k;
  *nd = static_cast<int>(e.size());
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)
