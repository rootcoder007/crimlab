// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Self-test of the network kernels: closed forms on paths, stars, cycles,
// complete graphs and two cliques joined by one edge; published values on
// Zachary's karate club; and, on seeded random graphs, brute force from
// the adjacency matrix: all-pairs distances and path counts for
// betweenness and closeness, peeling for cores, transitive closure for
// components, every triple for brokerage and every node set for the
// conductance.
#include <cmath>

#include "crimlab.h"
#include "tchk.h"

namespace {

constexpr int NX = 40;    // nodes at most
constexpr int EX = 800;   // CSR entries at most

struct gr {
  int n, dir;
  int off[NX + 1], adj[EX], ioff[NX + 1], iadj[EX];
  double wt[EX], iwt[EX];
};

// the CSR of an edge list (and of its reverse, for a directed graph)
int mk(gr &g, int n, const int *a, const int *b, const double *w, int m,
       int dir) {
  int iw[NX], nl, nd;
  g.n = n;
  g.dir = dir;
  int rc = clab_ntcsr(a, b, w, m, n, dir, g.off, g.adj, g.wt, iw, &nl, &nd);
  if (dir && rc == CLAB_OK) {
    rc = clab_ntcsr(b, a, w, m, n, 1, g.ioff, g.iadj, g.iwt, iw, &nl, &nd);
  }
  return rc;
}

// edges of named families
int path(int n, int *a, int *b) {
  for (int i = 0; i + 1 < n; ++i) {
    a[i] = i;
    b[i] = i + 1;
  }
  return n - 1;
}
int star(int k, int *a, int *b) {
  for (int i = 0; i < k; ++i) {
    a[i] = 0;
    b[i] = i + 1;
  }
  return k;
}
int cyc(int n, int *a, int *b) {
  for (int i = 0; i < n; ++i) {
    a[i] = i;
    b[i] = (i + 1) % n;
  }
  return n;
}
int kn(int n, int *a, int *b, int base = 0) {
  int m = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      a[m] = base + i;
      b[m++] = base + j;
    }
  }
  return m;
}
// two copies of K_k, nodes 0 .. k - 1 and k .. 2k - 1, joined by 0 - k
int bar(int k, int *a, int *b) {
  int m = kn(k, a, b);
  m += kn(k, a + m, b + m, k);
  a[m] = 0;
  b[m++] = k;
  return m;
}

// Zachary (1977), doi:10.1086/jar.33.4.3629752: the 78 ties among 34
// members, numbered from 1 as in the paper
const int zk[78][2] = {
    {1, 2},   {1, 3},   {1, 4},   {1, 5},   {1, 6},   {1, 7},   {1, 8},
    {1, 9},   {1, 11},  {1, 12},  {1, 13},  {1, 14},  {1, 18},  {1, 20},
    {1, 22},  {1, 32},  {2, 3},   {2, 4},   {2, 8},   {2, 14},  {2, 18},
    {2, 20},  {2, 22},  {2, 31},  {3, 4},   {3, 8},   {3, 9},   {3, 10},
    {3, 14},  {3, 28},  {3, 29},  {3, 33},  {4, 8},   {4, 13},  {4, 14},
    {5, 7},   {5, 11},  {6, 7},   {6, 11},  {6, 17},  {7, 17},  {9, 31},
    {9, 33},  {9, 34},  {10, 34}, {14, 34}, {15, 33}, {15, 34}, {16, 33},
    {16, 34}, {19, 33}, {19, 34}, {20, 34}, {21, 33}, {21, 34}, {23, 33},
    {23, 34}, {24, 26}, {24, 28}, {24, 30}, {24, 33}, {24, 34}, {25, 26},
    {25, 28}, {25, 32}, {26, 32}, {27, 30}, {27, 34}, {28, 34}, {29, 32},
    {29, 34}, {30, 33}, {30, 34}, {31, 33}, {31, 34}, {32, 33}, {32, 34},
    {33, 34}};

// all-pairs hop distances (-1: unreachable) and path counts of a graph
struct ap {
  int d[NX][NX];
  double s[NX][NX];
};
void allp(const gr &g, ap &r) {
  const int n = g.n;
  bool e[NX][NX] = {};
  for (int i = 0; i < n; ++i) {
    for (int p = g.off[i]; p < g.off[i + 1]; ++p) e[i][g.adj[p]] = true;
  }
  // Floyd-Warshall, then counts layer by layer: sigma(s, t) sums
  // sigma(s, u) over u one step nearer with an edge u -> t
  const int big = 1 << 20;
  int d[NX][NX];
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) d[i][j] = i == j ? 0 : e[i][j] ? 1 : big;
  }
  for (int k = 0; k < n; ++k) {
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        if (d[i][k] + d[k][j] < d[i][j]) d[i][j] = d[i][k] + d[k][j];
      }
    }
  }
  for (int s = 0; s < n; ++s) {
    for (int t = 0; t < n; ++t) {
      r.d[s][t] = d[s][t] >= big ? -1 : d[s][t];
      r.s[s][t] = s == t ? 1.0 : 0.0;
    }
    for (int l = 1; l < n; ++l) {
      for (int t = 0; t < n; ++t) {
        if (r.d[s][t] != l) continue;
        for (int u = 0; u < n; ++u) {
          if (e[u][t] && r.d[s][u] == l - 1) r.s[s][t] += r.s[s][u];
        }
      }
    }
  }
}

void bfbtw(const gr &g, const ap &r, double *bc) {
  for (int v = 0; v < g.n; ++v) {
    double b = 0.0;
    for (int s = 0; s < g.n; ++s) {
      for (int t = 0; t < g.n; ++t) {
        if (s == t || s == v || t == v || r.d[s][t] < 0) continue;
        if (r.d[s][v] >= 0 && r.d[v][t] >= 0 &&
            r.d[s][v] + r.d[v][t] == r.d[s][t]) {
          b += r.s[s][v] * r.s[v][t] / r.s[s][t];
        }
      }
    }
    bc[v] = g.dir ? b : b / 2;
  }
}

double btw(const gr &g, double *bc, int norm = 0) {
  double w[2 * NX];
  int iw[2 * NX];
  return clab_ntbtw(g.n, g.off, g.adj, g.dir, norm, bc, w, iw);
}

struct cl {
  double c[NX], nc[NX], h[NX];
  int r[NX];
};
cl cls(const gr &g) {
  cl o{};
  int iw[2 * NX];
  clab_ntcls(g.n, g.off, g.adj, o.c, o.nc, o.h, o.r, iw);
  return o;
}

struct ev {
  double x[NX], val;
  int it, conv;
};
ev evc(const gr &g, const double *wt = nullptr, double tol = 1e-14) {
  ev o{};
  double w[4 * NX];
  int iw[2 * NX];
  clab_ntevc(g.n, g.off, g.adj, wt, tol, 100000, o.x, &o.val, &o.it, &o.conv,
             w, iw);
  return o;
}

void core(const gr &g, int *c) {
  int iw[3 * NX];
  clab_ntcore(g.n, g.off, g.adj, c, iw);
}

double mod(const gr &g, const int *m, double res = 1.0,
           const double *wt = nullptr) {
  double q, w[2 * NX];
  clab_ntmod(g.n, g.off, g.adj, wt, m, res, &q, w);
  return q;
}

struct lv {
  int m[NX], lev;
  double q;
};
lv lvn(const gr &g, unsigned seed, int use, const double *wt = nullptr,
       double res = 1.0) {
  lv o{};
  static double w[2 * EX + 3 * NX];
  static int iw[2 * EX + 9 * NX + 3];
  double len, ilen;
  clab_ntlvnw(g.n, g.off[g.n], &len, &ilen);
  clab_ntlvn(g.n, g.off, g.adj, wt, res, seed, use, o.m, &o.q, &o.lev, w, iw);
  return o;
}

struct ch {
  double o[4];
  int set[NX], it, conv, rc;
};
ch chg(const gr &g, const double *wt = nullptr) {
  ch c{};
  double w[4 * NX];
  int iw[NX];
  c.rc = clab_ntchg(g.n, g.off, g.adj, wt, 1e-13, 200000, c.o, c.set, &c.it,
                    &c.conv, w, iw);
  return c;
}

// the least conductance over every proper node set, and the conductance
// of one set
double phi(const gr &g, const double *wt, const int *in) {
  double vol = 0, tot = 0, cut = 0;
  for (int i = 0; i < g.n; ++i) {
    for (int p = g.off[i]; p < g.off[i + 1]; ++p) {
      const double w = wt ? wt[p] : 1.0;
      tot += w;
      if (in[i]) vol += w;
      if (in[i] && !in[g.adj[p]]) cut += w;
    }
  }
  return cut / (vol < tot - vol ? vol : tot - vol);
}
double hmin(const gr &g, const double *wt) {
  double h = INFINITY;
  int in[NX];
  for (long s = 1; s < (1L << g.n) - 1; ++s) {
    for (int i = 0; i < g.n; ++i) in[i] = (s >> i) & 1;
    const double f = phi(g, wt, in);
    if (f < h) h = f;
  }
  return h;
}

void csr(tk &t) {
  // a loop, a repeat in the other direction and an isolated node
  const int a[5] = {0, 1, 2, 2, 1}, b[5] = {1, 2, 2, 1, 0};
  const double w[5] = {1, 2, 9, 3, 0.5};
  int off[5], adj[10], iw[4], nl, nd;
  double wt[10];
  tk_ieq(t, "net1", clab_ntcsr(a, b, w, 5, 4, 0, off, adj, wt, iw, &nl, &nd),
         CLAB_OK);
  tk_ieq(t, "net2", nl, 1);
  tk_ieq(t, "net3", nd, 2);
  const int eo[5] = {0, 1, 3, 4, 4}, ea[4] = {1, 0, 2, 1};
  const double ew[4] = {1.5, 1.5, 5, 5};
  bool ok = true;
  for (int i = 0; i < 5; ++i) ok = ok && off[i] == eo[i];
  for (int p = 0; p < 4; ++p) ok = ok && adj[p] == ea[p] && wt[p] == ew[p];
  tk_true(t, "net4", ok);
  // directed: 1 -> 0 and 0 -> 1 stay apart; 2 -> 1 twice merges
  const int a2[4] = {0, 1, 2, 2}, b2[4] = {1, 0, 1, 1};
  tk_ieq(t, "net5",
         clab_ntcsr(a2, b2, nullptr, 4, 3, 1, off, adj, wt, iw, &nl, &nd),
         CLAB_OK);
  tk_true(t, "net6", nl == 0 && nd == 1 && off[3] == 3 && wt[2] == 2);
  // the guards
  const int bad[1] = {4};
  const double neg[1] = {-1};
  tk_ieq(t, "net7", clab_ntcsr(bad, b, nullptr, 1, 4, 0, off, adj, wt, iw,
                               &nl, &nd), CLAB_EARG);
  tk_ieq(t, "net8", clab_ntcsr(a, b, neg, 1, 4, 0, off, adj, wt, iw, &nl, &nd),
         CLAB_EDOM);
  tk_ieq(t, "net9", clab_ntcsr(a, b, nullptr, 1, 0, 0, off, adj, wt, iw, &nl,
                               &nd), CLAB_EARG);
  tk_ieq(t, "net10", clab_ntcsr(a, b, nullptr, -1, 4, 0, off, adj, wt, iw,
                                &nl, &nd), CLAB_EARG);
  tk_ieq(t, "net11", clab_ntcsr(a, b, nullptr, (1 << 30) + 1, 4, 0, off, adj,
                                wt, iw, &nl, &nd), CLAB_EARG);
}

void degree(tk &t) {
  gr g;
  const int a[3] = {0, 0, 1}, b[3] = {1, 2, 0};
  const double w[3] = {2, 3, 4};
  mk(g, 4, a, b, w, 3, 1);
  int d[4], di[4];
  double s[4], si[4];
  clab_ntdeg(4, g.off, g.wt, d, s);
  clab_ntdeg(4, g.ioff, g.iwt, di, si);
  tk_true(t, "net20", d[0] == 2 && d[1] == 1 && d[2] == 0 && d[3] == 0);
  tk_true(t, "net21", di[0] == 1 && di[1] == 1 && di[2] == 1 && di[3] == 0);
  tk_true(t, "net22", s[0] == 5 && s[1] == 4 && si[2] == 3 && si[0] == 4);
  // handshake: degrees sum to twice the edges (Lean: handshake)
  int za[78], zb[78];
  for (int k = 0; k < 78; ++k) {
    za[k] = zk[k][0] - 1;
    zb[k] = zk[k][1] - 1;
  }
  mk(g, 34, za, zb, nullptr, 78, 0);
  int zd[34];
  double zs[34];
  clab_ntdeg(34, g.off, nullptr, zd, zs);
  int sum = 0;
  for (int i = 0; i < 34; ++i) sum += zd[i];
  tk_ieq(t, "net23", sum, 156);
  // the published degrees of members 1, 34, 33 and 12
  tk_true(t, "net24", zd[0] == 16 && zd[33] == 17 && zd[32] == 12 &&
                          zd[11] == 1);
  int off2[2] = {0, -1};
  tk_ieq(t, "net25", clab_ntdeg(1, off2, nullptr, zd, zs), CLAB_EARG);
  tk_ieq(t, "net26", clab_ntdeg(0, off2, nullptr, zd, zs), CLAB_EARG);
  const double bw[2] = {1, NAN};
  int o3[3] = {0, 1, 2};
  tk_ieq(t, "net28", clab_ntdeg(2, o3, bw, zd, zs), CLAB_EDOM);
}

void between(tk &t) {
  int a[EX], b[EX];
  gr g;
  double bc[NX];
  // path: node i separates i (n - 1 - i) pairs
  mk(g, 7, a, b, nullptr, path(7, a, b), 0);
  btw(g, bc);
  bool ok = true;
  for (int i = 0; i < 7; ++i) ok = ok && bc[i] == i * (6.0 - i);
  tk_true(t, "net30", ok);
  // star: the centre separates every pair of its k leaves
  mk(g, 6, a, b, nullptr, star(5, a, b), 0);
  btw(g, bc);
  tk_true(t, "net31", bc[0] == 10 && bc[3] == 0);
  btw(g, bc, 1);
  tk_eq(t, "net32", bc[0], 1.0, 1e-15);
  // cycles: (n - 1)(n - 3) / 8 for odd n, (n - 2)^2 / 8 for even n
  mk(g, 9, a, b, nullptr, cyc(9, a, b), 0);
  btw(g, bc);
  tk_eq(t, "net33", bc[4], 8.0 * 6 / 8, 1e-14);
  mk(g, 10, a, b, nullptr, cyc(10, a, b), 0);
  btw(g, bc);
  tk_eq(t, "net34", bc[7], 64.0 / 8, 1e-14);
  // complete: nobody is between
  mk(g, 6, a, b, nullptr, kn(6, a, b), 0);
  btw(g, bc);
  tk_eq(t, "net35", bc[2], 0.0, 0);
  // the karate club: 231.0714 for member 1, as tabulated in the
  // literature, and every member against the brute force
  int m = 0;
  for (; m < 78; ++m) {
    a[m] = zk[m][0] - 1;
    b[m] = zk[m][1] - 1;
  }
  mk(g, 34, a, b, nullptr, 78, 0);
  btw(g, bc);
  tk_eq(t, "net36", bc[0], 3235.0 / 14, 1e-12);
  static ap kp;
  double kb[NX];
  allp(g, kp);
  bfbtw(g, kp, kb);
  tk_veq(t, "net37", bc, kb, 34, 1e-13);
  // a directed path: only forward pairs pass through
  mk(g, 4, a, b, nullptr, path(4, a, b), 1);
  btw(g, bc);
  tk_true(t, "net38", bc[1] == 2 && bc[2] == 2 && bc[0] == 0);
  btw(g, bc, 1);
  tk_eq(t, "net39", bc[1], 2.0 / 6, 1e-15);
  mk(g, 2, a, b, nullptr, path(2, a, b), 0);
  btw(g, bc, 1);
  tk_nan(t, "net40", bc[0]);
  int bad[2] = {1, 0};
  tk_ieq(t, "net41", clab_ntbtw(1, bad, a, 0, 0, bc, nullptr, nullptr),
         CLAB_EARG);
  int o2[2] = {0, 1}, a1[1] = {3};
  tk_ieq(t, "net42", clab_ntbtw(1, o2, a1, 0, 0, bc, nullptr, nullptr),
         CLAB_EARG);
}

void closeness(tk &t) {
  int a[EX], b[EX];
  gr g;
  mk(g, 8, a, b, nullptr, path(7, a, b), 0);  // node 7 is isolated
  cl c = cls(g);
  // node i of a path of 7: distances sum to i(i + 1)/2 + (6 - i)(7 - i)/2
  bool ok = true;
  for (int i = 0; i < 7; ++i) {
    const double s = i * (i + 1) / 2.0 + (6 - i) * (7 - i) / 2.0;
    ok = ok && c.r[i] == 6 && std::fabs(c.c[i] - 1 / s) < 1e-16 &&
         std::fabs(c.nc[i] - 6 / s) < 1e-15;
  }
  tk_true(t, "net50", ok);
  // harmonic, node 0: 1 + 1/2 + ... + 1/6
  tk_eq(t, "net51", c.h[0], 49.0 / 20, 1e-15);
  tk_true(t, "net52", c.r[7] == 0 && std::isnan(c.c[7]) &&
                          std::isnan(c.nc[7]) && c.h[7] == 0);
  mk(g, 6, a, b, nullptr, star(5, a, b), 0);
  c = cls(g);
  tk_true(t, "net53", c.c[0] == 0.2 && c.c[2] == 1.0 / 9 && c.h[2] == 3);
  // a directed path reaches forward only
  mk(g, 4, a, b, nullptr, path(4, a, b), 1);
  c = cls(g);
  tk_true(t, "net54", c.r[0] == 3 && c.r[3] == 0 && c.c[1] == 1.0 / 3 &&
                          c.h[2] == 1);
  int bad[2] = {1, 0};
  tk_ieq(t, "net55", clab_ntcls(1, bad, a, c.c, c.nc, c.h, c.r, nullptr),
         CLAB_EARG);
}

void eigen(tk &t) {
  int a[EX], b[EX];
  gr g;
  const double pi = 3.14159265358979323846;
  // path of 6: sin(pi (i + 1) / 7) scaled, eigenvalue 2 cos(pi / 7)
  mk(g, 6, a, b, nullptr, path(6, a, b), 0);
  ev e = evc(g);
  double want[NX];
  for (int i = 0; i < 6; ++i) want[i] = std::sin(pi * (i + 1) / 7) / std::sin(3 * pi / 7);
  tk_veq(t, "net60", e.x, want, 6, 1e-12);
  tk_eq(t, "net61", e.val, 2 * std::cos(pi / 7), 1e-13);
  tk_true(t, "net62", e.conv == 1 && e.it > 1);
  // a star is bipartite: the shift makes it converge; centre 1, leaves
  // 1 / sqrt(k), eigenvalue sqrt(k)
  mk(g, 5, a, b, nullptr, star(4, a, b), 0);
  e = evc(g);
  tk_true(t, "net63", e.conv == 1);
  tk_eq(t, "net64", e.x[3], 0.5, 1e-12);
  tk_eq(t, "net65", e.val, 2.0, 1e-13);
  // a cycle of even length (bipartite) is flat
  mk(g, 8, a, b, nullptr, cyc(8, a, b), 0);
  e = evc(g);
  tk_true(t, "net66", e.conv == 1 && e.x[5] == 1 && std::fabs(e.val - 2) < 1e-14);
  // disconnected: K4 (eigenvalue 3) next to a path of 3 (sqrt 2) and an
  // isolated node; only K4 keeps weight
  int m = kn(4, a, b);
  a[m] = 4;
  b[m++] = 5;
  a[m] = 5;
  b[m++] = 6;
  mk(g, 8, a, b, nullptr, m, 0);
  e = evc(g);
  tk_true(t, "net67", e.x[0] == 1 && e.x[3] == 1 && e.x[4] == 0 &&
                          e.x[5] == 0 && e.x[7] == 0);
  tk_eq(t, "net68", e.val, 3.0, 1e-13);
  // two equal components: each scaled to 1
  m = kn(3, a, b);
  m += kn(3, a + m, b + m, 3);
  mk(g, 6, a, b, nullptr, m, 0);
  e = evc(g);
  tk_true(t, "net69", e.x[0] == 1 && e.x[4] == 1);
  // no edges: every node 1, eigenvalue 0
  mk(g, 3, a, b, nullptr, 0, 0);
  e = evc(g);
  tk_true(t, "net70", e.x[1] == 1 && e.val == 0);
  // weights scale the matrix: the leading eigenvector is unchanged when
  // every weight is doubled
  mk(g, 6, a, b, nullptr, path(6, a, b), 0);
  double w2[EX];
  for (int p = 0; p < g.off[6]; ++p) w2[p] = 2;
  e = evc(g, w2);
  tk_veq(t, "net71", e.x, want, 6, 1e-12);
  tk_eq(t, "net72", e.val, 4 * std::cos(pi / 7), 1e-13);
  // maxit reached: conv 0
  double x[6], val, wk[24];
  int it, cv, iw[12];
  clab_ntevc(6, g.off, g.adj, nullptr, 1e-14, 2, x, &val, &it, &cv, wk, iw);
  tk_true(t, "net73", cv == 0 && it == 2);
  tk_ieq(t, "net74", clab_ntevc(6, g.off, g.adj, nullptr, 0, 2, x, &val, &it,
                                &cv, wk, iw), CLAB_EDOM);
  tk_ieq(t, "net75", clab_ntevc(6, g.off, g.adj, nullptr, 1e-9, 0, x, &val,
                                &it, &cv, wk, iw), CLAB_EARG);
  // the karate club: residual of A x = val x on the converged vector
  for (m = 0; m < 78; ++m) {
    a[m] = zk[m][0] - 1;
    b[m] = zk[m][1] - 1;
  }
  mk(g, 34, a, b, nullptr, 78, 0);
  e = evc(g);
  double r = 0;
  for (int i = 0; i < 34; ++i) {
    double s = 0;
    for (int p = g.off[i]; p < g.off[i + 1]; ++p) s += e.x[g.adj[p]];
    r = std::fmax(r, std::fabs(s - e.val * e.x[i]));
  }
  tk_lt(t, "net76", r, 1e-11);
  tk_true(t, "net77", e.x[33] == 1);
}

void cores(tk &t) {
  int a[EX], b[EX], c[NX];
  gr g;
  // K4 with a path of two hanging off node 0, and an isolated node
  int m = kn(4, a, b);
  a[m] = 0;
  b[m++] = 4;
  a[m] = 4;
  b[m++] = 5;
  mk(g, 7, a, b, nullptr, m, 0);
  core(g, c);
  tk_true(t, "net80", c[0] == 3 && c[3] == 3 && c[4] == 1 && c[5] == 1 &&
                          c[6] == 0);
  mk(g, 9, a, b, nullptr, cyc(9, a, b), 0);
  core(g, c);
  tk_true(t, "net81", c[0] == 2 && c[8] == 2);
  int bad[2] = {1, 0};
  tk_ieq(t, "net82", clab_ntcore(1, bad, a, c, nullptr), CLAB_EARG);
}

void comp(tk &t) {
  int a[EX], b[EX], c[NX], nc, iw[6 * NX];
  gr g;
  // 0 -> 1 -> 2 -> 0, 2 -> 3, 4 alone, 5 -> 3
  const int da[5] = {0, 1, 2, 2, 5}, db[5] = {1, 2, 0, 3, 3};
  mk(g, 6, da, db, nullptr, 5, 1);
  clab_ntwcc(6, g.off, g.adj, g.ioff, g.iadj, c, &nc, iw);
  tk_true(t, "net90", nc == 2 && c[0] == 0 && c[3] == 0 && c[5] == 0 &&
                          c[4] == 1);
  // following out-edges only, 3 cannot reach 5's side
  clab_ntwcc(6, g.off, g.adj, nullptr, nullptr, c, &nc, iw);
  tk_true(t, "net91", nc == 3 && c[5] == 2 && c[4] == 1);
  clab_ntscc(6, g.off, g.adj, c, &nc, iw);
  tk_true(t, "net92", nc == 4 && c[0] == 0 && c[1] == 0 && c[2] == 0 &&
                          c[3] == 1 && c[4] == 2 && c[5] == 3);
  // a long directed cycle needs no recursion
  mk(g, NX, a, b, nullptr, cyc(NX, a, b), 1);
  clab_ntscc(NX, g.off, g.adj, c, &nc, iw);
  tk_true(t, "net93", nc == 1 && c[NX - 1] == 0);
  int bad[2] = {1, 0};
  tk_ieq(t, "net94", clab_ntscc(1, bad, a, c, &nc, iw), CLAB_EARG);
  tk_ieq(t, "net95", clab_ntwcc(1, bad, a, nullptr, nullptr, c, &nc, iw),
         CLAB_EARG);
  tk_ieq(t, "net96", clab_ntwcc(6, g.off, g.adj, bad, a, c, &nc, iw),
         CLAB_EARG);
}

void modular(tk &t) {
  int a[EX], b[EX];
  gr g;
  // two K5 joined by one edge: Q of the two cliques is 20/21 - 1/2
  mk(g, 10, a, b, nullptr, bar(5, a, b), 0);
  int two[10], one[10], sing[10];
  for (int i = 0; i < 10; ++i) {
    two[i] = i / 5;
    one[i] = 0;
    sing[i] = i;
  }
  tk_eq(t, "net100", mod(g, two), 20.0 / 21 - 0.5, 1e-15);
  tk_eq(t, "net101", mod(g, one), 0.0, 1e-15);
  // singletons: minus the sum of squared degree shares
  tk_eq(t, "net102", mod(g, sing), -(8 * 16.0 + 2 * 25.0) / (42.0 * 42),
        1e-15);
  // resolution 0 counts only the weight inside
  tk_eq(t, "net103", mod(g, two, 0.0), 40.0 / 42, 1e-15);
  // Louvain finds the two cliques, in node order and under any seed
  lv l = lvn(g, 0, 0);
  bool ok = true;
  for (int i = 0; i < 10; ++i) ok = ok && l.m[i] == two[i];
  tk_true(t, "net104", ok);
  tk_eq(t, "net105", l.q, 20.0 / 21 - 0.5, 1e-15);
  for (unsigned s = 1; s < 6; ++s) {
    l = lvn(g, s, 1);
    ok = true;
    for (int i = 0; i < 10; ++i) ok = ok && l.m[i] == two[i];
    tk_true(t, "net106", ok && l.lev >= 1);
  }
  // no edges: everyone alone, Q NaN
  mk(g, 3, a, b, nullptr, 0, 0);
  l = lvn(g, 0, 0);
  tk_true(t, "net107", l.m[2] == 2 && l.lev == 0 && std::isnan(l.q));
  // the karate club: the best partition has Q = 0.4197896 (Brandes et al.
  // 2008, doi:10.1109/TKDE.2007.190689); Louvain is a heuristic, and some
  // visiting orders stop at a local optimum near 0.39
  int m = 0;
  for (; m < 78; ++m) {
    a[m] = zk[m][0] - 1;
    b[m] = zk[m][1] - 1;
  }
  mk(g, 34, a, b, nullptr, 78, 0);
  for (unsigned s = 0; s < 8; ++s) {
    l = lvn(g, s, s > 0);
    tk_le(t, "net108", l.q, 0.4197896 + 1e-7);
    tk_ge(t, "net109", l.q, 0.38);
    tk_eq(t, "net110", l.q, mod(g, l.m), 1e-15);
    // the same seed, the same partition
    lv l2 = lvn(g, s, s > 0);
    ok = true;
    for (int i = 0; i < 34; ++i) ok = ok && l.m[i] == l2.m[i];
    tk_true(t, "net111", ok);
    // numbered by first node
    int hi = -1;
    for (int i = 0; i < 34; ++i) {
      ok = ok && l.m[i] <= hi + 1;
      hi = l.m[i] > hi ? l.m[i] : hi;
    }
    tk_true(t, "net112", ok);
  }
  // some seed reaches the optimum
  double qb = 0;
  for (unsigned s = 1; s < 40; ++s) qb = std::fmax(qb, lvn(g, s, 1).q);
  tk_eq(t, "net121", qb, 0.4197896, 1e-6);
  // weights: a heavy edge keeps its pair together
  mk(g, 4, a, b, nullptr, path(4, a, b), 0);
  double w[EX];
  for (int p = 0; p < g.off[4]; ++p) w[p] = 1;
  w[2] = w[3] = 10;  // the edge 1 - 2, held in rows 1 and 2
  l = lvn(g, 0, 0, w);
  tk_true(t, "net113", l.m[1] == l.m[2]);
  int bm[4] = {0, 0, 4, 0};
  double q, wk[8];
  tk_ieq(t, "net114", clab_ntmod(4, g.off, g.adj, nullptr, bm, 1, &q, wk),
         CLAB_EARG);
  bm[2] = 0;
  tk_ieq(t, "net115", clab_ntmod(4, g.off, g.adj, nullptr, bm, -1, &q, wk),
         CLAB_EDOM);
  int bad[2] = {1, 0};
  tk_ieq(t, "net116", clab_ntmod(1, bad, a, nullptr, bm, 1, &q, wk),
         CLAB_EARG);
  int mm[4], lev;
  tk_ieq(t, "net117", clab_ntlvn(4, g.off, g.adj, nullptr, 0, 0, 0, mm, &q,
                                 &lev, nullptr, nullptr), CLAB_EDOM);
  tk_ieq(t, "net118", clab_ntlvn(1, bad, a, nullptr, 1, 0, 0, mm, &q, &lev,
                                 nullptr, nullptr), CLAB_EARG);
  double len, il;
  tk_ieq(t, "net119", clab_ntlvnw(0, 0, &len, &il), CLAB_EARG);
  clab_ntlvnw(4, 6, &len, &il);
  tk_true(t, "net120", len == 24 && il == 51);
}

void broker(tk &t) {
  // two-paths through node 0 of every kind, groups A (0, 1, 2), B (3, 4),
  // C (5) and D (6); the counts by brute force over every ordered triple
  const int a[9] = {1, 0, 3, 0, 5, 0, 1, 6, 0}, b[9] = {0, 2, 0, 4, 0, 6, 2, 0, 5};
  const int grp[7] = {0, 0, 0, 1, 1, 2, 3};
  gr g;
  mk(g, 7, a, b, nullptr, 9, 1);
  double cnt[35];
  int iw[7];
  clab_ntbrk(7, g.off, g.adj, g.ioff, g.iadj, grp, cnt, iw);
  // brute force over every ordered triple
  bool e[7][7] = {};
  for (int k = 0; k < 9; ++k) e[a[k]][b[k]] = true;
  double bf[35] = {};
  for (int v = 0; v < 7; ++v) {
    for (int i = 0; i < 7; ++i) {
      for (int j = 0; j < 7; ++j) {
        if (i == j || i == v || j == v || !e[i][v] || !e[v][j] || e[i][j]) continue;
        const int gi = grp[i], gv = grp[v], gj = grp[j];
        int r = 4;
        if (gi == gv && gv == gj) r = 0;
        else if (gi == gj) r = 1;
        else if (gi == gv) r = 2;
        else if (gv == gj) r = 3;
        bf[r * 7 + v] += 1;
      }
    }
  }
  tk_veq(t, "net130", cnt, bf, 35, 0);
  // 1 -> 0 -> 2 is closed by 1 -> 2, so 0 coordinates nothing; 1 -> 0 -> 4
  // makes 0 a representative of A
  tk_true(t, "net131", cnt[0] == 0 && cnt[14] >= 1);
  int bad[2] = {1, 0};
  tk_ieq(t, "net132", clab_ntbrk(1, bad, a, g.ioff, g.iadj, grp, cnt, iw),
         CLAB_EARG);
  tk_ieq(t, "net133", clab_ntbrk(7, g.off, g.adj, bad, a, grp, cnt, iw),
         CLAB_EARG);
}

void cheeger(tk &t) {
  int a[EX], b[EX];
  gr g;
  const double pi = 3.14159265358979323846;
  // closed forms of lambda2: path 1 - cos(pi / (n - 1)), cycle
  // 1 - cos(2 pi / n), complete n / (n - 1), star 1
  mk(g, 8, a, b, nullptr, path(8, a, b), 0);
  ch c = chg(g);
  tk_eq(t, "net140", c.o[0], 1 - std::cos(pi / 7), 1e-11);
  tk_true(t, "net141", c.conv == 1 && c.rc == CLAB_OK);
  mk(g, 9, a, b, nullptr, cyc(9, a, b), 0);
  c = chg(g);
  tk_eq(t, "net142", c.o[0], 1 - std::cos(2 * pi / 9), 1e-11);
  mk(g, 6, a, b, nullptr, kn(6, a, b), 0);
  c = chg(g);
  tk_eq(t, "net143", c.o[0], 6.0 / 5, 1e-11);
  mk(g, 6, a, b, nullptr, star(5, a, b), 0);
  c = chg(g);
  tk_eq(t, "net144", c.o[0], 1.0, 1e-11);
  // two cliques: the sweep finds one clique, the cut of one edge
  mk(g, 10, a, b, nullptr, bar(5, a, b), 0);
  c = chg(g);
  tk_eq(t, "net145", c.o[3], 1.0 / 21, 1e-15);
  tk_true(t, "net146", c.set[0] == c.set[4] && c.set[5] == c.set[9] &&
                           c.set[0] != c.set[5]);
  // seeded random connected graphs: lambda2 / 2 <= h <= sweep <= sqrt(2
  // lambda2), h by every node set
  rng r = {77};
  for (int rep = 0; rep < 40; ++rep) {
    const int n = rsam(r, 3, 11);
    int m = path(n, a, b);  // a spanning path keeps it connected
    for (int i = 0; i < n; ++i) {
      for (int j = i + 2; j < n; ++j) {
        if (runi(r) < 0.3) {
          a[m] = i;
          b[m++] = j;
        }
      }
    }
    double w[EX];
    for (int k = 0; k < m; ++k) w[k] = rep % 2 ? runi(r, 0.5, 3) : 1;
    mk(g, n, a, b, w, m, 0);
    c = chg(g, g.wt);
    const double h = hmin(g, g.wt);
    tk_le(t, "net147", c.o[1], h + 1e-12);
    tk_le(t, "net148", h, c.o[3] + 1e-12);
    tk_le(t, "net149", c.o[3], c.o[2] + 1e-9);
    tk_eq(t, "net150", phi(g, g.wt, c.set), c.o[3], 1e-12);
  }
  // not connected, too small, bad arguments
  mk(g, 4, a, b, nullptr, path(3, a, b), 0);
  tk_ieq(t, "net151", chg(g).rc, CLAB_EDOM);
  mk(g, 1, a, b, nullptr, 0, 0);
  tk_ieq(t, "net152", chg(g).rc, CLAB_EDOM);
  int bad[2] = {1, 0};
  double o[4], wk[8];
  int set[2], it, cv, iw[2];
  tk_ieq(t, "net153", clab_ntchg(1, bad, a, nullptr, 1e-9, 9, o, set, &it, &cv,
                                 wk, iw), CLAB_EARG);
  // maxit reached
  mk(g, 8, a, b, nullptr, path(8, a, b), 0);
  double wk8[32];
  int s8[8], i8[8];
  clab_ntchg(8, g.off, g.adj, nullptr, 1e-13, 3, o, s8, &it, &cv, wk8, i8);
  tk_true(t, "net154", cv == 0 && it == 3);
}

void cooff(tk &t) {
  // incidents 0: offenders 0 1 2; 1: 1 0 (again) 1; 2: 3 alone; 3: 2 1
  const int inc[9] = {0, 0, 0, 1, 1, 1, 2, 3, 3},
            off[9] = {0, 1, 2, 1, 0, 1, 3, 2, 1};
  double np, il;
  int iw[64];
  tk_ieq(t, "net160", clab_coofw(inc, off, 9, 4, 4, &np, &il, iw), CLAB_OK);
  tk_eq(t, "net161", np, 3 + 1 + 0 + 1, 0);
  tk_eq(t, "net162", il, 5 + 18 + 4 + 16, 0);
  int pa[5], pb[5], n;
  double pw[5];
  clab_coofg(inc, off, 9, 4, 4, 5, pa, pb, pw, &n, iw);
  tk_ieq(t, "net163", n, 3);
  tk_true(t, "net164", pa[0] == 0 && pb[0] == 1 && pw[0] == 2 && pa[1] == 0 &&
                           pb[1] == 2 && pw[1] == 1 && pa[2] == 1 &&
                           pb[2] == 2 && pw[2] == 2);
  // too small a capacity and bad codes
  tk_ieq(t, "net165", clab_coofg(inc, off, 9, 4, 4, 2, pa, pb, pw, &n, iw),
         CLAB_EARG);
  tk_ieq(t, "net166", clab_coofw(inc, off, 9, 3, 4, &np, &il, iw), CLAB_EARG);
  tk_ieq(t, "net167", clab_coofg(inc, off, 9, 4, 3, 5, pa, pb, pw, &n, iw),
         CLAB_EARG);
  tk_ieq(t, "net168", clab_coofw(inc, off, -1, 4, 4, &np, &il, iw), CLAB_EARG);
  // one incident of 32769 offenders has more pairs than the table holds
  const int big = 32769;
  static int bi[big], bo[big], bw[2 + 2 * big];
  for (int r = 0; r < big; ++r) {
    bi[r] = 0;
    bo[r] = r;
  }
  tk_ieq(t, "net169", clab_coofw(bi, bo, big, 1, big, &np, &il, bw), CLAB_EARG);
}

// seeded random graphs against the brute-force references
void rgraph(tk &t, int big) {
  rng r = {2024};
  int a[EX], b[EX];
  gr g;
  const int reps = big ? 400 : 60;
  for (int rep = 0; rep < reps; ++rep) {
    const int n = rsam(r, 1, big ? 14 : 10), dir = rep % 3 == 0;
    const double pe = runi(r, 0.05, 0.6);
    int m = 0;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        if ((dir ? i != j : i < j) && runi(r) < pe) {
          a[m] = i;
          b[m++] = j;
        }
      }
    }
    mk(g, n, a, b, nullptr, m, dir);
    ap p;
    allp(g, p);
    double bc[NX], bf[NX];
    btw(g, bc);
    bfbtw(g, p, bf);
    tk_veq(t, "net170", bc, bf, n, 1e-13);
    cl c = cls(g);
    bool ok = true;
    for (int s = 0; s < n; ++s) {
      double sum = 0, h = 0;
      int rc = 0;
      for (int u = 0; u < n; ++u) {
        if (u == s || p.d[s][u] < 0) continue;
        ++rc;
        sum += p.d[s][u];
        h += 1.0 / p.d[s][u];
      }
      ok = ok && c.r[s] == rc && std::fabs(c.h[s] - h) < 1e-13 &&
           (rc == 0 ? std::isnan(c.c[s]) : c.c[s] == 1 / sum);
    }
    tk_true(t, "net171", ok);
    // components by transitive closure, numbered by first node
    int sc[NX], wc[NX], nc, iw[6 * NX];
    clab_ntscc(n, g.off, g.adj, sc, &nc, iw);
    clab_ntwcc(n, g.off, g.adj, dir ? g.ioff : nullptr,
               dir ? g.iadj : nullptr, wc, &nc, iw);
    int es[NX], ew[NX], ns = 0, nw = 0;
    for (int i = 0; i < n; ++i) es[i] = ew[i] = -1;
    for (int i = 0; i < n; ++i) {
      if (es[i] < 0) {
        for (int j = i; j < n; ++j) {
          if (es[j] < 0 && p.d[i][j] >= 0 && p.d[j][i] >= 0) es[j] = ns;
        }
        ++ns;
      }
    }
    // weak: reach in the symmetric closure
    gr u;
    int ua[EX], ub[EX];
    for (int k = 0; k < m; ++k) {
      ua[k] = a[k];
      ub[k] = b[k];
    }
    mk(u, n, ua, ub, nullptr, m, 0);
    ap pu;
    allp(u, pu);
    for (int i = 0; i < n; ++i) {
      if (ew[i] < 0) {
        for (int j = i; j < n; ++j) {
          if (pu.d[i][j] >= 0) ew[j] = nw;
        }
        ++nw;
      }
    }
    ok = nc == nw;
    for (int i = 0; i < n; ++i) ok = ok && sc[i] == es[i] && wc[i] == ew[i];
    tk_true(t, "net172", ok);
    // cores by peeling the undirected graph
    int cr[NX];
    core(u, cr);
    int bfc[NX], alive[NX];
    for (int i = 0; i < n; ++i) bfc[i] = 0;
    for (int k = 1; k < n; ++k) {
      for (int i = 0; i < n; ++i) alive[i] = 1;
      bool ch = true;
      while (ch) {
        ch = false;
        for (int i = 0; i < n; ++i) {
          if (!alive[i]) continue;
          int d = 0;
          for (int q = u.off[i]; q < u.off[i + 1]; ++q) d += alive[u.adj[q]];
          if (d < k) {
            alive[i] = 0;
            ch = true;
          }
        }
      }
      for (int i = 0; i < n; ++i) {
        if (alive[i]) bfc[i] = k;
      }
    }
    ok = true;
    for (int i = 0; i < n; ++i) ok = ok && cr[i] == bfc[i];
    tk_true(t, "net173", ok);
    // Louvain: its Q is the modularity of its partition, never below the
    // singletons' (moves only raise Q)
    if (u.off[n] > 0) {
      lv l = lvn(u, rep, rep % 2);
      int sg[NX];
      for (int i = 0; i < n; ++i) sg[i] = i;
      tk_eq(t, "net174", l.q, mod(u, l.m), 1e-14);
      tk_ge(t, "net175", l.q + 1e-14, mod(u, sg));
    }
  }
}

}  // namespace

void t_net(tk &t, int big) {
  csr(t);
  degree(t);
  between(t);
  closeness(t);
  eigen(t);
  cores(t);
  comp(t);
  modular(t);
  broker(t);
  cheeger(t);
  cooff(t);
  rgraph(t, big);
}
