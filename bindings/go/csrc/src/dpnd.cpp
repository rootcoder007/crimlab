// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The work of gmmpnl(), wsctst() and ccepnl() between arguments and
// results (dpnd.h): the design of difference and system GMM as Roodman
// (2009) builds it, the calls of the kernels of dpn.cpp, and the results.
// Every array whose size is a product of the panel's dimensions is
// counted and refused past a bound before it is made.
#include "dpnd.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdarg>
#include <cstddef>
#include <cstdio>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "crimlab.h"
#include "kutil.h"
#include "gla.h"
#include "pev.h"

namespace dpnd {

namespace {

using std::size_t;
inline size_t z(int i) { return static_cast<size_t>(i); }

// the largest array made, in elements (2^26), and the most instruments
// or regressors (their weight matrices are square in that number)
constexpr double cap = 67108864.0;
constexpr int mmax = 2000;

std::string fmt(const char *f, ...) {
  std::va_list ap, aq;
  va_start(ap, f);
  va_copy(aq, ap);
  std::vector<char> b(256);
  const int n = std::vsnprintf(b.data(), b.size(), f, ap);
  if (n >= static_cast<int>(b.size())) {
    b.resize(static_cast<size_t>(n) + 1);
    std::vsnprintf(b.data(), b.size(), f, aq);
  }
  va_end(aq);
  va_end(ap);
  return b.data();
}

// The kernels' design of the stacked equations and what the result
// reports of it.
struct dsg {
  int n = 0, k = 0, m = 0, ng = 0, nd = 0, nobs = 0, omin = 0, omax = 0;
  std::vector<int> go, eq, dgo, dt, sub, grp;
  std::vector<double> y, x, z, h, yd, xd, hc, hdd;
  std::vector<int> lrow;          // the row of y of each level row
  std::vector<double> ly, lx;     // the level rows: y and the regressors kept
  std::vector<std::string> xn, zn, drop, dum;
  std::string bad;
  double big = 0;
};

std::string pv(double v) {
  char b[32];
  std::snprintf(b, sizeof b, "%.15g", v);
  return b;
}

// The rows of a panel by unit, in increasing period: unit a's rows are
// ord[uo[a]] .. ord[uo[a + 1] - 1].
struct idx {
  const pnl &pn;
  std::vector<int> ord, uo;
  explicit idx(const pnl &p) : pn(p), ord(z(p.n)), uo(z(p.nu) + 1, 0) {
    std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(), [&](int a, int b) {
      return p.u[z(a)] != p.u[z(b)] ? p.u[z(a)] < p.u[z(b)] : p.t[z(a)] < p.t[z(b)];
    });
    for (int i = 0; i < p.n; ++i) ++uo[z(p.u[z(i)]) + 1];
    std::partial_sum(uo.begin(), uo.end(), uo.begin());
  }
  // the row of unit a in period s, -1 when there is none
  int at(int a, int s) const {
    if (s < 0) return -1;
    const auto b = ord.begin() + uo[z(a)], e = ord.begin() + uo[z(a) + 1];
    const auto it = std::lower_bound(b, e, s, [&](int r, int v) { return pn.t[z(r)] < v; });
    return it != e && pn.t[z(*it)] == s ? *it : -1;
  }
};

// The level rows, the transformed rows and their weights on the level rows
// (Roodman 2009, eq 21 and sec. 3.1).
struct rows {
  std::vector<int> lr, lu, ls, luo;     // level rows: row, unit, period
  std::vector<int> tu, tl, wo, wq;      // transformed rows: unit, label
  std::vector<double> wv;               // weights wv on level rows wq
  std::vector<int> tuo;                 // transformed rows of each unit
};

rows mkrows(const pnl &pn, const idx &ix, const spec &sp, double &big) {
  rows w;
  w.luo.assign(z(pn.nu) + 1, 0);
  w.tuo.assign(z(pn.nu) + 1, 0);
  for (int a = 0; a < pn.nu; ++a) {
    for (int q = ix.uo[z(a)]; q < ix.uo[z(a) + 1]; ++q) {
      const int r = ix.ord[z(q)], s = pn.t[z(r)];
      bool ok = true;
      for (int l = 1; l <= sp.ylag && ok; ++l) ok = ix.at(a, s - l) >= 0;
      if (!ok) continue;
      w.lr.push_back(r);
      w.lu.push_back(a);
      w.ls.push_back(s);
    }
    w.luo[z(a) + 1] = static_cast<int>(w.lr.size());
  }
  // forward deviations weigh every later row: count them first
  double nw = 0.0;
  for (int a = 0; a < pn.nu; ++a) {
    const double l = w.luo[z(a) + 1] - w.luo[z(a)];
    nw += sp.fod ? l * (l + 1) / 2 : 2 * l;
  }
  if (nw > cap) {
    big = nw;
    return w;
  }
  w.wo.push_back(0);
  for (int a = 0; a < pn.nu; ++a) {
    const int j0 = w.luo[z(a)], j1 = w.luo[z(a) + 1];
    for (int j = j0; j < j1; ++j) {
      if (!sp.fod) {
        if (j == j0 || w.ls[z(j)] != w.ls[z(j - 1)] + 1) continue;
        w.wq.insert(w.wq.end(), {j, j - 1});
        w.wv.insert(w.wv.end(), {1.0, -1.0});
        w.tl.push_back(w.ls[z(j)]);
      } else {
        if (j == j1 - 1) continue;
        const double t = j1 - 1 - j, c = std::sqrt(t / (t + 1));
        w.wq.push_back(j);
        w.wv.push_back(c);
        for (int h = j + 1; h < j1; ++h) {
          w.wq.push_back(h);
          w.wv.push_back(-c / t);
        }
        w.tl.push_back(w.ls[z(j)] + 1);
      }
      w.tu.push_back(a);
      w.wo.push_back(static_cast<int>(w.wq.size()));
    }
    w.tuo[z(a) + 1] = static_cast<int>(w.tu.size());
  }
  return w;
}

// The value of GMM-style variable v (-1: y) of unit a in period s; false
// when the unit has no row then.
bool gval(const pnl &pn, const idx &ix, const xm &xa, int v, int a, int s,
          double &out) {
  const int r = ix.at(a, s);
  if (r < 0) return false;
  out = v < 0 ? pn.y[r] : xa.v[z(v) * z(pn.n) + z(r)];
  return true;
}

// The modified Gram-Schmidt pass over the columns of x (n x k): keep[j]
// is false for a column that is 0 or within 1e-9 of the span of the kept
// ones before it, in squared norm.
std::vector<char> indep(const std::vector<double> &x, int n, int k) {
  std::vector<char> keep(z(k), 0);
  std::vector<std::vector<double>> qs;
  for (int j = 0; j < k; ++j) {
    std::vector<double> v(x.begin() + static_cast<std::ptrdiff_t>(z(j) * z(n)),
                          x.begin() + static_cast<std::ptrdiff_t>(z(j + 1) * z(n)));
    double n0 = 0.0;
    for (double e : v) n0 += e * e;
    if (!(n0 > 0.0)) continue;
    for (int pass = 0; pass < 2; ++pass) {
      for (const auto &q : qs) {
        double d = 0.0;
        for (int i = 0; i < n; ++i) d += q[z(i)] * v[z(i)];
        for (int i = 0; i < n; ++i) v[z(i)] -= d * q[z(i)];
      }
    }
    double n1 = 0.0;
    for (double e : v) n1 += e * e;
    if (!(n1 > 1e-9 * n0)) continue;
    for (double &e : v) e /= std::sqrt(n1);
    qs.push_back(std::move(v));
    keep[z(j)] = 1;
  }
  return keep;
}

// The candidate GMM-style instruments of variable v: lags a..b of the
// level for the transformed rows (one column per lag, or per label and lag
// uncollapsed), lag a - 1 of the difference for the level rows (one
// column, or one per period).
struct gcol {
  int v, a, b, base, nl, lbase;
};

// Codes of build(): 0 done; 1 no transformed equation; 2 a regressor of
// the user is collinear (bad); 3 too many instruments or arrays too large
// (big); 4 interrupted; 5 too many regressors (big).
int build(const pnl &pn, const xm &xa, const xm &zx, const spec &sp,
          dsg &d) {
  const idx ix(pn);
  const rows w = mkrows(pn, ix, sp, d.big);
  if (d.big > 0) return 3;
  const int nt = pn.nt, p = xa.p, nlr = static_cast<int>(w.lr.size());
  const int ntr = static_cast<int>(w.tu.size());
  if (ntr == 0) return 1;
  const int kc = sp.sys && sp.cons, ndum = sp.tdum ? nt : 0, k0 = sp.ylag + p;
  const int kall = k0 + kc + ndum, qz = zx.p, nc = kall + qz;
  if (kall > mmax) {
    d.big = kall;
    return 5;
  }
  // the stack: each unit's transformed rows, then its level rows
  std::vector<int> st, sl;  // per stacked row: transformed index, or -1; level index
  for (int a = 0; a < pn.nu; ++a) {
    const int t0 = w.tuo[z(a)], t1 = w.tuo[z(a) + 1];
    const int l0 = w.luo[z(a)], l1 = sp.sys ? w.luo[z(a) + 1] : l0;
    if (t1 == t0 && l1 == l0) continue;
    d.grp.push_back(a);
    d.go.push_back(static_cast<int>(st.size()));
    for (int r = t0; r < t1; ++r) {
      st.push_back(r);
      sl.push_back(-1);
    }
    for (int q = l0; q < l1; ++q) {
      st.push_back(-1);
      sl.push_back(q);
    }
    const int ob = sp.sys ? l1 - l0 : t1 - t0;
    d.nobs += ob;
    d.omin = d.ng == 0 ? ob : std::min(d.omin, ob);
    d.omax = std::max(d.omax, ob);
    ++d.ng;
  }
  const int ns = static_cast<int>(st.size());
  d.go.push_back(ns);
  if ((static_cast<double>(ns) + nlr) * nc > cap) {
    d.big = (static_cast<double>(ns) + nlr) * nc;
    return 3;
  }
  // the level values of every column: lags of y, x, constant, dummies, z
  std::vector<double> lv(z(nlr) * z(nc), 0.0);
  for (int q = 0; q < nlr; ++q) {
    const int r = w.lr[z(q)], a = w.lu[z(q)], s = w.ls[z(q)];
    for (int c = 0; c < sp.ylag; ++c) lv[z(c) * z(nlr) + z(q)] = pn.y[ix.at(a, s - 1 - c)];
    for (int c = 0; c < p; ++c) lv[z(sp.ylag + c) * z(nlr) + z(q)] = xa.v[z(c) * z(pn.n) + z(r)];
    if (kc) lv[z(k0) * z(nlr) + z(q)] = 1.0;
    if (ndum) lv[z(k0 + kc + s) * z(nlr) + z(q)] = 1.0;
    for (int c = 0; c < qz; ++c) lv[z(kall + c) * z(nlr) + z(q)] = zx.v[z(c) * z(pn.n) + z(r)];
  }
  auto sval = [&](int i, int c) {
    if (st[z(i)] < 0) return lv[z(c) * z(nlr) + z(sl[z(i)])];
    const int t = st[z(i)];
    double v = 0.0;
    for (int e = w.wo[z(t)]; e < w.wo[z(t) + 1]; ++e) {
      v += w.wv[z(e)] * lv[z(c) * z(nlr) + z(w.wq[z(e)])];
    }
    return v;
  };
  std::vector<double> sx(z(ns) * z(nc)), sy(z(ns));
  for (int i = 0; i < ns; ++i) {
    for (int c = 0; c < nc; ++c) sx[z(c) * z(ns) + z(i)] = sval(i, c);
    if (st[z(i)] < 0) {
      sy[z(i)] = pn.y[w.lr[z(sl[z(i)])]];
    } else {
      double v = 0.0;
      for (int e = w.wo[z(st[z(i)])]; e < w.wo[z(st[z(i)]) + 1]; ++e) {
        v += w.wv[z(e)] * pn.y[w.lr[z(w.wq[z(e)])]];
      }
      sy[z(i)] = v;
    }
  }
  // the names of the columns; the user's columns must be independent
  std::vector<std::string> nm;
  for (int l = 1; l <= sp.ylag; ++l) nm.push_back("L" + std::to_string(l) + ".y");
  nm.insert(nm.end(), xa.nm.begin(), xa.nm.end());
  if (kc) nm.emplace_back("(const)");
  for (int s = 0; s < ndum; ++s) nm.push_back("t" + pv(pn.tv[z(s)]));
  if (kutil::stop(static_cast<double>(ns) * kall * kall)) return 4;
  const std::vector<char> keep = indep(sx, ns, kall);
  std::vector<int> kx;
  for (int c = 0; c < kall; ++c) {
    if (keep[z(c)]) {
      kx.push_back(c);
      if (c >= k0 + kc) d.dum.push_back(nm[z(c)]);
    } else if (c < k0 + kc) {
      d.bad = nm[z(c)];
      return 2;
    } else {
      d.drop.push_back(nm[z(c)]);
    }
  }
  const int k = static_cast<int>(kx.size());
  // the GMM-style candidates, those that are not 0 throughout marked
  std::vector<gcol> gc;
  double ncand = 0.0;
  for (size_t j = 0; j < sp.gv.size(); ++j) {
    gcol g{sp.gv[j], sp.lo[j], std::min(sp.hi[j], nt - 1), 0, 0, 0};
    g.nl = g.b >= g.a ? g.b - g.a + 1 : 0;
    g.base = static_cast<int>(ncand);
    ncand += sp.coll ? g.nl : static_cast<double>(nt) * g.nl;
    g.lbase = static_cast<int>(std::min(ncand, cap));
    ncand += sp.sys ? (sp.coll ? 1 : nt) : 0;
    gc.push_back(g);
    if (ncand > cap) break;
  }
  if (ncand > cap) {
    d.big = ncand;
    return 3;
  }
  const int nca = static_cast<int>(ncand);
  std::vector<int> col(z(nca), -1);
  // pass 0 marks the candidates used, pass 1 fills z
  for (int pass = 0; pass < 2; ++pass) {
    for (int i = 0; i < ns; ++i) {
      for (const gcol &g : gc) {
        if (st[z(i)] >= 0) {
          const int t = st[z(i)], a = w.tu[z(t)], lab = w.tl[z(t)];
          for (int l = g.a; l <= g.b && l <= lab; ++l) {
            double v = 0.0;
            if (!gval(pn, ix, xa, g.v, a, lab - l, v) || v == 0.0) continue;
            const int c = g.base + (sp.coll ? l - g.a : lab * g.nl + l - g.a);
            if (pass == 0) col[z(c)] = 0;
            else d.z[z(col[z(c)]) * z(ns) + z(i)] = v;
          }
        } else {
          const int q = sl[z(i)], a = w.lu[z(q)], s = w.ls[z(q)];
          double v1 = 0.0, v0 = 0.0;
          if (!gval(pn, ix, xa, g.v, a, s - g.a + 1, v1) ||
              !gval(pn, ix, xa, g.v, a, s - g.a, v0) || v1 == v0) {
            continue;
          }
          const int c = g.lbase + (sp.coll ? 0 : s);
          if (pass == 0) col[z(c)] = 0;
          else d.z[z(col[z(c)]) * z(ns) + z(i)] = v1 - v0;
        }
      }
    }
    if (pass == 1) break;
    if (kutil::stop(static_cast<double>(ns) * nca)) return 4;
    // the IV-style columns: the x named, z, the constant and the dummies
    std::vector<int> ivc;
    for (int j : sp.iv) ivc.push_back(sp.ylag + j);
    for (int c = k0; c < k0 + kc; ++c) ivc.push_back(c);
    for (int c = 0; c < qz; ++c) ivc.push_back(kall + c);
    for (int c = k0 + kc; c < kall; ++c) {
      if (keep[z(c)]) ivc.push_back(c);
    }
    std::vector<int> ivk;
    for (int c : ivc) {
      bool nz = false;
      for (int i = 0; i < ns && !nz; ++i) nz = sx[z(c) * z(ns) + z(i)] != 0.0;
      if (nz) ivk.push_back(c);
    }
    int m = 0;
    for (const gcol &g : gc) {
      const std::string vn = g.v < 0 ? "y" : xa.nm[z(g.v)];
      for (int c = g.base; c < g.base + (sp.coll ? g.nl : nt * g.nl); ++c) {
        if (col[z(c)] < 0) continue;
        col[z(c)] = m++;
        const int o = c - g.base, l = g.a + (sp.coll ? o : o % g.nl);
        d.zn.push_back("L" + std::to_string(l) + "." + vn +
                       (sp.coll ? "" : "@" + pv(pn.tv[z(o / g.nl)])));
        d.sub.push_back(0);
      }
      for (int c = g.lbase; c < g.lbase + (sp.sys ? (sp.coll ? 1 : nt) : 0); ++c) {
        if (col[z(c)] < 0) continue;
        col[z(c)] = m++;
        d.zn.push_back("D.L" + std::to_string(g.a - 1) + "." + vn +
                       (sp.coll ? "" : "@" + pv(pn.tv[z(c - g.lbase)])));
        d.sub.push_back(1);
      }
    }
    const int mg = m;
    m += static_cast<int>(ivk.size());
    if (m > mmax || static_cast<double>(ns) * m > cap) {
      d.big = m;
      return 3;
    }
    d.m = m;
    d.z.assign(z(ns) * z(m), 0.0);
    for (size_t j = 0; j < ivk.size(); ++j) {
      const int c = ivk[j];
      d.zn.push_back(c >= kall ? zx.nm[z(c - kall)] : nm[z(c)]);
      d.sub.push_back(0);
      std::copy(sx.begin() + static_cast<std::ptrdiff_t>(z(c) * z(ns)),
                sx.begin() + static_cast<std::ptrdiff_t>(z(c + 1) * z(ns)),
                d.z.begin() + static_cast<std::ptrdiff_t>(z(mg + static_cast<int>(j)) * z(ns)));
    }
  }
  // H, group by group (Roodman 2009, eqs 22-23 and 26)
  double hs = 0.0;
  for (int a = 0; a < d.ng; ++a) {
    const double ni = d.go[z(a) + 1] - d.go[z(a)];
    hs += ni * ni;
  }
  if (hs > cap) {
    d.big = hs;
    return 3;
  }
  d.h.reserve(static_cast<size_t>(hs));
  for (int a = 0; a < d.ng; ++a) {
    const int u = d.grp[z(a)], t0 = w.tuo[z(u)], nta = w.tuo[z(u) + 1] - t0;
    const int l0 = w.luo[z(u)], nla = w.luo[z(u) + 1] - l0, ni = d.go[z(a) + 1] - d.go[z(a)];
    std::vector<double> hb(z(ni) * z(ni), 0.0), wm(z(nta) * z(nla), 0.0);
    for (int r = 0; r < nta; ++r) {
      for (int e = w.wo[z(t0 + r)]; e < w.wo[z(t0 + r) + 1]; ++e) {
        wm[z(w.wq[z(e)] - l0) * z(nta) + z(r)] = w.wv[z(e)];
      }
    }
    for (int r = 0; r < ni; ++r) hb[z(r) * z(ni) + z(r)] = 1.0;
    if (sp.hm != 2) {
      for (int r = 0; r < nta; ++r) {
        for (int c = 0; c < nta; ++c) {
          double s = 0.0;
          for (int j = 0; j < nla; ++j) s += wm[z(j) * z(nta) + z(r)] * wm[z(j) * z(nta) + z(c)];
          hb[z(c) * z(ni) + z(r)] = s;
        }
      }
      for (int r = 0; r < nta && sp.hm == 0 && sp.sys; ++r) {
        for (int j = 0; j < nla; ++j) {
          hb[z(nta + j) * z(ni) + z(r)] = hb[z(r) * z(ni) + z(nta + j)] = wm[z(j) * z(nta) + z(r)];
        }
      }
    }
    d.h.insert(d.h.end(), hb.begin(), hb.end());
    for (int r = 0; r < ni; ++r) d.eq.push_back(r < nta);
    if (sp.sig) {
      // the covariances, over sigma^2, of the stacked rows and the
      // differenced rows with independent homoskedastic errors (Roodman
      // 2009, n. 19): the transformed rows weigh the level rows by wm, a
      // level row is itself, a difference is +1 and -1 on two of them
      std::vector<int> jd;
      for (int j = 1; j < nla; ++j) {
        if (w.ls[z(l0 + j)] == w.ls[z(l0 + j - 1)] + 1) jd.push_back(j);
      }
      const int ndg = static_cast<int>(jd.size());
      for (int s = 0; s < ndg; ++s) {
        const int j = jd[z(s)];
        for (int r = 0; r < nta; ++r) {
          d.hc.push_back(wm[z(j) * z(nta) + z(r)] - wm[z(j - 1) * z(nta) + z(r)]);
        }
        for (int r = 0; r < ni - nta; ++r) d.hc.push_back(r == j ? 1.0 : r == j - 1 ? -1.0 : 0.0);
      }
      for (int s = 0; s < ndg; ++s) {
        for (int e = 0; e < ndg; ++e) {
          const int dj = jd[z(s)] - jd[z(e)];
          d.hdd.push_back(dj == 0 ? 2.0 : dj == 1 || dj == -1 ? -1.0 : 0.0);
        }
      }
    }
    // the differenced rows of the Arellano-Bond tests
    d.dgo.push_back(static_cast<int>(d.dt.size()));
    for (int j = l0 + 1; j < l0 + nla; ++j) {
      if (w.ls[z(j)] != w.ls[z(j - 1)] + 1) continue;
      d.dt.push_back(w.ls[z(j)]);
      d.yd.push_back(pn.y[w.lr[z(j)]] - pn.y[w.lr[z(j - 1)]]);
    }
  }
  d.nd = static_cast<int>(d.dt.size());
  d.dgo.push_back(d.nd);
  d.xd.assign(z(d.nd) * z(k), 0.0);
  int r = 0;
  for (int a = 0; a < d.ng; ++a) {
    const int u = d.grp[z(a)], l0 = w.luo[z(u)], l1 = w.luo[z(u) + 1];
    for (int j = l0 + 1; j < l1; ++j) {
      if (w.ls[z(j)] != w.ls[z(j - 1)] + 1) continue;
      for (int c = 0; c < k; ++c) {
        const size_t o = z(kx[z(c)]) * z(nlr);
        d.xd[z(c) * z(d.nd) + z(r)] = lv[o + z(j)] - lv[o + z(j - 1)];
      }
      ++r;
    }
  }
  d.lrow = w.lr;
  d.ly.resize(z(nlr));
  d.lx.assign(z(nlr) * z(k), 0.0);
  for (int q = 0; q < nlr; ++q) {
    d.ly[z(q)] = pn.y[w.lr[z(q)]];
    for (int c = 0; c < k; ++c) d.lx[z(c) * z(nlr) + z(q)] = lv[z(kx[z(c)]) * z(nlr) + z(q)];
  }
  d.n = ns;
  d.k = k;
  d.y = sy;
  d.x.assign(z(ns) * z(k), 0.0);
  for (int c = 0; c < k; ++c) {
    std::copy(sx.begin() + static_cast<std::ptrdiff_t>(z(kx[z(c)]) * z(ns)),
              sx.begin() + static_cast<std::ptrdiff_t>(z(kx[z(c)] + 1) * z(ns)),
              d.x.begin() + static_cast<std::ptrdiff_t>(z(c) * z(ns)));
    d.xn.push_back(nm[z(kx[z(c)])]);
  }
  return 0;
}


// A coefficient table with normal intervals, or t intervals on df > 0.
void ctab(tab &t, const std::vector<std::string> &nm, const double *b, const double *v, int k,
          double lv, double df) {
  const double q = df > 0 ? clab_qt(1 - (1 - lv) / 2, df, 1) : clab_qnorm(1 - (1 - lv) / 2);
  t = tab();
  for (int j = 0; j < k; ++j) {
    const double se = std::sqrt(v[z(j) * z(k) + z(j)]), s = b[j] / se;
    t.name.push_back(nm[z(j)]);
    t.est.push_back(b[j]);
    t.se.push_back(se);
    t.stat.push_back(s);
    t.p.push_back(df > 0 ? 2 * clab_pt(-std::fabs(s), df, 1) : 2 * clab_pnorm(-std::fabs(s)));
    t.lo.push_back(b[j] - q * se);
    t.hi.push_back(b[j] + q * se);
  }
}

// The joint Wald test of the coefficients other than the constant on
// dfm degrees of freedom: chi-squared, or with df2 > 0 F = W / dfm on dfm
// and df2; NaN when dfm is below 1.
void wtest(tab &t, const std::vector<std::string> &nm, const std::vector<double> &b,
           const std::vector<double> &v, int k, int dfm, double df2) {
  std::vector<int> j;
  for (int i = 0; i < k; ++i) {
    if (nm[z(i)] != "(const)") j.push_back(i);
  }
  const int q = static_cast<int>(j.size());
  std::vector<double> a(z(q) * z(q)), ai(z(q) * z(q));
  for (int r = 0; r < q; ++r) {
    for (int c = 0; c < q; ++c) a[z(c) * z(q) + z(r)] = v[z(j[z(c)]) * z(k) + z(j[z(r)])];
  }
  double w = NAN;
  if (q > 0 && dfm > 0 && pev::chol(a.data(), q)) {
    pev::cinv(a.data(), q, ai.data());
    w = 0.0;
    for (int r = 0; r < q; ++r) {
      for (int c = 0; c < q; ++c) w += b[z(j[z(r)])] * ai[z(c) * z(q) + z(r)] * b[z(j[z(c)])];
    }
  }
  t = tab();
  t.name.push_back(df2 > 0 ? "F" : "Wald");
  t.stat.push_back(df2 > 0 ? w / dfm : w);
  t.df.push_back(dfm);
  t.df2.push_back(df2 > 0 ? df2 : NAN);
  t.p.push_back(std::isnan(w)   ? NAN
                : df2 > 0       ? clab_pf(w / dfm, dfm, df2, 0)
                                : clab_pchisq(w, dfm, 0));
}

// The model degrees of freedom of the joint test as xtabond2 3.7.2
// computes them (Roodman 2009): the rank of the regressors in the levels
// (system GMM, whose transformed equation the program leaves out of the
// test) or in the transformed equation, less one for the constant of
// system GMM; without that constant, one less again when the vector of
// ones lies in their column space, which the program finds by comparing
// i'X_p (X'X)^-1 X_p'i with the rows of X_p, X the whole stack, to the
// relative tolerance of the machine epsilon times the groups times the
// periods.
int mdf(const dsg &d, bool sys, int cons, int nt) {
  std::vector<int> rr;
  for (int r = 0; r < d.n; ++r) {
    if (!sys || !d.eq[z(r)]) rr.push_back(r);
  }
  const int np = static_cast<int>(rr.size()), k = d.k;
  std::vector<double> px(z(np) * z(k)), xi(z(k), 0.0);
  for (int c = 0; c < k; ++c) {
    for (int i = 0; i < np; ++i) {
      const double e = d.x[z(c) * z(d.n) + z(rr[z(i)])];
      px[z(c) * z(np) + z(i)] = e;
      xi[z(c)] += e;
    }
  }
  int dfm = -cons;
  for (char c : indep(px, np, k)) dfm += c;
  if (cons) return dfm;
  std::vector<double> xx(z(k) * z(k), 0.0), xv(z(k) * z(k));
  for (int c = 0; c < k; ++c) {
    for (int e = 0; e < k; ++e) {
      double s = 0.0;
      for (int i = 0; i < d.n; ++i) s += d.x[z(c) * z(d.n) + z(i)] * d.x[z(e) * z(d.n) + z(i)];
      xx[z(e) * z(k) + z(c)] = s;
    }
  }
  if (!pev::chol(xx.data(), k)) return dfm;
  pev::cinv(xx.data(), k, xv.data());
  double q = 0.0;
  for (int c = 0; c < k; ++c) {
    for (int e = 0; e < k; ++e) q += xi[z(c)] * xv[z(e) * z(k) + z(c)] * xi[z(e)];
  }
  const double tol = DBL_EPSILON * static_cast<double>(d.ng) * nt;
  if (std::fabs(q - np) / (np + 1.0) < tol) --dfm;
  return dfm;
}

// The tests: AR(1..nar), Sargan, Hansen and, for system GMM, the
// difference test of the reported overidentification statistic (Sargan's
// with the conventional one-step variance, Hansen's otherwise).
void ttab(tab &t, const double *st, const double *ar, int nar, int m, int k, int msub, bool sys,
          bool sarg) {
  t = tab();
  for (int l = 0; l < nar; ++l) {
    t.name.push_back("AR(" + std::to_string(l + 1) + ")");
    t.stat.push_back(ar[l]);
    t.df.push_back(NAN);
    t.p.push_back(std::isnan(ar[l]) ? NAN : 2 * clab_pnorm(-std::fabs(ar[l])));
  }
  const double hx = sarg ? st[5] : st[2], dif = (sarg ? st[0] : st[1]) - hx;
  const char *nm[] = {"Sargan", "Hansen",
                      sarg ? "Sargan without level GMM" : "Hansen without level GMM",
                      sarg ? "Difference-in-Sargan" : "Difference-in-Hansen"};
  const double sv[] = {st[0], st[1], hx, dif};
  const double dv[] = {static_cast<double>(m - k), static_cast<double>(m - k),
                       static_cast<double>(m - msub - k), static_cast<double>(msub)};
  for (int j = 0; j < (sys ? 4 : 2); ++j) {
    t.name.push_back(nm[j]);
    t.stat.push_back(sv[j]);
    t.df.push_back(dv[j]);
    t.p.push_back(std::isnan(sv[j]) || dv[j] < 1 ? NAN : clab_pchisq(sv[j], dv[j], 0));
  }
}

}  // namespace

std::string gcheck(const spec &sp) {
  if (sp.wc < 0) return "";
  if (sp.steps == 1) return "`wc` corrects the two-step variance: give it with `steps = 2`.";
  if ((sp.wc != 0) != (sp.rob != 0)) {
    return fmt("`wc = %s` contradicts `robust = %s`: the two-step variance is Windmeijer's with "
               "`robust = TRUE` and uncorrected with `robust = FALSE`.",
               sp.wc ? "TRUE" : "FALSE", sp.rob ? "TRUE" : "FALSE");
  }
  return "";
}

std::string gfit(const pnl &pn, const xm &x, const xm &zx, spec sp, gout &o, int &code) {
  code = CLAB_OK;
  sp.sig = sp.steps == 1 && !sp.rob;
  dsg d;
  const int bc = build(pn, x, zx, sp, d);
  if (bc == 4) {
    code = CLAB_EINTR;
    return "";
  }
  if (bc == 1) return "no unit has the rows the transformed equation needs (two usable periods).";
  if (bc == 2) {
    return "the regressor " + d.bad +
           " is constant over time or collinear with the others in the stacked equations.";
  }
  if (bc == 3) {
    return fmt("the instruments or the arrays would be too large (%.0f; at most %d instruments): "
               "collapse them or narrow the lags.",
               d.big, mmax);
  }
  if (bc == 5) return fmt("the regressors would number %.0f, more than %d.", d.big, mmax);
  if (d.m < d.k) return fmt("the %d instruments are fewer than the %d coefficients.", d.m, d.k);
  // small, as xtabond2 3.7.2 computes it: the rows n of the equation
  // reported, the rows w of the transformed equation in system GMM with H
  // other than the identity (n otherwise), the groups G and the constant c
  // of system GMM
  int ntr = 0;
  for (int e : d.eq) ntr += e;
  const double wt = sp.sys && sp.hm != 2 ? ntr : d.nobs;
  if (sp.small && d.nobs <= d.k) {
    return fmt("`small = TRUE` needs more rows (%d) than coefficients (%d).", d.nobs, d.k);
  }
  if (sp.small && wt <= d.k) {
    return fmt("`small = TRUE` needs more rows in the transformed equation (%.0f) than "
               "coefficients (%d).",
               wt, d.k);
  }
  if (sp.small && !sp.sig && d.ng < 2) {
    return "`small = TRUE` with the robust or two-step variance needs two groups or more.";
  }
  const double fs = sp.small ? wt / (wt - d.k) : 1.0;
  int msub = 0;
  for (int f : d.sub) msub += f;
  std::vector<double> b(z(d.k)), v(z(d.k) * z(d.k)), st(6), av(z(sp.nar) + 1);
  const int rc = clab_gmmpnl(d.n, d.k, d.m, d.ng, d.go.data(), d.y.data(), d.x.data(), d.z.data(),
                             d.h.data(), d.eq.data(), d.nd, d.dgo.data(), d.dt.data(),
                             d.yd.data(), d.xd.data(), sp.sig ? d.hc.data() : nullptr,
                             sp.sig ? d.hdd.data() : nullptr,
                             sp.sys && msub ? d.sub.data() : nullptr, sp.steps, sp.rob, sp.nar,
                             sp.sig ? fs : 1.0, b.data(), v.data(), st.data(), av.data());
  if (rc == CLAB_EINTR || clab_intrd(0) != 0) {
    code = CLAB_EINTR;
    return "";
  }
  if (rc == CLAB_EDOM) {
    return "the instruments do not identify the coefficients, or the moments of the residuals "
           "are all 0.";
  }
  if (rc != CLAB_OK) return "the work arrays of the estimation would be too large.";
  if (d.m > d.ng) {
    o.note.push_back(std::to_string(d.m) + " instruments for " + std::to_string(d.ng) +
                     " groups: more instruments than groups weakens the Hansen test "
                     "(Roodman 2009, sec. 2.6).");
  }
  if (st[4] > 0) {
    o.note.push_back("a weight matrix was singular: a generalized inverse was used, and the "
                     "difference-in-Hansen statistic may be negative.");
  }
  // small (Roodman 2009, p. 123: "small requests small-sample corrections
  // to the covariance matrix estimate, resulting in t-test instead of
  // z-test statistics for the coefficients and an F test instead of a
  // Wald chi2 test for overall fit"), as xtabond2 3.7.2 computes it: the
  // conventional one-step variance times w / (w - k), any other times
  // (n - 1) / (n - k) G / (G - 1); sigma^2 times w / (w - k); the
  // residual degrees of freedom n - dfm - c for the conventional one-step
  // variance and G - c otherwise
  const int cons = sp.sys && sp.cons, dfm = mdf(d, sp.sys != 0, cons, pn.nt);
  double dfr = 0.0;
  if (sp.small) {
    const double n1 = d.nobs, g1 = d.ng, fv = sp.sig ? fs : (n1 - 1) / (n1 - d.k) * g1 / (g1 - 1);
    for (double &e : v) e *= fv;
    dfr = (sp.sig ? n1 - dfm : g1) - cons;
  }
  ctab(o.coef, d.xn, b.data(), v.data(), d.k, sp.level, dfr);
  // the statistic as xtabond2 3.7.2 prints it: b' V^-1 b over every
  // coefficient, the constant's included, on the same degrees of freedom
  {
    std::vector<std::string> nx(d.xn.size(), "x");
    wtest(o.xwald, nx, b, v, d.k, dfm, dfr);
  }
  const int nlr = static_cast<int>(d.lrow.size());
  o.resid.assign(z(pn.n), NAN);
  for (int q = 0; q < nlr; ++q) {
    double e = d.ly[z(q)];
    for (int c = 0; c < d.k; ++c) e -= d.lx[z(c) * z(nlr) + z(q)] * b[z(c)];
    o.resid[z(d.lrow[z(q)])] = e;
  }
  o.tstat = sp.small != 0;
  ttab(o.tests, st.data(), av.data(), sp.nar, d.m, d.k, msub, sp.sys != 0, sp.sig != 0);
  wtest(o.wald, d.xn, b, v, d.k, dfm, dfr);
  o.n = d.nobs;
  o.groups = d.ng;
  o.inst = d.m;
  o.omin = d.omin;
  o.oavg = static_cast<double>(d.nobs) / d.ng;
  o.omax = d.omax;
  o.steps = sp.steps;
  o.k = d.k;
  o.sigma2 = st[3] * fs;
  o.vcov = v;
  o.inst_names = d.zn;
  o.dropped = d.drop;
  o.dummies = d.dum;
  o.type = sp.sys ? "system" : "diff";
  o.trans = sp.fod ? "fod" : "fd";
  return "";
}

std::string wfit(const pnl &pn, const xm &x, int cons, int raw, wout &o, int &code) {
  code = CLAB_OK;
  if (x.p > mmax) return fmt("`x` must have at most %d columns.", mmax);
  std::vector<double> vb(z(x.p) * z(x.p)), r(6);
  o.b.assign(z(x.p), 0.0);
  const int rc = clab_wsctst(pn.n, x.p, pn.y, x.v.data(), pn.u.data(), pn.t.data(), cons, raw,
                             o.b.data(), vb.data(), r.data());
  if (rc == CLAB_EINTR || clab_intrd(0) != 0) {
    code = CLAB_EINTR;
    return "";
  }
  if (rc == CLAB_EDOM) {
    return "too few consecutive first differences, or units with them, for the test, or the "
           "differences of `x` are collinear.";
  }
  if (rc != CLAB_OK) return "the work arrays of the test would be too large.";
  const double s[] = {r[2], 1,    r[3] - 1, clab_pf(r[2], 1, r[3] - 1, 0), r[0],
                      r[1], r[4], r[3],     r[5]};
  std::copy(s, s + 9, o.stat);
  o.se.clear();
  for (int j = 0; j < x.p; ++j) o.se.push_back(std::sqrt(vb[z(j) * z(x.p) + z(j)]));
  return "";
}

namespace {

// The largest whole p with p^3 at most t: the integer part of t^(1/3).
int icbrt(int t) {
  int p = 0;
  while (static_cast<double>(p + 1) * (p + 1) * (p + 1) <= t) ++p;
  return p;
}

// The eigenvalues (decreasing) and vectors of the symmetric m x m matrix
// a; false when the decomposition fails.
bool eig(std::vector<double> a, int m, std::vector<double> &val, std::vector<double> &vec) {
  val.assign(z(m), 0.0);
  vec.assign(z(m) * z(m), 0.0);
  return gla::jacobi(a.data(), m, val.data(), vec.data());
}

// The number of factors the averages kept identify, as Juodis and Reese
// (2026, app. A.1) set out the two regularized procedures, on a balanced
// panel: out[0] the eigenvalue ratio of Ahn and Horenstein (2013) on the
// normalized proxies Zbar Sigma^-1/2, Sigma = sum_i (Z_i - Zbar)'(Z_i -
// Zbar) / (N T) (Juodis 2022; steps 1-3, r_max the averages less one);
// out[1] the eigenvalues of Zbar'Zbar / T at least 0.01 of the largest
// (Vogt et al. 2025; steps 1-2); out[2] the averages. NaN where undefined.
void nfac(const std::vector<double> &yy, const std::vector<double> &xx,
          const std::vector<double> &wv, int nu, int nt, const std::vector<int> &av, double *out) {
  std::vector<int> cs;
  for (size_t c = 0; c < av.size(); ++c) {
    if (av[c]) cs.push_back(static_cast<int>(c) - 1);
  }
  const int kz = static_cast<int>(cs.size());
  out[0] = out[1] = NAN;
  out[2] = kz;
  if (kz == 0) return;
  auto val = [&](int c, int i, int t) {
    return c < 0 ? yy[z(i) * z(nt) + z(t)] : xx[(z(c) * z(nu) + z(i)) * z(nt) + z(t)];
  };
  std::vector<double> zb(z(nt) * z(kz), 0.0), sg(z(kz) * z(kz), 0.0), m2(z(kz) * z(kz), 0.0);
  for (int c = 0; c < kz; ++c) {
    for (int i = 0; i < nu; ++i) {
      for (int t = 0; t < nt; ++t) zb[z(c) * z(nt) + z(t)] += wv[z(i)] * val(cs[z(c)], i, t);
    }
  }
  for (int i = 0; i < nu; ++i) {
    for (int t = 0; t < nt; ++t) {
      for (int a = 0; a < kz; ++a) {
        const double ea = val(cs[z(a)], i, t) - zb[z(a) * z(nt) + z(t)];
        for (int b = 0; b < kz; ++b) {
          const double eb = val(cs[z(b)], i, t) - zb[z(b) * z(nt) + z(t)];
          sg[z(b) * z(kz) + z(a)] += ea * eb / (static_cast<double>(nu) * nt);
        }
      }
    }
  }
  for (int a = 0; a < kz; ++a) {
    for (int b = 0; b < kz; ++b) {
      double s = 0.0;
      for (int t = 0; t < nt; ++t) s += zb[z(a) * z(nt) + z(t)] * zb[z(b) * z(nt) + z(t)];
      m2[z(b) * z(kz) + z(a)] = s / nt;
    }
  }
  std::vector<double> v1, e1, v2, e2;
  if (eig(m2, kz, v2, e2) && v2[0] > 0.0) {
    int r = 0;
    for (double e : v2) r += e >= 0.01 * v2[0];
    out[1] = r;
  }
  if (kz < 2 || !eig(sg, kz, v1, e1) || !(v1[z(kz) - 1] > 0.0)) return;
  // Sigma^-1/2 = V diag(1 / sqrt(val)) V', then Sigma^-1/2 M2 Sigma^-1/2
  std::vector<double> si(z(kz) * z(kz), 0.0), f(z(kz) * z(kz), 0.0), g(z(kz) * z(kz), 0.0);
  for (int e = 0; e < kz; ++e) {
    const double r = 1.0 / std::sqrt(v1[z(e)]);
    for (int a = 0; a < kz; ++a) {
      for (int b = 0; b < kz; ++b) {
        si[z(b) * z(kz) + z(a)] += r * e1[z(e) * z(kz) + z(a)] * e1[z(e) * z(kz) + z(b)];
      }
    }
  }
  for (int a = 0; a < kz; ++a) {
    for (int b = 0; b < kz; ++b) {
      double s = 0.0;
      for (int c = 0; c < kz; ++c) s += m2[z(c) * z(kz) + z(a)] * si[z(b) * z(kz) + z(c)];
      f[z(b) * z(kz) + z(a)] = s;
    }
  }
  for (int a = 0; a < kz; ++a) {
    for (int b = 0; b < kz; ++b) {
      double s = 0.0;
      for (int c = 0; c < kz; ++c) s += si[z(c) * z(kz) + z(a)] * f[z(b) * z(kz) + z(c)];
      g[z(b) * z(kz) + z(a)] = s;
    }
  }
  if (!eig(g, kz, v1, e1)) return;
  int best = 1;
  double top = -1.0;
  for (int r = 1; r < kz; ++r) {
    const double q = v1[z(r)] > 0.0 ? v1[z(r) - 1] / v1[z(r)] : HUGE_VAL;
    if (q > top) {
      top = q;
      best = r;
    }
  }
  out[0] = best;
}

// a b for k x k matrices, column-major
std::vector<double> mmul(const std::vector<double> &a, const std::vector<double> &b, int k) {
  std::vector<double> c(z(k) * z(k), 0.0);
  for (int j = 0; j < k; ++j) {
    for (int h = 0; h < k; ++h) {
      const double v = b[z(j) * z(k) + z(h)];
      for (int i = 0; i < k; ++i) c[z(j) * z(k) + z(i)] += a[z(h) * z(k) + z(i)] * v;
    }
  }
  return c;
}

// The variance of the pooled estimate b (equal weights, balanced) that
// Brown, Schmidt and Wooldridge (2026, eqs 36-40) give for fixed T with
// Psi = H: A^-1 B A^-1 / N, A = sum X_i'M X_i / N, B = sum s_i s_i' / N,
// s_i = X_i'M u_i - G'(I + K)[Psi (Psi'Psi)^-1 (x) M] q_i, u_i = y_i - X_i
// b, G = sum u_i (x) (X_i - Xbar) / N and q_i = vec(Z_i - Zbar), zero in
// the columns of the observed factors. With (B (x) A) vec(Q) = vec(A Q B')
// and (I + K) vec(W) = vec(W + W'), G' vec(S) is sum_m (X_m - Xbar)' S
// u_m / N. False when Psi'Psi or A is singular.
bool bswv(const std::vector<double> &yy, const std::vector<double> &xx,
          const std::vector<double> &d, int q, const std::vector<int> &av, int nu, int nt,
          int k, const std::vector<double> &b, std::vector<double> &vp) {
  std::vector<int> cs;  // the columns of H: -2 - c for factor c, -1 y, j for x
  for (int c = 0; c < q; ++c) cs.push_back(-2 - c);
  for (size_t c = 0; c < av.size(); ++c) {
    if (av[c]) cs.push_back(static_cast<int>(c) - 1);
  }
  const int hp = static_cast<int>(cs.size());
  auto xv = [&](int j, int i, int t) { return xx[(z(j) * z(nu) + z(i)) * z(nt) + z(t)]; };
  std::vector<double> yb(z(nt), 0.0), xb(z(nt) * z(k), 0.0), h(z(nt) * z(hp));
  for (int i = 0; i < nu; ++i) {
    for (int t = 0; t < nt; ++t) {
      yb[z(t)] += yy[z(i) * z(nt) + z(t)] / nu;
      for (int j = 0; j < k; ++j) xb[z(j) * z(nt) + z(t)] += xv(j, i, t) / nu;
    }
  }
  for (int c = 0; c < hp; ++c) {
    for (int t = 0; t < nt; ++t) {
      const int e = cs[z(c)];
      h[z(c) * z(nt) + z(t)] = e <= -2 ? d[z(-2 - e) * z(nt) + z(t)]
                               : e < 0 ? yb[z(t)]
                                       : xb[z(e) * z(nt) + z(t)];
    }
  }
  // P = H (H'H)^-1 (T x hp), M = I - P H'
  std::vector<double> hh(z(hp) * z(hp)), hi(z(hp) * z(hp)), pm(z(nt) * z(hp), 0.0),
      mm(z(nt) * z(nt));
  for (int a = 0; a < hp; ++a) {
    for (int c = 0; c < hp; ++c) {
      double s = 0.0;
      for (int t = 0; t < nt; ++t) s += h[z(a) * z(nt) + z(t)] * h[z(c) * z(nt) + z(t)];
      hh[z(c) * z(hp) + z(a)] = s;
    }
  }
  if (!pev::chol(hh.data(), hp)) return false;
  pev::cinv(hh.data(), hp, hi.data());
  for (int c = 0; c < hp; ++c) {
    for (int a = 0; a < hp; ++a) {
      const double e = hi[z(c) * z(hp) + z(a)];
      for (int t = 0; t < nt; ++t) pm[z(c) * z(nt) + z(t)] += h[z(a) * z(nt) + z(t)] * e;
    }
  }
  for (int s1 = 0; s1 < nt; ++s1) {
    for (int t = 0; t < nt; ++t) {
      double v = t == s1 ? 1.0 : 0.0;
      for (int c = 0; c < hp; ++c) v -= pm[z(c) * z(nt) + z(t)] * h[z(c) * z(nt) + z(s1)];
      mm[z(s1) * z(nt) + z(t)] = v;
    }
  }
  // the residuals, G (T x T x k as gm[(j T + s) T + t]) and A
  std::vector<double> u(z(nu) * z(nt)), gm(z(k) * z(nt) * z(nt), 0.0), am(z(k) * z(k), 0.0),
      mx(z(nt) * z(k));
  for (int i = 0; i < nu; ++i) {
    for (int t = 0; t < nt; ++t) {
      double e = yy[z(i) * z(nt) + z(t)];
      for (int j = 0; j < k; ++j) e -= xv(j, i, t) * b[z(j)];
      u[z(i) * z(nt) + z(t)] = e;
    }
    for (int j = 0; j < k; ++j) {
      for (int s1 = 0; s1 < nt; ++s1) {
        const double us = u[z(i) * z(nt) + z(s1)] / nu;
        for (int t = 0; t < nt; ++t) {
          gm[(z(j) * z(nt) + z(s1)) * z(nt) + z(t)] += (xv(j, i, t) - xb[z(j) * z(nt) + z(t)]) * us;
        }
      }
    }
  }
  std::vector<double> bm(z(k) * z(k), 0.0), si(z(k)), qm(z(nt) * z(hp)), w(z(nt) * z(nt)),
      mq(z(nt) * z(hp));
  for (int i = 0; i < nu; ++i) {
    if (kutil::stop(static_cast<double>(nt) * nt * (hp + k))) return false;
    for (int j = 0; j < k; ++j) {
      for (int t = 0; t < nt; ++t) {
        double s = 0.0;
        for (int r = 0; r < nt; ++r) s += mm[z(r) * z(nt) + z(t)] * xv(j, i, r);
        mx[z(j) * z(nt) + z(t)] = s;
      }
    }
    for (int j = 0; j < k; ++j) {
      for (int g = 0; g < k; ++g) {
        double s = 0.0;
        for (int t = 0; t < nt; ++t) s += mx[z(j) * z(nt) + z(t)] * xv(g, i, t);
        am[z(g) * z(k) + z(j)] += s / nu;
      }
    }
    // Q_i, then W = M Q_i P' and S = W + W'
    for (int c = 0; c < hp; ++c) {
      const int e = cs[z(c)];
      for (int t = 0; t < nt; ++t) {
        qm[z(c) * z(nt) + z(t)] = e <= -2 ? 0.0
                                  : e < 0 ? yy[z(i) * z(nt) + z(t)] - yb[z(t)]
                                          : xv(e, i, t) - xb[z(e) * z(nt) + z(t)];
      }
    }
    for (int c = 0; c < hp; ++c) {
      for (int t = 0; t < nt; ++t) {
        double s = 0.0;
        for (int r = 0; r < nt; ++r) s += mm[z(r) * z(nt) + z(t)] * qm[z(c) * z(nt) + z(r)];
        mq[z(c) * z(nt) + z(t)] = s;
      }
    }
    for (int s1 = 0; s1 < nt; ++s1) {
      for (int t = 0; t < nt; ++t) {
        double s = 0.0;
        for (int c = 0; c < hp; ++c) s += mq[z(c) * z(nt) + z(t)] * pm[z(c) * z(nt) + z(s1)];
        w[z(s1) * z(nt) + z(t)] = s;
      }
    }
    for (int j = 0; j < k; ++j) {
      double s = 0.0;
      for (int t = 0; t < nt; ++t) s += mx[z(j) * z(nt) + z(t)] * u[z(i) * z(nt) + z(t)];
      double g = 0.0;
      for (int s1 = 0; s1 < nt; ++s1) {
        for (int t = 0; t < nt; ++t) {
          g += (w[z(s1) * z(nt) + z(t)] + w[z(t) * z(nt) + z(s1)]) *
               gm[(z(j) * z(nt) + z(s1)) * z(nt) + z(t)];
        }
      }
      si[z(j)] = s - g;
    }
    for (int j = 0; j < k; ++j) {
      for (int g = 0; g < k; ++g) bm[z(g) * z(k) + z(j)] += si[z(j)] * si[z(g)] / nu;
    }
  }
  std::vector<double> ai(z(k) * z(k));
  if (!pev::chol(am.data(), k)) return false;
  pev::cinv(am.data(), k, ai.data());
  vp.assign(z(k) * z(k), 0.0);
  const std::vector<double> t1 = mmul(ai, bm, k), t3 = mmul(t1, ai, k);
  for (size_t e = 0; e < vp.size(); ++e) vp[e] = t3[e] / nu;
  return true;
}

}  // namespace

namespace {

// A fit of the dynamic kernel on the periods [t0, t1): its unit
// estimates (k x nu), mean group estimate and variance, pooled estimate
// and variance, the heterogeneity-bias statistic and the rows used.
struct dfit {
  std::vector<double> bi, bmg, vmg, bp, vp;
  double hb[2] = {0, 0};
  int n = 0;
};

// The regressors of the dynamic model (the lags 1..L of y, then x) and
// the rows of each unit's regression in [t0, t1): a row needs its unit
// observed then and in the L periods before, and p periods of averages
// before it, all within [t0, t1) (Chudik and Pesaran 2015, eq 24; each
// half of the jackknife as a panel of its own). need: the rows a unit
// must have; a unit with fewer is named in bad.
struct dset {
  int nu = 0, nt = 0, k = 0, ll = 0, p = 0, q = 0;
  const std::vector<double> *yy = nullptr, *xx = nullptr, *wv = nullptr, *d = nullptr;
  std::vector<int> av;
};

int drun(const dset &s, int t0, int t1, int need, int pool, dfit &f, int &bad, int &sht) {
  const int nu = s.nu, nt = s.nt, kd = s.ll + s.k;
  const std::vector<double> &yy = *s.yy;
  std::vector<double> xk(z(kd) * z(nu) * z(nt), NAN);
  std::vector<int> rm(z(nu) * z(nt), 0), av2(z(kd) + 1, 0);
  av2[0] = s.av[0];
  for (int j = 0; j < s.k; ++j) av2[z(s.ll + j) + 1] = s.av[z(j) + 1];
  f.n = 0;
  sht = -1;
  for (int i = 0; i < nu; ++i) {
    int ni = 0;
    for (int t = 0; t < nt; ++t) {
      const size_t r = z(i) * z(nt) + z(t);
      for (int j = 0; j < s.k; ++j) {
        const size_t c = z(i) * z(nt) + z(t);
        xk[z(s.ll + j) * z(nu) * z(nt) + c] = (*s.xx)[z(j) * z(nu) * z(nt) + c];
      }
      for (int l = 1; l <= s.ll && t - l >= 0; ++l) {
        xk[(z(l - 1) * z(nu) + z(i)) * z(nt) + z(t)] = yy[r - z(l)];
      }
      bool ok = !std::isnan(yy[r]) && t >= t0 && t < t1 && t - s.p >= t0 && t - s.ll >= t0;
      for (int l = 1; l <= s.ll && ok; ++l) ok = !std::isnan(yy[r - z(l)]);
      rm[r] = ok;
      ni += ok;
    }
    f.n += ni;
    if (ni < need && sht < 0) {
      sht = i;
      bad = ni;
    }
  }
  if (sht >= 0) return CLAB_EDOM;
  f.bi.assign(z(kd) * z(nu), 0.0);
  f.bmg.assign(z(kd), 0.0);
  f.vmg.assign(z(kd) * z(kd), 0.0);
  f.bp.assign(z(kd), 0.0);
  f.vp.assign(z(kd) * z(kd), 0.0);
  bad = -1;
  return clab_ccedyn(nu, nt, kd, yy.data(), xk.data(), s.wv->data(), s.d->data(), s.q, av2.data(),
                     s.p, rm.data(), pool, f.bi.data(), f.bmg.data(), f.vmg.data(), f.bp.data(),
                     f.vp.data(), f.hb, &bad);
}

// The mean of the unit estimates b (k x nu) and its variance (Pesaran
// 2006, eqs 53 and 58; Chudik and Pesaran 2015, eq 32).
void mgv(const std::vector<double> &b, int k, int nu, std::vector<double> &m,
         std::vector<double> &v) {
  m.assign(z(k), 0.0);
  v.assign(z(k) * z(k), 0.0);
  for (int i = 0; i < nu; ++i) {
    for (int j = 0; j < k; ++j) m[z(j)] += b[z(i) * z(k) + z(j)] / nu;
  }
  for (int i = 0; i < nu; ++i) {
    for (int j = 0; j < k; ++j) {
      for (int g = 0; g < k; ++g) {
        const double dj = b[z(i) * z(k) + z(j)] - m[z(j)], dg = b[z(i) * z(k) + z(g)] - m[z(g)];
        v[z(g) * z(k) + z(j)] += dj * dg / (nu * (nu - 1.0));
      }
    }
  }
}

// Recursive mean adjustment (Chudik and Pesaran 2015, sec. 4.1.2 and n.
// 13): each value less the mean of the unit's values in the periods
// before it; a unit's first period has none and drops, and so does the
// first period of the panel, which no unit keeps.
void rma(const std::vector<double> &a, int nu, int nt, int nc, std::vector<double> &o,
         const std::vector<double> &yy) {
  o.assign(z(nc) * z(nu) * z(nt - 1), NAN);
  for (int c = 0; c < nc; ++c) {
    for (int i = 0; i < nu; ++i) {
      double sm = 0.0;
      int cn = 0;
      for (int t = 0; t < nt; ++t) {
        if (std::isnan(yy[z(i) * z(nt) + z(t)])) continue;
        const double v = a[(z(c) * z(nu) + z(i)) * z(nt) + z(t)];
        if (cn > 0) o[(z(c) * z(nu) + z(i)) * z(nt - 1) + z(t - 1)] = v - sm / cn;
        sm += v;
        ++cn;
      }
    }
  }
}

}  // namespace

std::string cfit(const pnl &pn, const xm &x, const xm &fm, int cons, const std::vector<double> &w,
                 bool useall, const std::vector<std::string> &avg, const copt &co, double level,
                 cres &o, int &code) {
  code = CLAB_OK;
  const int n = pn.n, nu = pn.nu, nt = pn.nt, k = x.p, ll = co.ylag;
  if (k > mmax) return fmt("`x` must have at most %d columns.", mmax);
  const int q = cons + fm.p;
  // the averages kept: y, then the columns of x
  std::vector<int> av(z(k) + 1, useall ? 1 : 0);
  for (const std::string &v : avg) {
    const auto it = std::find(x.nm.begin(), x.nm.end(), v);
    if (v != "y" && it == x.nm.end()) {
      return "`avg` names " + v + ", which is neither y nor a column of `x`.";
    }
    av[v == "y" ? 0 : z(static_cast<int>(it - x.nm.begin()) + 1)] = 1;
  }
  int na = 0;
  for (int f : av) na += f;
  if (q + na == 0) {
    return "nothing is taken out: keep an average in `avg`, or give `cons` or `fac`.";
  }
  const int p = co.lags < 0 ? (ll > 0 ? icbrt(nt) : 0) : co.lags;
  if (p >= nt - 1) return fmt("`lags` must be below the periods less one (%d).", nt - 1);
  const bool dyn = ll > 0 || p > 0 || co.bias != 0;
  const int pool = co.pool < 0 ? (ll == 0 && co.bias == 0) : co.pool;
  if (pool && ll > 0) {
    return "`pool = TRUE` with `ylag` above 0: with a lagged dependent variable and "
           "heterogeneous slopes pooled estimators are inconsistent (Chudik and Pesaran 2015, n. "
           "5, after Pesaran and Smith 1995), and the theory of dynamic common correlated effects "
           "is that of the mean group estimator.";
  }
  if (pool && co.bias) {
    return "`pool = TRUE` with `bias`: the corrections of Chudik and Pesaran (2015, sec. 4.1) "
           "are of the mean group estimator.";
  }
  if (co.pvar && (!pool || dyn || !w.empty() || n != nu * nt)) {
    return "`pvar = \"bsw\"` is the variance of Brown, Schmidt and Wooldridge (2026, sec. 5.4) for "
           "the static pooled estimator with equal weights on a balanced panel: give it without "
           "`ylag`, `lags`, `bias`, `w` or `pool = FALSE`, and with every unit in every period.";
  }
  const double cells = static_cast<double>(nu) * nt;
  if (cells * (k + ll + 1) > cap) {
    return fmt("the units by periods (%.0f cells) are too many for the averages.", cells);
  }
  std::vector<int> tu(z(nu), 0), ut(z(nt), 0);
  for (int i = 0; i < n; ++i) {
    ++tu[z(pn.u[z(i)])];
    ++ut[z(pn.t[z(i)])];
  }
  for (int t = 0; t < nt; ++t) {
    if (ut[z(t)] < 2) {
      return "period " + pv(pn.tv[z(t)]) +
             " has one unit: its cross-section averages need two or more.";
    }
  }
  const int need = q + na * (p + 1) + 2 * (k + ll);
  for (int a = 0; a < nu && !dyn; ++a) {
    if (tu[z(a)] < need) {
      return fmt("unit %s has %d periods: its regression on the regressors, the averages and the "
                 "observed factors needs at least %d.",
                 pn.ulev[z(a)].c_str(), tu[z(a)], need);
    }
  }
  std::vector<double> wv(z(nu), 1.0 / nu), d(z(nt) * z(q), 1.0), yy(z(nu) * z(nt), NAN),
      xx(z(nu) * z(nt) * z(k), 0.0);
  std::vector<char> seen(z(nt), 0), useen(z(nu), 0);
  if (!w.empty()) {
    double s = 0.0;
    for (int i = 0; i < n; ++i) {
      const int a = pn.u[z(i)];
      if (!useen[z(a)]) {
        wv[z(a)] = w[z(i)];
        s += wv[z(a)];
        useen[z(a)] = 1;
      } else if (wv[z(a)] != w[z(i)]) {
        return "`w` must be constant within each unit.";
      }
    }
    for (double &e : wv) e /= s;
  }
  for (int i = 0; i < n; ++i) {
    const int a = pn.u[z(i)], t = pn.t[z(i)];
    yy[z(a) * z(nt) + z(t)] = pn.y[i];
    for (int j = 0; j < k; ++j) xx[(z(j) * z(nu) + z(a)) * z(nt) + z(t)] = x.v[z(j) * z(n) + z(i)];
    for (int j = 0; j < fm.p; ++j) {
      const double v = fm.v[z(j) * z(n) + z(i)];
      double &e = d[z(cons + j) * z(nt) + z(t)];
      if (seen[z(t)] && e != v) return "`fac` must be constant within each period.";
      e = v;
    }
    seen[z(t)] = 1;
  }
  dset ds;
  ds.nu = nu;
  ds.nt = nt;
  ds.k = k;
  ds.ll = ll;
  ds.p = p;
  ds.q = q;
  ds.yy = &yy;
  ds.xx = &xx;
  ds.wv = &wv;
  ds.d = &d;
  ds.av = av;
  const int kd = ll + k;
  int bad = -1, sh = -1, off = 0;
  const char *where = "";
  // the error of a fit, "" when it succeeded
  auto err = [&](int rc) -> std::string {
    if (rc == CLAB_EINTR || clab_intrd(0) != 0) {
      code = CLAB_EINTR;
      return "stop";
    }
    if (rc == CLAB_OK) return "";
    if (sh >= 0) {
      return fmt("unit %s has %d usable periods%s: its regression on the regressors, the lags of "
                 "y, the averages with their lags and the observed factors needs at least %d.",
                 pn.ulev[z(sh)].c_str(), bad, where, need);
    }
    if (rc == CLAB_EDOM && bad >= 0) {
      return "unit " + pn.ulev[z(bad)] +
             ": its regressors are collinear once the cross-section averages are taken out.";
    }
    if (rc == CLAB_EDOM && bad <= -2) {
      return "after the recursive mean adjustment period " + pv(pn.tv[z(-2 - bad + off)]) +
             " has fewer than two units.";
    }
    if (rc == CLAB_EDOM) return "the pooled matrix of the regressors is singular.";
    return "the work arrays of the estimation would be too large.";
  };
  dfit f0;
  std::string m = err(drun(ds, 0, nt, need, pool, f0, bad, sh));
  if (!m.empty()) return code == CLAB_EINTR ? "" : m;
  dfit fc = f0;  // the estimate reported
  if (co.bias == 1) {
    // the half-panel jackknife (eq 32), each half a panel of its own
    const int h = nt / 2;
    dfit fa, fb;
    where = " in the first half of the periods";
    m = err(drun(ds, 0, h, need, 0, fa, bad, sh));
    if (!m.empty()) return code == CLAB_EINTR ? "" : m;
    where = " in the second half of the periods";
    m = err(drun(ds, h, nt, need, 0, fb, bad, sh));
    if (!m.empty()) return code == CLAB_EINTR ? "" : m;
    for (size_t e = 0; e < fc.bi.size(); ++e) {
      fc.bi[e] = 2 * f0.bi[e] - (fa.bi[e] + fb.bi[e]) / 2;
    }
    mgv(fc.bi, kd, nu, fc.bmg, fc.vmg);
  } else if (co.bias == 2) {
    // recursive mean adjustment of y and x, then the same estimator
    std::vector<double> y2, x2, d2(z(nt - 1) * z(q));
    rma(yy, nu, nt, 1, y2, yy);
    rma(xx, nu, nt, k, x2, yy);
    for (int c = 0; c < q; ++c) {
      for (int t = 1; t < nt; ++t) d2[z(c) * z(nt - 1) + z(t - 1)] = d[z(c) * z(nt) + z(t)];
    }
    dset d2s = ds;
    d2s.nt = nt - 1;
    d2s.yy = &y2;
    d2s.xx = &x2;
    d2s.d = &d2;
    off = 1;
    where = " after the recursive mean adjustment";
    m = err(drun(d2s, 0, nt - 1, need, 0, fc, bad, sh));
    if (!m.empty()) return code == CLAB_EINTR ? "" : m;
  }
  std::vector<std::string> nm;
  for (int l = 1; l <= ll; ++l) nm.push_back("L" + std::to_string(l) + ".y");
  nm.insert(nm.end(), x.nm.begin(), x.nm.end());
  ctab(o.mg, nm, fc.bmg.data(), fc.vmg.data(), kd, level, 0);
  o.vmg = fc.vmg;
  o.bias = co.bias != 0;
  if (o.bias) {
    ctab(o.mg0, nm, f0.bmg.data(), f0.vmg.data(), kd, level, 0);
    o.vmg0 = f0.vmg;
  }
  o.pool = pool != 0;
  o.hb[0] = o.hb[1] = o.hb[2] = NAN;
  if (pool) {
    if (co.pvar && !bswv(yy, xx, d, q, av, nu, nt, k, f0.bp, f0.vp)) {
      if (clab_intrd(0) != 0) {
        code = CLAB_EINTR;
        return "";
      }
      return "`pvar = \"bsw\"`: the columns of the averages and observed factors, or the "
             "regressors left by them, are collinear.";
    }
    ctab(o.pooled, nm, f0.bp.data(), f0.vp.data(), kd, level, 0);
    o.vp = f0.vp;
    o.hb[0] = f0.hb[0];
    o.hb[1] = f0.hb[1];
    o.hb[2] = std::isnan(f0.hb[0]) ? NAN : clab_pchisq(f0.hb[0], f0.hb[1], 0);
  }
  o.nf[0] = o.nf[1] = NAN;
  o.nf[2] = na;
  if (n == nu * nt) nfac(yy, xx, wv, nu, nt, av, o.nf);
  o.bi.assign(z(nu) * z(kd), 0.0);
  for (int a = 0; a < nu; ++a) {
    for (int j = 0; j < kd; ++j) o.bi[z(j) * z(nu) + z(a)] = fc.bi[z(a) * z(kd) + z(j)];
  }
  o.nm = nm;
  o.w = wv;
  o.units = nu;
  o.periods = nt;
  o.k = kd;
  o.h = q + na * (p + 1);
  o.n = fc.n;
  o.lags = p;
  o.ylag = ll;
  return "";
}

}  // namespace dpnd
