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

std::string cfit(const pnl &pn, const xm &x, const xm &fm, int cons, const std::vector<double> &w,
                 bool useall, const std::vector<std::string> &avg, double level, cres &o,
                 int &code) {
  code = CLAB_OK;
  const int n = pn.n, nu = pn.nu, nt = pn.nt, k = x.p;
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
  const double cells = static_cast<double>(nu) * nt;
  if (cells * (k + 1) > cap) {
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
  const int need = q + na + 2 * k;
  for (int a = 0; a < nu; ++a) {
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
  std::vector<double> bi(z(k) * z(nu)), bmg(z(k)), vmg(z(k) * z(k)), bp(z(k)), vp(z(k) * z(k));
  int bad = -1;
  const int rc = clab_ccepnl(nu, nt, k, yy.data(), xx.data(), wv.data(), d.data(), q, av.data(),
                             bi.data(), bmg.data(), vmg.data(), bp.data(), vp.data(), &bad);
  if (rc == CLAB_EINTR || clab_intrd(0) != 0) {
    code = CLAB_EINTR;
    return "";
  }
  if (rc == CLAB_EDOM && bad >= 0) {
    return "unit " + pn.ulev[z(bad)] +
           ": its regressors are collinear once the cross-section averages are taken out.";
  }
  if (rc == CLAB_EDOM) return "the pooled matrix of the regressors is singular.";
  if (rc != CLAB_OK) return "the work arrays of the estimation would be too large.";
  ctab(o.mg, x.nm, bmg.data(), vmg.data(), k, level, 0);
  ctab(o.pooled, x.nm, bp.data(), vp.data(), k, level, 0);
  o.bi.assign(z(nu) * z(k), 0.0);
  for (int a = 0; a < nu; ++a) {
    for (int j = 0; j < k; ++j) o.bi[z(j) * z(nu) + z(a)] = bi[z(a) * z(k) + z(j)];
  }
  o.vmg = vmg;
  o.vp = vp;
  o.w = wv;
  o.units = nu;
  o.periods = nt;
  o.k = k;
  o.h = q + na;
  o.n = n;
  return "";
}

}  // namespace dpnd
