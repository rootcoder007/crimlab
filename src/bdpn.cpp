// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the dynamic and heterogeneous panels: difference
// and system GMM (gmmpnl), Wooldridge's test of serial correlation
// (wsctst) and the common correlated effects estimators (ccepnl). Each
// checks its arguments, builds the kernel's arrays (dpn.cpp) and builds
// the R result. Every array whose size is a product of the panel's
// dimensions is counted and refused past a bound before it is made, and
// the arrays are made inside try blocks, so that a failed allocation ends
// in an R error, never in an exception through .Call. Nothing calls back
// into the R evaluator.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "crimlab.h"
#include "kutil.h"
#include "pevr.h"
#include "rchk.h"
#include "spin.h"

extern "C" double Rf_pnorm5(double, double, double, int, int);
extern "C" double Rf_pchisq(double, double, int, int);
extern "C" double Rf_pf(double, double, double, int, int);

using clab::fail;
using pevr::cols;
using pevr::nlist;
using pevr::z;

namespace {

// the largest array made, in elements (2^26), and the most instruments
// or regressors (their weight matrices are square in that number)
constexpr double cap = 67108864.0;
constexpr int mmax = 2000;

// An interrupted kernel unwinds to the entry point (glue.cpp).
void intr(int code) {
  if (code == CLAB_EINTR || clab_intrd(0) != 0) throw clab::Intr{};
}

std::string pv(double v) {
  char b[32];
  std::snprintf(b, sizeof b, "%.15g", v);
  return b;
}

// The rows of a panel by unit, in increasing period: unit a's rows are
// ord[uo[a]] .. ord[uo[a + 1] - 1].
struct idx {
  const pevr::pan &pn;
  std::vector<int> ord, uo;
  explicit idx(const pevr::pan &p) : pn(p), ord(z(p.n)), uo(z(p.nu) + 1, 0) {
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

// What gmmpnl() was asked for.
struct spec {
  int ylag = 1, sys = 1, fod = 0, coll = 0, cons = 1, tdum = 1, hm = 0;
  std::vector<int> gv, lo, hi;  // the GMM-style variables (-1: y) and lags
  std::vector<int> iv;          // the columns of x that are IV-style
};

// The kernel's arrays and what the result reports of the design.
struct dsg {
  int n = 0, k = 0, m = 0, ng = 0, nd = 0, nobs = 0, omin = 0, omax = 0;
  std::vector<int> go, eq, dgo, dt, sub, grp;
  std::vector<double> y, x, z, h, yd, xd;
  std::vector<std::string> xn, zn, drop, dum;
  std::string bad;
  double big = 0;
};

// The level rows, the transformed rows and their weights on the level rows
// (Roodman 2009, eq 21 and sec. 3.1).
struct rows {
  std::vector<int> lr, lu, ls, luo;     // level rows: row, unit, period
  std::vector<int> tu, tl, wo, wq;      // transformed rows: unit, label
  std::vector<double> wv;               // weights wv on level rows wq
  std::vector<int> tuo;                 // transformed rows of each unit
};

rows mkrows(const pevr::pan &pn, const idx &ix, const spec &sp, double &big) {
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
bool gval(const pevr::pan &pn, const idx &ix, const pevr::mat &xm, int v, int a, int s,
          double &out) {
  const int r = ix.at(a, s);
  if (r < 0) return false;
  out = v < 0 ? pn.y[r] : xm.v[z(v) * z(pn.n) + z(r)];
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
int build(const pevr::pan &pn, const pevr::mat &xm, const pevr::mat &zx, const spec &sp,
          dsg &d) {
  const idx ix(pn);
  const rows w = mkrows(pn, ix, sp, d.big);
  if (d.big > 0) return 3;
  const int nt = pn.nt, p = xm.p, nlr = static_cast<int>(w.lr.size());
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
    for (int c = 0; c < p; ++c) lv[z(sp.ylag + c) * z(nlr) + z(q)] = xm.v[z(c) * z(pn.n) + z(r)];
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
  nm.insert(nm.end(), xm.nm.begin(), xm.nm.end());
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
            if (!gval(pn, ix, xm, g.v, a, lab - l, v) || v == 0.0) continue;
            const int c = g.base + (sp.coll ? l - g.a : lab * g.nl + l - g.a);
            if (pass == 0) col[z(c)] = 0;
            else d.z[z(col[z(c)]) * z(ns) + z(i)] = v;
          }
        } else {
          const int q = sl[z(i)], a = w.lu[z(q)], s = w.ls[z(q)];
          double v1 = 0.0, v0 = 0.0;
          if (!gval(pn, ix, xm, g.v, a, s - g.a + 1, v1) ||
              !gval(pn, ix, xm, g.v, a, s - g.a, v0) || v1 == v0) {
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
      const std::string vn = g.v < 0 ? "y" : xm.nm[z(g.v)];
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

// Columns of x named by s (character), each once.
std::vector<int> cnames(SEXP s, const char *an, const pevr::mat &xm) {
  if (TYPEOF(s) != STRSXP) fail("`%s` must be names of columns of `x`.", an);
  std::vector<int> o;
  for (R_xlen_t i = 0; i < XLENGTH(s); ++i) {
    if (STRING_ELT(s, i) == NA_STRING) fail("`%s` must not hold NA.", an);
    const std::string v = Rf_translateCharUTF8(STRING_ELT(s, i));
    const auto it = std::find(xm.nm.begin(), xm.nm.end(), v);
    if (it == xm.nm.end()) fail("`%s` names %s, which is not a column of `x`.", an, v.c_str());
    const int j = static_cast<int>(it - xm.nm.begin());
    if (std::find(o.begin(), o.end(), j) != o.end()) fail("`%s` names %s twice.", an, v.c_str());
    o.push_back(j);
  }
  return o;
}

// Lag limits, one value or one per GMM-style variable: whole numbers from
// lo, or Inf (big) when inf.
std::vector<int> lags(SEXP s, const char *an, int nv, int lo, bool inf) {
  if (!clab::isnum(s) || (XLENGTH(s) != 1 && XLENGTH(s) != nv)) {
    fail("`%s` must be one number or %d, one per GMM-style variable.", an, nv);
  }
  SEXP d = PROTECT(clab::dbl(s));
  std::vector<int> o;
  for (int i = 0; i < nv; ++i) {
    const double v = REAL(d)[XLENGTH(s) == 1 ? 0 : i];
    if (inf && v == R_PosInf) {
      o.push_back(1 << 30);
      continue;
    }
    if (!std::isfinite(v) || v != std::floor(v) || v < lo || v > 1e6) {
      fail("`%s` must hold whole numbers from %d%s.", an, lo, inf ? ", or Inf" : "");
    }
    o.push_back(static_cast<int>(v));
  }
  UNPROTECT(1);
  return o;
}

// A coefficient table with normal intervals.
SEXP cframe(const std::vector<std::string> &nm, const double *b, const double *v, int k,
            double lv) {
  const double q = Rf_qnorm5(1 - (1 - lv) / 2, 0.0, 1.0, 1, 0);
  cols co(7);
  SEXP ct = co.str("term", k);
  double *ce = co.dbl("est", k), *cs = co.dbl("se", k), *cz = co.dbl("z", k),
         *cp = co.dbl("p", k), *cl = co.dbl("lo", k), *ch = co.dbl("hi", k);
  for (int j = 0; j < k; ++j) {
    SET_STRING_ELT(ct, j, Rf_mkCharCE(nm[z(j)].c_str(), CE_UTF8));
    ce[j] = b[j];
    cs[j] = std::sqrt(v[z(j) * z(k) + z(j)]);
    cz[j] = ce[j] / cs[j];
    cp[j] = pevr::na(2 * Rf_pnorm5(-std::fabs(cz[j]), 0.0, 1.0, 1, 0));
    cl[j] = ce[j] - q * cs[j];
    ch[j] = ce[j] + q * cs[j];
  }
  return co.frame();
}

SEXP strs(const std::vector<std::string> &v) {
  SEXP s = PROTECT(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(v.size())));
  for (size_t i = 0; i < v.size(); ++i) {
    SET_STRING_ELT(s, static_cast<R_xlen_t>(i), Rf_mkCharCE(v[i].c_str(), CE_UTF8));
  }
  UNPROTECT(1);
  return s;
}

SEXP rmat(const std::vector<double> &v, int r, int c) {
  SEXP m = PROTECT(Rf_allocMatrix(REALSXP, r, c));
  std::copy(v.begin(), v.end(), REAL(m));
  UNPROTECT(1);
  return m;
}

// The options of gmmpnl() into sp.
void gopts(spec &sp, const pevr::mat &xm, SEXP ylag, SEXP gmm, SEXP iv, SEXP lmin, SEXP lmax,
           SEXP coll, SEXP type, SEXP trans, SEXP hmat, SEXP cons, SEXP tdum) {
  static const char *const tyo[] = {"system", "diff"};
  static const char *const tro[] = {"fd", "fod"};
  static const char *const hmo[] = {"ab", "dpd", "id"};
  sp.ylag = clab::whole(ylag, "ylag", 0);
  if (sp.ylag > 1000) fail("`ylag` must be at most 1000.");
  if (sp.ylag == 0 && xm.p == 0) fail("the model needs `ylag` above 0 or columns in `x`.");
  sp.sys = csp::pick(type, "type", tyo, 2) == 0;
  sp.fod = csp::pick(trans, "trans", tro, 2);
  sp.hm = csp::pick(hmat, "hmat", hmo, 3);
  sp.coll = csp::flag(coll, "collapse");
  sp.cons = csp::flag(cons, "cons");
  sp.tdum = csp::flag(tdum, "tdum");
  const std::vector<int> gx = Rf_isNull(gmm) ? std::vector<int>() : cnames(gmm, "gmm", xm);
  if (sp.ylag > 0) sp.gv.push_back(-1);
  sp.gv.insert(sp.gv.end(), gx.begin(), gx.end());
  if (Rf_isNull(iv)) {
    for (int j = 0; j < xm.p; ++j) {
      if (std::find(gx.begin(), gx.end(), j) == gx.end()) sp.iv.push_back(j);
    }
  } else {
    sp.iv = cnames(iv, "iv", xm);
    for (int j : sp.iv) {
      if (std::find(gx.begin(), gx.end(), j) != gx.end()) {
        fail("`iv` and `gmm` both name %s.", xm.nm[z(j)].c_str());
      }
    }
  }
  const int nv = static_cast<int>(sp.gv.size());
  if (nv > 0) {
    sp.lo = lags(lmin, "lagmin", nv, 1, false);
    sp.hi = lags(lmax, "lagmax", nv, 1, true);
    for (int i = 0; i < nv; ++i) {
      if (sp.hi[z(i)] < sp.lo[z(i)]) fail("`lagmax` must be at least `lagmin`.");
    }
  }
}

// The tests of a gmmpnl() result.
SEXP gtests(const double *st, const double *ar, int nar, int m, int k, int msub, bool sys) {
  const int nr = nar + 2 + (sys ? 2 : 0);
  cols c(4);
  SEXP tn = c.str("test", nr);
  double *s = c.dbl("stat", nr), *df = c.dbl("df", nr), *pp = c.dbl("p", nr);
  for (int l = 0; l < nar; ++l) {
    SET_STRING_ELT(tn, l, Rf_mkChar(("AR(" + std::to_string(l + 1) + ")").c_str()));
    s[l] = pevr::na(ar[l]);
    df[l] = NA_REAL;
    pp[l] = std::isnan(ar[l]) ? NA_REAL : 2 * Rf_pnorm5(-std::fabs(ar[l]), 0.0, 1.0, 1, 0);
  }
  const double hx = st[2], dif = st[1] - st[2];
  const char *nm[] = {"Sargan", "Hansen", "Hansen without level GMM", "Difference-in-Hansen"};
  const double sv[] = {st[0], st[1], hx, dif};
  const double dv[] = {static_cast<double>(m - k), static_cast<double>(m - k),
                       static_cast<double>(m - msub - k), static_cast<double>(msub)};
  for (int j = 0; j < nr - nar; ++j) {
    const int r = nar + j;
    SET_STRING_ELT(tn, r, Rf_mkChar(nm[j]));
    s[r] = pevr::na(sv[j]);
    df[r] = dv[j];
    pp[r] = std::isnan(sv[j]) || dv[j] < 1 ? NA_REAL : Rf_pchisq(sv[j], dv[j], 0, 0);
  }
  return c.frame();
}

}  // namespace

extern "C" {

//' @srrstats {G2.6} gmmpnl(), wsctst() and ccepnl() read the outcome, units
//'   and periods as plain vectors of one length; regressors, instruments
//'   and observed factors as a numeric matrix with one row per row.
//' @noRd
SEXP C_clab_r_gmmpnl(SEXP y, SEXP x, SEXP unit, SEXP time, SEXP ylag, SEXP gmm, SEXP iv,
                     SEXP zz, SEXP lmin, SEXP lmax, SEXP coll, SEXP type, SEXP trans, SEXP steps,
                     SEXP wc, SEXP hmat, SEXP cons, SEXP tdum, SEXP ar, SEXP level) {
  clab::who = "gmmpnl";
  int np = 0;
  pevr::pan pn = pevr::panel(y, unit, time, clab::ninf, &np);
  const pevr::mat xm = Rf_isNull(x) ? pevr::mat() : pevr::readm(x, "x", "x", pn.n, clab::ninf);
  const pevr::mat zx = Rf_isNull(zz) ? pevr::mat() : pevr::readm(zz, "z", "z", pn.n, clab::ninf);
  spec sp;
  gopts(sp, xm, ylag, gmm, iv, lmin, lmax, coll, type, trans, hmat, cons, tdum);
  const int ns = clab::whole(steps, "steps", 1);
  if (ns > 2) fail("`steps` must be 1 or 2.");
  const bool wcr = csp::flag(wc, "wc");
  const int nar = clab::whole(ar, "ar", 0);
  if (nar > 100) fail("`ar` must be at most 100.");
  const double lv = clab::one(level, "level", 0.0, 1.0);
  dsg d;
  int code = 0;
  try {
    code = build(pn, xm, zx, sp, d);
  } catch (const clab::Intr &) {
    throw;
  } catch (...) {
    code = 6;
  }
  if (code == 4) throw clab::Intr{};
  if (code == 1) fail("no unit has the rows the transformed equation needs (two usable periods).");
  if (code == 2) {
    fail("the regressor %s is constant over time or collinear with the others in the stacked "
         "equations.",
         d.bad.c_str());
  }
  if (code == 3) {
    fail("the instruments or the arrays would be too large (%.0f; at most %d instruments): "
         "collapse them or narrow the lags.",
         d.big, mmax);
  }
  if (code == 5) fail("the regressors would number %.0f, more than %d.", d.big, mmax);
  if (code == 6) fail("the design could not be built: out of memory.");
  if (d.m < d.k) fail("the %d instruments are fewer than the %d coefficients.", d.m, d.k);
  int msub = 0;
  for (int f : d.sub) msub += f;
  std::vector<double> b(z(d.k)), v(z(d.k) * z(d.k)), st(5), av(z(nar));
  const int rc = clab_gmmpnl(d.n, d.k, d.m, d.ng, d.go.data(), d.y.data(), d.x.data(), d.z.data(),
                             d.h.data(), d.eq.data(), d.nd, d.dgo.data(), d.dt.data(),
                             d.yd.data(), d.xd.data(), sp.sys && msub ? d.sub.data() : nullptr,
                             ns, wcr ? 1 : 0, nar, b.data(), v.data(), st.data(), av.data());
  intr(rc);
  if (rc == CLAB_EDOM) {
    fail("the instruments do not identify the coefficients, or the moments of the residuals "
         "are all 0.");
  }
  if (rc != CLAB_OK) fail("the estimation could not be run: out of memory.");
  std::vector<std::string> note;
  if (d.m > d.ng) {
    note.push_back(std::to_string(d.m) + " instruments for " + std::to_string(d.ng) +
                   " groups: more instruments than groups weakens the Hansen test "
                   "(Roodman 2009, sec. 2.6).");
  }
  if (st[4] > 0) {
    note.push_back("a weight matrix was singular: a generalized inverse was used, and the "
                   "difference-in-Hansen statistic may be negative.");
  }
  SEXP coef = PROTECT(cframe(d.xn, b.data(), v.data(), d.k, lv));
  SEXP te = PROTECT(gtests(st.data(), av.data(), nar, d.m, d.k, msub, sp.sys));
  cols fs(9);
  double *fn = fs.dbl("n", 1), *fg = fs.dbl("groups", 1), *fi = fs.dbl("inst", 1),
         *f0 = fs.dbl("omin", 1), *fa = fs.dbl("oavg", 1), *f1 = fs.dbl("omax", 1),
         *fst = fs.dbl("steps", 1), *fk = fs.dbl("k", 1), *fsg = fs.dbl("sigma2", 1);
  fn[0] = d.nobs;
  fg[0] = d.ng;
  fi[0] = d.m;
  f0[0] = d.omin;
  fa[0] = static_cast<double>(d.nobs) / d.ng;
  f1[0] = d.omax;
  fst[0] = ns;
  fk[0] = d.k;
  fsg[0] = pevr::na(st[3]);
  SEXP fit = PROTECT(fs.frame());
  SEXP vm = PROTECT(rmat(v, d.k, d.k));
  SEXP dn = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dn, 0, strs(d.xn));
  SET_VECTOR_ELT(dn, 1, strs(d.xn));
  Rf_setAttrib(vm, R_DimNamesSymbol, dn);
  SEXP in = PROTECT(strs(d.zn)), dr = PROTECT(strs(d.drop)), du = PROTECT(strs(d.dum));
  SEXP nt = PROTECT(strs(note));
  static const char *const tyo[] = {"system", "diff"};
  static const char *const tro[] = {"fd", "fod"};
  SEXP ty = PROTECT(Rf_mkString(tyo[sp.sys ? 0 : 1])), tr = PROTECT(Rf_mkString(tro[sp.fod]));
  SEXP lvs = PROTECT(Rf_ScalarReal(lv));
  np += 12;
  SEXP out = PROTECT(nlist({{"coef", coef}, {"tests", te}, {"fit", fit}, {"vcov", vm},
                            {"inst", in}, {"dropped", dr}, {"dummies", du}, {"note", nt},
                            {"type", ty}, {"trans", tr}, {"level", lvs}}));
  Rf_classgets(out, PROTECT(Rf_mkString("gmmpnl")));
  UNPROTECT(np + 2);
  return out;
}

SEXP C_clab_r_wsctst(SEXP y, SEXP x, SEXP unit, SEXP time, SEXP cons, SEXP ssc) {
  clab::who = "wsctst";
  static const char *const sso[] = {"cr1", "none"};
  int np = 0;
  pevr::pan pn = pevr::panel(y, unit, time, clab::ninf, &np);
  const pevr::mat xm = Rf_isNull(x) ? pevr::mat() : pevr::readm(x, "x", "x", pn.n, clab::ninf);
  const int c = csp::flag(cons, "cons");
  const int raw = csp::pick(ssc, "ssc", sso, 2);
  if (xm.p > mmax) fail("`x` must have at most %d columns.", mmax);
  std::vector<double> b(z(xm.p)), vb(z(xm.p) * z(xm.p)), o(6);
  const int rc = clab_wsctst(pn.n, xm.p, pn.y, xm.v.data(), pn.u.data(), pn.t.data(), c, raw,
                             b.data(), vb.data(), o.data());
  intr(rc);
  if (rc == CLAB_EDOM) {
    fail("too few consecutive first differences, or units with them, for the test, or the "
         "differences of `x` are collinear.");
  }
  if (rc != CLAB_OK) fail("the test could not be run: out of memory.");
  const double st[] = {o[2], 1, o[3] - 1, Rf_pf(o[2], 1, o[3] - 1, 0, 0), o[0], o[1], o[4],
                       o[3], o[5]};
  static const char *const sn[] = {"F", "df1", "df2", "p", "coef", "se", "pairs", "clusters",
                                   "ndiff"};
  SEXP sv = PROTECT(clab::nvec(9, sn, st));
  SEXP co = R_NilValue;
  if (xm.p > 0) {
    cols cf(3);
    SEXP tn = cf.str("term", xm.p);
    double *ce = cf.dbl("est", xm.p), *cs = cf.dbl("se", xm.p);
    for (int j = 0; j < xm.p; ++j) {
      SET_STRING_ELT(tn, j, Rf_mkCharCE(xm.nm[z(j)].c_str(), CE_UTF8));
      ce[j] = b[z(j)];
      cs[j] = std::sqrt(vb[z(j) * z(xm.p) + z(j)]);
    }
    co = cf.frame();
  }
  PROTECT(co);
  SEXP sv2 = PROTECT(Rf_mkString(sso[raw]));
  np += 3;
  SEXP out = PROTECT(nlist({{"stat", sv}, {"coef", co}, {"ssc", sv2}}));
  Rf_classgets(out, PROTECT(Rf_mkString("wsctst")));
  UNPROTECT(np + 2);
  return out;
}

SEXP C_clab_r_ccepnl(SEXP y, SEXP x, SEXP unit, SEXP time, SEXP fac, SEXP cons, SEXP w,
                     SEXP level) {
  clab::who = "ccepnl";
  int np = 0;
  pevr::pan pn = pevr::panel(y, unit, time, clab::ninf, &np);
  const int n = pn.n, nu = pn.nu, nt = pn.nt;
  if (static_cast<double>(nu) * nt != n) {
    fail("the panel must be balanced: every unit in every period (Pesaran 2006 sets the "
         "estimators out for such a panel).");
  }
  const pevr::mat xm = pevr::readm(x, "x", "x", n, clab::ninf);
  const int k = xm.p;
  if (k > mmax) fail("`x` must have at most %d columns.", mmax);
  const pevr::mat fm = Rf_isNull(fac) ? pevr::mat() : pevr::readm(fac, "fac", "f", n, clab::ninf);
  const int cn = csp::flag(cons, "cons"), q = cn + fm.p;
  if (q + 1 + 2 * k > nt) {
    fail("the %d periods must be at least the %d averages and observed factors and twice the "
         "%d regressors together.",
         nt, q + 1 + k, k);
  }
  const double lv = clab::one(level, "level", 0.0, 1.0);
  std::vector<double> wv(z(nu), 1.0 / nu), d(z(nt) * z(q), 1.0), yy(z(n)), xx(z(n) * z(k));
  std::vector<char> seen(z(nt), 0), useen(z(nu), 0);
  if (!Rf_isNull(w)) {
    clab::vec(w, "w");
    SEXP ws = PROTECT(clab::num(w, "w", n, 0.0, true));
    ++np;
    double s = 0.0;
    for (int i = 0; i < n; ++i) {
      const int a = pn.u[z(i)];
      if (!useen[z(a)]) {
        wv[z(a)] = REAL(ws)[i];
        s += wv[z(a)];
        useen[z(a)] = 1;
      } else if (wv[z(a)] != REAL(ws)[i]) {
        fail("`w` must be constant within each unit.");
      }
    }
    for (double &e : wv) e /= s;
  }
  for (int i = 0; i < n; ++i) {
    const int a = pn.u[z(i)], t = pn.t[z(i)];
    yy[z(a) * z(nt) + z(t)] = pn.y[i];
    for (int j = 0; j < k; ++j) {
      xx[(z(j) * z(nu) + z(a)) * z(nt) + z(t)] = xm.v[z(j) * z(n) + z(i)];
    }
    for (int j = 0; j < fm.p; ++j) {
      const double v = fm.v[z(j) * z(n) + z(i)];
      double &e = d[z(cn + j) * z(nt) + z(t)];
      if (seen[z(t)] && e != v) fail("`fac` must be constant within each period.");
      e = v;
    }
    seen[z(t)] = 1;
  }
  std::vector<double> bi(z(k) * z(nu)), bmg(z(k)), vmg(z(k) * z(k)), bp(z(k)), vp(z(k) * z(k));
  const int rc = clab_ccepnl(nu, nt, k, yy.data(), xx.data(), wv.data(), d.data(), q, bi.data(),
                             bmg.data(), vmg.data(), bp.data(), vp.data());
  intr(rc);
  if (rc == CLAB_EDOM) {
    fail("a unit's regressors are collinear once the cross-section averages are taken out, "
         "or the pooled matrix is singular.");
  }
  if (rc != CLAB_OK) fail("the estimation could not be run: out of memory.");
  SEXP mg = PROTECT(cframe(xm.nm, bmg.data(), vmg.data(), k, lv));
  SEXP pl = PROTECT(cframe(xm.nm, bp.data(), vp.data(), k, lv));
  std::vector<double> bt(z(nu) * z(k));
  for (int a = 0; a < nu; ++a) {
    for (int j = 0; j < k; ++j) bt[z(j) * z(nu) + z(a)] = bi[z(a) * z(k) + z(j)];
  }
  SEXP bm = PROTECT(rmat(bt, nu, k));
  SEXP dn = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dn, 0, pn.ulev);
  SET_VECTOR_ELT(dn, 1, strs(xm.nm));
  Rf_setAttrib(bm, R_DimNamesSymbol, dn);
  SEXP vm1 = PROTECT(rmat(vmg, k, k)), vm2 = PROTECT(rmat(vp, k, k));
  cols fs(4);
  double *fu = fs.dbl("units", 1), *ft = fs.dbl("periods", 1), *fk = fs.dbl("k", 1),
         *fh = fs.dbl("h", 1);
  fu[0] = nu;
  ft[0] = nt;
  fk[0] = k;
  fh[0] = q + 1 + k;
  SEXP fit = PROTECT(fs.frame());
  SEXP wr = PROTECT(Rf_allocVector(REALSXP, nu));
  std::copy(wv.begin(), wv.end(), REAL(wr));
  Rf_setAttrib(wr, R_NamesSymbol, pn.ulev);
  SEXP lvs = PROTECT(Rf_ScalarReal(lv));
  np += 9;
  SEXP out = PROTECT(nlist({{"mg", mg}, {"pooled", pl}, {"bi", bm}, {"vmg", vm1}, {"vp", vm2},
                            {"fit", fit}, {"w", wr}, {"level", lvs}}));
  Rf_classgets(out, PROTECT(Rf_mkString("ccepnl")));
  UNPROTECT(np + 2);
  return out;
}

}  // extern "C"
