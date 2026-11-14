// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Networks (src/bind.cpp): degree, betweenness, closeness, eigenvector,
// cores, components, Louvain communities, modularity, brokerage, the
// Cheeger bound and the co-offending graph.

// nlab(): labels (a vector) with no missing values.
bool nlab(Ctx &c, const Val &x, const char *nm, std::vector<std::string> &s) {
  int dim[2] = {-1, -1};
  if (x.t == Val::ARR && shape(x, 2, dim, 0) && dim[0] > 1 && dim[1] > 1) return notvec(c, nm);
  if (x.t == Val::NUL || !alabs(x, s)) return c.fail("`%s` must be labels (a vector) with no missing values.", nm);
  return true;
}

// Codes 0, 1, ... by first appearance; first[k] the position of code k's
// first appearance.
std::vector<int> fcode(const std::vector<std::string> &v, int &nd, std::vector<std::size_t> *first) {
  std::vector<int> cd(v.size());
  std::vector<std::size_t> ix(v.size());
  for (std::size_t i = 0; i < v.size(); ++i) ix[i] = i;
  std::stable_sort(ix.begin(), ix.end(), [&v](std::size_t a, std::size_t b) { return v[a] < v[b]; });
  // the first of each run of equal labels is its first appearance
  std::vector<std::size_t> head(v.size());
  for (std::size_t q = 0; q < ix.size(); ++q) {
    head[ix[q]] = q > 0 && v[ix[q]] == v[ix[q - 1]] ? head[ix[q - 1]] : ix[q];
  }
  int k = 0;
  if (first) first->clear();
  for (std::size_t i = 0; i < v.size(); ++i) {
    if (head[i] == i) {
      if (first) first->push_back(i);
      cd[i] = k++;
    } else {
      cd[i] = cd[head[i]];
    }
  }
  nd = k;
  return cd;
}

// A graph read from its labels: n nodes (labels lab), m edges a[k] - b[k].
struct Net {
  std::vector<std::string> lab;
  int n = 0, m = 0;
  std::vector<int> a, b;
};

bool net(Ctx &c, Net &g) {
  std::vector<std::string> f, t;
  CK(nlab(c, c["from"], "from", f));
  CK(nlab(c, c["to"], "to", t));
  if (t.size() != f.size()) {
    return c.fail("`to` must have the length of `from`, %s, not %s.", rfmt(static_cast<double>(f.size())).c_str(),
                  rfmt(static_cast<double>(t.size())).c_str());
  }
  if (2.0 * f.size() > INT_MAX) return c.fail("`from` is too long for crimlab: more than %d values.", INT_MAX);
  const std::size_t m = f.size();
  g.m = static_cast<int>(m);
  g.a.assign(m, 0);
  g.b.assign(m, 0);
  if (!c.has("nodes")) {
    std::vector<std::string> v(2 * m);
    for (std::size_t k = 0; k < m; ++k) {
      v[2 * k] = f[k];
      v[2 * k + 1] = t[k];
    }
    std::vector<std::size_t> fst;
    const std::vector<int> cd = fcode(v, g.n, &fst);
    for (std::size_t k = 0; k < m; ++k) {
      g.a[k] = cd[2 * k];
      g.b[k] = cd[2 * k + 1];
    }
    for (std::size_t i : fst) g.lab.push_back(v[i]);
  } else {
    std::vector<std::string> nd;
    CK(nlab(c, c["nodes"], "nodes", nd));
    if (nd.size() > static_cast<std::size_t>(INT_MAX)) return c.fail("`nodes` is too long for crimlab.");
    g.n = static_cast<int>(nd.size());
    std::vector<std::size_t> ix(nd.size());
    for (std::size_t i = 0; i < nd.size(); ++i) ix[i] = i;
    std::stable_sort(ix.begin(), ix.end(), [&nd](std::size_t p, std::size_t q) { return nd[p] < nd[q]; });
    // the first repeat in input order
    std::size_t rep = nd.size();
    for (std::size_t q = 1; q < ix.size(); ++q) {
      if (nd[ix[q]] == nd[ix[q - 1]]) rep = std::min(rep, ix[q]);
    }
    if (rep < nd.size()) return c.fail("`nodes` must not repeat a label: %s.", nd[rep].c_str());
    for (int e = 0; e < 2; ++e) {
      const std::vector<std::string> &x = e ? t : f;
      std::vector<int> &o = e ? g.b : g.a;
      for (std::size_t k = 0; k < m; ++k) {
        const auto it = std::lower_bound(ix.begin(), ix.end(), x[k],
                                         [&nd](std::size_t p, const std::string &s) { return nd[p] < s; });
        if (it == ix.end() || nd[*it] != x[k]) {
          return c.fail("`%s` names a node not in `nodes`: %s.", e ? "to" : "from", x[k].c_str());
        }
        o[k] = static_cast<int>(*it);
      }
    }
    g.lab = nd;
  }
  if (g.n == 0) return c.fail("the graph has no nodes: give edges or `nodes`.");
  return true;
}

// The edge weights: none, or one positive number per edge.
bool nwts(Ctx &c, const Net &g, std::vector<double> &w) {
  w.clear();
  if (!c.has("w")) return true;
  Dv v;
  CK(num(c, c["w"], "w", v, g.m, 0.0, true));
  w = v.v;
  return true;
}

// A graph in compressed sparse rows; rev reverses each edge.
struct Csr {
  std::vector<int> off, adj;
  std::vector<double> wt;
};
Csr csr(const Net &g, const std::vector<double> &w, bool dir, bool rev = false) {
  const std::size_t cap = static_cast<std::size_t>(g.m) * (dir ? 1 : 2);
  Csr o;
  o.off.assign(static_cast<std::size_t>(g.n) + 1, 0);
  o.adj.assign(cap + 1, 0);
  o.wt.assign(cap + 1, 0.0);
  std::vector<int> iw(static_cast<std::size_t>(g.n));
  int nl = 0, nd = 0;
  clab_ntcsr(rev ? g.b.data() : g.a.data(), rev ? g.a.data() : g.b.data(), w.empty() ? nullptr : w.data(), g.m, g.n,
             dir, o.off.data(), o.adj.data(), o.wt.data(), iw.data(), &nl, &nd);
  return o;
}

std::vector<int> ivec(double n) { return std::vector<int>(static_cast<std::size_t>(n) + 1); }
std::vector<double> dvec(double n) { return std::vector<double>(static_cast<std::size_t>(n) + 1); }

// One value per node, coded by first appearance in node order: from the
// names of an object matched to the labels, or an array in node order.
bool pernode(Ctx &c, const char *nm, const Net &g, std::vector<int> &code, std::vector<std::string> &val) {
  const Val &x = c[nm];
  std::vector<std::string> s;
  CK(nlab(c, x, nm, s));
  val.assign(static_cast<std::size_t>(g.n), std::string());
  if (x.t != Val::OBJ) {
    if (s.size() != static_cast<std::size_t>(g.n)) {
      return c.fail("`%s` must have one value per node (%d) or names.", nm, g.n);
    }
    val = s;
  } else {
    for (std::size_t i = 0; i < x.k.size(); ++i) {
      for (std::size_t j = 0; j < i; ++j) {
        if (x.k[i] == x.k[j]) return c.fail("`%s` must not repeat a name: %s.", nm, x.k[i].c_str());
      }
    }
    for (std::size_t i = 0; i < static_cast<std::size_t>(g.n); ++i) {
      const auto it = std::find(x.k.begin(), x.k.end(), g.lab[i]);
      if (it == x.k.end()) return c.fail("`%s` has no value for node %s.", nm, g.lab[i].c_str());
      val[i] = s[static_cast<std::size_t>(it - x.k.begin())];
    }
  }
  int nd = 0;
  code = fcode(val, nd, nullptr);
  return true;
}

// The node column, then the columns given.
Frame nframe(const Net &g) {
  Frame f;
  f.str("node", g.lab);
  return f;
}
std::vector<int> plus1(const int *x, int n) {
  std::vector<int> v(x, x + n);
  for (int &e : v) ++e;
  return v;
}

bool ntdeg(Ctx &c, Val &out) {
  Net g;
  CK(net(c, g));
  std::vector<double> w;
  CK(nwts(c, g, w));
  bool d = false;
  CK(flagd(c, "dir", false, d));
  const Csr o = csr(g, w, d);
  const std::size_t n = static_cast<std::size_t>(g.n);
  std::vector<int> od(n);
  std::vector<double> os(n);
  CK(kchk(c, clab_ntdeg(g.n, o.off.data(), o.wt.data(), od.data(), os.data())));
  Frame f = nframe(g);
  if (!d) {
    f.ints("deg", od.data(), n);
    f.num("wdeg", os);
    out = f.o;
    return true;
  }
  const Csr ic = csr(g, w, true, true);
  std::vector<int> id(n), td(n);
  std::vector<double> is(n), ts(n);
  CK(kchk(c, clab_ntdeg(g.n, ic.off.data(), ic.wt.data(), id.data(), is.data())));
  for (std::size_t i = 0; i < n; ++i) {
    td[i] = id[i] + od[i];
    ts[i] = is[i] + os[i];
  }
  f.ints("in", id.data(), n);
  f.ints("out", od.data(), n);
  f.ints("deg", td.data(), n);
  f.num("win", is);
  f.num("wout", os);
  f.num("wdeg", ts);
  out = f.o;
  return true;
}

bool ntbtw(Ctx &c, Val &out) {
  Net g;
  CK(net(c, g));
  bool d = false, nz = false;
  CK(flagd(c, "dir", false, d));
  CK(flagd(c, "norm", false, nz));
  const Csr o = csr(g, {}, d);
  std::vector<double> b(static_cast<std::size_t>(g.n)), wk = dvec(2.0 * g.n);
  std::vector<int> iw = ivec(2.0 * g.n);
  CK(kchk(c, clab_ntbtw(g.n, o.off.data(), o.adj.data(), d, nz, b.data(), wk.data(), iw.data())));
  Frame f = nframe(g);
  f.num("btw", b);
  out = f.o;
  return true;
}

bool ntcls(Ctx &c, Val &out) {
  Net g;
  CK(net(c, g));
  bool d = false;
  CK(flagd(c, "dir", false, d));
  const Csr o = csr(g, {}, d);
  const std::size_t n = static_cast<std::size_t>(g.n);
  std::vector<int> r(n), iw = ivec(2.0 * g.n);
  std::vector<double> cl(n), nc(n), hr(n);
  CK(kchk(c, clab_ntcls(g.n, o.off.data(), o.adj.data(), cl.data(), nc.data(), hr.data(), r.data(), iw.data())));
  Frame f = nframe(g);
  f.ints("reach", r.data(), n);
  f.num("cls", cl);
  f.num("ncls", nc);
  f.num("hrm", hr);
  out = f.o;
  return true;
}

bool ntevc(Ctx &c, Val &out) {
  Net g;
  CK(net(c, g));
  std::vector<double> w;
  CK(nwts(c, g, w));
  double tl = 0.0;
  int mi = 0;
  CK(oned(c, "tol", 1e-13, 0, 1, tl));
  CK(wholed(c, "maxit", 100000, 1, mi));
  const Csr o = csr(g, w, false);
  std::vector<double> x(static_cast<std::size_t>(g.n)), wk = dvec(4.0 * g.n);
  std::vector<int> iw = ivec(2.0 * g.n);
  double val = 0.0;
  int it = 0, cv = 0;
  CK(kchk(c, clab_ntevc(g.n, o.off.data(), o.adj.data(), o.wt.data(), tl, mi, x.data(), &val, &it, &cv, wk.data(),
                        iw.data())));
  Frame f = nframe(g);
  f.num("evc", x);
  out = vobj();
  put(out, "cent", f.o);
  put(out, "value", vnum(val));
  put(out, "iter", vnum(it));
  put(out, "conv", vbool(cv != 0));
  return true;
}

bool ntcore(Ctx &c, Val &out) {
  Net g;
  CK(net(c, g));
  const Csr o = csr(g, {}, false);
  std::vector<int> k(static_cast<std::size_t>(g.n)), iw = ivec(3.0 * g.n);
  CK(kchk(c, clab_ntcore(g.n, o.off.data(), o.adj.data(), k.data(), iw.data())));
  Frame f = nframe(g);
  f.ints("core", k.data(), k.size());
  out = f.o;
  return true;
}

bool ntcmp(Ctx &c, Val &out) {
  Net g;
  CK(net(c, g));
  bool d = false;
  CK(flagd(c, "dir", false, d));
  const Csr o = csr(g, {}, d);
  Csr ic;
  if (d) ic = csr(g, {}, true, true);
  std::vector<int> wc(static_cast<std::size_t>(g.n)), sc(static_cast<std::size_t>(g.n)), iw = ivec(6.0 * g.n);
  int nc = 0;
  CK(kchk(c, clab_ntwcc(g.n, o.off.data(), o.adj.data(), d ? ic.off.data() : nullptr, d ? ic.adj.data() : nullptr,
                        wc.data(), &nc, iw.data())));
  CK(kchk(c, clab_ntscc(g.n, o.off.data(), o.adj.data(), sc.data(), &nc, iw.data())));
  Frame f = nframe(g);
  const std::vector<int> w1 = plus1(wc.data(), g.n), s1 = plus1(sc.data(), g.n);
  f.ints("weak", w1.data(), w1.size());
  f.ints("strong", s1.data(), s1.size());
  out = f.o;
  return true;
}

bool ntlvn(Ctx &c, Val &out) {
  Net g;
  CK(net(c, g));
  std::vector<double> w;
  CK(nwts(c, g, w));
  Val d1;
  Dv rs;
  CK(num(c, arg(c, "res", d1, 1), "res", rs, 1, 0.0, true));
  const bool use = c.has("seed");
  int sd = 0;
  if (use) CK(whole(c, c["seed"], "seed", 0, sd));
  const Csr o = csr(g, w, false);
  double len = 0.0, il = 0.0, q = 0.0;
  CK(kchk(c, clab_ntlvnw(g.n, o.off[static_cast<std::size_t>(g.n)], &len, &il)));
  std::vector<int> m(static_cast<std::size_t>(g.n)), iw = ivec(il);
  std::vector<double> wk = dvec(len);
  int lev = 0;
  CK(kchk(c, clab_ntlvn(g.n, o.off.data(), o.adj.data(), o.wt.data(), rs[0], static_cast<unsigned>(sd), use, m.data(),
                        &q, &lev, wk.data(), iw.data())));
  Frame f = nframe(g);
  const std::vector<int> m1 = plus1(m.data(), g.n);
  f.ints("comm", m1.data(), m1.size());
  out = vobj();
  put(out, "memb", f.o);
  put(out, "q", vnum(q));
  put(out, "levels", vnum(lev));
  return true;
}

bool ntmod(Ctx &c, Val &out) {
  Net g;
  CK(net(c, g));
  std::vector<int> m;
  std::vector<std::string> val;
  CK(pernode(c, "memb", g, m, val));
  std::vector<double> w;
  CK(nwts(c, g, w));
  Val d1;
  Dv rs;
  CK(num(c, arg(c, "res", d1, 1), "res", rs, 1, 0.0));
  const Csr o = csr(g, w, false);
  std::vector<double> wk = dvec(2.0 * g.n);
  double q = 0.0;
  CK(kchk(c, clab_ntmod(g.n, o.off.data(), o.adj.data(), o.wt.data(), m.data(), rs[0], &q, wk.data())));
  out = vnum(q);
  return true;
}

bool ntbrk(Ctx &c, Val &out) {
  Net g;
  CK(net(c, g));
  std::vector<int> gr;
  std::vector<std::string> val;
  CK(pernode(c, "grp", g, gr, val));
  bool d = false;
  CK(flagd(c, "dir", false, d));
  const Csr o = csr(g, {}, d);
  Csr ic;
  if (d) ic = csr(g, {}, true, true);
  const Csr &i2 = d ? ic : o;
  const std::size_t n = static_cast<std::size_t>(g.n);
  std::vector<double> cnt(5 * n);
  std::vector<int> iw = ivec(g.n);
  CK(kchk(c, clab_ntbrk(g.n, o.off.data(), o.adj.data(), i2.off.data(), i2.adj.data(), gr.data(), cnt.data(),
                        iw.data())));
  Frame f = nframe(g);
  f.str("grp", val);
  static const char *const nm[] = {"coord", "itin", "repr", "gate", "liais"};
  std::vector<double> tot(n);
  for (int r = 0; r < 5; ++r) {
    f.num(nm[r], cnt.data() + r * n, n);
    for (std::size_t i = 0; i < n; ++i) tot[i] += cnt[r * n + i];
  }
  f.num("total", tot);
  out = f.o;
  return true;
}

bool ntchg(Ctx &c, Val &out) {
  Net g;
  CK(net(c, g));
  std::vector<double> w;
  CK(nwts(c, g, w));
  double tl = 0.0;
  int mi = 0;
  CK(oned(c, "tol", 1e-10, 0, 1, tl));
  CK(wholed(c, "maxit", 100000, 1, mi));
  if (g.n < 2) return c.fail("the graph needs at least two nodes.");
  const Csr o = csr(g, w, false);
  std::vector<int> iw = ivec(g.n), set = ivec(g.n);
  int nc = 0, it = 0, cv = 0;
  CK(kchk(c, clab_ntwcc(g.n, o.off.data(), o.adj.data(), nullptr, nullptr, set.data(), &nc, iw.data())));
  if (nc > 1) return c.fail("the graph must be connected, not in %d components (see ntcmp()).", nc);
  double b[4];
  std::vector<double> wk = dvec(4.0 * g.n);
  CK(kchk(c, clab_ntchg(g.n, o.off.data(), o.adj.data(), o.wt.data(), tl, mi, b, set.data(), &it, &cv, wk.data(),
                        iw.data())));
  std::vector<std::string> s;
  for (std::size_t i = 0; i < static_cast<std::size_t>(g.n); ++i) {
    if (set[i]) s.push_back(g.lab[i]);
  }
  static const char *const bn[] = {"lambda2", "lower", "upper", "sweep"};
  out = vobj();
  put(out, "bound", nvec(4, bn, b));
  put(out, "set", unbox(vstrs(s)));
  put(out, "iter", vnum(it));
  put(out, "conv", vbool(cv != 0));
  return true;
}

bool coofg(Ctx &c, Val &out) {
  std::vector<std::string> si, so;
  CK(nlab(c, c["inc"], "inc", si));
  CK(nlab(c, c["off"], "off", so));
  const std::size_t k = si.size();
  if (so.size() != k) {
    return c.fail("`off` must have the length of `inc`, %s, not %s.", rfmt(static_cast<double>(k)).c_str(),
                  rfmt(static_cast<double>(so.size())).c_str());
  }
  if (k == 0) return c.fail("`inc` must name at least one incident.");
  if (2.0 * k > INT_MAX) return c.fail("`inc` is too long for crimlab: more than %d values.", INT_MAX);
  int ni = 0, no = 0;
  std::vector<std::size_t> fo;
  const std::vector<int> ci = fcode(si, ni, nullptr), co = fcode(so, no, &fo);
  double np = 0.0, il = 0.0;
  std::vector<int> iw = ivec(static_cast<double>(ni) + 1 + static_cast<double>(k) + no);
  CK(kchk(c, clab_coofw(ci.data(), co.data(), static_cast<int>(k), ni, no, &np, &il, iw.data())));
  const int cap = static_cast<int>(np);
  std::vector<int> pa = ivec(cap), pb = ivec(cap), iw2 = ivec(il);
  std::vector<double> pw = dvec(cap);
  int n = 0;
  CK(kchk(c, clab_coofg(ci.data(), co.data(), static_cast<int>(k), ni, no, cap, pa.data(), pb.data(), pw.data(), &n,
                        iw2.data())));
  std::vector<std::string> f, t, nd;
  for (int p = 0; p < n; ++p) {
    f.push_back(so[fo[static_cast<std::size_t>(pa[static_cast<std::size_t>(p)])]]);
    t.push_back(so[fo[static_cast<std::size_t>(pb[static_cast<std::size_t>(p)])]]);
  }
  for (int j = 0; j < no; ++j) nd.push_back(so[fo[static_cast<std::size_t>(j)]]);
  Frame e;
  e.str("from", f);
  e.str("to", t);
  e.num("w", pw.data(), static_cast<std::size_t>(n));
  out = vobj();
  put(out, "edges", e.o);
  put(out, "nodes", unbox(vstrs(nd)));
  return true;
}
