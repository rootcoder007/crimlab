// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The explanation of a prediction function (src/xpl.cpp): R/xpl.R in C++.
// Outside R the function is the learners named by f, fitted to x and y by
// suplrn() with the settings of lctl, and called once on all the rows a
// method asks for.
// ponytail: the rows reach suplrn() as JSON values, some 100 bytes a
// cell; a design near the 2^26 values allowed needs gigabytes, call
// clab_stack() on the doubles if that size matters.

namespace xp {

inline std::size_t sz(double v) { return static_cast<std::size_t>(v); }

const char *const set[] = {"fam", "nfold", "fold", "disc", "alpha", "nlam", "nfin", "ntree", "mtry", "seed"};

// .xpck(): x a numeric matrix, finite (column-major into x, its names x1,
// x2, ...); f the names of learners; lctl their settings; y given.
inline bool ck(Ctx &c, int &n, int &p, std::vector<double> &x, std::vector<std::string> &nm) {
  const Val &o = c["x"];
  int dim[2] = {-1, -1};
  if (o.t != Val::ARR || o.a.empty() || o.a[0].t != Val::ARR || !shape(o, 2, dim, 0) || dim[0] < 1 || dim[1] < 1) {
    return ckstop(c, "`x` must be a numeric matrix with at least one row and one column");
  }
  n = dim[0];
  p = dim[1];
  x.assign(sz(n) * sz(p), 0.0);
  bool bad = false;
  fill(o, 2, dim, 0, 0, 1, x, bad, false);
  if (bad) return ckstop(c, "`x` must be finite numbers with no missing values");
  nm.clear();
  for (int j = 0; j < p; ++j) nm.push_back("x" + std::to_string(j + 1));
  const Val &f = c["f"];
  bool lrn = f.t == Val::STR || (f.t == Val::ARR && !f.a.empty());
  for (std::size_t i = 0; f.t == Val::ARR && i < f.a.size(); ++i) lrn = lrn && f.a[i].t == Val::STR;
  if (!lrn) return ckstop(c, "`f` must be a function or the names of learners of suplrn()");
  const Val &l = c["lctl"];
  bool ok = l.t == Val::NUL || (l.t == Val::ARR && l.a.empty()) || l.t == Val::OBJ;
  for (std::size_t i = 0; ok && l.t == Val::OBJ && i < l.k.size(); ++i) {
    bool in = false;
    for (const char *s : set) in = in || l.k[i] == s;
    for (std::size_t q = 0; q < i; ++q) in = in && l.k[q] != l.k[i];
    ok = in;
  }
  if (!ok) {
    return ckstop(c, "`lctl` must be a list named by settings of suplrn(): fam, nfold, fold, disc, alpha, nlam, "
                     "nfin, ntree, mtry, seed");
  }
  if (!c.has("y")) return ckstop(c, "`y` must be given to fit the learners named in `f`");
  return true;
}

// .xpf(): the predictions of the learners at the rows of m (nr x p).
inline bool pred(Ctx &c, const std::vector<double> &m, int nr, int p, std::vector<double> &out) {
  Val a = vobj();
  put(a, "x", c["x"]);
  put(a, "y", c["y"]);
  put(a, "lrns", c["f"]);
  put(a, "newx", vmat(m, nr, p));
  const Val &l = c["lctl"];
  for (std::size_t i = 0; l.t == Val::OBJ && i < l.k.size(); ++i) put(a, l.k[i].c_str(), l.a[i]);
  Ctx s{c.fn, &a, std::string(), {}};
  Val r;
  if (!suplrn(s, r)) {
    c.err = s.err;
    return false;
  }
  Dv d;
  readd(::col(r, "pred"), d);
  out = d.v;
  return true;
}

// The rows of a (na rows) over those of b (nb rows), p columns.
inline std::vector<double> rbind(const std::vector<double> &a, int na, const std::vector<double> &b, int nb, int p) {
  std::vector<double> m(sz(na + nb) * sz(p));
  for (int j = 0; j < p; ++j) {
    std::copy(a.begin() + static_cast<std::ptrdiff_t>(sz(j) * sz(na)),
              a.begin() + static_cast<std::ptrdiff_t>(sz(j + 1) * sz(na)),
              m.begin() + static_cast<std::ptrdiff_t>(sz(j) * sz(na + nb)));
    std::copy(b.begin() + static_cast<std::ptrdiff_t>(sz(j) * sz(nb)),
              b.begin() + static_cast<std::ptrdiff_t>(sz(j + 1) * sz(nb)),
              m.begin() + static_cast<std::ptrdiff_t>(sz(j) * sz(na + nb) + sz(na)));
  }
  return m;
}

// .xpbig()
inline bool big(Ctx &c, double rows, int p, const char *what) {
  if (rows * p > 67108864) {
    return ckstop(c, "`%s` asks for %.0f rows of %d columns to predict, more than 67108864 values", what, rows, p);
  }
  return true;
}

}  // namespace xp

bool xpimp(Ctx &c, Val &out) {
  int n = 0, p = 0;
  std::vector<double> x;
  std::vector<std::string> nm;
  CK(xp::ck(c, n, p, x, nm));
  Nv yv, rv, sd;
  CK(cknum(c, c["y"], "y", yv, false));
  if (yv.n() != xp::sz(n)) return ckstop(c, "`y` must have one value per row of `x`, %d", n);
  CK(cknumd(c, "nrep", 10, rv, 1, HUGE_VAL, false, true));
  int lo = 0;
  CK(ckarg(c, "loss", {"se", "ae"}, lo));
  bool rt = false;
  CK(flagd(c, "ratio", false, rt));
  CK(cknumd(c, "seed", 1, sd, 0));
  const int k = static_cast<int>(rv[0]);
  const double rows = static_cast<double>(n) * p * k;
  CK(xp::big(c, rows, p, "nrep"));
  std::vector<double> des(xp::sz(rows) * xp::sz(p));
  std::vector<int> perm(xp::sz(rows));
  CK(kchk(c, clab_xpperm(x.data(), n, p, k, static_cast<unsigned int>(std::fmod(sd[0], 4294967296.0)), des.data(),
                         perm.data())));
  std::vector<double> v;
  CK(xp::pred(c, xp::rbind(x, n, des, static_cast<int>(rows), p), n + static_cast<int>(rows), p, v));
  std::vector<double> rep(xp::sz(p) * xp::sz(k)), st(1 + 2 * xp::sz(p));
  const int rc = clab_xpimp(yv.p(), v.data(), v.data() + n, n, p, k, lo, rt ? 1 : 0, rep.data(), st.data());
  if (rc == CLAB_EDOM) return ckstop(c, "`ratio` needs a loss above 0 before permuting");
  CK(kchk(c, rc));
  std::vector<double> pm(perm.begin(), perm.end());
  for (double &e : pm) e += 1;
  const int d3[3] = {n, k, p};
  Frame f;
  f.str("term", nm);
  f.num("imp", st.data() + 1, xp::sz(p));
  f.num("sd", st.data() + 1 + p, xp::sz(p));
  out = vobj();
  put(out, "imp", f.o);
  put(out, "base", vnum(st[0]));
  put(out, "rep", vmat(rep, p, k));
  put(out, "perm", varr3(pm.data(), d3));
  put(out, "loss", vstr(lo == 0 ? "se" : "ae"));
  put(out, "ratio", vbool(rt));
  return true;
}

namespace xp {

// .xpv(): the column v of x by name or number, from 0.
inline bool vcol(Ctx &c, const std::vector<std::string> &nm, int &j) {
  const Val &o = c["v"];
  const Val &e = o.t == Val::ARR && o.a.size() == 1 ? o.a[0] : o;
  j = -1;
  if (e.t == Val::STR) {
    for (std::size_t q = 0; q < nm.size() && j < 0; ++q) {
      if (nm[q] == e.s) j = static_cast<int>(q);
    }
  } else if (e.t == Val::NUM && std::isfinite(e.n) && e.n == std::floor(e.n) && e.n >= 1 &&
             e.n <= static_cast<double>(nm.size())) {
    j = static_cast<int>(e.n) - 1;
  }
  if (j < 0) return ckstop(c, "`v` must name or number one column of `x`");
  return true;
}

// .xpg(): the grid of column j of x (n rows).
inline bool grid(Ctx &c, const std::vector<double> &x, int n, int j, std::vector<double> &g) {
  int ty = 0;
  CK(ckarg(c, "gtype", {"even", "quant"}, ty));
  Dv d;
  bool ok = true;
  if (!c.has("grid")) {
    d.v = {20};
  } else {
    ok = readd(c["grid"], d) == NUM && d.n() > 0;
    for (double e : d.v) ok = ok && std::isfinite(e);
  }
  if (ok && d.n() == 1) ok = d[0] >= 2 && d[0] == std::floor(d[0]) && d[0] <= 67108864.0;
  if (!ok) return ckstop(c, "`grid` must be a whole number of points from 2, or at least two values");
  if (d.n() > 1) {
    g = d.v;
    return true;
  }
  const int ng = static_cast<int>(d[0]);
  g.assign(sz(ng), 0.0);
  int m = 0;
  CK(kchk(c, clab_xpgrid(x.data() + sz(j) * sz(n), n, ng, ty, g.data(), &m)));
  g.resize(sz(m));
  return true;
}

}  // namespace xp

bool xppdp(Ctx &c, Val &out) {
  int n = 0, p = 0, j = 0;
  std::vector<double> x, g;
  std::vector<std::string> nm;
  CK(xp::ck(c, n, p, x, nm));
  CK(xp::vcol(c, nm, j));
  CK(xp::grid(c, x, n, j, g));
  const int ng = static_cast<int>(g.size());
  CK(xp::big(c, static_cast<double>(n) * ng, p, "grid"));
  std::vector<double> des(xp::sz(n) * xp::sz(ng) * xp::sz(p)), v, pd(xp::sz(ng));
  CK(kchk(c, clab_xpcpd(x.data(), n, p, j, g.data(), ng, des.data())));
  CK(xp::pred(c, des, n * ng, p, v));
  CK(kchk(c, clab_xppd(v.data(), n, ng, pd.data())));
  Frame f;
  f.num("grid", g);
  f.num("pd", pd);
  out = vobj();
  put(out, "pd", f.o);
  put(out, "term", vstr(nm[xp::sz(j)]));
  put(out, "n", vnum(n));
  return true;
}

namespace xp {

// .xpnew(): the rows of newx (by default the first row of x), with the p
// columns of x, finite.
inline bool newx(Ctx &c, const std::vector<double> &x, int n, int p, std::vector<double> &s, int &ns) {
  if (!c.has("newx")) {
    ns = 1;
    s.assign(sz(p), 0.0);
    for (int j = 0; j < p; ++j) s[sz(j)] = x[sz(j) * sz(n)];
    return true;
  }
  const Val &o = c["newx"];
  int dim[2] = {-1, -1};
  bool bad = false;
  if (o.t == Val::ARR && !o.a.empty() && o.a[0].t == Val::ARR && shape(o, 2, dim, 0) && dim[0] >= 1 && dim[1] == p) {
    s.assign(sz(dim[0]) * sz(p), 0.0);
    fill(o, 2, dim, 0, 0, 1, s, bad, false);
  } else {
    bad = true;
  }
  if (bad) {
    return ckstop(c, "`newx` must be a numeric matrix with the %d columns of `x`, no missing values", p);
  }
  ns = dim[0];
  return true;
}

}  // namespace xp

bool xpcp(Ctx &c, Val &out) {
  int n = 0, p = 0, j = 0, ns = 0;
  std::vector<double> x, s, g;
  std::vector<std::string> nm;
  CK(xp::ck(c, n, p, x, nm));
  CK(xp::vcol(c, nm, j));
  CK(xp::newx(c, x, n, p, s, ns));
  CK(xp::grid(c, x, n, j, g));
  const int ng = static_cast<int>(g.size());
  CK(xp::big(c, static_cast<double>(ns) * (ng + 1), p, "grid"));
  std::vector<double> des(xp::sz(ns) * xp::sz(ng) * xp::sz(p)), v;
  CK(kchk(c, clab_xpcpd(s.data(), ns, p, j, g.data(), ng, des.data())));
  CK(xp::pred(c, xp::rbind(des, ns * ng, s, ns, p), ns * (ng + 1), p, v));
  std::vector<double> pr(v.begin() + static_cast<std::ptrdiff_t>(ns * ng), v.end());
  out = vobj();
  put(out, "prof", vmat(v.data(), ns, ng));
  put(out, "grid", vauto(g));
  put(out, "pred", vauto(pr));
  put(out, "term", vstr(nm[xp::sz(j)]));
  return true;
}

bool xpale(Ctx &c, Val &out) {
  int n = 0, p = 0, j = 0;
  std::vector<double> x;
  std::vector<std::string> nm;
  CK(xp::ck(c, n, p, x, nm));
  CK(xp::vcol(c, nm, j));
  Nv kb;
  CK(cknumd(c, "nbin", 10, kb, 1, HUGE_VAL, false, true));
  if (kb[0] > 67108864.0) return ckstop(c, "`nbin` must be at most 67108864");
  CK(xp::big(c, 2.0 * n, p, "x"));
  const int K = static_cast<int>(kb[0]);
  std::vector<double> z(xp::sz(K) + 1), des(2 * xp::sz(n) * xp::sz(p)), v;
  std::vector<int> bin(xp::sz(n));
  int nz = 0;
  const int rc = clab_xpaled(x.data(), n, p, j, K, z.data(), &nz, bin.data(), des.data());
  if (rc == CLAB_EDOM) return ckstop(c, "`v` names a column with one value: it has no local effects");
  CK(kchk(c, rc));
  z.resize(xp::sz(nz));
  CK(xp::pred(c, des, 2 * n, p, v));
  const int nb = nz - 1;
  std::vector<double> eff(xp::sz(nb)), cnt(xp::sz(nb)), ale(xp::sz(nz));
  CK(kchk(c, clab_xpale(v.data(), n, bin.data(), nb, eff.data(), cnt.data(), ale.data())));
  Frame a, b;
  a.num("z", z);
  a.num("ale", ale);
  b.num("lo", z.data(), xp::sz(nb));
  b.num("hi", z.data() + 1, xp::sz(nb));
  b.num("n", cnt);
  b.num("eff", eff);
  out = vobj();
  put(out, "ale", a.o);
  put(out, "bins", b.o);
  put(out, "term", vstr(nm[xp::sz(j)]));
  return true;
}

bool xpshap(Ctx &c, Val &out) {
  int n = 0, p = 0, ns = 0;
  std::vector<double> x, s;
  std::vector<std::string> nm;
  CK(xp::ck(c, n, p, x, nm));
  CK(xp::newx(c, x, n, p, s, ns));
  Nv mv, sd;
  CK(cknumd(c, "nsamp", 200, mv, 1, HUGE_VAL, false, true));
  CK(cknumd(c, "seed", 1, sd, 0));
  const double rows = 2.0 * ns * p * mv[0];
  CK(xp::big(c, rows + ns + n, p, "nsamp"));
  const int m = static_cast<int>(mv[0]), q = static_cast<int>(rows);
  std::vector<double> des(xp::sz(rows) * xp::sz(p)), v;
  CK(kchk(c, clab_xpshd(x.data(), n, s.data(), ns, p, m, static_cast<unsigned int>(std::fmod(sd[0], 4294967296.0)),
                        des.data())));
  CK(xp::pred(c, xp::rbind(xp::rbind(des, q, s, ns, p), q + ns, x, n, p), q + ns + n, p, v));
  std::vector<double> phi(xp::sz(ns) * xp::sz(p)), se(phi.size());
  CK(kchk(c, clab_xpsh(v.data(), ns, p, m, phi.data(), se.data())));
  double b = 0;
  CK(kchk(c, clab_xppd(v.data() + q + ns, n, 1, &b)));
  std::vector<double> pr(v.begin() + q, v.begin() + q + ns);
  out = vobj();
  put(out, "phi", vmat(phi, ns, p));
  put(out, "se", vmat(se, ns, p));
  put(out, "pred", vauto(pr));
  put(out, "base", vnum(b));
  put(out, "nsamp", vnum(m));
  return true;
}
