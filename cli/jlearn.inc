// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Learners (src/rlearn.cpp): the elastic-net path, logistic regression and
// the random forest.

// rl::mat(): a numeric matrix (an array of rows) with at least one row,
// every value finite; column-major into v.
bool lmat(Ctx &c, const Val &x, const char *nm, int &n, int &p, std::vector<double> &v) {
  int dim[2] = {-1, -1};
  if (x.t != Val::ARR || x.a.empty() || !shape(x, 2, dim, 0) || dim[0] < 1) {
    return c.fail("`%s` must be a numeric matrix with at least one row.", nm);
  }
  n = dim[0];
  p = dim[1];
  v.assign(static_cast<std::size_t>(n) * static_cast<std::size_t>(p), 0.0);
  bool bad = false;
  fill(x, 2, dim, 0, 0, 1, v, bad, true);
  for (double e : v) {
    if (!std::isfinite(e)) return c.fail("`%s` must be finite numbers with no missing values.", nm);
  }
  return true;
}
// rl::folds(): whole numbers 1..K, n of them, every fold used, K >= 2.
bool lfolds(Ctx &c, const Val &f, const char *nm, int n, std::vector<int> &o, int &k) {
  Dv v;
  CK(num(c, f, nm, v, n, 1.0));
  double mx = 1.0;
  for (double e : v.v) {
    if (e != std::floor(e)) return c.fail("`%s` must be whole numbers.", nm);
    mx = std::fmax(mx, e);
  }
  if (mx < 2.0 || mx > n) return c.fail("`%s` must number at least two folds.", nm);
  k = static_cast<int>(mx);
  o.assign(v.n(), 0);
  std::vector<char> seen(static_cast<std::size_t>(k), 0);
  for (std::size_t i = 0; i < v.n(); ++i) {
    o[i] = static_cast<int>(v[i]) - 1;
    seen[static_cast<std::size_t>(o[i])] = 1;
  }
  for (int j = 0; j < k; ++j) {
    if (!seen[static_cast<std::size_t>(j)]) return c.fail("`%s` leaves fold %d empty.", nm, j + 1);
  }
  return true;
}
// rl::ychk(): 0 or 1 when bin, not constant.
bool ychk(Ctx &c, const Dv &y, bool bin, const char *nm) {
  bool one = false, zero = false, diff = false;
  for (double e : y.v) {
    if (bin && e != 0.0 && e != 1.0) return c.fail("`%s` must be 0 or 1 for a binary outcome.", nm);
    one = one || e == 1.0;
    zero = zero || e == 0.0;
    diff = diff || e != y.v[0];
  }
  if (!diff || (bin && !(one && zero))) return c.fail("`%s` is constant.", nm);
  return true;
}
bool lflag(Ctx &c, const char *nm, bool def, int &b) {
  std::vector<int> v = {def ? 1 : 0};
  if (c.has(nm)) CK(lgl(c, c[nm], nm, 1, v));
  b = v[0];
  return true;
}

bool enpath(Ctx &c, Val &out) {
  int n = 0, p = 0;
  std::vector<double> xv;
  CK(lmat(c, c["x"], "x", n, p, xv));
  if (n < 2 || p < 1) return c.fail("`x` must have at least two rows and one column.");
  Dv yv;
  CK(num(c, c["y"], "y", yv, n));
  static const char *const fams[] = {"gaussian", "binomial"};
  int f = 0;
  CK(pickd(c, "fam", fams, 2, 0, f));
  CK(ychk(c, yv, f == 1, "y"));
  double a = 0.0;
  CK(p01d(c, "alpha", 1, a));
  int nl = 0;
  double r = 0.0;
  std::vector<double> lam;
  if (!c.has("lambda")) {
    CK(wholed(c, "nlam", 100, 1, nl));
    r = n > p ? 1e-4 : 1e-2;
    if (c.has("lmr")) CK(one(c, c["lmr"], "lmr", 0, 1, r));
    lam.assign(static_cast<std::size_t>(nl), 0.0);
  } else {
    Dv l;
    CK(num(c, c["lambda"], "lambda", l, -1, 0.0));
    lam = l.v;
    nl = static_cast<int>(lam.size());
    for (std::size_t k = 1; k < lam.size(); ++k) {
      if (lam[k] > lam[k - 1]) return c.fail("`lambda` must be decreasing.");
    }
  }
  int nf = 0;
  std::vector<int> fd;
  if (c.has("fold")) CK(lfolds(c, c["fold"], "fold", n, fd, nf));
  Dv th;
  Val d1;
  CK(num(c, arg(c, "tol", d1, 1e-7), "tol", th, 1, 0.0, true));
  int mi = 0, sc = 0;
  CK(wholed(c, "maxit", 1e5, 1, mi));
  CK(lflag(c, "scr", true, sc));
  const double dp = p, dn = n;
  const std::size_t L = static_cast<std::size_t>(nl), P = static_cast<std::size_t>(p);
  std::vector<double> a0(L), beta(P * L), cv(2 * L);
  std::vector<double> w(static_cast<std::size_t>(dp * dp + 6 * dp + 5 * dn + dn * dp + nl * (dp + 1) + nf * (nl + 1.0)) + 1);
  std::vector<int> iw(static_cast<std::size_t>(2 * dn + dp) + 1);
  int best[2], nk = nl;
  const int rc = clab_enpath(xv.data(), n, p, yv.p(), f, a, lam.data(), nl, r, th[0], mi, sc, fd.empty() ? nullptr : fd.data(),
                             nf, a0.data(), beta.data(), cv.data(), best, &nk, w.data(), iw.data());
  if (rc == CLAB_EDOM) {
    return c.fail("no fit: no convergence within `maxit` passes, or `y` constant once a fold is held out.");
  }
  CK(kchk(c, rc));
  out = vobj();
  put(out, "lambda", vauto(std::vector<double>(lam.begin(), lam.begin() + nk)));
  put(out, "a0", vauto(std::vector<double>(a0.begin(), a0.begin() + nk)));
  put(out, "beta", vmat(beta, p, nk));
  if (!fd.empty()) {
    put(out, "cvm", vauto(std::vector<double>(cv.begin(), cv.begin() + nk)));
    put(out, "cvsd", vauto(std::vector<double>(cv.begin() + nk, cv.begin() + 2 * nk)));
    put(out, "lmin", vnum(lam[static_cast<std::size_t>(best[0])]));
    put(out, "l1se", vnum(lam[static_cast<std::size_t>(best[1])]));
  } else {
    for (const char *k : {"cvm", "cvsd", "lmin", "l1se"}) put(out, k, Val());
  }
  return true;
}

bool lgtreg(Ctx &c, Val &out) {
  int n = 0, p = 0;
  std::vector<double> xv;
  CK(lmat(c, c["x"], "x", n, p, xv));
  Dv yv;
  CK(num(c, c["y"], "y", yv, n));
  CK(ychk(c, yv, true, "y"));
  if (n <= p + 1) return c.fail("`x` needs more rows than columns plus one.");
  Dv th;
  Val d1;
  CK(num(c, arg(c, "tol", d1, 1e-10), "tol", th, 1, 0.0, true));
  int mi = 0;
  CK(wholed(c, "maxit", 50, 1, mi));
  const std::size_t q = static_cast<std::size_t>(p) + 1;
  std::vector<double> w(2 * q * q + 4 * q), cf(q), se(q);
  std::vector<int> iw(static_cast<std::size_t>(n));
  double dev = 0.0;
  const int rc = clab_lgtreg(xv.data(), n, p, yv.p(), th[0], mi, cf.data(), se.data(), &dev, w.data(), iw.data());
  if (rc == CLAB_EDOM) {
    return c.fail("no fit: collinear columns, or separated outcomes (no convergence within `maxit` steps).");
  }
  CK(kchk(c, rc));
  std::vector<std::string> nms = {"(Intercept)"};
  for (int j = 0; j < p; ++j) nms.push_back("x" + std::to_string(j + 1));
  out = vobj();
  put(out, "coef", nvec(nms, cf.data()));
  put(out, "se", nvec(nms, se.data()));
  put(out, "dev", vnum(dev));
  return true;
}

bool rfor(Ctx &c, Val &out) {
  int n = 0, p = 0, nt = 0, pt = 0;
  std::vector<double> xv, tv;
  CK(lmat(c, c["x"], "x", n, p, xv));
  if (n < 2) return c.fail("`x` must have at least two rows.");
  Dv yv;
  CK(num(c, c["y"], "y", yv, n));
  CK(lmat(c, c.has("newx") ? c["newx"] : c["x"], "newx", nt, pt, tv));
  if (pt != p) return c.fail("`newx` must have the %d columns of `x`.", p);
  int ntr = 0, mt = 0, mn = 0;
  CK(wholed(c, "ntree", 500, 1, ntr));
  mt = static_cast<int>(std::fmin(std::ceil(std::sqrt(p) + 20.0), p));
  if (c.has("mtry")) CK(whole(c, c["mtry"], "mtry", 1, mt));
  if (mt > p) return c.fail("`mtry` must be at most the %d columns of `x`.", p);
  CK(wholed(c, "minnode", 5, 1, mn));
  Dv fr, ib, sd;
  Val d1, d2, d3;
  CK(num(c, arg(c, "frac", d1, 1), "frac", fr, 1, 0.0, true));
  if (fr[0] > 1.0) return c.fail("`frac` must be at most 1.");
  int rp = 0, hn = 0;
  CK(lflag(c, "repl", true, rp));
  CK(lflag(c, "hon", false, hn));
  CK(num(c, arg(c, "imb", d2, 0), "imb", ib, 1, 0.0));
  if (ib[0] >= 0.5) return c.fail("`imb` must be below 0.5.");
  CK(num(c, arg(c, "seed", d3, 1), "seed", sd, 1, 0.0));
  if (std::llround(fr[0] * n) < (hn ? 2 : 1)) return c.fail("`frac` leaves too few rows in a bag.");
  std::vector<double> pr(static_cast<std::size_t>(nt)), ob(static_cast<std::size_t>(n));
  std::vector<double> w(static_cast<std::size_t>(lrn::fwork(n, p, nt)) + 1);
  std::vector<int> iw(static_cast<std::size_t>(lrn::fiwork(n, p, nt)) + static_cast<std::size_t>(n + nt) + 1);
  CK(kchk(c, clab_rfor(xv.data(), n, p, yv.p(), tv.data(), nt, ntr, mt, mn, fr[0], rp, hn, ib[0],
                       static_cast<unsigned int>(std::fmod(sd[0], 4294967296.0)), pr.data(), ob.data(), w.data(),
                       iw.data())));
  out = vobj();
  put(out, "pred", vauto(pr));
  put(out, "oob", vauto(ob));
  return true;
}

// Double machine learning and the super learner.

const char *const dmdl[] = {"plr", "pliv", "irm", "iivm", "aipw", "tmle", "did", "didcs", "ssm"};
const char *const dnn[][5] = {{"l", "m", "g"},     {"l", "m", "r", "g"}, {"g0", "g1", "m"},
                              {"g0", "g1", "m", "r0", "r1"}, {"g0", "g1", "m"}, {"g0", "g1", "m"},
                              {"g0", "m"},         {"g00", "g01", "g10", "g11", "m"}, {"pi", "m", "g0", "g1"}};
const char *const lch[] = {"glm", "lasso", "ridge", "enet", "forest", "stack"};

struct Lset {
  double alpha;
  int nlam, nfin;
  double tol;
  int scr, hon;
  unsigned int seed;
  int ntree, mtry;
};

clab_lrn lbase(int k, const Lset &s) {
  clab_lrn r = {k == 0 ? 0 : (k == 4 ? 2 : 1), k == 1 ? 1.0 : (k == 2 ? 0.0 : s.alpha), s.nlam, s.nfin, s.seed, s.ntree,
                s.mtry, 0, s.tol, s.scr, s.hon, nullptr, 0, 0, 0};
  return r;
}

// The text values of an argument (with a default when left out).
std::vector<Val> strv(Ctx &c, const char *nm, std::initializer_list<const char *> def) {
  std::vector<Val> e;
  if (!c.has(nm)) {
    for (const char *s : def) e.push_back(vstr(s));
    return e;
  }
  const Val &v = c[nm];
  if (v.t == Val::ARR) e = v.a;
  else e.push_back(v);
  return e;
}

bool lsubs(Ctx &c, const std::vector<Val> &stk, const char *nm, const Lset &s, clab_lrn *sub, int &ns) {
  const std::size_t m = stk.size();
  bool ok = m >= 1 && m <= 5;
  for (const Val &e : stk) ok = ok && e.t == Val::STR;
  if (!ok) return c.fail("`%s` must name from one to five learners.", nm);
  int seen[5] = {0, 0, 0, 0, 0};
  for (std::size_t j = 0; j < m; ++j) {
    int k = 0;
    CK(pick(c, stk[j], nm, lch, 5, k));
    if (seen[k]++ > 0) return c.fail("`%s` names \"%s\" twice.", nm, lch[k]);
    sub[j] = lbase(k, s);
  }
  ns = static_cast<int>(m);
  return true;
}

bool learner(Ctx &c, const Val &l, const char *nm, const Lset &s, const std::vector<Val> &stk, int nstk, int disc,
             clab_lrn *sub, clab_lrn &r) {
  int k = 0;
  CK(pick(c, l.t == Val::ARR && l.a.size() == 1 ? l.a[0] : l, nm, lch, 6, k));
  if (k < 5) {
    r = lbase(k, s);
    return true;
  }
  int ns = 0;
  CK(lsubs(c, stk, "stk", s, sub, ns));
  r = {3, s.alpha, s.nlam, s.nfin, s.seed, s.ntree, s.mtry, 0, s.tol, s.scr, s.hon, sub, ns, nstk, disc};
  return true;
}

// bin(): numbers, n of them, 0 or 1 and not constant.
bool lbin(Ctx &c, const char *nm, int n, Dv &v) {
  CK(num(c, c[nm], nm, v, n));
  return ychk(c, v, true, nm);
}

bool dmlx(Ctx &c, Val &out) {
  int md = 0;
  CK(pickd(c, "model", dmdl, 9, 0, md));
  int n = 0, p = 0;
  std::vector<double> xv;
  CK(lmat(c, c["x"], "x", n, p, xv));
  if (n < 4 || p < 1) return c.fail("`x` must have at least four rows and one column.");
  const std::size_t sn = static_cast<std::size_t>(n);
  Dv sv, yv, dv, zv;
  if (md == 8) CK(lbin(c, "s", n, sv));
  else if (md == 7) CK(lbin(c, "t", n, sv));
  if (md == 8) {
    const Val &y = c["y"];
    if (!isnum(y) || vlen(y) != sn) return c.fail("`y` must be numbers, one per row of `x`.");
    yv.v = dbl(y);
    for (std::size_t i = 0; i < sn; ++i) {
      if (sv[i] == 0.0) yv.v[i] = 0.0;
      else if (!std::isfinite(yv.v[i])) return c.fail("`y` must be finite where `s` is 1.");
    }
  } else {
    CK(num(c, c["y"], "y", yv, n));
  }
  if (md >= 2) CK(lbin(c, "d", n, dv));
  else CK(num(c, c["d"], "d", dv, n));
  if (md == 1) CK(num(c, c["z"], "z", zv, n));
  else if (md == 3) CK(lbin(c, "z", n, zv));
  int sc = 0;
  if (c.has("score")) {
    static const char *const s01[] = {"po", "iv"}, *const s23[] = {"ate", "atte"}, *const s67[] = {"obs", "exp"};
    if (md == 0 || md == 1) CK(pick(c, c["score"], "score", s01, 2, sc));
    else if (md == 2) CK(pick(c, c["score"], "score", s23, 2, sc));
    else if (md == 6 || md == 7) CK(pick(c, c["score"], "score", s67, 2, sc));
    else return c.fail("`score` must be NULL for model \"%s\".", dmdl[md]);
  }
  // folds: a vector, or a matrix with one column per repetition
  const Val &fo = c["fold"];
  int fdim[2] = {-1, -1};
  const bool fm = fo.t == Val::ARR && !fo.a.empty() && fo.a[0].t == Val::ARR && shape(fo, 2, fdim, 0);
  const int nrep = fm ? fdim[1] : 1;
  if (nrep < 1 || (fm ? fdim[0] : static_cast<int>(vlen(fo))) != n) return c.fail("`fold` must have one row per row of `x`.");
  std::vector<Val> fcols;
  if (fm) {
    for (int r = 0; r < nrep; ++r) {
      Val v = varr();
      for (const Val &row : fo.a) v.a.push_back(row.a[static_cast<std::size_t>(r)]);
      fcols.push_back(v);
    }
  } else {
    fcols.push_back(fo);
  }
  for (const Val &v : fcols) {
    if (!isnum(v)) return c.fail("`fold` must be finite numbers with no missing values.");
  }
  int nf = 0;
  std::vector<int> fv(sn * static_cast<std::size_t>(nrep));
  for (int r = 0; r < nrep; ++r) {
    std::vector<int> f;
    int k = 0;
    CK(lfolds(c, fcols[static_cast<std::size_t>(r)], "fold", n, f, k));
    if (r > 0 && k != nf) return c.fail("`fold` must number the same folds in every column.");
    nf = k;
    std::copy(f.begin(), f.end(), fv.begin() + static_cast<std::ptrdiff_t>(static_cast<std::size_t>(r) * sn));
  }
  int nc = 0, nc2 = 0, nrk = nrep;
  std::vector<int> cw;
  if (c.has("clu")) {
    if (md > 1) return c.fail("`clu` is for models \"plr\" and \"pliv\" only.");
    const Val &cu = c["clu"];
    int cd[2] = {-1, -1};
    const bool two = cu.t == Val::ARR && !cu.a.empty() && cu.a[0].t == Val::ARR && shape(cu, 2, cd, 0) && cd[1] == 2;
    if (two) {
      if (cd[0] != n || nrep != 2) {
        return c.fail("two-way `clu` needs one row per row of `x` and `fold` with two columns, the folds of each kind "
                      "of cluster.");
      }
      Val c1 = varr(), c2v = varr();
      for (const Val &row : cu.a) {
        c1.a.push_back(row.a[0]);
        c2v.a.push_back(row.a[1]);
      }
      std::vector<int> a, b;
      CK(lfolds(c, c1, "clu", n, a, nc));
      CK(lfolds(c, c2v, "clu", n, b, nc2));
      cw = a;
      cw.insert(cw.end(), b.begin(), b.end());
      nrk = 1;
    } else {
      CK(lfolds(c, cu, "clu", n, cw, nc));
    }
    std::vector<int> cf(static_cast<std::size_t>(nc > nc2 ? nc : nc2));
    for (int r = 0; r < nrep; ++r) {
      const int *cr = cw.data() + (two ? static_cast<std::size_t>(r) * sn : 0);
      const int ncr = two && r == 1 ? nc2 : nc;
      for (int q = 0; q < ncr; ++q) cf[static_cast<std::size_t>(q)] = -1;
      for (std::size_t i = 0; i < sn; ++i) {
        const int f = fv[static_cast<std::size_t>(r) * sn + i];
        int &e = cf[static_cast<std::size_t>(cr[i])];
        if (e >= 0 && e != f) return c.fail("`fold` must be constant within each cluster.");
        e = f;
      }
    }
  }
  double al = 0.0;
  CK(p01d(c, "alpha", 0.5, al));
  int nl = 0, ni = 0, nt = 0, mt = 0, ks = 0;
  CK(wholed(c, "nlam", 100, 1, nl));
  CK(wholed(c, "nfin", 5, 2, ni));
  Dv sdv, tlv;
  Val d1, d2;
  CK(num(c, arg(c, "seed", d1, 1), "seed", sdv, 1, 0.0));
  const unsigned int us = static_cast<unsigned int>(std::fmod(sdv[0], 4294967296.0));
  CK(wholed(c, "ntree", 500, 1, nt));
  if (c.has("mtry")) CK(whole(c, c["mtry"], "mtry", 1, mt));
  CK(num(c, arg(c, "tol", d2, 1e-7), "tol", tlv, 1, 0.0, true));
  int sr = 0, hn = 0, dc = 0;
  CK(lflag(c, "scr", true, sr));
  CK(lflag(c, "hon", false, hn));
  CK(wholed(c, "nstk", 5, 2, ks));
  CK(lflag(c, "disc", false, dc));
  clab_lrn sy[5], sdd[5];
  const Lset ls = {al, nl, ni, tlv[0], sr, hn, us, nt, mt};
  Lset lsd = ls;
  lsd.seed = us + 1U;
  const std::vector<Val> stk = strv(c, "stk", {"glm", "lasso", "forest"});
  const Val lrn0 = c.has("lrn") ? c["lrn"] : vstr("lasso");
  clab_lrn ly{}, ld{};
  CK(learner(c, lrn0, "lrn", ls, stk, ks, dc, sy, ly));
  CK(learner(c, c.has("lrnd") ? c["lrnd"] : lrn0, "lrnd", lsd, stk, ks, dc, sdd, ld));
  double tr = 0.0;
  if (!c.has("trim")) {
    tr = md == 4 ? 0.025 : (md == 5 ? 5.0 / std::sqrt(n) / std::log(n) : (md == 6 || md == 7 ? 0.01 : 1e-12));
  } else {
    Dv t;
    CK(num(c, c["trim"], "trim", t, 1, 0.0));
    tr = t[0];
    if (tr >= 0.5) return c.fail("`trim` must be below 0.5.");
  }
  const int k = clab_dmlxk(md, sc);
  const std::size_t sk = static_cast<std::size_t>(k), sr2 = static_cast<std::size_t>(nrk);
  std::vector<double> pr(sn * sk * sr2);
  int fit = 1;
  if (c.has("pred")) {
    const Val &pd = c["pred"];
    int dim3[3] = {-1, -1, -1}, dim2[2] = {-1, -1};
    std::vector<double> pv;
    bool ok = false;
    if (nrk == 1 && shape(pd, 2, dim2, 0) && pd.t == Val::ARR && !pd.a.empty()) {
      ok = dim2[0] == n && dim2[1] == k;
      if (ok) {
        pv.assign(sn * sk, 0.0);
        bool bad = false;
        fill(pd, 2, dim2, 0, 0, 1, pv, bad, true);
      }
    } else if (pd.t == Val::ARR && !pd.a.empty() && shape(pd, 3, dim3, 0)) {
      ok = dim3[0] == n && dim3[1] == k && dim3[2] == nrk;
      if (ok) {
        pv.assign(sn * sk * sr2, 0.0);
        bool bad = false;
        fill(pd, 3, dim3, 0, 0, 1, pv, bad, true);
      }
    }
    if (!ok) {
      return c.fail("`pred` must be an array of %d rows, %d columns and one slice per column of `fold`.", n, k);
    }
    for (double e : pv) {
      if (!std::isfinite(e)) return c.fail("`pred` must be finite.");
    }
    std::copy(pv.begin(), pv.end(), pr.begin());
    fit = 0;
  }
  long niw = 0;
  const std::size_t nw = static_cast<std::size_t>(clab_dmlxsz(md, n, p, nf, nrk, &ly, &ld, &niw));
  std::vector<double> w(nw + 1), est(2 + 2 * sr2), ps(2 * sn * sr2);
  std::vector<int> iw(static_cast<std::size_t>(niw) + 1);
  const int rc = clab_dmlx(md, sc, xv.data(), n, p, yv.p(), dv.p(), zv.n() ? zv.p() : nullptr, sv.n() ? sv.p() : nullptr,
                           fv.data(), nf, nrk, cw.empty() ? nullptr : cw.data(), nc, nc2, &ly, &ld, tr, fit, pr.data(),
                           est.data(), ps.data(), w.data(), iw.data());
  if (rc == CLAB_EDOM) return c.fail("no estimate: a learner failed, a training set was empty, or the score has no root.");
  CK(kchk(c, rc));
  Frame reps;
  reps.num("coef", est.data() + 2, sr2);
  reps.num("se", est.data() + 2 + sr2, sr2);
  const int d3[3] = {n, k, nrk};
  out = vobj();
  put(out, "coef", vnum(est[0]));
  put(out, "se", vnum(est[1]));
  put(out, "reps", reps.o);
  put(out, "psia", vmat(ps.data(), n, nrk));
  put(out, "psib", vmat(ps.data() + sn * sr2, n, nrk));
  put(out, "pred", varr3(pr.data(), d3));
  return true;
}

bool suplrn(Ctx &c, Val &out) {
  int n = 0, p = 0;
  std::vector<double> xv;
  CK(lmat(c, c["x"], "x", n, p, xv));
  if (n < 4 || p < 1) return c.fail("`x` must have at least four rows and one column.");
  static const char *const fams[] = {"gaussian", "binomial"};
  int fm = 0;
  CK(pickd(c, "fam", fams, 2, 0, fm));
  Dv yv;
  CK(num(c, c["y"], "y", yv, n));
  CK(ychk(c, yv, fm == 1, "y"));
  int nt = n;
  std::vector<double> tv = xv;
  if (c.has("newx")) {
    int pt = 0;
    CK(lmat(c, c["newx"], "newx", nt, pt, tv));
    if (pt != p) return c.fail("`newx` must have the %d columns of `x`.", p);
  }
  Dv sdv;
  Val d1;
  CK(num(c, arg(c, "seed", d1, 1), "seed", sdv, 1, 0.0));
  const unsigned int us = static_cast<unsigned int>(std::fmod(sdv[0], 4294967296.0));
  double al = 0.0;
  int nl = 0, ni = 0, ntr = 0, mt = 0, dc = 0;
  CK(p01d(c, "alpha", 0.5, al));
  CK(wholed(c, "nlam", 100, 1, nl));
  CK(wholed(c, "nfin", 5, 2, ni));
  CK(wholed(c, "ntree", 500, 1, ntr));
  if (c.has("mtry")) CK(whole(c, c["mtry"], "mtry", 1, mt));
  const Lset ls = {al, nl, ni, 1e-7, 1, 0, us, ntr, mt};
  CK(lflag(c, "disc", false, dc));
  clab_lrn sub[5];
  int ns = 0;
  const std::vector<Val> lr = strv(c, "lrns", {"glm", "lasso", "forest"});
  CK(lsubs(c, lr, "lrns", ls, sub, ns));
  int k = 0;
  std::vector<int> fv;
  const std::size_t sn = static_cast<std::size_t>(n);
  if (!c.has("fold")) {
    CK(wholed(c, "nfold", 5, 2, k));
    if (k > n) return c.fail("`nfold` must be at most the %d rows of `x`.", n);
    fv.assign(sn, 0);
    std::vector<int> wk(sn);
    CK(kchk(c, clab_kfold(n, k, us, fv.data(), wk.data())));
  } else {
    CK(lfolds(c, c["fold"], "fold", n, fv, k));
  }
  const clab_lrn l = {3, ls.alpha, ls.nlam, ls.nfin, us, ls.ntree, ls.mtry, 0, ls.tol, ls.scr, ls.hon, sub, ns, k, dc};
  long niw = 0;
  const std::size_t nw = static_cast<std::size_t>(clab_stksz(n, p, nt, &l, &niw));
  std::vector<double> w(nw + 1), oof(sn * static_cast<std::size_t>(ns)), risk(static_cast<std::size_t>(ns)),
      wt(static_cast<std::size_t>(ns)), pr(static_cast<std::size_t>(nt));
  std::vector<int> iw(static_cast<std::size_t>(niw) + 1);
  const int rc = clab_stack(xv.data(), n, p, yv.p(), fm, &l, fv.data(), tv.data(), nt, oof.data(), risk.data(), wt.data(),
                            pr.data(), w.data(), iw.data());
  if (rc == CLAB_EDOM) {
    return c.fail("a learner failed on a fold (a fold too small for its inner folds, or a design it cannot fit).");
  }
  CK(kchk(c, rc));
  std::vector<std::string> nm;
  for (int j = 0; j < ns; ++j) nm.push_back(lr[static_cast<std::size_t>(j)].s);
  std::vector<int> fo(sn);
  for (std::size_t i = 0; i < sn; ++i) fo[i] = fv[i] + 1;
  out = vobj();
  put(out, "wt", nvec(nm, wt.data()));
  put(out, "risk", nvec(nm, risk.data()));
  put(out, "fold", unbox(vints(fo.data(), sn)));
  put(out, "oof", vmat(oof, n, ns));
  put(out, "pred", vauto(pr));
  return true;
}
