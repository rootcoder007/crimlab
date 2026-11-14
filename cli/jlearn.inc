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
  int n, p;
  std::vector<double> xv;
  CK(lmat(c, c["x"], "x", n, p, xv));
  if (n < 2 || p < 1) return c.fail("`x` must have at least two rows and one column.");
  Dv yv;
  CK(num(c, c["y"], "y", yv, n));
  static const char *const fams[] = {"gaussian", "binomial"};
  int f;
  CK(pickd(c, "fam", fams, 2, 0, f));
  CK(ychk(c, yv, f == 1, "y"));
  double a;
  CK(p01d(c, "alpha", 1, a));
  int nl;
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
  int mi, sc;
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
  int n, p;
  std::vector<double> xv;
  CK(lmat(c, c["x"], "x", n, p, xv));
  Dv yv;
  CK(num(c, c["y"], "y", yv, n));
  CK(ychk(c, yv, true, "y"));
  if (n <= p + 1) return c.fail("`x` needs more rows than columns plus one.");
  Dv th;
  Val d1;
  CK(num(c, arg(c, "tol", d1, 1e-10), "tol", th, 1, 0.0, true));
  int mi;
  CK(wholed(c, "maxit", 50, 1, mi));
  const std::size_t q = static_cast<std::size_t>(p) + 1;
  std::vector<double> w(2 * q * q + 4 * q), cf(q), se(q);
  std::vector<int> iw(static_cast<std::size_t>(n));
  double dev;
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
  int n, p, nt, pt;
  std::vector<double> xv, tv;
  CK(lmat(c, c["x"], "x", n, p, xv));
  if (n < 2) return c.fail("`x` must have at least two rows.");
  Dv yv;
  CK(num(c, c["y"], "y", yv, n));
  CK(lmat(c, c.has("newx") ? c["newx"] : c["x"], "newx", nt, pt, tv));
  if (pt != p) return c.fail("`newx` must have the %d columns of `x`.", p);
  int ntr, mt, mn;
  CK(wholed(c, "ntree", 500, 1, ntr));
  mt = static_cast<int>(std::fmin(std::ceil(std::sqrt(p) + 20.0), p));
  if (c.has("mtry")) CK(whole(c, c["mtry"], "mtry", 1, mt));
  if (mt > p) return c.fail("`mtry` must be at most the %d columns of `x`.", p);
  CK(wholed(c, "minnode", 5, 1, mn));
  Dv fr, ib, sd;
  Val d1, d2, d3;
  CK(num(c, arg(c, "frac", d1, 1), "frac", fr, 1, 0.0, true));
  if (fr[0] > 1.0) return c.fail("`frac` must be at most 1.");
  int rp, hn;
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
