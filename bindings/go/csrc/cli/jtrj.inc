// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Group-based trajectory models for counts and for several outcomes.

// shares(): ng positive numbers that add up to 1.
bool shares(Ctx &c, Dv &v) {
  CK(num(c, c["prob"], "prob", v, -1, 0.0, true));
  if (std::fabs(nsum(v.p(), v.n()) - 1) > 1e-9) return c.fail("`prob` must sum to 1.");
  return true;
}

// The ZIP trajectory model's inputs.
struct Zin {
  std::vector<double> y, x, z, off, par;
  Dv prob;
  int n = 0, nt = 0, p = 0, q = 0, ng = 0;
  bool hoff = false;
};
bool zread(Ctx &c, Zin &a) {
  CK(nmat(c, c["y"], "y", a.n, a.nt, a.y, true));
  int r = 0, k = 0;
  CK(nmat(c, c["x"], "x", r, a.p, a.x));
  if (r != a.nt) return c.fail("`x` must have one row per period (column of `y`).");
  if (c.has("z")) {
    CK(nmat(c, c["z"], "z", r, a.q, a.z));
    if (r != a.nt) return c.fail("`z` must have one row per period (column of `y`).");
  }
  a.hoff = c.has("off");
  if (a.hoff) {
    CK(nmat(c, c["off"], "off", r, k, a.off));
    if (r != a.n || k != a.nt) return c.fail("`off` must have the dimensions of `y`.");
  }
  CK(shares(c, a.prob));
  a.ng = static_cast<int>(a.prob.n());
  CK(nmat(c, c["par"], "par", r, k, a.par));
  if (r != a.p + a.q || k != a.ng) {
    return c.fail("`par` must have one row per column of `x` and `z` and one column per group.");
  }
  return true;
}

bool zgll(Ctx &c, Val &out) {
  Zin a;
  CK(zread(c, a));
  std::vector<double> post(static_cast<std::size_t>(a.n) * a.ng), wk(2 * static_cast<std::size_t>(a.nt));
  double ll = 0;
  if (clab_zgll(a.y.data(), a.n, a.nt, a.hoff ? a.off.data() : nullptr, a.x.data(), a.p, a.z.data(), a.q, a.ng,
                a.par.data(), a.prob.p(), post.data(), &ll, wk.data()) != CLAB_OK) {
    return c.fail("`y` must hold whole counts from 0 (NA or a negative value is missing).");
  }
  out = vobj();
  put(out, "ll", vnum(ll));
  put(out, "post", vmat(post, a.n, a.ng));
  return true;
}

bool zgem(Ctx &c, Val &out) {
  Zin a;
  CK(zread(c, a));
  int mi = 0;
  double tl = 0.0;
  CK(wholed(c, "maxit", 500, 1, mi));
  CK(oned(c, "tol", 1e-8, 0.0, Inf, tl));
  const std::size_t pq = static_cast<std::size_t>(a.p + a.q);
  std::vector<double> post(static_cast<std::size_t>(a.n) * a.ng), wk(2 * pq * pq + 3 * pq + 7 * static_cast<std::size_t>(a.nt));
  std::vector<double> pm = a.par, pr = a.prob.v;
  double ll = 0;
  int it = 0;
  if (clab_zgem(a.y.data(), a.n, a.nt, a.hoff ? a.off.data() : nullptr, a.x.data(), a.p, a.z.data(), a.q, a.ng,
                pm.data(), pr.data(), post.data(), &ll, mi, tl, &it, wk.data()) != CLAB_OK) {
    return c.fail("`y` must hold whole counts from 0 (NA or a negative value is missing).");
  }
  out = vobj();
  put(out, "par", vmat(pm, static_cast<int>(pq), a.ng));
  put(out, "prob", vauto(pr));
  put(out, "post", vmat(post, a.n, a.ng));
  put(out, "ll", vnum(ll));
  put(out, "iter", vnum(it));
  return true;
}

// The multivariate trajectory model's inputs.
struct Min {
  std::vector<double> y, x, b, s;
  Dv prob;
  int dy[3], dx[2], db[3], ds[3], ng = 0;
};
bool mread(Ctx &c, Min &a) {
  CK(narr(c, c["y"], "y", 3, a.dy, a.y, true));
  CK(narr(c, c["x"], "x", 2, a.dx, a.x));
  CK(shares(c, a.prob));
  a.ng = static_cast<int>(a.prob.n());
  CK(narr(c, c["beta"], "beta", 3, a.db, a.b));
  CK(narr(c, c["sig"], "sig", 3, a.ds, a.s));
  const int nt = a.dy[1], r = a.dy[2], p = a.dx[1];
  if (a.dx[0] != nt || a.db[0] != p || a.db[1] != r || a.db[2] != a.ng || a.ds[0] != r || a.ds[1] != r ||
      a.ds[2] != a.ng) {
    return c.fail("the dimensions must agree: `y` n x nt x r, `x` nt x p, `beta` p x r x groups, `sig` r x r x "
                  "groups, `prob` one per group.");
  }
  return true;
}

bool mvtll(Ctx &c, Val &out) {
  Min a;
  CK(mread(c, a));
  const int n = a.dy[0], nt = a.dy[1], r = a.dy[2], p = a.dx[1];
  std::vector<double> post(static_cast<std::size_t>(n) * a.ng), wk(static_cast<std::size_t>(r) * r + r);
  double ll = 0;
  if (clab_mvtll(a.y.data(), n, nt, r, a.x.data(), p, a.ng, a.b.data(), a.s.data(), a.prob.p(), post.data(), &ll,
                 wk.data()) != CLAB_OK) {
    return c.fail("a covariance in `sig` is not positive definite.");
  }
  out = vobj();
  put(out, "ll", vnum(ll));
  put(out, "post", vmat(post, n, a.ng));
  return true;
}

bool mvtem(Ctx &c, Val &out) {
  Min a;
  CK(mread(c, a));
  const int n = a.dy[0], nt = a.dy[1], r = a.dy[2], p = a.dx[1];
  int mi = 0;
  double tl = 0.0;
  CK(wholed(c, "maxit", 500, 1, mi));
  CK(oned(c, "tol", 1e-8, 0.0, Inf, tl));
  std::vector<double> post(static_cast<std::size_t>(n) * a.ng),
      wk(static_cast<std::size_t>(p) * p + static_cast<std::size_t>(p) * r + static_cast<std::size_t>(r) * r + r);
  std::vector<double> pr = a.prob.v;
  double ll = 0;
  int it = 0;
  if (clab_mvtem(a.y.data(), n, nt, r, a.x.data(), p, a.ng, a.b.data(), a.s.data(), pr.data(), post.data(), &ll, mi,
                 tl, &it, wk.data()) != CLAB_OK) {
    return c.fail("a weighted design is singular, a residual variance is 0, or a covariance is not positive "
                  "definite.");
  }
  out = vobj();
  put(out, "beta", varr3(a.b.data(), a.db));
  put(out, "sig", varr3(a.s.data(), a.ds));
  put(out, "prob", vauto(pr));
  put(out, "post", vmat(post, n, a.ng));
  put(out, "ll", vnum(ll));
  put(out, "iter", vnum(it));
  return true;
}
