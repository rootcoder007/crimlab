// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Trace evidence: the two-level model of multivariate measurements,
// Scott-Knott grouping, lowess, the Chumbley score of tool marks and glass
// transfer.

bool tlvc(Ctx &c, Val &out) {
  int n = 0, p = 0;
  std::vector<double> xv;
  CK(nmat(c, c["x"], "x", n, p, xv));
  Lab g;
  CK(codes(c, c["item"], "item", static_cast<std::size_t>(n), g));
  if (g.g < 2 || n <= g.g) return c.fail("at least 2 items are needed, and more observations than items.");
  const std::size_t P = static_cast<std::size_t>(p);
  std::vector<double> mu(P), u(P * P), cv(P * P), wk(static_cast<std::size_t>(g.g) * P + 2 * static_cast<std::size_t>(g.g));
  if (clab_tlvc(xv.data(), n, p, g.code.data(), g.g, mu.data(), u.data(), cv.data(), wk.data()) != CLAB_OK) {
    return c.fail("the kernel refused the input.");
  }
  out = vobj();
  put(out, "mu", vauto(mu));
  put(out, "u", vmat(u, p, p));
  put(out, "c", vmat(cv, p, p));
  return true;
}

bool tlvlr(Ctx &c, Val &out) {
  const Val &bg = c["bg"];
  if ((bg.t != Val::OBJ && bg.t != Val::ARR) || bg.a.size() != 3) return c.fail("`bg` must be the list tlvc() returns.");
  Dv mu;
  CK(num(c, bg.a[0], "bg$mu", mu));
  const int p = static_cast<int>(mu.n());
  int r1 = 0, c1 = 0, r2 = 0, c2 = 0;
  std::vector<double> u, cv;
  CK(nmat(c, bg.a[1], "bg$u", r1, c1, u));
  CK(nmat(c, bg.a[2], "bg$c", r2, c2, cv));
  if (r1 != p || c1 != p || r2 != p || c2 != p) {
    return c.fail("`bg$u` and `bg$c` must be square with the length of `bg$mu`.");
  }
  Dv xa, xb;
  CK(num(c, c["xa"], "xa", xa, p));
  CK(num(c, c["xb"], "xb", xb, p));
  double ma = 0.0, mb = 0.0;
  CK(one(c, c["na"], "na", 0.0, Inf, ma));
  CK(one(c, c["nb"], "nb", 0.0, Inf, mb));
  std::vector<double> wk(static_cast<std::size_t>(p) * p + 2 * static_cast<std::size_t>(p));
  double o[2];
  if (clab_tlvlr(u.data(), cv.data(), mu.p(), p, xa.p(), ma, xb.p(), mb, o, wk.data()) != CLAB_OK) {
    return c.fail("a covariance matrix of the model is not positive definite.");
  }
  static const char *const nm[] = {"ln", "log10"};
  out = nvec(2, nm, o);
  return true;
}


bool skp(Ctx &c, Val &out) {
  Dv v;
  CK(num(c, c["x"], "x", v));
  const int n = static_cast<int>(v.n());
  double s = 0.0, a = 0.0;
  CK(one(c, c["s2"], "s2", 0.0, Inf, s));
  CK(oned(c, "alpha", 0.05, 0.0, 1.0, a));
  int ns = 0;
  CK(wholed(c, "nsim", 999, 1, ns));
  unsigned int sv;
  CK(seed(c, sv));
  std::vector<int> grp(v.n()), iw(3 * v.n());
  std::vector<double> wk(2 * v.n());
  int ng = 0;
  if (clab_skp(v.p(), n, s, a, ns, sv, grp.data(), &ng, iw.data(), wk.data()) != CLAB_OK) {
    return c.fail("the kernel refused the input.");
  }
  out = unbox(vints(grp.data(), grp.size()));
  return true;
}


bool lwss(Ctx &c, Val &out) {
  Args a(c, {"x", "y"});
  CK(a.byname(c));
  Dv xv, yv;
  CK(num(c, a[0], "x", xv));
  const std::size_t n = xv.n();
  CK(num(c, a[1], "y", yv, static_cast<long>(n)));
  double sp = 0.0;
  CK(oned(c, "f", 2.0 / 3.0, 0.0, Inf, sp));
  int it = 0;
  CK(wholed(c, "iter", 3, 0, it));
  std::vector<int> ord(n);
  for (std::size_t i = 0; i < n; ++i) ord[i] = static_cast<int>(i);
  const double *px = xv.p();
  std::stable_sort(ord.begin(), ord.end(), [px](int p, int q) { return px[p] < px[q]; });
  std::vector<double> sx(n), sy(n);
  for (std::size_t i = 0; i < n; ++i) {
    sx[i] = px[ord[i]];
    sy[i] = yv[static_cast<std::size_t>(ord[i])];
  }
  const double dl = c.has("delta") ? val1(c["delta"]) : 0.01 * (sx[n - 1] - sx[0]);
  if (!(dl >= 0) || !std::isfinite(dl)) return c.fail("`delta` must be one number from 0.");
  std::vector<double> fit(n), wk(3 * n), o(n);
  if (clab_lwss(sx.data(), sy.data(), static_cast<int>(n), sp, it, dl, fit.data(), wk.data()) != CLAB_OK) {
    return c.fail("the kernel refused the input.");
  }
  for (std::size_t i = 0; i < n; ++i) o[static_cast<std::size_t>(ord[i])] = fit[i];
  out = vauto(o);
  return true;
}

bool chmb(Ctx &c, Val &out) {
  Dv av, bv;
  CK(num(c, c["a"], "a", av));
  CK(num(c, c["b"], "b", bv));
  const int na = static_cast<int>(av.n()), nb = static_cast<int>(bv.n());
  int wo = 0, wv = 0;
  CK(whole(c, c["wopt"], "wopt", 2, wo));
  CK(whole(c, c["wval"], "wval", 2, wv));
  if (wv > wo) return c.fail("`wval` must be at most `wopt`.");
  double co = 0.0;
  CK(oned(c, "coarse", 0.25, 0.0, Inf, co));
  std::vector<double> wk(3 * av.n() + 3 * bv.n() + 7 * std::max(av.n(), bv.n()));
  int ns = 0, nd = 0;
  double u = 0, p = 0;
  const int rc = clab_chmb(av.p(), na, bv.p(), nb, wo, wv, co, &ns, &nd, &u, &p, wk.data());
  if (rc != CLAB_OK) {
    return c.fail("%s", rc == CLAB_EARG ? "a trimmed profile is shorter than `wopt`."
                                        : "every window is flat: the profiles carry no marks.");
  }
  const double o[4] = {static_cast<double>(ns), static_cast<double>(nd), u, p};
  static const char *const nm[] = {"nsame", "ndiff", "u", "p"};
  out = nvec(4, nm, o);
  return true;
}


bool tfrp(Ctx &c, Val &out) {
  Dv pv;
  CK(num(c, c["par"], "par", pv, 17));
  bool de = false;
  CK(flagd(c, "deffect", true, de));
  static const char *const to[] = {"nb", "nbtrunc", "unif"};
  int td = 0, ym = 0, ns = 0, q = 0;
  CK(pickd(c, "tdist", to, 3, 0, td));
  CK(wholed(c, "ymax", 20, 0, ym));
  CK(wholed(c, "nsim", 1000, 1, ns));
  unsigned int s;
  CK(seed(c, s));
  const double h = val1d(c, "hl", 0.02);
  CK(wholed(c, "nq", 256, 1, q));
  int wl = 0;
  if (clab_tfrw(pv.p(), de ? 1 : 0, td, ym, h, q, &wl) != CLAB_OK) {
    return c.fail("`par` is outside the model's domain, or `hl` is outside 1e-4 to 1, or `nq` above 100000.");
  }
  std::vector<double> wk(static_cast<std::size_t>(wl)), p(static_cast<std::size_t>(ym) + 1),
      se(static_cast<std::size_t>(ym) + 1), y(static_cast<std::size_t>(ym) + 1);
  double cut[2];
  if (clab_tfrp(pv.p(), de ? 1 : 0, td, ym, ns, s, h, q, p.data(), se.data(), cut, wk.data()) != CLAB_OK) {
    return c.fail("the kernel refused the input.");
  }
  for (int j = 0; j <= ym; ++j) y[static_cast<std::size_t>(j)] = j;
  Frame f;
  f.num("y", y);
  f.num("p", p);
  f.num("se", se);
  static const char *const cn[] = {"time", "grid"};
  out = vobj();
  put(out, "prob", f.o);
  put(out, "cut", nvec(2, cn, cut));
  return true;
}
