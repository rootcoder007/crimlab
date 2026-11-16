// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The split-population survival model (src/bspc.cpp).

// xcov(): covariates, null, an object of columns, a matrix of n rows or
// one vector; column-major into v, their names into nm.
bool spx(Ctx &c, const Val &x, std::size_t n, std::vector<double> &v, std::vector<std::string> &nm) {
  v.clear();
  nm.clear();
  if (x.t == Val::NUL) return true;
  int dim[2] = {-1, -1};
  if (x.t == Val::ARR && !x.a.empty() && x.a[0].t == Val::ARR) {
    if (!shape(x, 2, dim, 0) || !isnum(x.a[0]) || static_cast<std::size_t>(dim[0]) != n) {
      return c.fail("`x` must be a numeric matrix with one row per person, %d.", static_cast<int>(n));
    }
    v.assign(n * zz(dim[1]), 0.0);
    bool bad = false;
    fill(x, 2, dim, 0, 0, 1, v, bad, true);
    for (double e : v) {
      if (!std::isfinite(e)) return c.fail("`x` must hold finite numbers with no missing values.");
    }
    for (int j = 0; j < dim[1]; ++j) nm.push_back("x" + std::to_string(j + 1));
    return true;
  }
  std::vector<const Val *> cs;
  std::vector<std::string> names;
  if (x.t == Val::OBJ) {
    for (const Val &e : x.a) cs.push_back(&e);
    names = x.k;
  } else {
    cs.push_back(&x);
  }
  for (std::size_t j = 0; j < cs.size(); ++j) {
    Dv d;
    CK(num(c, *cs[j], "x", d, static_cast<std::int64_t>(n)));
    v.insert(v.end(), d.v.begin(), d.v.end());
    nm.push_back(j < names.size() && !names[j].empty() ? names[j] : "x" + std::to_string(j + 1));
  }
  return true;
}

bool spcure(Ctx &c, Val &out) {
  Dv tv;
  CK(num(c, c["time"], "time", tv, -1, 0.0, true));
  const std::size_t n = tv.n();
  std::vector<int> ev;
  CK(bin01(c, c["event"], "event", n, ev));
  int ne = 0;
  for (int e : ev) ne += e;
  if (ne < 1 || ne >= static_cast<int>(n)) return c.fail("`event` must hold failures and censored times both.");
  std::vector<double> xv;
  std::vector<std::string> xn;
  CK(spx(c, c["x"], n, xv, xn));
  const int p = static_cast<int>(xn.size());
  static const char *const dn[] = {"weibull", "exp", "lognormal"};
  int d = 0;
  CK(pickd(c, "dist", dn, 3, 0, d));
  double lv = 0.0;
  CK(oned(c, "level", 0.95, 0.0, 1.0, lv));
  int mi = 0;
  CK(wholed(c, "maxit", 500, 1, mi));
  std::vector<double> evd(ev.begin(), ev.end());
  const int q = p + 1, np = 2 * q + (d != 1);
  std::vector<double> par(zz(np)), cov(zz(np) * zz(np)), st(5);
  const int rc = clab_spcure(tv.p(), evd.data(), p ? xv.data() : nullptr, static_cast<int>(n), p, d, mi, par.data(),
                             cov.data(), st.data());
  if (rc == CLAB_EINTR) return kchk(c, rc);
  if (rc != CLAB_OK) return c.fail("the likelihood could not be maximised from its start.");
  const double z = qnorm(0.5 + 0.5 * lv);
  std::vector<std::string> part, term;
  std::vector<double> se, lo, hi, zv, pv;
  for (int k = 0; k < np; ++k) {
    const int j = k < 2 * q ? k % q : 0;
    part.push_back(k < q ? "fail" : k < 2 * q ? "time" : "shape");
    term.push_back(k >= 2 * q ? (d == 2 ? "log(sigma)" : "log(shape)") : j == 0 ? "(Intercept)" : xn[zz(j - 1)]);
    const double s = std::sqrt(cov[zz(k) * zz(np) + zz(k)]), e = par[zz(k)];
    se.push_back(na(s));
    lo.push_back(na(e - z * s));
    hi.push_back(na(e + z * s));
    zv.push_back(na(e / s));
    pv.push_back(std::isnan(s) ? NA : 2 * clab_pnorm(-std::fabs(e / s)));
  }
  Frame f;
  f.str("part", part);
  f.str("term", term);
  f.num("est", par);
  f.num("se", se);
  f.num("lo", lo);
  f.num("hi", hi);
  f.num("z", zv);
  f.num("p", pv);
  static const char *const sn[] = {"loglik", "aic", "iter", "conv", "pfail"};
  out = vobj();
  put(out, "coef", f.o);
  put(out, "fit", nvec(5, sn, st.data()));
  put(out, "dist", vstr(dn[d]));
  put(out, "n", vnum(static_cast<double>(n)));
  return true;
}
