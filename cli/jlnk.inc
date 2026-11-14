// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Crime linkage (src/blnk.cpp).

// mocols(): 0/1 behaviours, an array of rows or an object of columns; n
// rows, column-major into v; the number of columns.
bool lkmo(Ctx &c, const Val &mo, std::size_t n, std::vector<int> &v, int &q) {
  v.clear();
  q = 0;
  if (mo.t == Val::NUL) return true;
  auto bit = [](const Val &e, int &b) {
    if (e.t == Val::BOOL) b = e.b ? 1 : 0;
    else if (e.t == Val::NUM && (e.n == 0 || e.n == 1)) b = static_cast<int>(e.n);
    else return false;
    return true;
  };
  if (mo.t == Val::ARR && !mo.a.empty() && mo.a[0].t == Val::ARR) {
    const std::size_t k = mo.a[0].a.size();
    if (mo.a.size() != n) return c.fail("`mo` must be a 0/1 or logical matrix with one row per crime, %d.", static_cast<int>(n));
    v.assign(n * k, 0);
    for (std::size_t i = 0; i < n; ++i) {
      if (mo.a[i].t != Val::ARR || mo.a[i].a.size() != k) {
        return c.fail("`mo` must be a 0/1 or logical matrix with one row per crime, %d.", static_cast<int>(n));
      }
      for (std::size_t j = 0; j < k; ++j) {
        int b = 0;
        if (!bit(mo.a[i].a[j], b)) return c.fail("`mo` must hold 0 or 1 (or FALSE or TRUE), with no missing values.");
        v[j * n + i] = b;
      }
    }
    q = static_cast<int>(k);
    return true;
  }
  if (mo.t != Val::OBJ) return c.fail("`mo` must be a matrix, a data frame or a list of 0/1 columns.");
  for (const Val &e : mo.a) {
    std::vector<int> b;
    CK(bin01(c, e, "mo", n, b));
    v.insert(v.end(), b.begin(), b.end());
  }
  q = static_cast<int>(mo.a.size());
  return true;
}

bool lnkpr(Ctx &c, Val &out) {
  Pts p;
  CK(xy(c, PLANAR, p));
  const int n = p.n;
  if (n < 2) return c.fail("`x` must hold two crimes at least.");
  if (n > 3000) return c.fail("at most 3000 crimes: their pairs are kept in memory.");
  Dv tv;
  CK(num(c, c["t"], "t", tv, n));
  const Val &sv = c["series"];
  if (vlen(sv) != zz(n)) return c.fail("`series` must have one label per crime, %d.", n);
  std::vector<int> ser(zz(n), -1);
  {
    std::vector<std::string> s(zz(n));
    std::vector<char> has(zz(n), 0);
    std::vector<std::string> at;
    for (std::size_t i = 0; i < zz(n); ++i) {
      const Val &e = sv.t == Val::ARR || sv.t == Val::OBJ ? sv.a[i] : sv;
      if (e.t == Val::NUL) continue;
      if (!lab(e, s[i])) return c.fail("`series` must be labels (numbers, text or a factor), NA for an unsolved crime.");
      has[i] = 1;
      at.push_back(s[i]);
    }
    std::sort(at.begin(), at.end());
    at.erase(std::unique(at.begin(), at.end()), at.end());
    for (std::size_t i = 0; i < zz(n); ++i) {
      if (has[i]) ser[i] = static_cast<int>(std::lower_bound(at.begin(), at.end(), s[i]) - at.begin());
    }
  }
  std::vector<int> mv;
  int q = 0;
  CK(lkmo(c, c["mo"], zz(n), mv, q));
  const std::size_t m = zz(n) * zz(n - 1) / 2;
  std::vector<int> pi(m), pj(m), lk(m);
  std::vector<double> di(m), dtv(m), ja(m), pr(m);
  const int p3 = q > 0 ? 3 : 2;
  double cf[4], se[4], au[4], dev = 0.0;
  const int rc = clab_lnkpr(p.x.data(), p.y.data(), tv.p(), q ? mv.data() : nullptr, q, ser.data(), n, pi.data(),
                            pj.data(), di.data(), dtv.data(), ja.data(), lk.data(), pr.data(), cf, se, au, &dev);
  if (rc == CLAB_EINTR) return kchk(c, rc);
  if (rc != CLAB_OK) return c.fail("the solved crimes need linked and unlinked pairs that do not separate perfectly.");
  for (std::size_t k = 0; k < m; ++k) {
    ++pi[k];
    ++pj[k];
    ja[k] = na(ja[k]);
  }
  static const char *const tn[] = {"(Intercept)", "dist", "dt", "jac"};
  Frame cf2;
  std::vector<std::string> term;
  std::vector<double> es, ss, zv, pv;
  for (int j = 0; j <= p3; ++j) {
    term.push_back(tn[j]);
    es.push_back(cf[j]);
    ss.push_back(se[j]);
    zv.push_back(cf[j] / se[j]);
    pv.push_back(2 * clab_pnorm(-std::fabs(cf[j] / se[j])));
  }
  cf2.str("term", term);
  cf2.num("est", es);
  cf2.num("se", ss);
  cf2.num("z", zv);
  cf2.num("p", pv);
  for (double &e : au) e = na(e);
  static const char *const an[] = {"dist", "dt", "jac", "model"};
  Frame f;
  f.ints("i", pi.data(), m);
  f.ints("j", pj.data(), m);
  f.num("dist", di);
  f.num("dt", dtv);
  f.num("jac", ja);
  f.lgl("linked", lk.data(), m);
  f.num("prob", pr);
  out = vobj();
  put(out, "coef", cf2.o);
  put(out, "auc", nvec(4, an, au));
  put(out, "pairs", f.o);
  put(out, "dev", vnum(dev));
  put(out, "nmo", vnum(q));
  return true;
}
