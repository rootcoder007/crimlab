// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Internal benchmarking of officers (src/bibn.cpp).

// xcols(): the circumstances, null, an object of columns, a matrix of n
// rows or one vector; column-major into v; the number of columns.
bool ibx(Ctx &c, const Val &x, std::size_t n, std::vector<double> &v, int &p) {
  v.clear();
  p = 0;
  if (x.t == Val::NUL) return true;
  int dim[2] = {-1, -1};
  if (x.t == Val::ARR && !x.a.empty() && x.a[0].t == Val::ARR) {
    if (!shape(x, 2, dim, 0) || !isnum(x.a[0]) || static_cast<std::size_t>(dim[0]) != n) {
      return c.fail("`x` must be a numeric matrix with one row per stop, %d.", static_cast<int>(n));
    }
    v.assign(n * zz(dim[1]), 0.0);
    bool bad = false;
    fill(x, 2, dim, 0, 0, 1, v, bad, true);
    for (double e : v) {
      if (!std::isfinite(e)) return c.fail("`x` must hold finite numbers with no missing values.");
    }
    p = dim[1];
    return true;
  }
  std::vector<const Val *> cs;
  if (x.t == Val::OBJ) {
    for (const Val &e : x.a) cs.push_back(&e);
  } else {
    cs.push_back(&x);
  }
  for (const Val *e : cs) {
    Dv d;
    CK(num(c, *e, "x", d, static_cast<std::int64_t>(n)));
    v.insert(v.end(), d.v.begin(), d.v.end());
  }
  p = static_cast<int>(cs.size());
  return true;
}

bool ibench(Ctx &c, Val &out) {
  const std::size_t n = vlen(c["y"]);
  std::vector<int> yb;
  CK(bin01(c, c["y"], "y", n, yb));
  if (n < 2) return c.fail("`y` must hold two stops at least.");
  Lab o;
  CK(codes(c, c["off"], "off", n, o));
  if (o.g < 2) return c.fail("`off` must hold at least two officers.");
  std::vector<double> xv;
  int p = 0;
  CK(ibx(c, c["x"], n, xv, p));
  if (p > 8190) return c.fail("`x` must have at most 8190 columns.");
  int mn = 0;
  CK(wholed(c, "minn", 20, 1, mn));
  double fd = 0.0;
  CK(oned(c, "fdr", 0.05, 0.0, 1.0, fd));
  bool udr = false, ulf = false;
  CK(flagd(c, "dr", false, udr));
  CK(flagd(c, "local", false, ulf));
  std::vector<double> yd(yb.begin(), yb.end());
  const int no = o.g;
  std::vector<double> res(8 * zz(no)), nl3(3);
  std::vector<int> ok(zz(no));
  const int rc = clab_ibench(yd.data(), o.code.data(), p ? xv.data() : nullptr, static_cast<int>(n), p, no, mn, udr,
                             ulf, res.data(), ok.data(), nl3.data());
  if (rc == CLAB_EDOM && ulf) {
    return c.fail("the local false discovery rate could not be estimated: the null could not be fitted to the "
                  "middle of the officers' z statistics.");
  }
  CK(kchk(c, rc));
  Frame t;
  t.str("off", o.lev);
  static const char *const cn[] = {"n", "rate", "bench", "ess", "z", "p", "q", "lfdr"};
  for (int k = 0; k < (ulf ? 8 : 7); ++k) t.num(cn[k], nav(&res[zz(k) * zz(no)], zz(no)));
  std::vector<int> fl(zz(no));
  for (std::size_t j = 0; j < zz(no); ++j) {
    const double q = res[(ulf ? 7 : 6) * zz(no) + j];
    fl[j] = std::isnan(q) ? -1 : q < fd;
  }
  t.lgl("flag", fl.data(), fl.size());
  out = vobj();
  put(out, "tab", t.o);
  put(out, "fdr", vnum(fd));
  put(out, "minn", vnum(mn));
  put(out, "nx", vnum(p));
  put(out, "dr", vbool(udr));
  static const char *const nn[] = {"delta", "sigma", "p0"};
  put(out, "null", ulf ? nvec(3, nn, nl3.data()) : Val());
  return true;
}
