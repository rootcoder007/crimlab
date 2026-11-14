// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Hot spots on grids and in counts: local Gi* on a grid and the spatial
// point pattern test.

bool gstr(Ctx &c, Val &out) {
  int nr, nc;
  std::vector<double> cv, ev;
  CK(nmat(c, c["cnt"], "cnt", nr, nc, cv));
  if (!c.has("ex")) {
    ev.assign(cv.size(), 1.0);
  } else {
    int er, ec;
    CK(nmat(c, c["ex"], "ex", er, ec, ev));
    if (er != nr || ec != nc) return c.fail("`ex` must have the dimensions of `cnt`.");
  }
  if (nr > 46340 || nc > 46340) return c.fail("`cnt` must have at most 46340 rows and columns.");
  const double rd = val1d(c, "rad", 1.5);
  if (!(rd >= 0 && rd <= 1000)) return c.fail("`rad` must be one number from 0 to 1000.");
  int ns;
  CK(wholed(c, "nsim", 0, 0, ns));
  unsigned int s = 0;
  if (ns > 0) CK(seed(c, s));
  const std::size_t m = static_cast<std::size_t>(nr) * static_cast<std::size_t>(nc);
  const std::size_t q = 2 * static_cast<std::size_t>(std::floor(rd)) + 1;
  std::vector<double> o(7 * m), wk(m);
  std::vector<int> iw(2 * q * q + 3 * m);
  if (clab_gstr(cv.data(), ev.data(), nr, nc, rd, ns, s, o.data(), wk.data(), iw.data()) != CLAB_OK) {
    return c.fail("`cnt` must be whole numbers from 0, `ex` numbers from 0, with at least 2 cells of positive exposure.");
  }
  static const char *const nm[] = {"gz", "pz", "qz", "pp", "ze", "pe", "qe"};
  std::vector<double> rw(m), cl(m);
  for (std::size_t i = 0; i < m; ++i) {
    rw[i] = static_cast<double>(i % static_cast<std::size_t>(nr) + 1);
    cl[i] = static_cast<double>(i / static_cast<std::size_t>(nr) + 1);
  }
  Frame f;
  f.num("row", rw);
  f.num("col", cl);
  for (int j = 0; j < 7; ++j) f.num(nm[j], o.data() + static_cast<std::size_t>(j) * m, m);
  out = f.o;
  return true;
}


bool sptst(Ctx &c, Val &out) {
  Args a(c, {"a", "b"});
  a.v[0] = gvec(a.v[0]);
  a.v[1] = gvec(a.v[1]);
  CK(a.byname(c));
  if (!isnum(a[0]) || !isnum(a[1]) || a[0].t == Val::NUL || a[1].t == Val::NUL) {
    return c.fail("`a` and `b` must be counts.");
  }
  const std::vector<double> av = dbl(a[0]), bv = dbl(a[1]);
  const int n = static_cast<int>(av.size());
  if (n < 1 || bv.size() != av.size()) return c.fail("`a` and `b` must have the same length, from 1.");
  int nbt;
  CK(wholed(c, "nboot", 200, 1, nbt));
  unsigned int s;
  CK(seed(c, s));
  double al, sc;
  CK(oned(c, "alpha", 0.05, 0.0, 1.0, al));
  CK(oned(c, "scale", 100, 0.0, Inf, sc));
  const std::size_t N = av.size();
  std::vector<double> o(3 * N), wk(6 * N);
  std::vector<int> iw(4 * N);
  double g[8];
  if (clab_sptst(av.data(), bv.data(), n, nbt, s, al, sc, o.data(), g, wk.data(), iw.data()) != CLAB_OK) {
    return c.fail("`a` and `b` must be whole numbers from 0 (NA leaves a unit out).");
  }
  Frame f;
  f.num("diff", o.data(), N);
  f.num("p", o.data() + N, N);
  f.num("flag", o.data() + 2 * N, N);
  static const char *const gn[] = {"left", "empty", "elig", "flagged", "s", "z", "p", "nboot"};
  out = vobj();
  put(out, "units", f.o);
  put(out, "glob", nvec(8, gn, g));
  return true;
}
