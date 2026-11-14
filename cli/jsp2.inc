// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Spatial weights and the statistics built on them (src/bsp2.cpp), the
// directional mean, head-bang smoothing, the correlated random walk and
// risk-adjusted nearest-neighbour clustering (src/bsp3.cpp).

// Weights as compressed sparse rows, rows in increasing order of neighbour.
struct Wcsr {
  std::vector<int> rp, ci;
  std::vector<double> w;
};
bool wbuild(Ctx &c, int n, const std::vector<int> &ei, const std::vector<int> &ej, const std::vector<double> &ew,
            Wcsr &o) {
  const std::size_t m = ei.size();
  std::vector<int> ord(m);
  for (std::size_t q = 0; q < m; ++q) ord[q] = static_cast<int>(q);
  std::sort(ord.begin(), ord.end(), [&](int a, int b) {
    const std::size_t i = static_cast<std::size_t>(a), j = static_cast<std::size_t>(b);
    return ei[i] < ei[j] || (ei[i] == ei[j] && ej[i] < ej[j]);
  });
  o.rp.assign(static_cast<std::size_t>(n) + 1, 0);
  o.ci.assign(m, 0);
  o.w.assign(m, 0.0);
  for (std::size_t q = 0; q < m; ++q) {
    const std::size_t a = static_cast<std::size_t>(ord[q]);
    if (ei[a] == ej[a]) return c.fail("a unit cannot be its own neighbour in `w`.");
    if (q > 0) {
      const std::size_t b = static_cast<std::size_t>(ord[q - 1]);
      if (ei[a] == ei[b] && ej[a] == ej[b]) return c.fail("`w` lists a pair of units twice.");
    }
    o.ci[q] = ej[a];
    o.w[q] = ew[a];
    ++o.rp[static_cast<std::size_t>(ei[a]) + 1];
  }
  for (std::size_t i = 0; i < static_cast<std::size_t>(n); ++i) o.rp[i + 1] += o.rp[i];
  return true;
}

// wread(): the weights spwt() returns, for n units.
bool wread(Ctx &c, int n, Wcsr &o) {
  const Val &w = c["w"];
  const Val &nn = col(w, "n");
  if (nn.t == Val::NUL) return c.fail("`w` must be the list spwt() returns or an spdep listw.");
  if (val1(nn) != n) return c.fail("`w` is for %g units, not %d.", val1(nn), n);
  Dv i, j, v;
  CK(num(c, col(w, "i"), "w$i", i, -1, 1.0));
  CK(num(c, col(w, "j"), "w$j", j, static_cast<long>(i.n()), 1.0));
  CK(num(c, col(w, "w"), "w$w", v, static_cast<long>(i.n()), 0.0));
  std::vector<int> ei, ej;
  for (std::size_t q = 0; q < i.n(); ++q) {
    const double a = i[q], b = j[q];
    if (a > n || b > n || a != std::floor(a) || b != std::floor(b)) {
      return c.fail("`w$i` and `w$j` must be unit numbers from 1 to %d.", n);
    }
    ei.push_back(static_cast<int>(a) - 1);
    ej.push_back(static_cast<int>(b) - 1);
  }
  return wbuild(c, n, ei, ej, v.v, o);
}

bool spwt(Ctx &c, Val &out) {
  static const char *const to[] = {"knn", "band", "contig"};
  static const char *const so[] = {"row", "raw"};
  int tp = 0, rs = 0, me = 0;
  CK(pickd(c, "type", to, 3, 0, tp));
  CK(pickd(c, "style", so, 2, 0, rs));
  rs = rs == 0 ? 1 : 0;
  CK(metric(c, me));
  const double pw = val1d(c, "idw", 0);
  if (!(pw >= 0) || !std::isfinite(pw)) return c.fail("`idw` must be one number from 0.");
  const bool pts = c.has("x");
  Pts p;
  if (pts) CK(xy(c, PROJ, p));
  int n = 0;
  std::vector<int> rp, ci;
  std::vector<double> dv;
  if (tp == 2) {
    const Val &nb = c["nb"];
    if (nb.t != Val::ARR) return c.fail("`nb` must be a list of neighbour numbers, one element per unit.");
    if (nb.a.size() > static_cast<std::size_t>(INT_MAX)) return c.fail("`nb` is too long for crimlab.");
    n = static_cast<int>(nb.a.size());
    if (pts && p.n != n) return c.fail("`nb` must have one element per point.");
    if (!pts && pw > 0) return c.fail("inverse-distance weights need the points `x`.");
    std::vector<int> ei, ej;
    for (int i = 0; i < n; ++i) {
      const Val &a = nb.a[static_cast<std::size_t>(i)];
      if (!isnum(a)) return c.fail("`nb` must be a list of neighbour numbers, one element per unit.");
      const std::vector<double> av = dbl(a);
      for (double j : av) {
        if (j == 0 && av.size() == 1) break;
        if (!(j >= 1 && j <= n && j == std::floor(j))) return c.fail("`nb` must hold unit numbers from 1 to %d.", n);
        ei.push_back(i);
        ej.push_back(static_cast<int>(j) - 1);
      }
    }
    std::vector<double> ed(ei.size());
    for (std::size_t q = 0; q < ei.size(); ++q) {
      const std::size_t a = static_cast<std::size_t>(ei[q]), b = static_cast<std::size_t>(ej[q]);
      ed[q] = pts ? (me ? std::fabs(p.x[a] - p.x[b]) + std::fabs(p.y[a] - p.y[b]) : std::hypot(p.x[a] - p.x[b], p.y[a] - p.y[b]))
                  : -HUGE_VAL;
    }
    Wcsr o;
    CK(wbuild(c, n, ei, ej, ed, o));
    rp = o.rp;
    ci = o.ci;
    dv = o.w;
  } else {
    if (!pts) return c.fail("`x` must give the points for `type = \"%s\"`.", to[tp]);
    n = p.n;
    if (n < 2) return c.fail("at least 2 points are needed.");
    std::vector<int> iw(2 * static_cast<std::size_t>(n) + 2);
    if (tp == 0) {
      int kk = 0;
      CK(wholed(c, "k", 4, 1, kk));
      if (kk > n - 1) return c.fail("`k` must be below the number of points.");
      const std::size_t m = static_cast<std::size_t>(n) * static_cast<std::size_t>(kk);
      rp.assign(static_cast<std::size_t>(n) + 1, 0);
      ci.assign(m, 0);
      dv.assign(m, 0.0);
      CK(kc(c, clab_swknb(p.x.data(), p.y.data(), n, kk, me, ci.data(), dv.data(), iw.data()),
            "the kernel refused the input."));
      for (int i = 0; i <= n; ++i) rp[static_cast<std::size_t>(i)] = i * kk;
      for (int i = 0; i < n; ++i) {
        int *cc = ci.data() + static_cast<std::size_t>(i) * kk;
        double *e = dv.data() + static_cast<std::size_t>(i) * kk;
        for (int u = 1; u < kk; ++u) {
          const int j = cc[u];
          const double f = e[u];
          int v = u - 1;
          while (v >= 0 && cc[v] > j) {
            cc[v + 1] = cc[v];
            e[v + 1] = e[v];
            --v;
          }
          cc[v + 1] = j;
          e[v + 1] = f;
        }
      }
    } else {
      double dm = 0.0;
      CK(one(c, c["d"], "d", 0.0, Inf, dm));
      int nnz = 0;
      rp.assign(static_cast<std::size_t>(n) + 1, 0);
      CK(kc(c, clab_swband(p.x.data(), p.y.data(), n, dm, me, 0, rp.data(), nullptr, nullptr, &nnz, iw.data()),
            "the kernel refused the input."));
      ci.assign(static_cast<std::size_t>(nnz) + 1, 0);
      dv.assign(static_cast<std::size_t>(nnz) + 1, 0.0);
      CK(kc(c, clab_swband(p.x.data(), p.y.data(), n, dm, me, nnz, rp.data(), ci.data(), dv.data(), &nnz, iw.data()),
            "the kernel refused the input."));
      ci.resize(static_cast<std::size_t>(nnz));
      dv.resize(static_cast<std::size_t>(nnz));
    }
  }
  const std::size_t m = ci.size();
  std::vector<double> wv(m + 1);
  CK(kc(c, clab_swwt(rp.data(), dv.data(), n, pw, rs, wv.data()),
        "inverse-distance weights need distinct points: two neighbours coincide."));
  wv.resize(m);
  std::vector<int> iv(m), jv(m);
  for (std::size_t i = 0; i < static_cast<std::size_t>(n); ++i) {
    for (int q = rp[i]; q < rp[i + 1]; ++q) {
      iv[static_cast<std::size_t>(q)] = static_cast<int>(i) + 1;
      jv[static_cast<std::size_t>(q)] = ci[static_cast<std::size_t>(q)] + 1;
    }
  }
  if (tp == 2 && !pts) {
    for (double &e : dv) e = NA;
  }
  out = vobj();
  put(out, "n", vnum(n));
  put(out, "i", unbox(vints(iv.data(), m)));
  put(out, "j", unbox(vints(jv.data(), m)));
  put(out, "d", vauto(dv));
  put(out, "w", vauto(wv));
  return true;
}

// x, and the weights for its values: gearc(), gorg(), lmorn(), locgi().
bool xw(Ctx &c, double lo, int nmin, Dv &xv, Wcsr &w) {
  CK(num(c, c["x"], "x", xv, -1, lo));
  const int n = static_cast<int>(xv.n());
  if (n < nmin) return c.fail("at least %d values are needed.", nmin);
  return wread(c, n, w);
}

bool gearc(Ctx &c, Val &out) {
  Dv xv;
  Wcsr w;
  CK(xw(c, -HUGE_VAL, 4, xv, w));
  int ns = 0;
  CK(wholed(c, "nsim", 0, 0, ns));
  unsigned int s = 0;
  if (ns > 0) CK(seed(c, s));
  std::vector<double> wk(xv.n());
  double o[9];
  CK(kc(c, clab_gearc(xv.p(), static_cast<int>(xv.n()), w.rp.data(), w.ci.data(), w.w.data(), ns, s, wk.data(), o),
        "`x` is constant, or `w` has no positive weight: C is undefined."));
  static const char *const nm[] = {"c", "vn", "vr", "zn", "zr", "pn", "pr", "pperm", "s0"};
  out = nvec(9, nm, o);
  return true;
}

bool gorg(Ctx &c, Val &out) {
  Dv xv;
  Wcsr w;
  CK(xw(c, 0.0, 4, xv, w));
  int ns = 0;
  CK(wholed(c, "nsim", 0, 0, ns));
  unsigned int s = 0;
  if (ns > 0) CK(seed(c, s));
  std::vector<double> wk(xv.n());
  double o[6];
  CK(kc(c, clab_gorg(xv.p(), static_cast<int>(xv.n()), w.rp.data(), w.ci.data(), w.w.data(), ns, s, wk.data(), o),
        "G is undefined: fewer than 2 positive values of `x`, or no positive weight."));
  static const char *const nm[] = {"g", "e", "v", "z", "p", "pperm"};
  out = nvec(6, nm, o);
  return true;
}

bool lmorn(Ctx &c, Val &out) {
  Dv xv;
  Wcsr w;
  CK(xw(c, -HUGE_VAL, 3, xv, w));
  int ns = 0;
  CK(wholed(c, "nsim", 0, 0, ns));
  unsigned int s = 0;
  if (ns > 0) CK(seed(c, s));
  const std::size_t n = xv.n();
  std::vector<int> iw(n);
  std::vector<double> wk(n), o(7 * n);
  CK(kc(c, clab_lmorn(xv.p(), static_cast<int>(n), w.rp.data(), w.ci.data(), w.w.data(), ns, s, iw.data(), wk.data(),
                      o.data()),
        "`x` is constant: local Moran is undefined."));
  static const char *const nm[] = {"ii", "e", "v", "z", "p", "pperm", "quad"};
  Frame f;
  for (int j = 0; j < 7; ++j) f.num(nm[j], o.data() + static_cast<std::size_t>(j) * n, n);
  out = f.o;
  return true;
}

bool locgi(Ctx &c, Val &out) {
  Dv xv;
  Wcsr w;
  CK(xw(c, -HUGE_VAL, 3, xv, w));
  bool st = false;
  CK(flagd(c, "star", true, st));
  std::vector<double> o(xv.n());
  CK(kc(c, clab_locgi(xv.p(), static_cast<int>(xv.n()), w.rp.data(), w.ci.data(), w.w.data(), st ? 1 : 0, o.data()),
        "the kernel refused the input."));
  out = vauto(o);
  return true;
}

bool dirmn(Ctx &c, Val &out) {
  Args a(c, {"x", "y", "angle", "w"});
  CK(a.byname(c));
  if ((a[2].t == Val::NUL) == (a[0].t == Val::NUL)) {
    return c.fail("give the points `x` (with `origin`) or the angles `angle`, not both.");
  }
  std::vector<double> an;
  if (a[2].t != Val::NUL) {
    Dv d;
    CK(num(c, a[2], "angle", d));
    an = d.v;
  } else {
    Pts p;
    CK(xy(c, a[0], a[1], PROJ, p));
    Val d1 = varr();
    d1.a = {vnum(0), vnum(0)};
    Dv o;
    CK(num(c, c.has("origin") ? c["origin"] : d1, "origin", o, 2));
    double ox = o[0], oy = o[1];
    CK(fwd(c, p, &ox, &oy, 1));
    for (int i = 0; i < p.n; ++i) {
      const double dx = p.x[static_cast<std::size_t>(i)] - ox, dy = p.y[static_cast<std::size_t>(i)] - oy;
      if (dx == 0 && dy == 0) return c.fail("a point lies on `origin`: it has no direction.");
      an.push_back(std::atan2(dy, dx));
    }
  }
  const int n = static_cast<int>(an.size());
  Dv wv;
  if (a[3].t != Val::NUL) CK(num(c, a[3], "w", wv, n, 0.0));
  double o[5];
  CK(kc(c, clab_dirmn(an.data(), wv.n() ? wv.p() : nullptr, n, o), "`w` must have a positive sum."));
  double b = 90.0 - o[0] * 180.0 / 3.14159265358979323846;
  if (b < 0) b += 360.0;
  const double r[7] = {o[0], b, o[1], o[2], o[3], o[4], static_cast<double>(n)};
  static const char *const nm[] = {"theta", "bearing", "r", "var", "sd", "w", "n"};
  out = nvec(7, nm, r);
  return true;
}

bool hdbang(Ctx &c, Val &out) {
  Args a(c, {"x", "y", "v", "w"});
  CK(a.byname(c));
  Pts p;
  CK(xy(c, a[0], a[1], PROJ, p));
  const int n = p.n;
  if (n < 3) return c.fail("at least 3 points are needed.");
  Dv vv, wv;
  CK(num(c, a[2], "v", vv, n));
  if (a[3].t != Val::NUL) CK(num(c, a[3], "w", wv, n, 0.0, true));
  int k = 0, nt = 0, it = 0;
  CK(wholed(c, "nn", 8, 2, k));
  if (k > n - 1) return c.fail("`nn` must be below the number of points.");
  const double an = val1d(c, "ang", 135);
  if (!(an >= 90 && an <= 180)) return c.fail("`ang` must be one number from 90 to 180.");
  CK(wholed(c, "ntrip", 0, 0, nt));
  CK(wholed(c, "iter", 1, 0, it));
  const std::size_t T = static_cast<std::size_t>(k) * static_cast<std::size_t>(k - 1) / 2, N = static_cast<std::size_t>(n);
  std::vector<int> iw(2 * N + 2 + N * static_cast<std::size_t>(k) + 3 * T);
  std::vector<double> wk(N * static_cast<std::size_t>(k) + N + 5 * T), o(N);
  CK(kc(c, clab_hdbang(p.x.data(), p.y.data(), vv.p(), wv.n() ? wv.p() : nullptr, n, k, an, nt, it, iw.data(), wk.data(),
                       o.data()),
        "the kernel refused the input."));
  out = vauto(o);
  return true;
}

bool crwalk(Ctx &c, Val &out) {
  Pts p;
  CK(xy(c, PROJ, p));
  const int n = p.n;
  if (n < 3) return c.fail("at least 3 locations are needed.");
  int hh = 0, ns = 0;
  CK(wholed(c, "h", 1, 1, hh));
  CK(wholed(c, "nsim", 999, 1, ns));
  unsigned int sd = 0;
  CK(seed(c, sd));
  std::vector<double> st(static_cast<std::size_t>(n) - 1), tu(static_cast<std::size_t>(n) - 2);
  double f[5];
  CK(kc(c, clab_crwfit(p.x.data(), p.y.data(), n, st.data(), tu.data(), f),
        "two successive locations coincide: a repeated location has no heading; drop repeats first."));
  const std::size_t H = static_cast<std::size_t>(hh), S = static_cast<std::size_t>(ns);
  std::vector<double> obs(H), ex(H), kk(H), sx(H), sy(H), qx(H), qy(H);
  CK(kc(c, clab_crwmsd(p.x.data(), p.y.data(), n, f[0], f[1], f[2], hh, obs.data(), ex.data()),
        "the kernel refused the input."));
  std::vector<double> px(H * S), py(H * S);
  CK(kc(c, clab_crwsim(st.data(), n - 1, tu.data(), n - 2, p.x[static_cast<std::size_t>(n) - 1],
                       p.y[static_cast<std::size_t>(n) - 1], f[4], hh, ns, sd, px.data(), py.data()),
        "the kernel refused the input."));
  for (std::size_t k = 0; k < H; ++k) {
    double mx = 0, my = 0, vx = 0, vy = 0;
    const double *ax = px.data() + k * S, *ay = py.data() + k * S;
    for (std::size_t s = 0; s < S; ++s) {
      mx += ax[s];
      my += ay[s];
    }
    mx /= ns;
    my /= ns;
    for (std::size_t s = 0; s < S; ++s) {
      vx += (ax[s] - mx) * (ax[s] - mx);
      vy += (ay[s] - my) * (ay[s] - my);
    }
    kk[k] = static_cast<double>(k + 1);
    sx[k] = mx;
    sy[k] = my;
    qx[k] = ns > 1 ? std::sqrt(vx / (ns - 1)) : NaN;
    qy[k] = ns > 1 ? std::sqrt(vy / (ns - 1)) : NaN;
  }
  back(p, sx.data(), sy.data(), hh);
  const double fv[6] = {static_cast<double>(n - 1), f[0], f[1], f[2], f[3], f[4]};
  static const char *const fn[] = {"steps", "ml", "ml2", "c", "s", "head"};
  Frame m;
  m.num("k", kk);
  m.num("obs", obs);
  m.num("exp", ex);
  m.num("x", sx);
  m.num("y", sy);
  m.num("sdx", qx);
  m.num("sdy", qy);
  Val sim = vobj();
  put(sim, "x", vmat(px, ns, hh));
  put(sim, "y", vmat(py, ns, hh));
  out = vobj();
  put(out, "fit", nvec(6, fn, fv));
  put(out, "msd", m.o);
  put(out, "sim", sim);
  return true;
}

bool rnnhcl(Ctx &c, Val &out) {
  Pts p;
  CK(xy(c, PROJ, p));
  const int n = p.n;
  if (n < 2) return c.fail("at least 2 points are needed.");
  Pts b;
  CK(xy(c, c["base"], nul, p.ll ? RAW : PLANAR, b));
  std::vector<double> bx = b.x, by = b.y;
  CK(fwd(c, p, bx.data(), by.data(), b.n));
  Dv bv;
  if (c.has("bw")) CK(num(c, c["bw"], "bw", bv, b.n, 0.0));
  std::vector<double> ones(static_cast<std::size_t>(b.n), 1.0);
  if (bv.n()) ones = bv.v;
  double hb = 0.0, pv = 0.0, cc = 0.0;
  CK(one(c, c["h"], "h", 0.0, Inf, hb));
  CK(oned(c, "p", 0.5, 0.0, 1.0, pv));
  const double z = qnorm(pv);
  int mp = 0, kk = 0;
  CK(wholed(c, "minpts", 10, 2, mp));
  CK(wholed(c, "k", 1, 0, kk));
  if (kk > 2) return c.fail("`k` must be 0, 1 or 2.");
  CK(oned(c, "c", 1, 0.0, Inf, cc));
  const std::size_t N = static_cast<std::size_t>(n);
  std::vector<int> iw(4 * N), cl(N);
  std::vector<double> wk(3 * N), thr(N), cs(8 * N);
  int nc = 0;
  CK(kc(c, clab_rnnhcl(p.x.data(), p.y.data(), n, bx.data(), by.data(), ones.data(), b.n, hb, z, mp, kk, cc, iw.data(),
                       wk.data(), thr.data(), cl.data(), &nc, cs.data()),
        "`bw` must have a positive sum, and `p` must leave a positive threshold."));
  for (int &e : cl) e = e >= 0 ? e + 1 : INT_MIN;
  static const char *const cn[] = {"points", "mx", "my", "wx", "wy", "a", "b", "angle"};
  out = vobj();
  put(out, "cl", unbox(vints(cl.data(), N)));
  put(out, "clusters", ctab(p, cs, nc, 8, cn, 1));
  put(out, "thr", vauto(thr));
  return true;
}
