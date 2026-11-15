// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Air pollution (src/apol.cpp): dispersion and the health side, with the
// checks of R/apol.R and R/ahlth.R in their words.

// format() of one double, as R prints it with 7 significant digits.
inline std::string rfmt7(double x) {
  if (std::isnan(x)) return "NaN";
  if (std::isinf(x)) return x > 0 ? "Inf" : "-Inf";
  if (x == 0) return "0";
  char b[64];
  std::snprintf(b, sizeof b, "%.6e", x);
  const double v = std::strtod(b, nullptr);
  int d = 1;
  for (; d < 7; ++d) {
    std::snprintf(b, sizeof b, "%.*e", d - 1, x);
    if (std::strtod(b, nullptr) == v) break;
  }
  std::snprintf(b, sizeof b, "%.*e", d - 1, x);
  const int e = std::atoi(std::strchr(b, 'e') + 1);
  const int neg = x < 0;
  const int sci = neg + d + (d > 1) + (std::abs(e) >= 100 ? 5 : 4);
  const int rgt = std::max(0, d - e - 1);
  const int fix = neg + (e >= 0 ? e + 1 : 1) + (rgt > 0 ? rgt + 1 : 0);
  if (fix <= sci) std::snprintf(b, sizeof b, "%.*f", rgt, x);
  return b;
}

// .adnum(): the dispersion inputs.
inline bool adnum(Ctx &c, const Val &o, const char *w, Nv &d, double mn = -HUGE_VAL, bool positive = false,
                  bool scalar = true, bool whole = false) {
  d = Nv();
  bool lg = false;
  const char *cl = o.t == Val::NUL ? nullptr : jcls(o, lg);
  if (lg) {
    bool allna = o.t != Val::BOOL;
    for (const Val &e : o.a) allna = allna && e.t == Val::NUL;
    if (!allna) return ckstop(c, "`%s` must be numeric, not logical", w);
  }
  if (o.t == Val::NUL) return ckstop(c, "`%s` must be numeric, got NULL", w);
  if (cl != nullptr && !lg) return ckstop(c, "`%s` must be numeric, got %s", w, cl);
  Dv t;
  readd(o, t);
  d.v = t.v;
  if (scalar && d.n() != 1) return ckstop(c, "`%s` must be a single number", w);
  bool fin = d.n() > 0;
  for (double x : d.v) fin = fin && std::isfinite(x);
  if (!fin) return ckstop(c, "`%s` must be finite, with no missing values", w);
  for (double x : d.v) {
    if (positive && x <= 0) return ckstop(c, "`%s` must be positive", w);
  }
  for (double x : d.v) {
    if (x < mn) return ckstop(c, "`%s` must be at least %s", w, rfmt7(mn).c_str());
  }
  for (double x : d.v) {
    if (std::fabs(x) > 1e300) return ckstop(c, "`%s` must be below 1e300 in size", w);
  }
  for (double x : d.v) {
    if (whole && x != std::nearbyint(x)) return ckstop(c, "`%s` must be a whole number", w);
  }
  return true;
}

inline bool adnumd(Ctx &c, const char *nm, double def, Nv &d, double mn = -HUGE_VAL, bool positive = false,
                   bool whole = false) {
  if (!c.has(nm)) {
    d = Nv();
    d.v.push_back(def);
    return true;
  }
  return adnum(c, c[nm], nm, d, mn, positive, true, whole);
}

inline bool adcls(Ctx &c, int &k) {
  std::string s = "D";
  const Val &o = c["stab"];
  bool ok = true;
  if (o.t != Val::NUL) {
    const Val *v = o.t == Val::ARR && o.a.size() == 1 ? &o.a[0] : &o;
    ok = lab(*v, s);
  }
  for (char &ch : s) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
  if (!ok || s.size() != 1 || s[0] < 'A' || s[0] > 'F') {
    return ckstop(c, "`stab` must be one of the Pasquill classes \"A\" to \"F\"");
  }
  k = s[0] - 'A';
  return true;
}

inline bool adset(Ctx &c, int &u) {
  u = 0;
  const Val &o = c["set"];
  if (o.t == Val::NUL) return true;
  const Val *v = o.t == Val::ARR && o.a.size() == 1 ? &o.a[0] : &o;
  if (v->t != Val::STR || (v->s != "rural" && v->s != "urban")) {
    return ckstop(c, "`set` must be \"rural\" or \"urban\"");
  }
  u = v->s == "urban";
  return true;
}

inline bool adimg(Ctx &c, int &k) {
  Nv d;
  CK(adnumd(c, "nimg", 3, d, 0, false, true));
  if (d[0] > 1000) {
    return ckstop(c, "`nimg` must be at most 1,000: the image terms decay geometrically and more add only time");
  }
  k = static_cast<int>(d[0]);
  return true;
}

inline bool adlid(Ctx &c, double h, double &lid) {
  lid = NaN;
  if (!c.has("mixh")) return true;
  Nv d;
  CK(adnum(c, c["mixh"], "mixh", d, -HUGE_VAL, true));
  lid = d[0];
  if (h > lid) {
    c.warn.push_back(std::string(c.fn) + ": release height " + rfmt7(h) + " is above the mixing height " +
                     rfmt7(lid) + ": the plume is above the lid and the reflection terms do not describe it.");
  }
  return true;
}

// The receptors, n x 3 column-major: a matrix (rows), an object of three
// columns, an array of triples or a flat array of x, y, z in turn.
inline bool adrec(Ctx &c, std::vector<double> &r, int &n) {
  const Val &o = asmat(c["rcpt"]);
  r.clear();
  n = 0;
  auto bad = [&]() { return ckstop(c, "`rcpt` must give x, y, z for at least one receptor"); };
  if (o.t == Val::ARR && !o.a.empty() && o.a[0].t == Val::ARR) {
    std::vector<double> x, y, z;
    for (const Val &e : o.a) {
      Dv t;
      if (e.t != Val::ARR || readd(e, t) != NUM || t.n() != 3) return bad();
      x.push_back(t[0]);
      y.push_back(t[1]);
      z.push_back(t[2]);
    }
    r = x;
    r.insert(r.end(), y.begin(), y.end());
    r.insert(r.end(), z.begin(), z.end());
    n = static_cast<int>(x.size());
  } else {
    Dv t;
    if (readd(o, t) != NUM || t.n() == 0 || t.n() % 3 != 0) return bad();
    n = static_cast<int>(t.n() / 3);
    r.assign(3 * static_cast<std::size_t>(n), 0.0);
    for (int i = 0; i < n; ++i) {
      for (int k = 0; k < 3; ++k) r[static_cast<std::size_t>(k) * n + i] = t[3 * static_cast<std::size_t>(i) + k];
    }
  }
  for (double v : r) {
    if (std::isnan(v)) return bad();
  }
  return true;
}

inline void adhts(Ctx &c, const std::vector<double> &r, int n, double lid) {
  int lo = 0, hi = 0;
  for (int i = 0; i < n; ++i) {
    const double z = r[2 * static_cast<std::size_t>(n) + i];
    lo += z < 0;
    hi += !std::isnan(lid) && z > lid;
  }
  if (lo) {
    c.warn.push_back(std::string(c.fn) + ": " + std::to_string(lo) +
                     " receptor(s) below ground level (z < 0): the plume is reflected at the ground, so these "
                     "values are mirror images.");
  }
  if (hi) {
    c.warn.push_back(std::string(c.fn) + ": " + std::to_string(hi) + " receptor(s) above the mixing height " +
                     rfmt7(lid) +
                     ": the lid reflects the plume, so these values are not concentrations in the free "
                     "atmosphere.");
  }
}

bool pgsig(Ctx &c, Val &out) {
  int k = 0, u = 0;
  CK(adcls(c, k));
  CK(adset(c, u));
  Nv x;
  CK(adnum(c, c["x"], "x", x, 0, false, false));
  std::vector<double> sy(x.n()), sz(x.n());
  CK(kchk(c, clab_pgsig(x.p(), static_cast<int>(x.n()), k, u, sy.data(), sz.data())));
  out = vobj();
  put(out, "sigmay", vauto(sy));
  put(out, "sigmaz", vauto(sz));
  return true;
}

bool briggs(Ctx &c, Val &out) {
  int k = 0;
  CK(adcls(c, k));
  const double u = val1(c["u"]);
  if (std::isnan(u) || !std::isfinite(u) || u <= 0 || c["u"].t == Val::BOOL) {
    return ckstop(c, "`u` must be a single positive wind speed");
  }
  Nv x, di, ev, ts, ta, g, dz;
  CK(adnum(c, c["x"], "x", x, 0, false, false));
  CK(adnum(c, c["diam"], "diam", di, -HUGE_VAL, true));
  CK(adnum(c, c["exitvel"], "exitvel", ev, -HUGE_VAL, true));
  CK(adnum(c, c["stktmp"], "stktmp", ts, -HUGE_VAL, true));
  CK(adnum(c, c["ambtmp"], "ambtmp", ta, -HUGE_VAL, true));
  if (ts[0] < 150 || ta[0] < 150) {
    return ckstop(c,
                  "temperatures are in kelvin: stktmp %s and ambtmp %s look like degrees Celsius (ambient air is "
                  "about 290 K)",
                  rfmt7(ts[0]).c_str(), rfmt7(ta[0]).c_str());
  }
  CK(adnumd(c, "g", 9.80616, g, -HUGE_VAL, true));
  if (ts[0] < ta[0] - 0.5) {
    return ckstop(c,
                  "`stktmp` is below `ambtmp`: Briggs' formulas describe a buoyant plume, so the stack must be as "
                  "hot as the air or hotter (equal gives no rise)");
  }
  double lap = NaN;
  if (c.has("dthdz")) {
    CK(adnum(c, c["dthdz"], "dthdz", dz, -HUGE_VAL, true));
    lap = dz[0];
  }
  double o[3];
  std::vector<double> rise(x.n());
  CK(kchk(c, clab_briggs(x.p(), static_cast<int>(x.n()), u, di[0], ev[0], std::max(ts[0], ta[0]), ta[0], k, lap,
                         g[0], o, rise.data())));
  out = vobj();
  put(out, "flux", vnum(o[0]));
  put(out, "finalrise", vnum(o[1]));
  put(out, "xfinal", vnum(o[2]));
  put(out, "rise", vauto(rise));
  return true;
}

// gplume() and gspuff(): puff when t is given.
bool plume(Ctx &c, Val &out, bool puff) {
  Nv q, u, h, t;
  CK(adnum(c, c[puff ? "mass" : "q"], puff ? "mass" : "q", q, 0));
  CK(adnum(c, c["u"], "u", u, -HUGE_VAL, true));
  CK(adnum(c, c["h"], "h", h, 0));
  if (puff) CK(adnum(c, c["t"], "t", t, -HUGE_VAL, true));
  int ni = 0;
  CK(adimg(c, ni));
  double lid = NaN;
  CK(adlid(c, h[0], lid));
  std::vector<double> r;
  int n = 0;
  CK(adrec(c, r, n));
  adhts(c, r, n, lid);
  int k = 0, ur = 0;
  CK(adcls(c, k));
  CK(adset(c, ur));
  std::vector<double> v(static_cast<std::size_t>(n));
  if (puff) {
    CK(kchk(c, clab_gspuff(q[0], u[0], h[0], r.data(), n, t[0], k, ur, lid, ni, v.data())));
  } else {
    CK(kchk(c, clab_gplume(q[0], u[0], h[0], r.data(), n, k, ur, lid, ni, v.data())));
  }
  out = vauto(v);
  return true;
}

bool gplume(Ctx &c, Val &out) { return plume(c, out, false); }
bool gspuff(Ctx &c, Val &out) { return plume(c, out, true); }

bool plmfld(Ctx &c, Val &out) {
  Nv x, y, z;
  CK(adnum(c, c["x"], "x", x, -HUGE_VAL, false, false));
  CK(adnum(c, c["y"], "y", y, -HUGE_VAL, false, false));
  bool up = true;
  for (std::size_t i = 1; i < x.n(); ++i) up = up && x[i] > x[i - 1];
  for (std::size_t i = 1; i < y.n(); ++i) up = up && y[i] > y[i - 1];
  if (!up) return ckstop(c, "`x` and `y` must be strictly increasing receptor coordinates");
  CK(adnumd(c, "z", 0, z, 0));
  const std::size_t nx = x.n(), ny = y.n(), n = nx * ny;
  Val rc = varr();
  std::vector<double> flat;
  for (std::size_t j = 0; j < ny; ++j) {
    for (std::size_t i = 0; i < nx; ++i) {
      flat.push_back(x[i]);
      flat.push_back(y[j]);
      flat.push_back(z[0]);
    }
  }
  // gplume() on the grid, its checks in its own name
  Val in = *c.in;
  bool had = false;
  for (std::size_t i = 0; i < in.k.size(); ++i) {
    if (in.k[i] == "rcpt") had = true;
  }
  if (!had) put(in, "rcpt", vnums(flat));
  Ctx g{"gplume", &in, std::string(), {}};
  Val v;
  if (!plume(g, v, false)) {
    c.err = g.err;
    return false;
  }
  for (auto &w : g.warn) c.warn.push_back(w);
  std::vector<double> f(n);
  if (v.t == Val::NUM) f[0] = v.n;
  else for (std::size_t i = 0; i < n; ++i) f[i] = v.a[i].n;
  auto med = [](std::vector<double> d) {
    std::sort(d.begin(), d.end());
    const std::size_t m = d.size();
    return m % 2 ? d[m / 2] : (d[m / 2 - 1] + d[m / 2]) / 2;
  };
  std::vector<double> ddx, ddy;
  for (std::size_t i = 1; i < nx; ++i) ddx.push_back(x[i] - x[i - 1]);
  for (std::size_t i = 1; i < ny; ++i) ddy.push_back(y[i] - y[i - 1]);
  std::string st = "D";
  if (c.has("stab")) lab(c["stab"].t == Val::ARR ? c["stab"].a[0] : c["stab"], st);
  Nv h;
  CK(adnum(c, c["h"], "h", h, 0));
  out = vobj();
  put(out, "field", vmat(f, static_cast<int>(nx), static_cast<int>(ny)));
  put(out, "x", vauto(x.v));
  put(out, "y", vauto(y.v));
  put(out, "dx", ddx.empty() ? Val() : vnum(med(ddx)));
  put(out, "dy", ddy.empty() ? Val() : vnum(med(ddy)));
  put(out, "units", vstr("concentration"));
  put(out, "title", vstr("Gaussian plume, class " + st + ", h = " + rfmt7(h[0]) + " m"));
  return true;
}

// format(x, big.mark = ",", scientific = FALSE) of a whole number.
inline std::string bigm(double x) {
  char b[64];
  std::snprintf(b, sizeof b, "%.0f", std::fabs(x));
  std::string s = b, r;
  const int n = static_cast<int>(s.size());
  for (int i = 0; i < n; ++i) {
    r += s[static_cast<std::size_t>(i)];
    if ((n - i - 1) % 3 == 0 && i < n - 1) r += ',';
  }
  return x < 0 ? "-" + r : r;
}

bool advdif2d(Ctx &c, Val &out) {
  int nx = 0, ny = 0;
  std::vector<double> c0;
  {
    int d[2] = {-1, -1};
    const Val &o = c["c0"];
    if (o.t != Val::ARR || o.a.empty() || o.a[0].t != Val::ARR || !shape(o, 2, d, 0)) {
      return ckstop(c, "`c0` must be a numeric matrix (one row per x cell, one column per y cell)");
    }
    c0.assign(static_cast<std::size_t>(d[0]) * static_cast<std::size_t>(d[1]), 0.0);
    bool bad = false;
    fill(o, 2, d, 0, 0, 1, c0, bad, true);
    nx = d[0];
    ny = d[1];
  }
  for (double v : c0) {
    if (!std::isfinite(v)) return ckstop(c, "`c0` must be a finite numeric matrix");
  }
  for (double v : c0) {
    if (v < 0) return ckstop(c, "`c0` must be non-negative: it is a concentration field");
  }
  std::vector<double> src;
  if (c.has("source")) {
    int d[2] = {-1, -1};
    const Val &o = c["source"];
    if (!shape(o, 2, d, 0) || d[0] != nx || d[1] != ny) return ckstop(c, "`source` must have the dimensions of `c0`");
    src.assign(c0.size(), 0.0);
    bool bad = false;
    fill(o, 2, d, 0, 0, 1, src, bad, true);
    for (double v : src) {
      if (!std::isfinite(v)) return ckstop(c, "`source` must be a finite numeric matrix");
    }
    for (double v : src) {
      if (v < 0) return ckstop(c, "`source` must be non-negative: it is an emission rate, like `c0`");
    }
  }
  Nv u, v, kx, ky, dx, dy, dt, ns, mw;
  CK(adnum(c, c["u"], "u", u));
  CK(adnum(c, c["v"], "v", v));
  CK(adnum(c, c["kx"], "kx", kx, 0));
  CK(adnum(c, c["ky"], "ky", ky, 0));
  CK(adnum(c, c["dx"], "dx", dx, -HUGE_VAL, true));
  CK(adnum(c, c["dy"], "dy", dy, -HUGE_VAL, true));
  CK(adnum(c, c["dt"], "dt", dt, -HUGE_VAL, true));
  CK(adnum(c, c["nsteps"], "nsteps", ns, 0, false, true, true));
  auto nb = [&](double t, double *o) {
    o[0] = u[0] * t / dx[0];
    o[1] = v[0] * t / dy[0];
    o[2] = kx[0] * t / (dx[0] * dx[0]);
    o[3] = ky[0] * t / (dy[0] * dy[0]);
    o[4] = std::fabs(o[0]) + std::fabs(o[1]) + 2 * (o[2] + o[3]);
  };
  double b[5];
  double t = dt[0], nst = ns[0];
  nb(t, b);
  if (!std::isfinite(b[4])) return ckstop(c, "the stability number is not finite: check u, v, kx, ky, dx, dy and dt");
  const double sub = b[4] > 1 ? std::ceil(b[4]) : 1;
  const double work = static_cast<double>(nx) * ny * nst * sub;
  if (!std::isfinite(work)) return ckstop(c, "the amount of work is not finite: check nsteps and the grid");
  double maxw = 1e8;
  if (c.has("maxwork")) {
    const double m = val1(c["maxwork"]);
    const bool inf = c["maxwork"].t == Val::STR && c["maxwork"].s == "Inf";
    if ((std::isnan(m) && !inf) || (!inf && m <= 0)) {
      return ckstop(c, "`maxwork` must be a single positive number (Inf lifts the cap)");
    }
    maxw = inf ? HUGE_VAL : m;
  }
  if (work > maxw) {
    char sb[64];
    std::snprintf(sb, sizeof sb, " from a stability number of %.4g", b[4]);
    return ckstop(c,
                  "this run updates %s cells x steps x sub-steps (%d x %d cells, %s steps, %s sub-steps per step%s), "
                  "above maxwork = %s. Pass maxwork = Inf to run it anyway, or use a coarser grid, fewer steps or a "
                  "shorter total time",
                  bigm(work).c_str(), nx, ny, bigm(nst).c_str(), bigm(sub).c_str(), sub > 1 ? sb : "",
                  bigm(maxw).c_str());
  }
  if (nst > INT_MAX) return ckstop(c, "`nsteps` must be at most 2^31 - 1");
  if (b[4] > 1) {
    char w[512];
    std::snprintf(w, sizeof w,
                  "advdif2d: explicit scheme outside its stability bound (|u dt/dx| + |v dt/dy| + 2 (kx dt/dx^2 + ky "
                  "dt/dy^2) = %.3f > 1): each step is taken as %d sub-steps of dt/%d; pass a smaller dt to silence "
                  "this.",
                  b[4], static_cast<int>(sub), static_cast<int>(sub));
    c.warn.push_back(w);
    t = t / sub;
    nb(t, b);
    nst = nst * sub;
  }
  CK(kchk(c, clab_advdif(c0.data(), nx, ny, u[0], v[0], kx[0], ky[0], dx[0], dy[0], t, nst,
                         src.empty() ? nullptr : src.data())));
  long double s = 0;
  for (double e : c0) s += e;
  std::vector<double> xa, ya;
  for (int i = 1; i <= nx; ++i) xa.push_back(i * dx[0]);
  for (int j = 1; j <= ny; ++j) ya.push_back(j * dy[0]);
  out = vobj();
  put(out, "field", vmat(c0, nx, ny));
  put(out, "mass", vnum(static_cast<double>(s) * dx[0] * dy[0]));
  put(out, "cfl", vnum(std::fabs(b[0]) + std::fabs(b[1])));
  put(out, "difn", vnum(2 * (b[2] + b[3])));
  put(out, "stabn", vnum(b[4]));
  put(out, "stable", vbool(b[4] <= 1));
  put(out, "dt", vnum(t));
  put(out, "substeps", vnum(sub));
  put(out, "x", vauto(xa));
  put(out, "y", vauto(ya));
  put(out, "dx", vnum(dx[0]));
  put(out, "dy", vnum(dy[0]));
  put(out, "units", vstr("concentration"));
  put(out, "title", vstr("Advection-diffusion after " + bigm(nst) + " steps of dt = " + rfmt7(t)));
  return true;
}

bool lagrprt(Ctx &c, Val &out) {
  double seed = 0;
  if (c.has("seed")) {
    seed = val1(c["seed"]);
    if (c["seed"].t == Val::BOOL || std::isnan(seed) || !std::isfinite(seed) || seed != std::nearbyint(seed) ||
        seed < -2147483647.0 || seed >= 9007199254740992.0) {
      return ckstop(c, "`seed` must be a single whole number from -2147483647 to 2^53 - 1");
    }
  }
  Nv n, x0, y0, u, v, kx, ky, dt, ns;
  CK(adnum(c, c["n"], "n", n, -HUGE_VAL, true, true, true));
  if (n[0] > 1e8) return ckstop(c, "`n` must be at most 100 million particles");
  CK(adnum(c, c["x0"], "x0", x0));
  CK(adnum(c, c["y0"], "y0", y0));
  CK(adnum(c, c["u"], "u", u));
  CK(adnum(c, c["v"], "v", v));
  CK(adnum(c, c["kx"], "kx", kx, 0));
  CK(adnum(c, c["ky"], "ky", ky, 0));
  CK(adnum(c, c["dt"], "dt", dt, -HUGE_VAL, true));
  CK(adnum(c, c["nsteps"], "nsteps", ns, 0, false, true, true));
  double maxw = 1e8;
  if (c.has("maxwork")) {
    const double m = val1(c["maxwork"]);
    const bool inf = c["maxwork"].t == Val::STR && c["maxwork"].s == "Inf";
    if ((std::isnan(m) && !inf) || (!inf && m <= 0)) {
      return ckstop(c, "`maxwork` must be a single positive number (Inf lifts the cap)");
    }
    maxw = inf ? HUGE_VAL : m;
  }
  if (n[0] * ns[0] > maxw) {
    return ckstop(c, "this run moves %s particles x %s steps, above maxwork = %s; pass maxwork = Inf to run it anyway",
                  bigm(n[0]).c_str(), bigm(ns[0]).c_str(), bigm(maxw).c_str());
  }
  if (ns[0] > INT_MAX) return ckstop(c, "`nsteps` must be at most 2^31 - 1");
  Nv g;
  const bool hg = c.has("grid");
  if (hg) {
    if (vlen(c["grid"]) != 6) return ckstop(c, "`grid` must be c(xmin, xmax, ymin, ymax, nx, ny)");
    CK(adnum(c, c["grid"], "grid", g, -HUGE_VAL, false, false));
    if (g[1] <= g[0] || g[3] <= g[2] || g[4] < 1 || g[5] < 1 || g[4] != std::nearbyint(g[4]) ||
        g[5] != std::nearbyint(g[5])) {
      return ckstop(c, "`grid` needs xmin < xmax, ymin < ymax and whole positive nx, ny");
    }
    if (g[4] * g[5] > maxw) {
      return ckstop(c, "a %s x %s grid has %s cells, above maxwork = %s; pass maxwork = Inf to allocate it anyway",
                    bigm(g[4]).c_str(), bigm(g[5]).c_str(), bigm(g[4] * g[5]).c_str(), bigm(maxw).c_str());
    }
  }
  const int m = static_cast<int>(n[0]);
  std::vector<double> xs(static_cast<std::size_t>(m)), ys(xs);
  CK(kchk(c, clab_lagr(m, x0[0], y0[0], u[0], v[0], kx[0], ky[0], dt[0], ns[0], seed, xs.data(), ys.data())));
  double mean[2];
  std::vector<double> conc(hg ? static_cast<std::size_t>(g[4] * g[5]) : 0);
  CK(kchk(c, clab_laggrid(xs.data(), ys.data(), m, hg ? g.p() : nullptr, conc.data(), mean)));
  out = vobj();
  put(out, "x", vauto(xs));
  put(out, "y", vauto(ys));
  put(out, "mean_x", vnum(mean[0]));
  put(out, "meany", vnum(mean[1]));
  put(out, "grid", hg ? c["grid"] : Val());
  if (hg) put(out, "concentration", vmat(conc, static_cast<int>(g[4]), static_cast<int>(g[5])));
  return true;
}

// .polsc(): one finite number.
inline bool polsc(Ctx &c, const Val &o, const char *w, double &v, bool nonneg = true, bool positive = false) {
  bool lg = false;
  const char *cl = o.t == Val::NUL ? nullptr : jcls(o, lg);
  if (lg) {
    bool allna = o.t != Val::BOOL;
    for (const Val &e : o.a) allna = allna && e.t == Val::NUL;
    if (!allna) return ckstop(c, "`%s` must be numeric, not logical", w);
  }
  if (o.t == Val::NUL) return ckstop(c, "`%s` must be numeric, got NULL", w);
  if (cl != nullptr && !lg) return ckstop(c, "`%s` must be numeric, got %s", w, cl);
  Dv t;
  readd(o, t);
  if (t.n() != 1 || !std::isfinite(t[0])) return ckstop(c, "`%s` must be a single finite number", w);
  v = t[0];
  if (positive && v <= 0) return ckstop(c, "`%s` must be positive", w);
  if (nonneg && v < 0) return ckstop(c, "`%s` must be non-negative", w);
  return true;
}

inline bool polscd(Ctx &c, const char *nm, double def, double &v, bool nonneg = true) {
  if (!c.has(nm)) {
    v = def;
    return true;
  }
  return polsc(c, c[nm], nm, v, nonneg);
}

inline bool polexp(Ctx &c, Nv &z) {
  const Val &o = c["expo"];
  bool lg = false;
  const char *cl = jcls(o, lg);
  if (lg) return ckstop(c, "`expo` must be numeric, not logical");
  if (cl != nullptr) return ckstop(c, "`expo` must be numeric, got %s", cl);
  Dv t;
  readd(o, t);
  z.v = t.v;
  bool fin = z.n() > 0;
  for (double x : z.v) fin = fin && std::isfinite(x);
  if (!fin) return ckstop(c, "`expo` must be one or more finite concentrations, with no missing values");
  for (double x : z.v) {
    if (x < 0) return ckstop(c, "`expo` must be non-negative");
  }
  return true;
}

inline bool crfout(Ctx &c, const std::string &y, const std::vector<std::string> &av, const char *tail) {
  if (std::find(av.begin(), av.end(), y) != av.end()) return true;
  std::string a;
  for (std::size_t i = 0; i < av.size(); ++i) a += (i ? ", " : "") + av[i];
  return ckstop(c, "Unknown outcome '%s'. Available: %s%s", y.c_str(), a.c_str(), tail);
}

// The outcome key: one string; anything else is shown as R pastes it.
inline bool ykey(Ctx &c, std::string &y, bool &one) {
  y = "allcause";
  one = true;
  if (!c.has("y")) return true;
  const Val &o = c["y"];
  const Val *s = o.t == Val::ARR && o.a.size() == 1 ? &o.a[0] : &o;
  if (s->t == Val::STR) {
    y = s->s;
    return true;
  }
  one = false;
  std::vector<std::string> l;
  labs(o, l);
  y.clear();
  for (std::size_t i = 0; i < l.size(); ++i) y += (i ? "," : "") + l[i];
  return true;
}

// The outcome of the pollutant's CRF: the PM2.5 or the NO2 table.
inline bool crfchk(Ctx &c, bool pm, const std::string &y, bool one) {
  static const std::vector<std::string> a = {"allcause", "ihd", "stroke"}, b = {"allcause", "respiratory", "chiast"};
  const char *tail = pm ? "" : " or pass betaper10 explicitly";
  if (!one) return crfout(c, y + std::string(1, '\0'), pm ? a : b, tail);
  return crfout(c, y, pm ? a : b, tail);
}

inline double rmean(const std::vector<double> &v) { return nmean(v.data(), v.size()); }

// crfpm25() and crfno2() as one: the CRF object, from z.
bool crfv(Ctx &c, bool pm, const Nv &z, const std::string &y, double refc, bool hasb, double bpar, Val &out) {
  std::vector<double> rr(z.n()), lr(z.n()), lo(z.n(), NA), up(z.n(), NA);
  double par[3] = {0, 0, 0}, ci[2] = {NA, NA};
  int form = 0;
  Val ex = vobj();
  std::string cite;
  if (pm && y == "allcause") {
    par[0] = std::log(1.08);
    ci[0] = std::log(1.06);
    ci[1] = std::log(1.09);
  } else if (pm) {
    form = 1;
    const bool ihd = y == "ihd";
    par[0] = ihd ? 1.91 : 1.46;
    par[1] = ihd ? 0.14 : 0.13;
    par[2] = ihd ? 0.49 : 0.61;
  } else if (hasb) {
    par[0] = bpar;
  } else if (y == "allcause") {
    par[0] = std::log(1.02);
    ci[0] = std::log(1.01);
    ci[1] = std::log(1.04);
  } else if (y == "respiratory") {
    par[0] = std::log(1.03);
    ci[0] = std::log(1.00);
    ci[1] = std::log(1.05);
  } else {
    par[0] = std::log(1.05) * 10 / 4;
    ci[0] = std::log(1.02) * 10 / 4;
    ci[1] = std::log(1.07) * 10 / 4;
  }
  const int n = static_cast<int>(z.n());
  CK(kchk(c, clab_crf(z.p(), n, refc, form, par, rr.data(), lr.data())));
  if (form == 0 && !std::isnan(ci[0])) {
    std::vector<double> t(z.n());
    CK(kchk(c, clab_crf(z.p(), n, refc, 0, &ci[0], lo.data(), t.data())));
    CK(kchk(c, clab_crf(z.p(), n, refc, 0, &ci[1], up.data(), t.data())));
  }
  const bool one = n == 1;
  auto m = [&](const std::vector<double> &v) {
    if (one) return v[0];
    bool any = false;
    for (double e : v) any = any || std::isnan(e);
    return any ? NA : rmean(v);
  };
  if (form == 0) {
    put(ex, "betaper10", vnum(par[0]));
    put(ex, "betaper10ci", vnums(std::vector<double>{na(ci[0]), na(ci[1])}));
  } else {
    put(ex, "alpha", vnum(par[0]));
    put(ex, "gamma", vnum(par[1]));
    put(ex, "delta", vnum(par[2]));
  }
  if (pm) put(ex, "form", vstr(form ? "IER" : "log-linear"));
  put(ex, "outcome", vstr(y));
  put(ex, "rrperunit", one ? Val() : vnums(rr));
  if (pm) {
    cite = form ? "Burnett et al. (2014) EHP 122(4):397-403"
                : "Chen & Hoek (2020) Environ Int 143:105974; WHO (2021) Global AQ Guidelines";
  } else {
    cite = y == "chiast" && !hasb ? "Khreis et al. (2017) Environ Int 100:1-31; WHO (2021) Global AQ Guidelines"
                                  : "Huangfu & Atkinson (2020) Environ Int 144:105998; WHO (2021) Global AQ Guidelines";
  }
  out = vobj();
  put(out, "rr", vnum(m(rr)));
  put(out, "logrr", vnum(m(lr)));
  put(out, "rrlo", vnum(na(m(lo))));
  put(out, "rrup", vnum(na(m(up))));
  put(out, "refconc", vnum(refc));
  put(out, "expconc", vnum(m(z.v)));
  put(out, "pollutant", vstr(pm ? "PM2.5" : "NO2"));
  put(out, "citation", vstr(cite));
  put(out, "extra", ex);
  return true;
}

bool crfpm25(Ctx &c, Val &out) {
  std::string y;
  bool one = true;
  ykey(c, y, one);
  CK(crfchk(c, true, y, one));
  double refc = 0;
  CK(polscd(c, "refconc", 5.8, refc));
  Nv z;
  CK(polexp(c, z));
  return crfv(c, true, z, y, refc, false, 0, out);
}

bool crfno2(Ctx &c, Val &out) {
  std::string y;
  bool one = true;
  ykey(c, y, one);
  double b = 0;
  const bool hb = c.has("betaper10");
  if (!hb) {
    CK(crfchk(c, false, y, one));
  } else {
    CK(polsc(c, c["betaper10"], "betaper10", b, false));
  }
  double refc = 0;
  CK(polscd(c, "refconc", 10, refc));
  Nv z;
  CK(polexp(c, z));
  return crfv(c, false, z, y, refc, hb, b, out);
}

bool attrfrac(Ctx &c, Val &out) {
  double rr = 0, p = 0;
  CK(polsc(c, c["rr"], "rr", rr, false, true));
  CK(polsc(c, c["expoprev"], "expoprev", p));
  if (p > 1) return ckstop(c, "expoprev must be in [0,1]; got %s", rfmt7(p).c_str());
  if (c.has("rrlo") || c.has("rrup")) {
    if (!c.has("rrlo") || !c.has("rrup")) return ckstop(c, "rrlo and rrup go together");
    double lo = 0, hi = 0;
    CK(polsc(c, c["rrlo"], "rrlo", lo, false, true));
    CK(polsc(c, c["rrup"], "rrup", hi, false, true));
    if (lo > rr || hi < rr) return ckstop(c, "rr must lie between rrlo and rrup");
  }
  out = vnum(clab_levin(rr, p));
  return true;
}

inline bool flag01(Ctx &c, const char *nm, bool &b) {
  b = false;
  if (!c.has(nm)) return true;
  const Val &o = c[nm];
  const Val *s = o.t == Val::ARR && o.a.size() == 1 ? &o.a[0] : &o;
  if (s->t != Val::BOOL) return ckstop(c, "`%s` must be TRUE or FALSE", nm);
  b = s->b;
  return true;
}

inline bool polrate(Ctx &c, double r, bool rare, bool hi) {
  if (r > 1) {
    return ckstop(c,
                  "`brate` is a rate per person-year, between 0 and 1 (0.008 is 800 per 100,000); %s looks like a "
                  "rate per 100,000",
                  rfmt7(r).c_str());
  }
  if (r == 0) return true;
  if (r < 1e-11) {
    return ckstop(c, "`brate` %s per person-year is below one case per 100 billion person-years: not a usable rate",
                  rfmt7(r).c_str());
  }
  if (r < 1e-5 && !rare) {
    return ckstop(c,
                  "`brate` %s per person-year is below one case per 100,000 a year. A rate per 100,000 divided by "
                  "100,000 is the rate here (800 per 100,000 is 0.008); if the outcome really is this rare, say so "
                  "with rare = TRUE",
                  rfmt7(r).c_str());
  }
  if (r > 0.1 && !hi) {
    char b[64];
    std::snprintf(b, sizeof b, "%s", rfmt7(r * 1e5).c_str());
    std::string s = b;
    // format(big.mark = ",", scientific = FALSE)
    const std::size_t dot = s.find('.');
    std::string ip = s.substr(0, dot), fp = dot == std::string::npos ? "" : s.substr(dot);
    std::string w;
    for (std::size_t i = 0; i < ip.size(); ++i) {
      w += ip[i];
      if ((ip.size() - i - 1) % 3 == 0 && i + 1 < ip.size()) w += ',';
    }
    return ckstop(c,
                  "`brate` %s per person-year is %s per 100,000 a year, above a whole population's all-cause death "
                  "rate (an age-specific rate for the oldest ages can be this high). If it is a rate per 100,000, "
                  "divide it by 100,000; if the outcome really is this common, say so with hirate = TRUE",
                  rfmt7(r).c_str(), (w + fp).c_str());
  }
  return true;
}

bool mortdisp(Ctx &c, Val &out) {
  double dx = 0, pop = 0, r = 0, b = 0;
  bool ra = false, hi = false;
  CK(polsc(c, c["expdel"], "expdel", dx, false));
  CK(polsc(c, c["pop"], "pop", pop));
  CK(polsc(c, c["brate"], "brate", r));
  CK(flag01(c, "rare", ra));
  CK(flag01(c, "hirate", hi));
  CK(polrate(c, r, ra, hi));
  CK(polsc(c, c["betaperunit"], "betaperunit", b, false));
  out = vnum(clab_hif(dx, pop, r, b));
  return true;
}

inline std::string polkey(const Val &o) {
  const Val *s = o.t == Val::ARR && o.a.size() == 1 ? &o.a[0] : &o;
  if (o.t == Val::NUL) return "NO2";
  if (s->t != Val::STR) return "";
  std::string k;
  for (char ch : s->s) {
    const char l = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    if ((l >= 'a' && l <= 'z') || (l >= '0' && l <= '9')) k += l;
  }
  return k == "pm25" ? "PM2.5" : k == "no2" ? "NO2" : "";
}

// polbrdn() on given values: its checks, then the burden.
bool brdn(Ctx &c, const Val &pl, double em, double ep, double br, double pop, const std::string &y, bool yone,
          const Val &refv, bool ra, bool hi, Val &out) {
  const std::string pol = polkey(pl);
  if (pol.empty()) {
    std::vector<std::string> l;
    labs(pl, l);
    std::string a;
    for (std::size_t i = 0; i < l.size(); ++i) a += (i ? "," : "") + l[i];
    return ckstop(c, "Unknown pollutant '%s'. Use 'PM2.5' or 'NO2'", a.c_str());
  }
  CK(polrate(c, br, ra, hi));
  double refc = pol == "PM2.5" ? 5.8 : 10;
  if (refv.t != Val::NUL) CK(polsc(c, refv, "refconc", refc));
  Nv z;
  z.v = {em};
  Val crf;
  Ctx cc{pol == "PM2.5" ? "crfpm25" : "crfno2", c.in, std::string(), {}};
  if (!crfchk(cc, pol == "PM2.5", y, yone)) {
    c.err = cc.err;
    return false;
  }
  CK(crfv(c, pol == "PM2.5", z, y, refc, false, 0, crf));
  if (ep > 1) {
    Ctx ca{"attrfrac", c.in, std::string(), {}};
    ckstop(ca, "expoprev must be in [0,1]; got %s", rfmt7(ep).c_str());
    c.err = ca.err;
    return false;
  }
  const double rr = col(crf, "rr").n, lrr = col(crf, "logrr").n;
  const double paf = clab_levin(rr, ep), base = br * pop;
  Val ex = vobj();
  put(ex, "rr", vnum(rr));
  put(ex, "logrr", vnum(lrr));
  put(ex, "outcome", vstr(y));
  put(ex, "unit", vstr(y == "chiast" ? "incident cases" : "deaths"));
  out = vobj();
  put(out, "paf", vnum(paf));
  put(out, "attrcases", vnum(paf * base));
  put(out, "basecases", vnum(base));
  put(out, "population", vnum(pop));
  put(out, "brate", vnum(br));
  put(out, "expomean", vnum(em));
  put(out, "refconc", vnum(refc));
  put(out, "pollutant", vstr(pol));
  put(out, "citation", vstr(col(crf, "citation").s + "; Rothman et al. (2008) section 5"));
  put(out, "extra", ex);
  return true;
}

bool polbrdn(Ctx &c, Val &out) {
  if (polkey(c["poll"]).empty()) return brdn(c, c["poll"], 0, 0, 0, 0, "", true, Val(), false, false, out);
  double em = 0, ep = 0, br = 0, pop = 0;
  bool ra = false, hi = false;
  CK(polsc(c, c["expomean"], "expomean", em));
  CK(polsc(c, c["expoprev"], "expoprev", ep));
  CK(polsc(c, c["brate"], "brate", br));
  CK(flag01(c, "rare", ra));
  CK(flag01(c, "hirate", hi));
  CK(polrate(c, br, ra, hi));
  CK(polsc(c, c["pop"], "pop", pop));
  std::string y;
  bool one = true;
  ykey(c, y, one);
  return brdn(c, c["poll"], em, ep, br, pop, y, one, c["refconc"], ra, hi, out);
}

inline bool colnm(Ctx &c, const char *nm, const char *def, std::string &s) {
  s = def;
  if (!c.has(nm)) return true;
  const Val &o = c[nm];
  const Val *v = o.t == Val::ARR && o.a.size() == 1 ? &o.a[0] : &o;
  if (v->t != Val::STR || v->s.empty()) return ckstop(c, "`%s` must name one column of `data`", nm);
  s = v->s;
  return true;
}

bool polbrdarea(Ctx &c, Val &out) {
  std::string ac, ec, pc, bc;
  CK(colnm(c, "areacol", "fsa", ac));
  CK(colnm(c, "expcol", "exposure", ec));
  CK(colnm(c, "popcol", "population", pc));
  CK(colnm(c, "basratecol", "brate", bc));
  const Val &t = c["areatab"];
  std::string miss;
  for (const std::string &k : {ac, ec, pc, bc}) {
    if (t.t != Val::OBJ || col(t, k.c_str()).t == Val::NUL) miss += (miss.empty() ? "" : ", ") + k;
  }
  if (!miss.empty()) return ckstop(c, "areatab missing columns: %s", miss.c_str());
  const std::size_t n = vlen(col(t, ac.c_str()));
  if (n == 0) return ckstop(c, "areatab has no rows");
  std::vector<std::vector<double>> v(3);
  int q = 0;
  for (const std::string &k : {ec, pc, bc}) {
    const Val &cl = col(t, k.c_str());
    if (!isnum(cl)) return ckstop(c, "column `%s` must be numeric (a factor or character column is not a count)", k.c_str());
    v[static_cast<std::size_t>(q++)] = dbl(cl);
  }
  bool ra = false, hi = false;
  CK(flag01(c, "rare", ra));
  CK(flag01(c, "hirate", hi));
  std::string y;
  bool one = true;
  ykey(c, y, one);
  std::vector<double> rr(n), pf(n), at(n), bs(n);
  for (std::size_t i = 0; i < n; ++i) {
    Ctx r{"polbrdn", c.in, std::string(), {}};
    Val b;
    Val e = vnum(v[0][i]), p = vnum(v[1][i]), br = vnum(v[2][i]);
    double em = 0, pop = 0, rt = 0;
    bool ok = polsc(r, e, "expomean", em) && polsc(r, br, "brate", rt) && polrate(r, rt, ra, hi) &&
              polsc(r, p, "pop", pop) && brdn(r, c["poll"], em, 1, rt, pop, y, one, Val(), ra, hi, b);
    if (polkey(c["poll"]).empty()) ok = brdn(r, c["poll"], 0, 0, 0, 0, y, one, Val(), ra, hi, b);
    if (!ok) {
      std::string m = r.err;
      for (const char *pre : {"polbrdn: ", "crfpm25: ", "crfno2: "}) {
        if (m.rfind(pre, 0) == 0) m = m.substr(std::strlen(pre));
      }
      if (!m.empty() && m.back() == '.') m.pop_back();
      const std::pair<const char *, const char *> mp[] = {
          {"`expomean`", "`expcol`"}, {"`expo`", "`expcol`"}, {"`pop`", "`popcol`"}, {"`brate`", "`basratecol`"}};
      for (const auto &kv : mp) {
        for (std::size_t at2 = m.find(kv.first); at2 != std::string::npos; at2 = m.find(kv.first, at2)) {
          m.replace(at2, std::strlen(kv.first), kv.second);
          at2 += std::strlen(kv.second);
        }
      }
      return ckstop(c, "%s", m.c_str());
    }
    rr[i] = col(col(b, "extra"), "rr").n;
    pf[i] = col(b, "paf").n;
    at[i] = col(b, "attrcases").n;
    bs[i] = col(b, "basecases").n;
  }
  std::vector<std::size_t> o(n);
  for (std::size_t i = 0; i < n; ++i) o[i] = i;
  std::stable_sort(o.begin(), o.end(), [&](std::size_t a, std::size_t b) { return -at[a] < -at[b]; });
  auto pick = [&](const std::vector<double> &x) {
    std::vector<double> r;
    for (std::size_t i : o) r.push_back(x[i]);
    return vnums(r);
  };
  Val lab = varr();
  const Val &al = col(t, ac.c_str());
  for (std::size_t i : o) lab.a.push_back(al.t == Val::ARR ? al.a[i] : al);
  out = vobj();
  put(out, ac.c_str(), lab);
  put(out, ec.c_str(), pick(v[0]));
  put(out, pc.c_str(), pick(v[1]));
  put(out, bc.c_str(), pick(v[2]));
  put(out, "rr", pick(rr));
  put(out, "paf", pick(pf));
  put(out, "attrcases", pick(at));
  put(out, "basecases", pick(bs));
  return true;
}

bool expconcidx(Ctx &c, Val &out) {
  std::string ex, inc;
  if (!c.has("expo")) return ckstop(c, "`expo` must name one column of `data`");
  CK(colnm(c, "expo", "", ex));
  if (!c.has("income")) return ckstop(c, "`income` must name one column of `data`");
  CK(colnm(c, "income", "", inc));
  const Val &d = c["data"];
  std::string miss;
  for (const std::string &k : {ex, inc}) {
    if (d.t != Val::OBJ || col(d, k.c_str()).t == Val::NUL) miss += (miss.empty() ? "" : ", ") + k;
  }
  if (!miss.empty()) return ckstop(c, "data is missing columns: %s", miss.c_str());
  const Val &hc = col(d, ex.c_str()), &ic = col(d, inc.c_str());
  // complete rows: null (NA) in either column drops the row
  Val h2 = varr(), i2 = varr();
  const std::size_t n0 = std::max(vlen(hc), vlen(ic));
  for (std::size_t i = 0; i < n0; ++i) {
    const Val &a = hc.t == Val::ARR ? hc.a[i % hc.a.size()] : hc, &b = ic.t == Val::ARR ? ic.a[i % ic.a.size()] : ic;
    if (a.t == Val::NUL || b.t == Val::NUL) continue;
    h2.a.push_back(a);
    i2.a.push_back(b);
  }
  if (h2.a.size() < 2) return ckstop(c, "Need at least 2 rows to compute concentration index");
  Nv h, iv;
  CK(cknum(c, h2, ex.c_str(), h, false));
  CK(cknum(c, i2, inc.c_str(), iv, false));
  for (double x : h.v) {
    if (x < 0) return ckstop(c, "`%s` must be non-negative: an exposure is a concentration", ex.c_str());
  }
  std::vector<double> u = iv.v;
  std::sort(u.begin(), u.end());
  if (std::unique(u.begin(), u.end()) - u.begin() < 2) {
    return ckstop(c,
                  "`%s` must take at least two values: with one income level there is no income gradient and the "
                  "index is undefined",
                  inc.c_str());
  }
  if (!(rmean(h.v) > 0)) return ckstop(c, "Mean exposure is zero; concentration index undefined");
  const std::size_t n = h.n();
  double o[3];
  std::vector<double> cp(n), ce(n);
  CK(kchk(c, clab_concidx(h.p(), iv.p(), static_cast<int>(n), o, cp.data(), ce.data())));
  char b[256];
  if (o[0] < -0.01) {
    std::snprintf(b, sizeof b,
                  "CI = %.4f. Pro-poor exposure burden: lower-income individuals bear disproportionately higher "
                  "pollution.",
                  o[0]);
  } else if (o[0] > 0.01) {
    std::snprintf(b, sizeof b,
                  "CI = %.4f. Pro-rich exposure burden: higher-income individuals bear disproportionately higher "
                  "pollution.",
                  o[0]);
  } else {
    std::snprintf(b, sizeof b, "CI = %.4f. Exposure distributed approximately evenly across income.", o[0]);
  }
  Frame cv;
  cv.num("population", cp);
  cv.num("exposure", ce);
  Val e = vobj();
  put(e, "n", vnum(static_cast<double>(n)));
  put(e, "cov_h_R", vnum(o[1]));
  put(e, "curve", cv.o);
  out = vobj();
  put(out, "concidx", vnum(o[0]));
  put(out, "interpretation", vstr(b));
  put(out, "expomean", vnum(o[2]));
  put(out, "citation", vstr("Wagstaff, Paci & van Doorslaer (1991) Soc Sci Med 33(5):545-557"));
  put(out, "extra", e);
  return true;
}
