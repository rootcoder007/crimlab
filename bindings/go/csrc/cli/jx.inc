// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Helpers of the entries in the j*.inc files (included by jrun.cpp, in
// its unnamed namespace): the readers of the R binding (src/rchk.cpp,
// src/spin.cpp and the b*.cpp files) over JSON values, with the same
// checks in the same order and the same words, and the builders of
// results in the form R's values take in JSON: a data frame an object of
// columns, a named vector an object of numbers, a list an object, an
// unnamed vector of one value outside a data frame that value, a matrix
// an array of its rows.

constexpr double Inf = HUGE_VAL;

// The number of values: 0 for null, 1 for one value.
inline std::size_t vlen(const Val &o) {
  if (o.t == Val::NUL) return 0;
  return o.t == Val::ARR || o.t == Val::OBJ ? o.a.size() : 1;
}

// is.numeric(): a number, or an array or object of numbers and nulls.
inline bool isnum(const Val &o) {
  if (o.t == Val::NUM) return true;
  if (o.t != Val::ARR && o.t != Val::OBJ) return false;
  for (const Val &e : o.a) {
    if (e.t != Val::NUM && e.t != Val::NUL) return false;
  }
  return true;
}

// The values of a numeric argument, as as.double() gives them (null is NA).
inline std::vector<double> dbl(const Val &o) {
  Dv d;
  readd(o, d);
  for (std::size_t i = 0; i < d.n(); ++i) {
    if (std::isnan(d.v[i])) d.v[i] = NA;
  }
  return d.v;
}

// as.character(seq_len(n)) when nm is empty, else nm.
inline std::vector<std::string> nmsof(const std::vector<std::string> &nm, std::size_t n) { return grp(nm, n); }

// Results.

inline Val vobj() {
  Val v;
  v.t = Val::OBJ;
  return v;
}
// Integers; INT_MIN (R's NA_integer_) is null.
inline Val vints(const int *x, std::size_t n) {
  Val v = varr();
  for (std::size_t i = 0; i < n; ++i) v.a.push_back(x[i] == INT_MIN ? Val() : vnum(x[i]));
  return v;
}
// TRUE, FALSE; a negative value (NA) is null.
inline Val vlgls(const int *x, std::size_t n) {
  Val v = varr();
  for (std::size_t i = 0; i < n; ++i) v.a.push_back(x[i] < 0 ? Val() : vbool(x[i] != 0));
  return v;
}
// An unnamed vector outside a data frame: one value is that value.
inline Val unbox(Val v) {
  if (v.t == Val::ARR && v.a.size() == 1) return v.a[0];
  return v;
}
inline Val vauto(const std::vector<double> &x) { return unbox(vnums(x)); }
// A named vector of numbers.
inline Val nvec(const std::vector<std::string> &nm, const double *x) {
  Val o = vobj();
  for (std::size_t i = 0; i < nm.size(); ++i) put(o, nm[i].c_str(), vnum(x[i]));
  return o;
}
// A column-major nr x nc matrix as an array of its rows.
inline Val vmat(const double *x, int nr, int nc) {
  Val v = varr();
  for (int i = 0; i < nr; ++i) {
    Val r = varr();
    for (int j = 0; j < nc; ++j) r.a.push_back(vnum(x[i + static_cast<std::size_t>(nr) * j]));
    v.a.push_back(r);
  }
  return v;
}
inline Val vmat(const std::vector<double> &x, int nr, int nc) { return vmat(x.data(), nr, nc); }
inline Val vlmat(const int *x, int nr, int nc) {
  Val v = varr();
  for (int i = 0; i < nr; ++i) {
    Val r = varr();
    for (int j = 0; j < nc; ++j) {
      const int e = x[i + static_cast<std::size_t>(nr) * j];
      r.a.push_back(e < 0 ? Val() : vbool(e != 0));
    }
    v.a.push_back(r);
  }
  return v;
}
// A column-major d0 x d1 x d2 array as nested arrays, first index outermost.
inline Val varr3(const double *x, const int *d) {
  Val v = varr();
  for (int i = 0; i < d[0]; ++i) {
    Val r = varr();
    for (int j = 0; j < d[1]; ++j) {
      Val s = varr();
      for (int k = 0; k < d[2]; ++k) {
        s.a.push_back(vnum(x[i + static_cast<std::size_t>(d[0]) * (j + static_cast<std::size_t>(d[1]) * k)]));
      }
      r.a.push_back(s);
    }
    v.a.push_back(r);
  }
  return v;
}
inline Val vstrs(const std::vector<std::string> &x, const std::vector<int> &code) {
  Val v = varr();
  for (int k : code) v.a.push_back(k < 0 ? Val() : vstr(x[static_cast<std::size_t>(k)]));
  return v;
}

// Arguments.

// .flag(): TRUE or FALSE, one value.
inline bool flag(Ctx &c, const Val &o, const char *nm, bool &b) {
  const Val *e = o.t == Val::ARR && o.a.size() == 1 ? &o.a[0] : &o;
  if (e->t != Val::BOOL) return c.fail("`%s` must be TRUE or FALSE.", nm);
  b = e->b;
  return true;
}
// A flag with its default when left out.
inline bool flagd(Ctx &c, const char *nm, bool def, bool &b) {
  if (!c.has(nm)) {
    b = def;
    return true;
  }
  return flag(c, c[nm], nm, b);
}

// seed() (src/spin.cpp): one whole number from 0. In R a seed left out is
// drawn from R's generator (set.seed()); outside R there is none, so a
// function that draws needs it.
inline bool seed(Ctx &c, unsigned int &s) {
  if (!c.has("seed")) return c.fail("`seed` is needed outside R: one whole number from 0.");
  int v = 0;
  CK(whole(c, c["seed"], "seed", 0, v));
  s = static_cast<unsigned int>(v);
  return true;
}

// The argument nm, or the number def when left out.
inline const Val &arg(Ctx &c, const char *nm, Val &def, double x) {
  if (c.has(nm)) return c[nm];
  def = vnum(x);
  return def;
}
// whole(), one(), p01() and num() of an argument with an R default.
inline bool wholed(Ctx &c, const char *nm, double def, int lo, int &r) {
  Val d;
  return whole(c, arg(c, nm, d, def), nm, lo, r);
}
inline bool oned(Ctx &c, const char *nm, double def, double lo, double hi, double &v) {
  Val d;
  return one(c, arg(c, nm, d, def), nm, lo, hi, v);
}
inline bool p01d(Ctx &c, const char *nm, double def, double &v) {
  Val d;
  return p01(c, arg(c, nm, d, def), nm, v);
}
inline double val1d(Ctx &c, const char *nm, double def) { return c.has(nm) ? val1(c[nm]) : def; }

// .pick() with its default when left out.
inline bool pickd(Ctx &c, const char *nm, const char *const *opt, int n, int def, int &k) {
  if (!c.has(nm)) {
    k = def;
    return true;
  }
  return pick(c, c[nm], nm, opt, n, k);
}

// The labels of an argument (as.character()), false when one is missing
// or not a label.
inline bool alabs(const Val &o, std::vector<std::string> &s) {
  s.clear();
  std::string t;
  if (o.t == Val::ARR || o.t == Val::OBJ) {
    for (const Val &e : o.a) {
      if (!lab(e, t)) return false;
      s.push_back(t);
    }
    return true;
  }
  if (!lab(o, t)) return false;
  s.push_back(t);
  return true;
}

// codes() (src/spin.cpp): labels to codes from 0, by the sorted distinct
// labels (as.character(), compared byte by byte). A factor is read in R by
// its levels; JSON has none, so labels are always sorted.
struct Lab {
  std::vector<int> code;
  std::vector<std::string> lev;
  int g = 0;
};
inline bool codes(Ctx &c, const Val &o, const char *nm, std::size_t n, Lab &L) {
  if (vlen(o) != n) return c.fail("`%s` must have one label per value, %d.", nm, static_cast<int>(n));
  std::vector<std::string> s;
  if (!alabs(o, s)) {
    return c.fail("`%s` must be labels (numbers, text or a factor) with no missing values.", nm);
  }
  for (const Val &e : o.a) {
    if (e.t == Val::ARR || e.t == Val::OBJ) return notvec(c, nm);
  }
  std::vector<int> ord(n);
  for (std::size_t i = 0; i < n; ++i) ord[i] = static_cast<int>(i);
  std::stable_sort(ord.begin(), ord.end(), [&s](int a, int b) {
    return std::strcmp(s[static_cast<std::size_t>(a)].c_str(), s[static_cast<std::size_t>(b)].c_str()) < 0;
  });
  L.code.assign(n, 0);
  L.lev.clear();
  int g = -1;
  for (std::size_t q = 0; q < n; ++q) {
    const std::size_t i = static_cast<std::size_t>(ord[q]);
    if (q == 0 || s[static_cast<std::size_t>(ord[q - 1])] != s[i]) {
      ++g;
      L.lev.push_back(s[i]);
    }
    L.code[i] = g;
  }
  L.g = g + 1;
  return true;
}

// bin01() (src/spin.cpp): 0 or 1 outcomes, TRUE or FALSE or numbers, n.
inline bool bin01(Ctx &c, const Val &o, const char *nm, std::size_t n, std::vector<int> &b) {
  b.clear();
  const char *m = "`%s` must be 0 or 1 (or TRUE or FALSE), one per case.";
  if (vlen(o) != n) return c.fail(m, nm);
  auto one1 = [&b](const Val &e) {
    if (e.t == Val::BOOL) b.push_back(e.b ? 1 : 0);
    else if (e.t == Val::NUM && (e.n == 0 || e.n == 1)) b.push_back(static_cast<int>(e.n));
    else return false;
    return true;
  };
  if (o.t == Val::ARR || o.t == Val::OBJ) {
    for (const Val &e : o.a) {
      if (!one1(e)) return c.fail(m, nm);
    }
  } else if (!one1(o)) {
    return c.fail(m, nm);
  }
  return true;
}

// ref() (src/rchk.cpp): the label named, or the position from 1, among the
// n labels g; an and what name the argument and what it picks.
inline bool rref(Ctx &c, const Val &r, const std::vector<std::string> &g, const char *an, const char *what, int &i) {
  const Val *e = r.t == Val::ARR && r.a.size() == 1 ? &r.a[0] : &r;
  if (e->t == Val::STR) {
    for (std::size_t k = 0; k < g.size(); ++k) {
      if (g[k] == e->s) {
        i = static_cast<int>(k) + 1;
        return true;
      }
    }
    return c.fail("`%s` names no %s: %s.", an, what, e->s.c_str());
  }
  const double v = val1(r);
  const int n = static_cast<int>(g.size());
  if (std::isnan(v) || v != std::floor(v) || v < 1 || v > n) {
    return c.fail("`%s` must be one %s name or a position from 1 to %d.", an, what, n);
  }
  i = static_cast<int>(v);
  return true;
}

// A numeric array of nd dimensions, given as nested arrays with the first
// index outermost (a matrix as an array of its rows): its dimensions into
// dim and its values in R's order (column-major) into v; finite, or NA
// (null) when na. The words of narr() in the b*.cpp files.
inline bool shape(const Val &o, int nd, int *dim, int at) {
  if (at == nd) return o.t == Val::NUM || o.t == Val::NUL;
  if (o.t != Val::ARR || o.a.size() > static_cast<std::size_t>(INT_MAX)) return false;
  const int n = static_cast<int>(o.a.size());
  if (dim[at] < 0) dim[at] = n;
  if (dim[at] != n) return false;
  if (n == 0) {
    for (int j = at + 1; j < nd; ++j) dim[j] = dim[j] < 0 ? 0 : dim[j];
  }
  for (const Val &e : o.a) {
    if (!shape(e, nd, dim, at + 1)) return false;
  }
  return true;
}
inline void fill(const Val &o, int nd, const int *dim, int at, std::size_t off, std::size_t step, std::vector<double> &v,
          bool &bad, bool na) {
  if (at == nd) {
    const double x = o.t == Val::NUM ? o.n : NA;
    if (na ? std::isinf(x) : !std::isfinite(x)) bad = true;
    v[off] = x;
    return;
  }
  for (std::size_t i = 0; i < o.a.size(); ++i) {
    fill(o.a[i], nd, dim, at + 1, off + i * step, step * static_cast<std::size_t>(dim[at]), v, bad, na);
  }
}
inline bool narr(Ctx &c, const Val &o, const char *nm, int nd, int *dim, std::vector<double> &v, bool na = false) {
  for (int j = 0; j < nd; ++j) dim[j] = -1;
  if (!shape(o, nd, dim, 0)) {
    return c.fail(nd == 2 ? "`%s` must be a numeric matrix." : "`%s` must be a numeric array of 3 dimensions.", nm);
  }
  double tot = 1;
  for (int j = 0; j < nd; ++j) tot *= dim[j];
  if (tot > 4503599627370496.0) return c.fail("`%s` is too large.", nm);
  v.assign(static_cast<std::size_t>(tot), 0.0);
  bool bad = false;
  fill(o, nd, dim, 0, 0, 1, v, bad, na);
  if (bad) return c.fail(na ? "`%s` must hold finite numbers or NA." : "`%s` must hold finite numbers.", nm);
  return true;
}
inline bool nmat(Ctx &c, const Val &o, const char *nm, int &nr, int &nc, std::vector<double> &v, bool na = false) {
  int d[2];
  CK(narr(c, o, nm, 2, d, v, na));
  nr = d[0];
  nc = d[1];
  return true;
}

// A data frame (an object of equal-length columns) as the matrix
// as.matrix() makes of it; any other value as it is.
inline Val asmat(const Val &o) {
  if (o.t != Val::OBJ || o.a.empty()) return o;
  const std::size_t n = vlen(o.a[0]);
  for (const Val &e : o.a) {
    if (e.t != Val::ARR || e.a.size() != n) return o;
  }
  Val m = varr();
  for (std::size_t i = 0; i < n; ++i) {
    Val r = varr();
    for (const Val &e : o.a) r.a.push_back(e.a[i]);
    m.a.push_back(r);
  }
  return m;
}

// byname() (src/spin.cpp) on whole values: when two or more of the
// arguments given are objects of single values, all of one length with
// names that do not repeat, the later ones are put in the order of the
// first, which they must name alike.
inline bool vbyname(Ctx &c, std::vector<std::pair<const char *, Val *>> as) {
  std::vector<std::size_t> g;
  for (std::size_t j = 0; j < as.size(); ++j) {
    const Val &x = *as[j].second;
    if (x.t == Val::NUL) continue;
    if (x.t != Val::OBJ) return true;
    for (const Val &e : x.a) {
      if (e.t == Val::ARR || e.t == Val::OBJ) return true;
    }
    std::vector<std::string> s = x.k;
    std::sort(s.begin(), s.end());
    if (std::adjacent_find(s.begin(), s.end()) != s.end()) return true;
    if (!g.empty() && x.a.size() != as[g[0]].second->a.size()) return true;
    g.push_back(j);
  }
  if (g.size() < 2) return true;
  const Val &f = *as[g[0]].second;
  for (std::size_t q = 1; q < g.size(); ++q) {
    Val &d = *as[g[q]].second;
    Val o = vobj();
    for (std::size_t i = 0; i < f.k.size(); ++i) {
      const auto it = std::find(d.k.begin(), d.k.end(), f.k[i]);
      if (it == d.k.end()) {
        return c.fail("`%s` must name the same units as `%s`; `%s` has no %s.", as[g[q]].first, as[g[0]].first,
                      as[g[q]].first, f.k[i].c_str());
      }
      put(o, f.k[i].c_str(), d.a[static_cast<std::size_t>(it - d.k.begin())]);
    }
    d = o;
  }
  return true;
}

// The arguments named, copied, for vbyname().
struct Args {
  std::vector<Val> v;
  std::vector<const char *> nm;
  Args(Ctx &c, std::initializer_list<const char *> n) : nm(n) {
    for (const char *k : nm) v.push_back(c[k]);
  }
  bool byname(Ctx &c) {
    std::vector<std::pair<const char *, Val *>> as;
    for (std::size_t j = 0; j < v.size(); ++j) as.emplace_back(nm[j], &v[j]);
    return vbyname(c, as);
  }
  const Val &operator[](std::size_t j) const { return v[j]; }
};

// A column of a data frame by name; null when absent.
inline const Val &col(const Val &df, const char *nm) {
  if (df.t != Val::OBJ) return nul;
  for (std::size_t j = 0; j < df.k.size(); ++j) {
    if (df.k[j] == nm) return df.a[j];
  }
  return nul;
}

// A data frame of columns built in order.
struct Frame {
  Val o = vobj();
  void num(const char *nm, const double *x, std::size_t n) { put(o, nm, vnums(x, n)); }
  void num(const char *nm, const std::vector<double> &x) { put(o, nm, vnums(x)); }
  void ints(const char *nm, const int *x, std::size_t n) { put(o, nm, vints(x, n)); }
  void lgl(const char *nm, const int *x, std::size_t n) { put(o, nm, vlgls(x, n)); }
  void str(const char *nm, const std::vector<std::string> &x) { put(o, nm, vstrs(x)); }
  void add(const char *nm, Val v) { put(o, nm, std::move(v)); }
};

// The normal quantile at p (lower tail).
inline double qnorm(double p) { return clab_qnorm(p); }
// The normal quantile of a two-sided interval at level conf.
inline bool zconf(Ctx &c, const char *nm, double def, double &z) {
  double cf;
  CK(oned(c, nm, def, 0.0, 1.0, cf));
  z = qnorm(0.5 + 0.5 * cf);
  return true;
}

// Points (src/spin.cpp, xy()): x and y as two arrays, or x alone as an
// object of columns x and y (X and Y, or lon and lat) or an array of rows
// of two; with lonlat, longitude and latitude in degrees, projected by
// clab_eqproj about their own centre so that the planar kernels receive
// metres. RAW: the input is longitude and latitude, not projected; there
// an array of two numbers is one point. PLANAR: lon and lat are refused.
struct Pts {
  std::vector<double> x, y;
  int n = 0;
  bool ll = false;
  double lon0 = 0, lat0 = 0;
};
enum Mode { PROJ, PLANAR, RAW };
constexpr double earth = 6371008.8;
constexpr double deg = 3.14159265358979323846 / 180.0;

inline int find(const std::vector<std::string> &nm, const char *a) {
  for (std::size_t j = 0; j < nm.size(); ++j) {
    if (nm[j] == a) return static_cast<int>(j);
  }
  return -1;
}

inline bool xy(Ctx &c, const Val &x, const Val &y, Mode m, Pts &p) {
  p = Pts();
  if (m == RAW) p.ll = true;
  if (m == PROJ) CK(flagd(c, "lonlat", false, p.ll));
  const char *nx = p.ll ? "lon" : "x", *ny = p.ll ? "lat" : "y";
  Val cx = x, cy = y;
  int dim[2] = {-1, -1};
  if (y.t == Val::NUL && m == RAW && x.t == Val::ARR && isnum(x) && x.a.size() == 2) {
    cx = x.a[0];
    cy = x.a[1];
  } else if (y.t == Val::NUL) {
    std::vector<std::string> nm;
    std::size_t nc = 0;
    const bool lst = x.t == Val::OBJ;
    if (lst) {
      nm = x.k;
      nc = x.a.size();
    } else if (x.t == Val::ARR && !x.a.empty() && shape(x, 2, dim, 0)) {
      nc = static_cast<std::size_t>(dim[1]);
    } else {
      return c.fail("`y` is missing: give `x` and `y`, or a matrix or data frame of coordinates as `x`.");
    }
    int jx = -1, jy = -1;
    bool isll = false;
    if (!nm.empty()) {
      if (m != RAW) {
        jx = find(nm, "x");
        jy = find(nm, "y");
        if (jx < 0 || jy < 0) {
          jx = find(nm, "X");
          jy = find(nm, "Y");
        }
      }
      if (jx < 0 || jy < 0) {
        jx = find(nm, "lon");
        jy = find(nm, "lat");
        isll = true;
      }
    } else if (nc == 2) {
      jx = 0;
      jy = 1;
    }
    if (jx < 0 || jy < 0) {
      return c.fail(m == RAW ? "`x` must have columns named lon and lat, or be a matrix of exactly two unnamed columns."
                             : "`x` must have columns named x and y, or lon and lat, or be a matrix of exactly two "
                               "unnamed columns.");
    }
    if (isll && m == PLANAR) {
      return c.fail("columns lon and lat hold degrees: this function takes planar coordinates; project them first "
                    "with eqproj().");
    }
    if (isll && !p.ll) return c.fail("columns lon and lat hold degrees: set `lonlat = TRUE`.");
    if (lst) {
      cx = x.a[static_cast<std::size_t>(jx)];
      cy = x.a[static_cast<std::size_t>(jy)];
    } else {
      cx = varr();
      cy = varr();
      for (const Val &r : x.a) {
        cx.a.push_back(r.a[static_cast<std::size_t>(jx)]);
        cy.a.push_back(r.a[static_cast<std::size_t>(jy)]);
      }
    }
  }
  CK(vbyname(c, {{"x", &cx}, {"y", &cy}}));
  Dv vx, vy;
  CK(num(c, cx, nx, vx));
  CK(num(c, cy, ny, vy, static_cast<long>(vx.n())));
  p.n = static_cast<int>(vx.n());
  if (!p.ll) {
    p.x = vx.v;
    p.y = vy.v;
    return true;
  }
  double ss = 0.0, sc = 0.0, sl = 0.0;
  for (int i = 0; i < p.n; ++i) {
    const double lo = vx.v[static_cast<std::size_t>(i)], la = vy.v[static_cast<std::size_t>(i)];
    if (lo < -180.0 || lo > 180.0) return c.fail("`%s` must be longitudes from -180 to 180.", nx);
    if (la < -90.0 || la > 90.0) return c.fail("`%s` must be latitudes from -90 to 90.", ny);
    ss += std::sin(lo * deg);
    sc += std::cos(lo * deg);
    sl += la;
  }
  p.lon0 = std::atan2(ss, sc) / deg;
  p.lat0 = sl / p.n;
  if (m == RAW) {
    p.x = vx.v;
    p.y = vy.v;
    return true;
  }
  p.x.assign(vx.n(), 0.0);
  p.y.assign(vx.n(), 0.0);
  if (clab_eqproj(vx.p(), vy.p(), p.n, p.lon0, p.lat0, earth, p.x.data(), p.y.data()) != CLAB_OK) {
    return c.fail("the points reach the antipode of their centre: no equal-area projection about it holds them all.");
  }
  return true;
}
// The points of the arguments x and y of the call.
inline bool xy(Ctx &c, Mode m, Pts &p) { return xy(c, c["x"], c["y"], m, p); }

// fwd(): points given with the points, to their projection.
inline bool fwd(Ctx &c, const Pts &p, double *px, double *py, int n) {
  if (!p.ll) return true;
  for (int i = 0; i < n; ++i) {
    if (px[i] < -180.0 || px[i] > 180.0 || py[i] < -90.0 || py[i] > 90.0) {
      return c.fail("a longitude or latitude given with the points is out of range.");
    }
  }
  if (clab_eqproj(px, py, n, p.lon0, p.lat0, earth, px, py) != CLAB_OK) {
    return c.fail("a point given with the points lies at the antipode of their centre.");
  }
  return true;
}
// back(): results back to longitude and latitude.
inline void back(const Pts &p, double *px, double *py, int n) {
  if (p.ll) clab_eqinv(px, py, n, p.lon0, p.lat0, earth, px, py);
}

// win(): c(xmin, xmax, ymin, ymax).
inline bool win(Ctx &c, const Val &s, double *w) {
  Dv v;
  CK(num(c, s, "win", v, 4));
  for (int j = 0; j < 4; ++j) w[j] = v[static_cast<std::size_t>(j)];
  if (!(w[0] < w[1] && w[2] < w[3])) {
    return c.fail("`win` must be c(xmin, xmax, ymin, ymax) with xmin < xmax and ymin < ymax.");
  }
  return true;
}
inline bool inside(const double *w, double x, double y) { return x >= w[0] && x <= w[1] && y >= w[2] && y <= w[3]; }

// Weights: none, or n finite numbers from 0 with a positive sum.
inline bool pwts(Ctx &c, const Val &w, std::size_t n, std::vector<double> &out) {
  out.clear();
  if (w.t == Val::NUL) return true;
  Dv v;
  CK(num(c, w, "w", v, static_cast<long>(n), 0.0));
  bool any = false;
  for (double e : v.v) any = any || e > 0;
  if (!any) return c.fail("`w` must have a positive sum.");
  out = v.v;
  return true;
}
inline const double *ptr(const std::vector<double> &v) { return v.empty() ? nullptr : v.data(); }

// The helpers of src/bcol.h and src/bfx.h.

// bj::inr(): one number from lo to hi, both included (hi Inf: from lo up).
inline bool binr(Ctx &c, const Val &x, const char *nm, double lo, double hi, double &v) {
  v = val1(x);
  if (!(v >= lo && v <= hi) || !std::isfinite(v)) {
    if (std::isfinite(hi)) return c.fail("`%s` must be one number from %s to %s.", nm, rfmt(lo).c_str(), rfmt(hi).c_str());
    return c.fail("`%s` must be one number from %s.", nm, rfmt(lo).c_str());
  }
  return true;
}
inline bool binrd(Ctx &c, const char *nm, double def, double lo, double hi, double &v) {
  Val d;
  return binr(c, arg(c, nm, d, def), nm, lo, hi, v);
}
// A kernel's "no value" (NaN) as R's NA.
inline double na(double v) { return std::isnan(v) ? NA : v; }
// bfx::nna(): numbers or NA (true and false as 1 and 0), no infinite
// value, n of them or one recycled to n (n < 0: any length, at least one).
inline bool nna(Ctx &c, const Val &x, const char *nm, long n, std::vector<double> &v) {
  const char *m = "`%s` must be numbers or NA, one per case or one for all.";
  const std::size_t k = vlen(x);
  bool ok = k >= 1 && (n < 0 || k == static_cast<std::size_t>(n) || k == 1);
  std::vector<double> d;
  const bool arr = x.t == Val::ARR || x.t == Val::OBJ;
  for (std::size_t i = 0; ok && i < k; ++i) {
    const Val &u = arr ? x.a[i] : x;
    if (u.t == Val::NUM) d.push_back(u.n);
    else if (u.t == Val::BOOL) d.push_back(u.b ? 1.0 : 0.0);
    else if (u.t == Val::NUL) d.push_back(NA);
    else ok = false;
  }
  if (!ok) return c.fail(m, nm);
  const std::size_t out = n < 0 ? k : static_cast<std::size_t>(n);
  v.assign(out, 0.0);
  for (std::size_t i = 0; i < out; ++i) {
    v[i] = d[k == 1 ? 0 : i];
    if (std::isinf(v[i])) return c.fail(m, nm);
  }
  return true;
}
// bfx::mat(): a numeric matrix of finite numbers from lo with nc columns
// (any when nc < 0) and at least one row; a vector is one column. Values
// column-major into v, rows into nr.
inline bool bmat(Ctx &c, const Val &x, const char *nm, int nc, int &nr, std::vector<double> &v, double lo = -HUGE_VAL) {
  int dim[2] = {-1, -1};
  const bool m = x.t == Val::ARR && !x.a.empty() && x.a[0].t == Val::ARR && shape(x, 2, dim, 0);
  const bool vec = isnum(x);
  const std::size_t r = m ? static_cast<std::size_t>(dim[0]) : vlen(x), k = m ? static_cast<std::size_t>(dim[1]) : 1;
  if (!(m || vec) || r < 1 || (nc >= 0 && k != static_cast<std::size_t>(nc))) {
    if (nc >= 0) return c.fail("`%s` must be a numeric matrix of %d columns.", nm, nc);
    return c.fail("`%s` must be a numeric matrix.", nm);
  }
  if (m) {
    v.assign(r * k, 0.0);
    bool bad = false;
    fill(x, 2, dim, 0, 0, 1, v, bad, true);
  } else {
    v = dbl(x);
  }
  for (double e : v) {
    if (!std::isfinite(e) || e < lo) return c.fail("`%s` must hold finite numbers from %g.", nm, lo);
  }
  if (r > static_cast<std::size_t>(INT_MAX)) return c.fail("`%s` is too long for crimlab.", nm);
  nr = static_cast<int>(r);
  return true;
}
// Codes from 0 (-1 for NA) as the labels of a factor over lev.
inline Val fac(const int *v, std::size_t n, const char *const *lev) {
  Val o = varr();
  for (std::size_t i = 0; i < n; ++i) o.a.push_back(v[i] < 0 ? Val() : vstr(lev[v[i]]));
  return o;
}
// bj::perlev(): numbers from 0 given per label of lev, named ones matched
// by name (names of other labels ignored), unnamed ones in the order of lev.
inline bool perlev(Ctx &c, const Val &x, const char *nm, const std::vector<std::string> &lev, std::vector<double> &out) {
  Dv d;
  CK(num(c, x, nm, d, -1, 0.0));
  const std::size_t k = lev.size();
  out.assign(k, 0.0);
  if (d.nm.empty()) {
    if (d.n() != k) return c.fail("`%s` must have one value per label, %d, or be named by label.", nm, static_cast<int>(k));
    out = d.v;
    return true;
  }
  std::vector<int> seen(k, 0);
  for (std::size_t i = 0; i < d.n(); ++i) {
    const auto f = std::find(lev.begin(), lev.end(), d.nm[i]);
    if (f == lev.end()) continue;
    const std::size_t j = static_cast<std::size_t>(f - lev.begin());
    if (seen[j]++) return c.fail("`%s` names %s twice.", nm, d.nm[i].c_str());
    out[j] = d.v[i];
  }
  for (std::size_t j = 0; j < k; ++j) {
    if (!seen[j]) return c.fail("`%s` has no value for %s.", nm, lev[j].c_str());
  }
  return true;
}

// gvec() (src/spin.cpp): a one-row matrix or a one-row data frame (an
// object of single values or of one-value columns) as the vector it holds.
inline Val gvec(const Val &x) {
  if (x.t == Val::OBJ) {
    Val v = vobj();
    for (std::size_t j = 0; j < x.a.size(); ++j) {
      const Val &e = x.a[j];
      put(v, x.k[j].c_str(), e.t == Val::ARR && e.a.size() == 1 ? e.a[0] : e);
    }
    return v;
  }
  int dim[2] = {-1, -1};
  if (x.t == Val::ARR && shape(x, 2, dim, 0) && (dim[0] == 1 || dim[1] == 1)) {
    Val v = varr();
    for (const Val &r : x.a) {
      for (const Val &e : r.a) v.a.push_back(e);
    }
    return v;
  }
  return x;
}
