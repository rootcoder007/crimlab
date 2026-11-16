// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Exposure and crime (src/apx.cpp, src/cnt.cpp) for clab_jrun(): the R
// wrappers of R/apx.R and R/apx2.R in C++, with their checks in their
// order and words, their defaults, and their results written as
// dev/clref.R writes R's. A data frame is an object of columns: a column
// of text is a character column, one of numbers a numeric column, one of
// true and false a logical column. JSON has no dates, so a day column is
// numbers (smkiv() then needs `seas` for its seasons), and no factors, so
// labels are sorted byte by byte. The cluster bootstrap (se = "boot")
// draws what R draws after set.seed(seed): R's Mersenne-Twister and
// sample() (igh::Rng, src/igk.h); outside R there is no generator to draw
// from without a seed, so it needs one. The functions that take an
// estimator function (plcday, runfx, hetfx, ccross, jntx, dosebin,
// typefx) are in jax2.inc.

namespace ax {

using Cols = std::vector<std::string>;
using Ix = std::vector<std::size_t>;

inline std::size_t nr(const Val &d) { return d.a.empty() ? 0 : d.a[0].a.size(); }

// A value as as.numeric() reads it: null is NA.
inline double nv(const Val &e) {
  double x = 0.0;
  if (e.t == Val::NUM) return e.n;
  if (e.t == Val::BOOL) return e.b ? 1.0 : 0.0;
  if (snum(e, x)) return x;
  return NA;
}
// is.na(): null or NaN.
inline bool vna(const Val &e) {
  double x = 0.0;
  return e.t == Val::NUL || (e.t == Val::NUM && std::isnan(e.n)) || (snum(e, x) && std::isnan(x));
}
// as.character() of a value.
inline std::string vlab(const Val &e) {
  std::string s;
  return lab(e, s) ? s : std::string("NA");
}
// class() of a column: character, numeric, or logical (true, false, null).
inline const char *ccls(const Val &a) {
  bool num = false;
  double x = 0.0;
  for (const Val &e : a.a) {
    if (e.t == Val::STR && !snum(e, x)) return "character";
    if (e.t == Val::NUM || e.t == Val::STR) num = true;
  }
  return num ? "numeric" : "logical";
}
inline bool isnumc(const Val &a) { return std::strcmp(ccls(a), "numeric") == 0; }
// is.numeric() of an argument: numbers and nulls, not all null.
inline bool numv(const Val &o) {
  if (!isnum(o)) return false;
  if (o.t != Val::ARR && o.t != Val::OBJ) return true;
  for (const Val &e : o.a) {
    if (e.t != Val::NUL) return true;
  }
  return false;
}
inline std::vector<double> nums(const Val &o) {
  std::vector<double> v;
  if (o.t == Val::ARR || o.t == Val::OBJ) {
    for (const Val &e : o.a) v.push_back(nv(e));
  } else if (o.t != Val::NUL) {
    v.push_back(nv(o));
  }
  return v;
}
inline std::vector<double> colv(const Val &a) {
  std::vector<double> v;
  v.reserve(a.a.size());
  for (const Val &e : a.a) v.push_back(nv(e));
  return v;
}
inline Cols coll(const Val &a) {
  Cols v;
  v.reserve(a.a.size());
  for (const Val &e : a.a) v.push_back(vlab(e));
  return v;
}
// The values of an argument as one array.
inline Val asarr(const Val &o) {
  if (o.t == Val::ARR) return o;
  Val a = varr();
  if (o.t == Val::OBJ) a.a = o.a;
  else if (o.t != Val::NUL) a.a.push_back(o);
  return a;
}
inline std::string join(const Cols &v, const char *sep) {
  std::string s;
  for (std::size_t i = 0; i < v.size(); ++i) s += (i ? sep : "") + v[i];
  return s;
}
inline bool has(const Cols &v, const std::string &s) { return std::find(v.begin(), v.end(), s) != v.end(); }
// TRUE or FALSE, one value.
inline bool tf(const Val &o, bool &b) {
  const Val *e = o.t == Val::ARR && o.a.size() == 1 ? &o.a[0] : &o;
  if (e->t != Val::BOOL) return false;
  b = e->b;
  return true;
}
// R's %% of doubles (myfmod of src/main/arithmetic.c).
inline double rmod(double x, double y) {
  const double q = x / y;
  const long double t = static_cast<long double>(x) - std::floor(q) * static_cast<long double>(y);
  return static_cast<double>(t - std::floor(t / y) * y);
}
// mean() and sum() of doubles as R computes them: in long double, the
// mean with its second pass.
inline double lmean(const std::vector<double> &x) {
  long double s = 0.0;
  for (double v : x) s += v;
  s /= static_cast<long double>(x.size());
  if (std::isfinite(static_cast<double>(s))) {
    long double t = 0.0;
    for (double v : x) t += v - s;
    s += t / static_cast<long double>(x.size());
  }
  return static_cast<double>(s);
}
inline double lsum(const std::vector<double> &x) {
  long double s = 0.0;
  for (double v : x) s += v;
  return static_cast<double>(s);
}
// cov() of the columns of x (nr x p, column-major) over the rows ok, as
// R's stats computes it (cov.c): long double means with the second pass,
// long double sums of the products, divided by the rows less one.
inline std::vector<double> rcov(const std::vector<double> &x, std::size_t n, int p, const std::vector<char> &ok) {
  const std::size_t sp = zz(p);
  long double nobs = 0;
  for (char k : ok) nobs += k ? 1 : 0;
  std::vector<double> xm(sp), v(sp * sp);
  for (std::size_t i = 0; i < sp; ++i) {
    const double *xx = &x[i * n];
    long double sum = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
      if (ok[k]) sum += xx[k];
    }
    long double tmp = sum / nobs;
    if (std::isfinite(static_cast<double>(tmp))) {
      sum = 0.0;
      for (std::size_t k = 0; k < n; ++k) {
        if (ok[k]) sum += xx[k] - tmp;
      }
      tmp = tmp + sum / nobs;
    }
    xm[i] = static_cast<double>(tmp);
  }
  for (std::size_t i = 0; i < sp; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      const double *xx = &x[i * n], *yy = &x[j * n];
      long double sum = 0.0;
      for (std::size_t k = 0; k < n; ++k) {
        if (ok[k]) sum += static_cast<long double>(xx[k] - xm[i]) * (yy[k] - xm[j]);
      }
      v[i + j * sp] = v[j + i * sp] = static_cast<double>(sum / (nobs - 1));
    }
  }
  return v;
}
// The warning R gives, prefixed with the function.
inline void wrn(Ctx &c, const std::string &m) { c.warn.push_back(std::string(c.fn) + ": " + m); }

// .ckmarg(): left out, or all the choices, is the first; else one text
// matched exactly or by a prefix of one choice alone (ckarg of jsar.inc).
inline bool marg(Ctx &c, const char *nm, const std::vector<const char *> &ch, int &k) {
  return ckarg(c, nm, std::vector<std::string>(ch.begin(), ch.end()), k);
}

// Data frames: an object of columns (arrays of one length).

inline const char *vcls(const Val &o) {
  double x = 0.0;
  switch (o.t) {
    case Val::NUL: return "NULL";
    case Val::NUM: return "numeric";
    case Val::BOOL: return "logical";
    case Val::STR: return snum(o, x) ? "numeric" : "character";
    case Val::OBJ: return "list";
    default: break;
  }
  for (const Val &e : o.a) {
    if (e.t == Val::ARR || e.t == Val::OBJ) return "list";
  }
  bool lg = false;
  const char *cl = jcls(o, lg);
  return cl != nullptr ? cl : "numeric";
}
// .axdat(): a data frame with rows, each column an array.
inline bool dat(Ctx &c, const Val &o, const char *what, Val &d) {
  bool df = o.t == Val::OBJ;
  std::size_t n = 0;
  for (std::size_t j = 0; df && j < o.a.size(); ++j) {
    const Val &v = o.a[j];
    if (v.t == Val::OBJ || v.t == Val::NUL) df = false;
    if (v.t == Val::ARR) {
      for (const Val &e : v.a) df = df && e.t != Val::ARR && e.t != Val::OBJ;
    }
    const std::size_t len = v.t == Val::ARR ? v.a.size() : 1;
    if (j > 0 && len != n) df = false;
    n = len;
  }
  if (!df) return ckstop(c, "`%s` must be a data frame, got %s", what, vcls(o));
  if (n == 0) return ckstop(c, "`%s` has no rows", what);
  d = vobj();
  for (std::size_t j = 0; j < o.a.size(); ++j) put(d, o.k[j].c_str(), asarr(o.a[j]));
  return true;
}
inline const Val *fc(const Val &d, const std::string &nm) {
  for (std::size_t j = 0; j < d.k.size(); ++j) {
    if (d.k[j] == nm) return &d.a[j];
  }
  return nullptr;
}
inline void setc(Val &d, const std::string &nm, Val v) {
  for (std::size_t j = 0; j < d.k.size(); ++j) {
    if (d.k[j] == nm) {
      d.a[j] = std::move(v);
      return;
    }
  }
  put(d, nm.c_str(), std::move(v));
}
// d[r, , drop = FALSE]
inline Val rows(const Val &d, const Ix &r) {
  Val o = vobj();
  for (std::size_t j = 0; j < d.a.size(); ++j) {
    Val a = varr();
    a.a.reserve(r.size());
    for (std::size_t i : r) a.a.push_back(d.a[j].a[i]);
    put(o, d.k[j].c_str(), std::move(a));
  }
  return o;
}
inline Ix which(const std::vector<int> &k) {
  Ix r;
  for (std::size_t i = 0; i < k.size(); ++i) {
    if (k[i]) r.push_back(i);
  }
  return r;
}
inline Ix all(std::size_t n) {
  Ix r(n);
  for (std::size_t i = 0; i < n; ++i) r[i] = i;
  return r;
}

// Column names: the arguments joined as c() joins them.
struct Nm {
  bool null = true, chr = false, na = false, lst = false;
  Cols v;
};
inline void nmadd(Nm &m, const Val &e) {
  std::string s;
  if (e.t == Val::STR) {
    m.chr = true;
    m.v.push_back(e.s);
  } else if (e.t == Val::ARR || e.t == Val::OBJ) {
    m.lst = true;
  } else if (lab(e, s)) {
    m.v.push_back(s);
  } else {
    m.na = true;
    m.v.emplace_back();
  }
}
inline Nm nms(Ctx &c, std::initializer_list<const char *> as) {
  Nm m;
  for (const char *a : as) {
    const Val &o = c[a];
    if (o.t == Val::NUL) continue;
    m.null = false;
    if (o.t == Val::ARR || o.t == Val::OBJ) {
      for (const Val &e : o.a) nmadd(m, e);
    } else {
      nmadd(m, o);
    }
  }
  if (m.lst) m.chr = false;
  return m;
}
inline Nm nmof(const Cols &v) {
  Nm m;
  m.null = false;
  m.chr = true;
  m.v = v;
  return m;
}
// .axcols(): names of columns of d, nn of them (any number for none),
// present, given once, with no missing values; null: may be left out.
inline bool acols(Ctx &c, const Val &d, const Nm &m, const char *what, std::initializer_list<int> nn, bool null,
                  Cols &out) {
  out.clear();
  if (m.null) {
    if (null) return true;
    return ckstop(c, "`%s` must name columns of `data`, got NULL", what);
  }
  if (!m.chr || m.na || m.v.empty()) return ckstop(c, "`%s` must name columns of `data`", what);
  if (nn.size() && std::find(nn.begin(), nn.end(), static_cast<int>(m.v.size())) == nn.end()) {
    std::string s;
    for (int k : nn) s += (s.empty() ? "" : " or ") + std::to_string(k);
    return ckstop(c, "`%s` must name %s column(s)", what, s.c_str());
  }
  Cols miss;
  for (const std::string &s : m.v) {
    if (fc(d, s) == nullptr && !has(miss, s)) miss.push_back(s);
  }
  if (!miss.empty()) return ckstop(c, "`%s` names columns not in `data`: %s", what, join(miss, ", ").c_str());
  for (std::size_t i = 0; i < m.v.size(); ++i) {
    for (std::size_t j = 0; j < i; ++j) {
      if (m.v[i] == m.v[j]) return ckstop(c, "`%s` names a column twice", what);
    }
  }
  for (const std::string &s : m.v) {
    for (const Val &e : fc(d, s)->a) {
      if (vna(e)) return ckstop(c, "column `%s` has missing values", s.c_str());
    }
  }
  out = m.v;
  return true;
}
inline bool acol(Ctx &c, const Val &d, const char *nm, std::initializer_list<int> nn, bool null, Cols &out) {
  return acols(c, d, nms(c, {nm}), nm, nn, null, out);
}
// .axmat(): the columns k of d, numeric or logical and finite, as a
// column-major matrix.
inline bool amat(Ctx &c, const Val &d, const Cols &k, std::vector<double> &m) {
  const std::size_t n = nr(d);
  m.assign(n * k.size(), 0.0);
  for (std::size_t j = 0; j < k.size(); ++j) {
    const Val &a = *fc(d, k[j]);
    const char *cl = ccls(a);
    if (cl[0] == 'c') return ckstop(c, "column `%s` must be numeric, got %s", k[j].c_str(), cl);
    for (std::size_t i = 0; i < n; ++i) {
      const double x = nv(a.a[i]);
      if (!std::isfinite(x)) return ckstop(c, "column `%s` must be finite", k[j].c_str());
      m[j * n + i] = x;
    }
  }
  return true;
}
// The same, the columns known to be fine.
inline std::vector<double> matv(const Val &d, const Cols &k) {
  const std::size_t n = nr(d);
  std::vector<double> m(n * k.size());
  for (std::size_t j = 0; j < k.size(); ++j) {
    const Val &a = *fc(d, k[j]);
    for (std::size_t i = 0; i < n; ++i) m[j * n + i] = nv(a.a[i]);
  }
  return m;
}

// .axcode(): codes from 0 of labels, in order of first appearance.
struct Code {
  std::vector<int> c;
  int n = 0;
  Cols lev;
};
inline Code acode(const Val &a) {
  Code r;
  std::unordered_map<std::string, int> h;
  r.c.reserve(a.a.size());
  for (const Val &e : a.a) {
    const std::string s = vlab(e);
    const auto it = h.emplace(s, r.n);
    if (it.second) {
      r.lev.push_back(s);
      ++r.n;
    }
    r.c.push_back(it.first->second);
  }
  return r;
}
// Codes v (from 0, ng of them) of the rows r, recoded from 0 in order of
// first appearance among them, as .axcode() of the rows' labels.
inline std::vector<int> recode(const std::vector<int> &v, int ng, const Ix &r, int &g) {
  std::vector<int> map(zz(ng), -1), o(r.size());
  g = 0;
  for (std::size_t i = 0; i < r.size(); ++i) {
    int &m = map[zz(v[r[i]])];
    if (m < 0) m = g++;
    o[i] = m;
  }
  return o;
}
inline const double *pd(const std::vector<double> &v) {
  static const double z = 0.0;
  return v.empty() ? &z : v.data();
}

// The message of a kernel's failure, as the R entry point of its name
// words it (src/bapx.cpp).
inline bool kerr(Ctx &c, const char *who, int rc) {
  c.err = std::string(who) + ": " +
          (rc == CLAB_EDOM ? std::string("the fit failed: a regressor or an instrument does not vary within the fixed "
                                         "effects (they absorb it), or the columns are collinear.")
                           : std::string(clab_strerror(rc)));
  return false;
}

// Another function of the family, called as the R wrappers call it.
using Fn = bool (*)(Ctx &, Val &);
inline bool call(Ctx &c, const char *fn, Fn f, const Val &a, Val &out) {
  Ctx s{fn, &a, std::string(), {}};
  const bool ok = f(s, out);
  for (const std::string &w : s.warn) c.warn.push_back(w);
  if (!ok) c.err = s.err;
  return ok;
}

// .axlevel()
inline bool alevel(Ctx &c, double &lv) {
  lv = val1d(c, "level", 0.95);
  if (!(std::isfinite(lv) && lv > 0 && lv < 1)) {
    return ckstop(c, "`level` must be one number strictly between 0 and 1");
  }
  return true;
}

// The Poisson fit (.axppfit): the numbers of all the rows of a frame.
// .axfe(): the effects as codes, the sets after the second in fx.
struct Fe {
  bool two = false;
  Code f1, f2;
  std::vector<Code> fx;
};
inline Fe afe(const Val &d, const Cols &fe) {
  Fe a;
  a.f1 = acode(*fc(d, fe[0]));
  a.two = fe.size() >= 2;
  if (a.two) a.f2 = acode(*fc(d, fe[1]));
  for (std::size_t j = 2; j < fe.size(); ++j) a.fx.push_back(acode(*fc(d, fe[j])));
  return a;
}
// The effects of the rows r, coded from 0 over them as .axfe() of those
// rows codes them: f1 and f2, and the further sets as an n x mx matrix cx
// (column-major) with their counts nx.
struct Fr {
  std::vector<int> c1, c2, cx, nx;
  int n1 = 0, n2 = 0;
  const int *p2() const { return c2.empty() ? nullptr : c2.data(); }
  const int *px() const { return cx.empty() ? nullptr : cx.data(); }
  const int *pn() const { return nx.empty() ? nullptr : nx.data(); }
  int mx() const { return static_cast<int>(nx.size()); }
  // the doubles of scratch the kernels ask for each level of every set
  std::size_t ng() const {
    std::size_t g = zz(n1 + n2);
    for (int v : nx) g += zz(v);
    return g;
  }
};
inline Fr frows(const Fe &a, const Ix &r) {
  Fr f;
  f.c1 = recode(a.f1.c, a.f1.n, r, f.n1);
  if (a.two) f.c2 = recode(a.f2.c, a.f2.n, r, f.n2);
  for (const Code &x : a.fx) {
    int g = 0;
    const std::vector<int> v = recode(x.c, x.n, r, g);
    f.cx.insert(f.cx.end(), v.begin(), v.end());
    f.nx.push_back(g);
  }
  return f;
}
// .axfesep(): clab_fesep on the first two sets, then on each further set
// alone; keep, one per row of r. Negative for a code out of range.
inline int fesep(const std::vector<double> &y, const Fr &f, std::vector<int> &keep) {
  const std::size_t n = y.size();
  keep.assign(n, 0);
  std::vector<double> wk(f.ng() + 1);
  int rc = clab_fesep(pd(y), static_cast<int>(n), f.c1.data(), f.n1, f.p2(), f.n2, keep.data(), wk.data());
  std::vector<int> k2(n);
  for (std::size_t j = 0; rc >= 0 && j < f.nx.size(); ++j) {
    rc = clab_fesep(pd(y), static_cast<int>(n), &f.cx[j * n], f.nx[j], nullptr, 0, k2.data(), wk.data());
    for (std::size_t i = 0; i < n; ++i) keep[i] = keep[i] && k2[i];
  }
  return rc;
}

struct Pd {
  std::size_t n = 0;
  int p = 0;
  std::vector<double> y, x, off;
  bool hoff = false;
  Fe fe;
};
inline Pd pdat(const Val &d, const std::string &y, const Cols &x, const Cols &fe, const Cols &expo) {
  Pd a;
  a.n = nr(d);
  a.p = static_cast<int>(x.size());
  a.y = colv(*fc(d, y));
  a.x = matv(d, x);
  if (!expo.empty()) {
    a.hoff = true;
    for (double v : colv(*fc(d, expo[0]))) a.off.push_back(std::log(v));
  }
  a.fe = afe(d, fe);
  return a;
}
struct Pf {
  std::vector<double> beta, st, sco, brd;
};
// The fit on the rows r (the effects' codes made over them), with
// frequency weights w (nullptr or one per row of r).
inline int ppfit(const Pd &a, const Ix &r, const double *w, double tol, int maxit, Pf &f) {
  const std::size_t n = r.size(), sp = zz(a.p);
  if (n > static_cast<std::size_t>(INT_MAX)) return CLAB_EARG;
  std::vector<double> y(n), x(n * sp), off(a.hoff ? n : 0);
  for (std::size_t i = 0; i < n; ++i) {
    y[i] = a.y[r[i]];
    for (std::size_t j = 0; j < sp; ++j) x[j * n + i] = a.x[j * a.n + r[i]];
    if (a.hoff) off[i] = a.off[r[i]];
  }
  const Fr e = frows(a.fe, r);
  // st is 7 long, as src/bapx.cpp gives it: the kernel writes Pearson's
  // sum last
  f.beta.assign(sp, 0.0);
  f.st.assign(7, 0.0);
  f.sco.assign(n * sp, 0.0);
  f.brd.assign(sp * sp, 0.0);
  std::vector<double> vc(sp * sp, NA);
  std::vector<double> wk(n * (2 * sp + 8) + 3 * sp * sp + sp + 4 * e.ng() + 1);
  std::vector<int> iw(n + 1);
  return clab_ppmlfek(pd(y), pd(x), static_cast<int>(n), a.p, a.hoff ? off.data() : nullptr, w, e.c1.data(), e.n1,
                      e.p2(), e.n2, e.px(), e.pn(), e.mx(), nullptr, 0, nullptr, 0, maxit, 1e-10, tol,
                      f.beta.data(), vc.data(), f.sco.data(), f.brd.data(), f.st.data(), wk.data(), iw.data());
}
// .axfesep() on the rows r: keep, one per row of r.
inline int psep(const Pd &a, const Ix &r, std::vector<int> &keep) {
  const std::size_t n = r.size();
  std::vector<double> y(n);
  for (std::size_t i = 0; i < n; ++i) y[i] = a.y[r[i]];
  return fesep(y, frows(a.fe, r), keep);
}
// .axsep(): the rows kept after separation by the effects.
inline bool asep(Ctx &c, const std::vector<double> &y, const Val &d, const Cols &fe, std::vector<int> &k) {
  if (fesep(y, frows(afe(d, fe), all(y.size())), k) < 0) {
    c.err = std::string("ppmlfe: ") + clab_strerror(CLAB_EARG);
    return false;
  }
  std::size_t kept = 0;
  for (int v : k) kept += v ? 1 : 0;
  if (kept == 0) return ckstop(c, "every group of the fixed effects has no counts");
  if (kept < k.size()) {
    wrn(c, std::to_string(k.size() - kept) + " row(s) dropped: a group of `" + join(fe, "` or `") +
               "` with no counts at all (its effect would be minus infinity; it carries nothing for the other "
               "coefficients).");
  }
  return true;
}

// .axcl(): the clusters, refused for one, a warning under ten.
struct Cl {
  Code a, b;
  bool two = false;
  std::vector<double> g;
};
inline bool acl(Ctx &c, const Val &d, const Cols &clus, Cl &r) {
  r.a = acode(*fc(d, clus[0]));
  r.two = clus.size() == 2;
  if (r.two) r.b = acode(*fc(d, clus[1]));
  r.g = {static_cast<double>(r.a.n)};
  if (r.two) r.g.push_back(r.b.n);
  for (double g : r.g) {
    if (g < 2) return ckstop(c, "`clus` has one cluster: a cluster-robust variance needs two or more");
  }
  for (double g : r.g) {
    if (g < 10) {
      wrn(c, "fewer than ten clusters; cluster-robust errors are unreliable with so few (wcbx() is made for this "
             "case).");
      break;
    }
  }
  return true;
}
// .axnest(): each label of column f falls in one label of column cc.
inline bool nest(const Val &d, const std::string &f, const std::string &cc) {
  const Val &a = *fc(d, f), &b = *fc(d, cc);
  std::unordered_map<std::string, std::string> h;
  for (std::size_t i = 0; i < a.a.size(); ++i) {
    const std::string l = vlab(b.a[i]);
    const auto it = h.emplace(vlab(a.a[i]), l);
    if (!it.second && it.first->second != l) return false;
  }
  return true;
}
// .axtime()
inline bool atime(Ctx &c, const Val &d, std::vector<double> &tv, const char *se = "nw") {
  const Nm m = nms(c, {"time"});
  if (m.null) return ckstop(c, "`se = \"%s\"` needs `time`, the name of the day column", se);
  Cols t;
  CK(acols(c, d, m, "time", {1}, false, t));
  const Val &a = *fc(d, t[0]);
  if (!isnumc(a)) return ckstop(c, "column `%s` must hold numbers or dates", t[0].c_str());
  tv = colv(a);
  return true;
}
// .axnwlag(): the bandwidth of Newey and West (1994) for the Bartlett
// kernel from the rows' scores sco (n x p), the units u (codes, ng of
// them) and the days tv.
inline double nwlag(const std::vector<double> &sco, std::size_t n, int p, const std::vector<int> &u, int ng,
                    const std::vector<double> &tv) {
  // rowSums() and sum() add in long double
  std::vector<double> h(n);
  for (std::size_t i = 0; i < n; ++i) {
    long double s = 0.0;
    for (std::size_t j = 0; j < zz(p); ++j) s += sco[j * n + i];
    h[i] = static_cast<double>(s);
  }
  std::vector<int> tab(zz(ng), 0);
  for (int v : u) ++tab[zz(v)];
  const double tt = *std::max_element(tab.begin(), tab.end());
  const int nn = static_cast<int>(std::floor(4 * std::pow(tt / 100, 2.0 / 9)));
  // match(paste(u, tv - j), paste(u, tv)): the first row of a unit and day
  std::map<std::pair<int, double>, std::size_t> key;
  for (std::size_t i = 0; i < n; ++i) key.emplace(std::make_pair(u[i], tv[i]), i);
  std::vector<double> sg;
  for (int j = 0; j <= nn; ++j) {
    long double s = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      const auto it = key.find(std::make_pair(u[i], tv[i] - j));
      if (it != key.end()) s += static_cast<long double>(h[i]) * h[it->second];
    }
    sg.push_back(static_cast<double>(s));
  }
  long double a = 0.0, b = 0.0;
  for (int j = 1; j <= nn; ++j) {
    a += static_cast<long double>(j) * sg[zz(j)];
    b += sg[zz(j)];
  }
  const double s1 = 2 * static_cast<double>(a), s0 = sg[0] + 2 * static_cast<double>(b);
  if (!std::isfinite(s1 / s0)) return 0;
  return std::floor(1.1447 * std::pow(rpow(s1 / s0, 2), 1.0 / 3) * std::pow(tt, 1.0 / 3));
}

// .axvar(): the variance from the rows' scores as se asks (0 cluster, 1
// twoway, 2 hc1, 3 nw, 5 quasi, 6 dk), the small-sample factor by ssc (0 G, 1
// nk, 2 none), the two-way factors by twc (0 min, 1 each) and psd (1 eig:
// the negative eigenvalues of a two-way variance set to 0); the clusters
// counted, the lag and, for quasi, the dispersion psi from Pearson's sum
// pr. time and lag are the caller's.
struct Vr {
  std::vector<double> v, g;
  double lag = NA, psi = NA;
};
inline bool avar(Ctx &c, const Val &d, int p, const std::vector<double> &sco, const std::vector<double> &brd, int se,
                 const Cols &clus, const Cols &fe, int ssc, int twc, int psd, double pr, Vr &r) {
  const std::size_t n = nr(d);
  auto kk = [&](const Cols &cls) {
    double k = p + 1;
    for (const std::string &f : fe) {
      bool any = false;
      for (const std::string &cc : cls) any = any || nest(d, f, cc);
      if (!any) k += static_cast<double>(acode(*fc(d, f)).n) - 1;
    }
    return k;
  };
  std::vector<int> c1, c2;
  int n1 = 0, n2 = 0, kind = 0, lg = 0;
  std::vector<double> tv;
  double k = 0;
  if (se == 5) {
    // the model-based (X~'WX~)^-1 times the Pearson dispersion
    k = kk(Cols());
    r.psi = pr / (static_cast<double>(n) - k);
    r.v.resize(brd.size());
    for (std::size_t a = 0; a < brd.size(); ++a) r.v[a] = brd[a] * r.psi;
    r.g = {static_cast<double>(n)};
    r.lag = NA;
    return true;
  }
  if (se <= 1) {
    Cl cl;
    CK(acl(c, d, clus, cl));
    c1 = cl.a.c;
    n1 = cl.a.n;
    if (cl.two) {
      c2 = cl.b.c;
      n2 = cl.b.n;
    }
    r.g = cl.g;
    kind = se == 1 ? 1 : 0;
    k = kk(clus);
  } else {
    r.g = {static_cast<double>(n)};
    kind = se == 2 ? 2 : (se == 3 ? 3 : 4);
    k = kk(Cols());
    if (se == 6) {
      // .axdklag(): the bandwidth of .axnwlag() on the series of the scores
      // summed over the rows of each time
      CK(atime(c, d, tv, "dk"));
      std::vector<double> ut = tv;
      std::sort(ut.begin(), ut.end());
      ut.erase(std::unique(ut.begin(), ut.end()), ut.end());
      r.g = {static_cast<double>(ut.size())};
      if (c.has("lag")) {
        Nv q;
        CK(cknum(c, c["lag"], "lag", q, true, 0, Inf, false, true));
        lg = q[0] > INT_MAX ? INT_MAX : static_cast<int>(q[0]);
      } else {
        const std::size_t nt = ut.size();
        std::vector<double> h(nt * zz(p), 0.0);
        for (std::size_t i = 0; i < n; ++i) {
          const std::size_t t = static_cast<std::size_t>(std::lower_bound(ut.begin(), ut.end(), tv[i]) - ut.begin());
          for (std::size_t j = 0; j < zz(p); ++j) h[j * nt + t] += sco[j * n + i];
        }
        const double l = nwlag(h, nt, p, std::vector<int>(nt, 0), 1, ut);
        lg = l > INT_MAX ? INT_MAX : static_cast<int>(l);
      }
    }
    if (se == 3) {
      const Code u = acode(*fc(d, clus[0]));
      c1 = u.c;
      n1 = u.n;
      CK(atime(c, d, tv));
      const Val lv = c.has("lag") ? c["lag"] : vnum(nwlag(sco, n, p, u.c, u.n, tv));
      Nv q;
      CK(cknum(c, lv, "lag", q, true, 0, Inf, false, true));
      lg = q[0] > INT_MAX ? INT_MAX : static_cast<int>(q[0]);
    }
  }
  const double dn = static_cast<double>(n);
  auto ff = [&](double gg) { return ssc == 0 ? gg / (gg - 1) : ssc == 1 ? gg / (gg - 1) * (dn - 1) / (dn - k) : 1.0; };
  const double gm = *std::min_element(r.g.begin(), r.g.end());
  auto vcx = [&](int kd, const std::vector<int> &a1, int m1, const std::vector<int> &a2, int m2, double fac,
                 std::vector<double> &v) {
    v.assign(zz(p) * zz(p), 0.0);
    const int rc = clab_vcx(pd(sco), static_cast<int>(n), p, pd(brd), kd, a1.empty() ? nullptr : a1.data(), m1,
                            a2.empty() ? nullptr : a2.data(), m2, tv.empty() ? nullptr : tv.data(), lg, fac,
                            v.data());
    if (rc == CLAB_EDOM) {
      c.err = "vcx: the time column has a value that is not finite.";
      return false;
    }
    if (rc != CLAB_OK) {
      c.err = std::string("vcx: ") + clab_strerror(rc);
      return false;
    }
    return true;
  };
  if (se == 1 && twc == 1) {
    // Cameron, Gelbach and Miller (2011, p. 241): each of the three
    // components with its own factor, G, H and I (the intersections,
    // coded as .axcode(paste(c1, c2)) codes them)
    std::map<std::pair<int, int>, int> h;
    std::vector<int> c12(n);
    for (std::size_t i = 0; i < n; ++i) {
      c12[i] = h.emplace(std::make_pair(c1[i], c2[i]), static_cast<int>(h.size())).first->second;
    }
    const int n12 = static_cast<int>(h.size());
    std::vector<double> v1, v2, v3;
    const std::vector<int> none;
    CK(vcx(0, c1, n1, none, 0, ff(n1), v1));
    CK(vcx(0, c2, n2, none, 0, ff(n2), v2));
    CK(vcx(0, c12, n12, none, 0, ff(n12), v3));
    r.v.resize(v1.size());
    for (std::size_t a = 0; a < v1.size(); ++a) r.v[a] = v1[a] + v2[a] - v3[a];
  } else {
    CK(vcx(kind, c1, n1, c2, n2, ff(gm), r.v));
  }
  if (se == 1 && psd == 1) {
    // the eigenvalues below 0 set to 0 (Cameron, Gelbach and Miller 2011,
    // p. 241-242, after Politis 2007)
    svk::Mat m(p, p), e;
    m.v = r.v;
    std::vector<double> ev;
    svk::eigsym(m, ev, e);
    for (int a = 0; a < p; ++a) {
      for (int b = 0; b < p; ++b) {
        double s = 0.0;
        for (int l = 0; l < p; ++l) s += e(a, l) * std::max(ev[zz(l)], 0.0) * e(b, l);
        r.v[zz(b) * zz(p) + zz(a)] = s;
      }
    }
  }
  r.lag = se == 3 || se == 6 ? lg : NA;
  return true;
}
// sqrt(pmax(diag(v), 0))
inline std::vector<double> sdiag(const std::vector<double> &v, int p) {
  std::vector<double> s;
  for (std::size_t j = 0; j < zz(p); ++j) s.push_back(std::sqrt(std::max(v[j * zz(p) + j], 0.0)));
  return s;
}
// .axtab(): t on df degrees of freedom, the normal for df NaN; rr adds
// exp(est).
inline Val atab(const Cols &tm, const std::vector<double> &b, const std::vector<double> &se, double lv, double df,
                bool rr) {
  const bool t = !std::isnan(df);
  const double q = t ? clab_qt(0.5 + lv / 2, df, 1) : clab_qnorm(0.5 + lv / 2);
  std::vector<double> st, p, lo, hi, ex;
  for (std::size_t j = 0; j < b.size(); ++j) {
    const double s = b[j] / se[j];
    st.push_back(s);
    p.push_back(2 * (t ? clab_pt(-std::fabs(s), df, 1) : clab_pnorm(-std::fabs(s))));
    lo.push_back(b[j] - q * se[j]);
    hi.push_back(b[j] + q * se[j]);
    ex.push_back(std::exp(b[j]));
  }
  Frame f;
  f.str("term", tm);
  f.num("est", b);
  f.num("se", se);
  f.num("stat", st);
  f.num("p", p);
  f.num("lo", lo);
  f.num("hi", hi);
  if (rr) f.num("rr", ex);
  return f.o;
}

// The cluster bootstrap. nboot and seed as R reads them, nboot first.
inline bool bargs(Ctx &c, int &nb, bool &hs, double &sd) {
  Nv q;
  CK(cknumd(c, "nboot", 199, q, 2, Inf, false, true));
  nb = q[0] > INT_MAX ? INT_MAX : static_cast<int>(q[0]);
  hs = c.has("seed");
  if (hs) {
    CK(cknum(c, c["seed"], "seed", q, true, 0, Inf, false, true));
    sd = q[0];
  }
  return true;
}
// .axdraws(): how often each of g clusters is drawn in each draw,
// tabulate(sample.int(g, g, replace = TRUE)) after set.seed(seed), made a
// block of draws at a time (next(): k draws, g x k, column-major) so that
// the g x nboot counts are never held at once.
struct Adraw {
  igh::Rng rng;
  int g;
  std::size_t blk;  // draws per block: about 2^22 counts
  Adraw(std::uint32_t s, int gg)
      : rng(s), g(gg), blk(std::max<std::size_t>(1, (static_cast<std::size_t>(1) << 22) / zz(gg))) {}
  void next(std::size_t k, std::vector<int> &bm) {
    bm.assign(zz(g) * k, 0);
    for (std::size_t j = 0; j < k; ++j) {
      for (int i = 0; i < g; ++i) ++bm[j * zz(g) + static_cast<std::size_t>(rng.index(g))];
    }
  }
};
inline bool adseed(Ctx &c, bool hs, double sd, std::uint32_t &s) {
  if (!hs) return c.fail("`seed` is needed outside R: one whole number from 0.");
  if (sd > 2147483647.0) {
    c.err = "supplied seed is not a valid integer";
    return false;
  }
  s = static_cast<std::uint32_t>(static_cast<std::int32_t>(sd));
  return true;
}
// apply(x, 2, sd, na.rm = TRUE) of an n x p matrix.
inline std::vector<double> csd(const std::vector<double> &x, std::size_t n, int p) {
  std::vector<double> s;
  for (std::size_t j = 0; j < zz(p); ++j) {
    std::vector<double> col(x.begin() + static_cast<std::ptrdiff_t>(j * n),
                            x.begin() + static_cast<std::ptrdiff_t>((j + 1) * n));
    std::vector<char> ok(n);
    for (std::size_t i = 0; i < n; ++i) ok[i] = !std::isnan(col[i]);
    s.push_back(std::sqrt(rcov(col, n, 1, ok)[0]));
  }
  return s;
}

// .axboot(): fit(rows, weights, beta) on each draw of the clusters of
// column cl, false for a draw that cannot be fitted (left out). dr: nb x
// p, column-major, NA for those.
struct Bt {
  std::vector<double> dr, se;
  std::vector<char> ok;
  int nf = 0, g = 0;
};
template <class F>
bool aboot(Ctx &c, const Val &d, const std::string &cl, int nb, bool hs, double sd, int p, F fit, Bt &r) {
  const Code cc = acode(*fc(d, cl));
  if (cc.n < 2) return ckstop(c, "the bootstrap needs two or more clusters in `%s`", cl.c_str());
  std::uint32_t s0 = 0;
  CK(adseed(c, hs, sd, s0));
  Adraw ad(s0, cc.n);
  std::vector<int> bm;
  const std::size_t n = cc.c.size(), snb = zz(nb), sp = zz(p);
  r.g = cc.n;
  r.dr.assign(snb * sp, NA);
  Ix rw;
  std::vector<double> w, beta;
  for (std::size_t b0 = 0; b0 < snb; b0 += ad.blk) {
    const std::size_t k = std::min(ad.blk, snb - b0);
    ad.next(k, bm);
    for (std::size_t j = 0; j < k; ++j) {
      const std::size_t b = b0 + j;
      rw.clear();
      w.clear();
      for (std::size_t i = 0; i < n; ++i) {
        const int m = bm[j * zz(cc.n) + zz(cc.c[i])];
        if (m > 0) {
          rw.push_back(i);
          w.push_back(m);
        }
      }
      if (fit(rw, w, beta)) {
        for (std::size_t q = 0; q < sp; ++q) r.dr[q * snb + b] = beta[q];
      }
    }
  }
  r.ok.assign(snb, 1);
  r.nf = 0;
  for (std::size_t b = 0; b < snb; ++b) {
    for (std::size_t j = 0; j < sp; ++j) r.ok[b] = r.ok[b] && !std::isnan(r.dr[j * snb + b]);
    r.nf += r.ok[b] ? 0 : 1;
  }
  if (r.nf > 0) {
    wrn(c, std::to_string(r.nf) + " of " + std::to_string(nb) +
               " bootstrap draws could not be fitted and are left out.");
  }
  if (nb - r.nf < 2) return ckstop(c, "fewer than two bootstrap draws could be fitted");
  r.se = csd(r.dr, snb, p);
  return true;
}


// Two-stage least squares (.axivfit): the numbers of all the rows.
struct Id {
  std::size_t n = 0;
  int pe = 0, px = 0, q = 0;
  std::vector<double> y, xe, xx, zi;
  Fe fe;
};
inline Id idat(const Val &d, const std::string &y, const Cols &endo, const Cols &x, const Cols &inst, const Cols &fe) {
  Id a;
  a.n = nr(d);
  a.pe = static_cast<int>(endo.size());
  a.px = static_cast<int>(x.size());
  a.q = static_cast<int>(inst.size());
  a.y = colv(*fc(d, y));
  a.xe = matv(d, endo);
  a.xx = matv(d, x);
  a.zi = matv(d, inst);
  a.fe = afe(d, fe);
  return a;
}
struct If {
  std::vector<double> beta, vc, fst, st, res, fres, sco, brd;
};
inline std::vector<double> sub(const std::vector<double> &m, std::size_t n, int p, const Ix &r) {
  std::vector<double> o(r.size() * zz(p));
  for (std::size_t j = 0; j < zz(p); ++j) {
    for (std::size_t i = 0; i < r.size(); ++i) o[j * r.size() + i] = m[j * n + r[i]];
  }
  return o;
}
// The fit on the rows r, weights w (nullptr or one per row of r), the
// first-stage clusters c1 and c2 (codes of the rows r, or nullptr).
inline int ivfit(const Id &a, const Ix &r, const double *w, const int *c1, int n1, const int *c2, int n2, double tol,
                 If &f) {
  const std::size_t n = r.size(), p = zz(a.pe + a.px), pe = zz(a.pe);
  if (n > static_cast<std::size_t>(INT_MAX)) return CLAB_EARG;
  const std::vector<double> y = sub(a.y, a.n, 1, r), xe = sub(a.xe, a.n, a.pe, r), xx = sub(a.xx, a.n, a.px, r),
                            zi = sub(a.zi, a.n, a.q, r);
  const Fr e = frows(a.fe, r);
  // fst is pe x 4, as src/bapx.cpp gives it: the Wald statistic, the rk
  // statistic, the effective F and its degrees of freedom
  f.beta.assign(p, 0.0);
  f.vc.assign(p * p, NA);
  f.fst.assign(pe * 4, NA);
  f.st.assign(3, 0.0);
  f.res.assign(n, 0.0);
  f.fres.assign(n * pe, 0.0);
  f.sco.assign(n * p, 0.0);
  f.brd.assign(p * p, 0.0);
  double z[3];
  auto po = [&](std::vector<double> &v) { return v.empty() ? z : v.data(); };
  return clab_ivfek(pd(y), pd(xe), a.pe, pd(xx), a.px, pd(zi), a.q, static_cast<int>(n), w, e.c1.data(), e.n1,
                    e.p2(), e.n2, e.px(), e.pn(), e.mx(), c1, n1, c2, n2, tol, po(f.beta), po(f.vc), po(f.fst),
                    po(f.fres), po(f.res), po(f.sco), po(f.brd), po(f.st));
}
// qchisq(p, df, ncp) of the lower tail, as R's qnchisq() finds it: from
// Pearson's approximation, bounds by doubling and halving, then bisection
// to a relative width of 1e-13.
inline double qnchisq(double p, double df, double ncp) {
  const double eps = 1e-11, accu = 1e-13, big = std::numeric_limits<double>::max(),
               dmin = std::numeric_limits<double>::min(), deps = std::numeric_limits<double>::epsilon();
  const double b = ncp * ncp / (df + 3 * ncp), cc = (df + 3 * ncp) / (df + 2 * ncp), ff = (df + 2 * ncp) / (cc * cc);
  double ux = b + cc * clab_qchisq(p, ff, 1);
  if (ux < 0) ux = 1;
  const double ux0 = ux;
  double pp = std::fmin(1 - deps, p * (1 + eps));
  for (; ux < big && clab_pnchisq(ux, df, ncp, 1) < pp; ux *= 2) {
  }
  pp = p * (1 - eps);
  double lx = std::fmin(ux0, big);
  for (; lx > dmin && clab_pnchisq(lx, df, ncp, 1) > pp; lx *= 0.5) {
  }
  double nx = 0;
  do {
    nx = 0.5 * (lx + ux);
    if (clab_pnchisq(nx, df, ncp, 1) > p) ux = nx;
    else lx = nx;
  } while ((ux - lx) / nx > accu);
  return 0.5 * (ux + lx);
}
// .axfirst()
inline Val afirst(const Cols &endo, const std::vector<double> &fst, int q) {
  const std::size_t pe = endo.size();
  std::vector<double> wd, kp, ef, ke, cv, ni(pe, q);
  for (std::size_t j = 0; j < pe; ++j) {
    wd.push_back(fst[j]);
    kp.push_back(fst[pe + j]);
    ef.push_back(fst[2 * pe + j]);
    // the simplified test of Montiel Olea and Pflueger (2013, section
    // 2.2.2 and table 1): the upper 5% quantile of a noncentral
    // chi-squared on keff degrees of freedom with noncentrality x keff,
    // over keff, x = 1 / tau = 10
    const double k = fst[3 * pe + j];
    ke.push_back(k);
    cv.push_back(std::isfinite(k) && k > 0 ? qnchisq(0.95, k, 10 * k) / k : NA);
  }
  Frame f;
  f.str("endo", endo);
  f.num("wald", wd);
  f.num("kp", kp);
  f.num("eff", ef);
  f.num("keff", ke);
  f.num("cv", cv);
  f.num("ninst", ni);
  return f.o;
}
inline Val vs(const Cols &v) { return unbox(vstrs(v)); }
// The rows of a coefficient table whose term is in tm.
inline Val terms(const Val &tab, const Cols &tm) {
  const Val &t = col(tab, "term");
  Ix r;
  for (std::size_t i = 0; i < t.a.size(); ++i) {
    if (t.a[i].t == Val::STR && has(tm, t.a[i].s)) r.push_back(i);
  }
  return rows(tab, r);
}

const std::vector<const char *> sscs = {"G", "nk", "none"}, twcs = {"min", "each"}, psds = {"none", "eig"};

// The rank and the pivot of qr(x) (n x p, column-major, overwritten) with
// its default tolerance: LINPACK's dqrdc2 as R changes it, Householder
// reflections with limited pivoting, a column whose norm falls below tol
// of its first norm moved to the end.
inline int qrank(std::vector<double> &x, int n, int p, double tol, std::vector<int> &piv) {
  auto at = [&](int i, int j) -> double & { return x[zz(j) * zz(n) + zz(i)]; };
  auto nrm = [&](int i0, int j) {
    double s = 0.0;
    for (int i = i0; i < n; ++i) s += at(i, j) * at(i, j);
    return std::sqrt(s);
  };
  std::vector<double> qa(zz(p)), w1(zz(p)), w2(zz(p));
  piv.resize(zz(p));
  for (int j = 0; j < p; ++j) {
    piv[zz(j)] = j;
    qa[zz(j)] = w1[zz(j)] = nrm(0, j);
    w2[zz(j)] = qa[zz(j)] == 0 ? 1.0 : qa[zz(j)];
  }
  // l and k from 1 as in the Fortran
  int k = p + 1;
  for (int l = 1; l <= std::min(n, p); ++l) {
    const int c = l - 1;
    while (l < k && qa[zz(c)] < w2[zz(c)] * tol) {
      for (int i = 0; i < n; ++i) {
        const double t = at(i, c);
        for (int j = c + 1; j < p; ++j) at(i, j - 1) = at(i, j);
        at(i, p - 1) = t;
      }
      std::rotate(piv.begin() + c, piv.begin() + c + 1, piv.end());
      std::rotate(qa.begin() + c, qa.begin() + c + 1, qa.end());
      std::rotate(w1.begin() + c, w1.begin() + c + 1, w1.end());
      std::rotate(w2.begin() + c, w2.begin() + c + 1, w2.end());
      --k;
    }
    if (l == n) continue;
    double nx = nrm(c, c);
    if (nx == 0) continue;
    if (at(c, c) != 0) nx = std::copysign(nx, at(c, c));
    const double sc = 1.0 / nx;
    for (int i = c; i < n; ++i) at(i, c) *= sc;
    at(c, c) += 1.0;
    for (int j = c + 1; j < p; ++j) {
      double t = 0.0;
      for (int i = c; i < n; ++i) t += at(i, c) * at(i, j);
      t = -t / at(c, c);
      for (int i = c; i < n; ++i) at(i, j) += t * at(i, c);
      if (qa[zz(j)] == 0) continue;
      const double tt = std::max(1.0 - rpow(std::fabs(at(c, j)) / qa[zz(j)], 2), 0.0);
      if (std::fabs(tt) < 1e-6) {
        qa[zz(j)] = w1[zz(j)] = nrm(l, j);
      } else {
        qa[zz(j)] *= std::sqrt(tt);
      }
    }
    qa[zz(c)] = at(c, c);
    at(c, c) = -nx;
  }
  return std::min(k - 1, n);
}
// .axsst(): the regressors to drop by the check of Santos Silva and
// Tenreyro (2010, p. 311, steps 1 and 2), those qr() finds aliased among
// the regressors demeaned by the effects on the rows with positive counts
// (yv, one per row of d), in their order. A regressor the effects absorb
// on those rows is left by the demeaning as rounding noise, which the
// rank routine would judge against its own size and keep: it is set to
// zero first, where what is left is below 1e-7 of its spread before.
inline bool asst(Ctx &c, const Val &d, const std::vector<double> &yv, const Cols &x, const Cols &fe, double tol,
                 Cols &drop) {
  drop.clear();
  Ix pos;
  for (std::size_t i = 0; i < yv.size(); ++i) {
    if (yv[i] > 0) pos.push_back(i);
  }
  const Val dp = rows(d, pos);
  std::vector<double> m;
  CK(amat(c, dp, x, m));
  const Fr e = frows(afe(dp, fe), all(pos.size()));
  const int n = static_cast<int>(pos.size()), p = static_cast<int>(x.size());
  std::vector<double> s0(zz(p), 0.0);
  for (int j = 0; j < p; ++j) {
    const double *cj = m.data() + zz(j) * zz(n);
    double mu = 0;
    for (int i = 0; i < n; ++i) mu += cj[i];
    mu /= n;
    double ss = 0;
    for (int i = 0; i < n; ++i) ss += (cj[i] - mu) * (cj[i] - mu);
    s0[zz(j)] = std::sqrt(ss);
  }
  std::vector<double> wk(4 * e.ng() + 1);
  const int rc = clab_demean(m.data(), n, p, e.c1.data(), e.n1, e.p2(), e.n2, e.px(), e.pn(), e.mx(), tol,
                             wk.data());
  if (rc != CLAB_OK) {
    c.err = std::string("ppmlfe: ") + clab_strerror(rc);
    return false;
  }
  for (int j = 0; j < p; ++j) {
    double *cj = m.data() + zz(j) * zz(n);
    double ss = 0;
    for (int i = 0; i < n; ++i) ss += cj[i] * cj[i];
    if (std::sqrt(ss) <= 1e-7 * s0[zz(j)]) std::fill(cj, cj + n, 0.0);
  }
  std::vector<int> piv;
  const int rk = qrank(m, n, p, 1e-7, piv);
  std::vector<int> dr(piv.begin() + rk, piv.end());
  std::sort(dr.begin(), dr.end());
  for (int j : dr) drop.push_back(x[zz(j)]);
  return true;
}

bool ppmlfe(Ctx &c, Val &out) {
  static const std::vector<const char *> ses = {"cluster", "twoway", "hc1", "nw", "boot", "quasi"};
  int se = 0, ssc = 0, twc = 0, psd = 0, sep = 0;
  CK(marg(c, "se", ses, se));
  CK(marg(c, "ssc", sscs, ssc));
  CK(marg(c, "twc", twcs, twc));
  CK(marg(c, "psd", psds, psd));
  CK(marg(c, "sep", {"fe", "sst", "none"}, sep));
  Val d;
  CK(dat(c, c["data"], "data", d));
  Cols y, x, fe, expo, clus;
  CK(acol(c, d, "y", {1}, false, y));
  CK(acol(c, d, "x", {}, false, x));
  CK(acol(c, d, "fe", {}, false, fe));
  CK(acol(c, d, "expo", {1}, true, expo));
  Nm cn = nms(c, {"clus"});
  if (cn.null) cn = nmof(se == 1 ? fe : Cols{fe[0]});
  CK(acols(c, d, cn, "clus", {1, 2}, false, clus));
  if (se == 0 && clus.size() == 2) se = 1;
  if (se == 1 && clus.size() != 2) return ckstop(c, "`se = \"twoway\"` needs two `clus` columns");
  double lv = 0;
  CK(alevel(c, lv));
  Nv tl, mx;
  CK(cknumd(c, "tol", 1e-12, tl, -Inf, Inf, true));
  CK(cknumd(c, "maxit", 100, mx, 1, Inf, false, true));
  const double tol = tl[0];
  const int maxit = mx[0] > INT_MAX ? INT_MAX : static_cast<int>(mx[0]);
  std::vector<double> yv, t;
  CK(amat(c, d, y, yv));
  bool any = false;
  for (double v : yv) {
    if (v < 0) return ckstop(c, "column `%s` must hold counts, at least 0", y[0].c_str());
    any = any || v > 0;
  }
  if (!any) return ckstop(c, "column `%s` is all 0", y[0].c_str());
  std::vector<int> k(yv.size(), 1);
  if (sep != 2) CK(asep(c, yv, d, fe, k));
  const Ix kr = which(k);
  const Val ds = rows(d, kr);
  const std::size_t n = nr(ds);
  Cols xdrop;
  if (sep == 1) {
    std::vector<double> yk;
    for (std::size_t i : kr) yk.push_back(yv[i]);
    CK(asst(c, ds, yk, x, fe, tol, xdrop));
    if (xdrop.size() == x.size()) return ckstop(c, "every regressor is collinear on the rows with counts");
    if (!xdrop.empty()) {
      wrn(c, "regressor(s) dropped, collinear with the others on the rows with positive counts (Santos Silva and "
             "Tenreyro 2010): " + join(xdrop, ", ") + ".");
      Cols kx;
      for (const std::string &s : x) {
        if (!has(xdrop, s)) kx.push_back(s);
      }
      x = kx;
    }
  }
  CK(amat(c, ds, x, t));
  if (!expo.empty()) {
    CK(amat(c, ds, expo, t));
    for (double v : t) {
      if (v <= 0) return ckstop(c, "column `%s` must be positive", expo[0].c_str());
    }
  }
  const Pd a = pdat(ds, y[0], x, fe, expo);
  Pf r;
  const int rc = ppfit(a, all(n), nullptr, tol, maxit, r);
  if (rc != CLAB_OK) return kerr(c, "ppmlfe", rc);
  if (r.st[1] != 1) wrn(c, "the fit did not converge in " + rfmt(mx[0]) + " iterations.");
  const int p = static_cast<int>(x.size());
  std::vector<double> v, sev, g;
  double lg = NA, psi = NA;
  Val draws;
  if (se == 4) {
    int nb = 0;
    bool hs = false;
    double sd = 0;
    CK(bargs(c, nb, hs, sd));
    auto fit = [&](const Ix &rw, const std::vector<double> &w, std::vector<double> &beta) {
      Ix r2 = rw;
      std::vector<double> w2 = w;
      if (sep != 2) {
        std::vector<int> kk;
        if (psep(a, rw, kk) < 0) return false;
        r2.clear();
        w2.clear();
        for (std::size_t i = 0; i < rw.size(); ++i) {
          if (kk[i]) {
            r2.push_back(rw[i]);
            w2.push_back(w[i]);
          }
        }
      }
      Pf f;
      if (ppfit(a, r2, pd(w2), tol, maxit, f) != CLAB_OK) return false;
      beta = f.beta;
      return true;
    };
    Bt b;
    CK(aboot(c, ds, clus[0], nb, hs, sd, p, fit, b));
    v = rcov(b.dr, zz(nb), p, b.ok);
    sev = b.se;
    g = {static_cast<double>(b.g)};
    draws = vmat(b.dr, nb, p);
  } else {
    Vr vr;
    CK(avar(c, ds, p, r.sco, r.brd, se, clus, fe, ssc, twc, psd, r.st[6], vr));
    v = vr.v;
    sev = sdiag(v, p);
    g = vr.g;
    lg = vr.lag;
    psi = vr.psi;
  }
  out = vobj();
  put(out, "coef", atab(x, r.beta, sev, lv, NaN, true));
  put(out, "vcov", vmat(v, p, p));
  put(out, "n", vnum(static_cast<double>(n)));
  put(out, "nsep", vnum(static_cast<double>(yv.size() - n)));
  put(out, "xdrop", vs(xdrop));
  put(out, "ncl", vauto(g));
  put(out, "se", vstr(ses[zz(se)]));
  put(out, "ssc", vstr(sscs[zz(ssc)]));
  put(out, "lag", vnum(lg));
  put(out, "psi", vnum(psi));
  put(out, "iter", vnum(r.st[0]));
  put(out, "conv", vbool(r.st[1] == 1));
  put(out, "dev", vnum(r.st[2]));
  put(out, "loglik", vnum(r.st[3]));
  put(out, "draws", draws);
  put(out, "fe", vs(fe));
  put(out, "clus", vs(clus));
  put(out, "level", vnum(lv));
  return true;
}

bool ivfe(Ctx &c, Val &out) {
  static const std::vector<const char *> ses = {"cluster", "twoway", "hc1", "nw", "boot", "dk"};
  int se = 0, ssc = 0, twc = 0, psd = 0;
  CK(marg(c, "se", ses, se));
  CK(marg(c, "ssc", sscs, ssc));
  CK(marg(c, "twc", twcs, twc));
  CK(marg(c, "psd", psds, psd));
  Val d;
  CK(dat(c, c["data"], "data", d));
  Cols y, endo, inst, x, fe, clus;
  CK(acol(c, d, "y", {1}, false, y));
  CK(acol(c, d, "endo", {}, true, endo));
  CK(acol(c, d, "inst", {}, true, inst));
  CK(acol(c, d, "x", {}, true, x));
  CK(acol(c, d, "fe", {}, false, fe));
  Nm cn = nms(c, {"clus"});
  if (cn.null) cn = nmof(se == 1 ? fe : Cols{fe[0]});
  CK(acols(c, d, cn, "clus", {1, 2}, false, clus));
  if (se == 0 && clus.size() == 2) se = 1;
  if (se == 1 && clus.size() != 2) return ckstop(c, "`se = \"twoway\"` needs two `clus` columns");
  double lv = 0;
  CK(alevel(c, lv));
  Nv tl;
  CK(cknumd(c, "tol", 1e-12, tl, -Inf, Inf, true));
  const double tol = tl[0];
  if (inst.size() < endo.size()) {
    return ckstop(c, "`inst` must give at least as many instruments as `endo` has columns");
  }
  if (endo.empty() && !inst.empty()) return ckstop(c, "`inst` needs `endo`: the columns it instruments");
  if (endo.empty() && x.empty()) return ckstop(c, "give `endo` or `x`: there is nothing to estimate");
  Cols tm = endo;
  tm.insert(tm.end(), x.begin(), x.end());
  for (const std::string &s : tm) {
    if (has(inst, s)) return ckstop(c, "a column cannot be both a regressor and an excluded instrument");
  }
  {
    Cols ck = y;
    ck.insert(ck.end(), tm.begin(), tm.end());
    ck.insert(ck.end(), inst.begin(), inst.end());
    std::vector<double> t;
    CK(amat(c, d, ck, t));
  }
  const std::size_t n = nr(d);
  // .axfc(): the first stages' clusters
  std::vector<int> c1, c2;
  int n1 = 0, n2 = 0;
  if (se == 2) {
    for (std::size_t i = 0; i < n; ++i) c1.push_back(static_cast<int>(i));
    n1 = static_cast<int>(n);
  } else {
    const Code ca = acode(*fc(d, clus[0]));
    Code cb;
    if (se == 1) cb = acode(*fc(d, clus[1]));
    if (ca.n < 2 || (se == 1 && cb.n < 2)) {
      return ckstop(c, "`clus` has one cluster: a cluster-robust variance needs two or more");
    }
    c1 = ca.c;
    n1 = ca.n;
    if (se == 1) {
      c2 = cb.c;
      n2 = cb.n;
    }
  }
  const Id a = idat(d, y[0], endo, x, inst, fe);
  If r;
  const int rc = ivfit(a, all(n), nullptr, c1.data(), n1, se == 1 ? c2.data() : nullptr, n2, tol, r);
  if (rc != CLAB_OK) return kerr(c, "ivfe", rc);
  const int p = static_cast<int>(tm.size());
  std::vector<double> v, sev, g;
  double lg = NA, df = NaN;
  Val draws;
  if (se == 4) {
    int nb = 0;
    bool hs = false;
    double sd = 0;
    CK(bargs(c, nb, hs, sd));
    auto fit = [&](const Ix &rw, const std::vector<double> &w, std::vector<double> &beta) {
      If f;
      if (ivfit(a, rw, pd(w), nullptr, 0, nullptr, 0, tol, f) != CLAB_OK) return false;
      beta = f.beta;
      return true;
    };
    Bt b;
    CK(aboot(c, d, clus[0], nb, hs, sd, p, fit, b));
    v = rcov(b.dr, zz(nb), p, b.ok);
    sev = b.se;
    g = {static_cast<double>(b.g)};
    draws = vmat(b.dr, nb, p);
  } else {
    Vr vr;
    CK(avar(c, d, p, r.sco, r.brd, se == 5 ? 6 : se, clus, fe, ssc, twc, psd, NaN, vr));
    v = vr.v;
    sev = sdiag(v, p);
    g = vr.g;
    lg = vr.lag;
    if (se <= 1) df = *std::min_element(g.begin(), g.end()) - 1;
  }
  // the over-identification tests (Yamada 2026, section 7 and table 8,
  // item 6): Sargan and the cluster-robust Hansen J of the two-step
  // estimate, on the clusters of the first stages
  Val overid, gmm;
  if (inst.size() > endo.size() && !endo.empty()) {
    const Fr e = frows(a.fe, all(n));
    std::vector<double> ob(zz(p)), os(3);
    const int ro = clab_ivoid(pd(a.y), pd(a.xe), a.pe, pd(a.xx), a.px, pd(a.zi), a.q, static_cast<int>(n),
                              e.c1.data(), e.n1, e.p2(), e.n2, e.px(), e.pn(), e.mx(), c1.data(), n1,
                              se == 1 ? c2.data() : nullptr, n2, tol, ob.data(), os.data());
    if (ro != CLAB_OK) return kerr(c, "ivfe", ro);
    Frame f;
    f.str("test", {"sargan", "hansen"});
    f.num("stat", {os[0], os[1]});
    f.num("df", {os[2], os[2]});
    f.num("p", {clab_pchisq(os[0], os[2], 0), clab_pchisq(os[1], os[2], 0)});
    overid = f.o;
    gmm = vobj();
    for (std::size_t j = 0; j < tm.size(); ++j) put(gmm, tm[j].c_str(), vnum(ob[j]));
  }
  // the first stages in the exposure's units, per unit of each
  // instrument and over the exposure's mean
  Val fstab;
  if (!endo.empty()) {
    Frame f;
    Val fe1 = varr(), fi = varr(), fb = varr(), fs = varr(), frl = varr();
    for (const std::string &e : endo) {
      Cols xs = inst;
      xs.insert(xs.end(), x.begin(), x.end());
      Val ar = vobj(), s1;
      put(ar, "data", d);
      put(ar, "y", vstr(e));
      put(ar, "x", vstrs(xs));
      put(ar, "fe", vstrs(fe));
      put(ar, "clus", vstrs(clus));
      put(ar, "se", vstr(se == 4 ? "cluster" : ses[zz(se)]));
      put(ar, "ssc", vstr(sscs[zz(ssc)]));
      if (c.has("time")) put(ar, "time", c["time"]);
      if (c.has("lag")) put(ar, "lag", c["lag"]);
      put(ar, "level", vnum(lv));
      put(ar, "tol", vnum(tol));
      put(ar, "twc", vstr(twcs[zz(twc)]));
      put(ar, "psd", vstr(psds[zz(psd)]));
      CK(call(c, "ivfe", ivfe, ar, s1));
      const Val tb = terms(col(s1, "coef"), inst);
      const double mu = lmean(colv(*fc(d, e)));
      const Val &tt = col(tb, "term"), &te = col(tb, "est"), &ts = col(tb, "se");
      for (std::size_t i = 0; i < tt.a.size(); ++i) {
        fe1.a.push_back(vstr(e));
        fi.a.push_back(tt.a[i]);
        fb.a.push_back(te.a[i]);
        fs.a.push_back(ts.a[i]);
        frl.a.push_back(vnum(nv(te.a[i]) / mu));
      }
    }
    f.add("endo", fe1);
    f.add("inst", fi);
    f.add("est", fb);
    f.add("se", fs);
    f.add("rel", frl);
    fstab = f.o;
  }
  out = vobj();
  put(out, "coef", atab(tm, r.beta, sev, lv, df, false));
  put(out, "first", afirst(endo, r.fst, a.q));
  put(out, "fstab", fstab);
  put(out, "overid", overid);
  put(out, "gmm", gmm);
  put(out, "vcov", vmat(v, p, p));
  put(out, "n", vnum(static_cast<double>(n)));
  put(out, "ncl", vauto(g));
  put(out, "se", vstr(ses[zz(se)]));
  put(out, "ssc", vstr(sscs[zz(ssc)]));
  put(out, "lag", vnum(lg));
  put(out, "rss", vnum(r.st[2]));
  put(out, "draws", draws);
  put(out, "fe", vs(fe));
  put(out, "clus", vs(clus));
  put(out, "level", vnum(lv));
  return true;
}

// ppmlcf(family = "nb") (R/cfnb.R): on the rows rw of ds (repeated as
// drawn), the fixed effects as indicators (the first level present left
// out), the first stage by least squares of each endogenous column on an
// intercept, the instruments, the controls and the indicators, and the
// negative binomial fit of clab_cntreg() on the endogenous columns, the
// controls, the first-stage residuals, an intercept and the indicators.
struct Nbf {
  std::vector<double> beta, vb, vc, st;
  int k = 0;
};
inline bool nbcf(const Val &ds, const Ix &rw, const std::string &y, const Cols &endo, const Cols &inst, const Cols &x,
                 const Cols &fe, const std::vector<double> &off, const int *cl, int ncl, int maxit, Nbf &o) {
  const std::size_t m = rw.size();
  auto colr = [&](const std::string &nm) {
    const std::vector<double> v = colv(*fc(ds, nm));
    std::vector<double> r(m);
    for (std::size_t i = 0; i < m; ++i) r[i] = v[rw[i]];
    return r;
  };
  std::vector<double> dm;
  std::size_t nd = 0;
  for (const std::string &f : fe) {
    const Code cf = acode(*fc(ds, f));
    std::vector<char> pres(zz(cf.n), 0);
    for (std::size_t i = 0; i < m; ++i) pres[zz(cf.c[rw[i]])] = 1;
    bool first = true;
    for (int l = 0; l < cf.n; ++l) {
      if (!pres[zz(l)]) continue;
      if (first) {
        first = false;
        continue;
      }
      for (std::size_t i = 0; i < m; ++i) dm.push_back(cf.c[rw[i]] == l ? 1.0 : 0.0);
      ++nd;
    }
  }
  std::vector<double> z(m, 1.0);
  for (const std::string &v : inst) {
    const std::vector<double> c1 = colr(v);
    z.insert(z.end(), c1.begin(), c1.end());
  }
  for (const std::string &v : x) {
    const std::vector<double> c1 = colr(v);
    z.insert(z.end(), c1.begin(), c1.end());
  }
  z.insert(z.end(), dm.begin(), dm.end());
  const int pz = static_cast<int>(1 + inst.size() + x.size() + nd);
  std::vector<double> x2;
  for (const std::string &v : endo) {
    const std::vector<double> c1 = colr(v);
    x2.insert(x2.end(), c1.begin(), c1.end());
  }
  for (const std::string &v : x) {
    const std::vector<double> c1 = colr(v);
    x2.insert(x2.end(), c1.begin(), c1.end());
  }
  std::vector<double> bz(zz(pz)), xi(zz(pz) * zz(pz)), res(m), wk(zz(pz) * zz(pz) + zz(pz));
  for (const std::string &v : endo) {
    const std::vector<double> ev = colr(v);
    if (pev::ols(z.data(), ev.data(), static_cast<int>(m), pz, nullptr, bz.data(), xi.data(), res.data(), wk.data()) !=
        CLAB_OK) {
      return false;
    }
    x2.insert(x2.end(), res.begin(), res.end());
  }
  x2.insert(x2.end(), m, 1.0);
  x2.insert(x2.end(), dm.begin(), dm.end());
  const int k = static_cast<int>(2 * endo.size() + x.size() + 1 + nd);
  std::vector<double> yr = colr(y), of;
  for (std::size_t i = 0; i < m && !off.empty(); ++i) of.push_back(off[rw[i]]);
  std::vector<int> cr;
  for (std::size_t i = 0; i < m && cl; ++i) cr.push_back(cl[rw[i]]);
  o.k = k;
  o.beta.assign(zz(k), 0.0);
  o.vb.assign(zz(k) * zz(k), 0.0);
  o.vc.assign(zz(k) * zz(k), 0.0);
  o.st.assign(8, 0.0);
  std::vector<double> work(m * zz(k) + 6 * m + zz(k) + zz(k) * zz(k) + 2 * zz(ncl) * zz(k));
  std::vector<int> iw(zz(k)), al(zz(k));
  const int rc = clab_cntreg(yr.data(), x2.data(), static_cast<int>(m), k, of.empty() ? nullptr : of.data(),
                             cl ? cr.data() : nullptr, cl ? ncl : 0, 1, maxit, 1e-10, o.beta.data(), o.vb.data(),
                             o.vc.data(), o.st.data(), al.data(), work.data(), iw.data());
  if (rc != CLAB_OK) return false;
  for (int a : al) {
    if (a) return false;
  }
  return true;
}

bool ppmlcf(Ctx &c, Val &out) {
  static const std::vector<const char *> hows = {"cf", "plug"}, ses = {"boot", "naive"},
                                         fams = {"poisson", "nb"};
  int how = 0, se = 0, ssc = 0, fam = 0;
  CK(marg(c, "how", hows, how));
  CK(marg(c, "se", ses, se));
  CK(marg(c, "ssc", sscs, ssc));
  CK(marg(c, "family", fams, fam));
  if (fam == 1 && how != 0) return ckstop(c, "`family = \"nb\"` needs `how = \"cf\"`");
  Val d;
  CK(dat(c, c["data"], "data", d));
  Cols y, endo, inst, x, fe, expo, clus;
  CK(acol(c, d, "y", {1}, false, y));
  CK(acol(c, d, "endo", {}, false, endo));
  CK(acol(c, d, "inst", {}, false, inst));
  CK(acol(c, d, "x", {}, true, x));
  CK(acol(c, d, "fe", {1, 2}, false, fe));
  CK(acol(c, d, "expo", {1}, true, expo));
  Nm cn = nms(c, {"clus"});
  if (cn.null) cn = nmof(Cols{fe[0]});
  CK(acols(c, d, cn, "clus", {1}, false, clus));
  double lv = 0;
  CK(alevel(c, lv));
  Nv tl, mx;
  CK(cknumd(c, "tol", 1e-12, tl, -Inf, Inf, true));
  CK(cknumd(c, "maxit", 100, mx, 1, Inf, false, true));
  const double tol = tl[0];
  const int maxit = mx[0] > INT_MAX ? INT_MAX : static_cast<int>(mx[0]);
  if (inst.size() < endo.size()) {
    return ckstop(c, "`inst` must give at least as many instruments as `endo` has columns");
  }
  for (const std::string &s : inst) {
    if (has(endo, s) || has(x, s)) return ckstop(c, "a column cannot be both a regressor and an excluded instrument");
  }
  std::vector<double> yv, t;
  CK(amat(c, d, y, yv));
  bool any = false;
  for (double v : yv) {
    if (v < 0) return ckstop(c, "column `%s` must hold counts, at least 0", y[0].c_str());
    any = any || v > 0;
  }
  if (!any) return ckstop(c, "column `%s` is all 0", y[0].c_str());
  std::vector<int> k;
  CK(asep(c, yv, d, fe, k));
  const Val ds = rows(d, which(k));
  const std::size_t n = nr(ds);
  std::vector<double> off;
  if (!expo.empty()) {
    CK(amat(c, ds, expo, t));
    for (double v : t) {
      if (v <= 0) return ckstop(c, "column `%s` must be positive", expo[0].c_str());
    }
    for (double v : t) off.push_back(std::log(v));
  }
  const Code f1 = acode(*fc(ds, fe[0]));
  Code f2;
  if (fe.size() == 2) f2 = acode(*fc(ds, fe[1]));
  Cl cl;
  CK(acl(c, ds, clus, cl));
  Cols tm = endo;
  tm.insert(tm.end(), x.begin(), x.end());
  if (how == 0) {
    for (const std::string &e : endo) tm.push_back("vhat:" + e);
  }
  std::vector<double> me, mxx, mz;
  const int pe = static_cast<int>(endo.size()), q = static_cast<int>(inst.size());
  if (fam == 1) {
    // .axcfnbr(): the first-stage table of ivfe(), then the fit
    int nb = 0;
    bool hs = false;
    double sd = 0;
    if (se == 0) CK(bargs(c, nb, hs, sd));
    const Id a = idat(ds, y[0], endo, x, inst, fe);
    If r1;
    const int rc1 = ivfit(a, all(n), nullptr, cl.a.c.data(), cl.a.n, nullptr, 0, tol, r1);
    if (rc1 != CLAB_OK) return kerr(c, "ivfe", rc1);
    const std::size_t kp = tm.size();
    Nbf f0;
    if (!nbcf(ds, all(n), y[0], endo, inst, x, fe, off, se == 1 ? cl.a.c.data() : nullptr, cl.a.n, maxit, f0)) {
      return ckstop(c, "the negative binomial fit failed");
    }
    std::vector<double> b(f0.beta.begin(), f0.beta.begin() + static_cast<std::ptrdiff_t>(kp)), sev(kp);
    Val draws;
    int nfl = 0;
    if (se == 1) {
      const double g = cl.a.n, dn = static_cast<double>(n), kk = f0.k;
      const double fac = ssc == 0 ? g / (g - 1) : ssc == 1 ? g / (g - 1) * (dn - 1) / (dn - kk) : 1.0;
      for (std::size_t j = 0; j < kp; ++j) sev[j] = std::sqrt(fac * f0.vc[j * zz(f0.k) + j]);
    } else {
      auto fit = [&](const Ix &rw, const std::vector<double> &w, std::vector<double> &beta) {
        Ix rr;
        for (std::size_t i = 0; i < rw.size(); ++i) {
          for (int r = 0; r < static_cast<int>(w[i]); ++r) rr.push_back(rw[i]);
        }
        Nbf fb;
        if (!nbcf(ds, rr, y[0], endo, inst, x, fe, off, nullptr, 0, maxit, fb)) return false;
        beta.assign(fb.beta.begin(), fb.beta.begin() + static_cast<std::ptrdiff_t>(kp));
        return true;
      };
      Bt bt;
      CK(aboot(c, ds, clus[0], nb, hs, sd, static_cast<int>(kp), fit, bt));
      sev = bt.se;
      draws = vmat(bt.dr, nb, static_cast<int>(kp));
      nfl = bt.nf;
    }
    out = vobj();
    put(out, "coef", atab(tm, b, sev, lv, NaN, true));
    put(out, "first", afirst(endo, r1.fst, q));
    put(out, "draws", draws);
    put(out, "nfail", vnum(nfl));
    put(out, "n", vnum(static_cast<double>(n)));
    put(out, "nsep", vnum(static_cast<double>(yv.size() - n)));
    put(out, "ncl", vnum(cl.g[0]));
    put(out, "nboot", vnum(se == 0 ? nb : 0));
    put(out, "how", vstr("cf"));
    put(out, "se", vstr(ses[zz(se)]));
    put(out, "fe", vs(fe));
    put(out, "clus", vs(clus));
    put(out, "level", vnum(lv));
    put(out, "family", vstr("nb"));
    put(out, "theta", vnum(f0.st[2]));
    return true;
  }
  if (se == 1) {
    CK(amat(c, ds, endo, me));
    CK(amat(c, ds, x, mxx));
    CK(amat(c, ds, inst, mz));
    const Id a = idat(ds, y[0], endo, x, inst, fe);
    If r1;
    const int rc = ivfit(a, all(n), nullptr, cl.a.c.data(), cl.a.n, nullptr, 0, tol, r1);
    if (rc != CLAB_OK) return kerr(c, "ivfe", rc);
    Val e = ds;
    for (std::size_t j = 0; j < endo.size(); ++j) {
      std::vector<double> v(n);
      const std::vector<double> ev = colv(*fc(ds, endo[j]));
      for (std::size_t i = 0; i < n; ++i) v[i] = how == 0 ? r1.fres[j * n + i] : ev[i] - r1.fres[j * n + i];
      setc(e, how == 0 ? "vhat:" + endo[j] : endo[j], vnums(v));
    }
    Val ar = vobj(), pf;
    put(ar, "data", e);
    put(ar, "y", vstr(y[0]));
    put(ar, "x", vstrs(tm));
    put(ar, "fe", vstrs(fe));
    if (!expo.empty()) put(ar, "expo", vstr(expo[0]));
    put(ar, "clus", vstrs(clus));
    put(ar, "ssc", vstr(sscs[zz(ssc)]));
    put(ar, "level", vnum(lv));
    put(ar, "tol", vnum(tol));
    put(ar, "maxit", vnum(mx[0]));
    put(ar, "sep", vstr("none"));
    CK(call(c, "ppmlfe", ppmlfe, ar, pf));
    out = vobj();
    put(out, "coef", col(pf, "coef"));
    put(out, "first", afirst(endo, r1.fst, q));
    put(out, "draws", Val());
    put(out, "nfail", vnum(0));
    put(out, "n", vnum(static_cast<double>(n)));
    put(out, "nsep", vnum(static_cast<double>(yv.size() - n)));
    put(out, "ncl", vnum(cl.g[0]));
    put(out, "nboot", vnum(0));
  } else {
    int nb = 0;
    bool hs = false;
    double sd = 0;
    CK(bargs(c, nb, hs, sd));
    std::uint32_t s0 = 0;
    CK(adseed(c, hs, sd, s0));
    Adraw ad(s0, cl.a.n);
    std::vector<int> bm;
    CK(amat(c, ds, endo, me));
    CK(amat(c, ds, x, mxx));
    CK(amat(c, ds, inst, mz));
    const std::size_t kp = tm.size();
    // fst is pe x 4, as src/bapx.cpp gives it
    std::vector<double> beta(kp), fst(zz(pe) * 4, NA), draw(zz(nb) * kp), st(3), db;
    const std::vector<double> ys = colv(*fc(ds, y[0]));
    // the draws a block at a time, as R makes them: the full-sample fit
    // is the same in every block
    int nf = 0;
    for (std::size_t b0 = 0; b0 < zz(nb); b0 += ad.blk) {
      const std::size_t kb = std::min(ad.blk, zz(nb) - b0);
      ad.next(kb, bm);
      db.assign(kb * kp, NA);
      const int rc = clab_ppmlcf(pd(ys), pd(me), pe, pd(mxx), static_cast<int>(x.size()), pd(mz), q,
                                 static_cast<int>(n), off.empty() ? nullptr : off.data(), f1.c.data(), f1.n,
                                 fe.size() == 2 ? f2.c.data() : nullptr, f2.n, cl.a.c.data(), cl.a.n, bm.data(),
                                 static_cast<int>(kb), how, maxit, 1e-10, tol, beta.data(), fst.data(), db.data(),
                                 st.data());
      if (rc != CLAB_OK) return kerr(c, "ppmlcf", rc);
      for (std::size_t j = 0; j < kp; ++j) {
        for (std::size_t e = 0; e < kb; ++e) draw[j * zz(nb) + b0 + e] = db[j * kb + e];
      }
      nf += static_cast<int>(st[2]);
    }
    if (nf > 0) {
      wrn(c, std::to_string(nf) + " of " + std::to_string(nb) +
                 " bootstrap draws could not be fitted and are left out.");
    }
    if (nb - nf < 2) return ckstop(c, "fewer than two bootstrap draws could be fitted");
    const std::vector<double> sev = csd(draw, zz(nb), static_cast<int>(kp));
    out = vobj();
    put(out, "coef", atab(tm, beta, sev, lv, NaN, true));
    put(out, "first", afirst(endo, fst, q));
    put(out, "draws", vmat(draw, nb, static_cast<int>(kp)));
    put(out, "nfail", vnum(nf));
    put(out, "n", vnum(static_cast<double>(n)));
    put(out, "nsep", vnum(static_cast<double>(yv.size() - n)));
    put(out, "ncl", vnum(cl.g[0]));
    put(out, "nboot", vnum(nb));
  }
  put(out, "how", vstr(hows[zz(how)]));
  put(out, "se", vstr(ses[zz(se)]));
  put(out, "fe", vs(fe));
  put(out, "clus", vs(clus));
  put(out, "level", vnum(lv));
  return true;
}

// cut(): the labels R makes for breaks b (right-closed; the first closed
// on the left too with low), formatC(b, digits = d, width = 1) for the
// first d from 3 that tells the breaks apart; and the bin of x (from 0,
// -1 for none), as .bincode().
inline std::string fmtg(double x, int dig) {
  if (std::isnan(x)) return isna(x) ? "NA" : "NaN";
  if (std::isinf(x)) return x > 0 ? "Inf" : "-Inf";
  char b[64];
  std::snprintf(b, sizeof b, "%1.*g", dig, 0 + x);
  return b;
}
inline Cols cutlab(const std::vector<double> &b, bool low) {
  Cols ch;
  for (int dig = 3; dig <= 12; ++dig) {
    ch.clear();
    for (double x : b) ch.push_back(fmtg(x, dig));
    bool ok = true;
    for (std::size_t i = 0; i + 1 < ch.size(); ++i) ok = ok && ch[i] != ch[i + 1];
    if (ok) break;
  }
  Cols l;
  for (std::size_t i = 0; i + 1 < ch.size(); ++i) {
    l.push_back((i == 0 && low ? "[" : "(") + ch[i] + "," + ch[i + 1] + "]");
  }
  return l;
}
inline int bincode(double x, const std::vector<double> &b, bool low) {
  if (std::isnan(x)) return -1;
  std::size_t lo = 0, hi = b.size() - 1;
  if (x < b[lo] || b[hi] < x || (x == b[lo] && !low)) return -1;
  while (hi - lo >= 2) {
    const std::size_t nw = (hi + lo) / 2;
    if (x > b[nw]) lo = nw;
    else hi = nw;
  }
  return static_cast<int>(lo);
}
// is.unsorted(x, strictly = TRUE) is false
inline bool incr(const std::vector<double> &x) {
  for (std::size_t i = 0; i + 1 < x.size(); ++i) {
    if (!(x[i] < x[i + 1])) return false;
  }
  return true;
}
inline bool anyna(const std::vector<double> &x) {
  for (double v : x) {
    if (std::isnan(v)) return true;
  }
  return false;
}

// clab_windbin: n x na (nb - 1), column-major.
inline int wbin(Ctx &c, const std::vector<double> &dir, double width, double start, const std::vector<int> *ac, int na,
                int base, int edge, std::vector<double> &m, int &nc) {
  const double nbd = std::nearbyint(360 / width);
  if (nbd > kutil::maxcell) return (c.fail("`width` must divide 360 into at most %.0f sectors.", kutil::maxcell), CLAB_EARG);
  if (static_cast<double>(na) * (nbd - 1) > 2147483647.0) return (c.fail("the areas times the sectors of `width` must be at most 2147483647 columns."), CLAB_EARG);
  const int nb = static_cast<int>(nbd);
  nc = na * (nb - 1);
  m.assign(dir.size() * zz(nc), 0.0);
  const int rc = clab_windbin(pd(dir), static_cast<int>(dir.size()), width, start, ac ? ac->data() : nullptr, na,
                              base, edge, m.empty() ? nullptr : m.data());
  return rc == CLAB_OK ? rc : (kerr(c, "windiv", rc), rc);
}
// .axsec(): the sector (1..nb) of each direction.
inline double asec(double dir, double width, double start, bool paper) {
  const double a = rmod(dir - start, 360);
  if (std::isnan(a)) return NA;
  const double b = paper ? (a == 0 ? 0 : std::ceil(a / width) - 1) : std::floor(a / width);
  return std::min(b, std::nearbyint(360 / width) - 1) + 1;
}

bool windiv(Ctx &c, Val &out) {
  int edges = 0, avg = 0;
  CK(marg(c, "edges", {"halfopen", "paper"}, edges));
  CK(marg(c, "avg", {"vector", "mode", "mean"}, avg));
  const Val &dv = c["dir"];
  if (!numv(dv) || vlen(dv) == 0) return ckstop(c, "`dir` must be numeric directions in degrees");
  std::vector<double> dir = nums(dv);
  const std::size_t n = dir.size();
  bool share = false;
  if (c.has("share") && !tf(c["share"], share)) return ckstop(c, "`share` must be TRUE or FALSE");
  std::vector<double> spd;
  if (c.has("spd")) {
    spd = nums(c["spd"]);
    bool ok = numv(c["spd"]) && spd.size() == n;
    for (double v : spd) ok = ok && !(v < 0);
    if (!ok) return ckstop(c, "`spd` must give one speed from 0 per direction");
  }
  Cols area;
  if (c.has("area")) {
    const Val ar = asarr(c["area"]);
    bool ok = ar.a.size() == n;
    for (const Val &e : ar.a) ok = ok && !vna(e) && e.t != Val::ARR && e.t != Val::OBJ;
    if (!ok) return ckstop(c, "`area` must give one label per direction, none missing");
    area = coll(ar);
  }
  const bool sect = c.has("sect");
  double width = 90, start = 0, s1 = 0, s2 = 0;
  int nb = 0, base = 1;
  if (!sect) {
    Nv q;
    CK(cknumd(c, "width", 90, q, -Inf, Inf, true));
    width = q[0];
    CK(cknumd(c, "start", 0, q));
    start = q[0];
    nb = static_cast<int>(std::nearbyint(360 / width));
    if (width >= 360 || std::fabs(360 / width - nb) > 1e-9) {
      return ckstop(c, "`width` must divide 360 into two or more sectors");
    }
    CK(cknumd(c, "base", 1, q, 1, nb, false, true));
    base = static_cast<int>(q[0]);
  } else {
    const std::vector<double> s = nums(c["sect"]);
    bool ok = numv(c["sect"]) && s.size() == 2;
    for (double v : s) ok = ok && std::isfinite(v);
    if (!ok || rmod(s[0], 360) == rmod(s[1], 360)) {
      return ckstop(c, "`sect` must be two different finite directions in degrees");
    }
    s1 = s[0];
    s2 = s[1];
  }
  const bool paper = edges == 1;
  const double pi = 3.141592653589793238462643383279502884;
  if (c.has("by")) {
    const Val bv = asarr(c["by"]);
    bool ok = bv.a.size() == n;
    for (const Val &e : bv.a) ok = ok && !vna(e) && e.t != Val::ARR && e.t != Val::OBJ;
    if (!ok) return ckstop(c, "`by` must give one group per direction, none missing");
    const Code gk = acode(bv);
    const std::size_t ng = zz(gk.n);
    std::vector<char> ok2(n);
    for (std::size_t i = 0; i < n; ++i) ok2[i] = std::isfinite(dir[i]) && (spd.empty() || std::isfinite(spd[i]));
    // the groups' means of the values of the rows ok
    auto gmean = [&](const std::vector<double> &v) {
      std::vector<std::vector<double>> s(ng);
      for (std::size_t i = 0; i < n; ++i) {
        if (ok2[i]) s[zz(gk.c[i])].push_back(v[i]);
      }
      std::vector<double> m;
      for (const auto &g : s) m.push_back(g.empty() ? NaN : lmean(g));
      return m;
    };
    std::vector<double> sn(n), cs(n);
    for (std::size_t i = 0; i < n; ++i) {
      const double rd = dir[i] * pi / 180, wv = spd.empty() ? 1.0 : spd[i];
      sn[i] = std::sin(rd) * wv;
      cs[i] = std::cos(rd) * wv;
    }
    const std::vector<double> ms = gmean(sn), mc = gmean(cs);
    std::vector<double> vdir(ng);
    for (std::size_t g = 0; g < ng; ++g) vdir[g] = rmod(std::atan2(ms[g], mc[g]) * 180 / pi, 360);
    if (!area.empty()) {
      Cols ga(ng);
      std::vector<char> seen(ng, 0);
      for (std::size_t i = 0; i < n; ++i) {
        const std::size_t g = zz(gk.c[i]);
        if (seen[g] && ga[g] != area[i]) return ckstop(c, "`area` must be constant within each `by` group");
        if (!seen[g]) ga[g] = area[i];
        seen[g] = 1;
      }
      area = ga;
    }
    if (share) {
      if (sect) return ckstop(c, "`share` needs sectors by `width`, not `sect`");
      bool anyok = false;
      std::vector<double> d2(n);
      for (std::size_t i = 0; i < n; ++i) {
        d2[i] = ok2[i] ? dir[i] : NA;
        anyok = anyok || ok2[i];
      }
      if (!anyok) return ckstop(c, "`dir` must be numeric directions in degrees");
      std::vector<double> full;
      int k = 0;
      if (wbin(c, d2, width, start, nullptr, 1, base - 1, paper ? 1 : 0, full, k) != CLAB_OK) return false;
      std::vector<double> cnt(ng, 0.0), m(ng * zz(k), 0.0);
      for (std::size_t i = 0; i < n; ++i) cnt[zz(gk.c[i])] += ok2[i] ? 1 : 0;
      for (std::size_t j = 0; j < zz(k); ++j) {
        for (std::size_t i = 0; i < n; ++i) {
          const double v = full[j * n + i];
          m[j * ng + zz(gk.c[i])] += std::isnan(v) ? 0.0 : v;
        }
        for (std::size_t g = 0; g < ng; ++g) m[j * ng + g] /= cnt[g];
      }
      if (area.empty()) {
        out = vmat(m, static_cast<int>(ng), k);
        return true;
      }
      Cols lev = area;
      std::sort(lev.begin(), lev.end());
      lev.erase(std::unique(lev.begin(), lev.end()), lev.end());
      const std::size_t nl = lev.size();
      std::vector<double> mm(ng * zz(k) * nl, 0.0);
      for (std::size_t g = 0; g < ng; ++g) {
        const std::size_t ai = zz(static_cast<int>(std::lower_bound(lev.begin(), lev.end(), area[g]) - lev.begin()));
        for (std::size_t j = 0; j < zz(k); ++j) mm[(ai * zz(k) + j) * ng + g] = m[j * ng + g];
      }
      out = vmat(mm, static_cast<int>(ng), static_cast<int>(zz(k) * nl));
      return true;
    }
    if (avg == 0) {
      dir = vdir;
    } else if (avg == 2) {
      dir = gmean(dir);
    } else {
      if (sect) return ckstop(c, "`avg = \"mode\"` needs sectors by `width`, not `sect`");
      std::vector<std::vector<int>> tab(ng, std::vector<int>(zz(nb), 0));
      std::vector<char> hit(ng, 0);
      for (std::size_t i = 0; i < n; ++i) {
        if (!ok2[i]) continue;
        const double s = asec(dir[i], width, start, paper);
        if (std::isnan(s)) continue;
        ++tab[zz(gk.c[i])][zz(static_cast<int>(s) - 1)];
        hit[zz(gk.c[i])] = 1;
      }
      dir.assign(ng, NA);
      for (std::size_t g = 0; g < ng; ++g) {
        if (!hit[g]) continue;
        const int md = static_cast<int>(std::max_element(tab[g].begin(), tab[g].end()) - tab[g].begin()) + 1;
        dir[g] = rmod(start + (md - 0.5) * width, 360);
      }
    }
  }
  const std::size_t m = dir.size();
  if (sect) {
    std::vector<double> v(m);
    for (std::size_t i = 0; i < m; ++i) {
      const double a = rmod(dir[i] - s1, 360), b = rmod(s2 - s1, 360);
      v[i] = std::isnan(a) || std::isnan(b) ? NA : (a < b ? 1.0 : 0.0);
    }
    out = vmat(v, static_cast<int>(m), 1);
    return true;
  }
  std::vector<int> ac;
  int na = 1;
  if (!area.empty()) {
    Cols lev = area;
    std::sort(lev.begin(), lev.end());
    lev.erase(std::unique(lev.begin(), lev.end()), lev.end());
    for (const std::string &s : area) {
      ac.push_back(static_cast<int>(std::lower_bound(lev.begin(), lev.end(), s) - lev.begin()));
    }
    na = static_cast<int>(lev.size());
  }
  std::vector<double> w;
  int k = 0;
  if (wbin(c, dir, width, start, area.empty() ? nullptr : &ac, na, base - 1, paper ? 1 : 0, w, k) != CLAB_OK) {
    return false;
  }
  out = vmat(w, static_cast<int>(m), k);
  return true;
}

// invsn() by group (R's .axinby()): per sorted group the readings with
// both temperatures, their mean difference by sign and the inversion by
// rule (0 all, 1 any, 2 the mean above 0).
bool invby(Ctx &c, const std::vector<double> &up, int sg, int rule, Val &out) {
  const std::size_t n = up.size();
  const Val by = asarr(c["by"]);
  bool num = true, na = by.a.size() != n;
  for (const Val &e : by.a) {
    na = na || vna(e);
    num = num && e.t == Val::NUM;
  }
  if (na) return ckstop(c, "`by` must give the group of each reading, with no NA");
  std::vector<double> gv;
  Cols gs;
  for (const Val &e : by.a) {
    if (num) gv.push_back(nv(e));
    else gs.push_back(vlab(e));
  }
  std::vector<double> lv = gv;
  Cols ls = gs;
  std::sort(lv.begin(), lv.end());
  lv.erase(std::unique(lv.begin(), lv.end()), lv.end());
  std::sort(ls.begin(), ls.end());
  ls.erase(std::unique(ls.begin(), ls.end()), ls.end());
  const std::size_t ng = num ? lv.size() : ls.size();
  std::vector<std::vector<double>> u(ng);
  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t k =
        num ? static_cast<std::size_t>(std::lower_bound(lv.begin(), lv.end(), gv[i]) - lv.begin())
            : static_cast<std::size_t>(std::lower_bound(ls.begin(), ls.end(), gs[i]) - ls.begin());
    if (std::isfinite(up[i])) u[k].push_back(up[i]);
  }
  std::vector<double> nn(ng), df(ng);
  std::vector<int> inv(ng);
  for (std::size_t k = 0; k < ng; ++k) {
    nn[k] = static_cast<double>(u[k].size());
    if (u[k].empty()) {
      df[k] = NAN;
      inv[k] = INT_MIN;
      continue;
    }
    const double m = lmean(u[k]);
    bool all = true, any = false;
    for (double x : u[k]) {
      all = all && x > 0;
      any = any || x > 0;
    }
    df[k] = sg == 0 ? m : -m;
    inv[k] = rule == 0 ? all : (rule == 1 ? any : m > 0);
  }
  Frame f;
  if (num) f.num("by", lv);
  else f.str("by", ls);
  f.num("n", nn);
  Val dv = varr();
  for (double x : df) dv.a.push_back(std::isnan(x) ? Val() : vnum(x));
  f.add("diff", dv);
  f.ints("inv", inv.data(), ng);
  if (c.has("brk")) {
    const std::vector<double> br = nums(c["brk"]);
    if (!numv(c["brk"]) || br.size() < 2 || anyna(br) || !incr(br)) {
      return ckstop(c, "`brk` must be two or more increasing numbers");
    }
    const Cols l = cutlab(br, false);
    Val bin = varr();
    for (double x : df) {
      const int k = std::isnan(x) ? -1 : bincode(x, br, false);
      bin.a.push_back(k < 0 ? Val() : vstr(l[zz(k)]));
    }
    f.add("bin", bin);
  }
  out = f.o;
  return true;
}

bool invsn(Ctx &c, Val &out) {
  int sg = 0, rule = 0;
  CK(marg(c, "sign", {"upper", "ground"}, sg));
  CK(marg(c, "rule", {"all", "any", "mean"}, rule));
  const Val &lo = c["tlow"], &hi = c["tup"];
  const std::vector<double> a = nums(lo), b = nums(hi);
  if (!numv(lo) || !numv(hi) || a.size() != b.size() || a.empty()) {
    return ckstop(c, "`tlow` and `tup` must be numeric vectors of one length");
  }
  const std::size_t n = a.size();
  if (c.has("by")) {
    std::vector<double> up(n);
    for (std::size_t i = 0; i < n; ++i) up[i] = b[i] - a[i];
    return invby(c, up, sg, rule, out);
  }
  std::vector<double> df(n);
  std::vector<int> inv(n);
  for (std::size_t i = 0; i < n; ++i) {
    const double up = b[i] - a[i];
    df[i] = sg == 0 ? up : -up;
    inv[i] = std::isnan(up) ? INT_MIN : (up > 0 ? 1 : 0);
  }
  Frame f;
  f.num("diff", df);
  f.ints("inv", inv.data(), n);
  if (c.has("brk")) {
    const std::vector<double> br = nums(c["brk"]);
    if (!numv(c["brk"]) || br.size() < 2 || anyna(br) || !incr(br)) {
      return ckstop(c, "`brk` must be two or more increasing numbers");
    }
    const Cols l = cutlab(br, false);
    Val bin = varr();
    for (double x : df) {
      const int k = bincode(x, br, false);
      bin.a.push_back(k < 0 ? Val() : vstr(l[zz(k)]));
    }
    f.add("bin", bin);
  }
  out = f.o;
  return true;
}

bool expjn(Ctx &c, Val &out) {
  Val e, o;
  CK(dat(c, c["expo"], "expo", e));
  CK(dat(c, c["outc"], "outc", o));
  Cols by = {"area", "day"};
  if (c.has("by")) {
    const Nm m = nms(c, {"by"});
    bool dup = false;
    for (std::size_t i = 0; i < m.v.size(); ++i) {
      for (std::size_t j = 0; j < i; ++j) dup = dup || m.v[i] == m.v[j];
    }
    if (!m.chr || m.v.empty() || m.na || dup) return ckstop(c, "`by` must name the key columns");
    by = m.v;
  }
  for (int t = 0; t < 2; ++t) {
    const Val &d = t == 0 ? e : o;
    const char *nm = t == 0 ? "expo" : "outc";
    Cols miss;
    for (const std::string &s : by) {
      if (fc(d, s) == nullptr) miss.push_back(s);
    }
    if (!miss.empty()) return ckstop(c, "`%s` lacks the key columns: %s", nm, join(miss, ", ").c_str());
    for (const std::string &s : by) {
      for (const Val &v : fc(d, s)->a) {
        if (vna(v)) return ckstop(c, "`%s` has missing keys", nm);
      }
    }
  }
  Cols val, clash;
  for (const std::string &s : e.k) {
    if (!has(by, s) && !has(val, s)) val.push_back(s);
  }
  if (val.empty()) return ckstop(c, "`expo` has no columns besides the keys");
  for (const std::string &s : val) {
    if (fc(o, s) != nullptr) clash.push_back(s);
  }
  if (!clash.empty()) return ckstop(c, "columns in both tables besides the keys: %s", join(clash, ", ").c_str());
  auto keys = [&](const Val &d) {
    Cols k(nr(d));
    for (std::size_t j = 0; j < by.size(); ++j) {
      const Val &a = *fc(d, by[j]);
      for (std::size_t i = 0; i < k.size(); ++i) k[i] += (j ? "\r" : "") + vlab(a.a[i]);
    }
    return k;
  };
  const Cols ke = keys(e), ko = keys(o);
  std::unordered_map<std::string, std::size_t> he, ho;
  int de = 0, dq = 0;
  for (std::size_t i = 0; i < ke.size(); ++i) de += he.emplace(ke[i], i).second ? 0 : 1;
  for (std::size_t i = 0; i < ko.size(); ++i) dq += ho.emplace(ko[i], i).second ? 0 : 1;
  if (de > 0) return ckstop(c, "`expo` has %d duplicate key(s)", de);
  if (dq > 0) return ckstop(c, "`outc` has %d duplicate key(s)", dq);
  Ix hit, mh, miss;
  for (std::size_t i = 0; i < ko.size(); ++i) {
    const auto it = he.find(ko[i]);
    if (it == he.end()) {
      miss.push_back(i);
    } else {
      hit.push_back(i);
      mh.push_back(it->second);
    }
  }
  if (hit.empty()) return ckstop(c, "no outcome row has an exposure: check that the keys agree");
  Val dd = rows(o, hit);
  const Val em = rows(e, mh);
  for (const std::string &s : val) put(dd, s.c_str(), *fc(em, s));
  Val un = vobj();
  const Val om = rows(o, miss);
  for (const std::string &s : by) put(un, s.c_str(), *fc(om, s));
  int nexp = 0;
  for (const std::string &s : ke) nexp += ho.count(s) ? 0 : 1;
  out = vobj();
  put(out, "data", dd);
  put(out, "nout", vnum(static_cast<double>(miss.size())));
  put(out, "nexp", vnum(nexp));
  put(out, "unout", un);
  return true;
}

bool smkiv(Ctx &c, Val &out) {
  static const std::vector<const char *> adjs = {"raw", "own", "nbr"};
  int adj = 0, cf = 0;
  CK(marg(c, "adj", adjs, adj));
  CK(marg(c, "cf", {"mean", "median"}, cf));
  Val d;
  CK(dat(c, c["data"], "data", d));
  const std::size_t n = nr(d);
  Cols pud, pm, seas, by;
  CK(acols(c, d, nms(c, {"plume", "unit", "day"}), "plume, unit, day", {3}, false, pud));
  CK(acol(c, d, "pm", {1}, true, pm));
  CK(acol(c, d, "seas", {1}, true, seas));
  CK(acol(c, d, "by", {1}, true, by));
  std::string dens;
  if (c.has("dens")) {
    const Val &o = c["dens"];
    const Val *s = o.t == Val::ARR && o.a.size() == 1 ? &o.a[0] : &o;
    if (s->t != Val::STR || fc(d, s->s) == nullptr) return ckstop(c, "`dens` must name a column of `data`");
    dens = s->s;
  }
  std::vector<double> pl;
  CK(amat(c, d, Cols{pud[0]}, pl));
  for (double v : pl) {
    if (v != 0 && v != 1) return ckstop(c, "column `%s` must hold 0 or 1", pud[0].c_str());
  }
  Nv zq;
  CK(cknumd(c, "z", 1.64, zq, 0));
  const double z = zq[0];
  if (adj != 0 && pm.empty()) return ckstop(c, "`adj = \"%s\"` needs `pm`", adjs[zz(adj)]);
  Cols sv;
  if (!seas.empty()) {
    sv = coll(*fc(d, seas[0]));
  } else {
    // a day column outside R is numbers, never dates
    if (!pm.empty()) return ckstop(c, "give `seas`, or `day` as dates for the meteorological seasons");
    sv.assign(n, "all");
  }
  const Cols ul = coll(*fc(d, pud[1])), dl = coll(*fc(d, pud[2]));
  Cols key(n);
  for (std::size_t i = 0; i < n; ++i) key[i] = ul[i] + "\r" + sv[i];
  Cols nm = {"smoke"};
  std::vector<std::vector<double>> cols = {pl};
  if (!pm.empty()) {
    std::vector<double> pv;
    CK(amat(c, d, pm, pv));
    std::unordered_map<std::string, std::vector<double>> bg;
    std::unordered_map<std::string, int> nk;
    for (std::size_t i = 0; i < n; ++i) {
      int &m = nk[key[i]];
      if (pl[i] == 0) {
        ++m;
        bg[key[i]].push_back(pv[i]);
      }
    }
    for (const auto &kv : nk) {
      if (kv.second < 2) {
        return ckstop(c, "a place and season has fewer than two days without a plume for its background");
      }
    }
    std::unordered_map<std::string, std::array<double, 3>> st;
    for (auto &kv : bg) {
      std::vector<double> &v = kv.second;
      std::vector<char> ok(v.size(), 1);
      const double mu = lmean(v), sd = std::sqrt(rcov(v, v.size(), 1, ok)[0]);
      std::vector<double> s = v;
      std::sort(s.begin(), s.end());
      const std::size_t h = (s.size() + 1) / 2;
      const double md = s.size() % 2 ? s[h - 1] : lmean({s[h - 1], s[h]});
      st[kv.first] = {mu, sd, md};
    }
    std::vector<double> hi(n), sm(n), bgv(n), spm(n);
    for (std::size_t i = 0; i < n; ++i) {
      const auto &s = st[key[i]];
      hi[i] = pl[i] == 1 && pv[i] > s[0] + z * s[1] ? 1 : 0;
    }
    sm = pl;
    if (adj == 1) sm = hi;
    if (adj == 2) {
      if (!c.has("nbr")) return ckstop(c, "`adj = \"nbr\"` needs `nbr`, the neighbour of each row");
      Cols nb;
      CK(acol(c, d, "nbr", {1}, false, nb));
      const Cols nl = coll(*fc(d, nb[0]));
      std::unordered_map<std::string, std::size_t> h;
      for (std::size_t i = 0; i < n; ++i) h.emplace(ul[i] + "\r" + dl[i], i);
      for (std::size_t i = 0; i < n; ++i) {
        const auto it = h.find(nl[i] + "\r" + dl[i]);
        sm[i] = it == h.end() ? NA : hi[it->second];
      }
    }
    for (std::size_t i = 0; i < n; ++i) {
      const auto &s = st[key[i]];
      bgv[i] = cf == 0 ? s[0] : s[2];
      const double v = pv[i] - bgv[i];
      spm[i] = sm[i] == 1 ? (0 > v ? 0.0 : v) : 0.0;
    }
    cols[0] = sm;
    nm.push_back("bg");
    cols.push_back(bgv);
    nm.push_back("spm");
    cols.push_back(spm);
  }
  if (!dens.empty()) {
    const Val &a = *fc(d, dens);
    Cols lv;
    for (const Val &e : a.a) {
      if (!vna(e) && vlab(e) != "") lv.push_back(vlab(e));
    }
    std::sort(lv.begin(), lv.end());
    lv.erase(std::unique(lv.begin(), lv.end()), lv.end());
    for (const std::string &l : lv) {
      std::vector<double> v(n);
      for (std::size_t i = 0; i < n; ++i) v[i] = !vna(a.a[i]) && vlab(a.a[i]) == l ? 1.0 : 0.0;
      nm.push_back("dens:" + l);
      cols.push_back(v);
    }
  }
  Frame f;
  if (!by.empty()) {
    const Code g = acode(*fc(d, by[0]));
    f.str(by[0].c_str(), g.lev);
    for (std::size_t j = 0; j < nm.size(); ++j) {
      std::vector<std::vector<double>> s(zz(g.n));
      for (std::size_t i = 0; i < n; ++i) s[zz(g.c[i])].push_back(cols[j][i]);
      std::vector<double> v;
      for (const auto &x : s) v.push_back(nm[j] == "bg" || nm[j] == "spm" ? lmean(x) : lsum(x));
      f.num(nm[j].c_str(), v);
    }
  } else {
    for (std::size_t j = 0; j < nm.size(); ++j) f.num(nm[j].c_str(), cols[j]);
  }
  out = f.o;
  return true;
}

bool exprun(Ctx &c, Val &out) {
  Val d;
  CK(dat(c, c["data"], "data", d));
  Cols vud;
  CK(acols(c, d, nms(c, {"var", "unit", "day"}), "var, unit, day", {3}, false, vud));
  Nv q;
  CK(cknumd(c, "thr", 0, q));
  const double thr = q[0];
  CK(cknumd(c, "cap", 3, q, 1, Inf, false, true));
  const double cap = q[0];
  if (cap > kutil::maxcell) return c.fail("`cap` must be at most %.0f.", kutil::maxcell);
  CK(cknumd(c, "step", 1, q, -Inf, Inf, true));
  const double step = q[0];
  const Val &tc = *fc(d, vud[2]);
  if (!isnumc(tc)) return ckstop(c, "column `%s` must hold numbers or dates", vud[2].c_str());
  std::vector<double> e;
  CK(amat(c, d, Cols{vud[0]}, e));
  for (double &v : e) v = v > thr ? 1.0 : 0.0;
  const Code u = acode(*fc(d, vud[1]));
  const std::vector<double> tv = colv(tc);
  const std::size_t n = e.size();
  std::vector<int> run(n), ord(n);
  const int rc = clab_runlen(u.c.data(), tv.data(), e.data(), static_cast<int>(n), step, run.data(), ord.data());
  if (rc == CLAB_EDOM) {
    c.err = "exprun: a unit has two rows on one day: give one row per unit and day.";
    return false;
  }
  if (rc != CLAB_OK) {
    c.err = std::string("exprun: ") + clab_strerror(rc);
    return false;
  }
  for (int &v : run) {
    if (v < 0) v = INT_MIN;
  }
  Frame f;
  f.ints("run", run.data(), n);
  for (double k = 1; k <= cap; ++k) {
    const std::string nm = "run" + rfmt(k) + (k < cap ? "" : "p");
    std::vector<double> v(n);
    for (std::size_t i = 0; i < n; ++i) {
      v[i] = run[i] == INT_MIN ? NA : ((k < cap ? run[i] == k : run[i] >= k) ? 1.0 : 0.0);
    }
    f.num(nm.c_str(), v);
  }
  out = f.o;
  return true;
}

bool updwn(Ctx &c, Val &out) {
  static const std::vector<const char *> norms = {"mean", "none", "log", "pois"},
                                         ses = {"hc1", "cluster", "twoway", "nw", "boot"};
  int twin = 0, norm = 0, se = 0, ssc = 0;
  CK(marg(c, "twin", {"drop", "near", "keep"}, twin));
  CK(marg(c, "norm", norms, norm));
  CK(marg(c, "se", ses, se));
  CK(marg(c, "ssc", sscs, ssc));
  Val d, s;
  CK(dat(c, c["data"], "data", d));
  const std::size_t n = nr(d);
  Cols y, axy, day, wind, area, clus, wx;
  CK(acol(c, d, "y", {1}, false, y));
  CK(acols(c, d, nms(c, {"ax", "ay"}), "ax, ay", {2}, false, axy));
  CK(acol(c, d, "day", {1}, false, day));
  CK(acol(c, d, "wind", {1}, false, wind));
  CK(acol(c, d, "area", {1}, true, area));
  CK(acol(c, d, "clus", {1, 2}, true, clus));
  CK(acol(c, d, "wx", {}, true, wx));
  CK(dat(c, c["src"], "src", s));
  if (fc(s, "x0") == nullptr || fc(s, "y0") == nullptr) return ckstop(c, "`src` must have columns x0 and y0");
  const std::size_t ns = nr(s);
  std::vector<std::vector<double>> sv;
  for (const char *k : {"x0", "y0", "x1", "y1"}) {
    const Val *a = fc(s, k);
    if (a != nullptr && !isnumc(*a)) {
      bool alln = true;
      for (const Val &e : a->a) alln = alln && vna(e);
      if (!alln) return ckstop(c, "`src` column %s must be numeric", k);
    }
    sv.push_back(a != nullptr ? colv(*a) : std::vector<double>(ns, NA));
  }
  const std::vector<double> &s0 = sv[0], &t0 = sv[1], &s1 = sv[2], &t1 = sv[3];
  for (std::size_t j = 0; j < ns; ++j) {
    if (!std::isfinite(s0[j]) || !std::isfinite(t0[j])) return ckstop(c, "`src` x0 and y0 must be finite");
  }
  std::vector<char> line(ns);
  bool allseg = true;
  for (std::size_t j = 0; j < ns; ++j) {
    line[j] = std::isfinite(s1[j]) && std::isfinite(t1[j]);
    allseg = allseg && line[j];
  }
  Nv q;
  CK(cknumd(c, "tol", 60, q, 0, 90, true));
  const double tol = q[0];
  double band = Inf;
  if (c.has("band")) {
    Dv bd;
    const bool inf = readd(c["band"], bd) == NUM && bd.n() == 1 && bd[0] == Inf && c["band"].t != Val::OBJ;
    if (!inf) {
      CK(cknum(c, c["band"], "band", q, true, -Inf, Inf, true));
      band = q[0];
    }
  }
  std::vector<double> dist;
  if (c.has("dist")) {
    dist = nums(c["dist"]);
    if (!numv(c["dist"]) || dist.size() < 2 || anyna(dist) || dist[0] < 0 || !incr(dist)) {
      return ckstop(c, "`dist` must be two or more increasing distances from 0");
    }
  }
  bool agg = true;
  if (c.has("agg") && !tf(c["agg"], agg)) return ckstop(c, "`agg` must be TRUE or FALSE");
  if (area.empty() && (!agg || !allseg)) return ckstop(c, "`area` is needed for point sources and for agg = FALSE");
  std::vector<double> offs;
  if (c.has("offs")) {
    if (!allseg) return ckstop(c, "`offs` moves segments: every source must be a segment");
    offs = nums(c["offs"]);
    bool ok = numv(c["offs"]) && !offs.empty();
    for (double v : offs) ok = ok && std::isfinite(v);
    if (!ok) return ckstop(c, "`offs` must be finite numbers");
  }
  if (se == 2 && clus.size() != 2) return ckstop(c, "`se = \"twoway\"` needs two `clus` columns");
  double lv = 0;
  CK(alevel(c, lv));
  std::vector<double> yv, xy, wd, t;
  CK(amat(c, d, y, yv));
  for (double v : yv) {
    if (v < 0) return ckstop(c, "column `%s` must hold counts, at least 0", y[0].c_str());
  }
  CK(amat(c, d, axy, xy));
  CK(amat(c, d, wind, wd));
  if (!wx.empty()) CK(amat(c, d, wx, t));
  const Val &dc = *fc(d, day[0]);
  const Cols dv = coll(dc);
  const bool tnum = isnumc(dc);
  if (se == 3 && !tnum) return ckstop(c, "`se = \"nw\"` needs a day column of numbers or dates");
  Cols av;
  if (!area.empty()) {
    av = coll(*fc(d, area[0]));
    std::unordered_map<std::string, int> h;
    for (std::size_t i = 0; i < n; ++i) {
      if (!h.emplace(av[i] + "\r" + dv[i], 0).second) return ckstop(c, "an area has two rows on one day");
    }
  }
  static const char *const labs[] = {"downwind", "upwind", "along", "nosource", "twosources", "onsource", "missing"};
  // labs by dn: 1 downwind, 0 upwind, -1 along, ... -5 missing
  auto li = [](int dn) { return zz(dn == 1 ? 0 : 1 - dn); };
  const std::vector<double> axv(xy.begin(), xy.begin() + static_cast<std::ptrdiff_t>(n)),
      ayv(xy.begin() + static_cast<std::ptrdiff_t>(n), xy.end());
  struct K {
    std::vector<int> src, side, dn, row;
    std::vector<double> dst;
  };
  struct Run {
    Val fit, tab, coef;
    double n = 0;
  };
  // the kernel on sources [j0, j1), appended to k
  auto kern = [&](const std::vector<double> &a0, const std::vector<double> &b0, const std::vector<double> &a1,
                  const std::vector<double> &b1, std::size_t j0, std::size_t j1, int tw, K &k) {
    std::vector<int> sr(n), sd(n), dn(n);
    std::vector<double> ds(n);
    const int rc = clab_updwn(pd(axv), pd(ayv), pd(wd), static_cast<int>(n), &a0[j0], &b0[j0], &a1[j0], &b1[j0],
                              static_cast<int>(j1 - j0), band, tol, tw, sr.data(), sd.data(), dn.data(), ds.data());
    if (rc != CLAB_OK) return kerr(c, "updwn", rc);
    for (std::size_t i = 0; i < n; ++i) {
      k.src.push_back(sr[i] < 0 ? -1 : sr[i] + 1 + static_cast<int>(j0));
      k.side.push_back(sd[i]);
      k.dn.push_back(dn[i]);
      k.dst.push_back(ds[i]);
      k.row.push_back(static_cast<int>(i) + 1);
    }
    return true;
  };
  auto run = [&](const std::vector<double> &a0, const std::vector<double> &b0, const std::vector<double> &a1,
                 const std::vector<double> &b1, Run &res) {
    K k;
    std::vector<double> tb(7, 0.0);
    if (twin == 2) {
      for (std::size_t j = 0; j < ns; ++j) {
        if (!kern(a0, b0, a1, b1, j, j + 1, 0, k)) return false;
      }
      // an area near no source counts once in the table
      std::vector<char> cov(n, 0);
      for (std::size_t i = 0; i < k.dn.size(); ++i) {
        if (k.dn[i] != -2) {
          tb[li(k.dn[i])] += 1;
          cov[zz(k.row[i] - 1)] = 1;
        }
      }
      for (char v : cov) tb[3] += v ? 0 : 1;
    } else {
      if (!kern(a0, b0, a1, b1, 0, ns, twin == 1 ? 1 : 0, k)) return false;
      for (int v : k.dn) tb[li(v)] += 1;
    }
    Ix use;
    for (std::size_t i = 0; i < k.dn.size(); ++i) {
      if (k.dn[i] >= 0) use.push_back(i);
    }
    if (use.size() < 2) return ckstop(c, "fewer than two rows are upwind or downwind of a source");
    const Cols rl = dist.empty() ? Cols() : cutlab(dist, true);
    Ix u2;
    Cols ring;
    for (std::size_t i : use) {
      if (dist.empty()) {
        u2.push_back(i);
        ring.emplace_back();
        continue;
      }
      const int b = bincode(k.dst[i], dist, true);
      if (b < 0) continue;
      u2.push_back(i);
      ring.push_back(rl[zz(b)]);
    }
    if (u2.empty()) return ckstop(c, "no row falls in the rings of `dist`");
    const std::size_t m = u2.size();
    Cols g(m), cell(m), sdy(m);
    Ix rw(m);
    for (std::size_t i = 0; i < m; ++i) {
      const std::size_t a = u2[i];
      rw[i] = zz(k.row[a] - 1);
      const std::string sr = std::to_string(k.src[a]);
      g[i] = line[zz(k.src[a] - 1)] && agg ? "s\r" + sr + "\r" + std::to_string(k.side[a]) + "\r" + ring[i]
                                            : "a\r" + (area.empty() ? std::string() : av[rw[i]]) + "\r" + sr + "\r" +
                                                  ring[i];
      cell[i] = g[i] + "\r" + dv[rw[i]];
      sdy[i] = sr + "\r" + dv[rw[i]];
    }
    // whether a value differs within a cell
    auto same = [&](const Cols &v) {
      std::unordered_map<std::string, std::string> h;
      for (std::size_t i = 0; i < m; ++i) {
        const auto it = h.emplace(cell[i], v[i]);
        if (!it.second && it.first->second != v[i]) return true;
      }
      return false;
    };
    {
      Cols dn(m);
      for (std::size_t i = 0; i < m; ++i) dn[i] = std::to_string(k.dn[u2[i]]);
      if (same(dn)) {
        return ckstop(c, "the wind differs between areas of one side and day: give one direction per day, or agg = "
                         "FALSE");
      }
    }
    Cols cw = clus;
    cw.insert(cw.end(), wx.begin(), wx.end());
    for (const std::string &cc : cw) {
      const Cols l = coll(*fc(d, cc));
      Cols v(m);
      for (std::size_t i = 0; i < m; ++i) v[i] = l[rw[i]];
      if (same(v)) return ckstop(c, "column `%s` must be constant within each side (area) and day", cc.c_str());
    }
    // one row per cell, in order of first appearance; y summed
    std::unordered_map<std::string, std::size_t> hc;
    Ix o, oc(m);
    for (std::size_t i = 0; i < m; ++i) {
      const auto it = hc.emplace(cell[i], o.size());
      if (it.second) o.push_back(i);
      oc[i] = it.first->second;
    }
    const std::size_t ne = o.size();
    std::vector<double> ey(ne, 0.0), down(ne), rown(ne);
    for (std::size_t i = 0; i < m; ++i) ey[oc[i]] += yv[rw[i]];
    Cols eg(ne), esd(ne);
    for (std::size_t j = 0; j < ne; ++j) {
      eg[j] = g[o[j]];
      esd[j] = sdy[o[j]];
      down[j] = k.dn[u2[o[j]]];
      rown[j] = static_cast<double>(j + 1);
    }
    Val e = vobj();
    put(e, "y", vnums(ey));
    put(e, "g", vstrs(eg));
    put(e, "sd", vstrs(esd));
    put(e, "down", vnums(down));
    put(e, "row", vnums(rown));
    if (tnum) {
      std::vector<double> tv(ne);
      for (std::size_t j = 0; j < ne; ++j) tv[j] = nv(dc.a[rw[o[j]]]);
      setc(e, "t", vnums(tv));
    }
    for (const std::string &cc : clus) {
      Val a = varr();
      for (std::size_t j = 0; j < ne; ++j) a.a.push_back(fc(d, cc)->a[rw[o[j]]]);
      setc(e, cc, a);
    }
    Cols xs = {"down"};
    if (!dist.empty()) {
      xs.clear();
      for (const std::string &l : rl) {
        bool in = false;
        for (std::size_t j = 0; j < ne; ++j) in = in || ring[o[j]] == l;
        if (!in) continue;
        std::vector<double> v(ne);
        for (std::size_t j = 0; j < ne; ++j) v[j] = down[j] * (ring[o[j]] == l ? 1.0 : 0.0);
        setc(e, "down:" + l, vnums(v));
        xs.push_back("down:" + l);
      }
    }
    if (!wx.empty()) {
      // each side's own response to the weather, less one side per source
      Cols gl;
      std::vector<int> gs;
      for (std::size_t j = 0; j < ne; ++j) {
        if (!has(gl, eg[j])) {
          gl.push_back(eg[j]);
          gs.push_back(k.src[u2[o[j]]]);
        }
      }
      for (const std::string &w : wx) {
        const std::vector<double> wv0 = colv(*fc(d, w));
        for (std::size_t l = 0; l < gl.size(); ++l) {
          bool dup = false;
          for (std::size_t h = 0; h < l; ++h) dup = dup || gs[h] == gs[l];
          if (!dup) continue;
          std::vector<double> v(ne);
          for (std::size_t j = 0; j < ne; ++j) v[j] = wv0[rw[o[j]]] * (eg[j] == gl[l] ? 1.0 : 0.0);
          setc(e, w + ":side" + std::to_string(l + 1), vnums(v));
        }
      }
    }
    Cols xa = xs;
    for (const std::string &nm : e.k) {
      if (nm == "y" || nm == "g" || nm == "sd" || nm == "down" || nm == "row" || nm == "t" || has(clus, nm) ||
          has(xa, nm)) {
        continue;
      }
      xa.push_back(nm);
    }
    const Cols cl = se == 0 ? Cols{"row"} : clus.empty() ? Cols{"g"} : clus;
    if (norm != 3) {
      std::vector<double> yy = ey;
      if (norm == 0) {
        const double mu = lmean(yy);
        if (!(mu > 0)) return ckstop(c, "the counts are all 0");
        for (double &v : yy) v /= mu;
      } else if (norm == 2) {
        for (double &v : yy) v = std::log1p(v);
      }
      setc(e, "y", vnums(yy));
    }
    Val ar = vobj();
    put(ar, "data", e);
    put(ar, "y", vstr("y"));
    put(ar, "x", vstrs(xa));
    put(ar, "fe", vstrs({"g", "sd"}));
    put(ar, "clus", vstrs(cl));
    put(ar, "se", vstr(ses[zz(se)]));
    put(ar, "ssc", vstr(sscs[zz(ssc)]));
    if (se == 3) put(ar, "time", vstr("t"));
    for (const char *a : {"lag", "nboot", "seed"}) {
      if (c.has(a)) put(ar, a, c[a]);
    }
    put(ar, "level", vnum(lv));
    if (norm == 3) {
      CK(call(c, "ppmlfe", ppmlfe, ar, res.fit));
    } else {
      CK(call(c, "ivfe", ivfe, ar, res.fit));
    }
    Val tab = vobj();
    for (std::size_t j = 0; j < 7; ++j) put(tab, labs[j], vnum(tb[j]));
    res.tab = tab;
    res.n = static_cast<double>(ne);
    res.coef = terms(col(res.fit, "coef"), xs);
    return true;
  };
  Run r0;
  CK(run(s0, t0, s1, t1, r0));
  Val pl;
  if (!offs.empty()) {
    std::vector<double> nx(ns), ny(ns);
    for (std::size_t j = 0; j < ns; ++j) {
      const double dx = s1[j] - s0[j], dy = t1[j] - t0[j];
      const double len = std::sqrt(dx * dx + dy * dy);
      nx[j] = -dy / len;
      ny[j] = dx / len;
    }
    std::vector<double> po, pe, ps, plo, phi, pp, pn;
    for (double o : offs) {
      std::vector<double> a0(ns), b0(ns), a1(ns), b1(ns);
      for (std::size_t j = 0; j < ns; ++j) {
        a0[j] = s0[j] + o * nx[j];
        b0[j] = t0[j] + o * ny[j];
        a1[j] = s1[j] + o * nx[j];
        b1[j] = t1[j] + o * ny[j];
      }
      Run ro;
      po.push_back(o);
      // a moved segment with areas on one side only has nothing to compare
      if (!run(a0, b0, a1, b1, ro)) {
        c.err.clear();
        wrn(c, "no estimate for the segments moved by " + rfmt7(o) +
                   ": too few areas upwind and downwind of them.");
        for (auto *v : {&pe, &ps, &plo, &phi, &pp}) v->push_back(NA);
        pn.push_back(0);
        continue;
      }
      auto first = [&](const char *k) {
        const Val &a = col(ro.coef, k);
        return a.a.empty() ? NA : nv(a.a[0]);
      };
      pe.push_back(first("est"));
      ps.push_back(first("se"));
      plo.push_back(first("lo"));
      phi.push_back(first("hi"));
      pp.push_back(first("p"));
      pn.push_back(ro.n);
    }
    Frame f;
    f.num("off", po);
    f.num("est", pe);
    f.num("se", ps);
    f.num("lo", plo);
    f.num("hi", phi);
    f.num("p", pp);
    f.num("n", pn);
    pl = f.o;
  }
  out = vobj();
  put(out, "coef", r0.coef);
  put(out, "class", r0.tab);
  put(out, "n", vnum(r0.n));
  put(out, "fit", r0.fit);
  put(out, "placebo", pl);
  put(out, "norm", vstr(norms[zz(norm)]));
  put(out, "tol", vnum(tol));
  put(out, "band", vnum(band));
  return true;
}

}  // namespace ax
