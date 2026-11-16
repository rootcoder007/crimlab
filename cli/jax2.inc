// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Exposure and crime, continued: expgmm() of R/apx.R, whose design and fit
// are clab_expgmmd (src/apx.cpp), which the R function calls too, and the
// functions of R/apx.R and R/apx2.R that take an estimator function. In R
// `est` is a function and its arguments come through `...`; here `est`
// names one of ppmlfe, ivfe and ppmlcf (ppmlfe when left out) and its
// arguments are given beside the function's own, as in an R call. An
// argument the estimator does not take is refused as R refuses it.

namespace ax {

// expgmm()
bool expgmm(Ctx &c, Val &out) {
  static const std::vector<const char *> errs = {"mult", "add"}, steps = {"two", "one", "iter"},
                                         ses = {"cluster", "hc1"}, gsscs = {"G", "none"};
  int err = 0, step = 0, se = 0, ssc = 0;
  CK(marg(c, "err", errs, err));
  CK(marg(c, "step", steps, step));
  CK(marg(c, "se", ses, se));
  CK(marg(c, "ssc", gsscs, ssc));
  Val d;
  CK(dat(c, c["data"], "data", d));
  Cols y, endo, inst, x, fe, expo, clus;
  CK(acol(c, d, "y", {1}, false, y));
  CK(acol(c, d, "endo", {}, false, endo));
  CK(acol(c, d, "inst", {}, false, inst));
  CK(acol(c, d, "x", {}, true, x));
  CK(acol(c, d, "fe", {1, 2}, true, fe));
  CK(acol(c, d, "expo", {1}, true, expo));
  if (se == 0) {
    Nm cn = nms(c, {"clus"});
    if (cn.null && !fe.empty()) cn = nmof(Cols{fe[0]});
    if (cn.null) return ckstop(c, "`se = \"cluster\"` needs `clus` (or `fe`, its first column)");
    CK(acols(c, d, cn, "clus", {1}, false, clus));
  }
  double lv = 0;
  CK(alevel(c, lv));
  Nv tl, mx;
  CK(cknumd(c, "tol", 1e-10, tl, -HUGE_VAL, HUGE_VAL, true));
  CK(cknumd(c, "maxit", 100, mx, 1, INT_MAX, false, true));
  if (inst.size() < endo.size()) {
    return ckstop(c, "`inst` must give at least as many instruments as `endo` has columns");
  }
  for (const std::string &v : inst) {
    if (has(endo, v) || has(x, v)) return ckstop(c, "a column cannot be both a regressor and an excluded instrument");
  }
  // .axgmrows()
  std::vector<double> yv;
  CK(amat(c, d, y, yv));
  bool any = false;
  for (double v : yv) {
    if (v < 0) return ckstop(c, "column `%s` must hold counts, at least 0", y[0].c_str());
    any = any || v > 0;
  }
  if (!any) return ckstop(c, "column `%s` is all 0", y[0].c_str());
  std::vector<int> k(yv.size(), 1);
  if (!fe.empty()) CK(asep(c, yv, d, fe, k));
  const Ix kr = which(k);
  d = rows(d, kr);
  std::vector<double> yk, off;
  for (std::size_t i : kr) yk.push_back(yv[i]);
  if (!expo.empty()) {
    CK(amat(c, d, expo, off));
    for (double &v : off) {
      if (v <= 0) return ckstop(c, "column `%s` must be positive", expo[0].c_str());
      v = std::log(v);
    }
  }
  // .axgmdes()
  int nl = 0;
  std::vector<Code> fc2;
  for (const std::string &v : fe) {
    fc2.push_back(acode(*fc(d, v)));
    nl += fc2.back().n;
  }
  if (nl > 200) {
    return ckstop(c, "the fixed effects have %d levels; they enter as indicator columns, at most 200 (use ppmlcf() "
                     "for many)", nl);
  }
  Cols tm = endo, zn = inst;
  tm.insert(tm.end(), x.begin(), x.end());
  zn.insert(zn.end(), x.begin(), x.end());
  std::vector<double> xm, zm;
  CK(amat(c, d, tm, xm));
  CK(amat(c, d, zn, zm));
  const std::size_t n = yk.size();
  std::vector<int> cc;
  int nc = static_cast<int>(n);
  if (se == 0) {
    Cl cl;
    CK(acl(c, d, clus, cl));
    cc = cl.a.c;
    nc = cl.a.n;
  } else {
    for (std::size_t i = 0; i < n; ++i) cc.push_back(static_cast<int>(i));
  }
  const int nz = static_cast<int>(zn.size()) + (fe.empty() ? 1 : nl - static_cast<int>(fe.size()) + 1);
  if (nz > nc) {
    return ckstop(c, "%d moment conditions but %d clusters: the weight matrix would be singular (fewer instruments, "
                     "controls or effects, or more clusters)", nz, nc);
  }
  // .axgmfit()
  const int px = static_cast<int>(tm.size());
  const int kk = px + (fe.empty() ? 1 : nl - static_cast<int>(fe.size()) + 1);
  std::vector<double> beta(zz(kk)), vc(zz(kk) * zz(kk)), st(6);
  const int rc = clab_expgmmd(yk.data(), xm.data(), px, zm.data(), static_cast<int>(zn.size()), static_cast<int>(n),
                              fe.empty() ? nullptr : fc2[0].c.data(), fe.empty() ? 0 : fc2[0].n,
                              fe.size() < 2 ? nullptr : fc2[1].c.data(), fe.size() < 2 ? 0 : fc2[1].n,
                              off.empty() ? nullptr : off.data(), cc.data(), nc, err, step == 1 ? 1 : step == 0 ? 2 : 3,
                              static_cast<int>(mx[0]), tl[0], beta.data(), vc.data(), st.data());
  if (rc == CLAB_EDOM) {
    return c.fail("the fit failed: the instruments or the regressors are collinear, there are no more clusters than "
                  "moment conditions (the weight matrix is singular), or the mean overflows from the start values.");
  }
  if (rc != CLAB_OK) return c.fail("%s", clab_strerror(rc));
  if (st[3] != 1) wrn(c, "the fit did not converge in " + rfmt(mx[0]) + " iterations.");
  // .axgmtab()
  const double fac = ssc == 0 ? static_cast<double>(nc) / (nc - 1) : 1.0;
  std::vector<double> v(zz(px) * zz(px));
  for (std::size_t a = 0; a < zz(px); ++a) {
    for (std::size_t b = 0; b < zz(px); ++b) v[b * zz(px) + a] = vc[b * zz(kk) + a] * fac;
  }
  std::vector<double> b(beta.begin(), beta.begin() + px);
  Frame j;
  j.num("stat", {st[0]});
  j.num("df", {st[1]});
  j.num("p", {std::isnan(st[0]) ? NA : clab_pchisq(st[0], st[1], 0)});
  out = vobj();
  put(out, "coef", atab(tm, b, sdiag(v, px), lv, NaN, true));
  put(out, "j", j.o);
  put(out, "vcov", vmat(v, px, px));
  put(out, "n", vnum(static_cast<double>(n)));
  put(out, "nsep", vnum(static_cast<double>(yv.size() - n)));
  put(out, "ncl", vnum(nc));
  put(out, "err", vstr(errs[zz(err)]));
  put(out, "step", vstr(steps[zz(step)]));
  put(out, "se", vstr(ses[zz(se)]));
  put(out, "ssc", vstr(gsscs[zz(ssc)]));
  put(out, "iter", vnum(st[2]));
  put(out, "conv", vbool(st[3] == 1));
  put(out, "nw", vnum(st[4]));
  put(out, "fe", fe.empty() ? Val() : vs(fe));
  put(out, "clus", clus.empty() ? Val() : vs(clus));
  put(out, "level", vnum(lv));
  return true;
}

}  // namespace ax

namespace ax {

// The estimators `est` may name, and the arguments each takes.
const char *const ests[] = {"ppmlfe", "ivfe", "ppmlcf"};
const char *const eargs[] = {
    "data y x fe expo clus se ssc time lag nboot seed sep level tol maxit twc psd",
    "data y endo inst x fe clus se ssc time lag nboot seed level tol twc psd",
    "data y endo inst x fe expo clus how se ssc nboot seed level tol maxit"};
const Fn efns[] = {ppmlfe, ivfe, ppmlcf};

// .axest(): `est` left out is ppmlfe.
inline bool aest(Ctx &c, int &e) {
  e = 0;
  const Val &o = c["est"];
  if (o.t == Val::NUL) return true;
  if (o.t == Val::STR) {
    for (int i = 0; i < 3; ++i) {
      if (o.s == ests[i]) {
        e = i;
        return true;
      }
    }
  }
  return ckstop(c, "`est` must be an estimator function, such as ppmlfe");
}
// list(...): the arguments given that are not the function's own.
inline Val dots(Ctx &c, std::initializer_list<const char *> own) {
  Val a = vobj();
  for (std::size_t i = 0; i < c.in->k.size(); ++i) {
    bool mine = false;
    for (const char *o : own) mine = mine || c.in->k[i] == o;
    if (!mine) put(a, c.in->k[i].c_str(), c.in->a[i]);
  }
  return a;
}
// .clas(): a message of another function of the package as this one's.
inline std::string clas(const char *fn, const std::string &m) {
  const std::size_t k = m.find(": ");
  if (k == std::string::npos || k == 0) return m;
  for (std::size_t i = 0; i < k; ++i) {
    const char ch = m[i];
    if (!((ch >= 'a' && ch <= 'z') || (i > 0 && ch >= '0' && ch <= '9'))) return m;
  }
  return m.compare(0, k, fn) == 0 ? m : fn + m.substr(k);
}
// Another function of the package called as R calls it inside .clas().
inline bool ccall(Ctx &c, const char *fn, Fn f, const Val &a, Val &out) {
  Ctx s{fn, &a, std::string(), {}};
  const bool ok = f(s, out);
  for (const std::string &w : s.warn) c.warn.push_back(clas(c.fn, w));
  if (!ok) c.err = clas(c.fn, s.err);
  return ok;
}
// est(...) with the arguments a: one it does not take is refused as R
// refuses it.
inline bool ecall(Ctx &c, int e, const Val &a, Val &out) {
  const std::vector<std::string> ok = words(eargs[e]);
  std::string bad;
  int nb = 0;
  for (std::size_t i = 0; i < a.k.size(); ++i) {
    if (std::find(ok.begin(), ok.end(), a.k[i]) == ok.end()) {
      bad += (nb++ ? ", " : "") + a.k[i] + " = " + jdep(a.a[i]);
    }
  }
  if (nb > 0) {
    c.err = std::string(nb > 1 ? "unused arguments (" : "unused argument (") + bad + ")";
    return false;
  }
  return ccall(c, ests[e], efns[e], a, out);
}
// An argument object: data and the named values first, then the others.
inline Val eargs0(const Val &d, std::initializer_list<std::pair<const char *, Val>> first, const Val &rest) {
  Val a = vobj();
  put(a, "data", d);
  for (const auto &p : first) put(a, p.first, p.second);
  for (std::size_t i = 0; i < rest.k.size(); ++i) put(a, rest.k[i].c_str(), rest.a[i]);
  return a;
}
// c(v, w) of names: text or arrays of text.
inline Val cat(const Val &v, const Val &w) {
  Val r = varr();
  for (const Val *o : {&v, &w}) {
    if (o->t == Val::ARR) r.a.insert(r.a.end(), o->a.begin(), o->a.end());
    else if (o->t != Val::NUL) r.a.push_back(*o);
  }
  return r;
}
// The rows of a fit's coefficients whose term is t, in the order of t
// (match(): a term not there is a row of NA).
inline Val trows(const Val &cf, const Cols &t) {
  const Val &tm = col(cf, "term");
  Val o = vobj();
  for (std::size_t j = 0; j < cf.k.size(); ++j) {
    Val a = varr();
    for (const std::string &s : t) {
      std::size_t i = 0;
      while (i < tm.a.size() && !(tm.a[i].t == Val::STR && tm.a[i].s == s)) ++i;
      a.a.push_back(i < tm.a.size() ? cf.a[j].a[i] : Val());
    }
    put(o, cf.k[j].c_str(), a);
  }
  return o;
}
inline std::vector<double> cnum(const Val &f, const char *k) { return colv(col(f, k)); }
// .axvcov(): the covariance of a fit, of the terms t (in the order of
// its coefficients); from the bootstrap draws when it has none.
inline std::vector<double> avcov(const Val &r, const Cols &t) {
  const Val &tm = col(col(r, "coef"), "term");
  const std::size_t p = tm.a.size(), k = t.size();
  std::vector<double> v(p * p, NA), o(k * k, NA);
  const Val &vc = col(r, "vcov");
  if (vc.t == Val::ARR && vc.a.size() == p) {
    for (std::size_t i = 0; i < p; ++i) {
      for (std::size_t j = 0; j < p; ++j) v[j * p + i] = nv(vc.a[i].a[j]);
    }
  } else {
    const Val &dr = col(r, "draws");
    const std::size_t nb = dr.a.size();
    std::vector<double> x(nb * p);
    std::vector<char> ok(nb, 1);
    for (std::size_t b = 0; b < nb; ++b) {
      for (std::size_t j = 0; j < p; ++j) {
        x[j * nb + b] = nv(dr.a[b].a[j]);
        if (std::isnan(x[j * nb + b])) ok[b] = 0;
      }
    }
    v = rcov(x, nb, static_cast<int>(p), ok);
  }
  std::vector<std::size_t> ix;
  for (const std::string &s : t) {
    std::size_t i = 0;
    while (i < p && !(tm.a[i].t == Val::STR && tm.a[i].s == s)) ++i;
    ix.push_back(i);
  }
  for (std::size_t a = 0; a < k; ++a) {
    for (std::size_t b = 0; b < k; ++b) {
      if (ix[a] < p && ix[b] < p) o[b * k + a] = v[ix[b] * p + ix[a]];
    }
  }
  return o;
}
// .axwald(): R b = 0 for b (k), v (k x k) and r (q x k, column-major).
inline Val awald(Ctx &c, const std::vector<double> &b, const std::vector<double> &v, const std::vector<double> &r,
                 int q) {
  const std::size_t k = b.size(), uq = zz(q);
  std::vector<double> rb(uq, 0.0), rv(uq * k, 0.0), m(uq * uq, 0.0);
  for (std::size_t i = 0; i < uq; ++i) {
    for (std::size_t j = 0; j < k; ++j) rb[i] += r[j * uq + i] * b[j];
  }
  for (std::size_t i = 0; i < uq; ++i) {
    for (std::size_t j = 0; j < k; ++j) {
      for (std::size_t l = 0; l < k; ++l) rv[j * uq + i] += r[l * uq + i] * v[j * k + l];
    }
  }
  for (std::size_t i = 0; i < uq; ++i) {
    for (std::size_t j = 0; j < uq; ++j) {
      for (std::size_t l = 0; l < k; ++l) m[j * uq + i] += rv[l * uq + i] * r[l * uq + j];
    }
  }
  std::vector<double> mq = m;
  std::vector<int> piv;
  const int rk = qrank(mq, q, q, 1e-7, piv);
  Frame f;
  if (rk < q) {
    // .axwald() warns without a function name
    c.warn.push_back("the joint test of " + std::to_string(q) + " restrictions is not computed: the covariance has rank " +
                     std::to_string(rk) + " (too few clusters for it).");
    f.num("stat", {NA});
    f.num("df", {static_cast<double>(q)});
    f.num("p", {NA});
    return f.o;
  }
  // solve(m, rb) by elimination with partial pivoting, as LAPACK's dgesv
  std::vector<double> x = rb;
  for (std::size_t j = 0; j < uq; ++j) {
    std::size_t pv = j;
    for (std::size_t i = j + 1; i < uq; ++i) {
      if (std::fabs(m[j * uq + i]) > std::fabs(m[j * uq + pv])) pv = i;
    }
    if (pv != j) {
      for (std::size_t l = 0; l < uq; ++l) std::swap(m[l * uq + j], m[l * uq + pv]);
      std::swap(x[j], x[pv]);
    }
    for (std::size_t i = j + 1; i < uq; ++i) {
      const double t = m[j * uq + i] / m[j * uq + j];
      for (std::size_t l = j; l < uq; ++l) m[l * uq + i] -= t * m[l * uq + j];
      x[i] -= t * x[j];
    }
  }
  for (std::size_t j = uq; j-- > 0;) {
    for (std::size_t l = j + 1; l < uq; ++l) x[j] -= m[l * uq + j] * x[l];
    x[j] /= m[j * uq + j];
  }
  double w = 0.0;
  for (std::size_t i = 0; i < uq; ++i) w += rb[i] * x[i];
  f.num("stat", {w});
  f.num("df", {static_cast<double>(q)});
  f.num("p", {clab_pchisq(w, q, 0)});
  return f.o;
}
// quantile(x, p) of type 7.
inline std::vector<double> aqnt(const std::vector<double> &x, const std::vector<double> &p) {
  std::vector<double> q(p.size());
  clab_quant(x.data(), static_cast<int>(x.size()), p.data(), static_cast<int>(p.size()), q.data());
  return q;
}
// The rows of the coefficient table of a fit for term t, or R's refusal.
inline bool arow(Ctx &c, const Val &r, const std::string &t, const char *msg, std::size_t &i) {
  const Val &tm = col(col(r, "coef"), "term");
  for (i = 0; i < tm.a.size(); ++i) {
    if (tm.a[i].t == Val::STR && tm.a[i].s == t) return true;
  }
  return ckstop(c, msg, t.c_str());
}
inline double cell(const Val &r, const char *k, std::size_t i) { return nv(col(col(r, "coef"), k).a[i]); }
// .axrcs(): x and the k - 2 cubic terms of the restricted cubic spline at
// the knots t, column-major n x (k - 1); x alone for fewer than 3 knots.
inline std::vector<double> rcs(const std::vector<double> &x, const std::vector<double> &t, int &nc) {
  const std::size_t n = x.size(), k = t.size();
  nc = k < 3 ? 1 : static_cast<int>(k) - 1;
  std::vector<double> m(n * zz(nc));
  std::copy(x.begin(), x.end(), m.begin());
  if (k < 3) return m;
  const auto p = [](double u) { return std::pow(u > 0 ? u : 0.0, 3.0); };
  const double tk = t[k - 1], tk1 = t[k - 2], d2 = rpow(tk - t[0], 2);
  for (std::size_t j = 0; j + 2 < k; ++j) {
    for (std::size_t i = 0; i < n; ++i) {
      m[(j + 1) * n + i] = (p(x[i] - t[j]) - p(x[i] - tk1) * (tk - t[j]) / (tk - tk1) +
                            p(x[i] - tk) * (tk1 - t[j]) / (tk - tk1)) /
                           d2;
    }
  }
  return m;
}
// seq(a, b, length.out = n), as R makes it.
inline std::vector<double> seqn(double a, double b, int n) {
  std::vector<double> s;
  if (n == 1) return {a};
  if (a == b) return std::vector<double>(zz(n), a);
  const int n1 = n - 1;
  s.push_back(a);
  for (int i = 1; i < n1; ++i) s.push_back(a + i * ((b - a) / n1));
  s.push_back(b);
  return s;
}
// clab_dshift() of codes u and days t by k, the message of R's binding.
inline bool ashift(Ctx &c, const std::vector<int> &u, const std::vector<double> &t, double k, std::vector<int> &ix) {
  const std::size_t n = t.size();
  std::vector<int> ord(n);
  ix.assign(n, -1);
  const int rc = clab_dshift(u.data(), t.data(), static_cast<int>(n), k, ix.data(), ord.data());
  if (rc == CLAB_EDOM) return ckstop(c, "a unit has two rows on one day: give one row per unit and day");
  if (rc != CLAB_OK) return c.fail("%s", clab_strerror(rc));
  return true;
}
inline Val pick(const Val &a, const std::vector<int> &ix) {
  Val o = varr();
  for (int i : ix) o.a.push_back(i < 0 ? Val() : a.a[zz(i)]);
  return o;
}

}  // namespace ax

// runfx()
bool runfx(Ctx &c, Val &out) {
  using namespace ax;
  int e = 0;
  CK(aest(c, e));
  Val ra = vobj(), r;
  for (const char *k : {"data", "var", "unit", "day", "thr", "cap", "step"}) {
    if (c.has(k)) put(ra, k, c[k]);
  }
  CK(ccall(c, "exprun", exprun, ra, r));
  Val d;
  CK(dat(c, c["data"], "data", d));
  Cols rc(r.k.begin() + 1, r.k.end());
  for (std::size_t j = 1; j < r.k.size(); ++j) setc(d, r.k[j], r.a[j]);
  const Val a = eargs0(d, {{"y", c["y"]}, {"x", cat(vstrs(rc), c["ctrl"])}},
                       dots(c, {"data", "y", "var", "unit", "day", "thr", "cap", "step", "ctrl", "est"}));
  Val fit;
  CK(ecall(c, e, a, fit));
  Ix keep;
  const Val &tm = col(col(fit, "coef"), "term");
  for (std::size_t i = 0; i < tm.a.size(); ++i) {
    if (tm.a[i].t == Val::STR && ax::has(rc, tm.a[i].s)) keep.push_back(i);
  }
  out = vobj();
  put(out, "coef", rows(col(fit, "coef"), keep));
  put(out, "fit", fit);
  return true;
}

// hetfx()
bool hetfx(Ctx &c, Val &out) {
  using namespace ax;
  int how = 0;
  CK(marg(c, "how", {"inter", "split"}, how));
  Val d;
  CK(dat(c, c["data"], "data", d));
  Cols vm;
  CK(acols(c, d, nms(c, {"var", "mod"}), "var, mod", {2}, false, vm));
  int e = 0;
  CK(aest(c, e));
  const Cols mv = coll(*fc(d, vm[1]));
  Cols lv = mv;
  std::sort(lv.begin(), lv.end());
  lv.erase(std::unique(lv.begin(), lv.end()), lv.end());
  if (lv.size() < 2) return ckstop(c, "column `%s` has one level: nothing to compare", vm[1].c_str());
  Val a = dots(c, {"data", "var", "mod", "est", "how"});
  const std::size_t n = mv.size(), k = lv.size();
  std::vector<double> b, v(k * k, 0.0), se, lo, hi, p, nn;
  if (how == 1) {
    for (const std::string &l : lv) {
      Ix rw;
      for (std::size_t i = 0; i < n; ++i) {
        if (mv[i] == l) rw.push_back(i);
      }
      Val r;
      CK(ecall(c, e, eargs0(rows(d, rw), {}, a), r));
      std::size_t i = 0;
      CK(arow(c, r, vm[0], "the estimator returned no coefficient for `%s`", i));
      b.push_back(cell(r, "est", i));
      se.push_back(cell(r, "se", i));
      lo.push_back(cell(r, "lo", i));
      hi.push_back(cell(r, "hi", i));
      p.push_back(cell(r, "p", i));
      nn.push_back(nv(col(r, "n")));
    }
    for (std::size_t j = 0; j < k; ++j) v[j * k + j] = rpow(se[j], 2);
  } else {
    const std::vector<double> xv = colv(*fc(d, vm[0]));
    Cols nm;
    for (const std::string &l : lv) {
      std::vector<double> w(n);
      for (std::size_t i = 0; i < n; ++i) w[i] = xv[i] * (mv[i] == l ? 1.0 : 0.0);
      nm.push_back(vm[0] + ":" + l);
      setc(d, nm.back(), vnums(w));
    }
    // .axrepl(): var by its interactions in x and endo, each instrument by
    // its own
    for (std::size_t j = 0; j < a.k.size(); ++j) {
      if (a.k[j] != "x" && a.k[j] != "endo") continue;
      const Val ns = cat(a.a[j], Val());
      bool hv = false;
      Val o = varr();
      for (const Val &s : ns.a) {
        if (s.t == Val::STR && s.s == vm[0]) hv = true;
        else o.a.push_back(s);
      }
      if (hv) a.a[j] = cat(o, vstrs(nm));
    }
    for (std::size_t j = 0; j < a.k.size(); ++j) {
      if (a.k[j] != "inst") continue;
      Cols zi;
      for (const Val &z : cat(a.a[j], Val()).a) {
        const Val *zc = z.t == Val::STR ? fc(d, z.s) : nullptr;
        if (zc == nullptr) continue;
        const std::vector<double> zv = colv(*zc);
        for (const std::string &l : lv) {
          std::vector<double> w(n);
          for (std::size_t i = 0; i < n; ++i) w[i] = zv[i] * (mv[i] == l ? 1.0 : 0.0);
          zi.push_back(z.s + ":" + l);
          setc(d, zi.back(), vnums(w));
        }
      }
      a.a[j] = vstrs(zi);
    }
    Val r;
    CK(ecall(c, e, eargs0(d, {}, a), r));
    const Val cf = trows(col(r, "coef"), nm);
    b = cnum(cf, "est");
    se = cnum(cf, "se");
    lo = cnum(cf, "lo");
    hi = cnum(cf, "hi");
    p = cnum(cf, "p");
    for (const std::string &l : lv) nn.push_back(static_cast<double>(std::count(mv.begin(), mv.end(), l)));
    v = avcov(r, nm);
  }
  Frame f;
  f.str("level", lv);
  f.num("est", b);
  f.num("se", se);
  f.num("lo", lo);
  f.num("hi", hi);
  f.num("p", p);
  f.num("n", nn);
  const int q = static_cast<int>(k) - 1;
  std::vector<double> rr(zz(q) * k, 0.0);
  for (std::size_t i = 0; i < zz(q); ++i) {
    rr[i * zz(q) + i] = 1.0;
    rr[(i + 1) * zz(q) + i] = -1.0;
  }
  out = vobj();
  put(out, "tab", f.o);
  put(out, "test", awald(c, b, v, rr, q));
  put(out, "how", vstr(how == 0 ? "inter" : "split"));
  return true;
}

// ccross(): the day a number of days since 1970-01-01, as a date is.
bool ccross(Ctx &c, Val &out) {
  using namespace ax;
  int strat = 0;
  CK(marg(c, "strat", {"ymw", "ym"}, strat));
  Val d;
  CK(dat(c, c["data"], "data", d));
  Cols yud, ma, lagv;
  CK(acols(c, d, nms(c, {"y", "unit", "day"}), "y, unit, day", {3}, false, yud));
  CK(acol(c, d, "ma", {}, true, ma));
  bool dt = false;
  CK(aday(c, d, yud[2], dt));
  if (!isnumc(*fc(d, yud[2]))) return ckstop(c, "column `%s` must hold dates", yud[2].c_str());
  int e = 0;
  CK(aest(c, e));
  Nv q;
  CK(cknumd(c, "k", 2, q, 1, Inf, false, true));
  const double k = q[0];
  CK(acol(c, d, "lagv", {}, true, lagv));
  CK(cknumd(c, "lag", 0, q, 0, Inf, false, true));
  const double lag = q[0];
  const Val &ns = c["ns"];
  Cols nsn;
  std::vector<double> nsd;
  if (ns.t != Val::NUL) {
    bool ok = ns.t == Val::OBJ && numv(ns);
    for (std::size_t j = 0; ok && j < ns.a.size(); ++j) {
      const double v = nv(ns.a[j]);
      ok = !std::isnan(v) && v >= 1 && v == std::nearbyint(v);
      nsn.push_back(ns.k[j]);
      nsd.push_back(v);
    }
    if (!ok) return ckstop(c, "`ns` must be whole degrees of freedom named by column");
    Cols t;
    CK(acols(c, d, nmof(nsn), "names(ns)", {}, false, t));
  }
  // .axccma()
  const auto keep = [&](const Cols &cs) {
    Ix r;
    for (std::size_t i = 0; i < nr(d); ++i) {
      bool ok = true;
      for (const std::string &s : cs) ok = ok && !vna(fc(d, s)->a[i]);
      if (ok) r.push_back(i);
    }
    d = rows(d, r);
  };
  if (!lagv.empty() && lag > 0) {
    std::vector<int> ix;
    CK(ashift(c, acode(*fc(d, yud[1])).c, colv(*fc(d, yud[2])), -lag, ix));
    for (const std::string &v : lagv) setc(d, v, pick(*fc(d, v), ix));
    keep(lagv);
  }
  if (!ma.empty() && k > 1) {
    const std::vector<int> uc = acode(*fc(d, yud[1])).c;
    const std::vector<double> tv = colv(*fc(d, yud[2]));
    std::vector<std::vector<double>> s;
    for (const std::string &v : ma) s.push_back(colv(*fc(d, v)));
    for (double j = 1; j <= k - 1; ++j) {
      std::vector<int> ix;
      CK(ashift(c, uc, tv, -j, ix));
      for (std::size_t m = 0; m < ma.size(); ++m) {
        const std::vector<double> mv = colv(*fc(d, ma[m]));
        for (std::size_t i = 0; i < ix.size(); ++i) s[m][i] += ix[i] < 0 ? NA : mv[zz(ix[i])];
      }
    }
    for (std::size_t m = 0; m < ma.size(); ++m) {
      for (double &v : s[m]) v /= k;
      setc(d, ma[m], vnums(s[m]));
    }
    keep(ma);
  }
  // the strata: place, year and month, and weekday
  const std::size_t n = nr(d);
  const Cols ul = coll(*fc(d, yud[1]));
  const std::vector<double> dv = colv(*fc(d, yud[2]));
  Cols st(n), ym(n);
  for (std::size_t i = 0; i < n; ++i) {
    // civil date of a day number (H. Hinnant's days_from_civil inverse)
    const std::int64_t z = static_cast<std::int64_t>(std::floor(dv[i])) + 719468;
    const std::int64_t era = (z >= 0 ? z : z - 146096) / 146097;
    const std::int64_t doe = z - era * 146097;
    const std::int64_t yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    const std::int64_t doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    const std::int64_t mp = (5 * doy + 2) / 153;
    const std::int64_t mo = mp < 10 ? mp + 3 : mp - 9;
    const std::int64_t yr = yoe + era * 400 + (mo <= 2 ? 1 : 0);
    char b[48];
    std::snprintf(b, sizeof b, "%04lld-%02lld", static_cast<long long>(yr), static_cast<long long>(mo));
    ym[i] = b;
    const std::int64_t zd = static_cast<std::int64_t>(std::floor(dv[i]));
    const std::int64_t wd = ((zd + 4) % 7 + 7) % 7;
    st[i] = ul[i] + "\r" + ym[i] + "\r" + (strat == 0 ? std::to_string(wd) : std::string());
  }
  setc(d, ".strat", vstrs(st));
  setc(d, ".ym", vstrs(ym));
  Val a = dots(c, {"data", "y", "unit", "day", "strat", "ma", "k", "lagv", "lag", "ns", "est"});
  // .axccns()
  for (std::size_t j = 0; j < nsn.size(); ++j) {
    const Val &xc = *fc(d, nsn[j]);
    if (!isnumc(xc)) return ckstop(c, "column `%s` must be numeric", nsn[j].c_str());
    const std::vector<double> x = colv(xc);
    const int dg = static_cast<int>(nsd[j]);
    std::vector<double> pr;
    for (int i = 1; i < dg; ++i) pr.push_back(static_cast<double>(i) / dg);
    std::vector<double> kn = {*std::min_element(x.begin(), x.end())};
    for (double v : aqnt(x, pr)) kn.push_back(v);
    kn.push_back(*std::max_element(x.begin(), x.end()));
    int nc = 0;
    const std::vector<double> bm = rcs(x, kn, nc);
    Cols nm;
    for (int i = 0; i < nc; ++i) {
      nm.push_back(nsn[j] + ":ns" + std::to_string(i + 1));
      setc(d, nm.back(), vnums(&bm[zz(i) * n], n));
    }
    bool hx = false;
    for (std::size_t i = 0; i < a.k.size(); ++i) {
      if (a.k[i] == "x") {
        a.a[i] = cat(a.a[i], vstrs(nm));
        hx = true;
      }
    }
    if (!hx) put(a, "x", vstrs(nm));
  }
  if (!c.has("clus")) {
    bool hc = false;
    for (const std::string &s : a.k) hc = hc || s == "clus";
    if (!hc) put(a, "clus", vstr(".ym"));
  }
  return ecall(c, e, eargs0(d, {{"y", c["y"]}, {"fe", vstr(".strat")}}, a), out);
}

// jntx()
bool jntx(Ctx &c, Val &out) {
  using namespace ax;
  int form = 0;
  CK(marg(c, "form", {"spline", "linear"}, form));
  Val d;
  CK(dat(c, c["data"], "data", d));
  // .axjxck()
  Cols yx;
  CK(acols(c, d, nms(c, {"y", "x1", "x2"}), "y, x1, x2", {3}, false, yx));
  int e = 0;
  CK(aest(c, e));
  Nv q, p1, p2;
  CK(cknumd(c, "k1", 5, q, 3, 7, false, true));
  const int k1 = static_cast<int>(q[0]);
  CK(cknumd(c, "k2", 4, q, 3, 7, false, true));
  const int k2 = static_cast<int>(q[0]);
  if (c.has("p1")) {
    CK(cknum(c, c["p1"], "p1", p1, false, 0, 1));
  } else {
    p1.v = {0.1, 0.9};
  }
  if (c.has("p2")) {
    CK(cknum(c, c["p2"], "p2", p2, false, 0, 1));
  } else {
    p2.v = {0.05, 0.5, 0.95};
  }
  if (p1.n() != 2) return ckstop(c, "`p1` must be two quantiles");
  if (p2.n() < 2) return ckstop(c, "`p2` must be two or more quantiles");
  CK(cknumd(c, "ngrid", 15, q, 2, Inf, false, true));
  if (q[0] > 8192) return ckstop(c, "`ngrid` times `ngrid` must be at most 67108864");
  const int ng = static_cast<int>(q[0]);
  Nv rf;
  if (c.has("ref")) {
    CK(cknum(c, c["ref"], "ref", rf, false));
    if (rf.n() != 2) return ckstop(c, "`ref` must be two numbers, x1 and x2");
  }
  std::vector<double> a1, a2;
  CK(amat(c, d, {yx[1]}, a1));
  CK(amat(c, d, {yx[2]}, a2));
  static const double pk[][7] = {{0.1, 0.5, 0.9},
                                 {0.05, 0.35, 0.65, 0.95},
                                 {0.05, 0.275, 0.5, 0.725, 0.95},
                                 {0.05, 0.23, 0.41, 0.59, 0.77, 0.95},
                                 {0.025, 0.1833, 0.3417, 0.5, 0.6583, 0.8167, 0.975}};
  const auto knots = [&](const std::vector<double> &x, int kk) {
    return aqnt(x, std::vector<double>(pk[kk - 3], pk[kk - 3] + kk));
  };
  const std::vector<double> t1 = knots(a1, k1), t2 = knots(a2, k2);
  const auto dup = [](std::vector<double> t) {
    std::sort(t.begin(), t.end());
    return std::adjacent_find(t.begin(), t.end()) != t.end();
  };
  if (form == 0 && (dup(t1) || dup(t2))) {
    return ckstop(c, "an exposure has too few distinct values for its knots");
  }
  // .axjxbas(): column-major, rows of u and v
  const auto basis = [&](const std::vector<double> &u, const std::vector<double> &v, int &nc) {
    const std::size_t m = u.size();
    std::vector<double> o;
    if (form == 1) {
      nc = 3;
      o = u;
      o.insert(o.end(), v.begin(), v.end());
      for (std::size_t i = 0; i < m; ++i) o.push_back(u[i] * v[i]);
      return o;
    }
    int c1 = 0, c2 = 0;
    const std::vector<double> b1 = rcs(u, t1, c1), b2 = rcs(v, t2, c2);
    o = b1;
    o.insert(o.end(), b2.begin(), b2.end());
    for (int j = 0; j < c2; ++j) {
      for (int l = 0; l < c1; ++l) {
        for (std::size_t i = 0; i < m; ++i) o.push_back(b1[zz(l) * m + i] * b2[zz(j) * m + i]);
      }
    }
    nc = c1 + c2 + c1 * c2;
    return o;
  };
  int nb = 0;
  const std::vector<double> bm = basis(a1, a2, nb);
  const int nb1 = form == 1 ? 1 : k1 - 1, nb2 = form == 1 ? 1 : k2 - 1;
  const std::size_t n = a1.size();
  Cols nm;
  for (int j = 0; j < nb; ++j) {
    nm.push_back(".j" + std::to_string(j + 1));
    setc(d, nm.back(), vnums(&bm[zz(j) * n], n));
  }
  Val r;
  CK(ecall(c, e,
           eargs0(d, {{"y", vstr(yx[0])}, {"x", cat(vstrs(nm), c["ctrl"])}},
                  dots(c, {"data", "y", "x1", "x2", "form", "k1", "k2", "p1", "p2", "ref", "ngrid", "ctrl", "est"})),
           r));
  const std::vector<double> b = cnum(trows(col(r, "coef"), nm), "est");
  const std::vector<double> v = avcov(r, nm);
  std::vector<double> ref = rf.v;
  if (ref.empty()) ref = {clab_med(a1.data(), static_cast<int>(n)), clab_med(a2.data(), static_cast<int>(n))};
  // rel(): the change of the linear predictor and its standard error
  const auto rel = [&](const std::vector<double> &u, const std::vector<double> &w, const std::vector<double> &u0,
                       const std::vector<double> &w0, std::vector<double> &es, std::vector<double> &ss) {
    int nc = 0;
    const std::vector<double> x = basis(u, w, nc), x0 = basis(u0, w0, nc);
    const std::size_t m = u.size(), kk = b.size();
    es.assign(m, 0.0);
    ss.assign(m, 0.0);
    std::vector<double> db(kk);
    for (std::size_t i = 0; i < m; ++i) {
      for (std::size_t j = 0; j < kk; ++j) db[j] = x[j * m + i] - x0[j * m + i];
      for (std::size_t j = 0; j < kk; ++j) es[i] += db[j] * b[j];
      double s = 0.0;
      for (std::size_t j = 0; j < kk; ++j) {
        double t = 0.0;
        for (std::size_t l = 0; l < kk; ++l) t += db[l] * v[j * kk + l];
        s += t * db[j];
      }
      ss[i] = std::sqrt(s > 0 ? s : (std::isnan(s) ? s : 0.0));
    }
  };
  // .axjxout()
  const std::vector<double> q1 = aqnt(a1, p1.v), q2 = aqnt(a2, p2.v);
  const std::size_t m2 = q2.size();
  std::vector<double> ce, cs;
  rel(std::vector<double>(m2, q1[1]), q2, std::vector<double>(m2, q1[0]), q2, ce, cs);
  const double z = clab_qnorm(0.975);
  std::vector<double> clo, chi, crr;
  for (std::size_t i = 0; i < m2; ++i) {
    clo.push_back(ce[i] - z * cs[i]);
    chi.push_back(ce[i] + z * cs[i]);
    crr.push_back(std::exp(ce[i]));
  }
  Frame cf;
  cf.num("x2", q2);
  cf.num("p2", p2.v);
  cf.num("est", ce);
  cf.num("se", cs);
  cf.num("lo", clo);
  cf.num("hi", chi);
  cf.num("rr", crr);
  const auto mm1 = std::minmax_element(a1.begin(), a1.end()), mm2 = std::minmax_element(a2.begin(), a2.end());
  const std::vector<double> g1 = seqn(*mm1.first, *mm1.second, ng), g2 = seqn(*mm2.first, *mm2.second, ng);
  std::vector<double> gx, gy;
  for (double y2 : g2) {
    for (double x1 : g1) {
      gx.push_back(x1);
      gy.push_back(y2);
    }
  }
  std::vector<double> se2, ss2, sr;
  rel(gx, gy, std::vector<double>(gx.size(), ref[0]), std::vector<double>(gx.size(), ref[1]), se2, ss2);
  for (double x : se2) sr.push_back(std::exp(x));
  Frame sf;
  sf.num("x1", gx);
  sf.num("x2", gy);
  sf.num("est", se2);
  sf.num("se", ss2);
  sf.num("rr", sr);
  const int qn = nb - nb1 - nb2;
  std::vector<double> rr(zz(qn) * zz(nb), 0.0);
  for (int i = 0; i < qn; ++i) rr[zz(nb1 + nb2 + i) * zz(qn) + zz(i)] = 1.0;
  out = vobj();
  put(out, "cond", cf.o);
  put(out, "surf", sf.o);
  put(out, "test", awald(c, b, v, rr, qn));
  put(out, "fit", r);
  put(out, "form", vstr(form == 0 ? "spline" : "linear"));
  return true;
}

// dosebin()
bool dosebin(Ctx &c, Val &out) {
  using namespace ax;
  Val d;
  CK(dat(c, c["data"], "data", d));
  Cols yv;
  CK(acols(c, d, nms(c, {"y", "var"}), "y, var", {2}, false, yv));
  int e = 0;
  CK(aest(c, e));
  std::vector<double> xv;
  CK(amat(c, d, {yv[1]}, xv));
  // .axdbbin()
  std::vector<double> brk;
  Nv q;
  if (!c.has("brk")) {
    CK(cknumd(c, "nq", 5, q, 2, Inf, false, true));
    if (q[0] > kutil::maxcell) return ckstop(c, "`nq` must be at most 67108864");
    const int nq = static_cast<int>(q[0]);
    brk = aqnt(xv, seqn(0, 1, nq + 1));
    brk.erase(std::unique(brk.begin(), brk.end()), brk.end());
    if (brk.size() < 3) return ckstop(c, "`%s` gives fewer than two bins", yv[1].c_str());
  } else {
    const Val &o = c["brk"];
    brk = nums(o);
    if (!numv(o) || brk.size() < 3 || anyna(brk) || !incr(brk)) {
      return ckstop(c, "`brk` must be three or more increasing numbers");
    }
  }
  const Cols lv = cutlab(brk, true);
  const std::size_t n = xv.size(), nbin = lv.size();
  std::vector<int> bn(n);
  for (std::size_t i = 0; i < n; ++i) {
    bn[i] = bincode(xv[i], brk, true);
    if (bn[i] < 0) return ckstop(c, "some values of `%s` fall outside `brk`", yv[1].c_str());
  }
  CK(cknumd(c, "ref", 1, q, 1, static_cast<double>(nbin), false, true));
  const std::size_t ref = static_cast<std::size_t>(q[0]) - 1;
  Cols nm, use;
  std::vector<std::size_t> ui;
  for (std::size_t j = 0; j < nbin; ++j) {
    nm.push_back("bin:" + lv[j]);
    if (j == ref) continue;
    ui.push_back(j);
    use.push_back(nm[j]);
    std::vector<double> w(n);
    for (std::size_t i = 0; i < n; ++i) w[i] = bn[i] == static_cast<int>(j) ? 1.0 : 0.0;
    setc(d, nm[j], vnums(w));
  }
  Val r;
  CK(ecall(c, e,
           eargs0(d, {{"y", vstr(yv[0])}, {"x", cat(vstrs(use), c["ctrl"])}},
                  dots(c, {"data", "y", "var", "nq", "brk", "ref", "ctrl", "est"})),
           r));
  const Val cf = trows(col(r, "coef"), use);
  std::vector<double> mid(nbin), cnt(nbin), est(nbin, 0.0), se(nbin, NA), lo(nbin, NA), hi(nbin, NA), p(nbin, NA);
  for (std::size_t j = 0; j < nbin; ++j) {
    std::vector<double> in;
    for (std::size_t i = 0; i < n; ++i) {
      if (bn[i] == static_cast<int>(j)) in.push_back(xv[i]);
    }
    mid[j] = in.empty() ? NA : clab_med(in.data(), static_cast<int>(in.size()));
    cnt[j] = static_cast<double>(in.size());
  }
  const std::vector<double> ce = cnum(cf, "est"), cs = cnum(cf, "se"), cl = cnum(cf, "lo"), ch = cnum(cf, "hi"),
                            cp = cnum(cf, "p");
  for (std::size_t u = 0; u < ui.size(); ++u) {
    est[ui[u]] = ce[u];
    se[ui[u]] = cs[u];
    lo[ui[u]] = cl[u];
    hi[ui[u]] = ch[u];
    p[ui[u]] = cp[u];
  }
  Frame tb;
  tb.str("bin", lv);
  tb.num("mid", mid);
  tb.num("n", cnt);
  tb.num("est", est);
  tb.num("se", se);
  tb.num("lo", lo);
  tb.num("hi", hi);
  tb.num("p", p);
  // .axdbtst()
  const std::vector<double> v = avcov(r, use);
  const int k = static_cast<int>(ce.size());
  std::vector<double> id(zz(k) * zz(k), 0.0);
  for (std::size_t i = 0; i < zz(k); ++i) id[i * zz(k) + i] = 1.0;
  Val t1 = awald(c, ce, v, id, k);
  Cols tn = {"all"};
  std::vector<double> ts = {nv(col(t1, "stat").a[0])}, td = {static_cast<double>(k)}, tp = {nv(col(t1, "p").a[0])};
  if (k >= 2) {
    const int ql = k - 1;
    std::vector<double> rl(zz(ql) * zz(k), 0.0);
    for (int qq = 2; qq <= k; ++qq) {
      const std::size_t row = zz(qq - 2);
      rl[zz(qq - 1) * zz(ql) + row] = 1 / (mid[ui[zz(qq - 1)]] - mid[ref]);
      rl[row] = -1 / (mid[ui[0]] - mid[ref]);
    }
    const Val t2 = awald(c, ce, v, rl, ql);
    tn.push_back("linear");
    ts.push_back(nv(col(t2, "stat").a[0]));
    td.push_back(ql);
    tp.push_back(nv(col(t2, "p").a[0]));
  }
  Frame tt;
  tt.str("test", tn);
  tt.num("stat", ts);
  tt.num("df", td);
  tt.num("p", tp);
  out = vobj();
  put(out, "tab", tb.o);
  put(out, "test", tt.o);
  put(out, "fit", r);
  return true;
}

// typefx()
bool typefx(Ctx &c, Val &out) {
  using namespace ax;
  int adj = 0;
  static const std::vector<const char *> adjs = {"holm", "bonferroni", "BH", "none"};
  CK(marg(c, "adj", adjs, adj));
  Val d;
  CK(dat(c, c["data"], "data", d));
  Cols ys;
  CK(acol(c, d, "ys", {}, false, ys));
  const Val &tv = c["term"];
  if (tv.t != Val::STR && !(tv.t == Val::ARR && tv.a.size() == 1 && tv.a[0].t == Val::STR)) {
    return ckstop(c, "`term` must name one coefficient");
  }
  const std::string term = tv.t == Val::STR ? tv.s : tv.a[0].s;
  int e = 0;
  CK(aest(c, e));
  const Val a = dots(c, {"data", "ys", "term", "est", "adj"});
  std::vector<double> b, se, lo, hi, p, nn;
  for (const std::string &y : ys) {
    Val r;
    CK(ecall(c, e, eargs0(d, {{"y", vstr(y)}}, a), r));
    std::size_t i = 0;
    CK(arow(c, r, term, "the estimator returned no coefficient for `%s`", i));
    b.push_back(cell(r, "est", i));
    se.push_back(cell(r, "se", i));
    lo.push_back(cell(r, "lo", i));
    hi.push_back(cell(r, "hi", i));
    p.push_back(cell(r, "p", i));
    nn.push_back(nv(col(r, "n")));
  }
  // p.adjust()
  const std::size_t m = p.size();
  std::vector<double> pa = p;
  std::vector<std::size_t> o(m);
  for (std::size_t i = 0; i < m; ++i) o[i] = i;
  const double dm = static_cast<double>(m);
  if (adj == 0) {
    std::stable_sort(o.begin(), o.end(), [&](std::size_t x, std::size_t y) { return p[x] < p[y]; });
    double mx = -HUGE_VAL;
    for (std::size_t i = 0; i < m; ++i) {
      mx = std::max(mx, (dm - static_cast<double>(i)) * p[o[i]]);
      pa[o[i]] = std::min(1.0, mx);
    }
  } else if (adj == 1) {
    for (double &x : pa) x = std::min(1.0, dm * x);
  } else if (adj == 2) {
    std::stable_sort(o.begin(), o.end(), [&](std::size_t x, std::size_t y) { return p[x] > p[y]; });
    double mn = HUGE_VAL;
    for (std::size_t i = 0; i < m; ++i) {
      mn = std::min(mn, dm / (dm - static_cast<double>(i)) * p[o[i]]);
      pa[o[i]] = std::min(1.0, mn);
    }
  }
  Frame f;
  f.str("y", ys);
  f.num("est", b);
  f.num("se", se);
  f.num("lo", lo);
  f.num("hi", hi);
  f.num("p", p);
  f.num("padj", pa);
  f.num("n", nn);
  out = vobj();
  put(out, "tab", f.o);
  put(out, "adj", vstr(adjs[zz(adj)]));
  return true;
}

// plcday()
bool plcday(Ctx &c, Val &out) {
  using namespace ax;
  int sh = 0;
  CK(marg(c, "shift", {"both", "x", "z"}, sh));
  Val d;
  CK(dat(c, c["data"], "data", d));
  // .axpdck()
  Cols var, unit, day, with;
  CK(acol(c, d, "var", {1}, false, var));
  CK(acol(c, d, "unit", {1}, false, unit));
  CK(acol(c, d, "day", {1}, false, day));
  CK(acol(c, d, "with", {}, true, with));
  if (sh == 2 && with.empty()) return ckstop(c, "`shift = \"z\"` needs the instruments in `with`");
  int e = 0;
  CK(aest(c, e));
  std::vector<double> sf = {-365, -31, -7, 7, 31, 365};
  if (c.has("shifts")) {
    const Val &o = c["shifts"];
    sf = nums(o);
    bool ok = numv(o) && !sf.empty();
    for (double v : sf) ok = ok && std::isfinite(v);
    for (std::size_t i = 0; ok && i < sf.size(); ++i) {
      ok = sf[i] == std::nearbyint(sf[i]) && sf[i] != 0;
      for (std::size_t j = 0; ok && j < i; ++j) ok = sf[j] != sf[i];
    }
    if (!ok) return ckstop(c, "`shifts` must be distinct whole numbers of days other than 0");
  }
  if (!isnumc(*fc(d, day[0]))) return ckstop(c, "column `%s` must hold numbers or dates", day[0].c_str());
  const std::vector<double> tv = colv(*fc(d, day[0]));
  const std::vector<int> uc = acode(*fc(d, unit[0])).c;
  Cols mv = var;
  if (sh == 0) mv.insert(mv.end(), with.begin(), with.end());
  if (sh == 2) mv = with;
  const Val a = dots(c, {"data", "var", "unit", "day", "shifts", "est", "with", "shift"});
  std::vector<double> s0 = {0}, es, ss, ls, hs, ps, ns;
  s0.insert(s0.end(), sf.begin(), sf.end());
  for (double s : s0) {
    std::vector<int> ix;
    CK(ashift(c, uc, tv, s, ix));
    Ix kr;
    std::vector<int> ik;
    for (std::size_t i = 0; i < ix.size(); ++i) {
      if (ix[i] >= 0) {
        kr.push_back(i);
        ik.push_back(ix[i]);
      }
    }
    if (kr.size() < 2) return ckstop(c, "a shift of %s days leaves fewer than two rows", rfmt(s).c_str());
    Val ds = rows(d, kr);
    for (const std::string &v : mv) setc(ds, v, pick(*fc(d, v), ik));
    Val r;
    CK(ecall(c, e, eargs0(ds, {}, a), r));
    std::size_t i = 0;
    CK(arow(c, r, var[0], "`est` returned no coefficient for `%s`", i));
    es.push_back(cell(r, "est", i));
    ss.push_back(cell(r, "se", i));
    ls.push_back(cell(r, "lo", i));
    hs.push_back(cell(r, "hi", i));
    ps.push_back(cell(r, "p", i));
    ns.push_back(static_cast<double>(kr.size()));
  }
  Frame f;
  f.num("shift", s0);
  f.num("est", es);
  f.num("se", ss);
  f.num("lo", ls);
  f.num("hi", hs);
  f.num("p", ps);
  f.num("n", ns);
  out = vobj();
  put(out, "tab", f.o);
  put(out, "var", vstr(var[0]));
  put(out, "shift", vstr(sh == 0 ? "both" : sh == 1 ? "x" : "z"));
  return true;
}
