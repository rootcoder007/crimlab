// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of crimlab's time, contagion and operations methods:
// the Hawkes process (hwkfit, hwkres, hwksim) and series (trndrr, stpcp,
// cusrv, xcser). Each checks its arguments,
// builds the kernel's arrays, calls the kernels and builds the R result.
// Nothing calls back into the R evaluator.
#include <algorithm>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>

#include "crimlab.h"
#include "kutil.h"
#include "rchk.h"
#include "spin.h"

// R's maths library (declared as bstat.cpp does).
extern "C" double Rf_pnorm5(double, double, double, int, int);
extern "C" double Rf_pt(double, double, int, int);
extern "C" double Rf_qt(double, double, int, int);
extern "C" double Rf_pchisq(double, double, int, int);

using clab::fail;

namespace {

// Protects what it is given and unprotects it all when it goes.
struct prot {
  int n = 0;
  SEXP operator()(SEXP x) {
    ++n;
    return PROTECT(x);
  }
  ~prot() { UNPROTECT(n); }
};

// A named list of the given elements.
SEXP nlist(prot &p, std::initializer_list<std::pair<const char *, SEXP>> e) {
  SEXP v = p(Rf_allocVector(VECSXP, static_cast<R_xlen_t>(e.size())));
  std::vector<const char *> nm;
  R_xlen_t i = 0;
  for (const auto &x : e) {
    SET_VECTOR_ELT(v, i++, x.second);
    nm.push_back(x.first);
  }
  return clab::nlist(v, nm.data());
}

SEXP dvec(prot &p, const double *v, int n) {
  SEXP x = p(Rf_allocVector(REALSXP, n));
  for (int i = 0; i < n; ++i) REAL(x)[i] = std::isnan(v[i]) ? NA_REAL : v[i];
  return x;
}

SEXP strs(prot &p, const std::vector<std::string> &s) {
  SEXP x = p(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(s.size())));
  for (std::size_t i = 0; i < s.size(); ++i) SET_STRING_ELT(x, static_cast<R_xlen_t>(i), Rf_mkChar(s[i].c_str()));
  return x;
}

void setnames(prot &p, SEXP x, const std::vector<std::string> &s) {
  Rf_setAttrib(x, R_NamesSymbol, strs(p, s));
}

SEXP frame(prot &p, std::initializer_list<std::pair<const char *, SEXP>> e) {
  SEXP v = p(Rf_allocVector(VECSXP, static_cast<R_xlen_t>(e.size())));
  std::vector<const char *> nm;
  R_xlen_t i = 0;
  for (const auto &x : e) {
    SET_VECTOR_ELT(v, i++, x.second);
    nm.push_back(x.first);
  }
  return clab::frame(v, nm.data());
}

// Each element of a character vector among the options, as positions.
std::vector<int> picks(SEXP s, const char *nm, const char *const *opt, int n) {
  if (TYPEOF(s) != STRSXP || XLENGTH(s) < 1) {
    s = R_NilValue;  // refused below by pick()
  }
  std::vector<int> out;
  for (R_xlen_t i = 0; i < (Rf_isNull(s) ? 1 : XLENGTH(s)); ++i) {
    SEXP e = PROTECT(Rf_isNull(s) ? R_NilValue : Rf_ScalarString(STRING_ELT(s, i)));
    out.push_back(csp::pick(e, nm, opt, n));
    UNPROTECT(1);
  }
  return out;
}

double pos1(SEXP x, const char *nm) { return REAL(clab::num(x, nm, 1, 0.0, true))[0]; }

// The values of a numeric matrix: finite, at least lo (protected by the
// caller).
SEXP numat(SEXP x, const char *nm, double lo) {
  clab::unlab(x, nm);
  if (!clab::isnum(x)) fail("`%s` must hold numbers.", nm);
  SEXP d = clab::dbl(x);
  for (R_xlen_t i = 0; i < XLENGTH(d); ++i) {
    if (!std::isfinite(REAL(d)[i]) || REAL(d)[i] < lo) fail("`%s` must hold finite numbers from %g.", nm, lo);
  }
  return d;
}

// ------------------------------------------------------------ Hawkes

const char *const kinds[] = {"exp", "weibull", "lomax"};
const char *const bases[] = {"const", "cycle"};

int npar(int kind, int base) { return (base ? 4 : 1) + 1 + (kind ? 2 : 1); }

std::vector<std::string> pnames(int kind, int base) {
  std::vector<std::string> s;
  if (base) s = {"a0", "a1", "a2", "a3"};
  else s = {"a0"};
  s.emplace_back("eta");
  if (kind) {
    s.emplace_back("shape");
    s.emplace_back("scale");
  } else {
    s.emplace_back("beta");
  }
  return s;
}

// Event times: finite, sorted, in [0, hor]; tied ones warned about.
std::vector<double> evtimes(SEXP times, double hor, int nmin) {
  SEXP tv = PROTECT(clab::num(times, "times"));
  const int n = clab::ilen(static_cast<double>(XLENGTH(tv)), "times");
  std::vector<double> t(REAL(tv), REAL(tv) + n);
  UNPROTECT(1);
  if (n < nmin) fail("`times` must hold at least %d events, not %d.", nmin, n);
  for (int i = 0; i < n; ++i) {
    if (t[static_cast<std::size_t>(i)] < 0 || t[static_cast<std::size_t>(i)] > hor) {
      fail("`times` must lie between 0 and `hor`.");
    }
    if (i > 0 && t[static_cast<std::size_t>(i)] < t[static_cast<std::size_t>(i - 1)]) {
      fail("`times` must be sorted from earliest to latest.");
    }
  }
  return t;
}

void tiewarn(const std::vector<double> &t, const char *fn) {
  for (std::size_t i = 1; i < t.size(); ++i) {
    if (t[i] == t[i - 1]) {
      Rf_warningcall(R_NilValue,
                     "%s: some times are tied; events at one instant do not excite each "
                     "other. Times recorded to a resolution (days, say) are better spread "
                     "across it first (`res`).",
                     fn);
      return;
    }
  }
}

// The parameters of one model: n of them, finite.
std::vector<double> parv(SEXP par, const char *nm, int kind, int base) {
  const int np = npar(kind, base);
  SEXP v = PROTECT(clab::num(par, nm, np, R_NegInf));
  std::vector<double> th(REAL(v), REAL(v) + np);
  UNPROTECT(1);
  if (!(th[static_cast<std::size_t>(np - (kind ? 3 : 2))] >= 0)) fail("`%s`: eta must be at least 0.", nm);
  for (int j = np - (kind ? 2 : 1); j < np; ++j) {
    if (!(th[static_cast<std::size_t>(j)] > 0)) fail("`%s`: the kernel parameters must be positive.", nm);
  }
  return th;
}

}  // namespace

extern "C" {

//' @srrstats {TS1.2} hwkfit(), hwkres() and hwksim() take event times as a
//'   plain numeric vector in one unit, sorted, inside the window [0, hor].
//' @noRd
SEXP C_clab_r_hwkfit(SEXP times, SEXP hor, SEXP kernel, SEXP base, SEXP per, SEXP start,
                     SEXP res, SEXP seed, SEXP maxit) {
  clab::who = "hwkfit";
  prot p;
  const double h = pos1(hor, "hor");
  const double pr = pos1(per, "per");
  const int mx = clab::whole(maxit, "maxit", 1);
  std::vector<int> ks = picks(kernel, "kernel", kinds, 3), bs = picks(base, "base", bases, 2);
  std::vector<double> t;
  bool jit = false;
  if (!Rf_isNull(res)) {
    // spread each time across its interval [t, t + res), cut at hor
    const double r = pos1(res, "res");
    SEXP tv = p(clab::num(times, "times"));
    const int n = clab::ilen(static_cast<double>(XLENGTH(tv)), "times");
    t.assign(REAL(tv), REAL(tv) + n);
    kutil::rng g(csp::seed(seed));
    for (double &e : t) {
      if (e < 0 || e > h) fail("`times` must lie between 0 and `hor`.");
      e += g.unif() * std::fmin(r, h - e);
    }
    kutil::hsort(t.data(), n, [](double a, double b) { return a < b; });
    SEXP tj = p(Rf_allocVector(REALSXP, n));
    std::copy(t.begin(), t.end(), REAL(tj));
    t = evtimes(tj, h, 10);
    jit = true;
  } else {
    t = evtimes(times, h, 10);
    tiewarn(t, "hwkfit");
  }
  const int n = static_cast<int>(t.size());
  const int nc = static_cast<int>(ks.size() * bs.size());
  std::vector<double> st0;
  if (!Rf_isNull(start)) {
    if (nc > 1) fail("`start` needs one kernel and one baseline.");
    st0 = parv(start, "start", ks[0], bs[0]);
    if (!(st0[static_cast<std::size_t>(npar(ks[0], bs[0]) - (ks[0] ? 3 : 2))] < 1)) {
      fail("`start`: eta must be below 1.");
    }
  }
  // every combination; the one with the least AIC is the fit
  SEXP ck = p(Rf_allocVector(STRSXP, nc)), cb = p(Rf_allocVector(STRSXP, nc));
  SEXP cl = p(Rf_allocVector(REALSXP, nc)), ckp = p(Rf_allocVector(INTSXP, nc)),
       ca = p(Rf_allocVector(REALSXP, nc)), cbi = p(Rf_allocVector(REALSXP, nc)),
       ce = p(Rf_allocVector(REALSXP, nc)), cv = p(Rf_allocVector(LGLSXP, nc));
  int bi = -1;
  std::vector<double> bth, bcov, bst;
  int bk = 0, bb = 0;
  for (int c = 0; c < nc; ++c) {
    const int kd = ks[static_cast<std::size_t>(c) / bs.size()], bd = bs[static_cast<std::size_t>(c) % bs.size()];
    const int np = npar(kd, bd);
    std::vector<double> th(static_cast<std::size_t>(np)), cov(static_cast<std::size_t>(np * np)), st(3);
    const int rc = clab_hwkfit(t.data(), n, h, kd, bd, pr, st0.empty() ? nullptr : st0.data(), mx,
                               th.data(), cov.data(), st.data());
    if (rc != CLAB_OK) fail("the %s fit failed: the likelihood is not finite at any start.", kinds[kd]);
    SET_STRING_ELT(ck, c, Rf_mkChar(kinds[kd]));
    SET_STRING_ELT(cb, c, Rf_mkChar(bases[bd]));
    REAL(cl)[c] = st[0];
    INTEGER(ckp)[c] = np;
    REAL(ca)[c] = 2 * np - 2 * st[0];
    REAL(cbi)[c] = np * std::log(static_cast<double>(n)) - 2 * st[0];
    REAL(ce)[c] = th[static_cast<std::size_t>(np - (kd ? 3 : 2))];
    LOGICAL(cv)[c] = st[2] > 0;
    if (bi < 0 || REAL(ca)[c] < REAL(ca)[bi]) {
      bi = c;
      bth = th;
      bcov = cov;
      bst = st;
      bk = kd;
      bb = bd;
    }
  }
  const int np = npar(bk, bb);
  std::vector<std::string> pn = pnames(bk, bb);
  SEXP par = dvec(p, bth.data(), np);
  setnames(p, par, pn);
  std::vector<double> se(static_cast<std::size_t>(np));
  for (int j = 0; j < np; ++j) se[static_cast<std::size_t>(j)] = std::sqrt(bcov[static_cast<std::size_t>(j * np + j)]);
  SEXP sev = dvec(p, se.data(), np);
  setnames(p, sev, pn);
  SEXP cov = p(Rf_allocMatrix(REALSXP, np, np));
  for (int j = 0; j < np * np; ++j) REAL(cov)[j] = std::isnan(bcov[static_cast<std::size_t>(j)]) ? NA_REAL : bcov[static_cast<std::size_t>(j)];
  SEXP dn = p(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dn, 0, strs(p, pn));
  SET_VECTOR_ELT(dn, 1, strs(p, pn));
  Rf_setAttrib(cov, R_DimNamesSymbol, dn);
  // the branching arithmetic: offspring per event, cluster size, the
  // stationary rate (constant baseline) and the share triggered
  const double eta = bth[static_cast<std::size_t>(np - (bk ? 3 : 2))];
  const double br[4] = {eta, eta < 1 ? 1 / (1 - eta) : R_PosInf,
                        bb == 0 && eta < 1 ? std::exp(bth[0]) / (1 - eta) : NA_REAL, eta};
  static const char *const bn[] = {"eta", "size", "rate", "share"};
  SEXP branch = p(clab::nvec(4, bn, br));
  if (eta > 0.99) {
    Rf_warningcall(R_NilValue,
                   "hwkfit: the branching ratio is %.3f, near 1: the fit is near-critical "
                   "and probably mistakes a trend or a cycle for contagion (try base = "
                   "\"cycle\").",
                   eta);
  }
  SEXP tt = dvec(p, t.data(), n);
  SEXP cmp = frame(p, {{"kernel", ck}, {"base", cb}, {"loglik", cl}, {"k", ckp}, {"aic", ca},
                       {"bic", cbi}, {"eta", ce}, {"conv", cv}});
  SEXP out = nlist(p, {{"par", par},
                       {"se", sev},
                       {"cov", cov},
                       {"loglik", p(Rf_ScalarReal(bst[0]))},
                       {"aic", p(Rf_ScalarReal(REAL(ca)[bi]))},
                       {"bic", p(Rf_ScalarReal(REAL(cbi)[bi]))},
                       {"branch", branch},
                       {"conv", p(Rf_ScalarLogical(bst[2] > 0))},
                       {"iter", p(Rf_ScalarInteger(static_cast<int>(bst[1])))},
                       {"kernel", p(Rf_mkString(kinds[bk]))},
                       {"base", p(Rf_mkString(bases[bb]))},
                       {"per", p(Rf_ScalarReal(pr))},
                       {"n", p(Rf_ScalarInteger(n))},
                       {"hor", p(Rf_ScalarReal(h))},
                       {"jit", p(Rf_ScalarLogical(jit))},
                       {"times", tt},
                       {"cmp", cmp}});
  Rf_classgets(out, p(Rf_mkString("hwkfit")));
  return out;
}

SEXP C_clab_r_hwkres(SEXP times, SEXP hor, SEXP par, SEXP kernel, SEXP base, SEXP per) {
  clab::who = "hwkres";
  prot p;
  const double h = pos1(hor, "hor"), pr = pos1(per, "per");
  const int kd = csp::pick(kernel, "kernel", kinds, 3), bd = csp::pick(base, "base", bases, 2);
  std::vector<double> t = evtimes(times, h, 1);
  std::vector<double> th = parv(par, "par", kd, bd);
  const int n = static_cast<int>(t.size());
  std::vector<double> tau(static_cast<std::size_t>(n)), work(static_cast<std::size_t>(n < 100 ? 3 * (2 * n + 1) * (2 * n + 1) + n : n));
  double ks[2];
  // the arguments are checked: the kernel cannot refuse them
  clab_hwkres(t.data(), n, h, kd, bd, pr, th.data(), tau.data(), ks, work.data());
  std::vector<double> u(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) u[static_cast<std::size_t>(i)] = -std::expm1(-tau[static_cast<std::size_t>(i)]);
  SEXP out = nlist(p, {{"tau", dvec(p, tau.data(), n)},
                       {"u", dvec(p, u.data(), n)},
                       {"d", p(Rf_ScalarReal(ks[0]))},
                       {"p", p(Rf_ScalarReal(ks[1]))},
                       {"n", p(Rf_ScalarInteger(n))}});
  Rf_classgets(out, p(Rf_mkString("hwkres")));
  return out;
}

SEXP C_clab_r_hwksim(SEXP par, SEXP hor, SEXP kernel, SEXP base, SEXP per, SEXP win,
                     SEXP sigma, SEXP seed, SEXP nmax) {
  clab::who = "hwksim";
  prot p;
  const double h = pos1(hor, "hor"), pr = pos1(per, "per");
  const int kd = csp::pick(kernel, "kernel", kinds, 3), bd = csp::pick(base, "base", bases, 2);
  std::vector<double> th = parv(par, "par", kd, bd);
  const int cap = clab::whole(nmax, "nmax", 1);
  double w[4], sg = 0;
  const bool sp = !Rf_isNull(win);
  if (sp) {
    csp::win(win, w);
    sg = pos1(sigma, "sigma");
    // the baseline is per unit area: the window's rate is a0 + log(area)
    th[0] += std::log((w[1] - w[0]) * (w[3] - w[2]));
  }
  const std::size_t sc = static_cast<std::size_t>(cap);
  std::vector<double> t(sc), x(sp ? sc : 0), y(sp ? sc : 0), work(3 * sc);
  std::vector<int> gen(sc), iwork(2 * sc);
  int n = 0;
  const int rc = clab_hwksim(h, kd, bd, pr, th.data(), sp ? w : nullptr, sg, csp::seed(seed), cap,
                             t.data(), sp ? x.data() : nullptr, sp ? y.data() : nullptr, gen.data(),
                             &n, work.data(), iwork.data());
  if (rc != CLAB_OK) {
    fail("more than `nmax` = %d events: raise it, or lower eta (at 1 or above the process "
         "explodes).",
         cap);
  }
  SEXP gv = p(Rf_allocVector(INTSXP, n));
  std::copy(gen.begin(), gen.begin() + n, INTEGER(gv));
  if (!sp) return frame(p, {{"t", dvec(p, t.data(), n)}, {"gen", gv}});
  return frame(p, {{"t", dvec(p, t.data(), n)},
                   {"x", dvec(p, x.data(), n)},
                   {"y", dvec(p, y.data(), n)},
                   {"gen", gv}});
}

}  // extern "C"

namespace {

// ------------------------------------------------------------ series

const char *const alts[] = {"both", "up", "down"};

// Periods: NULL for 1..n, else n finite numbers, strictly increasing.
std::vector<double> periods(SEXP x, int n) {
  std::vector<double> v(static_cast<std::size_t>(n));
  if (Rf_isNull(x)) {
    for (int i = 0; i < n; ++i) v[static_cast<std::size_t>(i)] = i + 1;
    return v;
  }
  SEXP xv = PROTECT(clab::num(x, "x", n));
  for (int i = 0; i < n; ++i) {
    v[static_cast<std::size_t>(i)] = REAL(xv)[i];
    if (i > 0 && !(v[static_cast<std::size_t>(i)] > v[static_cast<std::size_t>(i - 1)])) {
      fail("`x` must be strictly increasing.");
    }
  }
  UNPROTECT(1);
  return v;
}

double zlev(SEXP level) { return Rf_qnorm5(0.5 + 0.5 * clab::one(level, "level", 0.0, 1.0), 0.0, 1.0, 1, 0); }

SEXP lgl1(prot &p, double v) { return p(Rf_ScalarLogical(v != 0)); }

}  // namespace

extern "C" {

SEXP C_clab_r_trndrr(SEXP y, SEXP x, SEXP expo, SEXP seas, SEXP alt, SEXP level) {
  clab::who = "trndrr";
  prot p;
  SEXP yv = p(clab::num(y, "y", -1, 0.0));
  const int n = clab::ilen(static_cast<double>(XLENGTH(yv)), "y");
  if (n < 3) fail("`y` must hold at least 3 periods, not %d.", n);
  const double *yy = REAL(yv);
  std::vector<double> xv = periods(x, n);
  const int a = csp::pick(alt, "alt", alts, 3);
  const double z = zlev(level), lv = REAL(clab::dbl(level))[0];
  std::vector<double> off(static_cast<std::size_t>(n), 0.0);
  if (!Rf_isNull(expo)) {
    SEXP e = p(clab::num(expo, "expo", n, 0.0, true));
    for (int i = 0; i < n; ++i) off[static_cast<std::size_t>(i)] = std::log(REAL(e)[i]);
  }
  // the Poisson trend: log mu = b0 + b1 (x - mean x) + log expo
  double xm = 0, tot = 0;
  for (int i = 0; i < n; ++i) {
    xm += xv[static_cast<std::size_t>(i)] / n;
    tot += yy[i];
  }
  if (!(tot > 0)) fail("`y` has no events: every count is 0.");
  int lo = n, hi = -1;
  for (int i = 0; i < n; ++i) {
    if (yy[i] > 0) {
      lo = std::min(lo, i);
      hi = std::max(hi, i);
    }
  }
  std::vector<double> d(static_cast<std::size_t>(2 * n)), mu(static_cast<std::size_t>(n), NA_REAL);
  for (int i = 0; i < n; ++i) {
    d[static_cast<std::size_t>(i)] = 1;
    d[static_cast<std::size_t>(n + i)] = xv[static_cast<std::size_t>(i)] - xm;
  }
  double rr[8];
  for (double &v : rr) v = NA_REAL;
  rr[7] = 0;
  if (lo == hi && (lo == 0 || hi == n - 1)) {
    // all the events in the first or the last period: no finite trend
    rr[0] = lo == 0 ? 0.0 : R_PosInf;
    Rf_warningcall(R_NilValue,
                   "trndrr: every event falls in the %s period: the trend has no finite "
                   "estimate.",
                   lo == 0 ? "first" : "last");
  } else {
    double beta[2], vb[4], vc[4], st[8];
    int al[2], iw[2];
    std::vector<double> w(static_cast<std::size_t>(2 * n + 6 * n + 2 + 4));
    const int rc = clab_cntreg(yy, d.data(), n, 2, off.data(), nullptr, 0, 0, 25, 1e-8, beta, vb, vc, st, al, w.data(), iw);
    if (rc != CLAB_OK) fail("the Poisson fit failed: the deviance left the finite range.");
    double pe = 0;
    for (int i = 0; i < n; ++i) {
      const double m = std::exp(beta[0] + beta[1] * d[static_cast<std::size_t>(n + i)] + off[static_cast<std::size_t>(i)]);
      mu[static_cast<std::size_t>(i)] = m;
      pe += (yy[i] - m) * (yy[i] - m) / m;
    }
    // the quasi-Poisson interval when the dispersion is above 1.5
    const double disp = pe / (n - 2), b = beta[1];
    double se = std::sqrt(vb[3]), q = z, pv;
    const bool over = disp > 1.5;
    if (over) {
      se *= std::sqrt(disp);
      q = Rf_qt(0.5 + 0.5 * lv, n - 2, 1, 0);
      pv = 2 * Rf_pt(-std::fabs(b / se), n - 2, 1, 0);
    } else {
      pv = 2 * Rf_pnorm5(-std::fabs(b / se), 0.0, 1.0, 1, 0);
    }
    const double r0[8] = {std::exp(b), std::exp(b - q * se), std::exp(b + q * se), pv, b, se, disp,
                          over ? 1.0 : 0.0};
    for (int j = 0; j < 8; ++j) rr[j] = r0[j];
  }
  SEXP rrf = frame(p, {{"rr", dvec(p, rr, 1)},
                       {"lo", dvec(p, rr + 1, 1)},
                       {"hi", dvec(p, rr + 2, 1)},
                       {"p", dvec(p, rr + 3, 1)},
                       {"b", dvec(p, rr + 4, 1)},
                       {"se", dvec(p, rr + 5, 1)},
                       {"disp", dvec(p, rr + 6, 1)},
                       {"over", lgl1(p, rr[7])}});
  // Mann-Kendall and Theil-Sen
  double st[10];
  std::vector<double> w(static_cast<std::size_t>(n) * static_cast<std::size_t>(n + 1) / 2);
  std::vector<int> iw(static_cast<std::size_t>(n));
  clab_trndrr(xv.data(), yy, n, a, z, st, w.data(), iw.data());
  SEXP mk = frame(p, {{"s", dvec(p, st, 1)},
                      {"var", dvec(p, st + 1, 1)},
                      {"tau", dvec(p, st + 2, 1)},
                      {"p", dvec(p, st + 3, 1)},
                      {"exact", lgl1(p, n <= 8)}});
  SEXP sen = frame(p, {{"slope", dvec(p, st + 4, 1)},
                       {"int", dvec(p, st + 5, 1)},
                       {"lo", dvec(p, st + 6, 1)},
                       {"hi", dvec(p, st + 7, 1)},
                       {"clamp", lgl1(p, st[8] + st[9])}});
  // period on period
  std::vector<double> chg(static_cast<std::size_t>(n), NA_REAL), pct(static_cast<std::size_t>(n), NA_REAL);
  for (int i = 1; i < n; ++i) {
    chg[static_cast<std::size_t>(i)] = yy[i] - yy[i - 1];
    pct[static_cast<std::size_t>(i)] = yy[i - 1] > 0 ? 100 * (yy[i] - yy[i - 1]) / yy[i - 1] : NA_REAL;
  }
  SEXP yoy = frame(p, {{"x", dvec(p, xv.data(), n)},
                       {"y", dvec(p, yy, n)},
                       {"chg", dvec(p, chg.data(), n)},
                       {"pct", dvec(p, pct.data(), n)}});
  // seasons: the chi-square test of equal totals
  SEXP sv = R_NilValue;
  if (!Rf_isNull(seas)) {
    csp::lab l = csp::codes(seas, "seas", n);
    p.n += 2;
    if (l.g < 2) fail("`seas` has one season: the test needs two or more.");
    std::vector<double> nn(static_cast<std::size_t>(l.g), 0.0), ex(static_cast<std::size_t>(l.g), tot / l.g);
    for (int i = 0; i < n; ++i) nn[static_cast<std::size_t>(INTEGER(l.code)[i])] += yy[i];
    double x2 = 0;
    for (int g = 0; g < l.g; ++g) x2 += (nn[static_cast<std::size_t>(g)] - ex[static_cast<std::size_t>(g)]) * (nn[static_cast<std::size_t>(g)] - ex[static_cast<std::size_t>(g)]) / ex[static_cast<std::size_t>(g)];
    SEXP tab = frame(p, {{"season", l.lev}, {"n", dvec(p, nn.data(), l.g)}, {"expd", dvec(p, ex.data(), l.g)}});
    sv = nlist(p, {{"tab", tab},
                   {"stat", p(Rf_ScalarReal(x2))},
                   {"df", p(Rf_ScalarInteger(l.g - 1))},
                   {"p", p(Rf_ScalarReal(Rf_pchisq(x2, l.g - 1, 0, 0)))}});
  }
  SEXP out = nlist(p, {{"rr", rrf},
                       {"mk", mk},
                       {"sen", sen},
                       {"yoy", yoy},
                       {"seas", sv},
                       {"fitted", dvec(p, mu.data(), n)},
                       {"n", p(Rf_ScalarInteger(n))},
                       {"alt", p(Rf_mkString(alts[a]))},
                       {"level", p(Rf_ScalarReal(lv))}});
  Rf_classgets(out, p(Rf_mkString("trndrr")));
  return out;
}

SEXP C_clab_r_stpcp(SEXP y, SEXP minseg, SEXP nperm, SEXP seed) {
  clab::who = "stpcp";
  prot p;
  SEXP yv = p(clab::num(y, "y"));
  const int n = clab::ilen(static_cast<double>(XLENGTH(yv)), "y");
  const int m = clab::whole(minseg, "minseg", 1), np = clab::whole(nperm, "nperm", 1);
  if (n < 2 * m) fail("`y` must hold at least 2 `minseg` = %d values, not %d.", 2 * m, n);
  double st[12];
  std::vector<double> w(static_cast<std::size_t>(n));
  std::vector<int> iw(static_cast<std::size_t>(n));
  clab_stpcp(REAL(yv), n, m, np, csp::seed(seed), st, w.data(), iw.data());
  SEXP ts = p(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(ts, 0, Rf_mkChar("step"));
  SET_STRING_ELT(ts, 1, Rf_mkChar("pettitt"));
  SEXP k = p(Rf_allocVector(INTSXP, 2));
  INTEGER(k)[0] = static_cast<int>(st[0]);
  INTEGER(k)[1] = static_cast<int>(st[7]);
  const double be[2] = {st[1], st[10]}, af[2] = {st[2], st[11]}, sa[2] = {st[3], st[8]},
               pv[2] = {st[4], st[9]};
  SEXP tab = frame(p, {{"test", ts},
                       {"k", k},
                       {"before", dvec(p, be, 2)},
                       {"after", dvec(p, af, 2)},
                       {"stat", dvec(p, sa, 2)},
                       {"p", dvec(p, pv, 2)}});
  SEXP out = nlist(p, {{"tab", tab},
                       {"nperm", p(Rf_ScalarReal(st[5]))},
                       {"exact", lgl1(p, st[6])},
                       {"n", p(Rf_ScalarInteger(n))}});
  Rf_classgets(out, p(Rf_mkString("stpcp")));
  return out;
}

SEXP C_clab_r_cusrv(SEXP y, SEXP expd, SEXP start, SEXP k, SEXP h, SEXP trans, SEXP reset) {
  clab::who = "cusrv";
  prot p;
  static const char *const tr[] = {"std", "rossi", "anscombe", "none"};
  SEXP yv = p(clab::num(y, "y", -1, 0.0));
  const int n = clab::ilen(static_cast<double>(XLENGTH(yv)), "y");
  const double kk = REAL(clab::num(k, "k", 1))[0], hh = pos1(h, "h");
  const int t = csp::pick(trans, "trans", tr, 4);
  const bool rs = csp::flag(reset, "reset");
  int s0 = 0;
  std::vector<double> m(static_cast<std::size_t>(n));
  if (Rf_isNull(expd)) {
    // the mean of the counts before start
    if (Rf_isNull(start)) fail("give `expd`, or `start`, the first period watched.");
    s0 = clab::whole(start, "start", 2) - 1;
    if (s0 >= n) fail("`start` must be at most the number of counts, %d.", n);
    double b = 0;
    for (int i = 0; i < s0; ++i) b += REAL(yv)[i] / s0;
    if (!(b > 0)) fail("the counts before `start` are all 0: no expected count.");
    for (double &v : m) v = b;
  } else {
    SEXP e = p(clab::num(expd, "expd", n, 0.0, true));
    std::copy(REAL(e), REAL(e) + n, m.begin());
  }
  const int nw = n - s0;
  std::vector<double> z(static_cast<std::size_t>(nw)), s(static_cast<std::size_t>(nw));
  SEXP al = p(Rf_allocVector(LGLSXP, nw));
  clab_cusrv(REAL(yv) + s0, m.data() + s0, nw, kk, hh, t, rs ? 1 : 0, z.data(), s.data(), LOGICAL(al));
  SEXP ix = p(Rf_allocVector(INTSXP, nw));
  for (int i = 0; i < nw; ++i) INTEGER(ix)[i] = s0 + i + 1;
  return frame(p, {{"i", ix},
                   {"y", dvec(p, REAL(yv) + s0, nw)},
                   {"expd", dvec(p, m.data() + s0, nw)},
                   {"z", dvec(p, z.data(), nw)},
                   {"s", dvec(p, s.data(), nw)},
                   {"alarm", al}});
}

SEXP C_clab_r_xcser(SEXP x, SEXP y, SEXP lagmax, SEXP dif, SEXP level) {
  clab::who = "xcser";
  prot p;
  SEXP xv = p(clab::num(x, "x")), yv = p(clab::num(y, "y"));
  SEXP xn = Rf_getAttrib(x, R_NamesSymbol), yn = Rf_getAttrib(y, R_NamesSymbol);
  std::vector<double> a, b;
  SEXP common = R_NilValue;
  if (!Rf_isNull(xn) && !Rf_isNull(yn)) {
    // the periods both series have, in the order of x
    std::vector<int> ix;
    for (R_xlen_t i = 0; i < XLENGTH(xv); ++i) {
      for (R_xlen_t j = 0; j < XLENGTH(yv); ++j) {
        if (clab::seql(STRING_ELT(xn, i), STRING_ELT(yn, j))) {
          a.push_back(REAL(xv)[i]);
          b.push_back(REAL(yv)[j]);
          ix.push_back(static_cast<int>(i));
          break;
        }
      }
    }
    common = p(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(ix.size())));
    for (std::size_t i = 0; i < ix.size(); ++i) SET_STRING_ELT(common, static_cast<R_xlen_t>(i), STRING_ELT(xn, ix[i]));
  } else {
    if (XLENGTH(xv) != XLENGTH(yv)) fail("`x` and `y` must have one length, or names to match them by.");
    a.assign(REAL(xv), REAL(xv) + XLENGTH(xv));
    b.assign(REAL(yv), REAL(yv) + XLENGTH(yv));
  }
  if (csp::flag(dif, "dif")) {
    for (std::size_t i = a.empty() ? 0 : a.size() - 1; i > 0; --i) {
      a[i] -= a[i - 1];
      b[i] -= b[i - 1];
    }
    if (!a.empty()) {
      a.erase(a.begin());
      b.erase(b.begin());
    }
  }
  const int n = static_cast<int>(a.size());
  if (n < 3) fail("the series have %d periods in common: at least 3 are needed.", n);
  const double z = zlev(level);
  const int lag = Rf_isNull(lagmax) ? std::max(0, std::min(n - 1, static_cast<int>(std::floor(10 * std::log10(n / 2.0)))))
                                    : clab::whole(lagmax, "lagmax", 0);
  if (lag >= n) fail("`lagmax` must be below the number of periods, %d.", n);
  std::vector<double> r(static_cast<std::size_t>(2 * lag + 1));
  if (clab_xcser(a.data(), b.data(), n, lag, r.data()) != CLAB_OK) fail("a series is constant: no correlation.");
  const double r0 = r[static_cast<std::size_t>(lag)];
  const double fz = std::atanh(r0), sz = n > 3 ? 1 / std::sqrt(n - 3.0) : NA_REAL;
  const double tt = r0 * std::sqrt((n - 2) / (1 - r0 * r0));
  const double c[4] = {r0, n > 3 ? std::tanh(fz - z * sz) : NA_REAL, n > 3 ? std::tanh(fz + z * sz) : NA_REAL,
                       std::fabs(r0) < 1 ? 2 * Rf_pt(-std::fabs(tt), n - 2, 1, 0) : 0.0};
  SEXP lg = p(Rf_allocVector(INTSXP, 2 * lag + 1));
  for (int k = -lag; k <= lag; ++k) INTEGER(lg)[lag + k] = k;
  SEXP cor = frame(p, {{"r", dvec(p, c, 1)},
                       {"lo", dvec(p, c + 1, 1)},
                       {"hi", dvec(p, c + 2, 1)},
                       {"p", dvec(p, c + 3, 1)},
                       {"n", p(Rf_ScalarInteger(n))}});
  SEXP out = nlist(p, {{"cor", cor},
                       {"ccf", frame(p, {{"lag", lg}, {"r", dvec(p, r.data(), 2 * lag + 1)}})},
                       {"common", common},
                       {"level", p(Rf_ScalarReal(REAL(clab::dbl(level))[0]))}});
  Rf_classgets(out, p(Rf_mkString("xcser")));
  return out;
}

}  // extern "C"

extern "C" {

SEXP C_clab_r_kmrec(SEXP time, SEXP event, SEXP grp, SEXP unit, SEXP end, SEXP level) {
  clab::who = "kmrec";
  prot p;
  SEXP tv = p(clab::num(time, "time", -1, 0.0));
  const int n0 = clab::ilen(static_cast<double>(XLENGTH(tv)), "time");
  const double z = zlev(level);
  std::vector<double> t;
  std::vector<int> ev, gr;
  SEXP lev;
  int g = 1;
  csp::lab gl{R_NilValue, R_NilValue, 1};
  if (!Rf_isNull(grp)) {
    gl = csp::codes(grp, "grp", n0);
    p.n += 2;
    g = gl.g;
  }
  if (Rf_isNull(unit)) {
    // durations
    t.assign(REAL(tv), REAL(tv) + n0);
    if (Rf_isNull(event)) {
      ev.assign(static_cast<std::size_t>(n0), 1);
    } else {
      SEXP e = p(csp::bin01(event, "event", n0));
      ev.assign(INTEGER(e), INTEGER(e) + n0);
    }
    for (int i = 0; i < n0; ++i) gr.push_back(Rf_isNull(grp) ? 0 : INTEGER(gl.code)[i]);
  } else {
    // recurrence: the gaps between a unit's events, the last one censored
    // at end
    if (!Rf_isNull(event)) fail("`event` is not read with `unit`: the gaps are the events.");
    csp::lab u = csp::codes(unit, "unit", n0);
    p.n += 2;
    double e = R_PosInf;
    if (!Rf_isNull(end)) {
      e = REAL(clab::num(end, "end", 1))[0];
      for (int i = 0; i < n0; ++i) {
        if (REAL(tv)[i] > e) fail("`end` must be at or after every time.");
      }
    }
    std::vector<int> ix(static_cast<std::size_t>(n0));
    for (int i = 0; i < n0; ++i) ix[static_cast<std::size_t>(i)] = i;
    const int *uc = INTEGER(u.code);
    const double *tp = REAL(tv);
    kutil::hsort(ix.data(), n0, [&](int a, int b) {
      return uc[a] < uc[b] || (uc[a] == uc[b] && (tp[a] < tp[b] || (tp[a] == tp[b] && a < b)));
    });
    for (int i = 0; i < n0; ++i) {
      const int a = ix[static_cast<std::size_t>(i)];
      const int gg = Rf_isNull(grp) ? 0 : INTEGER(gl.code)[a];
      if (i > 0 && uc[ix[static_cast<std::size_t>(i - 1)]] == uc[a]) {
        const int b = ix[static_cast<std::size_t>(i - 1)];
        if (!Rf_isNull(grp) && INTEGER(gl.code)[b] != gg) fail("`grp` must be the same for every time of a unit.");
        t.push_back(tp[a] - tp[b]);
        ev.push_back(1);
        gr.push_back(gg);
      }
      const bool last = i == n0 - 1 || uc[ix[static_cast<std::size_t>(i + 1)]] != uc[a];
      if (last && std::isfinite(e)) {
        t.push_back(e - tp[a]);
        ev.push_back(0);
        gr.push_back(gg);
      }
    }
    if (t.empty()) fail("no unit has two times and `end` is not given: there are no gaps.");
  }
  lev = Rf_isNull(grp) ? p(Rf_mkString("all")) : gl.lev;
  const int n = static_cast<int>(t.size());
  const std::size_t sn = static_cast<std::size_t>(n);
  std::vector<int> go(sn), iw(sn);
  std::vector<double> tt(sn), nr(sn), nd(sn), nc(sn), s(sn), se(sn), lo(sn), hi(sn),
      med(2 * static_cast<std::size_t>(g)), w(3 * static_cast<std::size_t>(g) + static_cast<std::size_t>(g) * static_cast<std::size_t>(g));
  double lr[2];
  int m = 0;
  clab_kmrec(t.data(), ev.data(), gr.data(), n, g, z, &m, go.data(), tt.data(), nr.data(), nd.data(),
             nc.data(), s.data(), se.data(), lo.data(), hi.data(), med.data(), lr, w.data(), iw.data());
  SEXP gk = p(Rf_allocVector(STRSXP, m));
  for (int i = 0; i < m; ++i) SET_STRING_ELT(gk, i, STRING_ELT(lev, go[static_cast<std::size_t>(i)]));
  SEXP km = frame(p, {{"grp", gk},
                      {"time", dvec(p, tt.data(), m)},
                      {"nrisk", dvec(p, nr.data(), m)},
                      {"nevent", dvec(p, nd.data(), m)},
                      {"ncens", dvec(p, nc.data(), m)},
                      {"surv", dvec(p, s.data(), m)},
                      {"se", dvec(p, se.data(), m)},
                      {"lo", dvec(p, lo.data(), m)},
                      {"hi", dvec(p, hi.data(), m)}});
  std::vector<double> gn(static_cast<std::size_t>(g), 0.0), ge(static_cast<std::size_t>(g), 0.0);
  for (int i = 0; i < n; ++i) {
    gn[static_cast<std::size_t>(gr[static_cast<std::size_t>(i)])] += 1;
    ge[static_cast<std::size_t>(gr[static_cast<std::size_t>(i)])] += ev[static_cast<std::size_t>(i)];
  }
  SEXP md = frame(p, {{"grp", lev},
                      {"n", dvec(p, gn.data(), g)},
                      {"events", dvec(p, ge.data(), g)},
                      {"median", dvec(p, med.data(), g)}});
  SEXP lrv = R_NilValue;
  if (g > 1) {
    const double pv = std::isnan(lr[0]) ? NA_REAL : Rf_pchisq(lr[0], lr[1], 0, 0);
    lrv = nlist(p, {{"stat", dvec(p, lr, 1)}, {"df", p(Rf_ScalarInteger(g - 1))}, {"p", p(Rf_ScalarReal(pv))}});
  }
  SEXP out = nlist(p, {{"km", km},
                       {"med", md},
                       {"lr", lrv},
                       {"n", p(Rf_ScalarInteger(n))},
                       {"level", p(Rf_ScalarReal(REAL(clab::dbl(level))[0]))}});
  Rf_classgets(out, p(Rf_mkString("kmrec")));
  return out;
}

}  // extern "C"

extern "C" {

SEXP C_clab_r_hypq(SEXP lam, SEXP pref, SEXP mu, SEXP queue, SEXP tt) {
  clab::who = "hypq";
  prot p;
  static const char *const qs[] = {"loss", "wait"};
  SEXP lv = p(clab::num(lam, "lam", -1, 0.0, true));
  const int m = clab::ilen(static_cast<double>(XLENGTH(lv)), "lam");
  const double u = pos1(mu, "mu");
  const int w = csp::pick(queue, "queue", qs, 2);
  // the orders: an m x n matrix of unit numbers, or a list of m vectors
  int n = -1;
  std::vector<int> pr;
  if (TYPEOF(pref) == VECSXP) {
    if (XLENGTH(pref) != m) fail("`pref` must hold one order per atom, %d.", m);
    n = static_cast<int>(XLENGTH(VECTOR_ELT(pref, 0)));
    pr.assign(static_cast<std::size_t>(m * n), 0);
    for (int j = 0; j < m; ++j) {
      SEXP v = VECTOR_ELT(pref, j);
      if (XLENGTH(v) != n) fail("`pref`: every order must list all %d units.", n);
      SEXP d = p(clab::num(v, "pref"));
      for (int r = 0; r < n; ++r) pr[static_cast<std::size_t>(j + r * m)] = static_cast<int>(REAL(d)[r]) - 1;
    }
  } else {
    if (!Rf_isMatrix(pref) || Rf_nrows(pref) != m) fail("`pref` must be a matrix with one row per atom, %d, or a list.", m);
    n = Rf_ncols(pref);
    SEXP d = p(numat(pref, "pref", R_NegInf));
    pr.resize(static_cast<std::size_t>(m * n));
    for (int i = 0; i < m * n; ++i) pr[static_cast<std::size_t>(i)] = static_cast<int>(REAL(d)[i]) - 1;
  }
  for (std::size_t i = 0; i < pr.size(); ++i) {
    if (pr[i] < 0 || pr[i] >= n) fail("`pref` must hold unit numbers from 1 to %d.", n);
  }
  if (n > 10) fail("at most 10 units: the model has 2^n states.");
  const std::size_t s = static_cast<std::size_t>(1) << n;
  std::vector<double> pst(s), wl(static_cast<std::size_t>(n)), dp(static_cast<std::size_t>(n * m)), work(s * s + s);
  double st[3];
  const int rc = clab_hypq(REAL(lv), pr.data(), m, n, u, w, pst.data(), wl.data(), dp.data(), st, work.data());
  if (rc != CLAB_OK) {
    fail(w ? "every order must list each unit once, and the calls must be fewer than the units "
             "can serve (rho below 1) for a queue."
           : "every order must list each unit once.");
  }
  SEXP bs = p(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(s)));
  for (std::size_t b = 0; b < s; ++b) {
    std::string k;
    for (int v = 0; v < n; ++v) k += ((b >> v) & 1u) ? '1' : '0';
    SET_STRING_ELT(bs, static_cast<R_xlen_t>(b), Rf_mkChar(k.c_str()));
  }
  SEXP dm = p(Rf_allocMatrix(REALSXP, n, m));
  std::copy(dp.begin(), dp.end(), REAL(dm));
  double big = 0;
  for (int j = 0; j < m; ++j) big += REAL(lv)[j];
  SEXP trv = R_NilValue;
  if (!Rf_isNull(tt)) {
    if (!Rf_isMatrix(tt) || Rf_nrows(tt) != n || Rf_ncols(tt) != m) {
      fail("`tt` must be a matrix of travel times, one row per unit and one column per atom.");
    }
    SEXP tv = p(numat(tt, "tt", 0.0));
    std::vector<double> at(static_cast<std::size_t>(m));
    double num = 0, den = 0;
    for (int j = 0; j < m; ++j) {
      double a = 0, b = 0;
      for (int v = 0; v < n; ++v) {
        a += dp[static_cast<std::size_t>(v + j * n)] * REAL(tv)[v + j * n];
        b += dp[static_cast<std::size_t>(v + j * n)];
      }
      at[static_cast<std::size_t>(j)] = a / b;
      num += a;
      den += b;
    }
    trv = nlist(p, {{"mean", p(Rf_ScalarReal(num / den))}, {"atom", dvec(p, at.data(), m)}});
  }
  SEXP out = nlist(p, {{"state", frame(p, {{"busy", bs}, {"p", dvec(p, pst.data(), static_cast<int>(s))}})},
                       {"wl", dvec(p, wl.data(), n)},
                       {"disp", dm},
                       {"queue", p(Rf_mkString(qs[w]))},
                       {"pfull", p(Rf_ScalarReal(st[0]))},
                       {"lq", p(Rf_ScalarReal(st[1]))},
                       {"wq", p(Rf_ScalarReal(st[1] / big))},
                       {"rho", p(Rf_ScalarReal(st[2]))},
                       {"travel", trv}});
  Rf_classgets(out, p(Rf_mkString("hypq")));
  return out;
}

SEXP C_clab_r_sqrtlw(SEXP area, SEXP units, SEXP busy, SEXP metric, SEXP speed) {
  clab::who = "sqrtlw";
  prot p;
  static const char *const ms[] = {"euclid", "grid"};
  const double a = pos1(area, "area"), un = pos1(units, "units");
  const double b = REAL(clab::num(busy, "busy", 1, 0.0))[0];
  if (!(b < 1)) fail("`busy` must be below 1: some units must be free.");
  const int mt = csp::pick(metric, "metric", ms, 2);
  double d, c;
  clab_sqrtlw(a, un, b, mt, &d, &c);
  const double tm = Rf_isNull(speed) ? NA_REAL : d / pos1(speed, "speed");
  return nlist(p, {{"dist", p(Rf_ScalarReal(d))},
                   {"const", p(Rf_ScalarReal(c))},
                   {"avail", p(Rf_ScalarReal(un * (1 - b)))},
                   {"time", p(Rf_ScalarReal(tm))}});
}

}  // extern "C"

