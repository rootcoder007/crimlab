// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of crimlab's time, contagion and operations methods:
// the Hawkes process (hwkfit, hwkres, hwksim). Each checks its arguments,
// builds the kernel's arrays, calls the kernels and builds the R result.
// Nothing calls back into the R evaluator.
#include <cmath>
#include <cstring>
#include <string>
#include <vector>

#include "crimlab.h"
#include "kutil.h"
#include "rchk.h"
#include "spin.h"

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
