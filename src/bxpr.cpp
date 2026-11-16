// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the experiments, the trials and the survey tests
// (src/rnd.cpp and the kernels after it). The R wrappers check the
// arguments and pass doubles, integers, codes from 0 and a seed that is a
// whole number; these call the kernels and return their arrays. Nothing
// calls back into the R evaluator.
#include <cmath>
#include <vector>

#include "crimlab.h"
#include "pev.h"
#include "pevr.h"
#include "rchk.h"
#include "spin.h"

namespace {

void kc(int rc, const char *dom = nullptr) {
  clab::ichk(rc);
  if (rc == CLAB_EDOM && dom) clab::fail("%s", dom);
  if (rc != CLAB_OK) clab::fail("%s", clab_strerror(rc));
}

SEXP lst(int n, const char *const *nm, SEXP *v) {
  SEXP out = PROTECT(Rf_allocVector(VECSXP, n));
  SEXP s = PROTECT(Rf_allocVector(STRSXP, n));
  for (int i = 0; i < n; ++i) {
    SET_VECTOR_ELT(out, i, v[i]);
    SET_STRING_ELT(s, i, Rf_mkChar(nm[i]));
  }
  Rf_setAttrib(out, R_NamesSymbol, s);
  UNPROTECT(2);
  return out;
}

unsigned int sd(SEXP seed) { return static_cast<unsigned int>(Rf_asReal(seed)); }

}  // namespace

extern "C" {

SEXP C_clab_r_rndalc(SEXP cnt, SEXP ratio, SEXP block, SEXP seed) {
  clab::who = "rndalc";
  const int ns = Rf_length(cnt), k = Rf_length(ratio), nb = Rf_length(block);
  R_xlen_t tot = 0;
  for (int s = 0; s < ns; ++s) tot += INTEGER(cnt)[s];
  SEXP v[2];
  v[0] = PROTECT(Rf_allocVector(INTSXP, tot));
  v[1] = PROTECT(Rf_allocVector(INTSXP, tot));
  kc(clab_rndalc(ns, INTEGER(cnt), k, INTEGER(ratio), INTEGER(block), nb, sd(seed), INTEGER(v[0]),
                 INTEGER(v[1])));
  static const char *const nms[] = {"arm", "block"};
  SEXP out = lst(2, nms, v);
  UNPROTECT(2);
  return out;
}

SEXP C_clab_r_blkprm(SEXP y, SEXP d, SEXP blk, SEXP nb, SEXP stat, SEXP alt, SEXP nperm, SEXP exact,
                     SEXP seed) {
  clab::who = "blkprm";
  const int n = Rf_length(y), ex = Rf_asLogical(exact);
  // the arrangements: nperm draws, or the exact count the wrapper checked
  const int nd = Rf_asInteger(nperm);
  SEXP v[2];
  v[0] = PROTECT(Rf_allocVector(REALSXP, 4));
  v[1] = PROTECT(Rf_allocVector(REALSXP, nd));
  kc(clab_blkprm(REAL(y), INTEGER(d), INTEGER(blk), n, Rf_asInteger(nb), Rf_asInteger(stat),
                 Rf_asInteger(alt), nd, ex, ex ? 0u : sd(seed), nd, REAL(v[0]), REAL(v[1])));
  static const char *const nms[] = {"st", "dist"};
  SEXP out = lst(2, nms, v);
  UNPROTECT(2);
  return out;
}

SEXP C_clab_r_pwranv(SEXP k, SEXP n, SEXP f, SEXP alpha, SEXP power) {
  clab::who = "pwranv";
  const int m = Rf_length(f);
  SEXP o = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(m) * 6));
  kc(clab_pwranv(Rf_asInteger(k), m, Rf_isNull(n) ? nullptr : REAL(n), REAL(f), Rf_asReal(alpha),
                 Rf_asReal(power), REAL(o)));
  UNPROTECT(1);
  return o;
}

SEXP C_clab_r_pwrsim(SEXP fam, SEXP par, SEXP n0, SEXP n1, SEXP eff, SEXP test, SEXP alt, SEXP alpha,
                     SEXP nsim, SEXP nperm, SEXP level, SEXP seed) {
  clab::who = "pwrsim";
  const int m = Rf_length(eff);
  SEXP o = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(m) * 5));
  kc(clab_pwrsim(Rf_asInteger(fam), Rf_asReal(par), m, INTEGER(n0), INTEGER(n1), REAL(eff),
                 Rf_asInteger(test), Rf_asInteger(alt), Rf_asReal(alpha), Rf_asInteger(nsim),
                 Rf_asInteger(nperm), Rf_asReal(level), sd(seed), REAL(o)));
  UNPROTECT(1);
  return o;
}

SEXP C_clab_r_falsctl(SEXP x, SEXP d, SEXP y, SEXP ctl, SEXP hc, SEXP level, SEXP adj, SEXP nperm,
                      SEXP seed) {
  clab::who = "falsctl";
  const int n = Rf_length(d), p = Rf_ncols(x), m = Rf_ncols(y), np = Rf_asInteger(nperm);
  SEXP o = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(m) * 8));
  kc(clab_falsctl(REAL(x), n, p, REAL(d), REAL(y), m, INTEGER(ctl), Rf_asLogical(hc), Rf_asReal(level),
                  Rf_asInteger(adj), np, np > 0 ? sd(seed) : 0u, REAL(o)),
     "the covariates are collinear, or `d` does not vary once they are taken out.");
  UNPROTECT(1);
  return o;
}

SEXP C_clab_r_plntfx(SEXP x, SEXP y, SEXP d, SEXP eff, SEXP nrep, SEXP hc, SEXP level, SEXP seed) {
  clab::who = "plntfx";
  const int n = Rf_length(y), p = Rf_ncols(x), ne = Rf_length(eff);
  SEXP o = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(ne) * 5));
  kc(clab_plntfx(REAL(x), n, p, REAL(y), REAL(d), REAL(eff), ne, Rf_asInteger(nrep), Rf_asLogical(hc),
                 Rf_asReal(level), sd(seed), REAL(o)),
     "the covariates are collinear, or `d` does not vary once they are taken out.");
  UNPROTECT(1);
  return o;
}

SEXP C_clab_r_specrv(SEXP x, SEXP nx, SEXP m, SEXP y, SEXP d, SEXP hc, SEXP level, SEXP boot, SEXP nsim,
                     SEXP seed) {
  clab::who = "specrv";
  const int n = Rf_length(d), ny = Rf_ncols(y), mm = Rf_asInteger(m), ns = Rf_asInteger(nsim);
  const R_xlen_t K = static_cast<R_xlen_t>(ny) << mm;
  SEXP v[3];
  v[0] = PROTECT(Rf_allocVector(REALSXP, K * 6));
  v[1] = PROTECT(Rf_allocVector(REALSXP, 8));
  v[2] = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(ns) * 3));
  kc(clab_specrv(REAL(x), n, Rf_asInteger(nx), mm, REAL(y), ny, REAL(d), Rf_asLogical(hc), Rf_asReal(level),
                 Rf_asLogical(boot), ns, sd(seed), REAL(v[0]), REAL(v[1]), REAL(v[2])),
     "the covariates are collinear, or `d` does not vary once they are taken out.");
  static const char *const nms[] = {"spec", "st", "dist"};
  SEXP out = lst(3, nms, v);
  UNPROTECT(3);
  return out;
}

SEXP C_clab_r_gensw(SEXP x, SEXP s, SEXP target, SEXP trim, SEXP bins) {
  clab::who = "gensw";
  const int n = Rf_length(s), p = Rf_ncols(x);
  SEXP v[5];
  v[0] = PROTECT(Rf_allocVector(REALSXP, n));
  v[1] = PROTECT(Rf_allocVector(REALSXP, n));
  v[2] = PROTECT(Rf_allocVector(REALSXP, p + 1));
  v[3] = PROTECT(Rf_allocVector(REALSXP, 6));
  v[4] = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(p) * 5));
  kc(clab_gensw(REAL(x), n, p, INTEGER(s), Rf_asLogical(target), Rf_asReal(trim), Rf_asInteger(bins),
                REAL(v[0]), REAL(v[1]), REAL(v[2]), REAL(v[3]), REAL(v[4])),
     "the sampling model cannot be fitted: the covariates separate the trial from the rest, or are "
     "collinear.");
  static const char *const nms[] = {"ps", "w", "coef", "st", "bal"};
  SEXP out = lst(5, nms, v);
  UNPROTECT(5);
  return out;
}

SEXP C_clab_r_trnsp(SEXP x, SEXP y, SEXP d, SEXP s, SEXP gen, SEXP trim, SEXP nboot, SEXP seed) {
  clab::who = "trnsp";
  const int n = Rf_length(y), p = Rf_ncols(x), nb = Rf_asInteger(nboot);
  SEXP v[6];
  v[0] = PROTECT(Rf_allocVector(REALSXP, 2));
  v[1] = PROTECT(Rf_allocVector(REALSXP, 2));
  v[2] = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(nb) * 2));
  v[3] = PROTECT(Rf_allocVector(REALSXP, n));
  v[4] = PROTECT(Rf_allocVector(REALSXP, n));
  v[5] = PROTECT(Rf_allocVector(REALSXP, 1));
  // gen: TRUE to generalise, which clab_gensw() numbers 0
  kc(clab_trnsp(REAL(x), n, p, REAL(y), INTEGER(d), INTEGER(s), Rf_asLogical(gen) ? 0 : 1, Rf_asReal(trim), nb,
                nb > 0 ? sd(seed) : 0u, REAL(v[0]), REAL(v[1]), REAL(v[2]), REAL(v[3]), REAL(v[4]), REAL(v[5])),
     "a model cannot be fitted: an arm of the trial has too few rows for the covariates, or the covariates "
     "separate the trial from the rest.");
  static const char *const nms[] = {"est", "se", "draws", "ps", "w", "nfail"};
  SEXP out = lst(6, nms, v);
  UNPROTECT(6);
  return out;
}

SEXP C_clab_r_rcttml(SEXP x, SEXP y, SEXP d, SEXP bin, SEXP inter, SEXP gx, SEXP g, SEXP gbound, SEXP bounds,
                     SEXP qbound) {
  clab::who = "rcttml";
  const int n = Rf_length(y), p = Rf_ncols(x);
  SEXP v[2];
  v[0] = PROTECT(Rf_allocVector(REALSXP, 6));
  v[1] = PROTECT(Rf_allocVector(REALSXP, 3));
  kc(clab_rcttml(REAL(x), n, p, REAL(y), INTEGER(d), Rf_asLogical(bin), Rf_asLogical(inter), Rf_asLogical(gx),
                 Rf_asReal(g), Rf_asReal(gbound), REAL(bounds)[0], REAL(bounds)[1], Rf_asReal(qbound), REAL(v[0]),
                 REAL(v[1])),
     "the outcome or treatment model cannot be fitted: the covariates are collinear, or separate the outcome or "
     "the arms.");
  static const char *const nms[] = {"out", "aux"};
  SEXP out = lst(2, nms, v);
  UNPROTECT(2);
  return out;
}

SEXP C_clab_r_rscchi(SEXP x, SEXP y, SEXP w, SEXP strata, SEXP psu, SEXP fpc) {
  clab::who = "rscchi";
  using pevr::z;
  const int n = clab::xlen(x, "x");
  if (n < 2) clab::fail("`x` must have two values at least.");
  int np = 0;
  csp::lab lx = csp::codes(x, "x", n);
  csp::lab ly = csp::codes(y, "y", n);
  np += 4;
  if (lx.g < 2 || ly.g < 2) clab::fail("`x` and `y` must each hold two labels or more.");
  std::vector<double> wv(z(n), 1.0);
  if (!Rf_isNull(w)) {
    SEXP wx = PROTECT(clab::num(w, "w", n, 0.0, true));
    ++np;
    wv.assign(REAL(wx), REAL(wx) + n);
  }
  // strata and first-stage units, as svymx() reads them
  std::vector<int> st(z(n), 0), ps(z(n));
  int nst = 1;
  if (!Rf_isNull(strata)) {
    csp::lab l = csp::codes(strata, "strata", n);
    st.assign(INTEGER(l.code), INTEGER(l.code) + n);
    nst = l.g;
    UNPROTECT(2);
  }
  std::vector<int> pu(z(n));
  int npu = n;
  if (Rf_isNull(psu)) {
    for (int i = 0; i < n; ++i) pu[z(i)] = i;
  } else {
    csp::lab l = csp::codes(psu, "psu", n);
    pu.assign(INTEGER(l.code), INTEGER(l.code) + n);
    npu = l.g;
    UNPROTECT(2);
  }
  std::vector<int> iw(z(nst) * z(npu));
  const int nps = pev::pairs(st.data(), pu.data(), n, nst, npu, ps.data(), iw.data());
  std::vector<double> nh(z(nst), 0.0);
  {
    std::vector<int> sp(z(nps), -1);
    for (int i = 0; i < n; ++i) sp[z(ps[z(i)])] = st[z(i)];
    for (int a : sp) nh[z(a)] += 1.0;
  }
  for (double v : nh) {
    if (v < 2) clab::fail("a stratum has one first-stage unit: its variance cannot be estimated.");
  }
  std::vector<double> fv;
  if (!Rf_isNull(fpc)) {
    SEXP fx = PROTECT(clab::num(fpc, "fpc", n, 0.0, true));
    ++np;
    fv.assign(z(nst), -1.0);
    for (int i = 0; i < n; ++i) {
      double &f = fv[z(st[z(i)])];
      if (f >= 0 && f != REAL(fx)[i]) clab::fail("`fpc` must be constant within each stratum.");
      f = REAL(fx)[i];
    }
    for (int h = 0; h < nst; ++h) {
      if (fv[z(h)] < nh[z(h)]) clab::fail("`fpc` must be at least the first-stage units sampled in its stratum.");
    }
  }
  const int K = (lx.g - 1) * (ly.g - 1);
  SEXP v[8];
  v[0] = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(lx.g) * ly.g));
  v[1] = lx.lev;
  v[2] = ly.lev;
  v[3] = PROTECT(Rf_allocVector(REALSXP, 16));
  v[4] = PROTECT(Rf_allocVector(REALSXP, K));
  v[5] = PROTECT(Rf_allocVector(REALSXP, 3));
  v[6] = PROTECT(Rf_ScalarInteger(nst));
  v[7] = PROTECT(Rf_ScalarInteger(nps));
  np += 6;
  kc(clab_rscchi(INTEGER(lx.code), INTEGER(ly.code), lx.g, ly.g, n, wv.data(), st.data(), nst, ps.data(), nps,
                 fv.empty() ? nullptr : fv.data(), REAL(v[0]), REAL(v[3]), REAL(v[4]), REAL(v[5])),
     "the design effects cannot be found: a margin of the table, or its variance, is 0.");
  static const char *const nms[] = {"tab", "xlev", "ylev", "test", "deff", "aux", "nstrata", "npsu"};
  SEXP out = lst(8, nms, v);
  UNPROTECT(np);
  return out;
}

SEXP C_clab_r_mipmm(SEXP x, SEXP tg, SEXP m, SEXP k, SEXP maxit, SEXP ridge, SEXP seed, SEXP fy, SEXP fx,
                    SEXP level) {
  clab::who = "mipmm";
  const int n = Rf_nrows(x), q = Rf_ncols(x), nt = Rf_length(tg), mm = Rf_asInteger(m), f = Rf_asInteger(fy),
            nfx = Rf_length(fx);
  R_xlen_t nmis = 0;
  for (int j = 0; j < nt; ++j) {
    const double *c = REAL(x) + static_cast<R_xlen_t>(INTEGER(tg)[j]) * n;
    for (int i = 0; i < n; ++i) nmis += std::isnan(c[i]) ? 1 : 0;
  }
  SEXP v[2];
  v[0] = PROTECT(Rf_allocVector(REALSXP, nmis * mm));
  v[1] = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(f >= 0 ? nfx + 1 : nt) * 8));
  kc(clab_mipmm(REAL(x), n, q, INTEGER(tg), nt, mm, Rf_asInteger(k), Rf_asInteger(maxit), Rf_asReal(ridge),
                sd(seed), f, INTEGER(fx), nfx, Rf_asReal(level), REAL(v[0]), REAL(v[1])),
     "a column of `y` has too few observed values for its predictors, or they are collinear.");
  static const char *const nms[] = {"imp", "pool"};
  SEXP out = lst(2, nms, v);
  UNPROTECT(2);
  return out;
}

}  // extern "C"
