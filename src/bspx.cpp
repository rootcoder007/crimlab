// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the spatial crime methods: near repeats (src/knx.cpp)
// and the rest that follow. Each reads its points through csp::xy()
// (spin.h), checks the other arguments, calls the kernels and builds the R
// result. Nothing here calls back into the R evaluator.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <vector>

#include "crimlab.h"
#include "rchk.h"
#include "spin.h"

// R's distribution functions under their library names (Rmath.h declares
// them so only in newer versions)
extern "C" double Rf_ppois(double, double, int, int);
extern "C" double Rf_pnorm5(double, double, double, int, int);
extern "C" double Rf_pchisq(double, double, int, int);
extern "C" double Rf_qchisq(double, double, int, int);

using clab::fail;
using clab::num;

namespace {

SEXP real(R_xlen_t n) { return Rf_allocVector(REALSXP, n); }

void kc(int rc, const char *msg) {
  if (rc != CLAB_OK) fail("%s", msg);
}

// a data frame of k double columns of length n, named nm; col[j] gives
// each column's storage
SEXP dframe(int k, int n, const char *const *nm, double **col) {
  SEXP v = PROTECT(Rf_allocVector(VECSXP, k));
  for (int j = 0; j < k; ++j) {
    SET_VECTOR_ELT(v, j, real(n));
    col[j] = REAL(VECTOR_ELT(v, j));
  }
  SEXP f = clab::frame(v, nm);
  UNPROTECT(1);
  return f;
}

// band edges: positive, finite and increasing (returned protected)
SEXP edges(SEXP e, const char *nm) {
  SEXP v = PROTECT(num(e, nm, -1, 0.0, true));
  const R_xlen_t k = XLENGTH(v);
  for (R_xlen_t q = 1; q < k; ++q) {
    if (!(REAL(v)[q] > REAL(v)[q - 1])) fail("`%s` must increase.", nm);
  }
  return v;
}

// a classed list
SEXP klass(SEXP x, const char *cl) {
  Rf_classgets(x, PROTECT(Rf_mkString(cl)));
  UNPROTECT(1);
  return x;
}

}  // namespace

extern "C" {

// ---------------------------------------------------------- near repeats

SEXP C_clab_r_knxnr(SEXP x, SEXP y, SEXP t, SEXP ds, SEXP dt, SEXP nsim, SEXP seed, SEXP lonlat) {
  clab::who = "knxnr";
  csp::pts p = csp::xy(x, y, lonlat);
  const int n = p.n;
  if (n < 3) fail("at least 3 incidents are needed.");
  SEXP tv = PROTECT(num(t, "t", n));
  SEXP sv = edges(ds, "ds");
  SEXP bv = edges(dt, "dt");
  const int ks = static_cast<int>(XLENGTH(sv)), kt = static_cast<int>(XLENGTH(bv));
  const int ns = clab::whole(nsim, "nsim", 0);
  const unsigned int sd = ns > 0 ? csp::seed(seed) : 0u;
  const int m = ks * kt;
  SEXP o = PROTECT(real(3 * static_cast<R_xlen_t>(m) + 1));
  kc(clab_knxnr(p.x, p.y, REAL(tv), n, REAL(sv), ks, REAL(bv), kt, ns, sd, REAL(o)),
     "the kernel refused the input.");
  static const char *const tn[] = {"dlo", "dhi", "tlo", "thi", "obs", "exp", "ratio", "p", "ppois"};
  double *c[9];
  SEXP tab = PROTECT(dframe(9, m, tn, c));
  const double *ob = REAL(o), *ex = REAL(o) + m, *pp = REAL(o) + 2 * m;
  for (int b = 0; b < kt; ++b) {
    for (int s = 0; s < ks; ++s) {
      const int q = s + b * ks;
      c[0][q] = s ? REAL(sv)[s - 1] : 0.0;
      c[1][q] = REAL(sv)[s];
      c[2][q] = b ? REAL(bv)[b - 1] : 0.0;
      c[3][q] = REAL(bv)[b];
      c[4][q] = ob[q];
      c[5][q] = ex[q];
      c[6][q] = ex[q] > 0 ? ob[q] / ex[q] : std::numeric_limits<double>::quiet_NaN();
      c[7][q] = pp[q];
      // P(X >= obs) for X Poisson with the expected count as its mean
      c[8][q] = ob[q] > 0 ? Rf_ppois(ob[q] - 1, ex[q], 0, 0) : 1.0;
    }
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out, 0, tab);
  SET_VECTOR_ELT(out, 1, Rf_ScalarInteger(n));
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(REAL(o)[3 * m]));
  SET_VECTOR_ELT(out, 3, Rf_ScalarInteger(ns));
  static const char *const on[] = {"tab", "n", "pairs", "nsim"};
  clab::nlist(out, on);
  UNPROTECT(7);
  return klass(out, "knxnr");
}

}  // extern "C"
