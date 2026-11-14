// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry point of repeat victimisation (rptvic, src/rpv.cpp): it
// checks the arguments, calls the kernel and builds the R result. Nothing
// calls back into the R evaluator.
#include <cmath>
#include <vector>

#include "bcol.h"
#include "rchk.h"

using clab::fail;

extern "C" {

SEXP C_clab_r_rptvic(SEXP id, SEXP t, SEXP pop, SEXP brk) {
  clab::who = "rptvic";
  int np = 0;
  const R_xlen_t xn = XLENGTH(id);
  if (xn < 1) fail("`id` must hold one incident at least.");
  csp::lab v = bj::labs(id, "id", xn);
  np += 2;
  const int n = bj::len(v.code, "id");
  SEXP tv = R_NilValue;
  if (!Rf_isNull(t)) {
    tv = PROTECT(bj::nums(t, "t", n));
    ++np;
  }
  double pp = NAN;
  if (!Rf_isNull(pop)) {
    pp = bj::inr(pop, "pop", 1.0, R_PosInf);
    if (pp < v.g) fail("`pop` must be at least the number of victims, %d.", v.g);
  }
  SEXP bv = R_NilValue;
  int nb = 0;
  if (!Rf_isNull(brk)) {
    if (Rf_isNull(t)) fail("`brk` needs the times `t`.");
    bv = PROTECT(bj::nums(brk, "brk"));
    ++np;
    nb = bj::len(bv, "brk");
    if (nb < 2) fail("`brk` must hold two breaks at least.");
    for (int j = 1; j < nb; ++j) {
      if (!(REAL(bv)[j] > REAL(bv)[j - 1])) fail("`brk` must increase.");
    }
  }
  std::vector<double> tab(3 * (static_cast<std::size_t>(n) + 1)), st(10), gap(3),
      lag(3 * static_cast<std::size_t>(nb > 1 ? nb - 1 : 0) + 1);
  int K = 0;
  bj::kc(clab_rptvic(INTEGER(v.code), Rf_isNull(tv) ? nullptr : REAL(tv), n, v.g, pp,
                     Rf_isNull(bv) ? nullptr : REAL(bv), nb, tab.data(), &K, st.data(), gap.data(),
                     lag.data()),
         "the kernel refused the input.");
  const bool hp = !std::isnan(pp);
  const int r0 = hp ? 0 : 1, nr = K + 1 - r0;
  const std::size_t R = static_cast<std::size_t>(K) + 1;
  bj::cols f(4, nr);
  int *kk = f.itg("k");
  double *vi = f.dbl("victims"), *in = f.dbl("incidents");
  double *ex = hp ? f.dbl("expected") : nullptr;
  for (int r = 0; r < nr; ++r) {
    const std::size_t q = static_cast<std::size_t>(r + r0);
    kk[r] = r + r0;
    vi[r] = tab[q];
    in[r] = tab[R + q];
    if (ex) ex[r] = tab[2 * R + q];
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  ++np;
  SET_VECTOR_ELT(out, 0, f.done());
  static const char *const sn[] = {"incidents", "victims", "repeats", "rshare", "conc",
                                   "incid", "prev", "rvict", "erep", "ratio"};
  for (double &x : st) x = bj::na(x);
  SET_VECTOR_ELT(out, 1, clab::nvec(10, sn, st.data()));
  int ncol = 1;
  if (!Rf_isNull(tv)) {
    for (double &x : gap) x = bj::na(x);
    static const char *const gn[] = {"n", "mean", "median"};
    SET_VECTOR_ELT(out, 2, clab::nvec(3, gn, gap.data()));
  }
  if (nb > 1) {
    const int m = nb - 1;
    bj::cols l(5, m);
    ++ncol;
    double *lo = l.dbl("lo"), *hi = l.dbl("hi"), *c = l.dbl("n"), *sh = l.dbl("share"), *cu = l.dbl("cum");
    for (int j = 0; j < m; ++j) {
      lo[j] = REAL(bv)[j];
      hi[j] = REAL(bv)[j + 1];
      c[j] = lag[static_cast<std::size_t>(j)];
      sh[j] = bj::na(lag[static_cast<std::size_t>(m + j)]);
      cu[j] = bj::na(lag[static_cast<std::size_t>(2 * m + j)]);
    }
    SET_VECTOR_ELT(out, 3, l.done());
  }
  static const char *const on[] = {"tab", "stat", "gap", "lag"};
  bj::named(out, on);
  Rf_classgets(out, PROTECT(Rf_mkString("rptvic")));
  UNPROTECT(np + ncol + 1);
  return out;
}

}  // extern "C"
