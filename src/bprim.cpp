// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry point of the summary statistics (src/prim.cpp): R/prim.R
// has checked every argument and passes doubles; op names the kernel.
// Nothing calls back into the R evaluator.
#include <vector>

#include "bcol.h"
#include "rchk.h"

extern "C" {

SEXP C_clab_r_prim(SEXP op, SEXP x, SEXP y, SEXP a) {
  const int k = Rf_asInteger(op);
  const int n = static_cast<int>(clab::xlen(x, "x"));
  const double *px = REAL(x);
  const double *pa = Rf_isNull(a) ? nullptr : REAL(a);
  SEXP out = R_NilValue;
  auto vec = [&out](R_xlen_t m, const char *nm = "x") {
    out = PROTECT(clab::xalloc(REALSXP, m, nm));
    return REAL(out);
  };
  switch (k) {
    case 0: *vec(1) = clab_mean(px, n); break;
    case 1: *vec(1) = clab_var(px, n, static_cast<int>(pa[0])); break;
    case 2: *vec(1) = clab_cor(px, REAL(y), n); break;
    case 3: *vec(1) = clab_edist(px, REAL(y), n); break;
    case 4: clab_dnorm(px, n, pa[0], pa[1], static_cast<int>(pa[2]), vec(n)); break;
    case 5: clab_moms(px, n, vec(4)); break;
    case 6: clab_quant(px, n, REAL(y), static_cast<int>(clab::xlen(y, "y")), vec(clab::xlen(y, "y"))); break;
    case 7: *vec(1) = clab_mad(px, n, pa[0]); break;
    case 8: *vec(1) = clab_trim(px, n, pa[0]); break;
    case 9: *vec(1) = clab_wins(px, n, pa[0]); break;
    case 10: clab_wtd(px, REAL(y), n, vec(2)); break;
    case 11: *vec(1) = clab_spear(px, REAL(y), n); break;
    case 12: clab_midrk(px, n, vec(n)); break;
    case 13: {
      const int p = static_cast<int>(pa[0]), r = p ? n / p : 0;
      clab_covm(px, r, p, vec(static_cast<R_xlen_t>(p) * p));
      break;
    }
    case 14: {
      const int b = static_cast<int>(pa[0]);
      clab::who = "cbootmean";
      double *o = vec(b, "B");
      clab::kchk(clab_bmean(px, n, b, pa[1], o), "bootstrap");
      break;
    }
    case 15: clab_ipw(px, REAL(y), n, pa[0], pa[1], vec(n)); break;
    case 16:
      clab::who = "cunif";
      clab_unif(pa[1], static_cast<int>(pa[0]), vec(static_cast<R_xlen_t>(pa[0]), "n"));
      break;
    case 18: {
      double *o = vec(3 * static_cast<R_xlen_t>(n));
      clab::kchk(clab_shrs(px, n, pa[0], pa[1], o, o + n, o + 2 * n), "shrs");
      break;
    }
    case 19: clab_ks(px, n, REAL(y), static_cast<int>(clab::xlen(y, "y")), vec(3)); break;
    case 20: {
      const int b = static_cast<int>(pa[0]), ny = static_cast<int>(clab::xlen(y, "y"));
      std::vector<double> p(static_cast<std::size_t>(b) + 2), q(p.size());
      const int nbin = clab_psibin(px, n, REAL(y), ny, b, p.data(), q.data());
      double r[2];
      clab_psi(p.data(), q.data(), nbin, pa[1], r);
      double *o = vec(2);
      o[0] = r[0];
      o[1] = r[1];
      break;
    }
    case 21: clab_chig(px, REAL(y), n, vec(3)); break;
    case 22: clab_chih(px, REAL(y), n, vec(4)); break;
    case 23: clab_ysum(px, REAL(y), n, pa[0], vec(8)); break;
    default: {
      double *o = vec(n);
      for (int i = 0; i < n; ++i) o[i] = clab_pgamma(px[i], pa[0], 1);
    }
  }
  UNPROTECT(1);
  return out;
}

}  // extern "C"
