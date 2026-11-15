// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of search theory for search and rescue (src/sar.cpp).
// The R wrappers (R/sar.R) check the arguments; these call the kernels and
// build the results. Nothing calls back into the R evaluator.
#include <vector>

#include "bcol.h"
#include "rchk.h"

namespace {

int n32(SEXP x) { return clab::ilen(static_cast<double>(XLENGTH(x)), nullptr); }

void kc(int rc) {
  clab::ichk(rc);
  if (rc != CLAB_OK) clab::fail("%s", clab_strerror(rc));
}

SEXP dvec(const std::vector<double> &v) {
  SEXP r = Rf_allocVector(REALSXP, static_cast<R_xlen_t>(v.size()));
  for (std::size_t i = 0; i < v.size(); ++i) REAL(r)[i] = v[i];
  return r;
}

}  // namespace

extern "C" {

SEXP C_clab_r_sarlr(SEXP x, SEXP W, SEXP model, SEXP m, SEXP navsd) {
  clab::who = "ltrlrng";
  const int n = n32(x);
  SEXP p = PROTECT(Rf_allocVector(REALSXP, n));
  kc(clab_sarlr(REAL(x), n, Rf_asReal(W), Rf_asInteger(model), Rf_asReal(m), Rf_asReal(navsd), REAL(p)));
  UNPROTECT(1);
  return p;
}

SEXP C_clab_r_sarsw(SEXP x, SEXP p) {
  clab::who = "swpwdth";
  double w = 0;
  kc(clab_sarsw(REAL(x), REAL(p), n32(x), &w));
  return Rf_ScalarReal(w);
}

SEXP C_clab_r_sarpod(SEXP cv, SEXP model, SEXP m) {
  clab::who = "srchpod";
  const int n = n32(cv);
  SEXP p = PROTECT(Rf_allocVector(REALSXP, n));
  kc(clab_sarpod(REAL(cv), n, Rf_asInteger(model), Rf_asReal(m), REAL(p)));
  UNPROTECT(1);
  return p;
}

SEXP C_clab_r_sarpsw(SEXP W, SEXP S, SEXP model, SEXP m, SEXP navsd, SEXP npts, SEXP ntrk) {
  clab::who = "pswppod";
  double p = 0;
  kc(clab_sarpsw(Rf_asReal(W), Rf_asReal(S), Rf_asInteger(model), Rf_asReal(m), Rf_asReal(navsd),
                 Rf_asInteger(npts), Rf_asInteger(ntrk), &p));
  return Rf_ScalarReal(p);
}

SEXP C_clab_r_sarscs(SEXP poc, SEXP D) {
  clab::who = "srchscs";
  const int na = n32(poc), ns = Rf_nrows(D);
  std::vector<double> pa(static_cast<std::size_t>(na)), up(pa), cp(pa);
  double pos = 0;
  int cert = 0;
  kc(clab_sarscs(REAL(poc), na, REAL(D), ns, &pos, pa.data(), up.data(), cp.data(), &cert));
  for (double &v : up) v = bj::na(v);
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 5));
  SET_VECTOR_ELT(r, 0, Rf_ScalarReal(pos));
  SET_VECTOR_ELT(r, 1, dvec(pa));
  SET_VECTOR_ELT(r, 2, dvec(up));
  SET_VECTOR_ELT(r, 3, dvec(cp));
  SET_VECTOR_ELT(r, 4, Rf_ScalarLogical(cert));
  UNPROTECT(1);
  return r;
}

SEXP C_clab_r_sarmap(SEXP xe, SEXP ye, SEXP datum, SEXP cen, SEXP sigma, SEXP end, SEXP nodes, SEXP wts) {
  clab::who = "datumpmap";
  const int nxe = n32(xe), nye = n32(ye);
  SEXP g = PROTECT(Rf_allocMatrix(REALSXP, nye - 1, nxe - 1));
  kc(clab_sarmap(REAL(xe), nxe, REAL(ye), nye, Rf_asInteger(datum), REAL(cen), Rf_asReal(sigma), REAL(end),
                 Rf_asInteger(nodes), Rf_isNull(wts) ? nullptr : REAL(wts), REAL(g)));
  UNPROTECT(1);
  return g;
}

SEXP C_clab_r_sarcrc(SEXP sigma, SEXP W, SEXP effort) {
  clab::who = "optcrcsrch";
  SEXP r = PROTECT(Rf_allocVector(REALSXP, 3));
  kc(clab_sarcrc(Rf_asReal(sigma), Rf_asReal(W), Rf_asReal(effort), REAL(r)));
  UNPROTECT(1);
  return r;
}

SEXP C_clab_r_saropt(SEXP poc, SEXP area, SEXP W, SEXP effort, SEXP tol) {
  clab::who = "opteffort";
  const int n = n32(poc);
  std::vector<double> z(static_cast<std::size_t>(n)), cv(z), pd(z);
  double pos = 0, lam = 0;
  kc(clab_saropt(REAL(poc), REAL(area), n, Rf_asReal(W), Rf_asReal(effort), Rf_asReal(tol), z.data(), cv.data(),
                 pd.data(), &pos, &lam));
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 5));
  SET_VECTOR_ELT(r, 0, dvec(z));
  SET_VECTOR_ELT(r, 1, dvec(cv));
  SET_VECTOR_ELT(r, 2, dvec(pd));
  SET_VECTOR_ELT(r, 3, Rf_ScalarReal(pos));
  SET_VECTOR_ELT(r, 4, Rf_ScalarReal(lam));
  UNPROTECT(1);
  return r;
}

SEXP C_clab_r_sarsq(SEXP z, SEXP model, SEXP m) {
  clab::who = "optsqr";
  const int k = n32(z);
  std::vector<double> o(static_cast<std::size_t>(k) * 8);
  kc(clab_sarsq(REAL(z), k, Rf_asInteger(model), Rf_asReal(m), o.data()));
  SEXP r = PROTECT(Rf_allocMatrix(REALSXP, k, 7));
  for (int i = 0; i < k; ++i) {
    for (int j = 0; j < 7; ++j) REAL(r)[static_cast<std::size_t>(j) * k + i] = o[static_cast<std::size_t>(i) * 8 + j];
  }
  UNPROTECT(1);
  return r;
}

SEXP C_clab_r_sarfig(void) {
  clab::who = "srchfig";
  const int n = clab_sarfig(-1, nullptr);
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 4));
  for (int k = 0; k < 4; ++k) SET_VECTOR_ELT(r, k, Rf_allocVector(STRSXP, n));
  const char *row[4];
  for (int i = 0; i < n; ++i) {
    clab_sarfig(i, row);
    for (int k = 0; k < 4; ++k) SET_STRING_ELT(VECTOR_ELT(r, k), i, Rf_mkCharCE(row[k], CE_UTF8));
  }
  UNPROTECT(1);
  return r;
}

}  // extern "C"
