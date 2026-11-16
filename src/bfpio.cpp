// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of what the footprint reads from the machine
// (src/fpio.cpp). The R wrappers (R/fp.R) check the arguments; these call
// the readers and build the results. Nothing calls back into the R
// evaluator.
#include <cmath>
#include <cstddef>
#include <cstring>
#include <vector>

#include "bcol.h"
#include "rchk.h"

namespace {

// One string argument as UTF-8 text, "" for NULL, an empty vector or NA.
const char *s8(SEXP x) {
  if (!Rf_isString(x) || XLENGTH(x) < 1 || STRING_ELT(x, 0) == NA_STRING) return "";
  return Rf_translateCharUTF8(STRING_ELT(x, 0));
}

// One path argument in the native encoding, as R's file functions take it.
const char *fpath(SEXP x) {
  if (!Rf_isString(x) || XLENGTH(x) < 1 || STRING_ELT(x, 0) == NA_STRING) return "";
  return Rf_translateChar(STRING_ELT(x, 0));
}

void ign(const char *var, SEXP v) {
  Rf_warningcall(R_NilValue,
                 "dtctloc: %s='%s' is not a location co2inttab() knows (code, alpha-3 or name); it is "
                 "ignored.",
                 var, Rf_translateChar(STRING_ELT(v, 0)));
}

}  // namespace

extern "C" {

SEXP C_clab_r_fploc(SEXP env, SEXP cod, SEXP tz, SEXP lc) {
  clab::who = "dtctloc";
  const char *loc = nullptr;
  std::size_t len = 0;
  int ig = 0;
  std::vector<char> how(256);
  int rc = clab_fploc(s8(env), s8(cod), s8(tz), s8(lc), &loc, how.data(), how.size(), &len, &ig);
  if (rc == CLAB_OK && len >= how.size()) {
    how.resize(len + 1);
    rc = clab_fploc(s8(env), s8(cod), s8(tz), s8(lc), &loc, how.data(), how.size(), &len, &ig);
  }
  if (rc != CLAB_OK) clab::fail("%s", clab_strerror(rc));
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(r, 0, Rf_ScalarString(Rf_mkCharCE(loc, CE_UTF8)));
  SET_VECTOR_ELT(r, 1, Rf_ScalarString(Rf_mkCharCE(how.data(), CE_UTF8)));
  static const char *const nm[] = {"location", "method"};
  clab::nlist(r, nm);
  if ((ig & 1) != 0) ign("CRIMLAB_LOCATION", env);
  if ((ig & 2) != 0) ign("CODECARBON_COUNTRY_ISO_CODE", cod);
  UNPROTECT(1);
  return r;
}

SEXP C_clab_r_cizone(SEXP key) {
  int zone = 0;
  const int i = clab_cizone(Rf_translateCharUTF8(STRING_ELT(key, 0)), &zone);
  SEXP r = PROTECT(Rf_allocVector(INTSXP, 2));
  INTEGER(r)[0] = i < 0 ? NA_INTEGER : i + 1;
  INTEGER(r)[1] = zone;
  UNPROTECT(1);
  return r;
}

SEXP C_clab_r_rapldom(SEXP dir) {
  clab::who = "co2fp";
  const char *d = fpath(dir);
  std::size_t len = 0;
  int n = 0;
  std::vector<char> b(4096);
  int rc = clab_rapldom(d, b.data(), b.size(), &len, &n);
  while (rc == CLAB_OK && len > b.size()) {
    b.resize(len);
    rc = clab_rapldom(d, b.data(), b.size(), &len, &n);
  }
  if (rc != CLAB_OK) clab::fail("%s", clab_strerror(rc));
  SEXP r = PROTECT(Rf_allocVector(STRSXP, n));
  std::size_t at = 0;
  for (int i = 0; i < n; ++i) {
    SET_STRING_ELT(r, i, Rf_mkCharCE(b.data() + at, CE_NATIVE));
    at += std::strlen(b.data() + at) + 1;
  }
  UNPROTECT(1);
  return r;
}

SEXP C_clab_r_raplrd(SEXP dom) {
  const R_xlen_t n = Rf_isString(dom) ? XLENGTH(dom) : 0;
  if (n == 0) return R_NilValue;
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP uj = Rf_allocVector(REALSXP, n);
  SET_VECTOR_ELT(r, 0, uj);
  SEXP mx = Rf_allocVector(REALSXP, n);
  SET_VECTOR_ELT(r, 1, mx);
  for (R_xlen_t i = 0; i < n; ++i) {
    clab_raplrd(Rf_translateChar(STRING_ELT(dom, i)), REAL(uj) + i, REAL(mx) + i);
    if (std::isnan(REAL(uj)[i])) {
      UNPROTECT(1);
      return R_NilValue;
    }
    REAL(mx)[i] = bj::na(REAL(mx)[i]);
  }
  static const char *const nm[] = {"uj", "max"};
  clab::nlist(r, nm);
  UNPROTECT(1);
  return r;
}

SEXP C_clab_r_cpuacc(SEXP path) {
  double v[2] = {0.0, 0.0};
  if (clab_cpuacc(fpath(path), v) != CLAB_OK) return R_NilValue;
  static const char *const nm[] = {"idle", "total"};
  return clab::nvec(2, nm, v);
}

SEXP C_clab_r_memgb(void) {
  const double g = clab_memgb();
  return Rf_ScalarReal(std::isnan(g) ? NA_REAL : g);
}

SEXP C_clab_r_fpnote(SEXP flags, SEXP machine) {
  const int fl = Rf_asInteger(flags), mc = Rf_asLogical(machine) == TRUE ? 1 : 0;
  int n = 0;
  for (int i = 0; i < 5; ++i) n += clab_fpnote(i, fl, mc) != nullptr ? 1 : 0;
  SEXP r = PROTECT(Rf_allocVector(STRSXP, n));
  int k = 0;
  for (int i = 0; i < 5; ++i) {
    const char *t = clab_fpnote(i, fl, mc);
    if (t != nullptr) SET_STRING_ELT(r, k++, Rf_mkChar(t));
  }
  UNPROTECT(1);
  return r;
}

}  // extern "C"
