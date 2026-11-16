// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry point of the Indigenous health, published-table, census and
// Calls to Action functions and of polex() (igh.h): the R arguments made
// values of igh::V, the function run, its result made an R object again
// (attributes and class included). Every C++ object is gone before an
// error or a warning is raised. Nothing calls back into the R evaluator.
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <exception>
#include <string>
#include <vector>

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#include "crimlab.h"
#include "igh.h"
#include "rchk.h"

namespace {

using igh::V;

char msg[2048];

std::vector<std::string> strs(SEXP s) {
  std::vector<std::string> o;
  if (Rf_isNull(s) || TYPEOF(s) != STRSXP) return o;
  for (R_xlen_t j = 0; j < XLENGTH(s); ++j) {
    SEXP e = STRING_ELT(s, j);
    o.push_back(e == NA_STRING ? "NA" : Rf_translateCharUTF8(e));
  }
  return o;
}

V fromr(SEXP x) {
  V v;
  switch (TYPEOF(x)) {
    case NILSXP:
      return v;
    case CLOSXP:
    case BUILTINSXP:
    case SPECIALSXP:
      v.fn = true;
      return v;
    case ENVSXP:
      v.env = true;
      return v;
    case REALSXP:
      v.t = V::DBL;
      v.d.assign(REAL(x), REAL(x) + XLENGTH(x));
      break;
    case INTSXP:
      v.t = Rf_isFactor(x) ? V::FAC : V::INT;
      v.i.assign(INTEGER(x), INTEGER(x) + XLENGTH(x));
      if (v.t == V::FAC) v.s = strs(Rf_getAttrib(x, R_LevelsSymbol));
      break;
    case LGLSXP:
      v.t = V::LGL;
      v.i.assign(LOGICAL(x), LOGICAL(x) + XLENGTH(x));
      break;
    case STRSXP:
      v.t = V::STR;
      for (R_xlen_t j = 0; j < XLENGTH(x); ++j) {
        SEXP e = STRING_ELT(x, j);
        v.na.push_back(e == NA_STRING);
        v.s.push_back(e == NA_STRING ? std::string() : std::string(Rf_translateCharUTF8(e)));
      }
      break;
    case VECSXP:
      v.t = V::LIST;
      v.df = Rf_inherits(x, "data.frame");
      for (R_xlen_t j = 0; j < XLENGTH(x); ++j) v.a.push_back(fromr(VECTOR_ELT(x, j)));
      if (v.df) {
        SEXP rn = Rf_getAttrib(x, R_RowNamesSymbol);
        v.nr = static_cast<int>(XLENGTH(x) ? Rf_xlength(VECTOR_ELT(x, 0)) : Rf_xlength(rn));
      }
      break;
    default:
      v.odd = true;
      return v;
  }
  v.nm = strs(Rf_getAttrib(x, R_NamesSymbol));
  v.cls = strs(Rf_getAttrib(x, R_ClassSymbol));
  return v;
}

SEXP mkstrs(const std::vector<std::string> &s, const std::vector<char> *na) {
  SEXP o = PROTECT(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(s.size())));
  for (std::size_t j = 0; j < s.size(); ++j) {
    const bool m = na && j < na->size() && (*na)[j];
    SET_STRING_ELT(o, static_cast<R_xlen_t>(j), m ? NA_STRING : Rf_mkCharCE(s[j].c_str(), CE_UTF8));
  }
  UNPROTECT(1);
  return o;
}

SEXP tor(const V &v) {
  SEXP o;
  switch (v.t) {
    case V::NUL:
      return R_NilValue;
    case V::DBL:
      o = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(v.d.size())));
      if (!v.d.empty()) std::memcpy(REAL(o), v.d.data(), v.d.size() * sizeof(double));
      break;
    case V::INT:
    case V::FAC:
      o = PROTECT(Rf_allocVector(INTSXP, static_cast<R_xlen_t>(v.i.size())));
      for (std::size_t j = 0; j < v.i.size(); ++j) INTEGER(o)[j] = v.i[j];
      if (v.t == V::FAC) Rf_setAttrib(o, R_LevelsSymbol, mkstrs(v.s, nullptr));
      break;
    case V::LGL:
      o = PROTECT(Rf_allocVector(LGLSXP, static_cast<R_xlen_t>(v.i.size())));
      for (std::size_t j = 0; j < v.i.size(); ++j) LOGICAL(o)[j] = v.i[j];
      break;
    case V::STR:
      o = PROTECT(mkstrs(v.s, &v.na));
      break;
    default:
      o = PROTECT(Rf_allocVector(VECSXP, static_cast<R_xlen_t>(v.a.size())));
      for (std::size_t j = 0; j < v.a.size(); ++j) SET_VECTOR_ELT(o, static_cast<R_xlen_t>(j), tor(v.a[j]));
      break;
  }
  if (!v.nm.empty()) Rf_setAttrib(o, R_NamesSymbol, mkstrs(v.nm, nullptr));
  if (v.t != V::LIST && v.nr >= 0 && v.nc >= 0) {
    SEXP d = PROTECT(Rf_allocVector(INTSXP, 2));
    INTEGER(d)[0] = v.nr;
    INTEGER(d)[1] = v.nc;
    Rf_setAttrib(o, R_DimSymbol, d);
    SEXP dn = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(dn, 0, v.rn.empty() ? R_NilValue : mkstrs(v.rn, nullptr));
    SET_VECTOR_ELT(dn, 1, v.cn.empty() ? R_NilValue : mkstrs(v.cn, nullptr));
    Rf_setAttrib(o, R_DimNamesSymbol, dn);
    UNPROTECT(2);
  }
  if (v.t == V::LIST && v.df) {
    const int n = v.a.empty() ? (v.nr > 0 ? v.nr : 0) : static_cast<int>(v.a[0].len());
    // the kept rows' numbers (row names are always the input's positions)
    SEXP rn;
    if (!v.rn.empty()) {
      rn = PROTECT(Rf_allocVector(INTSXP, n));
      for (int j = 0; j < n; ++j) INTEGER(rn)[j] = std::atoi(v.rn[static_cast<std::size_t>(j)].c_str());
    } else {
      rn = PROTECT(Rf_allocVector(INTSXP, 2));
      INTEGER(rn)[0] = NA_INTEGER;
      INTEGER(rn)[1] = -n;
    }
    Rf_setAttrib(o, R_RowNamesSymbol, rn);
    UNPROTECT(1);
  }
  for (const auto &a : v.at) Rf_setAttrib(o, Rf_install(a.first.c_str()), tor(a.second));
  if (!v.cls.empty()) Rf_setAttrib(o, R_ClassSymbol, mkstrs(v.cls, nullptr));
  UNPROTECT(1);
  return o;
}

}  // namespace

extern "C" {

SEXP C_clab_r_igh(SEXP fn, SEXP args) {
  if (TYPEOF(fn) != STRSXP || clab::xlen(fn, "fn") != 1) Rf_error("igh: bad function name");
  const char *f = CHAR(STRING_ELT(fn, 0));
  int k = 0;
  for (; k < igh::nfns; ++k) {
    if (std::strcmp(igh::fns[k].name, f) == 0) break;
  }
  if (k == igh::nfns) Rf_error("igh: no function %s", f);
  clab::who = igh::fns[k].name;
  // a NULL seed draws one from R's generator, as the session's stream would
  int seed = 0;
  bool draw = false;
  SEXP an = Rf_getAttrib(args, R_NamesSymbol);
  for (R_xlen_t j = 0; j < clab::xlen(args, "args"); ++j) {
    if (std::strcmp(CHAR(STRING_ELT(an, j)), "seed") == 0 && Rf_isNull(VECTOR_ELT(args, j))) {
      GetRNGstate();
      seed = static_cast<int>(std::floor(unif_rand() * 2147483647.0));
      PutRNGstate();
      draw = true;
    }
  }
  SEXP res = R_NilValue, wv = R_NilValue;
  int code = CLAB_OK;
  {
    std::vector<std::string> nm;
    std::vector<V> av;
    for (R_xlen_t j = 0; j < clab::xlen(args, "args"); ++j) {
      nm.push_back(CHAR(STRING_ELT(an, j)));
      if (draw && nm.back() == "seed") {
        V s;
        s.t = V::DBL;
        s.d.push_back(seed);
        av.push_back(s);
      } else {
        av.push_back(fromr(VECTOR_ELT(args, j)));
      }
      av.back().given = true;
    }
    V out;
    std::string err;
    std::vector<std::string> warn;
    try {
      code = igh::call(f, nm, av, out, err, warn);
    } catch (const std::exception &e) {
      code = CLAB_EARG;
      err = e.what();
    }
    if (code == CLAB_OK) {
      res = PROTECT(tor(out));
      wv = PROTECT(mkstrs(warn, nullptr));
    } else {
      std::snprintf(msg, sizeof msg, "%s", err.c_str());
    }
  }
  if (code != CLAB_OK) clab::fail("%s", msg);
  for (R_xlen_t j = 0; j < XLENGTH(wv); ++j) {
    Rf_warningcall(R_NilValue, "%s", Rf_translateCharUTF8(STRING_ELT(wv, j)));
  }
  UNPROTECT(2);
  return res;
}

}  // extern "C"
