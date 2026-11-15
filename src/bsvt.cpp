// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry point of the scaling, ideal-point and count-dynamics kernels
// (svk.h): R's arguments into the node tree, the kernel's result back into
// R lists, vectors, matrices and arrays (NaN as NA). Nothing calls back
// into the R evaluator.
#include <cmath>
#include <cstdio>
#include <string>

#include "crimlab.h"
#include "rchk.h"
#include "svk.h"

namespace {

using svk::Node;

bool fromr(SEXP x, Node &n, std::string &err, const std::string &nm) {
  switch (TYPEOF(x)) {
    case NILSXP:
      n.k = Node::NUL;
      return true;
    case REALSXP:
    case INTSXP:
    case LGLSXP: {
      const R_xlen_t len = XLENGTH(x);
      n.k = TYPEOF(x) == REALSXP ? Node::NUM : TYPEOF(x) == INTSXP ? Node::INT : Node::LGL;
      if (Rf_inherits(x, "factor")) {
        err = "`" + nm + "` must not be a factor.";
        return false;
      }
      n.d.resize(static_cast<std::size_t>(len));
      for (R_xlen_t i = 0; i < len; ++i) {
        double v;
        if (TYPEOF(x) == REALSXP) {
          v = REAL(x)[i];
        } else {
          const int e = TYPEOF(x) == INTSXP ? INTEGER(x)[i] : LOGICAL(x)[i];
          v = e == NA_INTEGER ? NAN : static_cast<double>(e);
        }
        n.d[static_cast<std::size_t>(i)] = v;
      }
      SEXP dm = Rf_getAttrib(x, R_DimSymbol);
      if (!Rf_isNull(dm)) {
        for (R_xlen_t i = 0; i < XLENGTH(dm); ++i) n.dim.push_back(INTEGER(dm)[i]);
      }
      return true;
    }
    case STRSXP:
      n.k = Node::STR;
      for (R_xlen_t i = 0; i < XLENGTH(x); ++i) {
        if (STRING_ELT(x, i) == NA_STRING) {
          err = "`" + nm + "` must not hold NA.";
          return false;
        }
        n.s.emplace_back(Rf_translateCharUTF8(STRING_ELT(x, i)));
      }
      return true;
    case VECSXP: {
      n.k = Node::LIST;
      SEXP nms = Rf_getAttrib(x, R_NamesSymbol);
      for (R_xlen_t i = 0; i < XLENGTH(x); ++i) {
        Node e;
        if (!fromr(VECTOR_ELT(x, i), e, err, nm)) return false;
        if (!Rf_isNull(nms)) n.nm.emplace_back(Rf_translateCharUTF8(STRING_ELT(nms, i)));
        n.l.push_back(std::move(e));
      }
      return true;
    }
    default:
      err = "`" + nm + "` must be numbers, strings or a list of them.";
      return false;
  }
}

// The result, protected once on the caller's count.
SEXP tor(const Node &n) {
  SEXP o = R_NilValue;
  switch (n.k) {
    case Node::NUL:
      return R_NilValue;
    case Node::NUM:
      o = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(n.d.size())));
      for (std::size_t i = 0; i < n.d.size(); ++i) REAL(o)[i] = std::isnan(n.d[i]) ? NA_REAL : n.d[i];
      break;
    case Node::INT:
      o = PROTECT(Rf_allocVector(INTSXP, static_cast<R_xlen_t>(n.d.size())));
      for (std::size_t i = 0; i < n.d.size(); ++i) {
        INTEGER(o)[i] = std::isnan(n.d[i]) ? NA_INTEGER : static_cast<int>(n.d[i]);
      }
      break;
    case Node::LGL:
      o = PROTECT(Rf_allocVector(LGLSXP, static_cast<R_xlen_t>(n.d.size())));
      for (std::size_t i = 0; i < n.d.size(); ++i) LOGICAL(o)[i] = std::isnan(n.d[i]) ? NA_LOGICAL : n.d[i] != 0;
      break;
    case Node::STR:
      o = PROTECT(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(n.s.size())));
      for (std::size_t i = 0; i < n.s.size(); ++i) SET_STRING_ELT(o, static_cast<R_xlen_t>(i), Rf_mkCharCE(n.s[i].c_str(), CE_UTF8));
      break;
    case Node::LIST: {
      o = PROTECT(Rf_allocVector(VECSXP, static_cast<R_xlen_t>(n.l.size())));
      for (std::size_t i = 0; i < n.l.size(); ++i) {
        SEXP e = PROTECT(tor(n.l[i]));
        SET_VECTOR_ELT(o, static_cast<R_xlen_t>(i), e);
        UNPROTECT(1);
      }
      if (!n.nm.empty()) {
        SEXP nm = PROTECT(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(n.nm.size())));
        for (std::size_t i = 0; i < n.nm.size(); ++i) SET_STRING_ELT(nm, static_cast<R_xlen_t>(i), Rf_mkChar(n.nm[i].c_str()));
        Rf_setAttrib(o, R_NamesSymbol, nm);
        UNPROTECT(1);
      }
      UNPROTECT(1);
      return o;
    }
  }
  if (!n.dim.empty()) {
    SEXP dm = PROTECT(Rf_allocVector(INTSXP, static_cast<R_xlen_t>(n.dim.size())));
    for (std::size_t i = 0; i < n.dim.size(); ++i) INTEGER(dm)[i] = n.dim[i];
    Rf_setAttrib(o, R_DimSymbol, dm);
    UNPROTECT(1);
  }
  UNPROTECT(1);
  return o;
}

// The kernel's run with every C++ object gone before an R error.
int go(const char *name, SEXP args, SEXP &res, char *msg, std::size_t cap) {
  std::string err;
  int rc = 0;
  {
    Node a = svk::list(), out;
    SEXP nms = Rf_getAttrib(args, R_NamesSymbol);
    for (R_xlen_t i = 0; i < XLENGTH(args) && rc == 0; ++i) {
      Node e;
      const std::string nm = Rf_translateCharUTF8(STRING_ELT(nms, i));
      if (!fromr(VECTOR_ELT(args, i), e, err, nm)) {
        rc = 2;
        break;
      }
      a.add(nm.c_str(), std::move(e));
    }
    if (rc == 0 && !svk::run(name, a, out, err)) rc = err.empty() ? 1 : 2;
    if (rc == 0) res = PROTECT(tor(out));
  }
  std::snprintf(msg, cap, "%s", err.c_str());
  return rc;
}

}  // namespace

extern "C" {

SEXP C_clab_r_svk(SEXP fn, SEXP args) {
  static char who[64];
  std::snprintf(who, sizeof who, "%s", Rf_translateCharUTF8(STRING_ELT(fn, 0)));
  clab::who = who;
  char msg[1024];
  SEXP res = R_NilValue;
  const int rc = go(who, args, res, msg, sizeof msg);
  if (rc == 1) clab::ichk(CLAB_EINTR);
  if (rc == 2) clab::fail("%s", msg);
  UNPROTECT(1);
  return res;
}

}  // extern "C"
