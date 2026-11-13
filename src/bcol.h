// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Helpers of the R entry points of the custody, court and published-table
// families (bflw.cpp, bmnd.cpp, bcnc.cpp, bjdg.cpp, blcr.cpp, bdkf.cpp,
// bsev.cpp, bpub.cpp): the result builders and the reading of labels and
// numbers they share. Header only.
#ifndef CLAB_BCOL_H
#define CLAB_BCOL_H

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>

#include "crimlab.h"
#include "rchk.h"
#include "spin.h"

namespace bj {

// The columns of a data frame of n rows (at most 40). The list stays
// protected until the entry point's final UNPROTECT, which counts it once.
struct cols {
  SEXP v;
  const char *nm[40];
  int k = 0, n;
  cols(int ncol, int nr) : n(nr) { v = PROTECT(Rf_allocVector(VECSXP, ncol)); }
  void put(const char *name, SEXP c) {
    SET_VECTOR_ELT(v, k, c);
    nm[k++] = name;
  }
  double *dbl(const char *name) {
    SEXP c = Rf_allocVector(REALSXP, n);
    put(name, c);
    return REAL(c);
  }
  int *itg(const char *name) {
    SEXP c = Rf_allocVector(INTSXP, n);
    put(name, c);
    return INTEGER(c);
  }
  int *lgl(const char *name) {
    SEXP c = Rf_allocVector(LGLSXP, n);
    put(name, c);
    return LOGICAL(c);
  }
  // the frame of the k columns added
  SEXP done() {
    SEXP w = PROTECT(Rf_allocVector(VECSXP, k));
    for (int j = 0; j < k; ++j) SET_VECTOR_ELT(w, j, VECTOR_ELT(v, j));
    SEXP f = clab::frame(w, nm);
    UNPROTECT(1);
    return f;
  }
};

// A named list of the k elements of the protected list x.
inline SEXP named(SEXP x, const char *const *nm) { return clab::nlist(x, nm); }

inline int len(SEXP x, const char *nm) { return clab::ilen(static_cast<double>(XLENGTH(x)), nm); }

// Labels to codes from 0 (a factor by its levels that occur, anything else
// by its sorted values), refusing a labelled vector and a matrix. Protects
// two.
inline csp::lab labs(SEXP x, const char *nm, R_xlen_t n) {
  clab::unlab(x, nm);
  clab::vec(x, nm);
  return csp::codes(x, nm, n);
}

// Finite numbers from lo (above lo when open), at least one, refusing a
// matrix; n of them unless n < 0. Not protected.
inline SEXP nums(SEXP x, const char *nm, R_xlen_t n = -1, double lo = clab::ninf,
                 bool open = false) {
  return clab::num(x, nm, n, lo, open);
}

// One number from lo to hi, both included (hi Inf: from lo up).
inline double inr(SEXP x, const char *nm, double lo, double hi) {
  clab::unlab(x, nm);
  const double v = clab::val1(x);
  if (!(v >= lo && v <= hi) || !std::isfinite(v)) {
    char a[48], b[48];
    clab::fmt(lo, a, sizeof a);
    clab::fmt(hi, b, sizeof b);
    if (std::isfinite(hi)) clab::fail("`%s` must be one number from %s to %s.", nm, a, b);
    clab::fail("`%s` must be one number from %s.", nm, a);
  }
  return v;
}

// NaN as R's NA_real_ (a kernel's "no value").
inline double na(double v) { return std::isnan(v) ? NA_REAL : v; }

// A kernel's code: OK, or stop with msg.
inline void kc(int rc, const char *msg) {
  if (rc != CLAB_OK) clab::fail("%s", msg);
}

// One label of a character vector as a C string (UTF-8).
inline const char *str(SEXP s, R_xlen_t i) { return Rf_translateCharUTF8(STRING_ELT(s, i)); }

// Numbers from 0 given per label of lev (k labels) into out: named ones
// matched by name (names of labels not in lev are ignored, a label without
// a value or with two is refused), unnamed ones in the order of lev.
inline void perlev(SEXP v, const char *nm, SEXP lev, int k, double *out) {
  SEXP d = PROTECT(nums(v, nm, -1, 0.0));
  SEXP vn = Rf_getAttrib(v, R_NamesSymbol);
  const R_xlen_t m = XLENGTH(d);
  if (Rf_isNull(vn)) {
    if (m != k) clab::fail("`%s` must have one value per label, %d, or be named by label.", nm, k);
    for (int j = 0; j < k; ++j) out[j] = REAL(d)[j];
  } else {
    // the first label named twice (dup) or left without a value (none);
    // the map is gone before any error is raised
    R_xlen_t dup = -1;
    int none = -1;
    {
      std::unordered_map<std::string, int> at;
      for (int j = 0; j < k; ++j) at.emplace(str(lev, j), j);
      std::vector<int> seen(static_cast<std::size_t>(k), 0);
      for (R_xlen_t i = 0; i < m && dup < 0; ++i) {
        if (STRING_ELT(vn, i) == NA_STRING) continue;
        auto f = at.find(str(vn, i));
        if (f == at.end()) continue;
        if (seen[static_cast<std::size_t>(f->second)]++) dup = i;
        out[f->second] = REAL(d)[i];
      }
      for (int j = 0; j < k && none < 0 && dup < 0; ++j) {
        if (!seen[static_cast<std::size_t>(j)]) none = j;
      }
    }
    if (dup >= 0) clab::fail("`%s` names %s twice.", nm, str(vn, dup));
    if (none >= 0) clab::fail("`%s` has no value for %s.", nm, str(lev, none));
  }
  UNPROTECT(1);
}

}  // namespace bj

#endif
