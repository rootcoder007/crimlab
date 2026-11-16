// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The reader of a numeric matrix argument (bmat.h).
#include "bmat.h"

#include <cstddef>
#include <string>

#include "rchk.h"

using clab::fail;

namespace bmat {

namespace {
inline std::size_t z(int i) { return static_cast<std::size_t>(i); }
}  // namespace

// The columns of a numeric matrix (or one vector) of n rows, finite and
// at least lo, with their names (pre1, pre2, ... where unnamed).
mat readm(SEXP x, const char *an, const char *pre, int n, double lo) {
  mat m;
  if (!clab::isnum(x)) fail("`%s` must be a numeric matrix or vector, one row per row of `y`.", an);
  const bool ism = Rf_isMatrix(x);
  const int nr = ism ? Rf_nrows(x) : static_cast<int>(XLENGTH(x));
  m.p = ism ? Rf_ncols(x) : 1;
  if (nr != n || m.p < 1) fail("`%s` must have one row per row of `y`, %d.", an, n);
  // num() takes a vector and refuses a matrix: it is given the values
  // alone, after the check for value labels that it would make on x
  clab::unlab(x, an);
  SEXP xd = PROTECT(clab::dbl(x));
  SEXP xv = PROTECT(clab::dcopy(xd));
  SEXP v = PROTECT(clab::num(xv, an, static_cast<R_xlen_t>(n) * m.p, lo));
  m.v.assign(REAL(v), REAL(v) + z(n) * z(m.p));
  SEXP dn = ism ? Rf_getAttrib(x, R_DimNamesSymbol) : R_NilValue;
  SEXP cn = Rf_isNull(dn) ? R_NilValue : VECTOR_ELT(dn, 1);
  for (int j = 0; j < m.p; ++j) {
    const bool named = !Rf_isNull(cn) && STRING_ELT(cn, j) != NA_STRING &&
                       CHAR(STRING_ELT(cn, j))[0] != '\0';
    m.nm.push_back(named ? std::string(Rf_translateCharUTF8(STRING_ELT(cn, j)))
                         : pre + std::to_string(j + 1));
  }
  UNPROTECT(3);
  return m;
}

}  // namespace bmat
