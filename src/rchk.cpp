// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The argument checks and result builders of the R binding. Each check
// keeps the wording and the order of tests of the R helper it replaced
// (.num, .one, .p01, .whole, .lgl, .ref); the R expression is quoted
// above each one.
#include "rchk.h"

#include <R_ext/Parse.h>
#define R_NO_REMAP_RMATH
#include <Rmath.h>

#include <climits>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>

#include "crimlab.h"

namespace clab {

const char *who = "crimlab";

//' @srrstats {G5.2a} Every message goes through fail() or kchk(), which put
//'   the raising function's name (clab::who) in front: no two functions
//'   raise the same message. errs() in tests/testthat/helper.R checks it.
//' @noRd
void fail(const char *f, ...) {
  // The message is formatted into R's transient memory (first 256 bytes,
  // then exactly the size needed) and handed to R as stop() hands it: R
  // truncates a long one the same way.
  // (f is never null; saying so keeps gcc's -fsanitize=undefined
  // instrumentation from inventing a null-format path)
  if (f == nullptr) f = "";
  std::size_t len = 256;
  char *buf = R_alloc(len, 1);
  std::va_list ap;
  va_start(ap, f);
  const int n = std::vsnprintf(buf, len, f, ap);
  va_end(ap);
  if (n >= 0 && static_cast<std::size_t>(n) >= len) {
    len = static_cast<std::size_t>(n) + 1;
    buf = R_alloc(len, 1);
    va_start(ap, f);
    std::vsnprintf(buf, len, f, ap);
    va_end(ap);
  }
  Rf_errorcall(R_NilValue, "%s: %s", who, buf);
}

SEXP lang(const char *src) {
  ParseStatus st;
  SEXP s = PROTECT(Rf_mkString(src));
  SEXP e = PROTECT(R_ParseVector(s, 1, &st, R_NilValue));
  SEXP c = VECTOR_ELT(e, 0);
  UNPROTECT(2);
  return c;
}

void kchk(int code, const char *src) {
  if (code == CLAB_OK) return;
  // Rf_error() gives the call of the R function that made the .Call, here
  // the exported one, as before. Where the old .Call sat in a helper, src
  // is that helper's call, parsed only on this error path.
  if (src == nullptr) Rf_error("%s: %s", who, clab_strerror(code));
  Rf_errorcall(PROTECT(lang(src)), "%s: %s", who, clab_strerror(code));
}

void fmt(double x, char *buf, int n) {
  // paste() turns a double into text with as.character(), which is
  // coerceVector(): 15 significant digits, "1e+05" for 100000.
  SEXP v = PROTECT(Rf_ScalarReal(x));
  SEXP s = PROTECT(Rf_coerceVector(v, STRSXP));
  std::snprintf(buf, static_cast<std::size_t>(n), "%s", CHAR(STRING_ELT(s, 0)));
  UNPROTECT(2);
}

void lfmt(double x, char *buf, int n) {
  // length() is an integer up to INT_MAX and a double above it
  if (x <= INT_MAX) {
    std::snprintf(buf, static_cast<std::size_t>(n), "%d", static_cast<int>(x));
  } else {
    fmt(x, buf, n);
  }
}

int ilen(double n, const char *nm) {
  if (n > INT_MAX) {
    if (nm == nullptr) Rf_error("%s: %s", who, clab_strerror(CLAB_EARG));
    fail("`%s` is too long for crimlab: more than %d values.", nm, INT_MAX);
  }
  return static_cast<int>(n);
}

double rd(double x) {
  // a volatile store rounds x to double: a product is never fused into the
  // sum that follows it, whatever the compiler and the target.
  volatile double v = x;
  return v;
}

//' @srrstats {G2.1} Types are asserted: isnum() is is.numeric(), and lgl(),
//'   the genotype and name checks of bind.cpp test their types likewise.
//' @srrstats {G2.6} A classed number (one with units, say) is read as its
//'   data; dates, times and factors are refused (test-bind.R).
//' @srrstats {EA2.6} Attributes and classes do not change results (G2.6).
//' @noRd
bool isnum(SEXP x) {
  // is.numeric(x): integer or double and not a factor. It dispatches, and
  // base's methods for Date, POSIXt and difftime say FALSE; another class
  // is read as its data (is.numeric() is TRUE for it unless the class has
  // its own method, which the binding does not consult).
  const int t = TYPEOF(x);
  if (t != INTSXP && t != REALSXP) return false;
  if (!Rf_isObject(x)) return true;
  return !Rf_inherits(x, "factor") && !Rf_inherits(x, "Date") &&
         !Rf_inherits(x, "POSIXt") && !Rf_inherits(x, "difftime");
}

bool anyna(SEXP x) {
  // anyNA(x): NaN counts; in a list, an element that is one NA value
  const R_xlen_t n = Rf_xlength(x);
  switch (TYPEOF(x)) {
    case REALSXP: {
      const double *p = REAL_RO(x);
      for (R_xlen_t i = 0; i < n; ++i) if (std::isnan(p[i])) return true;
      return false;
    }
    case INTSXP: {
      const int *p = INTEGER_RO(x);
      for (R_xlen_t i = 0; i < n; ++i) if (p[i] == NA_INTEGER) return true;
      return false;
    }
    case LGLSXP: {
      const int *p = LOGICAL_RO(x);
      for (R_xlen_t i = 0; i < n; ++i) if (p[i] == NA_LOGICAL) return true;
      return false;
    }
    case STRSXP:
      for (R_xlen_t i = 0; i < n; ++i) if (STRING_ELT(x, i) == NA_STRING) return true;
      return false;
    case VECSXP:
      for (R_xlen_t i = 0; i < n; ++i) {
        SEXP e = VECTOR_ELT(x, i);
        if (Rf_isVectorAtomic(e) && XLENGTH(e) == 1 && anyna(e)) return true;
      }
      return false;
    default:
      return false;
  }
}

//' @srrstats {G2.4} Inputs are converted explicitly once checked: dbl()
//'   (as.double()), whole() and ref() (as.integer()), labs()
//'   (as.character()).
//' @srrstats {G2.4b} dbl() is as.double() for integer and double input.
//' @noRd
SEXP dbl(SEXP x) {
  // as.double(x) for integer or double x; only values are read from it
  if (TYPEOF(x) == REALSXP) return x;
  const R_xlen_t n = XLENGTH(x);
  SEXP v = PROTECT(Rf_allocVector(REALSXP, n));
  const int *p = INTEGER_RO(x);
  double *o = REAL(v);
  for (R_xlen_t i = 0; i < n; ++i) o[i] = p[i] == NA_INTEGER ? NA_REAL : p[i];
  UNPROTECT(1);
  return v;
}

SEXP dcopy(SEXP x) {
  const R_xlen_t n = XLENGTH(x);
  SEXP v = PROTECT(Rf_allocVector(REALSXP, n));
  if (n > 0) std::memcpy(REAL(v), REAL_RO(x), static_cast<std::size_t>(n) * sizeof(double));
  UNPROTECT(1);
  return v;
}

//' @srrstats {G2.0} Lengths are asserted: num() with len, one(), p01() and
//'   whole() take one value, and paired vectors must match.
//' @srrstats {G2.13} Missing values are checked first: num() and the other
//'   checks refuse NA, naming the argument.
//' @srrstats {G2.15} No computation sees a missing value: every argument is
//'   checked before a sum, a mean or a kernel reads it.
//' @srrstats {G2.16} NaN, Inf and -Inf are refused here and again by the
//'   kernels (pos() and nonneg() in core.cpp call std::isfinite(), and
//'   range tests fail for NaN), which return CLAB_EDOM, not a number.
//' @noRd
//' @srrstats {G2.4e} A labelled vector (haven's labelled class, or any
//'   vector with a labels attribute) holds category codes: unlab() refuses
//'   it rather than read the codes as numbers or as group labels.
//' @noRd
void unlab(SEXP x, const char *nm) {
  if (!Rf_isNull(Rf_getAttrib(x, Rf_install("labels")))) {
    fail("`%s` has value labels: its numbers are category codes. Decode it "
         "first (haven::as_factor(), say).", nm);
  }
}

SEXP num(SEXP x, const char *nm, R_xlen_t len, double lo, bool open) {
  unlab(x, nm);
  // if (!is.numeric(x) || length(x) < 1L || anyNA(x) || !all(is.finite(x)))
  //   stop("`", nm, "` must be finite numbers with no missing values.")
  // if (!is.null(len) && length(x) != len)
  //   stop("`", nm, "` must have length ", len, ", not ", length(x), ".")
  // if (any(if (open) x <= lo else x < lo))
  //   stop("`", nm, "` must be ", "greater than " or "at least ", lo, ".")
  // as.double(x)
  if (!isnum(x) || XLENGTH(x) < 1) {
    fail("`%s` must be finite numbers with no missing values.", nm);
  }
  vec(x, nm);
  SEXP v = PROTECT(dbl(x));
  const double *p = REAL_RO(v);
  const R_xlen_t n = XLENGTH(v);
  for (R_xlen_t i = 0; i < n; ++i) {
    if (!std::isfinite(p[i])) {
      fail("`%s` must be finite numbers with no missing values.", nm);
    }
  }
  if (len >= 0 && n != len) {
    char a[32], b[32];
    lfmt(static_cast<double>(len), a, sizeof a);
    lfmt(static_cast<double>(n), b, sizeof b);
    fail("`%s` must have length %s, not %s.", nm, a, b);
  }
  for (R_xlen_t i = 0; i < n; ++i) {
    if (open ? p[i] <= lo : p[i] < lo) {
      char a[48];
      fmt(lo, a, sizeof a);
      fail("`%s` must be %s%s.", nm, open ? "greater than " : "at least ", a);
    }
  }
  UNPROTECT(1);
  return v;
}

//' @srrstats {G2.2} vec() refuses a matrix with more than one row and
//'   column where a vector is expected; num(), lgl() and the vector checks
//'   of bind.cpp call it (test-srr.R).
//' @noRd
void vec(SEXP x, const char *nm) {
  SEXP d = Rf_getAttrib(x, R_DimSymbol);
  if (TYPEOF(d) != INTSXP) return;
  int big = 0;
  for (R_xlen_t i = 0; i < XLENGTH(d); ++i) big += INTEGER_RO(d)[i] > 1;
  if (big > 1) {
    fail("`%s` must be a vector, not a matrix with more than one row and "
         "column.", nm);
  }
}

double val1(SEXP x) {
  if (!isnum(x) || XLENGTH(x) != 1) return NA_REAL;
  if (TYPEOF(x) == REALSXP) return REAL_RO(x)[0];
  const int v = INTEGER_RO(x)[0];
  return v == NA_INTEGER ? NA_REAL : v;
}

double one(SEXP x, const char *nm, double lo, double hi) {
  // !is.numeric(x) || length(x) != 1L || is.na(x) || x <= lo || x >= hi
  unlab(x, nm);
  const double v = val1(x);
  if (std::isnan(v) || v <= lo || v >= hi) {
    char a[48], b[48];
    fmt(lo, a, sizeof a);
    fmt(hi, b, sizeof b);
    fail("`%s` must be one number strictly between %s and %s.", nm, a, b);
  }
  return v;
}

double p01(SEXP x, const char *nm) {
  // !is.numeric(x) || length(x) != 1L || is.na(x) || x < 0 || x > 1
  unlab(x, nm);
  const double v = val1(x);
  if (std::isnan(v) || v < 0 || v > 1) {
    fail("`%s` must be one number from 0 to 1.", nm);
  }
  return v;
}

//' @srrstats {G2.4a} whole() and ref() return integers, as as.integer().
//' @noRd
int whole(SEXP x, const char *nm, int lo) {
  // !is.numeric(x) || length(x) != 1L || is.na(x) || x != round(x) ||
  //   x < lo || x > .Machine[["integer.max"]]
  // x == round(x) holds exactly when x is a whole number, as floor() says.
  unlab(x, nm);
  const double v = val1(x);
  if (std::isnan(v) || v != std::floor(v) || v < lo || v > INT_MAX) {
    fail("`%s` must be one whole number from %d.", nm, lo);
  }
  return static_cast<int>(v);
}

void lgl(SEXP x, const char *nm, R_xlen_t n) {
  // !is.logical(x) || length(x) != n || anyNA(x)
  unlab(x, nm);
  if (TYPEOF(x) != LGLSXP || XLENGTH(x) != n || anyna(x)) {
    fail("`%s` must be TRUE or FALSE, one per unit, with no missing values.",
         nm);
  }
  vec(x, nm);
}

int ref(SEXP r, SEXP nms, int n, const char *an, const char *what) {
  // if (is.character(ref) && length(ref) == 1L && !is.na(ref)) {
  //   i <- match(ref, nms)
  //   if (is.na(i)) stop("`ref` names no group: ", ref, ".")
  // }
  // is.numeric(ref) && length(ref) == 1L && !is.na(ref) &&
  //   ref == round(ref) && ref >= 1 && ref <= n
  // A factor is read by its label (its code would be a position).
  unlab(r, an);
  if (Rf_isFactor(r) && XLENGTH(r) == 1) r = labs(r);
  PROTECT(r);
  if (TYPEOF(r) == STRSXP && XLENGTH(r) == 1 && STRING_ELT(r, 0) != NA_STRING) {
    SEXP c = STRING_ELT(r, 0);
    const R_xlen_t m = Rf_xlength(nms);
    for (R_xlen_t i = 0; i < m; ++i) {
      if (seql(c, STRING_ELT(nms, i))) {
        UNPROTECT(1);
        return static_cast<int>(i + 1);
      }
    }
    fail("`%s` names no %s: %s.", an, what, Rf_translateChar(c));
  }
  UNPROTECT(1);
  const double v = val1(r);
  if (std::isnan(v) || v != std::floor(v) || v < 1 || v > n) {
    fail("`%s` must be one %s name or a position from 1 to %d.", an, what, n);
  }
  return static_cast<int>(v);
}

SEXP seqs(R_xlen_t n) {
  // as.character(seq_len(n)): integers (doubles past INT_MAX) coerced as
  // as.character() coerces them; R defers the conversion of a vector
  // without attributes, so no string is made until one is read
  const bool big = n > INT_MAX;
  SEXP v = PROTECT(Rf_allocVector(big ? REALSXP : INTSXP, n));
  for (R_xlen_t i = 0; i < n; ++i) {
    if (big) REAL(v)[i] = static_cast<double>(i + 1); else INTEGER(v)[i] = static_cast<int>(i + 1);  // NOLINT
  }
  SEXP s = Rf_coerceVector(v, STRSXP);
  UNPROTECT(1);
  return s;
}

double nsum(const double *x, R_xlen_t n) {
  // Neumaier's compensated sum, left to right; the compensation is added
  // only to a finite total (an infinite or NaN total is returned as is)
  double s = 0.0, c = 0.0;
  for (R_xlen_t i = 0; i < n; ++i) {
    const double t = s + x[i];
    c += std::fabs(s) >= std::fabs(x[i]) ? (s - t) + x[i] : (x[i] - t) + s;
    s = t;
  }
  return std::isfinite(s) ? s + c : s;
}

double nmean(const double *x, R_xlen_t n) {
  // mean(x): the sum over n, then, for a finite mean, the mean of the
  // residuals added back, as R does. R rescales a sum that overflows; the
  // terms averaged here are logarithms (below 2000 in size) of at most
  // INT_MAX ratios, so the sum cannot overflow and that branch is not
  // needed.
  // The first sum is plain, left to right: the correction pass that
  // follows compensates it, as in R.
  const double dn = static_cast<double>(n);
  double s = 0.0;
  for (R_xlen_t i = 0; i < n; ++i) s += x[i];
  s /= dn;
  if (!std::isfinite(s)) return s;
  double t = 0.0;
  for (R_xlen_t i = 0; i < n; ++i) t += x[i] - s;
  return s + t / dn;
}

double rpow(double x, double y) {
  // R's ^ on doubles: y == 2 squares, anything else is R_pow()
  return y == 2.0 ? x * x : R_pow(x, y);
}

Key key(SEXP c) {
  if (c == NA_STRING) return {nullptr, 0};
  if (Rf_getCharCE(c) == CE_BYTES) return {CHAR(c), 1};
  return {Rf_translateCharUTF8(c), 0};
}

int kcmp(const Key &a, const Key &b) {
  if (a.s == nullptr || b.s == nullptr) {
    return (a.s != nullptr) - (b.s != nullptr);
  }
  if (a.b != b.b) return a.b - b.b;
  return std::strcmp(a.s, b.s);
}

bool seql(SEXP a, SEXP b) {
  // R's Seql(): the same CHARSXP, or the same text in UTF-8; NA equals only
  // NA; a "bytes" string equals only a "bytes" string of the same bytes
  if (a == b) return true;
  if (a == NA_STRING || b == NA_STRING) return false;
  // two strings of one encoding are one CHARSXP when they are equal
  if (Rf_getCharCE(a) == Rf_getCharCE(b)) return false;
  const void *v = vmaxget();
  const bool eq = kcmp(key(a), key(b)) == 0;
  vmaxset(v);
  return eq;
}

SEXP labs(SEXP x) {
  // as.character(x): a factor gives its labels (levels(x)[x]); an atomic
  // vector or a list is coerced as as.character() coerces it. Date, POSIXt
  // and difftime were turned into text by their own methods; they are
  // refused (R_NilValue) rather than formatted differently.
  if (Rf_isFactor(x)) {
    SEXP lv = Rf_getAttrib(x, R_LevelsSymbol);
    const R_xlen_t n = XLENGTH(x), m = Rf_xlength(lv);
    const int *p = INTEGER_RO(x);
    SEXP out = PROTECT(Rf_allocVector(STRSXP, n));
    for (R_xlen_t i = 0; i < n; ++i) {
      const int k = p[i];
      SET_STRING_ELT(out, i, k == NA_INTEGER || k < 1 || k > m ? NA_STRING : STRING_ELT(lv, k - 1));
    }
    UNPROTECT(1);
    return out;
  }
  if (Rf_isObject(x) && (Rf_inherits(x, "Date") || Rf_inherits(x, "POSIXt") ||
                         Rf_inherits(x, "difftime"))) {
    return R_NilValue;
  }
  if (TYPEOF(x) == STRSXP) return x;
  return Rf_coerceVector(x, STRSXP);
}

SEXP frame(SEXP cols, const char *const *nm, bool late) {
  PROTECT(cols);
  const R_xlen_t k = XLENGTH(cols);
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, k));
  for (R_xlen_t i = 0; i < k; ++i) SET_STRING_ELT(nms, i, Rf_mkChar(nm[i]));
  Rf_setAttrib(cols, R_NamesSymbol, nms);
  // every column has passed ilen(), so -n fits an int
  const int n = static_cast<int>(XLENGTH(VECTOR_ELT(cols, 0)));
  SEXP rn = PROTECT(Rf_allocVector(INTSXP, n > 0 ? 2 : 0));
  if (n > 0) {
    INTEGER(rn)[0] = NA_INTEGER;
    INTEGER(rn)[1] = -n;
  }
  SEXP cl = PROTECT(Rf_mkString("data.frame"));
  if (late) Rf_setAttrib(cols, R_RowNamesSymbol, rn);
  Rf_classgets(cols, cl);
  if (!late) Rf_setAttrib(cols, R_RowNamesSymbol, rn);
  UNPROTECT(4);
  return cols;
}

SEXP nvec(int n, const char *const *nm, const double *v) {
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, n));
  for (int i = 0; i < n; ++i) {
    REAL(out)[i] = v[i];
    SET_STRING_ELT(nms, i, Rf_mkChar(nm[i]));
  }
  Rf_setAttrib(out, R_NamesSymbol, nms);
  UNPROTECT(2);
  return out;
}

SEXP nlist(SEXP x, const char *const *nm) {
  PROTECT(x);
  const R_xlen_t n = XLENGTH(x);
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, n));
  for (R_xlen_t i = 0; i < n; ++i) SET_STRING_ELT(nms, i, Rf_mkChar(nm[i]));
  Rf_setAttrib(x, R_NamesSymbol, nms);
  UNPROTECT(2);
  return x;
}

}  // namespace clab
