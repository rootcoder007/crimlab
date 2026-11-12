// SPDX-License-Identifier: AGPL-3.0-or-later
//
// Helpers of the R binding (bind.cpp): the argument checks the R wrappers
// made, word for word, and the builders of R results. Every check stops
// with Rf_errorcall(R_NilValue, ...), the condition stop(..., call. = FALSE)
// gives. Nothing here calls back into the R evaluator.
#ifndef CLAB_RCHK_H
#define CLAB_RCHK_H

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#include <limits>

namespace clab {

constexpr double ninf = -std::numeric_limits<double>::infinity();

// stop(<message>, call. = FALSE).
[[noreturn]] void fail(const char *fmt, ...);
// The R call in src, parsed (for conditions rebuilt on error paths).
SEXP lang(const char *src);
// A kernel's error code, raised with the call the old R wrapper made, src.
void kchk(int code, const char *src);
// as.character(x) for one double, into a buffer of n bytes.
void fmt(double x, char *buf, int n);
// A length as R prints it, into a buffer of n bytes.
void lfmt(R_xlen_t x, char *buf, int n);
// A length the kernels can take: at most INT_MAX.
int ilen(R_xlen_t n, const char *nm);
// x, rounded to a double: no fused multiply-add across it.
double rd(double x);

// is.numeric(x).
bool isnum(SEXP x);
// anyNA(x) for an atomic vector or a list.
bool anyna(SEXP x);
// as.double(x) of a numeric x: x itself when it is double already.
SEXP dbl(SEXP x);
// A fresh double copy of a double vector without its attributes.
SEXP dcopy(SEXP x);

// .num(): finite numbers, of length len (unless len < 0), all at least lo
// (above lo when open). Returns as.double(x).
SEXP num(SEXP x, const char *nm, R_xlen_t len = -1, double lo = ninf,
         bool open = false);
// x[1] as a double when is.numeric(x) and length(x) == 1, NA otherwise.
double val1(SEXP x);
// .one(): one number strictly between lo and hi.
double one(SEXP x, const char *nm, double lo, double hi);
// .p01(): one number from 0 to 1.
double p01(SEXP x, const char *nm);
// .whole(): one whole number from lo up to .Machine$integer.max.
int whole(SEXP x, const char *nm, int lo);
// .lgl(): TRUE or FALSE, n of them, no NA.
void lgl(SEXP x, const char *nm, R_xlen_t n);
// .ref(): the group named or numbered by ref, one-based.
int ref(SEXP r, SEXP nms, int n);
// as.character(seq_len(n)).
SEXP seqs(R_xlen_t n);
// .grp(): names(x), or "1", "2", ... for each value of x.
SEXP grp(SEXP x);

// Sums in a fixed order, compensated (Neumaier) while finite: the same on
// every platform. R's sum() accumulates in long double where the platform
// has one, so the old values could differ from these in the last bits.
double nsum(const double *x, R_xlen_t n);
// mean(): the sum over n, then one correction pass, as R does (in double).
double nmean(const double *x, R_xlen_t n);
// 10^y as R's ^ gives it (y == 2 is x * x, otherwise R_pow).
double rpow(double x, double y);

// A string as a sort key: the text in UTF-8 (bytes for "bytes" strings,
// which equal only "bytes" strings), NA before all others.
struct Key {
  const char *s;  // nullptr for NA
  int b;          // 1 for a "bytes" string
};
Key key(SEXP c);
int kcmp(const Key &a, const Key &b);
// match() and == on two strings.
bool seql(SEXP a, SEXP b);
// as.character(x) of atomic or list x; R_NilValue for a classed object the
// old code converted only through an S3 method (Date, POSIXt, difftime).
SEXP labs(SEXP x);

// A data frame from the list of columns cols, named nm. Attributes are set
// in the order data.frame() sets them (names, class, row.names), or in the
// order a column added later leaves them (names, row.names, class) when
// late; row names are the compact c(NA, -n).
SEXP frame(SEXP cols, const char *const *nm, bool late = false);
// A named double vector.
SEXP nvec(int n, const char *const *nm, const double *v);
// A named list of the n elements already in x.
SEXP nlist(SEXP x, const char *const *nm);

}  // namespace clab

#endif
