// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Helpers of the policy-evaluation R entry points (bdid.cpp, bstg.cpp,
// bsyn.cpp, bwcb.cpp): reading a unit-by-period panel, adjustment
// columns, adoption periods, and building R results. Every check stops
// through clab::fail(), so messages carry the function's name. Nothing
// here calls back into the R evaluator.
#ifndef CLAB_PEVR_H
#define CLAB_PEVR_H

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#include <cstddef>
#include <initializer_list>
#include <string>
#include <utility>
#include <vector>

namespace pevr {

inline std::size_t z(int i) { return static_cast<std::size_t>(i); }

// NA for NaN.
double na(double v);

// Columns of a result, all protected until the frame is built.
struct cols {
  std::vector<std::string> nm;
  SEXP v = R_NilValue;
  int k = 0;
  explicit cols(int ncol);
  double *dbl(const char *name, int n);
  SEXP str(const char *name, int n);
  // the data frame; unprotects the list
  SEXP frame();
};

// A named list of the given elements (protected by the caller).
SEXP nlist(std::initializer_list<std::pair<const char *, SEXP>> e);

// One column of a design: its values and its term label.
struct term {
  std::vector<double> x;
  std::string nm;
};

// The adjustment columns: a data frame or list of columns, a numeric
// matrix (its columns), or one vector, n rows. Labels (factor, text,
// logical) become indicators of every level but the first; numbers enter
// as they are.
void adjcols(SEXP adj, int n, std::vector<term> &out);

// A panel: the outcome y (n), the unit codes u (0..nu-1, labels ulev), the
// period codes t (0..nt-1) and the periods' values tv, sorted.
struct pan {
  int n = 0;
  const double *y = nullptr;
  std::vector<int> u, t;
  int nu = 0, nt = 0;
  SEXP ulev = R_NilValue;
  std::vector<double> tv;
};

// The panel of y (finite numbers, at least lo), unit (labels) and time
// (finite numbers), one row per unit and period. Protects what it keeps
// and adds the count to *np.
pan panel(SEXP y, SEXP unit, SEXP time, double lo, int *np);

// Each unit's first treated period (a code; nt for never) from the 0 or 1
// treatment d of each row, which must stay 1 once it is 1.
std::vector<int> firsts(const pan &pn, const int *d);

// Each unit's first treated period as a value from the column first (one
// per row, constant within a unit): NA, 0 and Inf mean never (Inf).
std::vector<double> firstcol(SEXP first, const pan &pn);

}  // namespace pevr

#endif
