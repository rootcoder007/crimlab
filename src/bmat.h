// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The reader of a numeric matrix argument with one row per row of `y`,
// for the R entry points of the environmental-crime and dynamic-panel
// families (benv.cpp, bdpn.cpp). R-side only: kernels never include it.
#ifndef CLAB_BMAT_H
#define CLAB_BMAT_H

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#include <string>
#include <vector>

#include "dpnd.h"

namespace bmat {

// The columns of a numeric matrix (or one vector) of n rows, finite and
// at least lo, with their names (pre1, pre2, ... where unnamed).
using mat = dpnd::xm;
mat readm(SEXP x, const char *an, const char *pre, int n, double lo);

}  // namespace bmat

#endif
