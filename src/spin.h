// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Helpers of the spatial R entry points (bsp.cpp): reading coordinates in
// their accepted forms, longitude and latitude by request only, seeds,
// choices among named options. Every check stops through clab::fail(), so
// messages carry the function's name. Nothing here calls back into the R
// evaluator.
#ifndef CLAB_SPIN_H
#define CLAB_SPIN_H

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

namespace csp {

// Points read by xy(): projected coordinates (planar units), n of them, in
// input order, and how they were obtained.
struct pts {
  SEXP v;          // list(x, y) of doubles, protected by the caller
  const double *x;
  const double *y;
  int n;
  bool ll;         // the input was longitude and latitude in degrees
  double lon0;     // the centre of the projection when ll
  double lat0;
};

// The coordinates of a spatial function's first argument(s). x is a
// numeric vector with y another of the same length, or (y NULL) a matrix
// or data frame of two coordinate columns, matched by name: "x" and "y",
// or "lon" and "lat"; an unnamed two-column matrix is read as x then y.
// Other names, missing or non-finite values, and lon/lat columns without
// lonlat = TRUE are refused. ll says whether the coordinates are
// longitude and latitude in degrees (from lonlat = TRUE): they must then
// lie in [-180, 180] x [-90, 90], and are projected by clab_eqproj about
// their own centre (the circular mean of the longitudes, the mean of the
// latitudes), so the planar kernels receive metres. The result's v is
// protected on return: the caller unprotects it. Mode proj: as above. planar: the function takes projected coordinates
// only (lon/lat columns are refused; lonlat is not read). raw: the input
// is longitude and latitude (lonlat is not read), range-checked, not
// projected; columns named lon and lat or two unnamed ones, and with y
// NULL a plain vector of two numbers is one point. The centre lon0, lat0
// is set whenever the input is lon/lat.
enum mode { proj, planar, raw };
pts xy(SEXP x, SEXP y, SEXP lonlat, mode m = proj);

// Forward to the projection of p (clab_eqproj about its centre) for the n
// points (px, py) in place, when p.ll; nothing otherwise.
void fwd(const pts &p, double *px, double *py, int n);

// Back to longitude and latitude (clab_eqinv) for the n points (px, py) in
// place, when p.ll; nothing otherwise.
void back(const pts &p, double *px, double *py, int n);

// Per-unit vectors matched by name, the rule of the rest of the package
// (bind.cpp): when each of the k arguments v given as a plain vector (NULL
// ones are skipped) has names, all of one length, and no name repeats
// within one, the later ones are put in the order of the first by name,
// and names that differ are an error listing them; otherwise (some
// unnamed, a repeated name, a matrix or data frame) positions are used.
// Replaced vectors are protected: the result is how many to unprotect.
int byname(SEXP *v, const char *const *nm, int k);

// A per-group or per-unit argument in its other forms, as a named numeric
// vector: a list of single numbers or a one-row data frame, or a one-row
// or one-column numeric matrix (named by its column or row names). Any
// other input is returned as it is. The result is not protected.
SEXP gvec(SEXP x, const char *nm);

// TRUE or FALSE, one value.
bool flag(SEXP x, const char *nm);

// The seed for a kernel's generator: from R's generator (set.seed()) when
// s is NULL, otherwise one whole number from 0.
unsigned int seed(SEXP s);

// The position of the one string in s among the n options of opt.
int pick(SEXP s, const char *nm, const char *const *opt, int n);

// A rectangle c(xmin, xmax, ymin, ymax), finite, with xmin < xmax and
// ymin < ymax, into w.
void win(SEXP s, double *w);

// TRUE when the point (x, y) lies in the rectangle w (closed).
inline bool inside(const double *w, double x, double y) {
  return x >= w[0] && x <= w[1] && y >= w[2] && y <= w[3];
}

}  // namespace csp

#endif
