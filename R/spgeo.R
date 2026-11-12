#' Great-circle distance
#'
#' The distance along the sphere between pairs of points given in
#' longitude and latitude. The central angle is computed as the `atan2` of
#' the cross and dot products of the two unit vectors, written so that
#' nothing cancels: accurate to a few units in the last place from
#' coincident points to antipodes (the spherical case of Vincenty 1975). The
#' Earth is taken as a sphere of mean radius `r`; distances on the
#' ellipsoid differ by up to about 0.5 per cent.
#'
#' @param from,to The two ends of each pair, in degrees: a matrix or data
#'   frame with columns `lon` and `lat` (an unnamed two-column matrix is
#'   read as longitude then latitude), or one point as `c(lon, lat)`. Both
#'   have the same number of points, or one of them has one point, which is
#'   paired with every point of the other. Longitudes from -180 to 180,
#'   latitudes from -90 to 90; no missing values.
#' @param r The radius of the sphere, positive; the distances come back in
#'   its unit. The default is the mean Earth radius in metres.
#' @return A numeric vector of distances, one per pair, in the order of the
#'   input.
#' @references Vincenty, T. (1975). Direct and inverse solutions of
#'   geodesics on the ellipsoid with application of nested equations.
#'   *Survey Review*, 23(176), 88-93. \doi{10.1179/sre.1975.23.176.88}
#' @seealso [eqproj()] to project longitude and latitude to metres.
#' @examples
#' # Toronto to Montreal, metres
#' hvdist(c(-79.38, 43.65), c(-73.57, 45.50))
#' pts <- data.frame(lon = c(0, 90, 180), lat = c(0, 0, 0))
#' hvdist(pts, c(0, 0))
#' @export
hvdist <- function(from, to, r = 6371008.8) {
  .Call(C_clab_r_hvdist, from, to, r)
}

#' Equal-area projection of longitude and latitude, and its inverse
#'
#' `eqproj()` maps longitude and latitude to planar metres by the Lambert
#' azimuthal equal-area projection of the sphere about a centre; `eqinv()`
#' maps them back. Areas are exact, so intensities, nearest-neighbour
#' expectations and K functions computed on the result are right; lengths
#' are scaled by `cos(c / 2)` along the radius from the centre and its
#' inverse across it (`c` the central angle), so within 100 km of the
#' centre every length is right to 3.1e-5 relative. Every planar function of
#' the package takes the projected coordinates; those with a `lonlat`
#' argument call this projection themselves.
#'
#' @param x,y The points: `x` and `y` numeric vectors, or `x` a matrix or
#'   data frame with columns named `lon` and `lat` for `eqproj()` and `x`
#'   and `y` for `eqinv()` (an unnamed two-column matrix is read in that
#'   order) and `y` `NULL`. Longitudes from -180 to 180, latitudes from -90
#'   to 90; no missing values.
#' @param lon0,lat0 The centre of the projection in degrees. For
#'   `eqproj()`, `NULL` (both) takes the circular mean of the longitudes and
#'   the mean of the latitudes.
#' @inheritParams hvdist
#' @return A data frame with one row per point, in input order: `x` and `y`
#'   (east and north, in the unit of `r`) from `eqproj()`, with the centre
#'   used as the attribute `centre`; `lon` and `lat` from `eqinv()`.
#' @references Snyder, J. P. (1987). *Map Projections: A Working Manual*.
#'   U.S. Geological Survey Professional Paper 1395, chapter 24.
#'   \doi{10.3133/pp1395}
#' @examples
#' p <- eqproj(c(-79.40, -79.38, -79.35), c(43.64, 43.66, 43.65))
#' p
#' a <- attr(p, "centre")
#' eqinv(p, lon0 = a[1], lat0 = a[2])
#' @export
eqproj <- function(x, y = NULL, lon0 = NULL, lat0 = NULL, r = 6371008.8) {
  .Call(C_clab_r_eqproj, x, y, lon0, lat0, r)
}

#' @rdname eqproj
#' @export
eqinv <- function(x, y = NULL, lon0, lat0, r = 6371008.8) {
  .Call(C_clab_r_eqinv, x, y, lon0, lat0, r)
}

#' Mean centre and dispersion of points
#'
#' The (weighted) mean centre of a point pattern, the standard deviations of
#' the two coordinates, the standard distance deviation (the root mean
#' squared distance from the mean centre), and the geometric and harmonic
#' means of the coordinates (Bachi 1957; Neft 1962). The mean centre is the
#' point that minimises the sum of squared distances to the points. The
#' dispersions divide by the total weight less `k`: `k = 1` makes the
#' squared values unbiased for the variances; `k = 2` is the convention of
#' some crime-mapping texts; `k = 0` divides by the total weight.
#'
#' @param x,y The points: `x` and `y` numeric vectors of equal length; or
#'   `x` a matrix or data frame with columns named `x` and `y` (or `X` and
#'   `Y`, as [sf::st_coordinates()] names them, or `lon` and `lat` with
#'   `lonlat = TRUE`) and `y` `NULL`; or `x` an sf or sfc object of POINT
#'   geometries, whose coordinate reference system must agree with
#'   `lonlat`. An unnamed two-column matrix is read as `x` then `y`; other
#'   columns are never taken by position. Vectors given separately (`x`,
#'   `y`, and `w` or `t` where a function takes them) are matched by name
#'   when every one is named and no name repeats, as elsewhere in the
#'   package; otherwise by position. No missing or infinite values.
#' @param w `NULL` for equal weights, or one non-negative weight per point
#'   (counts at a location, say), with a positive sum.
#' @param k The divisor offset: 0, 1 or 2; the total weight must exceed it.
#' @param lonlat `FALSE` (the default) when `x` and `y` are planar
#'   coordinates in one length unit; `TRUE` when they are longitude and
#'   latitude in degrees, which are then projected to metres by [eqproj()]
#'   about their own centre, and centres are returned in degrees. It is
#'   never guessed.
#' @return A named numeric vector: `w` (the total weight), `mx`, `my` (the
#'   mean centre, in degrees when `lonlat = TRUE`), `sx`, `sy` (standard
#'   deviations of the coordinates), `sdd` (the standard distance
#'   deviation), `gx`, `gy` (geometric means) and `hx`, `hy` (harmonic
#'   means). The last four are `NaN` for a coordinate with a non-positive
#'   value under positive weight, and for `lonlat = TRUE` (they depend on the
#'   origin of the coordinates). Lengths are in the unit of the coordinates,
#'   or metres.
#' @references Bachi, R. (1957). Statistical analysis of geographical
#'   series. *Bulletin de l'Institut International de Statistique*, 36(2),
#'   229-240.
#'
#'   Neft, D. S. (1966). *Statistical Analysis for Areal Distributions*.
#'   Regional Science Research Institute, Philadelphia.
#' @examples
#' set.seed(1)
#' p <- cbind(x = rnorm(50, 100, 10), y = rnorm(50, 200, 5))
#' cgstat(p)
#' cgstat(p[, 1], p[, 2], w = rep(1:2, 25))
#' cgstat(p, k = 0)  # dispersions over n, not n - 1
#' @export
cgstat <- function(x, y = NULL, w = NULL, k = 1, lonlat = FALSE) {
  .Call(C_clab_r_cgstat, x, y, w, k, lonlat)
}

#' Median centre of points
#'
#' The median of the `x` coordinates and the median of the `y`
#' coordinates (the mean of the two middle values for an even number of
#' points). It depends on the orientation of the axes: rotating the points
#' does not rotate it. [webpt()] gives the rotation-invariant centre of
#' minimum distance.
#'
#' @inheritParams cgstat
#' @return A named numeric vector `c(x, y)`, in degrees when
#'   `lonlat = TRUE` (the medians are taken of the projected coordinates).
#' @references Neft, D. S. (1966). *Statistical Analysis for Areal
#'   Distributions*. Regional Science Research Institute, Philadelphia.
#' @examples
#' cgmedn(c(1, 2, 3, 10), c(5, 1, 2, 2))
#' @export
cgmedn <- function(x, y = NULL, lonlat = FALSE) {
  .Call(C_clab_r_cgmedn, x, y, lonlat)
}

#' Centre of minimum distance
#'
#' The point `c` that minimises the weighted sum of distances
#' `sum(w * |p - c|)` to the points: the Weber point or spatial median,
#' which, unlike the median centre, does not depend on the axes. It is
#' found by Weiszfeld's iteration with the modification of Vardi and Zhang
#' (2000), which keeps every step defined and the sum non-increasing when
#' an iterate lands on a data point. The optimality gap certifies the
#' answer: it is 0 exactly at a minimiser. The minimiser is unique unless
#' all the points lie on one line with the weight split evenly.
#'
#' @inheritParams cgstat
#' @param start `NULL` to start at the mean centre, or `c(x, y)` (in
#'   degrees when `lonlat = TRUE`).
#' @param tol Stop when a step is at most `tol` times the root mean squared
#'   distance to the mean centre: one positive number.
#' @param maxit The largest number of steps: one whole number from 0.
#' @return A named numeric vector: `x`, `y` (the centre, in degrees when
#'   `lonlat = TRUE`), `obj` (the weighted sum of distances), `gap` (the
#'   optimality gap, 0 at a minimiser), `iter` (steps taken) and `code`: 2
#'   at an exact minimiser, 1 when the step fell below `tol`, 0 when
#'   `maxit` was reached.
#' @references Weiszfeld, E. and Plastria, F. (2009). On the point for
#'   which the sum of the distances to n given points is minimum. *Annals of
#'   Operations Research*, 167, 7-41 (translation of the 1937 paper).
#'   \doi{10.1007/s10479-008-0352-z}
#'
#'   Vardi, Y. and Zhang, C.-H. (2000). The multivariate L1-median and
#'   associated data depth. *Proceedings of the National Academy of
#'   Sciences*, 97(4), 1423-1426. \doi{10.1073/pnas.97.4.1423}
#' @examples
#' webpt(c(0, 1, 0, 10), c(0, 0, 1, 10))
#' @export
webpt <- function(x, y = NULL, w = NULL, start = NULL, tol = 1e-10,
                  maxit = 1000, lonlat = FALSE) {
  .Call(C_clab_r_webpt, x, y, w, start, tol, maxit, lonlat)
}

#' Standard deviational ellipse
#'
#' The ellipse whose axes lie along the eigenvectors of the weighted
#' covariance matrix of the points about their mean centre, with semi-axes
#' `c` times the square roots of its eigenvalues (Wang, Shi and Miao 2015).
#' It is the ellipse of the rotated-axes formulas of Lefever (1926) and
#' Ebdon (1988) once the same divisor is used; texts differ only in the
#' divisor (`k`) and in the multiplier (some multiply by `sqrt(2)`). The
#' curve of the standard deviation by direction that Lefever drew is not
#' this ellipse (Gong 2002); the ellipse is returned. For bivariate normal
#' points the ellipse with `c = sqrt(-2 * log(1 - p))` covers a share `p` of
#' them: `c = 1` covers 39.3 per cent, `c = 2` 86.5 per cent.
#'
#' @inheritParams cgstat
#' @param c The multiplier of the standard deviations: one positive number.
#' @return A named numeric vector: `mx`, `my` (the mean centre, in degrees
#'   when `lonlat = TRUE`), `a` and `b` (the semi-major and semi-minor axes,
#'   `a >= b`, in the unit of the coordinates or metres), `angle` (of the
#'   major axis, in radians anticlockwise from the x axis (east), in
#'   `[0, pi)`; 0 when `a == b`) and `area` (`pi * a * b`).
#' @references Wang, B., Shi, W. and Miao, Z. (2015). Confidence analysis
#'   of standard deviational ellipse and its extension into higher
#'   dimensional Euclidean space. *PLoS ONE*, 10(3), e0118537.
#'   \doi{10.1371/journal.pone.0118537}
#'
#'   Gong, J. (2002). Clarifying the standard deviational ellipse.
#'   *Geographical Analysis*, 34(2), 155-167.
#'   \doi{10.1111/j.1538-4632.2002.tb01082.x}
#' @examples
#' set.seed(2)
#' x <- rnorm(200)
#' sdell(x, 0.5 * x + rnorm(200, sd = 0.3), c = 2)
#' @export
sdell <- function(x, y = NULL, w = NULL, k = 1, c = 1, lonlat = FALSE) {
  .Call(C_clab_r_sdell, x, y, w, k, c, lonlat)
}

#' Convex hull of points
#'
#' The smallest convex polygon holding every point, by Andrew's monotone
#' chain with an exact orientation predicate (Shewchuk 1997), so
#' near-collinear points never make a wrong turn. Points on an edge are not
#' vertices; of coincident points the first is reported.
#'
#' @inheritParams cgstat
#' @return A list: `idx`, the row numbers of the vertices, anticlockwise
#'   from the point with the lowest `x` (then lowest `y`); `area` and
#'   `perim`, the hull's area and perimeter (in the square and the unit of
#'   the coordinates, or square metres and metres when `lonlat = TRUE`). One
#'   distinct point gives one vertex and zeros; collinear points give their
#'   two ends, area 0 and twice the segment's length. Each coordinate must
#'   be 0 or of magnitude from 1e-140 to 1e150, where the predicate is
#'   exact.
#' @references Andrew, A. M. (1979). Another efficient algorithm for convex
#'   hulls in two dimensions. *Information Processing Letters*, 9(5),
#'   216-219. \doi{10.1016/0020-0190(79)90072-3}
#'
#'   Shewchuk, J. R. (1997). Adaptive precision floating-point arithmetic
#'   and fast robust geometric predicates. *Discrete and Computational
#'   Geometry*, 18(3), 305-363. \doi{10.1007/PL00009321}
#' @examples
#' cvhull(c(0, 1, 1, 0, 0.5), c(0, 0, 1, 1, 0.5))
#' @export
cvhull <- function(x, y = NULL, lonlat = FALSE) {
  .Call(C_clab_r_cvhull, x, y, lonlat)
}
