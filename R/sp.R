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

#' Nearest-neighbour index
#'
#' The mean distance from each point to its k-th nearest neighbour,
#' against its expectation under complete spatial randomness (CSR) in the
#' rectangle `win`: the nearest-neighbour index of Clark and Evans (1954)
#' and its k-order extension (Thompson 1956), with Donnelly's (1978)
#' edge-corrected expectation for the first order and, when `nsim > 0`, a
#' Monte Carlo test against `nsim` patterns of as many uniform points in
#' the rectangle. The Monte Carlo p-values have exact size; the z
#' statistics ignore that nearest-neighbour distances are dependent and,
#' uncorrected, the edge. The index tests only the k-th neighbour
#' distances, not complete spatial randomness as a whole.
#'
#' @param x,y The points in planar coordinates (project longitude and
#'   latitude with [eqproj()] first): `x` and `y` numeric vectors of equal
#'   length, or `x` a matrix or data frame with columns named `x` and `y`
#'   (an unnamed two-column matrix is read as `x` then `y`). At least two
#'   points, all inside `win`; no missing values.
#' @param win The study rectangle `c(xmin, xmax, ymin, ymax)`, in the unit
#'   of the coordinates. The intensity is the number of points over its
#'   area.
#' @param kmax The highest neighbour order: one whole number from 1 to the
#'   number of points less 1.
#' @param nsim The number of simulated patterns for the Monte Carlo test:
#'   one whole number from 0 (no test).
#' @param seed `NULL` to draw the simulations' seed from R's generator (so
#'   [set.seed()] fixes them), or one whole number from 0.
#' @param dist `TRUE` to return every point's neighbour distances.
#' @param metric The distance: `"euclid"` or `"manhattan"` (`|dx| + |dy|`,
#'   a street grid aligned with the axes), as in [spwt()] and [jtsurf()].
#'   Under Manhattan distance the expectations use the area of the
#'   Manhattan ball, `2 r^2`: the mean k-th distance is
#'   `sqrt(pi / 2)` times the Euclidean expectation and `E R_k^2 = k / (2
#'   lambda)`; no edge-corrected (Donnelly) expectation is published for it,
#'   so `ce` holds `NaN` except `z`, about the uncorrected expectation with
#'   `Var R_1 = (4 - pi) / (8 lambda)`. The Monte Carlo test is unchanged.
#' @return A list. `order`, a data frame with one row per order `k`:
#'   `mean` (the mean k-th neighbour distance), `exp` (its expectation under
#'   CSR without edge effects, `c_k / sqrt(lambda)`), `ratio`, `pclu` and
#'   `preg` (the Monte Carlo p-values against clustering and against
#'   regularity, `NaN` without simulations) and `sd` (the standard deviation
#'   of one k-th distance under CSR). `ce`, a named vector for the first
#'   order: `donnelly` (the edge-corrected expectation), `ratio` (the index
#'   against it), `z` and `zc` (the Clark-Evans z statistics about the
#'   uncorrected and the corrected expectation). `dist`, `NULL` or a matrix
#'   with one row per point (input order) and one column per order.
#'   Distances are in the unit of the coordinates.
#' @references Clark, P. J. and Evans, F. C. (1954). Distance to nearest
#'   neighbor as a measure of spatial relationships in populations.
#'   *Ecology*, 35(4), 445-453. \doi{10.2307/1931034}
#'
#'   Thompson, H. R. (1956). Distribution of distance to nth neighbour in a
#'   population of randomly distributed individuals. *Ecology*, 37(2),
#'   391-394. \doi{10.2307/1933159}
#'
#'   Donnelly, K. P. (1978). Simulations to determine the variance and edge
#'   effect of total nearest-neighbour distance. In I. Hodder (ed.),
#'   *Simulation Methods in Archaeology*, 91-95. Cambridge University
#'   Press.
#' @examples
#' set.seed(1)
#' x <- runif(200)
#' y <- runif(200)
#' nnidx(x, y, win = c(0, 1, 0, 1), kmax = 3, nsim = 99)
#' @export
nnidx <- function(x, y = NULL, win, kmax = 1, nsim = 0, seed = NULL,
                  dist = FALSE, metric = "euclid") {
  .Call(C_clab_r_nnidx, x, y, win, kmax, nsim, seed, dist, metric)
}

#' Ripley's K function
#'
#' Ripley's K function of a point pattern in the rectangle `win` at radii
#' `r`, and Besag's `L = sqrt(K / pi)`. Under complete spatial randomness
#' `K(r) = pi r^2` and `L(r) = r`. Four estimators: `"none"` (no edge
#' correction), `"border"` (only points at least `r` from the edge are
#' centres), `"modif"` (the modified border estimator, with the area of the
#' eroded window; `NaN` once `r` reaches half a side) and `"iso"` (Ripley's
#' isotropic correction, each pair weighted by the inverse of the share of
#' the circle through the other point that lies in the rectangle, computed
#' exactly). The normalisation is `n (n - 1)`, as in spatstat.
#'
#' @inheritParams nnidx
#' @param r The radii: non-negative numbers in strictly increasing order,
#'   in the unit of the coordinates.
#' @param corr The edge correction: `"iso"`, `"border"`, `"modif"` or
#'   `"none"`.
#' @param metric The distance: `"euclid"` or `"manhattan"` (`|dx| + |dy|`,
#'   a street grid aligned with the axes), as in [spwt()] and [jtsurf()].
#'   Under Manhattan distance `K(r) = 2 r^2` and `L = sqrt(K / 2)` under
#'   randomness; the border and modified border
#'   corrections carry over (a Manhattan ball of radius `r` lies in the
#'   rectangle exactly when its centre is at least `r` from every edge),
#'   while the isotropic correction, the share of a circle, has no
#'   Manhattan form and is refused.
#' @return A data frame with one row per radius: `r`, `k` (in the square of
#'   the coordinates' unit) and `l` (`sqrt(k / pi)`, or `sqrt(k / 2)` under
#'   Manhattan distance).
#' @references Ripley, B. D. (1976). The second-order analysis of
#'   stationary point processes. *Journal of Applied Probability*, 13(2),
#'   255-266. \doi{10.2307/3212829}
#'
#'   Ripley, B. D. (1977). Modelling spatial patterns. *Journal of the
#'   Royal Statistical Society, Series B*, 39(2), 172-212.
#'   \doi{10.1111/j.2517-6161.1977.tb01615.x}
#' @examples
#' set.seed(1)
#' ripk(runif(100), runif(100), win = c(0, 1, 0, 1), r = seq(0, 0.25, 0.05))
#' @export
ripk <- function(x, y = NULL, win, r, corr = "iso", metric = "euclid") {
  .Call(C_clab_r_ripk, x, y, win, r, corr, metric)
}

#' Envelopes of Ripley's L function and a test of complete spatial
#' randomness
#'
#' `L(r) - r` for the data, its pointwise minimum, maximum and mean over
#' `nsim` patterns of as many uniform points in the rectangle, and a global
#' Monte Carlo test on the maximum absolute deviation
#' `T = max_r |L(r) - r|`, whose p-value `(1 + #{T_sim >= T}) / (nsim + 1)`
#' has exact size (Ripley 1977; Diggle 2013). The pointwise envelopes are
#' not a test of the whole curve: each radius alone falls outside them with
#' probability `2 / (nsim + 1)` under randomness.
#'
#' @inheritParams ripk
#' @inheritParams nnidx
#' @param nsim The number of simulated patterns: one whole number from 1.
#' @details Euclidean distance only; [ripk()] takes `metric`.
#' @return A list. `env`, a data frame with one row per radius: `r`, `obs`
#'   (`L(r) - r` of the data), `lo`, `hi` (the pointwise envelope) and
#'   `mean`. `test`, a named vector: `t` (the maximum absolute deviation)
#'   and `p` (its Monte Carlo p-value).
#' @references Diggle, P. J. (2013). *Statistical Analysis of Spatial and
#'   Spatio-Temporal Point Patterns* (3rd ed.), chapter 6. CRC Press.
#'
#'   Ripley, B. D. (1977). Modelling spatial patterns. *Journal of the
#'   Royal Statistical Society, Series B*, 39(2), 172-212.
#'   \doi{10.1111/j.2517-6161.1977.tb01615.x}
#' @examples
#' set.seed(1)
#' ripenv(runif(100), runif(100), win = c(0, 1, 0, 1),
#'        r = seq(0.01, 0.25, 0.02), nsim = 39)$test
#' @export
ripenv <- function(x, y = NULL, win, r, corr = "iso", nsim = 99,
                   seed = NULL) {
  .Call(C_clab_r_ripenv, x, y, win, r, corr, nsim, seed)
}

#' Fuzzy mode
#'
#' For each point, the number of points (itself and coincident ones
#' included) within distance `r` of it. With `r = 0` it counts coincident
#' points, the mode of exact locations; with `r > 0` it is the fuzzy mode
#' of the crime-mapping literature, a fixed-radius scan centred on each
#' incident. Counts of nearby points overlap by design: each point counts
#' its own neighbourhood.
#'
#' @inheritParams cgstat
#' @param r The radius: one number from 0, in the unit of the coordinates
#'   (metres when `lonlat = TRUE`). Distances up to and including `r`
#'   count.
#' @param count `NULL` for points, or for zones: one whole count from 0
#'   per zone (`x`, `y` the zones' centroids); the result is then each
#'   zone's total count within `r` of it, its own included.
#' @param metric The distance: `"euclid"` or `"manhattan"` (`|dx| + |dy|`),
#'   as in [spwt()] and [jtsurf()].
#' @return An integer vector, one count per point, in input order (with
#'   `count`, a numeric vector of zone totals).
#' @references Eck, J. E., Chainey, S., Cameron, J. G., Leitner, M. and
#'   Wilson, R. E. (2005). *Mapping Crime: Understanding Hot Spots*.
#'   National Institute of Justice, NCJ 209393.
#' @examples
#' fzmode(c(0, 0, 1, 5), c(0, 0, 0, 5), r = 1)
#' @export
fzmode <- function(x, y = NULL, r, count = NULL, metric = "euclid",
                   lonlat = FALSE) {
  .Call(C_clab_r_fzmode, x, y, r, count, metric, lonlat)
}

#' K-means partition of points
#'
#' Partitions the points into `k` groups by Lloyd's iterations: each point
#' joins its nearest centre (the lowest-numbered on ties), then each centre
#' moves to its points' mean (`upd = "mean"`, K-means, minimising the
#' within-group sum of squares) or to their centre of minimum distance
#' (`upd = "weber"`, a spatial K-medians minimising the sum of distances),
#' until no point changes group or after `maxit` moves. An emptied group
#' keeps its centre. The result depends on the start: k-means++ seeding
#' (Arthur and Vassilvitskii 2007), the centres of the busiest grid cells
#' at least a given distance apart, or centres given.
#'
#' @inheritParams cgstat
#' @param k The number of groups: one whole number from 1 to the number of
#'   points.
#' @param init The start: `"pp"` (k-means++, drawn from `seed`), `"cells"`
#'   (square cells of side `cell[1]` from the lower-left corner of the
#'   points' bounding box, taken by decreasing count, a cell's centre
#'   becoming a seed when it is at least `cell[2]` from every earlier one),
#'   or the `k` starting centres as a matrix or data frame of coordinates
#'   (degrees when `lonlat = TRUE`).
#' @param cell For `init = "cells"`: `c(side, separation)`, the side
#'   positive and the separation from 0, in the unit of the coordinates.
#' @param upd The centre update: `"mean"` or `"weber"`.
#' @param maxit The largest number of iterations: one whole number from 0.
#' @param count `NULL` for points, or for zones: one whole count from 0
#'   per zone, the points `x`, `y` being the zones' centroids. A zone is
#'   taken as `count` incidents at its centroid, which these hot-spot rules
#'   treat exactly as coincident points. Zones with count 0 get `NA`.
#' @inheritParams nnidx
#' @return A list. `cl`, the group of each point (1 to `k`, input order).
#'   `cen`, a data frame with one row per group: `x`, `y` (its centre, in
#'   degrees when `lonlat = TRUE`), `n` (points), `ss` (sum of squared
#'   distances to the centre) and `sd` (sum of distances). `fit`, a named
#'   vector: `ss` and `sd` (totals), `iter` (moves made) and `conv` (1 if no
#'   point changed group at the end).
#' @references Lloyd, S. (1982). Least squares quantization in PCM. *IEEE
#'   Transactions on Information Theory*, 28(2), 129-137.
#'   \doi{10.1109/TIT.1982.1056489}
#'
#'   Arthur, D. and Vassilvitskii, S. (2007). k-means++: the advantages of
#'   careful seeding. *Proceedings of the Eighteenth Annual ACM-SIAM
#'   Symposium on Discrete Algorithms*, 1027-1035.
#' @examples
#' set.seed(1)
#' x <- c(rnorm(30, 0), rnorm(30, 5))
#' y <- c(rnorm(30, 0), rnorm(30, 5))
#' kmpart(x, y, k = 2, seed = 1)$cen
#' @export
kmpart <- function(x, y = NULL, k, init = "pp", cell = NULL, upd = "mean",
                   maxit = 100, seed = NULL, count = NULL, lonlat = FALSE) {
  .Call(C_clab_r_kmpart, x, y, k, init, cell, upd, maxit, seed, count,
        lonlat)
}

#' Nearest-neighbour hierarchical clustering
#'
#' Groups points that lie closer together than a threshold distance, then
#' groups the groups, order by order (after Johnson 1967, as used for crime
#' hot spots). At the first order two points are linked when their distance
#' is below the threshold: `fixd`, or the mean nearest-neighbour distance
#' expected under complete spatial randomness in an area `area` plus `z`
#' standard errors, `0.5 sqrt(A / m) + z sqrt((4 - pi) A / (4 pi)) / m` for
#' `m` units, with `z = qnorm(p)` (`p = 0.5` gives the expectation itself,
#' smaller `p` a stricter threshold). Points with a linked neighbour are
#' taken in decreasing number of links; each one not yet clustered takes
#' every unclustered point linked to it; groups of at least `minpts` points
#' are kept. Their centres of minimum distance become the units of the next
#' order, with the threshold recomputed for their number, while at least
#' four remain, up to `maxord` orders. The iterative reassignment some
#' descriptions add after the first pass is not done: its rule is not stated
#' unambiguously.
#'
#' @inheritParams cgstat
#' @param area The study area for the random-expectation threshold, in the
#'   square of the coordinates' unit (square metres when `lonlat = TRUE`).
#'   Give `area` or `fixd`, not both.
#' @param fixd A fixed threshold distance instead, positive.
#' @param p The one-tailed probability that a pair is as close by chance:
#'   one number strictly between 0 and 1.
#' @param minpts,minhi The fewest points (order 1) and member units (higher
#'   orders) a cluster keeps: whole numbers from 2.
#' @param maxord The highest order: one whole number from 1 to 20.
#' @param c The multiplier of each cluster's standard deviational ellipse
#'   (see [sdell()]).
#' @param count `NULL` for points, or for zones: one whole count from 0
#'   per zone, the points `x`, `y` being the zones' centroids. A zone is
#'   taken as `count` incidents at its centroid, which these hot-spot rules
#'   treat exactly as coincident points. Zones with count 0 get `NA`; the
#'   random-expectation threshold counts incidents, not zones.
#' @return A list. `cl`, an integer matrix with one row per point (input
#'   order) and one column per order: the cluster (a row of `clusters`) the
#'   point belongs to at that order, or `NA`. `clusters`, a data frame with
#'   one row per cluster: `order`, `points`, `units` (member units),
#'   `mx`, `my` (mean centre), `wx`, `wy` (centre of minimum distance; both
#'   in degrees when `lonlat = TRUE`), `a`, `b`, `angle` (the ellipse of the
#'   member units, `NaN` when too few for the divisor) and `thresh` (the
#'   threshold distance used).
#' @references Johnson, S. C. (1967). Hierarchical clustering schemes.
#'   *Psychometrika*, 32(3), 241-254. \doi{10.1007/BF02289588}
#'
#'   Eck, J. E., Chainey, S., Cameron, J. G., Leitner, M. and Wilson, R. E.
#'   (2005). *Mapping Crime: Understanding Hot Spots*. National Institute
#'   of Justice, NCJ 209393.
#' @examples
#' set.seed(1)
#' x <- c(rnorm(40, 0, 0.3), rnorm(40, 6, 0.3), runif(40, -5, 10))
#' y <- c(rnorm(40, 0, 0.3), rnorm(40, 6, 0.3), runif(40, -5, 10))
#' nnhclu(x, y, area = 225, p = 0.05, minpts = 5)$clusters
#' @export
nnhclu <- function(x, y = NULL, area = NULL, fixd = NULL, p = 0.5,
                   minpts = 10, minhi = minpts, maxord = 3, k = 1, c = 1,
                   count = NULL, lonlat = FALSE) {
  .Call(C_clab_r_nnhclu, x, y, area, fixd, p, minpts, minhi, maxord, k, c,
        count, lonlat)
}

#' Search-circle hot-spot areas
#'
#' Hot-spot areas from search circles (Block and Block 1995): circles of
#' radius `r` are centred on the nodes of an `nx` by `ny` grid spanning the
#' rectangle `win` (corners included), the points in each are counted, the
#' circles holding at least two points are ranked by count and the first
#' `top` are kept; kept circles that share a point are merged, repeatedly,
#' and a hot-spot area is the set of points in a merged group, kept if it
#' has at least `minpts`. Descriptions differ on the circle's radius and on
#' the grid; both are the caller's choice here.
#'
#' @inheritParams nnidx
#' @param win The rectangle the grid spans, `c(xmin, xmax, ymin, ymax)`;
#'   points may lie outside it.
#' @param nx,ny The numbers of grid nodes across and up: whole numbers from
#'   2 to 10000.
#' @param r The search radius, positive, in the unit of the coordinates.
#' @param top The number of busiest circles kept: one whole number from 1.
#' @param minpts The fewest points a hot-spot area keeps: a whole number
#'   from 2.
#' @inheritParams nnhclu
#' @return A list. `cl`, the hot-spot area of each point (zone) (a row of
#'   `clusters`) or `NA`, in input order. `clusters`, a data frame with one
#'   row per area, by size: `points`, `mx`, `my` (mean centre), `wx`, `wy`
#'   (centre of minimum distance), `a`, `b` and `angle` (the standard
#'   deviational ellipse of its points).
#' @references Block, R. L. and Block, C. R. (1995). Space, place and
#'   crime: hot spot areas and hot places of liquor-related crime. In J. E.
#'   Eck and D. Weisburd (eds.), *Crime and Place*, Crime Prevention Studies
#'   4, 145-183. Criminal Justice Press.
#' @examples
#' set.seed(1)
#' x <- c(rnorm(30, 2, 0.2), runif(60, 0, 10))
#' y <- c(rnorm(30, 2, 0.2), runif(60, 0, 10))
#' stacel(x, y, win = c(0, 10, 0, 10), r = 0.75, top = 5, minpts = 5)$clusters
#' @export
stacel <- function(x, y = NULL, win, nx = 20, ny = 20, r, top = 25,
                   minpts = 5, k = 1, c = 1, count = NULL) {
  .Call(C_clab_r_stacel, x, y, win, nx, ny, r, top, minpts, k, c, count)
}

#' Mantel's space-time test
#'
#' The Pearson correlation, over all pairs of events, between their
#' spatial dissimilarity (the distance, or `1 / (distance + cs)` when
#' `cs > 0`, Mantel's reciprocal transform) and their temporal one
#' (`|t_i - t_j|`, or `1 / (|t_i - t_j| + ct)`), with a permutation test
#' that permutes the times over the fixed locations (Mantel 1967). Drawing
#' distances and intervals independently instead, as some descriptions do,
#' does not keep the dependence among the pairs of one distance matrix and
#' gives no test of stated size. Memory grows with the number of events,
#' not pairs; time with pairs times permutations.
#'
#' @inheritParams cgstat
#' @param t The time of each event: finite numbers (for dates, `as.numeric()`
#'   of them, in days or seconds).
#' @param cs,ct The constants of the reciprocal transforms: numbers from 0
#'   (0 uses the raw distance or interval).
#' @param nsim The number of permutations: one whole number from 0.
#' @inheritParams nnidx
#' @param metric The spatial distance: `"euclid"` or `"manhattan"`; the
#'   test permutes times and needs no expectation, so both are exact.
#' @return A named numeric vector: `r` (the correlation), `pairs`, and the
#'   permutation p-values `pup` (`(1 + #{r_perm >= r}) / (nsim + 1)`, large
#'   r: events close in space are close in time when both are raw
#'   distances), `plo` (`r_perm <= r`) and `p2` (`|r_perm| >= |r|`); `NaN`
#'   without permutations. 3 to 46340 events.
#' @references Mantel, N. (1967). The detection of disease clustering and a
#'   generalized regression approach. *Cancer Research*, 27(2), 209-220.
#' @examples
#' set.seed(1)
#' x <- runif(60)
#' y <- runif(60)
#' mantst(x, y, t = x + rnorm(60, sd = 0.1), nsim = 199, seed = 1)
#' @export
mantst <- function(x, y = NULL, t, cs = 0, ct = 0, nsim = 999, seed = NULL,
                   metric = "euclid", lonlat = FALSE) {
  .Call(C_clab_r_mantst, x, y, t, cs, ct, nsim, seed, metric, lonlat)
}

#' Spatial-temporal moving average
#'
#' For each event, the mean centre of the `span` events centred on it in
#' time order (ties in input order), the window shortened at the ends of
#' the series. It traces how the centre of a series of incidents moves
#' over time. It is a descriptive device of the crime-mapping literature
#' with no original paper to cite; the mean centre is that of [cgstat()].
#'
#' @inheritParams mantst
#' @param span The window: one odd whole number from 1.
#' @return A data frame with one row per event, in input order: `mx`, `my`
#'   (the window's mean centre, in degrees when `lonlat = TRUE`) and `n`
#'   (events in the window).
#' @examples
#' stmavg(x = 1:7, y = c(1, 1, 2, 3, 5, 8, 13), t = 7:1, span = 3)
#' @export
stmavg <- function(x, y = NULL, t, span = 5, lonlat = FALSE) {
  .Call(C_clab_r_stmavg, x, y, t, span, lonlat)
}

#' Distance-decay functions for journey to crime
#'
#' The decay of the chance of offending with distance from an offender's
#' base, in the forms of Brantingham and Brantingham (1981), Rossmo (2000)
#' and Levine (2009): `"linear"`, `f = max(0, A + B d)`, `par = c(A, B)`;
#' `"negexp"`, `f = A exp(-b d)`, `par = c(A, b)` with `b >= 0`;
#' `"normal"`, `f = A dnorm(d, mu, s)`, `par = c(A, mu, s)`; `"lognormal"`,
#' `f = A dlnorm(d, mu, s)`, `par = c(A, mu, s)` (the lognormal density of
#' `d`; some texts use `log d^2` and `1 / d^2`, which is not a density);
#' `"trunc"`, the truncated negative exponential, `f = P d / dp` up to the
#' peak `dp` and `P exp(-C (d - dp))` beyond, `par = c(P, dp, C)`; and
#' `"table"`, values `tab` at distances `0, delta, 2 delta, ...`, linear in
#' between and 0 beyond, `par = delta` (as [jtcal()] tabulates).
#'
#' @param d Distances: numbers from 0.
#' @param type The form: one of the names above.
#' @param par The parameters of the form, of the length given above.
#' @param tab For `type = "table"`: at least two values from 0.
#' @return A numeric vector of `f(d)`, one per distance.
#' @references Brantingham, P. J. and Brantingham, P. L. (1981).
#'   *Environmental Criminology*. Sage.
#'
#'   Levine, N. (2009). Introduction to the special issue on Bayesian
#'   journey-to-crime modelling. *Journal of Investigative Psychology and
#'   Offender Profiling*, 6(3), 167-185. \doi{10.1002/jip.107}
#' @examples
#' jtdecay(c(0, 0.5, 1, 2, 5), type = "trunc", par = c(1, 0.5, 0.8))
#' @export
jtdecay <- function(d, type, par, tab = NULL) {
  .Call(C_clab_r_jtdecay, d, type, par, tab)
}

#' Journey-to-crime surface
#'
#' Ranks the cells of a grid as candidate bases of an offender from the
#' locations of a series of incidents. `mode = "sum"` gives the classic
#' aggregate surface, the sum over incidents of the decay of their distance
#' to the cell (Canter et al. 2000; Levine 2009); `mode = "lik"` gives the
#' sum of the logarithms, the log-likelihood of the cell as base when the
#' decay function is the planar density of an incident about the base and
#' the incidents are independent (O'Leary 2009), so that `prob` is the
#' posterior over cells under a uniform prior. The sum surface is not a
#' likelihood.
#'
#' What this does not license: the surface ranks places to search under
#' the assumptions that the series is one offender's, that the offender
#' travelled from one fixed base inside the grid, and that the decay
#' function, which describes offenders in general, fits this one. It is a
#' search priority, not evidence that any person lives anywhere. Accuracy on
#' solved series in the same place ([jtacc()]) is the only guide to its
#' worth.
#'
#' @inheritParams nnidx
#' @param win The grid's rectangle `c(xmin, xmax, ymin, ymax)`.
#' @param nx,ny The numbers of cells across and up: whole numbers from 1,
#'   with at most 1e9 cells.
#' @inheritParams jtdecay
#' @param metric `"euclid"` (straight-line distance) or `"manhattan"`
#'   (`|dx| + |dy|`, a street grid aligned with the axes).
#' @param mode `"sum"` or `"lik"`.
#' @return A list. `raw` and `prob`, matrices of `nx` rows (across) and
#'   `ny` columns (up), cell `[i, j]` centred at `(xmin + (i - 1/2) dx,
#'   ymin + (j - 1/2) dy)`: the surface and its normalised version (summing
#'   to 1). `peak`, a named vector: `cell` (the first cell of maximum value,
#'   numbered down the columns of the matrix), `x` and `y` (its centre).
#' @references Canter, D., Coffey, T., Huntley, M. and Missen, C. (2000).
#'   Predicting serial killers' home base using a decision support system.
#'   *Journal of Quantitative Criminology*, 16(4), 457-478.
#'   \doi{10.1023/A:1007551316253}
#'
#'   O'Leary, M. (2009). The mathematics of geographic profiling. *Journal
#'   of Investigative Psychology and Offender Profiling*, 6(3), 253-265.
#'   \doi{10.1002/jip.111}
#' @examples
#' s <- jtsurf(c(2, 3, 4), c(2, 4, 3), win = c(0, 6, 0, 6), nx = 30,
#'             ny = 30, type = "negexp", par = c(1, 1), mode = "lik")
#' s$peak
#' @export
jtsurf <- function(x, y = NULL, win, nx, ny, type, par, tab = NULL,
                   metric = "euclid", mode = "sum") {
  .Call(C_clab_r_jtsurf, x, y, win, nx, ny, type, par, tab, metric, mode)
}

#' Empirical distance decay from solved trips
#'
#' A kernel density estimate of trip lengths from known offender trips,
#' to calibrate a decay function (Levine 2009; Silverman 1986). With
#' `refl = TRUE` the kernels are reflected at 0, so no mass falls below 0
#' and the estimate integrates to 1 over `[0, Inf)`. Tabulated at evenly
#' spaced distances it is a `"table"` decay for [jtdecay()] and [jtsurf()].
#'
#' @param d Trip lengths: numbers from 0.
#' @param t The distances to evaluate at (from 0 when `refl = TRUE`).
#' @param h The bandwidth: one positive number.
#' @param w `NULL` for equal weights, or one weight from 0 per trip, with a
#'   positive sum.
#' @param kern The kernel: `"normal"`, `"quartic"`, `"triangular"`,
#'   `"uniform"` or `"laplace"` (double exponential).
#' @param refl `TRUE` to reflect at 0.
#' @return A numeric vector of densities, one per value of `t`.
#' @references Silverman, B. W. (1986). *Density Estimation for Statistics
#'   and Data Analysis*. Chapman and Hall. \doi{10.1201/9781315140919}
#'
#'   Levine, N. (2009). Introduction to the special issue on Bayesian
#'   journey-to-crime modelling. *Journal of Investigative Psychology and
#'   Offender Profiling*, 6(3), 167-185. \doi{10.1002/jip.107}
#' @examples
#' set.seed(1)
#' trips <- rexp(200, 1 / 2)
#' jtcal(trips, t = seq(0, 10, 0.5), h = 0.5)
#' @export
jtcal <- function(d, t, h, w = NULL, kern = "normal", refl = TRUE) {
  .Call(C_clab_r_jtcal, d, t, h, w, kern, refl)
}

#' Accuracy of a journey-to-crime surface on a solved series
#'
#' Given a surface over the grid of [jtsurf()] and the offender's known
#' base: the surface's value in the base's cell, the hit score (the share
#' of the grid searched, in decreasing order of the surface, up to and
#' including the base's cell), in the best and the worst order of ties
#' (Rossmo 2000; Canter et al. 2000), and the error distance from the peak
#' cell's centre to the base.
#'
#' @param prob The surface: a matrix of values from 0, `nx` rows by `ny`
#'   columns, as [jtsurf()] returns.
#' @param win The grid's rectangle, as given to [jtsurf()].
#' @param base The known base, `c(x, y)`, inside `win`.
#' @return A named numeric vector: `p` (the value in the base's cell),
#'   `hitbest` and `hitworst` (shares from 0 to 1), `err` (the error
#'   distance), `peak` and `cell` (the peak's and the base's cell, numbered
#'   down the columns).
#' @references Rossmo, D. K. (2000). *Geographic Profiling*. CRC Press.
#' @examples
#' s <- jtsurf(c(2, 3, 4), c(2, 4, 3), win = c(0, 6, 0, 6), nx = 30,
#'             ny = 30, type = "negexp", par = c(1, 1), mode = "lik")
#' jtacc(s$prob, win = c(0, 6, 0, 6), base = c(3.2, 2.9))
#' @export
jtacc <- function(prob, win, base) {
  .Call(C_clab_r_jtacc, prob, win, base)
}

#' Bayesian journey to crime
#'
#' Combines a journey-to-crime surface over origin zones with an
#' origin-destination matrix of solved offences (Levine 2009; Levine and
#' Block 2011): `P(Jtc)`, the surface; `P(O)`, the origins of all solved
#' offences; `P(O|Jtc)`, the origins of offences in the zones where the
#' series happened; their product, Bayes' numerator with the surface as
#' prior; and the risk, the product over `P(O)`, renormalised. Each block
#' sums to 1. The matrix describes past solved offences and inherits their
#' selection.
#'
#' @param pj The journey-to-crime value of each origin zone: numbers from
#'   0, not all 0.
#' @param od The origin-destination counts: a matrix with one row per origin
#'   zone and one column per destination zone, numbers from 0.
#' @param dest The destination zone (column of `od`) of each incident of
#'   the series: whole numbers from 1; repeats count.
#' @return A data frame with one row per origin zone: `pj`, `po`, `pod`,
#'   `prod` and `risk` (0 where `po` is 0).
#' @references Levine, N. and Block, R. (2011). Bayesian journey-to-crime
#'   estimation: an improvement in geographic profiling methodology. *The
#'   Professional Geographer*, 63(2), 213-229.
#'   \doi{10.1080/00330124.2010.547152}
#' @examples
#' od <- matrix(c(5, 1, 0, 2, 4, 1, 0, 1, 6), 3)
#' bjtc(pj = c(0.5, 0.3, 0.2), od = od, dest = c(1, 2, 2))
#' @export
bjtc <- function(pj, od, dest) {
  .Call(C_clab_r_bjtc, pj, od, dest)
}
