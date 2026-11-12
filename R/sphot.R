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
#' @return An integer vector, one count per point, in input order.
#' @references Eck, J. E., Chainey, S., Cameron, J. G., Leitner, M. and
#'   Wilson, R. E. (2005). *Mapping Crime: Understanding Hot Spots*.
#'   National Institute of Justice, NCJ 209393.
#' @examples
#' fzmode(c(0, 0, 1, 5), c(0, 0, 0, 5), r = 1)
#' @export
fzmode <- function(x, y = NULL, r, lonlat = FALSE) {
  .Call(C_clab_r_fzmode, x, y, r, lonlat)
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
                   maxit = 100, seed = NULL, lonlat = FALSE) {
  .Call(C_clab_r_kmpart, x, y, k, init, cell, upd, maxit, seed, lonlat)
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
                   lonlat = FALSE) {
  .Call(C_clab_r_nnhclu, x, y, area, fixd, p, minpts, minhi, maxord, k, c,
        lonlat)
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
#' @return A list. `cl`, the hot-spot area of each point (a row of
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
                   minpts = 5, k = 1, c = 1) {
  .Call(C_clab_r_stacel, x, y, win, nx, ny, r, top, minpts, k, c)
}
