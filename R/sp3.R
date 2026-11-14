#' Head-bang smoothing
#'
#' Smooths values at points or areal units (rates by area centroid, say)
#' by Hansen's (1991) head-bang rule, which removes isolated peaks and
#' troughs while keeping edges between regions of different level. For
#' each unit `i` the rule takes its `nn` nearest neighbours and every pair
#' `(j, k)` of them whose angle at `i` is at least `ang` degrees, so that
#' `i` lies roughly between `j` and `k` (a triple); with `ntrip > 0` only the
#' `ntrip` triples whose ends are nearest (smallest `|p_j - p_i| +
#' |p_k - p_i|`) are kept. Each triple gives a low value `min(v_j, v_k)` and
#' a high value `max(v_j, v_k)`; the low screen is the median of the lows,
#' the high screen the median of the highs, and the new value is the median
#' of the low screen, `v_i` and the high screen, so `v_i` moves only if it
#' lies outside both. With weights `w` (populations at risk, say) the
#' screens are weighted medians, each end weighted by its own unit's
#' weight, after Mungiole, Pickle and Simonson (1999): the smallest value
#' whose cumulative weight reaches half the total. A unit with no triple
#' (on the outer edge of the points) keeps its value; Hansen's extrapolated
#' end triples are not formed. Each of `iter` passes starts from the
#' previous pass's values.
#'
#' @inheritParams cgstat
#' @param v The values to smooth: finite numbers, one per point.
#' @param w `NULL`, or positive weights, one per point.
#' @param nn The number of nearest neighbours searched: from 2 to the
#'   number of points less 1.
#' @param ang The smallest angle at `i` of a triple, in degrees, from 90 to
#'   180 (135 by convention).
#' @param ntrip The most triples used per unit, nearest first; 0 for all.
#' @param iter The number of passes: a whole number from 0.
#' @return A numeric vector of smoothed values, one per point, in input
#'   order.
#' @references Hansen, K. M. (1991). Head-banging: robust smoothing in the
#'   plane. *IEEE Transactions on Geoscience and Remote Sensing*, 29(3),
#'   369-378.
#'   \doi{10.1109/36.79427}
#'
#'   Mungiole, M., Pickle, L. W. and Simonson, K. H. (1999). Application of
#'   a weighted head-banging algorithm to mortality data maps. *Statistics
#'   in Medicine*, 18(23), 3201-3209.
#'   \doi{10.1002/(SICI)1097-0258(19991215)18:23<3201::AID-SIM310>3.0.CO;2-U}
#' @examples
#' g <- expand.grid(x = 1:6, y = 1:6)
#' v <- rep(1, 36)
#' v[15] <- 50
#' hdbang(g, v = v)[15]
#' @export
hdbang <- function(x, y = NULL, v, w = NULL, nn = 8, ang = 135, ntrip = 0,
                   iter = 1, lonlat = FALSE) {
  .Call(C_clab_r_hdbang, x, y, v, w, nn, ang, ntrip, iter, lonlat)
}

#' Correlated random walk of a series of events
#'
#' Treats the successive locations of a series (incidents attributed to
#' one offender, in time order) as a correlated random walk (Kareiva and
#' Shigesada 1983): steps of independent lengths and changes of heading
#' (turning angles) independent of each other and of the steps. The fit
#' gives the mean and mean squared step and the mean cosine and sine of
#' the turns; the expected squared displacement after `k` moves is
#' Kareiva and Shigesada's for symmetric turning (mean sine 0), `k E l^2 +
#' 2 (E l)^2 c / (1 - c) (k - (1 - c^k) / (1 - c))`, next to the observed
#' mean over every start. The walk is then continued `nsim` times for `h`
#' steps from the last location and heading, each turn and step drawn
#' independently and uniformly from the observed ones (an empirical
#' correlated random walk), and the simulated positions summarise where
#' the next events may fall. It describes movement; it does not identify
#' anyone, and a series with few events gives a rough walk.
#'
#' @inheritParams cgstat
#' @param h The number of steps ahead: a whole number from 1.
#' @param nsim The number of simulated continuations: a whole number from 1.
#' @inheritParams nnidx
#' @return A list. `fit`, a named vector: `steps`, `ml` (mean step), `ml2`
#'   (mean squared step), `c` and `s` (mean cosine and sine of the turning
#'   angles) and `head` (the last heading, radians anticlockwise from
#'   east). `msd`, a data frame with one row per step ahead `k`: `obs` (the
#'   observed mean squared displacement over `k` moves, `NA` beyond the
#'   series), `exp` (its expectation; it assumes `s` near 0), and `x`, `y`,
#'   `sdx`, `sdy` (the mean and standard deviations of the simulated
#'   positions; the means in degrees when `lonlat = TRUE`, the rest in
#'   metres). `sim`, the simulated positions: matrices `x` and `y`, one row
#'   per walk and one column per step (planar, metres when `lonlat = TRUE`).
#' @references Kareiva, P. M. and Shigesada, N. (1983). Analyzing insect
#'   movement as a correlated random walk. *Oecologia*, 56(2-3), 234-238.
#'   \doi{10.1007/BF00379695}
#' @examples
#' set.seed(1)
#' a <- cumsum(rnorm(20, 0, 0.5))
#' p <- cbind(x = cumsum(cos(a)), y = cumsum(sin(a)))
#' r <- crwalk(p, h = 3, nsim = 200, seed = 1)
#' r$fit
#' r$msd
#' @export
crwalk <- function(x, y = NULL, h = 1, nsim = 999, seed = NULL,
                   lonlat = FALSE) {
  .Call(C_clab_r_crwalk, x, y, h, nsim, seed, lonlat)
}

#' Risk-adjusted nearest-neighbour hierarchical clustering
#'
#' Nearest-neighbour clustering (as [nnhclu()], first order) with the
#' threshold distance adjusted for a baseline population at risk, so that
#' incidents dense only because people are dense there are not called a
#' cluster. In [nnhclu()] the threshold is the mean nearest-neighbour
#' distance under complete spatial randomness plus `z` standard errors,
#' `0.5 sqrt(A / n) + z sqrt((4 - pi) A / (4 pi)) / n`, which is
#' `f / sqrt(lambda)` with `lambda = n / A` and
#' `f = 0.5 + z sqrt((4 - pi) / (4 pi n))`. Here the uniform intensity is
#' replaced by the intensity the `n` incidents would have if they followed
#' the baseline: `lambda_i = n b(p_i) / B`, where `b` is the normal-kernel
#' density (bandwidth `h`) of the baseline units with sizes `bw` and `B` is
#' their total, and two incidents `i` and `j` are linked when
#' `|p_i - p_j| < f / sqrt((lambda_i + lambda_j) / 2)`. Taking the mean of
#' the two intensities is this function's reading of "the threshold at the
#' pair's baseline density" (the published description does not say how
#' the two ends combine); a uniform baseline gives back [nnhclu()]'s
#' threshold. Grouping is [nnhclu()]'s first order; higher orders are not
#' formed.
#'
#' @inheritParams cgstat
#' @param base The baseline units: coordinates as for `x` (zone centroids,
#'   or points), in the same system as `x` (degrees when `lonlat = TRUE`).
#' @param bw `NULL` for size 1 each, or the baseline sizes (populations),
#'   numbers from 0 with a positive sum.
#' @param h The bandwidth of the baseline's normal kernel, positive, in the
#'   unit of the coordinates (metres when `lonlat = TRUE`).
#' @param minpts The fewest points a cluster keeps: a whole number from 2.
#' @inheritParams nnhclu
#' @return A list: `cl`, each point's cluster (a row of `clusters`) or `NA`;
#'   `clusters`, a data frame of `points`, `mx`, `my`, `wx`, `wy`, `a`, `b`
#'   and `angle` as in [nnhclu()]; `thr`, each point's own threshold
#'   `f / sqrt(lambda_i)`.
#' @references Johnson, S. C. (1967). Hierarchical clustering schemes.
#'   *Psychometrika*, 32(3), 241-254. \doi{10.1007/BF02289588}
#'
#'   Clark, P. J. and Evans, F. C. (1954). Distance to nearest neighbor as
#'   a measure of spatial relationships in populations. *Ecology*, 35(4),
#'   445-453. \doi{10.2307/1931034}
#' @examples
#' set.seed(1)
#' town <- cbind(x = rnorm(400, 0, 1), y = rnorm(400, 0, 1))
#' inc <- rbind(town[1:60, ], cbind(x = rnorm(15, 3, 0.1),
#'                                  y = rnorm(15, 3, 0.1)))
#' r <- rnnhcl(inc, base = town, h = 0.5, minpts = 5)
#' r$clusters
#' @export
rnnhcl <- function(x, y = NULL, base, bw = NULL, h, p = 0.5, minpts = 10,
                   k = 1, c = 1, lonlat = FALSE) {
  .Call(C_clab_r_rnnhcl, x, y, base, bw, h, p, minpts, k, c, lonlat)
}
