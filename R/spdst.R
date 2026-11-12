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
                  dist = FALSE) {
  .Call(C_clab_r_nnidx, x, y, win, kmax, nsim, seed, dist)
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
#' @return A data frame with one row per radius: `r`, `k` (in the square of
#'   the coordinates' unit) and `l`.
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
ripk <- function(x, y = NULL, win, r, corr = "iso") {
  .Call(C_clab_r_ripk, x, y, win, r, corr)
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
