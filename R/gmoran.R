#' Global Moran's I with its inference
#'
#' Moran's (1950) I of values on units with spatial weights `w`: with `z =
#' x - mean(x)`, `I = n / S0 sum_ij w_ij z_i z_j / sum z_i^2`, its
#' expectation `-1 / (n - 1)` and its variance under normality and under
#' randomisation (Cliff and Ord 1981), the z statistics and their upper
#' normal tails (positive autocorrelation: neighbours alike), and, with
#' `nsim > 0`, the permutation p-value `(1 + #{I* >= I}) / (nsim + 1)`, the
#' values shuffled over the units from `seed`. These are the statistics of
#' `spdep::moran.test()` and `spdep::moran.mc()`. Given a matrix or data
#' frame, each column is tested on the same weights: one row per category
#' or year, a sweep of a variable across a polygon surface. [lmorn()] gives
#' the local statistics, which add up to `S0 I` with row-standardised
#' weights.
#'
#' @param x The values: finite numbers, one per unit; or a matrix or data
#'   frame with one column per variable.
#' @param w The weights: [spwt()]'s list, an spdep `listw` object, or an `n
#'   x n` matrix of weights from 0 (zero diagonal).
#' @inheritParams gstr
#' @inheritParams nnidx
#' @return A data frame with one row per variable: `var` (its name), `i`,
#'   `e`, `vn` and `vr`, `zn` and `zr`, `pn` and `pr`, `pperm` (`NaN`
#'   without permutations) and `s0` (the sum of the weights); `NA` but `s0`
#'   for a constant variable.
#' @references Moran, P. A. P. (1950). Notes on continuous stochastic
#'   phenomena. *Biometrika*, 37(1-2), 17-23. \doi{10.2307/2332142}
#'
#'   Cliff, A. D. and Ord, J. K. (1981). *Spatial Processes: Models and
#'   Applications*. Pion, London.
#' @examples
#' set.seed(1)
#' p <- cbind(x = runif(50), y = runif(50))
#' w <- spwt(p, k = 4)
#' gmoran(p[, 1] + rnorm(50, sd = 0.2), w, nsim = 99, seed = 1)
#' # a sweep over three yearly counts
#' gmoran(data.frame(y1 = rpois(50, 3), y2 = rpois(50, 3 + 4 * p[, 1]),
#'                   y3 = rpois(50, 9 * p[, 1])), w)
#' @export
gmoran <- function(x, w, nsim = 0, seed = NULL) {
  if (is.data.frame(x)) x <- as.matrix(x)
  .Call(C_clab_r_gmoran, x, w, nsim, seed)
}
