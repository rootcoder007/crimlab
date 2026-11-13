#' Bivariate Moran's I
#'
#' Whether high values of one variable sit next to high values of another
#' (Wartenberg 1985; Anselin, Syabri and Smirnov 2002): with `zx` and `zy`
#' the centred values, `I = n / S0 sum_ij w_ij zx_i zy_j / sqrt(sum zx^2
#' sum zy^2)`, which is Moran's I ([gmoran()]) when `y` is `x` (proved in
#' Lean: see `proof("bvmor")`) and, with row-standardised weights,
#' `spdep::moran_bv()`'s statistic. Inference is by permutation: `y`'s
#' values are shuffled over the units `nsim` times from `seed` while `x`
#' stays, so under the null its expectation is 0. The statistic mixes the
#' two variables' correlation at a unit with their spatial association; it
#' is not a test of spillover from `x` to `y`.
#'
#' @param x,y The two variables: finite numbers, one per unit, neither
#'   constant.
#' @inheritParams gmoran
#' @param nsim The number of permutations: a whole number from 0.
#' @return A named numeric vector: `i`, `mean` and `sd` (of the permuted
#'   statistics), `z` (`(i - mean) / sd`), `p` (`(1 + #{I* >= I}) / (nsim +
#'   1)`) and `ptwo` (`(1 + min(#{I* >= I}, #{I* <= I})) / (nsim + 1)`);
#'   `NaN` without permutations.
#' @references Wartenberg, D. (1985). Multivariate spatial correlation: a
#'   method for exploratory geographical analysis. *Geographical Analysis*,
#'   17(4), 263-283. \doi{10.1111/j.1538-4632.1985.tb00849.x}
#'
#'   Anselin, L., Syabri, I. and Smirnov, O. (2002). Visualizing
#'   multivariate spatial correlation with dynamically linked windows. In
#'   L. Anselin and S. Rey (eds), *New Tools for Spatial Data Analysis*.
#'   Center for Spatially Integrated Social Science, Santa Barbara.
#' @examples
#' set.seed(1)
#' p <- cbind(x = runif(60), y = runif(60))
#' a <- p[, 1] + rnorm(60, sd = 0.3)
#' b <- p[, 1] + rnorm(60, sd = 0.3)
#' bvmor(a, b, spwt(p, k = 5), nsim = 199, seed = 1)
#' @export
bvmor <- function(x, y, w, nsim = 999, seed = NULL) {
  .Call(C_clab_r_bvmor, x, y, w, nsim, seed)
}
