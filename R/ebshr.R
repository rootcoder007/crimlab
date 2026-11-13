#' Empirical-Bayes shrinkage of hot-spot counts
#'
#' The places with the highest counts in one period tend to fall in the
#' next by chance alone: a count `y = theta + e` is a true rate plus noise
#' of variance `v` (for Poisson counts, the mean count, the default). The
#' linear shrinkage estimate `(1 - B) y + B ybar` toward the weighted mean
#' has expected loss `(1 - B)^2 S_e + B^2 S_theta` when the noise has
#' weighted mean 0 and is uncorrelated with the rates, which is smallest at
#' `B* = S_e / (S_e + S_theta)`, where it never exceeds the loss of the raw
#' counts (James and Stein 1961; Efron and Morris 1975). `S_theta` is
#' estimated as the excess of the weighted scatter of `y` over `v` (0 if
#' negative). The predicted fall of a place is `B (y - ybar)`: the size of
#' the regression to the mean a hot-spot programme should expect without
#' doing anything.
#'
#' @param y The first-period counts (or rates): finite numbers, at least
#'   two.
#' @param v `NULL` for the weighted mean of `y` (Poisson noise), or the
#'   noise variance, a number from 0.
#' @param w `NULL`, or one weight from 0 per place, positive sum.
#' @param b `NULL` for `B*`, or a shrinkage factor from 0 to 1.
#' @return A list of class `ebshr`: `fit`, a named vector: `mean`, `total`
#'   (the weighted scatter of `y`), `noise` (`v`), `signal` (the estimated
#'   variance of the rates), `bstar` and `b` (the factor used); `est`, a
#'   data frame with one row per place: `y`, `shrunk` and `fall`.
#' @references Efron, B. and Morris, C. (1975). Data analysis using Stein's
#'   estimator and its generalizations. *Journal of the American
#'   Statistical Association*, 70(350), 311-319.
#'   \doi{10.1080/01621459.1975.10479864}
#' @examples
#' y <- c(40, 12, 9, 25, 7, 31, 5, 18)
#' s <- ebshr(y)
#' s
#' s$est[which.max(y), ]
#' @export
ebshr <- function(y, v = NULL, w = NULL, b = NULL) {
  .Call(C_clab_r_ebshr, y, v, w, b)
}

#' @export
print.ebshr <- function(x, ...) {
  f <- x[["fit"]]
  cat("ebshr: shrinkage ", format(signif(f[["b"]], 4)), " (B* ",
      format(signif(f[["bstar"]], 4)), ") toward ",
      format(signif(f[["mean"]], 4)), "; signal variance ",
      format(signif(f[["signal"]], 4)), ", noise ",
      format(signif(f[["noise"]], 4)), "\n", sep = "")
  print(x[["est"]], row.names = FALSE)
  invisible(x)
}
