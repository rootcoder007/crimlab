#' Bounds on sentence effects under monotone response and selection
#'
#' Asks what can be said about the effect of a harsher sentence (a longer
#' term, custody against probation) on an outcome such as reconviction
#' when who received which sentence was not random. With no assumption the
#' mean outcome everyone would have under sentence `t` is known only to an
#' interval (see [sntbnd()] for two sentences). Manski and Pepper (2000)
#' narrow it with two assumptions about ordered treatments:
#'
#' * monotone treatment response (MTR): each person's outcome does not fall
#'   as the sentence gets harsher (code `y` so that this holds); then
#'   those who got a milder sentence `u < t` would have at least their
#'   observed outcome under `t`, and those who got a harsher one at most
#'   theirs;
#' * monotone treatment selection (MTS): those who received harsher
#'   sentences would, on average, have had outcomes at least as high under
#'   any one sentence (the courts sentence higher-risk people more
#'   harshly).
#'
#' With `P(u)` the share and `m(u)` the mean outcome of those who received
#' sentence `u`, and outcomes between `ylo` and `yhi`, the bound on the
#' mean under `t` is, under MTR, from `sum_{u <= t} m(u) P(u) + ylo P(z >
#' t)` to `yhi P(z < t) + sum_{u >= t} m(u) P(u)`; under MTS, from `ylo
#' P(z < t) + m(t) P(z >= t)` to `m(t) P(z <= t) + yhi P(z > t)`; and under
#' both, from `sum_{u < t} m(u) P(u) + m(t) P(z >= t)` to `m(t) P(z <= t) +
#' sum_{u > t} m(u) P(u)` (their Proposition 2). The effect of `t` against
#' the first (mildest) sentence is bounded by the difference of the
#' bounds, its lower end at least 0 under MTR; these differences contain
#' the effect but need not be its sharpest bounds.
#'
#' @param y The outcome of each person: numbers between `ylo` and `yhi`
#'   (0 or 1 for reconviction).
#' @param z The sentence of each person: labels in their order of severity
#'   (a factor's levels, or the sorted values), two at least.
#' @param ylo,yhi The least and the greatest possible outcome: one number
#'   each, `ylo` below `yhi`.
#' @param assume `"none"`, `"mtr"`, `"mts"` or `"both"`.
#' @param w `NULL` for equal weights, or one weight from 0 per person.
#' @return An object of class `mtrbnd`, a list: `bnd`, one row per sentence
#'   (`z`, `share`, `mean`, `lower` and `upper` of the mean outcome under
#'   it, and `atelo`, `atehi`, the bounds on its effect against the first
#'   sentence, `NA` for that one); `assume`; and `range`.
#' @references Manski, C. F. and Pepper, J. V. (2000). Monotone
#'   instrumental variables: with an application to the returns to
#'   schooling. *Econometrica*, 68(4), 997-1010.
#'   \doi{10.1111/1468-0262.00144}
#' @examples
#' set.seed(13)
#' z <- factor(sample(c("fine", "probation", "custody"), 800, TRUE),
#'             levels = c("fine", "probation", "custody"))
#' y <- rbinom(800, 1, c(0.2, 0.35, 0.5)[as.integer(z)])
#' mtrbnd(y, z, assume = "both")
#' mtrbnd(y, z, assume = "none")
#' @export
mtrbnd <- function(y, z, ylo = 0, yhi = 1, assume = "both", w = NULL) {
  .Call(C_clab_r_mtrbnd, y, z, ylo, yhi, assume, w)
}

#' @export
print.mtrbnd <- function(x, ...) {
  cat("mtrbnd: bounds under ", x[["assume"]], "\n", sep = "")
  print(x[["bnd"]], row.names = FALSE)
  invisible(x)
}
