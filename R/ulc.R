#' The U-learner of conditional average treatment effects
#'
#' Estimates how the effect of a 0/1 treatment varies with the covariates
#' (Kunzel et al. 2019, sec. 2 and appendix): the outcome regression
#' `m(x) = E[y | x]` and the propensity `e(x) = P(d = 1 | x)` are learned
#' out of fold, as in [dmlx()], each fold predicted by learners fitted on
#' the others; the pseudo-outcome `u = (y - m) / (d - e)` has expectation
#' `tau(x)` given `x` when the treatment is unconfounded given `x`; and
#' `tau` is the learner `lrnt` fitted to `u` on all the rows. The mean of
#' `tau` over the rows is reported as `ate`.
#'
#' The denominator `d - e` is near 0 for a treated row with a propensity
#' near 1 or an untreated row with one near 0, where `u` explodes: it is
#' moved to `clip` (or `-clip`) when nearer 0 than that. The learners, the
#' folds, the clipping and the final fit are all compiled.
#'
#' @inheritParams dmlx
#' @param y The outcome: one finite number per row of `x`.
#' @param d The treatment: 0 or 1, both present.
#' @param fold The fold of each row, whole numbers from 1 to the number of
#'   folds (at least 2), each used.
#' @param lrn The learner of the outcome regression `m`: `"glm"`,
#'   `"lasso"`, `"ridge"`, `"enet"`, `"forest"` or `"stack"`, as in
#'   [dmlx()].
#' @param lrnd The learner of the propensity `e`, from the same choices (a
#'   probability: logistic for `"glm"`).
#' @param lrnt The learner that fits `tau` to the pseudo-outcome, from the
#'   same choices; `"glm"` gives an effect linear in the covariates.
#' @param clip The smallest size of `d - e`: one number above 0 and below
#'   0.5.
#' @param pred `NULL` to learn `m` and `e`, or their predictions: a matrix
#'   with one row per observation and two columns, `m` then `e`.
#' @param seed The seed of the inner folds and of the forests: one number,
#'   at least 0.
#' @return An object of class `ulcate`, a list: `tau`, the effect at each
#'   row; `ate`, its mean; `u`, the pseudo-outcome; `pred`, the
#'   out-of-fold predictions (columns `m` and `e`).
#' @references
#'   Kunzel, S. R., Sekhon, J. S., Bickel, P. J. and Yu, B. (2019).
#'   Metalearners for estimating heterogeneous treatment effects using
#'   machine learning. *Proceedings of the National Academy of Sciences*,
#'   116(10), 4156-4165. \doi{10.1073/pnas.1804597116}
#'
#'   Nie, X. and Wager, S. (2021). Quasi-oracle estimation of heterogeneous
#'   treatment effects. *Biometrika*, 108(2), 299-319.
#'   \doi{10.1093/biomet/asaa076}
#' @examples
#' set.seed(6)
#' n <- 400
#' x <- matrix(rnorm(n * 2), n)
#' d <- rbinom(n, 1, plogis(0.5 * x[, 1]))
#' y <- x[, 2] + d * (1 + x[, 1]) + rnorm(n)
#' fold <- sample(rep(1:4, length.out = n))
#' r <- ulcate(y, d, x, fold, lrn = "glm", lrnd = "glm")
#' r
#' @export
ulcate <- function(y, d, x, fold, lrn = "lasso", lrnd = lrn, lrnt = "glm",
                   clip = 0.01, alpha = 0.5, nlam = 100, nfin = 5, tol = 1e-7,
                   scr = TRUE, hon = FALSE, ntree = 500, mtry = NULL,
                   pred = NULL, stk = c("glm", "lasso", "forest"), nstk = 5,
                   disc = FALSE, seed = 1) {
  .Call(C_clab_r_ulcate, y, d, x, fold, lrn, lrnd, lrnt, clip, alpha, nlam,
        nfin, tol, scr, hon, ntree, mtry, pred, stk, nstk, disc, seed)
}

#' @rdname ulcate
#' @param ... Not used by `print()`.
#' @export
print.ulcate <- function(x, ...) {
  cat("ulcate: U-learner of the effect of d, ", length(x[["tau"]]),
      " rows; mean effect ", format(x[["ate"]], digits = 4), "\n",
      "tau from ", format(min(x[["tau"]]), digits = 4), " to ",
      format(max(x[["tau"]]), digits = 4), "\n",
      sep = "")
  invisible(x)
}
