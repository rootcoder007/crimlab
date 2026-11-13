#' Elastic net by coordinate descent, with cross-validation
#'
#' Fits the lasso, ridge and elastic net path by cyclical coordinate
#' descent, for a numeric (gaussian) or a 0/1 (binomial) outcome, and
#' cross-validates it over folds the caller supplies. The objective, the
#' standardisation of the columns and the scale of `lambda` are those of
#' 'glmnet', so that with the same `lambda` and the same folds the two give
#' the same coefficients, up to their convergence thresholds: the gaussian
#' family minimises `(1 / 2n) |y - a0 - x b|^2 + lambda pen(b)`, the
#' binomial family the mean negative log-likelihood plus `lambda pen(b)`,
#' with `pen(b) = sum((1 - alpha) / 2 (s b)^2 + alpha |s b|)`, `s` the
#' population standard deviations of the columns. For the gaussian family
#' the outcome is scaled to unit variance while fitting, as 'glmnet' does.
#' The intercept is never penalised; a constant column keeps a zero
#' coefficient.
#'
#' Speed against accuracy: `tol` and `maxit` set how far the descent
#' goes. The default `tol = 1e-7` is glmnet's; `tol = 1e-20, maxit = 1e7`
#' is the exact setting, which reproduces glmnet run at the same threshold
#' to 1e-8 in every coefficient (tests). Between the two, the penalised
#' objective at each penalty differs by less than 1e-6 of its value
#' (tested on the gaussian and binomial paths). `nlam` and `lmr` set the
#' length and the end of a path chosen here, and `scr` the screening of a
#' binomial path.
#'
#' Each coordinate step is the minimiser of its one-dimensional objective,
#' by soft-thresholding, the proximal map of the absolute value: see
#' `proof("enpath")`.
#'
#' The binomial fit is a proximal Newton method: iteratively reweighted
#' least squares around coordinate descent over the active columns, with
#' the optimality conditions checked over every column before a `lambda`
#' is left. The gaussian fit updates the
#' gradient through the column inner products, formed only for the
#' columns that enter the model.
#'
#' @param x The predictors: a numeric matrix, one row per observation,
#'   with no missing values.
#' @param y The outcome: one finite number per row of `x`; 0 or 1 for the
#'   binomial family.
#' @param fam `"gaussian"` or `"binomial"`.
#' @param alpha The mixing parameter, from 0 (ridge) to 1 (lasso).
#' @param lambda `NULL`, or the penalties to fit: non-negative numbers in
#'   decreasing order.
#' @param nlam When `lambda` is `NULL`, the number of penalties on the
#'   path: one whole number, at least 1.
#' @param lmr When `lambda` is `NULL`, the smallest penalty as a share of
#'   the largest: one number strictly between 0 and 1, or `NULL` for 1e-4
#'   when `x` has more rows than columns and 0.01 otherwise. The path runs
#'   from the smallest penalty that zeroes every coefficient down to `lmr`
#'   times it, evenly spaced on the log scale, and stops early as 'glmnet'
#'   does: from the fifth penalty on, once the share of the deviance
#'   explained grows by less than 1e-5 of itself or exceeds 0.999.
#' @param fold `NULL`, or the fold of each row: whole numbers from 1 to the
#'   number of folds (at least 2), each fold used.
#' @param tol The convergence threshold: a pass that moves no standardised
#'   coefficient by more than `sqrt(tol)` (weighted by the curvature, for
#'   the binomial family) ends the descent. One positive number.
#' @param maxit The most passes over the columns, in all: one whole number.
#' @param scr `TRUE` to screen the columns of a binomial path by the
#'   sequential strong rule of Tibshirani et al. (2012) before each penalty
#'   (faster; the optimality conditions are still checked over every
#'   column, so the solution does not change), `FALSE` to cycle over the
#'   active columns only. Ignored for the gaussian family.
#' @return A list: `lambda`, the penalties fitted; `a0`, the intercepts;
#'   `beta`, a matrix of coefficients with one column per penalty and one
#'   row per column of `x`; and, when `fold` is given, `cvm` and `cvsd`, the
#'   mean held-out loss (squared error, or the binomial deviance with
#'   probabilities kept within 1e-5 of 0 and 1) and its standard error over
#'   folds, computed as 'glmnet' computes them, `lmin`, the penalty at the
#'   minimum, and `l1se`, the largest penalty within one standard error of
#'   it; these four are `NULL` without folds.
#' @references Friedman, J., Hastie, T. and Tibshirani, R. (2010).
#'   Regularization paths for generalized linear models via coordinate
#'   descent. *Journal of Statistical Software*, 33(1), 1-22.
#'   \doi{10.18637/jss.v033.i01}
#'
#'   Tibshirani, R., Bien, J., Friedman, J., Hastie, T., Simon, N., Taylor,
#'   J. and Tibshirani, R. J. (2012). Strong rules for discarding predictors
#'   in lasso-type problems. *Journal of the Royal Statistical Society:
#'   Series B*, 74(2), 245-266. \doi{10.1111/j.1467-9868.2011.01004.x}
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(200 * 5), 200)
#' y <- x[, 1] - 0.5 * x[, 2] + rnorm(200)
#' f <- enpath(x, y, nlam = 20, fold = rep(1:5, 40))
#' f$lmin
#' f$beta[, f$lambda == f$lmin]
#' @export
enpath <- function(x, y, fam = "gaussian", alpha = 1, lambda = NULL,
                  nlam = 100, lmr = NULL, fold = NULL, tol = 1e-7,
                  maxit = 1e5, scr = TRUE) {
  .Call(C_clab_r_enpath, x, y, fam, alpha, lambda, nlam, lmr, fold, tol, maxit,
        scr)
}

#' Logistic regression
#'
#' Maximum likelihood for a 0/1 outcome and a logit link, by Newton's
#' method on the design with centred columns, the method `glm()` uses with
#' `family = binomial()` when the columns are not collinear. Standard
#' errors come from the inverse Fisher information.
#'
#' @param x The predictors: a numeric matrix with more rows than columns
#'   plus one, with no missing values; an intercept is always added.
#' @param y The outcome: 0 or 1, one per row of `x`, both present.
#' @param tol Newton steps stop when none moves a coefficient by more than
#'   `tol (1 + |coefficient|)`: one positive number.
#' @param maxit The most Newton steps: one whole number.
#' @return A list: `coef` and `se`, named numeric vectors with the
#'   intercept first, and `dev`, the residual deviance.
#' @references McCullagh, P. and Nelder, J. A. (1989). *Generalized Linear
#'   Models* (2nd ed.). Chapman and Hall. \doi{10.1007/978-1-4899-3242-6}
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(300 * 2), 300)
#' y <- rbinom(300, 1, plogis(0.5 + x[, 1] - x[, 2]))
#' lgtreg(x, y)$coef
#' @export
lgtreg <- function(x, y, tol = 1e-10, maxit = 50) {
  .Call(C_clab_r_lgtreg, x, y, tol, maxit)
}
