#' Age-period-cohort model by the intrinsic estimator
#'
#' Asks how much of the change in a rate (arrests, imprisonment,
#' victimisation per head) over the years belongs to age, to the period
#' and to the generation (cohort) people were born into. From a table of
#' counts by age group (rows) and period (columns), with equal widths, and
#' the population at risk in each cell, the Poisson log-linear model
#' \deqn{\log E(y_{ij}) = \log(expo_{ij}) + \mu + \alpha_i + \beta_j +
#' \gamma_{j - i + A}}{log E(y_ij) = log(expo_ij) + mu + alpha_i +
#' beta_j + gamma_(j - i + A)}
#' has effects of age, period and cohort, each summing to 0. Because cohort
#' = period - age, the three linear trends cannot be told apart: the
#' design has one null direction, and the data fit equally well along it.
#' The intrinsic estimator (Yang, Fu and Land 2004) picks the solution with
#' no component along that direction, the one a principal-components
#' (minimum-norm) fit gives; its covariance is the model-based one
#' projected the same way. The model is fitted by the iteratively
#' reweighted least squares of `glm()` (as [cntreg()] fits).
#'
#' The constraint is a choice, like any other that identifies the model:
#' the age, period and cohort trends it reports are one of many that fit
#' the same counts, and differ from another constraint's by a linear trend
#' in each. The curvatures (deviations from the trends) are identified
#' whatever the constraint.
#'
#' @param y The counts: a numeric matrix with one row per age group (youngest
#'   first) and one column per period (earliest first), at least three of
#'   each, numbers from 0, not all 0. Row and column names label the
#'   effects.
#' @param expo The population at risk of each cell: a matrix of the same
#'   dimensions, positive.
#' @param level The confidence level of the Wald intervals: one number
#'   strictly between 0 and 1.
#' @param maxit The steps of the fit: one whole number from 1.
#' @return An object of class `apcie`, a list: `coef`, one row for the
#'   intercept and for every age group, period and cohort (`term`, `level`,
#'   `est`, `se`, `lo`, `hi`; cohorts numbered from 1, the oldest, the
#'   oldest age group in the first period); `fit`, a named vector
#'   (`deviance`, `df`, `loglik`, `iter`, `conv`); and `dims`.
#' @references Yang, Y., Fu, W. J. and Land, K. C. (2004). A methodological
#'   comparison of age-period-cohort models: the intrinsic estimator and
#'   conventional generalized linear models. *Sociological Methodology*,
#'   34(1), 75-110. \doi{10.1111/j.0081-1750.2004.00148.x}
#' @examples
#' # arrests by five-year age group (rows) and five-year period (columns)
#' set.seed(11)
#' ag <- c(-1, 0.5, 0.2, -0.3, -0.6)
#' pe <- c(-0.2, 0, 0.1, 0.1)
#' co <- seq(-0.4, 0.4, length.out = 8)
#' ck <- outer(1:5, 1:4, function(i, j) j - i + 5)
#' eta <- outer(ag, pe, "+") + matrix(co[ck], 5)
#' pop <- matrix(1e5, 5, 4)
#' y <- matrix(rpois(20, pop * exp(-6 + eta)), 5,
#'             dimnames = list(c("15-19", "20-24", "25-29", "30-34", "35-39"),
#'                             c(2000, 2005, 2010, 2015)))
#' apcie(y, pop)
#' @export
apcie <- function(y, expo, level = 0.95, maxit = 100) {
  .Call(C_clab_r_apcie, y, expo, level, maxit)
}

#' @export
print.apcie <- function(x, ...) {
  f <- x[["fit"]]
  d <- x[["dims"]]
  cat(
    "apcie: ",
    d[["ages"]],
    " ages, ",
    d[["periods"]],
    " periods, ",
    d[["cohorts"]],
    " cohorts; deviance ",
    format(signif(f[["deviance"]], 5)),
    " on ",
    f[["df"]],
    " df\n",
    sep = ""
  )
  print(x[["coef"]], row.names = FALSE)
  invisible(x)
}
