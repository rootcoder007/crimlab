#' The parts of a fitted regression
#'
#' The fits of [ppmlfe()], [ivfe()], [ppmlcf()], [expgmm()], [fepnl()],
#' [gmmpnl()], [didfe()], [rdloc()] and [sarml()] are lists, and the
#' usual accessors read them: `coef()` gives the estimates by term,
#' `confint()` their intervals, `vcov()` the covariance matrix of the
#' estimates and `nobs()` the rows the fit used.
#'
#' The intervals are those the fit computed, at the fit's own level (the
#' bias-corrected interval of [rdloc()], the bootstrap interval of
#' [ppmlcf()]): they are read from the fit, not rebuilt from the
#' standard errors, so `confint()` takes no other level; fit again with
#' the `level` wanted. [sarml()] stores no interval, and its `confint()`
#' is the normal one, estimate plus or minus the quantile times the
#' standard error, at any `level`.
#'
#' [ppmlcf()] (bootstrap) and [rdloc()] (one estimate by two methods)
#' store no covariance matrix, so they have no `vcov()`. The matrix of
#' [sarml()] has a row for `rho` and one for `lambda` whichever the model
#' holds (`NaN` for the one it does not), as its help page says.
#'
#' @section How a model is given:
#' By the names of columns of `data`, or by vectors, each in the argument
#' of its role: the outcome, the regressors, the instruments, the fixed
#' effects, the clusters. There is no formula interface: a formula has no
#' place for those roles, and the fixed effects are absorbed, not
#' expanded into columns of a model matrix. Regressors are numeric
#' columns; a column of labels is taken only as a fixed effect or a
#' cluster. Missing and infinite values are refused, each in its own
#' words, and no option drops them: the caller subsets first.
#'
#' @section What a fit keeps:
#' The terms keep the names of the input columns, in `coef()`,
#' `confint()` and `vcov()`. Row names are not kept: a fit has one row
#' per term, not per observation. Fitted values need the fixed effects,
#' which are absorbed and never estimated, so the fits do not return
#' them; [fepnl()] and [gmmpnl()] return their residuals (`resid`).
#' Collinear regressors, and regressors the fixed effects absorb, stop
#' the fit with a message that says so. A count fit that has not
#' converged warns, and records it (`conv`, `iter`); `tol` and `maxit`
#' set the stopping rule.
#'
#' @section Speed:
#' The time of a fit grows in proportion to the number of rows: ten
#' times the rows take about ten times as long (measured for [ivfe()],
#' [ppmlfe()], [fepnl()] and [didfe()] from a thousand to a hundred
#' thousand rows).
#'
#' @param object A fit of one of the functions above.
#' @param parm `NULL` (every term) or the terms wanted: names of terms
#'   or their positions.
#' @param level The level of the intervals: `NULL` (the fit's own; 0.95
#'   for [sarml()]) or one number between 0 and 1. For every fit but
#'   that of [sarml()] it must be the fit's level.
#' @param ... Not used.
#' @return `coef()`: a named numeric vector, one estimate per term (per
#'   method for [rdloc()]). `confint()`: a matrix with a row per term
#'   and the columns `lower` and `upper`, with the level as the
#'   attribute `"level"`. `vcov()`: a square matrix named by the terms.
#'   `nobs()`: one whole number.
#' @examples
#' set.seed(1)
#' d <- data.frame(u = rep(1:20, each = 6), x = rnorm(120))
#' d$y <- 1 + 0.5 * d$x + rnorm(120)
#' f <- ivfe(d, "y", x = "x", fe = "u")
#' coef(f)
#' confint(f)
#' vcov(f)
#' nobs(f)
#' @name racc
#' @rawNamespace importFrom(stats,coef)
#' @rawNamespace importFrom(stats,confint)
#' @rawNamespace importFrom(stats,nobs)
#' @rawNamespace importFrom(stats,vcov)
NULL

#' The table of estimates of a fit and the names of its rows: the terms
#' (the methods for rdloc(), whose table has no column `term`).
#' @noRd
.racf <- function(x, fn) {
  cf <- x[["coef"]]
  if (!is.data.frame(cf) || is.null(cf[["est"]])) {
    .ckstop(fn, "`object` is not a complete fit (it has no table `coef`)")
  }
  nm <- if (is.null(cf[["term"]])) cf[["method"]] else cf[["term"]]
  list(cf = cf, nm = as.character(nm))
}

#' coef() of a regression fit: the estimates named by term.
#' @noRd
.racoef <- function(object, ...) {
  a <- .racf(object, "coef")
  stats::setNames(as.numeric(a$cf[["est"]]), a$nm)
}

#' The rows `parm` asks for, as positions.
#' @noRd
.raparm <- function(parm, nm) {
  if (is.null(parm)) {
    return(seq_along(nm))
  }
  i <- if (is.character(parm)) match(parm, nm) else parm
  if (!is.numeric(i) || !length(i) || anyNA(i) || any(i < 1) ||
        any(i > length(nm)) || any(i != round(i))) {
    .ckstop("confint", "`parm` must be NULL, or terms of the fit: %s",
            toString(nm))
  }
  as.integer(i)
}

#' confint() of a regression fit: the intervals the fit computed, or the
#'   normal ones when it stored none.
#' @noRd
.raci <- function(object, parm = NULL, level = NULL, ...) {
  a <- .racf(object, "confint")
  i <- .raparm(parm, a$nm)
  fl <- object[["level"]]
  if (!is.null(level) && !.pnum1(level, .Machine$double.xmin, 1 - 1e-15)) {
    .ckstop("confint", "`level` must be NULL or one number between 0 and 1")
  }
  if (is.null(a$cf[["lo"]])) {
    lv <- if (is.null(level)) 0.95 else level
    z <- stats::qnorm(0.5 + lv / 2)
    lo <- a$cf[["est"]] - z * a$cf[["se"]]
    hi <- a$cf[["est"]] + z * a$cf[["se"]]
  } else {
    if (!is.null(level) && !isTRUE(all.equal(level, fl))) {
      .ckstop("confint", paste0("the intervals of this fit are at level ",
                                "%s; fit again with the `level` wanted"),
              format(fl))
    }
    lv <- fl
    lo <- a$cf[["lo"]]
    hi <- a$cf[["hi"]]
  }
  m <- cbind(lower = lo, upper = hi)
  rownames(m) <- a$nm
  structure(m[i, , drop = FALSE], level = lv)
}

#' vcov() of a regression fit, named by term.
#' @noRd
.ravc <- function(object, ...) {
  a <- .racf(object, "vcov")
  v <- object[["vcov"]]
  if (!is.matrix(v)) {
    .ckstop("vcov", "this fit stores no covariance matrix")
  }
  if (nrow(v) == length(a$nm) && is.null(rownames(v))) {
    dimnames(v) <- list(a$nm, a$nm)
  }
  v
}

#' nobs() of a regression fit.
#' @noRd
.ranobs <- function(object, ...) {
  n <- object[["n"]]
  if (is.null(n)) {
    n <- object[["fit"]][["n"]]
  }
  if (!is.numeric(n) || length(n) != 1L) {
    .ckstop("nobs", "this fit records no count of rows")
  }
  as.integer(n)
}

#' @rdname racc
#' @export
coef.ppmlfe <- .racoef
#' @rdname racc
#' @export
coef.ivfe <- .racoef
#' @rdname racc
#' @export
coef.ppmlcf <- .racoef
#' @rdname racc
#' @export
coef.expgmm <- .racoef
#' @rdname racc
#' @export
coef.fepnl <- .racoef
#' @rdname racc
#' @export
coef.gmmpnl <- .racoef
#' @rdname racc
#' @export
coef.didfe <- .racoef
#' @rdname racc
#' @export
coef.rdloc <- .racoef
#' @rdname racc
#' @export
coef.sarml <- .racoef

#' @rdname racc
#' @export
confint.ppmlfe <- .raci
#' @rdname racc
#' @export
confint.ivfe <- .raci
#' @rdname racc
#' @export
confint.ppmlcf <- .raci
#' @rdname racc
#' @export
confint.expgmm <- .raci
#' @rdname racc
#' @export
confint.fepnl <- .raci
#' @rdname racc
#' @export
confint.gmmpnl <- .raci
#' @rdname racc
#' @export
confint.didfe <- .raci
#' @rdname racc
#' @export
confint.rdloc <- .raci
#' @rdname racc
#' @export
confint.sarml <- .raci

#' @rdname racc
#' @export
vcov.ppmlfe <- .ravc
#' @rdname racc
#' @export
vcov.ivfe <- .ravc
#' @rdname racc
#' @export
vcov.expgmm <- .ravc
#' @rdname racc
#' @export
vcov.fepnl <- .ravc
#' @rdname racc
#' @export
vcov.gmmpnl <- .ravc
#' @rdname racc
#' @export
vcov.didfe <- .ravc
#' @rdname racc
#' @export
vcov.sarml <- .ravc

#' @rdname racc
#' @export
nobs.ppmlfe <- .ranobs
#' @rdname racc
#' @export
nobs.ivfe <- .ranobs
#' @rdname racc
#' @export
nobs.ppmlcf <- .ranobs
#' @rdname racc
#' @export
nobs.expgmm <- .ranobs
#' @rdname racc
#' @export
nobs.fepnl <- .ranobs
#' @rdname racc
#' @export
nobs.gmmpnl <- .ranobs
#' @rdname racc
#' @export
nobs.didfe <- .ranobs
#' @rdname racc
#' @export
nobs.rdloc <- .ranobs
#' @rdname racc
#' @export
nobs.sarml <- .ranobs
