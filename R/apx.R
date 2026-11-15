# Exposure and crime: counts by place and day against a daily exposure
# (air pollution, heat, noise), with the estimators and instruments of
# Bondy, Roth and Sager (2020) and Herrnstadt, Heyes, Muehlegger and
# Saberian (2019, 2021). The numerics are in src/apx.cpp and src/cnt.cpp;
# these wrappers check the arguments and build the results.

# A data frame argument.
.axdat <- function(data, fn, what = "data") {
  if (!is.data.frame(data)) {
    .ckstop(fn, "`%s` must be a data frame, got %s", what, class(data)[1L])
  }
  if (!nrow(data)) .ckstop(fn, "`%s` has no rows", what)
  data
}

# Column names: character, present in d, no missing values in them.
.axcols <- function(d, cols, what, fn, n = NULL, null = FALSE) {
  if (is.null(cols)) {
    if (null) return(NULL)
    .ckstop(fn, "`%s` must name columns of `data`, got NULL", what)
  }
  if (!is.character(cols) || anyNA(cols) || !length(cols)) {
    .ckstop(fn, "`%s` must name columns of `data`", what)
  }
  if (!is.null(n) && !length(cols) %in% n) {
    .ckstop(fn, "`%s` must name %s column(s)", what, paste(n, collapse = " or "))
  }
  miss <- setdiff(cols, names(d))
  if (length(miss)) {
    .ckstop(fn, "`%s` names columns not in `data`: %s", what,
            paste(miss, collapse = ", "))
  }
  if (anyDuplicated(cols)) .ckstop(fn, "`%s` names a column twice", what)
  for (k in cols) {
    if (anyNA(d[[k]])) .ckstop(fn, "column `%s` has missing values", k)
  }
  cols
}

# Numeric columns as a matrix (n x 0 for none), finite.
.axmat <- function(d, cols, fn) {
  if (!length(cols)) return(matrix(numeric(0), nrow(d), 0L))
  for (k in cols) {
    v <- d[[k]]
    if (!is.numeric(v) && !is.logical(v)) {
      .ckstop(fn, "column `%s` must be numeric, got %s", k, class(v)[1L])
    }
    if (!all(is.finite(v))) .ckstop(fn, "column `%s` must be finite", k)
  }
  m <- matrix(as.numeric(unlist(d[cols], use.names = FALSE)), nrow(d))
  colnames(m) <- cols
  m
}

# Codes from 0 of a label column, in order of first appearance.
.axcode <- function(v) {
  if (is.factor(v)) v <- as.character(v)
  u <- unique(v)
  list(c = match(v, u) - 1L, n = length(u))
}

# The effects (one or two columns) as codes.
.axfe <- function(d, fe) {
  a <- .axcode(d[[fe[1L]]])
  b <- if (length(fe) == 2L) .axcode(d[[fe[2L]]]) else list(c = NULL, n = 0L)
  list(c1 = a$c, n1 = a$n, c2 = b$c, n2 = b$n)
}

# The clusters (one or two columns) as codes; refuses one cluster, warns
# under ten, as didfe() does.
.axcl <- function(d, clus, fn) {
  a <- .axcode(d[[clus[1L]]])
  b <- if (length(clus) == 2L) .axcode(d[[clus[2L]]]) else list(c = NULL, n = 0L)
  g <- c(a$n, if (length(clus) == 2L) b$n)
  if (any(g < 2L)) {
    .ckstop(fn, paste0("`clus` has one cluster: a cluster-robust variance ",
                       "needs two or more"))
  }
  if (any(g < 10L)) {
    warning(fn, ": fewer than ten clusters; cluster-robust errors are ",
            "unreliable with so few (wcbx() is made for this case).",
            call. = FALSE)
  }
  list(c1 = a$c, n1 = a$n, c2 = b$c, n2 = b$n, g = g)
}

.axlevel <- function(level, fn) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
        level <= 0 || level >= 1) {
    .ckstop(fn, "`level` must be one number strictly between 0 and 1")
  }
  level
}

# The coefficient table: est, se, stat, p and the interval; t on df
# degrees of freedom, or the normal when df is NULL.
.axtab <- function(term, est, se, level, df = NULL) {
  st <- est / se
  if (is.null(df)) {
    q <- stats::qnorm(0.5 + level / 2)
    p <- 2 * stats::pnorm(-abs(st))
  } else {
    q <- stats::qt(0.5 + level / 2, df)
    p <- 2 * stats::pt(-abs(st), df)
  }
  data.frame(term = term, est = est, se = se, stat = st, p = p,
             lo = est - q * se, hi = est + q * se, row.names = NULL)
}

# The rows kept after separation by the effects, with a warning naming
# how many rows and groups went.
.axsep <- function(y, d, fe, fn) {
  f <- .axfe(d, fe)
  k <- .Call(C_clab_r_fesep, y, f$c1, f$n1, f$c2, f$n2)
  if (!any(k)) .ckstop(fn, "every group of the fixed effects has no counts")
  if (!all(k)) {
    warning(fn, ": ", sum(!k), " row(s) dropped: a group of `",
            paste(fe, collapse = "` or `"), "` with no counts at all ",
            "(its effect would be minus infinity; it carries nothing for the ",
            "other coefficients).", call. = FALSE)
  }
  k
}

#' Poisson pseudo-maximum likelihood with absorbed fixed effects
#'
#' Fits a count model with one or two sets of fixed effects, such as place
#' and day, absorbed rather than estimated as columns, an optional exposure
#' (offset) and cluster-robust standard errors, one-way or two-way. This is
#' the estimator of Bondy, Roth and Sager (2020): daily counts of offences
#' by ward, ward fixed effects and time effects, and ward population as the
#' offset. It serves any count by place and period against any exposure.
#'
#' The model is `E(y | x) = exp(x'b + a_i + g_t + log(expo))`. Poisson
#' pseudo-maximum likelihood (Santos Silva and Tenreyro 2006) needs only
#' that mean to be right: the counts need not be Poisson, their variance
#' need not equal their mean, and zero counts are kept (a regression of
#' the log of the count would have to drop them). The fixed effects take
#' out everything about a place that does not change from day to day, and
#' everything about a day that is common to all places; the coefficients
#' are then identified from changes within a place over days, relative to
#' the other places on the same day. What is left that moves with the
#' exposure from day to day in one place and not the others (traffic, an
#' event, the weather if it is not a regressor) still biases them: that is
#' what the instruments of [ivfe()] and [ppmlcf()] are for.
#'
#' The fit is iteratively reweighted least squares in which each step
#' demeans the working response and the regressors by the fixed effects
#' (alternating projections, until the group means left are below `tol` of
#' each column's largest value), starting from `mu = (y + mean(y)) / 2`,
#' until the relative change of the deviance is below 1e-10. A regressor
#' that does not vary within the fixed effects is refused.
#'
#' Separation. A place (or a day) whose counts are all zero has a fixed
#' effect of minus infinity and says nothing about the coefficients; its
#' rows are dropped with a warning, the rule of 'fixest' and of Stata's
#' ppmlhdfe. Separation by the regressors (a regressor whose value
#' predicts zeros exactly) is not detected; neither paper discusses it.
#'
#' Standard errors are cluster-robust: `(X'WX)^-1 (sum_g s_g s_g') (X'WX)^-1`
#' with the scores `s_g` of each cluster at the final fit, `X` the demeaned
#' regressors, times `G / (G - 1)`; two-way clustering (Cameron, Gelbach and
#' Miller 2011) is `V1 + V2 - V12`, each term with the factor of `G` the
#' smaller number of clusters. This equals `fixest::fepois()` with
#' `ssc(K.adj = FALSE, G.adj = TRUE, G.df = "min")`. Bondy, Roth and Sager
#' cluster by ward and year-month (their table 2).
#'
#' @param data A data frame, one row per place and period.
#' @param y The name of the count column: numbers from 0, not all 0.
#' @param x The names of the regressors: the exposure and any controls
#'   (numeric columns). For `print()`, an object of class `ppmlfe`.
#' @param fe The names of one or two fixed-effect columns (labels), for
#'   example `c("area", "day")`. For several time effects at once (day of
#'   week and year-month), give their combination as one column.
#' @param expo `NULL` or the name of a column of positive exposures
#'   (population): its logarithm is the offset.
#' @param clus `NULL` (the first fixed effect) or the names of one or two
#'   cluster columns.
#' @param level The confidence level of the intervals.
#' @param tol The convergence tolerance of the demeaning.
#' @param maxit The most iterations of the fit.
#' @param object An object of class `ppmlfe`.
#' @param ... Not used.
#' @return An object of class `ppmlfe`, a list: `coef`, a data frame with
#'   one row per regressor (`term`, `est` on the log scale, `se`, `stat`
#'   (z), `p`, `lo`, `hi`, and `rr`, `exp(est)`, the rate ratio of one more
#'   unit); `vcov`, the covariance matrix; `n`, the rows used; `nsep`, the
#'   rows dropped for separation; `ncl`, the number of clusters of each
#'   clustering; `iter`, `conv`, `dev` and `loglik` of the fit; `fe`,
#'   `clus` and `level`. `summary()` returns `coef`, ready to plot.
#' @references
#'   Bondy, M., Roth, S. and Sager, L. (2020). Crime is in the air: the
#'   contemporaneous relationship between air pollution and crime.
#'   *Journal of the Association of Environmental and Resource Economists*,
#'   7(3), 555-585. \doi{10.1086/707127}
#'
#'   Santos Silva, J. M. C. and Tenreyro, S. (2006). The log of gravity.
#'   *Review of Economics and Statistics*, 88(4), 641-658.
#'   \doi{10.1162/rest.88.4.641}
#'
#'   Cameron, A. C., Gelbach, J. B. and Miller, D. L. (2011). Robust
#'   inference with multiway clustering. *Journal of Business and Economic
#'   Statistics*, 29(2), 238-249. \doi{10.1198/jbes.2010.07136}
#' @seealso [ivfe()], [ppmlcf()], [plcday()], [expjn()]
#' @examples
#' # simulated counts in 12 areas over 40 days
#' set.seed(1)
#' d <- expand.grid(day = 1:40, area = 1:12)
#' d$pop <- rep(round(runif(12, 1000, 5000)), each = 40)
#' d$pm <- rnorm(nrow(d), 20, 5)
#' d$temp <- rnorm(nrow(d), 15, 4)
#' d$crimes <- rpois(nrow(d), d$pop / 1000 *
#'   exp(0.02 * d$pm + 0.01 * d$temp + rnorm(12, 0, 0.3)[d$area]))
#' f <- ppmlfe(d, "crimes", c("pm", "temp"), fe = c("area", "day"),
#'             expo = "pop")
#' f
#' # two-way clustering, by area and by day
#' summary(ppmlfe(d, "crimes", c("pm", "temp"), fe = c("area", "day"),
#'                expo = "pop", clus = c("area", "day")))
#' @export
ppmlfe <- function(
  data,
  y,
  x,
  fe,
  expo = NULL,
  clus = NULL,
  level = 0.95,
  tol = 1e-12,
  maxit = 100
) {
  f <- "ppmlfe"
  d <- .axdat(data, f)
  .axcols(d, y, "y", f, n = 1L)
  .axcols(d, x, "x", f)
  .axcols(d, fe, "fe", f, n = 1:2)
  .axcols(d, expo, "expo", f, n = 1L, null = TRUE)
  if (is.null(clus)) clus <- fe[1L]
  .axcols(d, clus, "clus", f, n = 1:2)
  level <- .axlevel(level, f)
  tol <- .cknum(tol, "tol", f, positive = TRUE)
  maxit <- .cknum(maxit, "maxit", f, min = 1, whole = TRUE)
  yv <- .axmat(d, y, f)[, 1L]
  if (any(yv < 0)) .ckstop(f, "column `%s` must hold counts, at least 0", y)
  if (!any(yv > 0)) .ckstop(f, "column `%s` is all 0", y)
  k <- .axsep(yv, d, fe, f)
  d <- d[k, , drop = FALSE]
  yv <- yv[k]
  X <- .axmat(d, x, f)
  off <- NULL
  if (!is.null(expo)) {
    e <- .axmat(d, expo, f)[, 1L]
    if (any(e <= 0)) .ckstop(f, "column `%s` must be positive", expo)
    off <- log(e)
  }
  fx <- .axfe(d, fe)
  cl <- .axcl(d, clus, f)
  r <- .Call(C_clab_r_ppmlfe, yv, X, off, fx$c1, fx$n1, fx$c2, fx$n2,
             cl$c1, cl$n1, cl$c2, cl$n2, tol, as.integer(maxit))
  if (r$st[2L] != 1) {
    warning(f, ": the fit did not converge in ", maxit, " iterations.",
            call. = FALSE)
  }
  v <- r$vc
  dimnames(v) <- list(x, x)
  tab <- .axtab(x, r$beta, sqrt(pmax(diag(v), 0)), level)
  tab$rr <- exp(tab$est)
  structure(
    list(coef = tab, vcov = v, n = length(yv), nsep = sum(!k), ncl = cl$g,
         iter = r$st[1L], conv = r$st[2L] == 1, dev = r$st[3L],
         loglik = r$st[4L], fe = fe, clus = clus, level = level),
    class = "ppmlfe"
  )
}

#' @rdname ppmlfe
#' @export
print.ppmlfe <- function(x, ...) {
  cat("ppmlfe: Poisson pseudo-likelihood, ", x[["n"]], " rows, fixed effects ",
      paste(x[["fe"]], collapse = " and "), "; clusters ",
      paste0(x[["clus"]], " (", x[["ncl"]], ")", collapse = " and "), "\n",
      sep = "")
  if (x[["nsep"]] > 0) {
    cat(x[["nsep"]], " row(s) dropped for groups with no counts\n", sep = "")
  }
  print(x[["coef"]], row.names = FALSE, digits = 4)
  invisible(x)
}

#' @rdname ppmlfe
#' @export
summary.ppmlfe <- function(object, ...) object[["coef"]]
