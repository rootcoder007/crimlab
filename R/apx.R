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

#' Two-stage least squares with absorbed fixed effects
#'
#' Estimates the effect of an exposure that is measured with error or
#' shares unmeasured causes with the outcome, by instrumental variables:
#' each endogenous column is first regressed on the excluded instruments
#' and the controls (the first stage), and the outcome on the fitted values
#' and the controls (the second stage), with one or two sets of fixed
#' effects absorbed throughout. This is the linear estimator of
#' Herrnstadt, Heyes, Muehlegger and Saberian (2019): the log of daily
#' assaults near each monitor on ozone, instrumented by a westerly wind,
#' with monitor and time fixed effects and errors clustered by monitor and
#' by date. Without `endo` and `inst` it is least squares with the fixed
#' effects.
#'
#' What an instrument must satisfy, in plain words. Relevance: it moves the
#' exposure, after the fixed effects and controls; the first-stage
#' statistic `fstat` measures this (a common rule of thumb asks for more
#' than 10; Bondy, Roth and Sager report 13 to 25, Herrnstadt and others
#' 94 to 100). Exclusion: it changes the outcome only through the
#' exposure. Wind direction moves pollution, but if it also moves where
#' people are (a wind off a lake on a warm day) the exclusion fails; no
#' statistic can show that it holds. With effects that differ between
#' places, two-stage least squares estimates a weighted average of them
#' over the places and days the instrument moves, and only if the
#' instrument pushes the exposure the same way everywhere (monotonicity;
#' Bondy, Roth and Sager let wind act differently by region for that
#' reason: see [windiv()]).
#'
#' Standard errors are cluster-robust on the second-stage design with the
#' structural residuals (the outcome less the actual, not the fitted,
#' endogenous columns times their coefficients), times `G / (G - 1)`;
#' two-way as in [ppmlfe()]. `fstat` is the cluster-robust Wald statistic
#' of the excluded instruments in each first stage divided by their
#' number, with the same clusters. These equal `fixest::feols()` with
#' `ssc(K.adj = FALSE, G.adj = TRUE, G.df = "min")` and its `fitstat`
#' `"ivwald1"`. Tests use the t distribution on `G - 1` degrees of
#' freedom.
#'
#' @param data A data frame, one row per place and period.
#' @param y The name of the outcome column (numeric).
#' @param endo `NULL` or the names of the endogenous columns (the
#'   exposure).
#' @param inst `NULL` or the names of the excluded instruments, at least as
#'   many as `endo` (from [windiv()] or [invsn()], say).
#' @param x `NULL` or the names of the exogenous controls. For `print()`,
#'   an object of class `ivfe`.
#' @param fe,clus,level,tol As in [ppmlfe()].
#' @param object An object of class `ivfe`.
#' @param ... Not used.
#' @return An object of class `ivfe`, a list: `coef` (as in [ppmlfe()],
#'   without `rr`); `first`, one row per endogenous column (`endo`,
#'   `fstat`, `ninst`); `vcov`; `n`; `ncl`; `rss`; `fe`, `clus` and
#'   `level`. `summary()` returns `coef`.
#' @references
#'   Herrnstadt, E., Heyes, A., Muehlegger, E. and Saberian, S. (2019). Air
#'   pollution as a cause of violent crime. Working paper.
#'
#'   Bondy, M., Roth, S. and Sager, L. (2020). Crime is in the air: the
#'   contemporaneous relationship between air pollution and crime.
#'   *Journal of the Association of Environmental and Resource Economists*,
#'   7(3), 555-585. \doi{10.1086/707127}
#' @seealso [ppmlcf()] for counts, [windiv()], [invsn()], [plcday()]
#' @examples
#' # simulated: an unmeasured daily shock u raises both the exposure and
#' # the outcome; the wind moves the exposure only
#' set.seed(2)
#' d <- expand.grid(day = 1:60, area = 1:12)
#' d$wind <- runif(nrow(d), 0, 360)
#' w <- windiv(d$wind, width = 90)
#' d <- cbind(d, w)
#' u <- rnorm(nrow(d))
#' d$oz <- 3 * w[, 1] - 2 * w[, 2] + u + rnorm(nrow(d))
#' d$out <- 0.5 * d$oz + u + rnorm(nrow(d))
#' ivfe(d, "out", endo = "oz", inst = colnames(w), fe = c("area", "day"),
#'      clus = c("area", "day"))
#' @export
ivfe <- function(
  data,
  y,
  endo = NULL,
  inst = NULL,
  x = NULL,
  fe,
  clus = NULL,
  level = 0.95,
  tol = 1e-12
) {
  f <- "ivfe"
  d <- .axdat(data, f)
  .axcols(d, y, "y", f, n = 1L)
  .axcols(d, endo, "endo", f, null = TRUE)
  .axcols(d, inst, "inst", f, null = TRUE)
  .axcols(d, x, "x", f, null = TRUE)
  .axcols(d, fe, "fe", f, n = 1:2)
  if (is.null(clus)) clus <- fe[1L]
  .axcols(d, clus, "clus", f, n = 1:2)
  level <- .axlevel(level, f)
  tol <- .cknum(tol, "tol", f, positive = TRUE)
  if (length(inst) < length(endo)) {
    .ckstop(f, "`inst` must give at least as many instruments as `endo` has columns")
  }
  if (!length(endo) && length(inst)) {
    .ckstop(f, "`inst` needs `endo`: the columns it instruments")
  }
  if (!length(endo) && !length(x)) .ckstop(f, "give `endo` or `x`: there is nothing to estimate")
  if (length(intersect(c(endo, x), inst))) {
    .ckstop(f, "a column cannot be both a regressor and an excluded instrument")
  }
  yv <- .axmat(d, y, f)[, 1L]
  fx <- .axfe(d, fe)
  cl <- .axcl(d, clus, f)
  r <- .Call(C_clab_r_ivfe, yv, .axmat(d, endo, f), .axmat(d, x, f),
             .axmat(d, inst, f), fx$c1, fx$n1, fx$c2, fx$n2, cl$c1, cl$n1,
             cl$c2, cl$n2, tol)
  tm <- c(endo, x)
  v <- r$vc
  dimnames(v) <- list(tm, tm)
  g <- min(cl$g)
  structure(
    list(coef = .axtab(tm, r$beta, sqrt(pmax(diag(v), 0)), level, df = g - 1),
         first = data.frame(endo = as.character(endo), fstat = r$fst,
                            ninst = rep(length(inst), length(endo))),
         vcov = v, n = length(yv), ncl = cl$g, rss = r$st[3L], fe = fe,
         clus = clus, level = level),
    class = "ivfe"
  )
}

#' @rdname ivfe
#' @export
print.ivfe <- function(x, ...) {
  iv <- nrow(x[["first"]]) > 0L
  cat("ivfe: ", if (iv) "two-stage least squares" else "least squares", ", ",
      x[["n"]], " rows, fixed effects ", paste(x[["fe"]], collapse = " and "),
      "; clusters ", paste0(x[["clus"]], " (", x[["ncl"]], ")", collapse = " and "),
      "\n", sep = "")
  print(x[["coef"]], row.names = FALSE, digits = 4)
  if (iv) {
    cat("first stage, cluster-robust Wald F of the excluded instruments:\n")
    print(x[["first"]], row.names = FALSE, digits = 4)
  }
  invisible(x)
}

#' @rdname ivfe
#' @export
summary.ivfe <- function(object, ...) object[["coef"]]

# The bootstrap draws: how often each of g clusters is drawn, b draws
# (g x b), from R's generator; seed, when given, is used and the
# generator's state put back.
.axdraw <- function(g, b, seed) {
  if (!is.null(seed)) {
    env <- globalenv()
    had <- exists(".Random.seed", envir = env, inherits = FALSE)
    old <- if (had) get(".Random.seed", envir = env, inherits = FALSE)
    on.exit(if (had) {
      assign(".Random.seed", old, envir = env)
    } else {
      rm(".Random.seed", envir = env)
    })
    set.seed(seed)
  }
  vapply(seq_len(b), function(i) {
    tabulate(sample.int(g, g, replace = TRUE), g)
  }, integer(g))
}

#' The control function for a count outcome with fixed effects
#'
#' Instruments an exposure in the count model of [ppmlfe()]: the first
#' stages regress each endogenous column on the excluded instruments and
#' the controls with the fixed effects (as in [ivfe()]), and their
#' residuals are added to the Poisson pseudo-likelihood fit with the same
#' fixed effects and offset. The residual takes up the part of the
#' exposure that moves with unmeasured causes of the outcome, so the
#' coefficient on the exposure is that of its part the instruments move.
#' Its own coefficient tests whether the exposure was endogenous at all.
#'
#' Bondy, Roth and Sager (2020) instrument pollution in their Poisson
#' model by wind direction and by inversions (their equations 2 and 3);
#' their text writes the second stage with the fitted pollution in place
#' of pollution, and does not say how its standard errors account for the
#' first stage. With a log link, putting the fitted values in place is not
#' consistent in general; the control function is, when the first-stage
#' residual enters the mean as a linear index (Wooldridge 2014). That is
#' the method here, and it is a choice made for crimlab, not a
#' reproduction of the papers' computations.
#'
#' Standard errors come from a cluster bootstrap that refits both stages:
#' `nboot` draws of the clusters with replacement (the first `clus`
#' column), made on the R side with R's generator (so `set.seed()` or
#' `seed` reproduce them), each draw weighting a cluster by how often it
#' was drawn. The interval is the estimate plus or minus the normal
#' quantile times the bootstrap standard error. Separation is handled as
#' in [ppmlfe()], in the full sample and in each draw; a draw that cannot
#' be fitted (an instrument without variation in it, say) is left out and
#' counted in `nfail`. The first-stage statistic is that of [ivfe()],
#' clustered by `clus`.
#'
#' @param data,y,fe,expo,level,tol,maxit As in [ppmlfe()].
#' @param endo The names of the endogenous columns.
#' @param inst The names of the excluded instruments, at least as many.
#' @param x `NULL` or the names of the exogenous controls. For `print()`,
#'   an object of class `ppmlcf`.
#' @param clus `NULL` (the first fixed effect) or the name of one cluster
#'   column: the unit the bootstrap draws.
#' @param nboot The number of bootstrap draws, at least 2.
#' @param seed `NULL` (the generator's current state) or one whole number.
#' @param object An object of class `ppmlcf`.
#' @param ... Not used.
#' @return An object of class `ppmlcf`, a list: `coef` (as in [ppmlfe()],
#'   the endogenous columns, the controls, then the residual terms
#'   `"vhat:<endo>"`); `first` (as in [ivfe()]); `draws`, the bootstrap
#'   estimates (one row per draw); `nfail`; `n`; `nsep`; `ncl`; `nboot`;
#'   `fe`, `clus` and `level`. `summary()` returns `coef`.
#' @references
#'   Bondy, M., Roth, S. and Sager, L. (2020). Crime is in the air: the
#'   contemporaneous relationship between air pollution and crime.
#'   *Journal of the Association of Environmental and Resource Economists*,
#'   7(3), 555-585. \doi{10.1086/707127}
#'
#'   Wooldridge, J. M. (2014). Quasi-maximum likelihood estimation and
#'   testing for nonlinear models with endogenous explanatory variables.
#'   *Journal of Econometrics*, 182(1), 226-234.
#' @seealso [ppmlfe()], [ivfe()], [windiv()], [invsn()]
#' @examples
#' # simulated counts; the shock u raises pollution and crime together
#' set.seed(3)
#' d <- expand.grid(day = 1:50, area = 1:15)
#' d$wind <- runif(nrow(d), 0, 360)
#' w <- windiv(d$wind, width = 90)
#' d <- cbind(d, w)
#' u <- rnorm(nrow(d), 0, 0.5)
#' d$pm <- 2 * w[, 1] - w[, 3] + u + rnorm(nrow(d))
#' d$crimes <- rpois(nrow(d), exp(1 + 0.1 * d$pm + 0.5 * u))
#' ppmlcf(d, "crimes", endo = "pm", inst = colnames(w), fe = c("area", "day"),
#'        nboot = 49, seed = 1)
#' @export
ppmlcf <- function(
  data,
  y,
  endo,
  inst,
  x = NULL,
  fe,
  expo = NULL,
  clus = NULL,
  nboot = 199,
  seed = NULL,
  level = 0.95,
  tol = 1e-12,
  maxit = 100
) {
  f <- "ppmlcf"
  d <- .axdat(data, f)
  .axcols(d, y, "y", f, n = 1L)
  .axcols(d, endo, "endo", f)
  .axcols(d, inst, "inst", f)
  .axcols(d, x, "x", f, null = TRUE)
  .axcols(d, fe, "fe", f, n = 1:2)
  .axcols(d, expo, "expo", f, n = 1L, null = TRUE)
  if (is.null(clus)) clus <- fe[1L]
  .axcols(d, clus, "clus", f, n = 1L)
  level <- .axlevel(level, f)
  tol <- .cknum(tol, "tol", f, positive = TRUE)
  maxit <- .cknum(maxit, "maxit", f, min = 1, whole = TRUE)
  nboot <- .cknum(nboot, "nboot", f, min = 2, whole = TRUE)
  if (!is.null(seed)) seed <- .cknum(seed, "seed", f, min = 0, whole = TRUE)
  if (length(inst) < length(endo)) {
    .ckstop(f, "`inst` must give at least as many instruments as `endo` has columns")
  }
  if (length(intersect(c(endo, x), inst))) {
    .ckstop(f, "a column cannot be both a regressor and an excluded instrument")
  }
  yv <- .axmat(d, y, f)[, 1L]
  if (any(yv < 0)) .ckstop(f, "column `%s` must hold counts, at least 0", y)
  if (!any(yv > 0)) .ckstop(f, "column `%s` is all 0", y)
  k <- .axsep(yv, d, fe, f)
  d <- d[k, , drop = FALSE]
  yv <- yv[k]
  off <- NULL
  if (!is.null(expo)) {
    e <- .axmat(d, expo, f)[, 1L]
    if (any(e <= 0)) .ckstop(f, "column `%s` must be positive", expo)
    off <- log(e)
  }
  fx <- .axfe(d, fe)
  cl <- .axcl(d, clus, f)
  bm <- .axdraw(cl$n1, nboot, seed)
  storage.mode(bm) <- "integer"
  r <- .Call(C_clab_r_ppmlcf, yv, .axmat(d, endo, f), .axmat(d, x, f),
             .axmat(d, inst, f), off, fx$c1, fx$n1, fx$c2, fx$n2, cl$c1,
             cl$n1, bm, tol, as.integer(maxit))
  tm <- c(endo, x, paste0("vhat:", endo))
  dr <- r$draw
  colnames(dr) <- tm
  nfail <- as.integer(r$st[3L])
  if (nfail > 0L) {
    warning(f, ": ", nfail, " of ", nboot, " bootstrap draws could not be ",
            "fitted and are left out.", call. = FALSE)
  }
  if (nboot - nfail < 2L) .ckstop(f, "fewer than two bootstrap draws could be fitted")
  se <- apply(dr, 2L, stats::sd, na.rm = TRUE)
  tab <- .axtab(tm, r$beta, se, level)
  tab$rr <- exp(tab$est)
  structure(
    list(coef = tab,
         first = data.frame(endo = endo, fstat = r$fst,
                            ninst = rep(length(inst), length(endo))),
         draws = dr, nfail = nfail, n = length(yv), nsep = sum(!k),
         ncl = cl$g, nboot = nboot, fe = fe, clus = clus, level = level),
    class = "ppmlcf"
  )
}

#' @rdname ppmlcf
#' @export
print.ppmlcf <- function(x, ...) {
  cat("ppmlcf: Poisson control function, ", x[["n"]], " rows, fixed effects ",
      paste(x[["fe"]], collapse = " and "), "; ", x[["nboot"]],
      " bootstrap draws of ", x[["ncl"]], " clusters (", x[["clus"]], ")\n",
      sep = "")
  print(x[["coef"]], row.names = FALSE, digits = 4)
  cat("first stage, cluster-robust Wald F of the excluded instruments:\n")
  print(x[["first"]], row.names = FALSE, digits = 4)
  invisible(x)
}

#' @rdname ppmlcf
#' @export
summary.ppmlcf <- function(object, ...) object[["coef"]]

#' Wind-direction instruments
#'
#' Turns daily wind directions into indicators of the sector the wind
#' blows from, one sector left out, optionally one set per area so that
#' the same wind may act differently in different parts of a city. These
#' are the instruments of Bondy, Roth and Sager (2020), who use four
#' quadrants (0-90, 91-180, 181-270 degrees; the fourth left out) in each
#' of five regions of London, 15 instruments; and, with `sect`, the single
#' westerly-wind indicator of Herrnstadt, Heyes, Muehlegger and Saberian
#' (2019), wind from between 225 and 285 degrees.
#'
#' Directions are where the wind blows from, in degrees clockwise from
#' north (the meteorological convention); 360 is 0. Sector `b` is `[start
#' + (b - 1) width, start + b width)`, so a direction on a boundary goes to
#' the sector it opens. A daily direction should be the direction of the
#' day's mean wind vector, not the mean of the angles (1 and 359 degrees
#' average to north, not south): `(atan2(mean(sin(r)), mean(cos(r))) * 180
#' / pi) %% 360` for the hourly directions `r` in radians. Bondy, Roth and
#' Sager use instead the share of the day's hours in each quadrant; the
#' daily mean of these indicators over hourly directions gives that share.
#'
#' The indicators are instruments only if the wind direction affects the
#' outcome through the exposure alone (see [ivfe()]).
#'
#' @param dir Wind directions in degrees (numbers; `NA` gives a row of
#'   `NA`).
#' @param width The width of a sector in degrees; it must divide 360.
#' @param start Where the first sector begins, in degrees.
#' @param area `NULL` or the area of each direction (labels): one set of
#'   indicators per area, in sorted order of the labels.
#' @param base The sector left out, counted from 1.
#' @param sect `NULL`, or two numbers `c(from, to)` in degrees: then one
#'   indicator, 1 when the wind blows from between `from` (included) and
#'   `to` (excluded), going clockwise; `width`, `start` and `base` are not
#'   used.
#' @return A numeric matrix, one row per direction, columns named
#'   `"w<from>_<to>"`, with `":<area>"` added when `area` is given.
#' @references
#'   Bondy, M., Roth, S. and Sager, L. (2020). Crime is in the air: the
#'   contemporaneous relationship between air pollution and crime.
#'   *Journal of the Association of Environmental and Resource Economists*,
#'   7(3), 555-585. \doi{10.1086/707127}
#'
#'   Herrnstadt, E., Heyes, A., Muehlegger, E. and Saberian, S. (2019). Air
#'   pollution as a cause of violent crime. Working paper.
#' @seealso [ivfe()], [ppmlcf()], [invsn()]
#' @examples
#' windiv(c(10, 100, 200, 300))
#' # quadrants by region, as in London
#' windiv(c(10, 100, 200, 300), area = c("east", "east", "west", "west"))
#' # one westerly sector
#' windiv(c(230, 300), sect = c(225, 285))
#' @export
windiv <- function(dir, width = 90, start = 0, area = NULL, base = 1, sect = NULL) {
  f <- "windiv"
  if (!is.numeric(dir) || !length(dir)) .ckstop(f, "`dir` must be numeric directions in degrees")
  dir <- as.numeric(dir)
  if (!is.null(sect)) {
    if (!is.numeric(sect) || length(sect) != 2L || !all(is.finite(sect)) ||
          sect[1L] %% 360 == sect[2L] %% 360) {
      .ckstop(f, "`sect` must be two different finite directions in degrees")
    }
    m <- matrix(as.numeric((dir - sect[1L]) %% 360 < (sect[2L] - sect[1L]) %% 360))
    colnames(m) <- paste0("w", format(sect[1L]), "_", format(sect[2L]))
    return(m)
  }
  width <- .cknum(width, "width", f, positive = TRUE)
  start <- .cknum(start, "start", f)
  nb <- round(360 / width)
  if (width >= 360 || abs(360 / width - nb) > 1e-9) {
    .ckstop(f, "`width` must divide 360 into two or more sectors")
  }
  base <- .cknum(base, "base", f, min = 1, max = nb, whole = TRUE)
  ac <- NULL
  na <- 1L
  lev <- NULL
  if (!is.null(area)) {
    if (length(area) != length(dir) || anyNA(area)) {
      .ckstop(f, "`area` must give one label per direction, none missing")
    }
    lev <- if (is.factor(area)) levels(droplevels(area)) else sort(unique(as.character(area)))
    ac <- match(as.character(area), lev) - 1L
    na <- length(lev)
  }
  m <- .Call(C_clab_r_windbin, dir, width, start, ac, na, as.integer(base - 1))
  lo <- (start + (seq_len(nb) - 1) * width) %% 360
  nm <- paste0("w", format(lo, trim = TRUE), "_", format((lo + width) %% 360, trim = TRUE))[-base]
  colnames(m) <- if (is.null(lev)) nm else paste0(rep(nm, na), ":", rep(lev, each = nb - 1L))
  m
}

#' Temperature inversions
#'
#' The instrument of Bondy, Roth and Sager (2020, p. 566-567): "the
#' continuous difference between temperatures at the ground level and
#' temperatures at circa 750 m altitude (925 hectopascals [hPa] pressure
#' level)". On most days the air cools with height; in an inversion it is
#' warmer above than at the ground, and pollutants are trapped near the
#' ground. The paper does not say which temperature is subtracted from
#' which; here `diff` is the upper temperature less the ground one, so it
#' is positive in an inversion, and `inv` is 1 then.
#'
#' Inversions cover large areas (the paper's satellite readings are on a
#' grid of about 110 by 70 km), so they vary little between places in a
#' city; Bondy, Roth and Sager use them at the borough and city level. As
#' an instrument the inversion must affect the outcome only through the
#' exposure; cold, still days may change what people do as well, so
#' weather controls matter (see [ivfe()]).
#'
#' @param tlow Temperatures at the ground (numbers, `NA` allowed).
#' @param tup Temperatures aloft, at the same places and times, in the same
#'   unit.
#' @return A data frame: `diff`, `tup - tlow`, and `inv`, 1 when `diff` is
#'   above 0 and 0 otherwise (`NA` where a temperature is missing).
#' @references
#'   Bondy, M., Roth, S. and Sager, L. (2020). Crime is in the air: the
#'   contemporaneous relationship between air pollution and crime.
#'   *Journal of the Association of Environmental and Resource Economists*,
#'   7(3), 555-585. \doi{10.1086/707127}
#' @seealso [windiv()], [ivfe()], [ppmlcf()]
#' @examples
#' invsn(c(10, 5, 2), c(8, 6, 2))
#' @export
invsn <- function(tlow, tup) {
  f <- "invsn"
  if (!is.numeric(tlow) || !is.numeric(tup) || length(tlow) != length(tup) ||
        !length(tlow)) {
    .ckstop(f, "`tlow` and `tup` must be numeric vectors of one length")
  }
  dd <- as.numeric(tup) - as.numeric(tlow)
  data.frame(diff = dd, inv = as.integer(dd > 0))
}

#' Irrelevant-day placebo
#'
#' Refits an estimator with the exposure taken from other days of the same
#' place, a week, a month or a year before and after, and returns the
#' table of coefficients, the same day first. This is the placebo of Bondy,
#' Roth and Sager (2020, table A1: 7, 31 and 365 days back and forward) and
#' of Herrnstadt, Heyes, Muehlegger and Saberian (2019, table 9: 100 days).
#'
#' What it can show: if pollution on irrelevant days predicts crime as
#' well as pollution on the day does, the estimate is picking up something
#' slow (a season, a trend, a place that is polluted and violent for other
#' reasons) that the fixed effects did not remove. What it cannot show: a
#' cause on the same day that moves with the exposure (a heat wave, an
#' event) passes this test, and so does a wrong instrument; a null placebo
#' is necessary for, not proof of, a causal reading. Shifted exposure is
#' correlated with same-day exposure when pollution persists, so a small
#' placebo coefficient is expected even then.
#'
#' The shift is by the day's value: row `(place, day)` gets the exposure
#' of `(place, day + shift)` and is dropped when that row does not exist;
#' `n` reports the rows left. For the instrumental-variable estimators,
#' Bondy, Roth and Sager do not say whether the instruments are shifted
#' with the exposure; name them in `with` to shift them too.
#'
#' @param data A data frame, one row per place and day.
#' @param var The name of the exposure column to shift.
#' @param unit,day The names of the place and day columns; days are
#'   numbers or dates (one per unit at most).
#' @param shifts The shifts in days, whole numbers other than 0; negative
#'   is the past.
#' @param est The estimator: [ppmlfe()], [ivfe()] or [ppmlcf()], or
#'   another function with a `data` argument returning a list with `coef`
#'   holding a row whose `term` is `var`.
#' @param with `NULL` or the names of other columns to shift with `var`.
#' @param ... The other arguments of `est` (`y`, `x`, `fe`, ...).
#' @param object,x An object of class `plcday` (for the methods).
#' @return An object of class `plcday`, a list: `tab`, a data frame with
#'   one row per shift, 0 (the same day) first (`shift`, `est`, `se`,
#'   `lo`, `hi`, `p`, `n`); `var` and `level`. `summary()` returns `tab`.
#' @references
#'   Bondy, M., Roth, S. and Sager, L. (2020). Crime is in the air: the
#'   contemporaneous relationship between air pollution and crime.
#'   *Journal of the Association of Environmental and Resource Economists*,
#'   7(3), 555-585. \doi{10.1086/707127}
#'
#'   Herrnstadt, E., Heyes, A., Muehlegger, E. and Saberian, S. (2019). Air
#'   pollution as a cause of violent crime. Working paper.
#' @seealso [ppmlfe()], [ivfe()], [ppmlcf()]
#' @examples
#' # simulated: crime responds to same-day pollution only
#' set.seed(4)
#' d <- expand.grid(day = 1:80, area = 1:12)
#' d$pm <- rnorm(nrow(d), 20, 5)
#' d$crimes <- rpois(nrow(d), exp(0.5 + 0.03 * d$pm))
#' plcday(d, "pm", "area", "day", shifts = c(-7, 7), y = "crimes", x = "pm",
#'        fe = c("area", "day"))
#' @export
plcday <- function(
  data,
  var,
  unit,
  day,
  shifts = c(-365, -31, -7, 7, 31, 365),
  est = ppmlfe,
  with = NULL,
  ...
) {
  f <- "plcday"
  d <- .axdat(data, f)
  .axcols(d, var, "var", f, n = 1L)
  .axcols(d, unit, "unit", f, n = 1L)
  .axcols(d, day, "day", f, n = 1L)
  .axcols(d, with, "with", f, null = TRUE)
  if (!is.function(est)) .ckstop(f, "`est` must be an estimator function, such as ppmlfe")
  if (!is.numeric(shifts) || !length(shifts) || !all(is.finite(shifts)) ||
        any(shifts != round(shifts)) || any(shifts == 0) || anyDuplicated(shifts)) {
    .ckstop(f, "`shifts` must be distinct whole numbers of days other than 0")
  }
  tv <- d[[day]]
  if (!is.numeric(tv) && !inherits(tv, "Date")) {
    .ckstop(f, "column `%s` must hold numbers or dates", day)
  }
  tv <- as.numeric(tv)
  uc <- .axcode(d[[unit]])$c
  rows <- lapply(c(0, shifts), function(s) {
    ix <- .Call(C_clab_r_dshift, uc, tv, s)
    k <- !is.na(ix)
    if (sum(k) < 2L) .ckstop(f, "a shift of %s days leaves fewer than two rows", format(s))
    ds <- d[k, , drop = FALSE]
    for (v in c(var, with)) ds[[v]] <- d[[v]][ix[k]]
    r <- est(data = ds, ...)
    cf <- r[["coef"]]
    if (!is.data.frame(cf) || !var %in% cf[["term"]]) {
      .ckstop(f, "`est` returned no coefficient for `%s`", var)
    }
    cf <- cf[cf[["term"]] == var, , drop = FALSE]
    data.frame(shift = s, est = cf[["est"]], se = cf[["se"]], lo = cf[["lo"]],
               hi = cf[["hi"]], p = cf[["p"]], n = sum(k))
  })
  structure(list(tab = do.call(rbind, rows), var = var), class = "plcday")
}

#' @rdname plcday
#' @export
print.plcday <- function(x, ...) {
  cat("plcday: the coefficient of ", x[["var"]],
      " with the exposure of other days (shift 0 is the same day)\n", sep = "")
  print(x[["tab"]], row.names = FALSE, digits = 4)
  invisible(x)
}

#' @rdname plcday
#' @export
summary.plcday <- function(object, ...) object[["tab"]]

#' The upwind and downwind design around a road or a point source
#'
#' Compares counts on the downwind side of a pollution source with counts
#' on its upwind side on the same day, the design of Herrnstadt, Heyes,
#' Muehlegger and Saberian (2021) for the interstates of Chicago (first as
#' Herrnstadt and Muehlegger 2015). The upwind side shares the day's
#' weather and everything else about the day near that road, so it is the
#' control for the downwind side; the side fixed effects take out what
#' differs between the two sides on every day.
#'
#' What a user must supply, for any city and any source: planar
#' coordinates (projected, in one unit such as metres) of each area and of
#' each source; the daily wind direction at each row (where the wind blows
#' from, degrees clockwise from north, the direction of the day's mean
#' wind vector: see [windiv()]); and the count of each area and day. A
#' source is a straight segment (`x1`, `y1` given; split a bending road
#' into straight segments, as the paper splits I-90) or a point (`x1`,
#' `y1` `NA`: a stack, a port, an airport).
#'
#' The classification (crimlab's kernel, following the paper's rules). An
#' area belongs to the one source within `band` of it; areas within `band`
#' of two sources are dropped, as the paper drops crimes within a mile of
#' two interstates. For a segment, the side of the area is set by which
#' side of the line it lies on; the area is downwind on a day when the
#' direction the wind blows to is within `tol` degrees of the line
#' orthogonal to the road pointing to its side, upwind when within `tol`
#' of the opposite direction, and the day is dropped for that road
#' otherwise (wind along the road). The paper's main rule is a band of one
#' mile and 60 degrees; its robustness checks use 30 to 90 degrees and a
#' quarter mile to one mile. For a point source the same rule applies to
#' the direction from the source to the area (downwind within `tol` of
#' where the wind blows to): this extends the design and is not in the
#' papers.
#'
#' The estimate is the paper's equation (3): `y_ist = a_is + g_it + b
#' downwind_ist + e_ist`, with route-side fixed effects `a_is` and
#' route-date fixed effects `g_it`; with `agg = TRUE` the counts of the
#' areas on one side of one segment are summed by day first, as the paper
#' counts crimes by interstate side and day. Without aggregation (and for
#' point sources) the first effect is the area. With `norm = TRUE` the
#' outcome is divided by its mean over the rows used, so that `b` is a
#' share of the mean count, as in the paper; `fam = "ppml"` fits the count
#' model of [ppmlfe()] instead (the paper's appendix table A.1). The paper
#' reports heteroskedasticity-robust errors (here `clus = NULL`: each row
#' its own cluster, with the factor `n / (n - 1)`, not Stata's), and
#' clustered by side and month or year in its appendix (give such a
#' column, constant within a side and day, as `clus`).
#'
#' Placebos. Run the design on an outcome the exposure should not affect
#' (the paper uses property crime); and, with `offs`, move every segment
#' sideways by each offset and re-estimate, the paper's placebo
#' interstates drawn parallel to I-290 a mile apart: the effect should
#' peak at offset 0.
#'
#' @param data A data frame, one row per area and day.
#' @param y The name of the count column.
#' @param ax,ay The names of the columns of the areas' planar coordinates.
#' @param day The name of the day column.
#' @param wind The name of the wind-direction column (degrees, blowing
#'   from).
#' @param src A data frame of sources: columns `x0`, `y0`, and `x1`, `y1`
#'   (`NA` for a point source; the columns may be left out when every
#'   source is a point).
#' @param area `NULL` or the name of the area column (needed for point
#'   sources and for `agg = FALSE`).
#' @param tol The angle in degrees, above 0 and at most 90.
#' @param band The largest distance of an area from its source, in the
#'   coordinates' unit.
#' @param agg `TRUE` to sum the areas on each side of each segment by day.
#' @param fam `"ols"` (least squares, the paper's main model) or `"ppml"`.
#' @param norm `TRUE` to divide the outcome by its mean (least squares
#'   only).
#' @param clus `NULL` (robust by row) or the name of a cluster column,
#'   constant within each side (area) and day.
#' @param offs `NULL` or offsets (in the coordinates' unit, positive to the
#'   left of each segment) for the placebo segments; segments only.
#' @param level The confidence level.
#' @param object,x An object of class `updwn` (for the methods).
#' @param ... Not used.
#' @return An object of class `updwn`, a list: `coef`, the row of the
#'   downwind indicator (as in [ivfe()] or [ppmlfe()]); `class`, the rows
#'   by outcome of the classification (`downwind`, `upwind`, `along`
#'   (wind along the road or crosswind), `nosource`, `twosources`,
#'   `onsource`, `missing`); `n`, the rows of the fit; `fit`, the fit;
#'   `placebo`, `NULL` or one row per offset (`off`, `est`, `se`, `lo`,
#'   `hi`, `p`, `n`). `summary()` returns `coef`, with the placebo rows
#'   when there are any.
#' @references
#'   Herrnstadt, E., Heyes, A., Muehlegger, E. and Saberian, S. (2021). Air
#'   pollution and criminal activity: microgeographic evidence from
#'   Chicago. *American Economic Journal: Applied Economics*, 13(4),
#'   70-100. \doi{10.1257/app.20190091}
#'
#'   Herrnstadt, E. and Muehlegger, E. (2015). Air pollution and criminal
#'   activity: evidence from Chicago microdata. NBER Working Paper 21787.
#'   \doi{10.3386/w21787}
#' @seealso [ppmlfe()], [ivfe()], [gplume()] for the concentration
#'   downwind of a stack
#' @examples
#' # simulated: a road along the x axis, areas on both sides, crime up 10%
#' # on the downwind side
#' set.seed(5)
#' a <- data.frame(area = 1:20, ax = runif(20, 0, 5000),
#'                 ay = rep(c(-800, 800), 10))
#' d <- merge(a, data.frame(day = 1:120))
#' d$wind <- rep(runif(120, 0, 360), each = 1)[d$day]
#' blowsn <- cos((d$wind + 180) * pi / 180)  # +1 when the wind blows north
#' dn <- (d$ay > 0 & blowsn > 0.5) | (d$ay < 0 & blowsn < -0.5)
#' d$crimes <- rpois(nrow(d), 3 * (1 + 0.1 * dn))
#' road <- data.frame(x0 = 0, y0 = 0, x1 = 5000, y1 = 0)
#' updwn(d, "crimes", "ax", "ay", "day", "wind", road, band = 1600,
#'       offs = c(-1600, 0, 1600))
#' @export
updwn <- function(
  data,
  y,
  ax,
  ay,
  day,
  wind,
  src,
  area = NULL,
  tol = 60,
  band = Inf,
  agg = TRUE,
  fam = "ols",
  norm = TRUE,
  clus = NULL,
  offs = NULL,
  level = 0.95
) {
  f <- "updwn"
  d <- .axdat(data, f)
  .axcols(d, y, "y", f, n = 1L)
  .axcols(d, c(ax, ay), "ax, ay", f, n = 2L)
  .axcols(d, day, "day", f, n = 1L)
  .axcols(d, wind, "wind", f, n = 1L)
  .axcols(d, area, "area", f, n = 1L, null = TRUE)
  .axcols(d, clus, "clus", f, n = 1L, null = TRUE)
  s <- .axdat(src, f, "src")
  if (!all(c("x0", "y0") %in% names(s))) .ckstop(f, "`src` must have columns x0 and y0")
  for (k in c("x1", "y1")) if (is.null(s[[k]])) s[[k]] <- NA_real_
  for (k in c("x0", "y0", "x1", "y1")) {
    if (!is.numeric(s[[k]]) && !all(is.na(s[[k]]))) .ckstop(f, "`src` column %s must be numeric", k)
  }
  s0 <- as.numeric(s$x0)
  t0 <- as.numeric(s$y0)
  s1 <- as.numeric(s$x1)
  t1 <- as.numeric(s$y1)
  if (!all(is.finite(s0) & is.finite(t0))) .ckstop(f, "`src` x0 and y0 must be finite")
  line <- is.finite(s1) & is.finite(t1)
  tol <- .cknum(tol, "tol", f, positive = TRUE, min = 0, max = 90)
  if (!identical(band, Inf)) band <- .cknum(band, "band", f, positive = TRUE)
  if (!is.logical(agg) || length(agg) != 1L || is.na(agg)) .ckstop(f, "`agg` must be TRUE or FALSE")
  if (!is.logical(norm) || length(norm) != 1L || is.na(norm)) .ckstop(f, "`norm` must be TRUE or FALSE")
  if (!is.character(fam) || length(fam) != 1L || !fam %in% c("ols", "ppml")) {
    .ckstop(f, "`fam` must be \"ols\" or \"ppml\"")
  }
  if (is.null(area) && (!agg || any(!line))) {
    .ckstop(f, "`area` is needed for point sources and for agg = FALSE")
  }
  if (!is.null(offs)) {
    if (any(!line)) .ckstop(f, "`offs` moves segments: every source must be a segment")
    if (!is.numeric(offs) || !length(offs) || !all(is.finite(offs))) {
      .ckstop(f, "`offs` must be finite numbers")
    }
  }
  level <- .axlevel(level, f)
  yv <- .axmat(d, y, f)[, 1L]
  if (any(yv < 0)) .ckstop(f, "column `%s` must hold counts, at least 0", y)
  xy <- .axmat(d, c(ax, ay), f)
  wd <- .axmat(d, wind, f)[, 1L]
  dv <- as.character(d[[day]])
  av <- if (!is.null(area)) as.character(d[[area]])
  if (!is.null(area) && anyDuplicated(paste(av, dv, sep = "\r"))) {
    .ckstop(f, "an area has two rows on one day")
  }
  labs <- c(`1` = "downwind", `0` = "upwind", `-1` = "along", `-2` = "nosource",
            `-3` = "twosources", `-4` = "onsource", `-5` = "missing")
  run <- function(a0, b0, a1, b1) {
    k <- .Call(C_clab_r_updwn, xy[, 1L], xy[, 2L], wd, a0, b0, a1, b1, band, tol)
    tb <- table(factor(labs[as.character(k$dn)], levels = labs))
    use <- k$dn >= 0L
    if (sum(use) < 2L) .ckstop(f, "fewer than two rows are upwind or downwind of a source")
    ln <- line[k$src[use]]
    g <- ifelse(ln & agg, paste("s", k$src[use], k$side[use], sep = "\r"),
                paste("a", if (is.null(area)) "" else av[use], sep = "\r"))
    cell <- paste(g, dv[use], sep = "\r")
    sd <- paste(k$src[use], dv[use], sep = "\r")
    dn <- k$dn[use]
    if (any(tapply(dn, cell, function(v) length(unique(v))) > 1L)) {
      .ckstop(f, "the wind differs between areas of one side and day: give one direction per day, or agg = FALSE")
    }
    cv <- if (!is.null(clus)) as.character(d[[clus]])[use]
    if (!is.null(clus) && any(tapply(cv, cell, function(v) length(unique(v))) > 1L)) {
      .ckstop(f, "column `%s` must be constant within each side (area) and day", clus)
    }
    o <- !duplicated(cell)
    e <- data.frame(y = as.numeric(rowsum(yv[use], cell, reorder = FALSE)),
                    g = g[o], sd = sd[o], down = as.numeric(dn[o]),
                    cl = if (is.null(clus)) seq_len(sum(o)) else cv[o])
    if (fam == "ols") {
      if (norm) {
        m <- mean(e$y)
        if (!(m > 0)) .ckstop(f, "the counts are all 0")
        e$y <- e$y / m
      }
      r <- ivfe(e, "y", x = "down", fe = c("g", "sd"), clus = "cl", level = level)
    } else {
      r <- ppmlfe(e, "y", "down", fe = c("g", "sd"), clus = "cl", level = level)
    }
    list(fit = r, tab = tb, n = nrow(e))
  }
  r0 <- run(s0, t0, s1, t1)
  pl <- NULL
  if (!is.null(offs)) {
    len <- sqrt((s1 - s0)^2 + (t1 - t0)^2)
    nx <- -(t1 - t0) / len
    ny <- (s1 - s0) / len
    pl <- do.call(rbind, lapply(offs, function(o) {
      ro <- run(s0 + o * nx, t0 + o * ny, s1 + o * nx, t1 + o * ny)
      c0 <- ro$fit$coef
      data.frame(off = o, est = c0$est, se = c0$se, lo = c0$lo, hi = c0$hi,
                 p = c0$p, n = ro$n)
    }))
  }
  structure(
    list(coef = r0$fit$coef, class = r0$tab, n = r0$n, fit = r0$fit,
         placebo = pl, fam = fam, tol = tol, band = band),
    class = "updwn"
  )
}

#' @rdname updwn
#' @export
print.updwn <- function(x, ...) {
  cat("updwn: downwind against upwind, ", x[["n"]], " side (area) days, ",
      if (x[["fam"]] == "ols") "least squares" else "Poisson pseudo-likelihood",
      ", tolerance ", format(x[["tol"]]), " degrees\n", sep = "")
  print(x[["coef"]], row.names = FALSE, digits = 4)
  cat("rows by classification:\n")
  print(x[["class"]])
  if (!is.null(x[["placebo"]])) {
    cat("placebo segments moved sideways:\n")
    print(x[["placebo"]], row.names = FALSE, digits = 4)
  }
  invisible(x)
}

#' @rdname updwn
#' @export
summary.updwn <- function(object, ...) {
  cf <- object[["coef"]]
  if (is.null(object[["placebo"]])) return(cf)
  cbind(off = c(NA, object[["placebo"]][["off"]]),
        rbind(cf[c("est", "se", "lo", "hi", "p")],
              object[["placebo"]][c("est", "se", "lo", "hi", "p")]))
}

#' Join an exposure table to counts by place and day
#'
#' Makes the analysis table of the estimators of this family: the rows of
#' the outcome table (place, day, count, ...) with the exposure of the same
#' place and day joined on. Exposure may come from monitors (assigned to
#' areas beforehand, for example by inverse squared distance from the three
#' nearest, as Bondy, Roth and Sager do) or from crimlab's dispersion
#' functions ([gplume()], [plmfld()]).
#'
#' Duplicate keys in either table are refused: a place with two exposures
#' on one day has no single value to join. Outcome rows without an
#' exposure are left out of `data` and counted in `nout` (their keys in
#' `unout`); exposure rows without an outcome are counted in `nexp`. Keys
#' are compared as text, so a day given as a date in one table must be a
#' date in the other.
#'
#' @param expo The exposure table: a data frame with the key columns and
#'   one or more value columns.
#' @param outc The outcome table: a data frame with the key columns.
#' @param by The names of the key columns, in both tables.
#' @param object,x An object of class `expjn`.
#' @param ... Not used.
#' @return An object of class `expjn`, a list: `data`, the joined table
#'   (the outcome table's columns, then the exposure's other columns);
#'   `nout`, `nexp`; `unout`, the keys of the outcome rows left out.
#'   `summary()` returns the counts as a data frame.
#' @seealso [ppmlfe()], [plmfld()]
#' @examples
#' ex <- data.frame(area = c("a", "a", "b"), day = c(1, 2, 1), pm = c(10, 12, 8))
#' oc <- data.frame(area = c("a", "a", "b", "b"), day = c(1, 2, 1, 2),
#'                  crimes = c(3, 1, 0, 2))
#' j <- expjn(ex, oc)
#' j
#' j$data
#' @export
expjn <- function(expo, outc, by = c("area", "day")) {
  f <- "expjn"
  e <- .axdat(expo, f, "expo")
  o <- .axdat(outc, f, "outc")
  if (!is.character(by) || !length(by) || anyNA(by) || anyDuplicated(by)) {
    .ckstop(f, "`by` must name the key columns")
  }
  for (nm in c("expo", "outc")) {
    t <- if (nm == "expo") e else o
    miss <- setdiff(by, names(t))
    if (length(miss)) .ckstop(f, "`%s` lacks the key columns: %s", nm, paste(miss, collapse = ", "))
    if (anyNA(t[by])) .ckstop(f, "`%s` has missing keys", nm)
  }
  val <- setdiff(names(e), by)
  if (!length(val)) .ckstop(f, "`expo` has no columns besides the keys")
  clash <- intersect(val, names(o))
  if (length(clash)) {
    .ckstop(f, "columns in both tables besides the keys: %s", paste(clash, collapse = ", "))
  }
  key <- function(t) do.call(paste, c(lapply(t[by], as.character), sep = "\r"))
  ke <- key(e)
  ko <- key(o)
  if (anyDuplicated(ke)) .ckstop(f, "`expo` has %d duplicate key(s)", sum(duplicated(ke)))
  if (anyDuplicated(ko)) .ckstop(f, "`outc` has %d duplicate key(s)", sum(duplicated(ko)))
  m <- match(ko, ke)
  hit <- !is.na(m)
  if (!any(hit)) .ckstop(f, "no outcome row has an exposure: check that the keys agree")
  dd <- cbind(o[hit, , drop = FALSE], e[m[hit], val, drop = FALSE])
  rownames(dd) <- NULL
  un <- o[!hit, by, drop = FALSE]
  rownames(un) <- NULL
  structure(list(data = dd, nout = sum(!hit), nexp = sum(!ke %in% ko), unout = un),
            class = "expjn")
}

#' @rdname expjn
#' @export
print.expjn <- function(x, ...) {
  cat("expjn: ", nrow(x[["data"]]), " rows joined; ", x[["nout"]],
      " outcome row(s) without an exposure left out; ", x[["nexp"]],
      " exposure row(s) without an outcome\n", sep = "")
  invisible(x)
}

#' @rdname expjn
#' @export
summary.expjn <- function(object, ...) {
  data.frame(joined = nrow(object[["data"]]), nout = object[["nout"]],
             nexp = object[["nexp"]])
}
