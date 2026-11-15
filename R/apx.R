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

# Whether each level of f falls in one level of c (the effect is nested
# in the clusters).
.axnest <- function(f, c) {
  u <- unique(data.frame(a = as.character(f), b = as.character(c)))
  !anyDuplicated(u$a)
}

# A time column as numbers (numbers or dates).
.axtime <- function(d, time, fn) {
  if (is.null(time)) .ckstop(fn, "`se = \"nw\"` needs `time`, the name of the day column")
  .axcols(d, time, "time", fn, n = 1L)
  tv <- d[[time]]
  if (!is.numeric(tv) && !inherits(tv, "Date")) {
    .ckstop(fn, "column `%s` must hold numbers or dates", time)
  }
  as.numeric(tv)
}

# The variance of a fit from its rows' scores (r$sco, r$brd), as `se`
# asks: cluster (kind 0), twoway (1), hc1 (2), nw (3); the small-sample
# factor by `ssc`. Returns the matrix, the cluster counts and the lag.
.axvar <- function(r, d, se, clus, time, lag, ssc, fe, fn) {
  n <- nrow(r$sco)
  p <- ncol(r$sco)
  kk <- function(cls) {
    k <- p + 1
    for (f in fe) {
      if (!any(vapply(cls, function(cc) .axnest(d[[f]], d[[cc]]), NA))) {
        k <- k + length(unique(d[[f]])) - 1
      }
    }
    k
  }
  c1 <- c2 <- tv <- NULL
  n1 <- n2 <- 0L
  lg <- 0L
  if (se %in% c("cluster", "twoway")) {
    cl <- .axcl(d, clus, fn)
    c1 <- cl$c1
    n1 <- cl$n1
    c2 <- cl$c2
    n2 <- cl$n2
    g <- cl$g
    kind <- if (se == "twoway") 1L else 0L
    k <- kk(clus)
  } else {
    g <- n
    kind <- if (se == "hc1") 2L else 3L
    k <- kk(character())
    if (se == "nw") {
      u <- .axcode(d[[clus[1L]]])
      c1 <- u$c
      n1 <- u$n
      tv <- .axtime(d, time, fn)
      if (is.null(lag)) {
        tt <- max(tabulate(u$c + 1L, u$n))
        lag <- floor(4 * (tt / 100)^(2 / 9))
      }
      lg <- .cknum(lag, "lag", fn, min = 0, whole = TRUE)
    }
  }
  gm <- min(g)
  fac <- switch(ssc, G = gm / (gm - 1), nk = gm / (gm - 1) * (n - 1) / (n - k), none = 1)
  v <- .Call(C_clab_r_vcx, r$sco, r$brd, kind, c1, n1, c2, n2, tv, as.integer(lg), fac)
  list(v = v, g = g, lag = if (se == "nw") lg else NA_integer_, k = k)
}

# The cluster bootstrap: nboot draws of the clusters of column cl from R's
# generator; fit(rows, w) returns the coefficients of the rows drawn with
# their weights, or an error for a draw that cannot be fitted (left out).
.axboot <- function(d, cl, nboot, seed, p, fit, fn) {
  cc <- .axcode(d[[cl]])
  if (cc$n < 2L) .ckstop(fn, "the bootstrap needs two or more clusters in `%s`", cl)
  bm <- .axdraw(cc$n, nboot, seed)
  dr <- t(vapply(seq_len(nboot), function(b) {
    m <- bm[cc$c + 1L, b]
    k <- m > 0
    tryCatch(fit(which(k), m[k]), error = function(e) rep(NA_real_, p))
  }, numeric(p)))
  nf <- sum(!stats::complete.cases(dr))
  if (nf > 0L) {
    warning(fn, ": ", nf, " of ", nboot, " bootstrap draws could not be fitted and are ",
            "left out.", call. = FALSE)
  }
  if (nboot - nf < 2L) .ckstop(fn, "fewer than two bootstrap draws could be fitted")
  list(draws = dr, se = apply(dr, 2L, stats::sd, na.rm = TRUE), nfail = nf, g = cc$n)
}

# The Poisson fit on the rows of d (codes made here), with weights w.
.axppfit <- function(d, y, x, fe, expo, w, tol, maxit, fn) {
  fx <- .axfe(d, fe)
  off <- if (!is.null(expo)) log(d[[expo]])
  .Call(C_clab_r_ppmlfe, as.numeric(d[[y]]), .axmat(d, x, fn), off, w, fx$c1, fx$n1,
        fx$c2, fx$n2, NULL, 0L, NULL, 0L, tol, as.integer(maxit))
}

.axsen <- function(se) {
  c(cluster = "clustered", twoway = "two-way clustered", hc1 = "heteroskedasticity-robust",
    nw = "Newey-West", boot = "cluster bootstrap")[[se]]
}

#' Poisson pseudo-maximum likelihood with absorbed fixed effects
#'
#' Fits a count model with one or two sets of fixed effects, such as place
#' and day, absorbed rather than estimated as columns, an optional exposure
#' (offset) and robust standard errors of several kinds. This is the
#' estimator of Bondy, Roth and Sager (2020): daily counts of offences by
#' ward, ward fixed effects and time effects, and ward population as the
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
#' Separation (`sep`). A place (or a day) whose counts are all zero has a
#' fixed effect of minus infinity and says nothing about the coefficients;
#' `sep = "fe"` (the default) drops its rows with a warning, the rule of
#' 'fixest' and of Stata's ppmlhdfe; `sep = "none"` keeps them (the fit
#' then converges slowly towards the boundary). Separation by the
#' regressors (a regressor whose value predicts zeros exactly) is not
#' detected: the iterative rectifier of Correia, Guimaraes and Zylkin that
#' finds it is not offered, as its paper was not at hand to build it from.
#' Neither Bondy, Roth and Sager nor Herrnstadt and others discuss
#' separation.
#'
#' Standard errors (`se`), all from the rows' scores `s_i` at the final fit
#' and the bread `(X'WX)^-1` of the demeaned regressors `X`:
#' `"cluster"` (the default; `sum_g s_g s_g'` over the clusters of the
#' first `clus` column); `"twoway"` (Cameron, Gelbach and Miller 2011, `V1
#' + V2 - V12`; Bondy, Roth and Sager cluster by ward and year-month, table
#' 2; Herrnstadt and others 2019 by monitor and date; given two `clus`
#' columns, `"cluster"` means this too); `"hc1"` (each row alone, the
#' robust errors of Herrnstadt and others 2021); `"nw"` (Newey-West within
#' each unit of the first `clus` column over the days of `time`, Bartlett
#' weights `1 - d / (lag + 1)` for days `d` apart up to `lag`, by default
#' `floor(4 (T / 100)^(2/9))` for `T` days (Newey and West 1994); the
#' appendix of Herrnstadt and others 2021); `"boot"` (a cluster bootstrap
#' over the first `clus` column, `nboot` draws from R's generator, each
#' refitting the model). The small-sample factor `ssc`: `"G"` (the default)
#' `G / (G - 1)` with `G` the (smaller) number of clusters, or `n / (n - 1)`
#' for `"hc1"` and `"nw"`; `"nk"` that times `(n - 1) / (n - K)`, `K` the
#' coefficients, an intercept and the levels less one of each fixed effect
#' not nested in a cluster (fixest's default for clustered errors); `"none"`
#' no factor. `"cluster"` with `"G"` equals `fixest::fepois()` with
#' `ssc(K.adj = FALSE, G.adj = TRUE, G.df = "min")`.
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
#' @param clus `NULL` (the fixed effects: the first, or both for
#'   `se = "twoway"`) or the names of one or two cluster columns.
#' @param se The standard errors: `"cluster"`, `"twoway"`, `"hc1"`, `"nw"`
#'   or `"boot"`.
#' @param ssc The small-sample factor: `"G"`, `"nk"` or `"none"`.
#' @param time The name of the day column (numbers or dates), for
#'   `se = "nw"`.
#' @param lag `NULL` (the rule above) or the Newey-West lag in days.
#' @param nboot,seed The bootstrap draws (at least 2) and `NULL` or a whole
#'   number for R's generator, for `se = "boot"`.
#' @param sep `"fe"` or `"none"`: see Separation.
#' @param level The confidence level of the intervals.
#' @param tol The convergence tolerance of the demeaning.
#' @param maxit The most iterations of the fit.
#' @param object An object of class `ppmlfe`.
#' @param ... Not used.
#' @return An object of class `ppmlfe`, a list: `coef`, a data frame with
#'   one row per regressor (`term`, `est` on the log scale, `se`, `stat`
#'   (z), `p`, `lo`, `hi`, and `rr`, `exp(est)`, the rate ratio of one more
#'   unit); `vcov`; `n`, the rows used; `nsep`, the rows dropped for
#'   separation; `ncl`, the number of clusters (rows for `"hc1"` and
#'   `"nw"`); `se`, `ssc`, `lag`; `iter`, `conv`, `dev` and `loglik` of the
#'   fit; `draws` for the bootstrap; `fe`, `clus` and `level`. `summary()`
#'   returns `coef`, ready to plot.
#' @references
#'   Bondy, M., Roth, S. and Sager, L. (2020). Crime is in the air: the
#'   contemporaneous relationship between air pollution and crime.
#'   *Journal of the Association of Environmental and Resource Economists*,
#'   7(3), 555-585. \doi{10.1086/707127}
#'
#'   Herrnstadt, E., Heyes, A., Muehlegger, E. and Saberian, S. (2021). Air
#'   pollution and criminal activity: microgeographic evidence from
#'   Chicago. *American Economic Journal: Applied Economics*, 13(4),
#'   70-100. \doi{10.1257/app.20190091}
#'
#'   Santos Silva, J. M. C. and Tenreyro, S. (2006). The log of gravity.
#'   *Review of Economics and Statistics*, 88(4), 641-658.
#'   \doi{10.1162/rest.88.4.641}
#'
#'   Cameron, A. C., Gelbach, J. B. and Miller, D. L. (2011). Robust
#'   inference with multiway clustering. *Journal of Business and Economic
#'   Statistics*, 29(2), 238-249. \doi{10.1198/jbes.2010.07136}
#'
#'   Newey, W. K. and West, K. D. (1994). Automatic lag selection in
#'   covariance matrix estimation. *Review of Economic Studies*, 61(4),
#'   631-653. \doi{10.2307/2297912}
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
#' # two-way clustering by area and day; Newey-West within areas
#' summary(ppmlfe(d, "crimes", c("pm", "temp"), fe = c("area", "day"),
#'                expo = "pop", se = "twoway"))
#' summary(ppmlfe(d, "crimes", c("pm", "temp"), fe = c("area", "day"),
#'                expo = "pop", se = "nw", time = "day", lag = 3))
#' @export
ppmlfe <- function(
  data,
  y,
  x,
  fe,
  expo = NULL,
  clus = NULL,
  se = c("cluster", "twoway", "hc1", "nw", "boot"),
  ssc = c("G", "nk", "none"),
  time = NULL,
  lag = NULL,
  nboot = 199,
  seed = NULL,
  sep = c("fe", "none"),
  level = 0.95,
  tol = 1e-12,
  maxit = 100
) {
  f <- "ppmlfe"
  se <- match.arg(se)
  ssc <- match.arg(ssc)
  sep <- match.arg(sep)
  d <- .axdat(data, f)
  .axcols(d, y, "y", f, n = 1L)
  .axcols(d, x, "x", f)
  .axcols(d, fe, "fe", f, n = 1:2)
  .axcols(d, expo, "expo", f, n = 1L, null = TRUE)
  if (is.null(clus)) clus <- if (se == "twoway") fe else fe[1L]
  .axcols(d, clus, "clus", f, n = 1:2)
  if (se == "cluster" && length(clus) == 2L) se <- "twoway"
  if (se == "twoway" && length(clus) != 2L) .ckstop(f, "`se = \"twoway\"` needs two `clus` columns")
  level <- .axlevel(level, f)
  tol <- .cknum(tol, "tol", f, positive = TRUE)
  maxit <- .cknum(maxit, "maxit", f, min = 1, whole = TRUE)
  yv <- .axmat(d, y, f)[, 1L]
  if (any(yv < 0)) .ckstop(f, "column `%s` must hold counts, at least 0", y)
  if (!any(yv > 0)) .ckstop(f, "column `%s` is all 0", y)
  k <- if (sep == "fe") .axsep(yv, d, fe, f) else rep(TRUE, length(yv))
  d <- d[k, , drop = FALSE]
  .axmat(d, x, f)
  if (!is.null(expo) && any(.axmat(d, expo, f)[, 1L] <= 0)) {
    .ckstop(f, "column `%s` must be positive", expo)
  }
  r <- .axppfit(d, y, x, fe, expo, NULL, tol, maxit, f)
  if (r$st[2L] != 1) {
    warning(f, ": the fit did not converge in ", maxit, " iterations.", call. = FALSE)
  }
  dr <- NULL
  if (se == "boot") {
    nboot <- .cknum(nboot, "nboot", f, min = 2, whole = TRUE)
    if (!is.null(seed)) seed <- .cknum(seed, "seed", f, min = 0, whole = TRUE)
    b <- .axboot(d, clus[1L], nboot, seed, length(x), function(rows, w) {
      ds <- d[rows, , drop = FALSE]
      kk <- if (sep == "fe") {
        fx <- .axfe(ds, fe)
        .Call(C_clab_r_fesep, as.numeric(ds[[y]]), fx$c1, fx$n1, fx$c2, fx$n2)
      } else {
        rep(TRUE, nrow(ds))
      }
      .axppfit(ds[kk, , drop = FALSE], y, x, fe, expo, as.numeric(w[kk]), tol, maxit, f)$beta
    }, f)
    v <- stats::cov(b$draws, use = "complete.obs")
    sev <- b$se
    g <- b$g
    lg <- NA_integer_
    dr <- b$draws
    colnames(dr) <- x
  } else {
    vv <- .axvar(r, d, se, clus, time, lag, ssc, fe, f)
    v <- vv$v
    sev <- sqrt(pmax(diag(v), 0))
    g <- vv$g
    lg <- vv$lag
  }
  dimnames(v) <- list(x, x)
  tab <- .axtab(x, r$beta, sev, level)
  tab$rr <- exp(tab$est)
  structure(
    list(coef = tab, vcov = v, n = nrow(d), nsep = sum(!k), ncl = g, se = se, ssc = ssc,
         lag = lg, iter = r$st[1L], conv = r$st[2L] == 1, dev = r$st[3L], loglik = r$st[4L],
         draws = dr, fe = fe, clus = clus, level = level),
    class = "ppmlfe"
  )
}

#' @rdname ppmlfe
#' @export
print.ppmlfe <- function(x, ...) {
  cat("ppmlfe: Poisson pseudo-likelihood, ", x[["n"]], " rows, fixed effects ",
      paste(x[["fe"]], collapse = " and "), "; ", .axsen(x[["se"]]), " errors (",
      paste0(x[["clus"]], collapse = ", "), "; ", paste(x[["ncl"]], collapse = " and "),
      ")\n", sep = "")
  if (x[["nsep"]] > 0) {
    cat(x[["nsep"]], " row(s) dropped for groups with no counts\n", sep = "")
  }
  print(x[["coef"]], row.names = FALSE, digits = 4)
  invisible(x)
}

#' @rdname ppmlfe
#' @export
summary.ppmlfe <- function(object, ...) object[["coef"]]

# The ivfe kernel on the rows of d with weights w; first-stage clusters
# fc (codes) as the variance kind asks.
.axivfit <- function(d, y, endo, inst, x, fe, w, fc, tol, fn) {
  fx <- .axfe(d, fe)
  .Call(C_clab_r_ivfe, as.numeric(d[[y]]), .axmat(d, endo, fn), .axmat(d, x, fn),
        .axmat(d, inst, fn), w, fx$c1, fx$n1, fx$c2, fx$n2, fc$c1, fc$n1, fc$c2, fc$n2, tol)
}

# The clusters of the first-stage statistics: those of the variance for
# "cluster" and "twoway", each row for "hc1", the first clus column for
# "nw" and "boot".
.axfc <- function(d, se, clus, fn) {
  if (se == "hc1") return(list(c1 = seq_len(nrow(d)) - 1L, n1 = nrow(d), c2 = NULL, n2 = 0L))
  a <- .axcode(d[[clus[1L]]])
  b <- if (se == "twoway") .axcode(d[[clus[2L]]]) else list(c = NULL, n = 0L)
  if (a$n < 2L || (se == "twoway" && b$n < 2L)) {
    .ckstop(fn, paste0("`clus` has one cluster: a cluster-robust variance ",
                       "needs two or more"))
  }
  list(c1 = a$c, n1 = a$n, c2 = b$c, n2 = b$n)
}

.axfirst <- function(endo, fst, q) {
  m <- matrix(fst, ncol = 3L)
  data.frame(endo = as.character(endo), wald = m[, 1L], kp = m[, 2L], eff = m[, 3L],
             ninst = rep(q, length(endo)))
}

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
#' exposure, after the fixed effects and controls. Exclusion: it changes
#' the outcome only through the exposure. Wind direction moves pollution,
#' but if it also moves where people are (a wind off a lake on a warm day)
#' the exclusion fails; no statistic can show that it holds. With effects
#' that differ between places, two-stage least squares estimates a
#' weighted average of them over the places and days the instrument moves,
#' and only if the instrument pushes the exposure the same way everywhere
#' (monotonicity; Bondy, Roth and Sager let wind act differently by region
#' for that reason: see [windiv()]).
#'
#' First-stage strength, all three in `first`: `wald`, the robust Wald
#' statistic of the excluded instruments divided by their number (the
#' "F-test" Bondy, Roth and Sager report, as far as their text says);
#' `kp`, the Kleibergen-Paap rk Wald F of Herrnstadt and others (2019,
#' table 7), which for one endogenous column is the same robust Wald F (it
#' is `NA` for more than one: the rank statistic is not computed); and
#' `eff`, the effective F of Montiel Olea and Pflueger (2013), `b'Qb /
#' tr(V Q)` with `b` the instruments' first-stage coefficients, `V` their
#' robust covariance and `Q` their cross products after the controls and
#' the fixed effects (one endogenous column only). A common rule of thumb
#' asks for more than 10; Montiel Olea and Pflueger give critical values
#' for the effective F. The first-stage statistics use the clusters of
#' `se` (each row for `"hc1"`, the first `clus` column for `"nw"` and
#' `"boot"`).
#'
#' Standard errors (`se`, `ssc`, `time`, `lag`, `nboot`, `seed`) are as in
#' [ppmlfe()], on the second-stage design with the structural residuals
#' (the outcome less the actual, not the fitted, endogenous columns times
#' their coefficients). With `"cluster"` and `ssc = "G"` they equal
#' `fixest::feols()` with `ssc(K.adj = FALSE, G.adj = TRUE, G.df = "min")`.
#' Tests use the t distribution on `G - 1` degrees of freedom for
#' clustered errors and the normal otherwise.
#'
#' @param data A data frame, one row per place and period.
#' @param y The name of the outcome column (numeric).
#' @param endo `NULL` or the names of the endogenous columns (the
#'   exposure).
#' @param inst `NULL` or the names of the excluded instruments, at least as
#'   many as `endo` (from [windiv()] or [invsn()], say).
#' @param x `NULL` or the names of the exogenous controls. For `print()`,
#'   an object of class `ivfe`.
#' @param fe,clus,se,ssc,time,lag,nboot,seed,level,tol As in [ppmlfe()].
#' @param object An object of class `ivfe`.
#' @param ... Not used.
#' @return An object of class `ivfe`, a list: `coef` (as in [ppmlfe()],
#'   without `rr`); `first`, one row per endogenous column (`endo`, `wald`,
#'   `kp`, `eff`, `ninst`); `vcov`; `n`; `ncl`; `se`, `ssc`, `lag`; `rss`;
#'   `draws` for the bootstrap; `fe`, `clus` and `level`. `summary()`
#'   returns `coef`.
#' @references
#'   Herrnstadt, E., Heyes, A., Muehlegger, E. and Saberian, S. (2019). Air
#'   pollution as a cause of violent crime. Working paper.
#'
#'   Bondy, M., Roth, S. and Sager, L. (2020). Crime is in the air: the
#'   contemporaneous relationship between air pollution and crime.
#'   *Journal of the Association of Environmental and Resource Economists*,
#'   7(3), 555-585. \doi{10.1086/707127}
#'
#'   Kleibergen, F. and Paap, R. (2006). Generalized reduced rank tests
#'   using the singular value decomposition. *Journal of Econometrics*,
#'   133(1), 97-126. \doi{10.1016/j.jeconom.2005.02.011}
#'
#'   Montiel Olea, J. L. and Pflueger, C. (2013). A robust test for weak
#'   instruments. *Journal of Business and Economic Statistics*, 31(3),
#'   358-369. \doi{10.1080/00401706.2013.806694}
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
#' f <- ivfe(d, "out", endo = "oz", inst = colnames(w), fe = c("area", "day"),
#'           se = "twoway")
#' f
#' f$first
#' @export
ivfe <- function(
  data,
  y,
  endo = NULL,
  inst = NULL,
  x = NULL,
  fe,
  clus = NULL,
  se = c("cluster", "twoway", "hc1", "nw", "boot"),
  ssc = c("G", "nk", "none"),
  time = NULL,
  lag = NULL,
  nboot = 199,
  seed = NULL,
  level = 0.95,
  tol = 1e-12
) {
  f <- "ivfe"
  se <- match.arg(se)
  ssc <- match.arg(ssc)
  d <- .axdat(data, f)
  .axcols(d, y, "y", f, n = 1L)
  .axcols(d, endo, "endo", f, null = TRUE)
  .axcols(d, inst, "inst", f, null = TRUE)
  .axcols(d, x, "x", f, null = TRUE)
  .axcols(d, fe, "fe", f, n = 1:2)
  if (is.null(clus)) clus <- if (se == "twoway") fe else fe[1L]
  .axcols(d, clus, "clus", f, n = 1:2)
  if (se == "cluster" && length(clus) == 2L) se <- "twoway"
  if (se == "twoway" && length(clus) != 2L) .ckstop(f, "`se = \"twoway\"` needs two `clus` columns")
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
  .axmat(d, c(y, endo, x, inst), f)
  fc <- .axfc(d, se, clus, f)
  r <- .axivfit(d, y, endo, inst, x, fe, NULL, fc, tol, f)
  tm <- c(endo, x)
  dr <- NULL
  df <- NULL
  if (se == "boot") {
    nboot <- .cknum(nboot, "nboot", f, min = 2, whole = TRUE)
    if (!is.null(seed)) seed <- .cknum(seed, "seed", f, min = 0, whole = TRUE)
    b <- .axboot(d, clus[1L], nboot, seed, length(tm), function(rows, w) {
      .axivfit(d[rows, , drop = FALSE], y, endo, inst, x, fe, as.numeric(w),
               list(c1 = NULL, n1 = 0L, c2 = NULL, n2 = 0L), tol, f)$beta
    }, f)
    v <- stats::cov(b$draws, use = "complete.obs")
    sev <- b$se
    g <- b$g
    lg <- NA_integer_
    dr <- b$draws
    colnames(dr) <- tm
  } else {
    vv <- .axvar(r, d, se, clus, time, lag, ssc, fe, f)
    v <- vv$v
    sev <- sqrt(pmax(diag(v), 0))
    g <- vv$g
    lg <- vv$lag
    if (se %in% c("cluster", "twoway")) df <- min(g) - 1
  }
  dimnames(v) <- list(tm, tm)
  structure(
    list(coef = .axtab(tm, r$beta, sev, level, df = df),
         first = .axfirst(endo, r$fst, length(inst)), vcov = v, n = nrow(d), ncl = g,
         se = se, ssc = ssc, lag = lg, rss = r$st[3L], draws = dr, fe = fe, clus = clus,
         level = level),
    class = "ivfe"
  )
}

#' @rdname ivfe
#' @export
print.ivfe <- function(x, ...) {
  iv <- nrow(x[["first"]]) > 0L
  cat("ivfe: ", if (iv) "two-stage least squares" else "least squares", ", ",
      x[["n"]], " rows, fixed effects ", paste(x[["fe"]], collapse = " and "), "; ",
      .axsen(x[["se"]]), " errors (", paste0(x[["clus"]], collapse = ", "), "; ",
      paste(x[["ncl"]], collapse = " and "), ")\n", sep = "")
  print(x[["coef"]], row.names = FALSE, digits = 4)
  if (iv) {
    cat("first stage: robust Wald F, Kleibergen-Paap rk Wald F, effective F:\n")
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

#' Instrumental variables for a count outcome with fixed effects
#'
#' Instruments an exposure in the count model of [ppmlfe()], in one of two
#' ways (`how`). The first stages regress each endogenous column on the
#' excluded instruments and the controls with the fixed effects (as in
#' [ivfe()]). Then:
#'
#' `how = "cf"` (the default), the control function: the first-stage
#' residuals are added to the Poisson pseudo-likelihood fit with the same
#' fixed effects and offset. The residual takes up the part of the exposure
#' that moves with unmeasured causes of the outcome; its own coefficient
#' tests whether the exposure was endogenous at all. With a log link this
#' is consistent when the first-stage residual enters the mean as a linear
#' index (Wooldridge 2014), and it is the default for that reason.
#'
#' `how = "plug"`, the fitted exposure (the exposure less its first-stage
#' residual) in place of the exposure in the Poisson fit: the second stage
#' as Bondy, Roth and Sager (2020, equation 3, p. 565) write it, with
#' "AQI-hat". With a log link the plug-in is not consistent in general
#' (the mean of `exp()` of a fitted value is not `exp()` of the mean); it
#' is offered to reproduce the paper's specification.
#'
#' Standard errors (`se`): `"boot"` (the default), a cluster bootstrap that
#' refits both stages on each of `nboot` draws of the clusters of `clus`,
#' made on the R side with R's generator (so `set.seed()` or `seed`
#' reproduce them), each draw weighting a cluster by how often it was
#' drawn; `"naive"`, the clustered errors of the second-stage Poisson fit
#' alone, which treat the first-stage residual or fitted value as data and
#' are therefore wrong (too small, as a rule): they are given only to
#' reproduce papers that report them. Bondy, Roth and Sager do not say how
#' their second-stage errors account for the first stage. The exponential
#' generalised method of moments with multiplicative errors (Mullahy 1997;
#' Windmeijer and Santos Silva 1997) is not offered: its papers were not at
#' hand to build it from.
#'
#' Separation is handled as in [ppmlfe()] (groups of zeros dropped), in the
#' full sample and in each draw; a draw that cannot be fitted (an
#' instrument without variation in it, say) is left out and counted in
#' `nfail`. The first-stage table is that of [ivfe()], clustered by `clus`.
#'
#' @param data,y,fe,expo,level,tol,maxit As in [ppmlfe()].
#' @param endo The names of the endogenous columns.
#' @param inst The names of the excluded instruments, at least as many.
#' @param x `NULL` or the names of the exogenous controls. For `print()`,
#'   an object of class `ppmlcf`.
#' @param clus `NULL` (the first fixed effect) or the name of one cluster
#'   column: the unit the bootstrap draws and the naive errors cluster by.
#' @param how `"cf"` or `"plug"`.
#' @param se `"boot"` or `"naive"`.
#' @param ssc The small-sample factor of the naive errors, as in
#'   [ppmlfe()].
#' @param nboot The number of bootstrap draws, at least 2.
#' @param seed `NULL` (the generator's current state) or one whole number.
#' @param object An object of class `ppmlcf`.
#' @param ... Not used.
#' @return An object of class `ppmlcf`, a list: `coef` (as in [ppmlfe()],
#'   the endogenous columns, the controls, then for `"cf"` the residual
#'   terms `"vhat:<endo>"`); `first` (as in [ivfe()]); `draws`, the
#'   bootstrap estimates (one row per draw; `NULL` for `"naive"`); `nfail`;
#'   `n`; `nsep`; `ncl`; `nboot`; `how`, `se`; `fe`, `clus` and `level`.
#'   `summary()` returns `coef`.
#' @references
#'   Bondy, M., Roth, S. and Sager, L. (2020). Crime is in the air: the
#'   contemporaneous relationship between air pollution and crime.
#'   *Journal of the Association of Environmental and Resource Economists*,
#'   7(3), 555-585. \doi{10.1086/707127}
#'
#'   Wooldridge, J. M. (2014). Quasi-maximum likelihood estimation and
#'   testing for nonlinear models with endogenous explanatory variables.
#'   *Journal of Econometrics*, 182(1), 226-234.
#'   \doi{10.1016/j.jeconom.2014.04.020}
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
#' # the paper's plug-in, with its naive errors
#' summary(ppmlcf(d, "crimes", endo = "pm", inst = colnames(w),
#'                fe = c("area", "day"), how = "plug", se = "naive"))
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
  how = c("cf", "plug"),
  se = c("boot", "naive"),
  ssc = c("G", "nk", "none"),
  nboot = 199,
  seed = NULL,
  level = 0.95,
  tol = 1e-12,
  maxit = 100
) {
  f <- "ppmlcf"
  how <- match.arg(how)
  se <- match.arg(se)
  ssc <- match.arg(ssc)
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
  tm <- c(endo, x, if (how == "cf") paste0("vhat:", endo))
  if (se == "naive") {
    r1 <- .axivfit(d, y, endo, inst, x, fe, NULL,
                   list(c1 = cl$c1, n1 = cl$n1, c2 = NULL, n2 = 0L), tol, f)
    fr <- r1$fres
    e <- d
    if (how == "cf") {
      for (j in seq_along(endo)) e[[paste0("vhat:", endo[j])]] <- fr[, j]
    } else {
      for (j in seq_along(endo)) e[[endo[j]]] <- e[[endo[j]]] - fr[, j]
    }
    p <- ppmlfe(e, y, tm, fe, expo = expo, clus = clus, ssc = ssc, level = level,
                tol = tol, maxit = maxit, sep = "none")
    return(structure(
      list(coef = p$coef, first = .axfirst(endo, r1$fst, length(inst)), draws = NULL,
           nfail = 0L, n = length(yv), nsep = sum(!k), ncl = cl$g, nboot = 0L, how = how,
           se = se, fe = fe, clus = clus, level = level),
      class = "ppmlcf"
    ))
  }
  nboot <- .cknum(nboot, "nboot", f, min = 2, whole = TRUE)
  if (!is.null(seed)) seed <- .cknum(seed, "seed", f, min = 0, whole = TRUE)
  bm <- .axdraw(cl$n1, nboot, seed)
  storage.mode(bm) <- "integer"
  r <- .Call(C_clab_r_ppmlcf, yv, .axmat(d, endo, f), .axmat(d, x, f),
             .axmat(d, inst, f), off, fx$c1, fx$n1, fx$c2, fx$n2, cl$c1,
             cl$n1, bm, if (how == "plug") 1L else 0L, tol, as.integer(maxit))
  dr <- r$draw
  colnames(dr) <- tm
  nfail <- as.integer(r$st[3L])
  if (nfail > 0L) {
    warning(f, ": ", nfail, " of ", nboot, " bootstrap draws could not be ",
            "fitted and are left out.", call. = FALSE)
  }
  if (nboot - nfail < 2L) .ckstop(f, "fewer than two bootstrap draws could be fitted")
  sev <- apply(dr, 2L, stats::sd, na.rm = TRUE)
  tab <- .axtab(tm, r$beta, sev, level)
  tab$rr <- exp(tab$est)
  structure(
    list(coef = tab, first = .axfirst(endo, r$fst, length(inst)), draws = dr, nfail = nfail,
         n = length(yv), nsep = sum(!k), ncl = cl$g, nboot = nboot, how = how, se = se,
         fe = fe, clus = clus, level = level),
    class = "ppmlcf"
  )
}

#' @rdname ppmlcf
#' @export
print.ppmlcf <- function(x, ...) {
  cat("ppmlcf: Poisson ", if (x[["how"]] == "cf") "control function" else "plug-in",
      ", ", x[["n"]], " rows, fixed effects ", paste(x[["fe"]], collapse = " and "), "; ",
      if (x[["se"]] == "boot") {
        paste0(x[["nboot"]], " bootstrap draws of ", x[["ncl"]], " clusters (", x[["clus"]], ")")
      } else {
        "naive second-stage errors (they ignore the first stage)"
      }, "\n", sep = "")
  print(x[["coef"]], row.names = FALSE, digits = 4)
  cat("first stage: robust Wald F, Kleibergen-Paap rk Wald F, effective F:\n")
  print(x[["first"]], row.names = FALSE, digits = 4)
  invisible(x)
}

#' @rdname ppmlcf
#' @export
summary.ppmlcf <- function(object, ...) object[["coef"]]

# The sector (1..nb) of each direction, by the kernel's edges.
.axsec <- function(dir, width, start, edges) {
  a <- (dir - start) %% 360
  b <- if (edges == "paper") ifelse(a == 0, 0, ceiling(a / width) - 1) else floor(a / width)
  pmin(b, round(360 / width) - 1) + 1
}

#' Wind-direction instruments
#'
#' Turns wind directions into indicators of the sector the wind blows
#' from, one sector left out, optionally one set per area so that the same
#' wind may act differently in different parts of a city; from hourly
#' directions it first makes one per day.
#'
#' The choices, and which paper made which:
#'
#' Sectors (`width`, `start`, `edges`). Bondy, Roth and Sager (2020,
#' footnote 18) use quadrants "0-90, 91-180, 181-270" degrees with the
#' fourth left out, whole degrees closed at the top: `edges = "paper"`.
#' The default `"halfopen"` puts a boundary direction in the sector it
#' opens, `[lo, hi)`, which is the same for directions that are not whole
#' numbers. Bondy, Roth and Sager let each of five regions of London have
#' its own set (`area`), 15 instruments; their appendix also ranks four or
#' eight directions by each ward's long-run pollution (table A3), which is
#' not built here.
#'
#' One sector (`sect`). Herrnstadt, Heyes, Muehlegger and Saberian (2019,
#' footnote 24) use one indicator, wind from between 225 and 285 degrees
#' (60 degrees); their table 5 widens it to 75 and 90 degrees without
#' saying where the wider sectors start, so give `sect` as you choose it
#' (for example `c(217.5, 292.5)` for 75 degrees about the same centre).
#'
#' From hours to days (`by`, `avg`, `share`). With `by`, `dir` holds
#' several directions per day (hours) and each group gives one row. `avg =
#' "vector"` (the default) takes the direction of the mean wind vector,
#' `atan2(mean(sin), mean(cos))`, weighted by `spd` when given: Herrnstadt
#' and Muehlegger (2015, footnote 7) point out that 1 and 359 degrees
#' average to north, not south. `avg = "mode"` takes the sector most hours
#' fall in (ties: the first), the "predominant wind direction" of Bondy,
#' Roth and Sager's figure A1. `avg = "mean"` is the plain mean of the
#' angles, wrong across north and given only for comparison. `share =
#' TRUE` instead gives the share of each group's hours in each sector,
#' Bondy, Roth and Sager's instrument ("the share of hours in the 24 hours
#' of each day", footnote 18); `avg` is then not used.
#'
#' Directions are where the wind blows from, in degrees clockwise from
#' north; 360 is 0. The indicators are instruments only if the wind
#' direction affects the outcome through the exposure alone (see
#' [ivfe()]).
#'
#' @param dir Wind directions in degrees (numbers; `NA` gives a row of
#'   `NA`, and is left out of a group's average).
#' @param width The width of a sector in degrees; it must divide 360.
#' @param start Where the first sector begins, in degrees.
#' @param area `NULL` or the area of each direction (labels; constant
#'   within a `by` group): one set of indicators per area, in sorted order
#'   of the labels.
#' @param base The sector left out, counted from 1.
#' @param sect `NULL`, or two numbers `c(from, to)` in degrees: then one
#'   indicator, 1 when the wind blows from between `from` (included) and
#'   `to` (excluded), going clockwise; `width`, `start`, `base` and `edges`
#'   are not used.
#' @param edges `"halfopen"` or `"paper"`.
#' @param by `NULL`, or a group (day, or place and day) for each direction.
#' @param avg `"vector"`, `"mode"` or `"mean"`: how a group's directions
#'   become one.
#' @param share `TRUE` for the shares of each group's directions in each
#'   sector.
#' @param spd `NULL` or wind speeds, weights of the vector average.
#' @return A numeric matrix, one row per direction (or per `by` group,
#'   named by it), columns named `"w<from>_<to>"`, with `":<area>"` added
#'   when `area` is given. With `by`, the attribute `"dir"` holds the
#'   group's direction (for `share`, its vector average).
#' @references
#'   Bondy, M., Roth, S. and Sager, L. (2020). Crime is in the air: the
#'   contemporaneous relationship between air pollution and crime.
#'   *Journal of the Association of Environmental and Resource Economists*,
#'   7(3), 555-585. \doi{10.1086/707127}
#'
#'   Herrnstadt, E., Heyes, A., Muehlegger, E. and Saberian, S. (2019). Air
#'   pollution as a cause of violent crime. Working paper.
#'
#'   Herrnstadt, E. and Muehlegger, E. (2015). Air pollution and criminal
#'   activity: evidence from Chicago microdata. NBER Working Paper 21787.
#'   \doi{10.3386/w21787}
#' @seealso [ivfe()], [ppmlcf()], [invsn()]
#' @examples
#' windiv(c(10, 100, 200, 300))
#' # the paper's whole-degree quadrants: 90 is in the first
#' windiv(c(90, 91), edges = "paper")
#' # quadrants by region, as in London
#' windiv(c(10, 100, 200, 300), area = c("east", "east", "west", "west"))
#' # one westerly sector
#' windiv(c(230, 300), sect = c(225, 285))
#' # hours to days: 350 and 10 degrees average to north
#' windiv(c(350, 10, 100, 120), by = c(1, 1, 2, 2))
#' windiv(c(350, 10, 100, 120), by = c(1, 1, 2, 2), share = TRUE)
#' @export
windiv <- function(
  dir,
  width = 90,
  start = 0,
  area = NULL,
  base = 1,
  sect = NULL,
  edges = c("halfopen", "paper"),
  by = NULL,
  avg = c("vector", "mode", "mean"),
  share = FALSE,
  spd = NULL
) {
  f <- "windiv"
  edges <- match.arg(edges)
  avg <- match.arg(avg)
  if (!is.numeric(dir) || !length(dir)) .ckstop(f, "`dir` must be numeric directions in degrees")
  dir <- as.numeric(dir)
  if (!is.logical(share) || length(share) != 1L || is.na(share)) {
    .ckstop(f, "`share` must be TRUE or FALSE")
  }
  if (!is.null(spd) && (!is.numeric(spd) || length(spd) != length(dir) || any(spd < 0, na.rm = TRUE))) {
    .ckstop(f, "`spd` must give one speed from 0 per direction")
  }
  if (!is.null(area) && (length(area) != length(dir) || anyNA(area))) {
    .ckstop(f, "`area` must give one label per direction, none missing")
  }
  if (is.null(sect)) {
    width <- .cknum(width, "width", f, positive = TRUE)
    start <- .cknum(start, "start", f)
    nb <- round(360 / width)
    if (width >= 360 || abs(360 / width - nb) > 1e-9) {
      .ckstop(f, "`width` must divide 360 into two or more sectors")
    }
    base <- .cknum(base, "base", f, min = 1, max = nb, whole = TRUE)
  } else if (!is.numeric(sect) || length(sect) != 2L || !all(is.finite(sect)) ||
               sect[1L] %% 360 == sect[2L] %% 360) {
    .ckstop(f, "`sect` must be two different finite directions in degrees")
  }
  gl <- NULL
  if (!is.null(by)) {
    if (length(by) != length(dir) || anyNA(by)) {
      .ckstop(f, "`by` must give one group per direction, none missing")
    }
    gk <- .axcode(by)
    gl <- unique(as.character(by))
    rd <- dir * pi / 180
    wv <- if (is.null(spd)) rep(1, length(dir)) else as.numeric(spd)
    ok <- is.finite(dir) & is.finite(wv)
    vs <- function(fun) as.numeric(tapply(ifelse(ok, fun(rd) * wv, NA), gk$c, mean, na.rm = TRUE))
    vdir <- (atan2(vs(sin), vs(cos)) * 180 / pi) %% 360
    if (!is.null(area)) {
      if (any(tapply(as.character(area), gk$c, function(v) length(unique(v))) > 1L)) {
        .ckstop(f, "`area` must be constant within each `by` group")
      }
      area <- as.character(area)[!duplicated(gk$c)]
    }
    if (share) {
      if (!is.null(sect)) .ckstop(f, "`share` needs sectors by `width`, not `sect`")
      full <- windiv(ifelse(ok, dir, NA), width, start, NULL, base = base, edges = edges)
      lev <- if (is.null(area)) NULL else sort(unique(area))
      m <- rowsum(ifelse(is.na(full), 0, full), gk$c, reorder = TRUE) /
        as.numeric(tabulate(gk$c[ok] + 1L, gk$n))
      if (!is.null(lev)) {
        ai <- match(area, lev)
        mm <- matrix(0, nrow(m), ncol(m) * length(lev))
        for (r in seq_len(nrow(m))) mm[r, (ai[r] - 1L) * ncol(m) + seq_len(ncol(m))] <- m[r, ]
        colnames(mm) <- paste0(rep(colnames(full), length(lev)), ":", rep(lev, each = ncol(m)))
        m <- mm
      }
      rownames(m) <- gl
      attr(m, "dir") <- vdir
      return(m)
    }
    dir <- switch(avg,
      vector = vdir,
      mean = as.numeric(tapply(ifelse(ok, dir, NA), gk$c, mean, na.rm = TRUE)),
      mode = {
        if (!is.null(sect)) .ckstop(f, "`avg = \"mode\"` needs sectors by `width`, not `sect`")
        sc <- .axsec(dir, width, start, edges)
        md <- as.numeric(tapply(ifelse(ok, sc, NA), gk$c, function(v) {
          v <- v[!is.na(v)]
          if (!length(v)) NA else as.numeric(names(which.max(table(factor(v, levels = seq_len(nb))))))
        }))
        # the middle of the sector, which falls in it under either edges
        (start + (md - 0.5) * width) %% 360
      })
  }
  if (!is.null(sect)) {
    m <- matrix(as.numeric((dir - sect[1L]) %% 360 < (sect[2L] - sect[1L]) %% 360))
    colnames(m) <- paste0("w", format(sect[1L]), "_", format(sect[2L]))
  } else {
    ac <- NULL
    na <- 1L
    lev <- NULL
    if (!is.null(area)) {
      lev <- if (is.factor(area)) levels(droplevels(area)) else sort(unique(as.character(area)))
      ac <- match(as.character(area), lev) - 1L
      na <- length(lev)
    }
    m <- .Call(C_clab_r_windbin, dir, width, start, ac, na, as.integer(base - 1),
               if (edges == "paper") 1L else 0L)
    lo <- (start + (seq_len(nb) - 1) * width) %% 360
    nm <- paste0("w", format(lo, trim = TRUE), "_",
                 format((lo + width) %% 360, trim = TRUE))[-base]
    colnames(m) <- if (is.null(lev)) nm else paste0(rep(nm, na), ":", rep(lev, each = nb - 1L))
  }
  if (!is.null(gl)) {
    rownames(m) <- gl
    attr(m, "dir") <- dir
  }
  m
}

#' Temperature inversions
#'
#' The instrument of Bondy, Roth and Sager (2020, p. 566-567): "the
#' continuous difference between temperatures at the ground level and
#' temperatures at circa 750 m altitude (925 hectopascals [hPa] pressure
#' level)". On most days the air cools with height; in an inversion it is
#' warmer above than at the ground, and pollutants are trapped near the
#' ground.
#'
#' The choices. The paper does not say which temperature is subtracted
#' from which: `sign = "upper"` (the default) gives the upper temperature
#' less the ground one, positive in an inversion; `sign = "ground"` the
#' reverse, as the paper's words run. `inv` is 1 in an inversion (upper
#' warmer than ground) whatever the sign; the indicator, rather than the
#' continuous difference, is the other common form of the instrument. With
#' `brk`, `bin` gives the strength in bins of the difference (by the
#' chosen sign), for an instrument that lets weak and strong inversions act
#' differently.
#'
#' Cautions. Inversions cover large areas (the paper's satellite readings
#' are on a grid of about 110 by 70 km), so they vary little between places
#' in a city; Bondy, Roth and Sager use them at the borough and city level
#' and report a weaker first stage there than for wind (their table 5). As
#' an instrument the inversion must affect the outcome only through the
#' exposure; inversions come with cold, still, clear weather, which may
#' change what people do as well, so weather controls matter (see
#' [ivfe()]).
#'
#' @param tlow Temperatures at the ground (numbers, `NA` allowed).
#' @param tup Temperatures aloft, at the same places and times, in the same
#'   unit.
#' @param sign `"upper"` or `"ground"`: which temperature comes first.
#' @param brk `NULL` or the breaks of the bins of the difference (as in
#'   [cut()], right-closed).
#' @return A data frame: `diff`, the difference by `sign`; `inv`, 1 when
#'   `tup` is above `tlow` and 0 otherwise (`NA` where a temperature is
#'   missing); and with `brk`, `bin`, a factor.
#' @references
#'   Bondy, M., Roth, S. and Sager, L. (2020). Crime is in the air: the
#'   contemporaneous relationship between air pollution and crime.
#'   *Journal of the Association of Environmental and Resource Economists*,
#'   7(3), 555-585. \doi{10.1086/707127}
#' @seealso [windiv()], [ivfe()], [ppmlcf()]
#' @examples
#' invsn(c(10, 5, 2), c(8, 6, 2))
#' invsn(c(10, 5, 2), c(8, 6, 2), sign = "ground", brk = c(-Inf, -1, 0, Inf))
#' @export
invsn <- function(tlow, tup, sign = c("upper", "ground"), brk = NULL) {
  f <- "invsn"
  sign <- match.arg(sign)
  if (!is.numeric(tlow) || !is.numeric(tup) || length(tlow) != length(tup) ||
        !length(tlow)) {
    .ckstop(f, "`tlow` and `tup` must be numeric vectors of one length")
  }
  up <- as.numeric(tup) - as.numeric(tlow)
  r <- data.frame(diff = if (sign == "upper") up else -up, inv = as.integer(up > 0))
  if (!is.null(brk)) {
    if (!is.numeric(brk) || length(brk) < 2L || anyNA(brk) || is.unsorted(brk, strictly = TRUE)) {
      .ckstop(f, "`brk` must be two or more increasing numbers")
    }
    r$bin <- cut(r$diff, brk)
  }
  r
}

#' Irrelevant-day placebo
#'
#' Refits an estimator with the exposure taken from other days of the same
#' place, leads and lags in one call, and returns one table of
#' coefficients, the same day first. This is the placebo of Bondy, Roth
#' and Sager (2020, table A1 and footnote 27: 7, 31 and 365 days back and
#' forward) and of Herrnstadt, Heyes, Muehlegger and Saberian (2019, table
#' 9: 100 days); `shifts = c(-365, -100, -31, -7, 7, 31, 100, 365)` runs
#' both.
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
#' The shift is by the day's value: row `(place, day)` gets the value of
#' `(place, day + shift)` and is dropped when that row does not exist;
#' `n` reports the rows left. For the instrumental-variable estimators the
#' papers do not say what is moved: `shift = "both"` (the default) moves
#' the exposure `var` and the columns in `with` (the instruments) together,
#' a placebo of the whole design on another day; `"x"` moves the exposure
#' alone and keeps the same-day instruments; `"z"` moves the instruments
#' alone (a reduced-form placebo).
#'
#' @param data A data frame, one row per place and day.
#' @param var The name of the exposure column.
#' @param unit,day The names of the place and day columns; days are
#'   numbers or dates (one per unit at most).
#' @param shifts The leads and lags in days, whole numbers other than 0;
#'   negative is the past.
#' @param est The estimator: [ppmlfe()], [ivfe()] or [ppmlcf()], or
#'   another function with a `data` argument returning a list with `coef`
#'   holding a row whose `term` is `var`.
#' @param with `NULL` or the names of other columns to move (the
#'   instruments).
#' @param shift `"both"`, `"x"` or `"z"`: what is moved.
#' @param ... The other arguments of `est` (`y`, `x`, `fe`, ...).
#' @param object,x An object of class `plcday` (for the methods).
#' @return An object of class `plcday`, a list: `tab`, a data frame with
#'   one row per shift, 0 (the same day) first (`shift`, `est`, `se`,
#'   `lo`, `hi`, `p`, `n`); `var`, `shift`. `summary()` returns `tab`.
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
  shift = c("both", "x", "z"),
  ...
) {
  f <- "plcday"
  shift <- match.arg(shift)
  d <- .axdat(data, f)
  .axcols(d, var, "var", f, n = 1L)
  .axcols(d, unit, "unit", f, n = 1L)
  .axcols(d, day, "day", f, n = 1L)
  .axcols(d, with, "with", f, null = TRUE)
  if (shift == "z" && !length(with)) .ckstop(f, "`shift = \"z\"` needs the instruments in `with`")
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
  mv <- switch(shift, both = c(var, with), x = var, z = with)
  rows <- lapply(c(0, shifts), function(s) {
    ix <- .Call(C_clab_r_dshift, uc, tv, s)
    k <- !is.na(ix)
    if (sum(k) < 2L) .ckstop(f, "a shift of %s days leaves fewer than two rows", format(s))
    ds <- d[k, , drop = FALSE]
    for (v in mv) ds[[v]] <- d[[v]][ix[k]]
    r <- est(data = ds, ...)
    cf <- r[["coef"]]
    if (!is.data.frame(cf) || !var %in% cf[["term"]]) {
      .ckstop(f, "`est` returned no coefficient for `%s`", var)
    }
    cf <- cf[cf[["term"]] == var, , drop = FALSE]
    data.frame(shift = s, est = cf[["est"]], se = cf[["se"]], lo = cf[["lo"]],
               hi = cf[["hi"]], p = cf[["p"]], n = sum(k))
  })
  structure(list(tab = do.call(rbind, rows), var = var, shift = shift), class = "plcday")
}

#' @rdname plcday
#' @export
print.plcday <- function(x, ...) {
  cat("plcday: the coefficient of ", x[["var"]],
      " with values of other days (shift 0 is the same day; moved: ",
      x[["shift"]], ")\n", sep = "")
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
#' area belongs to a source within `band` of it. For a segment, the side of
#' the area is set by which side of the line it lies on; the area is
#' downwind on a day when the direction the wind blows to is within `tol`
#' degrees of the line orthogonal to the road pointing to its side, upwind
#' when within `tol` of the opposite direction, and the day is dropped for
#' that road otherwise (wind along the road). For a point source the same
#' rule applies to the direction from the source to the area: this extends
#' the design and is not in the papers.
#'
#' The choices, with the paper's own first:
#'
#' * `tol`: 60 degrees in the paper's main results (p. 12); 30 to 90 in its
#'   robustness table (pp. 15-16).
#' * `band`: one mile (1609 m) in the paper (footnote 23); a quarter mile
#'   to one mile in its robustness table. `dist` (breaks of distance) gives
#'   one downwind coefficient per ring instead, the paper's finding that
#'   the effect declines with distance.
#' * `twin`, areas within `band` of two sources: `"drop"` as the paper
#'   (p. 11-12); `"near"` gives them to the nearer source; `"keep"` counts
#'   them for each source.
#' * `norm`, the outcome: `"mean"`, counts divided by their mean over the
#'   rows used, least squares (the paper's main tables: the coefficient is
#'   a share of the mean count); `"none"`, raw counts, least squares;
#'   `"log"`, `log(1 + y)`, least squares; `"pois"`, the Poisson
#'   pseudo-likelihood of [ppmlfe()] (the paper's appendix table A.1).
#' * `wx`: weather columns interacted with the side effects, the paper's
#'   equation (4), p. 14, which lets each side respond to the weather in
#'   its own way (one side's interactions are left out: the route-date
#'   effects absorb the weather itself).
#' * `se`: `"hc1"` (the default), robust errors as the paper reports;
#'   `"cluster"` (by `clus`, such as side-month or side-year as in its
#'   appendix table A.2), `"twoway"`, `"nw"` (Newey-West within each side
#'   over days, its appendix too) or `"boot"`, as in [ppmlfe()]. Stata's
#'   robust errors use `n / (n - K)`; here that is `ssc = "nk"`.
#'
#' The estimate is the paper's equation (3), p. 12: `y_ist = a_is + g_it +
#' b downwind_ist + e_ist`, with route-side fixed effects `a_is` and
#' route-date fixed effects `g_it`; with `agg = TRUE` the counts of the
#' areas on one side of one segment are summed by day first, as the paper
#' counts crimes by interstate side and day. Without aggregation (and for
#' point sources) the first effect is the area.
#'
#' Placebos. Run the design on an outcome the exposure should not affect
#' (the paper uses property crime); and, with `offs`, move every segment
#' sideways by each offset and re-estimate, the paper's placebo
#' interstates drawn parallel to I-290 a mile apart (pp. 17-18): the effect
#' should peak at offset 0.
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
#' @param dist `NULL` or increasing breaks of distance from 0 for rings.
#' @param twin `"drop"`, `"near"` or `"keep"`.
#' @param agg `TRUE` to sum the areas on each side of each segment by day.
#' @param norm `"mean"`, `"none"`, `"log"` or `"pois"`.
#' @param wx `NULL` or the names of weather columns, constant within a side
#'   and day.
#' @param se,ssc,lag,nboot,seed As in [ppmlfe()]; the units of `"nw"` are
#'   the sides (areas) and the time is `day`.
#' @param clus `NULL` (for `"cluster"`, the sides) or the names of cluster
#'   columns, constant within each side (area) and day.
#' @param offs `NULL` or offsets (in the coordinates' unit, positive to the
#'   left of each segment) for the placebo segments; segments only.
#' @param level The confidence level.
#' @param object,x An object of class `updwn` (for the methods).
#' @param ... Not used.
#' @return An object of class `updwn`, a list: `coef`, the rows of the
#'   downwind terms (one, or one per ring); `class`, the rows by outcome of
#'   the classification (`downwind`, `upwind`, `along` (wind along the road
#'   or crosswind), `nosource`, `twosources`, `onsource`, `missing`); `n`,
#'   the rows of the fit; `fit`, the fit; `placebo`, `NULL` or one row per
#'   offset (`off`, `est`, `se`, `lo`, `hi`, `p`, `n`); `norm`, `tol`,
#'   `band`. `summary()` returns `coef`, with the placebo rows when there
#'   are any.
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
#' d$wind <- runif(120, 0, 360)[d$day]
#' blowsn <- cos((d$wind + 180) * pi / 180)  # +1 when the wind blows north
#' dn <- (d$ay > 0 & blowsn > 0.5) | (d$ay < 0 & blowsn < -0.5)
#' d$crimes <- rpois(nrow(d), 3 * (1 + 0.1 * dn))
#' road <- data.frame(x0 = 0, y0 = 0, x1 = 5000, y1 = 0)
#' updwn(d, "crimes", "ax", "ay", "day", "wind", road, band = 1600,
#'       offs = c(-1600, 0, 1600))
#' # the Poisson model and two distance rings
#' summary(updwn(d, "crimes", "ax", "ay", "day", "wind", road, band = 1600,
#'               norm = "pois", dist = c(0, 500, 1600)))
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
  dist = NULL,
  twin = c("drop", "near", "keep"),
  agg = TRUE,
  norm = c("mean", "none", "log", "pois"),
  wx = NULL,
  se = c("hc1", "cluster", "twoway", "nw", "boot"),
  ssc = c("G", "nk", "none"),
  clus = NULL,
  lag = NULL,
  nboot = 199,
  seed = NULL,
  offs = NULL,
  level = 0.95
) {
  f <- "updwn"
  twin <- match.arg(twin)
  norm <- match.arg(norm)
  se <- match.arg(se)
  ssc <- match.arg(ssc)
  d <- .axdat(data, f)
  .axcols(d, y, "y", f, n = 1L)
  .axcols(d, c(ax, ay), "ax, ay", f, n = 2L)
  .axcols(d, day, "day", f, n = 1L)
  .axcols(d, wind, "wind", f, n = 1L)
  .axcols(d, area, "area", f, n = 1L, null = TRUE)
  .axcols(d, clus, "clus", f, n = 1:2, null = TRUE)
  .axcols(d, wx, "wx", f, null = TRUE)
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
  if (!is.null(dist)) {
    if (!is.numeric(dist) || length(dist) < 2L || anyNA(dist) || dist[1L] < 0 ||
          is.unsorted(dist, strictly = TRUE)) {
      .ckstop(f, "`dist` must be two or more increasing distances from 0")
    }
  }
  if (!is.logical(agg) || length(agg) != 1L || is.na(agg)) .ckstop(f, "`agg` must be TRUE or FALSE")
  if (is.null(area) && (!agg || any(!line))) {
    .ckstop(f, "`area` is needed for point sources and for agg = FALSE")
  }
  if (!is.null(offs)) {
    if (any(!line)) .ckstop(f, "`offs` moves segments: every source must be a segment")
    if (!is.numeric(offs) || !length(offs) || !all(is.finite(offs))) {
      .ckstop(f, "`offs` must be finite numbers")
    }
  }
  if (se == "twoway" && length(clus) != 2L) .ckstop(f, "`se = \"twoway\"` needs two `clus` columns")
  level <- .axlevel(level, f)
  yv <- .axmat(d, y, f)[, 1L]
  if (any(yv < 0)) .ckstop(f, "column `%s` must hold counts, at least 0", y)
  xy <- .axmat(d, c(ax, ay), f)
  wd <- .axmat(d, wind, f)[, 1L]
  if (length(wx)) .axmat(d, wx, f)
  dv <- as.character(d[[day]])
  tnum <- if (is.numeric(d[[day]]) || inherits(d[[day]], "Date")) as.numeric(d[[day]])
  if (se == "nw" && is.null(tnum)) .ckstop(f, "`se = \"nw\"` needs a day column of numbers or dates")
  av <- if (!is.null(area)) as.character(d[[area]])
  if (!is.null(area) && anyDuplicated(paste(av, dv, sep = "\r"))) {
    .ckstop(f, "an area has two rows on one day")
  }
  labs <- c(`1` = "downwind", `0` = "upwind", `-1` = "along", `-2` = "nosource",
            `-3` = "twosources", `-4` = "onsource", `-5` = "missing")
  same <- function(v, cell) any(tapply(v, cell, function(u) length(unique(u))) > 1L)
  run <- function(a0, b0, a1, b1) {
    if (twin == "keep") {
      parts <- lapply(seq_along(a0), function(j) {
        k <- .Call(C_clab_r_updwn, xy[, 1L], xy[, 2L], wd, a0[j], b0[j], a1[j], b1[j], band,
                   tol, 0L)
        k$src[!is.na(k$src)] <- j
        k$row <- seq_len(nrow(xy))
        k
      })
      k <- list(src = unlist(lapply(parts, `[[`, "src")), side = unlist(lapply(parts, `[[`, "side")),
                dn = unlist(lapply(parts, `[[`, "dn")), dst = unlist(lapply(parts, `[[`, "dst")),
                row = unlist(lapply(parts, `[[`, "row")))
      # an area near no source counts once in the table
      hit <- k$dn != -2L
      none <- setdiff(seq_len(nrow(xy)), k$row[hit])
      tb <- table(factor(labs[as.character(c(k$dn[hit], rep(-2L, length(none))))], levels = labs))
    } else {
      k <- .Call(C_clab_r_updwn, xy[, 1L], xy[, 2L], wd, a0, b0, a1, b1, band, tol,
                 if (twin == "near") 1L else 0L)
      k$row <- seq_len(nrow(xy))
      tb <- table(factor(labs[as.character(k$dn)], levels = labs))
    }
    use <- k$dn >= 0L
    if (sum(use) < 2L) .ckstop(f, "fewer than two rows are upwind or downwind of a source")
    rw <- k$row[use]
    ring <- if (is.null(dist)) rep("", sum(use)) else {
      rg <- cut(k$dst[use], dist, include.lowest = TRUE)
      as.character(rg)
    }
    keep <- !is.na(ring)
    rw <- rw[keep]
    ring <- ring[keep]
    sr <- k$src[use][keep]
    sd0 <- k$side[use][keep]
    dn <- k$dn[use][keep]
    if (!length(rw)) .ckstop(f, "no row falls in the rings of `dist`")
    ln <- line[sr]
    g <- ifelse(ln & agg, paste("s", sr, sd0, ring, sep = "\r"),
                paste("a", if (is.null(area)) "" else av[rw], sr, ring, sep = "\r"))
    cell <- paste(g, dv[rw], sep = "\r")
    sdy <- paste(sr, dv[rw], sep = "\r")
    if (same(dn, cell)) {
      .ckstop(f, "the wind differs between areas of one side and day: give one direction per day, or agg = FALSE")
    }
    for (cc in c(clus, wx)) {
      if (same(as.character(d[[cc]])[rw], cell)) {
        .ckstop(f, "column `%s` must be constant within each side (area) and day", cc)
      }
    }
    o <- !duplicated(cell)
    e <- data.frame(y = as.numeric(rowsum(yv[rw], cell, reorder = FALSE)),
                    g = g[o], sd = sdy[o], down = as.numeric(dn[o]), row = seq_len(sum(o)))
    if (!is.null(tnum)) e$t <- tnum[rw][o]
    for (cc in clus) e[[cc]] <- d[[cc]][rw][o]
    xs <- "down"
    if (!is.null(dist)) {
      rl <- ring[o]
      lv <- levels(cut(dist[-1L], dist, include.lowest = TRUE))
      xs <- character()
      for (l in lv[lv %in% rl]) {
        nm <- paste0("down:", l)
        e[[nm]] <- e$down * (rl == l)
        xs <- c(xs, nm)
      }
    }
    if (length(wx)) {
      # each side's own response to the weather, less one side per source
      # (the route-date effects absorb the weather of each route)
      gl <- unique(e$g)
      gs <- sr[o][match(gl, e$g)]
      for (w in wx) {
        wv <- as.numeric(d[[w]][rw][o])
        for (l in gl[duplicated(gs)]) {
          nm <- paste0(w, ":side", match(l, gl))
          e[[nm]] <- wv * (e$g == l)
        }
      }
    }
    xa <- setdiff(names(e), c("y", "g", "sd", "down", "row", "t", clus))
    xa <- c(xs, setdiff(xa, xs))
    cl <- if (se == "hc1") "row" else if (is.null(clus)) "g" else clus
    tcol <- if (se == "nw") "t" else NULL
    if (norm == "pois") {
      r <- ppmlfe(e, "y", xa, fe = c("g", "sd"), clus = cl, se = se, ssc = ssc, time = tcol,
                  lag = lag, nboot = nboot, seed = seed, level = level)
    } else {
      e$y <- switch(norm, mean = {
        m <- mean(e$y)
        if (!(m > 0)) .ckstop(f, "the counts are all 0")
        e$y / m
      }, none = e$y, log = log1p(e$y))
      r <- ivfe(e, "y", x = xa, fe = c("g", "sd"), clus = cl, se = se, ssc = ssc, time = tcol,
                lag = lag, nboot = nboot, seed = seed, level = level)
    }
    list(fit = r, tab = tb, n = nrow(e), coef = r$coef[r$coef$term %in% xs, , drop = FALSE])
  }
  r0 <- run(s0, t0, s1, t1)
  pl <- NULL
  if (!is.null(offs)) {
    len <- sqrt((s1 - s0)^2 + (t1 - t0)^2)
    nx <- -(t1 - t0) / len
    ny <- (s1 - s0) / len
    pl <- do.call(rbind, lapply(offs, function(o) {
      ro <- run(s0 + o * nx, t0 + o * ny, s1 + o * nx, t1 + o * ny)
      c0 <- ro$coef[1L, ]
      data.frame(off = o, est = c0$est, se = c0$se, lo = c0$lo, hi = c0$hi,
                 p = c0$p, n = ro$n)
    }))
  }
  structure(
    list(coef = r0$coef, class = r0$tab, n = r0$n, fit = r0$fit, placebo = pl,
         norm = norm, tol = tol, band = band),
    class = "updwn"
  )
}

#' @rdname updwn
#' @export
print.updwn <- function(x, ...) {
  cat("updwn: downwind against upwind, ", x[["n"]], " side (area) days, outcome ",
      c(mean = "over its mean", none = "as counts", log = "log(1 + y)",
        pois = "by Poisson")[[x[["norm"]]]],
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
  cbind(off = c(rep(NA, nrow(cf)), object[["placebo"]][["off"]]),
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
