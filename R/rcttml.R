#' Targeted maximum likelihood for a randomised trial
#'
#' The average treatment effect of a randomised trial with adjustment for
#' baseline covariates by targeted maximum likelihood (Moore and van der
#' Laan 2009): more precise than the difference in means when the
#' covariates predict the outcome, and consistent whatever the outcome
#' model, since the chance of treatment is known by design.
#'
#' The outcome regression `Q(d, x)` is fitted on an intercept, `d` and `x`
#' (and `d` times each covariate with `inter = TRUE`): least squares for
#' `fam = "gaussian"`, logistic regression for `"binomial"` (an outcome of
#' 0 and 1). The chance of treatment `g` is the trial's share treated
#' (`g = NULL`), a known value (`g`), or a logistic regression on `x`
#' (`gx = TRUE`), bounded to `[gbound, 1 - gbound]`. `Q` is then
#' fluctuated along the clever covariate `H = d / g - (1 - d) / (1 - g)`
#' by a logistic regression with offset `logit(Q)`; a continuous outcome
#' is first put on `[0, 1]` by `bounds` (its range by default) and `Q`
#' kept in `[qbound, 1 - qbound]` (Gruber and van der Laan 2010). The
#' effect is the mean of `Q*(1, x) - Q*(0, x)`, its standard error that of
#' the influence curve, `sd(IC) / sqrt(n)`. The table also gives the
#' plug-in of the unfluctuated `Q` (`gcomp`, without a standard error) and
#' the difference in means (`unadj`, with the unpooled standard error);
#' intervals and p-values are normal.
#'
#' @param data A data frame, one row per participant of the trial.
#' @param y The name of the outcome column (numbers; 0 or 1 for
#'   `"binomial"`).
#' @param d The name of the treatment column: 0 or 1 (or `FALSE`, `TRUE`),
#'   each arm with two participants or more.
#' @param x `NULL` or the names of the baseline covariate columns
#'   (numbers); for `print()`, an object of class `rcttml`.
#' @param fam `"gaussian"` or `"binomial"`.
#' @param inter `TRUE` to let the effect vary with the covariates in `Q`.
#' @param g `NULL` (the share treated) or the known chance of treatment:
#'   one number strictly between 0 and 1.
#' @param gx `TRUE` to fit the chance of treatment on `x`.
#' @param gbound The bound of `g`, from 0 to below 0.5.
#' @param bounds `NULL` (the range of `y`) or two numbers, below and above
#'   every outcome, that put a continuous outcome on `[0, 1]`.
#' @param qbound The bound of the scaled `Q`, above 0 and below 0.5.
#' @param level The level of the intervals: one number strictly between 0
#'   and 1.
#' @return An object of class `rcttml`, a list: `tab`, a data frame with
#'   rows `tmle`, `gcomp` and `unadj` (`method`, `est`, `se`, `lo`, `hi`,
#'   `p`); `eps`, the fluctuation; `g`, the smallest and largest chance of
#'   treatment used; and `fam`, `bounds`, `n`, `level`.
#' @references
#'   Moore, K. L. and van der Laan, M. J. (2009). Covariate adjustment in
#'   randomized trials with binary outcomes: targeted maximum likelihood
#'   estimation. *Statistics in Medicine*, 28(1), 39-64.
#'   \doi{10.1002/sim.3445}
#'
#'   Gruber, S. and van der Laan, M. J. (2010). A targeted maximum
#'   likelihood estimator of a causal effect on a bounded continuous
#'   outcome. *The International Journal of Biostatistics*, 6(1), 26.
#'   \doi{10.2202/1557-4679.1260}
#'
#'   van der Laan, M. J. and Rubin, D. (2006). Targeted maximum likelihood
#'   learning. *The International Journal of Biostatistics*, 2(1), 11.
#'   \doi{10.2202/1557-4679.1043}
#' @examples
#' set.seed(10)
#' n <- 300
#' d <- data.frame(prior = rpois(n, 2), age = runif(n, 18, 50),
#'                 t = rbinom(n, 1, 0.5))
#' d$days <- 30 - 6 * d$t + 3 * d$prior - 0.1 * d$age + rnorm(n, sd = 4)
#' rcttml(d, "days", "t", c("prior", "age"))
#' d$rearrest <- rbinom(n, 1, plogis(-1 - 0.6 * d$t + 0.4 * d$prior))
#' rcttml(d, "rearrest", "t", "prior", fam = "binomial")$tab
#' @export
rcttml <- function(
  data,
  y,
  d,
  x = NULL,
  fam = c("gaussian", "binomial"),
  inter = FALSE,
  g = NULL,
  gx = FALSE,
  gbound = 0.01,
  bounds = NULL,
  qbound = 0.005,
  level = 0.95
) {
  f <- "rcttml"
  dd <- .axdat(data, f)
  .axcols(dd, y, "y", f, n = 1L)
  .axcols(dd, d, "d", f, n = 1L)
  .axcols(dd, x, "x", f, null = TRUE)
  if (anyDuplicated(c(y, d, x))) {
    .ckstop(f, "`y`, `d` and `x` must name different columns")
  }
  fam <- .ckmarg("fam", f)
  yv <- .axmat(dd, y, f)[, 1L]
  dv <- .axmat(dd, d, f)[, 1L]
  xm <- .axmat(dd, x, f)
  if (!all(dv %in% c(0, 1))) .ckstop(f, "column `%s` must hold 0 or 1", d)
  if (sum(dv) < 2 || sum(1 - dv) < 2) {
    .ckstop(f, "column `%s` must have two participants or more in each arm",
            d)
  }
  if (fam == "binomial" && !all(yv %in% c(0, 1))) {
    .ckstop(f, "with `fam = \"binomial\"`, column `%s` must hold 0 or 1", y)
  }
  inter <- .ckflag(inter, "inter", f)
  gx <- .ckflag(gx, "gx", f)
  if (!is.null(g)) {
    g <- .xpopen(g, "g", f)
    if (gx) .ckstop(f, "give `g` or `gx`, not both")
  }
  if (gx && is.null(x)) .ckstop(f, "`gx` needs covariates `x`")
  gbound <- .cknum(gbound, "gbound", f, min = 0)
  if (gbound >= 0.5) .ckstop(f, "`gbound` must be below 0.5")
  if (fam == "gaussian") {
    if (is.null(bounds)) bounds <- range(yv)
    bounds <- .cknum(bounds, "bounds", f, scalar = FALSE)
    if (length(bounds) != 2L || bounds[1L] >= bounds[2L] ||
          bounds[1L] > min(yv) || bounds[2L] < max(yv)) {
      .ckstop(f, "`bounds` must be two increasing numbers around `y`")
    }
  } else {
    bounds <- c(0, 1)
  }
  qbound <- .cknum(qbound, "qbound", f, positive = TRUE)
  if (qbound >= 0.5) .ckstop(f, "`qbound` must be below 0.5")
  level <- .xpopen(level, "level", f)
  k <- .Call(C_clab_r_rcttml, xm, yv, as.integer(dv), fam == "binomial",
             inter, gx, if (is.null(g)) NA_real_ else g, gbound, bounds,
             qbound)
  est <- k$out[1:3]
  se <- k$out[4:6]
  se[is.nan(se)] <- NA_real_
  z <- stats::qnorm(0.5 + level / 2)
  structure(
    list(tab = data.frame(method = c("tmle", "gcomp", "unadj"), est = est,
                          se = se, lo = est - z * se, hi = est + z * se,
                          p = 2 * stats::pnorm(-abs(est / se)),
                          stringsAsFactors = FALSE),
         eps = k$aux[1L], g = k$aux[2:3], fam = fam, bounds = bounds,
         n = nrow(dd), level = level),
    class = "rcttml"
  )
}

#' @rdname rcttml
#' @param ... Not used by `print()`.
#' @export
print.rcttml <- function(x, ...) {
  cat("rcttml: ", x[["fam"]], " outcome, ", x[["n"]],
      " participants; targeted maximum likelihood, g-computation and the ",
      "difference in means\n", sep = "")
  print(x[["tab"]], row.names = FALSE, digits = 4)
  invisible(x)
}
