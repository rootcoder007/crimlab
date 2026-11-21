# Intervals for a proportion estimated from a complex survey, and the
# release quality of an estimate. The command line and the other
# languages compute the same in C++ (cli/jval.inc).

#' Confidence interval of a survey proportion
#'
#' An interval for a proportion estimated from a complex survey, from the
#' estimate and its design-based standard error (those of [svymx()], say).
#' The usual interval, the estimate plus and minus a multiple of its
#' standard error (Wald), is too narrow for a small sample or a proportion
#' near 0 or 1, and can leave the range 0 to 1; statistical agencies ask
#' for one of the three others here instead.
#'
#' The first two put the design into an effective sample size, the size
#' of a simple random sample with the same variance,
#' `ne = p (1 - p) / se^2`, corrected for the degrees of freedom of the
#' design (Korn and Graubard 1998):
#' `ne* = ne (t(n - 1) / t(df))^2`, with `t()` the quantile of Student's t
#' at the level. With `trunc`, the default, `ne*` is never above `n`, as
#' Korn and Graubard and SAS have it; `survey::svyciprop()` does not
#' truncate, and `trunc = FALSE` gives its interval. An estimate of
#' exactly 0 or 1, or a standard error of 0, has no effective size of its
#' own and takes `n`.
#'
#' `"beta"`: the Clopper-Pearson interval at `ne*` with `p ne*` events, by
#' beta quantiles (Korn and Graubard 1998): the "modified Clopper-Pearson"
#' interval. `"wilson"`: the Wilson score interval at `ne*` with the
#' normal quantile: the "modified Wilson" interval (Dean and Pagano 2015).
#' `"logit"`: the Wald interval of `log(p / (1 - p))`, with standard error
#' `se / (p (1 - p))` and Student's t on `df`, taken back to a proportion;
#' it does not exist at 0 or 1. `"wald"`: `p` plus and minus the t
#' quantile times `se`, the interval [svymx()] gives.
#'
#' Health Canada's guides to the Canadian surveys of alcohol and drug use
#' recommend, in this order, the modified Wilson, the modified
#' Clopper-Pearson and the logit interval, and advise against the Wald
#' interval (Health Canada 2022, section 7.3).
#'
#' Errors: `est` outside 0 to 1; `se` negative; `n` not a whole number of
#' at least 2; `df` below 1; `level` outside (0, 1); `method` not one of
#' the four; `trunc` not `TRUE` or `FALSE`; arguments that are not numeric
#' or of different lengths (a length of 1 is repeated).
#'
#' @param est The estimated proportion, from 0 to 1.
#' @param se Its standard error.
#' @param n The number of respondents of the sample (unweighted).
#' @param df The degrees of freedom of the design: first-stage units less
#'   strata, as [svymx()] reports; `n - 1` by default (a simple random
#'   sample).
#' @param method `"beta"` (the default), `"wilson"`, `"logit"` or
#'   `"wald"`.
#' @param trunc `TRUE` (the default) keeps the effective sample size at
#'   or below `n`.
#' @inheritParams ppmlfe
#' @return A data frame with a row per estimate: `est`, `lo`, `hi` and
#'   `neff`, the effective sample size `ne*` (for every method: it is what
#'   the design is worth). `NA` where an input is missing, and for the
#'   logit interval at 0 or 1.
#' @references
#'   Korn, E. L. and Graubard, B. I. (1998). Confidence intervals for
#'   proportions with small expected number of positive counts estimated
#'   from survey data. *Survey Methodology*, 24(2), 193-201.
#'
#'   Dean, N. and Pagano, M. (2015). Evaluating confidence interval
#'   methods for binomial proportions in clustered surveys. *Journal of
#'   Survey Statistics and Methodology*, 3(4), 484-503.
#'   \doi{10.1093/jssam/smv024}
#'
#'   Health Canada (2022). *Canadian Postsecondary Education Alcohol and
#'   Drug Use Survey, 2021-22: Public Use Microdata File user guide*.
#'   Section 7.3, release guidelines.
#' @seealso [svymx()] for the estimate and its standard error;
#'   [svyqual()] for the quality level of an estimate.
#' @examples
#' # 12% of 400 respondents, standard error 2.1 points, 38 design df
#' svypci(0.12, 0.021, 400, df = 38)
#' for (m in c("beta", "wilson", "logit", "wald")) {
#'   print(svypci(0.12, 0.021, 400, df = 38, method = m))
#' }
#' # near 0 the Wald interval leaves the range; the others do not
#' svypci(0.01, 0.008, 400, df = 38, method = "wald")
#' svypci(0.01, 0.008, 400, df = 38, method = "wilson")
#' try(svypci(1.2, 0.1, 400))
#' @export
svypci <- function(
  est,
  se,
  n,
  df = n - 1,
  method = c("beta", "wilson", "logit", "wald"),
  level = 0.95,
  trunc = TRUE
) {
  fn <- "svypci"
  method <- .ckmarg("method", fn)
  trunc <- .ckflag(trunc, "trunc", fn)
  v <- .svyck(est, se, n, df, level, fn)
  est <- v[["est"]]
  se <- v[["se"]]
  n <- v[["n"]]
  a <- 1 - level
  tq <- stats::qt(1 - a / 2, v[["df"]])
  own <- !is.na(est) & !is.na(se) & est > 0 & est < 1 & se > 0
  ne <- est * (1 - est) / se^2 * (stats::qt(1 - a / 2, n - 1) / tq)^2
  ne <- ifelse(own, if (trunc) pmin(n, ne) else ne, n)
  ne[is.na(est) | is.na(se) | is.na(n) | is.na(tq)] <- NA_real_
  ci <- .svyint(method, est, se, ne, tq, a)
  data.frame(est = est, lo = ci[[1L]], hi = ci[[2L]], neff = ne)
}

#' The arguments of svypci() checked, each of one length.
#' @noRd
.svyck <- function(est, se, n, df, level, fn) {
  num <- function(v, what) {
    if (!is.numeric(v) || !length(v)) {
      .ckstop(fn, "`%s` must be numeric", what)
    }
    as.numeric(v)
  }
  est <- num(est, "est")
  se <- num(se, "se")
  n <- num(n, "n")
  df <- num(df, "df")
  if (any(est < 0 | est > 1, na.rm = TRUE)) {
    .ckstop(fn, "`est` must be proportions, from 0 to 1")
  }
  if (any(se < 0, na.rm = TRUE)) {
    .ckstop(fn, "`se` must be numbers from 0")
  }
  if (any(n < 2 | n != floor(n), na.rm = TRUE)) {
    .ckstop(fn, "`n` must be whole numbers of at least 2")
  }
  if (any(df < 1, na.rm = TRUE)) {
    .ckstop(fn, "`df` must be at least 1")
  }
  if (
    !is.numeric(level) ||
      length(level) != 1L ||
      !is.finite(level) ||
      level <= 0 ||
      level >= 1
  ) {
    .ckstop(fn, "`level` must be one number in (0, 1)")
  }
  len <- c(length(est), length(se), length(n), length(df))
  k <- max(len)
  if (!all(len %in% c(1L, k))) {
    .ckstop(fn, "the arguments must have one length, or length 1")
  }
  list(
    est = rep_len(est, k),
    se = rep_len(se, k),
    n = rep_len(n, k),
    df = rep_len(df, k)
  )
}

#' The lower and the upper limits of svypci() by `method`, at the
#' effective sample size `ne`.
#' @noRd
.svyint <- function(method, est, se, ne, tq, a) {
  x <- est * ne
  switch(
    method,
    beta = list(
      suppressWarnings(
        ifelse(x > 0, stats::qbeta(a / 2, x, ne - x + 1), 0)
      ),
      suppressWarnings(
        ifelse(x < ne, stats::qbeta(1 - a / 2, x + 1, ne - x), 1)
      )
    ),
    wilson = {
      z <- stats::qnorm(1 - a / 2)
      mid <- est + z^2 / (2 * ne)
      half <- z * sqrt(est * (1 - est) / ne + z^2 / (4 * ne^2))
      list((mid - half) / (1 + z^2 / ne), (mid + half) / (1 + z^2 / ne))
    },
    logit = {
      in01 <- !is.na(est) & est > 0 & est < 1
      l <- ifelse(in01, log(est / (1 - est)), NA_real_)
      h <- tq * se / (est * (1 - est))
      list(1 / (1 + exp(-(l - h))), 1 / (1 + exp(-(l + h))))
    },
    wald = list(est - tq * se, est + tq * se)
  )
}

#' Release quality of a survey estimate
#'
#' Grades an estimate by the rule of Health Canada's guides to the
#' Canadian surveys of alcohol and drug use (Health Canada 2022, Table
#' 1): by the unweighted number of observations behind it and by its
#' coefficient of variation. `"acceptable"`: at least `nmin` observations
#' with the characteristic (and, for a proportion, at least `dmin` in the
#' domain) and a coefficient of variation of at most `cvlo` per cent.
#' `"marginal"`: the same counts and a coefficient above `cvlo` up to
#' `cvhi`; to be released with a warning. `"unacceptable"`: fewer
#' observations, or a coefficient above `cvhi`; not to be released.
#'
#' The thresholds are arguments: other agencies set others, and Statistics
#' Canada's own have changed over the years.
#'
#' Errors: `cv` not numeric or negative; `n` or `d` not numeric or
#' negative; a threshold that is not one number from 0, or `cvlo` above
#' `cvhi`; arguments of different lengths (a length of 1 is repeated).
#'
#' @param cv The coefficient of variation of the estimate, in per cent (as
#'   [svymx()] reports it).
#' @param n The unweighted number of observations with the characteristic
#'   (the numerator of a proportion; for a mean or a total, the
#'   observations that contribute).
#' @param d For a proportion, the unweighted number of observations in the
#'   domain (the denominator); `NULL` (the default) for a mean or a
#'   total.
#' @param cvlo,cvhi The coefficients of variation, in per cent, up to
#'   which an estimate is acceptable and marginal.
#' @param nmin,dmin The fewest observations: with the characteristic, and
#'   in the domain.
#' @return A character vector: `"acceptable"`, `"marginal"` or
#'   `"unacceptable"`; `NA` where the counts suffice and `cv` is missing.
#' @references
#'   Health Canada (2022). *Canadian Postsecondary Education Alcohol and
#'   Drug Use Survey, 2021-22: Public Use Microdata File user guide*.
#'   Section 7.3 and Table 1, quality level guidelines for estimates.
#' @seealso [svymx()] and [svypci()].
#' @examples
#' svyqual(cv = c(8, 20, 40, 8, 8), n = c(120, 60, 30, 6, 40),
#'         d = c(900, 500, 300, 200, 45))
#' # a mean or a total: no denominator
#' svyqual(12, 35)
#' try(svyqual(-1, 10))
#' @export
svyqual <- function(
  cv,
  n,
  d = NULL,
  cvlo = 16.6,
  cvhi = 33.3,
  nmin = 10,
  dmin = 50
) {
  fn <- "svyqual"
  num <- function(v, what) {
    if (!is.numeric(v) || !length(v)) {
      .ckstop(fn, "`%s` must be numeric", what)
    }
    if (any(v < 0, na.rm = TRUE)) {
      .ckstop(fn, "`%s` must be numbers from 0", what)
    }
    as.numeric(v)
  }
  cv <- num(cv, "cv")
  n <- num(n, "n")
  if (!is.null(d)) {
    d <- num(d, "d")
  }
  for (k in c("cvlo", "cvhi", "nmin", "dmin")) {
    v <- get(k)
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 0) {
      .ckstop(fn, "`%s` must be one number from 0", k)
    }
  }
  if (cvlo > cvhi) {
    .ckstop(fn, "`cvlo` must not be above `cvhi`")
  }
  len <- c(length(cv), length(n), if (!is.null(d)) length(d))
  k <- max(len)
  if (!all(len %in% c(1L, k))) {
    .ckstop(fn, "the arguments must have one length, or length 1")
  }
  cv <- rep_len(cv, k)
  n <- rep_len(n, k)
  few <- n < nmin
  if (!is.null(d)) {
    few <- few | rep_len(d, k) < dmin
  }
  out <- ifelse(
    cv <= cvlo,
    "acceptable",
    ifelse(cv <= cvhi, "marginal", "unacceptable")
  )
  out[!is.na(few) & few] <- "unacceptable"
  out[is.na(few) & !(!is.na(cv) & cv > cvhi)] <- NA_character_
  out
}
