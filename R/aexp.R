#' Exposure-response estimate with a bootstrap sensitivity interval
#'
#' Fits the partially linear exposure-response model (the effect of an
#' exposure on an outcome, the confounders partialled out by cross-fitted
#' regressions), or any estimator with the same interface, then re-fits it
#' on \code{nboot} non-parametric resamples of the rows for a percentile
#' interval (Efron and Tibshirani 1993, chapter 13). The bootstrap interval
#' does not lean on the estimator's own standard error, so it is a check on
#' it: when the two disagree, the analytic one should not be trusted alone.
#'
#' In environmental criminology the exposure is a pollutant or heat and the
#' outcome a count of offences or victimisations by place and period, with
#' the neighbourhood's population and deprivation as confounders: the
#' estimate is the change in the outcome per unit of exposure.
#'
#' @details
#' The built-in estimator (\code{estr = NULL}) is the partially linear
#' regression \eqn{y = \theta d + g(x) + e} fitted by double machine
#' learning with cross-fitting (the DML2 estimator of Chernozhukov et al.
#' 2018): the rows with every column present are split at random into five
#' folds, the outcome and the exposure are each regressed on the confounders
#' by ordinary least squares on four folds and predicted on the fifth, and
#' \eqn{\theta} is the regression of the outcome's residuals on the
#' exposure's, with the heteroskedasticity-robust standard error. The folds
#' are drawn from R's generator seeded with 1 for every fit, the stream then
#' restored, so the estimate does not depend on the session. It refuses a
#' column that is not numeric, fewer than ten complete rows, and an exposure
#' with no variation left after the confounders are partialled out.
#'
#' Each resample draws \code{nrow(data)} rows with replacement and calls
#' \code{estr(resample, y, expo, confs)}; a resample on which the estimator
#' fails is skipped, and the run stops if fewer than 10 succeed. The interval
#' is the 2.5\% and 97.5\% type-7 quantiles of the replicate estimates, and
#' \code{seboot} their standard deviation. Resampling uses R's generator
#' seeded with \code{rstate} and restores R's own random stream afterwards.
#' The cost is \code{nboot + 1} estimator calls; the built-in estimator runs
#' in the compiled core.
#'
#' Refusals: column names that are not strings, an \code{nboot} that is not a
#' single non-negative whole number, an estimator that is not a function, and
#' the estimator's own errors on the full sample.
#'
#' @param data A data frame, one row per unit (an area and period, a person).
#' @param y,expo,confs Column names: the outcome and the exposure (single
#'   strings) and the confounders (one or more strings).
#' @param nboot Number of resamples, a non-negative whole number (default
#'   100); at least 10 must succeed.
#' @param rstate Seed for the resampling (default 42), a whole number from
#'   -2147483647 to \eqn{2^{53} - 1}.
#' @param estr \code{NULL} (the default) for the built-in cross-fitted
#'   partially linear regression, or a function (or the name of one)
#'   \code{(data, outcome, treatment, covariates)} returning a list with
#'   \code{ate} and \code{se} (and optionally \code{method}).
#' @param x A result of \code{exprspplr()}.
#' @param ... Not used.
#' @return A list of class \code{exprspplr} with \code{ate} (the
#'   full-sample estimate), \code{seana} (the estimator's standard error),
#'   \code{seboot} (the standard deviation of the replicates), \code{cilobs}
#'   and \code{ciupbs} (the 95\% percentile interval), \code{nboot} (the
#'   number of successful resamples) and \code{method} (the estimator's
#'   method with \code{" + percentile bootstrap"}).
#' @references Efron, B. and Tibshirani, R. J. (1993). An Introduction to the
#'   Bootstrap. Chapman and Hall.
#'
#'   Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C.,
#'   Newey, W. and Robins, J. (2018). Double/debiased machine learning for
#'   treatment and structural parameters. The Econometrics Journal 21(1),
#'   C1-C68. \doi{10.1111/ectj.12097}
#' @seealso \code{\link{dmlx}} for double machine learning with other
#'   learners and models.
#' @examples
#' # simulated areas: exposure driven by x1, the outcome by the exposure and
#' # x2
#' set.seed(1)
#' d <- data.frame(x1 = rnorm(200), x2 = rnorm(200))
#' d$exposure <- 20 + 2 * d$x1 + rnorm(200)
#' d$asthma <- 5 + 0.3 * d$exposure + d$x2 + rnorm(200)
#' r <- exprspplr(d, y = "asthma", expo = "exposure", confs = c("x1", "x2"),
#'                nboot = 20)
#' r$ate
#' c(r$cilobs, r$ciupbs)
#'
#' # Any estimator with the same interface: here an ordinary regression.
#' ols <- function(data, outcome, treatment, covariates) {
#'   f <- stats::reformulate(c(treatment, covariates), outcome)
#'   cf <- summary(stats::lm(f, data = data))$coefficients
#'   list(ate = cf[treatment, 1], se = cf[treatment, 2], method = "OLS")
#' }
#' exprspplr(d, "asthma", "exposure", c("x1", "x2"), nboot = 20,
#'           estr = ols)$method
#'
#' # Fewer than 10 resamples cannot give an interval.
#' try(exprspplr(d, "asthma", "exposure", c("x1", "x2"), nboot = 5))
#' @export
exprspplr <- function(
  data,
  y,
  expo,
  confs,
  nboot = 100L,
  rstate = 42L,
  estr = NULL
) {
  f <- "exprspplr"
  .colnm(y, "y", f)
  .colnm(expo, "expo", f)
  .colnm(confs, "confs", f, several = TRUE)
  if (is.character(estr) && length(estr) == 1L && !is.na(estr)) {
    nm <- estr
    estr <- tryCatch(match.fun(estr), error = function(e) {
      .ckstop(f, "`estr` names no function: \"%s\"", nm)
    })
  }
  if (!is.null(estr) && !is.function(estr)) {
    .ckstop(
      f,
      "`estr` must be a function (one like the built-in estimator), or NULL"
    )
  }
  rstate <- .ckseed(rstate, "rstate", f)
  nboot <- .polsc(nboot, "nboot", f)
  if (nboot != round(nboot) || nboot > 1e6) {
    .ckstop(f, "`nboot` must be a whole number of resamples, at most 1000000")
  }
  nboot <- as.integer(nboot)
  n <- nrow(data)
  if (is.null(estr)) {
    r <- .plrfit(data, y, expo, confs, nboot, rstate, f)
    base <- r$base
    boot <- r$boot
  } else {
    base <- tryCatch(estr(data, y, expo, confs), error = function(e) {
      .ckstop(
        f,
        paste(
          "`estr` failed on the full sample, called as",
          "estr(data, y, expo, confs): %s"
        ),
        sub("\\.$", "", conditionMessage(e))
      )
    })
    boot <- .clseed(rstate, {
      out <- numeric(0)
      for (b in seq_len(nboot)) {
        idx <- sample.int(n, n, replace = TRUE)
        rb <- tryCatch(
          estr(data[idx, , drop = FALSE], y, expo, confs),
          error = function(e) NULL
        )
        if (!is.null(rb)) out <- c(out, as.numeric(rb$ate))
      }
      out
    })
  }
  if (length(boot) < 10L) {
    .ckstop(
      f,
      paste0(
        "Only %d successful bootstrap fits; need >=10. Inspect the ",
        "model's convergence on your resamples"
      ),
      length(boot)
    )
  }
  s <- .Call(C_clab_r_bootci, as.numeric(boot), c(0.025, 0.975))
  structure(
    list(
      ate = as.numeric(base$ate),
      seana = as.numeric(base$se),
      seboot = s[1L],
      cilobs = s[2L],
      ciupbs = s[3L],
      nboot = length(boot),
      method = paste0(
        if (is.null(base$method)) "exposure-response" else base$method,
        " + percentile bootstrap"
      )
    ),
    class = c("exprspplr", "list")
  )
}

#' @rdname exprspplr
#' @export
print.exprspplr <- function(x, ...) {
  cat(
    "Exposure-response estimate: ", format(x$ate, digits = 4), "\n",
    "  Analytic standard error: ", format(x$seana, digits = 4), "\n",
    "  Bootstrap standard error: ", format(x$seboot, digits = 4),
    " (", x$nboot, " resamples)\n",
    "  95% percentile interval: ", format(x$cilobs, digits = 4), " to ",
    format(x$ciupbs, digits = 4), "\n",
    "  Method: ", x$method, "\n",
    sep = ""
  )
  invisible(x)
}

# The built-in estimator on the full sample and the resamples: the rows and
# folds drawn here with R's generator, the fits in the core.
.plrfit <- function(data, y, expo, confs, nboot, rstate, fn) {
  cols <- c(y, expo, confs)
  miss <- setdiff(cols, names(data))
  if (length(miss)) {
    .ckstop(fn, "data is missing columns: %s", toString(miss))
  }
  keep <- stats::complete.cases(data[, cols, drop = FALSE])
  for (cl in cols) {
    if (!is.numeric(data[[cl]][keep])) {
      .ckstop(
        fn,
        paste0(
          "column `%s` must be numeric (a factor or character column is ",
          "not a code)"
        ),
        cl
      )
    }
  }
  n <- nrow(data)
  if (sum(keep) < 10L) {
    .ckstop(
      fn,
      "need at least 10 complete rows for 5 folds; have %d",
      sum(keep)
    )
  }
  yv <- as.numeric(data[[y]])
  dv <- as.numeric(data[[expo]])
  xm <- matrix(
    as.numeric(unlist(data[confs], use.names = FALSE)),
    nrow = n
  )
  # the folds plrcf() draws: seed 1 for every fit, so one set per size
  fd <- new.env(parent = emptyenv())
  folds <- function(m) {
    k <- as.character(m)
    if (is.null(fd[[k]])) {
      fd[[k]] <- .clseed(1L, sample(rep_len(seq_len(5L), m)))
    }
    fd[[k]]
  }
  rows0 <- which(keep)
  idx <- .clseed(rstate, lapply(seq_len(nboot), function(b) {
    sample.int(n, n, replace = TRUE)
  }))
  rb <- lapply(idx, function(i) i[keep[i]])
  ok <- lengths(rb) >= 10L
  smp <- c(list(rows0), rb[ok])
  fit <- .Call(
    C_clab_r_plr, yv, dv, xm, as.integer(unlist(smp)),
    as.integer(lengths(smp)),
    as.integer(unlist(lapply(lengths(smp), folds)))
  )
  if (is.na(fit[[1L]][1L])) {
    .ckstop(
      fn,
      paste(
        "the treatment has no variation left after partialling out the",
        "covariates"
      )
    )
  }
  est <- fit[[1L]][-1L]
  list(
    base = list(
      ate = fit[[1L]][1L],
      se = fit[[2L]][1L],
      method = "partially linear regression, DML2 cross-fitting, OLS nuisances"
    ),
    boot = est[!is.na(est)]
  )
}
