#' Multiple imputation by predictive mean matching
#'
#' Fills the missing values of numeric columns `m` times by predictive
#' mean matching, each completed data set from its own chain of chained
#' equations, and pools an analysis of the completed data sets by Rubin's
#' rules: the step after [mcartst()] has said how the values are missing.
#' Predictive mean matching imputes only values that were observed, so
#' an imputed count stays a whole number and an imputed age stays in its
#' range (van Buuren 2018, sec. 3.4).
#'
#' Each chain starts the missing values of every column of `y` at draws
#' from its observed values and, `maxit` times, imputes each column in
#' turn from all the others (the other columns of `y` as last imputed, and
#' the complete columns `x`): the least squares fit on its observed rows,
#' `ridge` times the diagonal of the cross-products added to it; a draw of
#' the residual scale (`rss` over a chi-square on the residual degrees of
#' freedom) and of the coefficients about their estimate; and, for each
#' missing row, a donor drawn with equal chances among the `k` observed
#' rows whose fitted values (from the estimate) are nearest its own (from
#' the draw), whose observed value it takes (van Buuren 2018, sec. 4.5;
#' type 1 matching). The draws come from the package's generator.
#'
#' The analysis is the least squares fit of `fy` on an intercept and `fx`
#' (classical variances), or, without `fy`, the mean of each column of
#' `y` (variance `s^2 / n`). Pooled by Rubin (1987): `est` the mean of the
#' `m` estimates, `se` the square root of `T = U + (1 + 1 / m) B` (`U` the
#' mean variance, `B` the variance between imputations), `df` that of
#' Barnard and Rubin (1999) with the complete-data degrees of freedom of
#' the analysis, Student's t interval and p-value, `riv` the relative
#' increase in variance `(1 + 1 / m) B / U` and `fmi` the fraction of
#' missing information, as the `mice` package pools.
#'
#' @inheritParams cemtch
#' @inheritParams rndalc
#' @inheritParams plntfx
#' @param y The names of the numeric columns to impute (missing values as
#'   `NA`; observed values finite).
#' @param x `NULL` or the names of complete numeric columns that predict
#'   them; for `print()`, an object of class `mipmm`.
#' @param m The imputations: a whole number from 2 to 1000.
#' @param k The donors to draw from: a whole number from 1 to 10000.
#' @param maxit The rounds of each chain: a whole number from 1 to 1000.
#' @param ridge The ridge, as a share of the diagonal: a number from 0.
#' @param fy `NULL` or the outcome of the analysis: a column of `y` or `x`.
#' @param fx `NULL` or its predictors: columns of `y` or `x` other than
#'   `fy`.
#' @return An object of class `mipmm`, a list: `imp`, for each column of
#'   `y` a matrix of its imputed values, one row per missing value (named
#'   by its row of `data`) and one column per imputation; `pool`, a data
#'   frame (`term`, `est`, `se`, `df`, `lo`, `hi`, `p`, `riv`, `fmi`);
#'   `nmis`, the missing values of each column; and `m`, `k`, `maxit`,
#'   `seed`.
#' @references
#'   van Buuren, S. (2018). *Flexible Imputation of Missing Data*, 2nd ed.
#'   Chapman and Hall/CRC. \doi{10.1201/9780429492259}
#'
#'   Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in
#'   Surveys*. Wiley. \doi{10.1002/9780470316696}
#'
#'   Barnard, J. and Rubin, D. B. (1999). Small-sample degrees of freedom
#'   with multiple imputation. *Biometrika*, 86(4), 948-955.
#'   \doi{10.1093/biomet/86.4.948}
#'
#'   Little, R. J. A. (1988). Missing-data adjustments in large surveys.
#'   *Journal of Business and Economic Statistics*, 6(3), 287-296.
#'   \doi{10.1080/07350015.1988.10509663}
#' @examples
#' set.seed(11)
#' n <- 200
#' d <- data.frame(age = round(runif(n, 18, 60)), male = rbinom(n, 1, 0.5))
#' d$priors <- rpois(n, exp(-1 + 0.04 * d$age))
#' d$days <- round(20 + 3 * d$priors + 5 * d$male + rnorm(n, sd = 5))
#' d$priors[sample(n, 40)] <- NA
#' d$days[sample(n, 30)] <- NA
#' r <- mipmm(d, c("priors", "days"), c("age", "male"), fy = "days",
#'            fx = c("priors", "male"), seed = 1)
#' r
#' head(r$imp$priors)
#' @export
mipmm <- function(
  data,
  y,
  x = NULL,
  m = 5,
  k = 5,
  maxit = 5,
  ridge = 1e-5,
  fy = NULL,
  fx = NULL,
  level = 0.95,
  seed = NULL
) {
  f <- "mipmm"
  dd <- .axdat(data, f)
  .mpcols(dd, y, x, f)
  ym <- vapply(y, function(v) .mpval(dd[[v]], v, f), numeric(nrow(dd)))
  ym <- matrix(ym, nrow(dd))
  xm <- .axmat(dd, x, f)
  m <- .cknum(m, "m", f, min = 2, max = 1000, whole = TRUE)
  k <- .cknum(k, "k", f, min = 1, max = 1e4, whole = TRUE)
  maxit <- .cknum(maxit, "maxit", f, min = 1, max = 1000, whole = TRUE)
  ridge <- .cknum(ridge, "ridge", f, min = 0)
  a <- c(y, x)
  .mpfml(fy, fx, a, f)
  level <- .xpopen(level, "level", f)
  seed <- .xpseed(seed, f)
  o <- .Call(C_clab_r_mipmm, cbind(ym, xm), seq_along(y) - 1L, as.integer(m),
             as.integer(k), as.integer(maxit), ridge, seed,
             if (is.null(fy)) -1L else match(fy, a) - 1L,
             match(fx, a) - 1L, level)
  term <- if (is.null(fy)) y else c("(Intercept)", fx)
  structure(
    list(imp = .mpimp(o$imp, ym, y, m), pool = .mppool(o$pool, term),
         nmis = stats::setNames(colSums(is.na(ym)), y), m = m, k = k,
         maxit = maxit, seed = seed),
    class = "mipmm"
  )
}

#' @rdname mipmm
#' @inheritParams cemtch
#' @export
print.mipmm <- function(x, ...) {
  cat("mipmm: ", sum(x[["nmis"]]), " missing value(s) in ",
      length(x[["nmis"]]), " column(s), ", x[["m"]],
      " imputations by predictive mean matching\n", sep = "")
  print(x[["pool"]], row.names = FALSE, digits = 4)
  invisible(x)
}

# A column to impute: numbers, finite where observed.
.mpval <- function(v, nm, f) {
  if (!is.numeric(v) && !is.logical(v)) {
    .ckstop(f, "column `%s` must be numeric, got %s", nm, class(v)[1L])
  }
  v <- as.numeric(v)
  if (any(!is.na(v) & !is.finite(v)) || any(is.nan(v))) {
    .ckstop(f, "column `%s` must be finite where it is observed", nm)
  }
  v
}

# The columns y to impute and x that predict them, checked.
.mpcols <- function(dd, y, x, f) {
  if (is.null(y)) .ckstop(f, "`y` must name columns of `data`, got NULL")
  if (!is.character(y) || anyNA(y) || !length(y)) {
    .ckstop(f, "`y` must name columns of `data`")
  }
  miss <- setdiff(y, names(dd))
  if (length(miss)) {
    .ckstop(f, "`y` names columns not in `data`: %s", toString(miss))
  }
  if (anyDuplicated(y)) .ckstop(f, "`y` names a column twice")
  .axcols(dd, x, "x", f, null = TRUE)
  if (any(y %in% x)) .ckstop(f, "`y` and `x` must name different columns")
}

# The analysis checked: fy one of the columns a, fx others of them.
.mpfml <- function(fy, fx, a, f) {
  if (!is.null(fy) && (!is.character(fy) || length(fy) != 1L ||
                         !fy %in% a)) {
    .ckstop(f, "`fy` must name one column of `y` or `x`")
  }
  if (!is.null(fx) && (is.null(fy) || !is.character(fx) || anyNA(fx) ||
                         !all(fx %in% setdiff(a, fy)) ||
                         anyDuplicated(fx))) {
    .ckstop(f, paste("`fx` must name columns of `y` or `x` other than `fy`,",
                     "each once"))
  }
}

# The imputed values of each column of y, one matrix each.
.mpimp <- function(v, ym, y, m) {
  imp <- list()
  at <- 0
  for (j in seq_along(y)) {
    r <- which(is.na(ym[, j]))
    imp[[y[j]]] <- matrix(v[at + seq_len(length(r) * m)], length(r), m,
                          dimnames = list(r, NULL))
    at <- at + length(r) * m
  }
  imp
}

# The pooled table from the kernel's pooled values p.
.mppool <- function(p, term) {
  np <- length(term)
  data.frame(term = term, est = p[seq_len(np)], se = p[np + seq_len(np)],
             df = p[2L * np + seq_len(np)], lo = p[3L * np + seq_len(np)],
             hi = p[4L * np + seq_len(np)], p = p[5L * np + seq_len(np)],
             riv = p[6L * np + seq_len(np)], fmi = p[7L * np + seq_len(np)],
             stringsAsFactors = FALSE)
}
