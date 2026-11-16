#' The generalised ordered logit
#'
#' An ordered outcome (a risk level, a charge's seriousness, a sentence
#' band) regressed on covariates without forcing each covariate to act the
#' same way at every cut of the scale. For an outcome with levels
#' `1, ..., M`, each threshold `j < M` has its own intercept and slopes,
#' `P(y > j) = F(a_j + x b_j)` with `F` the logistic distribution function
#' (Peterson and Harrell 1990; Williams 2006, the parameterisation of
#' 'gologit2'); a covariate named in `par` keeps one slope for all the
#' thresholds, as in the proportional-odds model, so that `par = x` is that
#' model. O'Connell and Laniyonu (2025) find racial disparity in a risk
#' tool concentrated at the low-to-medium cut and absent at the
#' medium-to-high one, which a proportional-odds fit averages away.
#'
#' The parameters are estimated by maximum likelihood (quasi-Newton steps
#' with the analytic gradient, in the compiled core), with standard errors
#' from the inverse of the Hessian. A free covariate's slopes may cross,
#' so that a category's fitted probability is negative for some `x`; the
#' fit refuses such points within the data. `wald` tests, for each free
#' covariate, that its slopes are equal at every threshold (the test
#' 'gologit2' uses to choose the parallel ones), and `lr` is the
#' likelihood-ratio test of the fit against the proportional-odds model.
#'
#' @inheritParams cemtch
#' @inheritParams spcure
#' @param y The name of the ordered outcome column: an ordered factor or a
#'   factor (its levels in order), numbers (in increasing order) or text
#'   (with `lev`).
#' @param x The names of the covariate columns: numeric or logical, at
#'   least one.
#' @param par `NULL`, or the names of covariates in `x` held parallel (one
#'   slope for every threshold).
#' @param lev `NULL`, or the levels of the outcome from lowest to highest,
#'   each present in the data; needed for text.
#' @return An object of class `gologt`, a list: `coef`, one row per
#'   parameter (`term`, `thr` the threshold, as `"low|med"`, or `"all"` for
#'   a parallel slope; `est`, `se`, `z`, `p`, `lo`, `hi`); `wald`, one row
#'   per covariate (`term`, `stat`, `df`, `p`; `NA` for a parallel one);
#'   `lr`, the test against proportional odds (`stat`, `df`, `p`; `NA` when
#'   every slope is parallel); `fit` (`loglik`, `aic`, `iter`, `conv`, 1
#'   when converged); `vcov`; `lev`; and `n`.
#' @references
#'   Peterson, B. and Harrell, F. E. (1990). Partial proportional odds
#'   models for ordinal response variables. *Journal of the Royal
#'   Statistical Society Series C*, 39(2), 205-217. \doi{10.2307/2347760}
#'
#'   Williams, R. (2006). Generalized ordered logit/partial proportional
#'   odds models for ordinal dependent variables. *The Stata Journal*,
#'   6(1), 58-82. \doi{10.1177/1536867X0600600104}
#'
#'   O'Connell, M. and Laniyonu, A. (2025). Threshold-specific
#'   cumulative-logit models for actuarial-risk audit. *Race and Justice*,
#'   15(3), 428-453.
#' @seealso [gocoef()] for the slopes by threshold.
#' @examples
#' set.seed(7)
#' n <- 400
#' d <- data.frame(race = rbinom(n, 1, 0.4), age = rnorm(n))
#' e1 <- 0.5 + 0.8 * d$race - 0.5 * d$age
#' e2 <- -1 + 0.1 * d$race - 0.5 * d$age
#' u <- runif(n)
#' d$risk <- factor(1 + (u < plogis(e1)) + (u < plogis(e2)), 1:3,
#'                  c("low", "med", "high"))
#' r <- gologt(d, "risk", c("race", "age"), par = "age")
#' r
#' @export
gologt <- function(data, y, x, par = NULL, lev = NULL, level = 0.95,
                   maxit = 500) {
  f <- "gologt"
  dd <- .axdat(data, f)
  .axcols(dd, y, "y", f, n = 1L)
  .axcols(dd, x, "x", f)
  xm <- .axmat(dd, x, f)
  .gopar(par, x, f)
  lv <- .golev(dd[[y]], lev, y, f)
  yc <- match(as.character(dd[[y]]), as.character(lv)) - 1L
  level <- .cknum(level, "level", f)
  if (level <= 0 || level >= 1) {
    .ckstop(f, "`level` must be strictly between 0 and 1")
  }
  maxit <- .cknum(maxit, "maxit", f, min = 1, whole = TRUE)
  M <- length(lv)
  pv <- as.integer(x %in% par)
  k <- .Call(C_clab_r_gologt, yc, M, xm, pv, min(maxit, 2147483647), level)
  lb <- paste(lv[-M], lv[-1L], sep = "|")
  term <- c(rep("(Intercept)", M - 1L),
            unlist(lapply(seq_along(x), function(j) {
              rep(x[j], if (pv[j] == 1L) 1L else M - 1L)
            })))
  thr <- c(lb, unlist(lapply(pv, function(v) if (v == 1L) "all" else lb)))
  cf <- k$coef
  nm <- paste(term, thr, sep = ":")
  v <- k$vcov
  dimnames(v) <- list(nm, nm)
  fit <- k$fit
  n <- length(yc)
  structure(
    list(coef = data.frame(term = term, thr = thr, est = cf[, 1L],
                           se = cf[, 2L], z = cf[, 3L], p = cf[, 4L],
                           lo = cf[, 5L], hi = cf[, 6L]),
         wald = data.frame(term = x, stat = k$wald[1L, ], df = k$wald[2L, ],
                           p = k$wald[3L, ]),
         lr = c(stat = k$lr[1L], df = k$lr[2L], p = k$lr[3L]),
         fit = c(loglik = fit[1L], aic = -2 * fit[1L] + 2 * nrow(cf),
                 iter = fit[2L], conv = fit[3L]),
         vcov = v, lev = as.character(lv), n = n),
    class = "gologt"
  )
}

#' @rdname gologt
#' @inheritParams cemtch
#' @export
print.gologt <- function(x, ...) {
  f <- x[["fit"]]
  cat("gologt: generalised ordered logit, ", x[["n"]], " rows, ",
      length(x[["lev"]]), " levels; log-likelihood ",
      format(signif(f[["loglik"]], 6)),
      if (f[["conv"]] != 1) " (not converged)", "\n", sep = "")
  print(x[["coef"]], row.names = FALSE)
  l <- x[["lr"]]
  if (!is.na(l[["stat"]])) {
    cat("against proportional odds: LR ", format(signif(l[["stat"]], 4)),
        " on ", l[["df"]], " df, p ", format(signif(l[["p"]], 3)), "\n",
        sep = "")
  }
  invisible(x)
}

#' The slopes of a generalised ordered logit by threshold
#'
#' The coefficients of [gologt()] laid out as a table: one row per
#' threshold of the outcome, one column per term, so that a covariate's
#' slopes can be read across the cuts of the scale (a parallel covariate
#' repeats its one slope), with their standard errors and, for the free
#' covariates, the Wald test that the slopes are equal.
#'
#' @param object A result of [gologt()].
#' @param v `NULL` for every term, or the names of terms: `"(Intercept)"`
#'   and the covariates.
#' @return A list: `est` and `se`, matrices with one row per threshold
#'   (named as `"low|med"`) and one column per term of `v`; `wald`, the rows
#'   of the fit's Wald tests for the covariates of `v`.
#' @seealso [gologt()].
#' @examples
#' set.seed(8)
#' n <- 300
#' d <- data.frame(race = rbinom(n, 1, 0.4), age = rnorm(n))
#' u <- runif(n)
#' d$risk <- 1 + (u < plogis(0.5 + 0.8 * d$race)) +
#'   (u < plogis(-1 + 0.1 * d$race - 0.3 * d$age))
#' gocoef(gologt(d, "risk", c("race", "age")), "race")
#' @export
gocoef <- function(object, v = NULL) {
  f <- "gocoef"
  if (!inherits(object, "gologt")) {
    .ckstop(f, "`object` must be a result of gologt()")
  }
  cf <- object[["coef"]]
  tm <- unique(cf[["term"]])
  if (is.null(v)) v <- tm
  if (!is.character(v) || !length(v) || anyNA(v) || !all(v %in% tm)) {
    .ckstop(f, "`v` must name terms of the fit: %s", toString(tm))
  }
  lb <- cf[["thr"]][cf[["term"]] == "(Intercept)"]
  tab <- function(col) {
    m <- vapply(v, function(t) {
      r <- cf[[col]][cf[["term"]] == t]
      if (length(r) == 1L) rep(r, length(lb)) else r
    }, numeric(length(lb)))
    matrix(m, length(lb), dimnames = list(lb, v))
  }
  w <- object[["wald"]]
  w <- w[w[["term"]] %in% v, , drop = FALSE]
  rownames(w) <- NULL
  list(est = tab("est"), se = tab("se"), wald = w)
}

# The levels of the ordered outcome v, lowest first: lev, the levels of a
# factor, or the sorted values of numbers; each present, at least three.
.golev <- function(v, lev, y, f) {
  if (!is.null(lev)) {
    if (!(is.character(lev) || is.numeric(lev)) || anyNA(lev) ||
          anyDuplicated(lev)) {
      .ckstop(f, "`lev` must be distinct levels with no missing values")
    }
    lv <- lev
  } else if (is.factor(v)) {
    lv <- levels(v)
  } else if (is.numeric(v) || is.logical(v)) {
    lv <- sort(unique(as.numeric(v)))
  } else {
    .ckstop(f, "`lev` must give the order of the levels of column `%s`, text",
            y)
  }
  s <- as.character(v)
  l <- as.character(lv)
  if (!all(s %in% l)) {
    .ckstop(f, "column `%s` has values not in its levels: %s", y,
            toString(utils::head(unique(s[!s %in% l]), 5L)))
  }
  if (!all(l %in% s)) {
    .ckstop(f, "column `%s` has no rows at level(s) %s", y,
            toString(l[!l %in% s]))
  }
  if (length(lv) < 3L) {
    .ckstop(f, "column `%s` must have at least three levels", y)
  }
  lv
}

# The parallel terms: NULL, or columns of x, each once.
.gopar <- function(par, x, f) {
  if (!is.null(par) && (!is.character(par) || anyNA(par) ||
                          !all(par %in% x) || anyDuplicated(par))) {
    .ckstop(f, "`par` must name columns of `x`, each once")
  }
}
