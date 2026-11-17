#' Long-run effects from a common correlated effects fit
#'
#' The cross-sectionally augmented autoregressive distributed lag
#' (CS-ARDL) long-run effects of Chudik et al. (2016) as Ditzen (2021,
#' sec. 3.2) sets them out, from the unit estimates of [ccepnl()].
#'
#' A dynamic panel `y_it = mu_i + sum_l lambda_l,i y_i,t-l + sum_l
#' beta_l,i x_i,t-l + averages + e_it` fitted unit by unit (`ccepnl()`
#' with `ylag` lags of `y` and the lags of `x` as columns of its `x`) has
#' the long-run effect of `x` on the steady state of `y`
#' `theta_i = sum_l beta_l,i / (1 - sum_l lambda_l,i)` (Ditzen 2021, eq. 3
#' and sec. 3.2). Each unit's `theta_i` is computed from its estimates;
#' the mean group estimate is their mean and its variance the mean
#' group variance `sum_i (theta_i - mean)^2 / (N (N - 1))` (Pesaran and
#' Smith 1995, as Chudik et al. 2016 apply it). The adjustment speed
#' `phi_i = 1 - sum_l lambda_l,i` of the error-correction form (Ditzen
#' 2021, sec. 3.1) is reported the same way as `"phi"`.
#'
#' The cross-sectionally augmented distributed lag (CS-DL) estimator of
#' the same paper needs no function of its own: `ccepnl()` with `x` and
#' its first differences `x_t - x_t-1`, ..., lagged up to
#' `floor(T^(1/3)) - 1`, as regressors and `ylag = 0` estimates `theta` as the
#' coefficient of `x` (Ditzen 2021, sec. 3.3); see the examples. Ditzen (2021,
#' sec. 3.3) also lags the averages of `x` up to `floor(T^(1/3))` and keeps the
#' average of `y` unlagged, where `ccepnl(lags = )` lags both.
#'
#' The exponent of cross-sectional dependence of Bailey, Kapetanios and
#' Pesaran (2016), the other topic of Ditzen (2021), is not built here
#' (its bias correction is not set out in Ditzen 2021);
#' Pesaran's CD test is [pcdtst()].
#'
#' @param object An object of class `ccepnl`, or a data frame or matrix
#'   of unit estimates with one named column per term (as its `bi`).
#' @param x The terms of one regressor and its lags (columns of `bi`), or
#'   a named list of such sets, one long-run effect each.
#' @param ylag `NULL` (the columns named `L1.y`, `L2.y`, ... of `bi`) or
#'   the names of the columns holding the lags of `y`.
#' @inheritParams ppmlfe
#' @return A data frame, one row per long-run effect and a last row
#'   `"phi"`: `term`, `est` (the mean group estimate), `se`, `z`, `p`,
#'   `lo`, `hi`; with the unit values as the attribute `units` (a matrix,
#'   units by rows of the table).
#' @references
#'   Ditzen, J. (2021). Estimating long-run effects and the exponent of
#'   cross-sectional dependence: an update to xtdcce2. *The Stata Journal*,
#'   21(3), 687-707. \doi{10.1177/1536867X211045560}
#'
#'   Chudik, A., Mohaddes, K., Pesaran, M. H. and Raissi, M. (2016).
#'   Long-run effects in large heterogeneous panel data models with
#'   cross-sectionally correlated errors. In *Essays in Honor of Aman
#'   Ullah* (Advances in Econometrics, vol. 36), 85-135. Emerald.
#'
#'   Pesaran, M. H. and Smith, R. (1995). Estimating long-run
#'   relationships from dynamic heterogeneous panels. *Journal of
#'   Econometrics*, 68, 79-113. \doi{10.1016/0304-4076(94)01644-F}
#'
#'   Bailey, N., Kapetanios, G. and Pesaran, M. H. (2016). Exponent of
#'   cross-sectional dependence: estimation and inference. *Journal of
#'   Applied Econometrics*, 31, 929-960. \doi{10.1002/jae.2476}
#' @examples
#' set.seed(3)
#' tt <- 30
#' d <- expand.grid(year = 1:tt, prov = paste0("p", 1:25))
#' f <- rnorm(tt)[d$year]
#' g <- rnorm(25, 1)[as.integer(substring(d$prov, 2))]
#' d$x <- g * f + rnorm(nrow(d))
#' d$xl <- ave(d$x, d$prov, FUN = function(v) c(NA, v[-length(v)]))
#' d$y <- 0
#' for (i in which(d$year > 1)) {
#'   d$y[i] <- 0.4 * d$y[i - 1] + 0.5 * d$x[i] + 0.1 * d$xl[i] +
#'     g[i] * f[i] + rnorm(1)
#' }
#' k <- d$year > 1
#' # CS-ARDL: theta = (0.5 + 0.1) / (1 - 0.4) = 1
#' r <- ccepnl(d$y[k], d[k, c("x", "xl")], d$prov[k], d$year[k], ylag = 1)
#' lrcce(r, c("x", "xl"))
#' # CS-DL: x and its first difference; the coefficient of x is theta
#' d$dx <- d$x - d$xl
#' ccepnl(d$y[k], d[k, c("x", "dx")], d$prov[k], d$year[k], lags = 3)$mg
#' @export
lrcce <- function(object, x, ylag = NULL, level = 0.95) {
  f <- "lrcce"
  b <- .lrcb(object, f)
  a <- .lrcx(b, x, ylag, f)
  x <- a$x
  ylag <- a$ylag
  level <- .axlevel(level, f)
  phi <- 1 - rowSums(b[, ylag, drop = FALSE])
  u <- cbind(vapply(x, function(s) rowSums(b[, s, drop = FALSE]) / phi,
                    numeric(nrow(b))), phi = phi)
  if (!all(is.finite(u))) .ckstop(f, "a unit's lags of `y` sum to 1")
  n <- nrow(u)
  se <- sqrt(colSums(sweep(u, 2L, colMeans(u))^2) / (n * (n - 1)))
  tab <- .axtab(colnames(u), colMeans(u), se, level)
  names(tab)[names(tab) == "stat"] <- "z"
  attr(tab, "units") <- u
  tab
}

# The unit estimates of lrcce()'s object, checked.
.lrcb <- function(object, f) {
  b <- if (inherits(object, "ccepnl")) object[["bi"]] else object
  if (is.data.frame(b)) b <- as.matrix(b)
  if (!is.matrix(b) || !is.numeric(b) || is.null(colnames(b))) {
    .ckstop(f, "`object` must be a ccepnl fit or a matrix of unit %s",
            "estimates with named columns")
  }
  if (nrow(b) < 2L) .ckstop(f, "`object` needs two or more units")
  b
}

# The terms x (a list, named by their first columns when unnamed) and
# the lags of y (by default the columns L<k>.y) of lrcce(), checked
# against the unit estimates b: list(x, ylag).
.lrcx <- function(b, x, ylag, f) {
  if (!is.list(x)) x <- list(x)
  if (is.null(ylag)) ylag <- grep("^L[0-9]+\\.y$", colnames(b), value = TRUE)
  for (s in c(x, list(ylag))) {
    if (!is.character(s) || anyNA(s) || !all(s %in% colnames(b))) {
      .ckstop(f, "`x` and `ylag` must name columns of the unit estimates")
    }
  }
  if (is.null(names(x))) names(x) <- vapply(x, `[`, "", 1L)
  if (!length(ylag)) .ckstop(f, "the fit has no lags of `y`")
  list(x = x, ylag = ylag)
}
