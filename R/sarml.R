#' Spatial regression: lag, error and combined models
#'
#' Fits the spatial lag model `y = rho W y + X beta + e`, the spatial error
#' model `y = X beta + u` with `u = lambda W u + e`, or both at once (`"sac"`,
#' one weights matrix for both), with normal errors, by maximum likelihood
#' (Ord 1975; Anselin 1988): the likelihood is concentrated in `rho` or
#' `lambda` and maximised by Brent's method on `(-1 / m, 1 / m)`, `m` the
#' largest row sum of `W` (which bounds its spectral radius), with `log |I
#' - rho W|` from a dense LU factorisation, so the cost grows as `n^3` per
#' evaluation (at most 5000 units). The combined model is maximised over
#' both parameters by Nelder and Mead's simplex. Standard errors come from
#' the inverse of the expected information matrix (Anselin 1988, ch. 6);
#' these are the estimates of `spatialreg::lagsarlm()`,
#' `errorsarlm()` and `sacsarlm()`. `durbin = TRUE` adds the lags `W x` of
#' the varying columns (the spatial Durbin model, or the Durbin error
#' model). `method = "2sls"` fits the lag model by two-stage least squares
#' with the instruments `W x` and `W^2 x` (Kelejian and Prucha 1998), as
#' `spatialreg::stsls()` does, without a likelihood; its standard errors
#' divide the residual sum of squares by `n` less the coefficients.
#'
#' Next to each coefficient is ordinary least squares on the same columns
#' and its usual standard error: leaving out a real spatial lag biases
#' those coefficients, leaving out a spatial error leaves them unbiased but
#' their standard errors wrong. In a lag model a coefficient is not the
#' effect of its variable: see [spill()].
#'
#' @param y The outcome: finite numbers, one per unit.
#' @param x The covariates: a numeric vector, matrix or data frame with
#'   one row per unit (no intercept column when `int = TRUE`).
#' @inheritParams gmoran
#' @param model `"lag"`, `"error"` or `"sac"`.
#' @param durbin `TRUE` to add the lagged covariates.
#' @param method `"ml"` or `"2sls"` (the lag model only).
#' @param int `TRUE` to add an intercept.
#' @return A list of class `sarml`: `coef`, a data frame of `term` (the
#'   covariates, `lag.` before a lagged one, then `rho` and/or `lambda`),
#'   `est`, `se`, `z`, `p` (two-sided normal), `ols` and `olsse`; `fit`,
#'   a named vector: `loglik` and `aic` (`NaN` for two-stage least
#'   squares), `s2` (the error variance, the residual sum of squares over
#'   `n`), `llols` (the log-likelihood of ordinary least squares) and `n`;
#'   `vcov`, the covariance matrix of the coefficients, rows and columns for
#'   every covariate, `rho` and `lambda` (`NaN` for the one not in the
#'   model); `model`, `method`, `durbin`, `nlag` (the lagged covariates) and
#'   `w`.
#' @references Anselin, L. (1988). *Spatial Econometrics: Methods and
#'   Models*. Kluwer, Dordrecht. \doi{10.1007/978-94-015-7799-1}
#'
#'   Kelejian, H. H. and Prucha, I. R. (1998). A generalized spatial
#'   two-stage least squares procedure for estimating a spatial
#'   autoregressive model with autoregressive disturbances. *Journal of
#'   Real Estate Finance and Economics*, 17(1), 99-121.
#'   \doi{10.1023/A:1007707430416}
#' @examples
#' set.seed(1)
#' p <- cbind(x = runif(80), y = runif(80))
#' w <- spwt(p, k = 5)
#' W <- matrix(0, 80, 80)
#' W[cbind(w$i, w$j)] <- w$w
#' x <- rnorm(80)
#' y <- solve(diag(80) - 0.5 * W, 1 + 2 * x + rnorm(80))
#' f <- sarml(y, x, w)
#' f
#' spill(f)
#' @export
sarml <- function(y, x, w, model = "lag", durbin = FALSE, method = "ml",
                  int = TRUE) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x))) x <- matrix(x, dimnames = list(NULL, "x"))
  .Call(C_clab_r_sarml, y, x, w, model, durbin, method, int)
}

#' @export
print.sarml <- function(x, ...) {
  f <- x[["fit"]]
  cat("sarml: ", x[["model"]], if (x[["durbin"]]) " Durbin", " model by ",
      x[["method"]], ", ", f[["n"]], " units", sep = "")
  if (x[["method"]] == "ml") {
    cat(", log-likelihood ", format(signif(f[["loglik"]], 6)), " (OLS ",
        format(signif(f[["llols"]], 6)), ")", sep = "")
  }
  cat("\n")
  print(x[["coef"]], row.names = FALSE)
  invisible(x)
}

#' Direct, indirect and total effects of a spatial lag model
#'
#' In a model with a spatial lag of the outcome, `y = rho W y + X beta + W X
#' theta + e`, a change in covariate `k` at every unit moves the outcomes
#' by `S_k = (I - rho W)^-1 (beta_k I + theta_k W)`: the mean of its
#' diagonal is the direct effect, the mean row sum the total, and their
#' difference the indirect effect, the spillover to other units (LeSage
#' and Pace 2009). An analysis that reads `beta_k` alone can report "no
#' effect" where the effect is all spillover. With row-standardised
#' weights the total is `(beta_k + theta_k) / (1 - rho)`.
#' Standard errors are by the delta method in `rho`,
#' `beta_k` and `theta_k`, with the exact derivatives `d (I - rho W)^-1 /
#' d rho = (I - rho W)^-1 W (I - rho W)^-1`. The effects are those of
#' `spatialreg::impacts()`, which gives simulated rather than delta-method
#' standard errors.
#'
#' @param fit A [sarml()] fit of a lag or combined model (`"lag"` or
#'   `"sac"`, with or without `durbin`); or `NULL`, with `w`, `rho`, `b`
#'   (and `th`) given.
#' @inheritParams gmoran
#' @param rho The spatial lag coefficient.
#' @param b The coefficients of the covariates (named, or `x1`, `x2`, ...).
#' @param th `NULL`, or the coefficients of their lags, as many as `b`.
#' @return A data frame with one row per covariate: `term`, `direct`,
#'   `indirect`, `total` and their standard errors `sed`, `sei`, `set` (`NaN`
#'   without a fit).
#' @references LeSage, J. and Pace, R. K. (2009). *Introduction to Spatial
#'   Econometrics*. Chapman and Hall/CRC, Boca Raton.
#'   \doi{10.1201/9781420064254}
#' @examples
#' set.seed(2)
#' p <- cbind(x = runif(40), y = runif(40))
#' w <- spwt(p, k = 4)
#' spill(w = w, rho = 0.4, b = c(gent = 0.1), th = 0.3)
#' @export
spill <- function(fit = NULL, w = NULL, rho = NULL, b = NULL, th = NULL) {
  if (is.null(fit)) {
    if (is.null(w) || is.null(rho) || is.null(b)) {
      stop("spill: give a sarml() fit, or `w`, `rho` and `b`.", call. = FALSE)
    }
    nm <- names(b)
    if (is.null(nm)) nm <- paste0("x", seq_along(b))
    return(.Call(C_clab_r_spill, w, rho, b, th, NULL, nm))
  }
  if (!inherits(fit, "sarml") || !fit[["model"]] %in% c("lag", "sac")) {
    stop("spill: `fit` must be a sarml() fit of a lag or combined model; ",
         "in an error model the coefficients are the effects.", call. = FALSE)
  }
  tm <- fit[["coef"]][["term"]]
  v <- fit[["vcov"]]
  nl <- fit[["nlag"]]
  lg <- if (nl > 0L) utils::tail(rownames(v)[seq_len(nrow(v) - 2L)], nl)
  xs <- setdiff(rownames(v)[seq_len(nrow(v) - 2L)], c("(Intercept)", lg))
  if (nl > 0L) xs <- sub("^lag[.]", "", lg)
  est <- stats::setNames(fit[["coef"]][["est"]], tm)
  id <- c("rho", xs, lg)
  .Call(C_clab_r_spill, fit[["w"]], est[["rho"]], unname(est[xs]),
        if (nl > 0L) unname(est[lg]), v[id, id], xs)
}
