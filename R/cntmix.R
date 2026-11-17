#' Count regression with a random intercept per cluster
#'
#' Fits the log-linear model of a count with a normal random intercept
#' for each cluster, \eqn{\log E[y_{ij} | u_j] = x_{ij}'\beta + u_j}{log
#' E[y_ij | u_j] = x_ij' beta + u_j}, \eqn{u_j \sim N(0, \sigma^2)}{u_j ~
#' N(0, sigma^2)}, with a Poisson or a negative binomial (variance `mu +
#' mu^2 / theta`) response, by maximum likelihood with the Laplace
#' approximation of each cluster's integral (Breslow and Clayton 1993), as
#' 'glmmTMB' and 'lme4' fit it. It is the main model of Ruhela (2026):
#' placement volatility on alert complexity with a random intercept per
#' regional cluster.
#'
#' The algorithm: for given `beta`, `log sigma` (and `log theta`), each
#' cluster's mode of the random intercept is found by Newton's method on
#' its one-dimensional log-density, all clusters at once (`O(n)` per step,
#' quadratic convergence); the Laplace approximation of the marginal
#' log-likelihood is then a sum over clusters. The outer parameters are
#' maximised by BFGS from the Poisson fit without random effects, and
#' their standard errors are from the inverse of the numerically
#' differentiated Hessian at the maximum.
#'
#' When the negative binomial's `theta` runs off to the Poisson limit (above
#' `1e6`), its likelihood is flat in that direction and the Hessian is
#' singular; the Poisson model is then reported, as it is the limit, with
#' `theta` noted (`fam` says which). Fewer than two clusters, or clusters of
#' one row each, are refused: the variance of the intercepts is then not
#' identified.
#'
#' A coefficient that cannot be estimated is `NA` with its standard error,
#' test and interval, and `notes` says why, never a standard error of 0:
#' a column aliased with the others (a constant, or a combination of other
#' columns, found by the QR decomposition of the design before the fit,
#' which drops it as [stats::lm()] does); a 0/1 column all of whose counts
#' on one side are 0 (separation: its rate ratio is 0 or infinite, and so
#' is the intercept's when the zeros are where the column is 0); and an
#' estimate whose variance from the Hessian is not a positive number.
#'
#' @param y The counts: whole numbers from 0.
#' @param x The covariates: a data frame (factors and text enter as
#'   indicators of each level but the first) or a numeric matrix, one row
#'   per count; an intercept is added.
#' @param grp The cluster of each count: labels, at least two clusters, one
#'   of them with two rows or more.
#' @param fam `"nb"` or `"poisson"`.
#' @param w `NULL` or a weight from 0 for each count.
#' @inheritParams cntreg
#' @return An object of class `cntmix`, a list: `coef`, a data frame
#'   (`term`, `est`, `se`, `z`, `p`, `irr`, `lo`, `hi`; `NA` for a term
#'   that cannot be estimated); `sigma`, the standard deviation of the
#'   intercepts; `theta` (`NA` for Poisson); `loglik`, `aic`; `fam`, the
#'   family reported; `n`, `ngrp`; `conv`, the optimiser's code (0 when it
#'   converged); `level`; `notes`, why a term has no estimate (empty when
#'   every term has one).
#' @references Breslow, N. E. and Clayton, D. G. (1993). Approximate
#'   inference in generalized linear mixed models. *Journal of the
#'   American Statistical Association*, 88(421), 9-25.
#'   \doi{10.1080/01621459.1993.10594284}
#'
#'   Brooks, M. E., Kristensen, K., van Benthem, K. J., Magnusson, A.,
#'   Berg, C. W., Nielsen, A., Skaug, H. J., Maechler, M. and Bolker, B.
#'   M. (2017). glmmTMB balances speed and flexibility among packages for
#'   zero-inflated generalized linear mixed modeling. *The R Journal*,
#'   9(2), 378-400. \doi{10.32614/RJ-2017-066}
#' @examples
#' set.seed(2)
#' g <- rep(1:12, each = 40)
#' d <- rbinom(480, 1, 0.3)
#' y <- rpois(480, exp(-1 + 0.4 * d + rnorm(12, 0, 0.5)[g]))
#' cntmix(y, data.frame(d = d), g, fam = "poisson")
#' @export
cntmix <- function(
  y,
  x,
  grp,
  fam = c("nb", "poisson"),
  w = NULL,
  level = 0.95
) {
  fn <- "cntmix"
  fam <- .ckmarg("fam", fn)
  y <- .cknum(y, "y", fn, scalar = FALSE, min = 0, whole = TRUE)
  X <- .cmx(x, length(y), fn)
  w <- if (is.null(w)) {
    rep(1, length(y))
  } else {
    .cknum(w, "w", fn, scalar = FALSE, min = 0)
  }
  level <- .xpopen(level, "level", fn)
  g <- .cmgrp(grp, length(y), fn)
  if (length(w) != length(y)) {
    .ckstop(fn, "`w` must have one value per count")
  }
  if (all(y == 0)) {
    .ckstop(fn, "`y` is 0 throughout: there is no rate to model")
  }
  e <- .cmest(X, y, w)
  Xk <- X[, e$keep, drop = FALSE]
  b0 <- stats::glm.fit(Xk, y, w, family = stats::poisson())$coefficients
  b0[is.na(b0)] <- 0
  fit <- .cmfit(y, Xk, g, w, fam, b0)
  if (fam == "nb" && exp(fit$par[length(fit$par)]) > 1e6) {
    th <- exp(fit$par[length(fit$par)])
    fit <- .cmfit(y, Xk, g, w, "poisson", fit$par[seq_len(ncol(Xk))])
    fit$fam <- sprintf(
      "poisson (negative binomial theta %.3g: the Poisson limit)",
      th
    )
  }
  .cmout(fit, X, e, y, g, level)
}

#' @export
print.cntmix <- function(x, ...) {
  cat(
    "cntmix: ",
    x$fam,
    " counts with a random intercept per cluster\n",
    x$n,
    " counts in ",
    x$ngrp,
    " clusters; sigma ",
    format(x$sigma, digits = 4),
    if (!is.na(x$theta)) paste0(", theta ", format(x$theta, digits = 4)),
    "; log-likelihood ",
    format(x$loglik, digits = 7),
    if (x$conv != 0) " (the optimiser did not converge)",
    "\n",
    sep = ""
  )
  print(x$coef, row.names = FALSE, digits = 4)
  if (length(x$notes)) {
    cat(paste0("note: ", x$notes, "\n"), sep = "")
  }
  invisible(x)
}

# The design matrix with an intercept, its columns named.
.cmx <- function(x, n, fn) {
  if (is.matrix(x) && is.numeric(x)) {
    x <- as.data.frame(x)
  }
  if (!is.data.frame(x) || nrow(x) != n || anyNA(x)) {
    .ckstop(
      fn,
      "`x` must be a data frame or numeric matrix, one complete row per count"
    )
  }
  x[] <- lapply(x, function(v) {
    if (is.character(v) || is.logical(v)) factor(v) else v
  })
  X <- stats::model.matrix(
    ~.,
    x,
    contrasts.arg = lapply(
      Filter(is.factor, x),
      function(v) "contr.treatment"
    )
  )
  X[, colSums(abs(X)) > 0 | colnames(X) == "(Intercept)", drop = FALSE]
}

# The clusters as integers 1..G, refusing what cannot identify sigma.
.cmgrp <- function(grp, n, fn) {
  if (length(grp) != n || anyNA(grp)) {
    .ckstop(fn, "`grp` must give a cluster for every count")
  }
  g <- as.integer(factor(grp))
  tb <- tabulate(g)
  if (length(tb) < 2L || all(tb < 2L)) {
    .ckstop(fn, "`grp` needs two clusters or more, one with two rows or more")
  }
  g
}

# Log-density of each row and its first two derivatives in the linear
# predictor eta.
.cmll <- function(y, eta, th) {
  mu <- exp(eta)
  if (is.na(th)) {
    return(list(l = y * eta - mu - lgamma(y + 1), d1 = y - mu, d2 = -mu))
  }
  list(
    l = lgamma(y + th) -
      lgamma(th) -
      lgamma(y + 1) +
      th * log(th / (th + mu)) +
      y * (eta - log(th + mu)),
    d1 = th * (y - mu) / (th + mu),
    d2 = -th * mu * (y + th) / (th + mu)^2
  )
}

# Laplace approximation of the marginal log-likelihood at par = (beta,
# log sigma[, log theta]). Each cluster's log-density in its intercept is
# concave, so Newton's method from 0 with steps capped at 1 converges; the
# fixed start keeps the value a function of par alone, as the finite
# differences of the optimiser need. A point where it overflows is given
# a very low value.
.cmlap <- function(par, y, X, g, w, fam) {
  p <- ncol(X)
  s2 <- exp(2 * par[p + 1L])
  th <- if (fam == "nb") exp(par[p + 2L]) else NA
  eta0 <- drop(X %*% par[seq_len(p)])
  u <- rep(0, max(g))
  for (it in 1:200) {
    d <- .cmll(y, eta0 + u[g], th)
    gr <- drop(rowsum(w * d$d1, g)) - u / s2
    h <- drop(rowsum(w * d$d2, g)) - 1 / s2
    u <- u - pmax(pmin(gr / h, 1), -1)
    if (!all(is.finite(u)) || max(abs(gr)) < 1e-10) break
  }
  d <- .cmll(y, eta0 + u[g], th)
  h <- drop(rowsum(w * d$d2, g)) - 1 / s2
  ll <- sum(w * d$l) + sum(-u^2 / (2 * s2) - 0.5 * log(s2) - 0.5 * log(-h))
  if (is.finite(ll)) ll else -1e300
}

# BFGS on the outer parameters, then the Hessian.
.cmfit <- function(y, X, g, w, fam, b0) {
  par <- c(b0, 0, if (fam == "nb") 0)
  f <- function(p) -.cmlap(p, y, X, g, w, fam)
  o <- stats::optim(
    par,
    f,
    method = "BFGS",
    control = list(maxit = 1000, reltol = 1e-12)
  )
  hs <- stats::optimHess(o$par, f)
  list(
    par = o$par,
    ll = -o$value,
    conv = o$convergence,
    hess = hs,
    fam = fam,
    nb = fam == "nb"
  )
}

# The columns of the design that can be estimated (e$keep, by the QR
# decomposition with lm()'s tolerance) and the terms that have no finite
# estimate (e$na: those aliased, and those a 0/1 column separates), with
# a note on each (e$notes). Rows of weight 0 tell nothing.
.cmest <- function(X, y, w) {
  q <- qr(X)
  keep <- sort(q$pivot[seq_len(q$rank)])
  na <- colnames(X)[-keep]
  notes <- if (length(na)) {
    sprintf(
      "%s: aliased with the other columns (constant, or collinear), dropped",
      toString(na)
    )
  }
  on <- w > 0
  for (j in keep[colnames(X)[keep] != "(Intercept)"]) {
    v <- X[on, j]
    if (!all(v %in% c(0, 1)) || all(v == v[1L])) next
    z <- c(all(y[on][v == 0] == 0), all(y[on][v == 1] == 0))
    if (any(z)) {
      na <- c(na, colnames(X)[j], if (z[1L]) "(Intercept)")
      notes <- c(notes, sprintf(
        "%s: every count where it is %d is 0, so its rate ratio is %s",
        colnames(X)[j], if (z[1L]) 0L else 1L,
        if (z[1L]) "infinite" else "0"
      ))
    }
  }
  list(keep = keep, na = unique(na), notes = notes)
}

# The result object: the coefficients of the columns e$keep from the fit,
# NA for the terms e$na and for a variance that is not a positive number.
.cmout <- function(fit, X, e, y, g, level) {
  p <- length(e$keep)
  v <- tryCatch(solve(fit$hess), error = function(err) {
    matrix(NA_real_, nrow(fit$hess), ncol(fit$hess))
  })
  d <- diag(v)[seq_len(p)]
  est <- se <- rep(NA_real_, ncol(X))
  est[e$keep] <- fit$par[seq_len(p)]
  se[e$keep[is.finite(d) & d > 0]] <- sqrt(d[is.finite(d) & d > 0])
  bad <- colnames(X) %in% e$na
  est[bad] <- se[bad] <- NA_real_
  hs <- colnames(X)[!bad & is.na(se)]
  notes <- c(e$notes, if (length(hs)) {
    sprintf("%s: the Hessian gives no positive variance, so no standard error",
            toString(hs))
  })
  z <- stats::qnorm(1 - (1 - level) / 2)
  cf <- data.frame(
    term = colnames(X),
    est = est,
    se = se,
    z = est / se,
    p = 2 * stats::pnorm(-abs(est / se)),
    irr = exp(est),
    lo = exp(est - z * se),
    hi = exp(est + z * se)
  )
  k <- length(fit$par)
  structure(
    list(
      coef = cf,
      sigma = exp(fit$par[p + 1L]),
      theta = if (fit$nb) exp(fit$par[k]) else NA_real_,
      loglik = fit$ll,
      aic = 2 * k - 2 * fit$ll,
      fam = fit$fam,
      n = length(y),
      ngrp = max(g),
      conv = fit$conv,
      level = level,
      notes = as.character(notes)
    ),
    class = "cntmix"
  )
}
