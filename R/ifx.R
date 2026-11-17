#' Generalized synthetic control with interactive fixed effects
#'
#' The effect of a treatment adopted by some units of a balanced panel at
#' possibly different times, estimated as the generalized synthetic
#' control of Xu (2017): a model of the untreated outcome fitted to the
#' units never treated, then carried to each treated unit through its
#' periods before adoption.
#'
#' The untreated outcome of unit i at time t is
#' `y_it = x_it'b + mu + alpha_i + xi_t + lambda_i'f_t + e_it`, the
#' interactive fixed effects model of Bai (2009): unit and time effects
#' plus `r` common factors `f_t` with unit loadings `lambda_i`, which let
#' units follow different trends (a difference in differences assumes
#' `r = 0`). On the never-treated units the estimator of Bai (2009)
#' alternates the factors, `sqrt(T)` times the first `r` eigenvectors of
#' `W W'` for the doubly demeaned residuals `W` net of covariates
#' (normalised `F'F / T = I`), and `b`, by least squares of the doubly
#' demeaned outcome less the common component on the doubly demeaned
#' covariates, until `b` changes by less than `tol`. Each treated unit's
#' `alpha_i` and `lambda_i` come from least squares of its residuals on
#' `(1, f_t)` over its periods before adoption, which must number more
#' than `r + 1`; its counterfactual `y0` is the model's value at every
#' period, and its effect at a period is `y - y0`. `att` averages the
#' effects of the treated periods; `eff` averages them by periods since
#' adoption (event time `e`, 0 the first treated period; the effects
#' before adoption are the fit's in-sample gaps).
#'
#' `r = NULL` chooses the number of factors from 0 to `rmax` by the
#' cross-validation of Xu (2017, sec. 3.1, algorithm 1): each treated
#' unit's pre-adoption periods are left out one at a time, the loadings
#' fitted on the rest and the left-out outcome predicted; the `r` with the
#' smallest mean squared prediction error is kept (`cv`). The search stops
#' at the first `r` some treated unit has too few periods for.
#'
#' The standard errors are those of a bootstrap of units, the "simple
#' nonparametric bootstrap" Xu (2017, sec. 3.2, p. 63) finds valid when
#' the treated units are many: each draw takes as many units as the panel
#' has, with replacement, and refits with the chosen `r`; draws with no
#' treated unit or with no more control units than `r` fail and are left
#' out. With few treated units Xu (2017)
#' prefers a parametric bootstrap that resamples residuals, which is not
#' built here. Rows of a unit drawn twice enter twice.
#'
#' Barber (2026) estimates the effect of the New York City congestion
#' charge on fine particles, traffic speeds and arrests with the
#' generalized synthetic control, and the spatial shift of arrests with a
#' difference in differences (see [didfe()]).
#'
#' @param data A data frame, one row per unit and time (balanced: every
#'   unit at every time once).
#' @param y The name of the outcome column (numbers).
#' @param d The name of the treatment column: 0 or 1 (or `FALSE`, `TRUE`),
#'   staying 1 once it turns 1; never-treated units are the controls.
#' @param unit,time The names of the unit column (labels) and the time
#'   column (numbers or dates; ordered by value).
#' @param x `NULL` or the names of covariate columns (numbers).
#' @param r `NULL` (cross-validation) or the number of factors, a whole
#'   number from 0.
#' @param rmax The largest number of factors cross-validation tries.
#' @param se `"boot"` (the bootstrap of units) or `"none"`.
#' @param tol The convergence tolerance of `b`.
#' @param maxit The most alternations of Bai's estimator.
#' @param object An object of class `ifesc`.
#' @inheritParams ppmlcf
#' @inheritParams ppmlfe
#' @return An object of class `ifesc`, a list: `att`, a one-row data frame
#'   (`est`, `se`, `stat`, `p`, `lo`, `hi`); `eff`, one row per event time
#'   `e` (`n` treated units observed then, and the same columns); `r`; `cv`,
#'   a data frame of `r` and `mspe` (`NULL` when `r` was given); `beta`,
#'   the covariates' coefficients; `cf`, one row per treated unit and time
#'   (`unit`, `time`, `e`, `y`, `y0`, `eff`); `factors`, `T x r`; `draws`,
#'   the bootstrap's `att` then `eff` per draw, a row each (`NULL` for
#'   `se = "none"`);
#'   `ntr`, `nco`, `nt`; `iter`, `conv`; `nboot`, `nfail`, `level`.
#'   `summary()` returns `eff`.
#' @references
#'   Bai, J. (2009). Panel data models with interactive fixed effects.
#'   *Econometrica*, 77(4), 1229-1279. \doi{10.3982/ECTA6135}
#'
#'   Xu, Y. (2017). Generalized synthetic control method: causal inference
#'   with interactive fixed effects models. *Political Analysis*, 25(1),
#'   57-76. \doi{10.1017/pan.2016.2}
#'
#'   Barber, T. J. (2026). *The New York City congestion charge: air
#'   quality, traffic, and crime*. MS thesis, Cornell University.
#' @examples
#' set.seed(1)
#' p <- expand.grid(id = 1:30, yr = 1:12)
#' f <- cumsum(rnorm(12))
#' l <- rnorm(30)
#' p$tr <- as.integer(p$id <= 6 & p$yr >= 9)
#' p$y <- l[p$id] * f[p$yr] + 2 * p$tr + rnorm(nrow(p), 0, 0.3)
#' r <- ifesc(p, "y", "tr", "id", "yr", r = 1, nboot = 50, seed = 1)
#' r
#' @export
ifesc <- function(
  data,
  y,
  d,
  unit,
  time,
  x = NULL,
  r = NULL,
  rmax = 5,
  se = c("boot", "none"),
  nboot = 200,
  seed = NULL,
  level = 0.95,
  tol = 1e-8,
  maxit = 1000
) {
  f <- "ifesc"
  se <- .ckmarg("se", f)
  dd <- .axdat(data, f)
  pn <- .ifpan(dd, y, d, unit, time, x, f)
  if (!is.null(r)) {
    r <- .cknum(r, "r", f, min = 0, max = .Machine$integer.max, whole = TRUE)
  }
  rmax <- .cknum(rmax, "rmax", f, min = 0, max = 50, whole = TRUE)
  level <- .axlevel(level, f)
  tol <- .cknum(tol, "tol", f, positive = TRUE)
  maxit <- .cknum(maxit, "maxit", f, min = 1, max = .Machine$integer.max,
                  whole = TRUE)
  bm <- NULL
  if (se == "boot") {
    nboot <- .cknum(nboot, "nboot", f, min = 2, max = 1e7, whole = TRUE)
    if (!is.null(seed)) seed <- .axseed(seed, f)
    bm <- .ifdraws(ncol(pn$y), nboot, seed)
  }
  k <- .Call(C_clab_r_ifesc, pn$y, pn$x, pn$a, if (is.null(r)) -1L else
               as.integer(r), as.integer(rmax), as.integer(maxit), tol, bm)
  .ifout(k, pn, x, r, se, if (se == "boot") nboot else 0L, level, f)
}

#' @rdname ifesc
#' @export
print.ifesc <- function(x, ...) {
  cat("ifesc: generalized synthetic control, ", x[["ntr"]], " treated and ",
      x[["nco"]], " control units, ", x[["nt"]], " times, ", x[["r"]],
      " factor(s)", if (is.null(x[["cv"]])) "" else " (cross-validated)",
      "\n", sep = "")
  print(x[["att"]], row.names = FALSE, digits = 4)
  if (!x[["conv"]]) cat("Bai's estimator did not converge.\n")
  invisible(x)
}

#' @rdname ifesc
#' @export
summary.ifesc <- function(object, ...) object[["eff"]]

# The balanced panel: y (T x N), x (T x N x p), each unit's adoption
# period a (0-based; T for never treated), the units and the times.
.ifpan <- function(d, y, tr, unit, time, x, f) {
  .axcols(d, y, "y", f, n = 1L)
  .axcols(d, tr, "d", f, n = 1L)
  .axcols(d, unit, "unit", f, n = 1L)
  .axcols(d, time, "time", f, n = 1L)
  x <- .axcols(d, x, "x", f, null = TRUE)
  yv <- .axmat(d, y, f)[, 1L]
  xm <- .axmat(d, x, f)
  dv <- .axmat(d, tr, f)[, 1L]
  if (!all(dv %in% c(0, 1))) .ckstop(f, "column `%s` must hold 0 or 1", tr)
  tv <- d[[time]]
  if (!is.numeric(tv) && !inherits(tv, "Date")) {
    .ckstop(f, "column `%s` must hold numbers or dates", time)
  }
  uu <- unique(as.character(d[[unit]]))
  tt <- sort(unique(tv))
  iu <- match(as.character(d[[unit]]), uu)
  it <- match(tv, tt)
  nt <- length(tt)
  if (nrow(d) != length(uu) * nt || anyDuplicated(cbind(iu, it))) {
    .ckstop(f, "the panel must be balanced: each unit once at each time")
  }
  o <- order(iu, it)
  dm <- matrix(dv[o], nt)
  a <- apply(dm, 2L, function(v) if (any(v == 1)) which.max(v) - 1L else nt)
  if (any(apply(dm, 2L, is.unsorted))) {
    .ckstop(f, "the treatment must stay 1 once it turns 1")
  }
  if (any(a == 0L)) {
    .ckstop(f, "a treated unit is treated at the first time: it has no %s",
            "periods before adoption")
  }
  .ifck(a, nt, f)
  list(y = matrix(yv[o], nt), x = array(xm[o, ], c(nt, length(uu), ncol(xm))),
       a = as.integer(a), units = uu, times = tt)
}

.ifck <- function(a, nt, f) {
  if (all(a == nt)) .ckstop(f, "no unit is treated")
  if (sum(a == nt) < 2L) {
    .ckstop(f, "the panel needs two or more never-treated units")
  }
}

# The bootstrap counts of the units, N x nboot.
.ifdraws <- function(n, nboot, seed) {
  acc <- new.env(parent = emptyenv())
  acc$bm <- matrix(0L, n, nboot)
  .axdraws(n, nboot, seed, function(bm, i) {
    acc$bm[, i + seq_len(ncol(bm))] <- bm
  })
  storage.mode(acc$bm) <- "integer"
  acc$bm
}

# The result of ifesc() from the kernel's list k.
.ifout <- function(k, pn, x, r, se, nboot, level, f) {
  st <- k$st
  ne <- st[8L]
  e <- st[7L] + seq_len(ne) - 1L
  tr <- which(pn$a < length(pn$times))
  nt <- length(pn$times)
  ee <- outer(seq_len(nt) - 1L, pn$a[tr], "-")
  dr <- NULL
  sev <- rep(NA_real_, ne + 1L)
  if (se == "boot") {
    dr <- k$bef
    nf <- st[6L]
    if (nf > 0L) {
      warning(f, ": ", nf, " of ", nboot, " bootstrap draws could not be ",
              "fitted and are left out.", call. = FALSE)
    }
    if (nboot - nf < 2L) {
      .ckstop(f, "fewer than two bootstrap draws could be fitted")
    }
    sev <- apply(dr, 2L, stats::sd, na.rm = TRUE)
  }
  tb <- .axtab(c("att", e), k$ef, sev, level)
  tb$term <- NULL
  eff <- cbind(data.frame(e = e, n = as.vector(table(factor(ee, e)))),
               tb[-1L, , drop = FALSE])
  rownames(eff) <- NULL
  y0 <- k$y0[, tr, drop = FALSE]
  yy <- pn$y[, tr, drop = FALSE]
  cf <- data.frame(unit = rep(pn$units[tr], each = nt),
                   time = rep(pn$times, length(tr)), e = as.vector(ee),
                   y = as.vector(yy), y0 = as.vector(y0),
                   eff = as.vector(yy - y0))
  structure(
    list(att = tb[1L, , drop = FALSE], eff = eff, r = as.integer(st[1L]),
         cv = if (is.null(r)) .ifcv(k$cv) else NULL,
         beta = stats::setNames(k$beta, x), cf = cf,
         factors = matrix(k$fac, nt), draws = dr, ntr = as.integer(st[4L]),
         nco = as.integer(st[5L]), nt = nt, iter = as.integer(st[2L]),
         conv = st[3L] == 1, nboot = as.integer(nboot),
         nfail = if (se == "boot") as.integer(st[6L]) else 0L,
         level = level),
    class = "ifesc"
  )
}

.ifcv <- function(cv) {
  k <- !is.na(cv)
  data.frame(r = which(k) - 1L, mspe = cv[k])
}
