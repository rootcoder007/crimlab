#' Transport a trial's effect to a target population
#'
#' The average treatment effect of a randomised trial carried over to a
#' target population whose covariates differ from the trial's, by the
#' outcome model and by inverse odds weighting (Dahabreh et al. 2020;
#' Westreich et al. 2017). The data hold the trial's units (`s = 1`), with
#' outcome and treatment, and a sample of the target (`s = 0`), with the
#' covariates only; `y` and `d` are not read outside the trial and may be
#' `NA` there.
#'
#' The target is the rows outside the trial (`target = "transport"`) or
#' all rows, trial included (`"generalize"`). The outcome model fits the
#' outcome on an intercept and `x` by least squares in each arm of the
#' trial and averages the difference of the two predictions over the
#' target (`om`). Inverse odds weighting gives each trial unit the
#' weight of [gensw()] (`(1 - ps) / ps` to transport, `1 / ps` to
#' generalise, the sampling score `ps` trimmed to `[trim, 1 - trim]`) and
#' takes the difference of the arms' weighted means, the weights
#' normalised within each arm (`iow`; Tipton and Hartman 2023, eq. 3.10).
#' The first is right when the outcome model is, the second when the
#' sampling model is. Standard errors come from `nboot` bootstrap draws of
#' the rows within the trial and within the rest, a draw whose fit fails
#' left out (`nfail`), with normal intervals at `level`.
#'
#' @inheritParams gensw
#' @inheritParams plntfx
#' @param y The name of the outcome column (numbers; read in the trial).
#' @param d The name of the treatment column: 0 or 1 in the trial, both
#'   present.
#' @param x The names of the covariate columns (numbers), complete in
#'   every row; for `print()`, an object of class `trnsp`.
#' @param target `"transport"` or `"generalize"`.
#' @param trim The bound of the sampling scores, from 0 to below 0.5.
#' @param nboot The bootstrap draws: a whole number from 0 (no standard
#'   errors) to 100000.
#' @param seed `NULL` (drawn from R's generator, so [set.seed()] fixes it)
#'   or a whole number from 0 to 2147483647; not used without draws.
#' @return An object of class `trnsp`, a list: `tab`, a data frame with
#'   rows `om` and `iow` (`method`, `est`, `se`, `lo`, `hi`); `ps` and `w`,
#'   the sampling scores and weights of every row; `draws` (a matrix,
#'   `NA` for a draw left out); `nfail`; and `target`, `nboot`, `level`,
#'   `seed`.
#' @references
#'   Dahabreh, I. J., Robertson, S. E., Steingrimsson, J. A., Stuart, E. A.
#'   and Hernan, M. A. (2020). Extending inferences from a randomized trial
#'   to a new target population. *Statistics in Medicine*, 39(14),
#'   1999-2014. \doi{10.1002/sim.8426}
#'
#'   Westreich, D., Edwards, J. K., Lesko, C. R., Stuart, E. and Cole, S.
#'   R. (2017). Transportability of trial results using inverse odds of
#'   sampling weights. *American Journal of Epidemiology*, 186(8),
#'   1010-1014. \doi{10.1093/aje/kwx164}
#'
#'   Tipton, E. and Hartman, E. (2023). Generalizability and
#'   transportability. In Zubizarreta, J. R., Stuart, E. A., Small, D. S.
#'   and Rosenbaum, P. R. (eds), *Handbook of Matching and Weighting
#'   Adjustments for Causal Inference*, ch. 3. Chapman and Hall/CRC.
#'   \doi{10.1201/9781003102670}
#' @examples
#' set.seed(9)
#' n <- 600
#' d <- data.frame(age = rnorm(n, 35, 10))
#' d$s <- rbinom(n, 1, plogis(-0.05 * (d$age - 35)))
#' d$t <- ifelse(d$s == 1, rbinom(n, 1, 0.5), NA)
#' d$y <- ifelse(d$s == 1, 0.5 * d$t + 0.02 * d$t * d$age + rnorm(n), NA)
#' trnsp(d, "y", "t", "s", "age", nboot = 100, seed = 1)
#' @export
trnsp <- function(
  data,
  y,
  d,
  s,
  x,
  target = c("transport", "generalize"),
  trim = 0,
  nboot = 200,
  level = 0.95,
  seed = NULL
) {
  f <- "trnsp"
  dd <- .axdat(data, f)
  .axcols(dd, s, "s", f, n = 1L)
  .axcols(dd, x, "x", f)
  .trcol(dd, y, "y", f)
  .trcol(dd, d, "d", f)
  if (anyDuplicated(c(y, d, s, x))) {
    .ckstop(f, "`y`, `d`, `s` and `x` must name different columns")
  }
  target <- .ckmarg("target", f)
  sv <- .axmat(dd, s, f)[, 1L]
  xm <- .axmat(dd, x, f)
  if (!all(sv %in% c(0, 1))) .ckstop(f, "column `%s` must hold 0 or 1", s)
  if (all(sv == 1) || all(sv == 0)) {
    .ckstop(f, "column `%s` must mark units in and out of the trial", s)
  }
  tr <- sv == 1
  yv <- .trval(dd, y, tr, f)
  dv <- .trval(dd, d, tr, f)
  if (!all(dv[tr] %in% c(0, 1))) {
    .ckstop(f, "column `%s` must hold 0 or 1 in the trial", d)
  }
  if (all(dv[tr] == 1) || all(dv[tr] == 0)) {
    .ckstop(f, "column `%s` must hold both arms in the trial", d)
  }
  trim <- .cknum(trim, "trim", f, min = 0)
  if (trim >= 0.5) .ckstop(f, "`trim` must be below 0.5")
  nboot <- .cknum(nboot, "nboot", f, min = 0, max = 1e5, whole = TRUE)
  level <- .xpopen(level, "level", f)
  seed <- if (nboot > 0) .xpseed(seed, f) else NA_real_
  k <- .Call(C_clab_r_trnsp, xm, yv, as.integer(dv), as.integer(sv),
             target == "generalize", trim, as.integer(nboot), seed)
  se <- k$se
  se[is.nan(se)] <- NA_real_
  z <- stats::qnorm(0.5 + level / 2)
  dr <- matrix(k$draws, nboot, 2L, dimnames = list(NULL, c("om", "iow")))
  dr[is.nan(dr)] <- NA_real_
  structure(
    list(tab = data.frame(method = c("om", "iow"), est = k$est, se = se,
                          lo = k$est - z * se, hi = k$est + z * se,
                          stringsAsFactors = FALSE),
         ps = k$ps, w = k$w, draws = dr, nfail = k$nfail, target = target,
         nboot = nboot, level = level, seed = seed),
    class = "trnsp"
  )
}

#' @rdname trnsp
#' @inheritParams cemtch
#' @export
print.trnsp <- function(x, ...) {
  cat("trnsp: the trial's effect ", if (x[["target"]] == "transport")
        "transported to the rows outside it" else
        "generalised to all rows", ", ", x[["nboot"]],
      " bootstrap draws (", x[["nfail"]], " failed)\n", sep = "")
  print(x[["tab"]], row.names = FALSE, digits = 4)
  invisible(x)
}

#' A column name for a variable read in the trial only: the checks of
#' .axcols() but missing values.
#' @noRd
.trcol <- function(d, v, what, f) {
  if (is.null(v)) .ckstop(f, "`%s` must name columns of `data`, got NULL", what)
  if (!is.character(v) || anyNA(v) || !length(v)) {
    .ckstop(f, "`%s` must name columns of `data`", what)
  }
  if (length(v) != 1L) .ckstop(f, "`%s` must name 1 column(s)", what)
  if (!v %in% names(d)) {
    .ckstop(f, "`%s` names columns not in `data`: %s", what, v)
  }
  v
}

#' Its values: numbers, finite in the trial (rows keep), 0 elsewhere.
#' @noRd
.trval <- function(d, v, keep, f) {
  z <- d[[v]]
  if (!is.numeric(z) && !is.logical(z)) {
    .ckstop(f, "column `%s` must be numeric, got %s", v, class(z)[1L])
  }
  z <- as.numeric(z)
  if (!all(is.finite(z[keep]))) {
    .ckstop(f, "column `%s` must be finite in the trial", v)
  }
  z[!keep] <- 0
  z
}
