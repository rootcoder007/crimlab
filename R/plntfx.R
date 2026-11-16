#' Does the analysis recover a planted effect?
#'
#' Plants effects of known size in the real data and checks that the
#' analysis finds them: the outcome is shuffled across rows, which removes
#' whatever association it had with the treatment and the covariates
#' while keeping its distribution, an effect `eff` is added to the treated
#' (`y* = shuffled y + eff d`), and the analysis is run again, `nrep`
#' times for each effect. An analysis that is biased, or whose intervals
#' are too narrow, shows it here on data where the truth is known, a
#' positive control in the sense of Lipsitch, Tchetgen Tchetgen and Cohen
#' (2010).
#'
#' The analysis is the coefficient of `d` in the least squares fit of the
#' outcome on an intercept, `x` and `d`, with a classical or HC1 standard
#' error and Student's t interval at `level`, as [falsctl()] fits it. For
#' each effect: `mean` and `sd` of the estimates, `bias = mean - eff`,
#' `rmse` (about `eff`), `cover` (the share of intervals holding `eff`),
#' `rej` (the share rejecting 0 at `1 - level`), and `ok`, whether `|bias|`
#' is at most `tol`. By default `tol` is three Monte Carlo standard errors
#' of the mean, `3 sd / sqrt(nrep)`. One shuffle serves every effect, and
#' the shuffles come from the package's generator.
#'
#' @inheritParams cemtch
#' @inheritParams rndalc
#' @inheritParams ifesc
#' @param d The name of the treatment column (numbers; it must vary).
#' @param x `NULL` or the names of the covariate columns (numbers); for
#'   `print()`, an object of class `plntfx`.
#' @param eff The effects to plant, in the outcome's units: finite numbers.
#' @param nrep The shuffles per effect: a whole number from 2 to 1000000.
#' @param se `"hc1"` or `"iid"`.
#' @param level The level of the intervals: one number strictly between 0
#'   and 1.
#' @param tol `NULL` (three Monte Carlo standard errors) or the largest
#'   `|bias|` accepted: one number from 0.
#' @return An object of class `plntfx`, a list: `tab`, a data frame with
#'   one row per effect (`eff`, `mean`, `bias`, `sd`, `rmse`, `cover`,
#'   `rej`, `tol`, `ok`); and `nrep`, `se`, `level` and `seed`.
#' @references
#'   Lipsitch, M., Tchetgen Tchetgen, E. and Cohen, T. (2010). Negative
#'   controls: a tool for detecting confounding and bias in observational
#'   studies. *Epidemiology*, 21(3), 383-388.
#'   \doi{10.1097/EDE.0b013e3181d61eeb}
#'
#'   Morris, T. P., White, I. R. and Crowther, M. J. (2019). Using
#'   simulation studies to evaluate statistical methods. *Statistics in
#'   Medicine*, 38(11), 2074-2102. \doi{10.1002/sim.8086}
#' @examples
#' set.seed(6)
#' d <- data.frame(age = runif(200, 18, 60), t = rbinom(200, 1, 0.4))
#' d$y <- 0.02 * d$age + rnorm(200)
#' plntfx(d, "y", "t", x = "age", eff = c(0, 0.25, 0.5), nrep = 200,
#'        seed = 1)
#' @export
plntfx <- function(
  data,
  y,
  d,
  x = NULL,
  eff,
  nrep = 200,
  se = c("hc1", "iid"),
  level = 0.95,
  tol = NULL,
  seed = NULL
) {
  f <- "plntfx"
  dd <- .axdat(data, f)
  .axcols(dd, y, "y", f, n = 1L)
  .axcols(dd, d, "d", f, n = 1L)
  .axcols(dd, x, "x", f, null = TRUE)
  if (d == y || any(c(y, d) %in% x)) {
    .ckstop(f, "`y`, `d` and `x` must name different columns")
  }
  se <- .ckmarg("se", f)
  yv <- .axmat(dd, y, f)[, 1L]
  dv <- .axmat(dd, d, f)[, 1L]
  xm <- .axmat(dd, x, f)
  if (nrow(dd) <= length(x) + 2L) {
    .ckstop(f, "`data` needs more rows than covariates plus two")
  }
  if (all(dv == dv[1L])) .ckstop(f, "column `%s` must vary", d)
  eff <- .cknum(eff, "eff", f, scalar = FALSE)
  nrep <- .cknum(nrep, "nrep", f, min = 2, max = 1e6, whole = TRUE)
  level <- .xpopen(level, "level", f)
  if (!is.null(tol)) tol <- .cknum(tol, "tol", f, min = 0)
  seed <- .xpseed(seed, f)
  o <- .Call(C_clab_r_plntfx, xm, yv, dv, eff, as.integer(nrep), se == "hc1",
             level, seed)
  m <- length(eff)
  mn <- o[seq_len(m)]
  sdv <- o[m + seq_len(m)]
  tl <- if (is.null(tol)) 3 * sdv / sqrt(nrep) else rep(tol, m)
  tab <- data.frame(eff = eff, mean = mn, bias = mn - eff, sd = sdv,
                    rmse = o[2L * m + seq_len(m)],
                    cover = o[3L * m + seq_len(m)],
                    rej = o[4L * m + seq_len(m)], tol = tl,
                    ok = abs(mn - eff) <= tl)
  structure(list(tab = tab, nrep = nrep, se = se, level = level, seed = seed),
            class = "plntfx")
}

#' @rdname plntfx
#' @inheritParams cemtch
#' @export
print.plntfx <- function(x, ...) {
  t <- x[["tab"]]
  cat("plntfx: ", sum(t$ok), " of ", nrow(t),
      " planted effect(s) recovered within tolerance over ",
      format(x[["nrep"]], scientific = FALSE), " shuffles\n", sep = "")
  print(t, row.names = FALSE, digits = 4)
  invisible(x)
}
