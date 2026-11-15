# Regression discontinuity: local polynomial fits on either side of a
# cutoff with kernel weights (the fits are clab_ivfe of src/apx.cpp with
# weights and one group), the bandwidth of Imbens and Kalyanaraman (2012),
# the robust interval of Calonico, Cattaneo and Titiunik (2014) with the
# pilot bandwidth equal to the main one, placebo cutoffs, and the density
# test of McCrary (2008) (src/rdd.cpp).

# The kernel weights of u = distance / h: 0 outside (-1, 1).
.rdkw <- function(u, kern) {
  switch(
    kern,
    tri = pmax(1 - abs(u), 0),
    uni = ifelse(abs(u) <= 1, 0.5, 0),
    epa = pmax(0.75 * (1 - u^2), 0)
  )
}

# One local polynomial fit of order p at 0 of the centred running
# variable rv with bandwidth h: the jump in y (sharp) or its ratio to the
# jump in tv (fuzzy, the treatment instrumented by the side), and its
# variance from the rows' scores, heteroskedasticity-robust or clustered
# by cl (codes of all rows, or NULL).
.rdfit <- function(yv, rv, tv, h, kern, p, cl, ssc, fn) {
  w <- .rdkw(rv / h, kern)
  k <- w > 0
  side <- rv >= 0
  if (sum(k & side) < p + 2L || sum(k & !side) < p + 2L) {
    .ckstop(fn, paste0("fewer than %d rows on a side within the bandwidth ",
                       "%g: widen it"), p + 2L, h)
  }
  r <- rv[k]
  dd <- as.numeric(side[k])
  pm <- do.call(cbind, lapply(seq_len(p), function(j) {
    cbind(r^j, dd * r^j)
  }))
  n <- sum(k)
  if (is.null(tv)) {
    xe <- matrix(numeric(0), n, 0L)
    xx <- cbind(dd, pm)
    zi <- matrix(numeric(0), n, 0L)
  } else {
    xe <- matrix(tv[k], n, 1L)
    xx <- pm
    zi <- matrix(dd, n, 1L)
  }
  f0 <- integer(n)
  fit <- .Call(C_clab_r_ivfe, as.numeric(yv[k]), xe, xx, zi, as.numeric(w[k]),
               f0, 1L, NULL, 0L, NULL, NULL, NULL, 0L, NULL, 0L, 1e-12)
  if (is.null(cl)) {
    fac <- if (ssc == "G") n / (n - 1) else 1
    v <- .Call(C_clab_r_vcx, fit$sco, fit$brd, 2L, NULL, 0L, NULL, 0L, NULL,
               0L, fac)
    g <- n
  } else {
    cc <- .axcode(cl[k])
    if (cc$n < 2L) .ckstop(fn, "one cluster within the bandwidth")
    fac <- if (ssc == "G") cc$n / (cc$n - 1) else 1
    v <- .Call(C_clab_r_vcx, fit$sco, fit$brd, 0L, cc$c, cc$n, NULL, 0L,
               NULL, 0L, fac)
    g <- cc$n
  }
  list(est = fit$beta[1L], se = sqrt(max(v[1L, 1L], 0)), n = n,
       nl = sum(!side[k]), nr = sum(side[k]), g = g)
}

#' Regression discontinuity by local polynomials
#'
#' Estimates the jump at a cutoff in a running variable: the outcome's
#' (sharp design) or its ratio to the jump in the treatment (fuzzy
#' design), by local-linear fits on either side with kernel weights,
#' with the bandwidth of Imbens and Kalyanaraman (2012) or one given,
#' conventional and robust intervals, clustered errors and placebo
#' cutoffs. Carrillo, Lopez-Luzuriaga and Malik (2018, equation 6, p. 8)
#' use the sharp design in space: crime counts by distance to a
#' restriction border, a line on each side, within a window.
#'
#' The fit. With `r` the running variable less `cut` and `D = 1(r >= 0)`,
#' weighted least squares of `y` on `1, D, r, D r` with weights `K(r /
#' h)`; the coefficient of `D` is the jump, the difference of the two
#' local-linear intercepts. Fuzzy (`trt` given): two-stage least squares
#' of `y` on the treatment, instrumented by `D`, with `1, r, D r` and the
#' same weights, which is the ratio of the jumps in `y` and in the
#' treatment. Kernels (`kern`, each on `[-1, 1]`): `"tri"` (the default,
#' `1 - |u|`, the edge kernel Imbens and Kalyanaraman recommend, p. 941),
#' `"uni"` (a constant: a plain window, as Carrillo et al.), `"epa"`
#' (`3/4 (1 - u^2)`).
#'
#' The bandwidth (`bw`): `"ik"` (the default), the algorithm of Imbens and
#' Kalyanaraman (2012, section 4.2, steps 1-3, equations 11-17, p.
#' 941-943) for the outcome: a pilot `h1 = 1.84 S_X N^(-1/5)` for the
#' density and the variances at the cutoff, a cubic with a jump for the
#' third derivative, local quadratics within `h2` for the second
#' derivatives, the regularisation terms `2160 sigma^2 / (N2 h2^4)`, and `h
#' = C_K ((sigma-^2 + sigma+^2) / (f (dm^2 + r+ + r-)))^(1/5) N^(-1/5)`
#' with `C_K = (C2 / (4 C1))^(1/5)` from the kernel's moments (their Lemma
#' 3.1; 3.4375 for the triangle, 2.70 for the uniform on `[-1, 1]`, which
#' is their 5.40 for the uniform on `[-1/2, 1/2]`). For a fuzzy design it
#' is the outcome's bandwidth, the simplification the paper offers (p.
#' 946); its fuzzy formula (equation 24) is not built, as the paper does
#' not say at which bandwidth its preliminary `tau` is estimated. `"user"`:
#' the bandwidth `h`. The selector of Calonico, Cattaneo and Titiunik
#' (2014, section 4) is not built.
#'
#' Intervals (`ci`): both are returned in `coef`, the one named first.
#' `"robust"` (the default): Calonico, Cattaneo and Titiunik (2014) bias-
#' correct the local-linear estimate with a local-quadratic one at a pilot
#' bandwidth `b` and widen the variance accordingly; with `b = h` the
#' bias-corrected estimate is the local-quadratic estimate and its robust
#' variance is that estimate's variance (their Remark 7, p. 2305), which
#' is what is computed: the row `"robust"` gives the local-quadratic
#' estimate, its standard error and interval (`b != h` is not built).
#' `"conventional"`: the local-linear estimate with its own standard
#' error, which ignores the smoothing bias. Standard errors (`se`) from
#' the rows' scores: `"hc1"`, heteroskedasticity-robust, or `"cluster"`
#' by `clus`; `ssc = "G"` multiplies by `G / (G - 1)` (`G` the clusters or
#' the rows), `"none"` does not. The normal gives the intervals.
#'
#' Placebo cutoffs (`plac`): at each false cutoff the same sharp fit of
#' `y` (the same `h` and kernel, conventional errors) on the rows of that
#' side of the true cutoff only, so that the true jump cannot enter; jumps
#' found there speak against the design. For the density of the running
#' variable at the cutoff see [rddens()].
#'
#' @param data A data frame.
#' @param y The name of the outcome column.
#' @param run The name of the running-variable column.
#' @param cut The cutoff (treated from `cut` on).
#' @param trt `NULL` (sharp) or the name of the treatment column (fuzzy).
#' @param bw `"ik"` or `"user"`.
#' @param h The bandwidth, for `bw = "user"`.
#' @param kern `"tri"`, `"uni"` or `"epa"`.
#' @param ci `"robust"` or `"conventional"`: the interval named first.
#' @param se `"hc1"` or `"cluster"`.
#' @param clus The name of the cluster column, for `se = "cluster"`.
#' @param ssc `"G"` or `"none"`.
#' @param plac `NULL` or false cutoffs, none equal to `cut`.
#' @param level The confidence level.
#' @param x An object of class `rdloc`.
#' @param object An object of class `rdloc`.
#' @param ... Not used.
#' @return An object of class `rdloc`, a list: `coef`, two rows (`method`,
#'   `est`, `se`, `stat`, `p`, `lo`, `hi`), the one of `ci` first; `h`;
#'   `n`, `nl`, `nr`, the rows within the bandwidth (all, left, right);
#'   `ncl`; `ik`, the steps of the bandwidth (`h1`, `f`, `vl`, `vr`, `m3`,
#'   `h2l`, `h2r`, `m2l`, `m2r`, `ck`; `NULL` for `"user"`); `placebo`
#'   (`cut`, `est`, `se`, `p`, `n`; `NULL` without `plac`); `cut`, `kern`,
#'   `fuzzy`, `se`, `level`. `summary()` returns `coef`.
#' @references
#'   Imbens, G. and Kalyanaraman, K. (2012). Optimal bandwidth choice for
#'   the regression discontinuity estimator. *Review of Economic Studies*,
#'   79(3), 933-959. \doi{10.1093/restud/rdr043}
#'
#'   Calonico, S., Cattaneo, M. D. and Titiunik, R. (2014). Robust
#'   nonparametric confidence intervals for regression-discontinuity
#'   designs. *Econometrica*, 82(6), 2295-2326. \doi{10.3982/ECTA11757}
#'
#'   Carrillo, P. E., Lopez-Luzuriaga, A. and Malik, A. S. (2018).
#'   Pollution or crime: the effect of driving restrictions on criminal
#'   activity. *Journal of Public Economics*, 164, 50-69.
#'   \doi{10.1016/j.jpubeco.2018.05.007}
#' @seealso [rddens()], [ivfe()]
#' @examples
#' set.seed(5)
#' d <- data.frame(r = runif(2000, -1, 1))
#' d$y <- 1 + 0.8 * d$r + 0.5 * (d$r >= 0) + rnorm(2000, 0, 0.3)
#' f <- rdloc(d, "y", "r")
#' f
#' # a given bandwidth, the uniform kernel and two placebo cutoffs
#' summary(rdloc(d, "y", "r", bw = "user", h = 0.3, kern = "uni",
#'               plac = c(-0.5, 0.5)))
#' # fuzzy: treatment taken up more often above the cutoff
#' d$t <- as.numeric(runif(2000) < 0.2 + 0.6 * (d$r >= 0))
#' d$y2 <- 1 + 0.8 * d$r + 0.7 * d$t + rnorm(2000, 0, 0.3)
#' rdloc(d, "y2", "r", trt = "t")$coef
#' @export
rdloc <- function(
  data,
  y,
  run,
  cut = 0,
  trt = NULL,
  bw = c("ik", "user"),
  h = NULL,
  kern = c("tri", "uni", "epa"),
  ci = c("robust", "conventional"),
  se = c("hc1", "cluster"),
  clus = NULL,
  ssc = c("G", "none"),
  plac = NULL,
  level = 0.95
) {
  f <- "rdloc"
  bw <- match.arg(bw)
  kern <- match.arg(kern)
  ci <- match.arg(ci)
  se <- match.arg(se)
  ssc <- match.arg(ssc)
  d <- .axdat(data, f)
  .axcols(d, c(y, run), "y, run", f, n = 2L)
  .axcols(d, trt, "trt", f, n = 1L, null = TRUE)
  yv <- .axmat(d, y, f)[, 1L]
  rv <- .axmat(d, run, f)[, 1L]
  tv <- if (is.null(trt)) NULL else .axmat(d, trt, f)[, 1L]
  cut <- .cknum(cut, "cut", f)
  level <- .axlevel(level, f)
  if (!any(rv < cut) || !any(rv >= cut)) {
    .ckstop(f, "`cut` must have rows on both sides")
  }
  cl <- NULL
  if (se == "cluster") {
    if (is.null(clus)) .ckstop(f, "`se = \"cluster\"` needs `clus`")
    .axcols(d, clus, "clus", f, n = 1L)
    cl <- d[[clus]]
  }
  kc <- match(kern, c("tri", "uni", "epa")) - 1L
  ik <- NULL
  if (bw == "user") {
    if (is.null(h)) .ckstop(f, "`bw = \"user\"` needs `h`")
    h <- .cknum(h, "h", f, positive = TRUE)
  } else {
    if (!is.null(h)) .ckstop(f, "give `h` with `bw = \"user\"`")
    b <- .Call(C_clab_r_ikbw, as.numeric(yv), as.numeric(rv), cut, kc)
    h <- b$h
    ik <- stats::setNames(b$st, c("h1", "f", "vl", "vr", "m3", "h2l", "h2r",
                                  "m2l", "m2r", "ck"))
  }
  r1 <- .rdfit(yv, rv - cut, tv, h, kern, 1L, cl, ssc, f)
  r2 <- .rdfit(yv, rv - cut, tv, h, kern, 2L, cl, ssc, f)
  tab <- .axtab(c("conventional", "robust"), c(r1$est, r2$est),
                c(r1$se, r2$se), level)
  names(tab)[1L] <- "method"
  if (ci == "robust") tab <- tab[2:1, , drop = FALSE]
  rownames(tab) <- NULL
  pl <- NULL
  if (!is.null(plac)) {
    plac <- .cknum(plac, "plac", f, scalar = FALSE)
    if (any(plac == cut)) .ckstop(f, "`plac` must not hold `cut` itself")
    pl <- do.call(rbind, lapply(plac, function(p0) {
      s <- if (p0 < cut) rv < cut else rv >= cut
      q <- .rdfit(yv[s], rv[s] - p0, NULL, h, kern, 1L,
                  if (is.null(cl)) NULL else cl[s], ssc, f)
      data.frame(cut = p0, est = q$est, se = q$se,
                 p = 2 * stats::pnorm(-abs(q$est / q$se)), n = q$n)
    }))
  }
  structure(
    list(coef = tab, h = h, n = r1$n, nl = r1$nl, nr = r1$nr, ncl = r1$g,
         ik = ik, placebo = pl, cut = cut, kern = kern,
         fuzzy = !is.null(trt), se = se, level = level),
    class = "rdloc"
  )
}

#' @rdname rdloc
#' @export
print.rdloc <- function(x, ...) {
  cat("rdloc: ", if (x[["fuzzy"]]) "fuzzy" else "sharp",
      " regression discontinuity at ", format(x[["cut"]]), ", ",
      x[["kern"]], " kernel, bandwidth ", format(x[["h"]], digits = 4),
      if (is.null(x[["ik"]])) " (given)" else " (Imbens-Kalyanaraman)",
      "; ", x[["nl"]], " + ", x[["nr"]], " rows within it\n", sep = "")
  print(x[["coef"]], row.names = FALSE, digits = 4)
  if (!is.null(x[["placebo"]])) {
    cat("placebo cutoffs:\n")
    print(x[["placebo"]], row.names = FALSE, digits = 4)
  }
  invisible(x)
}

#' @rdname rdloc
#' @export
summary.rdloc <- function(object, ...) object[["coef"]]

#' The density test of McCrary at a cutoff
#'
#' Tests whether the density of a running variable jumps at the cutoff,
#' which manipulation of the running variable would cause, by the
#' two-step estimator of McCrary (2008): a finely binned histogram, then
#' local-linear smoothing of the bin heights on each side with the
#' triangle kernel.
#'
#' The steps, as the paper states them. The bins: `g(R) = floor((R - c) /
#' b) b + b / 2 + c` (equation 2, p. 702), so no bin straddles the cutoff,
#' on a grid of width `b` over the support, each height the share of the
#' rows in the bin over `b`; the binsize `b` by default `2 s n^(-1/2)`
#' with `s` the standard deviation of the running variable (p. 705, step
#' 1). The bandwidth `h` by default: a quartic in the bin midpoints on
#' each side, `3.348 (s2 (b - a) / sum f''(X_j)^2)^(1/5)` with `s2` the
#' quartic's mean squared error and `b - a` the side's length, averaged
#' over the two sides (p. 705, step 2; the paper calls it a guide for
#' subjective choice and halves it to undersmooth in its simulations, p.
#' 702-703). The estimate: `theta = ln f+ - ln f-` from the two local-linear
#' intercepts, with the standard error `sqrt((1 / (n h)) (24 / 5) (1 / f+ +
#' 1 / f-))` (equation 5); the test uses the normal.
#'
#' @param r The running variable (numbers).
#' @param cut The cutoff.
#' @param b `NULL` (the rule above) or the binsize.
#' @param h `NULL` (the rule above) or the bandwidth.
#' @param level The confidence level.
#' @param x An object of class `rddens`.
#' @param object An object of class `rddens`.
#' @param ... Not used.
#' @return An object of class `rddens`, a list: `tab`, one row (`theta`,
#'   `se`, `stat`, `p`, `lo`, `hi`, `fl`, `fr`, the densities left and
#'   right); `b`, `h`, `nbin`, `n`, `cut`, `level`. `summary()` returns
#'   `tab`.
#' @references
#'   McCrary, J. (2008). Manipulation of the running variable in the
#'   regression discontinuity design: a density test. *Journal of
#'   Econometrics*, 142(2), 698-714. \doi{10.1016/j.jeconom.2007.05.005}
#' @seealso [rdloc()]
#' @examples
#' set.seed(6)
#' rddens(runif(3000, -1, 1))
#' # bunching just above the cutoff
#' rddens(c(runif(3000, -1, 1), runif(300, 0, 0.1)))
#' @export
rddens <- function(r, cut = 0, b = NULL, h = NULL, level = 0.95) {
  f <- "rddens"
  if (!is.numeric(r) || anyNA(r) || !all(is.finite(r))) {
    .ckstop(f, "`r` must hold finite numbers")
  }
  cut <- .cknum(cut, "cut", f)
  b <- if (is.null(b)) 0 else .cknum(b, "b", f, positive = TRUE)
  h <- if (is.null(h)) 0 else .cknum(h, "h", f, positive = TRUE)
  level <- .axlevel(level, f)
  s <- .Call(C_clab_r_mccr, as.numeric(r), cut, b, h)
  tab <- .axtab("theta", s[1L], s[2L], level)[, -1L]
  tab$fl <- s[3L]
  tab$fr <- s[4L]
  names(tab)[1L] <- "theta"
  structure(
    list(tab = tab, b = s[5L], h = s[6L], nbin = as.integer(s[7L]),
         n = length(r), cut = cut, level = level),
    class = "rddens"
  )
}

#' @rdname rddens
#' @export
print.rddens <- function(x, ...) {
  cat("rddens: McCrary density test at ", format(x[["cut"]]), ", ",
      x[["n"]], " rows, binsize ", format(x[["b"]], digits = 4),
      ", bandwidth ", format(x[["h"]], digits = 4), "\n", sep = "")
  print(x[["tab"]], row.names = FALSE, digits = 4)
  invisible(x)
}

#' @rdname rddens
#' @export
summary.rddens <- function(object, ...) object[["tab"]]
