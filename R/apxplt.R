# Plot methods of the exposure-and-crime family (R/apx.R, R/apx2.R). They
# draw from what each result holds, with the drawing kit of R/plot.R, and
# return invisibly what they drew.

# ---- the methods -------------------------------------------------------------

#' @rdname ppmlfe
#' @inheritParams plotlrnz
#' @param terms `NULL` (every coefficient) or the names of the terms to
#'   draw.
#' @details `plot()` draws the rate ratio `exp(est)` of each term with its
#'   interval on a log scale, and returns the rows drawn. The fits of
#'   [ccross()] are of this class (or of [ppmlcf()]'s), so their spline
#'   terms can be left out with `terms`.
#' @export
plot.ppmlfe <- function(x, ..., main = NULL, col = NULL, pal = NULL,
                        terms = NULL) {
  .pck(x, "ppmlfe", "coef")
  invisible(.axcoef(.axterm(x$coef, terms), TRUE,
                    .pmain(main, "Poisson fixed-effects estimates"), col,
                    pal, ...))
}

#' @rdname ivfe
#' @inheritParams plotlrnz
#' @param terms `NULL` (every coefficient) or the names of the terms to
#'   draw.
#' @details `plot()` draws the coefficients with their intervals and, with
#'   endogenous columns, the first-stage coefficient of each instrument
#'   (`fstab`, with normal intervals at `level`) and the first-stage F
#'   statistics beside the 5% critical value of the effective F; it
#'   returns a list of what it drew (`coef`, `fstab`, `first`).
#' @export
plot.ivfe <- function(x, ..., main = NULL, col = NULL, pal = NULL,
                      terms = NULL) {
  .pck(x, "ivfe", c("coef", "first"))
  iv <- nrow(x$first) > 0L && !is.null(x$fstab)
  op <- .ppanels(if (iv) 3L else 1L)
  on.exit(graphics::par(op), add = TRUE)
  m <- if (iv) "Two-stage least-squares estimates" else
    "Least-squares estimates"
  out <- list(coef = .axcoef(.axterm(x$coef, terms), FALSE, .pmain(main, m),
                             col, pal, ...))
  if (iv) {
    out$fstab <- .axfstab(x$fstab, x$level, col, pal, ...)
    out$first <- .axfst(x$first, col, pal, ...)
  }
  invisible(out)
}

#' @rdname ppmlcf
#' @inheritParams plotlrnz
#' @param terms `NULL` (every coefficient) or the names of the terms to
#'   draw.
#' @details `plot()` draws the rate ratios with their intervals on a log
#'   scale and the first-stage F statistics beside the 5% critical value
#'   of the effective F; it returns a list of what it drew (`coef`,
#'   `first`).
#' @export
plot.ppmlcf <- function(x, ..., main = NULL, col = NULL, pal = NULL,
                        terms = NULL) {
  .pck(x, "ppmlcf", c("coef", "first"))
  op <- .ppanels(2L)
  on.exit(graphics::par(op), add = TRUE)
  m <- if (identical(x$how, "plug")) "Poisson plug-in estimates" else
    "Poisson control-function estimates"
  out <- list(coef = .axcoef(.axterm(x$coef, terms), TRUE, .pmain(main, m),
                             col, pal, ...),
              first = .axfst(x$first, col, pal, ...))
  invisible(out)
}

#' @rdname expgmm
#' @inheritParams plotlrnz
#' @param terms `NULL` (every coefficient) or the names of the terms to
#'   draw.
#' @details `plot()` draws the rate ratio `exp(est)` of each term with its
#'   interval on a log scale, and returns the rows drawn.
#' @export
plot.expgmm <- function(x, ..., main = NULL, col = NULL, pal = NULL,
                        terms = NULL) {
  .pck(x, "expgmm", "coef")
  invisible(.axcoef(.axterm(x$coef, terms), TRUE,
                    .pmain(main, "Exponential-mean GMM estimates"), col, pal,
                    ...))
}

#' @rdname plcday
#' @inheritParams plotlrnz
#' @details `plot()` draws the coefficient of each shift with its interval,
#'   from the earliest day to the latest, the same day in its own colour
#'   and shape, and returns the rows drawn.
#' @export
plot.plcday <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "plcday", "tab")
  t <- x$tab[order(x$tab$shift), , drop = FALSE]
  s <- t$shift
  lab <- paste(abs(s), ifelse(s < 0, "days before", "days after"))
  lab[s == 0] <- "same day"
  invisible(.pdots(lab, t$est, t$lo, t$hi,
                   .pmain(main, "Placebo: the exposure of other days"),
                   paste("Coefficient of", x$var), col, pal, ref = 0,
                   grp = ifelse(s == 0, "same day", "other days"), ...))
}

#' @rdname updwn
#' @inheritParams plotlrnz
#' @details `plot()` draws the downwind difference with its interval (by
#'   distance, at the midpoints of the rings, with `dist`), the
#'   estimates of the placebo segments by offset with the source itself
#'   marked (with `offs`), and the rows of each classification; it returns
#'   a list of what it drew (`coef`, `placebo`, `class`). The object holds
#'   the difference of the two sides, not each side's level, so the sides
#'   are not drawn apart.
#' @export
plot.updwn <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "updwn", c("coef", "class"))
  pl <- x$placebo
  op <- .ppanels(2L + !is.null(pl))
  on.exit(graphics::par(op), add = TRUE)
  cf <- x$coef
  m <- .pmain(main, "Downwind against upwind")
  yl <- "Downwind difference"
  out <- list()
  if (nrow(cf) > 1L) {
    out$coef <- .axprof(.axring(cf$term)$mid, cf$est, cf$lo, cf$hi, m,
                        "Distance from the source (ring midpoint)", yl, col,
                        pal, ...)
  } else {
    out$coef <- .pdots(cf$term, cf$est, cf$lo, cf$hi, m, yl, col, pal,
                       ref = 0, ...)
  }
  if (!is.null(pl)) {
    out$placebo <- .axprof(pl$off, pl$est, pl$lo, pl$hi,
                           "Placebo segments moved sideways", "Offset", yl,
                           col, pal, mark = 0,
                           key = c("moved", "the source itself"), ...)
  }
  k <- x$class
  out$class <- .pbars(as.numeric(k), names(k), "Rows by classification",
                      "Rows", col, pal, ...)
  invisible(out)
}

#' @rdname exprun
#' @inheritParams plotlrnz
#' @details `plot()` draws the effect of each day of a run with its
#'   interval (the last day standing for that day and later), and returns
#'   the rows drawn.
#' @export
plot.runfx <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "runfx", "coef")
  cf <- x$coef
  k <- .axrun(cf$term)
  lab <- ifelse(grepl("p$", cf$term), paste0(k, "+"), k)
  invisible(.axprof(k, cf$est, cf$lo, cf$hi,
                    .pmain(main, "Effect by day of a run of exposed days"),
                    "Day of the run", "Coefficient", col, pal, at = k,
                    lab = lab, ...))
}

#' @rdname hetfx
#' @inheritParams plotlrnz
#' @details `plot()` draws the effect of each level with its interval, the
#'   p-value of the test of equal effects in the axis title, and returns
#'   the rows drawn.
#' @export
plot.hetfx <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "hetfx", c("tab", "test"))
  t <- x$tab
  xl <- paste0("Coefficient (equal effects: p ",
               signif(x$test$p, 2), ")")
  invisible(.pdots(t$level, t$est, t$lo, t$hi,
                   .pmain(main, "Effect by level of the moderator"), xl, col,
                   pal, ref = 0, ...))
}

#' @rdname jntx
#' @param main,col,pal The title, colour and palette, as in [plotm]; `pal`
#'   shades the surface (by default `"diverging"`).
#' @details `plot()` draws the rate ratio of `x1` (`p1[2]` against
#'   `p1[1]`) at each quantile of `x2` with its interval on a log scale,
#'   and the surface of rate ratios against the reference point; it
#'   returns a list of what it drew (`cond`, `surf`).
#' @export
plot.jntx <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "jntx", c("cond", "surf"))
  op <- .ppanels(2L)
  on.exit(graphics::par(op), add = TRUE)
  cd <- x$cond
  lab <- paste0("x2 at ", 100 * cd$p2, "% (", signif(cd$x2, 3), ")")
  out <- list(cond = .pdots(lab, exp(cd$est), exp(cd$lo), exp(cd$hi),
                            .pmain(main, "Effect of x1 at levels of x2"),
                            "Rate ratio of x1 (log scale)", col, pal,
                            ref = 1, log = "x", ...))
  out$surf <- .axsurf(x$surf, pal, ...)
  invisible(out)
}

#' @rdname dosebin
#' @inheritParams plotlrnz
#' @details `plot()` draws the effect of each bin with its interval at the
#'   bin's median exposure, the reference bin (0) marked, the p-values of
#'   the two tests in the axis title, and returns the rows drawn.
#' @export
plot.dosebin <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "dosebin", c("tab", "test"))
  t <- x$tab
  tt <- x$test
  xl <- paste0("Exposure, median of each bin (",
               paste(tt$test, "p", signif(tt$p, 2),
                     collapse = "; "), ")")
  invisible(.axprof(t$mid, t$est, t$lo, t$hi,
                    .pmain(main, "Effect of each bin against the reference"),
                    xl, "Coefficient", col, pal,
                    mark = t$mid[is.na(t$se) & t$est == 0],
                    key = c("bin", "reference bin"), ...))
}

#' @rdname typefx
#' @inheritParams plotlrnz
#' @details `plot()` draws the effect on each outcome with its interval,
#'   the adjusted p-value beside the outcome's name, and returns the rows
#'   drawn.
#' @export
plot.typefx <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "typefx", "tab")
  t <- x$tab
  lab <- paste0(t$y, " (adjusted p ", signif(t$padj, 2), ")")
  invisible(.pdots(lab, t$est, t$lo, t$hi,
                   .pmain(main, "Effect by outcome"), "Coefficient", col,
                   pal, ref = 0, ...))
}

# ---- helpers -----------------------------------------------------------------

# normal intervals at `level` about estimates with standard errors
.axci <- function(est, se, level) {
  z <- stats::qnorm(0.5 + level / 2)
  list(lo = est - z * se, hi = est + z * se)
}

# the rings of an updwn() fit from its term names "down:(a,b]": the ends
# and the midpoint of each (cut() prints the breaks to three digits, so
# the positions are drawn to that precision)
.axring <- function(term) {
  v <- gsub("^[^:]*:[[(]|[])]$", "", term)
  e <- matrix(as.numeric(unlist(strsplit(v, ",", fixed = TRUE))), ncol = 2L,
              byrow = TRUE)
  data.frame(lo = e[, 1L], hi = e[, 2L], mid = (e[, 1L] + e[, 2L]) / 2)
}

# the day of a run of each runfx() term: "run2" is 2, "run3p" is 3
.axrun <- function(term) as.numeric(gsub("^run|p$", "", term))

# the coefficient rows drawn: every row, or the rows of `terms`
.axterm <- function(cf, terms) {
  if (is.null(terms)) {
    return(cf)
  }
  if (!is.character(terms) || !length(terms) || anyNA(terms) ||
        !all(terms %in% cf$term)) {
    .ckstop("plot", "`terms` must be NULL or names of terms of the fit: %s",
            toString(cf$term))
  }
  cf[cf$term %in% terms, , drop = FALSE]
}

# coefficients as a forest chart: rate ratios on a log scale for the
# Poisson and exponential-mean fits, the coefficients for least squares
.axcoef <- function(cf, rr, main, col = NULL, pal = NULL, ...) {
  if (rr) {
    .pdots(cf$term, exp(cf$est), exp(cf$lo), exp(cf$hi), main,
           "Rate ratio (log scale)", col, pal, ref = 1, log = "x", ...)
  } else {
    .pdots(cf$term, cf$est, cf$lo, cf$hi, main, "Coefficient", col, pal,
           ref = 0, ...)
  }
}

# the strength of the first stage of each endogenous column: the robust
# Wald F, the Kleibergen-Paap F and the effective F beside the 5% critical
# value of Montiel Olea and Pflueger
.axfst <- function(fs, col = NULL, pal = NULL, ...) {
  k <- c("Wald F", "Kleibergen-Paap F", "effective F", "5% critical value")
  .pdots(rep(fs$endo, 4L), c(fs$wald, fs$kp, fs$eff, fs$cv), NA, NA,
         "First-stage strength", "F statistic", col, pal,
         grp = rep(k, each = nrow(fs)), ...)
}

# the first-stage coefficient of each excluded instrument, with normal
# intervals at the fit's level
.axfstab <- function(fb, level, col = NULL, pal = NULL, ...) {
  ci <- .axci(fb$est, fb$se, level)
  .pdots(paste0(fb$inst, " (", fb$endo, ")"), fb$est, ci$lo, ci$hi,
         "First-stage coefficients", "Change in the exposure", col, pal,
         ref = 0, ...)
}

# estimates at positions on a numeric axis (days of a run, bin medians,
# distances, offsets) with their intervals, joined by a dotted line; the
# positions in `mark` (the reference, or the true design) are drawn as
# squares in the second colour, named in the legend by `key`; `at` and
# `lab` give the ticks of the axis; returns the rows drawn
.axprof <- function(x, est, lo, hi, main, xlab, ylab, col = NULL, pal = NULL,
                    mark = NULL, key = NULL, at = NULL, lab = NULL, ...) {
  cl <- .pcols(2L, col, pal)
  op <- .ppar()
  on.exit(graphics::par(op), add = TRUE)
  .pfrm(.prange(x), .prange(est, lo, hi, 0), main, xlab, ylab,
        xaxt = if (is.null(at)) "s" else "n", ...)
  if (!is.null(at)) graphics::axis(1, at = at, labels = lab)
  graphics::abline(h = 0, lty = 2, col = "#636363")
  m <- x %in% mark
  o <- order(x)
  graphics::lines(x[o], est[o], col = cl[1L], lty = 3)
  w <- is.finite(lo) & is.finite(hi)
  graphics::segments(x[w], lo[w], x[w], hi[w], col = cl[m + 1L][w],
                     lwd = 1.6)
  graphics::points(x, est, pch = ifelse(m, 15, 19), col = cl[m + 1L],
                   cex = 1.2)
  if (any(m) && !is.null(key)) .pkey("topright", key, cl, pch = c(19, 15))
  invisible(data.frame(x = x, est = est, lo = lo, hi = hi, mark = m))
}

# the surface of a jntx() fit: the log rate ratio shaded on a scale
# centred at 0, with contour lines labelled by the rate ratio
.axsurf <- function(s, pal = NULL, ...) {
  u <- sort(unique(s$x1))
  v <- sort(unique(s$x2))
  z <- matrix(log(s$rr[order(s$x2, s$x1)]), length(u), length(v))
  sh <- blpal(.palnm(if (is.null(pal)) "diverging" else pal, "plot"), 9L)
  r <- max(abs(.prange(z)))
  op <- .ppar()
  on.exit(graphics::par(op), add = TRUE)
  graphics::image(u, v, z, col = sh, zlim = c(-r, r), main = "",
                  xlab = "First exposure (x1)", ylab = "Second exposure (x2)",
                  ...)
  graphics::contour(u, v, exp(z), add = TRUE, col = "#636363", nlevels = 6,
                    labcex = 0.7)
  graphics::title(main = "Rate ratio against the reference point",
                  font.main = 1, cex.main = 1.05)
  graphics::box(col = "#8C8C8C")
  invisible(s)
}
