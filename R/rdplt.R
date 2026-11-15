# Plot methods of the regression-discontinuity results (R/rdloc.R), and
# the pieces they draw: rdloc() and rddens() keep them in the attribute
# "plt" of their result, computed in R from the data of the fit, so that
# the plot draws without refitting and the elements of the result stay
# as they are.

#' @rdname rdloc
#' @inheritParams plotlrnz
#' @return The attribute `"plt"` of the result holds what `plot()` draws:
#'   `bins`, the outcome's means in 20 equal bins of the running variable
#'   on each side (`x` the midpoint, `y`, `n`, `side`), and `fit`, the
#'   ends of the local-linear line of each side over the rows within the
#'   bandwidth, with the kernel's weights (`x`, `y`, `side`).
#' @details `plot()` draws the binned means and the local-linear line of
#'   each side, whose gap at the cutoff is the conventional estimate of a
#'   sharp design (for a fuzzy one, the jump in the outcome, the numerator
#'   of the ratio), and with `plac` the jump at each placebo cutoff with
#'   its normal interval; it returns a list of what it drew (`bins`,
#'   `fit`, `placebo`).
#' @export
plot.rdloc <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "rdloc", c("coef", "cut"))
  p <- .rdpget(x, "rdloc")
  pl <- x$placebo
  op <- .ppanels(1L + !is.null(pl))
  on.exit(graphics::par(op), add = TRUE)
  out <- .rddraw(p, x$cut, .pmain(main, "Regression discontinuity"),
                 "Running variable", "Outcome: bin means and fits", col, pal,
                 ...)
  if (!is.null(pl)) {
    ci <- .axci(pl$est, pl$se, x$level)
    out$placebo <- .pdots(paste("cutoff", signif(pl$cut, 4)), pl$est, ci$lo,
                          ci$hi, "Jumps at placebo cutoffs", "Jump", col,
                          pal, ref = 0, ...)
  }
  invisible(out)
}

#' @rdname rddens
#' @inheritParams plotlrnz
#' @return The attribute `"plt"` of the result holds what `plot()` draws:
#'   `bins`, the first-step histogram (`x` the bin midpoint, `y` the
#'   height, `side`), and `fit`, the local-linear smooth of the heights on
#'   each side at 40 points up to the cutoff (`x`, `y`, `side`), whose
#'   values at the cutoff are `fl` and `fr`.
#' @details `plot()` draws the bin heights and the smooth of each side,
#'   the log jump and its p-value in the axis title, and returns a list of
#'   what it drew (`bins`, `fit`).
#' @export
plot.rddens <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "rddens", c("tab", "cut"))
  p <- .rdpget(x, "rddens")
  t <- x$tab
  xl <- paste0("Running variable (log jump ", signif(t$theta, 3), ", p ",
               signif(t$p, 2), ")")
  invisible(.rddraw(p, x$cut, .pmain(main, "Density at the cutoff"), xl,
                    "Density: bin heights and fits", col, pal, ...))
}

# the weighted least-squares line y = a + b u: c(a, b), NA when the
# weights rest on fewer than two distinct values of u
.rdwls <- function(u, y, w) {
  s0 <- sum(w)
  s1 <- sum(w * u)
  s2 <- sum(w * u^2)
  t0 <- sum(w * y)
  t1 <- sum(w * u * y)
  dn <- s0 * s2 - s1^2
  if (dn <= 1e-12 * s0 * s2) {
    return(c(NA_real_, NA_real_))
  }
  c((s2 * t0 - s1 * t1) / dn, (s0 * t1 - s1 * t0) / dn)
}

# the means of y in nb equal bins of u over lim (one side of the cutoff,
# so no bin straddles it): the bins' midpoints x, means y and rows n; the
# empty bins left out
.rdbin <- function(u, y, lim, nb, side) {
  if (lim[2L] <= lim[1L]) {
    return(data.frame(x = lim[1L], y = mean(y), n = length(y), side = side))
  }
  e <- seq(lim[1L], lim[2L], length.out = nb + 1L)
  g <- findInterval(u, e, rightmost.closed = TRUE, all.inside = TRUE)
  n <- tabulate(g, nb)
  m <- as.numeric(tapply(y, factor(g, levels = seq_len(nb)), mean))
  k <- n > 0L
  data.frame(x = ((e[-1L] + e[-(nb + 1L)]) / 2)[k], y = m[k], n = n[k],
             side = rep(side, sum(k)))
}

# what plot.rdloc() draws: the outcome's means in nb bins on each side,
# and the local-linear line of each side over the rows the kernel weighs
# (its intercepts differ by the conventional sharp estimate)
.rdplt <- function(yv, rv, cut, h, kern, nb = 20L) {
  u <- rv - cut
  lf <- u < 0
  bins <- rbind(.rdbin(u[lf], yv[lf], c(min(u), 0), nb, "left"),
                .rdbin(u[!lf], yv[!lf], c(0, max(u)), nb, "right"))
  bins$x <- bins$x + cut
  w <- .rdkw(u / h, kern)
  fit <- do.call(rbind, lapply(c("left", "right"), function(s) {
    k <- w > 0 & (if (s == "left") lf else !lf)
    ab <- .rdwls(u[k], yv[k], w[k])
    e <- if (s == "left") c(min(u[k]), 0) else c(0, max(u[k]))
    data.frame(x = cut + e, y = ab[1L] + ab[2L] * e, side = s)
  }))
  list(bins = bins, fit = fit)
}

# what plot.rddens() draws: McCrary's first-step histogram (the share of
# the rows in each bin of width b, over b, at the bins' midpoints) and
# the local-linear smooth of the heights on each side, triangle kernel of
# bandwidth h, at ng points from the side's last bin to the cutoff
.rddplt <- function(r, cut, b, h, ng = 40L) {
  j <- floor((r - cut) / b)
  lo <- min(j)
  m <- as.integer(max(j) - lo + 1)
  x <- (lo + seq_len(m) - 1) * b + b / 2 + cut
  bins <- data.frame(x = x, y = tabulate(as.integer(j - lo + 1), m) /
                       (length(r) * b),
                     side = ifelse(x < cut, "left", "right"))
  fit <- do.call(rbind, lapply(c("left", "right"), function(s) {
    k <- bins$side == s
    if (!any(k)) {
      return(NULL)
    }
    xs <- bins$x[k] - cut
    g <- seq(if (s == "left") min(xs) else 0, if (s == "left") 0 else max(xs),
             length.out = ng)
    f <- vapply(g, function(g0) {
      .rdwls(xs - g0, bins$y[k], pmax(1 - abs((xs - g0) / h), 0))[1L]
    }, 0)
    data.frame(x = g + cut, y = f, side = s)
  }))
  list(bins = bins, fit = fit[is.finite(fit$y), , drop = FALSE])
}

# the pieces kept by the fit, or a refusal for a result without them
.rdpget <- function(x, cls) {
  p <- attr(x, "plt")
  if (!is.list(p) || is.null(p$bins) || is.null(p$fit)) {
    .ckstop("plot", "`x` lacks the binned data of the plot: refit it with %s()",
            cls)
  }
  p
}

# the bins as points (open below the cutoff, filled from it on) and the
# fitted lines of the two sides, the cutoff dashed; returns both
.rddraw <- function(p, cut, main, xlab, ylab, col = NULL, pal = NULL, ...) {
  cl <- .pcols(2L, col, pal)
  b <- p$bins
  f <- p$fit
  op <- .ppar()
  on.exit(graphics::par(op), add = TRUE)
  .pfrm(.prange(b$x, f$x), .prange(b$y, f$y), main, xlab, ylab, ...)
  graphics::abline(v = cut, lty = 2, col = "#636363")
  s <- (b$side == "right") + 1L
  graphics::points(b$x, b$y, pch = c(1, 19)[s], col = cl[s])
  for (k in 1:2) {
    w <- f$side == c("left", "right")[k]
    graphics::lines(f$x[w], f$y[w], col = cl[k], lty = k, lwd = 2)
  }
  .pkey("topleft", c("below the cutoff", "from the cutoff on"), cl,
        lty = 1:2, pch = c(1, 19))
  invisible(list(bins = b, fit = f))
}
