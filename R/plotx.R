#' Lorenz and funnel plots
#'
#' `plotlrnz()` draws the Lorenz curve of a total spread over units: the
#' cumulative share of the total held by the cumulative share of units,
#' smallest first, against the line of equality, with the Gini coefficient
#' (from [cconc()]) in the title. `plotfnl()` draws a funnel plot
#' (Spiegelhalter 2005): each area's ratio of observed to expected counts
#' against its expected count, with exact Poisson control limits around
#' `target`; areas outside the outer limits are drawn in a second colour
#' and a second shape, and labelled when `area` is given.
#'
#' @details The funnel limits at expected count \eqn{E} and level \eqn{p} are
#'   the Poisson quantiles \eqn{q_{(1-p)/2}} and \eqn{q_{(1+p)/2}} of a count
#'   with mean `target` \eqn{E}, divided by \eqn{E}; they are drawn through
#'   the areas' own expected counts. An area with either count missing is
#'   dropped whole.
#'
#' @param x For `plotlrnz()`, non-negative values, one per unit, or a data
#'   frame with columns `pop` and `share` (a Lorenz curve already computed,
#'   such as `lorenz` of [cconc()]).
#' @param obs,expd Observed and expected counts, one per area; `expd`
#'   positive.
#' @param area `NULL` or a label per area.
#' @param lvls The levels of the control limits: probabilities strictly
#'   between 0 and 1.
#' @param target The ratio expected under no excess: one positive number.
#' @param ... Passed to the base-graphics call that draws the frame.
#' @param main `NULL` (a default title) or one string.
#' @param col `NULL` (the palette) or colours.
#' @param pal `NULL` (the session palette, see [plotopt()]) or a palette
#'   name for [blpal()].
#' @return Invisibly, the data drawn: for `plotlrnz()` the curve (`pop`,
#'   `share`) from the origin, with the Gini as attribute `gini` (`NA` for
#'   a curve given as a frame); for `plotfnl()` a data frame with `area`,
#'   `observed`, `expected`, `ratio` and `outside` (beyond the outer
#'   limits).
#' @references Spiegelhalter, D. J. (2005). Funnel plots for comparing
#'   institutional performance. *Statistics in Medicine*, 24(8), 1185-1202.
#'   \doi{10.1002/sim.1970}
#' @examples
#' plotlrnz(c(0, 0, 1, 1, 2, 3, 8, 25))
#' plotfnl(obs = c(4, 20, 50, 2), expd = c(5, 22, 40, 20),
#'         area = c("a", "b", "c", "d"))
#' @export
plotlrnz <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  if (is.data.frame(x)) {
    d <- .plz(x)
    g <- NA_real_
  } else {
    r <- .pcc(x, "plotlrnz")
    d <- r[["lorenz"]]
    g <- r[["stat"]][["gini"]]
  }
  cl <- .pcols(1L, col, pal, "plotlrnz")
  op <- .ppar()
  on.exit(graphics::par(op), add = TRUE)
  .pfrm(c(0, 1), c(0, 1), .pmain(main, if (is.na(g)) "Lorenz curve" else
    sprintf("Lorenz curve (Gini %.3f)", g), "plotlrnz"),
    "Cumulative share of units, smallest first",
    "Cumulative share of the total", ...)
  graphics::abline(0, 1, col = "#636363", lty = 2)
  graphics::lines(d[["pop"]], d[["share"]], col = cl, lwd = 2.2)
  attr(d, "gini") <- g
  invisible(d)
}

#' @rdname plotlrnz
#' @export
plotfnl <- function(obs, expd, area = NULL, lvls = c(0.95, 0.998),
                    target = 1, ..., main = NULL, col = NULL, pal = NULL) {
  k <- .pfchk(obs, expd, area, lvls, target)
  obs <- k$obs
  expd <- k$expd
  keep <- k$keep
  lv <- k$lv
  # the limits are drawn through the areas' own expected counts
  e <- sort(unique(expd))
  lim <- function(p, at) {
    a <- (1 - p) / 2
    cbind(stats::qpois(a, target * at) / at,
          stats::qpois(1 - a, target * at) / at)
  }
  ratio <- obs / expd
  ob <- lim(max(lv), expd)
  out <- ratio < ob[, 1L] | ratio > ob[, 2L]
  cl <- .pcols(2L, col, if (is.null(pal) && is.null(col)) "diverging" else
    pal, "plotfnl")
  cl <- c(cl[1L], cl[length(cl)])
  op <- .ppar()
  on.exit(graphics::par(op), add = TRUE)
  bands <- lapply(lv, lim, at = e)
  .pfrm(.prange(e), .prange(ratio, unlist(bands), target),
          .pmain(main, "Funnel plot", "plotfnl"), "Expected count",
          "Observed / expected", ...)
  graphics::abline(h = target, col = "#636363")
  for (i in seq_along(lv)) {
    graphics::matlines(e, bands[[i]], col = "#636363",
                       lty = if (i == length(lv)) 1 else 2)
  }
  graphics::points(expd, ratio, pch = ifelse(out, 17, 19),
                   col = cl[out + 1L], cex = 1.3)
  if (!is.null(area) && any(out)) {
    graphics::text(expd[out], ratio[out], as.character(area[keep])[out],
                   pos = 3, cex = 0.75)
  }
  invisible(data.frame(
    area = if (is.null(area)) NA_character_ else as.character(area[keep]),
    observed = obs, expected = expd, ratio = ratio, outside = out,
    stringsAsFactors = FALSE
  ))
}

#' Distributions by group: ridgelines, distribution functions, violins,
#' hexagonal bins and concentration curves
#'
#' Plots of many values, split by a group. `plotrdg()` stacks a kernel
#' density per group (a ridgeline plot); `plotcdfs()` overlays the
#' empirical distribution functions; `plotvln()` draws a violin per group
#' with its median and quartiles; `plothex()` counts two variables in
#' hexagonal bins, shaded by count; `plotgain()` overlays the
#' concentration (Lorenz) curves of the groups with their Gini
#' coefficients.
#'
#' @param x Numbers (non-negative for `plotgain()`).
#' @param y For `plothex()`, numbers, one per value of `x`.
#' @param grp `NULL` (one group) or a group per value of `x`.
#' @param bw The kernel bandwidth: one positive number or a rule known to
#'   [stats::density()] (`"nrd0"`, `"nrd"`, `"ucv"`, `"bcv"`, `"sj"`,
#'   `"sj-ste"`, `"sj-dpi"`).
#' @param bins The number of hexagons across: one whole number from 2 to
#'   500.
#' @inheritParams plotlrnz
#' @return Invisibly, a data frame of what was drawn: `grp`, `x` and `y`
#'   (the density, or the share at or below), for `plotvln()` with the
#'   quartiles as attribute `summary`; for `plothex()` the hexagon centres
#'   `x`, `y` and `count`, by row of hexagons and then column; for
#'   `plotgain()` `grp`, `pop` and `share`, with the Gini per group as
#'   attribute `gini`. Values that are not finite are left out.
#' @examples
#' set.seed(1)
#' x <- c(rnorm(100), rnorm(100, 2))
#' g <- rep(c("a", "b"), each = 100)
#' plotrdg(x, g)
#' plotcdfs(x, g)
#' plotvln(x, g)
#' plothex(x, x + rnorm(200), bins = 10)
#' plotgain(rexp(200), g)
#' @export
plotrdg <- function(x, grp = NULL, bw = "nrd0", ..., main = NULL,
                    col = NULL, pal = NULL) {
  d <- .pxg(x, grp, "plotrdg")
  .pbw(bw, "plotrdg")
  lv <- levels(d$grp)
  out <- do.call(rbind, lapply(lv, function(g) {
    v <- d$x[d$grp == g]
    if (length(v) < 2L) {
      return(NULL)
    }
    k <- .pdens(v, bw)
    data.frame(grp = g, x = k$x, y = k$y, stringsAsFactors = FALSE)
  }))
  if (is.null(out)) {
    stop("plotrdg: no group has two values.", call. = FALSE)
  }
  cl <- .pcols(length(lv), col, pal, "plotrdg")
  op <- .ppar(mar = c(4.2, .plmar(lv), 3, 1))
  on.exit(graphics::par(op), add = TRUE)
  top <- max(out$y)
  step <- 0.7 * top
  .pfrm(range(out$x), c(0, step * (length(lv) - 1L) + top),
          .pmain(main, "Distribution by group", "plotrdg"), "Value",
          yaxt = "n", ...)
  for (i in rev(seq_along(lv))) {
    s <- out[out$grp == lv[i], ]
    b <- step * (i - 1L)
    graphics::polygon(c(s$x, rev(s$x)), c(b + s$y, rep(b, nrow(s))),
                      col = grDevices::adjustcolor(cl[i], 0.5),
                      border = "white")
    graphics::lines(s$x, b + s$y, col = cl[i], lwd = 1.8, lty = i)
  }
  graphics::axis(2, at = step * (seq_along(lv) - 1L), labels = lv,
                 tick = FALSE)
  invisible(out)
}

#' @rdname plotrdg
#' @export
plotcdfs <- function(x, grp = NULL, ..., main = NULL, col = NULL,
                     pal = NULL) {
  d <- .pxg(x, grp, "plotcdfs")
  lv <- levels(d$grp)
  cl <- .pcols(length(lv), col, pal, "plotcdfs")
  op <- .ppar()
  on.exit(graphics::par(op), add = TRUE)
  .pfrm(range(d$x), c(0, 1),
          .pmain(main, "Empirical distribution functions", "plotcdfs"),
          "Value", "Share at or below", ...)
  out <- do.call(rbind, lapply(seq_along(lv), function(i) {
    v <- sort(d$x[d$grp == lv[i]])
    p <- seq_along(v) / length(v)
    graphics::lines(c(v[1L], v), c(0, p), type = "s", col = cl[i], lty = i,
                    lwd = 2)
    data.frame(grp = lv[i], x = v, y = p, stringsAsFactors = FALSE)
  }))
  if (length(lv) > 1L) .pkey("bottomright", lv, cl, lty = seq_along(lv))
  invisible(out)
}

#' @rdname plotrdg
#' @export
plotvln <- function(x, grp = NULL, bw = "nrd0", ..., main = NULL,
                    col = NULL, pal = NULL) {
  d <- .pxg(x, grp, "plotvln")
  .pbw(bw, "plotvln")
  lv <- levels(d$grp)
  cl <- .pcols(length(lv), col, pal, "plotvln")
  op <- .ppar()
  on.exit(graphics::par(op), add = TRUE)
  .pfrm(c(0.5, length(lv) + 0.5), range(d$x),
          .pmain(main, "Violins: density, median and quartiles", "plotvln"),
          "", "Value", xaxt = "n", ...)
  out <- NULL
  sm <- NULL
  for (i in seq_along(lv)) {
    v <- d$x[d$grp == lv[i]]
    q <- stats::quantile(v, c(0.25, 0.5, 0.75), names = FALSE)
    sm <- rbind(sm, data.frame(grp = lv[i], median = q[2L], q1 = q[1L],
                               q3 = q[3L], stringsAsFactors = FALSE))
    if (length(v) >= 2L) {
      k <- .pdens(v, bw)
      w <- 0.42 * k$y / max(k$y)
      graphics::polygon(c(i - w, rev(i + w)), c(k$x, rev(k$x)),
                        col = grDevices::adjustcolor(cl[i], 0.5),
                        border = cl[i])
      out <- rbind(out, data.frame(grp = lv[i], x = k$x, y = k$y,
                                   stringsAsFactors = FALSE))
    }
    graphics::segments(i, q[1L], i, q[3L], lwd = 5, col = "#252525")
    graphics::points(i, q[2L], pch = 21, bg = "white", col = "#252525",
                     cex = 1.3)
  }
  graphics::axis(1, at = seq_along(lv), labels = lv, tick = FALSE)
  attr(out, "summary") <- sm
  invisible(out)
}

#' @rdname plotrdg
#' @export
plothex <- function(x, y, bins = 30L, ..., main = NULL, col = NULL,
                    pal = NULL) {
  .pxg(x, NULL, "plothex")
  if (!is.numeric(y) || length(y) != length(x)) {
    stop("plothex: `y` must be numbers, one per value of `x`.",
         call. = FALSE)
  }
  ok <- is.finite(x) & is.finite(y)
  x <- as.numeric(x[ok])
  y <- as.numeric(y[ok])
  if (length(x) < 2L) {
    stop("plothex: needs at least two points with both values.",
         call. = FALSE)
  }
  if (!.pnum1(bins, 2, 500, TRUE)) {
    stop("plothex: `bins` must be one whole number from 2 to 500.",
         call. = FALSE)
  }
  rx <- .prange(x)
  ry <- .prange(y)
  out <- .phex(x, y, bins, rx, ry)
  sh <- if (is.null(col)) {
    blpal(.palnm(if (is.null(pal)) "sequential" else pal, "plothex"), 9L)
  } else {
    .pcols(9L, col, NULL, "plothex")
  }
  lvl <- pmin(9L, 1L + floor(8 * (out$count - 1) /
                               max(1, max(out$count) - 1)))
  op <- .ppar()
  on.exit(graphics::par(op), add = TRUE)
  .pfrm(rx, ry, .pmain(main, sprintf("Counts in hexagons (most: %d)",
                                       max(out$count)), "plothex"),
          "x", "y", ...)
  a <- seq(pi / 6, 2 * pi, by = pi / 3)
  hx <- cos(a) / sqrt(3) / bins * diff(rx)
  hy <- sin(a) / sqrt(3) / bins * diff(ry)
  for (i in seq_len(nrow(out))) {
    graphics::polygon(out$x[i] + hx, out$y[i] + hy, col = sh[lvl[i]],
                      border = NA)
  }
  invisible(out)
}

#' @rdname plotrdg
#' @export
plotgain <- function(x, grp = NULL, ..., main = NULL, col = NULL,
                     pal = NULL) {
  d <- .pxg(x, grp, "plotgain")
  if (any(d$x < 0)) {
    stop("plotgain: `x` must not be negative (a share of a total).",
         call. = FALSE)
  }
  lv <- levels(d$grp)
  cs <- lapply(lv, function(g) .pcc(d$x[d$grp == g], "plotgain"))
  gini <- stats::setNames(vapply(cs, function(r) r$stat[["gini"]], 1), lv)
  cl <- .pcols(length(lv), col, pal, "plotgain")
  op <- .ppar()
  on.exit(graphics::par(op), add = TRUE)
  .pfrm(c(0, 1), c(0, 1), .pmain(main, "Concentration curves",
                                   "plotgain"),
          "Cumulative share of units, smallest first",
          "Cumulative share of the total", ...)
  graphics::abline(0, 1, col = "#636363", lty = 2)
  out <- do.call(rbind, lapply(seq_along(lv), function(i) {
    l <- cs[[i]]$lorenz
    graphics::lines(l$pop, l$share, col = cl[i], lty = i, lwd = 2.2)
    data.frame(grp = lv[i], pop = l$pop, share = l$share,
               stringsAsFactors = FALSE)
  }))
  .pkey("topleft", sprintf("%s (Gini %.3f)", lv, gini), cl,
        lty = seq_along(lv))
  attr(out, "gini") <- gini
  invisible(out)
}

#' Many series over time: calendar, horizon, small-multiple and bump charts
#'
#' `plotcal()` draws daily values as a calendar: one row of weeks per year,
#' each day a cell shaded by its value. `plothrz()` draws a horizon chart:
#' each series' distance from its own median, folded into `bands` bands of
#' darkening shade, blue above and red below, one row per series.
#' `plotmult()` draws each series in its own panel on a shared scale (small
#' multiples). `plotbump()` draws the rank of each series in each period
#' (a bump chart).
#'
#' @param x For `plotcal()`, dates (or text `as.Date()` reads), one per
#'   value of `y`.
#' @param y Numbers.
#' @param prd The period of each value: numbers, dates, a factor or text
#'   (in order of first appearance).
#' @param grp `NULL` (one series) or the series of each value; required by
#'   `plotbump()`.
#' @param bands The number of bands: one whole number from 1 to 6.
#' @param desc `TRUE` to rank the largest value first.
#' @inheritParams plotlrnz
#' @return Invisibly, a data frame of what was drawn: for `plotcal()` the
#'   `date`, `year`, `week`, `wday` (1 for Monday) and `y`; for the others
#'   `grp`, `prd` and `y`, with `band` (`plothrz()`) or `rank`
#'   (`plotbump()`).
#' @examples
#' days <- seq(as.Date("2024-01-01"), as.Date("2024-03-31"), by = "day")
#' plotcal(days, seq_along(days) %% 7)
#' d <- expand.grid(year = 2015:2020, area = c("n", "s", "e"))
#' d$v <- c(1:6, 6:1, rep(3, 6))
#' plothrz(d$v, d$year, d$area, bands = 2)
#' plotmult(d$v, d$year, d$area)
#' plotbump(d$v, d$year, d$area)
#' @export
plotcal <- function(x, y, ..., main = NULL, col = NULL, pal = NULL) {
  # before R 4.3 a number needs the origin given; text and dates ignore it
  dt <- tryCatch(as.Date(x, origin = "1970-01-01"), error = function(e) NULL)
  if (is.null(dt) || !is.numeric(y) || length(dt) != length(y) ||
        anyNA(dt[!is.na(x)])) {
    stop("plotcal: `x` must be dates, one per value of the numbers `y`.",
         call. = FALSE)
  }
  ok <- !is.na(dt)
  d <- data.frame(date = dt[ok], y = as.numeric(y[ok]))
  d$year <- as.integer(format(d$date, "%Y"))
  d$wday <- as.integer(format(d$date, "%u"))
  j1 <- as.Date(paste0(d$year, "-01-01"))
  d$week <- as.integer((as.numeric(d$date - j1) +
                          as.integer(format(j1, "%u")) - 1) %/% 7)
  yrs <- sort(unique(d$year))
  sh <- if (is.null(col)) {
    blpal(.palnm(if (is.null(pal)) "sequential" else pal, "plotcal"), 9L)
  } else {
    .pcols(9L, col, NULL, "plotcal")
  }
  r <- .prange(d$y)
  lvl <- 1L + pmin(8L, floor(8 * (d$y - r[1L]) / diff(r)))
  op <- .ppar()
  on.exit(graphics::par(op), add = TRUE)
  graphics::plot(c(0, 54), c(0, 8 * length(yrs)), type = "n", axes = FALSE,
                 xlab = "", ylab = "", main = "", ...)
  for (i in seq_along(yrs)) {
    s <- d$year == yrs[i]
    top <- 8 * (length(yrs) - i + 1L)
    graphics::rect(d$week[s], top - d$wday[s], d$week[s] + 0.92,
                   top - d$wday[s] + 0.92,
                   col = ifelse(is.na(lvl[s]), "#E0E0E0", sh[lvl[s]]),
                   border = NA)
    graphics::text(-0.5, top - 4, yrs[i], adj = 1, xpd = NA, cex = 0.85)
  }
  graphics::axis(1, at = seq(0, 52, by = 13),
                 labels = c("Jan", "Apr", "Jul", "Oct", "Dec"), tick = FALSE)
  graphics::title(main = .pmain(main, sprintf(
    "Daily values, %s to %s (darker is higher)", format(r[1L]),
    format(r[2L])
  ), "plotcal"), font.main = 1, cex.main = 1.05)
  invisible(d[, c("date", "year", "week", "wday", "y")])
}

#' @rdname plotcal
#' @export
plothrz <- function(y, prd, grp = NULL, bands = 3L, ..., main = NULL,
                    col = NULL, pal = NULL) {
  d <- .pser(y, prd, grp, "plothrz")
  if (!.pnum1(bands, 1, 6, TRUE)) {
    stop("plothrz: `bands` must be one whole number from 1 to 6.",
         call. = FALSE)
  }
  lv <- levels(d$grp)
  d$dev <- d$y - stats::ave(d$y, d$grp, FUN = function(v) {
    stats::median(v, na.rm = TRUE)
  })
  top <- max(abs(d$dev), 0, na.rm = TRUE)
  if (top == 0) top <- 1
  bw <- top / bands
  up <- blpal("sequential", bands + 1L)[-1L]
  dn <- grDevices::colorRampPalette(c("#FBD5CC", "#B2182B"))(bands)
  op <- .ppar(mar = c(4.2, .plmar(lv), 3, 1))
  on.exit(graphics::par(op), add = TRUE)
  graphics::plot(range(d$pos), c(0, length(lv)), type = "n", yaxt = "n",
                 xlab = "Period", ylab = "", main = "", ...)
  for (i in seq_along(lv)) {
    s <- d[d$grp == lv[i] & is.finite(d$dev), ]
    b0 <- length(lv) - i
    for (b in seq_len(bands)) {
      v <- pmin(pmax(abs(s$dev) - (b - 1) * bw, 0), bw) / bw * 0.9
      graphics::polygon(c(s$pos, rev(s$pos)),
                        b0 + c(ifelse(s$dev > 0, v, 0), rep(0, nrow(s))),
                        border = NA, col = up[b])
      graphics::polygon(c(s$pos, rev(s$pos)),
                        b0 + c(ifelse(s$dev < 0, v, 0), rep(0, nrow(s))),
                        border = NA, col = dn[b], density = 40)
    }
  }
  graphics::axis(2, at = length(lv) - seq_along(lv) + 0.45, labels = lv,
                 tick = FALSE)
  graphics::title(main = .pmain(main, sprintf(paste0(
    "Distance from each series' median, %d bands (above solid, below ",
    "hatched)"), bands), "plothrz"), font.main = 1, cex.main = 1)
  d$band <- pmin(bands, ceiling(abs(d$dev) / bw))
  invisible(d[, c("grp", "prd", "y", "band")])
}

#' @rdname plotcal
#' @export
plotmult <- function(y, prd, grp = NULL, ..., main = NULL, col = NULL,
                     pal = NULL) {
  d <- .pser(y, prd, grp, "plotmult")
  lv <- levels(d$grp)
  cl <- .pcols(length(lv), col, pal, "plotmult")
  main <- .pmain(main, "", "plotmult")
  op <- .ppanels(length(lv))
  on.exit(graphics::par(op), add = TRUE)
  yl <- .prange(d$y)
  for (i in seq_along(lv)) {
    s <- d[d$grp == lv[i], ]
    .pfrm(range(d$pos), yl, lv[i], ...)
    graphics::lines(s$pos, s$y, col = cl[i], lwd = 2, type = "o", pch = 19)
  }
  if (nzchar(main)) graphics::mtext(main, side = 3, line = -1.2, outer = TRUE)
  invisible(d[, c("grp", "prd", "y")])
}

#' @rdname plotcal
#' @export
plotbump <- function(y, prd, grp, desc = TRUE, ..., main = NULL, col = NULL,
                     pal = NULL) {
  d <- .pser(y, prd, grp, "plotbump")
  if (!isTRUE(desc) && !isFALSE(desc)) {
    stop("plotbump: `desc` must be TRUE or FALSE.", call. = FALSE)
  }
  d$rank <- stats::ave(if (desc) -d$y else d$y, d$pos, FUN = function(v) {
    rank(v, ties.method = "min")
  })
  lv <- levels(d$grp)
  cl <- .pcols(length(lv), col, pal, "plotbump")
  op <- .ppar(mar = c(4.2, 4.6, 3, 6))
  on.exit(graphics::par(op), add = TRUE)
  nr <- max(d$rank, na.rm = TRUE)
  .pfrm(range(d$pos), c(nr + 0.5, 0.5),
          .pmain(main, "Ranks over time", "plotbump"), "Period", "Rank",
          yaxt = "n", ...)
  graphics::axis(2, at = seq_len(nr), tick = FALSE)
  for (i in seq_along(lv)) {
    s <- d[d$grp == lv[i], ]
    graphics::lines(s$pos, s$rank, col = cl[i], lwd = 3, lty = i)
    graphics::points(s$pos, s$rank, pch = 20 + (i - 1L) %% 5L, bg = cl[i],
                     col = "white", cex = 1.4)
    graphics::text(s$pos[nrow(s)], s$rank[nrow(s)], lv[i], pos = 4,
                   xpd = NA, cex = 0.8, col = cl[i])
  }
  invisible(d[, c("grp", "prd", "y", "rank")])
}

#' Event studies for several outcomes in one grid
#'
#' Draws the event-study path of each result in a named list, one panel
#' each, as `plot()` draws one: estimates by time since adoption with their
#' intervals, before adoption in open circles and after in filled ones.
#'
#' @param x A named list of results of [evtst()] or [stgdid()].
#' @param ... Passed to the base-graphics call that draws each panel.
#' @param main `NULL` or one string, an overall title.
#' @param col,pal As for [plotlrnz()].
#' @return Invisibly, a data frame of the points drawn: `panel`, `t`, `est`,
#'   `lo` and `hi`.
#' @examples
#' d <- expand.grid(year = 2015:2024, town = LETTERS[1:12])
#' i <- match(d$town, LETTERS)
#' first <- rep(c(2019, 2021, NA), each = 4)
#' d$first <- first[i]
#' d$rate <- 30 + i + 0.5 * (d$year - 2015) -
#'   4 * (!is.na(d$first) & d$year >= d$first) + sin(seq_len(120))
#' d$cases <- d$rate / 2 + cos(seq_len(120))
#' plotevgrid(list(
#'   Rate = evtst(d$rate, d$town, d$year, d$first, leads = 3, lags = 3),
#'   Cases = evtst(d$cases, d$town, d$year, d$first, leads = 3, lags = 3)
#' ))
#' @export
plotevgrid <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  if (!.pevok(x)) {
    stop("plotevgrid: `x` must be a named list of evtst or stgdid results.",
         call. = FALSE)
  }
  main <- .pmain(main, "", "plotevgrid")
  .pcols(2L, col, pal, "plotevgrid")
  op <- .ppanels(length(x))
  on.exit(graphics::par(op), add = TRUE)
  out <- do.call(rbind, lapply(names(x), function(k) {
    e <- .pev(x[[k]])
    r <- .ppath(e$evtm, e$est, e$lo, e$hi, k, "Periods since adoption",
                "Effect", col, pal, ...)
    cbind(panel = k, r, stringsAsFactors = FALSE)
  }))
  if (nzchar(main)) graphics::mtext(main, side = 3, line = -1.2, outer = TRUE)
  invisible(out)
}

# ---- checks and computations for the helpers --------------------------------

.pxg <- function(x, grp, fn) {
  if (!is.numeric(x) || is.object(x)) {
    stop(fn, ": `x` must be a numeric vector.", call. = FALSE)
  }
  if (is.null(grp)) grp <- rep("all", length(x))
  if (length(grp) != length(x)) {
    stop(fn, ": `grp` must have one entry per value of `x`.", call. = FALSE)
  }
  d <- data.frame(x = as.numeric(x), grp = as.character(grp),
                  stringsAsFactors = FALSE)
  d <- d[is.finite(d$x) & !is.na(d$grp), , drop = FALSE]
  if (nrow(d) < 2L) {
    stop(fn, ": `x` needs at least two finite values.", call. = FALSE)
  }
  d$grp <- factor(d$grp, levels = unique(d$grp))
  d
}

.pbw <- function(bw, fn) {
  rules <- c("nrd0", "nrd", "ucv", "bcv", "sj", "sj-ste", "sj-dpi")
  ok <- if (is.character(bw)) {
    length(bw) == 1L && !is.na(bw) && bw %in% rules
  } else {
    is.numeric(bw) && length(bw) == 1L && is.finite(bw) && bw > 0
  }
  if (!ok) {
    stop(fn, ": `bw` must be one positive number or one of ",
         toString(rules), ".", call. = FALSE)
  }
  invisible(bw)
}

# stats::density(x, bw, adjust, kernel, weights, n, from, to, cut) by the
# compiled kernel of the command line's plotrdg() and plotvln(): a list of
# the points x, the densities y, the bandwidth bw, the points n of the FFT
# grid and warn, TRUE when a ucv or bcv minimum is at one end of its range.
# from and to NA: not given. An unknown rule or kernel is an error.
.dens <- function(x, bw = "nrd0", adjust = 1, kernel = "gaussian",
                  weights = NULL, n = 512L, from = NA, to = NA, cut = 3) {
  rs <- c("nrd0", "nrd", "ucv", "bcv", "sj-ste", "sj-dpi", "sj")
  ks <- c("gaussian", "epanechnikov", "rectangular", "triangular",
          "biweight", "cosine", "optcosine")
  rule <- if (is.character(bw)) c(1:6, 5L)[match(tolower(bw), rs)] else 0L
  if (!is.null(weights)) weights <- as.double(weights)
  .Call(C_clab_r_dens, as.double(x), weights,
        as.integer(c(rule, match(kernel, ks) - 1L, n)),
        as.double(c(if (is.character(bw)) 0 else bw, adjust, cut, from, to)))
}

# the density of one group for plotrdg() and plotvln(): the compiled kernel,
# so R, the command line and the other languages draw the same curve on
# every version of R, with the warning bw.ucv() and bw.bcv() give and a
# refusal in the words of stats::density(), under no other name
.pdens <- function(v, bw) {
  k <- tryCatch(.dens(v, bw = bw), error = function(e) {
    stop(sub("^dens: ", "", conditionMessage(e)), call. = FALSE)
  })
  if (k$warn) {
    warning("minimum occurred at one end of the range", call. = FALSE)
  }
  k
}

# the Lorenz curve and Gini of non-negative values, from the compiled core
.pcc <- function(x, fn) {
  ok <- is.numeric(x) && !is.object(x)
  if (!ok || !length(x) || any(!is.finite(x) | x < 0) || sum(x) <= 0) {
    stop(fn, ": `x` must be finite values from 0 with a positive total.",
         call. = FALSE)
  }
  cconc(as.numeric(x), nsim = 0)
}

.plz <- function(d) {
  if (!all(c("pop", "share") %in% names(d))) {
    stop("plotlrnz: a frame `x` needs the columns pop and share.",
         call. = FALSE)
  }
  d <- d[order(d$pop, d$share), c("pop", "share"), drop = FALSE]
  v <- c(d$pop, d$share)
  if (!nrow(d) || anyNA(v) || any(v < 0 | v > 1) || any(diff(d$share) < 0)) {
    stop("plotlrnz: `x` must be a Lorenz curve: shares in [0, 1], ",
         "non-decreasing.", call. = FALSE)
  }
  if (d$pop[1L] != 0 || d$share[1L] != 0) {
    d <- rbind(data.frame(pop = 0, share = 0), d)
  }
  rownames(d) <- NULL
  d
}

.phex <- function(x, y, bins, rx, ry) {
  sx <- (x - rx[1L]) / diff(rx) * bins
  sy <- (y - ry[1L]) / diff(ry) * bins
  h <- sqrt(3) / 2
  r1 <- round(sy / h)
  c1 <- round(sx - (r1 %% 2) / 2)
  r2 <- floor(sy / h)
  r2 <- ifelse(r2 == r1, r2 + 1, r2)
  c2 <- round(sx - (r2 %% 2) / 2)
  d1 <- (sx - c1 - (r1 %% 2) / 2)^2 + (sy - r1 * h)^2
  d2 <- (sx - c2 - (r2 %% 2) / 2)^2 + (sy - r2 * h)^2
  rr <- ifelse(d1 <= d2, r1, r2)
  cc <- ifelse(d1 <= d2, c1, c2)
  # the hexagons in order of row, then column, the same in every locale
  key <- paste(rr, cc)
  f <- which(!duplicated(key))
  f <- f[order(rr[f], cc[f])]
  data.frame(x = rx[1L] + (cc[f] + (rr[f] %% 2) / 2) / bins * diff(rx),
             y = ry[1L] + rr[f] * h / bins * diff(ry),
             count = tabulate(match(key, key[f]), length(f)))
}

.pser <- function(y, prd, grp, fn) {
  if (!is.numeric(y) || is.object(y)) {
    stop(fn, ": `y` must be a numeric vector.", call. = FALSE)
  }
  if (is.null(grp)) grp <- rep("series", length(y))
  if (length(prd) != length(y) || length(grp) != length(y)) {
    stop(fn, ": `prd` and `grp` must have one entry per value of `y`.",
         call. = FALSE)
  }
  d <- data.frame(grp = as.character(grp), prd = prd, y = as.numeric(y),
                  stringsAsFactors = FALSE)
  d <- d[!is.na(d$grp) & !is.na(d$prd), , drop = FALSE]
  if (!nrow(d)) {
    stop(fn, ": no value has both a period and a series.", call. = FALSE)
  }
  d$grp <- factor(d$grp, levels = unique(d$grp))
  d$pos <- if (is.numeric(d$prd) || inherits(d$prd, "Date")) {
    as.numeric(d$prd)
  } else if (is.factor(d$prd)) {
    as.integer(d$prd)
  } else {
    match(d$prd, unique(d$prd))
  }
  d[order(d$grp, d$pos), , drop = FALSE]
}

# the checked inputs of plotfnl(): areas missing either count dropped
.pfchk <- function(obs, expd, area, lvls, target) {
  for (a in list(list(obs, "obs"), list(expd, "expd"), list(lvls, "lvls"))) {
    if (!is.numeric(a[[1L]])) {
      stop("plotfnl: `", a[[2L]], "` must be numbers.", call. = FALSE)
    }
  }
  n <- length(obs)
  if (n != length(expd) || !n) {
    stop("plotfnl: `obs` and `expd` must be the same length, from 1.",
         call. = FALSE)
  }
  if (!is.null(area) && length(area) != n) {
    stop("plotfnl: `area` must have one label per area.", call. = FALSE)
  }
  keep <- !is.na(obs) & !is.na(expd)
  obs <- as.numeric(obs[keep])
  expd <- as.numeric(expd[keep])
  if (any(c(!length(obs), !is.finite(c(obs, expd)), obs < 0, expd <= 0))) {
    stop("plotfnl: the counts must be finite, `obs` from 0 and `expd` ",
         "positive.", call. = FALSE)
  }
  if (any(c(!length(lvls), is.na(lvls), lvls <= 0, lvls >= 1))) {
    stop("plotfnl: `lvls` must be probabilities strictly between 0 and 1.",
         call. = FALSE)
  }
  if (!.pnum1(target, .Machine$double.xmin)) {
    stop("plotfnl: `target` must be one positive number.", call. = FALSE)
  }
  list(obs = obs, expd = expd, keep = keep, lv = sort(as.numeric(lvls)))
}

# a named list of event studies
.pevok <- function(x) {
  ok <- c(is.list(x), !is.object(x), length(x) > 0L, !is.null(names(x)))
  all(ok) && all(nzchar(names(x)), vapply(x, inherits, NA,
                                          c("evtst", "stgdid")))
}
