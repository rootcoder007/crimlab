#' Plot styling: palettes and session options
#'
#' Every plot in the package draws with base graphics and one small set of
#' style choices, which these functions expose. `blpal()` returns
#' colour-blind-safe palettes: `"okabe"` (the Okabe-Ito categorical set),
#' `"diverging"` (blue, grey, orange-red), `"sequential"` (light to dark
#' blue) and `"grey"`. `plotopt()` sets, for the session, the defaults every
#' plot uses: `palette`, `grid`, `cex`, `las`, `bg` and `panel`. Each plot
#' also takes `main`, `col` (which overrides the palette), `pal` and `...`
#' (passed to the base-graphics call), so one plot can change without
#' changing the session.
#'
#' Colour never carries meaning alone: series that differ in colour differ
#' in line type or point shape too, so a plot reads the same in grey print.
#'
#' @details `"okabe"` is the eight-colour set of Okabe and Ito (2008),
#'   distinguishable under the common forms of colour-vision deficiency.
#'   `"diverging"` has three anchors (blue `#2166AC`, grey `#8C8C8C`,
#'   orange-red `#D6604D`) and is used for changes and residuals;
#'   `"sequential"` has five blues from light to dark, and `"grey"` five
#'   greys, black first. With `n` given, `"okabe"` and `"grey"` repeat
#'   their colours and `"diverging"` and `"sequential"` are interpolated.
#'
#'   The session options are kept in `options(crimlab.plot)`; the defaults
#'   are `palette = "grey"`, `grid = FALSE`, `cex = 1`, `las = 1`,
#'   `bg = "white"` and `panel = NA` (the fill of the plot region; `NA` for
#'   none): black and grey on white, as a journal prints it. For colour,
#'   `plotopt(palette = "okabe")`; for a tinted panel behind white grid
#'   lines, `plotopt(grid = TRUE, panel =` `"#EEF1F6")`. Setting some options
#'   keeps the others; an unusable value set by hand with `options()` falls
#'   back to its default.
#'
#' @param name The palette: one of `"okabe"` (the default), `"diverging"`,
#'   `"sequential"` and `"grey"`.
#' @param n `NULL` (the palette's own colours) or the number of colours
#'   wanted: one whole number from 1 to 67108864 (2^26).
#' @param ... For `plotopt()`, named options to set (see Details). With
#'   none, the options in force are returned.
#' @return `blpal()`: a character vector of colours. `plotopt()`: the named
#'   list of the six options in force, invisibly when setting.
#' @references Okabe, M. and Ito, K. (2008). Color Universal Design (CUD):
#'   how to make figures and presentations that are friendly to colorblind
#'   people. <https://jfly.uni-koeln.de/color/>
#' @srrstats {EA5.0} All plots share one base-graphics style set through
#'   plotopt(), with colour-blind-safe palettes from blpal().
#' @srrstats {EA5.0b} The default palette is black and grey, with Okabe-Ito
#'   for colour; series differ by line type or point shape as well as colour.
#' @examples
#' blpal("okabe", 3)
#' old <- plotopt()
#' plotopt(palette = "okabe", grid = TRUE)
#' plotopt(palette = old$palette, grid = old$grid)
#' plot(1:8, pch = 19, cex = 3, col = blpal())
#' @export
blpal <- function(name = "okabe", n = NULL) {
  base <- .pals[[.palnm(name, "blpal")]]
  if (is.null(n)) {
    return(base)
  }
  if (!.pnum1(n, 1, Inf, TRUE)) {
    stop("blpal: `n` must be one positive whole number.", call. = FALSE)
  }
  if (n > 2^26) stop("blpal: `n` must be at most 67108864.", call. = FALSE)
  if (name %in% c("okabe", "grey")) {
    return(rep_len(base, n))
  }
  if (n == 1) {
    return(base[1L])
  }
  grDevices::colorRampPalette(base)(n)
}

#' @rdname blpal
#' @export
plotopt <- function(...) {
  cur <- .popt()
  new <- list(...)
  if (!length(new)) {
    return(cur)
  }
  nm <- names(new)
  if (is.null(nm) || !all(nm %in% names(cur))) {
    stop("plotopt: the options are ", toString(names(cur)), ".",
         call. = FALSE)
  }
  for (k in nm) {
    if (k == "palette") .palnm(new[[k]], "plotopt")
    if (!.pok(k, new[[k]])) {
      stop("plotopt: `", k, "` must be ", .pmsg[[k]], ".", call. = FALSE)
    }
  }
  cur[nm] <- new
  options(crimlab.plot = cur)
  invisible(cur)
}

#' @rdname blpal
#' @details `yoypals()` lists the palettes of the year-over-year change
#'   tables: `"diverging"` (orange-red for a rise judged bad, teal for one
#'   judged good), `"safe"` (orange and blue, safe for red-green colour
#'   blindness) and `"mono"` (black and grey).
#' @return `yoypals()`: the names of the change-table palettes.
#' @export
yoypals <- function() names(.YOY_PALETTES)

# the colours of each change-table palette: text for a bad, good and flat
# change, and the backgrounds of bad and good rows
.YOY_PALETTES <- list(
  diverging = list(bad = "#c2492f", good = "#1b7f6b", flat = "#6b6b6b",
                   badbg = "#fbeae6", goodbg = "#e6f4f1"),
  safe = list(bad = "#b8621b", good = "#18567d", flat = "#6b6b6b",
              badbg = "#fdf0e4", goodbg = "#e7eff6"),
  mono = list(bad = "#1a1a1a", good = "#1a1a1a", flat = "#6b6b6b",
              badbg = "#f0f0f0", goodbg = "#f7f7f7")
)

.pals <- list(
  okabe = c("#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7",
            "#56B4E9", "#F0E442", "#000000"),
  diverging = c("#2166AC", "#8C8C8C", "#D6604D"),
  sequential = c("#C6DBEF", "#6BAED6", "#2171B5", "#08306B", "#041E42"),
  grey = c("#000000", "#636363", "#969696", "#3D3D3D", "#B0B0B0")
)

.pmsg <- list(grid = "TRUE or FALSE", cex = "one positive number",
              las = "0, 1, 2 or 3", bg = "one colour",
              panel = "one colour, or NA")

#' The name of a palette checked.
#' @noRd
.palnm <- function(name, fn) {
  if (!is.character(name) || length(name) != 1L ||
        !name %in% names(.pals)) {
    stop(fn, ": the palette must be one of ", toString(names(.pals)), ".",
         call. = FALSE)
  }
  name
}

#' one colour: a name, a hex code or a palette index; an index past the
#' range of an integer, or text that is not valid in its encoding, is
#' refused before col2rgb() warns or fails on it
#' @noRd
.iscol <- function(v) {
  if (length(v) != 1L || is.na(v) || !(is.numeric(v) || is.character(v))) {
    return(FALSE)
  }
  if (is.numeric(v) && !(v >= 0 && v <= .Machine$integer.max)) {
    return(FALSE)
  }
  if (is.character(v) && !validEnc(v)) {
    return(FALSE)
  }
  !is.null(tryCatch(grDevices::col2rgb(v), error = function(e) NULL))
}

#' Whether the value `v` is a valid setting of the plot option `k`.
#' @noRd
.pok <- function(k, v) {
  switch(
    k,
    palette = is.character(v) && length(v) == 1L && v %in% names(.pals),
    grid = isTRUE(v) || isFALSE(v),
    cex = .pnum1(v, .Machine$double.xmin),
    las = .pnum1(v, 0, 3, TRUE),
    bg = .iscol(v),
    panel = (length(v) == 1L && is.na(v)) || .iscol(v)
  )
}

#' one finite number in \[lo, hi\], a whole number when `int`
#' @noRd
.pnum1 <- function(v, lo = -Inf, hi = Inf, int = FALSE) {
  is.numeric(v) && length(v) == 1L && is.null(dim(v)) &&
    all(c(is.finite(v), v >= lo, v <= hi, !int | v == round(v)))
}

#' the options in force: the defaults, overridden by usable set values
#' @noRd
.popt <- function() {
  d <- list(palette = "grey", grid = FALSE, cex = 1, las = 1, bg = "white",
            panel = NA)
  o <- getOption("crimlab.plot")
  for (k in intersect(names(o), names(d))) {
    if (.pok(k, o[[k]])) d[[k]] <- o[[k]]
  }
  d
}

# ---- the drawing kit every plot uses ----------------------------------------

#' n colours: `col` as given (recycled) or the palette's
#' @noRd
.pcols <- function(n, col = NULL, pal = NULL, fn = "plot") {
  if (!is.null(col)) {
    if (!is.atomic(col) || !length(col) || !all(vapply(col, .iscol, NA))) {
      stop(fn, ": `col` must be colours.", call. = FALSE)
    }
    return(rep_len(col, n))
  }
  blpal(.palnm(if (is.null(pal)) .popt()$palette else pal, fn), n)
}

#' The title of a plot: `main` when given (one string), else the default.
#' @noRd
.pmain <- function(main, def, fn = "plot") {
  if (is.null(main)) {
    return(def)
  }
  if (!is.character(main) || length(main) != 1L || is.na(main)) {
    stop(fn, ": `main` must be one string or NULL.", call. = FALSE)
  }
  .cktxt(main, "main", fn)
}

#' sets the session look and the margins; returns the old values to restore
#' @noRd
.ppar <- function(mar = c(4.2, 4.6, 3, 1), ...) {
  o <- .popt()
  graphics::par(bg = o$bg, las = o$las, cex = o$cex, mar = mar,
                mgp = c(2.6, 0.7, 0), ...)
}

#' the left margin that fits row labels
#' @noRd
.plmar <- function(lab) {
  max(4.6, min(18, 0.55 * max(nchar(as.character(lab)), 1L) + 1.5))
}

#' an empty frame: panel fill, grid, title and axis titles
#' @noRd
.pfrm <- function(xlim, ylim, main, xlab = "", ylab = "", ...) {
  o <- .popt()
  graphics::plot(xlim, ylim, type = "n", main = "", xlab = xlab,
                 ylab = ylab, ...)
  u <- graphics::par("usr")
  lx <- graphics::par("xlog")
  ly <- graphics::par("ylog")
  if (!is.na(o$panel)) {
    graphics::rect(if (lx) 10^u[1L] else u[1L], if (ly) 10^u[3L] else u[3L],
                   if (lx) 10^u[2L] else u[2L], if (ly) 10^u[4L] else u[4L],
                   col = o$panel, border = NA)
  }
  if (o$grid) {
    gc <- if (is.na(o$panel)) "#DDDDDD" else "white"
    graphics::abline(h = graphics::axTicks(2), v = graphics::axTicks(1),
                     col = gc, lwd = 0.8)
  }
  graphics::box()
  graphics::title(main = main, font.main = 2, cex.main = 1.05)
  invisible(NULL)
}

#' The range of the finite values given, widened when it is a point, and 0 to
#'   1 when there is none.
#' @noRd
.prange <- function(...) {
  r <- suppressWarnings(range(c(...), finite = TRUE))
  if (!all(is.finite(r))) r <- c(0, 1)
  if (diff(r) == 0) r <- r + c(-0.5, 0.5)
  r
}

#' a grid of n panels; returns the old par to restore
#' @noRd
.ppanels <- function(n) {
  nc <- ceiling(sqrt(n))
  # setting mfrow resets cex: keep the values before it, mfrow first so that
  # restoring it does not undo the cex restored after it
  op <- graphics::par(c("mfrow", "bg", "las", "cex", "mar", "mgp"))
  graphics::par(mfrow = c(ceiling(n / nc), nc))
  .ppar(mar = c(4, 4.2, 2.6, 0.8))
  op
}

#' a legend that pairs every colour with a line type or point shape
#' @noRd
.pkey <- function(where, lab, col, lty = NA, pch = NA, ...) {
  a <- list(where, legend = lab, col = col, bty = "n", cex = 0.85, ...)
  if (!all(is.na(lty))) a$lty <- lty
  if (!all(is.na(pch))) a$pch <- pch
  do.call(graphics::legend, a)
}

#' point shapes, one per group, that also tell the groups apart in grey
#' @noRd
.ppch <- function(g) c(19, 17, 15, 1, 2, 0)[(as.integer(g) - 1L) %% 6L + 1L]

#' a dot-and-interval (forest) chart, one row per label from the top; `grp`
#' gives each group its colour and shape (several groups may share a row),
#' named in a key unless `key` is FALSE (the row labels then say it);
#' returns the points drawn
#' @noRd
.pdots <- function(lab, est, lo, hi, main, xlab, col = NULL, pal = NULL,
                   ref = NULL, grp = NULL, log = "", key = TRUE, ...) {
  lab <- as.character(lab)
  n <- length(lab)
  lo <- rep_len(as.numeric(lo), n)
  hi <- rep_len(as.numeric(hi), n)
  g <- factor(if (is.null(grp)) rep("", n) else grp,
              levels = unique(if (is.null(grp)) "" else grp))
  cl <- .pcols(nlevels(g), col, pal)
  u <- unique(lab)
  k <- length(u)
  y <- k + 1L - match(lab, u)
  op <- .ppar(mar = c(4.2, .plmar(u), 3, 1))
  on.exit(graphics::par(op), add = TRUE)
  v <- c(est, lo, hi, ref)
  if (log == "x") v <- v[v > 0]
  # rows of room above the first row for the key, so that it covers no
  # point: its height in inches against that of the plot region
  top <- 0
  if (!is.null(grp) && key) {
    kh <- nlevels(g) * graphics::par("cin")[2L] * graphics::par("cex") * 0.95
    ph <- graphics::par("fin")[2L] - sum(graphics::par("mai")[c(1L, 3L)])
    if (ph > kh) top <- kh * (k + 1) / (ph - kh)
  }
  .pfrm(.prange(v), c(0.5, k + 0.5 + top), main, xlab, "", yaxt = "n",
          log = log, ...)
  # every row keeps its label (axis() would drop one that comes close to
  # its neighbour), in letters small enough for the height of a row
  rh <- graphics::par("pin")[2L] / (k + top)
  ch <- graphics::par("cin")[2L] * graphics::par("cex")
  graphics::axis(2, at = rev(seq_len(k)), labels = u, tick = FALSE, las = 1,
                 gap.axis = -1, cex.axis = max(0.5, min(1, 0.85 * rh / ch)))
  if (!is.null(ref)) graphics::abline(v = ref, lty = 2, col = "#636363")
  w <- is.finite(lo) & is.finite(hi)
  .pint(lo[w], y[w], hi[w], y[w], cl[g][w], 2)
  graphics::points(est, y, pch = .ppch(g), col = cl[g], cex = 1.2)
  if (!is.null(grp) && key) {
    .pkey("topright", levels(g), cl, lty = if (any(w)) 1 else NA,
          pch = .ppch(seq_len(nlevels(g))))
  }
  d <- data.frame(label = lab, est = est, lo = lo, hi = hi,
                  stringsAsFactors = FALSE)
  if (!is.null(grp)) d$group <- as.character(g)
  invisible(d)
}

#' one of the values a plot option takes
#' @noRd
.pany <- function(v, ok, nm) {
  if (!is.character(v) || length(v) != 1L || is.na(v) || !v %in% ok) {
    stop("plot: `", nm, "` must be one of ", toString(dQuote(ok, FALSE)),
         ".", call. = FALSE)
  }
  v
}

#' a histogram (or a density with its 5th and 95th percentiles dashed) of
#' `z` with a red line at `obs`, and lines of text above the bars
#' @noRd
.phist <- function(z, obs, main, xlab, type, fill, note = NULL, ...) {
  op <- .ppar()
  on.exit(graphics::par(op), add = TRUE)
  if (type == "density" && length(unique(z)) > 1L) {
    k <- stats::density(z)
    .pfrm(.prange(k$x, obs), .prange(0, k$y), main, xlab, "Density", ...)
    graphics::lines(k$x, k$y, col = if (fill == "#D5D5D5") "#000000" else fill,
                    lwd = 2)
    graphics::abline(v = stats::quantile(z, c(0.05, 0.95)), lty = 2,
                     col = "#636363")
  } else {
    h <- graphics::hist(z, breaks = min(30L, max(5L, length(z) %/% 3L)),
                        plot = FALSE)
    top <- max(h$counts) * if (is.null(note)) 1.05 else 1.45
    .pfrm(.prange(h$breaks, obs), c(0, top), main, xlab, "Count", ...)
    nb <- length(h$breaks)
    graphics::rect(h$breaks[-nb], 0, h$breaks[-1L], h$counts, col = fill,
                   border = "white")
  }
  graphics::abline(v = obs, col = "#C0392B", lwd = 2)
  if (!is.null(note)) {
    graphics::legend("topleft", legend = note, bty = "n", cex = 0.8)
  }
  graphics::box()
  invisible(NULL)
}

#' the fill of histogram bars: grey, or the colour the caller asked for
#' @noRd
.pfill <- function(col, pal) {
  if (is.null(col) && is.null(pal)) "#D5D5D5" else .pcols(1L, col, pal)
}

#' intervals as bars with a cap at each end; one of no length gets no cap
#' @noRd
.pint <- function(x0, y0, x1, y1, col, lwd) {
  col <- rep_len(col, length(x0))
  z <- x0 == x1 & y0 == y1
  graphics::segments(x0[z], y0[z], x1[z], y1[z], col = col[z], lwd = lwd)
  if (!all(z)) {
    graphics::arrows(x0[!z], y0[!z], x1[!z], y1[!z], length = 0.03,
                     angle = 90, code = 3, col = col[!z], lwd = lwd)
  }
}

#' estimates over time with intervals: before the reference in one colour and
#' open circles, after in another and filled; returns the rows drawn
#' @noRd
.ppath <- function(t, est, lo, hi, main, xlab, ylab, col = NULL, pal = NULL,
                   ref = 0, at = -0.5, ..., style = "bars") {
  .pany(style, c("bars", "stata"), "style")
  cl <- .pcols(2L, col, pal)
  op <- .ppar()
  on.exit(graphics::par(op), add = TRUE)
  .pfrm(.prange(t), .prange(est, lo, hi, ref), main, xlab, ylab, ...)
  graphics::abline(h = ref, lty = 2, col = "#636363")
  graphics::abline(v = at, lty = 3, col = "#636363")
  post <- t > at
  w <- is.finite(lo) & is.finite(hi) & hi > lo
  if (style == "stata") {
    # hollow points joined by a line, their intervals dashed, in one colour
    o <- order(t)
    graphics::par(lty = 2)
    .pint(t[w], lo[w], t[w], hi[w], cl[1L], 1)
    graphics::par(lty = 1)
    graphics::lines(t[o], est[o], col = cl[1L], lwd = 1.4)
    graphics::points(t, est, pch = 21, bg = .popt()$bg, col = cl[1L],
                     cex = 1.1)
    return(invisible(data.frame(t = t, est = est, lo = lo, hi = hi)))
  }
  .pint(t[w], lo[w], t[w], hi[w], cl[post + 1L][w], 1.6)
  graphics::points(t, est, pch = ifelse(post, 19, 1), col = cl[post + 1L],
                   cex = 1.2, lwd = 1.6)
  .pkey("bottomleft", c("before", "after"), cl, pch = c(1, 19))
  invisible(data.frame(t = t, est = est, lo = lo, hi = hi))
}

#' the frame of a plot over `x`, and `x` as numbers. Periods written as text
#' are the numbers they spell ("2021") or else their order ("2020/21"), and
#' the axis shows them as they were written, each at its place
#' @noRd
.pxfrm <- function(x, ylim, main, xlab, ylab, ...) {
  if (is.numeric(x)) {
    .pfrm(.prange(x), ylim, main, xlab, ylab, ...)
    return(x)
  }
  lab <- unique(as.character(x))
  at <- suppressWarnings(as.numeric(lab))
  if (anyNA(at)) at <- seq_along(lab)
  x <- at[match(as.character(x), lab)]
  .pfrm(.prange(x), ylim, main, xlab, ylab, xaxt = "n", ...)
  graphics::axis(1, at = at, labels = lab)
  x
}

#' one line per group, each with its own colour, line type and point shape,
#' and a band (a step band for steps) when `lo` and `hi` are given; returns
#' the points drawn
#' @noRd
.pgrp <- function(x, y, g, main, xlab, ylab, col = NULL, pal = NULL,
                  lo = NULL, hi = NULL, type = "o", ref = NULL,
                  where = "topright", vref = NULL, ...) {
  g <- factor(g, levels = unique(g))
  k <- nlevels(g)
  cl <- .pcols(k, col, pal)
  op <- .ppar()
  on.exit(graphics::par(op), add = TRUE)
  x0 <- x
  x <- .pxfrm(x, .prange(y, lo, hi, ref), main, xlab, ylab, ...)
  if (!is.null(ref)) graphics::abline(h = ref, lty = 2, col = "#636363")
  if (!is.null(vref)) graphics::abline(v = vref, lty = 3, col = "#636363")
  # steps and plain lines carry no point marks
  st <- type != "o"
  for (i in seq_len(k)) {
    s <- which(g == levels(g)[i])
    if (!is.null(lo) && type == "s") {
      graphics::lines(x[s], lo[s], type = "s", col = cl[i], lty = 3)
      graphics::lines(x[s], hi[s], type = "s", col = cl[i], lty = 3)
    } else if (!is.null(lo)) {
      w <- s[is.finite(lo[s]) & is.finite(hi[s])]
      graphics::polygon(c(x[w], rev(x[w])), c(lo[w], rev(hi[w])),
                        col = grDevices::adjustcolor(cl[i], 0.2),
                        border = NA)
    }
    graphics::lines(x[s], y[s], type = type, col = cl[i], lty = i,
                    pch = if (st) NA else .ppch(i), lwd = 2)
  }
  if (k > 1L) {
    .pkey(where, levels(g), cl, lty = seq_len(k),
          pch = if (st) NA else .ppch(seq_len(k)))
  }
  d <- data.frame(x = x0, series = as.character(g), y = y,
                  stringsAsFactors = FALSE)
  if (!is.null(lo)) {
    d$lo <- lo
    d$hi <- hi
  }
  invisible(d)
}

#' named series on a shared x
#' @noRd
.plines <- function(x, ys, main, xlab, ylab, col = NULL, pal = NULL,
                    lo = NULL, hi = NULL, ...) {
  ys <- as.list(ys)
  n <- length(x)
  .pgrp(rep(x, length(ys)), unlist(ys, use.names = FALSE),
        rep(names(ys), each = n), main, xlab, ylab, col, pal,
        lo = if (!is.null(lo)) unlist(lo, use.names = FALSE),
        hi = if (!is.null(hi)) unlist(hi, use.names = FALSE), ...)
}

#' horizontal bars, first row at the top; returns the rows drawn
#' @noRd
.pbars <- function(v, lab, main, xlab, col = NULL, pal = NULL, ref = NULL,
                   ...) {
  cl <- .pcols(1L, col, pal)
  op <- .ppar(mar = c(4.2, .plmar(lab), 3, 1))
  on.exit(graphics::par(op), add = TRUE)
  k <- length(v)
  y <- rev(seq_len(k))
  .pfrm(.prange(0, v, ref), c(0.4, k + 0.6), main, xlab, yaxt = "n", ...)
  graphics::axis(2, at = y, labels = lab, tick = FALSE, las = 1)
  w <- is.finite(v)
  graphics::rect(0, y[w] - 0.35, v[w], y[w] + 0.35, col = cl, border = NA)
  if (!all(w)) graphics::text(0, y[!w], "no finite value", pos = 4, cex = 0.75)
  if (!is.null(ref)) graphics::abline(v = ref, lty = 2, col = "#636363")
  invisible(data.frame(label = lab, value = v, stringsAsFactors = FALSE))
}

#' a matrix as coloured cells with its numbers printed in them, so the
#' shading is never the only reading; returns the long data
#' @noRd
.pmat <- function(m, main, xlab = "", ylab = "", div = FALSE, mid = 0,
                  lab = m, pal = NULL, ...) {
  m <- as.matrix(m)
  rl <- if (is.null(rownames(m))) seq_len(nrow(m)) else rownames(m)
  cn <- if (is.null(colnames(m))) seq_len(ncol(m)) else colnames(m)
  pn <- if (!is.null(pal)) pal else if (div) "diverging" else "sequential"
  sh <- blpal(.palnm(pn, "plot"), 9L)
  r <- .prange(m)
  if (div) r <- mid + c(-1, 1) * max(abs(r - mid))
  op <- .ppar(mar = c(4.2, .plmar(rl), 3, 1))
  on.exit(graphics::par(op), add = TRUE)
  # image() draws the first row at the bottom: flip so it reads as printed
  rv <- rev(seq_len(nrow(m)))
  f <- t(m[rv, , drop = FALSE])
  graphics::image(seq_len(ncol(m)), seq_len(nrow(m)), f, col = sh, zlim = r,
                  axes = FALSE, main = "", xlab = xlab, ylab = ylab, ...)
  graphics::title(main = main, font.main = 2, cex.main = 1.05)
  graphics::axis(1, at = seq_len(ncol(m)), labels = cn, tick = FALSE)
  graphics::axis(2, at = rv, labels = rl, tick = FALSE)
  graphics::box()
  z <- (f - r[1L]) / diff(r)
  tc <- ifelse(is.finite(z) & (z > 0.7 | (div & z < 0.3)), "white", "black")
  graphics::text(rep(seq_len(ncol(m)), nrow(m)),
                 rep(seq_len(nrow(m)), each = ncol(m)),
                 format(t(as.matrix(lab)[rv, , drop = FALSE]), digits = 2),
                 col = tc, cex = 0.8)
  invisible(data.frame(row = rep(rl, ncol(m)), col = rep(cn, each = nrow(m)),
                       value = as.vector(m), stringsAsFactors = FALSE))
}

#' a surface on a grid, shaded and with contour lines
#' @noRd
.psurf <- function(x, y, z, main, xlab = "x", ylab = "y", pal = NULL, ...) {
  op <- .ppar()
  on.exit(graphics::par(op), add = TRUE)
  sh <- blpal(.palnm(if (is.null(pal)) "sequential" else pal, "plot"), 9L)
  graphics::image(x, y, z, col = sh, main = "", xlab = xlab, ylab = ylab,
                  asp = 1, ...)
  graphics::contour(x, y, z, add = TRUE, drawlabels = FALSE,
                    col = "#636363", nlevels = 6)
  graphics::title(main = main, font.main = 2, cex.main = 1.05)
  invisible(list(x = x, y = y, z = z))
}

#' the result object holds the parts a plot needs, or the plot refuses it
#' @noRd
.pck <- function(x, cls, need) {
  ok <- is.list(x) && all(need %in% names(x)) &&
    all(vapply(need, function(k) !is.null(x[[k]]), NA))
  if (!ok) {
    stop("plot: `x` is not a complete ", cls, " result (it needs ",
         toString(need), ").", call. = FALSE)
  }
  invisible(x)
}
