#' Choropleth maps in base graphics
#'
#' `chomap()` shades each area of a map by the class of its value, with no
#' spatial package: an area is a matrix of coordinates. The class breaks,
#' the class of each value and the bounding box are computed by the
#' compiled kernel; R only draws. Areas are drawn with `polypath()`, so a
#' ring inside another is a hole.
#'
#' @details With `brk` a number of classes, `meth = "quantile"` puts the
#'   breaks at R's type-7 quantiles of the finite values (so the classes
#'   hold about equal counts; breaks that repeat are merged, which leaves
#'   fewer classes) and `meth = "equal"` cuts their range into classes of
#'   equal width, as `seq(min, max, length.out = brk + 1)`. A value falls
#'   in the class whose interval holds it, intervals closed on the left
#'   and the last closed on both sides, as `findInterval(v, breaks,
#'   rightmost.closed = TRUE, all.inside = TRUE)` gives. When every finite
#'   value is the same, the one class runs from 0.5 below it to 0.5 above.
#'   Values that are missing or not finite are left unshaded.
#'
#'   The colours are a sequential ramp (light for low classes, dark for
#'   high), so the order of the classes reads in grey print too.
#'
#' @param polys A list of areas, each a numeric matrix of two columns (x
#'   and y) with at least three points; the rings of one area are separated
#'   by a row of `NA`. Names, when given, label the areas.
#' @param v Numbers, one per area.
#' @param brk The number of classes (a whole number from 1 to 9), or the
#'   breaks themselves (increasing numbers covering the finite values of
#'   `v`).
#' @param meth How a number of classes is cut: `"quantile"` or `"equal"`.
#' @param ... Passed to [graphics::plot.window()].
#' @param main `NULL` (a default title) or one string.
#' @param col `NULL` (the palette) or colours, recycled over the classes.
#' @param pal `NULL` (`"sequential"`) or a palette name for [blpal()].
#' @return Invisibly, a list: `areas`, a data frame with `area` (the label
#'   or the position of each area), `v`, `cls` (its class from 1, `NA` when
#'   unshaded) and `col` (its fill, `NA` when unshaded); `breaks`; and
#'   `bbox`, the box drawn (`xmin`, `xmax`, `ymin`, `ymax`).
#' @examples
#' sq <- function(x0, y0) cbind(x0 + c(0, 1, 1, 0), y0 + c(0, 0, 1, 1))
#' p <- list(a = sq(0, 0), b = sq(1, 0), c = sq(0, 1), d = sq(1, 1))
#' chomap(p, c(12, 3, 7, 25), brk = 2)
#' chomap(p, c(12, 3, 7, 25), brk = c(0, 5, 10, 30), main = "Rates")
#' @export
chomap <- function(polys, v, brk = 5, meth = "quantile", ..., main = NULL,
                   col = NULL, pal = NULL) {
  f <- "chomap"
  p <- .mpoly(polys, f)
  if (!is.numeric(v) || length(v) != length(polys)) {
    .ckstop(f, "`v` must be numbers, one per area")
  }
  v <- as.numeric(v)
  b <- .mbrk(v, brk, meth, f)
  fl <- .mfill(b, col, pal, f)
  bb <- .Call(C_clab_r_mapbb, f, p$x, p$y)
  mn <- .pmain(main, "Choropleth", f)
  op <- .ppar(mar = c(0.5, 0.5, 3, 0.5))
  on.exit(graphics::par(op), add = TRUE)
  .mdraw(polys, fl$col, b$br, fl$pal, bb, mn, TRUE, ...)
  invisible(list(
    areas = data.frame(area = p$lab, v = v, cls = b$cls, col = fl$col,
                       stringsAsFactors = FALSE),
    breaks = b$br, bbox = bb
  ))
}

# The areas checked: their coordinates in one vector each (an NA between
# areas, for the box) and their labels.
.mpoly <- function(polys, fn) {
  if (!is.list(polys) || is.data.frame(polys) || !length(polys)) {
    .ckstop(fn, "`polys` must be a non-empty list of areas")
  }
  for (i in seq_along(polys)) {
    m <- polys[[i]]
    ok <- is.matrix(m) && is.numeric(m) && ncol(m) == 2L &&
      !any(is.infinite(m)) && identical(is.na(m[, 1L]), is.na(m[, 2L])) &&
      sum(!is.na(m[, 1L])) >= 3L
    if (!ok) {
      .ckstop(fn, paste0("`polys[[%d]]` must be a two-column numeric ",
                         "matrix of at least three points (rows of NA ",
                         "between rings)"), i)
    }
  }
  list(x = unlist(lapply(polys, function(m) c(m[, 1L], NA)), FALSE, FALSE),
       y = unlist(lapply(polys, function(m) c(m[, 2L], NA)), FALSE, FALSE),
       lab = if (is.null(names(polys))) seq_along(polys) else names(polys))
}

# The breaks and the class of each value (NA when not finite), from the
# compiled kernel.
.mbrk <- function(v, brk, meth, fn) {
  if (!is.character(meth) || length(meth) != 1L ||
        !meth %in% c("quantile", "equal")) {
    .ckstop(fn, "`meth` must be \"quantile\" or \"equal\"")
  }
  fv <- v[is.finite(v)]
  if (!length(fv)) .ckstop(fn, "`v` has no finite value")
  one <- length(brk) == 1L
  ok <- if (one) {
    .pnum1(brk, 1, 9, TRUE)
  } else {
    is.numeric(brk) && length(brk) >= 2L && all(is.finite(brk)) &&
      all(diff(brk) > 0) && brk[1L] <= min(fv) && brk[length(brk)] >= max(fv)
  }
  if (!ok) {
    .ckstop(fn, paste0("`brk` must be a whole number of classes from 1 to ",
                       "9, or increasing breaks covering `v`"))
  }
  m <- if (!one) 2L else if (meth == "equal") 0L else 1L
  r <- .Call(C_clab_r_mapbrk, fn, v, if (one) as.integer(brk) else 0L, m,
             if (one) NULL else as.numeric(brk))
  cls <- r$cls
  cls[cls == 0L] <- NA_integer_
  list(br = r$br, cls = cls)
}

# The colour of each class and of each value.
.mfill <- function(b, col, pal, fn) {
  nc <- length(b$br) - 1L
  pl <- if (is.null(col)) {
    blpal(.palnm(if (is.null(pal)) "sequential" else pal, fn), nc)
  } else {
    .pcols(nc, col, NULL, fn)
  }
  list(pal = pl, col = pl[b$cls])
}

# An empty map frame: the box at an aspect of 1 and the title.
.mfrm <- function(bb, main, ...) {
  graphics::plot.new()
  graphics::plot.window(bb[1:2], bb[3:4], asp = 1, ...)
  graphics::title(main = main, font.main = 1, cex.main = 1.05)
}

.mdraw <- function(polys, fill, br, pl, bb, main, key, ...) {
  .mfrm(bb, main, ...)
  for (i in seq_along(polys)) {
    graphics::polypath(polys[[i]][, 1L], polys[[i]][, 2L], col = fill[i],
                       border = "#636363", rule = "evenodd")
  }
  if (key) {
    s <- formatC(signif(br, 3), format = "g")
    graphics::legend("bottomright", legend = paste(s[-length(s)], "-", s[-1L]),
                     fill = pl, bty = "n", cex = 0.8)
  }
}

#' Small multiples of a choropleth map by period
#'
#' `smmap()` draws one choropleth map per period in a grid of panels, all
#' with the same class breaks and colours, so that a change between
#' periods is a change of shade. The areas, breaks and classes are those
#' of [chomap()], the breaks computed from the values of every period
#' together; the legend is drawn in the last panel.
#'
#' @inheritParams chomap
#' @param v A numeric matrix (or a data frame of numeric columns), one row
#'   per area and one column per period; the column names label the
#'   periods.
#' @return Invisibly, a list: `areas`, a data frame with one row per area
#'   and period (the periods in order, the areas within each) holding
#'   `area`, `prd` (the period's label), `v`, `cls` and `col`; `breaks`;
#'   and `bbox`, as for [chomap()].
#' @examples
#' sq <- function(x0, y0) cbind(x0 + c(0, 1, 1, 0), y0 + c(0, 0, 1, 1))
#' p <- list(sq(0, 0), sq(1, 0), sq(0, 1), sq(1, 1))
#' v <- cbind(`2022` = c(3, 8, 5, 1), `2023` = c(4, 9, 2, 1),
#'            `2024` = c(6, 12, 2, 0))
#' smmap(p, v, brk = 3)
#' @export
smmap <- function(polys, v, brk = 5, meth = "quantile", ..., main = NULL,
                  col = NULL, pal = NULL) {
  f <- "smmap"
  p <- .mpoly(polys, f)
  if (is.data.frame(v) && length(v) && all(vapply(v, is.numeric, NA))) {
    v <- as.matrix(v)
  }
  if (!is.matrix(v) || !is.numeric(v) || nrow(v) != length(polys) ||
        !ncol(v)) {
    .ckstop(f, paste("`v` must be a numeric matrix, one row per area and",
                     "one column per period"))
  }
  prd <- colnames(v)
  if (is.null(prd)) prd <- as.character(seq_len(ncol(v)))
  vv <- as.numeric(v)
  b <- .mbrk(vv, brk, meth, f)
  fl <- .mfill(b, col, pal, f)
  bb <- .Call(C_clab_r_mapbb, f, p$x, p$y)
  mn <- .pmain(main, "Small multiples by period", f)
  op <- .ppanels(ncol(v))
  on.exit(graphics::par(op), add = TRUE)
  # restored first, so that the margins of the panels go before op's
  oo <- graphics::par(oma = c(0, 0, 2, 0), mar = c(0.5, 0.5, 2, 0.5))
  on.exit(graphics::par(oo), add = TRUE, after = FALSE)
  na <- length(polys)
  for (j in seq_along(prd)) {
    i <- (j - 1L) * na + seq_len(na)
    .mdraw(polys, fl$col[i], b$br, fl$pal, bb, prd[j], j == length(prd),
           ...)
  }
  graphics::mtext(mn, outer = TRUE, line = 0.5, cex = 1.05)
  invisible(list(
    areas = data.frame(area = rep(p$lab, length(prd)),
                       prd = rep(prd, each = na), v = vv, cls = b$cls,
                       col = fl$col, stringsAsFactors = FALSE),
    breaks = b$br, bbox = bb
  ))
}

#' Point maps, with clusters
#'
#' `ptmap()` draws points on a map at an aspect of 1, optionally coloured
#' by cluster: a result of [dbscn()] or a cluster number per point. Each
#' cluster has its own colour and point shape (so the clusters read in
#' grey print too) and its number at its centre, the mean of its points;
#' noise is drawn as grey crosses. The box and the cluster centres are
#' computed by the compiled kernel; R only draws. `plot()` of a [dbscn()]
#' result draws it through `ptmap()`, given the points it was found in.
#'
#' @param x The points: a numeric matrix or data frame of two columns (x
#'   and y), or the x coordinates with `y`. For `plot()`, a result of
#'   [dbscn()].
#' @param y `NULL` or the y coordinates, one per value of `x`. For
#'   `plot()`, the points the clusters were found in, as `x` of `ptmap()`.
#' @param cl `NULL`, a result of [dbscn()], or the cluster of each point: a
#'   whole number from 1, `NA` for noise.
#' @param ... Passed to [graphics::plot.window()].
#' @param main `NULL` (a default title) or one string.
#' @param col `NULL` (the palette) or colours, one per cluster.
#' @param pal `NULL` (the session palette, see [plotopt()]) or a palette
#'   name for [blpal()].
#' @return Invisibly, a list: `pts`, a data frame of the points (`x`, `y`,
#'   `cl` the cluster or `NA`, `col` and `pch` as drawn); `cent`, one row
#'   per cluster number from 1 to the largest (`cl`, the centre `x` and
#'   `y`, `NaN` for a number no point has, and `n` its points); and
#'   `bbox`, the box drawn (`xmin`, `xmax`, `ymin`, `ymax`).
#' @examples
#' set.seed(1)
#' p <- rbind(cbind(rnorm(30, 0, 0.2), rnorm(30, 0, 0.2)),
#'            cbind(rnorm(30, 3, 0.2), rnorm(30, 1, 0.2)),
#'            cbind(runif(10, -2, 5), runif(10, -2, 3)))
#' ptmap(p)
#' d <- dbscn(p, eps = 0.3, minpts = 5)
#' plot(d, p)
#' @export
ptmap <- function(x, y = NULL, cl = NULL, ..., main = NULL, col = NULL,
                  pal = NULL) {
  f <- "ptmap"
  p <- .mxy(x, y, f)
  n <- length(p$x)
  k <- 0L
  if (is.null(cl)) {
    cl <- rep(NA_integer_, n)
    pc <- rep(.pcols(1L, col, pal, f), n)
    ph <- rep(16L, n)
  } else {
    cl <- .mcl(cl, n, f)
    k <- max(0L, cl, na.rm = TRUE)
    pl <- .pcols(max(k, 1L), col, pal, f)
    on <- !is.na(cl)
    pc <- rep("#969696", n)
    pc[on] <- pl[cl[on]]
    ph <- rep(4L, n)
    ph[on] <- .mpch[(cl[on] - 1L) %% 8L + 1L]
  }
  bb <- .Call(C_clab_r_mapbb, f, p$x, p$y)
  ce <- .Call(C_clab_r_ptcen, p$x, p$y, ifelse(is.na(cl), 0L, cl), k)
  mn <- .pmain(main, if (k) sprintf("Points in %d clusters", k) else "Points",
               f)
  op <- .ppar(mar = c(0.5, 0.5, 3, 0.5))
  on.exit(graphics::par(op), add = TRUE)
  .mfrm(bb, mn, ...)
  graphics::points(p$x, p$y, col = pc, pch = ph)
  sh <- ce$n > 0L
  if (any(sh)) graphics::text(ce$x[sh], ce$y[sh], which(sh), font = 2)
  invisible(list(
    pts = data.frame(x = p$x, y = p$y, cl = cl, col = pc, pch = ph,
                     stringsAsFactors = FALSE),
    cent = data.frame(cl = seq_len(k), x = ce$x, y = ce$y, n = ce$n),
    bbox = bb
  ))
}

#' @rdname ptmap
#' @export
plot.dbscn <- function(x, y, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "dbscn", "cl")
  if (missing(y)) {
    stop("plot: `y` must be the points the clusters were found in.",
         call. = FALSE)
  }
  ptmap(y, cl = x, ..., main = .pmain(main, "DBSCAN clusters", "plot"),
        col = col, pal = pal)
}

# the point shapes of the clusters, in turn
.mpch <- c(16L, 17L, 15L, 18L, 1L, 2L, 0L, 5L)

# The points checked: x and y, or x a matrix or data frame of two columns.
.mxy <- function(x, y, fn) {
  if (is.null(y)) {
    ok <- (is.matrix(x) && is.numeric(x) && ncol(x) == 2L) ||
      (is.data.frame(x) && length(x) == 2L &&
         all(vapply(x, is.numeric, NA)))
    if (!ok) {
      .ckstop(fn, paste("`x` must be a two-column matrix or data frame of",
                        "coordinates, or numbers with `y`"))
    }
    y <- x[, 2L]
    x <- x[, 1L]
  } else if (!is.numeric(x) || !is.numeric(y) || length(x) != length(y)) {
    .ckstop(fn, "`x` and `y` must be numbers of the same length")
  }
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (!length(x) || !all(is.finite(c(x, y)))) {
    .ckstop(fn, "the coordinates must be finite, at least one point")
  }
  list(x = x, y = y)
}

# The cluster of each point: a dbscn result's, or whole numbers from 1 and
# NA for noise.
.mcl <- function(cl, n, fn) {
  if (inherits(cl, "dbscn")) cl <- cl[["cl"]]
  ok <- is.numeric(cl) && length(cl) == n &&
    all(is.na(cl) | (cl >= 1 & cl == round(cl) & cl <= .Machine$integer.max))
  if (!ok) {
    .ckstop(fn, paste("`cl` must be a dbscn result or a cluster per point:",
                      "whole numbers from 1, NA for noise"))
  }
  as.integer(cl)
}

#' Proportional symbol maps
#'
#' `prpmap()` draws a circle at each point (a district's centre, say)
#' whose area is proportional to its value: the radius is `rmax` times the
#' square root of the value over the largest value, so the largest has
#' radius `rmax`. The circles are drawn largest first, so that the smaller
#' stay visible over them, and filled half-transparent. The radii, the
#' drawing order and the box (widened by `rmax` on every side so the
#' circles fit) are computed by the compiled kernel; R only draws.
#'
#' @inheritParams ptmap
#' @param x The points: a numeric matrix or data frame of two columns (x
#'   and y), or the x coordinates with `y`.
#' @param y `NULL` or the y coordinates, one per value of `x`.
#' @param v Numbers from 0, one per point; `NA` draws no circle (not all
#'   missing).
#' @param rmax `NULL` (a fifteenth of the longer side of the points' box)
#'   or the radius of the largest circle, in the units of the coordinates.
#' @param col `NULL` (the palette) or a colour.
#' @return Invisibly, a list: `pts`, a data frame of the points (`x`, `y`,
#'   `v`, `r` the radius, `NA` when not drawn, and `ord` the place in the
#'   drawing order, `NA` when not drawn); `rmax`; and `bbox`, the box drawn
#'   (`xmin`, `xmax`, `ymin`, `ymax`).
#' @examples
#' x <- c(0, 4, 1, 3, 2)
#' y <- c(0, 1, 3, 4, 2)
#' prpmap(x, y, c(120, 30, 5, 60, 0))
#' prpmap(cbind(x, y), v = c(120, 30, NA, 60, 10), rmax = 1)
#' @export
prpmap <- function(x, y = NULL, v, rmax = NULL, ..., main = NULL,
                   col = NULL, pal = NULL) {
  f <- "prpmap"
  p <- .mxy(x, y, f)
  n <- length(p$x)
  ok <- is.numeric(v) && length(v) == n && !all(is.na(v)) &&
    !any(!is.na(v) & !(is.finite(v) & v >= 0))
  if (!ok) {
    .ckstop(f, "`v` must be numbers from 0, one per point, not all missing")
  }
  if (!is.null(rmax) && !.pnum1(rmax, .Machine$double.xmin)) {
    .ckstop(f, "`rmax` must be one positive number")
  }
  v <- as.numeric(v)
  s <- .Call(C_clab_r_mapsym, p$x, p$y, v,
             if (is.null(rmax)) NA_real_ else as.numeric(rmax))
  cl <- .pcols(1L, col, pal, f)
  mn <- .pmain(main, sprintf("Proportional symbols (largest %s)",
                             format(max(v, na.rm = TRUE))), f)
  op <- .ppar(mar = c(0.5, 0.5, 3, 0.5))
  on.exit(graphics::par(op), add = TRUE)
  .mfrm(s$bb, mn, ...)
  o <- s$ord
  graphics::symbols(p$x[o], p$y[o], circles = s$r[o], inches = FALSE,
                    add = TRUE, fg = cl,
                    bg = grDevices::adjustcolor(cl, alpha.f = 0.45))
  dr <- rep(NA_integer_, n)
  dr[o] <- seq_along(o)
  r <- s$r
  r[is.na(v)] <- NA_real_
  invisible(list(
    pts = data.frame(x = p$x, y = p$y, v = v, r = r, ord = dr),
    rmax = s$rmax, bbox = s$bb
  ))
}
