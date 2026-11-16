sq <- function(x0, y0) cbind(x0 + c(0, 1, 1, 0), y0 + c(0, 0, 1, 1))

test_that("chomap() classes the areas by quantile, equal or given breaks", {
  p <- list(a = sq(0, 0), b = sq(1, 0), c = sq(0, 1), d = sq(1, 1),
            e = rbind(sq(2, 0), NA, sq(2, 1)))
  v <- c(3, 10, 4, 7, NA)
  r <- drw(chomap(p, v, brk = 3))
  br <- unique(stats::quantile(v, seq(0, 1, length.out = 4), na.rm = TRUE,
                               names = FALSE))
  expect_equal(r$breaks, br, tolerance = 1e-15)
  cls <- findInterval(v, br, rightmost.closed = TRUE, all.inside = TRUE)
  expect_identical(r$areas$cls, cls)
  expect_identical(r$areas$area, names(p))
  expect_identical(r$areas$v, v)
  expect_identical(r$areas$col, blpal("sequential", length(br) - 1L)[cls])
  expect_identical(r$bbox, c(0, 3, 0, 2))
  e <- drw(chomap(p, v, brk = 2, meth = "equal"))
  expect_identical(e$breaks, seq(3, 10, length.out = 3))
  expect_identical(e$areas$cls, c(1L, 2L, 1L, 2L, NA))
  g <- drw(chomap(p, v, brk = c(0, 5, 10)))
  expect_identical(g$areas$cls, findInterval(v, c(0, 5, 10),
                                             rightmost.closed = TRUE,
                                             all.inside = TRUE))
  # repeated quantile breaks are merged: 1, 1, 1, 2 in two classes
  q <- drw(chomap(p[1:4], c(1, 1, 1, 2), brk = 2))
  expect_identical(q$breaks, c(1, 2))
  expect_identical(q$areas$cls, rep(1L, 4))
  one <- drw(chomap(p[1], 2))
  expect_identical(one$breaks, c(1.5, 2.5))
  expect_identical(one$bbox, c(0, 1, 0, 1))
  u <- drw(chomap(unname(p[1:2]), c(1, 2), brk = 2,
                  col = c("red", "blue")))
  expect_identical(u$areas$area, 1:2)
  expect_identical(u$areas$col, c("red", "blue"))
  k <- keeps(chomap(p, v, main = "Rates", pal = "grey"))
  expect_false(k$visible)
  expect_true(k$par)
})

test_that("chomap() refuses what it cannot draw", {
  p <- list(sq(0, 0), sq(1, 0))
  expect_error(chomap(list(), 1),
               "chomap: `polys` must be a non-empty list of areas.",
               fixed = TRUE)
  expect_error(chomap(data.frame(x = 1:3, y = 1:3), 1),
               "chomap: `polys` must be a non-empty list", fixed = TRUE)
  msg <- paste("chomap: `polys[[2]]` must be a two-column numeric matrix of",
               "at least three points (rows of NA between rings).")
  expect_error(chomap(list(sq(0, 0), sq(0, 0)[1:2, ]), 1:2), msg,
               fixed = TRUE)
  expect_error(chomap(list(sq(0, 0), sq(0, 0) * Inf), 1:2), msg,
               fixed = TRUE)
  expect_error(chomap(list(sq(0, 0), rbind(sq(0, 0), c(1, NA))), 1:2), msg,
               fixed = TRUE)
  expect_error(chomap(p, 1), "chomap: `v` must be numbers, one per area.",
               fixed = TRUE)
  expect_error(chomap(p, c("a", "b")), "chomap: `v` must be numbers",
               fixed = TRUE)
  expect_error(chomap(p, c(NA, Inf)), "chomap: `v` has no finite value.",
               fixed = TRUE)
  expect_error(chomap(p, 1:2, meth = "jenks"),
               "chomap: `meth` must be \"quantile\" or \"equal\".",
               fixed = TRUE)
  msg <- paste("chomap: `brk` must be a whole number of classes from 1 to 9,",
               "or increasing breaks covering `v`.")
  expect_error(chomap(p, 1:2, brk = 10), msg, fixed = TRUE)
  expect_error(chomap(p, 1:2, brk = 2.5), msg, fixed = TRUE)
  expect_error(chomap(p, 1:2, brk = c(1.5, 3)), msg, fixed = TRUE)
  expect_error(chomap(p, 1:2, brk = c(2, 1)), msg, fixed = TRUE)
  expect_error(chomap(p, 1:2, main = 1),
               "chomap: `main` must be one string or NULL.", fixed = TRUE)
  expect_error(chomap(p, 1:2, pal = "neon"),
               "chomap: the palette must be one of", fixed = TRUE)
  expect_error(chomap(p, 1:2, col = "nocolour"),
               "chomap: `col` must be colours.", fixed = TRUE)
})

test_that("smmap() shares one set of breaks over the periods", {
  p <- list(sq(0, 0), sq(1, 0), sq(2, 0))
  v <- cbind(`2023` = c(1, 5, 9), `2024` = c(2, NA, 12))
  r <- drw(smmap(p, v, brk = 2, meth = "equal"))
  br <- seq(1, 12, length.out = 3)
  expect_identical(r$breaks, br)
  expect_identical(r$areas$prd, rep(c("2023", "2024"), each = 3))
  expect_identical(r$areas$area, rep(1:3, 2))
  expect_identical(r$areas$v, as.vector(v))
  cls <- findInterval(as.vector(v), br, rightmost.closed = TRUE,
                      all.inside = TRUE)
  expect_identical(r$areas$cls, cls)
  expect_identical(r$areas$col, blpal("sequential", 2)[cls])
  expect_identical(r$bbox, c(0, 3, 0, 1))
  q <- drw(smmap(p, v, brk = 4))
  qb <- unique(stats::quantile(v, seq(0, 1, length.out = 5), na.rm = TRUE,
                               names = FALSE))
  expect_equal(q$breaks, qb, tolerance = 1e-15)
  d <- drw(smmap(p, as.data.frame(v)))
  expect_identical(d$areas$prd, rep(c("2023", "2024"), each = 3))
  u <- drw(smmap(p, unname(v)))
  expect_identical(u$areas$prd, rep(c("1", "2"), each = 3))
  k <- keeps(smmap(p, v, main = "By year"))
  expect_false(k$visible)
  expect_true(k$par)
})

test_that("smmap() refuses what it cannot draw", {
  p <- list(sq(0, 0), sq(1, 0))
  msg <- paste("smmap: `v` must be a numeric matrix, one row per area and",
               "one column per period.")
  expect_error(smmap(p, 1:2), msg, fixed = TRUE)
  expect_error(smmap(p, matrix(1, 3, 2)), msg, fixed = TRUE)
  expect_error(smmap(p, matrix("a", 2, 2)), msg, fixed = TRUE)
  expect_error(smmap(p, data.frame(a = c("x", "y"))), msg, fixed = TRUE)
  expect_error(smmap(list(), matrix(1)),
               "smmap: `polys` must be a non-empty list of areas.",
               fixed = TRUE)
  expect_error(smmap(p, matrix(NA_real_, 2, 2)),
               "smmap: `v` has no finite value.", fixed = TRUE)
  expect_error(smmap(p, matrix(1:4, 2), brk = 0),
               "smmap: `brk` must be a whole number", fixed = TRUE)
  expect_error(smmap(p, matrix(1:4, 2), main = c("a", "b")),
               "smmap: `main` must be one string or NULL.", fixed = TRUE)
})

test_that("ptmap() draws points by cluster, with the clusters' centres", {
  x <- c(0, 0.2, 0.1, 5, 5.2, 9)
  y <- c(0, 0.1, 0.2, 5, 5.1, 0)
  cl <- c(1L, 1L, 1L, 2L, 2L, NA)
  r <- drw(ptmap(x, y, cl = cl))
  expect_identical(r$pts$x, x)
  expect_identical(r$pts$cl, cl)
  expect_identical(r$cent$cl, 1:2)
  expect_identical(r$cent$n, c(3L, 2L))
  expect_equal(r$cent$x, as.vector(tapply(x, cl, mean)), tolerance = 1e-14)
  expect_equal(r$cent$y, as.vector(tapply(y, cl, mean)), tolerance = 1e-14)
  expect_identical(r$bbox, c(0, 9, 0, 5.1))
  pl <- blpal(plotopt()$palette, 2)
  expect_identical(r$pts$col, c(rep(pl[1], 3), rep(pl[2], 2), "#969696"))
  expect_identical(r$pts$pch, c(16L, 16L, 16L, 17L, 17L, 4L))
  # a number no point has: its centre is NaN
  s <- drw(ptmap(x, y, cl = c(1, 1, 1, 3, 3, NA), col = c("red", "blue")))
  expect_identical(s$cent$n, c(3L, 0L, 2L))
  expect_true(is.nan(s$cent$x[2]))
  expect_identical(s$pts$col[4], "red")
  m <- drw(ptmap(cbind(x, y)))
  expect_true(all(is.na(m$pts$cl)))
  expect_identical(nrow(m$cent), 0L)
  expect_identical(m$pts$col, rep(blpal(plotopt()$palette, 1), 6))
  d <- drw(ptmap(data.frame(a = x, b = y), cl = cl))
  expect_identical(d$pts$y, y)
  one <- drw(ptmap(2, 3))
  expect_identical(one$bbox, c(1.5, 2.5, 2.5, 3.5))
  k <- keeps(ptmap(x, y, cl = cl, main = "Clusters"))
  expect_false(k$visible)
  expect_true(k$par)
})

test_that("plot() of a dbscn result draws its clusters", {
  set.seed(1)
  p <- rbind(cbind(stats::rnorm(20, 0, 0.1), stats::rnorm(20, 0, 0.1)),
             cbind(stats::rnorm(20, 3, 0.1), stats::rnorm(20, 3, 0.1)),
             c(9, 9))
  d <- dbscn(p, eps = 0.5, minpts = 4)
  g <- drw(plot(d, p))
  expect_identical(g$pts$cl, d$cl)
  expect_identical(g$cent$n, d$sizes)
  for (j in seq_along(d$sizes)) {
    expect_equal(g$cent$x[j], mean(p[which(d$cl == j), 1]),
                 tolerance = 1e-14)
  }
  k <- keeps(plot(d, p))
  expect_false(k$visible)
  expect_true(k$par)
  expect_error(drw(plot(d)),
               "plot: `y` must be the points the clusters were found in.",
               fixed = TRUE)
  expect_error(drw(plot(d, p, main = 1)),
               "plot: `main` must be one string or NULL.", fixed = TRUE)
  expect_error(drw(plot(d, p[1:5, ])), "ptmap: `cl` must be a dbscn result",
               fixed = TRUE)
})

test_that("ptmap() refuses points it cannot draw", {
  expect_error(ptmap("a"),
               paste("ptmap: `x` must be a two-column matrix or data frame",
                     "of coordinates, or numbers with `y`."), fixed = TRUE)
  expect_error(ptmap(1:3), "ptmap: `x` must be a two-column matrix",
               fixed = TRUE)
  expect_error(ptmap(1:3, 1:2),
               "ptmap: `x` and `y` must be numbers of the same length.",
               fixed = TRUE)
  expect_error(ptmap(c(1, NA), c(1, 2)),
               "ptmap: the coordinates must be finite, at least one point.",
               fixed = TRUE)
  expect_error(ptmap(numeric(0), numeric(0)),
               "ptmap: the coordinates must be finite", fixed = TRUE)
  msg <- paste("ptmap: `cl` must be a dbscn result or a cluster per point:",
               "whole numbers from 1, NA for noise.")
  expect_error(ptmap(1:2, 1:2, cl = 1), msg, fixed = TRUE)
  expect_error(ptmap(1:2, 1:2, cl = c(0, 1)), msg, fixed = TRUE)
  expect_error(ptmap(1:2, 1:2, cl = c(1.5, 1)), msg, fixed = TRUE)
  expect_error(ptmap(1:2, 1:2, cl = c("a", "b")), msg, fixed = TRUE)
  expect_error(ptmap(1:2, 1:2, col = "nocolour"),
               "ptmap: `col` must be colours.", fixed = TRUE)
})

test_that("prpmap() scales circle areas to the values", {
  x <- c(0, 3, 0, 1)
  y <- c(0, 0, 4, 1)
  v <- c(4, 1, NA, 4)
  r <- drw(prpmap(x, y, v))
  # the box is 3 by 4: the largest radius is 4 / 15
  rm <- 4 / 15
  expect_equal(r$rmax, rm, tolerance = 1e-15)
  expect_equal(r$pts$r, rm * sqrt(v / 4), tolerance = 1e-15)
  expect_identical(r$pts$ord, c(1L, 3L, NA, 2L))
  expect_identical(r$pts$ord[order(-v, na.last = NA)], 1:3)
  expect_equal(r$bbox, c(-rm, 3 + rm, -rm, 4 + rm), tolerance = 1e-15)
  g <- drw(prpmap(cbind(x, y), v = v, rmax = 2))
  expect_identical(g$pts$r, 2 * sqrt(v / 4))
  expect_identical(g$bbox, c(-2, 5, -2, 6))
  z <- drw(prpmap(x, y, c(0, 0, 0, 0)))
  expect_identical(z$pts$r, rep(0, 4))
  one <- drw(prpmap(1, 1, 5))
  expect_equal(one$rmax, 1 / 15, tolerance = 1e-15)
  k <- keeps(prpmap(x, y, v, main = "Counts", col = "darkred"))
  expect_false(k$visible)
  expect_true(k$par)
})

test_that("prpmap() refuses what it cannot draw", {
  x <- c(0, 3, 0, 1)
  msg <- "prpmap: `v` must be numbers from 0, one per point, not all missing."
  expect_error(prpmap(x, x, c(1, -1, 2, 3)), msg, fixed = TRUE)
  expect_error(prpmap(x, x, rep(NA_real_, 4)), msg, fixed = TRUE)
  expect_error(prpmap(x, x, c(1, 2)), msg, fixed = TRUE)
  expect_error(prpmap(x, x, c(1, Inf, 2, 3)), msg, fixed = TRUE)
  expect_error(prpmap(x, x, letters[1:4]), msg, fixed = TRUE)
  expect_error(prpmap(x, x, 1:4, rmax = 0),
               "prpmap: `rmax` must be one positive number.", fixed = TRUE)
  expect_error(prpmap(x, x, 1:4, rmax = c(1, 2)),
               "prpmap: `rmax` must be one positive number.", fixed = TRUE)
  expect_error(prpmap(x, 1:3, 1:4),
               "prpmap: `x` and `y` must be numbers of the same length.",
               fixed = TRUE)
  expect_error(prpmap(x, x, 1:4, col = "nocolour"),
               "prpmap: `col` must be colours.", fixed = TRUE)
})
