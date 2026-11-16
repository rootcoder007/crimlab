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
