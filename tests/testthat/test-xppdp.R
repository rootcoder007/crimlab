test_that("xppdp() averages the function over the rows at each grid value", {
  set.seed(43)
  n <- 50
  x <- cbind(a = rnorm(n), b = runif(n))
  f <- function(m) m[, "a"]^2 + 2 * m[, "b"] + m[, "a"] * m[, "b"]
  pd <- function(g) {
    vapply(g, function(v) mean(v^2 + 2 * x[, "b"] + v * x[, "b"]), 0)
  }
  r <- xppdp(x, f, "a", grid = 7)
  g <- seq(min(x[, "a"]), max(x[, "a"]), length.out = 7)
  expect_equal(r$pd$grid, g, tolerance = 1e-15)
  expect_equal(r$pd$pd, pd(g), tolerance = 1e-13)
  expect_identical(r$term, "a")
  expect_identical(r$n, n)
  q <- xppdp(x, f, 1, grid = 6, gtype = "quant")
  g <- unique(stats::quantile(x[, "a"], seq(0, 1, length.out = 6),
                              names = FALSE))
  expect_equal(q$pd$grid, g, tolerance = 1e-15)
  expect_equal(q$pd$pd, pd(g), tolerance = 1e-13)
  u <- xppdp(x, f, "a", grid = c(2, -1, 0))
  expect_identical(u$pd$grid, c(2, -1, 0))
  expect_equal(u$pd$pd, pd(c(2, -1, 0)), tolerance = 1e-13)
  # repeated quantiles are dropped
  z <- cbind(a = c(0, 0, 0, 0, 1), b = 1:5)
  expect_identical(xppdp(z, f, "a", grid = 5, gtype = "quant")$pd$grid,
                   c(0, 1))
})

test_that("xppdp() of a least-squares learner is its slope", {
  set.seed(44)
  x <- cbind(u = rnorm(30), v = rnorm(30))
  y <- 1 - 3 * x[, 1] + x[, 2] + rnorm(30)
  cf <- stats::lm.fit(cbind(1, x), y)$coefficients
  r <- xppdp(x, "glm", "u", y = y, grid = 4, lctl = list(nfold = 2))
  g <- r$pd$grid
  want <- cf[[1L]] + cf[[2L]] * g + cf[[3L]] * mean(x[, 2])
  expect_equal(r$pd$pd, want, tolerance = 1e-9)
})

test_that("xppdp() refuses bad input in words", {
  x <- cbind(a = c(1, 2, 3), b = c(3, 1, 2))
  f <- function(m) m[, 1]
  e <- function(expr, m) expect_error(expr, m, fixed = TRUE)
  e(xppdp(x, f, "z"), "xppdp: `v` must name or number one column of `x`.")
  e(xppdp(x, f, 3), "xppdp: `v` must name or number one column of `x`.")
  e(xppdp(x, f, "a", grid = 1),
    paste("xppdp: `grid` must be a whole number of points from 2, or at",
          "least two values."))
  e(xppdp(x, f, "a", grid = c(1, NA)),
    paste("xppdp: `grid` must be a whole number of points from 2, or at",
          "least two values."))
  e(xppdp(x, f, "a", gtype = "odd"),
    "xppdp: `gtype` must be one of \"even\", \"quant\", not \"odd\".")
})
