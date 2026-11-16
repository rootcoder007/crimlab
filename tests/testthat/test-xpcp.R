test_that("xpcp() predicts each row with one input moved over the grid", {
  set.seed(45)
  x <- cbind(a = rnorm(40), b = runif(40))
  f <- function(m) exp(m[, "a"]) * m[, "b"]
  r <- xpcp(x, f, "a", newx = x[2:4, ], grid = 5)
  g <- seq(min(x[, "a"]), max(x[, "a"]), length.out = 5)
  expect_equal(r$grid, g, tolerance = 1e-15)
  expect_equal(r$prof, outer(x[2:4, "b"], exp(g)), tolerance = 1e-14)
  expect_equal(r$pred, f(x[2:4, ]), tolerance = 1e-15)
  expect_identical(r$term, "a")
  # the default row is the first; the mean profile over x is xppdp()
  expect_equal(xpcp(x, f, 2, grid = c(0, 1))$prof,
               matrix(exp(x[1L, "a"]) * c(0, 1), 1L), tolerance = 1e-15)
  a <- xpcp(x, f, "a", newx = x, grid = 4)
  expect_equal(colMeans(a$prof), xppdp(x, f, "a", grid = 4)$pd$pd,
               tolerance = 1e-14)
})

test_that("xpcp() refuses bad input in words", {
  x <- cbind(a = c(1, 2, 3), b = c(3, 1, 2))
  f <- function(m) m[, 1]
  m <- paste("xpcp: `newx` must be a numeric matrix with the 2 columns of",
             "`x`, no missing values.")
  expect_error(xpcp(x, f, "a", newx = x[, 1, drop = FALSE]), m, fixed = TRUE)
  expect_error(xpcp(x, f, "a", newx = c(1, 2)), m, fixed = TRUE)
  expect_error(xpcp(x, f, "a", newx = cbind(1, NA)), m, fixed = TRUE)
  expect_error(xpcp(x, f, "c"),
               "xpcp: `v` must name or number one column of `x`.",
               fixed = TRUE)
})
