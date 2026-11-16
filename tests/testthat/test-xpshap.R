test_that("xpshap() of a linear function against one reference row", {
  x <- cbind(a = 1, b = 2, c = -1)
  b <- c(2, -1, 0.5)
  f <- function(m) drop(m %*% b)
  newx <- rbind(c(3, 0, 1), c(-1, 4, 2))
  r <- xpshap(x, f, newx = newx, nsamp = 7)
  want <- t(t(sweep(newx, 2L, x[1L, ])) * b)
  expect_equal(unname(r$phi), want, tolerance = 1e-14)
  expect_equal(unname(r$se), matrix(0, 2, 3), tolerance = 1e-14)
  expect_identical(colnames(r$phi), c("a", "b", "c"))
  expect_equal(r$pred, f(newx), tolerance = 1e-15)
  expect_equal(r$base, f(x), tolerance = 1e-15)
})

test_that("xpshap() samples the Shapley values of a game with interactions", {
  # exact values by enumerating the coalitions: v(S) = f(newx on S, w off)
  x <- cbind(a = 0.5, b = -1, c = 2)
  newx <- cbind(a = 2, b = 1, c = -1)
  f <- function(m) m[, 1] * m[, 2] + exp(m[, 3] / 2) * m[, 1] + m[, 2]^2
  v <- function(s) {
    z <- x
    z[1L, s] <- newx[1L, s]
    f(z)
  }
  exact <- vapply(1:3, function(j) {
    o <- setdiff(1:3, j)
    ss <- list(integer(0), o[1L], o[2L], o)
    sum(vapply(ss, function(s) {
      k <- length(s)
      factorial(k) * factorial(2 - k) / 6 * (v(c(s, j)) - v(s))
    }, 0))
  }, 0)
  r <- xpshap(x, f, newx = newx, nsamp = 4000, seed = 2)
  expect_true(all(abs(r$phi[1L, ] - exact) <= 5 * r$se[1L, ] + 1e-12))
  expect_true(all(r$se > 0))
  # the exact values add up to the prediction less the reference
  expect_equal(sum(exact), unname(f(newx) - f(x)), tolerance = 1e-13)
  expect_identical(xpshap(x, f, newx = newx, nsamp = 5, seed = 2)$phi,
                   xpshap(x, f, newx = newx, nsamp = 5, seed = 2)$phi)
})

test_that("xpshap() refuses bad input in words", {
  x <- cbind(a = c(1, 2), b = c(3, 1))
  f <- function(m) m[, 1]
  expect_error(xpshap(x, f, nsamp = 0), "xpshap: `nsamp` must be at least 1.",
               fixed = TRUE)
  expect_error(xpshap(x, f, newx = 1),
               paste("xpshap: `newx` must be a numeric matrix with the 2",
                     "columns of `x`, no missing values."), fixed = TRUE)
  expect_error(xpshap(x, f, seed = NA),
               "xpshap: `seed` must be finite, with no missing values.",
               fixed = TRUE)
})
