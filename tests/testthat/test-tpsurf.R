test_that("tpsurf() solves the thin-plate system", {
  set.seed(5)
  x <- stats::runif(15)
  y <- stats::runif(15)
  z <- sin(3 * x) + y^2
  phi <- function(r) ifelse(r > 0, r^2 * log(r), 0)
  k <- phi(as.matrix(stats::dist(cbind(x, y))))
  tm <- cbind(1, x, y)
  lam <- 0.01
  a <- rbind(cbind(k + lam * diag(15), tm), cbind(t(tm), matrix(0, 3, 3)))
  s <- solve(a, c(z, 0, 0, 0))
  t1 <- tpsurf(x, y, z, lam = lam, nx = 4, ny = 3)
  expect_s3_class(t1, "tpsurf")
  expect_equal(t1$a, unname(s[1:15]), tolerance = 1e-8)
  expect_equal(unname(t1$beta), unname(s[16:18]), tolerance = 1e-8)
  expect_named(t1$beta, c("const", "x", "y"))
  fit <- as.vector(k %*% s[1:15] + tm %*% s[16:18])
  expect_equal(t1$fitted, fit, tolerance = 1e-8)
  expect_equal(t1$resid, z - t1$fitted, tolerance = 1e-15)
  gx <- seq(min(x), max(x), length.out = 4)
  gy <- seq(min(y), max(y), length.out = 3)
  expect_identical(t1$x, gx)
  expect_identical(t1$y, gy)
  g <- expand.grid(x = gx, y = gy)
  d <- sqrt(outer(g$x, x, "-")^2 + outer(g$y, y, "-")^2)
  sz <- phi(d) %*% s[1:15] + cbind(1, g$x, g$y) %*% s[16:18]
  expect_identical(dim(t1$z), c(4L, 3L))
  expect_equal(as.vector(t1$z), as.vector(sz), tolerance = 1e-8)
  # with no smoothing the surface passes through the values
  t0 <- tpsurf(x, y, z, nx = 5)
  expect_equal(t0$fitted, z, tolerance = 1e-8)
  expect_identical(dim(t0$z), c(5L, 5L))
  # three points: a plane, 1 + x + 2 y, and no bending
  p3 <- tpsurf(c(0, 1, 0), c(0, 0, 1), c(1, 2, 3), nx = 2)
  expect_equal(unname(p3$beta), c(1, 1, 2), tolerance = 1e-12)
  expect_equal(p3$a, c(0, 0, 0), tolerance = 1e-12)
  expect_equal(p3$z[2, 2], 4, tolerance = 1e-12)
  expect_output(print(t1), "tpsurf: thin-plate spline of 15 points")
  pd <- drw(plot(t1, main = "Surface"))
  expect_identical(pd$z, t1$z)
  k <- keeps(plot(t1))
  expect_false(k$visible)
  expect_true(k$par)
})

test_that("tpsurf() refuses what it cannot fit", {
  x <- c(0, 1, 0, 1)
  y <- c(0, 0, 1, 1)
  expect_error(tpsurf(1:3, 1:2, 1:3),
               "tpsurf: `x`, `y` and `z` must be numbers of the same length.",
               fixed = TRUE)
  expect_error(tpsurf(x, y, letters[1:4]),
               "tpsurf: `x`, `y` and `z` must be numbers", fixed = TRUE)
  msg <- "tpsurf: `x`, `y` and `z` must be finite, at least three points."
  expect_error(tpsurf(1:2, 1:2, 1:2), msg, fixed = TRUE)
  expect_error(tpsurf(x, y, c(1, NA, 2, 3)), msg, fixed = TRUE)
  expect_error(tpsurf(x, y, 1:4, lam = -1),
               "tpsurf: `lam` must be one number from 0.", fixed = TRUE)
  expect_error(tpsurf(x, y, 1:4, nx = 1),
               "tpsurf: `nx` and `ny` must be whole numbers from 2 to 1000.",
               fixed = TRUE)
  expect_error(tpsurf(x, y, 1:4, ny = 2.5), "tpsurf: `nx` and `ny` must be",
               fixed = TRUE)
  sing <- paste("tpsurf: the system is singular: the points lie on one line,",
                "or repeat with `lam` 0.")
  expect_error(tpsurf(c(0, 1, 2), c(0, 1, 2), c(1, 2, 3)), sing,
               fixed = TRUE)
  expect_error(tpsurf(c(0, 1, 0, 0), c(0, 0, 1, 0), c(1, 2, 3, 4)), sing,
               fixed = TRUE)
  expect_error(drw(plot(t1 <- tpsurf(x, y, 1:4), main = 2)),
               "plot: `main` must be one string or NULL.", fixed = TRUE)
})
