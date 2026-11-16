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
