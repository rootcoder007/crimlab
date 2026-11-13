# Empirical-Bayes shrinkage of hot-spot counts. The kernel's optimality is
# checked in src/t_hot.cpp; here the
# binding against the formulas written out.

test_that("shrinkage follows its formulas", {
  y <- c(40, 12, 9, 25, 7, 31, 5, 18)
  s <- ebshr(y)
  expect_s3_class(s, "ebshr")
  expect_named(s, c("fit", "est"))
  m <- mean(y)
  tot <- mean((y - m)^2)
  b <- m / (m + tot - m)
  expect_equal(unname(s$fit), c(m, tot, m, tot - m, b, b), tolerance = 1e-15)
  expect_named(s$fit, c("mean", "total", "noise", "signal", "bstar", "b"))
  expect_equal(s$est$shrunk, (1 - b) * y + b * m, tolerance = 1e-15)
  expect_equal(s$est$fall, b * (y - m), tolerance = 1e-14)
  w <- c(2, 1, 1, 3, 1, 1, 2, 1)
  sw <- ebshr(y, v = 10, w = w, b = 0.25)
  mw <- sum(w * y) / sum(w)
  tw <- sum(w * (y - mw)^2) / sum(w)
  expect_equal(unname(sw$fit), c(mw, tw, 10, tw - 10, 10 / tw, 0.25),
               tolerance = 1e-14)
  # all noise: everything shrinks to the mean
  expect_identical(ebshr(c(3, 4, 5), v = 100)$fit[["bstar"]], 1)
  expect_identical(ebshr(c(3, 3, 3), v = 0)$fit[["bstar"]], 0)
  expect_output(print(s), "ebshr: shrinkage")
})

test_that("ebshr() refuses bad input in words", {
  errs("ebshr", list(
    list("at least 2 places", 1),
    list("negative: give the noise variance", c(-1, -2)),
    list("`w` must have a positive sum", c(1, 2), w = c(0, 0)),
    list("`w` must have a positive sum", c(1, 2), v = 1, w = c(0, 0)),
    list("`v`", c(1, 2), v = -1),
    list("`b`", c(1, 2), b = 2)
  ))
})
