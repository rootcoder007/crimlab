# Bounds under monotone response and selection: the kernel is checked by
# hand in src/t_mtr.cpp; here the binding against the propositions of
# Manski and Pepper written out in R.

test_that("the bounds follow the formulas for every assumption", {
  set.seed(13)
  z <- factor(
    sample(c("fine", "probation", "custody"), 300, TRUE),
    levels = c("fine", "probation", "custody")
  )
  y <- rbinom(300, 1, c(0.2, 0.35, 0.5)[as.integer(z)])
  P <- as.numeric(table(z)) / 300
  m <- as.numeric(tapply(y, z, mean))
  k <- 3
  bd <- function(t, as) {
    lt <- seq_len(k) < t
    gt <- seq_len(k) > t
    switch(
      as,
      none = c(m[t] * P[t], m[t] * P[t] + 1 - P[t]),
      mtr = c(sum((m * P)[!gt]), sum(P[lt]) + sum((m * P)[!lt])),
      mts = c(m[t] * sum(P[!lt]), m[t] * sum(P[!gt]) + sum(P[gt])),
      both = c(
        sum((m * P)[lt]) + m[t] * sum(P[!lt]),
        m[t] * sum(P[!gt]) + sum((m * P)[gt])
      )
    )
  }
  for (as in c("none", "mtr", "mts", "both")) {
    r <- mtrbnd(y, z, assume = as)
    b <- sapply(1:3, bd, as = as)
    expect_equal(r$bnd$lower, b[1, ], tolerance = 1e-14)
    expect_equal(r$bnd$upper, b[2, ], tolerance = 1e-14)
    lo <- b[1, ] - b[2, 1]
    if (as %in% c("mtr", "both")) {
      lo <- pmax(lo, 0)
    }
    expect_equal(r$bnd$atelo, c(NA, lo[-1]), tolerance = 1e-14)
    expect_equal(r$bnd$atehi, c(NA, (b[2, ] - b[1, 1])[-1]), tolerance = 1e-14)
  }
  expect_identical(r$bnd$z, c("fine", "probation", "custody"))
  # the assumptions narrow the worst-case interval
  w0 <- mtrbnd(y, z, assume = "none")$bnd
  w3 <- mtrbnd(y, z)$bnd
  expect_true(all(w3$upper - w3$lower <= w0$upper - w0$lower + 1e-15))
  expect_output(print(w3 <- mtrbnd(y, z)), "bounds under both")
})

test_that("weights and other ranges", {
  r <- mtrbnd(
    c(2, 5, 3, 9),
    c(1, 1, 2, 2),
    ylo = 0,
    yhi = 10,
    w = c(1, 3, 1, 1),
    assume = "none"
  )
  expect_equal(r$bnd$mean, c((2 + 15) / 4, 6))
  expect_equal(r$bnd$share, c(4, 2) / 6)
  expect_equal(r$bnd$lower[1], (17 / 4) * (4 / 6))
})

test_that("mtrbnd() refuses bad input in words", {
  errs(
    "mtrbnd",
    list(
      list("`y` must lie between", c(0, 2), c(1, 2)),
      list("`z` must hold two treatments", c(0, 1), c(1, 1)),
      list("`ylo` must be below", c(0, 1), c(1, 2), ylo = 1, yhi = 1),
      list("`assume`", c(0, 1), c(1, 2), assume = "mono")
    )
  )
})
