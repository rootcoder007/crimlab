# The weighted displacement quotient against its formula in R.

test_that("the quotient and the expected counts follow Bowers and Johnson", {
  a <- c(100, 60)
  b <- c(50, 40)
  cc <- c(200, 180)
  r <- wdqx(a, b, cc)
  sm <- a[2] / cc[2] - a[1] / cc[1]
  bdm <- b[2] / cc[2] - b[1] / cc[1]
  expect_named(r, c("sm", "bdm", "wdq", "ge", "ne", "ea", "eb", "tne"))
  expect_equal(r[["wdq"]], bdm / sm, tolerance = 1e-14)
  expect_gt(r[["wdq"]], 0)
  expect_identical(r[["ge"]], 40)
  expect_equal(r[["ne"]], -sm, tolerance = 1e-15)
  expect_equal(r[["tne"]], (90 - 60) + (45 - 40), tolerance = 1e-13)
  # crime rising in the buffer is displacement: a negative quotient
  expect_lt(wdqx(a, c(50, 70), cc)[["wdq"]], 0)
  expect_true(is.na(wdqx(c(10, 10), b, c(10, 10))[["wdq"]]))
})

test_that("the worked example of Guerette (2009, Appendix D) is reproduced", {
  # Guerette (2009), Analyzing Crime Displacement and Diffusion, Appendix
  # D: GE = 207 - 110 = 97, NE = .326, TNE = 130.5
  r <- wdqx(a = c(207, 110), b = c(178, 157), c = c(308, 318))
  expect_identical(r[["ge"]], 97)
  expect_equal(round(r[["ne"]], 3), 0.326)
  expect_equal(round(r[["tne"]], 1), 130.5)
  # the guide prints WDQ = .411; its equation 3 on these counts gives 0.258
  expect_equal(
    r[["wdq"]],
    (157 / 318 - 178 / 308) / (110 / 318 - 207 / 308),
    tolerance = 1e-14
  )
  expect_equal(round(r[["wdq"]], 3), 0.258)
})

test_that("wdqx() refuses bad input in words", {
  errs(
    "wdqx",
    list(
      list("`a` must have length 2", 1, c(1, 2), c(1, 2)),
      list("`c` must be greater than 0", c(1, 2), c(1, 2), c(0, 2)),
      list("`b` must be at least 0", c(1, 2), c(-1, 2), c(1, 2))
    )
  )
})
