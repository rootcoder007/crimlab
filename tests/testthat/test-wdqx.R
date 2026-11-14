# The weighted displacement quotient against its formula in R.

test_that("the quotient and the expected counts follow Bowers and Johnson", {
  a <- c(100, 60)
  b <- c(50, 40)
  cc <- c(200, 180)
  r <- wdqx(a, b, cc)
  sm <- a[2] / cc[2] - a[1] / cc[1]
  bdm <- b[2] / cc[2] - b[1] / cc[1]
  expect_named(r, c("sm", "bdm", "wdq", "ea", "eb", "net"))
  expect_equal(r[["wdq"]], bdm / sm, tolerance = 1e-14)
  expect_gt(r[["wdq"]], 0)
  expect_equal(r[["net"]], (90 - 60) + (45 - 40), tolerance = 1e-13)
  # crime rising in the buffer is displacement: a negative quotient
  expect_lt(wdqx(a, c(50, 70), cc)[["wdq"]], 0)
  expect_true(is.na(wdqx(c(10, 10), b, c(10, 10))[["wdq"]]))
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
