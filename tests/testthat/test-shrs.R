# Shares with Wilson intervals (R/shrs.R). The values were computed once
# by rmoriebricklayer 1.1.1's shrs() on the same simulated tables.

s <- data.frame(division = c("North", "South", "East", "North", "East"),
                yr = c(1, 1, 1, 2, 2), stops = c(412, 77, 3, 20, NA))

test_that("the shares and intervals equal rmoriebricklayer's", {
  a <- shrs(s, stops, by = c("division", "yr"))
  expect_identical(a$division, c("East", "North", "South", "East", "North"))
  expect_identical(a$count, c(3, 412, 77, 0, 20))
  expect_identical(a$share, c(0.5859375, 80.46875, 15.0390625, 0, 3.90625))
  expect_equal(a$lower, c(0.19946757921860914, 76.813211410143509,
                          12.203768918306174, 0, 2.5427079865239395),
               tolerance = 1e-14)
  expect_equal(a$upper, c(1.708378068394411, 83.670488506743311,
                          18.39506258731727, 0.7446975722882736,
                          5.956310087929312), tolerance = 1e-14)
  b <- shrs(c(a = 0, b = 7, c = 3), conflvl = 0.9)
  expect_identical(b$label, c("a", "b", "c"))
  expect_equal(b$lower, c(0, 44.169979586983523, 12.687658390319791),
               tolerance = 1e-14)
  expect_equal(b$upper, c(21.294197008340682, 87.312341609680203,
                          55.830020413016477), tolerance = 1e-14)
  d <- shrs(s, "stops", total = 10000)
  expect_equal(c(d$lower, d$upper), c(4.7049863664017622, 5.569481327350168),
               tolerance = 1e-14)
  expect_s3_class(d, "shrs")
  expect_identical(attr(d, "share")$total, 10000)
  expect_output(print(d), "Wilson")
})

test_that("the ends are exact and the shares of a full grouping sum to 100", {
  u <- unclass(shrs(c(a = 0, b = 10)))
  expect_identical(u$lower[1L], 0)
  expect_identical(u$upper[2L], 100)
  expect_equal(sum(shrs(c(3, 7))$share), 100, tolerance = testthat_tolerance())
  expect_false("label" %in% names(shrs(c(3, 7))))
  expect_identical(shrs(s, "stops", total = "stops")$total, 512)
})

test_that("the refusals are those of rmoriebricklayer", {
  d <- data.frame(n = c(1, 2), g = c("a", "b"))
  expect_error(shrs(d, count = "n", by = "zz"), "`by` column not found: zz")
  expect_error(shrs(data.frame(n = c("a", "b")), count = "n"),
               "`n` must be numeric")
  expect_error(shrs(d, count = "n", total = 0), "the total must be positive")
  expect_error(shrs(d, count = "nn"), "`count` column not found: nn")
  expect_error(shrs(d, count = "n", total = "zz"), "`total` column not found")
  expect_error(shrs(d, count = "n", total = -1), "a column name or one")
  expect_error(shrs(d), "`count` must name the column")
  expect_error(shrs(list()), "`x` must be finite counts")
  expect_error(shrs(c(NA, NA)), "`x` must be finite counts")
  expect_error(shrs(c(3, 7), conflvl = 2), "strictly between 0 and 1")
  expect_error(shrs(c(5, 2), total = 4), "larger than the total")
  expect_error(shrs(c(-1, 2)), "non-negative")
  expect_error(shrs(c(1, 2), by = "g"), "unused argument")
})
