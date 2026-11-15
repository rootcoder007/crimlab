# The drift battery (R/drift.R). The values were computed once by
# rmoriebricklayer 1.1.1 on the same simulated input and printed with
# sprintf("%.17g"); the tolerance is 1e-10.

x <- c(3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5)
y <- c(2, 7, 1, 8, 2, 8, 1, 8, 2, 8, 4, 5)
v <- function(r) as.vector(unclass(r))
t <- 1e-10

test_that("drftks() equals rmoriebricklayer and stats::ks.test()", {
  expect_equal(v(drftks(x, y)),
               c(0.32575757575757569, 0.57633448967605139, 5.7391304347826084),
               tolerance = t)
  expect_equal(v(drftks(x, x)), c(0, 1, 5.5), tolerance = t)
  expect_equal(v(drftks(1, 2)), c(1, 0.69937419913101551, 0.5), tolerance = t)
  expect_equal(v(drftks(1:300, (1:300) + 40)),
               c(0.13333333333333341, 0.0096558989010344324, 150),
               tolerance = t)
  a <- c(0.3, 1.2, -0.5, 2.2, 0.9)
  b <- c(1.1, 2.5, 0.4, 3.1)
  expect_equal(drftks(a, b)[["statistic"]],
               unname(suppressWarnings(stats::ks.test(a, b))$statistic),
               tolerance = testthat_tolerance())
  r <- drftks(c(x, NA), y)
  expect_identical(attr(r, "lakner"), list(nx = 11L, ny = 12L))
  expect_s3_class(r, "drftks")
  expect_output(print(r), "inputs: nx = 11")
  expect_type(r[1:2], "double")
  expect_false(is.object(r[1:2]))
  expect_error(drftks(c(1, Inf), 1), "finite")
  expect_error(drftks(NA, 1), "at least one non-missing value")
  expect_error(drftks("a", 1), "must be numbers")
})

test_that("drftpsi() equals rmoriebricklayer", {
  expect_equal(v(drftpsi(x, y)), c(4.0220529665518514, 0.18532727412781119),
               tolerance = t)
  expect_equal(v(drftpsi(x, y, bins = 4, eps = 0.01)),
               c(0.70485843918983426, 0.084303292908706529), tolerance = t)
  expect_equal(v(drftpsi(c(1, 1, 1), c(9, 9, 9))),
               c(27.630993484907432, 0.69314718055994529), tolerance = t)
  expect_identical(v(drftpsi(x, x)), c(0, 0))
  expect_error(drftpsi(x, y, bins = 1), "at least 2")
  expect_error(drftpsi(x, y, bins = 2e6), "at most 1,000,000")
  expect_error(drftpsi(x, y, eps = 1), "below one")
})

test_that("drftchisq() equals rmoriebricklayer and stats::chisq.test()", {
  r <- drftchisq(c(a = 40, b = 60, c = 100), c(a = 50, b = 30, c = 20))
  expect_equal(v(r), c(126, 2, 4.3596100000630809e-28, -63), tolerance = t)
  expect_equal(r[["statistic"]], unname(stats::chisq.test(
    c(40, 60, 100), p = c(0.5, 0.3, 0.2))$statistic),
               tolerance = testthat_tolerance())
  expect_identical(v(drftchisq(c(a = 5), c(a = 2))), c(0, 0, NA, NA))
  expect_identical(v(drftchisq(c("a", "a", "b", "d"), c("a", "a", "b", "b"))),
                   c(Inf, 1, 0, -Inf))
  expect_equal(v(drftchisq(c("a", "a", "b", "b"), c("a", "b"))), c(0, 1, 1, 0),
               tolerance = testthat_tolerance())
  expect_error(drftchisq(c(1, 2), c(a = 1)), "must be named")
  expect_error(drftchisq(c(a = -1), c(a = 1)), "non-negative counts")
  expect_error(drftchisq(c(a = 0), c(a = 1)), "`obs` has no counts")
  expect_error(drftchisq(c(a = 1), c(a = 0)), "`expd` has no counts")
})

test_that("drfthomog() equals rmoriebricklayer on each branch", {
  set.seed(1)
  a <- sample(c("x", "y", "z"), 300, TRUE)
  b <- sample(c("x", "y", "z"), 300, TRUE, prob = c(0.5, 0.3, 0.2))
  r <- drfthomog(a, b)
  expect_equal(v(r), c(10.89951341174428, 2, 0.004297350083618914,
                       -5.4497567058721401), tolerance = t)
  expect_identical(attr(r, "method"), "chi-square approximation")
  s <- .Random.seed
  r <- drfthomog(c(a = 3, b = 1, c = 2), c(a = 1, b = 4, c = 1))
  expect_identical(.Random.seed, s)
  expect_equal(v(r), c(3.1333333333333333, 2, 0.35736426357364265,
                       -1.0289996712571428), tolerance = t)
  expect_identical(attr(r, "method"),
                   "Monte Carlo permutation (small expected counts)")
  r <- drfthomog(c(a = 4), c(a = 5))
  expect_identical(v(r), c(0, 0, NA, NA))
  expect_identical(attr(r, "method"), "inapplicable: one category")
  r <- drfthomog(c(a = 1, b = 1), c(c = 1, d = 1))
  expect_identical(v(r), c(4, 3, NA, NA))
  expect_match(attr(r, "method"), "4 categories for 4 observations")
  expect_equal(v(drfthomog(c(a = 30, b = 20), c(a = 25, b = 25))),
               c(1.0101010101010099, 1, 0.31487864133641991,
                 -1.1555679799890286), tolerance = t)
  expect_output(print(r), "method: inapplicable")
  expect_error(drfthomog(c(1, 2), c(a = 3)), "must be named")
  expect_error(drfthomog(c(a = 0), c(a = 3)), "at least one observation")
  expect_error(drfthomog(c(a = 3, b = 1), c(a = 1, b = 4), seed = 0.5),
               "one whole number")
})
