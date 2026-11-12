test_that("gtmp is the product rule at theta 0", {
  r <- gtmp(c(a = 0.1, b = 0.2), c(0.3, NA), hom = c(FALSE, TRUE))
  expect_equal(r[["locus"]][["mp"]], c(2 * 0.1 * 0.3, 0.2^2))
  expect_equal(r[["log10mp"]], log10(0.06 * 0.04))
  expect_equal(r[["log10lr"]], -r[["log10mp"]])
  expect_equal(r[["locus"]][["locus"]], c("a", "b"))
  expect_equal(gtmp(0.5, hom = TRUE)[["locus"]][["mp"]], 0.25)
})

test_that("gtmp matches Balding-Nichols by hand", {
  th <- 0.03
  p <- 0.1
  q <- 0.3
  d <- (1 + th) * (1 + 2 * th)
  hom <- (2 * th + (1 - th) * p) * (3 * th + (1 - th) * p) / d
  het <- 2 * (th + (1 - th) * p) * (th + (1 - th) * q) / d
  r <- gtmp(c(p, p), c(NA, q), hom = c(TRUE, FALSE), theta = th)
  expect_equal(r[["locus"]][["mp"]], c(hom, het), tolerance = 1e-14)
  # shared ancestry makes a match less surprising
  expect_gt(hom, p^2)
  expect_gt(het, 2 * p * q)
  # 400 loci would underflow as a product; the logarithm does not
  big <- gtmp(rep(1e-3, 400), hom = TRUE)
  expect_equal(big[["log10mp"]], -2400)
})

test_that("cllr is 1 for an uninformative system and 0 for a perfect one", {
  expect_equal(cllr(c(1, 1, 1), c(1, 1)), 1)
  expect_lt(cllr(c(1e12, 1e15), c(1e-12, 1e-15)), 1e-10)
  expect_gt(cllr(c(0.01, 0.1), c(100, 10)), 1)
  s <- c(50, 200, 8)
  d <- c(0.01, 0.2, 1.5)
  expect_equal(cllr(s, d),
               0.5 * (mean(log2(1 + 1 / s)) + mean(log2(1 + d))),
               tolerance = 1e-14)
})

test_that("lrpost applies Bayes rule", {
  expect_equal(lrpost(1e6, 1 / 5e6), (1e6 / (5e6 - 1)) / (1 + 1e6 / (5e6 - 1)))
  expect_equal(lrpost(1, 0.3), 0.3)
  expect_equal(lrpost(c(1, 3), c(0.5, 0.25)), c(0.5, 0.5))
  expect_equal(lrpost(1e308, 0.5), 1)
})

test_that("bad forensic input stops with the argument named", {
  expect_error(gtmp(0), "`p`")
  expect_error(gtmp(1.2, hom = TRUE), "at most 1")
  expect_error(gtmp(0.2), "`q`")
  expect_error(gtmp(0.6, 0.6), "`q`")
  expect_error(gtmp(0.2, "x"), "`q`")
  expect_error(gtmp(c(0.2, 0.2), 0.1), "`q`")
  expect_error(gtmp(0.2, hom = NA), "`hom`")
  expect_error(gtmp(0.2, hom = TRUE, theta = 1), "`theta`")
  expect_error(cllr(c(1, 0), 1), "`same`")
  expect_error(cllr(1, -1), "`diff`")
  expect_error(lrpost(0, 0.5), "`lr`")
  expect_error(lrpost(2, 1), "`prior`")
  expect_error(lrpost(2, c(0.1, 0.2)), "`prior`")
})
