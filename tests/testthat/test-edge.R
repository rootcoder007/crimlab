# Failure modes recorded for other packages computing the same quantities
# (in the author's notes; the file is named in each comment). Each
# test shows crimlab does not share the failure. The kernels' side of each
# lives in src/t_*.cpp under ids edgeNN; here what R forms itself.

test_that("evidence weights are summed as logarithms and never underflow", {
  # comparison.md: a likelihood ratio formed as exp(a) / exp(b) is 0 / 0
  # for an atypical item. Here 500 loci each worth 1e-3 combine to 1e-1500,
  # far below the smallest double, and the logarithm is still exact.
  p <- rep(1e-3, 500)
  r <- kinidx(cbind(seq_along(p), seq_along(p)), cbind(seq_along(p),
                                                       seq_along(p)),
              cbind(p, p), rel = c(0.9999, 0.0001, 0))
  expect_equal(r[["log10idx"]], 500 * log10(0.9999 + 0.0001 / 1e-3),
               tolerance = 1e-10)
  # the posterior of a ratio of 1e300 at a prior of 1e-300 is one half
  expect_equal(lrpost(1e300, 1e-300), 0.5)
  expect_equal(lrpost(1e-300, 1 - 1e-16), 1 / (1 + 1e-16 / (1 - 1e-16) *
                                                  1e300), tolerance = 1e-10)
})

test_that("allele frequencies outside their range are refused, not used", {
  # DNAtools.md: frequencies that do not form a distribution were used
  # silently. crimlab stops on any frequency that cannot be one.
  expect_error(mixlr(c(a = 0.5, b = 0.5 + 1e-9), "a", 1, "a", 1), "at most 1")
  expect_error(gtmp(0.5, 0.5 + 1e-9), "`q`")
  expect_error(kinidx(c(1, 2), c(1, 2), c(0, 0.2)), "`freq`")
  expect_error(patidx(c(1, 2), c(1, 2), c(0.1, 1.01), mother = c(1, 1)),
               "`freq`")
})

test_that("single-source match statistics are the closed forms", {
  # forensim.md lists these as the checks it could not run there
  p <- 0.07
  q <- 0.11
  expect_equal(gtmp(p, q)[["log10lr"]], -log10(2 * p * q), tolerance = 1e-10)
  expect_equal(gtmp(p, hom = TRUE)[["log10lr"]], -log10(p^2), tolerance = 1e-10)
})

test_that("a group never seen in daylight is no warning", {
  expect_no_warning(r <- vodtst(c(a = 50, b = 0), c(a = 40, b = 3)))
})

test_that("the cross-entropy plot is finite for a subnormal ratio", {
  want <- 0.5 * (310 * log2(10) + 1)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_equal(eceplt(1e-310, 1, lo = 0)[["ece"]], want, tolerance = 1e-12)
})

test_that("orrisk refuses an odds ratio whose reciprocal overflows", {
  expect_error(orrisk(1e-320, 0.5, 0.5), "reciprocal")
  expect_error(.Call(C_clab_orrisk, 1e-320, 0.5, 0.5),
               "outside the domain")
})
