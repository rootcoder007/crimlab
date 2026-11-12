# Failure modes recorded for other packages computing the same quantities
# (dossiers in the author's notes; the file is named in each comment). Each
# test shows csls does not share the failure.

test_that("evidence weights are summed as logarithms and never underflow", {
  # comparison.md: a likelihood ratio formed as exp(a) / exp(b) is 0 / 0
  # for an atypical item. Here 500 loci each worth 1e-3 combine to 1e-1500,
  # far below the smallest double, and the logarithm is still exact.
  p <- rep(1e-3, 500)
  r <- kinidx(cbind(seq_along(p), seq_along(p)), cbind(seq_along(p),
                                                       seq_along(p)),
              cbind(p, p), rel = c(0.9999, 0.0001, 0))
  expect_true(all(is.finite(r[["locus"]][["idx"]])))
  expect_equal(r[["log10idx"]], 500 * log10(0.9999 + 0.0001 / 1e-3),
               tolerance = 1e-10)
  expect_equal(gtmp(rep(1e-2, 500), hom = TRUE)[["log10mp"]], -2000,
               tolerance = 1e-10)
  # the posterior of a ratio of 1e300 at a prior of 1e-300 is one half
  expect_equal(lrpost(1e300, 1e-300), 0.5)
  expect_equal(lrpost(1e-300, 1 - 1e-16), 1 / (1 + 1e-16 / (1 - 1e-16) *
                                                  1e300), tolerance = 1e-10)
  # the cost of extreme ratios is finite: log2(1e300) bits, no Inf, no NaN
  expect_equal(cllr(1e-300, 1e-300), 0.5 * log2(1 + 1e300), tolerance = 1e-14)
  expect_equal(cllr(1, 1e300), 0.5 * (1 + log2(1 + 1e300)), tolerance = 1e-14)
})

test_that("allele frequencies outside their range are refused, not used", {
  # DNAtools.md: frequencies that do not form a distribution were used
  # silently. csls stops on any frequency that cannot be one.
  expect_error(mixlr(c(a = 0.5, b = 0.5 + 1e-9), "a", 1, "a", 1), "at most 1")
  expect_error(gtmp(0.5, 0.5 + 1e-9), "`q`")
  expect_error(kinidx(c(1, 2), c(1, 2), c(0, 0.2)), "`freq`")
  expect_error(patidx(c(1, 2), c(1, 2), c(0.1, 1.01), mother = c(1, 1)),
               "`freq`")
  # a heterozygote of the only two alleles at a locus is allowed
  expect_equal(gtmp(0.5, 0.5)[["locus"]][["mp"]], 0.5)
})

test_that("single-source match statistics are the closed forms", {
  # forensim.md lists these as the checks it could not run there
  p <- 0.07
  q <- 0.11
  expect_equal(gtmp(p, q)[["log10lr"]], -log10(2 * p * q), tolerance = 1e-10)
  expect_equal(gtmp(p, hom = TRUE)[["log10lr"]], -log10(p^2), tolerance = 1e-10)
})

test_that("a group never seen in daylight gives a ratio of 0, not a fit", {
  # searchlight.md: a near-separated logistic fit was labelled estimable.
  # The veil-of-darkness ratio here is a closed form in the counts: a zero
  # gives exactly 0 and a lopsided table a finite ratio, with no warning.
  expect_no_warning(r <- vodtst(c(a = 50, b = 0), c(a = 40, b = 3)))
  expect_equal(r[["ratio"]], c(1, 0), tolerance = 1e-10)
  expect_equal(vodtst(c(1, 1e6), c(1e6, 1))[["ratio"]], c(1, 1e12),
               tolerance = 1e-10)
})

# Two defects these tests found in csls itself, now fixed in the kernels.

test_that("cllr is finite for a subnormal same-source ratio", {
  # log1p(1 / same) overflowed once 1 / same passed the largest double; the
  # cost itself, 0.5 * (log2(1 + 1e310) + 1), is about 515.4 bits
  want <- 0.5 * (310 * log2(10) + 1)
  expect_equal(cllr(1e-310, 1), want, tolerance = 1e-12)
  expect_equal(cllrmin(1e-310, 1)[["cllr"]], want, tolerance = 1e-12)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_equal(eceplt(1e-310, 1, lo = 0)[["ece"]], want, tolerance = 1e-12)
  # both branches of the rewritten term agree where they meet
  expect_equal(cllr(1 - 1e-12, 1), cllr(1 + 1e-12, 1), tolerance = 1e-12)
})

test_that("orrisk stays in [0, 1] for a tiny odds ratio", {
  # OR - 1 rounded to -1 and 1 + (OR - 1) b cancelled as b neared 1: the
  # exposed risk came back above 1 (1e-15) or Inf (1e-17)
  for (or in c(1e-8, 1e-15, 1e-17, 1e-300)) {
    for (q in c(0.01, 0.5, 0.99)) {
      for (s in c(0.02, 0.5, 0.98)) {
        r <- orrisk(or, q, s)
        a <- r[["exposed"]]
        b <- r[["unexposed"]]
        expect_gt(a, 0)
        expect_lt(a, 1)
        expect_gt(b, 0)
        # the unexposed risk can be within rounding of 1
        expect_lte(b, 1)
        expect_equal(s * a + (1 - s) * b, q, tolerance = 1e-12)
      }
    }
  }
  # the mirror image: a tiny ratio is the reciprocal with the groups exchanged
  lo <- orrisk(1e-6, 0.3, 0.2)
  hi <- orrisk(1e6, 0.3, 0.8)
  expect_equal(lo[["exposed"]], hi[["unexposed"]], tolerance = 1e-14)
  expect_equal(lo[["unexposed"]], hi[["exposed"]], tolerance = 1e-14)
  expect_error(orrisk(1e-320, 0.5, 0.5), "reciprocal")
  expect_error(.Call(C_csls_orrisk, 1e-320, 0.5, 0.5),
               "outside the domain")
})
