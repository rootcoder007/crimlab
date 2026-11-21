# svypci() and svyqual() (R/svyci.R). The interval references were
# computed once by survey 4.4 (svyciprop() on the api data, a one-stage
# cluster sample of 183 schools in 15 districts) and printed with
# sprintf("%.15g"); the tolerance is 1e-9.

p <- 0.87431693989071035
s <- 0.02035947724462743

test_that("svypci() equals svyciprop() of survey, beta and logit", {
  b <- svypci(p, s, 183, df = 14, trunc = FALSE)
  expect_named(b, c("est", "lo", "hi", "neff"))
  expect_equal(
    c(b[["lo"]], b[["hi"]]),
    c(0.823702210847501, 0.914704492126865),
    tolerance = 1e-9
  )
  l <- svypci(p, s, 183, df = 14, method = "logit")
  expect_equal(
    c(l[["lo"]], l[["hi"]]),
    c(0.823798892379033, 0.911900077021222),
    tolerance = 1e-9
  )
  # the effective size is p (1 - p) / se^2, corrected for the design's df
  ne <- p * (1 - p) / s^2 * (stats::qt(0.975, 182) / stats::qt(0.975, 14))^2
  expect_equal(b[["neff"]], ne, tolerance = 1e-12)
  expect_gt(ne, 183)
  # truncated at the sample size, as Korn and Graubard have it
  t <- svypci(p, s, 183, df = 14)
  expect_identical(t[["neff"]], 183)
  expect_identical(t[["lo"]], stats::qbeta(0.025, p * 183, 183 - p * 183 + 1))
  expect_identical(t[["hi"]], stats::qbeta(0.975, p * 183 + 1, 183 - p * 183))
  expect_lt(t[["lo"]], b[["lo"]])
})

test_that("svypci() gives the Wilson and Wald intervals written out", {
  ne <- 0.12 * 0.88 / 0.021^2 * (stats::qt(0.975, 399) / stats::qt(0.975, 38))^2
  z <- stats::qnorm(0.975)
  w <- svypci(0.12, 0.021, 400, df = 38, method = "wilson")
  expect_equal(w[["neff"]], ne, tolerance = 1e-12)
  mid <- 0.12 + z^2 / (2 * ne)
  half <- z * sqrt(0.12 * 0.88 / ne + z^2 / (4 * ne^2))
  expect_equal(
    c(w[["lo"]], w[["hi"]]),
    c(mid - half, mid + half) / (1 + z^2 / ne),
    tolerance = 1e-12
  )
  # a simple random sample: the effective size is the sample size, and the
  # interval is the Wilson score interval of the count
  srs <- svypci(0.3, sqrt(0.3 * 0.7 / 200), 200, method = "wilson")
  expect_equal(srs[["neff"]], 200, tolerance = 1e-12)
  pt <- stats::prop.test(60, 200, correct = FALSE)$conf.int
  expect_equal(c(srs[["lo"]], srs[["hi"]]), as.vector(pt), tolerance = 1e-9)
  # and the beta interval is Clopper and Pearson's
  cp <- stats::binom.test(60, 200)$conf.int
  bb <- svypci(0.3, sqrt(0.3 * 0.7 / 200), 200)
  expect_equal(c(bb[["lo"]], bb[["hi"]]), as.vector(cp), tolerance = 1e-9)
  a <- svypci(0.12, 0.021, 400, df = 38, method = "wald")
  expect_identical(
    c(a[["lo"]], a[["hi"]]),
    0.12 + c(-1, 1) * stats::qt(0.975, 38) * 0.021
  )
  # near 0 the Wald interval leaves the range and the others do not
  expect_lt(svypci(0.01, 0.008, 400, df = 38, method = "wald")[["lo"]], 0)
  for (m in c("beta", "wilson", "logit")) {
    expect_gte(svypci(0.01, 0.008, 400, df = 38, method = m)[["lo"]], 0)
  }
  # a wider interval at a higher level
  expect_lt(
    svypci(0.12, 0.021, 400, df = 38, level = 0.99)[["lo"]],
    svypci(0.12, 0.021, 400, df = 38)[["lo"]]
  )
})

test_that("svypci() at 0, at 1, with no variance and with missing values", {
  e <- svypci(c(0, 1, 0.5, NA), 0.02, 100, df = 20)
  expect_identical(e[["neff"]][1:2], c(100, 100))
  expect_identical(e[["lo"]][1L], 0)
  expect_identical(e[["hi"]][1L], stats::qbeta(0.975, 1, 100))
  expect_identical(e[["hi"]][2L], 1)
  expect_identical(e[["lo"]][2L], stats::qbeta(0.025, 100, 1))
  expect_true(all(is.na(e[4L, ])))
  # no standard error: the sample size stands in
  expect_identical(svypci(0.5, 0, 100)[["neff"]], 100)
  lg <- svypci(c(0, 1, 0.5), 0.02, 100, df = 20, method = "logit")
  expect_identical(is.na(lg[["lo"]]), c(TRUE, TRUE, FALSE))
  expect_identical(is.na(lg[["hi"]]), c(TRUE, TRUE, FALSE))
  wl <- svypci(c(0, 1), 0.02, 100, df = 20, method = "wilson")
  expect_identical(wl[["lo"]][1L], 0)
  expect_equal(wl[["hi"]][2L], 1, tolerance = 1e-12)
  expect_true(is.na(svypci(0.5, NA_real_, 100)[["lo"]]))
  # one row per estimate, the others repeated
  expect_identical(nrow(svypci(c(0.1, 0.2, 0.3), 0.02, 100)), 3L)
})

test_that("svypci() refuses bad input", {
  errs(
    "svypci",
    list(
      list("`est` must be numeric", "a", 0.1, 100),
      list("`se` must be numeric", 0.1, NULL, 100),
      list("`n` must be numeric", 0.1, 0.1, "x"),
      list("`df` must be numeric", 0.1, 0.1, 100, df = "x"),
      list("`est` must be proportions, from 0 to 1", 1.2, 0.1, 100),
      list("`se` must be numbers from 0", 0.2, -0.1, 100),
      list("`n` must be whole numbers of at least 2", 0.2, 0.1, 1),
      list("`n` must be whole numbers of at least 2", 0.2, 0.1, 10.5),
      list("`df` must be at least 1", 0.2, 0.1, 100, df = 0),
      list("`level` must be one number in", 0.2, 0.1, 100, level = 1),
      list("`level` must be one number in", 0.2, 0.1, 100, level = "a"),
      list("`trunc` must be TRUE or FALSE", 0.2, 0.1, 100, trunc = NA),
      list(
        "the arguments must have one length, or length 1",
        c(0.1, 0.2),
        c(0.1, 0.2, 0.3),
        100
      ),
      list("`method` must be one of", 0.2, 0.1, 100, method = "exact")
    )
  )
})

test_that("svyqual() grades an estimate by the guide's Table 1", {
  expect_identical(
    svyqual(
      cv = c(8, 20, 40, 8, 8, 16.6, 16.7, 33.3, 33.4),
      n = c(120, 60, 30, 6, 40, 10, 10, 10, 10),
      d = c(900, 500, 300, 200, 45, 50, 50, 50, 50)
    ),
    c(
      "acceptable",
      "marginal",
      "unacceptable",
      "unacceptable",
      "unacceptable",
      "acceptable",
      "marginal",
      "marginal",
      "unacceptable"
    )
  )
  # a mean or a total: the numerator alone
  expect_identical(
    svyqual(c(12, 12), c(35, 9)),
    c("acceptable", "unacceptable")
  )
  # what cannot be known is missing; too few or too variable is still said
  expect_identical(
    svyqual(c(NA, NA, 50, 10), c(50, 5, NA, NA)),
    c(NA, "unacceptable", "unacceptable", NA)
  )
  expect_identical(svyqual(10, 50, c(NA, 40)), c(NA, "unacceptable"))
  # other thresholds
  expect_identical(svyqual(20, 30, cvlo = 25), "acceptable")
  expect_identical(svyqual(20, 30, nmin = 40), "unacceptable")
  expect_identical(svyqual(20, 30, 60, dmin = 100), "unacceptable")
  expect_identical(svyqual(30, 30, cvhi = 25, cvlo = 10), "unacceptable")
  errs(
    "svyqual",
    list(
      list("`cv` must be numeric", "a", 10),
      list("`cv` must be numbers from 0", -1, 10),
      list("`n` must be numbers from 0", 10, -1),
      list("`d` must be numeric", 10, 10, "a"),
      list("`cvlo` must be one number from 0", 10, 10, cvlo = -1),
      list("`dmin` must be one number from 0", 10, 10, dmin = c(1, 2)),
      list("`cvlo` must not be above `cvhi`", 10, 10, cvlo = 40),
      list(
        "the arguments must have one length, or length 1",
        c(1, 2),
        c(1, 2, 3)
      )
    )
  )
})
