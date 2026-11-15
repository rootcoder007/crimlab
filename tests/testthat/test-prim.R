# The summary statistics (R/prim.R). The first part is the test file of
# the same functions in rmoriebricklayer, carried over; the last tests
# hold values computed once by rmoriebricklayer 1.1.1 on the same input.

# Statistical kernels in the compiled core.
#
# Anchors outside the package: base R's own implementations
# (stats::sd, var, quantile type 7, median, mad, mean(trim=),
# weighted.mean, cor(method="spearman"), rank, cov, dnorm(log=),
# pgamma, pexp), closed-form moment definitions written out by hand, and
# for the Hawkes likelihood an independent recomputation of the same
# recursion in R.

test_that("csdev and cdistm agree with base R", {
  set.seed(1)
  x <- stats::rnorm(80)
  y <- stats::rnorm(80)
  expect_equal(
    csdev(x),
    stats::sd(x),
    tolerance = testthat::testthat_tolerance()
  )
  expect_identical(csdev(1:10), stats::sd(1:10))
  # ddof = 0 is the population figure
  expect_equal(csdev(x, ddof = 0), sqrt(mean((x - mean(x))^2)),
    tolerance = 1e-12)
  expect_equal(
    csdev(x)^2,
    stats::var(x),
    tolerance = testthat::testthat_tolerance()
  )
  # a constant vector has no spread
  expect_identical(csdev(rep(4, 10)), 0)
  # a single observation has no sample standard deviation to report
  expect_true(is.na(csdev(5)))

  expect_identical(cdistm(c(0, 0), c(3, 4)), 5)
  expect_identical(cdistm(x, y), sqrt(sum((x - y)^2)))
  # a point is zero from itself, and the metric is symmetric
  expect_identical(cdistm(x, x), 0)
  expect_identical(cdistm(x, y), cdistm(y, x))
  # and it satisfies the triangle inequality
  z <- stats::rnorm(80)
  expect_lte(cdistm(x, z), cdistm(x, y) + cdistm(y, z))
  expect_error(cdistm(1:3, 1:4), "same length")
})

test_that("cnormlpdf is dnorm on the log scale", {
  x <- seq(-4, 4, by = 0.25)
  expect_equal(
    cnormlpdf(x),
    stats::dnorm(x, log = TRUE),
    tolerance = testthat::testthat_tolerance()
  )
  expect_equal(
    cnormlpdf(x, 0.3, 1.7),
    stats::dnorm(x, 0.3, 1.7, log = TRUE),
    tolerance = testthat::testthat_tolerance()
  )
  # exponentiating recovers the density
  expect_equal(
    exp(cnormlpdf(x, -1, 2)),
    stats::dnorm(x, -1, 2),
    tolerance = testthat::testthat_tolerance()
  )
  # the peak is -log(sd sqrt(2 pi))
  expect_equal(
    cnormlpdf(5, mean = 5, sd = 2),
    -log(2 * sqrt(2 * pi)),
    tolerance = testthat::testthat_tolerance()
  )
  # it stays finite where the density itself underflows to zero
  expect_identical(stats::dnorm(50), 0)
  expect_true(is.finite(cnormlpdf(50)))
  expect_identical(cnormlpdf(50), stats::dnorm(50, log = TRUE))
  # a non-positive scale is not a distribution: refused since 0.5.13, as
  # cnormpdf() does (it returned NaN)
  expect_error(cnormlpdf(0, 0, 0), "one positive finite number")
  expect_error(cnormlpdf(0, 0, -1), "one positive finite number")
})

test_that("cmom computes all four moments in one pass", {
  set.seed(2)
  x <- stats::rnorm(200, mean = 3, sd = 2)
  m <- cmom(x)
  expect_named(m, c("mean", "variance", "skewness", "kurtosis"))
  expect_identical(m[["mean"]], mean(x))
  expect_equal(m[["variance"]], stats::var(x), tolerance = 1e-12)
  # the shape statistics are the sample-moment definitions
  cm <- function(k) mean((x - mean(x))^k)
  expect_equal(
    m[["skewness"]],
    cm(3) / cm(2)^1.5,
    tolerance = testthat::testthat_tolerance()
  )
  expect_equal(
    m[["kurtosis"]],
    cm(4) / cm(2)^2 - 3,
    tolerance = testthat::testthat_tolerance()
  )

  # a symmetric sample has no skew
  expect_identical(cmom(c(-2, -1, 0, 1, 2))[["skewness"]], 0)
  # kurtosis is reported as EXCESS, so a normal sample sits near zero and
  # a heavy-tailed one above it
  expect_lt(abs(cmom(stats::rnorm(5000))[["kurtosis"]]), 0.4)
  expect_gt(cmom(c(rep(0, 50), -20, 20))[["kurtosis"]], 5)
  # a right tail gives positive skew
  expect_gt(cmom(c(rep(1, 20), 50))[["skewness"]], 0)
  expect_lt(cmom(c(rep(1, 20), -50))[["skewness"]], 0)
  # the one-pass recurrence is numerically stable at a large offset,
  # where a naive sum-of-squares loses the variance entirely
  big <- 1e9 + c(1, 2, 3, 4, 5)
  expect_identical(cmom(big)[["variance"]], stats::var(big))

  # a shape statistic that is not defined is NaN rather than a guess
  expect_true(is.nan(cmom(c(1, 2))[["skewness"]]))
  expect_true(is.nan(cmom(c(1, 2, 3))[["kurtosis"]]))
  expect_true(is.nan(cmom(1)[["variance"]]))
  # a constant vector has no scale to standardise by
  expect_true(is.nan(cmom(rep(2, 10))[["skewness"]]))
  # NA propagates rather than being silently dropped
  expect_true(all(is.nan(cmom(c(1, 2, NA, 4)))))
})

test_that("cqnt is the type-7 quantile", {
  set.seed(3)
  x <- stats::rnorm(97)
  pr <- c(0, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1)
  expect_equal(
    unname(cqnt(x, pr)),
    unname(stats::quantile(x, pr, type = 7)), tolerance = 1e-12)
  expect_named(cqnt(x, c(0.25, 0.5)), c("25%", "50%"))
  # the extremes are the sample extremes
  expect_identical(unname(cqnt(x, 0)), min(x))
  expect_identical(unname(cqnt(x, 1)), max(x))
  # a tiny sample still interpolates the way base R does
  expect_equal(unname(cqnt(c(1, 2), 0.5)), 1.5)
  expect_identical(unname(cqnt(5, c(0, 0.5, 1))), rep(5, 3))
  # quantiles are non-decreasing in p
  expect_true(all(diff(cqnt(x, seq(0, 1, by = 0.05))) >= 0))
  expect_error(cqnt(x, 1.5), "\\[0, 1\\]")
  expect_error(cqnt(x, -0.1), "\\[0, 1\\]")
  expect_error(cqnt(x, numeric(0)), "empty")
  expect_error(cqnt(x, NA), "\\[0, 1\\]")
})

test_that("the robust summaries agree with base R and resist outliers", {
  x <- c(2, 4, 4, 4, 5, 5, 7, 9)
  expect_identical(cmed(x), stats::median(x))
  expect_identical(cmad(x), stats::mad(x))
  expect_identical(cmad(x, const = 1), stats::mad(x, constant = 1))
  expect_identical(ciqr(x), unname(stats::IQR(x)))
  set.seed(4)
  y <- stats::rnorm(150)
  expect_identical(cmed(y), stats::median(y))
  expect_identical(cmad(y), stats::mad(y))
  expect_identical(ciqr(y), unname(stats::IQR(y)))
  # an even-length sample averages the two middle values
  expect_equal(cmed(c(1, 2, 3, 4)), 2.5)

  # the fences are drawn k IQRs beyond the quartiles
  f <- ctkyfnc(x)
  q <- stats::quantile(x, c(0.25, 0.75), type = 7)
  expect_identical(f[["lower"]], unname(q[1] - 1.5 * (q[2] - q[1])))
  expect_identical(f[["upper"]], unname(q[2] + 1.5 * (q[2] - q[1])))
  expect_named(f, c("lower", "upper"))
  # a wider k gives a wider interval
  f3 <- ctkyfnc(x, k = 3)
  expect_lt(f3[["lower"]], f[["lower"]])
  expect_gt(f3[["upper"]], f[["upper"]])
  # and the fences identify the planted outlier
  wild <- c(x, 1000)
  fw <- ctkyfnc(wild)
  expect_identical(wild[wild < fw[["lower"]] | wild > fw[["upper"]]], 1000)

  # robustness: one wild value moves the mean far and the median barely
  expect_gt(abs(mean(wild) - mean(x)), 100)
  expect_lte(abs(cmed(wild) - cmed(x)), 1)
  # NA propagates
  expect_true(is.na(cmed(c(1, NA))))
  expect_true(is.na(cmad(c(1, NA))))
})

test_that("trimmed and winsorized means behave as documented", {
  x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 100)
  for (tr in c(0, 0.05, 0.1, 0.2, 0.3, 0.4)) {
    expect_identical(ctrmmean(x, tr), mean(x, trim = tr))
  }
  # trim = 0 is the plain mean; trim = 0.5 is the median
  expect_identical(ctrmmean(x, 0), mean(x))
  expect_identical(ctrmmean(x, 0.5), cmed(x))
  expect_identical(cwinsmean(x, 0.5), cmed(x))

  # the winsorized mean pulls the tails in rather than dropping them, so
  # it stays below the raw mean when there is one high outlier
  expect_lt(cwinsmean(x, 0.1), mean(x))
  # computed by hand: replace the 100 with 9 and the 1 with 2. On this
  # sample the two coincide, because the trimmed pair happens to average
  # to the same value as the pulled-in pair
  expect_identical(
    cwinsmean(x, 0.1),
    mean(c(2, 2, 3, 4, 5, 6, 7, 8, 9, 9))
  )
  expect_identical(ctrmmean(x, 0.1), cwinsmean(x, 0.1))
  # where the tails are asymmetric they separate: with two high outliers
  # winsorizing keeps some of their weight and trimming discards it
  asym <- c(1, 2, 3, 4, 5, 6, 7, 8, 500, 900)
  expect_gt(cwinsmean(asym, 0.1), ctrmmean(asym, 0.1))
  # on a symmetric sample with no outliers all three agree
  s <- c(1, 2, 3, 4, 5)
  expect_identical(ctrmmean(s, 0.2), mean(s))
  expect_identical(cwinsmean(s, 0.2), mean(s))
  # trimming more can only move the estimate toward the median, though
  # successive trims can land on the same value
  tms <- vapply(
    c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
    function(t) ctrmmean(x, t),
    0
  )
  expect_true(all(diff(tms) <= 0))
  expect_lt(tms[length(tms)], tms[1])

  # a trim outside [0, 0.5] is refused since 0.5.13 (it returned NaN)
  expect_error(ctrmmean(x, 0.6), "from 0 to 0.5")
  expect_error(ctrmmean(x, -0.1), "from 0 to 0.5")
  expect_true(is.na(ctrmmean(c(1, NA), 0.1)))
})

test_that("cwtd matches base R and reduces correctly", {
  x <- c(10, 20, 30, 40)
  w <- c(1, 1, 2, 4)
  r <- cwtd(x, w)
  expect_named(r, c("mean", "variance"))
  expect_identical(r[["mean"]], stats::weighted.mean(x, w))
  # equal weights recover the unweighted moments exactly
  eq <- cwtd(x, rep(1, 4))
  expect_identical(eq[["mean"]], mean(x))
  expect_identical(eq[["variance"]], stats::var(x))
  expect_identical(cwtd(x, rep(3, 4))[["mean"]], mean(x))
  # the weighted variance is the reliability-weight estimator
  sw <- sum(w)
  sw2 <- sum(w^2)
  mu <- sum(w * x) / sw
  expect_identical(r[["variance"]], sum(w * (x - mu)^2) / (sw - sw2 / sw))
  # scaling every weight leaves both moments unchanged
  expect_identical(cwtd(x, w * 7)[["mean"]], r[["mean"]])
  expect_identical(cwtd(x, w * 7)[["variance"]], r[["variance"]])
  # a dominant weight pulls the mean onto that observation
  expect_equal(cwtd(x, c(1, 1, 1, 1e6))[["mean"]], 40, tolerance = 1e-3)
  # a zero weight is the same as dropping the observation
  expect_identical(cwtd(x, c(1, 1, 1, 0))[["mean"]], mean(x[1:3]))

  expect_error(cwtd(1:3, 1:4), "same length")
  expect_error(cwtd(x, c(-1, 1, 1, 1)), "non-negative")
  # all-zero or infinite weights have no mean: refused, not NaN
  expect_error(cwtd(x, rep(0, 4)), "not all zero")
  expect_error(cwtd(x, c(1, Inf, 1, 1)), "finite")
})

test_that("Spearman's rho and the midranks match stats::cor and rank", {
  set.seed(5)
  x <- stats::rnorm(120)
  y <- stats::rnorm(120)
  expect_equal(
    cspr(x, y),
    stats::cor(x, y, method = "spearman"),
    tolerance = testthat::testthat_tolerance()
  )
  # ties get midranks, which is what makes the agreement hold
  xt <- c(1, 2, 2, 2, 5, 5, 7)
  yt <- c(3, 1, 1, 4, 4, 9, 2)
  expect_identical(
    cspr(xt, yt),
    stats::cor(xt, yt, method = "spearman")
  )
  expect_identical(cmidrk(xt), rank(xt))
  expect_identical(cmidrk(yt), rank(yt))
  expect_identical(cmidrk(c(5, 5, 5)), rep(2, 3))
  expect_identical(cmidrk(c(3, 1, 2)), c(3, 1, 2))

  # a monotone relationship is rho = 1 even where it is far from linear
  m <- 1:20
  expect_identical(cspr(m, m^3), 1)
  expect_identical(cspr(m, -exp(m / 5)), -1)
  # invariant to any order-preserving transformation of either variable
  expect_identical(cspr(x, y), cspr(exp(x), y))
  expect_identical(cspr(x, y), cspr(x, y^3 * 2))
  # bounded, and zero for a constant column since there is no ordering
  expect_lte(abs(cspr(x, y)), 1)
  expect_true(is.nan(cspr(x, rep(1, 120))))
  expect_error(cspr(1:3, 1:4), "same length")
})

test_that("ccovm matches stats::cov and keeps its names", {
  set.seed(6)
  X <- cbind(a = stats::rnorm(50), b = stats::rnorm(50), c = stats::rnorm(50))
  expect_equal(
    ccovm(X),
    stats::cov(X),
    tolerance = testthat::testthat_tolerance()
  )
  expect_identical(dimnames(ccovm(X)), list(colnames(X), colnames(X)))
  # the diagonal is the column variances and the matrix is symmetric
  expect_equal(
    diag(ccovm(X)),
    apply(X, 2, stats::var),
    tolerance = testthat::testthat_tolerance()
  )
  expect_identical(ccovm(X), t(ccovm(X)))
  # the off-diagonal is the pairwise covariance
  expect_equal(
    ccovm(X)["a", "b"],
    stats::cov(X[, "a"], X[, "b"]),
    tolerance = testthat::testthat_tolerance()
  )
  # data frames and integer matrices are accepted
  df <- data.frame(u = 1:5, v = c(2, 1, 4, 3, 6))
  expect_identical(ccovm(df), stats::cov(df))
  expect_equal(
    ccovm(cbind(1:5, 5:1)),
    stats::cov(cbind(1:5, 5:1)),
    ignore_attr = TRUE
  )
  # a single column is a 1x1 variance
  expect_identical(as.numeric(ccovm(cbind(z = 1:6))), stats::var(1:6))

  expect_error(ccovm(1:5), "matrix or data frame")
  expect_error(
    ccovm(data.frame(a = 1:3, b = letters[1:3])),
    "must be numeric"
  )
})

test_that("bootstrap replicate means are reproducible and unbiased", {
  set.seed(7)
  x <- stats::rnorm(60, mean = 5)
  b <- cbootmean(x, B = 3000, seed = 11)
  expect_length(b, 3000L)
  # the same seed reproduces the replicates exactly
  expect_identical(cbootmean(x, 3000, seed = 11), b)
  expect_false(identical(cbootmean(x, 3000, seed = 12), b))
  # the replicates centre on the sample mean, and their spread estimates
  # the standard error. The bootstrap SE is itself a random quantity, so
  # this is a loose agreement by nature, not an identity
  expect_equal(mean(b), mean(x), tolerance = 0.02)
  expect_equal(stats::sd(b), stats::sd(x) / sqrt(length(x)), tolerance = 0.1)
  # every replicate is a mean of values drawn from x, so all lie inside
  # the sample range
  expect_true(all(b >= min(x) & b <= max(x)))
  # the core uses its own generator, so R's stream is untouched
  set.seed(99)
  before <- stats::runif(1)
  set.seed(99)
  invisible(cbootmean(x, 50, seed = 3))
  expect_identical(stats::runif(1), before)
  # a constant sample bootstraps to that constant
  expect_identical(unique(cbootmean(rep(2, 10), 20, seed = 1)), 2)
  expect_error(cbootmean(x, 0), "whole number from 1")
})

test_that("IPW weights invert the clamped propensity score", {
  treat <- c(1, 0, 1, 0)
  e <- c(0.5, 0.25, 0.02, 0.9)
  expect_equal(
    cipw(treat, e),
    c(1 / 0.5, 1 / 0.75, 1 / 0.02, 1 / 0.1),
    tolerance = testthat::testthat_tolerance()
  )
  # a balanced score gives weight 2 to either arm
  expect_identical(cipw(c(1, 0), c(0.5, 0.5)), c(2, 2))
  # the clamp bites before the inversion, so it bounds the weight
  expect_identical(cipw(1, 0.001, trimlo = 0.05), 1 / 0.05)
  expect_equal(
    cipw(0, 0.999, trimhi = 0.95),
    1 / 0.05,
    tolerance = testthat::testthat_tolerance()
  )
  # a tighter clamp can only shrink an extreme weight
  expect_lt(
    cipw(1, 0.02, trimlo = 0.10),
    cipw(1, 0.02, trimlo = 0.01)
  )
  # logical indicators work
  expect_identical(
    cipw(c(TRUE, FALSE), c(0.4, 0.4)),
    c(1 / 0.4, 1 / 0.6)
  )
  # every weight is at least 1, since a probability cannot exceed one
  expect_true(all(cipw(treat, e) >= 1))
  expect_error(cipw(1:3, 1:4), "same length")
  expect_error(cipw(1, 0.5, trimlo = 0), "trimlo")
  expect_error(
    cipw(1, 0.5, trimlo = 0.9, trimhi = 0.1),
    "trimlo"
  )
})

test_that("cgamcdf is the unit-rate gamma distribution function", {
  q <- c(0.01, 0.5, 1, 2, 4, 12, 40)
  for (a in c(0.5, 1, 2, 3.5, 10)) {
    expect_equal(
      cgamcdf(a, q),
      stats::pgamma(q, shape = a),
      tolerance = testthat::testthat_tolerance()
    )
  }
  # shape 1 is the exponential
  expect_equal(
    cgamcdf(1, q),
    stats::pexp(q),
    tolerance = testthat::testthat_tolerance()
  )
  # it is a distribution function: 0 at the origin, rising to 1
  expect_identical(cgamcdf(2, 0), 0)
  expect_identical(cgamcdf(2, 1e4), 1)
  expect_true(all(diff(cgamcdf(3, seq(0, 20, by = 0.5))) >= 0))
  expect_true(all(cgamcdf(3, q) >= 0 & cgamcdf(3, q) <= 1))
  # the chi-square tail built from it matches stats::pchisq
  for (df in c(1, 2, 5, 8)) {
    for (stat in c(0.5, 2, 6, 15)) {
      expect_equal(
        1 - cgamcdf(df / 2, stat / 2),
        stats::pchisq(stat, df, lower.tail = FALSE),
        tolerance = testthat::testthat_tolerance()
      )
    }
  }
  expect_error(cgamcdf(0, 1), "one positive")
  expect_error(cgamcdf(c(1, 2), 1), "one positive")
})

test_that("core stats match base R", {
  x <- c(2, 4, 4, 4, 5, 5, 7, 9)
  y <- (seq_along(x))^2
  expect_identical(cmn(x), mean(x))
  expect_identical(cvari(x), var(x))
  expect_identical(ccorr(x, y), cor(x, y))
  expect_equal(
    cnormpdf(c(-1, 0, 1, 2.5), mean = 0.3, sd = 1.7),
    dnorm(c(-1, 0, 1, 2.5), mean = 0.3, sd = 1.7),
    tolerance = testthat::testthat_tolerance()
  )
})

test_that("core stats handle degenerate input the way R's contract implies", {
  expect_true(is.nan(cvari(1))) # sample variance needs n >= 2
  expect_true(is.nan(ccorr(c(1, 1, 1), c(2, 3, 4)))) # zero variance
  # sd must be > 0: refused, not a silent NaN
  expect_error(cnormpdf(0, sd = 0), "positive")
})

test_that("ccorr rejects length-mismatched input", {
  expect_error(ccorr(1:3, 1:4), "same length")
})

test_that("the values equal those of rmoriebricklayer", {
  # rmoriebricklayer 1.1.1 on these inputs, printed with sprintf("%.17g")
  x <- c(2, 4, 4, 4, 5, 5, 7, 9, 1.5, 3.25)
  y <- c(1, 3, 2, 5, 4, 6, 8, 7, 9, 10)
  t <- 1e-15
  expect_equal(cmn(x), 4.4749999999999996, tolerance = t)
  expect_equal(cvari(x), 4.9506944444444452, tolerance = t)
  expect_equal(ccorr(x, y), 0.17936924302140742, tolerance = t)
  expect_equal(csdev(x, 0), 2.1108351427811693, tolerance = t)
  expect_equal(cdistm(x, y), 10.761621625015442, tolerance = t)
  expect_equal(cnormpdf(c(-1, 0, 2.5), 0.3, 1.7),
               c(0.17517771124055995, 0.23104616994892685,
                 0.10157711386142985), tolerance = t)
  expect_equal(cnormlpdf(c(-1, 0, 50), 0.3, 1.7),
               c(-1.741954327519438, -1.4651377185228982,
                 -428.80077785693129), tolerance = t)
  expect_equal(unname(cmom(x)),
               c(4.4749999999999996, 4.9506944444444434,
                 0.69524854265814851, -0.086964115068321066), tolerance = t)
  expect_equal(unname(cqnt(x, c(0, 0.1, 0.33, 0.5, 0.9, 1))),
               c(1.5, 1.95, 3.9775, 4, 7.1999999999999993, 9), tolerance = t)
  expect_identical(cmed(x), 4)
  expect_equal(cmad(x), 1.4825999999999999, tolerance = t)
  expect_equal(ciqr(x), 1.5625, tolerance = t)
  expect_equal(unname(ctkyfnc(x, 2)), c(0.3125, 8.125), tolerance = t)
  expect_equal(ctrmmean(c(x, 100), 0.2), 4.6071428571428568, tolerance = t)
  expect_equal(cwinsmean(c(x, 100), 0.2), 4.7954545454545459, tolerance = t)
  expect_equal(unname(cwtd(x, y)), c(4.6727272727272728, 6.2080018939393939),
               tolerance = t)
  expect_equal(cspr(x, y), 0.067700320038633008, tolerance = t)
  expect_identical(cmidrk(x), c(2, 5, 5, 5, 7.5, 7.5, 9, 10, 1, 3))
  expect_equal(as.vector(ccovm(cbind(x, y))),
               c(4.9506944444444452, 1.2083333333333333, 1.2083333333333333,
                 9.1666666666666661), tolerance = t)
  expect_equal(cbootmean(x, 6, 7),
               c(4.1500000000000004, 3.9500000000000002, 4.5,
                 5.0499999999999998, 4.3250000000000002, 4.4000000000000004),
               tolerance = t)
  expect_equal(cipw(c(1, 0, 1, 0), c(0.5, 0.25, 0.002, 0.995)),
               c(2, 1.3333333333333333, 100, 99.999999999999915),
               tolerance = t)
  # crimlabs own gamma kernel, another algorithm than rmoriebricklayers
  expect_equal(cgamcdf(3.5, c(0.5, 1, 4, 12, 0, -1)),
               c(0.0051714634834844359, 0.040159631269898174,
                 0.66740609740059931, 0.99886064882105241, 0, 0),
               tolerance = 1e-12)
  expect_identical(cunif(4, 42),
                   c(0.74156487877182331, 0.1599103928769201,
                     0.27860113025513866, 0.34419071652363753))
  expect_identical(cunif(2, 0), c(0.88331080821364261, 0.43152799704850997))
})

test_that("the refusals name the argument", {
  expect_error(cmn("a"), "cmn: `x` must be numbers")
  expect_error(csdev(1:3, -1), "`ddof` must be one whole number from 0")
  expect_error(cmad(1:3, -1), "`const` must be one number from 0")
  expect_error(ctkyfnc(1:3, NA), "`k` must be one number")
  expect_error(cipw(c(1, 2), c(0.5, 0.5)), "coded 0/1")
  expect_error(cipw(c(1, 0), c(0.5, 1.5)), "`ps` must lie")
  expect_error(cunif(1.5, 1), "`n` must be one whole number")
  expect_error(cunif(2, 2^53), "below 2\\^53")
  expect_error(ccovm(matrix(letters[1:4], 2)), "numeric matrix")
  expect_identical(cgamcdf(2, c(NA, Inf)), c(NA, 1))
  expect_identical(length(cbootmean(numeric(), 3)), 3L)
  expect_true(all(is.nan(cbootmean(numeric(), 3))))
  expect_identical(cmn(NULL), NaN)
})
