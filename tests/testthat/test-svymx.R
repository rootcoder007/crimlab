# svymx(): every statistic, the domains, a contrast and the finite
# population correction against stored values of survey 4.5
# (dev/ref-svy.R); simple random samples by hand; the error against the
# Monte Carlo spread of a stratified cluster sample; and the refusals.

d <- svdata()
sv <- function(...) svymx(d$y, d$w, d$region, d$block, ...)

test_that("the estimates and errors are survey's", {
  want <- list(mean = c(2.88288615553608, 0.136620796168663),
               total = c(53271.0480033886, 2479.77406151403),
               var = c(3.09475888701346, 0.274279470672564))
  for (s in names(want)) {
    r <- sv(stat = s)$est
    expect_equal(c(r$est, r$se), want[[s]], tolerance = 1e-12, info = s)
  }
  r <- sv(stat = "ratio", den = d$x)$est
  expect_equal(c(r$est, r$se), c(0.483989952420717, 0.0292923583424013),
               tolerance = 1e-12)
  p <- svymx(d$v, d$w, d$region, d$block, stat = "prop")$est
  expect_equal(c(p$est, p$se), c(0.285362043550148, 0.0344583306160338),
               tolerance = 1e-12)
  f <- sv(fpc = d$N)$est
  expect_equal(f$se, 0.125328778175561, tolerance = 1e-12)
  m <- sv()
  expect_identical(m$est$df, 15)
  expect_equal(m$est$hi, m$est$est + qt(0.975, 15) * m$est$se,
               tolerance = testthat_tolerance())
  expect_equal(m$est$cv, 100 * m$est$se / m$est$est,
               tolerance = testthat_tolerance())
  expect_output(print(m), "18 first-stage units")
})

test_that("domains and their differences are survey's", {
  r <- sv(by = d$region)
  expect_equal(r$est$est, c(3.215615336862, 3.0416635533909,
                            2.39660807681759), tolerance = 1e-12)
  expect_equal(r$est$se, c(0.209414859731169, 0.265763241998262,
                           0.228242262459905), tolerance = 1e-12)
  expect_equal(r$cmp$diff[1], -0.173951783471097, tolerance = 1e-12)
  expect_equal(r$cmp$se[1], 0.338355854498857, tolerance = 1e-12)
  expect_equal(r$cmp$p, 2 * pt(-abs(r$cmp$stat), 15),
               tolerance = testthat_tolerance())
  s <- sv(by = d$region, ref = "s")
  expect_identical(s$cmp$dom, c("e", "n"))
  expect_output(print(s), "diff")
})

test_that("simple random samples: Bessel, the correction and a census", {
  y <- c(3, 7, 1, 9, 4, 6)
  r <- svymx(y)
  expect_equal(r$est$se, sd(y) / sqrt(6), tolerance = 1e-14)
  expect_equal(svymx(y, stat = "var")$est$est, var(y), tolerance = 1e-14)
  expect_equal(svymx(y, fpc = rep(20, 6))$est$se,
               sqrt((1 - 6 / 20) * var(y) / 6), tolerance = 1e-14)
  expect_identical(svymx(y, fpc = rep(6, 6))$est$se, 0)
  g <- svymx(y, w = c(1, 2, 1, 2, 1, 2))$des
  expect_equal(g$effn, 81 / 15, tolerance = testthat_tolerance())
  expect_equal(g$deffw, 6 * 15 / 81, tolerance = testthat_tolerance())
})

test_that("the error fits the spread of a stratified cluster sample", {
  set.seed(8)
  pop <- data.frame(h = rep(1:4, each = 500), c = rep(1:200, each = 10))
  pop$y <- rnorm(200)[pop$c] + pop$h + rnorm(2000)
  est <- se <- numeric(300)
  for (k in seq_along(est)) {
    cs <- unlist(lapply(1:4, function(h) {
      sample(unique(pop$c[pop$h == h]), 10)
    }))
    s <- pop[pop$c %in% cs, ]
    r <- svymx(s$y, rep(5, nrow(s)), s$h, s$c, fpc = rep(50, nrow(s)))
    est[k] <- r$est$est
    se[k] <- r$est$se
  }
  expect_lt(abs(mean(est) - mean(pop$y)), 4 * sd(est) / sqrt(300))
  expect_lt(abs(mean(se) / sd(est) - 1), 0.1)
})

test_that("svymx() refuses what it cannot estimate", {
  expect_error(svymx(1), "two values")
  expect_error(svymx(c(0, 2), stat = "prop"), "0 and 1")
  expect_error(sv(stat = "ratio"), "needs `den`")
  expect_error(sv(den = d$x), "denominator")
  expect_error(sv(stat = "ratio", den = d$x * 0), "sums to 0")
  expect_error(svymx(d$y, d$w, d$block, d$block), "one first-stage unit")
  expect_error(sv(fpc = seq_len(180)), "constant within")
  expect_error(sv(fpc = rep(2, 180)), "at least the first-stage")
  expect_error(svymx(c(1, 2, 3), by = c("a", "a", "b"), stat = "var"),
               "one unit")
  expect_error(svymx(d$y, -d$w), "`w`")
})
