test_that("pwranv() is the noncentral F upper tail at the critical value", {
  r <- pwranv(4, n = c(10, 25), f = c(0.1, 0.3), alpha = 0.05)
  expect_named(r, c("k", "n", "f", "df1", "df2", "ncp", "fcrit",
                    "power"))
  expect_identical(r$n, c(10, 25, 10, 25))
  expect_identical(r$f, c(0.1, 0.1, 0.3, 0.3))
  d2 <- 4 * (r$n - 1)
  l <- 4 * r$n * r$f^2
  fc <- stats::qf(0.05, 3, d2, lower.tail = FALSE)
  expect_identical(r$df1, rep(3, 4))
  expect_identical(r$df2, d2)
  expect_equal(r$ncp, l, tolerance = 1e-15)
  expect_equal(r$fcrit, fc, tolerance = 1e-10)
  expect_equal(r$power, stats::pf(fc, 3, d2, ncp = l, lower.tail = FALSE),
               tolerance = 1e-8)
})

test_that("pwranv() agrees with power.anova.test()", {
  mu <- c(0, 0, 0.5)
  f <- sqrt(mean((mu - mean(mu))^2))
  want <- stats::power.anova.test(groups = 3, n = 20, between.var = var(mu),
                                  within.var = 1)$power
  expect_equal(pwranv(3, n = 20, f = f)$power, want, tolerance = 1e-8)
  # with no effect the power is the level
  expect_equal(pwranv(3, n = 20, f = 0, alpha = 0.1)$power, 0.1,
               tolerance = 1e-10)
})

test_that("pwranv() finds the smallest group size that reaches the power", {
  mu <- c(0, 0, 0.5)
  f <- sqrt(mean((mu - mean(mu))^2))
  nc <- stats::power.anova.test(groups = 3, between.var = var(mu),
                                within.var = 1, power = 0.8)$n
  r <- pwranv(3, f = c(f, 0), power = 0.8)
  expect_identical(r$n, c(ceiling(nc), NA))
  pw <- pwranv(3, n = ceiling(nc) - c(1, 0), f = f)$power
  expect_lt(pw[1L], 0.8)
  expect_gte(pw[2L], 0.8)
  expect_equal(r$power[1L], pw[2L], tolerance = 1e-15)
  expect_identical(r$df1, c(2, 2))
  expect_true(all(is.na(r[2L, c("df2", "ncp", "fcrit", "power")])))
})

test_that("pwranv() refuses sizes, levels and powers it cannot use", {
  expect_error(pwranv(1, n = 10, f = 0.2),
               "pwranv: `k` must be between 2 and 1e+06.", fixed = TRUE)
  expect_error(pwranv(3, f = 0.25), "pwranv: give one of `n` and `power`.",
               fixed = TRUE)
  expect_error(pwranv(3, n = 10, f = 0.25, power = 0.8),
               "pwranv: give one of `n` and `power`.", fixed = TRUE)
  expect_error(pwranv(3, n = 1, f = 0.25), "pwranv: `n` must be above 1.",
               fixed = TRUE)
  expect_error(pwranv(3, n = 10, f = -0.1),
               "pwranv: `f` must be at least 0.", fixed = TRUE)
  expect_error(pwranv(3, n = 10, f = 0.2, alpha = 1),
               "pwranv: `alpha` must be one number strictly between 0 and 1.",
               fixed = TRUE)
  expect_error(pwranv(3, f = 0.2, power = 0.01),
               "pwranv: `power` must be above `alpha` and below 1.",
               fixed = TRUE)
})
