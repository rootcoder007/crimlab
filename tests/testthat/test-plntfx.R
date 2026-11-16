test_that("plntfx() recovers an effect exactly from a constant outcome", {
  set.seed(51)
  d <- data.frame(t = rbinom(40, 1, 0.5), age = runif(40, 18, 60), y = 2)
  r <- plntfx(d, "y", "t", x = "age", eff = c(0, 1.5), nrep = 20, seed = 2)
  tab <- r$tab
  expect_identical(names(tab), c("eff", "mean", "bias", "sd", "rmse",
                                 "cover", "rej", "tol", "ok"))
  # shuffling a constant changes nothing: every estimate is eff
  expect_equal(tab$mean, c(0, 1.5), tolerance = 1e-10)
  expect_lt(max(abs(tab$bias)), 1e-10)
  expect_lt(max(tab$sd, tab$rmse), 1e-10)
})

test_that("plntfx()'s estimates centre on the planted effect", {
  set.seed(52)
  d <- data.frame(t = rbinom(150, 1, 0.4), y = rexp(150))
  r <- plntfx(d, "y", "t", eff = c(0, 0.5), nrep = 400, se = "iid",
              seed = 7)
  tab <- r$tab
  expect_equal(tab$bias, tab$mean - tab$eff, tolerance = 1e-15)
  expect_true(all(abs(tab$bias) < 4 * tab$sd / sqrt(400)))
  expect_true(all(abs(tab$cover - 0.95) < 4 * sqrt(0.95 * 0.05 / 400)))
  expect_lt(abs(tab$rej[1L] - 0.05), 4 * sqrt(0.95 * 0.05 / 400))
  # the root mean squared error splits into bias and spread
  expect_equal(tab$rmse^2, tab$bias^2 + tab$sd^2 * 399 / 400,
               tolerance = 1e-10)
  expect_equal(tab$tol, 3 * tab$sd / sqrt(400), tolerance = 1e-15)
  expect_identical(tab$ok, abs(tab$bias) <= tab$tol)
  # a planted effect moves every estimate by the same amount: the spread
  # is the same for both effects, the shuffles being shared
  expect_equal(tab$sd[1L], tab$sd[2L], tolerance = 1e-10)
  expect_equal(tab$mean[2L] - tab$mean[1L], 0.5, tolerance = 1e-10)
})

test_that("plntfx() keeps its seed and takes a tolerance", {
  set.seed(53)
  d <- data.frame(t = rbinom(60, 1, 0.5), y = rnorm(60))
  r <- plntfx(d, "y", "t", eff = 1, nrep = 30, tol = 0, seed = 3)
  expect_identical(r$tab$tol, 0)
  expect_identical(r$tab$ok, abs(r$tab$bias) <= 0)
  expect_identical(plntfx(d, "y", "t", eff = 1, nrep = 30, tol = 0,
                          seed = 3), r)
  set.seed(2)
  a <- plntfx(d, "y", "t", eff = 1, nrep = 10)
  expect_identical(plntfx(d, "y", "t", eff = 1, nrep = 10, seed = a$seed), a)
  out <- utils::capture.output(v <- print(r))
  expect_identical(out[1L], paste0(
    "plntfx: ", sum(r$tab$ok), " of 1 planted effect(s) recovered within ",
    "tolerance over 30 shuffles"
  ))
  expect_identical(v, r)
})

test_that("plntfx() refuses what it cannot plant", {
  d <- data.frame(t = c(0, 1, 0, 1, 1, 0), y = 1:6, z = 1)
  expect_error(plntfx(d, "y", "y", eff = 1, seed = 1),
               "plntfx: `y`, `d` and `x` must name different columns.",
               fixed = TRUE)
  expect_error(plntfx(d, "y", "z", eff = 1, seed = 1),
               "plntfx: column `z` must vary.", fixed = TRUE)
  expect_error(plntfx(d, "y", "t", eff = 1, nrep = 1, seed = 1),
               "plntfx: `nrep` must be between 2 and 1e+06.", fixed = TRUE)
  expect_error(plntfx(d, "y", "t", eff = 1, tol = -1, seed = 1),
               "plntfx: `tol` must be at least 0.", fixed = TRUE)
  expect_error(plntfx(d, "y", "t", x = "z", eff = 1, seed = 1),
               paste("plntfx: the covariates are collinear, or `d` does not",
                     "vary once they are taken out."), fixed = TRUE)
})
