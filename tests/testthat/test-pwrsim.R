test_that("pwrsim() with no effect rejects at about the level", {
  r <- pwrsim(20, eff = 0, nsim = 4000, seed = 3)
  t <- r$tab
  expect_identical(names(t), c("n", "n1", "eff", "power", "se", "lo", "hi",
                               "rej"))
  expect_identical(t$power, t$rej / 4000)
  expect_lt(abs(t$power - 0.05), 4 * sqrt(0.05 * 0.95 / 4000))
  # the Monte Carlo standard error and Wilson's interval, by hand
  p <- t$power
  z <- stats::qnorm(0.975)
  mid <- (p + z^2 / 8000) / (1 + z^2 / 4000)
  hw <- z * sqrt(p * (1 - p) / 4000 + z^2 / (4 * 4000^2)) / (1 + z^2 / 4000)
  expect_equal(t$se, sqrt(p * (1 - p) / 4000), tolerance = 1e-14)
  expect_equal(c(t$lo, t$hi), c(mid - hw, mid + hw), tolerance = 1e-12)
})

test_that("pwrsim()'s Welch test has the power power.t.test() gives", {
  want <- stats::power.t.test(n = 64, delta = 0.5, sd = 2)$power
  r <- pwrsim(64, eff = 0.5, sd = 2, nsim = 3000, seed = 8)$tab
  expect_lt(abs(r$power - want), 4 * sqrt(want * (1 - want) / 3000) + 0.005)
  g <- pwrsim(64, eff = 1, alt = "greater", nsim = 2000, seed = 8)$tab
  w1 <- stats::power.t.test(n = 64, delta = 1, sd = 1,
                            alternative = "one.sided")$power
  expect_lt(abs(g$power - w1), 4 * sqrt(w1 * (1 - w1) / 2000) + 0.005)
  l <- pwrsim(64, eff = 1, alt = "less", nsim = 2000, seed = 8)$tab
  expect_identical(l$rej, 0)
})

test_that("pwrsim() draws binary and count outcomes", {
  # every treated unit 1, every control 0: always rejected
  b <- pwrsim(10, eff = 1, fam = "binary", p0 = 0, nsim = 50, seed = 1)$tab
  expect_identical(b$power, 1)
  expect_identical(b$rej, 50)
  # no outcome at all: never rejected
  z <- pwrsim(10, eff = 0, fam = "binary", p0 = 0, nsim = 50, seed = 1)$tab
  expect_identical(z$rej, 0)
  q <- pwrsim(30, eff = 0, fam = "poisson", mu0 = 3, nsim = 3000,
              seed = 5)$tab
  expect_lt(abs(q$power - 0.05), 4 * sqrt(0.05 * 0.95 / 3000) + 0.005)
})

test_that("pwrsim()'s permutation test holds its level", {
  r <- pwrsim(15, eff = 0, test = "perm", nperm = 99, nsim = 1500,
              seed = 2)$tab
  expect_lt(abs(r$power - 0.05), 4 * sqrt(0.05 * 0.95 / 1500) + 0.005)
  s <- pwrsim(15, eff = 3, test = "perm", nperm = 99, nsim = 50,
              seed = 2)$tab
  expect_identical(s$power, 1)
})

test_that("pwrsim() lays out its grid and keeps each design's stream", {
  r <- pwrsim(c(10, 20), eff = c(0, 0.5), ratio = 1.5, nsim = 200,
              seed = 4)
  t <- r$tab
  expect_identical(t$n, c(10, 20, 10, 20))
  expect_identical(t$n1, c(15, 30, 15, 30))
  expect_identical(t$eff, c(0, 0, 0.5, 0.5))
  # the first design alone has the same stream
  one <- pwrsim(10, eff = 0, ratio = 1.5, nsim = 200, seed = 4)$tab
  expect_identical(one$rej, t$rej[1L])
  set.seed(9)
  a <- pwrsim(10, eff = 0.2, nsim = 100)
  expect_identical(pwrsim(10, eff = 0.2, nsim = 100, seed = a$seed), a)
  out <- utils::capture.output(v <- print(r))
  expect_identical(out[1L], paste("pwrsim: normal outcome, welch test",
                                  "(two.sided) at alpha 0.05, 200 trials",
                                  "per design"))
  expect_identical(v, r)
})

test_that("pwrsim() refuses designs it cannot draw", {
  expect_error(pwrsim(1, eff = 0, seed = 1),
               "pwrsim: `n` must be between 2 and 1e+06.", fixed = TRUE)
  expect_error(pwrsim(10, eff = 0.6, fam = "binary", p0 = 0.5, seed = 1),
               "pwrsim: `p0 + eff` must be from 0 to 1.", fixed = TRUE)
  expect_error(pwrsim(10, eff = -2, fam = "poisson", mu0 = 1, seed = 1),
               "pwrsim: `mu0 + eff` must be at least 0.", fixed = TRUE)
  expect_error(pwrsim(10, eff = 0, sd = 0, seed = 1),
               "pwrsim: `sd` must be positive.", fixed = TRUE)
  expect_error(pwrsim(10, eff = 0, alpha = 0, seed = 1),
               "pwrsim: `alpha` must be one number strictly between 0 and 1.",
               fixed = TRUE)
  expect_error(pwrsim(10, eff = 0, test = "perm", nperm = 9, seed = 1),
               paste("pwrsim: `nperm` is too few for `alpha`: the least",
                     "p-value, 1 / (nperm + 1), is above it."), fixed = TRUE)
  expect_error(pwrsim(2, eff = 0, ratio = 0.5, seed = 1),
               paste("pwrsim: `ratio * n` must round to treated arms of 2",
                     "to 1000000."), fixed = TRUE)
  expect_error(pwrsim(10, eff = 0, fam = "gamma", seed = 1),
               paste("pwrsim: `fam` must be one of \"normal\", \"binary\",",
                     "\"poisson\", not \"gamma\"."), fixed = TRUE)
})
