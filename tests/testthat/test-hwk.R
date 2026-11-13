# hwkfit(), hwkres(), hwksim(): the log-likelihood against a direct sum
# written here, the fit against stats::optim() on that sum, the residuals'
# test against stats::ks.test(), recovery of known parameters, the
# simulation's mean count, and the refusals.

# the log-likelihood, written out: every pair, every kernel
llr <- function(t, hor, th, kernel = "exp") {
  mu <- exp(th[1])
  eta <- th[2]
  g <- switch(kernel,
    exp = function(d) th[3] * exp(-th[3] * d),
    weibull = function(d) dweibull(d, th[3], th[4]),
    lomax = function(d) th[3] / th[4] * (1 + d / th[4])^(-th[3] - 1)
  )
  cg <- switch(kernel,
    exp = function(d) pexp(d, th[3]),
    weibull = function(d) pweibull(d, th[3], th[4]),
    lomax = function(d) 1 - (1 + d / th[4])^(-th[3])
  )
  lam <- vapply(seq_along(t), function(i) {
    e <- t[t < t[i]]
    mu + eta * sum(g(t[i] - e))
  }, 0)
  unname(sum(log(lam)) - mu * hor - eta * sum(cg(hor - t)))
}

set.seed(11)
ev <- hwksim(c(a0 = log(0.6), eta = 0.45, beta = 1.2), 300, seed = 4)$t

test_that("the log-likelihood is the direct sum, every kernel", {
  for (k in c("exp", "weibull", "lomax")) {
    f <- hwkfit(ev, 300, kernel = k)
    expect_equal(f$loglik, llr(ev, 300, f$par, k), tolerance = 1e-10)
    expect_identical(f$kernel, k)
  }
})

test_that("the fit is the maximum stats::optim() finds", {
  f <- hwkfit(ev, 300)
  o <- optim(c(log(0.5), 0.3, 1), function(p) -llr(ev, 300, p),
             method = "L-BFGS-B", lower = c(-5, 1e-6, 1e-3),
             upper = c(5, 0.999, 50), control = list(factr = 1, pgtol = 0))
  expect_equal(unname(f$par), o$par, tolerance = 1e-5)
  expect_equal(f$loglik, -o$value, tolerance = 1e-10)
  expect_true(f$conv)
  # the covariance: the inverse of optimHess() of the same sum
  h <- optimHess(f$par, function(p) -llr(ev, 300, p))
  expect_equal(unname(f$cov), unname(solve(h)), tolerance = 1e-4)
  expect_equal(f$aic, 6 - 2 * f$loglik, tolerance = testthat_tolerance())
  expect_equal(f$bic, 3 * log(length(ev)) - 2 * f$loglik,
               tolerance = testthat_tolerance())
})

test_that("known parameters come back from a long simulation", {
  x <- hwksim(c(log(1), 0.6, 2), 5000, seed = 9)
  f <- hwkfit(x$t, 5000)
  z <- (f$par - c(0, 0.6, 2)) / f$se
  expect_true(all(abs(z) < 4))
  expect_equal(f$branch[["size"]], 1 / (1 - f$par[["eta"]]),
               tolerance = testthat_tolerance())
  expect_equal(f$branch[["rate"]], exp(f$par[["a0"]]) / (1 - f$par[["eta"]]),
               tolerance = testthat_tolerance())
  expect_equal(f$branch[["share"]], f$par[["eta"]],
               tolerance = testthat_tolerance())
  # the simulated count: mu hor / (1 - eta), to a few per cent
  expect_equal(nrow(x), 5000 / 0.4, tolerance = 0.05)
})

test_that("kernels and baselines are compared, the least AIC kept", {
  f <- hwkfit(ev, 300, kernel = c("exp", "lomax"), base = c("const", "cycle"),
              per = 50)
  expect_identical(nrow(f$cmp), 4L)
  expect_identical(f$aic, min(f$cmp$aic))
  expect_identical(f$cmp$k, c(3L, 6L, 4L, 7L))
  expect_output(print(f), "branching ratio")
  expect_output(print(hwkfit(ev, 300, maxit = 1)), "stopped before")
  expect_named(hwkfit(ev, 300, base = "cycle", per = 60)$par,
               c("a0", "a1", "a2", "a3", "eta", "beta"))
})

test_that("a start and a resolution are used", {
  f <- hwkfit(ev, 300, kernel = "weibull", start = c(0, 0.3, 1, 1))
  expect_named(f$par, c("a0", "eta", "shape", "scale"))
  d <- floor(ev)
  expect_warning(hwkfit(d, 300), "tied")
  j <- hwkfit(d, 300, res = 1, seed = 2)
  expect_true(j$jit)
  expect_false(anyDuplicated(j$times) > 0)
  expect_true(all(j$times >= d & j$times < d + 1))
  expect_identical(j$times, hwkfit(d, 300, res = 1, seed = 2)$times)
})

test_that("the residuals' test is ks.test()'s", {
  f <- hwkfit(ev[1:80], ev[80])
  r <- hwkres(f)
  k <- ks.test(r$u, "punif")
  expect_equal(r$d, unname(k$statistic), tolerance = 1e-14)
  expect_equal(r$p, k$p.value, tolerance = 1e-12)
  r2 <- hwkres(ev, 300, f$par)
  k2 <- ks.test(r2$u, "punif")
  expect_equal(r2$p, k2$p.value, tolerance = 1e-12)
  expect_equal(r2$u, 1 - exp(-r2$tau), tolerance = testthat_tolerance())
  # eta = 0: the gaps times the rate
  r0 <- hwkres(ev, 300, c(log(2), 0, 1))
  expect_equal(r0$tau, 2 * diff(c(0, ev)), tolerance = 1e-12)
  expect_output(print(r), "Kolmogorov-Smirnov")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_identical(plot(r), r)
})

test_that("simulation is seeded, in time and in space", {
  a <- hwksim(c(log(0.2), 0.5, 1), 100, seed = 3)
  expect_identical(a, hwksim(c(log(0.2), 0.5, 1), 100, seed = 3))
  expect_false(is.unsorted(a$t))
  set.seed(1)
  b <- hwksim(c(log(0.2), 0.5, 1), 100)
  set.seed(1)
  expect_identical(b, hwksim(c(log(0.2), 0.5, 1), 100))
  s <- hwksim(c(log(0.05), 0.5, 1), 100, win = c(0, 4, 0, 5), sigma = 0.2,
              seed = 1)
  expect_named(s, c("t", "x", "y", "gen"))
  expect_true(all(s$x >= 0 & s$x <= 4 & s$y >= 0 & s$y <= 5))
  # the background is per unit area: 0.05 x 20 x 100 = 100 on average
  expect_equal(sum(s$gen == 0), 100, tolerance = 0.3)
  f <- hwkfit(ev, 300)
  expect_named(hwksim(f, 50, seed = 1), c("t", "gen"))
})

test_that("bad input is refused", {
  errs("hwkfit", list(
    list("at least 10", c(1, 2), 3),
    list("between 0 and `hor`", ev, 100),
    list("between 0 and `hor`", c(-1, ev), 300, res = 1),
    list("sorted", rev(ev), 300),
    list("`kernel` must be one of", ev, 300, kernel = 1),
    list("`kernel` must be one of", ev, 300, kernel = "power"),
    list("`base` must be one of", ev, 300, base = NA_character_),
    list("one kernel and one baseline", ev, 300, kernel = c("exp", "lomax"),
         start = c(0, 0.5, 1)),
    list("`start`: eta must be at least 0", ev, 300, start = c(0, -1, 1)),
    list("eta must be below 1", ev, 300, start = c(0, 1, 1)),
    list("kernel parameters must be positive", ev, 300, start = c(0, 0.5, 0)),
    list("`start` must have length 3", ev, 300, start = c(0, 0.5)),
    list("`hor` must be", ev, -1),
    list("`maxit`", ev, 300, maxit = 0)
  ))
  errs("hwksim", list(
    list("more than `nmax`", c(log(5), 0.5, 1), 100, nmax = 10),
    list("`sigma`", c(0, 0.5, 1), 10, win = c(0, 1, 0, 1))
  ))
  errs("hwkres", list(list("`par` must have length 4", ev, 300, c(0, 0.5, 1),
                           kernel = "lomax")))
})

test_that("a near-critical fit is warned about", {
  # a rate that climbs steeply: a constant baseline reads it as contagion
  t <- sqrt(seq(1, 4e4, length.out = 400))
  expect_warning(f <- hwkfit(t, 200), "near 1")
  expect_gt(f$branch[["eta"]], 0.99)
})
