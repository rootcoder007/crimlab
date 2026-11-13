# detrsp(), dtrid(), inspg(), insphs() and colarr(): against optima, ranks,
# least squares, integration and odds ratios computed here in base R, and
# the refusals.

test_that("detrsp() picks the largest optimal level", {
  set.seed(8)
  x <- sample(0:15)
  ben <- 4 * sqrt(x) + rnorm(16, sd = 0.3)
  sanc <- x^1.4
  p <- c(0, 0.05, 0.2, 0.5, 1, 3)
  out <- detrsp(x, ben, sanc, p)
  ref <- vapply(p, function(q) {
    u <- ben - q * sanc
    c(max(x[u >= max(u) - 1e-12]), max(u))
  }, c(0, 0))
  expect_identical(out$p, p)
  expect_identical(out$xopt, ref[1, ])
  expect_identical(out$value, ref[2, ])
  expect_true(all(diff(out$xopt) <= 0))
})

test_that("dtrid() gives the rank, the identified effects and the fit", {
  a <- dtrid(p = c(0.1, 0.3, 0.5), s = c(2, 2, 2), c = c(30, 30, 30))
  expect_identical(a$rank, 2L)
  expect_identical(a$ident,
                   c(certainty = TRUE, severity = FALSE, celerity = FALSE))
  expect_identical(a$const,
                   c(certainty = FALSE, severity = TRUE, celerity = TRUE))
  expect_null(a$coef)
  expect_null(a$cost)
  set.seed(2)
  p <- runif(40, 0.1, 0.9)
  s <- runif(40, 1, 10)
  c <- runif(40, 0, 12)
  y <- 2 - 0.8 * log(p) - 0.8 * log(s) + 0.08 * c + rnorm(40, sd = 0.1)
  b <- dtrid(p, s, c, y, log = TRUE, rho = 0.1)
  ref <- unname(coef(lm(y ~ log(p) + log(s) + c)))
  expect_identical(b$rank, 4L)
  expect_identical(names(b$coef), c("(Intercept)", "log_p", "log_s", "c"))
  expect_equal(unname(b$coef), ref, tolerance = 1e-12)
  expect_equal(b$disc, -ref[4] / ref[2], tolerance = 1e-12)
  expect_equal(b$cost, p * s * exp(-0.1 * c), tolerance = 1e-15)
  q <- dtrid(p, s, c, y)
  expect_equal(unname(q$coef), unname(coef(lm(y ~ p + s + c))),
               tolerance = 1e-12)
  expect_identical(names(q$coef), c("(Intercept)", "p", "s", "c"))
  expect_true(is.na(q$disc))
  # a dimension that copies another: neither is identified
  d <- dtrid(p, 2 * p, c)
  expect_identical(d$rank, 3L)
  expect_identical(d$ident,
                   c(certainty = FALSE, severity = FALSE, celerity = TRUE))
  z <- dtrid(c(1, 1, 2, 3), c(1, 2, 3, 4), 1:4,
             c(1, 2, 3, 4), log = TRUE)
  expect_true(is.na(z$disc))
})

test_that("inspg() integrates the cycle and keeps its invariant", {
  out <- inspg(1.5, c(0.5, 0.3, 0.2), steps = 300, dt = 0.02)
  f <- function(x) 1.5 * x * (x[c(2, 3, 1)] - x[c(3, 1, 2)])
  x <- c(0.5, 0.3, 0.2)
  tr <- matrix(0, 301, 3)
  tr[1, ] <- x
  for (k in 1:300) {
    k1 <- f(x)
    k2 <- f(x + 0.01 * k1)
    k3 <- f(x + 0.01 * k2)
    k4 <- f(x + 0.02 * k3)
    x <- x + 0.02 / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    tr[k + 1, ] <- x
  }
  expect_equal(unname(out$traj), tr, tolerance = 1e-13)
  expect_identical(colnames(out$traj), c("x1", "x2", "x3"))
  v <- tr[, 1] * tr[, 2] * tr[, 3]
  expect_equal(out$inv, v, tolerance = 1e-13)
  expect_equal(out$drift, max(abs(v - v[1])), tolerance = 1e-6)
  expect_lt(out$drift, 1e-8)
  expect_equal(out$dvdt, 0, tolerance = 1e-17)
  expect_equal(unname(out$field), f(c(0.5, 0.3, 0.2)), tolerance = 1e-15)
  expect_identical(unname(out$fixed), rep(1 / 3, 3))
})

test_that("insphs() runs the replicator steps on the grid", {
  temp <- c(0.5, 1.2, 1.7)
  cost <- c(0.1, 0.6)
  out <- insphs(temp, cost, steps = 40, rate = 0.3, x0 = c(1, 2, 1))
  one <- function(tt, g) {
    a <- rbind(c(1, 0, 1), c(tt, 0, tt - 1), rep(1 - g, 3))
    x <- c(1, 2, 1) / 4
    for (k in 1:40) {
      f <- drop(a %*% x)
      x <- pmax(0, x + 0.3 * x * (f - sum(x * f)))
      x <- x / sum(x)
    }
    x
  }
  ref <- array(0, c(3, 3, 2))
  for (i in 1:3) for (j in 1:2) ref[, i, j] <- one(temp[i], cost[j])
  expect_equal(unname(out$coop), ref[1, , ], tolerance = 1e-13)
  expect_equal(unname(out$crime), ref[2, , ], tolerance = 1e-13)
  expect_equal(unname(out$insp), ref[3, , ], tolerance = 1e-13)
  expect_identical(dimnames(out$crime), list(c("0.5", "1.2", "1.7"),
                                             c("0.1", "0.6")))
  expect_identical(out$temp, temp)
  expect_identical(dim(insphs(steps = 2)$crime), c(20L, 20L))
})

test_that("colarr() gives the odds ratio among arrestees", {
  a <- c(0.3, 0.6, 0.1)
  e <- 0.2
  bg <- c(0.5, 0.1, 0.01)
  out <- colarr(a, e, bg)
  cl <- cbind(a * e, a * (1 - e), (1 - a) * e, (1 - a) * (1 - e) * bg)
  cl <- cl / rowSums(cl)
  or <- cl[, 1] * cl[, 4] / (cl[, 2] * cl[, 3])
  expect_equal(out$or, or, tolerance = 1e-14)
  expect_equal(out$or, bg, tolerance = 1e-14)
  expect_equal(out$lor, log(or), tolerance = 1e-13)
  expect_equal(as.matrix(out[c("c11", "c10", "c01", "c00")]), cl,
               tolerance = 1e-15, ignore_attr = TRUE)
  expect_identical(out$e, rep(0.2, 3))
})

test_that("the deterrence functions refuse bad input", {
  errs("detrsp", list(
    list("`sanc` must increase with `x`", 1:3, 1:3, c(1, 3, 2), 1),
    list("must not repeat a level", c(1, 1, 2), 1:3, 1:3, 1),
    list("`p` must be at least 0", 1:3, 1:3, 1:3, -1),
    list("`ben` must have length 3", 1:3, 1:2, 1:3, 1)
  ))
  errs("dtrid", list(
    list("`p` and `s` must be positive with `log = TRUE`", c(0, 1), 1:2,
         1:2, log = TRUE),
    list("the design has rank 2 of 4", 1:3, c(2, 2, 2), c(1, 1, 1), 1:3),
    list("`log` must be TRUE or FALSE", 1:2, 1:2, 1:2, log = NA),
    list("`rho` must be one number", 1:2, 1:2, 1:2, rho = NA),
    list("`y` must have length 2", 1:2, 1:2, 1:2, 1)
  ))
  errs("inspg", list(
    list("`a` must not be 0", 0, c(0.5, 0.3, 0.2)),
    list("`a` must be one number", NA, c(0.5, 0.3, 0.2)),
    list("`x` must be three shares summing to 1", 1, c(0.5, 0.5, 0.5)),
    list("`x` must have length 3", 1, c(0.5, 0.5)),
    list("`steps` must be at most 1e7", 1, c(0.5, 0.3, 0.2), steps = 2e7),
    list("`dt` must be one number", 1, c(0.5, 0.3, 0.2), dt = 0)
  ))
  errs("insphs", list(
    list("`x0` must have a positive sum", x0 = c(0, 0, 0)),
    list("`rate` must be one number", rate = 0),
    list("make more than 1e7 cells", temp = seq_len(4000),
         cost = seq_len(4000)),
    list("`steps` must be one whole number", steps = 1.5)
  ))
  errs("colarr", list(
    list("`e` must be one number or one per scenario, 3", c(0.1, 0.2, 0.3),
         c(0.1, 0.2), 0.5),
    list("`a` and `e` must be strictly between 0 and 1", 1, 0.5, 0.5),
    list("`bg` must be finite", 0.5, 0.5, NA)
  ))
})
