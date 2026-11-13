# wcbx() and wystep(): the wild cluster bootstrap against every sign
# pattern computed here with lm(), the permutations against every
# assignment refitted with didfe(), the power's size against its level,
# the stepdown against algorithm 2.8 written out here, and the refusals.

wpanel <- function(seed, nu = 8) {
  set.seed(seed)
  d <- expand.grid(year = 2015:2024, town = sprintf("t%02d", seq_len(nu)))
  i <- as.integer(d$town)
  d$on <- as.integer(i <= 3 & d$year >= 2019)
  d$x <- rnorm(10 * nu)
  d$rate <- 30 + rnorm(nu)[i] + 0.5 * (d$year - 2015) - 1.5 * d$on +
    0.3 * d$x + rnorm(10 * nu)
  d
}

# The restricted wild bootstrap p-value of b0 over every sign pattern, by
# lm() on the dummy regression: K counts the year effects, the regressors
# and the constant (the units are nested in the clusters).
wcr <- function(d, b0, adj = FALSE) {
  f <- if (adj) rate ~ on + x + factor(town) + factor(year) else
    rate ~ on + factor(town) + factor(year)
  x <- model.matrix(f, d)
  g <- as.integer(d$town)
  ng <- max(g)
  k <- 1 + 9 + 1 + adj
  tt <- function(y) {
    m <- lm.fit(x, y)
    b <- solve(crossprod(x))
    s <- rowsum(x * m$residuals, g)
    v <- (b %*% crossprod(s) %*% b)[2, 2] * ng / (ng - 1) * 79 / (80 - k)
    (m$coefficients[2] - b0) / sqrt(v)
  }
  y0 <- d$rate - b0 * d$on
  r <- lm.fit(x[, -2], y0)
  t0 <- abs(tt(d$rate))
  s <- as.matrix(expand.grid(rep(list(c(-1, 1)), ng)))
  ts <- apply(s, 1, function(v) {
    abs(tt(r$fitted.values + b0 * d$on + r$residuals * v[g]))
  })
  c(t0, mean(ts >= t0 - 1e-12), mean(ts > t0 + 1e-12))
}

test_that("the wild bootstrap enumerates every sign pattern", {
  d <- wpanel(1)
  r <- wcbx(d$rate, d$town, d$year, d$on)
  w <- wcr(d, 0)
  expect_true(r$enum)
  expect_identical(r$draws, 256L)
  expect_equal(abs(r$res$stat), unname(w[1]), tolerance = 1e-10)
  expect_equal(r$res$p, unname(w[2]), tolerance = 1e-14)
  f <- didfe(d$rate, d$town, d$year, d$on)
  expect_equal(r$res$stat, f$coef$stat, tolerance = 1e-10)
  rs <- wcbx(d$rate, d$town, d$year, d$on, ties = "strict")
  expect_equal(rs$res$p, unname(w[3]), tolerance = 1e-14)
  # the interval's ends are where the p-value crosses 5%
  for (e in c(r$res$lo, r$res$hi)) {
    expect_lte(wcr(d, e + 1e-6 * sign(e - r$res$est))[2], 0.05)
    expect_gt(wcr(d, e - 1e-6 * sign(e - r$res$est))[2], 0.05)
  }
  ra <- wcbx(d$rate, d$town, d$year, d$on, adj = d["x"])
  expect_equal(ra$res$p, unname(wcr(d, 0, adj = TRUE)[2]),
               tolerance = 1e-14)
  expect_output(print(r), "all enumerated")
})

test_that("random draws are seeded; few clusters give the whole line", {
  d <- wpanel(2)
  a <- wcbx(d$rate, d$town, d$year, d$on, wtype = "webb", reps = 99,
            seed = 1)
  b <- wcbx(d$rate, d$town, d$year, d$on, wtype = "webb", reps = 99,
            seed = 1)
  expect_identical(a, b)
  expect_false(a$enum)
  set.seed(3)
  c1 <- wcbx(d$rate, d$town, d$year, d$on, reps = 99)
  set.seed(3)
  c2 <- wcbx(d$rate, d$town, d$year, d$on, reps = 99)
  expect_identical(c1, c2)
  expect_false(c1$enum)
  # four clusters: the smallest p-value is 2/16 > 5%
  reg <- (as.integer(d$town) + 1) %/% 2
  r <- wcbx(d$rate, d$town, d$year, d$on, clus = reg)
  expect_identical(c(r$res$lo, r$res$hi), c(-Inf, Inf))
  # an effect far beyond any null: the interval is bounded
  d$rate <- d$rate - 50 * d$on
  r2 <- wcbx(d$rate, d$town, d$year, d$on, wtype = "webb", reps = 199,
             seed = 2)
  expect_true(all(is.finite(c(r2$res$lo, r2$res$hi))))
})

test_that("the permutations are every assignment of the dates", {
  d <- wpanel(4)
  r <- wcbx(d$rate, d$town, d$year, d$on, meth = "perm")
  # choose(8, 3) = 56 assignments of the 2019 date
  cmb <- combn(8, 3)
  ts <- apply(cmb, 2, function(s) {
    tr <- as.integer(as.integer(d$town) %in% s & d$year >= 2019)
    didfe(d$rate, d$town, d$year, tr)$coef$stat
  })
  expect_identical(r$draws, 56L)
  expect_true(r$enum)
  expect_equal(r$res$p, mean(abs(ts) >= abs(r$res$stat) - 1e-12),
               tolerance = 1e-14)
  re <- wcbx(d$rate, d$town, d$year, d$on, meth = "perm", stat = "est",
             ties = "strict")
  es <- apply(cmb, 2, function(s) {
    tr <- as.integer(as.integer(d$town) %in% s & d$year >= 2019)
    didfe(d$rate, d$town, d$year, tr)$coef$est
  })
  expect_equal(re$res$p, mean(abs(es) > abs(re$res$est) + 1e-12),
               tolerance = 1e-14)
  # two outcomes: the stepdown over the same assignments
  y2 <- data.frame(a = d$rate, b = d$x)
  r2 <- wcbx(y2, d$town, d$year, d$on, meth = "perm")
  expect_identical(r2$res$outcome, c("a", "b"))
  tb <- apply(cmb, 2, function(s) {
    tr <- as.integer(as.integer(d$town) %in% s & d$year >= 2019)
    didfe(d$x, d$town, d$year, tr)$coef$stat
  })
  w <- wystep(r2$res$stat, cbind(ts, tb))
  expect_equal(r2$res$pstep, w$padj, tolerance = 1e-14)
  r3 <- wcbx(cbind(d$rate, d$x), d$town, d$year, d$on, meth = "perm")
  expect_identical(r3$res$outcome, c("y1", "y2"))
  # random assignments when there are more than reps
  r4 <- wcbx(d$rate, d$town, d$year, d$on, meth = "perm", reps = 20,
             seed = 1)
  expect_false(r4$enum)
  expect_equal(r4$res$p * 21, round(r4$res$p * 21), tolerance = 1e-12)
})

test_that("the power's size is near its level; the effect is detected", {
  d <- wpanel(5, nu = 30)
  d$rate <- d$rate + 1.5 * d$on
  r <- wcbx(d$rate, d$town, d$year, meth = "power", ntrt = 10,
            adopt = 2018:2021, eff = c(0, 0.5, 3), reps = 400, seed = 1)
  expect_lt(abs(r$size - 0.05), 4 * sqrt(0.05 * 0.95 / 400))
  expect_identical(r$res$power[3], 1)
  p <- r$res$power
  want <- 0.5 + (0.8 - p[2]) / (p[3] - p[2]) * 2.5
  expect_equal(r$mde, if (p[2] >= 0.8) 0.5 else want,
               tolerance = testthat_tolerance())
  r0 <- wcbx(d$rate, d$town, d$year, meth = "power", ntrt = 10,
             adopt = 2020, reps = 20, seed = 2, power = 0.01)
  expect_identical(nrow(r0$res), 5L)
  expect_identical(r0$mde, 0)
  expect_output(print(r0), "minimum detectable effect")
})

test_that("the stepdown is algorithm 2.8", {
  set.seed(6)
  null <- matrix(rnorm(3000), ncol = 3)
  null[, 2] <- null[, 2] + 0.5 * null[, 1]
  obs <- c(2.5, 0.4, -2.1)
  r <- wystep(obs, null)
  o <- order(-abs(obs))
  q <- abs(null[, o])
  q[, 2] <- pmax(q[, 2], q[, 3])
  q[, 1] <- pmax(q[, 1], q[, 2])
  p <- cummax(colMeans(sweep(q, 2, abs(obs[o]) - 1e-12, ">=")))
  expect_equal(r$padj[o], p, tolerance = 1e-14)
  expect_equal(r$p, colMeans(abs(null) >= rep(abs(obs), each = 1000)),
               tolerance = 1e-14)
  nn <- rbind(null, NA)
  expect_identical(wystep(obs, nn), r)
  expect_identical(wystep(1, matrix(c(0.5, 2), 2))$padj, 0.5)
})

test_that("wcbx() and wystep() refuse what they cannot test", {
  d <- wpanel(7)
  y <- d$rate
  u <- d$town
  t <- d$year
  expect_error(wcbx(y, u, t), "`trt` is needed")
  expect_error(wcbx(cbind(y, y), u, t, d$on), "several outcomes")
  expect_error(wcbx(list(), u, t, d$on), "no outcome")
  expect_error(wcbx(y, u, t, d$on, clus = rep(1, 80)), "one cluster")
  expect_error(wcbx(y, u, t, d$on, adj = data.frame(a = as.integer(u))),
               "absorb")
  expect_error(wcbx(y, u, t, d$on, meth = "perm",
                    clus = (as.integer(u) + 1) %/% 2), "same period")
  expect_error(wcbx(y, u, t, d$on, meth = "perm",
                    adj = data.frame(a = as.integer(u))), "absorb")
  expect_error(wcbx(y, u, t, meth = "power", ntrt = 8), "untreated")
  expect_error(wcbx(y, u, t, meth = "power"), "adoption periods")
  expect_error(wcbx(y, u, t, meth = "power", adopt = 2015), "after the first")
  expect_error(wcbx(y, u, t, meth = "power", adopt = numeric()), "adopt")
  expect_error(wcbx(y, u, t, meth = "power", adopt = 2019,
                    eff = numeric()), "eff")
  expect_error(wcbx(y, u, t, d$on, reps = 3), "reps")
  expect_error(wystep(c(1, NA), matrix(0, 2, 2)), "obs")
  expect_error(wystep(1:2, matrix(0, 2, 3)), "one column per")
  expect_error(wystep(1, matrix(NA_real_, 2, 1)), "no draw")
})
