# clrilr(), preslr(), tpdmc() and tpbhm(): against log-ratios, normal
# densities, the fundamental matrix and the conjugate posterior computed
# here in base R, and the refusals.

test_that("clrilr() gives the centred and isometric log-ratios", {
  set.seed(6)
  x <- matrix(rgamma(40, 2), 10, 4,
              dimnames = list(NULL, c("a", "b", "c", "d")))
  x[3, 2] <- 0
  ps <- 1e-6
  cl <- (x + ps) / rowSums(x + ps)
  l <- log(cl)
  expect_equal(clrilr(x), l - rowMeans(l), tolerance = 1e-14)
  il <- sapply(1:3, function(i) {
    sqrt(i / (i + 1)) * (rowMeans(l[, 1:i, drop = FALSE]) - l[, i + 1])
  })
  colnames(il) <- paste0("ilr", 1:3)
  expect_equal(clrilr(x, "ilr"), il, tolerance = 1e-13)
  d <- as.data.frame(x[, 1:2])
  expect_identical(colnames(clrilr(d, ps = 0.5)), c("a", "b"))
  expect_null(dimnames(clrilr(unname(x)))[[2]])
})

test_that("preslr() is the ratio of two normal likelihoods", {
  x <- c(1.9, 2.3, 2.1, 1.7)
  s2 <- c(0.4, 0.5, 0.4, 0.6)
  out <- preslr(x, 2, 0.3, 1, s2)
  l1 <- sum(dnorm(x, 2, 0.3, log = TRUE))
  l2 <- sum(dnorm(x, 1, s2, log = TRUE))
  expect_equal(out$ll1, l1, tolerance = 1e-14)
  expect_equal(out$ll2, l2, tolerance = 1e-14)
  expect_equal(out$llr, l1 - l2, tolerance = 1e-14)
  expect_equal(out$lr, exp(l1 - l2), tolerance = 1e-13)
  expect_equal(out$l10, (l1 - l2) / log(10), tolerance = 1e-14)
  band <- c("no support either way", "weak support", "moderate support",
            "moderately strong support", "strong support",
            "very strong support", "extremely strong support")
  b <- findInterval(abs(out$l10), c(0, 1, 2, 3, 4, 6), left.open = TRUE)
  expect_identical(out$verbal, paste(band[b + 1], "for H1"))
  expect_identical(preslr(1, 1, 1, 0, 1)$verbal, "weak support for H1")
  expect_identical(preslr(0, 1, 1, 0, 1)$verbal, "weak support for H2")
  expect_identical(preslr(0, 1, 1, -1, 1)$verbal, "no support either way")
})

test_that("tpdmc() solves the absorbing chain", {
  k <- 4
  pres <- 0.3
  prog <- 0.6 * (1 - pres)
  mum <- 0.4 * pres
  p <- matrix(0, k + 2, k + 2)
  for (i in 1:k) {
    p[i, i] <- 1 - prog - mum
    p[i, i + 1] <- prog
    p[i, k + 2] <- p[i, k + 2] + mum
  }
  p[k + 1, k + 1] <- 1
  p[k + 2, k + 2] <- 1
  n <- solve(diag(k) - p[1:k, 1:k])
  b <- n %*% p[1:k, k + 1:2]
  st <- c("fresh", "bloat", "active", "advanced")
  out <- tpdmc(pres, 0.6, 0.4)
  expect_equal(unname(out$P), p, tolerance = 1e-15)
  expect_identical(dimnames(out$P)[[1]], c(st, "skeletal", "mummified"))
  expect_equal(unname(out$N), n, tolerance = 1e-14)
  expect_equal(unname(out$absorb), unname(b), tolerance = 1e-14)
  expect_identical(colnames(out$absorb), c("skeletal", "mummified"))
  expect_equal(unname(out$time), rowSums(n), tolerance = 1e-14)
  expect_identical(names(out$time), st)
  expect_null(out$path)
  # no preservation: nothing mummifies
  expect_identical(unname(tpdmc(0)$absorb[, "mummified"]), rep(0, 4))
  # a drawn path starts where asked and ends absorbed, or after steps
  s <- tpdmc(pres, 0.6, 0.4, start = "bloat", steps = 500, seed = 2)
  expect_identical(s$start, 2L)
  expect_identical(s$path[1], "bloat")
  expect_true(s$path[length(s$path)] %in% c("skeletal", "mummified"))
  expect_identical(s, tpdmc(pres, 0.6, 0.4, start = 2, steps = 500, seed = 2))
  short <- tpdmc(states = c("a", "b"), steps = 1, seed = 3)$path
  expect_identical(length(short), 2L)
  set.seed(1)
  expect_length(tpdmc(0.5, steps = 3)$path[1], 1)
})

test_that("tpbhm() is the conjugate normal posterior", {
  set.seed(7)
  n <- 40
  x <- cbind(lime = rbinom(n, 1, 0.5), temp = rnorm(n, 12, 3))
  y <- 0.3 + 0.2 * x[, 1] - 0.01 * x[, 2] + rnorm(n, sd = 0.05)
  g <- rep(c("a", "b", "c", "d"), 10)
  m0 <- c(0, 0.2, 0)
  s0 <- c(10, 0.1, 1)
  out <- tpbhm(y, x, g, pm = m0, psd = s0, level = 0.9)
  xm <- cbind(1, x)
  s2 <- sum(lm.fit(xm, y)$residuals^2) / (n - 3)
  sig <- solve(crossprod(xm) / s2 + diag(1 / s0^2))
  mu <- as.vector(sig %*% (crossprod(xm, y) / s2 + m0 / s0^2))
  sd <- unname(sqrt(diag(sig)))
  expect_identical(out$coef$term, c("(Intercept)", "lime", "temp"))
  expect_equal(out$coef$mean, mu, tolerance = 1e-10)
  expect_equal(out$coef$sd, sd, tolerance = 1e-10)
  expect_equal(out$coef$lo, mu - qnorm(0.95) * sd, tolerance = 1e-10)
  expect_equal(out$coef$hi, mu + qnorm(0.95) * sd, tolerance = 1e-10)
  expect_equal(out$coef$ppos, pnorm(mu / sd), tolerance = 1e-10)
  expect_identical(out$coef$prior, s0)
  expect_equal(out$sigma, sqrt(s2), tolerance = 1e-12)
  f <- as.vector(xm %*% mu)
  expect_equal(out$fitted, f, tolerance = 1e-12)
  gm <- c(tapply(y - f, g, mean))
  nj <- c(tapply(y - f, g, length))
  tau2 <- max(0, var(gm) - mean(s2 / nj))
  lam <- tau2 / (tau2 + s2 / nj)
  expect_equal(out$tau2, tau2, tolerance = 1e-10)
  expect_identical(out$groups$group, names(gm))
  expect_identical(out$groups$n, as.integer(nj))
  expect_equal(out$groups$raw, unname(gm), tolerance = 1e-10)
  expect_equal(out$groups$shrink, unname(lam), tolerance = 1e-8)
  expect_equal(out$groups$pooled, unname(lam * gm), tolerance = 1e-8)
  # a data frame, no groups, one prior for all
  d <- tpbhm(y, as.data.frame(x))
  expect_null(d$groups)
  expect_true(is.na(d$tau2))
  expect_identical(d$coef$prior, rep(10, 3))
  expect_identical(tpbhm(y, unname(x))$coef$term, c("(Intercept)", "x1", "x2"))
})

test_that("the taphonomy functions refuse bad input", {
  errs("clrilr", list(
    list("at least two parts", matrix(1:3)),
    list("`type` must be one of", matrix(1:4, 2), "x"),
    list("`ps` must be one number from 0", matrix(1:4, 2), ps = -1),
    list("a part is 0 after the pseudocount", matrix(0:3, 2), ps = 0),
    list("`x` must hold finite numbers from 0", matrix(-1:2, 2))
  ))
  errs("preslr", list(
    list("`s1` must be one number or one per value, 2", 1:2, 0, 1:3, 0, 1),
    list("`s2` must be greater than 0", 1, 0, 1, 0, 0),
    list("`m1` must be finite", 1, NA, 1, 0, 1)
  ))
  errs("tpdmc", list(
    list("`pres` must be one number from 0 to 1", 2),
    list("`decay` must be above 0", decay = 0),
    list("`decay` must be one number from 0 to 1", decay = 2),
    list("`states` must name the transient stages", states = 1:2),
    list("`states` must name each stage once", states = c("a", "a")),
    list("`start` must be a transient stage", start = "skeletal"),
    list("`start` must be a transient stage", start = 5),
    list("`start` must be a transient stage", start = TRUE),
    list("`steps` must be one whole number from 0", steps = -1),
    list("nothing leaves a stage", 1, mum = 0)
  ))
  x <- cbind(1:6, c(2, 1, 4, 3, 6, 5))
  errs("tpbhm", list(
    list("`x` must have one row per case, 3", 1:3, x),
    list("at least two cases", 1, 1),
    list("`pm` and `psd` must be one number or one per term, 3", 1:6, x,
         pm = c(1, 2)),
    list("`grp` must hold at least two groups", 1:6, x, grp = rep(1, 6)),
    list("`x` has dependent columns", 1:6, cbind(1:6, 2 * (1:6))),
    list("`level` must be one number", 1:6, x, level = 0)
  ))
})
