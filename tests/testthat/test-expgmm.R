# The exponential-mean GMM (expgmm(), R/apx.R). The reference values were
# computed once with base R alone, by an independent route: the moment
# conditions written out as matrices (X the regressors and the day
# indicator columns, Z the instruments, the control and the same
# columns), each step's criterion minimised by stats::nlminb() with its
# analytic gradient (rel.tol = 1e-15), S summed over the area clusters
# with rowsum(), and the variances (G'S^-1 G)^-1 and the one-step
# sandwich times 40 / 39. The tolerance 1e-6 allows for nlminb()'s
# stopping rule, which is looser than the kernel's Gauss-Newton one. The
# panel is simulated; no number here comes from the papers.

simgm <- function() {
  set.seed(20261116)
  nu <- 40L
  nt <- 10L
  n <- nu * nt
  d <- data.frame(
    unit = rep(seq_len(nu), each = nt),
    day = rep(seq_len(nt), nu)
  )
  d$z1 <- rnorm(n)
  d$z2 <- rnorm(n)
  d$w <- rnorm(n)
  v <- rnorm(n, 0, 0.5)
  d$x <- 0.6 * d$z1 - 0.4 * d$z2 + 0.3 * d$w + v + rnorm(n, 0, 0.5)
  g <- rnorm(nt, 0, 0.3)
  d$y <- rpois(n, exp(0.5 + 0.2 * d$x + 0.1 * d$w + g[d$day] + v))
  d
}

# the reference route, kept here so that the numbers can be regenerated
refgm <- function(d, err) {
  dd <- outer(d$day, 1:10, "==") + 0
  xm <- cbind(d$x, d$w, dd)
  zm <- cbind(d$z1, d$z2, d$w, dd)
  y <- d$y
  u <- function(b) {
    e <- drop(xm %*% b)
    if (err == "mult") y * exp(-e) - 1 else y - exp(e)
  }
  du <- function(b) {
    e <- drop(xm %*% b)
    (if (err == "mult") -y * exp(-e) else -exp(e)) * xm
  }
  fit <- function(w, b0) {
    q <- function(b) {
      g <- crossprod(zm, u(b))
      drop(t(g) %*% w %*% g)
    }
    gr <- function(b) {
      drop(2 * t(crossprod(zm, du(b))) %*% w %*% crossprod(zm, u(b)))
    }
    stats::nlminb(b0, q, gr, control = list(rel.tol = 1e-15,
                                            x.tol = 1e-15,
                                            eval.max = 5000,
                                            iter.max = 5000))$par
  }
  s <- function(b) crossprod(rowsum(zm * u(b), d$unit))
  w1 <- solve(crossprod(zm))
  b1 <- fit(w1, c(0, 0, log(tapply(y, d$day, mean))))
  w2 <- solve(s(b1))
  b2 <- fit(w2, b1)
  g2 <- crossprod(zm, du(b2))
  g1 <- crossprod(zm, du(b1))
  a <- solve(t(g1) %*% w1 %*% g1)
  v1 <- a %*% t(g1) %*% w1 %*% s(b1) %*% w1 %*% g1 %*% a * 40 / 39
  m2 <- crossprod(zm, u(b2))
  list(b1 = b1[1:2], se1 = sqrt(diag(v1))[1:2], b2 = b2[1:2],
       se2 = sqrt(diag(solve(t(g2) %*% w2 %*% g2) * 40 / 39))[1:2],
       j = drop(t(m2) %*% w2 %*% m2))
}

test_that("expgmm equals the moment conditions solved by nlminb()", {
  d <- simgm()
  f2 <- expgmm(d, "y", endo = "x", inst = c("z1", "z2"), x = "w",
               fe = "day", clus = "unit")
  expect_equal(f2$coef$est, c(0.170645698911901, 0.139008628563222),
               tolerance = 1e-6)
  expect_equal(f2$coef$se, c(0.0717823318437703, 0.058829333521988),
               tolerance = 1e-6)
  expect_equal(f2$j$stat, 0.0078947701478867, tolerance = 1e-5)
  expect_identical(f2$j$df, 1)
  expect_true(f2$conv)
  f1 <- expgmm(d, "y", endo = "x", inst = c("z1", "z2"), x = "w",
               fe = "day", clus = "unit", step = "one")
  expect_equal(f1$coef$est, c(0.171074762194083, 0.140454225882242),
               tolerance = 1e-6)
  expect_equal(f1$coef$se, c(0.0720024335929021, 0.0610972019663146),
               tolerance = 1e-6)
  expect_true(is.na(f1$j$stat))
  fa <- expgmm(d, "y", endo = "x", inst = c("z1", "z2"), x = "w",
               fe = "day", clus = "unit", err = "add")
  expect_equal(fa$coef$est, c(0.153736000672564, 0.135616003978279),
               tolerance = 1e-6)
  expect_equal(fa$coef$se, c(0.0627368114734188, 0.0605336690824264),
               tolerance = 1e-6)
  expect_equal(fa$j$stat, 0.578064017487062, tolerance = 1e-5)
  fi <- expgmm(d, "y", endo = "x", inst = c("z1", "z2"), x = "w",
               fe = "day", clus = "unit", step = "iter")
  expect_gte(fi$nw, 2L)
  expect_equal(fi$coef$est, f2$coef$est, tolerance = 1e-2)
  expect_output(print(f2), "Hansen J")
  expect_identical(summary(f2), f2$coef)
})

test_that("expgmm without endogeneity is the Poisson fit", {
  # z = x: the multiplicative moments differ from Poisson's, the additive
  # ones are its score equations, so the additive fit is ppmlfe()'s
  d <- simgm()
  d$xx <- d$x
  fa <- expgmm(d, "y", endo = "x", inst = "xx", x = "w", fe = "day",
               clus = "unit", err = "add")
  fp <- ppmlfe(d, "y", c("x", "w"), fe = "day")
  expect_equal(fa$coef$est, fp$coef$est, tolerance = 1e-7)
})

test_that("expgmm refuses what it cannot fit", {
  d <- simgm()
  expect_error(expgmm(d, "y", endo = "x", inst = character()), "expgmm: ")
  expect_error(expgmm(d, "y", endo = "x", inst = "z1"),
               "needs `clus`")
  expect_error(expgmm(d, "y", endo = "x", inst = c("z1", "z2"),
                      fe = "unit"), "moment conditions but 40 clusters")
  expect_error(expgmm(d, "y", endo = "x", inst = "x", clus = "unit"),
               "both a regressor")
  d$y[d$day == 1] <- 0
  expect_warning(f <- expgmm(d, "y", endo = "x", inst = c("z1", "z2"),
                             fe = "day", clus = "unit"), "dropped")
  expect_identical(f$nsep, 40L)
})
