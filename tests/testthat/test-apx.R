# Exposure and crime (R/apx.R). The reference values were computed once,
# on the simulated panel below, with fixest 0.14.2: fepois() and feols()
# with ssc(K.adj = FALSE, G.adj = TRUE, G.df = "min"), fixef.tol = 3e-12,
# glm.tol = 1e-12; the first-stage statistic is fitstat(, "ivwald1"). The
# control function's reference is feols() for the first stage, then
# fepois() with its residual, and for the bootstrap the same two fits with
# weights equal to the cluster draws below. The panel is simulated; no
# number here comes from the papers.

simpan <- function() {
  set.seed(20261115)
  nu <- 15L
  nt <- 30L
  n <- nu * nt
  d <- data.frame(unit = rep(seq_len(nu), each = nt), day = rep(seq_len(nt), nu))
  a <- rnorm(nu, 0, 0.5)
  g <- rnorm(nt, 0, 0.3)
  d$z1 <- rnorm(n)
  d$z2 <- rnorm(n)
  d$w <- rnorm(n)
  e <- rnorm(n)
  d$x <- 0.8 * d$z1 - 0.5 * d$z2 + 0.3 * d$w + a[d$unit] + 0.5 * e + rnorm(n, 0, 0.5)
  d$pop <- round(runif(nu, 500, 5000))[d$unit]
  d$y <- rpois(n, d$pop / 1000 * exp(0.2 * d$x + 0.1 * d$w + a[d$unit] + g[d$day] + 0.3 * e))
  d$yl <- 1 + 0.5 * d$x + 0.2 * d$w + a[d$unit] + g[d$day] + e + rnorm(n)
  d
}

test_that("ppmlfe equals fepois with two effects, an offset and one- and two-way clusters", {
  d <- simpan()
  f <- ppmlfe(d, "y", c("x", "w"), fe = c("unit", "day"), expo = "pop")
  expect_equal(f$coef$est, c(0.29195662923680077, 0.073073645297524289), tolerance = 1e-8)
  expect_equal(f$coef$se, c(0.02029437190826942, 0.01940571968849104), tolerance = 1e-6)
  expect_true(f$conv)
  expect_identical(f$n, 450L)
  f2 <- ppmlfe(d, "y", c("x", "w"), fe = c("unit", "day"), expo = "pop",
               clus = c("unit", "day"))
  expect_equal(f2$coef$se, c(0.028958692481084434, 0.019570056579545823), tolerance = 1e-6)
  expect_equal(f2$coef$est, f$coef$est, tolerance = 1e-12)
  f1 <- ppmlfe(d, "y", c("x", "w"), fe = "unit", expo = "pop")
  expect_equal(f1$coef$est, c(0.29019645982478903, 0.063439369545109175), tolerance = 1e-8)
  expect_equal(f1$coef$se, c(0.023525811320309342, 0.029059218018322456), tolerance = 1e-6)
  expect_equal(f$coef$rr, exp(f$coef$est))
  expect_identical(summary(f), f$coef)
  expect_output(print(f), "ppmlfe: Poisson pseudo-likelihood, 450 rows")
})

test_that("ppmlfe drops a unit of zeros, as fepois does, and nothing else changes", {
  d <- simpan()
  z <- data.frame(unit = 16L, day = 1:30, z1 = 0, z2 = 0, w = seq(-1, 1, length.out = 30),
                  x = seq(0, 2, length.out = 30), pop = 1000, y = 0, yl = 0)
  d0 <- rbind(d, z)
  expect_warning(f <- ppmlfe(d0, "y", c("x", "w"), fe = c("unit", "day"), expo = "pop"),
                 "30 row\\(s\\) dropped")
  expect_equal(f$coef$est, c(0.29195662923680077, 0.073073645297524289), tolerance = 1e-8)
  expect_identical(f$nsep, 30L)
  expect_identical(f$n, 450L)
})
