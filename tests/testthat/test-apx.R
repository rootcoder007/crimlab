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

test_that("ivfe equals feols IV: coefficients, clustered errors, first-stage Wald F", {
  d <- simpan()
  f <- ivfe(d, "yl", endo = "x", inst = c("z1", "z2"), x = "w", fe = c("unit", "day"))
  expect_equal(f$coef$est, c(0.53281354358983513, 0.21211428905136548), tolerance = 1e-8)
  expect_equal(f$coef$se, c(0.04844890651382984, 0.043953255267726352), tolerance = 1e-6)
  expect_equal(f$first$fstat, 420.99264832346375, tolerance = 1e-6)
  f2 <- ivfe(d, "yl", endo = "x", inst = c("z1", "z2"), x = "w", fe = c("unit", "day"),
             clus = c("unit", "day"))
  expect_equal(f2$coef$se, c(0.070034681564906318, 0.062053418838878639), tolerance = 1e-6)
  fj <- ivfe(d, "yl", endo = "x", inst = "z1", x = "w", fe = c("unit", "day"))
  expect_equal(fj$coef$est, c(0.42225638810732974, 0.2386127220186717), tolerance = 1e-8)
  expect_equal(fj$coef$se, c(0.072590551828609379, 0.044684989103406354), tolerance = 1e-6)
  expect_output(print(f), "two-stage least squares")
})

test_that("ivfe without instruments is least squares with the effects (didfe's fit)", {
  d <- simpan()
  d$trt <- as.integer(d$unit <= 7 & d$day > 15)
  a <- ivfe(d, "yl", x = "trt", fe = c("unit", "day"))
  b <- didfe(d$yl, d$unit, d$day, d$trt)
  expect_equal(a$coef$est, b$coef$est[1L], tolerance = 1e-10)
})

test_that("ppmlcf equals the two fixest fits and their cluster bootstrap", {
  d <- simpan()
  f <- ppmlcf(d, "y", endo = "x", inst = c("z1", "z2"), x = "w", fe = c("unit", "day"),
              expo = "pop", nboot = 30, seed = 7)
  expect_equal(f$coef$est, c(0.21298786456069305, 0.09126882860013609, 0.26288312058292385),
               tolerance = 1e-8)
  expect_equal(unname(f$draws[1L, ]),
               c(0.23993633897694167, 0.08229024562148296, 0.2806087297677553), tolerance = 1e-8)
  expect_equal(f$coef$se, c(0.03202807019148101, 0.031920108481147909, 0.051975977465311039),
               tolerance = 1e-6)
  expect_identical(f$coef$term, c("x", "w", "vhat:x"))
  expect_identical(f$nfail, 0L)
  # the seed leaves the generator as it was
  set.seed(1)
  u1 <- runif(1)
  set.seed(1)
  ppmlcf(d, "y", endo = "x", inst = "z1", fe = c("unit", "day"), nboot = 2, seed = 3)
  expect_identical(runif(1), u1)
})

test_that("the refusals: no variation, too few clusters, bad columns", {
  d <- simpan()
  d$zc <- 1
  expect_error(ivfe(d, "yl", endo = "x", inst = "zc", fe = c("unit", "day")),
               "does not vary within the fixed effects")
  d$uz <- d$unit
  expect_error(ppmlfe(d, "y", "uz", fe = c("unit", "day")),
               "does not vary within the fixed effects")
  d$one <- 1
  expect_error(ppmlfe(d, "y", "x", fe = c("unit", "day"), clus = "one"), "one cluster")
  expect_warning(ppmlfe(d[d$unit <= 5, ], "y", "x", fe = c("unit", "day")),
                 "fewer than ten clusters")
  expect_error(ppmlfe(d, "y", "nope", fe = "unit"), "not in `data`")
  expect_error(ivfe(d, "yl", endo = c("x", "w"), inst = "z1", fe = "unit"),
               "at least as many instruments")
  d$y[1] <- -1
  expect_error(ppmlfe(d, "y", "x", fe = "unit"), "at least 0")
})

test_that("windiv cuts sectors, leaves one out and interacts with areas", {
  m <- windiv(c(10, 100, 200, 300, 360))
  expect_identical(colnames(m), c("w90_180", "w180_270", "w270_0"))
  expect_identical(unname(m[2, ]), c(1, 0, 0))
  expect_identical(unname(m[5, ]), c(0, 0, 0))
  m2 <- windiv(c(10, 100), area = c("b", "a"))
  expect_identical(colnames(m2)[c(1, 4)], c("w90_180:a", "w90_180:b"))
  expect_identical(unname(m2[2, 1]), 1)
  m3 <- windiv(c(230, 300, 20), sect = c(225, 285))
  expect_identical(unname(m3[, 1]), c(1, 0, 0))
  expect_true(all(is.na(windiv(NA_real_)[1, ])))
  expect_error(windiv(1, width = 70), "must divide 360")
})

test_that("invsn is the temperature aloft less the ground", {
  r <- invsn(c(10, 5, 2), c(8, 6, 2))
  expect_identical(r$diff, c(-2, 1, 0))
  expect_identical(r$inv, c(0L, 1L, 0L))
})

test_that("plcday refits on shifted days and the same day equals the plain fit", {
  d <- simpan()
  p <- plcday(d, "x", "unit", "day", shifts = c(-7, 7), y = "y", x = c("x", "w"),
              fe = c("unit", "day"), expo = "pop")
  f <- ppmlfe(d, "y", c("x", "w"), fe = c("unit", "day"), expo = "pop")
  expect_identical(p$tab$shift, c(0, -7, 7))
  expect_equal(p$tab$est[1L], f$coef$est[1L])
  expect_identical(p$tab$n, c(450L, 345L, 345L))
  # by hand: the exposure of a week later
  d7 <- d[d$day <= 23, ]
  d7$x <- d$x[match(paste(d7$unit, d7$day + 7), paste(d$unit, d$day))]
  f7 <- ppmlfe(d7, "y", c("x", "w"), fe = c("unit", "day"), expo = "pop")
  expect_equal(p$tab$est[3L], f7$coef$est[1L])
  d2 <- rbind(d, d[1, ])
  expect_error(plcday(d2, "x", "unit", "day", y = "y", x = "x", fe = "unit"),
               "two rows on one day")
})

test_that("updwn classifies sides and estimates the downwind difference", {
  set.seed(5)
  a <- data.frame(area = 1:20, ax = runif(20, 0, 5000), ay = rep(c(-800, 800), 10))
  d <- merge(a, data.frame(day = 1:150))
  wnd <- runif(150, 0, 360)
  d$wind <- wnd[d$day]
  blown <- cos((d$wind + 180) * pi / 180)
  dn <- (d$ay > 0 & blown > 0.5) | (d$ay < 0 & blown < -0.5)
  d$crimes <- rpois(nrow(d), 3 * (1 + 0.2 * dn))
  road <- data.frame(x0 = 0, y0 = 0, x1 = 5000, y1 = 0)
  r <- updwn(d, "crimes", "ax", "ay", "day", "wind", road, band = 1600)
  # 60 degrees of the north-south line: the days with |cos| >= 0.5
  ndays <- sum(abs(cos(wnd * pi / 180)) >= 0.5)
  expect_identical(r$n, 2L * ndays)
  expect_identical(unname(r$class[["along"]]), 20L * (150L - ndays))
  # by hand: the sides' daily sums, side and day effects, least squares
  k <- abs(cos(d$wind * pi / 180)) >= 0.5
  s <- aggregate(crimes ~ day + north, data.frame(d[k, ], north = d$ay[k] > 0), sum)
  s$down <- as.numeric((s$north & cos((wnd[s$day] + 180) * pi / 180) > 0) |
                         (!s$north & cos((wnd[s$day] + 180) * pi / 180) < 0))
  s$crimes <- s$crimes / mean(s$crimes)
  s$cl <- seq_len(nrow(s))
  h <- ivfe(s, "crimes", x = "down", fe = c("north", "day"), clus = "cl")
  expect_equal(r$coef$est, h$coef$est, tolerance = 1e-10)
  rp <- updwn(d, "crimes", "ax", "ay", "day", "wind", road, band = 1600,
              offs = c(-1600, 0, 1600))
  expect_identical(rp$placebo$off, c(-1600, 0, 1600))
  expect_equal(rp$placebo$est[2L], r$coef$est)
  expect_error(updwn(d, "crimes", "ax", "ay", "day", "wind", road, tol = 100),
               "between 0 and 90")
})
