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
  d <- data.frame(
    unit = rep(seq_len(nu), each = nt),
    day = rep(seq_len(nt), nu)
  )
  a <- rnorm(nu, 0, 0.5)
  g <- rnorm(nt, 0, 0.3)
  d$z1 <- rnorm(n)
  d$z2 <- rnorm(n)
  d$w <- rnorm(n)
  e <- rnorm(n)
  d$x <- 0.8 * d$z1 - 0.5 * d$z2 + 0.3 * d$w + a[d$unit] + 0.5 * e +
    rnorm(n, 0, 0.5)
  d$pop <- round(runif(nu, 500, 5000))[d$unit]
  d$y <- rpois(
    n,
    d$pop / 1000 * exp(0.2 * d$x + 0.1 * d$w + a[d$unit] + g[d$day] + 0.3 * e)
  )
  d$yl <- 1 + 0.5 * d$x + 0.2 * d$w + a[d$unit] + g[d$day] + e + rnorm(n)
  d
}

test_that("ppmlfe equals fepois: two effects, an offset, two clusterings", {
  d <- simpan()
  f <- ppmlfe(d, "y", c("x", "w"), fe = c("unit", "day"), expo = "pop")
  expect_equal(
    f$coef$est,
    c(0.29195662923680077, 0.073073645297524289),
    tolerance = 1e-8
  )
  expect_equal(
    f$coef$se,
    c(0.02029437190826942, 0.01940571968849104),
    tolerance = 1e-6
  )
  expect_true(f$conv)
  expect_identical(f$n, 450L)
  f2 <- ppmlfe(d, "y", c("x", "w"), fe = c("unit", "day"), expo = "pop",
               clus = c("unit", "day"))
  expect_equal(
    f2$coef$se,
    c(0.028958692481084434, 0.019570056579545823),
    tolerance = 1e-6
  )
  expect_equal(f2$coef$est, f$coef$est, tolerance = 1e-12)
  f1 <- ppmlfe(d, "y", c("x", "w"), fe = "unit", expo = "pop")
  expect_equal(
    f1$coef$est,
    c(0.29019645982478903, 0.063439369545109175),
    tolerance = 1e-8
  )
  expect_equal(
    f1$coef$se,
    c(0.023525811320309342, 0.029059218018322456),
    tolerance = 1e-6
  )
  expect_equal(f$coef$rr, exp(f$coef$est))
  expect_identical(summary(f), f$coef)
  expect_output(print(f), "ppmlfe: Poisson pseudo-likelihood, 450 rows")
})

test_that("ppmlfe drops a unit of zeros, as fepois does", {
  d <- simpan()
  z <- data.frame(
    unit = 16L,
    day = 1:30,
    z1 = 0,
    z2 = 0,
    w = seq(-1, 1, length.out = 30),
    x = seq(0, 2, length.out = 30),
    pop = 1000,
    y = 0,
    yl = 0
  )
  d0 <- rbind(d, z)
  expect_warning(
    f <- ppmlfe(d0, "y", c("x", "w"), fe = c("unit", "day"), expo = "pop"),
    "30 row\\(s\\) dropped"
  )
  expect_equal(
    f$coef$est,
    c(0.29195662923680077, 0.073073645297524289),
    tolerance = 1e-8
  )
  expect_identical(f$nsep, 30L)
  expect_identical(f$n, 450L)
})

test_that("ivfe equals feols IV: estimates, errors, first-stage F", {
  d <- simpan()
  f <- ivfe(
    d,
    "yl",
    endo = "x",
    inst = c("z1", "z2"),
    x = "w",
    fe = c("unit", "day")
  )
  expect_equal(
    f$coef$est,
    c(0.53281354358983513, 0.21211428905136548),
    tolerance = 1e-8
  )
  expect_equal(
    f$coef$se,
    c(0.04844890651382984, 0.043953255267726352),
    tolerance = 1e-6
  )
  expect_equal(f$first$wald, 420.99264832346375, tolerance = 1e-6)
  # the effective F of Montiel Olea and Pflueger, by hand from fixest's
  # first stage (its clustered covariance and the instruments' cross
  # products after the control and the effects)
  expect_equal(f$first$eff, 536.97395780637555, tolerance = 1e-6)
  expect_identical(f$first$kp, f$first$wald)
  # the effective degrees of freedom of the simplified test (Montiel Olea
  # and Pflueger 2013, equation 4, x = 10), by base R matrix algebra on
  # the same Q and clustered V; the critical value is their Patnaik
  # formula, 23.11 for keff = 1 in their table 1
  expect_equal(f$first$keff, 1.51825656086039, tolerance = 1e-6)
  expect_equal(f$first$cv, stats::qchisq(0.95, f$first$keff,
                                         ncp = 10 * f$first$keff) /
                 f$first$keff, tolerance = 1e-12)
  # the first stage in units (Yamada 2026, item 5): fixest's first-stage
  # coefficients over the mean of x
  expect_equal(f$fstab$est, c(0.73993728634395717, -0.5556931692314554),
               tolerance = 1e-8)
  expect_equal(f$fstab$rel, f$fstab$est / mean(d$x))
  f2 <- ivfe(
    d,
    "yl",
    endo = "x",
    inst = c("z1", "z2"),
    x = "w",
    fe = c("unit", "day"),
    clus = c("unit", "day")
  )
  expect_equal(
    f2$coef$se,
    c(0.070034681564906318, 0.062053418838878639),
    tolerance = 1e-6
  )
  fj <- ivfe(d, "yl", endo = "x", inst = "z1", x = "w", fe = c("unit", "day"))
  expect_equal(
    fj$coef$est,
    c(0.42225638810732974, 0.2386127220186717),
    tolerance = 1e-8
  )
  expect_equal(
    fj$coef$se,
    c(0.072590551828609379, 0.044684989103406354),
    tolerance = 1e-6
  )
  expect_output(print(f), "two-stage least squares")
})

test_that("ivfe without instruments is didfe's least squares", {
  d <- simpan()
  d$trt <- as.integer(d$unit <= 7 & d$day > 15)
  a <- ivfe(d, "yl", x = "trt", fe = c("unit", "day"))
  b <- didfe(d$yl, d$unit, d$day, d$trt)
  expect_equal(a$coef$est, b$coef$est[1L], tolerance = 1e-10)
})

test_that("ppmlcf equals the two fixest fits and their cluster bootstrap", {
  d <- simpan()
  f <- ppmlcf(
    d,
    "y",
    endo = "x",
    inst = c("z1", "z2"),
    x = "w",
    fe = c("unit", "day"),
    expo = "pop",
    nboot = 30,
    seed = 7
  )
  expect_equal(
    f$coef$est,
    c(0.21298786456069305, 0.09126882860013609, 0.26288312058292385),
    tolerance = 1e-8
  )
  expect_equal(
    unname(f$draws[1L, ]),
    c(0.23993633897694167, 0.08229024562148296, 0.2806087297677553),
    tolerance = 1e-8
  )
  expect_equal(
    f$coef$se,
    c(0.03202807019148101, 0.031920108481147909, 0.051975977465311039),
    tolerance = 1e-6
  )
  expect_identical(f$coef$term, c("x", "w", "vhat:x"))
  expect_identical(f$nfail, 0L)
  # the seed leaves the generator as it was
  set.seed(1)
  u1 <- runif(1)
  set.seed(1)
  ppmlcf(
    d,
    "y",
    endo = "x",
    inst = "z1",
    fe = c("unit", "day"),
    nboot = 2,
    seed = 3
  )
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
  expect_error(
    ppmlfe(d, "y", "x", fe = c("unit", "day"), clus = "one"),
    "one cluster"
  )
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

test_that("plcday refits on shifted days; shift 0 is the plain fit", {
  d <- simpan()
  p <- plcday(
    d,
    "x",
    "unit",
    "day",
    shifts = c(-7, 7),
    y = "y",
    x = c("x", "w"),
    fe = c("unit", "day"),
    expo = "pop"
  )
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

# Simulated areas about a road along the x axis from 0 to 5000: 20 of
# them 800 north or south of it, 20 more spread up to 4000 away (so
# that the placebo roads moved sideways have areas on both sides).
simud <- function() {
  set.seed(5)
  a <- data.frame(area = 1:40, ax = runif(40, 0, 5000),
                  ay = c(rep(c(-800, 800), 10), runif(20, -4000, 4000)))
  d <- merge(a, data.frame(day = 1:150))
  wnd <- runif(150, 0, 360)
  d$wind <- wnd[d$day]
  d$temp <- rnorm(150, 20, 5)[d$day]
  blown <- cos((d$wind + 180) * pi / 180)
  dn <- (d$ay > 0 & blown > 0.5) | (d$ay < 0 & blown < -0.5)
  d$crimes <- rpois(nrow(d), 3 * (1 + 0.2 * dn))
  list(d = d, wnd = wnd)
}

# The sides' daily sums by hand for the areas within band of the road.
udhand <- function(d, wnd, band = 1600) {
  k <- abs(cos(d$wind * pi / 180)) >= 0.5 & abs(d$ay) <= band
  s <- aggregate(crimes ~ day + north,
                 data.frame(d[k, ], north = d$ay[k] > 0), sum)
  bl <- cos((wnd[s$day] + 180) * pi / 180)
  s$down <- as.numeric((s$north & bl > 0) | (!s$north & bl < 0))
  s$cl <- seq_len(nrow(s))
  s
}

test_that("updwn classifies sides and estimates the downwind difference", {
  u <- simud()
  d <- u$d
  road <- data.frame(x0 = 0, y0 = 0, x1 = 5000, y1 = 0)
  r <- updwn(d, "crimes", "ax", "ay", "day", "wind", road, band = 1600)
  # 60 degrees of the north-south line: the days with |cos| >= 0.5
  ndays <- sum(abs(cos(u$wnd * pi / 180)) >= 0.5)
  inb <- sum(abs(d$ay[d$day == 1]) <= 1600)
  expect_identical(r$n, 2L * ndays)
  expect_identical(unname(r$class[["along"]]), inb * (150L - ndays))
  s <- udhand(d, u$wnd)
  s$crimes <- s$crimes / mean(s$crimes)
  h <- ivfe(s, "crimes", x = "down", fe = c("north", "day"), clus = "cl")
  expect_equal(r$coef$est, h$coef$est, tolerance = 1e-10)
  rp <- updwn(d, "crimes", "ax", "ay", "day", "wind", road, band = 1600,
              offs = c(-1600, 0, 1600))
  expect_identical(rp$placebo$off, c(-1600, 0, 1600))
  expect_equal(rp$placebo$est[2L], r$coef$est)
  expect_true(all(is.finite(rp$placebo$est)))
  expect_error(updwn(d, "crimes", "ax", "ay", "day", "wind", road, tol = 100),
               "between 0 and 90")
})

test_that("expjn joins, counts misses and refuses duplicate keys", {
  ex <- data.frame(area = c("a", "a", "b"), day = c(1, 2, 1), pm = c(10, 12, 8))
  oc <- data.frame(
    area = c("a", "a", "b", "b"),
    day = c(1, 2, 1, 2),
    crimes = c(3, 1, 0, 2)
  )
  j <- expjn(ex, oc)
  expect_identical(j$data$pm, c(10, 12, 8))
  expect_identical(j$nout, 1L)
  expect_identical(j$nexp, 0L)
  expect_identical(j$unout$area, "b")
  expect_error(expjn(rbind(ex, ex[1, ]), oc), "`expo` has 1 duplicate key")
  expect_error(expjn(ex, rbind(oc, oc[2, ])), "`outc` has 1 duplicate key")
  expect_output(print(j), "1 outcome row\\(s\\) without an exposure")
})

# Second round: the options. References from the same fixest session
# (fepois with ssc() at its default for "nk"; vcov = "hetero" with
# ssc(K.adj = FALSE, G.adj = FALSE) for "hc1" with "none"; NW(2) ~ unit +
# day with no adjustment, which equals the panel Newey-West computed by
# hand; the plug-in as fepois on x less feols' residual, and the cluster
# bootstraps with the draws of set.seed(7) as above).

test_that("ppmlfe's standard errors: nk, hc1, nw, boot", {
  d <- simpan()
  a <- function(...) {
    ppmlfe(d, "y", c("x", "w"), fe = c("unit", "day"), expo = "pop", ...)
  }
  expect_equal(
    a(ssc = "nk")$coef$se,
    c(0.02103345634587318, 0.020112441014339042),
    tolerance = 1e-6
  )
  h0 <- a(se = "hc1", ssc = "none")
  expect_equal(
    h0$coef$se,
    c(0.024527219446829767, 0.027257790635697002),
    tolerance = 1e-6
  )
  expect_equal(
    a(se = "hc1")$coef$se,
    h0$coef$se * sqrt(450 / 449),
    tolerance = 1e-12
  )
  nw <- a(se = "nw", ssc = "none", time = "day", lag = 2)
  expect_equal(
    nw$coef$se,
    c(0.023987274516971947, 0.027108298235302562),
    tolerance = 1e-6
  )
  expect_equal(
    a(se = "nw", ssc = "none", time = "day", lag = 0)$coef$se,
    h0$coef$se,
    tolerance = 1e-12
  )
  expect_identical(nw$lag, 2L)
  expect_error(a(se = "nw"), "needs `time`")
  b <- a(se = "boot", nboot = 30, seed = 7)
  expect_equal(
    b$coef$se,
    c(0.022722352431319974, 0.028702646182253613),
    tolerance = 1e-6
  )
  expect_identical(dim(b$draws), c(30L, 2L))
  expect_equal(
    a(ssc = "none")$coef$se,
    a()$coef$se * sqrt(14 / 15),
    tolerance = 1e-12
  )
  expect_error(a(se = "twoway", clus = "unit"), "needs two `clus` columns")
  expect_error(a(se = "robust"), "should be one of")
})

test_that("ppmlfe sep = \"none\" keeps a unit of zeros", {
  d <- simpan()
  z <- data.frame(
    unit = 16L,
    day = 1:30,
    z1 = 0,
    z2 = 0,
    w = seq(-1, 1, length.out = 30),
    x = seq(0, 2, length.out = 30),
    pop = 1000,
    y = 0,
    yl = 0
  )
  f <- suppressWarnings(ppmlfe(rbind(d, z), "y", c("x", "w"),
                               fe = c("unit", "day"), expo = "pop",
                               sep = "none", maxit = 200))
  expect_identical(f$n, 480L)
  expect_identical(f$nsep, 0L)
  expect_equal(
    f$coef$est,
    c(0.29195662923680077, 0.073073645297524289),
    tolerance = 1e-4
  )
})

test_that("ivfe's other errors agree with ppmlfe's machinery", {
  d <- simpan()
  f <- ivfe(
    d,
    "yl",
    endo = "x",
    inst = c("z1", "z2"),
    x = "w",
    fe = c("unit", "day"),
    se = "hc1"
  )
  expect_identical(f$se, "hc1")
  expect_true(all(is.finite(f$first$wald)))
  d$id <- seq_len(nrow(d))
  g <- ivfe(
    d,
    "yl",
    endo = "x",
    inst = c("z1", "z2"),
    x = "w",
    fe = c("unit", "day"),
    clus = "id"
  )
  expect_equal(f$coef$se, g$coef$se, tolerance = 1e-12)
  b <- ivfe(
    d,
    "yl",
    endo = "x",
    inst = c("z1", "z2"),
    x = "w",
    fe = c("unit", "day"),
    se = "boot",
    nboot = 5,
    seed = 1
  )
  expect_identical(nrow(b$draws), 5L)
})

test_that("ppmlcf: the plug-in of the paper, naive and bootstrap errors", {
  d <- simpan()
  a <- function(...) ppmlcf(d, "y", endo = "x", inst = c("z1", "z2"), x = "w",
                            fe = c("unit", "day"), expo = "pop", ...)
  p <- a(how = "plug", nboot = 30, seed = 7)
  expect_equal(
    p$coef$est,
    c(0.23021881035819208, 0.085582594321049033),
    tolerance = 1e-8
  )
  expect_equal(
    p$coef$se,
    c(0.03480437819221463, 0.033705894253545116),
    tolerance = 1e-6
  )
  expect_identical(p$coef$term, c("x", "w"))
  n <- a(how = "plug", se = "naive")
  expect_equal(n$coef$est, p$coef$est, tolerance = 1e-8)
  expect_equal(
    n$coef$se,
    c(0.033905920356246934, 0.02426990187010988),
    tolerance = 1e-6
  )
  cn <- a(se = "naive")
  expect_equal(
    cn$coef$est,
    c(0.21298786456069305, 0.09126882860013609, 0.26288312058292385),
    tolerance = 1e-8
  )
  expect_equal(p$first$eff, 536.97395780637555, tolerance = 1e-6)
  expect_error(a(how = "gmm"), "should be one of")
})

test_that("windiv: the paper's edges, hours to days", {
  expect_identical(
    unname(windiv(c(90, 91, 180, 181), edges = "paper")[, 1]),
    c(0, 1, 1, 0)
  )
  expect_identical(unname(windiv(c(90, 91, 180, 181))[, 1]), c(1, 1, 0, 0))
  v <- windiv(c(350, 10, 100, 120), by = c(1, 1, 2, 2))
  # 350 and 10 degrees average to north (0, or 360 less a rounding)
  expect_lt(min(attr(v, "dir")[1L], 360 - attr(v, "dir")[1L]), 1e-9)
  expect_equal(attr(v, "dir")[2L], 110, tolerance = 1e-12)
  expect_identical(rownames(v), c("1", "2"))
  expect_identical(unname(v[2, ]), c(1, 0, 0))
  m <- windiv(c(350, 10, 100, 120), by = c(1, 1, 2, 2), avg = "mean")
  expect_equal(attr(m, "dir"), c(180, 110))
  o <- windiv(c(100, 110, 200), by = c(1, 1, 1), avg = "mode")
  expect_identical(unname(o[1, ]), c(1, 0, 0))
  s <- windiv(c(100, 200, 200, 10), by = c(1, 1, 1, 1), share = TRUE)
  expect_equal(unname(s[1, ]), c(0.25, 0.5, 0))
  sa <- windiv(
    c(100, 200, 10, 10),
    by = c(1, 1, 2, 2),
    area = c("a", "a", "b", "b"),
    share = TRUE
  )
  expect_identical(dim(sa), c(2L, 6L))
  expect_equal(unname(sa[1, 1:3]), c(0.5, 0.5, 0))
  sp <- windiv(c(0, 90), by = c(1, 1), spd = c(1, 3))
  expect_equal(attr(sp, "dir"), atan2(3, 1) * 180 / pi, tolerance = 1e-12)
})

test_that("invsn: the sign, the indicator and bins", {
  r <- invsn(c(10, 5, 2), c(8, 6, 2), sign = "ground", brk = c(-Inf, 0, Inf))
  expect_identical(r$diff, c(2, -1, 0))
  expect_identical(r$inv, c(0L, 1L, 0L))
  expect_identical(as.integer(r$bin), c(2L, 1L, 1L))
  expect_error(invsn(1, 2, brk = 3), "two or more increasing")
  g <- attr(invsn(c(10, 5, NA), c(8, 6, 2)), "diag")
  expect_identical(g$n, 2L)
  expect_equal(g$zero, 0.5)
  expect_identical(c(g$min, g$max), c(-2, 1))
})

test_that("plcday moves the exposure, the instruments or both", {
  d <- simpan()
  a <- function(s) {
    plcday(d, "x", "unit", "day", shifts = 7, est = ivfe,
           with = c("z1", "z2"), shift = s, y = "yl", endo = "x",
           inst = c("z1", "z2"), x = "w", fe = c("unit", "day"))
  }
  b <- a("both")$tab
  x <- a("x")$tab
  z <- a("z")$tab
  expect_equal(b$est[1L], x$est[1L])
  expect_equal(b$est[1L], z$est[1L])
  d7 <- d[d$day <= 23, ]
  k <- match(paste(d7$unit, d7$day + 7), paste(d$unit, d$day))
  d7$z1 <- d$z1[k]
  d7$z2 <- d$z2[k]
  h <- ivfe(
    d7,
    "yl",
    endo = "x",
    inst = c("z1", "z2"),
    x = "w",
    fe = c("unit", "day")
  )
  expect_equal(z$est[2L], h$coef$est[1L])
  expect_false(isTRUE(all.equal(x$est[2L], b$est[2L])))
  expect_error(
    plcday(d, "x", "unit", "day", shift = "z", y = "y", x = "x", fe = "unit"),
    "needs the instruments"
  )
})

test_that("updwn: the outcome, the rings, two roads and the weather", {
  u <- simud()
  d <- u$d
  road <- data.frame(x0 = 0, y0 = 0, x1 = 5000, y1 = 0)
  f <- function(...) {
    updwn(d, "crimes", "ax", "ay", "day", "wind", road, band = 1600, ...)
  }
  s <- udhand(d, u$wnd)
  s$ly <- log1p(s$crimes)
  h <- function(y, est = ivfe) {
    est(s, y, x = "down", fe = c("north", "day"), clus = "cl")$coef$est
  }
  expect_equal(f(norm = "none")$coef$est, h("crimes"), tolerance = 1e-10)
  expect_equal(f(norm = "log")$coef$est, h("ly"), tolerance = 1e-10)
  expect_equal(f(norm = "pois")$coef$est,
               ppmlfe(s, "crimes", "down", fe = c("north", "day"),
                      clus = "cl")$coef$est, tolerance = 1e-8)
  r2 <- f(dist = c(0, 500, 1600))
  expect_identical(r2$coef$term[1L], "down:[0,500]")
  expect_identical(nrow(r2$coef), 2L)
  rw <- f(wx = "temp")
  expect_identical(rw$coef$term, "down")
  expect_true("temp:side2" %in% rw$fit$coef$term)
  expect_identical(f(se = "nw", lag = 1)$fit$se, "nw")
  # two roads end to end, one 1600 north of the other: areas near the gap
  # are within 1000 of both
  set.seed(13)
  a <- data.frame(area = 1:44,
                  ax = c(runif(20, 0, 5000), runif(20, 6000, 11000),
                         rep(5500, 4)),
                  ay = c(rep(c(-800, 800), 10), 1600 + rep(c(-800, 800), 10),
                         rep(800, 4)))
  dd <- merge(a, data.frame(day = 1:150))
  dd$wind <- u$wnd[dd$day]
  dd$crimes <- rpois(nrow(dd), 3)
  two <- data.frame(x0 = c(0, 6000), y0 = c(0, 1600), x1 = c(5000, 11000),
                    y1 = c(0, 1600))
  g <- function(tw) {
    updwn(dd, "crimes", "ax", "ay", "day", "wind", two, band = 1000,
          twin = tw)
  }
  r0 <- g("drop")
  rk <- g("keep")
  rn <- g("near")
  expect_identical(unname(r0$class[["twosources"]]), 4L * 150L)
  expect_identical(unname(rk$class[["twosources"]]), 0L)
  expect_identical(rk$n, r0$n)
  expect_identical(rn$n, r0$n)
  expect_false(isTRUE(all.equal(rk$coef$est, r0$coef$est)))
})

test_that("two-way factors each and the eigenvalue fix (CGM 2011)", {
  # twc = "each": V1 + V2 - V12, each one-way variance with its own G / (G
  # - 1), the intersections being the rows here
  d <- simpan()
  d$ud <- paste(d$unit, d$day)
  a <- function(...) {
    ppmlfe(d, "y", c("x", "w"), fe = c("unit", "day"), expo = "pop", ...)
  }
  fe2 <- a(clus = c("unit", "day"), twc = "each")
  v <- a(clus = "unit")$vcov + a(clus = "day")$vcov - a(clus = "ud")$vcov
  expect_equal(fe2$vcov, v, tolerance = 1e-10)
  fm <- a(clus = c("unit", "day"))
  fp <- a(clus = c("unit", "day"), psd = "eig")
  if (all(eigen(fm$vcov)$values >= 0)) {
    expect_equal(fp$vcov, fm$vcov, tolerance = 1e-10)
  }
  expect_true(all(eigen(fp$vcov)$values >= -1e-14))
})

test_that("the Newey-West (1994) automatic lag, by hand", {
  # one unit, 5 days, scores 1, -1, 2, 0, 1: n = [4 (5 / 100)^(2/9)] = 2,
  # sigma_0 = 7, sigma_1 = -3, sigma_2 = 4, s1 = 2 (-3 + 8) = 10, s0 = 7 +
  # 2 (1) = 9
  g <- crimlab:::.axnwlag(matrix(c(1, -1, 2, 0, 1)), rep(0L, 5), 1:5)
  expect_identical(g, floor(1.1447 * ((10 / 9)^2)^(1 / 3) * 5^(1 / 3)))
})

test_that("the Kleibergen-Paap rk statistic for two endogenous columns", {
  # the reference follows the paper's own route (Kleibergen and Paap
  # 2006, equations 11, 12, 19 and 22) in base R: the endogenous columns
  # and the instruments less w and the two sets of effects by lm(),
  # Theta = G Pi F' with F = (Z'Z)^(1/2) and G = (Y'Y)^(-1/2), svd(),
  # Lambda from U22 and V22, Omega from the unit-clustered covariance of
  # vec(Pi) times 15 / 14; the kernel takes a shorter route (the
  # orthonormal bases of the singular vectors), so agreement checks both
  d <- simpan()
  set.seed(7)
  d$z3 <- rnorm(nrow(d))
  d$x2 <- 0.4 * d$z2 + 0.5 * d$z3 + 0.2 * d$z1 + rnorm(nrow(d))
  f <- ivfe(d, "yl", endo = c("x", "x2"), inst = c("z1", "z2", "z3"),
            x = "w", fe = c("unit", "day"))
  expect_equal(f$first$kp, rep(79.3812891173097, 2), tolerance = 1e-6)
  expect_true(all(is.na(f$first$eff)))
})

test_that("sep = sst drops a regressor collinear on the positive counts", {
  # Santos Silva and Tenreyro (2010, p. 310-311): a dummy that is 0
  # wherever the count is positive and 1 on some zeros has no estimate;
  # their steps drop it and refit, which is the fit without it
  d <- simpan()
  z <- which(d$y == 0)
  d$bad <- 0
  d$bad[z[seq_len(min(5L, length(z)))]] <- 1
  skip_if(length(z) == 0L)
  expect_warning(
    f <- ppmlfe(d, "y", c("x", "bad", "w"), fe = c("unit", "day"),
                expo = "pop", sep = "sst"),
    "Santos Silva and Tenreyro"
  )
  expect_identical(f$xdrop, "bad")
  g <- ppmlfe(d, "y", c("x", "w"), fe = c("unit", "day"), expo = "pop")
  expect_equal(f$coef, g$coef, tolerance = 1e-12)
  expect_identical(g$xdrop, character())
})

test_that("ivfe gives Sargan and the clustered Hansen J (Yamada 2026)", {
  # the reference, base R matrix algebra: every column less its unit and
  # day means by lm() on the dummies, Z = (z1, z2, w), X = (x, w), 2SLS,
  # then S = crossprod(rowsum(Z * u, unit)) at its residuals, the two-step
  # estimate and J = u2'Z S^-1 Z'u2; Sargan n u'Pz u / u'u
  d <- simpan()
  f <- ivfe(d, "yl", endo = "x", inst = c("z1", "z2"), x = "w",
            fe = c("unit", "day"))
  expect_equal(f$overid$stat, c(4.94072255039924, 5.30458827826541),
               tolerance = 1e-8)
  expect_identical(f$overid$df, c(1, 1))
  expect_equal(unname(f$gmm), c(0.563607320395195, 0.213473897850748),
               tolerance = 1e-8)
  expect_output(print(f), "Hansen J")
  f1 <- ivfe(d, "yl", endo = "x", inst = "z1", x = "w", fe = c("unit", "day"))
  expect_null(f1$overid)
})

test_that("se = quasi is glm()'s quasi-Poisson with dummy columns", {
  # the reference: glm(y ~ x + w + factor(unit) + factor(day) +
  # offset(log(pop)), family = quasipoisson), epsilon = 1e-14; its
  # dispersion is Pearson's chi-squared over the residual degrees of
  # freedom, as Armstrong, Gasparrini and Tobias (2014, p. 4) estimate it
  d <- simpan()
  f <- ppmlfe(d, "y", c("x", "w"), fe = c("unit", "day"), expo = "pop",
              se = "quasi")
  expect_equal(f$coef$se, c(0.0250069325985391, 0.0308754583727149),
               tolerance = 1e-6)
  expect_equal(f$psi, 1.20631966625669, tolerance = 1e-7)
})

test_that("a third set of effects equals its indicator columns", {
  # the references: glm(poisson) and lm() with all three sets as factor
  # columns (epsilon = 1e-14); absorbing changes only the arithmetic
  d <- simpan()
  set.seed(11)
  d$h <- sample(1:4, nrow(d), TRUE)
  f <- ppmlfe(d, "y", c("x", "w"), fe = c("unit", "day", "h"), expo = "pop")
  expect_equal(f$coef$est, c(0.296241030216946, 0.0679892442511943),
               tolerance = 1e-8)
  g <- ivfe(d, "yl", x = c("x", "w"), fe = c("unit", "day", "h"))
  expect_equal(g$coef$est, c(0.809239042169907, 0.144096303385526),
               tolerance = 1e-8)
  # a group of the third set without counts is dropped
  d$y[d$h == 4] <- 0
  expect_warning(f4 <- ppmlfe(d, "y", "x", fe = c("unit", "day", "h")),
                 "dropped")
  expect_identical(f4$nsep, sum(d$h == 4))
})
