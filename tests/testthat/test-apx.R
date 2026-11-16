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
  expect_equal(f$coef$rr, exp(f$coef$est), tolerance = 1e-12)
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
  expect_equal(f$fstab$rel, f$fstab$est / mean(d$x), tolerance = 1e-12)
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
  expect_equal(p$tab$est[1L], f$coef$est[1L], tolerance = 1e-12)
  expect_identical(p$tab$n, c(450L, 345L, 345L))
  # by hand: the exposure of a week later
  d7 <- d[d$day <= 23, ]
  d7$x <- d$x[match(paste(d7$unit, d7$day + 7), paste(d$unit, d$day))]
  f7 <- ppmlfe(d7, "y", c("x", "w"), fe = c("unit", "day"), expo = "pop")
  expect_equal(p$tab$est[3L], f7$coef$est[1L], tolerance = 1e-12)
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
  expect_equal(rp$placebo$est[2L], r$coef$est, tolerance = 1e-12)
  expect_true(all(is.finite(rp$placebo$est)))
  expect_error(updwn(d, "crimes", "ax", "ay", "day", "wind", road, tol = 100),
               "between 0 and 90")
})

test_that("expjn joins, counts misses and refuses duplicate keys", {
  ex <- data.frame(area = c("a", "a", "b"), day = c(1, 2, 1), pm = c(10, 12, 8),
                   stringsAsFactors = FALSE)
  oc <- data.frame(
    area = c("a", "a", "b", "b"),
    day = c(1, 2, 1, 2),
    crimes = c(3, 1, 0, 2),
    stringsAsFactors = FALSE
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
  expect_error(a(se = "robust"), "^ppmlfe: `se` must be one of \"cluster\"")
  expect_error(a(ssc = NA), "^ppmlfe: `ssc` must be one of .*, not NA")
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
  a <- function(...) {
    ppmlcf(d, "y", endo = "x", inst = c("z1", "z2"), x = "w",
           fe = c("unit", "day"), expo = "pop", ...)
  }
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
  expect_error(a(how = "gmm"), "^ppmlcf: `how` must be one of")
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
  expect_equal(attr(m, "dir"), c(180, 110), tolerance = 1e-12)
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

test_that("invsn by day: Batkeyev and DeRemer's every-reading rule", {
  # three-hourly readings of three days: day 1 inverted at every reading,
  # day 2 at two of three, day 3 at one reading with two missing; by hand
  lo <- c(1, 2, 3, 1, 2, 3, 1, NA, NA)
  hi <- c(4, 5, 6, 4, 1, 6, 2, 5, NA)
  d <- rep(c("d1", "d2", "d3"), each = 3)
  r <- invsn(lo, hi, by = d)
  expect_identical(r$by, c("d1", "d2", "d3"))
  expect_identical(r$n, c(3L, 3L, 1L))
  expect_identical(r$inv, c(1L, 0L, 1L))
  expect_equal(r$diff, c(3, (3 - 1 + 3) / 3, 1), tolerance = 1e-15)
  expect_identical(invsn(lo, hi, by = d, rule = "any")$inv, c(1L, 1L, 1L))
  expect_identical(invsn(lo, hi, by = d, rule = "mean")$inv, c(1L, 1L, 1L))
  g <- invsn(lo, hi, sign = "ground", by = d, brk = c(-Inf, -2, Inf))
  expect_equal(g$diff, -r$diff, tolerance = 1e-15)
  expect_identical(as.character(g$bin),
                   c("(-Inf,-2]", "(-2, Inf]", "(-2, Inf]"))
  dg <- attr(r, "diag")
  expect_identical(dg$n, 3L)
  expect_equal(dg$zero, 1 / 3, tolerance = 1e-15)
  # a day without a reading is NA
  e <- invsn(c(1, NA), c(2, 3), by = c(1, 2))
  expect_identical(e$inv, c(1L, NA))
  expect_identical(e$n, c(1L, 0L))
  expect_error(invsn(1, 2, by = c(1, 2)),
               "invsn: `by` must give the group of each reading, with no NA.",
               fixed = TRUE)
  expect_error(invsn(1, 2, by = 1, rule = "some"), "invsn: `rule` must be",
               fixed = TRUE)
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
  expect_equal(b$est[1L], x$est[1L], tolerance = 1e-12)
  expect_equal(b$est[1L], z$est[1L], tolerance = 1e-12)
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
  expect_equal(z$est[2L], h$coef$est[1L], tolerance = 1e-12)
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
  # One unit, 5 days, scores 1, -1, 2, 0, 1. The selection parameter is
  # the whole part of 4 (5 / 100)^(2/9), which is 2. The autocovariances
  # at lags 0, 1 and 2 are 7, -3 and 4, so s1 is twice (-3 + 8), which is
  # 10, and s0 is 7 plus twice 1, which is 9.
  g <- .axnwlag(matrix(c(1, -1, 2, 0, 1)), rep(0L, 5), 1:5)
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

test_that("sep = \"sst\" drops a regressor the effects absorb where y > 0", {
  # On the rows with positive counts z1 is a constant of the area and z2 a
  # constant of the area plus one of the day, so each is a combination of
  # the effects there; both move only where the count is zero, the case
  # of Santos Silva and Tenreyro (2010). Demeaning leaves rounding noise
  # of them, not exact zeros.
  set.seed(1)
  d <- expand.grid(area = letters[1:8], day = 1:30, stringsAsFactors = FALSE)
  d$pm <- rnorm(nrow(d))
  az <- c(a = 0.3, b = 1.7, c = -0.4, d = 2.2, e = 0.9, f = -1.1, g = 0.05,
          h = 3.3)
  d$y <- rpois(nrow(d), exp(0.2 * d$pm + 0.5))
  d$z1 <- az[d$area] + (d$y == 0)
  d$z2 <- az[d$area] + sin(d$day) + (d$y == 0)
  expect_gt(sum(d$y == 0), 0)
  expect_warning(
    r1 <- ppmlfe(d, "y", c("pm", "z1"), fe = "area", sep = "sst", se = "hc1"),
    "regressor\\(s\\) dropped"
  )
  expect_identical(r1$xdrop, "z1")
  expect_identical(r1$coef$term, "pm")
  expect_warning(
    r2 <- ppmlfe(d, "y", c("pm", "z2"), fe = c("area", "day"), sep = "sst",
                 se = "hc1"),
    "regressor\\(s\\) dropped"
  )
  expect_identical(r2$xdrop, "z2")
  # a regressor with its own variation after the effects stays
  r3 <- ppmlfe(d, "y", "pm", fe = c("area", "day"), sep = "sst", se = "hc1")
  expect_identical(r3$xdrop, character())
  expect_equal(r2$coef$est, r3$coef$est, tolerance = 1e-10)
})

# Third round: the refusals, the printed summaries and the edge cases of
# the helpers, each reached through the public functions where it can be.

test_that("the refusals of the column and data checks", {
  d <- simpan()
  d$s <- "a"
  d$one <- 1
  d$dc <- as.character(d$day)
  dn <- d
  dn$x[3] <- NA
  di <- d
  di$x[1] <- Inf
  d0 <- d
  d0$y <- 0
  dp <- d
  dp$pop[1] <- 0
  errs(
    "ppmlfe",
    list(
      list("`data` must be a data frame, got list", data = list(a = 1),
           y = "y", x = "x", fe = "unit"),
      list("`data` has no rows", data = d[0, ], y = "y", x = "x",
           fe = "unit"),
      list("`y` must name columns of `data`, got NULL", data = d, y = NULL,
           x = "x", fe = "unit"),
      list("`y` must name 1 column\\(s\\)", data = d, y = c("y", "yl"),
           x = "x", fe = "unit"),
      list("`x` names a column twice", data = d, y = "y", x = c("x", "x"),
           fe = "unit"),
      list("column `x` has missing values", data = dn, y = "y", x = "x",
           fe = "unit"),
      list("column `s` must be numeric, got character", data = d, y = "y",
           x = "s", fe = "unit"),
      list("column `x` must be finite", data = di, y = "y", x = "x",
           fe = "unit"),
      list("`level` must be one number strictly between 0 and 1", data = d,
           y = "y", x = "x", fe = "unit", level = 1),
      list("column `dc` must hold numbers or dates", data = d, y = "y",
           x = "x", fe = "unit", se = "nw", time = "dc"),
      list("the bootstrap needs two or more clusters in `one`", data = d,
           y = "y", x = "x", fe = "unit", se = "boot", clus = "one",
           nboot = 2),
      list("column `y` is all 0", data = d0, y = "y", x = "x", fe = "unit"),
      list("column `pop` must be positive", data = dp, y = "y", x = "x",
           fe = "unit", expo = "pop")
    )
  )
  # a regressor that is 0 wherever the count is positive is all there is
  d$bad <- as.numeric(d$y == 0)
  expect_error(
    ppmlfe(d, "y", "bad", fe = c("unit", "day"), sep = "sst"),
    "every regressor is collinear on the rows with counts"
  )
})

test_that("ivfe, ppmlcf and expgmm refuse a model they cannot fit", {
  d <- simpan()
  d$one <- 1
  d$y0 <- 0
  d$yn <- d$y
  d$yn[1] <- -1
  d$p0 <- d$pop
  d$p0[1] <- 0
  errs(
    "ivfe",
    list(
      list("`se = \"twoway\"` needs two `clus` columns", data = d, y = "yl",
           endo = "x", inst = "z1", fe = "unit", se = "twoway",
           clus = "unit"),
      list("`inst` needs `endo`: the columns it instruments", data = d,
           y = "yl", inst = "z1", x = "w", fe = "unit"),
      list("give `endo` or `x`: there is nothing to estimate", data = d,
           y = "yl", fe = "unit"),
      list("a column cannot be both a regressor and an excluded instrument",
           data = d, y = "yl", endo = "x", inst = c("z1", "w"), x = "w",
           fe = "unit"),
      list("`clus` has one cluster: a cluster-robust variance needs two",
           data = d, y = "yl", endo = "x", inst = "z1", fe = "unit",
           clus = "one")
    )
  )
  errs(
    "ppmlcf",
    list(
      list("at least as many instruments as `endo` has columns", data = d,
           y = "y", endo = c("x", "w"), inst = "z1", fe = "unit"),
      list("a column cannot be both a regressor and an excluded instrument",
           data = d, y = "y", endo = "x", inst = "z1", x = "z1",
           fe = "unit"),
      list("column `yn` must hold counts, at least 0", data = d, y = "yn",
           endo = "x", inst = "z1", fe = "unit"),
      list("column `y0` is all 0", data = d, y = "y0", endo = "x",
           inst = "z1", fe = "unit"),
      list("column `p0` must be positive", data = d, y = "y", endo = "x",
           inst = "z1", fe = "unit", expo = "p0")
    )
  )
  d$yp <- d$y + 1
  d$id <- seq_len(nrow(d))
  errs(
    "expgmm",
    list(
      list("at least as many instruments as `endo` has columns", data = d,
           y = "y", endo = c("x", "w"), inst = "z1", clus = "unit"),
      list("column `yn` must hold counts, at least 0", data = d, y = "yn",
           endo = "x", inst = "z1", clus = "unit"),
      list("column `y0` is all 0", data = d, y = "y0", endo = "x",
           inst = "z1", clus = "unit"),
      list("column `p0` must be positive", data = d, y = "y", endo = "x",
           inst = "z1", clus = "unit", expo = "p0"),
      # every count positive: no row is separated, all 450 levels enter
      list("the fixed effects have 450 levels; they enter as indicator",
           data = d, y = "yp", endo = "x", inst = "z1", fe = "id",
           se = "hc1")
    )
  )
})

test_that("the printed summaries and summary() of the estimators", {
  d <- simpan()
  iv <- ivfe(d, "yl", endo = "x", inst = "z1", fe = c("unit", "day"))
  expect_identical(summary(iv), iv$coef)
  cf <- ppmlcf(d, "y", endo = "x", inst = "z1", fe = c("unit", "day"),
               nboot = 3, seed = 1)
  expect_identical(summary(cf), cf$coef)
  expect_output(
    print(cf),
    paste0("ppmlcf: Poisson control function, 450 rows, fixed effects ",
           "unit and day; 3 bootstrap draws of 15 clusters \\(unit\\)")
  )
  expect_output(print(cf), "first stage: robust Wald F")
  pn <- ppmlcf(d, "y", endo = "x", inst = "z1", fe = c("unit", "day"),
               how = "plug", se = "naive")
  expect_output(
    print(pn),
    "Poisson plug-in, .*naive second-stage errors \\(they ignore"
  )
  # a unit of zeros, dropped and reported
  z <- data.frame(unit = 16L, day = 1:30, z1 = 0, z2 = 0, w = 0, x = 0,
                  pop = 1000, y = 0, yl = 0)
  expect_warning(
    f <- ppmlfe(rbind(d, z), "y", "x", fe = c("unit", "day")),
    "30 row\\(s\\) dropped"
  )
  expect_output(print(f), "30 row\\(s\\) dropped for groups with no counts")
})

test_that("ppmlfe and expgmm warn when the fit does not converge", {
  d <- simpan()
  expect_warning(
    ppmlfe(d, "y", c("x", "w"), fe = c("unit", "day"), maxit = 1),
    "the fit did not converge in 1 iterations"
  )
  # iterated GMM: the second weight is already past maxit = 1
  expect_warning(
    g <- expgmm(d, "y", endo = "x", inst = c("z1", "z2"), x = "w",
                fe = "day", se = "hc1", step = "iter", maxit = 1),
    "the fit did not converge in 1 iterations"
  )
  expect_false(g$conv)
})

test_that("expgmm: robust errors, an offset of 0, a dropped day", {
  d <- simpan()
  d$one <- 1
  f0 <- expgmm(d, "y", endo = "x", inst = c("z1", "z2"), x = "w",
               fe = "day", se = "hc1")
  f1 <- expgmm(d, "y", endo = "x", inst = c("z1", "z2"), x = "w",
               fe = "day", se = "hc1", expo = "one")
  # an exposure of 1 is an offset of 0: the same fit
  expect_equal(f1$coef$est, f0$coef$est, tolerance = 1e-10)
  expect_output(print(f0), "heteroskedasticity-robust errors")
  d$y[d$day == 1] <- 0
  expect_warning(
    g <- expgmm(d, "y", endo = "x", inst = c("z1", "z2"), x = "w",
                fe = "day", se = "hc1", expo = "pop"),
    "15 row\\(s\\) dropped"
  )
  expect_identical(g$nsep, 15L)
  expect_output(print(g), "15 row\\(s\\) dropped for groups with no counts")
  # the regressors collinear: the kernel's refusal
  d <- simpan()
  d$x2 <- d$x
  expect_error(
    expgmm(d, "y", endo = "x", inst = c("z1", "z2"), x = "x2", fe = "day",
           se = "hc1"),
    "the fit failed: the instruments or the regressors are collinear"
  )
})

test_that("factor effects, two-way nk factors each, the automatic NW lag", {
  d <- simpan()
  a <- function(...) {
    ppmlfe(d, "y", c("x", "w"), expo = "pop", ...)
  }
  d$uf <- factor(d$unit)
  expect_equal(a(fe = c("uf", "day"))$coef$est,
               a(fe = c("unit", "day"))$coef$est, tolerance = 1e-12)
  # both effects are nested in a cluster column, so k = 2 + 1 and each of
  # the three components takes (n - 1) / (n - k) = 449 / 447 more
  ge <- a(fe = c("unit", "day"), clus = c("unit", "day"), twc = "each")
  ne <- a(fe = c("unit", "day"), clus = c("unit", "day"), twc = "each",
          ssc = "nk")
  expect_equal(ne$vcov, ge$vcov * 449 / 447, tolerance = 1e-10)
  nw <- a(fe = c("unit", "day"), se = "nw", time = "day", ssc = "none")
  expect_gte(nw$lag, 0L)
  expect_equal(
    a(fe = c("unit", "day"), se = "nw", time = "day", ssc = "none",
      lag = nw$lag)$coef$se,
    nw$coef$se,
    tolerance = 1e-12
  )
  # scores of zero leave no ratio: the lag is 0
  expect_identical(.axnwlag(matrix(0, 5L, 1L), rep(0L, 5), 1:5), 0)
})

test_that("the bootstrap: draws refitted with sep = none, failed draws", {
  d <- simpan()
  # every count positive: no group is separated in any draw
  d$y <- d$y + 1
  b1 <- ppmlfe(d, "y", "x", fe = c("unit", "day"), se = "boot", nboot = 3,
               seed = 2)
  b0 <- ppmlfe(d, "y", "x", fe = c("unit", "day"), se = "boot", nboot = 3,
               seed = 2, sep = "none")
  expect_equal(b0$draws, b1$draws, tolerance = 1e-12)
  # a fit that fails on the first draw only; the sum of the weights of a
  # draw is 450 rows whatever is drawn (15 clusters of 30 rows)
  st <- new.env()
  st$i <- 0L
  fit <- function(rows, w) {
    st$i <- st$i + 1L
    if (st$i == 1L) stop("no fit")
    as.numeric(sum(w))
  }
  expect_warning(
    b <- .axboot(d, "unit", 4L, 1, 1L, fit, "ppmlfe"),
    "ppmlfe: 1 of 4 bootstrap draws could not be fitted and are left out"
  )
  expect_identical(b$nfail, 1L)
  expect_true(is.na(b$draws[1L, 1L]))
  expect_equal(b$draws[-1L, 1L], rep(450, 3), tolerance = 1e-12)
  expect_equal(b$se, 0, tolerance = 1e-12)
  expect_error(
    suppressWarnings(
      .axboot(d, "unit", 3L, 1, 1L, function(rows, w) stop("no"), "ppmlfe")
    ),
    "fewer than two bootstrap draws could be fitted"
  )
})

test_that("ppmlcf leaves out the draws of a single cluster", {
  d <- simpan()
  d <- d[d$unit <= 2L, ]
  quiet <- function(expr) {
    st <- new.env()
    st$w <- character()
    r <- withCallingHandlers(expr, warning = function(e) {
      st$w <- c(st$w, conditionMessage(e))
      invokeRestart("muffleWarning")
    })
    list(r = r, w = st$w)
  }
  # with two clusters a draw of one of them twice has one unit, and the
  # day effects then absorb every column: that draw cannot be fitted
  bm <- .axdraw(2L, 20L, 1)
  nf <- sum(colSums(bm > 0L) == 1L)
  o <- quiet(ppmlcf(d, "y", endo = "x", inst = "z1",
                    fe = c("unit", "day"), nboot = 20, seed = 1))
  expect_identical(o$r$nfail, as.integer(nf))
  expect_true(any(grepl(
    paste(nf, "of 20 bootstrap draws could not be fitted"), o$w,
    fixed = TRUE
  )))
  # a seed whose two draws both take one cluster only
  s <- Find(function(k) all(colSums(.axdraw(2L, 2L, k) > 0L) == 1L), 1:200)
  expect_error(
    quiet(ppmlcf(d, "y", endo = "x", inst = "z1", fe = c("unit", "day"),
                 nboot = 2, seed = s)),
    "fewer than two bootstrap draws could be fitted"
  )
})

test_that("the seed of a draw leaves no generator where there was none", {
  env <- globalenv()
  had <- exists(".Random.seed", envir = env, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = env, inherits = FALSE)
  on.exit(if (had) assign(".Random.seed", old, envir = env), add = TRUE)
  if (had) rm(".Random.seed", envir = env)
  m <- .axdraw(3L, 2L, 1)
  expect_false(exists(".Random.seed", envir = env, inherits = FALSE))
  expect_identical(colSums(m), c(3, 3))
})

test_that("the Kleibergen-Paap rk statistic with two-way clusters", {
  # with the rows as the second clustering the two-way meat is that of
  # the first: V(unit) + V(row) - V(unit and row), the last two the same
  # partition, and the factor is that of the 15 units either way
  d <- simpan()
  set.seed(7)
  d$z3 <- rnorm(nrow(d))
  d$x2 <- 0.4 * d$z2 + 0.5 * d$z3 + 0.2 * d$z1 + rnorm(nrow(d))
  d$id <- seq_len(nrow(d))
  a <- function(cl) {
    ivfe(d, "yl", endo = c("x", "x2"), inst = c("z1", "z2", "z3"),
         x = "w", fe = c("unit", "day"), clus = cl)$first$kp
  }
  expect_equal(a(c("unit", "id")), a("unit"), tolerance = 1e-8)
})

test_that("invsn and plcday: the refusals, print and summary", {
  expect_error(invsn("a", 1), "must be numeric vectors of one length")
  d <- simpan()
  d$dc <- as.character(d$day)
  # an estimator by hand: the mean of the exposure on the rows it gets
  mx <- function(data, ...) {
    list(coef = data.frame(term = "x", est = mean(data$x), se = 1, lo = 0,
                           hi = 1, p = 0.5))
  }
  p <- plcday(d, "x", "unit", "day", shifts = 7, est = mx)
  expect_equal(p$tab$est[2L], mean(d$x[d$day > 7]), tolerance = 1e-12)
  expect_identical(summary(p), p$tab)
  expect_output(print(p), "plcday: the coefficient of x with values of other")
  nox <- function(data, ...) list(coef = data.frame(term = "w"))
  errs(
    "plcday",
    list(
      list("`est` must be an estimator function", data = d, var = "x",
           unit = "unit", day = "day", est = "ppmlfe"),
      list("`shifts` must be distinct whole numbers of days other than 0",
           data = d, var = "x", unit = "unit", day = "day", shifts = 0,
           est = mx),
      list("`shifts` must be distinct whole numbers of days other than 0",
           data = d, var = "x", unit = "unit", day = "day", shifts = 1.5,
           est = mx),
      list("column `dc` must hold numbers or dates", data = d, var = "x",
           unit = "unit", day = "dc", est = mx),
      list("a shift of 100 days leaves fewer than two rows", data = d,
           var = "x", unit = "unit", day = "day", shifts = 100, est = mx),
      list("`est` returned no coefficient for `x`", data = d, var = "x",
           unit = "unit", day = "day", shifts = 7, est = nox)
    )
  )
})

test_that("updwn: print, summary, the placebo without areas, clusters", {
  u <- simud()
  d <- u$d
  road <- data.frame(x0 = 0, y0 = 0, x1 = 5000, y1 = 0)
  f <- function(...) {
    updwn(d, "crimes", "ax", "ay", "day", "wind", road, band = 1600, ...)
  }
  r <- f()
  expect_identical(summary(r), r$coef)
  expect_output(print(r), "updwn: downwind against upwind, .* over its mean")
  # a road moved 50 km has no area near it: no estimate, a warning
  expect_warning(rp <- f(offs = c(0, 50000)),
                 "no estimate for the segments moved by 50000")
  expect_true(is.na(rp$placebo$est[2L]))
  expect_identical(rp$placebo$n[2L], 0L)
  expect_output(print(rp), "placebo segments moved sideways")
  s <- summary(rp)
  expect_identical(s$off, c(NA, 0, 50000))
  expect_equal(s$est[1:2], c(r$coef$est, r$coef$est), tolerance = 1e-12)
  # the day as the cluster: one per day with the wind across the road
  rc <- f(se = "cluster", clus = "day")
  expect_identical(rc$fit$ncl, sum(abs(cos(u$wnd * pi / 180)) >= 0.5))
  expect_equal(rc$coef$est, r$coef$est, tolerance = 1e-12)
})

test_that("updwn refuses what it cannot classify or fit", {
  u <- simud()
  d <- u$d
  road <- data.frame(x0 = 0, y0 = 0, x1 = 5000, y1 = 0)
  pt <- data.frame(x0 = 2500, y0 = 0)
  dn <- d
  dn$crimes[1] <- -1
  dc <- d
  dc$dc <- as.character(dc$day)
  d0 <- d
  d0$crimes <- 0
  # half the areas of each side with the wind turned round
  dw <- d
  dw$wind <- dw$wind + 180 * (dw$area %% 4 < 2)
  a <- function(...) {
    c(list(y = "crimes", ax = "ax", ay = "ay", day = "day", wind = "wind"),
      list(...))
  }
  errs(
    "updwn",
    list(
      c("`src` must have columns x0 and y0", a(data = d, src = road[1])),
      c("`src` column x0 must be numeric",
        a(data = d, src = data.frame(x0 = "a", y0 = 0,
                                     stringsAsFactors = FALSE))),
      c("`src` x0 and y0 must be finite",
        a(data = d, src = data.frame(x0 = NA_real_, y0 = 0))),
      c("`dist` must be two or more increasing distances from 0",
        a(data = d, src = road, dist = 5)),
      c("`agg` must be TRUE or FALSE", a(data = d, src = road, agg = NA)),
      c("`area` is needed for point sources and for agg = FALSE",
        a(data = d, src = pt)),
      c("`offs` moves segments: every source must be a segment",
        a(data = d, src = pt, area = "area", offs = 1)),
      c("`offs` must be finite numbers", a(data = d, src = road, offs = Inf)),
      c("`se = \"twoway\"` needs two `clus` columns",
        a(data = d, src = road, se = "twoway")),
      c("column `crimes` must hold counts, at least 0", a(data = dn,
                                                          src = road)),
      c("`se = \"nw\"` needs a day column of numbers or dates",
        utils::modifyList(a(data = dc, src = road, se = "nw"),
                          list(day = "dc"))),
      c("an area has two rows on one day",
        a(data = rbind(d, d[1, ]), src = road, area = "area")),
      c("fewer than two rows are upwind or downwind of a source",
        a(data = d, src = road, band = 1)),
      c("no row falls in the rings of `dist`",
        a(data = d, src = road, band = 1600, dist = c(0, 1e-3))),
      c("the wind differs between areas of one side and day",
        a(data = dw, src = road, band = 1600)),
      c("column `area` must be constant within each side \\(area\\) and day",
        a(data = d, src = road, band = 1600, clus = "area")),
      c("the counts are all 0", a(data = d0, src = road, band = 1600))
    )
  )
})

test_that("the kernel of updwn: missing values, on the source", {
  # an area on the road and one with no wind; then a point source and an
  # area on it, and one beside it with the wind along the line between
  k <- .Call(C_clab_r_updwn, c(2500, 100), c(0, 800), c(0, NA_real_), 0, 0,
             5000, 0, Inf, 60, 0L)
  expect_identical(k$dn, c(-4L, -5L))
  k <- .Call(C_clab_r_updwn, c(0, 100), c(0, 0), c(0, 0), 0, 0, NA_real_,
             NA_real_, Inf, 60, 0L)
  expect_identical(k$dn, c(-4L, -1L))
})

test_that("windiv: the refusals, modes, missing groups and factor areas", {
  errs(
    "windiv",
    list(
      list("`dir` must be numeric directions in degrees", dir = "a"),
      list("`share` must be TRUE or FALSE", dir = 1, share = NA),
      list("`spd` must give one speed from 0 per direction", dir = 1,
           spd = -1),
      list("`area` must give one label per direction, none missing",
           dir = c(1, 2), area = "a"),
      list("`sect` must be two different finite directions in degrees",
           dir = 1, sect = c(10, 370)),
      list("`by` must give one group per direction, none missing",
           dir = c(1, 2), by = 1),
      list("`area` must be constant within each `by` group", dir = c(1, 2),
           by = c(1, 1), area = c("a", "b")),
      list("`avg = \"mode\"` needs sectors by `width`, not `sect`",
           dir = c(1, 2), by = c(1, 1), avg = "mode", sect = c(0, 90)),
      list("`share` needs sectors by `width`, not `sect`", dir = c(1, 2),
           by = c(1, 1), share = TRUE, sect = c(0, 90))
    )
  )
  # the mode under the paper's edges: 90 is in the first sector, (0, 90]
  o <- windiv(c(90, 90, 200), by = c(1, 1, 1), avg = "mode", edges = "paper")
  expect_identical(unname(o[1, ]), c(0, 0, 0))
  h <- windiv(c(90, 90, 200), by = c(1, 1, 1), avg = "mode")
  expect_identical(unname(h[1, ]), c(1, 0, 0))
  # a group without a direction has no mode
  m <- windiv(c(NA, 100), by = c(1, 2), avg = "mode")
  expect_true(all(is.na(m[1, ])))
  expect_identical(unname(m[2, ]), c(1, 0, 0))
  # a factor's levels keep their order, unused ones dropped
  fa <- factor(c("b", "a"), levels = c("c", "b", "a"))
  f <- windiv(c(10, 100), area = fa)
  expect_identical(colnames(f)[c(1, 4)], c("w90_180:b", "w90_180:a"))
  expect_identical(unname(f[2, 4]), 1)
})

test_that("expjn: the refusals of the keys and summary()", {
  ex <- data.frame(area = c("a", "a", "b"), day = c(1, 2, 1),
                   pm = c(10, 12, 8), stringsAsFactors = FALSE)
  oc <- data.frame(area = c("a", "a", "b", "b"), day = c(1, 2, 1, 2),
                   crimes = c(3, 1, 0, 2), stringsAsFactors = FALSE)
  j <- expjn(ex, oc)
  expect_identical(summary(j),
                   data.frame(joined = 3L, nout = 1L, nexp = 0L))
  en <- ex
  en$day[1] <- NA
  oo <- oc
  oo$pm <- 1
  errs(
    "expjn",
    list(
      list("`by` must name the key columns", expo = ex, outc = oc,
           by = character()),
      list("`expo` lacks the key columns: day", expo = ex[c("area", "pm")],
           outc = oc),
      list("`expo` has missing keys", expo = en, outc = oc),
      list("`expo` has no columns besides the keys",
           expo = ex[c("area", "day")], outc = oc),
      list("columns in both tables besides the keys: pm", expo = ex,
           outc = oo),
      list("no outcome row has an exposure: check that the keys agree",
           expo = transform(ex, day = day + 10), outc = oc)
    )
  )
})

test_that("the kernel of exprun: a missing exposure ends a run", {
  # exprun() refuses missing values before it; the kernel marks the row
  r <- .Call(C_clab_r_runlen, rep(0L, 4), c(1, 2, 3, 4), c(1, NA, 1, 1), 1)
  expect_identical(r, c(1L, NA, 1L, 2L))
})

test_that("windiv() refuses more sectors than its bound", {
  # 360 / 2^27 divides 360 exactly, into 2^27 sectors
  expect_error(windiv(c(10, 100), width = 360 / 2^27),
               "windiv: `width` must divide 360 into at most 67108864 sectors.",
               fixed = TRUE)
})

test_that("the bootstrap made a draw at a time equals the stored draws", {
  # one block of every draw is the clusters x nboot matrix the bootstrap
  # used to hold; blocks of one draw are what a large one now makes
  d <- simpan()
  d$y <- d$y + 1
  boot <- function(blk) {
    old <- options(crimlab.drawblock = blk)
    on.exit(options(old))
    ppmlfe(d, "y", "x", fe = c("unit", "day"), se = "boot", nboot = 7,
           seed = 2)
  }
  expect_identical(boot(1), boot(1e9))
  expect_identical(boot(3), boot(1e9))
  cf <- function(blk) {
    old <- options(crimlab.drawblock = blk)
    on.exit(options(old))
    suppressWarnings(ppmlcf(d, "y", endo = "x", inst = "z1",
                            fe = c("unit", "day"), nboot = 7, seed = 1))
  }
  expect_identical(cf(1), cf(1e9))
  expect_identical(cf(2), cf(1e9))
})

test_that("every choice argument is refused in the package's words", {
  cs <- list(
    ppmlfe = c("se", "ssc", "twc", "psd", "sep"),
    ivfe = c("se", "ssc", "twc", "psd"), ppmlcf = c("how", "se", "ssc"),
    expgmm = c("err", "step", "se", "ssc"), windiv = c("edges", "avg"),
    invsn = "sign", plcday = "shift", updwn = c("twin", "norm", "se", "ssc"),
    smkiv = c("adj", "cf"), hetfx = "how", ccross = "strat", jntx = "form",
    typefx = "adj", rdloc = c("bw", "kern", "ci", "se", "ssc")
  )
  for (fn in names(cs)) {
    for (a in cs[[fn]]) {
      for (v in list("zz", NA, 1)) {
        m <- paste0("^", fn, ": `", a, "` must be one of ")
        expect_error(do.call(fn, stats::setNames(list(v), a)), m,
                     info = paste(fn, a))
      }
    }
  }
})

test_that("an iteration count past an integer is refused without a warning", {
  m <- "`maxit` must be between 1 and 2147483647"
  expect_error(expect_warning(ppmlfe(simpan(), "y", "x", fe = "unit",
                                     maxit = 1e300), NA), m)
})

test_that("a bootstrap seed past an integer is refused in words", {
  d <- simpan()
  m <- "`seed` must be between 0 and 2147483647"
  expect_error(ppmlcf(d, "y", endo = "x", inst = c("z1", "z2"), x = "w",
                      fe = c("unit", "day"), expo = "pop", se = "boot",
                      nboot = 9, seed = 1e300), m)
})

test_that("ivfe(se = \"dk\") is Driscoll and Kraay's variance written out", {
  # least squares with unit and day effects: the demeaned regressor x~
  # (residuals of x on the effects), the residuals e, the scores x~ e
  # summed by day, h_t, and the Bartlett-weighted sum of their products
  # up to lag days apart, over (x~'x~)^2
  set.seed(61)
  d <- expand.grid(u = 1:12, t = 1:15)
  d$x <- rnorm(nrow(d)) + rnorm(15)[d$t]
  d$y <- 0.5 * d$x + rnorm(15)[d$t] + rnorm(nrow(d))
  r <- ivfe(d, "y", x = "x", fe = c("u", "t"), se = "dk", time = "t",
            lag = 2, ssc = "none")
  xt <- stats::residuals(stats::lm(x ~ factor(u) + factor(t), d))
  e <- stats::residuals(stats::lm(y ~ x + factor(u) + factor(t), d))
  h <- as.numeric(tapply(xt * e, d$t, sum))
  m <- sum(h^2)
  for (l in 1:2) m <- m + 2 * (1 - l / 3) * sum(h[-(1:l)] * h[1:(15 - l)])
  expect_equal(r$coef$se, sqrt(m) / sum(xt^2), tolerance = 1e-9)
  expect_identical(r$lag, 2L)
  # the small-sample factor of the times, T / (T - 1)
  g <- ivfe(d, "y", x = "x", fe = c("u", "t"), se = "dk", time = "t",
            lag = 2)
  expect_equal(g$coef$se, r$coef$se * sqrt(15 / 14), tolerance = 1e-12)
  expect_error(ivfe(d, "y", x = "x", fe = c("u", "t"), se = "dk"),
               "ivfe: `se = \"dk\"` needs `time`, the name of the day column.",
               fixed = TRUE)
})

test_that("ppmlcf(family = \"nb\") is the negative binomial control function", {
  # cities over three years: lead use instrumented by distance to a
  # refinery, homicides negative binomial with the population as exposure
  # (Feigenbaum and Muller 2016). Written out here: the first stage by
  # lm(), the second-stage likelihood maximised by optim() over the
  # coefficients and log theta; 1e-5 is the optimiser's agreement
  set.seed(62)
  n <- 240
  d <- data.frame(city = rep(1:80, 3), yr = rep(1:3, each = 80))
  dist <- runif(80, 0, 5)[d$city]
  u <- rnorm(80)[d$city]
  d$dist <- dist
  d$lead <- 1.5 - 0.3 * dist + 0.5 * u + rnorm(n, 0, 0.3)
  d$pop <- exp(rnorm(80, 10, 0.5))[d$city]
  mu <- d$pop * exp(-9 + 0.4 * d$lead + 0.1 * d$yr + 0.3 * u)
  d$hom <- stats::rnbinom(n, size = 3, mu = mu)
  r <- ppmlcf(d, "hom", "lead", "dist", fe = "yr", expo = "pop",
              clus = "city", family = "nb", se = "naive")
  vh <- stats::residuals(stats::lm(lead ~ dist + factor(yr), d))
  xm <- cbind(d$lead, vh, 1, d$yr == 2, d$yr == 3)
  nll <- function(q) {
    b <- q[1:5]
    th <- exp(q[6])
    m <- exp(drop(xm %*% b) + log(d$pop))
    -sum(lgamma(d$hom + th) - lgamma(th) - lgamma(d$hom + 1) +
           th * log(th / (th + m)) + d$hom * log(m / (th + m)))
  }
  o <- stats::optim(c(0, 0, -9, 0, 0, 0), nll, method = "BFGS",
                    control = list(maxit = 1000, reltol = 1e-14))
  o <- stats::optim(o$par, nll, method = "BFGS",
                    control = list(maxit = 1000, reltol = 1e-15))
  expect_identical(r$coef$term, c("lead", "vhat:lead"))
  expect_equal(r$coef$est, o$par[1:2], tolerance = 1e-5)
  expect_equal(r$theta, exp(o$par[6]), tolerance = 1e-4)
  expect_identical(r$family, "nb")
  b <- ppmlcf(d, "hom", "lead", "dist", fe = "yr", expo = "pop",
              clus = "city", family = "nb", nboot = 5, seed = 1)
  expect_identical(dim(b$draws), c(5L, 2L))
  expect_error(ppmlcf(d, "hom", "lead", "dist", fe = "yr", family = "nb",
                      how = "plug"),
               "ppmlcf: `family = \"nb\"` needs `how = \"cf\"`.", fixed = TRUE)
test_that("columns of very different scale are not taken for collinear", {
  # a share next to money: the rank check judges each column against its
  # own size; glm() and lm() with the indicators are the reference
  set.seed(1)
  d <- expand.grid(day = 1:20, area = 1:12)
  d$rate <- rnorm(nrow(d), 0.05, 0.01)
  d$gdp <- rnorm(nrow(d), 5e9, 1e9)
  d$crimes <- rpois(nrow(d), exp(1 + 2 * d$rate))
  g <- glm(crimes ~ rate + gdp + factor(area) + factor(day), poisson, d,
           control = glm.control(epsilon = 1e-12, maxit = 100))
  f <- ppmlfe(d, "crimes", c("rate", "gdp"), c("area", "day"), tol = 1e-12)
  expect_equal(f$coef$est, unname(coef(g)[c("rate", "gdp")]),
               tolerance = 1e-7)
  d$y <- 0.3 * d$rate + rnorm(nrow(d))
  h <- ivfe(d, "y", x = c("rate", "gdp"), fe = c("area", "day"))
  l <- lm(y ~ rate + gdp + factor(area) + factor(day), d)
  expect_equal(h$coef$est, unname(coef(l)[c("rate", "gdp")]),
               tolerance = 1e-8)
  # an instrument in other units: two-stage least squares does not change
  d$z <- d$rate + rnorm(nrow(d), 0, 0.01)
  d$zb <- d$z * 1e9
  a <- ivfe(d, "y", endo = "rate", inst = "z", fe = c("area", "day"))
  b <- ivfe(d, "y", endo = "rate", inst = "zb", fe = c("area", "day"))
  expect_equal(b$coef$est, a$coef$est, tolerance = 1e-9)
  expect_equal(b$coef$se, a$coef$se, tolerance = 1e-9)
})

test_that("windiv() names its sectors by their edges, formatted once", {
  # the names of sectors that are each the next one's lower edge, against
  # the edges formatted by hand
  for (w in c(90, 22.5, 360 / 7, 360 / 2^10)) {
    for (s in c(0, 10)) {
      nb <- round(360 / w)
      lo <- (s + (seq_len(nb) - 1) * w) %% 360
      a <- format(lo, trim = TRUE)
      b <- format(c(lo[-1L], lo[1L]), trim = TRUE)
      m <- windiv(c(5, 100), width = w, start = s, base = 2)
      expect_identical(colnames(m), paste0("w", a, "_", b)[-2L])
    }
  }
  expect_identical(colnames(windiv(1, width = 22.5))[1:2],
                   c("w22.5_45.0", "w45.0_67.5"))
})
