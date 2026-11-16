# The designs of R/apx2.R, on simulated data, checked against the same
# model fitted by hand with ppmlfe() or ivfe(). No number here comes from
# the papers.

test_that("exprun counts runs and runfx fits the day-in-run indicators", {
  d <- data.frame(area = c(1, 1, 1, 1, 1, 1, 2, 2), day = c(1:6, 1, 2),
                  smoke = c(0, 1, 1, 1, 1, 0, 1, 1))
  r <- exprun(d, "smoke", "area", "day")
  expect_identical(r$run, c(0L, 1L, 2L, 3L, 4L, 0L, 1L, 2L))
  expect_identical(r$run3p, c(0, 0, 0, 1, 1, 0, 0, 0))
  g <- d[-3, ]
  expect_identical(exprun(g, "smoke", "area", "day")$run[3L], 1L)
  expect_error(exprun(rbind(d, d[1, ]), "smoke", "area", "day"),
               "two rows on one day")
  set.seed(7)
  s <- expand.grid(day = 1:80, area = 1:12)
  s$smoke <- rbinom(nrow(s), 1, 0.3)
  rr <- exprun(s, "smoke", "area", "day")
  s$crimes <- rpois(nrow(s), exp(1 + 0.2 * rr$run1))
  f <- runfx(s, "crimes", "smoke", "area", "day", fe = c("area", "day"))
  h <- ppmlfe(cbind(s, rr[-1L]), "crimes", c("run1", "run2", "run3p"),
              fe = c("area", "day"))
  # the hand fit builds its own design: the iterative fits agree to their
  # convergence, not to the last digit
  expect_equal(f$coef$est, h$coef$est, tolerance = 1e-10)
})

test_that("smkiv: raw, own and neighbour adjustment, density, counts", {
  d <- data.frame(place = rep(c("a", "b"), each = 6),
                  day = rep(as.Date("2024-07-01") + 0:5, 2),
                  plume = c(0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1),
                  pm = c(5, 6, 7, 20, 6, 6, 5, 7, 6, 30, 6, 6),
                  near = rep(c("b", "a"), each = 6),
                  dens = c("", "", "", "heavy", "light", "", "", "", "",
                           "light", "", "medium"),
                  mon = "07", stringsAsFactors = FALSE)
  expect_identical(smkiv(d, "plume", "place", "day")$smoke, d$plume)
  o <- smkiv(d, "plume", "place", "day", pm = "pm", adj = "own")
  expect_identical(o$smoke, c(0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0))
  expect_equal(o$spm[4L], 20 - mean(c(5, 6, 7, 6)), tolerance = 1e-12)
  n <- smkiv(d, "plume", "place", "day", pm = "pm", adj = "nbr", nbr = "near")
  expect_identical(n$smoke[4L], 1)
  expect_identical(n$smoke[10L], 1)
  ds <- smkiv(d, "plume", "place", "day", dens = "dens")
  expect_named(ds, c("smoke", "dens:heavy", "dens:light", "dens:medium"))
  m <- smkiv(d, "plume", "place", "day", pm = "pm", adj = "own", by = "place")
  expect_identical(m$smoke, c(1, 1))
  expect_error(smkiv(d, "plume", "place", "day", adj = "own"), "needs `pm`")
})

test_that("hetfx: interactions equal a hand fit, split equals subsets", {
  set.seed(8)
  d <- expand.grid(day = 1:60, area = 1:20)
  d$poor <- ifelse(d$area <= 10, "poor", "rich")
  d$pm <- rnorm(nrow(d), 20, 5)
  d$crimes <- rpois(nrow(d), exp(0.5 + 0.04 * d$pm * (d$poor == "poor")))
  h <- hetfx(d, "pm", "poor", y = "crimes", x = "pm", fe = c("area", "day"))
  d$a <- d$pm * (d$poor == "poor")
  d$b <- d$pm * (d$poor == "rich")
  f <- ppmlfe(d, "crimes", c("a", "b"), fe = c("area", "day"))
  # the hand fit builds its own design: the iterative fits agree to their
  # convergence, not to the last digit
  expect_equal(h$tab$est, f$coef$est, tolerance = 1e-10)
  dd <- f$coef$est[1L] - f$coef$est[2L]
  v <- f$vcov[1, 1] + f$vcov[2, 2] - 2 * f$vcov[1, 2]
  expect_equal(h$test$stat, dd^2 / v, tolerance = 1e-10)
  sp <- hetfx(d, "pm", "poor", how = "split", y = "crimes", x = "pm",
              fe = c("area", "day"))
  s1 <- ppmlfe(d[d$poor == "poor", ], "crimes", "pm", fe = c("area", "day"))
  expect_equal(sp$tab$est[1L], s1$coef$est, tolerance = 1e-10)
  expect_identical(sp$test$df, 1L)
})

test_that("ccross is Poisson with place-year-month-weekday strata", {
  set.seed(9)
  d <- expand.grid(day = as.Date("2023-01-01") + 0:364, area = 1:4)
  d$pm <- rnorm(nrow(d), 15, 4)
  d$crimes <- rpois(nrow(d), exp(1 + 0.02 * d$pm))
  r <- ccross(d, "crimes", "area", "day", x = "pm")
  d$s <- paste(d$area, format(d$day, "%Y-%m"), as.POSIXlt(d$day)$wday)
  d$ym <- format(d$day, "%Y-%m")
  h <- ppmlfe(d, "crimes", "pm", fe = "s", clus = "ym")
  # the hand fit builds its own design: the iterative fits agree to their
  # convergence, not to the last digit
  expect_equal(r$coef$est, h$coef$est, tolerance = 1e-10)
  expect_equal(r$coef$se, h$coef$se, tolerance = 1e-10)
  m <- ccross(d, "crimes", "area", "day", ma = "pm", x = "pm")
  expect_identical(m$n, nrow(d) - 4L)
  expect_error(ccross(transform(d, day = 1), "crimes", "area", "day"),
               "must hold dates")
})

test_that("jntx: the linear form equals a hand fit with the product", {
  set.seed(10)
  d <- expand.grid(day = 1:150, city = 1:10)
  d$temp <- rnorm(nrow(d), 25, 4)
  d$o3 <- 0.6 * d$temp + rnorm(nrow(d), 30, 6)
  d$crimes <- rpois(nrow(d), exp(1 + 0.02 * d$temp))
  j <- jntx(d, "crimes", "temp", "o3", form = "linear", fe = "city")
  d$tx <- d$temp * d$o3
  h <- ppmlfe(d, "crimes", c("temp", "o3", "tx"), fe = "city")
  # the hand fit builds its own design: the iterative fits agree to their
  # convergence, not to the last digit
  expect_equal(j$fit$coef$est, h$coef$est, tolerance = 1e-10)
  q1 <- stats::quantile(d$temp, c(0.1, 0.9), names = FALSE)
  q2 <- stats::quantile(d$o3, 0.5, names = FALSE)
  b <- h$coef$est
  expect_equal(j$cond$est[2L], (q1[2] - q1[1]) * (b[1] + b[3] * q2),
               tolerance = 1e-10)
  expect_equal(j$test$stat, (b[3] / h$coef$se[3])^2, tolerance = 1e-10)
  s <- jntx(d, "crimes", "temp", "o3", fe = "city")
  expect_identical(ncol(s$fit$vcov), 4L + 3L + 12L)
  expect_identical(nrow(s$surf), 225L)
})

test_that("dosebin: quintiles, the reference and the tests", {
  set.seed(11)
  d <- expand.grid(day = 1:100, area = 1:10)
  d$aqi <- rgamma(nrow(d), 4, 0.1)
  d$crimes <- rpois(nrow(d), exp(1 + 0.3 * (d$aqi > 50)))
  r <- dosebin(d, "crimes", "aqi", fe = c("area", "day"))
  expect_identical(nrow(r$tab), 5L)
  expect_identical(r$tab$est[1L], 0)
  expect_identical(r$test$test, c("all", "linear"))
  expect_identical(r$test$df, c(4L, 3L))
  b <- dosebin(d, "crimes", "aqi", brk = c(0, 50, 100, Inf),
               fe = c("area", "day"))
  expect_identical(nrow(b$tab), 3L)
})

test_that("typefx: one fit per outcome and the Holm correction", {
  set.seed(12)
  d <- expand.grid(day = 1:60, area = 1:12)
  d$pm <- rnorm(nrow(d), 20, 5)
  d$assault <- rpois(nrow(d), exp(0.5 + 0.03 * d$pm))
  d$theft <- rpois(nrow(d), 2)
  t <- typefx(d, c("assault", "theft"), "pm", x = "pm", fe = c("area", "day"))
  a <- ppmlfe(d, "assault", "pm", fe = c("area", "day"))
  expect_equal(t$tab$est[1L], a$coef$est, tolerance = 1e-12)
  expect_identical(t$tab$padj, stats::p.adjust(t$tab$p, "holm"))
  expect_error(typefx(d, "assault", "pm", adj = "romano"), "should be one of")
})

test_that("smkiv cf = childs is the three-year month median (Childs 2022)", {
  # one place, three January days in each of 2019-2021; the background
  # of a year is the median of the smoke-free days of that year and the
  # years either side (Childs et al. 2022, equation 2), by hand:
  # 2019: 10, 12, 20, 22 -> 16; 2020: 10, 12, 20, 22, 30 -> 20;
  # 2021: 20, 22, 30 -> 22
  d <- data.frame(
    area = 1,
    day = as.Date(c("2019-01-01", "2019-01-02", "2019-01-03",
                    "2020-01-01", "2020-01-02", "2020-01-03",
                    "2021-01-01", "2021-01-02", "2021-01-03")),
    pm = c(10, 12, 40, 20, 22, 50, 30, 60, 15),
    plume = c(0, 0, 1, 0, 0, 1, 0, 1, 1)
  )
  r <- smkiv(d, "plume", "area", "day", pm = "pm", cf = "childs")
  expect_equal(r$bg, rep(c(16, 20, 22), each = 3), tolerance = 1e-12)
  expect_equal(r$spm, c(0, 0, 24, 0, 0, 30, 0, 38, 0), tolerance = 1e-12)
  expect_error(smkiv(transform(d, day = 1:9), "plume", "area", "day",
                     pm = "pm", cf = "childs", seas = "area"), "as dates")
})

test_that("ccross: single-day lags and spline controls (Han et al. 2026)", {
  set.seed(31)
  d <- expand.grid(day = as.Date("2023-01-01") + 0:180, area = 1:3)
  d$pm <- rnorm(nrow(d), 15, 4)
  d$temp <- rnorm(nrow(d), 10, 5)
  d$crimes <- rpois(nrow(d), exp(1 + 0.02 * d$pm))
  f <- ccross(d, "crimes", "area", "day", lagv = "pm", lag = 1, x = "pm")
  # the same fit by hand: pm of the day before, the first day dropped
  e <- d[order(d$area, d$day), ]
  e$pm <- stats::ave(e$pm, e$area, FUN = function(v) c(NA, v[-length(v)]))
  e <- e[!is.na(e$pm), ]
  e$.strat <- paste(e$area, format(e$day, "%Y-%m"), as.POSIXlt(e$day)$wday)
  e$.ym <- format(e$day, "%Y-%m")
  g <- ppmlfe(e, "crimes", "pm", fe = ".strat", clus = ".ym")
  expect_equal(f$coef$est, g$coef$est, tolerance = 1e-10)
  h <- ccross(d, "crimes", "area", "day", x = "pm", ns = c(temp = 4))
  expect_identical(h$coef$term, c("pm", paste0("temp:ns", 1:4)))
  expect_error(ccross(d, "crimes", "area", "day", x = "pm", ns = 3),
               "named by column")
})

# The refusals and printed summaries of the designs, and the branches
# their usual estimators do not reach: an estimator without vcov (the
# draws of ppmlcf()), instruments interacted by level, a one-degree
# spline.

test_that("smkiv refuses a bad plume, density, season and background", {
  d <- data.frame(place = rep(c("a", "b"), each = 6),
                  day = rep(as.Date("2024-07-01") + 0:5, 2),
                  plume = c(0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1),
                  pm = c(5, 6, 7, 20, 6, 6, 5, 7, 6, 30, 6, 6))
  d2 <- d
  d2$plume[1] <- 2
  dn <- d
  dn$day <- seq_len(nrow(dn))
  ds <- d
  ds$plume[1:5] <- 1
  errs(
    "smkiv",
    list(
      list("`dens` must name a column of `data`", data = d, plume = "plume",
           unit = "place", day = "day", dens = 1),
      list("column `plume` must hold 0 or 1", data = d2, plume = "plume",
           unit = "place", day = "day"),
      list("give `seas`, or `day` as dates for the meteorological seasons",
           data = dn, plume = "plume", unit = "place", day = "day",
           pm = "pm"),
      list("a place and season has fewer than two days without a plume",
           data = ds, plume = "plume", unit = "place", day = "day",
           pm = "pm", adj = "own"),
      list("`adj = \"nbr\"` needs `nbr`, the neighbour of each row",
           data = d, plume = "plume", unit = "place", day = "day", pm = "pm",
           adj = "nbr")
    )
  )
})

test_that("smkiv cf = childs: a month without a smoke-free day is NA", {
  d <- data.frame(
    area = 1,
    day = as.Date(c("2019-01-01", "2019-01-02", "2019-01-03",
                    "2020-01-01", "2020-01-02", "2020-01-03",
                    "2021-01-01", "2021-01-02", "2021-01-03", "2020-02-01")),
    pm = c(10, 12, 40, 20, 22, 50, 30, 60, 15, 25),
    plume = c(0, 0, 1, 0, 0, 1, 0, 1, 1, 1)
  )
  expect_warning(
    r <- smkiv(d, "plume", "area", "day", pm = "pm", cf = "childs"),
    "1 row\\(s\\) have no smoke-free days in their place, month and three"
  )
  # the January rows as in the hand computation above; February has none
  expect_equal(r$bg[1:9], rep(c(16, 20, 22), each = 3), tolerance = 1e-12)
  expect_true(is.na(r$bg[10L]))
})

test_that("the estimator checks and the one-level moderator", {
  set.seed(14)
  d <- expand.grid(day = 1:20, area = 1:4)
  d$day2 <- as.Date("2024-01-01") + d$day
  d$pm <- rnorm(nrow(d), 20, 5)
  d$crimes <- rpois(nrow(d), 3)
  d$one <- "x"
  d$dc <- as.character(d$day)
  d$tc <- "a"
  ee <- "`est` must be an estimator function, such as ppmlfe"
  errs("runfx", list(list(ee, data = d, y = "crimes", var = "pm",
                          unit = "area", day = "day", est = 1)))
  errs("hetfx", list(
    list(ee, data = d, var = "pm", mod = "area", est = 1),
    list("column `one` has one level: nothing to compare", data = d,
         var = "pm", mod = "one")
  ))
  errs("ccross", list(
    list(ee, data = d, y = "crimes", unit = "area", day = "day2", est = 1),
    list("column `tc` must be numeric", data = d, y = "crimes",
         unit = "area", day = "day2", x = "pm", ns = c(tc = 3))
  ))
  errs("jntx", list(
    list(ee, data = d, y = "crimes", x1 = "pm", x2 = "day", est = 1),
    list("`p1` must be two quantiles", data = d, y = "crimes", x1 = "pm",
         x2 = "day", p1 = 0.5)
  ))
  errs("dosebin", list(list(ee, data = d, y = "crimes", var = "pm",
                            est = 1)))
  errs("typefx", list(
    list("`term` must name one coefficient", data = d, ys = "crimes",
         term = c("a", "b")),
    list(ee, data = d, ys = "crimes", term = "pm", est = 1),
    list("the estimator returned no coefficient for `nope`", data = d,
         ys = "crimes", term = "nope", x = "pm", fe = "area", se = "hc1")
  ))
  errs("exprun", list(list("column `dc` must hold numbers or dates",
                           data = d, var = "pm", unit = "area",
                           day = "dc")))
})

test_that("the printed summaries of the designs", {
  set.seed(15)
  d <- expand.grid(day = 1:40, area = 1:10)
  d$smoke <- rbinom(nrow(d), 1, 0.3)
  d$pm <- rnorm(nrow(d), 20, 5)
  d$o3 <- rnorm(nrow(d), 30, 5)
  d$g <- ifelse(d$area <= 5, "a", "b")
  d$crimes <- rpois(nrow(d), exp(0.5 + 0.02 * d$pm))
  fe <- c("area", "day")
  expect_output(print(runfx(d, "crimes", "smoke", "area", "day", fe = fe)),
                "runfx: the effect by day of a run of exposed days")
  h <- hetfx(d, "pm", "g", y = "crimes", x = "pm", fe = fe)
  expect_output(print(h), "hetfx: the effect by level \\(interactions\\)")
  expect_output(print(h), "Wald test of equal effects: .* on 1 df")
  j <- jntx(d, "crimes", "pm", "o3", form = "linear", fe = fe)
  expect_output(print(j), "jntx: the effect of x1 \\(linear surface\\)")
  expect_output(print(j), "Wald test of no interaction: .* on 1 df")
  b <- dosebin(d, "crimes", "pm", nq = 3, fe = fe)
  expect_output(print(b), "dosebin: the effect of each bin against the")
  t <- typefx(d, "crimes", "pm", x = "pm", fe = fe)
  expect_output(print(t), "typefx: one fit per outcome, p-values adjusted by")
})

test_that("hetfx with ppmlcf: instruments by level, covariance from draws", {
  set.seed(16)
  d <- expand.grid(day = 1:30, area = 1:12)
  d$z <- rnorm(nrow(d))
  d$x <- 0.8 * d$z + rnorm(nrow(d))
  d$g <- ifelse(d$area <= 6, "a", "b")
  d$crimes <- rpois(nrow(d), exp(0.5 + 0.2 * d$x))
  fe <- c("area", "day")
  h <- hetfx(d, "x", "g", est = ppmlcf, y = "crimes", endo = "x",
             inst = "z", fe = fe, nboot = 5, seed = 1)
  # the same fit by hand: the exposure and the instrument times each level
  e <- d
  for (l in c("a", "b")) {
    e[[paste0("x:", l)]] <- e$x * (e$g == l)
    e[[paste0("z:", l)]] <- e$z * (e$g == l)
  }
  p <- ppmlcf(e, "crimes", endo = c("x:a", "x:b"), inst = c("z:a", "z:b"),
              fe = fe, nboot = 5, seed = 1)
  expect_equal(h$tab$est, p$coef$est[1:2], tolerance = 1e-12)
  v <- stats::cov(p$draws, use = "complete.obs")[1:2, 1:2]
  w <- (p$coef$est[1] - p$coef$est[2])^2 / (v[1, 1] + v[2, 2] - 2 * v[1, 2])
  expect_equal(h$test$stat, w, tolerance = 1e-10)
})

test_that("ccross with one degree of freedom is the linear term", {
  set.seed(17)
  d <- expand.grid(day = as.Date("2023-01-01") + 0:120, area = 1:3)
  d$pm <- rnorm(nrow(d), 15, 4)
  d$temp <- rnorm(nrow(d), 10, 5)
  d$crimes <- rpois(nrow(d), exp(1 + 0.02 * d$pm))
  a <- ccross(d, "crimes", "area", "day", x = "pm", ns = c(temp = 1))
  b <- ccross(d, "crimes", "area", "day", x = c("pm", "temp"))
  expect_identical(a$coef$term, c("pm", "temp:ns1"))
  expect_equal(a$coef$est, b$coef$est, tolerance = 1e-12)
})

test_that("jntx and dosebin refuse exposures that give no bins", {
  set.seed(18)
  d <- expand.grid(day = 1:30, area = 1:4)
  d$two <- rep(c(1, 2), length.out = nrow(d))
  d$aqi <- rgamma(nrow(d), 4, 0.1)
  d$crimes <- rpois(nrow(d), 3)
  expect_error(jntx(d, "crimes", "two", "aqi", fe = "area"),
               "an exposure has too few distinct values for its knots")
  errs("dosebin", list(
    list("`brk` must be three or more increasing numbers", data = d,
         y = "crimes", var = "aqi", brk = c(0, 1)),
    list("some values of `aqi` fall outside `brk`", data = d, y = "crimes",
         var = "aqi", brk = c(0, 1, 2)),
    list("`two` gives fewer than two bins", data = d, y = "crimes",
         var = "two", nq = 5)
  ))
})

test_that("dosebin refuses a constant exposure: fewer than two bins", {
  # With one distinct value the quantile breaks are a single number,
  # which cut() would read as a number of intervals (7 here): the
  # refusal has to come before cut().
  set.seed(19)
  d <- expand.grid(day = 1:30, area = 1:4)
  d$k <- 7
  d$crimes <- rpois(nrow(d), 3)
  expect_error(dosebin(d, "crimes", "k", fe = "area", se = "hc1"),
               "dosebin: `k` gives fewer than two bins")
})

test_that("counts that size the designs are refused past their bounds", {
  d <- data.frame(area = c(1, 1, 1, 1, 1, 1, 2, 2), day = c(1:6, 1, 2),
                  smoke = c(0, 1, 1, 1, 1, 0, 1, 1))
  expect_error(exprun(d, "smoke", "area", "day", cap = 1e9),
               "exprun: `cap` must be at most 67108864.", fixed = TRUE)
  set.seed(10)
  g <- expand.grid(day = 1:20, city = 1:3)
  g$temp <- rnorm(nrow(g), 25, 4)
  g$o3 <- rnorm(nrow(g), 45, 6)
  g$crimes <- rpois(nrow(g), 3)
  expect_error(jntx(g, "crimes", "temp", "o3", ngrid = 1e5),
               "jntx: `ngrid` times `ngrid` must be at most 67108864.",
               fixed = TRUE)
  g$aqi <- rgamma(nrow(g), 4, 0.1)
  expect_error(dosebin(g, "crimes", "aqi", nq = 1e9),
               "dosebin: `nq` must be at most 67108864.", fixed = TRUE)
})
