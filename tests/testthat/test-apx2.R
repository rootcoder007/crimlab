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
  expect_equal(f$coef$est, h$coef$est)
})

test_that("smkiv: raw, own and neighbour adjustment, density, counts", {
  d <- data.frame(place = rep(c("a", "b"), each = 6),
                  day = rep(as.Date("2024-07-01") + 0:5, 2),
                  plume = c(0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1),
                  pm = c(5, 6, 7, 20, 6, 6, 5, 7, 6, 30, 6, 6),
                  near = rep(c("b", "a"), each = 6),
                  dens = c("", "", "", "heavy", "light", "", "", "", "",
                           "light", "", "medium"),
                  mon = "07")
  expect_identical(smkiv(d, "plume", "place", "day")$smoke, d$plume)
  o <- smkiv(d, "plume", "place", "day", pm = "pm", adj = "own")
  expect_identical(o$smoke, c(0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0))
  expect_equal(o$spm[4L], 20 - mean(c(5, 6, 7, 6)))
  n <- smkiv(d, "plume", "place", "day", pm = "pm", adj = "nbr", nbr = "near")
  expect_identical(n$smoke[4L], 1)
  expect_identical(n$smoke[10L], 1)
  ds <- smkiv(d, "plume", "place", "day", dens = "dens")
  expect_identical(names(ds), c("smoke", "dens:heavy", "dens:light",
                                "dens:medium"))
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
  expect_equal(h$tab$est, f$coef$est)
  dd <- f$coef$est[1L] - f$coef$est[2L]
  v <- f$vcov[1, 1] + f$vcov[2, 2] - 2 * f$vcov[1, 2]
  expect_equal(h$test$stat, dd^2 / v)
  sp <- hetfx(d, "pm", "poor", how = "split", y = "crimes", x = "pm",
              fe = c("area", "day"))
  s1 <- ppmlfe(d[d$poor == "poor", ], "crimes", "pm", fe = c("area", "day"))
  expect_equal(sp$tab$est[1L], s1$coef$est)
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
  expect_equal(r$coef$est, h$coef$est)
  expect_equal(r$coef$se, h$coef$se)
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
  expect_equal(j$fit$coef$est, h$coef$est)
  q1 <- stats::quantile(d$temp, c(0.1, 0.9), names = FALSE)
  q2 <- stats::quantile(d$o3, 0.5, names = FALSE)
  b <- h$coef$est
  expect_equal(j$cond$est[2L], (q1[2] - q1[1]) * (b[1] + b[3] * q2))
  expect_equal(j$test$stat, (b[3] / h$coef$se[3])^2)
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
  expect_equal(t$tab$est[1L], a$coef$est)
  expect_equal(t$tab$padj, stats::p.adjust(t$tab$p, "holm"))
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
  expect_equal(r$bg, rep(c(16, 20, 22), each = 3))
  expect_equal(r$spm, c(0, 0, 24, 0, 0, 30, 0, 38, 0))
  expect_error(smkiv(transform(d, day = 1:9), "plume", "area", "day",
                     pm = "pm", cf = "childs", seas = "area"), "as dates")
})
