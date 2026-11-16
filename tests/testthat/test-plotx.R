test_that("Lorenz and funnel plots", {
  x <- c(0, 0, 1, 1, 2, 3, 8, 25)
  d <- drw(plotlrnz(x))
  r <- cconc(x, nsim = 0)
  expect_identical(d$share, r$lorenz$share)
  expect_identical(attr(d, "gini"), r$stat[["gini"]])
  hand <- data.frame(pop = c(0.5, 1), share = c(0.2, 1))
  h <- drw(plotlrnz(hand, main = "given"))
  expect_identical(h$share, c(0, 0.2, 1))
  expect_true(is.na(attr(h, "gini")))
  expect_identical(nrow(drw(plotlrnz(r$lorenz))), nrow(r$lorenz))
  expect_error(plotlrnz(data.frame(a = 1)), "^plotlrnz: a frame `x` needs")
  expect_error(plotlrnz(data.frame(pop = c(0.5, 1), share = c(0.6, 0.2))),
               "^plotlrnz: `x` must be a Lorenz curve")
  expect_error(plotlrnz(c(-1, 2)), "^plotlrnz: `x` must be finite values")
  f <- drw(plotfnl(obs = c(4, 20, 50, 2), expd = c(5, 22, 40, 20),
                   area = c("a", "b", "c", "d")))
  expect_identical(f$ratio, c(4, 20, 50, 2) / c(5, 22, 40, 20))
  e <- c(5, 22, 40, 20)
  lo <- stats::qpois(0.001, e) / e
  hi <- stats::qpois(0.999, e) / e
  expect_identical(f$outside, f$ratio < lo | f$ratio > hi)
  expect_true(f$outside[4L])
  expect_true(all(is.na(drw(plotfnl(c(4, 20), c(5, 22), col = "red"))$area)))
  g <- drw(plotfnl(c(4, NA, 9), c(5, 10, NA), area = c("p", "q", "r")))
  expect_identical(g$area, "p")
  expect_identical(c(g$observed, g$expected), c(4, 5))
  expect_error(plotfnl("a", 1), "^plotfnl: `obs` must be numbers")
  expect_error(plotfnl(c(1, 2), c(1, 2, 3)), "^plotfnl: `obs` and `expd`")
  expect_error(plotfnl(1, 1, area = c("a", "b")), "one label per area")
  expect_error(plotfnl(1, 0), "^plotfnl: the counts must be finite")
  expect_error(plotfnl(1, 1, lvls = 1), "^plotfnl: `lvls` must be")
  expect_error(plotfnl(1, 1, target = 0), "^plotfnl: `target` must be")
})

test_that("distributions by group", {
  set.seed(1)
  x <- c(stats::rnorm(100), stats::rnorm(100, 2))
  g <- rep(c("a", "b"), each = 100)
  r <- drw(plotrdg(x, g))
  expect_identical(unique(r$grp), c("a", "b"))
  expect_equal(sum(r$y[r$grp == "a"]) * diff(r$x[1:2]), 1, tolerance = 0.02)
  expect_identical(r$y[r$grp == "a"], .dens(x[1:100])$y)
  # a refusal of the kernel comes in the words of stats::density(), and the
  # warning of bw.ucv() when the minimum is at an end of its range
  expect_error(drw(plotrdg(c(1, 1, 1, 1), bw = "ucv")),
               "^data are constant in bandwidth calculation$")
  dv <- c(2.1, 3.4, 1.9, 5.2, 4.4, 3.3, 2.8, 6.1, 3.9, 4.7, 1.2, 3.0)
  expect_warning(drw(plotrdg(dv, bw = "ucv")),
                 "minimum occurred at one end of the range", fixed = TRUE)
  expect_s3_class(drw(plotrdg(c(x, 5), c(g, "lone"), main = "Ridges")),
                  "data.frame")
  e <- drw(plotcdfs(x, g))
  expect_identical(e$y[e$grp == "a"][100L], 1)
  expect_identical(e$x[e$grp == "b"], sort(x[101:200]))
  v <- drw(plotvln(c(x, 9), c(g, "c")))
  expect_identical(attr(v, "summary")$median[1L], stats::median(x[1:100]))
  h <- drw(plothex(x, x + stats::rnorm(200), bins = 10))
  expect_identical(sum(h$count), 200L)
  expect_identical(sum(drw(plothex(c(1, 1, 1), c(2, 2, 2)))$count), 3L)
  expect_identical(sum(drw(plothex(x, x, col = "black", bins = 5))$count),
                   200L)
  w <- stats::rexp(200)
  k <- drw(plotgain(w, g))
  expect_identical(attr(k, "gini")[["a"]],
                   cconc(w[1:100], nsim = 0)$stat[["gini"]])
  expect_identical(utils::tail(k$share[k$grp == "a"], 1), 1)
  expect_error(plotrdg("a"), "^plotrdg: `x` must be a numeric vector")
  expect_error(plotrdg(1:3, 1:2), "^plotrdg: `grp` must have one entry")
  expect_error(plotrdg(1:2, c("a", "b")), "^plotrdg: no group has two")
  expect_error(plotrdg(1:3, bw = -1), "^plotrdg: `bw` must be one positive")
  expect_error(plotvln(1:3, bw = "x"), "^plotvln: `bw` must be")
  expect_error(plotcdfs(c(1, NA)), "^plotcdfs: `x` needs at least two")
  expect_error(plothex(1:3, "a"), "^plothex: `y` must be numbers")
  expect_error(plothex(c(1, NA, 3), c(NA, 2, 4)), "at least two points")
  expect_error(plothex(1:3, 1:3, bins = 1), "^plothex: `bins` must be")
  expect_error(plotgain(c(-1, 2, 3)), "^plotgain: `x` must not be negative")
})

test_that("series over time", {
  days <- seq(as.Date("2023-12-25"), as.Date("2024-03-01"), by = "day")
  cal <- drw(plotcal(days, c(NA, seq_along(days)[-1L])))
  expect_identical(nrow(cal), length(days))
  expect_true(all(cal$wday %in% 1:7))
  expect_identical(cal$y[-1L], as.numeric(seq_along(days)[-1L]))
  expect_s3_class(drw(plotcal(as.character(days[1:3]), c(1, 1, 1),
                              col = "black")), "data.frame")
  expect_error(plotcal("not a date", 1), "^plotcal: `x` must be dates")
  d <- expand.grid(year = 2015:2020, area = c("n", "s", "e"))
  d$v <- c(1:6, 6:1, rep(3, 6))
  hz <- drw(plothrz(d$v, d$year, d$area, bands = 2))
  expect_true(all(hz$band %in% 0:2))
  expect_identical(hz$band[hz$grp == "e"], rep(0, 6))
  expect_s3_class(drw(plothrz(rep(1, 3), 1:3)), "data.frame")
  expect_error(plothrz(1:3, 1:3, bands = 9), "^plothrz: `bands` must be")
  sm <- drw(plotmult(d$v, d$year, d$area, main = "Small"))
  expect_identical(nrow(sm), 18L)
  expect_identical(sm$y, as.numeric(d$v))
  q <- drw(plotmult(1:4, factor(c("q1", "q2", "q3", "q4"))))
  expect_identical(q$y, as.numeric(1:4))
  expect_s3_class(drw(plotmult(1:4, c("q1", "q2", "q3", "q4"))),
                  "data.frame")
  b <- drw(plotbump(d$v, d$year, d$area))
  expect_identical(b$rank[b$grp == "n" & b$prd == 2015], 3)
  expect_identical(b$rank[b$grp == "s" & b$prd == 2015], 1)
  asc <- drw(plotbump(d$v, d$year, d$area, desc = FALSE))
  expect_identical(asc$rank[asc$grp == "n" & asc$prd == 2015], 1)
  expect_error(plotbump(1, 1, 1, desc = NA), "^plotbump: `desc` must be")
  expect_error(plotmult("a", 1, 1), "^plotmult: `y` must be a numeric")
  expect_error(plotmult(1:3, 1:2), "^plotmult: `prd` and `grp` must")
  expect_error(plotmult(1:2, c(NA, NA)), "^plotmult: no value has both")
})

test_that("event studies, alone and in a grid", {
  d <- did()
  e <- evtst(d$rate, d$town, d$year, d$first, leads = 3, lags = 3)
  p <- drw(plot(e))
  expect_identical(p$est, e$est$est)
  expect_identical(p$lo, e$est$lo)
  s <- stgdid(d$rate, d$town, d$year, d$first, reps = 19, seed = 1)
  q <- drw(plot(s, main = "Staggered"))
  expect_identical(q$t, s$dynamic$evtm)
  expect_identical(q$hi, s$dynamic$hi)
  g <- drw(plotevgrid(list(Event = e, Staggered = s), main = "Two"))
  expect_identical(unique(g$panel), c("Event", "Staggered"))
  expect_identical(g$est[g$panel == "Event"], e$est$est)
  expect_s3_class(drw(plotevgrid(list(One = e))), "data.frame")
  expect_error(plotevgrid(list(e)), "^plotevgrid: `x` must be a named list")
})

test_that("colours past an integer and a matrix of bins are refused quietly", {
  x <- c(1, 2, 3, 4)
  for (cl in list(1e300, Inf, -1, "#zz", sum)) {
    expect_error(expect_warning(plothex(x, x, col = cl), NA),
                 "^plothex: `col` must be colours\\.$")
  }
  expect_error(expect_warning(plothex(x, x, bins = matrix(5)), NA),
               "^plothex: `bins` must be one whole number from 2 to 500\\.$")
})

test_that("a title that is not valid UTF-8 is refused before the device", {
  bad <- rawToChar(as.raw(c(0xff, 0xfe, 0x41)))
  expect_error(plothex(c(1, 2), c(1, 2), main = bad),
               "^plothex: `main` holds text that is not valid")
})

test_that("the compiled density agrees with stats::density()", {
  # the kernel is the density() of R 4.3 and later (old.coords = FALSE):
  # before that R placed the grid, and worded some refusals, differently
  skip_if(getRversion() < "4.3.0", "stats::density() differs before R 4.3")
  dens <- .dens
  # the kernel's own radix-2 FFT and R's mixed-radix one round differently:
  # the two agree to about 1e-14 relative, so 1e-12 leaves room. ucv and
  # bcv minimise a criterion with optimize()'s steps: rounding in the
  # criterion (fused multiply-adds or not) moves the minimum found by about
  # 1e-12 relative on 1000 values, so 1e-10 for those
  near <- function(r, g, tol = 1e-12) {
    expect_equal(g$bw, r$bw, tolerance = tol)
    expect_equal(g$x, r$x, tolerance = tol)
    expect_equal(g$y, r$y, tolerance = tol)
  }
  set.seed(1)
  sets <- list(faithful$eruptions, as.numeric(precip),
               rep(c(1, 2, 2, 3, 5), 20), c(0.5, 1.5, 4),
               c(stats::rnorm(700), stats::rexp(300)))
  rules <- c("nrd0", "nrd", "ucv", "bcv", "sj", "sj-ste", "SJ-dpi")
  for (x in sets) {
    for (r in rules) {
      near(suppressWarnings(stats::density(x, bw = r)), dens(x, r),
           if (r %in% c("ucv", "bcv")) 1e-10 else 1e-12)
    }
  }
  # bw.ucv() and bw.bcv() warn when the minimum is at one end of the range
  expect_warning(stats::bw.ucv(c(0.5, 1.5, 4)), "one end of the range")
  expect_true(dens(c(0.5, 1.5, 4), "ucv")$warn)
  expect_false(dens(faithful$eruptions, "ucv")$warn)
  x <- faithful$eruptions
  ks <- c("gaussian", "epanechnikov", "rectangular", "triangular",
          "biweight", "cosine", "optcosine")
  for (k in ks) {
    for (n in c(100L, 512L, 1000L)) {
      args <- list(x, kernel = k, n = n, adjust = 0.7, cut = 2)
      near(do.call(stats::density, args), do.call(dens, args))
    }
  }
  # past 512 points the grid is the next power of 2
  expect_identical(dens(x, n = 1000L)$n, 1024)
  expect_identical(dens(x, n = 100L)$n, 512)
  w <- seq_along(x) / sum(seq_along(x))
  near(stats::density(x, bw = 0.3, weights = w, from = 1, to = 6),
       dens(x, 0.3, weights = w, from = 1, to = 6))
  # a rule ignores the weights, with R's warning
  near(suppressWarnings(stats::density(x, bw = "ucv", weights = w)),
       dens(x, "ucv", weights = w), 1e-10)
  # values that are not finite are dropped, and the mass is scaled by the
  # share of those kept
  xi <- c(x, Inf, -Inf)
  near(stats::density(xi, bw = 0.3), dens(xi, 0.3))
  wi <- c(w, 0.1, 0.2)
  near(suppressWarnings(stats::density(xi, bw = 0.3, weights = wi)),
       dens(xi, 0.3, weights = wi))
  # from far above to: the grid runs down, and R's approx() sorts it
  near(stats::density(x, from = 5, to = 1), dens(x, from = 5, to = 1))
  near(stats::density(x, n = 1L), dens(x, n = 1L))
  # values within one step of the grid below its first point or above its
  # last share their weight with that point alone, as R's BinDist() does
  g <- seq(0, 10, by = 0.001)
  near(stats::density(g, bw = 0.1, from = 5, to = 6),
       dens(g, 0.1, from = 5, to = 6))
  # constant data and a zero IQR fall back as bw.nrd0() does
  for (v in list(rep(3, 5), rep(0, 5), c(rep(1, 10), 2))) {
    near(stats::density(v), dens(v))
  }
  near(stats::density(1e6 + x, bw = 1e-3), dens(1e6 + x, 1e-3))
})

test_that("the compiled density refuses what stats::density() refuses", {
  # the kernel is the density() of R 4.3 and later (old.coords = FALSE):
  # before that R placed the grid, and worded some refusals, differently
  skip_if(getRversion() < "4.3.0", "stats::density() differs before R 4.3")
  dens <- .dens
  x <- faithful$eruptions
  same <- function(r, g) {
    m <- tryCatch(r, error = conditionMessage)
    expect_error(g, paste0("dens: ", m), fixed = TRUE)
  }
  same(stats::density(1, bw = "nrd0"), dens(1, "nrd0"))
  same(stats::density(rep(2, 5), bw = "ucv"), dens(rep(2, 5), "ucv"))
  same(stats::density(rep(2, 5), bw = "sj"), dens(rep(2, 5), "sj"))
  same(stats::density(rep(2, 600), bw = "sj"), dens(rep(2, 600), "sj"))
  # past 500 values R bins constant data into no bin and fails later
  same(stats::density(rep(2, 600), bw = "ucv"), dens(rep(2, 600), "ucv"))
  same(stats::density(rep(2, 600), bw = "bcv"), dens(rep(2, 600), "bcv"))
  same(stats::density(c(rep(0, 10), 1, 2), bw = "sj"),
       dens(c(rep(0, 10), 1, 2), "sj"))
  same(stats::density(x, bw = Inf), dens(x, Inf))
  same(stats::density(x, bw = -1), dens(x, -1))
  same(stats::density(x, bw = 0.3, adjust = 0), dens(x, 0.3, adjust = 0))
  same(stats::density(c(rep(1, 10), 2), bw = "nrd"),
       dens(c(rep(1, 10), 2), "nrd"))
  same(stats::density(c(-1.7e308, 0), bw = 1e308), dens(c(-1.7e308, 0), 1e308))
  same(stats::density(c(0, 1.7e308), bw = 1e308), dens(c(0, 1.7e308), 1e308))
  same(stats::density(x, bw = 0.3, from = Inf), dens(x, 0.3, from = Inf))
  same(stats::density(x, bw = 0.3, to = -Inf), dens(x, 0.3, to = -Inf))
  # what R refuses in its own words, the kernel refuses as out of range
  arg <- "dens: a size or an index is out of range"
  expect_error(stats::density(x, kernel = "box"), "should be one of")
  expect_error(dens(x, kernel = "box"), arg, fixed = TRUE)
  expect_error(stats::density(x, bw = "box"), "unknown bandwidth rule")
  expect_error(dens(x, "box"), arg, fixed = TRUE)
  expect_error(dens(x, n = 0L), arg, fixed = TRUE)
  expect_error(dens(numeric(0), 1), arg, fixed = TRUE)
  expect_error(dens(c(Inf, NaN), 1), arg, fixed = TRUE)
  expect_error(stats::density(x, weights = c(-1, rep(1, 271))),
               "'weights' must not be negative")
  expect_error(dens(x, 0.3, weights = c(-1, rep(1, 271))), arg, fixed = TRUE)
  expect_error(stats::density(x, weights = c(Inf, rep(1, 271))),
               "'weights' must all be finite")
  expect_error(dens(x, 0.3, weights = c(Inf, rep(1, 271))), arg, fixed = TRUE)
  expect_error(stats::density(x, weights = 1), "unequal length")
  expect_error(dens(x, 0.3, weights = 1), arg, fixed = TRUE)
})
