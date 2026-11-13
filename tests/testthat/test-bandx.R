# bandx(), pbbnd() and yoypx(): band labels read by hand, values and
# expansions by hand, the envelope formulas, and binom.test(), p.adjust()
# and qbeta() for the year-over-year tests.

test_that("bandx() reads the usual band labels", {
  lab <- c(
    "1-5",
    "6 to 15",
    "16 through 30",
    "31 – 60",
    "Over 60",
    "60+",
    "60 and over",
    "more than 1,000",
    "10 or less",
    "Under 10",
    "<10",
    "less than 2.5",
    "7",
    "about",
    "5-1",
    "1 to over 5",
    "1 2"
  )
  b <- bandx(lab, fac = 2, floor = 0)$bands
  expect_identical(
    b$lower,
    c(1, 6, 16, 31, 61, 60, 60, 1001, NA, NA, NA, NA, 7, NA, NA, NA, NA)
  )
  expect_identical(
    b$upper,
    c(5, 15, 30, 60, NA, NA, NA, NA, 10, 9, 9, 2.5, 7, NA, NA, NA, NA)
  )
  expect_identical(
    b$openup,
    c(rep(FALSE, 4), TRUE, TRUE, TRUE, TRUE, rep(FALSE, 9))
  )
  expect_identical(b$openlo[9:12], rep(TRUE, 4))
  # an open top band ends at twice its start; an open bottom at the floor
  expect_identical(b$value[5], (61 + 122) / 2)
  expect_identical(b$value[10], 4.5)
  expect_identical(b$value[14], NA_real_)
  # 60+ from 61 without clsup, over 60 from 60 without intsca
  expect_identical(bandx("60+", clsup = FALSE)$bands$lower, 61)
  expect_identical(bandx("over 60", intsca = FALSE)$bands$lower, 60)
  expect_identical(bandx(factor("2-4"))$bands$value, 3)
  expect_identical(bandx(c(NA, ""))$bands$lower, c(NA_real_, NA_real_))
})

test_that("bandx() places values, expands counts and sweeps the cap", {
  b <- c("1-5", "6-15", "16-30", "Over 30")
  n <- c(10, 6, 3, 1)
  r <- bandx(b, n, rule = "lower")
  expect_identical(r$vals, rep(c(1, 6, 16, 31), n))
  g <- bandx(b, rule = "geometric", cap = 50)$bands
  expect_equal(
    g$value,
    sqrt(c(1 * 5, 6 * 15, 16 * 30, 31 * 50)),
    tolerance = 1e-15
  )
  expect_identical(
    bandx(b, rule = "upper", cap = 50)$bands$value,
    c(5, 15, 30, 50)
  )
  s <- bandx(b, n, stat = mean, caps = c(40, 60, 100))
  mids <- function(cp) c(3, 10.5, 23, (31 + cp) / 2)
  want <- vapply(c(40, 60, 100), function(cp) mean(rep(mids(cp), n)), 0)
  expect_equal(s$sens$value, want, tolerance = 1e-15)
  expect_equal(s$span, max(want) - min(want), tolerance = 1e-14)
  expect_equal(s$relspan, s$span / median(want), tolerance = 1e-14)
  d <- bandx(b, n, stat = median)
  expect_identical(nrow(d$sens), 12L)
  expect_identical(d$sens$cap[1], round(31 * 1.05, 6))
  z <- bandx(c("1-2", "Over 2"), c(1, 1), stat = function(v) 0, caps = c(5, 6))
  expect_identical(c(z$span, z$relspan), c(0, NA))
  expect_output(print(s), "over the caps tried")
  expect_output(print(r), "4 bands")
  # unreadable bands: an error, or left out
  expect_identical(
    bandx(c("1-2", "lots"), c(2, 5), drop = TRUE)$vals,
    c(1.5, 1.5)
  )
})

test_that("bandx() refuses bad input in words", {
  errs(
    "bandx",
    list(
      list("`x` must be band labels", 1:3),
      list("`clsup` must be TRUE or FALSE", "1-2", clsup = NA),
      list("`rule` must be one of", "1-2", rule = "mean"),
      list("`fac` must be one number", "1-2", fac = 0),
      list("wrong side of its other end", "Over 60", cap = 10),
      list("`cnts` must be numbers from 0", "1-2", cnts = -1),
      list("cannot be read: lots", c("1-2", "lots"), c(1, 1)),
      list("`stat` must be a function", "1-2", 1, stat = "mean"),
      list("needs `cnts` and an open top band", "1-2", 1, stat = mean),
      list("`caps` must be numbers", "Over 2", 1, stat = mean, caps = "a"),
      list("must give one number", "Over 2", 1, stat = range)
    )
  )
})

test_that("pbbnd() bounds rounded and suppressed cells and their change", {
  r <- pbbnd(c("1,240", " 35", "<5", "x", NA, ""), rnd = 5, suplim = 10)
  b <- r$bnd
  expect_identical(b$value, c(1240, 35, NA, NA, NA, NA))
  expect_identical(b$lower, c(1237.5, 32.5, 0, 0, NA, NA))
  expect_identical(b$upper, c(1242.5, 37.5, 4, 9, NA, NA))
  expect_identical(
    b$status,
    c("rounded", "rounded", "suppressed", "suppressed", "missing", "missing")
  )
  expect_null(r$chg)
  q <- pbbnd(
    c(120, 15),
    rnd = 5,
    kind = "random",
    prev = c(100, 20),
    pop = 5e4,
    prepop = c(4e4, 5e4),
    per = 1e4
  )
  expect_identical(q$bnd$lower, c(116, 11))
  c1 <- q$chg
  expect_identical(c1$chglo, c(116 - 104, 11 - 24))
  expect_identical(c1$chgup, c(124 - 96, 19 - 16))
  expect_equal(c1$pctlo, 100 * (c(116, 11) / c(104, 24) - 1), tolerance = 1e-15)
  expect_equal(c1$rateup, 1e4 * c(124, 19) / 5e4, tolerance = 1e-15)
  expect_equal(c1$preratelo, 1e4 * c(96, 16) / c(4e4, 5e4), tolerance = 1e-15)
  e <- pbbnd(c(3, 4), prev = factor(c("2", "<3")))$chg
  expect_identical(e$chglo, c(1, 2))
  expect_identical(e$pctup, c(50, NA))
  expect_named(
    pbbnd(1, prev = 1, pop = 10)$chg,
    c("chglo", "chgup", "pctlo", "pctup", "ratelo", "rateup")
  )
  expect_output(print(q), "what the published cells allow")
})

test_that("pbbnd() refuses bad input in words", {
  errs(
    "pbbnd",
    list(
      list("`supmar` must be text", 1, supmar = 1),
      list("`x` must be numbers or text cells", list(1)),
      list("`x` must be finite", Inf),
      list(
        "cannot read these cells as numbers or suppression marks: a, b",
        c("a", "1", "b")
      ),
      list("`x` must hold counts from 0", -1),
      list("`rnd` must be above 0", 1, rnd = 0),
      list("`kind` must be one of", 1, kind = "up"),
      list("a suppressed cell has no limit", "x"),
      list("`prev` must have one cell per cell", 1:2, prev = 1),
      list("`prev` must hold counts from 0", 1, prev = -1),
      list("a suppressed cell of `prev`", 1, prev = "x"),
      list(
        "`pop` must be one number or one per cell, 2",
        1:2,
        prev = 1:2,
        pop = 1:3
      ),
      list(
        "`prepop` must be one number or one per cell",
        1:2,
        prev = 1:2,
        prepop = 1:3
      ),
      list("`per` must be one number", 1, prev = 1, per = 0)
    )
  )
})

test_that("yoypx() gives changes with exact tests and adjusted p-values", {
  d <- data.frame(
    year = c(2019:2023, 2019:2023, 2021),
    svc = c(rep("A", 5), rep("B", 5), "A"),
    n = c(120, 131, 98, 105, 140, 30, 28, 41, 25, 0, 2)
  )
  r <- yoypx(d$n, d$year, d$svc, adj = "holm")
  expect_s3_class(r, "yoypx")
  a <- r[r$group == "A", ]
  expect_identical(a$value, c(120, 131, 100, 105, 140))
  expect_identical(a$prev, c(NA, 120, 131, 100, 105))
  expect_equal(a$pctchg[2], 100 * 11 / 120, tolerance = 1e-15)
  ok <- !is.na(r$prev) & r$value + r$prev > 0
  p <- mapply(
    function(v, b) binom.test(v, v + b)$p.value,
    r$value[ok],
    r$prev[ok]
  )
  expect_equal(r$p[ok], p, tolerance = 1e-13)
  expect_equal(r$padj[ok], p.adjust(p, "holm"), tolerance = 1e-13)
  expect_identical(r$signif[ok], p.adjust(p, "holm") < 0.05)
  # the exact interval of the ratio, as a percent change
  v <- 131
  b <- 120
  lo <- qbeta(0.025, v, b + 1)
  hi <- qbeta(0.975, v + 1, b)
  expect_equal(a$pctlo[2], 100 * (lo / (1 - lo) - 1), tolerance = 1e-13)
  expect_equal(a$pctup[2], 100 * (hi / (1 - hi) - 1), tolerance = 1e-13)
  z <- r[r$group == "B", ]
  expect_identical(z$flag, c("no earlier period", NA, NA, NA, NA))
  expect_identical(
    z$pctup[5],
    100 *
      (qbeta(0.975, 1, 25) /
        (1 - qbeta(0.975, 1, 25)) -
        1)
  )
  expect_identical(z$pctlo[5], -100)
  f <- yoypx(c(0, 5, 0, 0), 1:4, minbase = 3)
  expect_identical(
    f$flag,
    c("no earlier period", "earlier period is 0", NA, "earlier period is 0")
  )
  expect_identical(f$p[4], NA_real_)
  expect_identical(yoypx(c(2, 5), 1:2, minbase = 3)$flag[2], "small base")
  expect_identical(yoypx(c(4, 5), 1:2)$p[1], NA_real_)
  expect_output(print(r), "change on the earlier period")
})

test_that("every adjustment and the rounding envelope agree with base R", {
  set.seed(3)
  x <- rpois(12, 40)
  for (m in c("none", "bonferroni", "hochberg", "BH", "BY")) {
    r <- yoypx(x, 1:12, adj = m)
    ok <- !is.na(r$p)
    expect_equal(r$padj[ok], p.adjust(r$p[ok], m), tolerance = 1e-14, info = m)
  }
  # an even split: the p-value is 1
  expect_identical(yoypx(c(7, 7), 1:2)$p[2], 1)
  e <- yoypx(c(100, 120), c(2020, 2021), rnd = 10)
  expect_identical(c(e$envlo[2], e$envup[2]), c(10, 30))
  expect_equal(e$envplo[2], 100 * (115 / 105 - 1), tolerance = 1e-14)
  expect_identical(e$cpctlo[2], min(e$pctlo[2], e$envplo[2]))
  expect_identical(e$cpctup[2], max(e$pctup[2], e$envpup[2]))
  expect_identical(e$cpctlo[1], NA_real_)
  er <- yoypx(c(100, 120), c(2020, 2021), rnd = 10, kind = "random")
  expect_identical(er$envlo[2], 120 - 9 - 109)
  # other units: means of a cell's rows, points for percentages
  k <- yoypx(c(1.5, 2.5, 4), c(1, 1, 2), units = "continuous")
  expect_identical(k$value, c(2, 4))
  expect_false("p" %in% names(k))
  pp <- yoypx(c(10, 12.5), c(1, 3), lag = 2, units = "percent")
  expect_identical(pp$ppchg, c(NA, 2.5))
})

test_that("yoypx() refuses bad input in words", {
  errs(
    "yoypx",
    list(
      list("`units` must be one of", 1, 1, units = "n"),
      list("`x` must be at least 0", -1, 1),
      list("`prd` must have length 2", 1:2, 1),
      list("`lag` must be one number", 1, 1, lag = 0),
      list("`minbase` must be one number from 0", 1, 1, minbase = -1),
      list("`rnd` must be above 0", 1, 1, rnd = 0),
      list("`rnd` goes with counts", 1, 1, units = "continuous", rnd = 5),
      list("`adj` must be one of", 1, 1, adj = "fdr"),
      list("`alpha` must be one number", 1, 1, alpha = 1),
      list("`by` must have one label per value", 1:2, 1:2, by = 1)
    )
  )
})
