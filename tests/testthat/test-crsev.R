# crsev() and zcomp(): against the sums and standardisations done here in
# base R, the worked example of rmorie's CrimeSeverityIndex(), and the
# refusals.

test_that("crsev() sums weighted counts per unit and period", {
  set.seed(11)
  d <- expand.grid(
    off = c("a", "b", "c"),
    unit = c("u1", "u2", "u3"),
    time = 2019:2021,
    stringsAsFactors = FALSE
  )
  d$x <- rpois(nrow(d), 20)
  d$pop <- c(u1 = 5000, u2 = 7000, u3 = 9000)[d$unit]
  w <- c(b = 3, a = 0.5, c = 40)
  r <- crsev(d$x, d$off, w, d$unit, d$time, d$pop, per = 1e3, base = 4)
  key <- paste(d$unit, d$time)
  sev <- c(tapply(d$x * w[d$off], key, sum))
  cnt <- c(tapply(d$x, key, sum))
  pop <- c(tapply(d$pop, key, unique))
  k <- paste(r$tab$unit, r$tab$time)
  expect_identical(k, sort(unique(key)))
  expect_equal(r$tab$sev, unname(sev[k]), tolerance = 1e-15)
  expect_equal(r$tab$count, unname(cnt[k]), tolerance = 1e-15)
  expect_equal(r$tab$srate, unname(sev[k] / pop[k] * 1e3), tolerance = 1e-14)
  expect_equal(r$tab$rate, unname(cnt[k] / pop[k] * 1e3), tolerance = 1e-14)
  expect_equal(r$tab$idx, 100 * r$tab$srate / r$tab$srate[4], tolerance = 1e-14)
  expect_identical(r$w, data.frame(off = c("a", "b", "c"), w = c(0.5, 3, 40)))
  expect_identical(r$base, 4L)
  expect_output(print(r), "index 100 at row 4")
})

test_that("crsev() reproduces the worked example and its options", {
  # rmorie::CrimeSeverityIndex(rbind(c(10, 2), c(8, 3)), c(1, 50),
  # c(1000, 1000))$index is 100 and 143.6364 (158 / 110 * 100)
  r <- crsev(
    c(10, 2, 8, 3),
    c("t", "r", "t", "r"),
    c(t = 1, r = 50),
    time = c(1, 1, 2, 2),
    pop = rep(1000, 4),
    base = 1
  )
  expect_equal(r$tab$idx, c(100, 15800 / 110), tolerance = 1e-15)
  expect_named(r$tab, c("time", "count", "sev", "pop", "rate", "srate", "idx"))
  # weights as rate times sentence, and unnamed weights in label order
  s <- crsev(
    c(10, 2, 8, 3),
    c("t", "r", "t", "r"),
    c(r = 5, t = 0.5),
    time = c(1, 1, 2, 2),
    sent = c(10, 2)
  )
  expect_identical(s$w$w, c(50, 1))
  expect_identical(s$tab$sev, c(110, 158))
  expect_named(s$tab, c("time", "count", "sev"))
  expect_true(is.na(s$base))
  expect_output(print(s), "^crsev: severity-weighted counts\n")
  # no unit or time: one cell; rows of one offence add up
  o <- crsev(c(1, 2, 3), c("a", "a", "b"), c(2, 1), base = 1)
  expect_identical(o$tab, data.frame(count = 6, sev = 9, idx = 100))
  # a zero base gives no index
  z <- crsev(c(0, 4), c("a", "a"), 1, time = 1:2, base = 1)
  expect_identical(z$tab$idx, c(NA_real_, NA_real_))
})

test_that("crsev() refuses bad input in words", {
  errs(
    "crsev",
    list(
      list("`x` must be at least 0", -1, "a", 1),
      list("one label per value", 1:2, "a", 1),
      list("`w` has no value for b", 1:2, c("a", "b"), c(a = 1)),
      list(
        "`w` names a twice",
        1:2,
        c("a", "b"),
        stats::setNames(c(1, 2, 1), c("a", "a", "b"))
      ),
      list("one value per label, 2", 1:2, c("a", "b"), 1),
      list("`w` must be at least 0", 1, "a", -1),
      list("`sent` has no value for a", 1, "a", 1, sent = c(b = 1)),
      list(
        "`unit` must be a vector",
        1:4,
        rep("a", 4),
        1,
        unit = matrix(1:4, 2)
      ),
      list("same on every row", 1:2, c("a", "b"), 1:2, pop = c(10, 20)),
      list("`pop` must be greater than 0", 1, "a", 1, pop = 0),
      list("`per` must be one number", 1, "a", 1, per = 0),
      list("`base` must be one whole number", 1, "a", 1, base = 0.5),
      list(
        "`base` must be a row of the result, from 1 to 1",
        1,
        "a",
        1,
        base = 2
      )
    )
  )
})

test_that("zcomp() standardises each category and sums with weights", {
  set.seed(5)
  d <- expand.grid(
    unit = letters[1:6],
    cat = c("x", "y", "z"),
    stringsAsFactors = FALSE
  )
  d$n <- rpois(nrow(d), 12)
  d <- d[-3, ]
  r <- zcomp(d$n, d$unit, d$cat, w = c(z = 1, x = 2, y = 0.5))
  m <- tapply(d$n, list(d$unit, d$cat), sum)
  m[is.na(m)] <- 0
  z <- scale(m)
  expect_equal(r$z, z, ignore_attr = TRUE, tolerance = 1e-14)
  expect_identical(dimnames(r$z), dimnames(m))
  comp <- drop(z %*% c(2, 0.5, 1))
  expect_equal(r$tab$comp, unname(comp), tolerance = 1e-14)
  expect_identical(r$tab$rank, as.integer(rank(-comp, ties.method = "min")))
  expect_equal(r$cat$mean, unname(colMeans(m)), tolerance = 1e-15)
  expect_output(print(r), "composite of 3 standardised categories")
  # a category that does not vary adds 0; ties share the lower rank
  e <- zcomp(
    c(1, 2, 2, 7, 7, 7),
    rep(c("a", "b", "c"), 2),
    rep(c("p", "q"), each = 3)
  )
  expect_identical(e$cat$sd[2], 1)
  expect_identical(e$tab$rank, c(3L, 1L, 1L))
  errs(
    "zcomp",
    list(
      list("at least two units", 1, "a", "p"),
      list("more than 1e8 cells", rep(1, 10001), 1:10001, 1:10001),
      list(
        "`w` has no value for q",
        1:2,
        c("a", "b"),
        c("p", "q"),
        w = c(p = 1)
      ),
      list("`x` must be at least 0", -1:0, c("a", "b"), c("p", "p"))
    )
  )
})
