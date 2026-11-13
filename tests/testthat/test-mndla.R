# mndla(): against the classes and counts worked out here in base R, the
# exact interval of binom.test(), and the refusals.

test_that("mndla() gives the three denominators by group", {
  set.seed(6)
  d <- data.frame(
    id = sample(1:60, 200, TRUE),
    yr = sample(2020:2022, 200, TRUE),
    days = rgeom(200, 0.08) + 1
  )
  r <- mndla(d$days, d$id, d$yr)
  tb <- r$tab
  expect_identical(unique(tb$group), c("2020", "2021", "2022", "all"))
  expect_identical(tb$den[1:3], c("stays", "people", "cumul"))
  q <- d$days > 15
  expect_identical(r$cls, as.integer(q))
  for (y in c(2020:2022, NA)) {
    s <- if (is.na(y)) rep(TRUE, 200) else d$yr == y
    g <- if (is.na(y)) "all" else as.character(y)
    row <- tb[tb$group == g, ]
    expect_equal(
      row$n,
      c(sum(s), rep(length(unique(d$id[s])), 2)),
      tolerance = 0
    )
    cum <- tapply(d$days[s], d$id[s], sum)
    expect_equal(
      row$k,
      c(sum(q & s), length(unique(d$id[s & q])), sum(cum > 15)),
      tolerance = 0
    )
  }
  b <- binom.test(tb$k[1], tb$n[1])$conf.int
  expect_equal(c(tb$lo[1], tb$hi[1]), as.numeric(b), tolerance = 1e-13)
  expect_identical(r$rule, "duration")
  expect_output(print(r), "duration")
})

test_that("mndla() applies Rules 43 and 44, alerts, eligibility and contact", {
  r <- mndla(
    c(10, 20, 20, 30),
    hours = c(1, 2, 3, 0.5),
    miss = c(100, 100, 100, 99)
  )
  expect_identical(r$cls, c(1L, 2L, 0L, 0L))
  expect_identical(r$tab$sol, 1)
  expect_identical(r$tab$tor, 1)
  expect_identical(r$tab$group, "all")
  expect_output(print(r), "Rules 43 and 44")
  a <- mndla(
    c(5, 16, 20, 3),
    alert = data.frame(
      mh = c("Yes", "no", "No", " yes "),
      sr = c(TRUE, FALSE, FALSE, FALSE)
    )
  )
  expect_identical(a$alerts, c(2L, 0L, 0L, 1L))
  expect_identical(a$cls, c(1L, 1L, 1L, 0L))
  m <- mndla(c(5, 16), alert = matrix(c(1, 0, 1, 0), 2), amin = 2)
  expect_identical(m$cls, c(1L, 1L))
  # matrices of every kind, and whole numbers
  for (al in list(
    matrix(c(TRUE, FALSE, TRUE, FALSE), 2),
    matrix(c("yes", "no", "Yes", "no"), 2),
    matrix(c(1L, 0L, 1L, 0L), 2)
  )) {
    expect_identical(mndla(c(5, 6), alert = al)$alerts, c(2L, 0L))
  }
  expect_identical(mndla(c(16, 16), elig = c(1L, 0L))$cls, c(1L, 0L))
  e <- mndla(c(16, 20, 30), elig = c(1, 0, 1), excl = c(FALSE, FALSE, TRUE))
  expect_identical(e$cls, c(1L, 0L, 0L))
  l <- mndla(c(16, 20), alert = list(c("yes", "no")), amin = 1)
  expect_identical(l$cls, c(1L, 1L))
  # aggregate rows with counts: no interval for fractional counts
  w <- mndla(c(10, 20), w = c(2.5, 1))
  expect_identical(w$tab$n, 3.5)
  expect_true(is.na(w$tab$lo))
  z <- mndla(c(1, 2), level = 0.9)$tab
  expect_identical(c(z$k, z$lo), c(0, 0))
  f <- mndla(c(20, 30))$tab
  expect_identical(f$hi, 1)
})

test_that("mndla() refuses bad input in words", {
  errs(
    "mndla",
    list(
      list("`days` must be at least 0", -1),
      list("`thr` must be one number from 0", 1, thr = -1),
      list("`hours` and `miss` go together", 1, hours = 1),
      list("`miss` must be percentages", 1, hours = 1, miss = 101),
      list(
        "`alert` goes with the duration rule",
        1,
        hours = 1,
        miss = 100,
        alert = TRUE
      ),
      list(
        "`alert` must have one row per stay",
        1:2,
        alert = matrix(TRUE, 3, 2)
      ),
      list("`alert` must have one value per stay", 1:2, alert = TRUE),
      list("`alert` must be TRUE or FALSE", 1:2, alert = c("yes", "maybe")),
      list("`alert` must be TRUE or FALSE", 1, alert = "yesss"),
      list("`alert` must be TRUE or FALSE", 1, alert = factor("yes")),
      list("`elig` must be TRUE or FALSE", 1, elig = NA),
      list("`excl` must be TRUE or FALSE", 1, excl = NA_real_),
      list("`amin` must be one whole number", 1, amin = 0),
      list("`by` must have one label per value", 1:2, by = 1),
      list("`w` must have length 2", 1:2, w = 1),
      list("`level` must be one number", 1, level = 2)
    )
  )
})
