# agemix(), cumrsk() and jpipe(): against cov.wt(), the life-table
# products written out, windows summed one by one, and the refusals.

test_that("agemix() gives the aggregate curve and the covariance of types", {
  ages <- 12:30
  m <- cbind(
    al = dnorm(ages, 17, 2),
    lcp = rep(0.08, 19),
    late = dnorm(ages, 26, 3)
  )
  w <- c(0.6, 0.15, 0.25)
  r <- agemix(m, w, ages)
  expect_equal(r$curve$agg, drop(m %*% w), tolerance = 1e-15)
  expect_named(r$curve, c("age", "agg", "al", "lcp", "late"))
  cw <- cov.wt(t(m), wt = w, method = "ML")$cov
  expect_equal(r$cov, cw, ignore_attr = TRUE, tolerance = 1e-14)
  expect_identical(
    dimnames(r$cov),
    list(as.character(ages), as.character(ages))
  )
  # one type: the covariance is 0 whatever the curve
  one <- agemix(m, c(1, 0, 0))
  expect_true(all(one$cov == 0))
  expect_identical(one$curve$age, as.numeric(1:19))
  expect_null(r$panel)
  # unnamed columns are named by position
  u <- agemix(data.frame(1:3, 3:1)[, 1:2, drop = FALSE] * 1, c(0.5, 0.5))
  expect_named(u$curve, c("age", "agg", "X1.3", "X3.1"))
  v <- agemix(matrix(1:4, 2) + 0, c(0.5, 0.5))
  expect_identical(names(v$curve)[3:4], c("type1", "type2"))
  expect_output(print(r), "mixture of types")
})

test_that("agemix() tests one type on a panel", {
  set.seed(4)
  a <- rpois(60, 2)
  b <- a + rpois(60, 1)
  r <- agemix(panel = cbind(a, b), nperm = 199, seed = 2)$panel
  expect_equal(r[["cov"]], mean(a * b) - mean(a) * mean(b), tolerance = 1e-14)
  expect_identical(r[["p"]], 1 / 200)
  expect_identical(r[["n"]], 60)
  q <- agemix(panel = data.frame(a, sample(b)), nperm = 0)$panel
  expect_true(is.na(q[["p"]]))
  set.seed(3)
  x <- agemix(panel = cbind(a, rev(b)), nperm = 99)$panel
  set.seed(3)
  expect_identical(agemix(panel = cbind(a, rev(b)), nperm = 99)$panel, x)
  expect_output(print(agemix(panel = cbind(a, b), nperm = 0)), "cov")
})

test_that("agemix() refuses bad input in words", {
  errs(
    "agemix",
    list(
      list("give `curves` and `shares`, or `panel`"),
      list("`curves` must be a numeric matrix", "a", 1),
      list("`curves` must have at least one column", list(), 1),
      list("`curves` must hold finite numbers", matrix(c(1, NA), 1), c(1, 0)),
      list(
        "`curves` may have at most 30 types",
        matrix(0, 1, 31),
        rep(1 / 31, 31)
      ),
      list("`shares` must have length 2", matrix(1, 2, 2), 1),
      list("`shares` must sum to 1", matrix(1, 2, 2), c(0.5, 0.4)),
      list("`ages` must have length 2", matrix(1, 2, 2), c(0.5, 0.5), ages = 1),
      list("`panel` must have two columns", panel = matrix(1, 3, 3)),
      list("at least three individuals", panel = matrix(1, 2, 2)),
      list(
        "`nperm` must be one whole number",
        panel = matrix(1, 3, 2),
        nperm = -1
      )
    )
  )
})

test_that("cumrsk() is one minus the product of the complements", {
  age <- 18:25
  ev <- c(120, 110, 90, 80, 60, 50, 40, 30)
  pop <- rep(10000, 8)
  r <- cumrsk(age, ev, pop)
  nr <- numeric(8)
  s <- 1
  cr <- v <- numeric(8)
  vs <- 0
  for (i in 1:8) {
    nr[i] <- pop[i] * (1 - if (i > 1) cr[i - 1] else 0)
    h <- ev[i] / nr[i]
    s <- s * (1 - h)
    cr[i] <- 1 - s
    vs <- vs + h / (nr[i] * (1 - h))
    v[i] <- vs
  }
  expect_equal(r$cumrisk, cr, tolerance = 1e-14)
  expect_equal(r$atrisk, nr, tolerance = 1e-14)
  z <- qnorm(0.975)
  expect_equal(r$lo, 1 - (1 - cr) * exp(z * sqrt(v)), tolerance = 1e-14)
  expect_equal(r$hi, 1 - (1 - cr) * exp(-z * sqrt(v)), tolerance = 1e-14)
  expect_true(all(r$lo < r$cumrisk & r$cumrisk < r$hi))
  d <- cumrsk(age, ev, pop, dth = rep(0.001, 8))
  expect_true(all(is.na(d$lo)))
  expect_lt(d$cumrisk[8], r$cumrisk[8])
  a <- cumrsk(18:19, c(5, 5), c(100, 100), atrisk = c(50, 40))
  expect_identical(a$hazard, c(0.1, 0.125))
  expect_s3_class(r, "cumrsk")
  expect_output(print(r), "first event by age")
  errs(
    "cumrsk",
    list(
      list("`age` must be increasing", c(19, 18), 1:2, c(10, 10)),
      list("`ev` must have length 2", 18:19, 1, c(10, 10)),
      list("`pop` must be greater than 0", 18, 1, 0),
      list("`dth` must be death rates", 18, 1, 10, dth = 1),
      list("more first events", 18, 20, 10),
      list("`level` must be one number", 18, 1, 10, level = 1)
    )
  )
})

test_that("jpipe() divides lagged window sums", {
  num <- c(5, 7, 6, 9, 8, 10, 4, 6)
  den <- c(20, 22, 19, 25, 24, 26, 18, 21)
  g <- rep(c("b", "a"), each = 4)
  r <- jpipe(num, den, win = 2, lag = 1, grp = g)
  win <- function(x, e, w) if (e < w) NA else sum(x[(e - w + 1):e])
  a <- 5:8
  want <- vapply(1:4, function(e) win(num[a], e, 2), 0)
  wd <- vapply(1:4, function(e) if (e - 1 < 2) NA else win(den[a], e - 1, 2), 0)
  expect_identical(r$group, rep(c("a", "b"), each = 4))
  expect_identical(r$num[1:4], want)
  expect_identical(r$den[1:4], wd)
  expect_equal(r$rate[1:4], want / wd, tolerance = 1e-15)
  expect_identical(r$index, rep(1:4, 2))
  o <- jpipe(num, den, win = 3)
  expect_identical(o$group, rep("all", 8))
  expect_equal(o$rate[3], 18 / 61, tolerance = 1e-15)
  expect_output(print(o), "lagged rates")
  errs(
    "jpipe",
    list(
      list("`num` must be at least 0", -1, 1),
      list("`den` must have length 2", 1:2, 1),
      list("`win` must be one whole number from 1", 1, 1, win = 0),
      list("`lag` must be one whole number from 0", 1, 1, lag = -1),
      list("`grp` must have one label per value", 1:2, 1:2, grp = 1)
    )
  )
})
