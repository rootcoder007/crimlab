# Values of the corrections kernels are checked in C++ (src/t_corr.cpp);
# here what R forms itself, the labels, the messages and rmorie parity
# (rmorie's values stored in ref/corr.txt by dev/ref/corr.R, version named
# there).

test_that("incap labels its groups and adds the aggregate row", {
  lam <- c(a = 2, b = 10, c = 0.5)
  r <- incap(lam, c(0.1, 0.3, 1), c(1, 2.5, 0))
  expect_identical(r[["group"]], c("a", "b", "c"))
  r <- incap(c(low = 2, high = 10), 0.1, 1, share = c(0.8, 0.2))
  expect_identical(r[["group"]], c("low", "high", "all"))
  expect_true(is.na(r[["marg"]][3]))
})

test_that("sntbnd names its terms and divides the observed rates in R", {
  y <- c(1, 0, 1, 1, 0, 0)
  z <- c("a", "a", "b", "b", "b", "a")
  w <- c(1, 2, 1, 0.5, 1, 1.5)
  r <- sntbnd(y, z, treat = "b", w = w)
  expect_identical(r[["term"]], c("a", "b", "effect"))
  expect_equal(r[["observed"]], c(1 / 4.5, 1.5 / 2.5, 1.5 / 2.5 - 1 / 4.5),
               tolerance = 1e-15)
  # the treatment picked by name flips the effect
  s <- sntbnd(y, z, treat = "a", w = w)
  expect_identical(s[["term"]], c("b", "a", "effect"))
})

test_that("corrections functions stop on bad input, naming the argument", {
  errs("incap", list(
    list("`lam`", 0, 0.1, 1),
    list("`q`", 1, 0, 1),
    list("at most 1", 1, 1.5, 1),
    list("`s`", 1, 0.5, -1),
    list("`q` must have length 1 or 3", c(1, 2, 3), c(0.1, 0.2), 1),
    list("`s` must have length 1 or 3", c(1, 2, 3), 0.1, c(1, 2)),
    list("sum to 1", c(1, 2), 0.1, 1, share = c(0.5, 0.6)),
    list("`share`", c(1, 2), 0.1, 1, share = c(1.5, -0.5))
  ))
  errs("incapr", list(
    list("never increase", c(1, 2), 0, 1),
    list("`s`", c(2, 1), 0, 0),
    list("`t0`", c(2, 1), -1, 1),
    list("`t0`", c(2, 1), 0.5, 1),
    list("length at least", c(2, 1), 1, 1),
    list("`r`", c(2, 1), 0, 1, r = 2),
    list("`lam`", c(2, NA), 0, 1)
  ))
  tt <- c(TRUE, TRUE)
  ff <- c(FALSE, FALSE)
  errs("sellab", list(
    list("`rel`", c(1, 0), c(TRUE, NA), c(TRUE, TRUE)),
    list("`rule`", c(1, 0), c(TRUE, TRUE), c(1, 0)),
    list("0 or 1", c(1, 2), c(TRUE, TRUE), c(TRUE, TRUE)),
    list("0 or 1", c(1, NA), c(TRUE, TRUE), c(TRUE, TRUE)),
    list("`y`", "a", TRUE, TRUE),
    list("`w`", c(1, 0), tt, tt, w = c(1, -1)),
    list("released nobody", c(1, 0), ff, tt),
    list("releases nobody", c(1, 0), tt, ff)
  ))
  errs("losbnd", list(list("`done`", numeric(0), 1), list("`open`", 1, -1)))
  errs("sntbnd", list(
    list("0 or 1", c(0, 2), c("a", "b")),
    list("exactly two", c(0, 1), c("a", "a")),
    list("`z`", c(0, 1), c("a", NA)),
    list("`treat`", c(0, 1), c("a", "b"), treat = "c"),
    list("positive total", c(0, 1), c("a", "b"), "b", w = c(0, 0))
  ))
  errs("hzdsel", list(
    list("`s`", 1, 0.5, 0.1, c(1, 1), c(1, 1)),
    list("`h`", 0.5, 1.5, 0.1, c(1, 1), c(1, 1)),
    list("below `h`", 0.5, 0.5, 0.6, c(1, 1), c(1, 1)),
    list("`surh`", 0.5, 0.5, 0.1, c(0, 1), c(1, 1)),
    list("`surl`", 0.5, 0.5, 0.1, c(1, 1), c(1, 1.2)),
    list("`surh`", 0.5, 0.5, 0.1, 1, c(1, 1)),
    list("`surl`", 0.5, 0.4, 0.2, c(0.9, 0.8), c(0.9, 1.2)),
    list("`surh`", 0.5, 0.4, 0.2, c(1.2, 0.8), c(0.9, 0.8))
  ))
})

test_that("corrections results agree with rmorie", {
  ref <- refget("corr")
  # from rmorie::morie_incapacitation(c(2, 10), 0.1, 1,
  #   shares = c(0.8, 0.2))
  o <- ref$inc
  d <- incap(c(2, 10), 0.1, 1, share = c(0.8, 0.2))
  expect_equal(d[["rate"]][1:2], o$rate, tolerance = 1e-14)
  expect_equal(d[["prev"]][1:2], o$prev, tolerance = 1e-14)
  expect_equal(d[["marg"]][1:2], o$marg, tolerance = 1e-14)
  expect_equal(d[["prev"]][3], o$agg, tolerance = 1e-14)

  # from rmorie::morie_incapacitation_career(lam, 2, 3, 0.25): prevented,
  # prevented_net, upper, lower, later
  lam <- c(12, 10, 8, 6, 5, 4, 3, 2, 2, 1)
  d <- incapr(lam, 2, 3, 0.25)
  expect_equal(unname(d[c("prev", "net", "upper", "lower", "later")]),
               ref$car, tolerance = 1e-14)
  # rmorie returns NA for `later` when the path is exactly t0 + s + 1 long
  # (morie_incapacitation_career(c(5, 4, 3), 0, 2)), although the sum
  # needs ages up to t0 + s only; crimlab computes it
  expect_true(is.na(ref$car3))
  expect_equal(incapr(c(5, 4, 3), 0, 2)[["later"]], 7, tolerance = 1e-15)

  set.seed(4)
  risk <- runif(500)
  y <- rbinom(500, 1, risk)
  # from rmorie::morie_selective_labels(y, risk < 0.7, risk < 0.8)
  o <- ref$sel
  expect_equal(c(sum(risk), sum(y)), o$inp, tolerance = 1e-12,
               info = "inputs changed")
  d <- sellab(y, risk < 0.7, risk < 0.8)
  expect_equal(d[["obs"]], o$obs, tolerance = 1e-14)
  expect_equal(d[["unobs"]], o$unobs, tolerance = 1e-14)
  expect_equal(unname(d[c("lower", "upper")]), o$bounds, tolerance = 1e-14)
  expect_equal(d[["naive"]], o$naive, tolerance = 1e-14)

  # rmorie::morie_backlog_censoring(c(30, 45, 60, 90, 120),
  #   c(100, 150, 200))
  d <- losbnd(c(30, 45, 60, 90, 120), c(100, 150, 200))
  expect_equal(d[["lower"]], ref$los$lower, tolerance = 1e-14)
  expect_equal(d[["bias"]], ref$los$bias, tolerance = 1e-14)

  set.seed(1)
  z <- sample(c("community", "custody"), 500, TRUE, c(0.6, 0.4))
  y <- rbinom(500, 1, ifelse(z == "custody", 0.55, 0.35))
  # from rmorie::morie_sentence_effect_bounds(y, z)
  o <- ref$snt
  expect_equal(c(sum(z == "custody"), sum(y)), o$inp, tolerance = 1e-12,
               info = "inputs changed")
  d <- sntbnd(y, z, treat = "custody")
  expect_equal(d[["lower"]][3], o$lower, tolerance = 1e-14)
  expect_equal(d[["upper"]][3], o$upper, tolerance = 1e-14)
  expect_equal(d[["observed"]][3], o$naive, tolerance = 1e-14)

  # from rmorie::morie_hazard_selection(0.5, 0.5, 0.1, c(0.5, 0.8),
  #   c(0.9, 0.9))
  d <- hzdsel(0.5, 0.5, 0.1, c(0.5, 0.8), c(0.9, 0.9))
  expect_equal(unname(d[c("haz0", "haz1")]), ref$hz$haz, tolerance = 1e-14)
  expect_equal(d[["ratio"]], ref$hz$ratio, tolerance = 1e-14)
})
