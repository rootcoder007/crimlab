test_that("incap matches the renewal cycle by hand", {
  lam <- c(a = 2, b = 10, c = 0.5)
  q <- c(0.1, 0.3, 1)
  s <- c(1, 2.5, 0)
  r <- incap(lam, q, s)
  # one cycle: 1 / q crimes in a free spell of 1 / (lam q) years, then s
  # years inside
  cyc <- (1 / q) / (1 / (lam * q) + s)
  expect_equal(r[["rate"]], unname(cyc), tolerance = 1e-14)
  # the share prevented is the time inside over the cycle length
  expect_equal(r[["prev"]], unname(s / (1 / (lam * q) + s)), tolerance = 1e-14)
  # one more year, by difference
  r1 <- incap(lam, q, s + 1)
  expect_equal(r[["marg"]], r[["rate"]] - r1[["rate"]], tolerance = 1e-12)
  expect_equal(r[["group"]], c("a", "b", "c"))
  expect_equal(r[["prev"]][3], 0)
  expect_true(all(r[["prev"]] < 1))
})

test_that("incap: the documented monotone claims hold on a grid", {
  for (lam in c(0.2, 1, 5, 40)) for (q in c(0.01, 0.2, 1)) {
    s <- seq(0, 30, by = 0.5)
    r <- incap(rep(lam, length(s)), q, s)
    expect_true(all(diff(r[["rate"]]) < 0))
    # the marginal year buys less each time (tested, not a theorem)
    expect_true(all(diff(r[["marg"]]) < 0))
    expect_true(all(r[["marg"]] > 0))
  }
  # surer sentences lower the rate
  r <- incap(rep(4, 5), c(0.1, 0.2, 0.4, 0.8, 1), 2)
  expect_true(all(diff(r[["rate"]]) < 0))
  # a uniform sentence prevents a larger share of the high-rate group
  expect_true(all(diff(incap(c(0.5, 2, 10, 50), 0.2, 1)[["prev"]]) > 0))
})

test_that("incap aggregates over groups", {
  r <- incap(c(low = 2, high = 10), 0.1, 1, share = c(0.8, 0.2))
  expect_equal(r[["group"]], c("low", "high", "all"))
  free <- 0.8 * 2 + 0.2 * 10
  inc <- 0.8 * 2 / 1.2 + 0.2 * 10 / 2
  expect_equal(r[["lam"]][3], free, tolerance = 1e-15)
  expect_equal(r[["rate"]][3], inc, tolerance = 1e-15)
  expect_equal(r[["prev"]][3], 1 - inc / free, tolerance = 1e-15)
  expect_true(is.na(r[["marg"]][3]))
})

test_that("incapr sums the path and its bounds", {
  lam <- c(12, 10, 8, 6, 5, 4, 3, 2, 2, 1)
  r <- incapr(lam, t0 = 2, s = 3, r = 0.25)
  expect_equal(r[["prev"]], 8 + 6 + 5)
  expect_equal(r[["net"]], 0.75 * 19)
  expect_equal(r[["upper"]], 3 * 8)
  expect_equal(r[["lower"]], 3 * 4)
  expect_equal(r[["later"]], 6 + 5 + 4)
  # the shortest path that holds the sentence still has a later value
  r <- incapr(c(5, 4, 3), t0 = 0, s = 2)
  expect_equal(unname(r), c(9, 9, 10, 6, 7))
  # constant rate: the bounds meet
  r <- incapr(rep(3, 6), t0 = 1, s = 4, r = 1)
  expect_equal(unname(r), c(12, 0, 12, 12, 12))
  set.seed(2)
  for (k in 1:200) {
    lam <- sort(rexp(sample(3:20, 1)), decreasing = TRUE)
    s <- sample(seq_len(length(lam) - 1), 1)
    t0 <- sample(0:(length(lam) - 1 - s), 1)
    rr <- runif(1)
    r <- incapr(lam, t0, s, rr)
    expect_lte(r[["lower"]], r[["prev"]] * (1 + 1e-12))
    expect_lte(r[["prev"]], r[["upper"]] * (1 + 1e-12))
    expect_lte(r[["later"]], r[["prev"]] * (1 + 1e-12))
    expect_lte(r[["net"]], r[["upper"]] * (1 + 1e-12))
  }
})

test_that("sellab bounds are the extremes over every filling of the unknown", {
  y <- c(1, 0, 1, 0, NA, NA, NA, 0)
  rel <- c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE)
  rule <- c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE)
  w <- c(1, 2, 1, 1, 0.5, 1.5, 1, 1)
  r <- sellab(y, rel, rule, w)
  expect_equal(r[["obs"]], (1 + 1) / (1 + 2 + 1 + 1 + 1), tolerance = 1e-15)
  # enumerate the outcomes of the two detained defendants the rule releases
  rates <- apply(expand.grid(0:1, 0:1), 1, function(f) {
    yy <- y
    yy[5:6] <- f
    sum((w * yy)[rule]) / sum(w[rule])
  })
  expect_equal(r[["lower"]], min(rates), tolerance = 1e-15)
  expect_equal(r[["upper"]], max(rates), tolerance = 1e-15)
  expect_equal(r[["unobs"]], 2 / 6, tolerance = 1e-15)
  # Csls.Select.unobserved_width
  expect_equal(r[["width"]], r[["unobs"]], tolerance = 1e-15)
  expect_equal(r[["naive"]], 1 / 4, tolerance = 1e-15)
  # a rule inside the judge's releases is identified
  r <- sellab(y, rel, c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE))
  expect_equal(r[["width"]], 0)
  expect_equal(r[["lower"]], 1 / 3)
  expect_equal(r[["naive"]], r[["lower"]])
  # the naive score sits at an end of the interval when the observed part
  # never fails or always fails
  rel <- c(TRUE, TRUE, FALSE, FALSE)
  rule <- c(TRUE, FALSE, TRUE, FALSE)
  r <- sellab(c(0, 0, NA, NA), rel, rule)
  expect_equal(unname(r[c("lower", "upper", "naive")]), c(0, 0.5, 0))
  r <- sellab(c(1, 1, NA, NA), rel, rule)
  expect_equal(unname(r[c("lower", "upper", "naive")]), c(0.5, 1, 1))
  # a rule that releases only detained defendants
  r <- sellab(c(1, NA), c(TRUE, FALSE), c(FALSE, TRUE))
  expect_equal(unname(r[c("lower", "upper")]), c(0, 1))
  expect_true(is.nan(r[["naive"]]))
  # Csls.Select.naive_between: the naive score never leaves the interval
  set.seed(9)
  for (k in 1:300) {
    n <- sample(2:40, 1)
    rel <- runif(n) < 0.6
    rule <- runif(n) < 0.6
    if (!any(rel) || !any(rule)) next
    y <- ifelse(rel, rbinom(n, 1, 0.4), NA)
    w <- runif(n)
    r <- sellab(y, rel, rule, w)
    if (is.nan(r[["naive"]])) next
    expect_gte(r[["naive"]], r[["lower"]] - 1e-12)
    expect_lte(r[["naive"]], r[["upper"]] + 1e-12)
  }
})

test_that("losbnd: the bound by hand, and never above a completed cohort", {
  r <- losbnd(c(30, 45, 60, 90, 120), c(100, 150, 200))
  expect_equal(r[["mean"]], 69)
  expect_equal(r[["age"]], 150)
  expect_equal(r[["lower"]], (345 + 450) / 8)
  # Csls.Censor.lower_bound_sub
  expect_equal(r[["bias"]], r[["lower"]] - r[["mean"]], tolerance = 1e-12)
  set.seed(8)
  for (k in 1:200) {
    t <- rexp(sample(1:30, 1), 1 / 50)
    a <- rexp(sample(1:30, 1), 1 / 80)
    r <- losbnd(t, a)
    # any completion u >= a has a mean at least the bound
    u <- a + rexp(length(a), 1 / 40)
    tm <- mean(c(t, u))
    expect_gte(tm, r[["lower"]] - 1e-12)
    expect_gte(tm - r[["mean"]], r[["bias"]] - 1e-12)
  }
  # open spells younger than the completed mean: the bound lies below it
  expect_lt(losbnd(c(100, 200), 10)[["bias"]], 0)
})

test_that("sntbnd bounds are the extremes over every counterfactual filling", {
  y <- c(1, 0, 1, 1, 0, 0)
  z <- c("a", "a", "b", "b", "b", "a")
  w <- c(1, 2, 1, 0.5, 1, 1.5)
  r <- sntbnd(y, z, w = w)
  expect_equal(r[["term"]], c("a", "b", "effect"))
  W <- sum(w)
  expect_equal(r[["share"]][1:2], c(4.5, 2.5) / W, tolerance = 1e-15)
  # brute force: fill each person's missing outcome with 0 or 1
  fills <- as.matrix(expand.grid(rep(list(0:1), 6)))
  eff <- apply(fills, 1, function(f) {
    y1 <- ifelse(z == "b", y, f)
    y0 <- ifelse(z == "a", y, f)
    c(sum(w * y0) / W, sum(w * y1) / W, sum(w * (y1 - y0)) / W)
  })
  expect_equal(r[["lower"]], apply(eff, 1, min), tolerance = 1e-14)
  expect_equal(r[["upper"]], apply(eff, 1, max), tolerance = 1e-14)
  # Csls.Sentence.ate_width_one, ate_contains_zero
  expect_equal(r[["upper"]][3] - r[["lower"]][3], 1, tolerance = 1e-15)
  expect_true(r[["lower"]][3] <= 0 && r[["upper"]][3] >= 0)
  expect_equal(r[["observed"]], c(1 / 4.5, 1.5 / 2.5, 1.5 / 2.5 - 1 / 4.5),
               tolerance = 1e-15)
  # the treatment picked by name flips the effect
  s <- sntbnd(y, z, treat = "a", w = w)
  expect_equal(s[["term"]], c("b", "a", "effect"))
  expect_equal(s[["lower"]][3], -r[["upper"]][3], tolerance = 1e-15)
  # an arm with zero weight leaves its observed rate undefined
  e <- sntbnd(c(1, 0), c("a", "b"), w = c(1, 0))
  expect_true(is.nan(e[["observed"]][3]))
})

test_that("hzdsel matches survivors counted by hand", {
  r <- hzdsel(s = 0.5, h = 0.5, l = 0.1, surh = c(0.5, 0.8), surl = c(0.9, 0.9))
  expect_equal(r[["share0"]], 5 / 14, tolerance = 1e-15)
  expect_equal(r[["share1"]], 8 / 17, tolerance = 1e-15)
  expect_equal(r[["haz0"]], 17 / 70, tolerance = 1e-15)
  expect_equal(r[["haz1"]], 49 / 170, tolerance = 1e-15)
  expect_equal(r[["ratio"]], 343 / 289, tolerance = 1e-14)
  # a cohort of 10,000: survivors and their reconvictions
  set.seed(6)
  for (k in 1:100) {
    s <- runif(1)
    h <- runif(1, 0.2, 1)
    l <- runif(1, 0, h)
    sh <- runif(2, 0.05, 1)
    sl <- runif(2, 0.05, 1)
    r <- hzdsel(s, h, l, sh, sl)
    hi <- 1e4 * s * sh
    lo <- 1e4 * (1 - s) * sl
    hz <- (hi * h + lo * l) / (hi + lo)
    expect_equal(unname(r[c("haz0", "haz1")]), hz, tolerance = 1e-12)
    # Csls.Hazard.hr2_gt_one_of_depletion, same low-risk survival
    r2 <- hzdsel(s, h, l, sort(sh), c(sl[1], sl[1]))
    if (sh[1] != sh[2]) expect_gt(r2[["ratio"]], 1)
  }
})

test_that("corrections functions stop on bad input, naming the argument", {
  expect_error(incap(0, 0.1, 1), "`lam`")
  expect_error(incap(1, 0, 1), "`q`")
  expect_error(incap(1, 1.5, 1), "at most 1")
  expect_error(incap(1, 0.5, -1), "`s`")
  expect_error(incap(c(1, 2, 3), c(0.1, 0.2), 1), "`q` must have length 1 or 3")
  expect_error(incap(c(1, 2, 3), 0.1, c(1, 2)), "`s` must have length 1 or 3")
  expect_error(incap(c(1, 2), 0.1, 1, share = c(0.5, 0.6)), "sum to 1")
  expect_error(incap(c(1, 2), 0.1, 1, share = c(1.5, -0.5)), "`share`")
  expect_error(incapr(c(1, 2), 0, 1), "never increase")
  expect_error(incapr(c(2, 1), 0, 0), "`s`")
  expect_error(incapr(c(2, 1), -1, 1), "`t0`")
  expect_error(incapr(c(2, 1), 0.5, 1), "`t0`")
  expect_error(incapr(c(2, 1), 1, 1), "length at least")
  expect_error(incapr(c(2, 1), 0, 1, r = 2), "`r`")
  expect_error(incapr(c(2, NA), 0, 1), "`lam`")
  expect_error(sellab(c(1, 0), c(TRUE, NA), c(TRUE, TRUE)), "`rel`")
  expect_error(sellab(c(1, 0), c(TRUE, TRUE), c(1, 0)), "`rule`")
  expect_error(sellab(c(1, 2), c(TRUE, TRUE), c(TRUE, TRUE)), "0 or 1")
  expect_error(sellab(c(1, NA), c(TRUE, TRUE), c(TRUE, TRUE)), "0 or 1")
  expect_error(sellab("a", TRUE, TRUE), "`y`")
  tt <- c(TRUE, TRUE)
  ff <- c(FALSE, FALSE)
  expect_error(sellab(c(1, 0), tt, tt, w = c(1, -1)), "`w`")
  expect_error(sellab(c(1, 0), ff, tt), "released nobody")
  expect_error(sellab(c(1, 0), tt, ff), "releases nobody")
  expect_error(losbnd(numeric(0), 1), "`done`")
  expect_error(losbnd(1, -1), "`open`")
  expect_error(sntbnd(c(0, 2), c("a", "b")), "0 or 1")
  expect_error(sntbnd(c(0, 1), c("a", "a")), "exactly two")
  expect_error(sntbnd(c(0, 1), c("a", NA)), "`z`")
  expect_error(sntbnd(c(0, 1), c("a", "b"), treat = "c"), "`treat`")
  expect_error(sntbnd(c(0, 1), c("a", "b"), w = c(0, 0)), "positive total")
  expect_error(hzdsel(1, 0.5, 0.1, c(1, 1), c(1, 1)), "`s`")
  expect_error(hzdsel(0.5, 1.5, 0.1, c(1, 1), c(1, 1)), "`h`")
  expect_error(hzdsel(0.5, 0.5, 0.6, c(1, 1), c(1, 1)), "below `h`")
  expect_error(hzdsel(0.5, 0.5, 0.1, c(0, 1), c(1, 1)), "`surh`")
  expect_error(hzdsel(0.5, 0.5, 0.1, c(1, 1), c(1, 1.2)), "`surl`")
  expect_error(hzdsel(0.5, 0.5, 0.1, 1, c(1, 1)), "`surh`")
})

test_that("corrections results agree with rmorie", {
  skip_if_not_installed("rmorie")
  o <- rmorie::morie_incapacitation(c(2, 10), 0.1, 1, shares = c(0.8, 0.2))
  d <- incap(c(2, 10), 0.1, 1, share = c(0.8, 0.2))
  expect_equal(d[["rate"]][1:2], o[["rate"]], tolerance = 1e-14)
  expect_equal(d[["prev"]][1:2], o[["prevented_share"]], tolerance = 1e-14)
  expect_equal(d[["marg"]][1:2], o[["marginal_prevention"]], tolerance = 1e-14)
  agg <- attr(o, "aggregate")
  expect_equal(d[["prev"]][3], unname(agg[["prevented_share"]]),
               tolerance = 1e-14)

  lam <- c(12, 10, 8, 6, 5, 4, 3, 2, 2, 1)
  o <- rmorie::morie_incapacitation_career(lam, 2, 3, 0.25)
  d <- incapr(lam, 2, 3, 0.25)
  expect_equal(unname(d[c("prev", "net", "upper", "lower", "later")]),
               c(o[["prevented"]], o[["prevented_net"]], o[["upper"]],
                 o[["lower"]], o[["later"]]), tolerance = 1e-14)
  # rmorie returns NA for `later` when the path is exactly t0 + s + 1 long,
  # although the sum needs ages up to t0 + s only; csls computes it
  o <- rmorie::morie_incapacitation_career(c(5, 4, 3), 0, 2)
  expect_true(is.na(o[["later"]]))
  expect_equal(incapr(c(5, 4, 3), 0, 2)[["later"]], 7)

  set.seed(4)
  risk <- runif(500)
  y <- rbinom(500, 1, risk)
  o <- rmorie::morie_selective_labels(y, risk < 0.7, risk < 0.8)
  d <- sellab(y, risk < 0.7, risk < 0.8)
  expect_equal(d[["obs"]], o[["observed_rate"]], tolerance = 1e-14)
  expect_equal(d[["unobs"]], o[["rule_share_unobserved"]], tolerance = 1e-14)
  expect_equal(unname(d[c("lower", "upper")]), unname(o[["bounds"]]),
               tolerance = 1e-14)
  expect_equal(d[["naive"]], o[["naive_rate"]], tolerance = 1e-14)

  o <- rmorie::morie_backlog_censoring(c(30, 45, 60, 90, 120), c(100, 150, 200))
  d <- losbnd(c(30, 45, 60, 90, 120), c(100, 150, 200))
  expect_equal(d[["lower"]], o[["lower_bound"]], tolerance = 1e-14)
  expect_equal(d[["bias"]], o[["bias_lower"]], tolerance = 1e-14)

  set.seed(1)
  z <- sample(c("community", "custody"), 500, TRUE, c(0.6, 0.4))
  y <- rbinom(500, 1, ifelse(z == "custody", 0.55, 0.35))
  o <- rmorie::morie_sentence_effect_bounds(y, z)
  d <- sntbnd(y, z)
  expect_equal(d[["lower"]][3], o[["ate_bounds"]][["lower"]], tolerance = 1e-14)
  expect_equal(d[["upper"]][3], o[["ate_bounds"]][["upper"]], tolerance = 1e-14)
  expect_equal(d[["observed"]][3], o[["naive_difference"]], tolerance = 1e-14)

  o <- rmorie::morie_hazard_selection(0.5, 0.5, 0.1, c(0.5, 0.8), c(0.9, 0.9))
  d <- hzdsel(0.5, 0.5, 0.1, c(0.5, 0.8), c(0.9, 0.9))
  expect_equal(unname(d[c("haz0", "haz1")]), unname(o[["period2_hazard"]]),
               tolerance = 1e-14)
  expect_equal(d[["ratio"]], o[["period2_hazard_ratio"]], tolerance = 1e-14)
})
