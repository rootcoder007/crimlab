test_that("clrsmp decomposes the toy Simpson case by hand", {
  r <- clrsmp(c(minor = 90, serious = 20), c(100, 50), c(19, 60), c(20, 140))
  p <- r[["parts"]]
  tot <- r[["total"]]
  expect_identical(p[["part"]], c("minor", "serious"))
  expect_equal(p[["ratea"]], c(0.9, 0.4), tolerance = 1e-15)
  expect_equal(p[["rateb"]], c(0.95, 60 / 140), tolerance = 1e-15)
  expect_equal(p[["sharea"]], c(2 / 3, 1 / 3), tolerance = 1e-15)
  expect_equal(p[["shareb"]], c(1 / 8, 7 / 8), tolerance = 1e-15)
  expect_equal(unname(tot[c("ratea", "rateb")]), c(110 / 150, 79 / 160),
               tolerance = 1e-15)
  # by hand: rate terms (wa + wb) / 2 (rb - ra)
  re <- c((2 / 3 + 1 / 8) / 2 * 0.05, (1 / 3 + 7 / 8) / 2 * (60 / 140 - 0.4))
  expect_equal(p[["rateeff"]], re, tolerance = 1e-12)
  # mix terms centred on the mean of the aggregates
  cen <- (110 / 150 + 79 / 160) / 2
  me <- ((c(0.9, 0.4) + c(0.95, 60 / 140)) / 2 - cen) *
    (c(1 / 8, 7 / 8) - c(2 / 3, 1 / 3))
  expect_equal(p[["mixeff"]], me, tolerance = 1e-12)
  expect_equal(unname(tot["change"]), 79 / 160 - 110 / 150, tolerance = 1e-15)
  expect_equal(unname(tot[["up"]]), 2, tolerance = 1e-15)
  expect_equal(unname(tot[["down"]]), 0, tolerance = 1e-15)
  expect_equal(unname(tot[["simpson"]]), 2, tolerance = 1e-15)
  # Csls.Clearance.change_ge_mix: every type up, so the fall is all mix
  expect_gt(tot[["rateeff"]], 0)
  expect_lt(tot[["mixeff"]], tot[["change"]] + 1e-15)
})

test_that("clrsmp terms add up and the mix total ignores the centre", {
  set.seed(11)
  for (k in 1:200) {
    n <- sample(2:6, 1)
    da <- runif(n, 1, 500)
    db <- runif(n, 1, 500)
    na <- da * runif(n)
    nb <- db * runif(n)
    r <- clrsmp(na, da, nb, db)
    tot <- r[["total"]]
    p <- r[["parts"]]
    chg <- sum(nb) / sum(db) - sum(na) / sum(da)
    expect_equal(unname(tot[["change"]]), chg, tolerance = 1e-12)
    expect_equal(tot[["rateeff"]] + tot[["mixeff"]], chg, tolerance = 1e-12)
    expect_equal(tot[["rateeff"]], sum(p[["rateeff"]]), tolerance = 1e-12)
    # the uncentred mix effect has the same total (Csls.Clearance.kitagawa)
    unc <- sum((na / da + nb / db) / 2 * (db / sum(db) - da / sum(da)))
    expect_equal(tot[["mixeff"]], unc, tolerance = 1e-12)
    # Csls.Clearance.agg_eq_weighted
    expect_equal(tot[["ratea"]], sum(da / sum(da) * na / da), tolerance = 1e-12)
    # Csls.Clearance.change_ge_mix whenever no rate falls
    if (all(nb / db >= na / da)) {
      expect_gte(tot[["change"]] - tot[["mixeff"]], -1e-12)
    }
  }
})

test_that("clrsmp flags the three kinds", {
  # every rate down, aggregate up
  r <- clrsmp(c(19, 60), c(20, 140), c(90, 20), c(100, 50))[["total"]]
  expect_equal(unname(r[c("up", "down", "simpson")]), c(0, 2, 2),
               tolerance = 1e-15)
  # one part up, one down, the aggregate and the rate effect both fall
  r <- clrsmp(c(80, 10), c(100, 100), c(85 * 0.5, 0.05 * 150),
              c(50, 150))[["total"]]
  re <- (0.5 + 0.25) / 2 * (0.85 - 0.8) + (0.5 + 0.75) / 2 * (0.05 - 0.1)
  expect_equal(r[["rateeff"]], re, tolerance = 1e-12)
  expect_equal(unname(r[c("up", "down", "simpson")]), c(1, 1, 0),
               tolerance = 1e-15)
  r <- clrsmp(c(80, 10), c(100, 100), c(85, 9.6), c(100, 100))[["total"]]
  expect_equal(unname(r[["simpson"]]), 0, tolerance = 1e-15)
  expect_equal(r[["mixeff"]], 0, tolerance = 1e-15)
  # aggregate falls, rate effect rises, not every part up: kind 1
  r <- clrsmp(c(80, 10), c(100, 100), c(45, 7.2), c(50, 150))[["total"]]
  expect_gt(r[["rateeff"]], 0)
  expect_lt(r[["change"]], 0)
  expect_equal(unname(r[c("up", "down", "simpson")]), c(1, 1, 1),
               tolerance = 1e-15)
  # no change in the counts: nothing moves
  r <- clrsmp(c(5, 6), c(10, 10), c(5, 6), c(10, 10))[["total"]]
  expect_equal(unname(r[c("change", "rateeff", "mixeff", "simpson")]),
               c(0, 0, 0, 0), tolerance = 1e-15)
  # rates above one (clearances of earlier years' offences) are allowed
  expect_equal(clrsmp(12, 10, 5, 10)[["total"]][["ratea"]], 1.2)
})

test_that("expbnd bounds hold for every admissible exposure and are reached", {
  r <- expbnd(c(white = 100, black = 300), c(1000, 1000), gamma = 1.5)
  expect_equal(r[["rate"]], c(0.1, 0.3))
  expect_equal(r[["lo"]], c(0.1, 0.3) / 1.5, tolerance = 1e-15)
  expect_equal(r[["hi"]], c(0.1, 0.3) * 1.5, tolerance = 1e-15)
  expect_equal(r[["ratio"]], c(1, 3), tolerance = 1e-15)
  expect_equal(r[["rlo"]], c(1, 3) / 2.25, tolerance = 1e-15)
  expect_equal(r[["rhi"]], c(1, 3) * 2.25, tolerance = 1e-15)
  expect_identical(r[["sign"]], c(0L, 1L))
  # ref by name; a group below the band
  r <- expbnd(c(a = 100, b = 300), c(1000, 1000), gamma = 1.5, ref = "b")
  expect_identical(r[["sign"]], c(-1L, 0L))
  # the band closes at gamma = 1, and a wide band leaves the sign open
  expect_identical(expbnd(c(2, 3), c(1, 1), 1)[["sign"]], c(0L, 1L))
  expect_identical(expbnd(c(1, 3), c(1, 1), 2)[["sign"]], c(0L, 0L))
  # random exposures inside the band give ratios inside [rlo, rhi]; the
  # ends come from putting the two groups at opposite ends of their bands
  set.seed(3)
  y <- c(40, 90, 15)
  m <- c(200, 250, 90)
  g <- 1.7
  b <- expbnd(y, m, g, ref = 2)
  for (k in 1:500) {
    e <- m * exp(runif(3, -log(g), log(g)))
    tr <- (y / e) / (y[2] / e[2])
    expect_true(all(tr >= b[["rlo"]] * (1 - 1e-12)))
    expect_true(all(tr <= b[["rhi"]] * (1 + 1e-12)))
    tr1 <- y / e
    expect_true(all(tr1 >= b[["lo"]] * (1 - 1e-12)))
    expect_true(all(tr1 <= b[["hi"]] * (1 + 1e-12)))
  }
  lo1 <- (y[1] / (g * m[1])) / (y[2] / (m[2] / g))
  hi1 <- (y[1] / (m[1] / g)) / (y[2] / (g * m[2]))
  expect_equal(b[["rlo"]][1], lo1, tolerance = 1e-14)
  expect_equal(b[["rhi"]][1], hi1, tolerance = 1e-14)
})

test_that("misrd matches the 2 x 2 tables enumerated cell by cell", {
  cells <- function(r1, r0, p, se, sp) {
    # P(A = a, recorded = s), and the outcome risk given A
    pa <- c(p, 1 - p)
    rec1 <- c(se, 1 - sp)
    risk <- c(r1, r0)
    m1 <- sum(pa * rec1 * risk) / sum(pa * rec1)
    m0 <- sum(pa * (1 - rec1) * risk) / sum(pa * (1 - rec1))
    m1 - m0
  }
  g <- expand.grid(r1 = c(0, 0.3, 0.9), r0 = c(0.1, 0.5), p = c(0.05, 0.4, 0.9),
                   se = c(0, 0.3, 0.8, 1), sp = c(0.2, 0.7, 1))
  for (i in seq_len(nrow(g))) {
    r1 <- g[["r1"]][i]
    r0 <- g[["r0"]][i]
    p <- g[["p"]][i]
    se <- g[["se"]][i]
    sp <- g[["sp"]][i]
    d1 <- p * se + (1 - p) * (1 - sp)
    if (d1 <= 0 || d1 >= 1) next
    x <- misrd(r1, r0, p, se, sp)
    expect_equal(x[["measured"]], cells(r1, r0, p, se, sp), tolerance = 1e-12)
    # Csls.Misclass.binary_attenuation
    expect_equal(x[["measured"]], x[["factor"]] * (r1 - r0),
                 tolerance = 1e-12)
    expect_equal(x[["d1"]] + x[["d0"]], 1, tolerance = 1e-15)
    if (se + sp >= 1) {
      # never reversed (binary_sign) and never larger: the documented claim
      expect_gte(x[["factor"]], -1e-15)
      expect_lte(x[["factor"]], 1 + 1e-12)
    } else {
      expect_lt(x[["factor"]], 0)
    }
  }
  # perfect measurement changes nothing
  expect_equal(misrd(0.3, 0.1, 0.4, 1, 1)[["measured"]], 0.2, tolerance = 1e-15)
})

test_that("rtmhs equals the statistic on the data with its swapped copy", {
  x1 <- c(10, 2, 7, 4, 9)
  x2 <- c(6, 3, 8, 5, 2)
  r <- rtmhs(x1, x2, cut = 6)
  expect_equal(r[["nsel"]], 3, tolerance = 1e-15)
  expect_equal(r[["change"]], mean(c(6 - 10, 8 - 7, 2 - 9)), tolerance = 1e-15)
  expect_equal(r[["mirror"]], 7 - 8, tolerance = 1e-15)
  # selection on x1 > cut in the stacked data and its period-swapped copy
  a <- c(x1, x2)
  b <- c(x2, x1)
  expect_equal(r[["sym"]], mean((b - a)[a > 6]), tolerance = 1e-15)
  expect_equal(r[["lowsym"]], mean((b - a)[a < 6]), tolerance = 1e-15)
  expect_equal(r[["low"]], mean(c(3 - 2, 5 - 4)), tolerance = 1e-15)
  expect_equal(r[["excess"]], r[["change"]] - r[["sym"]], tolerance = 1e-15)
  # weights: a place of weight 2 counts as two copies
  rw <- rtmhs(x1, x2, 6, w = c(2, 1, 1, 1, 1))
  rr <- rtmhs(c(x1, 10), c(x2, 6), 6)
  expect_equal(unname(rw[-1]), unname(rr[-1]), tolerance = 1e-14)
  # nothing above the cut in period 2, nothing below in either
  expect_true(is.nan(rtmhs(c(5, 1), c(1, 1), 4)[["mirror"]]))
  expect_true(is.nan(rtmhs(c(5, 6), c(7, 8), 4)[["low"]]))
})

test_that("rtmhs: Csls.RegMean signs hold on any data", {
  set.seed(19)
  for (k in 1:300) {
    n <- sample(5:60, 1)
    x1 <- rpois(n, 5) * runif(1, 0.5, 3)
    x2 <- rpois(n, 5) + sample(-3:6, 1)
    w <- runif(n)
    cut <- quantile(x1, runif(1, 0.1, 0.8), names = FALSE)
    if (sum(w[x1 > cut]) <= 0) next
    r <- rtmhs(x1, x2, cut, w)
    expect_lte(r[["sym"]], 1e-12)
    if (!is.nan(r[["lowsym"]])) expect_gte(r[["lowsym"]], -1e-12)
  }
  # under exchangeable periods with no treatment the selected places fall:
  # 400 replications, each the mean over 2,000 places; the mean change is
  # several standard errors below 0
  set.seed(7)
  ch <- replicate(400, {
    mu <- rgamma(2000, 2, 0.5)
    x1 <- rpois(2000, mu)
    x2 <- rpois(2000, mu)
    rtmhs(x1, x2, quantile(x1, 0.9, names = FALSE))[["change"]]
  })
  expect_lt(mean(ch) + 5 * sd(ch) / sqrt(400), 0)
})

test_that("policing functions stop on bad input, naming the argument", {
  expect_error(clrsmp(c(1, 1), c(1, 0), c(1, 1), c(1, 1)), "`reca`")
  expect_error(clrsmp(c(1, -1), c(1, 1), c(1, 1), c(1, 1)), "`clra`")
  expect_error(clrsmp(c(1, 1), c(1, 1), c(1, NA), c(1, 1)), "`clrb`")
  expect_error(clrsmp(c(1, 1), c(1, 1), c(1, 1), c(1, 1, 1)), "length 2")
  expect_error(clrsmp(c(1, 1), c(1, 1), c(1, 1), "a"), "`recb`")
  expect_error(expbnd(c(1, 0), c(1, 1), 1), "`y`")
  expect_error(expbnd(c(1, 1), c(1, -1), 1), "`m`")
  expect_error(expbnd(c(1, 1), c(1, 1), 0.9), "`gamma`")
  expect_error(expbnd(c(1, 1), c(1, 1), c(1, 2)), "`gamma`")
  expect_error(expbnd(c(1, 1), c(1, 1), 1, ref = 3), "`ref`")
  expect_error(misrd(1.1, 0.1, 0.5, 0.9, 0.9), "`r1`")
  expect_error(misrd(0.1, NA, 0.5, 0.9, 0.9), "`r0`")
  expect_error(misrd(0.1, 0.1, 1, 0.9, 0.9), "`prev`")
  expect_error(misrd(0.1, 0.1, 0.5, -0.1, 0.9), "`se`")
  expect_error(misrd(0.1, 0.1, 0.5, 0.9, c(1, 1)), "`sp`")
  expect_error(misrd(0.1, 0.1, 0.5, 0, 1), "recorded group empty")
  expect_error(misrd(0.1, 0.1, 0.5, 1, 0), "recorded group empty")
  expect_error(rtmhs(c(1, 2), c(1, 2, 3), 0), "length 2")
  expect_error(rtmhs(c(1, Inf), c(1, 2), 0), "`x1`")
  expect_error(rtmhs(c(1, 2), c(1, 2), NA), "`cut`")
  expect_error(rtmhs(c(1, 2), c(1, 2), 0, w = c(1, -1)), "`w`")
  expect_error(rtmhs(c(1, 2), c(1, 2), 5), "exceeds `cut`")
  expect_error(rtmhs(c(1, 2), c(1, 2), 1, w = c(1, 0)), "exceeds `cut`")
})

test_that("policing results agree with rmorie", {
  skip_if_not_installed("rmorie")
  y <- c(A = 300, B = 100, C = 50)
  m <- c(A = 1000, B = 1000, C = 400)
  o <- rmorie::morie_disparity_exposure_bounds(y, m, gamma = 1.5,
                                               reference = "B")
  d <- expbnd(y, m, 1.5, ref = "B")
  expect_equal(d[["rate"]], o[["proxy_rate"]], tolerance = 1e-14)
  expect_equal(d[["lo"]], o[["rate_lower"]], tolerance = 1e-14)
  expect_equal(d[["hi"]], o[["rate_upper"]], tolerance = 1e-14)
  expect_equal(d[["rlo"]], o[["ratio_lower"]], tolerance = 1e-14)
  expect_equal(d[["rhi"]], o[["ratio_upper"]], tolerance = 1e-14)
  expect_identical(d[["sign"]] != 0, o[["direction_identified"]])

  o <- rmorie::morie_misclassified_rd(0.3, 0.1, 0.4, 0.8, 0.9)
  d <- misrd(0.3, 0.1, 0.4, 0.8, 0.9)
  expect_equal(d[["measured"]], o[["measured_rd"]], tolerance = 1e-14)
  expect_equal(d[["factor"]], o[["factor"]], tolerance = 1e-14)

  set.seed(5)
  mu <- rgamma(300, 2, 0.5)
  x1 <- rpois(300, mu)
  x2 <- rpois(300, mu)
  cut <- quantile(x1, 0.8, names = FALSE)
  o <- rmorie::morie_regression_to_mean(x1, x2, cut)
  d <- rtmhs(x1, x2, cut)
  expect_equal(d[["change"]], o[["selected_change"]], tolerance = 1e-14)
  expect_equal(d[["mirror"]], o[["mirror_change"]], tolerance = 1e-14)
  expect_equal(d[["sym"]], o[["symmetrised_change"]], tolerance = 1e-14)
  expect_equal(d[["low"]], o[["low_selected_change"]], tolerance = 1e-14)
  expect_equal(d[["lowsym"]], o[["low_symmetrised_change"]], tolerance = 1e-14)
})
