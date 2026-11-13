# grprr(), dastd(), sirex(): intervals against stats::poisson.test(),
# binom.test() and prop.test(), Newcombe's and Fay and Feuer's formulas by
# hand, Marshall's shrinkage by hand, zero and huge counts, the input forms
# and the refusals.

test_that("rates: exact and approximate intervals against poisson.test()", {
  r <- grprr(
    c(a = 12, b = 30, c = 0),
    c(a = 4000, b = 5000, c = 800),
    per = 1000
  )
  q <- r$rates
  expect_identical(q$group, c("a", "b", "c"))
  for (i in 1:3) {
    pt <- poisson.test(q$x[i], q$n[i])$conf.int * 1000
    expect_equal(c(q$elo[i], q$ehi[i]), as.numeric(pt), tolerance = 1e-12)
  }
  expect_equal(q$lo[1], 3 * exp(-qnorm(0.975) / sqrt(12)), tolerance = 1e-12)
  expect_true(is.na(q$lo[3]))
  # the ratio: conditional on the total, as the two-sample poisson.test()
  m <- r$cmp
  pt <- poisson.test(c(30, 12), c(5000, 4000))
  expect_equal(m$ratio[1], unname(pt$estimate), tolerance = 1e-12)
  expect_equal(
    c(m$elo[1], m$ehi[1]),
    as.numeric(pt$conf.int),
    tolerance = 1e-10
  )
  s <- sqrt(1 / 30 + 1 / 12)
  expect_equal(m$hi[1], 2 * exp(qnorm(0.975) * s), tolerance = 1e-12)
  d <- 1000 * qnorm(0.975) * sqrt(30 / 5000^2 + 12 / 4000^2)
  expect_equal(m$dlo[1], 3 - d, tolerance = 1e-12)
  # a zero count: ratio 0, its exact interval from 0
  expect_identical(m$ratio[2], 0)
  expect_identical(m$elo[2], 0)
  expect_true(is.na(m$lo[2]))
  # representation: share of events over share of exposure
  expect_equal(q$rep[2], (30 / 42) / (5000 / 9800), tolerance = 1e-14)
})

test_that("proportions: Wilson, Clopper-Pearson, Katz and Newcombe", {
  r <- grprr(c(15, 40), c(60, 80), type = "prop", ref = 2)
  q <- r$rates
  for (i in 1:2) {
    w <- prop.test(q$x[i], q$n[i], correct = FALSE)$conf.int
    expect_equal(c(q$lo[i], q$hi[i]), as.numeric(w), tolerance = 1e-12)
    b <- binom.test(q$x[i], q$n[i])$conf.int
    expect_equal(c(q$elo[i], q$ehi[i]), as.numeric(b), tolerance = 1e-12)
  }
  m <- r$cmp
  expect_identical(r$ref, "2")
  expect_equal(m$ratio, 0.25 / 0.5, tolerance = 1e-14)
  s <- sqrt(1 / 15 - 1 / 60 + 1 / 40 - 1 / 80)
  expect_equal(m$lo, 0.5 * exp(-qnorm(0.975) * s), tolerance = 1e-12)
  expect_true(is.na(m$elo))
  # Newcombe (1998), method 10, by hand from the Wilson limits
  l1 <- q$lo[1]
  u1 <- q$hi[1]
  l2 <- q$lo[2]
  u2 <- q$hi[2]
  expect_equal(
    m$dlo,
    -0.25 - sqrt((0.25 - l1)^2 + (u2 - 0.5)^2),
    tolerance = 1e-12
  )
  expect_equal(
    m$dhi,
    -0.25 + sqrt((u1 - 0.25)^2 + (0.5 - l2)^2),
    tolerance = 1e-12
  )
  # one row per case gives the same
  y <- c(rep(1, 15), rep(0, 45), rep(1, 40), rep(0, 40))
  g <- rep(1:2, c(60, 80))
  r2 <- grprr(y, grp = g, type = "prop", ref = 2)
  expect_identical(r2$cmp, r$cmp)
  expect_output(print(r2), "proportions")
})

test_that("swapping the reference inverts the ratio", {
  a <- grprr(c(5, 9), c(100, 120))
  b <- grprr(c(5, 9), c(100, 120), ref = 2)
  expect_equal(b$cmp$ratio, 1 / a$cmp$ratio, tolerance = 1e-14)
  expect_equal(b$cmp$lo, 1 / a$cmp$hi, tolerance = 1e-12)
  expect_equal(b$cmp$elo, 1 / a$cmp$ehi, tolerance = 1e-10)
  expect_equal(b$cmp$diff, -a$cmp$diff, tolerance = 1e-14)
})

test_that("huge counts and one group", {
  r <- grprr(c(3e9, 1e9), c(1e12, 1e12))
  expect_equal(r$cmp$ratio, 1 / 3, tolerance = 1e-14)
  expect_lt(r$cmp$hi - r$cmp$lo, 1e-4)
  one <- grprr(5, 10)
  expect_identical(nrow(one$cmp), 0L)
  expect_output(print(one), "grprr")
})

test_that("direct standardisation: Fay and Feuer by hand", {
  d <- data.frame(
    age = rep(c("y", "m", "o"), 2),
    grp = rep(c("A", "B"), each = 3),
    x = c(30, 22, 5, 41, 30, 12),
    n = c(5000, 9000, 6000, 3000, 8000, 7000)
  )
  std <- c(y = 30, m = 45, o = 25)
  r <- dastd(d$x, d$n, d$age, d$grp, std = std, per = 1e4)
  w <- std[d$age] / 100 / d$n
  for (g in c("A", "B")) {
    i <- d$grp == g
    y <- sum(w[i] * d$x[i])
    v <- sum(w[i]^2 * d$x[i])
    wm <- max(w[i])
    k <- r$rates$group == g
    expect_equal(r$rates$rate[k], 1e4 * y, tolerance = 1e-12)
    expect_equal(
      r$rates$lo[k],
      1e4 * v / (2 * y) * qchisq(0.025, 2 * y^2 / v),
      tolerance = 1e-10
    )
    expect_equal(
      r$rates$hi[k],
      1e4 *
        (v + wm^2) /
        (2 * (y + wm)) *
        qchisq(0.975, 2 * (y + wm)^2 / (v + wm^2)),
      tolerance = 1e-10
    )
  }
  # unnamed weights follow the strata's sorted labels: m, o, y
  r2 <- dastd(d$x, d$n, d$age, d$grp, std = c(45, 25, 30), per = 1e4)
  expect_identical(r2$rates, r$rates)
  a <- r$rates$rate
  expect_equal(r$cmp$ratio, a[2] / a[1], tolerance = 1e-14)
  # a standardised rate lies between its smallest and largest stratum rate
  sr <- d$x / d$n * 1e4
  expect_gte(a[1], min(sr[1:3]))
  expect_lte(a[1], max(sr[1:3]))
  expect_output(print(r), "against A")
  # a group with no events: rate 0, ratio 0 with no log-normal interval
  z <- dastd(c(0, 0, 0, 3, 4, 5), d$n, d$age, d$grp, std = std)
  expect_identical(z$cmp$ratio, Inf)
  expect_true(is.na(z$cmp$lo))
  p <- structure(list(x = d$x, n = d$n, strat = d$age, grp = d$grp, std = std),
                 class = "clpre")
  expect_identical(dastd(p, per = 1e4)$rates, r$rates)
  one <- dastd(d$x[1:3], d$n[1:3], d$age[1:3], std = std)
  expect_identical(one$rates$group, "all")
})

test_that("indirect standardisation and shrinkage by hand", {
  d <- data.frame(
    area = rep(c("N", "S", "E", "W"), each = 2),
    age = rep(c("young", "old"), 4),
    x = c(12, 3, 30, 9, 4, 1, 20, 2),
    n = c(4000, 2000, 8000, 3000, 1500, 900, 3000, 2500)
  )
  r <- sirex(d$x, d$n, d$age, d$area)
  s <- r$sir
  rt <- tapply(d$x, d$age, sum) / tapply(d$n, d$age, sum)
  e <- tapply(d$n * rt[d$age], d$area, sum)
  o <- tapply(d$x, d$area, sum)
  expect_equal(s$expd, as.numeric(e[s$group]), tolerance = 1e-12)
  expect_equal(s$sir, as.numeric((o / e)[s$group]), tolerance = 1e-12)
  i <- 1
  expect_equal(
    c(s$lo[i], s$hi[i]),
    as.numeric(poisson.test(s$obs[i], s$expd[i])$conf.int),
    tolerance = 1e-12
  )
  expect_equal(s$flo[i], qpois(0.025, s$expd[i]) / s$expd[i], tolerance = 1e-14)
  # Marshall's moment estimates
  m <- sum(o) / sum(e)
  v <- sum(e * (o / e - m)^2) / sum(e) - m / mean(e)
  expect_equal(r$prior[["mean"]], m, tolerance = 1e-12)
  if (v > 0) {
    eb <- (o + m^2 / v) / (e + m / v)
    expect_equal(s$eb, as.numeric(eb[s$group]), tolerance = 1e-12)
  }
  # no extra variation: every ratio shrunk to the mean
  flat <- sirex(c(10, 10), c(1000, 1000), c(1, 1), c("a", "b"))
  expect_identical(flat$sir$shr, c(1, 1))
  expect_identical(flat$prior[["var"]], 0)
  # given stratum rates
  r2 <- sirex(d$x, d$n, d$age, d$area, rates = c(young = 0.004, old = 0.002))
  expect_equal(
    r2$sir$expd[r2$sir$group == "N"],
    4000 * 0.004 + 2000 * 0.002,
    tolerance = 1e-14
  )
  expect_output(print(r), "funnel")
  p <- structure(list(x = d$x, n = d$n, strat = d$age, grp = d$area),
                 class = "clpre")
  expect_identical(sirex(p)$sir, r$sir)
})

test_that("refusals", {
  errs(
    "grprr",
    list(
      list("`x` must be at least 0", c(-1, 2), c(1, 2)),
      list("`n` must be greater than 0", c(1, 2), c(0, 2)),
      list("`n` must give the exposure", c(1, 2)),
      list("more events than trials", c(5, 2), c(3, 4), type = "prop"),
      list("`type` must be one of", c(1, 2), c(3, 4), type = "odds"),
      list("`per` must be one number", c(1, 2), c(3, 4), per = 0),
      list("`grp` must be labels", c(1, 2), c(3, 4), c("a", NA)),
      list("`ref` names no group", c(1, 2), c(3, 4), ref = "q")
    )
  )
  errs(
    "dastd",
    list(
      list("one row per stratum", c(1, 2), c(1, 1), c("a", "a"), std = 1),
      list(
        "every group needs a row",
        c(1, 2),
        c(1, 1),
        c("a", "b"),
        c("g", "h"),
        std = c(1, 1)
      ),
      list(
        "`std` must have one value per stratum",
        c(1, 2),
        c(1, 1),
        c("a", "b"),
        std = 1
      ),
      list(
        "`std` has no value for stratum b",
        c(1, 2),
        c(1, 1),
        c("a", "b"),
        std = c(a = 1, c = 1)
      ),
      list(
        "`std` names stratum a twice",
        c(1, 2),
        c(1, 1),
        c("a", "b"),
        std = stats::setNames(c(1, 1, 2), c("a", "a", "b"))
      ),
      list("positive sum", c(1, 2), c(1, 1), c("a", "b"), std = c(0, 0))
    )
  )
  errs(
    "sirex",
    list(
      list("`grp` must give the group", c(1, 2), c(1, 1), c("a", "b"), NULL),
      list("positive expected count", c(0, 0), c(1, 1), c("a", "b"), c(1, 2))
    )
  )
})

# Parity with the source functions (helper-port.R): live when rmorie and
# rmoriebricklayer are installed, recorded values otherwise.

test_that("dastd() reproduces agestd() and rategap()", {
  d <- pin$dsr
  b <- portref("agestd")
  r <- dastd(d$x, d$n, d$age, d$grp, std = c(y = 30, m = 45, o = 25), per = 1e4)
  expect_equal(r$rates$rate, b$rate, tolerance = 1e-12)
  expect_equal(r$rates$lo, b$lower, tolerance = 1e-10)
  expect_equal(r$rates$hi, b$upper, tolerance = 1e-10)
  expect_equal(
    unlist(
      r$cmp[c("ratio", "lo", "hi", "diff", "dlo", "dhi")],
      use.names = FALSE
    ),
    unlist(
      b[c("ratio", "ratiolo", "ratioup", "diff", "diflo", "difup")],
      use.names = FALSE
    ),
    tolerance = 1e-10
  )
})

test_that("sirex() reproduces stdinc(), ebrates(), fnllim() and expcnts()", {
  b <- portref("stdinc")
  r <- sirex(pin$obs, pin$expd, rep(1, 4), 1:4, rates = 1)$sir
  expect_equal(r$lo, b$lower, tolerance = 1e-12)
  expect_equal(r$hi, b$upper, tolerance = 1e-12)
  expect_equal(r$eb, b$eb, tolerance = 1e-12)
  expect_equal(r$shr, b$shrinkage, tolerance = 1e-12)
  # fnllim() sorts by the expected count
  o <- order(pin$expd)
  expect_equal(r$flo[o], b$flo, tolerance = 1e-14)
  expect_equal(r$fhi[o], b$fhi, tolerance = 1e-14)
  d <- pin$dsr
  expect_equal(
    sirex(d$x, d$n, d$age, d$grp)$sir$expd,
    portref("expcnts"),
    tolerance = 1e-12
  )
})

test_that("grprr() reproduces rtcalc() and ratechg()", {
  b <- portref("rtcalc")
  x <- pin$chg
  r <- grprr(sum(x$count), sum(x$pop), per = 1000)$rates
  expect_equal(
    c(r$rate, r$elo, r$ehi),
    c(b$rate, b$lower, b$upper),
    tolerance = 1e-12
  )
  k <- grprr(x$count, x$pop, x$year)$cmp
  expect_equal(k$ratio, b$rr, tolerance = 1e-14)
  expect_equal(
    100 * (c(k$elo, k$ehi) - 1),
    c(b$pctlo, b$pctup),
    tolerance = 1e-10
  )
})

test_that("grprr() reproduces the risk ratio and difference of rmorie", {
  b <- portref("riskratio")
  t <- pin$tab
  r <- grprr(
    c(t[["a"]], t[["c"]]),
    c(t[["a"]] + t[["b"]], t[["c"]] + t[["d"]]),
    type = "prop",
    ref = 2
  )$cmp
  expect_equal(c(r$ratio, r$lo, r$hi), c(b$rr, b$lo, b$hi), tolerance = 1e-12)
  expect_equal(
    c(r$ratio, r$lo, r$hi),
    c(b$rr2, b$lo2, b$hi2),
    tolerance = 1e-12
  )
  expect_equal(
    c(r$diff, r$dlo, r$dhi),
    c(b$rd, b$rdlo, b$rdhi),
    tolerance = 1e-12
  )
  # disagreement: risk_difference() gives the Wald interval, which
  # under-covers near 0 and 1; crimlab gives Newcombe's
  expect_equal(r$diff, b$wrd, tolerance = 1e-14)
  wald <- -0.25 + c(-1, 1) * qnorm(0.975) * sqrt(0.25 * 0.75 / 60 + 0.25 / 80)
  expect_equal(c(b$wlo, b$whi), wald, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(r$dlo, b$wlo)))
  # disagreement: with a zero cell risk_ratio() adds 1/2 to every cell;
  # crimlab reports the ratio 0 and no log-normal interval
  z <- grprr(c(0, 5), c(10, 10), type = "prop", ref = 2)$cmp
  expect_identical(z$ratio, 0)
  expect_true(is.na(z$lo))
  expect_equal(b$zrr, (0.5 / 11) / (5.5 / 11), tolerance = 1e-14)
})

test_that("grprr() reproduces use-of-force and SIU rates", {
  b <- portref("uofdisp")
  u <- pin$uof
  # rmorie's default baseline is the largest group; crimlab's is the first
  r <- grprr(u$force, grp = u$race, type = "prop", ref = b$baseline)
  q <- r$rates[match(b$category, r$rates$group), ]
  expect_equal(q$rate, b$rate, tolerance = 1e-14)
  expect_equal(q$lo, b$lo, tolerance = 1e-12)
  expect_equal(q$hi, b$hi, tolerance = 1e-12)
  rr <- c(r$cmp$ratio, 1)[match(b$category, c(r$cmp$group, b$baseline))]
  expect_equal(rr, b$rr, tolerance = 1e-14)
  s <- portref("siusvc")
  v <- pin$siu
  g <- grprr(
    as.numeric(v$charges_recommended),
    grp = v$police_service,
    type = "prop"
  )
  expect_equal(
    g$rates$n[match(names(s$counts), g$rates$group)],
    as.numeric(s$counts),
    tolerance = 0
  )
  shown <- vapply(s$rows, function(e) e[[5]], "")
  ours <- sprintf(
    "%.1f%%",
    100 *
      g$rates$rate[match(
        vapply(s$rows, function(e) e[[1]], ""),
        g$rates$group
      )]
  )
  expect_identical(ours, shown)
  t <- portref("tpsuof")
  expect_equal(grprr(t$n, 250)$rates$rate, t$rate, tolerance = 1e-15)
})
