# stkflw(), lostay(), bklog(), dktflw(), lftab() and trnmat(): against the
# identities written out, quantile(), t.test(), dlnorm() and dexp(),
# cut(), chisq.test(), prop.trend.test(), eigen() and table().

test_that("stkflw() fills in the stock and flow from any two measures", {
  p <- c(5200, 4800, 4300)
  d <- c(310000, 318000, 330000)
  r <- stkflw(
    people = p,
    days = d,
    prd = 2021:2023,
    expo = 1e6,
    rel = c(5000, 4700, 4400),
    early = 182.5
  )
  f <- r$flow
  expect_identical(f$period, c("2021", "2022", "2023"))
  expect_equal(f$adp, d / 365, tolerance = 1e-15)
  expect_equal(f$alos, d / p, tolerance = 1e-15)
  expect_equal(f$peochg, 100 * (p / p[1] - 1), tolerance = 1e-14)
  expect_equal(
    f$aloschg,
    100 * ((d / p) / (d[1] / p[1]) - 1),
    tolerance = 1e-13
  )
  # days = people x alos, so the log changes add up
  expect_equal(
    f$pshare[-1],
    log(p[-1] / p[1]) / log(d[-1] / d[1]),
    tolerance = 1e-12
  )
  expect_equal(f$stkrate, 1e5 * d / 365 / 1e6, tolerance = 1e-14)
  expect_equal(f$cut, c(5000, 4700, 4400) * 182.5 / 365, tolerance = 1e-15)
  expect_equal(r$mean[["alos"]], mean(d / p), tolerance = 1e-14)
  expect_identical(r$per, 1e5)
  expect_output(print(r), "against the first period")
  # every pair gives the same table
  a <- stkflw(adp = d / 365, alos = d / p)$flow
  expect_equal(a$people, p, tolerance = 1e-12)
  b <- stkflw(people = p, alos = d / p)$flow
  expect_equal(b$days, d, tolerance = 1e-12)
  q <- stkflw(days = d, alos = d / p, base = "previous")$flow
  expect_equal(q$people, p, tolerance = 1e-12)
  expect_identical(q$peochg[1], NA_real_)
  expect_true(is.na(r$flow$pshare[1]))
  # headcounts give the average daily population
  h <- stkflw(people = 10, adp = list(c(8, 9, 10, 9)))$flow
  expect_identical(h$adp, 9)
  expect_identical(h$days, 9 * 365)
  # periods from dates
  dt <- stkflw(
    people = 10,
    days = 300,
    from = "2024-01-01",
    to = as.Date("2024-01-31")
  )
  expect_identical(dt$flow$t, 31)
  expect_named(
    dt$flow,
    c(
      "period",
      "t",
      "people",
      "days",
      "adp",
      "alos",
      "peochg",
      "aloschg",
      "dayschg",
      "adpchg",
      "pshare"
    )
  )
  expect_true(is.na(dt$per))
})

test_that("stkflw() refuses bad input in words", {
  errs(
    "stkflw",
    list(
      list("give two of", people = 1),
      list("not `days` with `adp`", days = 1, adp = 1),
      list("`people` must be greater than 0", people = 0, days = 1),
      list("`adp` must have one value per period, 2", people = 1:2, adp = 1),
      list(
        "`t` must be one number or one per period",
        people = 1:2,
        days = 1:2,
        t = 1:3
      ),
      list(
        "`expo` must be one number or one per period",
        people = 1:2,
        days = 1:2,
        expo = 1:3
      ),
      list("`base` must be one of", people = 1, days = 1, base = "last"),
      list("`rel` and `early` go together", people = 1, days = 1, rel = 1),
      list(
        "`early` must be one number from 0",
        people = 1,
        days = 1,
        rel = 1,
        early = -1
      ),
      list(
        "`prd` must have one label per period",
        people = 1:2,
        days = 1:2,
        prd = 1
      ),
      list(
        "`prd` must be labels",
        people = 1,
        days = 1,
        prd = as.Date("2024-01-01")
      ),
      list(
        "`from` and `to` must both be dates",
        people = 1,
        days = 1,
        from = 1,
        to = 2
      ),
      list(
        "`from` holds something that is not a date",
        people = 1,
        days = 1,
        from = "soon",
        to = "2024-01-01"
      ),
      list(
        "`to` must not precede `from`",
        people = 1,
        days = 1,
        from = "2024-02-01",
        to = "2024-01-01"
      )
    )
  )
})

test_that("lostay() summarises, bands and fits the stays", {
  set.seed(3)
  s <- round(rlnorm(300, 2.3, 1.1)) + 1
  fl <- s > 40 & runif(300) < 0.5
  r <- lostay(s, breaks = c(0, 5, 15, 31, 61, Inf), thr = c(15, 60), flag = fl)
  st <- r$stay
  expect_equal(st[["mean"]], mean(s), tolerance = 1e-14)
  expect_equal(st[["sd"]], sd(s), tolerance = 1e-13)
  expect_equal(
    unname(st[c("q1", "median", "q3")]),
    unname(quantile(s, c(0.25, 0.5, 0.75))),
    tolerance = 1e-15
  )
  expect_equal(
    unname(st[c("lo", "hi")]),
    as.numeric(t.test(s)$conf.int),
    tolerance = 1e-12
  )
  tb <- table(cut(s, c(0, 5, 15, 31, 61, Inf)))
  expect_equal(r$bands$n, as.numeric(tb), tolerance = 0)
  expect_equal(r$bands$share, as.numeric(tb) / 300, tolerance = 1e-15)
  fb <- tapply(fl, cut(s, c(0, 5, 15, 31, 61, Inf)), mean)
  expect_equal(r$bands$flagged, as.numeric(fb), tolerance = 1e-15)
  expect_equal(r$thr$atmost, c(mean(s <= 15), mean(s <= 60)), tolerance = 1e-15)
  # the fits by maximum likelihood
  ml <- mean(log(s))
  sl <- sqrt(mean((log(s) - ml)^2))
  aic <- c(
    4 - 2 * sum(dlnorm(s, ml, sl, log = TRUE)),
    2 - 2 * sum(dexp(s, 1 / mean(s), log = TRUE))
  )
  expect_equal(r$fit$aic[1:2], aic, tolerance = 1e-12)
  al <- 300 / sum(log(s / min(s)))
  pll <- sum(log(al) + al * log(min(s)) - (al + 1) * log(s))
  expect_equal(r$fit$aic[3], 4 - 2 * pll, tolerance = 1e-12)
  expect_identical(r$fit$best, r$fit$aic == min(r$fit$aic))
  expect_equal(unname(r$par), c(ml, sl, al), tolerance = 1e-13)
  expect_output(print(r), "stays, mean")
  # frequencies count as repeated stays
  w <- lostay(c(3, 9, 20), w = c(4, 2, 1))
  e <- lostay(rep(c(3, 9, 20), c(4, 2, 1)))
  expect_equal(w$stay, e$stay, tolerance = 1e-14)
  expect_null(w$bands)
  one <- lostay(5)
  expect_true(is.na(one$stay[["lo"]]))
  expect_true(all(is.na(one$fit$aic)))
  expect_identical(one$fit$best, c(FALSE, FALSE, FALSE))
})

test_that("lostay() refuses bad input in words", {
  errs(
    "lostay",
    list(
      list("`days` must be at least 0", -1),
      list("`w` must be whole-number counts", 1, w = 0.5),
      list("`w` must have a positive total", 1, w = 0),
      list("`breaks` must be at least two band edges", 1, breaks = 1),
      list("`breaks` must be increasing", 1, breaks = c(2, 1)),
      list("`thr` must be at least 0", 1, thr = -1),
      list("`flag` must be TRUE or FALSE", 1:2, flag = NA),
      list("`level` must be one number", 1, level = 0)
    )
  )
})

test_that("bklog() holds Little's law", {
  a <- c(0, 2, 3, 7, 9, 11)
  d <- c(5, 4, 12, NA, 15, 30)
  r <- bklog(a, d, horizon = 20, target = 2)
  e <- r$est
  occ <- sum(pmin(ifelse(is.na(d), 20, d), 20) - a)
  expect_equal(e[["occup"]], occ, tolerance = 1e-15)
  expect_equal(e[["backlog"]], occ / 20, tolerance = 1e-15)
  expect_identical(e[["pending"]], 2)
  expect_equal(e[["wait"]], mean((d - a)[1:5][-4]), tolerance = 1e-15)
  expect_true(is.na(e[["resid"]]))
  expect_equal(r$need, 2 / (6 / 20), tolerance = 1e-15)
  # nothing pending: L = lambda W exactly
  f <- bklog(c(0, 1, 2), c(3, 2, 6))
  expect_identical(f$horizon, 6)
  expect_equal(f$est[["resid"]], 0, tolerance = 1e-15)
  expect_null(f$need)
  expect_output(print(f), "Little's law")
  errs(
    "bklog",
    list(
      list("`arrive` must be at least 0", -1, 1),
      list("`dispose` must be numbers", 1, "a"),
      list("`dispose` must be finite or NA", 1, Inf),
      list("precedes its arrival \\(case 1\\)", 5, 1),
      list("inside the horizon", 5, 6, horizon = 2),
      list("`target` must be one number from 0", 1, 2, target = -1)
    )
  )
})

test_that("dktflw() gives clearance, rolling clearance and Little's law", {
  rc <- c(120, 130, 125, 140, 90)
  ds <- c(110, 135, 120, 150, 95)
  pd <- c(400, 395, 400, 390, 385)
  r <- dktflw(
    rc,
    ds,
    pd,
    prd = 2019:2023,
    win = 2,
    span = 91,
    mtime = c(300, 290, 310, 280, 300),
    shares = rbind(c(0.5, 0.3, 0.2))[rep(1, 5), ],
    floors = c(0, 90, 180)
  )
  f <- r$flow
  expect_equal(f$clear, ds / rc, tolerance = 1e-15)
  lo <- qbeta(0.025, ds, rc + 1)
  hi <- qbeta(0.975, ds + 1, rc)
  expect_equal(f$lo, lo / (1 - lo), tolerance = 1e-14)
  expect_equal(f$hi, hi / (1 - hi), tolerance = 1e-14)
  expect_equal(
    f$rclear[2:5],
    (ds[1:4] + ds[2:5]) / (rc[1:4] + rc[2:5]),
    tolerance = 1e-15
  )
  expect_identical(f$pchg, c(NA, diff(pd)))
  pace <- c(ds[1], (ds[1:4] + ds[2:5]) / 2)
  expect_equal(f$toclear, pd / pace, tolerance = 1e-15)
  expect_equal(f$frate, rc / 91, tolerance = 1e-15)
  expect_equal(f$impl, pd / (rc / 91), tolerance = 1e-15)
  expect_equal(f$ratio, f$impl / c(300, 290, 310, 280, 300), tolerance = 1e-15)
  expect_equal(f$floorbd, rep(0.3 * 90 + 0.2 * 180, 5), tolerance = 1e-14)
  expect_null(r$mix)
  # without the pending caseload: the cumulative net
  n <- dktflw(rc, ds)$flow
  expect_identical(n$cumnet, cumsum(rc - ds))
  expect_true(all(is.na(n$rclear)))
  z <- dktflw(c(0, 5, 0), c(3, 0, 0))$flow
  expect_true(is.infinite(z$hi[1]))
  expect_identical(z$lo[2], 0)
  expect_true(is.na(z$lo[3]))
  expect_output(print(r), "case flow")
})

test_that("dktflw() tests the mix of outcomes over periods", {
  m <- matrix(
    c(50, 55, 60, 70, 30, 28, 25, 20, 20, 17, 15, 10),
    4,
    dimnames = list(2020:2023, c("guilty", "withdrawn", "other"))
  )
  r <- dktflw(mix = m, level = 0.9)$mix
  s <- r$shares
  expect_identical(s$period, rep(as.character(2020:2023), 3))
  expect_equal(s$share, as.vector(m / rowSums(m)), tolerance = 1e-15)
  a <- qchisq(1 - 0.1 / 3, 1)
  ni <- as.vector(m)
  nn <- rep(rowSums(m), 3)
  hf <- sqrt(a * (a + 4 * ni * (nn - ni) / nn))
  expect_equal(
    s$lo,
    (a + 2 * ni - hf) / (2 * (nn + a)),
    tolerance = 1e-14,
    ignore_attr = TRUE
  )
  ct <- chisq.test(m, correct = FALSE)
  expect_equal(r$homog[["chisq"]], unname(ct$statistic), tolerance = 1e-12)
  expect_equal(r$homog[["p"]], ct$p.value, tolerance = 1e-12)
  expect_equal(
    r$homog[["cramerv"]],
    sqrt(
      unname(ct$statistic) /
        (sum(m) * 2)
    ),
    tolerance = 1e-13
  )
  tr <- vapply(
    1:3,
    function(j) {
      unname(suppressWarnings(prop.trend.test(m[, j], rowSums(m)))$statistic)
    },
    0
  )
  expect_equal(r$trend$stat, tr, tolerance = 1e-12)
  pv <- pchisq(tr, 1, lower.tail = FALSE)
  expect_equal(r$trend$holm, p.adjust(pv, "holm"), tolerance = 1e-13)
  expect_equal(
    r$trend$change,
    unname(m[4, ] / sum(m[4, ]) - m[1, ] / sum(m[1, ])),
    tolerance = 1e-15
  )
  # unnamed, with scores, and an outcome never seen
  u <- dktflw(mix = cbind(c(5, 6, 9), c(5, 4, 1), 0), scores = c(1, 2, 4))$mix
  expect_identical(u$trend$outcome, c("outcome1", "outcome2", "outcome3"))
  expect_true(is.na(u$trend$stat[3]))
  expect_true(is.na(u$trend$holm[3]))
  expect_equal(
    u$trend$stat[1],
    unname(
      prop.trend.test(
        c(5, 6, 9),
        c(10, 10, 10),
        c(1, 2, 4)
      )$statistic
    ),
    tolerance = 1e-12
  )
  expect_output(print(dktflw(mix = m)), "withdrawn")
  one <- dktflw(mix = matrix(c(3, 4), 1))$mix
  expect_true(is.na(one$homog[["p"]]))
})

test_that("dktflw() refuses bad input in words", {
  errs(
    "dktflw",
    list(
      list("give `recv` and `disp`, or `mix`"),
      list("`disp` must come with `recv`", recv = 1),
      list(
        "`win` must be one whole number from 1",
        recv = 1,
        disp = 1,
        win = 0
      ),
      list("`span` needs `pend`", recv = 1, disp = 1, span = 30),
      list(
        "`span` must be one number or one per period",
        recv = 1:2,
        disp = 1:2,
        pend = 1:2,
        span = 1:3
      ),
      list("`mtime` needs `span` and `pend`", recv = 1, disp = 1, mtime = 3),
      list("`shares` must be a numeric matrix", recv = 1, disp = 1, shares = 1),
      list(
        "`floors` must be increasing",
        recv = 1,
        disp = 1,
        shares = matrix(c(0.5, 0.5), 1),
        floors = c(5, 1)
      ),
      list(
        "`shares` must be numbers from 0",
        recv = 1,
        disp = 1,
        shares = matrix(c(-0.5, 1.5), 1),
        floors = c(0, 1)
      ),
      list(
        "must sum to 1",
        recv = 1,
        disp = 1,
        shares = matrix(c(0.5, 0.2), 1),
        floors = c(0, 1)
      ),
      list("`mix` must be a numeric matrix", mix = 1:3),
      list("`mix` must have a period and an outcome", mix = matrix(0, 0, 2)),
      list("every period with a case", mix = matrix(c(1, 0, 1, 0), 2)),
      list("`scores` must have length 2", mix = matrix(1, 2, 2), scores = 1),
      list("`level` must be one number", recv = 1, disp = 1, level = 1)
    )
  )
})

test_that("lftab() builds the actuarial life table", {
  br <- c(0, 3, 6, 12, 24, Inf)
  d <- c(40, 30, 20, 10, 5)
  w <- c(2, 2, 3, 3, 4)
  r <- lftab(br, d, w, ci = "log")
  tb <- r$tab
  ent <- sum(d + w) - c(0, cumsum(d + w))[1:5]
  ar <- ent - w / 2
  q <- d / ar
  s <- cumprod(1 - q)
  gw <- cumsum(d / (ar * (ar - d)))
  expect_equal(tb$entering, ent, tolerance = 1e-15)
  expect_equal(tb$surv, s, tolerance = 1e-14)
  expect_equal(tb$se, s * sqrt(gw), tolerance = 1e-14)
  z <- qnorm(0.975)
  expect_equal(tb$lo, s * exp(-z * sqrt(gw)), tolerance = 1e-14)
  j <- which(s <= 0.5)[1]
  s0 <- c(1, s)[j]
  expect_equal(
    r$median,
    br[j] + (br[j + 1] - br[j]) * (s0 - 0.5) / (s0 - s[j]),
    tolerance = 1e-14
  )
  expect_true(is.na(tb$hazard[5]))
  expect_identical(r$ci, "log")
  expect_output(print(r), "life table")
  # band labels, and one censored count for every band
  l <- lftab(c("0-3", "3-6", "Over 6"), c(10, 5, 5), cens = 1, ci = "plain")
  expect_identical(l$tab$upper, c(3, 6, Inf))
  expect_identical(l$tab$cens, c(1, 1, 1))
  u <- lftab(c("Under 3", "3-6"), c(1, 1))
  expect_identical(u$tab$lower, c(0, 3))
  expect_identical(u$median, 3)
  expect_identical(lftab(c(0, 1, 2), c(1, 1), cens = 10)$median, Inf)
  errs(
    "lftab",
    list(
      list("these bands cannot be read: soon", c("0-3", "soon"), 1:2),
      list("`breaks` must be at least two band edges", 1, 1),
      list("`events` must have length 2", c(0, 1, 2), 1),
      list(
        "`cens` must be one number or one per band",
        c(0, 1, 2),
        1:2,
        cens = 1:3
      ),
      list("`ci` must be one of", c(0, 1), 1, ci = "exact"),
      list("must be contiguous", c(0, 2, 1), c(1, 1))
    )
  )
})

test_that("trnmat() counts moves and finds the stationary distribution", {
  set.seed(5)
  id <- rep(1:40, each = 6)
  st <- sample(c("N", "E", "W", "C"), 240, TRUE, prob = 4:1)
  tm <- rep(6:1, 40)
  r <- trnmat(st, id, tm)
  o <- order(id, tm)
  s2 <- st[o]
  i2 <- id[o]
  k <- which(i2[-1] == i2[-240])
  tb <- table(
    factor(s2[k], sort(unique(st))),
    factor(s2[k + 1], sort(unique(st)))
  )
  expect_equal(r$count, unclass(tb), ignore_attr = TRUE, tolerance = 0)
  expect_equal(
    r$prob,
    unclass(tb / rowSums(tb)),
    ignore_attr = TRUE,
    tolerance = 1e-15
  )
  e <- eigen(t(r$prob))
  v <- Re(e$vectors[, which.min(abs(e$values - 1))])
  expect_equal(unname(r$pi), v / sum(v), tolerance = 1e-12)
  expect_named(r$pi, sort(unique(st)))
  expect_equal(r$summ[["stay"]], sum(diag(tb)) / sum(tb), tolerance = 1e-15)
  off <- tb[row(tb) != col(tb)]
  p <- off[off > 0] / sum(off)
  expect_equal(
    r$summ[["theil"]],
    sum(p * log(p * length(p))),
    tolerance = 1e-13
  )
  expect_output(print(r), "moves")
  # row order without times; no moves at all
  n <- trnmat(c("a", "b", "a"), c(1, 1, 1))
  expect_identical(n$count[1, 2], 1)
  z <- trnmat(c("a", "b"), c(1, 2))
  expect_identical(z$summ[["moves"]], 0)
  expect_true(all(is.na(z$pi)))
  errs(
    "trnmat",
    list(
      list("`id` must have one label per value", c("a", "b"), 1),
      list("`time` must have length 2", c("a", "b"), c(1, 1), time = 1),
      list("more than 5000 states", 1:5001, rep(1, 5001))
    )
  )
})
