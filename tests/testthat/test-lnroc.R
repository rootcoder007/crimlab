# Eyewitness lineups. The kernel is checked by hand in src/t_lnr.cpp;
# here the binding against the rates counted in base R.

lnsim <- function(n = 300, seed = 3) {
  set.seed(seed)
  tp <- rep(0:1, each = n / 2)
  proc <- rep(c("sim", "seq"), n / 2)
  pick <- ifelse(
    tp == 1,
    sample(c("suspect", "filler", "none"), n, TRUE, c(0.5, 0.15, 0.35)),
    sample(c("suspect", "filler", "none"), n, TRUE, c(0.1, 0.3, 0.6))
  )
  conf <- ifelse(pick == "none", NA, sample(1:5, n, TRUE))
  data.frame(tp, proc, pick, conf)
}

test_that("ROC points, diagnosticity and the partial area follow the counts", {
  d <- lnsim()
  r <- lnroc(d$tp, d$pick, d$conf, proc = d$proc)
  expect_s3_class(r, "lnroc")
  lev <- sort(unique(d$conf[d$pick != "none"]), decreasing = TRUE)
  for (p in c("seq", "sim")) {
    s <- d[d$proc == p, ]
    hr <- vapply(
      lev,
      function(c) {
        sum(s$tp == 1 & s$pick == "suspect" & s$conf >= c, na.rm = TRUE)
      },
      1
    ) /
      sum(s$tp == 1)
    fa <- vapply(
      lev,
      function(c) {
        sum(s$tp == 0 & s$pick == "suspect" & s$conf >= c, na.rm = TRUE)
      },
      1
    ) /
      sum(s$tp == 0)
    k <- r$roc$proc == p
    expect_equal(r$roc$hr[k], hr, tolerance = 1e-14)
    expect_equal(r$roc$far[k], fa, tolerance = 1e-14)
    expect_equal(r$roc$dr[k], hr / fa, tolerance = 1e-14)
    # the trapezoids from (0, 0) up to the common false alarm rate
    x <- c(0, fa)
    y <- c(0, hr)
    f <- r$far
    keep <- x <= f
    xs <- c(x[keep], f)
    ys <- c(y[keep], approx(x, y, f, ties = "ordered", rule = 2)$y)
    if (max(x) < f) {
      ys[length(ys)] <- y[length(y)]
    }
    a <- sum(diff(xs) * (head(ys, -1) + tail(ys, -1)) / 2)
    expect_equal(r$pauc$pauc[r$pauc$proc == p], a, tolerance = 1e-12)
    hb <- diff(c(0, hr))
    fb <- diff(c(0, fa))
    expect_equal(
      r$cac$acc[r$cac$proc == p],
      ifelse(hb + fb > 0, hb / (hb + fb), NA),
      tolerance = 1e-14
    )
  }
  expect_identical(
    r$far,
    min(tapply(seq_len(nrow(r$roc)), r$roc$proc, function(i) max(r$roc$far[i])))
  )
  expect_true(all(is.na(r$pauc$se)))
  expect_output(print(r), "partial area")
})

test_that("the bootstrap compares procedures; size counts absent picks", {
  d <- lnsim(200, 4)
  r <- lnroc(
    d$tp,
    d$pick,
    d$conf,
    proc = d$proc,
    nboot = 200,
    seed = 2,
    far = 0.1
  )
  expect_identical(r$far, 0.1)
  expect_true(all(r$pauc$se > 0))
  expect_true(is.na(r$pauc$p[1]))
  expect_gte(r$pauc$p[2], 0)
  expect_lte(r$pauc$p[2], 1)
  expect_identical(
    lnroc(d$tp, d$pick, d$conf, proc = d$proc, nboot = 50, seed = 2),
    lnroc(d$tp, d$pick, d$conf, proc = d$proc, nboot = 50, seed = 2)
  )
  s <- lnroc(d$tp, d$pick, d$conf, size = 6)
  c1 <- max(d$conf, na.rm = TRUE)
  expect_equal(
    s$roc$far[1],
    sum(d$tp == 0 & d$pick != "none" & d$conf >= c1, na.rm = TRUE) /
      (6 * sum(d$tp == 0)),
    tolerance = 1e-14
  )
  expect_null(s$roc$proc)
})

test_that("lineup fairness by hand and against the binomial", {
  ch <- c(rep(3, 15), rep(1, 5), rep(2, 6), rep(4, 4), rep(5, 7), rep(6, 3))
  r <- lnsize(ch, k = 6, susp = 3)
  expect_s3_class(r, "lnsize")
  p <- as.numeric(table(factor(ch, 1:6))) / 40
  expect_identical(r$count, as.numeric(table(factor(ch, 1:6))))
  expect_equal(r$stat[["etr"]], 1 / sum(p^2), tolerance = 1e-14)
  expect_equal(
    r$stat[["esize"]],
    6 - sum(abs(p * 40 - 40 / 6)) / (2 * 40 / 6),
    tolerance = 1e-14
  )
  expect_equal(
    r$stat[["pbin"]],
    pbinom(14, 40, 1 / 6, lower.tail = FALSE),
    tolerance = 1e-12
  )
  expect_identical(r$stat[["fsize"]], 40 / 15)
  h <- qnorm(0.975) * sqrt(4 * (sum(p^3) - sum(p^2)^2) / 40)
  expect_equal(r$stat[["etrlo"]], 1 / (sum(p^2) + h), tolerance = 1e-12)
  expect_equal(r$stat[["etrhi"]], 1 / (sum(p^2) - h), tolerance = 1e-12)
  # no member left unchosen: both effective sizes agree
  expect_identical(r$stat[["esize"]], r$stat[["esizeb"]])
  expect_output(print(r), "40 mock witnesses")
  expect_identical(lnsize(c(1, 2, 2))$stat[["k"]], 2)
})

# The mock witnesses' choices of a table of counts, 0 for "not present".
mock <- function(cnt, none = 0) c(rep(seq_along(cnt), cnt), rep(0, none))

test_that("the worked examples of Tredoux are reproduced", {
  # Tredoux (1996), Table 4.2 (p. 70): lineup, E_a (adjusted for null
  # foils), E_b (not adjusted); the suspect is member 4
  t42 <- list(
    g = list(c(0, 25, 5, 25, 3, 2), 2.83, 3.00),
    h = list(c(10, 10, 9, 10, 11, 10), 5.90, 5.90),
    j = list(c(7, 7, 7, 24, 8, 7), 4.60, 4.60),
    k = list(c(12, 6, 9, 13, 14, 6), 5.10, 5.10),
    l = list(c(6, 19, 3, 20, 8, 10), 4.45, 4.45)
  )
  for (x in t42) {
    s <- lnsize(mock(x[[1]]), k = 6, susp = 4)$stat
    expect_equal(round(s[["esize"]], 2), x[[2]], tolerance = 1e-12)
    expect_equal(round(s[["esizeb"]], 2), x[[3]], tolerance = 1e-12)
  }
  # lineup i (1, 0, 12, 12, 0, 11): E_b 3.17 as printed; E_a is printed
  # 3.17 too, but the formula gives 4 - 16 / 18 = 3.11
  s <- lnsize(mock(c(1, 0, 12, 12, 0, 11)), k = 6, susp = 4)$stat
  expect_equal(round(s[["esizeb"]], 2), 3.17)
  expect_equal(s[["esize"]], 4 - 16 / 18, tolerance = 1e-14)
  # Table 4.1 (p. 67): functional size with "not present" choices
  t41 <- list(
    d = list(c(5, 3, 6, 30, 9, 3), 4, 2),
    e = list(c(8, 9, 8, 10, 9, 9), 7, 6),
    f = list(c(1, 3, 20, 10, 21, 4), 1, 6)
  )
  for (x in t41) {
    s <- lnsize(mock(x[[1]], x[[2]]), k = 6, susp = 4)$stat
    expect_equal(s[["fsize"]], x[[3]], tolerance = 1e-14)
    expect_identical(s[["n"]], 60)
    expect_identical(s[["nm"]], 60 - x[[2]])
  }
  # Table 6.6 (p. 144): I .7326, E' 3.740; 95% interval of I .7015 to
  # .7637, of E' to 4.232 (the lower end printed 3.503 is 1 / (1 - .7015)
  # = 3.350); the effective size of equation 5 is 3.65 (p. 144)
  s <- lnsize(mock(c(3, 24, 30, 33, 10)), k = 5, susp = 4)$stat
  expect_equal(round(1 - 1 / s[["etr"]], 4), 0.7326)
  expect_equal(round(s[["etr"]], 3), 3.740)
  expect_equal(round(1 - 1 / s[["etrlo"]], 4), 0.7015)
  expect_equal(round(1 - 1 / s[["etrhi"]], 4), 0.7637)
  # the upper end printed 4.232 is 1 / (1 - .7637), from the rounded end
  # of I; unrounded it is 4.231
  expect_equal(round(1 / (1 - round(1 - 1 / s[["etrhi"]], 4)), 3), 4.232)
  expect_equal(round(s[["etrhi"]], 3), 4.231)
  expect_equal(round(s[["etrlo"]], 3), 3.350)
  expect_equal(round(s[["esizeb"]], 2), 3.65)
  # Tredoux (1999), Table 1 (p. S10): 21 mock witnesses, the suspect
  # second; proportion .38, effective size 4.71 (functional size printed
  # 2.61; 21 / 8 is 2.625)
  s <- lnsize(mock(c(2, 8, 2, 3, 3, 3)), k = 6, susp = 2)$stat
  expect_equal(round(s[["psusp"]], 2), 0.38)
  expect_equal(round(s[["esize"]], 2), 4.71)
  expect_identical(s[["fsize"]], 21 / 8)
})

test_that("the partial areas of Mickes, Flowe and Wixted are reproduced", {
  # Table 3 (p. 372), Experiment 1a: hit and false alarm rates at
  # confidence 100 down to 10, made from 150 target-present and 175
  # target-absent lineups of six members (counts rounded to the printed
  # rates; the false alarm rate is a filler identification over 6 lineups)
  conf <- seq(100, 10, by = -10)
  sim <- list(
    hr = c(0.073, 0.147, 0.260, 0.427, 0.487, 0.553, 0.620, 0.640, 0.640,
           0.640),
    far = c(0.000, 0.006, 0.019, 0.029, 0.037, 0.049, 0.053, 0.057, 0.057,
            0.057)
  )
  sq <- list(
    hr = c(0.047, 0.107, 0.213, 0.247, 0.247, 0.267, 0.293, 0.320, 0.327,
           0.333),
    far = c(0.002, 0.003, 0.013, 0.017, 0.022, 0.026, 0.029, 0.037, 0.046,
            0.049)
  )
  lineups <- function(r, proc) {
    nh <- diff(c(0, round(r$hr * 150)))
    nf <- diff(c(0, round(r$far * 6 * 175)))
    tp <- c(rep(1, 150), rep(0, 175))
    pick <- c(
      rep("suspect", sum(nh)), rep("none", 150 - sum(nh)),
      rep("filler", sum(nf)), rep("none", 175 - sum(nf))
    )
    cf <- c(
      rep(conf, nh), rep(NA, 150 - sum(nh)),
      rep(conf, nf), rep(NA, 175 - sum(nf))
    )
    data.frame(tp, pick, cf, proc)
  }
  d <- rbind(lineups(sim, "sim"), lineups(sq, "seq"))
  # the range: 0 to a false alarm rate slightly above the simultaneous
  # procedure's largest, 0.0575 here (0.345 on their scale)
  r <- lnroc(d$tp, d$pick, d$cf, proc = d$proc, size = 6, far = 0.0575)
  k <- r$roc$proc == "sim"
  expect_equal(round(r$roc$hr[k], 3), sim$hr, tolerance = 1e-12)
  expect_equal(round(r$roc$far[k], 3), sim$far, tolerance = 1e-12)
  expect_equal(round(r$roc$far[!k], 3), sq$far, tolerance = 1e-12)
  # p. 371: partial areas 0.13 (simultaneous) and 0.09 (sequential)
  pa <- setNames(r$pauc$pauc, r$pauc$proc)
  expect_equal(round(6 * pa[["sim"]], 2), 0.13)
  expect_equal(round(6 * pa[["seq"]], 2), 0.09)
  # past its last point (0.049) the sequential curve runs to (1 / 6, 1)
  x0 <- tail(r$roc$far[!k], 1)
  y0 <- tail(r$roc$hr[!k], 1)
  y1 <- y0 + (1 - y0) * (0.0575 - x0) / (1 / 6 - x0)
  x <- c(0, r$roc$far[!k], 0.0575)
  y <- c(0, r$roc$hr[!k], y1)
  expect_equal(
    pa[["seq"]],
    sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2),
    tolerance = 1e-12
  )
})

test_that("lnroc() and lnsize() refuse bad input in words", {
  errs(
    "lnroc",
    list(
      list("`pick` must be", c(1, 0), c("suspect", "guess"), c(1, 1)),
      list(
        "`conf` must be a finite number",
        c(1, 0),
        c("suspect", "none"),
        c(NA, 1)
      ),
      list(
        "target-present and target-absent",
        c(1, 1),
        c("suspect", "none"),
        c(1, NA)
      ),
      list("no witness made", c(1, 0), c("none", "none"), c(NA, NA)),
      list("`size`", c(1, 0), c("suspect", "none"), c(1, NA), size = 1)
    )
  )
  errs(
    "lnsize",
    list(
      list("`choice` must be whole", c(1, 2.5)),
      list("`choice` must be at least 0", c(1, -1)),
      list("`choice` must name a lineup member", c(0, 0)),
      list("position 7 of a lineup of 6", c(1, 7), k = 6),
      list("`susp` must be a position", c(1, 2), susp = 3)
    )
  )
})

test_that("the bootstrap errors do not depend on how many replicates exist", {
  # the replicates are drawn twice from the seed (sums, then squared
  # deviations) instead of being kept: a run is reproducible and its
  # errors finite
  d <- lnsim()
  a <- lnroc(d$tp, d$pick, d$conf, proc = d$proc, nboot = 40, seed = 2)
  b <- lnroc(d$tp, d$pick, d$conf, proc = d$proc, nboot = 40, seed = 2)
  expect_identical(a, b)
  expect_true(all(is.finite(a$pauc$se)))
})

test_that("lnsize() refuses a lineup larger than its bound", {
  expect_error(lnsize(c(1, 2, 2), k = 1e9),
               "lnsize: `k` must be at most 67108864.", fixed = TRUE)
})
