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
  se <- sqrt(4 * (sum(p^3) - sum(p^2)^2) / 40) / sum(p^2)^2
  expect_equal(
    r$stat[["etrhi"]] - r$stat[["etr"]],
    qnorm(0.975) * se,
    tolerance = 1e-12
  )
  expect_output(print(r), "40 mock witnesses")
  expect_identical(lnsize(c(1, 2, 2))$stat[["k"]], 2)
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
      list("position 7 of a lineup of 6", c(1, 7), k = 6),
      list("`susp` must be a position", c(1, 2), susp = 3)
    )
  )
})
