# trndrr(), stpcp(), cusrv(), xcser(): the rate ratio against glm(), the
# Mann-Kendall p-value against cor.test(), Theil-Sen and Sen's interval
# written out here, the seasons against chisq.test(), the step scan
# against every ordering enumerated here, CUSUM by its recursion, and the
# cross-correlations against ccf() and cor.test().

perms <- function(v) {
  if (length(v) <= 1L) return(list(v))
  out <- list()
  for (i in seq_along(v)) {
    for (p in perms(v[-i])) out[[length(out) + 1L]] <- c(v[i], p)
  }
  out
}

test_that("the rate ratio is glm()'s, Poisson and quasi-Poisson", {
  y <- c(210, 198, 190, 181, 160, 158, 140)
  x <- 2017:2023
  e <- c(10, 10.2, 10.3, 10.5, 10.4, 10.6, 10.9)
  r <- trndrr(y, x, expo = e)
  g <- summary(glm(y ~ I(x - mean(x)) + offset(log(e)), family = poisson))
  co <- coef(g)[2, ]
  expect_equal(r$rr$b, co[[1]], tolerance = 1e-10)
  expect_equal(r$rr$se, co[[2]], tolerance = 1e-10)
  expect_equal(r$rr$p, co[[4]], tolerance = 1e-8)
  ci <- exp(co[[1]] + c(-1, 1) * qnorm(0.975) * co[[2]])
  expect_equal(c(r$rr$lo, r$rr$hi), ci, tolerance = 1e-10)
  expect_false(r$rr$over)
  w <- c(10, 50, 5, 60, 8, 70)
  q <- trndrr(w)
  # the Pearson dispersion at the fitted counts; glm(family =
  # quasipoisson) takes it from the last iteration's working weights,
  # which differ in the sixth digit
  gp <- glm(w ~ I(seq_along(w) - 3.5), family = poisson)
  d <- sum(residuals(gp, type = "pearson")^2) / 4
  co <- coef(summary(gp))[2, ]
  se <- co[[2]] * sqrt(d)
  expect_true(q$rr$over)
  expect_equal(q$rr$disp, d, tolerance = 1e-10)
  expect_equal(q$rr$se, se, tolerance = 1e-10)
  expect_equal(q$rr$p, 2 * pt(-abs(co[[1]] / se), 4), tolerance = 1e-8)
  expect_equal(q$rr$hi, exp(co[[1]] + qt(0.975, 4) * se), tolerance = 1e-10)
  gq <- glm(w ~ seq_along(w), family = quasipoisson)
  expect_equal(q$rr$disp, summary(gq)$dispersion, tolerance = 1e-5)
  expect_equal(q$fitted, unname(fitted(gp)), tolerance = 1e-10)
})

test_that("Mann-Kendall is cor.test()'s Kendall test without ties", {
  y <- c(402, 377, 390, 368, 331, 340, 305)
  r <- trndrr(y)
  k <- cor.test(seq_along(y), y, method = "kendall", exact = TRUE)
  expect_equal(r$mk$tau, unname(k$estimate), tolerance = 1e-14)
  expect_equal(r$mk$p, k$p.value, tolerance = 1e-12)
  expect_true(r$mk$exact)
  ku <- cor.test(seq_along(y), y, method = "kendall", exact = TRUE,
                 alternative = "less")
  expect_equal(trndrr(y, alt = "down")$mk$p, ku$p.value, tolerance = 1e-12)
})

test_that("past 8 periods: the normal approximation with ties, by hand", {
  y <- c(5, 7, 7, 6, 9, 8, 8, 8, 11, 10, 12, 12)
  r <- trndrr(y)
  n <- length(y)
  s <- sum(sapply(1:(n - 1), function(i) sum(sign(y[(i + 1):n] - y[i]))))
  tt <- table(y)
  v <- (n * (n - 1) * (2 * n + 5) - sum(tt * (tt - 1) * (2 * tt + 5))) / 18
  expect_equal(r$mk$s, s)
  expect_equal(r$mk$var, v, tolerance = 1e-14)
  expect_equal(r$mk$p, 2 * pnorm(-(s - 1) / sqrt(v)), tolerance = 1e-12)
  expect_false(r$mk$exact)
  # Theil-Sen and Sen's interval
  sl <- outer(y, y, "-") / outer(seq_len(n), seq_len(n), "-")
  sl <- sort(sl[lower.tri(sl)])
  m <- length(sl)
  ca <- qnorm(0.975) * sqrt(v)
  expect_equal(r$sen$slope, median(sl))
  expect_equal(r$sen$int, median(y - median(sl) * seq_len(n)))
  expect_equal(r$sen$lo, sl[floor((m - ca) / 2)])
  expect_equal(r$sen$hi, sl[ceiling((m + ca) / 2) + 1])
  expect_false(r$sen$clamp)
})

test_that("year-over-year changes and seasons", {
  y <- c(5, 9, 7, 4, 6, 10, 8, 5)
  r <- trndrr(y, seas = rep(c("q1", "q2", "q3", "q4"), 2))
  expect_equal(r$yoy$chg, c(NA, diff(y)))
  expect_equal(r$yoy$pct, c(NA, 100 * diff(y) / y[-8]))
  tot <- tapply(y, rep(1:4, 2), sum)
  ch <- chisq.test(tot)
  expect_equal(r$seas$stat, unname(ch$statistic), tolerance = 1e-14)
  expect_equal(r$seas$p, ch$p.value, tolerance = 1e-14)
  expect_identical(r$seas$df, 3L)
  expect_identical(r$seas$tab$season, c("q1", "q2", "q3", "q4"))
  expect_true(is.na(trndrr(c(0, 1, 2))$yoy$pct[2]))
  expect_output(print(r), "seasons: chi-square")
  expect_output(print(trndrr(c(10, 50, 5, 60, 8, 70))), "overdispersed")
})

test_that("events in one end period: no finite trend", {
  expect_warning(r <- trndrr(c(0, 0, 5)), "last period")
  expect_identical(r$rr$rr, Inf)
  expect_true(is.na(r$rr$lo))
  expect_warning(r0 <- trndrr(c(4, 0, 0)), "first period")
  expect_identical(r0$rr$rr, 0)
  # in the middle the estimate is finite
  expect_true(is.finite(trndrr(c(0, 3, 0))$rr$rr))
})

test_that("the step scan's p-value counts every ordering", {
  y <- c(100, 104, 98, 60, 63, 58)
  r <- stpcp(y, x = 2018:2023)
  tk <- function(v, m = 2) {
    n <- length(v)
    max(sapply(m:(n - m), function(k) {
      k * (n - k) / n * (mean(v[1:k]) - mean(v[-(1:k)]))^2
    }))
  }
  all <- vapply(perms(y), tk, 0)
  expect_equal(r$tab$stat[1], tk(y))
  expect_equal(r$tab$p[1], mean(all >= tk(y) - 1e-9), tolerance = 1e-14)
  expect_identical(r$tab$at[1], 2020L)
  expect_true(r$exact)
  expect_identical(r$nperm, 720)
  # Pettitt by its definition
  n <- 6
  u <- sapply(1:(n - 1), function(t) {
    sum(sign(outer(y[(t + 1):n], y[1:t], "-")))
  })
  kk <- max(abs(u))
  expect_equal(r$tab$stat[2], kk)
  expect_equal(r$tab$p[2], min(1, 2 * exp(-6 * kk^2 / (n^3 + n^2))))
  expect_identical(r$tab$k[2], which.max(abs(u)))
  expect_output(print(r), "exact over 720")
  m <- stpcp(c(y, y + 1), nperm = 99, seed = 3)
  expect_identical(m, stpcp(c(y, y + 1), nperm = 99, seed = 3))
  expect_false(m$exact)
  expect_output(print(m), "from 99")
})

test_that("CUSUM follows its recursion", {
  y <- c(10, 12, 9, 11, 10, 13, 19, 22, 25, 12)
  r <- cusrv(y, start = 6)
  m <- mean(y[1:5])
  z <- (y[6:10] - m) / sqrt(m)
  s <- Reduce(function(a, b) max(0, a + b - 1.04), z, accumulate = TRUE,
              init = 0)[-1]
  expect_equal(r$z, z)
  expect_equal(r$s, s)
  expect_identical(r$alarm, s >= 2.26)
  expect_identical(r$i, 6:10)
  e <- cusrv(y, expd = rep(11, 10), trans = "none", k = 12, h = 5)
  expect_equal(e$s, Reduce(function(a, b) max(0, a + b - 12), y,
                           accumulate = TRUE, init = 0)[-1])
})

test_that("the correlations are cor.test()'s and ccf()'s", {
  x <- c(5, 3, 8, 6, 9, 12, 10, 14, 13, 15)
  y <- c(2, 4, 3, 7, 6, 8, 11, 9, 13, 12)
  r <- xcser(x, y, lagmax = 3)
  k <- cor.test(x, y)
  expect_equal(r$cor$r, unname(k$estimate), tolerance = 1e-14)
  expect_equal(c(r$cor$lo, r$cor$hi), as.vector(k$conf.int), tolerance = 1e-14)
  expect_equal(r$cor$p, k$p.value, tolerance = 1e-12)
  cc <- ccf(x, y, lag.max = 3, plot = FALSE)
  expect_equal(r$ccf$r, as.vector(cc$acf), tolerance = 1e-14)
  expect_identical(r$ccf$lag, -3:3)
  d <- ccf(x, y, plot = FALSE)
  expect_identical(nrow(xcser(x, y)$ccf), length(d$lag))
  dd <- xcser(x, y, dif = TRUE)
  expect_equal(dd$cor$r, cor(diff(x), diff(y)), tolerance = 1e-14)
  # by name: the years both have
  a <- c("2019" = 120, "2020" = 98, "2021" = 105, "2022" = 130, "2023" = 141)
  b <- c("2023" = 70, "2018" = 50, "2019" = 61, "2020" = 44, "2021" = 52,
         "2022" = 66)
  n <- xcser(a, b)
  expect_identical(n$common, names(a))
  expect_equal(n$cor$r, cor(a, b[names(a)]), tolerance = 1e-14)
  expect_output(print(n), "5 periods")
  three <- xcser(1:3, c(1, 3, 2))
  expect_true(is.na(three$cor$lo))
  expect_lt(xcser(1:4, 2 * (1:4))$cor$p, 1e-10)
})

test_that("bad input is refused", {
  errs("trndrr", list(
    list("at least 3 periods", c(1, 2)),
    list("no events", c(0, 0, 0)),
    list("strictly increasing", c(1, 2, 3), c(1, 3, 2)),
    list("`expo`", c(1, 2, 3), expo = c(1, 0, 1)),
    list("one season", c(1, 2, 3), seas = c("a", "a", "a")),
    list("`alt` must be one of", c(1, 2, 3), alt = "two.sided"),
    list("at least 0", c(1, -2, 3))
  ))
  errs("stpcp", list(
    list("at least 2 `minseg`", c(1, 2, 3)),
    list("one period per value", 1:6, x = 1:5),
    list("`nperm`", 1:6, nperm = 0)
  ))
  errs("cusrv", list(
    list("give `expd`", c(1, 2)),
    list("at most the number", c(1, 2), start = 3),
    list("all 0", c(0, 0, 3), start = 3),
    list("`expd`", c(1, 2), expd = c(1, 0)),
    list("`trans` must be one of", c(1, 2), start = 2, trans = "log")
  ))
  errs("xcser", list(
    list("one length", 1:4, 1:5),
    list("at least 3", c(a = 1, b = 2), c(a = 1, b = 2)),
    list("below the number", 1:4, c(2, 1, 4, 3), lagmax = 4),
    list("constant", 1:4, c(1, 1, 1, 1))
  ))
})
