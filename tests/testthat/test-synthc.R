# synthc() and sdidx(): the weights against the optimality (KKT)
# conditions of their convex problems, written out here; exact recovery of
# a convex combination; the placebos and the jackknife recomputed from the
# weights; recovery of a known effect with the errors against the Monte
# Carlo spread; and the refusals.

scpanel <- function(seed, nu = 10, ntr = 1, eff = -3) {
  set.seed(seed)
  d <- expand.grid(year = 2010:2024, town = sprintf("u%02d", seq_len(nu)))
  i <- as.integer(d$town)
  d$first <- ifelse(i <= ntr, 2019, NA)
  d$rate <- 20 + rnorm(nu)[i] + (d$year - 2010) * runif(nu, 0.2, 0.8)[i] +
    rnorm(15 * nu, 0, 0.3) + eff * (i <= ntr & d$year >= 2019)
  d
}

# The gradient of (1/m) ||A w - b||^2 + z2 ||w||^2 (columns of A - b
# centred when icpt) at w, and the KKT check on the simplex: equal on the
# support, no smaller off it.
kkt <- function(a, b, w, z2 = 0, icpt = FALSE, tol = 1e-8) {
  d <- a - b
  if (icpt) d <- sweep(d, 2, colMeans(d))
  g <- 2 * crossprod(d, d %*% w) / nrow(a) + 2 * z2 * w
  on <- w > 1e-10
  s <- max(abs(g))
  expect_equal(sum(w), 1, tolerance = 1e-12)
  expect_true(all(w >= 0))
  expect_lt(diff(range(g[on])), tol * s)
  expect_true(all(g[!on] >= min(g[on]) - tol * s))
}

test_that("the synthetic control's weights solve their problem", {
  d <- scpanel(1)
  s <- synthc(d$rate, d$town, d$year, d$first)
  y <- matrix(d$rate, 15)
  kkt(y[1:9, -1], y[1:9, 1], s$weights$w)
  expect_equal(s$path$synth, drop(y[, -1] %*% s$weights$w),
               tolerance = 1e-12)
  expect_equal(s$fit$rmspepre, sqrt(mean(s$path$gap[1:9]^2)),
               tolerance = testthat_tolerance())
  expect_equal(s$fit$p, (1 + sum(s$plac$ratio >= s$fit$ratio)) / 10,
               tolerance = testthat_tolerance())
  # a placebo is the same method on that unit without the treated one
  k <- d$town != "u01"
  p3 <- synthc(d$rate[k], d$town[k], d$year[k],
               ifelse(d$town[k] == "u03", 2019, NA), plac = FALSE)
  expect_equal(s$plac$ratio[s$plac$unit == "u03"], p3$fit$ratio,
               tolerance = 1e-12)
  expect_true(is.na(p3$fit$p))
  expect_null(p3$plac)
  expect_output(print(s), "synthc: unit u01")
})

test_that("a convex combination of donors is recovered exactly", {
  d <- scpanel(2)
  y <- matrix(d$rate, 15)
  y[, 1] <- 0.3 * y[, 2] + 0.7 * y[, 5] - 2 * (2010:2024 >= 2019)
  d$rate <- as.vector(y)
  s <- synthc(d$rate, d$town, d$year, d$first)
  expect_equal(s$weights$w, c(0.3, 0, 0, 0.7, 0, 0, 0, 0, 0),
               tolerance = 1e-9)
  expect_equal(s$path$gap, -2 * (2010:2024 >= 2019), tolerance = 1e-9)
})

test_that("the synthetic DiD weights solve their problems; tau from them", {
  d <- scpanel(3, nu = 12, ntr = 3)
  r <- sdidx(d$rate, d$town, d$year, d$first)
  y <- matrix(d$rate, 15)
  y0 <- y[, 4:12]
  y1 <- rowMeans(y[, 1:3])
  sg <- sd(as.vector(diff(y0[1:9, ])))
  kkt(y0[1:9, ], y1[1:9], r$unitw$w, z2 = sqrt(3 * 6) * sg^2, icpt = TRUE)
  kkt(t(y0[1:9, ]), colMeans(y0[10:15, ]), r$timew$w, z2 = (1e-6 * sg)^2,
      icpt = TRUE, tol = 1e-6)
  w <- r$unitw$w
  l <- r$timew$w
  ch <- function(v) mean(v[10:15]) - sum(l * v[1:9])
  tau <- ch(y1) - sum(w * apply(y0, 2, ch))
  expect_equal(r$est$est, tau, tolerance = 1e-12)
  expect_identical(r$meth, "jackknife")
  # the jackknife with the weights held
  u <- c(vapply(1:9, function(j) {
    ww <- w[-j] / sum(w[-j])
    ch(y1) - sum(ww * apply(y0[, -j], 2, ch))
  }, 0), vapply(1:3, function(j) {
    ch(rowMeans(y[, setdiff(1:3, j)])) - sum(w * apply(y0, 2, ch))
  }, 0))
  expect_equal(r$est$se, sqrt(11 / 12 * sum((u - mean(u))^2)),
               tolerance = 1e-12)
  expect_equal(r$path$synth, drop(y0 %*% w), tolerance = 1e-12)
  expect_output(print(r), "jackknife")
})

test_that("one treated unit takes the placebo error, seeded", {
  d <- scpanel(4, nu = 8)
  a <- sdidx(d$rate, d$town, d$year, d$first, reps = 50, seed = 3)
  b <- sdidx(d$rate, d$town, d$year, d$first, reps = 50, seed = 3)
  expect_identical(a, b)
  expect_identical(a$meth, "placebo")
  set.seed(5)
  c1 <- sdidx(d$rate, d$town, d$year, d$first, reps = 20)
  set.seed(5)
  c2 <- sdidx(d$rate, d$town, d$year, d$first, reps = 20)
  expect_identical(c1, c2)
  # comparison units flat before adoption: no noise, no ridge, all the
  # weight on one unit, and without it no jackknife
  e <- scpanel(6, nu = 5, ntr = 2)
  k <- as.integer(e$town) > 2 & e$year < 2019
  e$rate[k] <- 10 + as.integer(e$town[k])
  r <- sdidx(e$rate, e$town, e$year, e$first)
  expect_equal(max(r$unitw$w), 1, tolerance = testthat_tolerance())
  expect_true(is.na(r$est$se))
  expect_true(is.na(r$est$lo))
})

test_that("the effect is recovered; the errors fit the spread", {
  # unit and period effects, unit trends and unit noise: the model of
  # Arkhangelsky et al. with an effect of -2 on five of twenty units
  est <- se <- numeric(100)
  for (k in seq_along(est)) {
    set.seed(500 + k)
    d <- expand.grid(year = 2010:2024, town = sprintf("u%02d", 1:20))
    i <- as.integer(d$town)
    d$first <- ifelse(i <= 5, 2019, NA)
    d$rate <- 20 + rnorm(20)[i] + rnorm(15)[d$year - 2009] +
      (d$year - 2010) * runif(20, 0.2, 0.8)[i] + rnorm(300) -
      2 * (i <= 5 & d$year >= 2019)
    r <- sdidx(d$rate, d$town, d$year, d$first)
    est[k] <- r$est$est
    se[k] <- r$est$se
  }
  expect_lt(abs(mean(est) + 2), 4 * sd(est) / sqrt(100))
  # 100 draws: the spread's own error is about 7%; the jackknife runs a
  # little conservative
  expect_gt(mean(se) / sd(est), 0.85)
  expect_lt(mean(se) / sd(est), 1.3)
  gp <- numeric(40)
  for (k in seq_along(gp)) {
    d <- scpanel(300 + k, nu = 12, eff = -2)
    gp[k] <- mean(synthc(d$rate, d$town, d$year, d$first,
                         plac = FALSE)$path$gap[10:15])
  }
  expect_lt(abs(mean(gp) + 2), 4 * sd(gp) / sqrt(40))
})

test_that("synthc() and sdidx() refuse what they cannot fit", {
  d <- scpanel(7)
  y <- d$rate
  u <- d$town
  t <- d$year
  f <- d$first
  expect_error(synthc(y[-1], u[-1], t[-1], f[-1]), "balanced")
  expect_error(synthc(y, u, t, ifelse(u %in% c("u01", "u02"), 2019, NA)),
               "one unit must be treated")
  expect_error(sdidx(y, u, t, c(2019, 2020)[match(u, c("u01", "u02"))]),
               "one adoption period")
  expect_error(sdidx(y, u, t, rep(NA_real_, 150)), "no unit is treated")
  expect_error(sdidx(y, u, t, ifelse(u == "u10", NA, 2019)),
               "two units never treated")
  expect_error(sdidx(y, u, t, ifelse(u == "u01", 2010, NA)),
               "after the first")
  expect_error(sdidx(y, u, t, ifelse(u == "u01", 2019.5, NA)),
               "after the first")
  expect_error(sdidx(y, u, t, f, reps = 1), "reps")
  expect_error(synthc(y, u, t, f, plac = NA), "plac")
  y2 <- y
  y2[u != "u01"] <- 0
  expect_error(synthc(y2, u, t, f), NA)
})
