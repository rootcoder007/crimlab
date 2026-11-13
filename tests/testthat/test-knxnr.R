# The Knox test of near repeats. The kernel is checked in src/t_knx.cpp
# (brute force, and the expectation as the mean over every order of the
# times); here the binding against the table written out in R.

knox <- function(p, t, ds, dt) {
  d <- as.matrix(dist(p))
  e <- abs(outer(t, t, "-"))
  up <- upper.tri(d)
  s <- findInterval(d[up], c(0, ds), left.open = TRUE)
  s[d[up] == 0] <- 1L
  b <- findInterval(e[up], c(0, dt), left.open = TRUE)
  b[e[up] == 0] <- 1L
  s[s > length(ds)] <- 0L
  b[b > length(dt)] <- 0L
  ok <- s > 0 & b > 0
  obs <- tabulate(s[ok] + (b[ok] - 1L) * length(ds), length(ds) * length(dt))
  np <- length(t) * (length(t) - 1) / 2
  ex <- as.vector(outer(tabulate(s[s > 0], length(ds)),
                        tabulate(b[b > 0], length(dt))) / np)
  list(obs = obs, exp = ex, pairs = sum(s > 0))
}

test_that("the Knox table and its expectation follow the definition", {
  set.seed(51)
  n <- 200
  p <- cbind(x = round(runif(n, 0, 1000)), y = round(runif(n, 0, 1000)))
  t <- round(runif(n, 0, 120))
  ds <- c(50, 100, 200)
  dt <- c(7, 14, 30)
  k <- knxnr(p, t = t, ds = ds, dt = dt, nsim = 0)
  r <- knox(p, t, ds, dt)
  expect_s3_class(k, "knxnr")
  expect_named(k, c("tab", "n", "pairs", "nsim"))
  expect_identical(k$tab$obs, as.numeric(r$obs))
  expect_equal(k$tab$exp, r$exp, tolerance = 1e-14)
  expect_identical(k$pairs, as.numeric(r$pairs))
  expect_identical(k$tab$dhi, rep(ds, 3))
  expect_identical(k$tab$tlo, rep(c(0, dt[1:2]), each = 3))
  expect_equal(k$tab$ratio, ifelse(r$exp > 0, r$obs / r$exp, NaN),
               tolerance = 1e-14)
  expect_true(all(is.nan(k$tab$p)))
  expect_equal(k$tab$ppois,
               ifelse(r$obs > 0,
                      ppois(r$obs - 1, r$exp, lower.tail = FALSE), 1),
               tolerance = 1e-14)
  expect_output(print(k), "knxnr: 200 incidents")
})

test_that("permutation p-values find a near-repeat cluster", {
  set.seed(52)
  n <- 150
  p <- cbind(x = runif(n, 0, 5000), y = runif(n, 0, 5000))
  t <- runif(n, 0, 365)
  # twenty repeats a few metres and days after twenty incidents
  q <- p[1:20, ] + runif(40, -20, 20)
  p2 <- rbind(p, q)
  t2 <- c(t, t[1:20] + runif(20, 0, 5))
  k <- knxnr(p2, t = t2, ds = c(50, 500), dt = c(7, 60), nsim = 199,
             seed = 3)
  expect_identical(k, knxnr(p2, t = t2, ds = c(50, 500), dt = c(7, 60),
                            nsim = 199, seed = 3))
  expect_lt(k$tab$p[1], 0.01)
  expect_gt(k$tab$ratio[1], 5)
  expect_true(all(k$tab$p >= 1 / 200 & k$tab$p <= 1))
  expect_true(all(k$tab$p[k$tab$obs == 0] == 1))
  set.seed(4)
  a <- knxnr(p2, t = t2, ds = c(50, 500), dt = c(7, 60), nsim = 19)
  set.seed(4)
  expect_identical(a, knxnr(p2, t = t2, ds = c(50, 500), dt = c(7, 60),
                            nsim = 19))
})

test_that("dates, date-times and longitude and latitude are read", {
  set.seed(53)
  p <- cbind(x = runif(40, 0, 3000), y = runif(40, 0, 3000))
  d <- as.Date("2024-01-01") + sample(0:200, 40, replace = TRUE)
  k <- knxnr(p, t = as.numeric(d), ds = 500, dt = 10, nsim = 0)
  expect_identical(knxnr(p, t = d, ds = 500, dt = 10, nsim = 0), k)
  tm <- as.POSIXct(d) + 3600
  expect_identical(knxnr(p, t = tm, ds = 500, dt = 10, nsim = 0)$tab$obs,
               k$tab$obs)
  ll <- cbind(lon = -79.4 + p[, 1] / 1e5, lat = 43.7 + p[, 2] / 1e5)
  q <- eqproj(ll)
  expect_identical(knxnr(ll, t = d, ds = 500, dt = 10, nsim = 0,
                         lonlat = TRUE)$tab$obs,
                   knxnr(q, t = d, ds = 500, dt = 10, nsim = 0)$tab$obs)
})

test_that("knxnr() refuses bad input in words", {
  p <- cbind(x = 1:5, y = c(2, 1, 4, 3, 5))
  errs("knxnr", list(
    list("at least 3 incidents", p[1:2, ], t = 1:2, ds = 1, dt = 1),
    list("`ds` must increase", p, t = 1:5, ds = c(2, 1), dt = 1),
    list("`dt` must be finite numbers", p, t = 1:5, ds = 1,
         dt = numeric(0)),
    list("`t`", p, t = 1:4, ds = 1, dt = 1),
    list("`dt`", p, t = 1:5, ds = 1, dt = 0),
    list("`nsim`", p, t = 1:5, ds = 1, dt = 1, nsim = -1)
  ))
})
