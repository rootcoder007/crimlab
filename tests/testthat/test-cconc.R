# cconc() and hilltl(): against the mean-difference form of the Gini
# coefficient, sorting and cumulative sums, an exact enumeration of the
# Poisson null, lm(), a zeta function summed directly with optimize(), and
# the refusals.

gmd <- function(x) sum(abs(outer(x, x, "-"))) / (2 * length(x)^2 * mean(x))

test_that("cconc() gives the Gini, its split, top shares and the curve", {
  set.seed(8)
  x <- rnbinom(300, mu = 1.5, size = 0.5)
  r <- cconc(x, fracs = c(0.05, 0.1, 0.5), nsim = 0)
  s <- r$stat
  expect_equal(s[["gini"]], gmd(x), tolerance = 1e-13)
  expect_equal(s[["gpos"]], gmd(x[x > 0]), tolerance = 1e-13)
  expect_equal(s[["zero"]], mean(x == 0), tolerance = 1e-15)
  expect_equal(
    s[["gini"]],
    s[["zero"]] + (1 - s[["zero"]]) * s[["gpos"]],
    tolerance = 1e-14
  )
  expect_equal(s[["var"]], var(x), tolerance = 1e-13)
  expect_equal(s[["disp"]], var(x) / mean(x), tolerance = 1e-13)
  expect_equal(s[["ivar"]], max(var(x) - mean(x), 0), tolerance = 1e-13)
  expect_equal(
    s[["cvb"]],
    sqrt(mean((x - mean(x))^2) / 2) / mean(x),
    tolerance = 1e-13
  )
  expect_equal(s[["nullb"]], 1 / sqrt(2 * mean(x)), tolerance = 1e-15)
  expect_equal(s[["znull"]], exp(-mean(x)), tolerance = 1e-15)
  expect_lte(s[["gini"]], s[["cvb"]])
  d <- sort(x, decreasing = TRUE)
  k <- ceiling(c(0.05, 0.1, 0.5) * 300)
  expect_identical(r$top$units, as.numeric(k))
  expect_equal(r$top$share, cumsum(d)[k] / sum(x), tolerance = 1e-14)
  expect_equal(
    r$lorenz$share,
    c(0, cumsum(sort(x)) / sum(x)),
    tolerance = 1e-14
  )
  expect_equal(r$lorenz$pop, (0:300) / 300, tolerance = 1e-15)
  expect_null(r$null)
  expect_null(r$by)
  expect_output(print(r), "Gini 0\\.")
})

test_that("weights count as repeated units", {
  v <- c(1, 2, 4, 9)
  w <- c(50, 20, 6, 2)
  r <- cconc(v, w = w, fracs = c(0.1, 0.25), nsim = 0)
  e <- cconc(rep(v, w), fracs = c(0.1, 0.25), nsim = 0)
  expect_equal(r$stat, e$stat, tolerance = 1e-13)
  expect_equal(r$top, e$top, tolerance = 1e-14)
  # the Gini by group, on the units of each group
  g <- c("a", "a", "b", "b")
  b <- cconc(v, w = w, by = g, nsim = 0)$by
  expect_identical(b$group, c("a", "b"))
  expect_equal(b$gini[1], gmd(rep(v[1:2], w[1:2])), tolerance = 1e-13)
  expect_identical(b$units, c(70, 8))
  z <- cconc(c(0, 0, 3), by = c("p", "p", "q"), nsim = 0)$by
  expect_identical(z$gini, c(NA, 0))
  expect_identical(z$total, c(0, 3))
})

test_that("the Poisson null matches an exact enumeration", {
  # E[Gini] of three Poisson(1) counts (0 when all are 0), summed exactly
  k <- 0:25
  p <- dpois(k, 1)
  g <- expand.grid(a = k, b = k, c = k)
  pr <- p[g$a + 1] * p[g$b + 1] * p[g$c + 1]
  gi <- apply(as.matrix(g), 1, function(v) if (sum(v)) gmd(v) else 0)
  eg <- sum(pr * gi)
  sg <- sqrt(sum(pr * gi^2) - eg^2)
  r <- cconc(c(0, 1, 2), nsim = 20000, seed = 3)
  expect_lt(abs(r$null[["mean"]] - eg), 4 * sg / sqrt(20000))
  expect_equal(
    r$null[["excess"]],
    r$stat[["gini"]] - r$null[["mean"]],
    tolerance = 1e-15
  )
  expect_true(r$null[["p"]] > 0 && r$null[["p"]] <= 1)
  expect_identical(
    cconc(c(0, 1, 2), nsim = 50, seed = 3)$null,
    cconc(c(0, 1, 2), nsim = 50, seed = 3)$null
  )
  set.seed(1)
  a <- cconc(c(0, 1, 2), nsim = 50)$null
  set.seed(1)
  expect_identical(cconc(c(0, 1, 2), nsim = 50)$null, a)
  expect_output(print(r), "simulated null")
})

test_that("places of events give counts and the growth of distinct places", {
  set.seed(6)
  pl <- sample(letters, 120, TRUE, prob = 26:1)
  r <- cconc(pl, by = rep(1:2, 60), nsim = 0)
  tab <- table(pl)
  expect_identical(r$places, names(tab))
  expect_equal(r$stat[["gini"]], gmd(as.numeric(tab)), tolerance = 1e-13)
  cv <- r$grow$curve
  dist <- cumsum(!duplicated(pl))
  expect_identical(cv$dist, as.numeric(dist))
  m <- r$grow$mhat
  expect_equal(m * sum(1 / (m + 0:119)), dist[120], tolerance = 1e-13)
  expect_equal(cv$expd, m * cumsum(1 / (m + 0:119)), tolerance = 1e-12)
  expect_equal(cv$lo, m * log((m + 1:120) / m), tolerance = 1e-12)
  expect_equal(cv$hi, 1 + m * log((m + 0:119) / m), tolerance = 1e-12)
  h <- 60:120
  expect_equal(
    r$grow$slope,
    unname(coef(lm(log(dist[h]) ~ log(h)))[2]),
    tolerance = 1e-12
  )
  # the Gini within each group, of the places counted in that group
  g1 <- table(pl[c(TRUE, FALSE)])
  expect_equal(r$by$gini[1], gmd(as.numeric(g1)), tolerance = 1e-13)
  # every event at a new place: no finite M
  a <- cconc(c("a", "b", "c"), nsim = 0)$grow
  expect_identical(a$mhat, Inf)
  expect_true(all(is.na(a$curve$expd)) && is.na(a$slope))
})

test_that("the urban-scaling exponent is least squares on the logs", {
  set.seed(9)
  pop <- round(runif(40, 1e3, 1e5))
  x <- rpois(40, 2e-3 * pop^1.1)
  x[3] <- 0
  r <- cconc(x, pop = pop, nsim = 0)$scale
  f <- summary(lm(log(x[-3]) ~ log(pop[-3])))
  expect_equal(r[["beta"]], coef(f)[2, 1], tolerance = 1e-12)
  expect_equal(r[["se"]], coef(f)[2, 2], tolerance = 1e-10)
  expect_equal(r[["pref"]], exp(coef(f)[1, 1]), tolerance = 1e-11)
  expect_equal(r[["r2"]], f$r.squared, tolerance = 1e-12)
  expect_identical(r[["used"]], 39)
  few <- cconc(c(1, 0, 2), pop = 1:3, nsim = 0)$scale
  expect_identical(unname(few), c(NA, NA, NA, NA, 2))
})

test_that("cconc() refuses bad input in words", {
  errs(
    "cconc",
    list(
      list("`x` must be at least 0", c(1, -1)),
      list("every count is 0", c(0, 0), nsim = 0),
      list("at least two events", "a"),
      list("`w` goes with counts", c("a", "b"), w = 1:2),
      list("`pop` goes with counts", c("a", "b"), pop = 1:2),
      list("`w` must be greater than 0", 1:2, w = c(1, 0)),
      list("`fracs` must be shares from 0 to 1", 1:2, fracs = 2),
      list("`nsim` must be one whole number", 1:2, nsim = -1),
      list("whole-number weights", 1:2, w = c(1, 1.5), nsim = 10),
      list("mean count is above 1e8", c(1e9, 1e9), nsim = 1, seed = 1),
      list("`by` must have one label per value", 1:2, by = 1, nsim = 0)
    )
  )
})

# the Hurwitz zeta summed term by term to 1e6, smallest terms first, with
# the integral beyond and the half term at the cut. Largest first, the sum
# is right to 1e-12 only where R adds in extended precision (not on arm64).
zsum <- function(s, q) {
  k <- 999999:0
  sum((q + k)^-s) + (q + 1e6)^(1 - s) / (s - 1) + (q + 1e6)^-s / 2
}

test_that("the zeta kernel agrees with closed forms and direct sums", {
  z <- function(s, q) .Call(C_clab_r_hzta, s, q)
  expect_equal(z(c(2, 4), 1), c(pi^2 / 6, pi^4 / 90), tolerance = 1e-15)
  expect_equal(z(2, 0.5), pi^2 / 2, tolerance = 1e-15)
  for (s in c(1.5, 2.3, 7)) {
    for (q in c(0.3, 1, 4.5, 60)) {
      expect_equal(z(s, q), zsum(s, q), tolerance = 1e-12,
                   info = paste(s, q))
    }
  }
  expect_error(z(1, 1), "greater than 1")
})

test_that("hilltl() matches the estimators written out", {
  set.seed(12)
  y <- 2 * runif(300)^(-1 / 1.7)
  r <- hilltl(y, xmin = 3, discrete = FALSE)
  tl <- y[y >= 3]
  a <- 1 + length(tl) / sum(log(tl / 3))
  expect_equal(r$alpha, a, tolerance = 1e-14)
  expect_equal(r$se, (a - 1) / sqrt(length(tl)), tolerance = 1e-14)
  u <- sort(tl)
  f <- 1 - (u / 3)^(1 - a)
  n <- length(u)
  expect_equal(
    r$ks,
    max(abs(seq_len(n) / n - f), abs((seq_len(n) - 1) / n - f)),
    tolerance = 1e-14
  )
  expect_identical(r$ntl, n)
  expect_true(r$reliable)
  expect_output(print(r), "continuous Hill estimator")
  # counts: the closed form with xmin - 1/2
  k <- floor(runif(500)^(-1 / 1.4))
  ap <- hilltl(k, xmin = 2, approx = TRUE)
  tk <- k[k >= 2]
  expect_equal(ap$alpha, 1 + length(tk) / sum(log(tk / 1.5)), tolerance = 1e-14)
})

test_that("the discrete estimate maximises the zeta likelihood", {
  set.seed(13)
  k <- floor(runif(400)^(-1 / 1.5))
  r <- hilltl(k, xmin = 1)
  ls <- sum(log(k))
  nll <- function(a) length(k) * log(zsum(a, 1)) + a * ls
  o <- optimize(nll, c(1.0001, 25), tol = 1e-10)$minimum
  expect_equal(r$alpha, o, tolerance = 1e-7)
  h <- 1e-4
  d2 <- (nll(r$alpha + h) - 2 * nll(r$alpha) + nll(r$alpha - h)) / h^2
  expect_equal(r$se, sqrt(1 / d2), tolerance = 1e-4)
  u <- sort(unique(k))
  f <- 1 - vapply(u + 1, function(v) zsum(r$alpha, v), 0) / zsum(r$alpha, 1)
  expect_equal(
    r$ks,
    max(abs(cumsum(table(k)) / length(k) - f)),
    tolerance = 1e-9
  )
  expect_identical(r$xmin, 1)
  expect_match(r$method, "^discrete maximum likelihood$")
  # the default threshold is the 50th largest value
  expect_identical(hilltl(k)$xmin, sort(k, decreasing = TRUE)[50])
  expect_identical(hilltl(c(3, 5, 9))$xmin, 3)
})

test_that("hilltl() says when it cannot fit", {
  few <- hilltl(c(1, 2, 3, 50), xmin = 10)
  expect_true(is.na(few$alpha) && is.na(few$ks) && !few$reliable)
  expect_identical(few$method, "too few tail observations")
  eq <- hilltl(c(1, 5, 5, 5), xmin = 5)
  expect_match(eq$method, "every tail value is the same")
  edge <- hilltl(c(1001, rep(1000, 999)), xmin = 1000)
  expect_match(edge$method, "at an end of the search")
  expect_false(edge$reliable)
  expect_output(print(edge), "not reliable")
  errs(
    "hilltl",
    list(
      list("no values above 0", c(0, 0)),
      list("`discrete` must be TRUE or FALSE", 1:3, discrete = NA),
      list("`xmin` must be one number", 1:3, xmin = 0),
      list("`mintl` must be one whole number", 1:3, mintl = 0)
    )
  )
})
