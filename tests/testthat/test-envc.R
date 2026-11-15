# The environmental-crime family: lmdix() against the LMDI written out
# here on hand cases and simulated chains, fepnl() against lm() with unit
# and period dummies and a cluster sandwich written out here, its VIFs,
# slopes and Jarque-Bera statistic by hand, pcdtst() against the pairwise
# correlations of base R, and the refusals.

# a simulated chain: cases = (cases / P) (P / S) S over 5 units, 4 years
envsim <- function(seed = 4) {
  set.seed(seed)
  d <- expand.grid(year = 2016:2019, unit = paste0("u", 1:5))
  d$P <- runif(20, 1, 10)
  d$S <- runif(20, 5, 50)
  d$cases <- rpois(20, 20) + 1
  d
}

# L(a, b) ln(f1 / f0), by hand
lm1 <- function(y0, y1, f0, f1) {
  l <- ifelse(y0 == y1, y0, (y1 - y0) / (log(y1) - log(y0)))
  l * log(f1 / f0)
}

test_that("a two-factor change splits by the log mean", {
  r <- lmdix(c(6, 12, 10, 20), cbind(P = c(3, 4, 5, 5)), c("a", "a", "b", "b"),
             c(1, 2, 1, 2))
  d <- r$dec
  expect_identical(r$factors, c("y/P", "P"))
  l <- 6 / log(2)
  expect_equal(d[["y/P"]], c(l * log(1.5), 10), tolerance = 1e-12)
  expect_equal(d[["P"]], c(l * log(4 / 3), 0), tolerance = 1e-12)
  expect_identical(d$dy, c(6, 10))
  expect_identical(d$from, c(1, 1))
  expect_identical(d$to, c(2, 2))
})

test_that("temporal contributions add up to the change, interval by interval", {
  d <- envsim()
  r <- lmdix(d$cases, d[c("P", "S")], d$unit, d$year,
             at = c(2016, 2018, 2019), fnm = c("rate", "dens", "area"))
  e <- r$dec
  expect_identical(nrow(e), 10L)
  expect_equal(rowSums(e[c("rate", "dens", "area")]), e$dy,
               tolerance = 1e-12)
  a <- d[d$year == 2016, ]
  b <- d[d$year == 2018, ]
  expect_equal(e$rate[e$from == 2016],
               lm1(a$cases, b$cases, a$cases / a$P, b$cases / b$P),
               tolerance = 1e-12)
  expect_equal(e$dens[e$from == 2016], lm1(a$cases, b$cases, a$P / a$S,
                                           b$P / b$S), tolerance = 1e-12)
  # the factor values are the rates of the chain
  expect_equal(r$fac$rate, d$cases / d$P, tolerance = 1e-15)
  expect_equal(r$fac$area, d$S, tolerance = 1e-15)
  expect_null(r$ref)
})

test_that("the spatial gap is to the region of mean count and quantities", {
  d <- envsim()
  r <- lmdix(d$cases, d[c("P", "S")], d$unit, d$year, type = "space",
             at = 2017)
  e <- r$dec
  s <- d[d$year == 2017, ]
  m <- colMeans(s[c("cases", "P", "S")])
  expect_equal(e$ref, rep(m[["cases"]], 5), tolerance = 1e-14)
  expect_equal(e$dy, s$cases - m[["cases"]], tolerance = 1e-12)
  expect_equal(rowSums(e[r$factors]), e$dy, tolerance = 1e-12)
  expect_equal(e[["P/S"]], lm1(rep(m[["cases"]], 5), s$cases,
                               m[["P"]] / m[["S"]], s$P / s$S),
               tolerance = 1e-12)
  expect_equal(unlist(r$ref[1, r$factors]), c(m[["cases"]] / m[["P"]],
                                       m[["P"]] / m[["S"]], m[["S"]]),
               tolerance = 1e-14, ignore_attr = TRUE)
  g <- lmdix(d$cases, d[c("P", "S")], d$unit, d$year, type = "space",
             ref = "gmean")
  expect_equal(g$ref$y[1], exp(mean(log(d$cases[d$year == 2016]))),
               tolerance = 1e-14)
  expect_identical(nrow(g$dec), 20L)
})

test_that("zeros give NA, or are replaced by `zero`", {
  y <- c(0, 5, 4, 8)
  q <- cbind(P = c(2, 3, 4, 4))
  u <- c("a", "a", "b", "b")
  r <- lmdix(y, q, u, c(1, 2, 1, 2))
  expect_true(all(is.na(r$dec[1, r$factors])))
  expect_equal(r$dec[["y/P"]][2], 4, tolerance = 1e-12)
  z <- lmdix(y, q, u, c(1, 2, 1, 2), zero = 1e-8)
  expect_equal(sum(z$dec[1, z$factors]), 5 - 1e-8, tolerance = 1e-12)
  # the replaced zero's factor takes most, not all, of the change
  expect_gt(z$dec[["y/P"]][1], 4.5)
})

# CR1 by hand: (X'X)^-1 (sum_g s_g s_g') (X'X)^-1 G / (G - 1) (n - 1) /
# (n - K) on the dummy regression
envcr1 <- function(m, cl, k) {
  x <- model.matrix(m)
  r <- residuals(m)
  b <- solve(crossprod(x))
  s <- rowsum(x * r, cl)
  g <- length(unique(cl))
  n <- nrow(x)
  b %*% crossprod(s) %*% b * g / (g - 1) * (n - 1) / (n - k)
}

envpan <- function(seed = 5) {
  set.seed(seed)
  d <- expand.grid(year = 2016:2023, prov = paste0("p", 1:12))
  d$cper <- runif(96, 0.005, 0.04)
  d$ur <- runif(96, 0.01, 0.03)
  d$eipa <- rnorm(96, 5)
  d$cases <- 50 + 9000 * d$cper + 1500 * d$ur - 4e5 * d$cper * d$ur +
    as.integer(substring(d$prov, 2)) + (d$year - 2016) + rnorm(96, sd = 5)
  d
}

test_that("fepnl() is lm() with dummies and its error CR1", {
  d <- envpan()
  r <- fepnl(d$cases, d[c("cper", "eipa", "ur")], d$prov, d$year,
             int = c("cper", "ur"))
  m <- lm(cases ~ cper + eipa + ur + I(cper * ur) + factor(prov) +
            factor(year), d)
  expect_equal(r$coef$est, unname(coef(m)[2:5]), tolerance = 1e-9)
  expect_identical(r$coef$term, c("cper", "eipa", "ur", "cper:ur"))
  # the province effects are nested in the clusters: K = 4 + 1 + 7
  v <- envcr1(m, d$prov, 12)[2:5, 2:5]
  expect_equal(unname(r$vcov), unname(v), tolerance = 1e-8)
  expect_identical(r$fit$k, 12)
  tq <- qt(0.975, 11)
  expect_equal(r$coef$lo, r$coef$est - tq * r$coef$se, tolerance = 1e-14)
  expect_equal(r$coef$p, 2 * pt(-abs(r$coef$est / r$coef$se), 11),
               tolerance = 1e-12)
  expect_equal(r$fit$r2, summary(m)$r.squared, tolerance = 1e-10)
  expect_equal(r$resid, unname(residuals(m)), tolerance = 1e-8)
  x <- cbind(d$cper, d$eipa, d$ur, d$cper * d$ur)
  expect_equal(r$coef$beta, r$coef$est * apply(x, 2, sd) / sd(d$cases),
               tolerance = 1e-12)
  vf <- vapply(1:4, function(j) {
    1 / (1 - summary(lm(x[, j] ~ x[, -j]))$r.squared)
  }, 0)
  expect_equal(r$coef$vif, vf, tolerance = 1e-9)
  s <- r$slopes
  mu <- mean(d$ur) + c(-1, 0, 1) * sd(d$ur)
  b <- r$coef$est
  expect_equal(s$at, mu, tolerance = 1e-14)
  expect_equal(s$slope, b[1] + b[4] * mu, tolerance = 1e-12)
  vc <- r$vcov
  expect_equal(s$se, sqrt(vc[1, 1] + mu^2 * vc[4, 4] + 2 * mu * vc[1, 4]),
               tolerance = 1e-12)
  e <- residuals(m) - mean(residuals(m))
  sk <- mean(e^3) / mean(e^2)^1.5
  ku <- mean(e^4) / mean(e^2)^2
  jb <- 96 / 6 * (sk^2 + (ku - 3)^2 / 4)
  expect_equal(r$jb[["jb"]], jb, tolerance = 1e-6)
  expect_equal(r$jb[["p"]], pchisq(jb, 2, lower.tail = FALSE),
               tolerance = 1e-6)
})

test_that("fepnl() options: no factor, slopes where asked, no interaction", {
  d <- envpan()
  r <- fepnl(d$cases, d[c("cper", "ur")], d$prov, d$year,
             int = c("cper", "ur"), at = c(0.01, 0.02))
  n <- fepnl(d$cases, d[c("cper", "ur")], d$prov, d$year,
             int = c("cper", "ur"), ssc = "none")
  f <- 12 / 11 * 95 / (96 - 11)
  expect_equal(n$vcov * f, r$vcov, tolerance = 1e-12)
  expect_identical(r$slopes$at, c(0.01, 0.02))
  p <- fepnl(d$cases, d$cper, d$prov, d$year)
  expect_null(p$slopes)
  expect_identical(p$coef$term, "x1")
  expect_identical(p$coef$vif, 1)
})

test_that("pcdtst() is the CD of the pairwise correlations", {
  set.seed(6)
  d <- expand.grid(year = 1:9, prov = paste0("p", 1:6))
  d$e <- rnorm(54) + rnorm(9)[d$year]
  d <- d[-c(3, 20), ]
  r <- pcdtst(d$e, d$prov, d$year)
  w <- tapply(d$e, list(d$prov, d$year), identity)
  rho <- cor(t(w), use = "pairwise.complete.obs")
  tij <- crossprod(t(!is.na(w)))
  k <- lower.tri(rho)
  cd <- sum(sqrt(tij[k]) * rho[k]) / sqrt(sum(k))
  expect_equal(r[["cd"]], cd, tolerance = 1e-12)
  expect_equal(r[["p"]], 2 * pnorm(-abs(cd)), tolerance = 1e-12)
  expect_equal(r[["rho"]], mean(rho[k]), tolerance = 1e-12)
  expect_equal(r[["absrho"]], mean(abs(rho[k])), tolerance = 1e-12)
  expect_identical(r[["pairs"]], 15)
})

test_that("the family's refusals name the function", {
  d <- envsim()
  q <- d[c("P", "S")]
  expect_error(lmdix(-d$cases, q, d$unit, d$year), "lmdix: .*at least 0")
  expect_error(lmdix(d$cases, q[-1, ], d$unit, d$year), "one row per row")
  expect_error(lmdix(d$cases, q, d$unit, d$year, at = 2030), "periods of")
  expect_error(lmdix(d$cases, q, d$unit, d$year, at = 2016), "at least 2")
  expect_error(lmdix(d$cases, q, d$unit, d$year, zero = -1), "`zero`")
  expect_error(lmdix(d$cases, q, d$unit, d$year, fnm = "a"), "3 names")
  expect_error(lmdix(d$cases, q, d$unit, d$year, type = "x"), "`type`")
  p <- envpan()
  x <- p[c("cper", "ur")]
  expect_error(fepnl(p$cases, x, p$prov, p$year, int = "cper"),
               "fepnl: .*two names")
  expect_error(fepnl(p$cases, x, p$prov, p$year, int = c("cper", "z")),
               "not a column")
  expect_error(fepnl(p$cases, x, p$prov, p$year, at = 1), "give `int`")
  expect_error(fepnl(p$cases, cbind(x, k = 1), p$prov, p$year),
               "absorbed")
  expect_warning(fepnl(p$cases, x, p$prov, p$year,
                       clus = rep(1:4, each = 24)), "ten clusters")
  expect_error(pcdtst(rep(1, 4), c(1, 1, 2, 2), c(1, 2, 1, 2)),
               "pcdtst: .*undefined")
})

test_that("print and plot return their input and what they drew", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  d <- envsim()
  r <- lmdix(d$cases, d[c("P", "S")], d$unit, d$year)
  expect_identical(capture.output(x <- print(r))[1],
                   "lmdix: time decomposition, 5 rows, 3 factors")
  expect_identical(x, r)
  p <- plot(r)
  expect_equal(p$value, unname(colSums(r$dec[r$factors])),
               tolerance = 1e-14)
  expect_identical(plot(r, units = "u1")$value,
                   unname(unlist(r$dec[1, r$factors])))
  e <- envpan()
  f <- fepnl(e$cases, e[c("cper", "ur")], e$prov, e$year,
             int = c("cper", "ur"))
  expect_output(y <- print(f), "simple slopes")
  expect_identical(y, f)
  g <- plot(f)
  expect_identical(g$coef$est, f$coef$est)
  expect_identical(g$slopes$est, f$slopes$slope)
  expect_error(plot(structure(list(), class = "fepnl")), "not a complete")
})
