# didfe(): the least-squares fit against lm() with the unit and period
# dummies and a cluster sandwich written out here, the clustered errors and
# the Poisson fit against 'fixest' (values stored), the negative binomial
# against MASS::glm.nb(), the decomposition against 'bacondecomp' (values
# stored), the matching by hand, recovery of a known effect with errors
# against the Monte Carlo spread, and the refusals.

pd <- pepanel()

# CR1 by hand: (X'X)^-1 (sum_g s_g s_g') (X'X)^-1 G / (G - 1) (n - 1) /
# (n - K) on the dummy regression
cr1 <- function(m, cl, k) {
  x <- model.matrix(m)
  r <- residuals(m)
  w <- if (is.null(weights(m))) 1 else weights(m)
  b <- solve(crossprod(x * sqrt(w)))
  s <- rowsum(x * r * w, cl)
  g <- length(unique(cl))
  n <- nrow(x)
  b %*% crossprod(s) %*% b * g / (g - 1) * (n - 1) / (n - k)
}

test_that("the estimate is lm()'s with dummies, its error CR1", {
  r <- didfe(pd$rate, pd$town, pd$year, pd$on, adj = pd["x"])
  m <- lm(rate ~ on + x + factor(town) + factor(year), pd)
  expect_equal(r$coef$est, unname(coef(m)[2:3]), tolerance = 1e-10)
  # the unit effects are nested in the clusters: K = 2 + 1 + 9
  v <- cr1(m, pd$town, 12)
  expect_equal(r$coef$se, sqrt(diag(v))[2:3], tolerance = 1e-10,
               ignore_attr = TRUE)
  expect_identical(r$fit$k, 12)
  expect_identical(r$fit$ncl, 30)
  expect_equal(unname(r$vcov), unname(v[2:3, 2:3]), tolerance = 1e-10)
  tq <- qt(0.975, 29)
  expect_equal(r$coef$lo, r$coef$est - tq * r$coef$se, tolerance = 1e-14)
  expect_equal(r$coef$p, 2 * pt(-abs(r$coef$est / r$coef$se), 29),
               tolerance = 1e-12)
  # the comparison mean: the never-treated towns before 2018
  nv <- is.na(pd$first) & pd$year < 2018
  expect_equal(r$fit$ctlmean, mean(pd$rate[nv]), tolerance = 1e-14)
  expect_equal(r$fit$releff, r$coef$est[1] / mean(pd$rate[nv]),
               tolerance = 1e-14)
})

test_that("the clustered errors are fixest's, one-way and two-way", {
  # fixest 0.14.2, dev/ref-pe.R: feols(rate ~ on + x | town + year,
  # cluster = ~region) and cluster = ~region + year
  expect_warning(
    r1 <- didfe(pd$rate, pd$town, pd$year, pd$on, adj = pd["x"],
                clus = pd$region),
    NA
  )
  expect_equal(r1$coef$se, c(0.308144194504073, 1.07609354297877),
               tolerance = 1e-10)
  r2 <- didfe(pd$rate, pd$town, pd$year, pd$on, adj = pd["x"],
              clus = pd[c("region", "year")])
  expect_equal(r2$coef$se, c(0.546417553269131, 0.958445746214596),
               tolerance = 1e-10)
  expect_identical(r2$fit$ncl, 10)
  # weighted, fixest 0.14.2: feols(rate ~ on | town + year, weights = ~w)
  r3 <- didfe(pd$rate, pd$town, pd$year, pd$on, w = pd$w)
  expect_equal(c(r3$coef$est, r3$coef$se),
               c(-3.1261554539309, 0.482271447513932), tolerance = 1e-10)
})

test_that("an unbalanced panel, period effects only, and zero weights", {
  keep <- !(pd$town %in% c(3, 7) & pd$year %in% c(2016, 2021)) &
    pd$town != 30 | pd$year > 2019
  u <- pd[keep, ]
  r <- didfe(u$rate, u$town, u$year, u$on)
  m <- lm(rate ~ on + factor(town) + factor(year), u)
  expect_equal(r$coef$est, unname(coef(m)[2]), tolerance = 1e-9)
  expect_equal(r$coef$se, sqrt(cr1(m, u$town, 11)[2, 2]), tolerance = 1e-9)
  rt <- didfe(pd$rate, pd$town, pd$year, pd$on, fe = "time")
  mt <- lm(rate ~ on + factor(year), pd)
  expect_equal(rt$coef$est, unname(coef(mt)[2]), tolerance = 1e-10)
  # K = 1 + 1 + 9 + 29 (the units are not effects now; clusters by town)
  expect_identical(rt$fit$k, 11)
  w0 <- as.numeric(pd$town != 5)
  r0 <- didfe(pd$rate, pd$town, pd$year, pd$on, w = w0)
  r5 <- didfe(pd$rate[pd$town != 5], pd$town[pd$town != 5],
              pd$year[pd$town != 5], pd$on[pd$town != 5])
  expect_equal(r0$coef, r5$coef, tolerance = 1e-12)
  expect_identical(r0$fit$n, 290)
})

test_that("the Poisson fit is fixest's, the negative binomial glm.nb()'s", {
  # fixest 0.14.2: fepois(cnt ~ on | town + year, offset = ~log(pop),
  # cluster = ~town), and cluster = ~region + year without the offset
  r <- didfe(pd$cnt, pd$town, pd$year, pd$on, fam = "pois", expo = pd$pop)
  expect_equal(c(r$coef$est, r$coef$se),
               c(0.277910537133839, 0.0997874163832055), tolerance = 1e-7)
  expect_equal(r$coef$rr, exp(r$coef$est), tolerance = testthat_tolerance())
  expect_equal(r$fit$releff, exp(r$coef$est) - 1,
               tolerance = testthat_tolerance())
  r2 <- didfe(pd$cnt, pd$town, pd$year, pd$on, fam = "pois",
              clus = list(pd$region, pd$year))
  expect_equal(c(r2$coef$est, r2$coef$se),
               c(0.277910537133839, 0.0729701185450552), tolerance = 1e-7)
  z <- qnorm(0.975)
  expect_equal(r$coef$rrhi, exp(r$coef$est + z * r$coef$se),
               tolerance = testthat_tolerance())
  skip_if_not_installed("MASS")
  # counts with real overdispersion (size 3): theta is identified, so it can
  # be compared; on the Poisson counts above the likelihood is flat in theta
  # and its maximiser, in the tens of thousands, differs between platforms
  set.seed(7)
  od <- pd
  od$cnt <- rnbinom(nrow(od), mu = exp(2 + od$town / 50 + 0.3 * od$on),
                    size = 3)
  nb <- MASS::glm.nb(
    cnt ~ on + factor(town) + factor(year) + offset(log(pop)), od
  )
  rn <- didfe(od$cnt, od$town, od$year, od$on, fam = "nb", expo = od$pop)
  expect_equal(rn$coef$est, unname(coef(nb)[2]), tolerance = 1e-6)
  expect_equal(rn$fit$theta, nb$theta, tolerance = 1e-5)
  # no overdispersion to find: theta runs off, and the estimate is the
  # Poisson one
  rp <- suppressWarnings(
    didfe(pd$cnt, pd$town, pd$year, pd$on, fam = "nb", expo = pd$pop)
  )
  expect_gt(rp$fit$theta, 1000)
  expect_equal(rp$coef$est, r$coef$est, tolerance = 1e-4)
})

test_that("units and periods without counts are dropped", {
  y <- pd$cnt
  y[pd$town == 4] <- 0
  y[pd$year == 2016] <- 0
  r <- didfe(y, pd$town, pd$year, pd$on, fam = "pois")
  k <- pd$town != 4 & pd$year != 2016
  r2 <- didfe(y[k], pd$town[k], pd$year[k], pd$on[k], fam = "pois")
  expect_equal(r$coef, r2$coef, tolerance = testthat_tolerance())
  expect_identical(r$fit$dropped, 1)
  expect_identical(r$fit$n, 261)
})

test_that("the decomposition is bacondecomp's and sums to the estimate", {
  r <- didfe(pd$rate, pd$town, pd$year, pd$on, bacon = TRUE)
  b <- r$bacon
  expect_equal(sum(b$weight), 1, tolerance = 1e-14)
  expect_equal(sum(b$weight * b$est), r$coef$est, tolerance = 1e-12)
  # bacondecomp 0.1.1, dev/ref-pe.R, sorted by type, treated, control
  o <- order(b$type, b$treated, b$control)
  expect_equal(
    b$est[o],
    c(-2.3628249, -4.1439955, -2.8326409, 0.2160198, 1.2809375,
      -0.1383249, -5.4264955, -4.4146535, -2.5180267),
    tolerance = 1e-7
  )
  expect_equal(
    b$weight[o],
    c(0.03157895, 0.06315789, 0.05263158, 0.05263158, 0.06315789,
      0.03157895, 0.22105263, 0.26315789, 0.22105263),
    tolerance = 1e-7
  )
  expect_identical(unique(b$type[o]),
                   c("earlier vs later", "later vs earlier",
                     "treated vs never"))
  expect_true(all(is.infinite(b$control[b$type == "treated vs never"])))
  expect_output(print(r), "Goodman-Bacon")
})

test_that("leaving out each treated unit refits without it", {
  r <- didfe(pd$rate, pd$town, pd$year, pd$on, loo = TRUE)
  ever <- sort(unique(as.character(pd$town[pd$on == 1])))
  expect_identical(sort(r$loo$unit), ever)
  k <- pd$town != 1
  r1 <- didfe(pd$rate[k], pd$town[k], pd$year[k], pd$on[k])
  expect_equal(r$loo$est[r$loo$unit == "1"], r1$coef$est, tolerance = 1e-12)
  expect_equal(r$loo$se[r$loo$unit == "1"], r1$coef$se, tolerance = 1e-12)
  # one treated unit: without it nothing is treated
  one <- pd[pd$town %in% c(1, 4, 5, 9, 10), ]
  ro <- didfe(one$rate, one$town, one$year, one$on * (one$town == 1),
              loo = TRUE)
  expect_true(is.na(ro$loo$est))
})

test_that("matching takes the nearest comparison on the pre-period", {
  r <- didfe(pd$rate, pd$town, pd$year, pd$on, mtch = 2)
  # by hand: outcomes 2015-2017, pooled within-group covariance
  wide <- reshape(pd[pd$year < 2018, c("town", "year", "rate")],
                  idvar = "town", timevar = "year", direction = "wide")
  f <- as.matrix(wide[order(wide$town), -1])
  tr <- !is.na(tapply(pd$first, pd$town, max))
  s <- ((sum(tr) - 1) * cov(f[tr, ]) + (sum(!tr) - 1) * cov(f[!tr, ])) /
    (nrow(f) - 2)
  dd <- vapply(which(tr), function(i) {
    d <- sqrt(mahalanobis(f[!tr, ], f[i, ], s))
    which(!tr)[order(d)[1:2]]
  }, integer(2))
  want <- paste(rep(which(tr), each = 2), as.vector(dd))
  expect_setequal(paste(r$match$treated, r$match$control), want)
  w <- numeric(30)
  w[tr] <- 1
  for (j in as.vector(dd)) w[j] <- w[j] + 0.5
  pw <- pd
  pw$w <- w[pd$town]
  m <- lm(rate ~ on + factor(town) + factor(year), pw[pw$w > 0, ],
          weights = w)
  expect_equal(r$coef$est, unname(coef(m)[2]), tolerance = 1e-10)
  # within regions and a caliper: fewer matches
  ex <- rep(1:2, each = 150)
  re <- didfe(pd$rate, pd$town, pd$year, pd$on, mtch = 1, exact = ex,
              calip = 1)
  expect_true(all(re$match$dist <= max(r$match$dist)))
  expect_true(all((as.numeric(re$match$treated) > 15) ==
                    (as.numeric(re$match$control) > 15)))
})

test_that("the effect is recovered and its error matches the spread", {
  set.seed(11)
  est <- se <- numeric(300)
  for (k in seq_along(est)) {
    d <- expand.grid(t = 1:6, u = 1:24)
    tr <- d$u <= 12 & d$t >= 4
    y <- rnorm(24)[d$u] + 0.3 * d$t + 1.5 * tr +
      rnorm(24, 0, 0.5)[d$u] * (d$t - 3.5) / 3 + rnorm(144)
    r <- didfe(y, d$u, d$t, tr)
    est[k] <- r$coef$est
    se[k] <- r$coef$se
  }
  # the mean within 4 Monte Carlo errors of the effect; the mean error
  # within 10% of the spread (300 draws: the sd's own error is about 4%)
  expect_lt(abs(mean(est) - 1.5), 4 * sd(est) / sqrt(300))
  expect_lt(abs(mean(se) / sd(est) - 1), 0.1)
})

test_that("didfe() refuses what it cannot fit", {
  u <- pd$town
  t <- pd$year
  y <- pd$rate
  d <- pd$on
  expect_error(didfe(y, u, t, rep(0, 300)), "1 in some rows and 0")
  expect_error(didfe(y, u, t, d, expo = pd$pop), "exposure for counts")
  expect_error(didfe(-y, u, t, d, fam = "pois"), "at least 0")
  expect_error(didfe(y, u, t, d, bacon = TRUE, fam = "nb"), "fam = \"ols\"")
  expect_error(didfe(y, u, t, d, bacon = TRUE, w = pd$w), "plain two-way")
  expect_error(didfe(y, u, t, d, exact = u), "give `mtch` too")
  expect_error(didfe(y, u, t, d, calip = 1), "give `mtch` too")
  expect_error(didfe(y[-1], u[-1], t[-1], d[-1], bacon = TRUE), "balanced")
  expect_error(didfe(y, u, t, pd$on * (u > 1) + (u == 1) * (t == 2016),
                     bacon = TRUE), "stay 1")
  expect_error(didfe(y, u, t, as.integer(t >= 2015 & u <= 15), mtch = 1),
               "period before the first adoption")
  expect_error(didfe(y, u, t, as.integer(t >= 2018 & u <= 15) +
                       as.integer(t >= 2020 & u > 15), mtch = 1),
               "never treated")
  expect_error(didfe(y, u, t, d, mtch = 1, exact = t), "constant within")
  expect_error(didfe(y, u, t, d, mtch = 1, exact = u), "no treated unit")
  yy <- y
  yy[t == 2016] <- 2 * yy[t == 2015]
  yy[t == 2017] <- 3 * yy[t == 2015]
  expect_error(didfe(yy, u, t, d, mtch = 1), "collinear across units")
  expect_error(didfe(y, u, t, d, adj = data.frame(a = u)),
               "fixed effects absorb")
  expect_error(didfe(y, u, t, d, clus = rep(1, 300)), "one cluster")
  expect_error(didfe(y, u, t, d, clus = list(u, t, u)), "list of two")
  expect_error(didfe(y, u, t, d, clus = list(u, rep(1, 300))), "one cluster")
  expect_warning(didfe(y, u, t, d, clus = pd$region %% 5), "fewer than ten")
  t2 <- t
  t2[2] <- t2[1]
  expect_error(didfe(y, u, t2, d), "repeat a pair")
  expect_error(didfe(y[1:3], u[1:3], t[1:3], d[1:3]), "at least four")
  expect_error(didfe(y[1:10], u[1:10], t[1:10], d[1:10]), "two units")
  expect_error(didfe(y, u, t, d, adj = matrix(1, 2, 2)), "one row per row")
  expect_error(didfe(y, u, t, d, adj = matrix(NA_real_, 300, 1)),
               "finite numbers")
  expect_error(didfe(y, u, t, d, adj = list(a = 1:2)), "one value per row")
  expect_error(didfe(y, u, t, d, fe = "unit"), "fe")
  expect_error(didfe(y, u, t, d, level = 1), "level")
})

test_that("labelled controls enter as indicators; print shows the fit", {
  f <- ifelse(pd$x > 0, "hi", "lo")
  r <- didfe(pd$rate, pd$town, pd$year, pd$on,
             adj = list(f, z = pd$x, matrix(pd$x^2)))
  expect_identical(r$coef$term[1:2], c("trt", "adj1=lo"))
  r2 <- didfe(pd$rate, pd$town, pd$year, pd$on,
              adj = cbind(pd$x, b = 2 * pd$x^2))
  expect_identical(r2$coef$term, c("trt", "adj1", "b"))
  # treated from the first period: no comparison mean
  r4 <- didfe(pd$rate, pd$town, pd$year,
              as.integer(pd$town <= 5 | pd$on == 1))
  expect_true(is.na(r4$fit$ctlmean))
  r3 <- didfe(pd$rate, pd$town, pd$year, pd$on, adj = pd$x)
  expect_identical(r3$coef$term, c("trt", "adj"))
  expect_output(print(r3), "didfe: ols, 300 rows, 30 units")
})
