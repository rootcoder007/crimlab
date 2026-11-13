# jdgiv(), jmono() and mmshr(): against the covariance form of the Wald
# ratio, a population enumerated by type, chisq.test(), prop.test(),
# lm.wfit(), p.adjust(), and sums done here.

test_that("jdgiv() gives the Wald ratio and its defier sensitivity", {
  set.seed(7)
  z <- rbinom(500, 1, 0.5)
  d <- rbinom(500, 1, ifelse(z == 1, 0.55, 0.3))
  y <- rnorm(500, 1 + 0.8 * d)
  r <- jdgiv(z, d, y, dsh = c(0, 0.1), deff = c(-1, 1))
  cv <- function(a, b) mean(a * b) - mean(a) * mean(b)
  expect_equal(r$est[["wald"]], cv(z, y) / cv(z, d), tolerance = 1e-12)
  fs <- mean(d[z == 1]) - mean(d[z == 0])
  expect_equal(r$est[["fs"]], fs, tolerance = 1e-14)
  expect_equal(r$shares[["always"]], mean(d[z == 0]), tolerance = 1e-14)
  s <- r$sens
  expect_identical(s$dsh, c(0, 0.1, 0, 0.1))
  expect_equal(
    s$ceff,
    (r$est[["wald"]] * fs + s$dsh * s$deff) /
      (fs + s$dsh),
    tolerance = 1e-13
  )
  expect_identical(r$type, "obs")
  expect_output(print(r), "observed")
  # weights, and logical input
  w <- rep(c(2, 1), 250)
  rw <- jdgiv(z == 1, d == 1, y, w)
  expect_equal(
    rw$est[["itt"]],
    weighted.mean(y[z == 1], w[z == 1]) -
      weighted.mean(y[z == 0], w[z == 0]),
    tolerance = 1e-13
  )
  # no first stage: no ratio
  n0 <- jdgiv(c(1, 0, 1, 0), c(1, 1, 1, 1), 1:4)
  expect_true(is.na(n0$est[["wald"]]))
  expect_true(all(is.na(n0$sens$ceff)))
})

test_that("jdgiv() decomposes a known population exactly", {
  set.seed(8)
  n <- 400
  d0 <- rbinom(n, 1, 0.3)
  d1 <- ifelse(runif(n) < 0.9, pmax(d0, rbinom(n, 1, 0.5)), 1 - d0)
  y0 <- rnorm(n)
  y1 <- y0 + rnorm(n, 1)
  r <- jdgiv(d = cbind(d0, d1), y = data.frame(y0, y1))
  ty <- rep("never", n)
  ty[d0 == 1] <- "always"
  ty[d0 == 1 & d1 == 0] <- "defier"
  ty[d0 == 0 & d1 == 1] <- "complier"
  for (k in c("complier", "defier", "always", "never")) {
    expect_equal(r$shares[[k]], mean(ty == k), tolerance = 1e-14)
    expect_equal(r$effects[[k]], mean((y1 - y0)[ty == k]), tolerance = 1e-13)
  }
  sc <- r$shares
  ef <- r$effects
  # ITT and first stage as the type decomposition says
  expect_equal(
    r$est[["itt"]],
    sc[["complier"]] * ef[["complier"]] - sc[["defier"]] * ef[["defier"]],
    tolerance = 1e-13
  )
  expect_equal(
    r$est[["fs"]],
    sc[["complier"]] - sc[["defier"]],
    tolerance = 1e-14
  )
  expect_equal(r$est[["ate"]], mean(y1 - y0), tolerance = 1e-13)
  # no defiers: the ratio is the compliers' effect
  m <- jdgiv(
    d = cbind(c(0, 0, 1, 0), c(1, 1, 1, 0)),
    y = cbind(c(0, 0, 1, 2), c(3, 1, 1, 2))
  )
  expect_equal(m$est[["wald"]], m$effects[["complier"]], tolerance = 1e-15)
  lg <- jdgiv(
    d = cbind(c(FALSE, FALSE, TRUE, FALSE), c(TRUE, TRUE, TRUE, FALSE)),
    y = cbind(c(0L, 0L, 1L, 2L), c(3L, 1L, 1L, 2L))
  )
  expect_identical(lg$est, m$est)
  expect_true(is.na(m$effects[["defier"]]))
  expect_output(print(m), "population")
})

test_that("jdgiv() refuses bad input in words", {
  errs(
    "jdgiv",
    list(
      list("`z` must be 0 or 1", c(2, 0), c(1, 0), 1:2),
      list("`d` must have one value|`d` must be 0 or 1", c(1, 0), 1, 1:2),
      list("take both values", c(1, 1), c(1, 0), 1:2),
      list("`dsh` must be shares from 0 to 1", c(1, 0), c(1, 0), 1:2, dsh = 2),
      list("`deff` must be finite", c(1, 0), c(1, 0), 1:2, deff = NA),
      list("`w` must be at least 0", c(1, 0), c(1, 0), 1:2, w = c(-1, 1)),
      list("`d` must have two columns", NULL, 1:3, cbind(1, 2)),
      list("`y` must have two columns", NULL, cbind(0, 1), list(1)),
      list("`y` must hold numbers", NULL, cbind(0, 1), cbind("a", "b")),
      list("positive total", NULL, cbind(0, 1), cbind(1, 2), w = 0)
    )
  )
})

test_that("jmono() gives the rates, the homogeneity test and the slope test", {
  set.seed(9)
  j <- sample(paste0("J", 1:5), 600, TRUE)
  len <- c(J1 = 0.2, J2 = 0.3, J3 = 0.35, J4 = 0.4, J5 = 0.6)
  d <- rbinom(600, 1, len[j])
  y <- rbinom(600, 1, 0.3 + 0.2 * d)
  r <- jmono(j, d, y)
  rt <- c(tapply(d, j, mean))
  expect_identical(r$judges$judge, names(sort(rt)))
  expect_equal(r$judges$rate, unname(sort(rt)), tolerance = 1e-14)
  ym <- c(tapply(y, j, mean))
  expect_equal(r$judges$outcome, unname(ym[names(sort(rt))]), tolerance = 1e-14)
  ct <- suppressWarnings(chisq.test(table(j, d), correct = FALSE))
  expect_equal(r$homog[["chisq"]], unname(ct$statistic), tolerance = 1e-12)
  expect_equal(r$homog[["p"]], ct$p.value, tolerance = 1e-12)
  expect_identical(r$homog[["df"]], 4)
  # every pair, the more lenient first
  pr <- combn(5, 2)
  rr <- r$judges$rate
  oo <- r$judges$outcome
  s <- r$slope
  expect_equal(s$dp, rr[pr[2, ]] - rr[pr[1, ]], tolerance = 1e-14)
  expect_equal(s$dy, oo[pr[2, ]] - oo[pr[1, ]], tolerance = 1e-14)
  expect_identical(s$viol, abs(s$dy) > s$dp + 1e-12)
  expect_identical(s$j, r$judges$judge[pr[1, ]])
  b <- jmono(j, d, y, lo = -1, hi = 2)$slope
  expect_equal(b$bound, 3 * b$dp, tolerance = 1e-14)
  expect_output(print(r), "slope test")
  # weights move the rates, not the test
  w <- runif(600)
  rw <- jmono(j, d, w = w)
  expect_equal(
    rw$judges$rate[rw$judges$judge == "J1"],
    weighted.mean(d[j == "J1"], w[j == "J1"]),
    tolerance = 1e-14
  )
  expect_identical(rw$homog[["chisq"]], r$homog[["chisq"]])
  expect_null(rw$slope)
  one <- jmono(c("a", "b"), c(1, 1))
  expect_true(is.na(one$homog[["chisq"]]))
  expect_true(is.na(one$homog[["p"]]))
})

test_that("jmono() compares judges within groups", {
  set.seed(10)
  j <- sample(paste0("J", 1:4), 2000, TRUE)
  g <- sample(c("a", "b", "c"), 2000, TRUE)
  len <- c(J1 = 0.2, J2 = 0.35, J3 = 0.45, J4 = 0.6)
  d <- rbinom(2000, 1, len[j] + ifelse(g == "c" & j == "J1", 0.7, 0))
  r <- jmono(j, d, grp = g, minn = 100, adj = "BH")
  p <- r$sub$pairs
  jt <- c(tapply(d, j, mean))
  for (k in seq_len(nrow(p))) {
    s <- j == p$stricter[k] & g == p$group[k]
    l <- j == p$lenient[k] & g == p$group[k]
    expect_gt(jt[[p$stricter[k]]], jt[[p$lenient[k]]])
    pt <- prop.test(c(sum(d[l]), sum(d[s])), c(sum(l), sum(s)), correct = FALSE)
    expect_equal(p$z[k]^2, unname(pt$statistic), tolerance = 1e-10)
    expect_equal(p$grates[k], mean(d[s]), tolerance = 1e-14)
  }
  expect_equal(p$p, pnorm(-p$z), tolerance = 1e-15)
  expect_equal(p$padj, p.adjust(p$p, "BH"), tolerance = 1e-14)
  expect_identical(p$reversed, p$grates < p$gratel)
  sm <- r$sub$summary
  expect_identical(sm[["pairs"]], as.numeric(nrow(p)))
  expect_identical(sm[["signif"]], as.numeric(sum(p$signif)))
  expect_gt(sm[["refuted"]], 0)
  # the slope of each group's rates on the overall ones
  sl <- r$sub$slopes
  for (k in seq_len(nrow(sl))) {
    s <- g == sl$group[k]
    cr <- c(tapply(d[s], j[s], mean))
    cn <- c(tapply(d[s], j[s], length))
    f <- lm.wfit(cbind(1, jt[names(cr)]), cr, cn)
    expect_equal(sl$slope[k], unname(f$coefficients[2]), tolerance = 1e-12)
  }
  expect_output(print(r), "refuted")
  # no pair reaches minn: an empty table
  e <- jmono(j, d, grp = g, minn = 1e4)$sub
  expect_identical(nrow(e$pairs), 0L)
  expect_true(is.na(e$summary[["share"]]))
})

test_that("jmono() refuses bad input in words", {
  errs(
    "jmono",
    list(
      list("`d` must be 0 or 1", c("a", "b"), c(2, 0)),
      list("`judge` must have one label per value", "a", c(1, 0)),
      list("`minn` must be one number from 0", c("a", "b"), c(1, 0), minn = -1),
      list("`adj` must be one of", c("a", "b"), c(1, 0), adj = "x"),
      list("`alpha` must be one number", c("a", "b"), c(1, 0), alpha = 0),
      list("positive total weight", c("a", "b"), c(1, 0), w = c(0, 1)),
      list(
        "`lo` must be at most the smallest outcome",
        c("a", "b"),
        c(1, 0),
        y = 1:2,
        lo = 2
      ),
      list(
        "`hi` must be at least the largest outcome",
        c("a", "b"),
        c(1, 0),
        y = 1:2,
        hi = 1
      ),
      list("`grp` must have one label per value", c("a", "b"), c(1, 0), grp = 1)
    )
  )
})

test_that("jmono() refuses more pairs or cells than it can hold", {
  j <- seq_len(10001)
  d <- rep(0:1, length.out = 10001)
  errs("jmono", list(
    list("too many judges for every pair", j, d, y = d),
    list("more than 1e8 cells", j, d, grp = j)
  ))
})

test_that("mmshr() gives the shares against the minimum", {
  s <- c(365, 400, 300, 365, 730, 365)
  g <- c("a", "a", "b", "b", "b", "c")
  r <- mmshr(s, 365, g)
  expect_s3_class(r, "mmshr")
  expect_identical(r$group, c("a", "b", "c", "all"))
  e <- s - 365
  expect_equal(r$at[4], mean(e == 0), tolerance = 1e-15)
  expect_equal(r$below[2], 1 / 3, tolerance = 1e-15)
  expect_identical(r$mexcess, c(35, 365, 0, 200))
  expect_identical(r$n, c(2, 3, 1, 6))
  a <- mmshr(c(10, 12), c(10, 11))
  expect_identical(a$group, "all")
  expect_identical(a$above, 0.5)
  expect_output(print(r), "against the minimum")
  errs(
    "mmshr",
    list(
      list("`sent` must be finite", NA, 1),
      list("`min` must be one number or one per case, 2", 1:2, 1:3),
      list("`grp` must have one label per value", 1:2, 1, grp = 1)
    )
  )
})
