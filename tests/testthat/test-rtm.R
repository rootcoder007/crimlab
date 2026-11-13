# crchyp(), hitscr() and rtmfit(): against distances, counts and Poisson
# fits done here in base R, and the refusals.

test_that("crchyp() draws the circle on the farthest pair", {
  set.seed(3)
  x <- runif(9, 0, 10)
  y <- runif(9, 0, 10)
  d <- as.matrix(dist(cbind(x, y)))
  d[lower.tri(d, diag = TRUE)] <- -1
  ij <- which(d == max(d), arr.ind = TRUE)[1L, ]
  ij <- sort(unname(ij))
  cc <- c(x = mean(x[ij]), y = mean(y[ij]))
  r <- max(d) / 2
  home <- rbind(cc + c(r / 2, 0), cc + c(0, 2 * r))
  out <- crchyp(x, y, home)
  expect_named(out, c("centre", "r", "pair", "share", "mar"))
  expect_identical(out$pair, ij)
  expect_equal(out$centre, cc, tolerance = 1e-15)
  expect_equal(out$r, r, tolerance = 1e-14)
  inside <- sqrt((x - cc[1])^2 + (y - cc[2])^2) <= r + 1e-12 * max(1, r)
  expect_identical(out$share, mean(inside))
  expect_identical(out$mar, c(TRUE, FALSE))
  # one home as c(x, y); no home
  expect_identical(crchyp(c(0, 4, 2), c(0, 0, 1), c(2, -1))$mar, TRUE)
  expect_null(crchyp(c(0, 4), c(0, 0))$mar)
  expect_identical(crchyp(c(0, 4), c(0, 0))$share, 1)
})

test_that("hitscr() counts the cells searched before the home", {
  set.seed(4)
  s <- round(runif(40), 1)
  h <- c(3, 17, 40)
  out <- hitscr(s, h)
  ab <- vapply(h, function(k) sum(s > s[k]), 1L)
  ti <- vapply(h, function(k) sum(s == s[k]), 1L)
  expect_identical(out$home, as.integer(h))
  expect_identical(out$above, ab)
  expect_identical(out$ties, ti)
  expect_identical(out$rank, ab + 1L)
  expect_identical(out$hit, 100 * ab / 40)
  expect_identical(out$hitt, 100 * (ab + ti) / 40)
})

test_that("rtmfit() is the Poisson regression of the counts on the layers", {
  set.seed(5)
  n <- 60
  l1 <- rbinom(n, 1, 0.4)
  l2 <- rpois(n, 2)
  y <- rpois(n, exp(0.2 + 0.7 * l1 + 0.25 * l2))
  g <- glm(y ~ l1 + l2, family = poisson(),
           control = glm.control(epsilon = 1e-14, maxit = 100))
  b <- unname(coef(g))
  # the information at the estimate (vcov() of glm() takes the weights of
  # the step before the last)
  xm <- cbind(1, l1, l2)
  se <- sqrt(diag(solve(crossprod(xm, xm * fitted(g)))))
  se <- unname(se)
  out <- rtmfit(y, data.frame(bars = l1, stops = l2), level = 0.9)
  expect_identical(out$coef$term, c("(Intercept)", "bars", "stops"))
  expect_equal(out$coef$est, b, tolerance = 1e-9)
  expect_equal(out$coef$se, se, tolerance = 1e-9)
  expect_equal(out$coef$z, b / se, tolerance = 1e-9)
  expect_equal(out$coef$p, 2 * pnorm(-abs(b / se)), tolerance = 1e-8)
  expect_equal(out$coef$rrv, exp(b), tolerance = 1e-9)
  expect_equal(out$coef$lo, exp(b - qnorm(0.95) * se), tolerance = 1e-9)
  expect_equal(out$coef$hi, exp(b + qnorm(0.95) * se), tolerance = 1e-9)
  rr <- exp(b[2] * l1 + b[3] * l2)
  expect_equal(out$cells$rr, rr, tolerance = 1e-9)
  expect_equal(out$cells$score, rr / min(rr), tolerance = 1e-9)
  expect_equal(out$cells$fitted, unname(fitted(g)), tolerance = 1e-9)
  expect_identical(out$cells$y, as.numeric(y))
  expect_equal(out$dev, deviance(g), tolerance = 1e-9)
  expect_equal(out$aic, AIC(g), tolerance = 1e-9)
  expect_true(out$conv)
  expect_identical(out$level, 0.9)
  expect_identical(dimnames(out$lay)[[2]], c("bars", "stops"))
  # a matrix without names, a plain vector, an aliased layer
  m <- rtmfit(y, cbind(l1, l2))
  expect_identical(m$coef$term, c("(Intercept)", "l1", "l2"))
  v <- rtmfit(y, l1)
  expect_identical(v$coef$term, c("(Intercept)", "layer1"))
  expect_equal(v$coef$est,
               unname(coef(glm(y ~ l1, family = poisson(),
                               control = glm.control(1e-14, 100)))),
               tolerance = 1e-9)
  a <- rtmfit(y, list(l1, 2 * l1))
  expect_identical(a$coef$term, c("(Intercept)", "layer1", "layer2"))
  expect_true(is.na(a$coef$est[3]))
  expect_equal(a$cells$rr, exp(v$coef$est[2] * l1), tolerance = 1e-12)
})

test_that("rtmfit() builds proximity, density and kernel layers", {
  cell <- as.matrix(expand.grid(x = 1:5, y = 1:5))
  f1 <- rbind(c(1, 1), c(4, 4), c(4.5, 4))
  f2 <- rbind(c(2, 5))
  y <- c(3, 1, 0, 0, 0, 2, 1, 0, 1, 0, 0, 1, 2, 2, 0, 1, 0, 2, 4, 3, 1, 0,
         1, 2, 1)
  dd <- function(f) {
    sqrt(outer(cell[, 1], f[, 1], "-")^2 + outer(cell[, 2], f[, 2], "-")^2)
  }
  dens <- function(f) rowSums(dd(f) <= 1.5)
  fit <- function(l) {
    unname(coef(glm(y ~ ., data = as.data.frame(l), family = poisson(),
                    control = glm.control(1e-14, 100))))
  }
  out <- rtmfit(y, list(a = f1, b = f2), cell = cell, rad = 1.5,
                kind = "dens")
  expect_identical(unname(out$lay[, "a"]), as.numeric(dens(f1)))
  expect_identical(unname(out$lay[, "b"]), as.numeric(dens(f2)))
  expect_equal(out$coef$est, fit(cbind(dens(f1), dens(f2))),
               tolerance = 1e-9)
  p <- rtmfit(y, list(a = f1), cell = cell, rad = 1.5)
  expect_identical(unname(p$lay[, 1]), as.numeric(dens(f1) > 0))
  k <- rtmfit(y, f1, cell = cell, rad = 1, kind = "kde")
  kd <- rowSums(exp(-0.5 * dd(f1)^2))
  expect_equal(unname(k$lay[, 1]), as.numeric(scale(kd)), tolerance = 1e-13)
  expect_identical(k$coef$term, c("(Intercept)", "layer1"))
  expect_equal(k$coef$est, fit(scale(kd)), tolerance = 1e-9)
})

test_that("crchyp(), hitscr() and rtmfit() refuse bad input", {
  errs("crchyp", list(
    list("at least two sites", 1, 1),
    list("`y` must have length", c(1, 2), 1),
    list("`home` must be a numeric matrix of 2 columns", c(1, 2), c(1, 2),
         c(1, 2, 3)),
    list("`x` must be finite", c(1, NA), c(1, 2))
  ))
  errs("hitscr", list(
    list("`home` must be cell numbers from 1 to 3", c(1, 2, 3), 4),
    list("`home` must be cell numbers", c(1, 2, 3), 1.5),
    list("`s` must be finite", c(1, Inf), 1)
  ))
  y <- c(1, 0, 2)
  errs("rtmfit", list(
    list("at least one count above 0", c(0, 0), c(1, 2)),
    list("`level` must be one number", y, 1:3, level = 1),
    list("`lay` must hold at least one layer", y, list()),
    list("`lay` must have length 3", y, list(1:2)),
    list("`lay` must have one row per cell, 3", y, matrix(1, 2, 1)),
    list("`lay` must be a numeric matrix", y, "a"),
    list("`cell` must have one row per count, 3", y, list(cbind(1, 1)),
         cell = cbind(1:2, 1:2), rad = 1),
    list("`rad` must be one number", y, list(cbind(1, 1)),
         cell = cbind(1:3, 1:3), rad = 0),
    list("`kind` must be one of", y, list(cbind(1, 1)),
         cell = cbind(1:3, 1:3), rad = 1, kind = "x"),
    list("`lay` must hold at least one layer", y, list(),
         cell = cbind(1:3, 1:3), rad = 1),
    list("`lay` must be a numeric matrix of 2 columns", y, list(1:3),
         cell = cbind(1:3, 1:3), rad = 1),
    list("needs at least two cells", 1, list(cbind(1, 1)),
         cell = cbind(1, 1), rad = 1, kind = "kde")
  ))
})
