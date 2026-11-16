test_that("lrcce() is the mean group of the units' long-run effects", {
  # unit estimates written out: theta_i = (b0 + b1) / (1 - l1 - l2)
  b <- cbind(L1.y = c(0.3, 0.4, 0.5, 0.2), L2.y = c(0.1, 0, 0.1, 0.2),
             x = c(0.5, 0.6, 0.4, 0.5), xl = c(0.1, 0, 0.2, 0.1),
             z = c(1, 2, 3, 4))
  r <- lrcce(b, list(x = c("x", "xl"), z = "z"))
  th <- (b[, "x"] + b[, "xl"]) / (1 - b[, "L1.y"] - b[, "L2.y"])
  tz <- b[, "z"] / (1 - b[, "L1.y"] - b[, "L2.y"])
  ph <- 1 - b[, "L1.y"] - b[, "L2.y"]
  expect_identical(r$term, c("x", "z", "phi"))
  expect_equal(r$est, c(mean(th), mean(tz), mean(ph)), tolerance = 1e-14)
  expect_equal(r$se, c(sd(th), sd(tz), sd(ph)) / 2, tolerance = 1e-14)
  expect_equal(r$lo, r$est - stats::qnorm(0.975) * r$se, tolerance = 1e-14)
  expect_equal(unname(attr(r, "units")[, 1L]), unname(th), tolerance = 1e-14)
  # a data frame and a ccepnl fit's bi give the same
  expect_identical(lrcce(as.data.frame(b), list(x = c("x", "xl"), z = "z")),
                   r)
})

test_that("lrcce() reads a ccepnl fit and refuses what it cannot use", {
  set.seed(5)
  d <- expand.grid(year = 1:20, prov = paste0("p", 1:12))
  d$x <- rnorm(nrow(d))
  d$y <- 0
  for (i in which(d$year > 1)) d$y[i] <- 0.3 * d$y[i - 1] + d$x[i] + rnorm(1)
  k <- d$year > 1
  f <- ccepnl(d$y[k], d[k, "x", drop = FALSE], d$prov[k], d$year[k],
              ylag = 1)
  r <- lrcce(f, "x")
  expect_equal(r$est[1L], mean(f$bi[, "x"] / (1 - f$bi[, "L1.y"])),
               tolerance = 1e-14)
  expect_error(lrcce(f, "w"),
               "lrcce: `x` and `ylag` must name columns of the unit estimates.",
               fixed = TRUE)
  expect_error(lrcce(f$bi[, "x", drop = FALSE], "x"),
               "lrcce: the fit has no lags of `y`.", fixed = TRUE)
  expect_error(lrcce(1:3, "x"),
               paste("lrcce: `object` must be a ccepnl fit or a matrix of",
                     "unit estimates with named columns."), fixed = TRUE)
  expect_error(lrcce(cbind(L1.y = c(1, 0.5), x = c(1, 1)), "x"),
               "lrcce: a unit's lags of `y` sum to 1.", fixed = TRUE)
})
