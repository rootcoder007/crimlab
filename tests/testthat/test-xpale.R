# The algorithm of 'ALEPlot' (Apley and Zhu 2020) written out with
# quantile(), cut() and tapply().
aleref <- function(f, x, j, K) {
  z <- c(min(x[, j]), stats::quantile(x[, j], seq(1 / K, 1, length.out = K),
                                      type = 1, names = FALSE))
  z <- unique(z)
  K <- length(z) - 1L
  b <- as.numeric(cut(x[, j], breaks = z, include.lowest = TRUE))
  x1 <- x
  x2 <- x
  x1[, j] <- z[b]
  x2[, j] <- z[b + 1L]
  d <- as.numeric(tapply(f(x2) - f(x1), b, mean))
  a <- c(0, cumsum(d))
  w <- as.numeric(table(b))
  list(z = z, ale = a - sum((a[1:K] + a[2:(K + 1)]) / 2 * w) / sum(w),
       n = w, eff = d)
}

test_that("xpale() accumulates the local effects as ALEPlot does", {
  set.seed(46)
  x <- cbind(a = rnorm(80), b = rnorm(80))
  x[, "b"] <- x[, "b"] + x[, "a"]
  f <- function(m) m[, "a"] * m[, "b"] + m[, "a"]^2
  for (K in c(1, 4, 10)) {
    r <- xpale(x, f, "a", nbin = K)
    w <- aleref(f, x, 1, K)
    expect_equal(r$ale$z, w$z, tolerance = 1e-15)
    expect_equal(r$ale$ale, w$ale, tolerance = 1e-12)
    expect_equal(r$bins$n, w$n)
    expect_equal(r$bins$eff, w$eff, tolerance = 1e-12)
    expect_equal(r$bins$lo, w$z[-length(w$z)])
  }
  # a linear term: the effect is the slope times the centred edges
  g <- function(m) 3 * m[, "a"] + sin(m[, "b"])
  r <- xpale(x, g, "a", nbin = 5)
  z <- r$ale$z
  a <- 3 * (z - z[1L])
  expect_equal(r$ale$ale,
               a - sum((a[-length(a)] + a[-1L]) / 2 * r$bins$n) / 80,
               tolerance = 1e-12)
  # tied values: repeated edges dropped
  t <- cbind(a = rep(1:4, c(10, 1, 1, 8)), b = 1:20)
  w <- aleref(f, t, 1, 10)
  r <- xpale(t, f, 1, nbin = 10)
  expect_equal(r$ale$z, w$z)
  expect_equal(r$ale$ale, w$ale, tolerance = 1e-12)
})

test_that("xpale() refuses bad input in words", {
  x <- cbind(a = c(1, 1, 1), b = c(3, 1, 2))
  f <- function(m) m[, 1]
  expect_error(xpale(x, f, "a"),
               paste("xpale: `v` names a column with one value: it has no",
                     "local effects."), fixed = TRUE)
  expect_error(xpale(x, f, "b", nbin = 0), "xpale: `nbin` must be at least 1.",
               fixed = TRUE)
  expect_error(xpale(x, f, "b", nbin = 1.5),
               "xpale: `nbin` must be a whole number.", fixed = TRUE)
})
