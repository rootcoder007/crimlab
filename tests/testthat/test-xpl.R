test_that("xpimp() measures the loss after each returned shuffle", {
  set.seed(41)
  n <- 60
  x <- cbind(a = rnorm(n), b = runif(n), c = rnorm(n))
  b <- c(1.5, -2, 0)
  f <- function(m) drop(m %*% b)
  y <- f(x) + rnorm(n, sd = 0.3)
  for (loss in c("se", "ae")) {
    for (ratio in c(FALSE, TRUE)) {
      r <- xpimp(x, f, y, nrep = 4, loss = loss, ratio = ratio, seed = 3)
      l <- function(e) if (loss == "se") mean(e^2) else mean(abs(e))
      base <- l(y - f(x))
      want <- matrix(0, 3, 4)
      for (j in 1:3) {
        for (k in 1:4) {
          z <- x
          z[, j] <- x[r$perm[, k, j], j]
          want[j, k] <- if (ratio) l(y - f(z)) / base else l(y - f(z)) - base
        }
      }
      expect_equal(r$base, base, tolerance = 1e-14)
      expect_equal(unname(r$rep), want, tolerance = 1e-13)
      expect_equal(r$imp$imp, rowMeans(want), tolerance = 1e-13)
      expect_equal(r$imp$sd, apply(want, 1L, sd), tolerance = 1e-12)
      expect_identical(r$imp$term, c("a", "b", "c"))
    }
  }
  # a column the function ignores has no importance
  expect_equal(r$rep["c", ], rep(1, 4), tolerance = 1e-14)
  # each shuffle is a permutation of the rows
  for (j in 1:3) {
    for (k in 1:4) expect_identical(sort(r$perm[, k, j]), seq_len(n))
  }
  expect_identical(xpimp(x, f, y, nrep = 4, seed = 3)$perm,
                   xpimp(x, f, y, nrep = 4, seed = 3)$perm)
  expect_false(identical(xpimp(x, f, y, nrep = 4, seed = 3)$perm,
                         xpimp(x, f, y, nrep = 4, seed = 4)$perm))
  expect_true(is.na(xpimp(x, f, y, nrep = 1)$imp$sd[1L]))
})

test_that("xpimp() fits the learners named in f once, by least squares", {
  set.seed(42)
  n <- 40
  x <- cbind(u = rnorm(n), v = rnorm(n))
  y <- 1 + 2 * x[, 1] - x[, 2] + rnorm(n)
  cf <- stats::lm.fit(cbind(1, x), y)$coefficients
  f <- function(m) drop(cbind(1, m) %*% cf)
  a <- xpimp(x, "glm", y, nrep = 3, seed = 9, lctl = list(nfold = 2))
  b <- xpimp(x, f, y, nrep = 3, seed = 9)
  expect_equal(a$rep, b$rep, tolerance = 1e-9)
  expect_equal(a$base, b$base, tolerance = 1e-9)
})

test_that("xpimp() refuses bad input in words", {
  x <- cbind(a = 1:5, b = c(2, 1, 4, 3, 5))
  f <- function(m) m[, 1]
  y <- 1:5
  e <- function(expr, m) expect_error(expr, m, fixed = TRUE)
  e(xpimp(1:5, f, y),
    "xpimp: `x` must be a numeric matrix with at least one row and one column.")
  e(xpimp(x + c(NA, 0, 0, 0, 0), f, y),
    "xpimp: `x` must be finite numbers with no missing values.")
  e(xpimp(x, 1, y),
    "xpimp: `f` must be a function or the names of learners of suplrn().")
  e(xpimp(x, f, y, lctl = list(ntree = 5)),
    "xpimp: `lctl` is for learners named in `f`, not for a function.")
  e(xpimp(x, "glm", y, lctl = list(trees = 5)),
    paste("xpimp: `lctl` must be a list named by settings of suplrn(): fam,",
          "nfold, fold, disc, alpha, nlam, nfin, ntree, mtry, seed."))
  e(xpimp(x, "glm", NULL),
    "xpimp: `y` must be given to fit the learners named in `f`.")
  e(xpimp(x, f, 1:4), "xpimp: `y` must have one value per row of `x`, 5.")
  e(xpimp(x, f, y, nrep = 0), "xpimp: `nrep` must be at least 1.")
  e(xpimp(x, f, y, loss = "l2"),
    "xpimp: `loss` must be one of \"se\", \"ae\", not \"l2\".")
  e(xpimp(x, f, y, ratio = NA), "xpimp: `ratio` must be TRUE or FALSE.")
  e(xpimp(x, f, y, seed = -1), "xpimp: `seed` must be at least 0.")
  e(xpimp(x, function(m) 1, y),
    "xpimp: `f` must return a finite number for each of the 105 rows.")
  e(xpimp(x, f, as.numeric(1:5), ratio = TRUE),
    "xpimp: `ratio` needs a loss above 0 before permuting.")
  e(xpimp(x, f, y, nrep = 4e6),
    paste("xpimp: `nrep` asks for 40000000 rows of 2 columns to predict, more",
          "than 67108864 values."))
})

test_that("xppdp() averages the function over the rows at each grid value", {
  set.seed(43)
  n <- 50
  x <- cbind(a = rnorm(n), b = runif(n))
  f <- function(m) m[, "a"]^2 + 2 * m[, "b"] + m[, "a"] * m[, "b"]
  pd <- function(g) {
    vapply(g, function(v) mean(v^2 + 2 * x[, "b"] + v * x[, "b"]), 0)
  }
  r <- xppdp(x, f, "a", grid = 7)
  g <- seq(min(x[, "a"]), max(x[, "a"]), length.out = 7)
  expect_equal(r$pd$grid, g, tolerance = 1e-15)
  expect_equal(r$pd$pd, pd(g), tolerance = 1e-13)
  expect_identical(r$term, "a")
  expect_identical(r$n, as.integer(n))
  q <- xppdp(x, f, 1, grid = 6, gtype = "quant")
  g <- unique(stats::quantile(x[, "a"], seq(0, 1, length.out = 6),
                              names = FALSE))
  expect_equal(q$pd$grid, g, tolerance = 1e-15)
  expect_equal(q$pd$pd, pd(g), tolerance = 1e-13)
  u <- xppdp(x, f, "a", grid = c(2, -1, 0))
  expect_identical(u$pd$grid, c(2, -1, 0))
  expect_equal(u$pd$pd, pd(c(2, -1, 0)), tolerance = 1e-13)
  # repeated quantiles are dropped
  z <- cbind(a = c(0, 0, 0, 0, 1), b = 1:5)
  expect_identical(xppdp(z, f, "a", grid = 5, gtype = "quant")$pd$grid,
                   c(0, 1))
})

test_that("xppdp() of a least-squares learner is its slope", {
  set.seed(44)
  x <- cbind(u = rnorm(30), v = rnorm(30))
  y <- 1 - 3 * x[, 1] + x[, 2] + rnorm(30)
  cf <- stats::lm.fit(cbind(1, x), y)$coefficients
  r <- xppdp(x, "glm", "u", y = y, grid = 4, lctl = list(nfold = 2))
  g <- r$pd$grid
  want <- cf[[1L]] + cf[[2L]] * g + cf[[3L]] * mean(x[, 2])
  expect_equal(r$pd$pd, want, tolerance = 1e-9)
})

test_that("xppdp() refuses bad input in words", {
  x <- cbind(a = c(1, 2, 3), b = c(3, 1, 2))
  f <- function(m) m[, 1]
  e <- function(expr, m) expect_error(expr, m, fixed = TRUE)
  e(xppdp(x, f, "z"), "xppdp: `v` must name or number one column of `x`.")
  e(xppdp(x, f, 3), "xppdp: `v` must name or number one column of `x`.")
  e(xppdp(x, f, "a", grid = 1),
    paste("xppdp: `grid` must be a whole number of points from 2, or at",
          "least two values."))
  e(xppdp(x, f, "a", grid = c(1, NA)),
    paste("xppdp: `grid` must be a whole number of points from 2, or at",
          "least two values."))
  e(xppdp(x, f, "a", gtype = "odd"),
    "xppdp: `gtype` must be one of \"even\", \"quant\", not \"odd\".")
})

test_that("xpcp() predicts each row with one input moved over the grid", {
  set.seed(45)
  x <- cbind(a = rnorm(40), b = runif(40))
  f <- function(m) exp(m[, "a"]) * m[, "b"]
  r <- xpcp(x, f, "a", newx = x[2:4, ], grid = 5)
  g <- seq(min(x[, "a"]), max(x[, "a"]), length.out = 5)
  expect_equal(r$grid, g, tolerance = 1e-15)
  expect_equal(r$prof, outer(x[2:4, "b"], exp(g)), tolerance = 1e-14)
  expect_equal(r$pred, f(x[2:4, ]), tolerance = 1e-15)
  expect_identical(r$term, "a")
  # the default row is the first; the mean profile over x is xppdp()
  expect_equal(xpcp(x, f, 2, grid = c(0, 1))$prof,
               matrix(exp(x[1L, "a"]) * c(0, 1), 1L), tolerance = 1e-15)
  a <- xpcp(x, f, "a", newx = x, grid = 4)
  expect_equal(colMeans(a$prof), xppdp(x, f, "a", grid = 4)$pd$pd,
               tolerance = 1e-14)
})

test_that("xpcp() refuses bad input in words", {
  x <- cbind(a = c(1, 2, 3), b = c(3, 1, 2))
  f <- function(m) m[, 1]
  m <- paste("xpcp: `newx` must be a numeric matrix with the 2 columns of",
             "`x`, no missing values.")
  expect_error(xpcp(x, f, "a", newx = x[, 1, drop = FALSE]), m, fixed = TRUE)
  expect_error(xpcp(x, f, "a", newx = c(1, 2)), m, fixed = TRUE)
  expect_error(xpcp(x, f, "a", newx = cbind(1, NA)), m, fixed = TRUE)
  expect_error(xpcp(x, f, "c"),
               "xpcp: `v` must name or number one column of `x`.",
               fixed = TRUE)
})

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
    expect_identical(r$bins$n, w$n)
    expect_equal(r$bins$eff, w$eff, tolerance = 1e-12)
    expect_identical(r$bins$lo, w$z[-length(w$z)])
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
  expect_identical(r$ale$z, w$z)
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

test_that("xpshap() of a linear function against one reference row", {
  x <- cbind(a = 1, b = 2, c = -1)
  b <- c(2, -1, 0.5)
  f <- function(m) drop(m %*% b)
  newx <- rbind(c(3, 0, 1), c(-1, 4, 2))
  r <- xpshap(x, f, newx = newx, nsamp = 7)
  want <- t(t(sweep(newx, 2L, x[1L, ])) * b)
  expect_equal(unname(r$phi), want, tolerance = 1e-14)
  expect_equal(unname(r$se), matrix(0, 2, 3), tolerance = 1e-14)
  expect_identical(colnames(r$phi), c("a", "b", "c"))
  expect_equal(r$pred, f(newx), tolerance = 1e-15)
  expect_equal(r$base, f(x), tolerance = 1e-15)
})

test_that("xpshap() samples the Shapley values of a game with interactions", {
  # exact values by enumerating the coalitions: v(S) = f(newx on S, w off)
  x <- cbind(a = 0.5, b = -1, c = 2)
  newx <- cbind(a = 2, b = 1, c = -1)
  f <- function(m) m[, 1] * m[, 2] + exp(m[, 3] / 2) * m[, 1] + m[, 2]^2
  v <- function(s) {
    z <- x
    z[1L, s] <- newx[1L, s]
    f(z)
  }
  exact <- vapply(1:3, function(j) {
    o <- setdiff(1:3, j)
    ss <- list(integer(0), o[1L], o[2L], o)
    sum(vapply(ss, function(s) {
      k <- length(s)
      factorial(k) * factorial(2 - k) / 6 * (v(c(s, j)) - v(s))
    }, 0))
  }, 0)
  r <- xpshap(x, f, newx = newx, nsamp = 4000, seed = 2)
  expect_true(all(abs(r$phi[1L, ] - exact) <= 5 * r$se[1L, ] + 1e-12))
  expect_true(all(r$se > 0))
  # the exact values add up to the prediction less the reference
  expect_equal(sum(exact), unname(f(newx) - f(x)), tolerance = 1e-13)
  expect_identical(xpshap(x, f, newx = newx, nsamp = 5, seed = 2)$phi,
                   xpshap(x, f, newx = newx, nsamp = 5, seed = 2)$phi)
})

test_that("xpshap() refuses bad input in words", {
  x <- cbind(a = c(1, 2), b = c(3, 1))
  f <- function(m) m[, 1]
  expect_error(xpshap(x, f, nsamp = 0), "xpshap: `nsamp` must be at least 1.",
               fixed = TRUE)
  expect_error(xpshap(x, f, newx = 1),
               paste("xpshap: `newx` must be a numeric matrix with the 2",
                     "columns of `x`, no missing values."), fixed = TRUE)
  expect_error(xpshap(x, f, seed = NA),
               "xpshap: `seed` must be finite, with no missing values.",
               fixed = TRUE)
})
