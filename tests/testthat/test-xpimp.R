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
