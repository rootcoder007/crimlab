test_that("ulcate() from given nuisances is least squares on the clipped U", {
  set.seed(51)
  n <- 120
  x <- matrix(rnorm(n * 2), n)
  d <- rbinom(n, 1, 0.5)
  y <- x[, 1] + d * (1 + x[, 2]) + rnorm(n)
  m <- y + rnorm(n, sd = 0.5)
  e <- c(0.995, 0.002, runif(n - 2, 0.1, 0.9))
  d[1:2] <- c(1, 0)
  r <- ulcate(y, d, x, fold = rep(1:2, length.out = n), lrnt = "glm",
              pred = cbind(m, e), clip = 0.01)
  den <- d - e
  lo <- abs(den) < 0.01
  den[lo] <- ifelse(den[lo] >= 0, 0.01, -0.01)
  u <- (y - m) / den
  expect_equal(r$u, u, tolerance = 1e-14)
  expect_equal(den[1:2], c(0.01, -0.01))
  tau <- drop(cbind(1, x) %*% stats::lm.fit(cbind(1, x), u)$coefficients)
  expect_equal(r$tau, tau, tolerance = 1e-10)
  expect_equal(r$ate, mean(tau), tolerance = 1e-12)
  expect_identical(unname(r$pred), unname(cbind(m, e)))
  expect_identical(colnames(r$pred), c("m", "e"))
  expect_s3_class(r, "ulcate")
})

test_that("ulcate() learns m and e out of fold", {
  set.seed(52)
  n <- 150
  x <- matrix(rnorm(n * 2), n)
  d <- rbinom(n, 1, plogis(0.4 * x[, 1]))
  y <- x[, 2] + d * (1 + x[, 1]) + rnorm(n)
  fold <- sample(rep(1:3, length.out = n))
  r <- ulcate(y, d, x, fold, lrn = "glm", lrnd = "glm", lrnt = "glm")
  m <- e <- numeric(n)
  for (k in 1:3) {
    tr <- fold != k
    xt <- cbind(1, x[tr, ])
    b <- stats::lm.fit(xt, y[tr])$coefficients
    g <- stats::glm.fit(xt, d[tr], family = stats::binomial(),
                        control = list(epsilon = 1e-14, maxit = 100))
    m[!tr] <- drop(cbind(1, x[!tr, ]) %*% b)
    e[!tr] <- stats::plogis(drop(cbind(1, x[!tr, ]) %*% g$coefficients))
  }
  expect_equal(unname(r$pred[, "m"]), m, tolerance = 1e-10)
  expect_equal(unname(r$pred[, "e"]), e, tolerance = 1e-8)
  u <- (y - m) / (d - e)
  expect_equal(r$u, u, tolerance = 1e-7)
  # tau is the least-squares fit of u
  b <- stats::lm.fit(cbind(1, x), r$u)$coefficients
  expect_equal(r$tau, drop(cbind(1, x) %*% b), tolerance = 1e-9)
})

test_that("ulcate() refuses bad input in words", {
  x <- matrix(rnorm(20), 10)
  y <- rnorm(10)
  d <- rep(0:1, 5)
  f <- rep(1:2, 5)
  e <- function(expr, m) expect_error(expr, m, fixed = TRUE)
  e(ulcate(y, rep(2, 10), x, f),
    "ulcate: `d` must be 0 or 1 for a binary outcome.")
  e(ulcate(y, rep(1, 10), x, f), "ulcate: `d` is constant.")
  e(ulcate(y, d, x, f, clip = 0.5), "ulcate: `clip` must be below 0.5.")
  e(ulcate(y, d, x, f, pred = matrix(0.5, 10, 3)),
    "ulcate: `pred` must be a matrix of 10 rows and 2 columns.")
  e(ulcate(y, d, x, f, pred = matrix(c(NA, rep(0.5, 19)), 10)),
    "ulcate: `pred` must be finite.")
  e(ulcate(y, d, x[1:3, , drop = FALSE], f[1:3]),
    "ulcate: `x` must have at least four rows and one column.")
})
