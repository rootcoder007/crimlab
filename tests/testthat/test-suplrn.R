# Non-negative least squares from scratch: every subset of the columns of
# z, its least squares fit kept when every coefficient is positive, and the
# feasible fit of least residual sum of squares chosen (the solution of
# NNLS is the least squares fit on its own support).
nnref <- function(z, y) {
  k <- ncol(z)
  best <- rep(0, k)
  rss <- sum(y^2)
  for (m in seq_len(2^k - 1)) {
    s <- which(bitwAnd(m, 2^(seq_len(k) - 1)) > 0)
    b <- qr.solve(z[, s, drop = FALSE], y)
    if (all(b > 0)) {
      r <- sum((y - z[, s, drop = FALSE] %*% b)^2)
      if (r < rss) {
        rss <- r
        best <- rep(0, k)
        best[s] <- b
      }
    }
  }
  best
}

dat <- function(n = 200, bin = FALSE, seed = 1) {
  set.seed(seed)
  x <- matrix(rnorm(n * 4), n)
  e <- sin(2 * x[, 1]) + x[, 2] - 0.5 * x[, 3]
  y <- if (bin) rbinom(n, 1, plogis(e)) else e + rnorm(n, sd = 0.5)
  list(x = x, y = y)
}

test_that("the weights are the scaled NNLS fit on the cross-fitted columns", {
  for (bin in c(FALSE, TRUE)) {
    d <- dat(bin = bin)
    s <- suplrn(d$x, d$y, if (bin) "binomial" else "gaussian",
                lrns = c("glm", "ridge", "forest"), ntree = 30, nlam = 20)
    expect_true(all(s$wt >= 0), info = bin)
    expect_equal(sum(s$wt), 1, tolerance = 1e-15)
    w <- nnref(s$oof, d$y)
    expect_equal(unname(s$wt), w / sum(w), tolerance = 1e-9, info = bin)
    expect_equal(s$risk, colMeans((d$y - s$oof)^2), tolerance = 1e-14)
    expect_named(s$wt, c("glm", "ridge", "forest"))
    expect_identical(colnames(s$oof), c("glm", "ridge", "forest"))
    if (bin) expect_true(all(s$pred >= 0 & s$pred <= 1))
  }
})

test_that("the weights stay NNLS when a learner enters and another leaves", {
  # these data take the active set method through a step back (a learner's
  # weight reaching 0 as another enters)
  set.seed(1)
  x <- matrix(rnorm(360), 120)
  y <- x[, 1] + 0.5 * x[, 2]^2 + rnorm(120)
  s <- suplrn(x, y, lrns = c("glm", "lasso", "ridge", "forest"), nlam = 10,
              ntree = 20, nfold = 3, seed = 1)
  w <- nnref(s$oof, y)
  expect_equal(unname(s$wt), w / sum(w), tolerance = 1e-9)
})

test_that("each fold's predictions come from fits on the other folds", {
  d <- dat()
  s <- suplrn(d$x, d$y, lrns = "glm", nfold = 4)
  for (k in 1:4) {
    i <- s$fold == k
    f <- lm.fit(cbind(1, d$x[!i, ]), d$y[!i])
    expect_equal(s$oof[i, 1], drop(cbind(1, d$x[i, ]) %*% f$coefficients),
                 tolerance = 1e-12, info = k)
  }
  expect_identical(tabulate(s$fold), rep(50L, 4))
})

test_that("with one learner the stack is that learner", {
  d <- dat()
  nx <- matrix(rnorm(40), 10)
  s <- suplrn(d$x, d$y, lrns = "forest", ntree = 20, newx = nx, seed = 3)
  expect_identical(unname(s$wt), 1)
  expect_identical(s$pred, rfor(d$x, d$y, nx, ntree = 20, seed = 3)$pred)
  b <- dat(bin = TRUE)
  s <- suplrn(b$x, b$y, "binomial", lrns = "glm", newx = nx)
  expect_equal(s$pred, drop(plogis(cbind(1, nx) %*% lgtreg(b$x, b$y)$coef)),
               tolerance = 1e-12)
})

test_that("the discrete stack picks the learner of least risk", {
  d <- dat()
  s <- suplrn(d$x, d$y, lrns = c("glm", "forest"), ntree = 30, disc = TRUE)
  expect_identical(unname(s$wt), as.numeric(seq_along(s$risk) ==
                                              which.min(s$risk)))
  # noise: no learner beats the mean, every NNLS weight is 0, and the
  # stack falls back on the learner of least risk
  set.seed(4)
  x <- matrix(rnorm(60), 30)
  y <- rnorm(30)
  y <- y - mean(y)
  s <- suplrn(x, y, lrns = c("glm", "ridge"), nlam = 10, nfold = 3)
  expect_identical(unname(s$wt), as.numeric(seq_along(s$risk) ==
                                              which.min(s$risk)))
})

test_that("the folds and the result depend on the seed alone", {
  d <- dat()
  a <- suplrn(d$x, d$y, lrns = c("glm", "lasso"), nlam = 10, seed = 5)
  expect_identical(suplrn(d$x, d$y, lrns = c("glm", "lasso"), nlam = 10,
                          seed = 5), a)
  expect_identical(suplrn(d$x, d$y, lrns = c("glm", "lasso"), nlam = 10,
                          fold = a$fold, seed = 5), a)
  expect_true(all(table(a$fold) %in% 40L))
  expect_false(identical(suplrn(d$x, d$y, lrns = "glm", seed = 6)$fold,
                         a$fold))
})

test_that("suplrn refuses bad input", {
  d <- dat(n = 20)
  errs("suplrn", list(
    list("at least four rows", matrix(1:3, 3), 1:3),
    list("`fam` must be one of", d$x, d$y, fam = "pois"),
    list("names \"glm\" twice", d$x, d$y, lrns = c("glm", "glm")),
    list("from one to five", d$x, d$y, lrns = character()),
    list("`lrns` must be one of", d$x, d$y, lrns = "stack"),
    list("`nfold` must be at most", d$x, d$y, nfold = 21),
    list("`newx` must have the 4 columns", d$x, d$y, newx = matrix(1, 2, 3)),
    list("a learner failed on a fold", d$x[, c(1, 1)], d$y, lrns = "glm")
  ))
})

test_that("dmlx takes the stack as a nuisance learner", {
  set.seed(1)
  n <- 300
  x <- matrix(rnorm(n * 3), n)
  d <- x[, 1] + rnorm(n)
  y <- 0.5 * d + x[, 1] - x[, 2] + rnorm(n)
  fold <- rep(1:3, 100)
  a <- dmlx(y, d, x, "plr", fold, lrn = "glm", lrnd = "glm")
  s <- dmlx(y, d, x, "plr", fold, lrn = "stack", lrnd = "stack", stk = "glm")
  expect_identical(s$coef, a$coef)
  expect_identical(s$se, a$se)
  m <- dmlx(y, d, x, "plr", fold, lrn = "stack", stk = c("glm", "forest"),
            ntree = 20, nstk = 3)
  expect_lt(abs(m$coef - a$coef), 3 * a$se)
  t <- as.numeric(d > 0)
  i <- dmlx(y, t, x, "irm", fold, lrn = "stack", lrnd = "stack",
            stk = c("glm", "lasso"), nlam = 10, disc = TRUE)
  expect_true(is.finite(i$coef) && is.finite(i$se))
  errs("dmlx", list(
    list("names \"glm\" twice", y, d, x, "plr", fold, lrn = "stack",
         stk = c("glm", "glm")),
    list("no estimate", y[1:6], d[1:6], x[1:6, ], "plr", rep(1:2, 3),
         lrn = "stack", stk = "glm", nstk = 4)
  ))
})
