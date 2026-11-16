mp_dat <- function(n = 120, seed = 101) {
  set.seed(seed)
  d <- data.frame(age = round(runif(n, 18, 60)), male = rbinom(n, 1, 0.5))
  d$priors <- rpois(n, exp(-1 + 0.04 * d$age))
  d$days <- round(20 + 3 * d$priors + 5 * d$male + rnorm(n, sd = 5))
  d$priors[sample(n, 25)] <- NA
  d$days[sample(n, 20)] <- NA
  d
}

# Rubin's rules written out, with the degrees of freedom of Barnard and
# Rubin (1999) as the mice package computes them.
mp_pool <- function(q, u, m, dfc) {
  qb <- rowMeans(q)
  ub <- rowMeans(u)
  b <- apply(q, 1L, stats::var)
  t <- ub + (1 + 1 / m) * b
  lam <- (1 + 1 / m) * b / t
  tmp <- (1 - lam) * (1 + dfc) * dfc
  df <- (m - 1) * tmp / ((dfc + 3) * (m - 1) + lam^2 * tmp)
  riv <- (1 + 1 / m) * b / ub
  list(est = unname(qb), se = unname(sqrt(t)), df = unname(df),
       riv = unname(riv), fmi = unname((riv + 2 / (df + 3)) / (riv + 1)))
}

test_that("mipmm() takes the nearest donor's value when the fit is exact", {
  d <- data.frame(x = c(1, 2, 4, 7, 11, 16, 22, 29, 37, 46, 3.2, 20))
  d$y <- 2 * d$x
  d$y[11:12] <- NA
  r <- mipmm(d, "y", "x", m = 3, k = 1, ridge = 0, seed = 1)
  # fitted 6.4 is nearest 8 (x = 4); fitted 40 is nearest 44 (x = 22)
  expect_identical(r$imp$y, matrix(c(8, 44), 2L, 3L,
                                   dimnames = list(c("11", "12"), NULL)))
  expect_identical(r$nmis, c(y = 2))
})

test_that("mipmm() imputes observed values and pools a fit by Rubin's rules", {
  d <- mp_dat()
  r <- mipmm(d, c("priors", "days"), c("age", "male"), m = 6, fy = "days",
             fx = c("priors", "male"), seed = 3)
  expect_identical(names(r$imp), c("priors", "days"))
  expect_identical(dim(r$imp$priors), c(25L, 6L))
  expect_identical(rownames(r$imp$days), as.character(which(is.na(d$days))))
  expect_true(all(r$imp$priors %in% d$priors[!is.na(d$priors)]))
  expect_true(all(r$imp$days %in% d$days[!is.na(d$days)]))
  fits <- lapply(seq_len(6), function(j) {
    e <- d
    e$priors[is.na(d$priors)] <- r$imp$priors[, j]
    e$days[is.na(d$days)] <- r$imp$days[, j]
    stats::lm(days ~ priors + male, e)
  })
  q <- sapply(fits, stats::coef)
  u <- sapply(fits, function(f) diag(stats::vcov(f)))
  w <- mp_pool(q, u, 6, nrow(d) - 3)
  p <- r$pool
  expect_identical(p$term, c("(Intercept)", "priors", "male"))
  expect_equal(p$est, w$est, tolerance = 1e-10)
  expect_equal(p$se, w$se, tolerance = 1e-10)
  expect_equal(p$df, w$df, tolerance = 1e-10)
  expect_equal(p$riv, w$riv, tolerance = 1e-10)
  expect_equal(p$fmi, w$fmi, tolerance = 1e-10)
  ct <- stats::qt(0.975, w$df)
  expect_equal(p$lo, w$est - ct * w$se, tolerance = 1e-10)
  expect_equal(p$p, 2 * stats::pt(-abs(w$est / w$se), w$df),
               tolerance = 1e-8)
})

test_that("mipmm() pools the means of the imputed columns by default", {
  d <- mp_dat(seed = 102)
  r <- mipmm(d, c("priors", "days"), "age", m = 4, k = 3, maxit = 2,
             seed = 4)
  cm <- sapply(1:4, function(j) {
    c(mean(replace(d$priors, is.na(d$priors), r$imp$priors[, j])),
      mean(replace(d$days, is.na(d$days), r$imp$days[, j])))
  })
  cv <- sapply(1:4, function(j) {
    c(stats::var(replace(d$priors, is.na(d$priors), r$imp$priors[, j])),
      stats::var(replace(d$days, is.na(d$days), r$imp$days[, j]))) /
      nrow(d)
  })
  w <- mp_pool(cm, cv, 4, nrow(d) - 1)
  expect_identical(r$pool$term, c("priors", "days"))
  expect_equal(r$pool$est, w$est, tolerance = 1e-12)
  expect_equal(r$pool$se, w$se, tolerance = 1e-10)
  expect_equal(r$pool$df, w$df, tolerance = 1e-10)
  expect_identical(mipmm(d, c("priors", "days"), "age", m = 4, k = 3,
                         maxit = 2, seed = 4), r)
  set.seed(6)
  a <- mipmm(d, "priors", "age", m = 2)
  expect_identical(mipmm(d, "priors", "age", m = 2, seed = a$seed), a)
  out <- utils::capture.output(v <- print(r))
  expect_identical(out[1L], paste("mipmm: 45 missing value(s) in 2",
                                  "column(s), 4 imputations by predictive",
                                  "mean matching"))
  expect_identical(v, r)
})

test_that("mipmm() refuses columns and analyses it cannot use", {
  d <- mp_dat()
  expect_error(mipmm(d, c("priors", "age"), "age"),
               "mipmm: `y` and `x` must name different columns.",
               fixed = TRUE)
  expect_error(mipmm(d, "priors", "days"),
               "mipmm: column `days` has missing values.", fixed = TRUE)
  e <- d
  e$days[1L] <- Inf
  expect_error(mipmm(e, "days", "age"),
               "mipmm: column `days` must be finite where it is observed.",
               fixed = TRUE)
  expect_error(mipmm(d, "priors", "age", fy = "male"),
               "mipmm: `fy` must name one column of `y` or `x`.",
               fixed = TRUE)
  expect_error(mipmm(d, "priors", "age", fx = "age"),
               paste("mipmm: `fx` must name columns of `y` or `x` other than",
                     "`fy`, each once."), fixed = TRUE)
  expect_error(mipmm(d, "priors", "age", m = 1),
               "mipmm: `m` must be between 2 and 1000.", fixed = TRUE)
  f <- data.frame(a = c(1, NA, NA, NA), b = 1:4)
  expect_error(mipmm(f, "a", "b", seed = 1),
               paste("mipmm: a column of `y` has too few observed values for",
                     "its predictors, or they are collinear."), fixed = TRUE)
})
