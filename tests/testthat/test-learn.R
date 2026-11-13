# The learners against glmnet and glm(), with reference values from
# dev/ref-learn.R; the kernels' own checks are in src/t_learn.cpp.
ref <- readRDS(test_path("ref-learn.rds"))

test_that("enpath agrees with glmnet on the same lambdas and folds", {
  for (f in ref$fits) {
    yy <- if (f$fam == "gaussian") ref$y else ref$yb
    r <- enpath(ref$x, yy, f$fam, f$alpha, lambda = f$lambda,
               fold = ref$fold, tol = 1e-20, maxit = 1e7)
    info <- paste(f$fam, f$alpha)
    # coefficients and losses to 1e-8: both converged with threshold 1e-20
    expect_lt(max(abs(r$beta - f$beta)), 1e-8, label = info)
    expect_lt(max(abs(r$a0 - f$a0)), 1e-8, label = info)
    expect_lt(max(abs(r$cvm - f$cvm)), 1e-8, label = info)
    expect_lt(max(abs(r$cvsd - f$cvsd)), 1e-8, label = info)
    expect_identical(r$lmin, f$lmin, info = info)
    expect_identical(r$l1se, f$l1se, info = info)
    expect_identical(rownames(r$beta), colnames(ref$x))
    # the path from the data is glmnet's
    a <- enpath(ref$x, yy, f$fam, f$alpha, nlam = 15, lmr = 0.02)
    # (both stop early by the same rule)
    expect_equal(a$lambda, f$lambda, tolerance = 1e-12, info = info)
  }
})

test_that("lgtreg agrees with glm", {
  r <- lgtreg(ref$x, ref$yb)
  expect_equal(unname(r$coef), ref$glm$coef, tolerance = 1e-10)
  # glm() takes its covariance from the weights of its last iteration
  expect_equal(unname(r$se), ref$glm$se, tolerance = 1e-7)
  expect_equal(r$dev, ref$glm$dev, tolerance = 1e-12)
  expect_named(r$coef, c("(Intercept)", colnames(ref$x)))
  expect_identical(names(lgtreg(unname(ref$x), ref$yb)$se)[2:3],
                   c("x1", "x2"))
})

test_that("the learners stop on bad input, naming the argument", {
  x <- ref$x[1:40, ]
  y <- ref$y[1:40]
  yb <- ref$yb[1:40]
  errs("enpath", list(
    list("`x` must be a numeric matrix", as.data.frame(x), y),
    list("`x` must be a numeric matrix", x[, 1], y),
    list("`x` must be finite", replace(x, 3, NA), y),
    list("two rows and one column", x[1, , drop = FALSE], y[1]),
    list("`y` must have length 40", x, y[-1]),
    list("`fam` must be one of", x, y, "poisson"),
    list("`fam` must be one of", x, y, NA_character_),
    list("`y` must be 0 or 1", x, y, "binomial"),
    list("`y` is constant", x, rep(1, 40)),
    list("`y` is constant", x, rep(1, 40), "binomial"),
    list("`alpha`", x, y, alpha = 2),
    list("`nlam`", x, y, nlam = 0),
    list("`lmr`", x, y, lmr = 1),
    list("`lambda` must be decreasing", x, y, lambda = c(0.1, 0.2)),
    list("`lambda` must be at least 0", x, y, lambda = -1),
    list("`fold` must be whole", x, y, fold = rep(1.5, 40)),
    list("at least two folds", x, y, fold = rep(1, 40)),
    list("leaves fold 2 empty", x, y, fold = rep(c(1, 3), 20)),
    list("`tol`", x, y, tol = 0),
    list("`maxit`", x, y, maxit = 0.5),
    list("no fit", x, y, maxit = 1, tol = 1e-30)
  ))
  sep <- as.numeric(x[, 1] > mean(x[, 1]))
  errs("lgtreg", list(
    list("`x` must be a numeric matrix", "a", yb),
    list("`y` must be 0 or 1", x, y),
    list("`y` is constant", x, rep(0, 40)),
    list("more rows than columns", x[1:9, ], rep(0:1, length.out = 9)),
    list("`tol`", x, yb, tol = -1),
    list("`maxit`", x, yb, maxit = NA),
    list("no fit", x[, 1, drop = FALSE], sep)
  ))
})
