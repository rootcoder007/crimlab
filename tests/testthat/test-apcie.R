# The intrinsic estimator against glm() with sum-to-zero contrasts, its
# solution projected off the null vector of the design found by svd().

test_that("the estimates are the glm() solution off the null vector", {
  set.seed(11)
  A <- 5
  P <- 4
  ag <- c(-1, 0.5, 0.2, -0.3, -0.6)
  pe <- c(-0.2, 0, 0.1, 0.1)
  co <- seq(-0.4, 0.4, length.out = A + P - 1)
  ck <- outer(1:A, 1:P, function(i, j) j - i + A)
  pop <- matrix(round(runif(A * P, 5e4, 2e5)), A)
  y <- matrix(
    rpois(A * P, pop * exp(-6 + outer(ag, pe, "+") + matrix(co[ck], A))),
    A
  )
  r <- apcie(y, pop)
  expect_s3_class(r, "apcie")
  d <- data.frame(
    y = c(y),
    e = c(pop),
    a = factor(row(y)),
    p = factor(col(y)),
    k = factor(c(ck))
  )
  X <- model.matrix(
    ~ a + p + k,
    d,
    contrasts.arg = list(a = "contr.sum", p = "contr.sum", k = "contr.sum")
  )
  g <- glm(
    d$y ~ X - 1,
    family = poisson,
    offset = log(d$e),
    control = glm.control(epsilon = 1e-12, maxit = 100)
  )
  expect_true(g$converged)
  b <- coef(g)
  b[is.na(b)] <- 0
  V <- matrix(0, ncol(X), ncol(X))
  ok <- !is.na(coef(g))
  V[ok, ok] <- vcov(g, complete = FALSE)
  sv <- svd(X)
  b0 <- sv$v[, which.min(sv$d)]
  Pm <- diag(ncol(X)) - tcrossprod(b0)
  bie <- drop(Pm %*% b)
  vie <- Pm %*% V %*% Pm
  L <- function(nl) rbind(diag(nl - 1), -1)
  blk <- function(cols, nl) {
    m <- matrix(0, nl, ncol(X))
    m[, cols] <- L(nl)
    m
  }
  M <- rbind(
    c(1, rep(0, ncol(X) - 1)),
    blk(2:A, A),
    blk(A + 1:(P - 1), P),
    blk(A + P - 1 + 1:(A + P - 2), A + P - 1)
  )
  expect_equal(r$coef$est, drop(M %*% bie), tolerance = 1e-6)
  expect_equal(r$coef$se, sqrt(diag(M %*% vie %*% t(M))), tolerance = 1e-6)
  expect_equal(r$fit[["deviance"]], deviance(g), tolerance = 1e-8)
  expect_identical(r$fit[["df"]], A * P - (ncol(X) - 1))
  expect_identical(r$coef$term[1:3], c("(Intercept)", "age", "age"))
  expect_identical(sum(r$coef$term == "cohort"), as.integer(A + P - 1))
  expect_output(print(r), "5 ages, 4 periods, 8 cohorts")
})

test_that("dimnames label the age groups and periods", {
  y <- matrix(
    c(5, 8, 3, 9, 12, 6, 4, 10, 7),
    3,
    dimnames = list(c("a1", "a2", "a3"), c("p1", "p2", "p3"))
  )
  r <- apcie(y, matrix(100, 3, 3))
  expect_identical(r$coef$level[2:7], c("a1", "a2", "a3", "p1", "p2", "p3"))
})

test_that("apcie() refuses bad input in words", {
  errs(
    "apcie",
    list(
      list("`y` must be a numeric matrix", 1:9, matrix(1, 3, 3)),
      list("three ages and three periods", matrix(1, 2, 3), matrix(1, 2, 3)),
      list(
        "`expo` must be a matrix of the dimensions",
        matrix(1, 3, 3),
        matrix(1, 3, 2)
      ),
      list("`expo` must be greater than 0", matrix(1, 3, 3), matrix(0, 3, 3)),
      list("no events", matrix(0, 3, 3), matrix(1, 3, 3))
    )
  )
})
