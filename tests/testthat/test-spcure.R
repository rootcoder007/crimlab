# The split-population model. The kernel is checked against a closed
# form in src/t_spc.cpp; here against the likelihood written out in R and
# maximised by optim(), with optimHess() for the standard errors.

spsim <- function() {
  set.seed(10)
  n <- 500
  prior <- rbinom(n, 1, 0.4)
  fail <- rbinom(n, 1, plogis(-0.5 + 1.2 * prior))
  tt <- ifelse(fail == 1, rweibull(n, 1.3, 300 * exp(-0.3 * prior)), Inf)
  cens <- runif(n, 365, 1500)
  list(t = pmin(tt, cens), e = as.numeric(tt <= cens), x = prior)
}

spll <- function(u, t, e, x, dist) {
  d <- plogis(u[1] + u[2] * x)
  eta <- u[3] + u[4] * x
  if (dist == "lognormal") {
    s <- exp(u[5])
    lf <- dlnorm(t, eta, s, log = TRUE)
    S <- plnorm(t, eta, s, lower.tail = FALSE)
  } else {
    k <- if (dist == "exp") 1 else exp(u[5])
    lf <- dweibull(t, k, exp(eta), log = TRUE)
    S <- pweibull(t, k, exp(eta), lower.tail = FALSE)
  }
  -sum(ifelse(e == 1, log(d) + lf, log(1 - d + d * S)))
}

test_that("the fit maximises the split-population likelihood", {
  s <- spsim()
  for (dist in c("weibull", "exp", "lognormal")) {
    r <- spcure(s$t, s$e, data.frame(prior = s$x), dist = dist)
    expect_s3_class(r, "spcure")
    expect_identical(r$fit[["conv"]], 1)
    u <- r$coef$est
    o <- optim(
      u,
      spll,
      t = s$t,
      e = s$e,
      x = s$x,
      dist = dist,
      method = "BFGS",
      control = list(reltol = 1e-14, maxit = 1000)
    )
    expect_equal(
      r$fit[["loglik"]],
      -spll(u, s$t, s$e, s$x, dist),
      tolerance = 1e-10
    )
    expect_lte(-r$fit[["loglik"]], o$value + 1e-6)
    expect_equal(u, o$par, tolerance = 1e-3)
    h <- optimHess(u, spll, t = s$t, e = s$e, x = s$x, dist = dist)
    expect_equal(r$coef$se, sqrt(diag(solve(h))), tolerance = 1e-3)
  }
  expect_identical(r$coef$part, c("fail", "fail", "time", "time", "shape"))
  expect_identical(
    r$coef$term,
    c("(Intercept)", "prior", "(Intercept)", "prior", "log(sigma)")
  )
  expect_output(print(r), "share who would ever fail")
})

test_that("spcure() refuses bad input in words", {
  errs(
    "spcure",
    list(
      list("`time`", c(0, 1), c(1, 0)),
      list("failures and censored times both", c(1, 2), c(1, 1)),
      list("`dist`", c(1, 2), c(1, 0), dist = "gamma"),
      list(
        "`x` must be a numeric matrix",
        c(1, 2),
        c(1, 0),
        x = matrix(1, 3, 1)
      )
    )
  )
})

test_that("a covariate matrix is the data frame of its columns", {
  set.seed(1)
  n <- 200
  x <- rnorm(n)
  t <- rexp(n)
  ev <- rbinom(n, 1, 0.6)
  a <- spcure(t, ev, data.frame(prior = x))
  expect_identical(spcure(t, ev, cbind(prior = x))$coef, a$coef)
  # an unnamed column is named by its position
  expect_identical(spcure(t, ev, matrix(x))$coef$term,
                   c("(Intercept)", "x1", "(Intercept)", "x1", "log(shape)"))
  expect_error(spcure(t, ev, matrix(c(x[-1], Inf))),
               "spcure: `x` must hold finite numbers with no missing values.",
               fixed = TRUE)
})
