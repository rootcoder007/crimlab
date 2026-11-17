# The count model with a random intercept per cluster (cntmix()), checked
# against the exact marginal likelihood and against repeated rows.

test_that("cntmix() finds the maximum of the marginal likelihood", {
  set.seed(3)
  g <- rep(1:4, each = 6)
  d <- rep(0:1, 12)
  y <- rpois(24, exp(-0.3 + 0.5 * d + c(-0.6, 0, 0.4, 0.8)[g]))
  f <- cntmix(y, data.frame(d = d), g, fam = "poisson")
  # the exact marginal likelihood by quadrature, at (b0, b1, sigma)
  ex <- function(b0, b1, s) {
    sum(vapply(
      1:4,
      function(j) {
        i <- g == j
        log(
          stats::integrate(
            function(u) {
              vapply(
                u,
                function(v) {
                  exp(sum(stats::dpois(
                    y[i],
                    exp(b0 + b1 * d[i] + v),
                    log = TRUE
                  )))
                },
                0
              ) *
                stats::dnorm(u, 0, s)
            },
            -8 * s,
            8 * s
          )$value
        )
      },
      0
    ))
  }
  b <- f$coef$est
  l0 <- ex(b[1], b[2], f$sigma)
  # the Laplace value is close to the exact one at the estimate
  expect_lt(abs(f$loglik - l0), 0.05)
  # and moving the treatment coefficient lowers the exact likelihood
  expect_lt(ex(b[1], b[2] + 0.1, f$sigma), l0)
  expect_lt(ex(b[1], b[2] - 0.1, f$sigma), l0)
  # relabelling the clusters changes nothing
  f2 <- cntmix(y, data.frame(d = d), letters[5 - g], fam = "poisson")
  expect_equal(f2$coef$est, f$coef$est, tolerance = 1e-6)
  errs(
    "cntmix",
    list(
      list("two clusters or more", y, data.frame(d = d), rep(1, 24)),
      list("0 throughout", 0 * y, data.frame(d = d), g),
      list("one complete row", y, data.frame(d = d[-1]), g),
      list("must be one of", y, data.frame(d = d), g, fam = "x")
    )
  )
})

test_that("cntmix() takes weights and a matrix, and says when theta runs off", {
  set.seed(1)
  g <- rep(1:6, each = 10)
  d <- rep(0:1, 30)
  # binomial counts are less dispersed than Poisson: theta runs off
  y <- stats::rbinom(
    60,
    3,
    stats::plogis(-0.5 + 0.4 * d + c(-0.3, 0, 0.2, 0.4, -0.1, 0.1)[g])
  )
  x <- data.frame(d = d)
  fp <- cntmix(y, x, g, fam = "poisson")
  # integer weights are the rows repeated
  w <- rep(1:3, 20)
  i <- rep(seq_along(y), w)
  fw <- cntmix(y, x, g, fam = "poisson", w = w)
  fr <- cntmix(y[i], data.frame(d = d[i]), g[i], fam = "poisson")
  expect_equal(
    fw$coef[c("est", "se")],
    fr$coef[c("est", "se")],
    tolerance = 1e-6
  )
  expect_equal(fw$sigma, fr$sigma, tolerance = 1e-6)
  expect_equal(fw$loglik, fr$loglik, tolerance = 1e-8)
  # a numeric matrix is its data frame
  expect_equal(cntmix(y, cbind(d = d), g, fam = "poisson"), fp)
  # the Poisson limit: the Poisson fit, theta named in fam
  fn <- cntmix(y, x, g, fam = "nb")
  expect_match(fn$fam, "^poisson \\(negative binomial theta .*Poisson limit\\)")
  expect_gt(as.numeric(sub(".*theta ([^:]+):.*", "\\1", fn$fam)), 1e6)
  expect_identical(fn$theta, NA_real_)
  expect_equal(fn$coef$est, fp$coef$est, tolerance = 1e-4)
  expect_equal(fn$coef$se, fp$coef$se, tolerance = 1e-3)
  # where the likelihood overflows, a very low value
  expect_identical(
    .cmlap(c(800, 0, 0), y, cbind(1, d), g, rep(1, 60), "poisson"),
    -1e300
  )
  errs(
    "cntmix",
    list(
      list("one value per count", y, x, g, w = 1:3),
      list("a cluster for every count", y, x, g[-1]),
      list("a cluster for every count", y, x, replace(g, 1, NA))
    )
  )
  prt(fn, "^cntmix: poisson \\(negative binomial theta .*\n60 counts in 6 ")
  expect_false(any(grepl("theta [0-9.]+;", capture.output(print(fn)))))
  # an overdispersed sample keeps the negative binomial and shows theta
  yn <- stats::rnbinom(60, size = 1, mu = exp(0.5 + 0.4 * d))
  fb <- cntmix(yn, x, g, fam = "nb")
  expect_identical(fb$fam, "nb")
  prt(fb, "^cntmix: nb counts with a random intercept per cluster\n60 counts")
  expect_output(print(fb), paste0("theta ", format(fb$theta, digits = 4)))
  fb$conv <- 1L
  expect_output(print(fb), "\\(the optimiser did not converge\\)")
  # a singular Hessian (as when a column is given twice): no standard
  # error, the estimates as they were
  fit <- .cmfit
  local_mocked_bindings(.cmfit = function(...) {
    f <- fit(...)
    f$hess[] <- 0
    f
  })
  fs <- cntmix(y, x, g, fam = "poisson")
  expect_true(all(is.na(fs$coef[c("se", "z", "p", "lo", "hi")])))
  expect_identical(fs$coef$est, fp$coef$est)
})

test_that("cntmix() reports an aliased column as NA, with the reason", {
  set.seed(5)
  n <- 200
  g <- sample(10, n, TRUE)
  x <- rnorm(n)
  y <- rpois(n, exp(0.2 + 0.3 * x + rnorm(10, 0, 0.5)[g]))
  f <- cntmix(y, data.frame(x, k = 5, x2 = 2 * x), g, fam = "poisson")
  a <- cntmix(y, data.frame(x), g, fam = "poisson")
  expect_identical(f$coef$term, c("(Intercept)", "x", "k", "x2"))
  expect_identical(anyNA(f$coef[1:2, ]), FALSE)
  expect_identical(all(is.na(f$coef[3:4, -1])), TRUE)
  expect_equal(f$coef[1:2, ], a$coef, tolerance = 1e-6)
  expect_identical(f$notes, paste("k, x2: aliased with the other columns",
                                  "(constant, or collinear), dropped"))
  expect_identical(a$notes, character(0))
  expect_output(print(f), "note: k, x2: aliased")
  errs("cntmix", list(
    list("`level` must be one number strictly between 0 and 1", y,
         data.frame(x), g, level = 1),
    list("`level` must be one number strictly between 0 and 1", y,
         data.frame(x), g, level = 0)
  ))
})

test_that("cntmix() says when a 0/1 column separates the counts", {
  set.seed(6)
  n <- 200
  g <- sample(8, n, TRUE)
  d <- rbinom(n, 1, 0.3)
  f <- cntmix(rpois(n, 2) * (1 - d), data.frame(d), g, fam = "poisson")
  expect_identical(is.na(f$coef$se), c(FALSE, TRUE))
  expect_identical(f$notes,
                   "d: every count where it is 1 is 0, so its rate ratio is 0")
  f <- cntmix(rpois(n, 2) * d, data.frame(d), g, fam = "poisson")
  expect_identical(is.na(f$coef$est), c(TRUE, TRUE))
  expect_match(f$notes, "where it is 0 is 0, so its rate ratio is infinite")
})

test_that("a variance from the Hessian that is not positive is NA", {
  X <- cbind("(Intercept)" = 1, a = c(0, 1, 0, 1))
  e <- list(keep = 1:2, na = character(0), notes = NULL)
  fit <- list(par = c(0.1, 0.2, 0), hess = diag(c(1, -1, 1)), ll = -1,
              conv = 0, fam = "poisson", nb = FALSE)
  o <- .cmout(fit, X, e, rep(1, 4), c(1, 1, 2, 2), 0.95)
  expect_identical(is.na(o$coef$se), c(FALSE, TRUE))
  expect_identical(o$notes, paste("a: the Hessian gives no positive",
                                  "variance, so no standard error"))
  fit$hess <- matrix(0, 3, 3)
  o <- .cmout(fit, X, e, rep(1, 4), c(1, 1, 2, 2), 0.95)
  expect_identical(is.na(o$coef$se), c(TRUE, TRUE))
})
