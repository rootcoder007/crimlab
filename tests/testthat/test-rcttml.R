rt_dat <- function(n = 200, seed = 91) {
  set.seed(seed)
  d <- data.frame(prior = rpois(n, 2), age = runif(n, 18, 50),
                  t = rbinom(n, 1, 0.5))
  d$y <- 10 - 2 * d$t + 1.5 * d$prior - 0.05 * d$age + rnorm(n, sd = 2)
  d$b <- rbinom(n, 1, stats::plogis(-0.5 - 0.7 * d$t + 0.4 * d$prior))
  d
}

# The fluctuation written out with glm(): the clever covariate, the
# offset, and the influence curve.
rt_hand <- function(d, q1, q0, yv, g, a, b) {
  ys <- (yv - a) / (b - a)
  qb <- 0.005
  q1 <- pmin(pmax((q1 - a) / (b - a), qb), 1 - qb)
  q0 <- pmin(pmax((q0 - a) / (b - a), qb), 1 - qb)
  h <- d$t / g - (1 - d$t) / (1 - g)
  qd <- ifelse(d$t == 1, q1, q0)
  fit <- suppressWarnings(stats::glm(
    ys ~ -1 + h + offset(stats::qlogis(qd)), stats::quasibinomial(),
    control = stats::glm.control(epsilon = 1e-14, maxit = 100)
  ))
  eps <- unname(stats::coef(fit))
  s1 <- stats::plogis(stats::qlogis(q1) + eps / g)
  s0 <- stats::plogis(stats::qlogis(q0) - eps / (1 - g))
  psi <- mean(s1 - s0)
  ic <- h * (ys - ifelse(d$t == 1, s1, s0)) + s1 - s0 - psi
  c(est = psi * (b - a), se = sqrt(stats::var(ic) / nrow(d)) * (b - a),
    eps = eps)
}

test_that("rcttml() fluctuates lm()'s Q along the clever covariate", {
  d <- rt_dat()
  r <- rcttml(d, "y", "t", c("prior", "age"))
  fit <- stats::lm(y ~ t + prior + age, d)
  q1 <- stats::predict(fit, transform(d, t = 1))
  q0 <- stats::predict(fit, transform(d, t = 0))
  w <- rt_hand(d, q1, q0, d$y, mean(d$t), min(d$y), max(d$y))
  expect_identical(r$tab$method, c("tmle", "gcomp", "unadj"))
  expect_equal(r$tab$est[1L], w[["est"]], tolerance = 1e-7)
  expect_equal(r$tab$se[1L], w[["se"]], tolerance = 1e-7)
  expect_equal(r$eps, w[["eps"]], tolerance = 1e-6)
  expect_equal(r$tab$est[2L], unname(stats::coef(fit)[["t"]]),
               tolerance = 1e-10)
  expect_true(is.na(r$tab$se[2L]))
  m1 <- d$y[d$t == 1]
  m0 <- d$y[d$t == 0]
  expect_equal(r$tab$est[3L], mean(m1) - mean(m0), tolerance = 1e-12)
  expect_equal(r$tab$se[3L],
               sqrt(stats::var(m1) / length(m1) + stats::var(m0) / length(m0)),
               tolerance = 1e-12)
  z <- stats::qnorm(0.975)
  expect_equal(r$tab$lo, r$tab$est - z * r$tab$se, tolerance = 1e-15)
  expect_equal(r$tab$p[1L], 2 * stats::pnorm(-abs(w[["est"]] / w[["se"]])),
               tolerance = 1e-6)
  expect_equal(r$g, rep(mean(d$t), 2), tolerance = 1e-15)
})

test_that("rcttml() takes interactions, a known g and a fitted g", {
  d <- rt_dat(seed = 92)
  fit <- stats::lm(y ~ t * prior, d)
  q1 <- stats::predict(fit, transform(d, t = 1))
  q0 <- stats::predict(fit, transform(d, t = 0))
  r <- rcttml(d, "y", "t", "prior", inter = TRUE, g = 0.5,
              bounds = c(-20, 40))
  w <- rt_hand(d, q1, q0, d$y, 0.5, -20, 40)
  expect_equal(r$tab$est[1L], w[["est"]], tolerance = 1e-7)
  expect_equal(r$tab$est[2L], mean(q1 - q0), tolerance = 1e-10)
  expect_identical(r$bounds, c(-20, 40))
  gf <- stats::fitted(stats::glm(t ~ prior + age, stats::binomial(), d))
  x <- rcttml(d, "y", "t", c("prior", "age"), gx = TRUE)
  fq <- stats::lm(y ~ t + prior + age, d)
  wx <- rt_hand(d, stats::predict(fq, transform(d, t = 1)),
                stats::predict(fq, transform(d, t = 0)), d$y, unname(gf),
                min(d$y), max(d$y))
  expect_equal(x$tab$est[1L], wx[["est"]], tolerance = 1e-6)
  expect_equal(x$g, range(gf), tolerance = 1e-7)
})

test_that("rcttml() fits a binary outcome by logistic regression", {
  d <- rt_dat(seed = 93)
  fit <- stats::glm(b ~ t + prior, stats::binomial(), d)
  q1 <- stats::predict(fit, transform(d, t = 1), type = "response")
  q0 <- stats::predict(fit, transform(d, t = 0), type = "response")
  r <- rcttml(d, "b", "t", "prior", fam = "binomial")
  w <- rt_hand(d, q1, q0, d$b, mean(d$t), 0, 1)
  expect_equal(r$tab$est[1L], w[["est"]], tolerance = 1e-6)
  expect_equal(r$tab$se[1L], w[["se"]], tolerance = 1e-6)
  expect_equal(r$tab$est[2L], mean(q1 - q0), tolerance = 1e-7)
  out <- utils::capture.output(v <- print(r))
  expect_identical(out[1L], paste("rcttml: binomial outcome, 200",
                                  "participants; targeted maximum",
                                  "likelihood, g-computation and the",
                                  "difference in means"))
  expect_identical(v, r)
})

test_that("rcttml() refuses trials it cannot adjust", {
  d <- rt_dat()
  expect_error(rcttml(d, "y", "prior"),
               "rcttml: column `prior` must hold 0 or 1.", fixed = TRUE)
  expect_error(rcttml(d, "y", "t", fam = "binomial"),
               paste("rcttml: with `fam = \"binomial\"`, column `y` must hold",
                     "0 or 1."), fixed = TRUE)
  expect_error(rcttml(d, "y", "t", "prior", g = 0.5, gx = TRUE),
               "rcttml: give `g` or `gx`, not both.", fixed = TRUE)
  expect_error(rcttml(d, "y", "t", gx = TRUE),
               "rcttml: `gx` needs covariates `x`.", fixed = TRUE)
  expect_error(rcttml(d, "y", "t", bounds = c(0, 1)),
               "rcttml: `bounds` must be two increasing numbers around `y`.",
               fixed = TRUE)
  expect_error(rcttml(d[1:3, ], "y", "t"),
               paste("rcttml: column `t` must have two participants or more",
                     "in each arm."), fixed = TRUE)
  e <- d
  e$p2 <- 2 * e$prior
  expect_error(rcttml(e, "y", "t", c("prior", "p2")),
               paste("rcttml: the outcome or treatment model cannot be",
                     "fitted: the covariates are collinear, or separate the",
                     "outcome or the arms."), fixed = TRUE)
})
