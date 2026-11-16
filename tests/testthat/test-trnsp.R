tp_dat <- function(n = 500, seed = 81) {
  set.seed(seed)
  d <- data.frame(age = rnorm(n, 35, 10), prior = rpois(n, 1.5))
  d$s <- rbinom(n, 1, plogis(-0.04 * (d$age - 35) + 0.2 * d$prior - 0.3))
  d$t <- ifelse(d$s == 1, rbinom(n, 1, 0.5), NA)
  d$y <- ifelse(d$s == 1, 0.5 * d$t + 0.03 * d$t * d$age + 0.2 * d$prior +
                  rnorm(n), NA)
  d
}

test_that("trnsp()'s outcome model is lm() by arm averaged on the target", {
  d <- tp_dat()
  r <- trnsp(d, "y", "t", "s", c("age", "prior"), nboot = 0)
  tr <- d[d$s == 1, ]
  f1 <- stats::lm(y ~ age + prior, tr[tr$t == 1, ])
  f0 <- stats::lm(y ~ age + prior, tr[tr$t == 0, ])
  tg <- d[d$s == 0, ]
  om <- mean(stats::predict(f1, tg) - stats::predict(f0, tg))
  expect_equal(r$tab$est[1L], om, tolerance = 1e-10)
  g <- trnsp(d, "y", "t", "s", c("age", "prior"), target = "generalize",
             nboot = 0)
  expect_equal(g$tab$est[1L],
               mean(stats::predict(f1, d) - stats::predict(f0, d)),
               tolerance = 1e-10)
  expect_identical(r$tab$method, c("om", "iow"))
  expect_true(all(is.na(r$tab$se)))
  expect_identical(r$seed, NA_real_)
})

test_that("trnsp() weights the trial by the inverse odds of glm()'s score", {
  d <- tp_dat(seed = 82)
  r <- trnsp(d, "y", "t", "s", c("age", "prior"), nboot = 0)
  ps <- unname(stats::fitted(stats::glm(s ~ age + prior, stats::binomial(),
                                        d)))
  w <- (1 - ps) / ps
  tr <- d$s == 1
  m1 <- tr & d$t %in% 1
  m0 <- tr & d$t %in% 0
  iow <- sum(w[m1] * d$y[m1]) / sum(w[m1]) - sum(w[m0] * d$y[m0]) / sum(w[m0])
  expect_equal(r$tab$est[2L], iow, tolerance = 1e-7)
  expect_equal(r$ps, ps, tolerance = 1e-7)
  expect_equal(r$w, ifelse(tr, w, 0), tolerance = 1e-6)
  g <- trnsp(d, "y", "t", "s", c("age", "prior"), target = "generalize",
             nboot = 0)
  wg <- 1 / ps
  ig <- sum(wg[m1] * d$y[m1]) / sum(wg[m1]) -
    sum(wg[m0] * d$y[m0]) / sum(wg[m0])
  expect_equal(g$tab$est[2L], ig, tolerance = 1e-7)
})

test_that("trnsp()'s errors are the spread of its bootstrap draws", {
  d <- tp_dat(seed = 83)
  r <- trnsp(d, "y", "t", "s", "age", nboot = 60, seed = 4)
  expect_identical(dim(r$draws), c(60L, 2L))
  expect_identical(r$nfail, sum(is.na(r$draws[, 1L])) + 0)
  expect_equal(r$tab$se, unname(apply(r$draws, 2L, stats::sd, na.rm = TRUE)),
               tolerance = 1e-12)
  z <- stats::qnorm(0.975)
  expect_equal(r$tab$lo, r$tab$est - z * r$tab$se, tolerance = 1e-15)
  expect_identical(trnsp(d, "y", "t", "s", "age", nboot = 60, seed = 4), r)
  set.seed(5)
  a <- trnsp(d, "y", "t", "s", "age", nboot = 5)
  expect_identical(trnsp(d, "y", "t", "s", "age", nboot = 5, seed = a$seed),
                   a)
  out <- utils::capture.output(v <- print(r))
  expect_identical(out[1L], paste("trnsp: the trial's effect transported to",
                                  "the rows outside it, 60 bootstrap draws",
                                  paste0("(", r$nfail, " failed)")))
  expect_identical(v, r)
})

test_that("trnsp() refuses trials and targets it cannot use", {
  d <- tp_dat()
  expect_error(trnsp(d, "y", "t", "s", "y"),
               "trnsp: column `y` has missing values.", fixed = TRUE)
  expect_error(trnsp(d, "y", "t", "s", "s"),
               "trnsp: `y`, `d`, `s` and `x` must name different columns.",
               fixed = TRUE)
  expect_error(trnsp(d, "w", "t", "s", "age"),
               "trnsp: `y` names columns not in `data`: w.", fixed = TRUE)
  e <- d
  e$y[which(e$s == 1)[1L]] <- NA
  expect_error(trnsp(e, "y", "t", "s", "age"),
               "trnsp: column `y` must be finite in the trial.", fixed = TRUE)
  e <- d
  e$t[e$s == 1] <- 1
  expect_error(trnsp(e, "y", "t", "s", "age"),
               "trnsp: column `t` must hold both arms in the trial.",
               fixed = TRUE)
  e <- d
  e$s <- 1
  expect_error(trnsp(e, "y", "t", "s", "age"),
               "trnsp: column `s` must mark units in and out of the trial.",
               fixed = TRUE)
  expect_error(trnsp(d, "y", "t", "s", "age", nboot = -1),
               "trnsp: `nboot` must be between 0 and 1e+05.", fixed = TRUE)
})
