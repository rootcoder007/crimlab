gw_dat <- function(n = 400, seed = 71) {
  set.seed(seed)
  d <- data.frame(age = rnorm(n, 35, 10), prior = rpois(n, 1.5))
  d$s <- rbinom(n, 1, plogis(-0.5 + 0.04 * (d$age - 35) + 0.3 * d$prior))
  d
}

test_that("gensw() fits the sampling score of glm() and weights by it", {
  d <- gw_dat()
  fit <- stats::glm(s ~ age + prior, stats::binomial(), d)
  ps <- unname(stats::fitted(fit))
  g <- gensw(d, "s", c("age", "prior"))
  expect_equal(g$coef, stats::coef(fit), tolerance = 1e-7)
  expect_equal(g$ps, ps, tolerance = 1e-7)
  expect_equal(g$w, ifelse(d$s == 1, 1 / g$ps, 0), tolerance = 1e-15)
  t <- gensw(d, "s", c("age", "prior"), target = "transport")
  expect_equal(t$w, ifelse(d$s == 1, (1 - t$ps) / t$ps, 0),
               tolerance = 1e-15)
  expect_identical(t$target, "transport")
})

test_that("gensw()'s diagnostics are their definitions written out", {
  d <- gw_dat(seed = 72)
  g <- gensw(d, "s", c("age", "prior"), target = "transport", bins = 5)
  ps <- g$ps
  tr <- d$s == 1
  tg <- d$s == 0
  l <- stats::qlogis(ps)
  expect_identical(g$diag[["ntrial"]], sum(tr) + 0)
  expect_identical(g$diag[["ntarget"]], sum(tg) + 0)
  expect_equal(g$diag[["smd"]], (mean(l[tr]) - mean(l[tg])) / stats::sd(l[tg]),
               tolerance = 1e-12)
  expect_equal(g$diag[["dmean"]], mean(ps[tr]) - mean(ps[tg]),
               tolerance = 1e-12)
  b <- pmin(floor((ps - min(ps)) / (max(ps) - min(ps)) * 5), 4)
  f1 <- tabulate(b[tr] + 1, 5) / sum(tr)
  ft <- tabulate(b[tg] + 1, 5) / sum(tg)
  expect_equal(g$diag[["gindex"]], sum(sqrt(f1 * ft)), tolerance = 1e-12)
  w <- g$w[tr]
  expect_equal(g$diag[["ess"]], sum(w)^2 / sum(w^2), tolerance = 1e-12)
  bal <- g$bal
  expect_identical(bal$var, c("age", "prior"))
  for (j in 1:2) {
    v <- d[[bal$var[j]]]
    sdt <- stats::sd(v[tg])
    expect_equal(bal$target[j], mean(v[tg]), tolerance = 1e-12)
    expect_equal(bal$trial[j], mean(v[tr]), tolerance = 1e-12)
    expect_equal(bal$wtrial[j], sum(w * v[tr]) / sum(w), tolerance = 1e-12)
    expect_equal(bal$smd[j], (mean(v[tr]) - mean(v[tg])) / sdt,
                 tolerance = 1e-12)
    expect_equal(bal$wsmd[j], (sum(w * v[tr]) / sum(w) - mean(v[tg])) / sdt,
                 tolerance = 1e-12)
  }
  # generalising, the target is every row
  a <- gensw(d, "s", c("age", "prior"))
  expect_identical(a$diag[["ntarget"]], nrow(d) + 0)
  expect_equal(a$bal$target, colMeans(d[c("age", "prior")]),
               tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("gensw() trims the scores and prints its diagnostics", {
  d <- gw_dat(seed = 73)
  g <- gensw(d, "s", c("age", "prior"), trim = 0.3)
  expect_gte(min(g$ps), 0.3)
  expect_lte(max(g$ps), 0.7)
  out <- utils::capture.output(v <- print(g))
  expect_match(out[1L], paste0("^gensw: ", sum(d$s), " trial units to ",
                               "generalize to ", nrow(d), "; logit SMD "))
  expect_identical(v, g)
})

test_that("gensw() refuses what it cannot weight", {
  d <- gw_dat()
  expect_error(gensw(d, "age", "prior"),
               "gensw: column `age` must hold 0 or 1.", fixed = TRUE)
  d1 <- transform(d, s = 1)
  expect_error(gensw(d1, "s", "age"),
               "gensw: column `s` must mark units in and out of the trial.",
               fixed = TRUE)
  expect_error(gensw(d, "s", "s"),
               "gensw: `s` and `x` must name different columns.", fixed = TRUE)
  expect_error(gensw(d, "s", "age", trim = 0.5),
               "gensw: `trim` must be below 0.5.", fixed = TRUE)
  sep <- data.frame(s = rep(0:1, each = 10), x = 1:20)
  expect_error(gensw(sep, "s", "x"),
               paste("gensw: the sampling model cannot be fitted: the",
                     "covariates separate the trial from the rest, or are",
                     "collinear."), fixed = TRUE)
})
