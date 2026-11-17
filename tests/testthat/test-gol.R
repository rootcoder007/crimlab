# The generalised ordered logit written out: P(y > j) = plogis(a_j + x b_j),
# b_j one slope for a parallel column, maximised by optim().
golref <- function(yc, X, par) {
  M <- max(yc)
  J <- M - 1L
  k <- ncol(X)
  nb <- ifelse(par, 1L, J)
  ll <- function(u) {
    a <- u[seq_len(J)]
    b <- split(u[-seq_len(J)], rep(seq_len(k), nb))
    eta <- vapply(seq_len(J), function(j) {
      a[j] + drop(X %*% vapply(b, function(s) s[min(j, length(s))], 0))
    }, numeric(nrow(X)))
    up <- cbind(1, stats::plogis(eta), 0)
    pr <- up[cbind(seq_along(yc), yc)] - up[cbind(seq_along(yc), yc + 1L)]
    if (any(pr <= 0)) return(-1e10)
    sum(log(pr))
  }
  st <- c(stats::qlogis(1 - cumsum(tabulate(yc, M))[-M] / length(yc)),
          rep(0, sum(nb)))
  o <- stats::optim(st, ll, method = "BFGS", hessian = TRUE,
                    control = list(fnscale = -1, reltol = 1e-14,
                                   maxit = 5000))
  list(est = o$par, ll = o$value, se = sqrt(diag(solve(-o$hessian))))
}

sim <- function(n, seed) {
  set.seed(seed)
  d <- data.frame(race = rbinom(n, 1, 0.4), age = rnorm(n))
  u <- runif(n)
  d$y <- 1 + (u < plogis(0.6 + 0.9 * d$race - 0.5 * d$age)) +
    (u < plogis(-0.8 + 0.1 * d$race - 0.5 * d$age))
  d
}

test_that("gologt() maximises the likelihood written out", {
  d <- sim(500, 61)
  X <- cbind(d$race, d$age)
  for (par in list(NULL, "age", c("race", "age"))) {
    r <- gologt(d, "y", c("race", "age"), par = par)
    w <- golref(d$y, X, c("race", "age") %in% par)
    expect_equal(r$coef$est, w$est, tolerance = 1e-5)
    expect_equal(r$fit[["loglik"]], w$ll, tolerance = 1e-9)
    expect_equal(r$coef$se, w$se, tolerance = 1e-3)
    expect_identical(r$fit[["conv"]], 1)
  }
  r <- gologt(d, "y", c("race", "age"), par = "age")
  expect_identical(r$coef$term,
                   c("(Intercept)", "(Intercept)", "race", "race", "age"))
  expect_identical(r$coef$thr, c("1|2", "2|3", "1|2", "2|3", "all"))
  expect_identical(r$lev, c("1", "2", "3"))
  z <- stats::qnorm(0.975)
  expect_equal(r$coef$lo, r$coef$est - z * r$coef$se, tolerance = 1e-14)
  expect_equal(r$coef$p, 2 * stats::pnorm(-abs(r$coef$est / r$coef$se)),
               tolerance = 1e-12)
  # the likelihood-ratio test against proportional odds
  po <- gologt(d, "y", c("race", "age"), par = c("race", "age"))
  lr <- 2 * (r$fit[["loglik"]] - po$fit[["loglik"]])
  expect_equal(r$lr[["stat"]], lr, tolerance = 1e-12)
  expect_identical(r$lr[["df"]], 1)
  expect_equal(r$lr[["p"]], stats::pchisq(lr, 1, lower.tail = FALSE),
               tolerance = 1e-12)
  expect_true(is.na(po$lr[["stat"]]))
  # the Wald test that race has one slope
  b <- r$coef$est[3:4]
  v <- r$vcov[3:4, 3:4]
  wd <- (b[1] - b[2])^2 / (v[1, 1] + v[2, 2] - 2 * v[1, 2])
  expect_equal(r$wald$stat[1L], wd, tolerance = 1e-10)
  expect_true(is.na(r$wald$stat[2L]))
  expect_identical(r$fit[["aic"]], -2 * r$fit[["loglik"]] + 10)
})

test_that("gologt() reads the order of the levels", {
  d <- sim(300, 62)
  a <- gologt(d, "y", "race")
  d$f <- factor(c("lo", "mid", "hi")[d$y], c("lo", "mid", "hi"))
  b <- gologt(d, "f", "race")
  expect_identical(b$coef$est, a$coef$est)
  expect_identical(b$coef$thr[1:2], c("lo|mid", "mid|hi"))
  d$s <- as.character(d$f)
  expect_identical(
    gologt(d, "s", "race", lev = c("lo", "mid", "hi"))$coef$est, a$coef$est
  )
})

test_that("gologt() refuses bad input in words", {
  d <- sim(60, 63)
  d$s <- c("a", "b", "c")[d$y]
  e <- function(expr, m) expect_error(expr, m, fixed = TRUE)
  e(gologt(d, "s", "race"),
    "gologt: `lev` must give the order of the levels of column `s`, text.")
  e(gologt(d, "s", "race", lev = c("a", "b")),
    "gologt: column `s` has values not in its levels: c.")
  e(gologt(d, "s", "race", lev = c("a", "b", "c", "d")),
    "gologt: column `s` has no rows at level(s) d.")
  e(gologt(transform(d, y = pmin(y, 2)), "y", "race"),
    "gologt: column `y` must have at least three levels.")
  e(gologt(d, "y", "race", par = "age"),
    "gologt: `par` must name columns of `x`, each once.")
  e(gologt(d, "y", "race", level = 1),
    "gologt: `level` must be strictly between 0 and 1.")
  e(gologt(d, "y", "race", maxit = 0),
    "gologt: `maxit` must be between 1 and 2147483647.")
  e(gologt(d, "y", "zz"), "gologt: `x` names columns not in `data`: zz.")
})

test_that("gocoef() lays the coefficients out by threshold", {
  set.seed(64)
  n <- 300
  d <- data.frame(race = rbinom(n, 1, 0.4), age = rnorm(n))
  u <- runif(n)
  d$y <- 1 + (u < plogis(0.5 + 0.8 * d$race)) +
    (u < plogis(-1 + 0.1 * d$race - 0.3 * d$age))
  r <- gologt(d, "y", c("race", "age"), par = "age")
  g <- gocoef(r)
  cf <- r$coef
  expect_identical(dimnames(g$est),
                   list(c("1|2", "2|3"), c("(Intercept)", "race", "age")))
  expect_identical(g$est[, "race"], setNames(cf$est[3:4], c("1|2", "2|3")))
  expect_identical(unname(g$est[, "age"]), rep(cf$est[5L], 2))
  expect_identical(unname(g$se[, "(Intercept)"]), cf$se[1:2])
  expect_identical(g$wald, r$wald)
  h <- gocoef(r, "race")
  expect_identical(colnames(h$est), "race")
  expect_identical(h$wald$term, "race")
  expect_error(gocoef(r, "sex"),
               paste("gocoef: `v` must name terms of the fit: (Intercept),",
                     "race, age."),
               fixed = TRUE)
  expect_error(gocoef(list()), "gocoef: `object` must be a result of gologt().",
               fixed = TRUE)
})
