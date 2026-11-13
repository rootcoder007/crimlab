# cntreg(): Poisson and negative binomial fits against stats::glm() and
# MASS::glm.nb() (the same iterations, so agreement to rounding), the
# cluster sandwich against sandwich::vcovCL(), the exchangeable GEE against
# geepack, closed forms, adversarial counts, the input forms and the
# refusals.

set.seed(3)
cn <- local({
  id <- rep(1:8, each = 10)
  x <- rbinom(80, 1, 0.5)
  f <- factor(sample(c("a", "b", "c"), 80, TRUE))
  u <- rnorm(8, 0, 0.4)[id]
  y <- rnbinom(80, mu = exp(1 + 0.3 * x + u), size = 3)
  e <- runif(80, 1, 3)
  data.frame(y, x, f, id, e)
})

pick <- function(r, fit, what = "est") r$coef[[what]][r$coef$fit == fit]

test_that("the Poisson fit is glm()'s", {
  r <- cntreg(cn$y, cn$x, adj = cn["f"], expo = cn$e)
  g <- glm(y ~ factor(x) + f + offset(log(e)), family = poisson, data = cn)
  co <- unname(coef(summary(g)))
  expect_equal(pick(r, "pois"), co[, 1], tolerance = 1e-10)
  expect_equal(pick(r, "pois", "se"), co[, 2], tolerance = 1e-10)
  expect_equal(pick(r, "pois", "p"), co[, 4], tolerance = 1e-8)
  expect_equal(r$fit$aic[1], AIC(g), tolerance = 1e-10)
  ci <- exp(co[2, 1] + c(-1, 1) * qnorm(0.975) * co[2, 2])
  expect_equal(c(r$irr$lo[1], r$irr$hi[1]), ci, tolerance = 1e-10)
  expect_identical(r$irr$term, c("1", "1"))
  expect_identical(r$coef$term[1:4], c("intercept", "1", "f=b", "f=c"))
})

test_that("the negative binomial fit is MASS::glm.nb()'s", {
  skip_if_not_installed("MASS")
  r <- cntreg(cn$y, cn$x, adj = cn["f"], expo = cn$e)
  nb <- MASS::glm.nb(y ~ factor(x) + f + offset(log(e)), data = cn)
  co <- unname(coef(summary(nb)))
  expect_equal(pick(r, "nb"), co[, 1], tolerance = 1e-8)
  expect_equal(pick(r, "nb", "se"), co[, 2], tolerance = 1e-8)
  expect_equal(r$fit$theta[2], nb$theta, tolerance = 1e-8)
  expect_equal(r$fit$setheta[2], nb$SE.theta, tolerance = 1e-6)
  expect_equal(r$fit$aic[2], AIC(nb), tolerance = 1e-8)
  expect_equal(r$fit$loglik[2], as.numeric(logLik(nb)), tolerance = 1e-8)
})

test_that("the cluster sandwich is sandwich::vcovCL()'s", {
  skip_if_not_installed("sandwich")
  skip_if_not_installed("MASS")
  r <- cntreg(cn$y, cn$x, adj = cn["f"], expo = cn$e, clus = cn$id)
  g <- glm(y ~ factor(x) + f + offset(log(e)), family = poisson, data = cn)
  nb <- MASS::glm.nb(y ~ factor(x) + f + offset(log(e)), data = cn)
  expect_equal(
    pick(r, "poiscl", "se"),
    unname(sqrt(diag(sandwich::vcovCL(g, cluster = cn$id)))),
    tolerance = 1e-9
  )
  expect_equal(
    pick(r, "nbcl", "se"),
    unname(sqrt(diag(sandwich::vcovCL(nb, cluster = cn$id)))),
    tolerance = 1e-7
  )
  # cadj = FALSE drops G / (G - 1): the independence GEE variance
  r0 <- cntreg(
    cn$y,
    cn$x,
    adj = cn["f"],
    expo = cn$e,
    clus = cn$id,
    cadj = FALSE
  )
  expect_equal(
    pick(r0, "poiscl", "se"),
    pick(r, "poiscl", "se") * sqrt(7 / 8),
    tolerance = 1e-12
  )
})

test_that("the exchangeable GEE is geepack's", {
  skip_if_not_installed("geepack")
  r <- cntreg(cn$y, cn$x, adj = cn["f"], expo = cn$e, clus = cn$id, gee = TRUE)
  g <- geepack::geeglm(
    y ~ factor(x) + f + offset(log(e)),
    family = poisson,
    id = id,
    data = cn,
    corstr = "exchangeable",
    control = geepack::geese.control(epsilon = 1e-12, maxit = 100)
  )
  expect_equal(pick(r, "gee"), unname(coef(g)), tolerance = 1e-8)
  expect_equal(
    pick(r, "gee", "se"),
    unname(sqrt(diag(g$geese$vbeta))),
    tolerance = 1e-7
  )
  expect_equal(r$fit$alpha[3], unname(g$geese$alpha), tolerance = 1e-7)
  expect_true(r$fit$conv[3])
})

test_that("two groups: the ratio of the pooled rates, by hand", {
  r <- cntreg(c(3, 5, 10, 2), c("a", "a", "b", "b"), expo = c(2, 3, 4, 1))
  # to glm()'s stopping rule (the deviance settled to 1e-8): its estimate
  # and its variance, taken at the weights of the last step
  expect_equal(r$irr$irr[1], (12 / 5) / (8 / 5), tolerance = 1e-9)
  se <- sqrt(1 / 8 + 1 / 12)
  expect_equal(r$irr$lo[1], 1.5 * exp(-qnorm(0.975) * se), tolerance = 1e-5)
  # the reference swapped inverts the ratio
  r2 <- cntreg(
    c(3, 5, 10, 2),
    c("a", "a", "b", "b"),
    expo = c(2, 3, 4, 1),
    ref = "b"
  )
  expect_equal(r2$irr$irr[1], 1 / 1.5, tolerance = 1e-9)
  expect_identical(r2$ref, "b")
  expect_identical(r2$irr$term[1], "a")
})

test_that("adversarial counts", {
  # separation: the reference group has no counts
  r <- suppressWarnings(cntreg(c(0, 0, 4, 6), c(1, 1, 2, 2)))
  expect_true(r$irr$sep[1])
  expect_gt(r$irr$irr[1], 1e6)
  expect_output(print(r), "not estimable")
  # Poisson data: theta runs large (MASS::glm.nb() stops at its iteration
  # limit too) and the two fits agree
  set.seed(1)
  y <- rpois(400, 5)
  g <- rep(1:2, 200)
  r <- cntreg(y, g)
  expect_equal(r$irr$irr[2], r$irr$irr[1], tolerance = 1e-3)
  expect_gt(r$fit$theta[2], 50)
  expect_false(r$fit$conv[2])
  # huge counts
  r <- cntreg(c(2e9, 3e9, 1e9, 4e9), c(0, 1, 0, 1))
  expect_equal(r$irr$irr[1], 7 / 3, tolerance = 1e-10)
  # a column the group determines is aliased, the ratio kept
  r <- cntreg(cn$y, cn$x, adj = data.frame(xx = cn$x * 2, f = cn$f))
  expect_true(is.na(pick(r, "pois")[3]))
  expect_identical(r$fit$rank[1], 4)
  # counts beyond the double range: the fit leaves the finite range
  expect_error(cntreg(c(1.7e308, 1e308, 1e-300, 5), c(1, 2, 1, 2)),
               "the fit failed")
  # counts equal to their group means: no finite theta, the negative
  # binomial rows NA
  expect_warning(r <- cntreg(c(4, 4, 6, 6), c(1, 1, 2, 2)), "no finite")
  expect_true(is.na(r$irr$irr[2]))
  # a named adjustment matrix names its terms
  r <- cntreg(cn$y, cn$x, adj = cbind(w = cn$e))
  expect_identical(r$coef$term[3], "w")
  # whole numbers are not required
  r <- cntreg(c(1.5, 2.5, 4, 6.5), c(0, 0, 1, 1))
  expect_equal(r$irr$irr[1], 10.5 / 4, tolerance = 1e-10)
})

test_that("the input forms agree", {
  a <- cntreg(cn$y, cn$f, adj = cn$x, expo = cn$e)
  b <- cntreg(
    cn$y,
    as.character(cn$f),
    adj = matrix(cn$x, ncol = 1),
    expo = cn$e
  )
  expect_equal(a$irr$irr, b$irr$irr, tolerance = 1e-14)
  expect_identical(b$coef$term[4], "adj1")
  # named vectors are matched by name
  y <- stats::setNames(cn$y, seq_len(80))
  gg <- stats::setNames(cn$f, seq_len(80))
  o <- 80:1
  expect_equal(
    cntreg(y[o], gg)$irr$irr,
    cntreg(cn$y, cn$f)$irr$irr,
    tolerance = 1e-12
  )
  # logical adjustment columns are labels
  d <- cntreg(cn$y, cn$f, adj = list(big = cn$x == 1))
  expect_identical(d$coef$term[4], "big=TRUE")
  # a clpre list fills the arguments by name; arguments given override it
  p <- structure(list(y = cn$y, grp = cn$f, expo = cn$e), class = "clpre")
  expect_identical(cntreg(p)$irr, cntreg(cn$y, cn$f, expo = cn$e)$irr)
  expect_identical(cntreg(p, ref = "b")$ref, "b")
  expect_output(print(cntreg(cn$y, cn$x, clus = cn$id)), "8 clusters")
})

test_that("refusals", {
  errs(
    "cntreg",
    list(
      list("`y` has no events", c(0, 0), c(1, 2)),
      list("`y` must be at least 0", c(-1, 2), c(1, 2)),
      list("`y` must be finite", c(NA, 2), c(1, 2)),
      list("`grp` must give the group", c(1, 2)),
      list("`grp` has one group", c(1, 2), c(1, 1)),
      list("`grp` must be labels", c(1, 2), c(1, NA)),
      list("`expo` must be greater than 0", c(1, 2), c(1, 2), NULL, c(0, 1)),
      list("`clus` has one cluster", c(1, 2), c(1, 2), NULL, NULL, c(1, 1)),
      list("`gee = TRUE` needs `clus`", c(1, 2), c(1, 2), gee = TRUE),
      list("`level` must be one number", c(1, 2), c(1, 2), level = 1),
      list("`ref` names no group", c(1, 2), c(1, 2), ref = "z"),
      list("`adj` column adj1 must be a vector", c(1, 2), c(1, 2), list(1:3)),
      list(
        "has value labels",
        c(1, 2),
        c(1, 2),
        list(structure(c(1, 2), labels = c(a = 1)))
      ),
      list("`adj` must have one row", c(1, 2), c(1, 2), matrix(1, 3, 2)),
      list("`adj` must hold finite", c(1, 2), c(1, 2), matrix(NA_real_, 2, 2))
    )
  )
})
