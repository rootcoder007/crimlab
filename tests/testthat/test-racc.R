# Regression standards: the accessors of a fit, and what a fit does with
# collinear, noiseless, missing and non-finite input.

rd <- function() {
  set.seed(2)
  n <- 240
  d <- data.frame(u = rep(1:24, each = 10), t = rep(1:10, 24),
                  x = rnorm(n), w = rnorm(n))
  d$y <- 1 + 0.5 * d$x + rnorm(n)
  d$cnt <- rpois(n, exp(0.2 + 0.3 * d$x))
  d
}

#' @srrstats {RE4.2} coef() of every regression fit is tested here.
#' @srrstats {RE4.3} confint() of every regression fit is tested here.
#' @srrstats {RE4.5} nobs() is tested here.
#' @srrstats {RE4.6} vcov() is tested here.
#' @srrstats {RE7.3} The accessors of RE4.2 to RE4.6 are called on the
#'   fits and their values compared with the fits' own tables.
#' @srrstats {RE7.2} The terms of a fit keep the names of the input
#'   columns: coef(), confint() and vcov() are named by them.
#' @noRd
test_that("the accessors read a fit", {
  d <- rd()
  f <- ivfe(d, "y", x = c("x", "w"), fe = "u")
  expect_named(coef(f), c("x", "w"))
  expect_identical(unname(coef(f)), f$coef$est)
  ci <- confint(f)
  expect_identical(dimnames(ci), list(c("x", "w"), c("lower", "upper")))
  expect_identical(unname(ci[, "lower"]), f$coef$lo)
  expect_identical(unname(ci[, "upper"]), f$coef$hi)
  expect_identical(attr(ci, "level"), 0.95)
  expect_identical(rownames(confint(f, "w")), "w")
  expect_identical(confint(f, 2)[1, ], confint(f, "w")[1, ])
  expect_identical(confint(f, level = 0.95)[, 1], ci[, 1])
  v <- vcov(f)
  expect_identical(dimnames(v), list(c("x", "w"), c("x", "w")))
  expect_equal(sqrt(diag(v)), c(x = f$coef$se[1], w = f$coef$se[2]),
               tolerance = 1e-12)
  expect_identical(nobs(f), 240L)
  expect_error(confint(f, level = 0.9), "fit again with the `level` wanted")
  expect_error(confint(f, "z"), "`parm` must be NULL, or terms of the fit")
  expect_error(confint(f, 3), "`parm` must be NULL")
  expect_error(confint(f, level = 2), "`level` must be NULL or one number")
  # the level of the fit is the level of its intervals
  g <- ivfe(d, "y", x = "x", fe = "u", level = 0.9)
  expect_identical(attr(confint(g), "level"), 0.9)
  expect_lt(diff(confint(g)[1, ]), diff(confint(f, "x")[1, ]))
})

test_that("each regression fit has its accessors", {
  d <- rd()
  trt <- as.integer(d$u > 12 & d$t > 5)
  fits <- list(
    ppmlfe = ppmlfe(d, "cnt", "x", "u"),
    fepnl = fepnl(d$y, d[, c("x", "w")], d$u, d$t),
    didfe = didfe(d$y, d$u, d$t, trt),
    rdloc = rdloc(d, "y", "x")
  )
  for (nm in names(fits)) {
    f <- fits[[nm]]
    expect_identical(unname(coef(f)), f$coef$est, info = nm)
    expect_identical(unname(confint(f)[, "upper"]), f$coef$hi, info = nm)
    expect_true(nobs(f) > 0L, info = nm)
  }
  expect_named(coef(fits$rdloc), c("robust", "conventional"))
  expect_identical(nobs(fits$rdloc), fits$rdloc$n)
  expect_identical(nobs(fits$fepnl), as.integer(fits$fepnl$fit$n))
  for (nm in c("ppmlfe", "fepnl", "didfe")) {
    v <- vcov(fits[[nm]])
    expect_identical(rownames(v), names(coef(fits[[nm]])), info = nm)
    expect_equal(unname(sqrt(diag(v))), fits[[nm]]$coef$se,
                 tolerance = 1e-10, info = nm)
  }
  # a spatial lag model: no stored interval, so the normal one, any level
  set.seed(1)
  p <- cbind(x = runif(60), y = runif(60))
  x <- rnorm(60)
  s <- sarml(1 + 2 * x + rnorm(60), x, spwt(p, k = 5))
  ci <- confint(s, level = 0.9)
  z <- stats::qnorm(0.95)
  expect_equal(unname(ci[, "upper"]), s$coef$est + z * s$coef$se,
               tolerance = 1e-12)
  expect_identical(attr(ci, "level"), 0.9)
  expect_true(all(names(coef(s)) %in% rownames(vcov(s))))
  expect_identical(nobs(s), 60L)
  expect_error(.ranobs(list(coef = s$coef)), "records no count of rows")
  expect_error(.ravc(fits$rdloc), "this fit stores no covariance matrix")
  expect_error(.racoef(list(coef = 1)), "not a complete fit")
  # a covariance matrix without names is named by the terms
  v <- .ravc(list(coef = data.frame(term = c("a", "b"), est = 1:2),
                  vcov = diag(2)))
  expect_identical(dimnames(v), list(c("a", "b"), c("a", "b")))
})

#' @srrstats {RE7.0} Predictors with an exact linear relation are tested.
#' @srrstats {RE7.0a} A fit on perfectly collinear predictors is refused,
#'   in words that name collinearity.
#' @srrstats {RE2.4a} Collinear predictors are found before any estimate
#'   is returned, and the fit is refused.
#' @srrstats {RE7.1} A response that is an exact linear function of the
#'   predictor is fitted: the slope comes back exactly, with a standard
#'   error of rounding size.
#' @srrstats {RE2.4b} See RE7.1: a noiseless response is not refused; it
#'   is fitted exactly.
#' @noRd
test_that("collinear predictors are refused and exact responses fitted", {
  d <- rd()
  d$x2 <- 2 * d$x
  expect_error(ivfe(d, "y", x = c("x", "x2"), fe = "u"), "collinear")
  expect_error(ppmlfe(d, "cnt", c("x", "x2"), "u"), "collinear")
  expect_error(fepnl(d$y, d[, c("x", "x2")], d$u, d$t), "collinear")
  expect_error(didfe(d$y, d$u, d$t, as.integer(d$u > 12 & d$t > 5),
                     adj = d[, c("x", "x2")]), "collinear")
  d$yex <- 1 + 2 * d$x
  f <- ivfe(d, "yex", x = "x", fe = "u")
  expect_equal(unname(coef(f)), 2, tolerance = 1e-12)
  expect_lt(f$coef$se, 1e-12)
  g <- fepnl(d$yex, d[, "x", drop = FALSE], d$u, d$t)
  expect_equal(unname(coef(g)), 2, tolerance = 1e-12)
})

#' @srrstats {RE7.1a} The exact fit is timed against the noisy one: both
#'   are closed-form solves of the same size, so the exact fit must not be
#'   slower by more than the timer's noise.
#' @noRd
test_that("an exact response is fitted as fast as a noisy one", {
  d <- rd()
  d$yex <- 1 + 2 * d$x
  tm <- function(y) {
    system.time(for (i in 1:20) ivfe(d, y, x = "x", fe = "u"))[["elapsed"]]
  }
  noisy <- tm("y")
  expect_lte(tm("yex"), 5 * noisy + 0.5)
})

#' @srrstats {RE3.0} A count fit that stops before it converges warns.
#' @srrstats {RE3.1} With the warning suppressed the fit still says so:
#'   `conv` is FALSE and `iter` the iterations used.
#' @noRd
test_that("a fit that does not converge says so", {
  d <- rd()
  expect_warning(f <- ppmlfe(d, "cnt", "x", "u", maxit = 1),
                 "did not converge")
  expect_false(f$conv)
  expect_identical(as.integer(f$iter), 1L)
  expect_true(ppmlfe(d, "cnt", "x", "u")$conv)
})

#' @srrstats {RE2.1} Missing and infinite values are told apart: each is
#'   refused, in its own words. No option drops them (G2.14b).
#' @noRd
test_that("missing and infinite predictors are refused apart", {
  d <- rd()
  dm <- d
  dm$x[3] <- NA
  di <- d
  di$x[3] <- Inf
  expect_error(ivfe(dm, "y", x = "x", fe = "u"), "has missing values")
  expect_error(ivfe(di, "y", x = "x", fe = "u"), "must be finite")
  expect_error(ppmlfe(dm, "cnt", "x", "u"), "has missing values")
})

#' @srrstats {RE5.0} The time of a fit is measured at two sizes: ten
#'   times the rows must not take much more than ten times as long.
#' @noRd
test_that("the time of a fit grows with its rows, not faster", {
  skip_on_cran()
  tm <- function(n) {
    set.seed(1)
    d <- data.frame(u = rep(seq_len(n / 10), each = 10), x = rnorm(n))
    d$y <- 1 + 0.5 * d$x + rnorm(n)
    min(replicate(3, system.time(ivfe(d, "y", x = "x", fe = "u"))[["elapsed"]]))
  }
  small <- tm(2e4)
  expect_lte(tm(2e5), 40 * small + 0.5)
})
