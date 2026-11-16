# exprspplr(): the values of the reference implementation (rmoriebricklayer
# 1.1.1 with its default estimator plrcf(), computed once on these simulated
# inputs and printed with 17 significant digits), the custom estimator, R's
# stream left as it was, and the refusals.

ex_data <- function() {
  set.seed(1)
  d <- data.frame(x1 = round(rnorm(40), 3), x2 = round(rnorm(40), 3))
  d$exposure <- round(20 + 2 * d$x1 + rnorm(40), 3)
  d$asthma <- round(5 + 0.3 * d$exposure + d$x2 + rnorm(40), 3)
  d
}
ex_nm <- c("ate", "seana", "seboot", "cilobs", "ciupbs", "nboot", "method")
plr_m <- paste(
  "partially linear regression, DML2 cross-fitting, OLS nuisances +",
  "percentile bootstrap"
)

test_that("the built-in estimator gives the reference values", {
  d <- ex_data()
  expect_identical(d$asthma[c(1, 40)],
                   c(9.7829999999999995, 12.683999999999999))
  ref <- list(
    list(list(d, "asthma", "exposure", c("x1", "x2"), nboot = 20),
         c(0.44513547578785756, 0.14184699663664993, 0.16610723972828578,
           0.273793985819547, 0.81523421115678296), 20L),
    list(local({
      d2 <- d
      d2$x2[c(3, 7)] <- NA
      d2$asthma[11] <- NA
      list(d2, "asthma", "exposure", c("x1", "x2"), nboot = 15, rstate = 7)
    }),
    c(0.53391097611688754, 0.15059007181552186, 0.19265165227437603,
      0.33319233895369105, 0.94033052304400866), 15L),
    list(list(d, "asthma", "exposure", "x1", nboot = 12, rstate = 2^40 + 5),
         c(0.58332257873519111, 0.15469125985579976, 0.17853927683154383,
           0.30011209404619626, 0.83967914602003324), 12L),
    # collinear confounders: the least squares fit by QR
    list(local({
      d3 <- d
      d3$x3 <- 2 * d3$x1
      list(d3, "asthma", "exposure", c("x1", "x3", "x2"), nboot = 10,
           rstate = 3)
    }),
    c(0.44513547578785562, 0.14184699663665004, 0.1559221415359725,
      0.29864825518678056, 0.77007051448994146), 10L)
  )
  for (r in ref) {
    e <- do.call(exprspplr, r[[1L]])
    expect_s3_class(e, "exprspplr")
    expect_named(e, ex_nm)
    expect_equal(unlist(e[1:5]), stats::setNames(r[[2L]], ex_nm[1:5]),
                 tolerance = 1e-10)
    expect_identical(e$nboot, r[[3L]])
    expect_identical(e$method, plr_m)
  }
})

test_that("a custom estimator is resampled as in the reference", {
  d <- ex_data()
  ols <- function(data, outcome, treatment, covariates) {
    f <- stats::reformulate(c(treatment, covariates), outcome)
    cf <- summary(stats::lm(f, data = data))$coefficients
    list(ate = cf[treatment, 1], se = cf[treatment, 2], method = "OLS")
  }
  e <- exprspplr(d, "asthma", "exposure", c("x1", "x2"), nboot = 11,
                 estr = ols)
  expect_equal(
    unlist(e[1:5]),
    stats::setNames(c(0.51406800897354532, 0.17999475126763784,
                      0.1262353154591874, 0.35925705674114017,
                      0.70019079788589778), ex_nm[1:5]),
    tolerance = 1e-10
  )
  expect_identical(e$method, "OLS + percentile bootstrap")
  expect_output(print(e), "percentile interval")
})

test_that("R's random stream is left as it was", {
  d <- ex_data()
  set.seed(99)
  s <- .Random.seed
  exprspplr(d, "asthma", "exposure", "x1", nboot = 10)
  expect_identical(.Random.seed, s)
})

test_that("exprspplr() refuses what it cannot use", {
  d <- ex_data()
  expect_error(
    exprspplr(d, "asthma", "exposure", "x1", nboot = 5),
    paste0("^exprspplr: Only 5 successful bootstrap fits; need >=10\\. ",
           "Inspect the model's convergence on your resamples\\.$")
  )
  expect_error(exprspplr(d, "asthma", "exposure", "zz", nboot = 10),
               "^exprspplr: data is missing columns: zz\\.$")
  expect_error(
    exprspplr(d, "asthma", "exposure", "x1", nboot = 2.5),
    "^exprspplr: `nboot` must be a whole number of resamples, at most"
  )
  expect_error(
    exprspplr(d, "asthma", "exposure", "x1", rstate = NA),
    "^exprspplr: `rstate` must be a single whole number from -2147483647"
  )
  expect_error(exprspplr(d[1:8, ], "asthma", "exposure", "x1", nboot = 10),
               "need at least 10 complete rows for 5 folds; have 8\\.$")
  expect_error(exprspplr(d, "asthma", "x1", "x1", nboot = 10),
               "the treatment has no variation left")
  expect_error(exprspplr(d, 1, "exposure", "x1"),
               "^exprspplr: `y` must name one column of `data`\\.$")
  expect_error(exprspplr(d, "asthma", "exposure", "x1", estr = 3),
               "`estr` must be a function")
  expect_error(exprspplr(d, "asthma", "exposure", "x1", estr = "no.such.fn"),
               "`estr` names no function")
  d$x1 <- as.character(d$x1)
  expect_error(exprspplr(d, "asthma", "exposure", "x1"),
               "column `x1` must be numeric")
})

test_that("exprspplr() reports an estimator that fails on the full sample", {
  d <- ex_data()
  bad <- function(data, outcome, treatment, covariates) stop("no model.")
  expect_error(
    exprspplr(d, "asthma", "exposure", "x1", nboot = 10, estr = bad),
    paste0("^exprspplr: `estr` failed on the full sample, called as ",
           "estr\\(data, y, expo, confs\\): no model\\.$")
  )
})

test_that("exprspplr() refuses a list for a data frame in words", {
  expect_error(exprspplr(list(y = 1:3, e = 1:3), "y", "e", NULL),
               "^exprspplr: `data` must be a data frame\\.$")
})
