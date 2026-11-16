# dstfit(): the fits of the reference implementation (helper-dfref.R). The
# closed forms agree to rounding; the fits by the Nelder-Mead and BFGS
# minimisers follow optim()'s algorithms but not R's own density code, so
# they agree to the minimisers' tolerance (1e-6 relative); a t whose degrees
# of freedom run off to infinity (a flat likelihood) is compared on its
# location and scale only.

df_check <- function(r, ref, tol = 1e-6, flat = character(0)) {
  expect_s3_class(r, "dstfit")
  expect_named(
    r,
    c("table", "fits", "boot", "hill", "x", "nmiss", "discrete", "conflvl")
  )
  t <- r$table
  expect_named(
    t,
    c("rank", "family", "npar", "loglik", "aic", "bic", "ks", "ad",
      "converged", "params")
  )
  rt <- ref[[1L]]
  expect_identical(t$family, rt[[2L]])
  expect_equal(t$npar, rt[[3L]], tolerance = testthat_tolerance())
  ok <- !t$family %in% flat
  for (j in 4:8) {
    expect_equal(t[[j]][ok], rt[[j]][ok], tolerance = tol, info = names(t)[j])
  }
  expect_identical(t$params[ok], rt[[10L]][ok])
  for (nm in setdiff(names(r$fits), flat)) {
    expect_equal(unname(r$fits[[nm]]), ref[[2L]][[match(nm, t$family)]],
                 tolerance = tol, info = nm)
  }
  b <- r$boot
  rb <- ref[[3L]]
  expect_identical(b$family, rb[[1L]])
  expect_identical(b$param, rb[[2L]])
  keep <- !b$family %in% flat
  for (j in 4:5) {
    expect_equal(
      as.numeric(b[[j]][keep]), as.numeric(rb[[j]][keep]),
      tolerance = 1e-5
    )
  }
  if (is.null(ref[[4L]])) {
    expect_null(r$hill)
  } else {
    expect_equal(unname(unlist(r$hill)), unlist(ref[[4L]]), tolerance = 1e-12)
  }
  expect_equal(r$nmiss, ref[[5L]], tolerance = testthat_tolerance())
  expect_identical(r$discrete, ref[[6L]])
}

test_that("continuous data give the reference fits", {
  set.seed(1)
  g <- round(rgamma(60, shape = 2, rate = 0.5), 3)
  df_check(dstfit(g, nboot = 15), dfref$g, flat = "t")
  df_check(dstfit(g, nboot = 0, fams = c("pareto", "exp", "logn"), k = 12,
                  conflvl = 0.9), dfref$k, tol = 1e-12)
  set.seed(3)
  pr <- round(rbeta(40, 2, 5), 4)
  df_check(dstfit(pr, nboot = 0), dfref$b, flat = "t")
  set.seed(4)
  tt <- round(rt(45, 3) * 2 + 1, 3)
  df_check(dstfit(tt, nboot = 0, fams = c("normal", "t")), dfref$t)
})

test_that("counts give the reference fits", {
  set.seed(2)
  cn <- ifelse(runif(50) < 0.3, 0, rpois(50, 4))
  cn[c(3, 9)] <- NA
  set.seed(11)
  s <- .Random.seed
  r <- dstfit(cn, nboot = 10, seed = 7)
  expect_identical(.Random.seed, s)
  df_check(r, dfref$c, flat = c("t", "zinb"))
  expect_equal(r$fits$negbin, c(size = 0.90630779612651746, mu = 2.5),
               tolerance = 1e-8)
  expect_true(is.na(r$table$ad[r$table$family == "poisson"]))
})

test_that("dstfit() refuses what it cannot fit", {
  expect_error(dstfit(c(-1, 2, 3, 4, 5), fams = "gamma"),
               "^dstfit: `gamma` needs every value above 0\\.$")
  expect_error(dstfit(c(1, 2, 3)),
               "^dstfit: `x` needs at least five values that are not missing")
  expect_error(dstfit(rep(2, 6)), "every value of `x` is the same")
  expect_error(dstfit(c(1:5, Inf)), "^dstfit: `x` has an infinite value\\.$")
  expect_error(dstfit("a"), "^dstfit: `x` must be a numeric vector\\.$")
  expect_error(dstfit(1:10, fams = "cauchy"), "`fams` must be one or more of")
  expect_error(dstfit(1:10, nboot = -1), "`nboot` must be one whole number")
  expect_error(dstfit(1:20, nboot = 0, k = 30),
               "^dstfit: `k` must be one whole number from 2 to 19\\.$")
  expect_output(print(dstfit(1:10 + 0.5, nboot = 0)), "maximum likelihood")
})

test_that("dstfit() refuses a confidence level outside (0, 1)", {
  expect_error(
    dstfit(c(1.2, 2.3, 3.1, 4.8, 2.2), conflvl = 1),
    "dstfit: `conflvl` must be a single number strictly between 0 and 1.",
    fixed = TRUE
  )
})

test_that("the discrete fits of huge whole numbers return", {
  r <- dstfit(c(3, 1, 4, 1, 5, 9, 2, 6) * 1e50, nboot = 0)
  expect_s3_class(r, "dstfit")
  k <- r$table$ks[r$table$family == "pois"]
  expect_true(length(k) == 0L || is.na(k) || (k >= 0 && k <= 1))
})
