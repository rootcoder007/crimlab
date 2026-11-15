# The plot methods of the regression discontinuity (R/rdplt.R), checked
# as those of the exposure-and-crime family are in test-apxplt.R, with
# keeps() and chk() of helper-plot.R. The data are simulated and small.

test_that("the positions the discontinuity plots compute, by hand", {
  # a line through (-1, 1), (0, 2), (1, 3); with weights 1, 2, 1 on
  # 0, 1, 4: s1 = 0, so a = 6 / 4 and b = 4 / 2
  expect_equal(.rdwls(c(-1, 0, 1), c(1, 2, 3), c(1, 1, 1)), c(2, 1),
               tolerance = 1e-12)
  expect_equal(.rdwls(c(-1, 0, 1), c(0, 1, 4), c(1, 2, 1)), c(1.5, 2))
  expect_identical(.rdwls(c(1, 1), c(1, 2), c(1, 1)), c(NA_real_, NA_real_))
  # two bins of [-1, 0]: -1 and -0.9 in the first, -0.4 in the second
  b <- .rdbin(c(-1, -0.9, -0.4), c(1, 3, 5), c(-1, 0), 2L, "left")
  expect_identical(b$x, c(-0.75, -0.25))
  expect_identical(b$y, c(2, 5))
  expect_identical(b$n, c(2L, 1L))
  # an empty bin is left out
  expect_identical(.rdbin(c(-1, -0.9), c(1, 3), c(-1, 0), 2L, "left")$n, 2L)
})

test_that("rdloc keeps its bins and lines and draws them", {
  set.seed(20261120)
  d <- data.frame(r = runif(600, -1, 1))
  d$y <- 1 + 0.8 * d$r + 0.5 * (d$r >= 0) + rnorm(600, 0, 0.3)
  f <- rdloc(d, "y", "r", plac = c(-0.5, 0.5))
  # the attribute leaves the elements as they were
  expect_named(f, c("coef", "h", "n", "nl", "nr", "ncl", "ik", "placebo",
                    "cut", "kern", "fuzzy", "se", "level"))
  a <- attr(f, "plt")
  expect_identical(sum(a$bins$n), 600L)
  # the gap of the two lines at the cutoff is the conventional estimate
  ft <- a$fit
  gap <- ft$y[ft$side == "right"][1L] - ft$y[ft$side == "left"][2L]
  expect_equal(gap, f$coef$est[f$coef$method == "conventional"],
               tolerance = 1e-8)
  p <- chk(f)
  expect_named(p, c("bins", "fit", "placebo"))
  expect_identical(p$bins, a$bins)
  z <- stats::qnorm(0.975)
  expect_equal(p$placebo$lo, f$placebo$est - z * f$placebo$se,
               tolerance = 1e-14)
  g <- f
  attr(g, "plt") <- NULL
  expect_error(drw(plot(g)), "plot: `x` lacks the binned data")
})

test_that("rddens keeps its histogram and smooths and draws them", {
  set.seed(6)
  x <- rddens(runif(800, -1, 1))
  expect_named(x, c("tab", "b", "h", "nbin", "n", "cut", "level"))
  a <- attr(x, "plt")
  # the heights are shares over the binsize: they sum to 1 / b
  expect_equal(sum(a$bins$y) * x$b, 1, tolerance = 1e-12)
  # the smooths end at the densities of the fit
  ft <- a$fit
  expect_equal(utils::tail(ft$y[ft$side == "left"], 1L), x$tab$fl,
               tolerance = 1e-8)
  expect_equal(ft$y[ft$side == "right"][1L], x$tab$fr, tolerance = 1e-8)
  p <- chk(x)
  expect_named(p, c("bins", "fit"))
  expect_identical(p$fit, ft)
})
