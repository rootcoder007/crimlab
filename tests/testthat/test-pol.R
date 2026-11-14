# Values of the policing kernels are checked in C++ (src/t_pol.cpp); here
# the labels, the messages and rmorie parity (rmorie's values stored in
# ref/pol.txt by dev/ref/pol.R, version named there).

test_that("clrsmp labels its parts", {
  r <- clrsmp(c(minor = 90, serious = 20), c(100, 50), c(19, 60), c(20, 140))
  p <- r[["parts"]]
  expect_identical(p[["part"]], c("minor", "serious"))
})

test_that("policing functions stop on bad input, naming the argument", {
  errs("clrsmp", list(
    list("`reca`", c(1, 1), c(1, 0), c(1, 1), c(1, 1)),
    list("`clra`", c(1, -1), c(1, 1), c(1, 1), c(1, 1)),
    list("`clrb`", c(1, 1), c(1, 1), c(1, NA), c(1, 1)),
    list("length 2", c(1, 1), c(1, 1), c(1, 1), c(1, 1, 1)),
    list("`recb`", c(1, 1), c(1, 1), c(1, 1), "a")
  ))
  errs("expbnd", list(
    list("`y`", c(1, 0), c(1, 1), 1),
    list("`m`", c(1, 1), c(1, -1), 1),
    list("`gamma`", c(1, 1), c(1, 1), 0.9),
    list("`gamma`", c(1, 1), c(1, 1), c(1, 2)),
    list("`ref`", c(1, 1), c(1, 1), 1, ref = 3)
  ))
  errs("misrd", list(
    list("`r1`", 1.1, 0.1, 0.5, 0.9, 0.9),
    list("`r0`", 0.1, NA, 0.5, 0.9, 0.9),
    list("`prev`", 0.1, 0.1, 1, 0.9, 0.9),
    list("`se`", 0.1, 0.1, 0.5, -0.1, 0.9),
    list("`sp`", 0.1, 0.1, 0.5, 0.9, c(1, 1)),
    list("recorded group empty", 0.1, 0.1, 0.5, 0, 1),
    list("recorded group empty", 0.1, 0.1, 0.5, 1, 0)
  ))
  errs("rtmhs", list(
    list("length 2", c(1, 2), c(1, 2, 3), 0),
    list("`x1`", c(1, Inf), c(1, 2), 0),
    list("`cut`", c(1, 2), c(1, 2), NA),
    list("`w`", c(1, 2), c(1, 2), 0, w = c(1, -1)),
    list("exceeds `cut`", c(1, 2), c(1, 2), 5),
    list("exceeds `cut`", c(1, 2), c(1, 2), 1, w = c(1, 0))
  ))
})

test_that("policing results agree with rmorie", {
  ref <- refget("pol")
  y <- c(A = 300, B = 100, C = 50)
  m <- c(A = 1000, B = 1000, C = 400)
  # rmorie::morie_disparity_exposure_bounds(y, m, gamma = 1.5,
  #   reference = "B")
  o <- ref$exb
  d <- expbnd(y, m, 1.5, ref = "B")
  expect_equal(d[["rate"]], o$rate, tolerance = 1e-14)
  expect_equal(d[["lo"]], o$lo, tolerance = 1e-14)
  expect_equal(d[["hi"]], o$hi, tolerance = 1e-14)
  expect_equal(d[["rlo"]], o$rlo, tolerance = 1e-14)
  expect_equal(d[["rhi"]], o$rhi, tolerance = 1e-14)
  expect_identical(d[["sign"]] != 0, o$dir)

  # from rmorie::morie_misclassified_rd(0.3, 0.1, 0.4, 0.8, 0.9)
  d <- misrd(0.3, 0.1, 0.4, 0.8, 0.9)
  expect_equal(d[["measured"]], ref$mis$measured, tolerance = 1e-14)
  expect_equal(d[["factor"]], ref$mis$factor, tolerance = 1e-14)

  set.seed(5)
  mu <- rgamma(300, 2, 0.5)
  x1 <- rpois(300, mu)
  x2 <- rpois(300, mu)
  cut <- quantile(x1, 0.8, names = FALSE)
  # from rmorie::morie_regression_to_mean(x1, x2, cut)
  o <- ref$rtm
  expect_equal(c(sum(mu), sum(x1), sum(x2)), o$inp, tolerance = 1e-12,
               info = "inputs changed")
  d <- rtmhs(x1, x2, cut)
  expect_equal(d[["change"]], o$change, tolerance = 1e-14)
  expect_equal(d[["mirror"]], o$mirror, tolerance = 1e-14)
  expect_equal(d[["sym"]], o$sym, tolerance = 1e-14)
  expect_equal(d[["low"]], o$low, tolerance = 1e-14)
  expect_equal(d[["lowsym"]], o$lowsym, tolerance = 1e-14)
})
