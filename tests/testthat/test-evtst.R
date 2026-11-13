# evtst(): the binned event study against lm() and 'fixest' (values
# stored), the pre-trend test, the dropped window, the cohort trends by
# hand, and the refusals.

pd <- pepanel()

test_that("the leads and lags are fixest's, with its errors and F", {
  e <- evtst(pd$rate, pd$town, pd$year, pd$first, leads = 3, lags = 3)
  # fixest 0.14.2, dev/ref-pe.R: feols(rate ~ i(k, ref = c(-1, -1000)) |
  # town + year, cluster = ~town), k binned to [-3, 3]
  k <- e$est$evtm != -1
  expect_equal(e$est$est[k], c(0.200389386024178, 0.915200963302268,
                               -1.08276313044967, -2.72861583195067,
                               -3.67281342680414, -6.17992763628308),
               tolerance = 1e-9)
  expect_equal(e$est$se[k], c(0.383804035157557, 0.38143170242839,
                              0.347148973356146, 0.319943822552955,
                              0.377726395808488, 0.654282120428571),
               tolerance = 1e-9)
  expect_equal(e$est$est[!k], 0)
  expect_true(is.na(e$est$se[!k]))
  # fixest's wald() of the two leads: F = 4.14298537836381
  expect_equal(e$pre$f, 4.14298537836381, tolerance = 1e-9)
  expect_equal(e$pre$df1, 2)
  expect_equal(e$pre$df2, 29)
  expect_equal(e$pre$p, pf(e$pre$f, 2, 29, lower.tail = FALSE))
  expect_equal(dim(e$vcov), c(6L, 6L))
  expect_output(print(e), "leads jointly zero")
})

test_that("rows outside the window are dropped without binning", {
  e <- evtst(pd$rate, pd$town, pd$year, pd$first, leads = 2, lags = 2,
             bin = FALSE, ref = 0)
  k <- ifelse(is.na(pd$first), NA, pd$year - pd$first)
  keep <- is.na(k) | (k >= -2 & k <= 2)
  d <- pd[keep, ]
  kk <- k[keep]
  m <- lm(rate ~ I(kk %in% -2) + I(kk %in% -1) + I(kk %in% 1) +
            I(kk %in% 2) + factor(town) + factor(year), d)
  expect_equal(e$est$est[-3], unname(coef(m)[2:5]), tolerance = 1e-9)
  expect_equal(e$n, sum(keep))
})

test_that("the cohort trends are the moving means, scaled at -1", {
  e <- evtst(pd$rate, pd$town, pd$year, pd$first, leads = 2, lags = 2,
             ma = 2)
  m <- tapply(pd$rate[pd$first %in% 2020], pd$year[pd$first %in% 2020], mean)
  a <- stats::filter(m, c(0.5, 0.5), sides = 1)
  v <- a[names(m) %in% 2018:2022] / a[names(m) == "2019"]
  expect_equal(e$coh$value[e$coh$cohort == 2020], as.vector(v),
               tolerance = 1e-14)
  expect_identical(unique(e$coh$cohort), c(2018, 2020, 2022))
  e1 <- evtst(pd$rate, pd$town, pd$year, pd$first, leads = 1, lags = 1,
              norm = FALSE)
  expect_equal(e1$coh$value[1], mean(pd$rate[pd$first %in% 2018 &
                                                pd$year == 2017]))
  # leads = 1 with ref = -1: no lead to test
  expect_true(is.na(e1$pre$f))
  expect_output(print(e1), "evtst")
})

test_that("evtst() refuses what it cannot fit", {
  y <- pd$rate
  u <- pd$town
  t <- pd$year
  f <- pd$first
  expect_error(evtst(y, u, t, f, ref = 5), "`ref` must be")
  expect_error(evtst(y, u, t + 0.5, f), "whole numbers")
  expect_error(evtst(y, u, t, rep(NA_real_, 300)), "no unit is ever treated")
  expect_error(evtst(y, u, t, rep(1, 300) * t), "constant within")
  expect_error(evtst(y, u, t, -Inf * (u > 0)), "or NA, 0 or Inf")
  expect_error(evtst(y, u, t, letters[1:3]), "one per row")
  expect_error(evtst(y, u, t, f, leads = 9, bin = FALSE), "event time -9")
  expect_error(evtst(y, u, t, f, w = pd$w * 0 + (u > 3)), NA)
})
