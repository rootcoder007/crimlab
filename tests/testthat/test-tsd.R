# Counts over time and across categories (src/tsd.cpp) against the values
# of the rmorie functions they port. The ARIMA fits take the steps of the
# BFGS optimiser of optim(); the ARIMA(1,1,1) fit equals the reference to
# 1e-13 and is compared at 1e-10. The seasonal fit differs from it by
# 1e-7 in its coefficients, 6e-10 (relative) in its forecasts and 2e-8 in
# the errors of those forecasts, so it is compared at 1e-8 and its errors
# at 1e-7.

test_that("arimfc() equals the reference forecast", {
  y <- c(52, 48, 61, 55, 47, 50, 66, 70, 58, 49, 45, 53, 57, 51, 63, 60, 49,
    52, 71, 74, 61, 50, 47, 55, 59, 54, 66, 62, 51, 55, 73, 77, 63, 54, 49, 58,
    60, 57, 70, 64)
  r <- arimfc(y, 6)
  # rmorie 1.5.0: R$morie_tps_arima_forecast(), the same input,
  # its values printed with 17 significant digits
  ref <- c(58.579828274521411, 60.799864215173187, 59.890564689732336,
        60.26300263626586, 60.110456615799329, 60.172937597065506)
  expect_equal(unname(c(r$fc)), c(ref), tolerance = 1e-10)
  ref <- c(8.4923465829595219, 14.373636721695604, 17.697741289767642,
        20.756784480260087, 23.321592697261575, 25.667387510573167)
  expect_equal(unname(c(r$se)), c(ref), tolerance = 1e-10)
  ref <- 283.87014784803432
  expect_equal(unname(c(r$aic)), c(ref), tolerance = 1e-08)
  ref <- 288.93678621037611
  expect_equal(unname(c(r$bic)), c(ref), tolerance = 1e-08)
})


test_that("sarfc() equals the reference hold-out forecast", {
  y <- c(52, 48, 61, 55, 47, 50, 66, 70, 58, 49, 45, 53, 57, 51, 63, 60, 49,
    52, 71, 74, 61, 50, 47, 55, 59, 54, 66, 62, 51, 55, 73, 77, 63, 54, 49, 58,
    60, 57, 70, 64)
  r <- sarfc(y, 4)
  # rmorie 1.5.0: R$morie_tps_sarima_forecast(), the same input,
  # its values printed with 17 significant digits
  ref <- c(61.633732308124713, 56.553996127378085, 68.770543598695468,
        64.63414379219671)
  expect_equal(unname(c(r$fc)), c(ref), tolerance = 1e-8)
  ref <- 78.309663098001835
  expect_equal(unname(c(r$aic)), c(ref), tolerance = 1e-08)
  ref <- 1.5631415193793825
  # the search stops within about 2e-7 of the optimum, and where it stops
  # differs between systems
  expect_equal(unname(c(r$mape)), c(ref), tolerance = 1e-6)
  ref <- 1.093355168837483
  expect_equal(unname(c(r$rmse)), c(ref), tolerance = 1e-6)
})


test_that("ousim() equals the reference simulation", {
  y <- c(52, 48, 61, 55, 47, 50, 66, 70, 58, 49, 45, 53, 57, 51, 63, 60, 49,
    52, 71, 74, 61, 50, 47, 55, 59, 54, 66, 62, 51, 55, 73, 77, 63, 54, 49, 58,
    60, 57, 70, 64)
  r <- ousim(y, 3, 10)
  # rmorie 1.5.0: R$morie_tps_langevin_simulate(), the same input,
  # its values printed with 17 significant digits
  ref <- 0.60476815820474661
  expect_equal(unname(c(r$theta)), c(ref), tolerance = 1e-10)
  ref <- 58.252367038686891
  expect_equal(unname(c(r$mu)), c(ref), tolerance = 1e-10)
  ref <- 7.6830239230095758
  expect_equal(unname(c(r$sigma)), c(ref), tolerance = 1e-10)
  ref <- c(64, 64, 64, 71.057121146829701, 56.185425039833866,
        63.313938871379015, 68.175512122477343, 60.541448943603505,
        59.437504200552091, 73.787369613033491, 58.429817439330357,
        74.228368662438015, 63.910460795110154, 68.347845857168778,
        82.13494284052824, 49.817975862714327, 60.100480960996258,
        66.667210442703308, 59.804849202005464, 56.798878512623624,
        41.168570546455072, 40.115791627890154, 67.820364507831911,
        49.14439503218388, 37.398379622692133, 60.71309914362768,
        63.984981259954637, 64.571023885118294, 55.917622084439763,
        58.541471893188607)
  expect_equal(unname(c(r$paths)), c(ref), tolerance = 1e-10)
  ref <- c(56.180007065314648, 58.541471893188607, 63.968068685925324)
  expect_equal(unname(c(r$q)), c(ref), tolerance = 1e-10)
})


test_that("fpdens() equals the reference density", {
  y <- c(52, 48, 61, 55, 47, 50, 66, 70, 58, 49, 45, 53, 57, 51, 63, 60, 49,
    52, 71, 74, 61, 50, 47, 55, 59, 54, 66, 62, 51, 55, 73, 77, 63, 54, 49, 58,
    60, 57, 70, 64)
  r <- fpdens(y, 16, 20)
  # rmorie 1.5.0: R$morie_tps_fokker_planck_grid(), the same input,
  # its values printed with 17 significant digits
  ref <- c(0, 7.7666666666666666, 15.533333333333333, 23.300000000000001,
        31.066666666666666, 38.833333333333336, 46.600000000000001,
        54.366666666666667, 62.133333333333333, 69.900000000000006,
        77.666666666666671, 85.433333333333337, 93.200000000000003,
        100.96666666666667, 108.73333333333333, 116.5)
  expect_equal(unname(c(r$grid)), c(ref), tolerance = 1e-10)
  ref <- c(0, 0, 0, 0, 0, 0, 0.0050731645512773422, 0.034271422215353689,
        0.05701804208933807, 0.031119481711771927, 0.0012732542391259186,
        0, 0, 0, 0, 0)
  expect_equal(unname(c(r$density)), c(ref), tolerance = 1e-10)
  ref <- 62.133333333333333
  expect_equal(unname(c(r$peak)), c(ref), tolerance = 1e-10)
  ref <- 48.802880740914446
  expect_equal(unname(c(r$statvar)), c(ref), tolerance = 1e-10)
})


test_that("catcor() equals the reference correlations (to its 3 decimals)", {
  set.seed(16)
  u <- list(b = sample(letters[1:6], 40, TRUE),
    a = sample(letters[1:6], 30, TRUE), c = sample(letters[2:7], 25, TRUE))
  r <- catcor(u)
  # rmorie 1.5.0: R$morie_tps_category_correlation_matrix(), the same input,
  # its values printed with 17 significant digits
  ref <- c(1, 0.46100000000000002, -0.040000000000000001, 0.46100000000000002,
        1, -0.085000000000000006, -0.040000000000000001,
        -0.085000000000000006, 1)
  # rmorie rounds the correlations to 3 decimals; crimlab does not
  expect_equal(round(unname(c(r$r)), 3), ref, tolerance = 1e-12)
})


test_that("lvfit() equals the reference rates", {
  y <- c(410, 455, 430, 470, 520, 498, 515)
  r <- lvfit(y, 2017:2023)
  # rmorie 1.5.0: { utils::assignInNamespace(), the same input,
  # its values printed with 17 significant digits
  ref <- 0.11125715482779004
  expect_equal(unname(c(r$alpha)), c(ref), tolerance = 1e-10)
  ref <- 0.00024452121940173635
  expect_equal(unname(c(r$beta)), c(ref), tolerance = 1e-10)
  ref <- 0.0039665880110530544
  expect_equal(unname(c(r$gamma)), c(ref), tolerance = 1e-10)
  ref <- 8.4395489596873491e-06
  expect_equal(unname(c(r$delta)), c(ref), tolerance = 1e-10)
  ref <- 299.09375565989762
  expect_equal(unname(c(r$period)), c(ref), tolerance = 1e-10)
  ref <- c(410, 455, 431.66666666666663, 451.66666666666663, 473.33333333333331,
        496, 510.99999999999994)
  expect_equal(unname(c(r$pred)), c(ref), tolerance = 1e-10)
})

# Behaviour of the count-dynamics functions.

test_that("arimfc() and sarfc() forecast and refuse short series", {
  set.seed(1)
  y <- rpois(60, 50 + 10 * sin(2 * pi * (1:60) / 12))
  f <- arimfc(y, h = 6)
  expect_length(f$fc, 6)
  expect_true(all(diff(f$se) >= 0))
  expect_equal(f$bic, f$aic + 3 * (log(60) - 2),
               tolerance = testthat_tolerance())
  s <- sarfc(y, h = 6)
  expect_length(s$actual, 6)
  expect_equal(s$rmse, sqrt(mean((s$actual - s$fc)^2)),
               tolerance = testthat_tolerance())
  expect_error(arimfc(1:5), "24 values")
  expect_error(sarfc(y, order = c(1, 0, 0), seasonal = c(0, 0, 0, 12)),
               "difference")
  expect_error(sarfc(y[1:30]), "36 values")
  # 60 - 40 = 20 training values, fewer than the 24 needed
  expect_error(sarfc(y, h = 40), "training values")
})

test_that("catcor() counts by unit and correlates the categories", {
  u <- list(b = c("x", "y", "y"), a = c("x", "x", "y", "z"))
  f <- catcor(u)
  expect_identical(f$cats, c("a", "b"))
  expect_equal(unname(f$counts), cbind(c(2, 1, 1), c(1, 2, 0)),
               tolerance = testthat_tolerance())
  expect_equal(f$r[1, 2], cor(c(2, 1, 1), c(1, 2, 0)),
               tolerance = testthat_tolerance())
  expect_equal(unname(diag(f$r)), c(1, 1), tolerance = testthat_tolerance())
  expect_error(catcor(list(a = "x")), "two categories")
  expect_error(catcor(1:3), "named list")
})

test_that("ousim() and fpdens() share the Ornstein-Uhlenbeck fit", {
  set.seed(2)
  y <- rpois(100, 20)
  s <- ousim(y, npath = 4, days = 10)
  d <- fpdens(y, ngrid = 32, nstep = 50)
  expect_equal(s$theta, d$theta, tolerance = testthat_tolerance())
  expect_equal(s$halflife, log(2) / s$theta, tolerance = testthat_tolerance())
  expect_identical(dim(s$paths), c(4L, 10L))
  expect_equal(s$paths[, 1], rep(y[100], 4), tolerance = testthat_tolerance())
  expect_equal(sum(d$density) * diff(d$grid[1:2]), 1,
               tolerance = testthat_tolerance())
  expect_equal(d$statvar, d$sigma^2 / (2 * d$theta),
               tolerance = testthat_tolerance())
  expect_error(ousim(1:2), "three values")
})

test_that("lvfit() follows its formulas", {
  y <- c(410, 455, 430, 470, 520, 498, 515)
  f <- lvfit(y)
  expect_equal(f$alpha, median(diff(log(y))) + 0.05,
               tolerance = testthat_tolerance())
  expect_equal(f$period, 2 * pi / sqrt(f$alpha * f$gamma),
               tolerance = testthat_tolerance())
  expect_equal(f$year, 1:7, tolerance = testthat_tolerance())
  expect_error(lvfit(1:3), "five years")
  expect_error(lvfit(y, 1:2), "one value per count")
  expect_output(print(f), "lvfit:")
  expect_s3_class(summary(f), "summary.svfit")
})

# The fits below are compared with stats::arima() on the training values
# (the default method, CSS for the start, then exact maximum likelihood,
# with the same BFGS steps) and predict() for the forecasts.
ts_y <- c(52, 48, 61, 55, 47, 50, 66, 70, 58, 49, 45, 53, 57, 51, 63, 60, 49,
  52, 71, 74, 61, 50, 47, 55, 59, 54, 66, 62, 51, 55, 73, 77, 63, 54, 49, 58,
  60, 57, 70, 64)
ts_cmp <- function(y, order, seasonal, tol, h = 4L) {
  r <- sarfc(y, h, order = order, seasonal = seasonal)
  f <- stats::arima(y[seq_len(length(y) - h)], order = order,
                    seasonal = list(order = seasonal[1:3],
                                    period = seasonal[4]))
  expect_equal(unname(r$coef), unname(stats::coef(f)), tolerance = tol)
  expect_equal(unname(r$fc), as.numeric(stats::predict(f, h)$pred),
               tolerance = tol)
  expect_equal(r$aic, f$aic, tolerance = tol)
  r
}

test_that("sarfc() with two AR terms and a seasonal AR term is arima()", {
  # Both fits end where a quasi-Newton search stops, and where it stops
  # depends on the last digits of the system's exp() and log(): on Linux
  # the forecasts agree to 3e-12, on macOS and Windows to 4e-8. The bound
  # is that of a search with a numerical gradient, as in cli/tol.tsv.
  ts_cmp(ts_y, c(2, 1, 0), c(1, 0, 0, 4), 1e-6)
  ts_cmp(ts_y, c(2, 1, 0), c(1, 0, 0, 3), 1e-6)
  ts_cmp(ts_y, c(3, 1, 1), c(0, 0, 0, 4), 1e-6)
})

test_that("sarfc() fits a random walk and an AR(1) without CSS steps", {
  # no ARMA term: the forecast is the last training value
  r <- ts_cmp(ts_y, c(0, 1, 0), c(0, 0, 0, 12), 1e-12)
  expect_identical(r$fc, rep(58, 4))
  expect_length(r$coef, 0L)
  # one AR term, a state of one: the variance 1 / (1 - phi^2)
  ts_cmp(ts_y, c(1, 1, 0), c(0, 0, 0, 12), 1e-12)
})

test_that("an MA start with a root inside the unit circle is inverted", {
  set.seed(29)
  w <- round(rnorm(40, 50, 5))
  # the CSS start of the twice-differenced series has a root of modulus
  # 0.989 (stats::arima(method = "CSS")): it is moved outside first
  css <- stats::coef(stats::arima(w[1:36], c(0, 2, 2), method = "CSS"))
  expect_lt(min(Mod(polyroot(c(1, css)))), 1)
  # 6e-9 from arima() in the coefficients and the forecasts, near the
  # boundary of invertibility (an MA root of modulus 1.000002)
  ts_cmp(w, c(0, 2, 2), c(0, 0, 0, 12), 1e-7)
})

test_that("sarfc() refuses as arima() does", {
  ef <- function(call, msg) expect_error(call, msg, fixed = TRUE)
  # a constant series: the CSS objective at the start is log(0)
  ef(sarfc(rep(5, 40), order = c(1, 1, 0), seasonal = c(0, 0, 0, 12)),
     "sarfc: initial value in 'vmmin' is not finite")
  # differences that grow by 10% a step: an AR coefficient of 1.1
  g <- cumsum(1.1^(1:40))
  ef(sarfc(g, order = c(1, 1, 0), seasonal = c(0, 0, 0, 12)),
     "sarfc: non-stationary AR part from CSS")
  ef(sarfc(ts_y, 4, order = c(2, 1, 0), seasonal = c(1, 0, 0, 12)),
     "sarfc: non-stationary seasonal AR part from CSS")
  # differences of 1e154: the sum of squares is finite at the start and
  # overflows a step of 1e-3 away
  d <- sqrt(0.999 * .Machine$double.xmax / 34)
  ef(sarfc(cumsum(rep(d, 40)), 4, order = c(1, 1, 0),
           seasonal = c(0, 0, 0, 12)),
     "sarfc: non-finite finite-difference value [1]")
  ef(sarfc(1:40, order = c(1, 1)),
     "sarfc: `order` must have 3 whole numbers from 0.")
  ef(sarfc(1:40, order = c(1, 1.5, 0)),
     "sarfc: `order` must have 3 whole numbers from 0.")
  ef(sarfc(1:40, seasonal = c(0, 1, 1, 0)),
     "sarfc: the period, the fourth value of `seasonal`, must be at least 1.")
  ef(catcor(list(a = "x", b = "x")),
     "catcor: `units` must name two units at least.")
  ef(ousim(1:10, npath = 1e4, days = 1e5),
     "ousim: `npath` times the steps must be at most 1e8.")
  ef(lvfit(c(1, 2, -1, 3, 4)), "lvfit: `y` must hold counts from 0.")
})

test_that("BUG: an AR(1) fit at the unit root is refused", {
  # constant differences of 1e154: the likelihood drives the AR
  # coefficient to 1, where arima() stops ("system is exactly singular")
  # and a state of more than one gives "not stationary"
  d <- sqrt(0.999 * .Machine$double.xmax / 38)
  expect_error(
    sarfc(cumsum(rep(d, 40)), 4, order = c(1, 1, 0),
          seasonal = c(0, 0, 0, 12)),
    "sarfc: the fitted ARMA part is not stationary.", fixed = TRUE
  )
})

test_that("a horizon or a grid past its bound is refused before it is made", {
  y <- c(3, 5, 4, 6, 8, 7, 9, 11, 10, 12, 14, 13, 15, 17, 16, 18, 20, 19,
         21, 23, 22, 24, 26, 25)
  expect_error(arimfc(y, h = 1e9),
               "arimfc: `h` must be at most 67108864.", fixed = TRUE)
  expect_error(fpdens(y, ngrid = 1e9),
               "fpdens: `ngrid` must be at most 67108864.", fixed = TRUE)
})
