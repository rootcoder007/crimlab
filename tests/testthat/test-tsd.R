# Counts over time and across categories (src/tsd.cpp) against the values
# of the rmorie functions they port. The ARIMA fits stop where the BFGS
# optimiser of optim() stops (relative tolerance 1.5e-8 on the objective),
# so their forecasts are compared at 1e-4.

test_that("arimfc() equals the reference forecast", {
  y <- c(52, 48, 61, 55, 47, 50, 66, 70, 58, 49, 45, 53, 57, 51, 63, 60, 49,
    52, 71, 74, 61, 50, 47, 55, 59, 54, 66, 62, 51, 55, 73, 77, 63, 54, 49, 58,
    60, 57, 70, 64)
  r <- arimfc(y, 6)
  # rmorie 1.5.0: R$morie_tps_arima_forecast(), the same input,
  # its values printed with 17 significant digits
  ref <- c(58.579828274521411, 60.799864215173187, 59.890564689732336,
        60.26300263626586, 60.110456615799329, 60.172937597065506)
  expect_equal(unname(c(r$fc)), c(ref), tolerance = 0.0001)
  ref <- c(8.4923465829595219, 14.373636721695604, 17.697741289767642,
        20.756784480260087, 23.321592697261575, 25.667387510573167)
  expect_equal(unname(c(r$se)), c(ref), tolerance = 0.0001)
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
  expect_equal(unname(c(r$fc)), c(ref), tolerance = 0.0001)
  ref <- 78.309663098001835
  expect_equal(unname(c(r$aic)), c(ref), tolerance = 1e-08)
  ref <- 1.5631415193793825
  expect_equal(unname(c(r$mape)), c(ref), tolerance = 0.0001)
  ref <- 1.093355168837483
  expect_equal(unname(c(r$rmse)), c(ref), tolerance = 0.0001)
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
        58.541471893188607 )
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
        100.96666666666667, 108.73333333333333, 116.5 )
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
  expect_equal(unname(c(r$r)), c(ref), tolerance = 0.0005)
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
  expect_equal(f$bic, f$aic + 3 * (log(60) - 2))
  s <- sarfc(y, h = 6)
  expect_length(s$actual, 6)
  expect_equal(s$rmse, sqrt(mean((s$actual - s$fc)^2)))
  expect_error(arimfc(1:5), "24 values")
  expect_error(sarfc(y, order = c(1, 0, 0), seasonal = c(0, 0, 0, 12)),
               "difference")
  expect_error(sarfc(y[1:30]), "36 values")
  expect_error(sarfc(y, h = 30), "training values")
})

test_that("catcor() counts by unit and correlates the categories", {
  u <- list(b = c("x", "y", "y"), a = c("x", "x", "y", "z"))
  f <- catcor(u)
  expect_identical(f$cats, c("a", "b"))
  expect_equal(unname(f$counts), cbind(c(2, 1, 1), c(1, 2, 0)))
  expect_equal(f$r[1, 2], cor(c(2, 1, 1), c(1, 2, 0)))
  expect_equal(unname(diag(f$r)), c(1, 1))
  expect_error(catcor(list(a = "x")), "two categories")
  expect_error(catcor(1:3), "named list")
})

test_that("ousim() and fpdens() share the Ornstein-Uhlenbeck fit", {
  set.seed(2)
  y <- rpois(100, 20)
  s <- ousim(y, npath = 4, days = 10)
  d <- fpdens(y, ngrid = 32, nstep = 50)
  expect_equal(s$theta, d$theta)
  expect_equal(s$halflife, log(2) / s$theta)
  expect_identical(dim(s$paths), c(4L, 10L))
  expect_equal(s$paths[, 1], rep(y[100], 4))
  expect_equal(sum(d$density) * diff(d$grid[1:2]), 1)
  expect_equal(d$statvar, d$sigma^2 / (2 * d$theta))
  expect_error(ousim(1:2), "three values")
})

test_that("lvfit() follows its formulas", {
  y <- c(410, 455, 430, 470, 520, 498, 515)
  f <- lvfit(y)
  expect_equal(f$alpha, median(diff(log(y))) + 0.05)
  expect_equal(f$period, 2 * pi / sqrt(f$alpha * f$gamma))
  expect_equal(f$year, 1:7)
  expect_error(lvfit(1:3), "five years")
  expect_error(lvfit(y, 1:2), "one value per count")
  expect_output(print(f), "lvfit:")
  expect_s3_class(summary(f), "summary.svfit")
})
