# The simulation of the alert complexity estimators (cxsim()).

test_that("cxsim() reports bias, spread and coverage", {
  s <- cxsim(reps = 2, n = 300, ncl = 6, seed = 1)
  expect_identical(s$res$est, c("crude", "mixed", "matched"))
  expect_true(all(s$res$ok <= 2))
  expect_true(is.na(s$res$cover[1]))
  s2 <- cxsim(reps = 2, n = 300, ncl = 6, conf = 1, seed = 1)
  expect_identical(s2$conf, 1)
  errs("cxsim", list(
    list("`reps` must be between 1 and 1e\\+07", reps = 0),
    list("`reps` must be between 1 and 1e\\+07", reps = 2^53),
    list("`n` must be between 50 and 1e\\+07", n = 1e308),
    list("`beta` must be between -20 and 20", beta = 1e308),
    list("`sdcl` must be between 0 and 10", sdcl = 1e308)
  ))
  # a 1 x 1 matrix is a number
  f <- cxsim(reps = 2, n = 200, ncl = 5, beta = matrix(log(1.3)), seed = 1)
  expect_identical(f, cxsim(reps = 2, n = 200, ncl = 5, seed = 1))
})

test_that("cxsim() counts failed fits and prints its settings", {
  # every count model failing: those estimates are missing and not
  # counted, the crude ratio still given
  local_mocked_bindings(cntmix = function(...) stop("no fit"))
  f <- cxsim(reps = 2, n = 200, ncl = 5, seed = 1)
  expect_identical(f$res$ok, c(2L, 0L, 0L))
  expect_true(all(is.finite(f$draws$b[f$draws$est == "crude"])))
  expect_true(all(is.na(f$draws[f$draws$est != "crude", c("b", "lo", "hi")])))
  prt(f, "^cxsim: 2 tables of 200 person-years, log rate ratio 0.2624\n")
  prt(
    cxsim(reps = 2, n = 200, ncl = 5, conf = 0.5, seed = 1),
    paste0(
      "^cxsim: 2 tables of 200 person-years, log rate ratio 0.2624, ",
      "unobserved confounder 0.5\n"
    )
  )
})

test_that("cxsim() leaves the session's random stream as it was", {
  set.seed(12)
  s <- .Random.seed
  a <- cxsim(reps = 1, n = 100, ncl = 5, b0 = 0, seed = 2^40)
  expect_identical(.Random.seed, s)
  expect_identical(cxsim(reps = 1, n = 100, ncl = 5, b0 = 0, seed = 2^40), a)
})
