# The simulation of the alert complexity estimators (cxsim()).

test_that("cxsim() reports bias, spread and coverage", {
  s <- cxsim(reps = 2, n = 300, ncl = 6, seed = 1)
  expect_identical(s$res$est, c("crude", "mixed", "matched"))
  expect_true(all(s$res$ok <= 2))
  expect_true(is.na(s$res$cover[1]))
  s2 <- cxsim(reps = 2, n = 300, ncl = 6, conf = 1, seed = 1)
  expect_identical(s2$conf, 1)
  errs("cxsim", list(list("at least 1", reps = 0)))
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
