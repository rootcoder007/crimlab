test_that("gocoef() lays the coefficients out by threshold", {
  set.seed(64)
  n <- 300
  d <- data.frame(race = rbinom(n, 1, 0.4), age = rnorm(n))
  u <- runif(n)
  d$y <- 1 + (u < plogis(0.5 + 0.8 * d$race)) +
    (u < plogis(-1 + 0.1 * d$race - 0.3 * d$age))
  r <- gologt(d, "y", c("race", "age"), par = "age")
  g <- gocoef(r)
  cf <- r$coef
  expect_identical(dimnames(g$est),
                   list(c("1|2", "2|3"), c("(Intercept)", "race", "age")))
  expect_identical(g$est[, "race"], setNames(cf$est[3:4], c("1|2", "2|3")))
  expect_identical(unname(g$est[, "age"]), rep(cf$est[5L], 2))
  expect_identical(unname(g$se[, "(Intercept)"]), cf$se[1:2])
  expect_identical(g$wald, r$wald)
  h <- gocoef(r, "race")
  expect_identical(colnames(h$est), "race")
  expect_identical(h$wald$term, "race")
  expect_error(gocoef(r, "sex"),
               paste("gocoef: `v` must name terms of the fit: (Intercept),",
                     "race, age."),
               fixed = TRUE)
  expect_error(gocoef(list()), "gocoef: `object` must be a result of gologt().",
               fixed = TRUE)
})
