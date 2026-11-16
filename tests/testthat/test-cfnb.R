# The negative binomial control function of ppmlcf(family = "nb")
# (R/cfnb.R).

test_that("ppmlcf(family = \"nb\") is the negative binomial control function", {
  # cities over three years: lead use instrumented by distance to a
  # refinery, homicides negative binomial with the population as exposure
  # (Feigenbaum and Muller 2016). Written out here: the first stage by
  # lm(), the second-stage likelihood maximised by optim() over the
  # coefficients and log theta; 1e-5 is the optimiser's agreement
  set.seed(62)
  n <- 240
  d <- data.frame(city = rep(1:80, 3), yr = rep(1:3, each = 80))
  dist <- runif(80, 0, 5)[d$city]
  u <- rnorm(80)[d$city]
  d$dist <- dist
  d$lead <- 1.5 - 0.3 * dist + 0.5 * u + rnorm(n, 0, 0.3)
  d$pop <- exp(rnorm(80, 10, 0.5))[d$city]
  mu <- d$pop * exp(-9 + 0.4 * d$lead + 0.1 * d$yr + 0.3 * u)
  d$hom <- stats::rnbinom(n, size = 3, mu = mu)
  r <- ppmlcf(d, "hom", "lead", "dist", fe = "yr", expo = "pop",
              clus = "city", family = "nb", se = "naive")
  vh <- stats::residuals(stats::lm(lead ~ dist + factor(yr), d))
  xm <- cbind(d$lead, vh, 1, d$yr == 2, d$yr == 3)
  nll <- function(q) {
    b <- q[1:5]
    th <- exp(q[6])
    m <- exp(drop(xm %*% b) + log(d$pop))
    -sum(lgamma(d$hom + th) - lgamma(th) - lgamma(d$hom + 1) +
           th * log(th / (th + m)) + d$hom * log(m / (th + m)))
  }
  o <- stats::optim(c(0, 0, -9, 0, 0, 0), nll, method = "BFGS",
                    control = list(maxit = 1000, reltol = 1e-14))
  o <- stats::optim(o$par, nll, method = "BFGS",
                    control = list(maxit = 1000, reltol = 1e-15))
  expect_identical(r$coef$term, c("lead", "vhat:lead"))
  expect_equal(r$coef$est, o$par[1:2], tolerance = 1e-5)
  expect_equal(r$theta, exp(o$par[6]), tolerance = 1e-4)
  expect_identical(r$family, "nb")
  b <- ppmlcf(d, "hom", "lead", "dist", fe = "yr", expo = "pop",
              clus = "city", family = "nb", nboot = 5, seed = 1)
  expect_identical(dim(b$draws), c(5L, 2L))
  expect_error(ppmlcf(d, "hom", "lead", "dist", fe = "yr", family = "nb",
                      how = "plug"),
               "ppmlcf: `family = \"nb\"` needs `how = \"cf\"`.", fixed = TRUE)
})
