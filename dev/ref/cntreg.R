# Reference values of tests/testthat/test-cntreg.R: MASS::glm.nb(),
# sandwich::vcovCL() and geepack::geeglm() on the counts the tests build.
# From the package root, with MASS, sandwich and geepack in the library:
#   Rscript dev/ref/cntreg.R
# the packages it needs (dev/ref/all.R skips it without them)
need <- c("MASS", "sandwich", "geepack")
source(file.path("dev", "ref", "put.R"))

set.seed(3)
cn <- local({
  id <- rep(1:8, each = 10)
  x <- rbinom(80, 1, 0.5)
  f <- factor(sample(c("a", "b", "c"), 80, TRUE))
  u <- rnorm(8, 0, 0.4)[id]
  y <- rnbinom(80, mu = exp(1 + 0.3 * x + u), size = 3)
  e <- runif(80, 1, 3)
  data.frame(y, x, f, id, e)
})
g <- glm(y ~ factor(x) + f + offset(log(e)), family = poisson, data = cn)
nb <- MASS::glm.nb(y ~ factor(x) + f + offset(log(e)), data = cn)
co <- unname(coef(summary(nb)))
ge <- geepack::geeglm(
  y ~ factor(x) + f + offset(log(e)),
  family = poisson,
  id = id,
  data = cn,
  corstr = "exchangeable",
  control = geepack::geese.control(epsilon = 1e-12, maxit = 100)
)
put(list(
  inp = sum(cn$y, cn$x, as.integer(cn$f), cn$e),
  nb = list(
    est = co[, 1],
    se = co[, 2],
    theta = nb$theta,
    setheta = nb$SE.theta,
    aic = AIC(nb),
    loglik = as.numeric(logLik(nb))
  ),
  cl = list(
    pois = unname(sqrt(diag(sandwich::vcovCL(g, cluster = cn$id)))),
    nb = unname(sqrt(diag(sandwich::vcovCL(nb, cluster = cn$id))))
  ),
  gee = list(
    est = unname(coef(ge)),
    se = unname(sqrt(diag(ge$geese$vbeta))),
    alpha = unname(ge$geese$alpha)
  )
), "cntreg", need)
