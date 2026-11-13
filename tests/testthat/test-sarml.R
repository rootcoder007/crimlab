# Global and bivariate Moran's I, spatial regression and its effects. The
# kernels are checked in src/t_sreg.cpp; here the binding against spdep
# and spatialreg (stored, versions named) and the formulas written out.

sdat <- function() {
  set.seed(91)
  n <- 60
  p <- cbind(x = runif(n), y = runif(n))
  w <- spwt(p, k = 4)
  v <- p[, 1] + p[, 2] + rnorm(n, sd = 0.3)
  u <- p[, 1] + rnorm(n, sd = 0.3)
  W <- matrix(0, n, n)
  W[cbind(w$i, w$j)] <- w$w
  x1 <- rnorm(n)
  x2 <- runif(n)
  y <- as.vector(solve(diag(n) - 0.45 * W, 1 + 2 * x1 - x2 + rnorm(n)))
  list(p = p, w = w, v = v, u = u, W = W, x = cbind(x1, x2), y = y)
}

test_that("global Moran's I agrees with spdep::moran.test()", {
  d <- sdat()
  g <- gmoran(d$v, d$w)
  # spdep 1.4.2: moran.test(v, nb2listw(knn2nb(knearneigh(p, k = 4)))),
  # and with randomisation = FALSE
  expect_equal(g$i, 0.565255385168907, tolerance = 1e-13)
  expect_equal(g$e, -0.0169491525423729, tolerance = 1e-13)
  expect_equal(g$vr, 0.00670334086876704, tolerance = 1e-12)
  expect_equal(g$zr, 7.11098780707048, tolerance = 1e-12)
  expect_equal(g$pr, 5.76076597911219e-13, tolerance = 1e-10)
  expect_equal(g$vn, 0.00662436835090727, tolerance = 1e-12)
  expect_equal(g$zn, 7.15324910546279, tolerance = 1e-12)
  expect_equal(g$pn, 4.23738043256574e-13, tolerance = 1e-10)
  expect_identical(g$s0, 60)
  expect_true(is.nan(g$pperm))
  expect_identical(gmoran(d$v, d$W)[, -1], g[, -1])
  # a sweep: one row per column, each as alone
  m <- gmoran(data.frame(a = d$v, b = d$u, k = 1), d$w, nsim = 99, seed = 4)
  expect_identical(m$var, c("a", "b", "k"))
  expect_equal(m[1, 2:9], g[, 2:9], ignore_attr = TRUE)
  expect_true(all(is.na(m[3, 2:10])))
  expect_identical(m$pperm[1], 0.01)
  set.seed(1)
  a <- gmoran(d$u, d$w, nsim = 49)
  set.seed(1)
  expect_identical(gmoran(d$u, d$w, nsim = 49), a)
  expect_identical(gmoran(cbind(d$u, d$v), d$w)$var, c("x1", "x2"))
})

test_that("bivariate Moran agrees with spdep::moran_bv() and with Moran", {
  d <- sdat()
  b <- bvmor(d$u, d$v, d$w, nsim = 199, seed = 2)
  # spdep 1.4.2: moran_bv(u, v, lw, nsim = 9)$t0
  expect_equal(b[["i"]], 0.251663587577159, tolerance = 1e-13)
  expect_named(b, c("i", "mean", "sd", "z", "p", "ptwo"))
  expect_equal(bvmor(d$v, d$v, d$w, nsim = 0)[["i"]], gmoran(d$v, d$w)$i,
               tolerance = 1e-14)
  expect_true(all(is.nan(bvmor(d$u, d$v, d$w, nsim = 0)[-1])))
  expect_lt(abs(b[["mean"]]), 4 * b[["sd"]] / sqrt(199))
  expect_lt(b[["p"]], 0.05)
  expect_equal(b[["z"]], (b[["i"]] - b[["mean"]]) / b[["sd"]],
               tolerance = 1e-15)
})

test_that("maximum likelihood agrees with spatialreg", {
  d <- sdat()
  # spatialreg 1.4.3: lagsarlm(y ~ x1 + x2, method = "eigen"): coef,
  # s2, LL; then the standard errors of beta and rho
  f <- sarml(d$y, d$x, d$w)
  expect_s3_class(f, "sarml")
  cf <- stats::setNames(f$coef$est, f$coef$term)
  se <- stats::setNames(f$coef$se, f$coef$term)
  expect_equal(unname(cf[c("rho", "(Intercept)", "x1", "x2")]),
               c(0.606896909014245, 0.924452296102293, 1.92210970645926,
                 -1.02603714241967), tolerance = 1e-7)
  expect_equal(f$fit[["s2"]], 0.838767263794935, tolerance = 1e-8)
  expect_equal(f$fit[["loglik"]], -82.7377286381847, tolerance = 1e-10)
  expect_equal(unname(se[c("(Intercept)", "x1", "x2", "rho")]),
               c(0.249009387318348, 0.106779808361906, 0.425568117868853,
                 0.0659103883769072), tolerance = 1e-6)
  # the error model, errorsarlm()
  e <- sarml(d$y, d$x, d$w, model = "error")
  ce <- stats::setNames(e$coef$est, e$coef$term)
  expect_equal(unname(ce[c("lambda", "(Intercept)", "x1", "x2")]),
               c(0.69630097288685, 1.39926889033533, 1.73624510356828,
                 -1.19886497187181), tolerance = 1e-7)
  expect_equal(e$fit[["loglik"]], -89.3493461408994, tolerance = 1e-10)
  expect_equal(e$coef$se, c(0.487425170797838, 0.107743919578138,
                            0.462513590327871, 0.0860956829751044),
               tolerance = 1e-6)
  # the Durbin lag model, lagsarlm(Durbin = TRUE)
  s <- sarml(d$y, d$x, d$w, durbin = TRUE)
  expect_identical(s$coef$term, c("(Intercept)", "x1", "x2", "lag.x1",
                                  "lag.x2", "rho"))
  expect_equal(s$coef$est, c(0.440746223227439, 1.92912193088266,
                             -1.14638660059137, 0.0983546307148609,
                             1.12140967062863, 0.593741589310601),
               tolerance = 1e-7)
  expect_equal(s$fit[["loglik"]], -81.5677250947261, tolerance = 1e-10)
  expect_equal(s$coef$se, c(0.413950462439741, 0.111749703104606,
                            0.425905405868202, 0.411825839293856,
                            0.733730003847498, 0.101695915617604),
               tolerance = 1e-6)
  # sacsarlm(): the optimum (two optimisers: 1e-5), and the standard
  # errors from the expected information
  a <- sarml(d$y, d$x, d$w, model = "sac")
  expect_equal(a$coef$est, c(0.999938119237316, 1.92183719713386,
                             -1.12640663489178, 0.57544521359594,
                             0.13969130258284), tolerance = 1e-5)
  expect_equal(a$fit[["loglik"]], -82.5469720464817, tolerance = 1e-10)
  expect_equal(a$coef$se, c(0.264418141937921, 0.105223420163633,
                            0.431219359969637, 0.0945683525034107,
                            0.233850154788925), tolerance = 1e-5)
  # the OLS columns are lm()'s
  o <- summary(lm(d$y ~ d$x))$coefficients
  expect_equal(f$coef$ols[1:3], unname(o[, 1]), tolerance = 1e-12)
  expect_equal(f$coef$olsse[1:3], unname(o[, 2]), tolerance = 1e-12)
  expect_equal(f$fit[["llols"]], as.numeric(logLik(lm(d$y ~ d$x))),
               tolerance = 1e-12)
  expect_output(print(f), "sarml: lag model by ml")
  expect_output(print(s), "lag Durbin model")
})

test_that("two-stage least squares agrees with spatialreg::stsls()", {
  d <- sdat()
  t2 <- sarml(d$y, data.frame(d$x), d$w, method = "2sls")
  # spatialreg 1.4.3: stsls(y ~ x1 + x2): coefficients, standard errors
  expect_equal(t2$coef$est, c(0.930389649997745, 1.92042562474309,
                              -1.03775656923328, 0.62734927957874),
               tolerance = 1e-12)
  expect_equal(t2$coef$se, c(0.253255579989461, 0.109837579907715,
                             0.442238227293904, 0.111268901328071),
               tolerance = 1e-12)
  expect_equal(t2$fit[["s2"]], 0.832592540249328, tolerance = 1e-12)
  expect_true(is.nan(t2$fit[["loglik"]]))
  expect_output(print(t2), "by 2sls")
  t3 <- sarml(d$y, d$x, d$w, method = "2sls", durbin = TRUE)
  expect_identical(t3$coef$term[4:5], c("lag.x1", "lag.x2"))
  expect_identical(sarml(d$y, d$x[, 1], d$w)$coef$term,
                   c("(Intercept)", "x", "rho"))
})

test_that("effects agree with spatialreg::impacts() and the delta method", {
  d <- sdat()
  s <- sarml(d$y, d$x, d$w, durbin = TRUE)
  e <- spill(s)
  # spatialreg 1.4.3: impacts(lagsarlm(Durbin = TRUE), listw = lw); exact
  # at spatialreg's coefficients, to the optimiser's agreement at ours
  ref <- c(2.16925240324877, -1.06129099209691, 2.8213558112544,
           0.99981059118569, 4.99060821450317, -0.0614804009112224)
  expect_identical(e$term, c("x1", "x2"))
  expect_equal(c(e$direct, e$indirect, e$total), ref, tolerance = 1e-7)
  r <- spill(w = d$w, rho = 0.593741589310601,
             b = c(x1 = 1.92912193088266, x2 = -1.14638660059137),
             th = c(0.0983546307148609, 1.12140967062863))
  expect_equal(c(r$direct, r$indirect, r$total), ref, tolerance = 1e-13)
  l <- spill(w = d$w, rho = 0.606896909014245,
             b = c(x1 = 1.92210970645926, x2 = -1.02603714241967))
  expect_equal(c(l$direct, l$indirect, l$total),
               c(2.15753580218812, -1.15170942725383, 2.73204596538163,
                 -1.45838743015513, 4.88958176756975, -2.61009685740896),
               tolerance = 1e-13)
  # row-standardised: the total is (beta + theta) / (1 - rho)
  cf <- stats::setNames(s$coef$est, s$coef$term)
  expect_equal(e$total[1], (cf[["x1"]] + cf[["lag.x1"]]) / (1 - cf[["rho"]]),
               tolerance = 1e-12)
  # delta method against a numerical gradient of the effects
  eff <- function(th) {
    S <- solve(diag(60) - th[1] * d$W)
    M <- S %*% (th[2] * diag(60) + th[3] * d$W)
    c(mean(diag(M)), mean(rowSums(M)) - mean(diag(M)), mean(rowSums(M)))
  }
  id <- c("rho", "x1", "lag.x1")
  th <- cf[id]
  h <- 1e-6
  J <- vapply(1:3, function(k) {
    e1 <- th
    e2 <- th
    e1[k] <- e1[k] + h
    e2[k] <- e2[k] - h
    (eff(e1) - eff(e2)) / (2 * h)
  }, numeric(3))
  V <- s$vcov[id, id]
  expect_equal(c(e$sed[1], e$sei[1], e$set[1]), sqrt(diag(J %*% V %*% t(J))),
               tolerance = 1e-7)
  # given values, without a fit: no standard errors
  g <- spill(w = d$w, rho = 0.4, b = c(gent = 0.1), th = 0.3)
  expect_identical(g$term, "gent")
  expect_equal(g$total, 0.4 / 0.6, tolerance = 1e-12)
  expect_true(is.nan(g$sed))
  expect_identical(spill(w = d$W, rho = 0.4, b = 2)$term, "x1")
  lw <- structure(list(neighbours = lapply(split(d$w$j, d$w$i), as.integer),
                       weights = split(d$w$w, d$w$i)), class = "listw")
  expect_identical(spill(w = lw, rho = 0.4, b = 2),
                   spill(w = d$W, rho = 0.4, b = 2))
  expect_equal(spill(w = d$W, rho = 0.4, b = 2)$total, 2 / 0.6,
               tolerance = 1e-12)
  sc <- spill(sarml(d$y, d$x, d$w, model = "sac"))
  expect_identical(nrow(sc), 2L)
})

test_that("the spatial functions refuse bad input in words", {
  d <- sdat()
  nbs <- lapply(1:4, function(i) as.integer(i %% 4 + 1))
  lw <- structure(list(neighbours = nbs,
                       weights = as.list(rep(1, 4))), class = "listw")
  expect_identical(gmoran(c(1, 3, 2, 5), lw)$s0, 4)
  errs("gmoran", list(
    list("at least 4 units", 1:3, d$w),
    list("no positive weight", c(1, 2, 3, 5), matrix(0, 4, 4)),
    list("must be a 4 x 4 matrix", 1:4, matrix(0, 3, 3)),
    list("its own neighbour", 1:4, diag(4)),
    list("must be the list spwt", 1:4, list(1)),
    list("is for 60 units", 1:4, d$w),
    list("lists a pair of units twice", 1:4,
         list(n = 4, i = c(1, 1), j = c(2, 2), w = c(1, 1))),
    list("unit numbers from 1 to 4", 1:4,
         list(n = 4, i = 1.5, j = 2, w = 1)),
    list("one row of neighbours", 1:4,
         structure(list(neighbours = list(), weights = list()),
                   class = "listw")),
    list("integer neighbour lists", 1:4,
         structure(list(neighbours = as.list(c(2, 1, 4, 3)),
                        weights = as.list(rep(1, 4))), class = "listw")),
    list("one weight per neighbour", 1:4,
         structure(list(neighbours = as.list(c(2L, 1L, 4L, 3L)),
                        weights = as.list(rep("a", 4))), class = "listw")),
    list("outside 1 to 4", 1:4,
         structure(list(neighbours = as.list(c(2L, 1L, 9L, 3L)),
                        weights = as.list(rep(1, 4))), class = "listw")),
    list("finite weights from 0", 1:4,
         structure(list(neighbours = as.list(c(2L, 1L, 4L, 3L)),
                        weights = as.list(c(-1, 1, 1, 1))), class = "listw"))
  ))
  errs("bvmor", list(
    list("at least 3 units", 1:2, 1:2, d$w),
    list("constant", rep(1, 60), d$v, d$w)
  ))
  errs("sarml", list(
    list("at least 4 units", 1:3, 1:3, d$w),
    list("one row per value", d$y, 1:3, d$w),
    list("lag model only", d$y, d$x, d$w, model = "error",
         method = "2sls"),
    list("a column of `x` that varies", d$y, rep(1, 60), d$w,
         method = "2sls"),
    list("too many columns", d$y[1:4], matrix(rnorm(16), 4), d$W[1:4, 1:4]),
    list("cannot be fitted", d$y, cbind(d$x, d$x[, 1]), d$w),
    list("at most 5000 units", rep(1, 5001), rep(1, 5001), d$w)
  ))
  errs("spill", list(
    list("give a sarml", NULL),
    list("in an error model", sarml(d$y, d$x, d$w, model = "error")),
    list("`b` must be finite numbers", NULL, d$w, 0.4, numeric(0)),
    list("at least 2 units", NULL, matrix(0, 1, 1), 0.4, 1),
    list("must be the list spwt", NULL, list(1), 0.4, 1),
    list("singular", NULL, d$W, 1, 1),
    list("at most 5000 units", NULL, list(n = 5001, i = 1, j = 2, w = 1),
         0.4, 1)
  ))
})
