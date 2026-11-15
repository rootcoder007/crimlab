# Spatial epidemiology against the values of rmorie and of surveillance,
# fields, MASS, stats, pscl, deSolve, survival and spdep, computed once by
# dev/ref/sep.R and stored in ref/sep.txt.

sin_ <- function() {
  set.seed(51)
  obs <- rpois(70, 6) + c(rep(0, 66), 12, 0, 9, 0)
  rec <- cbind(x = runif(6), y = runif(6))
  src <- cbind(x = runif(5), y = runif(5))
  ll <- cbind(
    lon = -79.4 + runif(4, -0.1, 0.1),
    lat = 43.7 + runif(4, -0.1, 0.1)
  )
  ls <- cbind(
    lon = -79.4 + runif(3, -0.1, 0.1),
    lat = 43.7 + runif(3, -0.1, 0.1)
  )
  ne <- 40
  ex <- cbind(x1 = rnorm(ne), x2 = runif(ne))
  ee <- runif(ne, 3, 8)
  ey <- rpois(ne, ee * exp(0.2 + 0.4 * ex[, 1] - 0.5 * ex[, 2])) *
    rbinom(ne, 1, 0.8)
  i <- 0:29
  j <- 0:35
  P <- 1 * (abs(outer(1:12, 1:12, "-")) == 1)
  P[1, 12] <- P[12, 1] <- 1
  P[3, 9] <- P[9, 3] <- 1
  list(
    obs = obs,
    rec = rec,
    src = src,
    sw = c(2, 1, 3, 1, 0.5),
    ll = ll,
    ls = ls,
    A = rbind(
      c(0, 1, 1, 0, 0),
      c(1, 0, 1, 1, 0),
      c(1, 1, 0, 1, 1),
      c(0, 1, 1, 0, 1),
      c(0, 0, 1, 1, 0)
    ),
    ex = ex,
    ee = ee,
    ey = ey,
    hop = rbind(
      c(0, 1, 2, -1),
      c(1, 0, 1, -1),
      c(2, 1, 0, -1),
      c(-1, -1, -1, 0)
    ),
    W3 = rbind(c(0, 1, 0), c(1, 0, 1), c(0, 1, 0)),
    fx = matrix(sin(i)),
    ft = 0.5 + ((i * 7) %% 11) / 5 * exp(-0.5 * sin(i) - 0.3 * (i %% 3)),
    fe = as.numeric(i %% 5 != 0),
    fr = i %% 3 + 1,
    cx = matrix(sin(j)),
    ct = ifelse(j %% 4 == 0, 5, 0.2 + ((j * 7) %% 11) / 8),
    ce = as.numeric(!(j %% 4 == 0 | j %% 9 == 1)),
    cr = j %% 3 + 1,
    P = P,
    fv = sin(1:12 / 2) + 0.3 * cos(1:12),
    fz = cos(1:12 / 3)
  )
}

# rmorie's great-circle distances are in km on a radius of 6371 km,
# crimlab's in metres on the mean radius 6371008.8 m: a rmorie km is km
# metres
km <- 6371008.8 / 6371

test_that("the functions agree with rmorie", {
  k <- sin_()
  ref <- refget("sep")$rm
  expect_equal(
    c(
      sum(k$obs),
      sum(k$rec),
      sum(k$src),
      sum(k$ll),
      sum(k$ls),
      sum(k$ex),
      sum(k$ee),
      sum(k$ey)
    ),
    ref$inp,
    tolerance = 1e-12,
    info = "inputs changed"
  )
  b <- bayout(k$obs)
  expect_equal(b$ub, ref$bay$upperbound, tolerance = 0)
  expect_identical(b$alarm, ref$bay$alarm)
  expect_identical(b$t, 7:70)
  b2 <- bayout(k$obs, 10, 2, 2, acty = FALSE, alpha = 0.1, at = 60:70)
  expect_equal(b2$ub, ref$bay2$upperbound, tolerance = 0)
  expect_identical(b2$alarm, ref$bay2$alarm)
  be <- bufexp(k$rec, src = k$src, r = 0.4, w = k$sw)
  expect_equal(be$count, unname(ref$buf$count), tolerance = 0)
  expect_equal(be$wsum, ref$buf$weighted, tolerance = 1e-14)
  expect_equal(be$near, ref$buf$nearest, tolerance = 1e-14)
  bl <- bufexp(k$ll, src = k$ls, r = 8 * km, lonlat = TRUE)
  expect_equal(bl$count, unname(ref$bufll$count), tolerance = 0)
  expect_equal(bl$near, ref$bufll$nearest * km, tolerance = 1e-12)
  expect_equal(
    kerexp(k$rec, src = k$src, h = 0.3, w = k$sw),
    ref$kg,
    tolerance = 1e-13
  )
  expect_equal(
    kerexp(k$rec, src = k$src, h = 0.5, kern = "quartic"),
    ref$kq,
    tolerance = 1e-13
  )
  expect_equal(lerprc(k$A, 0.7, 2), ref$ler, tolerance = 1e-15)
  y <- bym2sc(k$A, 2, 0.4)
  expect_equal(y$scale, ref$bym$scaling_factor, tolerance = 1e-12)
  expect_equal(y$ginv, ref$bym$generalized_inverse, tolerance = 1e-12)
  expect_equal(y$cov, ref$bym$covariance, tolerance = 1e-12)
  p <- ecolrg(k$ey, k$ee, k$ex)
  expect_equal(p$coef, ref$pois$coefficients, tolerance = 1e-12)
  expect_equal(p$fitted, ref$pois$fitted, tolerance = 1e-12)
  expect_equal(p$loglik, ref$pois$loglik, tolerance = 1e-12)
  # the same iterations as rmorie: theta and the coefficients to rounding
  n <- ecolrg(k$ey, k$ee, k$ex, "negbin")
  expect_equal(n$coef, ref$nb$coefficients, tolerance = 1e-10)
  expect_equal(n$theta, ref$nb$theta, tolerance = 1e-10)
  expect_equal(n$loglik, ref$nb$loglik, tolerance = 1e-10)
  z <- ecolrg(k$ey, k$ee, k$ex, "zip", z = k$ex[, 1])
  expect_equal(z$coef, ref$zip$coefficients, tolerance = 1e-10)
  expect_equal(z$zcoef, ref$zip$zero_coefficients, tolerance = 1e-10)
  expect_equal(z$loglik, ref$zip$loglik, tolerance = 1e-10)
  expect_equal(spbias(k$hop, c(1, 0.5, 0.25), -5)$bias, ref$bias, tolerance = 0)
  s <- spsir(
    c(5, 0, 0),
    k$W3,
    c(1000, 800, 1200),
    beta = 0.5,
    gamma = 0.2,
    tmax = 20,
    dt = 0.5
  )
  expect_equal(s$attack, ref$sir$attack, tolerance = 1e-13)
  expect_equal(s$S, ref$sir$S, tolerance = 1e-13)
  expect_equal(s$I, ref$sir$I, tolerance = 1e-13)
  expect_equal(s$R, ref$sir$R, tolerance = 1e-13)
  # Newton-Raphson to a step of 1e-10, as rmorie: the estimates agree to
  # about that, the standard errors (from a numerical Hessian) less closely
  f <- sfrail(k$ft, k$fe, k$fx, k$fr, k$W3)
  expect_equal(f$pll, ref$fr$statistic, tolerance = 1e-10)
  expect_equal(f$shape, ref$fr$shape, tolerance = 1e-9)
  expect_equal(f$coef, ref$fr$coefficients, tolerance = 1e-9)
  expect_equal(f$frail, ref$fr$frailties, tolerance = 1e-9)
  expect_equal(f$se, ref$fr$se, tolerance = 1e-6)
  expect_identical(f$iter, ref$fr$iterations)
  cu <- sfcure(k$ct, k$ce, k$cx, k$cr, k$W3)
  expect_equal(cu$pll, ref$cr$statistic, tolerance = 1e-10)
  expect_equal(
    c(cu$shape, cu$lscale, cu$coef),
    c(ref$cr$shape, ref$cr$log_scale, ref$cr$coefficients),
    tolerance = 1e-9
  )
  expect_equal(cu$effect, ref$cr$regional_effects, tolerance = 1e-9)
  expect_equal(cu$se, ref$cr$se, tolerance = 1e-6)
  expect_equal(cu$cure, ref$cr$cure_fraction, tolerance = 1e-9)
  # eigenvector indices: crimlab from 1, rmorie's spatial_filter() from 0
  e1 <- evfilt(k$fv, k$P)
  expect_identical(e1$sel, ref$sf1$selected + 1L)
  expect_equal(e1$coef, ref$sf1$coefficients, tolerance = 1e-10)
  expect_equal(
    c(e1$imoran, e1$aic, e1$r2, e1$fr2),
    c(ref$sf1$residual_moran, ref$sf1$aic, ref$sf1$r2, ref$sf1$filter_r2),
    tolerance = 1e-10
  )
  expect_equal(e1$resid, ref$sf1$residuals, tolerance = 1e-10)
  e2 <- evfilt(k$fv, k$P, k$fz, "aic")
  expect_identical(e2$sel, ref$sf2$selected + 1L)
  expect_equal(e2$coef, ref$sf2$coefficients, tolerance = 1e-10)
  expect_equal(
    c(e2$imoran, e2$aic, e2$r2, e2$fr2),
    c(ref$sf2$residual_moran, ref$sf2$aic, ref$sf2$r2, ref$sf2$filter_r2),
    tolerance = 1e-10
  )
  e3 <- evfilt(k$fv, k$P, method = "r2", tol = 0.01)
  expect_identical(e3$sel, ref$sf3$selected + 1L)
  expect_equal(e3$r2, ref$sf3$r2, tolerance = 1e-10)
  expect_equal(e1$mi, ref$mi, tolerance = 1e-12)
})

test_that("the functions agree with the reference packages", {
  k <- sin_()
  ref <- refget("sep")$oth
  b <- bayout(k$obs)
  expect_equal(b$ub, ref$bay, tolerance = 0)
  expect_identical(b$alarm, ref$bayal)
  # rdist.earth() takes the arc cosine of a dot product, which loses about
  # 1e-11 of a distance of 10 km to rounding
  bl <- bufexp(k$ll, src = k$ls, r = 1, lonlat = TRUE)
  expect_equal(bl$near, ref$near * 1000, tolerance = 1e-9)
  y <- bym2sc(k$A, 2, 0.4)
  expect_equal(y$ginv, ref$ginv, tolerance = 1e-12)
  expect_equal(y$scale, ref$scale, tolerance = 1e-12)
  p <- ecolrg(k$ey, k$ee, k$ex)
  expect_equal(p$coef, ref$pois, tolerance = 1e-10)
  expect_equal(p$loglik, ref$poll, tolerance = 1e-12)
  # glm.nb() alternates theta.ml() and glm.fit() to its own tolerances
  n <- ecolrg(k$ey, k$ee, k$ex, "negbin")
  expect_equal(n$coef, ref$nb, tolerance = 1e-6)
  expect_equal(n$theta, ref$theta, tolerance = 1e-5)
  # zeroinfl() maximises the same likelihood by BFGS
  z <- ecolrg(k$ey, k$ee, k$ex, "zip", z = k$ex[, 1])
  expect_equal(c(z$coef, z$zcoef), ref$zip, tolerance = 1e-5)
  expect_equal(z$loglik, ref$zll, tolerance = 1e-8)
  s <- spsir(
    c(5, 0, 0),
    k$W3,
    c(1000, 800, 1200),
    beta = 0.5,
    gamma = 0.2,
    tmax = 20,
    dt = 0.5
  )
  expect_equal(cbind(s$S, s$I, s$R), ref$rk4, tolerance = 1e-12)
  # one area without a penalty: the Weibull model of survreg(), its
  # coefficient and intercept on the proportional-hazards scale
  f1 <- sfrail(k$ft, k$fe, k$fx, rep(1, 30), matrix(0, 1, 1), tau = 0)
  expect_equal(
    c(f1$shape, f1$coef, f1$frail),
    ref$weib[c(1, 3, 2)],
    tolerance = 1e-7
  )
  expect_equal(c(f1$shape, f1$coef, f1$frail), ref$fr1, tolerance = 1e-9)
  e1 <- evfilt(k$fv, k$P)
  expect_equal(e1$imoran, ref$moran, tolerance = 1e-10)
  expect_equal(e1$coef, ref$lmcoef, tolerance = 1e-10)
})

test_that("rmorie's own checks hold", {
  u <- (1:100) / 101
  obs <- as.integer(5 + 5 * u[1:30])
  r <- bayout(obs, w = 5)
  for (q in seq_along(r$t)) {
    t <- r$t[q]
    expect_equal(r$ub[q], qnbinom(0.95, sum(obs[(t - 5):(t - 1)]) + 0.5, 5 / 6),
                 tolerance = testthat_tolerance())
  }
  A <- rbind(c(0, 1, 1, 0), c(1, 0, 1, 0), c(1, 1, 0, 1), c(0, 0, 1, 0))
  expect_equal(rowSums(lerprc(A, 1)), rep(0, 4), tolerance = 1e-15)
  b <- bym2sc(A, 2, 0.4)
  Q <- diag(rowSums(A)) - A
  expect_equal(Q %*% b$ginv %*% Q, Q, tolerance = 1e-12)
  expect_equal(exp(mean(log(diag(b$ginv) / b$scale))), 1, tolerance = 1e-12)
  k <- kerexp(
    rbind(c(0.2, 0.3)),
    src = rbind(c(0.25, 0.35), c(0.6, 0.8)),
    h = 0.5,
    kern = "quartic"
  )
  d <- sqrt(c(0.05^2 + 0.05^2, 0.4^2 + 0.5^2)) / 0.5
  expect_equal(
    k,
    sum(ifelse(d < 1, 3 / pi * (1 - d^2)^2, 0)) / 0.25,
    tolerance = 1e-13
  )
  W <- rbind(c(0, 1, 0, 0), c(1, 0, 1, 0), c(0, 1, 0, 1), c(0, 0, 1, 0))
  N <- c(1000, 800, 1200, 500)
  s <- spsir(c(5, 0, 0, 0), W, N, beta = 0.5, gamma = 0.2, tmax = 40, dt = 0.25)
  expect_equal(s$S[161, ] + s$I[161, ] + s$R[161, ], N, tolerance = 1e-12)
  iso <- spsir(
    c(5, 0, 0, 0),
    NULL,
    N,
    beta = 0.5,
    gamma = 0.2,
    tmax = 40,
    dt = 0.25
  )
  expect_equal(iso$R[161, 3], 0, tolerance = testthat_tolerance())
  C <- 0.9 * diag(4) + 0.1 * W / rowSums(W)
  f <- function(y) {
    S <- y[1:4]
    I <- y[5:8]
    l <- as.vector(C %*% (I / N))
    c(-0.5 * S * l, 0.5 * S * l - 0.2 * I, 0.2 * I)
  }
  y0 <- c(N - c(5, 0, 0, 0), c(5, 0, 0, 0), rep(0, 4))
  h <- 0.25
  k1 <- f(y0)
  k2 <- f(y0 + h / 2 * k1)
  k3 <- f(y0 + h / 2 * k2)
  k4 <- f(y0 + h * k3)
  y1 <- y0 + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
  s1 <- spsir(c(5, 0, 0, 0), W, N, beta = 0.5, gamma = 0.2, tmax = h, dt = h)
  expect_equal(c(s1$S[2, ], s1$I[2, ], s1$R[2, ]), y1, tolerance = 1e-12)
  W3 <- rbind(c(0, 1, 0), c(1, 0, 1), c(0, 1, 0))
  i <- 0:29
  reg <- i %% 3
  x <- matrix(sin(i))
  t <- 0.5 + ((i * 7) %% 11) / 5 * exp(-0.5 * x[, 1] - 0.3 * reg)
  ev <- as.numeric(i %% 5 != 0)
  r <- sfrail(t, ev, x, reg + 1, W3)
  eta <- x[, 1] * r$coef + r$frail[reg + 1]
  H <- t^r$shape * exp(eta)
  expect_lt(abs(sum((ev - H) * x[, 1])), 1e-7)
  Q3 <- diag(rowSums(W3)) - W3
  gu <- vapply(0:2, function(k) sum((ev - H)[reg == k]), 0) -
    as.vector(Q3 %*% r$frail)
  expect_lt(max(abs(gu)), 1e-7)
  j <- 0:35
  x2 <- matrix(sin(j))
  t2 <- ifelse(j %% 4 == 0, 5, 0.2 + ((j * 7) %% 11) / 8)
  ev2 <- as.numeric(!(j %% 4 == 0 | j %% 9 == 1))
  cr <- sfcure(t2, ev2, x2, j %% 3 + 1, W3)
  th <- exp(x2[, 1] * cr$coef + cr$effect[j %% 3 + 1])
  G <- exp(cr$lscale + cr$shape * log(t2))
  expect_lt(abs(sum((ev2 - th * (1 - exp(-G))) * x2[, 1])), 1e-7)
  expect_equal(cr$cure, exp(-th), tolerance = 1e-12)
  # filters recompute: the residual Moran's I of the selected vectors' fit
  Wp <- 1 * (abs(outer(1:12, 1:12, "-")) == 1)
  yp <- c(1, 2, 2.5, 4, 3.5, 3, 1.5, 1, 2, 3.5, 4, 2.5)
  e <- evfilt(yp, Wp, method = "aic")
  d <- e$resid - mean(e$resid)
  expect_equal(
    e$imoran,
    12 / sum(Wp) * sum(Wp * outer(d, d)) / sum(d^2),
    tolerance = 1e-12
  )
})

test_that("results print and summarise", {
  k <- sin_()
  expect_output(print(ecolrg(k$ey, k$ee, k$ex)), "ecolrg: 40 counts")
  expect_named(
    ecolrg(k$ey, k$ee, as.data.frame(k$ex), "negbin"),
    c("coef", "theta", "fitted", "loglik")
  )
  expect_output(
    print(spsir(c(1, 0), tmax = 1, dt = 0.5)),
    "spsir: 2 patches to time 1, attack rate"
  )
  f <- sfrail(k$ft, k$fe, k$fx, k$fr, k$W3)
  expect_output(print(f), "sfrail: Weibull frailty model")
  sf <- summary(f)
  expect_identical(sf$term, c("logshape", "x1", "frail1", "frail2", "frail3"))
  expect_equal(sf$z, sf$est / sf$se, tolerance = testthat_tolerance())
  c0 <- sfcure(k$ct, k$ce, NULL, k$cr, k$W3)
  expect_length(c0$coef, 0)
  expect_output(print(c0), "sfcure: promotion-time cure model")
  expect_identical(
    summary(c0)$term,
    c("logshape", "lscale", "effect1", "effect2", "effect3")
  )
  expect_output(print(evfilt(k$fv, k$P)), "eigenvectors selected")
  expect_length(evfilt(k$fv, k$P, maxv = 0)$sel, 0)
  expect_length(evfilt(k$fv, k$P, cand = "all", method = "press")$mi, 11)
  expect_named(
    evfilt(k$fv, k$P, as.data.frame(k$fz), cand = "negative"),
    c("sel", "coef", "imoran", "aic", "r2", "fr2", "resid", "mi")
  )
  lw <- spwt(nb = list(2, c(1, 3), 2), type = "contig", style = "raw")
  expect_equal(lerprc(lw, 0.5), lerprc(k$W3, 0.5),
               tolerance = testthat_tolerance())
  expect_equal(bym2sc(lw)$scale, bym2sc(k$W3)$scale,
               tolerance = testthat_tolerance())
})

test_that("spatial epidemiology refuses bad input in words", {
  k <- sin_()
  errs(
    "bayout",
    list(
      list("`x` must be a numeric vector of counts", "a"),
      list("`x` must be counts from 0, or NA", c(1, -1)),
      list("the reference window needs at least 7", 1:3),
      list("`at` must be times from 1 to 8", 1:8, 52, 0, 6, TRUE, 0.05, 9),
      list("time 2 is too early", 1:8, 52, 0, 2, TRUE, 0.05, 2),
      list("no observed count", c(rep(NA, 6), 1), 52, 0, 6),
      list(
        "`alpha` must be one number strictly between 0 and 1",
        1:8,
        52,
        0,
        6,
        TRUE,
        1
      )
    )
  )
  errs(
    "bufexp",
    list(
      list("`r` must be at least 0", k$rec, NULL, k$src, -1),
      list("`w` must have length 5", k$rec, NULL, k$src, 1, 1:2),
      list("`lonlat` must be TRUE or FALSE", k$rec, NULL, k$src, 1, NULL, NA)
    )
  )
  errs(
    "kerexp",
    list(
      list(
        "`h` must be one number strictly between 0 and Inf",
        k$rec,
        NULL,
        k$src,
        0
      ),
      list("`kern` must be one of", k$rec, NULL, k$src, 1, NULL, "epan")
    )
  )
  errs(
    "lerprc",
    list(
      list("`rho` must be one number from 0 to 1", k$A, 2),
      list("`w` must be the list spwt\\(\\) returns", 1, 0.5),
      list(
        "`w` must be a 5 x 5 matrix",
        structure(k$A[, 1:4], dim = c(5L, 4L)),
        0.5
      )
    )
  )
  errs(
    "bym2sc",
    list(
      list("at least 2 areas", matrix(0, 1, 1)),
      list(
        "every area needs a neighbour",
        rbind(c(0, 1, 0), c(1, 0, 0), c(0, 0, 0))
      )
    )
  )
  errs(
    "ecolrg",
    list(
      list("`e` must be greater than 0", 1:3, c(1, 0, 1)),
      list("`family` must be one of", 1:3, c(1, 1, 1), NULL, "gamma"),
      list(
        "more values than the model has coefficients",
        1:2,
        c(1, 1),
        cbind(1:2)
      ),
      list("singular", 1:6, rep(1, 6), cbind(1:6, 2 * (1:6))),
      list(
        "`maxit` must be one whole number from 1",
        1:3,
        c(1, 1, 1),
        NULL,
        "poisson",
        NULL,
        1e-12,
        0
      )
    )
  )
  errs(
    "spbias",
    list(
      list("`d` must be a numeric matrix of hop distances", 1:3, 1),
      list("`d` must hold whole numbers from 0, or -1", rbind(c(0, 1.5)), 1),
      list("`d` must hold whole numbers from 0, or -1", rbind(c(0, -2)), 1)
    )
  )
  errs(
    "spsir",
    list(
      list("`i0` must not exceed `pop`", c(5, 1), NULL, c(3, 3)),
      list(
        "`kappa` must be one number from 0 to 1",
        1,
        NULL,
        NULL,
        0.5,
        0.2,
        2
      ),
      list("at most 1e6 steps", 1, NULL, NULL, 0.5, 0.2, 0.1, 1e7, 1)
    )
  )
  errs(
    "sfrail",
    list(
      list(
        "`reg` must be area numbers from 1 to 3",
        k$ft,
        k$fe,
        k$fx,
        k$fr + 1,
        k$W3
      ),
      list("`event` must be 0 or 1", k$ft, k$fe + 2, k$fx, k$fr, k$W3),
      list("`x` must have 30 rows", k$ft, k$fe, 1:3, k$fr, k$W3),
      list(
        "the Hessian is singular",
        k$ft,
        k$fe,
        rep(1, 30),
        rep(1, 30),
        matrix(0, 1, 1),
        0
      )
    )
  )
  errs(
    "sfcure",
    list(list(
      "`tol` must be one number strictly between 0 and 1",
      k$ct,
      k$ce,
      k$cx,
      k$cr,
      k$W3,
      1,
      2
    ))
  )
  errs(
    "evfilt",
    list(
      list("`v` must have at least 3 values", 1:2, diag(2)),
      list("`method` must be one of", k$fv, k$P, NULL, "bic"),
      list("`cand` must be one of", k$fv, k$P, NULL, "moran", "some"),
      list("no links", k$fv, matrix(0, 12, 12)),
      list(
        "`maxv` must be one whole number from 0",
        k$fv,
        k$P,
        NULL,
        "moran",
        "positive",
        0.1,
        -1
      )
    )
  )
})

test_that("data frames of covariates are read as matrices", {
  k <- sin_()
  a <- ecolrg(k$ey, k$ee, k$ex, "zip", z = data.frame(z = k$ex[, 1]))
  b <- ecolrg(k$ey, k$ee, k$ex, "zip", z = k$ex[, 1, drop = FALSE])
  expect_equal(a$zcoef, b$zcoef, tolerance = 1e-12)
  f1 <- sfrail(k$ft, k$fe, data.frame(x = k$fx[, 1]), k$fr, k$W3)
  f2 <- sfrail(k$ft, k$fe, k$fx, k$fr, k$W3)
  expect_equal(f1$coef, f2$coef, tolerance = 1e-12)
  c1 <- sfcure(k$ct, k$ce, data.frame(x = k$cx[, 1]), k$cr, k$W3)
  c2 <- sfcure(k$ct, k$ce, k$cx, k$cr, k$W3)
  expect_equal(c1$coef, c2$coef, tolerance = 1e-12)
})

test_that("bayout() bounds large counts and refuses counts past any bound", {
  # six weeks of 5000: p^s underflows, so the search from 0 never moved
  # and the call ran for ever; the bound is now R's quantile
  r <- bayout(rep(5000, 20))
  expect_equal(r$ub, rep(qnbinom(0.95, 30000.5, 6 / 7), 14), tolerance = 0)
  expect_error(bayout(rep(1e300, 10)),
               paste("bayout: the counts of a reference window are too large",
                     "for the negative binomial bound."), fixed = TRUE)
})

test_that("evfilt() on values whose sums overflow adds no eigenvector", {
  # every criterion is NaN: no candidate is best, and the selection used
  # to add an eigenvector that does not exist
  k <- sin_()
  e <- evfilt(k$fv * 1e300, k$P)
  expect_length(e$sel, 0L)
})
